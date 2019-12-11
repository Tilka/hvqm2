#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hvqm2dec.h>

#define MAXWIDTH 640
#define MAXHEIGHT 480
#define HVQM_DATASIZE_MAX 30000
#define OUTPUT_DIR "output"

#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-parameter"

_Static_assert(sizeof(HVQM2Frame) == 0x34, "");

typedef struct BitBuffer
{
    u32 bit;
    u32 value;
    u32 const *pos;
#ifdef NATIVE
    u32 size;
#endif
} BitBuffer;

static u16 read16(void const *src)
{
    u8 const *ptr = (u8 const *)src;
    u16 value = 0;
    value |= *ptr++;
    value <<= 8;
    value |= *ptr;
    return value;
}

static u32 read32(void const *src)
{
    u8 const *ptr = (u8 const *)src;
    u32 value = 0;
    value |= *ptr++;
    value <<= 8;
    value |= *ptr++;
    value <<= 8;
    value |= *ptr++;
    value <<= 8;
    value |= *ptr;
    return value;
}

typedef struct Tree
{
    s16 root;
    u16 array[2][0x200];
} Tree;

static struct
{
    u8 nestbuf[HVQM2_NESTSIZE_L * HVQM2_NESTSIZE_S];
    BitBuffer macroblock_buf;
    BitBuffer movevector_buf;
    BitBuffer basisnum_buf[2];
    BitBuffer basisnumrun_buf[2];
    BitBuffer scale_buf[3];
    BitBuffer dcval_buf[3];
    BitBuffer dcrun_buf[3];
    void const *fixvl[3];
    Tree dcval_tree;
    Tree basisnum_tree;
    Tree basisnumrun_tree;
    Tree scale_tree;
    Tree movevector_tree;
    u16 Huff_nodeno;
    s16 dc_max;
    s16 dc_min;
    u8 *basis[3];
    u8 *dcbuf[3];
    struct wcode
    {
        u8 *basis_prev_line;
        u8 *dcbuf_prev_line;
        u8 *basis_curr_line;
        u8 *dcbuf_curr_line;
        u8 *basis_next_line;
        u8 *dcbuf_next_line;
        u8 basis_next;
        u8 dcbuf_next;
        u8 basis_curr;
        u8 dcbuf_curr;
        u8 dcbuf_prev;
    } u_wcode, v_wcode, y0wcode, y1wcode;
    u32 divT[0x200];
    u8 clipT[0x300];
    void (*ColorConv)(u32 *outbuf, u16 const *pix_y, u16 const *pix_u, u16 const *pix_v);
    s32 fb_width;
    u32 h_sampling_rate;
    u32 v_sampling_rate;
    u32 mcu_h_pix;
    u32 mcu_v_pix;
    u32 mcu411;
    u32 next_macroblk_line;
    u32 im_width;
    u32 im_height;
    u32 lum_hblocks;
    u32 lum_vblocks;
    s32 lum_totalblocks;
    u32 col_hblocks;
    u32 col_vblocks;
    s32 col_totalblocks;
    u8 *nest;
    u32 nest_w;
    u32 nest_h;
    u8 yshift;
    u8 pix_alpha;
    u16 pix_y[4*16];
    u16 pix_u[16];
    u16 pix_v[16];
} global;

static void dumpPlane(char const *path, u8 const *data, int plane_idx)
{
    FILE *f = fopen(path, "wb+");
    if (!f)
    {
        fprintf(stderr, "failed opening '%s'", path);
        exit(-1);
    }
    u32 hblocks = plane_idx ? global.col_hblocks : global.lum_hblocks;
    u32 vblocks = plane_idx ? global.col_vblocks : global.lum_vblocks;
    fprintf(f, "P5\n%u %u\n255\n", hblocks, vblocks);
    for (u32 i = 0; i < vblocks; ++i)
    {
        for (u32 j = 0; j < hblocks; ++j)
        {
            fputc(*data++, f);
        }
    }
    fclose(f);
}

__attribute__((unused))
static void dumpPlanes(u32 frame, u32 frame_type)
{
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        char path[128];
        snprintf(path, 128, OUTPUT_DIR"/basis_%u_%c_%c.ppm", frame, "iph"[frame_type], "yuv"[plane_idx]);
        dumpPlane(path, global.basis[plane_idx], plane_idx);
        snprintf(path, 128, OUTPUT_DIR "/dcbuf_%u_%c_%c.ppm", frame, "iph"[frame_type], "yuv"[plane_idx]);
        dumpPlane(path, global.dcbuf[plane_idx], plane_idx);
    }
}

static u32 getBit(BitBuffer *buf)
{
    u32 bit = buf->bit;
    if (bit == 0)
    {
#ifdef NATIVE
        if (!buf->size)
        {
            fprintf(stderr, "error: BitBuffer overread\n");
            //exit(1);
        }
        buf->size -= 4;
#endif
        buf->value = read32(buf->pos);
        buf->pos += 1;
        bit = 0x80000000;
    }
    buf->bit = bit >> 1;
    return buf->value & bit ? 1 : 0;
}

static u32 getByte(BitBuffer *buf)
{
    u32 result = 0;
    for (int i = 7; i >= 0; --i)
        result |= getBit(buf) << i;
    return result;
}

// done
// HVQM4: _readTree
static s16 readHuffTree(BitBuffer *src, Tree *dst)
{
    if (getBit(src) == 0)
    {
        // leaf node
        return getByte(src);
    }
    else
    {
        // recurse
        s16 pos = global.Huff_nodeno++;
        // read the 0 side of the tree
        dst->array[0][pos] = readHuffTree(src, dst);
        // read the 1 side of the tree
        dst->array[1][pos] = readHuffTree(src, dst);
        return pos;
    }
}

static void setCode(BitBuffer *dst, void const *src)
{
    u32 size = read32(src);
    dst->pos = size ? src + 4 : NULL;
    dst->bit = 0;
#ifdef NATIVE
    dst->size = size;
#endif
}

static void readTree(BitBuffer *buf, Tree *tree)
{
    global.Huff_nodeno = 0x100;
    if (buf->pos)
        tree->root = readHuffTree(buf, tree);
}

// done
static s16 decodeHuff(BitBuffer *buf, Tree *tree)
{
    s16 pos = tree->root;
    while (pos >= 0x100)
        pos = tree->array[getBit(buf)][pos];
    return pos;
}

static s16 decodeHuff2(BitBuffer *buf, Tree *tree)
{
    s16 pos = tree->root;
    while (pos >= 0x100)
        pos = tree->array[getBit(buf)][pos];
    return tree->array[0][pos];
}

// done
// HVQM4: kinda like decodeSOvfSym
static s16 decodeDC(BitBuffer *buf)
{
    s16 value = decodeHuff2(buf, &global.dcval_tree);
    s16 sum = value;
    // '=='!
    if (value == global.dc_min || value == global.dc_max)
    {
        do
        {
            value = decodeHuff2(buf, &global.dcval_tree);
            sum += value;
        } while (value <= global.dc_min || value >= global.dc_max);
    }
    return sum;
}

// done
static u8 getDeltaBN(u8 *rle, BitBuffer *val_buf, BitBuffer *rle_buf)
{
    if (*rle)
    {
        --*rle;
        return 0;
    }
    u8 value = decodeHuff(val_buf, &global.basisnum_tree);
    if (value == 0)
        *rle = decodeHuff(rle_buf, &global.basisnumrun_tree);
    return value;
}

// TODO: the original decoder is making a trade-off here,
// we can improve the interpolation quality by using additional memory
static void WeightImBlock(u16 *pix, u8 curr, u8 top, u8 bottom, u8 left, u8 right)
{
    u32 const tmb = top - bottom;
    u32 const lmr = left - right;
    u32 const vph = tmb + lmr;
    u32 const vmh = tmb - lmr;

    u32 const c2 = curr * 2;
    u32 const c8 = (curr * 8) + 4;

    u32 const tpl = top    + left  - c2;
    u32 const tpr = top    + right - c2;
    u32 const bpr = bottom + right - c2;
    u32 const bpl = bottom + left  - c2;

    u32 const tml = top    - left;
    u32 const tmr = top    - right;
    u32 const bmr = bottom - right;
    u32 const bml = bottom - left;

    pix[0] = (c8 + vph + tpl) / 8;
    pix[1] = (c8 + vph + tml) / 8;
    pix[2] = (c8 + vmh + tmr) / 8;
    pix[3] = (c8 + vmh + tpr) / 8;

    pix[4] = (c8 + vph - tml) / 8;
    pix[5] = (c8 - bpr      ) / 8;
    pix[6] = (c8 - bpl      ) / 8;
    pix[7] = (c8 + vmh - tmr) / 8;

    pix[8]  = (c8 - vmh - bml) / 8;
    pix[9]  = (c8 - tpr      ) / 8;
    pix[10] = (c8 - tpl      ) / 8;
    pix[11] = (c8 - vph - bmr) / 8;

    pix[12] = (c8 - vmh + bpl) / 8;
    pix[13] = (c8 - vmh + bml) / 8;
    pix[14] = (c8 - vph + bmr) / 8;
    pix[15] = (c8 - vph + bpr) / 8;
}

// HVQM4: kinda like IpicBlockDec
static void decBlockCPU(u16 *pix, struct wcode *wcode, u32 plane_idx)
{
    if (wcode->basis_curr == 0)
    {
        //putchar('W');
        u8 curr   =  wcode->dcbuf_curr;
        u8 right  =  wcode->basis_next      == 0 ?  wcode->dcbuf_next      : curr;
        u8 top    = *wcode->basis_prev_line == 0 ? *wcode->dcbuf_prev_line : curr;
        u8 bottom = *wcode->basis_next_line == 0 ? *wcode->dcbuf_next_line : curr;
        u8 left   =  wcode->dcbuf_prev;
        WeightImBlock(pix, curr, top, bottom, left, right);
        wcode->dcbuf_prev = wcode->dcbuf_curr;
    }
    else
    {
        if (wcode->basis_curr == 8)
        {
            // HVQM4: OrgBlock (but on block type 6)
            //putchar('O');
            for (u32 i = 0; i < 16; ++i)
                *pix++ = *(u8*)global.fixvl[plane_idx]++;
        }
        else
        {
            // HVQM4: IntraAotBlock
            //putchar('A');
            for (u32 i = 0; i < 16; ++i)
                pix[i] = wcode->dcbuf_curr;
            for (u32 i = 0; i < wcode->basis_curr; ++i)
            {
                // aka F4
                u16 bits = read16(global.fixvl[plane_idx]);
                global.fixvl[plane_idx] += 2;
                u32 x_stride = ((bits >> 0) & 1) + 1;
                u32 y_stride = ((bits >> 1) & 1) + 1;
                u32 nest_pos_x, nest_pos_y;
                if (global.yshift == 8)
                {
                    nest_pos_y = ((bits >> 8) & 0x1F);
                    nest_pos_x = ((bits >> 2) & 0x3F);
                }
                else
                {
                    nest_pos_y = ((bits >> 7) & 0x3F);
                    nest_pos_x = ((bits >> 2) & 0x1F);
                }
                s16 tmp[4][4];
                u8 const *nest = global.nest + nest_pos_y * global.nest_w + nest_pos_x;
                s32 sum = 0;
                for (int y = 0; y < 4; ++y)
                {
                    for (int x = 0; x < 4; ++x)
                    {
                        u8 nest_value = nest[y * y_stride * global.nest_w + x * x_stride];
                        tmp[y][x] = nest_value;
                        sum += nest_value;
                    }
                }
                // FIXME?
                s32 mean = (sum + 8) >> 4;
                s16 max = 0;
                for (int y = 0; y < 4; ++y)
                {
                    for (int x = 0; x < 4; ++x)
                    {
                        tmp[y][x] -= mean;
                        s16 value = tmp[y][x];
                        if (value < 0)
                            value = -value;
                        if (max < value)
                            max = value;
                    }
                }
                s32 beta_q = bits >> 13;
                s32 scale = decodeHuff2(&global.scale_buf[plane_idx], &global.scale_tree);
                s32 div = global.divT[max];
                s32 bar = (beta_q + scale) * div;
                u16 *out = pix;
                for (int y = 0; y < 4; ++y)
                    for (int x = 0; x < 4; ++x)
                        *out++ += (tmp[y][x] * bar + 512) >> 10;
            }
            for (int j = 0; j < 16; ++j)
                pix[j] = 0;
        }
        wcode->dcbuf_prev = wcode->dcbuf_next;
    }
    wcode->basis_prev_line++;
    wcode->dcbuf_prev_line++;
    wcode->basis_next_line++;
    wcode->dcbuf_next_line++;
}

static u32 yuv2rgba(s16 y, s16 u, s16 v)
{
    // ITU-R formula in 10.6-bit fixed-point
    u -= 128;
    v -= 128;
    s32 y6 = y << 6;
    s16 r16 = 0x4020 + v * 90;
    u8 r = global.clipT[(y6 + r16) >> 6];
    s16 g16 = 0x4020 - v * 46 - u * 22;
    u8 g = global.clipT[(y6 + g16) >> 6];
    s16 b16 = 0x4020 + u * 113;
    u8 b = global.clipT[(y6 + b16) >> 6];
    u8 a = global.pix_alpha;
    return a << 24 | b << 16 | g << 8 | r;
}

static void ColorConv422(u32 *outbuf, u16 const *pix_y, u16 const *pix_u, u16 const *pix_v)
{
    u16 const *pix_y0 = pix_y;
    u16 const *pix_y1 = pix_y + 16;
    for (u32 i = 0; i < 4; ++i)
    {
        u32 *out = outbuf;
        for (u32 j = 0; j < 2; ++j)
        {
            *out++ = yuv2rgba(*pix_y0++, *pix_u,   *pix_v);
            *out++ = yuv2rgba(*pix_y0++, *pix_u++, *pix_v++);
        }
        for (u32 j = 0; j < 2; ++j)
        {
            *out++ = yuv2rgba(*pix_y1++, *pix_u,   *pix_v);
            *out++ = yuv2rgba(*pix_y1++, *pix_u++, *pix_v++);
        }
        outbuf += global.fb_width;
    }
}

static void ColorConv411(u32 *outbuf, u16 const *pix_y, u16 const *pix_u, u16 const *pix_v)
{
    u16 const *pix_y0 = pix_y;
    u16 const *pix_y1 = pix_y + 16;
    for (u32 k = 0; k < 2; ++k)
    {
        for (u32 i = 0; i < 2; ++i)
        {
            for (u32 n = 0; n < 2; ++n)
            {
                u32 *out = outbuf;
                u16 const *u = pix_u;
                u16 const *v = pix_v;
                for (u32 j = 0; j < 2; ++j)
                {
                    *out++ = yuv2rgba(*pix_y0++, *u,   *v);
                    *out++ = yuv2rgba(*pix_y0++, *u++, *v++);
                }
                for (u32 j = 0; j < 2; ++j)
                {
                    *out++ = yuv2rgba(*pix_y1++, *u,   *v);
                    *out++ = yuv2rgba(*pix_y1++, *u++, *v++);
                }
                outbuf += global.fb_width;
            }
            pix_u += 4;
            pix_v += 4;
        }
        pix_y0 += 16;
        pix_y1 += 16;
    }
}

static void advance_wcode(struct wcode *wcode, u32 amount)
{
    wcode->basis_prev_line += amount;
    wcode->dcbuf_prev_line += amount;
    wcode->basis_curr_line += amount;
    wcode->dcbuf_curr_line += amount;
    wcode->basis_next_line += amount;
    wcode->dcbuf_next_line += amount;
}

static void update_wcode(struct wcode *wcode)
{
    wcode->basis_curr = wcode->basis_next;
    wcode->dcbuf_curr = wcode->dcbuf_next;
    wcode->basis_curr_line++;
    wcode->dcbuf_curr_line++;
    wcode->basis_next = *wcode->basis_curr_line;
    wcode->dcbuf_next = *wcode->dcbuf_curr_line;
}

// HVQM4: kinda like IpicLineDec
static void decLine(u32 *outbuf)
{
    global.u_wcode.basis_next = *global.u_wcode.basis_curr_line;
    global.u_wcode.dcbuf_prev = *global.u_wcode.dcbuf_curr_line;
    global.u_wcode.dcbuf_next = *global.u_wcode.dcbuf_curr_line;
    global.v_wcode.basis_next = *global.v_wcode.basis_curr_line;
    global.v_wcode.dcbuf_prev = *global.v_wcode.dcbuf_curr_line;
    global.v_wcode.dcbuf_next = *global.v_wcode.dcbuf_curr_line;
    global.y0wcode.basis_next = *global.y0wcode.basis_curr_line;
    global.y0wcode.dcbuf_prev = *global.y0wcode.dcbuf_curr_line;
    global.y0wcode.dcbuf_next = *global.y0wcode.dcbuf_curr_line;
    if (global.mcu411)
    {
        global.y1wcode.basis_next = *global.y1wcode.basis_curr_line;
        global.y1wcode.dcbuf_prev = *global.y1wcode.dcbuf_curr_line;
        global.y1wcode.dcbuf_next = *global.y1wcode.dcbuf_curr_line;
    }
    for (u32 i = 0; i < global.col_hblocks - 1; ++i)
    {
        if (global.y0wcode.basis_next == 0x80)
        {
            advance_wcode(&global.y0wcode, 2);
            global.y0wcode.basis_next = *global.y0wcode.basis_curr_line;
            global.y0wcode.dcbuf_prev = *global.y0wcode.dcbuf_curr_line;
            global.y0wcode.dcbuf_next = *global.y0wcode.dcbuf_curr_line;
            if (global.mcu411)
            {
                advance_wcode(&global.y1wcode, 2);
                global.y1wcode.basis_next = *global.y1wcode.basis_curr_line;
                global.y1wcode.dcbuf_prev = *global.y1wcode.dcbuf_curr_line;
                global.y1wcode.dcbuf_next = *global.y1wcode.dcbuf_curr_line;
            }
            advance_wcode(&global.u_wcode, 1);
            global.u_wcode.basis_next = *global.u_wcode.basis_curr_line;
            global.u_wcode.dcbuf_prev = *global.u_wcode.dcbuf_curr_line;
            global.u_wcode.dcbuf_next = *global.u_wcode.dcbuf_curr_line;
            advance_wcode(&global.v_wcode, 1);
            global.v_wcode.basis_next = *global.v_wcode.basis_curr_line;
            global.v_wcode.dcbuf_prev = *global.v_wcode.dcbuf_curr_line;
            global.v_wcode.dcbuf_next = *global.v_wcode.dcbuf_curr_line;
        }
        else
        {
            update_wcode(&global.y0wcode);
            decBlockCPU(global.pix_y + 0*16, &global.y0wcode, 0);
            update_wcode(&global.y0wcode);
            decBlockCPU(global.pix_y + 1*16, &global.y0wcode, 0);
            if (global.mcu411)
            {
                update_wcode(&global.y1wcode);
                decBlockCPU(global.pix_y + 2*16, &global.y1wcode, 0);
                update_wcode(&global.y1wcode);
                decBlockCPU(global.pix_y + 3*16, &global.y1wcode, 0);
            }
            update_wcode(&global.u_wcode);
            decBlockCPU(global.pix_u, &global.u_wcode, 1);
            update_wcode(&global.v_wcode);
            decBlockCPU(global.pix_v, &global.v_wcode, 2);
            global.ColorConv(outbuf, global.pix_y, global.pix_u, global.pix_v);
        }
        outbuf += global.mcu_h_pix;
    }
    if (global.y0wcode.basis_next == 0x80)
    {
        advance_wcode(&global.y0wcode, 2);
        if (global.mcu411)
            advance_wcode(&global.y1wcode, 2);
        advance_wcode(&global.u_wcode, 1);
        advance_wcode(&global.v_wcode, 1);
    }
    else
    {
        update_wcode(&global.y0wcode);
        decBlockCPU(global.pix_y + 0*16, &global.y0wcode, 0);
        global.y0wcode.basis_curr = global.y0wcode.basis_next;
        global.y0wcode.dcbuf_curr = global.y0wcode.dcbuf_next;
        global.y0wcode.basis_curr_line++;
        global.y0wcode.dcbuf_curr_line++;
        decBlockCPU(global.pix_y + 1*16, &global.y0wcode, 0);
        if (global.mcu411)
        {
            update_wcode(&global.y1wcode);
            decBlockCPU(global.pix_y + 2*16, &global.y1wcode, 0);
            global.y1wcode.basis_curr = global.y1wcode.basis_next;
            global.y1wcode.dcbuf_curr = global.y1wcode.dcbuf_next;
            global.y1wcode.basis_curr_line++;
            global.y1wcode.dcbuf_curr_line++;
            decBlockCPU(global.pix_y + 3*16, &global.y1wcode, 0);
        }
        global.u_wcode.basis_curr = global.u_wcode.basis_next;
        global.u_wcode.dcbuf_curr = global.u_wcode.dcbuf_next;
        global.u_wcode.basis_curr_line++;
        global.u_wcode.dcbuf_curr_line++;
        decBlockCPU(global.pix_u, &global.u_wcode, 1);
        global.v_wcode.basis_curr = global.v_wcode.basis_next;
        global.v_wcode.dcbuf_curr = global.v_wcode.dcbuf_next;
        global.v_wcode.basis_curr_line++;
        global.v_wcode.dcbuf_curr_line++;
        decBlockCPU(global.pix_v, &global.v_wcode, 2);
        global.ColorConv(outbuf, global.pix_y, global.pix_u, global.pix_v);
    }
}

// HVQM4: kinda like IpicPlaneDec
static void decFrame(u32 *outbuf)
{
    // first line hast prev = curr
    global.u_wcode.basis_prev_line = global.basis[1];
    global.u_wcode.dcbuf_prev_line = global.dcbuf[1];
    global.u_wcode.basis_curr_line = global.basis[1];
    global.u_wcode.dcbuf_curr_line = global.dcbuf[1];
    global.u_wcode.basis_next_line = global.basis[1] + global.col_hblocks;
    global.u_wcode.dcbuf_next_line = global.dcbuf[1] + global.col_hblocks;
    global.v_wcode.basis_prev_line = global.basis[2];
    global.v_wcode.dcbuf_prev_line = global.dcbuf[2];
    global.v_wcode.basis_curr_line = global.basis[2];
    global.v_wcode.dcbuf_curr_line = global.dcbuf[2];
    global.v_wcode.basis_next_line = global.basis[2] + global.col_hblocks;
    global.v_wcode.dcbuf_next_line = global.dcbuf[2] + global.col_hblocks;
    global.y0wcode.basis_prev_line = global.basis[0];
    global.y0wcode.dcbuf_prev_line = global.dcbuf[0];
    global.y0wcode.basis_curr_line = global.basis[0];
    global.y0wcode.dcbuf_curr_line = global.dcbuf[0];
    global.y0wcode.basis_next_line = global.basis[0] + global.lum_hblocks;
    global.y0wcode.dcbuf_next_line = global.dcbuf[0] + global.lum_hblocks;
    if (global.mcu411)
    {
        global.y1wcode.basis_prev_line = global.basis[0];
        global.y1wcode.dcbuf_prev_line = global.dcbuf[0];
        global.y1wcode.basis_curr_line = global.basis[0] + global.lum_hblocks;
        global.y1wcode.dcbuf_curr_line = global.dcbuf[0] + global.lum_hblocks;
        global.y1wcode.basis_next_line = global.basis[0] + global.lum_hblocks + global.lum_hblocks;
        global.y1wcode.dcbuf_next_line = global.dcbuf[0] + global.lum_hblocks + global.lum_hblocks;
    }
    decLine(outbuf);

    // advance
    outbuf += global.mcu_v_pix;
    global.u_wcode.basis_prev_line = global.basis[1];
    global.u_wcode.dcbuf_prev_line = global.dcbuf[1];
    global.v_wcode.basis_prev_line = global.basis[2];
    global.v_wcode.dcbuf_prev_line = global.dcbuf[2];
    if (global.mcu411)
    {
        global.y0wcode.basis_curr_line += global.lum_hblocks;
        global.y0wcode.dcbuf_curr_line += global.lum_hblocks;
        global.y0wcode.basis_next_line += global.lum_hblocks;
        global.y0wcode.dcbuf_next_line += global.lum_hblocks;
        advance_wcode(&global.y1wcode, global.lum_hblocks);
    }
    else
    {
        global.y0wcode.basis_prev_line = global.basis[0];
        global.y0wcode.dcbuf_prev_line = global.dcbuf[0];
    }

    // middle lines
    for (u32 line = 1; line < global.col_vblocks - 1; ++line)
    {
        decLine(outbuf);
        outbuf += global.mcu_v_pix;
        if (global.mcu411)
        {
            advance_wcode(&global.y0wcode, global.lum_hblocks);
            advance_wcode(&global.y1wcode, global.lum_hblocks);
        }
    }

    // last line has next = curr
    global.u_wcode.basis_next_line = global.u_wcode.basis_curr_line;
    global.u_wcode.dcbuf_next_line = global.u_wcode.dcbuf_curr_line;
    global.v_wcode.basis_next_line = global.v_wcode.basis_curr_line;
    global.v_wcode.dcbuf_next_line = global.v_wcode.dcbuf_curr_line;
    if (global.mcu411)
    {
        global.y1wcode.basis_next_line = global.y1wcode.basis_curr_line;
        global.y1wcode.dcbuf_next_line = global.y1wcode.dcbuf_curr_line;
    }
    else
    {
        global.y0wcode.basis_next_line = global.y0wcode.basis_curr_line;
        global.y0wcode.dcbuf_next_line = global.y0wcode.dcbuf_curr_line;
    }
    decLine(outbuf);
}

// done
static void my_hvqm2Init2(u8 alpha)
{
    global.pix_alpha = alpha;
    for (s32 i = -0x100; i < 0x200; ++i)
    {
        s32 j = i < 0 ? 0 : (i < 0x100 ? i : 255);
        global.clipT[0x100 + i] = j;
    }
    global.divT[0] = 0;
    for (u32 i = 1; i < 0x200; ++i)
        global.divT[i] = 0x1000 / i;
}

// done
static u32 my_hvqm2Setup2(HVQM2Header *header, u32 outbufWidth)
{
    global.im_width = header->width;
    global.im_height = header->height;
    if (outbufWidth)
        global.fb_width = outbufWidth;
    else
        global.fb_width = header->width;
    global.h_sampling_rate = header->h_sampling_rate;
    global.v_sampling_rate = header->v_sampling_rate;
    global.mcu_h_pix = header->h_sampling_rate * 4;
    global.mcu_v_pix = header->v_sampling_rate * 4 * global.fb_width;
    global.lum_hblocks = (s32)header->width >> 2;
    global.lum_vblocks = (s32)header->height >> 2;
    global.lum_totalblocks = global.lum_hblocks * global.lum_vblocks;
    global.mcu411 = header->v_sampling_rate == 2;
    global.col_hblocks = (s32)header->width >> 3;
    global.col_vblocks = (s32)header->height >> (global.mcu411 ? 3 : 2);
    global.col_totalblocks = global.col_hblocks * global.col_vblocks;
    global.next_macroblk_line = global.fb_width * 8;
    global.ColorConv = global.mcu411 ? ColorConv411 : ColorConv422;
    global.yshift = header->y_shiftnum;
    if (global.yshift == 8)
    {
        global.nest_w = 70;
        global.nest_h = 38;
    }
    else
    {
        global.nest_w = 38;
        global.nest_h = 70;
    }
    global.dc_min = 0xFF80 << header->video_quantize_shift;
    global.dc_max = 0x007F << header->video_quantize_shift;
    for (int i = 0; i < 0x100; ++i)
    {
        global.scale_tree.array[0][i] = (s8)i << 3;
        global.dcval_tree.array[0][i] = (s8)i << header->video_quantize_shift;
    }
    return header->total_frames;
}

// done, works
static void Ipic_BasisNumDec()
{
    u8 rle = 0;
    for (int i = 0; i < global.lum_totalblocks; ++i)
    {
        if (rle)
        {
            global.basis[0][i] = 0;
            --rle;
        }
        else
        {
            u8 value = decodeHuff(&global.basisnum_buf[0], &global.basisnum_tree);
            if (value == 0)
                rle = decodeHuff(&global.basisnumrun_buf[0], &global.basisnumrun_tree);
            global.basis[0][i] = value;
        }
    }
    rle = 0;
    for (int i = 0; i < global.col_totalblocks; ++i)
    {
        if (rle)
        {
            global.basis[1][i] = 0;
            global.basis[2][i] = 0;
            --rle;
        }
        else
        {
            u8 value = decodeHuff(&global.basisnum_buf[1], &global.basisnum_tree);
            if (value == 0)
                rle = decodeHuff(&global.basisnumrun_buf[1], &global.basisnumrun_tree);
            global.basis[1][i] = value & 0xF;
            global.basis[2][i] = value >> 4;
        }
    }
}

// TODO: verify type of `rle`
static s16 getDeltaDC(int plane_idx, int *rle)
{
    if (*rle == 0)
    {
        s16 delta = decodeDC(&global.dcval_buf[plane_idx]);
        if (delta == 0)
            *rle = decodeHuff(&global.dcrun_buf[plane_idx], &global.basisnumrun_tree);
        return delta;
    }
    else
    {
        --*rle;
        return 0;
    }
}

static u32 mean(u32 a, u32 b)
{
    // HVQM4: (a + b + 1) / 2
    return (a + b) / 2;
}

static void IpicDcvDec()
{
    u8 dcvalY = 0;
    u8 dcvalU = 0;
    u8 dcvalV = 0;
    int rleY = 0;
    int rleU = 0;
    int rleV = 0;
    u8 *dcbufY = global.dcbuf[0];
    u8 *dcbufU = global.dcbuf[1];
    u8 *dcbufV = global.dcbuf[2];
    u8 *dcbufY_prev = global.dcbuf[0];
    u8 *dcbufU_prev = global.dcbuf[1];
    u8 *dcbufV_prev = global.dcbuf[2];
    for (u32 h = 0; h < global.col_hblocks; ++h)
    {
        dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY;
        dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY;
        dcvalU += getDeltaDC(1, &rleU); *dcbufU++ = dcvalU;
        dcvalV += getDeltaDC(2, &rleV); *dcbufV++ = dcvalV;
    }
    dcvalY = *dcbufY_prev++;
    dcvalU = *dcbufU_prev++;
    dcvalV = *dcbufV_prev++;
    if (global.mcu411)
    {
        for (u32 h = 0; h < global.col_hblocks; ++h)
        {
            dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY; dcvalY = mean(dcvalY, *dcbufY_prev++);
            dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY; dcvalY = mean(dcvalY, *dcbufY_prev++);
        }
    }
    for (u32 v = 0; v < global.col_vblocks - 1; ++v)
    {
        dcvalY = *(dcbufY_prev - 1);
        dcvalU = *(dcbufU_prev - 1);
        dcvalV = *(dcbufV_prev - 1);
        for (u32 h = 0; h < global.col_hblocks; ++h)
        {
            dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY; dcvalY = mean(dcvalY, *dcbufY_prev++);
            dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY; dcvalY = mean(dcvalY, *dcbufY_prev++);
            dcvalU += getDeltaDC(1, &rleU); *dcbufU++ = dcvalU; dcvalU = mean(dcvalU, *dcbufU_prev++);
            dcvalV += getDeltaDC(2, &rleV); *dcbufV++ = dcvalV; dcvalV = mean(dcvalV, *dcbufV_prev++);
        }
        dcvalY = *(dcbufY_prev - 1);
        if (global.mcu411)
        {
            for (u32 h = 0; h < global.col_hblocks; ++h)
            {
                dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY; dcvalY = mean(dcvalY, *dcbufY_prev++);
                dcvalY += getDeltaDC(0, &rleY); *dcbufY++ = dcvalY; dcvalY = mean(dcvalY, *dcbufY_prev++);
            }
        }
    }
}

static void MakeNest(u16 nest_start_x, u16 nest_start_y)
{
    // names from HVQM4
    u32 h_nest_blocks, v_nest_blocks, h_mirror, v_mirror, h_empty, v_empty;
    if (global.lum_hblocks < global.nest_w)
    {
        h_nest_blocks = global.lum_hblocks;
        h_mirror = global.lum_hblocks < global.nest_w - global.lum_hblocks ? global.lum_hblocks : global.nest_w - global.lum_hblocks;
        h_empty = global.nest_w - (h_nest_blocks + h_mirror);
    }
    else
    {
        h_nest_blocks = global.nest_w;
        h_mirror = 0;
        h_empty = 0;
    }
    if (global.lum_vblocks < global.nest_h)
    {
        v_nest_blocks = global.lum_vblocks;
        v_mirror = global.lum_vblocks < global.nest_h - global.lum_vblocks ? global.lum_vblocks : global.nest_h - global.lum_vblocks;
        v_empty = global.nest_w - (v_nest_blocks + v_mirror);
    }
    else
    {
        v_nest_blocks = global.nest_h;
        v_mirror = 0;
        v_empty = 0;
    }
    uint8_t *nest = global.nest;
    u8 const *ptr = &global.dcbuf[0][global.lum_hblocks * nest_start_y + nest_start_x];
    for (u32 i = 0; i < v_nest_blocks; ++i)
    {
        u8 const *p = ptr;
        for (u32 j = 0; j < h_nest_blocks; ++j)
            *nest++ = *p++;
        for (u32 j = 0; j < h_mirror; ++j)
            *nest++ = *--p;
        for (u32 j = 0; j < h_empty; ++j)
            *nest++ = 0;
        ptr += global.lum_hblocks;
    }

    uint8_t const *nest2 = nest - global.nest_w;
    for (u32 i = 0; i < v_mirror; ++i)
    {
        for (u32 j = 0; j < global.nest_w; ++j)
            *nest++ = nest2[j];
        nest2 -= global.nest_w;
    }

    for (u32 i = 0; i < v_empty; ++i)
        for (u32 j = 0; j < global.nest_w; ++j)
            *nest++ = 0;
}

static void HVQMDecodeIpic(HVQM2KeyFrame const *keyframe, void const *code)
{
    for (int i = 0; i < 3; ++i)
        setCode(&global.dcrun_buf[i], code + read32(&keyframe->dcrun_offset[i]));
    Ipic_BasisNumDec();
    IpicDcvDec();
    MakeNest(read16(&keyframe->nest_start_x), read16(&keyframe->nest_start_y));
}

static void HVQMDecodePpic(HVQM2PredictFrame const *predict, void const *code, void *outbuf, void const *previm)
{
    setCode(&global.movevector_buf, code + read32(&predict->movevector_offset));
    readTree(&global.movevector_buf, &global.movevector_tree);

    setCode(&global.macroblock_buf, code + read32(&predict->macroblock_offset));

    global.y0wcode.basis_curr_line = global.basis[0];
    global.y0wcode.dcbuf_curr_line = global.dcbuf[0];
    global.u_wcode.basis_curr_line = global.basis[1];
    global.u_wcode.dcbuf_curr_line = global.dcbuf[1];
    global.v_wcode.basis_curr_line = global.basis[2];
    global.v_wcode.dcbuf_curr_line = global.dcbuf[2];

    global.y0wcode.basis_next_line = global.basis[0] + global.lum_hblocks;
    global.y0wcode.dcbuf_next_line = global.dcbuf[0] + global.lum_hblocks;
    if (!global.mcu411)
    {
        global.u_wcode.basis_next_line = global.basis[1] + global.col_hblocks;
        global.u_wcode.dcbuf_next_line = global.dcbuf[1] + global.col_hblocks;
        global.v_wcode.basis_next_line = global.basis[2] + global.col_hblocks;
        global.v_wcode.dcbuf_next_line = global.dcbuf[2] + global.col_hblocks;
    }

    u8 rle_lum = 0;
    u8 rle_col = 0;
    u8 dc_y = 0, dc_u = 0, dc_v = 0;
    s8 mv_x = 0, mv_y = 0;
    void *outbuf_line = outbuf;
    for (u32 y = 0; y < global.im_height; y += 8)
    {
        void *outbuf_pos = outbuf_line;
        for (u32 x = 0; x < global.im_width; x += 8)
        {
            if (getBit(&global.macroblock_buf))
            {
                if (getBit(&global.macroblock_buf))
                {
                    dc_y += decodeDC(&global.dcval_buf[0]);
                    dc_u += decodeDC(&global.dcval_buf[1]);
                    dc_v += decodeDC(&global.dcval_buf[2]);
                    *global.y0wcode.basis_curr_line++ = 0;
                    *global.y0wcode.basis_curr_line++ = 0;
                    *global.y0wcode.dcbuf_curr_line++ = dc_y;
                    *global.y0wcode.dcbuf_curr_line++ = dc_y;
                    *global.y0wcode.basis_next_line++ = 0;
                    *global.y0wcode.basis_next_line++ = 0;
                    *global.y0wcode.dcbuf_next_line++ = dc_y;
                    *global.y0wcode.dcbuf_next_line++ = dc_y;
                    *global.u_wcode.basis_curr_line++ = 0;
                    *global.u_wcode.dcbuf_curr_line++ = dc_u;
                    *global.v_wcode.basis_curr_line++ = 0;
                    *global.v_wcode.dcbuf_curr_line++ = dc_v;
                    if (!global.mcu411)
                    {
                        *global.u_wcode.basis_next_line++ = 0;
                        *global.u_wcode.dcbuf_next_line++ = dc_u;
                        *global.v_wcode.basis_next_line++ = 0;
                        *global.v_wcode.dcbuf_next_line++ = dc_v;
                    }
                }
                else
                {
                    *global.y0wcode.basis_curr_line++ = getDeltaBN(&rle_lum, &global.basisnum_buf[0], &global.basisnumrun_buf[0]);
                    *global.y0wcode.basis_curr_line++ = getDeltaBN(&rle_lum, &global.basisnum_buf[0], &global.basisnumrun_buf[0]);
                    *global.y0wcode.dcbuf_curr_line++ = (dc_y += decodeDC(&global.dcval_buf[0]));
                    *global.y0wcode.dcbuf_curr_line++ = (dc_y += decodeDC(&global.dcval_buf[0]));
                    *global.y0wcode.basis_next_line++ = getDeltaBN(&rle_lum, &global.basisnum_buf[0], &global.basisnumrun_buf[0]);
                    *global.y0wcode.basis_next_line++ = getDeltaBN(&rle_lum, &global.basisnum_buf[0], &global.basisnumrun_buf[0]);
                    *global.y0wcode.dcbuf_next_line++ = (dc_y += decodeDC(&global.dcval_buf[0]));
                    *global.y0wcode.dcbuf_next_line++ = (dc_y += decodeDC(&global.dcval_buf[0]));
                    u8 deltaBN = getDeltaBN(&rle_col, &global.basisnum_buf[1], &global.basisnumrun_buf[1]);
                    *global.u_wcode.basis_curr_line++ = (deltaBN >> 0) & 0xF;
                    *global.v_wcode.basis_curr_line++ = (deltaBN >> 4) & 0xF;
                    *global.u_wcode.dcbuf_curr_line++ = (dc_u += decodeDC(&global.dcval_buf[1]));
                    *global.v_wcode.dcbuf_curr_line++ = (dc_v += decodeDC(&global.dcval_buf[2]));
                    if (!global.mcu411)
                    {
                        deltaBN = getDeltaBN(&rle_col, &global.basisnum_buf[1], &global.basisnumrun_buf[1]);
                        *global.u_wcode.basis_next_line++ = (deltaBN >> 0) & 0xF;
                        *global.v_wcode.basis_next_line++ = (deltaBN >> 4) & 0xF;
                        *global.u_wcode.dcbuf_next_line++ = (dc_u += decodeDC(&global.dcval_buf[1]));
                        *global.v_wcode.dcbuf_next_line++ = (dc_v += decodeDC(&global.dcval_buf[2]));
                    }
                }
            }
            else
            {
                mv_x += decodeHuff(&global.movevector_buf, &global.movevector_tree);
                mv_y += decodeHuff(&global.movevector_buf, &global.movevector_tree);

                u32 pos_y = y + mv_y;
                u32 pos_x = x + mv_x;
                u32 const *src = previm + pos_x + pos_y * global.fb_width;
                u32 *dst = outbuf_pos;
                for (u32 i = 0; i < 8; ++i)
                {
                    for (u32 j = 0; j < 8; ++j)
                        dst[j] = src[j];
                    src += global.fb_width;
                    dst += global.fb_width;
                }

                *global.y0wcode.basis_curr_line++ = 0x80;
                *global.y0wcode.basis_curr_line++ = 0x80;
                global.y0wcode.dcbuf_curr_line += 2;
                *global.y0wcode.basis_next_line++ = 0x80;
                *global.y0wcode.basis_next_line++ = 0x80;
                global.y0wcode.dcbuf_next_line += 2;

                *global.u_wcode.basis_curr_line++ = 0x80;
                global.u_wcode.dcbuf_curr_line++;
                *global.v_wcode.basis_curr_line++ = 0x80;
                global.v_wcode.dcbuf_curr_line++;
                if (!global.mcu411)
                {
                    *global.u_wcode.basis_next_line++ = 0x80;
                    global.u_wcode.dcbuf_next_line++;
                    *global.v_wcode.basis_next_line++ = 0x80;
                    global.v_wcode.dcbuf_next_line++;
                }
            }
            outbuf_pos += 32;
        }
        global.y0wcode.basis_curr_line += global.lum_hblocks;
        global.y0wcode.dcbuf_curr_line += global.lum_hblocks;
        global.y0wcode.basis_next_line += global.lum_hblocks;
        global.y0wcode.dcbuf_next_line += global.lum_hblocks;
        if (!global.mcu411)
        {
            global.u_wcode.basis_curr_line += global.col_hblocks;
            global.u_wcode.dcbuf_curr_line += global.col_hblocks;
            global.u_wcode.basis_next_line += global.col_hblocks;
            global.u_wcode.dcbuf_next_line += global.col_hblocks;
            global.v_wcode.basis_curr_line += global.col_hblocks;
            global.v_wcode.dcbuf_curr_line += global.col_hblocks;
            global.v_wcode.basis_next_line += global.col_hblocks;
            global.v_wcode.dcbuf_next_line += global.col_hblocks;
        }
        outbuf_line += global.next_macroblk_line;
    }
}

static void HVQMDecodeHpic(u32 *outbuf, u32 const *previm)
{
    for (u32 i = 0; i < global.im_height; ++i)
    {
        for (u32 j = 0; j < global.im_width; ++j)
            outbuf[j] = previm[j];
        outbuf += global.fb_width;
        previm += global.fb_width;
    }
}

static void my_hvqm2Decode2(void const *code, u32 format, u32 *outbuf, u32 *previm, void *workbuf)
{
    if (format == HVQM2_VIDEO_HOLD)
        return HVQMDecodeHpic(outbuf, previm);
    global.basis[0] = workbuf; workbuf += global.lum_totalblocks;
    global.dcbuf[0] = workbuf; workbuf += global.lum_totalblocks;
    global.basis[1] = workbuf; workbuf += global.col_totalblocks;
    global.dcbuf[1] = workbuf; workbuf += global.col_totalblocks;
    global.basis[2] = workbuf; workbuf += global.col_totalblocks;
    global.dcbuf[2] = workbuf;
    HVQM2Frame const *frame = code;
    for (int i = 0; i < 2; ++i)
    {
        setCode(&global.basisnum_buf[i], code + read32(&frame->basisnum_offset[i]));
        setCode(&global.basisnumrun_buf[i], code + read32(&frame->basnumrn_offset[i]));
    }
    for (int i = 0; i < 3; ++i)
    {
        setCode(&global.scale_buf[i], code + read32(&frame->scale_offset[i]));
        setCode(&global.dcval_buf[i], code + read32(&frame->dcval_offset[i]));
        global.fixvl[i] = code + read32(&frame->fixvl_offset[i]);
    }
    readTree(&global.basisnum_buf[0], &global.basisnum_tree);
    readTree(&global.basisnumrun_buf[0], &global.basisnumrun_tree);
    readTree(&global.scale_buf[0], &global.scale_tree);
    readTree(&global.dcval_buf[0], &global.dcval_tree);
    global.nest = &global.nestbuf[0];

    if (format == HVQM2_VIDEO_KEYFRAME)
    {
        HVQM2KeyFrame const *keyframe = code + sizeof(HVQM2Frame);
        HVQMDecodeIpic(keyframe, code);
        decFrame(outbuf);
    }
    else
    {
        HVQM2PredictFrame const *predict = code + sizeof(HVQM2Frame);
        HVQMDecodePpic(predict, code, outbuf, previm);
        //decFrame(outbuf);
    }
}

static void dumpRGB(HVQM2Header const *header, u32 frame, u32 frame_type, u8 const *outbuf)
{
    char path[128];
    snprintf(path, sizeof(path), OUTPUT_DIR "/video_%i_%c.ppm", frame, "iph"[frame_type]);
    FILE *outfile = fopen(path, "wb+");
    fprintf(outfile, "P6\n%u %u\n255\n", header->width, header->height);
    u8 temp[MAXWIDTH*MAXHEIGHT*3];
    // remove alpha
    for (u32 i = 0; i < header->width * header->height; ++i)
    {
        temp[3*i+0] = outbuf[4*i+0];
        temp[3*i+1] = outbuf[4*i+1];
        temp[3*i+2] = outbuf[4*i+2];
    }
    fwrite(temp, header->width * header->height * 3, 1, outfile);
    fclose(outfile);
}

int main(int argc, char **argv)
{
    HVQM2Header header;
    HVQM2Record record;
    u8 buffer[HVQM_DATASIZE_MAX];
    _Static_assert(sizeof(HVQM2Header) <= HVQM_DATASIZE_MAX, "");
    _Static_assert(sizeof(HVQM2Record) <= HVQM_DATASIZE_MAX, "");
    u8 outbuf[2][MAXWIDTH*MAXHEIGHT*4];
    u16 workbuf[(MAXWIDTH/8)*(MAXHEIGHT/4)*4];

    if (argc != 2)
    {
        fprintf(stderr, "usage: %s <video.hvqm>\n", argv[0]);
        return -1;
    }
    FILE *f = fopen(argv[1], "rb");
    if (!f)
    {
        fprintf(stderr, "error while opening input file\n");
        return -1;
    };
    fread(buffer, 1, sizeof(HVQM2Header), f);
    memcpy(header.file_version, buffer, 16);
    header.file_size             = read32(&buffer[16]);
    header.width                 = read16(&buffer[20]);
    header.height                = read16(&buffer[22]);
    header.h_sampling_rate       = buffer[24];
    header.v_sampling_rate       = buffer[25];
    header.y_shiftnum            = buffer[26];
    header.video_quantize_shift  = buffer[27];
    header.total_frames          = read32(&buffer[28]);
    header.usec_per_frame        = read32(&buffer[32]);
    header.max_frame_size        = read32(&buffer[36]);
    header.max_sp_packets        = read32(&buffer[40]);
    header.audio_format          = buffer[44];
    header.channels              = buffer[45];
    header.sample_bits           = buffer[46];
    header.audio_quantize_step   = buffer[47];
    header.total_audio_records   = read32(&buffer[48]);
    header.samples_per_sec       = read32(&buffer[52]);
    header.max_audio_record_size = read32(&buffer[56]);
    if (memcmp(header.file_version, "HVQM2 1.0\0\0\0\0\0\0", sizeof(header.file_version)) != 0)
    {
        fprintf(stderr, "that doesn't look like an HVQM2 1.0 file\n");
        return -1;
    }
    if (header.width > MAXWIDTH || header.height > MAXHEIGHT)
    {
        fprintf(stderr, "video too large\n");
        return -1;
    }

    u8 alpha = 0xFF;
    u32 outbufWidth = header.width;
    my_hvqm2Init2(alpha);
    my_hvqm2Setup2(&header, outbufWidth);
#ifndef NATIVE
    hvqm2Init2(alpha);
    hvqm2Setup2(&header, outbufWidth);
#endif

    u32 curr = 0;
    for (u32 frame = 0; frame < header.total_frames; ++frame)
    {
        fread(buffer, 1, sizeof(HVQM2Record), f);
        record.type   = read16(&buffer[0]);
        record.format = read16(&buffer[2]);
        record.size   = read32(&buffer[4]);
        if (record.size > HVQM_DATASIZE_MAX)
        {
            fprintf(stderr, "record too large\n");
            return -1;
        }
        fread(buffer, 1, record.size, f);
        if (record.type != HVQM2_VIDEO)
            continue;
        switch (record.format)
        {
            case HVQM2_VIDEO_KEYFRAME: putchar('I'); break;
            case HVQM2_VIDEO_PREDICT:  putchar('P'); break;
            case HVQM2_VIDEO_HOLD:     putchar('H'); break;
            default:
                fprintf(stderr, "invalid frame format\n");
                exit(-1);
        }
        fflush(stdout);
#ifdef NATIVE
        my_hvqm2Decode2(buffer, record.format, (u32*)outbuf[curr], (u32*)outbuf[1 - curr], workbuf);
#else
        hvqm2Decode2(buffer, record.format, (u32*)outbuf[curr], (u32*)outbuf[1 - curr], workbuf);
#endif

        if (record.format == HVQM2_VIDEO_KEYFRAME)
        {
            u8 *wb = (u8*)workbuf;
            global.basis[0] = wb; wb += global.lum_totalblocks;
            global.dcbuf[0] = wb; wb += global.lum_totalblocks;
            global.basis[1] = wb; wb += global.col_totalblocks;
            global.dcbuf[1] = wb; wb += global.col_totalblocks;
            global.basis[2] = wb; wb += global.col_totalblocks;
            global.dcbuf[2] = wb;
            dumpPlanes(frame, record.format);
            dumpRGB(&header, frame, record.format, outbuf[curr]);
        }

        curr ^= 1;
    }
    puts("");
}
