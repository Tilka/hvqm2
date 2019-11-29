#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hvqm2dec.h>

#define MAXWIDTH 640
#define MAXHEIGHT 480
#define HVQM_DATASIZE_MAX 30000

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
    u32 fixvl_arr[3];
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
        u8 const *basis_prev;
        u8 const *dcbuf_prev;
        u8 const *basis_curr;
        u8 const *dcbuf_curr;
        u8 const *basis_next;
        u8 const *dcbuf_next;
        u8 unk[5];
    } u_wcode, v_wcode, y0wcode, y1wcode;
    u32 divT[0x200];
    u8 clipT[0x300];
    void (*ColorConv)();
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
    //
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
static void dumpPlanes(char const *dir)
{
    static u32 frame = 0;
    ++frame;
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        char path[128];
        snprintf(path, 128, "%s/basis_%u_%c.ppm", dir, frame, "yuv"[plane_idx]);
        dumpPlane(path, global.basis[plane_idx], plane_idx);
        snprintf(path, 128, "%s/dcbuf_%u_%c.ppm", dir, frame, "yuv"[plane_idx]);
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
            fprintf(stderr, "error: BitBuffer overread\n");
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
// only used by P frames
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

// HVQM4: kinda like IpicBlockDec
static void decBlockCPU()
{
    // TODO
}

static void ColorConv422()
{
    // TODO
}

static void ColorConv411()
{
    // TODO
}

// HVQM4: kinda like IpicLineDec
static void decLine(u32 *outbuf)
{
    // TODO
    decBlockCPU();
}

static void advance_y_wcode(u32 amount)
{
    global.y0wcode.basis_prev += amount;
    global.y0wcode.dcbuf_prev += amount;
    global.y0wcode.basis_curr += amount;
    global.y0wcode.dcbuf_curr += amount;
    global.y0wcode.basis_next += amount;
    global.y0wcode.dcbuf_next += amount;
    global.y1wcode.basis_prev += amount;
    global.y1wcode.dcbuf_prev += amount;
    global.y1wcode.basis_curr += amount;
    global.y1wcode.dcbuf_curr += amount;
    global.y1wcode.basis_next += amount;
    global.y1wcode.dcbuf_next += amount;
}

// HVQM4: kinda like IpicPlaneDec
static void decFrame(u32 *outbuf)
{
    // first line
    global.u_wcode.basis_prev = global.basis[1];
    global.u_wcode.dcbuf_prev = global.dcbuf[1];
    global.u_wcode.basis_curr = global.basis[1];
    global.u_wcode.dcbuf_curr = global.dcbuf[1];
    global.u_wcode.basis_next = global.basis[1] + global.col_hblocks;
    global.u_wcode.dcbuf_next = global.dcbuf[1] + global.col_hblocks;
    global.v_wcode.basis_prev = global.basis[2];
    global.v_wcode.dcbuf_prev = global.dcbuf[2];
    global.v_wcode.basis_curr = global.basis[2];
    global.v_wcode.dcbuf_curr = global.dcbuf[2];
    global.v_wcode.basis_next = global.basis[2] + global.col_hblocks;
    global.v_wcode.dcbuf_next = global.dcbuf[2] + global.col_hblocks;
    global.y0wcode.basis_prev = global.basis[0];
    global.y0wcode.dcbuf_prev = global.dcbuf[0];
    global.y0wcode.basis_curr = global.basis[0];
    global.y0wcode.dcbuf_curr = global.dcbuf[0];
    global.y0wcode.basis_next = global.basis[0] + global.lum_hblocks;
    global.y0wcode.dcbuf_next = global.dcbuf[0] + global.lum_hblocks;
    if (global.mcu411)
    {
        global.y1wcode.basis_prev = global.basis[0];
        global.y1wcode.dcbuf_prev = global.dcbuf[0];
        global.y1wcode.basis_curr = global.basis[0] + global.lum_hblocks;
        global.y1wcode.dcbuf_curr = global.dcbuf[0] + global.lum_hblocks;
        global.y1wcode.basis_next = global.basis[0] + global.lum_hblocks + global.lum_hblocks;
        global.y1wcode.dcbuf_next = global.dcbuf[0] + global.lum_hblocks + global.lum_hblocks;
    }
    decLine(outbuf);

    // advance
    outbuf += global.mcu_v_pix;
    global.u_wcode.basis_prev = global.basis[1];
    global.u_wcode.dcbuf_prev = global.dcbuf[1];
    global.v_wcode.basis_prev = global.basis[2];
    global.v_wcode.dcbuf_prev = global.dcbuf[2];
    if (global.mcu411)
    {
        advance_y_wcode(global.lum_hblocks);
    }
    else
    {
        global.y0wcode.basis_prev = global.basis[0];
        global.y0wcode.dcbuf_prev = global.dcbuf[0];
    }

    // middle lines
    for (u32 line = 1; line < global.col_vblocks - 1; ++line)
    {
        decLine(outbuf);
        outbuf += global.mcu_v_pix;
        if (global.mcu411)
        {
            advance_y_wcode(global.lum_hblocks);
        }
    }

    // advance
    global.u_wcode.basis_next = global.u_wcode.basis_curr;
    global.u_wcode.dcbuf_next = global.u_wcode.dcbuf_curr;
    global.v_wcode.basis_next = global.v_wcode.basis_curr;
    global.v_wcode.dcbuf_next = global.v_wcode.dcbuf_curr;
    if (global.mcu411)
    {
        global.y1wcode.basis_next = global.y1wcode.basis_curr;
        global.y1wcode.dcbuf_next = global.y1wcode.dcbuf_curr;
    }
    else
    {
        global.y0wcode.basis_next = global.y0wcode.basis_curr;
        global.y0wcode.dcbuf_next = global.y0wcode.dcbuf_curr;
    }

    // last line
    decLine(outbuf);
}

// done
static void my_hvqm2Init2(u8 alpha)
{
    global.pix_alpha = alpha;
    for (s32 i = -0x100; i < 0x100; ++i)
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
    global.lum_hblocks = header->width / 4;
    global.lum_vblocks = header->height / 4;
    global.lum_totalblocks = global.lum_hblocks * global.lum_vblocks;
    global.mcu411 = header->v_sampling_rate == 2;
    global.col_hblocks = header->width / 8;
    global.col_vblocks = header->height / (global.mcu411 ? 8 : 4);
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
    printf("img: %ux%u\n", global.im_width, global.im_height);
    printf("lum: %ux%u = %u\n", global.lum_hblocks, global.lum_vblocks, global.lum_totalblocks);
    printf("col: %ux%u = %u\n", global.col_hblocks, global.col_vblocks, global.col_totalblocks);
    printf("mcu411: %u\n", global.mcu411);
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

static void MakeNest()
{
    // TODO
}

static void HVQMDecodeIpic(HVQM2KeyFrame const *keyframe, void const *code)
{
    for (int i = 0; i < 3; ++i)
        setCode(&global.dcrun_buf[i], code + read32(&keyframe->dcrun_offset[i]));
    Ipic_BasisNumDec();
    IpicDcvDec();
    MakeNest();
}

static void HVQMDecodePpic(HVQM2PredictFrame const *predict, void const *code)
{
    // TODO
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
        global.fixvl_arr[i] = read32(code + read32(&frame->fixvl_offset[i]));
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
    }
    else
    {
        HVQM2PredictFrame const *predict = code + sizeof(HVQM2Frame);
        HVQMDecodePpic(predict, code);
    }
    decFrame(outbuf);
}

static void dumpRGB(HVQM2Header const *header, int frame, u8 const *outbuf)
{
    char path[128];
    snprintf(path, sizeof(path), "output/video_rgb_%i.ppm", frame);
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
    u32 i_frame = 0;
    for (u32 frame = 0; frame < header.total_frames;)
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
            case HVQM2_VIDEO_KEYFRAME: putchar('I'); ++i_frame; break;
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
            dumpPlanes("/tmp/output");
            dumpRGB(&header, i_frame, outbuf[curr]);
        }

        curr ^= 1;
        ++frame;
    }
    puts("");
}
