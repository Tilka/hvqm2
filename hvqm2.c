#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hvqm2dec.h>

#define MAXWIDTH 640
#define MAXHEIGHT 480
#define HVQM_DATASIZE_MAX 30000

#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-parameter"

typedef struct BitBuffer
{
    u32 bit;
    u32 value;
    u32 const *pos;
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
    u16 root;
    u16 array[2][0x100];
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
    u16 dc_max;
    u16 dc_min;
    u8 *basisY;
    u8 *basisU;
    u8 *basisV;
    void *dcbufY;
    void *dcbufU;
    void *dcbufV;
    // some data that sgi calls wcode
    u32 divT[0x200];
    u8 clipT[0x300];
    void (*ColorConv)();
    s32 outbufWidth;
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

static u32 getBit(BitBuffer *buf)
{
    u32 bit = buf->bit;
    if (bit == 0)
    {
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
    return tree->array[0][pos];
}

// done
// HVQM4: kinda like decodeSOvfSym
static s16 decodeDC(BitBuffer *buf)
{
    s16 sum = 0;
    s16 value;
    do
    {
        value = decodeHuff(buf, &global.dcval_tree);
        sum += value;
    } while (value == global.dc_min || value == global.dc_max);
    return sum;
}

// done
// HVQM4: kinda like getDeltaDC
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
static void decBlockCPU(u32 a0, u8 const *ptr, u32 a2)
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
static void decLine()
{
    // TODO
}

// HVQM4: kinda like IpicPlaneDec
static void decFrame()
{
    // TODO
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
        global.outbufWidth = outbufWidth;
    else
        global.outbufWidth = header->width;
    global.h_sampling_rate = header->h_sampling_rate;
    global.v_sampling_rate = header->v_sampling_rate;
    global.mcu_h_pix = header->h_sampling_rate * 4;
    global.mcu_v_pix = header->v_sampling_rate * 4 * outbufWidth;
    global.lum_hblocks = header->width >> 2;
    global.lum_vblocks = header->height >> 2;
    global.lum_totalblocks = global.lum_hblocks * global.lum_vblocks;
    global.col_hblocks = header->width >> 3;
    global.col_vblocks = header->height >> 3;
    global.col_totalblocks = global.col_hblocks * global.col_vblocks;
    global.next_macroblk_line = global.outbufWidth * 8;
    global.mcu411 = (header->v_sampling_rate ^ 2) < 1;
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

// done
static void Ipic_BasisNumDec()
{
    u32 rle = 0;
    for (int i = 0; i < global.lum_totalblocks; ++i)
    {
        if (rle)
        {
            global.basisY[i] = 0;
            --rle;
        }
        else
        {
            u8 value = decodeHuff(&global.basisnum_buf[0], &global.basisnum_tree);
            if (value == 0)
                rle = decodeHuff(&global.basisnumrun_buf[0], &global.basisnumrun_tree);
            global.basisY[i] = value;
        }
    }
    rle = 0;
    for (int i = 0; i < global.col_totalblocks; ++i)
    {
        if (rle)
        {
            global.basisU[i] = 0;
            global.basisV[i] = 0;
            --rle;
        }
        else
        {
            u8 value = decodeHuff(&global.basisnum_buf[1], &global.basisnum_tree);
            if (value == 0)
                rle = decodeHuff(&global.basisnumrun_buf[1], &global.basisnumrun_tree);
            global.basisU[i] = value & 0xF;
            global.basisV[i] = value >> 4;
        }
    }
}

static void IpicDcvDec()
{
    // TODO
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
        outbuf += global.outbufWidth;
        previm += global.outbufWidth;
    }
}

static void my_hvqm2Decode2(void const *code, u32 format, u32 *outbuf, u32 *previm, void *workbuf_)
{
    if (format == HVQM2_VIDEO_HOLD)
        return HVQMDecodeHpic(outbuf, previm);
    u8 *workbuf = workbuf_;
    global.basisY = workbuf; workbuf += global.lum_totalblocks;
    global.dcbufY = workbuf; workbuf += global.lum_totalblocks;
    global.basisU = workbuf; workbuf += global.col_totalblocks;
    global.dcbufU = workbuf; workbuf += global.col_totalblocks;
    global.basisV = workbuf; workbuf += global.col_totalblocks;
    global.dcbufV = workbuf;
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
    // TODO
}

int main(int argc, char **argv)
{
    HVQM2Header header;
    HVQM2Record record;
    u8 buffer[HVQM_DATASIZE_MAX];
    _Static_assert(sizeof(HVQM2Header) < HVQM_DATASIZE_MAX, "");
    _Static_assert(sizeof(HVQM2Record) < HVQM_DATASIZE_MAX, "");
    u8 outbuf[2][MAXWIDTH*MAXHEIGHT*4];
    u16 workbuf[(MAXWIDTH/8)*(MAXHEIGHT/4)*4];

    if (argc != 2)
    {
        fprintf(stderr, "usage: %s <video.hvqm>", argv[0]);
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
#ifdef NATIVE
    my_hvqm2Init2(alpha);
    my_hvqm2Setup2(&header, outbufWidth);
#else
    hvqm2Init2(alpha);
    hvqm2Setup2(&header, outbufWidth);
#endif

    u32 curr = 0;
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
            case HVQM2_VIDEO_KEYFRAME: putchar('I'); break;
            case HVQM2_VIDEO_PREDICT:  putchar('P'); break;
            case HVQM2_VIDEO_HOLD:     putchar('H'); break;
            default:
                fprintf(stderr, "invalid frame format");
                exit(-1);
        }
        fflush(stdout);
#ifdef NATIVE
        my_hvqm2Decode2(buffer, record.format, (u32*)outbuf[curr], (u32*)outbuf[1 - curr], workbuf);
#else
        hvqm2Decode2(buffer, record.format, (u32*)outbuf[curr], (u32*)outbuf[1 - curr], workbuf);
#endif

        char path[128];
        snprintf(path, sizeof(path), "output/video_rgb_%i.ppm", frame);
        FILE *outfile = fopen(path, "wb+");
        fprintf(outfile, "P6\n%u %u\n255\n", header.width, header.height);
        u8 temp[MAXWIDTH*MAXHEIGHT*3];
        // remove alpha
        for (u32 i = 0; i < header.width * header.height; ++i)
        {
            temp[3*i+0] = outbuf[curr][4*i+0];
            temp[3*i+1] = outbuf[curr][4*i+1];
            temp[3*i+2] = outbuf[curr][4*i+2];
        }
        fwrite(temp, header.width * header.height * 3, 1, outfile);
        fclose(outfile);

        curr ^= 1;
        ++frame;
    }
    puts("");
}
