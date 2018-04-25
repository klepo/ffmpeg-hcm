/*
 * HCM encoder/decoder
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "libavutil/samplefmt.h"
#include "libavutil/internal.h"
#include "libavutil/opt.h"
#include "avcodec.h"
#include "internal.h"

#include <stdio.h>
#include <complex.h>

typedef _Fcomplex floatC;

#define COMPLEX_R(r) _FCOMPLEX_(r, 0.0f)

typedef struct HCMContext {
    AVClass *class;

    int period;
    int harmonics;
    int mos;

    int oSize;
    int bSize;
    //uint64_t stride;

    float *b;
    floatC *e;
    floatC *bE;
    floatC *bE_1;

    floatC *sCTmp1;
    floatC *sCTmp2;

    uint64_t frame;
    uint64_t frameOffset;
    uint64_t step;
    uint64_t stepOffset;

    int lastFrameDecoded;
} HCMContext;

floatC caddf(floatC n1, floatC n2)
{
    return _FCOMPLEX_(crealf(n1) + crealf(n2), cimagf(n1) + cimagf(n2));
}

/*floatC csubf(floatC n1, floatC n2)
{
    return _FCOMPLEX_(crealf(n1) - crealf(n2), cimagf(n1) - cimagf(n2));
}*/

floatC cmulf(floatC n1, floatC n2)
{
    return _FCOMPLEX_(crealf(n1) * crealf(n2) - cimagf(n1) * cimagf(n2), cimagf(n1) * crealf(n2) + crealf(n1) * cimagf(n2));
}

static void generateE(int period, int ih, int h, int bSize, floatC *e)
{
    floatC i = _FCOMPLEX_(0.0f, -1.0f);
    for (int x = 0; x < bSize; x++) {
        int hx = ih * bSize + x;
        e[hx] = cexpf(cmulf(i, COMPLEX_R((2.0f * (float)(M_PI) / ((float)(period) / (h))) * (float)(x))));
    }
}

static void generateBE(int ih, int bSize, int oSize, const float *b, const floatC *e, floatC *bE, floatC *bE_1, int normalize)
{
    for (int x = 0; x < bSize; x++) {
        int hx = ih * bSize + x;
        bE[hx] = cmulf(COMPLEX_R(b[x]), e[hx]);
        bE_1[hx] = cmulf(COMPLEX_R(b[(x + oSize) % (bSize - 1)]), e[ih * bSize + ((x + oSize) % (bSize - 1))]);
        if (normalize) {
            bE[hx] = cmulf(bE[hx], COMPLEX_R(2.0f / (float)(oSize)));
            bE_1[hx] = cmulf(bE_1[hx], COMPLEX_R(2.0f / (float)(oSize)));
        }
    }
}

static void triangular(int oSize, float *w)
{
    for (int x = 0; x < oSize; x++) {
        w[x] = (float)(x) / oSize;
    }
    for (int x = oSize; x < 2 * oSize + 1; x++) {
        w[x] = 2.0f - (float)(x) / oSize;
    }
}

/*static void hann(int oSize, float *w)
{
    for (int x = 0; x < 2 * oSize + 1; x++) {
        w[x] = (float)(pow(sin(M_PI * x / (2 * oSize)), 2));
    }
}*/

static void generateFunctions(int bSize, int oSize, int period, int harmonics, float *b, floatC *e, floatC *bE, floatC *bE_1, int normalize)
{
    // Generate basis function (window)
    triangular(oSize, b);  // Triangular window
    //hann(oSize, b);        // Hann window

    // Generate complex exponential functions
    for (int h = 0; h < harmonics; h++) {
        generateE(period, h, h + 1, bSize, e);
        generateBE(h, bSize, oSize, b, e, bE, bE_1, normalize);
    }
}

static void initHCMContext(HCMContext *hCMContext, int normalize)
{
    hCMContext->oSize = hCMContext->period * hCMContext->mos;
    hCMContext->bSize = hCMContext->oSize * 2 + 1;
    //hCMContext->stride = hCMContext->harmonics * 2;
    hCMContext->b = (float *)malloc(sizeof(float) * (size_t)hCMContext->bSize);
    hCMContext->e = (floatC *)malloc(sizeof(floatC) * (size_t)(hCMContext->harmonics * hCMContext->bSize));
    hCMContext->bE = (floatC *)malloc(sizeof(floatC) * (size_t)(hCMContext->harmonics * hCMContext->bSize));
    hCMContext->bE_1 = (floatC *)malloc(sizeof(floatC) * (size_t)(hCMContext->harmonics * hCMContext->bSize));

    generateFunctions(hCMContext->bSize, hCMContext->oSize, hCMContext->period, hCMContext->harmonics, hCMContext->b, hCMContext->e, hCMContext->bE, hCMContext->bE_1, normalize);

    hCMContext->sCTmp1 = (floatC *)malloc(sizeof(floatC) * (size_t)hCMContext->harmonics);
    hCMContext->sCTmp2 = (floatC *)malloc(sizeof(floatC) * (size_t)hCMContext->harmonics);
    // Set zeros
    for (int ih = 0; ih < hCMContext->harmonics; ih++) {
        hCMContext->sCTmp1[ih] = COMPLEX_R(0.0f);
        hCMContext->sCTmp2[ih] = COMPLEX_R(0.0f);
    }

    hCMContext->lastFrameDecoded = 0;
}

static void freeHCMContext(HCMContext *hCMContext)
{
    free(hCMContext->b);
    free(hCMContext->e);
    free(hCMContext->bE);
    free(hCMContext->bE_1);
    free(hCMContext->sCTmp2);
    free(hCMContext->sCTmp1);
}

static av_cold int hcm_encode_init(AVCodecContext *avctx)
{
    HCMContext *hCMContext = avctx->priv_data;

    initHCMContext(hCMContext, 1);

    hCMContext->frame = 0;
    hCMContext->stepOffset = 0;

    avctx->sample_rate = 48000;
    avctx->channels = 1;
    avctx->sample_fmt = avctx->codec->sample_fmts[0];
    //avctx->frame_size = 10;
    avctx->bits_per_coded_sample = 32;
    avctx->block_align = 2 * hCMContext->harmonics * avctx->channels * avctx->bits_per_coded_sample / 8;
    avctx->bit_rate = avctx->block_align * 8LL * avctx->sample_rate;
    avctx->channel_layout = AV_CH_LAYOUT_MONO;

    av_log(NULL, AV_LOG_INFO, "avctx->sample_rate: %d\n", avctx->sample_rate);
    av_log(NULL, AV_LOG_INFO, "avctx->channels: %d\n", avctx->channels);
    av_log(NULL, AV_LOG_INFO, "avctx->sample_fmt: %d\n", avctx->sample_fmt);
    av_log(NULL, AV_LOG_INFO, "avctx->frame_size: %d\n", avctx->frame_size);
    av_log(NULL, AV_LOG_INFO, "avctx->bits_per_coded_sample: %d\n", avctx->bits_per_coded_sample);
    av_log(NULL, AV_LOG_INFO, "avctx->block_align: %d\n", avctx->block_align);
    av_log(NULL, AV_LOG_INFO, "avctx->bit_rate: %d\n", avctx->bit_rate);
    av_log(NULL, AV_LOG_INFO, "avctx->channel_layout: %d\n", avctx->channel_layout);

    return 0;
}

static int hcm_encode_frame(AVCodecContext *avctx, AVPacket *avpkt, const AVFrame *frame, int *got_packet_ptr)
{
    HCMContext *hCMContext = avctx->priv_data;

    // avctx            codec context
    // avpkt            output AVPacket (may contain a user-provided buffer)
    // frame            input AVFrame containing the raw data to be encoded
    // got_packet_ptr   encoder sets to 0 or 1 to indicate that a non-empty packet was returned in avpkt.
    // return 0 on success, negative error code on failure

    int steps = frame->nb_samples;
    int maxOutputFrames = (int)(floorf((float)(steps) / hCMContext->oSize));
    int64_t minOutputBytes = hCMContext->harmonics * (int64_t)sizeof(floatC);
    int64_t maxOutputBytes = maxOutputFrames * minOutputBytes;

    int ret;
    if ((ret = ff_alloc_packet2(avctx, avpkt, maxOutputBytes, minOutputBytes)) < 0)
        return ret;

    float *dst = (float *)avpkt->data;
    const float *src = (const float *)frame->data[0];

    int frames = 0;
    // For every step (input AVFrame)
    for (int step = 0; step < steps; step++) {
        // Compute local index
        int stepLocal = (int)((hCMContext->stepOffset + (uint64_t)step) % (uint64_t)(hCMContext->bSize - 1));
        int savingFlag = ((stepLocal + 1) % hCMContext->oSize == 0) ? 1 : 0;
        int oddFrameFlag = ((hCMContext->frame + 1) % 2 == 0) ? 1 : 0;

        //For every harmonics
        for (int ih = 0; ih < hCMContext->harmonics; ih++) {
            // Correlation step
            hCMContext->sCTmp1[ih] = caddf(hCMContext->sCTmp1[ih], cmulf(hCMContext->bE[ih * hCMContext->bSize + stepLocal], COMPLEX_R(src[step])));
            hCMContext->sCTmp2[ih] = caddf(hCMContext->sCTmp2[ih], cmulf(hCMContext->bE_1[ih * hCMContext->bSize + stepLocal], COMPLEX_R(src[step])));
        }

        if (savingFlag) {
            // Select accumulated value
            floatC *dataC = oddFrameFlag ? hCMContext->sCTmp1 : hCMContext->sCTmp2;

            // Drop first "half" frame
            if (hCMContext->frame > 0) {
                memcpy(dst + frames * 2 * hCMContext->harmonics, dataC, sizeof(floatC) * (size_t)hCMContext->harmonics);
                frames++;
            }

            // Set zeros
            for (int ih = 0; ih < hCMContext->harmonics; ih++) {
                dataC[ih] = COMPLEX_R(0);
            }

            // Increment global number of frames
            hCMContext->frame++;
        }
    }
    // Increment global number of steps
    hCMContext->stepOffset += (uint64_t)steps;

    avpkt->size = frames * hCMContext->harmonics * (int64_t)sizeof(floatC);
    *got_packet_ptr = 1;
    return 0;
}

static av_cold int hcm_encode_close(AVCodecContext *avctx)
{
    freeHCMContext((HCMContext *)avctx->priv_data);
    return 0;
}

static av_cold int hcm_decode_init(AVCodecContext *avctx)
{
    HCMContext *hCMContext = avctx->priv_data;
    initHCMContext(hCMContext, 0);

    hCMContext->step = 0;
    hCMContext->frameOffset = 0;

    avctx->sample_fmt = avctx->codec->sample_fmts[0];

    /*av_log(NULL, AV_LOG_INFO, "avctx->sample_rate: %d\n", avctx->sample_rate);
    av_log(NULL, AV_LOG_INFO, "avctx->channels: %d\n", avctx->channels);
    av_log(NULL, AV_LOG_INFO, "avctx->sample_fmt: %d\n", avctx->sample_fmt);
    av_log(NULL, AV_LOG_INFO, "avctx->frame_size: %d\n", avctx->frame_size);
    av_log(NULL, AV_LOG_INFO, "avctx->bits_per_coded_sample: %d\n", avctx->bits_per_coded_sample);
    av_log(NULL, AV_LOG_INFO, "avctx->block_align: %d\n", avctx->block_align);
    av_log(NULL, AV_LOG_INFO, "avctx->bit_rate: %d\n", avctx->bit_rate);
    av_log(NULL, AV_LOG_INFO, "avctx->channel_layout: %d\n", avctx->channel_layout);*/

    return 0;
}

static int hcm_decode_frame(AVCodecContext *avctx, void *data, int *got_frame_ptr, AVPacket *avpkt)
{
    HCMContext *hCMContext = avctx->priv_data;

    // avctx            codec context
    // data             output data (AVFrame)
    // got_frame_ptr    output data size
    // avpkt            input AVPacket
    AVFrame *frame = data;

    av_log(NULL, AV_LOG_INFO, "avpkt->size: %d\n", avpkt->size);

    if (avpkt->size == 0 && hCMContext->lastFrameDecoded)
        return 0;

    int frames = avpkt->size / (int)sizeof(floatC) / hCMContext->harmonics;

    // Because of last coefficients duplication
    if (avpkt->size == 0) {
        frames = 2;
        hCMContext->lastFrameDecoded = 1;
    }
    int maxOutputSteps = (frames) * hCMContext->oSize;

    int ret;
    frame->nb_samples = maxOutputSteps;
    if ((ret = ff_get_buffer(avctx, frame, 0)) < 0)
        return ret;

    float *dst = (float *)frame->data[0];
    const float *src = (const float *)avpkt->data;

    int steps = 0;
    // For every frame (input AVPacket)
    for (int frame = 0; frame < frames; frame++) {
        uint64_t framesOffsetGlobal = hCMContext->frameOffset + (uint64_t)frame;
        int frameOffset = frame * hCMContext->harmonics * 2;

        // Decode steps
        for (int stepLocal = 0; stepLocal < hCMContext->oSize; stepLocal++) {
            dst[steps] = 0;

            // For every hamornics
            for (int ih = 0; ih < hCMContext->harmonics; ih++) {

                if (stepLocal == 0) {
                    int iHC = 2 * ih;

                    // Save last coefficients
                    hCMContext->sCTmp1[ih] = hCMContext->sCTmp2[ih];

                    // Copy first coefficients
                    if (framesOffsetGlobal == 0)
                        hCMContext->sCTmp1[ih] = conjf(_FCOMPLEX_(src[iHC], src[iHC + 1]));

                    // Don't load last two coefficients (duplicate)
                    if (avpkt->size != 0) {
                        int fiHC = frameOffset + iHC;
                        // Read coefficient
                        hCMContext->sCTmp2[ih] = conjf(_FCOMPLEX_(src[fiHC], src[fiHC + 1]));
                    }
                }

                // Compute new point value
                int sH = ih * hCMContext->bSize + stepLocal;
                dst[steps] += crealf(cmulf(hCMContext->sCTmp2[ih], hCMContext->bE[sH])) + crealf(cmulf(hCMContext->sCTmp1[ih], hCMContext->bE_1[sH]));
            }
            steps++;
            hCMContext->step++;
        }
    }

    hCMContext->frameOffset += (uint64_t)frames;

    frame->nb_samples = steps;
    *got_frame_ptr = 1;
    return avpkt->size;
}

static av_cold int hcm_decode_close(AVCodecContext *avctx)
{
    freeHCMContext((HCMContext *)avctx->priv_data);
    return 0;
}

static const AVOption options[] = {
    {"period", "Period of input signal", offsetof(HCMContext, period), AV_OPT_TYPE_INT, {.i64 = 15}, 1, 10000, AV_OPT_FLAG_DECODING_PARAM | AV_OPT_FLAG_ENCODING_PARAM | AV_OPT_FLAG_AUDIO_PARAM },
    {"harmonics", "Multiple of harmonic frequency", offsetof(HCMContext, harmonics), AV_OPT_TYPE_INT, {.i64 = 2}, 1, 20, AV_OPT_FLAG_DECODING_PARAM | AV_OPT_FLAG_ENCODING_PARAM | AV_OPT_FLAG_AUDIO_PARAM },
    {"mos", "Multiple of overlap size", offsetof(HCMContext, mos), AV_OPT_TYPE_INT, {.i64 = 1}, 1, 100, AV_OPT_FLAG_DECODING_PARAM | AV_OPT_FLAG_ENCODING_PARAM | AV_OPT_FLAG_AUDIO_PARAM },
    {NULL}
};

static const AVClass hcm_encoder_class = {
    .class_name = "HCM encoder",
    .item_name  = av_default_item_name,
    .option     = options,
    .version    = LIBAVUTIL_VERSION_INT,
};

AVCodec ff_hcm_encoder = {
    .name           = "hcm",
    .long_name      = NULL_IF_CONFIG_SMALL("HCM (HIFU compression method)"),
    .type           = AVMEDIA_TYPE_AUDIO,
    .id             = AV_CODEC_ID_HCM,
    .priv_data_size = sizeof(HCMContext),
    .init           = hcm_encode_init,
    .close          = hcm_encode_close,
    .encode2        = hcm_encode_frame,
    .capabilities   = /*AV_CODEC_CAP_DR1 | */AV_CODEC_CAP_VARIABLE_FRAME_SIZE,
    .sample_fmts    = (const enum AVSampleFormat[]) { AV_SAMPLE_FMT_FLT,
                                                      AV_SAMPLE_FMT_NONE },
    .priv_class     = &hcm_encoder_class,
};

static const AVClass hcm_decoder_class = {
    .class_name = "HCM decoder",
    .item_name  = av_default_item_name,
    .option     = options,
    .version    = LIBAVUTIL_VERSION_INT,
};

AVCodec ff_hcm_decoder = {
    .name           = "hcm",
    .long_name      = NULL_IF_CONFIG_SMALL("HCM (HIFU compression method)"),
    .type           = AVMEDIA_TYPE_AUDIO,
    .id             = AV_CODEC_ID_HCM,
    .priv_data_size = sizeof(HCMContext),
    .init           = hcm_decode_init,
    .close          = hcm_decode_close,
    .decode         = hcm_decode_frame,
    .capabilities   = AV_CODEC_CAP_DR1 | AV_CODEC_CAP_DELAY,
    .sample_fmts    = (const enum AVSampleFormat[]) { AV_SAMPLE_FMT_FLT,
                                                      AV_SAMPLE_FMT_NONE },
    .priv_class     = &hcm_decoder_class,
};
