#include "DoViBaker.h"

//#include <stdlib.h>
//#include <math.h>
//#include <strstream>
#include <array>

//////////////////////////////
// Code
//////////////////////////////

template<bool chromaSubsampling, bool quarterResolutionEl>
DoViBaker<chromaSubsampling, quarterResolutionEl>::DoViBaker(PClip _blChild, PClip _elChild, const char* rpuPath, bool _qnd, IScriptEnvironment* env)
  : GenericVideoFilter(_blChild), elChild(_elChild), qnd(_qnd)
{
	int bits_per_pixel = vi.BitsPerComponent();
	if (bits_per_pixel != 16) {
		env->ThrowError("DoViBaker: Video must be 16bit!");
	}
	vi.pixel_type = VideoInfo::CS_RGBP16;

	doviProc = new DoViProcessor(rpuPath, env);

  CPU_FLAG = env->GetCPUFlags();
}

template<bool chromaSubsampling, bool quarterResolutionEl>
DoViBaker<chromaSubsampling, quarterResolutionEl>::~DoViBaker()
{
}

template<bool chromaSubsampling, bool quarterResolutionEl>
template<int vertLen, int nD>
inline void DoViBaker<chromaSubsampling, quarterResolutionEl>::upsampleVert(PVideoFrame& mez, const PVideoFrame& src, const int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env)
{
	const int srcHeight = src->GetHeight(plane);
	const int srcWidth = src->GetRowSize(plane) / sizeof(uint16_t);
	const int srcPitch = src->GetPitch(plane) / sizeof(uint16_t);
	const uint16_t* srcPb = (const uint16_t*)src->GetReadPtr(plane);

	const int dstHeight = mez->GetHeight(plane);
	const int dstWidth = mez->GetRowSize(plane) / sizeof(uint16_t);
	const int dstPitch = mez->GetPitch(plane) / sizeof(uint16_t);
	uint16_t* dstPeven = (uint16_t*)mez->GetWritePtr(plane);
	uint16_t* dstPodd = dstPeven + dstPitch;

	std::array<const uint16_t*, vertLen> srcP;
	std::array<uint16_t, vertLen> value;
	auto& srcP0 = srcP[nD];

	for (int h0 = 0; h0 < srcHeight; h0++) {
		for (int i = 0; i < nD; i++) {
			int factor = max(h0 + Dn0p[i], 0);
			srcP[i] = srcPb + factor * srcPitch;
		}
		srcP0 = srcPb + h0 * srcPitch;
		for (int i = nD + 1; i < vertLen; i++) {
			int factor = min(h0 + Dn0p[i], srcHeight - 1);
			srcP[i] = srcPb + factor * srcPitch;
		}

		for (int w = 0; w < srcWidth; w++) {
			for (int i = 0; i < vertLen; i++) {
				value[i] = srcP[i][w];
			}
			dstPeven[w] = evenUpscaler(&value[0], nD);
			dstPodd[w] = oddUpscaler(&value[0], nD);
		}

		dstPeven += 2 * dstPitch;
		dstPodd += 2 * dstPitch;
	}
}

template<bool chromaSubsampling, bool quarterResolutionEl>
void DoViBaker<chromaSubsampling, quarterResolutionEl>::upsampleHorz(PVideoFrame& dst, const PVideoFrame& mez, int plane, IScriptEnvironment* env)
{
	const int srcHeight = mez->GetHeight(plane);
	const int srcWidth = mez->GetRowSize(plane) / sizeof(uint16_t);
	const int srcPitch = mez->GetPitch(plane) / sizeof(uint16_t);
	const uint16_t* srcP = (const uint16_t*)mez->GetReadPtr(plane);

	const int dstHeight = dst->GetHeight(plane);
	const int dstWidth = dst->GetRowSize(plane) / sizeof(uint16_t);
	const int dstPitch = dst->GetPitch(plane) / sizeof(uint16_t);
	uint16_t* dstP = (uint16_t*)dst->GetWritePtr(plane);

	static const std::array<int, 8> Dn0p{ -3,-2,-1, 0, 1, 2, 3, 4 };
	static const int nD = 3;
	static const int pD = Dn0p.size() - nD - 1;
	std::array<uint16_t, Dn0p.size()> value;
	
	for (int h = 0; h < srcHeight; h++) {
		for (int w = nD; w < srcWidth - pD; w++) {
			dstP[2 * w] = doviProc->upsampleHorzEven(&srcP[w - nD], nD);
			dstP[2 * w + 1] = doviProc->upsampleHorzOdd(&srcP[w - nD], nD);
		}
		for (int w = 0; w < nD; w++) {
			for (int i = 0; i < nD; i++) {
				int wd = max(w + Dn0p[i], 0);
				value[i] = srcP[wd];
			}
			std::copy_n(&srcP[w], pD + 1, &value[nD]);
			dstP[2 * w] = doviProc->upsampleHorzEven(&value[0], nD);
			dstP[2 * w + 1] = doviProc->upsampleHorzOdd(&value[0], nD);
		}
		for (int w = srcWidth - pD; w < srcWidth; w++) {
			for (int i = nD + 1; i < Dn0p.size(); i++) {
				int wd = min(w + Dn0p[i], srcWidth - 1);
				value[i] = srcP[wd];
			}
			std::copy_n(&srcP[w - nD], nD + 1, &value[0]);
			dstP[2 * w] = doviProc->upsampleHorzEven(&value[0], nD);
			dstP[2 * w + 1] = doviProc->upsampleHorzOdd(&value[0], nD);
		}
		srcP += srcPitch;
		dstP += dstPitch;
	}
}

template<bool chromaSubsampling, bool quarterResolutionEl>
void DoViBaker<chromaSubsampling, quarterResolutionEl>::upsampleEl(PVideoFrame& dst, const PVideoFrame& src, VideoInfo dstVi, IScriptEnvironment* env)
{
	dstVi.width /= 2;
	PVideoFrame mez = env->NewVideoFrame(dstVi);

	upsampleVert<5,2>(mez, src, PLANAR_Y, { -2,-1, 0, 1, 2 }, &DoViProcessor::upsampleElYvertEven, &DoViProcessor::upsampleElYvertOdd, env);
	upsampleVert<3,1>(mez, src, PLANAR_U, { -1, 0, 1 }, &DoViProcessor::upsampleElUVvertEven, &DoViProcessor::upsampleElUVvertOdd, env);
	upsampleVert<3,1>(mez, src, PLANAR_V, { -1, 0, 1 }, &DoViProcessor::upsampleElUVvertEven, &DoViProcessor::upsampleElUVvertOdd, env);
	
	upsampleHorz(dst, mez, PLANAR_Y, env);
	upsampleHorz(dst, mez, PLANAR_U, env);
	upsampleHorz(dst, mez, PLANAR_V, env);
}

template<bool chromaSubsampling, bool quarterResolutionEl>
void DoViBaker<chromaSubsampling, quarterResolutionEl>::to444(PVideoFrame& dst, const PVideoFrame& src, VideoInfo dstVi, IScriptEnvironment* env)
{
	dstVi.width /= 2;
	PVideoFrame mez = env->NewVideoFrame(dstVi);
	
	upsampleVert<7,3>(mez, src, PLANAR_U, { -3,-2,-1, 0, 1, 2, 3 }, &DoViProcessor::upsampleElUVvertEven, &DoViProcessor::upsampleElUVvertOdd, env);
	upsampleVert<7,3>(mez, src, PLANAR_V, { -3,-2,-1, 0, 1, 2, 3 }, &DoViProcessor::upsampleElUVvertEven, &DoViProcessor::upsampleElUVvertOdd, env);

	upsampleHorz(dst, mez, PLANAR_U, env);
	upsampleHorz(dst, mez, PLANAR_V, env);
}

/*PVideoFrame DoViBaker::GetFrame(int n, IScriptEnvironment* env) {
	PVideoFrame src = elChild->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrame(vi);
	
	static const int srcChannel[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };

	if (vi.IsPlanar()) {

		for (int i = 0; i < 3; i++) {
			
			const int src_pitch = src->GetPitch(srcChannel[i]);
			const int src_width = src->GetRowSize(srcChannel[i]);
			const int src_height = src->GetHeight(srcChannel[i]);
			const unsigned char* srcp = src->GetReadPtr(srcChannel[i]);

			const int dst_pitch = dst->GetPitch(srcChannel[i]);
			const int dst_width = dst->GetRowSize(srcChannel[i]);
			const int dst_height = dst->GetHeight(srcChannel[i]);
			unsigned char* dstp = dst->GetWritePtr(srcChannel[i]);

			for (int h = 0; h < src_height; h++) {
				for (int w = 0; w < src_width; w++) {
					*(dstp + w) = *(srcp + w);
				}
				srcp += src_pitch;
				dstp += dst_pitch;
			}
		}
		return dst;
	}
}*/

template<bool chromaSubsampling, bool quarterResolutionEl>
void DoViBaker<chromaSubsampling, quarterResolutionEl>::applyDovi(PVideoFrame& dst, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env) {

	const int blSrcHeightY = blSrc->GetHeight();
	const int blSrcWidthY = blSrc->GetRowSize() / sizeof(uint16_t);
	const int blSrcPitchY = blSrc->GetPitch() / sizeof(uint16_t);

	const int elSrcHeightY = elSrc->GetHeight();
	const int elSrcWidthY = elSrc->GetRowSize() / sizeof(uint16_t);
	const int elSrcPitchY = elSrc->GetPitch() / sizeof(uint16_t);

	const int dstHeightY = dst->GetHeight();
	const int dstWidthY = dst->GetRowSize() / sizeof(uint16_t);
	const int dstPitchY = dst->GetPitch() / sizeof(uint16_t);

	std::array<const uint16_t*, chromaSubsampling + 1> blSrcYp;
	std::array<const uint16_t*, chromaSubsampling + 1> elSrcYp;
	std::array<uint16_t*, chromaSubsampling + 1> dstYp;

	blSrcYp[0] = (const uint16_t*)blSrc->GetReadPtr();
	elSrcYp[0] = (const uint16_t*)elSrc->GetReadPtr();
	dstYp[0] = (uint16_t*)dst->GetWritePtr();
	if (chromaSubsampling) {
		blSrcYp[1] = blSrcYp[0] + blSrcPitchY;
		elSrcYp[1] = elSrcYp[0] + elSrcPitchY;
		dstYp[1] = dstYp[0] + dstPitchY;
	}

	const int blSrcHeightUV = blSrc->GetHeight(PLANAR_U);
	const int blSrcWidthUV = blSrc->GetRowSize(PLANAR_U) / sizeof(uint16_t);
	const int blSrcPitchUV = blSrc->GetPitch(PLANAR_U) / sizeof(uint16_t);

	const int elSrcHeightUV = elSrc->GetHeight(PLANAR_U);
	const int elSrcWidthUV = elSrc->GetRowSize(PLANAR_U) / sizeof(uint16_t);
	const int elSrcPitchUV = elSrc->GetPitch(PLANAR_U) / sizeof(uint16_t);

	const int dstHeightUV = dst->GetHeight(PLANAR_U);
	const int dstWidthUV = dst->GetRowSize(PLANAR_U) / sizeof(uint16_t);
	const int dstPitchUV = dst->GetPitch(PLANAR_U) / sizeof(uint16_t);

	const uint16_t* blSrcUp = (const uint16_t*)blSrc->GetReadPtr(PLANAR_U);
	const uint16_t* elSrcUp = (const uint16_t*)elSrc->GetReadPtr(PLANAR_U);
	uint16_t* dstUp = (uint16_t*)dst->GetWritePtr(PLANAR_U);

	const uint16_t* blSrcVp = (const uint16_t*)blSrc->GetReadPtr(PLANAR_V);
	const uint16_t* elSrcVp = (const uint16_t*)elSrc->GetReadPtr(PLANAR_V);
	uint16_t* dstVp = (uint16_t*)dst->GetWritePtr(PLANAR_V);

	for (int huv = 0; huv < blSrcHeightUV; huv++) {
		if (chromaSubsampling) {
			int wuv = 0;
			for (int j = 0; j < chromaSubsampling + 1; j++) {
				for (int i = 0; i < chromaSubsampling + 1; i++) {
					const int w = (chromaSubsampling + 1) * wuv + i;
					dstYp[j][w] = doviProc->processSampleY(blSrcYp[j][w], elSrcYp[j][w]);
				}
			}
			int mmrBlY1 = 3 * blSrcYp[0][2 * wuv] + blSrcYp[0][2 * wuv + 1] + 2;
			int mmrBlY2 = 3 * blSrcYp[1][2 * wuv] + blSrcYp[1][2 * wuv + 1] + 2;
			uint16_t mmrBlY = ((mmrBlY1 >> 2) + (mmrBlY2 >> 2) + 1) >> 1;

			dstUp[wuv] = doviProc->processSampleU(blSrcUp[wuv], elSrcUp[wuv], mmrBlY, blSrcUp[wuv], blSrcVp[wuv]);
			dstVp[wuv] = doviProc->processSampleV(blSrcVp[wuv], elSrcVp[wuv], mmrBlY, blSrcUp[wuv], blSrcVp[wuv]);
		}

		for (int wuv = chromaSubsampling; wuv < blSrcWidthUV - chromaSubsampling; wuv++) {
			for (int j = 0; j < chromaSubsampling + 1; j++) {
				for (int i = 0; i < chromaSubsampling + 1; i++) {
					const int w = (chromaSubsampling + 1) * wuv + i;
					dstYp[j][w] = doviProc->processSampleY(blSrcYp[j][w], elSrcYp[j][w]);
				}
			}
			uint16_t mmrBlY;
			if (chromaSubsampling) {
				int mmrBlY1 = blSrcYp[0][2 * wuv - 1] + 2 * blSrcYp[0][2 * wuv] + blSrcYp[0][2 * wuv + 1] + 2;
				int mmrBlY2 = blSrcYp[1][2 * wuv - 1] + 2 * blSrcYp[1][2 * wuv] + blSrcYp[1][2 * wuv + 1] + 2;
				mmrBlY = ((mmrBlY1 >> 2) + (mmrBlY2 >> 2) + 1) >> 1;
			}
			else {
				mmrBlY = blSrcYp[0][wuv];
			}
			dstUp[wuv] = doviProc->processSampleU(blSrcUp[wuv], elSrcUp[wuv], mmrBlY, blSrcUp[wuv], blSrcVp[wuv]);
			dstVp[wuv] = doviProc->processSampleV(blSrcVp[wuv], elSrcVp[wuv], mmrBlY, blSrcUp[wuv], blSrcVp[wuv]);
		}

		if (chromaSubsampling) {
			int wuv = blSrcWidthUV - chromaSubsampling;
			for (int j = 0; j < chromaSubsampling + 1; j++) {
				for (int i = 0; i < chromaSubsampling + 1; i++) {
					const int w = (chromaSubsampling + 1) * wuv + i;
					dstYp[j][w] = doviProc->processSampleY(blSrcYp[j][w], elSrcYp[j][w]);
				}
			}
			int mmrBlY1 = blSrcYp[0][2 * wuv - 1] + 3 * blSrcYp[0][2 * wuv] + 2;
			int mmrBlY2 = blSrcYp[1][2 * wuv - 1] + 3 * blSrcYp[1][2 * wuv] + 2;
			uint16_t mmrBlY = ((mmrBlY1 >> 2) + (mmrBlY2 >> 2) + 1) >> 1;

			dstUp[wuv] = doviProc->processSampleU(blSrcUp[wuv], elSrcUp[wuv], mmrBlY, blSrcUp[wuv], blSrcVp[wuv]);
			dstVp[wuv] = doviProc->processSampleV(blSrcVp[wuv], elSrcVp[wuv], mmrBlY, blSrcUp[wuv], blSrcVp[wuv]);
		}

		for (int i = 0; i < chromaSubsampling + 1; i++) {
			blSrcYp[i] += blSrcPitchY * (chromaSubsampling + 1);
			elSrcYp[i] += elSrcPitchY * (chromaSubsampling + 1);
			dstYp[i] += dstPitchY * (chromaSubsampling + 1);
		}
		blSrcUp += blSrcPitchUV;
		blSrcVp += blSrcPitchUV;
		elSrcUp += elSrcPitchUV;
		elSrcVp += elSrcPitchUV;
		dstUp += dstPitchUV;
		dstVp += dstPitchUV;
	}
}

template<bool chromaSubsampling, bool quarterResolutionEl>
void DoViBaker<chromaSubsampling, quarterResolutionEl>::doAllQuickAndDirty(PVideoFrame& dst, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env) {
	const int blSrcHeightY = blSrc->GetHeight(PLANAR_Y);
	const int blSrcWidthY = blSrc->GetRowSize(PLANAR_Y) / sizeof(uint16_t);
	const int blSrcPitchY = blSrc->GetPitch(PLANAR_Y) / sizeof(uint16_t);

	const int elSrcHeightY = elSrc->GetHeight(PLANAR_Y);
	const int elSrcWidthY = elSrc->GetRowSize(PLANAR_Y) / sizeof(uint16_t);
	const int elSrcPitchY = elSrc->GetPitch(PLANAR_Y) / sizeof(uint16_t);

	const int dstHeight = dst->GetHeight(PLANAR_R);
	const int dstWidth = dst->GetRowSize(PLANAR_R) / sizeof(uint16_t);
	const int dstPitch = dst->GetPitch(PLANAR_R) / sizeof(uint16_t);

	const int blSrcHeightUV = blSrc->GetHeight(PLANAR_U);
	const int blSrcWidthUV = blSrc->GetRowSize(PLANAR_U) / sizeof(uint16_t);
	const int blSrcPitchUV = blSrc->GetPitch(PLANAR_U) / sizeof(uint16_t);

	const int elSrcHeightUV = elSrc->GetHeight(PLANAR_U);
	const int elSrcWidthUV = elSrc->GetRowSize(PLANAR_U) / sizeof(uint16_t);
	const int elSrcPitchUV = elSrc->GetPitch(PLANAR_U) / sizeof(uint16_t);


	std::array<const uint16_t*, (chromaSubsampling + 1) * (quarterResolutionEl + 1)> blSrcYp;
	std::array<const uint16_t*, chromaSubsampling + 1> elSrcYp;
	std::array<uint16_t*, (chromaSubsampling + 1)* (quarterResolutionEl + 1)> dstRp;
	blSrcYp[0] = (const uint16_t*)blSrc->GetReadPtr(PLANAR_Y);
	elSrcYp[0] = (const uint16_t*)elSrc->GetReadPtr(PLANAR_Y);
	dstRp[0] = (uint16_t*)dst->GetWritePtr(PLANAR_R);

	std::array<const uint16_t*, quarterResolutionEl + 1> blSrcUp;
	std::array<const uint16_t*, 1> elSrcUp;
	std::array<uint16_t*, (chromaSubsampling + 1)* (quarterResolutionEl + 1)> dstGp;
	blSrcUp[0] = (const uint16_t*)blSrc->GetReadPtr(PLANAR_U);
	elSrcUp[0] = (const uint16_t*)elSrc->GetReadPtr(PLANAR_U);
	dstGp[0] = (uint16_t*)dst->GetWritePtr(PLANAR_G);

	std::array<const uint16_t*, quarterResolutionEl + 1> blSrcVp;
	std::array<const uint16_t*, 1> elSrcVp;
	std::array<uint16_t*, (chromaSubsampling + 1)* (quarterResolutionEl + 1)> dstBp;
	blSrcVp[0] = (const uint16_t*)blSrc->GetReadPtr(PLANAR_V);
	elSrcVp[0] = (const uint16_t*)elSrc->GetReadPtr(PLANAR_V);
	dstBp[0] = (uint16_t*)dst->GetWritePtr(PLANAR_B);


	for (int i = 1; i < elSrcVp.size(); i++) {
		elSrcUp[i] = elSrcUp[i - 1] + elSrcPitchUV;
		elSrcVp[i] = elSrcVp[i - 1] + elSrcPitchUV;
	}
	for (int i = 1; i < elSrcYp.size(); i++) {
		elSrcYp[i] = elSrcYp[i - 1] + elSrcPitchY;
	}
	for (int i = 1; i < blSrcVp.size(); i++) {
		blSrcUp[i] = blSrcUp[i - 1] + blSrcPitchUV;
		blSrcVp[i] = blSrcVp[i - 1] + blSrcPitchUV;
	}
	for (int i = 1; i < dstBp.size(); i++) {
		blSrcYp[i] = blSrcYp[i - 1] + blSrcPitchY;
		dstRp[i] = dstRp[i - 1] + dstPitch;
		dstGp[i] = dstGp[i - 1] + dstPitch;
		dstBp[i] = dstBp[i - 1] + dstPitch;
	}

	for (int heluv = 0; heluv < elSrcHeightUV; heluv++) {
		for (int weluv = 0; weluv < elSrcWidthUV; weluv++) {

			const uint16_t& elu = elSrcUp[0][weluv];
			const uint16_t& elv = elSrcVp[0][weluv];

			for (int hDbluv = 0; hDbluv < quarterResolutionEl + 1; hDbluv++) {
				for (int wDbluv = 0; wDbluv < quarterResolutionEl + 1; wDbluv++) {

					int wbluv = (weluv << quarterResolutionEl) + wDbluv;
					const uint16_t& blu = blSrcUp[hDbluv][wbluv];
					const uint16_t& blv = blSrcVp[hDbluv][wbluv];

					int hDbluvy = hDbluv << chromaSubsampling;
					int wbluvy = wbluv << chromaSubsampling;
					const uint16_t& mmrbly = blSrcYp[hDbluvy][wbluvy];

					const uint16_t& u = doviProc->processSampleU(blu, elu, mmrbly, blu, blv);
					const uint16_t& v = doviProc->processSampleV(blv, elv, mmrbly, blu, blv);

					for (int hDbly = 0; hDbly < chromaSubsampling + 1; hDbly++) {
						for (int wDbly = 0; wDbly < chromaSubsampling + 1; wDbly++) {

							int hDDbly = hDbluvy + hDbly;
							int wbly = wbluvy + wDbly;
							const uint16_t& bly = blSrcYp[hDDbly][wbly];

							int hDely = hDDbly >> quarterResolutionEl;
							int wely = wbly >> quarterResolutionEl;
							const uint16_t& ely = elSrcYp[hDely][wely];

							const uint16_t& y = doviProc->processSampleY(bly, ely);
							doviProc->sample2rgb(dstRp[hDDbly][wbly], dstGp[hDDbly][wbly], dstBp[hDDbly][wbly], y, u, v);
						}
					}
				}
			}
		}

		for (int i = 0; i < elSrcVp.size(); i++) {
			elSrcVp[i] += elSrcPitchUV;
			elSrcUp[i] += elSrcPitchUV;
		}
		for (int i = 0; i < chromaSubsampling + 1; i++) {
			elSrcYp[i] += elSrcPitchY * (chromaSubsampling + 1);
		}
		for (int i = 0; i < quarterResolutionEl + 1; i++) {
			blSrcVp[i] += blSrcPitchUV * (quarterResolutionEl + 1);
			blSrcUp[i] += blSrcPitchUV * (quarterResolutionEl + 1);
		}
		for (int i = 0; i < (quarterResolutionEl + 1) * (chromaSubsampling + 1); i++) {
			blSrcYp[i] += blSrcPitchY * (quarterResolutionEl + 1) * (chromaSubsampling + 1);
			dstRp[i] += dstPitch * (quarterResolutionEl + 1) * (chromaSubsampling + 1);
			dstGp[i] += dstPitch * (quarterResolutionEl + 1) * (chromaSubsampling + 1);
			dstBp[i] += dstPitch * (quarterResolutionEl + 1) * (chromaSubsampling + 1);
		}
	}
}

template<bool chromaSubsampling, bool quarterResolutionEl>
void DoViBaker<chromaSubsampling, quarterResolutionEl>::convert2rgb(PVideoFrame& dst, const PVideoFrame& srcY, const PVideoFrame& srcUV)
{
	const int srcHeightY = srcY->GetHeight(PLANAR_Y);
	const int srcWidthY = srcY->GetRowSize(PLANAR_Y) / sizeof(uint16_t);
	const int srcPitchY = srcY->GetPitch(PLANAR_Y) / sizeof(uint16_t);

	const int dstHeight = dst->GetHeight(PLANAR_R);
	const int dstWidth = dst->GetRowSize(PLANAR_R) / sizeof(uint16_t);
	const int dstPitch = dst->GetPitch(PLANAR_R) / sizeof(uint16_t);

	const uint16_t* srcYp = (const uint16_t*)srcY->GetReadPtr(PLANAR_Y);
	uint16_t* dstRp = (uint16_t*)dst->GetWritePtr(PLANAR_R);

	const int srcHeightUV = srcUV->GetHeight(PLANAR_U);
	const int srcWidthUV = srcUV->GetRowSize(PLANAR_U) / sizeof(uint16_t);
	const int srcPitchUV = srcUV->GetPitch(PLANAR_U) / sizeof(uint16_t);

	const uint16_t* srcUp = (const uint16_t*)srcUV->GetReadPtr(PLANAR_U);
	uint16_t* dstGp = (uint16_t*)dst->GetWritePtr(PLANAR_G);

	const uint16_t* srcVp = (const uint16_t*)srcUV->GetReadPtr(PLANAR_V);
	uint16_t* dstBp = (uint16_t*)dst->GetWritePtr(PLANAR_B);

	for (int huv = 0; huv < srcHeightUV; huv++) {
		for (int wuv = 0; wuv < srcWidthUV; wuv++) {
			doviProc->sample2rgb(dstRp[wuv], dstGp[wuv], dstBp[wuv], srcYp[wuv], srcUp[wuv], srcVp[wuv]);
		}

		srcYp += srcPitchY;
		srcUp += srcPitchUV;
		srcVp += srcPitchUV;

		dstRp += dstPitch;
		dstGp += dstPitch;
		dstBp += dstPitch;
	}
}

template<bool chromaSubsampling, bool quarterResolutionEl>
PVideoFrame DoViBaker<chromaSubsampling, quarterResolutionEl>::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame blSrc = child->GetFrame(n, env);
	PVideoFrame elSrc = elChild ? elChild->GetFrame(n, env) : blSrc;
	PVideoFrame dst = env->NewVideoFrameP(vi, &blSrc);

	doviProc->intializeFrame(n, env);
	env->propSetInt(env->getFramePropsRW(dst), "_dovi_max_content_light_level", doviProc->getMaxContentLightLevel(), 0);

	bool skipElProcessing = false;
	if (!elChild || !doviProc->isFEL()) {
		skipElProcessing = true;
		doviProc->forceDisableElProcessing();
	}

	if (qnd) {
		doAllQuickAndDirty(dst, blSrc, elSrc, env);
	}
	else {
		PVideoFrame mez = env->NewVideoFrame(child->GetVideoInfo());
		if (quarterResolutionEl && !skipElProcessing) {
			PVideoFrame elUpSrc = env->NewVideoFrame(child->GetVideoInfo());
			upsampleEl(elUpSrc, elSrc, child->GetVideoInfo(), env);
			applyDovi(mez, blSrc, elUpSrc, env);
		}
		else {
			applyDovi(mez, blSrc, elSrc, env);
		}
		if (chromaSubsampling) {
			VideoInfo vi444 = child->GetVideoInfo();
			vi444.pixel_type = VideoInfo::CS_YUV444P16;
			PVideoFrame mez444 = env->NewVideoFrame(vi444);
			to444(mez444, mez, vi444, env);
			convert2rgb(dst, mez, mez444);
		}
		else {
			convert2rgb(dst, mez, mez);
		}
	}
	return dst;
}

// explicitly instantiate the template for the linker
template class DoViBaker<0, 0>;
template class DoViBaker<0, 1>;
template class DoViBaker<1, 0>;
template class DoViBaker<1, 1>;