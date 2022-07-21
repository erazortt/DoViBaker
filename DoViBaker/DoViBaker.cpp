#include "DoViBaker.h"

//#include <stdlib.h>
//#include <math.h>
//#include <strstream>
#include <array>

//////////////////////////////
// Code
//////////////////////////////

template<bool chromaSubsampling>
DoViBaker<chromaSubsampling>::DoViBaker(PClip _blChild, PClip _elChild, const char* rpuPath, bool _quarterResolutionEl, IScriptEnvironment* env)
  : GenericVideoFilter(_blChild), elChild(_elChild), quarterResolutionEl(_quarterResolutionEl)
{
	int bits_per_pixel = vi.BitsPerComponent();
	if (bits_per_pixel != 16) {
		env->ThrowError("DoViBaker: Video must be 16bit!");
	}

	doviProc = new DoViProcessor(rpuPath, env);

  CPU_FLAG = env->GetCPUFlags();
}

template<bool chromaSubsampling>
DoViBaker<chromaSubsampling>::~DoViBaker()
{
}

template<bool chromaSubsampling>
template<int vertLen>
inline void DoViBaker<chromaSubsampling>::upsampleElVert(PVideoFrame& mez, const PVideoFrame& src, const std::array<int, vertLen>& pD, const std::array<int, vertLen>& nD, int p, IScriptEnvironment* env)
{
	static const int planes[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
	typedef uint16_t(*upscaler_t)(const uint16_t* srcSamples, int idx0);
	static const upscaler_t vertEvenUpscaler[3] = { &DoViProcessor::upsampleElYvertEven, &DoViProcessor::upsampleElUVvertEven, &DoViProcessor::upsampleElUVvertEven };
	static const upscaler_t vertOddUpscaler[3] = { &DoViProcessor::upsampleElYvertOdd, &DoViProcessor::upsampleElUVvertOdd, &DoViProcessor::upsampleElUVvertOdd };

	const int srcHeight = src->GetHeight(planes[p]);
	const int srcWidth = src->GetRowSize(planes[p]) / sizeof(uint16_t);
	const int srcPitch = src->GetPitch(planes[p]) / sizeof(uint16_t);
	const uint16_t* srcPb = (const uint16_t*)src->GetReadPtr(planes[p]);

	const int dstHeight = mez->GetHeight(planes[p]);
	const int dstWidth = mez->GetRowSize(planes[p]) / sizeof(uint16_t);
	const int dstPitch = mez->GetPitch(planes[p]) / sizeof(uint16_t);
	uint16_t* dstPeven = (uint16_t*)mez->GetWritePtr(planes[p]);
	uint16_t* dstPodd = dstPeven + dstPitch;

	std::array<const uint16_t*, vertLen + 1 + vertLen> srcP; // nD.size() + 1 + pD.size() = vertLen + 1 + vertLen
	std::array<uint16_t, srcP.size()> value;
	auto& srcP0 = srcP[nD.size()];

	for (int h0 = 0; h0 < srcHeight; h0++) {
		for (int i = 0; i < nD.size(); i++) {
			int factor = max(h0 + nD[i], 0);
			srcP[i] = srcPb + factor * srcPitch;
		}
		srcP0 = srcPb + h0 * srcPitch;
		for (int i = nD.size() + 1, j = 0; i < srcP.size(); i++, j++) {
			int factor = min(h0 + pD[j], srcHeight - 1);
			srcP[i] = srcPb + factor * srcPitch;
		}

		for (int w = 0; w < srcWidth; w++) {
			for (int i = 0; i < srcP.size(); i++) {
				value[i] = srcP[i][w];
			}
			dstPeven[w] = (*vertEvenUpscaler[p])(&value[0], nD.size());
			dstPodd[w] = (*vertOddUpscaler[p])(&value[0], nD.size());
		}

		dstPeven += 2 * dstPitch;
		dstPodd += 2 * dstPitch;
	}
}

template<bool chromaSubsampling>
void DoViBaker<chromaSubsampling>::upsampleElHorz(PVideoFrame& dst, const PVideoFrame& mez, int p, IScriptEnvironment* env)
{
	static const int planes[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
	const int srcHeight = mez->GetHeight(planes[p]);
	const int srcWidth = mez->GetRowSize(planes[p]) / sizeof(uint16_t);
	const int srcPitch = mez->GetPitch(planes[p]) / sizeof(uint16_t);
	const uint16_t* srcP = (const uint16_t*)mez->GetReadPtr(planes[p]);

	const int dstHeight = dst->GetHeight(planes[p]);
	const int dstWidth = dst->GetRowSize(planes[p]) / sizeof(uint16_t);
	const int dstPitch = dst->GetPitch(planes[p]) / sizeof(uint16_t);
	uint16_t* dstP = (uint16_t*)dst->GetWritePtr(planes[p]);

	static const std::array<int, 3> nD{ -3,-2,-1 };
	static const std::array<int, 4> pD{ 1, 2, 3, 4 };
	std::array<int, nD.size() + 1 + pD.size()> w;
	std::array<uint16_t, w.size()> value;
	int& w0 = w[nD.size()];

	for (int h = 0; h < srcHeight; h++) {
		for (w0 = nD.size(); w0 < srcWidth - pD.size(); w0++) {
			dstP[2 * w0] = doviProc->upsampleHorzEven(&srcP[w0 - nD.size()], nD.size());
			dstP[2 * w0 + 1] = doviProc->upsampleHorzOdd(&srcP[w0 - nD.size()], nD.size());
		}
		for (w0 = 0; w0 < nD.size(); w0++) {
			for (int i = 0; i < nD.size(); i++) {
				w[i] = max(w0 + nD[i], 0);
			}
			for (int i = 0; i < nD.size(); i++) {
				value[i] = srcP[w[i]];
			}
			std::copy_n(&srcP[w0], pD.size() + 1, &value[nD.size()]);
			dstP[2 * w0] = doviProc->upsampleHorzEven(&value[0], nD.size());
			dstP[2 * w0 + 1] = doviProc->upsampleHorzOdd(&value[0], nD.size());
		}
		for (w0 = srcWidth - pD.size(); w0 < srcWidth; w0++) {
			for (int i = nD.size() + 1, j = 0; i < w.size(); i++, j++) {
				w[i] = min(w0 + pD[j], srcWidth - 1);
			}
			for (int i = nD.size() + 1; i < w.size(); i++) {
				value[i] = srcP[w[i]];
			}
			std::copy_n(&srcP[w0 - nD.size()], nD.size() + 1, value.begin());
			dstP[2 * w0] = doviProc->upsampleHorzEven(&value[0], nD.size());
			dstP[2 * w0 + 1] = doviProc->upsampleHorzOdd(&value[0], nD.size());
		}
		srcP += srcPitch;
		dstP += dstPitch;
	}
}

template<bool chromaSubsampling>
void DoViBaker<chromaSubsampling>::upsampleEl(PVideoFrame& dst, const PVideoFrame& src, IScriptEnvironment* env)
{
	VideoInfo elvi = elChild->GetVideoInfo();
	elvi.height *= 2;
	PVideoFrame mez = env->NewVideoFrame(elvi);

	upsampleElVert<2>(mez, src, { -2,-1 }, { 2, 1 }, 0, env);
	upsampleElVert<1>(mez, src, { -1 }, { 1 }, 1, env);
	upsampleElVert<1>(mez, src, { -1 }, { 1 }, 2, env);
	
	upsampleElHorz(dst, mez, 0, env);
	upsampleElHorz(dst, mez, 1, env);
	upsampleElHorz(dst, mez, 2, env);
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

template<bool chromaSubsampling>
void DoViBaker<chromaSubsampling>::applyDovi(PVideoFrame& dst, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env) {

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

	// This section of code deals with the U and V planes of planar formats (e.g. YV12)
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

template<bool chromaSubsampling>
PVideoFrame DoViBaker<chromaSubsampling>::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame blSrc = child->GetFrame(n, env);
	PVideoFrame elSrc = elChild->GetFrame(n, env);
	PVideoFrame elUpSrc = env->NewVideoFrame(child->GetVideoInfo());
	PVideoFrame dst = env->NewVideoFrameP(vi, &blSrc);

	doviProc->intializeFrame(n, env);
	env->propSetInt(env->getFramePropsRW(dst), "_dovi_max_content_light_level", doviProc->getMaxContentLightLevel(), 0);
	std::string propBaseYcc2rgbCoef = "_dovi_ycc_to_rgb_coef";
	for (int i = 0; i < 9; i++) {
		env->propSetFloat(env->getFramePropsRW(dst), propBaseYcc2rgbCoef.append(std::to_string(i)).c_str() , doviProc->getYcc2RgbCoef(i), 0);
	}
	std::string propBaseYcc2rgbOff = "_dovi_ycc_to_rgb_offset";
	for (int i = 0; i < 3; i++) {
		env->propSetFloat(env->getFramePropsRW(dst), propBaseYcc2rgbOff.append(std::to_string(i)).c_str(), doviProc->getYcc2RgbOff(i), 0);
	}

	if (quarterResolutionEl) {
		upsampleEl(elUpSrc, elSrc, env);
		applyDovi(dst, blSrc, elUpSrc, env);
	}
	else {
		applyDovi(dst, blSrc, elSrc, env);
	}
	return dst;
}

// explicitly instantiate the template for the linker
template class DoViBaker<0>;
template class DoViBaker<1>;