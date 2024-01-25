#include "DoViTonemap.h"
#include "DoViProcessor.h"

// explicitly instantiate the template for the linker
template class DoViTonemap<8>;
template class DoViTonemap<10>;
template class DoViTonemap<12>;
template class DoViTonemap<14>;
template class DoViTonemap<16>;

template<int signalBitDepth>
DoViTonemap<signalBitDepth>::DoViTonemap(
	PClip child,
	float targetMaxNits,
	float targetMinNits,
	float masterMaxNits,
	float masterMinNits,
	float lumScale_,
	IScriptEnvironment* env)
	: GenericVideoFilter(child)
	, targetMaxPq(DoViProcessor::nits2pq(targetMaxNits))
	, targetMinPq(DoViProcessor::nits2pq(targetMinNits))
	, masterMaxPq(DoViProcessor::nits2pq(masterMaxNits < 0 ? 10000 : masterMaxNits))
	, masterMinPq(DoViProcessor::nits2pq(masterMinNits < 0 ? 0 : masterMinNits))
	, lumScale(lumScale_ < 0 ? 1 : lumScale_)
	, dynamicMasterMaxPq(masterMaxNits < 0)
	, dynamicMasterMinPq(masterMinNits < 0)
	, dynamicLumScale(lumScale < 0) 
{
	doviEetf = new DoViEetf<signalBitDepth>();
	doviEetf->generateEETF(
		targetMaxPq,
		targetMinPq,
		masterMaxPq,
		masterMinPq,
		lumScale);
}

template<int signalBitDepth>
DoViTonemap<signalBitDepth>::~DoViTonemap() 
{
	delete doviEetf;
	doviEetf = 0x0;
}

template<int signalBitDepth>
PVideoFrame DoViTonemap<signalBitDepth>::GetFrame(int n, IScriptEnvironment* env) {

	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrameP(vi, &src);

	uint16_t maxPq = masterMaxPq;
	uint16_t minPq = masterMinPq;
	float scale = lumScale;

	try {
		if (dynamicMasterMaxPq)
			maxPq = env->propGetInt(env->getFramePropsRO(src), "_dovi_dynamic_max_pq", 0, 0);
		if (dynamicMasterMinPq)
			minPq = env->propGetInt(env->getFramePropsRO(src), "_dovi_dynamic_min_pq", 0, 0);
		if (dynamicLumScale)
			scale = env->propGetFloat(env->getFramePropsRO(src), "_dovi_dynamic_luminosity_scale", 0, 0);
	}	catch(...) {
		env->ThrowError("DoViTonemap: Expected frame property not available");
		return dst;
	}

	if (maxPq != masterMaxPq || minPq != masterMinPq || scale != lumScale) {
		masterMinPq = minPq;
		masterMaxPq = maxPq;
		lumScale = scale;
		doviEetf->generateEETF(
			targetMaxPq,
			targetMinPq,
			masterMaxPq,
			masterMinPq,
			lumScale);
	}

	applyTonemapRGB(dst, src);

	return dst;
}

template<int signalBitDepth>
void DoViTonemap<signalBitDepth>::applyTonemapRGB(PVideoFrame& dst, const PVideoFrame& src) const
{
	//apply tonemap using R'G'B' scaling
	const unsigned int height = src->GetHeight(PLANAR_R);
	const unsigned int width = src->GetRowSize(PLANAR_R) / sizeof(uint16_t);
	const uint16_t* srcP[3];
	int srcPitch[3];
	uint16_t* dstP[3];
	int dstPitch[3];

	srcP[0] = (const uint16_t*)src->GetReadPtr(PLANAR_R);
	srcP[1] = (const uint16_t*)src->GetReadPtr(PLANAR_G);
	srcP[2] = (const uint16_t*)src->GetReadPtr(PLANAR_B);
	srcPitch[0] = src->GetPitch(PLANAR_R) / sizeof(uint16_t);
	srcPitch[1] = src->GetPitch(PLANAR_G) / sizeof(uint16_t);
	srcPitch[2] = src->GetPitch(PLANAR_B) / sizeof(uint16_t);
	dstP[0] = (uint16_t*)dst->GetWritePtr(PLANAR_R);
	dstP[1] = (uint16_t*)dst->GetWritePtr(PLANAR_G);
	dstP[2] = (uint16_t*)dst->GetWritePtr(PLANAR_B);
	dstPitch[0] = dst->GetPitch(PLANAR_R) / sizeof(uint16_t);
	dstPitch[1] = dst->GetPitch(PLANAR_G) / sizeof(uint16_t);
	dstPitch[2] = dst->GetPitch(PLANAR_B) / sizeof(uint16_t);

	for (unsigned p = 0; p < 3; ++p)
	{
		for (unsigned h = 0; h < height; ++h)
		{
			for (unsigned w = 0; w < width; ++w) {
				dstP[p][w] = doviEetf->applyEETF(srcP[p][w]);
			}
			srcP[p] += srcPitch[p];
			dstP[p] += dstPitch[p];
		}
	}
}

/*template<int signalBitDepth>
void DoViTonemap<signalBitDepth>::applyTonemapMaxRGB(PVideoFrame& dst, const PVideoFrame& src) const
{
	//apply tonemap using maxRGB scaling
	const unsigned int height = src->GetHeight(PLANAR_R);
	const unsigned int width = src->GetRowSize(PLANAR_R) / sizeof(uint16_t);
	const uint16_t* srcP[3];
	int srcPitch[3];
	uint16_t* dstP[3];
	int dstPitch[3];

	srcP[0] = (const uint16_t*)src->GetReadPtr(PLANAR_R);
	srcP[1] = (const uint16_t*)src->GetReadPtr(PLANAR_G);
	srcP[2] = (const uint16_t*)src->GetReadPtr(PLANAR_B);
	srcPitch[0] = src->GetPitch(PLANAR_R) / sizeof(uint16_t);
	srcPitch[1] = src->GetPitch(PLANAR_G) / sizeof(uint16_t);
	srcPitch[2] = src->GetPitch(PLANAR_B) / sizeof(uint16_t);
	dstP[0] = (uint16_t*)dst->GetWritePtr(PLANAR_R);
	dstP[1] = (uint16_t*)dst->GetWritePtr(PLANAR_G);
	dstP[2] = (uint16_t*)dst->GetWritePtr(PLANAR_B);
	dstPitch[0] = dst->GetPitch(PLANAR_R) / sizeof(uint16_t);
	dstPitch[1] = dst->GetPitch(PLANAR_G) / sizeof(uint16_t);
	dstPitch[2] = dst->GetPitch(PLANAR_B) / sizeof(uint16_t);

	for (unsigned h = 0; h < height; ++h)
	{
		for (unsigned w = 0; w < width; ++w) {
			uint16_t m1 = std::max(std::max(srcP[0][w], srcP[1][w]), srcP[2][w]);
			float m2 = pq2signal(applyLut(signal2pq(m1)));
			float factor = m2 / m1;
			for (int i = 0; i < 3; i++)
				//dstP[i][w] = (nits2pq(factor * pq2nits(srcP[i][w] >> (16 - 12))) << 16 - 12);
				dstP[i][w] = pq2signal(nits2pq(factor * pq2nits(signal2pq(srcP[i][w]))));
		}

		for (unsigned p = 0; p < 3; ++p)
		{
			srcP[p] += srcPitch[p];
			dstP[p] += dstPitch[p];
		}
	}
}*/
