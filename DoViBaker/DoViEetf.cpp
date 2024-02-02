#include "DoViEetf.h"
#include "DoViProcessor.h"

// explicitly instantiate the template for the linker
template class DoViEetf<10>;
template class DoViEetf<12>;
template class DoViEetf<14>;
template class DoViEetf<16>;

template<int signalBitDepth>
DoViEetf<signalBitDepth>::DoViEetf(float kneeOffset_, bool normalizeOutput_)
	: kneeOffset(kneeOffset_),
		normalizeOutput(normalizeOutput_) {}

template<int signalBitDepth>
void DoViEetf<signalBitDepth>::generateEETF(
	uint16_t targetMaxPq,
	uint16_t targetMinPq,
	uint16_t masterMaxPq,
	uint16_t masterMinPq,
	float lumScale)
{
	// based on report ITU-R BT.2408-7 Annex 5 (was in ITU-R BT.2390 until revision 7)
	float masterMaxEp = DoViProcessor::EOTFinv(DoViProcessor::EOTF(masterMaxPq / 4095.0f) * lumScale);
	float masterMinEp = DoViProcessor::EOTFinv(DoViProcessor::EOTF(masterMinPq / 4095.0f) * lumScale);
	float targetMaxEp = targetMaxPq / 4095.0f;
	float targetMinEp = targetMinPq / 4095.0f;

	// this is not from the report
	// any mapping is unnecessary while inside the range of tragets capabilities
	masterMaxEp = fmaxf(masterMaxEp, targetMaxEp);
	masterMinEp = fminf(masterMinEp, targetMinEp);

	float maxLum = (targetMaxEp - masterMinEp) / (masterMaxEp - masterMinEp);
	float minLum = (targetMinEp - masterMinEp) / (masterMaxEp - masterMinEp);
	float kneeStart = (1 + kneeOffset) * maxLum - kneeOffset;
	float b = minLum;

	for (int inSignal = 0; inSignal < LUT_SIZE; inSignal++) {
		float ep = inSignal / float(LUT_SIZE - 1);
		float e = DoViProcessor::EOTF(ep) * lumScale;
		ep = DoViProcessor::EOTFinv(e);
		float e1 = (ep - masterMinEp) / (masterMaxEp - masterMinEp);

		// this following claming is not from the report
		// it serves to stop using the spline above where it is supposed to be used 
		// this works only in conjunction with the change above
		e1 = std::clamp(e1, 0.0f, 1.0f); 

		float e2 = (kneeStart < 1 && e1 > kneeStart) ? eetfSpline(e1, kneeStart, maxLum) : e1;

		// following code line is not like in the report
		// there the tapring factor is b*(1-e2)^4. this is however incorrect, 
		// since this increases the brightness also in the high end above KS.
		// using e1 instead of e2 works only together with the change above
		float e3 = e2 + b * powf(1 - e1, 4);
		
		float e4 = e3 * (masterMaxEp - masterMinEp) + masterMinEp;
		e4 = std::clamp(e4, 0.0f, 1.0f);
		if (normalizeOutput) {
			e4 /= targetMaxEp;
		}
		uint16_t outSignal = e4 * (LUT_SIZE - 1) + 0.5f;
		lut[inSignal] = outSignal;
	}
}

