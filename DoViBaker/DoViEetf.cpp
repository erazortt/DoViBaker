#include "DoViEetf.h"
#include "DoViProcessor.h"
#include <algorithm>

template<int signalBitDepth>
DoViEetf<signalBitDepth>::DoViEetf(float kneeOffset_, bool normalizeOutput_)
	: kneeOffset(kneeOffset_),
		normalizeOutput(normalizeOutput_) {}

// explicitly instantiate the template for the linker
template class DoViEetf<8>;
template class DoViEetf<10>;
template class DoViEetf<12>;
template class DoViEetf<14>;
template class DoViEetf<16>;

template<int signalBitDepth>
void DoViEetf<signalBitDepth>::generateEETF(
	uint16_t targetMaxPq,
	uint16_t targetMinPq,
	uint16_t masterMaxPq,
	uint16_t masterMinPq,
	float lumScale,
	bool limitedInput)
{
	// based on report ITU-R BT.2408-7 Annex 5 (was in ITU-R BT.2390 until revision 7)
	float masterMaxEp = DoViProcessor::EOTFinv(DoViProcessor::EOTF(masterMaxPq / 4095.0f) * lumScale);
	float masterMinEp = DoViProcessor::EOTFinv(DoViProcessor::EOTF(masterMinPq / 4095.0f) * lumScale);
	const float targetMaxEp = targetMaxPq / 4095.0f;
	const float targetMinEp = targetMinPq / 4095.0f;

	// This is not from the report.
	// Any mapping is unnecessary while inside the range of tragets capabilities
	masterMaxEp = fmaxf(masterMaxEp, targetMaxEp);
	masterMinEp = fminf(masterMinEp, targetMinEp);

	const float maxLum = (targetMaxEp - masterMinEp) / (masterMaxEp - masterMinEp);
	const float minLum = (targetMinEp - masterMinEp) / (masterMaxEp - masterMinEp);
	const float kneeStart = (1 + kneeOffset) * maxLum - kneeOffset;
	const float b = minLum;

	// we don't need unnecessarily big values of the taper power, thus b is limited
	const float taperPwr = 1 / std::max(b, 0.01f);

	int inBlack = 0;
	int inWhite = LUT_SIZE - 1;
	if (limitedInput) {
		inBlack = 16 << (signalBitDepth - 8);
		inWhite = (235 << (signalBitDepth - 8));
	}
	int inMinSignal = masterMinEp * (inWhite - inBlack) + inBlack + 0.5f;
	int inMaxSignal = masterMaxEp * (inWhite - inBlack) + inBlack + 0.5f;

	// limiting is needed since masterMaxEp can be above 1 due to lumScale
	inMinSignal = std::min(inMinSignal, LUT_SIZE);
	inMaxSignal = std::min(inMaxSignal, LUT_SIZE);

	uint16_t outBlack = 0;
	uint16_t outWhite = LUT_SIZE - 1;
	if (!normalizeOutput) {
		outBlack = targetMinEp * (LUT_SIZE - 1) + 0.5f;
		outWhite = targetMaxEp * (LUT_SIZE - 1) + 0.5f;
	}

	int inSignal = 0;
	for (; inSignal < inMinSignal; inSignal++) {
		// skip unnecessary processing where result is known
		lut[inSignal] = outBlack;
	}
	for (; inSignal < inMaxSignal; inSignal++) {
		float ep = (inSignal - inBlack) / float(inWhite - inBlack);
		float e = DoViProcessor::EOTF(ep) * lumScale;
		ep = DoViProcessor::EOTFinv(e);
		float e1 = (ep - masterMinEp) / (masterMaxEp - masterMinEp);

		// This following clamping is not from the report.
		// It serves to stop using the spline above where it is supposed to be used.
		// This works only in conjunction with the change above
		e1 = std::clamp(e1, 0.0f, 1.0f);

		const float e2 = (e1 > kneeStart) ? eetfSpline(e1, kneeStart, maxLum) : e1;

		// Following code line is not like in the report.
		// There the tapering factor is b*(1-e2)^4. This is however incorrect, 
		// since the dependence on e2 increases the brightness also in the high-end.
		// Using e1 instead of e2 works only together with all the changes above.
		// Additionally the fixed power of 4 is unnecessarily low and thus increases the 
		// low-end brightness too much. It can be set to 1/b, so long targetMinEp/targetMaxEp<0.5.
		const float taper = b * powf(1 - e1, taperPwr);
		const float e3 = e2 + taper;
		
		float e4 = e3 * (masterMaxEp - masterMinEp) + masterMinEp;

		// Following code line is not like in the report.
		// This clamping ensures that even in unexpected situations the result is not outside the value range.
		// (like too low kneeOffsets and/or taperPowers and/or too high luminosity factors)
		e4 = std::clamp(e4, targetMinEp, targetMaxEp);

		if (normalizeOutput) {
			e4 = (e4 - targetMinEp) / (targetMaxEp - targetMinEp);
		}
		uint16_t outSignal = e4 * (LUT_SIZE - 1) + 0.5f;
		lut[inSignal] = outSignal;
	}
	for (; inSignal < LUT_SIZE; inSignal++) {
		// skip unnecessary processing where result is known
		lut[inSignal] = outWhite;
	}
}

template void DoViEetf<8>::generateEETF(uint16_t targetMaxPq, uint16_t targetMinPq, uint16_t masterMaxPq, uint16_t masterMinPq, float lumScale, bool limitedInput);
template void DoViEetf<10>::generateEETF(uint16_t targetMaxPq, uint16_t targetMinPq, uint16_t masterMaxPq, uint16_t masterMinPq, float lumScale, bool limitedInput);
template void DoViEetf<12>::generateEETF(uint16_t targetMaxPq, uint16_t targetMinPq, uint16_t masterMaxPq, uint16_t masterMinPq, float lumScale, bool limitedInput);
template void DoViEetf<14>::generateEETF(uint16_t targetMaxPq, uint16_t targetMinPq, uint16_t masterMaxPq, uint16_t masterMinPq, float lumScale, bool limitedInput);
template void DoViEetf<16>::generateEETF(uint16_t targetMaxPq, uint16_t targetMinPq, uint16_t masterMaxPq, uint16_t masterMinPq, float lumScale, bool limitedInput);
