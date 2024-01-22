#include "DoViTonemap.h"

DoViTonemap::DoViTonemap(
	float targetMaxNits_,
	float targetMinNits_,
	float masterMaxNits_,
	float masterMinNits_,
	float lumScale_)
	: staticMasterMaxPq(masterMaxNits_ > -1)
	, staticMasterMinPq(masterMinNits_ > -1)
	, staticLumScale(lumScale_ > -1)
	, masterMaxPq(-1)
	, masterMinPq(-1)
	, targetMaxPq(nits2pq(targetMaxNits_))
	, targetMinPq(nits2pq(targetMinNits_))
	, lumScale(-1)
{
	if (staticMasterMaxPq) masterMaxPq = nits2pq(masterMaxNits_);
	if (staticMasterMinPq) masterMinPq = nits2pq(masterMinNits_);
	if (staticLumScale) lumScale = lumScale_;
	if (staticMasterMaxPq && staticMasterMinPq && staticLumScale)
		generateLut();
}

void DoViTonemap::generateLut()
{
	// based on ITU-R BT.2408-7 Annex 5
	float masterMaxEp = EOTFinv(EOTF(masterMaxPq / 4095.0) * lumScale);
	//masterMaxEp = fmaxf(fminf(masterMaxEp, 1), 0);
	float masterMinEp = EOTFinv(EOTF(masterMinPq / 4095.0) * lumScale);
	//masterMinEp = fmaxf(fminf(masterMinEp, 1), 0);
	float targetMaxEp = targetMaxPq / 4095.0;
	float targetMinEp = targetMinPq / 4095.0;
	masterMaxEp = fmaxf(masterMaxEp, targetMaxEp);

	float maxLum = (targetMaxEp - masterMinEp) / (masterMaxEp - masterMinEp);
	float minLum = (targetMinEp - masterMinEp) / (masterMaxEp - masterMinEp);
	float KS = 1.5 * maxLum - 0.5;
	float b = minLum;

	for (uint16_t inPq = 0; inPq < LUT_SIZE; inPq++) {
		float ep = inPq / float(LUT_SIZE - 1);
		float Y = EOTF(ep) * lumScale;
		ep = EOTFinv(Y);
		float e1 = (ep - masterMinEp) / (masterMaxEp - masterMinEp);
		e1 = fmaxf(fminf(e1, 1), 0);
		float e2 = (KS < 1 && e1 > KS) ? EETF(e1, KS, maxLum) : e1;
		float e3 = e2 + b * powf(1 - e2, 4);
		float e4 = e3 * (masterMaxEp - masterMinEp) + masterMinEp;
		e4 = fmaxf(fminf(e4, 1), 0);
		uint16_t outPq = e4 * (LUT_SIZE - 1) + 0.5;
		lut[inPq] = outPq;
	}
}

