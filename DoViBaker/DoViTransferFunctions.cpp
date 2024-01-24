#include "DoViTransferFunctions.h"

DoViTransferFunctions::DoViTransferFunctions(
	float targetMaxNits_,
	float targetMinNits_,
	float masterMaxNits_,
	float masterMinNits_,
	float lumScale_)
	: masterMaxPq(nits2pq(masterMaxNits_))
	, masterMinPq(nits2pq(masterMinNits_))
	, targetMaxPq(nits2pq(targetMaxNits_))
	, targetMinPq(nits2pq(targetMinNits_))
	, lumScale(lumScale_)
{
		generateEETF();
}

void DoViTransferFunctions::generateEETF()
{
	// based on the report ITU-R BT.2408-7 Annex 5
	float masterMaxEp = EOTFinv(EOTF(masterMaxPq / 4095.0) * lumScale);
	float masterMinEp = EOTFinv(EOTF(masterMinPq / 4095.0) * lumScale);
	float targetMaxEp = targetMaxPq / 4095.0;
	float targetMinEp = targetMinPq / 4095.0;

	// this is not from the report
	// any mapping is unnecessary while inside the range of tragets capabilities
	masterMaxEp = fmaxf(masterMaxEp, targetMaxEp);
	masterMinEp = fminf(masterMinEp, targetMinEp);

	float maxLum = (targetMaxEp - masterMinEp) / (masterMaxEp - masterMinEp);
	float minLum = (targetMinEp - masterMinEp) / (masterMaxEp - masterMinEp);
	float KS = 1.5 * maxLum - 0.5;
	float b = minLum;

	for (uint16_t inPq = 0; inPq < LUT_SIZE; inPq++) {
		float ep = inPq / float(LUT_SIZE - 1);
		float Y = EOTF(ep) * lumScale;
		ep = EOTFinv(Y);
		float e1 = (ep - masterMinEp) / (masterMaxEp - masterMinEp);

		// this following claming is not from the report
		// it serves to stop using the spline above where it is supposed to be used 
		// this works only in conjunction with the change above
		e1 = std::clamp(e1, 0.0f, 1.0f); 

		float e2 = (KS < 1 && e1 > KS) ? eetfSpline(e1, KS, maxLum) : e1;

		// following code line is not like in the report
		// there the tapring factor is b*(1-e2)^4. this is however incorrect, 
		// since this increases the brightness also in the high end above KS.
		// using e1 instead of e2 works only together with the change above
		float e3 = e2 + b * powf(1 - e1, 4); 
		
		float e4 = e3 * (masterMaxEp - masterMinEp) + masterMinEp;
		e4 = std::clamp(e4, 0.0f, 1.0f);
		uint16_t outPq = e4 * (LUT_SIZE - 1) + 0.5;
		lut[inPq] = outPq;
	}
}
