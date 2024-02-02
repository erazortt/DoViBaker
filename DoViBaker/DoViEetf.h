#pragma once
#include <algorithm>

template<int signalBitDepth>
class DoViEetf
{
public:
  DoViEetf(bool normalizeOutput);

  inline uint16_t applyEETF(uint16_t s) const { return lut[s]; };
  void generateEETF(
    uint16_t targetMaxPq,
    uint16_t targetMinPq,
    uint16_t masterMaxPq,
    uint16_t masterMinPq,
    float lumScale);

private:
  static inline constexpr float eetfSpline(float e1, float KS, float maxLum);

  // the LUT becomes unnecessarily big for higher bit depths, cut at 12 bits
  static constexpr int LUT_BITS = signalBitDepth;
  static constexpr int LUT_SIZE = 1 << LUT_BITS;

  const bool normalizeOutput;
  uint16_t lut[LUT_SIZE];
};

template<int signalBitDepth>
constexpr float DoViEetf<signalBitDepth>::eetfSpline(float e1, float KS, float maxLum)
{
  float t = (e1 - KS) / (1 - KS);
  float p = ((2*t-3)*t*t+1)*KS + (((t-2)*t+1)*(1-KS) + (-2*t+3)*t*maxLum) * t;
  return p;
}
