#pragma once
#include <cstdint>

template<int signalBitDepth>
class DoViEetf
{
public:
  DoViEetf(float kneeOffset, bool normalizeOutput);

  inline uint16_t applyEETF(uint16_t s) const { return lut[s]; };
  void generateEETF(
    uint16_t targetMaxPq,
    uint16_t targetMinPq,
    uint16_t masterMaxPq,
    uint16_t masterMinPq,
    float lumScale,
    bool limitedInput);

private:
  static inline constexpr float eetfSpline(float e1, float kS, float maxLum);

  static constexpr int LUT_SIZE = 1 << signalBitDepth;
  const float kneeOffset;
  const bool normalizeOutput;
  uint16_t lut[LUT_SIZE];
};

template<int signalBitDepth>
constexpr float DoViEetf<signalBitDepth>::eetfSpline(float e1, float kS, float maxLum)
{
  float t = (e1 - kS) / (1 - kS);
  float p = ((2*t-3)*t*t+1)*kS + (((t-2)*t+1)*(1-kS) + (-2*t+3)*t*maxLum) * t;
  return p;
}
