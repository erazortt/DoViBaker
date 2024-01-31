#pragma once
#include <algorithm>

template<int signalBitDepth>
class DoViEetf
{
public:
  DoViEetf(bool normalizeOutput);

  inline uint16_t applyEETF(uint16_t s) const { return lut2signal(lut[signal2lut(s)]); };
  void generateEETF(
    uint16_t targetMaxPq,
    uint16_t targetMinPq,
    uint16_t masterMaxPq,
    uint16_t masterMinPq,
    float lumScale);

private:
  static inline float eetfSpline(float e1, float KS, float maxLum);
  static inline uint16_t signal2lut(uint16_t signal);
  static inline uint16_t lut2signal(uint16_t coord);

  // the LUT becomes unnecessarily big for higher bit depths, cut at 12 bits
  static constexpr int LUT_BITS = (signalBitDepth < 12) ? signalBitDepth : 12;
  static constexpr int LUT_SIZE = 1 << LUT_BITS;

  const bool normalizeOutput;
  uint16_t lut[LUT_SIZE];
};

template<int signalBitDepth>
float DoViEetf<signalBitDepth>::eetfSpline(float e1, float KS, float maxLum)
{
  float t = (e1 - KS) / (1 - KS);
  float p = (2*t*t*t-3*t*t+1)*KS+(t*t*t-2*t*t+t)*(1-KS)+(-2*t*t*t+3*t*t)*maxLum;
  float e2 = (e1 < KS) ? e1 : p;
  return e2;
}

template<int signalBitDepth>
uint16_t DoViEetf<signalBitDepth>::signal2lut(uint16_t signal) {
  const uint_fast32_t s = signal;
  if (signalBitDepth == 8) {
    return signal;
  }
  else if (signalBitDepth == 10) {
    return signal;
  }
  else if (signalBitDepth == 12) {
    return signal;
  }
  else if (signalBitDepth == 14) {
    return (s * 16381 + 0x8000) >> 16;
  }
  else if (signalBitDepth == 16) {
    return (s * 65521 + 0x80000) >> 20;
  }
}

template<int signalBitDepth>
uint16_t DoViEetf<signalBitDepth>::lut2signal(uint16_t coord) {
  const uint_fast32_t c = coord;
  if (signalBitDepth == 8) {
    return coord;
  }
  else if (signalBitDepth == 10) {
    return coord;
  }
  else if (signalBitDepth == 12) {
    return coord;
  }
  else if (signalBitDepth == 14) {
    return (c * 16387 + 2049) >> 12;
  }
  else if (signalBitDepth == 16) {
    return (c * 65551 + 2055) >> 12;
  }
}