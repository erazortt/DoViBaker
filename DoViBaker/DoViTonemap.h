#pragma once
#include <algorithm>

class DoViTonemap
{
public:
  DoViTonemap(
    float targetMaxNits, 
    float targetMinNits,
    float masterMaxNits,
    float masterMinNits,
    float scale);

  static inline float EOTF(float ep);
  static inline float EOTFinv(float fd);
  static inline float pq2nits(uint16_t pq);
  static inline uint16_t nits2pq(float nits);
  static inline float EETF(float e1, float KS, float maxLum);
  inline uint16_t applyLut(uint16_t pq) const { return lut[pq]; };

protected:
  void generateLut();
  const uint16_t targetMaxPq;
  const uint16_t targetMinPq;
  uint16_t masterMaxPq;
  uint16_t masterMinPq;
  float lumScale;

private:
  static constexpr float m1 = 2610.0 / 4096 / 4;
  static constexpr float m2 = 2523.0 / 4096 * 128;
  static constexpr float c3 = 2392.0 / 4096 * 32;
  static constexpr float c2 = 2413.0 / 4096 * 32;
  static constexpr float c1 = c3 - c2 + 1;
  static constexpr int LUT_SIZE = 4096;

  uint16_t lut[LUT_SIZE];
};

float DoViTonemap::EOTF(float ep)
{
  const float epower = powf(ep, 1 / m2);
  const float num = std::max(epower - c1, 0.0f);
  const float denom = c2 - c3 * epower;
  return powf(num / denom, 1 / m1);
}

float DoViTonemap::EOTFinv(float Y)
{
  const float epower = powf(Y, m1);
  const float num = c1 + c2 * epower;
  const float denom = 1 + c3 * epower;
  return powf(num / denom, m2);
}

float DoViTonemap::pq2nits(uint16_t pq)
{
  const float ep = pq / 4095.0;
  return DoViTonemap::EOTF(ep) * 10000;
}

uint16_t DoViTonemap::nits2pq(float nits)
{
  const float Y = nits / 10000;
  return DoViTonemap::EOTFinv(Y) * 4095.0 + 0.5;
}

float DoViTonemap::EETF(float e1, float KS, float maxLum)
{
  float t = (e1 - KS) / (1 - KS);
  float p = (2*t*t*t-3*t*t+1)*KS+(t*t*t-2*t*t+t)*(1-KS)+(-2*t*t*t+3*t*t)*maxLum;
  float e2 = (e1 < KS) ? e1 : p;
  return e2;
}
