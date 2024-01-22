#pragma once
#include "avisynth.h"
#include "DoViTonemap.h"

template<int bitDepth>
class DoViTonemapper : public DoViTonemap, public GenericVideoFilter
{
public:
  DoViTonemapper(
    PClip child,
    float masterMaxNits,
    float masterMinNits,
    float targetMaxNits,
    float targetMinNits,
    float scale,
    IScriptEnvironment* env);

  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;
  static inline uint16_t signal2pq(uint16_t signal);
  static inline uint16_t pq2signal(uint16_t pq);

private:
  void applyTonemapRGB(PVideoFrame& dst, const PVideoFrame& src) const;
  //void applyTonemapYUV(PVideoFrame& dst, const PVideoFrame& src) const;
};

template<int bitDepth>
uint16_t DoViTonemapper<bitDepth>::signal2pq(uint16_t signal) {
  if (bitDepth < 12) {
    return signal << (12 - bitDepth);
  } else {
    return signal >> (bitDepth - 12);
  }
}

template<int bitDepth>
uint16_t DoViTonemapper<bitDepth>::pq2signal(uint16_t pq) {
  if (bitDepth < 12) {
    return pq >> (12 - bitDepth);
  } else {
    return pq << (bitDepth - 12);
  }
}