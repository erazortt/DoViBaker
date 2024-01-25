#pragma once
#include "avisynth.h"
#include "DoViTransferFunctions.h"

template<int bitDepth>
class DoViTonemap : public DoViTransferFunctions, public GenericVideoFilter
{
public:
  DoViTonemap(
    PClip child,
    float targetMaxNits,
    float targetMinNits,
    float masterMaxNits,
    float masterMinNits,
    float scale,
    IScriptEnvironment* env);

  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;
  static inline uint16_t signal2pq(uint16_t signal);
  static inline uint16_t pq2signal(uint16_t pq);

private:
  void applyTonemapRGB(PVideoFrame& dst, const PVideoFrame& src) const;
  //void applyTonemapYUV(PVideoFrame& dst, const PVideoFrame& src) const;

  bool dynamicMasterMaxPq;
  bool dynamicMasterMinPq;
  bool dynamicLumScale;
};

template<int bitDepth>
uint16_t DoViTonemap<bitDepth>::signal2pq(uint16_t signal) {
  const uint_fast32_t s = signal;
  if (bitDepth == 8) {
    return (s * 4111 + 135) >> 8;
  } 
  else if (bitDepth == 10) {
    return (s * 4099 + 513) >> 10;
  } 
  else if (bitDepth == 12){
    return signal;
  }
  else if (bitDepth == 14) {
    return (s * 16381 + 32768) >> 16;
  }
  else if (bitDepth == 16) {
    return (s * 4095 + 34815) >> 16;
  }
}

template<int bitDepth>
uint16_t DoViTonemap<bitDepth>::pq2signal(uint16_t pq) {
  const uint_fast32_t p = pq;
  if (bitDepth == 8) {
    return (p * 4081 + 32768) >> 16;
  }
  else if (bitDepth == 10) {
    return (p * 4093 + 8192) >> 14;
  }
  else if (bitDepth == 12) {
    return pq;
  }
  else if (bitDepth == 14) {
    return (p * 16387 + 2049) >> 12;
  }
  else if (bitDepth == 16) {
    return (p * 65551 + 2055) >> 12;
  }
}