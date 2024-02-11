#pragma once
#include "avisynth.h"
#include "DoViEetf.h"

template<int signalBitDepth>
class DoViTonemap : public GenericVideoFilter
{
public:
  DoViTonemap(
    PClip child,
    float targetMaxNits,
    float targetMinNits,
    float masterMaxNits,
    float masterMinNits,
    float lumScale,
    float kneeOffset,
    bool normalizeOutput,
    IScriptEnvironment* env);

  virtual ~DoViTonemap();

  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  void applyTonemapRGB(PVideoFrame& dst, const PVideoFrame& src) const;
  //void applyTonemapYUV(PVideoFrame& dst, const PVideoFrame& src) const;

  const uint16_t targetMaxPq;
  const uint16_t targetMinPq;
  uint16_t masterMaxPq;
  uint16_t masterMinPq;
  float lumScale;

  const bool dynamicMasterMaxPq;
  const bool dynamicMasterMinPq;
  const bool dynamicLumScale;

  DoViEetf<signalBitDepth>* doviEetf;
};
