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
    float scale,
    bool normalize,
    IScriptEnvironment* env);

  virtual ~DoViTonemap();

  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  void applyTonemapRGB(PVideoFrame& dst, const PVideoFrame& src) const;
  //void applyTonemapYUV(PVideoFrame& dst, const PVideoFrame& src) const;

  uint16_t targetMaxPq;
  uint16_t targetMinPq;
  uint16_t masterMaxPq;
  uint16_t masterMinPq;
  float lumScale;

  bool dynamicMasterMaxPq;
  bool dynamicMasterMinPq;
  bool dynamicLumScale;

  DoViEetf<signalBitDepth>* doviEetf;
};
