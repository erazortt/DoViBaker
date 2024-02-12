#pragma once

#include <string>
#include "DoViProcessor.h"

template<int quarterResolutionEl>
class DoViBaker : public GenericVideoFilter
{
public:
  DoViBaker(
    PClip blChild, 
    PClip elChild, 
    const char* rpuPath, 
    bool blClipChromaSubSampled, 
    bool elClipChromaSubSampled, 
    uint16_t desiredTrimPq,
    float targetMinLum,
    float targetMaxLum,
    bool qnd, 
    bool rgbProof, 
    bool nlqProof,
    int sourceProfile,
    bool outYUV,
    IScriptEnvironment* env);
  virtual ~DoViBaker();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;
  int __stdcall SetCacheHints(int cachehints, int frame_range) override
  {
      return cachehints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
  }

private:
  void upscaleEl(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  void upsampleChroma(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  //void upsampleElChroma(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  //void upsampleBlChroma(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);

  template<int blChromaSubsampling, int elChromaSubsampling>
  void doAllQuickAndDirty(PVideoFrame& rgb, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env) const;

  template<int chromaSubsampling>
  void applyDovi(PVideoFrame& dst, const PVideoFrame& blSrcY, const PVideoFrame& blSrcUV, const PVideoFrame& elSrcY, const PVideoFrame& elSrcUV, IScriptEnvironment* env) const;
  void convert2rgb(PVideoFrame& rgb, const PVideoFrame& y, const PVideoFrame& uv) const;
  void applyTrim(PVideoFrame& dst, const PVideoFrame& src) const;

  typedef uint16_t(*upscaler_t)(const uint16_t* srcSamples, int idx0);
  template<int vertLen, int nD>
  void upsampleVert(PVideoFrame& dst, const PVideoFrame& src, int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env);
  template<int vertLen, int nD>
  void upsampleHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env);
  //void upsampleHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, IScriptEnvironment* env);

  PClip elChild;
  DoViProcessor* doviProc;
  
  const bool qnd;
  const bool outYUV;
  const bool blClipChromaSubSampled;
  const bool elClipChromaSubSampled;
};
