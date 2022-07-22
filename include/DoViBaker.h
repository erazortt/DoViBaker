
#ifndef __FILTERED_EWA_RESIZE_H
#define __FILTERED_EWA_RESIZE_H

#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

//#include <stdint.h>
#include "DoViProcessor.h"

template<bool chromaSubsampling, bool quarterResolutionEl>
class DoViBaker : public GenericVideoFilter
{
public:
  DoViBaker(PClip _blChild, PClip _elChild, const char* rpuPath, bool qnd, IScriptEnvironment* env);
  virtual ~DoViBaker();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  void upsampleEl(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  void to444(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  void applyDovi(PVideoFrame& dst, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env);
  void convert2rgb(PVideoFrame& rgb, const PVideoFrame& y, const PVideoFrame& uv);
  void doAllQuickAndDirty(PVideoFrame& rgb, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env);

  typedef uint16_t(*upscaler_t)(const uint16_t* srcSamples, int idx0);
  template<int vertLen, int nD>
  void upsampleVert(PVideoFrame& dst, const PVideoFrame& src, int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env);
  void upsampleHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, IScriptEnvironment* env);

  PClip elChild;
  int CPU_FLAG;
  DoViProcessor* doviProc;
  bool qnd;
};

#endif
