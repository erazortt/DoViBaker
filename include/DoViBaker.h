
#ifndef __FILTERED_EWA_RESIZE_H
#define __FILTERED_EWA_RESIZE_H

#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

//#include <stdint.h>
#include "DoViProcessor.h"

template<bool chromaSubsampling>
class DoViBaker : public GenericVideoFilter
{
public:
  DoViBaker(PClip _blChild, PClip _elChild, const char* rpuPath, bool halfResolutionEl, IScriptEnvironment* env);
  virtual ~DoViBaker();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  void upsampleEl(PVideoFrame& dst, const PVideoFrame& el, IScriptEnvironment* env);
  void applyDovi(PVideoFrame& dst, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env);

  template<int vertLen>
  void upsampleElVert(PVideoFrame& dst, const PVideoFrame& src, const std::array<int, vertLen>& pD, const std::array<int, vertLen>& nD, int plane, IScriptEnvironment* env);
  void upsampleElHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, IScriptEnvironment* env);

  PClip elChild;
  int CPU_FLAG;
  DoViProcessor* doviProc;
  bool halfResolutionEl;
};

#endif
