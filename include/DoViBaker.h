#pragma once
#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

//#include <stdint.h>
#include "DoViProcessor.h"
#include "lut.h"

template<bool chromaSubsampling, bool quarterResolutionEl>
class DoViBaker : public GenericVideoFilter
{
public:
  DoViBaker(PClip _blChild, PClip _elChild, const char* rpuPath, bool qnd, std::vector<std::pair<uint16_t, std::string>> &cubes, IScriptEnvironment* env);
  virtual ~DoViBaker();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  void upscaleEl(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  void upsample(PVideoFrame& dst, const PVideoFrame& el, VideoInfo dstVi, IScriptEnvironment* env);
  void applyDovi(PVideoFrame& dst, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env);
  void applyLut(PVideoFrame& dst, const PVideoFrame& src);

  inline float toFloat(uint16_t in) const { return float(in) / ((1UL << 16) - 1); }
  inline uint16_t fromFloat(float in) const { return in * ((1UL << 16) - 1); }
  inline void sample2rgb(uint16_t& r, uint16_t& g, uint16_t& b, const uint16_t& y, const uint16_t& u, const uint16_t& v) const;
  void convert2rgb(PVideoFrame& rgb, const PVideoFrame& y, const PVideoFrame& uv) const;
  void doAllQuickAndDirty(PVideoFrame& rgb, const PVideoFrame& blSrc, const PVideoFrame& elSrc, IScriptEnvironment* env) const;

  typedef uint16_t(*upscaler_t)(const uint16_t* srcSamples, int idx0);
  template<int vertLen, int nD>
  void upsampleVert(PVideoFrame& dst, const PVideoFrame& src, int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env);
  void upsampleHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, IScriptEnvironment* env);

  PClip elChild;
  int CPU_FLAG;
  DoViProcessor* doviProc;
  const bool qnd;
  std::vector<std::pair<uint16_t, std::unique_ptr<timecube::Lut>>> luts;
  const timecube::Lut* current_frame_lut;
};
