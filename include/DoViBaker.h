#pragma once
#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

#include "DoViProcessor.h"
#include "lut.h"

#include <string>

template<int quarterResolutionEl>
class DoViBaker : public GenericVideoFilter
{
public:
  DoViBaker(
    PClip _blChild, 
    PClip _elChild, 
    const char* rpuPath, 
    bool blClipChromaSubSampled, 
    bool elClipChromaSubSampled, 
    std::vector<std::pair<uint16_t, std::string>> &cubes, 
    bool qnd, 
    bool rgbProof, 
    bool nlqProof, 
    IScriptEnvironment* env);
  virtual ~DoViBaker();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

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
  void applyLut(PVideoFrame& dst, const PVideoFrame& src) const;

  typedef uint16_t(*upscaler_t)(const uint16_t* srcSamples, int idx0);
  template<int vertLen, int nD>
  void upsampleVert(PVideoFrame& dst, const PVideoFrame& src, int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env);
  template<int vertLen, int nD>
  void upsampleHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, const std::array<int, vertLen>& Dn0p, const upscaler_t evenUpscaler, const upscaler_t oddUpscaler, IScriptEnvironment* env);
  //void upsampleHorz(PVideoFrame& dst, const PVideoFrame& src, int plane, IScriptEnvironment* env);

  PClip elChild;
  int CPU_FLAG;
  DoViProcessor* doviProc;
  const bool qnd;
  const bool blClipChromaSubSampled;
  const bool elClipChromaSubSampled;
  std::vector<std::pair<uint16_t, std::unique_ptr<timecube::Lut>>> luts;
  const timecube::Lut* current_frame_lut;
};
