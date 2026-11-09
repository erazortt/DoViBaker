#pragma once

#include <array>
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
    IScriptEnvironment* env);
  virtual ~DoViBaker();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;
  int __stdcall SetCacheHints(int cachehints, int frame_range) override
  {
      return cachehints == CACHE_GET_MTMODE ? MT_MULTI_INSTANCE : 0;
  }

private:
  /*
  * these upsampling functions are not following the paper, but should be correct when assuming top-left chroma location
  */
  static inline constexpr uint16_t upsampleLumaEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleLumaOdd(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleChromaEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleChromaOdd(const uint16_t* srcSamples, int idx0);

  /*
  * these are the original upsampling functions from the paper and seem to assume center-left chroma location which is actually incorrect for hdr sources
  static inline constexpr uint16_t upsampleHorzEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleHorzOdd(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleBlVertEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleBlVertOdd(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElYvertEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElYvertOdd(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElUVvertEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElUVvertOdd(const uint16_t* srcSamples, int idx0);
  */

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
  const bool blClipChromaSubSampled;
  const bool elClipChromaSubSampled;
};

/*
*  these upsampling functions are not following the paper, but should be correct when assuming top-left chroma location
*/
template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleChromaEven(const uint16_t* y, int n)
{
  return y[n];
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleChromaOdd(const uint16_t* y, int n)
{
  // this is spline16 at positions 0.5, 1.5
  auto val = (-307 * y[n - 1] + 2355 * y[n] + 2355 * y[n + 1] - 307 * y[n + 2] + 2048) >> 12;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleLumaEven(const uint16_t* y, int n)
{
  // this matches quite exactly spline16 at positions -1.75, -0.75, 0.25, 1.25
  auto val = (-3 * y[n - 2] + 29 * y[n - 1] + 111 * y[n] - 9 * y[n + 1] + 64) >> 7;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleLumaOdd(const uint16_t* y, int n)
{
  // this matches quite exactly spline16 at positions -1.25, -0.25, 0.75, 1.75
  auto val = (-9 * y[n - 1] + 111 * y[n] + 29 * y[n + 1] - 3 * y[n + 2] + 64) >> 7;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

/*
* these are the original upsampling functions from the paper and seem to assume center-left chroma location which is actually incorrect for hdr sources

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleHorzEven(const uint16_t* y, int n)
{
  return y[n];
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleHorzOdd(const uint16_t* y, int n)
{
  // this matches quite exactly spline64 at positions 0.5, 1.5, 2.5, 3.5
  auto val = (22 * y[n - 3] + 94 * y[n - 2] - 524 * y[n - 1] + 2456 * y[n] + 2456 * y[n + 1] - 524 * y[n + 2] +
    94 * y[n + 3] + 22 * y[n + 4] + 2048) >> 12;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleBlVertEven(const uint16_t* y, int n)
{
  // this matches quite exactly spline36 at positions -2.75, -1.75, -0.75, 0.25, 1.25, 2.25
  auto val = (2 * y[n - 3] - 12 * y[n - 2] + 65 * y[n - 1] + 222 * y[n] - 25 * y[n + 1] + 4 * y[n + 2] + 128) >> 8;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleBlVertOdd(const uint16_t* y, int n)
{
  // this matches quite exactly spline36 at positions -2.25, -1.25, -0.25, 0.75, 1.75, 2.75
  auto val = (4 * y[n - 2] - 25 * y[n - 1] + 222 * y[n] + 65 * y[n + 1] - 12 * y[n + 2] + 2 * y[n + 3] + 128) >> 8;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleElYvertEven(const uint16_t* y, int n)
{
  // this matches quite exactly spline16 at positions -1.75, -0.75, 0.25, 1.25
  auto val = (-3 * y[n - 2] + 29 * y[n - 1] + 111 * y[n] - 9 * y[n + 1] + 64) >> 7;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleElYvertOdd(const uint16_t* y, int n)
{
  // this matches quite exactly spline16 at positions -1.25, -0.25, 0.75, 1.75
  auto val = (-9 * y[n - 1] + 111 * y[n] + 29 * y[n + 1] - 3 * y[n + 2] + 64) >> 7;
  return static_cast<uint16_t>(std::clamp(val, 0, 0xFFFF));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleElUVvertEven(const uint16_t* y, int n)
{
  auto val = (64 * y[n - 1] + 192 * y[n] + 128) >> 8;
  return static_cast<uint16_t>(min(0xFFFF, val));
}

template<int quarterResolutionEl>
constexpr uint16_t DoViBaker<quarterResolutionEl>::upsampleElUVvertOdd(const uint16_t* y, int n)
{
  auto val = (192 * y[n] + 64 * y[n + 1] + 128) >> 8;
  return static_cast<uint16_t>(min(0xFFFF, val));
}
*/
