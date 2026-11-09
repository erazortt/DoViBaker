#pragma once

#include <algorithm>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <climits>
#include <cmath>
#endif

// powf/fminf/fmaxf are used unqualified, the C header declares them globally.
#include <math.h>

#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)
#include "dovi/rpu_parser.h"

class DoViProcessor {
public:
  DoViProcessor(const char* rpuPath, IScriptEnvironment* env, uint8_t blContainerBits, uint8_t elContainerBits);
  virtual ~DoViProcessor();
  bool wasCreationSuccessful() const { return successfulCreation; }
  void setRgbProof(bool set = true) { rgbProof = set; }
  void setNlqProof(bool set = true) { nlqProof = set; }
  inline void setTrim(uint16_t trimPq, float targetMinNits, float targetMaxNits);

  bool intializeFrame(int frame, IScriptEnvironment* env, const uint8_t* rpubuf, size_t rpusize);
  inline int getClipLength() const { return rpus->len; }
  inline bool isIntegratedRpu() const { return !rpus; }
  inline bool isSceneChange() const { return scene_refresh_flag; }
  inline bool isLimitedRangeOutput() const { return !signal_full_range_flag; }
  inline bool elProcessingEnabled() const { return !disable_residual_flag; }
  inline bool trimProcessingEnabled() const { return !skipTrim; }
  inline uint16_t getNlqOffset(int cmp) const { return nlq_offset[cmp] << (outContainerBitDepth - el_bit_depth); }
  inline uint16_t getDynamicMinPq() const { return dynamic_min_pq; }
  inline uint16_t getDynamicMaxPq() const { return dynamic_max_pq; }
  inline uint16_t getDynamicMaxContentLightLevel() const { return dynamic_max_content_light_level; }
  inline uint16_t getStaticMaxPq() const { return static_max_pq; }
  inline uint16_t getStaticMaxContentLightLevel() const { return static_max_content_light_level; }
  inline uint16_t getStaticMaxAvgContentLightLevel() const { return static_max_avg_content_light_level; }
  inline uint16_t getStaticMasterDisplayMaxLuminance() const { return static_master_display_max_luminance; }
  inline uint16_t getStaticMasterDisplayMinLuminance() const { return static_master_display_min_luminance; }
  const std::vector<uint16_t>& getAvailableTrimPqs() const { return availableTrimPqs; }

  static inline float EOTF(float ep);
  static inline float EOTFinv(float fd);
  static inline float pq2nits(uint16_t pq);
  static inline uint16_t nits2pq(float nits);

  inline uint16_t processSampleY(uint16_t bl, uint16_t el) const;
  inline uint16_t processSampleU(uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const;
  inline uint16_t processSampleV(uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const;

  inline void sample2rgb(uint16_t& r, uint16_t& g, uint16_t& b, const uint16_t& y, const uint16_t& u, const uint16_t& v) const;
  void processTrim(uint16_t& ro, uint16_t& go, uint16_t& bo, const uint16_t& ri, const uint16_t& gi, const uint16_t& bi) const;

  static constexpr uint8_t outContainerBitDepth = 16;
private:
  static inline constexpr uint16_t ClipUint16(int value);
  void showMessage(const char* message, IScriptEnvironment* env);
  uint16_t processSample(int cmp, uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const;
  int getPivotIndex(int cmp, uint16_t sample) const;
  uint16_t polynompialMapping(int cmp, int pivot_idx, uint16_t sample) const;
  uint16_t mmrMapping(int cmp, int pivot_idx, uint64_t sampleY, uint64_t sampleU, uint64_t sampleV) const;
  int16_t nonLinearInverseQuantization(int cmp, uint16_t sample) const;
  uint16_t signalReconstruction(uint16_t v, int16_t r) const;
  void prepareTrimCoef();

  static constexpr float m1 = 2610.0 / 4096 / 4;
  static constexpr float m2 = 2523.0 / 4096 * 128;
  static constexpr float c3 = 2392.0 / 4096 * 32;
  static constexpr float c2 = 2413.0 / 4096 * 32;
  static constexpr float c1 = c3 - c2 + 1;

  const DoviRpuOpaqueList* rpus;

  bool successfulCreation;
  bool rgbProof;
  bool nlqProof;
  bool skipTrim;
  bool trimInfoMissing;

  const uint8_t blContainerBitDepth;
  const uint8_t elContainerBitDepth;
  uint8_t bl_bit_depth;
  uint8_t el_bit_depth;
  uint8_t out_bit_depth;
  uint8_t coeff_log2_denom;
  bool is_fel;
  bool disable_residual_flag;
  bool scene_refresh_flag;
  bool signal_full_range_flag;

  uint16_t dynamic_max_pq;
  uint16_t dynamic_min_pq;
  uint16_t dynamic_avg_pq;
  uint16_t dynamic_max_content_light_level;
  uint16_t static_max_pq;
  uint16_t static_max_content_light_level;
  uint16_t static_max_avg_content_light_level;
  uint16_t static_master_display_max_luminance;
  uint16_t static_master_display_min_luminance;
  int16_t ycc_to_rgb_coef[9];
  uint32_t ycc_to_rgb_offset[3];

  static const uint16_t ycc_to_rgb_coef_scale_shifts = 13;
  static const uint16_t ycc_to_rgb_offset_scale_shifts = (28-outContainerBitDepth);
  static const uint16_t rgb_to_lms_coef_scale_shifts = 14;

  uint8_t num_pivots_minus1[3];
  std::vector<std::vector<uint16_t>> pivot_value;

  std::vector<std::vector<uint8_t>> mapping_idc;
  std::vector<std::vector<uint8_t>> poly_order;
  std::vector<std::vector<std::vector<int32_t>>> fp_poly_coef;

  std::vector<std::vector<uint8_t>> mmr_order;
  std::vector<std::vector<int64_t>> fp_mmr_const;
  std::vector<std::vector<std::vector<std::vector<int64_t>>>> fp_mmr_coef;

  uint16_t nlq_offset[3];
  uint32_t fp_hdr_in_max[3];
  uint32_t fp_linear_deadzone_slope[3];
  uint32_t fp_linear_deadzone_threshold[3];

  uint16_t desiredTrimPq;
  float targetMaxNits;
  float targetMinNits;
  std::vector<uint16_t> availableTrimPqs;
  struct TrimCoefficients {
    uint16_t slope;
    uint16_t offset;
    uint16_t power;
    uint16_t chroma_weight;
    uint16_t saturation_gain;
    uint16_t tone_detail;

    float maxNits;
    float minNits;
    float ccc[3];
    float goP[3];
    float cS[2];
  } trim;
};

void DoViProcessor::setTrim(uint16_t trimPq, float targetMinNits, float targetMaxNits)
{
  desiredTrimPq = trimPq;
  this->targetMinNits = std::clamp(targetMinNits, 0.00001f, 5.0f);
  this->targetMaxNits = std::clamp(targetMaxNits, 5.0f, 10000.0f);
}

float DoViProcessor::EOTF(float ep)
{
  const float epower = powf(ep, 1 / m2);
  const float num = std::max(epower - c1, 0.0f);
  const float denom = c2 - c3 * epower;
  return powf(num / denom, 1 / m1);
}

float DoViProcessor::EOTFinv(float Y)
{
  const float epower = powf(Y, m1);
  const float num = c1 + c2 * epower;
  const float denom = 1 + c3 * epower;
  return powf(num / denom, m2);
}

float DoViProcessor::pq2nits(uint16_t pq)
{
  const float ep = pq / 4095.0f;
  return EOTF(ep) * 10000;
}

uint16_t DoViProcessor::nits2pq(float nits)
{
  const float Y = nits / 10000;
  return EOTFinv(Y) * 4095.0f + 0.5f;
}

constexpr uint16_t DoViProcessor::ClipUint16(int value)
{
  return static_cast<uint16_t>(std::clamp(value, 0, 0xFFFF));
}

uint16_t DoViProcessor::processSampleY(uint16_t bl, uint16_t el) const {
  return processSample(0, bl, el, 0, 0, 0);
}

uint16_t DoViProcessor::processSampleU(uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const {
  return processSample(1, bl, el, mmrBlY, mmrBlU, mmrBlV);
}

uint16_t DoViProcessor::processSampleV(uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const {
  return processSample(2, bl, el, mmrBlY, mmrBlU, mmrBlV);
}

void DoViProcessor::sample2rgb(uint16_t& r, uint16_t& g, uint16_t& b, const uint16_t& y, const uint16_t& u, const uint16_t& v) const
{
  int yf = y - ycc_to_rgb_offset[0];
  int uf = u - ycc_to_rgb_offset[1];
  int vf = v - ycc_to_rgb_offset[2];
  r = ClipUint16((ycc_to_rgb_coef[0] * yf + ycc_to_rgb_coef[1] * uf + ycc_to_rgb_coef[2] * vf) >> ycc_to_rgb_coef_scale_shifts);
  g = ClipUint16((ycc_to_rgb_coef[3] * yf + ycc_to_rgb_coef[4] * uf + ycc_to_rgb_coef[5] * vf) >> ycc_to_rgb_coef_scale_shifts);
  b = ClipUint16((ycc_to_rgb_coef[6] * yf + ycc_to_rgb_coef[7] * uf + ycc_to_rgb_coef[8] * vf) >> ycc_to_rgb_coef_scale_shifts);
}


// see also:
// https://code.videolan.org/videolan/libplacebo/-/blob/775a9325a23e26443b562b104c1fe949b99aa3c8/src/colorspace.c
// https://github.com/test-full-band/tfb-video/blob/master/core/src/main/java/band/full/video/dolby/VdrDmDataPayload.java
// https://ffmpeg.org/doxygen/trunk/dovi__meta_8h_source.html
// https://ffmpeg.org/doxygen/trunk/dovi__rpu_8c_source.html
