#pragma once
#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

#include "rpu_parser.h"
#include <vector>
#include <windows.h>
#include <string>

typedef DoviRpuOpaqueList* (*f_dovi_parse_rpu_bin_file)(const char* path);
typedef void (*f_dovi_rpu_list_free)(DoviRpuOpaqueList* ptr);
typedef const char* (*f_dovi_rpu_get_error)(const DoviRpuOpaque* ptr);
typedef const DoviRpuDataHeader* (*f_dovi_rpu_get_header)(const DoviRpuOpaque* ptr);
typedef void (*f_dovi_rpu_free_header)(const DoviRpuDataHeader* ptr);
typedef const DoviRpuDataNlq* (*f_dovi_rpu_get_data_nlq)(const DoviRpuOpaque* ptr);
typedef void (*f_dovi_rpu_free_data_nlq)(const DoviRpuDataNlq* ptr);
typedef const DoviRpuDataMapping* (*f_dovi_rpu_get_data_mapping)(const DoviRpuOpaque* ptr);
typedef void (*f_dovi_rpu_free_data_mapping)(const DoviRpuDataMapping* ptr);
typedef const DoviVdrDmData* (*f_dovi_rpu_get_vdr_dm_data)(const DoviRpuOpaque* ptr);
typedef void (*f_dovi_rpu_free_vdr_dm_data)(const DoviVdrDmData* ptr);

class DoViProcessor {
public:
  DoViProcessor(const char* rpuPath, IScriptEnvironment* env);
  virtual ~DoViProcessor();
  void intializeFrame(int frame, IScriptEnvironment* env);
  uint16_t getMaxContentLightLevel() const { return max_content_light_level; }
  float getYcc2RgbCoef(int i) const { return ycc_to_rgb_coef[i]; }
  float getYcc2RgbOff(int i) const { return ycc_to_rgb_offset[i]; }

  static inline constexpr uint16_t Clip3(uint16_t lower, uint16_t upper, int value);
  static inline constexpr uint16_t upsampleHorzEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleHorzOdd(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElYvertEven(const uint16_t *srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElYvertOdd(const uint16_t *srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElUVvertEven(const uint16_t* srcSamples, int idx0);
  static inline constexpr uint16_t upsampleElUVvertOdd(const uint16_t* srcSamples, int idx0);

  inline uint16_t processSampleY(uint16_t bl, uint16_t el) const;
  inline uint16_t processSampleU(uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const;
  inline uint16_t processSampleV(uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const;


private:
  void showMessage(const char* message, IScriptEnvironment* env);
  uint16_t processSample(int cmp, uint16_t bl, uint16_t el, uint16_t mmrBlY, uint16_t mmrBlU, uint16_t mmrBlV) const;
  int getPivotIndex(int cmp, uint16_t sample) const;
  uint16_t polynompialMapping(int cmp, int pivot_idx, uint16_t sample) const;
  uint16_t mmrMapping(int cmp, int pivot_idx, uint16_t sampleY, uint16_t sampleU, uint16_t sampleV) const;
  int16_t nonLinearInverseQuantization(int cmp, uint16_t sample) const;
  uint16_t signalReconstruction(uint16_t v, int16_t r) const;

  HINSTANCE doviLib;
  DoviRpuOpaqueList* rpus;

  f_dovi_parse_rpu_bin_file dovi_parse_rpu_bin_file;
  f_dovi_rpu_list_free dovi_rpu_list_free;
  f_dovi_rpu_get_header dovi_rpu_get_header;
  f_dovi_rpu_free_header dovi_rpu_free_header;
  f_dovi_rpu_get_data_nlq dovi_rpu_get_data_nlq;
  f_dovi_rpu_free_data_nlq dovi_rpu_free_data_nlq;
  f_dovi_rpu_get_vdr_dm_data dovi_rpu_get_vdr_dm_data;
  f_dovi_rpu_free_vdr_dm_data dovi_rpu_free_vdr_dm_data;
  f_dovi_rpu_get_data_mapping dovi_rpu_get_data_mapping;
  f_dovi_rpu_free_data_mapping dovi_rpu_free_data_mapping;
  f_dovi_rpu_get_error dovi_rpu_get_error;

  uint8_t bl_bit_depth;
  uint8_t el_bit_depth;
  uint8_t out_bit_depth;
  uint8_t coeff_log2_denom;
  bool disable_residual_flag;
  uint8_t nlq_method_idc;

  uint16_t max_content_light_level;
  int16_t ycc_to_rgb_coef[8];
  uint32_t ycc_to_rgb_offset[3];

  uint8_t num_pivots_minus1[3];
  std::vector<std::vector<uint16_t>> pivot_value;

  std::vector<std::vector<uint8_t>> mapping_idc;
  std::vector<std::vector<uint8_t>> poly_order;
  std::vector<std::vector<std::vector<int32_t>>> fp_poly_coef;

  std::vector<std::vector<uint8_t>> mmr_order;
  std::vector<std::vector<int32_t>> fp_mmr_const;
  std::vector<std::vector<std::vector<std::vector<int32_t>>>> fp_mmr_coef;

  uint16_t nlq_offset[3];
  uint32_t fp_hdr_in_max[3];
  uint32_t fp_linear_deadzone_slope[3];
  uint32_t fp_linear_deadzone_threshold[3];
};


constexpr uint16_t DoViProcessor::Clip3(uint16_t lower, uint16_t upper, int value)
{
  return max(min(value, upper), lower);
}

constexpr uint16_t DoViProcessor::upsampleHorzEven(const uint16_t* y, int n)
{
  return y[n];
}

constexpr uint16_t DoViProcessor::upsampleHorzOdd(const uint16_t* y, int n)
{
  auto val = (22 * y[n - 3] + 94 * y[n - 2] - 524 * y[n - 1] + 2456 * y[n] + 2456 * y[n + 1] - 524 * y[n + 2] +
    94 * y[n + 3] + 22 * y[n + 4] + 2048) >> 12;
  return Clip3(0, 0xFFFF, val);
}

constexpr uint16_t DoViProcessor::upsampleElYvertEven(const uint16_t* y, int n)
{
  auto val = (-3 * y[n - 2] + 29 * y[n - 1] + 111 * y[n] - 9 * y[n + 1] + 64) >> 7;
  return Clip3(0, 0xFFFF, val);
}

constexpr uint16_t DoViProcessor::upsampleElYvertOdd(const uint16_t* y, int n)
{
  auto val = (-9 * y[n - 1] + 111 * y[n] + 29 * y[n + 1] - 3 * y[n + 2] + 64) >> 7;
  return Clip3(0, 0xFFFF, val);
}

constexpr uint16_t DoViProcessor::upsampleElUVvertEven(const uint16_t* y, int n)
{
  auto val = (64 * y[n - 1] + 192 * y[n] + 128) >> 8;
  return min(0xFFFF, val);
}

constexpr uint16_t DoViProcessor::upsampleElUVvertOdd(const uint16_t* y, int n)
{
  auto val = (192 * y[n] + 64 * y[n + 1] + 128) >> 8;
  return min(0xFFFF, val);
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
