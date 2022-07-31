#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

#include "DoViBaker.h"
#include <vector>
#include <string>
#include <sstream>

AVSValue __cdecl Create_RealDoViBaker(
  PClip blclip, 
  PClip elclip, 
  const char* rpuPath, 
  std::string cubeFiles, 
  std::string nits, 
  std::string cubesBasePath,
  bool qnd,
  bool rgbProof,
  bool nlqProof,
  const AVSValue* args, 
  IScriptEnvironment* env)
{
  if (!blclip->GetVideoInfo().HasVideo() || (elclip && !elclip->GetVideoInfo().HasVideo())) {
    env->ThrowError("DoViBaker: Clip not available");
  }

  if (!blclip->GetVideoInfo().IsYUV() || (elclip && !elclip->GetVideoInfo().IsYUV())) {
    env->ThrowError("DoViBaker: Clip must be in YUV format");
  }

  int blClipChromaSubSampled = -1;
  if (blclip->GetVideoInfo().Is420()) {
    blClipChromaSubSampled = 1;
  }
  if (blclip->GetVideoInfo().Is444()) {
    blClipChromaSubSampled = 0;
  }
  if (blClipChromaSubSampled < 0) {
    env->ThrowError("DoViBaker: Only 444 and 420 subsampling allowed");
  }

  int elClipChromaSubSampled = blClipChromaSubSampled;
  if (elclip) {
    elClipChromaSubSampled = -1;
    if (elclip && elclip->GetVideoInfo().Is420()) {
      elClipChromaSubSampled = 1;
    }
    if (elclip && elclip->GetVideoInfo().Is444()) {
      elClipChromaSubSampled = 0;
    }
    if (elclip && elClipChromaSubSampled < 0) {
      env->ThrowError("DoViBaker: Only 444 and 420 subsampling allowed");
    }
  }

  int quarterResolutionEl = 0;
  if (elclip) {
    quarterResolutionEl = -1;
    if ((blclip->GetVideoInfo().width == elclip->GetVideoInfo().width) && (blclip->GetVideoInfo().height == elclip->GetVideoInfo().height)) {
      quarterResolutionEl = 0;
    }
    if ((blclip->GetVideoInfo().width == 2 * elclip->GetVideoInfo().width) && (blclip->GetVideoInfo().height == 2 * elclip->GetVideoInfo().height)) {
      quarterResolutionEl = 1;
    }
    if (quarterResolutionEl < 0) {
      env->ThrowError("DoViBaker: Enhancement Layer must either be same size or quarter size as Base Layer");
    }
  }
  std::stringstream ssCubeFiles(cubeFiles);
  std::vector<std::string> cubesList;
  std::string segment;
  while (std::getline(ssCubeFiles, segment, ';'))
  {
    segment.insert(0, cubesBasePath);
    cubesList.push_back(segment);
  }
  std::stringstream ssNits(nits);
  std::vector<uint16_t> nitsList;
  while (std::getline(ssNits, segment, ';'))
  {    
    nitsList.push_back(std::atoi(segment.c_str()));
  }
  std::vector<std::pair<uint16_t, std::string>> cubeNitsPairs;
  if (cubesList.size() > 0) {
    if (cubesList.size() <= nitsList.size()) {
      env->ThrowError("DoViBaker: List of LUTs must be one entry longer then the list of nits.");
    }
    cubeNitsPairs.push_back(std::pair(0, cubesList[0]));
    for (int i = 0; i < nitsList.size(); i++) {
      cubeNitsPairs.push_back(std::pair(nitsList[i], cubesList[i + 1]));
    }
  }
  
  if (quarterResolutionEl == 0) {
    return new DoViBaker<false>(blclip, elclip, rpuPath, blClipChromaSubSampled, elClipChromaSubSampled, cubeNitsPairs, qnd, rgbProof, nlqProof, env);
  }
  if (quarterResolutionEl == 1) {
    return new DoViBaker<true>(blclip, elclip, rpuPath, blClipChromaSubSampled, elClipChromaSubSampled, cubeNitsPairs, qnd, rgbProof, nlqProof, env);
  }
}

AVSValue __cdecl Create_DoViBaker(AVSValue args, void* user_data, IScriptEnvironment* env)
{
  //args.ArraySize()
  return Create_RealDoViBaker(
    args[0].AsClip(), 
    args[1].AsClip(), 
    args[2].AsString(), 
    args[3].AsString(""), 
    args[4].AsString(""), 
    args[5].AsString(""),
    args[6].AsBool(false),
    args[7].AsBool(false),
    args[8].AsBool(false),
    &args, env);
}

const AVS_Linkage *AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
  AVS_linkage = vectors;

  env->AddFunction("DoViBaker", "c[el]c[rpu]s[cubes]s[mclls]s[cubes_basepath]s[qnd]b[rgbProof]b[nlqProof]b", Create_DoViBaker, 0);

  return "Hey it is just a spectrogram!";
}

void ypp2ycc(uint16_t* ycc, float y, float u, float v) {
  static const uint16_t scale = 1 << DoViProcessor::containerBitDepth;
  static const uint16_t bias = 16 << (DoViProcessor::containerBitDepth - 8);
  static const uint16_t ltop = scale - (21 << (DoViProcessor::containerBitDepth - 8));
  static const uint16_t ctop = scale - (16 << (DoViProcessor::containerBitDepth - 8));

  ycc[0] = y * (ltop - bias) + bias;
  ycc[1] = (u + 0.5) * (ctop - bias) + bias;
  ycc[2] = (v + 0.5) * (ctop - bias) + bias;
}

inline uint16_t normalizeSample(uint16_t &sample) {
  // we just check for 8bit precicion deviations, assuming smaller differences to be not visible
  return (sample >>= (DoViProcessor::containerBitDepth - 8));
}

void normalizeRgb(uint16_t* rgb) {
  normalizeSample(rgb[0]);
  normalizeSample(rgb[1]);
  normalizeSample(rgb[2]);
}

bool checkElProcessing(const DoViProcessor &dovi) {
  uint16_t yuv[3];
  ypp2ycc(yuv, 0.5000, 0.0000, 0.0000);
  uint16_t inGrey = yuv[0];
  ypp2ycc(yuv, 1.0000, 0.0000, 0.0000);
  uint16_t elHi = yuv[0];
  ypp2ycc(yuv, 0.0000, 0.0000, 0.0000);
  uint16_t elLo = yuv[0];
  uint16_t outHi = dovi.processSampleY(inGrey, elHi);
  uint16_t outLo = dovi.processSampleY(inGrey, elLo);
  return outHi != outLo;
}

bool checkMatrix(const DoViProcessor& dovi) {
  uint16_t yuv[3];
  uint16_t rgb[3];

  static const uint16_t maxNormRGB = 255;
  static const uint16_t halfNormRGB = maxNormRGB >> 1;

  ypp2ycc(yuv, 1.0000, 0.0000, 0.0000);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  if (rgb[0] != maxNormRGB || rgb[1] != maxNormRGB || rgb[2] != maxNormRGB) return true;
  
  ypp2ycc(yuv, 0.0000, 0.0000, 0.0000);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  if (rgb[0] != 0 || rgb[1] != 0 || rgb[2] != 0) return true;

  ypp2ycc(yuv, 0.5000, 0.0000, 0.0000);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  if (rgb[0] != halfNormRGB || rgb[1] != halfNormRGB || rgb[2] != halfNormRGB) return true;

  ypp2ycc(yuv, 0.2627 / 2, -0.1396 / 2, 0.5000 / 2);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  if (rgb[0] != halfNormRGB || rgb[1] != 0 || rgb[2] != 0) return true;

  ypp2ycc(yuv, 0.6780 / 2, -0.3604 / 2, -0.4598 / 2);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  if (rgb[0] != 0 || rgb[1] != halfNormRGB || rgb[2] != 0) return true;
  
  ypp2ycc(yuv, 0.0593 / 2, 0.5000 / 2, -0.0402 / 2);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  if (rgb[0] != 0 || rgb[1] != 0 || rgb[2] != halfNormRGB) return true;

  return false;
}

bool checkNonIdentityMapping(const DoViProcessor& dovi) {
  uint16_t yuv[3];
  uint16_t ely = dovi.getNlqOffset(0);
  for (int i = 0; i <= 10; i++) {
    ypp2ycc(yuv, float(i) / 10.0, 0.0000, 0.0000);
    uint16_t bly = yuv[0];
    uint16_t y = dovi.processSampleY(bly, ely);
    if (normalizeSample(y) != normalizeSample(bly)) {
      return true;
    }
  }
  return false;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    printf("DoViAnalyzer: needing path to RPU.bin file\n");
    return 1;
  }
  
  DoViProcessor dovi(argv[1], NULL);
  if (!dovi.wasCreationSuccessful()) {
    return 1;
  }

  /*
  printf("Self test\n");
  printf("2081 pq = %i\n", DoViProcessor::pq2nits(2081));
  for (int i = 0; i <= 100; i += 10) {
    std::string out(std::to_string(i));
    out += "%: ";
    out += std::to_string(DoViProcessor::pq2nits(4095 * i * 0.01));
    out += "\n";
    printf(out.c_str());
  }*/

  int length = dovi.getClipLength();
  printf("clip length: %i\n", length);
  int clip_max_pq = 0;
  bool elMixing = false;
  bool unusualMatrix = false;
  bool nonIdentityMapping = false;
  for (int i = 0; i < length; i++) {
    dovi.intializeFrame(i, NULL);
    int frame_max_pq = dovi.getMaxPq();
    if (frame_max_pq > clip_max_pq) {
      clip_max_pq = frame_max_pq;
    }
    unusualMatrix |= checkMatrix(dovi);
    nonIdentityMapping |= checkNonIdentityMapping(dovi);
    elMixing |= checkElProcessing(dovi);
  }
  //printf("clip max pq: %i\n", clip_max_pq);
  printf("overall max cll: %i\n", DoViProcessor::pq2nits(clip_max_pq));
  printf("unusual color matrix: %i\n", unusualMatrix);
  printf("mapping non-identity: %i\n", nonIdentityMapping);
  printf("el-clip processing: %i\n", elMixing);

  /*
  uint16_t yuv[3];
  uint16_t rgb[3];

  ypp2ycc(yuv, 0.5000, 0.0000, 0.0000);
  uint16_t inGrey = yuv[0];
  ypp2ycc(yuv, 1.0000, 0.0000, 0.0000);
  uint16_t elHi = yuv[0];
  ypp2ycc(yuv, 0.0000, 0.0000, 0.0000);
  uint16_t elLo = yuv[0];
  uint16_t outHi = dovi.processSampleY(inGrey, elHi);
  uint16_t outLo = dovi.processSampleY(inGrey, elLo);
  printf("transfere: %f\n",float(outHi)/float(outLo));

  ypp2ycc(yuv, 1.0000, 0.0000, 0.0000);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  printf("white = %i %i %i\n", rgb[0], rgb[1], rgb[2]);

  ypp2ycc(yuv, 0.0000, 0.0000, 0.0000);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  printf("black = %i %i %i\n", rgb[0], rgb[1], rgb[2]);

  ypp2ycc(yuv, 0.5000, 0.0000, 0.0000);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  printf("50%% grey = %i %i %i\n", rgb[0], rgb[1], rgb[2]);

  ypp2ycc(yuv, 0.2627/2, -0.1396/2, 0.5000/2);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  printf("50%% red = %i %i %i\n", rgb[0], rgb[1], rgb[2]);

  ypp2ycc(yuv, 0.6780/2, -0.3604/2, -0.4598/2);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  printf("50%% green = %i %i %i\n", rgb[0], rgb[1], rgb[2]);

  ypp2ycc(yuv, 0.0593/2, 0.5000/2, -0.0402/2);
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  normalizeRgb(rgb);
  printf("50%% blue = %i %i %i\n", rgb[0], rgb[1], rgb[2]);
  */
}
