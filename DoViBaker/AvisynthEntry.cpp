#include <sstream>

#include "def.h"
#include "DoViBaker.h"
#include "DoViCubes.h"
#include "DoViTonemap.h"
#include "DoViStatsFileLoader.h"

AVSValue __cdecl Create_RealDoViBaker(
  PClip blclip,
  PClip elclip,
  const char* rpuPath,
  uint16_t desiredTrimPq,
  float targetMaxNits,
  float targetMinNits,
  bool qnd,
  bool rgbProof,
  bool nlqProof,
  int sourceProfile,
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

  if (sourceProfile != 7 && sourceProfile != 8) {
      env->ThrowError("DoViBaker: sourceProfile must be either 7 or 8.");
  }

  if (quarterResolutionEl == 0) {
    return new DoViBaker<false>(blclip, elclip, rpuPath, blClipChromaSubSampled, elClipChromaSubSampled, desiredTrimPq, targetMinNits, targetMaxNits, qnd, rgbProof, nlqProof, sourceProfile, env);
  }
  if (quarterResolutionEl == 1) {
    return new DoViBaker<true>(blclip, elclip, rpuPath, blClipChromaSubSampled, elClipChromaSubSampled, desiredTrimPq, targetMinNits, targetMaxNits, qnd, rgbProof, nlqProof, sourceProfile, env);
  }
}

AVSValue __cdecl Create_DoViBaker(AVSValue args, void* user_data, IScriptEnvironment* env)
{
  auto elClip = (args[1].Defined() ? args[1].AsClip() : nullptr);

  return Create_RealDoViBaker(
    args[0].AsClip(),
    elClip,
    args[2].AsString(""),
    args[3].AsInt(0),
    args[4].AsFloat(100),
    args[5].AsFloat(0),
    args[6].AsBool(false),
    args[7].AsBool(false),
    args[8].AsBool(false),
    args[9].AsInt(7),
    env);
}

AVSValue __cdecl Create_RealDoViCubes(
  PClip clip,
  std::string cubeFiles,
  std::string nits,
  std::string cubesBasePath,
  bool fullrange,
  const AVSValue* args,
  IScriptEnvironment* env)
{
  if (!clip->GetVideoInfo().IsPlanarRGB())
  {
    env->ThrowError("DoViTonemap: input must be planar RGB");
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
  return new DoViCubes(clip, cubeNitsPairs, fullrange, env);
}

AVSValue __cdecl Create_DoViCubes(AVSValue args, void* user_data, IScriptEnvironment* env)
{
  return Create_RealDoViCubes(
    args[0].AsClip(),
    args[1].AsString(""),
    args[2].AsString(""),
    args[3].AsString(""),
    args[4].AsBool(true),
    &args, env);
}

AVSValue __cdecl Create_RealDoViStatsFileLoader(
  PClip clip,
  std::string maxPqFile,
  std::string sceneCutFile,
  const AVSValue* args,
  IScriptEnvironment* env)
{
  return new DoViStatsFileLoader(clip, maxPqFile, sceneCutFile, env);
}

AVSValue __cdecl Create_DoViStatsFileLoader(AVSValue args, void* user_data, IScriptEnvironment* env)
{
  return Create_RealDoViStatsFileLoader(
    args[0].AsClip(),
    args[1].AsString(""),
    args[2].AsString(""),
    &args, env);
}

AVSValue __cdecl Create_RealDoViTonemap(
  PClip clip,
  float targetMaxNits,
  float targetMinNits,
  float masterMaxNits,
  float masterMinNits,
  float lumScale,
  float kneeOffset,
  bool normalizeOutput,
  const AVSValue* args,
  IScriptEnvironment* env)
{
  if (!clip->GetVideoInfo().IsPlanarRGB())
  {
    env->ThrowError("DoViTonemap: input must be planar RGB");
  }
  if (kneeOffset < 0.5 || kneeOffset > 2){
    env->ThrowError("DoViTonemap: valid range of knee offset is [0.5; 2.0]");
  }
  if (targetMaxNits < 0 || targetMinNits < 0) {
    env->ThrowError("DoViTonemap: target capabilities must be given explicitly");
  }

  switch (clip->GetVideoInfo().BitsPerComponent())
  {
  case 10: return new DoViTonemap<10>(clip, targetMaxNits, targetMinNits, masterMaxNits, masterMinNits, lumScale, kneeOffset, normalizeOutput, env); break;
  case 12: return new DoViTonemap<12>(clip, targetMaxNits, targetMinNits, masterMaxNits, masterMinNits, lumScale, kneeOffset, normalizeOutput, env); break;
  case 14: return new DoViTonemap<14>(clip, targetMaxNits, targetMinNits, masterMaxNits, masterMinNits, lumScale, kneeOffset, normalizeOutput, env); break;
  case 16: return new DoViTonemap<16>(clip, targetMaxNits, targetMinNits, masterMaxNits, masterMinNits, lumScale, kneeOffset, normalizeOutput, env); break;
  default:
    env->ThrowError("DoViTonemap: input bit depth not compatible");
    break;
  }
  return 0x0;
}


AVSValue __cdecl Create_DoViTonemap(AVSValue args, void* user_data, IScriptEnvironment* env)
{
  return Create_RealDoViTonemap(
    args[0].AsClip(),
    args[1].AsFloat(1000),
    args[2].AsFloat(0),
    args[3].AsFloat(-1),
    args[4].AsFloat(-1),
    args[5].AsFloat(1),
    args[6].AsFloat(0.75),
    args[7].AsBool(false),
    &args, env);
}

#ifdef DOVI_BAKER
const AVS_Linkage *AVS_linkage = nullptr;
extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
  AVS_linkage = vectors;

  env->AddFunction("DoViBaker", "c[el]c[rpu]s[trimPq]i[targetMaxNits]f[targetMinNits]f[qnd]b[rgbProof]b[nlqProof]b[sourceProfile]i", Create_DoViBaker, 0);
  env->AddFunction("DoViTonemap", "c[targetMaxNits]f[targetMinNits]f[masterMaxNits]f[masterMinNits]f[lumScale]f[kneeOffset]f[normalizeOutput]b", Create_DoViTonemap, 0);
  env->AddFunction("DoViCubes", "c[cubes]s[mclls]s[cubes_basepath]s[fullrange]b", Create_DoViCubes, 0);
  env->AddFunction("DoViStatsFileLoader", "c[statsFile]s[sceneCutsFile]s", Create_DoViStatsFileLoader, 0);

  return "Hey it is just a spectrogram!";
}
#endif //DOVI_BAKER

void ypp2ycc(uint16_t* ycc, float y, float u, float v) {
  //YPrPb to YCrCb
  static const uint16_t scale = 1 << DoViProcessor::outContainerBitDepth;
  static const uint16_t bias = 16 << (DoViProcessor::outContainerBitDepth - 8);
  static const uint16_t ltop = scale - (21 << (DoViProcessor::outContainerBitDepth - 8));
  static const uint16_t ctop = scale - (16 << (DoViProcessor::outContainerBitDepth - 8));

  ycc[0] = y * (ltop - bias) + bias;
  ycc[1] = (u + 0.5) * (ctop - bias) + bias;
  ycc[2] = (v + 0.5) * (ctop - bias) + bias;
}

inline uint16_t to8bits(uint16_t sample) {
  // we just check for 8bit precicion deviations, assuming smaller differences to be not visible
  return (sample >> (DoViProcessor::outContainerBitDepth - 8));
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

uint16_t checkMatrix(const DoViProcessor& dovi) {
  uint16_t yuv[3];
  uint16_t rgb[3];
  uint16_t diffBits = 0;

  static const uint16_t maxNormRGB = 255 << (DoViProcessor::outContainerBitDepth - 8);
  static const uint16_t halfNormRGB = maxNormRGB >> 1;

  ypp2ycc(yuv, 1.0000, 0.0000, 0.0000); // 10000nits white
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  diffBits |= std::abs(rgb[0] - maxNormRGB);
  diffBits |= std::abs(rgb[1] - maxNormRGB);
  diffBits |= std::abs(rgb[2] - maxNormRGB);

  ypp2ycc(yuv, 0.0000, 0.0000, 0.0000); // black
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  diffBits |= rgb[0];
  diffBits |= rgb[1];
  diffBits |= rgb[2];

  ypp2ycc(yuv, 0.5000, 0.0000, 0.0000); // 92nits white
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  diffBits |= std::abs(rgb[0] - halfNormRGB);
  diffBits |= std::abs(rgb[1] - halfNormRGB);
  diffBits |= std::abs(rgb[2] - halfNormRGB);

  ypp2ycc(yuv, 0.2627 / 2, -0.1396 / 2, 0.5000 / 2); // half intensity red
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  diffBits |= std::abs(rgb[0] - halfNormRGB);
  diffBits |= rgb[1];
  diffBits |= rgb[2];

  ypp2ycc(yuv, 0.6780 / 2, -0.3604 / 2, -0.4598 / 2); // half intensity blue
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  diffBits |= rgb[0];
  diffBits |= std::abs(rgb[1] - halfNormRGB);
  diffBits |= rgb[2];

  ypp2ycc(yuv, 0.0593 / 2, 0.5000 / 2, -0.0402 / 2); // half intensity green
  dovi.sample2rgb(rgb[0], rgb[1], rgb[2], yuv[0], yuv[1], yuv[2]);
  diffBits |= rgb[0];
  diffBits |= rgb[1];
  diffBits |= std::abs(rgb[2] - halfNormRGB);

  return diffBits;
}

uint16_t checkNonIdentityMapping(const DoViProcessor& dovi) {
  uint16_t yuv[3];
  uint16_t ely = dovi.getNlqOffset(0);
  uint16_t diffBits = 0;
  for (int i = 0; i <= 10; i++) {
    ypp2ycc(yuv, float(i) / 10.0, 0.0000, 0.0000);
    uint16_t bly = yuv[0];
    int y = dovi.processSampleY(bly, ely);
    diffBits |= std::abs(y - bly);
  }
  return diffBits;
}

double Spline16Filter(double value) {
  value = fabs(value);

  if (value < 1.0) {
    return ((value - 9.0 / 5.0) * value - 1.0 / 5.0) * value + 1.0;
  }
  else if (value < 2.0) {
    return ((-1.0 / 3.0 * (value - 1.0) + 4.0 / 5.0) * (value - 1.0) - 7.0 / 15.0) * (value - 1.0);
  }
  return 0.0;
}

double Spline36Filter(double value) {
  value = fabs(value);

  if (value < 1.0) {
    return ((13.0 / 11.0 * (value)-453.0 / 209.0) * (value)-3.0 / 209.0) * (value)+1.0;
  }
  else if (value < 2.0) {
    return ((-6.0 / 11.0 * (value - 1.0) + 270.0 / 209.0) * (value - 1.0) - 156.0 / 209.0) * (value - 1.0);
  }
  else if (value < 3.0) {
    return  ((1.0 / 11.0 * (value - 2.0) - 45.0 / 209.0) * (value - 2.0) + 26.0 / 209.0) * (value - 2.0);
  }
  return 0.0;
}

double Spline64Filter(double value) {
  value = fabs(value);

  if (value < 1.0) {
    return ((49.0 / 41.0 * (value)-6387.0 / 2911.0) * (value)-3.0 / 2911.0) * (value)+1.0;
  }
  else if (value < 2.0) {
    return ((-24.0 / 41.0 * (value - 1.0) + 4032.0 / 2911.0) * (value - 1.0) - 2328.0 / 2911.0) * (value - 1.0);
  }
  else if (value < 3.0) {
    return ((6.0 / 41.0 * (value - 2.0) - 1008.0 / 2911.0) * (value - 2.0) + 582.0 / 2911.0) * (value - 2.0);
  }
  else if (value < 4.0) {
    return ((-1.0 / 41.0 * (value - 3.0) + 168.0 / 2911.0) * (value - 3.0) - 97.0 / 2911.0) * (value - 3.0);
  }
  return 0.0;
}

//#define DOVI_ANALYZER
#ifdef DOVI_ANALYZER
int main(int argc, char* argv[])
{
  bool selfTest = false;
  bool showSplines = false;
  bool showNitsTable = false;
  bool showTonemap = false;
  if (argc < 2) {
    printf("DoViAnalyzer: provide path to RPU.bin file\n");
    return 1;
  }
  if (argv[1][0] == '-') {
    selfTest = true;
    printf("Self test\n");
    if (strcmp(argv[1], "-showSplines") == 0)
      showSplines = true;
    else if (strcmp(argv[1], "-showNitsTable") == 0)
      showNitsTable = true;
    else if (strcmp(argv[1], "-showTonemap") == 0)
      showTonemap = true;
    else {
      printf("DoViAnalyzer: self test not regognized\n");
      return 1;
    }
  }
  if (showSplines) {
    printf("Spline64 coefficients\n");
    for (int i = 0; i < 4; i++) {
      float val = i + 0.5;
      printf("%f %i\n", val, int(Spline64Filter(val) * (1 << 12)));
    }
    printf("Spline16 coefficients\n");
    for (int i = 0; i < 2; i++) {
      float val = i + 0.5;
      printf("%f %i\n", val, int(Spline16Filter(val) * (1 << 12)));
    }
    printf("Spline36 coefficients\n");
    for (int i = -2; i < 4; i++) {
      float val = i - 0.25;
      printf("%f %i\n", val, int(Spline36Filter(val) * (1 << 8)));
    }
    printf("Spline16 coefficients\n");
    for (int i = -1; i < 3; i++) {
      float val = i - 0.25;
      printf("%f %i\n", val, int(Spline16Filter(val) * (1 << 7)));
    }
  }

  if (showNitsTable) {
    printf("100 nits = %i pq\n", DoViProcessor::nits2pq(100));
    printf("600 nits = %i pq\n", DoViProcessor::nits2pq(600));
    printf("1000 nits = %i pq\n", DoViProcessor::nits2pq(1000));
    printf("2081 pq = %f nits\n", DoViProcessor::pq2nits(2081));
    printf("2851 pq = %f nits\n", DoViProcessor::pq2nits(2851));
    printf("3079 pq = %f nits\n", DoViProcessor::pq2nits(3079));
    for (int i = 0; i <= 120; i += 10) {
      printf("%i%%: %f\n", i, DoViProcessor::pq2nits(4095 * i * 0.01));
    }
  }

  if (showTonemap) {
    uint16_t targetMaxPq = DoViProcessor::nits2pq(1000);
    uint16_t targetMinPq = DoViProcessor::nits2pq(2);
    uint16_t masterMaxPq = DoViProcessor::nits2pq(2000);
    uint16_t masterMinPq = DoViProcessor::nits2pq(1);
    const int bitDepth = 12;
    DoViEetf<bitDepth> tonemap(0.75, false);
    tonemap.generateEETF(targetMaxPq, targetMinPq, masterMaxPq, masterMinPq, 1.0, false);
    for (int i = 0; i <= 255; i++) {
      uint16_t signal = i * ((1 << bitDepth) - 1) / 255.0 + 0.5;
      uint16_t signalPq = i * 4095 / 255.0 + 0.5;
      uint16_t mapped = tonemap.applyEETF(signal);
      uint16_t mappedPq = mapped * 4095.0 / ((1 << bitDepth) - 1) + 0.5;
      printf("%i %f %f %i\n",signal, DoViProcessor::pq2nits(signalPq), DoViProcessor::pq2nits(mappedPq), mapped);
    }
  }

  if (selfTest) { return 1; }

  FILE* fpSceneChange = NULL;
  if (argc > 2) {
    fpSceneChange = fopen(argv[2], "w");
  }

  DoViProcessor dovi(argv[1], NULL, DoViProcessor::outContainerBitDepth, DoViProcessor::outContainerBitDepth, 0);
  if (!dovi.wasCreationSuccessful()) {
    return 1;
  }

  int length = dovi.getClipLength();
  printf("clip length: %i\n", length);
  dovi.setTrim(1, 5, 5); //enable trimming so that we get all possible trims
  int clip_max_pq = 0;
  bool elMixing = false;
  uint16_t unusualMatrix = 0;
  uint16_t nonIdentityMapping = 0;
  bool limitedRangeOutput = false;
  std::vector<uint16_t> trimPq;
  for (int i = 0; i < length; i++) {
    if (!dovi.intializeFrame(i, NULL, 0, 0)) {
      return 1;
    }
    int frame_max_pq = dovi.getDynamicMaxPq();
    if (frame_max_pq > clip_max_pq) {
      clip_max_pq = frame_max_pq;
    }
    unusualMatrix |= checkMatrix(dovi);
    nonIdentityMapping |= checkNonIdentityMapping(dovi);
    elMixing |= checkElProcessing(dovi);
    if(fpSceneChange && dovi.isSceneChange()){
      fputs((std::to_string(i)+"\n").c_str(), fpSceneChange);
    }
    if (dovi.isLimitedRangeOutput()) {
      limitedRangeOutput = true;
    }
    for (int i = 0; i < dovi.getAvailableTrimPqs().size(); i++) {
      if (std::find(trimPq.begin(), trimPq.end(), dovi.getAvailableTrimPqs()[i]) == trimPq.end()) {
        trimPq.push_back(dovi.getAvailableTrimPqs()[i]);
      }
    }
  }

  if (fpSceneChange) {
    fclose(fpSceneChange);
  }

  //printf("clip max pq: %i\n", clip_max_pq);
  printf("overall max cll: %i\n", int(DoViProcessor::pq2nits(clip_max_pq) + 0.5));
  printf("color matrix deviation: %i\n", to8bits(unusualMatrix));
  printf("mapping deviation: %i\n", to8bits(nonIdentityMapping));
  if(elMixing) printf("el-clip processing: enabled\n");
  else printf("el-clip processing: disabled\n");
  printf("available trims: ");
  for (int i = 0; i < trimPq.size(); i++) {
    printf("%i nits (trimPq: %i)\n", int(DoViProcessor::pq2nits(trimPq[i]) + 0.5), trimPq[i]);
    printf("                 ");
  }
  (trimPq.size() > 0) ? printf("\n") : printf("none\n");
  if (limitedRangeOutput) printf("output is limited range!\n");
  else printf("output is full range\n");
}
#endif // DOVI_ANALYZER