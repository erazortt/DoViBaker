#include <windows.h>

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

int main()
{
	DoViProcessor dovi("Z:/rpu.bin", NULL);
  dovi.intializeFrame(1, NULL);

  printf((std::string("2081 pq = ") + std::to_string(DoViProcessor::pq2nits(2081))).c_str());
  for (int i = 0; i < 100; i += 10) {
    std::string out(std::to_string(i));
    out += "%: ";
    out += std::to_string(DoViProcessor::pq2nits(4095 * i * 0.01));
    out += "\n";
    printf(out.c_str());
  }
}
