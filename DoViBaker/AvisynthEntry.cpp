#include <windows.h>

#pragma warning(push)
#pragma warning(disable: 4512 4244 4100 693)
#include "avisynth.h"
#pragma warning(pop)

#include "DoViBaker.h"
#include <vector>

AVSValue __cdecl Create_RealDoViBaker(PClip blclip, PClip elclip, const char* rpuPath, const AVSValue* args, IScriptEnvironment* env)
{
  if (!blclip->GetVideoInfo().HasVideo() || !elclip->GetVideoInfo().HasVideo()) {
    env->ThrowError("DoViBaker: Clip not available");
  }

  if (!blclip->GetVideoInfo().IsYUV() || !elclip->GetVideoInfo().IsYUV()) {
    env->ThrowError("DoViBaker: Clip must be in YUV format");
  }

  int chromaSubSampling = -1;
  if (blclip->GetVideoInfo().Is420() && elclip->GetVideoInfo().Is420()) {
    chromaSubSampling = 1;
  }
  if (blclip->GetVideoInfo().Is444() && elclip->GetVideoInfo().Is444()) {
    chromaSubSampling = 0;
  }
  if (chromaSubSampling<0) {
    env->ThrowError("DoViBaker: Base Layer and Enhancement Layer clips must have both the same chroma subsampling, either 444 or 420");
  }

  int quarterResolutionEl = -1;
  if ((blclip->GetVideoInfo().width == elclip->GetVideoInfo().width) && (blclip->GetVideoInfo().height == elclip->GetVideoInfo().height)) {
    quarterResolutionEl = 0;
  }
  if ((blclip->GetVideoInfo().width == 2*elclip->GetVideoInfo().width) && (blclip->GetVideoInfo().height == 2*elclip->GetVideoInfo().height)) {
    quarterResolutionEl = 1;
  }
  if (quarterResolutionEl < 0) {
    env->ThrowError("DoViBaker: Enhancement Layer must either be same size or quarter size as Base Layer");
  }

  if (chromaSubSampling == 0) {
    return new DoViBaker<false>(blclip, elclip, rpuPath, quarterResolutionEl, env);
  }
  if (chromaSubSampling == 1) {
    return new DoViBaker<true>(blclip, elclip, rpuPath, quarterResolutionEl, env);
  }
}

AVSValue __cdecl Create_DoViBaker(AVSValue args, void* user_data, IScriptEnvironment* env)
{

  return Create_RealDoViBaker(args[0].AsClip(), args[1].AsClip(), args[2].AsString(), &args, env);
}

const AVS_Linkage *AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
  AVS_linkage = vectors;

  env->AddFunction("DoViBaker", "c[el]c[rpu]s", Create_DoViBaker, 0);

  return "Hey it is just a spectrogram!";
}

int main()
{
	DoViProcessor dovi("Z:/rpu.bin", NULL);
  dovi.intializeFrame(1, NULL);

}
