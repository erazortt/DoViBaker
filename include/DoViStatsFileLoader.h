#pragma once
#include "avisynth.h"
#include <vector>
#include <string>
#include <tuple>

class DoViStatsFileLoader : public GenericVideoFilter
{
public:
  DoViStatsFileLoader(
    PClip child,
    std::string maxPqFile,
    std::string sceneCutFile,
    IScriptEnvironment* env);
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  std::vector<std::tuple<uint32_t, uint16_t, uint16_t, float>> sceneMaxSignal;
  uint32_t currentScene;
  uint32_t previousFrame;
  uint16_t staticMaxPq;
  uint16_t staticMaxCll;
};

