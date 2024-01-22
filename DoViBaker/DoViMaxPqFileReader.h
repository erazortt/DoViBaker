#pragma once
#include "avisynth.h"
#include <vector>
#include <string>

class DoViMaxPqFileReader : public GenericVideoFilter
{
public:
  DoViMaxPqFileReader(
    PClip child,
    std::string maxPqFile,
    std::string sceneCutFile,
    IScriptEnvironment* env);
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  std::vector<std::tuple<uint32_t, uint16_t, uint16_t, uint8_t>> sceneMaxSignal;
  uint32_t currentScene;
  uint32_t previousFrame;
  uint16_t staticMaxPq;
  uint16_t staticMaxCll;
};

