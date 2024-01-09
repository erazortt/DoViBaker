#pragma once
#include "avisynth.h"
#include <vector>
#include <string>

class DoViSceneFileReader : public GenericVideoFilter
{
public:
  DoViSceneFileReader(
    PClip child,
    std::string sceneFile,
    IScriptEnvironment* env);
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  std::vector<std::pair<uint32_t, uint16_t>> sceneMaxSignal;
  uint32_t currentScene;
  uint32_t previousFrame;
};

