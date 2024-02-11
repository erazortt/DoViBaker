#pragma once
#include "avisynth.h"
#include "timecube/timecube.h"
#include <vector>
#include <string>

class DoViCubes: public GenericVideoFilter
{
public:
  DoViCubes(
    PClip child,
    std::vector<std::pair<uint16_t, std::string>>& cubes,
    bool fullrange,
    IScriptEnvironment* env);
  virtual ~DoViCubes();
  PVideoFrame GetFrame(int n, IScriptEnvironment* env) override;

private:
  struct TimecubeLutFree {
    void operator()(timecube_lut* ptr) { timecube_lut_free(ptr); }
  };

  struct TimecubeFilterFree {
    void operator()(timecube_filter* ptr) { timecube_filter_free(ptr); }
  };

  void applyLut(PVideoFrame& dst, const PVideoFrame& src) const;

  bool fullrange;
  std::vector<std::pair<uint16_t, timecube_filter*>> luts;
  const timecube_filter* currentFrameLut;
};
