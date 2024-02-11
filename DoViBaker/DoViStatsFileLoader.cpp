#include "DoViStatsFileLoader.h"
#include "DoViProcessor.h"
#include <fstream>
#include <sstream>
#include <deque>
#include <algorithm>

DoViStatsFileLoader::DoViStatsFileLoader(
	PClip child,
	std::string maxPqFile,
	std::string sceneCutFile,
	IScriptEnvironment* env)
	: GenericVideoFilter(child)
	, currentScene(0)
	, previousFrame(0)
	, staticMaxPq(0)
	, staticMaxCll(0)
{
	uint32_t frame = 0, isLastFrameInScene, frameMaxPq, frameMinPq, firstFrameNextScene = 0;
	uint16_t sceneMaxPq = 0, sceneMinPq = -1;
	float scale;
	std::deque<float> sceneScales;
	std::ifstream fpStats, fpSceneCut;
	std::string line, segment;

	fpStats.open(maxPqFile, std::ifstream::in);
	if (!fpStats.is_open()) {
		env->ThrowError((std::string("DoViMaxMqFileReader: cannot find stats file ") + maxPqFile).c_str());
	}
	if (!sceneCutFile.empty()) {
		fpSceneCut.open(sceneCutFile, std::ifstream::in);
		if (!fpSceneCut.is_open()) {
			env->ThrowError((std::string("DoViMaxMqFileReader: cannot find scene cut file ") + sceneCutFile).c_str());
		}
		while (firstFrameNextScene == 0) {
			if (!(fpSceneCut >> firstFrameNextScene)) {
				env->ThrowError((std::string("DoViMaxMqFileReader: error reading scene cut file ") + sceneCutFile).c_str());
			}
		}
	} 

	while (std::getline(fpStats, line))
	{
		std::istringstream ssline(line);
		if (std::getline(ssline, segment, ' ')) {
			frame = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading frame number from stats file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			isLastFrameInScene = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading scene change from stats file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			frameMaxPq = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading maxPq from stats file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			frameMinPq = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading minPq from stats file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			scale = std::atof(segment.c_str());
			sceneScales.push_back(scale);
		}

		if (frameMaxPq > sceneMaxPq) sceneMaxPq = frameMaxPq;
		if (sceneMaxPq > staticMaxPq) staticMaxPq = sceneMaxPq;
		if (frameMinPq < sceneMinPq) sceneMinPq = frameMinPq;
		if (fpSceneCut.is_open()) {
			if (firstFrameNextScene != frame + 1) continue;
		} else if (!isLastFrameInScene) continue;

		float sceneScaleMedian = 1;
		if (!sceneScales.empty()) {
			std::sort(sceneScales.begin(), sceneScales.end());
			sceneScaleMedian = sceneScales.at(sceneScales.size() / 2);
			sceneScales.clear();
		}

		sceneMaxSignal.push_back(std::tuple(frame + 1, sceneMaxPq, sceneMinPq, sceneScaleMedian));
		sceneMaxPq = 0;
		sceneMinPq = -1;
		if (fpSceneCut.is_open()) {
			if (!(fpSceneCut >> firstFrameNextScene)) {
				firstFrameNextScene = child->GetVideoInfo().num_frames;
			}
		}
	}
	float sceneScaleMedian = 1;
	if (!sceneScales.empty()) {
		std::sort(sceneScales.begin(), sceneScales.end());
		sceneScaleMedian = sceneScales.at(sceneScales.size() / 2);
	}

	sceneMaxSignal.push_back(std::tuple(frame + 1, sceneMaxPq, sceneMinPq, sceneScaleMedian));
	staticMaxCll = DoViProcessor::pq2nits(staticMaxPq) + 0.5f;

	if (child->GetVideoInfo().num_frames != frame + 1) {
		env->ThrowError((std::string("DoViMaxMqFileReader: clip length does not match stats file ") + maxPqFile).c_str());
	}

	fpStats.close();
	if(fpSceneCut.is_open()) fpSceneCut.close();
}

PVideoFrame DoViStatsFileLoader::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);

	if (previousFrame > n)
		currentScene = 0;
	previousFrame = n;

	while (std::get<0>(sceneMaxSignal.at(currentScene)) <= n) {
		currentScene++;
	}

	uint16_t maxPq = std::get<1>(sceneMaxSignal.at(currentScene));
	uint16_t maxCll = DoViProcessor::pq2nits(maxPq) + 0.5f;
	uint16_t minPq = std::get<2>(sceneMaxSignal.at(currentScene));
	float scale = std::get<3>(sceneMaxSignal.at(currentScene));
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_min_pq", minPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_max_pq", maxPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_max_content_light_level", maxCll, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_static_max_pq", staticMaxPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_static_max_content_light_level", staticMaxCll, 0);
	env->propSetFloat(env->getFramePropsRW(src), "_dovi_dynamic_luminosity_scale", scale, 0);

	env->propSetInt(env->getFramePropsRW(src), "_SceneChangeNext", std::get<0>(sceneMaxSignal.at(currentScene)) == n + 1, 0);
	if (currentScene > 0) {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", std::get<0>(sceneMaxSignal.at(currentScene - 1)) == n, 0);
	} else {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", 0, 0);
	}

	return src;
}