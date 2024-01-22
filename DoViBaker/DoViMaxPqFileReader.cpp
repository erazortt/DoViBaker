#include "DoViMaxPqFileReader.h"
#include "DoViTonemap.h"
#include <fstream>
#include <sstream>
#include <deque>
#include <algorithm>

DoViMaxPqFileReader::DoViMaxPqFileReader(
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
	uint32_t frame = 0, isLastFrameInScene, pq, firstFrameNextScene = 0;
	uint16_t maxPq = 0;
	uint8_t scale;
	std::deque<uint8_t> sceneScales;
	std::ifstream fpMaxPq, fpSceneCut;
	std::string line, segment;

	fpMaxPq.open(maxPqFile);
	if (!fpMaxPq.is_open()) {
		env->ThrowError((std::string("DoViMaxMqFileReader: cannot find maxPq file ") + maxPqFile).c_str());
	}
	if (!sceneCutFile.empty()) {
		fpSceneCut.open(sceneCutFile);
		if (!fpSceneCut.is_open()) {
			env->ThrowError((std::string("DoViMaxMqFileReader: cannot find scene cut file ") + sceneCutFile).c_str());
		}
		if (!(fpSceneCut >> firstFrameNextScene)) {
			env->ThrowError((std::string("DoViMaxMqFileReader: error reading scene cut file ") + sceneCutFile).c_str());
		}
	} 

	while (std::getline(fpMaxPq, line))
	{
		std::istringstream ssline(line);
		if (std::getline(ssline, segment, ' ')) {
			frame = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading maxPq file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			isLastFrameInScene = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading maxPq file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			pq = std::atoi(segment.c_str());
		} else env->ThrowError((std::string("DoViMaxMqFileReader: error reading maxPq file ") + sceneCutFile).c_str());
		if (std::getline(ssline, segment, ' ')) {
			scale = std::atoi(segment.c_str());
			sceneScales.push_back(scale);
		}

		if (pq > maxPq) maxPq = pq;
		if (maxPq > staticMaxPq) staticMaxPq = maxPq;
		if (fpSceneCut.is_open()) {
			if (firstFrameNextScene != frame + 1) continue;
		} else if (!isLastFrameInScene) continue;

		uint8_t sceneScaleMedian = -1;
		if (!sceneScales.empty()) {
			std::sort(sceneScales.begin(), sceneScales.end());
			sceneScaleMedian = sceneScales.at(sceneScales.size() / 2);
		}

		uint16_t maxCll = DoViTonemap::pq2nits(maxPq) + 0.5;
		sceneMaxSignal.push_back(std::tuple(frame + 1, maxPq, maxCll, sceneScaleMedian));
		maxPq = 0;
		if (fpSceneCut.is_open()) {
			if (!(fpSceneCut >> firstFrameNextScene)) {
				firstFrameNextScene = child->GetVideoInfo().num_frames;
			}
		}
	}
	uint8_t sceneScaleMedian = -1;
	if (!sceneScales.empty()) {
		std::sort(sceneScales.begin(), sceneScales.end());
		sceneScaleMedian = sceneScales.at(sceneScales.size() / 2);
	}

	uint16_t maxCll = DoViTonemap::pq2nits(maxPq) + 0.5;
	sceneMaxSignal.push_back(std::tuple(frame + 1, maxPq, maxCll, sceneScaleMedian));
	staticMaxCll = DoViTonemap::pq2nits(staticMaxPq) + 0.5;

	if (child->GetVideoInfo().num_frames != frame + 1) {
		env->ThrowError((std::string("DoViMaxMqFileReader: clip length does not match maxPq file ") + maxPqFile).c_str());
	}

	fpMaxPq.close();
	if(fpSceneCut.is_open()) fpSceneCut.close();
}

PVideoFrame DoViMaxPqFileReader::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);

	if (previousFrame > n)
		currentScene = 0;
	previousFrame = n;

	while (std::get<0>(sceneMaxSignal.at(currentScene)) <= n) {
		currentScene++;
	}

	uint16_t maxPq = std::get<1>(sceneMaxSignal.at(currentScene));
	uint16_t maxCll = std::get<2>(sceneMaxSignal.at(currentScene));
	float scale = std::get<3>(sceneMaxSignal.at(currentScene));
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_max_pq", maxPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_max_content_light_level", maxCll, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_static_max_pq", staticMaxPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_static_max_content_light_level", staticMaxCll, 0);
	if (scale > -1) {
		env->propSetFloat(env->getFramePropsRW(src), "_dovi_dynamic_luminosity_scale", scale/10, 0);
	}
	env->propSetInt(env->getFramePropsRW(src), "_SceneChangeNext", std::get<0>(sceneMaxSignal.at(currentScene)) == n + 1, 0);
	if (currentScene > 0) {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", std::get<0>(sceneMaxSignal.at(currentScene - 1)) == n, 0);
	} else {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", 0, 0);
	}

	return src;
}