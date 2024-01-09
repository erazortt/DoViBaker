#include "DoViMaxPqFileReader.h"
#include "DoViProcessor.h"

DoViMaxPqFileReader::DoViMaxPqFileReader(
	PClip child,
	std::string sceneCutFile,
	std::string maxPqFile,
	IScriptEnvironment* env)
	: GenericVideoFilter(child)
	, currentScene(0)
	, previousFrame(0)
{
	uint32_t frame, firstFrameNextScene, pq;
	uint16_t maxPq = 0;
	uint16_t staticMaxPq = 0;
	FILE *fpSceneCut, *fpMaxPq;

	fpSceneCut = fopen(sceneCutFile.c_str(), "r");
	fpMaxPq = fopen(maxPqFile.c_str(), "r");
	if (!fpSceneCut) {
		env->ThrowError((std::string("DoViMaxMqFileReader: cannot find scene cut file ") + sceneCutFile).c_str());
	}
	if (!fpMaxPq) {
		env->ThrowError((std::string("DoViMaxMqFileReader: cannot find maxPq file ") + maxPqFile).c_str());
	}
	if (fscanf(fpSceneCut, "%i\n", &firstFrameNextScene) != 1) {
		env->ThrowError((std::string("DoViMaxMqFileReader: error reading scene cut file ") + sceneCutFile).c_str());
	}
	while (fscanf(fpMaxPq, "%i %i\n", &frame, &pq) == 2) {
		if (pq > maxPq) maxPq = pq;
		if (maxPq > staticMaxPq) staticMaxPq = maxPq;
		if (firstFrameNextScene != frame + 1) continue;
		uint16_t maxCll = DoViProcessor::pq2nits(maxPq);
		sceneMaxSignal.push_back(std::tuple(frame + 1, maxPq, maxCll));
		maxPq = 0;
		if (fscanf(fpSceneCut, "%i\n", &firstFrameNextScene) != 1) {
			firstFrameNextScene = child->GetVideoInfo().num_frames;
		}
	}
	uint16_t maxCll = DoViProcessor::pq2nits(maxPq);
	sceneMaxSignal.push_back(std::tuple(frame + 1, maxPq, maxCll));
	staticMaxCll = DoViProcessor::pq2nits(staticMaxPq);

	if (child->GetVideoInfo().num_frames != frame + 1) {
		env->ThrowError((std::string("DoViMaxMqFileReader: clip length does not match maxPq file ") + maxPqFile).c_str());
	}

	fclose(fpMaxPq);
	fclose(fpSceneCut);
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
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_max_pq", maxPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_dynamic_max_content_light_level", maxCll, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_static_max_content_light_level", staticMaxCll, 0);
	env->propSetInt(env->getFramePropsRW(src), "_SceneChangeNext", std::get<0>(sceneMaxSignal.at(currentScene)) == n + 1, 0);
	if (currentScene > 0) {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", std::get<0>(sceneMaxSignal.at(currentScene - 1)) == n, 0);
	} else {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", 0, 0);
	}

	return src;
}