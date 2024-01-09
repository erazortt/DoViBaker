#include "DoViSceneFileReader.h"
#include "DoViProcessor.h"

DoViSceneFileReader::DoViSceneFileReader(
	PClip child,
	std::string sceneCllFile,
	IScriptEnvironment* env)
	: GenericVideoFilter(child)
	, currentScene(0)
	, previousFrame(0)
{
	uint32_t frame, lastFrameInScene, pq;
	uint16_t maxPq = 0;
	uint16_t staticMaxPq = 0;
	FILE* fp;

	fp = fopen(sceneCllFile.c_str(), "r");
	if (!fp) {
		env->ThrowError((std::string("DoViSceneFileReader: cannot find scenes file ") + sceneCllFile).c_str());
	}
	while (fscanf(fp, "%i %i %i\n", &frame, &lastFrameInScene, &pq) == 3) {
		if (pq > maxPq) maxPq = pq;
		if (maxPq > staticMaxPq) staticMaxPq = maxPq;
		if (!lastFrameInScene) continue;
		uint16_t maxCll = DoViProcessor::pq2nits(maxPq);
		sceneMaxSignal.push_back(std::tuple(frame + 1, maxPq, maxCll));
		maxPq = 0;
	}
	uint16_t maxCll = DoViProcessor::pq2nits(maxPq);
	sceneMaxSignal.push_back(std::tuple(frame + 1, maxPq, maxCll));
	staticMaxCll = DoViProcessor::pq2nits(staticMaxPq);
}

PVideoFrame DoViSceneFileReader::GetFrame(int n, IScriptEnvironment* env)
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