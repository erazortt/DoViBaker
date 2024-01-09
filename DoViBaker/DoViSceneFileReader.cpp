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
	FILE* fp;

	fp = fopen(sceneCllFile.c_str(), "r");
	if (!fp) {
		env->ThrowError((std::string("DoViSceneFileReader: cannot find scenes file ") + sceneCllFile).c_str());
	}
	while (fscanf(fp, "%i %i %i\n", &frame, &lastFrameInScene, &pq) == 3) {
		if (pq > maxPq) maxPq = pq;
		if (!lastFrameInScene) continue;
		sceneMaxSignal.push_back(std::pair(frame + 1, maxPq));
		maxPq = 0;
	}
	sceneMaxSignal.push_back(std::pair(frame + 1, maxPq));
}

PVideoFrame DoViSceneFileReader::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);

	uint16_t maxCll;
	if (previousFrame > n)
		currentScene = 0;
	previousFrame = n;

	while (sceneMaxSignal.at(currentScene).first <= n) {
		currentScene++;
	}

	uint16_t maxPq = sceneMaxSignal.at(currentScene).second;
  uint16_t nits = DoViProcessor::pq2nits(maxPq);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_max_pq", maxPq, 0);
	env->propSetInt(env->getFramePropsRW(src), "_dovi_max_content_light_level", nits, 0);
	env->propSetInt(env->getFramePropsRW(src), "_SceneChangeNext", sceneMaxSignal.at(currentScene).first == n + 1, 0);
	if (currentScene > 0) {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", sceneMaxSignal.at(currentScene - 1).first == n, 0);
	} else {
		env->propSetInt(env->getFramePropsRW(src), "_SceneChangePrev", 0, 0);
	}

	return src;
}