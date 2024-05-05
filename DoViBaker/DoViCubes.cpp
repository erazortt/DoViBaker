#include "DoViCubes.h"
#include <climits>
#include <filesystem>

AVS_FORCEINLINE void* aligned_malloc(size_t size, size_t align)
{
	void* result = [&]() {
#ifdef _MSC_VER 
		return _aligned_malloc(size, align);
#else 
		if (posix_memalign(&result, align, size))
			return result = nullptr;
		else
			return result;
#endif
		}();

		return result;
}

AVS_FORCEINLINE void aligned_free(void* ptr)
{
#ifdef _MSC_VER 
	_aligned_free(ptr);
#else 
	free(ptr);
#endif
}

DoViCubes::DoViCubes(
  PClip child,
  std::vector<std::pair<uint16_t, std::string>>& cubes,
  bool fullrange,
  IScriptEnvironment* env)
  : GenericVideoFilter(child)
{
	int lutMaxCpuCaps = INT_MAX;

	timecube_filter_params params{};
	params.width = vi.width;
	params.height = vi.height;
	params.src_type = TIMECUBE_PIXEL_WORD;
	params.src_depth = vi.BitsPerComponent();
	params.src_range = TIMECUBE_RANGE_FULL;
	params.dst_type = TIMECUBE_PIXEL_WORD;
	params.dst_depth = vi.BitsPerComponent();
	params.dst_range = fullrange ? TIMECUBE_RANGE_FULL : TIMECUBE_RANGE_LIMITED;
	params.interp = TIMECUBE_INTERP_TETRA;
	params.cpu = static_cast<timecube_cpu_type_e>(lutMaxCpuCaps);

	for (int i = 0; i < cubes.size(); i++) {
		auto cube_path = cubes[i].second;
		if (!std::filesystem::exists(std::filesystem::path(cube_path))) {
			env->ThrowError((std::string("DoViCubes: cannot find cube file ") + cube_path).c_str());
		}
		std::unique_ptr<timecube_lut, TimecubeLutFree> cube{ timecube_lut_from_file(cube_path.c_str()) };
		if (!cube) {
			env->ThrowError((std::string("DoViCubes: error reading LUT from file ") + cube_path).c_str());
		}

		timecube_filter* lut = timecube_filter_create(cube.get(), &params);
		if (!lut) {
			env->ThrowError((std::string("DoViCubes: error creating LUT from file ") + cube_path).c_str());
		}

		luts.push_back(std::pair(cubes[i].first, lut));
	}
}

DoViCubes::~DoViCubes() {}

void DoViCubes::applyLut(PVideoFrame& dst, const PVideoFrame& src) const
{
	const void* src_p[3];
	ptrdiff_t src_stride[3];
	void* dst_p[3];
	ptrdiff_t dst_stride[3];

	src_p[0] = src->GetReadPtr(PLANAR_R);
	src_p[1] = src->GetReadPtr(PLANAR_G);
	src_p[2] = src->GetReadPtr(PLANAR_B);
	src_stride[0] = src->GetPitch(PLANAR_R);
	src_stride[1] = src->GetPitch(PLANAR_G);
	src_stride[2] = src->GetPitch(PLANAR_B);
	dst_p[0] = dst->GetWritePtr(PLANAR_R);
	dst_p[1] = dst->GetWritePtr(PLANAR_G);
	dst_p[2] = dst->GetWritePtr(PLANAR_B);
	dst_stride[0] = dst->GetPitch(PLANAR_R);
	dst_stride[1] = dst->GetPitch(PLANAR_G);
	dst_stride[2] = dst->GetPitch(PLANAR_B);

	std::unique_ptr<void, decltype(&aligned_free)> tmp{ nullptr, aligned_free };
	tmp.reset(aligned_malloc(timecube_filter_get_tmp_size(currentFrameLut), 64));

	timecube_filter_apply(currentFrameLut, src_p, src_stride, dst_p, dst_stride, tmp.get());
}

PVideoFrame DoViCubes::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst = env->NewVideoFrameP(vi, &src);

	uint16_t maxCll = 0;
	if (env->propNumElements(env->getFramePropsRO(src), "_dovi_dynamic_max_content_light_level") > -1) {
		maxCll = env->propGetInt(env->getFramePropsRO(src), "_dovi_dynamic_max_content_light_level", 0, 0);
	}	else env->ThrowError("DoViCubes: Expected frame property not available");
	
	currentFrameLut = luts[luts.size() - 1].second;
	for (int i = 1; i < luts.size(); i++) {
		if (maxCll <= luts[i].first) {
			currentFrameLut = luts[i - 1].second;
			break;
		}
	}
	applyLut(dst, src);

	return dst;
}