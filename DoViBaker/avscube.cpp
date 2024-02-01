#include "avisynth.h"
#include <array>
#include <filesystem>
#include <io.h>
#include <timecube/timecube.h>

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

class AVSCube : public GenericVideoFilter
{
	struct TimecubeLutFree {
		void operator()(timecube_lut* ptr) { timecube_lut_free(ptr); }
	};

	struct TimecubeFilterFree {
		void operator()(timecube_filter* ptr) { timecube_filter_free(ptr); }
	};

	std::unique_ptr<timecube_filter, TimecubeFilterFree> m_lut;
	char cube_path[1024];
	int cpu;
	bool fullrange;
	bool has_at_least_v8;

public:
	AVSCube(PClip _child, char* _cube_path, int _cpu, bool _fullrange, IScriptEnvironment* env) : GenericVideoFilter(_child)
	{
		if (!vi.IsPlanarRGB())
		{
			env->ThrowError("AVSCube: input must be planar RGB");
		}
		has_at_least_v8 = true;
		try { env->CheckVersion(8); }
		catch (const AvisynthError&) { has_at_least_v8 = false; }

		timecube_filter_params params{};
		params.width = vi.width;
		params.height = vi.height;
		params.src_type = TIMECUBE_PIXEL_WORD;
		params.src_depth = vi.BitsPerComponent();
		params.src_range = TIMECUBE_RANGE_FULL;
		params.dst_type = TIMECUBE_PIXEL_WORD;
		params.dst_depth = vi.BitsPerComponent();
		params.dst_range = _fullrange ? TIMECUBE_RANGE_FULL : TIMECUBE_RANGE_LIMITED;
		params.interp = TIMECUBE_INTERP_TETRA;
		params.cpu = static_cast<timecube_cpu_type_e>(_cpu);

		strcpy_s(cube_path, 1024, _cube_path);
		cpu = _cpu;
		fullrange = _fullrange;
		if (_access(cube_path, 0))
			env->ThrowError("AVSCube: cannot open cube file");

		std::unique_ptr<timecube_lut, TimecubeLutFree> cube{ timecube_lut_from_file(cube_path) };
		if (!cube)
			throw std::runtime_error{ "AVSCube: error reading LUT from file" };

		m_lut.reset(timecube_filter_create(cube.get(), &params));
		if (!m_lut)
			throw std::runtime_error{ "AVSCube: error creating LUT" };
	}

	AVSCube::~AVSCube(void)
	{
	}

	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};

template <class T>
T* incr(T* ptr, ptrdiff_t count)
{
	return (T*)((const char*)ptr + count);
}

PVideoFrame __stdcall AVSCube::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);
	PVideoFrame dst;
	if (has_at_least_v8)
		dst = env->NewVideoFrameP(vi, &src);
	else
		dst = env->NewVideoFrame(vi);

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
	tmp.reset(aligned_malloc(timecube_filter_get_tmp_size(m_lut.get()), 64));

	timecube_filter_apply(m_lut.get(), src_p, src_stride, dst_p, dst_stride, tmp.get());

	return dst;
}

AVSValue __cdecl Create_AVSCube(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	return new AVSCube(args[0].AsClip(),
		(char*)args[1].AsString(""),
		args[2].AsInt(INT_MAX),
		args[3].AsBool(true),
		env);
}

const AVS_Linkage* AVS_linkage = 0;
extern "C" __declspec(dllexport) const char* __stdcall
AvisynthPluginInit3(IScriptEnvironment * env, AVS_Linkage * vectors)
{
	AVS_linkage = vectors;
	env->AddFunction("AVSCube", "c[cube]s[cpu]i[fullrange]b", Create_AVSCube, 0);
	return 0;
}
