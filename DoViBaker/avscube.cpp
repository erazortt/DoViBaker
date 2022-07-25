#include "avisynth.h"
#include "cube.h"
#include "lut.h"
#include "VSHelper.h"
#include <io.h>

class AVSCube : public GenericVideoFilter
{
	std::unique_ptr<timecube::Lut> m_lut;
	char cube_path[1024];
	int cpu;
	bool fullrange;
	bool has_at_least_v8;

public:
	AVSCube(PClip _child, char *_cube_path, int _cpu, bool _fullrange, IScriptEnvironment* env) : GenericVideoFilter(_child)
	{
		if (vi.pixel_type != VideoInfo::CS_RGBP16)
		{
			env->ThrowError("AVSCube: input must be CS_RGBP16");
		}
		has_at_least_v8 = true;
		try { env->CheckVersion(8); }
		catch (const AvisynthError&) { has_at_least_v8 = false; }

		strcpy_s(cube_path, 1024, _cube_path);
		cpu = _cpu;
		fullrange = _fullrange;
		if (_access(cube_path, 0))
			env->ThrowError("AVSCube: cannot open cube file");
		timecube::Cube cube = timecube::read_cube_from_file(cube_path);
		m_lut = timecube::create_lut_impl(cube, cpu);
	}

	AVSCube::~AVSCube(void)
	{
	}

	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
};

template <class T>
T *incr(T *ptr, ptrdiff_t count)
{
	return (T *)((const char *)ptr + count);
}

PVideoFrame __stdcall AVSCube::GetFrame(int n, IScriptEnvironment* env)
{
	PVideoFrame src = child->GetFrame(n, env);
	//	PVideoFrame dst = env->NewVideoFrame(vi);
	PVideoFrame dst;
	if (has_at_least_v8)
		dst = env->NewVideoFrameP(vi, &src);
	else
		dst = env->NewVideoFrame(vi);
	unsigned int width = vi.width;
	unsigned int height = vi.height;

	std::unique_ptr<float, decltype(&vs_aligned_free)> tmp_buf{ nullptr, vs_aligned_free };
	unsigned aligned_width = width % 8 ? (width - width % 8) + 8 : width;

	const void *src_p[3];
	ptrdiff_t src_stride[3];
	void *dst_p[3];
	ptrdiff_t dst_stride[3];
	float *tmp[3] = { 0 };

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

	tmp_buf.reset(vs_aligned_malloc<float>(aligned_width * 3 * sizeof(float), 32));
	if (!tmp_buf)
		throw std::bad_alloc{};

	tmp[0] = tmp_buf.get();
	tmp[1] = tmp_buf.get() + aligned_width;
	tmp[2] = tmp_buf.get() + aligned_width * 2;

	timecube::PixelFormat format;
	format.type = (timecube::PixelType) 1;
	format.depth = 16;
	format.fullrange = fullrange;

	for (unsigned i = 0; i < height; ++i)
	{
		m_lut->to_float(src_p, tmp, format, width);
		m_lut->process(tmp, tmp, width);
		m_lut->from_float(tmp, dst_p, format, width);

		for (unsigned p = 0; p < 3; ++p)
		{
			src_p[p] = incr(src_p[p], src_stride[p]);
			dst_p[p] = incr(dst_p[p], dst_stride[p]);
		}
	}

	return dst;
}

AVSValue __cdecl Create_AVSCube(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	return new AVSCube(args[0].AsClip(),
		(char *)args[1].AsString(""),
		args[2].AsInt(INT_MAX),
		args[3].AsBool(true),
		env);
}

const AVS_Linkage *AVS_linkage = 0;
extern "C" __declspec(dllexport) const char* __stdcall
AvisynthPluginInit3(IScriptEnvironment* env, AVS_Linkage* vectors)
{
	AVS_linkage = vectors;
	env->AddFunction("Cube", "c[cube]s[cpu]i[fullrange]b", Create_AVSCube, 0);
	return 0;
}
