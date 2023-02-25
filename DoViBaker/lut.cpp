#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "cube.h"
#include "lut.h"
#include "lut_x86.h"

namespace timecube {
namespace {

template <class T>
struct ToFloat {
	float scale;
	int offset;

	explicit ToFloat(const PixelFormat &format) :
		scale{},
		offset{}
	{
		assert(format.type == PixelType::BYTE || format.type == PixelType::WORD);

		if (format.fullrange) {
			scale = 1.0f / ((1UL << format.depth) - 1);
			offset = 0;
		} else {
			assert(format.depth >= 8);
			scale = 1.0f / (219UL << (format.depth - 8));
			offset = 16 << (format.depth - 8);
		}
	}

	float operator()(T x)
	{
		return (static_cast<int32_t>(x) - offset) * scale;
	}
};

template <class T>
struct FromFloat {
	float scale;
	int offset;
	unsigned maxval;

	explicit FromFloat(const PixelFormat &format) :
		scale{},
		offset{},
		maxval{ (1U << format.depth) - 1 }
	{
		assert(format.type == PixelType::BYTE || format.type == PixelType::WORD);

		if (format.fullrange) {
			scale = static_cast<float>((1UL << format.depth) - 1);
			offset = 0;
		} else {
			assert(format.depth >= 8);
			scale = static_cast<float>(219UL << (format.depth - 8));
			offset = 16 << (format.depth - 8);
		}
	}

	T operator()(float x)
	{
		long val = std::lrint(x * scale) + offset;
		return static_cast<T>(std::min(std::max(val, 0L), static_cast<long>(maxval)));
	}
};


struct Vector3 : public std::array<float, 3> {
	Vector3() = default;

	Vector3(float a, float b, float c) : std::array<float, 3>{ { a, b, c } } {}
};

Vector3 operator+(const Vector3 &lhs, const Vector3 &rhs)
{
	return{ lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2] };
}

Vector3 operator*(float lhs, const Vector3 &rhs)
{
	return{ lhs * rhs[0], lhs * rhs[1], lhs * rhs[2] };
}

template <class T>
T interp(T lo, T hi, float dist)
{
	return (1.0f - dist) * lo + dist * hi;
}

Vector3 trilinear_interp(const Vector3 tri[2][2][2], float dist_x, float dist_y, float dist_z)
{
	Vector3 tmp0 = interp(tri[0][0][0], tri[1][0][0], dist_x);
	Vector3 tmp1 = interp(tri[0][0][1], tri[1][0][1], dist_x);
	Vector3 tmp2 = interp(tri[0][1][0], tri[1][1][0], dist_x);
	Vector3 tmp3 = interp(tri[0][1][1], tri[1][1][1], dist_x);

	tmp0 = interp(tmp0, tmp2, dist_y);
	tmp1 = interp(tmp1, tmp3, dist_y);

	tmp0 = interp(tmp0, tmp1, dist_z);
	return tmp0;
}


class Lut1D : public Lut {
	std::vector<float> m_lut[3];
	float m_scale[3];
	float m_offset[3];
public:
	explicit Lut1D(const Cube &cube) :
		m_scale{},
		m_offset{}
	{
		for (unsigned i = 0; i < 3; ++i) {
			m_lut[i].resize(cube.n);
			m_scale[i] = 1.0f / (cube.domain_max[i] - cube.domain_min[i]);
			m_offset[i] = cube.domain_min[i] * m_scale[i];
		}

		for (size_t i = 0; i < cube.n; ++i) {
			m_lut[0][i] = cube.lut[i * 3 + 0];
			m_lut[1][i] = cube.lut[i * 3 + 1];
			m_lut[2][i] = cube.lut[i * 3 + 2];
		}
	}

	void process(const float * const src[3], float * const dst[3], unsigned width) const override
	{
		uint_least32_t lut_max = static_cast<uint_least32_t>(m_lut[0].size() - 1);
		float lut_clamp = std::nextafterf(static_cast<float>(lut_max), -INFINITY);

		for (unsigned p = 0; p < 3; ++p) {
			const float *src_p = src[p];
			float *dst_p = dst[p];

			for (unsigned i = 0; i < width; ++i) {
				float x, d;
				uint_least32_t idx;

				x = src_p[i];
				x = (x * m_scale[p] + m_offset[p]) * lut_max;

				x = std::min(std::max(x, 0.0f), lut_clamp);
				idx = static_cast<uint_least32_t>(x);
				d = x - idx;

				x = interp(m_lut[p][idx], m_lut[p][idx + 1], d);
				dst_p[i] = x;
			}
		}
	}
};

class Lut3D : public Lut {
	std::vector<Vector3> m_lut;
	uint_least32_t m_dim;
	float m_scale[3];
	float m_offset[3];
public:
	explicit Lut3D(const Cube &cube) :
		m_dim{ cube.n },
		m_scale{},
		m_offset{}
	{
		for (unsigned i = 0; i < 3; ++i) {
			m_scale[i] = 1.0f / (cube.domain_max[i] - cube.domain_min[i]);
			m_offset[i] = cube.domain_min[i] * m_scale[i];
		}

		m_lut.resize(m_dim * m_dim * m_dim);

		for (size_t i = 0; i < m_lut.size(); ++i) {
			m_lut[i][0] = cube.lut[i * 3 + 0];
			m_lut[i][1] = cube.lut[i * 3 + 1];
			m_lut[i][2] = cube.lut[i * 3 + 2];
		}
	}

	void process(const float * const src[3], float * const dst[3], unsigned width) const override
	{
		const float *src_r = src[0];
		const float *src_g = src[1];
		const float *src_b = src[2];
		float *dst_r = dst[0];
		float *dst_g = dst[1];
		float *dst_b = dst[2];

		uint_least32_t lut_max = m_dim - 1;
		float lut_clamp = std::nextafter(static_cast<float>(lut_max), -INFINITY);

		for (unsigned i = 0; i < width; ++i) {
			float r, g, b;
			float dist_r, dist_g, dist_b;
			uint_least32_t idx_r, idx_g, idx_b;

			Vector3 tri[2][2][2];
			Vector3 interp_result;

			r = src_r[i];
			g = src_g[i];
			b = src_b[i];

			r = (r * m_scale[0] + m_offset[0]) * lut_max;
			g = (g * m_scale[1] + m_offset[1]) * lut_max;
			b = (b * m_scale[2] + m_offset[2]) * lut_max;

			r = std::min(std::max(r, 0.0f), lut_clamp);
			g = std::min(std::max(g, 0.0f), lut_clamp);
			b = std::min(std::max(b, 0.0f), lut_clamp);

			idx_r = static_cast<uint_least32_t>(r);
			idx_g = static_cast<uint_least32_t>(g);
			idx_b = static_cast<uint_least32_t>(b);

			dist_r = r - idx_r;
			dist_g = g - idx_g;
			dist_b = b - idx_b;

			tri[0][0][0] = m_lut[(idx_r + 0) + (idx_g + 0) * m_dim + (idx_b + 0) * m_dim * m_dim];
			tri[0][0][1] = m_lut[(idx_r + 1) + (idx_g + 0) * m_dim + (idx_b + 0) * m_dim * m_dim];
			tri[0][1][0] = m_lut[(idx_r + 0) + (idx_g + 1) * m_dim + (idx_b + 0) * m_dim * m_dim];
			tri[0][1][1] = m_lut[(idx_r + 1) + (idx_g + 1) * m_dim + (idx_b + 0) * m_dim * m_dim];
			tri[1][0][0] = m_lut[(idx_r + 0) + (idx_g + 0) * m_dim + (idx_b + 1) * m_dim * m_dim];
			tri[1][0][1] = m_lut[(idx_r + 1) + (idx_g + 0) * m_dim + (idx_b + 1) * m_dim * m_dim];
			tri[1][1][0] = m_lut[(idx_r + 0) + (idx_g + 1) * m_dim + (idx_b + 1) * m_dim * m_dim];
			tri[1][1][1] = m_lut[(idx_r + 1) + (idx_g + 1) * m_dim + (idx_b + 1) * m_dim * m_dim];

			interp_result = trilinear_interp(tri, dist_b, dist_g, dist_r);
			r = interp_result[0];
			g = interp_result[1];
			b = interp_result[2];

			dst_r[i] = r;
			dst_g[i] = g;
			dst_b[i] = b;
		}
	}
};

} // namespace



void Lut::to_float(const void * const src[3], float * const dst[3], const PixelFormat &format, unsigned width) const
{
	switch (format.type) {
	case PixelType::BYTE:
		std::transform(static_cast<const uint8_t *>(src[0]), static_cast<const uint8_t *>(src[0]) + width, dst[0], ToFloat<uint8_t>{ format });
		std::transform(static_cast<const uint8_t *>(src[1]), static_cast<const uint8_t *>(src[1]) + width, dst[1], ToFloat<uint8_t>{ format });
		std::transform(static_cast<const uint8_t *>(src[2]), static_cast<const uint8_t *>(src[2]) + width, dst[2], ToFloat<uint8_t>{ format });
		break;
	case PixelType::WORD:
		std::transform(static_cast<const uint16_t *>(src[0]), static_cast<const uint16_t *>(src[0]) + width, dst[0], ToFloat<uint16_t>{ format });
		std::transform(static_cast<const uint16_t *>(src[1]), static_cast<const uint16_t *>(src[1]) + width, dst[1], ToFloat<uint16_t>{ format });
		std::transform(static_cast<const uint16_t *>(src[2]), static_cast<const uint16_t *>(src[2]) + width, dst[2], ToFloat<uint16_t>{ format });
		break;
	case PixelType::HALF:
		throw std::runtime_error{ "half precision not implemented" };
	case PixelType::FLOAT:
		std::copy_n(static_cast<const float *>(src[0]), width, dst[0]);
		std::copy_n(static_cast<const float *>(src[1]), width, dst[1]);
		std::copy_n(static_cast<const float *>(src[2]), width, dst[2]);
		break;
	default:
		throw std::logic_error{ "bad pixel type" };
	}
}

void Lut::from_float(const float * const src[3], void * const dst[3], const PixelFormat &format, unsigned width) const
{
	switch (format.type) {
	case PixelType::BYTE:
		std::transform(src[0], src[0] + width, static_cast<uint8_t *>(dst[0]), FromFloat<uint8_t>{ format });
		std::transform(src[1], src[1] + width, static_cast<uint8_t *>(dst[1]), FromFloat<uint8_t>{ format });
		std::transform(src[2], src[2] + width, static_cast<uint8_t *>(dst[2]), FromFloat<uint8_t>{ format });
		break;
	case PixelType::WORD:
		std::transform(src[0], src[0] + width, static_cast<uint16_t *>(dst[0]), FromFloat<uint16_t>{ format });
		std::transform(src[1], src[1] + width, static_cast<uint16_t *>(dst[1]), FromFloat<uint16_t>{ format });
		std::transform(src[2], src[2] + width, static_cast<uint16_t *>(dst[2]), FromFloat<uint16_t>{ format });
		break;
	case PixelType::HALF:
		throw std::runtime_error{ "half precision not implemented" };
	case PixelType::FLOAT:
		std::copy_n(src[0], width, static_cast<float *>(dst[0]));
		std::copy_n(src[1], width, static_cast<float *>(dst[1]));
		std::copy_n(src[2], width, static_cast<float *>(dst[2]));
		break;
	default:
		throw std::logic_error{ "bad pixel type" };
	}
}

bool Lut::supports_half() const
{
	return false;
}

std::unique_ptr<Lut> create_lut_impl(const Cube &cube, int simd)
{
	std::unique_ptr<Lut> ret;

#ifdef CUBE_X86
	ret = create_lut_impl_x86(cube, simd);
#endif

	if (!ret) {
		if (cube.is_3d)
			ret = std::unique_ptr<Lut>(new Lut3D{ cube });
		else
			ret = std::unique_ptr<Lut>(new Lut1D{ cube });
	}

	return ret;
}

} // namespace timecube
