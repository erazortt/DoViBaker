#ifdef CUBE_X86

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <smmintrin.h>
#include "cube.h"
#include "lut.h"
#include "lut_x86.h"

#if defined(_MSC_VER)
  #define FORCE_INLINE __forceinline
#elif defined(__GNUC__)
  #define FORCE_INLINE __attribute__((always_inline))
#else
  #define FORCE_INLINE
#endif

namespace timecube {
namespace {

template <class T>
struct AlignedAllocator {
	typedef T value_type;

	AlignedAllocator() = default;

	template <class U>
	AlignedAllocator(const AlignedAllocator<U> &) noexcept {}

	T *allocate(size_t n) const
	{
		void *ptr = _mm_malloc(n * sizeof(T), 16);
		if (!ptr)
			throw std::bad_alloc{};
		return static_cast<T *>(ptr);
	}

	void deallocate(void *ptr, size_t) const noexcept
	{
		_mm_free(ptr);
	}

	bool operator==(const AlignedAllocator &) const noexcept { return true; }
	bool operator!=(const AlignedAllocator &) const noexcept { return false; }
};


void byte_to_float(const uint8_t *src, float *dst, float scale, unsigned offset, unsigned width)
{
	__m128 scale_ps = _mm_set_ps1(scale);
	__m128i offset_epu8 = _mm_set1_epi8(offset);

	for (unsigned i = 0; i < (width - width % 16); i += 16) {
		__m128i lolo, lohi, hilo, hihi;
		__m128i lo, hi;
		__m128i x;
		__m128 y;

		x = _mm_load_si128((const __m128i *)(src + i));
		x = _mm_subs_epu8(x, offset_epu8);

		lo = _mm_unpacklo_epi8(x, _mm_setzero_si128());
		hi = _mm_unpackhi_epi8(x, _mm_setzero_si128());

		lolo = _mm_unpacklo_epi16(lo, _mm_setzero_si128());
		lohi = _mm_unpackhi_epi16(lo, _mm_setzero_si128());
		hilo = _mm_unpacklo_epi16(hi, _mm_setzero_si128());
		hihi = _mm_unpackhi_epi16(hi, _mm_setzero_si128());

		y = _mm_cvtepi32_ps(lolo);
		y = _mm_mul_ps(y, scale_ps);
		_mm_store_ps(dst + i + 0, y);

		y = _mm_cvtepi32_ps(lohi);
		y = _mm_mul_ps(y, scale_ps);
		_mm_store_ps(dst + i + 4, y);

		y = _mm_cvtepi32_ps(hilo);
		y = _mm_mul_ps(y, scale_ps);
		_mm_store_ps(dst + i + 8, y);

		y = _mm_cvtepi32_ps(hihi);
		y = _mm_mul_ps(y, scale_ps);
		_mm_store_ps(dst + i + 12, y);
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		dst[i] = (src[i] - offset) * scale;
	}
}

void word_to_float(const uint16_t *src, float *dst, float scale, unsigned offset, unsigned width)
{
	__m128 scale_ps = _mm_set_ps1(scale);
	__m128i offset_epu16 = _mm_set1_epi16(offset);

	for (unsigned i = 0; i < (width - width % 8); i += 8) {
		__m128i lo, hi;
		__m128i x;
		__m128 y;

		x = _mm_load_si128((const __m128i *)(src + i));
		x = _mm_subs_epu16(x, offset_epu16);

		lo = _mm_unpacklo_epi16(x, _mm_setzero_si128());
		hi = _mm_unpackhi_epi16(x, _mm_setzero_si128());

		y = _mm_cvtepi32_ps(lo);
		y = _mm_mul_ps(y, scale_ps);
		_mm_store_ps(dst + i + 0, y);

		y = _mm_cvtepi32_ps(hi);
		y = _mm_mul_ps(y, scale_ps);
		_mm_store_ps(dst + i + 4, y);
	}
	for (unsigned i = width - width % 8; i < width; ++i) {
		dst[i] = (src[i] - offset) * scale;
	}
}

void float_to_byte(const float *src, uint8_t *dst, float scale, unsigned offset, unsigned width)
{
	__m128 scale_ps = _mm_set_ps1(scale);
	__m128i offset_epu8 = _mm_set1_epi8(offset);

	for (unsigned i = 0; i < (width - width % 16); i += 16) {
		__m128i lolo, lohi, hilo, hihi;
		__m128i lo, hi;
		__m128 x;
		__m128i y;

		x = _mm_load_ps(src + i + 0);
		x = _mm_mul_ps(x, scale_ps);
		lolo = _mm_cvtps_epi32(x);

		x = _mm_load_ps(src + i + 4);
		x = _mm_mul_ps(x, scale_ps);
		lohi = _mm_cvtps_epi32(x);

		x = _mm_load_ps(src + i + 8);
		x = _mm_mul_ps(x, scale_ps);
		hilo = _mm_cvtps_epi32(x);

		x = _mm_load_ps(src + i + 12);
		x = _mm_mul_ps(x, scale_ps);
		hihi = _mm_cvtps_epi32(x);

		lo = _mm_packus_epi32(lolo, lohi);
		hi = _mm_packus_epi32(hilo, hihi);

		y = _mm_packus_epi16(lo, hi);
		y = _mm_adds_epu8(y, offset_epu8);
		_mm_store_si128((__m128i *)(dst + i), y);
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		float x = src[i] * scale;
		int y = _mm_cvt_ss2si(_mm_set_ss(x)) + offset;
		y = std::min(std::max(y, 0), static_cast<int>(UINT8_MAX));
		dst[i] = static_cast<uint8_t>(y);
	}
}

void float_to_word(const float *src, uint16_t *dst, unsigned depth, float scale, unsigned offset, unsigned width)
{
	__m128 scale_ps = _mm_set_ps1(scale);
	__m128i offset_epu16 = _mm_set1_epi16(offset);

	for (unsigned i = 0; i < (width - width % 8); i += 8) {
		__m128i lo, hi;
		__m128 x;
		__m128i y;

		x = _mm_load_ps(src + i + 0);
		x = _mm_mul_ps(x, scale_ps);
		lo = _mm_cvtps_epi32(x);

		x = _mm_load_ps(src + i + 4);
		x = _mm_mul_ps(x, scale_ps);
		hi = _mm_cvtps_epi32(x);

		y = _mm_packus_epi32(lo, hi);
		y = _mm_adds_epu16(y, offset_epu16);
		y = _mm_min_epu16(y, _mm_set1_epi16((1U << depth) - 1));
		_mm_store_si128((__m128i *)(dst + i), y);
	}
	for (unsigned i = width - width % 8; i < width; ++i) {
		float x = src[i] * scale;
		int y = _mm_cvt_ss2si(_mm_set_ss(x)) + offset;
		y = std::min(std::max(y, 0), static_cast<int>((1U << depth) - 1));
		dst[i] = static_cast<uint16_t>(y);
	}
}


static inline FORCE_INLINE __m128 mm_interp_ps(__m128 lo, __m128 hi, __m128 dist)
{
	__m128 x;

	// (1 - x) * a == a - x * a
	x = _mm_mul_ps(dist, lo);
	x = _mm_sub_ps(lo, x);
	// (a - x * a) + x * b
	x = _mm_add_ps(_mm_mul_ps(dist, hi), x);

	return x;
}

static inline FORCE_INLINE __m128i lut3d_calculate_index(const __m128 &r, const __m128 &g, const __m128 &b, const __m128i &stride_g, const __m128i &stride_b)
{
	__m128i idx_r, idx_g, idx_b;

	idx_r = _mm_cvttps_epi32(r);
	idx_r = _mm_slli_epi32(idx_r, 4); // 16 byte entries.

	idx_g = _mm_cvttps_epi32(g);
	idx_g = _mm_mullo_epi32(idx_g, stride_g);

	idx_b = _mm_cvttps_epi32(b);
	idx_b = _mm_mullo_epi32(idx_b, stride_b);

	return _mm_add_epi32(_mm_add_epi32(idx_r, idx_g), idx_b);
}

// Performs trilinear interpolation on one pixel.
// Returns [R G B x].
static inline FORCE_INLINE __m128 lut3d_trilinear_interp(const void *lut, ptrdiff_t stride_g, ptrdiff_t stride_b, ptrdiff_t idx,
                                                         __m128 r, __m128 g, __m128 b)
{
#define LUT_OFFSET(x) reinterpret_cast<const float *>(static_cast<const unsigned char *>(lut) + (x))
	__m128 r0g0b0, r1g0b0, r0g1b0, r1g1b0, r0g0b1, r1g0b1, r0g1b1, r1g1b1;

	r0g0b0 = _mm_load_ps(LUT_OFFSET(idx));
	r1g0b0 = _mm_load_ps(LUT_OFFSET(idx + 16));
	r0g1b0 = _mm_load_ps(LUT_OFFSET(idx + stride_g));
	r1g1b0 = _mm_load_ps(LUT_OFFSET(idx + stride_g + 16));
	r0g0b1 = _mm_load_ps(LUT_OFFSET(idx + stride_b));
	r1g0b1 = _mm_load_ps(LUT_OFFSET(idx + stride_b + 16));
	r0g1b1 = _mm_load_ps(LUT_OFFSET(idx + stride_g + stride_b));
	r1g1b1 = _mm_load_ps(LUT_OFFSET(idx + stride_g + stride_b + 16));

	r0g0b0 = mm_interp_ps(r0g0b0, r1g0b0, r);
	r0g1b0 = mm_interp_ps(r0g1b0, r1g1b0, r);
	r0g0b1 = mm_interp_ps(r0g0b1, r1g0b1, r);
	r0g1b1 = mm_interp_ps(r0g1b1, r1g1b1, r);

	r0g0b0 = mm_interp_ps(r0g0b0, r0g1b0, g);
	r0g0b1 = mm_interp_ps(r0g0b1, r0g1b1, g);

	r0g0b0 = mm_interp_ps(r0g0b0, r0g0b1, b);
	return r0g0b0;
#undef LUT_OFFSET
}

class Lut3D_SSE41 final : public Lut {
	std::vector<float, AlignedAllocator<float>> m_lut;
	uint_least32_t m_dim;
	float m_scale[3];
	float m_offset[3];
public:
	explicit Lut3D_SSE41(const Cube &cube) :
		m_dim{ cube.n },
		m_scale{},
		m_offset{}
	{
		for (unsigned i = 0; i < 3; ++i) {
			m_scale[i] = (m_dim - 1) / (cube.domain_max[i] - cube.domain_min[i]);
			m_offset[i] = cube.domain_min[i] * m_scale[i];
		}

		// Pad each LUT entry to 16 bytes.
		m_lut.resize(m_dim * m_dim * m_dim * 4);

		for (size_t i = 0; i < m_lut.size() / 4; ++i) {
			m_lut[i * 4 + 0] = cube.lut[i * 3 + 0];
			m_lut[i * 4 + 1] = cube.lut[i * 3 + 1];
			m_lut[i * 4 + 2] = cube.lut[i * 3 + 2];
		}
	}

	void to_float(const void * const src[3], float * const dst[3], const PixelFormat &format, unsigned width) const override
	{
		if (format.type == PixelType::BYTE || format.type == PixelType::WORD) {
			float scale;
			unsigned offset;

			if (format.fullrange) {
				scale = 1.0f / ((1UL << format.depth) - 1);
				offset = 0;
			} else {
				scale = 1.0f / (219UL << (format.depth - 8));
				offset = 16 << (format.depth - 8);
			}

			if (format.type == PixelType::BYTE) {
				byte_to_float(static_cast<const uint8_t *>(src[0]), dst[0], scale, offset, width);
				byte_to_float(static_cast<const uint8_t *>(src[1]), dst[1], scale, offset, width);
				byte_to_float(static_cast<const uint8_t *>(src[2]), dst[2], scale, offset, width);
			} else {
				word_to_float(static_cast<const uint16_t *>(src[0]), dst[0], scale, offset, width);
				word_to_float(static_cast<const uint16_t *>(src[1]), dst[1], scale, offset, width);
				word_to_float(static_cast<const uint16_t *>(src[2]), dst[2], scale, offset, width);
			}
		} else {
			Lut::to_float(src, dst, format, width);
		}
	}

	void from_float(const float * const src[3], void * const dst[3], const PixelFormat &format, unsigned width) const override
	{
		if (format.type == PixelType::BYTE || format.type == PixelType::WORD) {
			float scale;
			unsigned offset;

			if (format.fullrange) {
				scale = static_cast<float>((1UL << format.depth) - 1);
				offset = 0;
			} else {
				scale = static_cast<float>(219UL << (format.depth - 8));
				offset = 16 << (format.depth - 8);
			}

			if (format.type == PixelType::BYTE) {
				float_to_byte(src[0], static_cast<uint8_t *>(dst[0]), scale, offset, width);
				float_to_byte(src[1], static_cast<uint8_t *>(dst[1]), scale, offset, width);
				float_to_byte(src[2], static_cast<uint8_t *>(dst[2]), scale, offset, width);
			} else {
				float_to_word(src[0], static_cast<uint16_t *>(dst[0]), format.depth, scale, offset, width);
				float_to_word(src[1], static_cast<uint16_t *>(dst[1]), format.depth, scale, offset, width);
				float_to_word(src[2], static_cast<uint16_t *>(dst[2]), format.depth, scale, offset, width);
			}
		} else {
			Lut::from_float(src, dst, format, width);
		}
	}

	void process(const float * const src[3], float * const dst[3], unsigned width) const override
	{
		const float *lut = m_lut.data();
		uint32_t lut_stride_g = m_dim * sizeof(float) * 4;
		uint32_t lut_stride_b = m_dim * m_dim * sizeof(float) * 4;

		const float *src_r = src[0];
		const float *src_g = src[1];
		const float *src_b = src[2];
		float *dst_r = dst[0];
		float *dst_g = dst[1];
		float *dst_b = dst[2];

		const __m128 scale_r = _mm_set_ps1(m_scale[0]);
		const __m128 scale_g = _mm_set_ps1(m_scale[1]);
		const __m128 scale_b = _mm_set_ps1(m_scale[2]);
		const __m128 offset_r = _mm_set_ps1(m_offset[0]);
		const __m128 offset_g = _mm_set_ps1(m_offset[1]);
		const __m128 offset_b = _mm_set_ps1(m_offset[2]);

		const __m128 lut_max = _mm_set_ps1(std::nextafter(static_cast<float>(m_dim - 1), -INFINITY));
		const __m128i lut_stride_g_epi32 = _mm_set1_epi32(lut_stride_g);
		const __m128i lut_stride_b_epi32 = _mm_set1_epi32(lut_stride_b);

		for (unsigned i = 0; i < width; i += 4) {
			__m128 r = _mm_load_ps(src_r + i);
			__m128 g = _mm_load_ps(src_g + i);
			__m128 b = _mm_load_ps(src_b + i);

			__m128 result0, result1, result2, result3;
			__m128 rtmp, gtmp, btmp;
			__m128i idx;

			size_t idx_scalar;

			// Input domain remapping.
			r = _mm_add_ps(_mm_mul_ps(scale_r, r), offset_r);
			g = _mm_add_ps(_mm_mul_ps(scale_g, g), offset_g);
			b = _mm_add_ps(_mm_mul_ps(scale_b, b), offset_b);

			r = _mm_max_ps(r, _mm_setzero_ps());
			r = _mm_min_ps(r, lut_max);

			g = _mm_max_ps(g, _mm_setzero_ps());
			g = _mm_min_ps(g, lut_max);

			b = _mm_max_ps(b, _mm_setzero_ps());
			b = _mm_min_ps(b, lut_max);

			idx = lut3d_calculate_index(r, g, b, lut_stride_g_epi32, lut_stride_b_epi32);

			// Cube distances.
			r = _mm_sub_ps(r, _mm_floor_ps(r));
			g = _mm_sub_ps(g, _mm_floor_ps(g));
			b = _mm_sub_ps(b, _mm_floor_ps(b));

			// Interpolation.
#if SIZE_MAX >= UINT64_MAX
  #define EXTRACT_EVEN(out, x, idx) _mm_extract_epi64((x), (idx) / 2)
  #define EXTRACT_ODD(out, x, idx) ((out) >> 32)
#else
  #define EXTRACT_EVEN(out, x, idx) _mm_extract_epi32((x), (idx))
  #define EXTRACT_ODD(out, x, idx) _mm_extract_epi32((x), (idx))
#endif
			rtmp = _mm_shuffle_ps(r, r, _MM_SHUFFLE(0, 0, 0, 0));
			gtmp = _mm_shuffle_ps(g, g, _MM_SHUFFLE(0, 0, 0, 0));
			btmp = _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 0, 0, 0));
			idx_scalar = EXTRACT_EVEN(idx_scalar, idx, 0);
			result0 = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b, idx_scalar & 0xFFFFFFFFU, rtmp, gtmp, btmp);

			rtmp = _mm_shuffle_ps(r, r, _MM_SHUFFLE(1, 1, 1, 1));
			gtmp = _mm_shuffle_ps(g, g, _MM_SHUFFLE(1, 1, 1, 1));
			btmp = _mm_shuffle_ps(b, b, _MM_SHUFFLE(1, 1, 1, 1));
			idx_scalar = EXTRACT_ODD(idx_scalar, idx, 1);
			result1 = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b, idx_scalar, rtmp, gtmp, btmp);

			rtmp = _mm_shuffle_ps(r, r, _MM_SHUFFLE(2, 2, 2, 2));
			gtmp = _mm_shuffle_ps(g, g, _MM_SHUFFLE(2, 2, 2, 2));
			btmp = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2, 2, 2, 2));
			idx_scalar = EXTRACT_EVEN(idx_scalar, idx, 2);
			result2 = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b, idx_scalar & 0xFFFFFFFFU, rtmp, gtmp, btmp);

			rtmp = _mm_shuffle_ps(r, r, _MM_SHUFFLE(3, 3, 3, 3));
			gtmp = _mm_shuffle_ps(g, g, _MM_SHUFFLE(3, 3, 3, 3));
			btmp = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 3, 3, 3));
			idx_scalar = EXTRACT_ODD(idx_scalar, idx, 3);
			result3 = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b, idx_scalar, rtmp, gtmp, btmp);
#undef EXTRACT_ODD
#undef EXTRACT_EVEN

			_MM_TRANSPOSE4_PS(result0, result1, result2, result3);
			r = result0;
			g = result1;
			b = result2;

			if (i + 4 > width) {
				alignas(16) float rbuf[4];
				alignas(16) float gbuf[4];
				alignas(16) float bbuf[4];
				_mm_store_ps(rbuf, r);
				_mm_store_ps(gbuf, g);
				_mm_store_ps(bbuf, b);

				for (unsigned ii = i; ii < width; ++ii) {
					dst_r[ii] = rbuf[ii - i];
					dst_g[ii] = gbuf[ii - i];
					dst_b[ii] = bbuf[ii - i];
				}
			} else {
				_mm_store_ps(dst_r + i, r);
				_mm_store_ps(dst_g + i, g);
				_mm_store_ps(dst_b + i, b);
			}
		}
	}
};

} // namespace


std::unique_ptr<Lut> create_lut_impl_sse41(const Cube &cube)
{
	return cube.is_3d ? std::unique_ptr<Lut>(new Lut3D_SSE41{ cube }) : nullptr;
}

} // namespace timecube

#endif // CUBE_X86
