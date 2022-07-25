#ifdef CUBE_X86

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <immintrin.h>
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
		void *ptr = _mm_malloc(n * sizeof(T), 64);
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


void byte_to_float(const uint8_t *src, float *dst, float scale, float offset, unsigned width)
{
	const __m512 scale_ps = _mm512_set1_ps(scale);
	const __m512 offset_ps = _mm512_set1_ps(offset);

	for (unsigned i = 0; i < width - width % 16; i += 16) {
		__m512i x = _mm512_cvtepu8_epi32(_mm_load_si128((const __m128i *)(src + i)));
		__m512 y = _mm512_cvtepi32_ps(x);

		y = _mm512_fmadd_ps(scale_ps, y, offset_ps);
		_mm512_store_ps(dst + i, y);
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		dst[i] = _mm_cvtss_f32(_mm_fmadd_ss(_mm_set_ss(src[i]), _mm_set_ss(scale), _mm_set_ss(offset)));
	}
}

void word_to_float(const uint16_t *src, float *dst, float scale, float offset, unsigned width)
{
	const __m512 scale_ps = _mm512_set1_ps(scale);
	const __m512 offset_ps = _mm512_set1_ps(offset);

	for (unsigned i = 0; i < width - width % 16; i += 16) {
		__m512i x = _mm512_cvtepu16_epi32(_mm256_load_si256((const __m256i *)(src + i)));
		__m512 y = _mm512_cvtepi32_ps(x);

		y = _mm512_fmadd_ps(scale_ps, y, offset_ps);
		_mm512_store_ps(dst + i, y);
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		dst[i] = _mm_cvtss_f32(_mm_fmadd_ss(_mm_set_ss(src[i]), _mm_set_ss(scale), _mm_set_ss(offset)));
	}
}

void half_to_float(const uint16_t *src, float *dst, unsigned width)
{
	for (unsigned i = 0; i < width - width % 16; i += 16) {
		__m256i x = _mm256_load_si256((const __m256i *)(src + i));
		__m512 y = _mm512_cvtph_ps(x);
		_mm512_store_ps(dst + i, y);
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		dst[i] = _mm_cvtss_f32(_mm_cvtph_ps(_mm_set1_epi16(src[i])));
	}
}

void float_to_byte(const float *src, uint8_t *dst, float scale, float offset, unsigned width)
{
	const __m512 scale_ps = _mm512_set1_ps(scale);
	const __m512 offset_ps = _mm512_set1_ps(offset);

	for (unsigned i = 0; i < width - width % 16; i += 16) {
		__m512 x;
		__m512i xi;

		x = _mm512_load_ps(src + i);
		x = _mm512_fmadd_ps(scale_ps, x, offset_ps);
		xi = _mm512_cvtps_epi32(x);
		_mm_store_si128((__m128i *)(dst + i), _mm512_cvtusepi32_epi8(xi));
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		float x = _mm_cvtss_f32(_mm_fmadd_ss(_mm_set_ss(src[i]), _mm_set_ss(scale), _mm_set_ss(offset)));
		int y = _mm_cvt_ss2si(_mm_set_ss(x));
		y = std::min(std::max(y, 0), static_cast<int>(UINT8_MAX));
		dst[i] = static_cast<uint8_t>(y);
	}
}

void float_to_word(const float *src, uint16_t *dst, unsigned depth, float scale, float offset, unsigned width)
{
	const __m512 scale_ps = _mm512_set1_ps(scale);
	const __m512 offset_ps = _mm512_set1_ps(offset);

	for (unsigned i = 0; i < width - width % 16; i += 16) {
		__m512 x;
		__m512i xi;

		x = _mm512_load_ps(src + i);
		x = _mm512_fmadd_ps(scale_ps, x, offset_ps);
		xi = _mm512_cvtps_epi32(x);
		_mm256_store_si256((__m256i *)(dst + i), _mm256_min_epu16(_mm512_cvtusepi32_epi16(xi), _mm256_set1_epi16((1U << depth) - 1)));
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		float x = _mm_cvtss_f32(_mm_fmadd_ss(_mm_set_ss(src[i]), _mm_set_ss(scale), _mm_set_ss(offset)));
		int y = _mm_cvt_ss2si(_mm_set_ss(x));
		y = std::min(std::max(y, 0), static_cast<int>((1U << depth) - 1));
		dst[i] = static_cast<uint16_t>(y);
	}
}

void float_to_half(const float *src, uint16_t *dst, unsigned width)
{
	for (unsigned i = 0; i < width - width % 16; i += 16) {
		__m512 x = _mm512_load_ps(src + i);
		__m256i y = _mm512_cvtps_ph(x, 0);
		_mm256_store_si256((__m256i *)(dst + i), y);
	}
	for (unsigned i = width - width % 16; i < width; ++i) {
		int x = _mm_extract_epi16(_mm_cvtps_ph(_mm_set_ss(src[i]), 0), 0);
		dst[i] = static_cast<uint16_t>(x);
	}
}


static inline FORCE_INLINE __m512 mm512_interp_ps(__m512 lo, __m512 hi, __m512 dist)
{
	__m512 x;

	// (1 - x) * a == -x * a + a
	x = _mm512_fnmadd_ps(dist, lo, lo);
	// (-x * a + a) + x * b
	x = _mm512_fmadd_ps(dist, hi, x);

	return x;
}

// Transpose 2x2 128-bit elements within upper and lower 256-bit lanes.
static inline FORCE_INLINE void mm512_transpose2_ps128(__m512 &a, __m512 &b)
{
	__m512 tmp0 = a;
	__m512 tmp1 = b;
	a = _mm512_shuffle_f32x4(tmp0, tmp1, _MM_SHUFFLE(2, 0, 2, 0));
	b = _mm512_shuffle_f32x4(tmp0, tmp1, _MM_SHUFFLE(3, 1, 3, 1));
}

static inline FORCE_INLINE __m512i lut3d_calculate_index(const __m512 &r, const __m512 &g, const __m512 &b, const __m512i &stride_g, const __m512i &stride_b)
{
	__m512i idx_r, idx_g, idx_b;

	idx_r = _mm512_cvttps_epi32(r);
	idx_r = _mm512_slli_epi32(idx_r, 4); // 16 byte entries.

	idx_g = _mm512_cvttps_epi32(g);
	idx_g = _mm512_mullo_epi32(idx_g, stride_g);

	idx_b = _mm512_cvttps_epi32(b);
	idx_b = _mm512_mullo_epi32(idx_b, stride_b);

	return _mm512_add_epi32(_mm512_add_epi32(idx_r, idx_g), idx_b);
}


// Performs trilinear interpolation on four pixels.
// Returns [R0 G0 B0 xx R1 G1 B1 xx R2 G2 B2 x R3 G3 B3 x].
static inline FORCE_INLINE __m512 lut3d_trilinear_interp(const void *lut, ptrdiff_t stride_g, ptrdiff_t stride_b,
                                                         ptrdiff_t idx_lolo, ptrdiff_t idx_lohi, ptrdiff_t idx_hilo, ptrdiff_t idx_hihi,
                                                         __m512 r, __m512 g, __m512 b)
{
#define LUT_OFFSET(x) reinterpret_cast<const float *>(static_cast<const unsigned char *>(lut) + (x))
	__m512 g_lo = _mm512_shuffle_f32x4(g, g, _MM_SHUFFLE(1, 1, 0, 0));
	__m512 b_lo = _mm512_shuffle_f32x4(b, b, _MM_SHUFFLE(1, 1, 0, 0));
	__m512 g_hi = _mm512_shuffle_f32x4(g, g, _MM_SHUFFLE(3, 3, 2, 2));
	__m512 b_hi = _mm512_shuffle_f32x4(b, b, _MM_SHUFFLE(3, 3, 2, 2));

	__m512 g0b0_a, g0b1_a, g1b0_a, g1b1_a;
	__m512 g0b0_b, g0b1_b, g1b0_b, g1b1_b;

	g0b0_a = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_lolo))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_lohi)), 1);
	g1b0_a = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_lolo + stride_g))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_lohi + stride_g)), 1);
	g0b1_a = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_lolo + stride_b))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_lohi + stride_b)), 1);
	g1b1_a = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_lolo + stride_b + stride_g))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_lohi + stride_b + stride_g)), 1);

	g0b0_a = mm512_interp_ps(g0b0_a, g1b0_a, g_lo);
	g0b1_a = mm512_interp_ps(g0b1_a, g1b1_a, g_lo);

	g0b0_a = mm512_interp_ps(g0b0_a, g0b1_a, b_lo);

	g0b0_b = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_hilo))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_hihi)), 1);
	g1b0_b = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_hilo + stride_g))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_hihi + stride_g)), 1);
	g0b1_b = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_hilo + stride_b))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_hihi + stride_b)), 1);
	g1b1_b = _mm512_insertf32x8(_mm512_castps256_ps512(_mm256_loadu_ps(LUT_OFFSET(idx_hilo + stride_b + stride_g))),
	                            _mm256_loadu_ps(LUT_OFFSET(idx_hihi + stride_b + stride_g)), 1);

	g0b0_b = mm512_interp_ps(g0b0_b, g1b0_b, g_hi);
	g0b1_b = mm512_interp_ps(g0b1_b, g1b1_b, g_hi);

	g0b0_b = mm512_interp_ps(g0b0_b, g0b1_b, b_hi);

	mm512_transpose2_ps128(g0b0_a, g0b0_b);
	g0b0_a = mm512_interp_ps(g0b0_a, g0b0_b, r);

	return g0b0_a;
#undef LUT_OFFSET
}

// Converts packed [R0 G0 B0 xx R4 G4 B4 xx R8 G8 B8 xx Rc Gc Bc xx] ... to [R1 R2 ...] [G1 G2 ...] [B1 B2 ...].
static inline FORCE_INLINE void lut3d_unpack_result(const __m512 &result048c, const __m512 &result159d, const __m512 &result26ae, const __m512 &result37bf,
                                                    __m512 &r, __m512 &g, __m512 &b)
{
	__m512 t0 = _mm512_shuffle_ps(result048c, result159d, 0x44);
	__m512 t1 = _mm512_shuffle_ps(result26ae, result37bf, 0x44);
	__m512 t2 = _mm512_shuffle_ps(result048c, result159d, 0xEE);
	__m512 t3 = _mm512_shuffle_ps(result26ae, result37bf, 0xEE);

	__m512 tt0 = _mm512_shuffle_ps(t0, t1, 0x88); // r0 r1 r2 r3 | r4 r5 r6 r7
	__m512 tt1 = _mm512_shuffle_ps(t0, t1, 0xDD); // g0 g1 g2 g3 | g4 g5 g6 g7
	__m512 tt2 = _mm512_shuffle_ps(t2, t3, 0x88); // b0 b1 b2 b3 | b4 b5 b6 b7
	// __m512 tt3 = _mm512_shuffle_ps(t2, t3, 0xDD);

	r = tt0;
	g = tt1;
	b = tt2;
}

class Lut3D_AVX512 final : public Lut {
	std::vector<float, AlignedAllocator<float>> m_lut;
	uint_least32_t m_dim;
	float m_scale[3];
	float m_offset[3];
public:
	explicit Lut3D_AVX512(const Cube &cube) :
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
			float offset;

			if (format.fullrange) {
				scale = 1.0f / ((1UL << format.depth) - 1);
				offset = 0;
			} else {
				scale = 1.0f / (219UL << (format.depth - 8));
				offset = -static_cast<float>(16UL << (format.depth - 8)) * scale;
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
		} else if (format.type == PixelType::HALF) {
			half_to_float(static_cast<const uint16_t *>(src[0]), dst[0], width);
			half_to_float(static_cast<const uint16_t *>(src[1]), dst[1], width);
			half_to_float(static_cast<const uint16_t *>(src[2]), dst[2], width);
		} else {
			Lut::to_float(src, dst, format, width);
		}
	}

	void from_float(const float * const src[3], void * const dst[3], const PixelFormat &format, unsigned width) const override
	{
		if (format.type == PixelType::BYTE || format.type == PixelType::WORD) {
			float scale;
			float offset;

			if (format.fullrange) {
				scale = static_cast<float>((1UL << format.depth) - 1);
				offset = 0;
			} else {
				scale = static_cast<float>(219UL << (format.depth - 8));
				offset = static_cast<float>(16 << (format.depth - 8));
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
		} else if (format.type == PixelType::HALF) {
			float_to_half(src[0], static_cast<uint16_t *>(dst[0]), width);
			float_to_half(src[1], static_cast<uint16_t *>(dst[1]), width);
			float_to_half(src[2], static_cast<uint16_t *>(dst[2]), width);
		} else {
			Lut::from_float(src, dst, format, width);
		}
	}

	bool supports_half() const override { return true; }

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

		const __m512 scale_r = _mm512_set1_ps(m_scale[0]);
		const __m512 scale_g = _mm512_set1_ps(m_scale[1]);
		const __m512 scale_b = _mm512_set1_ps(m_scale[2]);
		const __m512 offset_r = _mm512_set1_ps(m_offset[0]);
		const __m512 offset_g = _mm512_set1_ps(m_offset[1]);
		const __m512 offset_b = _mm512_set1_ps(m_offset[2]);

		const __m512 lut_max = _mm512_set1_ps(std::nextafter(static_cast<float>(m_dim - 1), -INFINITY));
		const __m512i lut_stride_g_epi32 = _mm512_set1_epi32(lut_stride_g);
		const __m512i lut_stride_b_epi32 = _mm512_set1_epi32(lut_stride_b);

		for (unsigned i = 0; i < width; i += 16) {
			__m512 r = _mm512_load_ps(src_r + i);
			__m512 g = _mm512_load_ps(src_g + i);
			__m512 b = _mm512_load_ps(src_b + i);

			__m512 result048c, result159d, result26ae, result37bf;
			__m512 rtmp, gtmp, btmp;
			__m512i idx;
			__m128i idx_lolo, idx_lohi, idx_hilo, idx_hihi;

			size_t idx_scalar_lolo, idx_scalar_lohi, idx_scalar_hilo, idx_scalar_hihi;

			// Input domain remapping.
			r = _mm512_fmadd_ps(r, scale_r, offset_r);
			g = _mm512_fmadd_ps(g, scale_g, offset_g);
			b = _mm512_fmadd_ps(b, scale_b, offset_b);

			r = _mm512_max_ps(r, _mm512_setzero_ps());
			r = _mm512_min_ps(r, lut_max);

			g = _mm512_max_ps(g, _mm512_setzero_ps());
			g = _mm512_min_ps(g, lut_max);

			b = _mm512_max_ps(b, _mm512_setzero_ps());
			b = _mm512_min_ps(b, lut_max);

			// Base offset.
			idx = lut3d_calculate_index(r, g, b, lut_stride_g_epi32, lut_stride_b_epi32);
			idx_lolo = _mm512_castsi512_si128(idx);
			idx_lohi = _mm512_extracti32x4_epi32(idx, 1);
			idx_hilo = _mm512_extracti32x4_epi32(idx, 2);
			idx_hihi = _mm512_extracti32x4_epi32(idx, 3);

			// Cube distances.
			r = _mm512_sub_ps(r, _mm512_roundscale_ps(r, 1));
			g = _mm512_sub_ps(g, _mm512_roundscale_ps(g, 1));
			b = _mm512_sub_ps(b, _mm512_roundscale_ps(b, 1));

			// Interpolation.
#if SIZE_MAX >= UINT64_MAX
  #define EXTRACT_EVEN(out, x, idx) _mm_extract_epi64((x), (idx) / 2)
  #define EXTRACT_ODD(out, x, idx) ((out) >> 32)
#else
  #define EXTRACT_EVEN(out, x, idx) _mm_extract_epi32((x), (idx))
  #define EXTRACT_ODD(out, x, idx) _mm_extract_epi32((x), (idx))
#endif
			rtmp = _mm512_permute_ps(r, _MM_SHUFFLE(0, 0, 0, 0));
			gtmp = _mm512_permute_ps(g, _MM_SHUFFLE(0, 0, 0, 0));
			btmp = _mm512_permute_ps(b, _MM_SHUFFLE(0, 0, 0, 0));
			idx_scalar_lolo = EXTRACT_EVEN(idx_scalar_lolo, idx_lolo, 0);
			idx_scalar_lohi = EXTRACT_EVEN(idx_scalar_lohi, idx_lohi, 0);
			idx_scalar_hilo = EXTRACT_EVEN(idx_scalar_hilo, idx_hilo, 0);
			idx_scalar_hihi = EXTRACT_EVEN(idx_scalar_hihi, idx_hihi, 0);
			result048c = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b,
			                                    idx_scalar_lolo & 0xFFFFFFFFU, idx_scalar_lohi & 0xFFFFFFFFU, idx_scalar_hilo & 0xFFFFFFFFU, idx_scalar_hihi & 0xFFFFFFFFU,
			                                    rtmp, gtmp, btmp);

			rtmp = _mm512_permute_ps(r, _MM_SHUFFLE(1, 1, 1, 1));
			gtmp = _mm512_permute_ps(g, _MM_SHUFFLE(1, 1, 1, 1));
			btmp = _mm512_permute_ps(b, _MM_SHUFFLE(1, 1, 1, 1));
			idx_scalar_lolo = EXTRACT_ODD(idx_scalar_lolo, idx_lolo, 1);
			idx_scalar_lohi = EXTRACT_ODD(idx_scalar_lohi, idx_lohi, 1);
			idx_scalar_hilo = EXTRACT_ODD(idx_scalar_hilo, idx_hilo, 1);
			idx_scalar_hihi = EXTRACT_ODD(idx_scalar_hihi, idx_hihi, 1);
			result159d = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b,
			                                    idx_scalar_lolo, idx_scalar_lohi, idx_scalar_hilo, idx_scalar_hihi,
			                                    rtmp, gtmp, btmp);

			rtmp = _mm512_permute_ps(r, _MM_SHUFFLE(2, 2, 2, 2));
			gtmp = _mm512_permute_ps(g, _MM_SHUFFLE(2, 2, 2, 2));
			btmp = _mm512_permute_ps(b, _MM_SHUFFLE(2, 2, 2, 2));
			idx_scalar_lolo = EXTRACT_EVEN(idx_scalar_lolo, idx_lolo, 2);
			idx_scalar_lohi = EXTRACT_EVEN(idx_scalar_lohi, idx_lohi, 2);
			idx_scalar_hilo = EXTRACT_EVEN(idx_scalar_hilo, idx_hilo, 2);
			idx_scalar_hihi = EXTRACT_EVEN(idx_scalar_hihi, idx_hihi, 2);
			result26ae = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b,
			                                    idx_scalar_lolo & 0xFFFFFFFFU, idx_scalar_lohi & 0xFFFFFFFFU, idx_scalar_hilo & 0xFFFFFFFFU, idx_scalar_hihi & 0xFFFFFFFFU,
			                                    rtmp, gtmp, btmp);

			rtmp = _mm512_permute_ps(r, _MM_SHUFFLE(3, 3, 3, 3));
			gtmp = _mm512_permute_ps(g, _MM_SHUFFLE(3, 3, 3, 3));
			btmp = _mm512_permute_ps(b, _MM_SHUFFLE(3, 3, 3, 3));
			idx_scalar_lolo = EXTRACT_ODD(idx_scalar_lolo, idx_lolo, 3);
			idx_scalar_lohi = EXTRACT_ODD(idx_scalar_lohi, idx_lohi, 3);
			idx_scalar_hilo = EXTRACT_ODD(idx_scalar_hilo, idx_hilo, 3);
			idx_scalar_hihi = EXTRACT_ODD(idx_scalar_hihi, idx_hihi, 3);
			result37bf = lut3d_trilinear_interp(lut, lut_stride_g, lut_stride_b,
			                                    idx_scalar_lolo, idx_scalar_lohi, idx_scalar_hilo, idx_scalar_hihi,
			                                    rtmp, gtmp, btmp);
#undef EXTRACT_ODD
#undef EXTRACT_EVEN
			lut3d_unpack_result(result048c, result159d, result26ae, result37bf, r, g, b);

			if (i + 16 > width) {
				__mmask16 mask = 0xFFFFU >> (i + 16 - width);
				_mm512_mask_store_ps(dst_r + i, mask, r);
				_mm512_mask_store_ps(dst_g + i, mask, g);
				_mm512_mask_store_ps(dst_b + i, mask, b);
			} else {
				_mm512_store_ps(dst_r + i, r);
				_mm512_store_ps(dst_g + i, g);
				_mm512_store_ps(dst_b + i, b);
			}
		}
	}
};

} // namespace


std::unique_ptr<Lut> create_lut_impl_avx512(const Cube &cube)
{
	return cube.is_3d ? std::unique_ptr<Lut>(new Lut3D_AVX512{ cube }) : nullptr;
}

} // namespace timecube

#endif // CUBE_X86
