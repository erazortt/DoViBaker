#pragma once

#ifndef TIMECUBE_LUT_H_
#define TIMECUBE_LUT_H_

#include <memory>

struct timecube_filter {
protected:
	~timecube_filter() = default;
};

namespace timecube {

struct Cube;

enum class PixelType {
	BYTE,
	WORD,
	HALF,
	FLOAT,
};

struct PixelFormat {
	PixelType type;
	unsigned depth;
	bool fullrange;
};

class Lut : public ::timecube_filter {
public:
	virtual ~Lut() = default;

	virtual void to_float(const void * const src[3], float * const dst[3], const PixelFormat &format, unsigned width) const;

	virtual void from_float(const float * const src[3], void * const dst[3], const PixelFormat &format, unsigned width) const;

	virtual bool supports_half() const;

	virtual void process(const float * const src[3], float * const dst[3], unsigned width) const = 0;
};

std::unique_ptr<Lut> create_lut_impl(const Cube &cube, int simd);

} // namespace timecube

#endif // TIMECUBE_LUT_H_
