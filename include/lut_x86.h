#pragma once

#ifdef CUBE_X86

#ifndef TIMECUBE_LUT_X86_H_
#define TIMECUBE_LUT_X86_H_

#include <memory>

namespace timecube {

struct Cube;
class Lut;

std::unique_ptr<Lut> create_lut_impl_sse41(const Cube &cube);
std::unique_ptr<Lut> create_lut_impl_avx2(const Cube &cube);
std::unique_ptr<Lut> create_lut_impl_avx512(const Cube &cube);

std::unique_ptr<Lut> create_lut_impl_x86(const Cube &cube, int simd);

} // namespace timecube
#endif // TIMECUBE_LUT_X86_H_

#endif // CUBE_X86
