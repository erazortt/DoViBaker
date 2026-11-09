// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA, or visit
// http://www.gnu.org/copyleft/gpl.html .
//
// Linking Avisynth statically or dynamically with other modules is making a
// combined work based on Avisynth.  Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Avisynth give you
// permission to link Avisynth with independent modules that communicate with
// Avisynth solely through the interfaces defined in avisynth.h, regardless of the license
// terms of these independent modules, and to copy and distribute the
// resulting combined work under terms of your choice, provided that
// every copy of the combined work is accompanied by a complete copy of
// the source code of Avisynth (the version of Avisynth used to produce the
// combined work), being distributed under the terms of the GNU General
// Public License plus this exception.  An independent module is a module
// which is not derived from or based on Avisynth, such as 3rd-party filters,
// import and export plugins, or graphical user interfaces.

#ifndef AVSCORE_CPUID_H
#define AVSCORE_CPUID_H

#include <cstdint>
#include <cstddef>

// For GetCPUFlags/GetCPUFlagsEx.
// These are backwards-compatible with those in VirtualDub ending with SSE4_2
// For emulation see https://software.intel.com/en-us/articles/intel-software-development-emulator

// enum changed to constexpr since we need 64-bit constants and enum types are limited to int size
#ifndef AVS_CPU_FLAG_CONSTANT

#if (defined(__cplusplus) && __cplusplus >= 201103L) || \
    (defined(_MSC_VER) && _MSC_VER >= 1700) // MSVC 2012+ supports basic constexpr
    // Use constexpr for C++11 and later (resolves unused variable warnings)
#define AVS_CPU_FLAG_CONSTANT constexpr
#else
    // Fallback to static const for pre-C++11 legacy compilers
#define AVS_CPU_FLAG_CONSTANT static const
#endif

#endif // CPU_FLAG_CONSTANT

// Intel/AMD x86/x86-64 flags

AVS_CPU_FLAG_CONSTANT int64_t CPUF_FORCE = 0x01;    // N/A
AVS_CPU_FLAG_CONSTANT int64_t CPUF_FPU = 0x02;    // 386/486DX
AVS_CPU_FLAG_CONSTANT int64_t CPUF_MMX = 0x04;    // P55C, K6, PII
AVS_CPU_FLAG_CONSTANT int64_t CPUF_INTEGER_SSE = 0x08;    // PIII, Athlon
AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSE = 0x10;    // PIII, Athlon XP/MP
AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSE2 = 0x20;    // PIV, K8
AVS_CPU_FLAG_CONSTANT int64_t CPUF_3DNOW = 0x40;    // K6-2
AVS_CPU_FLAG_CONSTANT int64_t CPUF_3DNOW_EXT = 0x80;    // Athlon
AVS_CPU_FLAG_CONSTANT int64_t CPUF_X86_64 = 0xA0;    // Hammer

AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSE3 = 0x100;   // PIV+, K8 Venice
AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSSE3 = 0x200;   // Core 2
AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSE4 = 0x400;
AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSE4_1 = 0x400;   // Penryn, Wolfdale, Yorkfield
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX = 0x800;   // Sandy Bridge, Bulldozer
AVS_CPU_FLAG_CONSTANT int64_t CPUF_SSE4_2 = 0x1000;  // Nehalem

AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX2 = 0x2000;  // Haswell
AVS_CPU_FLAG_CONSTANT int64_t CPUF_FMA3 = 0x4000;
AVS_CPU_FLAG_CONSTANT int64_t CPUF_F16C = 0x8000;
AVS_CPU_FLAG_CONSTANT int64_t CPUF_MOVBE = 0x10000; // Big Endian move
AVS_CPU_FLAG_CONSTANT int64_t CPUF_POPCNT = 0x20000;
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AES = 0x40000;
AVS_CPU_FLAG_CONSTANT int64_t CPUF_FMA4 = 0x80000;

// AVX-512
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512F     = 0x00100000; // F Foundation.
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512DQ    = 0x00200000; // DQ (Double/Quad granular) Instructions
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512PF    = 0x00400000; // PF Prefetch
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512ER    = 0x00800000; // ER Exponential and Reciprocal
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512CD    = 0x01000000; // CD Conflict Detection
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512BW    = 0x02000000; // BW (Byte/Word granular) Instructions
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512VL    = 0x04000000; // VL (128/256 Vector Length) Extensions
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512IFMA  = 0x08000000; // IFMA integer 52 bit
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512VBMI  = 0x10000000; // VBMI, byte/word shuffling, sign/zero extension, and general pixel manipulation

// Group feature flags for convenience: checking a single flag for "base" and "fast" AVX512 feature sets.
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512_BASE = 0x20000000; // F, CD, BW, DQ, VL all set.
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512_FAST = 0x40000000; // Base + VNNI, VBMI, VBMI2, BITALG, VPOPCNTDQ. Spec detection logic excludes older/throttling models that also have these features.

// Last 32-bit flag reserved for future use:
// AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX10       = 0x80000000LL; // AVX10 as one flag, version query needed in distinct function.

// Flags exceeding the 32-bit limit (0xFFFFFFFF):
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512VNNI         = 0x00100000000LL; // VNNI, accumulated dot product on 8/16 bit integers
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512VBMI2        = 0x00200000000LL; // VBMI2: Byte/word load, store, & concatenation with shift for unaligned memory and packed data re-arrangement.
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512BITALG       = 0x00400000000LL; // BITALG, Bit Manipulation Instructions
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512VPOPCNTDQ    = 0x00800000000LL; // VPOPCNTDQ, Vector Population Count Double/Quadword
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512FP16         = 0x01000000000LL; // FP16, C++23 std::float16_t, S1E5M10, Limited range, higher precision (~3.3 decimal digits)
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512BF16         = 0x02000000000LL; // BF16, C++23 std::bfloat16_t, S1E8M7, Wide range (same as float); low precision (~2.3 decimal digits)

// We intentionally don't pollute the flag space with AVX512 crypto (VAES, VPCLMULQDQ, GFNI) and deprecated (VP2INTERSECT, 4VNNIW, 4FMAPS) features

// ARMv8-A flags (Values shared with X86 flags, usage guarded by platform macros)
AVS_CPU_FLAG_CONSTANT int64_t CPUF_ARM_NEON = 0x01; // NEON flag, minimum for aarch64
AVS_CPU_FLAG_CONSTANT int64_t CPUF_ARM_DOTPROD = 0x02; // Dot Product
AVS_CPU_FLAG_CONSTANT int64_t CPUF_ARM_SVE2 = 0x04; // SVE2
AVS_CPU_FLAG_CONSTANT int64_t CPUF_ARM_I8MM = 0x08; // I8MM
AVS_CPU_FLAG_CONSTANT int64_t CPUF_ARM_SVE2_1 = 0x10; // SVE2.1

#ifdef BUILDING_AVSCORE

// composite flags for feature group flags.
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512_BASE_ALL =
CPUF_AVX512F | CPUF_AVX512CD | CPUF_AVX512BW | CPUF_AVX512DQ | CPUF_AVX512VL |
CPUF_AVX512_BASE; // and the base single feature flag itself

AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512_FAST_ALL =
CPUF_AVX512_BASE_ALL |
CPUF_AVX512VNNI | CPUF_AVX512VBMI | CPUF_AVX512VBMI2 | CPUF_AVX512BITALG | CPUF_AVX512VPOPCNTDQ |
CPUF_AVX512_FAST; // and the FAST single feature flag itself

// Mask for all AVX-512 features listed here
AVS_CPU_FLAG_CONSTANT int64_t CPUF_AVX512_MASK =
CPUF_AVX512_FAST_ALL | CPUF_AVX512IFMA | CPUF_AVX512BF16 | CPUF_AVX512FP16 | CPUF_AVX512PF | CPUF_AVX512ER;

int GetCPUFlags();
int64_t GetCPUFlagsEx(); // more CPU flags, 32-bit was not enough
size_t GetL2CacheSize(); // in bytes
#endif

#endif // AVSCORE_CPUID_H
