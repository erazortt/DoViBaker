#pragma once

/*
 * This header selects *which* of the four products of this repository a
 * translation unit belongs to. The four products are:
 *
 *   DOVI_PRODUCT_BAKER     ->  DoViBaker   (Avisynth plugin: BL + EL + RPU -> PQ)
 *   DOVI_PRODUCT_ANALYZER  ->  DoViAnalyzer.exe (RPU.bin analyzer, console app)
 *   DOVI_PRODUCT_LUTGEN    ->  DoViLutGen.exe   (LUT generator, console app)
 *   DOVI_PRODUCT_AVSCUBE   ->  AVSCube     (Avisynth plugin: applies one LUT)
 *
 * The CMake build passes the selection on the compiler command line for every
 * target (e.g. -DDOVI_TARGET=DOVI_PRODUCT_ANALYZER), so all four products are
 * built in one go and nothing has to be edited here.
 *
 * The legacy MSVC project (DoViBaker/DoViBaker.vcxproj) passes nothing; then
 * the selection is taken from the "manual switch" below: uncomment exactly one
 * line and rebuild.
 */

#define DOVI_PRODUCT_BAKER    1
#define DOVI_PRODUCT_ANALYZER 2
#define DOVI_PRODUCT_LUTGEN   3
#define DOVI_PRODUCT_AVSCUBE  4

// ------- manual switch (only read when the build system sets nothing) -------
#if !defined(DOVI_TARGET)
#  define DOVI_TARGET DOVI_PRODUCT_BAKER
//#define DOVI_TARGET DOVI_PRODUCT_ANALYZER
//#define DOVI_TARGET DOVI_PRODUCT_LUTGEN
//#define DOVI_TARGET DOVI_PRODUCT_AVSCUBE
#endif
// ---------------------------------------------------------------------------

// Translate the selected product into the macro used by the sources.
#if   DOVI_TARGET == DOVI_PRODUCT_BAKER
#  define DOVI_BAKER
#elif DOVI_TARGET == DOVI_PRODUCT_ANALYZER
#  define DOVI_ANALYZER
#elif DOVI_TARGET == DOVI_PRODUCT_LUTGEN
#  define DOVI_LUTGEN
#elif DOVI_TARGET == DOVI_PRODUCT_AVSCUBE
#  define AVS_CUBE
#else
#  error "def.h: DOVI_TARGET must be DOVI_PRODUCT_BAKER, DOVI_PRODUCT_ANALYZER, DOVI_PRODUCT_LUTGEN or DOVI_PRODUCT_AVSCUBE"
#endif

// Exactly one product per compilation. Defining the product macros directly
// (the old way: -DDOVI_LUTGEN on the command line) keeps working, as long as
// it stays unique.
#if defined(DOVI_BAKER) && (defined(DOVI_ANALYZER) || defined(DOVI_LUTGEN) || defined(AVS_CUBE))
#  error "def.h: DOVI_BAKER conflicts with another product - select exactly one product (see comments in def.h)"
#elif defined(DOVI_ANALYZER) && (defined(DOVI_LUTGEN) || defined(AVS_CUBE))
#  error "def.h: DOVI_ANALYZER conflicts with another product - select exactly one product (see comments in def.h)"
#elif defined(DOVI_LUTGEN) && defined(AVS_CUBE)
#  error "def.h: DOVI_LUTGEN conflicts with another product - select exactly one product (see comments in def.h)"
#endif
