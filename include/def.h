#pragma once

// uncomment exactly one of the following lines to select which application or library to build
#define DOVI_BAKER
//#define DOVI_ANALYZER
//#define DOVI_LUTGEN
//#define AVS_CUBE

//#ifdef DOVI_BAKER
//#define DOVI_ANALYZER
//#endif

#ifdef DOVI_ANALYZER
#define DOVI_BAKER
#endif

#ifdef DOVI_LUTGEN
#define DOVI_BAKER
#endif

#ifdef AVS_CUBE
#define DOVI_ANALYZER
#endif
