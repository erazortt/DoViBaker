#pragma once

/*
 * Export of the Avisynth+ plugin entry point.
 *
 * The host looks for "AvisynthPluginInit3" on Windows and for
 * "avisynth_plugin_init3" on POSIX systems, and the plugins of this project
 * are built with hidden visibility there, so the entry point needs the
 * default visibility explicitly.
 */

#include "avisynth.h"

#ifdef AVS_POSIX
#  define DOVI_PLUGIN_EXPORT __attribute__((visibility("default")))
#  define DOVI_PLUGIN_INIT   avisynth_plugin_init3
#else
#  define DOVI_PLUGIN_EXPORT __declspec(dllexport)
#  define DOVI_PLUGIN_INIT   AvisynthPluginInit3
#endif
