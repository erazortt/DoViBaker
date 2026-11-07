# Common compile settings for every product of this repository, collected in
# the interface target dovibaker::settings.

option(DOVIBAKER_MSVC_STATIC_RUNTIME
    "Link the MSVC runtime statically (/MT), so the plugins need no VC redist" ON)
option(DOVIBAKER_FAST_MATH
    "Use fast (not fully standard conformant) floating point math in optimized builds" ${MSVC})

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
if(NOT MSVC)
    # The plugins are loaded by the host, the static helpers are linked into them.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(MSVC AND DOVIBAKER_MSVC_STATIC_RUNTIME)
    # Matches the historic MSVC project settings (Release: /MT, Debug: /MTd).
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

add_library(dovibaker_settings INTERFACE)
add_library(dovibaker::settings ALIAS dovibaker_settings)

target_compile_features(dovibaker_settings INTERFACE cxx_std_20)

target_include_directories(dovibaker_settings INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_definitions(dovibaker_settings INTERFACE
    $<$<PLATFORM_ID:Windows>:_CRT_SECURE_NO_WARNINGS>
    $<$<PLATFORM_ID:Windows>:NOMINMAX>
    $<$<PLATFORM_ID:Windows>:WIN32_LEAN_AND_MEAN>)

if(MSVC)
    target_compile_options(dovibaker_settings INTERFACE
        /utf-8            # the sources are plain UTF-8, without BOM
        /W3
        /MP               # parallel compilation
        /EHsc
        $<$<AND:$<NOT:$<CONFIG:Debug>>,$<BOOL:${DOVIBAKER_FAST_MATH}>>:/fp:fast>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(dovibaker_settings INTERFACE
        $<$<AND:$<NOT:$<CONFIG:Debug>>,$<BOOL:${DOVIBAKER_FAST_MATH}>>:-ffast-math>
        -Wno-unused-variable
        -Wno-unused-but-set-variable)
endif()

# Put every product into one directory, so the plugins, the tools and the
# runtime library of libdovi end up side by side.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
if(NOT MSVC)
    # Plugins (*.so) and the console tools belong together.
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

# Solution folders for the Visual Studio generator.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
