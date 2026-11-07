# Builds the timecube LUT library (and its graphengine helper) from the
# "timecube" submodule as static libraries which are linked into the plugins.
#
# Defines the static targets `timecube` and `graphengine`.

set(dovi_timecube_root ${CMAKE_CURRENT_SOURCE_DIR}/timecube)

if(NOT EXISTS ${dovi_timecube_root}/timecube/timecube.cpp)
    message(FATAL_ERROR
        "The timecube submodule is missing. Fetch it with:\n"
        "    git submodule update --init --recursive")
endif()

set(dovi_timecube_sources
    ${dovi_timecube_root}/timecube/cube.cpp
    ${dovi_timecube_root}/timecube/lut.cpp
    ${dovi_timecube_root}/timecube/timecube.cpp)

set(dovi_graphengine_sources
    ${dovi_timecube_root}/graphengine/graphengine/cpuinfo.cpp
    ${dovi_timecube_root}/graphengine/graphengine/graph.cpp
    ${dovi_timecube_root}/graphengine/graphengine/node.cpp)

# timecube only has x86 assembly paths, the plain C++ fallback is used elsewhere.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|i[3-6]86)$")
    list(APPEND dovi_timecube_sources
        ${dovi_timecube_root}/timecube/x86/lut_x86.cpp
        ${dovi_timecube_root}/timecube/x86/lut_sse41.cpp
        ${dovi_timecube_root}/timecube/x86/lut_avx2.cpp
        ${dovi_timecube_root}/timecube/x86/lut_avx512.cpp)
    list(APPEND dovi_graphengine_sources
        ${dovi_timecube_root}/graphengine/graphengine/x86/cpuinfo_x86.cpp)
    set(dovi_timecube_x86 TRUE)
else()
    set(dovi_timecube_x86 FALSE)
    message(STATUS "timecube: no x86 optimizations for '${CMAKE_SYSTEM_PROCESSOR}'")
endif()

add_library(graphengine STATIC EXCLUDE_FROM_ALL ${dovi_graphengine_sources})
add_library(timecube STATIC EXCLUDE_FROM_ALL ${dovi_timecube_sources})
add_library(dovibaker::timecube ALIAS timecube)
add_library(dovibaker::graphengine ALIAS graphengine)

foreach(dovi_tc_target graphengine timecube)
    target_include_directories(${dovi_tc_target} PRIVATE
        ${dovi_timecube_root}/graphengine/graphengine
        ${dovi_timecube_root}/graphengine/graphengine/x86
        ${dovi_timecube_root}/graphengine/include
        ${dovi_timecube_root}/timecube
        ${dovi_timecube_root}/timecube/x86)
    target_compile_features(${dovi_tc_target} PRIVATE cxx_std_17)
    # graphengine is embedded, so its symbols get their own namespace.
    target_compile_definitions(${dovi_tc_target} PRIVATE GRAPHENGINE_IMPL_NAMESPACE=timecube)
    if(NOT MSVC)
        set_target_properties(${dovi_tc_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
    set_target_properties(${dovi_tc_target} PROPERTIES FOLDER thirdparty)
endforeach()

target_include_directories(timecube PRIVATE ${dovi_timecube_root}/graphengine/include)
# timecube exposes a C API only, graphengine stays an implementation detail.
target_link_libraries(timecube PRIVATE graphengine)

if(MINGW)
    # timecube.cpp includes <Windows.h> with exactly that spelling, which a case
    # sensitive MinGW sysroot does not provide.
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("#include <Windows.h>\nint main() { return 0; }"
        DOVI_MINGW_WINDOWS_H_FOUND)
    if(NOT DOVI_MINGW_WINDOWS_H_FOUND)
        set(dovi_mingw_case_fixups ${CMAKE_BINARY_DIR}/_mingw_include_fixups)
        file(WRITE ${dovi_mingw_case_fixups}/Windows.h "#include_next <windows.h>\n")
        target_include_directories(timecube PRIVATE ${dovi_mingw_case_fixups})
        message(STATUS "timecube: added a spelling fixup header for <Windows.h>")
    endif()
endif()

if(dovi_timecube_x86)
    # CUBE_X86 switches on the SIMD paths of timecube (see timecube/Makefile).
    target_compile_definitions(timecube PRIVATE CUBE_X86)

    if(MSVC)
        # /bigobj: the SIMD lookup tables exceed the default object file format.
        set_source_files_properties(${dovi_timecube_root}/timecube/x86/lut_sse41.cpp
            PROPERTIES COMPILE_OPTIONS "/EHsc;/bigobj")
        set_source_files_properties(${dovi_timecube_root}/timecube/x86/lut_avx2.cpp
            PROPERTIES COMPILE_OPTIONS "/EHsc;/bigobj;/arch:AVX2")
        set_source_files_properties(${dovi_timecube_root}/timecube/x86/lut_avx512.cpp
            PROPERTIES COMPILE_OPTIONS "/EHsc;/bigobj;/arch:AVX512")
    else()
        set_source_files_properties(${dovi_timecube_root}/timecube/x86/lut_sse41.cpp
            PROPERTIES COMPILE_OPTIONS "-msse4.1")
        set_source_files_properties(${dovi_timecube_root}/timecube/x86/lut_avx2.cpp
            PROPERTIES COMPILE_OPTIONS "-mf16c;-mavx2;-mfma")
        set_source_files_properties(${dovi_timecube_root}/timecube/x86/lut_avx512.cpp
            PROPERTIES COMPILE_OPTIONS "-mf16c;-mfma;-mavx512f;-mavx512cd;-mavx512vl;-mavx512bw;-mavx512dq")
    endif()
endif()
