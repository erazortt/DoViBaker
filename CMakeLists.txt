cmake_minimum_required(VERSION 3.16)

project(DoViBaker LANGUAGES CXX)

add_library(libdovibaker SHARED
    DoViBaker/AvisynthEntry.cpp
    DoViBaker/cube.cpp
    DoViBaker/DoViBaker.cpp
    DoViBaker/DoViProcessor.cpp
    DoViBaker/lut.cpp
    DoViBaker/lut_avx2.cpp
    DoViBaker/lut_avx512.cpp
    DoViBaker/lut_sse41.cpp
    DoViBaker/lut_x86.cpp
)

target_include_directories(libdovibaker PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
if (build_type STREQUAL debug)
    target_compile_definitions(libdovibaker PRIVATE DEBUG_BUILD)
else (build_type STREQUAL release)
    target_compile_definitions(libdovibaker PRIVATE RELEASE_BUILD)
endif ()

message(STATUS "Build type - ${CMAKE_BUILD_TYPE}")

target_compile_features(libdovibaker PRIVATE cxx_std_17)

set_source_files_properties(DoViBaker/lut_sse41.cpp PROPERTIES COMPILE_OPTIONS "-mfpmath=sse;-msse4.1")
set_source_files_properties(DoViBaker/lut_avx2.cpp.cpp PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
set_source_files_properties(DoViBaker/lut_avx512.cpp PROPERTIES COMPILE_OPTIONS "-mavx512f;-mavx512bw;-mavx512dq;-mavx512vl;-mfma")

target_link_libraries(libdovibaker dovi)

find_package (Git)

if (GIT_FOUND)
    execute_process (COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
        OUTPUT_VARIABLE ver
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set_target_properties(libdovibaker PROPERTIES OUTPUT_NAME "dovibaker.${ver}")
else ()
    message (STATUS "GIT not found")
endif ()

include(GNUInstallDirs)

INSTALL(TARGETS libdovibaker LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/avisynth")

# uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
