cmake_minimum_required(VERSION 3.18)

# Multi-config generators (Visual Studio, Xcode) have no single build type.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

project(DoViBaker LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(GNUInstallDirs)
include(CTest)            # -> BUILD_TESTING, the smoke tests of the products
include(DoviVersion)      # -> DOVIBAKER_VERSION
include(DoviSettings)     # -> dovibaker::settings interface target

# ---------------------------------------------------------------------------
# What to build. All four products are built in one go by default.
# ---------------------------------------------------------------------------
option(DOVIBAKER_BUILD_ALL "Build all products (plugins and tools)" ON)
option(DOVIBAKER_BUILD_BAKER "Build the DoViBaker Avisynth plugin" ${DOVIBAKER_BUILD_ALL})
option(DOVIBAKER_BUILD_ANALYZER "Build the DoviAnalyzer console tool" ${DOVIBAKER_BUILD_ALL})
option(DOVIBAKER_BUILD_LUTGEN "Build the DoViLutGen console tool" ${DOVIBAKER_BUILD_ALL})
option(DOVIBAKER_BUILD_AVSCUBE "Build the AVSCube Avisynth plugin" ${DOVIBAKER_BUILD_ALL})

if(NOT DOVIBAKER_BUILD_BAKER AND NOT DOVIBAKER_BUILD_ANALYZER
   AND NOT DOVIBAKER_BUILD_LUTGEN AND NOT DOVIBAKER_BUILD_AVSCUBE)
    message(FATAL_ERROR "Nothing to build - enable at least one DOVIBAKER_BUILD_* option")
endif()

foreach(product IN ITEMS BAKER ANALYZER LUTGEN AVSCUBE)
    if(DOVIBAKER_BUILD_${product})
        message(STATUS "DoViBaker: will build ${product}")
    endif()
endforeach()

# ---------------------------------------------------------------------------
# The dependencies, only the ones the selected products need.
# ---------------------------------------------------------------------------
if(DOVIBAKER_BUILD_BAKER OR DOVIBAKER_BUILD_AVSCUBE)
    include(DoviTimecube)     # -> timecube, graphengine
endif()

if(DOVIBAKER_BUILD_BAKER OR DOVIBAKER_BUILD_ANALYZER)
    include(DoviLibdovi)      # -> dovi::dovi
endif()

# install rules (defines dovibaker_install_plugin/tool used below)
include(DoviInstall)

add_subdirectory(DoViBaker)

dovibaker_version_report()
