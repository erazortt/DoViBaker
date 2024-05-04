cmake_minimum_required(VERSION 3.28)

if (NOT CMAKE_GENERATOR MATCHES "Visual Studio")
    if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
    endif()
endif()

project(DoViBaker LANGUAGES CXX)

add_library(DoViBaker SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/AvisynthEntry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/DoViBaker.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/DoViCubes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/DoViEetf.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/DoViProcessor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/DoViStatsFileLoader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DoViBaker/DoViTonemap.cpp
)

target_include_directories(DoViBaker PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (UNIX)
    target_include_directories(DoViBaker PRIVATE
        /usr/local/include/avisynth
        /usr/local/include
    )
endif()

if (NOT CMAKE_GENERATOR MATCHES "Visual Studio")
    string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
    if (build_type STREQUAL debug)
        target_compile_definitions(DoViBaker PRIVATE DEBUG_BUILD)
    else (build_type STREQUAL release)
        target_compile_definitions(DoViBaker PRIVATE RELEASE_BUILD)
    endif()

    message(STATUS "Build type - ${CMAKE_BUILD_TYPE}")
endif()

target_compile_features(DoViBaker PRIVATE cxx_std_20)

# graphengine
add_library(graphengine STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/graphengine/cpuinfo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/graphengine/graph.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/graphengine/node.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/graphengine/x86/cpuinfo_x86.cpp
)

target_include_directories(graphengine PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/graphengine
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/graphengine/x86
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/include
)

if (UNIX)
    target_compile_options(graphengine PRIVATE "-fPIC")
endif()
# end graphengine

# timecube
add_library(timecube STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/cube.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/lut.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/timecube.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_x86.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_sse41.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_avx2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_avx512.cpp
)

target_include_directories(timecube PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86
    ${CMAKE_CURRENT_SOURCE_DIR}/timecube/graphengine/include
)

if (MSVC)
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_avx2.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_avx512.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX512")
else()
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_sse41.cpp PROPERTIES COMPILE_OPTIONS "-mfpmath=sse;-msse4.1")
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_avx2.cpp PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/timecube/timecube/x86/lut_avx512.cpp PROPERTIES COMPILE_OPTIONS "-mavx512f;-mavx512bw;-mavx512dq;-mavx512vl;-mfma")
    
    target_compile_options(timecube PRIVATE "-fPIC")
endif()

target_link_libraries(timecube PRIVATE graphengine)
#end timecube

find_library(dovi dovi)

target_link_libraries(DoViBaker PRIVATE
    ${dovi}
    timecube
)

if (WIN32)
    target_link_libraries(DoViBaker PRIVATE
    ws2_32
    bcrypt
    userenv
    shlwapi
    ntdll
    )
endif()

if (UNIX)
    find_package (Git)

    if (GIT_FOUND)
        execute_process (COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
            OUTPUT_VARIABLE ver
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set_target_properties(DoViBaker PROPERTIES OUTPUT_NAME "dovibaker.${ver}")
    else ()
        message (STATUS "GIT not found")
    endif ()

    include(GNUInstallDirs)

    INSTALL(TARGETS DoViBaker LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/avisynth")

    # uninstall target
    if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    endif()
endif()
