# The four products of this repository. Every product tells the sources which
# one it is by defining DOVI_TARGET (see include/def.h), therefore all of them
# can be built from a single build tree - no editing of files needed.

if(WIN32)
    # Historic naming of the plugin files: DoViBaker_x64.dll / DoViBaker_x86.dll
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(dovi_arch_suffix "_x64")
    else()
        set(dovi_arch_suffix "_x86")
    endif()
    set(dovi_baker_name "DoViBaker${dovi_arch_suffix}")
    set(dovi_avscube_name "AVSCube${dovi_arch_suffix}")
else()
    # Naming of the historic POSIX install: dovibaker.<version>.so
    set(dovi_baker_name "dovibaker.${DOVIBAKER_VERSION}")
    set(dovi_avscube_name "avscube.${DOVIBAKER_VERSION}")
endif()

# ---------------------------------------------------------------------------
# DoViBaker: the Avisynth plugin which bakes BL + EL + RPU into a PQ stream.
# It registers DoViBaker, DoViTonemap, DoViCubes and DoViStatsFileLoader.
# ---------------------------------------------------------------------------
if(DOVIBAKER_BUILD_BAKER)
    add_library(DoViBaker SHARED
        AvisynthEntry.cpp
        DoViBaker.cpp
        DoViCubes.cpp
        DoViEetf.cpp
        DoViProcessor.cpp
        DoViStatsFileLoader.cpp
        DoViTonemap.cpp)

    target_compile_definitions(DoViBaker PRIVATE DOVI_TARGET=DOVI_PRODUCT_BAKER)
    target_link_libraries(DoViBaker PRIVATE
        dovibaker::settings
        dovibaker::timecube
        dovi::dovi)

    set_target_properties(DoViBaker PROPERTIES
        OUTPUT_NAME ${dovi_baker_name}
        PREFIX ""          # plugins carry the name of the product, no "lib"
        CXX_EXTENSIONS OFF
        FOLDER plugins)

    dovibaker_deploy_libdovi(DoViBaker)
    dovibaker_install_plugin(DoViBaker)
endif()

# ---------------------------------------------------------------------------
# DoViAnalyzer: console tool, analyzes an RPU.bin file
# ---------------------------------------------------------------------------
if(DOVIBAKER_BUILD_ANALYZER)
    add_executable(DoViAnalyzer
        AvisynthEntry.cpp
        DoViEetf.cpp
        DoViProcessor.cpp)

    target_compile_definitions(DoViAnalyzer PRIVATE DOVI_TARGET=DOVI_PRODUCT_ANALYZER)
    target_link_libraries(DoViAnalyzer PRIVATE
        dovibaker::settings
        dovi::dovi)

    set_target_properties(DoViAnalyzer PROPERTIES FOLDER tools)

    dovibaker_deploy_libdovi(DoViAnalyzer)
    dovibaker_install_tool(DoViAnalyzer)
endif()

# ---------------------------------------------------------------------------
# DoViLutGen: console tool, creates the PQ/HLG/SDR LUTs. Stands on its own.
# ---------------------------------------------------------------------------
if(DOVIBAKER_BUILD_LUTGEN)
    add_executable(DoViLutGen DoViLutGen.cpp)

    target_compile_definitions(DoViLutGen PRIVATE DOVI_TARGET=DOVI_PRODUCT_LUTGEN)
    target_link_libraries(DoViLutGen PRIVATE dovibaker::settings)

    set_target_properties(DoViLutGen PROPERTIES FOLDER tools)

    dovibaker_install_tool(DoViLutGen)
endif()

# ---------------------------------------------------------------------------
# AVSCube: Avisynth plugin which applies one LUT file
# ---------------------------------------------------------------------------
if(DOVIBAKER_BUILD_AVSCUBE)
    add_library(AVSCube SHARED AvsCube.cpp)

    target_compile_definitions(AVSCube PRIVATE DOVI_TARGET=DOVI_PRODUCT_AVSCUBE)
    target_link_libraries(AVSCube PRIVATE
        dovibaker::settings
        dovibaker::timecube)

    set_target_properties(AVSCube PROPERTIES
        OUTPUT_NAME ${dovi_avscube_name}
        PREFIX ""          # plugins carry the name of the product, no "lib"
        CXX_EXTENSIONS OFF
        FOLDER plugins)

    dovibaker_install_plugin(AVSCube)
endif()

# ---------------------------------------------------------------------------
# Smoke tests: `ctest --test-dir build` after the build.
# ---------------------------------------------------------------------------
if(BUILD_TESTING)
    if(TARGET DoViLutGen)
        add_test(NAME DoViLutGen.pq2hlg
            COMMAND DoViLutGen ${CMAKE_CURRENT_BINARY_DIR}/test_pq2hlg.cube
                    -s 5 -i 0 -o 0)
        # every BT.709 SDR output format: 2 and 4 convert scene referred, 3 and 5
        # display referred, 2 and 3 map colors, 4 and 5 clip them
        foreach(dovi_sdr_format IN ITEMS 2 3 4 5)
            add_test(NAME DoViLutGen.pq2sdr709_format${dovi_sdr_format}
                COMMAND DoViLutGen ${CMAKE_CURRENT_BINARY_DIR}/test_pq2sdr709_o${dovi_sdr_format}.cube
                        -s 5 -i 0 -o ${dovi_sdr_format})
        endforeach()
        add_test(NAME DoViLutGen.hlg2sdr709
            COMMAND DoViLutGen ${CMAKE_CURRENT_BINARY_DIR}/test_hlg2sdr709.cube
                    -s 5 -i 2 -o 2)
    endif()

    set(dovi_test_rpu ${CMAKE_SOURCE_DIR}/dovi_tool/assets/hevc_tests/regular_rpu.bin)
    if(TARGET DoViAnalyzer AND EXISTS ${dovi_test_rpu})
        add_test(NAME DoViAnalyzer.rpu
            COMMAND DoViAnalyzer ${dovi_test_rpu})
    endif()
endif()

# The sources are shared by the products, keep them grouped in the solution.
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
    AvisynthEntry.cpp AvsCube.cpp DoViBaker.cpp DoViCubes.cpp DoViEetf.cpp
    DoViLutGen.cpp DoViProcessor.cpp DoViStatsFileLoader.cpp DoViTonemap.cpp)
