# Determines DOVIBAKER_VERSION from the git tags (falls back to the version
# given in project() when the tree is not a git checkout).

set(DOVIBAKER_VERSION "0.0.0")

find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE DOVIBAKER_VERSION_TAG
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE DOVIBAKER_GIT_RESULT)
    if(DOVIBAKER_GIT_RESULT EQUAL 0 AND DOVIBAKER_VERSION_TAG)
        string(REGEX REPLACE "^v" "" DOVIBAKER_VERSION "${DOVIBAKER_VERSION_TAG}")
    endif()
endif()

set(DOVIBAKER_VERSION "${DOVIBAKER_VERSION}" CACHE STRING
    "Version used for the plugin file name" FORCE)

function(dovibaker_version_report)
    message(STATUS "DoViBaker ${DOVIBAKER_VERSION}: build type '${CMAKE_BUILD_TYPE}', generator '${CMAKE_GENERATOR}'")
endfunction()
