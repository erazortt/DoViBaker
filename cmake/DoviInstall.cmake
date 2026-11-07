# Install rules for the products: the plugins go into the Avisynth plugin
# directory, the console tools into bin/.

function(dovibaker_install_plugin target)
    if(WIN32)
        install(TARGETS ${target} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
        install(TARGETS ${target} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/avisynth)
    endif()
endfunction()

function(dovibaker_install_tool target)
    install(TARGETS ${target} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endfunction()

# The scripts are helpers for Avisynth, not build products.
option(DOVIBAKER_INSTALL_SCRIPTS "Install the helper scripts from scripts/" OFF)
if(DOVIBAKER_INSTALL_SCRIPTS)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/
        DESTINATION ${CMAKE_INSTALL_DATADIR}/dovibaker/scripts
        FILES_MATCHING PATTERN "*.avsi" PATTERN "*.avs")
endif()

# "make uninstall" / "cmake --build . --target uninstall"
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
