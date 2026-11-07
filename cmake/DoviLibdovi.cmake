# Finds or builds libdovi, the C API of the dovi_tool submodule, and defines
# the imported target dovi::dovi.
#
# Resolution order:
#   1. DOVIBAKER_LIBDOVI_PREFIX (the --prefix of `cargo cinstall`, for instance)
#   2. pkg-config / find_library of an installed libdovi
#   3. cargo build of the dovi_tool submodule (DOVIBAKER_BUILD_LIBDOVI, default)
#
# The header is taken from the found installation; the copy vendored in
# include/dovi is the fallback.

set(DOVIBAKER_LIBDOVI_PREFIX "" CACHE PATH
    "Prefix of an already installed libdovi (contains include/ and lib/)")
option(DOVIBAKER_BUILD_LIBDOVI
    "Build libdovi from the dovi_tool submodule with cargo" ON)
set(DOVIBAKER_LIBDOVI_CARGO_ARGS "" CACHE STRING
    "Extra arguments passed to cargo when building libdovi")
set(DOVIBAKER_LIBDOVI_IMPLIB "" CACHE FILEPATH
    "Import library to link against when libdovi is built by cargo (empty = auto)")

set(dovi_LIBDOVI_VENDORED_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(dovi_LIBDOVI_SUBMODULE ${CMAKE_CURRENT_SOURCE_DIR}/dovi_tool/dolby_vision)

function(dovibaker_define_libdovi implib runtime header)
    add_library(dovi::dovi SHARED IMPORTED GLOBAL)
    set_target_properties(dovi::dovi PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${header}"
        IMPORTED_LOCATION "${runtime}"
        IMPORTED_IMPLIB "${implib}")

    # libdovi is Rust; its C API pulls in these system libraries on Windows.
    if(WIN32)
        set_property(TARGET dovi::dovi APPEND PROPERTY INTERFACE_LINK_LIBRARIES
            ws2_32 bcrypt userenv shlwapi ntdll)
    endif()
endfunction()

# On Windows find_library finds the import library, the runtime is elsewhere.
function(dovibaker_dovi_runtime out implib search_dir)
    set(${out} ${implib} PARENT_SCOPE)
    if(WIN32 AND implib MATCHES "\.lib$")
        find_file(dovi_DLL NAMES dovi.dll dolby_vision.dll
            PATHS ${search_dir} ${search_dir}/../bin NO_DEFAULT_PATH)
        if(dovi_DLL)
            set(${out} ${dovi_DLL} PARENT_SCOPE)
        endif()
    endif()
endfunction()

# Puts the libdovi runtime next to a product, so that the plugins find it at
# runtime without any installation step.
function(dovibaker_deploy_libdovi product)
    if(NOT TARGET dovi::dovi)
        return()
    endif()
    get_target_property(dovi_runtime dovi::dovi IMPORTED_LOCATION)
    # No EXISTS check here: when libdovi is built by cargo the file shows up
    # only during the build, the copy runs after it.
    if(NOT dovi_runtime OR dovi_runtime STREQUAL "dovi_runtime-NOTFOUND")
        return()
    endif()
    # Cargo has to run before the product is linked. Imported targets can carry
    # the dependency, but not every generator honours it, so add it as well.
    if(DOVI_CARGO_TARGET)
        add_dependencies(${product} ${DOVI_CARGO_TARGET})
    endif()

    if(dovi_runtime MATCHES "\.lib$")
        # only an import library is available, the runtime has to be installed
        # by the libdovi package itself
        return()
    endif()

    if(WIN32)
        add_custom_command(TARGET ${product} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${dovi_runtime} $<TARGET_FILE_DIR:${product}>
            COMMENT "Copying libdovi next to ${product}")
    else()
        # libdovi is installed as a plain library, reachable from both bin/
        # and the plugin directory.
        set_target_properties(${product} PROPERTIES
            INSTALL_RPATH "$ORIGIN:$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
        install(FILES ${dovi_runtime} DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
endfunction()

# ---------------------------------------------------------------------------
# 1. an installation chosen by the user
if(DOVIBAKER_LIBDOVI_PREFIX)
    find_path(dovi_INCLUDE_DIR
        NAMES dovi/rpu_parser.h libdovi/rpu_parser.h rpu_parser.h
        HINTS ${DOVIBAKER_LIBDOVI_PREFIX}/include)
    find_library(dovi_LIBRARY NAMES dovi dovi.dll libdovi libdovi.dll.a
        HINTS ${DOVIBAKER_LIBDOVI_PREFIX}/lib ${DOVIBAKER_LIBDOVI_PREFIX}/bin)
    if(NOT dovi_LIBRARY)
        message(FATAL_ERROR
            "No libdovi found below DOVIBAKER_LIBDOVI_PREFIX=${DOVIBAKER_LIBDOVI_PREFIX}")
    endif()
    if(NOT dovi_INCLUDE_DIR)
        set(dovi_INCLUDE_DIR ${dovi_LIBDOVI_VENDORED_INCLUDE})
    endif()
    dovibaker_dovi_runtime(dovi_RUNTIME ${dovi_LIBRARY} ${DOVIBAKER_LIBDOVI_PREFIX}/lib)
    dovibaker_define_libdovi(${dovi_LIBRARY} ${dovi_RUNTIME} ${dovi_INCLUDE_DIR})
    message(STATUS "libdovi: using ${dovi_LIBRARY}")
    return()
endif()

# 2. an installation of the system
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(DOVI_PC QUIET dovi)
endif()

find_library(dovi_LIBRARY NAMES dovi dovi.dll libdovi libdovi.dll.a
    HINTS ${DOVI_PC_LIBRARY_DIRS})
find_path(dovi_INCLUDE_DIR NAMES dovi/rpu_parser.h libdovi/rpu_parser.h
    HINTS ${DOVI_PC_INCLUDE_DIRS} ${dovi_LIBDOVI_VENDORED_INCLUDE})

if(dovi_LIBRARY)
    if(NOT dovi_INCLUDE_DIR)
        set(dovi_INCLUDE_DIR ${dovi_LIBDOVI_VENDORED_INCLUDE})
    endif()
    get_filename_component(dovi_LIB_DIR ${dovi_LIBRARY} DIRECTORY)
    dovibaker_dovi_runtime(dovi_RUNTIME ${dovi_LIBRARY} ${dovi_LIB_DIR})
    dovibaker_define_libdovi(${dovi_LIBRARY} ${dovi_RUNTIME} ${dovi_INCLUDE_DIR})
    message(STATUS "libdovi: using ${dovi_LIBRARY}")
    return()
endif()

# 3. build the submodule with cargo
if(NOT DOVIBAKER_BUILD_LIBDOVI)
    message(FATAL_ERROR
        "No libdovi found. Install it (see dovi_tool/dolby_vision/README.md,"
        " `cargo cinstall --release`) and set DOVIBAKER_LIBDOVI_PREFIX, or"
        " build it from the submodule with DOVIBAKER_BUILD_LIBDOVI=ON.")
endif()

if(NOT EXISTS ${dovi_LIBDOVI_SUBMODULE}/Cargo.toml)
    message(FATAL_ERROR
        "The dovi_tool submodule is missing. Fetch it with:\n"
        "    git submodule update --init --recursive")
endif()

find_program(CARGO_EXECUTABLE cargo)
find_program(RUSTC_EXECUTABLE rustc)
if(NOT CARGO_EXECUTABLE OR NOT RUSTC_EXECUTABLE)
    message(FATAL_ERROR
        "cargo/rustc not found - install Rust (https://rust-lang.org/tools/install/)"
        " to build libdovi from the submodule, or set DOVIBAKER_LIBDOVI_PREFIX"
        " to an existing installation.")
endif()

# The Rust triple decides where cargo puts its artifacts (and which target has
# to be installed: `rustup target add <triple>`). It is only needed to build
# libdovi here, cargo builds the host platform by default.
execute_process(COMMAND ${RUSTC_EXECUTABLE} -vV
    OUTPUT_VARIABLE rustc_version_output ERROR_QUIET)
string(REGEX MATCH "host: ([a-z0-9_-]+)" _dovi_rust_host "${rustc_version_output}")
set(dovi_RUST_HOST_TRIPLE "${CMAKE_MATCH_1}")
if(NOT dovi_RUST_HOST_TRIPLE)
    message(FATAL_ERROR "Could not determine the Rust host triple from 'rustc -vV'")
endif()

if(NOT DOVIBAKER_LIBDOVI_RUST_TRIPLE)
    if(MSVC)
        set(dovi_GUESS_TRIPLE "pc-windows-msvc")
    elseif(MINGW)
        set(dovi_GUESS_TRIPLE "pc-windows-gnu")
    elseif(APPLE)
        set(dovi_GUESS_TRIPLE "apple-darwin")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(dovi_GUESS_TRIPLE "unknown-linux-gnu")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Android")
        set(dovi_GUESS_TRIPLE "linux-android")
    endif()

    if(dovi_GUESS_TRIPLE)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
            set(dovi_GUESS_TRIPLE "x86_64-${dovi_GUESS_TRIPLE}")
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(i[3-6]86|x86)$")
            set(dovi_GUESS_TRIPLE "i686-${dovi_GUESS_TRIPLE}")
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
            set(dovi_GUESS_TRIPLE "aarch64-${dovi_GUESS_TRIPLE}")
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm$")
            set(dovi_GUESS_TRIPLE "armv7-${dovi_GUESS_TRIPLE}")
        else()
            set(dovi_GUESS_TRIPLE "")
        endif()
    endif()

    if(dovi_GUESS_TRIPLE)
        set(DOVIBAKER_LIBDOVI_RUST_TRIPLE "${dovi_GUESS_TRIPLE}" CACHE STRING
            "Rust target triple used to build libdovi from the submodule" FORCE)
    else()
        set(DOVIBAKER_LIBDOVI_RUST_TRIPLE "${dovi_RUST_HOST_TRIPLE}" CACHE STRING
            "Rust target triple used to build libdovi from the submodule" FORCE)
    endif()
endif()

if(NOT DOVIBAKER_LIBDOVI_RUST_TRIPLE STREQUAL dovi_RUST_HOST_TRIPLE)
    message(STATUS "libdovi: cargo target '${DOVIBAKER_LIBDOVI_RUST_TRIPLE}'"
        " (host is '${dovi_RUST_HOST_TRIPLE}')"
        " - install it with: rustup target add ${DOVIBAKER_LIBDOVI_RUST_TRIPLE}")
endif()

set(dovi_CARGO_TARGET_DIR ${CMAKE_BINARY_DIR}/_libdovi)
set(dovi_LIBDOVI_DIR ${dovi_CARGO_TARGET_DIR}/${DOVIBAKER_LIBDOVI_RUST_TRIPLE}/release)

# The library target of the dolby_vision package is named after the package,
# cargo-c is the only thing which renames it to "dovi".
set(dovi_LIBDOVI_CRATE_NAME dolby_vision)

set(dovi_LIBDOVI_SONAME_FLAG "")
if(WIN32)
    set(dovi_RUNTIME ${dovi_LIBDOVI_DIR}/${dovi_LIBDOVI_CRATE_NAME}.dll)
    if(MSVC)
        set(dovi_IMPLIB ${dovi_LIBDOVI_DIR}/${dovi_LIBDOVI_CRATE_NAME}.dll.lib)
    else()
        set(dovi_IMPLIB ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.dll.a)
    endif()
elseif(APPLE)
    set(dovi_RUNTIME ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.dylib)
    set(dovi_IMPLIB ${dovi_RUNTIME})
else()
    set(dovi_RUNTIME ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.so)
    set(dovi_IMPLIB ${dovi_RUNTIME})
    # Without a SONAME the plugins would record the absolute build path.
    set(dovi_LIBDOVI_SONAME_FLAG "-Clink-arg=-Wl,-soname,lib${dovi_LIBDOVI_CRATE_NAME}.so")
endif()

if(DOVIBAKER_LIBDOVI_IMPLIB)
    set(dovi_IMPLIB ${DOVIBAKER_LIBDOVI_IMPLIB})
endif()

separate_arguments(dovi_CARGO_EXTRA_ARGS NATIVE_COMMAND "${DOVIBAKER_LIBDOVI_CARGO_ARGS}")

# On ELF platforms the import "library" is the shared object itself.
set(dovi_CARGO_BYPRODUCTS ${dovi_RUNTIME})
if(NOT dovi_IMPLIB STREQUAL dovi_RUNTIME)
    list(APPEND dovi_CARGO_BYPRODUCTS ${dovi_IMPLIB})
endif()

add_custom_target(dovi_cargo
    COMMAND ${CMAKE_COMMAND} -E env RUSTFLAGS=${dovi_LIBDOVI_SONAME_FLAG}
        ${CARGO_EXECUTABLE}
            rustc
            --manifest-path ${dovi_LIBDOVI_SUBMODULE}/Cargo.toml
            --target-dir ${dovi_CARGO_TARGET_DIR}
            --target ${DOVIBAKER_LIBDOVI_RUST_TRIPLE}
            --release --locked --features capi
            --crate-type cdylib
            ${dovi_CARGO_EXTRA_ARGS}
    BYPRODUCTS ${dovi_CARGO_BYPRODUCTS}
    COMMENT "Building libdovi with cargo (${DOVIBAKER_LIBDOVI_RUST_TRIPLE}, release)"
    VERBATIM)
set_target_properties(dovi_cargo PROPERTIES FOLDER thirdparty)
set(DOVI_CARGO_TARGET dovi_cargo CACHE INTERNAL "cargo target which builds libdovi")

dovibaker_define_libdovi(${dovi_IMPLIB} ${dovi_RUNTIME} ${dovi_LIBDOVI_VENDORED_INCLUDE})
add_dependencies(dovi::dovi dovi_cargo)
message(STATUS "libdovi: building ${dovi_RUNTIME} from the submodule (cargo, release)")
