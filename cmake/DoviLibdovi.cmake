# Finds or builds libdovi, the C API of the dovi_tool submodule, and defines
# the imported target dovi::dovi.
#
# Link mode (DOVIBAKER_LIBDOVI_LINK, declared in the top level CMakeLists.txt):
#   static   cargo builds a static library and its objects end up inside every
#            product. No runtime library is built, copied, installed or looked
#            up, so a plugin can never fail to load over a missing
#            dolby_vision.dll. This is the default.
#   shared   cargo builds a library (dolby_vision.dll / libdolby_vision.so)
#            which has to travel with the plugins, see dovibaker_deploy_libdovi.
#
# Resolution order, the same for both modes:
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
    "Import library to link against in shared mode when libdovi is built by cargo (empty = auto)")

set(dovi_LIBDOVI_VENDORED_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(dovi_LIBDOVI_SUBMODULE ${CMAKE_CURRENT_SOURCE_DIR}/dovi_tool/dolby_vision)

# cargo-c renames the library to "dovi", a plain cargo build keeps the name of
# the crate ("dolby_vision"). Names with a dot are searched as file names, which
# is how the static archives are told apart from the import libraries.
set(dovi_LIBDOVI_CRATE_NAME dolby_vision)
set(dovi_LIBDOVI_STATIC_NAMES dovi.lib libdovi.a dolby_vision.lib libdolby_vision.a)
set(dovi_LIBDOVI_SHARED_NAMES dovi dovi.dll libdovi libdovi.dll.a dovi.dll.lib)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# An MSVC ".lib" is either a real static library or the import library of a DLL,
# and only the first one can be linked statically. The file names tell them
# apart, and that is what the searches below rely on: cargo-c installs "dovi.lib"
# (static) next to "dovi.dll.lib" (importing the DLL). This catches the extra
# case of the Windows SDK style import libraries, which start with the marker
# 00 00 FF FF. Archives, the COFF import libraries of link.exe included, start
# with "!<arch>" and are taken at face value.
function(dovibaker_is_pe_import_library out path)
    set(${out} FALSE PARENT_SCOPE)
    if(NOT path MATCHES "\.lib$")
        return()
    endif()
    file(READ "${path}" dovi_marker LIMIT 4 HEX)
    if(dovi_marker STREQUAL "0000ffff")
        set(${out} TRUE PARENT_SCOPE)
    endif()
endfunction()

# A Rust static library does not carry the system libraries its std uses, the
# final link has to name them, and which ones they are differs per target. Ask
# rustc: building a static library prints a "native-static-libs" note.
#   dovibaker_rust_native_libs(<out> [<rust triple>])
function(dovibaker_rust_native_libs out)
    set(${out} "" PARENT_SCOPE)
    if(NOT RUSTC_EXECUTABLE)
        find_program(RUSTC_EXECUTABLE rustc QUIET)
    endif()
    if(NOT RUSTC_EXECUTABLE)
        return()
    endif()

    set(dovi_probe_args "")
    if(ARGN)
        list(APPEND dovi_probe_args --target ${ARGN})
    endif()

    set(dovi_probe_dir ${CMAKE_BINARY_DIR}/_libdovi_probe)
    file(WRITE ${dovi_probe_dir}/lib.rs "// asks rustc which system libraries a static library needs\n")
    execute_process(COMMAND ${RUSTC_EXECUTABLE}
        --crate-type staticlib ${dovi_probe_args} --print native-static-libs
        ${dovi_probe_dir}/lib.rs -o ${dovi_probe_dir}/probe.a
        ERROR_VARIABLE dovi_probe_out
        RESULT_VARIABLE dovi_probe_result
        OUTPUT_QUIET ERROR_STRIP_TRAILING_WHITESPACE)
    if(NOT dovi_probe_result EQUAL 0)
        message(STATUS "libdovi: 'rustc --print native-static-libs' failed, the"
            " built in list of system libraries is used")
        return()
    endif()

    string(REGEX MATCH "native-static-libs:[^\r\n]*" dovi_probe_line "${dovi_probe_out}")
    string(REGEX REPLACE "native-static-libs:" "" dovi_probe_list "${dovi_probe_line}")
    string(REPLACE "\"" "" dovi_probe_list "${dovi_probe_list}")
    # The note is separated by spaces and a CMake list is not: without this the
    # whole line is one library name, which the GNU driver happens to split and
    # the MSVC linker reports as LNK1181 for a single, long file name.
    string(REGEX REPLACE "[ \t]+" ";" dovi_probe_list "${dovi_probe_list}")

    set(dovi_probe_libs "")
    foreach(dovi_lib IN LISTS dovi_probe_list)
        string(STRIP "${dovi_lib}" dovi_lib)
        # The MSVC targets list raw linker switches too ("/defaultlib:msvcrt.lib").
        # They are not libraries, and the CRT is picked in DoviSettings.cmake.
        if(dovi_lib MATCHES "^[-/\\\\]" AND NOT dovi_lib MATCHES "^-l")
            continue()
        endif()
        string(REGEX REPLACE "^-l" "" dovi_lib "${dovi_lib}")
        string(REGEX REPLACE "^lib" "" dovi_lib "${dovi_lib}")
        string(REGEX REPLACE "\.(lib|a|so)$" "" dovi_lib "${dovi_lib}")
        # plain names only, so that nothing reaches the link line quoted
        if(dovi_lib MATCHES "^[A-Za-z0-9_+-]+$")
            list(APPEND dovi_probe_libs ${dovi_lib})
        endif()
    endforeach()

    if(dovi_probe_libs)
        list(REMOVE_DUPLICATES dovi_probe_libs)
        set(${out} "${dovi_probe_libs}" PARENT_SCOPE)
    endif()
endfunction()

# The dovi_* entry points of a Rust archive have default visibility and
# -fvisibility=hidden does not reach them. Left alone, every plugin which
# embeds libdovi would export them, and two of them in one Avisynth+ process
# would interpose each other's copies. PE exports only what is
# __declspec(dllexport)'ed, so this concerns ELF and Mach-O only.
function(dovibaker_keep_archive_symbols_private target archive)
    if(WIN32)
        return()
    endif()
    include(CheckLinkerFlag)
    get_filename_component(dovi_archive_name ${archive} NAME)

    check_linker_flag(CXX "-Wl,--exclude-libs,ALL" DOVI_LINKER_EXCLUDE_LIBS)
    if(DOVI_LINKER_EXCLUDE_LIBS)
        set_property(TARGET ${target} APPEND PROPERTY INTERFACE_LINK_OPTIONS
            "-Wl,--exclude-libs,ALL")
        message(STATUS "libdovi: keeping the symbols of ${dovi_archive_name}"
            " inside the products")
        return()
    endif()

    check_linker_flag(CXX "-Wl,-hidden-l${dovi_archive_name}" DOVI_LINKER_HIDDEN_L)
    if(DOVI_LINKER_HIDDEN_L)
        set_property(TARGET ${target} APPEND PROPERTY INTERFACE_LINK_OPTIONS
            "-Wl,-hidden-l${dovi_archive_name}")
        message(STATUS "libdovi: keeping the symbols of ${dovi_archive_name}"
            " inside the products")
        return()
    endif()

    message(WARNING "libdovi: this linker cannot hide the symbols of a static"
        " archive, the products will export the dovi_* entry points")
endfunction()

# dovibaker_define_libdovi(<STATIC|SHARED> <implib|archive> <runtime> <header>)
function(dovibaker_define_libdovi type implib runtime header)
    if(type STREQUAL "STATIC")
        add_library(dovi::dovi STATIC IMPORTED GLOBAL)
        set_target_properties(dovi::dovi PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${header}"
            IMPORTED_LOCATION "${implib}")
    else()
        add_library(dovi::dovi SHARED IMPORTED GLOBAL)
        set_target_properties(dovi::dovi PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${header}"
            IMPORTED_LOCATION "${runtime}"
            IMPORTED_IMPLIB "${implib}")
    endif()

    # libdovi is Rust; its C API pulls in these system libraries on Windows.
    if(WIN32)
        set_property(TARGET dovi::dovi APPEND PROPERTY INTERFACE_LINK_LIBRARIES
            ws2_32 bcrypt userenv shlwapi ntdll)
    endif()

    # The static library carries Rust's std, which names even more of them.
    if(type STREQUAL "STATIC" AND dovi_NATIVE_LIBS)
        # A name with a space in it is a parsed note which was not split: the
        # linker would report it as one long file name (LNK1181 on MSVC).
        foreach(dovi_link_lib IN LISTS dovi_NATIVE_LIBS)
            if(dovi_link_lib MATCHES "[ \t]")
                message(FATAL_ERROR
                    "libdovi: the system libraries of the Rust static library were"
                    " not split into a list ('${dovi_link_lib}'). Full list:"
                    " ${dovi_NATIVE_LIBS}")
            endif()
        endforeach()
        set_property(TARGET dovi::dovi APPEND PROPERTY INTERFACE_LINK_LIBRARIES
            ${dovi_NATIVE_LIBS})
        get_target_property(dovi_link_libs dovi::dovi INTERFACE_LINK_LIBRARIES)
        list(REMOVE_DUPLICATES dovi_link_libs)
        set_target_properties(dovi::dovi PROPERTIES
            INTERFACE_LINK_LIBRARIES "${dovi_link_libs}")
    endif()

    if(type STREQUAL "STATIC")
        dovibaker_keep_archive_symbols_private(dovi::dovi ${implib})
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
# runtime without any installation step. In static mode there is nothing to
# ship, only the cargo build has to have run before the product is linked.
function(dovibaker_deploy_libdovi product)
    if(NOT TARGET dovi::dovi)
        return()
    endif()
    # Cargo has to run before the product is linked. Imported targets can carry
    # the dependency, but not every generator honours it, so add it as well.
    if(DOVI_CARGO_TARGET)
        add_dependencies(${product} ${DOVI_CARGO_TARGET})
    endif()
    if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
        return()
    endif()

    get_target_property(dovi_runtime dovi::dovi IMPORTED_LOCATION)
    # No EXISTS check here: when libdovi is built by cargo the file shows up
    # only during the build, the copy runs after it.
    if(NOT dovi_runtime OR dovi_runtime STREQUAL "dovi_runtime-NOTFOUND")
        return()
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

# The system libraries a static libdovi needs, for the modes which do not ask
# rustc about them.
function(dovibaker_native_libs_fallback)
    if(NOT dovi_NATIVE_LIBS AND WIN32)
        set(dovi_NATIVE_LIBS ws2_32 bcrypt userenv shlwapi ntdll PARENT_SCOPE)
    endif()
endfunction()

# ---------------------------------------------------------------------------
# 1. an installation chosen by the user
if(DOVIBAKER_LIBDOVI_PREFIX)
    find_path(dovi_INCLUDE_DIR
        NAMES dovi/rpu_parser.h libdovi/rpu_parser.h rpu_parser.h
        HINTS ${DOVIBAKER_LIBDOVI_PREFIX}/include)

    if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
        find_library(dovi_LIBRARY NAMES ${dovi_LIBDOVI_STATIC_NAMES}
            HINTS ${DOVIBAKER_LIBDOVI_PREFIX}/lib ${DOVIBAKER_LIBDOVI_PREFIX}/lib64
                  ${DOVIBAKER_LIBDOVI_PREFIX}/bin)
        if(dovi_LIBRARY AND WIN32)
            dovibaker_is_pe_import_library(dovi_IS_IMPORT_LIB ${dovi_LIBRARY})
            if(dovi_IS_IMPORT_LIB)
                message(FATAL_ERROR
                    "${dovi_LIBRARY} is the import library of the libdovi DLL, not a"
                    " static library. `cargo cinstall` installs both; configure with"
                    " -DDOVIBAKER_LIBDOVI_LINK=shared to use the DLL.")
            endif()
        endif()
    else()
        find_library(dovi_LIBRARY NAMES ${dovi_LIBDOVI_SHARED_NAMES}
            HINTS ${DOVIBAKER_LIBDOVI_PREFIX}/lib ${DOVIBAKER_LIBDOVI_PREFIX}/bin)
    endif()

    if(NOT dovi_LIBRARY)
        message(FATAL_ERROR
            "No ${DOVIBAKER_LIBDOVI_LINK} libdovi found below"
            " DOVIBAKER_LIBDOVI_PREFIX=${DOVIBAKER_LIBDOVI_PREFIX}")
    endif()
    if(NOT dovi_INCLUDE_DIR)
        set(dovi_INCLUDE_DIR ${dovi_LIBDOVI_VENDORED_INCLUDE})
    endif()

    if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
        dovibaker_rust_native_libs(dovi_NATIVE_LIBS "${DOVIBAKER_LIBDOVI_RUST_TRIPLE}")
        dovibaker_native_libs_fallback()
        dovibaker_define_libdovi(STATIC ${dovi_LIBRARY} ${dovi_LIBRARY} ${dovi_INCLUDE_DIR})
    else()
        dovibaker_dovi_runtime(dovi_RUNTIME ${dovi_LIBRARY} ${DOVIBAKER_LIBDOVI_PREFIX}/lib)
        dovibaker_define_libdovi(SHARED ${dovi_LIBRARY} ${dovi_RUNTIME} ${dovi_INCLUDE_DIR})
    endif()
    message(STATUS "libdovi: using ${dovi_LIBRARY} (${DOVIBAKER_LIBDOVI_LINK})")
    return()
endif()

# 2. an installation of the system
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(DOVI_PC QUIET dovi)
endif()

if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
    find_library(dovi_LIBRARY NAMES ${dovi_LIBDOVI_STATIC_NAMES}
        HINTS ${DOVI_PC_LIBRARY_DIRS})
else()
    find_library(dovi_LIBRARY NAMES ${dovi_LIBDOVI_SHARED_NAMES}
        HINTS ${DOVI_PC_LIBRARY_DIRS})
endif()
find_path(dovi_INCLUDE_DIR NAMES dovi/rpu_parser.h libdovi/rpu_parser.h
    HINTS ${DOVI_PC_INCLUDE_DIRS} ${dovi_LIBDOVI_VENDORED_INCLUDE})

if(dovi_LIBRARY)
    if(NOT dovi_INCLUDE_DIR)
        set(dovi_INCLUDE_DIR ${dovi_LIBDOVI_VENDORED_INCLUDE})
    endif()

    if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
        if(WIN32)
            dovibaker_is_pe_import_library(dovi_IS_IMPORT_LIB ${dovi_LIBRARY})
            if(dovi_IS_IMPORT_LIB)
                message(FATAL_ERROR
                    "${dovi_LIBRARY} is the import library of the libdovi DLL, not a"
                    " static library. Configure with -DDOVIBAKER_LIBDOVI_LINK=shared"
                    " to use the DLL, or install the static library too.")
            endif()
        endif()
        dovibaker_rust_native_libs(dovi_NATIVE_LIBS "${DOVIBAKER_LIBDOVI_RUST_TRIPLE}")
        dovibaker_native_libs_fallback()
        dovibaker_define_libdovi(STATIC ${dovi_LIBRARY} ${dovi_LIBRARY} ${dovi_INCLUDE_DIR})
    else()
        get_filename_component(dovi_LIB_DIR ${dovi_LIBRARY} DIRECTORY)
        dovibaker_dovi_runtime(dovi_RUNTIME ${dovi_LIBRARY} ${dovi_LIB_DIR})
        dovibaker_define_libdovi(SHARED ${dovi_LIBRARY} ${dovi_RUNTIME} ${dovi_INCLUDE_DIR})
    endif()
    message(STATUS "libdovi: using ${dovi_LIBRARY} (${DOVIBAKER_LIBDOVI_LINK})")
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

# crate type, artifact names and the rustc arguments of the final crate
set(dovi_LIBDOVI_SONAME_FLAG "")
set(dovi_CARGO_RUSTC_ARGS "")
if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
    set(dovi_CARGO_CRATE_TYPE staticlib)
    if(MSVC)
        set(dovi_ARCHIVE ${dovi_LIBDOVI_DIR}/${dovi_LIBDOVI_CRATE_NAME}.lib)
    else()
        set(dovi_ARCHIVE ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.a)
    endif()
    set(dovi_IMPLIB ${dovi_ARCHIVE})
    set(dovi_RUNTIME ${dovi_ARCHIVE})
    if(DOVIBAKER_LIBDOVI_IMPLIB)
        message(WARNING "libdovi: DOVIBAKER_LIBDOVI_IMPLIB is only used by"
            " DOVIBAKER_LIBDOVI_LINK=shared, ignoring it")
    endif()
    # The archive carries the symbol and debug tables of the whole crate; the
    # products are stripped on their own, so keep the archive small. Arguments
    # for rustc go after "--", in RUSTFLAGS they would hit the build scripts.
    set(dovi_CARGO_RUSTC_ARGS -- -Cstrip=symbols)
elseif(WIN32)
    set(dovi_CARGO_CRATE_TYPE cdylib)
    set(dovi_RUNTIME ${dovi_LIBDOVI_DIR}/${dovi_LIBDOVI_CRATE_NAME}.dll)
    if(MSVC)
        set(dovi_IMPLIB ${dovi_LIBDOVI_DIR}/${dovi_LIBDOVI_CRATE_NAME}.dll.lib)
    else()
        set(dovi_IMPLIB ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.dll.a)
    endif()
elseif(APPLE)
    set(dovi_CARGO_CRATE_TYPE cdylib)
    set(dovi_RUNTIME ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.dylib)
    set(dovi_IMPLIB ${dovi_RUNTIME})
else()
    set(dovi_CARGO_CRATE_TYPE cdylib)
    set(dovi_RUNTIME ${dovi_LIBDOVI_DIR}/lib${dovi_LIBDOVI_CRATE_NAME}.so)
    set(dovi_IMPLIB ${dovi_RUNTIME})
    # Without a SONAME the plugins would record the absolute build path.
    set(dovi_LIBDOVI_SONAME_FLAG "-Clink-arg=-Wl,-soname,lib${dovi_LIBDOVI_CRATE_NAME}.so")
endif()

if(DOVIBAKER_LIBDOVI_IMPLIB AND NOT DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
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
            --crate-type ${dovi_CARGO_CRATE_TYPE}
            ${dovi_CARGO_EXTRA_ARGS}
            ${dovi_CARGO_RUSTC_ARGS}
    BYPRODUCTS ${dovi_CARGO_BYPRODUCTS}
    COMMENT "Building libdovi with cargo (${DOVIBAKER_LIBDOVI_RUST_TRIPLE}, release, ${dovi_CARGO_CRATE_TYPE})"
    VERBATIM)
set_target_properties(dovi_cargo PROPERTIES FOLDER thirdparty)
set(DOVI_CARGO_TARGET dovi_cargo CACHE INTERNAL "cargo target which builds libdovi")

if(DOVIBAKER_LIBDOVI_LINK STREQUAL "static")
    dovibaker_rust_native_libs(dovi_NATIVE_LIBS "${DOVIBAKER_LIBDOVI_RUST_TRIPLE}")
    dovibaker_native_libs_fallback()
    dovibaker_define_libdovi(STATIC ${dovi_IMPLIB} ${dovi_RUNTIME} ${dovi_LIBDOVI_VENDORED_INCLUDE})
    message(STATUS "libdovi: building ${dovi_ARCHIVE} from the submodule"
        " (cargo, release, static - no runtime library to ship)")
else()
    dovibaker_define_libdovi(SHARED ${dovi_IMPLIB} ${dovi_RUNTIME} ${dovi_LIBDOVI_VENDORED_INCLUDE})
    message(STATUS "libdovi: building ${dovi_RUNTIME} from the submodule (cargo, release, shared)")
endif()
add_dependencies(dovi::dovi dovi_cargo)
