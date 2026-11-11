# CMakeLists.txt -- top-level build system for ORIGAMI
#
# The C sources live in code/ (which also has the original Makefile). This file
# sits one level above that directory and builds the same targets:
#   all       : vozinit voz1b1 voztie origamitag origamihalo origamicen origamienv
#   origami   : origamitag origamihalo origamicen origamienv
#   voboz     : vozinit voz1b1 voztie
#
# Typical use (out-of-source build, from the repository root):
#   cmake -S . -B build
#   cmake --build build                      # or: cmake --build build --target origami
#
# Options:
#   -DORIGAMI_USE_SYSTEM_QHULL=ON  Link against an installed Qhull (QHULL_ROOT
#                                  may be given as a hint) instead of building
#                                  the Qhull 2012.1 tarball shipped in the repo.
#   -DORIGAMI_INT_IDS=ON           Use readfiles_int.c (4-byte Gadget particle
#                                  IDs) instead of readfiles.c (8-byte IDs).
#   -DORIGAMI_BUILD_TESTS=ON       Build the C unit tests in tests/; run them
#                                  with `ctest --test-dir build` afterwards.
#   -DOpenMP_ROOT=$(brew --prefix libomp)
#                                  This may be needed on OS X to ensure use of
#                                  an OpenMP installation from homebrew
#
# Note: the Makefile in code/ hard-codes a Qhull location under /users/bfalck;
# that path is intentionally not used here.

cmake_minimum_required(VERSION 3.18)
project(origami LANGUAGES C)

# Directory holding the C sources and the original Makefile
set(ORIGAMI_CODE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/code")

# Match the Makefile's -O2 optimisation level when no build type is requested.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_C_FLAGS_RELEASE "-O2")

option(ORIGAMI_USE_SYSTEM_QHULL "Use an installed Qhull instead of the bundled 2012.1 tarball" OFF)
option(ORIGAMI_INT_IDS "Read Gadget particle IDs as 4-byte integers (readfiles_int.c)" OFF)
option(ORIGAMI_BUILD_TESTS "Build the C unit tests in tests/ (run with ctest)" OFF)
set(QHULL_ROOT "" CACHE PATH "Optional hint: prefix of an installed Qhull (with ORIGAMI_USE_SYSTEM_QHULL=ON)")

find_package(OpenMP REQUIRED)
find_library(M_LIBRARY m REQUIRED)

# ---------------------------------------------------------------------------
# Qhull
# ---------------------------------------------------------------------------
if(ORIGAMI_USE_SYSTEM_QHULL)
  set(_qhull_hints "")
  if(QHULL_ROOT)
    set(_qhull_hints HINTS "${QHULL_ROOT}" "${QHULL_ROOT}/include" "${QHULL_ROOT}/lib")
  endif()

  find_path(QHULL_INCLUDE_DIR qhull_a.h
    ${_qhull_hints}
    PATH_SUFFIXES libqhull qhull include/libqhull include/qhull src/libqhull)
  find_library(QHULL_LIBRARY NAMES qhull qhullstatic
    ${_qhull_hints}
    PATH_SUFFIXES lib)

  if(NOT QHULL_INCLUDE_DIR OR NOT QHULL_LIBRARY)
    message(FATAL_ERROR
      "Could not find Qhull (qhull_a.h and libqhull). Set QHULL_ROOT to its install prefix.")
  endif()

  add_library(qhull_ext INTERFACE)
  target_include_directories(qhull_ext INTERFACE "${QHULL_INCLUDE_DIR}")
  target_link_libraries(qhull_ext INTERFACE "${QHULL_LIBRARY}")
  set(QHULL_TARGET qhull_ext)
else()
  # Unpack the Qhull 2012.1 source tarball shipped in the repository root
  # into the build tree, then compile libqhull directly. This avoids relying
  # on the Qhull install paths hard-coded in the Makefile.
  set(QHULL_TARBALL "${CMAKE_CURRENT_SOURCE_DIR}/qhull-2012.1-src.tgz")
  set(QHULL_SRC_ROOT "${CMAKE_BINARY_DIR}/_deps/qhull-2012.1")
  set(QHULL_LIBDIR "${QHULL_SRC_ROOT}/src/libqhull")

  if(NOT EXISTS "${QHULL_LIBDIR}/qhull_a.h")
    if(NOT EXISTS "${QHULL_TARBALL}")
      message(FATAL_ERROR "Qhull tarball not found: ${QHULL_TARBALL}")
    endif()
    message(STATUS "Unpacking bundled Qhull: ${QHULL_TARBALL}")
    file(ARCHIVE_EXTRACT INPUT "${QHULL_TARBALL}" DESTINATION "${CMAKE_BINARY_DIR}/_deps")
  endif()

  # libqhull source list (from qhull's own CMakeLists.txt, version 2012.1)
  set(QHULL_SOURCES
    global.c stat.c geom2.c poly2.c merge.c libqhull.c geom.c poly.c
    qset.c mem.c random.c usermem.c userprintf.c io.c user.c
    rboxlib.c userprintf_rbox.c)
  list(TRANSFORM QHULL_SOURCES PREPEND "${QHULL_LIBDIR}/")

  add_library(qhull STATIC ${QHULL_SOURCES})
  target_include_directories(qhull PUBLIC "${QHULL_LIBDIR}")
  # Qhull 2012.1 is old C; silence warnings that originate in its own sources.
  target_compile_options(qhull PRIVATE -w)
  set(QHULL_TARGET qhull)
endif()

# ---------------------------------------------------------------------------
# Shared particle-reading code
# ---------------------------------------------------------------------------
if(ORIGAMI_INT_IDS)
  set(READFILES_SRC "${ORIGAMI_CODE_DIR}/readfiles_int.c")
else()
  set(READFILES_SRC "${ORIGAMI_CODE_DIR}/readfiles.c")
endif()

add_library(origami_readfiles STATIC ${READFILES_SRC})
target_include_directories(origami_readfiles PUBLIC "${ORIGAMI_CODE_DIR}")
target_link_libraries(origami_readfiles PUBLIC OpenMP::OpenMP_C)

# Helper: one executable per Makefile target, linked with the shared
# particle reader, OpenMP (-fopenmp) and libm. Source files are given
# relative to code/.
function(origami_add_tool name)
  set(_srcs "")
  foreach(_f IN LISTS ARGN)
    list(APPEND _srcs "${ORIGAMI_CODE_DIR}/${_f}")
  endforeach()
  add_executable(${name} ${_srcs})
  target_include_directories(${name} PRIVATE "${ORIGAMI_CODE_DIR}")
  target_link_libraries(${name} PRIVATE origami_readfiles OpenMP::OpenMP_C ${M_LIBRARY})
  install(TARGETS ${name} RUNTIME DESTINATION bin)
endfunction()

# ---------------------------------------------------------------------------
# VOBOZ tessellation tools (voz1b1 needs Qhull)
# ---------------------------------------------------------------------------
origami_add_tool(vozinit vozinit.c)
# vozinit.c includes qhull_a.h (the Makefile passes QINC for it), so it needs
# the Qhull include path even though it does not call Qhull routines.
target_link_libraries(vozinit PRIVATE ${QHULL_TARGET})

origami_add_tool(voz1b1 voz1b1.c vozutil.c)
target_link_libraries(voz1b1 PRIVATE ${QHULL_TARGET})

origami_add_tool(voztie voztie.c)

# ---------------------------------------------------------------------------
# ORIGAMI morphology / halo / environment tools
# ---------------------------------------------------------------------------
origami_add_tool(origamitag origamitag.c)
origami_add_tool(origamihalo origamihalo.c)
origami_add_tool(origamienv origamienv.c)
origami_add_tool(origamicen origamicen.c)

# ---------------------------------------------------------------------------
# Convenience targets mirroring the Makefile (executables are already part of
# the default "all" target, so these only select subsets)
# ---------------------------------------------------------------------------
add_custom_target(origami
  DEPENDS origamitag origamihalo origamicen origamienv
  COMMENT "Building ORIGAMI tools")

add_custom_target(voboz
  DEPENDS vozinit voz1b1 voztie
  COMMENT "Building VOBOZ tessellation tools")

# ---------------------------------------------------------------------------
# Unit tests (off by default). Build with -DORIGAMI_BUILD_TESTS=ON, then run
# `ctest --test-dir build --output-on-failure`.
# ---------------------------------------------------------------------------
if(ORIGAMI_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
