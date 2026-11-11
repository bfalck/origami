# CMakeLists.txt for the ORIGAMI C unit tests (tests/).
#
# Enabled from the top-level build with -DORIGAMI_BUILD_TESTS=ON, then run with
#   ctest --test-dir build --output-on-failure
#
# The tests need only a C compiler, OpenMP (the tools use it) and libm. They
# do not need Qhull, so they never touch the Delaunay/VOBOZ code.
#
# Each test that covers a tool's helpers #includes that tool's .c file, which
# gives the test direct access to its static-free helper functions.

set(ORIGAMI_EXAMPLE_PARAMS "${ORIGAMI_CODE_DIR}/parameters_example.txt")

# Builds one test executable and registers it with CTest.
#   name   : test name, also the source file name (tests/<name>.c)
#   reader : which readfiles implementation to link against
function(origami_add_unit_test name reader)
  add_executable(${name} ${name}.c ${reader})
  target_include_directories(${name} PRIVATE
    "${ORIGAMI_CODE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(${name} PRIVATE
    "ORIGAMI_EXAMPLE_PARAMS=\"${ORIGAMI_EXAMPLE_PARAMS}\"")
  target_link_libraries(${name} PRIVATE OpenMP::OpenMP_C ${M_LIBRARY})
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# origamitag.c uses whichever reader the build selected. The readfiles tests
# and origamidiv.c always use the 8-byte reader (readfiles.c): readfiles_int.c
# has no `nbuf` global and no divread(), so it does not match the current
# parameter-file format or origamidiv.c.
origami_add_unit_test(test_origamitag ${READFILES_SRC})
origami_add_unit_test(test_readfiles "${ORIGAMI_CODE_DIR}/readfiles.c")
origami_add_unit_test(test_origamidiv "${ORIGAMI_CODE_DIR}/readfiles.c")
