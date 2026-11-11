/*-----------------------------------------------------------------------------
 * test_origamidiv.c -- unit tests for the helper routines in code/origamidiv.c.
 *
 * origamidiv.c has its own isneg() and par() (separate from origamitag.c), so
 * they are tested here on their own. Its main() is renamed while including
 * the file so this translation unit can supply the test main().
 *---------------------------------------------------------------------------*/
#include <limits.h>
#include "test_harness.h"

#define main origamidiv_main
#include "origamidiv.c"
#undef main

static void test_isneg(void) {
  CHECK_INT_EQ(isneg(-3), 1);
  CHECK_INT_EQ(isneg(INT_MIN), 1);
  CHECK_INT_EQ(isneg(0), 0);
  CHECK_INT_EQ(isneg(42), 0);
}

static void test_par_axis_ordering(void) {
  const int ng = 4;
  CHECK_INT_EQ(par(0, 0, 0, ng), 0);
  CHECK_INT_EQ(par(1, 0, 0, ng), 1);
  CHECK_INT_EQ(par(0, 1, 0, ng), ng);
  CHECK_INT_EQ(par(0, 0, 1, ng), ng * ng);
  CHECK_INT_EQ(par(ng - 1, ng - 1, ng - 1, ng), ng * ng * ng - 1);
}

int main(void) {
  RUN_TEST(test_isneg);
  RUN_TEST(test_par_axis_ordering);
  TEST_SUMMARY("test_origamidiv");
}
