/*-----------------------------------------------------------------------------
 * test_origamitag.c -- unit tests for the helper routines in code/origamitag.c.
 *
 * origamitag.c defines its own main(). We rename it while including the file,
 * so this translation unit can supply the test main() instead.
 *---------------------------------------------------------------------------*/
#include <limits.h>
#include "test_harness.h"

#define main origamitag_main
#include "origamitag.c"
#undef main

static void test_isneg(void) {
  CHECK_INT_EQ(isneg(-1), 1);
  CHECK_INT_EQ(isneg(INT_MIN), 1);
  CHECK_INT_EQ(isneg(0), 0);
  CHECK_INT_EQ(isneg(1), 0);
  CHECK_INT_EQ(isneg(INT_MAX), 0);
}

static void test_par_axis_ordering(void) {
  const int ng = 8;
  CHECK_INT_EQ(par(0, 0, 0, ng), 0);
  CHECK_INT_EQ(par(1, 0, 0, ng), 1);             /* x changes fastest */
  CHECK_INT_EQ(par(0, 1, 0, ng), ng);            /* then y */
  CHECK_INT_EQ(par(0, 0, 1, ng), ng * ng);       /* then z */
  CHECK_INT_EQ(par(ng - 1, ng - 1, ng - 1, ng), ng * ng * ng - 1);
}

static void test_par_is_bijection(void) {
  enum { NG = 5, NTOT = NG * NG * NG };
  int seen[NTOT];
  int i, j, k, n_unique = 0;

  memset(seen, 0, sizeof(seen));
  for (k = 0; k < NG; k++)
    for (j = 0; j < NG; j++)
      for (i = 0; i < NG; i++) {
        int idx = par(i, j, k, NG);
        CHECK(idx >= 0 && idx < NTOT);
        if (idx >= 0 && idx < NTOT) {
          CHECK(seen[idx] == 0);   /* no two grid points share an index */
          seen[idx] = 1;
        }
      }
  for (i = 0; i < NTOT; i++) n_unique += seen[i];
  CHECK_INT_EQ(n_unique, NTOT);    /* every index is hit exactly once */
}

static void test_goodmod_periodic_wrap(void) {
  CHECK_INT_EQ(goodmod(2, 8), 2);    /* in range: unchanged */
  CHECK_INT_EQ(goodmod(8, 8), 0);    /* one period above wraps to 0 */
  CHECK_INT_EQ(goodmod(-1, 8), 7);   /* one step below wraps to the top */
}

static void test_max_macro(void) {
  CHECK_INT_EQ(max(3, 7), 7);
  CHECK_INT_EQ(max(7, 3), 7);
  CHECK_INT_EQ(max(-2, -5), -2);
}

int main(void) {
  RUN_TEST(test_isneg);
  RUN_TEST(test_par_axis_ordering);
  RUN_TEST(test_par_is_bijection);
  RUN_TEST(test_goodmod_periodic_wrap);
  RUN_TEST(test_max_macro);
  TEST_SUMMARY("test_origamitag");
}
