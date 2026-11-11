/*-----------------------------------------------------------------------------
 * test_harness.h -- a minimal, dependency-free C unit-test harness.
 *
 * Usage:
 *   #include "test_harness.h"
 *   static void test_something(void) { CHECK_INT_EQ(f(2), 4); }
 *   int main(void) { RUN_TEST(test_something); TEST_SUMMARY("name"); }
 *
 * Each CHECK_* records a failure and keeps going, so one run reports every
 * failing check. TEST_SUMMARY returns 0 only when every check passed, which
 * is what CTest uses to decide pass or fail.
 *---------------------------------------------------------------------------*/
#ifndef ORIGAMI_TEST_HARNESS_H
#define ORIGAMI_TEST_HARNESS_H

#include <math.h>
#include <stdio.h>
#include <string.h>

static int tests_run = 0;
static int tests_failed = 0;

/* Check that a boolean expression is true. */
#define CHECK(cond) do {                                                  \
    tests_run++;                                                          \
    if (!(cond)) {                                                        \
      tests_failed++;                                                     \
      fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);   \
    }                                                                     \
  } while (0)

/* Check that two integers are equal. */
#define CHECK_INT_EQ(actual, expected) do {                               \
    long _a = (long)(actual);                                             \
    long _e = (long)(expected);                                           \
    tests_run++;                                                          \
    if (_a != _e) {                                                       \
      tests_failed++;                                                     \
      fprintf(stderr, "  FAIL %s:%d: %s = %ld, expected %ld\n",           \
              __FILE__, __LINE__, #actual, _a, _e);                       \
    }                                                                     \
  } while (0)

/* Check that two floating-point values agree to within tol. */
#define CHECK_NEAR(actual, expected, tol) do {                            \
    double _a = (double)(actual);                                         \
    double _e = (double)(expected);                                       \
    tests_run++;                                                          \
    if (fabs(_a - _e) > (tol)) {                                          \
      tests_failed++;                                                     \
      fprintf(stderr, "  FAIL %s:%d: %s = %.9g, expected %.9g (tol %g)\n", \
              __FILE__, __LINE__, #actual, _a, _e, (double)(tol));        \
    }                                                                     \
  } while (0)

/* Check that two C strings are identical. */
#define CHECK_STR_EQ(actual, expected) do {                               \
    const char *_a = (actual);                                            \
    const char *_e = (expected);                                          \
    tests_run++;                                                          \
    if (_a == NULL || strcmp(_a, _e) != 0) {                              \
      tests_failed++;                                                     \
      fprintf(stderr, "  FAIL %s:%d: %s = \"%s\", expected \"%s\"\n",     \
              __FILE__, __LINE__, #actual, _a ? _a : "(null)", _e);       \
    }                                                                     \
  } while (0)

/* Run one test function and report its status. */
#define RUN_TEST(fn) do {                                                 \
    int _before = tests_failed;                                           \
    fprintf(stderr, "[ RUN  ] %s\n", #fn);                                \
    fn();                                                                 \
    fprintf(stderr, "[ %s ] %s\n",                                        \
            (tests_failed == _before) ? " OK " : "FAIL", #fn);            \
  } while (0)

/* Print totals and return the process exit code. */
#define TEST_SUMMARY(suite) do {                                          \
    fprintf(stderr, "%s: %d checks, %d failed\n",                         \
            (suite), tests_run, tests_failed);                            \
    return (tests_failed != 0) ? 1 : 0;                                   \
  } while (0)

#endif /* ORIGAMI_TEST_HARNESS_H */
