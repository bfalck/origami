/*-----------------------------------------------------------------------------
 * test_readfiles.c -- unit tests for parameter-file parsing and the binary
 * position reader in code/readfiles.c (or readfiles_int.c when built with
 * -DORIGAMI_INT_IDS=ON; both provide read_parameters() and posread()).
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "test_harness.h"
#include "params.h"   /* extern declarations of the global parameters */

/* Declared in readfiles.c / readfiles_int.c. */
void read_parameters(char *filename);
int posread(char *posfile, float ***p, float fact);

/* A parameter file laid out exactly like parameters_example.txt, with values
 * chosen so that none of them match the defaults. The reader skips one
 * header line and a few comment/blank lines between groups, so the layout
 * matters. */
static const char *fixture_params =
  "# General parameters\n"
  "posfile\t/data/test.pos\n"
  "outdir\t/out/\n"
  "taglabel\tn32\n"
  "boxsize\t50.5\n"
  "np1d\t32\n"
  "nsplit\t2\n"
  "numfiles\t0\n"
  "\n"
  "# Used by origamidiv_tag\n"
  "nbuf\t4\n"
  "\n"
  "# Used by vozinit\n"
  "buffer\t0.25\n"
  "ndiv\t3\n"
  "\n"
  "# Used by origamihalo\n"
  "volcut\t0.5\n"
  "npmin\t7\n"
  "halolabel\th32\n";

static void test_read_parameters_fixture(void) {
  char path[] = "/tmp/origami_params_XXXXXX";
  int fd = mkstemp(path);
  FILE *f;

  CHECK(fd >= 0);
  if (fd < 0) return;
  f = fdopen(fd, "w");
  CHECK(f != NULL);
  if (f == NULL) { close(fd); unlink(path); return; }
  fputs(fixture_params, f);
  fclose(f);

  read_parameters(path);
  unlink(path);

  CHECK_STR_EQ(posfile, "/data/test.pos");
  CHECK_STR_EQ(outdir, "/out/");
  CHECK_STR_EQ(taglabel, "n32");
  CHECK_NEAR(boxsize, 50.5, 1e-4);
  CHECK_INT_EQ(np1d, 32);
  CHECK_INT_EQ(nsplit, 2);
  CHECK_INT_EQ(numfiles, 0);
  CHECK_INT_EQ(nbuf, 4);
  CHECK_NEAR(buffer, 0.25, 1e-6);
  CHECK_INT_EQ(numdiv, 3);
  CHECK_NEAR(volcut, 0.5, 1e-6);
  CHECK_INT_EQ(npmin, 7);
  CHECK_STR_EQ(halolabel, "h32");
}

static void test_read_parameters_shipped_example(void) {
  /* The example file shipped with the code must always parse. */
  read_parameters(ORIGAMI_EXAMPLE_PARAMS);

  CHECK_STR_EQ(taglabel, "n256");
  CHECK_NEAR(boxsize, 100000.0, 1e-2);
  CHECK_INT_EQ(np1d, 256);
  CHECK_INT_EQ(nsplit, 1);
  CHECK_INT_EQ(numfiles, 0);
  CHECK_INT_EQ(nbuf, 64);
  CHECK_NEAR(buffer, 0.1, 1e-6);
  CHECK_INT_EQ(numdiv, 2);
  CHECK_NEAR(volcut, 0.01, 1e-6);
  CHECK_INT_EQ(npmin, 20);
  CHECK_STR_EQ(halolabel, "n256d100");
}

static void test_posread_binary_roundtrip(void) {
  /* Format: int count, then all x, all y, all z as 4-byte floats. */
  const int np = 4;
  const float x[4] = {0.0f, 1.5f, -2.25f, 10.0f};
  const float y[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  const float z[4] = {-1.0f, 0.5f, 7.0f, 100.0f};
  const float fact = 2.0f;
  char path[] = "/tmp/origami_pos_XXXXXX";
  int fd = mkstemp(path);
  FILE *f;
  float **p = NULL;
  int n, i;

  CHECK(fd >= 0);
  if (fd < 0) return;
  f = fdopen(fd, "wb");
  CHECK(f != NULL);
  if (f == NULL) { close(fd); unlink(path); return; }
  fwrite(&np, sizeof(int), 1, f);
  fwrite(x, sizeof(float), np, f);
  fwrite(y, sizeof(float), np, f);
  fwrite(z, sizeof(float), np, f);
  fclose(f);

  n = posread(path, &p, fact);
  unlink(path);

  CHECK_INT_EQ(n, np);
  if (n == np && p != NULL) {
    for (i = 0; i < np; i++) {
      CHECK_NEAR(p[i][0], x[i] * fact, 1e-6);   /* scaled by fact */
      CHECK_NEAR(p[i][1], y[i] * fact, 1e-6);
      CHECK_NEAR(p[i][2], z[i] * fact, 1e-6);
    }
  }
  if (p != NULL) {
    for (i = 0; i < np; i++) free(p[i]);
    free(p);
  }
}

int main(void) {
  RUN_TEST(test_read_parameters_fixture);
  RUN_TEST(test_read_parameters_shipped_example);
  RUN_TEST(test_posread_binary_roundtrip);
  TEST_SUMMARY("test_readfiles");
}
