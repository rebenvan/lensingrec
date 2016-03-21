#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "def.h"

void write_bias()
{
  int ii;
  int bins;
	char filename[128];
  FILE *pfile;

  bins = BIN*mtp;
  if (METHOD == 1 || METHOD == 3) {
    sprintf(filename, "%s/result/info_bias1.dat", path);
    pfile = fopen(filename, "wb");
    assert(pfile != NULL);
    fwrite(&bins, sizeof(int), (size_t)1, pfile);
    for (ii = 0; ii < BIN; ii++)
      fwrite(bis1[ii], sizeof(float), (size_t)mtp, pfile);
    fwrite(&bins, sizeof(int), (size_t)1, pfile);
    fclose(pfile);
  }
  if (METHOD == 2 || METHOD == 3) {
    sprintf(filename, "%s/result/info_bias2.dat", path);
    pfile = fopen(filename, "wb");
    assert(pfile != NULL);
    fwrite(&bins, sizeof(int), (size_t)1, pfile);
    for (ii = 0; ii < BIN; ii++)
      fwrite(bis2[ii], sizeof(float), (size_t)mtp, pfile);
    fwrite(&bins, sizeof(int), (size_t)1, pfile);
    fclose(pfile);
  }
}

void write_kappa()
{
	char filename[128];
  FILE *pfile;

  if (METHOD == 1 || METHOD == 3) {
    sprintf(filename, "%s/result/info_kappa1.dat", path);
    pfile = fopen(filename, "wb");
    assert(pfile != NULL);
    fwrite(&mtp, sizeof(int), (size_t)1, pfile);
    fwrite(kap1, sizeof(float), (size_t)mtp, pfile);
    fwrite(&mtp, sizeof(int), (size_t)1, pfile);
    fclose(pfile);
  }

  if (METHOD == 2 || METHOD == 3) {
    sprintf(filename, "%s/result/info_kappa2.dat", path);
    pfile = fopen(filename, "wb");
    assert(pfile != NULL);
    fwrite(&mtp, sizeof(int), (size_t)1, pfile);
    fwrite(kap2, sizeof(float), (size_t)mtp, pfile);
    fwrite(&mtp, sizeof(int), (size_t)1, pfile);
    fclose(pfile);
  }
}
