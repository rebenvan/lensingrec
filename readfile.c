#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "def.h"

void show(int mtp, float *ps);

void read_data()
{
  char filename[128];
  FILE *pf;
  int ref;

  int ii, jj, ll;
  int ind, bin;

  fprintf(stderr, "Read source galaxy power spectrum...");
  cmb = read_cl(-1, -1);
  fprintf(stderr, "mtp=%d", mtp);
  show(mtp, cmb);
  fprintf(stderr, "done!\n");

  fprintf(stderr, "Read image galaxy power spectrum...\n");
  ind = 0;
  for (ii = 0; ii < BIN; ii++) {
    for (jj = ii; jj < BIN; jj++) {
      fprintf(stderr, "...%d%d\t", ii, jj);
      cab[ind] = read_cl(ii, jj);
      show(mtp, cab[ind]);
      ind++;
    }
  }
  assert(ind == BIN2);
  fprintf(stderr, "done!\n");
	ind = 0;
	for (ii = 0; ii < BIN; ii++) {
		for (jj = ii; jj < BIN; jj++) {
			for (ll = 0; ll < mtp; ll++) {
				cab[ind][ll] = cmb[ll]*(1.0+(ii+jj+2)/2.0/BIN);
			}
			ind++;
		}
	}

  fprintf(stderr, "Read average number density...");
  sprintf(filename, "%s/ave.dat", path);
  pf = fopen(filename, "rb");
  assert(pf != NULL);
  ref = fread(&bin,sizeof(int),(size_t)1,pf);
  assert(bin == BIN);
  ref = fread(&ave,sizeof(float),(size_t)BIN,pf);
  ref = fread(&bin,sizeof(int),(size_t)1,pf);
  assert(bin == BIN);
  fclose(pf);
	show(BIN, ave);
  fprintf(stderr, "done!\n");

  fprintf(stderr, "Read lominosity slope...\n");
	sprintf(filename, "%s/slp.dat", path);
  pf = fopen(filename, "rb");
  assert(pf != NULL);
  ref = fread(&bin,sizeof(int),(size_t)1,pf);
  assert(bin == BIN);
  ref = fread(&slp,sizeof(float),(size_t)BIN,pf);
  ref = fread(&bin,sizeof(int),(size_t)1,pf);
  assert(bin == BIN);
  fclose(pf);
	show(BIN, slp);
  fprintf(stderr, "done!\n");
}

void show(int mtp, float *ps)
{
  int ii;
  float min, max;

  min = ps[0]; max = ps[0];
  for (ii = 1; ii < mtp; ii++) {
    if (ps[ii] < min) min = ps[ii];
    if (ps[ii] > max) max = ps[ii];
  }
  fprintf(stderr, "\t%.2e|%.2e\n", min, max);
}
