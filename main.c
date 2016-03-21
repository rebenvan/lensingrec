#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "def.h"

int main(int argc, char *argv[])
{
  int ii, ll;
  char filename[128];

  float min, max, ave;

  fprintf(stderr, "\nReading basic data from file...\n");

  read_data();	// Read the dark mattter, source galaxies and image galaxies

  fprintf(stderr, "completed!\n");

  if (METHOD == 1 || METHOD == 3) { // Minimal Variance Linear Estimator
    fprintf(stderr, "\nUsing method 1\n");
    int iter;
    int conv, lastconv;
    float *kt;
    iter = 0; conv = -1;
    kt = (float *)malloc(mtp*sizeof(float));
    assert(kt != NULL);
    do {
      fprintf(stderr, "round:%d\n", iter+1);
      lastconv = conv;

      cal_bias(iter);
      cal_weight(iter);
      cal_kappa(iter);

      if (iter > 1) {
	conv = 0;
	for (ll = 1; ll <= mtp; ll++) {
	  if (fabs(*(kap1+ll-1)-*(kt+ll-1))/(*(kt+ll-1)))
	    conv++;
	}
      }

      for (ll = 1; ll <= mtp; ll++)
	*(kt+ll-1) = *(kap1+ll-1);

      //fprintf(stderr, "\tconv=%d/%d", conv, mtp);

      iter++;

    } while(conv > lastconv || conv < 0.95*mtp);
    if (kalm != NULL) free(kalm);
    free(wgt[0]);

    fprintf(stderr, "completed\n");
  }

  if (METHOD == 2 || METHOD == 3) { // Direct Determination
    fprintf(stderr, "\nUsing method 2\n");

    cal_kappa2();

    fprintf(stderr, "bias:\n");
    for (ii = 0; ii < BIN; ii++) {
      min = *bis2[ii]; max = *bis2[ii]; ave = 0.F;
      for (ll = 1; ll <= mtp; ll++) {
	ave += *(bis2[ii]+ll-1);
	if (*(bis2[ii]+ll-1) < min) min = *(bis2[ii]+ll-1);
	if (*(bis2[ii]+ll-1) > max) max = *(bis2[ii]+ll-1);
      }
      fprintf(stderr, "\t\t%f|%f|%f\n", min, ave/mtp, max);
    }

    fprintf(stderr, "kappa:\n");
    min = *kap2; max = *kap2;
    for (ll = 1; ll <= mtp; ll++) {
      if (*(kap2+ll-1) < min) min = *(kap2+ll-1);
      if (*(kap2+ll-1) > max) max = *(kap2+ll-1);
    }
    fprintf(stderr, "\t\t%.2e|%.2e\n", min, max);

    fprintf(stderr, "completed\n");
  }

  //write_file("bias");
  //write_file("kappa");
  write_bias();
	write_kappa();

  free(cmb);
  for (ii = 0; ii < BIN2; ii++)
    free(cab[ii]);
  free(bis1[0]); free(kap1);
  free(bis2[0]); free(kap2);

  return 0;
}
