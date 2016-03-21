#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "healpix_utils.h"

void cal_darkmatter_deltat();
void cal_source_deltat();
void cal_image_deltat();

long Order;
int Plane;
int Bins;
char path[128];

void main(int argc, char *argv[])
{
	int ii;
	int dk=0,sg=0,ig=0;
	if (argc >= 4) {
		Plane = atoi(argv[1]);
		Order = atol(argv[2]);
		Bins = atoi(argv[3]);
		sprintf(path, "/huawei/osv1/FAN/kappadata/sphere_%04d", Plane);
		for (ii = 4; ii < argc; ii++) {
			if (strcmp(argv[ii], "dk") == 0) dk = 1;
			if (strcmp(argv[ii], "sg") == 0) sg = 1;
			if (strcmp(argv[ii], "ig") == 0) ig = 1;
		}
	}
	else {
		fprintf(stderr, "Please enter the args:\n");
		fprintf(stderr, "arg 1: Plane\n");
		fprintf(stderr, "arg 2: Order\n");
		fprintf(stderr, "arg 3: Bins\n");
		fprintf(stderr, "arg 4-6: dk/sg/ig if needed\n");
		return;
	}

	if (dk == 1)
		cal_darkmatter_deltat();

	if (sg == 1)
		cal_source_deltat();

	if (ig == 1)
		cal_image_deltat();

	return;
}

void cal_darkmatter_deltat()
{
	long npix, pix, dummy;
  float *Mapin;
  char filename[256];
  FILE *fpr;
  int ref;

  long ii, jj;
  float *dpix, tmp, ave, min, max;

  npix = order2npix(13L);
  Mapin = (float *)malloc(sizeof(float)*npix);
  assert(Mapin != NULL);

  fprintf(stderr, "Read the healpix point of the dark matter...\n");
  sprintf(filename, "/huawei/osv1/wcl/RAYTRACING/L500/lensplanes/lensplane%04d.dat", Plane);
  fpr = fopen(filename, "rb");
  assert(fpr != NULL);
  ref = fread(&dummy, sizeof(long), 1, fpr);
  assert(dummy == npix);
  ref = fread(Mapin, sizeof(float), npix, fpr);
  ref = fread(&dummy, sizeof(long), 1, fpr);
  assert(dummy == npix);
  fclose(fpr);
  fprintf(stderr, "Read successfully!\n\n");

	npix = order2npix(Order);
  assert(Order <= 13L);
  pix = 1L<<(2*(13L-Order));
  dpix = (float *)calloc(npix, sizeof(float));
  assert(dpix != NULL);
  ave = 0.0F;
  for (ii = 0; ii < npix; ii++) {
    tmp = 0.0F;
    for (jj = 0; jj < pix; jj++)
      tmp += Mapin[ii*pix+jj];
    ave += tmp;
    dummy = nest2ring(ii, Order);
    dpix[dummy] = tmp;
  }
  ave /= npix;
  for (ii = 0; ii < npix; ii++)
    dpix[ii] = (dpix[ii]-ave)/ave;
  min = dpix[0]; max = dpix[0];
  for (ii = 1L; ii < npix; ii++) {
    if (dpix[ii] < min) min = dpix[ii];
    if (dpix[ii] > max) max = dpix[ii];
  }
  fprintf(stderr, "\tAll  dpix=(%f ~ %f)\n", min, max);

  sprintf(filename, "%s/deltat_%ld/h_darkmatter_deltat.dat", path, Order);
  fpr = fopen(filename, "wb");
  fwrite(&npix, sizeof(long), 1, fpr);
  fwrite(dpix, sizeof(float), npix, fpr);
  fwrite(&npix, sizeof(long), 1, fpr);
  fclose(fpr);

  free(Mapin); free(dpix);
  return;
}

void cal_source_deltat()
{
	long ii, jj;
	FILE *fpr;
	char filename[256];
	int ref;

	long Ngals;
	float *gals;
	float *dpix;
	double theta, phi;
	long inpix, pix[4];
	double wgt[4];
	float ave;
	float min, max;
	long npix=12L*(1L<<(2*Order));

	fprintf(stderr, "Read source galaxy from file...\n");
	sprintf(filename, "%s/galaxies.bin", path);
	fpr = fopen(filename, "rb");
	assert(fpr != NULL);
	fread(&Ngals, sizeof(long), (size_t)1, fpr);
	gals = (float *)malloc(Ngals*5*sizeof(float));
	assert(gals != NULL);
	fread(gals, sizeof(float), (size_t)Ngals*5, fpr);
	fread(&Ngals, sizeof(long), (size_t)1, fpr);
	fclose(fpr);
	fprintf(stderr, "Read completed!\n");

	fprintf(stderr, "The range of the sources...\n");
  fprintf(stderr, "Ngals=\t%ld\n", Ngals);
  min = gals[0]; max = gals[0];
  for (ii = 1L; ii < Ngals; ii++) {
    if (gals[ii*5] < min) min = gals[ii*5];
    if (gals[ii*5] > max) max = gals[ii*5];
  }
  fprintf(stderr, "\tTheta    :(%f ~ %f)\n", min, max);
  min = gals[1]; max = gals[1];
  for (ii = 1L; ii < Ngals; ii++) {
    if (gals[ii*5+1] < min) min = gals[ii*5+1];
    if (gals[ii*5+1] > max) max = gals[ii*5+1];
  }
  fprintf(stderr, "\tPhi      :(%f ~ %f)\n", min, max);
  min = gals[2]; max = gals[2];
  for (ii = 1L; ii < Ngals; ii++) {
    if (gals[ii*5+2] < min) min = gals[ii*5+2];
    if (gals[ii*5+2] > max) max = gals[ii*5+2];
  }
  fprintf(stderr, "\tLum(Ls)  :(%.2e ~ %.2e)\n", min, max);
  min = gals[3]; max = gals[3];
  for (ii = 1L; ii < Ngals; ii++) {
    if (gals[ii*5+3] < min) min = gals[ii*5+3];
    if (gals[ii*5+3] > max) max = gals[ii*5+3];
  }
  fprintf(stderr, "\tDis(Mpc) :(%.2e ~ %.2e)\n", min, max);
	min = gals[4]; max = gals[4];
  for (ii = 1L; ii < Ngals; ii++) {
    if (gals[ii*5+4] < min) min = gals[ii*5+4];
    if (gals[ii*5+4] > max) max = gals[ii*5+4];
  }
  fprintf(stderr, "\tMass(Lm) :(%.2e ~ %.2e)\n", min, max);
  fprintf(stderr, "\n");

	fprintf(stderr, "Allocate the galaxies to healpix pixels...\n");
  dpix = (float *)calloc(npix, sizeof(float));
  assert(dpix != NULL);
  for (ii = 0L; ii < Ngals; ii++) {
    theta = (double)gals[ii*5];
    phi = (double)gals[ii*5+1];
		inpix = ang2ring(theta, phi, Order);
		dpix[inpix] += 1.0F;
    //get_interpol(theta, phi, pix, wgt, Order);
    //for (jj = 0; jj < 4; jj++)
      //dpix[pix[jj]] += (float)wgt[jj];
  }
  ave = 1.F*Ngals/npix;
  for (ii = 0; ii < npix; ii++)
    dpix[ii] = (dpix[ii]-ave)/ave;
  min = dpix[0]; max = dpix[0];
  for (ii = 1L; ii < npix; ii++) {
    if (dpix[ii] < min) min = dpix[ii];
    if (dpix[ii] > max) max = dpix[ii];
  }
  fprintf(stderr, "\tAll  dpix=(%f ~ %f)\n", min, max);
  gals = NULL;

	fprintf(stderr, "Write healpix overdensity to file...\n");
  sprintf(filename, "%s/deltat_%ld/h_source_deltat.dat", path, Order);
  fpr = fopen(filename, "wb");
  assert(fpr != NULL);
  fwrite(&npix, sizeof(long), (size_t)1, fpr);
  fwrite(dpix, sizeof(float), (size_t)npix, fpr);
  fwrite(&npix, sizeof(long), (size_t)1, fpr);
  fclose(fpr);
  fprintf(stderr, "Write completed.\n");
  free(dpix);

	return;
}

void cal_image_deltat()
{
	long ii, jj, kk;
	FILE *fpr = NULL;
	char filename[256];
	int ref;

	long npix=12L*(1L<<(2*Order));
	long NumGals;
	float *imgs;
	long *srt;
	float flux[Bins];
	long nimgs[Bins], snimgs[Bins];
	float ave[Bins], ave_all;
	float min, max;
	double theta, phi;
	long inpix, pix[4];
	double wgt[4];
	float *dpix_all, *dpix[Bins];

	fprintf(stderr, "Read images galaxies from file...\n");
	sprintf(filename, "%s/images.bin", path);
	fpr = fopen(filename, "rb");
  assert(fpr != NULL);
  ref = fread(&NumGals, sizeof(long), (size_t)1, fpr);
  imgs = (float *)malloc(NumGals*4*sizeof(float));
  assert(imgs != NULL);
  ref = fread(imgs, sizeof(float), (size_t)(NumGals*4), fpr);
  ref = fread(&NumGals, sizeof(long), (size_t)1, fpr);
  fclose(fpr);
	sprintf(filename, "%s/images_sorting.bin", path);
	fpr = fopen(filename, "rb");
	assert(fpr != NULL);
	ref = fread(&NumGals, sizeof(long), (size_t)1, fpr);
	srt = (long *)malloc(NumGals*sizeof(long));
	assert(srt != NULL);
	ref = fread(srt, sizeof(long), (size_t)NumGals, fpr);
	ref = fread(&NumGals, sizeof(long), (size_t)1, fpr);
  fclose(fpr);
	fprintf(stderr, "\n");

	fprintf(stderr, "The range of the images...\n");
	fprintf(stderr, "NumGals=\t%ld\n", NumGals);
  min = imgs[0]; max = imgs[0];
  for (ii = 1L; ii < NumGals; ii++) {
    if (imgs[ii*4] < min) min = imgs[ii*4];
    if (imgs[ii*4] > max) max = imgs[ii*4];
  }
  fprintf(stderr, "\tTheta:(%f ~ %f)\n", min, max);
  min = imgs[1]; max = imgs[1];
  for (ii = 1L; ii < NumGals; ii++) {
    if (imgs[ii*4+1] < min) min = imgs[ii*4+1];
    if (imgs[ii*4+1] > max) max = imgs[ii*4+1];
  }
  fprintf(stderr, "\tPhi  :(%f ~ %f)\n", min, max);
  min = imgs[2]; max = imgs[2];
  for (ii = 1L; ii < NumGals; ii++) {
    if (imgs[ii*4+2] < min) min = imgs[ii*4+2];
    if (imgs[ii*4+2] > max) max = imgs[ii*4+2];
  }
  fprintf(stderr, "\tFlux :(%.2e ~ %.2e)\n", min, max);
  min = imgs[3]; max = imgs[3];
  for (ii = 1L; ii < NumGals; ii++) {
    if (imgs[ii*4+3] < min) min = imgs[ii*4+3];
    if (imgs[ii*4+3] > max) max = imgs[ii*4+3];
  }
  fprintf(stderr, "\tMag  :(%.5f ~ %.5f)\n", min, max);
	fprintf(stderr, "\n");

	fprintf(stderr, "Divide images to different flux bins...\n");
	/*for (ii = 0; ii < Bins; ii++)
    nimgs[ii] = 0;
  for (ii = 0; ii < NumGals; ii++) {
    if (imgs[ii*4+2] > flux[0]) {
      for (jj = 1; jj < Bins; jj++) {
        if (imgs[ii*4+2] <= flux[jj]) {
          nimgs[jj-1]++;
          break;
        }
      }
			if (jj == Bins)
      	nimgs[Bins-1]++;
    }
  }*/
	for (ii = 0; ii < Bins; ii++) {
		nimgs[ii] = (long)(NumGals/Bins);
		if (ii == 0) snimgs[ii] = nimgs[ii];
		else snimgs[ii] = snimgs[ii-1] + nimgs[ii];
	}
	if (snimgs[Bins-1] != NumGals) {
		nimgs[Bins-1] += NumGals-snimgs[Bins-1];
		snimgs[Bins-1] = NumGals;
	}
	ave_all = 1.0F*NumGals/npix;
  for (ii = 0; ii < Bins; ii++)
    ave[ii] = 1.0F*nimgs[ii]/npix;
	/*for (ii = 0; ii < Bins-1; ii++)
    fprintf(stderr, "\tbin%d:(%-5.0f ~ %-5.0f),nimgs=%-9ld,ave=%4.3f\n", ii, flux[ii], flux[ii+1], nimgs[ii], ave[ii]);
  fprintf(stderr, "\tbin%d:(%-5.0f ~ max  ),nimgs=%-9ld,ave=%4.3f\n", ii, flux[ii], nimgs[ii], ave[ii]);*/
	for (ii = 0; ii < Bins; ii++) {
		if (ii == 0) min = imgs[srt[0]*4+2];
		else min = imgs[srt[snimgs[ii-1]]*4+2];
		max = imgs[srt[snimgs[ii]-1]*4+2];
		fprintf(stderr, "\tbin%d:(%-5.0f ~ %-5.0f),nimgs=%-9ld,ave=%4.3f\n", ii, min, max, nimgs[ii], ave[ii]);
	}
	fprintf(stderr, "\n");

  fprintf(stderr, "Write average overdensity to file...\n");
  sprintf(filename, "%s/ave_%ld.dat", path, Order);
  fpr = fopen(filename, "wb");
  fwrite(&Bins, sizeof(int), (size_t)1, fpr);
  fwrite(ave, sizeof(float), (size_t)Bins, fpr);
  fwrite(&Bins, sizeof(int), (size_t)1, fpr);
  fclose(fpr);
  fprintf(stderr, "\n");

	fprintf(stderr, "Distribute images to healpix points...\n");
  dpix_all = (float *)calloc(npix, sizeof(float));
  assert(dpix_all != NULL);
	for (ii = 0; ii < Bins; ii++) {
    dpix[ii] = (float *)calloc(npix, sizeof(float));
    assert(dpix[ii] != NULL);
  }
  /*for (ii = 0; ii < NumGals; ii++) {
    theta = (double)imgs[ii*4];
    phi = (double)imgs[ii*4+1];
    get_interpol(theta, phi, pix, wgt, Order);
    for (kk = 0; kk < 4; kk++)
      dpix_all[pix[kk]] += (float)wgt[kk];
		if (imgs[ii*4+2] > flux[0]) {
			for (jj = 1; jj < Bins; jj++) {
				if (imgs[ii*4+2] <= flux[jj])
					break;
			}
			for (kk = 0; kk < 4; kk++)
				dpix[jj-1][pix[kk]] += (float)wgt[kk];
		}
  }*/
	for (ii = 0; ii < NumGals; ii++) {
		theta = (double)imgs[srt[ii]*4];
		phi = (double)imgs[srt[ii]*4+1];
		inpix = ang2ring(theta, phi, Order);
		dpix_all[inpix] += 1.0F;
		//get_interpol(theta, phi, pix, wgt, Order);
    //for (kk = 0; kk < 4; kk++)
      //dpix_all[pix[kk]] += (float)wgt[kk];
		for (jj = 0; jj < Bins; jj++) {
			if (ii < snimgs[jj]) {
				dpix[jj][inpix] += 1.0F;
				//for (kk = 0; kk < 4; kk++)
					//dpix[jj][pix[kk]] += (float)wgt[kk];
				break;
			}
		}
	}
  for (ii = 0; ii < npix; ii++)
    dpix_all[ii] = (dpix_all[ii]-ave_all)/ave_all;
  min = dpix_all[0]; max = dpix_all[0];
  for (ii = 0; ii < npix; ii++) {
    if (dpix_all[ii] < min) min = dpix_all[ii];
    if (dpix_all[ii] > max) max = dpix_all[ii];
  }
  fprintf(stderr, "\tAll  dpix:(%.3f ~ %.3f)\n", min, max);
	for (ii = 0; ii < Bins; ii++) {
  	for (jj = 0; jj < npix; jj++)
    	dpix[ii][jj] = (dpix[ii][jj]-ave[ii])/ave[ii];
  	min = dpix[ii][0]; max = dpix[ii][0];
  	for (jj = 0; jj < npix; jj++) {
    	if (dpix[ii][jj] < min) min = dpix[ii][jj];
    	if (dpix[ii][jj] > max) max = dpix[ii][jj];
  	}
  	fprintf(stderr, "\tBin%d dpix:(%.3f ~ %.3f)\n", ii, min, max);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "Write healpix data to file...\n");
  sprintf(filename, "%s/deltat_%ld/h_image_deltat.dat", path, Order);
  fpr = fopen(filename, "wb");
  fwrite(&npix, sizeof(long), (size_t)1, fpr);
  fwrite(dpix_all, sizeof(float), (size_t)npix, fpr);
  fwrite(&npix, sizeof(long), (size_t)1, fpr);
  fclose(fpr);
  for (ii = 0; ii < Bins; ii++) {
    sprintf(filename, "%s/deltat_%ld/h_image_deltat%ld.dat", path, Order, ii);
    fpr = fopen(filename, "wb");
    fwrite(&npix, sizeof(long), (size_t)1, fpr);
    fwrite(dpix[ii], sizeof(float), (size_t)npix, fpr);
    fwrite(&npix, sizeof(long), (size_t)1, fpr);
    fclose(fpr);
  }
  fprintf(stderr, "\n");

	free(imgs);
  free(dpix_all);
  for (ii = 0; ii < Bins; ii++)
    free(dpix[ii]);

	fprintf(stderr, "All done!\n");
	return;
}
