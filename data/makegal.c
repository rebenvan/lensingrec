#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fitsio.h>
#include "healpix_utils.h"

void read_source_galaxy(long *Numgals, float **Gals);
void read_image_galaxy(long Numgals, float *Gals, float **Imgs);
void readsource(long *Numgals, float **Gals);
void readimage(long Numgals, float *Gals, float **Imgs);
void sort_images(long Nimgs, float *imgs, long *imgs_dex);

int Plane;
char path[128];

int main(int argc, char *argv[])
{
	long Numgals;
	float *Gals, *Imgs;

	if (argc == 2) {
		Plane = atoi(argv[1]);
		sprintf(path, "/huawei/osv1/FAN/kappadata/sphere_%04d", Plane);
	}
	else {
		fprintf(stderr, "Please enter the Plane.\n");
		return 0;
	}

	fprintf(stderr, "\n\nFor source...\n");
	fprintf(stderr, "----------------------------------------\n\n");
	read_source_galaxy(&Numgals, &Gals);
	fprintf(stderr, "\n----------------------------------------\n");

	fprintf(stderr, "\n\nFor images...\n");
	fprintf(stderr, "----------------------------------------\n\n");
	read_image_galaxy(Numgals, Gals, &Imgs);
	fprintf(stderr, "\n----------------------------------------\n");

	free(Gals); free(Imgs);

	return 0;
}

void read_source_galaxy(long *Numgals, float **Gals)
{
	long ii;
	char filename[128];
	FILE *fpr;
	int ref;

	long Ngals;
	float *gals;
	float min, max;

	// ----- Read the source galaxy data from file -----
	fprintf(stderr, "Read the source galaxies...\n");

	readsource(Numgals, Gals);
	Ngals = *Numgals;
	gals = *Gals;

	sprintf(filename, "%s/galaxies.bin", path);
	fpr = fopen(filename, "wb");
	assert(fpr != NULL);
	fwrite(&Ngals, sizeof(long), (size_t)1, fpr);
	fwrite(gals, sizeof(float), (size_t)Ngals*5, fpr);
	fwrite(&Ngals, sizeof(long), (size_t)1, fpr);
	fclose(fpr);
	fprintf(stderr, "Read succefully!\n\n");

	// ----- Show the range of the galaxies -----
	fprintf(stderr, "The range of the galaxies...\n");
	min = *gals; max = *gals;
	for (ii = 1L; ii < Ngals; ii++) {
		if (*(gals+ii*5) < min) min = *(gals+ii*5);
		if (*(gals+ii*5) > max) max = *(gals+ii*5);
	}
	fprintf(stderr, "\tTheta:(%f, %f)\n", min, max);
	min = *(gals+1); max = *(gals+1);
	for (ii = 1L; ii < Ngals; ii++) {
		if (*(gals+ii*5+1) < min) min = *(gals+ii*5+1);
		if (*(gals+ii*5+1) > max) max = *(gals+ii*5+1);
	}
	fprintf(stderr, "\tPhi  :(%f, %f)\n", min, max);
	min = *(gals+2); max = *(gals+2);
	for (ii = 1L; ii < Ngals; ii++) {
		if (*(gals+ii*5+2) < min) min = *(gals+ii*5+2);
		if (*(gals+ii*5+2) > max) max = *(gals+ii*5+2);
	}
	fprintf(stderr, "\tLum  :(%.2e, %.2e)\n", min, max);
	min = *(gals+3); max = *(gals+3);
	for (ii = 1L; ii < Ngals; ii++) {
		if (*(gals+ii*5+3) < min) min = *(gals+ii*5+3);
		if (*(gals+ii*5+3) > max) max = *(gals+ii*5+3);
	}
	fprintf(stderr, "\tDist :(%.2e, %.2e)\n", min, max);
	min = *(gals+4); max = *(gals+4);
  for (ii = 1L; ii < Ngals; ii++) {
    if (*(gals+ii*5+4) < min) min = *(gals+ii*5+4);
    if (*(gals+ii*5+4) > max) max = *(gals+ii*5+4);
  }
  fprintf(stderr, "\tMass :(%.2e, %.2e)\n\n", min, max);

	gals = NULL;
	return;
}

void read_image_galaxy(long Numgals, float *Gals, float **Imgs)
{
	long ii;
	char filename[128];
	FILE *fpr;
	int ref;

	long Nimgs;
	float *imgs;
	long *imgs_dex;
	float min, max;

	// ---- Read the image galaxy from file -----
	fprintf(stderr, "Read the image galaxies...\n");

	readimage(Numgals, Gals, Imgs);
	Nimgs = Numgals;
	imgs = *Imgs;

	sprintf(filename, "%s/images.bin", path);
	fpr = fopen(filename, "wb");
	assert(fpr != NULL);
	fwrite(&Nimgs, sizeof(long), (size_t)1, fpr);
	fwrite(imgs, sizeof(float), (size_t)Nimgs*4, fpr);
	fwrite(&Nimgs, sizeof(long), (size_t)1, fpr);
	fclose(fpr);
	fprintf(stderr, "Read succefully!\n\n");

	// ----- Show the range of image galaxies, kappa and flux -----
	fprintf(stderr, "The range of the images...\n");
	min = imgs[0]; max = imgs[0];
	for (ii = 1L; ii < Nimgs; ii++) {
		if (imgs[ii*4] < min) min = imgs[ii*4];
		if (imgs[ii*4] > max) max = imgs[ii*4];
	}
	fprintf(stderr, "\tTheta:(%f, %f)\n", min, max);
	min = imgs[1]; max = imgs[1];
	for (ii = 1L; ii < Nimgs; ii++) {
		if (imgs[ii*4+1] < min) min = imgs[ii*4+1];
		if (imgs[ii*4+1] > max) max = imgs[ii*4+1];
	}
	fprintf(stderr, "\tPhi  :(%f, %f)\n", min, max);
	min = imgs[2]; max = imgs[2];
	for (ii = 1L; ii < Nimgs; ii++) {
		if (imgs[ii*4+2] < min) min = imgs[ii*4+2];
		if (imgs[ii*4+2] > max) max = imgs[ii*4+2];
	}
	fprintf(stderr, "\tFlux :(%.2e, %.2e)\n", min, max);
	min = imgs[3]; max = imgs[3];
	for (ii = 1L; ii < Nimgs; ii++) {
		if (imgs[ii*4+3] < min) min = imgs[ii*4+3];
		if (imgs[ii*4+3] > max) max = imgs[ii*4+3];
	}
	fprintf(stderr, "\tMag  :(%f, %f)\n\n", min, max);

	// ----- Sorting the images in flux -----
	fprintf(stderr, "Sorting all the images with flux...\n");
	imgs_dex = (long *)malloc(Nimgs*sizeof(long));
	assert(imgs_dex != NULL);

	sort_images(Nimgs, imgs, imgs_dex);

	sprintf(filename, "%s/images_sorting.bin", path);
	fpr = fopen(filename, "wb");
	assert(fpr != NULL);
	fwrite(&Nimgs, sizeof(long), (size_t)1, fpr);
	fwrite(imgs_dex, sizeof(long), (size_t)Nimgs, fpr);
	fwrite(&Nimgs, sizeof(long), (size_t)1, fpr);
	fclose(fpr);
	fprintf(stderr, "Sort Completed!\n\n");

	imgs = NULL;
	free(imgs_dex);
	return;
}

void readsource(long *Numgals, float **Gals)
{
	char filename[256];
	FILE *fpr;
	int ref;
	typedef struct {
		long index;
		int type;
		float pos[3];
		float totalmass;
		float stellarmass;
		float luminosity;
		float axialRatio;
		float posAngle;
		float nvec[3];
	} lcGalaxyDemo;
	lcGalaxyDemo *sourgals;
	long TotNumGals;
	double vec[3], theta, phi;
	long ii;

	sprintf(filename, "/huawei/osv1/wcl/RAYTRACING/L500/lensplanes/sourceGalaxy/sourceData/sourceGalaxy%04d.bin", Plane);
	fpr = fopen(filename, "rb");
	assert(fpr != NULL);
	ref = fread(&TotNumGals, sizeof(long), (size_t)1, fpr);
	fprintf(stderr, "\tTotNumGals=%ld\n", TotNumGals);
	sourgals = (lcGalaxyDemo *)malloc(TotNumGals*sizeof(lcGalaxyDemo));
	assert(sourgals != NULL);
	ref = fread(sourgals, sizeof(lcGalaxyDemo), (size_t)TotNumGals, fpr);
	ref = fread(&TotNumGals, sizeof(long), (size_t)1, fpr);
	fclose(fpr);

	*Gals = (float *)malloc(TotNumGals*5*sizeof(float));
	assert(*Gals != NULL);
	for (ii = 0; ii < TotNumGals; ii++) {
		vec[0] = (double)(sourgals+ii)->pos[0];
		vec[1] = (double)(sourgals+ii)->pos[1];
		vec[2] = (double)(sourgals+ii)->pos[2];
		vec2ang(vec, &theta, &phi);
		*(*Gals+ii*5) = (float)theta;	// 0 to PI
		*(*Gals+ii*5+1) = (float)phi;	// 0 to 2*PI
		*(*Gals+ii*5+2) = (float)1.0e10*((sourgals+ii)->luminosity);	//lum of sun
		*(*Gals+ii*5+3) = (float)sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);	// Mpc
		*(*Gals+ii*5+4) = (float)1.0e8*(sourgals+ii)->totalmass;	// mass of sun
	}

	free(sourgals);
	*Numgals = TotNumGals;

	return;
}

void readimage(long Numgals, float *Gals, float **Imgs)
{
	fitsfile *fpr;
	FILE *fp;
	long firstrow, firstelem;
	float floatnull;
	int anynull;
	long TotNumGals, BaseNumGalsInChunk;
	int round;
	int tfields = 7;
	int colnum;
	char filename[256];
	char *ttype[] = {"index","ra","dec","A00","A01","A10","A11"};
	int hdutype;
	int nkeys, keypos;
	char card[FLEN_CARD];
	long *index;
	double *ra, *dec, *A00, *A01, *A10, *A11;
	int status;

	double theta, phi, lum, flux;
	double mu, mag;
	double ds = 1.0/206265.0;	// distance from sun to earth(pc)
	double pc = 3.08567758e18; // distance unit(cm)
	double Ls = 3.846e33; // flux(erg/s)
	double ms = -26.74;	// magnitude of sun(v)
	double jy = 1.0e-23; // unit of Jansky(erg/s/cm2/Hz)
	double fr = 1.420406e9;	// frequency of the luminosity(Hz)
	long ii;

	status = 0;
	sprintf(filename, "/huawei/osv1/wcl/RAYTRACING/L500/outputs/galImage%04d.0000.fit", Plane);
	fits_open_file(&fpr, filename, READONLY, &status);
	if (status)
		fits_report_error(stderr, status);
	fits_movabs_hdu(fpr, BINARY_TBL, &hdutype, &status);
	if (status)
		fits_report_error(stderr, status);
	fits_get_hdrpos(fpr, &nkeys, &keypos, &status);
	if (status)
		fits_report_error(stderr, status);
	//fprintf(stderr, "\tHeader, listing for HDU #%d:\n", BINARY_TBL);
	for (ii = 0; ii < nkeys; ii++) {
		fits_read_record(fpr, ii, card, &status);
		if (status)
			fits_report_error(stderr, status);
		//fprintf(stderr, "\t%s\n", card);
	}

	fits_get_num_rows(fpr, &TotNumGals, &status);
	if (status)
		fits_report_error(stderr, status);
	fits_get_rowsize(fpr, &BaseNumGalsInChunk, &status);
	if (status)
		fits_report_error(stderr, status);
	round = TotNumGals/BaseNumGalsInChunk;
	if (round*BaseNumGalsInChunk != TotNumGals)
		round++;

	fprintf(stderr, "\tTotNumGals=%ld\n", TotNumGals);
	if (TotNumGals != Numgals) {
		fprintf(stderr, "Total source galaxies:%ld\n", Numgals);
		fprintf(stderr, "Total image galaxies: %ld\n", TotNumGals);
		TotNumGals = (TotNumGals<Numgals)?TotNumGals:Numgals;
	}
	index = (long *)malloc(TotNumGals*sizeof(long));
	ra = (double *)malloc(TotNumGals*sizeof(double));
	dec = (double *)malloc(TotNumGals*sizeof(double));
	A00 = (double *)malloc(TotNumGals*sizeof(double));
	A01 = (double *)malloc(TotNumGals*sizeof(double));
	A10 = (double *)malloc(TotNumGals*sizeof(double));
	A11 = (double *)malloc(TotNumGals*sizeof(double));
	assert(index != NULL && ra != NULL && dec != NULL);
	assert(A00 != NULL && A01 != NULL && A10 != NULL && A11 != NULL);

	firstelem = 1;
	firstrow = 1;
	for (ii = 1; ii <= round; ii++) {
		if (ii == round)
			BaseNumGalsInChunk = TotNumGals-firstrow+1;

		fits_get_colnum(fpr, CASEINSEN, ttype[0], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TLONG, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, index+firstrow-1, &anynull, &status);
		if (status)
			fits_report_error(stderr, status);

		fits_get_colnum(fpr, CASEINSEN, ttype[1], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TDOUBLE, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, ra+firstrow-1, &anynull, &status);
		if (status)
			fits_report_error(stderr, status);

		fits_get_colnum(fpr, CASEINSEN, ttype[2], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TDOUBLE, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, dec+firstrow-1, &anynull, &status);
		if (status)
			fits_report_error(stderr, status);

		fits_get_colnum(fpr, CASEINSEN, ttype[3], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TDOUBLE, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, A00+firstrow-1, &anynull, &status);
		if (status)
			fits_report_error(stderr, status);

		fits_get_colnum(fpr, CASEINSEN, ttype[4], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TDOUBLE, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, A01+firstrow-1, &anynull, &status);
		if (status)
			fits_report_error(stderr, status);

		fits_get_colnum(fpr, CASEINSEN, ttype[5], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TDOUBLE, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, A10+firstrow-1, &anynull, &status); 
		if (status)
			fits_report_error(stderr, status);

		fits_get_colnum(fpr, CASEINSEN, ttype[6], &colnum, &status);
		if (status)
			fits_report_error(stderr, status);
		fits_read_col(fpr, TDOUBLE, colnum, firstrow, firstelem, BaseNumGalsInChunk, &floatnull, A11+firstrow-1, &anynull, &status); 
		if (status)
			fits_report_error(stderr, status);

		firstrow += BaseNumGalsInChunk;
	}

	fits_close_file(fpr, &status);
	if (status)
		fits_report_error(stderr, status);

	*Imgs = (float *)malloc(TotNumGals*4*sizeof(float));
	assert(*Imgs != NULL);
	for (ii = 0L; ii < TotNumGals; ii++) {
		radec2ang(&theta, &phi, ra[ii], dec[ii]);
		*(*Imgs+ii*4) = (float)theta;
		*(*Imgs+ii*4+1) = (float)phi;
		mu = fabs(1.0/(A00[ii]*A11[ii]-A01[ii]*A10[ii]));
		flux = 1.0e6*Gals[index[ii]*5+2]*mu*Ls/(4.0*M_PI*pow(Gals[index[ii]*5+3]*1.0e6*pc,2))/(jy*fr);
		*(*Imgs+ii*4+2) = (float)flux;
		mag = ms+(-2.5)*log10(Gals[index[ii]*5+2]/pow(1.0e6*Gals[index[ii]*5+3],2)/(1.0/pow(ds,2)));
		*(*Imgs+ii*4+3) = (float)mag;
	}

	free(index); free(ra); free(dec);
	free(A00); free(A01); free(A10); free(A11);
	return;
}

void sort_images(long Nimgs, float *imgs, long *imgs_dex)
{
	long ii, Need;

	struct farray {
    struct farray *last;
    struct farray *father;
    long index;
    struct farray *next;
  };
  struct farray *first, *this, *new, *old;
	char filename[256];
	FILE *fpr;
  time_t mtime;

	fprintf(stderr, "Sorting images...\n");
  mtime = -time(NULL);
  first = (struct farray *)malloc(sizeof(struct farray));
  assert(first != NULL);
  first->index = 0L; first->father = NULL; first->last = NULL; first->next = NULL;
  for (ii = 1L; ii < Nimgs; ii++) { // Sorting the images to farray
    this = first;
    for ( ; ; ) {
      if (*(imgs+ii*4+2) >= *(imgs+(this->index)*4+2)) {
        if (this->next == NULL) {
          new = (struct farray *)malloc(sizeof(struct farray));
          assert(new != NULL);
          new->father = this;
          this->next = new;
          this = new;
          this->index = ii;
          this->last = NULL;
          this->next = NULL;
          break;
        }
        else {
          this = this->next;
        }
      }
      else {
        if (this->last == NULL) {
          new = (struct farray *)malloc(sizeof(struct farray));
          assert(new != NULL);
          new->father = this;
          this->last = new;
          this = new;
          this->index = ii;
          this->last = NULL;
          this->next = NULL;
          break;
        }
        else {
          this = this->last;
        }
      }
    }
    if ((ii+1)%100000 == 0)
      fprintf(stderr, "\r%ld   \t...\t   %.2f%%", ii+1, 100.F*(ii+1)/Nimgs);
  }
  fprintf(stderr, "\r%ld   \t...\t   100.00%%\n", Nimgs);
  mtime += time(NULL);
  fprintf(stderr, "Sorting completed!\n\tTime used:%lds\n\n", mtime);

	mtime = -time(NULL);
  Need = 0;
  this = first;
  do { // Read the index of sorted images to imgs_dex
    if (this->last != NULL) {
      this = this->last;
    }
    else {
      *(imgs_dex+Need) = this->index;
      Need++;
      if (Need%100000 == 0)
        fprintf(stderr, "\r%ld   \t...\t   %.2f%%", Need, 100.F*Need/Nimgs);
      if (this->next != NULL) {
        this = this->next;
      }
      else {
        for ( ; ; ) {
          if (this->father == NULL)
            break;
          if (this->father->next == this) {
            old = this;
            this = this->father;
            this->next = NULL;
            free(old);
          }
          else {
            old = this;
            this = this->father;
            this->last = NULL;
            free(old);
            break;
          }
        }
      }
    }
  } while (this != NULL && Need < Nimgs);
  fprintf(stderr, "\r%ld   \t...\t   100.00%%\n", Nimgs);
  free(first);
  mtime += time(NULL);
  fprintf(stderr, "Read sorting completed!\n\tTime used:%lds\n", mtime);

	return;
}
