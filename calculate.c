#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "def.h"

void make_njab(int ll, float ans[BIN1], float *nja, float *njb);

float *read_cl(int dex1, int dex2)
{
  int ll, ll2;

  char filename[128];
  FILE *pf;
  int ref;

  float tmp, *ps;

  if (dex1 == -1 || dex2 == -1)
    sprintf(filename, "%s/alm_cl/h_source_cl.dat", path);
  else
    sprintf(filename, "%s/alm_cl/h_image%d%d_cl.dat", path, dex1, dex2);
  pf = fopen(filename, "rb");
  assert(pf != NULL);
  ref = fread(&ll,sizeof(int),(size_t)1,pf);
  ps = (float *)malloc((ll-1)*sizeof(float));
  assert(ps != NULL);
  ref = fread(&tmp, sizeof(float), (size_t)1, pf);
  ref = fread(ps, sizeof(float), (size_t)(ll-1), pf);
  ref = fread(&ll2,sizeof(int),(size_t)1,pf);
  assert(ll2 == ll);
  fclose(pf);

  if (dex1 == -1 || dex2 == -1)
    mtp = ll-1;

  return ps;
}

float *read_alm(int dex)
{
  int Nalm;

  char filename[128];
  FILE *pf;
  int ref;

  int *index;
  float *alm;

  sprintf(filename, "%s/alm_cl/h_image%d%d_alm.dat", path, dex, dex);
  pf = fopen(filename, "rb");
  assert(pf != NULL);
  ref = fread(&Nalm,sizeof(int),(size_t)1,pf);
  assert(Nalm == (mtp+1)*(mtp+2)/2);
  index = (int *)malloc(Nalm*sizeof(int));
  alm = (float *)malloc(Nalm*2*sizeof(float));
  assert(index != NULL && alm != NULL);
  ref = fread(index,sizeof(int),(size_t)Nalm,pf);
  ref = fread(alm,sizeof(float),(size_t)Nalm*2,pf);
  ref = fread(&Nalm,sizeof(int),(size_t)1,pf);
  assert(Nalm == (mtp+1)*(mtp+2)/2);
  fclose(pf);

  assert(index[Nalm-1] == (mtp+1)*(mtp+1));
  free(index);

  return alm;
}

float *make_power_spectrum(float *alm)
{
  int mm, ll;
  int Nalm, pp;
  float *cl;
  float alm_[2];

  cl = (float *)malloc(mtp*sizeof(float));
  assert(cl != NULL);

  Nalm = (mtp+1)*(mtp+2)/2;

  pp = 1;
  for (ll = 1; ll <= mtp; ll++) {
    *(cl+ll-1) = 0.F;
    *(cl+ll-1) += (*(alm+pp))*(*(alm+pp))+(*(alm+Nalm+pp))*(*(alm+Nalm+pp));
    pp++;
    for (mm = 1; mm <= ll; mm++) {
      *(cl+ll-1) += (*(alm+pp))*(*(alm+pp))+(*(alm+Nalm+pp))*(*(alm+Nalm+pp));
      alm_[1] = *(alm+pp);
      alm_[2] = pow(-1, mm)*(*(alm+Nalm+pp));
      *(cl+ll-1) += alm_[1]*alm_[1]+alm_[2]*alm_[2];
      pp++;
    }
    *(cl+ll-1) /= (2*ll+1);
    *(cl+ll-1) *= ll*(ll+1)/(2*M_PI);
  }
  assert(pp == Nalm);
  return cl;
}

void cal_bias(int iter)
{
  int ii, jj, ll;
  int pp;
  int ind[10];

  char filename[1024];
  FILE *pf;

  int Nalm;
  float *alm;
  float *cij[BIN];

  float min, max, ave;

  fprintf(stderr, "Compute the galaxy bias...\n");

  pp = 0;
  for (ii = 0; ii < BIN; ii++) {
    ind[ii] = pp;
    for (jj = ii; jj < BIN; jj++)
      pp++;
  }

  if (iter == 0) {
    bis1[0] = (float *)malloc(mtp*BIN*sizeof(float));
    assert(bis1[0] != NULL);
    for (ii = 1; ii < BIN; ii++)
      bis1[ii] = bis1[0]+mtp*ii;
  }

  for (ii = 0; ii < BIN; ii++) {
    if (iter == 0)
      cij[ii] = cab[ind[ii]];
    else {
      alm = read_alm(ii);
      pp = 1;
      Nalm = (mtp+1)*(mtp+2)/2;
      for (ll = 1; ll <= mtp; ll++) {
	for (jj = 0; jj <= ll; jj++) {
	  *(alm+pp) -= slp[ii]*(*(kalm+pp));
	  *(alm+Nalm+pp) -= slp[ii]*(*(kalm+Nalm+pp));
	  pp++;
	}
      }
      cij[ii] = make_power_spectrum(alm);
      free(alm);
    }
  }

  for (ii = 0; ii < BIN; ii++) {
    for (ll = 1; ll <= mtp; ll++)
      *(bis1[ii]+ll-1) = sqrt((*(cij[ii]+ll-1))/(*(cmb+ll-1)));
  }

  for (ii = 0; ii < BIN; ii++) {
    ave = 0.F;
    min = *bis1[ii]; max = *bis1[ii];
    for (ll = 1; ll <= mtp; ll++) {
      ave += *(bis1[ii]+ll-1);
      if (*(bis1[ii]+ll-1) < min) min = *(bis1[ii]+ll-1);
      if (*(bis1[ii]+ll-1) > max) max = *(bis1[ii]+ll-1);
    }
    fprintf(stderr, "\t\t%f|%f|%f\n", min, ave/mtp, max);
  }

  if (iter > 0)
    for (ii = 0; ii < BIN; ii++)
      free(cij[ii]);
}

void cal_weight(int iter)
{
  int ii, ll;
  float Snb2, Snbg, Sng2;
  float lambda1, lambda2;

  fprintf(stderr, "Compute the weight...\n");

  if (iter == 0) {
    wgt[0] = (float *)malloc(mtp*BIN*sizeof(float));
    assert(wgt[0] != NULL);
    for (ii = 0; ii < BIN; ii++)
      wgt[ii] = wgt[0]+mtp*ii;
  }

  for (ll = 1; ll <= mtp; ll++) {
    Snb2 = 0.0; Snbg = 0.0; Sng2 = 0.0;
    for (ii = 0; ii < BIN; ii++) {
      Snb2 += ave[ii]*pow(*(bis1[ii]+ll-1),2.0);
      Snbg += ave[ii]*(*(bis1[ii]+ll-1))*slp[ii];
      Sng2 += ave[ii]*pow(slp[ii],2.0);
    }
    lambda1 = -2*Snb2/(pow(Snbg,2.0)-Snb2*Sng2);
    lambda2 =  2*Snbg/(pow(Snbg,2.0)-Snb2*Sng2);
    for (ii = 0; ii < BIN; ii++)
      *(wgt[ii]+ll-1) = (ave[ii]/2)*(lambda1*slp[ii]+lambda2*(*(bis1[ii]+ll-1)));
  }
}

void cal_kappa(int iter)
{
  int ii, mm, ll;
  int ind;

  int Nalm;
  float *alm;
  float min, max;

  fprintf(stderr, "Compute the kappa...\n");

  if (kalm != NULL) free(kalm);
  kalm = (float *)calloc((mtp+1)*(mtp+2), sizeof(float));
  assert(kalm != NULL);

  for (ii = 0; ii < BIN; ii++) {

    alm = read_alm(ii);

    ind = 1;
    Nalm = (mtp+1)*(mtp+2)/2;
    for (ll = 1; ll <= mtp; ll++) {
      for (mm = 0; mm <= ll; mm++) {
				*(kalm+ind) += *(wgt[ii]+ll-1)*(*(alm+ind));
				*(kalm+Nalm+ind) += *(wgt[ii]+ll-1)*(*(alm+Nalm+ind));
				ind++;
      }
    }
  }
  free(alm);

  if (kap1 != NULL) free(kap1);
  kap1 = make_power_spectrum(kalm);

  min = *kap1; max = *kap1;
  for (ll = 2; ll <= mtp; ll++) {
    if (*(kap1+ll-1) < min) min = *(kap1+ll-1);
    if (*(kap1+ll-1) > max) max = *(kap1+ll-1);
  }
  fprintf(stderr, "\t\t%.2e|%.2e\n", min, max);
  fprintf(stderr, "Done!\n");
}

void cal_kappa2()
{
  int ii, jj, ll;
  int round, ind;

  int *index;
  float *alm;

  float *nja, *njb;
  float delx[BIN1];
  float ans[BIN1];
  float err, lasterr;

  nja = (float *)calloc(BIN1*BIN1, sizeof(float));
  njb = (float *)calloc(BIN1, sizeof(float));
  assert(nja != NULL && njb != NULL);

  bis2[0] = (float *)calloc(mtp*BIN, sizeof(float));
  for (ii = 0; ii < BIN; ii++)
    bis2[ii] = bis2[0]+mtp*ii;
  kap2 = (float *)calloc(mtp, sizeof(float));
  assert(bis2[0] != NULL && kap2 != NULL);

  for (ll = 1; ll <= mtp; ll++) {

    round = 0; lasterr = 1.F;
    for (ii = 0; ii < BIN; ii++) // Initial value
      ans[ii] = sqrt(cmb[ll-1]);
    ans[BIN] = 1.0e-3*cmb[ll-1];

    do {
      if (round > 1)
				lasterr = err;

      make_njab(ll, ans, nja, njb);

      for (ii = 0; ii < BIN1; ii++)
       	delx[ii] = *(njb+ii) / *(nja+BIN1*ii+ii);
      for (ii = 0; ii < BIN1; ii++)
       	ans[ii] += delx[ii];
      ans[BIN] = fabs(ans[BIN]);
      err = 0.F;
      ind = 0;
      for (ii = 0; ii < BIN; ii++) {
       	for (jj = ii; jj < BIN; jj++) {
         	err += pow(*(cab[ind]+ll-1)-ans[ii]*ans[jj]-ans[BIN]*slp[ii]*slp[jj], 2.0);
	  			ind++;
				}
    	}

      round++;

    } while (err < lasterr || (err > 1.e-10F && round < 10));

    for (ii = 0; ii < BIN; ii++)
      *(bis2[ii]+ll-1) = ans[ii] - (err>lasterr ? delx[ii] : 0);
    *(kap2+ll-1) = ans[BIN] - (err>lasterr ? delx[BIN] : 0);
  }

  free(nja); free(njb);
}

void make_njab(int ll, float ans[BIN1], float *nja, float *njb)
{
  int ii, jj, kk;
  int ind;
  float *ja, *jaT, *jb;//ja[55][11], jaT[11][55], jb[55]
  float div[BIN1];

  ja = (float *)calloc(BIN1*BIN2, sizeof(float));
  jaT = (float *)calloc(BIN2*BIN1, sizeof(float));
  jb = (float *)calloc(BIN2, sizeof(float));

  ind = 0;
  for (ii = 0; ii < BIN; ii++) { // write ja and jb
    for (jj = ii; jj < BIN; jj++) {
      if (ii == jj) {
				*(ja+ind*BIN1+ii) = 2*ans[ii];
      }
      else {
				*(ja+ind*BIN1+ii) = ans[jj];
				*(ja+ind*BIN1+jj) = ans[ii];
      }
      *(ja+ind*BIN1+BIN) = slp[ii]*slp[jj];
      *(jb+ind) = -(ans[ii]*ans[jj]+ans[BIN]*slp[ii]*slp[jj]-(*(cab[ind]+ll-1)));
      ind++;
    }
  }
  assert(ind == BIN2);

  for (ii = 0; ii < BIN2; ii++)	// write jaT
    for (jj = 0; jj < BIN1; jj++)
      *(jaT+BIN2*jj+ii) = *(ja+BIN1*ii+jj);

  for (ii = 0; ii < BIN1; ii++) { // write nja=jaT*ja
    for (jj = 0; jj < BIN1; jj++) {
      *(nja+BIN1*ii+jj) = 0.F;
      for (kk = 0; kk < BIN2; kk++)
				*(nja+BIN1*ii+jj) += *(jaT+BIN2*ii+kk)*(*(ja+BIN1*kk+jj));
    }
  }

  for (ii = 0; ii < BIN1; ii++) { // write njb=jaT*jb
    *(njb+ii) = 0.F;
    for (jj = 0; jj < BIN2; jj++)
      *(njb+ii) += *(jaT+BIN2*ii+jj)*(*(jb+jj));
  }

  for (ii = 0; ii < BIN1; ii++) {
    for (jj = ii+1; jj < BIN1; jj++) {
      div[jj] = *(nja+BIN1*jj+ii) / *(nja+BIN1*ii+ii);
      *(nja+BIN1*jj+ii) = 0.F;
    }
    for (jj = ii+1; jj < BIN1; jj++)
      for (kk = ii+1; kk < BIN1; kk++)
        *(nja+BIN1*jj+kk) -= *(nja+BIN1*ii+kk)*div[jj];
    for (jj = ii+1; jj < BIN1; jj++)
      *(njb+jj) -= *(njb+ii)*div[jj];
  }

  for (ii = BIN; ii >= 0; ii--) {
    for (jj = 0; jj < ii; jj++) {
      div[jj] = *(nja+BIN1*jj+ii) / *(nja+BIN1*ii+ii);
      *(nja+BIN1*jj+ii) = 0.F;
      *(njb+jj) -= *(njb+ii)*div[jj];
    }
  }

  free(ja); free(jaT); free(jb);
}
