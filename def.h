#define ORDER 12L
#define BIN 10
#define BIN1 11
#define BIN2 55

#define METHOD 3

extern int mtp;
extern float ave[BIN];
extern float slp[BIN];
extern float *cmb;
extern float *cab[BIN2];
extern float *bis1[BIN];
extern float *bis2[BIN];
extern float *wgt[BIN];
extern float *kalm;
extern float *kap1;
extern float *kap2;
extern char path[128];

void read_data();
float *read_cl(int dex1, int dex2);
float *read_alm(int dex);
float *make_power_spectrum(float *alm);
void cal_bias(int iter);
void cal_weight(int iter);
void cal_kappa(int iter);
void cal_kappa2();
void write_bias();
void write_kappa();
