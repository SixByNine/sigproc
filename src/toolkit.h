/* ************************************************************** */
/* Routines to convert doubles to floats                          */
/* ************************************************************** */

void TKconvertFloat1(double *x,float *ox,int n);
void TKconvertFloat2(double *x,double *y,float *ox,float *oy,int n);


/* ************************************************************** */
/* Basic statistics                                               */
/* ************************************************************** */

float  TKfindMin_f(float *x,int n);
float  TKfindMedian_f(float *val,int count);
float  TKfindRMS_f(float *x,int n);
float  TKfindRMS_d(double *x,int n);
float  TKfindMax_f(float *x,int n);
float  TKmean_f(float *x,int n);

double TKvariance_d(double *x,int n);
double TKrange_d(double *x,int n);
double TKfindMin_d(double *x,int n);
double TKfindMax_d(double *x,int n);
double TKmean_d(double *x,int n);
double TKfindMin_d(double *x,int n);
double TKsign_d(double a,double b);
double TKretMax_d(double a,double b);
double TKretMin_d(double a,double b);
float TKretMax_f(float a,float b);
float TKretMin_f(float a,float b);
int TKretMin_i(int a,int b);

void meanrms(float* vals, int l, float* m ,float* r);

/* ************************************************************** */
/* Sorting                                                        */
/* ************************************************************** */

void TKsort_f(float *val,int nobs);
void TKsort_d(double *val,int nobs);
void TKsort_2f(float *val,float *val2,int nobs);
void quicksort_index(float* array, int* index, int npts);
void quicksort_inplace(float* array, int npts);



/* ************************************************************** */
/* Routines that modify a data array                              */
/* ************************************************************** */

void TKzeromean_d(int n,double *y);

/* ************************************************************** */
/* Random number routines                                         */
/* ************************************************************** */

double TKranDev(long *seed);
double TKgaussDev(long *seed);
long TKsetSeed();
void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
double genrand_real1(void);

