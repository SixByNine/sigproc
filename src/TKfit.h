void TKleastSquares_svd(double *x,double *y,double *sig,int n,double *p,double *e,int nf,double **cvm, double *chisq, void (*fitFuncs)(double, double [], int),int weight);
void TKremovePoly_f(float *px,float *py,int n,int m);
void TKremovePoly_d(double *px,double *py,int n,int m);
void TKleastSquares_svd_noErr(double *x,double *y,int n,double *p,int nf, void (*fitFuncs)(double, double [], int));     
void TKfitPoly(double x,double *v,int m);
void TKfitPolyOrigin(double x,double *v,int m);
