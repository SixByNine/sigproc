#include <math.h>

/*  Compute pulsar phase and frequency at time mjd0+mjd1/86400.  
    Note that it would really be better to save the polyco MJDs 
    as two doubles. */

/* Hardwiring max sets of coefficients to 16 - this may need to change. */
/* Hardwiring max coefficients to 12 - this may need to change. */

#define MAXSETS 16
#define MAXCOEFF 12

void phcalc(double mjd0,double mjd1,double *phase,double *psrfreq,double *rphase,double *psr_f0,double *poly_tmid,double *coeff,int *num_coeff) /*includefile*/
{
  double dtmin,mjd;
  int i,j;
  int show;
  double nblk;
  static int icurr;

  show = 0;

  nblk = (poly_tmid[1]-poly_tmid[0])*1440.;
  mjd = mjd0+mjd1/86400.;

  icurr = -1;
  *psrfreq = psr_f0[0];                     /* Default psrfreq */
  for (j=0;j<MAXSETS;j++) {
    dtmin = (mjd-poly_tmid[j])*1440.;  /* Time from center of this set*/
    if (fabs(dtmin) <= nblk/2.) {
      if(show) {
	printf("Set %1d, %1d coeffs works: %f %f\n",j,num_coeff[j],poly_tmid[j],poly_tmid[j+1]);
	for(i=0;i<num_coeff[j];i++) 
	  printf("%1d %e\n",i,coeff[j*num_coeff[j]+i]);
      }
      *psrfreq = 0.;                    /* Compute psrfreq and phase from */
      *phase = coeff[j*num_coeff[j]+num_coeff[j]-1];      /* the polynomial coeffs. */
      if (show) printf("phase = %21.15e   :%21.15e\n",*phase,coeff[j*num_coeff[j]+num_coeff[j]-1]);  
      for(i=num_coeff[j]-1;i>0;--i) {
        *psrfreq = dtmin*(*psrfreq) + i*coeff[j*num_coeff[j]+i];
        *phase = dtmin*(*phase) + coeff[j*num_coeff[j]+i-1];
	if (show) printf("phase = %21.15e   :%21.15e\n",*phase,coeff[j*num_coeff[j]+i-1]); 
      }

      *psrfreq = psr_f0[j] + *psrfreq/60.;  /* Add in the DC terms and scale */
      *phase += rphase[j]+dtmin*60.*psr_f0[j];
      if (show) printf("phase = %21.15e   f0: %21.15e\n",*phase,psr_f0[j]);
      *phase -= floor(*phase);
      if ((*phase < 0.) || (*phase > 1.))
        { printf("phase = %21.15f\n",*phase); exit(1); }
      icurr = j;
      break;
    }
  }
  if (icurr == -1) {
    printf("MJD %9.3f out of range (%9.3f to %9.3f)\n",
         (mjd),poly_tmid[0]-nblk/2880.,poly_tmid[MAXSETS-1]+nblk/2880.);
    *phase = -999.;
    printf("isets = %d\n",MAXSETS);
    exit(1);
  }

}
