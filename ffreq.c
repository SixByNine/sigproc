#include<math.h>
double ffreq(double tsamp, int npf, int fold, int k) /*includefile*/
{
/*
  Returns the fluctuation frequency (Hz) of bin k in the amplitude
  spectrum having npf points. Tsamp is the sampling interval
  of the corresponding time domain data (seconds), whilst fold
  refers to 1 plus the number of harmonic sums that have produced
  the present spectrum. e.g fold=1 -> refers to the raw spectrum
  fold=5 refers to 4 harmonic sums (16 harmonics).
*/
  return((double)k/(2.0*tsamp*(double)npf*pow(2.0,(double)(fold-1))));
}
int fbin(double tsamp, int npf, int fold, double freq) /*includefile*/
/*
  Mathematical inverse of the above routine - returns bin number for
  given frequency
*/
{
  return((int)(freq*(2.0*tsamp)*(double)npf*pow(2.0,(double)(fold-1))));
}
