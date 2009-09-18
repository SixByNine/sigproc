#include <stdio.h>
#include <math.h>
float select(unsigned long, unsigned long, float *);
/*
 Routine to whiten a data stream of length npts by subtracting a running 
 median so that it has zero mean and then normalizing the data to unit rms.
 Used to flatten sloping baselines in amplitude spectra of noisy time series.
 Has the advantage over similar algorithms using running mean which can be
 biased by large spectral features. 
 
 dlorimer@atnf.csiro.au - 2001/10/12 - initial version for use in FIND
*/
void whiten_(float *data, int npts, int nsegments) 
{
  int i,j,k;
  unsigned long nwrk;
  float *work, *median, sumsq, rms;

  /* allocate space to store medians and data segments */
  nwrk=npts/nsegments;
  median=(float *) malloc(nsegments*sizeof(float));
  work=(float *) malloc(nwrk*sizeof(float));

  /* find the median of each segment */
  for (i=0; i<nsegments; i++) {
    k=0;
    for (j=i*nwrk; j<(i+1)*nwrk; j++) work[k++]=data[j];
    median[i]=select(nwrk/2,nwrk,work-1);
  }

  /* subtract running median and sum squares every npts/1024 samples */
  sumsq=0.0;  j=npts/1024;
  for (i=0; i<npts; i++) {
    if (data[i] != 0.0) data[i]-=median[i/nwrk];
    if (i%j == 0) sumsq+=data[i]*data[i];
  }

  /* normalize to unit rms */
  rms=sqrt(sumsq/1024.0);  for (i=0; i<npts; i++) data[i]/=rms;

  /* job done - free up space and return */
  free(median);  free(work);
}
