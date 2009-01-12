#include <stdio.h>
#include "sigproc.h"
void submean(float *profile, int nbins) /*includefile*/
{
  float *prcopy, mean;
  int i,j;

  /* assign space and make a copy of the original profile and find peak bin */
  prcopy=(float *) malloc(nbins * sizeof(float));
  mean=0.0; 
  for (i=0;i<nbins;i++) {
    prcopy[i]=profile[i];
    if (prcopy[i]>mean) {
      j=i;
      mean=prcopy[i];
    }
  }

  /* shift copy to center based on peak bin: j */
  shift_prof(prcopy,nbins,j+1-nbins);

  /* finally calculate mean from off-pulse data */
  mean=0.0;
  for (i=0;i<nbins/20;i++) mean+=prcopy[i];
  for (i=nbins-1-nbins/20;i<nbins;i++) mean+=prcopy[i];
  mean/=nbins/10;

  /* now subtract this from original profile */
  for (i=0;i<nbins;i++) profile[i]-=mean;

  /* discard the copy */
  free(prcopy);
}
