#include <stdio.h>
#include "sigproc.h"
void submedian(float *profile, int nbins) /*includefile*/
{
  float *prcopy, median;
  int i,j;

  /* assign space and make a copy of the original profile */
  prcopy=(float *) malloc(nbins * sizeof(float));
  for (i=0;i<nbins;i++) prcopy[i]=profile[i];

  /* find median (don't worry about even/odd numbers) */
  median=nrselect(nbins/2,nbins,prcopy-1);

  /* now subtract this from original profile */
  for (i=0;i<nbins;i++) profile[i]-=median;

  /* discard the copy */
  free(prcopy);
}

void subcal(float *profile, int nbins) /*includefile*/
{
  int i;
  float num,sum,sub;
  num=sum=0.0;
  for (i=nbins-nbins/3;i<nbins;i++) {
    sum+=profile[i];
    num+=1.0;
  }
  sub=sum/num;
  for (i=0;i<nbins;i++) profile[i]-=sub;
}
