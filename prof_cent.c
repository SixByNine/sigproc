/* centres profiles in the entire dataset */
#include <stdio.h>
#include "sigproc.h"
void prof_cent(float *profile, int nbins, int nchans) /*includefile*/
{
  int i,j,k,l,norg;
  float *temp;
  void centprof();

  temp=vector(1,nbins);
  k=l=0;
  for (i=1;i<=nchans;i++) {
    for (j=1;j<=nbins;j++) {
      k++;
      temp[j]=profile[k];
    }
    centprof(temp,nbins);
    for (j=1;j<=nbins;j++) {
      l++;
      profile[l]=temp[j];
    }
  }
  free_vector(temp,1,nbins);
}

/* centres a profile -- ie shifts it in phase so that peak is nbins/2 */
void centprof(float *profile, int nbins)
{
  int i,j;
  float prmax;
  
  prmax=profile[1];
  j=1;
  for (i=2;i<=nbins;i++) {
    if (profile[i] > prmax) {
      prmax=profile[i];
      j=i;
    }
  }

  shift_prof(profile,nbins,nbins/2-j);
}

