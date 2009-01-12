#include <stdio.h>
#include "sigproc.h"
void prof_adds(float *profile, int *nbins, int nchans, int nifs, int nadd)/*includefile*/
{
  int i,n=0;
  float sum=0.0;
  for (i=0;i<*nbins*nchans*nifs;i++) {
    sum+=profile[i];
    if ((i+1)%nadd==0) {
      profile[n]=sum/(float)nadd;
      sum=0.0;
      n++;
    }
  }
  *nbins/=nadd;
}

