/* sums IFs for the entire data set, sets nifs to 1 */
#include <stdio.h>
#include "sigproc.h"
void prof_sumifs(float *profile, int nbins, int nchans, int *nifs) /*includefile*/
{
  int i,j,k;
  float *tmpprof;

  tmpprof = (float *) malloc(nchans * nbins * (*nifs) * sizeof(float));

  for(i=0;i<nchans*nbins*(*nifs);i++) tmpprof[i] = profile[i];
  k=0;
  for(i=0;i<nchans;i++) 
    for(j=0;j<nbins;j++) {
      profile[k]=tmpprof[i*nbins*(*nifs)+j]+tmpprof[i*nbins*(*nifs)+nbins+j];
      k++;
    }
  free(tmpprof);
  *nifs = 1;

}
