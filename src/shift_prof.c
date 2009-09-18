#include <stdio.h>
/* 
   shift a profile located in the memory position *profile by an integer
   number of phase bins = shift. nbins is the total number of bins. 
*/
void shift_prof(float *profile, int nbins, int ishift) /* includefile */
{
  float *prcopy;
  int i,j;

  /* assign space and make a copy of the original profile */
  prcopy=(float *) malloc(nbins * sizeof(float));
  for (i=0;i<nbins;i++) prcopy[i]=profile[i];

  for (i=0;i<nbins;i++) {
    j = i + ishift;
    j = (j > nbins-1) ? j - nbins : j;
    j = (j < 0) ? j + nbins : j;
    profile[j]=prcopy[i];
  }
  free(prcopy);
}
