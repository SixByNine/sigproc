/* scale profile so each bin occupies an integer value between 0 and 65535 */
#include <stdio.h>
#include "sigproc.h"
void scale_prof(float *profile, int nbins, unsigned long *iprofile, float *scale, float *offset) /*includefile*/
{
  int j;
  float denominator;
  *offset=vmin(profile,nbins);
  denominator=(vmax(profile,nbins)-*offset);
  if (denominator!=0.0)
    *scale=65535.0/denominator;
  else
    *scale=0.0;
  for (j=0;j<nbins;j++) {
    iprofile[j]=(profile[j]-(*offset))*(*scale);
  }
}
