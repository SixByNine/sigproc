/* subtracts a baseline from the entire data set */
#include <stdio.h>
#include "sigproc.h"
void prof_sbas(char *srcname, float *profile, int nbins, int nchans, int nifs) /*includefile*/
{
  int i,j;
  for (i=0;i<nchans*nifs;i++) {
    if (strings_equal(srcname,"cal")) {
      subcal(profile,nbins);
    } else {
      submedian(profile,nbins);
    }
    profile+=nbins;
  }
}
