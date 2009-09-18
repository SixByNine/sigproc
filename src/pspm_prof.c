#include <stdio.h>
#include "sigproc.h"
/* read in header and profiles from timing data returns a pointer to the
   profiles which are stacked in order of decreasing radio frequency */
float *pspm_prof(FILE *input, int nbins, int nchans, int *table)/*includefile*/
{
  int swap_bytes,i,j,k,l;
  unsigned long *rawdat;
  float *profile;
 
  /* PSPM/BPP data are big endian -> must swap bytes if using little endian */
  swap_bytes=little_endian();
 
  /* allocate space for raw and ordered data */
  rawdat = (unsigned long *) malloc(nchans * nbins * sizeof(long));
  profile = (float *) malloc(nbins * nchans * sizeof(float));

  /* read in the profiles */
  fread(rawdat,nbins*nchans*sizeof(unsigned long),1,input);

  /* stack profiles in frequency order, swapping if necessary */
  for (i=0;i<nchans;i++) {
    k=table[i]*nbins;
    l=i*nbins;
    for (j=0;j<nbins;j++) {
      if (swap_bytes) swap_long(&rawdat[k+j]);
      profile[l+j]=rawdat[k+j];
    }
  }
  free(rawdat);
  return(profile);

  for (i=0;i<nchans;i++) {
    puts("#START");
    k=i*nbins;
    for (j=0;j<nbins;j++) {
      printf("%f %f\n",(float)j,profile[k+j]);
    }
    puts("#STOP");
  }
  exit(0);
}
