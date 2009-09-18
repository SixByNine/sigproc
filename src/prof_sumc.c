/* collapses nchannels down to a number of sub-bands */
#include <stdio.h>
void prof_sumc(float *profile, int nbins, int nbands, int *nchans, int nifs, int *ignore) /* includefile */
{
  int i,c,n,b,cpb,nc,m;
  float *sub_band;

  /* calculate number of channels in each sub-band (integer multiples only!) */
  nc=*nchans;
  cpb=nc/nbands;
  if (nc != cpb*nbands) 
    error_message("nchannels not integer multiple of number of sub-bands!");

  /* create working array to hold one sub-band */
  sub_band = (float *) malloc(nbins * sizeof(float));
  for (n=0;n<nbins;n++) sub_band[n]=0.0;
  m=0;
  
  for (b=0;b<nbands;b++) {
    for (i=0;i<nifs;i++) {
      /* sum over all channels to form the sub-band */
      for (c=cpb*b; c<cpb*(b+1); c++) {
	if (!ignore[c])
	  for (n=0; n<nbins; n++) sub_band[n]+=profile[c*nbins*nifs+i*nbins+n];
      }
      /* copy subband to profile and clear for next sum */
      for (n=0;n<nbins;n++) {
	profile[m]=sub_band[n];
	sub_band[n]=0.0;
	m++;
      } 
    }
  }
  free(sub_band);
  *nchans=nbands;
}



