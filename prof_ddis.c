/* dedispserses profiles with respect to arbitrary sky frequency (MHz) */
#include <stdio.h>
#include "sigproc.h"
int fftshift;
void prof_ddis(float *profile, int nbins, int nchans, int nbands, int nifs, double *chanfreq,  double period, double dm, double reference_frequency, float jyf1, float jyf2) /* includefile */
{
  int b,i,j,c,cpb,ishift;
  double freq,fshift;
  float jyf[2];
  cpb=nchans/nbands;
  for (b=0;b<nbands;b++) {
    for (c=b*cpb;c<(b+1)*cpb;c++) {
      /* choose reference frequency for dedispersion */
      if (reference_frequency<=0.0) 
	freq=chanfreq[b*cpb];
      else 
	freq=reference_frequency;
      /* calculate phase or bin shift for this channel */
      fshift=dmdelay(freq,chanfreq[c],dm)/period;
      ishift=dmdelay(freq,chanfreq[c],dm)/(period/nbins);
      /* check for wrapping */
      while (ishift < -nbins) ishift+=nbins;
      while (ishift > +nbins) ishift-=nbins;
      if (nbands>1)
	chanfreq[b]=chanfreq[b*cpb];
      else
	chanfreq[c]=freq;
      /* HACK! HACK! HACK!  Replace pol 0, bin 0 for Sept.1 2003 */
      /*      profile[0] = profile[1];  */
      jyf[0] = jyf1;
      jyf[1] = jyf2;
      for (i=0; i<nifs; i++) {
	if (fftshift) 
	  fshift_prof(profile,nbins,fshift); /* frequency-domain method */
	else
	  shift_prof(profile,nbins,ishift); /* bad old time-domain way! */
	for(j=0;j<nbins;j++)
          profile[j] *= jyf[i];
	profile+=nbins;
      }
    }
  }
}
