#include <stdio.h>
/* 
   return a pointer to an array of filterbank channel frequencies given the
   center frequency fmid and the offset between each of the nchan channels
   IHS incorporated additional offset for WAPP data (wapp_off)
*/
double *chan_freqs(double fmid, double foff, int nchans, int wapp_off) /* includefile */
{
  int i; 
  double *chanfreq;
  chanfreq = (double *) malloc( nchans * sizeof(double) ); 
  for (i=0;i<nchans;i++)  {
    chanfreq[i]=fmid+(nchans/2-i)*foff-0.5*((double)wapp_off)*foff;
  }
  return (chanfreq);
}
