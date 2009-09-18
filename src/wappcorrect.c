#include <stdio.h>
#include <math.h>
/* subroutine to return correction to wapp_time (us) based on mjd */
double wappcorrect(double mjd) /*includefile*/
{
  double correction;

  /* assume no correction initially */
  correction=0.0;

  if ( (mjd >= 51829.0) && (mjd < 51834.0) ) correction=-0.08;
  if ( (mjd >= 51834.0) && (mjd < 51859.0) ) correction=-0.68;
  if ( (mjd >= 51859.0) && (mjd < 51969.0) ) correction=+0.04;

  if (correction != 0.0) {
    fprintf(stderr,"WARNING: correction %f us applied for MJD %.1f\n",
	    correction,mjd);
    fflush(stderr);
  }

  return(correction);
}
