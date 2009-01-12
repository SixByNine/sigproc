#include <stdio.h>
#include <math.h>
#define NMAX 8192  /* Need a 2-d array */
#define PI 3.1415926
#define TWOPI 6.2831852
/* 
   shift a profile located in the memory position *profile by the
   correct phase = fshift. nbins is the total number of bins. 
   Note that this is a frequency-domain rotation!
*/
void fshift_prof(float *profile, int nbins, double fshift) /* includefile */
{
  float profout[NMAX][2],*amp,*pha;
  float pha1;
  int i,j,nh;
  int forward=1,back=-1,complex=0,real=1;

  nh = nbins/2;
  amp=(float *) malloc((nh+1) * sizeof(float));
  pha=(float *) malloc((nh+1) * sizeof(float));

  cprofc(profile,nbins,amp,pha);
  pha1 = fshift*TWOPI;
  for(i=1;i<nh;i++)
    pha[i] = fmod((pha[i] + i*pha1),TWOPI);
  uncprofc(amp,pha,nbins,&profout[0][0]);
  for(i=0;i<nbins;i++)
    for(j=0;j<2;j++)
      profout[i][j] /= nbins;
  ffft_(&profout[0][0],&nbins,&back,&complex);
  for(i=0;i<nbins;i++)
    profile[i] = profout[i][0]/2.;

  free(amp);
  free(pha);
}
