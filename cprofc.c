#include <math.h>

#define NBINMAX 8192

void cprofc(float *prof, int nbins, float *amp, float *pha)
{

  int i,nh,n2;
  int forward = 1, back = -1, real = 1, complex = 0; 
  float temp[NBINMAX][2];

  nh = nbins/2;
  n2 = nbins*2;
  
  for(i=0;i<nh;i++) {
    temp[i][0] = prof[2*i];
    temp[i][1] = prof[2*i+1];
    temp[nh+i][0] = 0.;
    temp[nh+i][1] = 0.;  
  }

  ffft_(&temp[0][0],&nbins,&forward,&real);

  for(i=0;i<nh+1;i++) {
    amp[i] = sqrt(temp[i][0]*temp[i][0] + temp[i][1]*temp[i][1]);
    pha[i] = 0.;
    if(amp[i] > 0.) 
      pha[i] = atan2(temp[i][1],temp[i][0]); 
  }

}

void uncprofc(float *amp, float *pha, int nbins, float *c)
{

  int i,nh;
  
  nh = nbins/2;

  for(i=0;i<nh;i++) {
    c[2*i] = amp[i]*cos((double)pha[i]);
    c[2*i+1] = amp[i]*sin((double)pha[i]);
  }
  for(i=1;i<nh+1;i++) {  /* Important to set Nyquist component!!!  */
    c[2*nbins-2*i] = amp[i]*cos((double)pha[i]);
    c[2*nbins-2*i+1] = -amp[i]*sin((double)pha[i]);
  }

}
