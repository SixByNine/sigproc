#include <math.h>
#include <stdio.h>

void find_formspec(int npts,float * data){

  // Just the amplitudes version

  //printf("formspec n = %d\n",npts);

  if (0){
  for (int i=0;i<npts/2-1;i++){
    float re=data[2*i];
    float im=data[2*i+1];
    data[i]=sqrt(re*re+im*im);
  }
  }
  // "Improved" bin-interbin version

  if (1){
  float re_l = data[0];
  float im_l = data[1];

  for (int i=1;i<npts/2;i++){
    float re=data[2*i];
    float im=data[2*i+1];

    float amp=sqrt(re*re+im*im);
    float amp_diff = sqrt(0.5*(re-re_l)*(re-re_l) + 0.5*(im-im_l)*(im-im_l));

    if (amp>amp_diff)
      data[i]=amp;
    else
      data[i]=amp_diff;
    re_l = re;
    im_l = im;
  }
 }
}
