#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "header.h"

/**
 * I have extracted the sub-zero functionality to a seperate function,
 * for neatness and style. M.Keith 2007
 * BASED ON ORIGIANAL WORK BY R.Eatough 2007
 */
void subzero(float* block, int nSamplesRead){

     int i,j,k;
     float sum;

     k=0;
     for(i=0; i<nSamplesRead; i++){
       sum=0;
       for(j=0; j<nchans; j++){
       sum=sum+block[k++]; 
       }
       sum=sum/(float)nchans;
       k=k-nchans;
       for(j=0; j<nchans; j++){
       block[k]=block[k]-sum;
       k=k+1;
       }
     }
}
