/*
  FIND - a program to find periodic signals in a noisy time series. 
*/
#include "find.h"

main(int argc, char *argv[])
{
  FILE *fileptr;
  int i,numsamps,nfft,nspc,headersize,binmin;
  float *data, *amplitude, ar, ai, arl, ail, a1, a2;
  double frequency, pmax=9.99999;

  fileptr=open_file(argv[1],"rb");

  if (!(headersize=read_header(fileptr))) {
    error_message("could not read header parameters!");
    exit(1);
  }

  numsamps=nsamples(argv[1],headersize,nbits,nifs,nchans);
  nfft=(int)pow(2.0,(double)np2(numsamps));
  nspc=nfft/2;
  data=(float *) malloc(sizeof(float)*numsamps);
  
  printf("#Reading in %d samples from file: %s...\n",numsamps,argv[1]);
  if ((read_block(fileptr,nbits,data,numsamps)) != numsamps) {
    error_message("error reading data");
  }

  printf("#Performing 2^%d point FFT...\n",np2(numsamps));
  realft(data-1,(unsigned long) nfft,1);

  printf("#Forming amplitude spectrum...\n");
  amplitude=(float *) malloc(nspc*sizeof(float));
  binmin=fbin(tsamp,nspc,1,1.0/pmax);
  ar=ai=arl=ail=a1=a2=0.0;
  for (i=0;i<nspc;i++) {
    if (i<binmin) 
      amplitude[i]=0.0;
    else {
      ar=data[2*i];
      ai=data[2*i+1];
      a1=ar*ar+ai*ai;
      a2=0.5*((ar-arl)*(ar-arl)+(ai-ail)*(ai-ail));
      amplitude[i]=(a1>a2)?a1:a2;
      arl=ar;
      ail=ai;
    }
    /*printf("%f %f\n",ffreq(tsamp,nspc,1,i+1),amplitude[i]);*/
  }
  free(data);

}
