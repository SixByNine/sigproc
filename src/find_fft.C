
/*
  A fortran->C++ translation of find_fft.f
  M Bailes 13 Jan 2009

  dat is an array of the data to be ffted (floats)
  ndat is the dimension of the array (int)

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern "C" {
  int sglfft_(float *, float *, int *, int *,int *,int *);
  int realtr_(float*,float *, int *, int *);
};

void find_fft(int * ndat, float * dat){

  printf("ndat is %d\n",*ndat);
  int n2 = 0;

  while ((int)pow(2,n2)<*ndat)n2++;

  int nfft = (int)pow(2,n2-1);
  printf("nfft is %d\n",nfft);

  int nfft2 = nfft/2;
  int two = 2;

  sglfft_(dat, &dat[1],&nfft2,&nfft2,&nfft2,&two);
  realtr_(dat,&dat[1],&nfft2,&two);
  *ndat = nfft;
}

