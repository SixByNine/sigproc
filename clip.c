#include "dedisperse.h"
int nbands,nobits;
double userdm;
main (int argc, char *argv[])
{
  float block[32768],copy[32768],mean,sigma,sum,ssq,mnsq,median;
  int i,n,headersize=0;
  input=fopen(argv[1],"rb");
  output=stdout;

  if (!(headersize=read_header(input))) 
    error_message("could not read header parameters!");
  userdm=refdm;
  nobits=nbits;
  nbands=1;
  dedisperse_header();

  while((n=read_block(input,nbits,block,32768))) {
    for (i=0;i<n;i++) copy[i]=block[i];
    median=nrselect((unsigned long)n/2,(unsigned long)n,copy);
    sum=ssq=0.0;
    for (i=0; i<n; i++) {
      sum+=block[i];
      ssq+=block[i]*block[i];
    }
    mean=sum/(float)n;
    mnsq=ssq/(float)n;
    sigma=sqrt(mnsq-mean*mean);
    for (i=0; i<n; i++) if (abs(block[i]-median)>sigma) block[i]=median;
    fwrite(block,sizeof(float),n,output);
  }

}
