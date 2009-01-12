#include "dedisperse.h"
int nbands,nobits;
double userdm;
main (int argc, char *argv[])
{
  float block[32768],copy[32768],median,median0;
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
    if (median0==0.0) median0=median;
    for (i=0;i<n;i++) copy[i]=(block[i]-median)/median0;
    fwrite(copy,sizeof(float),n,output);
  }

}
