#include "filterbank.h"
int obits;
#define GULP 32768*16
main (int argc, char *argv[])
{
  float block[GULP],copy[GULP];
  int ixnc,kxnc,i,j,k,n,off,nadd,nread,headersize=0;

  input=fopen(argv[1],"rb");
  output=fopen(argv[2],"wb");
  nadd=atoi(argv[3]);

  for (i=0;i<4;i++) ifstream[i]='Y';

  if (!(headersize=read_header(input))) 
    error_message("could not read header parameters!");
  obits=32;
  tsamp*=nadd;
  filterbank_header(output);
  for (i=0;i<GULP;i++) block[i]=copy[i]=0.0;
  nread=GULP/nchans;

  while(n=read_block(input,nbits,block,nread)) {
    k=off=0;
    for (i=0;i<n/nchans;i++) {
      ixnc=i*nchans;
      kxnc=k*nchans;
      for (j=0;j<nchans;j++) copy[off+j]+=block[ixnc+j];
      k++;
      if (k==nadd) {
	off+=nchans;
	k=0;
      }
    }
    for (j=0;j<n/nadd;j++) copy[j]/=(float)nadd;
    fwrite(copy,sizeof(float),n/nadd,output);
    for (j=0;j<n/nadd;j++) copy[j]=0.0;
  }
}
