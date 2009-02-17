/* read a block of data into datablock */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigproc.h"
int remove_mean,nchans;
int read_block(FILE *input, int nbits, float *block, int nread) /*includefile*/
{
  int i,j,k,s1,s2,iread,s,c,offset;
  float mean;
  unsigned char *charblock;
  unsigned short *shortblock;
  long seed=0;

  /* decide how to read the data based on the number of bits per sample */
  switch(nbits) {
  case 1:
    /* read n/8 bytes into character block containing n 1-bit pairs */
    charblock=(unsigned char *) malloc(nread/8);
    iread=fread(charblock,1,nread/8,input);
    k=0;
    /* unpack 1-bit pairs into the datablock */
    for (i=0; i<iread; i++) {
      for (j=0;j<8;j++) {
	block[k++]=charblock[i]&1;
	charblock[i]>>=1;
      }
    }
    iread=k; /* this is the number of samples read in */
    free(charblock);
    break;
  case 4:
    /* read n/2 bytes into character block containing n 4-bit pairs */
    charblock=(unsigned char *) malloc(nread/2);
    iread=fread(charblock,1,nread/2,input);
    j=0;
    /* unpack 4-bit pairs into the datablock */
    for (i=0; i<iread; i++) {
      char2ints(charblock[i],&s1,&s2);
      block[j++]=(float) s1;
      block[j++]=(float) s2;
    }
    iread*=2; /* this is the number of samples read in */
    free(charblock);
    break;
  case 8:
    /* read n bytes into character block containing n 1-byte numbers */
    charblock=(unsigned char *) malloc(nread);
    iread=fread(charblock,1,nread,input);
    /* copy numbers into datablock */
    for (i=0; i<iread; i++) {
      block[i]=(float) charblock[i];
    }
    free(charblock);
    break;
  case 16:
    /* read 2*n bytes into short block containing n 2-byte numbers */
    shortblock=(unsigned short *) malloc(2*nread);
    iread=fread(shortblock,2,nread,input);
    /* copy numbers into datablock */
    for (i=0; i<iread; i++) {
      block[i]=(float) shortblock[i];
    }
    free(shortblock);
    break;
  case 32:
    /* read 4*n bytes into floating block containing n 4-byte numbers */
    iread=fread(block,4,nread,input); 
    break;
  default:
    error_message("read_block - nbits can only be 4,8,16 or 32!");
  }

  if (!remove_mean) return iread;
  /* subtract the mean of the channels from each individual channel */
  for (s=0; s<iread/nchans; s++) {
	mean=0.0;
	offset=s*nchans;
	for (c=0; c<nchans; c++) mean+=block[offset+c];
	mean/=(float)nchans;
	for (c=0; c<nchans; c++) block[offset+c]-=mean;
  }	
  /* return number of samples read */
  return iread;
}
