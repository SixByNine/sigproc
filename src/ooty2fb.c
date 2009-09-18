/* ooty2fb - converts OOTY search-mode data into "filterbank" data */
#include "filterbank.h"
#define BUFSZ 4096
void ooty2fb(FILE *input, FILE *output) /* includefile */
{
  unsigned char *gulp;
  float realtime,*sample,sum;
  char string[80];
  int offset,s,i,j,k,doit,idump=0,nread,opened=0;
  gulp = (unsigned char *) malloc(BUFSZ * sizeof(unsigned char));
  sample = (float *) malloc(8*BUFSZ * sizeof(float));
  while ( (nread=fread(gulp,1,BUFSZ,input)) != 0) {
    for (i=0; i<nread; i++) gulp[i] = ~ gulp[i];
    realtime=tsamp*idump;
    if ((doit=process(realtime,start_time,final_time))==-1) break;
    if (doit) {
      switch (obits) {
      case 32:
	for (s=0; s<nread/32; s++) {
	  offset=s*32;
	  k=0;
	  sum=0.0;
	  for (i=0; i<32; i++) {
	    for (j=1;j<=8;j++) {
	      sample[k]=gulp[offset+i]&1;
	      gulp[offset+i]>>=1;
	      if (swapout) swap_float(&sample[k]);
	      k++;
	    }
	  }
	  fwrite(sample,sizeof(float),256,output);
	  idump++;
	}
	break;
      case 1:
	fwrite(gulp,1,nread,output);
	idump+=nread/32;
	break;
      }

      /* open up logfile if need be */
      if (idump%1024 == 0) {
	if (!opened) {
	  open_log("filterbank.monitor");
	  opened=1;
	}
	sprintf(string,"time:%.1fs",realtime);
	update_log(string);
      }
    }

  }
  free(gulp);free(sample);
}
