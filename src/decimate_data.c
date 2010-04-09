#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "decimate.h"

void decimate_data(FILE *input, FILE *output) /*includefile*/
{ 
  char string[80];
  float *fblock,min,max,mid,realtime,fsaved[2];
  unsigned short *sblock;
  unsigned char  *cblock;
  int nsaved=0,ns,nsblk,opened=0,nout,i;

  nsblk=nchans*nifs*naddt;
  fblock=(float *) malloc(nsblk*sizeof(float));
  sblock=(unsigned short *) malloc(nsblk*sizeof(unsigned short));
  cblock=(unsigned char *) malloc(nsblk*sizeof(unsigned short));
  realtime=min=0.0;
  max=(float) pow(2.0,(double)obits) -1.0;

  mid=((float) pow(2.0,(double)nbits) -1.0)/2.0f;
  if(naddc > 1){
	  mid*=naddc;
  }

  /* main decimation loop */
  while ((ns=read_block(input,nbits,fblock,nsblk))>0) {
    add_channels(fblock,ns,naddc);
    add_samples(fblock,nifs,nchans/naddc,naddt);
    if (!opened) {
      /* open up logfile */
      open_log("decimate.monitor");
      opened=1;
    }
    nout=ns/naddt/naddc;
    switch (obits) {
    case 32:
      fwrite(fblock,sizeof(float),nout,output);
      break;
    case 16:
      float2short(fblock,nout,min,max,sblock);
      fwrite(sblock,sizeof(unsigned short),nout,output);
      break;
    case 8:
      float2char(fblock,nout,min,max,cblock);
      fwrite(cblock,sizeof(unsigned char),nout,output);
      break;
    case 4:
      if (nout==1) {
	/* must have at least two samples for four-bit packing save this one */
	fsaved[nsaved]=fblock[0];
	nsaved++;
	if (nsaved==2) {
	  /* we have 2 saved! write out */
	  float2four(fsaved,nsaved,min,max,cblock);
	  fwrite(cblock,sizeof(unsigned char),1,output);
	  nsaved=0;
	}
      } else {
	/* normal case */
	float2four(fblock,nout,min,max,cblock);
	fwrite(cblock,sizeof(unsigned char),nout/2,output);
      }
      break;
    case 2:
      float2two(fblock,nout,min,max,cblock);
      fwrite(cblock,sizeof(unsigned char),nout/4,output);
      break;
    case 1:
      for (i=0; i < nout; i++){
	      if(fblock[i] - mid > 0){
		      fblock[i]=1;
	      } else {
		      fblock[i]=0;
	      }
      }
      float2one(fblock,nout,0,1,cblock);
      fwrite(cblock,sizeof(unsigned char),nout/8,output);
      break;
    }
    realtime+=(float) tsamp * (float) ns/(float) nchans/(float) nifs;
    sprintf(string,"time:%.1fs",realtime);
    update_log(string);
  }
  update_log("finished");
  close_log();
}
