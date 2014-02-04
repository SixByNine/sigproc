#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sigproc.h"
#include "header.h"
#include <math.h>
#include <time.h>
#include <tempo2pred.h>

float MX_val;
void write_block(int nobits,int nsout,int nbands, FILE* output,float* outblock);

int main (int argc, char** argv){
   FILE* input;
   FILE* output;
   char pred_fname[1024];
   int res;
   float *block;
   float *profile;
   FILE* prof_file;
   T2Predictor pred;

   if(argc < 3){
	  fprintf(stderr,"%s [.fil] [t2pred.dat] [prof.asc]\n\n",argv[0]);
	  fprintf(stderr,"Pulsar insertion tool. M. Keith 2014.\n");
	  fprintf(stderr,"Use tempo2 -f [.par] -pred \"...\" to generate the predictor\n");
	  exit(0);
   }

   srand(time(NULL));

   input=fopen(argv[1],"r");
   strcpy(pred_fname,argv[2]);
   prof_file=fopen(argv[3],"r");

   output=stdout;

   res = T2Predictor_Read(&pred, pred_fname);
   if(res!=0){
	  fprintf(stderr,"error, could not read predictor %d\n",res);
   }
   int hdrsize = read_header(input);

   if(!hdrsize){
	  fprintf(stderr,"error, could not read sigproc header\n");
	  exit(1);
   }

   fclose(input);
   input=fopen(argv[1],"r");
   char* buf = malloc(sizeof(char)*hdrsize);
   fread(buf,1,hdrsize,input);
   fwrite(buf,1,hdrsize,output);

   int nprof=0;
   int MX_prof=512;
   profile=malloc(sizeof(float)*MX_prof);
   while(!feof(prof_file)){
	  if(nprof >= MX_prof){
		 MX_prof *= 2;
		 profile=realloc(profile,sizeof(float)*MX_prof);
	  }
	  fscanf(prof_file,"%f\n",profile+nprof);
	  nprof++;
   }

   block=malloc(sizeof(float)*nchans);
   long double mjd = (long double)tstart;
   long double tsamp_mjd = (long double)tsamp/86400.0L;
   long double phase;
   int pbin;
   int ch=0;
   float freq[nchans];
   for(ch = 0; ch < nchans; ch++){
	  freq[ch] = fch1 + foff*ch;
   }

   float A;
   switch(nbits){
	  case 16:
		 MX_val=pow(2,16)-1;
		 break;
	  case 8:
		 MX_val=pow(2,8)-1;
		 break;
	  case 4:
		 MX_val=pow(2,4)-1;
		 break;
	  case 2:
		 MX_val=pow(2,2)-1;
		 break;
	  case 1:
		 MX_val=pow(2,1)-1;
		 break;
   }
   fprintf(stderr,"m=%f t0=%f dt=%f dmjd=%f\n",(float)mjd,(float)tstart,(float)tsamp,(float)tsamp_mjd);
   while(!feof(input)){
	  read_block(input,nbits,block,nchans);
	  for(ch = 0; ch < nchans; ch++){
		 phase = T2Predictor_GetPhase(&pred,mjd,freq[ch]);
		 phase = phase-floor(phase);
		 pbin = floor(phase*nprof);
		 A = profile[pbin];
		 block[ch]+=A*(rand()/(float)RAND_MAX); // this makes it probabilistic that it will pass 1 for fractional amplitudes
	  }
	  write_block(nbits,1,nchans,output,block);
	  mjd+=tsamp_mjd;
   }


   return 0;
}


void write_block(int nobits,int nsout,int nbands, FILE* output,float* outblock){
   char *onebyte;
   short *twobyte;
   int nifs=1;
   int i=0;

   /* now write out samples and bat on */
   switch (nobits) {
	  case 1:
		 onebyte = (char *) malloc(nsout*nifs*nbands);
		 float2one(outblock,nsout*nifs*nbands,0,MX_val,onebyte);
		 fwrite(onebyte,sizeof(char),nsout*nifs*nbands/8,output);
		 break;

	  case 2:
		 onebyte = (char *) malloc(nsout*nifs*nbands);
		 float2two(outblock,nsout*nifs*nbands,0,MX_val,onebyte);
		 fwrite(onebyte,sizeof(char),nsout*nifs*nbands/4,output);
		 break;
	  case 4:
		 onebyte = (char *) malloc(nsout*nifs*nbands);
		 float2four(outblock,nsout*nifs*nbands,0,MX_val,onebyte);
		 fwrite(onebyte,sizeof(char),nsout*nifs*nbands/2,output);
		 break;
	  case 8:
#if SIGNED
			onebyte = (char *) malloc(nsout*nifs*nbands);
			for (i=0; i<nsout*nifs*nbands; i++) 
			   onebyte[i] = ((char) outblock[i])-128;
			fwrite(onebyte,sizeof(char),nsout*nifs*nbands,output);
#else
			onebyte = (unsigned char *) malloc(nsout*nifs*nbands);
			for (i=0; i<nsout*nifs*nbands; i++) 
			   onebyte[i] = ((unsigned char) outblock[i]);
			fwrite(onebyte,sizeof(unsigned char),nsout*nifs*nbands,output);
#endif
			break;
	  case 16:
			twobyte = (short *) malloc(nsout*nifs*nbands);
			for (i=0; i<nsout*nifs*nbands; i++) 
			   twobyte[i] = (short) outblock[i];
			fwrite(twobyte,sizeof(short),nsout*nifs*nbands,output);
			break;
	  case 32:
			fwrite(outblock,sizeof(float),nsout*nifs*nbands,output);
			break;
	  default:
			error_message("requested output number of bits can only be 8, 16 or 32");
			break;
   }



}
