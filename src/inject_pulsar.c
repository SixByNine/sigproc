#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sigproc.h"
#include "header.h"
#include "mjklog.h"
#include <math.h>
#include <time.h>
#include <tempo2pred.h>
#include "mjk_cmd.h"
#include "mjk_random.h"

float MX_val;
void write_block(int nobits,int nsout,int nbands, FILE* output,float* outblock);

int main (int argc, char** argv){
   FILE* input;
   FILE* output;
   char pred_fname[1024];
   char fil_fname[1024];
   char prof_fname[1024];
   int res,i;
   float *block;
   float *profile;
   float *grads;
   float spec_index=0;
   float in_snr=10; // S/N
   float ref_freq=1400.0;
   float frac;
   uint32_t seed;
   FILE* prof_file;
   T2Predictor pred;

   if(argc < 3){
	  fprintf(stderr,"%s [.fil] [t2pred.dat] [prof.asc]\n\n",argv[0]);
	  fprintf(stderr,"Pulsar insertion tool. M. Keith 2014.\n");
	  fprintf(stderr,"Use tempo2 -f [.par] -pred \"...\" to generate the predictor\n");
	  exit(0);
   }


   seed=time(NULL);

   in_snr=getF("--snr","-s",argc,argv,in_snr);
   ref_freq=getF("--freq","-f",argc,argv,1400.0);
   spec_index=getF("--sidx","-i",argc,argv,-1.5);
   seed=getI("--seed","-S",argc,argv,seed);


   input=fopen(argv[1],"r");
   strcpy(pred_fname,argv[2]);
   prof_file=fopen(argv[3],"r");

   output=stdout;

   res = T2Predictor_Read(&pred, pred_fname);
   if(res!=0){
	  fprintf(stderr,"error, could not read predictor %d\n",res);
   }
   int hdrsize = read_header(input);
   long long nsamp = nsamples(argv[1],hdrsize,nbits,nifs,nchans);
   

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
   grads=malloc(sizeof(float)*MX_prof);
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
   float sidx[nchans];
   for(ch = 0; ch < nchans; ch++){
	  freq[ch] = fch1 + foff*ch;
	  sidx[ch] = pow(freq[ch]/1400.0,spec_index);
   }

   float A;
   switch(nbits){
	  case 8:
		 MX_val=pow(2,8)-1;
		 A=24;
		 break;
	  case 4:
		 MX_val=pow(2,4)-1;
		 A=3;
		 break;
	  case 2:
		 MX_val=pow(2,2)-1;
		 A=1;
		 break;
	  case 1:
		 MX_val=pow(2,1)-1;
		 A=1;
		 break;
   }
   float sum=0;
   // normalise the profile
   for (i=0; i < nprof; i++)
	  sum+=profile[i];

   double scale=1.0;
   float psr_freq = T2Predictor_GetFrequency(&pred,mjd,freq[0]);
	
   sum/=nprof;
   scale = A*in_snr / sqrt(nchans*nsamp)/sum;
   logmsg("scale=%lg",scale);

   // normalise to pseudo-S/N in 1s.
   for (i=0; i < nprof; i++)
	  profile[i]=profile[i]*scale;

   for (i=0; i < nprof-1; i++){
	  grads[i] = profile[i+1]-profile[i];
   }
   grads[nprof-1]=profile[0]-profile[nprof-1];

   mjk_rand_t *rnd = mjk_rand_init(seed);

   logmsg("m=%f t0=%f dt=%f dmjd=%f",(float)mjd,(float)tstart,(float)tsamp,(float)tsamp_mjd);
   uint64_t count=0;
   while(!feof(input)){
	  read_block(input,nbits,block,nchans);

	  if(count%1024==0){
		 fprintf(stderr,"%ld samples,  %.1fs\r",count,count*tsamp);
	  }
//#pragma omp parallel for private(ch,phase,pbin,frac,A) shared (block,sidx,nchans,nprof) 
	  for(ch = 0; ch < nchans; ch++){
		 phase = T2Predictor_GetPhase(&pred,mjd,freq[ch]);
		 phase = phase-floor(phase);
		 pbin = floor(phase*nprof);
		 frac = phase*nprof-pbin;
		 //if(ch==0 && profile[pbin] > 0)logmsg("%Lf %d %f %f %f",phase,pbin,frac,profile[pbin],grads[pbin]);
		 A = profile[pbin] + frac*grads[pbin];
		 block[ch]+=sidx[ch]*A*(mjk_rand_double(rnd)); // this makes it probabilistic that it will pass 1 for fractional amplitudes
	  }
	  write_block(nbits,1,nchans,output,block);
	  mjd+=tsamp_mjd;
	  count++;
   }


   mjk_rand_free(rnd);
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
		 for (i=0; i<nsout*nifs*nbands; i++) {
			if (outblock[i] > MX_val) onebyte[i]=255;
			else onebyte[i] = ((unsigned char) outblock[i]);
		 }
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
