#include <complex.h>
#include <fftw3.h>
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
void write_block(int_fast32_t nobits,int nsout,int nbands, FILE* output,float* outblock);
struct convolve_plan {
   fftwf_plan fwda;
   fftwf_plan fwdb;
   fftwf_plan bwd;
   float* a;
   float* b;
   float _Complex* c;
   float _Complex* d;
   float* output;
   int_fast32_t npts;
   float scale;
};
struct convolve_plan *setup_convolve(int_fast32_t npts, float* a, float* b,float* output);
void convolve(struct convolve_plan *plan);
void free_convolve_plan(struct convolve_plan *plan);

double sinc(double x){
   if(x==0)return 1.0;
   else return sin(x)/x;
}

int main (int_fast32_t argc, char** argv){
   FILE* input;
   FILE* output;
   char pred_fname[1024];
   char fil_fname[1024];
   char prof_fname[1024];
   char subprof_fname[1024];
   int_fast32_t res,i;
   float *block;
   float *profile;
   float spec_index=0;
   float in_snr=10; // S/N
   float ref_freq=1400.0;
   float t_scat=0;
   float scat_idx=0;
   float frac;
   uint32_t seed;
   uint_fast32_t nsubpulse=5;
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
   t_scat=getF("--scatter-time","-c",argc,argv,0.0);
   scat_idx=getF("--scatter-index","-C",argc,argv,4.0);
   spec_index=getF("--sidx","-i",argc,argv,-1.5);
   seed=getI("--seed","-S",argc,argv,seed);
   //seed=getI("--dmsmear","-d",argc,argv,seed);
   strcpy(subprof_fname,getS("--subprof","-b",argc,argv,""));


   input=fopen(argv[1],"r");
   strcpy(pred_fname,argv[2]);
   prof_file=fopen(argv[3],"r");

   output=stdout;

   res = T2Predictor_Read(&pred, pred_fname);
   if(res!=0){
	  fprintf(stderr,"error, could not read predictor %d\n",res);
   }
   int_fast32_t hdrsize = read_header(input);
   const uint_fast32_t nchan_const = nchans;
   const uint64_t nsamp = nsamples(argv[1],hdrsize,nbits,nifs,nchan_const);


   if(!hdrsize){
	  fprintf(stderr,"error, could not read sigproc header\n");
	  exit(1);
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

   fclose(input);
   input=fopen(argv[1],"r");
   char* buf = malloc(sizeof(char)*hdrsize);
   fread(buf,1,hdrsize,input);
   fwrite(buf,1,hdrsize,output);

   uint_fast32_t n=0;
   while(!feof(prof_file)){
	  fscanf(prof_file,"%f\n",&frac); // dummy value
	  n++;
   }
   const uint_fast32_t nprof=n;

   profile=fftwf_alloc_real(nprof);
   float *ism_conv = fftwf_alloc_real(nprof);
   float *subpulse_map = fftwf_alloc_real(nprof);
   float *subpulse_profile = fftwf_alloc_real(nprof);
   float *unsmeared_prof = fftwf_alloc_real(nprof);
   float *smeared_prof = fftwf_alloc_real(nprof);
   logmsg("Initialising convolution plans...");
   struct convolve_plan *ism_conv_plan = setup_convolve(nprof,unsmeared_prof,ism_conv,smeared_prof);
   struct convolve_plan *subpulse_conv_plan = setup_convolve(nprof,subpulse_profile,subpulse_map,unsmeared_prof);

   logmsg("%x",ism_conv_plan);
   logmsg("Done.");


  fseek(prof_file,0,SEEK_SET);
   n=0;
   while(!feof(prof_file)){
	  fscanf(prof_file,"%f\n",profile+n);
	  subpulse_profile[n]=1.0;
	  n++;
   }
   fclose(prof_file);
   if(strlen(subprof_fname)){
	  prof_file=fopen(subprof_fname,"r");
	  n=0;
	  while(!feof(prof_file)){
		 fscanf(prof_file,"%f\n",subpulse_profile+n);
		 n++;
	  }
	  fclose(prof_file);
   }


   double sum=0;
   // normalise the profile
   for (i=0; i < nprof; i++)
	  sum+=profile[i];
   double scale=1.0;

   sum/=(double)nprof;
   scale = A*in_snr / sqrt(nchan_const*nsamp)/sum;
   logmsg("mean=%lg scale=%lg",sum,scale);
// normalise to pseudo-S/N in 1s.
   for (i=0; i < nprof; i++)
	  profile[i]=profile[i]*scale;


   sum=0;
   // normalise the subpulse profile
   for (i=0; i < nprof; i++)
	  sum+=subpulse_profile[i];
   double subpulse_scale;
   sum/=(double)nprof;
   // we also normalise by the number of subpulses so that the final subpulse profile will have an area of 1.
   scale = 1.0/sum/(float)nsubpulse;

   for (i=0; i < nprof; i++)
	  subpulse_profile[i]=subpulse_profile[i]*scale;


   block=malloc(sizeof(float)*nchan_const);
   long double mjd = (long double)tstart;
   long double tsamp_mjd = (long double)tsamp/86400.0L;
   long double phase;
   int_fast32_t pbin;
   uint_fast32_t ch=0;
   float freq[nchan_const];
   float sidx[nchan_const];
   float dmdly[nchan_const];
   uint_fast32_t ism_idx[nchan_const];
   float ism_models[nchan_const][nprof];
   uint_fast32_t nism=0;
   const uint_fast32_t NCOPY = nprof*sizeof(float);
   const double psr_freq = (double)T2Predictor_GetFrequency(&pred,mjd,freq[0]);
   logmsg("Pulsar period ~%.2lfms",1000.0/psr_freq);

   for(ch = 0; ch < nchan_const; ch++){
	  freq[ch] = fch1 + foff*ch;
	  sidx[ch] = pow(freq[ch]/ref_freq,spec_index);
   }

   // set up the ISM model
   {
	  int_fast32_t dm_bin=-1;
	  int_fast32_t scatter_bin=-1;
	  int_fast32_t sbin;
	  float* dm_conv = ism_conv; // re-use an existing convolution function
	  float* scatt_conv = unsmeared_prof;
	  float* out_conv = smeared_prof;

	  const float tscat0 = (t_scat*pow(fch1/ref_freq,-scat_idx));
	  logmsg("Tscatter@%.1fMHz = %gs",ref_freq,t_scat);
	  logmsg("Tscatter@%.1fMHz = %gs",fch1,tscat0);
	  for(ch = 0; ch < nchan_const; ch++){
		 phase = T2Predictor_GetPhase(&pred,mjd,freq[ch]) - T2Predictor_GetPhase(&pred,mjd,freq[ch]+foff);
		 pbin = fabs(phase)*nprof;
		 if (pbin<1)pbin=1;
		 if (pbin%2==0)pbin-=1;

		 sbin = (int_fast32_t)floor(psr_freq*(t_scat*pow(freq[ch]/ref_freq,-scat_idx))*nprof);

		 if (pbin!=dm_bin || sbin != scatter_bin){
			// make new ISM model
			logmsg("New ISM Model...");
			logmsg("DM smearing. df=%lg dt=%lLg phase bins=%d",foff,phase,pbin);
			sum=0;
			for(i=0;i<nprof;i++){
			   int_fast32_t x = (i-pbin/2+nprof)%nprof;
			   double y = 2.*M_PI*(i-pbin/2.0+0.5)/pbin;
			   if(i<pbin){
				  //logmsg("%d %g %g",x,y,sinc(y));
				  dm_conv[x] = sinc(y);
				  sum+=dm_conv[x];
			   }
			   else dm_conv[x]=0;
			}
			for(i=0;i<nprof;i++){
			   //logmsg("%g %g",dm_conv[i],sum);
			   dm_conv[i]/=sum;
			}


			logmsg("Scattering. Exponential with phase bins=%d idx=%.1f",sbin,scat_idx);
			if(sbin > 0){
			   sum=0;
			   for(i=0;i<nprof;i++){
				  double y = i/(double)sbin;
				  //logmsg("%d %g %g",i,y,exp(-y));
				  scatt_conv[i] = exp(-y);
				  sum+=scatt_conv[i];
			   }
			   for(i=0;i<nprof;i++){
				  //logmsg("%g %g",ism_conv[i],sum);
				  scatt_conv[i]/=sum;
			   }
			} else {
			   scatt_conv[0]=1.0;
			   for(i=1;i<nprof;i++){
				  scatt_conv[i]=0;
			   }
			}

			convolve(ism_conv_plan); // dm_conv*scatt_conv => out_conv

			memcpy(ism_models[nism],out_conv,NCOPY);
			nism++;
			dm_bin=pbin;
			scatter_bin=sbin;
		 }
		 ism_idx[ch]=nism-1;
	  }
   }

   float  grads[nism][nprof];
   float output_prof[nism][nprof];
   double poff[nchan_const];
   double p0 = T2Predictor_GetPhase(&pred,mjd,freq[0]);
   int_fast32_t prevbin[nchan_const];
   int_fast32_t dbin,ii;
   mjk_rand_t **rnd = malloc(sizeof(mjk_rand_t*)*nchan_const);
   for(ch = 0; ch < nchan_const; ch++){
	  rnd[ch] = mjk_rand_init(seed+ch);
	  poff[ch] = T2Predictor_GetPhase(&pred,mjd,freq[ch])-p0;
	  prevbin[ch]=-1;
   }

   logmsg("m=%f t0=%f dt=%f dmjd=%f",(float)mjd,(float)tstart,(float)tsamp,(float)tsamp_mjd);
   uint64_t count=0;
   int64_t ip0;
   int64_t prev_ip0=INT64_MAX;
   while(!feof(input)){
	  if(count%1024==0){
		 fprintf(stderr,"%ld samples,  %.1fs\r",count,count*tsamp);
	  }

	  read_block(input,nbits,block,nchan_const);
	  p0 = T2Predictor_GetPhase(&pred,mjd,freq[0]);
	  ip0 = (int64_t)floor(p0);
	  if (ip0!=prev_ip0){
		 //logmsg("new pulse");
		 // need a new temp profile.

		 for(i=0; i < nprof;i++){
			subpulse_map[i]=0;
		 }
		 for(n=0; n < nsubpulse;n++){
			i=floor(mjk_rand_double(rnd[0])*nprof);
			subpulse_map[i]+=exp(mjk_rand_gauss(rnd[0]))/1.6487212707;
			//subpulse_map[i]+=(mjk_rand_double(rnd[0]))*2;
		 }
		 convolve(subpulse_conv_plan);
		 for(i=0; i < nprof;i++){
			unsmeared_prof[i]*=profile[i];
		 }
		 for(i =0; i<nism; i++){
			memcpy(ism_conv,ism_models[i],NCOPY);
			convolve(ism_conv_plan);
			memcpy(output_prof[i],smeared_prof,NCOPY);
			for (n=0; n < nprof-1; n++){
			   grads[i][n] = smeared_prof[n+1]-smeared_prof[n];
			}
			grads[i][nprof-1]=smeared_prof[0]-smeared_prof[nprof-1];
		 }
		 prev_ip0=ip0;
	  }
	  for(ch = 0; ch < nchan_const; ch++){
		 phase = p0+poff[ch];
		 phase = phase-floor(phase);
		 pbin = floor(phase*nprof);
		 frac = phase*nprof-pbin;
		 A = output_prof[ism_idx[ch]][pbin] + frac*grads[ism_idx[ch]][pbin];
		 if(prevbin[ch]<0)prevbin[ch]=pbin;
		 dbin=pbin-prevbin[ch];
		 while(dbin<0)dbin+=nprof;
		 while(dbin>1){
			A += output_prof[ism_idx[ch]][prevbin[ch]];
			prevbin[ch]+=1;
			dbin=pbin-prevbin[ch];
			while(dbin<0)dbin+=nprof;
		 }
		 if (A>0){
			block[ch] += sidx[ch]*A + mjk_rand_double(rnd[0]);
		 }
		 prevbin[ch]=pbin;
	  }
	  write_block(nbits,1,nchan_const,output,block);
	  mjd+=tsamp_mjd;
	  count++;
   }


   for(ch = 0; ch < nchan_const; ch++){
	  mjk_rand_free(rnd[ch]);
   }
   free(rnd);
   fftwf_free(profile);

   return 0;
}


void write_block(int_fast32_t nobits,int nsout,int nbands, FILE* output,float* outblock){
   char *onebyte;
   int16_t *twobyte;
   int_fast32_t nifs=1;
   int_fast32_t i=0;

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
		 twobyte = (int16_t *) malloc(nsout*nifs*nbands);
		 for (i=0; i<nsout*nifs*nbands; i++) 
			twobyte[i] = (int16_t) outblock[i];
		 fwrite(twobyte,sizeof(int16_t),nsout*nifs*nbands,output);
		 break;
	  case 32:
		 fwrite(outblock,sizeof(float),nsout*nifs*nbands,output);
		 break;
	  default:
		 error_message("requested output number of bits can only be 8, 16 or 32");
		 break;
   }
}
struct convolve_plan *setup_convolve(int_fast32_t npts, float* a, float* b, float* output){
   struct convolve_plan *plan = calloc(1,sizeof(struct convolve_plan));
   plan->a=a;
   plan->b=b;
   plan->c = fftwf_alloc_complex(npts);
   plan->d = fftwf_alloc_complex(npts);
   plan->output = output;
   plan->fwda = fftwf_plan_dft_r2c_1d(npts,plan->a,plan->c,FFTW_MEASURE);
   plan->fwdb = fftwf_plan_dft_r2c_1d(npts,plan->b,plan->d,FFTW_MEASURE);
   plan->bwd  = fftwf_plan_dft_c2r_1d(npts,plan->c,plan->output,FFTW_MEASURE);
   plan->npts=npts;
   plan->scale=(float)npts;
   return plan;
}
void convolve(struct convolve_plan *plan){
   int_fast32_t i;
   fftwf_execute(plan->fwda);
   fftwf_execute(plan->fwdb);
   for(i=0; i < plan->npts;i++){
	  plan->c[i] *= plan->d[i]/plan->scale;
   }
   fftwf_execute(plan->bwd);
}
void free_convolve_plan(struct convolve_plan *plan);


