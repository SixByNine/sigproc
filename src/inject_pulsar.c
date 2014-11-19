#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


// standard headers
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

// external libraries
#include <fftw3.h>
#include <tempo2pred.h>

// sigproc/mjk headers
#include "sigproc.h"
#include "header.h"
#include "mjklog.h"
#include "mjk_cmd.h"
#include "mjk_random.h"

#ifndef fftwf_alloc_real
#define fftwf_alloc_real(A) (float*)fftwf_malloc(A*sizeof(float))
#endif

#ifndef fftwf_alloc_complex
#define fftwf_alloc_complex(A) (fftwf_complex*)fftwf_malloc(A*sizeof(fftwf_complex*))
#endif

#define M_PI        3.14159265358979323846264338327950288


void print_help(){
   fprintf(stderr,\
		 "\ninject_pulsar [options] --pred t2pred.dat --prof prof.asc file.fil > output.fil\n"\
		 "\n"\
		 "Inject a simulated pulsar into the input .fil file. Writes to stdout.\n"\
		 "Michael Keith (2014) - mkeith@pulsarastronomy.net\n"\
		 "\n"\
		 "INPUTS:\n"\
		 "   t2pred.dat - A tempo2 predictor file (generate with tempo2 -f x.par -pred ...\n"\
		 "   prof.asc   - Pulse profile in single-column text format. 2^n bins is best\n"\
		 "   file.fil   - Input filterbank file. Can be a real file, or from fast_fake\n"\
		 "\n"\
		 "OPTIONS:\n"\
		 "   --help, -h          This help text\n"\
		 "   --snr,-s            Target signal-to-noise ratio (phase average S/N). (def=15)\n"\
		 "   --subprof,-b        Profile for sub-profile structure. Same format as prof.asc\n"\
		 "   --nsub,-n           Number of sub-pulses per profile, over full pulse phase (def=5).\n"\
		 "   --sidx,-i           Spectral index of pulsar. (def=-1.5)\n"\
		 "   --scatter-time,-c   Scattering timescale at ref freq, s. (def=no scattering).\n"\
		 "   --scint-bw,-C       Scintilation bandwidth, MHz. Cannot use in conjunction with -c\n"\
		 "   --scatter-index,-X  Index of scattering. (def=4.0).\n"\
		 "   --freq,-f           Reference frequency for scattering/spectral index, MHz. (def=1500)\n"\
		 "   --pulse-sigma,-E    'sigma' for log-normal pulse intensity distribution. (def=0.2)\n"\
		 "   --seed,-S           Random seed for simulation. (def=time())\n"\
		 "\n"\
		 "Note that the scint bandwidth is derived from the scatter time, so only one can be specified.\n"\
		 "\n");

}


float MX_val;
void write_block(int_fast32_t nobits,int nsout,int nbands, FILE* output,float* outblock);
struct convolve_plan {
   fftwf_plan fwda;
   fftwf_plan fwdb;
   fftwf_plan bwd;
   float* a;
   float* b;
   float complex* c;
   float complex* d;
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

int main (int argc, char** argv){
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
   float in_snr=15; // S/N
   float ref_freq=1400.0;
   float t_scat=0;
   float scint_bw;
   float scat_idx=0;
   float pulse_energy_sigma;
   float frac;
   double sum;
   uint32_t seed;
   uint_fast32_t nsubpulse;
   uint_fast32_t n;
   FILE* prof_file;
   T2Predictor pred;
   mjk_clock_t *CLK_setup = init_clock();
   mjk_clock_t *CLK_process = init_clock();
   mjk_clock_t *CLK_inner = init_clock();
   bool help=false;
   float A;

   // This is just a timer for the runtime computation.
   start_clock(CLK_setup);

   // set up fftw to run in multi-threaded mode.
   fftwf_init_threads();
   //fftwf_plan_with_nthreads(omp_get_max_threads());
   fftwf_plan_with_nthreads(1);

   // read command line parameters
   in_snr=getF("--snr","-s",argc,argv,in_snr);
   ref_freq=getF("--freq","-f",argc,argv,1400.0);
   t_scat=getF("--scatter-time","-c",argc,argv,-1);
   scint_bw=getF("--scint-bw","-C",argc,argv,-1);
   scat_idx=getF("--scatter-index","-X",argc,argv,4.0);
   spec_index=getF("--sidx","-i",argc,argv,-1.5);
   seed=getI("--seed","-S",argc,argv,time(NULL));
   strcpy(subprof_fname,getS("--subprof","-b",argc,argv,""));
   strcpy(prof_fname,getS("--prof","-p",argc,argv,"prof.asc"));
   strcpy(pred_fname,getS("--pred","-P",argc,argv,"t2pred.dat"));
   nsubpulse=getI("--nsub","-n",argc,argv,5);
   pulse_energy_sigma=getF("--pulse-sigma","-E",argc,argv,0.2);
   help=getB("--help","-h",argc,argv,0);
   getArgs(&argc,argv);

   
   if (argc!=2)help=1;
   if (help){
	  print_help();
	  exit(1);
   }


   if(scint_bw > 0 && t_scat>0){
	  logmsg("Note: ignoring scint_bw, using scattering time instead");
	  scint_bw=-1;
   }
   if(scint_bw < 0 && t_scat<0){
	  scint_bw=0;
	  t_scat=0;
	  logmsg("Disable scattering/scintillation");
   }
   if(scint_bw < 0){
	  scint_bw = 1e-6/(2*M_PI*t_scat);
   }

   if(t_scat < 0){
	  t_scat = 1e-6/(2*M_PI*scint_bw);
   }

   // double-check by printing values back to user.


   logmsg("input .fil file \t= '%s'",argv[1]);
   logmsg("Predicor     \t= '%s'",pred_fname);
   logmsg("Profile      \t= '%s'",prof_fname);
   if(strlen(subprof_fname)>0){
	  logmsg("SubProfile   \t= '%s'",subprof_fname);
	  logmsg("Nsubpulse    \t= %d",nsubpulse);
	  logmsg("Pulse energy sig\t= %.1g",pulse_energy_sigma);
   }
   logmsg("Target S/N   \t= %.2f",in_snr);
   logmsg("Ref Freq     \t= %.2f MHz",ref_freq);
   logmsg("Spectrum     \t= f^(%.2f)",spec_index);
   logmsg("t_scatter    \t= (%.2g s) * f^(%.2f)",t_scat,-scat_idx);
   logmsg("scint_bw     \t= %.2g MHz",scint_bw);
   logmsg("Random Seed  \t= 0x%"PRIx64,seed);

   // try and open the input files.
   input=fopen(argv[1],"r");
   if(!input){
	  logerr("Could not open input .fil file '%s'",argv[1]);
	  print_help();
	  exit(2);
   }
   prof_file=fopen(prof_fname,"r");
   if(!prof_file){
	  logmsg("Could not open input profile file '%s'",prof_fname);
	  print_help();
	  exit(3);
   }


   output=stdout;

   // read the tempo2 predictor.
   res = T2Predictor_Read(&pred, pred_fname);
   if(res!=0){
	  fprintf(stderr,"error, could not read predictor %d\n",res);
	  exit(4);
   }

   int_fast32_t hdrsize = read_header(input);
   const uint_fast32_t nchan_const = nchans;
   const uint64_t nsamp = nsamples(argv[1],hdrsize,nbits,nifs,nchan_const);

   if(!hdrsize){
	  fprintf(stderr,"error, could not read sigproc header\n");
	  exit(255);
   }
   logmsg("Input Nsamples  \t= %"PRIu64,nsamp);
   logmsg("Input Tobs      \t= %f",nsamp*tsamp);
   logmsg("Input Nchans    \t= %"PRIuFAST32,nchan_const);

   fseek(input,0,SEEK_SET);
   // read/write the sigproc header.
   char* buf = malloc(sizeof(char)*hdrsize); // alloc buf
   fread(buf,1,hdrsize,input);
   fwrite(buf,1,hdrsize,output);
   free(buf);

   n=0;
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
   logmsg("Initialising convolution plans.");
   struct convolve_plan *ism_conv_plan = setup_convolve(nprof,unsmeared_prof,ism_conv,smeared_prof);
   struct convolve_plan *subpulse_conv_plan = setup_convolve(nprof,subpulse_profile,subpulse_map,unsmeared_prof);

   // reread the profile
   logmsg("read profile file");
   fseek(prof_file,0,SEEK_SET);
   n=0;
   while(!feof(prof_file)){
	  fscanf(prof_file,"%f\n",profile+n);
	  subpulse_profile[n]=1.0;
	  n++;
   }
   fclose(prof_file);
   if(strlen(subprof_fname)){

	  logmsg("read subprofile file");
	  prof_file=fopen(subprof_fname,"r");
	  n=0;
	  while(!feof(prof_file)){
		 fscanf(prof_file,"%f\n",subpulse_profile+n);
		 n++;
	  }
	  fclose(prof_file);
   }

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
   const double noiseAmp = A;


   sum=0;
   // normalise the profile
   for (i=0; i < nprof; i++)
	  sum+=profile[i];
   double scale=1.0;

   sum/=(double)nprof;
   scale = noiseAmp*in_snr / sqrt(nchan_const*nsamp)/sum;
   logmsg("mean=%lg scale=%lg",sum,scale);
   // normalise to pseudo-S/N in 1s.
   for (i=0; i < nprof; i++)
	  profile[i]=profile[i]*scale;


   sum=0;
   // normalise the subpulse profile
   for (i=0; i < nprof; i++)
	  sum+=subpulse_profile[i];
   sum/=(double)nprof;
   // we also normalise by the number of subpulses so that the final subpulse profile will have an area of 1.
   scale = 1.0/sum/(float)nsubpulse;

   for (i=0; i < nprof; i++)
	  subpulse_profile[i]=subpulse_profile[i]*scale;

   i=0;
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
   float **ism_models= (float**)malloc(sizeof(float*)*nchan_const);
   for (i=0; i < nchan_const; i++){
	  ism_models[i] = (float*)malloc(sizeof(float)*nprof);
   }
   uint_fast32_t nism=0;
   const uint_fast32_t NCOPY = nprof*sizeof(float);
   const double psr_freq = (double)T2Predictor_GetFrequency(&pred,mjd,freq[0]);
   double poff[nchan_const];

   const double pulse_energy_norm = exp(pow(pulse_energy_sigma,2)/2.0);
   double p0 = T2Predictor_GetPhase(&pred,mjd,fch1);
   int_fast32_t prevbin[nchan_const];
   int_fast32_t dbin,ii;


   logmsg("Pulsar period ~%.2lfms",1000.0/psr_freq);
   mjk_rand_t **rnd = malloc(sizeof(mjk_rand_t*)*nchan_const);
   rnd[0] = mjk_rand_init(seed);
   for(ch = 0; ch < nchan_const; ch++){
	  freq[ch] = fch1 + foff*ch;
	  poff[ch] = T2Predictor_GetPhase(&pred,mjd,freq[ch])-p0;
	  prevbin[ch]=-1;
	  sidx[ch] = pow(freq[ch]/ref_freq,spec_index);
   }

   if(t_scat > 0){
	  float complex *scint_model = fftwf_alloc_complex(nchan_const*2);
	  float *scint_out = fftwf_alloc_real(nchan_const*2);
	  fftwf_plan fft_plan = fftwf_plan_dft_c2r_1d(nchan_const*2,scint_model,scint_out,FFTW_MEASURE);
	  double t = 0;
	  for(ch = 0; ch < nchan_const*2; ch++){
		 A = exp(-t/t_scat);
		 scint_model[ch] = A*mjk_rand_gauss(rnd[0]) + I*A*mjk_rand_gauss(rnd[0]);
		 t += 1.0/fabs(foff*1e6*nchan_const*2);
	  }

	  fftwf_execute(fft_plan);

	  sum=0;
	  for(ch = 0; ch < nchan_const; ch++){
		 scint_out[ch]*=scint_out[ch];
		 sum+=scint_out[ch];
	  }
	  sum/=(double)nchan_const;
	  for(ch = 0; ch < nchan_const; ch++){
		 sidx[ch]*=scint_out[ch]/sum;
	  }

	  fftwf_destroy_plan(fft_plan);
	  fftwf_free(scint_out);
	  fftwf_free(scint_model);

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
	  const float tscatN = (t_scat*pow(freq[nchan_const-1]/ref_freq,-scat_idx));
	  logmsg("Tscatter@%.1fMHz = %gs",ref_freq,t_scat);
	  logmsg("Tscatter@%.1fMHz = %gs",freq[0],tscat0);
	  logmsg("Tscatter@%.1fMHz = %gs",freq[nchan_const-1],tscatN);
	  for(ch = 0; ch < nchan_const; ch++){
		 phase = T2Predictor_GetPhase(&pred,mjd,freq[ch]) - T2Predictor_GetPhase(&pred,mjd,freq[ch]+foff);
		 pbin = fabs(phase)*nprof;
		 if (pbin<1)pbin=1;
		 if (pbin%2==0)pbin-=1;

		 sbin = (int_fast32_t)floor(psr_freq*(t_scat*pow(freq[ch]/ref_freq,-scat_idx))*nprof);

		 if (pbin!=dm_bin || sbin != scatter_bin){
			// make new ISM model
			//logmsg("New ISM Model...");
			//logmsg("DM smearing. df=%lg dt=%lLg phase bins=%d",foff,phase,pbin);
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


			//logmsg("Scattering. Exponential with phase bins=%d idx=%.1f",sbin,scat_idx);
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
   logmsg("Generated %d ISM models",nism);

   float  grads[nism][nprof];
   float output_prof[nism][nprof];
   stop_clock(CLK_setup);
   start_clock(CLK_process);

   logmsg("Starting simulation");
   uint64_t count=0;
   int64_t ip0;
   int64_t prev_ip0=INT64_MAX;
   float rands[nchan_const];
   while(!feof(input)){
	  if(count%1024==0){
		 fprintf(stderr,"%ld samples,  %.1fs\r",count,count*tsamp);
	  }

	  for(ch = 0; ch < nchan_const; ch++){
		 rands[ch]=mjk_rand_double(rnd[0]);
	  }
	  read_block(input,nbits,block,nchan_const);
	  p0 = T2Predictor_GetPhase(&pred,mjd,freq[0]);
	  ip0 = (int64_t)floor(p0);
	  //pragma omp paralell private(i,n,ch,phase,frac,pbin,A,dbin) shared(ip0,prev_ip0)
	  {
		 if (ip0!=prev_ip0){
			//logmsg("new pulse");
			// need a new temp profile.

			//pragma omp  for 
			for(i=0; i < nprof;i++){
			   subpulse_map[i]=0;
			}
			//pragma omp single
			for(n=0; n < nsubpulse;n++){
			   i=floor(mjk_rand_double(rnd[0])*nprof);
			   subpulse_map[i]+=exp(mjk_rand_gauss(rnd[0])*pulse_energy_sigma)/pulse_energy_norm;
			   //subpulse_map[i]+=(mjk_rand_double(rnd[0]))*2;
			}

			//pragma omp single
			{
			   convolve(subpulse_conv_plan);
			}
			//pragma omp  for 
			for(i=0; i < nprof;i++){
			   unsmeared_prof[i]*=profile[i];
			}

			//pragma omp single
			{
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
		 }
		 start_clock(CLK_inner);
		 //pragma omp  for  schedule(dynamic,128)
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
			   block[ch] += floor(sidx[ch]*A+rands[ch]);
			}
			prevbin[ch]=pbin;

		 }

		 stop_clock(CLK_inner);
	  }
	  write_block(nbits,1,nchan_const,output,block);
	  mjd+=tsamp_mjd;
	  count++;
   }
   fprintf(stderr,"\n");

   stop_clock(CLK_process);
   logmsg("Simulation ended.");
   logmsg("Setup took   %.2f s",read_clock(CLK_setup));
   logmsg("Process took %.2f s",read_clock(CLK_process));
   logmsg("Inner took  %.2f s",read_clock(CLK_inner));
   logmsg("Total was  %0.2g times 'real' time",(read_clock(CLK_setup)+read_clock(CLK_process))/(count*tsamp));
   logmsg("Check last Random %"PRIx32,mjk_rand(rnd[0]));

   for (i=0; i < nchan_const; i++){
	  free(ism_models[i]);
   }
   free(ism_models);

   mjk_rand_free(rnd[0]);
   //for(ch = 0; ch < nchan_const; ch++){
   //}
   free(rnd);
   free_convolve_plan(ism_conv_plan);
   free_convolve_plan(subpulse_conv_plan);
   fftwf_free(profile);
   fclose(input);

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
void free_convolve_plan(struct convolve_plan *plan){
   fftwf_free(plan->c);
   fftwf_free(plan->d);
   fftwf_destroy_plan(plan->fwda);
   fftwf_destroy_plan(plan->fwdb);
   fftwf_destroy_plan(plan->bwd);
   free(plan);
}


