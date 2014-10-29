#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* global variables describing the data */
#include "header.h"

/* list of subroutines and functions */
#include "sigproc.h"
#include "mjklog.h"
#include "mjk_random.h"
#include "mjk_cmd.h"

void fastfake_help(){
      fprintf(stderr,\
		 "\nfast_fake [options] > output.fil\n"\
		 "\n"\
		 "Create a gausian noise fake filterbank file. Writes to stdout.\n"\
		 "Michael Keith (2014) - mkeith@pulsarastronomy.net\n"\
		 "\n"\
		 "OPTIONS:\n"\
		 "   --help, -h          This help text\n"\
		 "   --tobs,-T           Total observation time, s (def=270)\n"\
		 "   --tsamp,-t          Sampling interval, us (def=64)\n"\
		 "   --mjd,-m            MJD of first sample (def=56000.0)\n"\
		 "   --fch1,-F           Frequency of channel 1, MHz (def=1581.804688)\n"\
		 "   --foff,-f           Channel bandwidth, MHz (def=-0.390625)\n"\
		 "   --nbits,-b          Output number of bits, 1,2,4,8,32 (def=2)\n"\
		 "   --nchans,-c         Output number of channels (def=1024)\n"\
		 "   --seed,-S           Random seed (def=time())\n"\
		 "   --name,-s           Source name for header (def=FAKE)\n"\
		 "\n"\
		 "Default parameters make a HTRU-style data file.\n"\
		 "\n");
	  exit(1);
}



char inpfile[80], outfile[80];

FILE *output;


int main (int argc, char *argv[])
{
   char headerless;
   double obstime;

   char string[80];
   long seed=-1;
   
   float fval;
   unsigned char A,B,C,D;
   char help=0;
   uint_fast32_t i;


   /* set up default variables */
   strcpy(inpfile,"stdin");
   strcpy(outfile,"stdout");
   headerless=0;
   machine_id=telescope_id=0;
   machine_id=10;
   telescope_id=4;
   nchans=1024;
   nbits=2;
   tstart=56000.0;
   tsamp=64.0; // microseconds.. will be converted to seconds later
   fch1=1581.804688;
   foff=-0.390625;
   nifs=1;
   nbeams=1;
   ibeam=1;
   obstime=270.0;
   output=stdout;

   help=getB("--help","-h",argc,argv,0);
   obstime=getF("--tobs","-T",argc,argv,obstime);
   tsamp=getF("--tsamp","-t",argc,argv,tsamp);
   tstart=getF("--mjd","-m",argc,argv,tstart);
   fch1=getF("--fch1","-F",argc,argv,fch1);
   foff=getF("--foff","-f",argc,argv,foff);
   nbits=getI("--nbits","-b",argc,argv,nbits);
   nchans=getI("--nchans","-c",argc,argv,nchans);
   seed=getI("--seed","-S",argc,argv,seed);
   char test_mode=getB("--test","-0",argc,argv,0);
   strcpy(outfile,getS("--out","-o",argc,argv,"stdout"));
   strcpy(source_name,getS("--name","-s",argc,argv,"FAKE"));
   getArgs(&argc,argv);
   if (help || argc > 1){
	  for(i=1; i < argc; i++)logerr("Unknown argument '%s'",argv[i]);
	  fastfake_help();
   }

   logmsg("FASTFAKE - M.Keith 2014");
   

   time_t t0 = time(NULL);
   if (seed<0)seed=t0;

   mjk_rand_t *rnd = mjk_rand_init(seed);


   logmsg("tobs          = %lfs",obstime);
   logmsg("tsamp         = %lfus",tsamp);
   logmsg("mjdstart      = %lf",tstart);
   logmsg("freq chan 1   = %lfMHz",fch1);
   logmsg("freq offset   = %lfMHz",foff);
   logmsg("output nbits  = %d",nbits);
   logmsg("random seed   = %ld",seed);
   logmsg("output file   = '%s'",outfile);

   tsamp*=1e-6; // convert tsamp to us

   if(STREQ(outfile,"stdout")){
	  output=stdout;
   } else{
	  output=fopen(outfile,"w");
   }

   if (!headerless) {
	  logmsg("write header");
	  send_string("HEADER_START");
	  send_string("source_name");
	  send_string(source_name);
	  send_int("machine_id",machine_id);
	  send_int("telescope_id",telescope_id);
	  send_int("data_type",1);
	  send_double("fch1",fch1);
	  send_double("foff",foff);
	  send_int("nchans",nchans);
	  send_int("nbits",nbits);
	  send_int("nbeams",nbeams);
	  send_int("ibeam",ibeam);
	  send_double("tstart",tstart);
	  send_double("tsamp",tsamp);
	  send_int("nifs",nifs);
	  if (nbits==8){
		 send_char("signed",OSIGN);
	  }
	  send_string("HEADER_END");
   }

   int bits_per_sample = (nbits * nchans * nifs);
   if (bits_per_sample % 8 != 0){
	  logerr("bits per sample is not a multiple of 8");
	  exit(1);
   }


   int bytes_per_sample = bits_per_sample / 8.0;
   uint64_t nsamples = (uint64_t)(obstime/tsamp+0.5);

   logmsg("Generate %lld samples, %.3lf GiB",nsamples,nsamples*bytes_per_sample/pow(2,30));

   uint64_t onepercent = nsamples/100;
   int percent=0;
   unsigned char* buffer = (unsigned char*) malloc(sizeof(unsigned char)*bytes_per_sample);

   if (nbits==1 || nbits==2 || nbits==4 | nbits==8){
	  // integer samples
	  for(uint64_t samp = 0; samp < nsamples; samp++){
		 if (samp%onepercent ==0){
			double t1=(double)(time(NULL)-t0)+1e-3;
			double bytespersec = samp*bytes_per_sample/t1;
			fprintf(stderr,"Complete: % 3d%%. Sample: % 9lld Real time % 6.1lfs, Sim time % 6.1lfs. Speed % 4.2lfMiB/s\r",percent,samp,t1,(double)samp*tsamp,bytespersec/pow(2,20));
			fflush(stderr);
			percent+=1;
		 }
		 mjk_rand_gauss_atleast(rnd,nchans);
		 const int chanskip = 8/nbits;
//#pragma omp parallel for schedule(dynamic,1)
		 for(uint64_t chan = 0; chan < nchans; chan+=chanskip){
			switch(nbits){
			   case 1:
				  buffer[chan/8] = mjk_rand(rnd)&0xFF; 
				  break;
			   case 2:
				  fval = (mjk_rand_gauss(rnd)+1.5);
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  A = (unsigned char)(mjk_rand_gauss(rnd)+1.5)&0x8;
				  fval = (mjk_rand_gauss(rnd)+1.5);
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  B = ((unsigned char)(mjk_rand_gauss(rnd)+1.5)&0x8)<<2;
				  fval = (mjk_rand_gauss(rnd)+1.5);
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  C = ((unsigned char)(mjk_rand_gauss(rnd)+1.5)&0x8)<<4;
				  fval = (mjk_rand_gauss(rnd)+1.5);
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  D = ((unsigned char)(mjk_rand_gauss(rnd)+1.5)&0x8)<<6;
				  buffer[chan/4]=A|B|C|D;
				  break;
			   case 4:
				  fval = (mjk_rand_gauss(rnd)*24.0+64.0);
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  A = (unsigned char)(mjk_rand_gauss(rnd)*3+7.5)&0xF;
				  fval = (mjk_rand_gauss(rnd)*2.5+6);
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  B = ((unsigned char)(mjk_rand_gauss(rnd)*3+7.5)&0xF)<<4;
				  buffer[chan/2] = A|B;
				  break;
			   case 8:
				  fval = (mjk_rand_gauss(rnd)*24.0+96.0); // more headroom.
				  fval = fmax(fval,0);
				  fval = fmin(fval,255.0);
				  buffer[chan] = (unsigned char)fval;
				  break;
			}
		 }
		 if(!test_mode)fwrite(buffer,sizeof(unsigned char),bytes_per_sample,output);
	  }
   } else if (nbits==32){
	  // float samples
	  float* fbuf = (float*)buffer;
	  for(uint64_t samp = 0; samp < nsamples; samp++){
		 if (samp%onepercent ==0){
			double t1=(double)(time(NULL)-t0)+1e-3;
			double bytespersec = samp*bytes_per_sample/t1;
			fprintf(stderr,"Complete: % 3d%%. Sample: % 9lld Real time % 6.1lfs, Sim time % 6.1lfs. Speed % 4.2lfMiB/s\r",percent,samp,t1,(double)samp*tsamp,bytespersec/pow(2,20));
			fflush(stderr);
			percent+=1;
		 }
		 for(uint64_t chan = 0; chan < nchans; chan++){
			fbuf[chan] = mjk_rand_gauss(rnd);
		 }
		 if(!test_mode)fwrite(buffer,sizeof(unsigned char),bytes_per_sample,output);
	  }
   }


   fprintf(stderr,"\n");
   free(buffer);
   mjk_rand_free(rnd);
   logmsg("Done!");
   fclose(output);

   return 0;
}

