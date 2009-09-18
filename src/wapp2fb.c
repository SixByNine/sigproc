/* 
   wapp2fb - converts WAPP search-mode data into "filterbank" data 

   Created from earlier wapp2fbfm code (which itself was based on routines 
   written by Andy Dowd and Jeff Hagen) --- dunc@naic.edu --- Aug 23, 2000.

   Aug 18, 2004 - added ALFA demultiplexing mode
   Mar 24, 2005 - added Jeff's code to get positions if not in header
*/

#include "filterbank.h"
#include "wapp_header.h"
#include <string.h>
#include <fcntl.h>

#include "alfa_position.c"

/* 
   From version 2.6 onwards, FFTW can be substituted if version 3 of
   the libraries has been installed on the system. The compiler directive
   FFTW is used to switch on these libraries --- can give 30% speed-up
   for this particular application (real-to-real FFT)
*/
#ifdef FFTW
#include <fftw3.h>
#endif

/* time between correlator dumps in us */
#define WAPP_DEAD_TIME 0.34

int wapp_lagtrunc, wapp_file, wapp_lagformat, wapp_sum, wapp_level, wapp_flip;
int wapp_header_size, wapp_incfile_length, wapp_isalfa, wapp_number;
double wapp_obstime, alfa_ang, alfa_raj[7], alfa_dej[7];

void wapp2fb(FILE *input, FILE *output) /* includefile */
{
  FILE *bptr, *fpou, *alfa[2];
  int pixel[2];
  double pra, pdec;
  double bw, bandwidth, scale, power, hweight, tsamp_us, crate, lmst;
  double *lag, *sum, *acf, *window, jan1, days, epoch, ras,des,rahr,dede; 
  float zerolag,*block,smin,smax;
  int doit,i,j,k,two_nlags,nlags,stat,rec_size,idump,swap_bytes,ifnum,opened;
  int filesize,headersize,beam,utsecs,iymdf[4],rah,ram,ded,dem;
  unsigned char *cblock, zuc;
  unsigned short *sblock, zus; 
  unsigned long zul;
  char message[80], outfile[80];
  void *dump;
  static float realtime=0.0;

#ifdef FFTW
  fftw_plan fftplan;
#endif
  /* establish whether we need to swap bytes (WAPP is little endian) */
  swap_bytes=big_endian();

  /* initialise correlator parameters used below */
  nlags=nchans;  two_nlags=2*nlags;  bandwidth=foff*nlags;
  tsamp_us=tsamp*1.0e6;
  if (bandwidth<0.0) bandwidth *= -1.0;

#ifdef FFTW
  acf = fftw_malloc(sizeof(double) * two_nlags);
  lag = fftw_malloc(sizeof(double) * two_nlags);
  /* set up fftw table and acf array when computing power spectra */
  fftplan=fftw_plan_r2r_1d(two_nlags,acf,lag,FFTW_R2HC,FFTW_PATIENT);
#endif
#ifndef FFTW
  /* set up acf array when computing power spectra */
  acf = (double *) malloc(two_nlags * sizeof(double));
  lag = (double *) malloc(two_nlags * sizeof(double));
#endif

  if (compute_spectra) {
    /* ranges for scaling spectra */
    smin=0.0;smax=3.0;
  } else {
    /* ranges for scaling correlation functions */
    smin=-0.5;smax=1.0;
  }

  /* set up the weights for windowing of ACF to monimize FFT leakage */
  if (hanning) {
    /* Hanning window */
    hweight=0.50;
  } else if (hamming) {
    /* Hamming window */
    hweight=0.54;
  } else {
    /* no window (default) */
    hweight=1.00;
  }

  /* define the smoothing window to be applied base on the above weight */
  window = (double *) malloc(nlags * sizeof(double));
  for (j=0; j<nlags; j++) window[j]=(hweight+(1.0-hweight)*cos(PI*j/nlags));

  /* work out number of IFs to loop over */
  if (sumifs && (nifs>1)) {
    smin*=2.0;
    smax*=2.0;
    ifnum=2;
  } else {
    sumifs=0;
    ifnum=nifs;
  }

  /* calculate required record size for reading - i.e. number of bytes/dump */
  rec_size = nifs*nlags*(nbits/8);
  dump = malloc(rec_size); /* pointer to the correlator dump */

  /* correlator data rate */
  crate = 1.0/(tsamp_us-WAPP_DEAD_TIME); 
  /* scale factor to normalize correlation functions */
  if (bandwidth < 50.0) 
    bw=50.0; /* correct scaling for narrow-band use */
  else
    bw=bandwidth;

  scale = crate/bw;
  if (wapp_level==9) scale/=16.0; /* 9-level sampling */
  if (wapp_sum) scale/=2.0;  /* summed IFs (search mode) */
  scale*=pow(2.0,(double)wapp_lagtrunc); /* needed for truncation modes */

  /* now define a number of working arrays to store lags and spectra */
  block = (float *) malloc(nlags * sizeof(float));
  cblock = (unsigned char *) malloc(nlags * sizeof(unsigned char));
  sblock = (unsigned short *) malloc(nlags * sizeof(unsigned short));

  /* if the file is ALFA data --- do the demultiplexing to two files */
  if (wapp_isalfa) {

    angle_split(src_raj,&rah,&ram,&ras);
    rahr=(double)rah+(double)ram/60.0+(double)ras/3600.0;
    angle_split(src_dej,&ded,&dem,&des);
    if (ded>0)
      dede=(double)ded+(double)dem/60.0+(double)des/3600.0;
    else
      dede=(double)ded-(double)dem/60.0-(double)des/3600.0;
    /* calculate local sidereal time in hours */
    lmst=slaGmst(tstart)*12.0/4.0/atan(1.0)-4.4502051459439667;
    if (lmst<0.0) lmst+=24.0;
    slaDjcal(5,tstart,iymdf,&stat);
    slaCaldj(iymdf[0],1,1,&jan1,&stat);
    days=tstart-jan1+1.0;
    epoch=(double)iymdf[0]+days/365.25;
    utsecs=86400*(tstart-floor(tstart));

    pixel[0]=(wapp_number-1)*2;
    pixel[1]=pixel[0]+1;
    for (i=0; i<2; i++) {
      if (alfa_raj[pixel[i]] == 0.0) {
      alfa_position(rahr,dede,lmst,epoch,alfa_ang,0.0,0.0,pixel[i],&pra,&pdec);
      src_raj=h2hms(pra);
      src_dej=deg2dms(pdec);
      } else {
      src_raj=h2hms(alfa_raj[pixel[i]]);
      src_dej=deg2dms(alfa_dej[pixel[i]]);
      }
      sprintf(outfile,"%s_%.0f_%05d_%04d_%s_%d.fil",
	    project,floor(tstart),utsecs,
	    scan_number,source_name,pixel[i]);
      if (!file_exists(outfile)) {
        if (i==0) puts("opening output files for demultiplexed ALFA data...");
        alfa[i]=open_file(outfile,"wb");
        puts(outfile);
        filterbank_header(alfa[i]);
      } else {
        if (i==0) puts("appending output files for demultiplexed ALFA data...");
        puts(outfile);
        alfa[i]=open_file(outfile,"ab");
      }
    }
    beam=0;
  }

  if (headerfile) {
    /* write output ASCII header file */
    fpou=open_file("head","w");  
    fprintf(fpou,"Original WAPP file: %s\n",inpfile);
    fprintf(fpou,"Sample time (us): %f\n",tsamp_us);
    fprintf(fpou,"Observation time (s): %f\n",wapp_obstime);
    fprintf(fpou,"Time stamp (MJD): %18.12f\n",tstart);
    fprintf(fpou,"Number of samples/record: %d\n",512);
    fprintf(fpou,"Center freq (MHz): %f\n",fch1+(float)nlags*foff/2.0);
    fprintf(fpou,"Channel band (kHz): %f\n",bandwidth*1000.0/nlags);
    fprintf(fpou,"Number of channels/record: %d\n",nlags);
    fprintf(fpou,"Nifs: %d\n",ifnum);
    fprintf(fpou,"RA (J2000): %f\n",src_raj);
    fprintf(fpou,"DEC (J2000):  %f\n",src_dej);
    fprintf(fpou,"Gal l: %.4f\n",srcl);
    fprintf(fpou,"Gal b: %.4f\n",srcb); 
    fprintf(fpou,"Name: %s\n",source_name);
    fprintf(fpou,"Lagformat: %d\n",wapp_lagformat);
    fprintf(fpou,"Sum: %d\n",wapp_sum);
    fprintf(fpou,"Level: %d\n",wapp_level);
    fprintf(fpou,"AZ at start: %f\n",az_start);
    fprintf(fpou,"ZA at start: %f\n",za_start);
    fprintf(fpou,"AST at start: %f\n",ast0);
    fprintf(fpou,"LST at start: %f\n",lst0);
    fprintf(fpou,"Project ID: %s\n",project);
    fprintf(fpou,"Observers: %s\n",culprits);
    filesize=sizeof_file(inpfile);
    fprintf(fpou,"File size (bytes): %d\n",filesize);
    headersize=wapp_header_size+wapp_incfile_length;
    fprintf(fpou,"Data size (bytes): %d\n",filesize-headersize);
    fprintf(fpou,"Number of samples: %d\n",nsamples(inpfile,headersize,nbits,nifs,nchans));
    fclose(fpou);
  }

  /* initialise various counters and flags */
  opened=idump=i=j=0; 

  /* main loop reading data from infile until no more left to read */
  while( (stat=read(wapp_file,dump,rec_size)) == rec_size) {

    /* calculate elapsed time and determine whether we process this record */
    realtime += (float) tsamp;
    if ( (doit=process(realtime,start_time,final_time)) == -1) break;

    if (doit) {

      /* set ALFA beam output if necessary */
      if (wapp_isalfa) {
	output=alfa[beam];       /* set output file for this loop */
	beam=!(beam);            /* flip file for next iteration */
      }

      /* clear zerolag and blocksum arrays */
      zerolag=0.0;
      for (j=0; j<nlags; j++) block[j]=0.0;

      /* loop over the IFs */
      for (i=0; i<ifnum; i++) {
	if (ifstream[i]=='Y') {
	  if (zerolagdump) {
	    /* only interested in the zero lag term for each IF */
	    switch (nbits) {
	    case 8:
	      zuc = *(((unsigned char *)dump)+i*nlags);
	      zerolag+=zuc;
	      break;
	    case 16:
	      zus = *(((unsigned short *)dump)+i*nlags);
	      if (swap_bytes) swap_short(&zus); 
	      zerolag+=zus;
	      break;
	    case 32:
	      zul = *(((unsigned long *)dump)+i*nlags);
	      if (swap_bytes) swap_long(&zul); 
	      zerolag+=zul;
	      break;
	    }
	    /* write out the data checking IF number for summed mode */
	    if ( (sumifs && (i==1)) || (!sumifs) ) {
	      if (obits==32) {
		if (swapout) swap_float(&zerolag);
		fwrite(&zerolag,sizeof(float),1,output);
	      } else {
		sprintf(message,"cannot write %d bits in zerolag mode",obits);
		error_message(message);
	      }
	    }
	  } else {
	    /* fill lag array with scaled CFs */
	    for (j=0; j<nlags; j++) {
	      switch (nbits) {
	      case 8:
		zuc = *(((unsigned char *)dump)+j+i*nlags);
		lag[j] = scale * (double) zuc - 1.0;
		break;
	      case 16:
		zus = *(((unsigned short *)dump)+j+i*nlags);
		if (swap_bytes) swap_short(&zus);
		lag[j] = scale * (double) zus - 1.0;
		break;
	      case 32:
		zul = *(((unsigned long  *)dump)+j+i*nlags);
		if (swap_bytes) swap_long(&zul);
		lag[j] = scale * (double) zul - 1.0;
		break;
	      }
	    }
	    /* calculate power and correct for finite level quantization */
	    power = inv_cerf(lag[0]);
	    power = 0.1872721836/power/power;
	    if (i<2) {
	      if (do_vanvleck) { 
		if (wapp_level==3) {
		  /* apply standard 3-level van vleck correction */
		  vanvleck3lev(lag,nlags);
		} else if (wapp_level==9) {
		  /* apply 9-level van vleck correction */
		  vanvleck9lev(lag,nlags);
		}
	      }
	    }

	    if (compute_spectra) {
	      /* form windowed even ACF in array */
	      for(j=1; j<nlags; j++) {
		acf[j]=window[j]*lag[j]*power;
		acf[two_nlags-j]=acf[j];
	      }   
	      acf[nlags]=0.0;
	      acf[0]=lag[0]*power; 
	      /* FFT the ACF (which is real and even) -> real and even FFT */
#ifdef FFTW
	      fftw_execute(fftplan);
#endif
#ifndef FFTW
	      rfft(two_nlags,acf,lag);
#endif
	      /* if the band needs to be flipped --- do it here */
	      if (wapp_flip) {
		/* use acf as temporary array */
		for (j=0;j<nlags;j++) acf[j]=lag[j]; 
		k=nlags-1;
		for (j=0;j<nlags;j++) {
		  lag[k]=acf[j];
		  k--;
		}
	      }
  	      /* add lags to block array */
	      for (j=0; j<nlags; j++) block[j]+=lag[j];
	    } else {			
	      /* just copy correlation functions into block */
	      for (j=0; j<nlags; j++) block[j]=lag[j];	
	    }
	    /* write out data block checking IF number for summed mode */
	    if ( (sumifs && (i==1)) || (!sumifs) ) {
	      if (obits==32) {
		if (swapout) for (j=0; j<nlags; j++) swap_float(&block[j]);
		fwrite(block,sizeof(float),nlags,output);
	      } else if (obits==16) {
		float2short(block,nlags,smin,smax,sblock);
		if (swapout) for (j=0; j<nlags; j++) swap_short(&sblock[j]);
		fwrite(sblock,sizeof(unsigned short),nlags,output);
	      } else if (obits==8) {
		float2char(block,nlags,smin,smax,cblock);
		fwrite(cblock,sizeof(unsigned char),nlags,output);
	      } else if (obits==4) {
		float2four(block,nlags,smin,smax,cblock);
		fwrite(cblock,sizeof(unsigned char),nlags/2,output);
	      } else {
		sprintf(message,"cannot write %d bits in wapp2fb",obits);
		error_message(message);
	      }
	    }
	  } /* end of zerolagdump if */

	  if (!sumifs) { /* reset block and zerolag if not summing */
	    zerolag=0.0;
	    for (j=0; j<nlags; j++) block[j]=0.0;
	  }
	}          /* end of IFstream if */
      }            /* end of loop over IFs */
    }              /* end of processing if */

    /* increment dump counter and update logfile every 512 dumps */
    idump++;
    if (idump%512 == 0) {
      if (!opened) {
	/* open up logfile */
	open_log("filterbank.monitor");
	opened=1;
      }
      sprintf(message,"time:%.1fs",realtime);
      update_log(message);
    }
  }                /* end of main read loop*/
  if (wapp_isalfa) {
      fclose(alfa[0]);
      fclose(alfa[1]);
  }

    
  /* job done - free up remaining arrays */
  free(dump);free(block);free(sblock);free(cblock);free(window);
#ifdef FFTW
  fftw_destroy_plan(fftplan);
  fftw_free(acf); fftw_free(lag);
#endif
#ifndef FFTW
  free(acf); free(lag);
#endif
}
