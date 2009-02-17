/* getpulse.c - reads SIGPROC dedispersed time series and outputs/plots a pulse */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
#include "cpgplot.h"
FILE *input;
main(int argc, char *argv[])
{
  FILE *output;
  int numerate,i,j,k,l,nsperdmp,nsamps,indexing,plot,indexnow;
  int ifchan[16],frchan[4096],ifnum,chnum,ns,sample,ntotal;
  float time_of_pulse, width_of_pulse, time, time_start, time_end;
  float series_amp[100000],series_time[100000],ampmin,ampmax;
  char message[80],outfile[80],filename[80],timestring[80],pgdev[80];
  unsigned char c;
  unsigned short s;
  float f[8];

    if (argc<2) {
    puts("");
    puts("getpulse - make and/or plot a time series from dedispersed file\n"); 
    puts("usage: getpulse {filename} -{options}\n");
    puts("filename is the dedispersed data file\n");
    puts("options:\n");
    puts("-t time - time (in seconds) on which to center time series (REQUIRED)\n");
    puts("-w width - width (in seconds) of time series to plot (def=1)\n");
    puts("-c fchan - to only output frequency channel fchan (def=all)\n");
    puts("-i ifchan - to only output IF channel ifchan (def=all)\n");
    puts("-p pgdev - pgplot device (def=/xs)\n");
    puts("-numerate - to precede each dump with sample number (def=time)\n");
    puts("-noindex - do not precede each dump with time/number\n");
    puts("-plot - PGPLOT the results\n");
    puts("");
    exit(0);
   }

  /* zero IF and frequency channel arrays */
  for (i=0;i<16;i++)   ifchan[i]=0;
  for (i=0;i<4096;i++) frchan[i]=0;

  /* default case is to read from standard input */
  input=stdin;
  numerate=0;
  ampmin = 1.e6;
  ampmax = -1.e6;
  plot=0;
  indexnow=0;
  indexing=1;
  strcpy(pgdev,"/xs");

  /* parse command line if arguments were given */
  if (argc>1) {
    i=1;
    while (i<argc) {
        if (strings_equal(argv[i],"-i")) {
	i++;
	ifchan[atoi(argv[i])-1]=1;
      } else if (strings_equal(argv[i],"-c")) {
	i++;
	frchan[atoi(argv[i])-1]=1;
      } else if (strings_equal(argv[i],"-t")) {
        i++;
	time_of_pulse = atof(argv[i]);
        strcpy(timestring,argv[i]);
	fprintf(stderr,"centering on time %f s\n",time_of_pulse); 
      } else if (strings_equal(argv[i],"-w")) {
        i++;
        width_of_pulse = atof(argv[i]);
        fprintf(stderr,"with width %f s\n",width_of_pulse);
      } else if (strings_equal(argv[i],"-p")) {
        i++;
        strcpy(pgdev,argv[i]); 
      } else if (strings_equal(argv[i],"-numerate")) {
	numerate=1;
      } else if (strings_equal(argv[i],"-noindex")) {
	indexing=0;
      } else if (strings_equal(argv[i],"-plot")) {
	plot=1;
      } else if (file_exists(argv[i])) {
	input=open_file(argv[i],"rb");
        strcpy(filename,argv[i]);
      } else {
	sprintf(message,"unknown argument (%s) passed to getpulse",argv[i]);
	error_message(message);
      }
      i++;
    }
  }

  strcat(filename,".pulse.");
  strcat(filename,timestring);
  strcpy(outfile,filename);
  output=fopen(outfile,"w");
  time_start = time_of_pulse - width_of_pulse/2.;
  time_end = time_of_pulse + width_of_pulse/2.;
  sample = 0;
  /* try to read the header */
  if (!read_header(input)) error_message("error reading header\n");

  /* check what IF and frequency channels (if any the user has selected) */
  j=0;
  for (i=0; i<nifs; i++) if (ifchan[i]) j++;
  if (j==0) for (i=0; i<nifs; i++) ifchan[i]=1;
  j=0;
  for (i=0; i<nchans; i++) if (frchan[i]) j++;
  if (j==0) for (i=0; i<nchans; i++) frchan[i]=1;

  /* number of samples to read per dump */
  nsperdmp=nifs*nchans;
  /* initialize loop counters and flags */
  ifnum=chnum=nsamps=l=0;

  while (!feof(input)) {

    /* unpack the sample(s) if necessary */
    switch (nbits) {
    case 1:
      fread(&c,1,1,input);
      for (i=0;i<8;i++) {
	f[i]=c&1;
	c>>=1;
      }
      ns=8;
      break;
    case 4:
      fread(&c,1,1,input);
      char2ints(c,&j,&k);
      f[0]=(float) j;
      f[1]=(float) k;
      ns=2;
      break;
    case 8:
      fread(&c,nbits/8,1,input);
      f[0]=(float) c;
      ns=1;
      break;
    case 16:
      fread(&s,nbits/8,1,input);
      f[0]=(float) s;
      ns=1;
      break;
    case 32:
      fread(&f[0],nbits/8,1,input);
      ns=1;
      break;
    default:
      sprintf(message,"cannot read %d bits per sample...\n",nbits);
      error_message(message);
      break;
    }

    if (plot) {
       ntotal = (time_end-time_start)/tsamp;
       if (ntotal > 100000) {
	  fprintf(stderr,"Too many samples to plot!\n");
	  plot=0;
	}
    }
    for (i=0; i<ns; i++) {
      /* time stamp or index the data */
      time = (double) tsamp * (double) l;
      if (time > time_end) {
         if (plot) {
	       indexnow=indexnow-1;
	       cpgbeg(0,pgdev,1,1);
               cpgsvp(0.1,0.9,0.1,0.9);
               cpgswin(series_time[0],series_time[indexnow],ampmin-0.1*ampmax,ampmax+0.1*ampmax);
               cpgbox("bcnst",0.0,0,"bcnst",0,0.0);
               cpgline(indexnow,series_time,series_amp);
               cpgscf(2);
	       cpgmtxt("B",2.5,0.5,0.5,"Time (s)");
	       cpgmtxt("L",1.8,0.5,0.5,"Amplitude");
               cpgend();}
	    exit(0);
      }
      /* print sample if it is one of the ones selected */
      if (ifchan[ifnum] && frchan[chnum] && time >= time_start && time <= time_end) {
	if (indexing) {
           if (numerate) fprintf(output,"%d %f\n",l,f[i]);
           else fprintf(output,"%f %f\n",time,f[i]);
        }	
	else {
	   fprintf(output,"%f\n",f[i]);
	}
	if (plot) {
	   series_amp[indexnow]=f[i];
	   series_time[indexnow]=time;
           if (f[i] < ampmin) ampmin = f[i];
           if (f[i] > ampmax) ampmax = f[i];
           indexnow++;
	}
      } 
      nsamps++;
      chnum++;
      if (chnum==nchans) {
	chnum=0;
	ifnum++;
	if (ifnum==nifs) ifnum=0;
      }
      if ((nsamps%nsperdmp)==0) {
	nsamps=0;
	l++;
      }
      }
	}
    fclose(output);

}

