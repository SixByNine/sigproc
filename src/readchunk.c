/* reader.c - reads SIGPROC data files and formats human output */
/* 14/7/2004 - added byte option for converting time series to single bytes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
FILE *input;
main(int argc, char *argv[])
{
  int numerate,i,j,k,l,stream,nsperdmp,nsamps,indexing,indexnow;
  int ifchan[16],frchan[4096],ifnum,chnum,ns,charout,sample;
  float time_of_pulse, width_of_pulse, time, time_start, time_end;
  char message[80],byte;
  unsigned char c;
  unsigned short s;
  float f[8];

  /* zero IF and frequency channel arrays */
  for (i=0;i<16;i++)   ifchan[i]=0;
  for (i=0;i<4096;i++) frchan[i]=0;

  /* default case is to read from standard input */
  input=stdin;
  charout=numerate=stream=0;
  indexing=1;

  /* parse command line if arguments were given */
  if (argc>1) {
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"help")) {
	reader_help();
	exit(0);
      } else if (strings_equal(argv[1],"version")) {
	printf("PROGRAM: %s SIGPROC version: %.1f\n",argv[0],SIGPROC_VERSION);
	exit(0);
      } else if (strings_equal(argv[i],"-i")) {
	i++;
	ifchan[atoi(argv[i])-1]=1;
      } else if (strings_equal(argv[i],"-c")) {
	i++;
	frchan[atoi(argv[i])-1]=1;
      } else if (strings_equal(argv[i],"-t")) {
        i++;
	time_of_pulse = atof(argv[i]);
	fprintf(stderr,"centering on time %f s\n",time_of_pulse); 
      } else if (strings_equal(argv[i],"-w")) {
        i++;
        width_of_pulse = atof(argv[i]);
        fprintf(stderr,"with width %f s\n",width_of_pulse); 
      } else if (strings_equal(argv[i],"-numerate")) {
	numerate=1;
      } else if (strings_equal(argv[i],"-noindex")) {
	indexing=0;
      } else if (strings_equal(argv[i],"-stream")) {
	stream=1;
      } else if (strings_equal(argv[i],"-byte")) {
	charout=1;
      } else if (file_exists(argv[i])) {
	input=open_file(argv[i],"rb");
      } else {
	reader_help();
	sprintf(message,"unknown argument (%s) passed to reader",argv[i]);
	error_message(message);
      }
      i++;
    }
  }

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
  indexnow=1;
  if (stream && indexing) numerate=1;

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

    for (i=0; i<ns; i++) {
      if (charout) {
	byte=(char) f[i];
	putchar(byte);
      } else {
      /* time stamp or index the data if needed */
	time = (double) tsamp * (double) l;
	if (time > time_end) exit(0);
      if (indexnow && stream) {
	puts("#START");
      } else if (indexnow && indexing && time >= time_start && time <= time_end) {
	/* if (numerate) 
	  printf("%d ",l);
	else
	  printf("# %f\n",time); */
      }
      indexnow=0;
      /* print sample if it is one of the ones selected */
      if (ifchan[ifnum] && frchan[chnum] && time >= time_start && time <= time_end) {
	if (stream && numerate) printf("%d ",nsamps);
	printf("%d %d %f \n",l,chnum,f[i]);
      } 
      nsamps++;
      chnum++;
      if (chnum==nchans) {
	chnum=0;
	ifnum++;
	if (ifnum==nifs) ifnum=0;
      }
      /* put newline and terminator if in streaming mode */
      if (stream) {
	if ((nsamps%nchans)==0) puts("#STOP");
      } 
      /* put newline if this is the last sample of the dump */
      if ((nsamps%nsperdmp)==0) {
	nsamps=0;
	indexnow=1;
	l++;
      }
      }
    }
  }
}

