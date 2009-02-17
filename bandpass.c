#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
FILE *input;
main(int argc, char *argv[])
{
  char string[80];
  int increment=0,i,j,k,l,ndumps,ic,x,y,cube,dumpcount;
  double elapsed_time;
  unsigned char c;
  unsigned short s;
  float *f,*g,temp,fchan;
  double secsperdump;

  /* default is to read from standard input and to sum all dumps */
  input=stdin;
  ndumps=cube=dumpcount=0;
  secsperdump=elapsed_time=0.0;

  /* check command-line parameters */ 
  i=1;
  while (i<argc) {
    print_version(argv[0],argv[1]);
    if (strings_equal(argv[i],"-d")) {
      i++;
      ndumps=atoi(argv[i]);
    } else if (strings_equal(argv[i],"-t")) {
      i++;
      secsperdump=atof(argv[i]);
    } else if (help_required(argv[i])) {
      bandpass_help();
      exit(0);
    } else if (strings_equal(argv[i],"-cube")) {
      cube=1;
    } else if (file_exists(argv[i])) {
      input=open_file(argv[i],"rb");
    } else {
      bandpass_help();
      sprintf(string,"unknown argument (%s) passed to bandpass",argv[i]);
      error_message(string);
    }
    i++;
  }

  if (!read_header(input)) 
    error_message("could not read filterbank header...");

  /* set number of dumps to average over if user has supplied seconds */
  if (secsperdump > 0.0) 
    ndumps = (int) rint(secsperdump/tsamp);

  /* initialize buffer for storing bandpass */
  ic=nchans*nifs;
  f = (float *) malloc(sizeof(float)*ic);
  g = (float *) malloc(sizeof(float)*ic);
  for (i=0; i<ic; i++) f[i]=g[i]=0.0;
  k=l=0;
  while (!feof(input)) {
    switch (nbits) {
    case 4:
      fread(&c,1,1,input);
      char2ints(c,&x,&y);
      f[l]=(float) x;
      l++;
      f[l]=(float) y;
      break;
    case 8:
      fread(&c,nbits/8,1,input);
      f[l]=(float) c;
      break;
    case 16:
      fread(&s,nbits/8,1,input);
      f[l]=(float) s;
      break;
    case 32:
      fread(&temp,nbits/8,1,input);
      f[l]=temp;
      break;
    default:
      sprintf(string,"cannot read %d bits per sample...\n",nbits);
      error_message(string);
      break;
    }
    l++;
    if (l==ic) {
      for (i=0;i<ic;i++) g[i]+=f[i];
      k++;
      if (k==ndumps) {
	dumpcount++;
	if (!cube) puts("#START");
	for (i=0; i<nchans; i++) {
	  fchan=fch1+(foff*(double) i);
	  if (cube) printf("%f ",elapsed_time);
	  printf("%f ",fchan);
	  for (j=0; j<nifs; j++) {
	    l=j*nchans+i;
	    printf("%f ",g[l]/(float)ndumps);
	    g[l]=f[l]=0.0;
	  }
	  printf("\n");
	}
	if (!cube) puts("#STOP");
	k=0;
      }
      l=0;
      if (cube) {
	elapsed_time+=tsamp;
	if (dumpcount==nchans) exit(0);
      }
    }
  }

  /* write out final dump if necessary */
  if (!ndumps) {
    ndumps=k;
    for (i=0; i<nchans; i++) {
      fchan=fch1+(foff*(double) i);
      printf("%f ",fchan);
      for (j=0; j<nifs; j++) {
	l=j*nchans+i;
	printf("%f ",g[l]/(float)ndumps);
      }
      printf("\n");
    }
  }
}
