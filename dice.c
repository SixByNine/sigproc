/*
  DICE - dice up a filterbank file by removing unwanted channels
  Initial version (Feb 2006) adapted from splice.c
  Modified (Feb 23, 2006) to write out fch1 and foff wherever possible
  Modified (Nov 28, 2006) to forcefully put zeros in bad channels and
	write them when force=1 (currently hardwired). The reason for this
	is so that the data can be sub-banded later on. zapped channels are
	written as zeros.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
FILE *output;
main (int argc, char **argv)
{
  int i=1, j, k, channum, *numbt, kchans=0, ns, nbytes, *keep, simple;
  int force=1;
  double *frmhz;
  FILE *input,*keepfile;
  float *block;
  unsigned char *charblock;
  unsigned short *shortblock;

  output=stdout;
  /* print help if necessary */
  if (argc!=3 || help_required(argv[1])) {
    /*dice_help();*/
    puts("");
    puts("usage: dice filterbank_file keep_file");
    puts("");
    puts("keep file is the list of channels to keep in the output");
    puts("");
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }

  if (strings_equal(argv[1],"stdin")) {
	input = stdin;
  } else {
        /* open up file */
  	if (!file_exists(argv[1])) 
    		error_message("input file does not exist...");
  	input=open_file(argv[1],"rb");
  }

  /* read in header info */
  if (!read_header(input)) 
    error_message("problem reading header parameters");

  if (data_type != 1) 
    error_message("input data are not in filterbank format!");

  keep=(int *) malloc(nchans*sizeof(int));

  /* open up file to give channel numbers to keep */
  if (!file_exists(argv[2])) 
    error_message("keep file does not exist...");

  keepfile=open_file(argv[2],"r");
  frmhz = (double *) malloc(sizeof(double)*nchans);
  keep  = (int *) malloc(sizeof(int)*nchans);
  for (i=0;i<nchans;i++) keep[i]=0;
  k=0;
  while (fscanf(keepfile,"%d\n",&channum)==1) {
    keep[channum-1]=1;
    frmhz[k++]=fch1+(channum-1)*foff;
    kchans++;
  }

  /* do a check to see whether fch1 and foff are sufficient to describe
     the diced channels */
  fch1=frmhz[0];
  simple=1;
  for (i=0;i<kchans;i++) if (frmhz[i] != fch1+i*foff) simple=0;

  /* force the output to include the zapped channels */
  if (force) {
	kchans=nchans;
	simple=1;
  }

  /* broadcast header for output file */
  send_string("HEADER_START");
  send_int("machine_id",machine_id);
  send_int("telescope_id",telescope_id);
  send_int("data_type",1);

  if (simple) {
    send_double("fch1",fch1);
    send_double("foff",foff);
    send_int("nchans",kchans);
  } else {
    send_string("FREQUENCY_START");
    send_int("nchans",kchans);
    for (i=0; i<kchans; i++) send_double("fchannel",frmhz[i]);
    send_string("FREQUENCY_END");
  }

  if (!strings_equal(source_name,"")) {
    send_string("source_name");
    send_string(source_name);
  }
  send_coords(src_raj,src_dej,az_start,za_start);
  send_int("nbits",nbits);
  send_double("tstart",tstart);
  send_double("tsamp",tsamp);
  send_int("nifs",nifs);
  send_string("HEADER_END");

  block = (float *) malloc(nchans*sizeof(float));
  charblock = (unsigned char *) malloc(nchans*sizeof(unsigned char));
  shortblock = (unsigned short *) malloc(nchans*sizeof(unsigned short));

  while (read_block(input,nbits,block,nchans)==nchans) {
    k=0;
    if (force) {
      for (i=0; i<nchans; i++) {
	if (keep[i]) 
	  block[k++]=block[i];
	else
	  block[k++]=0;
      }	
    } else {
      for (i=0; i<nchans; i++) if (keep[i]) block[k++]=block[i];
    }
    switch (nbits) {
    case 8:
      for (i=0;i<k;i++) charblock[i]=block[i];
      fwrite(charblock,1,k,output);
      break;
    case 16:
      for (i=0;i<k;i++) shortblock[i]=block[i];
      fwrite(shortblock,2,k,output);
      break;
    default:
    error_message("dice - currently only works with 8 or 16 bit data");
    }
  }

}
