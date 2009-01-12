/*
  SPLICE - splice several filterbank files starting on same time stamp together
  for use to analyse data from cloned machines sampling several parts of the 
  band. e.g. multi-WAPPs. 
  Initial version (Nov 2002) required contiguous bands. 
  New version     (Jun 2003) requires only high-lo frequency ordering
  Modified        (Jul 2005) to write out to file if -o file option is given
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
FILE *output;
main (int argc, char **argv)
{
  int i=1, j, k, nfiles=0, *numbt, schans=0, nbytes, *nchan;
  FILE *input[32];
  char *block;
  double *stamp, *frch1, *froff, *frmhz;
  output=stdout;
  /* print help if necessary */
  if (argc<=1 || help_required(argv[1])) {
    /*splice_help();*/
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }

  /* open up files */
  while (i<argc) {
    if (file_exists(argv[i])) {
      input[nfiles]=open_file(argv[i],"rb");
      nfiles++;
    } else if (strings_equal(argv[i],"-o")) {
      output=open_file(argv[++i],"wb");
    }
    i++;
  }

  /* read in headers and check time stamps */
  stamp = (double *) malloc(nfiles*sizeof(double));
  frch1 = (double *) malloc(nfiles*sizeof(double));
  froff = (double *) malloc(nfiles*sizeof(double));
  numbt = (int *) malloc(nfiles*sizeof(int));
  nchan = (int *) malloc(nfiles*sizeof(int));
  for (i=0; i<nfiles; i++) {
    if (read_header(input[i])) {
      stamp[i]=tstart;
      frch1[i]=fch1;
      froff[i]=foff;
      numbt[i]=nbits;
      nchan[i]=nchans;
      schans+=nchans;
    } else {
      error_message("problem reading header parameters");
    }
    if (data_type != 1) 
      error_message("input data are not in filterbank format!");
    if (stamp[i] != stamp[0]) 
      error_message("start times in input files are not identical!");
    if (numbt[i] != numbt[0])
      error_message("number of bits per sample in input files not identical!");
    if (i>0) {
      if (frch1[i] > frch1[i-1]) 
	error_message("input files not ordered in descending frequency!");
    }

  }


  send_string("HEADER_START");
  send_int("machine_id",machine_id);
  send_int("telescope_id",telescope_id);
  send_int("data_type",1);

  send_string("FREQUENCY_START");
  send_int("nchans",schans);
  frmhz = (double *) malloc(sizeof(double)*schans);
  k=0;
  for (i=0; i<nfiles; i++) {
    for (j=0; j<nchans; j++) {
      frmhz[k]=frch1[i]+j*froff[i];
      send_double("fchannel",frmhz[k++]);
    }
  }
  send_string("FREQUENCY_END");
  

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

  nbytes = nchans*nbits/8;
  block = (char *) malloc(nbytes);
  while (1) {
    for (i=0; i<nfiles; i++) {
      if (feof(input[i])) exit(0);
      fread(block,nbytes,1,input[i]);
      fwrite(block,nbytes,1,output);
    }
  }
}
