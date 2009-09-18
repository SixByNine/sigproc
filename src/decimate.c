/*
  DECIMATE  - decimate filterbank data by adding channels and/or time samples
*/

#include "decimate.h"

main (int argc, char *argv[])
{
  int i, nc, headersize, headerless=0;
  char string[80];

  /* set up default global variables */
  obits=headerless=naddc=naddt=nsamp=0;
  input=stdin;
  strcpy(inpfile,"stdin");
  output=stdout;
  strcpy(outfile,"stdout");

  if (argc > 1) {
    /* check command-line parameters */ 
    print_version(argv[0],argv[1]);
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"-c")) {
	i++;
	naddc=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-t")) {
	i++;
	naddt=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-o")) {
	/* get and open file for output */
	output=fopen(argv[++i],"wb");
      } else if (strings_equal(argv[i],"-T")) {
	i++;
	nsamp=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-n")) {
	i++;
	obits=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-headerless")) {
	headerless=1;
      } else if (help_required(argv[1])) {
	decimate_help();
	exit(0);
      } else if (file_exists(argv[i])) {
	strcpy(inpfile,argv[i]);
	input=open_file(inpfile,"rb");
      } else {
	decimate_help();
	sprintf(string,"unknown argument (%s) passed to decimate",argv[i]);
	error_message(string);
      }
      i++;
    }
  }

  /* read in the header to establish what the input data are... */
  if ((headersize=read_header(input))) {
    if ( (nsamp > 0) && !strings_equal(inpfile,"stdin") ) {
      naddt=nsamples(inpfile,headersize,nbits,nifs,nchans)/nsamp;
      if (naddt%2) naddt--;
    }
    switch (data_type) {
    case 1:
      break;
    case 2:
      nchans=1;
      break;
    default:
      error_message("input data to decimate is not in filterbank format");
      break;
    }
    /* check number of time samples to add */
    if (naddt <= 1) naddt=1;
    /*if (naddt%2) error_message("time decimation must be a power of 2");*/
    /* check number of frequency channels to add (integer multiple) */
    if (naddc<1) naddc=nchans;
    nc=nchans/naddc;
    if ( (nc*naddc) != nchans ) 
      error_message("nchans must be integer multiple of decimation factor");
    if (obits == 0) obits=nbits;
    if (obits==1) error_message("output of 1-bit data will result in vastly reduced S/N!\nselect a higher output bit size with the -n option");
    /* all ok - broadcast the new header */
    if (!headerless) decimate_header();
  } else {
    error_message("input data file is of unknown origin!!!");
  }
  
  /* finally decimate and output the data */
  decimate_data(input,output);

  exit(0);
}
