/*
  FILTERBANK  - convert raw data from a pulsar machine into "filterbank data"
  a stream of samples: [s1,c1] [s1,c2]...[s1,cn] [s2,c1]...[s2,cn]..... etc
  where the indices s1, s2.... are individual sample readouts and c1, c2.... 
  are the individual channels. For each readout, there are "n" channels.
  Before writing out a stream of samples, filterbank sends a header string
  which is intercepted by other programs (dedisperse, fold etc) downstream.
  (N.B. The header can be omitted using the -headerfile command-line option)
*/
#include "filterbank.h"

int wapp_isalfa;

main (int argc, char *argv[])
{
  int i,nfiles,fileidx,fileidx_start,inputdata,opened=0;
  char message[80];
  int numsamps;
  int bpp_headersize = 32768;
  float f;
  int sample_end,fileidx_final,byte_end,bytefinal;
  double sample_final,scantime;
  int data_size,sample_skip,bytestart,ns,blocksize;
  double sample_start;
  
  /* check number of command-line arguments */
  if (argc<2) {
    filterbank_help();
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }
 
  /* print help if necessary */
  if (help_required(argv[1])) {
    filterbank_help();
    exit(0);
  }
 

  /* set up default global variables */
  hanning=hamming=zerolagdump=swapout=sumifs=headerless=headerfile=0;
  wapp_isalfa=invert_band=clip_threshold=headeronly=0;
  time_offset=start_time=final_time=0.0;
  obits=-1;
  do_vanvleck=compute_spectra=1;
  strcpy(ifstream,"XXXX");

  /* work out how many files are on the command line */
  i=1;
  nfiles=0;
  while(file_exists(argv[i])) {
	nfiles++;
	i++;
  }
  if (!nfiles) error_message("no input files supplied on command line!");
  fileidx=1;

  /* now parse any remaining command line parameters */
  if (argc>nfiles) {
    i=nfiles+1;
    while (i<argc) {
      if (strings_equal(argv[i],"-o")) {
	/* get and open file for output */
	strcpy(outfile,argv[++i]);
	output=fopen(outfile,"wb");
	opened=1;
      } else if (strings_equal(argv[i],"-c")) {
	/* get clip threshold (sigma) */
	clip_threshold=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-s")) {
	/* get starting time (s) */
	start_time=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-r")) {
	/* get time to read (s) this is adjusted below if skipping */
	final_time=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-n")) {
	/* output number of bits per sample to write */
	obits=atoi(argv[++i]);
      } else if (strings_equal(argv[i],"-dt")) {
	/* add a time offset in seconds to tstart */
	time_offset=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-i")) {
	/* flag IF stream to write */
	i++;
	if (atoi(argv[i])<1 || atoi(argv[i])>4) {
	  error_message("IFstream must lie between 1 and 4");
	}
	ifstream[atoi(argv[i])-1]='Y';
      } else if (strings_equal(argv[i],"-swapout")) {
	/* perform byte swapping on all output data */
	swapout=1;
      } else if (strings_equal(argv[i],"-floats")) {
	/* write data as floating point numbers */
	obits=32;
      } else if (strings_equal(argv[i],"-sumifs")) {
	/* sum IFs if necessary */
	sumifs=1;
      } else if (strings_equal(argv[i],"-zerolag")) {
	/* zerolagdump used for correlators e.g. WAPP */
	zerolagdump=1;
	obits=32;
      } else if (strings_equal(argv[i],"-rawcfs")) {
	/* write correlation functions only */
	compute_spectra=do_vanvleck=0;
      } else if (strings_equal(argv[i],"-corcfs")) {
	/* write corrected correlation functions */
	compute_spectra=0;
	do_vanvleck=1;
      } else if (strings_equal(argv[i],"-novanvleck")) {
	/* don't apply van vleck correction */
	do_vanvleck=0;
      } else if (strings_equal(argv[i],"-invert")) {
	/* invert the band after FFT */
	invert_band=1;
      } else if (strings_equal(argv[i],"-hamming")) {
	/* Hamming smoothing */
	hamming=1;
	hanning=0;
      } else if (strings_equal(argv[i],"-hanning")) {
	/* Hanning smoothing */
	hanning=1;
	hamming=0;
      } else if (strings_equal(argv[i],"-headerfile")) {
	/* no binary headers but write data to "head" file */
	headerless=headerfile=1;
      } else if (strings_equal(argv[i],"-headeronly")) {
	/* only binary header written */
	headeronly=1;
      } else {
	/* unknown argument passed down - stop! */
	filterbank_help();
	sprintf(message,"unknown argument (%s) passed to filterbank.",argv[i]);
	error_message(message);
      }
      i++;
    }
  }

  /* if no IF streams selected, flag all */
  if (strings_equal(ifstream,"XXXX")) strcpy(ifstream,"YYYY");

  /* adjust finish time if necessary */
  if (final_time > 0.0) final_time+=start_time;

  if (!opened) {
    /* no output file selected, use standard output */
    output=stdout;
    strcpy(outfile,"stdout");
  }


  fileidx_start = fileidx;
  /* For BCPM data:  Test to see if the first input file is BCPM type*/
  /* If the file type is not BCPM then close the file and continue normally */
  /* open up the first input file*/
  strcpy(inpfile,argv[fileidx]);
  input=open_file(inpfile,"rb");
  inputdata=typeof_inputdata(input,inpfile);
  
  /*Establish whether or not the file data type is BCPM*/
  if (inputdata == 5) {
    /* If the start time is not zero, find the file to start with */
    data_size = sizeof_file(inpfile)-bpp_headersize;
    scantime = (double)(data_size/((double)nchans/2)*(double)tsamp);
    ns = 512;
    blocksize = ns*(double)nchans/2;
    
    if(start_time) {
      /* Calculate the no. of time samples */
      /* Calculate the file corresponding to the start time */
      fileidx_start = ceil(start_time/scantime);
      sample_skip = floor(start_time/tsamp);
      sample_start = sample_skip - (double)((fileidx_start-1)*(scantime/(double)tsamp));
      bytestart = (sample_start*(double)nchans)/2.;

      if(bytestart<blocksize/2){
	fprintf(stderr,"bytestart is less than blocksize/2\n");
	bytestart+=data_size-blocksize/2;
	fileidx_start-=1;
	sample_skip -= ns;
	sample_start = sample_skip - (double)((fileidx_start-1)*(scantime/(double)tsamp));
      }

      fprintf(stderr,"Starting Time:  \n");
      fprintf(stderr,"          Start File #              %2d\n",fileidx_start);
      fprintf(stderr,"          Start Sample #  %12.3f\n",sample_start);
      fprintf(stderr,"          Start Byte #    %12d\n",bytestart);
      fprintf(stderr,"          Start Time (s)  %12.3f\n",sample_skip*(double)tsamp-(fileidx_start-1)*scantime);
      fprintf(stderr,"\nAdvancing start time by   %12.5f\n\n",sample_skip*(double)tsamp);


      fileidx = fileidx_start;
    }
    if(final_time) {
      sample_end = ceil(final_time/(double)tsamp);
      fileidx_final = ceil(final_time/scantime);
      sample_final = sample_end-(double)((fileidx_final-1)*scantime/tsamp);
      byte_end = sample_end*(double)nchans/2;
      bytefinal = (double)sample_final*(double)nchans/2;
      nfiles = fileidx_final;


    }
  }
  fclose(input);

  /* main loop around input files */
  while (fileidx <= nfiles) {

    /* open up input file */
    strcpy(inpfile,argv[fileidx]);
    input=open_file(inpfile,"rb");

    /* open the raw data file and establish its origin and header pars */
    if (!(inputdata=typeof_inputdata(input,inpfile)))
    error_message("input data file is of unknown origin!!!");
    /* check for timing data files - not processed here! */
    switch (inputdata) {
    case 2:
      error_message("input data (PSPM timing mode) not read by this program!");
      break;
    case 4:
      error_message("input data (WAPP timing mode) not read by this program!");
      break;
    case 6:
      error_message("input data (BPP timing mode) not read by this program!");
      break;
    }

    if (fileidx == fileidx_start) {
  	/* add on a time offset in seconds to the start time */
  	tstart+=time_offset/86400.0;
  	/* broadcast the header */
  	if (!wapp_isalfa) filterbank_header(output);
	if (headeronly) exit(0);
    }

    /* now actually convert the raw data into filterbank format */
    switch (inputdata) {
    case 1: 
      /* PSPM search-mode data */
      pspm2fb(input,output);
      break;
    case 3:
      /* WAPP fast-sampled data */
      wapp2fb(input,output);
      break;
    case 5:
      /* BPP fast-sampled data */
      bpp2fb(input,output);
      break;
    case 7:
      /* AOFTM fast-sampled data */
      aoftm2fb(inpfile,output);
      break;
    case 8:
      /* OOTY filterbank data */
      ooty2fb(input,output);
      break;
    case 9:
      /* SCAMP filterbank data */
      scamp2fb(input,output);
      break;
    case 10:
      /* GMRT filterbank data */
      gmrt2fb(input,output);
      break;
    case 11:
      /* Effelsberg pulsar2000 data */
      pulsar2k2fb(input,output);
      break;
    }
    fileidx++;
    fclose(input);
  }

  /* all done, update log, close all files and exit normally */
  update_log("finished");
  close_log();
  /*fclose(output);*/
  exit(0);
}

