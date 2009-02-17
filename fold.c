/*
  FOLD  - fold filterbank channels or time series data to produce integrated
  pulse profiles or single pulses.
*/
#include "fold.h"

main (int argc, char *argv[])
{
  /* local variables */
  double pfactor,newmjd=0.0;
  float sefd;
  int i,opened_input=0,opened_output=0,headersize=0;
  char string[80];

  /* set up default globals */
  baseline=ascii=multiple=1;
  npuls=binary=totalpower=accumulate=0;
  time_offset=acceleration=skip_time=read_time=0.0;
  asciipol=psrfits=stream=headerless=npulses=0;
  phase_start=folding_period=dump_time=tsamp_user=0.0;
  phase_finish=pfactor=1.0;
  jyfactor=sefd=userbase=0.0;
  nbins=0; /* this will get set in the folding routine if not set by user */
  strcpy(polyco_file,"");

  /* check the command line parameters */
  i=1;
  while (i<argc) {
    print_version(argv[0],argv[1]);
    if (strings_equal(argv[i],"-o")) {
      /* get and open file for output */
      i++;
      strcpy(outfile,argv[i]);
      output=fopen(outfile,"wb");
      opened_output=1;
    } else if (strings_equal(argv[i],"-m")) {
      multiple=atoi(argv[++i]);
    } else if (strings_equal(argv[i],"-p")) {
      /* get folding period */
      i++;
      if (file_exists(argv[i])) {
	strcpy(polyco_file,argv[i]);
	folding_period=-1.0;
      } else {
	folding_period=atof(argv[i]);
      }
    } else if (strings_equal(argv[i],"-dt")) {
      /* add a time offset in seconds to tstart */
      time_offset=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-mjd")) {
      /* change the start time completely! */
      newmjd=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-sk")) {
      /* skip the first skip_time seconds before folding */
      skip_time=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-re")) {
      /* read and fold only read_time seconds of data */
      read_time=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-a")) {
      /* get acceleration for folding */
      acceleration=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-d")) {
      /* get dumptime or number of pulses for subintegrations */
      i++;
      if (strcspn(".",argv[i])) {
	npulses=atoi(argv[i]);
      } else {
	dump_time=atof(argv[i]);
      }
    } else if (strings_equal(argv[i],"-t")) {
      /* get user-supplied sampling time */
      i++;
      tsamp_user=atof(argv[i]);
    } else if (strings_equal(argv[i],"-j")) {
      /* get user-supplied Jansky calibration factor */
      jyfactor=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-s")) {
      /* get user-supplied SEFD */
      sefd=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-b")) {
      /* get user-supplied baseline */
      baseline=0;
      userbase=atof(argv[++i]);
    } else if (strings_equal(argv[i],"-f")) {
      /* get period multiplication factor */
      i++;
      pfactor=atof(argv[i]);
    } else if (strings_equal(argv[i],"-l")) {
      /* get leading phase of pulse */
      i++;
      phase_start=atof(argv[i]);
      if ( (phase_start < 0.0) || (phase_start > 1.0) ) 
	error_message("start pulse phase out of range!");
    } else if (strings_equal(argv[i],"-r")) {
      /* get trailing phase of pulse */
      i++;
      phase_finish=atof(argv[i]);
      if ( (phase_finish < 0.0) || (phase_finish > 1.0) ) 
	error_message("final pulse phase out of range!");
    } else if (strings_equal(argv[i],"-n")) {
      /* get number of bins */
      i++;
      nbins=atoi(argv[i]);
    } else if (strings_equal(argv[i],"-ascii")) {
      /* write data as ASCII numbers */
      ascii=1;
    } else if (strings_equal(argv[i],"-totalpower")) {
      /* sum polarizations 1+2 before writing */
      totalpower=1;
    } else if (strings_equal(argv[i],"-epn")) {
      /* write data in EPN format */
      ascii=0;
    } else if (strings_equal(argv[i],"-bin")) {
      /* write data in SIGPROC binary format */
      binary=1;
    } else if (strings_equal(argv[i],"-acc")) {
      /* write out accumulated pulse profiles in subints */
      accumulate=1;
    } else if (strings_equal(argv[i],"-asciipol")) {
      /* write data as ASCII numbers for Jim's polarization code */
      asciipol=1;
    } else if (strings_equal(argv[i],"-psrfits")) {
      /* write data in PSRFITS format */
      ascii=0;
      psrfits=1;
#ifndef PSRFITS
      error_message("-psrfits option not supported in this compilation...\nConsult the SIGPROC manual for further information about PSRFITS.");
#endif
    } else if (strings_equal(argv[i],"-stream")) {
      /* write data as ASCII streams */
      stream=1;
    } else if (strings_equal(argv[i],"-sub")) {
      /* shorthand for -nobaseline -stream -d x */
      stream=1;
      baseline=0;
      i++;
      if (strcspn(".",argv[i])) {
	npulses=atoi(argv[i]);
      } else {
	dump_time=atof(argv[i]);
      }
    } else if (strings_equal(argv[i],"-nobaseline")) {
      /* processing correlation functions so don't subtract baseline */
      baseline=0;
    } else if (file_exists(argv[i])) {
      /* get and open file for input */
      strcpy(inpfile,argv[i]);
      input=open_file(inpfile,"rb");
      opened_input=1;
    } else if (help_required(argv[i])) {
      fold_help();
      exit(0);
    } else {
	/* unknown argument passed down - stop! */
	fold_help();
	sprintf(string,"unknown argument (%s) passed to %s",argv[i],argv[0]);
	error_message(string);
    }
    i++;
  }

  /* get appropriate calibration factor from SEFD and baseline */
  if (sefd != 0.0 && userbase != 0.0) jyfactor=sefd/userbase;

  /* multiply folding period by user-supplied factor */
  if (folding_period != -1.0) folding_period*=pfactor;

  /* check start and end phase of pulse */
  if (phase_start >= phase_finish) 
    error_message("silly pulse phases selected!");

  /* check npulses versus dump_time */
  if (npulses < 0) error_message("npulses < 0!");
  if ((npulses > 0) && (dump_time > 0.0)) 
    error_message("can't have npulses AND dumptime defined!");

  /* check for folding period still set to zero - if so, look for polyco.dat */
  if (folding_period == 0.0) {
    strcpy(polyco_file,"polyco.dat");
    if (file_exists(polyco_file)) {
      folding_period=-1.0;
    } else {
      error_message("folding period not specified and no polyco.dat found!");
    }
  }

  if (!opened_input) {
    /* no input file selected, use standard input */
    input=stdin;
    strcpy(inpfile,"stdin");
  }

  /* read in the header parameters from the input stream */
  if (!(headersize=read_header(input))) 
    error_message("could not read header parameters!");

  if (acceleration != 0.0) {
    tobs=tsamp*(double)nsamples(inpfile,headersize,nbits,nifs,nchans);
    if (tobs <= 0.0) error_message("could not get sensible observation time");
  }

  /* override the header */
  if (newmjd!=0.0) tstart=newmjd;

  if (!opened_output) {
    /* no output file selected, use standard output */
    output=stdout;
    strcpy(outfile,"stdout");
  }

  /* open the raw data file and establish its origin and header pars */
  switch(data_type) {
  case 1: 
  case 2:
  case 6:
    open_log("fold.monitor");
    folded_profiles=fold_data();
    break;
  default:
    error_message("input data is of unknown origin!!!");
  }

  if ((npulses==0.0) && (dump_time==0.0))
     write_profiles(folded_profiles,nbins,nchans,nifs,output);
  if (stream) fprintf(output,"#DONE\n");

  /* all done, update and close logfile */
  update_log("finished");
  close_log();
  i=0;
#ifdef PSRFITS
  if (psrfits) fits_close_file(fits,&i);
#endif
  exit(0);
}
