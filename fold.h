#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *input, *output, *logfile;
char  inpfile[80], outfile[80];

/* global variables describing the operating mode */
int binary, ascii, asciipol, totalpower, stream, headerless, 
  baseline, multiple, psrfits;

char polyco_file[80];

double folding_period, tsamp_user, dump_time, phase_start, phase_finish;
double acceleration,tobs;
int accumulate, nbins, npulses;
float *folded_profiles, jyfactor, userbase;

/* global variables describing the data */
#include "header.h"
double time_offset, skip_time, read_time;

/* list of subroutines and functions */
#include "sigproc.h"

#ifdef PSRFITS
/* Define CFITS file pointer - global- defined in fitsio.h */
#include "fitsio.h"
fitsfile *fits;
#endif
