/* filterbank.h - include file for filterbank and related routines */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* input and output files and logfile (filterbank.monitor) */
FILE *input, *output, *logfile;
char  inpfile[80], outfile[80];

/* global variables describing the data */
#include "header.h"
double time_offset;

/* global variables describing the operating mode */
float start_time, final_time, clip_threshold;

int obits, sumifs, headerless, headerfile, swapout, invert_band;
int compute_spectra, do_vanvleck, hanning, hamming, zerolagdump;
int headeronly;
char ifstream[8];

/* library of subroutines and functions */
#include "sigproc.h"
