/* find.h - include file for find program and related routines */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* input and output files and logfile (filterbank.monitor) */
FILE *input, *output, *logfile;
char  inpfile[80], outfile[80];

/* global variables describing the data */
#include "header.h"
#include "sigproc.h"
