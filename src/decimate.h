#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sigproc.h"
#include "header.h"
int nsamp,naddc,naddt,headerless,obits;
char inpfile[80], outfile[80];
FILE *input, *output, *logfile;
