#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

extern FILE *input, *output, *logfile;
extern char  inpfile[80], outfile[80], ignfile[80];

/* global variables describing the operating mode */
extern int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off;
double refrf,userdm,fcorrect;
float clipvalue,jyfactor,jyf1,jyf2;

/* global variables describing the data */
#include "header.h"

/* list of subroutines and functions */
#include "sigproc.h"
