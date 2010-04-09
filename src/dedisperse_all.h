#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern FILE *input, *output, *logfile;
extern char  inpfile[180], outfile[180], ignfile[180];

/* global variables describing the operating mode */
extern int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off;
extern double refrf,userdm,fcorrect;
extern float clipvalue,jyfactor,jyf1,jyf2;

/* global variables describing the data */
#include "header.h"

/* list of subroutines and functions */

#include "sigproc.h"



