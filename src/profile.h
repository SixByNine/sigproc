#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *input;
double period;
long int np;
int nbins,site;
float *profile, pmin, pmax, prng;
/* global variables describing the data */
#include "header.h"

/* list of subroutines and functions */
#include "sigproc.h"

