#ifndef DEDISPERSE_H
#define DEDISPERSE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern FILE *input, *output, *logfile;
extern char  inpfile[80], outfile[80], ignfile[80];

void subzero(float* block, int nSamplesRead);
/*
 * block is a float array of length nSamplesRead*nchans
 * clipAbove/Below is a float array of length nchans
 * clipmode is a flag set: bit 1 = low, bit 2 = high.
 */
void inputClip(float* block, int nSamplesRead,float* clipAbove, float* clipBelow,int clipmode);

void reverseBaselineNormalise(float* data, unsigned long long sample0, int nsamps, int nchans);

char doBaselineFile;
float *baselineCurrentMeans;
float *baselineCurrentVariances;
unsigned long long baselineValidStart;
unsigned long long baselineValidEnd;
char baselineFile[80];

/* global variables describing the operating mode */
int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off,doSubZero,doInputClip;
double refrf,userdm,fcorrect;
float clipvalue,jyfactor,jyf1,jyf2,*inputClipHigh,*inputClipLow;

/* global variables describing the data */
#include "header.h"

/* list of subroutines and functions */
#include "sigproc.h"
#endif
