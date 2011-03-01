/***
 * Normalize method, with threshold overload.
 * Baseline Function using exponential smoothing and threshold overload.
 * R Neil, 27,Jan,2011
***/
#ifndef __BASELINE_H__
#define __BASELINE_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void normalise(int n, float * d);

void normalise(int n, float * d, float threshold);

void find_baseline(int ndat, float * dat,float smooth_nsamp);

void find_baseline(int ndat, float * dat, float smooth_nsamp, float threshold);
#endif
