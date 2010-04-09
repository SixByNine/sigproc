/* Set of routines that are commonly used in tempo2 and/or its plugins.
 * These routines are mainly stand-alone functions and exist for 
 * float and double precision variables
 *
 * G. Hobbs: v2, 31 Dec 2008.  Complete rewrite of the routines
 *
 * NOTES: Related toolkits include:
 *  TKspectrum.h: contains routines for spectral estimation
 *  TKfit.h:      contains routines for fitting
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "toolkit.h"


/* ************************************************************** */
/* Routines to convert doubles to floats                          */
/* ************************************************************** */

void TKconvertFloat1(double *x,float *ox,int n)
{
  int i;
  for (i=0;i<n;i++)
    ox[i] = (float)x[i];
}

void TKconvertFloat2(double *x,double *y,float *ox,float *oy,int n)
{
  int i;
  for (i=0;i<n;i++)
    {
      ox[i] = (float)x[i];
      oy[i] = (float)y[i];
    }
}

/* ************************************************************** */
/* Basic statistics                                               */
/* ************************************************************** */

float TKfindMin_f(float *x,int n)
{
  float ret;
  int i;

  ret = x[0];
  for (i=0;i<n;i++)
    {
      if (x[i] < ret) ret = x[i];
    }
  return ret;
}

double TKfindMin_d(double *x,int n)
{
  double ret;
  int i;

  ret = x[0];
  for (i=0;i<n;i++)
    {
      if (x[i] < ret) ret = x[i];
    }
  return ret;
}


double TKsign_d(double a,double b)
{
  if (b >= 0)
    return fabs(a);
  else
    return -fabs(a);
}

float TKfindMedian_f(float *val,int count)   /* Find median in x array */
{
  int i;
  float median = 0.0;
  float newval[count];

  for (i=0;i<count;i++)
    newval[i]=val[i];
  TKsort_f(newval,count);
  if (count%2==0)
    median = (newval[count/2-1]+newval[count/2])/2.0;
  else
    median = newval[(count-1)/2];

  return median;

}

float TKfindRMS_f(float *x,int n)
{
  int i;
  float mean,sdev=0.0;
  mean = TKmean_f(x,n);
  for (i=0;i<n;i++)
    sdev += pow(x[i]-mean,2);
  sdev/=(n-1);
  sdev = sqrt(sdev);
  return sdev;
}

float TKfindRMS_d(double *x,int n)
{
  int i;
  double mean,sdev=0.0;
  mean = TKmean_d(x,n);
  for (i=0;i<n;i++)
    sdev += pow(x[i]-mean,2);
  sdev/=(n-1);
  sdev = sqrt(sdev);
  return sdev;
}

float TKfindMax_f(float *x,int n)
{
  float ret;
  int i;

  ret = x[0];
  for (i=0;i<n;i++)
    {
      if (x[i] > ret) ret = x[i];
    }
  return ret;
}

double TKretMax_d(double a,double b)
{
  if (a>b) return a;
  else return b;
}

double TKretMin_d(double a,double b)
{
  if (a<b) return a;
  else return b;
}

float TKretMax_f(float a,float b)
{
  if (a>b) return a;
  else return b;
}

float TKretMin_f(float a,float b)
{
  if (a<b) return a;
  else return b;
}


int TKretMin_i(int a,int b)
{
  if (a<b) return a;
  else return b;
}

double TKfindMax_d(double *x,int n)
{
  double ret;
  int i;

  ret = x[0];
  for (i=0;i<n;i++)
    {
      if (x[i] > ret) ret = x[i];
    }
  return ret;
}

double TKmean_d(double *x,int n)
{
  double mean=0.0;
  int i;
  for (i=0;i<n;i++)
    mean += x[i];

  return mean/(double)n;
}

float TKmean_f(float *x,int n)
{
  float mean=0.0;
  int i;
  for (i=0;i<n;i++)
    mean += x[i];

  return mean/(float)n;
}

double TKvariance_d(double *x,int n)
{
  double mean;
  double var=0.0;
  double ep=0.0,s;
  int i;

  mean = TKmean_d(x,n);
  for (i=0;i<n;i++)
    {
      s = x[i]-mean;
      ep  += s;      
      var += s*s;
    }
  var = (var-ep*ep/(double)n)/(double)(n-1);
  return var;
}

double TKrange_d(double *x,int n)
{
  return (TKfindMax_d(x,n)-TKfindMin_d(x,n));
}


/* ************************************************************** */
/* Sorting                                                        */
/* ************************************************************** */


void TKsort_f(float *val,int nobs)
{
  int i,changed=0;
  float store;

  do 
    {
      changed=0;
      for (i=0;i<nobs-1;i++)
	{
	  if (val[i] > val[i+1])
	    {
	      store=val[i+1];
	      val[i+1]=val[i];
	      val[i]=store;
	      changed=1;
	    }
	}
    }while (changed==1);
}

void TKsort_d(double *val,int nobs)
{
  int i,changed=0;
  double store;

  do 
    {
      changed=0;
      for (i=0;i<nobs-1;i++)
	{
	  if (val[i] > val[i+1])
	    {
	      store=val[i+1];
	      val[i+1]=val[i];
	      val[i]=store;
	      changed=1;
	    }
	}
    }while (changed==1);
}

void TKsort_2f(float *val,float *val2,int nobs)
{
  int i,changed=0;
  float store;

  do 
    {
      changed=0;
      for (i=0;i<nobs-1;i++)
	{
	  if (val[i] > val[i+1])
	    {
	      store=val2[i+1];
	      val2[i+1]=val2[i];
	      val2[i]=store;
	      store=val[i+1];
	      val[i+1]=val[i];
	      val[i]=store;
	      changed=1;
	    }
	}
    }while (changed==1);
}

/* ************************************************************** */
/* Routines that modify a data array                              */
/* ************************************************************** */

void TKzeromean_d(int n,double *y)
{
    int i;
    double ysum = 0.0;

    for (i = 0; i < n; i++)
	ysum += y[i];
    ysum /= n;
    for (i = 0; i < n; i++)
      y[i] -= ysum;
}

/* ************************************************************** */
/* Random number generation                                       */
/* ************************************************************** */
/* Returns random deviate between 0 and 1  */
double TKranDev(long *seed)
{
  double x;

  if (*seed < 0)
    {
      init_genrand(-(*seed));
      *seed=1;
    }
  x = genrand_real1();
  return x;
}


/* Returns Gaussian random deviates */
/* Based on the GNU scientific library code: gsl_ran_gaussian
 * This is a polar (Box-Mueller) method */

double TKgaussDev(long *seed)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      x = -1 + 2 * TKranDev(seed);
      y = -1 + 2 * TKranDev(seed);

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return y * sqrt (-2.0 * log (r2) / r2);
}


/* Sets up a random number seed based on the clock (based on the second counter) */
long TKsetSeed()
{
  long seed;
  seed = -time(NULL);
  return seed;
}

/* Random number routines -- should call the above functions */

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* Period parameters */  
#define RAND_N 624
#define RAND_M 397

static unsigned long rand_mt[RAND_N]; /* the array for the state vector  */
static int rand_mti=RAND_N+1; /* rand_mti==RAND_N+1 means rand_mt[RAND_N] is not initialized */

/* initializes rand_mt[RAND_N] with a seed */
void init_genrand(unsigned long s)
{
    rand_mt[0]= s & 0xffffffffUL;
    for (rand_mti=1; rand_mti<RAND_N; rand_mti++) {
        rand_mt[rand_mti] = 
	    (1812433253UL * (rand_mt[rand_mti-1] ^ (rand_mt[rand_mti-1] >> 30)) + rand_mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array rand_mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        rand_mt[rand_mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}


/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, 0x9908b0dfUL};
    /* mag01[x] = x * 0x9908b0dfUL  for x=0,1 */

    if (rand_mti >= RAND_N) { /* generate RAND_N words at one time */
        int kk;
        if (rand_mti == RAND_N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<RAND_N-RAND_M;kk++) {
            y = (rand_mt[kk]&0x80000000UL)|(rand_mt[kk+1]&0x7fffffffUL);
            rand_mt[kk] = rand_mt[kk+RAND_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<RAND_N-1;kk++) {
            y = (rand_mt[kk]&0x80000000UL)|(rand_mt[kk+1]&0x7fffffffUL);
            rand_mt[kk] = rand_mt[kk+(RAND_M-RAND_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }

        y = (rand_mt[RAND_N-1]&0x80000000UL)|(rand_mt[0]&0x7fffffffUL);
        rand_mt[RAND_N-1] = rand_mt[RAND_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        rand_mti = 0;
    }
  
    y = rand_mt[rand_mti++];
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}


/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
  double x;
  x = genrand_int32()*(1.0/4294967295.0);
  return x;
  //    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}


/**
 * Finds the mean and stddev of the given numbers, using a 'better' method.
 *
 * vals is 'l' floats.
 *
 * 'm' is the output mean, 'r' is the output stddev
 *
 */
void meanrms(float* vals, int l, float* m ,float* r){
        int navg;
        float old_mean,old_S,mean,S,rms;
	float sum;

        navg = 1;
	sum=vals[0];
        old_mean = vals[0];
        old_S = 0.0;

        for(int i=1; i < l; i++){
                if(vals[i]!=0){
                        navg++;
			sum += vals[i];
                        mean = sum/navg;
                        S = old_S + (vals[i] - old_mean)*(vals[i] - mean);

                        old_mean = mean;
                        old_S = S;
                }
        }
        if(navg>0){
                rms=sqrt(S/(navg-1));
        } else {
                rms=0;
                mean=0;
        }

        (*r)=rms;
        (*m)=mean;
        return;

}


int quicksort_index_partition(float* array, int* index, int top, int bottom);
void quicksort_index_wk(float* array, int* index, int top, int bottom){
        int middle;
        if (top < bottom)
        {
                middle = quicksort_index_partition(array,index, top, bottom);
                quicksort_index_wk(array,index, top, middle);
                quicksort_index_wk(array,index, middle+1, bottom);
        }
        return;
}

void quicksort_index(float* array, int* index, int npts){
	quicksort_index_wk(array,index,1,npts-1);
}

int quicksort_index_partition(float* array, int* index, int top, int bottom){
        float x = array[index[top]];
        int topidx = top - 1;
        int botidx = bottom + 1;
        int swap;
        do
        {
                do
                {
                        botidx --;
                }while (x > array[index[botidx]]);

                do
                {
                        topidx++;
                } while (x < array[index[topidx]]);

                if (topidx < botidx)
                {
                        swap = index[topidx];
                        index[topidx] = index[botidx];
                        index[botidx] = swap;
                }
        }while (topidx < botidx);
        return botidx;
}


int quicksort_inplace_partition(float* array, int top, int bottom);
void quicksort_inplace_wk(float* array, int top, int bottom){
        int middle;
        if (top < bottom)
        {
                middle = quicksort_inplace_partition(array, top, bottom);
                quicksort_inplace_wk(array, top, middle);
                quicksort_inplace_wk(array, middle+1, bottom);
        }
        return;
}

void quicksort_inplace(float* array, int npts){
	quicksort_inplace_wk(array,1,npts-1);
}

int quicksort_inplace_partition(float* array, int top, int bottom){
        float x = array[top];
        int topidx = top - 1;
        int botidx = bottom + 1;
        float swap;
        do
        {
                do
                {
                        botidx --;
                }while (x >array[botidx]);

                do
                {
                        topidx++;
                } while (x <array[topidx]);

                if (topidx < botidx)
                {
                        swap = array[topidx];
                        array[topidx] = array[botidx];
                        array[botidx] = swap;
                }
        }while (topidx < botidx);
        return botidx;
}
