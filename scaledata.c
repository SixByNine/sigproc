/* routines for scaling and descaling arrays */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void float2int(float *f, int n, int b, float min, float max, int *i) /*includefile*/ 
{
  int j;
  float delta,imax;
  imax=pow(2.0,(double) b)-1.0;
  delta=max-min;
  for (j=0; j<n; j++) {
    f[j] = (f[j]>max) ? (max): f[j];
    i[j] = (f[j]<min) ? (0.0): (int) rint((f[j]-min)*imax/delta);
  }
}

void int2float(int *i, int n, int b, float min, float max, float *f) /*includefile*/ 
{
  int j;
  float delta,imax;
  imax=pow(2.0,(double) b)-1.0;
  delta=max-min;
  for (j=0; j<n; j++) f[j] = min+((float) i[j])*delta/imax;
}

void float2four(float *f, int n, float min, float max, unsigned char *c) /*includefile*/
{
  int *i,j;
  i=(int *) malloc(sizeof(int)*n);
  float2int(f,n,4,min,max,i);
  for (j=0;j<n;j++) c[j]=0;
  for (j=0;j<n-1;j+=2) c[j/2]=charof2ints(i[j],i[j+1]);
  free(i);
}

void float2char(float *f, int n, float min, float max, unsigned char *c) /*includefile*/
{
  int *i,j;
  i=(int *) malloc(sizeof(int)*n);
  float2int(f,n,8,min,max,i);
  for (j=0;j<n;j++) c[j]=(unsigned char) i[j];
  free(i);
}

void float2short(float *f, int n, float min, float max, unsigned short *s) /*includefile*/
{
  int *i,j;
  i=(int *) malloc(sizeof(int)*n);
  float2int(f,n,16,min,max,i);
  for (j=0;j<n;j++) s[j]=(unsigned short) i[j];
  free(i);
}

