#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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

unsigned char charof8ints(int* i){
	int posn;
	unsigned char retVal=0;
	for(posn = 0; posn < 8; posn++){
		if(i[posn] > 0) retVal = retVal | (1 << posn);
	}
	return retVal;

}


void float2one(float *f, int n, float min, float max, unsigned char *c) /*includefile*/
{
	int *i,j;
	i=(int *) malloc(sizeof(int)*n);
	float2int(f,n,4,min,max,i);
	for (j=0;j<n;j++) c[j]=0;
	for (j=0;j<n-7;j+=8) c[j/8]=charof8ints(i+j);
	free(i);
}

unsigned char charof4ints (int i, int j, int k, int l) /* includefile */
{
	return ( (l<<6)+(k<<4)+(j<<2)+i );
}

void float2two(float *f, int n, float min, float max, unsigned char *c) /*includefile*/
{
	int *i,j;
	i=(int *) malloc(sizeof(int)*n);
	float2int(f,n,2,min,max,i);
	for (j=0;j<n;j++) c[j]=0;
	for (j=0;j<n-3;j+=4) c[j/4]=charof4ints(i[j],i[j+1],i[j+2],i[j+3]);
	free(i);
}


