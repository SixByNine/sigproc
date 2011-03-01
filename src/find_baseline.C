/***
 * Normalize method, with threshold overload.
 * Baseline Function using exponential smoothing and threshold overload.
 * R Neil, 27,Jan,2011
***/
#include "find_baseline.h"
#include <limits>
#define MAX_BL_ITTR 5
/**
 * Normalize method, returns a normalized timeseries
 * Overloaded Normalize method takes a threshold and skips all values above that 
 * threshold in sigma's.
**/
void 
normalise(int n, float * d){
    int i;
    double sum = 0.0;
    double sumSq = 0.0;
    double mean = 0.0;
    double meanSq = 0.0;
    double sigma = 0.0;
    /* 2 pass method one pass to calculate mean, second to normalize */
    for(i=0; i<n; ++i) {
        sum+=d[i];
        sumSq+=d[i]*d[i];
    }
    mean=sum/(n*1.0);
    meanSq=sumSq/(n*1.0);
    sigma=sqrt(meanSq-(mean*mean));
    for (i=0; i<n; ++i) { d[i]=(d[i]-mean)/sigma; }
    printf("Normalized: %d samples, mean %f, sigma %f\n",n,mean,sigma);
}
/**
 * Overloaded Normalize method takes threshold.
**/
void 
normalise(int n, float * d, float threshold){
    int i;
    int nsum=0;
    double sum = 0.0;
    double sumSq = 0.0;
    double mean = 0.0;
    double meanSq = 0.0;
    double sigma = 0.0;
    /* Pre-normalize to allow reasonable threshold cuttoffs */
    normalise(n,d);
    /* compute sigma without threshold spikes */
    for(i=0; i<n; ++i) {
        if (fabs(d[i])<threshold){
            sum+=d[i];
	    sumSq+=d[i]*d[i];
	    ++nsum;
        }
    }
    mean = sum/(1.0*nsum);
    meanSq = sumSq/(1.0*nsum);
    sigma=sqrt(meanSq-(mean*mean));
    for (i=0; i<n; ++i) { d[i]=(d[i]-mean)/sigma; }
    printf("Normalized: %d samples, %d skipped, mean %f, sigma %f\n",n,nsum,mean,sigma);
}




float find_baseline_i(int ndat, float * dat, float smooth_nsamp, float threshold) {

    int i;
    float frac;
    double cSum = 0.0;
    /* set the smoothing factor */
    double sf = 1.0/smooth_nsamp;
    /* normalize the timeseries */
    normalise(ndat, dat,threshold);

    printf("sm=%f\tsf=%f\n",smooth_nsamp,sf);
    float* smooth=(float*)malloc(sizeof(float)*ndat);

    cSum=dat[0];
    for(i=1; i<ndat; i++) {
	if (dat[i]<threshold) {
        	cSum = (sf*dat[i])+((1-sf)*cSum);
	}
        smooth[i] = cSum;
    }
    cSum=dat[ndat-1];
    for(i=ndat-2; i>=0; i--) {
	if (dat[i]<threshold) {
        	cSum = (sf*dat[i])+((1-sf)*cSum);
	}
	frac=(float)i/ndat;
        smooth[i] = frac*smooth[i]+(1.0-frac)*cSum;
    }

    for(i=0; i<ndat; i++) {
	    dat[i]-=smooth[i];
    }


    normalise(ndat, dat,threshold);

    float ssq=0;
    for(i=0; i<ndat; i++) {
	    ssq+=smooth[i]*smooth[i];
    }

    free(smooth);
    return ssq/(float)ndat;
}




/**
 * baseline method using exponential smoothing, smoothing factor set at 0.02 default.
 * ndat is the size of the timeseries
 * dat is the timeseries
 * overloaded method takes threshold which will cut out spikes above specified 
 * threshold in sigma's.
**/
void
find_baseline( int ndat, float * dat,float smooth_nsamp) {
	find_baseline(ndat,dat,smooth_nsamp,std::numeric_limits<float>::max());
}

char sw=1;
/** 
 * Overloaded baseline method takes threshold
**/
void
find_baseline( int ndat, float * dat, float smooth_nsamp, float threshold) {
	int count=0;
/*
 * TEST code that put in a huge cubic!
 	if(sw){
		for(int i=0; i<ndat; i++) {
			float x= (float)i/100000.0;
			dat[i]+=x - 100*x*x + 5*x*x*x;
		}
		sw=0;
	} else{*/
		while(count < MAX_BL_ITTR){
			float goodness = find_baseline_i(ndat,dat,smooth_nsamp,threshold);
			printf("%02d goodness = %f\n",count,goodness);
			count++;
			if(goodness < 0.01) break;
		}
/*	}*/
}

