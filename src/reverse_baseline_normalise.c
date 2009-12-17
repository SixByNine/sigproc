#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "dedisperse.h"

/*
 *Globals:
 *baselineCurrentMeans
 *baselineCurrentVariances
 *baselineValidStart
 *baselineValieEnd
 *baselineFile
 */


void reverseBaselineNormalise(float* data, unsigned long long sample0, int nsamps, int nchans){

	FILE *file;
	unsigned long long sample;
	int isamp,ichan;
	float empty;
	unsigned long long sampstart, sampend;
	char string[256];
	double inmean;

	switch(nbits){
		case 1:
			inmean=0.5;
			break;
		case 2:
			inmean=1.5;
			break;
		case 4:
			inmean=7.5;
			break;
		case 8:
			inmean=127.5;
			break;
		default:
			fprintf(stderr,"Can't use reverse baseline for nbits!=1,2,4,8");
			return;
	}


//	fprintf(stderr,"%lld (%lld -> %lld) %d %d\n",sample0,baselineValidStart,baselineValidEnd,nsamps,nchans);

	sample = sample0;
	if(baselineValidStart == baselineValidEnd){
		baselineCurrentMeans = malloc(sizeof(float)*nchans);
		baselineCurrentVariances = malloc(sizeof(float)*nchans);
	}

	for(isamp=0; isamp < nsamps; isamp++){
		if(sample < baselineValidStart || sample >= baselineValidEnd ){
			// we are outside of this baseline valid range...
			
			// Look for a valid range in the file.
			fprintf(stderr,"Out of valid range! %d\n",sample);
			file = fopen(baselineFile,"r");
			while(!feof(file)){
				fscanf(file,"%s\n",string);
				if(strcmp(string,"#START#")==0){
					fscanf(file,"#%lld %lld\n",&sampstart, &sampend);
//					fprintf(stderr,"trying: (%lld -> %lld)\n",sampstart,sampend);
					if(sample < sampend && sample >= sampstart){
//						fprintf(stderr,"found: (%lld -> %lld)\n",sampstart,sampend);
						baselineValidStart = sampstart;
						baselineValidEnd = sampend;

						for(ichan = 0; ichan < nchans; ichan++){
							fscanf(file,"%f %f %f\n",&empty,baselineCurrentMeans+ichan,baselineCurrentVariances+ichan);
						}
						break;
					}
				}
			}
			fclose(file);
			

		}
		for(ichan = 0; ichan < nchans; ichan++){
//			if(sample==0||sample==1)fprintf(stderr,"%d %f %f %f\n",ichan,data[ichan + isamp*nchans],baselineCurrentMeans[ichan],baselineCurrentVariances[ichan]);

			data[ichan + isamp*nchans] = (data[ichan + isamp*nchans] - inmean);
			data[ichan + isamp*nchans] = data[ichan + isamp*nchans] * sqrt(baselineCurrentVariances[ichan]) + baselineCurrentMeans[ichan];
		}
		sample++;
	}



}
