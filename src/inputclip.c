#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "header.h"
#include <stdio.h>
//#define INPUTCLIP_TEST_MODE

/**
 * Clip the data as it is read.
 * This is designed to solve the Nancay survey problem of -ve spikes.
 * M Keith 2007
 */
void inputClip(float* block, int nSamplesRead, float* clipBelow, float* clipAbove,int clipmode){

	int i,j,k;
	float sum;
	int** histogram;
	int histogramSize,intValue,counter,medianLimit;
	float median;


#ifdef INPUTCLIP_TEST_MODE
	FILE *histLogFile;
	histLogFile = fopen("inputClip.log","w");
	fprintf(histLogFile,"=====================\n");
#endif

	/*
	 * Here I am checking the status of bit 3 and 4.
	 * Either of these are set if we should use the median as the threashold...
	 * 00001100 = 12
	 */
	if(clipmode & 12){

		medianLimit = nSamplesRead/2;

		/*
		 * We can only do this if we have 4,8 or 16 bit input data!
		 * 1 bit could be done, but it's kinda nonsense...
		 */

		if(nbits > 16){
			fprintf(stderr,"Error: Median computation not supported for nbits > 16\n");
			exit(1);
		}

		/*
		 * To get the median, we need compute a histogram for each channel...
		 * We therefore need a nchansx2^nbits array.
		 * We can quicky do 2^n by left shifting by n.
		 */
		histogramSize = 1 << nbits;
		histogram = (int**)malloc(nchans * sizeof(int*));
		for(i = 0; i < nchans; i++){
			histogram[i] = (int*)malloc(histogramSize * sizeof(int));
		}

		/*
		 * Set the histogram to zero.
		 */
		for(i=0; i<nchans; i++){
			for(j=0; j<histogramSize; j++){ 
				histogram[i][j] = 0;
			}
		}


		/*
		 * Now to compute the histogram.
		 * Basicaly we take the closest integer (to prevent errors from float rounding)
		 * and increment the appropriate box;
		 */
		k = 0;
		for(i=0; i<nSamplesRead; i++){
			for(j=0; j<nchans; j++){
				intValue = (int)(block[k++] + 0.5);
				histogram[j][intValue]++;
			}
		}


#ifdef INPUTCLIP_TEST_MODE
		/** If we are testing, print out the histogram! **/
		for(i=0; i<nchans; i++){
			fprintf(histLogFile,"%d:\t",i);
			for(j = 0; j < histogramSize; j++){
				fprintf(histLogFile,"%d\t",histogram[i][j]);
			}
			fprintf(histLogFile,"\n");
		}
#endif


		/*
		 * Now we need to get the median for each channel!
		 */
		for(i=0; i<nchans; i++){
			counter = 0;
			for(j = 0; j < histogramSize; j++){
				counter += histogram[i][j];
				if(counter > medianLimit){

					median = ((float)j + 0.5 - (float)(counter - medianLimit)/(float)histogram[i][j]);
#ifdef INPUTCLIP_TEST_MODE
					fprintf(histLogFile,"%d: %f\n",i,median);
#endif

					if(clipmode & 4) clipBelow[i] = median;
					if(clipmode & 8) clipAbove[i] = median;
					break;
				}
			}
		}





		/*
		 * Free the histogram array.
		 */

		for(i = 0; i < nchans; i++){
			free(histogram[i]);
		}
		free(histogram);


	}

#ifdef INPUTCLIP_TEST_MODE
	counter = 0;
#endif


	k=0;
	for(i=0; i<nSamplesRead; i++){
		for(j=0; j<nchans; j++){
			if((clipmode & 10) && block[k] > clipAbove[j]){
				 block[k] = clipAbove[j];
#ifdef INPUTCLIP_TEST_MODE
				 counter++;
#endif

			}
			if((clipmode & 5) && block[k] < clipBelow[j]){
				 block[k] = clipBelow[j];
#ifdef INPUTCLIP_TEST_MODE
                                 counter++;
#endif

			}
			k++;
		}
	}


#ifdef INPUTCLIP_TEST_MODE
	fprintf(histLogFile,":: %d\n",counter);
	fclose(histLogFile);
	exit(255);
#endif


}
