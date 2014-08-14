#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#include "filterbank.h"
void newgmrt2fb(FILE *input, FILE *output) /* includefile*/
{
	unsigned char *buffer = (unsigned char*)malloc(nchans);
	int ich;
	int flipch;
	unsigned char temp;

	while (!feof(input)){
		fread(buffer,1,nchans,input);
		for(ich=0; ich< nchans; ich++){
			buffer[ich] &= gmrtzap[ich];
		}
		if(invert_band){

			for(ich=0; ich < nchans/2; ich++){
				flipch = nchans - ich-1;
				SWAP(buffer[ich],buffer[flipch]);
			}
		}
		fwrite(buffer,1,nchans,output);
	}



}




