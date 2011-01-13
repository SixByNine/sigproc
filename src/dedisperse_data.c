#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "dedisperse.h"
/* 
   orders incoming blocks of data into dedispersed sub-bands 
 */
void dedisperse_data(FILE *input, FILE *output) /*includefile*/
{ 
	char message[80];
	float *buff[2], *dedisp, realtime, nextbaseline, *offset, *tmpblk;
	int readnext=0,isamp,bnum,nsamp,i,j,s,c,b,indx,ns[2],soffset,ddidx;
	int ic,ixnb,ixnc,*ishift,maxshift,nsblk,nsmax,cpb,d,spb,x,nsout,nxb;
	int *ignore;
	unsigned long long totalsamp;

	totalsamp = 0;

	/* calculate table of shift values given the filterbank and DM */
	ishift=dmshift(fch1,foff,nchans,nbands,userdm,refrf,tsamp,frequency_table);
	maxshift=ishift[nchans-1];

	/* set the buffer size based on the maximum shift */
	nsblk=256*nchans; nsout=32*nchans;  
	/*nsblk=256*nchans; nsout=32768*nchans;*/
	nsmax=maxshift*nifs*nchans;
	if (nsmax>nsblk) nsblk=nsmax;
	nxb=nifs*nbands;

	/* define the buffers and initialize counters */
	dedisp =(float *) malloc(nxb*nsout*sizeof(float));
	offset =(float *) malloc(nxb*sizeof(float));
	tmpblk =(float *) malloc(nsout*sizeof(float));
	buff[0]=(float *) malloc(nsblk*sizeof(float));
	buff[1]=(float *) malloc(nsblk*sizeof(float));
	for (i=0;i<nxb;i++) offset[i]=0.0;
	d=bnum=isamp=0;
	ic=nchans*nifs;
	nextbaseline=realtime=0.0;

	/* zero any channels that are in the ignored list of channels */
	if (file_exists(ignfile)) {
		ignore=ignored_channels(ignfile,nchans);
	} else {
		ignore=(int *) malloc(nchans*sizeof(int));
		for (i=0;i<nchans;i++) ignore[i]=0;
	}

	/* number of channels per band to dedisperse (cpb) must be an integer */
	cpb=nchans/nbands; 
	if ((cpb*nbands) != nchans) error_message("silly sub-band selection!");

	/* main loop - keep going until no more data comes in */
	while (1) {

		/* read in the buffer to be processed if not done so already */
		if (!readnext) {
			sprintf(message,"time:%.1fs:DM:%.1fpc/cc",realtime,refdm);
			update_log(message);
			if ((ns[bnum] = read_block (input, nbits, buff[bnum], nsblk)) <= 0){
				/* There were no more samples to read */
				if (isamp)write_dedisp(dedisp,isamp,nifs,nbands,offset,output);
				return;
			} else {
				/* We managed to read some samples!
				 * There may be some things we want to do to the data as we read it in
				 * These are done here if they need to be done.
				 */

				if(ns[bnum]%nchans){
					/* If we get here then ns[bnum]/nchans is not integer!
					 * This is really bad, as we assume that we read in an
					 * integer number of samples.
					 */
					fprintf(stderr,"ERROR: Could not read an integer number of samples (nread/nchans is not integer).\nPerhaps the file is incomplete or the stream was terminated.\nUnfortunately dedispersion cannot continue.\n");
					exit(2);
				}

 				/* Do we want to clip the input channel data
				 * Note that we must do this BEFORE zerodming
                                 */
                                if(doInputClip){
                                        inputClip(buff[bnum],ns[bnum]/nchans,inputClipLow,inputClipHigh,doInputClip); /* It is important that ns[bnum]/nchans is an integer! */
                                }


				/* Do we want to do the 'zero dm' subtraction?
				 */
				if(doSubZero){
					subzero(buff[bnum],ns[bnum]/nchans); /* It is important that ns[bnum]/nchans is an integer! */
				}

				if(doBaselineFile){
					reverseBaselineNormalise(buff[bnum],totalsamp,ns[bnum]/nchans,nchans);
				}

			}

		}

		/* number of samples in this buffer */
		nsamp=ns[bnum]/ic;

		/* flag to signify whether next buffer has been read in (1=yes;0=no) */
		readnext=0;

		/* dedispersing loop over all samples in this buffer */
		for (s=0; s<nsamp; s++) {
			soffset=isamp*nxb;
			/* loop over the IFs */
			for (i=0; i<nifs; i++) {
				/* number of channels to skip within this IF */
				ixnc=i*nchans; 
				ixnb=i*nbands;
				for (b=0; b<nbands; b++) {
					/* calculate index of this sample */
					ddidx=soffset+ixnb+b;
					/* clear array element for storing dedispersed subband */
					dedisp[ddidx]=0.0; 
					/* loop over the channels in this subband */
					for (c=b*cpb;c<(b+1)*cpb;c++) {
						/* proceed only if selected channel # is not in ignore list */
						if (!ignore[c]) 
						{
							/* calculate index of sample to be added */
							indx=(s+ishift[c])*ic+ixnc+c;
							/* required sample will be in either this buffer or the next */
							if (indx<ns[bnum]) {
								dedisp[ddidx]+=buff[bnum][indx];
							} else {
								if (!readnext) {
									if ((ns[!bnum]=read_block(input,nbits,buff[!bnum],nsblk))<=0) {
									/* There are no more samples to read, nread <= 0 */
										if (isamp) {
											write_dedisp(dedisp,isamp,nifs,nbands,offset,output);
										}
										return;
									} else {
										/* We read some samples, as nread > 0.
										 * There may be some things we want to do to the data as we read it in
										 * These are done here if they need to be done.
										 */

										if(ns[!bnum]%nchans){
											/* If we get here then ns[bnum]/nchans is not integer!
											 * This is really bad, as we assume that we read in an
											 * integer number of samples.
											 */
											fprintf(stderr,"ERROR: Could not read an integer number of samples (nread/nchans is not integer).\nPerhaps the file is incomplete or the stream was terminated.\nUnfortunately dedispersion cannot continue.\n");
											exit(2);
										}

                                                                                /* Do we want to clip the input channel data
                                                                                 */
                                                                                if(doInputClip){
                                                                                        inputClip(buff[!bnum],ns[!bnum]/nchans,inputClipLow,inputClipHigh,doInputClip);
											/* It is important that ns[!bnum]/nchans is an integer! */
                                                                                }

										/* Do we want to do the 'zero dm' subtraction?
										 */
										if(doSubZero){
											subzero(buff[!bnum],ns[!bnum]/nchans); /* It is important that ns[!bnum]/nchans is an integer! */
										}
										if(doBaselineFile){
											reverseBaselineNormalise(buff[!bnum],totalsamp,ns[!bnum]/nchans,nchans);
										}



									}
									sprintf(message,"time:%.1fs:DM:%.1fpc/cc",realtime,refdm);
									update_log(message);
									readnext=1;
								}
								dedisp[ddidx]+=buff[!bnum][indx-ns[bnum]];
							}
						}
					}
				} /* end of loop over subbands */
			} /* end of loop over IFs */
			/* update number of samples dedispersed and elapsed time */
			isamp++; realtime+=tsamp;
			totalsamp++;
			if (isamp==nsout) {
				if (baseline) {
					for (i=0;i<nifs;i++) {
						ixnb=i*nbands;
						for (b=0;b<nbands;b++) {
							for (j=0;j<nsout;j++) tmpblk[j]=dedisp[j*nxb+ixnb+b];
							offset[ixnb+b]=nrselect(nsout/2,nsout,tmpblk-1);
							//fprintf(stderr,"%d %f\n",ixnb+b,offset[ixnb+b]);
						}
					}
				}
				write_dedisp(dedisp,nsout,nifs,nbands,offset,output);
				isamp=0;
			}
		} /* end of loop over samples */
		/* switch to next buffer */
		bnum=!bnum;
	} /* end of main loop */
}

/* subtract current offset from the dedisperse time samples and write */
void write_dedisp(float *dedisp, int nsout, int nifs, int nbands, float *offset, FILE *output)/*includefile*/
{
	int s,i,b,ixnb,sxib,n;
	static int first=1;
	static float *clipthreshold;
	float *temp,outliers,sample,sumsq;
	unsigned char *onebyte;
	short *twobyte;

	/* multiply outgoing data by Jansky calibration factor if supplied */
	if (jyfactor != 1.0) for (i=0;i<nsout*nifs*nbands;i++) dedisp[i]*=jyfactor;

	if (first) {
		/* allocate an array for holding blocks from a given subband */
		temp=malloc(sizeof(float)*nsout);
		/* allocate an array for saving the clipping threshold */
		clipthreshold=malloc(sizeof(float)*nbands*nifs);
		for (i=0;i<nifs*nbands;i++) clipthreshold[i]=0.0;
	}

	for (i=0;i<nifs;i++) {
		ixnb=i*nbands;
		for (b=0;b<nbands;b++) {
			if (first) {
				/* copy sub-band into temporary store for absolute value */
				for (s=0;s<nsout;s++) {
					sxib=s*nifs*nbands;
					temp[s]=fabs(dedisp[sxib+ixnb+b]-offset[ixnb+b]);
				}
				/* find the value below which 90% of the samples lie */
				outliers=nrselect(nsout/10,nsout,temp-1);
				n=0;
				sumsq=0.0;
				/* calculate sum of squares based on the inner 90% of samples */
				for (s=0;s<nsout;s++) {
					if (temp[s]<outliers) {
						sumsq+=temp[s]*temp[s];
						n++;
					}
				}
				/* now set the threshold based on the sum of squares */
				if (n) 
					clipthreshold[ixnb+b]=clipvalue*sqrt((double)sumsq/(double)n);
				else 
					clipping=0;
			}
			for (s=0;s<nsout;s++) {
				sxib=s*nifs*nbands;
				/* subtract off median value of this block */
				sample=dedisp[sxib+ixnb+b]-offset[ixnb+b];
				/* clip this sample if it exceeds the threshold */
				if (fabs(sample)>clipthreshold[ixnb+b] && clipping) sample=0.0;
				/* store the final produce and swap bytes if necessary */
				dedisp[sxib+ixnb+b]=sample;
				if (swapout) swap_float(&dedisp[sxib+ixnb+b]);
			}
		}
	}

	/* now write out samples and bat on */
	switch (nobits) {
		case 8:
#if SIGNED
			onebyte = (char *) malloc(nsout*nifs*nbands);
			for (i=0; i<nsout*nifs*nbands; i++) 
				onebyte[i] = ((char) dedisp[i])-128;
			fwrite(onebyte,sizeof(char),nsout*nifs*nbands,output);
#else
			onebyte = (unsigned char *) malloc(nsout*nifs*nbands);
			for (i=0; i<nsout*nifs*nbands; i++) 
				onebyte[i] = ((unsigned char) dedisp[i]);
			fwrite(onebyte,sizeof(unsigned char),nsout*nifs*nbands,output);
#endif
			break;
		case 16:
			twobyte = (short *) malloc(nsout*nifs*nbands);
			for (i=0; i<nsout*nifs*nbands; i++) 
				twobyte[i] = (short) dedisp[i];
			fwrite(twobyte,sizeof(short),nsout*nifs*nbands,output);
			break;
		case 32:
			fwrite(dedisp,sizeof(float),nsout*nifs*nbands,output);
			break;
		default:
			error_message("requested output number of bits can only be 8, 16 or 32");
			break;
	}

	if (first) {
		first=0;
		free(temp);
	}
}
