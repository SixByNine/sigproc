/* pspm2fb - converts PSPM search-mode data into "filterbank" data */

#include "filterbank.h"
#include "pspmhdr.h"
PSPM_SEARCH_HEADER pspm_search;

void pspm2fb(FILE *input, FILE *output) /* includefile */
{
  FILE *fpou;
  int np=0,ns,nc,c,s,nints,rawdata[PSPM_INT_BLOCK],idump,doit,swap_bytes;
  int i,opened=0,drift;
  float datablock[PSPM_REA_BLOCK],datablock2[PSPM_REA_BLOCK],sum,realtime;
  unsigned char charblock[PSPM_REA_BLOCK],sample;
  unsigned short shortblock[PSPM_REA_BLOCK];
  char string[80];

  /* establish whether this is drift-mode data */
  if (pspm_search.HEADER_TYPE == 0) {
	drift=1;
  } else {
	drift=0;
  }
  idump=0;

  /* establish whether we need to swap bytes (PSPM is big endian) */
  swap_bytes=little_endian();

  /* shorthand for number of samples and channels in a datablock */
  ns=PSPM_SAM_BLOCK;
  nc=PSPM_NCH_BLOCK;

  if (headerfile) {
    /* write output ASCII header file */
    fpou=open_file("head","w");  
    fprintf(fpou,"Original PSPM file: %s\n",inpfile);
    fprintf(fpou,"Sample time (us): %f\n",tsamp*1.0e6);
    fprintf(fpou,"Number of samples/record: %d\n",ns);
    fprintf(fpou,"Center freq (MHz): %f\n",fch1+(float)nchans*foff/2.0);
    fprintf(fpou,"Channel band (kHz): %f\n",fabs(foff)*1000.0);
    fprintf(fpou,"Number of channels/record: %d\n",nc);
    fprintf(fpou,"Time stamp (MJD): %18.12f\n",tstart);
    fprintf(fpou,"AZ at start: %f\n",az_start);
    fprintf(fpou,"ZA at start: %f\n",za_start);
    fprintf(fpou,"RA (J2000): %f\n",src_raj);
    fprintf(fpou,"DEC (J2000):  %f\n",src_dej);
    fclose(fpou);
  }


  /* main loop over each record of input data file*/
  while (!feof(input)) {
    /* this is the time at that start of the block */
    realtime=tsamp*idump;
    /* read in a record and unscramble the channels */
    nints=fread(rawdata,sizeof(int),PSPM_INT_BLOCK,input);
    if ( (doit=process(realtime,start_time,final_time)) == -1) break;
    if (doit) {
      /* about to process this record, update log */
      np++;
      if (np%10 == 0) {
	if (!opened) {
	  /* open up logfile */
	  open_log("filterbank.monitor");
	  opened=1;
	}
	sprintf(string,"time:%.1fs",realtime);
	update_log(string);
      }

      if (swap_bytes) for (s=0;s<nints;s++) swap_int(&rawdata[s]);

      /* unscramble the channels using Ingrid's C routine */
      pspm_decode(rawdata,datablock);

      /* if the -invert option was specified, flip the band */
      i=0;
      if (invert_band) {
	for(s=0;s<ns;s++) {
	  for (c=nc-1;c>=0;c--) {
	    datablock2[i++]=datablock[s*nc+c];
	  }
	}
        for(i=0;i<ns*nc;i++) 
	    datablock[i]=datablock2[i];
      }

      realtime+=(float) ns * (float) tsamp; 
      /* decide on how to write out data */
      if (obits==32) {
	/* user has requested floating point numbers in binary file */
	if (swapout) for (s=0;s<ns*nc;s++) swap_float(&datablock[s]);
	fwrite(datablock,sizeof(float),ns*nc,output);
      } else if (obits==16) {
	/* user has requested unsigned shorts in binary file */
	float2short(datablock,ns*nc,0.0,15.0,shortblock);
	if (swapout) for (s=0;s<ns*nc;s++) swap_short(&shortblock[s]);
	fwrite(shortblock,sizeof(unsigned short),ns*nc,output);
      } else if (obits==8) {
	/* user has requested unsigned chars in binary file */
	float2char(datablock,ns*nc,0.0,15.0,charblock);
	fwrite(charblock,sizeof(unsigned char),ns*nc,output);
      } else if (obits==4) {
	/* default is to  write data out packed into character format */
	float2four(datablock,ns*nc,0.0,15.0,charblock);
	fwrite(charblock,sizeof(unsigned char),ns*nc/2,output);
      } else if (obits==0) {
	/* special mode to write out in different order for old ddsp program */
	for (c=0;c<nc;c++) {
	  for(s=0;s<ns;s++) {
	    sample=(unsigned char)datablock[s*nc+c];
	    fwrite(&sample,sizeof(unsigned char),1,output);
	  }
	}
      } else { 
	error_message("unknown bit format for writing");
      }
    }
    /* update dumps read/processed */
    idump+=ns;
    /* break out if this is in drift-mode and we've read 9 beams */
    if (drift && (idump == 4718592)) break;  
  }
}
