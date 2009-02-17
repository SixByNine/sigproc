/*
  bpp2fb - converts BPP search-mode data into "filterbank" data 

  BPP search-mode data is read as unsigned chars each consisting of
  two 4-bit words. Reading the data in this way means that no byte
  swapping is required. Data are read in as blocks of BLOCKSIZE samples
  so that BLOCKSIZE/2 unsigned chars are read in each time. 
*/

#include "filterbank.h"
#include "bpphdr.h" 
BPP_SEARCH_HEADER bpp_search; 
void bpp2fb(FILE *input, FILE *output) /* includefile */
{
  FILE *fpou;
  int np,ns,nc,i,c,b,s,nchars,idump,doit,opened,nsblk,i1,i2,*chtab,blocksize;
  int channel;
  float *tempblock, *datablock,realtime=0.0;
  unsigned char *charblock,sample;
  static unsigned char *charcarry;
  unsigned short *shortblock;
  char string[80];
  int nshift[8] = {28,24,20,16,12,8,4,0};
  static int ncarryover;
  static int first = 1;
  static int fileidx_start=1,file_count=1;
  static double end_time;
  double scantime,sample_start,sample_final;
  int sample_skip,bytestart;
  int bytefinal,fileidx_final,sample_diff,byte_diff=0;
  int bpp_headersize = 32768;
  int rd_jnq=0;
  double sample_end,byte_end;

  np=idump=opened=0;
  ns=512;
  blocksize=ns*nchans;
  nc=nchans;

  if(first) {
    charcarry  = (unsigned char *) malloc(sizeof(unsigned char)*blocksize);
    ncarryover = 0;
    first = 0;
    scantime = (double)((double)bpp_search.file_size)/((double)(nchans/2))*bpp_search.samp_rate*1.e-6;

    if(start_time){
      fileidx_start = ceil(start_time/scantime);
      file_count = fileidx_start;
      sample_skip = floor(start_time/tsamp);
      sample_start = sample_skip - (double)((fileidx_start-1)*(scantime/(double)tsamp));
      bytestart = (sample_start*(double)nchans)/2.;

      if(bytestart<blocksize/2){
	bytestart+=(double)bpp_search.file_size-blocksize/2;
	fileidx_start-=1;
	file_count -=1;
	sample_skip -= ns;
	sample_start = sample_skip - (double)((fileidx_start-1)*(scantime/(double)tsamp));
      }

      realtime += sample_skip*(double)tsamp - scantime*(fileidx_start-1);

      if((rd_jnq =fseek(input,bytestart,SEEK_CUR)< 0)) 	fprintf(stderr,"Error seeking to data byte %d in file\n",bytestart);
    }
    
    if(final_time){
      sample_end = ceil(final_time/(double)tsamp);
      fileidx_final = ceil(final_time/scantime);
      sample_final = sample_end-(double)((fileidx_final-1)*scantime/tsamp);
      byte_end = ceil(final_time/(double)tsamp)*(double)nchans/2;
      bytefinal = (double)sample_final*(double)nchans/2;
      end_time = (double)(byte_end/((double)nchans/2)*(double)tsamp);

      
      fprintf(stderr,"Ending Time:   \n");
      fprintf(stderr,"          End File #                %2d\n",fileidx_final);
      fprintf(stderr,"          End Sample #    %12.3f\n",sample_final);
      fprintf(stderr,"          End Time (s)    %12.5f     (w.r.t. File # 1)\n",end_time);

      if(start_time) {
	end_time -= (double)((fileidx_start-1)*scantime);
	fprintf(stderr,"          End Time (s)    %12.5f     (w.r.t. File # %d)\n",end_time,fileidx_start);
      }

    }
    if (!final_time)    end_time = 1000*scantime;
  }

  tempblock  = (float *) malloc(sizeof(float)*blocksize);
  datablock  = (float *) malloc(sizeof(float)*blocksize);
  charblock  = (unsigned char *) malloc(sizeof(unsigned char)*blocksize);
  shortblock = (unsigned short *) malloc(sizeof(unsigned short)*blocksize);

  if (headerfile) {
    /* write output ASCII header file */
    fpou=open_file("head","w");  
    fprintf(fpou,"Original BPP file: %s\n",inpfile);
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

  chtab=bpp_chans(bpp_search.bandwidth,bpp_search.mb_start_address,
		  bpp_search.mb_end_address,bpp_search.mb_start_board,
		  bpp_search.mb_end_board,bpp_search.cb_id,
		  bpp_search.aib_los,bpp_search.dfb_sram_freqs,
		  bpp_search.rf_lo);


  /************************************************/
  /* main loop over each record of input data file*/
  /************************************************/
  while (!feof(input)&& realtime<=end_time) {
    /* read in a record */

    if(!ncarryover)  /* No data left from previous file */
      nchars=fread(charblock,sizeof(unsigned char),blocksize/2,input);

    if(ncarryover>0) { /* Add to partial block left from previous file */
      nchars=fread(charblock,sizeof(unsigned char),(blocksize/2-ncarryover),input);
      for(c=0;c<nchars;c++) 
	charblock[c+ncarryover] = charblock[c];
      for(c=0;c<ncarryover;c++) 
	charblock[c] = charcarry[c];
      ncarryover = 0;
      nchars = blocksize/2;
      fprintf(stderr,"          Starting at beginning of File # %d\n",file_count);
    }
    if(!ncarryover && nchars<blocksize/2) { /* Don't process, just keep */
      ncarryover = nchars;
      for(c=0;c<nchars;c++) 
	charcarry[c] = charblock[c];
      file_count++;
      if(final_time)      end_time = end_time - (double)tsamp*(int)(scantime/(double)tsamp);
      fprintf(stderr,"Advancing to file # %d\n",file_count);
      /*      if(final_time)   fprintf(stderr,"            End time =  %f     (w.r.t. file # %d)\n",end_time,file_count);*/
      /*      fprintf(stderr,"          Realtime is           %f\n",realtime);*/
    }
    else {
      /* decide whether to process it */
      if ( (doit=process(realtime,sample_skip*(double)tsamp,end_time)) == -1) {
	fprintf(stderr,"realtime at time of break = %f (s)\n ",realtime);
	break; 
      }
      doit = 1;
      if (doit) {
	/* about to process this record, update log */
	np++;
	if (np%10 == 0) {
	  if (!opened) {
	    open_log("filterbank.monitor");
	    opened=1;
	  }
	  sprintf(string,"time:%.1fs",realtime);
	  update_log(string);
	}
	
	/* unscramble the 4-bit data into the float tempblock */
	
	nsblk=0;
	for (c=0;c<nchars;c++) {
	  char2ints(charblock[c],&i1,&i2);
	  tempblock[nsblk++] = (float) i2;
	  tempblock[nsblk++] = (float) i1;
	}
	
	/* update wallclock time */
	realtime+=(float) (nsblk/nchans/nifs) * (float) tsamp; 
	
	if (sumifs) for (s=0;s<nsblk;s++) datablock[s]=0.0;
	s=i1=i2=0;
	/* loop over all samples, summing IFs 1+2 -> total power if neccessary */
	while (s<nsblk) {
	  for (i=0;i<nifs;i++) {
	    channel=0;
	    if (invert_band) channel=nchans-1;
	    for (c=0;c<nchans;c++) {
	      if (invert_band) 
	      if (sumifs) {
		if (i<2) datablock[i1+channel]+=tempblock[i2+chtab[i*nchans+c]];
	      } else {
		datablock[i1+i*nchans+channel]=tempblock[i2+chtab[i*nchans+c]];
	      }
	      if (invert_band)
		channel--;
	      else
		channel++;
	      s++;
	    }
	  }
	  
	  /* update internal counters */
	  i2+=nifs*nchans;
	  if (sumifs) {
	    i1+=nchans;
	  } else {
	    i1+=nifs*nchans;
	  }
	}
	/* divide by 2 to in sumif mode to allow for 4-bit packing */
	if (sumifs) {
	  nsblk/=nifs;
	  for (s=0;s<nsblk;s++) datablock[s]/=2.0;
	}
	/* decide on how to write out data */
	if (obits==32) {
	  /* user has requested floating point numbers in binary file */
	  if (swapout) for (s=0;s<nsblk;s++) swap_float(&datablock[s]);
	  fwrite(datablock,sizeof(float),nsblk,output);
	} else if (obits==16) {
	  /* user has requested unsigned shorts in binary file */
	  float2short(datablock,nsblk,0.0,15.0,shortblock);
	  if (swapout) for (s=0;s<nsblk;s++) swap_short(&shortblock[s]);
	  fwrite(shortblock,sizeof(unsigned short),nsblk,output);
	} else if (obits==8) {
	  /* user has requested unsigned chars in binary file */
	  float2char(datablock,nsblk,0.0,15.0,charblock);
	  fwrite(charblock,sizeof(unsigned char),nsblk,output);
	} else if (obits==4) {
	  /* default is to  write data out packed into character format */
	  float2four(datablock,nsblk,0.0,15.0,charblock);
	  fwrite(charblock,sizeof(unsigned char),nsblk/2,output);
	} 
	else { 
	  error_message("unknown bit format for writing");
	}
      }
      /* update dumps read/processed */
      idump+=ns;
    }
  }
}
