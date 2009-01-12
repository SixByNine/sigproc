/* 
   look at the raw data file and establish origin and whether it's readable 
   the following data types are currently recognizable:
   
   0 - sigproc binary profiles
   1 - pspm search
   2 - pspm timing
   3 - wapp search
   4 - wapp timing
   5 - bpp search
   6 - bpp timing
   7 - aoftm search (no timing mode for this machine)
   8 - ooty search  (no timing mode for this machine)
   9 - parkes 1-bit
  10 - GMRT data
  11 - PULSAR2000 search data (Effelsberg); new data format (27.08.2002 BK)

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sigproc.h"
#include "header.h"

#include "pspmhdr.h"
PSPM_TIMING_HEADER pspm_timing;
PSPM_SEARCH_HEADER pspm_search;

#include "wapp_header.h"
#include "key.h"
struct WAPP_HEADER *wapp;
struct WAPP_HEADER head;

#include "bpphdr.h"
BPP_TIMING_HEADER bpp_timing;
BPP_SEARCH_HEADER bpp_search;  

#include "scamp_header.h"
int scamp_block_size, scamp_rawdata;
SCAMP_HEADER schead;
char *get_scamp_header(char *header_parameter, int nbytes) {
  char *string;
  string=(char *) malloc(nbytes);
  strcpy(string,"");
  strncat(string,header_parameter,nbytes);
  return(string);
}

int invert_band, wapp_lagtrunc, psr2Ksubbaseline;
int wapp_file, wapp_lagformat, wapp_sum, wapp_level, wapp_flip, obits, *chtab;
int wapp_header_size,wapp_incfile_length,wapp_isalfa,wapp_number;
double wapp_obstime, alfa_ang, alfa_raj[7], alfa_dej[7];

int scamp_ignore[1024],scamp_chans;

int typeof_inputdata(FILE *fptr, char *filename) /* includefile */
{ 
  FILE *hdr, *ignore_file;
  char *string, *label, junk[14], cbuffer[100], cdum[100];
  int dayno,year,scan,hh,mm,ss,mo,dd,swap_bytes,rah,ram,ded,dem,i,j;
  double ras,des,a,b;
  struct HEADERKEY *key;
  struct HEADERP *h;
  int fsys;
  char *dummy;
  unsigned short gmrt[3];

  /* set machine and telescope ID to unknown type initially */
  machine_id=telescope_id=0;
  rewind(fptr);

  nbeams = 1;
  ibeam = 1;
  if (read_header(fptr)) {
    if (data_type == 3) return(0);
  }

  /* try reading as PULSAR2000 data		-> case 11 */
  rewind( fptr );
  fgets( cbuffer, 100, fptr );

  if ( strncmp(cbuffer, "@@@PULSAR2000@@@", 16) == 0 )  {

    machine_id=8;				// PULSAR2000
    telescope_id=8;				// Effelsberg
    data_type=1;				// "filterbank"

    fgets( cbuffer, 100, fptr );			// ---------------------------
    fgets( cbuffer, 100, fptr );			// scan_number (notused!)
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
    sscanf(cbuffer, "%s %s", cdum, source_name);	// source_name
    fgets( cbuffer, 100, fptr );			// date (not used!)
    fgets( cbuffer, 100, fptr );			// ---------------------------       
    
    // center frequency (= sky frequency) [MHz]
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
    sscanf(cbuffer, "%s %lf %s", cdum, &fch1, cdum);

    // number of filterbank channels (filterbanks: 8, 16, 32; PSE++: 30, 60)
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);    
    sscanf(cbuffer, "%s %i", cdum, &nchans);
    
    // filterbank channel bandwidth (->convert to negative, defination)
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);      
    sscanf(cbuffer, "%s %s %lf %s", cdum, cdum, &foff, cdum);
    if ( (foff > 0.6) && (foff < 0.7) )  foff = -1.0*(2.0/3.0);	// PSE++: -0.66666 MHz
    if ( (foff > 1.3) && (foff < 1.4) )	 foff = -1.0*(4.0/3.0);	// PSE++: -1.33333 MHz
    if (foff > 0.0)  foff *= (-1.0);	// use NB or BB filterbanks: 4 MHz or 60 MHz
    
    // we assume two polarisation channels for both PSE and FB data
    nifs = 2;		// !!! nifs = 2 !!!

    // set number of channels per filterbank
    nchans /= nifs;

    // calculate the center frequency of the first filterbank channel
    fch1 -= ((nchans/2)-0.5) * foff;

    // sample time [s]
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
    sscanf(cbuffer, "%s %lf %s", cdum, &tsamp, cdum);
    tsamp /= (1000.0*1000.0);

    // number of bits per time sample: 4, 8 or 16 bits
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
    sscanf(cbuffer, "%s %i %s", cdum, &nbits, cdum);    

    // time series with subtracted baseline?
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
    sscanf(cbuffer, "%s %s", cdum, cbuffer);
    if (strncmp(cbuffer,"substract",9) == 0) psr2Ksubbaseline=1;
    else				     psr2Ksubbaseline=0;

    fgets( cbuffer, 100, fptr ); 		// timing-reference (not important!)
    fgets( cbuffer, 100, fptr );		// ---------------------------    

    // check if we have the new or old data format!
    fgets( cbuffer, 100, fptr );
    if ( strncmp(cbuffer, "RA1950", 6) == 0 )  {
      // --> new data format
      fgets( cbuffer, 100, fptr ); 		// DEC1950 (not used!)
    
      // RA2000
      fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
      sscanf(cbuffer, "%s %i %i %i", cdum, &hh, &mm, &ss);
      src_raj  = (double)hh*10000.0 + (double)mm*100.0 + (double)ss;

      // DEC2000
      fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
      sscanf(cbuffer, "%s %i %i %i", cdum, &hh, &mm, &ss);
      src_dej  = (double)hh*10000.0 + (double)mm*100.0 + (double)ss;

      fgets( cbuffer, 100, fptr ); 		// gal. Long.  (not used!)
      fgets( cbuffer, 100, fptr ); 		// gal. Lat.   (not used!)

      // AZM(start) [degree]
      fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
      sscanf(cbuffer, "%s %lf", cdum, &az_start);

      // ELV(start) [degree]
      fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
      sscanf(cbuffer, "%s %lf", cdum, &za_start);

      fgets( cbuffer, 100, fptr );		// ---------------------------
      fgets( cbuffer, 100, fptr ); 		// start UTC (not used!)
      fgets( cbuffer, 100, fptr ); 		// start LST (not used!)   
    }
    
    // start MJDN
    fgets( cbuffer, 100, fptr );  eraseDP( cbuffer);
    sscanf(cbuffer, "%s %s %lf", cdum, cdum, &tstart);

    return 11;
  }		// end of: PULSAR2000 data (new format)
    
  /* 
     Try reading as Parkes fast-sampled data       -> case 9 
     N.B. key is to check first integer read in - this should be 48/24k.
     48k is the usual case where one filter system has been used
     24k is the rarer case where two systems were used but only one
         was read back from tape (functionality to read this added 4/7/04)
	 **** note that currently it is not possible to deal with the
              case where both systems have been read into the same file ****
     In any case, depending on whether the data were written by big/little 
     endian machines, the integer will need to be swapped - so try both 
     possibilities. Since the rest of the data are chars, no further byte 
     swapping required. This trick allows one to read the data regardless 
     of the origin of the tape reader.
  */
  rewind(fptr);
  scamp_block_size=0;
  fread(&junk,10,1,fptr); /* dummy pad to ensure correct read */
  fread(&schead,sizeof(schead),1,fptr);
  string=get_scamp_header((schead.telid),sizeof(schead.telid));
  if (strings_equal(string,"PARKES    ")) {
    scamp_rawdata=1;
    scamp_block_size=49152;
  } else {
    scamp_rawdata=0;
    rewind(fptr);
    fread(&i,4,1,fptr);
    j=i;
    swap_int(&j);
    if (i == 49152 || j == 49152) scamp_block_size=49152;
    if (i == 24576 || j == 24576) scamp_block_size=24576;
  }
  if (scamp_block_size) {
    if (!scamp_rawdata) {
      string=(char *) malloc(strlen(filename)+1);
      strcpy(string,"");
      strncat(string,filename,strlen(filename)-3);
      strcat(string,"hdr");
      hdr=open_file(string,"rb");
      fread(&junk,14,1,hdr); /* dummy pad to ensure correct read */
      fread(&schead,sizeof(schead),1,hdr);
      fclose(hdr);
    }
    string=get_scamp_header(schead.ra,sizeof(schead.ra));
    rah=atoi(strtok(string,":"));
    ram=atoi(strtok(NULL,":"));
    ras=atof(strtok(NULL,":"));
    src_raj=ras+ram*100.0+rah*10000.0;
    string=get_scamp_header(schead.dec,sizeof(schead.dec));
    ded=atoi(strtok(string,":"));
    dem=atoi(strtok(NULL,":"));
    des=atof(strtok(NULL,":"));
    src_dej=des+100.0*dem+10000.0*abs(ded);
    if (ded<0) src_dej*=-1.0;
    string=get_scamp_header(schead.mjd,sizeof(schead.mjd));
    tstart=atof(string);
    string=get_scamp_header(schead.ut,sizeof(schead.ut));
    tstart+=atof(strtok(string,":"))/24.0+atof(strtok(NULL,":"))/1440.0+
      atof(strtok(NULL,":"))/86400.0;
    string=get_scamp_header(schead.nsys,sizeof(schead.nsys));
    fsys=atoi(string)-1; /* get the index of the filter system to use */
    if (fsys<0) fsys=0;
    if (scamp_rawdata) fsys=0;
    string=get_scamp_header(schead.fch1[fsys],sizeof(schead.fch1[fsys]));
    fch1=atof(string);
    string=get_scamp_header(schead.chbw[fsys],sizeof(schead.chbw[fsys]));
    foff=atof(string);
    string=get_scamp_header(schead.nchan[fsys],sizeof(schead.nchan[fsys]));
    nchans=atoi(string);
    scamp_chans=nchans;
    /* read in list of segments to ignore */
    for (i=0;i<1024;i++) scamp_ignore[i]=0;
    if (file_exists("scamp.ignore")) {
      ignore_file=open_file("scamp.ignore","r");
      while (1) {
	fscanf(ignore_file,"%d",&i);
	if (feof(ignore_file)) break;
	scamp_ignore[i]=1;
	nchans-=8;
      }
      fclose(ignore_file);
    }
    raw_fch1 = fch1;
    raw_foff = foff;
    if (foff>0) {
      invert_band=1;
      fch1+=foff*(double)(nchans-1);
      foff*=-1.0;
    }
    

    string=get_scamp_header(schead.l,sizeof(schead.l));
    gal_l=atof(string);
    string=get_scamp_header(schead.b,sizeof(schead.b));
    gal_b=atof(string);
    
    string=get_scamp_header(schead.nbeam,sizeof(schead.nbeam));
    while(*string==' ')string++;
    nbeams=atoi(string);
    string=get_scamp_header(schead.ibeam,sizeof(schead.ibeam));
    while(*string==' ')string++;
    ibeam=atoi(string);
    string=get_scamp_header(schead.obstime,sizeof(schead.obstime));
    header_tobs = atof(string); 
    string=get_scamp_header(schead.tsamp[fsys],sizeof(schead.tsamp[fsys]));
    tsamp=atof(string)*1.0e-3;
    string=get_scamp_header(schead.telid,sizeof(schead.telid));
    telescope_id=5;
    if (!strncmp(string,"PARKES",6)) telescope_id=4;
    if (!strncmp(string,"JODRELL",7)) telescope_id=5;
    string=get_scamp_header(schead.psrname,sizeof(schead.psrname));
    strcpy(source_name,strtok(string," "));
    free(string);
    machine_id=6;
    data_type=1;
    nbits=1;
    nifs=1;
    return 9;
  }

  /* Try reading as PSPM/BPP timing-mode data */
  rewind(fptr);
  fread(&bpp_timing,BPP_HEADER_SIZE,1,fptr); 

  /* establish whether we need to swap bytes (PSPM/BPPs are big endian) */
  swap_bytes=little_endian();
  if (swap_bytes) {
    swap_long(&bpp_timing.HEADER_TYPE);
    swap_long(&bpp_timing.BACKEND_TYPE);
  }
  if (bpp_timing.HEADER_TYPE == 1 || bpp_timing.HEADER_TYPE == 3) {

    switch (bpp_timing.BACKEND_TYPE) {
    case 0:
      /* this is PSPM data */
      machine_id=1;
      telescope_id=1;
      break;
    case 2:
      /* this is NBPP data */
      machine_id=4;
      telescope_id=3;
      break;
    default:
      error_message("unknown backend in BPP/PSPM timing file?!?");
      break;
    }

    switch (machine_id) {
    case 1:
      /* this is PSPM data */
      rewind(fptr);
      fread(&pspm_timing,PSPM_HEADER_SIZE,1,fptr); 
      if (swap_bytes) {
	swap_long(&pspm_timing.scan_num);
	swap_long(&pspm_timing.tc);
	swap_double(&pspm_timing.tick_offset);
	swap_double(&pspm_timing.freq);
	swap_long(&pspm_timing.num_phase_bins);
	swap_long(&pspm_timing.num_chans);
	swap_long(&pspm_timing.num_periods);
	swap_double(&pspm_timing.psr_period);
	swap_double(&pspm_timing.psr_dm);
      }
      /* get number of seconds since midnight (tstart) on MJD mjd_obs */
      tstart=pspm_tstart(pspm_timing.scan_num, pspm_timing.start_time,
		       pspm_timing.tick_offset,&mjdobs);
      read_aoscan(pspm_timing.scan_num,&dayno,&year,&scan_number);
      return 2;
      break;
    case 4:
      /* this is BPP data */
      if (swap_bytes) {
	swap_long(&bpp_timing.scan_num);
	swap_double(&bpp_timing.tick_offset);
	swap_long(&bpp_timing.num_phase_bins);
	swap_long(&bpp_timing.num_chans);
	swap_long(&bpp_timing.num_periods);
	swap_double(&bpp_timing.psr_period);
	swap_double(&bpp_timing.psr_dm);
	swap_double(&bpp_timing.bandwidth);
	swap_double(&bpp_timing.rf_lo);
	swap_int(&bpp_timing.mb_start_address);
	swap_int(&bpp_timing.mb_end_address);
	swap_int(&bpp_timing.mb_start_board);
	swap_int(&bpp_timing.mb_end_board);
	for (i=0;i<MAXNUMCB;i++) 
	  swap_int(&bpp_timing.cb_id[i]);
	for (i=0;i<MAX_NUM_LO_BOARDS;i++) 
	  swap_double(&bpp_timing.aib_los[i]);
	for (i=0;i<FB_CHAN_PER_BRD;i++)
	  swap_float(&bpp_timing.dfb_sram_freqs[i]);
      }
      /* parse search date dayno:year */
      string=(char *)malloc(strlen(bpp_timing.date));
      strcpy(string,bpp_timing.date);
      dayno=atoi(strtok(string,":"));
      dayno=atoi(strtok(NULL,":"));
      year=atoi(strtok(NULL,":"))+1900;
      free(string);
      /* parse UT start time hh:mm:ss */
      string=(char *)malloc(strlen(bpp_timing.start_time));
      strcpy(string,bpp_timing.start_time);
      hh=atoi(strtok(string,":"));
      mm=atoi(strtok(NULL,":"));
      ss=atoi(strtok(NULL,":"));
      free(string);
      /* calculate MJD start time based on the above info */
      tstart=mjd(year,1,1)+(double)dayno-1.0;
      tstart+=(double)hh/24.0+(double)mm/1440.0+(double)ss/86400.0;
      mjdobs=floor(tstart);
      tstart=(tstart-mjdobs)*86400.0;
      /* 
	 calculate nchans from other header info as 
	 num_chans from NBPP seems to be incorrect! 
      */
      bpp_timing.num_chans = 
	(bpp_timing.mb_end_address/2-bpp_timing.mb_start_address/2+1)*
	(bpp_timing.mb_end_board-bpp_timing.mb_start_board+1)*4;
      return 6;
      break;
    default:
      error_message("coding in typeof_intputdata not yet done!");
      break;
    }

  }

  /* Try reading as PSPM search-mode data */
  rewind(fptr);
  fread(&pspm_search,PSPM_HEADER_SIZE,1,fptr);

  /* establish whether we need to swap bytes (PSPM is big endian) */
  swap_bytes=little_endian();
  if (swap_bytes) {
    swap_long(&pspm_search.HEADER_TYPE);
    swap_long(&pspm_search.BACKEND_TYPE);
  }

  if (pspm_search.BACKEND_TYPE == 0 && ((pspm_search.HEADER_TYPE==0) || (pspm_search.HEADER_TYPE==2))) {
    machine_id=1;

    /* perform byte swapping on header parameters of interest if necessary */
    if (swap_bytes) {
      swap_double(&pspm_search.samp_rate);
      swap_double(&pspm_search.chan_spacing);
      swap_double(&pspm_search.rf_freq);
      swap_long(&pspm_search.num_chans);
      swap_long(&pspm_search.bit_mode);
      swap_long(&pspm_search.scan_num);
      swap_long(&pspm_search.scan_file_number);
      swap_double(&pspm_search.tick_offset);
      swap_double(&pspm_search.user_az);
      swap_double(&pspm_search.user_za);
      swap_double(&pspm_search.user_ra);
      swap_double(&pspm_search.user_dec);
    }

    tsamp=pspm_search.samp_rate*1.0e-6;
    nchans=pspm_search.num_chans;
    nifs=1;
    nbits=pspm_search.bit_mode;

    pspm_search.chan_spacing=0.062; 
    fch1=pspm_search.rf_freq+nchans*pspm_search.chan_spacing/2;
    foff=-1.0*pspm_search.chan_spacing;

    /* get number of seconds since midnight (tstart) on MJD mjd_obs */
    tstart=pspm_tstart(pspm_search.scan_num, pspm_search.start_time,
		       pspm_search.tick_offset,&mjdobs);

    /* 
       if this is drift scan data - correct tstart to start of this file
       each file is 10 lots of pow(2,19) 80 us samples, but only the first
       9 are independent, hence offset is pow(2,19)*9*80e-6 = 377.48736 s 
       N.B. pspm_search.scan_file_number gives number starting at 1
    */
    if (pspm_search.scan_file_number > 1)  
      tstart+=(pspm_search.scan_file_number-1)*377.48736;

    tstart=mjdobs+tstart/86400.0;

    /* get telescope az and za */
    az_start=pspm_search.user_az;
    za_start=pspm_search.user_za;

    /* get source ra and dec */
    src_raj=pspm_search.user_ra;
    src_dej=pspm_search.user_dec;

    strcpy(source_name,pspm_search.psr_name);
    read_aoscan(pspm_search.scan_num,&dayno,&year,&scan_number);

    telescope_id=1;
    return 1;
  } 

  /* Try reading as BPP fast-sampled data        -> case 6 */
  rewind(fptr);
  fread(&bpp_search,BPP_HEADER_SIZE,1,fptr);
  /* establish whether we need to swap bytes (BPPs are big endian) */
  swap_bytes=little_endian();
  if (swap_bytes) swap_long(&bpp_search.BACKEND_TYPE);
if (bpp_search.BACKEND_TYPE == 0 || bpp_search.BACKEND_TYPE == 1 || bpp_search.BACKEND_TYPE == 4) {
    machine_id=4;
    telescope_id=6;
    nifs=1;
    if (swap_bytes) {
      swap_double(&bpp_search.samp_rate);
      swap_double(&bpp_search.bandwidth);
      swap_double(&bpp_search.rf_lo);
      swap_long(&bpp_search.bit_mode);
      swap_long(&bpp_search.num_chans);
      swap_long(&bpp_search.scan_file_number);
      swap_long(&bpp_search.file_size);
      swap_int(&bpp_search.mb_start_address);
      swap_int(&bpp_search.mb_end_address);
      swap_int(&bpp_search.mb_start_board);
      swap_int(&bpp_search.mb_end_board);
      for (i=0;i<MAXNUMCB;i++) 
	swap_int(&bpp_search.cb_id[i]);
      for (i=0;i<MAX_NUM_LO_BOARDS;i++) 
	swap_double(&bpp_search.aib_los[i]);
      for (i=0;i<FB_CHAN_PER_BRD;i++)
	swap_float(&bpp_search.dfb_sram_freqs[i]);
      swap_double(&bpp_search.ra_1950);
      swap_double(&bpp_search.dec_1950);
    }
    /* parse search date dayno:year */
    string=(char *)malloc(strlen(bpp_search.date));
    strcpy(string,bpp_search.date);
    dayno=atoi(strtok(string,":"));
    year=atoi(strtok(NULL,":"));
    free(string);
    /* parse UT start time hh:mm:ss */
    string=(char *)malloc(strlen(bpp_search.start_time));
    strcpy(string,bpp_search.start_time);
    hh=atoi(strtok(string,":"));
    mm=atoi(strtok(NULL,":"));
    ss=atoi(strtok(NULL,":"));
    free(string);
    /* calculate MJD start time based on the above info */
    tstart=mjd(year,1,1)+(double)dayno; /* day0=jan1 for BCPM!!! */
    tstart+=(double)hh/24.0+(double)mm/1440.0+(double)ss/86400.0;
    tstart+=1.0/86400.0; /* now add on the mysterious 1-s offset */
    tsamp=bpp_search.samp_rate*1.0e-6;
    nchans=bpp_search.num_chans;
    nbits=bpp_search.bit_mode;
    /* get source name and ra and dec */
    src_raj=bpp_search.ra_1950;
    src_dej=bpp_search.dec_1950;
    strcpy(source_name,bpp_search.target_name);

    src_raj=src_dej=az_start=za_start=0.0;
    chtab=bpp_chans(bpp_search.bandwidth,bpp_search.mb_start_address,
		    bpp_search.mb_end_address,bpp_search.mb_start_board,
		    bpp_search.mb_end_board,bpp_search.cb_id,
		    bpp_search.aib_los,bpp_search.dfb_sram_freqs,
		    bpp_search.rf_lo);
    if (foff>0.0) foff*=-1.0;
    return 5;
  }

  /* Try reading as WAPP data */
  rewind(fptr);
  wapp = &head;
  if ( (h = head_parse( filename )) != NULL) {
    fetch_hdrval(h,"wapp_time",&(wapp->wapp_time),sizeof(wapp->wapp_time));

    /* establish whether we need to swap bytes (WAPP is little endian) */
    swap_bytes=big_endian();

    fetch_hdrval(h,"src_name",&(wapp->src_name),sizeof(wapp->src_name));
    strcpy(source_name,wapp->src_name);

    fetch_hdrval(h,"obs_type",&(wapp->obs_type),sizeof(wapp->obs_type));

    fetch_hdrval(h,"obs_time",&(wapp->obs_time),sizeof(wapp->obs_time));
    if (swap_bytes) swap_double(&wapp->obs_time);
    wapp_obstime=wapp->obs_time;
    fetch_hdrval(h,"header_size",&(wapp->header_size),sizeof(wapp->header_size));
    if (swap_bytes) swap_long(&wapp->header_size);
    wapp_header_size=wapp->header_size;

    /* get the date and start time from the header and convert to MJD */
    fetch_hdrval(h,"obs_date",&(wapp->obs_date),sizeof(wapp->obs_date));
    sscanf(wapp->obs_date,"%4d%2d%2d",&year,&mo,&dd);
    fetch_hdrval(h,"start_time",&(wapp->start_time),sizeof(wapp->start_time));
    sscanf(wapp->start_time,"%d:%d:%d",&hh,&mm,&ss);
    tstart=mjd(year,mo,dd)+hh/24.0+mm/1440.0+ss/86400.0;

    /* for data between April 17 and May 8 inclusive, the start times are
       off by 0.5 days due to the ntp daemon not running... fix here.
       this also occured on May 17! hopefully will not happen again... */
    if ( ((tstart >= 52016.0) && (tstart <= 52039.0)) 
	 || (floor(tstart) == 52046.0)) {
      fprintf(stderr,"WARNING: MJD start time off by 0.5 days! fixed...\n");
      tstart-=0.5;
    }

    fetch_hdrval(h, "samp_time", &(wapp->samp_time),sizeof(wapp->samp_time));
    if (swap_bytes) {
      swap_double(&wapp->wapp_time);
      swap_double(&wapp->samp_time);
    }
    tsamp=(wappcorrect(tstart)+(wapp->wapp_time))*1.0e-6;   
    fetch_hdrval(h,"num_lags",&(wapp->num_lags),sizeof(wapp->num_lags));
    nchans=wapp->num_lags;
    if (swap_bytes) swap_int(&nchans);
    if (swap_bytes) swap_long(&wapp->num_lags);

    /* get scan number from the file name's last four digits */
    scan_number=atoi(&filename[strlen(filename)-4]);
    if (strstr(filename,"wapp1") != NULL) wapp_number=1;
    if (strstr(filename,"wapp2") != NULL) wapp_number=2;
    if (strstr(filename,"wapp3") != NULL) wapp_number=3;
    if (strstr(filename,"wapp4") != NULL) wapp_number=4;
    /*
    fetch_hdrval(h,"scan_number",&(wapp->scan_number),
		 sizeof(wapp->scan_number));
    if (swap_bytes) swap_long(&wapp->scan_number);
    scan_number=(int) wapp->scan_number;
    */

    wapp_isalfa=0;
    fetch_hdrval(h,"isalfa",&(wapp->isalfa),sizeof(wapp->isalfa));
    i=wapp->isalfa;
    if (swap_bytes) swap_int(&i);
    wapp_isalfa=i;

    alfa_ang=0.0;
    fetch_hdrval(h,"alfa_ang",&(wapp->alfa_ang),sizeof(wapp->alfa_ang));
    alfa_ang=wapp->alfa_ang;
    if (swap_bytes) swap_double(&alfa_ang);

    fetch_hdrval(h,"alfa_raj",&(wapp->alfa_raj),sizeof(wapp->alfa_raj));
    fetch_hdrval(h,"alfa_decj",&(wapp->alfa_decj),sizeof(wapp->alfa_decj));
    for (i=0;i<7;i++) {
      alfa_raj[i]=wapp->alfa_raj[i];
      alfa_dej[i]=wapp->alfa_decj[i];
      if (swap_bytes) {
	swap_double(&alfa_raj[i]);
	swap_double(&alfa_dej[i]);
      }
    }

    fetch_hdrval(h,"nifs",&(wapp->nifs),sizeof(wapp->nifs));
    nifs=wapp->nifs;
    if (swap_bytes) swap_int(&wapp->nifs);
    if (swap_bytes) swap_int(&nifs);

    /* set a global to show 3/9 level sampling */
    fetch_hdrval(h,"level",&(wapp->level),sizeof(wapp->level));
    if (swap_bytes) swap_int(&wapp->level); 
    wapp_level=(wapp->level == 1) ? 3 : 9;  

    /* determine number of bits */
    fetch_hdrval(h,"lagformat",&(wapp->lagformat),sizeof(wapp->lagformat));
    if (swap_bytes) swap_int(&wapp->lagformat);
    switch (wapp->lagformat) {
    case 0:
      nbits=16;
      break;
    case 1:
      nbits=32;
      break;
    case 3: /* timing mode data - not relevant, but needs to work! */
      break;
    case 4:
      nbits=8;
      break;
    default:
      error_message("lagformat variable in header should be 0, 1 or 4");
      break;
    }

    /* save lag truncation for use in wapp2fb later on */
    fetch_hdrval(h,"lagtrunc",&(wapp->lagtrunc),sizeof(wapp->lagtrunc));
    if (swap_bytes) swap_int(&wapp->lagtrunc);
    wapp_lagtrunc = wapp->lagtrunc;

    fetch_hdrval(h,"bandwidth",&(wapp->bandwidth),sizeof(wapp->bandwidth));
    fetch_hdrval(h,"freqinversion",&(wapp->freqinversion),
		 sizeof(wapp->freqinversion));
    if (swap_bytes) {
      swap_double(&wapp->bandwidth);
      swap_int(&wapp->freqinversion);
    }

    foff=-1.0*wapp->bandwidth/nchans;
    if (!wapp->freqinversion) {
      wapp_flip=1;
    } else {
      wapp_flip=0;
    }
    if (invert_band) wapp_flip=!wapp_flip; /* force user-specified inversion */

    fetch_hdrval(h,"src_ra",&(wapp->src_ra),sizeof(wapp->src_ra));
    fetch_hdrval(h,"src_dec",&(wapp->src_dec),sizeof(wapp->src_dec));
    fetch_hdrval(h,"start_az",&(wapp->start_az),sizeof(wapp->start_az));
    fetch_hdrval(h,"start_za",&(wapp->start_za),sizeof(wapp->start_za));
    if (swap_bytes) {
      swap_double(&wapp->src_ra);
      swap_double(&wapp->src_dec);
      swap_double(&wapp->start_az);
      swap_double(&wapp->start_za);
    }
    src_raj=wapp->src_ra;
    src_dej=wapp->src_dec;
    az_start=wapp->start_az;
    za_start=wapp->start_za;
    angle_split(src_raj,&rah,&ram,&ras);
    angle_split(src_dej,&ded,&dem,&des);
    cel2gal(rah,ram,ras,ded,dem,des,&srcl,&srcb);

    fetch_hdrval(h,"cent_freq",&(wapp->cent_freq),sizeof(wapp->cent_freq));
    if (swap_bytes) swap_double(&wapp->cent_freq);
    fch1=wapp->cent_freq-nchans*foff/2.0;

    fetch_hdrval(h,"lagformat",&(wapp->lagformat),sizeof(wapp->lagformat));
    fetch_hdrval(h,"sum",&(wapp->sum),sizeof(wapp->sum));
    wapp_lagformat=wapp->lagformat; 
    if (swap_bytes) {
      swap_int(&wapp->lagformat);
      swap_int(&wapp_lagformat);
      swap_int(&wapp->sum); 
    }
    wapp_sum=wapp->sum; 
    wapp_file=h->fd;
    machine_id=2;
    /* get additional header values */

    fetch_hdrval(h, "az", &(wapp->start_az),sizeof(wapp->start_az));
    fetch_hdrval(h, "za", &(wapp->start_za),sizeof(wapp->start_za));
    fetch_hdrval(h, "ast", &(wapp->start_ast),sizeof(wapp->start_ast));
    fetch_hdrval(h, "lst", &(wapp->start_lst),sizeof(wapp->start_lst));
    if (swap_bytes) {
      swap_double(&wapp->start_az);
      swap_double(&wapp->start_za);
      swap_double(&wapp->start_ast);
      swap_double(&wapp->start_lst);
      az_start=wapp->start_az;
      za_start=wapp->start_za;
      ast0=wapp->start_ast;
      lst0=wapp->start_lst;
    }
    fetch_hdrval(h,"project_id",&(wapp->project_id),sizeof(wapp->project_id));
    strcpy(project,wapp->project_id);

    fetch_hdrval(h,"observers",&(wapp->observers),sizeof(wapp->observers));
    strcpy(culprits,wapp->observers);

    fetch_hdrval(h,"psr_dm",&(wapp->psr_dm),sizeof(wapp->psr_dm));
    if (swap_bytes) swap_double(&(wapp->psr_dm)); 

    fetch_hdrval(h,"dumptime",&(wapp->dumptime),sizeof(wapp->dumptime));
    if (swap_bytes) swap_double(&(wapp->dumptime));

    telescope_id=1;

    /* get number of bins which will be non zero in folding mode */
    fetch_hdrval(h,"nbins",&(wapp->nbins),sizeof(wapp->nbins));
    if (swap_bytes) swap_int(&(wapp->nbins));
    fetch_hdrval(h,"isfolding",&(wapp->isfolding),sizeof(wapp->isfolding));
    if (swap_bytes) swap_int(&(wapp->isfolding));

    /* get polyco information */

    fetch_hdrval(h,"rphase",&(wapp->rphase),sizeof(wapp->rphase));
    if (swap_bytes) for (i=0;i<16;i++) swap_double(&wapp->rphase[i]);
    fetch_hdrval(h,"psr_f0",&(wapp->psr_f0),sizeof(wapp->psr_f0));
    if (swap_bytes) for (i=0;i<16;i++) swap_double(&wapp->psr_f0[i]);
    fetch_hdrval(h,"poly_tmid",&(wapp->poly_tmid),sizeof(wapp->poly_tmid));
    if (swap_bytes) for (i=0;i<16;i++) swap_double(&wapp->poly_tmid[i]);
    fetch_hdrval(h,"coeff",&(wapp->coeff),sizeof(wapp->coeff));
    if (swap_bytes) for (i=0;i<192;i++) swap_double(&wapp->coeff[i]);
    fetch_hdrval(h,"num_coeffs",&(wapp->num_coeffs),sizeof(wapp->num_coeffs));
    if (swap_bytes) for (i=0;i<16;i++) swap_int(&wapp->num_coeffs[i]);

    if (wapp->isfolding) {
      mjdobs=floor(tstart);
      tstart=(tstart-mjdobs)*86400.0;
      return 4; /*folding mode*/
    } else {
      return 3; /*fast-sampled mode*/
    }
  }


  /* Try reading as AOFTM fast-sampled data      -> case 5 */
  rewind(fptr);
 if (aoftm_read_header(filename)) {
    machine_id=3;
    telescope_id=1;
    data_type=1;
    nifs=1;
    obits=8;
    return 7;
  }

  /* Try reading as GMRT fast-sampled data -> case 10 */
  rewind(fptr);
  for (i=1;i<=3;i++) fread(&gmrt[i],2,1,fptr);
  /*fprintf(stderr,"%d %d %d\n",gmrt[1],gmrt[2],gmrt[3]);exit(0);*/
  if ((gmrt[2]==32768)) { 
  }
    rewind(fptr);
    machine_id=7;
    telescope_id=7;
    data_type=1;
    fch1=626.031250;
    foff=-0.0625;
    nchans=256;
    nbits=16;
    tsamp=0.000256;
    nifs=1;
    return(10);

  /* Try reading as OOTY fast-sampled data       -> case 8 */
  rewind(fptr);
  machine_id=5;
  telescope_id=2;
  data_type=1;
  foff=-0.03125;
  fch1=330.5+foff/2.0;
  nchans=256;
  nbits=1;
  tstart=0.0;
  tsamp=0.00025800387;
  nifs=1;
  return 8;

  /* unknown data type */
  return 0;
}
