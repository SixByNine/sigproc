#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "filterbank.h"

#define HEADERRECS 36
#define HEADERSIZE HEADERRECS*80 		/* FITS 2880-byte block */
#define AOFTMBLEN 524288	/* buffer size in bytes */
#define CHANSPERQ 256		/* n chans per freq quadrant */

/* data management variables */
unsigned short buffer[AOFTMBLEN/2];
unsigned short BRO[256], BRO2[256], FOLD[256], FOLD2[256];
unsigned short aoftm_data[1024];
double quadif[4];
int quadsign[4];
int nfft, nbufs, nbits, nchans, twobit; 

char aoftm_header[HEADERRECS][80]; 
double aoftm_if[4];
int aoftm_sgn[4], aoftm_pol[16];

char *headername (char *filename) /*includefile*/
{
  char *string, *hdrnam;
  string=(char *)malloc(strlen(filename));
  hdrnam=(char *)malloc(strlen(filename));
  string=strtok(filename,".");
  while (string != NULL) {
    if (strcmp(string,"dat")) {
      strcat(hdrnam,string);
      strcat(hdrnam,".");
    } else {
      strcat(hdrnam,"hdr");
    }
    string=strtok(NULL,".");
  }
  free(string);
  return(hdrnam);
}

int aoftm_read_header(char *filename) /*includefile*/
{
  FILE *out_hdr, *board_config, *instream;
  char outfile[100];

  /* header info variables */
  char aoftm_polcfg[33], aoftm_chord[9];
  double freq, bw, dt, dt_max, int_time;
  char utc_start[8];
  char caltype[6], tick, dest, source[9];
  char hdrfname[80], runtype[6], receiver[9], projid[9];
  int rah, ram, decd, decm, day, month, year, hr, min, sec;
  float ras, decs;
  double glon, glat;
  int iboard, iquad, icount, nchan_max, boardno, test_vec;
  int seek, equinox;
  
  double fmjd, az, za;
  int j, mjd;
  int nscan;

  strcpy(hdrfname,"");
  strncat(hdrfname,filename,strlen(filename)-4);
  strcat(hdrfname,".hdr");

  if ((instream = fopen(hdrfname, "r")) == NULL) {
    /*fprintf(stderr,"Couldn't open pipe for header read!\n");*/
    return(0);
  }
  if (fread(aoftm_header, 1, HEADERSIZE, instream) != HEADERSIZE) {
    fprintf(stderr, "Unable to read header!\n");
    return(0);
  }
  close(instream);
  az=za=0.0;
  for (j=0; j<HEADERRECS; j++) {
    if (strstr(aoftm_header[j], "MJD-OBS") != NULL) 
      sscanf(aoftm_header[j], "MJD-OBS = %lf", &fmjd);
    if (strstr(aoftm_header[j], "AZ-OBS") != NULL)
      sscanf(aoftm_header[j], "AZ-OBS = %lf", &az);
    if (strstr(aoftm_header[j], "EL-OBS") != NULL)
      sscanf(aoftm_header[j], "EL-OBS = %lf", &za);
    if (strstr(aoftm_header[j], "EQUINOX") != NULL)
      sscanf(aoftm_header[j], "EQUINIX = %d", &equinox);
    if (strstr(aoftm_header[j], "SOURCE") != NULL)
      sscanf(aoftm_header[j], "SOURCE  = '%9c'",&source);
    if (strstr(aoftm_header[j], "AZ-OBS") != NULL)
      sscanf(aoftm_header[j], "AZ-OBS = '%lf'",&az);
    if (strstr(aoftm_header[j], "EL-OBS") != NULL)
      sscanf(aoftm_header[j], "EL-OBS = '%lf'",&za);
    if (strstr(aoftm_header[j], "DATE-OBS") != NULL)
      sscanf(aoftm_header[j], "DATE-OBS= ' %2d/%2d/%2d", &day,&month,&year);
    if (strstr(aoftm_header[j], "TIME-OBS") != NULL)
      sscanf(aoftm_header[j], "TIME-OBS= ' %2d:%2d:%2d", &hr,&min,&sec);
    if (strstr(aoftm_header[j], "RA-OBS") != NULL)
      sscanf(aoftm_header[j], "RA-OBS  = ' %2d:%2d:%.1f",&rah,&ram,&ras);
    if (strstr(aoftm_header[j], "DEC-OBS") != NULL)
      sscanf(aoftm_header[j], "DEC-OBS = ' %3d:%2d:%.0f",&decd,&decm,&decs);
    if (strstr(aoftm_header[j], "NFFT") != NULL)
      sscanf(aoftm_header[j], "NFFT = %d", &nfft);
    if (strstr(aoftm_header[j], "NINT") != NULL)
      sscanf(aoftm_header[j], "NINT = %d", &icount);
    if (strstr(aoftm_header[j], "BOARD") != NULL)
      sscanf(aoftm_header[j], "BOARD = %d", &boardno);
    if (strstr(aoftm_header[j], "T-VECTOR") != NULL)
      sscanf(aoftm_header[j], "T-VECTOR = %d", &test_vec);
    if (strstr(aoftm_header[j], "BUFFERS") != NULL)
      sscanf(aoftm_header[j], "BUFFERS = %d", &nbufs);
    if (strstr(aoftm_header[j], "NBITS") != NULL)
      sscanf(aoftm_header[j], "NBITS = %d", &nbits);
    if (strstr(aoftm_header[j], "DELTA-T") != NULL)
      sscanf(aoftm_header[j], "DELTA-T = %lf", &dt);
    if (strstr(aoftm_header[j], "TOTAL-T") != NULL)
      sscanf(aoftm_header[j], "TOTAL-T = %lf", &int_time);
    if (strstr(aoftm_header[j], "FREQ") != NULL)
      sscanf(aoftm_header[j], "FREQ = %lf", &freq);
    if (strstr(aoftm_header[j], "BW") != NULL)
      sscanf(aoftm_header[j], "BW = %lf", &bw);
    if (strstr(aoftm_header[j], "RECEIVER") != NULL)
      sscanf(aoftm_header[j], "RECEIVER='%9c'",&receiver);
    if (strstr(aoftm_header[j], "PROJECT") != NULL)
      sscanf(aoftm_header[j], "PROJECT ='%9c'",&projid);
    if (strstr(aoftm_header[j], "POL_CFG") != NULL)
      sscanf(aoftm_header[j],"POL_CFG = '%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d'",&aoftm_pol[0],&aoftm_pol[1],&aoftm_pol[2],&aoftm_pol[3],&aoftm_pol[4],&aoftm_pol[5],&aoftm_pol[6],&aoftm_pol[7],&aoftm_pol[8],&aoftm_pol[9],&aoftm_pol[10],&aoftm_pol[11],&aoftm_pol[12],&aoftm_pol[13],&aoftm_pol[14],&aoftm_pol[15]);
    if (strstr(aoftm_header[j], "IF_CFG") != NULL)
      sscanf(aoftm_header[j],"IF_CFG = '%lf %lf %lf %lf",&aoftm_if[0], 
	     &aoftm_if[1], &aoftm_if[2], &aoftm_if[3]);
    if (strstr(aoftm_header[j], "CHN_ORDR") != NULL)
      sscanf(aoftm_header[j],"CHN_ORDR = '%d %d %d %d'",&aoftm_sgn[0],&aoftm_sgn[1],&aoftm_sgn[2],&aoftm_sgn[3]);
    if (strstr(aoftm_header[j], "10-SEC") != NULL)
      sscanf(aoftm_header[j], "10-SEC = %c", &tick);
    if (strstr(aoftm_header[j], "CAL-TYPE") != NULL)
      sscanf(aoftm_header[j], "CAL-TYPE = '%6s'", caltype);
    if (strstr(aoftm_header[j], "SCAN-NO") != NULL)
      sscanf(aoftm_header[j], "SCAN-NO = %d", &nscan);
    if (strstr(aoftm_header[j], "SEEK-BGN") != NULL)
      sscanf(aoftm_header[j], "SEEK-BGN = %d", &seek);
    if (strstr(aoftm_header[j], "DESTINAT") != NULL)
      sscanf(aoftm_header[j], "DESTINAT = ' %c", &dest);
  }

  mjd = (int) fmjd;
  twobit = (nbits == 2) ? 1 : 0;
  nchans = (twobit) ? 4*nfft : nfft;
  strcpy(source_name,"");
  strncat(source_name,source,9);
  if (ras<0.0) ras=0.0;
  if (decs<0.0) decs=0.0;
  src_raj=(double)rah*10000.0+(double)ram*100.0+(double)ras;
  src_dej=abs(decd)*10000.0+decm*100.0+decs;
  if (decd<0.0) src_dej*=-1.0;
  cel2gal(rah,ram,ras,decd,decm,decs,&glon,&glat);
  az_start=az;
  za_start=za;
  tsamp=dt*1.0e-6;
  tstart=fmjd;
  fch1=freq+5115.0/1024.0; /* to get to the top of the band */
  foff=bw/(double)nchans;
  return(1);
}
void create_bro(int nfft) /*includefile*/
{
  int i,j,k;
 
  k = rint( log((double)nfft)/log(2.0) ); 
  for (i=0; i<nfft; i++) {
    BRO[i] = 0;
    for (j=0; j<k; j++)
      BRO[i] += ((i >> j) & 1) << (k-j-1);
    if (i < nfft/2)
      FOLD[nfft/2 - i - 1] = BRO[i];
    else
      FOLD[3*nfft/2 - i - 1] = BRO[i];
  }
 
  if (nfft < 256) {
    for (i=0; i<nfft; i++) {
      FOLD[i] += 1;
      if (FOLD[i] > (unsigned short)(nfft - 1)) FOLD[i] = 0;
    }
  }
 
  for (i=0; i<nfft; i++) {
    FOLD2[FOLD[i]] = i;
    BRO2[BRO[i]] = i;
  }

}
void two_bit_reorder(unsigned short *spectrum, int nfft) /*includefile*/
{
  int i;
 
  for(i=0; i<nfft; i += 2){         
    aoftm_data[FOLD2[i]]          = ((*(spectrum + i/2) >> 6) & 0x0003);
    aoftm_data[FOLD2[i+1]]        = ((*(spectrum + i/2) >> 14) & 0x0003);

    aoftm_data[FOLD2[i]+nfft]     = ((*(spectrum + i/2) >> 4) & 0x0003);
    aoftm_data[FOLD2[i+1]+nfft]   = ((*(spectrum + i/2) >> 12) & 0x0003);

    aoftm_data[FOLD2[i]+nfft*2]   = ((*(spectrum + i/2) >> 2) & 0x0003);
    aoftm_data[FOLD2[i+1]+nfft*2] = ((*(spectrum + i/2) >> 10) & 0x0003);

    aoftm_data[FOLD2[i]+nfft*3]   = ((*(spectrum + i/2) >> 0) & 0x0003);
    aoftm_data[FOLD2[i+1]+nfft*3] = ((*(spectrum + i/2) >> 8) & 0x0003);
  }



}
void eight_bit_reorder(unsigned short *spectrum, int nfft) /*includefile*/
{
  int i;
 
  for (i=0; i<nfft; i += 2){
    aoftm_data[FOLD2[i]]   = (*(spectrum + i/2) & 0x00ff);
    aoftm_data[FOLD2[i+1]] = (*(spectrum + i/2) >> 8);
  }
}
void aoftm2fb(char *filename, FILE *output) /*includefile*/
{ 
  unsigned char data_out[1024],flip_out[1024];
  unsigned short *p, *p0;
  int isam,ibuf,index,rr,i,j,doit;
  unsigned long idx[4*CHANSPERQ];
  FILE *instream;
  float chanif[4*CHANSPERQ], quadbw, df, realtime;
  int nskip,iquad, ichan, iboard, swap_bytes, samples_per_buffer;
  char message[80];

  samples_per_buffer=4194304/nchans/nbits;

  nskip=0;
  twobit = (nbits == 2) ? 1 : 0;
  nchans = (twobit) ? 4*nfft : nfft;

/* rr is the number of AOFTMBLEN buffers of data to be read and stored on disk 
(un-bit-reversed but packed); */

  rr = (int) nbufs;     /* or buffer numbers */
  rr = (nbufs < rr) ? nbufs : rr;

  /* open up logfile */
  open_log("filterbank.monitor");

  create_bro(nfft);   
  swap_bytes=little_endian();

  instream=open_file(filename,"rb");

  quadif[0]=50.0;quadif[1]=52.5;quadif[2]=55.0;quadif[3]=57.5;
  quadsign[0]=quadsign[1]=quadsign[2]=quadsign[3]=-1;

  if (twobit) {
    quadbw = quadif[0]-quadif[1];
    if (quadbw < 0.0) quadbw = -quadbw;
    df = quadbw/(double)CHANSPERQ;
    for (iquad=0;iquad<4;iquad++) {
      for (j=0;j<CHANSPERQ;j++) {
	ichan = j + iquad*CHANSPERQ;
	chanif[ichan] = quadif[iquad]
	  + quadsign[iquad]*((double)j*df - (quadbw-df)/2.);
      }
    }
    indexx(1024,&(chanif[-1]),&(idx[-1]));
    for (ichan=0;ichan<nchans;ichan++) idx[ichan]--;
    
  }
  else {
    for (ichan=0;ichan<CHANSPERQ;ichan++) {
      if (quadsign[1] < 0.) 
	idx[ichan] = nchans-ichan-1; 
      else 
	idx[ichan] = ichan;
    }
  }
  
  /* loop over rr buffers, un-bit-reversing one sample at a time, and 
     writing packed data to disk */

  for (ibuf=0; ibuf<rr; ibuf++) {
    realtime=(float)tsamp*(float)samples_per_buffer*(float)ibuf;
    if ( (doit=process(realtime,start_time,final_time)) == -1) break;
      if (fread(buffer, 2, AOFTMBLEN/2, instream) != AOFTMBLEN/2) {
	/*fprintf(stderr, "End of data at buffer %d\nbuf", ibuf+1);*/
	rr = ibuf + 1;
	break;          /* buffers are 2 byte values */
      }
    if (doit) {
      sprintf(message,"time:%.1fs",realtime);
      update_log(message);
      if (swap_bytes) 
	for (isam=0;isam<AOFTMBLEN/2;isam++) swap_short(&buffer[isam]);
      
      for (isam = 0; isam<AOFTMBLEN/nfft; isam++) {
	p = &(buffer[isam*nfft/2]);
	if (twobit) {
	  two_bit_reorder(p,nfft);
	  for (ichan=0; ichan<nchans; ichan++) {
	    index = idx[ichan];
	    if (aoftm_data[index] == 1) 
	      aoftm_data[index] = 3;
	    else if (aoftm_data[index] == 0) 
	      aoftm_data[index] = 2;
	    else if (aoftm_data[index] == 3) 
	      aoftm_data[index] = 1;
	    else 
	      aoftm_data[index] = 0;
	    data_out[ichan] = (unsigned char) aoftm_data[index];
	  }
	} else {
	  eight_bit_reorder(p,nfft); 
	  for (ichan=0; ichan<nchans; ichan++) {
	    index = idx[ichan];
	    data_out[ichan] = (unsigned char) aoftm_data[index];
	  }
	}
	j=nchans;
	for (i=0;i<nchans;i++) flip_out[i]=data_out[--j]; 
	if (fwrite(flip_out, 1, nchans, output) != nchans) {
	  fprintf(stderr,"Couldn't write AOFTM data!\n");
	  close(instream);
	  exit(-1);
	}
      }
    }
  }
  close(instream);
}



