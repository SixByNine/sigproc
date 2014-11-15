#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fitsio.h>

#define MX_NCHAN 4096

// Routine to read the header from a PSRFITS file 
void readpsrfits_hd(char *filename,int *machine_id,int *telescope_id,int *data_type,char *source_name,
		    double *fch1,int *nchans,double *foff,int *nifs,double *tsamp,
		    double *src_raj,double *src_dej,double *tstart,int *nbits,int *ibeam)
{
  fitsfile *fp;
  int status;
  float obsbw,fc;
  char sobsbw[100],sfc[100],ras[100],decs[100];
  char name[100],telescope[100];
  double h,m,sec,chanbw,smjd,offs;
  int imjd;

  int ncol, NCHANS, anynul;
  unsigned char nulval = 0;
  long long row, firstelem, nelem;
  float ch_freq[MX_NCHAN];

  // Defaults that are not being set
  *machine_id=-1;
  *data_type = 0;
  *ibeam = 1;
  status=0;

  fits_open_file(&fp,filename,READONLY,&status);
  fits_report_error(stderr,status);
  fits_movabs_hdu( fp, 1, NULL, &status );
  fits_read_key(fp, TSTRING, "SRC_NAME", source_name, NULL, &status);
  if(strlen(source_name)==0){
	  strcpy(source_name,"UNKNOWN      ");
  }
  fits_report_error(stderr,status);
  if (fits_read_key(fp, TINT, "IBEAM", ibeam, NULL, &status))
    {
      *ibeam = 1;
      status=0;
    }
  fits_report_error(stderr,status);
  //fits_read_key(fp, TINT, "OBSNCHAN", nchans, NULL, &status);
  //fits_report_error(stderr,status);
  fits_read_key(fp, TFLOAT, "OBSBW", &obsbw, NULL, &status);
  //  if (obsbw > 0) obsbw = -obsbw; // We must invert the data
  // The bandwidth in FITS is always positive
  // In order to determine whether sigproc should be given a negative bandwidth
  // the sign of the CHANBW parameter needs to be determined
  // This is done below
  fits_report_error(stderr,status);
  fits_read_key(fp, TFLOAT, "OBSFREQ", &fc, NULL, &status);
  fits_report_error(stderr,status);
  fits_read_key(fp, TINT, "STT_IMJD", &imjd, NULL, &status);
  fits_report_error(stderr,status);
  fits_read_key(fp, TDOUBLE, "STT_SMJD", &smjd, NULL, &status);
  fits_report_error(stderr,status);
  fits_read_key(fp, TDOUBLE, "STT_OFFS", &offs, NULL, &status);
  fits_report_error(stderr,status);
  fits_read_key(fp, TSTRING, "TELESCOP", telescope, NULL, &status);
  fits_report_error(stderr,status);
  if (strcasecmp(telescope,"PARKES")==0)
    *telescope_id=4;
  else if (strcasecmp(telescope,"ARECIBO")==0)
    *telescope_id=1;
  else if (strcasecmp(telescope,"JODRELL")==0)
    *telescope_id=5;
  else if (strcasecmp(telescope,"GBT")==0)
    *telescope_id=6;
  else if (strcasecmp(telescope,"EFFELSBERG")==0)
    *telescope_id=8;
  else if (strcasecmp(telescope,"SRT")==0)
    *telescope_id=10;
  else 
    *telescope_id = -1;

  // Start time
  *tstart = imjd+(smjd+offs)/86400.0;

  fits_read_key(fp, TSTRING, "RA", ras, NULL, &status);
  sscanf(ras,"%lf:%lf:%lf",&h,&m,&sec);
  *src_raj  = (double)h*10000.0 + (double)m*100.0 + (double)sec;

  fits_read_key(fp, TSTRING, "DEC", decs, NULL, &status);
  sscanf(decs,"%lf:%lf:%lf",&h,&m,&sec);
  
  *src_dej  = (double)fabs(h)*10000.0 + (double)m*100.0 + (double)sec;
  if (h < 0) (*src_dej)*=-1;

  
  // Now get information from the subint header
  strcpy(name,"SUBINT");
  fits_movnam_hdu(fp, BINARY_TBL, name, 0, &status);
  fits_read_key(fp, TDOUBLE, "TBIN", tsamp, NULL, &status );
  fits_read_key(fp, TINT, "NPOL", nifs, NULL, &status );
  fits_read_key(fp, TINT, "NBITS", nbits, NULL, &status );
  fits_read_key(fp, TDOUBLE, "CHAN_BW", &chanbw, NULL, &status );
//  if (chanbw < 0) obsbw=-(obsbw);

  fits_movnam_hdu( fp, BINARY_TBL, "SUBINT", 0, &status );
  fits_read_key( fp, TINT, "NCHAN", &NCHANS, NULL, &status);
  fits_get_colnum( fp, CASESEN, "DAT_FREQ", &ncol, &status );

  firstelem = 1LL;
  nelem = NCHANS;
  *nchans = NCHANS;
  row = 1LL;

  fits_read_col( fp, TFLOAT, ncol, row, firstelem, nelem, &nulval,
		  ch_freq, &anynul, &status );

  //*fch1 = ch_freq[0] > ch_freq[nelem-1] ? ch_freq[0] : ch_freq[nelem-1];
  //*fch1 = fc-obsbw/2.0 + obsbw/(*nchans)/2.0;
  //*fch1 = fc - ((*nchans-1)/ 2.0)*chanbw;
  *fch1=ch_freq[0];
  *foff = chanbw;

  fprintf(stderr,"Native fch1=%f, foff=%f\n",*fch1,*foff);
  if(chanbw > 0){
	 fprintf(stderr,"Will flip the band!\n");
  }
  
  fits_close_file(fp,&status);
  fits_report_error(stderr,status);
}
