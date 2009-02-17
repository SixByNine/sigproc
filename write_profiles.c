#include "fold.h"
/* 
   write_profiles.c:

   this subroutine does some final user-defined processing (scaling, baseline
   subtraction and IF summing) before writing out the folded pulse profiles 
   in a number of different formats. The time stamp and folding period in the 
   profiles refers to the START of a given integration. This routine is 
   normally called from within fold_data.c to dump profiles after a pre-defined
   subintegration. 

   Modification history:

   April 2, 2003 - added IF summing and checked time stamps (drl@jb.man.ac.uk)
   June 22, 2004 - added psrfits option (drl@jb.man.ac.uk)
*/
#define MAX_BLKS 20
#define MAX_COEFF 15
double tsta, tsec, tadd, pfld0, window;
long int pulsecount;
void write_profiles(float *prof,int nbins, int nchan, int nifs, FILE *out)/*includefile*/
{
  int i,j,c,b,k,m,year,month,day;
  float *profile,scale,offset,rms,epoch;
  double usec,tbin,tres,psec;
  static int first=1;
  struct EPN epn;

  int sta=0;
  short int sj;
  float x,binned_freq[16384],binned_weight[16384],
    binned_offset[16384],binned_scale[16384],binned_data[16384];
  int bitpix=8, naxis=0, nrcvr, nrows, ncols, col;
  static int subint_cnt;
  char datestr[10],timestr[8];
  char Cstr16[16], Estr16[16], Istr16[16];
  char *ttype[20], *tform[20], *tunit[20];
  long naxes[4];
  int junk,rah,ram,ded,dem,subint_hdu,hh,mm;
  float ss;
  int last_scanhdr_hdu;
  double ras,des,dx;
  char rastr[80], destr[80], sra[80], sde[80];
  char date_time[24];
  char *pch[MAX_BLKS];
  char site[MAX_BLKS][2];
  short int nspan[MAX_BLKS];
  short int ncoeff[MAX_BLKS];
  double rfreq[MAX_BLKS];
  double rmjd[MAX_BLKS];
  double rphase[MAX_BLKS];
  double lgfiterr[MAX_BLKS];
  double f0[MAX_BLKS];
  double coeff[MAX_BLKS][MAX_COEFF];
  static double srcl, srcb;
  /* For Pulsar History BINTABLE */
 
  static char *PHtype[18] = {
    "DATE_PRO","PROC_CMD","POL_TYPE","NPOL    ","NBIN    ","NBIN_PRD","TBIN    ","CTR_FREQ",
    "NCHAN   ","CHAN_BW ","PAR_CORR","RM_CORR ","DEDISP  ","DDS_MTHD","SC_MTHD ","CAL_MTHD",
    "CAL_FILE","RFI_MTHD"
  };
  
  static char *PHform[18] = {
    "24A     ","80A     ","8A      ","1I      ","1I      ","1I      ","1D      ","1D      ",
    "1I      ","1D      ","1I      ","1I      ","1I      ","32A     ","32A     ","32A     ",
    "32A     ","32A     "
  };
  
  static char *PHunit[18] = {
    "        ","        ","        ","        ","        ","        ","s       ","MHz     ",
    "        ","MHz     ","        ","        ","        ","        ","        ","        ",
    "        ","        "
  };

  /* subtract baseline from outgoing profiles if requested */
  if (userbase != 0.0) for (i=0;i<multiple*nbins*nifs*nchans;i++) 
    prof[i]-=userbase;

  /* multiply outgoing profiles by Jansky calibration factor if supplied */
  if (jyfactor != 0.0) for (i=0;i<multiple*nbins*nifs*nchans;i++) 
    prof[i]*=jyfactor;

  /* sum first two polarizations together if requested */
  if (totalpower && (nifs > 1)) {
    for (i=0;i<multiple*nbins*nchans;i++) 
      prof[i]+=prof[i+multiple*nbins*nchans];
    nifs=1;
  }

  if (binary) {
    /* write out profiles in binary format */
    folding_period=pfld0*1000.0;
    tadd=tstart;
    tstart=tstart+tsta/86400.0;
    npuls=pulsecount;
    fold_header();
    tstart=tadd;
    for (i=0; i<nifs; i++) {
      for (c=0; c<nchan; c++) {
	for (b=0; b<nbins; b++) 
	  fwrite(&prof[i*nchan*nbins+c*nbins+b],sizeof(float),1,out);
      }
    }
  } else if (stream) {
    /* write out profiles as ASCII streams with START/STOP boundaries */
    k=0;
    for (i=0;i<nifs;i++) {
      for (c=nchan-1;c>=0;c--) {
	fprintf(out,"#START %d %f %f\n",nbins,tsta,fch1+foff*(float)c);
	for (m=0;m<multiple;m++) {
	  for (b=0;b<nbins;b++) {
	    fprintf(out,"%d %f\n",b+m*nbins,prof[i*nchan*nbins+c*nbins+b]);
	  }
	}
	fprintf(out,"#STOP\n");
      }
    }
  } else if (asciipol) {
    /* write profiles in format for Jim's polarization code */
    for (b=0;b<nbins;b++) 
      for (i=0;i<nifs;i++)
	for (c=0;c<nchan;c++) 
	  fprintf(out,"%d %d %d %f\n",b,i,c,prof[i*nchan*nbins+c*nbins+b]);
  } else if (ascii) {
    fprintf(output,"# %.1f %.7f %.10f %ld %.3f %.3f %d %c %d %s\n",
     floor(tstart),(tstart-floor(tstart))*86400.0+tsta,pfld0,pulsecount,fch1,refdm,nbins,tempo_site(telescope_id),1,source_name);
    for (b=0;b<nbins;b++) {
      fprintf(out,"%d",b+1);
      for (i=0;i<nifs;i++) {
	for (c=0;c<nchan;c++) fprintf(out," %f",prof[i*nchan*nbins+c*nbins+b]);
      }
      fprintf(out,"\n");
    }
#ifdef PSRFITS
  } else if (psrfits) {
    if (first) {
      first=0;
      /* write profile in PSRFITS format */
      fits_create_file(&fits, "stdout", &sta);
      fits_create_img(fits,bitpix,naxis,naxes,&sta);
      fits_write_date(fits,&sta);
      /* Get DateTime string - required in various BINTABLEs */
      fits_get_system_time(date_time,&junk,&sta);
      fits_update_key(fits,TSTRING,"HDRVER","1.19","Header version",&sta);
      fits_update_key(fits,TSTRING,"OBSERVER",culprits,
		      "Observer name(s)",&sta);
      fits_update_key(fits,TSTRING,"PROJID",project,"Project name",&sta);
      fits_update_key(fits,TSTRING,"TELESCOP",telescope_name(telescope_id),
		      "Telescope name", &sta);
      fits_update_key(fits,TSTRING,"BACKEND",backend_name(machine_id),
		      "Backend ID",&sta);
      if (nifs>1) 
	nrcvr=2;
      else
	nrcvr=1;
      fits_update_key(fits,TINT,"NRCVR",&nrcvr,
		      "Number of receiver channels (I)",&sta);
      fits_update_key(fits,TSTRING,"OBS_MODE", "PSR",
		      "(PSR, CAL, SEARCH)", &sta);
      fits_update_key(fits,TSTRING,"SRC_NAME", source_name,
		      "Source or scan ID", &sta);
      fits_update_key(fits,TSTRING,"COORD_MD", "J2000",
		      "Coordinate mode (J2000, Gal, Ecliptic, etc.)", &sta);
      
      angle_split(src_raj,&rah,&ram,&ras);
      if (ras<10.0) 
	sprintf(sra,"0%.3f",ras);
      else
	sprintf(sra,"%.3f",ras);
      sprintf(rastr,"%02d:%02d:%s",rah,ram,sra);
      angle_split(src_dej,&ded,&dem,&des);
      if (des<10.0) 
	sprintf(sde,"0%.3f",des);
      else
	sprintf(sde,"%.3f",des);
      sprintf(destr,"%02d:%02d:%s",ded,dem,sde);
      cel2gal(rah,ram,ras,ded,dem,des,&srcl,&srcb);

      fits_update_key(fits, TSTRING, "STT_CRD1", rastr,
		      "Start coord 1 (hh:mm:ss.sss or ddd.ddd)", &sta);
      fits_update_key(fits, TSTRING, "STT_CRD2", destr,
		      "Start coord 2 (-dd:mm:ss.sss or -dd.ddd)", &sta);
      fits_update_key(fits, TSTRING, "TRK_MODE", "TRACK",
		      "Track mode (TRACK, SCANGC, SCANLAT)", &sta);
      fits_update_key(fits, TSTRING, "STP_CRD1", rastr,
		      "Stop coord 1 (hh:mm:ss.sss or ddd.ddd)", &sta);
      fits_update_key(fits, TSTRING, "STP_CRD2", destr,
		      "Stop coord 2 (-dd:mm:ss.sss or -dd.ddd)", &sta);
      fits_update_key(fits, TSTRING, "CAL_MODE", "OFF",
		      "Cal mode (OFF, SYNC, EXT1, EXT2)", &sta);
      fits_create_tbl(fits,BINARY_TBL,0,18,PHtype,
		      PHform,PHunit,"HISTORY",&sta);
      pch[0] = date_time;
      fits_write_col(fits,  TSTRING, 1, 1, 1, 1, pch, &sta);
      pch[0] = "SIGPROC";
      fits_write_col(fits,  TSTRING, 2, 1, 1, 1, pch, &sta);
      pch[0] = "??";
      fits_write_col(fits,  TSTRING, 3, 1, 1, 1, pch, &sta);
      /* Nr of pols  - i.e. actually written out */
      sj = nifs;
      fits_write_col(fits, TSHORT, 4, 1, 1, 1, &sj, &sta );
      
      /* Nr of bins per product (0 for SEARCH mode) */
      sj = nbins;
      fits_write_col( fits, TSHORT, 5, 1, 1, 1, &sj, &sta );
      
      /* Nr of bins per period */
      fits_write_col( fits, TSHORT, 6, 1, 1, 1, &sj, &sta );
      
      /* Bin time */
      dx =  folding_period/sj;
      fits_write_col( fits, TDOUBLE, 7, 1, 1, 1, &dx, &sta );
      
      /* Centre freq. */
      dx = ((double)(nchans/2)-1.0)*foff+fch1;
      fits_write_col( fits, TDOUBLE, 8, 1, 1, 1, &dx, &sta );
      
      /* Number of channels */
      sj = nchans;
      fits_write_col( fits, TSHORT, 9, 1, 1, 1, &sj, &sta );
      
      /* Channel bandwidth */
      dx = foff;
      fits_write_col( fits, TDOUBLE, 10, 1, 1, 1, &dx, &sta );
      
      /* Create frequency array for later */
      /* Get freq. of first channel */
      if ((fch1==0.0) && (foff==0.0)) {
	for (i=0; i<nchans; i++) binned_freq[i]=(float) frequency_table[i];
      } else {
	for (i=0; i<nchans; i++) 
	  binned_freq[i]=fch1+i*foff;
      }
      
      for (i=0; i<16384; i++) {
	binned_scale[i]=1.0;
	binned_offset[i]=0.0;
	binned_weight[i]=1.0;
      }
      
      sj = 0;
      /* Parallactic angle correction applied */
      fits_write_col( fits, TSHORT, 11, 1, 1, 1, &sj, &sta );
      /* RM correction applied */
      fits_write_col( fits, TSHORT, 12, 1, 1, 1, &sj, &sta );
      /* Data dedispersed */
      fits_write_col( fits, TSHORT, 13, 1, 1, 1, &sj, &sta );
      /* Dedispersion method */
      pch[0] = "NONE";
      fits_write_col( fits, TSTRING, 14, 1, 1, 1, pch, &sta );
      /* Scattered power correction method */
      pch[0] = "NONE";
      fits_write_col( fits, TSTRING, 15, 1, 1, 1, pch, &sta );
      /* Calibration method */
      pch[0] = "NONE";
      fits_write_col( fits, TSTRING, 16, 1, 1, 1, pch, &sta );
      /* Name of calibration file */
      pch[0] = "NONE";
      fits_write_col( fits, TSTRING, 17, 1, 1, 1, pch, &sta );
      /* RFI excision method */
      pch[0] = "NONE";
      fits_write_col( fits, TSTRING, 18, 1, 1, 1, pch, &sta );
      /* Store last header hdu number */
      fits_get_hdu_num( fits, &last_scanhdr_hdu );

      subint_cnt=1;
      /* Add START TIME to primary HDU */
      /* Move to primary HDU */
      fits_movabs_hdu( fits, 1, NULL, &sta );
      cal(tstart,&epn.year,&epn.month,&epn.day);
      sprintf(datestr,"%4d-%02d-%02d",epn.year,epn.month,epn.day);
      fits_update_key( fits, TSTRING, "STT_DATE", datestr,
		       "Start UT date (YYYY-MM-DD)", &sta);
      uttime(tstart,&hh,&mm,&ss);
      sprintf(timestr,"%02d:%02d:%02d",hh,mm,ss);
      fits_update_key( fits, TSTRING, "STT_TIME", timestr,
		       "Start UT (hh:mm:ss)", &sta);
      dx=86400*(tstart-floor(tstart));
      fits_update_key( fits, TINT, "STT_SMJD", &dx,
		     "[s] Start time (sec past UTC 00h) (J)", &sta);
      /* Move to last created HDU in scan header */
      fits_movabs_hdu( fits, last_scanhdr_hdu, NULL, &sta );

      /* Create SUBINT BINTABLE */
      ttype[0] = "ISUBINT ";/* Subint number. If NAXIS=-1, 0 indicates EOD. */
      tform[0] = "1J      ";
      tunit[0] = "";
      ttype[1] = "INDEXVAL";    /* Optionally used if INT_TYPE != TIME */
      tform[1] = "1D      ";
      tunit[1] = "";
      ttype[2] = "TSUBINT ";    /* [s] Length of subintegration */
      tform[2] = "1D      ";
      tunit[2] = "";
      ttype[3] = "OFFS_SUB"; /* [s] Offset from Start UTC of subint centre */
      tform[3] = "1D      ";
      tunit[3] = "";
      ttype[4] = "LST_SUB ";    /* [s] LST at subint centre */
      tform[4] = "1D      ";
      tunit[4] = "";
      ttype[5] = "RA_SUB  ";    /* [turns] RA (J2000) at subint centre */
      tform[5] = "1D      ";
      tunit[5] = "";
      ttype[6] = "DEC_SUB ";    /* [turns] Dec (J2000) at subint centre */
      tform[6] = "1D      ";
      tunit[6] = "";
      ttype[7] = "GLON_SUB";    /* [deg] Gal longitude at subint centre */
      tform[7] = "1D      ";
      tunit[7] = "";
      ttype[8] = "GLAT_SUB";    /* [deg] Gal latitude at subint centre */
      tform[8] = "1D      ";
      tunit[8] = "";
      ttype[9] = "FD_ANG  ";    /* [deg] Feed angle at subint centre */
      tform[9] = "1E      ";
      tunit[9] = "";
      ttype[10] = "POS_ANG ";/*[deg] Position angle of feed at subint centre */
      tform[10] = "1E      ";
      tunit[10] = "";
      ttype[11] = "PAR_ANG ";    /* [deg] Parallactic angle at subint centre */
      tform[11] = "1E      ";
      tunit[11] = "";
      ttype[12] = "TEL_AZ  ";    /* [deg] Telescope azimuth at subint centre */
      tform[12] = "1E      ";
      tunit[12] = "";
      ttype[13] = "TEL_ZEN ";/*[deg] Telescope zenith angle at subint centre */
      tform[13] = "1E      ";
      tunit[13] = "";

      sprintf( Cstr16, "%dE", nchans );
      ttype[14] = "DAT_FREQ";
      tform[14] = Cstr16;
      tunit[14] = "";
      ttype[15] = "DAT_WTS ";
      tform[15] = Cstr16;
      tunit[15] = "";
      
      sprintf( Estr16, "%dE", nifs*nchans );
      ttype[16] = "DAT_OFFS";
      tform[16] = Estr16;
      tunit[16] = "";
      ttype[17] = "DAT_SCL ";
      tform[17] = Estr16;
      tunit[17] = "";
      
      sprintf( Istr16, "%dE", nifs*nchans*nbins );
      ttype[18] = "DATA    ";
      tform[18] = Istr16;
      tunit[18] = "Jy      ";
      
      nrows = 0; /* naxis2 - Let CFITSIO sort this out */
      ncols = 19; /* tfields */
      fits_create_tbl( fits, BINARY_TBL, nrows, ncols, 
		       ttype, tform, tunit, "SUBINT  ", &sta);
      
      /* Add dimensions of column 'ncols' = SUBINT Data */
      naxes[0] = nbins;
      naxes[1] = nchans;
      naxes[2] = nifs;
    
      fits_write_tdim( fits, ncols, 3, naxes, &sta );
      
      /* Add keywords */
      fits_update_key( fits, TSTRING, "INT_TYPE", "TIME",
		       "Time axis (TIME, BINPHSPERI, BINLNGASC, etc)", &sta);
      
      fits_update_key( fits, TSTRING, "INT_UNIT", "SEC",
		       "Unit of time axis (SEC, PHS (0-1), DEG)", &sta);
      
      fits_update_key( fits, TINT, "NCH_FILE", &nchans,
		       "Number of channels/sub-bands in this file (I)", &sta);
      j = 0;
      fits_update_key( fits, TINT, "NCH_STRT", &j,
		     "Start channel/sub-band number (0 to NCHAN-1) (I)", &sta);
    
      /* Store subint hdu number */
      fits_get_hdu_num( fits, &subint_hdu );
    }
    
    /* Write SUBINT BINTABLE columns */
    
    /* Fill in columns of table */
    col = 1;
    
    /* Subint number. If NAXIS=-1, 0 indicates EOD. */
    j=subint_cnt;
    fits_write_col( fits, TINT, col, subint_cnt, 1, 1, &j, &sta );
    col++;
    
    /* INDEXVAL - Optionally used if INT_TYPE != TIME */
    dx = 0.0;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;
    
    /* [s] Length of subint */
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dump_time, &sta );
    col++;
    
    /* [s] Offset from Start UTC of subint centre */
    dx = 0.0;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;
    
    /* [s] LST at subint centre */
    dx=0.0;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;
    
    /* [turns] RA (J2000) at subint centre */
    dx=src_raj/360.0;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;
    
    /* [turns] Dec (J2000) at subint centre */
    dx=src_dej/360.0;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;
    
    /* [deg] Gal longitude at subint centre */
    dx=srcl;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;

    /* [deg] Gal latitude at subint centre */
    dx=srcb;
    fits_write_col( fits, TDOUBLE, col, subint_cnt, 1, 1, &dx, &sta );
    col++;

    /* [deg] Feed angle at subint centre */
    x=0.0;
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
    col++;

    /* [deg] Parallactic angle at subint centre */

    /* [deg] Position angle of feed at subint centre */
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
    col++;

    /* [deg] Parallactic angle at subint centre */
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
    col++;

    /* [deg] Telescope azimuth at subint centre */
    x=(float) az_start;
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
    col++;

    /* [deg] Telescope zenith angle at subint centre */
    x=(float) za_start;
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, 1, &x, &sta );
    col++;

    /* Centre freq. for each channel - NCHAN floats */
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1,nchans,binned_freq,&sta );
    col++;

    /* Weights for each channel -  NCHAN floats */
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, nchans,
		    binned_weight, &sta );
    col++;

    /* Data offset for each channel - NCHAN*NPOL floats */
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, nchans*nifs,
		    binned_offset, &sta );
    col++;

    /* Data scale factor for each channel - NCHAN*NPOL floats */
    fits_write_col( fits, TFLOAT, col, subint_cnt, 1, nchans*nifs, 
		    binned_scale, &sta );
    col++;

    /* Subint data table - Dimensions of data table = (NBIN,NCHAN,NPOL) */
    for (i=0;i<nchans*nifs*nbins;i++) binned_data[i]=prof[i];

    fits_write_col(fits,TFLOAT,col, subint_cnt, 1, nbins*nchans*nifs, 
		   binned_data, &sta );

    subint_cnt++;
    if (sta) fits_report_error(stderr,sta);
#endif
  } else {
    /* EPN format requested - set up some general EPN variables */
    sprintf(epn.history,"%s %s fast-sampled data reduced using fold",
	    telescope_name(telescope_id),backend_name(machine_id));
    while (strlen(epn.history)<65) strcat(epn.history," ");
    strcpy(epn.jname,"");
    strcpy(epn.cname,"");
    epn.pbar=pfld0;
    epn.dm=refdm;
    epn.rm=0.0;
    strcpy(epn.catref,"none");
    strcpy(epn.bibref,"none");
    epn.raj=0.0;
    epn.dec=0.0;
    strcpy(epn.telname,telescope_name(telescope_id));
    epn.epoch=(float) floor(tstart);
    epn.opos=0.0;
    epn.paflag=' ';
    epn.timflag='U';
    epn.xtel=0.0;
    epn.ytel=0.0;
    epn.ztel=0.0;
    cal((double)epn.epoch,&epn.year,&epn.month,&epn.day);
    epn.scanno=0;
    epn.subscan=0;
    epn.npol=nifs;
    epn.nfreq=nchan;
    epn.nbins=nbins;
    epn.tbin=1.0e6*pfld0/(double)nbins;
    epn.nint=0;
    epn.ncal=0;
    epn.lcal=0;
    epn.tres=epn.tbin;
    epn.fluxflag='U';
    epn.navg=1;
    strcpy(epn.uf,"MHz ");
    epn.df=1000.0*fabs(foff);
    strcpy(epn.ud,"kHz ");
    if (epn.df>=10000.0) {
      epn.df/=1000.0;
      strcpy(epn.ud,"MHz ");
    }
    epn.tstart=(tstart-floor(tstart))*86400.0+tsta;
    epn.tstart*=1.0e6;
    epn.iprofile=(unsigned long *) malloc(epn.nbins*sizeof(long));
    profile = (float *) malloc(sizeof(float)*nbins);
    for (i=0;i<nifs;i++) {
      strcpy(epn.idfield,"I");
      for (c=0;c<nchan;c++) {
	for (b=0;b<nbins;b++) {
	  profile[b]=prof[i*nchan*nbins+c*nbins+b];
	}
	scale_prof(profile,nbins,epn.iprofile,&epn.scale,&epn.offset);
	epn.f0=fch1+foff*c;
	epn.nband=c+1;
	epn.papp=pfld0;
	epn.rms=0.0;
	epn.tres=tbin=1.0e6*pfld0*window/(double)nbins;
	if (c==0) write_epn_header(out,epn);
	write_epn_subheader(out,epn);
      }
    }
    fflush(out);
    free(profile);
    free(epn.iprofile);
  }
}
