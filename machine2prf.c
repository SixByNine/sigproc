/*
## help on machine2prf
## machine2prf - read and dedisperse pre-folded data from various machines
## 
## passed down:
## 
## FILE *input  - pointer to raw data
## FILE *output - pointer to output device/file
## 
## dedispersed profiles are converted into either ASCII or EPN format. This
## routine is called from the dedisperse program. Type "dedisperse help" for
## a synopsis of the various command-line options. Currently supported machines
## are: PSPM, WAPP and BPP. Reading of headers is done in "typeof_inputdata.c"
## 
## see also:
## 
## dedisperse.c typof_intputdata.c 
## 
## last modified: Mar 27, 2001 (dunc@naic.edu)
## help end
*/

#include "dedisperse.h"
#include "pspmhdr.h"
PSPM_TIMING_HEADER pspm_timing;
#include "wapp_header.h"
struct WAPP_HEADER *wapp;
struct WAPP_HEADER head;
#include "bpphdr.h"
BPP_TIMING_HEADER bpp_timing;
double fmid;
double period;
int nbins;
void machine2prf(FILE *input, FILE *output) /* includefile */
{
  int i,b,c,dayno,day=0,month=0,year=0,subscan=0,scanno=0,nprds,nchans0;
  int *chtab,*ignore;
  unsigned long *iprofile;
  char  *telname,*catref,*bibref,*paflag,*timflag,*srcname;
  double dm,rm,*chanfreq,pra,pde,df;
  double phase,psrfreq,phase0;
  float *profile,opos,offset=0.0,scale=1.0,rms=0.0;
  struct EPN epn;

  if (data_type == 3) {
    /* SIGPROC binary profile format */
    wapp_off=0;
    srcname=source_name;
    fmid=fch1+fabs(foff)/2.0+foff*((float)nchans)/2.0;
    foff=fabs(foff);
    profile = (float *) malloc(sizeof(float)*nifs*nchans*nbins);
    fread(profile,sizeof(float),nifs*nchans*nbins,input);
    mjdobs=floor(tstart);
    tstart=(tstart-mjdobs)*86400.0;
    phase=0.0;
    nprds=npuls;
    dm=0.0;
  } else {
    /* set up header variables and read in data depending on the machine */
    switch (machine_id) {
    case 1:
      /* PSPM data */
      wapp_off=0;
      fmid=pspm_timing.freq; 
      foff=0.062;
      df=7.68;
      nbins=pspm_timing.num_phase_bins;
      nchans=pspm_timing.num_chans;
      period=pspm_timing.psr_period;
      read_aoscan(pspm_timing.scan_num, &dayno, &year, &subscan);
      dm = pspm_timing.psr_dm;
      chtab=pspm_chans(nchans);
      profile=pspm_prof(input,nbins,nchans,chtab);
      nifs=1;
      nprds=pspm_timing.num_periods;
      srcname=pspm_timing.psr_name+3;
      pra=pspm_timing.psmon_ra;
      pde=pspm_timing.psmon_dec;
      break;
    case 2:
      /* WAPP data */
      fmid=wapp->cent_freq;
      if(fcorrect>0.) fmid = fcorrect;
      foff=-(1-2*wapp_inv)*wapp->bandwidth/wapp->num_lags;
      df=wapp->bandwidth;
      nbins=wapp->nbins;
      nchans=wapp->num_lags;
      nifs=wapp->nifs;
      period=1.0/wapp->psr_f0[0];
      dm = wapp->psr_dm;
      nprds=1;
      pra=wapp->src_ra;
      pde=wapp->src_dec;
/*      if(profnum2 > 1000) profnum2 = wapp_ndumps; */
      profile=wapp_prof(nbins,nchans,nifs,profnum1,profnum2);
      srcname=wapp->src_name;
/*      tstart += profnum1*wapp->dumptime; */
      if((fabs(period-0.040))>0.001) { /* not a cal scan */
	phcalc(mjdobs,tstart,&phase0,&psrfreq,wapp->rphase,wapp->psr_f0,
	       wapp->poly_tmid,wapp->coeff,wapp->num_coeffs);
	tstart += wapp->dumptime/2.;
      /* Actually midpoint, not start, of combined dumps */
        tstart += (profnum1+profnum2)/2.*wapp->dumptime;
        phcalc(mjdobs,tstart,&phase,&psrfreq,wapp->rphase,wapp->psr_f0,
               wapp->poly_tmid,wapp->coeff,wapp->num_coeffs); 
        period = 1./psrfreq;
/*	printf("Phases 1: %f %f\n",phase, phase0); */
        phase -= phase0;  /* first sample went into bin 0 */
/*	printf("Phases 2: %f %f\n",phase, phase0);  */
        while (phase < 0.) phase += 1.;
      }
      break;
    case 4:
      /* BPP data */
      nbins=bpp_timing.num_phase_bins;
      nchans=bpp_timing.num_chans;
      period=bpp_timing.psr_period;
      read_aoscan(bpp_timing.scan_num, &dayno, &year, &subscan);
      dm = bpp_timing.psr_dm; 
      chtab=bpp_chans(bpp_timing.bandwidth,bpp_timing.mb_start_address,
		      bpp_timing.mb_end_address,bpp_timing.mb_start_board,
		      bpp_timing.mb_end_board,bpp_timing.cb_id,
		      bpp_timing.aib_los,bpp_timing.dfb_sram_freqs,
		      bpp_timing.rf_lo);
      df=foff*nchans;
      profile=pspm_prof(input,nbins,nchans,chtab);
      nifs=1;
      nprds=bpp_timing.num_periods;
      srcname=bpp_timing.psr_name+3;
      pra=pde=0.0;
      break;
    default:
      error_message("unknown machine!");
      break;
    }
  }

  /* user can specify a DM to dedisperse at (default is to value in header) */
  if (usrdm) dm=userdm;

  /* user can specify number of bins in output profiles (decimation) */
  if (userbins==0) userbins=nbins;

  /* set up number of dedisperse subbands (default is 1) */
  if (nbands==0) nbands=nchans;

  /* save starting number of channels */
  nchans0=nchans;

  /* calculate sky frequency of each channel */
  /*  chanfreq=chan_freqs(fmid,fabs(foff),nchans); */
  chanfreq=chan_freqs(fmid,foff,nchans,wapp_off);

  /* zero any profiles that are in the ignored list of channels */
  if (file_exists(ignfile)) {
    ignore=ignored_channels(ignfile,nchans);
  } else {
    ignore=(int *) malloc(nchans*sizeof(int));
    for (i=0;i<nchans;i++) ignore[i]=0;
  }

  /* check to see if reference frequency has been flagged at the command line*/
  if (refrf == -1.0) refrf=fmid; /* use mid frequency if this is the case */

  /* dedisperse profiles if requested and scale by Jyfactors if two-pol data */
  if (dm>0.0) 
    prof_ddis(profile,nbins,nchans,nbands,nifs,chanfreq,period,dm,refrf,jyf1,jyf2);

  /* sum profiles to form subbands if requested */
  if (nbands<nchans) prof_sumc(profile,nbins,nbands,&nchans,nifs,ignore);

  /* decimate profiles if requested */
  if (nbins>userbins) prof_adds(profile,&nbins,nchans,nifs,nbins/userbins);
  
  /* sum IFs if requested -- this alters nifs to 1 */
  if (sumifs) prof_sumifs(profile,nbins,nchans,&nifs);

  /* subtract baseline (median of each profile) if requested */
  if (baseline) {
    prof_sbas(srcname,profile,nbins,nchans,nifs); 
  } else {
    /* no baseline subtraction but normalise by number of periods/channels
       to get units in terms of the digitization range (e.g. 0-15 PSPM) */
    //    for (i=0;i<nbins*nchans*nifs;i++) profile[i]/=(float)nchans0*nprds/nbands;
  }

  /* set up sampling time (us) for output profile (nbins may be different) */
  tsamp=1.0e6*period/nbins;

  /* multiply outgoing profiles by single Jansky calibration factor if supplied */
  if (jyfactor != 1.0) for (i=0;i<nbins*nifs*nchans;i++) profile[i]*=jyfactor;

  if (ascii && !stream) {
    /* ascii output format requested */
    for (i=0;i<nifs;i++) {
      for (c=0;c<nchans;c++) {
	if (!headerless) 
        fprintf(output,"# %.1f %.7f %.10f %d %.3f %.3f %d %c %d %s %.8f\n",
        (float)mjdobs,tstart,period,nprds,chanfreq[c],dm,nbins,tempo_site(telescope_id),subscan,
        srcname,phase);
	for (b=0;b<nbins;b++) 
	  fprintf(output,"%d %f\n",b,profile[c*nbins*nifs+i*nbins+b]);
      }
    }
  } else if (asciipol) {
    /* write profiles in format for Jim's polarization code */
    for (b=0;b<nbins;b++) 
      for (i=0;i<nifs;i++)
	for (c=0;c<nchans;c++) 
	 fprintf(output,"%d %d %d %f\n",b,i,c,profile[c*nbins*nifs+i*nbins+b]);
  } else if (stream) {
    for (i=0;i<nifs;i++) {
      for (c=nchans-1;c>=0;c--) {
	fprintf(output,"#START %d %f %f\n",nbins,tstart,chanfreq[c]);
	for (b=0;b<nbins;b++) {
	 fprintf(output,"%d %f\n",b,profile[c*nbins*nifs+i*nbins+b]);
	}
	fprintf(output,"#STOP\n");
      }
    }
    fprintf(output,"#DONE\n");
  } else {
    /* EPN format requested - set up some general EPN variables */
    sprintf(epn.history,"%s %s timing-mode data reduced using dedisperse",
	    telescope_name(telescope_id),backend_name(machine_id));
    while (strlen(epn.history)<65) strcat(epn.history," ");
    strcpy(epn.jname,srcname);
    strcpy(epn.cname,srcname);
    epn.pbar=period;
    epn.dm=dm;
    epn.rm=0.0;
    strcpy(epn.catref,"none");
    strcpy(epn.bibref,"none");
    epn.raj=pra;
    epn.dec=pde;
    strcpy(epn.telname,telescope_name(telescope_id));
    epn.epoch=(float)mjdobs;
    epn.opos=0.0;
    epn.paflag=' ';
    epn.timflag='U';
    epn.xtel=0.0;
    epn.ytel=0.0;
    epn.ztel=0.0;
    epn.day=0;
    epn.month=0;
    epn.year=year;
    epn.scanno=scanno;
    epn.subscan=subscan;
    epn.npol=nifs;
    epn.nfreq=nchans;
    epn.nbins=nbins;
    epn.tbin=1.0e6*period/(double)nbins;
    epn.nint=0;
    epn.ncal=0;
    epn.lcal=0;
    epn.tres=epn.tbin;
    epn.fluxflag='U';
    epn.navg=nchans0/nbands;
    strcpy(epn.uf,"MHz ");
    epn.df=df;
    df/=(double)nbands;
    epn.df=df*1000.0;
    strcpy(epn.ud,"kHz ");
    if (epn.df>=10000.0) {
      epn.df/=1000.0;
      strcpy(epn.ud,"MHz ");
    }
    epn.tstart=tstart;
    epn.iprofile=(unsigned long *) malloc(epn.nbins*sizeof(long));
    /* loop over IFs and channels writing an EPN record foreach profile */
    for (i=0;i<nifs;i++) {
      strcpy(epn.idfield,"I");
      for (c=0;c<nchans;c++) {
	epn.f0=fch1-foff*c;
	epn.nband=c+1;
	epn.papp=period;
	scale_prof(profile,nbins,epn.iprofile,&epn.scale,&epn.offset); 
	epn.rms=0.0;
	write_epn(output,epn);
	profile+=nbins;
      }
    }
  }
}
