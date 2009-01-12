#include "fold.h"
/* folds incoming blocks of data */
/* modified 28 May 2001 to fix phase counting [JMC]
   Former:
		turn = tsec/psec;
		phase = turn - floor(turn);
   This doesn't work for binaries with fastly changing periods!
   Instead of using integrated time, need to use differential
   time / period.
   New:

		turn += tsmp / psec;
		phase = turn - floor(turn);	
*/
double tsec, tsmp, tsta, tlas, psec, turn, phase, window, pfld, pfld0, tadd;
double dump_last; /* added by MAM Nov 2004 - to write out dumps at exact
		     time values, rather than to the nearest pulse */
long int pulsecount;

/* set override flags to be OFF! - used only for depolyco.c */
int poly_override=0;
double override_f0=0.0;

float  *fold_data() /* includefile */
{
  float *data, *prof, *temp, *base, *cntp, *cntt, *samps;
  char string[80];
  FILE *pfile;
  long int tempcount, bytes_to_skip;
  int ibin,lbin,i,s,c,n,nread,nsamps,newbaseline,pupd,prfidx,datidx,blocksize;
  struct POLYCO polyco;
  double kappa, dump_next, speed_of_light=299792458.0;
  /* 
     initialise folding variables (all double precision):
     tsec - elapsed time in seconds
     tsta - elapsed time at start of first bin 
     turn - number of pulse turns
     phase - phase of sample
     window - phase range selected by user (default=1.0)
     tsmp - sample time in seconds (including WAPP dead time)
     psec - folding period in seconds
     pfld - mean folding period over the whole scan in second
  */
  lbin=-1;
  newbaseline=1;
  pulsecount=tempcount=dump_last=0;
  tsta=tsec=tadd=turn=0.0;
  tsmp=tsamp;			
  window=phase_finish-phase_start;
  kappa=0.0;
  dump_next=dump_time; /* keep track of the time to dump next */

  /* add arbitrary time offset to start time */
  tstart+=time_offset/86400.0;

  /* skip forward in the time series if need be */
  if (skip_time > 0.0) {
    tstart+=skip_time/86400.0;
    bytes_to_skip=(long)(skip_time/tsamp)*(long)(nbits*nifs*nchans/8);
    /*fprintf(stderr,"%ld %lf %d %d %d\n",
      bytes_to_skip,skip_time,nbits,nifs,nchans);*/
    if (fseek(input,bytes_to_skip,SEEK_CUR)) 
      error_message("skip error - check time value of -sk option");
  }
  if (read_time==0.0) read_time=1.0e32;
  
  /* allow the user to override sample time in the header from command line */
  if (tsamp_user > 0.0) tsmp=tsamp_user*1.0e-6;

  /* folding: acceleration; a polyco file; a .top file or from command line */
  if (acceleration != 0.0) {
    /* normalization for midpoint as reference */
    kappa=folding_period/1000.0/(1.0+acceleration*tobs/2.0/speed_of_light);
    pfld=psec=kappa;
    pupd=1;
  } else if (folding_period == -1.0) {
    pfile=open_file(polyco_file,"r");
    if (!read_polycoset(pfile,&polyco)) {
      /* it's a .top file - read period directly */
      rewind(pfile);
      fscanf(pfile,"%lf",&folding_period);
      pfld=psec=folding_period/1000.0;
      pupd=0;
    } else {
      /* it's a polyco - get nearest set and prodceed */
      get_nearest_polyco(polyco_file,tstart,&polyco); 
      pfld=psec=polyco_period(tstart,polyco);
      pupd=1;
    }
    fclose(pfile);
  } else {  
    /* use period from command line */
    pfld=psec=folding_period/1000.0;
    pupd=0;
  }
  pfld0=psec;
  /* if user has not specified bins make them phase_window/tsmp */
  if (nbins<1) nbins=(int) (window*psec/tsmp);

  /* current EPN format has maximum bin restriction (1k) */
  if (!ascii && !stream && nbins>1024) nbins=1024;

  /* work out nearest number of pulses for time dumps */ 
 /* if (dump_time > 0.0) npulses = (int) rint(dump_time/psec); */ 

  /* allocate space for the data and profiles */
  blocksize=nifs*512*nchans;
  data=(float *) malloc(blocksize*sizeof(float));
  prof=(float *) malloc(nifs*nbins*nchans*sizeof(float));
  cntp=(float *) malloc(nbins*sizeof(float));
  temp=(float *) malloc(nifs*nbins*nchans*sizeof(float));
  cntt=(float *) malloc(nbins*sizeof(float));
  base=(float *) malloc(nifs*nchans*sizeof(float));

  /* initialize baseline array */
  for (i=0; i<nifs*nchans; i++) base[i]=0.0;

  /* initialize profile arrays */
  for (i=0; i<nbins*nchans*nifs; i++) prof[i]=temp[i]=0.0;

  /* initialize bin counts */
  for (i=0; i<nbins; i++) cntt[i]=cntp[i]=0.0;

  while ((nread=read_block(input,nbits,data,blocksize))>0
	 && (tsec<read_time) ) {

    /* number of samples to process in this block */
    nsamps=nread/nchans/nifs;

    if (newbaseline && baseline) {
      samps=(float *) malloc(nsamps*sizeof(float));
      for (i=0; i<nsamps; i++) samps[i]=0.0;
      /* calculate a base line to subtract from this subintegration */
      for (i=0; i<nifs; i++) {
	n=i*nchans;
	for (c=0; c<nchans; c++) {
	  for (s=0; s<nsamps; s++) samps[s]=data[c+n+s*nchans*nifs];
	  base[n+c]=nrselect(nsamps/2,nsamps,samps-1);
	}
      }
      newbaseline=0;
      free(samps);
    }
    /* update the period if necessary and send out log */
    if (pupd) {
      if (acceleration != 0.0) {
	pfld=psec=kappa*(1.0+acceleration*tsec/speed_of_light);
      } else {
	get_nearest_polyco(polyco_file,tstart+tsec/86400.0,&polyco); 
	pfld=psec=polyco_period(tstart+tsec/86400.0,polyco);
      }
    }
    sprintf(string,"time:%.1fsP(fold):%.10fs",tsec, psec);
    update_log(string);

    /* main folding loop */
    for (s=0; s<nsamps; s++) {

      /* calculate phase of each sample */
      /* turn=tsec/psec; */
      phase=turn-floor(turn);

      /* process this sample only if within selected window */
      if ( (phase >= phase_start) && (phase <= phase_finish) ) {
	
	ibin=nbins*((phase-phase_start)/window);

	if ( (ibin-lbin) < 0 ) pulsecount++;

	cntp[ibin]+=1.0;
	cntt[ibin]+=1.0;
	
	/* loop over IFs adding sample into appropriate phase bin */
	for (i=0; i<nifs; i++) {
	  n=i*nchans;
	  /* loop over channels */
	  for (c=0; c<nchans; c++) {
	    prfidx=n*nbins+c*nbins+ibin;
	    datidx=c+n+s*nchans*nifs;
	    prof[prfidx]+=data[datidx]-base[n+c];
	    temp[prfidx]+=data[datidx]-base[n+c];
	  }
	}
      }

      /* if user requested profile-dump mode - see if this is a new pulse */
      if (npulses > 0) {
	if ( (ibin-lbin) < 0 ) {
	  if ((pulsecount) && !(pulsecount%npulses)) {
	    /* normalize prof by counts */
	    norm_prof(temp,cntt,nbins,nifs,nchans);
	    /* 
	       dump out the profile and prepare period for start of next dump
	       added use of temporary variable to keep track of the pulses
	       integrated within this sub-integration (drl - jul 2, 2004)
	    */
	    pulsecount-=tempcount;
	    write_profiles(temp,nbins,nchans,nifs,output);
	    pulsecount+=tempcount;
	    tempcount=pulsecount;
            pfld0=psec;
	    if (accumulate) {
	      denorm_prof(temp,cntt,nbins,nifs,nchans);
	    } else {
	      /* initialize subint profile and subint counts */
	      for (i=0; i<nbins*nchans*nifs; i++) temp[i]=0.0;
	      for (i=0; i<nbins; i++) cntt[i]=0.0;
	    }
	    /* set flag to get new baseline for next subintegration */
	    newbaseline=1;
	    /* OLD WAY (wrong!) 
		update start time based on number of pulses added 
	    tsta+=psec*(double)npulses; */
	    /* NEW WAY (right!)
		update start time based on number of samples since start */
	    tsta+=tsec-dump_last;
	    dump_last = tsec;
	  }
	}
      }
  
      /* when dumping profiles every dump_time seconds, seek the
         closest time sample -- this was previously in error!! */
      if (dump_time > 0) {
	if (dump_next-tsec <= tsmp/2.0) {
	   norm_prof(temp,cntt,nbins,nifs,nchans);
	   write_profiles(temp,nbins,nchans,nifs,output);
	   pfld0=psec; 
	   if (accumulate) {
	     denorm_prof(temp,cntt,nbins,nifs,nchans);
	   } else {
	     /* initialize subint profile and subint counts */
	     for (i=0; i<nbins*nchans*nifs; i++) temp[i]=0.0;
	     for (i=0; i<nbins; i++) cntt[i]=0.0;
	   }
	   /* set flag to get new baseline for next subintegration */
	   newbaseline=1; 
	   /* update start time based on number of samples since start */
	   tsta+=tsec-dump_last;
	   dump_last = tsec;
	   dump_next+= dump_time;
        }
      } 

      lbin=ibin; /* keep track of last bin number */

      /* update elapsed time counter */
      tsec+=tsmp;
      turn += tsmp/psec;
    }
  }

  /* normalize prof by counts */
  norm_prof(prof,cntp,nbins,nifs,nchans);

  /* free arrays and return folded profile(s) */
  free(data); free(base); free(temp);
  tsta=0.0;
  return prof;
}
