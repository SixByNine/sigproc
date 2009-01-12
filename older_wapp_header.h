
/* 
  RCS: $Id: wapp_header.h,v 1.5 2001/10/17 12:38:13 wapp Exp wapp $
   wapp_head.h - header parameters for WAPP (Wideband Arecibo Pulsar Processor)
   This structure (WAPP_HEADER) is based on the PSPM-style headers but has
   been somewhat simplified to include just essential/relevant parameters. 
*/

#define HEADER_VERSION 4

/* lagformat for wapp_header */

#define INTLAGS   0   /* 16 bit integers - searching only    */
#define LONGLAGS  1   /* 32 bit integers - searching only    */
#define FLOATLAGS 2   /* 32 bit float ACF/CCFs folding only  */
#define FLOATSPEC 3   /* 32 bit float fftd ACFs folding only */

struct WAPP_HEADER {

    long header_version; /* some integer that increments with each revision */
    long header_size;    /* size (in bytes) of this header (nom =1024) */
    char obs_type[24];   /* what kind of observation is this */
                         /* PULSAR_SEARCH */
                         /* PULSAR_FOLDING */
                         /* SPECTRA_TOTALPOWER */
/* 
    The following are obtained from current telescope status display
    note that start AST/LST are for reference purposes only and should 
    not be taken as accurate time stamps. The time stamp can be derived
    from the obs_date/start_time variables further down in the structure.
*/
    double src_ra;       /* requested ra J2000 (10000*hr+100*min+sec) */
    double src_dec;      /* requested dec J2000 (10000*deg+100*min+sec) */
    double start_az;     /* telescope azimuth at start of scan (deg) */
    double start_za;     /* telescope zenith angle at start of scan (deg) */
    double start_ast;    /* AST at start of scan (sec) */
    double start_lst;    /* local siderial time at start of scan (sec) */
/*  
    In the following, anything not supplied/requested by the user
    is assumed to be calculated by WAPP when it writes the header
*/
    double cent_freq;    /* user-supplied band center frequency (MHz) */
    double obs_time;     /* user-requested length of this integration (s) */
    double samp_time;    /* user-requested sample time (us) */
    double wapp_time;    /* actual sample time (us) i.e. requested+dead time */
    double bandwidth;    /* total bandwidth (MHz) for this observation */

    long num_lags;       /* user-requested number of lags per dump per spect */
    long scan_number;    /* built by WAPP from year+daynumber+3-digit-number */

    char src_name[24];   /* user-supplied source name (usually pulsar name) */
    char obs_date[24];   /* built by WAPP from yyyymmdd */
    char start_time[24]; /* UT seconds after midnight (start on 1-sec tick) */
    char project_id[24]; /* user-supplied AO proposal number (XYYYY) */
    char observers[24];  /* observer(s) name(s) */

    int nifs;            /* user-requested: number of IFs to be recorded     */
    int level;           /* user-requested: 1 means 3-level; 2 mean 9-level  */
    int sum;             /* user-requested: 1 means that data is sum of IFs  */
    int freqinversion;   /* 1 band is inverted, else band is not inverted    */
    long long timeoff;   /* number of reads between obs start and snap block */
    int lagformat;       /* 0=16 bit uint lags , 1=32 bit uint lags          */
                         /* 2=32 bit float lags, 3=32 bit float spectra      */
    int lagtrunc;        /* if we truncate data (0 no trunc)                 */
                         /* for 16 bit lagmux modes, selects which 16 bits   */
                         /* of the 32 are included as data                   */
                         /* 0 is bits 15-0 1,16-1 2,17-2...7,22-7            */
    int firstchannel;    /* 0 when correlator channel a is first, 1 if b     */
    int nbins;           /* number of time bins for pulsar folding mode      */
                         /* doulbles as maxrecs for snap mode                */
    double dumptime;     /* folded integrations for this period of time      */
    double power_analog[2];   /* Power measured by Analog Detector           */
/*    
    In the following, pulsar-specific information is recorded for use 
    by folding programs e.g. the quick-look software. This is passed to 
    WAPP by psrcontrol at the start of the observation. 

    The apparent pulse phase and frequency at time "dt" minutes with
    respect to the start of the observation are then calculated as:

    phase = rphase + dt*60*f0 + coeff[0] + dt*coeff[1] + dt*dt*coeff[2] + ...
    freq(Hz) = f0 + (1/60)*(coeff[1] + 2*dt*coeff[2] + 3*dt*dt*coeff[3] + ...)

    where the C notation has been used (i.e. coeff[0] is first coefficient etc)
    for details, see TEMPO notes (http://www.naic.edu/~pulsar/docs/tempo.txt)
*/
    double psr_dm;         /* pulsar's dispersion measure (cm-3 pc) */
    double rphase[16];      /* reference phase of pulse (0-1) */
    double psr_f0[16];      /* pulse frequency at reference epoch (Hz) */
    double poly_tmid[16];   /* mid point of polyco in (MJD) modified Julian date
 */     
    double coeff[192];     /* polynomial coefs made by TEMPO, 16 sets of 12 */
    int num_coeffs[16];     /* number of coefficients */

     /* this pads out the header to 2048 bytes */
};
