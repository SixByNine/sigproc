/* 
  RCS: $Id: wapp_header.h,v 1.8 2004/02/17 13:13:25 wapp Exp wapp $
   wapp_head.h - header parameters for WAPP (Wideband Arecibo Pulsar Processor)
   This structure (WAPP_HEADER) is based on the PSPM-style headers but has
   been somewhat simplified to include just essential/relevant parameters. 
*/

#define HEADER_VERSION 8

/* lagformat for wapp_header */

#define INTLAGS   0   /* 16 bit integers - searching only      */
#define LONGLAGS  1   /* 32 bit integers - searching only      */
#define FLOATLAGS 2   /* 32 bit float ACF/CCFs folding only    */
#define FLOATSPEC 3   /* 32 bit float fftd ACFs folding only   */
#define BYTELAGS  4   /* 8 bit integers ACF/CCF searching only */
#define NIBBLAGS  5   /* 4 bit integers ACF/CCF searching only */
#define FOLD32    8   /* bit mask with FLOATLAGS and FLOATSPEC */
                      /* if on 32 bit folding is turned on     */

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
    int isfolding;       /* is folding selected                              */
    int isalfa;          /* is ALFA selected                                 */
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
    double poly_tmid[16];   /* mid point of polyco in (MJD) modified Julian date */     
    double coeff[192];     /* polynomial coefs made by TEMPO, 16 sets of 12 */
    int num_coeffs[16];     /* number of coefficients */
    char hostname[24];    /* ascii hostname of machine that took this data */
    int fold_bits;        /* 0 if 16 bits (old default) 1 if 32 bit folding */
    int iflo_flip;        /* consider entire iflo and determine flip */
    int attena;           /* first board parallel port value */
    int attenb;           /* second board parallel port value */

/* alfa info */

    double rfeed_offaz[7];  /* deg az rotated offset all alfa beam to center */
    double rfeed_offza[7];  /* deg za rotated offset all alfa beam to center */
    double prfeed_offaz;  /* deg az offset to center of alfa beam */
    double prfeed_offza;  /* deg za offset to center of alfa beam */
    double alfa_raj[7];   /* hr starting actual ra position of alfa beams */
    double alfa_decj[7];  /* deg starting actual dec position of alfa beams */
    double alfa_az[7];    /* deg starting actual az position of alfa beams */
    double alfa_za[7];    /* deg starting actual za postion of alfa beams */
    double alfa_ang;      /* deg alfa rotation angle  */
    double para_ang;      /* deg paralactic angle of center beam */
    char frontend[24];    /* receiver abbrev name */

/* add iflo data */

    double syn1;      /* upstairs synthesizer freq Hz */
    double synfrq[4]; /* downstairs synthesizers (Hz) */


    unsigned char prfeed;  /* centered alfa beam */
    unsigned char shcl;    /* true if receiver shutter closed */
    unsigned char sbshcl;   /* true if Sband receiver shutter closed */

    unsigned char rfnum;     /* position of the receiver selectror */
    unsigned char zmnormal;  /* transfer switch to reverse channels, true normal */
    unsigned char rfattn[2];  /* attenuator position */
    unsigned char ifnum;      /* if selector, 1/300 2/750, 3/1500, 4/10GHz1500, 5-thru */
    unsigned char ifattn[2];  /* IF attenuator positions */
    unsigned char fiber;      /* true fi fiber is chosen (always the case) */
    unsigned char ac2sw;      /* ac power to various instruments and other stuff */
    unsigned char if750nb;    /* narrow band 750 filter selected */

    unsigned char phbsig;     /* converter combiner signal phase adjust */
    unsigned char hybrid;     /* converter combiner hybrid  */
    unsigned char phblo;      /* convert combiner lo phase adjust */

/* downstairs */

    unsigned char xfnormal;   /* transfer switch true if normal */
    unsigned char noise;      /* noise on */
    unsigned char gain[2];    /* gain of downstairs amps */
    unsigned char inpfrq;     /* input distributor position */
    unsigned char mixer[4];   /* mixer source switches */
    unsigned char vlbainp;    /* vlba input switch position */
    unsigned char ampinp[4];  /* amplifier input source switches */
    unsigned char extinp[4];  /* external input selector switches */
    unsigned char syndest[4]; /* synthesizer destinations */
    unsigned char calsrc;     /* cal source bit */

    unsigned char vis30mhz;   /* greg 1 ch 0 */
    unsigned char pwrmet;     /* power meter input switch */
    unsigned char blank430;   /* 430 blanking on */
    unsigned char fill[6];    /* fill */
};
