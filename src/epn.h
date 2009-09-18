/*
  C structure containing EPN header information and data for a pulse profile
  for further details on the header parameters, see the infamous EPN paper:
  Lorimer et al. in A&A Suppl. Series, v.128, p.541-544 [astro-ph/9801097] 
*/
struct EPN {
  char history[72];     /* some (preferably meaningful) notes about the data */
  char jname[12];       /* Name of pulsar derived from its J2000 coordinates */
  char cname[12];       /* Common name of pulsar (either B- or J-name) */
  double pbar;          /* Barycentric period of pulsar (s) (present epoch) */
  double dm;            /* Dispersion measure (pc/cc) */
  double rm;            /* Rotation measure of pulsar (rad m^2) */
  char catref[6];       /* Tag showing which catalogue in use */
  char bibref[8];       /* Bibliographic reference for data */
  double raj;           /* Right ascension J2000 hhmmss.s */
  double dec;           /* Declination J2000 ddmmss.s */
  char telname[8];      /* Telescope name */
  float epoch;          /* modified Julian date of observation */
  float opos;           /* position angle of feed (degrees) */
  char paflag;          /* "A" signifies absolute PAs else undefined */
  char timflag;         /* "U" signifies UTC time "B" barycentric UTC */
  float xtel;           /* Topocentric rectangular position of telescope [m] */
  float ytel;           /* Topocentric rectangular position of telescope [m] */
  float ztel;           /* Topocentric rectangular position of telescope [m] */
  int day;              /* year (eg 1996) */
  int month;            /* month (1-12) */
  int year;             /* day (1-31) */
  int scanno;           /* sequence number of the observation */
  int subscan;          /* sub-sequence number of the observation */
  int npol;             /* number of polarisations observed */
  int nfreq;            /* number of frequency bands per polarisation */
  int nbins;            /* number of phase bins per frequency */
  double tbin;          /* sampling interval (us) */
  int nint;             /* number of integrated pulses per data block */
  int ncal;             /* bin number of start of cal signal */
  int lcal;             /* number of bins in the cal signal */
  double tres;          /* temporal resolution (us) */
  char fluxflag;        /* "F" signifies flux (mJy) calibrated data */
  char idfield[8];      /* Description of data stream I Q U V etc. etc. */
  int nband;            /* Ordinal number of data stream */
  int navg;             /* Number of streams averaged into this one */
  double f0;            /* Effective centre sky frequency of this stream */
  char uf[8];           /* String giving unit of f0 [default is GHz] */
  double df;            /* Effective bandwidth of this stream */
  char ud[8];           /* String giving unit of df [default is MHz] */
  double tstart;        /* Time of first bin wrt EPOCH [us] */
  float scale;          /* scale factor for the data */
  float offset;         /* offset to be added to the data */
  float rms;            /* rms of the data */
  double papp;          /* apparent period at time of first phase bin (us) */
  /* 
     the final entry in the structure is iprofile ---- the profile scaled to
     an integer value between 0 and 65535. To convert this to floating point
     use the scaling:  profile[i]=scale*(float)iprofile[i]+offset
  */
  unsigned long *iprofile; 
} ;
