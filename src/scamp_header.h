/* 
 scamp_header.h - C structure of header files for SCAMP filterbank data
 
 dlorimer@atnf.csiro.au - Oct 20, 2001 - adapted from sc_pmhead.inc
*/
typedef struct {

  char tape[6];       /* Tape label */
  char file[4];       /* File number on tape */
  char block[8];      /* Block counter in file */
  char date[8];       /* UT date 'yy/mm/dd' */
  char mjd[8];        /* MJD at 00h UT */
  char ut[16];        /* UT at file start ' hh:mm:ss.ssssss' */
  char lst[12];       /* LST at file start ' hh:mm:ss.ss' */
  char coord[2];      /* Coord system  04 = Galactic  05 = J2000 */
  
  char ra[16];        /* RA at file start ' hh:mm:ss.ssss' */
  char dec[16];       /* Dec at file start '-dd:mm:ss.sss' */
  char l[8];          /* Gal longitude at start 'ddd.dddd' */
  char b[8];          /* Gal latitude at file start  '+dd.dddd' */
  char fangle[8];     /* Feed angle (degrees) */
  char obstime[8];    /* Obs length (seconds of time) 'ssss.sss' */
  char comment[64];   /* Comment entered by user */
  char nfilter[2];    /* No. of filter systems (01 or 02) */
  char chbw[2][8];    /* Chan inc (MHz, -ve if inverted) 'bb.bbbb' */
  char nchan[2][4];   /* No. of channels in each filter system */
  char fch1[2][12];   /* RFof 1st chan center    'ffffff.fffff' */
  
  char tsamp[2][12];  /* Sample interval in us */
  char nsgrp[2][4];   /* Samples per group (=1 for pmdaq) */
  char ngblk[8];      /* Groups per block */
  char blocksec[8];   /* Seconds per tape block */
  char fdcntrl[2];    /* 0 = none, 1 = fixed FA, 2 = fixed PA or GPA */
  char dtype[2][1];   /* Data type code (2=pol, 3=norm, 5=dedisp) */
  
  char uthd[16];      /* UT of block start ' hh:mm:ss.ssssss' */
  char lsthd[12];     /* LST of block start ' hh:mm:ss.ss' */
  char rahd[16];      /* RA at block start ' hh:mm:ss.ssss' */
  char dechd[16];     /* Dec at block start '-dd:mm:ss.sss' */
  char lhd[8];        /* Gal long at block start 'ddd.dddd' */
  char bhd[8];        /* Gal lat at block start '+dd.dddd' */
  char zahd[8];       /* Zenith angle at block start 'ddd.dddd' */
  char azhd[8];       /* Azimuth at block start 'ddd.dddd' */
    
  char atten[4][4];   /* Attenuator settings (db) */
  char tpower[20][4]; /* Total powers (Jy) 'iiii' */
  char nblk_read[8];  /* Number of tape blocks in disk file */
  char scanl[8];      /* scan rate in deg/min in long '-r.rrrrr' */
  char scanb[8];      /* scan rate in deg/min in lat  '-r.rrrrr' */
  char nbeam[4];      /* total number of beams */
  char ibeam[4];      /* number of this beam */
  char psrname[16];   /* pulsar name including B/J */
  char obsfile[16];   /* config file name used */
  char nbits[2][2];   /* number of bits per sample */
  char dmdd[8];       /* DM for dedispersion */
  char nddch[2][4];   /* number of channels per dedispersed band */
  char move[2];       /* 'k0': On k(0-9), '01': grd '02': off '03:scanning'*/
  
  char pnterr[6];     /* pointing error in acmin mmm.mmm */
  char tree[2];       /* ' ' normal 'T' tre dedisp 'D' pdm dedisp */
  char nsys[2];       /* Filter systems disk file */
  char telid[10];     /* Telescope ID */
  char pangle[8];     /* Parallactic angle degrees */
  char nspbsw[8];     /* Number of samples per beam switch */
  char pcalcyc[4];    /* Cal cycle period in samples */
  
  char reserve[22];

} SCAMP_HEADER;
