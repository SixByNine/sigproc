struct POLYCO {
  char psrname[80];
  char date[80];
  double utc;
  double tmid;
  double dm;
  double doppler;
  double lfrms;
  double rphase;
  double f0;
  double fobs;
  int span;
  int obsno;
  int nc;
  double *coeff;
} ;
