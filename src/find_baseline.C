
/*
  A fortran->C++ translation of find_baseline.f
  M Bailes 13 Jan 2009

  dat is an array of the data to be baselined (floats)
  ndat is the dimension of the array (int)
  tsamp_dat is the sample time in seconds
  tsmooth is the smoothing time in seconds
  threshold is the number of sigma to avoid for running mean averaging

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//actually this returns the standard deviation and reports the mean
float find_rms(float * dat, int ndat, float * ave){
  double s=0.0;
  double ss=0.0;
  for (int j=0;j<ndat;j++){
    s=s+dat[j];
    ss=ss+dat[j]*dat[j];
  }

  double div = 1.0 / (double) ndat;
  *ave=s*div;
  return((float)sqrt(ss*div-((s*div)*(s*div))));
}

void find_baseline(float * dat, int ndat, double tsamp_dat, double tsmooth,
		   float threshold){

  int verbose = 1;
  if (tsamp_dat == 0.0){
    fprintf(stderr,"find_baseline error: tsamp_dat is zero\n");
    exit(-1);
  }

  int nrm = (int) (tsmooth/tsamp_dat);
  if (nrm>ndat){
    fprintf(stderr,"find_baseline error: running mean longer than data\n");
    exit(-2);
  }
//  fprintf(stderr,"nrm is %d\nndat is %d\n",nrm,ndat);

  //  if (verbose) cout << "ndat " << ndat << "n running mean " << nrm << endl;

  float * tdat = new float[nrm];
  if (tdat == NULL) {
    fprintf(stderr,
	    "find_baseline error: Error allocating %d floats for tdat\n",nrm);
    exit(-3);
  }
  int nrh = nrm/2;  // Half of running mean number

  // Copy nrm dat's into tdat
  for (int i=0;i<nrm;i++){
      tdat[i]=dat[i];
  }

  // RMS of the whole thing
  float ave_f;
  float rms_f = find_rms(dat,ndat,&ave_f);
  fprintf(stderr,"found good rms: %f\n",rms_f);
  if (rms_f==0.0){
    fprintf(stderr,"find_baseline error: rms of data is zero\n");
    exit(-4);
  }
  float inverse_rms_f = 1.0/rms_f;

  // Subtract mean to half of smoothing time
  for (int j=0;j<nrh;j++) dat[j]=(dat[j]-ave_f)*inverse_rms_f;
  //------dorunning mean (eliminating high points)
//  fprintf(stderr,"down here\n");
   
  int na=nrh;
  int nb=ndat-nrh;
  int ja=0;
  float ave = ave_f;
  float div = 1.0 / nrm;
  float threshold_r = threshold * rms_f;

  for (int j=na;j<nb;j++){
    dat[j]=dat[j]-ave;            // Subtract average
    if (ja == nrm-1) ja = 0;
    float al=tdat[ja];
    float an=dat[j+na];
    if (fabs(an-ave) > threshold_r) an=ave; 
    tdat[ja]=an;
    ave=ave+(an-al)*div;
    dat[j]=dat[j]*inverse_rms_f;
    ja++;
  }
  // do end bit
  for (int j=nb;j<ndat;j++)dat[j]=(dat[j]-ave)*inverse_rms_f;
}

