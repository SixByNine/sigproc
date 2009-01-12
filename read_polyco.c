#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sigproc.h"
#include "header.h"
int poly_override;
double override_f0;
int read_polycoset(FILE *polycofile, struct POLYCO *polyco) /*includefile*/
{
  char coefficient[25];
  int i,j;

  /* read in the first two lines of the polyco set (header info) */
  if((fscanf(polycofile,"%10s%9s%12lf%20lf%21lf%7lf%7lf",
      &polyco->psrname,&polyco->date,&polyco->utc,&polyco->tmid,&polyco->dm,
      &polyco->doppler,&polyco->lfrms))== -1) return(0);
  if ((fscanf(polycofile,"%20lf%18lf%5i%6i%5i%9lf\n",&polyco->rphase,
        &polyco->f0,&polyco->obsno,&polyco->span,&polyco->nc,
	&polyco->fobs)) == -1) return(0);
  if (poly_override) polyco->f0=override_f0;
  polyco->coeff = (double *) malloc((polyco->nc)*sizeof(double));

  i=0;
  while (i<polyco->nc) {
    /* look for coefficients of the form .....D... (TEMPO - f77 style) */
    if ((fscanf(polycofile,"%25s",coefficient))==-1) return(0);
    if ((strstr(coefficient,"D") != 0) || (strstr(coefficient,"E") != 0)) {
      /* replace the D/Es with an e and copy string to coefficient array */
      for (j=0; j<strlen(coefficient); j++) {
	if(coefficient[j]=='D') coefficient[j]='e';
	if(coefficient[j]=='E') coefficient[j]='e';
      }
      polyco->coeff[i]=atof(coefficient);
      i++;
    }
  }
  return(polyco->nc);
}

void get_nearest_polyco(char *filename, double mjd, struct POLYCO *polyco) /*includefile*/
{
  double diff,best,dt,tmid0,tmid1;
  int i,set,first=1,minutes;
  FILE *polycofile;
  char sname[80],pname[80];
  static double last=0.0;

  /* check to see when routine last called - return if still within a span */ 
  if (last==0.0) {
	last = mjd;
  } else {
        minutes = (mjd-last)*1440.0;
	/* 28/5/04 - bugfix - previous versions did not have factor of 2!!! */
	if (minutes >= polyco->span/2.0) {
		last = mjd;
	} else {
		return;
	}
  }

  polycofile=fopen(filename,"r");
  best=1.0e32;
  i=0;

  while (read_polycoset(polycofile,polyco) != 0) 
    {
      if (i==0) tmid0=polyco->tmid;
      if (i==1) tmid1=polyco->tmid;
      diff=fabs((polyco->tmid)-mjd);
      if (diff<best) {
	set=i;
	best=diff;
      }
      i++;
    }

  fclose(polycofile);

  /* check to see if the difference is within the spans of neighbouring sets */
  if (best > (tmid1-tmid0)) 
    fprintf(stderr,"WARNING: polyco spacing: %.0f s, nearest set: %.0f s\n",
	(tmid1-tmid0)*86400.0,best*86400.0);

  /* read in the best set */
  polycofile=fopen(filename,"r");
  for(i=0; i<=set; i++) read_polycoset(polycofile,polyco);
  fclose(polycofile);


  if (first) {
    /* strip off any B/Js from the source/polyco names */
    if (source_name[0] == 'J' || source_name[0] == 'B') 
      strcpy(sname,source_name+1);
    else
      strcpy(sname,source_name);

    if (polyco->psrname[0] == 'J' || polyco->psrname[0] == 'B') 
      strcpy(pname,polyco->psrname+1);
    else
      strcpy(pname,polyco->psrname);

    /*
    if (strings_equal(sname,"")) {
      fprintf(stderr,"WARNING: no source name in time series header!\n");
    } else if (!strings_equal(sname,pname)) {
      fprintf(stderr,
      "WARNING: source name (%s) and polyco name (%s) do not match!\n",
	      source_name,polyco->psrname);
    }
    */
    first=0;
  }

}

double polyco_period(double mjd, struct POLYCO polyco) /*includefile*/
{
  int i;
  double dt,freq,factor,period;

  dt=(mjd-polyco.tmid)*1440.0;
  factor=polyco.coeff[1];
  for (i=2; i<polyco.nc; i++) factor+=pow(dt,i-1)*i*polyco.coeff[i]; 
  freq=polyco.f0+(factor/60.0);
  period=1.0/freq;
  return (period);
}

double polyco_phase(double mjd, struct POLYCO polyco) /*includefile*/
{
  int i;
  double dt,phase;
  dt=(mjd-polyco.tmid)*1440.0;
  phase=polyco.rphase +(dt*60.0*polyco.f0)+polyco.coeff[0];
  for (i=1; i<polyco.nc; i++) phase+=pow(dt,i)*polyco.coeff[i]; 
  return(phase);
}
