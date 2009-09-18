#include <stdio.h>
#include "sigproc.h"
double *bandfactors(int nchans) /*includefile*/
{
  FILE *file;
  double *corrfactor,dtmp;
  int j;
  float temp;

  corrfactor = (double *) malloc(nchans * sizeof(double));
  /* read in bandpass correction factors from file "file" */
  file=open_file("bandfactors","r");
  for (j=0; j<nchans; j++) {
    fscanf(file,"%f %lf",&temp,&dtmp);
    corrfactor[j]=dtmp;
  }
  close(file);

  return(corrfactor);
}
