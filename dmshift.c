#include <stdio.h>
#include "sigproc.h"
int *dmshift(double f1, double df, int nchans, int nbands, double dm, double refrf, double tsamp, double frequency_table[]) /*includefile*/
{
  int i, cpb, *shift;
  double fi;

  shift = (int *) malloc(nchans * sizeof(int));
  fi=f1;
  cpb=nchans/nbands;
  if (frequency_table[0] != 0.0) f1=frequency_table[0];

  for (i=0; i<nchans; i++) {
    if (refrf > 0.0) f1=refrf;
    if (frequency_table[0] != 0.0) fi=frequency_table[i];
    shift[i]=(int)(dmdelay(fi,f1,dm)/tsamp);
    fi+=df;
    if (!((i+1)%cpb)) f1+=(double) cpb * df; 
  }
  return shift;
}
