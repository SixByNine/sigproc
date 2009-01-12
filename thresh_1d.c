/*
  This routine will apply the threshold test.

  */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "pulse.h"

void thresh_1d(ndata, realdata, ns, npulse, pulses, thresh, scrdsk)
int ndata;
char scrdsk[];
float realdata[], thresh;
int ns;
int *npulse;
Pulsus *pulses;
{
  int i, nd;
  double mean, rms, smean, srms, test;
  char *outfile;

  char dir1[132];

  strcpy(dir1,scrdsk);

  mean = realdata[0];
  rms = realdata[0]*realdata[0];

  for (i=1; i<ndata; i++)
    {
      mean += realdata[i];
      rms += realdata[i]*realdata[i];
    }

  smean = mean;
  mean /= (double)ndata;
  srms = rms;
  test = (double)rms - (double)ndata*pow((double)mean, 2.0);
  if (test < 0) printf("test = %f, rms = %f, mean = %f, ndata = %d\n",test,rms,mean,ndata);
  rms = (double)sqrt(((double)rms - (double)ndata*pow((double)mean, 2.0))/(double)(ndata-1));

  /*
    Run through the time series again, this time comparing each sample
    with the threshold.

    NOTE:  The important thing, when figuring out the index in the time
    series, is not how many times the time series has been smoothed, but
    how many times it's been *decimated*.  The way single_ch() is currently
    structured, the number of times the series has been decimated is (ns-1),
    except for ns==0.

    Also, adjust the mean and rms to remove pulses.  This prepares for
    the final pass through the series below.
    */

  *npulse = 0;

  if (ns)
    nd = (int)pow(2.0, (double)(ns-1));
  else
    nd = 1;


  for (i=0; i<ndata; i++)
    {
      float x;

      x = (realdata[i] - mean)/rms; 
      if ( x > thresh )
        {
          pulses[*npulse].index = i*nd + iindx;
          pulses[*npulse].amp = x;
          pulses[*npulse].mean = mean;
          pulses[*npulse].rms = rms;

          smean -= realdata[i];
          srms -= realdata[i]*realdata[i];

          (*npulse)++;
        }
    }
      mean = smean/(double)(ndata-*npulse);
      test = (double)srms - (double)(ndata-*npulse)*pow((double)mean, 2.0);
      if (test < 0) printf("test = %f, rms = %f, mean = %f, ndata = %d\n",test,rms,mean,ndata);
      rms = (double)sqrt(((double)srms - (double)(ndata-*npulse)*pow((double)mean, 2.0))/(double)(ndata-*npulse-1));

      *npulse = 0;
      for (i=0; i<ndata; i++)
        {
          float x;

          x = (realdata[i] - mean)/rms;

          if ( x > thresh )
            {
              pulses[*npulse].index = i*nd + iindx;
              pulses[*npulse].amp = x;
              pulses[*npulse].mean = mean;
              pulses[*npulse].rms = rms;

              smean -= realdata[i];
              srms -= realdata[i]*realdata[i];

              (*npulse)++;
            }
    }

  return;
}


