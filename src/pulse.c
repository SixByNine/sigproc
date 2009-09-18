/*
  This program will search for non-periodic pulses in single-byte
  dedispersed data. It applues a threshold test on a time series for
  a single DM channel.

  This program serves as a driver program, reading in data, analyzing
  each DM channel, and looping over all channels.

  v. 0.0  TJWL  29 Sept 1994
  v. 0.1  MAM   26 Jul  1996 - corrected pulse indexing
  v. 0.2  MAM   15 Sep  2001 - adapted to Parkes data format
			     - compressed output files

  The inputs to this subroutine include. 

   'thresh' - Number of standard deviations by which a candidate must
	      exceed the mean to be considered a detection.

   'nsmax' - The maximum number of times to smooth the data.

   'ndm' - The DM index of the searched time series.
	 
   'length' - Number of time samples to be read.

   'data' - Input data array in single-byte format.

   'loopsize' - is equal to ndatfft/8, with its max defined in pulse.h 
*/
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "pulse.h"

#define STDOUT       1
#define FOPENMAX	100

int pulse_(realdata, thresh, nsmax, ndm, length, data, loopsize, scrdsk1, lsd1)
int *nsmax, *length, *ndm, *loopsize, *lsd1;
char scrdsk1[];
float *thresh, realdata[];
float *data;

{
  int done;
  int ndata;                          /* number of floats to be read       */
  int ndata2;                         /* number of floats sent to single_ch*/
  int nover;                          /* number of samples to overlap      */
                                      /* since the largest boxcar filter is*/
                                      /* 2**(NSMAX-1) "bins" wide, nover   */
                                      /* will be set to 2**(NSMAX+1)       */
  int n, i, iloop, iover, iopt, npulse, nloops, ierr;
  float testdata;
  double mean, rms, smean, srms, test;
  float *overlap;
  char *file, *file2, *outfile; 
  FILE *fd;

  scrdsk1[*lsd1]=NULL;

  nover = (int)pow(2.0, (double)(*nsmax+1));
  if ((overlap=(float *)malloc( sizeof(float)*(unsigned)nover)) == NULL )
    exit( 3 );

  if ( *thresh < 0.0 )
    {
      if ( fprintf(stderr, "Thresh = 0!\n")  == EOF ) exit( EOF );
      exit( 1 );
    }

  ierr = open_files(scrdsk1);		/* open new temporary output file */

  n = *loopsize - nover;                   /* used for           */
                                           /* overlapping        */
  iover = 0;
  iindx = 0;
  iloop = 0;
  nloops = (*length-nover)/(*loopsize) + 1;
  ndata = *loopsize;

  for (iloop=0; iloop<nloops; iloop++) {

  ndata2 = ndata + iover; 

  for (i=iover; i<ndata2; i++) realdata[i] = data[i+iindx];
  if ( iloop > 1 ) for (i=0; i<nover; i++)  realdata[i] = overlap[i];

  for (i=n; i<ndata2; i++) overlap[i-n] = realdata[i];

      single_ch(ndata2, realdata, *nsmax, *thresh, *ndm, scrdsk1); 

      iover = nover;
      ndata = n;
      iindx += n;
  }
return 0;
}
