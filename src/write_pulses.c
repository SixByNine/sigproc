/*
  This routine will write out the results of any pulses found.

*/
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "pulse.h"

int write_pulses(ns, npulse, pulse, ndm, scrdsk)
int ns;
int npulse;
char scrdsk[];
Pulsus pulse[];
{
  int i;
  int nchar;
  char *outfile;
  FILE *fd;
  char dir1[132];

  strcpy(dir1,scrdsk);
  strcat(dir1, "/best_tmp"); 

  /*
    Compose the file name based on the DM channel and amount of 
    smoothing.
    */
 
  if ( (outfile=(char *)malloc( 200)) == NULL ) 
    {
      if ( fprintf(stderr, "malloc failed in write_pulses at\n") == EOF ) 
	exit( EOF );
      if ( fprintf(stderr, "%d  %d  %d\n", ndm, ns, npulse) == EOF ) 
	exit( EOF );
      return( 1 );
    }
      
  sprintf(outfile, dir1);
 
  /*
    First, open the output data file for writing.
    */

  if ( (fd=fopen(outfile, "a")) == NULL ) 
    {
      if ( fprintf(stderr, "File open failed in write_pulses for\n") == EOF )
	exit( EOF );
      if ( fprintf(stderr, "%s\n", outfile) == EOF ) exit( EOF );
      return( EOF );
    }

  /*
    Write results.
    */

  for (i=0; i<npulse; i++)
    if ( fprintf(fd, 
		 "%d  %d  %d  %f  %f  %f\n", 
		 ndm, ns, pulse[i].index, 
                 pulse[i].amp, pulse[i].mean, pulse[i].rms) 
	== EOF ) 
      {
	if ( fprintf(stderr, "File write failed in write_pulses for\n") == EOF ) 
	  exit( EOF );
	return( EOF );
      }

  /*
    Close file.
    */
 
  fclose(fd);

}
