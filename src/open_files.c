#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "pulse.h"

int open_files(char scrdsk[])
{

  char *outfile;
  FILE *fd;
  char dir1[132];

  strcpy(dir1,scrdsk);

  strcat(dir1,"/best_tmp"); 
 
  outfile=(char *)malloc(200);
  sprintf(outfile, dir1);
 
  /*
    Open output files for writing.
    */

  if ( (fd=fopen(outfile, "w")) == NULL ) 
    {
      if ( fprintf(stderr, "File open failed in write_pulses for\n") == EOF )
	exit( EOF );
      if ( fprintf(stderr, "%s\n", outfile) == EOF ) exit( EOF );
      return( EOF );
    }

  fclose(fd);

  if ( (fd=fopen(outfile, "w")) == NULL )
    {
      if ( fprintf(stderr, "File open failed in write_pulses for\n") == EOF )
        exit( EOF );
      if ( fprintf(stderr, "%s\n", outfile) == EOF ) exit( EOF );
      return( EOF );
    }

   fclose(fd);
}
