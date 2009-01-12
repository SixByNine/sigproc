#include <stdio.h>
#include "sigproc.h"
int *ignored_channels(char *filename, int nchans) /*includefile*/
{
  int i,idx,*ignore;
  FILE *ignfile;
  
  /* allocate space for ignore array and initialize */
  ignore = (int *) malloc(nchans * sizeof(int));
  for (i=0;i<nchans;i++) ignore[i]=0;

  /* read list of ignored channel numbers from file */
  ignfile=open_file(filename,"r");
  while (1) {
    fscanf(ignfile,"%d",&idx);
    if (feof(ignfile)) break;
    idx--;
    if ( (idx >= 0) && (idx < nchans) ) ignore[idx]=1;
  }
  close(ignfile);

  return(ignore);
}
