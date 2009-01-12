#include <math.h>
long long sizeof_file(char name[]) ;
long long nsamples(char *filename,int headersize, int nbits, int nifs, int nchans) /*includefile*/
{
  long long datasize,numsamps;
  datasize=sizeof_file(filename)-headersize;
  numsamps=(long long) (long double) (datasize)/ (((long double) nbits) / 8.0) 
		 /(long double) nifs/(long double) nchans;
  return(numsamps);
}
