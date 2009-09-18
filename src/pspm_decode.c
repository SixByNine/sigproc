#include <stdio.h>
#include <math.h>

/*  
Based on decode.f.  Will unpack one 512-sample (= 32kb) block of PSPM
search data and fill the 1-d array "tmparray" with:
sample0 channel0, sample0 channel1,... sample0 channel127, sample1 channel0,...

It might make more sense to keep nboards, maxsam, nchpb and nfreqs in a
separate structure, but they're not all in PSPM_SEARCH_HEADER and therefore
need to be hardwired somewhere.

Ingrid Stairs Nov. 3, 1999 
*/

void pspm_decode(int *rawdata, float *tmparray) /* includefile */
{

  int i,j,k,l;
  unsigned int packed;
  int idx1,idx2;
  int nboards = 16;
  int maxsam = 512;
  int nchpb = 8;
  int nfreqs = 128;
  int nbits[8] = {16,20,24,28,0,4,8,12};
 
  for(j=0;j<maxsam;j++) 
    for(k=0;k<nboards;k++) {
      idx1 = j*nboards+k;
      idx2 = nfreqs-nchpb*k;
      packed = rawdata[idx1];
      for(l=0;l<nchpb;l++)
	tmparray[nfreqs*j+idx2-l-1] = (float) ((packed >> nbits[l]) & 15);
    }
}
