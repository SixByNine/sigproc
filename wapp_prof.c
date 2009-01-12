#include <stdio.h>
#include <fcntl.h>
#include "sigproc.h"
int wapp_file;

/* read in profiles from wapp timing-mode data returns a pointer to them */
float *wapp_prof(int nbins, int nchans, int nifs, int np1, int np2) /*includefile*/
{
  int swap_bytes,i,ip,b,c,n,rec_size;
  float *profile, *dump; 
  double pavg;
 
  swap_bytes=big_endian();
  rec_size=nchans*nifs*nbins*sizeof(float);
  profile = (float *) malloc(rec_size); 
  dump = (float *) malloc(rec_size);

  for(n=0;n<nchans*nifs*nbins;n++) profile[n]=0.;

  for(ip=0;ip<np1;ip++) 
    if ((read(wapp_file,dump,rec_size)) != rec_size) 
      error_message("reading WAPP profile...");

  pavg = (double)(np2-np1+1.);
  for(ip=np1;ip<=np2;ip++) {
    if ((read(wapp_file,dump,rec_size)) == rec_size) {
      n=0;
      for (c=0;c<nchans;c++) {
	for (i=0;i<nifs;i++) {
	  for (b=0;b<nbins;b++) {
	    if (swap_bytes) swap_float(&dump[b*nifs*nchans+i*nchans+c]);
	    profile[n]+=dump[b*nifs*nchans+i*nchans+c]/pavg;
	    /*	printf("%d %d %d %d %f\n",c,i,b,n,profile[n]); */
	    n++;
	  }
	}
      }
    }
    else {
      if(ip == np1)
	error_message("reading WAPP profile...");
      else 
	break;
    }
  }


  free(dump);
  return(profile);
}
