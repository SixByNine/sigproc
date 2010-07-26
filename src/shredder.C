/*
  shredder  - shreds raw filterbank data into .tim files suitable for seek
  Each tim file represents an individual channel
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#define LONG64BIT unsigned long long
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include "gtools.h"
extern "C" {
#include "dedisperse_all.h"
};
#include "getDMtable.h"
//#include "gtools.h"


void inline_shredder_help(){
  fprintf(stderr,"shredder help\n");
  fprintf(stderr,"Usage: shredder filename [options]\n");
  fprintf(stderr,"-v                 verbose mode (prints status every DM trial)\n");
  fprintf(stderr,"-g gulpsize        number of samples to dedisp at once\n");
  fprintf(stderr,"-n Nsamptotal      Only do Nsamptotal samples\n");
  fprintf(stderr,"-b channel_begin   start channel (start at 0)\n");
  fprintf(stderr,"-e channel_end     end channel (end at nchan-1)\n");
  fprintf(stderr,"-s Nsamps          Skip Nsamp samples before starting\n");
  fprintf(stderr,"-s Nsamps          Skip Nsamp samples before starting\n");
  fprintf(stderr,"shredder uses OpenMP and 16 bit words to\n");
  fprintf(stderr,"create many dedispersed files at once in a highly\n");
  fprintf(stderr,"parallel manner. It is for use on 64 bit machines\n");
  fprintf(stderr,"but still works on 32 bit machines.\n");
  fprintf(stderr,"Currently tested on 96x1bit files, 512x1bit files, 1024x2bit files.\n");
#ifdef HAVE_OPENMP
  fprintf(stderr,"Compiled with OpenMP: Multi-threading enabled\n");
#else
  fprintf(stderr,"**SINGLE THREADED MODE**\nLink against OpenMP (-fopenmp with GNU on gcc > 4.2 for multi-threaded)\n");
#endif
}

FILE *input, *output, *outfileptr, *dmlogfileptr;
char  inpfile[180], outfile[180], ignfile[180], dmlogfilename[180];
char outfile_root[180];

/* ugly global variables describing the operating mode */
int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off;
double refrf,userdm,fcorrect;
float clipvalue,jyf1,jyf2;
int fftshift;
#include "wapp_header.h"
#include "key.h"
struct WAPP_HEADER *wapp;
struct WAPP_HEADER head;

/*
  tsamp in seconds, f0, df in MHz
  returns the DM for a given delay between the top and
  bottom frequency channel in samples.
*/

int main (int argc, char *argv[])
{
  /* local variables */
  char string[180];
  int i,nfiles=0,fileidx,sigproc,scan_number,subscan=0;
  int numsamps=0;
  unsigned char * rawdata;
  unsigned short int * unpacked; //, * times;
  size_t nbytesraw;
  int ibyte,j,k;
  unsigned char abyte;
  unsigned short int ** times;
  int nread;
  int ndm=0;
  float total_MBytes = 0;
  int counts;
  int verbose=0;
  int readsamp = 0;
  int nreadsamp = 0;
  int skip = 0;
  int nskip = 0;
  int ntoload;
  int ntodedisp;
  int maxdelay = 0;
  int appendable = 0;
  int ngulp; //max number of samples at a time
  int gulping = 0;
  int ngulpsize;
  int nsampleft; 
  int killing=0;
  int chan_begin=0;
  int chan_end=0;

  if(sizeof(LONG64BIT) != 8 ){
	  fprintf(stderr,"ERROR: sofware has been compiled with LONG64BIT as a datatype of %d bytes, needs to be 8\n",sizeof(LONG64BIT));
  }

  /* check number of command-line arguments and print help if necessary */
  if (argc<2) {
    inline_shredder_help();
    exit(0);
  }
  
  /* print help if necessary */
  if (strcmp(argv[1],"-h")==0) {
    inline_shredder_help();
    exit(0);
  }
  
  /* set up default globals */
  userbins=usrdm=asciipol=stream=clipping=swapout=headerless=0;
  sumifs=wapp_inv=wapp_off=barycentric=0;
  nobits=32;
  ascii=1;
  fftshift=1;
  profnum1 = 0;
  profnum2 = 1000;
  nbands=baseline=1;
  clipvalue=refrf=userdm=fcorrect=0.0;
  refdm=-1.0;
  output=NULL;
  strcpy(ignfile,"");

  // **************************************
  //          PARSE COMMAND LINE
  // **************************************
  i = 1;
  while (i<argc) {
    printf("argv[i] is %s\n",argv[i]);

    if (i==1) {
      if ((input=fopen(argv[i],"r"))!=NULL){
	printf("Opened file %s\n",argv[i]);
        strcpy(inpfile,argv[i]);
        nfiles++;
      }
      i++;
    }
    if (!strcmp(argv[i],"-v")) {
      verbose=1;
    }
    else if (!strcmp(argv[i],"-n")) {
      /* read only X samples */
      ntodedisp=atoi(argv[++i]);
      readsamp=1;
    }
    else if (!strcmp(argv[i],"-s")) {
      /* skip first X samples */
      nskip=atoi(argv[++i]);
      skip=1;
    }
    else if (!strcmp(argv[i],"-g")) {
      ngulpsize=atoi(argv[++i]);
      gulping=1;
    }
    else if (!strcmp(argv[i],"-b")) {
      chan_begin=atoi(argv[++i]);
      gulping=1;
    }
    else if (!strcmp(argv[i],"-e")) {
      chan_end=atoi(argv[++i]);
      gulping=1;
    }
    else {
      /* unknown argument passed down - stop! */
      inline_shredder_help();
      fprintf(stderr,"unknown argument (%s) passed to %s\n\n",argv[i],argv[0]);
      exit(1);
    }
    i++;
  }

  char tmp[180];
  char *tmp2;
  strcpy(tmp,inpfile); // These few lines strip the input file's
  tmp2 = basename(tmp);// name from it's directory tag.
  strcpy(outfile_root,tmp2);

  // **************************************
  //         VALUE PRECALCULATION
  // **************************************  
  /* read in the header to establish what the input data are... */
  sigproc=read_header(input);
  if (!sigproc) {
    fprintf(stderr,"Not sigproc data\n");
    exit(-1);
  }
  if (foff > 0.0) {
    fprintf(stderr,"dedisperse can't handle low->high frequency ordering!");
    exit(1);
  }

  if (fileidx == 1) {
    /* this is filterbank data */
    if (output!=stdout) output=fopen(outfile,"wb");
    if (output==NULL){
      perror("Outfile error 1\n");
      fprintf(stderr,"Error opening file %s\n",outfile);
      exit(-1);
    }
  }
  
  numsamps = nsamples(inpfile,sigproc,nbits,nifs,nchans);	/* get numsamps */
  if (chan_end==0) chan_end = nchans-1;

  // Sanity
  if (chan_end<chan_begin) {
    fprintf(stderr,"End channel %d less than start channel %d\n",
	    chan_end,chan_begin);
    exit(-1);
  }

  if (readsamp) numsamps=ntodedisp+maxdelay;
  
  if (gulping) {
    if (ngulpsize>numsamps) ngulpsize = numsamps-maxdelay;
    ntodedisp=ngulpsize;
  } else {
    ntodedisp=numsamps-maxdelay;
    ngulpsize=ntodedisp;
  }

  ntoload = ntodedisp + maxdelay; 
  
  nbytesraw = size_t (nchans) * size_t(ntoload * nbits/8);
  printf("Loading %U bytes of raw data\n",nbytesraw);
  if ((rawdata = (unsigned char *) malloc(nbytesraw))==(unsigned char *)NULL){
      fprintf(stderr,"Error allocating %U bytes of RAM for raw data\n",nbytesraw);
      exit(-1);
  }
  // skip either 0 or nskip samples
  fseek(input, size_t(nskip)*nchans*nbits/8, SEEK_CUR);
  nsampleft-=nskip;
  // some values used for unpacking
  int sampperbyte = (int)(8/nbits);
  int andvalue = (int)pow(2,nbits)-1;
  int nwholegulps = (numsamps - maxdelay)/ngulpsize;
  int nleft = numsamps - ngulpsize * nwholegulps;
  int ngulps = nwholegulps + 1;


  printf("ngulps %d\n",ngulps);

  // Start of main loop
  for (int igulp=0; igulp<ngulps;igulp++){
    if (igulp==nwholegulps) {
      ntoload = nleft;
      nbytesraw = size_t(ntoload)*nbits*nchans/8;
      ntodedisp = nleft - maxdelay;
    }

    //read gulp from file
    fprintf(stderr,"Gulp %d Loading %d samples, i. e. %u bytes of Raw Data\n",igulp, ntoload, nbytesraw);
    nread = fread(rawdata,nbytesraw,1,input);
    if (nread != 1){
      fprintf(stderr, "Failed to read %u nread = %d \n", nbytesraw, nread);
      exit(-1);
    }
  
     //shred file into channels 
      refdm=0;
      nobits=8;
      float orig_fch1 = fch1;
      int orig_nchans = nchans;
      for (int ichan=chan_begin;ichan<=chan_end;ichan++){
	sprintf(outfile,"%s.%4.4d.tim", outfile_root, ichan);
	if (igulp==0) {
	  outfileptr=fopen(outfile,"w");
	  if (outfileptr==NULL) {
            perror("Outfile error 2\n");
	    fprintf(stderr,"Error opening file %s\n",outfile);
	    exit(-1);
	  }
	  output=outfileptr;
	  nchans=1;
	  fch1=fch1+foff*ichan;
	  dedisperse_header();
	}
	fch1=orig_fch1;
	nchans=orig_nchans;
	if (igulp!=0) outfileptr=fopen(outfile,"a");
	ibyte=ichan/sampperbyte;
	int rotate=(ichan-ibyte*sampperbyte)*nbits; // Amount to rotate
	int stride = nchans/sampperbyte;
	unsigned char databytes[ntoload];
#pragma omp parallel for private (abyte,j)
	for (j=0;j<ntoload;j++){
	  abyte = rawdata[ibyte+j*stride];
	  databytes[j] = (unsigned char)((abyte>>rotate)&andvalue);
	}
	fwrite(databytes,ntoload,1,outfileptr);
	fclose(outfileptr);
      }
  }
  return(0); // exit normally
} // main
