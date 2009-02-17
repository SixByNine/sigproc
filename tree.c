/*
  TREE - dedisperses filterbank data using the tree algorithm. 
  
  Ramach, 09-MAY-2004, Green Bank.
  added to SIGPROC version 3.3 (drl-May 2005)
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>

int     DMMin, DMMax, FlipSwitch, FlipSwitch;
double  TSkip, TRead, DMMinv, DMMaxv;
char   *unfname;

#include "sigproc.h"
#include "header.h"
int nbands, nobits, userdm;
FILE *output;

void tree_help() {
  puts("");
  puts("tree - dedisperses filterbank data rapidly using the tree algorithm");
  puts("");
  puts("usage : tree <options> <UniqueID>");
  puts("");
  puts("-s skip - time length to skip (s; def=0)");
  puts("-r read - time length to read (s; def=largest 2^n*tsamp)");
  puts("-l dmlo - lower DM index to write (def=0)");
  puts("-u dmup - upper DM index to write (def=nchan-1)");
  puts("");
  exit(0);
}


void tree_parms(int argc, char *argv[]) {

  int   i;

  unfname    = (char *) calloc(120, sizeof(char));
  DMMin      = 0;
  DMMax      = 0;
  TSkip      = 0.0;
  TRead      = 0.0;
  FlipSwitch = 1;

  for (i=0; i<argc; i++) {
     if (strcmp(argv[i], "-l") == 0) {
      sscanf(argv[i+1], "%d", &DMMin);
    }

    if (strcmp(argv[i], "-u") == 0) {
      sscanf(argv[i+1], "%d", &DMMax);
    }

    if (strcmp(argv[i], "-s") == 0) {
      sscanf(argv[i+1], "%lf", &TSkip);
    }

    if (strcmp(argv[i], "-r") == 0) {
      sscanf(argv[i+1], "%lf", &TRead);
    }

    if (strcmp(argv[i], "-noflip") == 0) {
      FlipSwitch = 0;
      printf("\n");
      printf("THE BAND IS NOT FLIPPED!!\n");
    }
  }

  strcpy(unfname,argv[argc-1]);
  printf("\n");
  printf("Unique file name                  : %s\n",unfname);

  return;
}
/*  ======================================================================  */
/*  This function bit-reverses the given value "inval" with the number of   */
/*  bits, "nbits".    ----  R. Ramachandran, 10-Nov-97, nfra.               */
/*  ======================================================================  */


int bitrev(int inval,int nbits)
{
     int     ifact,k,i,ibitr;

     if(nbits <= 1)
     {
          ibitr = inval;
     }
     else
     {
          ifact = 1;
          for (i=1; i<(nbits); ++i)
               ifact  *= 2;
          k     = inval;
          ibitr = (1 & k) * ifact;

          for (i=2; i < (nbits+1); i++)
          {
               k     /= 2;
               ifact /= 2;
               ibitr += (1 & k) * ifact;
          }
     }
     return ibitr;
}

void AxisSwap(float Inbuf[],
	      float Outbuf[], 
	      int   nchans, 
	      int   NTSampInRead) {
  int    j1, j2, indx, jndx;

  for (j1=0; j1<NTSampInRead; j1++) {
    indx  = (j1 * nchans);
    for (j2=(nchans-1); j2>=0; j2--) {
      jndx = (j2 * NTSampInRead + j1);
      Outbuf[jndx]  = Inbuf[indx+j2];
    }
  }

  return;
}
void  FlipBand(float  Outbuf[], 
	       int    nchans, 
	       int    NTSampInRead) {

  int    indx, jndx, kndx, i, j;
  float *temp;

  temp  = (float *) calloc((NTSampInRead*nchans), sizeof(float));

  indx  = (nchans - 1);
  for (i=0; i<nchans; i++) {
    jndx = (indx - i) * NTSampInRead;
    kndx = (i * NTSampInRead);
    memcpy(&temp[jndx], &Outbuf[kndx], sizeof(float)*NTSampInRead);
  }
  memcpy(Outbuf, temp, (sizeof(float)*NTSampInRead * nchans));

  free(temp);

  return;
}
/*  ======================================================================  */
/*  This is a function to Taylor-dedisperse a data stream. It assumes that  */
/*  the arrangement of data stream is, all points in Chn.1, all points in   */
/*  Chn.2, and so forth.                                                    */
/*                     R. Ramachandran, 07-Nov-97, nfra.                    */
/*                                                                          */
/*  outbuf[]       : input array (short int), replaced by dedispersed data  */
/*                   at the output                                          */
/*  mlen           : dimension of outbuf[] (int)                            */
/*  nchn           : number of frequency channels (int)                     */
/*                                                                          */
/*  This programme has been debugged and finalised! The main programme used */
/*  for this purpose is "testtaylor.c", along with the function bitrev.c,   */
/*  which does the bit-reversal.  ----  R. Ramachandran, 10-Nov-97, nfra.   */
/*  ======================================================================  */

void taylor_flt(float outbuf[], int mlen, int nchn)
{
  float itemp;
  int   nsamp,npts,ndat1,nstages,nmem,nmem2,nsec1,nfin, i;
  int   istages,isec,ipair,ioff1,i1,i2,koff,ndelay,ndelay2;
  int   bitrev(int, int);

  /*  ======================================================================  */

  nsamp   = ((mlen/nchn) - (2*nchn));
  npts    = (nsamp + nchn);
  ndat1   = (nsamp + 2 * nchn);
  nstages = (int)(log((float)nchn) / 0.6931471 + 0.5);
  nmem    = 1;


  for (istages=0; istages<nstages; istages++) {
    nmem  *= 2;
    nsec1  = (nchn/nmem);
    nmem2  = (nmem - 2);

    for (isec=0; isec<nsec1; isec++) {
      ndelay = -1;
      koff   = (isec * nmem);

      for (ipair=0; ipair<(nmem2+1); ipair += 2) {
	ioff1   = (bitrev(ipair,istages+1)+koff) * ndat1;
	i2      = (bitrev(ipair+1,istages+1) + koff) * ndat1;
	ndelay++;
	ndelay2 = (ndelay + 1);
	nfin    = (npts + ioff1);
	for (i1=ioff1; i1<nfin; i1++) {
	  itemp      = (outbuf[i1] + outbuf[i2+ndelay]);
	  outbuf[i2] = (outbuf[i1] + outbuf[i2+ndelay2]);
	  outbuf[i1] = itemp;
	  i2++;
	}
      }
    }
  }

  return;
}


main(int argc, char** argv)
{
  int dum_int, iread, iw, indx, ITStart, namelen, HdrSize,
    NTOT, NSkip, NPtsToRead, SkiByte, log2NTOT, NDM,
    NOPFiles, MAX, NTSampInRead, i, NSampInRead, BytePerFrame,
    NByteInRead, ITOffset, NBitChan, *ibrev, NRead, IFiles, FOld,
    FNew, FSwitch, KeepTrack, NTReadOld, NTReadNew, NReadOld,
    NReadNew, OPFileSize, IOffset, SkipByte, oldper, newper ,NT_Files;

  long long TotalTrack;

  float  FMin, FMax, FCen, *Inbuf, *Outbuf, timecount;

  char  *filename, *filefull, *parfile;

  FILE  *fpar, *Fout[8192], *Fin;

  if (argc <= 1) 
    tree_help();
  else 
    tree_parms(argc, argv);

  filename   = (char *) calloc(120, sizeof(char));
  filefull   = (char *) calloc(120, sizeof(char));
  parfile    = (char *) calloc(120, sizeof(char));

  ITStart = 0; NT_Files = 0;

  /* To read a sample header from first time file */

  strcpy(filename, unfname);
  namelen    = strlen(filename);

  printf("First filename   : %s\n", filename);
  if ((fpar     = fopen(filename, "rb")) == NULL) {
    printf("ERROR opening file %s.\n", filename);
    exit(0);
  }
  HdrSize  = read_header(fpar);
  fclose(fpar);

  FMin     = fch1;
  FMax     = FMin + (float)(nchans - 1) * foff;
  FCen     = (FMin + FMax) / 2.0;

  printf("No. of frequency channels        : %d\n",nchans);
  printf("Beginning radio frequency        : %f MHz\n", FMin - (foff / 2.0));
  printf("Ending radio frequency           : %f MHz\n", FMax + (foff / 2.0));
  printf("Centre frq. of the whole band    : %f MHz\n", FCen);
  printf("\n");
  printf("input sampling interval       : %f sec\n", tsamp);

  NTOT = 0;
  strcpy(filename, unfname);
  namelen    = strlen(filename);
  if ((fpar       = fopen(filename, "rb")) == NULL) {
    printf("ERROR opening %s.\n", filename);
  }

  HdrSize    = read_header(fpar);
  fclose(fpar);

  NT_Files = (int)nsamples(filename, HdrSize, nbits, nifs, nchans);

  NTOT       += NT_Files;
  printf("File %5d   : %s with %d samples\n", i, filename, NT_Files);

  NSkip      = (int)(TSkip / tsamp);
  if (NSkip > NT_Files) {
    printf("Initial Skip-length longer than the first file!\n");
    exit(0);
  }

  NPtsToRead = (int)(TRead / tsamp);
  if (NPtsToRead==0) {
    TRead=NT_Files*tsamp;
    NPtsToRead=NT_Files;
  }
  if (NPtsToRead > NTOT) NPtsToRead = NTOT;
  SkipByte   = (NSkip * nbits * nchans / 8);

  printf("\n");
  printf("This data set contains %d number of samples\n", NTOT);
  printf("Number of bits per sample          : %d\n", nbits);
  printf("Time length to skip at the begin.  : %f\n", TSkip);
  printf("Bytes to skip at the beginning     : %d\n", SkipByte);
  printf("No. of time samples to skip        : %d\n", NSkip);
  printf("Time length to read after skipping : %f\n", TRead);
  printf("Samples to read after skipping     : %d\n", NPtsToRead);

  NTOT  -= NSkip;
  printf("Time samples after initial skip    : %d\n", NTOT);

  if (NPtsToRead > NTOT) {
    NPtsToRead  = NTOT;
    printf("Too many samples to read. Truncated to %d samples\n", NTOT);
  } else
    NTOT  = NPtsToRead;

  log2NTOT = (int)(log((double)NTOT) / log((double)2.0));
  NTOT     = (1 << log2NTOT);
  printf("After truncating to lower 2^n      : %d\n", NTOT);

  if (DMMin < 0) DMMin = 0;
  if (DMMax >= nchans) DMMax = (nchans - 1);
  if ((DMMin == 0) && (DMMax == 0)) {
    DMMin = 0; DMMax = (nchans - 1);
  } else if ((DMMin != 0) && (DMMax == 0)) 
    DMMax = (nchans - 1);
  printf("\n");
  printf("Minimum DM index = %d   Maximum = %d\n", DMMin, DMMax);
  DMMinv=(tsamp/8.3e3)*(double)(DMMin)*pow((FMax+FMin)/2.,3.)/fabs(FMax-FMin);
  DMMaxv=(tsamp/8.3e3)*(double)(DMMax)*pow((FMax+FMin)/2.,3.)/fabs(FMax-FMin);
  printf("Minimum DM value = %f   Maximum = %f\n", DMMinv,DMMaxv);
  NDM        = (DMMax - DMMin + 1);

  NOPFiles   = NDM;
  OPFileSize = sizeof(float) * NTOT;

  printf("\n");
  printf("Total number of output files    : %d\n", NOPFiles);
  printf("Size of each output file        : %d bytes\n", OPFileSize);

  for (i=0; i<NOPFiles; i++) {
    strcpy(filefull,unfname);
    sprintf(&filefull[strlen(filefull)-4], ".DM%.4d.tim", (i+DMMin));
    if((Fout[i] = fopen(filefull,"wb")) == NULL) {
      printf("ERROR opening output file %s.\n", filefull);
      exit(0);
    } else {
      nbands=1;
      nobits=32;
      refdm=(tsamp/8.3e3)*(double)(i+DMMin)*pow((FMax+FMin)/2.,3.)/
	fabs(FMax-FMin);
      output=Fout[i];
      dedisperse_header();
    }
  }
  printf("\n");

  /*
  strcpy(parfile, "\0"); strcpy(parfile, unfname);
  strncat(parfile, ".par", 4);
  fpar = fopen(parfile, "w");
  fprintf(fpar, "%s\n", unfname);
  fprintf(fpar, "%d  %d  %10.8f \n ", NTOT, nchans, tsamp);
  fprintf(fpar, "%d  %d  %d  %d\n", DMMin, DMMax, NDM, NOPFiles);
  fprintf(fpar, "%d\n", nchans);
  fprintf(fpar, "%f %f\n", FMin, FMax);
  fclose(fpar);
  */

  MAX          = (1 << 12);
  if (MAX > NTOT) MAX = NTOT;

  NTSampInRead = (MAX + (2 * nchans));
  NSampInRead  = (nchans * NTSampInRead);
  BytePerFrame = (nchans * nbits / 8);
  NByteInRead  = (NTSampInRead * BytePerFrame);

  ITOffset     = 2 * nchans;
  IOffset      = (-ITOffset * BytePerFrame);

  NBitChan = (int)(log((double)nchans) / 0.6931471);

  printf("Time samples to read in one read      : %d\n", NTSampInRead);
  printf("Bytes to read in one read             : %d bytes\n", NByteInRead);
  printf("Offset counter for each read          : %d time samples\n", IOffset);
  printf("\n");
  printf("No of bits to represent all channels  : %d\n", NBitChan);

  Inbuf    = (float *) calloc(NSampInRead, sizeof(float));
  Outbuf   = (float *) calloc(NSampInRead, sizeof(float));
  ibrev    = (int *)   calloc(nchans, sizeof(int));

  for (i=0; i<nchans; i++) {
    ibrev[i] = bitrev(i, NBitChan);
  }

  NRead    = (int)(NTOT / MAX);
  printf("Total number of reads in the run       : %d\n", NRead);

  IFiles   = 0;
  FOld     = 0;
  FNew     = 0;
  FSwitch  = 1;

  strcpy(filename, unfname);  namelen    = strlen(filename);
  printf("First filename   : %s\n", filename);
  if ((Fin = fopen(filename, "rb")) == NULL) {
    printf("ERROR opening file %s.\n", filename);
    exit(0);
  }
  HdrSize  = read_header(Fin);
  printf("Going into the main loop\n");

  fseek(Fin, SkipByte, SEEK_CUR);
  KeepTrack = (NT_Files - NSkip);
  TotalTrack = 0;

  newper   = 0.0;
  oldper   = newper;
  printf("\n\n");
  for (iread=0; iread<NRead; iread++) {

    newper = 100.0 * (((float)iread / (float)NRead));
    if (newper > oldper) {
      timecount = ((float)iread * NTSampInRead * tsamp);
      printf("\rProcessed : %3d%%    Current file : %s    Time from beg : %8.2f sec", 
	     (int)newper, filename, timecount); fflush(stdout);
      oldper = newper;
    }

    KeepTrack  -= NTSampInRead;
    TotalTrack += NTSampInRead;

    if (KeepTrack > 0) {
      read_block(Fin, nbits, Inbuf, NSampInRead);
      fseek(Fin, IOffset, SEEK_CUR);
      KeepTrack  += ITOffset;
      TotalTrack -= ITOffset;
    }

    if (KeepTrack > 0) FSwitch = 1;

    memset(&Outbuf[0], 0, (sizeof(float) * NSampInRead));

    AxisSwap(Inbuf, Outbuf, nchans, NTSampInRead);

    if (FlipSwitch == 0)  {
      printf("BEFORE..."); fflush(stdout);
      FlipBand(Outbuf, nchans, NTSampInRead);
      printf("AFTER!\n");
    }

    taylor_flt(Outbuf, NSampInRead, nchans);

    for (iw=DMMin; iw<=DMMax; iw++) {
      indx  = (ibrev[iw] * NTSampInRead);
      fwrite(&Outbuf[indx], sizeof(float), MAX, Fout[iw-DMMin]);
    }
  }

  newper    = 100.0 * ((float)iread / NRead);
  timecount = ((float)iread * NTSampInRead * tsamp);
  printf("\rProcessed : %3d%%    Current file : %s    Time from beg : %8.2f sec", 
	 (int)newper, filename, timecount); fflush(stdout);
  printf("\n\n");

  for (i=0; i<NOPFiles; i++) fclose(Fout[i]);

  exit(0);
}
