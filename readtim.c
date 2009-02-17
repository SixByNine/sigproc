#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "header.h"
#include "sigproc.h"

FILE *timfile;

void openpresto_(char *filename){
  int i;
  /* first fix filename */
  for (i=0; i<strlen(filename); i++) {
    if (filename[i]==' ') {
      filename[i]='\0';
      break;
    }
  }
  timfile=open_file(filename,"rb");
  nbits=32;
}

void readhd_(char filename[80], float *rdm, float *rtsmp, char ra[12], char dec[12], float *fref, double *mjdstart, char srcname[80], char telname[80])
{
  int i;
  int rah,ram,ded,dem;
  double ras,des;
  char decsign,sra[6],sde[6];
  /* first fix filename */
  for (i=0; i<strlen(filename); i++) {
    if (filename[i]==' ') {
      filename[i]='\0';
      break;
    }
  }

  /* now open up the file and read in the header */
  refdm=tsamp=0.0;
  timfile=open_file(filename,"rb");
  read_header(timfile);
  strcpy(telname,telescope_name(telescope_id));
  *fref=(float)fch1;
  *mjdstart=tstart;
  *rtsmp=(float) tsamp;
  *rdm=(float) refdm;
  angle_split(src_raj,&rah,&ram,&ras);
  strcpy(srcname,source_name);
  if (ras<10.0) {
    sprintf(sra,"0%.1f",ras);
  } else {
    sprintf(sra,"%.1f",ras);
  }
  angle_split(src_dej,&ded,&dem,&des);
  if (src_dej > 0.0) 
    decsign = '+';
  else 
    decsign = '-';
  if (des<10.0) {
    sprintf(sde,"0%.1f",des);
  } else {
    sprintf(sde,"%.1f",des);
  }
  sprintf(ra,"%02d:%02d:%s",rah,ram,sra);
  sprintf(dec,"%c%02d:%02d:%s",decsign,abs(ded),dem,sde);
  if (nbits != 32 && nbits != 8 && nbits != 16) {
    fprintf(stderr,"cannot read %d bit data..\n");
    exit(1);
  }
}

int skipsample_(int *nskip)
{
  long offset;
  if (feof(timfile)) {
    fclose(timfile);
    return 0;
  }
  offset = (long) (*nskip) * (nbits/8);
  return fseek(timfile,offset,SEEK_CUR);
}
int readsample_(float *sample) 
{
  unsigned short *twobytes;
  char *byte;

  if (feof(timfile)) {
    fclose(timfile);
    return 0;
  }
  switch (nbits) {
  case 32:
    fread(sample,nbits/8,1,timfile);
    break;
  case 16:
    fread(twobytes,nbits/8,1,timfile);
    *sample = (float) *twobytes;
    break;
  case 8:
    fread(byte,nbits/8,1,timfile);
    *sample = (float) *byte;
    break;
  }
  return(1);
}
