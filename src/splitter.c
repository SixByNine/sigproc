#include <stdio.h>
#include <stdlib.h>
#include "sigproc.h"
#define DEFBLOCK 32768
#define DEFBYTES 2000000000
main(int argc, char **argv)
{
  FILE *output;
  char stem[80], newfile[80];
  unsigned char block[DEFBLOCK];
  int blocksize, filesize, nfiles=1, nread, snip=0, snipped=0;
  if (argc < 2) {
    puts("usage: splitter outputstem (filesize/bytes)");
    exit(0);
  }
  strcpy(stem,argv[1]);
  sprintf(newfile,"%s.%d",stem,nfiles);
  output=fopen(newfile,"wb");
  fprintf(stderr,"opened output file... %s\n",newfile);
  if (strings_equal(argv[1],"snip")) snip=1;

  if (argc>2) 
    filesize=atoi(argv[2]);
  else 
    filesize=DEFBYTES;

  if (DEFBLOCK > filesize) 
    blocksize=filesize;
  else
    blocksize=DEFBLOCK;

  while (nread=fread(block,1,blocksize,stdin)) {
    fwrite(block,1,nread,output);
    fflush(output);
    if (sizeof_file(newfile)>=filesize) {
      if (snip && snipped) {
	continue;
      } else {
	fclose(output);
	sprintf(newfile,"%s.%d",stem,++nfiles);
	output=fopen(newfile,"wb");
	fprintf(stderr,"opening new file... %s\n",newfile);
	snipped=1;
      }
    }
  }
  fclose(output);
}
