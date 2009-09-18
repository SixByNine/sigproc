#include <stdio.h>
FILE *open_file(char *filename, char *descriptor) /* includefile */
{
  FILE *fopen(), *fptr;
  if ((fptr=fopen(filename,descriptor)) == NULL) {
    fprintf(stderr,"Error in opening file: %s\n",filename);
    exit(1);
  }
  return fptr;
}
