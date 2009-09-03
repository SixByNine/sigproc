#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

void normalise(int,float*);

float * Creadspec(char * filename, int * npf, double * rate)
{
  float * f;
  FILE * fptr = fopen(filename,"r");
  if (fptr==NULL){
    fprintf(stderr,"Error opening file %s\n",filename);
    exit(-1);
  }

  float dm,ac;
  fseek(fptr,8,SEEK_CUR);    // crazy fortran io
  fread(&dm,sizeof(float),1,fptr);
  fread(&ac,sizeof(float),1,fptr);
  fseek(fptr,8,SEEK_CUR);    // crazy fortran io

  fprintf(stderr,"\ndm %f\n ac %f\n",dm,ac);

  double tsamp; int fold;
  fseek(fptr,8,SEEK_CUR);    // crazy fortran io
  fread(&tsamp,sizeof(double),1,fptr);
  printf("tsamp is %lf\n",tsamp);
  *rate = 1.0 / (tsamp);
  fread(npf,sizeof(float),1,fptr);
  fread(&fold,sizeof(float),1,fptr);
  fseek(fptr,8,SEEK_CUR);    // crazy fortran io

  fprintf(stderr,"\tsamp %lf\n npf %d\nfold %d\n",tsamp,*npf,fold);

  f = (float *) malloc(sizeof(float)* *npf);
  if (f == NULL){
    fprintf(stderr,"Error allocating %d bytes in readspec\n",sizeof(float)* *npf);
    exit(-2);
  }

  fprintf(stderr,"Reading in %d floats from %s\n",*npf,filename);
  fseek(fptr,4,SEEK_CUR);    // crazy fortran io
  int nread = fread(f,sizeof(float),*npf,fptr);
  printf("Read %d floats from %s\n",nread,filename);
  normalise(nread,f);

  fseek(fptr,4,SEEK_CUR);    // crazy fortran io

  fclose(fptr);
  return(f);
}

