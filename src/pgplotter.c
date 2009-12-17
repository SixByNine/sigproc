#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/* 
   pgplotter - a simple tool to plot SIGPROC streamed output with PGPLOT
   now included as part of sigproc-3.0 onwards (14-Jan-2005 - drl@jb.man.ac.uk)

   usage: fold ..... -stream | pgplotter
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "sigproc.h"

int verbose=0;
int maxpts=4096;
int get_stream(float *x, float *y) 
{
  float a,b;
  int i,stream;
  char line[80], key[80];

  strcpy(line,"");
  i=stream=0;
  while (1) {
    fgets(line,80,stdin);
    sscanf(line,"%s",key);
    if (!strcmp(key,"#DONE")) return(0);
    if (!strcmp(key,"#STOP")) return(i);
    if (stream) {
      if (i >= maxpts) {
	maxpts+=4096;
	x = (float *) realloc(x,maxpts * sizeof(float));
	y = (float *) realloc(y,maxpts * sizeof(float));
      }
      sscanf(line,"%f %f",&x[i],&y[i]);
      i++;
    }
    if (!strcmp(key,"#START")) stream=1;
  }
}
void minmax(float *data, int npts, float *datamin, float *datamax)
{ 
  int i; 
  *datamin=data[0];
  *datamax=data[0];
  for (i=1; i<npts; i++) {
    *datamin = (data[i] < *datamin) ? data[i] : *datamin;
    *datamax = (data[i] > *datamax) ? data[i] : *datamax;
  }
}
main(int argc, char *argv[])
{
  float *x, *y, xmin, xmax, ymin, ymax;
  int i,npts,nx,ny;
  char device[80];

  x = (float *) malloc(maxpts * sizeof(float));
  y = (float *) malloc(maxpts * sizeof(float));
  xmin=xmax=ymin=ymax=0.0;
  nx=ny=1;

  if (argc>2) print_version(argv[0],argv[1]);

  strcpy(device,"/xs");

  if (argc>1) nx=atoi(argv[1]);
  if (argc>2) ny=atoi(argv[2]);
  if (argc>3) strcpy(device,argv[3]);

  cpgbeg(0,device,nx,ny);
  cpgsvp(0.15,0.85,0.15,0.85);
  cpgsch(1.5);
  cpgscf(2);
  while (npts=get_stream(x,y)) {
    if (verbose) printf("read %d points... (max %d)\n",npts,maxpts);
    minmax(x,npts,&xmin,&xmax);
    minmax(y,npts,&ymin,&ymax);
    if (verbose) printf("xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);
    cpgask(0);
    cpgpage();
    cpgswin(xmin,xmax,ymin,ymax);	
    cpgbox("bcnst",0.0,0,"bcnst",0.0,0);
    cpgline(npts,x,y);
  }
  cpgend();
}
