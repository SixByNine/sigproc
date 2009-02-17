/* 
   pgplotter - a simple tool to plot SIGPROC streamed output with PGPLOT
   now included as part of sigproc-3.0 onwards (14-Jan-2005 - drl@jb.man.ac.uk)

   usage: fold ..... -stream | pgplotter

   Modification history

   Apr 25, 2007 - Duncan.Lorimer@mail.wvu.edu - added -twod flag for gray scale
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "sigproc.h"
int noffset=0;
int flip=0;
int verbose=0;
int maxpts=32768*32;
float xxmin=0.0,yymin=1.0e32,zzmin=1.0e32,grmin=1.0e32;
float xxmax=0.0,yymax=-1.0e32,zzmax=-1.0e32,grmax=-1.0e32;
float pdm=0.0, offset[999];
int get_stream(float *x, float *y) 
{
  float xval,yval;
  int i,n,stream;
  char line[80], key[80];

  strcpy(line,"");
  i=stream=0;
  while (1) {
    fgets(line,80,stdin);
    sscanf(line,"%s",key);
    if (!strcmp(key,"#DONE")) return (0);
    if (!strcmp(key,"#STOP")) return (i);
    if (stream) {
      if (i >= maxpts) {
	maxpts+=4096;
	x = (float *) realloc(x,maxpts * sizeof(float));
	y = (float *) realloc(y,maxpts * sizeof(float));
      }
      sscanf(line,"%f %f",&x[i],&y[i]);
      i++;
    }
    if (!strcmp(key,"#START")) {
      stream=1;
      sscanf(line,"%s %d %f %f",key,&n,&xval,&yval);
      yval*=1000.0;
      if (xval>zzmax) zzmax=xval;
      if (xval<zzmin) zzmin=xval;
      if (yval>yymax) yymax=yval;
      if (yval<yymin) yymin=yval;
    }
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
  float *x, *y, xmin, xmax, ymin, ymax, zmin, zmax, *twod, tr[6], deltat,f;
  int i,npts,nlines,nx,ny,pggray,ntwod,first,itwod,subints;
  char device[80], title[80];

  x = (float *) malloc(maxpts * sizeof(float));
  y = (float *) malloc(maxpts * sizeof(float));
  xmin=xmax=ymin=ymax=zmin=zmax=0.0;
  nx=ny=first=1;
  nlines=pggray=itwod=0;
  tr[0]=0.0;
  tr[1]=1.0;
  tr[2]=0.0;
  tr[3]=1.0;
  tr[4]=0.0;
  tr[5]=1.0;
  strcpy(device,"/xs");
  strcpy(title,"");

  if (argc > 1) {
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"-twod")) {
	pggray=1;
      } else if (strings_equal(argv[i],"-flip")) {
	flip=1;
      } else if (strings_equal(argv[i],"-pgdev")) {
	strcpy(device,argv[++i]);
      } else if (strings_equal(argv[i],"-title")) {
	strcpy(title,argv[++i]);
      } else if (strings_equal(argv[i],"-zmin")) {
	zmin=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-zmax")) {
	zmax=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-xmin")) {
	xxmin=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-xmax")) {
	xxmax=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-dm")) {
	pdm=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-offset")) {
	offset[noffset++]=atof(argv[++i]);
      }
      i++;
    }
  }

  cpgbeg(0,device,nx,ny);
  cpgsch(1.5);
  cpgscf(2);
  cpgsvp(0.15,0.85,0.15,0.85);

  if (pggray) {
    /* accumulate data */
    while (npts=get_stream(x,y)) {
      if (first) {
	ntwod=npts;
	twod=(float *) malloc(sizeof(float)*ntwod*ntwod);
	first=0;
      }
      if (ntwod != npts) error_message("number of bins changed!");
      minmax(y,npts,&ymin,&ymax);
      if (ymin<grmin) grmin=ymin;
      if (ymax>grmax) grmax=ymax;
      for (i=0; i<npts; i++) {
	twod[itwod++]=y[i];
	if (itwod > ntwod*ntwod) {
		puti(itwod);
	puti(ntwod*ntwod);
		error_message("array overflow in pgplotter");
	}
      }
      nlines++;
    }
    cpgswin(1.0,(float)ntwod+1.0,1.0,(float)nlines+2.0);
    if (zmin == 0.0) zmin=grmin;
    if (zmax == 0.0) zmax=grmax;
    cpggray(twod,ntwod,nlines,1,ntwod,1,nlines,zmax,zmin,tr);
    //printf("Grey scale range: %f %f\n",zmin,zmax);
    cpgsch(1.0);
    if (yymin == yymax) {
      yymin = zzmin;
      yymax = zzmax;
      subints = 1;
    } else {
      yymin/=1000.0;
      yymax/=1000.0;
      subints = 0;
    }
    if ((xxmin == 0.0) && (xxmax == 0.0)) {
      cpgmtxt("B",2.8, 0.5, 0.5, "Sample number");
    } else {
      cpgmtxt("B",2.8, 0.5, 0.5, "Time (ms)");
    }
    if (xxmin == 0.0) xxmin=1.0;
    if (xxmax == 0.0) xxmax=(float) ntwod;
    cpgswin(xxmin,xxmax,yymin,yymax); 
    cpgbox("bcnst",0.0,0,"bvcnst",0.0,0);
    cpgmtxt("T",2.8, 0.5, 0.5, title);
    if (subints) 
      cpgmtxt("L",3.4, 0.5, 0.5, "Time (s)");
    else
      //cpgmtxt("L",3.4, 0.5, 0.5, "Frequency (GHz)");
      cpgmtxt("L",3.4, 0.5, 0.5, "Frequency (MHz)");
    if (pdm != 0.0) {
      cpgsci(0);
      cpgslw(4);
      for (i=0; i<noffset; i++) {
	cpgmove(offset[i],yymax); 
	for (f=yymax; f>yymin; f-=0.001) {
	  deltat=4.148808*(pow(f,-2.0)-pow(yymax,-2.0))*pdm;
	  cpgdraw(offset[i]+deltat,f);
	}
      }
    }
  } else {
    /* do a regular plot */
    while (npts=get_stream(x,y)) {
      if (verbose) printf("read %d points... (max %d)\n",npts,maxpts);
      minmax(x,npts,&xmin,&xmax);
      minmax(y,npts,&ymin,&ymax);
      if (verbose) printf("xmin %f xmax %f ymin %f ymax %f\n",xmin,xmax,ymin,ymax);
      cpgask(0);
      cpgpage();
      cpgswin(xmin,xmax,ymin,ymax);	
      //cpgslw(4);
      cpgbox("bcnst",0.0,0,"bcnst",0.0,0);
      cpgmtxt("L",3.4, 0.5, 0.5, "Flux density (aU)");
      cpgmtxt("B",2.8, 0.5, 0.5, "Bin number");
      cpgline(npts,x,y);
    }
  }

  cpgend();
  exit(0);
}
