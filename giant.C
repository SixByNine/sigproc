#include "cpgplot.h"
#include "dialog.h"
#include "string.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void normalise(int n,float * d){
    float sum=0.0;
    float sumsq=0.0;
    int noff=0;
    while (noff<n) {
	sum+=d[noff];
	sumsq+=d[noff]*d[noff];
	noff++;
      }
    float mean=sum/(float)n;
    float meansq=sumsq/(float)n;
    float sigma=sqrt(meansq-mean*mean);
    for (int i=0;i<n;i++) d[i]=(d[i]-mean)/sigma;
}

void getminmax(int n, float * d, float * min, float * max){
  int i;
  *min=*max=d[0];
  for (i=0;i<n;i++) if (*min>d[i])*min=d[i];
  for (i=0;i<n;i++) if (*max<d[i])*max=d[i];
}

void getminmaxes(int n, float * d, int nplot, float * ymin, float * ymax){

  int i,j;
  if (n<nplot){
    for (i=0;i<n;i++) ymin[i]=ymax[i]=d[i];

  }
  else
    {
     for (j=0;j<nplot;j++){
       ymin[j]=ymax[j]=d[j*(n/nplot)];
       for (i=j*(n/nplot);i<(j+1)*(n/nplot);i++){
         if (ymin[j]>d[i])ymin[j]=d[i];
         if (ymax[j]<d[i])ymax[j]=d[i];
       }
       }
    }
}

/* scrunches by a factor 2 */

void bscrunch(int npoints, float * d){
  int i;
  for (i=0;i<npoints/2;i++) d[i]=(d[2*i]+d[2*i+1])/sqrt(2.0);
}

void plotminmax(int npoints, float * data){
  int i;
  /* first of all plot everything */
  float min,max,diff;
  float * x;
  float xmn,xmx,ymn,ymx;
  float xp1,yp,xp2;
  int xstart = 0;
  int xend = npoints;
  int npointsnow=npoints;
  int nplotnow=npoints;
  if (nplotnow>npoints)nplotnow=npoints;

  x = (float*)malloc(sizeof(float)*nplotnow);
  getminmax(npoints, &data[xstart], &min, &max);
  diff = max-min;
    for (i=0;i<nplotnow;i++) x[i]=(float)(i)/(npoints-1)*npoints;
    xmn=0.0;
    xmx=(float)npoints;
    ymn = min-diff*0.05;
    ymx = max+diff*0.05;
    printf("xmin xmax ymin ymax %f %f %f %f\n",xmn,xmx,ymn,ymx);
    cpgswin(xmn,xmx,ymn,ymx);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgsci(1);
    cpgline(nplotnow,x,data);
}

/* plots min and max for every npoints/nplot points */

void plotminmaxeff(int npoints, int nplot, float * data){

  if (nplot>npoints) {
    plotminmax(npoints,data);
  }else{
  int i;
  /* first of all plot everything */
  float min,max,diff;
  float * x, *ymin, *ymax;
  float xmn,xmx,ymn,ymx;
  float xp1,yp,xp2;
  int xstart = 0;
  int xend = npoints;
  int npointsnow=npoints;
  int nplotnow=nplot;
  if (nplotnow>npoints)nplotnow=npoints;

  x = (float*)malloc(sizeof(float)*nplotnow);
  ymin = (float *)malloc(sizeof(float)*nplotnow);
  ymax = (float *)malloc(sizeof(float)*nplotnow);

  getminmax(npointsnow, &data[xstart], &min, &max);
  getminmaxes(npointsnow, &data[xstart], nplotnow, ymin, ymax);
  diff = max-min;
    for (i=0;i<nplotnow;i++) x[i]=(float)(i)/(nplotnow-1)*nplotnow;
    xmn=0.0;
    xmx=(float)nplotnow;
    ymn = min-diff*0.05;
    ymx = max+diff*0.05;
    printf("xmin xmax ymin ymax %f %f %f %f\n",xmn,xmx,ymn,ymx);
    cpgswin(xmn,xmx,ymn,ymx);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgsci(1);
    cpgline(nplotnow,x,ymin);
    cpgline(nplotnow,x,ymax);
  }
}



/* A routine to find min of each nth block */

void fiddle(int npoints, int nplot, float * data){

  int i;
  /* first of all plot everything */
  float min,max,diff;
  float * x, *ymin, *ymax;
  float xmn,xmx,ymn,ymx;
  float xp1,yp,xp2;
  char ans=' ';
  int xstart = 0;
  int xend = npoints;
  int npointsnow=npoints;
  int nplotnow=nplot;

  x = (float*)malloc(sizeof(float)*nplot);
  ymin = (float *)malloc(sizeof(float)*nplot);
  ymax = (float *)malloc(sizeof(float)*nplot);

  while (ans!='q'){
    if (nplotnow>npointsnow) nplotnow=npointsnow;
    getminmax(npointsnow, &data[xstart], &min, &max);
    getminmaxes(npointsnow, &data[xstart], nplotnow, ymin, ymax);
    diff = max-min;
    for (i=0;i<nplotnow;i++) x[i]=(float)i+1;
    xmn=0.0;
    xmx=(float)nplotnow-1.0;
    ymn = min-diff*0.05;
    ymx = max+diff*0.05;
    printf("xmin xmax ymin ymax %f %f %f %f\n",xmn,xmx,ymn,ymx);
    cpgswin(xmn,xmx,ymn,ymx);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgline(nplotnow,x,ymin);
    cpgline(nplotnow,x,ymax);
    cpgcurs(&xp1,&yp,&ans);
    if (ans!=' ' && ans!='q' && ans!='b'){
     cpgsci(2);
     cpgmove(xp1,min);
     cpgdraw(xp1,max);
     xp2=xp1;
     cpgcurs(&xp2,&yp,&ans);
     cpgmove(xp2,min);
     cpgdraw(xp2,max);
     xstart=xstart+xp1/nplotnow*npointsnow;
     xend=xstart+xp2/nplotnow*npointsnow;
     npointsnow=(xp2-xp1)*npointsnow/nplotnow;
    } else{
      if (ans==' '){
         xstart=0;
         xend=npoints;
         npointsnow=npoints;
         nplotnow=nplot;
      }
      if (ans=='b'){
	bscrunch(npoints,data);
	 npoints=npoints/2;
         xstart=0;
         xend=npoints;
         npointsnow=npoints;
         nplotnow=nplot;
      }
    }
    cpgeras();
    cpgsci(1);
  }
}

class labels 
{
public:
  float x[100];
  float y[100];
  char * ascii[100];
  int nlabels;
  labels(int nl);
};

labels::labels(int nl){
  nlabels=0;
}

void min_means_min(float * min, float * max){
  if (*min>*max) {
    float temp = *max;
    *max=*min;
    *min = temp;
  }
}

int load_dunc_data(FILE * fptr,int nstart,int nwanted){
  static long int ntotal;
}
#define MAXFILES 32
#include "dedisperse.h"
char inpfile[80];
main (int argc, char *argv[]) 
{
  int ntimglobal=0;  // number of time samples in original
  int ngulp_original=0;       // number of time samples to look at at once
  int nskipstart=0;       // number skipped at start

  int i,ntim,headersize[MAXFILES],noff=0,gulp;
  float *time_series[MAXFILES],sum=0.0,sumsq=0.0,mean,meansq,sigma;
  int MAXMARKERS = 1024;
  int nfiles = 0;
  FILE *inputfile[MAXFILES];
  char filename[MAXFILES][256];

  if (argc<2 || help_required(argv[1])) {
    giant_help();
    exit(0);
  }
  print_version(argv[0],argv[1]);
  i=1;
  while (i<argc) {
    if (strings_equal(argv[i],"-s"))       sscanf(argv[++i],"%d",&nskipstart);
    if (strings_equal(argv[i],"-n"))       sscanf(argv[++i],"%d",&ngulp_original);
    if (file_exists(argv[i]))          {
      inputfile[nfiles]=open_file(argv[i],"r");
      strcpy(filename[nfiles],argv[i]);
      nfiles++;
    }
    if (nfiles>MAXFILES) error_message("too many open files");
    i++;
  }

  int ntimglobal_smallest=0, nsamp;
  for (i=0; i<nfiles; i++) {

    if ((headersize[i]=read_header(inputfile[i]))) {
      if (nbits!=32) 
	error_message("giant currently only works for 32-bit data");

      nsamp = nsamples(filename[i],headersize[i],nbits,nifs,nchans);
      if (i == 0) {
	ntimglobal_smallest=nsamp;
      } else {
	ntimglobal= nsamp;
	if (ntimglobal < ntimglobal_smallest) ntimglobal_smallest = ntimglobal;
      }
      
      // Space for data (time_series)
      time_series[i]=(float *) malloc(nsamp*sizeof(float));
      if (time_series[i]==NULL){
	fprintf(stderr,"Error mallocing %d floats of %d size\n",nsamp,
		sizeof(float));
	exit(-1);
      }
      
      // Skip data
      fprintf(stderr,"Skipping %d bytes\n",nskipstart*sizeof(float));
      fseek(inputfile[i],nskipstart*sizeof(float),SEEK_CUR);
      
    }
  }
  puti(ntimglobal_smallest);
  if (ngulp_original==0) ngulp_original=ntimglobal_smallest;

  cpgbeg(0,"/xs",1,1);
  /* create the dialog */
  dialog * d = new dialog();

  /* add the "action" buttons */
  int QUIT         = d->addbutton(0.02,0.65,"Quit");
  int PLOT         = d->addbutton(0.02,0.55,"Plot");
  int NEXT       = d->addbutton(0.02,0.45,"Next");
  int BSCRUNCH         = d->addbutton(0.02,0.35,"Bscrunch");
  int RESET         = d->addbutton(0.02,0.25,"Reset");
  int GLOBALRESET         = d->addbutton(0.02,0.15,"Global Reset");

  /* add the plot regions */
  d->addplotregion(0.2,0.99,0.85,0.9);
  float deltay = 0.7/(float)nfiles;
  for (i=0; i<nfiles; i++) 
      d->addplotregion(0.2,0.99,0.8-deltay*(float)(i+1),0.8-deltay*(float)i);

  d->draw();

  float x,y;
  char ans;
  int button=-1; int plotno=-1;
  int NPIXELS = 1024;
  float * xaxis = new float[NPIXELS];
  float * ymaxes = new float[NPIXELS];
  float * ymins = new float[NPIXELS];

  int nmarkers=0;
  int * markers= new int[MAXMARKERS];
  int nfileptr=nskipstart;
  int nplot=ngulp_original;
  int nstart=0;  
  int ngulp=ngulp_original;

  bool zoneplot=false;
  int ngates=0;
  float xgate=0.0;

  button=NEXT;
  while (button!=QUIT){
    // Plot the zone
    // Entire file is white
    if (button!=NEXT)button=d->manage(&x,&y,&ans,&plotno);
    if (button==BSCRUNCH) {
      button = plotno = -1;
    }
    if (button==GLOBALRESET) {
      plotno = -1;
      nstart = 0;
      nplot=ngulp_original;
      ngulp=ngulp_original;
      button=PLOT;
      // Skip to end of skipped data
      for (i=0; i<nfiles; i++)
	fseek(inputfile[i],-(nfileptr-nskipstart)*sizeof(float),SEEK_CUR);
      nfileptr=nskipstart;
      zoneplot=false;
      button=NEXT;
    }
    if (button==NEXT) {
      ngulp=ngulp_original;
      nstart=0;
      nplot=ngulp_original;
      // Read the data
      int NActuallyRead;
      for (i=0; i<nfiles; i++) {

	NActuallyRead = fread(time_series[i],sizeof(float),ngulp,inputfile[i]);
	puti(ngulp);
	if (NActuallyRead!=ngulp){
	  fprintf(stderr,"Could not read %d floats from file\n",ngulp);
	  ngulp = NActuallyRead;
	}
	normalise(ngulp,time_series[i]);
      }
      nfileptr+=ngulp;
      button = PLOT;
      plotno= -1;
      zoneplot=true;
    }
    if (button==RESET) {
      button = plotno = -1;
      nstart=0;
      nplot=ngulp;
      button=PLOT;
      zoneplot=true;
    }
    if (plotno>0){
      if (ans=='D'){
	markers[nmarkers]=(x/NPIXELS)*nplot+nstart+nfileptr-ngulp;
	nmarkers++;
	zoneplot=true;
      }
      if (ans=='A'){
      d->plotregions[plotno].reset();
      cpgsci(2);
      cpgmove(x,-1000);
      cpgdraw(x,1000);
      if (ngates==0){
	xgate=x;
	ngates++;
      } else {
	min_means_min(&x,&xgate);
	printf("x %f xgate %f\n",x,xgate);
	nstart=(int)(x/NPIXELS*nplot)+nstart;
	nplot=(xgate-x)/(NPIXELS)*nplot;
	//if (nplot<NPIXELS) nplot=NPIXELS;
	ngates=0;
	button=PLOT;
	zoneplot=true;
	printf("nplot %d nstart %d\n",nplot,nstart);
      }
      }
      if (ans=='z'){
	if (NPIXELS>nplot) {
	  nstart+=(int)x;
	}else
	nstart=(int)(x/(float)NPIXELS*nplot)+nstart;
	printf("nstart %d\n",nstart);
	nplot/=4;
	printf("nplot %d\n",nplot);
	nstart-=nplot/2;
	printf("nstart %d\n",nstart);
	//if (nplot<NPIXELS){nplot=NPIXELS;}
	button=PLOT;
	zoneplot=true;
      }
    }
    if (button==PLOT){
      for (i=0; i<nfiles; i++) {
	d->plotregions[i+1].erase(0.025,0.0,0.025,0.0);
	d->plotregions[i+1].set(0.0,(float)NPIXELS,0.0,1.0);
	plotminmaxeff(nplot,NPIXELS,&time_series[i][nstart]);
      }
      button = plotno = -1;
      zoneplot=true;
      ngates=0;
    }
    if (zoneplot){
      d->plotregions[0].set(0.0,(float)(ntimglobal_smallest),0.0,1.0);
      cpgsci(0);
      cpgsfs(1);
      cpgrect(0.0,(float)(ntimglobal_smallest),0.0,1.0);
      d->plotregions[0].set(0.0,(float)(ntimglobal_smallest),0.1,1.0);
      cpgsci(1);
      cpgsfs(1);
      cpgrect(0.0,(float)(ntimglobal_smallest),0.0,1.0);
      // Green is accessible
      cpgsci(3);
      cpgrect((float)nskipstart,(float)(ntimglobal_smallest),0.1,1.0);
      // Blue is in RAM
      cpgsci(4);
      cpgrect((float)nfileptr-ngulp,(float)(nfileptr),0.1,1.0);
      // Red is plotted
      cpgsci(2);
      cpgrect((float)nfileptr-ngulp+nstart,(float)(nfileptr-ngulp+nstart+nplot),
	      0.1,1.0);
      cpgsch(2.0);
      cpgslw(2);
      cpgsci(6);
      for (int i=0;i<nmarkers;i++){
	float xx=markers[i],yy=0.10;
	printf("MARKER %d\n",markers[i]);
	cpgpt(1,&xx,&yy,17);
      }
      cpgsch(1.0);
      cpgslw(1);
      zoneplot=false;
    }
  }
  cpgend();

  printf("Normal execution completed\n");
  return(0);

}
