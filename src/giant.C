#include "find_baseline.h"
#include "cpgplot.h"
#include "dialog.h"
#include "string.h"
#include "gtools.h"
#include "libplotfil.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

extern "C" {
  int newoldsumhrm_(float *, float *, int *, int *, float *,
		    float *, float *, float *, float *);
};

extern "C" {
  int newnewsumhrm_(float *, int *, int *, float *,
		    float *, float *, float *, float *);
};

extern "C" {
#include "dedisperse.h"
};

float * Creadspec(char * filename, int * npf, double * rate);
void find_formspec(int n, float * d);
void find_fft(int * n, float * d);
void normalise(int n, float * d);
void normalise(int n, float * d, float threshold);
void mowlawn(int n, float * d, int * mown, int * nrejects, float threshold, int maxnscrunch);
void getminmax(int n, float * d, float * min, float * max);
void getminmaxes(int n, float * d, int nplot, float * ymin, float * ymax);
void bscrunch(int npoints, float * d);
void plotminmax(int npoints, float * data, float tstart, float delta);
void plotminmaxeff(int npoints, int nplot, float * data, float tstart, float delta);
void min_means_min(float * min, float * max);
int load_dunc_data(FILE * fptr,int nstart,int nwanted);
int load_killdata(int * killdata,int kchans,char * killfile);
//void automow(int n, float * d, int * mown, int * nrejects, float threshold);

void plotfil(char *currentfile, long long int Sskip, int Sread, int Sdec, float inpDM, int pgpID, bool dokill, char *killfile);
void formpdf(float * pdfs, int pdfmax, int ngulp, float * time_series);
void  zap_improbables(float pdfs[][100], float ** time_series, int nbeam,
		     int ngulp, int maxsigma, float thresh, int nbeammax);


   class labels 
   {
   public:
     float x[100];
     float y[100];
     char * ascii[100];
     int nlabels;
     labels(int nl);
   };

   #define MAXFILES 32
//   #define MAXSIGMA 10
//test!
   #define MAXSIGMA 100
//

   char inpfile[80];

/* sigproc global variables describing the operating mode */
int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, 
  baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off;
//double refrf,userdm,fcorrect;
//float clipvalue,jyf1,jyf2;

int main (int argc, char *argv[]) 
{
  int ntimglobal=0;  // number of time samples in original
  int ngulp_original=0;       // number of time samples to look at at once
  int nskipstart=0;       // number skipped at start
  int nrejects;
  int zapswitch=0;
  double tsamp_orig=0;

  //gsearch setup & defaults
  float Gsigmacut=7.0;
  float delta, tstart;
  vector<Gpulse> * Giant = new vector<Gpulse>[MAXFILES];
  bool Gsearched=false;

  int i,j,ntim,headersize[MAXFILES],noff=0,gulp,toppeak;
  float *time_series[MAXFILES],sum=0.0,sumsq=0.0,mean,meansq,sigma;
  int MAXMARKERS = 1024;
  int nfiles = 0;
  FILE *inputfile[MAXFILES];
  char filename[MAXFILES][256];
  int spectra=0;
  int powerspectra=0;
  double dmoffirstfile;
  char *killfile;
  bool dokill=false;
  bool ssigned=true;
  int topfold=-1;
  int topgiant=-1;


  if (argc<2 || help_required(argv[1])) {
    fprintf(stderr,"Usage: giant filenames\n\t(e.g.>>  giant *.tim)\n\n\t-s  N\tskip N samples\n\t-n  N\tread N samples\n\t-S read spectra instead of amplitudes\n\t-i interpret signed chars as unsigned\n\t-z make a zap list of bad time samples\n");
    exit(0);
  }
  print_version(argv[0],argv[1]);
  i=1;
  while (i<argc) {
    if (file_exists(argv[i]))          {
      inputfile[nfiles]=open_file(argv[i],"r");
      strcpy(filename[nfiles],argv[i]);
      nfiles++;
    }
    if (strings_equal(argv[i],"-s"))       sscanf(argv[++i],"%d",&nskipstart);
    if (strings_equal(argv[i],"-S"))       spectra=1;
    if (strings_equal(argv[i],"-i"))       ssigned=false;
    if (strings_equal(argv[i],"-n"))       sscanf(argv[++i],"%d",&ngulp_original);
    if (strings_equal(argv[i],"-c"))       sscanf(argv[++i],"%f",&Gsigmacut);
    if (strings_equal(argv[i],"-z"))       zapswitch=1;
    if (strings_equal(argv[i],"-k"))       {killfile=(char*)malloc(strlen(argv[++i])+1); strcpy(killfile,argv[i]);dokill=true;}
    if (nfiles>MAXFILES) error_message("too many open files");
    i++;
  }

  int ntimglobal_smallest=0, nsamp;
  for (i=0; i<nfiles; i++) {

    if (spectra){
      int npf; 
      double rate;
      time_series[i]=Creadspec(filename[i],&npf,&rate);
      tsamp = 1.0/(rate);
      //normalise(npf,time_series[i]);
      nsamp = ntimglobal = ntimglobal_smallest = npf;
    }
    else
    {
    if ((headersize[i]=read_header(inputfile[i]))) {
      if (i==0) dmoffirstfile = refdm;
      if (nbits!=8 && nbits!=32)
	    error_message("giant currently only works for 8- or 32-bit data");

      nsamp = nsamples(filename[i],headersize[i],nbits,nifs,nchans);
      if (i == 0) {
	ntimglobal_smallest=nsamp;
      } else {
	ntimglobal= nsamp;
	if (ntimglobal < ntimglobal_smallest) ntimglobal_smallest = ntimglobal;
      }
      
      // Space for data (time_series)
      time_series[i]=(float *) malloc((nsamp+2)*sizeof(float));
      if (time_series[i]==NULL){
	fprintf(stderr,"Error mallocing %d floats of %d size\n",nsamp,
		sizeof(float));
	exit(-1);
      }
      tsamp_orig = tsamp;
      
      // Skip data
      fprintf(stderr,"Skipping %d bytes\n",nskipstart*nbits/8);
      fseek(inputfile[i],nskipstart*nbits/8,SEEK_CUR);
      
    } // each file
    } // spectra or not
  }  // for (i...)
  puti(ntimglobal_smallest);
  if (ngulp_original==0) ngulp_original=ntimglobal_smallest;

  // sbates 2009 ; switch to make a .killtchan file for time samples > 3.5 sigma

  int ngulp=ngulp_original;
//  int nrejects_max=ngulp_original/100;
  int * mown = new int[ngulp_original];
  int nstart=0;
  if (zapswitch){
    float dummy;
    int NActuallyRead;
    char *buffer;
    buffer = new char[ngulp*nbits/8];
    for (i=0; i<nfiles; i++){
      NActuallyRead = fread(buffer,nbits/8,ngulp,inputfile[i]);
      if (nbits==32){
	  memcpy(time_series[i],buffer,sizeof(float)*ngulp);
      } else {
	    for (int j=0;j<NActuallyRead;j++){
	      if (ssigned) time_series[i][j]=(float)buffer[j];
	      if (!ssigned) time_series[i][j]=(float)((unsigned char)buffer[j]);
	    }
      }
	    puti(ngulp);
	    //      find_baseline(time_series[i],ngulp,tsamp,2.0,3.0);
	    mowlawn(ngulp,time_series[i], mown, &nrejects,3.5,5);
    }
    printf("%f\n",dummy);
    printf("Bad time samples found...\n");
    exit(0);
    }

  int pgpID = cpgbeg(0,"/xs",1,1);
  /* create the dialog */
  dialog * d = new dialog();

  /* add the "action" buttons */
  int QUIT         = d->addbutton(0.02,0.95,"Quit");
  int POWER        = d->addbutton(0.07,0.85,"POWER");
  int SMHRM        = d->addbutton(0.075,0.80,"SMHRM");
  int FFT          = d->addbutton(0.02,0.85,"FFT");
  int PLOT         = d->addbutton(0.02,0.80,"Plot");
  int NEXT         = d->addbutton(0.02,0.75,"Next");
  int ZAPPEAK      = d->addbutton(0.075,0.75,"ZapPeak");
  int RESET        = d->addbutton(0.02,0.70,"Reset");
  int GLOBALRESET  = d->addbutton(0.02,0.65,"Global Reset");
  int HALVEPLOT    = d->addbutton(0.02,0.60,"Halve Plot");
  int BASELINE     = d->addbutton(0.02,0.50,"Baseline");
  int ZAPCOMMON    = d->addbutton(0.02,0.45,"Zap Common");
  int SUBTRACTMEAN = d->addbutton(0.02,0.40,"ZAP Mean");
  int BSCRUNCH     = d->addbutton(0.02,0.35, "Bscrunch");
  int NORMALISE    = d->addbutton(0.02,0.30,"Normalise"); 
  int HISTOGRAM    = d->addbutton(0.02,0.25,"Histogram"); 
  int GSEARCH      = d->addbutton(0.02,0.20,"Find Giants");
  int MOWLAWN      = d->addbutton(0.08,0.70, "LAWN");
  int SEEFIL       = d->addbutton(0.02,0.15,"View Band");
  int FWRITE       = d->addbutton(0.02,0.05,"Write File");
 

  /* add the plot regions */
  d->addplotregion(0.2,0.99,0.98,0.99);
  float deltay = 0.9/(float)nfiles;
  for (i=0; i<nfiles; i++) 
      d->addplotregion(0.2,0.99,0.95-deltay*(float)(i+1),0.95-deltay*(float)i);

  d->draw();

  float x,y;
  char ans;
  int button=-1; int plotno=-1;
  int NPIXELS = 1024;
  float * xaxis = new float[NPIXELS];
  float * ymaxes = new float[NPIXELS];
  float * ymins = new float[NPIXELS];

  int scrunch=1;
  int nmarkers=0;
  int * markers= new int[MAXMARKERS];
  int nfileptr=nskipstart;
  int nplot=ngulp_original;
  //  int nstart=0;  
  //int ngulp=ngulp_original;
  //  int nrejects_max=ngulp_original/100;
  float trialperiod;
  int doperiod=-1;
  float xperiod;

  //  int * mown = new int[nrejects_max];

  bool zoneplot=false;
  int ngates=0;
  float xgate=0.0;

  button=NEXT;
  if (spectra) button = PLOT;
  while (button!=QUIT){
    // Plot the zone
    // Entire file is white
    if (button!=NEXT)button=d->manage(&x,&y,&ans,&plotno);
//    printf("manage x %f y %f plotno %d\n",x,y,plotno);
    if (button==BASELINE) {
	for (i=0; i<nfiles; i++){
	  find_baseline(time_series[i],ngulp,tsamp,2.0,3.0);
	}
	button = PLOT;
	zoneplot=false;
        plotno = -1;
    }
    if (button==FWRITE) {
      // reread first header and close it. Sets globals.
      fclose(inputfile[0]);
      inputfile[0]=open_file(argv[1],"r");
      headersize[0]=read_header(inputfile[0]);
      output = open_file("giant.tim","w");
      nobits=32;
      nbands=1;
      dedisperse_header();
      fprintf(stderr,"Opened file, writing data\n");
      fwrite(time_series[0],sizeof(float),ngulp,output);
      fclose(output);
      button = -1;
      zoneplot=false;
      plotno =-1;
    }
    if (button==BSCRUNCH) {
	for (i=0; i<nfiles; i++){
	    bscrunch(ngulp,time_series[i]);
	}
	tsamp*=2;
	scrunch*=2;
      	ngulp/=2;
	nplot/=2;
	button = PLOT;
	zoneplot=false;
	Gsearched=false;
        plotno = -1;
    }
    if (button==FFT) {
	for (i=0; i<nfiles; i++){
	  ngulp = ngulp_original;
	  find_fft(&ngulp,time_series[i]);
	// Zap DC spike
	  time_series[i][0]=0.0;
	  time_series[i][1]=0.0;
	}
	spectra = 1;
	nplot = ngulp;
	button = PLOT;
	Gsearched=false;
        plotno = -1;
    }
    if (button==POWER) {
	for (i=0; i<nfiles; i++){
	  find_formspec(ngulp,time_series[i]);
	}
	ngulp/=2;
	powerspectra = 1;
	nplot = ngulp;
	button = PLOT;
        plotno = -1;
    }
    if (button==SMHRM) {
        nfiles = 6;
	for (i=1; i<nfiles; i++){
	  time_series[i]=(float *) malloc((ngulp+2)*sizeof(float));
	  if (time_series[i]==NULL){
	    fprintf(stderr,"Error allocating memory\n");
	    exit(-1);
	  }
	}
	for (i=1;i<nfiles;i++) memcpy(time_series[i],time_series[0],
				      (ngulp+2)*sizeof(float));
	d->nplotregion=1;
        float deltay = 0.9/(float)nfiles;
        for (i=0; i<nfiles; i++) 
          d->addplotregion(0.2,0.99,0.95-deltay*(float)(i+1),
			   0.95-deltay*(float)i);
	cpgeras();
	d->draw();

	float * workspace = new float[ngulp];

	// Set up space for data, now actually sumhrm
	int one=1;
			newoldsumhrm_(&time_series[0][1],workspace,&ngulp,&one,
	//		newoldsumhrm_(&time_series[0][0],workspace,&ngulp,&one,
		   time_series[1],time_series[2],time_series[3],
		   time_series[4],time_series[5]);
		/*	newnewsumhrm_(time_series[0],&ngulp,&one,
		   time_series[1],time_series[2],time_series[3],
		   time_series[4],time_series[5]);*/
		for (int iff=2;iff<6;iff++){
		  for (int i=0;i<ngulp;i++){
		    time_series[iff][i]/=sqrt(pow(2.0,(float)(iff-1)));
		  }
		}
	delete [] workspace;
	button = PLOT;
        plotno = -1;
    }
    if (button==NORMALISE) {
	for (i=0; i<nfiles; i++){
	  normalise(ngulp,time_series[i],3.0);
	}
	button = PLOT;
	Gsearched=false;
        plotno = -1;
    }
    if (button==HISTOGRAM) {
      float pdfs[nfiles][MAXSIGMA];
      //  create pdfs for each beam
      for (int i=0;i<nfiles;i++)
	formpdf(pdfs[i],MAXSIGMA,ngulp,time_series[i]);

      for (int i=0;i<nfiles;i++){
	for (int j=0;j<MAXSIGMA; j++){
	  fprintf(stderr, "pdfs[%d][%2d]=%8.0f %f \%\n", i, j+1, pdfs[i][j], 100*pdfs[i][j]/ngulp);
	}
      }
        button = PLOT;
        plotno = -1;
    }
    if (button==HALVEPLOT) {
	nplot/=2;
	button = PLOT;
	zoneplot=true;
	Gsearched=false;
        plotno = -1;
    }
    if (button==GLOBALRESET) {
      plotno = -1;
      nstart = 0;
      scrunch=1;
      tsamp = tsamp_orig;
      nplot=ngulp_original;
      ngulp=ngulp_original;
      button=PLOT;
      // Skip to end of skipped data
      for (i=0; i<nfiles; i++){
	fseek(inputfile[i],-(nfileptr-nskipstart)*nbits/8,SEEK_CUR);
	Giant[i].clear();
      }
      nfileptr=nskipstart;
      zoneplot=false;
      Gsearched=false;
      doperiod=-1;
      button=NEXT;
    }
    if (button==SUBTRACTMEAN && nfiles>1) {
      plotno = -1;
      nstart = 0;
      nplot=ngulp_original;
      ngulp=ngulp_original;
      button=PLOT;
      // Skip to end of skipped data
	for (int jj=0;jj<ngulp;jj++){
	  float sum;
	  sum=0.0;
	  for (i=1;i<nfiles;i++){
	    sum+=time_series[i][jj];
	  }
	  time_series[0][jj]-=sum/(float(nfiles-1));
	}
	Gsearched=false;
    }
    if (button==ZAPCOMMON && nfiles>1) {
      plotno = -1;
      nstart = 0;
      nplot=ngulp_original;
      ngulp=ngulp_original;
      button=PLOT;
      float pdfs[nfiles][MAXSIGMA];
      //  create pdfs for each beam
      for (int i=0;i<nfiles;i++)
	formpdf(pdfs[i],MAXSIGMA,ngulp,time_series[i]);
      //  for each point in each beam, mask if improbable
      float thresh = 3.0;
      int nbeammax = 5;
      zap_improbables(pdfs,time_series,nfiles,ngulp,MAXSIGMA,thresh,nbeammax);
      // Skip to end of skipped data
      //for (int jj=0;jj<ngulp;jj++){
      //  float sum;
      //  sum=0.0;
      //  for (i=1;i<nfiles;i++){
      //    sum+=time_series[i][jj];
      //  }
      //  time_series[0][jj]-=sum/(float(nfiles-1));
      //}
      //Gsearched=false;
    }
    if (button==NEXT) {
      ngulp=ngulp_original;
      nstart=0;
      nplot=ngulp_original;
      // Read the data
      int NActuallyRead;
      //      unsigned char *buffer;
      char *buffer;
      buffer = new char[ngulp*nbits/8];
      //buffer = new char[ngulp*nbits/8];
      for (i=0; i<nfiles; i++) {
//	NActuallyRead = fread(time_series[i],sizeof(float),ngulp,inputfile[i]);
	NActuallyRead = fread(buffer,nbits/8,ngulp,inputfile[i]);
	if (nbits==32){
	  memcpy(time_series[i],buffer,sizeof(float)*ngulp);
	} else {
	    for (int j=0;j<NActuallyRead;j++){
	      if (ssigned) time_series[i][j]=(float)buffer[j];
	      if (!ssigned) time_series[i][j]=(float)((unsigned char)buffer[j]);
	    }
	}
	
	puti(ngulp);
	if (NActuallyRead!=ngulp){
	  fprintf(stderr,"Could not read %d floats from file\n",ngulp);
	  ngulp = NActuallyRead;
	}
	if(nfiles==1){
	  // Add fake pulsar here....
	  //	  for (int ii=0;ii<ngulp;ii++) time_series[i][ii]+= 10.0*pow(sin(float(ii*2.0*M_PI/60.0)),250.0);
	}
	//normalise(ngulp,time_series[i]);
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
      Gsearched=false;
      if (ans=='p'){
	doperiod=-1;
      }
    }
    if (plotno>0){
      if (ans=='p'){  // hit p on a plot to type in a period
	d->plotregions[plotno].reset();
	//plot the thing;
	fprintf(stderr,"Please enter a period in seconds: ");
	cin>>trialperiod;
	xperiod = x;
	doperiod=plotno;
	button=PLOT;
      }
      if (ans==','){  // subtract 0.000001 seconds from period
	d->plotregions[plotno].reset();
	trialperiod-=0.000001;
	fprintf(stderr,"Trial period is now %f\n",trialperiod);
	doperiod=plotno;
	button=PLOT;
      }
      if (ans=='.'){  // add 0.000001 seconds to period
	d->plotregions[plotno].reset();
	trialperiod+=0.000001;
	fprintf(stderr,"Trial period is now %f\n",trialperiod);
	doperiod=plotno;
	button=PLOT;
      }
      if (ans=='<'){  // subtract 0.001 seconds from period
	d->plotregions[plotno].reset();
	trialperiod-=0.001;
	fprintf(stderr,"Trial period is now %f\n",trialperiod);
	doperiod=plotno;
	button=PLOT;
      }
      if (ans=='>'){  // add 0.001 seconds to period
	d->plotregions[plotno].reset();
	trialperiod+=0.001;
	fprintf(stderr,"Trial period is now %f\n",trialperiod);
	doperiod=plotno;
	button=PLOT;
      }
      if (ans=='X'){  // right click two points on a plot to calculate and plot a period
	d->plotregions[plotno].reset();
	cpgsci(3);
	cpgmove(x,-1000);
	cpgdraw(x,1000);
	if (ngates==0){
	  xgate=x;
	  ngates++;
	} else {
	  min_means_min(&x,&xgate);
	  printf("Period from %f to %f is %f\n",x,xgate,xgate-x);  
	  doperiod=plotno;
	  xperiod = x;
	  trialperiod=xgate-x;
	  ngates=0;
	  button=PLOT;
	}
      }
      if (ans=='D'){
	markers[nmarkers]=(int)(x/NPIXELS)*nplot+nstart+nfileptr-ngulp;
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
//	  printf("x %f xgate %f tstart %f\n",x,xgate,tstart);
	  nstart=(int)((x-tstart)/delta)+nstart;
	  nplot=(int)((xgate-x)/delta);
	  //if (nplot<NPIXELS) nplot=NPIXELS;
	  ngates=0;
	  button=PLOT;
	  zoneplot=true;
//	  printf("nplot %d nstart %d\n",nplot,nstart);
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
    if (button==GSEARCH){
      float topSNR=0.0;
      for (i=0;i<nfiles;i++){
	Giant[i] = giantsearch(nplot, &time_series[i][nstart], Gsigmacut, 1.0, 30, 1, refdm);
	for (int j=0;j<Giant[i].size();j++){
	  if (spectra) {
	    printf("%s %f\t SNR: %8.2f\t%8.3f Hz\t\t%i\t%d \n",filename[i],Giant[i][j].amp,Giant[i][j].SNR,(float)(Giant[i][j].loc)/nplot/2.0/tsamp,Giant[i][j].width,Giant[i][j].loc);
	    if (Giant[i][j].SNR>topSNR){
	      topSNR=Giant[i][j].SNR; toppeak=j; topfold=i;
	    }
	  }
	  else
	    printf("%s %f\t%f\t%i\t\t%i\n",filename[i],Giant[i][j].amp,Giant[i][j].SNR,Giant[i][j].loc,Giant[i][j].width);
	    Giant[i][j].start+=nstart;
	    Giant[i][j].loc+=nstart;
	}
//	if (!spectra)printf("%i GIANT CANDIDATES FOUND in %s\n\n",Giant[i].size(),filename[i]);
//	if (spectra) printf("%d peaks in fold %d\n",Giant[i].size(), (int)pow(2,i-1));
      }
//     printf("Top Pulsar is in fold %d SNR %f Freq %f Period %f ms\n",
//	     topfold, Giant[topfold][toppeak].SNR,(float)(Giant[i][j].loc)/nplot/2.0/tsamp,1.0/((float)(Giant[i][j].loc)/nplot/2.0/tsamp*1000.0) );
      button = plotno = -1;
      button=PLOT;
      Gsearched=true;
    }
    if (button==MOWLAWN){

      for (i=0; i<nfiles; i++){
	 mowlawn(ngulp,time_series[i], mown, &nrejects,3.5,5);
      }
      printf("Lawn mown\n");
      fflush(stdout);
      //      for (i=0;i<nrejects;i++){
      //printf("mown[%d] = %d\n",i,mown[i]);
      //}
      button = PLOT;
      Gsearched=false;
      plotno = -1;
      /*
      // try to read in the fil file
      char fbankfile[100];
      FILE *infile;
      int sizeofheader, filesizing, Sread, Sskip, Sdec, nsamps,k;
      int *rawarchive, *archive;
      float *floatarchive, fchlast;
      long long int nsampsinfile;
      bool isokay;
      

      sprintf(fbankfile,"%s.fil",filename[0]);
      if (file_exists(fbankfile)) {
	printf("file %s found\n",fbankfile);
	infile = fopen (fbankfile, "rb");
	if (infile==NULL) {fputs ("File error\n",stderr); exit (1);}
	
	sizeofheader = read_header(infile);
	if (data_type == 1){ //i.e., a normal sigproc binary profile.
	  nsampsinfile=nsamples(fbankfile,sizeofheader,nbits,nifs,nchans);
	  // printf("nsamps = %i, headsize = %i\n",nsampsinfile,sizeofheader);
	  Sskip = nstart*scrunch;
	  Sread = scrunch*nplot;
	  Sdec = scrunch;
	  if (Sread) nsampsinfile = Sread;

	  // allocate memory for the raw data
	  rawarchive = (int*) malloc(nifs*nchans*nsampsinfile*nbits/8);
	  if (rawarchive == NULL) {fputs ("Memory error",stderr); exit(2);}
	  
	    long long int nfseek = Sskip*nchans*nifs*nbits/8;
	    if (fseek(infile,nfseek,SEEK_CUR)) {fprintf(stderr,"\n\nFSEEK FAILED trying to skip %lld bytes\n\n",nfseek);exit(0);}


	  // read the data form inputfile to rawarchive memory
	  filesizing = fread(rawarchive,1,nifs*nchans*nsampsinfile*nbits/8,infile);
	  if(filesizing!=nsampsinfile*nchans*nbits*nifs/8) {fputs("Read error",stderr); exit(3);}
	  fclose(infile);
	  
	  //printf("filesizing = %i\n",filesizing);
	  archive = (int*) malloc(filesizing*sizeof(int)*8);
	  if (archive == NULL) {fprintf (stderr,"\nCould not malloc %f bytes for wise data.\n",filesizing*sizeof(int)*8); exit(2);}
	  
	  //---CONVERT BITWISE DATA---
	    filesizing = 0;
	    for (i=0;i<nsampsinfile*nchans*nifs*nbits/(sizeof(int)*8);i++){
		int abyte = rawarchive[i];
		for (j=0;j<(sizeof(int)*8/nbits);j++){
		  int andvalue = (int) pow(2,nbits)-1;
		  int bitshift = (int)nbits;
		  archive[filesizing++]= andvalue & abyte;
		  abyte = abyte >> bitshift;
		}
	    }
	  floatarchive = (float*) malloc(filesizing*sizeof(float));
	  floatarchive = filint2float(archive,filesizing);
	  for (i=0; i< filesizing; i++){
	    //printf("float = %f\n",floatarchive[i]);
	  }
	  //  printf("length = %i\n",i);
	}
	else {
	  fprintf(stderr,"File %s not found!\n",fbankfile);
	  isokay = false;
	  return -2;
	}
	fchlast = fch1 + (foff*nchans);
      }
      
      printf("mean being put in is; %f\n",mean);
      // try to change the frequency channels to the mean at all bad time samples (from mowlawn)
      for (i=0; i< Sread;i++){
      	for (j=0; j<nrejects; j++){
      	  if (i+Sskip == mown[j]){
      	    for (k=0; k<nchans; k++){
      	      floatarchive[i*nchans+k] = mean; // some mean value - will want to do something clever with a gaussian dist.
      	    }
      	  }
      	}
      }
      //      for (i=0; i<nrejects;i++)
      //	{
      //	  for (j=0; j<nchans; j++)
      //	    {
      //	      floatarchive[mown[i]*nchans+j] = mean;   
      //	    }
      //	}
      nsamps = nsampsinfile;
      nsamps/=Sdec;
      Sskip/=Sdec;
      tsamp*=Sdec;
      
      //try and plot this thing
      //      nsamps = nsampsinfile;
      //if (isokay){
	//	for (int sdecnow=2;sdecnow<=Sdec;sdecnow*=2){
	//  filavg(&nFBsamps,nchans,floatarchive);
	//}
	float timeseriesindex[nsamps];
	for (i=0; i<nsamps; i++){
	  timeseriesindex[i]=(Sskip + float(i))*tsamp;
	}
      
    cpgopen("?");
    //!!!    fprintf(stderr,"***Click filterbank window to return to Giant plot\n");
    int colortable=7;
    float snrmax, snrmin;
    setcolortable(colortable);
    if (isokay){
	snrmax = filgetmax(floatarchive,nsampsinfile);
	snrmin = filgetmin(floatarchive,nsampsinfile);
	cpgsvp(0.1,0.9,0.1,0.9);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff, 0.0};
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgimag(floatarchive,nchans,nsampsinfile,1,nchans,1,nsampsinfile,snrmin,snrmax,tr);
    }
    cpgslw(3);
    cpgsci(0);
    cpgclos();
 
    button = -1;
      //Gsearched=false;
    plotno = -1;
    */
    }

    if (button==SEEFIL){
      //      for (i=0;i<nfiles;i++){
      char fbankfile[100];
      sprintf(fbankfile,"%s.fil",filename[0]);
      //      fprintf(stderr,"plotfil %s %d %d %d %d %f\n",fbankfile,nstart,nplot*scrunch,scrunch,refdm);
      plotfil(fbankfile,nstart*scrunch,scrunch*nplot,scrunch,(float)dmoffirstfile,pgpID,dokill,killfile);
      button = plotno = -1;
      button=PLOT;
      zoneplot=true;
      Gsearched=true;
    }
    if (button==ZAPPEAK){
      char fbankfile[100];
      button = plotno = -1;
      button=PLOT;
    }
    if (button==PLOT){
      for (i=0; i<nfiles; i++) {
	delta = tsamp;
	if (spectra) delta = 1.0 / (2.0*tsamp*ngulp);
	//	if (powerspectra) delta*=2.0;
	tstart = delta * nstart;
	d->plotregions[i+1].erase(0.025,0.1,0.125,0.0);
	//	printf("nplot %d delta %f tstart %f\n",nplot,delta,tstart);
	d->plotregions[i+1].set(tstart,tstart+nplot*delta,0.0,1.0);
	//	printf("trace nplot %d NPIXELS %d start %f delta %f\n",nplot,NPIXELS,tstart,delta);
	plotminmaxeff(nplot,NPIXELS,&time_series[i][nstart],tstart,delta);
	//	printf("trace2\n");
	if (Gsearched){
	  cpgsci(3);
	  for (int j=0;j<Giant[i].size();j++){
	    float xx=(Giant[i][j].loc)*(tsamp);
	    if (spectra) xx = ((float)Giant[i][j].loc)/tsamp/nplot/2.0; // Hz?
	    cpgpt(1,&xx,&Giant[i][j].amp,25);
	  }
	}
      }
      if (doperiod>=0){
	d->plotregions[doperiod].reset();
	cpgsci(7);
	for (float pstep=xperiod;pstep<10*nplot;pstep+=trialperiod){
	  cpgmove(pstep,-1000);
	  cpgdraw(pstep,1000);
	}
	for (float pstep=xperiod;pstep>0;pstep-=trialperiod){
	  cpgmove(pstep,-1000);
	  cpgdraw(pstep,1000);
	}
      }
      //      doperiod=-1;
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


void normalise(int n, float * d){
  printf("Normalising %d samples\n",n);
    double sum=0.0;
    double sumsq=0.0;
    int noff=0;
    while (noff<n) {
	sum+=d[noff];
	sumsq+=d[noff]*d[noff];
	noff++;
      }
    double mean=sum/(double)n;
    double meansq=sumsq/(double)n;
    double sigma=sqrt(meansq-mean*mean);
    for (int i=0;i<n;i++) d[i]=(d[i]-mean)/sigma;
    printf("mean, %f sigma! %f\n",mean,sigma);
}

void normalise(int n, float * d, float threshold){
    printf("Normalising %d samples\n",n);
    double sum=0.0;
    double sumsq=0.0;
    int noff=0;
    while (noff<n) {
	sum+=d[noff];
	sumsq+=d[noff]*d[noff];
	noff++;
      }
    double mean=sum/(double)n;
    double meansq=sumsq/(double)n;
    double sigma=sqrt(meansq-mean*mean);
    for (int i=0;i<n;i++) d[i]=(d[i]-mean)/sigma;
    printf("Mean, %f sigma! %f\n",mean,sigma);
    // Now recompute sigma without threshold spikes
    sum = 0.0;
    sumsq=0.0;
    int nsum=0;
    noff=0;
    while (noff<n) {
      if (d[noff]<threshold){
	sum+=d[noff];
	sumsq+=d[noff]*d[noff];
	nsum++;
      }
      noff++;
    }
    mean=sum/(double)nsum; // shouldn't this be nsum?!?!?!?!
    meansq=sumsq/(double)nsum; // and this one!!!???
    sigma=sqrt(meansq-mean*mean);
    for (int i=0;i<n;i++) d[i]=(d[i]-mean)/sigma;
    printf("new mean, %f new sigma! %f\n",mean,sigma);
}


int quicksort_inplace_partition(int* array, int top, int bottom);
void quicksort_inplace(int* array, int top, int bottom){
	int middle;
	if (top < bottom)
	{
		middle = quicksort_inplace_partition(array, top, bottom);
		quicksort_inplace(array, top, middle);
		quicksort_inplace(array, middle+1, bottom);
	}
	return;
}


int quicksort_inplace_partition(int* array, int top, int bottom){
	int x = array[top];
	int topidx = top - 1;
	int botidx = bottom + 1;
	int swap;
	do
	{
		do
		{
			botidx --;
		}while (x <array[botidx]);

		do
		{
			topidx++;
		} while (x >array[topidx]);

		if (topidx < botidx)
		{
			swap = array[topidx];
			array[topidx] = array[botidx];
			array[botidx] = swap;
		}
	}while (topidx < botidx);
	return botidx;
}


void mowlawn(int n, float * d, int * mown, int * nrejects, float threshold, int max_nscrunch){
    printf("Mowing lawn %d samples\n",n);
    double sum=0.0;
    double sumsq=0.0;
    int noff=0;

    double meankeep;
    double sigmakeep;

    int nscrunch=1;
    int current_n=n;

     *nrejects=0;
    while(nscrunch < max_nscrunch){

     // Compute Mean and Sigma of timeseries
     while (noff<current_n) {
	 sum+=d[noff];
	 sumsq+=d[noff]*d[noff];
	 noff++;
       }

     double mean=sum/(double)current_n;
     double meansq=sumsq/(double)current_n;
     double sigma=sqrt(meansq-mean*mean);

     // Normalise timeseries
     for (int i=0;i<current_n;i++) d[i]=(d[i]-mean)/sigma;

     if(nscrunch==1){
       printf("Mean, %f sigma! %f\n",mean,sigma);
       meankeep = mean;
       sigmakeep = sigma;
     }


     // Now recompute sigma without spikes > threshold
     sum = 0.0;
     sumsq=0.0;
     int nsum=0;
     noff=0;
     while (noff<current_n) {
       if (d[noff]<threshold){
	 sum+=d[noff];
	 sumsq+=d[noff]*d[noff];
	 nsum++;
       }
       noff++;
     }

     mean=sum/(double)nsum; //was previosly noff
     meansq=sumsq/(double)nsum; // was previously noff
     sigma=sqrt(meansq-mean*mean);
     for (int i=0;i<current_n;i++) d[i]=(d[i]-mean)/sigma;
     printf("new mean, %f new sigma, %f\n",mean,sigma);
     noff =0;
     int i=0;
    
     while (noff<current_n) {
       if (fabs(d[noff]) > threshold){
	 //printf("d[%d]=%f\n",noff,d[noff]);
	 d[noff]= mean;
	 for(int z = 0 ; z < nscrunch; z++){
	   mown[*nrejects] = noff*nscrunch + z;
	   *nrejects = *nrejects+1;
	 }
	 i++;
       }

       noff++;
    /*   if ( *nrejects > (n/100)-1 ){
	 printf("Too many rejects %d (out of %d) change threshold......\n",*nrejects,n/100-1);
	 exit(-1);
       }*/
     }



     bscrunch(current_n,d);
     nscrunch*=2;
     current_n /=2;

    }


    FILE *chantkill;
    if ((chantkill = fopen("tchan.kill","w")) == NULL){
      printf("Error opening file\n");
      exit(-4);
    }
    fprintf(chantkill,"#%lf\t%lf\n",meankeep,sigmakeep);

    
    quicksort_inplace(mown,0,*nrejects);
//    for (int i =0 ; i < *nrejects; i++){
//	    printf("%d\n",mown[i]);
  //  }
    int currzap=0;
    for (int i =0 ; i < n; i++){
      // Skip on, writing 1s, until we reach the next thing that was mown
      if(currzap >= *nrejects || i < mown[currzap]){
	fprintf(chantkill,"1\n",i);
      } else {
	// When we reach something to mow, write a zero
	fprintf(chantkill,"0\n",i);
	currzap++;
	// in case we mow the same sample more than once, we check that we have moved on
	// to a larger sample to mow, otherwise we try the next entry in the mown array.
	while(currzap < *nrejects && mown[currzap]<=i)currzap++;
      }
    }
    


    fclose(chantkill);

}


void getminmax(int n, float * d, float * min, float * max){
  int i;
  *min=*max=d[0];
  for (i=0;i<n;i++) if (*min>d[i])*min=d[i];
  for (i=0;i<n;i++) if (*max<d[i])*max=d[i];
}

void getminmaxes(int n, float * d, int nplot, float * ymin, float * ymax){

  // problems with int overflows in this routine for large n, nplot.
  int i,j;
  if (n<nplot){
    for (i=0;i<n;i++) ymin[i]=ymax[i]=d[i];
  }
  else
    {
     for (j=0;j<nplot;j++){
       int index = (int) ((double)j*(double) n/double(nplot));
       ymin[j]=ymax[j]=d[index];
       for (i=index;i<index+n/nplot;i++){
         if (ymin[j]>d[i])ymin[j]=d[i];
         if (ymax[j]<d[i])ymax[j]=d[i];
       }
       }
    }
}

/* scrunches by a factor 2 */

void bscrunch(int npoints, float * d){
  int i;
  for (i=0;i<npoints/2;i++) d[i]=(d[2*i]+d[2*i+1])/2.0;
}

void plotminmax(int npoints, float * data, float tstart, float delta){
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

  //  printf("plotminmax entered\n");

  x = (float*)malloc(sizeof(float)*nplotnow);
  if (x==NULL){
    fprintf(stderr,"Error allocating %d bytes in plotminmax\n",
	    sizeof(float)*nplotnow);
    exit(-1);
  }
  //  printf("getminmax\n");
  getminmax(npoints, &data[xstart], &min, &max);
  diff = max-min;
  for (i=0;i<nplotnow;i++) x[i]=tstart + (float)(i) * delta;
    xmn=tstart;
    xmx=xmn + delta * (float)nplotnow;
    ymn = min-diff*0.05;
    ymx = max+diff*0.05;
//    printf("plotminmax xmin xmax ymin ymax %f %f %f %f\n",xmn,xmx,ymn,ymx);
    cpgswin(xmn,xmx,ymn,ymx);
    cpgsch(0.5);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgsci(1);
    cpgsch(1.0);
    cpgline(nplotnow,x,data);
}

/* plots min and max for every npoints/nplot points */

void plotminmaxeff(int npoints, int nplot, float * data, float tstart,
		   float delta){

  //  fprintf(stderr,"plotminmaxeff\n");
  if (nplot>npoints) {
    //    fprintf(stderr,"plotminmaxeff short\n");
    plotminmax(npoints, data, tstart, delta); // fixed ??
  }else{
  int i;
  /* first of all plot everything */
  float min,max,diff;
  float * x, * ymin, * ymax;
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
  if (x==NULL || ymin ==NULL || ymax==NULL){
    fprintf(stderr,"Error allocating memory in plotminmaxeff\n");
    exit(-1);
  }

  //  fprintf(stderr,"bing\n");
  getminmax(npointsnow, &data[xstart], &min, &max);
  //  fprintf(stderr,"bing\n");
  getminmaxes(npointsnow, &data[xstart], nplotnow, ymin, ymax);
  //  fprintf(stderr,"bing\n");
  diff = max-min;
  for (i=0;i<nplotnow;i++) x[i]=tstart + delta * npoints * 
       (float)(i)/((float)nplotnow-1);
    xmn=tstart;
    xmx=tstart+(float)npoints*delta;
    ymn = min-diff*0.05;
    ymx = max+diff*0.05;
    //    printf("xmin xmax ymin ymax %f %f %f %f\n",xmn,xmx,ymn,ymx);
    cpgswin(xmn,xmx,ymn,ymx);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgsci(1);
    cpgline(nplotnow,x,ymin);
    cpgline(nplotnow,x,ymax);
  }
}

/*
void formpdf(float * pdf, int pdfmax, int ngulp, float * time_series){
  //negative values are ignored

    normalise(ngulp,time_series);
    //zero pdf
    for (int i=0;i<pdfmax;i++)
      pdf[i]=0.0;
    
    for(int i=0;i<ngulp;i++){
      int index = (int)time_series[i];
      if(index<0) index=0;
      if(index>pdfmax-1) index=pdfmax-1;
      pdf[index]++;
    }
}
*/

//test!
void formpdf(float * pdf, int pdfmax, int ngulp, float * time_series){
  //negative values are ignored
  int negval = 0;
    normalise(ngulp,time_series);
    //zero pdf
    for (int i=0;i<pdfmax;i++)
      pdf[i]=0.0;
    
    for(int i=0;i<ngulp;i++){
      float tentimes = 10*time_series[i];
      int index = (int)tentimes;
      //int index = (int)time_series[i];
      if(index<0){ 
	index=0;
	//	fprintf(stderr,"Negative value!");
	negval++;
      }
      if(index>pdfmax-1) index=pdfmax-1;
      pdf[index]++;
    }
    fprintf(stderr,"Number of negative values: %d\n", negval);
    fprintf(stderr,"ngulp = %d\n", ngulp);
}


void zap_improbables(float pdfs[][100], float ** time_series, int nbeam,
		     int ngulp, int maxsigma, float thresh,
		     int nbeammax){
  //  nbeam is the total number of beams, ngulp is the number of points in each gulp, maxsigma is the same as pdfmax, 
  //  thresh is the treshold value for sigma accepted, nbeammax is the maximum number of beams with high sigma accepted
  //  time_series[BEAMID][PT#]

  int zapcount = 0;
    // For each point in time
    for (int i=0;i<ngulp;i++){
      int suspect = 0;    // Innocent until proven guilty
      int nbad = 0;

      // For each beam of that time
      for (int ib=0;ib<nbeam;ib++){
	if(time_series[ib][i]>thresh) nbad++;
      }
      if (nbad>=nbeammax){
	for (int ib=0;ib<nbeam;ib++){
	  if(time_series[ib][i]>thresh){ 
	    time_series[ib][i]=0.0;	  
	    zapcount++;
	  }
	}
      }
    }
    fprintf(stderr, "Number of zapped: %d\n", zapcount);
}


/*
// A routine to find min of each nth block 
void fiddle(int npoints, int nplot, float * data){
  int i;
  // first of all plot everything 
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

*/


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


labels::labels(int nl){
  nlabels=0;
}



void plotfil(char *currentfile, long long int Sskip, int Sread, int Sdec, float inpDM, int pgpID,bool dokill,char *killfile){

  // Preserve global variables from tim file
  double temprefdm = refdm;
  double temptsamp = tsamp;
  int tempnbits = nbits;
  int tempnchans = nchans;
  int tempnifs = nifs;
  double tempfch1 = fch1;
  double tempfoff = foff;
  double temptstart = tstart;


    int   i=0,j=0,nFBsamps,filesizing,nsamps,sizeofheader;//,Sread=0,Sdec=1;
    long long int nsampsinfile;//Sskip=0,
    //    char  currentfile[100]; // for plotting the SNR vs DM plot
    float snrmax, snrmin,dmoff=0; // for plotting the SNR vs DM plot
    int   *archive,*rawarchive; // Filterbank colorplot files (also floatarchive)
    float *floatarchive,fchlast;
    FILE  *infile;
    bool isokay;
    //unsigned char *buffer;
    //    char *buffer;
    int colortable=7; //autoset to "pseudocolor"



//--------------READ FILTERBANK FILE--------------
    fprintf(stderr,"Looking for %s\t...\t",currentfile);
    
    if (file_exists(currentfile)) {
	fprintf(stderr,"file found!\n");
	
	infile = fopen (currentfile, "rb");
	if (infile==NULL) {fputs ("File error\n",stderr); exit (1);}
	
	sizeofheader = read_header(infile);
	if (data_type == 1){ //i.e., a normal sigproc binary profile.		    
	    nsampsinfile=nsamples(currentfile,sizeofheader,nbits,nifs,nchans);
	    if (Sskip > nsampsinfile) {fprintf(stderr,"\n\tERROR in quickgplot:\n\t\tSkipping %lld data samples will surpass whole fil file!\n\n",Sskip); exit(7);}
	    if (Sread) {
		nFBsamps = Sread;
	    } else {
		nFBsamps = nsampsinfile - Sskip;
	    }

	    // 8 throughout is the number of bits per byte; necessary because this is bitwise date
	    rawarchive = (int*) malloc(nifs*nchans*nFBsamps*nbits/8);
	    if (rawarchive == NULL) {fputs ("Memory error",stderr); exit(2);}
	    
	    long long int nfseek = Sskip*nchans*nifs*nbits/8;
	    if (fseek(infile,nfseek,SEEK_CUR)) {fprintf(stderr,"\n\nFSEEK FAILED trying to skip %lld bytes\n\n",nfseek);exit(0);}

	    // filesizing is the number of bytes read in from the raw data
	    filesizing = fread(rawarchive,1,nifs*nchans*nFBsamps*nbits/8,infile);
	    if(filesizing!=nFBsamps*nchans*nbits*nifs/8) {fputs("Read error",stderr); exit(3);}
	    fclose(infile);

	    archive = (int*) malloc(filesizing*sizeof(int)*8);
	    if (archive == NULL) {fprintf (stderr,"\nCould not malloc %f bytes for wise data.\n",filesizing*sizeof(int)*8); exit(2);}

	    //---CONVERT BITWISE DATA---
	    filesizing = 0;
	    for (i=0;i<nFBsamps*nchans*nifs*nbits/(sizeof(int)*8);i++){
		int abyte = rawarchive[i];
		for (j=0;j<(sizeof(int)*8/nbits);j++){
		  int andvalue = (int) pow(2,nbits)-1;
		    int bitshift = (int)nbits;
		    archive[filesizing++]= andvalue & abyte;
		    abyte = abyte >> bitshift;
		}
	    }
	    floatarchive = (float*) malloc(filesizing*sizeof(float));
	    floatarchive = filint2float(archive,filesizing);
	} else {fprintf(stderr,"File type %d is an unknown format for filterbank file\n",data_type); exit(6); }
	isokay = true;
    } else {
	fprintf(stderr,"File %s not found!\n",currentfile);
	isokay = false;
	return;
    }
    fchlast = fch1 + (foff*nchans);

    int *killchans = new int[nchans];
    if (dokill) {
	load_killdata(killchans,nchans,killfile);
    }


    nsamps = nFBsamps;
    if (isokay){
	for (int sdecnow=2;sdecnow<=Sdec;sdecnow*=2){
	    filavg(&nFBsamps,nchans,floatarchive);
	}
	
	if (dokill){
	    for (int ich=0;ich<nchans;ich++){
		if (killchans[ich] == 0){
//		    fprintf(stderr,"ZAPMASK CHANNEL %d\n",ich);
		    killchan(floatarchive,nFBsamps,ich,Sdec,nbits,nchans);
		}
	    }
	}
    }

    nsamps/=Sdec;
    Sskip/=Sdec;
    tsamp*=Sdec;

    float timeseriesindex[nsamps];
    for (i=0; i<nsamps; i++){
	timeseriesindex[i]=(Sskip + float(i))*tsamp;
    }

    

    cpgopen("?");
    //!!!    fprintf(stderr,"***Click filterbank window to return to Giant plot\n");

    setcolortable(colortable);
    if (isokay){
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff, 0.0};
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);
    }

    cpgslw(3);
    cpgsci(0);
    plotdm(Sdec,Sskip,dmoff,inpDM, nchans, tsamp, foff, fch1);
    cpgsci(1);
    

//    cpgsls(4); //dashed:2, dotted:4
//    cpgsci(4);
//    cpgslw(3);
//    plotlambda(Sdec,Sskip,dmoff,inpDM,nchans,tsamp,foff,fch1);

    cpgsci(1);
    fprintf(stderr,"\na = scrunchx2 in frequency (under development)\nk = kill current channel\nm = move dm curve\nh = kill lower 150 channels (under development)\nq = quit\n");
    float dumx,dumy;
    char dumchar;
    int favg=1;
    while(dumchar != 'q'){
      cpgcurs(&dumx,&dumy,&dumchar);
      if (dumchar == 'k'){
	fprintf(stderr,"Killing channel %d\n\n",(int)((fch1-dumy)/fabs(foff*favg))-1);
	killchan(floatarchive,nFBsamps,(int)((fch1-dumy)/fabs(foff*favg))-1,Sdec,nbits,nchans);
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff*favg, 0.0};
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);	
	cpgslw(3);
	cpgsci(0);
	plotdm(Sdec,Sskip,dmoff,inpDM, nchans, tsamp, foff*favg, fch1);
	cpgsci(1);
      }
      if (dumchar == 'h'){
	fprintf(stderr,"HITRUN data; killing lower 180 channels of band.\n\n");
	for (int j=0;j<150;j++){
	  killchan(floatarchive,nFBsamps,j,Sdec,nbits,nchans);
	}
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff*favg, 0.0};
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);
	cpgslw(3);
	cpgsci(0);
	plotdm(Sdec,Sskip,dmoff,inpDM, nchans, tsamp, foff*favg, fch1);
	cpgsci(1);
      }
      if (dumchar == 'm'){
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff*favg, 0.0};
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);	
	cpgslw(3);
	cpgsci(0);
	fprintf(stderr,"DM offset is now %f, dumx is %f\n\n",dumx - timeseriesindex[0],dumx);
	dmoff = dumx - timeseriesindex[0];
	plotdm(Sdec,Sskip,dmoff,inpDM, nchans, tsamp, foff*favg, fch1);
	cpgsci(1);
      }
      if (dumchar == 'a'){
	favg*=2;
//	int tempchans = nchans - 1;
	filavg(&nchans,nFBsamps,floatarchive);
//	nchans = tempchans;
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff*favg, 0.0};
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);
	cpgslw(3);
	cpgsci(0);
	plotdm(Sdec,Sskip,dmoff,inpDM, nchans, tsamp, foff*favg, fch1);
	cpgsci(1);
      }
    }
    cpgclos();
    cpgslct(pgpID);

    //restore original global variables
    temprefdm = refdm = temprefdm;
    tsamp =  temptsamp;  
    nbits =  tempnbits; 
    nchans = tempnchans;
    nifs =   tempnifs;
    fch1 =   tempfch1;
    foff =   tempfoff;
    tstart = temptstart;

    return;
}


int load_killdata(int * killdata,int kchans,char * killfile){
  FILE * kptr;
  char line[1];
  kptr = fopen(killfile,"r");
    if (kptr==NULL){
      fprintf(stderr,"Error opening file %s\n",killfile);
      exit(-2);
    }
    for (int i=0; i<kchans;i++) {
      if (fgets(line,20,kptr)!=NULL){  // Read in whole line                                                                     
        int nscanned = sscanf(line,"%d",&killdata[i]);
//	fprintf(stderr,"READ %d\n",killdata[i]);
        if (nscanned==0) {
          fprintf(stderr,"Could not scan %s as 1 or 0\n",line);
          exit(-1);
        }
      } else{
        fprintf(stderr,"Error reading %dth value from %s\n",i,killfile);
        exit(-1);
      }
    }
  fclose(kptr);
  return(0);
}
