#include "cpgplot.h"
#include "string.h"
#include <string>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
using namespace std;
/*extern "C" {
#include "dedisperse.h"
}*/


float filgetmax(float *data, int arraysize){
    float max=data[0]; 
    int i;
    for (i=0; i<arraysize; i++){
	if (data[i] > max){
	    max = data[i];
	}
    }
    return(max);
}


float filgetmin(float *data, int arraysize){
    float min=data[0];
    int i;
    for (i=0; i<arraysize; i++){
	if (data[i] < min){
	    min = data[i];
	}
    }
    return(min);
}


float* filint2float(int *array, int arraysize){
    float *floatarray = (float*) malloc(sizeof(float)*arraysize);
    char *temparray;
    for (int i=0;i<arraysize;i++){
	floatarray[i] = float(array[i]);
    }
    return(floatarray);
}


bool getkillfile(int * killmask,int nchans,char *killfile){
    FILE * kptr;
    char line[100];
    kptr = fopen(killfile,"r");
    if (kptr==NULL){
	fprintf(stderr,"GETKILLFILE: Error opening file %s\nWill NOT excise RFI channels.\n",killfile);
	return(false);
    }
    for (int i=0; i<nchans;i++) {
	if (fgets(line,20,kptr)!=NULL){  // Read in whole line
	    int nscanned = sscanf(line,"%d",&killmask[i]);
	    if (nscanned==0) {
		fprintf(stderr,"GETKILLFILE: Could not scan %s as 1 or 0\nWill NOT excise RFI channels.\n",line);
		return(false);
	    }
	} else {
	    fprintf(stderr,"GETKILLFILE: Error reading %dth value from %s\nWill NOT excise RFI channels.\n",i,killfile);
	    return(false);
	}
    }
    fclose(kptr);
    return(true);
}


void killchan(float *data, int ntime, int chan, float fillvalue, int filnchans){
  long long int bin=0;
  for (int t=0;t<ntime;t++){
      bin=(long long int)(t*filnchans)+chan;
      data[bin] = fillvalue;
  }
}

void filavg(int *ntime, int numchans, float *data){
    int averagedbin;
    float average;
    int tempntime = *ntime;
    if (tempntime%2!=0) tempntime--;
    for (int t=0;t<tempntime;t+=2){
        for (int ch=0;ch<numchans;ch++){
          average = (data[t*numchans+ch] + data[(t+1)*numchans+ch])/2;
          averagedbin = ((t/2)*numchans+ch);
          data[averagedbin] = average;
        }
    }
    *ntime = (tempntime)/2;
}

void filchanavg(int ntime, int *numchans, float *data){
    int averagedbin=0;
    float average;
    int tempchans = *numchans;
    for (int t=0; t<ntime; t++){
	for (int c=0; c<tempchans; c+=2){
	    average = (data[t*tempchans+c] + data[t*tempchans+c+1])/2;
	    data[averagedbin] = average;
	    averagedbin++;
	}
    }
    if ((tempchans % 2)!=0) *numchans = (*numchans-1)/2; //ignores last odd channel
    else *numchans = (*numchans)/2;
}

void plotlambda(long long int Sskip, float dmoff, float inpDM, int filnchans, double filtsamp, double filfoff, double filfch1){
    float dtime[filnchans];
    float freqarr[filnchans];
    for(int i=0;i>filnchans*(int)filfoff;i+=(int)filfoff){
	dtime[i/(int)filfoff] = (Sskip*filtsamp+dmoff+0.03) + 2400*(pow(filfch1+i,-1)-pow(filfch1,-1));
	freqarr[i/(int)filfoff] = filfch1+i;
    }
    cpgline(filnchans,dtime,freqarr);
}

void plotdm(long long int Sskip, float dmoff, float inpDM, int filnchans, double filtsamp, double filfoff, double filfch1){
    float dtime[filnchans];
    float freqarr[filnchans];
    double fstep;
    for(int i=0;i<=filnchans;i++){
	dtime[(int)(fstep/filfoff)] = (Sskip*filtsamp+dmoff) + (4150000*(pow(filfch1+fstep,-2)-pow(filfch1,-2))*inpDM)/(1000);
	freqarr[(int)(fstep/filfoff)] = filfch1+fstep;
	fstep+=filfoff;
    }
    cpgline(filnchans,dtime,freqarr);
}


// Task removes a frequency-dependent baseline.
void chanbaseline(float *data, int ntime, int filnchans){
    float sum = 0;
    float mean;
    for (int c=0; c<filnchans; c++){
	sum=0;
	for (int t=0; t<ntime; t++){
	    sum = sum + data[t*filnchans+c];
	}
	mean = sum/ntime;
	int t;
#pragma omp parallel for private(t)
	for (t=0; t<ntime; t++)
	    data[t*filnchans+c] = data[t*filnchans+c]/mean;
    }
    
}


//LIFTED FROM WILLEM'S PSRCHIVE CODE AND ALTERED SLIGHTLY
void setcolortable(int colortable)
{
  switch (colortable) {
    
      case 0: {    //greyscale
    float grey_l[] = { 0.0, 1.};
    float grey_r[] = { 0.0, 1.};
    float grey_g[] = { 0.0, 1.};
    float grey_b[] = { 0.0, 1.};
    cpgctab (grey_l, grey_r, grey_g, grey_b, 2, 1.0, 0.5);
    break;
      }
  
      case 1: {//inverse
    float grey_l[] = { 0.0, 1.0};
    float grey_r[] = { 1.0, 0.0};
    float grey_g[] = { 1.0, 0.0};
    float grey_b[] = { 1.0, 0.0};
    cpgctab (grey_l, grey_r, grey_g, grey_b, 2, 1.0, 0.5);
    break;
  }

      case 2: {//heat
    float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    cpgctab (heat_l, heat_r, heat_g, heat_b, 5, 1.0, 0.5);
    break;
  }
  
      case 3: {//cold
    float cool_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float cool_r[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float cool_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float cool_b[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    cpgctab (cool_l, cool_r, cool_g, cool_b, 5, 1.0, 0.5);
    break;
  }
  
      case 4: {//plasma
    float cool_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float cool_r[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float cool_g[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float cool_b[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    cpgctab (cool_l, cool_r, cool_g, cool_b, 5, 1.0, 0.5);
    break;
  }
  
      case 5: {//forest
    float cool_l[] = {0.0, 0.4, 0.7, 0.9, 1.0};
    float cool_r[] = {0.0, 0.0, 0.0, 0.3, 1.0};
    float cool_g[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float cool_b[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    cpgctab (cool_l, cool_r, cool_g, cool_b, 5, 1.0, 0.5);
    break;
  }
  
      case 6: {//alien glow
    float test_l[] = { 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
		       0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.};
    float test_r[] = { 0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
		       0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.};
    float test_g[] = { 0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
		       0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.};
    float test_b[] = { 0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
		       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.};
    cpgctab (test_l, test_r, test_g, test_b, 20, 1.0, 0.5);
    break;
  }
  
      case 7: { //pseudo
    float test_l[] = { -0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
    float test_r[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.};
    float test_g[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.};
    float test_b[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.};
    cpgctab (test_l, test_r, test_g, test_b, 9, 1.0, 0.5);
    break;
  }
  
  } //end switch loop
}





void plotfilhelp(){
    fprintf(stderr,"\na = scrunchx2 in frequency\n");
    fprintf(stderr,"k = kill channel at cursor position\n");
    fprintf(stderr,"m = start DM curve at x-projected cursor position\n");
    fprintf(stderr,"h = print this help menu\n");
    fprintf(stderr,"q = quit/return to main window\n");
    fprintf(stderr,"b = normalize frequency-dept bandpass ripples\n");
//    fprintf(stderr,"\n");
    return;
}




void plotfil(float *floatarchive, int nFBsamps, long long int Sskip, float inpDM, int pgpID, double filtsamp, int filnchans, double filfch1, double filfoff, int colortable, bool askdevice){

    int   i=0,j=0;
    float snrmax, snrmin,dmoff=0; // for plotting the SNR vs DM plot
    float fchlast;
    float timeseriesindex[nFBsamps];


    for (i=0; i<nFBsamps; i++){
	timeseriesindex[i]=(Sskip + float(i))*filtsamp;
    }
    fchlast = filfch1 + (filfoff*filnchans);

    if (askdevice){
	cpgopen("?");
    } else {
	cpgopen("/xs");
    }
    fprintf(stderr,"***Type 'h' for interactive-filterbank-window help!\n");

    setcolortable(colortable);
    snrmax = filgetmax(floatarchive,nFBsamps);
    snrmin = filgetmin(floatarchive,nFBsamps);
    cpgsvp(0.1,0.9,0.1,0.9);
    float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff, 0.0};
    cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);

    cpgslw(3);
    cpgsci(0);
    plotdm(Sskip,dmoff,inpDM, filnchans, filtsamp, filfoff, filfch1);
    cpgsci(1);
    

//    cpgsls(4); //dashed:2, dotted:4
//    cpgsci(4);
//    cpgslw(3);
//    plotlambda(Sdec,Sskip,dmoff,inpDM,filnchans,filtsamp,filfoff,filfch1);

    cpgsci(1);
    plotfilhelp();
    float dumx,dumy;
    char dumchar;
    int favg=1;
    while(dumchar != 'q'){
      cpgcurs(&dumx,&dumy,&dumchar);
      if (dumchar == 'k'){
	fprintf(stderr,"Killing channel %d\n\n",(int)((filfch1-dumy)/fabs(filfoff*favg))-1);
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	killchan(floatarchive,nFBsamps,(int)((filfch1-dumy)/fabs(filfoff*favg))-1,(snrmax+snrmin)/2,filnchans);
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff*favg, 0.0};
	cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);	
	cpgslw(3);
	cpgsci(0);
	plotdm(Sskip,dmoff,inpDM, filnchans, filtsamp, filfoff*favg, filfch1);
	cpgsci(1);
      }
/*      if (dumchar == 'h'){
	fprintf(stderr,"HITRUN data; killing lower 180 channels of band.\n\n");
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	for (int j=0;j<150;j++){
	    killchan(floatarchive,nFBsamps,j,(snrmax+snrmin)/2,filnchans);
	}
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);	float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff*favg, 0.0};
	cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);
	cpgslw(3);
	cpgsci(0);
	plotdm(Sskip,dmoff,inpDM, filnchans, filtsamp, filfoff*favg, filfch1);
	cpgsci(1);
	}*/
      if (dumchar == 'h'){
	  plotfilhelp();
      }
      if (dumchar == 'm'){
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff*favg, 0.0};
	cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);	
	cpgslw(3);
	cpgsci(0);
	fprintf(stderr,"DM offset is now %f, dumx is %f\n\n",dumx - timeseriesindex[0],dumx);
	dmoff = dumx - timeseriesindex[0];
	plotdm(Sskip,dmoff,inpDM, filnchans, filtsamp, filfoff*favg, filfch1);
	cpgsci(1);
      }
      if (dumchar == 'p'){
	      cpgclos();

	      snrmax = filgetmax(floatarchive,nFBsamps);
	      snrmin = filgetmin(floatarchive,nFBsamps);
	      cpgopen("/ps");
	      setcolortable(0);
	      cpgsvp(0.15,0.85,0.15,0.85);
	      cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	      cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	      cpglab("Time (s)", "Frequency (MHz)", "");
	      float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff*favg, 0.0};
	      cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);
	      cpgclos();
	      cpgopen("/xs");
	      setcolortable(colortable);
	      cpgsvp(0.1,0.9,0.1,0.9);
	      cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	      cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	      cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);

      }
      if (dumchar == 'a'){
	favg*=2;
//	int tempchans = filnchans;
	filchanavg(nFBsamps,&filnchans,floatarchive);
//	filnchans = tempchans;
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff*favg, 0.0};
	cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);
	cpgslw(3);
	cpgsci(0);
	plotdm(Sskip,dmoff,inpDM, filnchans, filtsamp, filfoff*favg, filfch1);
	cpgsci(1);
      }
      if (dumchar == 'b'){
	  chanbaseline(floatarchive,nFBsamps,filnchans);
	  snrmax = filgetmax(floatarchive,nFBsamps);
	  snrmin = filgetmin(floatarchive,nFBsamps);
	  cpgsvp(0.1,0.9,0.1,0.9);
	  cpgswin(timeseriesindex[0],timeseriesindex[nFBsamps-1],filfch1,fchlast);
	  cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	  float tr[] = {timeseriesindex[0], 0.0, filtsamp, filfch1, filfoff*favg, 0.0};
	  cpgimag(floatarchive,filnchans,nFBsamps,1,filnchans,1,nFBsamps,snrmin,snrmax,tr);
	  cpgslw(3);
	  cpgsci(0);
	  plotdm(Sskip,dmoff,inpDM, filnchans, filtsamp, filfoff*favg, filfch1);
	  cpgsci(1);
      }
    }
    cpgclos();
    cpgslct(pgpID);

    return;
}

void plotfil(float *floatarchive, int nFBsamps, long long int Sskip, float inpDM, int pgpID, double filtsamp, int filnchans, double filfch1, double filfoff, int colortable){
    plotfil(floatarchive,nFBsamps,Sskip,inpDM,pgpID,filtsamp,filnchans,filfch1,filfoff,colortable,false);
    return;
}
