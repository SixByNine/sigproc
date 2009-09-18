#include "cpgplot.h"
#include "string.h"
#include <string>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
using namespace std;
//extern "C" {
//#include "dedisperse.h"
//}


float filgetmax(float *data, int arraysize){
    float max=0; 
    int i;
    for (i=0; i<arraysize; i++){
	if (data[i] > max){
	    max = data[i];
	}
    }
    return(max);
}


float filgetmin(float *data, int arraysize){
    float min=0;
    int i;
    for (i=0; i<arraysize; i++){ if (data[i] < min) min = data[i]; }
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

void killchan(float *data, int ntime, int chan, int tscr, double filnbits, int filnchans){
  long long int bin=0;
  for (int t=0;t<ntime;t++){
      bin=(long long int)(t*filnchans)+chan;
    data[bin] = ((float)filnbits)/2;
  }
}

void filavg(int *ntime, int numchans, float *data){
    int averagedbin;
    float average;
    for (int t=0;t<*ntime;t+=2){
	for (int ch=0;ch<numchans;ch++){
	  average = (data[t*numchans+ch] + data[(t+1)*numchans+ch])/2;
	  averagedbin = ((t/2)*numchans+ch);
	  data[averagedbin] = average;
	}
    }
    if (numchans%2!=0) *ntime = (*ntime-1)/2; //ignores last odd channel
    else *ntime = (*ntime)/2;
}

//void filchanavg(int ntime, int *navchans, float *data){
//    int averagedbin;
//    float average;
//    for (int ch=0;ch<*navchans;ch+=2){
//      for (int t=0;t<ntime;t++){ //data is frequency-ordered
/*	  average = (data[ch*ntime+t] + data[(t+1)*numchans+ch])/2;
	  averagedbin = ((t/2)*numchans+ch);
	  data[averagedbin] = average;
	}
    }
    if (numchans%2!=0) *navg = (*navg-1)/2; //ignores last odd channel
    else *navg = (*navg)/2;
    }*/


void plotlambda(int Sdec, long long int Sskip, float dmoff, float inpDM, int filnchans, double filtsamp, double filfoff, double filfch1){
    float dtime[filnchans];
    float freqarr[filnchans];
    float tstep = filtsamp*Sdec;
    for(int i=0;i>filnchans*(int)filfoff;i+=(int)filfoff){
	dtime[i/(int)filfoff] = (Sskip*filtsamp+dmoff+0.03) + 2400*(pow(filfch1+i,-1)-pow(filfch1,-1));
	freqarr[i/(int)filfoff] = filfch1+i;
    }
    cpgline(filnchans,dtime,freqarr);
}

void plotdm(int Sdec, long long int Sskip, float dmoff, float inpDM, int filnchans, double filtsamp, double filfoff, double filfch1){
    float dtime[filnchans];
    float freqarr[filnchans];
    float tstep = filtsamp*Sdec;
    double fstep;
    for(int i=0;i<=filnchans;i++){
	dtime[(int)(fstep/filfoff)] = (Sskip*filtsamp+dmoff) + (4150000*(pow(filfch1+fstep,-2)-pow(filfch1,-2))*inpDM)/(1000);
	freqarr[(int)(fstep/filfoff)] = filfch1+fstep;
	fstep+=filfoff;
    }
    cpgline(filnchans,dtime,freqarr);
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



void printhelp(){
    fprintf(stderr,"\n\n****************************\n");
    fprintf(stderr,"***      QUICKGPLOT      ***\n");
    fprintf(stderr,"****************************\n");
    fprintf(stderr,"Sarah Burke\nSwinburne University\n(c)2008\n");
    fprintf(stderr,"****************************\n\n");
    fprintf(stderr,"Multibeam plotting program to show results from gsearching software\n");
    fprintf(stderr,"\tUSAGE:\n\t\tquickgplot <indatabase> -<options>\n\n");
    fprintf(stderr,"Inputs are as follows:\n");
    fprintf(stderr,"\t-s     number of samples to skip [0]\n");
    fprintf(stderr,"\t-r     number samples to read [all]\n");
    fprintf(stderr,"\t-dec   factor to scrunch data by [1, no scrunching]\n");
    fprintf(stderr,"\t-dm    DM to plot over waterfall as black line\n");
    fprintf(stderr,"\t-dmoff time (seconds) to offset the DM line from the first sample");
    fprintf(stderr,"\t-c     color table for plotting filterbank [def-heat]:\n");
    fprintf(stderr,"\t\t0 gray\n");
    fprintf(stderr,"\t\t1 inverse gray\n");
    fprintf(stderr,"\t\t2 heat\n");
    fprintf(stderr,"\t\t3 cold\n");
    fprintf(stderr,"\t\t4 plasma\n");
    fprintf(stderr,"\t\t5 forest\n");
    fprintf(stderr,"\t\t6 alien glow\n");
    fprintf(stderr,"\t\t7 pseudocolor\n");
    fprintf(stderr,"\tindatabase - base name of data (without extensions) as it comes out of gsearch\n");
    fprintf(stderr,"\t\t     code (no default). Program will not accept wildcards.\n");
}


//void plotfil(char *currentfile, long long int Sskip, int Sread, int Sdec, int inpDM){
// 
//}

/*
void plotfil(char *currentfile, long long int Sskip, int Sread, int Sdec, float inpDM){

    int   i=0,j=0,nFBsamps,filesizing,nsamps,sizeofheader;//,Sread=0,Sdec=1;
    long long int nsampsinfile;//Sskip=0,
    //    char  currentfile[100]; // for plotting the SNR vs DM plot
    float snrmax, snrmin,dmoff=0; // for plotting the SNR vs DM plot
    int   *archive,*rawarchive; // Filterbank colorplot files (also floatarchive)
    float *floatarchive,fchlast;
    FILE  *infile;
    bool isokay;
    unsigned char *buffer;
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
		  int andvalue = pow(2,nbits)-1;
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
    }
    fchlast = fch1 + (foff*nchans);

    nsamps = nFBsamps;
    if (isokay){
	for (int sdecnow=2;sdecnow<=Sdec;sdecnow*=2){
	    filavg(&nFBsamps,nchans,floatarchive);
	}
	killchan(floatarchive, nFBsamps, 6, Sdec);
	killchan(floatarchive, nFBsamps, 7, Sdec);
	killchan(floatarchive, nFBsamps, 8, Sdec);
	//	killchan(floatarchive, nFBsamps, 62, Sdec);
	//	killchan(floatarchive, nFBsamps, 63, Sdec);
	//	killchan(floatarchive, nFBsamps, 31, Sdec);
    }

    nsamps/=Sdec;
    Sskip/=Sdec;
    tsamp*=Sdec;

    float timeseriesindex[nsamps];
    for (i=0; i<nsamps; i++){
	timeseriesindex[i]=(Sskip + float(i))*tsamp;
    }

    

    cpgopen("115/xs");

    setcolortable(colortable);
    if (isokay){
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(timeseriesindex[0],timeseriesindex[nsamps-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff, 0.0};
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);
    }
    cpgslw(3);
    cpgsci(0);
    plotdm(Sdec,Sskip,dmoff,inpDM);
//    cpgsls(4); //dashed:2, dotted:4
//    cpgsci(4);
//    cpgslw(3);
//    plotlambda(Sdec,Sskip,dmoff,inpDM);

    cpgclos();

}

*/
