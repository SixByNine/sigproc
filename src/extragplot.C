#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "cpgplot.h"
#include "string.h"
#include <string>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
using namespace std;

// STILL TO FIX:

//   - When you're only reading in certain DMs or times, don't read/allocate
//     memory for whole file
//   - Make other plots
//   - Decide on SNR hist
//   - Add gaussian fit to SNR hist

void printhelp();
void timeavg(int n, float * d);
float getmax(float *data, int arraysize, int *index);
double normalise(int n, float * d, double * dataaverage);
int comparefunction(const void *a, const void *b);

//--------------------------------------------------------
//                        MAIN
//--------------------------------------------------------

int main (int argc, char *argv[]){

    //General purpose variables
    int j=0,ncands=0,nsubcands=0;
    long long int Sskip=0,Sread=0,i=1;
    char  textline[300];
    char *typedata;
    float *snrdata,*dmdata,plotdmlo=-1,plotdmhi=-1,snrlimit=6.5;
    int *peakdata,*boxcardata,boxplotnow,bestboxcardetected,bestpeakdetected;
    float bestdmdetected,bestsnrdetected=-9999,datamaxplotpoint[1];
    float datamin,datamax,snrmax=-9999,snrmin=9999,scalesnrlo=-1,scalesnrhi=-1;
    FILE  *infile;
    bool needfiles = true;

    //Top plot variables
    float *dmhist,*dmhistvalues;
    int* snrhist = new int[100];//do snrhist or no? split into RFI/pulses/gaussian?
    float *cdmdata,*cdmhist,*cdmhistvalues;

    //Data file's "header" lines
    char procfilename[100];
    char UTCdatetime[100];
    char ra[100],dec[100];
    float ctrfreq,bandwidth,tsamp,gthresh,dmlo,dmhi;
    int ndms;
    long long int nsamp;



//----------------PARSE INPUT LINE----------------
    while (i<argc) {
	if (fopen(argv[i],"r")!=NULL){
	    strcpy(textline,argv[i]);
	    infile=fopen(textline,"r");
	    if (infile==NULL){
		fprintf(stderr,"Error opening file %s\n",textline);
		exit(-1);
	    }
	    needfiles=false;
	    fclose(infile);
	}
	if (!strcmp(argv[i],"-s"))       sscanf(argv[++i],"%lld",&Sskip);
	if (!strcmp(argv[i],"-r"))       sscanf(argv[++i],"%d",&Sread);
	if (!strcmp(argv[i],"-dmlo"))    sscanf(argv[++i],"%f",&plotdmlo);
	if (!strcmp(argv[i],"-dmhi"))    sscanf(argv[++i],"%f",&plotdmhi);
	if (!strcmp(argv[i],"-scale")){
	    sscanf(argv[++i],"%f",&scalesnrlo);
	    sscanf(argv[++i],"%f",&scalesnrhi);
	}
	if (!strcmp(argv[i],"-h"))       {printhelp();exit(0);}
	i++;
    }

    if (i<2){
	printhelp();
	fprintf(stderr,"\n***ERROR in extragplot:\n\tNo values supplied on input line\n\n");
	exit(0);
    } else if (needfiles){
	fprintf(stderr,"\n***ERROR in extragplot:\n\tNo valid input data file given.\n\n");
	exit(0);
    }

    
//--------------READ DETECTION DATA FILE--------------
    long long int numlines = 0;
    infile = fopen(textline, "r");
    if (infile==NULL){
	fprintf(stderr,"***FATAL ERROR: File %s not found!\n",textline);
	exit(0);
    } else {
	while(fgets(textline, 300, infile)){
	    if(textline[0] == '#' || (textline[0] == 'f' && textline[1] == 'i')) continue;
	    else numlines++;
	}
	rewind(infile);
	numlines = numlines - 1;
	if ((snrdata = (float*) malloc (sizeof(float)*numlines)) == NULL){fprintf(stderr,"Error allocating RAM for SNR data\n");exit(-1);}
	if ((peakdata = (int*) malloc (sizeof(int)*numlines)) == NULL){fprintf(stderr,"Error allocating RAM for peak data\n");exit(-1);}
	if ((boxcardata = (int*) malloc (sizeof(int)*numlines)) == NULL){fprintf(stderr,"Error allocating RAM for boxcar data\n");exit(-1);}
	if ((dmdata = (float*) malloc (sizeof(float)*numlines)) == NULL){fprintf(stderr,"Error allocating RAM for dm data\n");exit(-1);}
	if ((typedata = (char*) malloc (sizeof(char)*numlines)) == NULL){fprintf(stderr,"Error allocating RAM for type data\n");exit(-1);}
	//detnumdata = (int*) malloc (sizeof(int)*numlines);


	// Read "header" info from file.
	fgets(textline, 300, infile);
	sscanf(textline,"%*s %s %lld %f %f %f %s %s %s %f %d %f %f",&procfilename[0],&nsamp,&tsamp,&ctrfreq,&bandwidth,&ra,&dec,&UTCdatetime,&gthresh,&ndms,&dmlo,&dmhi);
	dmhist = new float[ndms];
	dmhistvalues = new float[ndms];
	cdmhist = new float[ndms];
	cdmhistvalues = new float[ndms];


	// Do input checks
	if (Sread > nsamp || Sskip > nsamp){	   
	    fprintf(stderr,"WARNING in extragplot:\n\tSamples to skip/read exceeds samples in file.\n\tWill plot whole file (%lld samples).\n\n", nsamp);
	    Sread = nsamp - Sskip;
	    Sskip = 0;
	} else if (Sskip + Sread > nsamp){
	    fprintf(stderr,"WARNING in extragplot:\n\tSamples to skip+read exceeds samples in file.\n\tWill plot from skip (sample %lld) to end of file (sample %lld).\n\n",Sskip,nsamp);
	    Sread = nsamp - Sskip;
	} else if (Sread == 0){
	    Sread = nsamp;
	}

	// If plotdmlo/hi set wrong or unset, make plot range 1.1*search range size (or 0 to dmhi*1.1).
	if (plotdmlo == -1){
	    if (dmlo == 0) plotdmlo = dmlo;
	    else plotdmlo = dmlo-0.1*fabs(dmlo);
	}
	if (plotdmhi == -1) plotdmhi = dmhi*1.1;
	if (plotdmlo >= plotdmhi){
	    fprintf(stderr,"\nWARNING in extragplot:\n\tdmlo > dmhi. Will plot from minimum to maximum searched DM (%f - %f).\n\n",dmlo,dmhi);
	    if (dmlo == 0) plotdmlo = dmlo;
	    else plotdmlo = dmlo-0.1*fabs(dmlo);
	    plotdmhi = dmhi*1.1;
	}


	// READ IN DATA.
	i = 0;
	while(fgets(textline, 300, infile)){
	    if(textline[0] == '#') continue;
	    else {
		sscanf(textline,"%*s %*f %f %*lld %lld %*d %d %f %*d %c",&snrdata[i],&peakdata[i],&boxcardata[i],&dmdata[i],&typedata[i]); //,&detnumdata[i]);
		if (snrmax==-9999) snrmax = snrdata[i];
		if (bestsnrdetected==-9999) bestsnrdetected = snrdata[i];
		if (snrmin==9999 && snrdata[i] != 0) snrmin = snrdata[i];
		if (snrdata[i] > bestsnrdetected && typedata[i] == 'C'){
		    bestpeakdetected = peakdata[i];
		    bestdmdetected = dmdata[i];
		    bestboxcardetected = boxcardata[i];
		    bestsnrdetected = snrdata[i];
		}
		if (typedata[i] == 'c') nsubcands++;
		if (typedata[i] == 'G' || typedata[i] == 'R' || typedata[i] == 'C'){
		    if (snrdata[i] > snrmax){
			snrmax=snrdata[i];
		    }
		}
		if (snrdata[i] < snrmin && snrdata[i] !=0){
		    snrmin = snrdata[i];
		}
	    }
	    i++;
	}
	cdmdata = new float[nsubcands];

	// If scalesnrlo/hi set wrong or unset, make plot range sensible
	if (scalesnrhi == -1) {
	    for (i=0;i<numlines;i++){
		if (typedata[i] == 'C'){
		    if (snrdata[i]>scalesnrhi) scalesnrhi = snrdata[i];
		}
	    }
	}
	if (scalesnrlo == -1){
	    scalesnrlo = 0.1*scalesnrhi;
	    if (scalesnrlo < snrlimit) scalesnrlo = snrlimit;
	}
	if (scalesnrlo >= scalesnrhi){
	    scalesnrlo = snrlimit;
	    for (i=0;i<numlines;i++){
		if (typedata[i] == 'C'){
		    if (snrdata[i]>scalesnrhi) scalesnrhi = snrdata[i];
		}
	    }
	    if (scalesnrlo >= scalesnrhi) scalesnrlo = 0.5*scalesnrhi;
	    fprintf(stderr,"\nWARNING in extragplot:\n\tBad scale range. Will scale dots between %f and %f.\n\n",snrlimit,3*snrlimit);
	}
	if (bestsnrdetected<scalesnrlo) fprintf(stderr,"WARNING in extragplot:\n\tNo candidate points above SNR scale lower limit %f.\n\n",scalesnrlo);
    }
    fclose(infile);

    //cout<<"Plotting from "<<scalesnrlo<<" to "<<scalesnrhi<<"\n";
    //cout<<"Plotting from "<<plotdmlo<<" to "<<plotdmhi<<"\n";


//-----------------PLOT EVERYTHING----------------
//    cout<<Sskip*tsamp<<" and "<<Sread*tsamp<<" and "<<plotdmlo<<" and "<<plotdmhi<<"\n";
    sprintf(textline,"summary%s.ps/cps",procfilename);
    cpgopen(textline);
//    cpgopen("/xs");
//    cpgopen("?");


    //*****************************
    //***** DM vs. TIME PLOT ******
    //*****************************
    float xdummy[1];
    float ydummy[1];
    cpgsch(0.7);
    cpgsvp(0.05,0.99,0.05,0.5);
    fprintf(stderr,"Skipping %g seconds, reading %g seconds.\n",Sskip*tsamp,(Sskip+Sread)*tsamp);
    cpgswin((float)(Sskip*tsamp),(float)((Sskip+Sread)*tsamp),plotdmlo,plotdmhi);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgmtxt("L",2.1,0.5,0.5,"DM in pc/cm\\u3\\d");
    cpgmtxt("B",2.1,0.5,0.5,"Time (s)");
    float hiscale = 3;
    float loscale = 0.5;

    // Plot all detections (light gray)
/*    cpgsci(15);
    i=0;
    while (i<numlines){	if (typedata[i] == 'r' || typedata[i] == 'c' || typedata[i] == 'g'){
	    if (snrdata[i] >= scalesnrhi) cpgsch(hiscale);
	    else if (snrdata[i] <= scalesnrlo) cpgsch(loscale);
	    else cpgsch(loscale + hiscale*(snrdata[i]-scalesnrlo)/(scalesnrhi-scalesnrlo));
	    xdummy[0] = (float)peakdata[i]*tsamp;
	    cpgpt(1,xdummy,&dmdata[i],21);
	}
	i++;
	}*/

    // Plot RFI (dark gray) and "gaussian noise" (red)
    i=0;
    j=0;
    while (i<numlines){
	if (typedata[i] == 'C' || typedata[i] == 'c'){
	    i++;
	    if (typedata[i]=='c'){
		cdmdata[j]=dmdata[i];
		j++;
	    }
	    continue;
	} else if (typedata[i] == 'R'){
	    cpgsci(14);
	    boxplotnow = boxcardata[i];
	    i++;
	    continue;
	} else if (typedata[i] == 'G'){
	    cpgsci(2);
	    boxplotnow = boxcardata[i];
	    i++;
	    continue;
	} else if (boxplotnow == boxcardata[i]) {
	    if (snrdata[i] >= scalesnrhi) cpgsch(hiscale);
	    else if (snrdata[i] <= scalesnrlo) cpgsch(loscale);
	    else cpgsch(loscale + hiscale*(snrdata[i]-scalesnrlo)/(scalesnrhi-scalesnrlo));
	    xdummy[0] = (float)peakdata[i]*tsamp;
	    cpgpt(1,xdummy,&dmdata[i],21);
	}
	i++;
    }

    // Plot "candidates" (soft-blue); 4=blue, 11=softblue
    i=0;
    while (i<numlines){
	if (typedata[i] != 'C' && typedata[i] != 'c'){
	    i++;
	    continue;
	} else if (typedata[i] == 'C'){
	    ncands++;
	    cpgsci(11);
	    boxplotnow = boxcardata[i];
	    i++;
	    continue;
	} else if (boxplotnow == boxcardata[i]) {
	    if (snrdata[i] >= scalesnrhi) cpgsch(hiscale);
	    else if (snrdata[i] <= scalesnrlo) cpgsch(loscale);
	    else cpgsch(loscale + hiscale*(snrdata[i]-scalesnrlo)/(scalesnrhi-scalesnrlo));
	    xdummy[0] = (float)peakdata[i]*tsamp;
	    cpgpt(1,xdummy,&dmdata[i],21);
	}
	i++;
    }
    cpgsci(1);
    cpgsch(0.7);



    //*****************************
    //***** SNR vs. DM PLOT *******
    //*****************************
    cpgsch(0.7);
    cpgsvp(0.383,0.666,0.57,0.99); //plot mid
    //cpgsvp(0.716,0.99,0.57,0.99); //plot right
    cpgswin(plotdmlo,plotdmhi,snrmin,snrmax*1.1);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgmtxt("B",2.2,0.5,0.5,"DM in pc/cm\\u3\\d");
    cpgmtxt("L",2.1,0.5,0.5,"SNR");
    i=0;
    while (i<numlines){
	if (typedata[i] == 'C' || typedata[i] == 'c'){
	    i++;
	    continue;
	} else if (typedata[i] == 'R'){
	    cpgsci(14);
	    boxplotnow = boxcardata[i];
	    i++;
	    continue;
	} else if (typedata[i] == 'G'){
	    cpgsci(2);
	    boxplotnow = boxcardata[i];
	    i++;
	    continue;
	} else if (boxplotnow == boxcardata[i]) {
	    cpgpt(1,&dmdata[i],&snrdata[i],21);
	}
	i++;
    }

    i=0;
    while (i<numlines){
	if (typedata[i] != 'C' && typedata[i] != 'c'){
	    i++;
	    continue;
	} else if (typedata[i] == 'C'){
	    cpgsci(11);
	    boxplotnow = boxcardata[i];
	    i++;
	    continue;
	} else if (boxplotnow == boxcardata[i]) {
	    cpgpt(1,&dmdata[i],&snrdata[i],21);
	}
	i++;
    }
    cpgsci(1);
    cpgsch(0.7);



    //*****************************
    //******* N vs. DM PLOT *******
    //*****************************
    // sort the DM array and initialize counter
    for (i=0;i<ndms;i++){
	dmhist[i]=0;
	dmhistvalues[i]=0;
	cdmhist[i]=0;
	cdmhistvalues[i]=0;
    }
    int idm = 0,cidm=0;
    float ndmmax=1,ncdmmax=1;
    qsort(dmdata,numlines,sizeof(float),comparefunction);
    qsort(cdmdata,nsubcands,sizeof(float),comparefunction);
    dmhist[0] = 1;
    dmhistvalues[0]=dmdata[0];
    cdmhist[0] = 1;
    cdmhistvalues[0]=dmdata[0];
    for (i=0;i<numlines-1;i++){
	if (dmdata[i]<plotdmlo){
	    dmhistvalues[0] = dmdata[i];
	    continue;
	} else if (dmdata[i]>plotdmhi){
	    break;
	}
	if (dmdata[i]!=dmdata[i+1]){
	    idm++;
	    if (dmdata[i+1]>plotdmhi) break;
	    dmhistvalues[idm]=dmdata[i+1];
	    //dmhist[idm] = 0;
	} else {
	    dmhist[idm]++;
	}
	if (i<nsubcands-1){
	    if (cdmdata[i]<plotdmlo){
		cdmhistvalues[0] = cdmdata[i];
		continue;
	    } else if (cdmdata[i]>plotdmhi){
		break;
	    }
	    if (cdmdata[i]!=cdmdata[i+1]){
		cidm++;
		if (cdmdata[i+1]>plotdmhi) break;
		cdmhistvalues[cidm]=cdmdata[i+1];
		//cdmhist[cidm] = 1;
	    } else {
		cdmhist[cidm]++;
		if (cdmhist[cidm]>ndmmax) ndmmax = cdmhist[cidm];
	    }
	}
    }
    fprintf(stderr,"Found detections at %d unique DMs.\n",idm);
    fprintf(stderr,"Found candidates at %d unique DMs.\n",cidm);

    //Now actually make the plot...
    ndmmax = 1.1*ndmmax;
    cpgsvp(0.05,0.333,0.57,0.99); //plot left
    //cpgsvp(0.383,0.666,0.57,0.99); //plot mid
    cpgswin(plotdmlo,plotdmhi,0.0,ndmmax);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    //cpglab("Time (s)","DM (pc/cm3)","");
    cpgmtxt("B",2.2,0.5,0.5,"DM in pc/cm\\u3\\d");
    cpgmtxt("L",2.1,0.5,0.5,"N detections");
    cpgslw(3);
    cpgline(idm,dmhistvalues,dmhist);
    cpgsci(11);
    cpgslw(4);
    cpgline(cidm,cdmhistvalues,cdmhist);
    cpgsci(1);
    cpgslw(1);





    //*****************************
    //***** File/det Summary ******
    //*****************************
    cpgsvp(0.68,1.0,0.57,0.99);
    cpgswin(0,100,0,100);
    //cpgscf(2);
    cpgsch(0.9);
    cpgsah(1,180,0);
    cpgsls(4);
    //cpgsci(1); cpgarro(0,100,100,100); cpgsci(1);
    sprintf(textline,"\\frFile: \\fn%s",procfilename);
    cpgtext(0,94,textline);
    sprintf(textline,"\\frJ2000:\\fn %s %s",ra,dec);
    cpgtext(0,87, textline);
    sprintf(textline,"\\frTsamp\\fn = %.3f us",tsamp*1000);
    cpgtext(0,80, textline);
    sprintf(textline,"\\frf, BW\\fn = %.2f, %.2f MHz",ctrfreq,bandwidth);
    cpgtext(0,73, textline);

    cpgsci(1);cpgarro(0,70,100,70);cpgsci(1);
    cpgtext(0,65, "\\frSearch parameters:");
    sprintf(textline,"\\frNsigma\\fn = %.2f",gthresh);
    cpgtext(0,58, textline);
    sprintf(textline,"\\frDM range\\fn = %.1f,%.1f pc/cm\\u3",dmlo,dmhi);
    cpgtext(0,51, textline);

    cpgsci(1);cpgarro(0,40,100,40);cpgsci(1);
    cpgtext(0,35, "\\frBest detection:");
    sprintf(textline,"Det. at %f seconds",tsamp*bestpeakdetected);
    cpgtext(0,28, textline);
    sprintf(textline,"SNR = %.2f",bestsnrdetected);
    cpgtext(0,21, textline);
    sprintf(textline,"DM = %.2f",bestdmdetected);
    cpgtext(0,14, textline);
    sprintf(textline,"Boxcar = %d",bestboxcardetected);
    cpgtext(0,7, textline);
    
    
/*
    cpgtext(0.04,0.5, "Parameters at best SNR detection:");
    sprintf(textline,"SNR = %.2f",snrdata[maxind]);
    cpgtext(0.1,0.333, textline);
    sprintf(textline,"DM = %.2f",dmdata[maxind]);
    cpgtext(0.1,0.167, textline);
    sprintf(textline,"Boxcar = %.2f",boxcardata[maxind]);
    cpgtext(0.1,0.0, textline);
    sprintf(textline,"RA = %.2f",src_raj);
    cpgtext(0.48,0.333, textline);
    sprintf(textline,"Dec = %.2f",src_dej);
    cpgtext(0.48,0.167, textline);
    if (gal_b==0) sprintf(textline,"glat = not provided");
    else sprintf(textline,"glat = %.4f",gal_b);
    cpgtext(0.48,0.0, textline);
    cpgsci(1);
*/

    cpgclos();

    return(0);  
}






//---SUBFUNCTIONS---------------------------------------------
void printhelp(){
    fprintf(stderr,"\n\n****************************\n");
    fprintf(stderr,"***      EXTRAGPLOT      ***\n");
    fprintf(stderr,"****************************\n");
    fprintf(stderr,"Sarah Burke\nSwinburne University\n(c)2008\n");
    fprintf(stderr,"****************************\n\n");
    fprintf(stderr,"Multibeam plotting program to show results from gsearching software\n");
    fprintf(stderr,"\tUSAGE:\n\t\textragplot <filename> -<options>\n\n");
    fprintf(stderr,"NOTE: The in data file must be formatted as the Gresults output file from\n");
    fprintf(stderr,"a call to (dedisperse_all -G). Program will not accept wildcards.\n\n");
    fprintf(stderr,"Options are as follows:\n");
    fprintf(stderr,"\t-s     N   number of samples to skip [0]\n");
    fprintf(stderr,"\t-r     N   number samples to read [all]\n");
    fprintf(stderr,"\t-dmlo  N   lowest DM to plot [auto-choice]\n");
    fprintf(stderr,"\t-dmhi  N   highest DM to plot [auto-choice]\n");
    fprintf(stderr,"\t-scale N N lower, upper scale of SNR circles\n");
    fprintf(stderr,"\t-h         this help menu\n");
}

// This is a function to give qsort.
int comparefunction(const void *a, const void *b){
    return (*(int*)a - *(int*)b);
}


void timeavg(int n, float * d){
    int asize;
    float avg;
    for (int i=0;i<n-1;i+=2){
	avg = (d[i] + d[i+1]) / 2;
	asize = i/2;
	d[asize] = avg;
//	if (noff < 10){//	    fprintf (stderr,"averaged %g and %g and put %g into d[%d]\n",d[noff],d[noff+1],avg,asize);//	}//	scanf("%d",avg);
    }
    if (n-1%2==0) {
	asize = n/2;
	d[asize] = d[n];
    }
}

float getmax(float *data, int arraysize, int *index){
    float max=data[0];
    *index = 0;
    int i;
    for (i=0; i<arraysize; i++){
//	cout<<"index is "<<data[i]<<" \n"<<endl;
	if (data[i] > max){
//	    cout<<"\nFound "<<data[i]<<" larger than "<<max<<"\n";
	    max = data[i];
	    *index = i;
	}
    }
    return(max);
}


//  normalization + report of SIGMA and MEAN
double normalise(int n, float * d, double * dataaverage) {
	double sum=0.0;
	double sumsq=0.0;
	int noff=0;
	while (noff<n) {
		sum+=d[noff];
		sumsq+=d[noff]*d[noff];
		noff++;
	}
	double mean=sum/(double)n;
	*dataaverage = mean;
	double meansq=sumsq/(double)n;
	double sigma=sqrt(meansq-mean*mean);
	if (sigma==0) {
		fprintf(stderr, "Normalise::RMS of data is zero\n");
		return (sigma);
	}
	int i;
#pragma omp parallel for private(i)
	for (int i=0; i<n; i++)
		d[i]=(d[i]-mean)/sigma;
	return (sigma);
}


/*
double getrms(int n,float * d, double * dataaverage){
    double sum=0.0;
    double sumsq=0.0;
    for (int i=0;i<n;i++) {
	sum+=d[i];
	sumsq+=d[i]*d[i];
    }
    double mean=sum/(double)n;
    *dataaverage = mean;
    double meansq=sumsq/(double)n;
    double sigma=sqrt(meansq-mean*mean);
    if (sigma==0) {
	fprintf(stderr, "getrms::RMS of data is zero\n");
        return(sigma);
    }
    int i;
#pragma omp parallel for private(i)
    for (int i=0;i<n;i++) d[i]=(d[i]-mean);
//    fprintf(stderr,"IN GETRMS FUNCTION:\n\tRMS of data is: %g\n\tMean is: %g\n",sigma,mean);
    return(sigma);
    }*/
