#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "cpgplot.h"
#include "string.h"
#include "libplotfil.h"
#include <string>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
using namespace std;

// TO FIX:
//
//   - If the pulse is at the end of the file, make sure file reads from start-width to file end
//   - Channel scrunching if detection is faint
//   - Get correct pulse width

/*
Because I have no timeseries[0] (i.e. beam 1)? (YES.)
*/

void printbeam(int n, float min, float max, float *xdata, float *ydata, const char *beamid, int bestbeam, float *maxindex, float *maxdata, bool dataok, int start, int stop);
void printhelp();


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


//--------------------------------------------------------
//                        MAIN
//--------------------------------------------------------

extern "C" {
#include "dedisperse.h"
}
int main (int argc, char *argv[]){

    int   i=0,j=0,nFBsamps=0,filesizing,color,nsamp,nreadinsamp,sizeofheader,Sread=0,Sdec=1,datamaxindex;
    long long int Sskip=0,nsampsinfile,nsampsStart;
    char  currentfile[100], beam[2]; // for plotting the SNR vs DM plot
    float *snrplot,*dmplot, snrmax, snrmin, dmmax, tscrmax,maxsnrdetected=0; // for plotting the SNR vs DM plot
    int   *archive,*rawarchive,bestbeam=-1,maxind, *beamplot; // Filterbank colorplot files (also floatarchive)
    float *tscrplot, datamaxplotpoint[1],besttscr;
    float *floatarchive,datamin,datamax,fchlast;
    float *timeseries[16];                         // timeseries[0-12]:  beams 1 to D
    float bestfitdm = -1;
    FILE  *infile;                                 // timeseries[13-15]: dm0, dmlo, dmhi
    bool isokay[17];
    bool foundone = false;
    bool ISsigned = true,fsigned=false;
    bool dokill = false;
    char *killfile;
    unsigned char *buffer;
    long long int tempstart,temppeak;
    int tempwid;
    int latestend=0,maxpeakdetected=0;
    int colortable=7; //autoset to "pseudocolor"

//----------------PARSE INPUT LINE----------------
    while (i<argc) {
	if (strings_equal(argv[i],"-s"))       sscanf(argv[++i],"%lld",&Sskip);
	if (strings_equal(argv[i],"-r"))       sscanf(argv[++i],"%d",&Sread);
	if (strings_equal(argv[i],"-dec"))     sscanf(argv[++i],"%d",&Sdec);
	if (strings_equal(argv[i],"-c"))       sscanf(argv[++i],"%d",&colortable);
	if (strings_equal(argv[i],"-bestbeam"))sscanf(argv[++i],"%d",&bestbeam);
	if (strings_equal(argv[i],"-i"))       ISsigned = false;
	if (strings_equal(argv[i],"-f"))       fsigned=true;
	if (strings_equal(argv[i],"-k"))       {killfile=(char*)malloc(strlen(argv[++i])+1); strcpy(killfile,argv[i]);dokill=true;}
//	if (strings_equal(argv[i],"-dm"))      sscanf(argv[++i],"%f",&bestfitdm);
	if (strings_equal(argv[i],"-h"))       {printhelp();exit(0);}
	i++;
    }

    if (Sdec > Sread && Sread > 0){
	fprintf(stderr,"\n[***ERROR in quickgplot: Scrunch factor larger than samples read***]\n\n");
	exit(0);
    }
    if (i<2){
	printhelp();
	fprintf(stderr,"\n[***ERROR in quickgplot: No values supplied on input line***]\n\n");
	exit(0);
    } else if (argv[1][0]=='-'){
	printhelp();
	fprintf(stderr,"\n[***ERROR in quickgplot: Data string must precede options***]\n\n");
	exit(0);
    }
    
//-----------------READ TIMESERIES----------------
    for (i=0 ; i<16 ; i++){
	switch ( i ) {
	    case 9 :
		beam[0] = 'A'; break;
	    case 10 : 
		beam[0] = 'B'; break;
	    case 11 : 
		beam[0] = 'C'; break;
	    case 12 : 
		beam[0] = 'D'; break;
	    case 13 :
		sprintf(currentfile,"%s.dm0",argv[1]); break;
	    case 14 : 
		sprintf(currentfile,"%s.dmlo",argv[1]); break;
	    case 15 : 
		sprintf(currentfile,"%s.dmhi",argv[1]); break;
	    default : 
		sprintf(beam,"%d",i+1);
	}
	if (i<13) sprintf(currentfile,"%s.beam%s",argv[1],beam);
	fprintf(stderr,"Looking for %s\t...\t",currentfile);
	if (file_exists(currentfile)) {
	    fprintf(stderr,"file found!\n");
	    infile = fopen (currentfile, "r");
	    if (infile==NULL) {fputs ("File error\n",stderr); exit (1);}
	    
	    if (sizeofheader=read_header(infile)) {
		if (! fsigned){
		    if (isign > 0) {
			ISsigned=false;
			fprintf(stderr,"using signed header variable to set UNSIGNED\n");
		    }
		    if (isign < 0) {
			ISsigned=true;
			fprintf(stderr,"using signed header variable to set SIGNED\n");
		    }
		}
		nsampsinfile = nsamples(currentfile,sizeofheader,nbits,nifs,nchans);
		if (nbits!=32 && nbits!=8) {fprintf(stderr,"\n\n\tERROR in quickgplot.C:\n\t\tFile %s must be 8- or 32-bit data",currentfile); exit(7);}
		if (Sskip > nsampsinfile) {fprintf(stderr, "\n\tERROR in quickgplot:\n\t\tSkipping %d data samples will surpass whole data file!\n\n",Sskip); exit(7);}
		else nreadinsamp = nsampsinfile - Sskip;
		if (Sread){
		    if (Sread > nsampsinfile){
			fprintf(stderr, "\n\t**WARNING:\n\t\tRead samples (%d) > number of samples in file (%d). Setting Sread = %d\n\n",Sread,nsampsinfile,nsampsinfile);
			Sread = nsampsinfile;
		    } else if (Sread + Sskip > nsampsinfile){
			fprintf(stderr, "\n\t**WARNING:\n\t\tskipped + read samples (%d + %d) > number of samples in file (%d). Setting Sread = %d\n\n",Sskip,Sread,nsampsinfile,nsampsinfile-Sskip);
			Sread = nsampsinfile-Sskip;
		    }
		    nreadinsamp = Sread;
		}

		if (!foundone) { nsamp = nreadinsamp; }
		else { if (nreadinsamp < nsamp) nsamp = nreadinsamp; }

//		timeseries[i] = new float[nreadinsamp];
		timeseries[i] = (float*) malloc (sizeof(float)*nreadinsamp);
		if (timeseries[i] == NULL) {fputs ("Memory error",stderr); exit (2);}
		
		fseek(infile,Sskip*nbits/8,SEEK_CUR);
		
//		buffer = new unsigned char[nreadinsamp*nbits/8];
		buffer = (unsigned char*)malloc(sizeof(unsigned char)*(nreadinsamp*nbits/8));
		filesizing = fread(buffer,nbits/8,nreadinsamp,infile);
		if (filesizing != nreadinsamp) {fputs ("Read error",stderr); exit (3);}	
		
		if (nbits==32){
		    memcpy(timeseries[i],buffer,sizeof(float)*nreadinsamp);
		} else {
		    for (int j=0;j<nreadinsamp;j++){
			if (ISsigned) timeseries[i][j]=(float)((signed char)buffer[j]);
			if (!ISsigned) timeseries[i][j]=(float)buffer[j];
		    }
		}

		free(buffer);
		fclose(infile);
	    }
	    isokay[i] = true;
	    foundone = true;
	} else {
	    fprintf(stderr,"File not found!\n");
	    isokay[i] = false;
	}
    } // for all timeseries


//--------------READ FILTERBANK FILE--------------
    sprintf(currentfile,"%s.fil",argv[1]);
    fprintf(stderr,"Looking for %s\t...\t",currentfile);
    
    if (file_exists(currentfile)) {
	fprintf(stderr,"file found!\n");
	
	infile = fopen (currentfile, "rb");
	if (infile==NULL) {fputs ("File error\n",stderr); exit (1);}
	
	sizeofheader = read_header(infile);
	if (data_type == 1){ //i.e., a normal sigproc binary profile.		    
	    nsampsinfile=nsamples(currentfile,sizeofheader,nbits,nifs,nchans);
	    if (Sskip > nsampsinfile) {fprintf(stderr, "\n\tERROR in quickgplot:\n\t\tSkipping %lld data samples will surpass whole fil file!\n\n",Sskip); exit(7);}
	    else nFBsamps = nsampsinfile - Sskip;
	    if (Sread) {
		if (Sread > nsampsinfile){
		    fprintf(stderr, "\n\t**WARNING:\n\t\tRead samples (%d) > number of samples in file (%d). Setting Sread = %d\n\n",Sread,nsampsinfile,nsampsinfile);
		    Sread = nsampsinfile;
		} else if (Sread + Sskip > nsampsinfile){
		    fprintf(stderr, "\n\t**WARNING:\n\t\tskipped + read samples (%d + %d) > number of samples in file (%d). Setting Sread = %d\n\n",Sskip,Sread,nsampsinfile,nsampsinfile-Sskip);
		    Sread = nsampsinfile-Sskip;
		}
		nFBsamps = Sread;
	    }

	    // 8 throughout is the number of bits per byte; necessary because this is bitwise date
	    rawarchive = (int*) malloc(nifs*nchans*nFBsamps*nbits/8);
	    if (rawarchive == NULL) {fputs ("Memory error",stderr); exit(2);}
	    
	    long long int nfseek = Sskip*nchans*nifs*nbits/8;
	    if (fseek(infile,nfseek,SEEK_CUR)) {fprintf(stderr,"\n\nfseek failed to read %lld bytes from file.\n\n",nfseek);exit(0);}

	    // filesizing is the number of bytes read in from the raw data
	    filesizing = fread(rawarchive,1,nifs*nchans*nbits*nFBsamps/8,infile);
	    if(filesizing!=nFBsamps*nchans*nbits*nifs/8) {fputs("Read error",stderr); exit(3);}
	    fclose(infile);
	    
	    archive = (int*) malloc(filesizing*sizeof(int)*8);
	    if (archive == NULL) {fprintf (stderr,"\nCould not malloc %f bytes for wise data.\n",filesizing*sizeof(int)*8);exit(2);}
	    
	    //---GRAB BITWISE DATA---
	    filesizing = 0;
	    for (i=0;i<nFBsamps*nchans*nifs*nbits/(sizeof(int)*8);i++){
		int abyte = rawarchive[i];
		for (j=0;j<(sizeof(int)*8/nbits);j++){
		    int andvalue = pow(2,(int)nbits)-1;
		    int bitshift = (int)nbits;
		    archive[filesizing++]= andvalue & abyte;
		    abyte = abyte >> bitshift;
		}
	    }
	    floatarchive = (float*) malloc(filesizing*sizeof(float));
	    floatarchive = filint2float(archive,filesizing);
	    free(archive);
	    free(rawarchive);
	} else {fprintf(stderr,"File type %d is an unknown format for filterbank file\n",data_type); exit(6); }
	isokay[16] = true;
    } else {
	fprintf(stderr,"File not found!\n");
	isokay[16] = false;
    }
    nsampsStart=nsampsinfile;
    
//    cout<<filesizing<<" should be "<<nFBsamps*nchans<<"\n";
    
    if (isokay[16]){
	for (int sdecnow=2;sdecnow<=Sdec;sdecnow*=2){
	    filavg(&nFBsamps,nchans,floatarchive);
	}
    }
    
//    cout<<"checkpoint 3\n";

    fchlast = fch1+nchans*foff;

//    cout<<"checkpoint 4\n";

    
//--------------READ PULSE DATA FILE--------------
    int numlines = 0;
    sprintf(currentfile,"%s.PULSEDATA",argv[1]);
    infile = fopen(currentfile, "r");
    if (infile==NULL){
	fprintf(stderr,"***FATAL ERROR: File %s not found!\n",currentfile);
	exit(0);
    } else {
	char textline[300];
	while(fgets(textline, 300, infile)){
	    if(textline[0] == '#' || (textline[0] == 'f' && textline[1] == 'i')) continue;
	    else numlines++;
	}
	rewind(infile);
	dmplot = (float*) malloc (sizeof(float)*numlines);
	snrplot = (float*) malloc (sizeof(float)*numlines);
	beamplot = (int*) malloc (sizeof(int)*numlines);
	tscrplot = (float*) malloc (sizeof(float)*numlines);
	
	i = 0;
	while(fgets(textline, 300, infile)){
	    if(textline[0] == '#' || (textline[0] == 'f' && textline[1] == 'i')) continue;
	    else {
		sscanf(textline,"%*s %*f %f %lld %lld %d %f %f %d",&snrplot[i],&tempstart,&temppeak,&tempwid,&tscrplot[i],&dmplot[i],&beamplot[i]);
		if (tempstart<nsampsStart) nsampsStart = tempstart;
		if (tempstart+tempwid>latestend) latestend = tempstart+tempwid;
//		if (tempwid>widestwid) widestwid = tempwid;
		if (snrplot[i]>maxsnrdetected){
		    maxpeakdetected = temppeak;
		    maxsnrdetected=snrplot[i];
		    bestfitdm = dmplot[i];
		}
	    }
	    i++;
	}
	fclose(infile);
    }

//    cout<<"nsampsstart orig is "<<nsampsStart<<" but subtracting "<<Sskip<<"\n";
    nsampsStart = (nsampsStart-Sskip);
    nsampsStart/=Sdec;
    latestend = (latestend-Sskip);
    latestend/=Sdec;
    maxpeakdetected=(maxpeakdetected-Sskip)/Sdec;//not currently using maxpeakdetected
 
//  cout<<"everything read fine...\n";

//-------If NECESSARY, TIME-AVERAGE THE DATA--------
    int scrunchfacnow = 2;
    while (scrunchfacnow <= Sdec){ 
	for (i=0;i<16;i++){
	    if (isokay[i]){
		timeavg(2*nsamp/scrunchfacnow,timeseries[i]);
	    }
	}
	scrunchfacnow*=2;
    }
    nsamp/=Sdec;
    Sskip/=Sdec;
    tsamp*=Sdec;

//    cout<<"There are "<<nsamp<<" elements in timeseriesindex\n";
//--------NORMALISE DATA AND SET PLOT RANGES--------
    float timeseriesindex[nsamp];
    for (i=0; i<nsamp; i++){
	timeseriesindex[i]=((float)(Sskip*tsamp) + float(i)*tsamp);
    }
    double tempvariable;
    float plottemp[14];
    datamin = 9999;
    datamax = -9999;

    for (i=0;i<16;i++){
	if (isokay[i]){
	    plottemp[i+1] = normalise(nsamp,timeseries[i],&tempvariable);
	    snrmin = filgetmin(timeseries[i],nsamp);
	    snrmax = filgetmax(timeseries[i],nsamp);
	    if (snrmin < datamin) datamin = snrmin;
	    if (snrmax > datamax){
		datamax = snrmax;
	    }
	}
    }

    
    datamin -= 0.1*datamax;
    datamax += 0.1*datamax; 

//    cout<<"checkpoint 5\n";   

    snrmax = getmax(snrplot,numlines,&maxind);
    dmmax = filgetmax(dmplot,numlines);
    tscrmax = filgetmax(tscrplot,numlines);
    snrmax = 1.1*snrmax;            //leave space in plot
    dmmax = 1.1*dmmax;              //leave space in plot
    tscrmax = 1.1*tscrmax;          //leave space in plot
    besttscr = tscrplot[maxind];
    snrmin = filgetmin(snrplot,numlines);
    switch ( beamplot[maxind] ) {
	case 'A' :
	    bestbeam = 9; break;
	case 'B' : 
	    bestbeam = 10; break;
	case 'C' : 
	    bestbeam = 11; break;
	case 'D' : 
	    bestbeam = 12; break;
	default : 
	    bestbeam=beamplot[maxind]-1;
    }
    bestfitdm = dmplot[maxind];
    
//    cout<<"checkpoint 6\n";
//    cout<<"\nnsampsStart "<<nsampsStart<<"\nlatestend "<<latestend<<"\n";
//    cout<<"timseer "<<timeseries[bestbeam][nsampsStart]<<"\n";

    plottemp[0] = getmax(&timeseries[bestbeam][nsampsStart],latestend-nsampsStart,&datamaxindex);

//    cout<<"checkpoint 7\n";
//    cout<<"\ndatamaxindex "<<datamaxindex<<"\n";
//    cout<<"\ndatamaxplotpoint"<<timeseriesindex[datamaxindex+nsampsStart]<<"\n";

    datamaxplotpoint[0] = timeseriesindex[datamaxindex+nsampsStart];//fixed?




//-----------------PLOT EVERYTHING----------------

//    sprintf(currentfile,"/xs");
    sprintf(currentfile,"%s.png/png",argv[1]);
    cpgopen(currentfile);
//    cpgopen("/xs");
    //--- Best DM timeseries ---//
    cpgsvp(0.04,0.47,0.25,0.35);
    cpgswin(timeseriesindex[0],timeseriesindex[nsamp-1],datamin,datamax);
    cpgbox("BCST",0.0,0,"BCST",0.0,0);
    sprintf(currentfile,"%.1f",dmplot[maxind]);
    cpgmtxt("R",1.0,0.5,0.5,currentfile);
    if (isokay[bestbeam]){
	cpgsci(15);
	cpgline(nsamp,timeseriesindex,timeseries[bestbeam]);
	cpgsci(2);
	cpgpt(1,datamaxplotpoint,plottemp,12);
    }
    cpgsci(1);

    //--- One above DM timeseries ---//
    cpgsvp(0.04,0.47,0.35,0.45);
    cpgswin(timeseriesindex[0],timeseriesindex[nsamp-1],datamin,datamax);
    cpgbox("BCST",0.0,0,"BCST",0.0,0);
    cpgmtxt("R",1.0,0.5,0.5,"1up");
    if (isokay[15]){
	cpgsci(15);
	cpgline(nsamp,timeseriesindex,timeseries[15]);
    }
    cpgsci(1);

    //--- One below DM timeseries ---//
    cpgsvp(0.04,0.47,0.15,0.25);
    cpgswin(timeseriesindex[0],timeseriesindex[nsamp-1],datamin,datamax);
    cpgbox("BCST",0.0,0,"BCST",0.0,0);
    cpgmtxt("R",1.0,0.5,0.5,"1dn");
    if (isokay[14]){
	cpgsci(15);
	cpgline(nsamp,timeseriesindex,timeseries[14]);
    }
    cpgsci(1);

    //--- 0.00 DM timeseries ---//
    cpgsvp(0.04,0.47,0.05,0.15);
    cpgswin(timeseriesindex[0],timeseriesindex[nsamp-1],datamin,datamax);
    cpgbox("BCNST",0.0,0,"BCNVST",0.0,0);
    cpgmtxt("R",1.0,0.5,0.5,"0.00");
    if (isokay[13]){
	cpgsci(15);
	cpgline(nsamp,timeseriesindex,timeseries[13]);
    }
    cpgsci(1);

//    cout<<"plotted the timeseries...\n";
    cpgsch(1.0);
    //--- SNR v DM ---//
    cpgsvp(0.55,0.75,0.63,0.93);
    dmmax = 1.1*dmmax;
    cpgswin(0,dmmax,0,snrmax);
    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgmtxt("T",1,0.5,0.5,"SNR vs. DM");
    color=0;
    cpgsci(color);

    int prevbeam = -1;
    float tempfloatdm[1], tempfloatsnr[1];
    i=0;
    while (i<numlines){
	if (tscrplot[i] == besttscr){//!!!!!
	    if (beamplot[i]!=prevbeam){
		cpgsci(++color);
	    }
	    tempfloatdm[0] = dmplot[i];
	    tempfloatsnr[0] = snrplot[i];
	    cpgpt(1,tempfloatdm,tempfloatsnr,3); // 3-->asterisks
	    prevbeam = beamplot[i];
	}//!!!!!
	i++;
    }
    cpgsci(1);

//    cout<<"plotted SNRvDM...\n";
    //--- SNR v TSCR ---//
    cpgsvp(0.75,0.95,0.63,0.93);
    cpgswin(0,tscrmax,0,snrmax);
    cpgbox("BCNST",0.0,0,"BCST",0.0,0);
    cpgmtxt("T",2,0.5,0.5,"SNR vs.");
    cpgmtxt("T",1,0.5,0.5,"Boxcar width");
    color=0;
    cpgsci(color);

    prevbeam = -1;
    i=0;
    while (i<numlines){
	if (dmplot[i] == dmplot[maxind]){//!!!!!
	    if (beamplot[i]!=prevbeam){
		cpgsci(++color);
	    }
	    tempfloatdm[0] = tscrplot[i];
	    tempfloatsnr[0] = snrplot[i];
	    cpgpt(1,tempfloatdm,tempfloatsnr,3); // 3-->asterisks
	    prevbeam = beamplot[i];
	}//!!!!!
	i++;
    }
    cpgsci(1);

    cpgsch(1.0);
    //--- SNR v DM legend ---//
    cpgsvp(0.5,0.99,0.47,0.60);
    cpgswin(0.0,1.0,0.0,1.0);
//    cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
    cpgsci(15);
    cpgtext(0.04,0.85, "POSITIVE DETECTION IN BEAMS:");
    double tempvar = 0.1;
    color=0; prevbeam = -1; i=0;
    cpgsci(color);
    char beamholder[2];
    beamholder[1] = 0;//this little operation cleans out any text after beamholder[0]; to avoid annoying jibberish
    while (i<numlines){
	if (tscrplot[i] == besttscr){//!!!!!
	    if (beamplot[i]!=prevbeam){
		switch ( beamplot[i] ) {
		    case 10 :
			beamholder[0] = 'A'; break;
		    case 11 : 
			beamholder[0] = 'B'; break;
		    case 12 : 
			beamholder[0] = 'C'; break;
		    case 13 : 
			beamholder[0] = 'D'; break;
		    default : 
			sprintf(beamholder,"%d",beamplot[i]);
		}
		cpgsci(++color);
		cpgtext(tempvar,0.65,beamholder);
		tempvar+=0.03;
	    }
	    prevbeam = beamplot[i];
	}//!!!!!
	i++;
    }


    cpgsci(15);
    cpgtext(0.04,0.5, "Parameters at best SNR detection:");
    sprintf(currentfile,"SNR = %.2f",snrplot[maxind]);
    cpgtext(0.1,0.333, currentfile);
    sprintf(currentfile,"DM = %.2f",dmplot[maxind]);
    cpgtext(0.1,0.167, currentfile);
    sprintf(currentfile,"tscr = %.2f",tscrplot[maxind]);
    cpgtext(0.1,0.0, currentfile);
    sprintf(currentfile,"RA = %.2f",src_raj);
    cpgtext(0.48,0.333, currentfile);
    sprintf(currentfile,"Dec = %.2f",src_dej);
    cpgtext(0.48,0.167, currentfile);
    if (gal_b==0) sprintf(currentfile,"glat = not provided");
    else sprintf(currentfile,"glat = %.4f",gal_b);
    cpgtext(0.48,0.0, currentfile);
    cpgsci(1);

    //--- Greyscale freq v. time plot ---//
    // filterbank plot
    int favg = 1;

    setcolortable(colortable);
    if (isokay[16]){
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);

	// kill channels in killfile
	int *killchans = new int[nchans];
	if (dokill){
	    dokill = getkillfile(killchans,nchans,killfile);
	}
	if (dokill){
	    for (int ich=0;ich<nchans;ich++){
                if (killchans[ich] == 0){
                    killchan(floatarchive,nFBsamps,ich,(snrmax+snrmin)/2,nchans);
                }
            }
	}
	// Chan-baseline to clean up the plot
	chanbaseline(floatarchive,nFBsamps,nchans);
	if (nchans>100){
	    if (snrplot[maxind]<18){
		filchanavg(nFBsamps,&nchans,floatarchive);
		favg*=2;
		chanbaseline(floatarchive,nFBsamps,nchans);
		if (snrplot[maxind]<13){
		    filchanavg(nFBsamps,&nchans,floatarchive);
		    favg*=2;
		    chanbaseline(floatarchive,nFBsamps,nchans);
		    if (snrplot[maxind]<9){
			filchanavg(nFBsamps,&nchans,floatarchive);
			favg*=2;
			chanbaseline(floatarchive,nFBsamps,nchans);		    
		    }
		}
	    }
	}
	snrmax = filgetmax(floatarchive,nFBsamps);
	snrmin = filgetmin(floatarchive,nFBsamps);
	cpgsvp(0.55,0.95,0.05,0.45);
	cpgswin(timeseriesindex[0],timeseriesindex[nsamp-1],fch1,fchlast);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
//	cpgscir(2,7);
	float tr[] = {timeseriesindex[0], 0.0, tsamp, fch1, foff*favg, 0.0};
//	float tr[] = {timeseriesindex[0], 0.0, tsamp/Sdec, fch1, foff, 0.0};
	cpgimag(floatarchive,nchans,nFBsamps,1,nchans,1,nFBsamps,snrmin,snrmax,tr);
	cpgslw(3);
	cpgsci(0);
	plotdm(Sskip,0,bestfitdm, nchans, tsamp, foff*favg, fch1);
	cpgslw(1);
	cpgsci(1);
    }
   




//    cout<<"about to plot all the beams...\n";
    //--- And the 13 beams ---//
    cpgsvp(0.205,0.295,0.705,0.795);//1
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[0],"1",bestbeam,datamaxplotpoint,plottemp,isokay[0],nsampsStart,latestend);
    cpgsvp(0.3075,0.3975,0.635,0.725);//2 (from 4)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[1],"2",bestbeam,datamaxplotpoint,plottemp,isokay[1],nsampsStart,latestend);
    cpgsvp(0.205,0.295,0.57,0.66);//3
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[2],"3",bestbeam,datamaxplotpoint,plottemp,isokay[2],nsampsStart,latestend);
    cpgsvp(0.1025,0.1925,0.635,0.725);//4 (from 2)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[3],"4",bestbeam,datamaxplotpoint,plottemp,isokay[3],nsampsStart,latestend);
    cpgsvp(0.1025,0.1925,0.775,0.865);//5 (from 7)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[4],"5",bestbeam,datamaxplotpoint,plottemp,isokay[4],nsampsStart,latestend);
    cpgsvp(0.205,0.295,0.84,0.93);//6
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[5],"6",bestbeam,datamaxplotpoint,plottemp,isokay[5],nsampsStart,latestend);
    cpgsvp(0.3075,0.3975,0.775,0.865);//7 (from 5)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[6],"7",bestbeam,datamaxplotpoint,plottemp,isokay[6],nsampsStart,latestend);
    cpgsvp(0.41,0.5,0.705,0.795);//8 (from B)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[7],"8",bestbeam,datamaxplotpoint,plottemp,isokay[7],nsampsStart,latestend);
    cpgsvp(0.3075,0.3975,0.5,0.59);//9 (from A)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[8],"9",bestbeam,datamaxplotpoint,plottemp,isokay[8],nsampsStart,latestend);
    cpgsvp(0.1025,0.1925,0.5,0.59);//A (from 9)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[9],"A",bestbeam,datamaxplotpoint,plottemp,isokay[9],nsampsStart,latestend);
    cpgsvp(0.0,0.09,0.705,0.795);//B (from 8)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[10],"B",bestbeam,datamaxplotpoint,plottemp,isokay[10],nsampsStart,latestend);
    cpgsvp(0.1025,0.1925,0.91,1.0);//C (from D)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[11],"C",bestbeam,datamaxplotpoint,plottemp,isokay[11],nsampsStart,latestend);
    cpgsvp(0.3075,0.3975,0.91,1.0);//D (from C)
    printbeam(nsamp,datamin,datamax,timeseriesindex,timeseries[12],"D",bestbeam,datamaxplotpoint,plottemp,isokay[12],nsampsStart,latestend);


    cpgclos();

    return(0);  
}






//---SUBFUNCTIONS---------------------------------------------


void printbeam(int n, float min, float max, float *xdata, float *ydata, const char *beamid, int bestbeam, float *maxindex, float *maxdata, bool dataok, int start, int stop){
    int bestbeamref=-1;
    cpgswin(xdata[0],xdata[n-1],min,max);
    cpgbox("BC",0.0,0,"BC",0.0,0);
    cpgmtxt("B",1,0.5,0.5,beamid);
    if (dataok){
	cpgsci(15);
	cpgline(n,xdata,ydata);
	switch ( *beamid ) {
	    case 'A' :
		bestbeamref = 9; break;
	    case 'B' : 
		bestbeamref = 10; break;
	    case 'C' : 
		bestbeamref = 11; break;
	    case 'D' : 
		bestbeamref = 12; break;
	    default : 
		bestbeamref = *beamid-'0';
		bestbeamref-=1;
	}
	if (bestbeamref == bestbeam){
	    cpgsci(2);
	    cpgpt(1,maxindex,maxdata,12);
	} else {//?????
	    cpgsci(2);//?????
	    cpgsls(4);//?????
	    float detlinex[2] = {maxindex[0],maxindex[0]};//?????
	    float detliney[2] = {min,max};//?????
	    cpgline(2,detlinex,detliney);//?????
	    detlinex[0] = xdata[start];//?????
	    detlinex[1] = xdata[stop];//?????
	    detliney[0] = maxdata[bestbeamref+1];//?????
	    detliney[1] = maxdata[bestbeamref+1];//?????
	    cpgsls(1);
	    cpgline(2,detlinex,detliney);//?????
	    cpgsls(1);
	}
    }
    cpgsci(1);
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
    fprintf(stderr,"\t-s   number of samples to skip [0]\n");
    fprintf(stderr,"\t-r   number samples to read [all]\n");
    fprintf(stderr,"\t-dec factor to scrunch data by [1, no scrunching]\n");
    fprintf(stderr,"\t-c   color table for plotting filterbank [def-heat]:\n");
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
