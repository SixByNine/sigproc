#include "gtools.h"
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

// This is a comment

// Gpulse CLASS "METHODS"
Gpulse::Gpulse() {
}
Gpulse::Gpulse(float a, float s, int b, int l, int w, int t, float d) {
    amp = a;
    SNR = s;
    loc = l;
    width = w;
    start = b;
    tscrfac = t;
    dm = d;
    beam=-1; //if beam is -1, the beam ID has not yet been initiated
}
Gpulse::Gpulse(float a, float s, int b, int l, int w, int t, float d, int bm) {
    amp = a;
    SNR = s;
    loc = l;
    width = w;
    start = b;
    tscrfac = t;
    dm = d;
    beam = bm;
}
void Gpulse::put_pulse(float a, float s, int b, int l, int w, int t, float d) {
	amp = a;
	SNR = s;
	loc = l;
	width = w;
	start = b;
	tscrfac = t;
	dm = d;
}
void Gpulse::put_pulse(float a, float s, int b, int l, int w, int t, float d, int bm) {
	amp = a;
	SNR = s;
	loc = l;
	width = w;
	start = b;
	tscrfac = t;
	dm = d;
	beam = bm;
}

//**************************************************
//**************************************************
//***               SEARCH CODES                 ***
//**************************************************
//**************************************************

//State: to hold results of multiple DM trials
GPulseState::GPulseState(int ndms) {
    DMtrials = new vector<Gpulse>[ndms];
    NDMtrials = ndms;
}

//Searchforgiants: Runs findgiants search for a trial DM of a GPulseState's DM range.
//                 Also corrects for sampling index offsets caused by gulping
void GPulseState::searchforgiants(int itrial, int numbersamples, int offset,  float * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac) {
    DMtrials[itrial] = findgiants(numbersamples, data, nsigma, bintol, usertscrfac, dm,starttscrfac);
    if (offset>0){
	for (int i=0;i<DMtrials[itrial].size();i++){
	    DMtrials[itrial].at(i).start += offset;
	    DMtrials[itrial].at(i).loc += offset;
	}
    }
}
void GPulseState::searchforgiants(int itrial, int numbersamples, int offset, unsigned short int * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac) {
    DMtrials[itrial] = findgiants(numbersamples, data, nsigma, bintol, usertscrfac, dm,starttscrfac);
    if (offset>0){
	for (int i=0;i<DMtrials[itrial].size();i++){
	    DMtrials[itrial].at(i).start += offset;
	    DMtrials[itrial].at(i).loc += offset;
	}
    }
}

void GPulseState::selfdestruct(){
    for (int i=0;i<NDMtrials;i++){
	DMtrials[i].erase(DMtrials[i].begin(),DMtrials[i].end());
    }
}

//Givetimes: Does multi-DM coincidence matching, returns candidate pulse sample stamps.
//Results are written to GResults.txt unless otherwise specified. Results can be beam-tracked if desired.
//The BEAM ID should ONLY be specified if wanting to do A MULTIBEAM SEARCH (otherwise beamID=-1)
int* GPulseState::givetimes(int* ndetected, float sampletime, float flo,float fhi,float irrel,char* filetimestamp){
    return(givetimes(ndetected,sampletime,flo,fhi,irrel,filetimestamp,-1,"GResults.txt"));
}
int* GPulseState::givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, char* filetimestamp, int beamID){
    return(givetimes(ndetected, sampletime, flo, fhi, irrel,filetimestamp,beamID,"GResults.txt"));
}
int* GPulseState::givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, char* filetimestamp, int beamID, char* resultsfilename) {
    int nsinglebeamcands, delayinsamples, totdelay;
    float delayinms;
    vector<Gpulse> suspectvectorstorage;
    vector<Gpulse> SPvectorstorage;
    Gpulse gpulsestorage;


    FILE* SPfile;
    char SPfilename[200];
    
    for (int i=0; i<NDMtrials; i++) {
	suspectvectorstorage.insert(suspectvectorstorage.end(),DMtrials[i].begin(), DMtrials[i].end());
     }

    fprintf(stderr,"\n\nN candidates in this block before associating: %d\n",suspectvectorstorage.size());

    vector<Gpulse>* suspectarraystorage = new vector<Gpulse>[suspectvectorstorage.size()];
    suspectarraystorage = assoc_giants(suspectvectorstorage,&nsinglebeamcands,resultsfilename,filetimestamp,beamID,irrel);


    fprintf(stderr,"\n\nN candidates in this block before associating: %d\n",suspectvectorstorage.size());
    fprintf(stderr,"N candidates in this block after associating: %d\n\n",nsinglebeamcands);
    int* detectiontimestamps = new int[nsinglebeamcands*2];
//    if (beamID<0){ BEAMID
    for (int i=0; i<(nsinglebeamcands*2); i+=2){
	gpulsestorage = suspectarraystorage[i/2].at(0);
	SPvectorstorage = suspectarraystorage[i/2];
	if (flo<fhi)
	    delayinms = gpulsestorage.dm * 4.15 * (pow(flo/1000, -2) - pow(fhi/1000, -2));
	if (fhi<flo)
	    delayinms = gpulsestorage.dm * 4.15 * (pow(fhi/1000, -2) - pow(flo/1000, -2));
	delayinsamples = (int)(delayinms/(sampletime*1000))+1;
	fprintf(stderr,"Candidate %4d: DM %5.2f SNR %5.2f SCR %d",i/2,gpulsestorage.dm,gpulsestorage.SNR,gpulsestorage.tscrfac);//2009-04-30-02:57:52.fil
	fprintf(stderr,";  there were %d detections of this candidate.\n",SPvectorstorage.size());
	
//	    SPfilename contains a .pulse file for EACH CANDIDATE.
	sprintf(SPfilename,"%f_%f_%d_%d_%d_%d_%07.2f_%s_%d.pulse",gpulsestorage.amp,gpulsestorage.SNR,gpulsestorage.start,gpulsestorage.loc,gpulsestorage.width,gpulsestorage.tscrfac,gpulsestorage.dm,filetimestamp,beamID);
	SPfile = fopen(SPfilename,"w");
	for (int j=1;j<SPvectorstorage.size();j++){ 
//                                   filename  amp   snr  start peak wid tscr dm   beam
//		fprintf(stderr,"filenamedummy %8.5f %8.5f %12d %12d %12d %6d %8.2f %d\n",SPvectorstorage.at(j).amp,SPvectorstorage.at(j).SNR,SPvectorstorage.at(j).start,SPvectorstorage.at(j).loc,SPvectorstorage.at(j).width,SPvectorstorage.at(j).tscrfac,SPvectorstorage.at(j).dm,beamID);
	    fprintf(SPfile,"%s\t%8.5f %8.5f %12d %12d %12d %6d %8.2f %d\n",filetimestamp,SPvectorstorage.at(j).amp,SPvectorstorage.at(j).SNR,SPvectorstorage.at(j).start,SPvectorstorage.at(j).loc,SPvectorstorage.at(j).width,SPvectorstorage.at(j).tscrfac,SPvectorstorage.at(j).dm,beamID);
	}
	fclose(SPfile);
//          resultsfile is a SUMMARY FILE for all the detected candidates.
//	fprintf(resultsfile,"Candidate %4d: DM %5.2f SNR %5.2f SCR %4d STARTBIN %13d PEAK %13d\n",i/2,gpulsestorage.dm,gpulsestorage.SNR,gpulsestorage.tscrfac,gpulsestorage.start,gpulsestorage.loc);
	totdelay = delayinsamples+gpulsestorage.width;
	if (gpulsestorage.start-(totdelay)<0)
	    detectiontimestamps[i] = 0;
	else
	    detectiontimestamps[i] = gpulsestorage.start-(totdelay);
	detectiontimestamps[i+1] = 3*(totdelay);
    }
    *ndetected = nsinglebeamcands;
    suspectvectorstorage.erase(suspectvectorstorage.begin(),suspectvectorstorage.end());
    SPvectorstorage.erase(SPvectorstorage.begin(),SPvectorstorage.end()); //added to try to abate a possible memory leak
    
    return (detectiontimestamps);
//    } else { BEAMID
    //do something BEAMID
//    } BEAMID
}

/* WHAT ABOUT MULTIBEAM?!?!?
 if (multibeam){
 //	suspectvectorstorage.clear();
 //	for (int i=0;i<nsinglebeamcands-1;i++){
 //	    suspectvectorstorage.push_back(suspectarraystorage[i].at(0));
 //	}
 candidates = beamassoc_giants(suspectarraystorage, 13, 5.0);
 //	beamassoc_giants();
 } else {
 int whatever = 6;
 //return assoc_giants results
 }
 */




//**************************************************
//**************************************************
//***             GSEARCHING LOOPS               ***
//**************************************************
//**************************************************


//*******************************************
// WRAPPERS TO SET DATA TYPES, TIME SCRUNCH,
//  & NORMALIZATION FOR TIMESERIES SEARCH
//*******************************************
// for 32-bit (float) data
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol, int usertscrfac, float dm, int starttscrfac) {
	vector<Gpulse> mybigvector, giant;
	double mean=0;
	double sigma=0;
	int newnpts;
	for (int tscrfac=starttscrfac; tscrfac<=usertscrfac; tscrfac*=2) {
		newnpts = npts*starttscrfac/tscrfac;

                // First pass: subtract mean and get stdev
		if (tscrfac==starttscrfac){
		    sigma = submeanrms(newnpts, data, &mean);
		    sigma = getmowedsigma(newnpts,data,sigma,0); //USELESS IF DATA NOT BASELINED FIRST.
		} else {
		    sigma = timeavg(2*newnpts, data, &mean);
		    sigma = getmowedsigma(newnpts,data,sigma,0);
		}
//		cout<<"sigma is "<<sigma<<" at tscrfac "<<tscrfac<<"\n";
		giant = giantsearch(newnpts, data, nsigma*sigma, sigma, bintol/tscrfac, tscrfac, dm);
		mybigvector.insert(mybigvector.end(), giant.begin(), giant.end());
	}
	return (mybigvector);
}

//shortcut (start scr factor@1)
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol, int usertscrfac, float dm) {
	return findgiants(npts, data, nsigma, bintol, usertscrfac, dm, 1);
}

// for unsigned short int data (e.g. HITRUN survey dedisperser)
vector<Gpulse> findgiants(int npts, unsigned short int * data, float nsigma, float bintol,
			  int usertscrfac, float dm,int starttscrfac) {
	vector<Gpulse> mybigvector;
	float *convertedarray;
	convertedarray = (float *) malloc(npts*sizeof(float));
	for (int j=0; j<npts; j++) {
		convertedarray[j]=(float)data[j];
	}
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,dm);
	free(convertedarray);
	return (mybigvector);
}

// for double (64-bit) data 
vector<Gpulse> findgiants(int npts, double * data, float nsigma, float bintol,
		int usertscrfac, float dm) {
	vector<Gpulse> mybigvector;
	float *convertedarray;
	convertedarray = (float *) malloc(npts*sizeof(float));
	for (int j=0; j<npts; j++) {
		convertedarray[j]=(float)data[j];
	}
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,dm);
	free(convertedarray);
	return (mybigvector);
}

// for int data 
vector<Gpulse> findgiants(int npts, int * data, float nsigma, float bintol,
		int usertscrfac, float dm) {
	vector<Gpulse> mybigvector;
	float *convertedarray;
	convertedarray = (float *) malloc(npts*sizeof(float));
	for (int j=0; j<npts; j++) {
		convertedarray[j]=(float)data[j];
	}
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,dm);
	free(convertedarray);
	return (mybigvector);
}

// for 8-bit (byte) data
vector<Gpulse> findgiants(int npts, unsigned char * data, float nsigma,
		float bintol, int usertscrfac, float dm) {
	vector<Gpulse> mybigvector;
	float *convertedarray;
	convertedarray = (float *) malloc(npts*sizeof(float));
	for (int j=0; j<npts; j++) {
		convertedarray[j]=(float)data[j];
	}
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,dm);
	free(convertedarray);
	return (mybigvector);
}


//*******************************************
//        ACTUAL TIMESERIES SEARCH
//*******************************************
vector<Gpulse> giantsearch(int n, float * data, float thresh, double RMS, float bintol, int tscrfac, float dm) {
	int i;
	Gpulse pulse;
	bool detected = false, first = true;
	int xstart = 0, width = 0; //first bin of detected pulse
	int j, imax=0, ngiants=0, lastend=0; //counters
	float snr; //imax will store the location in d of the max peak of a pulse
	vector<Gpulse> giants;

	for (i=0; i<n; i++) {
		if (data[i]>thresh && detected==false)//when we get above the threshold, note where it began.
		{
			detected = true;
			if (first==false) {
				if ((i-lastend) > bintol) {
					pulse.put_pulse(data[imax], data[imax]/RMS, xstart*tscrfac,
							imax*tscrfac, (lastend-xstart)*tscrfac, tscrfac, dm);
					giants.push_back(pulse);
					xstart = i;//if the last pulse ended more than bintol bins ago, count this as a new pulse.
				}
			} else {
				first = false;
				xstart = i;
			}
		} else if (data[i]<thresh && detected==true) //when we go under the threshold, search the in-to-out width for the maximum.
		{
			imax = xstart;
			lastend = i;
			for (j=xstart; j<i; j++)
				if (data[imax]<data[j])
					imax=j;
			detected = false;
		} else if (i==n-1) {
			if (data[i]>thresh && detected==true) {
				pulse.put_pulse(data[imax], data[imax]/RMS, xstart*tscrfac,
						imax*tscrfac, (i-xstart)*tscrfac, tscrfac, dm);
				giants.push_back(pulse);
			} else if (first==false && data[i]<thresh) {
				pulse.put_pulse(data[imax], data[imax]/RMS, xstart*tscrfac,
						imax*tscrfac, (lastend-xstart)*tscrfac, tscrfac, dm);
				giants.push_back(pulse);
			}
		}
	}
	return (giants);
}



//**************************************************
//**************************************************
//***          OTHER USEFUL FUNCTIONS            ***
//**************************************************
//**************************************************

//**************************************************
//             NORMALIZATION ROUTINES
//**************************************************
// straight normalization
void normalise(float *data, int n) {
	double sum=0.0, sumsq=0.0;
	int i=0;
	while (i<n) {
		sum += data[i];
		sumsq += (data[i] * data[i]);
		i++;
	}
	double mean=sum/(double)n;
	double meansquares=sumsq/(double)n;
	double sigma= sqrt(meansquares - (mean * mean));
#pragma omp parallel for private(i)
	for (int i=0; i<n; i++)
		data[i]=(data[i]-mean)/sigma;
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



//**************************************************
// subtract mean from data, return 1sigma value
//**************************************************

double submeanrms(int n, float * d, double * dataaverage) {
    double sum=0.0;
    double sumsq=0.0;
    for (int i=0; i<n; i++) {
	sum+=d[i];
	sumsq+=d[i]*d[i];
    }
    double mean=sum/(double)n;
    *dataaverage = mean;

    double rmssquared=sumsq/(double)n;
    double sigma=sqrt(rmssquared-mean*mean);
    if (rmssquared==0) {
	fprintf(stderr, "submeanrms::RMS of data is zero\n");
	return (sigma);
    }
    int i;
#pragma omp parallel for private(i)
    for (i=0; i<n; i++)
	d[i]=(d[i]-mean);
    return (sigma);
}

double getrms(int n, unsigned short int * d, double * dataaverage) {
    double sum=0.0;
    double sumsq=0.0;
    for (int i=0; i<n; i++) {
	sum+=(float)d[i];
	sumsq+=(float)d[i]*(float)d[i];
    }
    double mean = sum/(double)n;
    *dataaverage = mean;

    double rmssquared=sumsq/(double)n;
    double sigma=sqrt(rmssquared-mean*mean);
    if (rmssquared==0) {
	fprintf(stderr, "getrms::RMS of data is zero\n");
    }
    return (sigma);
}

double getrms(int n, float * d, double * dataaverage) {
    double sum=0.0;
    double sumsq=0.0;
    for (int i=0; i<n; i++) {
	sum+=(float)d[i];
	sumsq+=(float)d[i]*(float)d[i];
    }
    double mean = sum/(double)n;
    *dataaverage = mean;

    double rmssquared=sumsq/(double)n;
    double sigma=sqrt(rmssquared-mean*mean);
    if (rmssquared==0) {
	fprintf(stderr, "getrms::RMS of data is zero\n");
    }
    return (sigma);
}

// Recalculate sigma without SNR>6 samples.
double getmowedsigma(int n, float * d, double unmowedsigma, double mean) {
    double sum=0.0;
    double sumsq=0.0;
    double sigma;
    for (int i=0; i<n; i++) {
	if (d[i]<(3*unmowedsigma)){
	    sum+=d[i];
	    sumsq+=d[i]*d[i];
	} else {
	    sum += mean;
	    sumsq+=mean*mean;
	}
    }
    sigma=sqrt((sumsq/(double)n) - (mean*mean));
    if (sigma==0) {
	fprintf(stderr, "getmowedsigma::RMS of data is zero\n");
    }
    return(sigma);
}

//**************************************************
//            TIME-SCRUNCHING ROUTINES
//**************************************************

// Time-avg by two & return SIGMA
double timeavg(int n, float * d, double *dataaverage) {
	int asize;
	float avg;
	double sum=0.0;
	double sumsq=0.0;
	for (int i=0; i<n-1; i+=2) {
		avg = (d[i] + d[i+1]) / 2;
		asize = i/2;
		d[asize] = avg;
		sum+=avg;
		sumsq+=avg*avg;
	}
	if ((n-1)%2 == 0) {
		asize = (n-1)/2;
		d[asize] = d[n-1];
	}
	double mean=sum/(double)asize;
	*dataaverage = mean;
	double rmssquared=sumsq/(double)asize;
	double sigma=sqrt(rmssquared-mean*mean);
	if (rmssquared==0) {
		fprintf(stderr, "timeavg::RMS of data is zero\n");
		return (sigma);
	}
	return (sigma);
}

// Time-avg by factor two, no sigma-return
void timeavg(int n, float * d) {
	int asize;
	float avg;
	for (int i=0; i<n-1; i+=2) {
		avg = (d[i] + d[i+1]) / 2;
		asize = i/2;
		d[asize] = avg;
	}
	if (n-1%2==0) {
		asize = n/2;
		d[asize] = d[n];
	}
}


//**************************************************
//      MIN/MAX  +  TYPE CONVERSION ROUTINES
//**************************************************
float getmax(float *data, int arraysize) {
	float max=0;
	int i;
	for (i=0; i<arraysize; i++) {
		if (data[i] > max) {
			max = data[i];
		}
	}
	return (max);
}

float getmax(float *data, int arraysize, int *index) {
	float max=0;
	int i;
	for (i=0; i<arraysize; i++) {
		if (data[i] > max) {
			max = data[i];
			*index = i;
		}
	}
	return (max);
}

float getmin(float *data, int arraysize) {
	float min=0;
	int i;
	for (i=0; i<arraysize; i++) {
		if (data[i] < min)
			min = data[i];
	}
	return (min);
}

float* int2float(int *array, int arraysize) {
	float *floatarray = (float*) malloc(sizeof(float)*arraysize);
	for (int i=0; i<arraysize; i++) {
		floatarray[i] = float(array[i]);
	}
	return (floatarray);
}




//**************************************************
//**************************************************
//***  BASELINE REMOVAL & ADVANCED RFI EXCISION  ***
//**************************************************
//**************************************************

//**************************************************
//      FIRST-PASS RFI EXCISION: assoc_giants
//**************************************************
vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands, char* resultsfilename,char* filetimestamp, int beamID) {
        vector<Gpulse> *cands = new vector<Gpulse>[uapulses.size()/2];
	cands = assoc_giants(uapulses, nsinglebeamcands, resultsfilename,filetimestamp,beamID,5.0); //DEFAULT IRREL DM = 5.0
	return (cands); //IS nsinglebeamcands IMPLEMENTED CORRECTLY HERE?
}

vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands,char* resultsfilename,char* filetimestamp,int beamID,float irrel) {
	int maxdm=0, scanrangelo, scanrangehi, ncandidates=0;
	int npulses = uapulses.size(),suslo,sushi,ninonecand,maxsnrindex=1;
	vector<Gpulse> *candidates;
	vector<Gpulse> onecand;
	FILE* resultsfile = fopen(resultsfilename,"a");
	Gpulse placeholder;
	bool newpulse;
	float maxsnr;
	char dettype;
	cout<<"Associating...\n";
	if (npulses == 0) {
		fprintf(stderr,"\nDid not find any giant pulses in this block.\n");
		*nsinglebeamcands=0;
		return (candidates);
	}

	candidates = new vector<Gpulse>[uapulses.size()/2];
	while (uapulses.size() > 0) {
		//	if (uapulses[0].dm == uapulses[npulses].dm){
		//	    fprintf(stderr,"Only found lone high-DM spikes.\n",irrel);
		//	    return(candidates);
		//	}
		placeholder.put_pulse(0, 0, 0, 0, 0, 0, 0);
		onecand.push_back(placeholder); //first element holds "best of" info
//		printf("\nSussing a pulse at DM %f\n",placeholder.dm);
		onecand.push_back(uapulses[0]); //put on first pulse from unassociated ones
		uapulses.erase(uapulses.begin()); //remove first element
		npulses--;

		scanrangelo = onecand[1].start;
		scanrangehi = onecand[1].start+onecand[1].width;
		maxsnr = onecand[1].SNR;
		maxsnrindex = 1;
//		printf("\n\n\nmaxsnr at first is %f\n\n\n",maxsnr);
		//printf("*start* %d %f %f\n",(int)onecand.size(), onecand[1].dm,onecand[1].SNR);

		newpulse = true;
		ninonecand=1; //begins with pulse already pushed on
		while (newpulse) {
//		    printf("-----------------Scanning range %d to %d---------------\n",scanrangelo,scanrangehi);
		    newpulse = false;
		    for (int i=1; i<npulses; i++) {
			suslo = uapulses[i].start;
			sushi = uapulses[i].start+uapulses[i].width;
			if ((scanrangelo>=suslo&&scanrangelo<=sushi)||(scanrangehi>=suslo&&scanrangehi<=sushi)||(suslo>=scanrangelo&&suslo<=scanrangehi)){
//			    printf("Found associated pulse: DM %f, SNR %f, loc %i, start %i\n",uapulses[i].dm,uapulses[i].SNR,uapulses[i].loc,uapulses[i].start);
			    if (suslo < scanrangelo) {
				scanrangelo = suslo;
			    	newpulse = true;
			    }
			    if (sushi>scanrangehi) {
			    	scanrangehi = sushi;
			    	newpulse = true;
			    }
			    if (uapulses[i].SNR>maxsnr) {
//				printf("newp(%f) is more than maxsnr(%f)\n",uapulses[i].SNR,maxsnr);
			    	maxsnr=uapulses[i].SNR;
			    	maxsnrindex=ninonecand+1;
//				printf("new index for SNR %f is %d\n",maxsnr,maxsnrindex);
			    }
			    onecand.push_back(uapulses[i]);
//			    fprintf(stderr,"loc: %d start: %d wid: %d\n",uapulses[i].loc,uapulses[i].start,uapulses[i].width);
			    uapulses.erase(uapulses.begin()+i); //remove associated element
//			    fprintf(stderr,"i,npulses,uap.size,ninonecand: %d,%d,%d,%d\n",i,npulses,uapulses.size(),ninonecand);
			    onecand[0].width = scanrangehi-scanrangelo;
			    ninonecand++;
			    npulses--;
			    i--;
//			    fprintf(stderr,"cand loc: %d start: %d wid: %d\n\n",onecand[ninonecand].loc,onecand[ninonecand].start,onecand[ninonecand].width);
			}
		    }
		}

		if (onecand.size()<3 || onecand[maxsnrindex].dm<irrel) {
		    printf("The burst is RFI. Npulses:%d DM:%f SNR:%f Loc'n:%d Wid:%d Scr:%d\n\n",onecand.size(),onecand[maxsnrindex].dm,onecand[maxsnrindex].SNR,onecand[maxsnrindex].loc,onecand[0].width,onecand[maxsnrindex].tscrfac);
		    if (onecand.size()<3){
			dettype = 'g';
			fprintf(resultsfile,"%s\t%8.5f %8.5f %12d %12d %12d %6d %8.2f %d %c %d\n",filetimestamp,onecand[maxsnrindex].amp,onecand[maxsnrindex].SNR,scanrangelo,onecand[maxsnrindex].loc,scanrangehi-scanrangelo,onecand[maxsnrindex].tscrfac,onecand[maxsnrindex].dm,beamID,'G',onecand[maxsnrindex].loc);
		    } else if (onecand[maxsnrindex].dm<irrel){
			//nrfi++;
			//detnum = nrfi; //!!! NEED TO MAKE NRFI/DETNUM/DETTYPE
			dettype = 'r';
			fprintf(resultsfile,"%s\t%8.5f %8.5f %12d %12d %12d %6d %8.2f %d %c %d\n",filetimestamp,onecand[maxsnrindex].amp,onecand[maxsnrindex].SNR,scanrangelo,onecand[maxsnrindex].loc,scanrangehi-scanrangelo,onecand[maxsnrindex].tscrfac,onecand[maxsnrindex].dm,beamID,'R',onecand[maxsnrindex].loc);
		    }
		    for (int j=0;j<onecand.size();j++){
			fprintf(resultsfile,"%s\t%8.5f %8.5f %12d %12d %12d %6d %8.2f %d %c %d\n",filetimestamp,onecand.at(j).amp,onecand.at(j).SNR,onecand.at(j).start,onecand.at(j).loc,onecand.at(j).width,onecand.at(j).tscrfac,onecand.at(j).dm,beamID,dettype,onecand[maxsnrindex].loc);
		    }
		    //RFIpulses.insert(RFIpulses.end(), onecand.begin()+1, onecand.end());
		} else {
		    onecand[0].put_pulse(onecand[maxsnrindex].amp,
					 onecand[maxsnrindex].SNR, scanrangelo,
					 onecand[maxsnrindex].loc, scanrangehi-scanrangelo,
					 onecand[maxsnrindex].tscrfac, onecand[maxsnrindex].dm);
		    candidates[ncandidates] = onecand;
//		    printf("GOOD! %4d %10d %5.2f %5.3f\n",onecand.size(),onecand[0].start,onecand[maxsnrindex].dm,onecand[maxsnrindex].SNR);
		    fprintf(resultsfile,"%s\t%8.5f %8.5f %12d %12d %12d %6d %8.2f %d %c %d\n",filetimestamp,onecand[maxsnrindex].amp,onecand[maxsnrindex].SNR,scanrangelo,onecand[maxsnrindex].loc,scanrangehi-scanrangelo,onecand[maxsnrindex].tscrfac,onecand[maxsnrindex].dm,beamID,'C',onecand[maxsnrindex].loc);
		    for (int j=0;j<onecand.size();j++){
			fprintf(resultsfile,"%s\t%8.5f %8.5f %12d %12d %12d %6d %8.2f %d %c %d\n",filetimestamp,onecand.at(j).amp,onecand.at(j).SNR,onecand.at(j).start,onecand.at(j).loc,onecand.at(j).width,onecand.at(j).tscrfac,onecand.at(j).dm,beamID,'c',onecand[maxsnrindex].loc);
		    }
		    cout<<"...done printing.\n";
		    ncandidates++;
		}
		onecand.erase(onecand.begin(),onecand.end());
	}
	fclose(resultsfile);

	*nsinglebeamcands = ncandidates;
	return (candidates);
}






//**************************************************
//      REMOVE BASELINE FROM DATA BEFORE SEARCH
//**************************************************
// This task to be run BEFORE call to searchforgiants or findgiants.
// Adapted from Matthew's find_baseline.C
void removebaseline(unsigned short int *indata, float* outdata, int ndat, int runmeansize, float thresh){
//NOTE    runmeansize is the size of the smoothing filter in number of samples (i.e. tsmooth/tsamp)

    float * tdata = new float[runmeansize];
    if (tdata == NULL) {
	fprintf(stderr,
		"find_baseline error: Error allocating %d floats for tdata\n",runmeansize);
	exit(-3);
    }
//    float * outdata = new float[ndat];
//    if (outdata == NULL) {
//	fprintf(stderr,
//		"find_baseline error: Error allocating %d floats for outdata\n",ndat);
//	exit(-3);
//    }


    // Copy runmeansize indata's into tdata
    // and indata into outdata to preserve indata
    // Also take STDev & mean of all the data
    double sum=0.0;
    double sumsq=0.0;
    for (int i=0;i<ndat;i++){
	if (i<runmeansize)
	    tdata[i]=(float)indata[i];
	outdata[i]=(float)indata[i];
	sum+=outdata[i];
	sumsq+=outdata[i]*outdata[i];
    }
    float mean_f = (float)(sum/(double)ndat);
    float rms_f=sqrt((sumsq/(double)ndat)-mean_f*mean_f);


    if (rms_f==0.0){
	fprintf(stderr,"find_baseline error: rms of data is zero\n");
	exit(-4);
    }
    float inverse_rms_f = 1.0/rms_f;
    
    // Subtract mean to half of smoothing time
    for (int j=0;j<(runmeansize/2);j++) outdata[j]=(outdata[j]-mean_f)*inverse_rms_f;
    //------dorunning mean (eliminating high points)
//  fprintf(stderr,"down here\n");
    
    int na=runmeansize/2;
    int nb=ndat-(runmeansize/2);
    int ja=0;
    float mean = mean_f;
    float div = 1.0 / runmeansize;
    float thresh_r = thresh * rms_f;
    
    for (int j=na;j<nb;j++){
	outdata[j]=outdata[j]-mean;            // Subtract mean
	if (ja == runmeansize-1) ja = 0;
	float al=tdata[ja];
	float an=outdata[j+na];
	if (fabs(an-mean) > thresh_r) an=mean; 
	tdata[ja]=an;
	mean=mean+(an-al)*div;
	outdata[j]=outdata[j]*inverse_rms_f;
	ja++;
    }
    // do end bit
    for (int j=nb;j<ndat;j++)outdata[j]=(outdata[j]-mean)*inverse_rms_f;
    
    delete[] tdata;
}




//**************************************************
//     second-pass RFI excision; beamassocgiants
//**************************************************

/*
THE NODE 14 SOFTWARE WILL DO THE FOLLOWING:

   - Talk to the individual 13 nodes via sockets. It will communicate with MY
     GSEARCH CODE or MIKE'S PROGRAM (inst.of calling for results, send call to
     socket connect). The benefit of communication with my SP GSEARCH program
     is that I can continue to use vectors instead of returning just a
     timestamp

   - Convert input from SOCKETS to something runnable with beam_assoc_giants,
     containing also beam info.

   - Communicate results to Willem's 13 x NODE14EAR sockets

   - Save stats of the pointing to disk as text

   - Save candidate database to disk
 */


/*vector<Gpulse>* beamassoc_giants(vector<Gpulse> *beams, int nbeams, int *nmultibeamcands) {
	vector<Gpulse> *cands;
	cands = beamassoc_giants(beams, nbeams, nmultibeamcands, 5.0); //DEFAULT IRREL DM = 5.0
	return (cands);
}

vector<Gpulse>* beamassoc_giants(vector<Gpulse> *beams, int nbeams, int *nmultibeamcands, float irrel) {
	int maxdm=0, scanrangelo, scanrangehi, ncandidates=0, ninonecand;
	int npulses = uapulses.size(),suslo,sushi,maxsnrindex=1;
	int bdetect[nbeams] = {0}; //which beams the pulse is detected in
	vector<Gpulse> *candidates;
	vector<Gpulse> onecand, RFIpulses;
	Gpulse placeholder;
	bool newpulse;
	float maxsnr;

	if (npulses == 0) {
		fprintf(stderr,"\nDid not find any giant pulses in this block.\n");
		*nmultibeamcands=0;
		return (candidates);
	}

	candidates = new vector<Gpulse>[uapulses.size()/2];
	while (uapulses.size() > 0) {
		placeholder.put_pulse(0, 0, 0, 0, 0, 0, 0);
		onecand.push_back(placeholder); //first element holds "best of" info
//		printf("\nSussing a pulse at DM %f\n",placeholder.dm);
		onecand.push_back(uapulses[0]); //put on first pulse from unassociated ones
		uapulses.erase(uapulses.begin()); //remove first element
		npulses--;

		scanrangelo = onecand[1].start;
		scanrangehi = onecand[1].start+onecand[1].width;
		maxsnr = onecand[1].SNR;
		maxsnrindex = 1;
//		printf("\n\n\nmaxsnr at first is %f\n\n\n",maxsnr);
		//printf("*start* %d %f %f\n",(int)onecand.size(), onecand[1].dm,onecand[1].SNR);

		newpulse = true;
		ninonecand=1; //begins with pulse already pushed on
		while (newpulse) {
//		    printf("-----------------Scanning range %d to %d---------------\n",scanrangelo,scanrangehi);
		    newpulse = false;
		    for (int i=1; i<npulses; i++) {
			suslo = uapulses[i].start;
			sushi = uapulses[i].start+uapulses[i].width;
			if ((scanrangelo>=suslo&&scanrangelo<=sushi)||(scanrangehi>=suslo&&scanrangehi<=sushi)||(suslo>=scanrangelo&&suslo<=scanrangehi)){
//			    printf("Found associated pulse: DM %f, SNR %f, loc %i, start %i\n",uapulses[i].dm,uapulses[i].SNR,uapulses[i].loc,uapulses[i].start);
			    if (suslo < scanrangelo) {
				scanrangelo = suslo;
			    	newpulse = true;
			    }
			    if (sushi>scanrangehi) {
			    	scanrangehi = sushi;
			    	newpulse = true;
			    }
			    if (uapulses[i].SNR>maxsnr) {
//				printf("newp(%f) is more than maxsnr(%f)\n",uapulses[i].SNR,maxsnr);
			    	maxsnr=uapulses[i].SNR;
			    	maxsnrindex=ninonecand+1;
//				printf("new index for SNR %f is %d\n",maxsnr,maxsnrindex);
			    }
			    onecand.push_back(uapulses[i]);
//			    fprintf(stderr,"loc: %d start: %d wid: %d\n",uapulses[i].loc,uapulses[i].start,uapulses[i].width);
			    uapulses.erase(uapulses.begin()+i); //remove associated element
//			    fprintf(stderr,"i,npulses,uap.size,ninonecand: %d,%d,%d,%d\n",i,npulses,uapulses.size(),ninonecand);
			    onecand[0].width = scanrangehi-scanrangelo;
			    ninonecand++;
			    npulses--;
			    i--;
//			    fprintf(stderr,"cand loc: %d start: %d wid: %d\n\n",onecand[ninonecand].loc,onecand[ninonecand].start,onecand[ninonecand].width);
			}
		    }
		}
//		printf("finished loop. Max SNR at index %d was %f\n",maxsnrindex,onecand[maxsnrindex].SNR);
		if (onecand.size()<=3 || onecand[maxsnrindex].dm<irrel) {
// 		    printf("The burst is RFI. %d %d %f %f\n\n",onecand.size(),maxsnrindex,onecand[maxsnrindex].dm,onecand[maxsnrindex].SNR);
		    RFIpulses.insert(RFIpulses.end(), onecand.begin()+1, onecand.end());
		} else {
		    onecand[0].put_pulse(onecand[maxsnrindex].amp,
					 onecand[maxsnrindex].SNR, scanrangelo,
					 onecand[maxsnrindex].loc, scanrangehi-scanrangelo,
					 onecand[maxsnrindex].tscrfac, onecand[maxsnrindex].dm);
		    candidates[ncandidates] = onecand;
		    ncandidates++;
		}
		onecand.clear();
	}
	*nmultibeamcands = ncandidates;
	return (candidates);
}

*/
