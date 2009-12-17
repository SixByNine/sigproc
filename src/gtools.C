#include <stdlib.h>
#include "gtools.h"

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

//state stuff
GPulseState::GPulseState(int ndms) {
    DMtrials = new vector<Gpulse>[ndms]; //!!!!!!!I THINK THIS IS WHERE THE MEMORY LEAK IS.
    NDMtrials = ndms;
}

//NOTE: searchforgiants works fine; it preserves previous runs of dmtrials.
void GPulseState::searchforgiants(int itrial, int numbersamples, float * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac) {
    DMtrials[itrial] = findgiants(numbersamples, data, nsigma, bintol, usertscrfac, dm,starttscrfac);
//    printf("in searchforgiants the dm before %f, found second giant at SNR %f\n",dm, DMtrials[itrial].begin()->SNR);
}

void GPulseState::searchforgiants(int itrial, int numbersamples, unsigned short int * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac) {
    DMtrials[itrial] = findgiants(numbersamples, data, nsigma, bintol, usertscrfac, dm,starttscrfac);
//    printf("in searchforgiants the dm before %f, found second giant at SNR %f\n",dm, DMtrials[itrial].begin()->SNR);
//    fprintf(stderr,"At DM %f found %d candidates\n",dm,DMtrials[itrial].size());
}


int* GPulseState::givetimes(int* ndetected, float sampletime, float flo,float fhi,float irrel) {
    return(givetimes(ndetected,sampletime,flo,fhi,irrel,-1));
}

int* GPulseState::givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, int beamID){
    return(givetimes(ndetected, sampletime, flo, fhi, irrel,beamID,"GResults.txt"));
}

int* GPulseState::givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, int beamID, char* resultsfilename) {
    //look for giants and return array of data to save. or specifications
    //about data to save
    // The BEAM ID should ONLY be specified if wanting to do A MULTIBEAM SEARCH.
    //---------------------------------------------
    int nsinglebeamcands, delayinsamples, totdelay;
    float delayinms;
    vector<Gpulse> suspectvectorstorage;
    Gpulse gpulsestorage;
    FILE* resultsfile = fopen(resultsfilename,"a");
    
    for (int i=0; i<NDMtrials-1; i++) {
	suspectvectorstorage.insert(suspectvectorstorage.end(),
				    DMtrials[i].begin(), DMtrials[i].end());
    }
    vector<Gpulse>* suspectarraystorage = assoc_giants(suspectvectorstorage,&nsinglebeamcands,irrel);
//	fprintf(stderr,"sampletime: %f flo:%f fhi:%f\n",sampletime,flo,fhi);
    fprintf(stderr,"\n\nN candidates in this block before associating: %d\n",suspectvectorstorage.size());
    fprintf(stderr,"N candidates in this block after associating: %d\n\n",nsinglebeamcands);
    fprintf(resultsfile,"\n\nN candidates in this block before associating: %d\n",suspectvectorstorage.size());
    fprintf(resultsfile,"N candidates in this block after associating: %d\n\n",nsinglebeamcands);
    int* timestamps = new int[nsinglebeamcands*2];
    if (beamID<0){
	for (int i=0; i<(nsinglebeamcands*2); i+=2) {
	    gpulsestorage = suspectarraystorage[i/2].at(0);
	    if (flo<fhi)
		delayinms = gpulsestorage.dm * 4.15 * (pow(flo/1000, -2) - pow(fhi/1000, -2));
	    if (fhi<flo)
		delayinms = gpulsestorage.dm * 4.15 * (pow(fhi/1000, -2) - pow(flo/1000, -2));
	    delayinsamples = (int)(delayinms/(sampletime*1000))+1;
	    fprintf(stderr,"Candidate %4d: DM %5.2f SNR %5.2f SCR %d\n",i/2,gpulsestorage.dm,gpulsestorage.SNR,gpulsestorage.tscrfac);
	    fprintf(resultsfile,"Candidate %4d: DM %5.2f SNR %5.2f SCR %4d STARTBIN %13d PEAK %13d\n",i/2,gpulsestorage.dm,gpulsestorage.SNR,gpulsestorage.tscrfac,gpulsestorage.start,gpulsestorage.loc);
	    totdelay = delayinsamples+gpulsestorage.width;
	    if (gpulsestorage.start-(totdelay)<0)
		timestamps[i] = 0;
	    else
		timestamps[i] = gpulsestorage.start-(totdelay);
	    timestamps[i+1] = 3*(totdelay);
	}
	*ndetected = nsinglebeamcands;
	suspectvectorstorage.clear();
	fclose(resultsfile);
	
	return (timestamps);
    } else {
	//do something
    }
    fclose(resultsfile);
    return(0);
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

//*******************************************
//       TIMESERIES SEARCH ALGORITHM
//*******************************************
vector<Gpulse> giantsearch(int n, float * data, float thresh, double RMS, float bintol, int tscrfac, float dm) {
	int i;
	Gpulse pulse;
	bool detected = false, first = true;
	int xstart = 0, width = 0; //first bin of detected pulse
	int j, imax, ngiants=0, lastend=0; //counters
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
//***             GSEARCHING LOOPS               ***
//**************************************************
//**************************************************

// for 32-bit (float) data
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol,
		int usertscrfac, float dm, int starttscrfac) {
	vector<Gpulse> mybigvector, giant;
	double mean=0;
	double sigma=0;
	for (int tscrfac=starttscrfac; tscrfac<=usertscrfac; tscrfac*=2) {
		int newnpts = npts*starttscrfac/tscrfac;

		if (tscrfac==starttscrfac)
			sigma = getrms(newnpts, data, &mean);
		else
			sigma = timeavg(2*newnpts, data, &mean);
		//	printf("sigma is %g\n",sigma);
		giant = giantsearch(newnpts, data, nsigma*sigma, sigma,	bintol/tscrfac, tscrfac, dm); //	    if (giant.size()<1) printf("NO GIANTS FOUND IN %s\n",filename[0]); //	    printf("Giant pulse candidates in file %s:\n------filen------\t---amp---\t---SNR---\t--peak bin--\t---width---\n",filename[i]);
		mybigvector.insert(mybigvector.end(), giant.begin(), giant.end());
	}
	return (mybigvector);
}

vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol,
		int usertscrfac, float dm) {
	return findgiants(npts, data, nsigma, bintol, usertscrfac, dm, 1);
}

// for unsigned short int data, as in the HITRUN survey dedisperser
vector<Gpulse> findgiants(int npts, unsigned short int * data, float nsigma, float bintol,
			  int usertscrfac, float dm,int starttscrfac) {
	vector<Gpulse> mybigvector;
	float *convertedarray;
	convertedarray = (float *) malloc(npts*sizeof(float));
	for (int j=0; j<npts; j++) {
		convertedarray[j]=(float)data[j];
	}
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,
			dm);
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
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,
			dm);
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
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,
			dm);
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
	mybigvector = findgiants(npts, convertedarray, nsigma, bintol, usertscrfac,
			dm);
	free(convertedarray);
	return (mybigvector);
}

//**************************************************
//**************************************************
//***             USEFUL FUNCTIONS               ***
//**************************************************
//**************************************************


//**************************************************
//  Normalizes an array d of n elements by the RMS
//**************************************************
void normalise(float *data, int arraysize) {
	double sum=0.0, sumsquares=0.0;
	int i=0;
	while (i<arraysize) {
		sum += data[i];
		sumsquares += (data[i] * data[i]);
		i++;
	}
	double mean=sum/(double)arraysize;
	double meansquares=sumsquares/(double)arraysize;
	double sigma= sqrt(meansquares - (mean * mean));
	for (i=0; i<arraysize; i++)
		data[i]=(data[i]-mean)/sigma;//    cout<<"normalization:\n"<<"mean == "<<mean<<"\nsigma = "<<sigma<<"\n";
}

//**************************************************
//  Normalizes an array d of n elements by the RMS
//      and gives SIGMA and DATA AVERAGE back
//**************************************************
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
	for (int i=0; i<n; i++)
		d[i]=(d[i]-mean)/sigma; //    fprintf(stderr,"IN NORMALIZE FUNCTION:\n\tRMS of data is: %g\n\tMean is: %g\n",sigma,mean);
	return (sigma);
}

//**************************************************
// subtract mean from data and return 1sigma value
//**************************************************
double getrms(int n, float * d, double * dataaverage) {
	double sum=0.0;
	double sumsq=0.0;
	for (int i=0; i<n; i++) {
		sum+=d[i];
		sumsq+=d[i]*d[i];
	}
	double mean=sum/(double)n;
	*dataaverage = mean;
	double meansq=sumsq/(double)n;
	double sigma=sqrt(meansq-mean*mean);
	if (sigma==0) {
		fprintf(stderr, "getrms::RMS of data is zero\n");
		return (sigma);
	}
	int i;
#pragma omp parallel for private(i)
	for (i=0; i<n; i++)
		d[i]=(d[i]-mean);
	//    fprintf(stderr,"IN GETRMS FUNCTION:\n\tRMS of data is: %g\n\tMean is: %g\n",sigma,mean);
	return (sigma);
}

//**************************************************
// Time-avg data by a factor of two & return RMS
//**************************************************
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
	//    printf("asize is %d\n",asize);
	double mean=sum/(double)asize;
	*dataaverage = mean;
	double meansq=sumsq/(double)asize;
	double sigma=sqrt(meansq-mean*mean);
	if (sigma==0) {
		fprintf(stderr, "timeavg::RMS of data is zero\n");
		return (sigma);
	}
	//    fprintf(stderr,"IN TIMEAVG FUNCTION:\n\tRMS of data is: %g\n\tMean is: %g\n",sigma,mean);
	return (sigma);
}

//**************************************************
//    Time-average the data by a factor of two.
//**************************************************
void timeavg(int n, float * d) {
	int asize;
	float avg;
	for (int i=0; i<n-1; i+=2) {
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
	char *temparray;
	for (int i=0; i<arraysize; i++) {
		floatarray[i] = float(array[i]);
	}
	return (floatarray);
}

//**************************************************
//      first-pass RFI excision; assoc_giants
//**************************************************

vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands) {
	vector<Gpulse> *cands;
	cands = assoc_giants(uapulses, nsinglebeamcands, 5.0); //DEFAULT IRREL DM = 5.0
	return (cands); //IS nsinglebeamcands IMPLEMENTED CORRECTLY HERE?
}

vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands,float irrel) {
	int maxdm=0, scanrangelo, scanrangehi, ncandidates=0;
	int npulses = uapulses.size(),suslo,sushi,ninonecand,maxsnrindex=1;
	vector<Gpulse> *candidates;
	vector<Gpulse> onecand, RFIpulses;
	Gpulse placeholder;
	bool newpulse;
	float maxsnr;

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
//		    printf("!!GOOD!! %4d %10d %5.2f %5.3f\n",onecand.size(),onecand[0].start,onecand[maxsnrindex].dm,onecand[maxsnrindex].SNR);
		    ncandidates++;
		}
		onecand.clear();
	}
	*nsinglebeamcands = ncandidates;
	return (candidates);
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
//		    printf("!!GOOD!! %4d %10d %5.2f %5.3f\n",onecand.size(),onecand[0].start,onecand[maxsnrindex].dm,onecand[maxsnrindex].SNR);
		    ncandidates++;
		}
		onecand.clear();
	}
	*nmultibeamcands = ncandidates;
	return (candidates);
}

*/
