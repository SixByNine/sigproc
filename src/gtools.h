#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;
// Set of tools to do various useful gsearch things
//    Sarah Burke
//    Swinburne University
//    May, 2008
//
//  This .h file contains:
//      -definition of the GPULSE class
//      -the giantsearch searching code
//      -some useful functions


//**************************************************
//**************************************************
//***            GIANT PULSE CLASS               ***
//**************************************************
//**************************************************
class Gpulse{
public: //temporarily all are public.
    float amp, SNR, dm;       // amplitude, significance
    int start,loc, width, tscrfac, beam; // integer bin start, bin location, and bin width of pulse
    //float zSNR, DM;     // SNR at 0 dispersion measure, and DM of detected pulse
    Gpulse ();
    Gpulse (float,float,int,int,int,int,float);
    Gpulse (float,float,int,int,int,int,float,int);
    void put_pulse(float,float,int,int,int,int,float,int);
    void put_pulse(float,float,int,int,int,int,float);
    //int rank () {return (SNR/zSNR);}
};




//**************************************************
//**************************************************
//***          function declarations             ***
//**************************************************
//**************************************************


class GPulseState{
 public:
    vector<Gpulse> *DMtrials;
    int NDMtrials;
    GPulseState(int ndms);
    void searchforgiants(int itrial, int numbersamples, int offset, float * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac);
    void searchforgiants(int itrial, int numbersamples, int offset, unsigned short int * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac);
    int* givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, char* filetimestamp);
    int* givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, char* filetimestamp, int beamID);
    int* givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, char* filetimestamp, int beamID, char* resultsfilename);
    void selfdestruct();
};
void normalise(float *data, int n);                  //normalize data
double normalise(int n,float * d, double * dataaverage);//normalize and report SIGMA, MEAN
double getrms(int n,float * d, double * dataaverage);//report SIGMA
double getrms(int n, unsigned short int * d, double * dataaverage);//report SIGMA
double submeanrms(int n,float * d, double * dataaverage);//subtract mean and report SIGMA
void timeavg(int n, float * d);                      //scrunch data by factor of 2
double timeavg(int n, float * d, double *dataaverage);//scrunch data by factor of 2,ret.RMS and subtract mean
float getmax(float *data, int arraysize);            //get max of data
float getmax(float *data, int arraysize, int *index);//get max of data and provide array index of maximum
float getmin(float *data, int arraysize);            //get min of data
float* int2float(int *array, int arraysize);         //convert int array to float array
void removebaseline(unsigned short int *indata, float* outdata, int ndat, int runmeansize, float thresh);
double getmowedsigma(int n, float * d, double unmowedsigma, double mean);
vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands,char* resultsfilename,char* filetimestamp,int beamID);
vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands,char* resultsfilename,char* filetimestamp,int beamID,float irrel);
vector<Gpulse>* beamassoc_giants(vector<Gpulse> *beams, int nbeams, int *nmultibeamcands);
vector<Gpulse>* beamassoc_giants(vector<Gpulse> *beams, int nbeams, int *nmultibeamcands, float irrel);


vector<Gpulse> giantsearch(int n, float * data, float thresh, double RMS, float bintol, int tscrfac, float dm);
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol, int usertscrfac, float dm, int starttscrfac);
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol, int usertscrfac, float dm);
vector<Gpulse> findgiants(int npts, unsigned short int * data, float nsigma, float bintol, int usertscrfac, float dm, int starttscrfac);
vector<Gpulse> findgiants(int npts, double * data, float nsigma, float bintol, int usertscrfac, float dm);
vector<Gpulse> findgiants(int npts, int * data, float nsigma, float bintol, int usertscrfac, float dm);
vector<Gpulse> findgiants(int npts, unsigned char * data, float nsigma, float bintol, int usertscrfac, float dm);





