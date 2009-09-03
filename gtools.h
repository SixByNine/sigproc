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
    int finishedcount;
    GPulseState(int ndms);
    void searchforgiants(int itrial, int numbersamples, float * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac);
    void searchforgiants(int itrial, int numbersamples, unsigned short int * data, float nsigma, float bintol, int usertscrfac, float dm,int starttscrfac);
    int* givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel);
    int* givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, int beamID);
    int* givetimes(int* ndetected, float sampletime, float flo, float fhi, float irrel, int beamID, char* resultsfilename);
};
void normalise(float *data, int arraysize);          //normalize data
double normalise(int n,float * d, double * dataaverage);//normalize and return RMS; provide dataaverage
double getrms(int n,float * d, double * dataaverage);//get RMS of data
void timeavg(int n, float * d);                      //scrunch data by factor of 2
double timeavg(int n, float * d, double *dataaverage);//scrunch data by factor of 2,ret.RMS and subtract mean
float getmax(float *data, int arraysize);            //get max of data
float getmax(float *data, int arraysize, int *index);//get max of data and provide array index of maximum
float getmin(float *data, int arraysize);            //get min of data
float* int2float(int *array, int arraysize);         //convert int array to float array
vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands);
vector<Gpulse>* assoc_giants(vector<Gpulse> uapulses, int *nsinglebeamcands, float irrel);
vector<Gpulse>* beamassoc_giants(vector<Gpulse> *beams, int nbeams, int *nmultibeamcands);
vector<Gpulse>* beamassoc_giants(vector<Gpulse> *beams, int nbeams, int *nmultibeamcands, float irrel);


vector<Gpulse> giantsearch(int n, float * data, float thresh, double RMS, float bintol, int tscrfac, float dm);
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol, int usertscrfac, float dm, int starttscrfac);
vector<Gpulse> findgiants(int npts, float * data, float nsigma, float bintol, int usertscrfac, float dm);
vector<Gpulse> findgiants(int npts, unsigned short int * data, float nsigma, float bintol, int usertscrfac, float dm, int starttscrfac);
vector<Gpulse> findgiants(int npts, double * data, float nsigma, float bintol, int usertscrfac, float dm);
vector<Gpulse> findgiants(int npts, int * data, float nsigma, float bintol, int usertscrfac, float dm);
vector<Gpulse> findgiants(int npts, unsigned char * data, float nsigma, float bintol, int usertscrfac, float dm);





