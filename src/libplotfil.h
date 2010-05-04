using namespace std;


void filavg(int *ntime, int numchans, float *data);
void filchanavg(int ntime, int *numchans, float *data);
void setcolortable(int colortable);
void printhelp();
void plotdm(long long int Sskip, float dmoff, float inpDM, int filnchans, double filtsamp, double filfoff, double filfch1);
void plotlambda(long long int Sskip, float dmoff, float inpDM,int filnchans, double filtsamp, double filfoff, double filfch1);
void killchan(float *data, int ntime, int chan, float fillvalue, int filnchans);
float filgetmax(float *data, int arraysize);
float filgetmin(float *data, int arraysize);
float* filint2float(int *array, int arraysize);
void chanbaseline(float *data, int ntime, int filnchans);
void plotfilhelp();
void plotfil(float *floatarchive, int nFBsamps, long long int Sskip, float inpDM, int pgpID, double filtsamp, int filnchans, double filfch1, double filfoff, int colortable,bool askdevice);
void plotfil(float *floatarchive, int nFBsamps, long long int Sskip, float inpDM, int pgpID, double filtsamp, int filnchans, double filfch1, double filfoff, int colortable);
bool getkillfile(int * killmask,int nchans,char *killfile);
