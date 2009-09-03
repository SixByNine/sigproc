using namespace std;


void filavg(int *ntime, int nchans, float *data);
//void filchanavg(int ntime, int *navchans, float *data);
void setcolortable(int colortable);
void printhelp();
void plotdm(int Sdec, long long int Sskip, float dmoff, float inpDM, int filnchans, double filtsamp, double filfoff, double filfch1);
void plotlambda(int Sdec, long long int Sskip, float dmoff, float inpDM,int filnchans, double filtsamp, double filfoff, double filfch1);
void killchan(float *data, int ntime, int chan, int tscr, double filnbits, int filnchans);
float filgetmax(float *data, int arraysize);
float filgetmin(float *data, int arraysize);
float* filint2float(int *array, int arraysize);
//void plotfil(char *currentfile, long long int Sskip, int Sread, int Sdec, float inpDM);
