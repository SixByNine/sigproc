// Command line program to call getDMtable.C
// and get an acceleration table.
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
//#include "getDMtable.h"

using namespace std;

void getacctable(float DM, double amax, double tsamp, double ti, double tint, double bw, double cfreq, int Nchan, double tol, int* Nacc, float* &acctable);

int main(int argc, char *argv[]){
    float DM=0;
    double accmax,tsamp,tintrinsic,ttotal,bw,cfreq,tolerance;
    int Nchan,Nacc;
    float* acctable;
    char* inthing;

    int i=1;
    while (i<argc) {
	sscanf(argv[i],"%2c",inthing);
	cout<<inthing<<"\n";
	if (inthing == "-d"){ sscanf(argv[++i],"%f",DM);cout<<"got DM";}         // pc/cc
	else if (inthing == "-a") sscanf(argv[++i],"%L",&accmax);     // s/s?
	else if (inthing == "-t") sscanf(argv[++i],"%L",&tsamp);      // us
	else if (inthing == "-w") sscanf(argv[++i],"%L",&tintrinsic); // us
	else if (inthing == "-T") sscanf(argv[++i],"%L",&ttotal);     // s?
	else if (inthing == "-c") sscanf(argv[++i],"%L",&bw);         // 
	else if (inthing == "-C") sscanf(argv[++i],"%L",&cfreq);      // 
	else if (inthing == "-Nc") sscanf(argv[++i],"%d",&Nchan);     
	else if (inthing == "-f") sscanf(argv[++i],"%L",&tolerance);  // 1.25
	else if (inthing == "-Na") sscanf(argv[++i],"%d",&Nacc);
//	if (argv[i] == "-") sscanf(argv[++i],"%d",&);
//	if (strings_equal(argv[i],"-n"))       sscanf(argv[++i],"%d",&ngulp_original);
	i++;
    }
    
    cout<<DM<<" "<<accmax<<" "<<tsamp<<" "<<tintrinsic<<" \n";
    getacctable (DM,accmax,tsamp,tintrinsic,ttotal,bw,cfreq,Nchan,tolerance,&Nacc,acctable);

    for (i=0;i<Nacc;i++){
	printf("%f\n",acctable[i]);
    }

    return(0);
}

void getacctable(float DM, double amax, double tsamp, double ti, double tint, double bw, double cfreq, int Nchan, double tol, int* Nacc, float* &acctable){
    
        *Nacc = 0;
        double d_a; 

  	cout <<" DM = "<< DM << " acc = ";
  	cout <<"tsamp "<<tsamp<<"\ntol "<<tol<<"\nti "<<ti<<"\ntint "<<tint<<"\n";
	
	double acc = 0.0;
	while (acc < amax){
	    double e = (tol*tol-1)*(24.0*3e8)*(24.0*3e8)*1e-12*(tsamp*tsamp+ti*ti)/pow(tint,4);
	    double f = (tol*tol-1)*(24.0*3e8)*(24.0*3e8)*1e-12*8.3*8.3*bw*bw/pow(cfreq,6)/pow(tint,4);
	    cout<< "e "<<e<<"\nf "<<f<<"\n";
	    d_a = sqrt(e + f*DM*DM);
	    cout << acc << " ";
	    acc+=2*d_a;
	    (*Nacc)++;
	    }
	cout << endl;
	cout << "Nacc is " << *Nacc << endl;

	acctable = new float[*Nacc];
	*Nacc = 0;
	acc = 0.0;
	while(acc < amax){
	    acctable[*Nacc]=acc;
	    double e = (tol*tol-1)*(24.0*3e8)*(24.0*3e8)*1e-12*(tsamp*tsamp+ti*ti)/pow(tint,4);
	    double f = (tol*tol-1)*(24.0*3e8)*(24.0*3e8)*1e-12*8.3*8.3*bw*bw/
		pow(cfreq,6)/pow(tint,4);
	    d_a = sqrt(e + f*DM*DM);
	    acc+=2*d_a;
	    (*Nacc)++;
	    }
}
