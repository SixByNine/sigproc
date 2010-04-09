// Creates a list of accelerations for each DM 
// Lina Levin
// 2008-06-25

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

void   getDMtable(float DMstart, float DMmax, double tsamp, double ti, double bw, double cfreq, int Nchan, double tol,int* Ndms,float* &DMtable){
        
        *Ndms = 0;
        float DM = DMstart;
	int j;
	
/*	cout << "DMstart is " << DMstart << " DMmax is " << DMmax << " Ndms is " << *Ndms << endl;

	cout << "tsamp is " << tsamp << endl;
	cout << "ti is " << ti << endl;
	cout << "bw is " << bw << endl;
	cout << "cfreq is " << cfreq << endl;
	cout << "Nchan is " << Nchan << endl;
	cout << "tol is " << tol << endl;
*/

	while(DM<=DMmax){
	    double oldDM = DM;
	    double t00 = sqrt(tsamp*tsamp + ti*ti + pow(8.3*bw*DM/pow(cfreq,3.0),2));
	    double a = 8.3*bw/pow(cfreq,3);
	    double b = 8.3*Nchan*bw/4/pow(cfreq,3);
	    double c = tol*tol*t00*t00 - tsamp*tsamp - ti*ti;
	    double newDM = (b*b*DM + sqrt(-a*a*b*b*DM*DM + a*a*c + b*b*c))/(a*a+b*b);
	    //DM = oldDM + 2.0*(newDM-oldDM);
	    DM = newDM;
//            cout << "DM is " << DM << endl;
	    (*Ndms)++;
	    }

	//cout << "Ndms before new is " << *Ndms << endl;

	DMtable = new float[*Ndms];
	*Ndms = 0;
	DM = DMstart;
	while(DM<=DMmax){
	    DMtable[*Ndms]=DM;
	    double oldDM = DM;
	    double t00 = sqrt(tsamp*tsamp + ti*ti + pow(8.3*bw*DM/pow(cfreq,3.0),2));
	    double a = 8.3*bw/pow(cfreq,3);
	    double b = 8.3*Nchan*bw/4/pow(cfreq,3);
	    double c = tol*tol*t00*t00 - tsamp*tsamp - ti*ti;
	    double newDM = (b*b*DM + sqrt(-a*a*b*b*DM*DM + a*a*c + b*b*c))/(a*a+b*b);
	    //DM = oldDM + 2.0*(newDM-oldDM);
	    DM = newDM;
	    (*Ndms)++;
	    }
	
	//cout << "Ndms after new float is " << *Ndms << endl;
}

void getacctable(float DM, double amax, double tsamp, double ti, double tint, double bw, double cfreq, int Nchan, double tol, int* Nacc, float* &acctable){
    
        *Nacc = 0;
        double d_a; 

  	//cout <<" DM = "<< DM << " acc = ";
	
	double acc = 0.0;
	while (acc < amax){
	    double e = (tol*tol-1)*(24.0*3e8)*(24.0*3e8)*1e-12*(tsamp*tsamp+ti*ti)/pow(tint,4);
	    double f = (tol*tol-1)*(24.0*3e8)*(24.0*3e8)*1e-12*8.3*8.3*bw*bw/
		pow(cfreq,6)/pow(tint,4);
	    d_a = sqrt(e + f*DM*DM);
	    //cout << acc << " ";
	    acc+=2*d_a;
	    (*Nacc)++;
	    }
	cout << endl;
	//cout << "Nacc is " << *Nacc << endl;

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
 




void   getDMtable(float DMmax, double tsamp, double ti, double bw, double cfreq, int Nchan, double tol,int* Ndms,float* &DMtable){
    getDMtable(0.0,DMmax,tsamp,ti,bw,cfreq,Nchan,tol,Ndms,DMtable);
}
