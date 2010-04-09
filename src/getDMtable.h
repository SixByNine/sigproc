// Creates a list of accelerations for each DM 
// Lina Levin
// 2008-06-25

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

void   getDMtable(float DMstart, float DMmax, double tsamp, double ti, double bw, double cfreq, int Nchan, double tol,int* Ndms,float* &DMtable);
void   getDMtable(float DMmax, double tsamp, double ti, double bw, double cfreq, int Nchan, double tol,int* Ndms,float* &DMtable);
void getacctable(float DM, double amax, double tsamp, double ti, double tint, double bw, double cfreq, int Nchan, double tol, int* Nacc, float* &acctable);

 


