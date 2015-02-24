#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cinttypes>
#include <vector>
#include <algorithm>
#include <valarray>
#include <cstring>

// sigproc headers
#include "sigproc.h"
#include "header.h"
#include "mjk_cmd.h"


#define DM_CONST_T 2.41e-4
#define DM_CONST   1.0e12/DM_CONST_T

extern "C" {
// horible sigproc globals at work.
FILE *input, *output, *logfile;
char  inpfile[180], outfile[180], ignfile[180];
int ascii, asciipol, stream, swapout, headerless, nbands, userbins, usrdm, baseline, clipping, sumifs, profnum1, profnum2, nobits, wapp_inv, wapp_off;
double refrf,userdm,fcorrect;
float clipvalue,jyfactor,jyf1,jyf2;
}


using namespace std;

template <typename T>
using vec3d = vector<vector<valarray<T>>>;

template <typename T>
uint32_t dedisperse(const size_t Nsamps, const size_t Nchans,
	  const T *raw, vector<float> freq,
	  const float tsamp, const float dmmax,
	  vector<float> &dms, vector<valarray<T>> &output);

int main(int argc, char** argv){

   FILE *infile;
   char fname[1024];

   strcpy(fname,getS("--file","-f",argc,argv,"in.fil"));
   const float dmmax=getF("--dmmax","-d",argc,argv,100);

   infile=fopen(fname,"r");
   const size_t hdrlen = read_header(infile);
   const size_t nsamp = nsamples(fname,hdrlen,nbits,nifs,nchans);

   vector<float> freq(nchans);
   for (size_t c = 0; c < nchans; c++){
	  freq[c] = 1e6*(c*foff + fch1);
   }

   vector<float> dms(0);
   vector<valarray<float>> dedispersed(0);
   float* raw;

   const uint32_t nsblk=10000;
   raw=static_cast<float*>(calloc(nsblk*nchans,sizeof(float)));

   uint64_t nread=0;
   double obsnT=0;

   while(!feof(infile)){
	  fseek(infile,hdrlen+nread,SEEK_SET);
	  printf("T=%lf\n",obsnT);
	  size_t read = fread(raw,sizeof(float),nchans*nsblk,infile);
	  printf("read %d floats\n",read);
	  uint32_t nout = dedisperse(nsblk,nchans,raw,freq,tsamp,dmmax,dms,dedispersed);
	  nbands=1;
	  nobits=32;
	  char ofname[1024];
	  for (uint32_t idm=0; idm < dms.size(); idm++){
		 sprintf(ofname,"%s-%05.2f.tim",fname,dms[idm]);
		 if(nread==0){
			userdm=dms[idm];
			output = fopen(ofname,"w");
			dedisperse_header();
		 } else {
			output = fopen(ofname,"a");
		 }
		 fwrite(&(dedispersed.at(idm)[0]),nout,sizeof(float),output);
		 fclose(output);
	  }
	  obsnT+=tsamp*read/nchans;
	  nread += read*sizeof(float);
   }



   fclose(infile);
   return 0;

}

float isq(float f1, float f2){
   return (DM_CONST)*(pow(f1,-2)-pow(f2,-2));
}

int32_t get_t1(int32_t dT, int32_t t0, float f0, float f2){
   float f1=(f2+f0)/2.0;
   float C = isq(f1,f0) / isq(f2,f0);
   return t0-dT*C;
}

float dm_to_dt(float dm, float f1, float f2, float tsamp){
   return (dm*isq(f1,f2)/tsamp);
}

float dt_to_dm(uint32_t dt, float f1, float f2, float tsamp){
   return dt*tsamp / isq(f1,f2);
}


template<typename T>
vec3d<T> *allocate3d(size_t dTmax,size_t  nsub,size_t nsamp, T def){
   printf("Allocate3d %dx%dx%d = %.1e Bytes\n",dTmax,nsub,nsamp,(double)(dTmax*nsub*nsamp*sizeof(T)));
   valarray<T> v(def,nsamp);
   vector<valarray<T>> vv(nsub,v);
   vec3d<T> * ret = new vec3d<T>(dTmax,vv);
   return ret;
}


template <typename T>
uint32_t dedisperse(const size_t Nsamps, const size_t Nchans,
	  const T *raw, vector<float> freq,
	  const float tsamp, const float dmmax,
	  vector<float> &dms, vector<valarray<T>> &output){

   valarray<T> v((T)0,Nsamps);
   vector<valarray<T> >* data = new vector<valarray<T> > ( Nchans, v);


   float df = freq[1] - freq[0];
   bool reverse_freqs=false;

   if (df < 0){
	  // reverse frequencies
	  df = -df;
	  reverse(freq.begin(),freq.end());
	  reverse_freqs=true;
   }

   float f0 = freq[0]-0.5*df;
   float f2 = freq[1]+0.5*df;



   //transpose
   printf("transpose");
   for (uint32_t ichan =0; ichan < Nchans; ichan++){
	  uint32_t iichan=ichan;
	  if(reverse_freqs)iichan=Nchans-ichan-1;
	  for (uint32_t isamp =0; isamp < Nsamps; isamp++){
		 data->at(ichan)[isamp] = raw[iichan + isamp*Nchans];
	  }
   }




   uint32_t dTmax = (uint32_t)round(dm_to_dt(dmmax,f0,f2,tsamp))+1;

   vec3d<T> *A_out;
   vec3d<T> *A_in = allocate3d(dTmax,Nchans,Nsamps,(T)0);
   printf("Round 1\n");
   for (uint32_t dT =0; dT < dTmax ; dT++){
	  for (uint32_t c =0; c < Nchans ; c++){
		 A_in->at(dT)[c] *=0;
		 for (uint32_t i =0; i < dT+1 ; i++){
			A_in->at(dT)[c] += data->at(c).shift(-i);
		 }
	  }
   }

   uint32_t Nsub=2;

   while (Nsub <= Nchans){
	  if(Nsub < Nchans){
		 f0=freq[0];
		 f2=freq[Nsub*2-1];
		 dTmax=(int)(ceil(dm_to_dt(dmmax,f0,f2,tsamp)))+1;
	  } else {
		 f0=freq[0];
		 f2=freq[Nchans-1];
		 dTmax = (int)(round(dm_to_dt(dmmax,f0,f2,tsamp)));
	  }
	  uint32_t Npairs = Nchans/Nsub;

	  A_out = allocate3d(dTmax,Npairs,Nsamps,(T)0);
	  printf("Round %d\n",Nsub);

	  for (uint32_t ipair = 0; ipair < Npairs; ipair++){
		 uint32_t i=ipair*2;
		 uint32_t j=ipair*2+1;
		 f0=freq[ipair*Nsub];
		 f2=freq[(ipair+1)*Nsub-1];
		 for (uint32_t dT = 0; dT < dTmax; dT++){
			int32_t t0 = 0;
			int32_t t1 = get_t1(dT,t0,f0,f2);
			int32_t t2 = -dT;
			int32_t d0 = t0-t1;
			int32_t d1 = t1-t2;
			if (d0 >= A_in->size()){
			   d0=A_in->size()-1;
			}
			if (d1 >= A_in->size()){
			   d1=A_in->size()-1;
			}
			A_out->at(dT)[ipair] += A_in->at(d0)[i] + A_in->at(d1)[j].shift(-t1);
		 }
	  }
	  A_in->swap(*A_out);
	  A_out->clear();
	  delete A_out;
	  Nsub*=2;
   }

   dms.clear();
   f0=freq[0];
   f2=freq[Nchans-1];
   uint32_t nout=Nsamps - dTmax;
   printf("%d %d %d\n",nout,Nsamps,dTmax);
   output.clear();
   for (uint32_t dT = 0; dT < dTmax; dT++){
	  dms.push_back(dt_to_dm(dT,f0,f2,tsamp));
	  output.push_back(A_in->at(dT)[0]);
   }
   return nout;
}


