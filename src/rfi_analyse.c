#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <complex.h>
#include <inttypes.h>
#include <stdio.h>
#include <TKlog.h>
#include <TKmatrix.h>
#include <math.h>
#include <fftw3.h>
#include <cpgplot.h>
#include "mjk_cmd.h"
#include "mjk_random.h"
#include "sigproc.h"
#include "header.h"


float camp(complex float x){
	return sqrt(crealf(x)*crealf(x) + cimagf(x)*cimagf(x));
}



int main(int argc, char** argv){
	float tr[6];

	const float ZAP=32;
	const uint64_t TSIZE=18;
	const uint64_t zapE=64;
	fftwf_init_threads();
	fftwf_plan_with_nthreads(omp_get_max_threads());

	logmsg("Open file '%s'",argv[1]);
	FILE* f = fopen(argv[1],"r");

	int hdr_bytes = read_header(f);
	const uint64_t nskip = hdr_bytes;
	const uint64_t nchan = nchans;
	logmsg("Nchan=%"PRIu64", tsamp=%f",nchan,tsamp);
	mjk_rand_t *random = mjk_rand_init(12345);

	rewind(f);
	FILE* of = fopen("clean.fil","w");
	uint8_t hdr[nskip];
	fread(hdr,1,nskip,f);
	fwrite(hdr,1,nskip,of);
	const uint64_t nsamp_per_block=round(pow(2,TSIZE));

	logmsg("Tblock = %f",nsamp_per_block*tsamp);

	mjk_clock_t *t_all = init_clock();
	start_clock(t_all);

	mjk_clock_t *t_read = init_clock();
	mjk_clock_t *t_trns= init_clock();
	mjk_clock_t *t_rms = init_clock();
	mjk_clock_t *t_fft = init_clock();
	mjk_clock_t *t_spec = init_clock();

	const uint64_t bytes_per_block = nchan*nsamp_per_block;

	uint8_t *buffer = calloc(bytes_per_block,1);
	float **data = malloc_2df(nchan,nsamp_per_block);
	float **clean = malloc_2df(nchan,nsamp_per_block);
	float *bpass = calloc(nchan,sizeof(float));
	float  *ch_var=NULL;
	float  *ch_mean=NULL;
	float  *ch_fft_n=NULL;
	float  *ch_fft_p=NULL;


	logmsg("Planning FFT - this will take a long time the first time it is run!");
	start_clock(t_fft);
	FILE * wisfile;
	if(wisfile=fopen("wisdom.txt","r")){
		fftwf_import_wisdom_from_file(wisfile);
		fclose(wisfile);
	}
	const int fftX=nsamp_per_block;
	const int fftY=nchan;
	const int fftXo=nsamp_per_block/2+1;

	float *X = fftwf_malloc(sizeof(float)*fftX);

	for (uint64_t i = 0; i < nsamp_per_block ; i++){
		X[i]=i;
	}
	float *tseries = fftwf_malloc(sizeof(float)*fftX);
	float complex *fseries = fftwf_malloc(sizeof(float complex)*fftXo);
	float *pseries = fftwf_malloc(sizeof(float)*fftXo);
	uint8_t *mask = malloc(sizeof(uint8_t)*fftXo);
	fftwf_plan fft_1d = fftwf_plan_dft_r2c_1d(fftX,tseries,fseries,FFTW_MEASURE|FFTW_DESTROY_INPUT);

	complex float * fftd = fftwf_malloc(sizeof(complex float)*(fftXo*fftY));
	fftwf_plan fft_plan = fftwf_plan_many_dft_r2c(
			1,&fftX,fftY,
			data[0] ,&fftX,1,fftX,
			fftd    ,&fftXo,1,fftXo,
			FFTW_MEASURE|FFTW_PRESERVE_INPUT);
	logmsg("Planning iFFT - this will take a long time the first time it is run!");
	fftwf_plan ifft_plan = fftwf_plan_many_dft_c2r(
			1,&fftX,fftY,
			fftd ,&fftXo,1,fftXo,
			clean[0] ,&fftX,1,fftX,
			FFTW_MEASURE|FFTW_PRESERVE_INPUT);

	if(!fft_plan){
		logmsg("Error - could not do FFT plan");
		exit(2);
	}

	wisfile=fopen("wisdom.txt","w");
	fftwf_export_wisdom_to_file(wisfile);
	fclose(wisfile);
	stop_clock(t_fft);
	logmsg("T(planFFT)= %.2lfs",read_clock(t_fft));
	reset_clock(t_fft);



	float min_var=1e9;
	float max_var=0;

	float min_fft_n=1e9;
	float max_fft_n=0;


	float min_fft_p=1e9;
	float max_fft_p=0;

	float min_mean=1e9;
	float max_mean=0;
	uint64_t nblocks=0;
	uint64_t totread=0;
	while(!feof(f)){
		nblocks++;
		ch_var = realloc(ch_var,nchan*nblocks*sizeof(float));
		ch_mean = realloc(ch_mean,nchan*nblocks*sizeof(float));
		ch_fft_n = realloc(ch_fft_n,nchan*nblocks*sizeof(float));
		ch_fft_p = realloc(ch_fft_p,nchan*nblocks*sizeof(float));
		start_clock(t_read);
		uint64_t read = fread(buffer,1,bytes_per_block,f);
		stop_clock(t_read);
		if (read!=bytes_per_block){
			nblocks--;
			break;
		}
		totread+=read;
		logmsg("read=%"PRIu64" bytes. T=%fs",read,totread*tsamp/(float)nchan);
		uint64_t offset = (nblocks-1)*nchan;
		start_clock(t_trns);
		// transpose with small blocks in order to increase cache efficiency.
#define BLK 8
#pragma omp parallel for schedule(static,2) shared(buffer,data)
		for (uint64_t j = 0; j < nchan ; j+=BLK){
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				for (uint64_t k = 0; k < BLK ; k++){
					data[j+k][i] = buffer[i*nchan+j+k];
				}
			}
		}

#pragma omp parallel for shared(data)
		for (uint64_t j = 0; j < nchan ; j++){
			if(j<zapE || (nchan-j) < zapE){ 
				for (uint64_t i = 0; i < nsamp_per_block ; i++){
					data[j][i]=ZAP;
				}
			}
		}

		if(nblocks==1){
#pragma omp parallel for shared(data,bpass)
			for (uint64_t j = 0; j < nchan ; j++){
				for (uint64_t i = 0; i < nsamp_per_block ; i++){
					bpass[j]+=data[j][i];
				}
				bpass[j]/=(float)nsamp_per_block;
				bpass[j]-=ZAP;
			}
		}
#pragma omp parallel for shared(data,bpass)
		for (uint64_t j = 0; j < nchan ; j++){
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				data[j][i]-=bpass[j];
			}
		}




		stop_clock(t_trns);

		start_clock(t_rms);
#pragma omp parallel for shared(data,ch_mean,ch_var)
		for (uint64_t j = 0; j < nchan ; j++){
			float mean=0;
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				mean+=data[j][i];
			}
			mean/=(float)nsamp_per_block;
			if(mean > ZAP+5 || mean < ZAP-5){
				logmsg("ZAP ch=%"PRIu64,j);
				for (uint64_t i = 0; i < nsamp_per_block ; i++){
					data[j][i]=ZAP;
				}
			}

			float ss=0;
			float x=0;
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				x = data[j][i]-mean;
				ss+=x*x;
			}
			float var=ss/(float)nsamp_per_block;
			if (var > 0){
				for (uint64_t i = 0; i < nsamp_per_block ; i++){
					float v = (data[j][i]-mean)/sqrt(var);
					if(v > 3 || v < -3){
						data[j][i]=mjk_rand_gauss(random)*sqrt(var)+mean;
					}
				}
			}

			ch_var[offset+j] = var;
			ch_mean[offset+j] = mean;

		}
		stop_clock(t_rms);

		for (uint64_t i = 0; i < nsamp_per_block ; i++){
			tseries[i]=0;
		}

		float tmean=0;
		float tvar=0;
		float max=0;
		float min=1e99;
		//#pragma omp parallel for shared(data,tseries)
		// NOT THREAD SAFE
		for (uint64_t j = 0; j < nchan ; j++){
			tmean+=ch_mean[offset+j];
			tvar+=ch_var[offset+j];
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				tseries[i]+=data[j][i];
				if(data[j][i]>max)max=data[j][i];
				if(data[j][i]<min)min=data[j][i];
			}
		}
		float ss=0;
		float mm=0;
		for (uint64_t i = 0; i < nsamp_per_block ; i++){
			float x=tseries[i]-tmean;
			mm+=tseries[i];
			ss+=x*x;
		}
		float rvar=ss/(float)nsamp_per_block;
		logmsg("var=%g tvar=%g",ss/(float)nsamp_per_block,tvar);
		logmsg("mean=%g tmean=%g",mm/(float)nsamp_per_block,tmean);
		cpgopen("3/xs");
		cpgsvp(0.1,0.9,0.1,0.9);
		cpgswin(0,fftX,tmean-sqrt(tvar)*30,tmean+sqrt(tvar)*30);
		cpgbox("ABN",0,0,"ABN",0,0);
		cpgline(fftX,X,tseries);
		cpgsci(2);
		cpgclos();
		tr[0] = 0.0 ;
		tr[1] = 1;
		tr[2] = 0;
		tr[3] = 0.5;
		tr[4] = 0;
		tr[5] = 1;

		logmsg("max=%g min=%g",max,min);

		cpgopen("4/xs");
		cpgsvp(0.1,0.9,0.1,0.9);
		cpgswin(0,nsamp_per_block,0,nchan);
		cpgbox("ABN",0,0,"ABN",0,0);
		cpggray(*data,nsamp_per_block,nchan,1,nsamp_per_block,1,nchan,tmean/(float)nchan+sqrt(rvar/(float)nchan),tmean/(float)nchan-sqrt(rvar/(float)nchan),tr);
		cpgclos();




		start_clock(t_fft);
		fftwf_execute(fft_1d);
		fftwf_execute(fft_plan);
		stop_clock(t_fft);

		{
			float T = sqrt(fftXo*tvar)*12;
			logmsg("Zap T=%.2e",T);

			float fx[fftXo];
			float fT[fftXo];
#pragma omp parallel for shared(fseries,pseries,mask)
			for (uint64_t i = 0; i < fftXo ; i++){
				mask[i]=1;
			}
#pragma omp parallel for shared(fseries,pseries,mask)
			for (uint64_t i = 0; i < fftXo ; i++){
				pseries[i]=camp(fseries[i]);
				fx[i]=i;
				float TT = T;
				if (i>512)TT=T/2.0;
				if(i>32){
					fT[i]=TT;
					if (pseries[i] > TT) {
						mask[i]=0;
					}
				} else fT[i]=0;
			}

			uint64_t nmask=0;
			for (uint64_t i = 0; i < fftXo ; i++){
				if (mask[i]==0){
					nmask++;
				}
			}
			logmsg("masked=%d (%.2f%%)",nmask,100*nmask/(float)fftXo);
			cpgopen("1/xs");
			cpgsvp(0.1,0.9,0.1,0.9);
			cpgswin(0,fftXo,0,T*10);
			cpgbox("ABN",0,0,"ABN",0,0);
			cpgline(fftXo,fx,pseries);
			cpgsci(2);
			cpgline(fftXo,fx,fT);
			cpgclos();

		}


		//		exit(1);

		start_clock(t_spec);

		//FILE* ff=fopen("plot","w");
#pragma omp parallel for shared(fftd,ch_mean,ch_fft_n,ch_fft_p)
		for (uint64_t j = 0; j < nchan ; j++){
			float var = ch_var[offset+j];
			float m=sqrt(var*fftXo/2.0);
			float T = sqrt(var*fftXo)*3;
			uint64_t n=0;
			float p=0;
			float complex *fftch = fftd + fftXo*j;
			for(uint64_t i = 1; i < fftXo; i++){
				if (camp(fftch[i]) > T) {
					n++;
					p+=camp(fftch[i]);
				}
				//	 if(j==512)fprintf(ff,"%f ",camp(fftch[i]));
				if(mask[i]==0){
					fftch[i]=m*(mjk_rand_gauss(random) + I*mjk_rand_gauss(random)); 
				}
				//	 if(j==512)fprintf(ff,"%f\n",camp(fftch[i]));
			}
			ch_fft_n[offset+j]=n;
			ch_fft_p[offset+j]=p;
		}
		// fclose(ff);

		logmsg("iFFT");
		fftwf_execute(ifft_plan);

#pragma omp parallel for schedule(static,2) shared(buffer,clean)
		for (uint64_t j = 0; j < nchan ; j+=BLK){
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				for (uint64_t k = 0; k < BLK ; k++){
					clean[j+k][i]/=(float)fftX;
					buffer[i*nchan+j+k] = round(clean[j+k][i]);
				}
			}

			if(j==512){
				cpgopen("2/xs");
				cpgsvp(0.1,0.9,0.1,0.9);
				cpgswin(0,fftX,ch_mean[j]-sqrt(ch_var[j])*10,ch_mean[j]+sqrt(ch_var[j])*10);
				cpgbox("ABN",0,0,"ABN",0,0);
				cpgline(fftX,X,data[j]);
				cpgsci(2);
				cpgline(fftX,X,clean[j]);
				cpgclos();

			}

		}
		fwrite(buffer,1,bytes_per_block,of);


		for (uint64_t i = 0; i < nsamp_per_block ; i++){
			tseries[i]=0;
		}

		tmean=0;
		tvar=0;
		max=0;
		min=1e99;
		//#pragma omp parallel for shared(clean,tseries)
		// NOT THREAD SAFE
		for (uint64_t j = 0; j < nchan ; j++){
			tmean+=ch_mean[offset+j];
			tvar+=ch_var[offset+j];
			for (uint64_t i = 0; i < nsamp_per_block ; i++){
				tseries[i]+=clean[j][i];
				if(clean[j][i]>max)max=clean[j][i];
				if(clean[j][i]<min)min=clean[j][i];
			}
		}
		ss=0;
		mm=0;
		for (uint64_t i = 0; i < nsamp_per_block ; i++){
			float x=tseries[i]-tmean;
			mm+=tseries[i];
			ss+=x*x;
		}
		rvar=ss/(float)nsamp_per_block;
		logmsg("var=%g tvar=%g",ss/(float)nsamp_per_block,tvar);
		logmsg("mean=%g tmean=%g",mm/(float)nsamp_per_block,tmean);
		cpgopen("5/xs");
		cpgsvp(0.1,0.9,0.1,0.9);
		cpgswin(0,fftX,tmean-sqrt(tvar)*30,tmean+sqrt(tvar)*30);
		cpgbox("ABN",0,0,"ABN",0,0);
		cpgline(fftX,X,tseries);
		cpgsci(2);
		cpgclos();
		tr[0] = 0.0 ;
		tr[1] = 1;
		tr[2] = 0;
		tr[3] = 0.5;
		tr[4] = 0;
		tr[5] = 1;

		logmsg("max=%g min=%g",max,min);

		cpgopen("6/xs");
		cpgsvp(0.1,0.9,0.1,0.9);
		cpgswin(0,nsamp_per_block,0,nchan);
		cpgbox("ABN",0,0,"ABN",0,0);
		cpggray(*clean,nsamp_per_block,nchan,1,nsamp_per_block,1,nchan,tmean/(float)nchan+sqrt(rvar/(float)nchan),tmean/(float)nchan-sqrt(rvar/(float)nchan),tr);
		cpgclos();




		stop_clock(t_spec);
		for (uint64_t j = 0; j < nchan ; j++){
			float mean=ch_mean[offset+j];
			if (mean > max_mean)max_mean=mean;
			if (mean < min_mean)min_mean=mean;

			float var=ch_var[offset+j];
			if (var > max_var)max_var=var;
			if (var < min_var)min_var=var;
			float fft_n=ch_fft_n[offset+j];
			if (fft_n > max_fft_n)max_fft_n=fft_n;
			if (fft_n < min_fft_n)min_fft_n=fft_n;
			float fft_p=ch_fft_p[offset+j];
			if (fft_p > max_fft_p)max_fft_p=fft_p;
			if (fft_p < min_fft_p)min_fft_p=fft_p;

		}
	}
	stop_clock(t_all);

	fclose(of);

	logmsg("T(all)  = %.2lfs",read_clock(t_all));
	logmsg("T(read) = %.2lfs",read_clock(t_read));
	logmsg("T(trans)= %.2lfs",read_clock(t_trns));
	logmsg("T(fft)  = %.2lfs",read_clock(t_fft));
	logmsg("T(fan)  = %.2lfs",read_clock(t_spec));
	logmsg("T(rms)  = %.2lfs",read_clock(t_rms));
	logmsg("T(rest) = %.2lfs",read_clock(t_all)-read_clock(t_read)-read_clock(t_trns)-read_clock(t_rms)-read_clock(t_fft)-read_clock(t_spec));


	tr[0] = -tsamp*nsamp_per_block*0.5;
	tr[2] = tsamp*nsamp_per_block;
	tr[1] = 0;
	tr[3] = 0.5;
	tr[5] = 0;
	tr[4] = 1;




	cpgopen("1/xs");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_mean,nchan,nblocks,1,nchan,1,nblocks,max_mean,min_mean,tr);
	cpgclos();

	cpgopen("2/xs");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_var,nchan,nblocks,1,nchan,1,nblocks,max_var,min_var,tr);
	cpgclos();

	cpgopen("3/xs");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_fft_n,nchan,nblocks,1,nchan,1,nblocks,max_fft_n,min_fft_n,tr);
	cpgclos();

	cpgopen("4/xs");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_fft_p,nchan,nblocks,1,nchan,1,nblocks,max_fft_p,min_fft_p,tr);
	cpgclos();



	cpgopen("mean.ps/vcps");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_mean,nchan,nblocks,1,nchan,1,nblocks,max_mean,min_mean,tr);
	cpgclos();

	cpgopen("var.ps/vcps");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_var,nchan,nblocks,1,nchan,1,nblocks,max_var,min_var,tr);
	cpgclos();


	cpgopen("fft_n.ps/vcps");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_fft_n,nchan,nblocks,1,nchan,1,nblocks,max_fft_n,min_fft_n,tr);
	cpgclos();

	cpgopen("fft_p.ps/vcps");
	cpgsvp(0.1,0.9,0.1,0.9);
	cpgswin(0,nblocks*tsamp*nsamp_per_block,0,nchan);
	cpgbox("ABN",600,10,"ABN",100,1);
	cpggray(ch_fft_p,nchan,nblocks,1,nchan,1,nblocks,max_fft_p,min_fft_p,tr);
	cpgclos();



	fclose(f);
	free(buffer);
	free_2df(data);

	return 0;

}


