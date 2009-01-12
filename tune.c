/*
   TUNE
   MK 2005 -- Based on fold.c 
   This will try and tune using stack algorithm
 */
#include "fold.h"
#include "cpgplot.h"
#include <string.h>

#define SPEED_OF_LIGHT  299792458.0

void     quikgray(float*dat,int nx,int ny,int nxx);


struct tuneparam{
	float maxPer;
	float minPer;
	float incPer;
	float maxAcc;
	float minAcc;
	float incAcc;
	float maxJer;
	float minJer;
	float incJer;
};

int main (int argc, char *argv[])
{
	/* local variables */
	FILE *bestfile;
	FILE *jreaperfile;
	FILE *tunelistfile;
	double fftperiod,fftsnr,fftdm;
	double pfactor,newmjd=0.0;
	float sefd,sub_time;
	int i,opened_input=0,opened_output=0,headersize=0;
	char string[80],format[40],jreaperfilename[40];
	char stemname[80],susname[40];
	int j,b,n,psnum,hrmnum;
	int nsubints,nbins150;
	float **subintArray;
	float **subintEnds;
	float *fortranSubintArray; 
	float dmcurve[1000],dmcurvidx[1000];
	int ndms; /* The number of elements used in the dmcurve array */
	float realDM;
	float perFactor,perMaxMax,perMaxMin, realPer;
	float accFactor,accMaxMax,accMaxMin, realAcc;
	float jerFactor,jerMaxMax,jerMaxMin, realJer;
	int use_jerk,use_accn;
	struct tuneparam param;
	float *fPtr, *sumProfile;
	float s,sMax,sSquared,sSquMax,rms;
	float mean,variance,snr,snrMax,snrMaxMax;
	float sum;
	int width,widthMax,widthMaxMax,widthLimit;
	int houghposn;
	int usepgplot; /*flag to use pgplot or to write data out to file */
	int usequikgray,usequadratic;
	float *snrMap,*snrMapPtr;
	float *houghplot,*finalsubint;
	float acc_offset,acc_offparam;
	int writejreaper,writelogfile;
	double deltaT,truePeriod; 
	float max,min,maxMax,minMax;
	float tr[6];
	float x1,x2,y1,y2,xleft,xright,ytop,ybottom,xscale,yscale,scale;
	float pf,af,jf;
	/* set up default globals */

	fftperiod=fftdm=fftsnr=-1.0;
	usepgplot=baseline=ascii=multiple=1;
	npuls=binary=totalpower=accumulate=0;
	time_offset=acceleration=skip_time=read_time=0.0;
	use_jerk=use_accn=0;
	asciipol=psrfits=stream=headerless=npulses=0;
	phase_start=folding_period=dump_time=tsamp_user=0.0;
	phase_finish=pfactor=1.0;
	jyfactor=sefd=userbase=0.0;
	nbins=0; /* this will get set in the folding routine if not set by user */
	acc_offset=0.0;
	strcpy(polyco_file,"");
	psnum = 0;
	usequadratic=usequikgray=writejreaper = 0;
	maxMax = -1000;
	nsubints = 128;
	nbins = 64;
	baseline=0;
	pf=af=jf=1.0;
	strcpy(format,"/XSERV");
	ndms = 0;
	/* check the command line parameters */
	i=1;
	while (i<argc) {
		print_version(argv[0],argv[1]);
		if (strings_equal(argv[i],"-o")) {
			/* get and open file for output */
			i++;
			strcpy(outfile,argv[i]);
			output=fopen(outfile,"wb");
			opened_output=1;
		}else if (strings_equal(argv[i],"-format")) {
			/* get output format*/
			i++;
			strcpy(format,argv[i]);
		}else if (strings_equal(argv[i],"-jreaper")) {
			/* set writejreaper*/
			writejreaper=1;
			i++;
			strcpy(jreaperfilename,argv[i]);

		}else if (strings_equal(argv[i],"-quikgray")) {
			/* set quikgray*/
			usequikgray=1;

		}else if (strings_equal(argv[i],"-quadgray")) {
			/* set quadratic*/
			usequadratic=1;
		} else if (strings_equal(argv[i],"-useaccn")) {
			use_accn=1;
		} else if (strings_equal(argv[i],"-usejerk")) {
			use_jerk=1;
		}else if (strings_equal(argv[i],"logfile")) {
			/* set writejreaper*/
			writejreaper=1;
			i++;
			strcpy(jreaperfilename,argv[i]);


		}else if (strings_equal(argv[i],"-bestfile")) {
			/* load dmcurv and params from 'best' output */
			i++;
			strcpy(string,argv[i]);
			bestfile = fopen(string,"r");
			fscanf(bestfile,"%s",inpfile);
			input=open_file(inpfile,"rb");
			opened_input=1;


			fscanf(bestfile,"%lf %lf %lf %d\n",&fftperiod,&fftsnr,&fftdm,&hrmnum);
			folding_period = fftperiod;
			realDM=fftdm;
			ndms = 0;
			while(!feof(bestfile)){
				fscanf(bestfile,"%d %f %f\n",&j,dmcurvidx+ndms,dmcurve+ndms);
				ndms++;
			}

			fclose(bestfile);

		} else if (strings_equal(argv[i],"-m")) {
			multiple=atoi(argv[++i]);
		} else if (strings_equal(argv[i],"-p")) {
			/* get folding period */
			i++;
			if (file_exists(argv[i])) {
				strcpy(polyco_file,argv[i]);
				folding_period=-1.0;
			} else {
				folding_period=atof(argv[i]);
			}
		} else if (strings_equal(argv[i],"-dt")) {
			/* add a time offset in seconds to tstart */
			time_offset=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-mjd")) {
			/* change the start time completely! */
			newmjd=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-sk")) {
			/* skip the first skip_time seconds before folding */
			skip_time=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-re")) {
			/* read and fold only read_time seconds of data */
			read_time=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-a")) {
			/* get acceleration for folding */
			acc_offset=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-d")) {
			/* get dumptime or number of pulses for subintegrations */
			i++;
			if (strcspn(".",argv[i])) {
				npulses=atoi(argv[i]);
			} else {
				dump_time=atof(argv[i]);
			}
		} else if (strings_equal(argv[i],"-t")) {
			/* get user-supplied sampling time */
			i++;
			tsamp_user=atof(argv[i]);
		} else if (strings_equal(argv[i],"-j")) {
			/* get user-supplied Jansky calibration factor */
			jyfactor=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-s")) {
			/* get user-supplied SEFD */
			sefd=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-b")) {
			/* get user-supplied baseline */
			baseline=0;
			userbase=atof(argv[++i]);
		} else if (strings_equal(argv[i],"-f")) {
			/* get period multiplication factor */
			i++;
			pfactor=atof(argv[i]);
		} else if (strings_equal(argv[i],"-l")) {
			/* get leading phase of pulse */
			i++;
			phase_start=atof(argv[i]);
			if ( (phase_start < 0.0) || (phase_start > 1.0) ) 
				error_message("start pulse phase out of range!");
		} else if (strings_equal(argv[i],"-r")) {
			/* get trailing phase of pulse */
			i++;
			phase_finish=atof(argv[i]);
			if ( (phase_finish < 0.0) || (phase_finish > 1.0) ) 
				error_message("final pulse phase out of range!");
		} else if (strings_equal(argv[i],"-n")) {
			/* get number of bins */
			i++;
			nbins=atoi(argv[i]);
		} else if (strings_equal(argv[i],"-sub")) {
			/* get number of subints */
			i++;
			nsubints=atoi(argv[i]);

		} else if (strings_equal(argv[i],"-pf")) {
			/* get period multiplication factor */
			i++;
			pf=atof(argv[i]);
		} else if (strings_equal(argv[i],"-af")) {
			/* get period multiplication factor */
			i++;
			af=atof(argv[i]);
		} else if (strings_equal(argv[i],"-jf")) {
			/* get period multiplication factor */
			i++;
			jf=atof(argv[i]);
		} else if (strings_equal(argv[i],"-ascii")) {
			/* write data as ASCII numbers */
			ascii=1;
		} else if (strings_equal(argv[i],"-totalpower")) {
			/* sum polarizations 1+2 before writing */
			totalpower=1;
		} else if (strings_equal(argv[i],"-epn")) {
			/* write data in EPN format */
			ascii=0;
		} else if (strings_equal(argv[i],"-bin")) {
			/* write data in SIGPROC binary format */
			binary=1;
		} else if (strings_equal(argv[i],"-acc")) {
			/* write out accumulated pulse profiles in subints */
			accumulate=1;
		} else if (strings_equal(argv[i],"-asciipol")) {
			/* write data as ASCII numbers for Jim's polarization code */
			asciipol=1;
		} else if (strings_equal(argv[i],"-psrfits")) {
			/* write data in PSRFITS format */
			ascii=0;
			psrfits=1;
#ifndef PSRFITS
			error_message("-psrfits option not supported in this compilation...\nConsult the SIGPROC manual for further information about PSRFITS.");
#endif
		} else if (strings_equal(argv[i],"-stream")) {
			/* write data as ASCII streams */
			stream=1;
		} else if (strings_equal(argv[i],"-sub")) {
			/* shorthand for -nobaseline -stream -d x */
			stream=1;
			baseline=0;
			i++;
			if (strcspn(".",argv[i])) {
				npulses=atoi(argv[i]);
			} else {
				dump_time=atof(argv[i]);
			}
		} else if (strings_equal(argv[i],"-nobaseline")) {
			/* processing correlation functions so don't subtract baseline */
			baseline=0;
		} else if (file_exists(argv[i])) {
			/* get and open file for input */
			strcpy(inpfile,argv[i]);
			input=open_file(inpfile,"rb");
			opened_input=1;
		} else if (help_required(argv[i])) {
			tune_help();
			exit(0);
		} else {
			/* unknown argument passed down - stop! */
			fold_help();
			sprintf(string,"unknown argument (%s) passed to %s",argv[i],argv[0]);
			error_message(string);
		}
		i++;
	}

	/* get appropriate calibration factor from SEFD and baseline */
	if (sefd != 0.0 && userbase != 0.0) jyfactor=sefd/userbase;

	/* multiply folding period by user-supplied factor */
	if (folding_period != -1.0) folding_period*=pfactor;

	/* check start and end phase of pulse */
	if (phase_start >= phase_finish) 
		error_message("silly pulse phases selected!");

	/* check npulses versus dump_time */
	if (npulses < 0) error_message("npulses < 0!");
	if ((npulses > 0) && (dump_time > 0.0)) 
		error_message("can't have npulses AND dumptime defined!");

	/* check for folding period still set to zero - if so, look for polyco.dat */
	if (folding_period == 0.0) {
		strcpy(polyco_file,"polyco.dat");
		if (file_exists(polyco_file)) {
			folding_period=-1.0;
		} else {
			error_message("folding period not specified and no polyco.dat found!");
		}
	}

	if (!opened_input) {
		/* no input file selected, use standard input */
		input=stdin;
		strcpy(inpfile,"stdin");
	}

	/* read in the header parameters from the input stream */
	if (!(headersize=read_header(input))) 
		error_message("could not read header parameters!");

	/*if (acceleration != 0.0) {*/
	tobs=tsamp*(double)nsamples(inpfile,headersize,nbits,nifs,nchans);
	if (tobs <= 0.0) error_message("could not get sensible observation time");
	/*} */

	/* override the header */
	if (newmjd!=0.0) tstart=newmjd;

	if (!opened_output) {
		/* no output file selected, use standard output */
		output=stdout;
		strcpy(outfile,"stdout");
	}

	/* open the raw data file and establish its origin and header pars */
	switch(data_type) {
		case 1: 
		case 2:
		case 6:
			open_log("fold.monitor");
			/*   folded_profiles=fold_data();*/
			break;
		default:
			error_message("input data is of unknown origin!!!");
	}


	if(nbins > (folding_period/1000.0)/tsamp) nbins = (int)((folding_period/1000.0)/tsamp +0.0001);
	/* let people know we are about to start */

	while(tobs/nsubints < folding_period/1000.0){
		nsubints /= 2;
	}

	nbins150 = (int)(1.6*nbins);
	sumProfile = malloc(sizeof(float)*nbins150);


	param.minPer = -0.5/pf;
	param.maxPer =  0.5/pf;
	param.incPer =  0.01/pf;

	/*
	   param.minPer =  -0;
	   param.maxPer =  0;
	   param.incPer =  2/(double)nsubints;
	 */

	param.minAcc = -0.01/af;
	param.maxAcc =  0.01/af;
	param.incAcc =  0.0005/af;
	if(!use_accn){
		param.minAcc =  -0.;
		param.maxAcc =  0.;
		param.incAcc =  0.02/(double)nsubints;
	}


	param.minJer = -0.001/jf;
	param.maxJer =  0.001/jf;
	param.incJer =  0.0006/jf;

	if(!use_jerk){
		param.minJer =  0.000;
		param.maxJer =  0.000;
		param.incJer =  0.02/(double)nsubints;
	}



	acc_offparam = acc_offset * tobs * nbins/((nsubints/2)*(nsubints/2)*SPEED_OF_LIGHT);

	param.minAcc+=acc_offparam;
	param.maxAcc+=acc_offparam;


	deltaT = param.maxPer * nsubints * (folding_period/1000.0) / nbins;
	perMaxMax = folding_period  * (1.0 + deltaT/tobs);
	deltaT = param.minPer * nsubints * (folding_period/1000.0) / nbins;
	perMaxMin = folding_period  * (1.0 + deltaT/tobs);
	deltaT = param.incPer * nsubints * (folding_period/1000.0) / nbins;
	perFactor = folding_period  * (deltaT/tobs);


	deltaT = param.maxAcc * (nsubints/2)* (nsubints/2) * (folding_period/1000.0) / nbins;
	accMaxMax =     SPEED_OF_LIGHT * deltaT/(tobs/2.0 * folding_period/1000.0);
	deltaT = param.minAcc * (nsubints/2)* (nsubints/2)* (folding_period/1000.0) / nbins;
	accMaxMin =     SPEED_OF_LIGHT * deltaT/(tobs/2.0 * folding_period/1000.0);
	deltaT = param.incAcc *(nsubints/2)* (nsubints/2)* (folding_period/1000.0) / nbins;
	accFactor =     SPEED_OF_LIGHT * deltaT/(tobs/2.0 * folding_period/1000.0);

/*
	deltaT = param.maxAcc * (double)(nsubints* nsubints) * (folding_period/1000.0) / nbins / 4.0 ;
        accMaxMax =   - deltaT/(tobs/2.0 * (folding_period/1000.0)*(folding_period/1000.0));
        deltaT = param.minAcc * nsubints * nsubints* (folding_period/1000.0) / nbins;
        accMaxMin =    - deltaT/(tobs * (folding_period/1000.0)*(folding_period/1000.0));
 
        deltaT = param.incAcc * nsubints * nsubints* (folding_period/1000.0) / nbins;
        accFactor =   - deltaT/(tobs * (folding_period/1000.0)*(folding_period/1000.0));

*/

	deltaT = param.maxJer *nsubints* nsubints * nsubints * (folding_period/1000.0) / nbins;
	jerMaxMax = SPEED_OF_LIGHT * deltaT/(tobs*tobs*folding_period/1000.0);
	deltaT = param.minJer * nsubints* nsubints *nsubints * (folding_period/1000.0) / nbins;
	jerMaxMin =  SPEED_OF_LIGHT * deltaT/(tobs*tobs*folding_period/1000.0);
	deltaT = param.incJer * nsubints* nsubints *nsubints * (folding_period/1000.0) / nbins;
	jerFactor =  SPEED_OF_LIGHT * deltaT/(tobs*tobs*folding_period/1000.0);


	printf("Num Per: %d (%f,%f,%e) ms\n",(int)((param.maxPer - param.minPer)/param.incPer +1),perMaxMin,perMaxMax,perFactor);
	printf("Num Acc: %d (%f,%f,%f) m/s/s\n",(int)((param.maxAcc - param.minAcc)/param.incAcc +1),accMaxMin,accMaxMax,accFactor);
	printf("Num Jer: %d (%g,%g,%g) m/s/s/s\n",(int)((param.maxJer - param.minJer)/param.incJer +1),jerMaxMin,jerMaxMax,jerFactor);

	printf("Total Trials: %d\n",(int)((param.maxPer - param.minPer)/param.incPer +1)*(int)((param.maxAcc - param.minAcc)/param.incAcc +1)*(int)((param.maxJer - param.minJer)/param.incJer +1));


	skip_time = 0.0;
	sub_time = tobs/nsubints;
	subintArray = malloc(sizeof(float*)*nsubints);
	fortranSubintArray = malloc(sizeof(float)*nsubints*nbins);

	i = (int)(sub_time / (folding_period/1000.0));

	sub_time = (folding_period/1000.0) * i;

	read_time = sub_time;


	max = -10000.0f;
	min = 10000.0f;
	printf("Generating %d subints:\nTobs: %lf",nsubints,tobs);


	for(j = 0;j<nsubints;j++){	
		/*	printf("%4.3lf %4.3lf %4.2lf %4.2lf %4.2f\n",tobs,tsamp,read_time,skip_time,sub_time);*/
		subintArray[j]=fold_data();

		/* 	printf("%4.3lf %4.3lf %4.2lf %4.2lf\n",tobs,tsamp,read_time,skip_time);*/
		/* 	write_profiles(subintArray[j],nbins,nchans,nifs,output);*/

		for(i = 0; i < nbins; i++){
			fortranSubintArray[i+nbins*j] = subintArray[j][i]*subintArray[j][i];
			if(subintArray[j][i] < min)min = subintArray[j][i];
			if(subintArray[j][i] > max) max = subintArray[j][i];
		}
		rewind(input);
		read_header(input);
		skip_time = (j+1)*sub_time;
		read_time = sub_time;
	}
	for(j = 0;j<nsubints;j++){
		for(i = 0; i < nbins; i++){
			subintArray[j][i] = (subintArray[j][i]-min)/(max - min);	
			fortranSubintArray[i+nbins*j] = (fortranSubintArray[i+nbins*j]-min)/(max - min);	
		}
	}

	max = -10000.0f;
	min = 10000.0f;



	variance = 0.000498;
	mean = 0.0;
	snrMap = malloc(sizeof(float)*(int)((param.maxAcc - param.minAcc)/param.incAcc +1)*(int)((param.maxPer - param.minPer)/param.incPer +1));

	houghplot = malloc(sizeof(float)*nbins150*(int)((param.maxPer - param.minPer)/param.incPer +1));
	finalsubint = malloc(sizeof(float)*nbins*nsubints);

	for(j = 0; j < nbins150*(int)((param.maxPer - param.minPer)/param.incPer +1) ;j++) houghplot[j]=0.0f;

	snrMapPtr = snrMap;
	snrMaxMax = -10000;
	widthMaxMax = 0;

	accMaxMax = 0;
	jerMaxMax = 0;

	printf("Tuning period\n");
	jerFactor = param.minJer;

	while(jerFactor <= param.maxJer){

		accFactor = param.minAcc;
		while(accFactor <= param.maxAcc){
			houghposn = -1;
			perMaxMax = 0;
			perFactor = param.minPer;

			while(perFactor <= param.maxPer){

				houghposn++;
				/*printf("test3 %f / %f : %d / %d \n",accFactor,param.maxAcc,houghposn,(int)((param.maxAcc - param.minAcc)/param.incAcc +1));

				/*for(j = 0; j < nbins;j++) sumProfile[j]=0.0f;*/

				for(j = 0;j<nsubints;j++){
					b = j - nsubints/2;
					fPtr = subintArray[j]+(int)(perFactor * (j) + accFactor*b*b +jerFactor*b*b*b + 0.5);
					while((fPtr) < subintArray[j]) fPtr += nbins;		
					for(i = 0; i < nbins150; i++){
						while((fPtr + i) >= subintArray[j]+nbins){
							fPtr -= nbins;
						}
						houghplot[i+houghposn*nbins150] += *(fPtr+i);
						/*printf("%f , %f\n",*(fPtr+i),sumProfile[j]);*/	

					}

				}
				sum = 0;
				for(j = 0; j < nbins150;j++){
					houghplot[j+houghposn*nbins150] /= nsubints;
					sum += houghplot[j+houghposn*nbins150];
				}
				mean = sum / nbins150;
				snrMax = -100000;


				width = 1;
				while(width < nbins/2){
					s = 0;
					sMax = -100000;
					for(i = 0; i < width; i++){
						s += houghplot[i+houghposn*nbins150];
						s -= mean;
					}

					for(i = 0; i < nbins; i++){
						s -= houghplot[i+houghposn*nbins150];
						s += houghplot[i+width+houghposn*nbins150];
						if(s > sMax){
							sMax = s;
						}
					}

					snrMax = sMax/(float)sqrt(width);

					if(snrMax > snrMaxMax){
						snrMaxMax   = snrMax;
						perMaxMax   = perFactor;
						accMaxMax   = accFactor;
						jerMaxMax   = jerFactor;
						widthMaxMax = width;
						maxMax = -10000000.0f;
						minMax = 100000000.0f;
						for(i = 0; i < nbins150; i++){
							sumProfile[i] = houghplot[i+houghposn*nbins150];
							if(sumProfile[i] > maxMax) maxMax = sumProfile[i];
							if(sumProfile[i] < minMax) minMax = sumProfile[i];

						}
					}
					width *= 2;
				}


				perFactor += param.incPer;
				snrMapPtr++;
			}


			accFactor += param.incAcc;

		}
		jerFactor += param.incJer;
	}


	/* convert to a real SNR value! */

	for(i = 0; i < nbins150; i++){

		sumProfile[i] = (sumProfile[i] - minMax)/(maxMax - minMax);
	}


	n = (int)(((float)nbins / 6.0) + 0.5); /* The number of 'quiet' bins to find */	


	if(n < 10) n = (int)(((float)nbins / 2.0) + 0.5);
	if(n < 4) n = (nbins-widthMaxMax);
	s = 0.0;
	sSquared = 0.0;

	for(i = 0; i < n; i++){
		s += sumProfile[i];
		sSquared += sumProfile[i]*sumProfile[i]; 
	}
	sMax = s;
	sSquMax = sSquared;
	j = 0;
	for(i = 0; i < nbins; i++){
		s -= sumProfile[i];
		s += sumProfile[i+n];
		sSquared -= sumProfile[i]*sumProfile[i]; 
		sSquared += sumProfile[i+n]*sumProfile[i+n];

		if(s < sMax){
			sMax = s;
			sSquMax = sSquared;
			j = i;
		}
	}
	mean = (sMax / n);

	variance = sSquMax/n - mean*mean;
	s = 0.0;

	for(i = 0; i < widthMaxMax; i++){
		s += sumProfile[i]-mean;
	}
	sMax = s;
	j = 0;
	for(i = 0; i < nbins; i++){
		s +=  sumProfile[i+widthMaxMax] - mean;
		s -= sumProfile[i]-mean;
		if(s > sMax){
			sMax = s;
			j = i;
		}

	}




	snr = sMax/(sqrt(variance*widthMaxMax));

	printf("Using %d/%d 'Quiet bins' for snr\n",n,nbins);

	printf("Drawing plot\n");

	/***** Generate the final subint profile for the plot. *****/


	for(j = 0;j<nsubints;j++){
		b = j - nsubints/2;
		fPtr = subintArray[j] + (int)(perMaxMax * (j) +accMaxMax*b*b +jerMaxMax*b*b*b + 0.5);
		while((fPtr) < subintArray[j]) fPtr += nbins;
		for(i = 0; i < nbins; i++){
			while((fPtr + i) >= subintArray[j]+nbins){
				fPtr -= nbins;
			}
			if(usequadratic){
				finalsubint[i+j*nbins] += *(fPtr+i) * *(fPtr+i);
			}else{
				finalsubint[i+j*nbins] += *(fPtr+i);
			}
			/*printf("%f , %f\n",*(fPtr+i),sumProfile[j]);*/

		}

	}



	/***** BEGIN generating plot data and drawing ******/
	deltaT = perMaxMax * nsubints * (folding_period/1000.0) / nbins;
	realPer = folding_period  * (1.0 + deltaT/tobs);



	deltaT = accMaxMax *(nsubints/2) * (nsubints/2) * (folding_period/1000.0) / nbins;
	realAcc =     SPEED_OF_LIGHT * deltaT/(tobs/2.0 * folding_period/1000.0);

	
	/*realAcc =  - deltaT/(tobs * (folding_period/1000.0)*(folding_period/1000.0));
	realAcc = deltaT;	
*/
	deltaT = jerMaxMax *nsubints*nsubints* nsubints * (folding_period/1000.0) / nbins;
	realJer = SPEED_OF_LIGHT * deltaT/(tobs*tobs*folding_period/1000.0);



	if(writejreaper){

		jreaperfile = fopen(jreaperfilename,"wa");
		fprintf(jreaperfile,inpfile);	
		fprintf(jreaperfile,"\n%lf %lf %lf\n",src_raj,src_dej,tstart);
		fprintf(jreaperfile,"%lf %lf %lf\n",fftperiod,fftsnr,fftdm);
		fprintf(jreaperfile,"%f %f %f %d %f\n", realPer,realAcc,realJer,widthMaxMax,snr);
		fprintf(jreaperfile,"#GID# %s\n",source_name);

	}	


	/**** PGPLOT stuff for hough test ****/

	max = -10000.0f;
	min = 10000.0f;

	for(i = 0; i < nbins150 * houghposn; i++){
		if(max < houghplot[i]) max = houghplot[i];
		if(min > houghplot[i]) min = houghplot[i];
	}





	max -= min;
	for(i = 0; i < nbins150 * houghposn; i++){
		houghplot[i] -= min;
		if(houghplot[i] < 0.75 * max) houghplot[i] =0;
	}
	min = 0.75 * max;
	/*if(max > maxMax){
	  accMaxMax = accFactor;
	  jerMaxMax = jerFactor;
	  maxMax = max;
	  }*/
	sprintf(string,format,psnum);
	psnum++;


	cpgbeg(0, string, 1, 1);
	cpgsch(0.8);
	/* Draw Hough Plot */
	cpgsvp(0.05,0.45,0.425,0.55);
	cpgswin(0.0, (float)nbins, 0, (float)houghposn);

	cpgqwin( &x1, &x2, &y1, &y2 );

	xscale = ( x2 - x1 ) / nbins;
	yscale = ( y2 - y1 ) / houghposn;
	scale = ( xscale < yscale ) ? xscale : yscale;


	xleft   = 0.5f * ( x1 + x2 - nbins * scale );
	xright  = 0.5f * ( x1 + x2 + nbins * scale );
	ybottom = 0.5f * ( y1 + y2 - houghposn * scale );
	ytop    = 0.5f * ( y1 + y2 + houghposn * scale );


	tr[0] = xleft - 0.5f * scale;
	tr[1] = scale;
	tr[2] = 0.0f;
	tr[3] = ybottom - 0.5f * scale;
	tr[4] = 0.0f;
	tr[5] = scale;



	cpglab("Bin", "Period step","");

	cpggray(houghplot,nbins150,houghposn,1,nbins,1,houghposn,max,min,tr);

	cpgbox("ABCN",0.0,0,"ABCN",0.0,0);


	if(writejreaper){
		fprintf(jreaperfile,"#HOUGH# %d %d\n",nbins,houghposn);
		for(i = 0; i < houghposn; i++){
			for(j = 0; j < nbins; j++){
				fprintf(jreaperfile,"%2.2f ",houghplot[j+nbins150*i]);
			}
		}
	}




	/* Draw best profile */

	cpgsvp(0.55,0.95,0.1,0.325);
	cpgswin(0.0, (float)nbins, 0.0, 1.0);
	cpgbox("ABCSN",0.0,0,"ABCN",0.0,0);

	cpgmove(0.0,sumProfile[nbins-1]);
	for(i = 0; i < nbins; i++){
		cpgdraw(i,sumProfile[i]);
		cpgmove(i,sumProfile[i]);

	}

	if(writejreaper){
		fprintf(jreaperfile,"\n#BESTPROFILE# %d\n",nbins);
		for(i = 0; i < nbins; i++){
			fprintf(jreaperfile,"%1.2f ",sumProfile[i]-minMax);

		}
	}

	/* Draw origianal Profile */
	maxMax = -100000;
	minMax = 100000;

	for(j = 0; j < nbins; j++){
		sumProfile[j] = 0;
	}
	for(i = 0; i < nsubints; i++){
		for(j = 0; j < nbins; j++){
			sumProfile[j] = sumProfile[j] + subintArray[i][j];
		}
	}
	for(j = 0; j < nbins; j++){
		if(sumProfile[j] > maxMax)maxMax = sumProfile[j];
		if(sumProfile[j] < minMax)minMax = sumProfile[j];
	}
	cpgsvp(0.05,0.45,0.1,0.325);

	for(i = 0; i < nbins; i++){

		sumProfile[i] = (sumProfile[i] - minMax)/(maxMax - minMax);

	}

	cpgswin(0.0, (float)nbins, 0.0, 1);
	cpgbox("ABCSN",0.0,0,"ABCN",0.0,0);
	cpgmove(0.0,sumProfile[nbins-1]);
	for(i = 0; i < nbins; i++){
		cpgdraw(i,sumProfile[i]);
		cpgmove(i,sumProfile[i]);

	}

	if(writejreaper){
		fprintf(jreaperfile,"\n#ORIGPROFILE# %d\n",nbins);
		for(i = 0; i < nbins; i++){
			fprintf(jreaperfile,"%1.2f ",sumProfile[i]-minMax);

		}
	}


	/* Draw DM curve if we have one */

	if(ndms > 0){
		/* we have some dms in the dm curve, so draw it*/
		maxMax = -100000;
		minMax = 100000;

		for(i = 0; i < ndms; i++){
			if(dmcurve[i] > maxMax)maxMax = dmcurve[i];
			if(dmcurve[i] < minMax)minMax = dmcurve[i];
		}                                                                                                             

		cpgsvp(0.05,0.45,0.65,0.825);
		cpgswin(dmcurvidx[0], dmcurvidx[ndms-1], 0.0, maxMax-minMax);
		cpgbox("ABCSN",0.0,0,"ABCN",0.0,0);


		cpgmove(dmcurvidx[0],dmcurve[0]-minMax);
		for(i = 0; i < ndms; i++){
			cpgdraw(dmcurvidx[i],dmcurve[i]-minMax);
			cpgmove(dmcurvidx[i],dmcurve[i]-minMax);
			/*cpgpt1(dmcurvidx[i],dmcurve[i]-minMax,-1);*/

		}                                                                                                              
		if(writejreaper){
			fprintf(jreaperfile,"\n#DMCURVEIDX# %d\n",ndms);
			for(i = 0; i < ndms; i++){
				fprintf(jreaperfile,"%f ",dmcurvidx[i]);

			}
			fprintf(jreaperfile,"\n#DMCURVE# %d\n",ndms);
			for(i = 0; i < ndms; i++){
				fprintf(jreaperfile,"%f ",dmcurve[i]-minMax);
			}
		}


	}



	/* Draw subints */


	max = -10000.0f;
	min = 10000.0f;

	for(i = 0; i < nbins * nsubints; i++){
		if(max < finalsubint[i]) max = finalsubint[i];
		if(min > finalsubint[i]) min = finalsubint[i];
	}
	for(i = 0; i < nbins * nsubints; i++){
		finalsubint[i] = (finalsubint[i]-min)/(max-min);
	}
	/*max = (max - min)/2.0;
	  for(i = 0; i < nbins * nsubints; i++){
	  fortranSubintArray[i] = fortranSubintArray[i] - min;
	  if(fortranSubintArray[i] > max) fortranSubintArray[i] = max;
	  }*/

	cpgsvp(0.55,0.95,0.45,0.825);

	if(!usequikgray){
		cpgswin(0.0, (float)nbins, 0, (float)nsubints);
		cpgqwin( &x1, &x2, &y1, &y2 );

		xscale = ( x2 - x1 ) / nbins;
		yscale = ( y2 - y1 ) / nsubints;
		scale = ( xscale < yscale ) ? xscale : yscale;


		xleft   = 0.5f * ( x1 + x2 - nbins * scale );
		xright  = 0.5f * ( x1 + x2 + nbins * scale );
		ybottom = 0.5f * ( y1 + y2 - nsubints * scale );
		ytop    = 0.5f * ( y1 + y2 + nsubints * scale );


		tr[0] = xleft - 0.5f * scale;
		tr[1] = scale;
		tr[2] = 0.0f;
		tr[3] = ybottom - 0.5f * scale;
		tr[4] = 0.0f;
		tr[5] = scale;




		cpgswin(0.0,((float)nbins),0.0,(float)nsubints);

		cpggray(finalsubint,nbins,nsubints,1,nbins,1,nsubints,1.0,0.0,tr);
	}else{	
		cpgswin(0.0,((float)nbins),0.0,(float)nsubints);

		quikgray(finalsubint,nbins,nsubints,nbins);
	}
	cpgbox("ABCN",0.0,0,"ABCN",0.0,0);

	cpgsch(0.9);

	if(writejreaper){
		fprintf(jreaperfile,"\n#SUBINTS# %d %d\n",nbins,nsubints);
		for(i = 0; i < nbins*nsubints; i++){
			fprintf(jreaperfile,"%1.2f ",finalsubint[i]);

		}
	}



	/* draw title */

	cpgsvp(0.0,1.0,0.0,1.0);

	sprintf(string,"Tune Params: Period: %f, Freq: %f, Width: %d, SNR: %f",                 realPer,1000./realPer,widthMaxMax,snr);
	cpgswin(0.0,1.0,0.0,1.0);
	cpgtext(0.1,0.98,string);
	sprintf(string,"RA: %6.2lf DEC: %6.2lf Infile: %s MJD: %9.3lf" ,                        src_raj,src_dej,inpfile,tstart);
	cpgtext(0.1,0.95,string);

	sprintf(string,"Acc: %f, Jer: %f",realAcc,realJer);
	cpgtext(0.1,0.9,string);
	if(fftperiod > 0){

		sprintf(string,"Seek params: Period %lf, DM: %lf, SNR:%lf",fftperiod,fftdm,fftsnr);
		cpgtext(0.1,0.925,string);

	}	


	cpgsch(0.8);                        
	sprintf(string,"Sub Integrations");
	cpgtext(0.55,0.85,string);
	sprintf(string,"DM Curve");
	cpgtext(0.05,0.85,string);
	sprintf(string,"Hough Plane");
	cpgtext(0.05,0.575,string);
	sprintf(string,"Pre Tune profile");
	cpgtext(0.05,0.35,string);
	sprintf(string,"Post Tune profile");
	cpgtext(0.55,0.35,string);


	cpgend();


	if(writejreaper){
		fclose(jreaperfile);
	}

	deltaT = perMaxMax * nsubints * (folding_period/1000.0) / nbins;
	realPer = folding_period  * (1.0 + deltaT/tobs);



	deltaT = accMaxMax *(nsubints/2)* (nsubints/2) * (folding_period/1000.0) / nbins;
	realAcc =     SPEED_OF_LIGHT * deltaT/(tobs/2.0 * folding_period/1000.0);

	deltaT = jerMaxMax * nsubints*nsubints* nsubints  * (folding_period/1000.0) / nbins;
	realJer = SPEED_OF_LIGHT * deltaT/(tobs*tobs*folding_period/1000.0);

	/* hack - originally started by MK
	   - fiddled by SJ and DRL to get the output we want
	   - these first lines extract the directory file name from the entire dm name */
	if(writejreaper){
		i=0;
		for (j=0; j<strlen(inpfile); j++) {
			if (inpfile[j] == '_') i++;
			if (i == 3) break;
		}
		strncpy(stemname,inpfile,j);
		strncpy(susname,jreaperfilename,10);
		strcat(susname,".ps");
		tunelistfile = fopen("tune.lis","a");
		fprintf(tunelistfile,"%s %s %s %8.6lf %6.2lf %6.2lf %d %8.6f %6.2f %6.2f %4d %6.2f %6.2f\n",source_name,stemname,susname,fftperiod,fftsnr,fftdm,hrmnum,realPer,snr,realDM,widthMaxMax,realAcc,realJer);
		fclose(tunelistfile);	
	}


	printf("\n=========\nPeriod: %f, Acc: %f, Jerk: %f\n",realPer,realAcc,realJer);
	printf("Width: %d, SNR: %f\nDone.\n",widthMaxMax,snr);
	/* all done, update and close logfile */
	update_log("finished");
	close_log();
	i=0;
#ifdef PSRFITS
	if (psrfits) fits_close_file(fits,&i);
#endif






	free(snrMap);
	free(sumProfile);
	free(subintArray);
	free(finalsubint);
	exit(0);
}

void     quikgray(float*dat,int nx,int ny,int nxx)
{
	/*
	   Coarse grey-scale plot using PG plot routines
	   Assumed that viewport and window defined outside routine
	   nx and ny are the areas of the array display, the x dimension is
	   repeated out to nxx, for multiple cycles.
	 */

	int ksym[8]={0,1,20,21,2,3,17,18};
	float s=0.;
	float ss=0.;
	float smin=1.e30;
	float smax=-1.e30;
	float aa,rms,x,y;
	int i,ii,j,k;

	for (i=0;i<nx*ny;i=i+1)
	{
		aa=dat[i];
		s=s+aa;
		ss=ss+aa*aa;
		if(aa > smax)smax=aa;
		if(aa < smin)smin=aa;
	}

	s=s/(float)(nx*ny);
	rms=sqrt(ss/(float)(nx*ny)-s*s);
	if(s+7*rms > smax) rms=(smax-s)/7.0;
	printf("max: %f min %f mean: %f rms: %f\n",smax,smin,s,rms);
	for (j=0;j<ny;j=j+1)
	{
		ii=0;
		for (i=0;i<nxx;i=i+1)
		{
			ii=ii+1;
			if(ii==nx) ii=0;
			k=(int)((dat[(int)(j*nx+ii)]-s)/rms);
			if (k>7) k=7;
			x=(float)i+1;y=(float)j+1;
			if (k>0) cpgpt(1,&x,&y,ksym[k]);
		}
	}
}

