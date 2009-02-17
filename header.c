/*
   header - show header parameters in a data file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "sigproc.h"
int wapp_header_size, wapp_incfile_length;
int nbins;
	double period;
main(int argc, char *argv[]) 
{
	FILE *fileptr;
	char filename[80],*telescope,*backend,*datatype,message[80],unit[16];
	int i,j,year,month,day,check,rah,ram,ded,dem;
	double ras,des,frac,tobs;
	char sra[6],sde[6],decsign;
	int raw,uth,utm,uts;
	long long numsamps,datasize,headersize;

	int writeobsdbline;

	fileptr=stdin;
	strcpy(filename,"stdin");
	strcpy(rawdatafile,"stdin");
	pulsarcentric=barycentric=0;

	writeobsdbline=0;

	if (argc>1) {
		print_version(argv[0],argv[1]);
		if (help_required(argv[1])) {
			header_help();
			exit(0);
		} else if (file_exists(argv[1])) {
			strcpy(filename,argv[1]);
			fileptr=open_file(filename,"rb");
		} else if (!file_exists(argv[1]) && (strncmp(argv[1],"-",1) !=0)) {
			sprintf(message,"Data file: %s not found...\n",argv[1]);
			error_message(message);
			exit(1);
		}
	}


	if (!(headersize=read_header(fileptr))) {
		rewind(fileptr);
		if ((raw=typeof_inputdata(fileptr,filename))) {
			data_type=0;
			switch (raw) {
				case 1:
					headersize=32768;
					break;
				case 5:
					headersize=32768;
					break;
				case 3:
					headersize=wapp_header_size+wapp_incfile_length;
					break;
				default:
					break;
			}
		} else {
			error_message("could not read header parameters!");
			exit(1);
		}
	}

	/* attempt to find number of bytes of data and number of samples */
	if (!strings_equal(filename,"stdin")) {
		datasize=sizeof_file(filename)-headersize;
		numsamps=nsamples(filename,headersize,nbits,nifs,nchans);
	} else if (!strings_equal(rawdatafile,"stdin")) {
		datasize=sizeof_file(rawdatafile)-headersize;
		numsamps=nsamples(rawdatafile,headersize,nbits,nifs,nchans);
	} else {
		datasize=numsamps=0;
	}

	telescope=telescope_name(telescope_id);
	backend=backend_name(machine_id);
	datatype=data_category(data_type);

	if (argc>2) {
		check=1;
		i=2;
	} else if ((argc>1) && strings_equal(filename,"stdin")) {
		check=1;
		i=1;
	} else {
		check=0;
	}

	angle_split(src_raj,&rah,&ram,&ras);
	if (ras<10.0) {
		sprintf(sra,"0%.1f",ras);
	} else {
		sprintf(sra,"%.1f",ras);
	}

	angle_split(src_dej,&ded,&dem,&des);
	if (src_dej > 0.0) 
		decsign = '+';
	else 
		decsign = '-';
	if (des<10.0) {
		sprintf(sde,"0%.1f",des);
	} else {
		sprintf(sde,"%.1f",des);
	}

	cal(tstart,&year,&month,&day);

	if (check) {
		/* check command-line parameters */ 
		while (i<argc) {
			if (strings_equal(argv[i],"-telescope")) {
				puts(telescope);
			} else if (strings_equal(argv[i],"-obsdb")) {
				writeobsdbline=1;
			} else if (strings_equal(argv[i],"-machine")) {
				puts(backend);
			} else if (strings_equal(argv[i],"-source_name")) {
				puts(source_name);
			} else if (strings_equal(argv[i],"-scan_number")) {
				puti(scan_number);
			} else if (strings_equal(argv[i],"-datatype")) {
				puts(datatype);
			} else if (strings_equal(argv[i],"-frame")) {
				if (pulsarcentric) 
					puts("pulsarcentric");
				else if (barycentric) 
					puts("barycentric");
				else 
					puts("topocentric");
			} else if (strings_equal(argv[i],"-radec")) {
				printf("%02d:%02d:%s ",rah,ram,sra);
				printf("%c%02d:%02d:%s\n",decsign,abs(ded),dem,sde);
			} else if (strings_equal(argv[i],"-barycentric")) {
				puti(barycentric);
			} else if (strings_equal(argv[i],"-pulsarcentric")) {
				puti(pulsarcentric);
			} else if (strings_equal(argv[i],"-data_type")) {
				puti(data_type);
			} else if (strings_equal(argv[i],"-headersize")) {
				printf("%d\n",headersize);
			} else if (strings_equal(argv[i],"-datasize")) {
				printf("%lld\n",datasize);
			} else if (strings_equal(argv[i],"-nsamples")) {
				printf("%lld\n",numsamps);
			} else if (strings_equal(argv[i],"-tobs")) {
				printf("%f\n",(double)numsamps*tsamp);
			} else if (strings_equal(argv[i],"-az_start")) {
				printf("%f\n",az_start);
			} else if (strings_equal(argv[i],"-za_start")) {
				printf("%f\n",za_start);
			} else if (strings_equal(argv[i],"-fch1")) {
				printf("%.3f\n",fch1);
			} else if (strings_equal(argv[i],"-bandwidth")) {
				printf("%.3f\n",fabs(foff)*(double)nchans);
			} else if (strings_equal(argv[i],"-fmid")) {
				printf("%.3f\n",fch1+foff*nchans/2);
			} else if (strings_equal(argv[i],"-foff")) {
				printf("%f\n",foff);
			} else if (strings_equal(argv[i],"-refdm")||strings_equal(argv[i],"-dm")) {
				printf("%f\n",refdm);
			} else if (strings_equal(argv[i],"-nchans")) {
				printf("%d\n",nchans);
			} else if (strings_equal(argv[i],"-tstart")) {
				printf("%.12f\n",tstart);
			} else if (strings_equal(argv[i],"-frequencies")) {
				for (j=0; j<nchans; j++) printf("%f\n",frequency_table[j]);
			} else if (strings_equal(argv[i],"-mjd")) {
				printf("%d\n",(int)floor(tstart));
			} else if (strings_equal(argv[i],"-date")) {
				printf("%4d/%02d/%02d\n",year,month,day);
			} else if (strings_equal(argv[i],"-utstart")) {
				frac=tstart-floor(tstart);
				uth=(int) floor(24.0*frac);
				frac-=(double)uth/24.0;
				utm=(int) floor(1440.0*frac);
				frac-=(double)utm/1440.0;
				uts=(int) floor(86400.0*frac);
				printf("%02d:%02d:%02d\n",uth,utm,uts);
			} else if (strings_equal(argv[i],"-tsamp")) {
				printf("%.5f\n",tsamp*1.0e6);
			} else if (strings_equal(argv[i],"-nbits")) {
				printf("%d\n",nbits);
			} else if (strings_equal(argv[i],"-nifs")) {
				printf("%d\n",nifs);
			} else if (strings_equal(argv[i],"-src_raj")) {
				printf("%02d:%02d:%s\n",rah,ram,sra);
			} else if (strings_equal(argv[i],"-src_dej")) {
				printf("%c%02d:%02d:%s\n",decsign,abs(ded),dem,sde);
			} else {
				header_help();
				sprintf(message,"unknown argument (%s) passed to header",argv[i]);
				error_message(message);
			}
			i++;
		}
		/* if we are doing a obs line do this... otherwise continue normaly
		 * MK 2006, for the MM survey bookkeeping.
		 */


		if(writeobsdbline){
			printf("%s ",source_name);
                        printf("%3.3lf %3.3lf ",gal_l,gal_b);
			printf("%s ",filename);
                        printf("%5.3lf %6.6lf ",header_tobs,tstart);
			printf("%6.3lf %6.3lf ",src_raj,src_dej);
                        printf("%5.3lf %3.3lf ",raw_fch1,raw_foff);
			printf("%d %d %5.3lf",nbeams,nchans,tsamp*1000);
			printf("\n");

		}

		exit(0);
	}


	/* no command-line flags were specified - display full output */

	printf("Data file                        : %s\n",filename);
	printf("Header size (bytes)              : %d\n",headersize);
	if (datasize) 
		printf("Data size (bytes)                : %lld\n",datasize);
	if (pulsarcentric) 
		printf("Data type                        : %s (pulsarcentric)\n",datatype);
	else if (barycentric) 
		printf("Data type                        : %s (barycentric)\n",datatype);
	else
		printf("Data type                        : %s (topocentric)\n",datatype);

	printf("Telescope                        : %s\n",telescope);
	printf("Datataking Machine               : %s\n",backend);


	if (!strings_equal(source_name,"")) 
		printf("Source Name                      : %s\n",source_name);
	if (src_raj != 0.0) 
		printf("Source RA (J2000)                : %02d:%02d:%s\n",rah,ram,sra);
	if (src_dej != 0.0)
		printf("Source DEC (J2000)               : %c%02d:%02d:%s\n",
				decsign,abs(ded),dem,sde);
	if ((az_start != 0.0) && (az_start != -1.0))
		printf("Start AZ (deg)                   : %f\n",az_start);
	if ((za_start != 0.0) && (za_start != -1.0))
		printf("Start ZA (deg)                   : %f\n",za_start);

	switch (data_type) {
		case 0:
		case 1:
			if ((fch1==0.0) && (foff==0.0)) {
				printf("Highest frequency channel (MHz)  : %f\n",
						frequency_table[0]);
				printf("Lowest frequency channel  (MHz)  : %f\n",
						frequency_table[nchans-1]);
			} else {
				printf("Frequency of channel 1 (MHz)     : %f\n",fch1);
				printf("Channel bandwidth      (MHz)     : %f\n",foff); 
				printf("Number of channels               : %d\n",nchans);
				printf("Number of beams                  : %d\n",nbeams);
				printf("Beam number                      : %d\n",ibeam); 
			}
			break;
		case 2:
			nchans=1;
			printf("Reference DM (pc/cc)             : %f\n",refdm);
			printf("Reference frequency    (MHz)     : %f\n",fch1);
			break;
		case 3:
			if (refdm > 0.0)
				printf("Reference DM (pc/cc)             : %f\n",refdm);
			printf("Frequency of channel 1 (MHz)     : %f\n",fch1);
			printf("Channel bandwidth      (MHz)     : %f\n",foff); 
			printf("Number of channels               : %d\n",nchans);
			printf("Number of phase bins             : %d\n",nbins);
			printf("Folding period  (s)              : %.12f\n",period);
			break;
		case 6:
			printf("Reference DM (pc/cc)             : %f\n",refdm);
			printf("Frequency of channel 1 (MHz)     : %f\n",fch1);
			printf("Channel bandwidth      (MHz)     : %f\n",foff); 
			printf("Number of channels               : %d\n",nchans);
			break;
	}

	printf("Time stamp of first sample (MJD) : %.12f\n",tstart);
	printf("Gregorian date (YYYY/MM/DD)      : %4d/%02d/%02d\n",year,month,day);

	if (data_type != 3) 
		printf("Sample time (us)                 : %.5f\n",tsamp*1.0e6); 

	if (datasize && data_type != 3) {
		printf("Number of samples                : %lld\n",numsamps);
		tobs=(double)numsamps*tsamp;
		strcpy(unit,"(seconds)   ");
		if (tobs>60.0) {
			tobs/=60.0;
			strcpy(unit,"(minutes)   ");
			if (tobs>60.0) {
				tobs/=60.0;
				strcpy(unit,"(hours)     ");
				if (tobs>24.0) {
					tobs/=24.0;
					strcpy(unit,"(days)      ");
				}
			}
		}
		printf("Observation length %s  : %.1f\n",unit,tobs);
	}
	printf("Number of bits per sample        : %d\n",nbits);
	printf("Number of IFs                    : %d\n",nifs);
	exit(0);
}
