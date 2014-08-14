#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*
 * zerodm.c Apply "zerodm" filter to 8-bit data!
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "sigproc.h"
#include <string.h>
#include "mjk_random.h"


void zerodm_help(){
	printf("zerodm: apply 'zerodm' filter to 8-bit filterbank data\n");
	printf("\nzerodm file.fil\n");
}

void zerodm(char* block, unsigned int blocksize, mjk_rand_t* random);



int wapp_header_size, wapp_incfile_length;
int nbins;
double period;
main(int argc, char *argv[]) 
{
	FILE *fileptr, *outfile;
	char filename[1024],*telescope,*backend,*datatype,message[80],unit[16];
	int i,j,year,month,day,check,rah,ram,ded,dem;
	double ras,des,frac,tobs;
	char sra[6],sde[6],decsign;
	int raw,uth,utm,uts;
	long long numsamps,datasize,headersize;
	double readsec,skipsec;


	fileptr=stdin;
	outfile=stdout;
	strcpy(filename,"stdin");
	strcpy(rawdatafile,"stdin");
	pulsarcentric=barycentric=0;


	if (argc>1) {
		print_version(argv[0],argv[1]);
		if (help_required(argv[1])) {
			zerodm_help();
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
	char force_read=0;

	for (i = 0; i < argc ; i++){
		if (strcmp(argv[i],"-s")==0){
			skipsec=atof(argv[++i]);
		}
		if (strcmp(argv[i],"-r")==0){
			readsec=atof(argv[++i]);
		}
		if (strcmp(argv[i],"-f")==0){
		   force_read=1;
		}

	}

	rewind(fileptr);


	char* block_array;
	unsigned long long int count;
	unsigned long long int update_count;
	unsigned long long int bytes_per_sample=(unsigned long long int)(nchans*nbits)/8;
	unsigned long long int bytes_to_read=datasize;
	unsigned long long int blocksize = bytes_per_sample;
	unsigned long long int numblocks = bytes_to_read / blocksize;
	unsigned long long int update_size = (unsigned long long int) (10.0 * (bytes_per_sample/tsamp));

	fprintf(stderr,"\n\n==============\n");
	fprintf(stderr,"Copying header (%d bytes)\n",headersize);

	block_array = (char*)malloc(headersize);
	count = fread(block_array,1,headersize,fileptr);
	if ( count != headersize ){
		fprintf(stderr,"Error! Could not read header %d/%d\n",count,headersize);
		exit(1);
	}
	count = fwrite(block_array,1,headersize,outfile);
	if ( count != headersize ){
		fprintf(stderr,"Error! Could not write header\n");
		exit(1);
	}
	free(block_array);

	fprintf(stderr,"Allocating memory...\n");
	block_array = (char*)malloc(blocksize);

	fprintf(stderr,"Setting up random seed...\n");
	mjk_rand_t * random = mjk_rand_init(time(NULL));


	fprintf(stderr,"Applying zerodm filter...\n");
	count = 0;
	update_count = update_size;
	while ( count < bytes_to_read ) { 
		int read = fread(block_array,1,blocksize,fileptr);
		if ( read < 1 ) {
			fprintf(stderr,"Error! Could not read enough data\n");
			exit(2);
		}

		zerodm(block_array,blocksize,random);


		read = fwrite(block_array,1,read,outfile);
		if ( read < 1 ) {
			fprintf(stderr,"Error! Could not write enough data\n");
			exit(2);
		}
		count += read;
		update_count += read;
		if ( update_count >= update_size ){
			update_count=0;
			fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			fprintf(stderr,"\t% 8.1f s",tsamp*count/(float)bytes_per_sample);
		}
	}
	fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	fprintf(stderr,"\t% 8.1f s",tsamp*count/(float)bytes_per_sample);

	fprintf(stderr,"\nDone\n");

	free(block_array);

	exit(0);
}
void zerodm(char* block, unsigned int blocksize, mjk_rand_t* random){
	
	unsigned int ich = 0;
	int isign;
	int sum=0;
	for(ich=0; ich < blocksize; ich++){
		sum+=block[ich];
	}
	float mean = (float)sum/(float)blocksize;
	int isub = (int)mean;
	float fsub = mean - isub;
	for(ich=0; ich < blocksize; ich++){
		int x = block[ich];
		x -= isub;
		if(mjk_rand_double(random) < fsub){
			x -= 1;
		}
		// set new mean to 64
		x+=64;
		if(x < 0)x=0;
		if(x > 255)x=255;
		block[ich] = (unsigned char)x;
	}
}







