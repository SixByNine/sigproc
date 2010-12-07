#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/*
   DECIMATE  - decimate filterbank data by adding channels and/or time samples
   */

#include "decimate.h"
#include <stdio.h>
#include <string.h>

unsigned long long hashdjb2(unsigned char *str, int len)
{
	unsigned long long hash_long = 5381;
	unsigned long long c;
	int i;

	i = 0;
	while (i < len){
		c = *str++;
		//hash_long = ((hash_long << 5) + hash_long) + c; /* hash * 33 + c */
		hash_long = hash_long * (unsigned long long)33 ^ c;
		i++;
	}
	return hash_long;
}

int main (int argc, char *argv[])
{
	int i,len, nc, headersize, headerless=0;
	char string[80];
	unsigned char dataSample[4096];
	int year,month,day,check,rah,ram,ded,dem;
	        double ras,des,frac,tobs;
		double design;
int blen;

	/* set up default global variables */
	headerless=0;
	input=stdin;
	strcpy(inpfile,"stdin");
	output=stdout;
	strcpy(outfile,"stdout");

	if (argc > 1) {
		/* check command-line parameters */ 
		print_version(argv[0],argv[1]);
		i=1;
		while (i<argc) {
			if (file_exists(argv[i])) {
				strcpy(inpfile,argv[i]);
				input=open_file(inpfile,"rb");
			} else {
				sprintf(string,"unknown argument (%s) passed to decimate",argv[i]);
				error_message(string);
			}
			i++;
		}
	}

	/* read in the header to establish what the input data are... */
	if ((headersize=read_header(input))) {
		//    if ( (nsamp > 0) && !strings_equal(inpfile,"stdin") ) {
		//	    nsamples=nsamples(inpfile,headersize,nbits,nifs,nchans)/nsamp;
		//  }
		nsamp = nsamples(inpfile,headersize,nbits,nifs,nchans);

		switch (data_type) {
			case 1:
				break;
			case 2:
				nchans=1;
				break;
			default:
				break;
		}
		/* check number of time samples to add */
	} else {
		error_message("input data file is of unknown origin!!!");
	}

	        fseek(input,headersize,SEEK_SET);
		                fseek(input,0,SEEK_CUR);


	fprintf(stderr,"First 10 bytes are:\n");
        len = fread(dataSample,1, 10, input);
	for(i = 0; i < len; i++){
		fprintf(stderr,"%02x ",(unsigned char)dataSample[i]);
	}
	fprintf(stderr,"\n");

	fseek(input,0,SEEK_SET);
	len = fread(dataSample,1, 4096, input);


	angle_split(src_raj,&rah,&ram,&ras);

	
	angle_split(src_dej,&ded,&dem,&des);
			

	cal(tstart,&year,&month,&day);
	frac=tstart-floor(tstart);
	int  uth=(int) floor(24.0*frac);
	frac-=(double)uth/24.0;
	int utm=(int) floor(1440.0*frac);
	frac-=(double)utm/1440.0;
	int uts=(int) floor(86400.0*frac);


	// write out the new header!

	fprintf(output,"<?xml version='1.0'?>\n");
	fprintf(output,"<psrxml version='1'>\n");
	fprintf(output,"\t<source_name>%s</source_name>\n",source_name);
	fprintf(output,"\t<day_of_observation units='MJD'>%d</day_of_observation>\n",(int)floor(tstart));
	fprintf(output,"\t<midnight_to_first_sample units='ns'>%lld</midnight_to_first_sample>\n",(unsigned long long)((tstart - floor(tstart))*8.64e13));
	fprintf(output,"\t<current_sample_interval units='us'>%lf</current_sample_interval>\n",tsamp*1e6);
	fprintf(output,"\t<number_of_samples>%d</number_of_samples>\n",nsamp);
	fprintf(output,"\t<actual_obs_time units='s'>%lf</actual_obs_time>\n",nsamp*tsamp);
	fprintf(output,"\t<centre_freq_first_channel units='MHz'>%lf</centre_freq_first_channel>\n",fch1);
	fprintf(output,"\t<channel_offset units='MHz'>%lf</channel_offset>\n",foff);
	fprintf(output,"\t<number_of_channels>%d</number_of_channels>\n",nchans);
	fprintf(output,"\t<utc>%d-%02d-%02dT%02d:%02d:%02dZ</utc>\n",year,month,day,uth,utm,uts);
	fprintf(output,"\t<reference_dm>%f</reference_dm>\n",refdm);

	fprintf(output,"\t<start_coordinate>\n");
	fprintf(output,"\t\t<coordinate>\n");
	fprintf(output,"\t\t\t<ra units='degrees'>%lf</ra>\n",rah*15.0 + ram/4.0+ras/240.0);
	if(src_dej < 0){
		design = -1;
	} else {
		design = 1;
	}
	fprintf(output,"\t\t\t<dec units='degrees'>%lf</dec>\n",design*(abs(ded) + dem/60.0 + des/3600.0));
	fprintf(output,"\t\t\t<position_epoch>J2000</position_epoch>\n");
	fprintf(output,"\t\t</coordinate>\n");
	fprintf(output,"\t</start_coordinate>\n");

	//fprintf(output,"\t<receiver href='' />\n",machine_id);
	if(ibeam == 0) ibeam = 1;
	fprintf(output,"\t<receiver_beam>%d</receiver_beam>\n",ibeam);
	if(nbeams == 0)nbeams = 1;
	fprintf(output,"\t<total_beams_recorded>%d</total_beams_recorded>\n",nbeams);
	//fprintf(output,"\t<backend href='' />\n",machine_id);
	fprintf(output,"\t<recorded_polarisations>II</recorded_polarisations>\n");
	if(telescope_id==4){
		fprintf(output,"\t<telescope>\n");
		fprintf(output,"\t\t<name>Parkes</name>\n");
		fprintf(output,"\t\t<x>-4554231.6</x>\n");
		fprintf(output,"\t\t<y>2816759.1</y>\n");
		fprintf(output,"\t\t<z>-3454036.1</z>\n");
		fprintf(output,"\t\t<sigproc_code>4</sigproc_code>\n");
		fprintf(output,"\t\t<tempo_code>7</tempo_code>\n");
		fprintf(output,"\t\t<pulsarhunter_code>PARKES</pulsarhunter_code>\n");
		fprintf(output,"\t</telescope>\n");
	} else{
		fprintf(output,"\t<telescope>\n");
		fprintf(output,"\t\t<name>%s</name>\n",telescope_name(telescope_id));
		fprintf(output,"\t\t<sigproc_code>%d</sigproc_code>\n",telescope_id);
		fprintf(output,"\t</telescope>\n");
	}
	fprintf(output,"\t<observing_programme>%s</observing_programme>\n",project);

	if(telescope_id==4){
		fprintf(output,"\t<telescope_identifying_string>Parkes</telescope_identifying_string>");
	}
	fprintf(output,"\t<data>\n");
	fprintf(output,"\t\t<filename>%s</filename>\n",inpfile);
	fprintf(output,"\t\t<data_uid type='djb2:4096'>%016llx</data_uid>\n",hashdjb2(dataSample,len));
	fprintf(output,"\t\t<endian>");
	if(nbits < 9){
		fprintf(output,"INDEPENDANT");
	} else {
		fprintf(stderr,"Endianness not known!\n");
	}
	fprintf(output,"</endian>\n");
	fprintf(output,"\t\t<header_length units='bytes'>%d</header_length>\n",headersize);

	blen=nchans;
	while(blen < 4096){
		blen*=2;
	}
	fprintf(output,"\t\t<block_size units='bytes'>%d</block_size>\n",blen);
	fprintf(output,"\t\t<block_header_length units='bytes'>0</block_header_length>\n");
	fprintf(output,"\t\t<bits_per_sample>%d</bits_per_sample>\n",nbits);
	fprintf(output,"\t\t<data_order>TFP</data_order>\n");
	fprintf(output,"\t\t<bit_order_first_sample_in>LSB</bit_order_first_sample_in>\n");
	if(nbits > 8){
                fprintf(output,"\t\t<signed>TRUE</signed>\n");
	} else {
		fprintf(output,"\t\t<signed>FALSE</signed>\n");
	}

	fprintf(output,"\t</data>\n");

	fprintf(output,"</psrxml>\n");

	fclose(output);
	exit(0);
}
