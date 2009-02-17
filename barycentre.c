/* 
   barycentre.c - refer a filterbank/timeseries file to the rest
   frame of the solar system barycentre by computing an appropriate
   set of polynomial coefficients for the data sampling interval,
   then keeping track of the difference between elapsed time (topo)
   versus the time computed from the coefficients (barycentric) and
   adding or removing time samples so that the two timescales stay
   fixed. Added samples are set to zero.

   Modification history:
   March 18, 2007, drl added -mypar option to override position in header
                   this is useful for observations where there may be a small
                   position offset between the telescope position and the true
                   pulsar position and the resulting difference can cause
                   the barycentric time series to drift. Also included refdm
                   in the calculation.
*/
#include "dedisperse.h"
char polyco_filename[80];
double fcent;
/* subroutine to call TEMPO to calculate a polyco.bar file for barycentering */
char *make_polycofile(char ra[],char dec[],char topo[], char site,
			    double mjdtopo, double tsamp)
{
  FILE *resid2,*parfile, *tzfile;
  float junk;
  double mjdbary;
  char *polycofilename;
  polycofilename=(char *) malloc(80);

  parfile=fopen("tssb.par","w");
  fprintf(parfile,"PSR 0000+00\n");
  fprintf(parfile,"RAJ %s\n",ra);
  fprintf(parfile,"DECJ %s\n",dec);
  fprintf(parfile,"F0 1.0\n");
  fprintf(parfile,"DM %f\n",refdm);
  fprintf(parfile,"PEPOCH %s\n",topo);
  fclose(parfile);
  tzfile=fopen("tz.in","w");
  fprintf(tzfile,"%c    2  30  9 %lf\n",site,fcent);
  fprintf(tzfile,"\n \n");
  fprintf(tzfile,"0000+00 60 9 12 %lf\n",fcent);
  fclose(tzfile);
  tzfile=fopen("runtempo.csh","w");
  fprintf(tzfile,"#!/bin/csh\n",site);
  fprintf(tzfile,"tempo -z -f tssb.par << EOD\n");
  fprintf(tzfile,"%f %f\n",mjdtopo-1.0,mjdtopo+1.0);
  fprintf(tzfile,"EOD");
  fclose(tzfile);
  system("csh runtempo.csh > /dev/null");
  system("mv polyco.dat polyco.bar");
  system("rm -f tssb.par tz.in tz.tmp");
  system("rm -f fort.22 tempo.lis runtempo.csh ");
  strcpy(polycofilename,"polyco.bar");
  return(polycofilename);
}

/* subroutine to call TEMPO to calculate the barycentric MJD */
double barycentric_time(char ra[],char dec[],char topo[], char site,
			    double mjdtopo)
{
  FILE *resid2,*parfile, *timfile;
  float junk;
  double mjdbary;
  parfile=fopen("tssb.par","w");
  fprintf(parfile,"PSR 0000+00\n");
  fprintf(parfile,"RAJ %s\n",ra);
  fprintf(parfile,"DECJ %s\n",dec);
  fprintf(parfile,"F0 1.0\n");
  fprintf(parfile,"DM 0.0\n");
  fprintf(parfile,"PEPOCH %s\n",topo);
  fclose(parfile);
  timfile=fopen("tssb.tim","w");
  fprintf(timfile,"%c    0  0000+00 %.3lf %.13f     1.00\n",site,fcent,mjdtopo);
  fprintf(timfile,"%c    0  0000+00 %.3lf %.13f     1.00\n",site,fcent,mjdtopo);
  fprintf(timfile,"%c    0  0000+00 %.3lf %.13f     1.00\n",site,fcent,mjdtopo);
  fclose(timfile);
  system("tempo tssb.tim > /dev/null");
  resid2=fopen("resid2.tmp","r");
  fread(&junk,4,1,resid2);
  fread(&mjdbary,8,1,resid2);
  fclose(resid2);
  system("rm -f tssb.tim tssb.par 0000+00.par");
  system("rm -f matrix.tmp tempo.lis resid2.tmp ");
  return(mjdbary);
}

char inpfile[80], outfile[80];
main (int argc, char *argv[]) 
{
  int drop,add,i,j,n,ntim,headersize,rah,ram,ded,dem;
  int ndropped=0,nadded=0,nbytes_per_sample,verbose=0,mypolyco=0;
  unsigned char *rawdata, *dummy;
  double mjd, elapsed_time, barycentre_time;
  double ras,des,mjdbary;
  char ra[80], dec[80], topo[80], sra[6], sde[6], site;
  char message[80], myparfile[80], line[80], key[80];
  struct POLYCO polyco;
  FILE *polycofile,*parfile;

  if (argc<2 || help_required(argv[1])) {
    barycentre_help();
    exit(0);
  }
  print_version(argv[0],argv[1]);
  if (!file_exists(argv[1]))
    error_message("input file does not exist!");

  strcpy(inpfile,argv[1]);
  input=open_file(inpfile,"r");
  strcpy(outfile,"stdout");
  output=stdout;
  strcpy(myparfile,"");

  i=2;
  while (i<argc) {
    if (strings_equal(argv[i],"-verbose")) 
      verbose=1;
    if (strings_equal(argv[i],"-mypolyco")) 
      mypolyco=1;
    if (strings_equal(argv[i],"-myparfile")) 
      strcpy(myparfile,argv[++i]);
    i++;
  }

  if ((headersize=read_header(input))) {

    /* calculate centre frequency for use in TEMPO files */
    fcent=fch1+(double)(nchans/2)*foff;
    /* parse the header parameters for RA */
    angle_split(src_raj,&rah,&ram,&ras);
    if (ras<10.0) {
      sprintf(sra,"0%.1f",ras);
    } else {
      sprintf(sra,"%.1f",ras);
    }
    sprintf(ra,"%02d:%02d:%s",rah,ram,sra);
    /* parse the header parameters for DEC */
    angle_split(src_dej,&ded,&dem,&des);
    if (des<10.0) {
      sprintf(sde,"0%.1f",des);
    } else {
      sprintf(sde,"%.1f",des);
    }
    sprintf(dec,"%02d:%02d:%s",ded,dem,sde);
    /* now call TEMPO to calculate the barycentric MJD */
    sprintf(topo,"%.12f",tstart);
    site=tempo_site(telescope_id);
    if (verbose) 
      fprintf(stderr,"Telescope: %s TEMPO site code: %c\n",
	      telescope_name(telescope_id),site);
    if (!strings_equal(myparfile,"")) {
      parfile=open_file(myparfile,"r");
      while (fgets(line,80,parfile) != NULL) {
	strcpy(key,strtok(line," "));
	if (strings_equal(key,"RAJ"))
	  strcpy(ra,strtok(NULL," "));
	if (strings_equal(key,"DECJ"))
	  strcpy(dec,strtok(NULL," "));
      }
    }

    if (mypolyco) {
      strcpy(polyco_filename,"polyco.bar");
      if (verbose)
      fprintf(stderr,"Using barycentric polyco file: %s\n",polyco_filename);
    } else {
      strcpy(polyco_filename,make_polycofile(ra,dec,topo,site,tstart,tsamp));
      if (verbose)
      fprintf(stderr,"Created barycentric polyco file: %s\n",polyco_filename);
    }
      
    barycentric=1;
    mjdbary=barycentric_time(ra,dec,topo,site,tstart);
    if (verbose) {
      fprintf(stderr,"Topocentric MJD %.12f\n",tstart);
      fprintf(stderr,"Barycentric MJD %.12f\n",mjdbary);
    }
    /* write out header with barycentric MJD if required */
    send_string("HEADER_START");
    send_int("telescope_id",telescope_id); 
    send_int("machine_id",machine_id);
    send_coords(src_raj,src_dej,az_start,za_start);
    send_int("data_type",data_type);
    send_int("barycentric",1);
    send_int("pulsarcentric",0);
    if (nchans==1) send_double("refdm",refdm);
    if (fch1 == 0.0) 
      send_double("fch1",frequency_table[0]);
    else
      send_double("fch1",fch1);
    send_int("nchans",nchans);
    if (nchans>1) send_double("foff",foff);
    send_int("nbits",nbits);  
    send_int("nifs",nifs);  
    send_double ("tstart",mjdbary); 
    send_double("tsamp",tsamp);
    send_int("nbeams",nbeams);
    send_int("ibeam",ibeam);
    send_string("HEADER_END");
    open_log("barycentre.monitor");

    ntim=nsamples(inpfile,headersize,nbits,nifs,nchans);
    mjd=tstart;
    polycofile=open_file(polyco_filename,"r");
    if (!read_polycoset(polycofile,&polyco)) {
      error_message("depolyco: error reading polyco file");
    } else {
      get_nearest_polyco(polyco_filename,mjd,&polyco);
    }

    i=n=drop=add=0;
    nbytes_per_sample=nchans*nbits*nifs/8;
    dummy=(char *) malloc(nbytes_per_sample);
    for (j=0;j<nbytes_per_sample;j++) dummy[j]=0;
    elapsed_time=barycentre_time=0.0;
    while (i<ntim) {
      n++;
      elapsed_time+=tsamp;
      barycentre_time+=tsamp*polyco_period(mjd,polyco);
      if (elapsed_time-barycentre_time>tsamp) {
	add=1;
	elapsed_time-=tsamp;
      } else if (barycentre_time-elapsed_time>tsamp) {
	drop=1;
	elapsed_time+=tsamp;
      }
      if (drop || add) {
	rawdata=(char *) malloc(n*nbytes_per_sample);
	fread(rawdata,1,n*nbytes_per_sample,input);
	if (drop) {
	  ndropped++;
	  fwrite(rawdata,1,(n-1)*nbytes_per_sample,output);
	} else {
	  nadded++;
	  fwrite(rawdata,1,n*nbytes_per_sample,output);
	  fwrite(dummy,1,nbytes_per_sample,output);
	}
	n=0;
	free(rawdata);
	drop=add=0;
        sprintf(message,"time:%.1fs",elapsed_time);
        update_log(message);
      }
      mjd+=tsamp/86400.0;
      get_nearest_polyco(polyco_filename,mjd,&polyco);
      i++;
    }
    if (n) {
      rawdata=(char *) malloc(n*nbytes_per_sample);
      fread(rawdata,1,n*nbytes_per_sample,input);
      fwrite(rawdata,1,n*nbytes_per_sample,output);
      free(rawdata);
    }
    free(dummy);
    update_log("finished");
    close_log("barycentre.monitor");
    if (verbose && nadded)
      fprintf(stderr,"added %d samples\n",nadded);
    if (verbose && ndropped)
      fprintf(stderr,"dropped %d samples\n",ndropped);
  }
}
