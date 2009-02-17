/* 
   depolyco.c - barycentre or take out binary motion 
   contained within a preprepared set of polynomial coefficients 
   in "polyco.dat". The polyco file is made by setting the pulsar
   period to be the data sampling time, tsamp. Code adapted from
   seek's resample.f (dlorimer@atnf.csiro.au; Oct 15, 2003)

   Revised version Nov 1, 2003 (drl@jb.man.ac.uk) runs TEMPO to
   get the barycentric MJD which is written out into the resulting
   resampled time series header. 

   Revised version Aug 20, 2004 (drl@jb.man.ac.yk) runs TEMPO to
   generate a polycofile for the barycentric case. Pulsarcentering
   is done by the user generating a .par file with F0=1/tsamp!

   This _could_ be extended further to run TEMPO in the pulsarcentric
   case also. However more pressing is to really develop code so that
   entire filterbank files can be barycentered/pulsarcentered.
*/
#include "dedisperse.h"
char polyco_filename[80];
int verbose;
/* subroutine to call TEMPO to calculate a polyco.bar file for barycentering */
void make_polycofile(char ra[],char dec[],char topo[], char site,
			    double mjdtopo, double tsamp)
{
  FILE *resid2,*parfile, *tzfile;
  float junk;
  double mjdbary;
  parfile=fopen("tssb.par","w");
  fprintf(parfile,"PSR 0000+00\n");
  fprintf(parfile,"RAJ %s\n",ra);
  fprintf(parfile,"DECJ %s\n",dec);
  fprintf(parfile,"F0 %f\n",1.0/tsamp);
  fprintf(parfile,"DM 0.0\n");
  fprintf(parfile,"PEPOCH %s\n",topo);
  fclose(parfile);
  tzfile=fopen("tz.in","w");
  fprintf(tzfile,"%c    2  30  9 1410\n",site);
  fprintf(tzfile,"\n \n");
  fprintf(tzfile,"0000+00 60 9 12 1410\n");
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
  fprintf(timfile,"%c    0  0000+00 9999.000 %.13f     1.00\n",site,mjdtopo);
  fprintf(timfile,"%c    0  0000+00 9999.000 %.13f     1.00\n",site,mjdtopo);
  fprintf(timfile,"%c    0  0000+00 9999.000 %.13f     1.00\n",site,mjdtopo);
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


int poly_override;
double override_f0;
/* the subroutine that does all the work... adapted from resample.f */
void depolyco(int ntim, double tstart, double tsamp, int singlebyte)
{
  double time,taut,next;
  static int n8k=8192;
  float series[n8k],resamp[n8k],polyfreq,inputfreq;
  unsigned char packed[n8k];
  char message[80];
  int i,j,n;
  struct POLYCO polyco;
  FILE *polycofile;
  
  polycofile=open_file(polyco_filename,"r");
  if (verbose) fprintf(stderr,"opened %s\n",polyco_filename);
  poly_override=1;
  override_f0=1.0/tsamp;
  if (!read_polycoset(polycofile,&polyco)) {
    error_message("depolyco: error reading polyco.dat...");
  } else {
    get_nearest_polyco(polyco_filename,tstart,&polyco);
    taut=polyco_period(tstart,polyco);
    fread(series,sizeof(float),n8k,input);
  }

  polyfreq=(float) polyco.f0;
  inputfreq=(float) 1.0/tsamp;
  if (verbose) fprintf(stderr,"polyfreq %f Hz\n",polyfreq);
  if (polyfreq!=inputfreq)
  error_message("requested polyco file needs to be remade with F0=1.0/tsamp!");

  /* initialize loop counters and buffers */
  for (i=0; i<n8k; i++) resamp[i]=0.0;
  i=j=n=1;
  time=0.0;
  next=taut;

  while (n<ntim) {

    if (next>tsamp*(double)n) {
      resamp[j-1]+=series[i-1]*((double)n*tsamp-(next-taut))/tsamp;
      n++;
      i++;
    }
    if (n==ntim) break;
    
    if (next<=tsamp*(double)n) {
      resamp[j-1]+=series[i-1]*(next-(double)(n-1)*tsamp)/tsamp;
      j++;
      if (j==n8k) {
	/* it's time to write out a block of resampled data */
	if (singlebyte) {
	  /* single byte output assumes that the time series has zero mean! */
	  for (j=0;j<n8k;j++) 
	    packed[j] = 128 + (unsigned char) (resamp[j]);
	  fwrite(packed,sizeof(char),n8k,output);
	} else {
	  fwrite(resamp,sizeof(float),n8k,output);
	}
	/* now clear the block ready for the next set of samples */
	for (j=0;j<n8k;j++) resamp[j]=0.0;
	j=1;
	/* check if we're still reading the best polyco set */
	get_nearest_polyco(polyco_filename,tstart+time/86400.0,&polyco);
      }
      time=next;
      taut=polyco_period(tstart+time/86400.0,polyco);
      next+=taut;
    }


    if (i==n8k) {
      /* it's time to read in a block of raw data */
      fread(series,sizeof(float),n8k,input);
      i=1;
      /* update the log every 8k samples */
      sprintf(message,"time:%.1fs",time);
      update_log(message);
    }
  }
}

char inpfile[80], outfile[80];
main (int argc, char *argv[]) 
{
  int i,ntim,headersize,rah,ram,ded,dem,headerless,singlebyte,notempo;
  double ras,des,mjdbary;
  char ra[80], dec[80], topo[80], sra[6], sde[6], site;
  float *series;
  if (argc<2 || help_required(argv[1])) {
    depolyco_help();
    exit(0);
  } else {
    print_version(argv[0],argv[1]);
  }

  if (!file_exists(argv[1]))
    error_message("time series does not exist!");

  strcpy(inpfile,argv[1]);
  strcpy(outfile,"stdout");
  output=stdout;
  verbose=headerless=singlebyte=0;
  strcpy(polyco_filename,"unknown");
  strcpy(ra,"");strcpy(dec,"");
  i=2;
  notempo=headerless=singlebyte=verbose=0;
  while (i<argc) {
    if (strings_equal(argv[i],"-notempo")) notempo=1;
    if (strings_equal(argv[i],"-headerless")) headerless=1;
    if (strings_equal(argv[i],"-singlebyte")) singlebyte=1;
    if (strings_equal(argv[i],"-verbose")) verbose=1;
    if (strings_equal(argv[i],"-raj")) strcpy(ra,argv[++i]);
    if (strings_equal(argv[i],"-decj")) strcpy(dec,argv[++i]);
    /* if the user supplies a polyco file -- it's pulsarcentric! */
    if (file_exists(argv[i])) {
      strcpy(polyco_filename,argv[i]);
      pulsarcentric=1;
    }
    i++;
  }

  input=open_file(inpfile,"r");

  if ((headersize=read_header(input))) {
    if (data_type != 2) error_message("input file is not a time series!");
    if (nbits != 32) error_message("depolyco can only read 32-bit data!");

    /* parse the header parameters for RA */
    angle_split(src_raj,&rah,&ram,&ras);
    if (ras<10.0) {
      sprintf(sra,"0%.1f",ras);
    } else {
      sprintf(sra,"%.1f",ras);
    }
    if (strings_equal(ra,"")) sprintf(ra,"%02d:%02d:%s",rah,ram,sra);

    /* parse the header parameters for DEC */
    angle_split(src_dej,&ded,&dem,&des);
    if (des<10.0) {
      sprintf(sde,"0%.1f",des);
    } else {
      sprintf(sde,"%.1f",des);
    }
    if (strings_equal(dec,"")) sprintf(dec,"%02d:%02d:%s",ded,dem,sde);

    if (!notempo) {
      /* now call TEMPO to calculate the barycentric MJD */
      sprintf(topo,"%.12f",tstart);
      site=tempo_site(telescope_id);
      if (strings_equal(polyco_filename,"unknown")) {
	make_polycofile(ra,dec,topo,site,tstart,tsamp);
	strcpy(polyco_filename,"polyco.bar");
	if (verbose) fprintf(stderr,"generated polyco.dat\n");
	barycentric=1;
      }
      if (verbose) fprintf(stderr,"RAJ %s DECJ %s SITE %c\n",ra,dec,site);
      if (verbose) fprintf(stderr,"Topocentric MJD %.13f\n",tstart);
      mjdbary=barycentric_time(ra,dec,topo,site,tstart);
      if (verbose) fprintf(stderr,"Barycentric MJD %.13f\n",mjdbary);
    } else {
      mjdbary=tstart; /* no barycentric time correction */
    }

    if (strings_equal(polyco_filename,"unknown")) 
      error_message("polyco file not created...");

    /* write out header with barycentric MJD if required */
    if (!headerless) {
      send_string("HEADER_START");
      if (!strings_equal(source_name,"")) {
	send_string("source_name");
	send_string(source_name);
      }
      send_int("telescope_id",telescope_id); 
      send_int("machine_id",machine_id);
      send_coords(src_raj,src_dej,az_start,za_start);
      send_int("data_type",2);
      send_int("barycentric",barycentric);
      send_int("pulsarcentric",pulsarcentric);
      send_double("refdm",refdm);
      if (fch1 == 0.0) 
	send_double("fch1",frequency_table[0]);
      else
	send_double("fch1",fch1);
      send_int("nchans",1);
      if (singlebyte) 
	send_int("nbits",8);  
      else
	send_int("nbits",32);  
      send_double ("tstart",mjdbary); 
      send_double("tsamp",tsamp);
      send_int("nifs",nifs);
      send_string("HEADER_END");
    }

    /* now call the routine that does the work */
    open_log("depolyco.monitor");
    ntim=nsamples(inpfile,headersize,nbits,nifs,nchans);
    depolyco(ntim,tstart,tsamp,singlebyte);
    update_log("finished");
    close_log("depolyco.monitor");
  }
}
