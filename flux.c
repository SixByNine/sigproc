/*
  profile - a program to produce quick ASCII/grey-scale plots of profiles
  to the standard output. To use, simply fold a time series file and pipe
  the output to this program. For profiles produced in regular prf format
  an ASCII histogram is produced showing a 64-bin version of the profile.
  For profiles produced with -stream option, a pseudo grey-scale output
  is produced showing the profiles collapsed down to one line per profile.

  Created April 9/10, 2002 - drl@jb.man.ac.uk
 */
#include "profile.h"
void getflux(float *profile, int nbins, float *smean, float *speak) 
{
  int i;
  float sum;
  sum=0.0;
  *speak=profile[0];
  for (i=0;i<nbins;i++) {
    sum+=profile[i];
    if (profile[i] > *speak) *speak=profile[i];
  }
  *smean=sum/(float)nbins;
}
void outflux(float smean, float speak, int mean, int peak, char *unit)
{
  if (strings_equal(unit,"mJy")) {
    smean*=1.0e3;
    speak*=1.0e3;
  } else if (strings_equal(unit,"uJy")) {
    smean*=1.0e6;
    speak*=1.0e6;
  } else if (!strings_equal(unit,"Jy")) {
    error_message("invalid unit specified for flux density: Jy, mJy or uJy");
  }
  if (mean) printf("Smean = %.3f %s ",smean,unit);
  if (peak) printf("Speak = %.3f %s ",speak,unit);
  if (mean || peak) puts("");
}
main(int argc, char **argv) 
{
  float tsec,fmhz,smean,speak;
  int i,dummy,peak=0,mean=0;
  char line[240], hash, message[80], unit[3];

  input=stdin;
  strcpy(unit,"mJy");
  if (argc > 1) {
    print_version(argv[0],argv[1]);
    if (help_required(argv[1])) {
      /*flux_help();*/
      exit(0);
    }
    i=1;
    while (i<argc) {
      if (file_exists(argv[i])) {
	input=open_file(argv[i],"r");
      } else if (strings_equal(argv[i],"-peak")) {
	peak=1;
      } else if (strings_equal(argv[i],"-mean")) {
	mean=1;
      } else if (strings_equal(argv[i],"-unit")) {
	strcpy(unit,argv[++i]);
      } else {
	sprintf(message,"command-line argument %s not recognized...",argv[i]);
	error_message(message);
      }
      i++;
    }
  }
  if ((peak==0) && (mean==0)) mean=1;
  fgets(line,sizeof(line),input);
  sscanf(line,"%c %lf %lf %lf %ld %lf %lf %d",&hash,
	 &mjdobs,&tstart,&period,&np,&fch1,&refdm,&nbins);
  if (nbins > 0) {
    while (!feof(input)) {
      profile=(float *) malloc(sizeof(float)*nbins);
      for (i=0; i<nbins; i++) fscanf(input,"%d %f",&dummy,&profile[i]);
      getflux(profile,nbins,&smean,&speak);
      outflux(smean,speak,mean,peak,unit);
      free(profile);
      nbins=0;
      fgets(line,sizeof(line),input);
      fgets(line,sizeof(line),input);
      sscanf(line,"%c %lf %lf %lf %ld %lf %lf %d",&hash,
	 &mjdobs,&tstart,&period,&np,&fch1,&refdm,&nbins);
      if (nbins<=0) break;
    }
  } else {
    while (!feof(input)) {
      sscanf(line,"#START %d %f %f",&nbins,&tsec,&fmhz);
      if (nbins <=0) break;
      profile=(float *) malloc(sizeof(float)*nbins);
      for (i=0; i<nbins; i++) fscanf(input,"%d %f",&dummy,&profile[i]);
      getflux(profile,nbins,&smean,&speak);
      outflux(smean,speak,mean,peak,unit);
      fgets(line,sizeof(line),input);
      fgets(line,sizeof(line),input);
      free(profile);
      fgets(line,sizeof(line),input);
      nbins=0;
    }
  }
}
