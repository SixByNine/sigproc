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
main(int argc, char **argv) 
{
  float tsec,fmhz,peak=1.0,sum,sumsq,sigma,sig2,mean,meansq,chi2,diff;
  float *pcopy,n;
  int hh,mm,i,j,row,dummy,mbins=64,obins,first=1,idx,nprof=1,display=0,rc2=0,flux=0;
  unsigned long *indx;
  char line[240], hash, gry[10], message[80];

  strcpy(gry,"   ,:o*@$#");
  input=stdin;

  if (argc > 1) {
    print_version(argv[0],argv[1]);
    if (help_required(argv[1])) {
      profile_help();
      exit(0);
    }
    i=1;
    while (i<argc) {
      if (file_exists(argv[i])) {
	input=open_file(argv[i],"r");
      } else if (strings_equal(argv[i],"-frequency")) {
	display=1;
      } else if (strings_equal(argv[i],"-chi2")) {
	rc2=1;
      } else if (strings_equal(argv[i],"-sum")) {
	flux=1;
      } else if (strings_equal(argv[i],"-p")) {
	peak=atof(argv[++i]);
      } else {
	sprintf(message,"command-line argument %s not recognized...",argv[i]);
	error_message(message);
      }
      i++;
    }
  }

  fgets(line,sizeof(line),input);
  sscanf(line,"%c %lf %lf %lf %ld %lf %lf %d",&hash,
	 &mjdobs,&tstart,&period,&np,&fch1,&refdm,&nbins);
  if (nbins > 0) {
    while (!feof(input)) {
      profile=(float *) malloc(sizeof(float)*nbins);
      for (i=0; i<nbins; i++) {
	fscanf(input,"%d %f",&dummy,&profile[i]);
      }
      if (rc2 || flux) {
	pcopy=(float *) malloc(sizeof(float)*nbins);
	for (i=0;i<nbins;i++) pcopy[i]=profile[i];
	indx=(unsigned long *) malloc(sizeof(unsigned long)*nbins);
	indexx((unsigned long)nbins,pcopy,indx);
	if (flux) {
	  n=sum=sumsq=0.0;
	  for (i=0;i<40;i++) {
	    sum+=profile[i];
	    sumsq+=profile[i]*profile[i];
	    n+=1.0;
	  }
	  mean=sum/n;
	  meansq=sumsq/n;
	  sigma=sqrt(meansq-mean*mean)/sqrt(n-1.0);
	  n=sum=sumsq=0.0;
	  for (i=43;i<52;i++) {
	    sum+=profile[i];
	    sumsq+=profile[i]*profile[i];
	    n+=1.0;
	  }
	  printf("%f %f\n",sum/n,sigma);
	  exit(0);
	}
	n=sum=sumsq=0.0;
	for (i=0;i<nbins/2;i++) {
	  j=indx[i];
	  sum+=profile[j];
	  sumsq+=profile[j]*profile[j];
	  n+=1.0;
	}
	mean=sum/n;
	meansq=sumsq/n;
	sigma=sqrt(meansq-mean*mean);
	sig2=sigma*sigma;
	sumsq=0.0;
	for (i=0;i<nbins;i++) {
	  diff=profile[i]-mean;
	  sumsq+=diff*diff;
	}
	//chi2=sumsq/sig2/(float)(nbins-1);
	chi2=sumsq/(float)(nbins-1);
	printf("%f %f\n",period,chi2);
	exit(0);
      }
      if (nbins>mbins) {
	add_channels(profile,nbins,nbins/mbins);
	nbins=mbins;
      }
      pmin=vmin(profile,nbins);
      pmax=vmax(profile,nbins)*peak;
      prng=pmax-pmin;
      for (i=0; i<nbins; i++) profile[i]=19.0*(profile[i]-pmin)/prng+1.0;
      for (row=20; row>=1; row--) {
	for (i=0; i<nbins; i++) {
	  if (profile[i]>=(float) row) 
	    printf("#");
	  else
	    printf(" ");
	}
	printf("\n");
      }
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
      fgets(line,sizeof(line),input);
      fgets(line,sizeof(line),input);
      if (nbins>mbins) {
	obins=nbins;
	add_channels(profile,nbins,nbins/mbins);
	nbins=mbins;
      }
      if (first) {
	pmin=vmin(profile,nbins);
	pmax=vmax(profile,nbins)*peak;
	prng=pmax-pmin;
	/*first=0;*/
      }
      if (display) {
	printf("|%04d|%11.6f|",nprof,fmhz);
      } else {
	hh=(int) tsec/3600.0;
	tsec-=(float) hh*3600.0;
	mm=(int) tsec/60.0;
	tsec-=(float) mm*60.0;
	printf("|%04d|%02d:%02d:%02d|",nprof,hh,mm,(int)tsec);
      }
      for (i=0; i<nbins; i++) {
	idx=(int) (9.0*(profile[i]-pmin)/prng);
	if (idx<0) idx=0;
	if (idx>9) idx=9;
	putchar(gry[idx]);
      }
      puts("|");
      free(profile);
      fgets(line,sizeof(line),input);
      nbins=0;
      nprof++;
    }
  }
}
