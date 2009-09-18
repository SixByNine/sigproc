/*
  snrdm - a program to produce the expected signal-to-noise versus DM curve 
  for a pulsar search given a set of observing parameters.

  Created July 25, 2002 - drl@jb.man.ac.uk
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sigproc.h"
main(int argc, char **argv) 
{
  char message[80];
  int i;
  float fmhz=430.0, bw=7.68, nch=128.0, dm=50.0, pms=150.0,
        wms=5.0, tms=0.08, tsc=0.0, mdm=0.0, nrm=100.0;
  float xdm, wef, snr, tdm, tdd, snrmax=0.0;
  FILE *input;
  if (argc > 1) {
    print_version(argv[0],argv[1]);
    if (help_required(argv[1])) {
      /*snrdm_help();*/
      exit(0);
    }
    i=1;
    while (i<argc) {
      if (file_exists(argv[i])) {
	input=open_file(argv[i],"r");
      } else if (strings_equal(argv[i],"-f")) {
        fmhz=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-b")) {
        bw=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-c")) {
        nch=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-d")) {
        dm=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-p")) {
        pms=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-w")) {
        wms=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-t")) {
        tms=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-s")) {
        tsc=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-m")) {
        mdm=atof(argv[++i]);
      } else if (strings_equal(argv[i],"-n")) {
        nrm=atof(argv[++i]);
      } else {
	sprintf(message,"command-line argument %s not recognized...",argv[i]);
	error_message(message);
      }
      i++;
    }
  }
  if (mdm == 0.0) mdm=2.0*dm;
  tdm=(float) dmdelay(fmhz-bw/nch,fmhz,dm)*1000.0;
  for (xdm=0.0;xdm<mdm;xdm+=0.1) {
    tdd=(float) dmdelay(fmhz-bw/2.0,fmhz+bw/2.0,abs(dm-xdm))*1000.0;
    wef=sqrt(wms*wms+tsc*tsc+tms*tms+tdm*tdm+tdd*tdd);
    snr=sqrt((pms-wef)/wef);
    snrmax=(snr>snrmax)?snr:snrmax;
  }
  printf("#Signal-to-noise/DM curve assuming following parameters:\n");
  fprintf(stderr,"#Pulse Period %.3f ms\n",pms);
  fprintf(stderr,"#Intrinsic Width %.3f ms\n",wms);
  fprintf(stderr,"#Scattering Time %.3f ms\n",tsc);
  fprintf(stderr,"#True DM %.3f pc/cc\n",dm);
  fprintf(stderr,"#Centre Frequency %.3f MHz\n",fmhz);
  fprintf(stderr,"#Bandwidth %.3f MHz\n",bw);
  fprintf(stderr,"#Number of Channels %d\n",nch);
  printf("#START\n");
  for (xdm=0.0;xdm<mdm;xdm+=0.1) {
    tdd=(float) dmdelay(fmhz-bw/2.0,fmhz+bw/2.0,abs(dm-xdm))*1000.0;
    wef=sqrt(wms*wms+tsc*tsc+tms*tms+tdm*tdm+tdd*tdd);
    if (wef>pms) 
      snr=0.0;
    else
      snr=nrm*sqrt((pms-wef)/wef)/snrmax;
    printf("%.2f %f\n",xdm,snr);
  }
  printf("#STOP\n");
  printf("#DONE\n");
}
