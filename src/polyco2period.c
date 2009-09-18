/*
  polyco2period-calculates period for MJD specified on command line 
*/
#include "fold.h"
#include <math.h>

void polyco2period_help()
{
  puts("");
  puts("polyco2period: returns period (sec) given mjd and polyco file");
  puts("        usage: polyco2period mjd [-p polyco file name]\n");
  puts("");
}

main (int argc, char *argv[])
{
  /* local variables */
  int i,opened_input=0,opened_output=0;
  char string[80];

  double psec, phase, offset;
  double polyco_period( double mjd, struct POLYCO polyco);
  struct POLYCO polyco;

  if (argc<2) {
    polyco2period_help();
    exit(0);
  }
	offset=0.0;

  strcpy(polyco_file,"polyco.dat");

  /* check the command line parameters */
  i=2;
  while (i<argc) {
    if (strings_equal(argv[i],"-p")) {
      i++;
      if (file_exists(argv[i])) {
	strcpy(polyco_file,argv[i]);
	folding_period=-1.0;
      }
    } else if (strings_equal(argv[i],"-o")) {
      i++;
      offset=atof(argv[i]);
    } else  {
	/* unknown argument passed down - stop! */
	polyco2period_help();
	sprintf(string,"unknown argument (%s) passed to polyco2period",argv[i]);
	error_message(string);
    }
    i++;
  }

  tstart=atof(argv[1])+offset/86400.0;
  get_nearest_polyco(polyco_file,tstart,&polyco); 
  psec=polyco_period(tstart,polyco);
  phase=fabs(polyco_phase(tstart,polyco));
  phase=phase-floor(phase);
  printf("%12.9f (sec) for mjd = %18.12f Phase: %.12f\n", psec,tstart,phase);
  exit(0);
}
