/* blanker - a program to blank out pulses from a time series */

#include "dedisperse.h"
char polyco_filename[80];
char inpfile[80], outfile[80];
main (int argc, char *argv[]) 
{
  int i,ntim,headersize,noff=0,gulp;
  float *time_series,sum=0.0,sumsq=0.0,mean,meansq,sigma;
  double phase_start=0.0,phase_finish=0.0,phase,turn=0.0,last=0.0,period,mjd;
  double period_constant=0.0; /* !=0.0 means no polyco */
  long seed=0;
  struct POLYCO polyco;
  FILE *polycofile;
  if (argc<2 || help_required(argv[1])) {
    blanker_help();
    exit(0);
  }
  print_version(argv[0],argv[1]);
  if (!file_exists(argv[1]))
    error_message("input file does not exist!");

  strcpy(inpfile,argv[1]);
  input=open_file(inpfile,"r");
  strcpy(outfile,"stdout");
  output=stdout;
  strcpy(polyco_filename,"polyco.dat");

  i=2;
  while (i<argc) {
    if (strings_equal(argv[i],"-p"))       strcpy(polyco_filename,argv[++i]);
    if (strings_equal(argv[i],"-s"))       phase_start=atof(argv[++i]);
    if (strings_equal(argv[i],"-f"))       phase_finish=atof(argv[++i]);
    if (strings_equal(argv[i],"-P"))       period_constant=atof(argv[++i]);
    i++;
  }
  if (phase_start==phase_finish)
    error_message("start and end phases must be different");
  if (phase_start>phase_finish)
    error_message("start phases must be less than end phase");

  if ((headersize=read_header(input))) {
    if (nbits!=32) 
      error_message("blanker currently only works for 32-bit data");

    send_string("HEADER_START");

    if (!strings_equal(source_name,"")) {
      send_string("source_name");
      send_string(source_name);
    }
    send_int("telescope_id",telescope_id); 
    send_int("machine_id",machine_id);
    send_coords(src_raj,src_dej,az_start,za_start);
    send_int("data_type",2);
    send_double("fch1",fch1);
    send_double("refdm",refdm);
    send_int("barycentric",barycentric);
    send_int("nchans",1);
    send_int("nbits",nbits);  
    send_double ("tstart",tstart); 
    send_double("tsamp",tsamp);
    send_int("nifs",nifs);
    send_string("HEADER_END");

    gulp=8192;
    ntim=nsamples(inpfile,headersize,nbits,nifs,nchans);
    if (gulp>ntim) gulp=ntim;

    if (period_constant!=0.0) {
      fprintf(stderr,"Using constant period\n");
      period = period_constant;
    }
    else{
      polycofile=open_file(polyco_filename,"r");
      if (!read_polycoset(polycofile,&polyco)) {
        error_message("blanker: error reading polyco file");
      } else {
        get_nearest_polyco(polyco_filename,tstart,&polyco);
        period=polyco_period(tstart,polyco);
      }
    }
    time_series=(float *) malloc(gulp*sizeof(float));
    fread(time_series,sizeof(float),gulp,input);

    /* calculate the off-pulse mean and rms from first 1k points */
    i=0;
    mjd=tstart;
    while (noff<1024) {
      turn += tsamp/period;
      if ((turn-last)>=1.0) {
	if (period_constant==0.0){
	get_nearest_polyco(polyco_filename,tstart,&polyco);
	period=polyco_period(tstart,polyco);
	}
	else period = period_constant;
	last=turn;
      }
      phase = turn-floor(turn);
      if ((phase<phase_start)||(phase>phase_finish)) {
	sum+=time_series[i];
	sumsq+=time_series[i]*time_series[i];
	noff++;
      }
      tstart += tsamp/86400.0;
      i++;
    }
    mean=sum/(float)gulp;
    meansq=sumsq/(float)gulp;
    sigma=sqrt(meansq-mean*mean);

    /* now run through the file blanking with gaussian noise */
    turn=0.0;
    last=0.0;
    tstart=mjd;
    i=0;
    while (!feof(input)) {
      turn += tsamp/period;
      if ((turn-last)>=1.0) {
	if (period_constant>0.0) {
	  period=period_constant;
	}
	else {
	  get_nearest_polyco(polyco_filename,tstart,&polyco);
	  period=polyco_period(tstart,polyco);
	}
	last=turn;
      }
      phase = turn-floor(turn);
      if ((phase>=phase_start)&&(phase<=phase_finish)) 
	time_series[i]=gasdev(&seed)*sigma+mean;
      tstart += tsamp/86400.0;
      i++;
      if (i==gulp) {
	fwrite(time_series,sizeof(float),gulp,output);
	gulp=fread(time_series,sizeof(float),gulp,input);
	i=0;
      }
    }

  }
}
