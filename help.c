/* help.c - on-line synopsis for each program called via "programname help" */
int help_required(char *string) /* includefile */
{
  if (strings_equal(string,"help")) return(1);
  if (strings_equal(string,"HELP")) return(1);
  if (strings_equal(string,"-help")) return(1);
  if (strings_equal(string,"-HELP")) return(1);
  if (strings_equal(string,"--help")) return(1);
  if (strings_equal(string,"--HELP")) return(1);
  if (strings_equal(string,"-h")) return(1);
  if (strings_equal(string,"-H")) return(1);
  if (strings_equal(string,"--h")) return(1);
  if (strings_equal(string,"--H")) return(1);
  return(0);
}
void giant_help() /*includefile*/
{
  puts("");
  puts("giant - interactive search for giant pulses in time series");
  puts("");
  puts("usage: giant {timfile1} {timfile2} .... {timfilen} -{options}");
  puts("");
  puts("timfile1-n  - the name of the time series to be read");
  puts("-s nsamp    - number of samples to skip (def=0)");
  puts("-n ngulp    - size of a gulp to view in samples (def=all)");
  puts("");
}
void blanker_help() /*includefile*/
{
  puts("");
  puts("blanker - blanks out pulses from a time series");
  puts("");
  puts("usage: blanker {timfile} -{options}");
  puts("");
  puts("timfile     - the name of the time series to be read");
  puts("-p polyco   - name of polyco file to use (def=polyco.dat)");
  puts("-s startphi - starting pulse phase to blank (0->1)");
  puts("-f finisphi - final pulse phase to blank (0->1)");
  puts("");
}
void depolyco_help() /*includefile*/
{
  puts("");
  puts("depolyco - resample a time series to either barycentric or pulsarcentric frames");
  puts("");
  puts("usage: depolyco {timfile} {polycofile} -{options}");
  puts("");
  puts("timfile     - the name of the time series to be read");
  puts("polycofile  - the name of the polyco file for the pulsarcentric case");
  puts("-singlebyte - write output as unsigned characters (def=floats)");
  puts("-raj        - use different RA (J2000; hh:mm:ss.s) than header");
  puts("-decj       - use different DEC (J2000; dd:mm:ss.s) than header");
  puts("-verbose    - write out TEMPO information to stderr (def=quiet)");
  puts("");
  puts("N.B. no polyco file implies barycentric correction to be applied");
  puts("");
}
void barycentre_help() /*includefile*/
{
  puts("");
  puts("barycentre - refer a datafile to a frame at rest wrt the solar system barycentre");
  puts("");
  puts("usage: barycentre inputfile -{options} > outputfile");
  puts("");
  puts("inputfile   - the name of the filterbank/time series file");
  puts("-mypolyco   - take user-defined polyco.bar file (def=create one)");
  puts("-verbose    - write out barycentre information to stderr (def=quiet)");
  puts("");
}
void profile_help() /*includefile*/
{
  puts("");
  puts("profile - produce ASCII or pseudo grey-scale displays of folded data");
  puts("");
  puts("usage: profile {filename} -{options}");
  puts("");
  puts("options:");
  puts("");
  puts("   filename - profile file (def=stdin)");
  puts("-p fraction - set max value to fraction of peak (def=1.0)");
  puts("-frequency  - label grey-scale profiles in frequency (def=time)");
  puts("");
}
void bandpass_help() /*includefile*/
{
  puts("");
  puts("bandpass - outputs the pass band from a filterbank file\n");
  puts("usage: bandpass {filename} -{options}\n");
  puts("options:\n");
  puts("   filename - filterbank data file (def=stdin)");
  puts("-d numdumps - number of dumps to average over (def=all)");
  puts("-t dumptime - number of seconds to average over (def=all)");
  puts("");
}
void decimate_help() /*includefile*/
{
  puts("");
  puts("decimate - reduce time and/or frequency resolution of filterbank data\n");
  puts("usage: decimate {filename} -{options}\n");
  puts("options:\n");
  puts("   filename - filterbank data file (def=stdin)");
  puts("-c numchans - number of channels to add (def=all)");
  puts("-t numsamps - number of time samples to add (def=none)");
  puts("-T numsamps - (alternative to -t) specify number of output time samples");
  puts("-n numbits  - specify output number of bits (def=input)");
  puts("-headerless - do not broadcast resulting header (def=broadcast)");
  puts("");
}
void dedisperse_help() /*includefile*/
{
  puts("");
  puts("dedisperse  - form time series from filterbank data or profile from folded data\n");
  puts("usage: dedisperse {filename} -{options}\n");
  puts("options:\n");
  puts("   filename - full name of the raw data file to be read (def=stdin)");
  puts("-d dm2ddisp - set DM value to dedisperse at (def=0.0)");
  puts("-b numbands - set output number of sub-bands (def=1)");
  puts("-B num_bits - set output number of bits (def=32)");
  puts("-o filename - output file name (def=stdout)");
  puts("-c minvalue - clip samples > minvalue*rms (def=noclip)");
  puts("-f reffreq  - dedisperse relative to refrf MHz (def=topofsubband)");
  puts("-F newfreq  - correct header value of centre frequency to newfreq MHz (def=header value)");
  puts("-n num_bins - set number of bins if input is profile (def=input)");
  puts("-i filename - read list of channels to ignore from a file (def=none)");
  puts("-p np1 np2  - add profile numbers np1 thru np2 if multiple WAPP dumps (def=all)");
  puts("-j Jyfactor - multiply dedispersed data by Jyfactor to convert to Jy");
  puts("-J Jyf1 Jyf2 - multiply dedispersed data by Jyf1 and Jyf2 to convert to Jy (use only for two-polarization data)");
  puts("-wappinvert - invert WAPP channel order (when using LSB data) (def=USB)");
  puts("-wappoffset - assume wapp fsky between two middle channels (for pre-52900 data ONLY)");
  puts("-rmean      - subtract the mean of channels from each sample before dedispersion (def=no)");
  puts("-swapout    - perform byte swapping on output data (def=native)");
  puts("-nobaseline - don't subtract baseline from the data (def=subtract)");
  puts("-sumifs     - sum 2 IFs when creating the final profile (def=don't)");
  puts("-headerless - write out data without any header info");
  puts("-epn        - write profiles in EPN format (def=ASCII)");
  puts("-asciipol   - write profiles in ASCII format for polarization package");
  puts("-stream     - write profiles as ASCII streams with START/STOP boundaries");
  puts("");
}
void fake_help() /*includefile*/
{
  puts("");
  puts("fake - produce fake filterbank format data for testing downstream code\n");
  puts("usage: fake -{options}\n");
  puts("options:\n");
  puts("-period   p - period of fake pulsar in ms (def=random)");
  puts("-width    w - pulse width in percent (def=4)");
  puts("-snrpeak  s - signal-to-noise ratio of single pulse (def=1.0)");
  puts("-dm       d - dispersion measure of fake pulsar (def=random)");
  puts("-nbits    b - number of bits per sample (def=4)");
  puts("-nchans   n - number of filterbank channels (def=128)");
  puts("-tsamp    t - sampling time in us (def=80)");
  puts("-tobs     t - observation time in s (def=10)");
  puts("-tstart   t - MJD time stamp of first sample (def=50000.0)");
  puts("-nifs     n - number of IFs (def=1)");
  puts("-fch1     f - frequency of channel 1 in MHz (def=433.968)");
  puts("-foff     f - channel bandwidth in MHz (def=0.062)");
  puts("-seed     s - seed for Numerical Recipes ran1 (def=seconds since midnight)");
  puts("-nosmear    - do not add in dispersion/sampling smearing (def=add)");
  puts("-swapout    - perform byte swapping on output data (def=native)");
  puts("-evenodd    - even channels=1 odd channels=0 (def=noise+signal)");
  puts("-headerless - do not write header info at start of file (def=header)");
  puts("");
  puts("binary options:\n");
  puts("-binary     - create binary system");
  puts("-bper       - orbital period in hours (def=10.0)");  
  puts("-becc       - eccentricity (def=0.0, circular)");
  puts("-binc       - inclination in degrees (def=90.0)");
  puts("-bomega     - longitude of periastron in degrees (def= 0.0)");
  puts("-bphase     - starting orbital phase (number between 0 and 1, def=0.0)");
  puts("-bpmass     - pulsar mass in solar units (def=1.4))");
  puts("-bcmass     - companion mass in solar units (def=5.0)");
  puts("");
}
void filterbank_help() /*includefile*/
{
  puts("");
  puts("filterbank - convert raw pulsar-machine data to filterbank format\n");
  puts("usage: filterbank <rawdatafile1> .... <rawdatafileN> -{options}\n");
  puts("rawdatafile - raw data file (recognized machines: WAPP, PSPM, OOTY)");
  puts("\noptions:\n");
  puts("-o filename - output file containing filterbank data (def=stdout)");
  puts("-s skiptime - skip the first skiptime (s) of data (def=0.0)");
  puts("-r readtime - read readtime (s) of data (def=all)");
  puts("-i IFstream - write IFstream (IFstream=1,2,3,4)");
  puts("-n nbits    - write n-bit numbers (def=input format)");
  puts("-c minvalue - clip DM=0 samples > mean+minvalue*sigma (def=noclip)");
  puts("-swapout    - perform byte swapping on output data (def=native)");
  puts("-floats     - write floating-point numbers (equal to -n 32)");
  puts("-sumifs     - sum IFs 1+2 to form total-power data");
  puts("-headerfile - write header parameters to an ASCII file (head)");
  puts("-headeronly - write ONLY binary header parameters");
  puts("\noptions for correlator (currently WAPP) data:\n");
  puts("-hamming    - apply Hamming window before FFT (def=nowindow)");
  puts("-hanning    - apply Hanning window before FFT (def=nowindow)");
  puts("-novanvleck - don't do van Vleck correction before FFT (def=doit)");
  puts("-invert     - invert the band after FFT (def=noinversion)");
  puts("-zerolag    - write just the zero-lag value for each IF");
  puts("-rawcfs     - write raw correlation functions (novanvleck)");
  puts("-corcfs     - write corrected correlation functions (vanvleck)");
  puts("");
}
void fold_help() /*includefile*/
{
  puts("");
  puts("fold - fold filterbank channels/time series data\n");

  puts("usage: fold {filename} -{options}\n");
  puts("options:\n");
  puts("   filename - full name of the raw data file to be read (def=stdin)");
  puts("-o out_file - output file for pulse profile data (def=stdout)");
  puts("-p fold_prd - period to fold (ms) or polyco file (def=polyco.dat)");
  puts("-a accelern - fold using constant acceleration (def=0 m/s/s)");
  puts("-f p_factor - multiply the period by p_factor (def=1.0)");
  puts("-m m_factor - output multiple profiles (STREAM only; def=1)");
  puts("-n num_bins - number of bins in folded profile(s) (def=window/tsamp)");
  puts("-d time/num - dump profiles every time s or num pulses (def=nodumps)");
  puts("-t samptime - hard-wire the sampling time (us) (def=header)");
  puts("-l phaseval - phase value (turns) of left edge of pulse (def=0.0)");
  puts("-r phaseval - phase value (turns) of right edge of pulse (def=1.0)");
  puts("-j Jyfactor - multiply all profiles by Jyfactor to convert to Jy");
  puts("-b baseline - subtract baseline from all profiles (def=autobase)");
  puts("-dt timeoff - add a time offset in seconds to tstart (def=0.0)");
  puts("-sk skiptim - skip the first skiptim s before folding (def=0.0)");
  puts("-re readtim - read and fold only readtim s of data (def=ALL)");
  puts("-ascii      - write profiles as ASCII numbers (this is the default)");
  puts("-epn        - write profiles in EPN format (def=ASCII)");
  puts("-acc        - write out accumulated profiles (def=subints)");
  puts("-bin        - write profiles in SIGPROC binary format (def=ASCII)");
  puts("-sub subint - shorthand for -nobaseline -stream -d subint.0");
  puts("-psrfits    - write profiles in PSRFITS format (def=ASCII)");
  puts("-totalpower - sum polarizations 1+2 before writing (def=nosumming)");
  puts("-asciipol   - write profiles in JMCs ASCII format for polarization");
  puts("-stream     - write profiles as ASCII streams with START/STOP bounds");
  puts("-nobaseline - don't subtract baseline from profiles (def=subtract)");
  puts("");
}

void tune_help() /*includefile*/
{ 
	puts("");
	puts("tune - fine tune a period in a time series by stacking sub-integration");
	puts("");
	puts("usage: tune {filename} -{options}\n");
	puts("options:\n");
	puts("   filename - full name of the raw data file to be read (def=stdin)");

	puts("-p fold_prd - center period to fold (ms) or polyco file (def=polyco.dat)");
	puts("-a accelern - center acceleration to fold (def=0 m/s/s)");
	puts("-f p_factor - multiply the period by p_factor (def=1.0)");
	puts("-n num_bins - max number of bins in folded profile(s) (def=64)");
	puts("-t samptime - hard-wire the sampling time (us) (def=header)");
/*	I don't think these work anymore MJK

	puts("-dt timeoff - add a time offset in seconds to tstart (def=0.0)");
	puts("-sk skiptim - skip the first skiptim s before folding (def=0.0)");
	puts("-re readtim - read and fold only readtim s of data (def=ALL)");*/
	puts("-sub subint - Number of subints to use (def=128)");
	puts("-quikgray   - plot data using pdm-style quikgray code (def=pggray)");
	puts("-useaccn    - Do an acceleration search *Experimental* (def=don't)");
        puts("-usejerk    - Do an acceleration search *Experimental* (def=don't)");
	puts("-pf pfactor - divide the period search range by this number (def=1)");
        puts("-af afactor - divide the accel search range by this number (def=1)");
        puts("-jf jfactor - divide the jerk search range by this number (def=1)");
	puts("-format fmt - set the output format, standard PGPLOT formats (def=/xserv)");
	puts("-jreaper f  - write out ascii based output to file f (def=don't)");
	puts("-bestfile f - read the 'best' summery file for candidate info (def=don't)");
	puts("");
}

void header_help() /*includefile*/
{
  puts("");
  puts("header  - examine header parameters of filterbank data\n");
  puts("usage: header {filename} -{options}\n");
  puts("filename is the filterbank data file (def=stdin)\n");
  puts("options:\n");
  puts("-telescope   - return telescope name");
  puts("-machine     - return datataking machine name");
  puts("-source_name - return source name");
  puts("-fch1        - return frequency of channel 1 in MHz");
  puts("-foff        - return channel bandwidth in MHz");
  puts("-nchans      - return number of channels");
  puts("-tstart      - return time stamp of first sample (MJD)");
  puts("-tsamp       - return sample time (us)");
  puts("-nbits       - return number of bits per sample");
  puts("-nifs        - return number of IF channels");
  puts("-headersize  - return header size in bytes");
  puts("-datasize    - return data size in bytes if known");
  puts("-nsamples    - return number of samples if known");
  puts("-tobs        - return length of observation if known (s)");
  puts("");
}
void reader_help() /*includefile*/
{
  puts("");
  puts("reader  - look at filterbank data in ASCII format\n");
  puts("usage: reader {filename} -{options}\n");
  puts("filename is the filterbank data file (def=stdin)\n");
  puts("options:\n");
  puts("-c        c - output only frequency channel c (1...nchans) (def=all)");
  puts("-i        i - output only IF channel i (1...nifs) (def=all)");
  puts("-numerate   - precede each dump with sample number (def=time)");
  puts("-noindex    - do not precede each dump with number/time");
  puts("-stream     - produce a stream of numbers with START/STOP boundaries");
  puts("");
}
