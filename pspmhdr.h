/*
 * timing mode header information.
 *
 * version 2 July 29, 1994
 */

#define	PSPM_TIMING_HEADER_SIZE	(1<<15)	/* DO NOT CHANGE THIS		    */

#define	PSPM_TIMING_HEADER_VERSION	2


typedef struct {

  double psr_period;
  double samp_rate;
  double psr_dm;
  double psmon_az;
  double psmon_za;
  double psmon_ast;
  double freq;                  /* Sky Frequency (MHz)                       */
  double tick_offset;           /* delay from 10sec tick to start time (usec)*/
  double ast_start;             /* AST start time -- seconds past midnight   */
                                /* including the tick_offset measurement.    */
  double length_of_integration; /* estimate that does NOT take into account  */
                                /* the doppler corrections to samp_rate.     */
				/* For true integration time, the polyco.dat */
				/* entry appended on the end of the scan     */
				/* must be consulted.			     */
  long header_version;
  long bit_mode;
  long num_phase_bins;
  unsigned long scan_num;
  long tc;
  long num_chans;
  long num_periods;             /* exact number of pulse periods integrated  */
  long sizeof_polyco;
  long psmon_wrap;
  long psmon_feed;
  long psmon_daynumber;
  long psmon_epoch;
  char psr_name[12];
  char psr_ra[12];
  char psr_dec[12];
  char date[12];
  char start_time[12];
  char even_word_boundary_filler[4];
  double psmon_ra;
  double psmon_dec;
  double psmon_tolerance;
  char filler[32544];
  long UPDATE_DONE;
  long HEADER_TYPE;

} PSPM_TIMING_HEADER;


#ifndef PSPM_HEADER_SIZE
#define PSPM_HEADER_SIZE	(1<<15)
#endif

#define	DRIFT_TYPE	0
#define	POINT_TYPE	2


typedef struct {

    double samp_rate;
    double pasmon_az;
    double pasmon_za;
    double user_az;
    double user_za;
    double pasmon_lmst;		/* local mean siderial time in seconds */
    double rf_freq;		/* Sky Frequency which is down converted to  */
				/* IF frequency PSPM_IF (MHz)		     */
    double tick_offset;
    double bw;
    double length_of_integration; /* some fixed number */
    long header_version;
    long scan_file_number;
    long bit_mode;
    unsigned long scan_num;
    long tc;
    long num_chans;
    long pasmon_wrap;
    long pasmon_feed;
    long pasmon_daynumber;
    long pasmon_ast;		/* hhmmss */
    char psr_name[12];
    char date[12];
    char start_time[12];
    long file_size;
    long tape_num;
    long tape_file_number;
    char obs_group[12];
    char even_word_boundary_filler[4];
    double user_ra;		/* J2000 (10000.*hr+100.*min+sec)	    */
    double user_dec;		/* J2000 (10000.*deg+100.*min+sec)	    */
    double chan_first_freq;	/* IF center frequency of first channel (MHz)*/
    double chan_spacing;	/* Spaceing between adjacent channel center */
				/* frequencies (MHz)			    */
    int	SIDEBAND;		/* sideband				    */
    char filler[32536];
    long BACKEND_TYPE;
    long UPDATE_DONE;
    long HEADER_TYPE;

} PSPM_SEARCH_HEADER;

#define PSPM_BYT_BLOCK 32768
#define PSPM_INT_BLOCK PSPM_BYT_BLOCK/4
#define PSPM_SAM_BLOCK 512
#define PSPM_NCH_BLOCK 128
#define PSPM_REA_BLOCK PSPM_SAM_BLOCK*PSPM_NCH_BLOCK
