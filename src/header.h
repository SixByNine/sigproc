/* global variables describing the data */
char rawdatafile[80], source_name[80];
int machine_id, telescope_id, data_type, nchans, nbits, nifs, scan_number,
  barycentric,pulsarcentric; /* these two added Aug 20, 2004 DRL */
double tstart,mjdobs,tsamp,fch1,foff,refdm,az_start,za_start,src_raj,src_dej;
double gal_l,gal_b,header_tobs,raw_fch1,raw_foff;
int nbeams, ibeam;
/* added 20 December 2000    JMC */
double srcl,srcb;
double ast0, lst0;
long wapp_scan_number;
char project[8];
char culprits[24];
double analog_power[2];

/* added frequency table for use with non-contiguous data */
double frequency_table[4096]; /* note limited number of channels */
long int npuls; /* added for binary pulse profile format */
