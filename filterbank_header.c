#include "filterbank.h"

void filterbank_header(FILE *outptr) /* includefile */
{
  int i,j;
  output=outptr;
  if (obits == -1) obits=nbits;
  /* go no further here if not interested in header parameters */
  if (headerless) return;
  /* broadcast the header parameters to the output stream */
  if (machine_id != 0) {
    send_string("HEADER_START");
    send_string("rawdatafile");
    send_string(inpfile);
    if (!strings_equal(source_name,"")) {
      send_string("source_name");
      send_string(source_name);
    }
    send_int("machine_id",machine_id);
    send_int("telescope_id",telescope_id);
    send_coords(src_raj,src_dej,az_start,za_start);
    if (zerolagdump) {
      /* time series data DM=0.0 */
      send_int("data_type",2);
      refdm=0.0;
      send_double("refdm",refdm);
      send_int("nchans",1);
    } else {
      /* filterbank data */
      send_int("data_type",1);
      send_double("fch1",fch1);
      send_double("foff",foff);
      send_int("nchans",nchans);
    }
    /* beam info */
    send_int("nbeams",nbeams);
    send_int("ibeam",ibeam);
    /* number of bits per sample */
    send_int("nbits",obits);
    /* start time and sample interval */
    send_double("tstart",tstart+(double)start_time/86400.0);
    send_double("tsamp",tsamp);
    if (sumifs) {
      send_int("nifs",1);
    } else {
      j=0;
      for (i=1;i<=nifs;i++) if (ifstream[i-1]=='Y') j++;
      if (j==0) error_message("no valid IF streams selected!");
      send_int("nifs",j);
    }
    send_string("HEADER_END");
  }
  
}
