#include "dedisperse.h"

void dedisperse_header() /* includefile */
{
  /* go no further here if not interested in header parameters */
  if (headerless) return;

  /* broadcast the header parameters to the output stream */
  send_string("HEADER_START");

  if (!strings_equal(source_name,"")) {
    send_string("source_name");
    send_string(source_name);
  }
  send_int("telescope_id",telescope_id); 
  send_int("machine_id",machine_id);
  send_coords(src_raj,src_dej,az_start,za_start);
  if (refdm == -1.0) refdm=userdm;
  if (nbands==1) {
    send_int("data_type",2);
    send_double("refdm",refdm);
  } else {
    send_int("data_type",6);
    send_double("refdm",refdm);
    send_double("foff",foff*(double)(nchans/nbands));
  }

  if ( (fch1 == 0.0) && (frequency_table[0] == 0.0) ) 
    error_message("help... missing frequency information!");

  if (fch1 == 0.0) 
    send_double("fch1",frequency_table[0]);
  else
    send_double("fch1",fch1);

  send_int("barycentric",barycentric);
  send_int("nchans",nbands);
  send_int("nbits",nobits);  
  send_double ("tstart",tstart); 
  send_double("tsamp",tsamp);
  send_int("nifs",nifs);
  send_string("HEADER_END");
}
