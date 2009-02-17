#include "decimate.h"
void decimate_header() /* includefile */
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
  if (nchans/naddc==1) {
    send_int("data_type",2);
    send_double("refdm",0.0);
  } else {
    send_int("data_type",1);
  }
  send_double("fch1",fch1);
  send_double("foff",foff*(double)naddc);
  send_int("nchans",nchans/naddc);
  send_int("nbits",obits);
  send_double ("tstart",tstart); 
  send_double("tsamp",tsamp*(double)naddt);
  send_int("nifs",nifs);
  send_int("barycentric",barycentric);
  send_string("HEADER_END");
}
