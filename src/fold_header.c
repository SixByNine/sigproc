#include "fold.h"

void fold_header() /* includefile */
{
  send_string("HEADER_START");
  send_int("data_type",3);
  if (!strings_equal(source_name,"")) {
    send_string("source_name");
    send_string(source_name);
  }
  send_int("machine_id",machine_id);
  send_int("telescope_id",telescope_id);
  send_coords(src_raj,src_dej,az_start,za_start);
  send_int("nbits",32);
  send_int("nifs",nifs);
  send_double("fch1",fch1);
  send_double("foff",foff);
  send_double("tstart",tstart);
  if (npuls>0) send_long("npuls",npuls);
  send_int("nchans",nchans);
  send_int("nbins",nbins);
  send_double("period",folding_period/1000.0);
  send_string("HEADER_END");
}
