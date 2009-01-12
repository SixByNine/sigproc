#include "dedisperse.h"
main(int argc, char **argv)
{
  FILE *tim;
  char c;
  float s;
  output=stdout;
  send_string("HEADER_START");
  send_int("telescope_id",5);
  send_int("machine_id",6);
  send_int("data_type",2);
  send_double("refdm",0.0);
  send_int("nchans",1);
  send_int("nbits",32);  
  send_double("tsamp",atof(argv[1]));
  send_int("nifs",1);
  send_string("HEADER_END");
  tim=fopen(argv[2],"rb");
  while (!feof(tim)) {
    fread(&c, 1, 1, tim);	
    s=c;
    fwrite(&s,4, 1, output);
  }
}
