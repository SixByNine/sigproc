#include "filterbank.h"
FILE *input, *output;

main(int argc, char **argv)
{
  int i, start_sample, block_size;
  long int offset, samples_to_read, bytes_to_read;
  char *block;
  fpos_t pos;

  if (argc<4) {
    puts("usage: extract filterbank_file start_sample samples_to_read");
    exit(0);
  }
   
  if (file_exists(argv[1])) {
    strcpy(inpfile,argv[1]);
    input=open_file(inpfile,"rb");
    output=stdout;
  } else {
    error_message("input file does not exist...");
  }
  start_sample=atol(argv[2]); samples_to_read=atol(argv[3]);
  
  if ((read_header(input)) && (data_type==1)) {
    obits=-1;
    for (i=0; i<nifs; i++) ifstream[i]='Y';
    block_size = nifs*nchans*nbits/8;
    block = (char *) malloc(block_size);
    start_time=(double)(start_sample-1)*tsamp;
    filterbank_header(output);
    offset=(start_sample-1);
    /* jumping to the start position */
    for (i=0; i<block_size; i++) {
      if (fseek(input,offset,SEEK_CUR) != 0) {
	error_message("error setting start position in file");
      }
    }
    /* copying the samples */
    for (i=0; i<samples_to_read; i++) {
      fread(block,1,block_size,input);
      fwrite(block,1,block_size,output);
    }
  } else {
    error_message("extract can currently only work with filterbank data!");
  }

}
