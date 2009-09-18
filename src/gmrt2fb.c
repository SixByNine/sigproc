/* 
   gmrt2fb - converts GMRT search-mode data into "filterbank" data 
   adapted from f77 code gmrt2fb. This version produces identical
   16-bit data to the f77 version but also has a new 8-bit option
   which subtracts a mean from the data and writes the result as an
   unsigned character. The resulting output is centred on 128. This
   appears to be just as good as the 16-bit mode. drl-July-20-2005
*/
#include "filterbank.h"
void gmrt2fb(FILE *input, FILE *output) /* includefile*/
{
  double mean,sum,num;
  short junk,result[256];
  unsigned short ur[256];
  int c,r,opened;
  char string[80];
  unsigned char uc[256];
  r=opened=0;
  c=256;
  sum=num=0.0;
  while (!ferror(input)) {
    if ( (fread(&junk,2,1,input)) != 1) return;
    result[c--]= -1*((~junk) & 32767);
    if (c==0) {
      r++;
      if (r>10) {
	for (c=1;c<=256;c++) {
	  ur[c-1]=result[c];
	  uc[c-1]=128+((double)result[c]-mean);
	}
	ur[0]=ur[1]=ur[254]=ur[255]=uc[0]=uc[1]=uc[254]=uc[255]=0;
	if (obits==16) 
	  fwrite(ur,sizeof(short),256,output);
	if (obits==8)
	  fwrite(uc,sizeof(char),256,output);
      } else {
	/* for some reason, the first 10 samples were not
           used in the f77 code, so repeat that here but
           use these samples to calculate a mean */
	for (c=1;c<=256;c++) {
	  sum+=result[c];
	  num+=1.0;
	}
	mean=sum/num;
      }
      if (r%1024 == 0) {
	if (!opened) {
	  open_log("filterbank.monitor");
	  opened=1;
	}
	sprintf(string,"time:%.1fs",r*tsamp);
	update_log(string);
      }
      c=256;
    }
  }
}
