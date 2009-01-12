/**********************************************************************************/
/* This is to convert combined 16bit data into sigproc filterbank format..        */
/* i.e. it flips the frequency ordering as well writes a header...                */
/*                  - Sarala 7 Nov 2006                                           */
/*                                                                                */
/* This is modified such that it convert both the 16bit and 8bit combined data    */
/*                  - Sarala 7 Feb 2007                                           */
/*                                                                                */
/* gcc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -o convertSig convertSig.c -lm  */
/* /astro/pulsar/sarala/sig_install/libsigproc_linux.a                            */
/*                                                                                */
/**********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include <time.h>
#include "sigproc.h"
#include "header.h"
#include <errno.h>

extern int errno;
FILE *output;

#define nch 256
#define maxsamp 80000

main(int argc, char** argv)
{
    int 	  c,arg=1,i,j,k,nstep;
    long long int fsizeusb,nsamp,nextract=0,nsampuse;
    char	  fnameusb[180];
    FILE	  *finfileusb,*flog;
    short int     *DATA1,*DATA4,temp;
    unsigned char *DATA2,*DATA3;
    char          timeBuffer[80];
    time_t        now ;


    if (argc <= 1)
    {
      printf("USAGE : convertSig -f filename -r raj -d decj -m MJD -s timeSample -n nbits -w channelWidth -c first_channelFreq -t No of Channels \n");
      printf("\n");
      printf("    filename: Combined Data filename \n");
      printf("\n");
      exit(0);
    }

      flog=fopen("convert.log","a");
      one_more_time:
      now = time(0);
      strcpy(timeBuffer, ctime(&now));
      printf("\n\nStarting on %s\n",timeBuffer);
      fprintf(flog,"----------------------------\n");
      fprintf(flog,"Starting on %s\n",timeBuffer);

     while((c = getopt(argc, argv, "f:r:d:m:s:n:w:c:t:")) != -1)
      {
        switch (c)
        {
	 case 'f':
         sscanf(optarg,"%s",fnameusb);
              if((finfileusb = fopen(fnameusb,"rb"))== NULL)
              {
                    printf("%s\n",fnameusb);
                    printf("Error opening file \n");
                    fprintf(flog,"Error opening file \n");
                  //  exit(0);
               }
               else
               {
	            printf("\n opened file: %s \n",fnameusb);
	            fprintf(flog,"\n opened file: %s \n",fnameusb);
               }
        }
        switch (c)
        {
          case 'r':
             sscanf(optarg,"%lf",&src_raj);
             printf("testing %lf\n",src_raj);
        }
        switch (c)
        {
          case 'd':
             sscanf(optarg,"%lf",&src_dej);
        }
        switch (c)
        {
          case 'm':
             sscanf(optarg,"%lf",&tstart);
        }
        switch (c)
        {
          case 's':
             sscanf(optarg,"%lf",&tsamp);
        }
        switch (c)
        {
          case 'n':
             sscanf(optarg,"%d",&nbits);
        }
        switch (c)
        {
          case 'w':
             sscanf(optarg,"%lf",&foff);
        }
        switch (c)
        {
          case 'c':
             sscanf(optarg,"%lf",&fch1);
        }
        switch (c)
        {
          case 't':
             sscanf(optarg,"%d",&nchans);
        }

      }

     printf("Filling the header for sigproc file\n");
     strcpy(source_name,"DEBUG");
     data_type=1;
     nifs=1;

/*   
     machine_id=1;    // this is the id used for analysing sigproc data at NCRA
     telescope_id=3;  // this Arecibo id is used for analysing sigproc data at NCRA for the time being..
*/

     machine_id=9;    // these are the id's set for GMRT in sigproc, as per Dunc Lorimer 
     telescope_id=7;  // these are the id's set for GMRT in sigproc, as per Dunc Lorimer 

     az_start=0.0;
     za_start=0.0;
     printf("machine_id %d ; telescope_id %d  \n",machine_id,telescope_id);

/******************************************/

     /* Find the file size */
      fseeko(finfileusb,0,SEEK_END);
      fsizeusb = ftello(finfileusb);
      fseeko(finfileusb,0,SEEK_SET);

      printf("File size : %lld\n",fsizeusb);
      fprintf(flog,"File size : %lld\n",fsizeusb);
      
     /* Move the file pointer */
              
      printf("Location of file ptr: %lld\n",ftello(finfileusb));

      output = fopen("convertSig.raw","w");
      printf("Outputfile: convertSig.raw\n");
      fprintf(flog,"Outputfile: convertSig.raw\n");

      nsamp = fsizeusb;

if (nbits == 16) { // for 16bit data, read and write it as short int
    DATA1 = malloc (2*nch*maxsamp * sizeof (short int *));
    DATA4 = malloc (2*nch*maxsamp * sizeof (short int *));
    nsamp = nsamp/1024; /* converting nsamp in units of time samples */
    printf("Number of bytes: %lld, Number of time samples:%lld\n",1024*nsamp,nsamp);
    fprintf(flog,"Number of bytes: %lld, Number of time samples:%lld\n",1024*nsamp,nsamp);
}
else {// for 8bit data, read and write it as unsigned char
   printf("no of bits %d \n",nbits); 
    DATA2 = malloc (2*nch*maxsamp * sizeof (unsigned char *));
    DATA3 = malloc (2*nch*maxsamp * sizeof (unsigned char *));
    nsamp = nsamp/512; /* converting nsamp in units of time samples */
    printf("Number of bytes: %lld, Number of time samples:%lld\n",512*nsamp,nsamp);
    fprintf(flog,"Number of bytes: %lld, Number of time samples:%lld\n",512*nsamp,nsamp);
}


/* write out the header for the *fil files for sigproc that is needed*/

    send_string("HEADER_START");
    send_string("source_name");
    send_string(source_name);
    send_int("machine_id",machine_id);
    send_int("telescope_id",telescope_id);
    send_int("data_type",data_type);
    send_coords(src_raj,src_dej,az_start,za_start);
    send_double("fch1",fch1);
    send_double("foff",foff);
    send_int("nchans",nchans);
    send_int("nbits",nbits);
    send_double("tstart",tstart);
    send_double("tsamp",tsamp);
    send_int("nifs",nifs);
    send_string("HEADER_END");

// End of the Sigproc header...

     nstep=5;
     nsampuse = nsamp/nstep;

     while(nsampuse > maxsamp)
         {
           nstep=nstep+10;
           nsampuse = nsamp/nstep;
           printf("Nsampuse: %lld\n",nsampuse);
          }
          
      printf("samples per block:%lld\n",nsampuse);
      fprintf(flog,"samples per block:%lld\n",nsampuse);
          
      for(i=0;i<nstep;i++)
         {
           if (nbits == 16) fread(DATA4,sizeof(short),2*nch*nsampuse,finfileusb);
           else fread(DATA2,sizeof(unsigned char),2*nch*nsampuse,finfileusb);
           for(k=0;k<nsampuse;k++)
                {
                  for(j=0;j<nchans;j++) {
                    if (nbits == 16) DATA1[k*nchans+j]=DATA4[k*nchans+(nchans-1-j)];
                    else DATA3[k*nchans+j]=DATA2[k*nchans+(nchans-1-j)];
                  }
                }

           if (nbits == 16) fwrite(DATA1,sizeof(short int),2*nch*nsampuse,output);
           else fwrite(DATA3,sizeof(unsigned char),2*nch*nsampuse,output);
           printf("STEP:%d out of %d steps writing %lld time samples\n",i+1,nstep,(i+1)*nsampuse);
           fprintf(flog,"STEP:%d out of %d steps writing %lld time samples\n",i+1,nstep,(i+1)*nsampuse);
         }
         fprintf(flog,"Finished\n\n");


       fclose(output);
       fclose(finfileusb);
       fclose(flog);
}
