/* pulsar2k2fb - converts PULSAR2000 search-mode data into "filterbank" data */

// 03/Apr/2002	bklein@MPIfR-Bonn.MPG.de

// 04/Jan/2003	add the rebin function	BK
// 04/Sep/2005	correct read-in for 16-bit values


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "filterbank.h"

int psr2Ksubbaseline;
#define BUFSZ		(4*1024)

#define	MIN(x,y)	(((x)<(y))?(x):(y))
#define	MAX(x,y)	(((x)>(y))?(x):(y))


inline  char pack( char i, char j );             // pack a byte in two 4 bits
inline  void unpack( char p, char *i, char *j ); // unpack 4 bits in one byte


void pulsar2k2fb(FILE *input, FILE *output) /* includefile */
{
  FILE			*fpou;
  int			fd[64];
  int			filesize = 0;
  unsigned char		*byteblock, chardata[64];
  unsigned short	*wordblock, shortdata[64];
  float			floatdata[64*2], min = 0.0, max = 255.0;
  float			if1, if2;
  float 		realtime = 0.0;
  char 			string[80];
  unsigned int 		doit, ndump = 0, nread, opened = 0;
  unsigned int		nsamp, i, s, c, smax = (sumifs ? nchans : nchans*nifs);


  /* write output ASCII header file */
  if (headerfile)  {
    /* write output ASCII header file */
    fpou=open_file("head","w");
    fprintf(fpou,"Original PULSAR2000 file: %s\n",inpfile);
    fprintf(fpou,"Sample time (us) : %f\n",tsamp*1.0e6);
    // fprintf(fpou,"Number of samples: %d\n",ns);
    fprintf(fpou,"Number of channels: %d\n",nchans);
    fclose(fpou);
  }

  /* get filesize and calculate the number of samples */
  filesize = sizeof_file("c01.data");
  switch (nbits)  {
    case  4: nsamp = filesize * 2;  break;
    case  8: nsamp = filesize;	    break;
    case 16: nsamp = filesize / 2;  break;
    default: break;
  }

  /* open 'nchans' channel files */
  for (c = 0; c < nchans*nifs; c++)  {
    // if we have FB data we must read the channels in reverse order
    // because channel 1 has the lowest frequency!!!
#if 0	// PSE and the filterbanks have the same channel order!!!  (04.10.2002  BK)
    if (foff < (-4.0))		// we have filterbank data
      sprintf(string, "c%02u.data", nchans*nifs-c);
    else			// we have PSE++ data
      sprintf(string, "c%02u.data", c+1);
#else
    sprintf(string, "c%02u.data", c+1);          
#endif
    fd[c] = open(string, O_RDONLY);
    if ( fd[c] < 0)  {
      fprintf(stderr, "\nERROR: can't access %s! [%s, %u]",
	string, __FILE__, __LINE__);
      exit(-1);
    }
  }

  // allocate dyn. memory for BUFSZ * nifs * nchans in bytes
  if (nbits > 8)
    wordblock = (unsigned short *)malloc( BUFSZ * nifs * nchans * sizeof(short) );
  else
    byteblock = (unsigned char  *)malloc( BUFSZ * nifs * nchans * sizeof(char) );
  if ( ((nbits >  8) && (wordblock == NULL))  ||
       ((nbits < 16) && (byteblock == NULL)) )  {
    fprintf(stderr, "\nERROR: can't allocate enough memory! [%s, %u]\n",
	__FILE__, __LINE__);
    fflush(stderr);
    exit(-1);
  }

  ndump = 0;
  do  {

    switch (nbits)  {
      case 16:		// not well tested !!!
	nread = (ndump < nsamp-BUFSZ) ? BUFSZ : nsamp-ndump;
	for (c = 0; c < nchans*nifs; c++)  {
	  if ( read( fd[c], &wordblock[c*BUFSZ], nread*2 ) != nread*2 )  {
	    fprintf(stderr,"\nERROR: reading PULSAR2000 data! [%s, %u]\n",
		__FILE__, __LINE__);
	    fflush(stderr);
	    exit(-1);
	  }
	}	// end of: channel-loop
	break;
      case  8:
	nread = (ndump < nsamp-BUFSZ) ? BUFSZ : nsamp-ndump;
	for (c = 0; c < nchans*nifs; c++)  {
	  if ( read( fd[c], &byteblock[c*BUFSZ], nread*1 ) != nread )  {
//	  if ( read( fd[c], &byteblock[c*BUFSZ], nread*1 ) == -1 )  {
	    fprintf(stderr,"\nERROR: reading PULSAR2000 data! [%s, %u]\n",
		__FILE__, __LINE__);
	    fprintf(stderr, "DEBUG: nchans: %u  nifs: %u  c: %u  nread: %u\n",
		nchans, nifs, c, nread);
	    fflush(stderr);
	    exit(-1);
	  }
	}	// end of: channel-loop
	break;
      case  4:		// not well tested !!!
	nread = (ndump < nsamp-BUFSZ) ? BUFSZ : nsamp-ndump;
	for (c = 0; c < nchans*nifs; c++)  {
	  if ( read( fd[c], &byteblock[c*BUFSZ], nread/2 ) != nread )  {
	    fprintf(stderr,"\nERROR: reading PULSAR2000 data! [%s, %u]\n",
		__FILE__, __LINE__);
	    fflush(stderr);
	    exit(-1);
	  }
	}	// end of: channel-loop
	break;
    }		// end of switch (nbits)

    /* process loop */
    for (i = 0; i < nread; i++) {
      /* decide whether to process it */
      realtime = tsamp * ndump;
      if ( (doit=process(realtime,start_time,final_time)) == -1)  break;

      // fix a problem with the first sample (also in PULSAR2000 ???)
      if (ndump == 0)  doit = 0;	// simply ignore it!

      if (doit)  {
	/* about to process this record, update log: MONITOR */
	if (ndump % 1024 == 0)  {
	  if (!opened)  {
	    /* open up logfile */
	    open_log("filterbank.monitor");
	    opened = 1;
	  }
	  sprintf(string,"time:%.1fs", realtime);
	  update_log(string);
	}

	/* convert data to float and sum IFs if 'sumifs' is set */
	// psr2Ksubbaseline > 0: baseline substracted !!!
	switch (nbits)  {
	  case 16:	// not well tested !!!
	    if (psr2Ksubbaseline)  {
	      if (sumifs && (nifs > 1))	// add the two IFs and scale it down
	        for (c = 0; c < nchans; c++)  {
		  if1 = (short)wordblock[i + c*BUFSZ];
		  if2 = (short)wordblock[i + (c+nchans)*BUFSZ];
		  floatdata[c] = (if1 + if2) / 2.0; 		// +32768.0;
		}
	      else		// only convert to float
	        for (c = 0; c < nchans*nifs; c++)
		  floatdata[c] = (short)wordblock[i + c*BUFSZ];	// +32768.0;
	    }
	    else  {	// samples with baseline
	      if (sumifs && (nifs > 1))	// add the two IFs and scale it down
	        for (c = 0; c < nchans; c++)  {
		  if1 = (unsigned short)wordblock[i + c*BUFSZ];
		  if2 = (unsigned short)wordblock[i + (c+nchans)*BUFSZ];
		  floatdata[c] = (if1 + if2) / 2.0;
		}
	      else		// only convert to float
	        for (c = 0; c < nchans*nifs; c++)
		  floatdata[c] = (unsigned short)wordblock[i + c*BUFSZ];
	    }
	  break;

	  case 8:
	    if (psr2Ksubbaseline)  {
	      if (sumifs && (nifs > 1))	// add the two IFs and scale it down
	        for (c = 0; c < nchans; c++)  {
		  if1 = (char)byteblock[i + c*BUFSZ];
		  if2 = (char)byteblock[i + (c+nchans)*BUFSZ];
		  floatdata[c] = (if1 + if2) / 2.0;		// + 128.0
		}
	      else		// only convert to float
	        for (c = 0; c < nchans*nifs; c++)
		  floatdata[c] = (char)byteblock[i + c*BUFSZ];	// + 128.0;
	    }
	    else  {	// samples with baseline
	      if (sumifs && (nifs > 1))	// add the two IFs and scale it down
	        for (c = 0; c < nchans; c++)  {
		  if1 = (unsigned char)byteblock[i + c*BUFSZ];
		  if2 = (unsigned char)byteblock[i + (c+nchans)*BUFSZ];
		  floatdata[c] = (if1 + if2) / 2.0;
		}
	      else		// only convert to float
	        for (c = 0; c < nchans*nifs; c++)
		  floatdata[c] = (unsigned char)byteblock[i + c*BUFSZ];
	    }
	  break;

	  case 4:
		// to do !!!!!
	  break;
	}	// end of: switch(nbits)

	//* decide on how to write out data
	switch (obits)  {
	case 32:	// 32-bit float-data
	  // user has requested floating point numbers in binary file
	  if (swapout) for(s=0;s<smax;s++) swap_float(&floatdata[s]);
	  fwrite( floatdata, sizeof(float), smax, output );
	  break;
	case 16:
	  // user has requested unsigned shorts in binary file
	  float2short( floatdata, smax, min, max, shortdata );
	  if (swapout) for(s=0;s<smax;s++) swap_short(&shortdata[s]);
	  fwrite( shortdata, sizeof(unsigned short), smax, output );
	  break;
	case 8:
	  // user has requested unsigned chars in binary file
	  float2char( floatdata, smax, min, max, chardata );
	  fwrite( chardata, sizeof(unsigned char), smax, output );
	  break;
	case 4:
	  // user has requested to write data out packed into character format
	  float2four( floatdata, smax, min, max, chardata );
	  fwrite( chardata, sizeof(unsigned char), smax/2, output );
	  break;
	default:
	  error_message("unknown bit format for writing");
	break;
        }	// end of switch
      }		// end of: if (doit)
      ndump++;
    }		// end of: process - loop

  } while (ndump < nsamp);

  if (nbits > 8)  free(wordblock);
  else		  free(byteblock);

  // close all channel files
  for (c = 0; c < nchans; c++)  close(fd[c]);
}


inline char pack( char i, char j )
{
  return ( (j << 4) + (i & 0x0F) );
}


inline void unpack( char p, char *i, char *j )
{
  *i = p & 0x0F;
  *j = (p & 0xF0) >> 4;
}

// ===================================================================================

