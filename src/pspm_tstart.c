/*
## help on pspm_tstart 
## double pspm_tstart - returns UT time stamp for the start of any PSPM scan
##
## variables passed down:
##
## long scan_num      - the scan number from the header
## char *start_time   - the time stamp string from the header (HH:MM:SS)
## double tick_offset - the tick offset from the header in us
##
## variables returned:
## 
## double mjdobs    - the MJD of the observations
##
## pspm_tstart returns the number of seconds since midnight for mjdobs
##
## this routine is designed to cope with the change in start time convention
## (AST -> UT) written in the headers before and after Jan 1, 2000. Pre y2k,
## all start times were in AST. Post y2k, all start times are in UT.
##
## Created on Mar 10, 2001 (dunc@naic.edu)
## help end
*/
#include <stdio.h>
#include "sigproc.h"
double pspm_tstart(unsigned long scan_num, char *start_time, double tick_offset, double *mjdobs)/*includefile*/
{
  int dayno,year,scan,hh,mm,ss,seconds;
  /* read in the aoscan to get dayno year and scan */
  read_aoscan(scan_num, &dayno, &year, &scan);

  /* get number of seconds since midnight */
  sscanf(start_time,"%d:%d:%d",&hh,&mm,&ss);
  seconds=hh*3600+mm*60+ss;

  /* 
     before y2k the start time was AST - add on 4 hours (14400 seconds) 
     to convert to UT, checking for crossing midnight and new year during 
     conversion and adjusting year, seconds and day number if necessary 
  */
  if (year<2000) {
    seconds+=14400;
    if (seconds>86400) {
      seconds-=86400;
      dayno++;
    }
    if (dayno>366) {
      dayno=1;
      year++;
    }
  }

  /* form mjd based on jan 1 of current year plus number of days since then */
  *mjdobs=mjd(year,1,1)+(double)(dayno-1);

  /* return the start time in seconds */
  return ((double)seconds + tick_offset*1.0e-6);
}
