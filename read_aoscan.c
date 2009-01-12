#include <stdio.h>
#include <stdlib.h>
/* returns day number, year (yyyy), and scan number from an aoscan number */
void read_aoscan(unsigned long aoscan, int *day, int *year, int *scan) /*includefile*/
{
  char cao[12],cscan[3],cyear[4],cday[3];
  int i,l,s;

  sprintf(cao,"%u",aoscan);
  l=strlen(cao);
  /* scan number */
  sprintf(cscan,"%c%c%c",cao[l-3],cao[l-2],cao[l-1]);  
  *scan=atoi(cscan);
  /* year */
  if (cao[l-5]=='9') {
    sprintf(cyear,"%c%c",cao[l-5],cao[l-4]);
    *year=1900+atoi(cyear);
    s=l-6;
  } else {
    sprintf(cyear,"%c%c%c%c",cao[l-7],cao[l-6],cao[l-5],cao[l-4]);
    *year=atoi(cyear);
    s=l-8;
  }
  /* day */
  for (i=0; i<=s; i++) {
	cday[i]=cao[i];
  }
  *day=atoi(cday);
  if (*year > 2010) *year=0;
}
