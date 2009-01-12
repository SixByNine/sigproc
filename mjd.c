/* some wrappers for slalib routine calls */
#include <stdio.h>
#include "slalib.h"
#include "slamac.h"
#include "sigproc.h"
double mjd(int year, int month, int day) /*includefile*/
{
  double djm;
  int j;
  char message[80];
  slaCaldj(year,month,day,&djm,&j);
  if (j==0) return (djm);
  sprintf(message,"slaCaldj could not process %d %d %d",year,month,day);
  error_message(message);
}
void cal(double djm, int *year, int *month, int *day) /*includefile*/
{
  int iymdf[4],j;
  char message[80];
  slaDjcal(1,djm,iymdf,&j);
  if (j==0) {
    *year=iymdf[0];
    *month=iymdf[1];
    *day=iymdf[2];
    return;
  } else {
    sprintf(message,"slaDjcal could not process MJD: %f",djm);
    error_message(message);
  }
}

