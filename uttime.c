#include <math.h>
/* return UT hours minutes and seconds given an MJD */
void uttime(double mjd, int *hh, int *mm, float *ss) /*includefile*/
{
  mjd=(mjd-floor(mjd))*24.0; /* get UT in hours */
  *hh=(int) mjd;
  mjd=(mjd-floor(mjd))*60.0; /* get remainder in minutes */
  *mm=(int) mjd;
  mjd=(mjd-floor(mjd))*60.0; /* get remainder in seconds */
  *ss=(float) mjd;
}
