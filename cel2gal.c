#include <math.h>
/* J2000 coordinates of galactic north pole and x-axis */
#define N_GAL_POLE_RA    (12 + (51 + (26.2754/60))/60)*15
#define N_GAL_POLE_DEC   (27 + (07 + (41.705/60))/60)
#define GAL_CENTER_DEC  -(28 + (56 + (10.219/60))/60)

void cel2gal(int rah, int ram, double ras,int decd,int decm,double decs,double *glon,double *glat) /*includefile*/
{
double ra, dec, gl, gb;
int isign=1;
double twopi, raddeg, spdec, cpdec, glong_off, ra_off;
double sra, cra, sdec, cdec;

    ra = (rah + (ram + (ras/60))/60)*15;
    if (decd != 0) isign = decd/abs(decd);
    dec = decd + isign*(decm + (decs/60))/60;

  twopi = 4*acos(0.0);
  raddeg = twopi/360.0;

  spdec = sin((90 - N_GAL_POLE_DEC)*raddeg);
  cpdec = cos((90 - N_GAL_POLE_DEC)*raddeg);
  glong_off = -asin(sin(GAL_CENTER_DEC*raddeg)/spdec)/raddeg;
  ra_off = 90 + N_GAL_POLE_RA;

  sra = sin((ra - ra_off)*raddeg);
  cra = cos((ra - ra_off)*raddeg);
  sdec = sin(dec*raddeg);
  cdec = cos(dec*raddeg);

  *glon = glong_off;
  *glon = glong_off + atan2((cdec*sra*cpdec + sdec*spdec), cdec*cra)/raddeg;
  if (*glon < 0)
    *glon = *glon + 360;
  *glat = asin(sdec*cpdec - cdec*sra*spdec)/raddeg;

}
  
