/* standalone file for computing alfa position offsets */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 * Cactus File @(#)deg_trig.c	1.1
 *         SID 1.1
 *        Date 02/04/98
 */

static char SccsId[] = "@(#)deg_trig.c	1.1\t02/04/98";

#include <math.h>

#define PI		3.14159265358979323846
#define DEG_TO_RAD	PI/180.0
#define RAD_TO_DEG	180.0/PI

double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

/* Calculates trigonometric fn's with 'angle' in degrees, not radians.  */

double deg_sin(angle)
double angle;
{
  double sin();

  return( sin(DEG_TO_RAD * angle) );
}

double deg_cos(angle)
double angle;
{
  double cos();

  return( cos(DEG_TO_RAD * angle) );
}

double deg_tan(angle)
double  angle;
{
  double tan();

  return( tan(DEG_TO_RAD * angle) );
}

double deg_asin(value)
double  value;
{
  double asin();

  return( RAD_TO_DEG*asin(value) );
}

double deg_acos(value)
double  value;
{
  double acos();

  return( RAD_TO_DEG*acos(value) );
}

double deg_atan(value)
double  value;
{
  double atan();

  return( RAD_TO_DEG*atan(value) );
}



/*
   precession and nutaton routines converted to 'C'
   These are taken from the not so latest SLALIB sources - may 1992
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "prec_proto.h"

double fmod(), sin(), cos(), sqrt(), atan2();
double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();
double sla_epco(), sla_epj2d();
double current_time();
static double eqeq = 0.0;
int sla_deuler ( char *, double, double, double, double rmat[3][3]);
int sla_prec ( double,double, double rmatp[3][3]);

#define max(a,b) (((a)>(b))?(a):(b))
   
/* rambo usage: precess( epoch_ra, epoch_dec, &curr_ra, &curr_dec, when) */

double precess( epoch_ra, epoch_dec, ret_ra, ret_dec, curtime, src )
double epoch_ra, epoch_dec, *ret_ra, *ret_dec, curtime;
int src;
{
  double curr_ra, curr_dec, current_time();
  double d2r = PI/180.0;
  double h2r = 2.0*PI/24.0;

  curr_ra = epoch_ra*h2r; 
  curr_dec = epoch_dec*d2r;

  if( src == 'J' ) {
    sla_preces("FK5",sla_epco("J","J",2000.0),2000.0,
                &curr_ra, &curr_dec);
  } else {
/*              Proper motions not supplied */

/*              Precess to B1950 */
                sla_preces("FK4",sla_epco("B","B",1950.0),1950.0,
                  &curr_ra, &curr_dec);

/*              Add E-terms to make FK4 position */
                sla_addet(curr_ra,curr_dec,1950.0,&curr_ra,&curr_dec); 

/*              Convert to J2000 FK5 without proper motion */
                sla_fk45z(curr_ra,curr_dec,
                    sla_epco("J","B",1950.0),&curr_ra,&curr_dec);
  }

/*           Apparent */

/*           Convert to geocentric apparent  */
  sla_map(curr_ra,curr_dec,0.0,0.0,0.0,0.0,2000.0,
    sla_epj2d(sla_epco("J","J",curtime)), &curr_ra, &curr_dec );

  curr_ra = curr_ra/h2r;
  curr_dec = curr_dec/d2r;

  *ret_ra = curr_ra;
  *ret_dec = curr_dec;

  return(eqeq);
}

/*
*     - - - - - - -
*      P R E C E S
*     - - - - - - -
*
*  Precession - either FK4 (Bessel-Newcomb, pre IAU 1976) or
*  FK5 (Fricke, post IAU 1976) as required.
*
*  Given:
*     SYSTEM     char   precession to be applied: 'FK4' or 'FK5'
*     EP0,EP1    dp     starting and ending epoch
*     RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP0
*
*  Returned:
*     RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP1
*
*  Called:    sla_DRANRM, sla_PREBN, sla_PREC, sla_DCS2C,
*             sla_DMXV, sla_DCC2S
*
*  Notes:
*
*     1)  Lowercase characters in SYSTEM are acceptable.
*
*     2)  The epochs are Besselian if SYSTEM='FK4' and Julian if 'FK5'.
*         For example, to precess coordinates in the old system from
*         equinox 1900.0 to 1950.0 the call would be:
*             CALL sla_PRECES ('FK4', 1900D0, 1950D0, RA, DC)
*
*     3)  This routine will NOT correctly convert between the old and
*         the new systems - for example conversion from B1950 to J2000.
*         For these purposes see sla_FK425, sla_FK524, sla_FK45Z and
*         sla_FK54Z.
*
*     4)  If an invalid SYSTEM is supplied, values of -99D0,-99D0 will
*         be returned for both RA and DC.
*
*  P.T.Wallace   Starlink   20 April 1990
*/
sla_preces (system, ep0, ep1, ra, dc)
char *system;
double ep0,ep1,*ra,*dc;
{

  double pm[3][3],v1[3],v2[3];
  double sla_dranrm();

/*  Convert to uppercase and validate SYSTEM */
    if( strcasecmp(system, "FK4")!=0 && strcasecmp(system, "FK5")!=0 )
      return;

/*     Generate appropriate precession matrix */
    if(strcasecmp(system, "FK4") == 0)
      sla_prebn(ep0,ep1,pm);
    else
      sla_prec(ep0,ep1,pm);

/*     Convert RA,Dec to x,y,z */
    sla_dcs2c(*ra,*dc,v1);

/*     Precess */
    sla_dmxv(pm,v1,v2);

/*     Back to RA,Dec */
    sla_dcc2s(v2,ra,dc);
    *ra=sla_dranrm(*ra);
}

/*
*     - - - -
*      N U T
*     - - - -
*
*  Form the matrix of nutation for a given date - IAU 1980 theory
*  (double precision)
*
*  References:
*     Final report of the IAU Working Group on Nutation,
*      chairman P.K.Seidelmann, 1980.
*     Kaplan,G.H., 1981, USNO circular no. 163, pA3-6.
*
*  Given:
*     DATE   dp         TDB (loosely ET) as Modified Julian Date
*                                           (=JD-2400000.5)
*  Returned:
*     RMATN  dp(3,3)    nutation matrix
*
*  The matrix is in the sense   V(true)  =  RMATN * V(mean) .
*
*  Called:   sla_NUTC, sla_DEULER
*
*  P.T.Wallace   Starlink   10 May 1990
*/

sla_nut (date, rmatn)
double date,rmatn[3][3];
{
  double dpsi,deps,eps0;

/*  Nutation components and mean obliquity */
  sla_nutc(date,&dpsi,&deps,&eps0);

/*  Rotation matrix */
  sla_deuler("xzx",eps0,-dpsi,-(eps0+deps),rmatn);
}

/*
*     - - - - -
*      N U T C
*     - - - - -
*
*  Nutation:  longitude & obliquity components and mean
*  obliquity - IAU 1980 theory (double precision)
*
*  Given:
*
*     DATE        dp    TDB (loosely ET) as Modified Julian Date
*                                            (JD-2400000.5)
*  Returned:
*
*     DPSI,DEPS   dp    nutation in longitude,obliquity
*     EPS0        dp    mean obliquity
*
*  References:
*     Final report of the IAU Working Group on Nutation,
*      chairman P.K.Seidelmann, 1980.
*     Kaplan,G.H., 1981, USNO circular no. 163, pA3-6.
*
*  P.T.Wallace   Starlink   September 1987
*/
sla_nutc (date, dpsi, deps, eps0)
double date,*dpsi,*deps,*eps0;
{
double t,el,el2,el3;
double elp,elp2;
double f,f2,f4;
double d,d2,d4;
double om,om2;
double dp,de;
double a;

/*  Turns to arc seconds */
static double t2as=1296000.0;
/*  Arc seconds to radians */
static double as2r=0.4848136811095359949e-05;
/*  Units of 0.0001 arcsec to radians */
static double u2r=0.4848136811095359949e-05/1.0e4;

/*  Interval between basic epoch J2000.0 and current epoch (JC) */
      t=(date-51544.5)/36525.0;

/*
*  FUNDAMENTAL ARGUMENTS in the FK5 reference system
*/

/*  Mean longitude of the moon minus mean longitude of the moon's perigee */
      el=as2r*(485866.733+(1325*t2as+715922.633
              +(31.310+0.064*t)*t)*t);

/*  Mean longitude of the sun minus mean longitude of the sun's perigee */
      elp=as2r*(1287099.804+(99.0*t2as+1292581.224 +(-0.577-0.012*t)*t)*t);

/*  Mean longitude of the moon minus mean longitude of the moon's node */
      f=as2r*(335778.877+(1342*t2as+295263.137 +(-13.257+0.011*t)*t)*t);

/*  Mean elongation of the moon from the sun */
      d=as2r*(1072261.307+(1236*t2as+1105601.328 +(-6.891+0.019*t)*t)*t);

/*  Longitude of the mean ascending node of the lunar orbit on the */
/*   ecliptic, measured from the mean equinox of date */
      om=as2r*(450160.280+(-5.0*t2as-482890.539 +(7.455+0.008*t)*t)*t);

/*  Multiples of arguments */
      el2=el+el;
      el3=el2+el;
      elp2=elp+elp;
      f2=f+f;
      f4=f2+f2;
      d2=d+d;
      d4=d2+d2;
      om2=om+om;


/*
*  SERIES FOR THE NUTATION
*/
      dp=0.0;
      de=0.0;

/*  106 */
      dp=dp+sin(elp+d);
/*  105*/
      dp=dp-sin(f2+d4+om2);
/*  104*/
      dp=dp+sin(el2+d2);
/*  103*/
      dp=dp-sin(el-f2+d2);
/*  102*/
      dp=dp-sin(el+elp-d2+om);
/*  101*/
      dp=dp-sin(-elp+f2+om);
/*  100*/
      dp=dp-sin(el-f2-d2);
/*  99*/
      dp=dp-sin(elp+d2);
/*  98*/
      dp=dp-sin(f2-d+om2);
/*  97*/
      dp=dp-sin(-f2+om);
/*  96*/
      dp=dp+sin(-el-elp+d2+om);
/*  95*/
      dp=dp+sin(elp+f2+om);
/*  94*/
      dp=dp-sin(el+f2-d2);
/*  93*/
      dp=dp+sin(el3+f2-d2+om2);
/*  92*/
      dp=dp+sin(f4-d2+om2);
/*  91*/
      dp=dp-sin(el+d2+om);
/*  90*/
      dp=dp-sin(el2+f2+d2+om2);
/*  89*/
      a=el2+f2-d2+om;
      dp=dp+sin(a);
      de=de-cos(a);
/*  88*/
      dp=dp+sin(el-elp-d2);
/*  87*/
      dp=dp+sin(-el+f4+om2);
/*  86*/
      a=-el2+f2+d4+om2;
      dp=dp-sin(a);
      de=de+cos(a);
/*  85*/
      a=el+f2+d2+om;
      dp=dp-sin(a);
      de=de+cos(a);
/*  84*/
      a=el+elp+f2-d2+om2;
      dp=dp+sin(a);
      de=de-cos(a);
/*  83*/
      dp=dp-sin(el2-d4);
/*  82*/
      a=-el+f2+d4+om2;
      dp=dp-2.0*sin(a);
      de=de+cos(a);
/*  81*/
      a=-el2+f2+d2+om2;
      dp=dp+sin(a);
      de=de-cos(a);
/*  80*/
      dp=dp-sin(el-d4);
/*  79*/
      a=-el+om2;
      dp=dp+sin(a);
      de=de-cos(a);
/*  78*/
      a=f2+d+om2;
      dp=dp+2.0*sin(a);
      de=de-cos(a);
/*  77*/
      dp=dp+2.0*sin(el3);
/*  76*/
      a=el+om2;
      dp=dp-2.0*sin(a);
      de=de+cos(a);
/*  75*/
      a=el2+om;
      dp=dp+2.0*sin(a);
      de=de-cos(a);
/*  74*/
      a=-el+f2-d2+om;
      dp=dp-2.0*sin(a);
      de=de+cos(a);
/*  73*/
      a=el+elp+f2+om2;
      dp=dp+2.0*sin(a);
      de=de-cos(a);
/*  72*/
      a=-elp+f2+d2+om2;
      dp=dp-3.0*sin(a);
      de=de+cos(a);
/*  71*/
      a=el3+f2+om2;
      dp=dp-3.0*sin(a);
      de=de+cos(a);
/*  70*/
      a=-el2+om;
      dp=dp-2.0*sin(a);
      de=de+cos(a);
/*  69*/
      a=-el-elp+f2+d2+om2;
      dp=dp-3.0*sin(a);
      de=de+cos(a);
/*  68*/
      a=el-elp+f2+om2;
      dp=dp-3.0*sin(a);
      de=de+cos(a);
/*  67*/
      dp=dp+3.0*sin(el+f2);
/*  66*/
      dp=dp-3.0*sin(el+elp);
/*  65*/
      dp=dp-4.0*sin(d);
/*  64*/
      dp=dp+4.0*sin(el-f2);
/*  63*/
      dp=dp-4.0*sin(elp-d2);
/*  62*/
      a=el2+f2+om;
      dp=dp-5.0*sin(a);
      de=de+3.0*cos(a);
/*  61*/
      dp=dp+5.0*sin(el-elp);
/*  60*/
      a=-d2+om;
      dp=dp-5.0*sin(a);
      de=de+3.0*cos(a);
/*  59*/
      a=el+f2-d2+om;
      dp=dp+6.0*sin(a);
      de=de-3.0*cos(a);
/*  58*/
      a=f2+d2+om;
      dp=dp-7.0*sin(a);
      de=de+3.0*cos(a);
/*  57*/
      a=d2+om;
      dp=dp-6.0*sin(a);
      de=de+3.0*cos(a);
/*  56*/
      a=el2+f2-d2+om2;
      dp=dp+6.0*sin(a);
      de=de-3.0*cos(a);
/*  55*/
      dp=dp+6.0*sin(el+d2);
/*  54;*/
      a=el+f2+d2+om2;
      dp=dp-8.0*sin(a);
      de=de+3.0*cos(a);
/*  53*/
      a=-elp+f2+om2;
      dp=dp-7.0*sin(a);
      de=de+3.0*cos(a);
/*  52*/
      a=elp+f2+om2;
      dp=dp+7.0*sin(a);
      de=de-3.0*cos(a);
/*  51*/
      dp=dp-7.0*sin(el+elp-d2);
/*  50*/
      a=-el+f2+d2+om;
      dp=dp-10.0*sin(a);
      de=de+5.0*cos(a);
/*  49*/
      a=el-d2+om;
      dp=dp-13.0*sin(a);
      de=de+7.0*cos(a);
/*  48*/
      a=-el+d2+om;
      dp=dp+16.0*sin(a);
      de=de-8.0*cos(a);
/*  47*/
      a=-el+f2+om;
      dp=dp+21.0*sin(a);
      de=de-10.0*cos(a);
/*  46*/
      dp=dp+26.0*sin(f2);
      de=de-cos(f2);
/*  45*/
      a=el2+f2+om2;
      dp=dp-31.0*sin(a);
      de=de+13.0*cos(a);
/*  44*/
      a=el+f2-d2+om2;
      dp=dp+29.0*sin(a);
      de=de-12.0*cos(a);
/*  43*/
      dp=dp+29.0*sin(el2);
      de=de-cos(el2);
/*  42*/
      a=f2+d2+om2;
      dp=dp-38.0*sin(a);
      de=de+16.0*cos(a);
/*  41*/
      a=el+f2+om;
      dp=dp-51.0*sin(a);
      de=de+27.0*cos(a);
/*  40*/
      a=-el+f2+d2+om2;
      dp=dp-59.0*sin(a);
      de=de+26.0*cos(a);
/*  39*/
      a=-el+om;
      dp=dp+(-58.0-0.1*t)*sin(a);
      de=de+32.0*cos(a);
/*  38*/
      a=el+om;
      dp=dp+(63.0+0.1*t)*sin(a);
      de=de-33.0*cos(a);
/*  37*/
      dp=dp+63.0*sin(d2);
      de=de-2.0*cos(d2);
/*  36*/
      a=-el+f2+om2;
      dp=dp+123.0*sin(a);
      de=de-53.0*cos(a);
/*  35*/
      a=el-d2;
      dp=dp-158.0*sin(a);
      de=de-cos(a);
/*  34*/
      a=el+f2+om2;
      dp=dp-301.0*sin(a);
      de=de+(129.0-0.1*t)*cos(a);
/*  33*/
      a=f2+om;
      dp=dp+(-386.0-0.4*t)*sin(a);
      de=de+200.0*cos(a);
/*  32*/
      dp=dp+(712.0+0.1*t)*sin(el);
      de=de-7.0*cos(el);
/*  31*/
      a=f2+om2;
      dp=dp+(-2274.0-0.2*t)*sin(a);
      de=de+(977.0-0.5*t)*cos(a);
/*  30*/
      dp=dp-sin(elp+f2-d2);
/*  29*/
      dp=dp+sin(-el+d+om);
/*  28*/
      dp=dp+sin(elp+om2);
/*  27*/
      dp=dp-sin(elp-f2+d2);
/*  26*/
      dp=dp+sin(-f2+d2+om);
/*  25*/
      dp=dp+sin(el2+elp-d2);
/*  24*/
      dp=dp-4.0*sin(el-d);
/*  23*/
      a=elp+f2-d2+om;
      dp=dp+4.0*sin(a);
      de=de-2.0*cos(a);
/*  22*/
      a=el2-d2+om;
      dp=dp+4.0*sin(a);
      de=de-2.0*cos(a);
/*  21*/
      a=-elp+f2-d2+om;
      dp=dp-5.0*sin(a);
      de=de+3.0*cos(a);
/*  20*/
      a=-el2+d2+om;
      dp=dp-6.0*sin(a);
      de=de+3.0*cos(a);
/*  19*/
      a=-elp+om;
      dp=dp-12.0*sin(a);
      de=de+6.0*cos(a);
/*  18*/
      a=elp2+f2-d2+om2;
      dp=dp+(-16.0+0.1*t)*sin(a);
      de=de+7.0*cos(a);
/*  17*/
      a=elp+om;
      dp=dp-15.0*sin(a);
      de=de+9.0*cos(a);
/*  16*/
      dp=dp+(17.0-0.1*t)*sin(elp2);
/*  15*/
      dp=dp-22.0*sin(f2-d2);
/*  14*/
      a=el2-d2;
      dp=dp+48.0*sin(a);
      de=de+cos(a);
/*  13*/
      a=f2-d2+om;
      dp=dp+(129.0+0.1*t)*sin(a);
      de=de-70.0*cos(a);
/*  12*/
      a=-elp+f2-d2+om2;
      dp=dp+(217.0-0.5*t)*sin(a);
      de=de+(-95.0+0.3*t)*cos(a);
/*  11*/
      a=elp+f2-d2+om2;
      dp=dp+(-517.0+1.2*t)*sin(a);
      de=de+(224.0-0.6*t)*cos(a);
/*  10*/
      dp=dp+(1426.0-3.4*t)*sin(elp);
      de=de+(54.0-0.1*t)*cos(elp);
/*  9*/
      a=f2-d2+om2;
      dp=dp+(-13187.0-1.6*t)*sin(a);
      de=de+(5736.0-3.1*t)*cos(a);
/*  8*/
      dp=dp+sin(el2-f2+om);
/*  7*/
      a=-elp2+f2-d2+om;
      dp=dp-2.0*sin(a);
      de=de+1.0*cos(a);
/*  6*/
      dp=dp-3.0*sin(el-elp-d);
/*  5*/
      a=-el2+f2+om2;
      dp=dp-3.0*sin(a);
      de=de+1.0*cos(a);
/*  4*/
      dp=dp+11.0*sin(el2-f2);
/*  3*/
      a=-el2+f2+om;
      dp=dp+46.0*sin(a);
      de=de-24.0*cos(a);
/*  2*/
      dp=dp+(2062.0+0.2*t)*sin(om2);
      de=de+(-895.0+0.5*t)*cos(om2);
/*  1*/
      dp=dp+(-171996.0-174.2*t)*sin(om);
      de=de+(92025.0+8.9*t)*cos(om);

/*  Convert results to radians*/
      *dpsi=dp*u2r;
      *deps=de*u2r;

/*  Mean obliquity*/
      *eps0=as2r*(84381.448+ (-46.8150+ (-0.00059+ 0.001813*t)*t)*t);
      eqeq = *dpsi * cos(*eps0 + *deps) * 13750.987; 
}

/*
*     - - - - -
*      P R E C
*     - - - - -
*
*  Form the matrix of precession between two epochs (IAU 1976, FK5)
*  (double precision)
*
*  Given:
*     EP0    dp         beginning epoch
*     EP1    dp         ending epoch
*
*  Returned:
*     RMATP  dp(3,3)    precession matrix
*
*  Notes:
*
*     1)  The epochs are TDB (loosely ET) Julian epochs.
*
*     2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
*
*  References:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*      equations (6) & (7), p283.
*     Kaplan,G.H., 1981. USNO circular no. 163, pA2.
*
*  Called:  sla_DEULER
*
*  P.T.Wallace   Starlink   12 April 1990
*/
sla_prec (ep0, ep1, rmatp)
double ep0,ep1,rmatp[3][3];
{

/*  Arc seconds to radians*/
static double as2r=0.4848136811095359949e-05;

double t0,t,tas2r,w,zeta,z,theta;



/*  Interval between basic epoch J2000.0 and beginning epoch (JC)*/
      t0 = (ep0-2000.0)/100.0;

/*  Interval over which precession required (JC)*/
      t = (ep1-ep0)/100.0;

/*  Euler angles*/
      tas2r = t*as2r;
      w = 2306.2181+(1.39656-0.000139*t0)*t0;

      z = (w+((1.09468+0.000066*t0)+0.018203*t)*t)*tas2r;

      zeta = (w+((0.30188-0.000344*t0)+0.017998*t)*t)*tas2r;

      theta = (2004.3109+(-0.85330-0.000217*t0)*t0);
      theta = (theta +((-0.42665-0.000217*t0)-0.041833*t)*t);
      theta = theta*tas2r;

/*
      theta = ((2004.3109+(-0.85330-0.000217*t0)*t0)
             +((-0.42665-0.000217*t0)-0.041833*t)*t)*tas2r;

      printf("theta %f\n", theta );
*/

/*  Rotation matrix*/
      sla_deuler("zyz",-zeta,theta,-z,rmatp);
}

/*
*     - - - - - - -
*      P R E N U T
*     - - - - - - -
*
*  Form the matrix of precession and nutation (IAU1976/FK5)
*  (double precision)
*
*  Given:
*     EPOCH   dp         Julian Epoch for mean coordinates
*     DATE    dp         Modified Julian Date (JD-2400000.5)
*                        for true coordinates
*
*  Returned:
*     RMATPN  dp(3,3)    combined precession/nutation matrix
*
*  Called:  sla_PREC, sla_EPJ, sla_NUT, sla_DMXM
*
*  Notes:
*
*  1)  The epoch and date are TDB (loosely ET).
*
*  2)  The matrix is in the sense   V(true)  =  RMATPN * V(mean)
*
*  P.T.Wallace   Starlink   April 1987
*/
sla_prenut (epoch, date, rmatpn)
double epoch,date,rmatpn[3][3];
{

double rmatp[3][3],rmatn[3][3],sla_epj();



/*  Precession*/
      sla_prec(epoch,sla_epj(date),rmatp);

/*  Nutation*/
      sla_nut(date,rmatn);

/*  Combine the matrices:  PN = N x P*/
      sla_dmxm(rmatn,rmatp,rmatpn);
}

/*
*     - - - - - - 
*      P R E B N
*     - - - - - -
*
*  Generate the matrix of precession between two epochs,
*  using the old, pre IAU 1976, Bessel-Newcomb model, in
*  Andoyer's formulation (double precision)
*
*  Given:
*     BEP0    dp         beginning Besselian epoch
*     BEP1    dp         ending Besselian epoch
*
*  Returned:
*     RMATP  dp(3,3)    precession matrix
*
*  The matrix is in the sense   V(BEP1)  =  RMATP * V(BEP0) .
*
*  Reference:
*     Smith et al 1989, Astron. J., 97, 1 p269.
*
*  Called:  sla_DEULER
*
*  P.T.Wallace   Starlink   12 April 1990
*/
sla_prebn (bep0, bep1, rmatp)
double bep0,bep1,rmatp[3][3];
{

/*  Arc seconds to radians*/
static double as2r=0.4848136811095359949e-05;

double bigt,t,tas2r,w,zeta,z,theta;



/*  Interval between basic epoch B1850.0 and beginning epoch in TC*/
      bigt = (bep0-1850.0)/100.0;

/*  Interval over which precession required, in tropical centuries*/
      t = (bep1-bep0)/100.0;

/*  Euler angles*/
      tas2r = t*as2r;
      w = 2303.5545+(1.39720+0.000060*bigt)*bigt;

      zeta = (w+(0.30240-0.000270*bigt+0.017995*t)*t)*tas2r;
      z = (w+(1.09480+0.000390*bigt+0.018325*t)*t)*tas2r;
      theta = (2005.112+(-0.8529-0.00037*bigt)*bigt+
             (-0.4265-0.00037*bigt-0.04180*t)*t)*tas2r;

/*  Rotation matrix*/
      sla_deuler("zyz",-zeta,theta,-z,rmatp);
}


/*
*     - - - - - - -
*      D R A N R M
*     - - - - - - -
*
*  Normalise angle into range 0-2 pi  (double precision)
*
*  Given:
*     ANGLE     dp      the angle in radians
*
*  The result is ANGLE expressed in the range 0-2 pi (double
*  precision).
*
*  P.T.Wallace   Starlink   December 1984
*/
double sla_dranrm(angle)
double angle;
{
  double ret, fmod();
  static double d2pi=6.283185307179586476925287;

  ret =fmod(angle,d2pi);
  if(ret < 0.0 )
     ret=ret+d2pi;
  return(ret);
}

/*
*     - - - - - -
*      D C S 2 C
*     - - - - - -
*
*  Spherical coordinates to direction cosines (double precision)
*
*  Given:
*     A,B       dp      spherical coordinates in radians
*                        (RA,Dec), (Long,Lat) etc
*
*  Returned:
*     V         dp(3)   x,y,z unit vector
*
*  The spherical coordinates are longitude (+ve anticlockwise
*  looking from the +ve latitude pole) and latitude.  The
*  Cartesian coordinates are right handed, with the x axis
*  at zero longitude and latitude, and the z axis at the
*  +ve latitude pole.
*
*  P.T.Wallace   Starlink   October 1984
*/
sla_dcs2c (a, b, v)
double a,b,v[3];
{
  double cosb;

  cosb=cos(b);

  v[0]=cos(a)*cosb;
  v[1]=sin(a)*cosb;
  v[2]=sin(b);
}

/*
*     - - - -
*      E P J
*     - - - -
*
*  Conversion of Modified Julian Date to Julian Epoch (double precision)
*
*  Given:
*     DATE     dp       Modified Julian Date (JD - 2400000.5)
*
*  The result is the Julian Epoch.
*
*  Reference:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*
*  P.T.Wallace   Starlink   February 1984
*/
double sla_epj(date)
double date;
{
      return( 2000.0 + (date-51544.5)/365.25 );
}

/*
*     - - - - - 
*      D M X M
*     - - - - -
* 
*  Product of two 3x3 matrices:
*
*      matrix C  =  matrix A  x  matrix B
*
*  (double precision)
*
*  Given:
*      A      dp(3,3)        matrix
*      B      dp(3,3)        matrix
*
*  Returned:
*      C      dp(3,3)        matrix result
*
*  To comply with the ANSI Fortran 77 standard, A, B and C must
*  be different arrays.  However, the routine is coded so as to
*  work properly on the VAX and many other systems even if this
*  rule is violated.
*
*  P.T.Wallace   Starlink   5 April 1990
*/
sla_dmxm (a, b, c)
double a[3][3],b[3][3],c[3][3];
{
      int i,j,k;
      double w,wm[3][3];

/*  Multiply into scratch matrix */
      for(i=0; i< 3; i++ ) {
        for(j=0; j< 3; j++ ) {
            w=0.0;
            for( k=0; k<3; k++ ) 
               w=w+a[k][i]*b[j][k];
            wm[j][i]=w;
         }
      }

/*  Return the result */
      for(j=0; j< 3; j++ )
        for(i=0; i< 3; i++ )
            c[j][i]=wm[j][i];
}

/*
*     - - - - - - -
*      D E U L E R
*     - - - - - - -
*
*  Form a rotation matrix from the Euler angles - three successive
*  rotations about specified Cartesian axes (double precision)
*
*  Given:
*    ORDER  c*(*)    specifies about which axes the rotations occur
*    PHI    dp       1st rotation (radians)
*    THETA  dp       2nd rotation (   "   )
*    PSI    dp       3rd rotation (   "   )
*
*  Returned:
*    RMAT   dp(3,3)  rotation matrix
*
*  A rotation is positive when the reference frame rotates
*  anticlockwise as seen looking towards the origin from the
*  positive region of the specified axis.
*
*  The characters of ORDER define which axes the three successive
*  rotations are about.  A typical value is 'ZXZ', indicating that
*  RMAT is to become the direction cosine matrix corresponding to
*  rotations of the reference frame through PHI radians about the
*  old Z-axis, followed by THETA radians about the resulting X-axis,
*  then PSI radians about the resulting Z-axis.
*
*  The axis names can be any of the following, in any order or
*  combination:  X, Y, Z, uppercase or lowercase, 1, 2, 3.  Normal
*  axis labelling/numbering conventions apply;  the xyz (=123)
*  triad is right-handed.  Thus, the 'ZXZ' example given above
*  could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ').  ORDER
*  is terminated by length or by the first unrecognised character.
*
*  Fewer than three rotations are acceptable, in which case the later
*  angle arguments are ignored.  Zero rotations produces a unit RMAT.
*
*  P.T.Wallace   Starlink   November 1988
*/
sla_deuler (order, phi, theta, psi, rmat)
char *order;
double phi,theta,psi,rmat[3][3];
{

      int j,i,l,n,k;
      double result[3][3],rotn[3][3],angle,s,c,w,wm[3][3];
      char axis;

/*  Initialise result matrix */
      for(j=0; j<3; j++ ) {
        for(i=0; i<3; i++ ) {
            if(i != j)
               result[j][i] = 0.0;
            else
               result[j][i] = 1.0;
         }
      }

/*  Establish length of axis string */
      l = strlen(order) -1;

/*  Look at each character of axis string until finished */
      for(n=0; n<3; n++ ) {
         if(n<=l) {

/*        Initialise rotation matrix for the current rotation */
            for(j=0; j<3; j++ ) {
               for(i=0; i<3; i++ ) {
                  if (i!=j)
                     rotn[j][i] = 0.0;
                  else
                     rotn[j][i] = 1.0;
               }
            }

/*        Pick up the appropriate Euler angle and take sine & cosine */
            if (n==0)
               angle = phi;
            else if (n==1)
               angle = theta;
            else
               angle = psi;
            s = sin(angle);
            c = cos(angle);

/*        Identify the axis */
            axis = order[n];
            if(axis=='X'|| axis == 'x' || axis=='1') {

/*           Matrix for x-rotation */
               rotn[1][1] = c;
               rotn[2][1] = s;
               rotn[1][2] = -s;
               rotn[2][2] = c;

            } else if (axis=='Y'|| axis=='y'|| axis=='2') {

/*           Matrix for y-rotation */
               rotn[0][0] = c;
               rotn[2][0] = -s;
               rotn[0][2] = s;
               rotn[2][2] = c;

            } else if(axis=='Z'|| axis=='z'||axis=='3') {

/*           Matrix for z-rotation */
               rotn[0][0] = c;
               rotn[1][0] = s;
               rotn[0][1] = -s;
               rotn[1][1] = c;
            } else {

/*           Unrecognised character - fake end of string */
               l = 0;

            }

/*        Apply the current rotation (matrix ROTN x matrix RESULT) */
            for(i=0; i<3; i++ ) {
              for(j=0; j<3; j++ ) {
                  w = 0.0;
                  for(k=0; k<3; k++ )
                     w = w+rotn[k][i]*result[j][k];
                  wm[j][i] = w;
               } 
            }
            for(j=0; j<3; j++ )
              for(i=0; i<3; i++ )
                  result[j][i] = wm[j][i];
         }
      }

/*  Copy the result */
      for(j=0; j<3; j++ )
        for(i=0; i<3; i++ )
            rmat[j][i] = result[j][i];
}

/*
*     - - - - - -
*      D C C 2 S
*     - - - - - -
*
*  Direction cosines to spherical coordinates (double precision)
*
*  Given:
*     V     d(3)   x,y,z vector
*
*  Returned:
*     A,B   d      spherical coordinates in radians
*
*  The spherical coordinates are longitude (+ve anticlockwise
*  looking from the +ve latitude pole) and latitude.  The
*  Cartesian coordinates are right handed, with the x axis
*  at zero longitude and latitude, and the z axis at the
*  +ve latitude pole.
*
*  If V is null, zero A and B are returned.
*  At either pole, zero A is returned.
*
*  P.T.Wallace   Starlink   July 1989
*/
sla_dcc2s (v, a, b)
double v[3], *a,*b;
{
  double x,y,z,r;

      x = v[0];
      y = v[1];
      z = v[2];
      r = sqrt(x*x+y*y);

      if (r == 0.0)
        *a = 0.0;
      else
        *a = atan2(y,x);

      if(z == 0.0)
         *b = 0.0;
      else
         *b = atan2(z,r);
}

/*
*     - - - - - 
*      D M X V
*     - - - - -
*
*  Performs the 3-D forward unitary transformation:
*
*     vector VB = matrix DM * vector VA
*
*  (double precision)
*
*  Given:
*     DM       dp(3,3)    matrix
*     VA       dp(3)      vector
*
*  Returned:
*     VB       dp(3)      result vector
*
*  P.T.Wallace   Starlink   March 1986
*/
sla_dmxv (dm, va, vb)
double dm[3][3],va[3],vb[3];
{
      int i,j;
      double w,vw[3];


/*  Matrix DM * vector VA -> vector VW */
      for(j=0; j<3; j++ ) {
         w=0.0;
         for(i=0; i<3; i++ ) 
            w=w+dm[i][j]*va[i];
         vw[j]=w;
      }

/*  Vector VW -> vector VB */
      for(j=0; j<3; j++ )
         vb[j]=vw[j];
}


/*+
*     - - - - - -
*      A D D E T
*     - - - - - -
*
*  Add the E-terms (elliptic component of annual aberration)
*  to a pre IAU 1976 mean place to conform to the old
*  catalogue convention (double precision)
*
*  Given:
*     RM,DM     dp     RA,Dec (radians) without E-terms
*     EQ        dp     Besselian epoch of mean equator and equinox
*
*  Returned:
*     RC,DC     dp     RA,Dec (radians) with E-terms included
*
*  Called:
*     sla_ETRMS, sla_DCS2C, sla_DCC2S, sla_DRANRM, sla_DRANGE
*
*  Explanation:
*     Most star positions from pre-1984 optical catalogues (or
*     derived from astrometry using such stars) embody the
*     E-terms.  If it is necessary to convert a formal mean
*     place (for example a pulsar timing position) to one
*     consistent with such a star catalogue, then the RA,Dec
*     should be adjusted using this routine.
*
*  Reference:
*     Explanatory Supplement to the Astronomical Ephemeris,
*     section 2D, page 48.
*
*  P.T.Wallace   Starlink   July 1986
*/
sla_addet (rm, dm, eq, rc, dc)
double rm,dm,eq, *rc, *dc;
{
      double a[3],v[3];
      int i;

/*  E-terms vector */
      sla_etrms(eq,a);

/*  Spherical to Cartesian */
      sla_dcs2c(rm,dm,v);

/*  Include the E-terms */
      for(i=0; i<3; i++ )
         v[i]=v[i]+a[i];

/*  Cartesian to spherical */
      sla_dcc2s(v,rc,dc);

/*  Bring RA into conventional range */
      *rc=sla_dranrm(*rc);
}

/*
*     - - - - - -
*      E T R M S
*     - - - - - -
*
*  Compute the E-terms (elliptic component of annual aberration)
*  vector (double precision)
*
*  Given:
*     EP      dp      Besselian epoch
*
*  Returned:
*     EV      dp(3)   E-terms as (dx,dy,dz)
*
*  Note the use of the J2000 aberration constant (20.49552 arcsec).
*  This is a reflection of the fact that the E-terms embodied in
*  existing star catalogues were computed from a variety of
*  aberration constants.  Rather than adopting one of the old
*  constants the latest value is used here.
*
*  References:
*     1  Smith, C.A. et al., 1989.  Astr.J. 97, 265.
*     2  Yallop, B.D. et al., 1989.  Astr.J. 97, 274.
*
*  P.T.Wallace   Starlink   10 April 1990
*/
sla_etrms (ep, ev)
double ep,ev[3];
{

/*  Arcseconds to radians */
static double as2r=0.4848136811095359949e-5;

double t,e,e0,p,ek,cp;

/*  Julian centuries since B1950 */
      t=(ep-1950.0)*1.00002135903e-2;

/*  Eccentricity */
      e=0.01673011-(0.00004193+0.000000126*t)*t;

/*  Mean obliquity */
      e0=(84404.836-(46.8495+(0.00319+0.00181*t)*t)*t)*as2r;

/*  Mean longitude of perihelion */
      p=(1015489.951+(6190.67+(1.65+0.012*t)*t)*t)*as2r;

/*  E-terms */
      ek=e*20.49552*as2r;
      cp=cos(p);
      ev[0]= ek*sin(p);
      ev[1]=-ek*cp*cos(e0);
      ev[2]=-ek*cp*sin(e0);
}
/*
*     - - - - - -
*      F K 4 5 Z
*     - - - - - -
*
*  Convert B1950.0 FK4 star data to J2000.0 FK5 assuming zero
*  proper motion in an inertial frame (double precision)
*
*  This routine converts stars from the old, Bessel-Newcomb, FK4
*  system to the new, IAU 1976, FK5, Fricke system, in such a
*  way that the FK5 proper motion is zero.  Because such a star
*  has, in general, a non-zero proper motion in the FK4 system,
*  the routine requires the epoch at which the position in the
*  FK4 system was determined.
*
*  The method is from Appendix 2 of ref 1, but using the constants
*  of ref 2.
*
*  Given:
*     R1950,D1950     dp    B1950.0 FK4 RA,Dec at epoch (rad)
*     BEPOCH          dp    Besselian epoch (e.g. 1979.3D0)
*
*  Returned:
*     R2000,D2000     dp    J2000.0 FK5 RA,Dec (rad)
*
*  Notes:
*
*  1)  The epoch BEPOCH is strictly speaking Besselian, but
*      if a Julian epoch is supplied the result will be
*      affected only to a negligible extent.
*
*  2)  Conversion from Besselian epoch 1950.0 to Julian epoch
*      2000.0 only is provided for.  Conversions involving other
*      epochs will require use of the appropriate precession,
*      proper motion, and E-terms routines before and/or
*      after FK425 is called.
*
*  3)  In the FK4 catalogue the proper motions of stars within
*      10 degrees of the poles do not embody the differential
*      E-term effect and should, strictly speaking, be handled
*      in a different manner from stars outside these regions.
*      However, given the general lack of homogeneity of the star
*      data available for routine astrometry, the difficulties of
*      handling positions that may have been determined from
*      astrometric fields spanning the polar and non-polar regions,
*      the likelihood that the differential E-terms effect was not
*      taken into account when allowing for proper motion in past
*      astrometry, and the undesirability of a discontinuity in
*      the algorithm, the decision has been made in this routine to
*      include the effect of differential E-terms on the proper
*      motions for all stars, whether polar or not.  At epoch 2000,
*      and measuring on the sky rather than in terms of dRA, the
*      errors resulting from this simplification are less than
*      1 milliarcsecond in position and 1 milliarcsecond per
*      century in proper motion.
*
*  References:
*
*     1  Aoki,S., et al, 1983.  Astron.Astrophys., 128, 263.
*
*     2  Smith, C.A. et al, 1989.  "The transformation of astrometric
*        catalog systems to the equinox J2000.0".  Astron.J. 97, 265.
*
*     3  Yallop, B.D. etal, 1989.  "Transformation of mean star places
*        from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
*        Astron.J. 97, 274.
*
*  Called:  sla_DCS2C, sla_EPJ, sla_EPB2D, sla_DCC2S, sla_DRANRM
*
*  P.T.Wallace   Starlink   10 April 1990
*/

sla_fk45z (r1950,d1950,bepoch,r2000,d2000)
double r1950,d1950,bepoch, *r2000,*d2000;
{
      static double d2pi=6.283185307179586476925287;
      double w;
      int i,j;

/*  Position and position+velocity vectors */
      double r0[3],a1[3],v1[3],v2[6];

/*  Radians per year to arcsec per century */
      static double pmf=100.0*60.0*60.0*360.0/6.283185307179586476925287;

/*  Functions */
      double sla_epj(),sla_epb2d(),sla_dranrm();

/* */
/*  CANONICAL CONSTANTS  (see references) */
/* */

/*  Vectors A and Adt, and matrix M (only half of which is needed here) */
     static double  a[3] = { -1.62557e-6,  -0.31919e-6, -0.13843e-6 };
     static double ad[3] = { +1.245e-3,    -1.580e-3,   -0.659e-3   };

     static double em[3][6] = {
                         { +0.999925678186902,
                            +0.011182059571766,
                            +0.004857946721186,
                            -0.000541652366951,
                            +0.237917612131583,
                            -0.436111276039270 
                         },{
                            -0.011182059642247,
                             +0.999937478448132,
                             -0.000027147426498,
                             -0.237968129744288,
                             -0.002660763319071,
                             +0.012259092261564 
                         } , {
                            -0.004857946558960,
                             -0.000027176441185,
                             +0.999988199738770,
                             +0.436227555856097,
                             -0.008537771074048,
                             +0.002119110818172 
                          }
                        };

/*  Spherical to Cartesian */
      sla_dcs2c(r1950,d1950,r0);

/*  Adjust vector A to give zero proper motion in FK5 */
      w=(bepoch-1950.0)/pmf;
      for(i=0; i<3; i++ )
         a1[i]=a[i]+w*ad[i];

/*  Remove e-terms */
      w=r0[0]*a1[0]+r0[1]*a1[1]+r0[2]*a1[2];
      for(i=0; i<3; i++ )
         v1[i]=r0[i]-a1[i]+w*r0[i];

/*  Convert position vector to Fricke system */
      for(i=0; i<6; i++ ) {
         w=0.0;
         for(j=0; j<3; j++ )
            w=w+em[j][i]*v1[j];
         v2[i]=w;
      }

/*  Allow for fictitious proper motion in FK4 */
      w=(sla_epj(sla_epb2d(bepoch))-2000.0)/pmf;
      for(i=0; i<3; i++ )
         v2[i]=v2[i]+w*v2[i+3];

/*  Revert to spherical coordinates */
      sla_dcc2s(v2,&w,d2000);
      *r2000=sla_dranrm(w);
}

/*
*     - - - -
*      M A P
*     - - - -
*
*  Transform star RA,Dec from mean place to geocentric apparent
*
*  The reference frames and timescales used are post IAU 1976.
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Given:
*     RM,DM    dp     mean RA,Dec (rad)
*     PR,PD    dp     proper motions:  RA,Dec changes per Julian year
*     PX       dp     parallax (arcsec)
*     RV       dp     radial velocity (km/sec, +ve if receding)
*     EQ       dp     epoch and equinox of star data (Julian)
*     DATE     dp     TDB for apparent place (JD-2400000.5)
*
*  Returned:
*     RA,DA    dp     apparent RA,Dec (rad)
*
*  Called:
*     sla_MAPPA       star-independent parameters
*     sla_MAPQK       quick mean to apparent
*
*  Notes:
*
*  1)  EQ is the Julian epoch specifying both the reference
*      frame and the epoch of the position - usually 2000.
*      For positions where the epoch and equinox are
*      different, use the routine sla_PM to apply proper
*      motion corrections before using this routine.
*
*  2)  The distinction between the required TDB and TDT is
*      always negligible.  Moreover, for all but the most
*      critical applications UTC is adequate.
*
*  3)  The proper motions in RA are dRA/dt rather than
*      cos(Dec)*dRA/dt.
*
*  4)  This routine may be wasteful for some applications
*      because it recomputes the Earth position/velocity and
*      the precession/nutation matrix each time, and because
*      it allows for parallax and proper motion.  Where
*      multiple transformations are to be carried out for one
*      epoch, a faster method is to call the sla_MAPPA routine
*      once and then either the sla_MAPQK routine (which includes
*      parallax and proper motion) or sla_MAPQKZ (which assumes
*      zero parallax and proper motion).
*
*  P.T.Wallace   Starlink   11 April 1990
*/
sla_map (rm, dm, pr, pd, px, rv, eq, date, ra, da)
double rm,dm,pr,pd,px,rv,eq,date, *ra, *da;
{
  double amprms[21];

/*  Star-independent parameters */
  sla_mappa(eq,date,amprms);

/*  Mean to apparent */
      sla_mapqk(rm,dm,pr,pd,px,rv,amprms,ra,da);
}

/*
*     - - - - - -
*      M A P P A
*     - - - - - -
*
*  Compute star-independent parameters in preparation for
*  conversions between mean place and geocentric apparent place.
*
*  The parameters produced by this routine are required in the
*  parallax, light deflection, aberration, and precession/nutation
*  parts of the mean/apparent transformations.
*
*  The reference frames and timescales used are post IAU 1976.
*
*  Given:
*     EQ       d      epoch of mean equinox to be used (Julian)
*     DATE     d      TDB (JD-2400000.5)
*
*  Returned:
*     AMPRMS   d(21)  star-independent mean-to-apparent parameters:
*             
*       (1)      time interval for proper motion (Julian years)
*       (2-4)    barycentric position of the Earth (AU)
*       (5-7)    heliocentric direction of the Earth (unit vector)
*       (8)      (grav rad Sun)*2/(Sun-Earth distance)
*       (9-11)   ABV: barycentric Earth velocity in units of c
*       (12)     sqrt(1-v**2) where v=modulus(ABV)
*       (13-21)  precession/nutation (3,3) matrix
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Notes:
*
*  1)  For DATE, the distinction between the required TDB and TDT
*      is always negligible.  Moreover, for all but the most
*      critical applications UTC is adequate.
*
*  2)  The accuracy of the routines using the parameters AMPRMS is
*      limited by the routine sla_EVP, used here to compute the
*      Earth position and velocity by the methods of Stumpff.
*      The maximum error in the resulting aberration corrections is
*      about 0.3 milliarcsecond.
*
*  3)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to
*      the mean equinox and equator of epoch EQ.
*
*  4)  The parameters AMPRMS produced by this routine are used by
*      sla_MAPQK and sla_MAPQKZ.
*
*  5)  This routine replaces earlier routine sla_MAPP, which has
*      been withdrawn.
*
*  Called:
*     sla_EPJ         MDJ to Julian epoch
*     sla_EVP         earth position & velocity
*     sla_DVN         normalise vector
*     sla_PRENUT      precession/nutation matrix
*
*  P.T.Wallace   Starlink   11 April 1990
*/
sla_mappa (eq, date, amprms)
double eq,date,amprms[21];
{

/*  Light time for 1 AU (sec) */
      double cr = 499.004782;

/*  Gravitational radius of the sun x 2 (2*mu/c**2, AU) */
      double gr2 = 2.0*9.87063e-9;
      int i;

      double ebd[3],ehd[3],eh[3],e,vn[3],vm;

/*  Time interval for proper motion correction */
      amprms[0] = sla_epj(date)-eq;

/*  Get Earth barycentric and heliocentric position and velocity */
      sla_evp(date,eq,ebd,&amprms[1],ehd,eh);

/*  Heliocentric direction of earth (normalised) and modulus */
      sla_dvn(eh,&amprms[4],&e);

/*  Light deflection parameter */
      amprms[7] = gr2/e;

/*  Aberration parameters */
      for(i=0; i<3; i++ )
         amprms[i+8] = ebd[i]*cr;

      sla_dvn(&amprms[8],vn,&vm);
      amprms[11] = sqrt(1.0-vm*vm);

/*  Precession/nutation matrix */
      sla_prenut(eq,date,(void *)&amprms[12]);
}

/*
*     - - - -
*      E V P 
*     - - - -
*
*  Barycentric and heliocentric velocity and position of the Earth
*
*  All arguments are double precision
*
*  Given:
*
*     DATE          TDB (loosely ET) as a Modified Julian Date
*                                         (JD-2400000.5)
*
*     DEQX          Julian Epoch (e.g. 2000.0D0) of mean equator and
*                   equinox of the vectors returned.  If DEQX .LE. 0D0,
*                   all vectors are referred to the mean equator and
*                   equinox (FK5) of date DATE.
*
*  Returned (all 3D Cartesian vectors):
*
*     DVB,DPB       barycentric velocity, position
*
*     DVH,DPH       heliocentric velocity, position
*
*  (Units are au/s for velocity and au for position)
*
*  Called:  sla_EPJ, sla_PREC
*
*  Accuracy:
*
*     The maximum deviations from the JPL DE96 ephemeris are as
*     follows:
*
*     barycentric velocity                  42  cm/s
*     barycentric position           0.000 046  au
*
*     heliocentric velocity                 42  cm/s
*     heliocentric position          0.000 011  au
*
*  This routine is adapted from the BARVEL and BARCOR
*  subroutines of P.Stumpff, which are described in
*  Astron. Astrophys. Suppl. Ser. 41, 1-8 (1980).  Most of the
*  changes are merely cosmetic and do not affect the results at
*  all.  However, some adjustments have been made so as to give
*  results that refer to the new (IAU 1976 'FK5') equinox
*  and precession, although the differences these changes make
*  relative to the results from Stumpff's original 'FK4' version
*  are smaller than the inherent accuracy of the algorithm.  One
*  minor shortcoming in the original routines that has NOT been
*  corrected is that better numerical accuracy could be achieved
*  if the various polynomial evaluations were nested.  Note also
*  that one of Stumpff's precession constants differs by 0.001 arcsec
*  from the value given in the Explanatory Supplement to the A.E.
*
*  P.T.Wallace   Starlink   27 November 1991
*/

sla_evp (date, deqx, dvb, dpb, dvh, dph)
double date,deqx,dvb[3],dpb[3],dvh[3],dph[3];
{
      int ideq,i,j,k;
      double t,tsq,a,pertl,
           pertld,pertr,pertrd,cosa,sina,esq,param,twoe,twog,
           phi,f,sinf,cosf,phid,psid,pertp,pertpd,tl,sinlm,coslm,
           sigma,b,plon,pomg,pecc,flatm,flat;

      double dt,dtsq,dlocal,dml,
                       deps,dparam,dpsi,d1pdro,drd,drld,dtl,dsinls,
                       dcosls,dxhd,dyhd,dzhd,dxbd,dybd,dzbd,dcosep,
                       dsinep,dyahd,dzahd,dyabd,dzabd,dr,
                       dxh,dyh,dzh,dxb,dyb,dzb,dyah,dzah,dyab,
                       dzab,depj,deqcor;

      double sn[4], forbel[7],sorbel[17],sinlp[4],coslp[4];
      /* EQUIVALENCE (SORBEL(1),E),(FORBEL(1),G) */

      double  dprema[3][3],w,vw[3];

      static double  dc2pi=6.2831853071796,cc2pi=6.283185;
      static double  ds2r=0.7272205216643e-04;

/* */
/*   Constants DCFEL(I,K) of fast changing elements */
/*                     I=1                I=2              I=3 */
     static double dcfel[8][3] = {
                  { 1.7400353e+00, 6.2833195099091e+02, 5.2796e-06 },
                  { 6.2565836e+00, 6.2830194572674e+02,-2.6180e-06 },
                  { 4.7199666e+00, 8.3997091449254e+03,-1.9780e-05 },
                  { 1.9636505e-01, 8.4334662911720e+03,-5.6044e-05 },
                  { 4.1547339e+00, 5.2993466764997e+01, 5.8845e-06 },
                  { 4.6524223e+00, 2.1354275911213e+01, 5.6797e-06 },
                  { 4.2620486e+00, 7.5025342197656e+00, 5.5317e-06 },
                  { 1.4740694e+00, 3.8377331909193e+00, 5.6093e-06 }
               };

/*
*   Constants DCEPS and CCSEL(I,K) of slowly changing elements
*                      I=1           I=2           I=3
*/
      static double dceps[3] = {  4.093198e-01,-2.271110e-04,-2.860401e-08 };
      static double ccsel[17][3] = {
                   { 1.675104E-02,-4.179579E-05,-1.260516E-07 },
                   { 2.220221E-01, 2.809917E-02, 1.852532E-05 },
                   { 1.589963E+00, 3.418075E-02, 1.430200E-05 },
                   { 2.994089E+00, 2.590824E-02, 4.155840E-06 },
                   { 8.155457E-01, 2.486352E-02, 6.836840E-06 },
                   { 1.735614E+00, 1.763719E-02, 6.370440E-06 },
                   { 1.968564E+00, 1.524020E-02,-2.517152E-06 },
                   { 1.282417E+00, 8.703393E-03, 2.289292E-05 },
                   { 2.280820E+00, 1.918010E-02, 4.484520E-06 },
                   { 4.833473E-02, 1.641773E-04,-4.654200E-07 },
                   { 5.589232E-02,-3.455092E-04,-7.388560E-07 },
                   { 4.634443E-02,-2.658234E-05, 7.757000E-08 },
                   { 8.997041E-03, 6.329728E-06,-1.939256E-09 },
                   { 2.284178E-02,-9.941590E-05, 6.787400E-08 },
                   { 4.350267E-02,-6.839749E-05,-2.714956E-07 },
                   { 1.348204E-02, 1.091504E-05, 6.903760E-07 },
                   { 3.106570E-02,-1.665665E-04,-1.590188E-07 }
                };
/* */
/*   Constants of the arguments of the short-period perturbations */
/*   by the planets:   DCARGS(I,K) */
/*                       I=1               I=2 */
      static double dcargs[15][2] = {
                   { 5.0974222e+00,-7.8604195454652e+02 },
                   { 3.9584962e+00,-5.7533848094674e+02 },
                   { 1.6338070e+00,-1.1506769618935e+03 },
                   { 2.5487111e+00,-3.9302097727326e+02 },
                   { 4.9255514e+00,-5.8849265665348e+02 },
                   { 1.3363463e+00,-5.5076098609303e+02 },
                   { 1.6072053e+00,-5.2237501616674e+02 },
                   { 1.3629480e+00,-1.1790629318198e+03 },
                   { 5.5657014e+00,-1.0977134971135e+03 },
                   { 5.0708205e+00,-1.5774000881978e+02 },
                   { 3.9318944e+00, 5.2963464780000e+01 },
                   { 4.8989497e+00, 3.9809289073258e+01 },
                   { 1.3097446e+00, 7.7540959633708e+01 },
                   { 3.5147141e+00, 7.9618578146517e+01 },
                   { 3.5413158e+00,-5.4868336758022e+02 }
               };

/* */
/*   Amplitudes CCAMPS(N,K) of the short-period perturbations */
/*           N=1          N=2          N=3          N=4          N=5 */
      static double ccamps[15][5] = {
       { -2.279594E-5, 1.407414E-5, 8.273188E-6, 1.340565E-5,-2.490817E-7 },
       { -3.494537E-5, 2.860401E-7, 1.289448E-7, 1.627237E-5,-1.823138E-7 },
       {  6.593466E-7, 1.322572E-5, 9.258695E-6,-4.674248E-7,-3.646275E-7 },
       {  1.140767E-5,-2.049792E-5,-4.747930E-6,-2.638763E-6,-1.245408E-7 },
       {  9.516893E-6,-2.748894E-6,-1.319381E-6,-4.549908E-6,-1.864821E-7 },
       {  7.310990E-6,-1.924710E-6,-8.772849E-7,-3.334143E-6,-1.745256E-7 },
       { -2.603449E-6, 7.359472E-6, 3.168357E-6, 1.119056E-6,-1.655307E-7 },
       { -3.228859E-6, 1.308997E-7, 1.013137E-7, 2.403899E-6,-3.736225E-7 },
       {  3.442177E-7, 2.671323E-6, 1.832858E-6,-2.394688E-7,-3.478444E-7 },
       {  8.702406E-6,-8.421214E-6,-1.372341E-6,-1.455234E-6,-4.998479E-8 },
       { -1.488378E-6,-1.251789E-5, 5.226868E-7,-2.049301E-7, 0.0E0 },
       { -8.043059E-6,-2.991300E-6, 1.473654E-7,-3.154542E-7, 0.0E0 },
       {  3.699128E-6,-3.316126E-6, 2.901257E-7, 3.407826E-7, 0.0E0 },
       {  2.550120E-6,-1.241123E-6, 9.901116E-8, 2.210482E-7, 0.0E0 },
       { -6.351059E-7, 2.341650E-6, 1.061492E-6, 2.878231E-7, 0.0E0 }
       };

/* */
/*   Constants of the secular perturbations in longitude */
/*   CCSEC3 and CCSEC(N,K) */
/*                      N=1           N=2           N=3 */
     static double ccsec3 = -7.757020E-08;
     static double ccsec[4][3]  = {
                  { 1.289600E-06, 5.550147E-01, 2.076942E+00 },
                  { 3.102810E-05, 4.035027E+00, 3.525565E-01 },
                  { 9.124190E-06, 9.990265E-01, 2.622706E+00 },
                  { 9.793240E-07, 5.508259E+00, 1.559103E+01 }
               };

/*   Sidereal rate DCSLD in longitude, rate CCSGD in mean anomaly */
      static double dcsld = 1.990987e-07;
      static double ccsgd = 1.990969E-07;

/*   Some constants used in the calculation of the lunar contribution */
      static double cckm = 3.122140E-05;
      static double ccmld = 2.661699E-06;
      static double ccfdi = 2.399485E-07;

/* */
/*   Constants DCARGM(I,K) of the arguments of the perturbations */
/*   of the motion of the Moon */
/*                       I=1               I=2 */
      static double dcargm[3][2] =  {
                    { 5.1679830e+00, 8.3286911095275e+03 },
                    { 5.4913150e+00,-7.2140632838100e+03 },
                    { 5.9598530e+00, 1.5542754389685e+04 }
              };

/* */
/*   Amplitudes CCAMPM(N,K) of the perturbations of the Moon */
/*            N=1          N=2           N=3           N=4 */
      static double ccampm[3][4] = {
         {  1.097594E-01, 2.896773E-07, 5.450474E-02, 1.438491E-07 },
         { -2.223581E-02, 5.083103E-08, 1.002548E-02,-2.291823E-08 },
         {  1.148966E-02, 5.658888E-08, 8.249439E-03, 4.063015E-08 }
      };

/* */
/*   CCPAMV(K)=A*M*DL/DT (planets), DC1MME=1-MASS(Earth+Moon) */
      static double ccpamv[4] = { 8.326827E-11,1.843484E-11,1.988712E-12,1.881276E-12};
      static double dc1mme = 0.99999696;

/*   CCPAM(K)=A*M(planets), CCIM=INCLINATION(Moon) */
      static double ccpam[4] = {4.960906E-3,2.727436E-3,8.392311E-4,1.556861E-3};
      static double ccim = 8.978749E-2;

/* */
/*   EXECUTION */
/*   --------- */

/*   Control parameter IDEQ, and time arguments */
      ideq = 0;
      if (deqx > 0.0) ideq=1;
      dt = (date-15019.5)/36525.0;
      t = dt;
      dtsq = dt*dt;
      tsq = dtsq;

/*   Values of all elements for the instant DATE */
      for( k=0; k<8; k++) {
         dlocal = fmod(dcfel[k][0]+dt*dcfel[k][1]+dtsq*dcfel[k][2], dc2pi);
         if (k == 0)
            dml = dlocal;
         else
            forbel[k-1] = dlocal;
      }
      deps = fmod(dceps[0]+dt*dceps[1]+dtsq*dceps[2], dc2pi);
      for( k=0; k<17; k++) {
         sorbel[k] = fmod(ccsel[k][0]+t*ccsel[k][1]+tsq*ccsel[k][2], cc2pi);
      }

/*   Secular perturbations in longitude */
      for( k=0; k<4; k++) {
         a = fmod(ccsec[k][1]+t*ccsec[k][2], cc2pi);
         sn[k] = sin(a);
      }

/*   Periodic perturbations of the EMB (Earth-Moon barycentre) */
      pertl =  ccsec[0][0]          *sn[0] +ccsec[1][0]*sn[1]+
              (ccsec[2][0]+t*ccsec3)*sn[2] +ccsec[3][0]*sn[3];
      pertld = 0.0;
      pertr = 0.0;
      pertrd = 0.0;
      for( k=0; k<15; k++) {
         a = fmod(dcargs[k][0]+dt*dcargs[k][1], dc2pi);
         cosa = cos(a);
         sina = sin(a);
         pertl = pertl + ccamps[k][0]*cosa+ccamps[k][1]*sina;
         pertr = pertr + ccamps[k][2]*cosa+ccamps[k][3]*sina;
         if(k<10) {
            pertld = pertld+
                     (ccamps[k][1]*cosa-ccamps[k][0]*sina)*ccamps[k][4];
            pertrd = pertrd+
                     (ccamps[k][3]*cosa-ccamps[k][2]*sina)*ccamps[k][4];
         }
      }

/*   Elliptic part of the motion of the EMB */
      /* EQUIVALENCE (SORBEL(1),E),(FORBEL(1),G) */
      esq = sorbel[0]*sorbel[0];
      dparam = 1.0-esq;
      param = dparam;
      twoe = sorbel[0]+sorbel[0];
      twog = forbel[0]+forbel[0];
      phi = twoe*((1.0-esq*0.125)*sin(forbel[0])+sorbel[0]*0.625*sin(twog)
                +esq*0.5416667*sin(forbel[0]+twog) );
      f = forbel[0]+phi;
      sinf = sin(f);
      cosf = cos(f);
      dpsi = dparam/(1.0+sorbel[0]*cosf);
      phid = twoe*ccsgd*((1.0+esq*1.5)*cosf+sorbel[0]*(1.25-sinf*sinf*0.5));
      psid = ccsgd*sorbel[0]*sinf/sqrt(param);

/*   Perturbed heliocentric motion of the EMB */
      d1pdro = 1.0+pertr;
      drd = d1pdro*(psid+dpsi*pertrd);
      drld = d1pdro*dpsi*(dcsld+phid+pertld);
      dtl = fmod(dml+phi+pertl, dc2pi);
      dsinls = sin(dtl);
      dcosls = cos(dtl);
      dxhd = drd*dcosls-drld*dsinls;
      dyhd = drd*dsinls+drld*dcosls;

/*   Influence of eccentricity, evection and variation on the */
/*   geocentric motion of the Moon */
      pertl = 0.0;
      pertld = 0.0;
      pertp = 0.0;
      pertpd = 0.0;
      for( k=0; k<3; k++) {
         a = fmod(dcargm[k][0]+dt*dcargm[k][1], dc2pi);
         sina = sin(a);
         cosa = cos(a);
         pertl = pertl +ccampm[k][0]*sina;
         pertld = pertld+ccampm[k][1]*cosa;
         pertp = pertp +ccampm[k][2]*cosa;
         pertpd = pertpd-ccampm[k][3]*sina;
      }

/*   Heliocentric motion of the Earth */
      tl = forbel[1]+pertl;
      sinlm = sin(tl);
      coslm = cos(tl);
      sigma = cckm/(1.0+pertp);
      a = sigma*(ccmld+pertld);
      b = sigma*pertpd;
      dxhd = dxhd+(a*sinlm)+(b*coslm);
      dyhd = dyhd-(a*coslm)+(b*sinlm);
      dzhd =     -(sigma*ccfdi*cos(forbel[2]));

/*   Barycentric motion of the Earth */
      dxbd = dxhd*dc1mme;
      dybd = dyhd*dc1mme;
      dzbd = dzhd*dc1mme;
      for( k=0; k<4; k++) {
         plon = forbel[k+3];
         pomg = sorbel[k+1];
         pecc = sorbel[k+9];
         tl = fmod(plon+2.0*pecc*sin(plon-pomg), cc2pi);
         sinlp[k] = sin(tl);
         coslp[k] = cos(tl);
         dxbd = dxbd+(ccpamv[k]*(sinlp[k]+pecc*sin(pomg)));
         dybd = dybd-(ccpamv[k]*(coslp[k]+pecc*cos(pomg)));
         dzbd = dzbd-(ccpamv[k]*sorbel[k+13]*cos(plon-sorbel[k+5]));
      }

/*   Transition to mean equator of date */
      dcosep = cos(deps);
      dsinep = sin(deps);
      dyahd = dcosep*dyhd-dsinep*dzhd;
      dzahd = dsinep*dyhd+dcosep*dzhd;
      dyabd = dcosep*dybd-dsinep*dzbd;
      dzabd = dsinep*dybd+dcosep*dzbd;

/*   Heliocentric coordinates of the Earth */
      dr = dpsi*d1pdro;
      flatm = ccim*sin(forbel[2]);
      a = sigma*cos(flatm);
      dxh = dr*dcosls-(a*coslm);
      dyh = dr*dsinls-(a*sinlm);
      dzh =          -(sigma*sin(flatm));

/*   Barycentric coordinates of the Earth */
      dxb = dxh*dc1mme;
      dyb = dyh*dc1mme;
      dzb = dzh*dc1mme;
      for( k=0; k<4; k++) {
         flat = sorbel[k+13]*sin(forbel[k+3]-sorbel[k+5]);
         a = ccpam[k]*(1.0-sorbel[k+9]*cos(forbel[k+3]-sorbel[k+1]));
         b = a*cos(flat);
         dxb = dxb-(b*coslp[k]);
         dyb = dyb-(b*sinlp[k]);
         dzb = dzb-(a*sin(flat));
      }

/*   Transition to mean equator of date */
      dyah = dcosep*dyh-dsinep*dzh;
      dzah = dsinep*dyh+dcosep*dzh;
      dyab = dcosep*dyb-dsinep*dzb;
      dzab = dsinep*dyb+dcosep*dzb;

/*   Copy result components into vectors, correcting for FK4 equinox */
      depj=sla_epj(date);
      deqcor = ds2r*(0.035+0.00085*(depj-1950.0));
      dvh[0] = dxhd-deqcor*dyahd;
      dvh[1] = dyahd+deqcor*dxhd;
      dvh[2] = dzahd;
      dvb[0] = dxbd-deqcor*dyabd;
      dvb[1] = dyabd+deqcor*dxbd;
      dvb[2] = dzabd;
      dph[0] = dxh-deqcor*dyah;
      dph[1] = dyah+deqcor*dxh;
      dph[2] = dzah;
      dpb[0] = dxb-deqcor*dyab;
      dpb[1] = dyab+deqcor*dxb;
      dpb[2] = dzab;

/*   Was precession to another equinox requested? */
      if (ideq!=0) {

/*     Yes: compute precession matrix from MJD DATE to Julian epoch DEQX */
         sla_prec(depj,deqx,dprema);

/*     Rotate DVH */
         for(j=0; j<3; j++ ) {
            w=0.0;
            for(i=0; i<3; i++ ) {
               w=w+dprema[i][j]*dvh[i];
            }
            vw[j]=w;
         }
         for(j=0; j<3; j++ ) {
            dvh[j]=vw[j];
         }

/*     Rotate DVB */
         for(j=0; j<3; j++ ) {
            w=0.0;
            for(i=0; i<3; i++ ) {
               w=w+dprema[i][j]*dvb[i];
            }
            vw[j]=w;
         }
         for(j=0; j<3; j++ ) {
            dvb[j]=vw[j];
         }

/*     Rotate DPH */
         for(j=0; j<3; j++ ) {
            w=0.0;
            for(i=0; i<3; i++ ) {
               w=w+dprema[i][j]*dph[i];
            }
            vw[j]=w;
         }
         for(j=0; j<3; j++ ) {
            dph[j]=vw[j];
         }

/*     Rotate DPB */
         for(j=0; j<3; j++ ) {
            w=0.0;
            for(i=0; i<3; i++ ) {
               w=w+dprema[i][j]*dpb[i];
            }
            vw[j]=w;
         }
         for(j=0; j<3; j++ ) {
            dpb[j]=vw[j];
         }
      }
}

/*
*     - - - - 
*      D V N
*     - - - -
*
*  Normalises a 3-vector also giving the modulus (double precision)
*
*  Given:
*     V       dp(3)      vector
*
*  Returned:
*     UV      dp(3)      unit vector in direction of V
*     VM      dp         modulus of V
*
*  If the modulus of V is zero, UV is set to zero as well
*
*  P.T.Wallace   Starlink   November 1984
*/
sla_dvn (v, uv, vm)
double v[3],uv[3], *vm;
{
  int i;
  double w1,w2;

/*  Modulus */
      w1=0.0;
      for(i=0; i<3; i++ ) {
         w2=v[i];
         w1=w1+w2*w2;
      }
      w1=sqrt(w1);
      *vm=w1;

/*  Normalise the vector */
      if(w1 <= 0.0)
         w1=1.0;
      for(i=0; i<3; i++ )
         uv[i]=v[i]/w1;
}

/*
*     - - - - - -
*      M A P Q K
*     - - - - - -
*
*  Quick mean to apparent place:  transform a star RA,Dec from
*  mean place to geocentric apparent place, given the
*  star-independent parameters.
*
*  Use of this routine is appropriate when efficiency is important
*  and where many star positions, all referred to the same equator
*  and equinox, are to be transformed for one epoch.  The
*  star-independent parameters can be obtained by calling the
*  sla_MAPPA routine.
*
*  If the parallax and proper motions are zero the sla_MAPQKZ
*  routine can be used instead.
*
*  The reference frames and timescales used are post IAU 1976.
*
*  Given:
*     RM,DM    d      mean RA,Dec (rad)
*     PR,PD    d      proper motions:  RA,Dec changes per Julian year
*     PX       d      parallax (arcsec)
*     RV       d      radial velocity (km/sec, +ve if receding)
*
*     AMPRMS   d(21)  star-independent mean-to-apparent parameters:
*             
*       (1)      time interval for proper motion (Julian years)
*       (2-4)    barycentric position of the Earth (AU)
*       (5-7)    heliocentric direction of the Earth (unit vector)
*       (8)      (grav rad Sun)*2/(Sun-Earth distance)
*       (9-11)   barycentric Earth velocity in units of c
*       (12)     sqrt(1-v**2) where v=modulus(ABV)
*       (13-21)  precession/nutation (3,3) matrix
*
*  Returned:
*     RA,DA    d      apparent RA,Dec (rad)
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Notes:
*
*  1)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to
*      the mean equinox and equator of epoch EQ.
*
*  2)  Within about 300 arcsec of the centre of the Sun the
*      gravitational deflection term is set to zero to avoid
*      overflow.  Otherwise no account is taken of the
*      impossibility of observing stars which lie behind
*      the Sun.
*
*  3)  This routine replaces earlier routine sla_MAPQ, which is
*      to be withdrawn.
*
*  Called:
*     sla_DCS2C       spherical to Cartesian
*     sla_DVDV        dot product
*     sla_DMXV        matrix x vector
*     sla_DCC2S       Cartesian to spherical
*     sla_DRANRM      normalise angle 0-2Pi
*
*  P.T.Wallace   Starlink   November 1988
*/
sla_mapqk (rm, dm, pr, pd, px, rv, amprms, ra, da)
double rm,dm,pr,pd,px,rv,amprms[21],*ra,*da;
{

/*  Arc seconds to radians */
      static double as2r=0.4848136811095359949e-5;

/*  Km/s to AU/year */
      double vf = 0.21094502;

      int i;

      double pmt,gr2e,ab1,eb[3],ehn[3],abv[3],
             q[3],pxr,w,em[3],p[3],pn[3],pde,pdep1,
             p1[3],p1dv,p1dvp1,p2[3],p3[3];

      double sla_dvdv(),sla_dranrm();

/*  Unpack scalar and vector parameters */
      pmt = amprms[0];
      gr2e = amprms[7];
      ab1 = amprms[11];
      for(i=0; i<3; i++ ) {
         eb[i] = amprms[i+1];
         ehn[i] = amprms[i+4];
         abv[i] = amprms[i+8];
      }

/*  Spherical to x,y,z */
      sla_dcs2c(rm,dm,q);

/*  Space motion (radians per year) */
      pxr = px*as2r;
      w = vf*rv*pxr;
      em[0] = -pr*q[1]-pd*cos(rm)*sin(dm)+w*q[0];
      em[1] =  pr*q[0]-pd*sin(rm)*sin(dm)+w*q[1];
      em[2] =          pd*cos(dm)        +w*q[2];

/*  Geocentric direction of star (normalised) */
      for(i=0; i<3; i++ )
         p[i] = q[i]+pmt*em[i]-pxr*eb[i];
      sla_dvn(p,pn,&w);

/*  Light deflection */
      pde = sla_dvdv(pn,ehn);
      pdep1 = 1.0+pde;
      w = gr2e/max(pdep1,1e-6);
      for(i=0; i<3; i++ )
         p1[i] = pn[i]+w*(ehn[i]-pde*pn[i]);

/*  Aberration */
      p1dv = sla_dvdv(p1,abv);
      p1dvp1 = p1dv+1.0;
      w = 1.0+p1dv/(ab1+1.0);
      for(i=0; i<3; i++ )
         p2[i] = (ab1*p1[i]+w*abv[i])/p1dvp1;

/*  Precession and nutation */
      sla_dmxv((void *)&amprms[12],p2,p3);

/*  Geocentric apparent RA,Dec */
      sla_dcc2s(p3,ra,da);
      *ra = sla_dranrm(*ra);
}

/*
*     - - - - -
*      E P C O
*     - - - - -
*
*  Convert an epoch into the appropriate form - 'B' or 'J'
*
*  Given:
*     K0    char    form of result:  'B'=Besselian, 'J'=Julian
*     K     char    form of given epoch:  'B' or 'J'
*     E     dp      epoch
*
*  Called:  sla_EPB, sla_EPJ2D, sla_EPJ, sla_EPB2D
*
*  Notes:
*
*     1) The result is always either equal to or very close to
*        the given epoch E.  The routine is required only in
*        applications where punctilious treatment of heterogeneous
*        mixtures of star positions is necessary.
*
*     2) K0 and K are not validated.  They are interpreted as follows:
*
*        o  If K0 and K are the same the result is E.
*        o  If K0 is 'B' and K isn't, the conversion is J to B.
*        o  In all other cases, the conversion is B to J.
*
*  P.T.Wallace   Starlink   8 April 1990
*/
double sla_epco (k0, k, e)
char *k0,*k;
double e;
{
  double ret, sla_epb(),sla_epj2d(),sla_epj(),sla_epb2d();

      if (strcmp(k,k0)==0)
         ret = e;
      else if (strcmp(k0, "B")==0 )
         ret=sla_epb(sla_epj2d(e));
      else
         ret=sla_epj(sla_epb2d(e));
      return(ret);
}

/*
*     - - - - - -
*      E P J 2 D
*     - - - - - -
*
*  Conversion of Julian Epoch to Modified Julian Date (double precision)
*
*  Given:
*     EPJ      dp       Julian Epoch
*
*  The result is the Modified Julian Date (JD - 2400000.5).
*
*  Reference:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*
*  P.T.Wallace   Starlink   February 1984
*/
double sla_epj2d (epj)
double epj;
{
      return( 51544.5 + (epj-2000.0)*365.25);
}

/*
*     - - - - - -
*      E P B 2 D
*     - - - - - -
*
*  Conversion of Besselian Epoch to Modified Julian Date
*  (double precision)
*
*  Given:
*     EPB      dp       Besselian Epoch
*
*  The result is the Modified Julian Date (JD - 2400000.5).
*
*  Reference:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*
*  P.T.Wallace   Starlink   February 1984
*/
double sla_epb2d(epb)
double epb;
{
      return( 15019.81352 + (epb-1900.0)*365.242198781 );
}

/*
*     - - - -
*      E P B
*     - - - -
*
*  Conversion of Modified Julian Date to Besselian Epoch
*  (double precision)
*
*  Given:
*     DATE     dp       Modified Julian Date (JD - 2400000.5)
*
*  The result is the Besselian Epoch.
*
*  Reference:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*
*  P.T.Wallace   Starlink   February 1984
*/
double sla_epb(date)
double date;
{
      return( 1900.0 + (date-15019.81352)/365.242198781 );
}

/*
*     - - - - -
*      D V D V
*     - - - - -
* 
*  Scalar product of two 3-vectors  (double precision)
*
*  Given:
*      VA      dp(3)     first vector
*      VB      dp(3)     second vector
*
*  The result is the scalar product VA.VB (double precision)
*
*  P.T.Wallace   Starlink   November 1984
*/
double sla_dvdv(va, vb)
double va[3],vb[3];
{
      return(va[0]*vb[0]+va[1]*vb[1]+va[2]*vb[2]);
}

double current_time()
{
  struct tm *gmt, *gmtime();
  time_t t;
  double dys, year, frac;

  time(&t);
  gmt = gmtime(&t);

  dys = (double)dysize(gmt->tm_year);
  year = 1900.0 + (double)gmt->tm_year;  
  frac = ((double)gmt->tm_yday + ((double)gmt->tm_hour +
         ((double)gmt->tm_min + (double)gmt->tm_sec/60.0)/60.0)/24.0)/dys;
  
  return( year + frac);
}


/*  gall -> B1950 

h_slaGe50(ret,p)
struct NODE *ret;
struct NODE *p[];
{
  ret->type = typ_int;
  if( p[0]->type != typ_double ||
      p[1]->type != typ_double ||
      p[2]->type != typ_double ||
      p[3]->type != typ_double ) {
    ret->value.ival = 1;
    return;
  } else
    ret->value.ival = 0;

  slaGe50(p[0]->value.fval, p[1]->value.fval, &p[2]->value.fval, &p[3]->value.fval );
}
*/

/*  gall -> J2000 

h_slaGaleq(ret,p)
struct NODE *ret;
struct NODE *p[];
{
  ret->type = typ_int;
  if( p[0]->type != typ_double ||
      p[1]->type != typ_double ||
      p[2]->type != typ_double ||
      p[3]->type != typ_double ) {
    ret->value.ival = 1;
    return;
  } else
    ret->value.ival = 0;

  slaGaleq(p[0]->value.fval, p[1]->value.fval, &p[2]->value.fval, &p[3]->value.fval );
}
*/

/*  J2000 -> gal 

h_slaEqgal(ret,p)
struct NODE *ret;
struct NODE *p[];
{
  ret->type = typ_int;
  if( p[0]->type != typ_double ||
      p[1]->type != typ_double ||
      p[2]->type != typ_double ||
      p[3]->type != typ_double ) {
    ret->value.ival = 1;
    return;
  } else
    ret->value.ival = 0;

  slaEqgal(p[0]->value.fval, p[1]->value.fval, &p[2]->value.fval, &p[3]->value.fval );
}
*/

/*  B1950 ->gall 

h_slaEg50(ret,p)
struct NODE *ret;
struct NODE *p[];
{
  ret->type = typ_int;
  if( p[0]->type != typ_double ||
      p[1]->type != typ_double ||
      p[2]->type != typ_double ||
      p[3]->type != typ_double ) {
    ret->value.ival = 1;
    return;
  } else
    ret->value.ival = 0;

  slaEg50(p[0]->value.fval, p[1]->value.fval, &p[2]->value.fval, &p[3]->value.fval );
}
*/


int isleap(year)
int year;
{
  return((year % 4 == 0 && year % 100 != 0) || year % 400 == 0);
}

int dysize(yr)
{
  return( isleap( yr ) ? 366 : 365 );
}



/*
 Calculates Parallactic angle (degrees) from
  TLST (Sidereal time, Hours)
  RA   (Hours)
  DEC  (Degrees)
  ALAT (Latitude of observatory, degrees)

    Written by D. T. Emerson.
*/

double para(tlst, ra, dec, alat )
double  tlst, ra, dec, alat;
{
  double  ha, har, rar, decr, alatr, tanp1, tanp2, angle;
  double sin(), cos(), atan2();
      
/* Calculate Hour Angle HA (hours)  */

  ha = tlst - ra;

/* Convert to consistent units (here radians)  */  

  har   = ha*PI/12.0;
  rar   = ra*PI/12.0;
  decr  = dec*PI/180.0;
  alatr = alat*PI/180.0;

/* Calculate TANP, i.e. TAN(Parallactic angle) using ATAN2  */

  tanp1 = cos(alatr)*sin(har);
  tanp2 = (sin(alatr)*cos(decr) - cos(alatr)*sin(decr)*cos(har));

/* Catch indeterminate case, where numerator and denominator = 0  */

  if(tanp1 == 0.0 && tanp2 == 0.0)
    angle = 0.0;

/* Use ATAN2 instead of ATAN, to allow for denominator possibly being zero,
   and to get signs of all 4 quadrants correct.  */

   else
     angle = atan2(tanp1,tanp2)*180.0/PI;

  return(angle);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

#define CTD (360.0/(PI*2.0))
#define CTR (PI*2.0/360.0)
#define AOLAT   18.353806         /* degrees north */

double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

/* returns the azimuth and za for a source at apparent ra, dec and lst
   ra (hours)
   dec (degrees)
   lst (hours)
*/

int compute_azza( ra, dec, lst, paz, pza  ) 
double ra, dec, lst, *paz, *pza;
{
  double ha, az, el, sindec, cosdec, sinlat, coslat,
                cosha, sinel, cosel, cosazel, cosaz;


/* Calculate hour angle and convert to degrees (use convention
   that hour angle ranges between +/- 12 hours)    */

  ha = lst - ra;
  if( ha > 12.0 )
    ha = ha - 24.0;
  else if( ha < -12.0 )
    ha = ha + 24.0;
  ha = ha * 15.0;

/* Compute trig quantities needed; make angle calculations double
   precision to cut down on round-off-induced overflows near the
   zenith and the meridian.    */

  sindec = sin(dec*CTR);
  cosdec = cos(dec*CTR);
  sinlat = sin(AOLAT*CTR);
  coslat = cos(AOLAT*CTR);
  cosha  = cos(ha*CTR);

/* Compute elevation (trap for blow-ups near the zenith):  */

  sinel = sindec * sinlat + cosdec * cosha * coslat;

  if (sinel >= -1.0 && sinel <= 1.0)
    el = asin(sinel)*CTD;
  else
    printf("elevation too high");

/* Compute azimuth (trap for blow-ups near HA = 0):    */

  cosel = cos(el*CTR);
  cosazel = sindec * coslat - cosdec * cosha * sinlat;
  cosaz = (1.0/cosel) * cosazel;
  if (cosaz >= -1.0 && cosaz <= 1.0)
    az = acos(cosaz)*CTD;
  else {
    if (dec <= AOLAT)
      az = 180.0;
    else
      az = 0.0;
  }

/* Resolve quadrant ambiguity of azimuth:    */

  if (ha > 0.0)
    az = 360.0 - az;

  *pza = 90.0 - el;

  az += 180.0; /* arecibo az is 180 away from azel dish */
  if ( az > 360.0 )
    az -= 360.0;

  *paz = az;
  /* printf("az %f dec %f ha %f\n", az, dec, ha ); */
  return 0;
}

/* 
  thanks to Avinash Deshspande, he computed these mixed up 
  equations
*/

compute_zaha( az, dec, pza, pha )
double az;
double dec;
double *pza;
double *pha;
{
  double sinaz, cosaz, sindec, cosdec, sinlat, coslat;
  double sinelp, sinelm, sinel, cosel, sinha;
  double el, ha, za, sq, zam, zap;

  sinaz  = sin(az*CTR);
  cosaz  = cos(az*CTR);
  sindec = sin(dec*CTR);
  cosdec = cos(dec*CTR);
  sinlat = sin(AOLAT*CTR);
  coslat = cos(AOLAT*CTR);

  sq = cosdec*cosdec - coslat*coslat*sinaz*sinaz;
  if( sq < 0.0 )
    return;
  sinelp = (sinlat*sindec + coslat * cosaz * 
    sqrt( sq ) ) /
    ( 1.0 - coslat*coslat * sinaz*sinaz );

  sinelm = (sinlat*sindec - coslat * cosaz * 
    sqrt( sq ) ) /
    ( 1.0 - coslat*coslat * sinaz*sinaz );

  zam = 90.0 - deg_asin( sinelm );
  zap = 90.0 - deg_asin( sinelp );

/* find a za in the right range */

  if( zam > 0.0 && zam < 20.0 ) {
    sinel = sinelm;
    za = zam;
  } else {
    sinel = sinelp;
    za = zap;
  }
    
  cosel = (sindec - sinel * sinlat ) / coslat / cosaz ;

  sinha =  -1.0 * cosel * sinaz / cosdec;
  ha = deg_asin( sinha );
  ha /= 15;

  *pza = za;
  *pha = ha;
}

/* compute radec offsets from az/za offsets and positions */

int off_hadec( daz, dza, az, za, ha, dec, pdha, pddec )
double daz;
double dza; 
double az; 
double za; 
double ha; 
double dec; 
double *pdha;
double *pddec;
{
  double cosza, coslat, cosaz, cosdec;
  double sinza, sinlat, sinaz;
  double cosha, sinha, ddec, hacd;

   daz = CTR*daz;
   dza = CTR*dza;
   az = CTR*az;
   za = CTR*za;
   ha =  ha*PI/12.0;
   dec = CTR*dec;
  
   cosza = cos(za);
   coslat = cos(AOLAT);
   cosaz = cos(az);
   sinza = sin(za);
   sinlat = sin(AOLAT);
   sinaz = sin(az);
   cosdec = cos(dec);
   cosha = cos(ha);
   cosdec = cos(dec);

   ddec = (dza * ( cosza *coslat * cosaz - sinza * sinlat ) - 
          daz * ( coslat * sinaz ) ) / cosdec; 

   *pddec = CTD*ddec;

   hacd = -1.0 * (cosha * cosaz + sinha * sinaz *sinlat ) * daz -
                 (cosha * sinaz * cosza - 
                   sinha * ( sinza * coslat + cosza * sinlat * cosaz )) * dza;

   *pdha = hacd / cosdec;
}


#include <stdlib.h>
#include <math.h>
/*
   compute apparent ra, dec from az, el, lst
   Zombeck 1982 p. 71

   az, el, ra, dec in degrees
   lst in hours

*/

double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

#define AOLAT   18.353806         /* degrees north */

/* az za in degrees lst in hours */

compute_radec( az, za, lst, pra, pdec )
double az, za, lst;
double *pra, *pdec;
{
  double sindec, sinaz, sinel, cosaz, cosel;
  double cosdec, dec, ha, cosha, sinlat, coslat;
  double el;

  az += 180.0;
  if( az > 360.0 )
   az -= 360.0;
  el = 90.0 - za;
  sinel = deg_sin(el);
  cosel = deg_cos(el);
  sinaz = deg_sin(az);
  cosaz = deg_cos(az);
  sinlat = deg_sin(AOLAT);
  coslat = deg_cos(AOLAT);

  sindec = sinel * sinlat + cosel * coslat * cosaz;

  if( sindec >= -1.0 && sindec <= 1.0 )
    dec = deg_asin( sindec );
  else
    dec = 90.0;

  cosdec = deg_cos( dec );
  cosha = ( sinel * coslat - cosel * cosaz * sinlat)/cosdec ;

  if( cosha >= -1.0 && cosha <= 1.0 )
    ha = deg_acos( cosha )/15.0;
  else
    ha = 0.0;

/* Resolve quadrant ambiguity */

  if( az < 180.0 )
    ha = -ha;

  *pra =  fmod(lst - ha, 24.0);
  if( *pra < 0.0 )
    *pra += 24.0;

  *pdec = dec;
}

struct ALFAOFF {
  double az;
  double za;
};

static struct ALFAOFF *alfaoff = NULL;

double roty( double, double, double);
double rotx( double, double, double);
int find_alfarot(double, int, double *, double * );
double find_alfaoff(int, int);
double  deg_sin(double);
double  deg_cos(double);
double  deg_tan(double);
double  deg_asin(double);
double  deg_acos(double);
double  deg_atan(double);
double  sin(), cos(), tan();

double find_alfaoff(beam, isza)
int beam;
int isza;
{
  int count;
  char *s;
  char *delim = " \t\n";
  FILE *fd;
  char line[80];

  if( !alfaoff ) {
    alfaoff = (struct ALFAOFF *)malloc( sizeof(struct ALFAOFF)*7);
    bzero( alfaoff, sizeof(struct ALFAOFF)*7);

    alfaoff[0].az=0.0;
    alfaoff[0].za=0.0;
    alfaoff[1].az=-164.529999999/3600.0;
    alfaoff[1].za=332.558085181/3600.0;
    alfaoff[2].az=-329.06/3600.0;
    alfaoff[2].za=0.0;
    alfaoff[3].az=-164.529999999/3600.0; 
    alfaoff[3].za=-332.558085181/3600.0;
    alfaoff[4].az=164.529999999/3600.0;
    alfaoff[4].za=-332.558085181/3600.0;
    alfaoff[5].az=329.06/3600.0;
    alfaoff[5].za=0.0;
    alfaoff[6].az=164.529999999/3600.0;
    alfaoff[6].za=332.558085181/3600.0;

  }

  if ( beam < 0 || beam > 6 )
    return 0.0;
  else if( isza )
    return alfaoff[beam].za;
  else
    return alfaoff[beam].az;
}

find_alfarot(angle, beam, paz, pza )
double angle;
int beam;
double *paz;
double *pza;
{
  double az, za;

  az  = find_alfaoff(beam, 0);
  za  = find_alfaoff(beam, 1);

 *paz = rotx( az, za, angle );
 *pza = roty( az, za, angle );
}

#define ALFARATIO (329.06/384.005)

/* negative angle associated with clockwize rotation */

double rotx( offx, offy, ang )
double offx, offy, ang;
{
  double ret, h, a1;

  offy *= ALFARATIO;

  h = sqrt( offx*offx + offy*offy );
  a1 = atan2( offy, offx);
  ret = h * cos( a1 - ang/180.0*PI );


  return ret;
}

/* adjust for the ellipse */
/* negative angle associated with clockwize rotation */

double roty( offx, offy, ang )
double offx, offy, ang;
{
  double ret, h, a1;

  offy *= ALFARATIO;
  h = sqrt( offx*offx + offy*offy );
  a1 = atan2( offy, offx);
  ret = h * sin( a1 - ang/180.0*PI );
  ret /= ALFARATIO;
  return ret;
}

/* given 
  J2000 ra/dec in hh.hhh ddd.ddd
  lst time in hours from midnight
  current epoch as in 2004.56789
  rotation angle
  alfa beam

  compute new ra/dec based on alfa beam offsets

*/

alfa_position( ra, dec, lst, epoch, angle, off1, off2, beam, pra, pdec )
double ra;
double dec;
double lst;
double epoch;
double angle;
double off1;
double off2;
int beam;
double *pra;
double *pdec;
{
  double rra, rdec, nra, ndec, bra, bdec, offaz, offza;
  double saz, sza, nsaz, nsza;

  rra = ra*PI*2.0/24.0;
  rdec = dec*PI*2.0/360.0;
  sla_preces("FK5",2000.0,epoch, &rra, &rdec );
  nra = rra * 24.0 / 2.0 / PI;
  ndec = rdec * 180.0 / PI;

  compute_azza( nra, ndec, lst, &saz, &sza );

  find_alfarot(angle, beam, &offaz, &offza );

  saz = saz - (offaz - off1)/deg_sin(sza);
  sza = sza - (offza - off2);
  compute_radec( saz, sza, lst, &bra, &bdec );

  rra = bra*2.0*PI/24.0;
  rdec = bdec*PI/180.0;
  sla_preces("FK5",epoch, 2000.0, &rra, &rdec );
  bra = rra * 24.0 / 2.0 / PI;
  bdec = rdec * 180.0 / PI;
  *pra = bra;
  *pdec = bdec;

  /* printf("%d     %13.9f     %13.9f\n", beam, bra, bdec );*/
}
/* given 
  J2000 ra/dec in hh.hhh ddd.ddd
  lst time in hours from midnight
  current epoch as in 2004.56789
  rotation angle (deg)
  alfa beam (0-6)

  compute new ra/dec based on alfa beam offsets

main(int argc, char **argv){
double ra;
double dec;
double lst;
double epoch;
double angle;
double off1;
double off2;
int beam;
double pra;
double pdec;
ra=atof(argv[1]);
dec=atof(argv[2]);
lst=atof(argv[3]);
epoch=atof(argv[4]);
angle=atof(argv[5]);
off1=0.0;
off2=0.0;
beam=atoi(argv[6]);
printf("%f %f %f %f %f %f %f %d\n",ra, dec, lst, epoch, angle, off1, off2,beam);
alfa_position( ra, dec, lst, epoch, angle, off1, off2, beam, &pra, &pdec );
printf("%f %f\n",pra,pdec);
}
*/
