#include "slalib.h"
#include "slamac.h"
double slaDranrm ( double angle ) /*includefile*/
/*
**  - - - - - - - - - -
**   s l a D r a n r m
**  - - - - - - - - - -
**
**  Normalize angle into range 0-2 pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the range 0-2 pi (double).
**
**  Defined in slamac.h:  D2PI, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double w;

   w = dmod ( angle, D2PI );
   return ( w >= 0.0 ) ? w : w + D2PI;
}
double slaGmst ( double ut1 ) /*includefile*/
/*
**  - - - - - - - -
**   s l a G m s t
**  - - - - - - - -
**
**  Conversion from Universal Time to Sidereal Time.
**
**  (double precision)
**
**  Given:
**    ut1    double     Universal Time (strictly UT1) expressed as
**                      Modified Julian Date (JD-2400000.5)
**
**  The result is the Greenwich Mean Sidereal Time (double
**  precision, radians).
**
**  The IAU 1982 expression (see page S15 of the 1984 Astronomical
**  Almanac) is used, but rearranged to reduce rounding errors.
**  This expression is always described as giving the GMST at
**  0 hours UT.  In fact, it gives the difference between the
**  GMST and the UT, which happens to equal the GMST (modulo
**  24 hours) at 0 hours UT each day.  In this routine, the
**  entire UT is used directly as the argument for the
**  standard formula, and the fractional part of the UT is
**  added separately;  note that the factor 1.0027379... does
**  not appear.
**
**  See also the routine slaGmsta, which delivers better numerical
**  precision by accepting the UT date and time as separate arguments.
**
**  Called:  slaDranrm
**
**  Defined in slamac.h:  D2PI, DS2R, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double tu;

/* Julian centuries from fundamental epoch J2000 to this UT */
   tu = ( ut1 - 51544.5 ) / 36525.0;

/* GMST at this UT */
   return slaDranrm ( dmod ( ut1, 1.0 ) * D2PI +
                       ( 24110.54841 +
                       ( 8640184.812866 +
                       ( 0.093104 - 6.2e-6 * tu ) * tu ) * tu ) * DS2R );
}
void slaCaldj ( int iy, int im, int id, double *djm, int *j ) /*includefile*/
/*
**  - - - - - - - - -
**   s l a C a l d j
**  - - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  (Includes century default feature:  use slaCldj for years
**   before 100AD.)
**
**  Given:
**     iy,im,id   int      year, month, day in Gregorian calendar
**
**  Returned:
**     *djm       double   Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j         int      status:
**                           0 = ok
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  Acceptable years are 00-49, interpreted as 2000-2049,
**                       50-99,     "       "  1950-1999,
**                       100 upwards, interpreted literally.
**
**  Called:  slaCldj
**
**  Last revision:   21 October 1993
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int ny;

/* Default century if appropriate */
   if ( ( iy >= 0 ) && ( iy <= 49 ) )
      ny = iy + 2000;
   else if ( ( iy >= 50 ) && ( iy <= 99 ) )
      ny = iy + 1900;
   else
      ny = iy;

/* Modified Julian Date */
   slaCldj ( ny, im, id, djm, j );
}

void slaCalyd(int iy, int im, int id, int *ny, int *nd, int *j )/*includefile*/
/*
**  - - - - - - - - -
**   s l a C a l y d
**  - - - - - - - - -
**
**  Gregorian calendar date to year and day in year (in a Julian
**  calendar aligned to the 20th/21st century Gregorian calendar).
**
**  (Includes century default feature:  use slaClyd for years
**   before 100AD.)
**
**  Given:
**     iy,im,id   int    year, month, day in Gregorian calendar
**                       (year may optionally omit the century)
**  Returned:
**     *ny        int    year (re-aligned Julian calendar)
**     *nd        int    day in year (1 = January 1st)
**     *j         int    status:
**                         0 = OK
**                         1 = bad year (before -4711)
**                         2 = bad month
**                         3 = bad day (but conversion performed)
**
**  Notes:
**
**  1  This routine exists to support the low-precision routines
**     slaEarth, slaMoon and slaEcor.
**
**  2  Between 1900 March 1 and 2100 February 28 it returns answers
**     which are consistent with the ordinary Gregorian calendar.
**     Outside this range there will be a discrepancy which increases
**     by one day for every non-leap century year.
**
**  3  Years in the range 50-99 are interpreted as 1950-1999, and
**     years in the range 00-49 are interpreted as 2000-2049.
**
**  Called:  slaClyd
**
**  Last revision:   22 September 1995
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i;

/* Default century if appropriate */
   if ( ( iy >= 0 ) && ( iy <= 49 ) )
      i = iy + 2000;
   else if ( ( iy >= 50 ) && ( iy <= 99 ) )
      i = iy + 1900;
   else
      i = iy;

/* Perform the conversion */
   slaClyd ( i, im, id, ny, nd, j );
}
void slaDjcal ( int ndp, double djm, int iymdf[4], int *j ) /*includefile*/
/*
**  - - - - - - - - -
**   s l a D j c a l
**  - - - - - - - - -
**
**  Modified Julian Date to Gregorian calendar, expressed
**  in a form convenient for formatting messages (namely
**  rounded to a specified precision, and with the fields
**  stored in a single array).
**
**  Given:
**     ndp      int       number of decimal places of days in fraction
**     djm      double    Modified Julian Date (JD-2400000.5)
**
**  Returned:
**     iymdf    int[4]    year, month, day, fraction in Gregorian calendar
**     *j       int       status:  nonzero = out of range
**
**  Any date after 4701BC March 1 is accepted.
**
**  Large ndp values risk internal overflows.  It is typically safe
**  to use up to ndp=4.
**
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**
**  Defined in slamac.h:  dmod
**
**  Last revision:   17 August 1999
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double fd, df, f, d;
   long jd, n4, nd10;

/* Validate */
   if ( ( djm <= -2395520.0 ) || ( djm >= 1.0e9 ) ) {
      *j = - 1;
      return;
   } else {

   /* Denominator of fraction */
      fd = pow ( 10.0, (double) gmax ( ndp, 0 ) );
      fd = dnint ( fd );

   /* Round date and express in units of fraction */
      df = djm * fd;
      df = dnint ( df );

   /* Separate day and fraction */
      f = dmod ( df, fd );
      if ( f < 0.0 ) f += fd;
      d = ( df - f ) / fd;

   /* Express day in Gregorian calendar */
      jd = (long) dnint ( d ) + 2400001L;
      n4 = 4L * ( jd + ( ( 2L * ( ( 4L * jd - 17918L ) / 146097L)
                                       * 3L ) / 4L + 1L ) / 2L - 37L );
      nd10 = 10L * ( ( ( n4 - 237L ) % 1461L ) / 4L ) + 5L;
      iymdf[0] = (int) ( ( n4 / 1461L ) - 4712L );
      iymdf[1] = (int) ( ( ( nd10 / 306L + 2L ) % 12L ) + 1L );
      iymdf[2] = (int) ( ( nd10 % 306L ) / 10L + 1L );
      iymdf[3] = (int) dnint ( f );
      *j = 0;
   }
}
void slaCldj ( int iy, int im, int id, double *djm, int *j ) /*includefile*/
/*
**  - - - - - - - -
**   s l a C l d j
**  - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  Given:
**     iy,im,id     int    year, month, day in Gregorian calendar
**
**  Returned:
**     *djm         double Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j           int    status:
**                           0 = OK
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  The year must be -4699 (i.e. 4700BC) or later.
**
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**
**  Last revision:   29 August 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long iyL, imL;

/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



/* Validate year */
   if ( iy < -4699 ) { *j = 1; return; }

/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *j = 2; return; }

/* Allow for leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;

/* Validate day */
   *j = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

/* Lengthen year and month numbers to avoid overflow */
   iyL = (long) iy;
   imL = (long) im;

/* Perform the conversion */
   *djm = (double)
        ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
        + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
        - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
        + (long) id - 2399904L );
}
void slaClyd(int iy, int im, int id, int *ny, int *nd,int*jstat)/*uncludefile*/
/*
**  - - - - - - - -
**   s l a C l y d
**  - - - - - - - -
**
**  Gregorian calendar to year and day in year (in a Julian calendar
**  aligned to the 20th/21st century Gregorian calendar).
**
**  Given:
**     iy,im,id     int    year, month, day in Gregorian calendar
**
**  Returned:
**     ny          int    year (re-aligned Julian calendar)
**     nd          int    day in year (1 = January 1st)
**     jstat       int    status:
**                          0 = OK
**                          1 = bad year (before -4711)
**                          2 = bad month
**                          3 = bad day (but conversion performed)
**
**  Notes:
**
**  1  This routine exists to support the low-precision routines
**     slaEarth, slaMoon and slaEcor.
**
**  2  Between 1900 March 1 and 2100 February 28 it returns answers
**     which are consistent with the ordinary Gregorian calendar.
**     Outside this range there will be a discrepancy which increases
**     by one day for every non-leap century year.
**
**  3  The essence of the algorithm is first to express the Gregorian
**     date as a Julian Day Number and then to convert this back to
**     a Julian calendar date, with day-in-year instead of month and
**     day.  See 12.92-1 and 12.95-1 in the reference.
**
**  Reference:  Explanatory Supplement to the Astronomical Almanac,
**              ed P.K.Seidelmann, University Science Books (1992),
**              p604-606.
**
**  Last revision:   26 November 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long i, j, k, l, n, iyL, imL;

/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



/* Validate year */
   if ( iy < -4711 ) { *jstat = 1; return; }

/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *jstat = 2; return; }

/* Allow for (Gregorian) leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;

/* Validate day */
   *jstat = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

/* Perform the conversion */
   iyL = (long) iy;
   imL = (long) im;
   i = ( 14 - imL ) /12L;
   k = iyL - i;
   j = ( 1461L * ( k + 4800L ) ) / 4L
     + ( 367L * ( imL - 2L + 12L * i ) ) / 12L
     - ( 3L * ( ( k + 4900L ) / 100L ) ) / 4L + (long) id - 30660L;
   k = ( j - 1L ) / 1461L;
   l = j - 1461L * k;
   n = ( l - 1L ) / 365L - l / 1461L;
   j = ( ( 80L * ( l - 365L * n + 30L ) ) / 2447L ) / 11L;
   i = n + j;
   *nd = 59 + (int) ( l -365L * i + ( ( 4L - n ) / 4L ) * ( 1L - j ) );
   *ny = (int) ( 4L * k + i ) - 4716;
}
