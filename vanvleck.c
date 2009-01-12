#include <math.h>
#define NO    0
#define YES   1
/*------------------------------------------------------------------------*
 * Van Vleck Correction for 9-level sampling/correlation
 *  Samples {-4,-3,-2,-1,0,1,2,3,4}
 * Uses Zerolag to adjust correction
 *   data_array -> Points into autocorrelation function of at least 'count' points
 * this routine takes the first value as the zerolag and corrects the remaining
 * count-1 points.  Zerolag is set to a normalized 1
 * NOTE - The available routine works on lags normaized to -16<rho<16, so
 *  I need to adjust the values before/after the fit
 * Coefficent ranges
 *   c1 
 *     all
 *   c2
 *     r0 > 4.5
 *     r0 < 2.1
 *     rest 
 * NOTE - correction is done INPLACE ! Original values are destroyed
 * As reported by M. Lewis -> the polynomial fits are OK, but could be improved
 *------------------------------------------------------------------------*/
void vanvleck9lev(double *rho,int npts) /* includefile */
{
  static double coef1[5] = { 1.105842267, -0.053258115, 0.011830276,-0.000916417, 0.000033479 };

  static double coef2rg4p5[5] = { 0.111705575, -0.066425925, 0.014844439, -0.001369796, 0.000044119 };
  static double coef2rl2p1[5] =  { 1.285303775, -1.472216011, 0.640885537, -0.123486209, 0.008817175 };
  static double coef2rother[5] = { 0.519701391, -0.451046837, 0.149153116, -0.021957940, 0.001212970 };

  static double coef3rg2p0[5] = { 1.244495105, -0.274900651, 0.022660239, -0.000760938, -1.993790548 };
  static double coef3rother[5] = { 1.249032787, 0.101951346, -0.126743165, 0.015221707, -2.625961708 };

  static double coef4rg3p15[5] = { 0.664003237, -0.403651682, 0.093057131, -0.008831547, 0.000291295 };
  static double coef4rother[5] = { 9.866677289, -12.858153787, 6.556692205, -1.519871179, 0.133591758 };

  static double coef5rg4p0[4] = { 0.033076469, -0.020621902, 0.001428681, 0.000033733};
  static double coef5rg2p2[4] = { 5.284269565, 6.571535249, -2.897741312, 0.443156543};
  static double coef5rother[4] = {-1.475903733, 1.158114934, -0.311659264, 0.028185170};

  double acoef[5],dtmp,zl,ro;
  int i;
  ro = *rho;
  zl = (*rho)*16;
/*  for(i=0; i<npts; i++) 
    *(rho+i) *= *(rho+i)/ro; */
 
  
  acoef[0] = ((((coef1[4]*zl + coef1[3])*zl + coef1[2])*zl +coef1[1])*zl +coef1[0]);

  if( zl > 4.50)
    acoef[1] =((((coef2rg4p5[4]*zl + coef2rg4p5[3])*zl + coef2rg4p5[2])*zl + coef2rg4p5[1])*zl + coef2rg4p5[0]);
  else if(zl < 2.10)  
    acoef[1] =((((coef2rl2p1[4]*zl + coef2rl2p1[3])*zl + coef2rl2p1[2])*zl + coef2rl2p1[1])*zl + coef2rl2p1[0]);
  else
    acoef[1] =((((coef2rother[4]*zl + coef2rother[3])*zl + coef2rother[2])*zl + coef2rother[1])*zl + coef2rother[0]);

  if( zl > 2.00)
    acoef[2] = coef3rg2p0[4]/zl + (((coef3rg2p0[3]*zl + coef3rg2p0[2])*zl + coef3rg2p0[1])*zl + coef3rg2p0[0]);
  else
    acoef[2] = coef3rother[4]/zl + (((coef3rother[3]*zl + coef3rother[2])*zl + coef3rother[1])*zl + coef3rother[0]);  
  
  if( zl > 3.15)
    acoef[3] = ((((coef4rg3p15[4]*zl + coef4rg3p15[3])*zl + coef4rg3p15[2])*zl + coef4rg3p15[1])*zl + coef4rg3p15[0]);
  else
    acoef[3] = ((((coef4rg3p15[4]*zl + coef4rother[3])*zl + coef4rother[2])*zl + coef4rother[1])*zl + coef4rother[0]);  

  if( zl > 4.00)
    acoef[4] =(((coef5rg4p0[3]*zl + coef5rg4p0[2])*zl + coef5rg4p0[1])*zl + coef5rg4p0[0]);
  else if(zl < 2.2)  
    acoef[4] =(((coef5rg2p2[3]*zl + coef5rg2p2[2])*zl + coef5rg2p2[1])*zl + coef5rg2p2[0]);
  else
    acoef[4] =(((coef5rother[3]*zl + coef5rother[2])*zl + coef5rother[1])*zl + coef5rother[0]);
  
 /* printf("(0)-> %g (1)-> %g (2)-> %g (3)-> %g (4)-> %g\n",acoef[0],acoef[1],acoef[2],acoef[3],acoef[4]);*/
  for(i=1; i<npts; i++)
  {
    dtmp = *(rho+i);
    *(rho+i) =((((acoef[4]*dtmp + acoef[3])*dtmp + acoef[2])*dtmp + acoef[1])*dtmp + acoef[0])*dtmp;
  }
  *(rho) = 1.0;
  return;
}
/*------------------------------------------------------------------------*
 * Van Vleck Correction for 3-level sampling/correlation
 *  Samples {-1,0,1}
 * Uses Zerolag to adjust correction
 *   data_array -> Points into autocorrelation function of at least 'count' points
 * this routine takes the first value as the zerolag and corrects the remaining
 * count-1 points.  Zerolag is set to a normalized 1
 * 
 * NOTE - correction is done INPLACE ! Original values are destroyed
 *------------------------------------------------------------------------*/
int vanvleck3lev(double *rho,int npts) /* includefile */
{
    static double lo_const[3][4] = {
         { 0.939134371719, -0.567722496249, 1.02542540932, 0.130740914912 },
         { -0.369374472755, -0.430065136734, -0.06309459132, -0.00253019992917},
         { 0.888607422108, -0.230608118885, 0.0586846424223, 0.002012775510695}
         };
    static double high_const[5][4] = {
         {-1.83332160595, 0.719551585882, 1.214003774444, 7.15276068378e-5},
         {1.28629698818, -1.45854382672, -0.239102591283, -0.00555197725185},
         {-7.93388279993, 1.91497870485, 0.351469403030, 0.00224706453982},
         {8.04241371651, -1.51590759772, -0.18532022393, -0.00342644824947},
         {-13.076435520, 0.769752851477, 0.396594438775, 0.0164354218208}
         };

    double lo_u[3],lo_h[3];
    double high_u[5],high_h[5];
    double lo_coefficient[3];
    double high_coefficient[5];
    double zho,zho_3;
    double temp_data;
    double temp_data_1;
    int  ichan,ico,flag_any_high;
/* Perform Lo correction on All data that is not flaged for high correction --*/
    zho=(double)rho[0];
    
    zho_3=zho*zho*zho;

    lo_u[0]=zho;
    lo_u[1]=zho_3-(61.0/512.0);
    lo_u[2]=zho-(63.0/128.0);

    lo_h[0]=zho*zho;
    lo_h[2]=zho_3*zho_3*zho;         /* zlag ^7 */
    lo_h[1]=zho*lo_h[2];               /* zlag ^8 */
/* determine lo-correct coefficents -*/
    for(ico=0; ico<3; ico++)
    {
         lo_coefficient[ico]=(lo_u[ico]*(lo_u[ico]*(lo_u[ico]*lo_const[ico][0]+lo_const[ico][1])+lo_const[ico][2])+lo_const[ico][3])/lo_h[ico];
    }
/* perform correction --*/
    for(ichan=1,flag_any_high=NO; ichan<npts; ichan++)
    {
        temp_data=(double)rho[ichan];
        if(fabs(temp_data) > 0.199)
        {
            if(flag_any_high==NO)
            {
                   high_u[0]=lo_h[2];                  /* zlag ^7 */
                   high_u[1]=zho-(63.0/128.0);
                   high_u[2]=zho*zho-(31.0/128.0);
                   high_u[3]=zho_3-(61.0/512.0);
                   high_u[4]=zho-(63.0/128.0);

                   high_h[0]=lo_h[1];               /* zlag ^8 */
                   high_h[1]=lo_h[1];               /* zlag ^8 */
                   high_h[2]=lo_h[1]*zho_3*zho;   /* zlag ^12 */
                   high_h[3]=lo_h[1]*lo_h[1]*zho;  /* zlag ^17 */
                   high_h[4]=high_h[3];             /* zlag ^17 */
                   for(ico=0; ico<5; ico++)
                   {
                      high_coefficient[ico]=(high_u[ico]*(high_u[ico]*(high_u[ico]*high_const[ico][0]+high_const[ico][1])+high_const[ico][2])+high_const[ico][3])/high_h[ico];
                   }
                   flag_any_high=YES;
              }
              temp_data_1=fabs(temp_data*temp_data*temp_data);
              rho[ichan]=(temp_data*(temp_data_1*(temp_data_1*(temp_data_1*(temp_data_1*high_coefficient[4]+high_coefficient[3])+high_coefficient[2])+high_coefficient[1])+high_coefficient[0]));
         } else
         {
              temp_data_1=temp_data*temp_data;
              rho[ichan]=(temp_data*(temp_data_1*(temp_data_1*lo_coefficient[2]+lo_coefficient[1])+lo_coefficient[0]));
         }
    }
    rho[0] = 1.0;
    return(0);
}

