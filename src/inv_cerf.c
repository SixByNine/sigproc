#include <math.h>
/* Approximation for Inverse Complementary Error Function */
double inv_cerf(double input) /* includefile */
{
    static double numerator_const[3] = {1.591863138,-2.442326820,0.37153461};
    static double denominator_const[3] = {1.467751692,-3.013136362,1.0};
    double temp_data,temp_data_srq,erf_data;
    erf_data = 1.0-input;
    temp_data = (erf_data*erf_data - 0.5625);
    temp_data_srq = temp_data*temp_data;
    temp_data = (erf_data)*(numerator_const[0]+(temp_data*numerator_const[1])+(temp_data_srq*numerator_const[2]))/
(denominator_const[0]+(temp_data*denominator_const[1])+(temp_data_srq*denominator_const[2]));
    return(temp_data);
}
