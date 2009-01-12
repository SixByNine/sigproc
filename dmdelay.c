/* 
   return the delay in seconds between two sky frequencies f1 and f2 (MHz) 
   N.B. constant of proportionality derived from e^2 pcm2 / (2pi c me) where 
   the elementary charge is assumed to be in Gaussian units: 4.8032068e-20
   and pc2m = 3.0856776e16 is the conversion between metres and parsecs 
*/
double dmdelay(double f1, double f2, double dm) /* includefile */
{
  return(4148.741601*((1.0/f1/f1)-(1.0/f2/f2))*dm);
}
