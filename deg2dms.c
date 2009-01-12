double deg2dms(double angle) /*includefile*/
{
  int deg,min;
  double sec,sign;

  if (angle<0.0) 
    sign=-1.0;
  else
    sign=+1.0;

  angle=fabs(angle);
  deg=(int)angle;
  angle-=(double)deg;
  angle*=60.0;
  min=(int)angle;
  angle-=(double)min;
  angle*=60.0;
  sec=angle;
  return(sign*((double)deg*10000.0+(double)min*100.0+sec));
}
double h2hms(double hours) /*includefile*/
{
  double hms;
  hms=deg2dms(hours);
  return(hms);
}

