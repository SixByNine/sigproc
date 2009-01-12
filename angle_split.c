void angle_split(double angle, int *dd, int *mm, double *ss) /*includefile*/
{
  int negative;
  if (angle<0.0) {
    angle*=-1.0;
    negative=1;
  } else {
    negative=0;
  }
  *dd=(int) (angle/10000.0);
  angle-=(double) (*dd)*10000.0;
  *mm=(int) (angle/100.0);
  *ss=angle-100.0*(*mm);
  if (negative) *dd = *dd * -1;
}

