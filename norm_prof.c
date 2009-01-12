/* normalize prof by counts */
void norm_prof(float *prof, float *cnt, int nbins, int nifs, int nchans) /*includefile*/
{
  int i,s,c,n;
  float tmp;

  for (s=0; s<nbins; s++) {
    tmp=1.0;
    if(cnt[s] != 0.0) tmp=cnt[s];
    for (i=0; i<nifs; i++) {
      n=i*nchans;
      for (c=0; c<nchans; c++) prof[(n+c)*nbins+s] /= tmp;
    }
  }
}
/* denormalize prof by counts */
void denorm_prof(float *prof, float *cnt, int nbins, int nifs, int nchans) /*includefile*/
{
  int i,s,c,n;
  float tmp;

  for (s=0; s<nbins; s++) {
    tmp=1.0;
    if(cnt[s] != 0.0) tmp=cnt[s];
    for (i=0; i<nifs; i++) {
      n=i*nchans;
      for (c=0; c<nchans; c++) prof[(n+c)*nbins+s] *= tmp;
    }
  }
}
