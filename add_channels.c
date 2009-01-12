void add_channels(float *data, int nsamples, int nadd) /* includefile */
{
  int i,j,k;
  float sum,n;

  if (nadd <= 1) return;

  j=k=0;
  sum=0.0;
  n=(float) nadd;

  for (i=0; i<nsamples; i++) {
    sum+=data[i];
    j++;
    if (j==nadd) {
      data[k]=sum/n;
      k++; 
      j=0;
      sum=0.0;
    }
  }

}
