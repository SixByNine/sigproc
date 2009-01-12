void add_samples(float *data, int nifs, int nchans, int nadd) /* includefile */
{
  int i,c,t,inc,nxc;
  float n;

  if (nadd <= 1) return;

  nxc=nifs*nchans;
  for (t=1; t<nadd; t++) 
    for (i=0; i<nifs; i++) {
      inc=i*nchans;
      for (c=0; c<nchans; c++) data[inc+c]+=data[t*nxc+inc+c];
    }

/*
  commented out for now as this reduces dynamic range when summing 1bit data
  n=(float) nadd;
  for (i=0; i<nxc; i++) data[i]/=n;
*/

}
