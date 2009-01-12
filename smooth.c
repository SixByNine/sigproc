/*
  This routine will smooth a time series by adding adjacent samples.

  v. 0.0  TJWL  27 Sept 1994
  v. 0.1  TJWL  30 Sept 1994 - minor bug fixes

  Input:
  int ndata;     the number of time series samples
  float data[];  the time series

  Output:
  int ndata;     the number of smoothed time series samples
  float data[];  the smoothed time series    

*/

void smooth(ndata, data)
int *ndata;
float data[];
{
  int i, npoints;

  npoints = *ndata - 1;

  for (i=0; i<npoints; i++)
    data[i] = (data[i] + data[i+1])/2.0;

  *ndata = npoints;

  return;
}

