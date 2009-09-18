/* returns the minimum value of a vector */
float vmin(float *vec, float n) /*includefile*/
{
  float minvalue; int i;
  minvalue=vec[0];
  for (i=1;i<n;i++) if (vec[i]<minvalue) minvalue=vec[i];
  return (minvalue);
}

