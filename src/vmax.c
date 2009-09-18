/* returns the maximum value of a vector */
float vmax(float *vec, float n) /*includefile*/
{
  float maxvalue; int i;
  maxvalue=vec[0];
  for (i=1;i<n;i++) if (vec[i]>maxvalue) maxvalue=vec[i];
  return (maxvalue);
}
