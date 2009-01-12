/* np2 - returns the nearest power of two number to the integer passed down */
#include<math.h>
int np2(int n) /*includefile*/
{
  int p=0;
  while( pow(2.0,(double)p) < n) p++;
  return(p-1);
}
