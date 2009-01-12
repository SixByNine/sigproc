#include <math.h>
int bit(int bitindex, unsigned char byte) /*includefile*/
{
  unsigned char mask;
  mask=(unsigned char) pow(2,bitindex-1);
  
  if (mask & byte) {
    return (0);
  } else {
    return (1);
  }
}
