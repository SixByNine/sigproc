#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/* 
   pack two integers into a single char containing two 4-bit words 
   first integer passed (i) takes the lower four bits of the char,
   whilst the second integer passed (j) takes the higher four bits.
*/
unsigned char charof2ints (int i, int j) /* includefile */
{
  return ( (j<<4)+i );
}
/*
  reverse operation of the above routine - recovers original packed
  integers from within a character. i - low ; j - high as above.
*/
#define HI4BITS 240
#define LO4BITS   15
void char2ints (unsigned char c, int *i, int *j) /* includefile */
{
  *i =  c & LO4BITS;
  *j = (c & HI4BITS) >> 4;
}

#define HI2BITS 192
#define UPMED2BITS 48
#define LOMED2BITS 12
#define LO2BITS 3

void char2fourints (unsigned char c, int *i, int *j, int *k, int *l) /* includefile */
{
  *i =  c & LO2BITS;
  *j = (c & LOMED2BITS) >> 2;
  *k = (c & UPMED2BITS) >> 4;
  *l = (c & HI2BITS) >> 6;
}
