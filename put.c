#include <stdio.h>
void puti(int i) /*includefile*/
{
  printf("%d\n",i);
}
void putf(float f) /*includefile*/
{
  printf("%.f\n",f);
}
void putd(double d) /*includefile*/
{
  printf("%.12lf\n",d);
}
void putu(unsigned long l) /*includefile*/
{
  printf("%ld\n",l);
}
void putl(long l) /*includefile*/
{
  printf("%ld\n",l);
}
void putld(long double d) /*includefile*/
{
  printf("%llf\n",d);
}
