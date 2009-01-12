#include <time.h>
#include <stdio.h>
#include <stdlib.h>
/* ssm - returns the approximate number of seconds since midnight */
int ssm(void) /*includefile*/
{
  char hours[8],minutes[8],seconds[8];
  struct tm *ptr;
  time_t lt;

  lt = time(NULL);
  ptr= localtime(&lt);
  strftime(hours,8,"%H",ptr);
  strftime(minutes,8,"%M",ptr);
  strftime(seconds,8,"%S",ptr);

  return (atoi(hours)*3600 + atoi(minutes)*60 + atoi(seconds));
}

