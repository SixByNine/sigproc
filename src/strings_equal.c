#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
int strings_equal (char *string1, char *string2) /* includefile */
{
  if (!strcmp(string1,string2)) {
    return 1;
  } else {
    return 0;
  }
}
