#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
void error_message(char *message) /*includefile */
{
  fprintf(stderr,"ERROR: %s\n",message);
  exit(1);
}
