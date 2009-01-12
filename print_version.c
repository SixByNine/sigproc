#include "version.h"
void print_version(char *program, char *argument) /*includefile*/
{
  if ( (strings_equal(argument,"version")) ||
       (strings_equal(argument,"-version"))) {
    printf("PROGRAM: %s is part of SIGPROC version: %.1f\n",program,SIGPROC_VERSION);
    exit(0);
  }
}
