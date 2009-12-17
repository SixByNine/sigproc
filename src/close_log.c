#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
FILE *logfile;
int logging_mode;
void close_log() /*includefile*/
{
  if (logging_mode) fclose(logfile);
}
