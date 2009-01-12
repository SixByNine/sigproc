/* 
   open_log.c - checks for monitor.running and enables logging if it exists 
   global variable logging mode is switched off otherwise so that update_log
   and close_log routines will not do anything if no logging is required.
*/
#include <stdio.h>
#include "sigproc.h"
FILE *logfile;
int logging_mode;
void open_log(char *filename) /*includefile*/
{
  if (file_exists("monitor.running")) 
    logging_mode=1;
  else
    logging_mode=0;
  if (logging_mode) logfile=open_file(filename,"w");
}
