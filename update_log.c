#include <stdio.h>
FILE *logfile;
char inpfile[80], outfile[80];
int logging_mode;
void update_log(char *string) /* includefile */
{
  if (logging_mode) {
    fprintf(logfile,"input %s status %s output %s\n",inpfile,string,outfile);
    fflush(logfile);
  }
}
