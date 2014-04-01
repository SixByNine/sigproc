
#include <stdio.h>
/* define some functions for log message 
 * M.Keith 2012 - let me know if this fails to compile anywhere.
 * mkeith@pulsarastronomy.net
 **/
#define LOG_OUTFILE stderr
#define WHERESTR  "[%s:%d] "
#define WHEREARG  __FILE__, __LINE__
#define ENDL "\n"
#define WHEREERR "******\nERROR [%s:%d] "
#define WHERETCHK "[%s:%d] T=%.2f s: "
#define _LOG(...) fprintf(LOG_OUTFILE,__VA_ARGS__)
#define logmsg(_fmt, ...) _LOG(WHERESTR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logdbg(_fmt, ...)  if(debugFlag)logmsg(_fmt,##__VA_ARGS__)
#define logerr(_fmt, ...) _LOG(WHEREERR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logtchk(_fmt, ...) if(tcheck)_LOG(WHERETCHK _fmt ENDL, WHEREARG,(clock()-timer_clk)/(float)CLOCKS_PER_SEC,##__VA_ARGS__)



