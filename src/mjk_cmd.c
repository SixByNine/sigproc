#include "mjk_cmd.h"
#include <stdlib.h>
#include <string.h>

#ifdef __MACH__
#include <sys/time.h>
//clock_gettime is not implemented on OSX
int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

void getArgs(int *argc, char** argv){
   int i;
   int n=1;
   for (i=1;i< *argc; i++){
	  if(argv[i][0]<5){
		 i+=argv[i][0];
		 continue;
	  }
	  argv[n]=argv[i];
	  n++;
   }
   *argc=n;
}
char* getS(char *lo, char* so, int argc, char** argv, char* val){
   int i;
   char* ret = val;
   for (i=1; i < argc; i++){
	  if (STREQ(argv[i],lo) || STREQ(argv[i],so)){
		 ret = *(argv+i+1);
		 argv[i][0]=1;
		 break;
	  }
   }
   return ret;
}

char getB(char *lo, char* so, int argc, char** argv, char val){
   int i;
   char ret = val;
   for (i=1; i < argc; i++){
	  if (STREQ(argv[i],lo) || STREQ(argv[i],so)){
		 ret = !val;
		 argv[i][0]=0;
		 break;
	  }
   }
   return ret;
}


double getF(char *lo, char* so, int argc, char** argv, double val){
   char* str = getS(lo,so,argc,argv,NULL);
   if (str==NULL) return val;
   else return atof(str);
}

int getI(char *lo, char* so, int argc, char** argv, int val){
   char* str = getS(lo,so,argc,argv,NULL);
   if (str==NULL) return val;
   else return atoi(str);
}





mjk_clock_t *init_clock(){
   mjk_clock_t *ret =  calloc(1,sizeof(mjk_clock_t));

#ifndef __MACH__
   if(sysconf(_SC_MONOTONIC_CLOCK)){
	  ret->spec =CLOCK_MONOTONIC;
   }else{
	  ret->spec =CLOCK_REALTIME;
   }
#endif
   return ret;
}

void start_clock(mjk_clock_t* clock){
   if(!clock->running){
	  clock_gettime(clock->spec,&clock->last_time);
	  clock->running=1;
   }
}
void stop_clock(mjk_clock_t* clock){
   if(clock->running){
	  long n = clock->last_time.tv_nsec;
	  time_t s = clock->last_time.tv_sec;
	  clock_gettime(clock->spec,&clock->last_time);
	  clock->running=0;
	  clock->n += (clock->last_time.tv_nsec-n);
	  clock->s += (clock->last_time.tv_sec-s);
   }
}
void reset_clock(mjk_clock_t* clock){
   clock->running=0;
   clock->n=0;
   clock->s=0;
   clock_gettime(clock->spec,&clock->last_time);
}
double read_clock(mjk_clock_t* clock){
   return clock->s + clock->n/1e9L;
}




