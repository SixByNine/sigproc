#include "mjk_cmd.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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
   return ret;
}

void start_clock(mjk_clock_t* clock){
   if(!clock->running){
	  clock->stt = omp_get_wtime();
	  clock->running=1;
   }
}
void stop_clock(mjk_clock_t* clock){
   if(clock->running){
	  clock->running=0;
	  clock->time += omp_get_wtime() - clock->stt;
   }
}
void reset_clock(mjk_clock_t* clock){
   clock->running=0;
   clock->time=0;
   clock->stt=0;
}
double read_clock(mjk_clock_t* clock){
   return clock->time;
}



char* trim_string(char* str){
	int l = strlen(str)-1;
	while(l > 0){
		if(str[l]==' ')str[l]='\0';
		else break;
		l-=1;
	}
	int i=0;
	while (i < l){
		if(str[0]==' ')str++;
		else break;
		i+=1;
	}
	return str;
}


