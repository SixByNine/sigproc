#include "mjk_cmd.h"
#include <stdlib.h>

char* getS(char *lo, char* so, int argc, char** argv, char* val){
   int i;
   char* ret = val;
   for (i=1; i < argc; i++){
	  if (STREQ(argv[i],lo) || STREQ(argv[i],so)){
		 ret = *(argv+i+1);
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

