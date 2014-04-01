#ifndef _mjk_cmd_h
#define _mjk_cmd_h
#define STREQ(a,b) (strcmp(a,b)==0)
char* getS(char *lo, char* so, int argc, char** argv, char* val);

char getB(char *lo, char* so, int argc, char** argv, char val);

double getF(char *lo, char* so, int argc, char** argv, double val);
int getI(char *lo, char* so, int argc, char** argv, int val);


#endif
