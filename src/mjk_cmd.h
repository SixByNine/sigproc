#include <omp.h>
#ifndef _mjk_cmd_h
#define _mjk_cmd_h
#define STREQ(a,b) (strcmp(a,b)==0)

#ifdef __cplusplus
extern "C" {
#endif



char* getS(char *lo, char* so, int argc, char** argv, char* val);

char getB(char *lo, char* so, int argc, char** argv, char val);

double getF(char *lo, char* so, int argc, char** argv, double val);
int getI(char *lo, char* so, int argc, char** argv, int val);

void getArgs(int *argc, char** argv);

typedef struct mjk_clock{
   char running;
   double time;
   double stt;
} mjk_clock_t;

mjk_clock_t *init_clock();
void start_clock(mjk_clock_t* clock);
void stop_clock(mjk_clock_t* clock);
void reset_clock(mjk_clock_t* clock);
double read_clock(mjk_clock_t* clock);

char* trim_string(char* str);

#ifdef __cplusplus
}
#endif


#endif
