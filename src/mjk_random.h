#include <inttypes.h>
#include <stdlib.h>
#ifndef _mjk_rand_h
#define _mjk_rand_h
#define MJK_RAND_R1279_SZ 2048
#define MJK_RAND_R1279_SZ1 2047
#ifdef __cplusplus
extern "C" {
#endif



   typedef struct mjk_rand {
	   uint32_t nthreadmax;
	   uint32_t buffer_len;
	   double *buffer;
	   int32_t next;
	   uint32_t rmax;
	   int32_t gauss_next;
	   float* gauss_buffer;
	   uint64_t* seed;
	   int32_t* ir;
   } mjk_rand_t;


   double mjk_rand_double(mjk_rand_t *state);
   uint32_t mjk_rand(mjk_rand_t *state);
   float mjk_rand_gauss(mjk_rand_t *state);
   mjk_rand_t *mjk_rand_init(uint64_t seed);
   void mjk_rand_free(mjk_rand_t *state);
void mjk_rand_gauss_atleast(mjk_rand_t *state, uint64_t n);
#ifdef __cplusplus
}
#endif
#endif //_mjk_rand_h
