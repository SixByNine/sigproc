#include <inttypes.h>
#ifndef _mjk_rand_h
#define _mjk_rand_h
#ifdef __cplusplus
extern "C" {
#endif
   // the actual C header.
   typedef enum mjk_rand_alg {
	  FALLBACK,GNU48
   } mjk_rand_alg_t;


   typedef struct mjk_rand {
	  uint32_t nthreadmax;
	  uint32_t buffer_len;
	  uint32_t *buffer;
	  int32_t next;
	  mjk_rand_alg_t alg;
	  void* alg_state;
	  uint32_t rmax;
	  double gauss_next;
	  uint32_t gauss_state;
   } mjk_rand_t;


double mjk_rand_double(mjk_rand_t *state);
uint32_t mjk_rand(mjk_rand_t *state);
double mjk_rand_gauss(mjk_rand_t *state);
mjk_rand_t *mjk_rand_init(uint64_t seed);
void mjk_rand_free(mjk_rand_t *state);
#ifdef __cplusplus
}
#endif
#endif //_mjk_rand_h
