#include "mjk_random.h"
#include <stdlib.h>
#include <math.h>


uint32_t mjk_fallback_rand(int64_t *state){
   *state = *state*1103515245 + 12345;
   return ((uint32_t)(*state/INT32_MAX)%(INT32_MAX/2));
}

void mjk_rand_fill_fallback(mjk_rand_t *state){
   uint64_t iblk,ibuf;
   int64_t astate1 = *((int64_t*)(state->alg_state));
   int64_t astateP;
   int64_t astateN;

   const uint64_t buf_chunk = state->buffer_len/state->nthreadmax;
   const uint64_t chunk=1;

#pragma omp parallel shared(state,astateN) private(iblk,ibuf,astateP)
   {
#pragma omp for schedule(dynamic,chunk)
	  for(iblk=0; iblk < state->nthreadmax; iblk++){
		 astateP = astate1 + 89235*iblk;
		 for(ibuf=0; ibuf < buf_chunk; ibuf++){
			state->buffer[ibuf + iblk*buf_chunk] = mjk_fallback_rand(&astateP);
		 }
		 if (iblk==state->nthreadmax-1){
			astateN = astateP;
		 }
	  }
   }
   *((int64_t*)state->alg_state) = astateN;
   state->next=state->buffer_len-1;
}

double mjk_rand_gauss(mjk_rand_t *state) {
   double S, Z, U1, U2, u,v;
   double fac;

   if (state->gauss_state)
	  Z = state->gauss_next;
   else
   {
	  do
	  {
		 U1 = mjk_rand_double(state);
		 U2 = mjk_rand_double(state);

		 u = 2. * U1 - 1.;
		 v = 2. * U2 - 1.;
		 S = u * u + v * v;
	  } while(S >= 1);

	  fac = sqrt (-2. * log(S) / S);
	  Z = u * fac;
	  state->gauss_next=v*fac;
   }

   state->gauss_state = 1-state->gauss_state;

   return Z;
}



double mjk_rand_double(mjk_rand_t *state){
   uint32_t r = mjk_rand(state);
   return (double)r/(double)state->rmax;
}
uint32_t mjk_rand(mjk_rand_t *state){
   if (state->next < 0){
	  switch (state->alg){
		 case FALLBACK:
		 default:
			mjk_rand_fill_fallback(state);
	  }
   } 

   uint32_t ret = state->buffer[state->next];
   state->next--;
   return ret;
}


mjk_rand_t *mjk_rand_init(uint64_t seed){
   mjk_rand_t *state = calloc(1,sizeof(mjk_rand_t));
   state->next=-1;
   state->nthreadmax=16;
   state->buffer_len=1024*1024;
   state->alg=FALLBACK;
   state->alg_state=calloc(1,sizeof(int64_t));
   state->buffer=calloc(state->buffer_len,sizeof(uint32_t));
   state->rmax=INT32_MAX/2;
   state->gauss_next=0;
   state->gauss_state=0;
   *((int64_t*)(state->alg_state)) = seed;
   return state;
}


void mjk_rand_free(mjk_rand_t *state){
   free(state->buffer);
   free(state->alg_state);
   free(state);
}
