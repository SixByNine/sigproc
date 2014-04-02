#include "mjk_random.h"
#include <stdlib.h>
#include <math.h>
#include "mjklog.h"
#include <time.h>
#define debugFlag 0

uint32_t mjk_fallback_rand(int64_t *state){
	*state = *state*1103515245 + 12345;
	return ((uint32_t)(*state/INT32_MAX)%(INT32_MAX/2));
}

void mjk_rand_fill_fallback(mjk_rand_t *state){
	time_t t0=time(NULL);
	logdbg("fill random buffer");
	uint64_t iblk,ibuf;
	int64_t astate1 = *((int64_t*)(state->alg_state));
	int64_t astateP;
	int64_t astateN;

	const uint64_t buf_chunk = state->buffer_len/state->nthreadmax;
	const uint64_t chunk=64;

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
	time_t t1=time(NULL);
	logdbg("buffer filled %ds",t1-t0);
}


int mjk_rand_gauss2_old(mjk_rand_t *state,double U1, double U2,double *o1,double *o2) {
	double S, Z, u,v;
	double fac;

	u = 2. * U1 - 1.;
	v = 2. * U2 - 1.;
	S = u * u + v * v;
	if (S < 1){
		fac = sqrt (-2. * log(S) / S);
		*o1 = u * fac;
		*o2 = v * fac;
		return 0;
	} else {
		return 1;
	}


}

#define TWO_PI 6.2831853071795864769252866

int mjk_rand_gauss2(mjk_rand_t *state,double rand1, double rand2,double *o1,double *o2) {

   if(rand1 < 1e-100) rand1 = 1e-100;
   rand1 = -2 * log(rand1);
   rand2 = rand2 * TWO_PI;

   *o1 = sqrt(rand1) * cos(rand2);
   *o2 = sqrt(rand1) * sin(rand2);
}

void mjk_rand_fill_fallback_gauss(mjk_rand_t *state){
   mjk_rand_fill_fallback(state);

   time_t t0=time(NULL);
   logdbg("fill gauss buffer");
   uint64_t i;
   double o1,o2,r1,r2;
   int64_t i1,i2;
   const uint32_t chunk=64;
   if (state->gauss_buffer==NULL){
	  state->gauss_buffer = (float*)calloc(state->buffer_len,sizeof(float));
   }
#pragma omp parallel shared(state) private(i,o1,o2,i1,i2,r1,r2)
#pragma omp for  schedule(dynamic,chunk)
   for (i=0; i < state->buffer_len/2; i+=2){
	  i1=state->buffer[i];
	  i2=state->buffer[i+1];
	  r1 = (double)i1/(double)state->rmax;
	  r2 = (double)i2/(double)state->rmax;
	  while(mjk_rand_gauss2(state,r1,r2,&o1,&o2)){
		 r1 = (double)mjk_fallback_rand(&i1)/(double)state->rmax;
		 r2 = (double)mjk_fallback_rand(&i1)/(double)state->rmax;
	  }
	  state->gauss_buffer[i]=o1;
	  state->gauss_buffer[i+1]=o2;
   }

   state->gauss_next=state->buffer_len/2-1;
   state->next=-1;

   time_t t1=time(NULL);
   logdbg("buffer filled %ds",t1-t0);
}


float mjk_rand_gauss(mjk_rand_t *state){
   float ret;
#pragma omp critical (mjk_rand_gauss_CRITICAL)
   {
	  if (state->gauss_next < 0){
		 switch (state->alg){
			case FALLBACK:
			default:
			   mjk_rand_fill_fallback_gauss(state);
		 }
	  } 
	  ret = state->gauss_buffer[state->gauss_next];
	  state->gauss_next--;
   }


   return ret;
}



double mjk_rand_double(mjk_rand_t *state){
   uint32_t r = mjk_rand(state);
   return (double)r/(double)state->rmax;
}

uint32_t mjk_rand(mjk_rand_t *state){
   uint32_t ret;
#pragma omp critical (mjk_rand_CRITICAL)
   {
	  if (state->next < 0){
		 switch (state->alg){
			case FALLBACK:
			default:
			   mjk_rand_fill_fallback(state);
		 }
	  } 
	  ret = state->buffer[state->next--];
   }
   return ret;
}

void mjk_rand_gauss_atleast(mjk_rand_t *state, uint64_t n){
   if(state->gauss_next < n){
	  mjk_rand_fill_fallback_gauss(state);
   }
}

mjk_rand_t *mjk_rand_init(uint64_t seed){
   mjk_rand_t *state = calloc(1,sizeof(mjk_rand_t));
   state->next=-1;
   state->nthreadmax=1024*1024;
   state->buffer_len=state->nthreadmax*64;
   state->alg=FALLBACK;
   state->alg_state=calloc(1,sizeof(int64_t));
   state->buffer=calloc(state->buffer_len,sizeof(uint32_t));
   state->rmax=INT32_MAX/2;
   state->gauss_next=-1;
   state->gauss_buffer=NULL;
   *((int64_t*)(state->alg_state)) = seed;
   return state;
}


void mjk_rand_free(mjk_rand_t *state){
   free(state->buffer);
   if(state->gauss_buffer!=NULL)free(state->gauss_buffer);
   free(state->alg_state);
   free(state);
}
