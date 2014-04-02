#include "mjk_random.h"
#include <stdlib.h>
#include <math.h>
#include "mjklog.h"
#include <time.h>
#include <string.h>
#define debugFlag 0

uint32_t mjk_fallback_rand(int64_t *state){
   *state = *state*1103515245 + 12345;
   return ((uint32_t)(*state/INT32_MAX)%(INT32_MAX/2));
}

void mjk_rand_fill(mjk_rand_t *state){
   time_t t0=time(NULL);
   logdbg("fill random buffer");
   uint64_t iblk,ibuf;

   unsigned short astateP[3];
   const uint64_t buf_chunk = state->buffer_len/state->nthreadmax;
   const uint64_t chunk=64;

#pragma omp parallel shared(state) private(iblk,ibuf,astateP)
   {
#pragma omp for schedule(dynamic,chunk)
	  for(iblk=0; iblk < state->nthreadmax; iblk++){
		 astateP[0]+=iblk;
		 astateP[1]+=2*iblk;
		 astateP[2]-=3*iblk;
		 for(ibuf=0; ibuf < buf_chunk; ibuf++){
			state->buffer[ibuf + iblk*buf_chunk] = erand48(astateP);
		 }
		 if (iblk==state->nthreadmax-1){
			state->xi[0]=astateP[0];
			state->xi[1]=astateP[1];
			state->xi[2]=astateP[2];
		 }
	  }
   }
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

void mjk_rand_gauss2(mjk_rand_t *state,double rand1, double rand2,double *o1,double *o2) {

   if(rand1 < 1e-100) rand1 = 1e-100;
   rand1 = -2 * log(rand1);
   rand2 = rand2 * TWO_PI;

   *o1 = sqrt(rand1) * cos(rand2);
   *o2 = sqrt(rand1) * sin(rand2);
}

void mjk_rand_fill_gauss(mjk_rand_t *state){
   mjk_rand_fill(state);

   time_t t0=time(NULL);
   logdbg("fill gauss buffer");
   uint64_t i;
   double o1,o2;
   double i1,i2;
   const uint32_t chunk=64;
   if (state->gauss_buffer==NULL){
	  state->gauss_buffer = (float*)calloc(state->buffer_len,sizeof(float));
   }
#pragma omp parallel shared(state) private(i,o1,o2,i1,i2)
#pragma omp for  schedule(dynamic,chunk)
   for (i=0; i < state->buffer_len; i+=2){
	  i1=state->buffer[i];
	  i2=state->buffer[i+1];
	  mjk_rand_gauss2(state,i1,i2,&o1,&o2);
	  state->gauss_buffer[i]=o1;
	  state->gauss_buffer[i+1]=o2;
   }

   state->gauss_next=state->buffer_len-1;
   state->next=-1;

   time_t t1=time(NULL);
   logdbg("buffer filled %ds",t1-t0);
}


float mjk_rand_gauss(mjk_rand_t *state){
   float ret;
   {
	  if (state->gauss_next < 0){
		 mjk_rand_fill_gauss(state);
	  } 
	  ret = state->gauss_buffer[state->gauss_next];
	  state->gauss_next--;
   }


   return ret;
}



uint32_t mjk_rand(mjk_rand_t *state){
   double r = mjk_rand_double(state);
   return (uint32_t)(r*state->rmax);
}

double mjk_rand_double(mjk_rand_t *state){
   double ret;
   {
	  if (state->next < 0){
		 mjk_rand_fill(state);
	  } 
	  ret = state->buffer[state->next--];
   }
   return ret;
}

void mjk_rand_gauss_atleast(mjk_rand_t *state, uint64_t n){
   if(state->gauss_next < n){
	  mjk_rand_fill_gauss(state);
   }
}

mjk_rand_t *mjk_rand_init(uint64_t seed){
   mjk_rand_t *state = calloc(1,sizeof(mjk_rand_t));
   state->next=-1;
   state->nthreadmax=1024*1024;
   state->buffer_len=state->nthreadmax*16;
   state->buffer=calloc(state->buffer_len,sizeof(double));
   state->rmax=INT32_MAX/2;
   state->gauss_next=-1;
   state->gauss_buffer=NULL;
   memcpy(state->xi,&seed,6);
   logdbg("%x %x %x 0x%016Lx",state->xi[0],state->xi[1],state->xi[2],seed);
   return state;
}


void mjk_rand_free(mjk_rand_t *state){
   free(state->buffer);
   if(state->gauss_buffer!=NULL)free(state->gauss_buffer);
   free(state);
}
