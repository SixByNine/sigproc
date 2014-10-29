#include "mjk_random.h"
#include "sigproc.h"
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


double a2281 (uint64_t seed[MJK_RAND_R1279_SZ], int *i){
   const int i1 = (*i - 2281) & MJK_RAND_R1279_SZ1;
   const int i2 = (*i - 1252) & MJK_RAND_R1279_SZ1;
   seed[*i] = seed[i1] + seed[i2];
   double ret = seed[*i] / (double)UINT64_MAX;
   //fprintf(stderr,"%lu %lu %lu %lg\n",seed[*i],seed[i1],seed[i2],ret);
   (*i) = (*i + 1)&MJK_RAND_R1279_SZ1;
   return ret;
}

double r1279 (uint64_t seed[MJK_RAND_R1279_SZ], int *i){
   const int i1 = (*i - 1279) & MJK_RAND_R1279_SZ1;
   const int i2 = (*i - 418) & MJK_RAND_R1279_SZ1;
   seed[*i] = seed[i1] * seed[i2]+1279;
   double ret = seed[*i] / (double)UINT64_MAX;
   //fprintf(stderr,"%lu %lu %lu %lg\n",seed[*i],seed[i1],seed[i2],ret);
   (*i) = (*i + 1)&MJK_RAND_R1279_SZ1;
   return ret;
}


void mjk_rand_fill(mjk_rand_t *state){
   time_t t0=time(NULL);
   logdbg("fill random buffer");
   uint64_t ibuf,iblk;
   const int blksize = state->buffer_len / state->nthreadmax;
#pragma omp parallel private(ibuf,iblk) shared(state)
#pragma omp for schedule(static,1)
   for(iblk=0; iblk < state->nthreadmax; iblk++){
	  for(ibuf=0; ibuf < blksize; ibuf++){
		 state->buffer[ibuf + iblk*blksize] = 
			a2281(
				  state->seed+iblk*MJK_RAND_R1279_SZ,
				  &(state->ir)[iblk]
				 );
		 //fprintf(stderr,"ttt %lg\n",state->buffer[ibuf]);
	  }
   }
   state->next=state->buffer_len-1;
   time_t t1=time(NULL);
   logdbg("buffer filled %ds",t1-t0);
}




#define TWO_PI 6.2831853071795864769252866

void mjk_rand_gauss2(double rand1, double rand2,double *o1,double *o2) {

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
#pragma omp parallel private(i,o1,o2,i1,i2)
#pragma omp for schedule(dynamic,chunk)
   for (i=0; i < state->buffer_len; i+=2){
	  i1=state->buffer[i];
	  i2=state->buffer[i+1];
	  mjk_rand_gauss2(i1,i2,&o1,&o2);
	  state->gauss_buffer[i]=o1;
	  state->gauss_buffer[i+1]=o2;
   }
   /*
	  for (i=1620; i < state->buffer_len; i++){ // =4096){
	  fprintf(stderr,"%d %g %lg\n",i,state->gauss_buffer[i],state->buffer[i]);
	  }
	  exit(1);
	  */
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
   long nseed = -(long)seed;
   mjk_rand_t *state = calloc(1,sizeof(mjk_rand_t));
   state->next=-1;
   state->nthreadmax=8;
   state->buffer_len=1024*16*state->nthreadmax;
   state->buffer=calloc(state->buffer_len,sizeof(double));
   state->rmax=INT32_MAX/2;
   state->gauss_next=-1;
   state->gauss_buffer=NULL;
   srand(seed);
   state->ir = calloc(state->nthreadmax,sizeof(int));
   state->seed = calloc(MJK_RAND_R1279_SZ*state->nthreadmax,sizeof(uint64_t));
   int i;
   uint64_t j,k;
   for(i=0; i < MJK_RAND_R1279_SZ*state->nthreadmax; i++){
	  j = (uint64_t) (UINT32_MAX * (double)nrran2(&nseed));
	  j = j << 8*sizeof(uint32_t);
	  k = (uint64_t) (UINT32_MAX * (double)nrran2(&nseed));
	  //fprintf(stderr,"j=%016llx k=%016llx r=%016llx  %g\n",j,k,j+k,nrran2(&nseed));
	  state->seed[i] = j+k;
   }
   mjk_rand_fill(state);
   mjk_rand_fill(state);
   return state;
}


void mjk_rand_free(mjk_rand_t *state){
   free(state->buffer);
   free(state->seed);
   if(state->gauss_buffer!=NULL)free(state->gauss_buffer);
   free(state);
}
