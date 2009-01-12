#include <stdio.h>
#include "sigproc.h"
/* 
## help on pspm_chans
## int *pspm_chans - returns pointer to a look-up table ordering PSPM channels
## 
## variables passed down:
## 
## int nchans - the number of PSPM channels
## 
## see Brian Cadwell's PhD thesis p. 118 for rationale behind the ordering
## 
## Last modified: Mar 27, 2001 (dunc@naic.edu)
## help end
*/
int *pspm_chans(int nchans) /* includefile */
{
  int i,j,k, *table;
  table=(int *) malloc(nchans*sizeof(int));
  i=0;
  for (j=nchans/4;j>=1;j--) for (k=4;k>0;k--) table[i++]=j*4-k;
  return (table);
}

#include "bpphdr.h"
static int dfb_chan_lookup[MAXREGS][NIBPERREG] = {
  {4, 0, 4, 0},
  {5, 1, 5, 1},
  {6, 2, 6, 2}, 
  {7, 3, 7, 3},
  {4, 0, 4, 0},
  {5, 1, 5, 1},
  {6, 2, 6, 2},
  {7, 3, 7, 3}
};
  
/* This takes care of byte swap in outreg_b */
static float sideband_lookup[MAXREGS][NIBPERREG] = {
  {-1.0, -1.0, +1.0, +1.0},
  {-1.0, -1.0, +1.0, +1.0},
  {-1.0, -1.0, +1.0, +1.0},
  {-1.0, -1.0, +1.0, +1.0},
  {+1.0, +1.0, -1.0, -1.0},
  {+1.0, +1.0, -1.0, -1.0},
  {+1.0, +1.0, -1.0, -1.0},
  {+1.0, +1.0, -1.0, -1.0}
};

double fch1,foff,fmid;
int nifs, nchans;

int  *bpp_chans(double bw, int mb_start_addr, int mb_end_addr, int mb_start_brd, int mb_end_brd, int *cb_id, double *aib_los, float *dfb_sram_freqs, double rf_lo) /* includefile */
{
  int i, n=0, dfb_chan, logical_brd, regid, bid, nibble, *table;
  double  f_aib, u_or_l, f_sram, fc;
  float *fmhz;
  unsigned long *nridx;
  double rf_lo_mhz;

  nchans = (mb_end_addr/2-mb_start_addr/2+1)*(mb_end_brd-mb_start_brd+1)*4;
  fmhz   = (float *) malloc(nchans*sizeof(float));
  table  = (int *)   malloc(nchans*sizeof(int));
  nridx  = (unsigned long *) malloc(nchans*sizeof(unsigned long));
  if (-1.e6<rf_lo && rf_lo<1.e6) 
    rf_lo_mhz = rf_lo;
  else
    rf_lo_mhz = rf_lo/1.e6;

  /* 
     Loop over (16-bit) regs per board. divide by 2's are to make them 
     word addresses instead of byte addresses so we can index with them.
     Normal modes will be regid = 0..3, 0..7, or 4..7 
  */

  for (regid=mb_start_addr/2; regid<=mb_end_addr/2; regid++)
    /* Loop over each board */
    for (bid=mb_start_brd;bid<=mb_end_brd;bid++) {
      /* Now find which LOGICAL CB we are reading */
      logical_brd = -1;
      for (i=0; i<MAXNUMCB; i++) {
        if (bid == cb_id[i]) {
	  logical_brd = i;
	  break;
        }
      }
      if (logical_brd == -1) error_message("bpp_chan - logical_brd not found");
      /* Assumes cabling so that LO0 feeds MF0,1 which feeds leftmost CB! */
      f_aib = aib_los[logical_brd];
      /* Loop over 4 nibbles per reg */
      for (nibble=0; nibble<4; nibble++) {
        dfb_chan = dfb_chan_lookup[regid][nibble];
        u_or_l = sideband_lookup[regid][nibble];
        f_sram = dfb_sram_freqs[dfb_chan];
        fc = f_aib + f_sram + u_or_l * bw/4.0;
	if (rf_lo_mhz<1.e4) /* below 10 GHz LSB; above 10 GHz USB */
	  fmhz[n++]=rf_lo_mhz+800-fc/1.0e6;
	else
	  fmhz[n++]=rf_lo_mhz+fc/1.0e6;
      }
    }

  /* produce lookup table which gives channels in order of descending freq */
  indexx(96,fmhz-1,nridx-1);
  if (nchans==192) {
    nifs=2;
    indexx(96,fmhz+96-1,nridx+96-1);
    for (i=96; i<192; i++) 
      nridx[i] += 96;
  }
  n=nchans;
  for (i=0;i<nchans;i++) 
    table[i]=nridx[--n]-1;    
  nchans/=nifs;
  fmid=0.5*(fmhz[table[0]]+fmhz[table[95]]);
  fch1=fmhz[table[0]];
  foff=fch1-fmhz[table[1]];
  free(fmhz);
  free(nridx);
  return(table);
}
