c=============================================================================
      subroutine fftwdata(llog)
c=============================================================================
c
c     FFTs the time series using a call to the FFTW routines. 
c     - Ralph Eatough 12/07
c    
c     llog    - i4  - llogical unit number for all but warning messages.
c      
c=============================================================================
c
      implicit none
      include 'seek.inc'
      integer n,llog  !,slun
      integer*8 plan
      integer FFTW_ESTIMATE
      parameter (FFTW_ESTIMATE=64)

      n=ntim
      write(llog,*) 'FFT: (fftw-3.1.2)...'

c      call import_wisdom_from_file(isuccess, slun)
      call sfftw_plan_dft_r2c_1d(plan,n,series,series,FFTW_ESTIMATE)
      call sfftw_execute(plan)
      call sfftw_destroy_plan(plan)

c      call sglfft(series,series(2),n,n,n,2)
c      call realtr(series,series(2),n,2)

c      ntim=n*2
      end
c      
c=============================================================================
