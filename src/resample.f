c==============================================================================
      subroutine resample(llog)

c==============================================================================
c
c   New routine to re-sample a time series o f ntim points contained
c   in the array series() to the rest frame of an object with a changing
c   velocity v.
c
c   Algorithm based on resample.f from sigproc-3.7 by Duncan Lorimer.
c
c   Code structure adapted from Haydon Knight's Rebin.ic.
c
c   This combination of code structure and algorithm appears to provide
c   correct acceleration warping to data sets 2^28 samples in length 
c   without floating point precision problems.
c
c   David Barnes, January 2007.   
c
c==============================================================================
      implicit none
      integer llog
      include 'seek.inc'
      include 'csamp.inc'

      real*8 refv

      real*8 c
      parameter(c=2.99792458e8)
      integer i

      integer samp_in, samp_out

      real*8 earthtime
      real*8 timethru,newpulsartime
      real*8 earlyweight,lateweight

      real*8 tau0,tint,taut,next,ac,ad

      refv = 0.0
      ac   = -1.0*refac
      ad   = -1.0*refad
      tint = ntim * tsamp
      tau0 = tsamp/(1. + ac * tint/2./c + ad*0.01*tint*tint/6./c)
      taut = tau0
      next = taut

      write(llog,*)'Re-sampling time series. AC=',refac,' m/s/s. AD=',
     &             refad,' cm/s/s/s'

c initialise to zero output series - possibly not necessary
      do i=1,ntim
         samp(i)=0.0
      enddo

      newpulsartime = 0.
      do samp_in=0,(ntim-1)

         earthtime = samp_in * tsamp
         timethru = newpulsartime 
         newpulsartime = newpulsartime + tau0*(1.+ac*earthtime/c +
     +        ad*0.01*earthtime*earthtime/c)

         samp_out = (timethru / tsamp)

         earlyweight = 1. - (timethru - samp_out*tsamp) / tsamp
         lateweight = 1. - earlyweight
         
         if (samp_out.ge.0.0.and.samp_out.lt.ntim) then
            samp(samp_out+1) = samp(samp_out+1) + earlyweight * 
     +           series(samp_in+1)
         endif
         if (samp_out.ge.-1.and.samp_out.lt.(ntim-1)) then
            samp(samp_out+1+1) = lateweight * series(samp_in+1)
         endif

      enddo

      do i=1,ntim
c uncomment following and do accel = 0 to verify output looks like input
c         if (mod(i, 2**22).eq.0) then
c           write (*,*) 'i, orig, new = ', i, series(i), samp(i)
c         endif
         series(i)=samp(i)
      enddo

      end
c==============================================================================
