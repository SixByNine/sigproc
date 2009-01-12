c==============================================================================
      subroutine resample(llog)
c==============================================================================
c
c   Routine to re-sample a time series of ntim points contained in the array
c   series() to the rest frame of an object with a changing velocity v.
c      
c   The required sampling time is given by: tau(t)=tau0(1.0+v/c)
c
c   Two approximations to v are available:
c
c     - constant acceleration: v=a*t
c     - first derivative:      v=a*t + 0.5*adot*t*t
c     
c   The acceleration terms are stored in the common real*4's refac and refad.
c
c     refac - reference acceleration in m/s/s
c     refad - reference a-dot in cm/s/s/s     (cm used for convenience)
c
c   The constant tau0 is normalised so that the final time series
c   contains the same number of samples as the original one. The array
c   samp() is used as a buffer to store the re-sampled time series.
c
c     Adapted from rebin2.f 98/07/12 (dunc@mpifr-bonn.mpg.de)
c
c==============================================================================
      implicit none
      integer llog
      include 'seek.inc'
      include 'csamp.inc'
      real*8 tint,tau0,time,taut,next,tav,nav,c,ac,ad
      parameter(c=2.99792458e8)
      integer i,j

      if (refac.eq.0.0.and.refad.eq.0.0) return
      
      nav=0.0
      tav=0.0
      do i=1,ntim
         samp(i)=0.0
      enddo
      i=1
      j=1
      time=0.0
      tint=real(ntim)*tsamp
      ac=refac
      ad=refad*0.01 ! convert from cm/s/s/s to m/s/s/s for calculation

      tau0=tsamp/(1.0+ac*tint/2.0/c+ad*tint*tint/6.0/c)

      write(llog,*)'Re-sampling time series. AC=',refac,' m/s/s. AD=',
     &             refad,' cm/s/s/s'
      
      taut=tau0
      next=taut
 10   continue
      if (next.gt.real(i)*tsamp) then
      samp(j)=samp(j)+series(i)*(real(i)*tsamp-(next-taut))/tsamp
      i=i+1
      endif
      if (next.le.real(i)*tsamp) then
         samp(j)=samp(j)+series(i)*(next-real(i-1)*tsamp)/tsamp
         j=j+1
         time=next
         taut=tau0*(1.0+ac*time/c+0.5*ad*time*time/c)
         tav=tav+taut
         nav=nav+1.0
         next=next+taut
      endif
      if (i.lt.ntim) goto 10
      do i=1,ntim
         series(i)=samp(i)
      enddo
      end
c==============================================================================
