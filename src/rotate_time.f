C @(#)find_spc.f        3.1 12/17/92
c========================================================================
      Subroutine rotate_time(cdat,nf,nf1)
c========================================================================
      
c  Forms the Raw Power spectrum from the complex FFT .  Takes complex
c     differences? possible that time could be saved by making a version of
C realtr that does this.

c  PAH MXB
c  adapted from - N D'Amico     FIND software V1.3      May 1989
C pah -91/11/11: store the complex spectrum and rotate centre time
C                by negating the alternate complex frequencies
c mxb 91/11/18:  no longer zero the first nf1, let find_norm do this.
c dcs 96/12/17:  remove spectral interpolation lines to test if this
C                gives better sensitivity
c IHS 99/09/14   Comment IN spectral interpolation again
c======================================================================

      implicit none

c      real cdat(*), dat(*)
      real cdat(*)
      integer nf, nf1

      real  arn, ain, ar, ai, aml, amc, amn, tempr, tempi
      integer nfh, j

      nfh=nf-1
      
c  Start forming amplitude spectrum
c  (first zero out unused points)
      
c$$$      do j=1,nf1-1
c$$$         dat(j)=0
c$$$      end do

      ar=cdat(3)
      ai=cdat(4)
      aml=0.0

c for scalar speed.
       do j=2,nf

c rotate the original time series by negating alternate complex no.s  
         if(mod(j,2).eq.1) then
            cdat(2*j+1)=-cdat(2*j+1)
            cdat(2*j+2)=-cdat(2*j+2)
         endif

         arn=cdat(2*j+1)
         ain=cdat(2*j+2)

         tempr = ar+arn
         tempi = ai+ain
         amn=(tempr*tempr + tempi*tempi)*0.5
         amc=(ar*ar+ai*ai) 

C substitute the amplitude of the max of the average of this point and
C the last, this point and the next, and the point itself
C         dat(j-1)=sqrt(amc)
c         dat(j-1)=sqrt(max(aml,amc,amn))
         aml=amn
         ar=arn
         ai=ain
      end do
      
c  complete last block
c      amc=(ar*ar+ai*ai) 
c     dat(nf)=sqrt(max(amc,aml))

c Kill counter-related birdies and Nyquist freq
CIHS 060200 Don't do this anymore - redundant with zap_birdies.f

C     call kill_bc(dat,nf,nf1)

      return
      end
c
