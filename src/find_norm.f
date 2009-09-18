C @(#)find_norm.f       3.1 12/17/92
c=======================================================================
      subroutine FIND_Norm(cdat,runave,nf,nf1,llog)
c=======================================================================
      
c     Normalizes the amplitude spectrum  and the complex spectrum 
c     Modified to reduce the amplitudes of high frequencies to give
c     roughly equal numbers of suspects per period interval.
c     AGL Aug 29 2000
      
c======================================================================
      implicit none
      
c      include seek.inc
c      include csamp.inc

c      real fftdat(*)
      real runave(*), cdat(*)

      integer nf,nf1,j,llog
      real g, div, snr, fac(nf/1024)
 
c      write(llog,*)nf, nf1

      snr=7.
      do j=1,nf/1024
c         fac(j)=snr/(snr+0.35*alog10(float(j*1024)/(0.01*nf)))
         fac(j)=1.0
      enddo

c      do j=1,nf1-1
c         fftdat(j)=0.0
c     enddo

c      do j=1,nf1-1    My lines RE. The low frequency bins should be zeroed out? 
c         cdat(j)=0.0
c      enddo


c       write(llog,*)'1000th element of running mean', runave(1000)

c       write(llog,*)'******', nf

c      do j=nf1,nf   Here we can run off the end of the array if there 
c                    are only nsm (1024*8) points in the running averages. RE  


      do j=nf1,nf 
c         g=1
         g=runave(j)
c         write(llog,*) j, nf1, nf, g
         if(j.le.0.01*nf) then
            div=1.0/g
         else
            div=fac((j-1)/1024 + 1)/g
         endif
c         fftdat(j)=fftdat(j)*div-1
         cdat(2*j-1)=cdat(2*j-1)*div
         cdat(2*j)=cdat(2*j)*div
      enddo
      
      return
      end


