      subroutine rebin(series,ntim,tsamp,rfac)
      implicit none
      real series(*),tsamp
      integer ntim,rfac
      integer i,j,k
      real sum
      
      j=0
      k=0
      sum=0.0
      do i=1,ntim
         sum=sum+series(i)
         j=j+1
         if (j.eq.rfac) then
            k=k+1
            series(k)=sum/real(rfac)
            sum=0.0
            j=0
         endif
      enddo
      ntim=k
      tsamp=tsamp*real(rfac)
      end
      
