      real function spcsnr(dat,n)
      implicit none
      integer n
      real dat(n)
      integer i,j,nb
      real rms,peak,sumsq
      peak=-1.0e32
      sumsq=0.0
      nb=n/10
      j=0
      do i=1,n
         if (i.le.nb.or.i.ge.n-nb) then
            sumsq=sumsq+dat(i)*dat(i)
            j=j+1
         else
            peak=max(peak,dat(i))
         endif
      enddo
      rms=sqrt(sumsq/real(j))
      spcsnr=0.0
      if (peak.gt.0.0) spcsnr=peak/rms
      end
