c==============================================================================
      subroutine ralphs_mean(data,ndat,nrun,ramea,ncal,llog)
c==============================================================================
c
c     Calculates the running mean for nrun points of the array data(1:ndat)
c
c     data    - r4  - data array
c     ndat    - i4  - number of data points
c     nrun    - i4  - number of points to average
c      
c     rmea    - r4  - array to store running mean
c     ncal    - i4  - number of means calculated
c
c     Created: November 1997 (dunc@mpifr-bonn.mpg.de)
c
c==============================================================================
c     
      implicit none
      integer ndat,nrun,ncal,llog
      real data(*),ramea(*)
c
c     Local variables
c      
      integer h,i,j,k,l,m
      real sum
c
c     Initialise...
c      
      j=0
      k=0
      l=0
      m=0
      sum=0.0

c
c     Main loop
c      
      do h=1,ndat
c      do h=nrun,ndat
         m=m+1
         do i=1,nrun
            if (data(m+i-1).ne.0.0) then
               l=l+1
               sum=sum+data(m+i-1)
            endif
            j=j+1
            if (j.eq.nrun) then
               k=k+1
               if (l.gt.0) then
                  ramea(k)=sum/real(l)
c                  ramea(k)=sum
               else
                  ramea(k)=0.0
               endif
               sum=0.0
               j=0
               l=0
c     i=i-(k*nrun)+1
            endif
         enddo
      enddo


c
c     Pass back number of calculations
c      
      ncal=k
      end
c
c     That's all folks!
c      
c==============================================================================
