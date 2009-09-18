c==============================================================================
      subroutine getfast_rmea(data,ndat,nrun,rafmea,ncal,llog)
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
      integer ndat,nrun,llog,ncal
      real data(*),rafmea(*)

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
c      do h=1,ndat
c         m=m+1
c         do i=1,nrun
          do i=nrun,2*nrun
            if (data(i).ne.0.0) then
               l=l+1
               sum=sum+data(i)
            endif
            j=j+1
            if (j.eq.nrun) then
               k=k+1
               if (l.gt.0) then
                  rafmea(k)=sum/real(l)
c                  rafmea(k)=sum
               else
                  rafmea(k)=0.0
               endif
            endif
         enddo
         
c         i=0
         l=0
c         j=0 ! new line
         
         do i=2,ndat
            if (data(i).ne.0.0) then
               l=l+1
               sum=sum-data(i-1)+data(i+nrun-1)
c     else
c     sum=0.0
            endif
            k=k+1                ! old line 
c            j=j+1                 ! new line
c            if (j.eq.nrun) then ! new line
c               k=k+1            ! new line
               if (l.gt.0) then
                  rafmea(k)=sum/real(nrun)   ! old line
c                  rafmea(k)=sum/real(l)       ! new line
               else
                  rafmea(k)=0.0
               endif
c              j=0              ! new line
c              l=0              ! new line
c          endif               ! new line
         enddo

c               j=0
c               l=0
c     i=i-(k*nrun)+1
c            endif
c         enddo
c      enddo


c
c     Pass back number of calculations
c      
      ncal=k
      end
c
c     That's all folks!
c      
c==============================================================================

