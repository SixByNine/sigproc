c==============================================================================
      subroutine getmeanrms(data,ndat,nrun,rmea,ncal,rrms)
c==============================================================================
c
c     Calculates the running mean and rms for nrun points of the array data(1:ndat)
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
      integer ndat,nrun,ncal
      real data(*),rmea(*), rrms(*)
      real*4 mean,rms,block(nrun)
c
c     Local variables
c      
      integer i,j,k,l,nrej
      integer m,n

c
c     Initialise...
c      
      j=0
      k=0
      l=0
      m=0
c      n=1
      nrej=0

c
c     Main loop
c      
      do i=1,ndat
         if (data(i).ne.0.0) then
            l=l+1
c            m=m+1
            block(l)=data(i)    
         endif
         j=j+1
         if (j.eq.nrun) then
            k=k+1
            if (l.gt.0) then
               call meanrms(block, l, mean, rms, nrej) ! may want to pass l instead of nrun
               rmea(k)= mean
               rrms(k) = rms
            else
               rmea(k)=0.0
               rrms(k) = 0.0
            endif
c     WRITE(90,*) rmea(k)  !, rrms(k)
            
            j=0
            l=0
c            m=0
            
c            do n=1,nrun
c               block(n)=0.0
c            enddo
c     WRITE(69,*) block(n)
c            n=0
         endif
      enddo
c     
c     Pass back number of calculations
c      
      ncal=k
      end
c
c      
c==============================================================================
