c==============================================================================
      subroutine getrmed(data,ndat,nrun,rmed,ncal)
c==============================================================================
c
c     Calculates the running median for nrun points of the array data(1:ndat)
c
c     data    - r4  - data array
c     ndat    - i4  - number of data points
c     nrun    - i4  - number of points to average
c      
c     rmed    - r4  - array to store running mean
c     ncal    - i4  - number of means calculated
c
c     Created: November 1997 (dunc@mpifr-bonn.mpg.de)
c
c==============================================================================
c     
      implicit none
      integer ndat,nrun,ncal
      real data(*),rmed(*)
c
c     Local variables
c      
      integer i,j,k,nmax
      parameter (nmax=8192)
      real tmp(nmax)
      integer idx(nmax)
      if (nrun.gt.nmax) stop 'getrmed: array size limit reached!'
c
c     Initialise...
c      
      j=0
      k=0
c
c     Main loop
c      
      do i=1,ndat
        j=j+1
        tmp(j)=data(i)
        if (j.eq.nrun) then
  	   call indexxf77(nrun,tmp,idx)
           k=k+1
	   rmed(k)=tmp(idx(nrun/2))
           j=0
        endif
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
      subroutine getrmea(data,ndat,nrun,rmea,ncal)
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
      integer ndat,nrun,ncal
      real data(*),rmea(*)
c
c     Local variables
c      
      integer i,j,k,l
      real sum
c
c     Initialise...
c      
      j=0
      k=0
      l=0
      sum=0.0
c
c     Main loop
c      
      do i=1,ndat
        if (data(i).ne.0.0) then
          l=l+1
          sum=sum+data(i)
        endif
        j=j+1
        if (j.eq.nrun) then
           k=k+1
           if (l.gt.0) then
             rmea(k)=sum/real(l)
           else
             rmea(k)=0.0
           endif
           sum=0.0
           j=0
           l=0
        endif
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
