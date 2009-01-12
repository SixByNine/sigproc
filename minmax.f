      subroutine minmax(data,n,dmin,dmax)
      integer i,n
      real data(n),dmin,dmax
      dmin=+1.0e32
      dmax=-1.0e32
      do i=1,n
      dmin=min(dmin,data(i))
      dmax=max(dmax,data(i))
      enddo
      end
