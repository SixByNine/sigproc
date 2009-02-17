c==============================================================================
      subroutine clock(yy,mm,dd,hh,mi,ss)
c==============================================================================
c
c     Returns the present time from the Sun-OS clock
c      
c     yy - i4 - year   
c     mm - i4 - month
c     dd - i4 - day  
c     hh - i4 - hour
c     mi - i4 - minute
c     ss - i4 - second
c      
c     Creation date: 98/06/08 (dlorimer@naic.edu)
c      
c==============================================================================
c      
      implicit none
      integer yy,mm,dd,hh,mi,ss
      
      integer iarray(3)
      
      call idate(iarray)
      mm=iarray(2)
      dd=iarray(1)
      yy=iarray(3)
      call itime(iarray)
      hh=iarray(1)
      mi=iarray(2)
      ss=iarray(3)
      end
c
c==============================================================================

