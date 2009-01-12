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
c     MKeith 2007/05/09: Updated to DATE_AND_TIME as this is the standard      
c     Creation date: 98/06/08 (dlorimer@naic.edu)
c      
c==============================================================================
c      
      implicit none
      integer yy,mm,dd,hh,mi,ss
      
      integer iarray(8)
      character (LEN = 10) tmp(3)      
c      call idate(iarray)
c      mm=iarray(2)
c      dd=iarray(1)
c      yy=iarray(3)
c      call itime(iarray)
c      hh=iarray(1)
c      mi=iarray(2)
c      ss=iarray(3)

       call DATE_AND_TIME(tmp(1),tmp(2),tmp(3),iarray)
       yy=iarray(1)
       mm=iarray(2)
       dd=iarray(3)

       hh=iarray(5)
       mi=iarray(6)
       ss=iarray(7)
      end
c
c==============================================================================

