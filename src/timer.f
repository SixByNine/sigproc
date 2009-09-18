c
c     These two routines can be called from any two points within a program
c     to start and stop a stop watch based on MJDs from the ship's clock.
c     Note that the logical unit number for the output is required (6=screen)
c
c     Creation date: 98/04/30 (dunc@mpifr-bonn.mpg.de)
c      
c==============================================================================
      subroutine timstart(llog)
c==============================================================================
c      
      integer llog
      double precision mjds,mjdf
      common /timer/ mjds,mjdf
      save
      call getmjd(mjds)
      write(llog,*) 'Timer is up and running...'
      end
c      
c==============================================================================
      subroutine timfinis(llog)
c==============================================================================
c      
      integer llog
      double precision mjds,mjdf
      common /timer/ mjds,mjdf
      save
      call getmjd(mjdf)
      write(llog,*) 'Timer clocked',
     &nint((mjdf-mjds)*86400.0),' s for this job.'
      end
c      
c==============================================================================
