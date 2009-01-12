c=============================================================================
      subroutine fftwdata(llog)
c=============================================================================
c
c     Dummy wrapper which gets compiles if FFTW library absent
c     - DRL Jan 26, 2007
c    
c     llog    - i4  - llogical unit number for all but warning messages.
c      
c=============================================================================
c
      integer llog
      write(llog,*) 'FFTW not included in the compile! To use this',
     & ' feature, do the following:'
      write(llog,*) '1: go to your sigproc source code directory'
      write(llog,*) '2: type "make clean"'
      write(llog,*) '3: edit "makefile.$OSTYPE so that LFFTW is',
     & ' defined correctly'
      write(llog,*) '4: type "make seek"'
      write(llog,*) 'N.B. You will need to ensure that FFTW is',
     & ' installed on your computer!'
      write(llog,*)
      stop
      end
c      
c=============================================================================
