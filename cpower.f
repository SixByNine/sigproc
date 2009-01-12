c==============================================================================
      real function cpower(npf,nharm,pcand)
c==============================================================================
c
c
c==============================================================================
      implicit none
      include 'seek.inc' 
      include 'csamp.inc'
      integer npf,nharm
      real*8 pcand,df0
      real re,im
      integer ifun,i
c
c     Initialization
c
      cpower=0.0        ! coherent power TBD
      re=0.0            ! real part of coherent sum TBD
      im=0.0            ! imaginary part of coherent sum TBD
      df0=2.0*real(npf)*tsamp/(pcand/1000.0) !position of fundamental
c
c     Do the vector addition and return the coherent power
c
      do i=1,nharm
         if (df0*i.lt.npf) then
            ifun=nint(df0*i)*2+1
            re=re+series(ifun)
            im=im+series(ifun+1)
         endif
      enddo
      cpower=re*re+im*im
      end
c==============================================================================
