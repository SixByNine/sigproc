c==========================================================================
      subroutine formspec(npf,nf1)
c==========================================================================
      implicit none
c
c     Forms the power spectrum from the real and imaginary
c     parts ffted array series(ntim) (see seek.inc)
c
      include 'seek.inc'
      include 'csamp.inc'
      integer npf,nf1
      integer i,j
      real arl,ail,a1,a2,ar,ai,anf

      anf=series(2)**2
      arl=0.0
      ail=0.0
      do i=1,2*nf1
         series(i)=0.0
      enddo

      npf=ntim/2
      do j=1,npf-1
        ar=series(2*j+1)
        ai=series(2*j+2)
        a1=ar**2+ai**2
        a2=((ar-arl)**2+(ai-ail)**2)/2.
        samp(j)=sqrt(max(a1,a2))
        arl=ar
        ail=ai
      enddo
      samp(npf)=0.0
c      ar=sqrt(anf)
c      a1=anf
c      a2=((ar-arl)**2 + ail**2)/2.
c      samp(npf)=sqrt(max(a1,a2))
      end
c==========================================================================
