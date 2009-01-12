c==============================================================================
      subroutine sprof(profile,nbins,shift)
c==============================================================================
c
c     shifts the profile in the array profile() with nbins bins by the 
c     number of bins passed down in shift. shift>0 means shift forward,
c     shift<0 means shift backwards.
c
      implicit none

      real profile(*)
      integer nbins, shift
c
c     local variables
c
      integer i,j
      real dummy

      if (shift.lt.0) shift=nbins+shift
      do i=1,shift
         dummy=profile(nbins)
         do j=nbins,2,-1
            profile(j)=profile(j-1)
         enddo
         profile(1)=dummy
      enddo

      end
