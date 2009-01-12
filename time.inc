c==============================================================================
c Definition of the raw time series.
c==============================================================================
c
      include 'npts.inc'
      integer ntim ! actual number of points in time series
      integer ntim2! actual number of points in time series
      integer noff ! offset number of time samples (for fakedata)
      integer iseg ! offset number of segments (for fakedata)
      real series  ! time series
      common /data/ series(npts)
      common /idat/ ntim,ntim2,noff,iseg
c
c==============================================================================

