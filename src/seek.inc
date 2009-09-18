c=============================================================================
c   Common Variable declarations for find (includes common time series/head)
c=============================================================================
c
      include 'vers.inc'
      include 'time.inc'
      include 'head.inc'
      include 'pulse.inc'
      real power   ! power spectrum
      common /four/ power(nfolds,npts/2)
c    
c=============================================================================
    
