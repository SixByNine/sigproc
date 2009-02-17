c==============================================================================
c Definition of the header variables in the time series files.
c==============================================================================
c
      real*8 tsamp   ! sampling time (seconds)
      common /head/ tsamp
      integer nchans,channo,maxchans  ! number of filterbank channels
      parameter(maxchans=4)
c     addition to specify number of folds - dbarnes 20060329
      include 'folds.inc'
      real skyfreq,fref
      common /fbank/ nchans,skyfreq(maxchans),channo,fref
      real refdm     ! reference DM for time series
      common /dmref/ refdm
      real refac     ! reference acceleration for time series (m/s/s)
      real refad     ! reference adot for time series (cm/s/s/s)
      common /accn/ refac,refad
      character*80 filename,header,facc
      integer lst
      common /datafile/ filename,header,lst,facc
      integer dmidx
      common /dmindex/ dmidx
      character*80 zapfile,maskfile(nfolds)
      integer nmasks
      common /zapfname/ zapfile,maskfile,nmasks
      real fbrute
      common /brutal/ fbrute
      integer rfac
      common /rebinfac/ rfac
      real    skp
      common /skipfac/ skp
      integer tsize
      character*12 ra,dec
      common /coords/ ra,dec
      common /transform/ tsize
      logical dumpraw
      common /rawspec/ dumpraw
      logical dumpfft
      common /fftdump/ dumpfft
      logical noise
      common /test/ noise
      double precision mjdstart
      common /timestamp/ mjdstart
      character*80 source_name,telname
      common /names/ source_name,telname
c
c==============================================================================
