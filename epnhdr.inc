c==============================================================================
c
c     Include file "epnhdr.inc"
c
c     This file defines the variables used in the EPN-header
c
c==============================================================================
c		Main Header   - 6 Lines 80 characters per line
c==============================================================================
c     version     - A8     - version number of EPN format.
c     counter     - I4     - the number of 80 character records in this block
c     history     - A72    - some (preferably meaningful) notes about the data
c------------------------------------------------------------------------------
c     jname       - A12    - Name of pulsar derived from its J2000 coordinates
c     cname       - A12    - Common name of pulsar (either B- or J-name)
c     pbar        - F16.11 - Barycentric period of pulsar (s) (present epoch)
c     dm          - F8.3   - Dispersion measure (pc/cc)
c     rm          - F10.3  - Rotation measure of pulsar (rad m^2)
c     catref      - A6     - Tag showing which catalogue in use 
c     bibref      - A8     - Bibliographic reference for data 
c------------------------------------------------------------------------------
c     rah         - I2     - Hours of right ascension (J2000)
c     ram         - I2     - Mins. of right ascension (J2000)
c     ras         - F6.3   - Secs. of right ascension (J2000)
c     ded         - I3     - Degrees of declination   (J2000)
c     dem         - I2     - Minutes of declination   (J2000)
c     des         - F6.3   - Seconds of declination   (J2000)
c     telname     - A8     - Telescope name
c     epoch       - F10.3  - modified Julian date of observation
c     opos        - F8.3   - position angle of feed (degrees)
c     paflag      - A1     - "A" signifies absolute PAs else undefined
c     timflag     - A1     - "U" signifies UTC time "B" barycentric UTC 
c------------------------------------------------------------------------------
c     xtel        - F17.5  - Topocentric rectangular position of telescope [m]
c     ytel        - F17.5  - Topocentric rectangular position of telescope [m]
c     ztel        - F17.5  - Topocentric rectangular position of telescope [m]
c------------------------------------------------------------------------------
c     cdy         - I4     - year (eg 1996)
c     cdm         - I2     - month (1-12)
c     cdd         - I2     - day (1-31)
c     scanno      - I4     - sequence number of the observation
c     subscan     - I4     - sub-sequence number of the observation
c     npol        - I2     - number of polarisations observed
c     nfreq       - I4     - number of frequency bands per polarisation
c     nbin        - I4     - number of phase bins per frequency
c     tbin        - F12.6  - sampling interval (us)
c     nint        - I6     - number of integrated pulses per data block
c     ncal        - I4     - bin number of start of cal signal
c     lcal        - I4     - number of bins in the cal signal
c     tres        - F12.6  - temporal resolution (us)
c     fluxflag    - A1     - "F" signifies flux (mJy) calibrated data
c==============================================================================
c     Block sub-header.... 2 Lines 80 Characters per Line
c==============================================================================
c     idfield     - A8     - Description of data stream I Q U V etc. etc.
c     nband       - I4     - ordinal number of data stream
c     navg        - I4     - number of streams averaged into current one
c     f0          - f12.8  - Effective centre sky frequency of this stream
c     f0u         - A8     - String giving unit of f0 [default is GHz]
c     df          - f12.6  - Effective bandwidth of this stream
c     dfu         - A8     - String giving unit of df [default is MHz]
c     tstart      - F17.5  - Time of first bin wrt EPOCH [us]
c------------------------------------------------------------------------------
c     scale       - E12.6  - Scale factor for the data
c     offset      - E12.6  - Offset to be added to the data
c     rms         - E12.6  - rms of the data stream
c     papp        - F16.12 - Apparent period at time of first phase bin [s]
c==============================================================================
c
c     EPN Filename
c
      character*80 filename
c
c     First line of main header
c
      character version*8, history*72
      integer*4 counter
      common/epn1/version,history,counter
c      
c     Second line of main header
c
      character jname*12, cname*12
      real*8    pbar, dm, rm
      character catref*6
      character bibref*8
      common/epn2/jname,cname,pbar,dm,rm,catref,bibref
c
c     Third line of main header
c
      integer*4 rah, ram, ded, dem
      real      ras, des
      character telname*8
      real*8    epoch , opos
      character paflag*1, timflag*1
      common/epn3/rah,ram,ded,dem,ras,des,telname,epoch,opos,
     &             paflag,timflag
c
c     Fourth line of main header
c
      real*8 xtel, ytel, ztel
      common/epn4/xtel,ytel,ztel
c
c     Fifth line of main header
c
      integer*4 cdy, cdm, cdd
      integer*4 scanno, subscan, npol, nfreq, nbin
      real*8    tbin, tres
      integer*4 nint ,ncal, lcal
      real      fcal
      character fluxflag*1
      common/epn5/cdy,cdm,cdd,scanno,subscan,npol,nfreq,nbin,tbin,
     &             nint,ncal,lcal,fcal,tres,fluxflag
c
c     Maximum blocks writing to sub-header...
c
      integer maxblk
      parameter(maxblk=8)
c
c     First line of block sub-header 
c
      character idfield*8
      integer*4 nband, navg
      real*8    f0, df
      character*8 f0u,dfu
      double precision tstart
      common/epns1/idfield(maxblk),nband(maxblk),navg(maxblk),
     & f0(maxblk),f0u(maxblk),df(maxblk),dfu(maxblk),tstart(maxblk)
c
c     Second line of block sub-header (required to scale/descale data)
c
      real*8    scale,offset,rms,papp
      common/epns2/offset(maxblk),scale(maxblk),rms(maxblk),papp(maxblk)
c
c     Data block 
c
      integer maxbin
      parameter(maxbin=16384)
      real*4 rawdata
      common/dblk/rawdata(maxblk,maxbin)
c
c     Reading or writing?
c
      integer readwri     ! = 1 write = -1 read
c
c     Record number for reading... and record length
c
      integer recno
c
c     Pad out the record to its maximum length
c
      logical padout
c
c     That's all folks!
c
c==============================================================================
