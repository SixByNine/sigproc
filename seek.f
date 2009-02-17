c=============================================================================
      program seek 
c=============================================================================
c
c     A program (formerly known as find) to seek periodic signals in data
c
c     Created: 97/11/21 (dunc@mpifr-bonn.mpg.de)
c
c     Modification history:
c
c     98/04/10 - (dunc) added acceleration capability (no longer need rbin)
c     98/04/30 - (dunc) source code overhauled. Now more user-friendly
c     98/11/20 - (dunc) added read and FFT capability for ".dis" files
c     99/07/08 - (dunc) added read spectrum file capability
c     01/02/15 - (dunc) added capability to read new data format ".tim" files
c     01/10/11 - (dunc) added -pmzap option and fixed call to zapit routine
c     02/03/01 - (dunc) added -pulse option to call Maura's single-pulse code
c     02/03/20 - (dunc) changed ordering of spectral zapping in dosearch.f
c     02/03/21 - (dunc) added -pzero option to pad with zeros if need be
c     02/03/21 - (dunc) make pmax a command-line parameter
c     05/04/07 - (dunc) changed name of program to SEEK! to appease RNM et al.
c     05/04/28 - (dunc) added ability to read masks for all 5 harmonic folds
c     06/04/05 - (dunc) added -mmzap option for Methanol multibeam survey
c     07/01/26 - (dunc) added Ralph Eatough's fftw wrapper using -fftw option
c     07/02/01 - (dunc) added Ralph Eatough's profile recon mods using -recon
c     07/02/22 - (dunc) added header output to .prd file via -head option
c     07/06/13 - (dunc) added -fftdump option to write out .fft files
c     07/03/13 - (dunc) can now read PRESTO .dat/inf time series format
c
c=============================================================================
      implicit none
      include 'seek.inc'
      logical dump,rspc,acsearch,tanalyse,pmzap,mmzap,pulse,
     &        append,pzero,fftw,recon,prdh
      integer oldw ! = 0 new mean, = 1 old mean, = 2 median
      character*80 sfile
      real accn,adot
      real*8 pmax
      integer llog       
      call seekin(llog,dump,rspc,pmzap,mmzap,sfile,pulse,append,pzero,
     &	  pmax,nofft,fftw,recon,oldw,prdh,spthresh,ncandsmax,nsmax)
      accn=refac
      adot=refad
      call timstart(llog)                    ! fire up the ship's clock
      if (.not.rspc) call readdat(llog,pzero)! read in the time series 
      if (pulse) then
        call baseline(llog)
        call singlepulse(llog,append,spthresh,ncandsmax,nsmax)
      endif
      if (.not.nofft) then
      if (accn.ne.0.0) refac=accn
      if (adot.ne.0.0) refad=adot
      acsearch=accn.ne.0.0.or.adot.ne.0.0
      tanalyse=(index(filename,'.ser').gt.0.0)
     &     .or.(index(filename,'.tim').gt.0.0)
     &     .or.(index(filename,'.dat').gt.0.0)
     &     .or.(index(filename,'.dis').gt.0.0)
      if (rspc) tanalyse=.false.
      if (tanalyse) then                     ! (time series analysis only)
         if (acsearch) call resample(llog)   ! re-sample time series
	 if (fftw) then
         	call fftwdata(llog)          ! fft the data using FFTW
	 else 
         	call fftdata(llog)           ! fft the data using SINGLETON
         endif
         if (dumpfft) then
            call dumpdat(filename(1:index(filename,'.tim')-1)//'.fft')
            stop
         endif
      endif                                  ! (standard analysis follows)
      call dosearch(llog,dump,rspc,oldw,pmzap,
     &           mmzap,recon,prdh,sfile,pmax)! search 
      endif
      call timfinis(llog)                    ! stop the clock
      end
c=============================================================================
