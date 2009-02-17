c=============================================================================
      subroutine seekin(llog,dump,rspc,pmzap,mmzap,sfile,pulse,app,
     & pzero,pmax,nofft,fftw,recon,oldw,prdh,spthresh,ncandsmax,nsmax)
c=============================================================================
c
c   Controls the command-line inputs to seek
c
      include 'seek.inc'
      integer narg,iargc,i,llog,p2,lun
      character*80 option,sfile
      logical dump,rspc,pmzap,mmzap,pulse,app,pzero,fftw,recon,prdh
      integer oldw
      real*8 pmax

      p2=nint(log10(real(npts))/log10(2.0))
      narg=iargc()
      call getarg(1,option)
      if (option.eq.'version'.or.option.eq.'-version') then
         write(*,'(a,a)') 'PROGRAM: seek ',version
         stop
      endif
      if (option.eq.'-hsums') then
      write(*,759) 'Harmonic sums are: '
      do i=1,nfolds
         write(*,762) foldvals(i)
      enddo
 759  format(a,$)
 762  format(1x,i4,2x,$)
	write(*,*)
         stop
      endif
      if (option.eq.'-maxft') then
         write(*,2) 'Maximum length of time series 2^',p2,' pts'
         stop
      endif
      if (narg.lt.1.or.option.eq.'-help'.or.option.eq.'help') then
         write(*,*)
         write(*,1) 'seek - searches for periodic and/or'//
     &              ' transient signals in a noisy time series'
         write(*,*)
         write(*,1)'usage: seek <INFILE> -{options}'
         write(*,*)
         write(*,1)'The input file may be a time series, or a set of'//
     &             'Fourier coefficients.'
         write(*,1)'The file extension MUST, however, be either'//
     &             ' ".tim" ".ser" ".dat" ".dis" or ".fft"'
         write(*,1)'In the latter case, the FFT stage is skipped.'
         write(*,*)
         write(*,1)'options:'
         write(*,*)
         write(*,1)'-A       - append output ASCII files'
         write(*,1)'-s       - dump spectra to ".spc" files'
         write(*,1)'-q       - quiet mode - all messages > INFILE.log'
         write(*,1)'-pmzap   - calls the PM survey routine pm_zapbrd'
         write(*,1)'-mmzap   - calls the PM survey routine mm_zapbrd'
         write(*,1)'-pulse   - calls Maura single pulse search'
         write(*,1)'-pzero   - pads out data with zeros (def=Gaussian)'
         write(*,1)'-maxft   - report max length of Fourier transform'
	 write(*,1)'-nofft   - turns off FFT search'
         write(*,1)'-dumpfft - dump Fourier coefficients to a .fft file'
	 write(*,1)'-fftw    - uses FFTW instead of SINGLETON routine'
	 write(*,1)'-hsums   - show harmonic summing used'
	 write(*,1)'-recon   - calculate and report reconstructed S/N'
	 write(*,1)'-submn   - mean subtraction (old method) to whiten spectrum'
	 write(*,1)'-submd   - median subtraction to whiten spectrum'
	 write(*,1)'-head    - adds header info to output .prd files'
         write(*,*)
         write(*,1)'-m[file] - mask birdies from file (def="mask")'
         write(*,1)'-z[file] - zap birdies from file (def="birdies")'
	 write(*,1)'-b[freq] - zap 10-sig+ spikes < freq (def=100 Hz)'
         write(*,1)'-c[cfac] - add every cfac samples before FFT'
         write(*,1)'-a[accn] - re-sample at constant accn (m/s/s)'
         write(*,1)'-d[adot] - re-sample at constant adot (cm/s/s/s)'
         write(*,1)'-D[dmvl] - change header DM to be dmvl'
         write(*,1)'-t[tlen] - fix transform length to 2**tlen'
         write(*,1)'-i[tsec] - ignore tsec seconds of data on reading'
         write(*,1)'-p[pmax] - set maximum period of seach (def=9.999s)'
	 write(*,1)'-T[spth] - set single-pulse search threshold (sigma)'
	 write(*,1)'-n[nmax] - maximum number of single-pulse candidates per DM channel'
	 write(*,1)'-w[smax] - number of times to smooth time series for single-pulse search'
         write(*,*)
 1       format(a)
 2       format(a,i2,a)
         stop
      endif

      pmax=9.99999999
      dmidx=-1
      spthresh=4
      ncandsmax=100000
      nsmax=10
      fbrute=0.0
      zapfile=' '
      maskfile(1)=' '
      nmasks=0
      rfac=1
      llog=6
      tsize=0
      skp=0.0
      dump=.false.
      pmzap=.false.
      mmzap=.false.
      pulse=.false.
      nofft=.false.
      dumpfft=.false.
      pzero=.false.
      fftw=.false.
      prdh=.false.
      dumpraw=.false. 
      app=.false.
      recon=.false.
      oldw= 0 
      facc='sequential'
      refac=0.0
      rspc=.false.
      sfile='fold0.spc'
      refdm=0.0 
      call getarg(1,filename)
      
      if (index(filename,'.ser').gt.0) then
         lst=index(filename,'.ser')-1
      else if (index(filename,'.tim').gt.0) then
         lst=index(filename,'.tim')-1
      else if (index(filename,'.dis').gt.0) then
         lst=index(filename,'.dis')-1
         dmidx=1
      else if (index(filename,'.fft').gt.0) then
         lst=index(filename,'.fft')-1
      else if (index(filename,'.dat').gt.0) then
         lst=index(filename,'.dat')-1
      else if (index(filename,'.spc').gt.0) then
         lst=index(filename,'.spc')-1
         rspc=.true.
         sfile=filename
      else
         stop 'file type not recognized! Type seek for help.'
      endif
      
      do i=2,narg
        call getarg(i,option)
        if (index(option,'-z').gt.0) then
          zapfile='birdies' ! default file
	  if (option(3:).ne.' ') zapfile=option(3:)
        else if (index(option,'-mmzap').gt.0) then
           mmzap=.true.
        else if (index(option,'-m').gt.0) then
          maskfile(1)='mask' ! default file
	  if (option(3:).ne.' ') maskfile(1)=option(3:)
          call glun(lun)
          open(unit=lun,file=maskfile(1),status='old')
          read(lun,'(a)') option
c
c         establish whether 1 or 5 mask files are being read
c
          if (option.eq.'maskfiles are as follows:') then
             do f=1,5
                read(lun,'(a)') maskfile(f)
             enddo
             nmasks=5 ! added option to mask all five folds drl 28/05/05
          else
             nmasks=1 ! this is the old method of masking unfolded spectrum
          endif
          close(unit=lun)
        else if (index(option,'-b').gt.0) then
	  fbrute=100.0 ! default value
          if (option(3:).ne.' ') read(option(3:),*) fbrute
        else if (index(option,'-i').gt.0) then
c          read(option(3:),*) dmidx
          read(option(3:),*) skp
        else if (index(option,'-D').gt.0) then
          read(option(3:),*) refdm
        else if (index(option,'-a').gt.0) then
          read(option(3:),*) refac
        else if (index(option,'-d').gt.0) then
          read(option(3:),*) refad
        else if (index(option,'-c').gt.0) then
          read(option(3:),*) rfac
c          if (mod(rfac,2).ne.0) then
c            write(*,*) 'WARNING - rebin factor must be a multiple of 2!'
c            write(*,*) '*** input value',rfac,' ignored...'
c            rfac=1
c          endif
        else if (index(option,'-t').gt.0) then
           read(option(3:),*) tsize
        else if (index(option,'-q').gt.0) then
           call glun(llog)
        else if (index(option,'-recon').gt.0) then
           recon=.true.
        else if (index(option,'-submn').gt.0) then
           oldw=1
        else if (index(option,'-submd').gt.0) then
           oldw=2
        else if (index(option,'-pmzap').gt.0) then
           pmzap=.true.
        else if (index(option,'-pulse').gt.0) then
           pulse=.true.
	else if (index(option,'-nofft').gt.0) then
	   nofft=.true.
	else if (index(option,'-fftdump').gt.0) then
c	   dumpfft=.true.
	else if (index(option,'-fftw').gt.0) then
	   fftw=.true.
	else if (index(option,'-head').gt.0) then
	   prdh=.true.
        else if (index(option,'-pzero').gt.0) then
           pzero=.true.
        else if (index(option,'-p').gt.0) then
          read(option(3:),*) pmax
	else if (index(option,'-T').gt.0) then
          read(option(3:),*) spthresh
	else if (index(option,'-w').gt.0) then
	  read(option(3:),*) nsmax
	else if (index(option,'-n').gt.0) then
	  read(option(3:),*) ncandsmax
        else if (index(option,'-s').gt.0) then
          dump=.true.
	  if (option(3:).ne.' ') dumpraw=.true.
	  sfile=' '
	  if (option(3:).ne.' ') sfile=option(3:)
        else if (index(option,'-A').gt.0) then
          facc='append'
          app=.true.
        else
          write(*,*) 'WARNING.. command line option ',
     &    option(1:index(option,' ')-1),' not recognized!!!'
        endif
      enddo
      
      if (llog.ne.6)
     &open(unit=llog,file=filename(1:lst)//'.log',status='unknown',
     &     access=facc)
      write(llog,*) 'SEEK: ',version
      end
