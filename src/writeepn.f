C***********************************************************************
      subroutine writeepn(profile,nbins,ppsr,refdm,sfreq,chbw,
     &tsamp,mjd,tst,epnfile,prfsnr,header)
C***********************************************************************
      implicit none
      include 'epnhdr.inc'
      integer nbins,i
      real profile(nbins),prfsnr,peak,refdm,mean,sfreq,chbw
      character*80 header
      logical first
      data first/.true./
      real*8 ppsr,tsamp,tst,mjd
      save
      real kwmax,smmax
      character*(*) epnfile
c     
c     set up the EPN variables...
c
      history=header
      jname=' '
      cname=' ' 

      pbar=ppsr
      dm=refdm
      rm=0.0

      catref=' '
      bibref=' '

      timflag=' '
      
      rah=0
      ram=0
      ras=0
      ded=0
      dem=0
      des=0

      telname=' '
      xtel=0.0
      ytel=0.0
      ztel=0.0

      epoch=int(mjd)
      if (epoch.gt.0.0) timflag='U'
      opos=0.0
      paflag=' '

c      call idate(cdm,cdd,cdy)
      cdy=cdy+1900   
      cdy=1998
      cdm=01
      cdd=01

      scanno=1
      subscan=1
      npol=1 
      nfreq=1

      nbin=nbins
      tbin=tsamp*1.0e6
      tres=tbin
      fluxflag='U'
      nint=0
      ncal=0
      lcal=0
      fcal=1.0
c
c     Sub-header variables...
c
      do i=1,npol
         rms(i)=0.0
         idfield(i)=epnfile
         nband(i)=i
         navg(i)=1
         f0(i)=sfreq/1000.0
         f0u(i)=' GHz'
         df(i)=chbw
         dfu(i)=' MHz'
         tstart(i) = tst
         if (timflag.eq.'U') tstart(i)=(mjd-epoch)*86400.0*1.0e6
         papp(i)=ppsr
      enddo

      mean=0.0
      do i=1,nbins
         rawdata(1,i)=profile(i)
         mean=mean+profile(i)
      enddo
      mean=mean/real(nbins)
      do i=1,nbins
         rawdata(1,i)=rawdata(1,i)-mean
      enddo

      peak=-1.0e32
      do i=1,nbins
         peak=max(rawdata(1,i),peak)
      enddo

      if (first) then
         recno=-1
         first=.false.
      else
         recno=0          ! just keep writing to the same file...
      endif
c      write(*,'('' Folded profile written to '',a)') epnfile
      filename=epnfile
      readwri=1        ! writing required
      padout=.false.   ! don't pad out the data 
c
c     Now call the interface routine...
c
      call rwepn(filename, readwri, recno, padout)
c
c     Calculate PRFSNR
c
c
      prfsnr=0.0
      if (peak.gt.0.0) prfsnr=peak/rms(1)
      return
      do i=1,nbins
         profile(i)=rawdata(1,i)
      enddo
      call smooth(profile,nbins,kwmax,prfsnr,smmax)
      return
      end


