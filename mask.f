      program mask
      implicit none
      include 'vers.inc'
      integer mp,mc,nbins
      parameter(mp=2**25,mc=20,nbins=128)
      integer binval(nbins)
      real fbin(mp),samp(mp),smin,smax,lo,hi,thresh,histmax,hmax
      real dm,ac
      real*8 tsamp,freq
      integer narg,iargc,fold,npf,i,lun
      character*80 comline,title,fname
      logical plotit

      thresh=0.0
      lo=0.0
      hi=0.0
      fold=0
      title='Raw Spectrum'
      plotit=.true.
      lun=-999

      narg=iargc()
      if (narg.gt.0) then
        call getarg(1,fname)
	do i=2,narg
	  call getarg(i,comline)
	  if (comline(1:2).eq.'-t') read(comline(3:),*) thresh
	  if (comline(1:2).eq.'-l') read(comline(3:),*) lo
	  if (comline(1:2).eq.'-h') read(comline(3:),*) hi
	  if (comline(1:2).eq.'-f') plotit=.false.
	enddo
      else
         write(*,*)
         write(*,*) 'MASK: ',version
         write(*,*) 'Produces a spectral mask for use by SEEK'
         write(*,*)
         write(*,*)'usage: mask <SPECTRUM_FILE> -{options}'
         write(*,*)
         write(*,*)'The spectrum file is produced by running find -s'
         write(*,*)'This program is usually run on fold1.spc'
         write(*,*)
         write(*,*)'Available options are:'
         write(*,*)
         write(*,*)'-f: write the spectral mask to file "mask.out"'
         write(*,*)
         write(*,*)'-l[f_hz]: lowest frequency to consider (def=0)'
         write(*,*)'-h[f_hz]: highest frequency to consider (def=Nyqst)'
         write(*,*)'-t[ampl]: amplitude threshold to mask'
         write(*,*)
         write(*,*)'Options given as: mask fold1.spc -t10'
         write(*,*) 
         stop
      endif

      call readspec(fname,fold,samp,dm,ac,tsamp,npf)
      write(*,*) 'Npts:',npf

      smax=-1.0e32
      smin=+1.0e32
      if (lo.eq.0.0) lo=real(freq(tsamp,npf,fold,1))
      if (hi.eq.0.0) hi=real(freq(tsamp,npf,fold,npf))
      do i=1,npf
         fbin(i)=real(freq(tsamp,npf,fold,i))
         if (fbin(i).lt.hi.and.fbin(i).gt.lo) then
            smax=max(smax,samp(i))
            smin=min(smin,samp(i))
         endif
      enddo

      if (thresh.gt.0.0) write(*,*) 'Masking S/N-threshold:',thresh
      write(*,*) 'Fmin:',lo,' Hz'
      write(*,*) 'Fmax:',hi,' Hz'

      if (.not.plotit) then
	call glun(lun)
	open(unit=lun,file='mask.out',status='unknown')
      endif
      do i=1,npf
         if (fbin(i).ge.lo.and.fbin(i).le.hi.and.abs(samp(i)).gt.thresh
     &      .and.thresh.ne.0.0) then
            if (lun.ne.-999) write(lun,*) i,fbin(i),samp(i)
            samp(i)=0.0
c            smax=max(smax,samp(i))
c            smin=min(smin,samp(i))
         endif
      enddo

      if (.not.plotit) then
	write(*,*) '"mask.out" file created'
	stop
      endif

      call pgbegin(0,'?',1,1)
      call pgscf(2)
      call pgsch(1.5)
      call pgvport(0.15,0.85,0.15,0.85)
      call pgadvance

c      call minmax(samp,npf,smin,smax)
      call pgwindow(lo,hi,0.0,smax)
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
      call pgline(npf,fbin,samp)

      call pgend
      stop
      call pgask(.true.)
      call pgadvance
      smin=0.0
c      smax=3.0e5
c      smax=4.0e4
      hmax=histmax(npf,samp,smin,smax,nbins)
      call histval(npf,samp,smin,smax,nbins,binval)
      do i=1,nbins
         write(55,*) i,binval(i)
      enddo
      call pgwindow(smin,smax,0.0,hmax)
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
      call pghist(npf,samp,smin,smax,nbins,1)
      
      call pgend
      end
