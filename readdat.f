c=============================================================================
	subroutine readdat(llog,pzero)
c=============================================================================
        implicit none 
	include 'seek.inc'
	real rtsamp,tmin,temp,sum,sumsq,mean,sigma,sample,tmpdm,firstone
	integer lun,i,j,llog,near2,imin,isec,nskip,isamp,ichan,nadd
	integer length,nsrec,nchan,nread,seed,norg,nmax,ihour
 	integer readsample,skipsample
	character*80 pfile,presto
	real f0,chbw
	logical filterbank,opened,pzero
	data opened/.false./
	save

	if (index(filename,'.tim').gt.0.or.
     &     index(filename,'.dat').gt.0) then
	   header=' '
	   refac=0.0
	   nchans=1
	   if (index(filename,'.tim').gt.0) then
	      write(llog,*) 'Working with SIGPROC time series data...'
	      call readhd(filename,tmpdm,rtsamp,ra,dec,fref,
     &         mjdstart,source_name,telname)
	   else
	      write(llog,*) 'Working with PRESTO time series data...'
	      call glun(lun)
	      open(unit=lun,
     &        file=filename(1:index(filename,'.dat')-1)//'.inf',
     &        status='old')
	      presto=''
	      refac=0
	      do while (presto(1:15).ne.' Number of bins')
		 read(lun,'(a)') presto
	      enddo
	      read(presto(index(presto,'=')+1:),*) ntim
	      
	      rewind(lun)
	      do while (presto(1:19).ne.' Width of each time')
		 read(lun,'(a)') presto
	      enddo
	      read(presto(index(presto,'=')+1:),*) rtsamp
	      
	      rewind(lun)
	      do while (presto(1:11).ne.' Dispersion')
		 read(lun,'(a)') presto
	      enddo
	      read(presto(index(presto,'=')+1:),*) refdm
	      close(lun)
	      
	      call openpresto(filename)
	   endif

	   if (refdm.eq.0.0) refdm=tmpdm
	   ntim=1
	   nread=1
	   nadd=0
	   sum=0.0
	   nskip=skp/rtsamp
	   if (nskip.lt.0) stop "silly skip length passed down..."
	   if (nskip.gt.0) then
		write(*,*) 'Skipping',nskip,' samples...'
	        i=skipsample(nskip)
		if (i.ne.0) stop "Error skipping..."
		nskip=0
 	   endif
	   if (tsize.eq.0) then
	      nmax=npts
	   else
	      nmax=2**tsize
	   endif
	   j=0
	   do while (ntim.lt.nmax.and.nread.ne.0) 
	      nread=readsample(sample)
	      j=j+1
	      if (j.eq.1) firstone=sample
	      sample=sample-firstone
	      if (j.gt.nskip) then
		 sum=sum+sample
		 nadd=nadd+1
		 if (nadd.eq.rfac) then
		    series(ntim)=sum/real(rfac)
		    ntim=ntim+1
		    sum=0.0
		    nadd=0
		 endif
	      endif
	   enddo
	   rtsamp=rtsamp*real(rfac)
	   write(llog,*) 'Read',ntim,' samples...'
c	   write(llog,'(x,a80)') header
	   write(llog,*) 'Reference DM:',refdm
	   write(llog,*) 'Sampling time:',rtsamp*1.0e6,' us'
	   if (rfac.gt.1)write(llog,*) '(every',rfac,
     &     ' samples were added during reading of time series)'
	   tsamp=rtsamp
	   ntim2 = ntim
	   goto 2
	endif

	if (dmidx.ne.-1) then
	   call rdhead(filename(1:lst)//'.hdr',pfile,f0,chbw,
     &                 rtsamp,nsrec,nchan)
	   refdm=real(dmidx-1)*(32.0/real(nchan-1))*(rtsamp/0.08)
	   rtsamp=rtsamp/1000.0
	   tsamp=rtsamp
	   call getddis(llog,filename,nchan,dmidx)
	   write(llog,*) 'Reference DM:',refdm,' Reference AC:',refac
	   write(llog,*) 'Sampling time:',rtsamp*1.0e6,' us'
	   nchans=1
	   goto 2
	endif

	filterbank=.false.
	if (filterbank) then
	   if (.not.opened) call rdfbtab(skyfreq,maxchans,nchans)
	   write(llog,*) 'Working with filterbank data...'
	else if (index(filename,'.fft').gt.0) then
	   write(llog,*) 'Working with pre-FFT-ed data...'
        else
	   write(llog,*) 'Working with time series data...'
	endif

	if (opened) goto 1
	call glun(lun)
	open(unit=lun,file=filename,status='old',!err=999)
     &  form='unformatted',err=999)
	opened=.true.
	ichan=0
        write(llog,'('' Opened input data file: '',a)')
     &  filename(1:lst+4)
	header=' '
	header=header(1:length(header))
	
	if (filterbank) then
	   read(lun) ntim,nchans,isamp
	   rtsamp=real(isamp)*1.0e-6
	   refdm=0.0
	   refac=0.0
	else
	   read(lun) ntim,rtsamp,refdm,refac
	   nchans=1
	endif
	if (ntim.gt.npts) then
	   write(*,*) 'WARNING - too many points in time series!'
	   write(*,*) '** reading in max possible:',npts/1024,' kpts'
	   ntim=npts
	endif
	
 	write(llog,*) 'Reading',ntim/1024,' kpt data file...',ntim 
	write(llog,*) 'Reference DM:',refdm
	write(llog,*) 'Sampling time:',rtsamp*1.0e6,' us'
	tsamp=rtsamp

 1	read(lun) (series(i),i=1,ntim)

	ichan=ichan+1
	if (filterbank) write(llog,*) 'Sky frequency of this channel:',
     &  skyfreq(ichan),' MHz'
	
	if (nchans.eq.1) then
          close(unit=lun)
	  opened=.false.
	endif

 2	continue
c 	if (rfac.gt.1) then
c	   write(llog,*) 'Compressing time series by a factor of',rfac
c	   call rebin(series,ntim,rtsamp,rfac)
c           tsamp=rtsamp
c	   write(llog,*) 'Data file now contains',ntim/1024,' kpts'
c  	   write(llog,*) 'Sampling time now:',rtsamp*1.0e6,' us'
c	endif

c	nskip=0
c	if (skp.gt.0) nskip=skp*rtsamp
c	if (nskip.gt.ntim) then
c	   write(*,*) 'WARNING - requested skip size too large!'
c	   write(*,*) '** no skipping carried out...'
c	   nskip=0
c	endif
c	
c	if (nskip.gt.0) then
c	   j=0
c	   do i=nskip+1,ntim
c	      j=j+1
c	      series(j)=series(i)
c	   enddo
c	   ntim=j
c	   write(llog,*) 'Skipped ',nskip,' points...'
c        endif

	norg=ntim
	temp=log10(real(ntim))/log10(2.0)
	near2=nint(temp)

	if (tsize.eq.0) then
	   write(llog,*) 'Nearest power of 2:',near2
	   tsize=near2
	   ntim=2**tsize
	else
	   write(llog,*) 'Requested transform power of 2:',tsize
	   ntim=2**tsize
	   if (ntim.gt.npts) then
	     write(*,*) 'WARNING - requested transform length too large!'
	     write(*,*) '** using max possible:',npts/1024,' kpts'
	     ntim=npts
	   endif
	endif
	
c
c	pad out rest of time series with either gaussian noise or zeros
c
	if (ntim.gt.norg) then
	sum=0.0
	sumsq=0.0
	do i=1,256
	   sum=sum+series(i)
	   sumsq=sumsq+series(i)*series(i)
	enddo	
	mean=sum/256.0
	sigma=sqrt(sumsq)/16.0
	   write(llog,*) 'Padding time series with additional zeros...'
	   do i=norg+1,ntim
	      series(i)=0.0
	   enddo
	endif

        tmin=rtsamp*real(ntim)/60.0
        imin=int(tmin)
        isec=60.0*(tmin-real(imin))
	if (imin.lt.60) then
	   write(llog,*) 'Data length:',imin,' min',isec,' sec'
	else
	   ihour=imin/60
	   imin=imin-(ihour*60)
	   write(llog,*) 'Data length:',ihour,' hr',imin,' min'
	endif
	return
	
 999    write(llog,*) 'Data file not found!'
        stop
        end
c==============================================================================
	subroutine dumpdat(outfile)
c==============================================================================
	implicit none
	character*80 outfile
	include 'seek.inc'
	integer lun,i
	write(*,'(a)') ' Dumping '//outfile(1:index(outfile,' ')-1)
	call glun(lun)

	open(unit=lun,file=outfile,status='unknown',
     &          form='unformatted')
	write(lun) ntim,real(tsamp),refdm,refac
        write(lun) (series(i),i=1,ntim)
	close(unit=lun)
	end
c==============================================================================
