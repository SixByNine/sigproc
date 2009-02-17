c program to determine the statistics of 1 bit filter data
c orginally SJ - giant pulse searching code
c adapted SJ - Feb 2006 for MMB survey
c
c
	program chaninfo

	implicit none

c PARAMETER DEFINITIONS
c nbuf=buffer size in bytes
c chmax=max channels
c nblkmax=max number of blocks allowed
	integer chmax,nbuf,nblkmax,hdrsize
	parameter(nbuf=49152,hdrsize=640,chmax=513,nblkmax=10000)

c VARIABLE DECLARATIONS
	byte bbuf(nbuf)
	integer b,i,j,k,bval,samp_blk,lf
	integer nchans,nchans8,nch,ngulp,ival,ibit,nsamp
	integer bit_offset(chmax),chsum(chmax)
	integer itab(-1:7,-128:127),nbyte,nargs
	real chave(chmax),threshold,chbw,frch1,freq
	real*8 rms,chrms(chmax)
	character fname*70,hdrfname*70,ans*1,string*70,hdr*(640)

	integer blockgood(nblkmax),ngood,nbad,lastngood,niter
	real blockmean(nblkmax),globmean,globrms
	real lastglobmean,lastglobrms
	real bbm,sample_amp,sigma,amp_expected,rms_expected
c	integer band_amp(chmax)

c read in the arguments
	nargs = IARGC()
        if(nargs.lt.1)then
	  write(*,*)
          write(*,'(''chaninfo - time series stats of 1-bit data'')')
	  write(*,*)
          write(*,'(''usage: chaninfo file (nblk) (cthr) (bthr)'')')
	  write(*,*)
          write(*,'(''file - name of raw .dat file '')')
          write(*,'(''nblk - (optional) blocks to read (def=all)'')')
          write(*,'(''cthr - (optional) threshold for bad channels'',
     +              '' (def=+/-0.2)'')')
          write(*,'(''bthr - (optional) threshold for bad blocks'',
     +              '' (def=3 sigma)'')')
	  write(*,*)
          stop
        endif
	call getarg(1,fname)
	open(unit=10,file=fname,status='old',form='unformatted',err=990)
        ngulp=0
	if(nargs.gt.1)then
	  call getarg(2,string)
	  read(string,*)ngulp
        endif
	if(ngulp.eq.0)ngulp=nblkmax
        threshold=0.2
	sigma=3.0
	if(nargs.gt.2)then
          call getarg(3,string)
          read(string,*)threshold
        endif
	if(nargs.eq.4)then
          call getarg(4,string)
          read(string,*)sigma
        endif

c read the header for useful info
	lf=index(fname,' ')-1
	hdrfname=fname(1:lf-4)//'.hdr'
	open(unit=11,file=hdrfname,status='old',form='unformatted',err=991)
	read(11,err=992)hdr
	close(11)
	read(hdr(222:224),'(i3)',err=992)nchans
	if(nchans.gt.chmax)then
	  write(*,'('' ERROR - nchans is greater than chmax !! '')')
	  stop
	endif
	read(hdr(205:212),'(f8.0)',err=992)chbw
	if(chbw.gt.0.)then
	  write(*,'('' WARNING - chbw is positive - flipping channels'')')
	endif
	read(hdr(242:253),'(f12.0)')frch1

	nchans8=nchans/8
c expected mean and rms of a sample
	amp_expected=nchans/2.0
	rms_expected=0.5*sqrt(nchans/2.0)
c samples per block
	samp_blk=nbuf*8/nchans

c setup the table
c For byte value = ival, if bit ibit of byte is set, then tab(ibit,ival)=1
        do ival=-128,127
          do ibit=0,7
            if(btest(ival,ibit))then
              itab(ibit,ival)=1
            else
              itab(ibit,ival)=0
            endif
          enddo
        enddo
	do i=1,nchans
	  bit_offset(i)=mod(i-1,8)
	  chsum(i)=0
	  chrms(i)=0.0
	enddo

c read in the data and compute the mean and rms of each channel
c compute the mean of each block of data
	nsamp=0
 15 	do b=1,ngulp
 	  read(10,end=20)bbuf
	  nbyte=0
	  bbm=0.0
c	  do j=1,nchans8
c	    band_amp(j)=0
c	  enddo
	  do i=1,samp_blk
	    nsamp=nsamp+1
	    sample_amp=0.0
	    nch=0
	    do j=1,nchans8
	      nbyte=nbyte+1
	      do k=1,8
	        nch=nch+1
	        bval=itab(bit_offset(nch),bbuf(nbyte))
	        chsum(nch)=chsum(nch)+bval
	        chrms(nch)=chrms(nch)+bval*bval
	        sample_amp=sample_amp+bval
c	        band_amp(j)=band_amp(j)+bval
	      enddo
	    enddo
	    bbm=bbm+sample_amp
	  enddo
c	  write(77,'(24i8)')(band_amp(j),j=1,nchans8)
	  blockmean(b)=bbm/samp_blk
	enddo

 20	b=b-1
	close(10)
	open(unit=33,file='chan.info',form='formatted',status='unknown')
	write(33,'('' Finished reading '',i4,'' blocks '')')b
	write(33,'('' Number of samples per channel = '',i8)')nsamp
	rms=0.5/sqrt(real(nsamp))
	write(33,'('' Expected rms per channel = '',f12.8)')rms
	write(33,*)

	write(33,'('' Ch       Freq          Sum      Ave        RMS'')')
	do i=1,nchans
	  freq=frch1+(i-1)*chbw
	  chave(i)=real(chsum(i))/nsamp
          chrms(i)=sqrt(chrms(i)/nsamp - chave(i)*chave(i))
	  write(33,'(i4,f12.4,i12,f10.6,f12.6)')i,freq,
     +              chsum(i),chave(i),chrms(i)
	enddo
	close(33)

c determine the bad channels and write out to bad.chans
c good channels to good.chan
c NOTE - these are inverted if chbw>0 for dice routine
	open(unit=11,file='bad.chans',form='formatted',status='unknown')
	open(unit=12,file='good.chans',form='formatted',status='unknown')
	if(chbw.lt.0.)then
          do i=1,nchans
            if(chave(i).gt.0.5+threshold .or.
     +         chave(i).lt.0.5-threshold)then
              write(*,'('' Bad channel '',i4,'' with average '',
     +          f10.6)')i,chave(i)
              write(11,*)i,chave(i)
            else
              write(12,*)i
	    endif
	  enddo
	else
          do i=nchans,1,-1
            if(chave(i).gt.0.5+threshold .or.
     +         chave(i).lt.0.5-threshold)then
              write(*,'('' Bad channel '',i4,'' with average '',
     +          f10.6)')nchans-i+1,chave(i)
              write(11,*)nchans-i+1,chave(i)
            else
              write(12,*)nchans-i+1
	    endif
	  enddo
	endif
        close(11)
        close(12)

c do the block statistics
	do i=1,nblkmax
	  blockgood(i)=1
	enddo
	lastglobmean=0.0
	lastglobrms=1.e6
	niter=0

c this is an iterative loop which marks as bad those blocks
c more than "sigma" from the mean. The loop continues until
c either the number of bad blocks is constant or for 10 trials
 100	globmean=0.0
	globrms=0.0
	lastngood=ngood
	ngood=0
	nbad=0
	do i=1,b
	  if(blockmean(i).gt.lastglobmean+sigma*lastglobrms .or. 
     +	     blockmean(i).lt.lastglobmean-sigma*lastglobrms) then
	     blockgood(i)=0
             nbad=nbad+1
	  else
	    ngood=ngood+1
	    globmean=globmean+blockmean(i)
	    globrms=globrms+blockmean(i)*blockmean(i)
	  endif
	enddo
	globmean=globmean/ngood
	globrms=sqrt(globrms/ngood - globmean*globmean)
	lastglobmean=globmean
	lastglobrms=globrms
	niter=niter+1
	if(niter.eq.10)goto 120
	if(lastngood.ne.ngood)goto 100

c write out the information
 120	continue
	write(*,*)
	write(*,'('' There are '',i8,'' good blocks'')')ngood
	write(*,'('' There are '',i8,'' bad blocks'')')nbad
	write(*,'('' Iterations : '',i4)')niter
	open(unit=44,file='block.info',form='formatted',status='unknown')
	write(44,'('' Block   Mean   Status '')')
	open(unit=45,file='bad.blocks',form='formatted',status='unknown')
	nsamp=1
	do i=1,b
	  write(44,*)i,blockmean(i),blockgood(i)
	  if(blockgood(i).eq.0)write(45,*)nsamp,nsamp+samp_blk-1
	  nsamp=nsamp+samp_blk
	enddo
	close(44)
	close(45)
	goto 999

 990    write(*,'('' ERROR - could not open data file !!! '')')
	stop
 991    write(*,'('' ERROR - could not open header file !!! '')')
	stop
 992    write(*,'('' ERROR - problem reading header file !!! '')')
	stop

 999	continue
	end
