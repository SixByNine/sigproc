c==============================================================================
      program gray
c==============================================================================
c
c  This program plots EPN records to the standard output in a pseudo-grey
c  scale style suitable for viewing on a simple terminal with no graphic
c  capabilities.
c
c  Last modified 98/04/28 - writes out ASCII array of the grey scale for plot
c
c  dunc@mpifr-bonn.mpg.de
c
c==============================================================================
      
      implicit none

      include 'epnhdr.inc'
      real tmp(maxbin)
      integer nepnrec,nrec,nscr,nplt,step,next,nbins,nch,ichan,lun,istat
      integer narg, lfil,i,j,nskip, nread,iargc,ln,nc,hh,mm,nborg,binmax
      double precision date,ss
      logical chflag(maxblk),first,centre,help
      data first/.true./
      real*8 t1,tsub
      real dmin,dmax,rmsfac,grey(80),prmax
      parameter (nscr=60,nc=10)
      character*1 gry(nc),cprof(maxblk)*80,option*80
c
c     Here are the characters in increasing order of intensity...
c      
       gry(1)=' '
       gry(2)=' '
       gry(3)=' '
       gry(4)=','
       gry(5)=':'
       gry(6)='o'
       gry(7)='*'
       gry(8)='@'
       gry(9)='$'
      gry(10)='#'
      rmsfac=0.0
c
c     Prompt user if silly/no options entered...
c      
      narg=iargc()
      help=.false.
      if (narg.lt.1) then
         call glun(lun)
         open(unit=lun,file='file',status='old',iostat=istat)
         if (istat.eq.0) then
            read(lun,'(a)') filename
            close(lun)
         else
            help=.true.
         endif
      endif
         
      if (help) then
      write(*,'('' GREY : A program to greyplot EPN data.'')')
         write(*,'('' No options were specified!'')')
         write(*,'('' usage: plotg [filename] <options>'')')
         write(*,'('' N.B. Input file must be in EPN format!!'')')
	 write(*,'('' -r: set start record number (def=1)'')')
	 write(*,'('' -s: set min=0 & max=n*rms (def=autoscale)'')')
	 write(*,'('' -c: centres profile to first record (optional)'')')
         write(*,'('' Comments/Bugs etc. -> dunc@mpifr-bonn.mpg.de'')')
         stop 
      else
         if (narg.ge.1) call getarg(1,filename)
         lfil=index(filename,'.ser')-1
         if (lfil.eq.0) lfil=index(filename,'.dis')-1
         if (lfil.gt.0) filename=filename(1:lfil)
         lfil=index(filename,' ')-1
	 if (index(filename,'.epn').eq.0) then
	   filename=filename(1:lfil)//'.epn'
	   lfil=lfil+4
         endif
         nrec=nepnrec(filename)
         if (nrec.eq.0) then 
           write(*,'('' EPN file: '',a,'' not found!'')')
     &     filename(1:lfil)
           stop
         endif
      endif
c
c     Read the first record to find out the number of channels...
c      
      padout=.false.
      readwri=-1
      recno=1
      call rwepn(filename, readwri, recno, padout)
      t1=tstart(1)
      recno=2
      call rwepn(filename, readwri, recno, padout)
      tsub=tstart(1)-t1
c
c     Set up default values...
c      
      nskip=0
      nread=0
      recno=0
      do i=1,maxblk
         chflag(i)=.false.
      enddo
      chflag(1)=.true.
      centre=.false.
      binmax=0
c
c     Read optional extras from the standard input if necessary...
c      
      if (narg.ge.2) then
         do i=2,narg
  	   call getarg(i,option)
           if (index(option,'-r').gt.0) then
              read(option(3:),'(i5)') recno
              recno=recno-1
              if (recno.lt.0) recno=0
           else if (index(option,'-s').gt.0) then
              read(option(3:),'(f8.2)')rmsfac
           else if (index(option,'-c').gt.0) then
              centre=.true.
           else if (index(option,'-S').gt.0) then
              read(option(3:),*) binmax
           endif
         enddo
      endif
c
c     Main loop
c
      call glun(lun)
      open(lun,file='grey',status='unknown')
      do while(recno.ne.-999)

        recno=recno+1
        call rwepn(filename, readwri, recno, padout)
        date=epoch+tstart(1)/8.64e10
        call dattim(date,hh,mm,ss)
        ln=index(cname,' ')-1
        nch=0

        if (first.and.centre) then
           binmax=0
           prmax=-1.0e32
           do i=1,nbin
              if (rawdata(1,i).gt.prmax) then
                 prmax=rawdata(1,i)
                 binmax=i
              endif
           enddo
           binmax=nbin/2-binmax
        endif
        
        do ichan=1,npol
        if (chflag(ichan)) then

        nch=nch+1
        do i=1,nbin
           tmp(i)=0.0
        enddo
        nplt=nscr
        step=nbin/nscr
        if (step.lt.1) step=1
        next=step
        j=1
        do i=1,nbin
           tmp(j)=tmp(j)+rawdata(ichan,i)
           if (i.eq.next) then
              j=j+1
              next=next+step
           endif
        enddo
        nbins=j-1

        if (nbins.gt.nscr) nbins=nscr

c        if (centre) call sprof(tmp,nbins,binmax)
        if (binmax.ne.0) call sprof(tmp,nbins,binmax)

        dmin=+1.0e32
        dmax=-1.0e32

        do i=1,nbins
           dmax=max(tmp(i),dmax)
           dmin=min(tmp(i),dmin)
        enddo

        if (rmsfac.gt.0.) then
          dmin=0.0
          dmax=rmsfac*rms(ichan)
        endif
        cprof(nch)=' '

        do i=1,nbins
          j=(tmp(i)-dmin)/(dmax-dmin)*nc+1
          if (j.gt.nc) j=nc
	  if (j.lt.1) j=1
          cprof(nch)(i:i)=gry(j)
          grey(i)=tmp(i)
        enddo

        nborg=nbins
        do while(nbins+nborg.lt.nscr)
           do i=1,nborg
              cprof(nch)(nbins+i:nbins+i)=cprof(nch)(i:i)
              grey(nbins+i)=grey(i)
           enddo
           nbins=nbins+nborg
        enddo

        endif
        enddo
        if (recno.ne.-999) then
        write(*,'('' |'',i4.4,''|'',i2.2,'':'',i2.2,'':'',i2.2,
     &  ''|'',a,''|'')')
     &   recno,hh,mm,int(ss),cprof(1)(1:nbins)
        if (first) then
C           write(lun,'(x,a72)') history
           write(lun,*)
           write(lun,*) nbins,1.0/real(nborg),real(tsub)/1.0e6/60.0
           first=.false.
        endif
        write(lun,*) (grey(i),i=1,nbins)
        endif
      enddo
      close(lun)
c      write(*,*)1.0e6*real(papp(1)**2.0)/real(date*86400.0)/real(nbins),
c     &           ' us error in period would cause a drift of one bin...'
      end
