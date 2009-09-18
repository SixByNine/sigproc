c==============================================================================
      program quickplot
c==============================================================================
c
c     Create postscript summary plot of output from the "quicklook" script 
c
      implicit none

      include 'epnhdr.f'
      integer nepnrec,i,j,nchan,nshift,nsubints,slen,nc,lun,istat
      real x(maxbin),y(maxbin),xmin,xmax,ymin,ymax,dat(maxbin,1024),mean
      real psnr,ssnr,smmax,fc,bw,f1,cb,pdm,ts,to,window,phasestart,trms
      integer kwmax,greyscale,indx(maxbin),npts,n90,ndms,nbits
      logical filexists,susfile,cent,summary
      real*8 pfold,pfind,susdm,sussn
      real z(maxbin),dmmodelrms,bestrms,pms,dc,testrms,bestdc
      character*80 telescope,machine,source,comline,date,mjd,fstring
      character*80 period,tstring,junk,dmc
c
c     Initialise some variables...
c
      padout=.false.
      readwri=-1
      recno=1
      inquire(file='subbands.epn',exist=filexists)
      if (.not.filexists) stop 'subbands.epn not found!'
      inquire(file='subints.epn',exist=filexists)
      if (.not.filexists) stop 'subints.epn not found!'
      inquire(file='asciiheader',exist=filexists)
      if (.not.filexists) stop 'asciiheader file not found!'
      inquire(file='timeseries',exist=filexists)
      if (.not.filexists) stop 'timeseries file not found!'
      inquire(file='phasestart',exist=filexists)
      if (filexists) then
         call glun(lun)
         open(lun,file='phasestart',status='old')
         read(lun,*) phasestart
         close(lun)
      else
         phasestart=0.0
      endif

      nchan=nepnrec('subbands.epn')
      if (nchan.lt.1) stop 'no profiles found in subbands.epn'
      call rwepn('subbands.epn', readwri, nchan+1, padout)

      nsubints=nepnrec('subints.epn')-1
      if (nsubints.lt.1) stop 'no profiles found in subints.epn'
      call rwepn('subbands.epn', readwri, 1, padout)
      pfold=papp(1)*1000.0

      call glun(lun)
      open(lun,file='asciiheader',status='old',err=1)
      read(lun,'(a)',err=2) telescope
      read(lun,'(a)',err=2) machine
      read(lun,'(a)',err=2) filename
      read(lun,'(a)',err=2) source
      read(lun,'(a)',err=2) mjd
      read(lun,'(a)',err=2) date
      read(lun,*,err=2) fc
      read(lun,*,err=2) bw
      read(lun,*,err=2) f1
      read(lun,*,err=2) cb
      read(lun,*,err=2) nc
      read(lun,*,err=2) nbits ! new
      read(lun,*,err=2) ts
      read(lun,*,err=2) to
      read(lun,'(a)',err=2) period
      read(lun,'(a)',err=2) dmc
      read(dmc,*) pdm
    
c      write(fstring,'(''F\\dOBS\\u: '',f6.3,'' GHz '',
c     &      ''Bandwidth:'',f5.1,'' MHz Nchans: '',i4,
c     &      ''  t\\dsamp\\u: '',f6.0,''\\gms'')') 
c     & fc/1000.0,bw,nc,ts
      close(lun)
      
      greyscale=0
      call getarg(1,comline)
      call getarg(2,junk)
	
      if (junk.eq.'greyscale'.or.junk.eq.'grayscale') greyscale=1
      if (junk.eq.'centre'.or.junk.eq.'center') then
         cent=.true.
      else
         cent=.false.
      endif
	cent=.true.
      call getarg(3,junk)
      if (junk.eq.'summary') then
	summary=.true.
      else
	summary=.false.
      endif
      if (comline.eq.'?') then
         call pgbegin(0,'?',1,1)
      else 
         if (summary) then
           call pgbegin(0,filename(1:slen(filename))//'.ps/ps',1,1)
	   write(*,'(a)') filename(1:slen(filename))//'.png'
         else
           call pgbegin(0,filename(1:slen(filename))//'.ps/vps',1,1)
	   write(*,'(a)') filename(1:slen(filename))//'.ps'
	 endif
      endif
      if (comline.eq.'greyscale'.or.comline.eq.'grayscale') greyscale=1
      inquire(file=comline,exist=filexists)
      if (filexists) then
         susfile=.true.
      else
         susfile=.false.
      endif
      call pgscf(2)
      call pgsch(0.9)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (summary) then
	   call pgvport(0.08,0.98,0.6,0.7)
	   call pgwindow(0.0,1.0,0.0,1.0)
	   call pgtext(0.0,1.1,filename(1:slen(filename))//
     &               '    P = '//period(1:slen(period))//' ms'//
     &               '   DM = '//dmc(1:slen(dmc))//' cm-3 pc'
     &                 )
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmin=+1.0e32
      xmax=-1.0e32
      ymin=+1.0e32
      ymax=-1.0e32
      call glun(lun)
      open(lun,file='timeseries',status='old')
      i=0
      istat=0
      do while(istat.eq.0) 
         i=i+1
         read(lun,*,iostat=istat) x(i),y(i)
         if (istat.eq.0) then
            xmin=min(xmin,x(i))
            xmax=max(xmax,x(i))
         endif
      enddo
      close(lun)
      npts=i-1
      n90=npts-npts/10
      call indexx(npts,y,indx)
      mean=0.0
      do i=1,n90
         mean=mean+y(indx(i))
      enddo
      mean=mean/real(n90)
      do i=1,npts
         y(i)=y(i)-mean
      enddo
      trms=0.0
      do i=1,n90
         trms=trms+(y(indx(i))*y(indx(i)))
      enddo
      trms=sqrt(trms/real(n90))
      do i=1,npts
         y(i)=y(i)/trms
         ymin=min(ymin,y(i))
         ymax=max(ymax,y(i))
      enddo
      if (summary) goto 111
      if (susfile) then
         call pgvport(0.125,0.48,0.67,0.79)
      else 
         call pgvport(0.125,0.925,0.67,0.79)
      endif
      call pgwindow(xmin,xmax,ymin,10.0)
      
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
      call pgline(npts,x,y)
      call pglabel('Telescope Time (s)','Flux (rms~1)',' ')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(fstring,'(''Folded '',i5,'' s of'',
     &      i4,'' x'',f7.1,''-kHz '',i2,
     &      ''-bit filterbank data t\\dsmp\\u:'',
     &      i5,''\\gms'')')
     &      int(xmax),nc,1000*abs(bw)/real(nc),nbits,int(ts)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 111  if (susfile) then
         xmin=+1.0e32
         xmax=-1.0e32
         ymin=+1.0e32
         ymax=-1.0e32
         call glun(lun)
         open(lun,file=comline,status='old',err=1)
         read(lun,'(a)',iostat=istat) junk
         if (junk(1:1).eq."#") then
	   read(junk(2:),*) pfind,sussn,susdm
	 else
	   read(junk,*) pfind,sussn,susdm
	 endif
         i=0
         do while (istat.eq.0) 
            i=i+1
            read(lun,*,iostat=istat) ndms,x(i),y(i)
         enddo
         close(lun)
         ndms=i-1
         bestrms=1.0e32
	 bestdc=0.01
         pms = real(pfold)
         do i=1,999
           dc  = real(i)/1000.
	   testrms=dmmodelrms(ndms,x,y,z,xmin,xmax,pdm,bw,fc,pms,dc)
           if (testrms.lt.bestrms) then
		bestdc=dc
		bestrms=testrms
	   endif
         enddo 
         testrms=dmmodelrms(ndms,x,y,z,xmin,xmax,pdm,bw,fc,pms,bestdc) 
         if (summary) then
           call pgvport(0.3,0.5,0.5,0.7)
	 else
           call pgvport(0.57,0.925,0.67,0.79)
	 endif
         ymin=y(1)
	 ymax=y(1)
	 do i=1,ndms
	   ymin=min(y(i),ymin)
	   ymax=max(y(i),ymax)
         enddo
         call pgwindow(xmin,xmax,ymin,ymax*1.1)
         call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
c         call pgsls(4)
c         call pgline(ndms,x,y)
c         call pgslw(4)
         call pgsls(1)
         call pgline(ndms,x,z)
         call pgslw(1)
         call pgpoint(ndms,x,y,17)
         call pglabel('Trial DM (cm\\u-3\\d pc)','S/N',' ')
	 call pgwindow(0.0,1.0,0.0,1.0)
         ssnr=ymax
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ymax=-1.0e32
      call rwepn('subints.epn', readwri, nsubints, padout)
      do i=1,nbin
         if (rawdata(1,i).gt.ymax) then
            ymax=rawdata(1,i)
            j=i
         endif
      enddo
      nshift=nbin/2-j
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmin=+1.0e32
      xmax=-1.0e32
      ymin=+1.0e32
      ymax=-1.0e32
      do i=1,nsubints
         call rwepn('subints.epn', readwri, i, padout)

         do j=1,nbin
            x(j)=real(j)*real(tbin)*1.0e-3
            xmin=min(xmin,x(j))
            xmax=max(xmax,x(j))
            y(j)=rawdata(1,j)
            ymin=min(ymin,y(j))
            ymax=max(ymax,y(j))
         enddo
         if (cent) call sprof(y,nbin,nshift)
         do j=1,nbin
            dat(j,i)=y(j)
         enddo
      enddo
      
      xmin=0.0
      xmax=1.0
      ymin=0.0
      ymax=real(nsubints)
      if (summary) then
        call pgvport(0.56,0.76,0.5,0.7)
      else
        call pgvport(0.125,0.925,0.28,0.43)
      endif
      if (greyscale.eq.1) then
         call grayscale(dat,maxbin,1024,nbin,nsubints,1) 
      else 
         call quickgrey(dat,maxbin,1024,nbin,nsubints,1) 
      endif
      window=real(nbin)*real(tbin)*1.0e-3/pfold
c      call pgwindow(phasestart,phasestart+window,real(nsubints+1),0.0)
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
      call pglabel(' ','Subintegration',' ')
      if (summary) 
     & call pglabel('Pulse phase (turns)','Subintegration',' ')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,nchan
         call rwepn('subbands.epn', readwri, i, padout)
         do j=1,nbin
            y(j)=rawdata(1,j)
         enddo
         if (cent) call sprof(y,nbin,nshift)
         do j=1,nbin
            dat(j,i)=y(j)
            ymin=min(ymin,y(j))
            ymax=max(ymax,y(j))
         enddo
      enddo
      xmin=0.0
      xmax=1.0
      ymin=0.0
      ymax=real(nchan)
      if (summary) then
        call pgvport(0.82,0.98,0.5,0.7)
      else
        call pgvport(0.125,0.925,0.45,0.6)
      endif
      call pgwindow(1.0,real(nbin),0.0,real(nchan+1))
      if (greyscale.eq.1) then
         call grayscale(dat,maxbin,1024,nbin,nchan,1)
      else
         call quickgrey(dat,maxbin,1024,nbin,nchan,1)
      endif
c      call pgwindow(0.0,window,0.0,real(nchan+1))
      if (summary) then
	call pglabel('Pulse phase (turns)','Frequency band',' ')
        call pgbox('bcnst',0.0,0,'bcnt',0.0,0)
      else
        call pglabel('','Frequency band',' ')
        call pgbox('bcst',0.0,0,'bcnt',0.0,0)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 222  call rwepn('subints.epn', readwri, nsubints+1, padout)
      ymin=+1.0e32
      ymax=-1.0e32
      do j=1,nbin
         x(j)=real(j-1)*real(tbin)*1.0e-3
         x(j)=j
         y(j)=rawdata(1,j)
         ymin=min(ymin,y(j))
         ymax=max(ymax,y(j))
      enddo
      xmin=0.0
      xmax=pms
      xmin=0.5
      xmax=real(nbin)+0.5
      if (summary) then
        call pgvport(0.08,0.23,0.5,0.7)
      else
        call pgvport(0.125,0.925,0.1,0.25)
      endif
      ymin=ymin-(ymax-ymin)*0.05
      ymax=ymax+(ymax-ymin)*0.05
      call pgwindow(xmin,xmax,ymin,ymax)
      call pgbox('bcnst',0.0,0,'bc',0.0,0)
      call pglabel('Pulse Phase (bins)','Flux density',' ')
      if (cent) call sprof(y,nbin,nshift)
      call pgmove(xmin,y(1))
      do i=1,nbin-1
         call pgdraw(x(i)+0.5,y(i))
         call pgdraw(x(i)+0.5,y(i+1))
      enddo
      call pgdraw(xmax,y(nbin))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.cent) call sprof(y,nbin,nshift)
      call smooth(y,nbin,kwmax,psnr,smmax)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (summary) then
        call pgwindow(0.0,1.0,0.0,1.0)
	write(junk,'(f6.1)') psnr
	call pgtext(0.7,0.8,junk)
	goto 333
      endif
      call pgvport(0.125,0.925,0.1,1.0)
      call pgvport(0.06,0.925,0.1,1.0)
      call pgwindow(0.0,1.0,0.0,1.0)
      call pgtext(0.0,0.9,telescope(1:slen(telescope))//
     &     ' '//machine(1:slen(machine))//': '//
     &     filename(1:slen(filename))//'.ps')
      call pgtext(0.0,0.875,'J2000 coords: '//source(1:slen(source))//
     &     ' MJD: '//mjd(1:slen(mjd))//' Date: '//date(1:len(date)))
      call pgtext(0.0,0.85,fstring)
      write(tstring,'(''P\\dfromSEEK\\u: '',f16.10,
     & '' ms  DM: '',f7.2,'' cm\\u-3\\d pc'',''  SEEK S/N: '',f6.1)') 
     &   pfind,pdm,ssnr
      call pgtext(0.0,0.825,tstring)
      write(tstring,'(''P\\doptimized\\u: '',f16.10,
     & '' ms N\\dbins\\u: '',i4,
     & '' DC: '',f4.1,''%   PROF S/N: '',f6.1)') 
     & pfold,nbin,bestdc*100,psnr
      call pgtext(0.0,0.8,tstring)
 333  call pgend
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      stop
 1    stop 'ERROR: opening asciiheader'
 2    stop 'ERROR: reading asciiheader'
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function dmmodelrms(ndms,x,y,z,xmin,xmax,pdm,bw,fc,pfold,dc)
      integer ndms
      real x(*),y(*),z(*),xmin,xmax,pdm,bw,fc,pfold,dc,weff
      real ssq,zmax,wint,ymax
      integer i
      zmax=0.0
      wint=dc*pfold
      xmin=0.0
      xmax=0.0	
      ymax=0.0
      do i=1,ndms
        weff=sqrt(wint**2.0+(8.3e6*abs(pdm-x(i))*bw/fc/fc/fc)**2.0)
        if(weff.gt.pfold) then
  	  z(i)=0.0
	  if (xmin.ne.0.0.and.xmax.eq.0.0) xmax=x(i)*1.1
	else
	  if (xmin.eq.0.0.and.xmax.eq.0.0) xmin=x(i)*0.9
	  z(i)=sqrt(pfold-weff)/sqrt(weff)
	endif
        ymax=max(ymax,y(i))
	zmax=max(zmax,z(i))
      enddo
      ssq=0.0
      do i=1,ndms
        z(i)=z(i)*ymax/zmax
        ssq=ssq+(y(i)-z(i))**2.0
      enddo
      dmmodelrms=sqrt(ssq/real(ndms))
      if (xmax.eq.xmin) then
	xmin=x(1)
	xmax=x(ndms)
      endif
      if (xmax.eq.0.0) xmax=x(ndms)*1.1
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function slen(string)
      character string*(*)
c     
c     OBTAIN THE LOCATION OF THE LAST NON-SPACE CHARACTER.
c     
      integer ilen
c     search for the first null
      ilen = len(string)
c     use the position of the first null
      do 1 i = ilen, 1, -1
c     
c     LENGTH FOUND.
c     
         if (string(i:i) .ne. char(32) .and.
     &       string(i:i).ne.char(0)) then
            slen = i
            return 
         end if
c     
c     STRING IS ALL SPACES OR ZERO LENGTH.
c     
    1 continue
      slen = 0
c     
c     END OF INTEGER FUNCTION LENGTH.
c     
      return 
      end
      subroutine quickgrey(dat,nxd,nyd,nx,ny,flip)

c  Coarse grey-scale plot using PG plot routines
c  Assumed that viewport and window defined outside routine
        integer nxd,nyd
        integer nsym,flip
        parameter(nsym=10)
	integer*4 ksym(nsym)
	real*4 dat(nxd,nyd)
	data ksym/1,20,21,2,3,14,15,17,16,18/

	s=0.
	ss=0.
	smin=1.e30
	smax=-smin
	do j=1,ny
	  do i=1,nx
	    aa=dat(i,j)
	    s=s+aa
	  enddo
	enddo
	s=s/(nx*ny) 
	do j=1,ny
	  do i=1,nx
            aa=dat(i,j)-s
	    ss=ss+aa**2
	    if(aa.gt.smax)smax=aa
	    if(aa.lt.smin)smin=aa
	  enddo
	enddo
	rms=sqrt(ss/(nx*ny))

        if (flip.eq.1) then
           call pgwindow(1.0,real(nx),real(ny+1),0.0)
        else 
           call pgwindow(1.0,real(nx),0.0,real(ny+1))
        endif
	do j=1,ny
	  do i=1,nx
	    k=min(int((dat(i,j)-s)/rms),nsym)
	    if(k.gt.0)then
	      x=i
	      y=j
	      call pgpoint(1,x,y,ksym(k))
	    endif
	  enddo
	enddo

	return
	end



c==============================================================================
C nicked from pgplot
C*GRGLUN -- get a Fortran logical unit number (Sun/Convex-UNIX)
C+
      SUBROUTINE GLUN(LUN)
      INTEGER LUN
C
C Get an unused Fortran logical unit number. 
C Returns a Logical Unit Number that is not currently opened.
C After GRGLUN is called, the unit should be opened to reserve
C the unit number for future calls.  Once a unit is closed, it
C becomes free and another call to GRGLUN could return the same
C number.  Also, GRGLUN will not return a number in the range 1-9
C as older software will often use these units without warning.
C
C Arguments:
C  LUN    : receives the logical unit number, or -1 on error.
C--
C 12-Feb-1989 [AFT/TJP].
C DRL adapted to subroutine GLUN for use with stand-alone software
C 16-Jul-1993 @ JB
c==============================================================================
      INTEGER I
      LOGICAL QOPEN
C---
      DO 10 I=99,10,-1
          INQUIRE (UNIT=I,  OPENED=QOPEN)
          IF (.NOT.QOPEN) THEN
              LUN = I
              RETURN
          END IF
   10 CONTINUE
C none left
      STOP 'RAN OUT OF LUNs!' 
      END
c==============================================================================
c
c     File: "rwepn.f"
c     
c     This file contains the generalised EPN format reading/writing routines
c     It has been tested on an HP system with f77 -c rwepn.f to produce the 
c     object file. It requires the accompanying include file "epnhdr.f"
c
c==============================================================================
       subroutine rwepn(filename2,readwri,recno,padout)
c==============================================================================
c
c     Subroutine to read or write data in EPN format. All parameters are
c     defined in the include file "epnhdr.f" which should lie in the
c     same directory as this file. Passed down variables to this routine
c     are:
c
c     filename - character*80  - name of this EPN file
c     readwri  - integer*4     - switch: -1 => read, +1 => write
c     recno    - integer*4     - record number to read
c
c     Usage:
c
c     To incorporate this routine into your pulsar data reduction
c     package, include the file "epnhdr.f" in your source code
c     and set up the appropriate variables before.
c
c     Features:
c
c     When writing is requested, the routine will keep the file
c     open for successive writes by default. To override this, pass
c     down recno as -1, the presently opened file will be closed
c     and a new one opened.
c
c     If an error occurs in reading, recno is returned as -999
c
c     By the very nature of the format, record lengths must be fixed.
c     So if you want to have different lengthed data streams in the
c     same archive you need to switch "padout" on. As its name suggests,
c     this pads out the binned data to its maximum length so that
c     the records are fixed. A draw back is of course that it uses
c     up more space. Therefore if it is known that all the records
c     are of the same lenght anyway (eg for a stream of single pulses)
c     then it is recommended to switch padout off.
c
c==============================================================================
c
c     Created by Dunc Lorimer (dunc@mpifr-bonn.mpg.de) November 1996....
c
c     12/11/1996   Initial version
c     15/11/1996   Implemented common blocks
c     21/05/1997   Changed to HEX storage
c     01/07/1997   Added tres and fluxflag to header
c
c==============================================================================
c
      implicit none
      character*(*) filename2
      include 'epnhdr.f' 
c
c     Data buffers...
c
      integer*4 intdata(maxbin) 
      character*4 chrdata(maxbin)
      real*4 tmp(maxbin)
c
c     The maximum length of an EPN file is the 6 line header, plus
c     maxblk sets of maxbin bin data streams. Each stream of maxbin
c     bins occupies (maxbin/20) sets of 80 character lines plus a 2
c     line header. The number of charcters per record is:
c
c     (6+maxblk*(maxbin/20+2))*80 
c
c     One can tune this depending on the size of the data-sets...
c     the default is maxblk=8, maxbin=4160 -> maxrec=66080
c     suitable smaller setting for e.g. a PC would be...
c                    maxblk=4, maxbin=1040 -> maxrec=17760.
c     N.B. Both maxbin and maxblk are specified in "epnhdr.f"
c
      integer maxrec
      parameter(maxrec=(6+maxblk*(maxbin/20+2))*80)
      character line(maxrec/80)*80, recrd*(maxrec)
c
c==============================================================================
c
c     Local variables
c
      integer lun, lfil, nlin, next, i, j, k, l, iras, ipol,
     &        wrecno, recln, maxint, lct, length
      real rras, pvno, vno
      character*1 dsign, lin1*12, line40*40
      logical reading,writing,filex,opened,first
      save
      data opened/.false./, first/.true./
      line40='----------------------------------------'
c
c==============================================================================
c
c     Sort out whether reading or writing
c
      reading=(readwri.eq.-1)
      writing=(.not.reading)
      if (readwri.ne.-1.and.readwri.ne.1) then
         write(*,'(''Unknown read/write option passed to rwepn!!'')')
         write(*,'(''Your "readwri" integer reads: '',i4)')readwri
         stop
      endif
c
c     If writing required and recno=-1, then close a previously
c     opened file (if one is open) and proceed
c
      if (writing.and.recno.eq.-1.and.opened) then
         close(unit=lun)
         opened=.false.
      endif
c
c     Get a free logical unit number
c
      if (.not.opened) call glun(lun)
c
c     Check whether the file exists or not....
c     
c     For reading: the routine will stop if not
c     For writing: the routine will append to the previous file
c
      lfil=length(filename2)
      filename=filename2(1:lfil)
      inquire(file=filename(1:lfil),exist=filex)
c
c     Open the input file
c
      if (reading) then
         if (filex) then
c
c          Read in the very first line, and work out the
c          required record length...(number of 80 chr lines)
c
           open(unit=lun,file=filename(1:lfil),status='old',
     &          access='direct',form='unformatted',recl=12)
           read(lun,rec=1) lin1
           read(lin1,'(8x,i4)') counter
           close(unit=lun)
           recln=counter*80
c
c          Now open as a direct access unformatted file
c
           open(unit=lun,file=filename(1:lfil),status='old',
     &          access='direct',form='unformatted',recl=recln)
           if (recno.lt.1) recno=1
           read(lun,rec=recno,err=999) recrd(1:recln)
           close(unit=lun)
         else
           stop 'EPN file does not exist!'
         endif
      endif
      if (writing) then
c
c       Establish out how many 'lines' will make up this record...
c
        if (padout) then
           recln=maxrec ! use maximum record length
           nlin=maxbin/20
           npol=maxblk
           counter=6+npol*(nlin+2)
        else
           nlin=nbin/20
           next=nbin-nlin*20
           if (next.gt.0) nlin=nlin+1
           counter=6+npol*(nlin+2)
           recln=counter*80
        endif
        if (opened) then
c
c          Increase the internal record counter
c
           wrecno=wrecno+1
        else
c
c          Reset the internal record counter...
c
           wrecno=1
           open(unit=lun,file=filename(1:lfil),status='unknown',
     &          access='direct',form='unformatted',recl=recln)
           opened=.true.
        endif
      endif
c
c     As from Version 5.0 data is scaled into hex format...
c     This gives an increas in dynamic range by over a factor 6
c
      pvno=5.0
      maxint=65535
c
c     As from Version 6.0, actual resolution of the data is
c     recorded, the factor fcal is no longer recorded, instead
c     the flag fluxflag is written into line 4 to signify whether
c     the data are flux calibrated or not. (1.7.1997)
c
      pvno=6.0
      lct=1
c
c     Line 1
c
      if (reading) then
        read(recrd((lct-1)*80+1:lct*80),11) version, counter, history
        read(version(4:),'(f5.2)') vno
        if (vno.lt.pvno) then
           if (first) then
              write(*,*) 'Reading an old EPN file - Version',vno
              first=.false.
           endif
           maxint=9999
        endif
      else
        vno=pvno
        write(version,'(a3,f5.2)') 'EPN',vno
        write(line(lct),11) version, counter, history
      endif
c
c     Un-comment this line for tracing...
c      write(*,11) version, counter, history
c
 11   format(a8,i4,a68)
c
c     Line 2
c
      lct=lct+1
      if (reading) then
        read(recrd((lct-1)*80+1:lct*80),22) 
     &  jname,cname,pbar,dm,rm,catref,bibref
      else
        write(line(lct),22) jname, cname, pbar, dm, rm, catref, bibref
      endif
c
c     Un-comment this line for tracing...
c      write(*,22) jname, cname, pbar, dm, rm, catref, bibref
c
 22   format(a12,a12,f16.11,f8.3,f10.3,a6,a8,8x)
c
c     Line 3
c
      lct=lct+1
      if (reading) then
        read(recrd((lct-1)*80+1:lct*80),33) 
     &               rah,ram,iras,rras,dsign,ded,dem,des, 
     &               telname, epoch, opos, paflag, timflag
        ras=float(iras)+rras
        if (dsign.eq.'-') ded=-1*ded
        if (dsign.eq.'-'.and.ded.eq.0) then
           dem=-1*dem
           des=-1.0*des
        endif
      else
        rras=ras-int(ras)
        iras=int(ras)
        dsign='+'
        if (index(jname,'-').gt.0) dsign='-'
        if (ded.lt.0) ded=-1*ded
        write(line(lct),33) rah, ram, iras, rras, dsign, ded, dem, des, 
     &               telname, epoch, opos, paflag, timflag
      endif
c
c     Un-comment this line for tracing... 
c      write(*,33) rah, ram, iras, rras, dsign, ded, dem, des, telname, 
c     &           epoch, opos, paflag, timflag
c
 33    format(3i2.2,f4.3,a1,2i2.2,f6.3,a8,f10.3,f8.3,a1,a1,31x)
c
c     Line 4 - (X,Y,Z) of telescope - NEW - implemented in V 6.0
c
      if (vno.ge.6.0) then
        lct=lct+1
        if (reading) then
           read(recrd((lct-1)*80+1:lct*80),44) xtel,ytel,ztel
        else
           write(line(lct),44) xtel,ytel,ztel
        endif
      endif
 44   format(3f17.5,29x) 
c
c     Line 5
c
      lct=lct+1
      if (reading) then
        if (vno.ge.6.0) then
        read(recrd((lct-1)*80+1:lct*80),55) 
     &               cdd,cdm,cdy,scanno,subscan,npol,nfreq, 
     &               nbin, tbin, tres, nint, ncal, lcal, fluxflag
        else
        read(recrd((lct-1)*80+1:lct*80),56) 
     &               cdd,cdm,cdy,scanno,subscan,npol,nfreq, 
     &               nbin, tbin, nint, ncal, lcal, fcal
        endif
      else
        write(line(lct),55) cdd, cdm, cdy, scanno, subscan, npol, nfreq, 
     &               nbin, tbin, tres, nint, ncal, lcal, fluxflag
      endif
c
c     Un-comment this line for tracing...
c      write(*,55) cdd, cdm, cdy, scanno, subscan, npol, nfreq, nbin,
c     &             tbin, tres, nint, ncal, lcal, fluxflag
c      write(*,56) cdd, cdm, cdy, scanno, subscan, npol, nfreq, nbin,
c     &             tbin, nint, ncal, lcal, fcal
c
 55   format(2i2.2,i4.4,2i4.4,i2,i4,i4,2f12.6,i6,i4,i4,a1,15x) ! Version 6
 56   format(2i2.2,i4.4,2i4.4,i2,i4,i4,f12.6,i6,i4,i4,f8.3,20x)
c
c     Line 6 - Blank at the moment
c
      if (vno.ge.6.0) then
        lct=lct+1
        write(line(lct),66) line40,line40
      endif
 66   format(2a40) 
c      
c     Loop around each polarisation...
c
      do ipol=1,npol
c
c       sub-header line 1
c
        lct=lct+1
        if (reading) then
          if (vno.ge.6.0) then
          read(recrd((lct-1)*80+1:lct*80),77) idfield(ipol),nband(ipol),
     &    navg(ipol),f0(ipol),f0u(ipol),df(ipol),dfu(ipol),tstart(ipol)
          else
          f0u(ipol)=' GHz'
          dfu(ipol)=' MHz'
          read(recrd((lct-1)*80+1:lct*80),78) idfield(ipol),nband(ipol),
     &    navg(ipol),f0(ipol),df(ipol),tstart(ipol)
          endif
        else
          if (f0u(ipol).eq.' ') f0u(ipol)=' GHz'  ! Default units
          if (dfu(ipol).eq.' ') dfu(ipol)=' MHz'  ! Default units
          write(line(lct),77) idfield(ipol),nband(ipol),
     &    navg(ipol),f0(ipol),f0u(ipol),df(ipol),dfu(ipol),tstart(ipol)
        endif
c
c       Un-comment this line for tracing...
c       write(*,77) idfield(ipol),nband(ipol),navg(ipol),f0(ipol), 
c     &              f0u(ipol),df(ipol),dfu(ipol),tstart(ipol)
c       write(*,78) idfield(ipol),nband(ipol),navg(ipol),f0(ipol), 
c     &              df(ipol),tstart(ipol)
c
 77     format(a8,i4,i4,f12.8,a8,f12.6,a8,f17.5,7x)
 78     format(a8,i4,i4,f12.8,f12.6,f17.5,23x)
        if (writing) then 
c
c         Zero the integer and character arrays
c
          do i=1,maxbin
            intdata(i)=0
	    chrdata(i)='0000'
          enddo
c
c         Scale the data
c     
          do k=1,nbin
             tmp(k)=rawdata(ipol,k)
          enddo
          call rawtfint(tmp,nbin,intdata,
     &         scale(ipol),offset(ipol),rms(ipol),readwri,maxint)
c
c         Convert to a hex string
c
          do k=1,nbin
             call b102hex(intdata(k),chrdata(k))
          enddo
        endif
c
c       Sub-header line 2
c
        lct=lct+1
        if (reading) then
          read(recrd((lct-1)*80+1:lct*80),88)  scale(ipol),offset(ipol),
     &    rms(ipol),papp(ipol)
        else
          write(line(lct),88) scale(ipol),offset(ipol),
     &                            rms(ipol),papp(ipol)
        endif
c
c       Un-comment this line for tracing...
c        write(*,88) scale(ipol), offset(ipol), rms(ipol), papp(ipol)
c
 88     format(3e12.6,f16.11,28x) 
c
c       Now we finally get to the data!!
c
c
c       Read the data (20 i4s per 80 column line)
c
        nlin=(counter-4)/npol-2
        if (vno.ge.6.0) nlin=(counter-6)/npol-2
        k=1
        do i=1,nlin
          lct=lct+1
          if (reading) then
c            read(recrd((lct-1)*80+1:lct*80),99)(chrdata(j),j=k,k+19)
            l=0
	    do j=k,k+19
              chrdata(j)=recrd((lct-1)*80+1+l:(lct-1)*80+4+l)
              call hex2b10(chrdata(j),intdata(j))
              l=l+4
	    enddo
          else
c            write(line(lct),99)(chrdata(j),j=k,k+19)
            l=1
	    do j=k,k+19
              line(lct)(l:l+3)=chrdata(j)
              l=l+4
	    enddo
          endif
c
c         Un-comment this line for tracing...
c          write(*,99)(chrdata(j),j=k,k+19)
c
          k=k+20
        enddo
c
c       Format for the data
c
 99     format(20a4)
        if (reading) then
c
c         Convert hex string to integer for versions from 5.0
c
          do k=1,nbin
             if (vno.ge.5.0) then
c                call hex2b10(chrdata(k),intdata(k))
             else
                read(chrdata(k),'(i4.4)') intdata(k)
             endif
          enddo
c
c         De-scale the data
c     
          call rawtfint(tmp,nbin,intdata,
     &         scale(ipol),offset(ipol),rms(ipol),readwri,maxint)
          do k=1,nbin
             rawdata(ipol,k)=tmp(k)
          enddo
        endif
      enddo
c
c     Write the record...
c
      if (writing) then
        do i=1,lct
          j=(i-1)*80+1
          k=j+79
          recrd(j:k)=line(i)
        enddo
        write(lun,rec=wrecno) recrd(1:k)
        return
      endif
c
c     Bail out!
c
      return
c
c     On error return this record number
c
 999  recno=-999
      end
c==============================================================================
      subroutine hex2b10(hexstr,b10no)
c==============================================================================

      implicit none

      integer b10no
      character*4 hexstr
      character*1 hex(16)
      integer i,j,dig
      data hex/'0','1','2','3','4','5','6','7','8','9',
     &         'A','B','C','D','E','F'/

      b10no=0
      dig=0

      do i=4,1,-1
         do j=1,16
            if (hexstr(i:i).eq.hex(j)) goto 5
         enddo
 5       continue
         b10no=b10no+int(16.0**float(dig)*float(j-1))
         dig=dig+1
      enddo

      end
c==============================================================================
      subroutine b102hex(b10no,hexstr)
c==============================================================================
c
c	Converts a base 10 number passed down as the integer "b10no"
c	to a hexadeximal string returned as the character*4 "hexstr"
c       N.B. maximum integer for 4 character hex string is 65535
c
      implicit none

      integer b10no
      character*4 hexstr
      character*1 hex(16)
      integer i,j,no,num,dig
      data hex/'0','1','2','3','4','5','6','7','8','9',
     &         'A','B','C','D','E','F'/

      no=b10no
      dig=0

      do i=4,1,-1
         do j=15,0,-1
            num=int(float(j)*16.0**float(i-1))
            if (num.le.no) goto 5
         enddo
 5       continue
         dig=dig+1
         hexstr(dig:dig)=hex(j+1)
         no=no-num
      enddo

      end
c==============================================================================
      subroutine rawtfint(raw,nbin,intdata,scale,offset,rms,dirn,maxint)
c==============================================================================
c
c     Routine to convert "raw" (i.e. floating point) to/from integer
c     format. In this case, the integer may vary between 0 and maxint.
c     Raw data is passed down from site software and stored as integers
c     in the EPN files.
c
c     The direction of the conversion is controlled by "dirn":
c
c     dirn = +1 : raw -> int
c     dirn = -1 : int -> raw
c
c==============================================================================
c    
      implicit none
c
c     Passed down variables
c
      real raw(*)
      real*8 scale, offset, rms
      integer nbin, intdata(*), dirn, maxint
c
c     Local variables
c
      real dmin, dmax, sumsq, pmax, sum, mean, ndiv
      integer i, ibmax, nshift
      if (dirn.eq.1) then
c
c       Scale the data:
c
c       Find minimum & maximum values of the data
c
        dmin = raw(1)
        dmax = raw(1)
        do i=1,nbin
          dmax=max(raw(i),dmax)
          dmin=min(raw(i),dmin)
        enddo
        offset=dmin
        scale=float(maxint)/(dmax-dmin)
c
c       Now do the scaling
c
        do i=1,nbin
          intdata(i)=int((raw(i)-offset)*scale)
        enddo
c
c-------------------------------------------------------------
c       If no rms has been supplied by the user, have a go at 
c       calculating it from the wings of the profile after 
c       shifting it so that its peak lies at bin number nbin/2
c-------------------------------------------------------------
c
        if (rms.eq.0.0) then
          pmax=-1.0e32
          ibmax=0
          do i=1,nbin
             if (raw(i).gt.pmax) then
                ibmax=i
                pmax=raw(i)
             endif
          enddo
          nshift=nbin/2-ibmax
          call sprof(raw,nbin,nshift)
          sumsq=0.0
          sum=0.0
          do i=1,nbin/15
            sumsq=sumsq+raw(i)*raw(i)
            sumsq=sumsq+raw(nbin-i+1)*raw(nbin-i+1)
            sum=sum+raw(i)
            sum=sum+raw(nbin-i+1)
          enddo
	  ndiv=max(1.0,float(2*nbin/15))
          rms=sqrt(sumsq/ndiv)
          mean=sum/ndiv
        endif
      else if (dirn.eq.-1) then
c
c       Prepare to de-scale the data:
c
        do i=1,nbin
          raw(i)=real(offset+dble(intdata(i))/scale)
        enddo
      else
c
c       Silly option was passed down...
c
        write(*,*) 'Invalid value for dirn passed to rawtfint : ',dirn
        stop
      endif
c
c     Job done!
c
      end
c==============================================================================
      integer function nepnrec(filename)
c==============================================================================
c
c     Function to find out the number of records in an EPN file.
c
c     This can in principle be done with an INQUIRE statement but
c     it is system dependent. The following approach uses a simple
c     binary search algorithm to find the last record which is far
c     quicker than reading in every record in very large EPN files.
c
c     if the file doesn't exist nepnrec is returned as 0
c     nepnrec = -1 if an error occurred whilst reading the file
c
c     Created 96/12/03 DRL@MPIfR
c     Modified October 1997 to return -1 on error
c
c
      implicit none
      character*(*) filename

      character lin1*12,recrd*80
      integer irec,istat,counter,recln,hirec,lorec,orec,lun
      logical filex
c
c     Check to see whether the EPN file exists...
c     return nepnrec=0 if it doesn't
c
      nepnrec=0
      inquire(file=filename,exist=filex)
      if (.not.filex) return 
c
c     It exists! Open it and find out the true record length..
c
      call glun(lun)
      nepnrec=-1 ! This will be the value of nepnrec if an error occurs
      irec=-1
      open(unit=lun,file=filename,status='old',access='direct',
     &     form='unformatted',recl=12,err=999)
      read(lun,rec=1,err=999)lin1
      read(lin1,'(8x,i4)',err=999)counter
      close(unit=lun)
      recln=counter*80
      if (recln.le.0) return ! Error
c
c     Find out crude lower and upper bounds for the last record
c
      open(unit=lun,file=filename,status='old',access='direct',
     &     form='unformatted',recl=recln,err=999)
      irec=1
      istat=0
      do while(.true.)
         read(lun,rec=irec,iostat=istat) recrd
         if (istat.ne.0) goto 998
         irec=irec*2
      enddo
 998  hirec=irec
      lorec=irec/2
c
c     Now search for the last record within these bounds iteratively
c
      istat=0
      orec=0
      do while(orec.ne.irec)
         orec=irec
         irec=lorec+(hirec-lorec)/2
         read(lun,rec=irec,iostat=istat) recrd
         if (istat.eq.0) then
            lorec=irec
         else
            hirec=irec
         endif
      enddo
 999  continue
      close(unit=lun)
      nepnrec=irec
      end
c============================================================================= 
      subroutine dattim(date,hh,mm,ss)
c============================================================================= 
c
c     Returns the time obtained by truncating the MJD passed down
c     as date
c
      implicit none 
      double precision date,dummy,ss
      integer hh,mm


      dummy=abs(date-int(date))

      dummy=dummy*24.0
      hh=int(dummy)

      dummy=abs(dummy-int(dummy))

      dummy=dummy*60.0
      mm=int(dummy)

      dummy=abs(dummy-int(dummy))

      dummy=dummy*60.0
      ss=dummy

      end
      subroutine smooth(pr,nbin,kwmax,snrmax,smmax)
c******************************************************************
c
c  convolves profile pr(nbin) with a boxcar of width kw.  it returns
c    the width kwmax which gave the highest s/n ratio snrmax, and the
c    corresponding pulse amplitude smmax.
c
c******************************************************************
c

      integer nbin,kwmax
      real*4 pr(*),rmsp,snrmax,smmax
c
      integer j,k,kw,nn,ja,jj
      real*4 s,wrk(43),al,an,sn,smax
c
      snrmax=0.
c---------------------------------------
c  remove baseline
c      ksm=nbin/2.5+0.5
c      smax=1.e30
c      do 10 j = 1,nbin
c        s=0.0
c        do 20 k = 1,ksm
c          s = s + pr(mod(j+k-1,nbin)+1)
c   20   continue
c        if(s.lt.smax) smax=s
c   10 continue
c      smax=smax/ksm
c      do 30 j = 1,nbin
c        pr(j) = pr(j) - smax
c   30 continue
c--------------------------------------
c      remove baseline and calc rmsp
      do i=1,2
      s=0
      k=0
      do j=1,nbin
         if (j.le.nbin/10.or.j.ge.nbin-nbin/10) then
            k=k+1
            s=s+pr(j)
         endif
      enddo
      s=s/real(k)
      do j=1,nbin
         pr(j)=pr(j)-s
      enddo
      enddo
      k=0
      rmsp=0.0
      do j=1,nbin
         if (j.le.nbin/10.or.j.ge.nbin-nbin/10) then
            k=k+1
            rmsp=rmsp+pr(j)*pr(j)
         endif
      enddo
      rmsp=sqrt(rmsp/real(k))
c
c
      do 40 nn=1,6
        kw=2**(nn-1)
        if(kw.gt.nbin/2) return
	s=0.0
	do 50 k=1,kw
	  s=s+pr(k)
	  wrk(k)=pr(k)
   50   continue
	ja=0
	smax=s
	do 60 j=2,nbin
	  ja=ja+1
	  if(ja.gt.kw) ja=ja-kw
	  al=wrk(ja)
	  jj=j+kw-1
	  if(jj.gt.nbin)jj=jj-nbin
	  an=pr(jj)
	  s=s+an-al
	  wrk(ja)=an
	  if(s.gt.smax) smax=s
   60   continue

        sn=smax/(rmsp*sqrt(kw*(1.+real(kw)/nbin)))
        if(sn.gt.snrmax) then
          snrmax=sn
          kwmax=kw
          smmax=smax/kw
        endif
   40 continue

      end
c==============================================================================
      subroutine sprof(profile,nbins,shift)
c==============================================================================
c
c     shifts the profile in the array profile() with nbins bins by the 
c     number of bins passed down in shift. shift>0 means shift forward,
c     shift<0 means shift backwards.
c
      implicit none

      real profile(*)
      integer nbins, shift
c
c     local variables
c
      integer i,j
      real dummy

      if (shift.lt.0) shift=nbins+shift
      do i=1,shift
         dummy=profile(nbins)
         do j=nbins,2,-1
            profile(j)=profile(j-1)
         enddo
         profile(1)=dummy
      enddo

      end
      subroutine grayscale(dat,nxd,nyd,nx,ny,flip)

      integer nxd,nyd
      integer nsym,flip
      parameter(nsym=10)
      real*4 dat(nxd,nyd),tr(6),xw,yw
      save
      
      s=0.
      ss=0.
      do j=1,ny
	 smin=1.0e30
         do i=1,nx
	    aa=dat(i,j)
	    smin=min(smin,aa)
	    s=s+aa
         enddo
	 do i=1,nx
	    dat(i,j)=dat(i,j)-smin
	 enddo
      enddo
      s=s/(nx*ny) 
      smin=1.e30
      smax=-smin
      do j=1,ny
         do i=1,nx
            aa=dat(i,j)-s
	    ss=ss+aa**2
	    smax=max(smax,dat(i,j))
	    smin=min(smin,dat(i,j))
	    smin=min(smin,dat(i,j))
         enddo
      enddo
      rms=sqrt(ss/(nx*ny))
      
      tr(1)=0.0
      tr(2)=1.0
      tr(3)=0.0
      tr(4)=0.0
      tr(5)=0.0
      tr(6)=1.0
      xmin=xw(tr,0.5,0.5)
      ymin=yw(tr,0.5,0.5)
      xmax=xw(tr,real(nx)+0.5,real(ny)+0.5)
      ymax=yw(tr,real(nx)+0.5,real(ny)+0.5)

      if (flip.eq.1) then
         call pgwindow(xmin,xmax,ymax,ymin)
      else
         call pgwindow(xmin,xmax,ymin,ymax)
      endif
c      call pggray(dat,nxd,nyd,1,nx,1,ny,smax,s+0.5*rms,tr)
      call pggray(dat,nxd,nyd,1,nx,1,ny,smax,smin,tr)
      
      end


      real function xw(tr,i,j)
      real tr(6),i,j
      xw=tr(1)+tr(2)*i+tr(3)*j
      end
      
      real function yw(tr,i,j)
      real tr(6),i,j
      yw=tr(4)+tr(5)*i+tr(6)*j
      end
      
C from numerical recipes
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      if (n.eq.1) then
	indx(1)=1
        return
      endif
      DO 11 J=1,N
         INDX(J)=J
 11   CONTINUE
      L=N/2+1
      IR=N
 10   CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         INDXT=INDX(L)
         Q=ARRIN(INDXT)
      ELSE
         INDXT=INDX(IR)
         Q=ARRIN(INDXT)
         INDX(IR)=INDX(1)
         IR=IR-1
         IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
 20   IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
         ENDIF
         IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         ENDIF
         GO TO 20
      ENDIF
      INDX(I)=INDXT
      GO TO 10
      END
c     DECK LENGTH
c     
c     
c     
c     
c     RETURNS THE LENGTH OF 'STRING' EXCLUDING ANY TRAILING SPACES.
c     
      integer function length(string)
      implicit none
      character string*(*)
c     
c     OBTAIN THE LOCATION OF THE LAST NON-SPACE CHARACTER.
c     
      integer ilen,i
c     search for the first null
      ilen = len(string)
c     use the position of the first null
      do 1 i = ilen, 1, -1
c     
c     LENGTH FOUND.
c     
         if (string(i:i) .ne. char(32) .and.
     &       string(i:i).ne.char(0)) then
            length = i
            return 
         end if
c     
c     STRING IS ALL SPACES OR ZERO LENGTH.
c     
    1 continue
      length = 0
c     
c     END OF INTEGER FUNCTION LENGTH.
c     
      return 
      end










