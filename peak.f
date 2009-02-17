c=============================================================================
      program peak 
c=============================================================================
c
c     A program to find giant pulses in a noisy time series.  
c
c     Created: 98/11/28 (dunc@naic.edu)
c     
c     Modification history:
c
c     05/11/14 - (dunc) ressurected for use within SIGPROC
c
c=============================================================================
      implicit none
      include 'seek.inc'
      include 'peak.inc'
      integer i,llog                         ! lun for file in quiet mode
      call peakin(llog)                      ! command line inputs
      call timstart(llog)                    ! fire up the ship's clock
      call readdat(llog,pzero)                     ! read in the time series
      call baseline(llog)    ! normalise time series
      do i=1,nloop
      call findpeak(llog,i)                  ! find peaks above snrmin
      if (pout.or.xout) call plotpeak(llog,i)! make summary plot if required
      if (i.lt.nloop) call halve(series,ntim)! halve time resolution if repeat
      tsamp=tsamp*2.0
      enddo
      call timfinis(llog)                    ! stop the clock
      end
c=============================================================================
c==============================================================================
      subroutine findpeak(llog,iloop)
c==============================================================================
      implicit none
      include 'seek.inc'
      include 'peak.inc'
      integer llog,i,j,k,l,m,n,o,p,nstart,nfinis,last,iloop
      real prfsnr,pulse(nbins),rav

      write(llog,*) 'Searching for peaks...',tsamp,ntim
      nfound=0
      last=-8
      l=ntim/nsave
      m=0
      n=0
      rav=0.0
      if (wepn) call system('rm -f peak.epn')
      i=0
      do while(i.lt.ntim)
         i=i+1
         if (pout.or.xout) then
            m=m+1
            rav=rav+series(i)
            if (m.eq.l) then
               n=n+1
               smoothed(n)=rav/real(l)
               m=0
               rav=0.0
            endif
         endif
         if (series(i).gt.snrmin.and.(i-last).gt.8) then
            nfound=nfound+1
            p=i
            do o=i,i+8
               if (series(o).gt.series(p)) then
                  p=o
               endif
            enddo
            i=p
            if (nfound.le.nsave) then
               pulidx(nfound)=i
               pulsnr(nfound)=series(i)
            endif
            last=i
            if (wepn) then
               nstart=max(1,i-nbins/2)
               nfinis=nstart+nbins-1
               k=0
               do j=nstart,nfinis
                  k=k+1
                  pulse(k)=series(j)
               enddo
               call writeepn(pulse,nbins,1.0,refdm,1.0,1.0,
     &         tsamp,0.0,tsamp*real(i)*1.0e6,'peak.epn',prfsnr,' pulse')
            endif
         endif
      enddo

      write(llog,*) 'Found',nfound,' peaks with S/N >',snrmin
      call glun(i)
      if (iloop.gt.1) facc='append'
      open(unit=i,file='peak.sum',status='unknown',access=facc)
      write(i,*) refdm,int(tsamp*1.0e6),nfound,snrmin
      close(unit=i)

      end
c==============================================================================
c=============================================================================
      subroutine plotpeak(llog,iloop) 
c=============================================================================
      implicit none
      include 'seek.inc'
      include 'peak.inc'
      integer llog,i,j,k,l,m,idx(nsave),lun,ndisp,iloop,t,lout
      real xmin,xmax,ymin,ymax,rtsamp,const,delt,frac,dm(nsave)
      logical lexist
      character*80 rank,csnr,text,user*8,comment,pwd,line

c      call system('rm -f pwd')
c      call system('/usr/bin/pwd > pwd')
c      open(unit=10,file='pwd',status='old')
c      read(10,'(a)') pwd
c      close(unit=10)
c      call system('rm -f pwd')

      call getenv('PWD',pwd)
      comment='Comments: __________________________________'
      ndisp=5
      if (nfound.lt.1) return
      if (ndisp.gt.nfound) ndisp=nfound
      
      rtsamp=real(tsamp)

      if (nloop.eq.1) then
         text='peak.ps/ps'
      else
         write(text,'(''peak'',i1,''.ps/ps'')') iloop
      endif 
      if (pout) then
         call pgbegin(0,text,1,1)
      else
         call pgbegin(0,'/xs',1,1)
      endif
      call pgvport(0.1,0.9,0.1,0.25)
      call pgscf(2)
      xmin=rtsamp
      xmax=real(ntim)*rtsamp


      ymin=1.e32
      ymax=-ymin
c      do i=1,nsave
c         ymin=min(ymin,smoothed(i))
c         ymax=max(ymax,smoothed(i))
c      enddo
c      do i=1,ntim
c         ymin=min(ymin,series(i))
c         ymax=max(ymax,series(i))
c      enddo
	ymin=-1.5
	ymax=+1.5
      call pgwindow(xmin,xmax,ymin,ymax*1.1)
      call pgbox('bcnst',0.0,0,'bc',0.0,0)
      call pglabel('Time (s)',' ',' ')

      const=rtsamp*real(ntim)/nsave
      call pgmove(xmin,smoothed(1))
      do i=2,nsave
         call pgdraw(real(i)*const,smoothed(i))
      enddo

      call pgsls(2)
      call pgmove(xmin,0.0)
      call pgdraw(xmax,0.0)
      call pgsls(1)

      call pgvport(0.1,0.9,0.275,0.425)
      ymin=1.e32
      ymax=-ymin
      do i=1,ntim
         ymin=min(ymin,series(i))
         ymax=max(ymax,series(i))
      enddo
      ymin=snrmin-2

      call pgbox('bc',0.0,0,'bc',0.0,0)
      call pgwindow(xmin,xmax,ymin,ymax*1.1)
      const=rtsamp
      call pgmove(xmin,series(1))
      do i=2,ntim
         call pgdraw(real(i)*const,series(i))
      enddo

      call pgsls(2)
      call pgmove(xmin,snrmin)
      call pgdraw(xmax,snrmin)
      call pgsls(1)

c      do i=1,nsave
c         ymin=min(ymin,smoothed(i))
c         ymax=max(ymax,smoothed(i))
c      enddo
c
c      call pgwindow(xmin,xmax,ymin,ymax*1.1)
c
c      const=rtsamp*real(ntim)/nsave
c      call pgmove(xmin,smoothed(1))
c      do i=2,nsave
c         call pgdraw(real(i)*const,smoothed(i))
c      enddo

      call pgsch(0.7)
      if (nfound.gt.1) then
         call indexx(nfound,pulsnr,idx)
      else
         idx(1)=1
      endif
      j=0
      do i=nfound,nfound-(ndisp-1),-1
         j=j+1
         write(rank,'(a,i1,a)') '(',j,')'
         call pgtext(real(pulidx(idx(i)))*rtsamp,pulsnr(idx(i)),rank)
      enddo

	call glun(lout)
	open(lout,file="peak.out",status="unknown")
      do i=nfound,1,-1
         write(lout,*) real(pulidx(idx(i)))*rtsamp+skp,pulsnr(idx(i))
      enddo
	close(lout)

      call pgsch(0.6)
      j=0
      delt=0.8/real(ndisp)
      do i=nfound,nfound-(ndisp-1),-1
         j=j+1
         call pgvport(0.1+delt*real(j-1),0.1+delt*real(j),0.45,0.6)
         call pgwindow(0.0,1.0,0.0,1.0)
         call pgbox('bc',0.0,0,'bc',0.0,0)
         write(rank,'(a,i1,a)') '(',j,')'
         call pgtext(0.05,0.9,rank)
         write(csnr,'(a,f4.1)') 'S/N = ',pulsnr(idx(i))
         call pgtext(0.2,0.9,csnr)
         k=max(1,pulidx(idx(i))-nbins/2)
         if (j.eq.1) then
            xmin=1.0
            xmax=real(nbins)
            ymin=1.e32
            ymax=-ymin
            do l=k,k+nbins-1
               ymin=min(ymin,series(l))
               ymax=max(ymax,series(l))
            enddo
         endif
         call pgwindow(xmin,xmax,ymin,ymax*1.15)
         call pgmove(xmin,series(k))
         m=1
         do l=k+1,k+nbins-1
            m=m+1
            call pgdraw(real(m),series(l))
         enddo
         call pgsls(2)
         call pgmove(xmin,0.0)
         call pgdraw(xmax,0.0)
         call pgsls(1)

         
      enddo

      call pgvport(0.6,0.9,0.675,0.9)
      inquire(file='pkdm.sum',exist=lexist)
      
      if (lexist) then
         call glun(lun)
         open(unit=lun,file='pkdm.sum',status='old')
         xmin=1.e32
         xmax=-xmin
         ymin=1.e32
         ymax=-ymin
         t=nint(tsamp*1.0e6)
         do i=1,nsave
            read(lun,'(a)',end=1) line
            if (iloop.eq.1) read(line,*) dm(i),idx(i)
            if (iloop.eq.2) read(line,*) dm(i),j,idx(i)
            if (iloop.eq.3) read(line,*) dm(i),j,j,idx(i)
            if (iloop.eq.4) read(line,*) dm(i),j,j,j,idx(i)
            if (iloop.eq.5) read(line,*) dm(i),j,j,j,j,idx(i)
            xmin=min(xmin,dm(i))
            xmax=max(xmax,dm(i))
            ymin=min(ymin,real(idx(i)))
            ymax=max(ymax,real(idx(i)))
         enddo
 1       j=i-1
         close(unit=lun)
         call pgwindow(xmin,xmax,ymin,ymax*1.1)
         call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
         call pglabel('DM (cm\\u-3\\d pc)','N\\dTOT\\u',' ')
         call pgmove(dm(1),real(idx(1)))
         call pgpoint(1,dm(1),real(idx(1)),17)
         do i=2,j
            call pgdraw(dm(i),real(idx(i)))
            call pgpoint(1,dm(i),real(idx(i)),17)
         enddo
      else
         ymin=0.0
         ymax=1.1
         xmin=snrmin
         xmax=pulsnr(idx(nfound))
         call pgwindow(xmin*0.9,xmax*1.1,ymin,ymax)
         call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
         call pglabel('S/N of detected pulse','N(>S/N)/N\\dTOT\\u',' ')
         call pgsls(2)
         call pgmove(snrmin,ymin)
         call pgdraw(snrmin,ymax)
         call pgsls(1)
         k=nfound
         do i=1,nfound
            frac=real(k)/real(nfound)
            call pgpoint(1,pulsnr(idx(i)),frac,17)
            k=k-1
         enddo
      endif
c
c     Header info...
c
      call pgvport(0.1,0.5,0.5,0.9)
      call pgwindow(0.0,1.0,0.0,1.0)
c      call pgscf(4)
      call pgsch(1.5)
      call pgtext(0.0,0.9,'PEAK... Single Pulse Search')
      call pgsch(0.8)
      call pgscf(2)
c      call pgtext(0.0,0.8,'File: '//filename(1:lst+1))
      call pgtext(0.0,0.8,'File: '//pwd)
      write(text,'(a,i4,6x,a,f6.1,a)')
     &    'N\\dTOT\\u: ',nfound,
     &    'DM: ',refdm,' cm\\u-3\\d pc'
      call pgtext(0.0,0.7,text)

      text=' '
      call gruser(user,l)
      call grdate(text,m)

      call pgtext(0.0,0.5,comment)
      call pgscf(3)
      call pgtext(0.0,0.35,'Hunter ID: '//user//' Date:  '//text)

      call pgend
      
      end 
c=============================================================================


      subroutine halve(data,npts)
      implicit none
      real data(*)
      integer i,j,npts
      real sum
      sum=0.0
      j=0
      do i=1,npts
         sum=sum+data(i)
         if (mod(i,2).eq.0) then
	    j=j+1
            data(j)=sum/2.0
            sum=0.0
         endif
      enddo
      npts=npts/2
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
	write(*,*) 'hello there'
      END
c=============================================================================
      subroutine peakin(llog) 
c=============================================================================
c
c   Controls the command-line inputs to peak
c
      implicit none
      include 'seek.inc'
      include 'peak.inc'
      integer narg,iargc,i,llog,p2
      character*80 option

      llog=6
      wepn=.false.
      pout=.false.
      xout=.false.
      snrmin=5.0
      nloop=1
      rfac=1
      tsize=19
      p2=nint(log10(real(npts))/log10(2.0))
      dmidx=-1
      facc='sequential'
      skp=0

      narg=iargc()
 1    format(a)
      if (narg.lt.1) then
         write(*,*)
         write(*,1)'peak - searches for giant pulses in time series'
         write(*,*)
         write(*,1)'usage: peak <INFILE> -{options}'
         write(*,*)
         write(*,1)'The input file may be a time series, or a set of'
         write(*,1)'set of filterbank channels. The file suffix MUST,'
         write(*,1)'however, be either ".tim", ".ser" or ".dis".'
         write(*,*)
         write(*,1)'options:'
         write(*,*)
         write(*,1)'-A: append output ASCII files if already exist'
         write(*,1)'-q: quiet mode - all standard messages > peak.log'
         write(*,1)'-e: write an EPN file "peak.epn" with the pulses'
         write(*,1)'-p: write a PostScript plot summarising the search'
         write(*,1)'-X: write an X-display plot summarising the search'
         write(*,*)
         write(*,1)'-l[loop]: analysis loop times, doubling tsamp'
         write(*,1)'-s[smin]: set S/N threshold to smin (def=5.0)'
         write(*,1)'-t[tlen]: fix time series length to 2**tlen points'
         write(*,1)'-i[tsec]: ignore tsec seconds of data on reading'
         write(*,*)
         stop
      endif

      call getarg(1,filename)
      
      if (index(filename,'.ser').gt.0) then
         lst=index(filename,'.ser')-1
      else if (index(filename,'.tim').gt.0) then
         lst=index(filename,'.tim')-1
      else if (index(filename,'.dat').gt.0) then
         lst=index(filename,'.dat')-1
      else if (index(filename,'.dis').gt.0) then
         lst=index(filename,'.dis')-1
         dmidx=1
      else
         stop 'file type not recognized! Type peak for help.'
      endif
      
      do i=2,narg
        call getarg(i,option)
        if (index(option,'-s').gt.0) then
          read(option(3:),*) snrmin
        else if (index(option,'-t').gt.0) then
           read(option(3:),*) tsize
        else if (index(option,'-i').gt.0) then
          read(option(3:),*) skp
        else if (index(option,'-l').gt.0) then
          read(option(3:),*) nloop
        else if (index(option,'-q').gt.0) then
          call glun(llog)
        else if (index(option,'-e').gt.0) then
          wepn=.true.
        else if (index(option,'-p').gt.0) then
          pout=.true.
        else if (index(option,'-X').gt.0) then
          xout=.true.
        else if (index(option,'-A').gt.0) then
          facc='append'
        else
          write(*,*) 'WARNING.. command line option ',
     &    option(1:index(option,' ')-1),' not recognized!!!'
        endif
      enddo
      
      if (llog.ne.6) open(unit=llog,file='peak.log',status='unknown',
     &    access=facc)
      write(llog,*)
      write(llog,*) 'PEAK: ',version
      end 
c=============================================================================
