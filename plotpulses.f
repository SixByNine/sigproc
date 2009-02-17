c==============================================================================
      program plotpulses
c==============================================================================
      implicit none
      character*80 stem(99),line,pgdev,filename,ch*1,command*240
      integer lun,lst(99),npulses,width,idx,fft,maxp,i,narg,sym,
     &	      nsnbins,ndms,maxdmhist,mindmhist,dmplot,idm,j,
     &        n_this_dm, nmax, nfiles, istat, pgband, pgcurs
      logical filexists,interactive
      parameter(maxp=10000000)
      real tsamp,dm(maxp),sn(maxp),mindm,maxdm,minsn,maxsn,mints,maxts,
     &     time(maxp),thresh,snmid(10000),maxsnhist,snhist(10000),yref,
     &     maxdmhst,dmhist(10000),dmval(10000),peakdm,width_ms,xref,
     &	   tmin,tmax,wmin,wmax,maxsnr,rndms,dmlast,dmmin,dmmax,x,y
c==============================================================================
      narg=iargc()
      if (narg.lt.1) then
	 write(*,*)
	 write(*,*) 'plotpulses - PGPLOT results of single pulse search'
	 write(*,*)
         write(*,*) 'usage: plotpulses filestem1 .... filestemn',
     &              ' -s1 minsn -s2 maxsn',
     & 		    ' -w1 wmin -w2 wmax -d1 dmmin -d2 dmmax', 
     &		    ' -t1 tmin -t2 tmax -p pgdev -d dmplot',
     &		    ' -n nmax -i'
	 write(*,*) ''
	 write(*,*) 'minsn - minimum s/n to display (def=5)'
	 write(*,*) 'maxsn - maximum s/n to display (def=all)'
	 write(*,*) 'pgdev - pgplot device (def=/xs)'
	 write(*,*) 'wmin - minimum width (in ms) to display (def=all)'
	 write(*,*) 'wmax - maximum width (in ms) to display (def=all)'
	 write(*,*) 'dmmin - minimum DM (def=all)'
	 write(*,*) 'dmmax - maximum DM (def=all)'
         write(*,*) 'tmin - minimum time (seconds) to display',
     &		    ' (def=all)'
	 write(*,*) 'tmax - maximum time (seconds) to display',
     &              ' (def=all)'
	 write(*,*) 'dmplot - 1 to plot DM and 2 to plot DM channel (def=1)'
         write(*,*) 'nmax - maximum number of events to plot per DM',
     &		    ' channel (def=all)'
         write(*,*) 'allbeams - will do a DM-t stack for all beams',
     &              ' in the current directory'
	 write(*,*)
         stop
      endif
      lun=20
      nfiles=0
      minsn=5
      dmplot=1
      idm=0
      maxsnr=10000
      tmin=-1e32
      tmax=1e32
      wmin=0
      wmax=100000
      nmax=100000
      dmmin = -1e32
      dmmax = 1e32
      pgdev='/xs'
      interactive=.false.
      i=1
      do while (i.le.narg)
         call getarg(i,line)
         inquire(file=line(1:index(line,' ')-1)//'.pls',exist=filexists)
         if (filexists) then
            nfiles=nfiles+1
            call getarg(i,stem(nfiles))
            lst(nfiles)=index(stem(nfiles),' ')-1
         elseif (line.eq.'allbeams') then
            istat=system("ls *.pls|awk -F. '{print $1}' > allbeams.lis")
            open(lun,file='allbeams.lis',iostat=istat)
            do while (istat.eq.0)
               nfiles=nfiles+1
               read(lun,'(a)',iostat=istat) stem(nfiles)
               lst(nfiles)=index(stem(nfiles),' ')-1
            enddo
            if (nfiles.eq.0) stop 'no .pls files found...'
            close(lun)
            nfiles=nfiles-1
            istat=system("rm allbeams.lis")
         elseif (line.eq.'-s1') then
            i=i+1
            call getarg(i,line)
            read(line,*) minsn
         elseif (line.eq.'-s2') then
            i=i+1
            call getarg(i,line)
            read (line,*) maxsnr
         elseif (line.eq.'-p') then
            i=i+1
            call getarg(i,pgdev)
         elseif (line.eq.'-w1') then
            i=i+1
            call getarg(i,line)
            read(line,*) wmin
         elseif (line.eq.'-w2') then
            i=i+1
            call getarg(i,line)
            read(line,*) wmax
         elseif (line.eq.'-t1') then
            i=i+1
            call getarg(i,line)
            read(line,*) tmin
         elseif (line.eq.'-t2') then
            i=i+1
            call getarg(i,line)
            read(line,*) tmax
         elseif (line.eq.'-d1') then
            i=i+1
            call getarg(i,line)
            read(line,*) dmmin
         elseif (line.eq.'-d2') then
            i=i+1
            call getarg(i,line)
            read(line,*) dmmax
         elseif (line.eq.'-d') then
            i=i+1
            call getarg(i,line)
            read(line,*) dmplot
         elseif (line.eq.'-n') then
            i=i+1
            call getarg(i,line)
            read(line,*) nmax
         elseif (line.eq.'-i') then
            interactive=.true.
         endif
         i=i+1
      enddo 

c==============================================================================
      call pgbegin(0,pgdev,1,nfiles)
      do j=1,nfiles
         call pgadvance
         i=0
         mints=1.0e32
         maxts=-1.0e32
         mindm=1.0e32
         maxdm=-1.0e32
         maxsn=0.0
         n_this_dm = 0
         filename=stem(j)
         open(lun,file=filename(1:lst(j))//'.pls',status='unknown')
         read(lun,'(a)') line
         read(line(7:),*) tsamp
         read(line(24:),*) thresh
         do while(.true.)
            i=i+1
            read(lun,*,err=1,end=1) dm(i),width,idx,sn(i),fft
            if (i.eq.1) dmlast=dm(i)
            n_this_dm = n_this_dm + 1
            if (dm(i).ne.dmlast) then
               idm=idm+1
               n_this_dm = 1
            endif
            dmlast=dm(i)
            if (dmplot.eq.2) dm(i)=idm
            width_ms=(tsamp*2.**width)*1.e-3
            time(i)=idx*tsamp/1.0e6
c            if (sn(i).ge.minsn) then
            if (width_ms.le.wmax.and.width_ms.ge.wmin.and.
     &           sn(i).ge.minsn.and.sn(i).le.maxsnr.
     &           and.time(i).ge.tmin.and.time(i).le.tmax.
     &           and.n_this_dm.le.nmax.
     &           and.dm(i).ge.dmmin.and.dm(i).le.dmmax) then
               mints=min(time(i),mints)
               maxts=max(time(i),maxts)
               mindm=min(dm(i),mindm)
               maxdm=max(dm(i),maxdm)
               maxsn=max(sn(i),maxsn)
               snhist(sn(i)-thresh)=snhist(sn(i)-thresh)+1
            else
               i=i-1
            endif
         enddo
 1       close(lun)
         if (dmmin.ne.-1.0e32) mindm=dmmin
         if (dmmax.ne.+1.0e32) maxdm=dmmax
         write(*,*) dmmin,dmmax,mindm,maxdm
         npulses=i-1
         nsnbins = maxsn-thresh
         maxsnhist = 0
         do i = 1, nsnbins
            snmid(i) = thresh+(i-1)+0.5
            if (snhist(i).ne.0) then
               snhist(i) = log10(snhist(i))
            else
               snhist(i) =-0.5
            endif
            if (snhist(i).gt.maxsnhist) maxsnhist = snhist(i)
         end do
         
         filename=stem(j)
         open(lun,file=filename(1:lst(j))//'.hst',status='unknown')
         
         i = 0
         do while(.true.)
            i=i+1
            read (lun,*,err=1,end=2) dmval(i), dmhist(i)
            if (dmplot.eq.2) dmval(i) = i
         enddo
 2       close(lun)
         
         ndms = i-1
         maxdmhist = 0
         mindmhist = 1000000
         do i = 1, ndms
            if (dmhist(i).gt.maxdmhist) then
               maxdmhist = dmhist(i)
               peakdm = dmval(i)
            endif
            if (dmhist(i).lt.mindmhist) then
               mindmhist = dmhist(i)
            endif
         end do
         write(line,100) npulses,maxsn,peakdm
 100     format(i8,' events. max S/N =',f5.1,'. Peak DM =',f7.1,
     &        ' cm\\u-3\\d pc')
c==============================================================================
         call pgscf(2)
         call pgswin(0.0,1.0,0.0,1.0)
         filename=stem(j)
         if (nfiles.eq.1) then
c            goto 888
            call pgtext(-0.02,1.0,'File: '//filename(1:lst(j))
     &           //' '//line)
            call pgsvp(0.06, 0.31, 0.6, 0.87)
            call pgswin(minsn,maxsn,log10(0.5),1.1*maxsnhist)
            call pgsch(0.8)
            call pgbox('bcnst',0.0,0,'bclnst',0.0,0)
            call pgmtxt('B', 2.5, 0.5, 0.5, 'Signal-to-noise ratio')
            call pgmtxt('L', 1.8, 0.5, 0.5, 'Number of Pulses')
            call pgsch(1.0)
            call pgbin(nsnbins,snmid,snhist,1)
            call pgsvp(0.39, 0.64, 0.6, 0.87)
            if (dmplot.eq.1) then
               call pgswin(mindm-0.5,maxdm+0.5,mindmhist-0.1*maxdmhist,
     &              maxdmhist+0.1*maxdmhist)
            else
               rndms=ndms
               call pgswin(0.0,rndms,mindmhist-0.1*maxdmhist,
     &              maxdmhist+0.1*maxdmhist)
            endif
            call pgsch(0.8)
            call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
            if (dmplot.eq.1) then
               call pgmtxt('B', 2.5, 0.5, 0.5, 'DM (cm\\u-3\\d pc)')
            else
               call pgmtxt('B',2.5,0.5,0.5,'DM Channel')
            endif
            call pgmtxt('L', 1.8, 0.5, 0.5, 'Number of Pulses')
            call pgsch(1.0)
            call pgbin(ndms,dmval,dmhist,1)
            call pgsvp(0.72, 0.97, 0.6, 0.87)
            if (dmplot.eq.1) then
               call pgswin(mindm-0.5,maxdm+0.5,minsn,maxsn)
            else
               call pgswin(0.0,rndms,minsn,maxsn)
            endif
            call pgsch(0.8)
            call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
            if (dmplot.eq.1) then
               call pgmtxt('B', 2.5, 0.5, 0.5, 'DM (cm\\u-3\\d pc)')
            else
               call pgmtxt('B',2.5,0.5,0.5,'DM Channel')
            endif
            call pgmtxt('L', 1.8, 0.5, 0.5, 'Signal-to-noise ratio')
            call pgpt(npulses,dm,sn,20)
         endif
 888     if (tmin.ne.-1.0e32) mints=tmin
         if (tmax.ne.+1.0e32) maxts=tmax
         if (dmplot.eq.1) then
            call pgswin(mints,maxts,mindm,maxdm)
         else
            call pgswin(mints,maxts,0.0,rndms)
         endif
         if (nfiles.eq.1) then
            call pgsvp(0.06, 0.97, 0.08, 0.52)
            call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
         else
            call pgsvp(0.1,1.0,0.1,0.9)
            call pgsch(1.5)
            call pgbox('bnstc',0.0,0,'bvcnst',0.0,0)
         endif
c         call pgsch(0.8)
         if (nfiles.eq.1) then
            call pgmtxt('B',2.5,0.5,0.5,'Time (s)')
            if (dmplot.eq.1) then
               call pgmtxt('L',1.8,0.5,0.5,'DM (cm\\u-3\\d pc)')
            else
               call pgmtxt('L',2.5,0.5,0.5,'DM Channel')
            endif
         endif
         do i=1,npulses
            sym = (sn(i)-minsn)*0.5+20+0.5
            if (sym.gt.25) sym = 25
            if (sn(i).gt.minsn) call pgpt1(time(i),dm(i),sym)
         enddo
         do while (interactive)
            i=pgcurs(xref,yref,ch)
            if (ch.eq.'X') stop
            write(command,'(a,f9.3,f12.3,a)') 
     &        'extract '//stem(j)(1:index(stem(j),' ')-1)//
     &                 ' ',xref,yref,'> /dev/null'
            i=system(command)
         enddo
c==============================================================================
      enddo
      call pgend
      end
c==============================================================================
