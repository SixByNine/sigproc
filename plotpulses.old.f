c==============================================================================
      program plotpulses
c==============================================================================
      implicit none
      character*80 stem,line
      integer lun,lst,npulses,width,idx,fft,maxp,i,narg,sym
      parameter(maxp=1000000)
      real tsamp,dm(maxp),sn(maxp),mindm,maxdm,minsn,maxsn,mints,maxts,
     &     time(maxp)
c==============================================================================
      narg=iargc()
      if (narg.lt.2) then
         write(*,*) 'usage: plotpulses filestem minsn'
         stop
      endif
      call getarg(1,stem)
      lst=index(stem,' ')-1
      call getarg(2,line)
      read(line,*) minsn
      i=0
      lun=20
      mints=1.0e32
      maxts=-1.0e32
      mindm=1.0e32
      maxdm=-1.0e32
      maxsn=0.0
c==============================================================================
      open(lun,file=stem(1:lst)//'.pls',status='unknown')
      read(lun,'(a)') line
      read(line(7:),*) tsamp
      do while(.true.)
         i=i+1
         read(lun,*,err=1,end=1) dm(i),width,idx,sn(i),fft
         time(i)=idx*tsamp/1.0e6
         mints=min(time(i),mints)
         maxts=max(time(i),maxts)
         mindm=min(dm(i),mindm)
         maxdm=max(dm(i),maxdm)
         maxsn=max(sn(i),maxsn)
      enddo
 1    close(lun)
      
      npulses=i-1
      write(*,*) npulses,' events, DM range:',mindm,' ->',maxdm,
     &           ' pc/cc .... max S/N',maxsn
c==============================================================================
      call pgbegin(0,'?',1,1)
      call pgscf(2)
      call pgswin(mints,maxts,mindm,maxdm)
      call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
      call pglab('Time (s)','DM (cm\\u-3\\d pc)',
     &           'Single pulse plane for '//stem)
      do i=1,npulses
         sym=20
         if (sn(i).gt.5) sym=21
         if (sn(i).gt.7) sym=22
         if (sn(i).gt.9) sym=23
         if (sn(i).gt.10) sym=17
         if (sn(i).gt.minsn) call pgpt1(time(i),dm(i),sym)
      enddo
      call pgend
      end
c==============================================================================
