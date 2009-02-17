c==============================================================================
      program best
c==============================================================================
c
c     Looks through the .prd files to find the best candidates
c     from the analysis... (dunc@mpifr-bonn.mpg.de - November 1997)
c
c     Modification history:
c
c     98/02/02 -> dunc@mpifr-bonn.mpg.de - Acceleration options added
c     98/06/12 -> dlorimer@naic.edu      - code overhauled and commented
c                                        - harmonic zapping simplified
c                                        - -p option implemented
c     98/07/01 -> dunc@mpifr-bonn.mpg.de - -z option implemented
c                                        - directory etc listed on plot      
c     98/07/11 -> dunc@mpifr-bonn.mpg.de - -L and -H options implemented
c     98/07/15 -> dunc@mpifr-bonn.mpg.de - -a option no longer necessary
c                                        - adot read in and 2-D plot
c     98/07/17 -> dunc@mpifr-bonn.mpg.de - ".fld" file with folding commands
c     98/09/26 -> dunc@mpifr-bonn.mpg.de - assumes DM search if only one list
c     98/11/24 -> dunc@naic.edu          - corrected DM/ACIDX output feature
c     99/04/21 -> dunc@naic.edu          - added -F and -C options to calculate
c                                        - whether a candidate gets DM smeared
c     02/10/22 -> drl@jb.man.ac.uk       - -z option done by default -i inverts
c                                        - harmonic checking only for <16*P
c                                        - -t option allows tracing candidates
c     07/01/31 -> drl                    - has David Barnes' mods for nfolds
c      
c==============================================================================
c
      implicit none

      include 'vers.inc'
      include 'folds.inc'
      
      character*80 filename,option,title*160,susfile,cline,prf,ps
      character*70 string(10),period,acceln,accdot,dmidx*4,dmval*8

      integer fold,nc,mc,nt,mt,lst,lun2,lun3,lsum,length,mg,ngx,ngy
      integer ngz,ndmi,ndmm
c      parameter(mt=5000,mc=mt*512,mg=256)
c this change: mt needs to be roughly N(DM) * N(AC) * NFOLDS
c              mc needs to be roughly mt * nrows per set, which
c                 eg. is limited to 5 in one of my scripts.  50 is safe-ish
      parameter(mt=200000,mc=mt*50,mg=256)

      real*8 par(mc),ratio, ratio2
      
      real snr(mc), dm(mt), ac(mt), ad(mt), x(mt), y(mt), ymin,
     &     dmdsp, dd(mt), snrmax, snrmin, prdmin, prdmax, ymax,
     &     parmin(3), parmax(3), grid(mg,mg), tr(6), xmin, xmax,
     &     xval(mg), yval(mg), gridmin, gridmax, deltax, deltay,
     &     xx, yy, fcentre, chbandw,zval(mg)

      integer trid(mc),cidx(mc),i,j,k,l,sidx(mc),nids(mc),ndm
      integer narg,iargc,fld(mc),nsus,nf,mode,ix,iy,gidx,iz,bb

      logical lview, prd, frq, plotit, autop, harm, test, hzap,
     &        ok, lac, ldm, lad, consider, laccn, smeared, trace,
     &     lmode5
c
c==============================================================================
c
c     un-comment these lines if floating point exceptions
c     occur (heaven forbid!)
c
c      external myhandler
c      integer ieee_handler
c
c  TRAP DIVISION BY ZERO ERRORS and NaN
c      i=ieee_handler("set","common",myhandler)
c==============================================================================
c
c     If no command line input given - tell user what to do...
c
         gridmin=+1.0e32
         gridmax=-1.0e32
         do k=1,mg
            do l=1,mg
               grid(k,l)=0.0
            enddo
         enddo
      lmode5 = .false.
      narg=iargc()
      call getarg(1,option)
      if (option.eq.'version'.or.option.eq.'-version') then
         write(*,'(a,a)') 'PROGRAM: best ',version
         stop
      endif
      if (narg.lt.1.or.option.eq.'-help'.or.option.eq.'help') call help
c
c     Get filename and initialise some flags and variables
c
      call getarg(1,filename)
      lst=length(filename)-4
      cline='best '//filename(1:length(filename))
      l=length(cline)
      prd=index(filename,'.prd').gt.0.or.index(filename,'.top').gt.0
      frq=.not.prd
      fold=0
      autop=.false.
      snrmin=8.0
      prdmin=0.0
      prdmax=1.0e32
      lview=.false.
      hzap=.true.
      trace=.false.
      do i=1,2
         parmin(i)=-1.0e32
         parmax(i)=+1.0e32
      enddo
      chbandw=0.0
      fcentre=1.0
c
c     Sort out command-line arguments
c
      if (narg.gt.1) then
        do i=2,narg
          call getarg(i,option)
          cline=cline(1:l)//' '//option(1:length(option))
          l=length(cline)
          if (index(option,'-f').gt.0) read(option(3:),*) fold
          if (index(option,'-s').gt.0) read(option(3:),*) snrmin
          if (index(option,'-l').gt.0) read(option(3:),*) prdmin
          if (index(option,'-h').gt.0) read(option(3:),*) prdmax 
          if (index(option,'-L').gt.0) read(option(3:),*) parmin(1)
          if (index(option,'-H').gt.0) read(option(3:),*) parmax(1) 
          if (index(option,'-F').gt.0) read(option(3:),*) fcentre
          if (index(option,'-C').gt.0) read(option(3:),*) chbandw
          if (index(option,'-v').gt.0) lview=.true.
          if (index(option,'-p').gt.0) autop=.true.
          if (index(option,'-z').gt.0) hzap=.true.
          if (index(option,'-i').gt.0) hzap=.false.
          if (index(option,'-t').gt.0) trace=.true.
        enddo
      endif
c
c     Write out the command line inputs to a filename.bst file
c     so that the action can be repeated at a later date......
c      
      call glun(lun2)
      open(unit=lun2,file=filename(1:lst)//'.bst',status='unknown')
      write(lun2,'(a)') cline(1:length(cline))
      close(unit=lun2)
c
c     Read in the candidates from the ".prd" or ".frq" file...
c
      write(*,'('' File: '',a)') filename(1:index(filename,' ')-1)
      nc=0
      nt=0
      if (fold.eq.0) then
c         write(*,*) 'Folds: 1-5'
c         do fold=1,5
         write(*,*) 'Folds: 1-', nfolds
         do fold=1,nfolds
            call readsus(filename, fold, mc, nc, mt, nt,
     &                   par, snr, trid, fld, dm, ac, ad)
         enddo
	 nf=nfolds
      else
            write(*,*) 'Fold:',fold
            call readsus(filename, fold, mc, nc, mt, nt,
     &                   par, snr, trid, fld, dm, ac, ad)
	 nf=1
      endif
c
c     Look at dm,ac+ad arrays to establish parameter space of search
c
      ldm=.false.
      lac=.false.
      lad=.false.
      mode=0
      if (nt.eq.nfolds.and.nf.eq.nfolds) then
         ldm=.true. ! assume DM search if only one list
         mode=1
      else
         do i=1,nt
            if (i.gt.1.and.dm(i).ne.dm(i-1)) ldm=.true.
            if (ac(i).ne.0.0) lac=.true.
            if (ad(i).ne.0.0) lad=.true.
         enddo
      endif

      if (ldm.and.(.not.lac).and.(.not.lad)) then
         write(*,*) '1-D DM search...'
         mode=1
      else if ((.not.ldm).and.(lac.or.lad)) then
         if (lac.and.(.not.lad)) then
            write(*,*) '1-D AC search...'
            mode=2
         else if (lad.and.(.not.lac)) then
            write(*,*) '1-D AD search...'
            mode=3
         else if (lad.and.lac) then
            write(*,*) '2-D AC+AD search...'
            mode=4
         endif
      else if(ldm.and.lac.and.(.not.lad)) then
         write(*,*) '2-D DM+AC search...'
         mode=5
      endif
      
      if (mode.eq.0) stop 'Unrecognizable parameter space! HELP!'
c
c     Go through list, leaving out candidates that do not fall in
c     selected ranges.
c
      j=0
      k=0
      if (mode.eq.5) then
         call gdim(dm,nt,ngx,parmin(1),parmax(1),xval)
         call gdim(ac,nt,ngy,parmin(2),parmax(2),yval)
      endif
      do i=1,nc
         consider=(par(i).ge.prdmin.and.par(i).le.prdmax)
         k=trid(i)
         if (mode.eq.1.and.(dm(k).lt.parmin(1).or.dm(k).gt.parmax(1)))
     &   consider=.false.         
         if (mode.eq.1.and.smeared(real(par(i)),dm(k),fcentre,chbandw))
     &   consider=.false.         
         if (mode.eq.2.and.(ac(k).lt.parmin(1).or.ac(k).gt.parmax(1)))
     &   consider=.false.         
         if (mode.eq.3.and.(ad(k).lt.parmin(1).or.ad(k).gt.parmax(1)))
     &   consider=.false.
         if (mode.eq.5) then
            ix=gidx(xval,ngx,dm(k))
            iy=gidx(yval,ngy,ac(k))
            grid(ix,iy)=max(grid(ix,iy),snr(i))
            gridmin=min(gridmin,grid(ix,iy))
            gridmax=max(gridmax,grid(ix,iy))
         endif
         if (consider) then
            j=j+1
            par(j)=par(i)
            snr(j)=snr(i)
            trid(j)=trid(i)
            fld(j)=fld(i)
         endif
      enddo
      nc=j
c      
c     Initialise arrays ready for cross-correlation analysis
c
      do i=1,nc
         cidx(i)=0
         sidx(i)=0
         nids(i)=0
      enddo
c
c     Establish ranges of DMs (or ACs) searched
c
      do i=1,2
         parmin(i)=1.0e32
         parmax(i)=0.0
      enddo
      dmdsp=dm(1)
      do i=1,nt
         if (mode.eq.1) then
            dd(i)=dm(i)
         else if (mode.eq.2) then
            dd(i)=ac(i)
         else if (mode.eq.3) then
            dd(i)=ad(i)
         endif
c         if (mode.eq.5) then
c            parmin(1)=min(parmin(1),dm(i))
c            parmax(1)=max(parmax(1),dm(i))
c            parmin(2)=min(parmin(2),ac(i))
c            parmax(2)=max(parmax(2),ac(i))
c         else
            parmin(1)=min(parmin(1),dd(i))
            parmax(1)=max(parmax(1),dd(i))
c         endif
      enddo

c
c     Display info concerning number of candidates read in
c
      if (mode.eq.1) then
        write(*,*) nt/nf,' DM group(s). ',nc,' candidates'
        write(*,*) 'DM range:',parmin(1),parmax(1),' pc/cc'
      else if (mode.eq.2) then
        write(*,*) nt/nf,' AC group(s). ',nc,' candidates'
        write(*,*) 'AC range:',parmin(1),parmax(1),' m/s/s'
      else if (mode.eq.3) then
        write(*,*) nt/nf,' AD group(s). ',nc,' candidates'
        write(*,*) 'AD range:',parmin(1),parmax(1),' cm/s/s/s'
      else if (mode.eq.4) then
        write(*,*) nt/nf,' AD/AC group(s). ',nc,' candidates'
        call gdim(ac,nt,ngx,parmin(1),parmax(1),xval)
        call gdim(ad,nt,ngy,parmin(2),parmax(2),yval)
        write(*,*) 'AC x AD grid dimensions:',ngx,' x',ngy
        write(*,*) 'AC range:',parmin(1),parmax(1),' m/s/s'
        write(*,*) 'AD range:',parmin(2),parmax(2),' cm/s/s/s'
      else if (mode.eq.5) then
         write(*,*) nt/nf,' DM/AC group(s). ',nc,' candidates'
         call gdim(dm,nt,ngx,parmin(1),parmax(1),xval)
         call gdim(ac,nt,ngy,parmin(2),parmax(2),yval)
         write(*,*) 'DM range:',parmin(1),parmax(1),' pc/cc'
         write(*,*) 'AC range:',parmin(2),parmax(2),' m/s/s'
         write(*,*) 'DM x AC grid dimensions:',ngx,' x',ngy
      endif

      
c     
c     Sort signals in terms of their signal-to-noise ratios
c
      call indexxf77(nc,snr,sidx)
c
c     Cross-correlate signals for integer and non-integer harmonics
c     when found, these are then flagged out using the array cidx() 
c
      if (hzap) then
         write(*,*)'Zapping integer+non-integer harmonics...'
      else
         write(*,*)'Zapping only the integer harmonics...'
      endif
c
c     This loop is not as horrendous as it looks... (well ok, it is!)
c      
      do i=nc,1,-1
         k=sidx(i)  ! index of outer loop in terms of descending S/N
         if (cidx(k).eq.0) then ! has this signal already been seen?
            do j=1,nc
               test=par(j).ne.0.0.and.par(k).ne.0.0.and.cidx(j).eq.0
               if (test) then 
                  ratio=par(k)/par(j)                    
                  if (ratio.lt.1.0) ratio=1.0/ratio      
                  ratio2=ratio
                  ratio=ratio-int(ratio)
                  harm=(ratio.gt.0.98.or.ratio.lt.0.02)
                  if (hzap) then ! non-integer harmonics clobbered here
                   harm=harm.or.(ratio.gt.0.16.and.ratio.lt.0.167)
                   harm=harm.or.(ratio.gt.0.19999.and.ratio.lt.0.20001)
                   harm=harm.or.(ratio.gt.0.24995.and.ratio.lt.0.25005)
                   harm=harm.or.(ratio.gt.0.29999.and.ratio.lt.0.30001)
                   harm=harm.or.(ratio.gt.0.3332.and.ratio.lt.0.3334)
                   harm=harm.or.(ratio.gt.0.39.and.ratio.lt.0.41)
                   harm=harm.or.(ratio.gt.0.19.and.ratio.lt.0.21)
                   harm=harm.or.(ratio.gt.0.76.and.ratio.lt.0.78)
                   harm=harm.or.(ratio.gt.0.4995.and.ratio.lt.0.5005)
                   harm=harm.or.(ratio.gt.0.5997.and.ratio.lt.0.6003)
                   harm=harm.or.(ratio.gt.0.6665.and.ratio.lt.0.6668)
                   harm=harm.or.(ratio.gt.0.7499.and.ratio.lt.0.7501)
                   harm=harm.or.(ratio.gt.0.7997.and.ratio.lt.0.8003)
                     harm=harm.or.(ratio.gt.0.249.and.ratio.lt.0.251)
                     harm=harm.or.(ratio.gt.0.332.and.ratio.lt.0.334)
                     harm=harm.or.(ratio.gt.0.499.and.ratio.lt.0.501)
                     harm=harm.or.(ratio.gt.0.599.and.ratio.lt.0.601)
                     harm=harm.or.(ratio.gt.0.665.and.ratio.lt.0.667)
                     harm=harm.or.(ratio.gt.0.749.and.ratio.lt.0.751)
                  endif
c
c     check to see if candidate is harmonically related out to 32nd
c
c                  write(*,*) par(k),par(j),harm,ratio,ratio2
                  if (harm.and.ratio2.lt.16.0) then
                     cidx(j)=k
                     if (snr(k).ge.snrmin.and.snr(j).ge.snrmin)
     &                    nids(k)=nids(k)+1
                     if (j.ne.k) then
                       cidx(j)=-1*cidx(j)
                       if (trace) 
     &                 write(*,*) par(j),' related to',par(k),ratio2
                    endif
                  endif
               endif
            enddo
         endif
      enddo
c
c     Column headers for screen output...
c
      write(string(1),'(a,a)')'   P (ms)      S/N    AC   ACID NIDs f',
     &                     '    P/Ptop    Ptop/P'
      write(string(2),'(a,a)')'   f (Hz)      S/N    AC   ACID NIDs f',
     &                     '    f/ftop    ftop/f'
      write(string(3),'(a,a)')'   P (ms)      S/N    DM   DMID NIDs f',
     &                     '    P/Ptop    Ptop/P'
      write(string(4),'(a,a)')'   f (Hz)      S/N    DN   DMID NIDs f',
     &                     '    f/ftop    ftop/f'
      write(string(5),'(a,a)')'   P (ms)      S/N    AD   ACID NIDs f',
     &                     '    P/Ptop    Ptop/P'
      write(string(6),'(a,a)')'   f (Hz)      S/N    AD   ACID NIDs f',
     &                     '    f/Ptop    ftop/f'
      write(string(7),'(a,a)')'   P (ms)      S/N    AC     AD   ACID',
     &                     ' NIDs f    P/Ptop    Ptop/P'
      write(string(8),'(a,a)')'   f (HZ)      S/N    AC     AD   ACID',
     &                     ' NIDs f    f/ftop    ftop/f'
      write(string(9),'(a,a)')'   P (ms)      S/N    DM      AC',
     &      '      AD  ACID NIDs f    P/Ptop    Ptop/P'
      write(string(10),'(a,a)')'   f (HZ)      S/N    DM    AC     AD',
     &                     ' ACID    NIDs f    f/ftop    ftop/f'

      if (mode.eq.1) then
        if (prd) write(*,*) string(3)
        if (frq) write(*,*) string(4)
        laccn=.false.
      else if (mode.eq.2) then
        if (prd) write(*,*) string(1)
        if (frq) write(*,*) string(2)
        laccn=.true.
      else if (mode.eq.3) then
        if (prd) write(*,*) string(5)
        if (frq) write(*,*) string(6)
        laccn=.true.
      else if (mode.eq.4) then
        if (prd) write(*,*) string(7)
        if (frq) write(*,*) string(8)
        laccn=.true.
      else if (mode.eq.5) then
        if (prd) write(*,*) string(9)
        if (frq) write(*,*) string(10)
        laccn=.true.
      endif
c
c     Print out results of cross correlation analysis.... ie list
c     only the signals that have not been flagged in the cidx()
c     array as being related to something else (i.e. cidx(i)=0.0)
c     this output goes out to a "*.lis" file as well as to the sceen.
c
      call glun(lun2)
      open(unit=lun2,file=filename(1:index(filename,' ')-5)//'.lis',
     &     status='unknown')
      call glun(lun3)
      open(unit=lun3,file=filename(1:index(filename,' ')-5)//'.fld',
     &     status='unknown')
      nsus=0
      do i=nc,1,-1
         j=sidx(i)
         if (cidx(j).gt.0.and.snr(j).ge.snrmin.and.
     &       par(j).le.prdmax.and.par(j).ge.prdmin) then
            ratio=par(j)/par(sidx(nc))
            if (mode.lt.4) then
            nsus=nsus+1
            write(*,1) par(j),snr(j),dd(trid(j)),trid(j),nids(j),
     &                 fld(j),ratio,1.0/ratio
            write(lun2,1) par(j),snr(j),dd(trid(j)),trid(j),nids(j),
     &                 fld(j),ratio,1.0/ratio,
     &                 filename(1:index(filename,' ')-5),nsus
            call wstr8(par(j),period)
            call wstr4(dd(trid(j)),acceln)
            write(dmidx,'(i4.4)') trid(j)
            write(dmval,'(f8.2)') dm(trid(j))
            write(prf,'(a,i3,a)') '_sus',nsus,'.prf'
            if (nsus.lt.100) write(prf,'(a,i2,a)') 'sus',nsus,'.prf'
            if (nsus.lt.10) write(prf,'(a,i1,a)') 'sus',nsus,'.prf'
            write(ps,'(a,i3,a)') '_sus',nsus,'.ps/ps'
            if (nsus.lt.100) write(ps,'(a,i2,a)') '_sus',nsus,'.ps/ps'
            if (nsus.lt.10) write(ps,'(a,i1,a)') '_sus',nsus,'.ps/ps'
            if (mode.eq.1) then
               write(lun3,'(a)')    'filterbank '//
     &              filename(1:index(filename,' ')-5)//'.dat'//
     &              '| dedisperse -d '//dmval//'| fold -p '//
     &              period(1:length(period))//' -n 128 > '//
     &              filename(1:index(filename,' ')-5)//prf//
     &              '; plot prof '//         
     &              filename(1:index(filename,' ')-5)//prf//
     &              '-nobox -centre -T"P: '//period(1:length(period))//
     &              ' ms DM: '//dmval//' cm-3 pc"'//
     &              ' -d'//filename(1:index(filename,' ')-5)//ps
            else
               write(lun3,'(a)') 'fold '//
     &              filename(1:index(filename,' ')-5)//'.tim'//
     &                        ' -p '//period(1:length(period))//
     &                        ' -a '//acceln(1:length(acceln))//
     &           ' > '//prf
            endif
            else if(mode.eq.4) then
            write(*,2) par(j),snr(j),ac(trid(j)),ad(trid(j)),trid(j),
     &                 nids(j),fld(j),ratio,1.0/ratio,
     &                 filename(1:index(filename,' ')-5)
            write(lun2,2) par(j),snr(j),ac(trid(j)),ad(trid(j)),trid(j),
     &                 nids(j),fld(j),ratio,1.0/ratio,
     &                 filename(1:index(filename,' ')-5)
            call wstr8(par(j),period)
            call wstr4(ac(trid(j)),acceln)
            call wstr4(ad(trid(j)),accdot)
            write(lun3,'(a)') 'fold '//
     &              filename(1:index(filename,' ')-5)//'.ser'//
     &                        ' -p'//period(1:length(period))//
     &                        ' -a'//acceln(1:length(acceln))//
     &                        ' -d'//accdot(1:length(accdot))//
     &                        ' -s8; grey'
         else if(mode.eq.5) then
            write(*,12) par(j),snr(j),dm(trid(j)),ac(trid(j)),
     &           ad(trid(j)),trid(j),
     &           nids(j),fld(j),ratio,1.0/ratio,
     &           filename(1:index(filename,' ')-5)
            write(lun2,12) par(j),snr(j),dm(trid(j)),ac(trid(j)),
     &           ad(trid(j)),trid(j),
     &           nids(j),fld(j),ratio,1.0/ratio,
     &           filename(1:index(filename,' ')-5)
            call wstr8(par(j),period)
            call wstr4(ac(trid(j)),acceln)
            call wstr4(ad(trid(j)),accdot)
            write(lun3,'(a)') 'fold '//
     &           filename(1:index(filename,' ')-5)//'.ser'//
     &           ' -p'//period(1:length(period))//
     &           ' -a'//acceln(1:length(acceln))//
     &           ' -d'//accdot(1:length(accdot))//
     &           ' -s8; grey'
         endif
      endif
      enddo
      close(unit=lun3)
      close(unit=lun2)

      if (mode.eq.5) then
         deltax=(parmax(1)-parmin(1))/real(ngx)
         deltay=(parmax(2)-parmin(2))/real(ngy)
         tr(1)=parmin(1)-0.5*deltax
         tr(2)=deltax
         tr(3)=0.0
         tr(4)=parmin(2)-0.5*deltay
         tr(5)=0.0
         tr(6)=deltay
         call pgbegin(0,'/xs',1,1)
         call pgvport(0.1,0.9,0.1,0.9)
         call pgwindow(parmin(1),parmax(1),parmin(2),parmax(2))
         call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
         call pgsch(1.0)
         call pgscf(2)
         call pglabel('DM (cm\\u-3\\d pc)',
     &        'Acceleration (m s\\u-2\\d)',' ')
         call pggray(grid,mg,mg,1,ngx,1,ngy,gridmax,gridmin,tr)
         call pgend
         stop
      endif


c
c     Stop here if nothing else to be done (default mode)
c
      if ((.not.autop).and.(.not.lview)) stop
c
c     Carry on either in interactive (view) mode where the
c     postscript files are saved on demand, or in automatic
c     mode (-p) which writes all the postscript output
c
      if (lview) write(*,*) 'Cycling through the list...'
      nsus=0
      do i=nc,1,-1
         j=sidx(i)
         if (cidx(j).gt.0.and.snr(j).ge.snrmin.and.
     &       par(j).le.prdmax.and.par(j).ge.prdmin) then
            if (mode.eq.2.and.lview.and.(.not.autop)) then
               if (prd) write(*,*) string(1)(1:24)
               if (frq) write(*,*) string(2)(1:24)
            else if (mode.eq.3.and.lview.and.(.not.autop)) then
               if (prd) write(*,*) string(5)(1:24)
               if (frq) write(*,*) string(6)(1:24)
            else if (mode.eq.4.and.lview.and.(.not.autop)) then
               if (prd) write(*,*) string(7)(1:32)
               if (frq) write(*,*) string(8)(1:32)
            else if (mode.eq.5.and.lview.and.(.not.autop)) then
               if (prd) write(*,*) string(9)(1:40)
               if (frq) write(*,*) string(10)(1:40)
            else if (.not.autop) then
               if (prd) write(*,*) string(3)
               if (frq) write(*,*) string(4)
c               stop 'HELP!'
            endif
            if (lview) then
               if (mode.lt.4) then
                  write(*,5) par(j),snr(j),dd(trid(j))
               else if(mode.eq.4) then
                  write(*,6) par(j),snr(j),ac(trid(j)),
     &                 ad(trid(j))
               else if(mode.eq.5) then
                  write(*,13) par(j),snr(j),dm(trid(j)),ac(trid(j)),
     &                 ad(trid(j))
               endif
            endif
            nsus=nsus+1
c
c           Clear the grid if in 2-D mode
c      
            if (mode.eq.4.or.mode.eq.5) then
               gridmin=+1.0e32
               gridmax=-1.0e32
               do k=1,mg
                  do l=1,mg
                     grid(k,l)=0.0
                  enddo
               enddo
            endif
c
c           Open a summary file (list of signal-to-noise versus DM/AC
c           for reading by other programs... only done in auto mode
c
            call glun(lsum)
            if (nsus.lt.10) then
               write(susfile,'(a,i1,a)')'sus00',nsus,'.sum'
            else if (nsus.lt.100) then
               write(susfile,'(a,i2,a)')'sus0',nsus,'.sum'
            else
               write(susfile,'(a,i3,a)')'sus',nsus,'.sum'
            endif
            open(unit=lsum,file=susfile,status='unknown')
c	MK 2006: Added the name of the dedispersed
c		file. This is for 'tune' used with MMB  
		
            write(dmval,'(f8.2)') dd(trid(j))

		bb=index(dmval,' ')
               do while (bb.eq.1)
                  dmval=dmval(2:8)
                  bb=index(dmval,' ')
               enddo			
		write(lsum,*)filename(1:index(filename,' ')-5)//
     &              '_dice_dm'//dmval(1:bb-1)//'.tim'

            if (mode.eq.4) then
               write(lsum,*) par(j),snr(j),ac(trid(j)),ad(trid(j))
               xx=ac(trid(j))
               yy=ad(trid(j))
            else if(mode.eq.5) then
               write(lsum,*) par(j),snr(j),dm(trid(j)),ac(trid(j)),
     &              ad(trid(j))
               xx=ac(trid(j))
               yy=ad(trid(j))
            else
               write(lsum,*) par(j),snr(j),dd(trid(j)),fld(j)
            endif
            

c     We need to reset the x and y arrays for mode 5
            if(mode.eq.5) then
               do ndmi = 0,5000
                  x(ndmi) = 0
                  y(ndmi) = 0
               enddo
            endif


            ndm=0
            ndmi = 0
            ndmm = 1
            snrmax=0.0

            do k=1,nc
               if (abs(cidx(k)).eq.j.and.fld(k).eq.fld(j))then
                  ok=.true.
                  if (mode.eq.5) then
                     l=trid(k)
                     ix=gidx(xval,ngx,dm(l))
                     iy=gidx(yval,ngy,ac(l))
                     grid(ix,iy)=grid(ix,iy)+snr(k)
                     write(lsum,*) ix,iy,tr(1)+real(ix)*tr(2),
     &               tr(4)+real(iy)*tr(6),grid(ix,iy)
                     gridmin=min(gridmin,grid(ix,iy))
                     gridmax=max(gridmax,grid(ix,iy))
                     deltax=(parmax(1)-parmin(1))/real(ngx)
                     deltay=(parmax(2)-parmin(2))/real(ngy)
                     tr(1)=parmin(1)-0.5*deltax
                     tr(2)=deltax
                     tr(3)=0.0
                     tr(4)=parmin(2)-0.5*deltay
                     tr(5)=0.0
                     tr(6)=deltay
                  else
                     if (ndm.gt.0) ok=dd(trid(k)).ne.x(ndm)
                     if (ok) then   
                        ndm=ndm+1
                        x(ndm)=dd(trid(k))
                        y(ndm)=snr(k)
                        write(lsum,*) ndm,x(ndm),y(ndm)
                        snrmax=max(snrmax,snr(k))
                     endif
                  endif
               endif
            enddo
            close(lsum)

c     test

          

c
c           For the main plot, write out a title
c
            write(title,3) par(j),snr(j),dd(trid(j))
            if (laccn) write(title,4) par(j),snr(j),dd(trid(j)),dmdsp
            if (mode.eq.4) write(title,7) par(j),snr(j),xx,yy
            if (mode.eq.5) write(title,8) par(j),snr(j),
     &           dd(trid(k)),ac(trid(j))
c
c           Sort out whether this plot will get displayed...
c           always the case in auto mode, on demand otherwise
c    
c
 

	    plotit=.true.
            option=' '
            if (lview) then
               write(*,'('' View this one [y] ''$)')
               read(*,'(a)')option
               if (option(1:1).ne.'n'.and.option(1:1).ne.'N') then 
                  call pgbegin(0,'/xs',1,1)
               else
                  plotit=.false.
               endif
            else if (autop) then
               call pgbegin(0,
     &         susfile(1:index(susfile,'.sum')-1)//'.ps/ps',1,1)
               plotit=.true.
            else
               plotit=.false.
            endif
c
c           Plot it
c

            lmode5 = .false.

	    if (plotit) then
 200           if(lmode5) then
                  mode = 1
                  call pgvport(0.12,0.85,0.50,0.85)
                  call pgsch(1.0)
                  call pgscf(2)
                  
               else if(mode.eq.5) then
                  call pgvport(0.12,0.40,0.12,0.40)
                  call pgsch(1.0)
                  call pgscf(2)
                  mode = 4
                  lmode5 = .true.
                
               else
                  call pgvport(0.15,0.85,0.15,0.85)
                  call pgsch(1.3)
                  call pgscf(2)
               endif
               if (mode.eq.4) then
                  xmin=tr(1)+tr(2)*0.5
                  xmax=tr(1)+tr(2)*(real(ngx)+0.5)
                  ymin=tr(4)+tr(6)*0.5
                  ymax=tr(4)+tr(6)*(real(ngy)+0.5)
               else
c                  xmin=parmin(1)
                  xmin=dd(1)
                  xmax=parmax(1)
                  ymin=0.0
                  ymax=snrmax*1.1
               endif
               call pgwindow(xmin,xmax,ymin,ymax)
               call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
               if (mode.eq.1) then
                  call pglabel('Trial DM (cm\\u-3\\d pc)',
     &                'Signal-to-Noise Ratio',' ')
               else if (mode.eq.2) then
                  call pglabel('Trial Acceleration (m s\\u-2\\d)',
     &                'Signal-to-Noise Ratio',' ')
               else if (mode.eq.3) then
                  call pglabel('Trial Adot (cm s\\u-3\\d)',
     &                'Signal-to-Noise Ratio',' ')
               else if (mode.eq.4) then
                  call pglabel('Trial Acceleration (m s\\u-2\\d)',
     &                'Trial Adot (cm s\\u-3\\d)',' ')
               else if (mode.eq.5) then
                  call pglabel('Trial DM (cm\\u-3\\d pc)',
     &                'Trial Acceleration (m s\\u-2\\d)',' ')
               endif
               if (mode.eq.2.or.mode.eq.3) then
                  call pgsls(2)
                  call pgmove(parmin(1),snrmin)
                  call pgdraw(parmax(1),snrmin)
               endif
               if (mode.lt.4) then
                  call pgsls(4)
                  call pgmove(0.0,0.0)
                  call pgdraw(0.0,snrmax*1.1)
                  call pgsls(1)
                  call pgline(ndm,x,y)
                  call pgpoint(ndm,x,y,17)
               else if (mode.eq.4.or.mode.eq.5) then
                  call pggray(grid,mg,mg,1,ngx,1,ngy,gridmax,gridmin,tr)
                  call pgpoint(1,xx,yy,27)
                  call pgsls(4)
                  call pgmove(xmin,0.0)
                  call pgdraw(xmax,0.0)
                  call pgmove(0.0,ymin)
                  call pgdraw(0.0,ymax)
                  call pgsls(1)
               endif
              
              

               if(lmode5.and.mode.eq.4) then
                  goto 200
               else if(lmode5.and.mode.eq.1) then
                  mode = 5
                  lmode5 = .false.
                  
               endif
               call pgsch(1.0)
               call pgwindow(0.0,1.0,0.0,1.0)
               call pgtext(-0.05,1.25,title)
               call getenv('PWD',title)
               title='Directory: '//title(1:index(title,' ')-1)
     &              //' File: '//filename(1:lst)
               call pgtext(-0.05,1.35,title)
               
               call pgiden
               call pgend
            endif
         endif
      enddo
c
c     Formats for various screen/plot outputs
c      
 1    format(f13.8,1x,f6.1,1x,f7.1,1x,2(i4.4,1x),
     &   i1,2(1x,f9.4),1x,a,1x,i4)
 2    format(f13.8,1x,f6.1,1x,2(f7.1,1x),2(i4.4,1x),i1,2(1x,f9.4),1x,a)
 3    format('P =',f9.4,' ms S/N = ',f6.1,' DM = ',f7.2,
     &' cm\\u-3\\d pc')
 4    format('P =',f9.4,' ms S/N = ',f6.1,' AC = ',f7.1,
     &' m s\\u-2\\d (DM = ',f7.2,')')
 5    format(f13.8,1x,f6.1,1x,f6.1,1x)
 6    format(f13.8,1x,f6.1,1x,2(f6.1,1x))
 7    format('P =',f9.4,' ms S/N = ',f6.1,' AC = ',f7.1,
     &' m s\\u-2\\d AD = ',f7.2,' cm s\\u-3\\d')
 8    format('P =',f9.4,' ms S/N = ',f6.1,' AC = ',f7.1,
     &' m s\\u-2\\d AD = ',f7.2,' cm s\\u-3\\d DM = ',f7.2,
     &' cm\\u-3\\d pc')

 12   format(f13.8,1x,f6.1,f6.1,1x,2(f7.1,1x),2(i4.4,1x),i1,2(1x,f9.4),
     &     1x,a)
 13   format(f13.8,1x,f6.1,1x,3(f6.1,1x))

      end
c==============================================================================
      subroutine help
      write(*,*)
      write(*,1) 'best - displays the best suspects from SEEK'
      write(*,*)
      write(*,1) 'usage: best <INFILE> -{options}'
      write(*,*)
      write(*,1)
     & 'The input file may be of the ".prd", ".top" or ".frq" variety.'
      write(*,*)
      write(*,1)'options:'
      write(*,*)
      write(*,1)'-v: view output selectively'
      write(*,1)'-p: produce postscript file output automatically'
      write(*,1)'-i: zap integer harmonics only'
      write(*,1)'-t: trace mode... show all harmonic relationships'
      write(*,*)
      write(*,1)'-f[fold]: read a specific harmonic fold (1-5;def=all)'
      write(*,1)'-s[smin]: set signal-to-noise threshold smin (def=8)'
      write(*,1)'-l[mini]: set minimum prd/frq to consider (ms;def=0)'
      write(*,1)'-h[maxi]: set maximum prd/frq to consider (ms;def=inf)'
      write(*,1)'-L[mini]: set minimum ac/dm to consider'
      write(*,1)'-H[maxi]: set maximum ac/dm to consider'
      write(*,1)'-F[fMHz]: set centre freq for DMsmear test (optional)'
      write(*,1)'-C[fMHz]: set channel band for DMsmear test (optional)'
      write(*,*)
 1    format(a)
      stop
      end
c==============================================================================
c      integer function myhandler(sig,code,context)
c==============================================================================
c      integer sig, code, context(5)
c      write(*,*)'ieee exception'
c      call abort()
c      end
c==============================================================================
      subroutine gdim(x,nt,ng,xmin,xmax,xval)
c==============================================================================
c
c     Gets the dimensions of a 1-D slice of parameter space
c
c     x(nt) - r4 - array containing signal-to-noise ratios --- passed
c     nt    - i4 - integer nunber of elements in the slice --- passed
c     ng    - i4 - number of distinct grid elements        --- return
c     xmin  - r4 - lowest grid element                     --- return
c     xmax  - r4 - highest grid element                    --- return
c     xval  - r4 - array containing valiues of ng elements --- return
c
c     Created: 98/07/15 (dunc@mpifr-bonn.mpg.de)
c      
c==============================================================================
c
      implicit none
      integer nt,ng
      real x(nt),xmin,xmax,xval(*)
      
      integer i,j
      logical newpnt

      ng=0
      xmax=-1.0e32
      xmin=+1.0e32

      do i=1,nt
         xmin=min(xmin,x(i))
         xmax=max(xmax,x(i))
         newpnt=.true.
         if (ng.gt.0) then
            do j=1,ng
               if (x(i).eq.xval(j)) newpnt=.false.
            enddo
         endif
         if (newpnt) then
            ng=ng+1
            xval(ng)=x(i)
         endif
      enddo

      end
c==============================================================================
      integer function gidx(xval,nx,x)
c==============================================================================
c
c     Returns the index "gidx" in the array xval(nx) of the value 'x'
c
c     xval(*) - r4 - array containing values of ng elements  --- passed
c     nx      - i4 - number of used elements in xval         --- passed
c     x       - r4 - value whose index in xval is required   --- passed
c
c     Created: 98/07/15 (dunc@mpifr-bonn.mpg.de)
c      
c==============================================================================
c
      implicit none
      integer nx
      real xval(*),x
      integer i
      gidx=0
      do i=1,nx
         if (xval(i).eq.x) then
            gidx=i
            return
         endif
      enddo
      stop 'HELP! Grid value not found! - see GIDX for more info...'
      end
c==============================================================================
      subroutine wstr4(param4,cstring)
      implicit none
      real*4 param4
      real*8 param8
      character*(*) cstring
      integer i,length,j
      write(cstring,*) param4
      goto 1
      entry wstr8(param8,cstring)
      write(cstring,*) param8
 1    j=0
      do i=1,length(cstring)
         if (cstring(i:i).ne.' '.and.j.eq.0) j=i
      enddo
      cstring=cstring(j:)
      end
c==============================================================================
      logical function smeared(period,dm,f0,chbw)
c
c     Returns true if smearing across a single filterbank channel
c     of width chbw (MHz) for a centre frequency f0 (MHz) is larger
c     than the period of the signal... ie if it's dispersion smeared
c
      implicit none
      real period,dm,f0,chbw
      smeared=(8.3e6*dm*chbw/f0/f0/f0).gt.period
      end
c==============================================================================
