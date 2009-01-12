ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program FFAmain
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     new version:
c     reads in tim files - mk 20/02/04
c     merged into sigproc-3.7 and produces a .prd file - drl 03/04/06
c     added search over 5 different values of rebinnin - mk 04/04/06
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'vers.inc'

      integer ndataMax, ndata
      parameter (NdataMax = 2**24)
      real time(ndataMax)

      integer iargc,i, ia
      real p1, p2, log2, refdm

      integer isamp, ip1, ip2, nrebin, dmidx, ndatanew
      character*80 infile, dir*80, arg1*80, arg2*80

      logical rescale, zeropad, append, searchw, dumpspec

      ip1 = 1000  ! msec
      ip2 = 20000 ! msec
      nrebin = 1                  ! factor of rebinning

      dir='./'

      rescale = .false.
      zeropad = .false.
      append  = .false.
      searchw = .false.
      dumpspec = .false.
      
c     Get input parameters (note change of ordering -- drl apr 06)

      ia=iargc()
      call getarg(1,arg1)
      if(ia.lt.1 .or. arg1(1:2).eq.'-h') go to 999

      write(*,'(a,a)') ' PROGRAM: FFA ',version
      write(*,'(a)') ' written by MK & DRL based on code by P.Mueller'
      if(arg1.eq.'-version') stop
      call timstart(6)
      call getarg(1,infile)

      if (ia.gt.1) then
      
         i=2
         do while(i.le.ia)
            call getarg(i,arg1)  
            
            if(arg1(1:3).eq.'-p1') then
               i=i+1
               call getarg(i,arg2)  
               if(arg2(1:1).eq.'-')then
                  write(*,*)'ERROR: Argument required for option ',
     +                 arg1(1:2)
                  STOP
               else
                  read(arg2,*) p1
                  ip1=int(p1)
               endif   
               
            else if(arg1(1:3).eq.'-p2') then
               i=i+1
               call getarg(i,arg2)  
               if(arg2(1:1).eq.'-')then
                  write(*,*)'ERROR: Argument required for option ',
     +                 arg1(1:2)
                  STOP
               else
                  read(arg2,*) p2
                  ip2=int(p2)
               endif   
               
            else if(arg1(1:3).eq.'-rb') then
               i=i+1
               call getarg(i,arg2)  
               if(arg2(1:1).eq.'-')then
                  write(*,*)'ERROR: Argument required for option ',
     $                 arg1(1:2)
                  STOP
               else
                  read(arg2,*) nrebin
               endif
               
            else if(arg1(1:2).eq.'-z') then
               zeropad=.true.
               
            else if(arg1(1:2).eq.'-A') then
               append=.true.
               
            else if(arg1(1:2).eq.'-r') then
               rescale=.true.

            else if(arg1(1:2).eq.'-w') then
               searchw=.true.

            else if(arg1(1:2).eq.'-s') then
               dumpspec=.true.
            endif
            
            i=i+1
         enddo
      endif
      
      if (searchw.and.(nrebin.gt.1)) then
         write(*,*) 'Cannot choose rebinning and search for width',
     +        ' simultaneously...'
         stop
      endif

      call readsp(infile, ndataMax,ndata,time,isamp,refdm)

      log2=log(real(ndata))/log(2.)
      write(*,*) '#samples = 2**',log2,' = ',ndata
      log2=int(log2)
      if (zeropad) then
         ndatanew = 2**(log2+1)
         if (ndatanew.lt.ndatamax) then
            write(*,*) 'ZEROPADDING!'
            do i=ndata+1, ndatanew
               time(i)=0.0
            enddo
            ndata=ndatanew
         else
            ndata=2**log2
         endif
      else
         ndata=2**log2
      endif
      write(*,*) '#samples using = ',ndata
         
c     .. dmdix not used...
      call runffa(isamp,ndata,time,ip1,ip2,nrebin,dmidx,infile,dir,
     +            rescale,refdm,append,dumpspec,searchw)

      call timfinis(6)
      stop
 999  continue

      write(*,*)
      write(*,'(a)') 
     +'ffa - perform a fast folding analysis on time series data'
      write(*,*)
      write(*,'(a)') 'usage: ffa {filename} -{options}'
      write(*,*)
      write(*,'(a)') 'options:'
      write(*,*)
      write(*,'(a,i5)')'-p1 val - specify lower period (ms), def= ',ip1
      write(*,'(a,i5)')'-p2 val - specify upper period (ms), def= ',ip2
      write(*,'(a,i5)')'-rb fac - rebin data by factor, def= ',nrebin
      write(*,'(a)')'-A      - append output to .prd file (def=newfile)'
      write(*,'(a)')'-z      - zero-pad data to nearest power of 2'
      write(*,'(a)')'-r      - rescale sampling time to fix PKS problem'
      write(*,'(a)')'-w      - search width (=tries several rebins)'
      write(*,*)
      write(*,'(a)')'NOTE:  -w and -rb are exclusive choices'
      write(*,*)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readsp(infile, ndataMax, ndata,time, isamp, refdm)
c
c     read time series in SIGPROC's ".tim" format
c     requires linking with the file "readtim.c"
c
      implicit none

      character infile*(*)
      integer ndata, isamp,ndataMax
      real time(ndataMax)
      real refdm,tsamp,sample
      integer  readsample, closefile

      ndata=0
      call readhd(infile,refdm,tsamp)
      do while(ndata.lt.ndataMax)
        if (readsample(sample).eq.0) goto 1
        ndata=ndata+1
        time(ndata)=sample
      enddo
      write(*,*) 'WARNING: Array size limit reached',ndatamax
 1    continue
      isamp=int(tsamp*1.0e6+0.01) ! fixing some rounding problem
      write(*,*) 'Read',ndata,' samples...'
      write(*,*) 'Reference DM:',refdm
      write(*,*) 'Sampling time:',isamp,' us'
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine runffa(isamp,ndata,timeIn,ip10,ip20,np,ind,fname,dir,
     +                  fixsamp,refdm,append,dumpspec,searchw)

c     using code written by Peter Mueller, MPIfR
c     modified by M. Kramer for PMSURV

      implicit none
      integer isamp, ndata, ip1, ip2, np, ind,ip10,ip20
      real timeIn(ndata),refdm
      character fname*(*), outname*80, dir*(*), name*20
      logical fixsamp,append,dumpspec,searchw, init

      integer ndataMax, iunit
      parameter (ndataMax=2**24)
      real time(ndataMax)

      integer nspec, nprof
      parameter (nspec=8192)
      parameter (nprof=8192)
      real pspec(nspec), signal(nprof), snspec(nspec)
      byte bspec(nspec)

      integer iimax, ist, ipn, ip, i, j, percent, level
      real xmax0, xmax, fp, pmax
      real*8 p
      character uprow*3, filename*30

      integer irebinnum, irebinmax, irebin
      parameter (irebinmax=5)
      integer irebinval(irebinmax)
      data irebinval/16,32,64,128,256/
      
      
      integer istart, ntopmax, ntop, delta, status, nused
      parameter (ntopmax=100)
      integer idx(ntopmax), next(ntopmax), bottom
      real fracdif, ptop(ntopmax,irebinmax), snrtop(ntopmax,irebinmax)
      real supertopp, supertopsnr

      integer nfold
      parameter (nfold=128)
      real fold(nfold), fold8(8,nfold), rescale
      byte pattern(2,nfold), bfold(nfold)

      data istart/1/   ! start searching for candidates at bin 1
      data ntop/10/   ! print ten but save three candidates
      data fracdif/0.002/ ! fractional difference between candidates
      data delta/2/    ! separation


      uprow=char(27)//'[A'     

      if (fixsamp) then ! fix sampling time for Parkes data
         rescale = 1.0140375
         isamp=20000
         write(*,*) '*** Sampling time fixed to: ', isamp
         write(*,*) '*** Later rescaling...'         
      else
         rescale = 1.0
c         write(*,*) isamp
      endif
      
      do irebin=1, irebinmax
         do i=1, ntop
            ptop(i,irebin)=0.0
            snrtop(i,irebin)=0.0
         enddo
      enddo         
      
      if (searchw) then
         irebinnum=irebinmax
      else
         irebinnum=1
         irebinval(irebinnum) = np
      endif
      
      supertopsnr=0.0
      
      do irebin=1, irebinnum

         init = .true.
         xmax0=-99999999.0
         xmax=-99999999.0
         ip1=ip10
         ip2=ip20

         np = irebinval(irebin)

         if ((ndata/np).gt.ndataMax) then
            write(*,*) 'Rebinning too small for arrays in RUNFFA!'
            write(*,*) ndata, np, ndataMax
            stop
         endif

         call rebin(np,ndata,timeIn,time)
         write(*,*) 'Rebinned by: ',np


c      do i=1, ndata/np
c         write(88,*) i, time(i)
c      enddo
c      write(*,*) ip1, ip2
      
         iimax = ndata / np
         fp = 1000.0 / float(isamp) / float(np)
         ip1 = int( float(ip1) * fp)
         ip2 = int(float(ip2) * fp) + nint(max(0.0, (fp - 1.0)))
c         print*, ip1, ip2, fp
         
         write(*,*)

         ist = 0
         ipn = ip2 - ip1 + 1
         level=1
         do ip = ip1, ip2
            percent=int(100.*float(ip-ip1)/float(ip2-ip1))
            if (percent.gt.level) then
               write(*,33) uprow, percent
 33            format(a,' Processed: ',i3,'%')
               level=level+1
            endif
            call ffa(ndata,time,iimax, ip, fp, ipn, p, ist, xmax, 
     +           nprof, signal,nspec,pspec,snspec,init)
            init = .false.
            if (xmax. gt. xmax0) then
               pmax = p
               xmax0 = xmax
            endif
         enddo
         nused = ist
         
         call supertopvals(snspec,nused,istart,ntop,fracdif,delta,idx,
     &        next, bottom, status)
         
         if (dumpspec) then
            write(filename,121) irebin
 121        format('FFAspec.',i2.2)
            write(*,*) 'Dumping spectrum...',filename
           call glun(iunit)
           open(iunit,file=filename,form='formatted',status='unknown')
           do i=1, nused
              write(iunit,*) pspec(i), snspec(i)
           enddo
           close(iunit)
         endif
         
         j=bottom      
         do i=1,ntop
            ptop(ntop-i+1,irebin)=pspec(idx(j))
            snrtop(ntop-i+1,irebin)=snspec(idx(j))
            j=next(j)
         enddo 
         if (snrtop(1,irebin).gt.supertopsnr) then
            supertopsnr=snrtop(1,irebin)
            supertopp=ptop(1,irebin)
         endif
         write(*,*) 'Best candidate this fold:',ptop(1,irebin),
     +        '  S/N: ',snrtop(1,irebin)

      enddo
c     
c     write out a .prd file for use by best
c     
      call glun(iunit)
      fname=fname(1:index(fname,'.tim')-1)//'FFA.prd'
      if (append) then
         open(iunit,file=fname,status='unknown',access='append')
      else
         open(iunit,file=fname,status='unknown')
      endif
      write(iunit,*) 'DM:',refdm,' AC:',0.0,' AD:',0.0
      do i=1, ntop
         write(iunit,'(5(f5.1,1x,f13.8,1x))') 
     &   (snrtop(i,irebin),ptop(i,irebin)*rescale,irebin=1,irebinmax)
      enddo
      close(iunit)
         
      write(*,*) 'Best suspect',supertopp*rescale,' ms'
      write(*,'(x,a,x,f5.1)') 'S/N:',supertopsnr

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ffa(ndim,time,mmax,ip,fp,ipn,pmax,ist,xmax,
     +               nprof,signal,nspec,pspec,snspec,init)

c
c* Copyright (C) 1995 by Peter M"uller, MPIfR.  peter@mpifr-bonn.mpg.de
c*
c* Permission to use, copy, modify, and distribute this software and its
c* documentation for any purpose and without fee is hereby granted, provided
c* that the above copyright notice appear in all copies and that both that
c* copyright notice and this permission notice appear in supporting
c* documentation.  This software is provided "as is" without express or
c* implied warranty.
c
c     modified for PMsurv, M.Kramer Dec, 2001
c
c
      implicit real*4 (a-h,o-z)
      
      integer ndim, nspec, nprof
      real time(ndim), pspec(nspec), signal(nprof)
      real snspec(nspec)

      
      integer npuls
      parameter (npuls=2**24)
      real puls(npuls), lp

      real*8 pmax, p, p0, p1, p2, p3, p4, p5, ps, dp, xp 
      real*8 x


      real lpfirst
      logical first, init
      save first, lpfirst


      data first/.true./

      snr = -999999999.0

      in = 1
      ip2 = ip/2
      ip10 = max(1.0, ip * 0.15)
      nmax = ip * 2**int(log(float(mmax)/float(ip)) / log(2.0) + 0.9)
      if(nmax .gt. NDIM) nmax = nmax/2
      lmax = nmax / (ip+1)
      lmax1 = lmax - 1
      isnr = max(1, int(float(lmax) / max(1.0, 1024.0 / float(ipn))))
      kp = nint(log(float(lmax+1)) / log(2.0))
      dp = dfloat(ip + 1) / dfloat(nmax)
      p1 = 1.d0 / dfloat(ip)
      p2 = p1 * p1
      p3 = p2 * p1
      p4 = p3 * p1
      p5 = p4 * p1
      ps = p1 + p2 + p3 + p4 + p5

      if (init) then
         do i=1, npuls
            puls(i)=0.0
         enddo
         do i=1, nspec
            snspec(i)=0.0
         enddo
      endif

      time(mmax) = 0.0
      do n = 0, lmax1, in
         x = dp * n
         xp = (((((p5*x + p4)*x + p3)*x + p2)*x - ps)*x + p1)*x
c        p = ip + x - dp * nint(xp/dp - 0.5*mod(ip,2))
         p = ip + x - xp
         do l = 0, kp-1
            if (mod(n, 2**l) .eq. 0) kp1 = kp - l 
         enddo
         do k = kp1, kp
            np = 2**(kp-k)
            joff = nmax - nmax / 2**(k-1)
            ioff = nmax - nmax / 2**max(k-2,0)
            ish = mod((n + np) / np / 2, ip)
            do i = 0, np-1
	       iip = i * ip
               i0 = iip + joff
               i1 = 2 * iip + ioff
               i2 = i1 + ip
               do j = 1, ip
                  j1 = j + ish
                  if(j1 .gt. ip) j1 = j1 - ip
                  ind = j + i0
                  jnd = j + i1
                  knd = j1 + i2
                  if(kp1 .eq. 1 .and. k.eq.1) then
                     puls(ind) = time(jnd) + time(knd)
                  else
                     puls(ind) = puls(jnd) + puls(knd)
                  endif
               enddo
            enddo
         enddo

         xmax0 = -999999999.0
         do lm = 1, ip
            if(xmax0 .lt. puls(lm+joff)) then
               xmax0 = puls(lm+joff)
               imax = lm
            endif
         enddo
         sx = 0.0
         sxx = 0.0
	 lp = 0.0
         do im = imax+ip10, imax+ip2, 1
	    lp = lp + 1.
            lm = mod(im, ip)
            if(lm .eq. 0) lm = ip
            sx = sx + puls(lm+joff)
            sxx = sxx + puls(lm+joff)**2
         enddo
         do im = imax-ip10, imax-ip2, -1
	    lp = lp + 1.
            lm = mod(im+ip, ip)
            if(lm .eq. 0) lm = ip
            sx = sx + puls(lm+joff)
            sxx = sxx + puls(lm+joff)**2
         enddo
         rms = sqrt((lp*sxx - sx*sx) / (lp*(lp-1)))

c       -- updated. may 6 - mk, drl
c         rms = sqrt((lp*sxx - sx*sx)) / lp
c         if (first) then
c            first=.false.
c            lpfirst=lp
c         else
c            rms=rms*sqrt(lp/lpfirst)
c         endif
c     ------------------------

         if(abs(rms) .lt. 0.001) rms = 1.0
c
         xmax0 = (xmax0 - sx/lp) / rms
c
         if(xmax0 .gt. snr) then
	    snr = xmax0
	    p0 = p/fp
	 endif
         
c         write(*,*) rms, p0, lp,ip
c         pause

         if(mod(n+1, isnr).eq.0 .or. n.eq.lmax1) then
            ist = ist + 1
            if (ist.gt.nspec) then
               write(*,*) 'Long period buffer overflow in FFA!'
               write(*,*) ist, nspec
            else
               pspec(ist)=p0
               snspec(ist)=snr
            endif           
            snr = -999999999.0 
         endif
         if(xmax0 .gt. xmax) then
            xmax = xmax0
            pmax = p/fp
            if (ip.gt.nprof) then
               write(*,*) 'Long period buffer#2 overflow in FFA!'
            else               
               do lm = 1, ip
                  signal(lm) = (puls(lm+joff) - sx/lp) / float(lmax)
               enddo
            endif
         endif
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rebin(np,nmax,timeIn,time)
c
c* original 1995 by Peter M"uller, MPIfR.  peter@mpifr-bonn.mpg.de

c  modified for PMsurv, M.Kramer, Dec 2001
c
      implicit none

      integer np, nmax
      real timeIn(nmax), time(nmax/np)

      integer ind, n, i, newnmax
      real sum, ssum


      if(np .le. 1) then
         do i=1, nmax
            time(i)=timeIn(i)
         enddo
         return
      endif

c      print*,'REBIN: ',np, nmax

      newnmax=nmax/np

      ssum=0

      ind = 0
      do n = 1, nmax, np
         sum = 0.0
         do i = 1, np
            sum = sum + timeIn(n + i - 1)
         enddo
         ind = ind + 1
         
         if (ind.le.newnmax) then
            time(ind) = sum
            ssum=ssum+sum
         else
            write(*,*) 'REBIN: Array exceeded! ', ind
         endif

c     CHECK for offsets from zero!
c        xx = ind / (float(nmax)/float(np))
c        write(27,*) xx, sum

      enddo
      
c     - subtract non-zero mean -
      ssum=ssum/(float(nmax)/float(np))
      do i=1, newnmax
         time(i)=time(i)-ssum
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @(#)supertopvals.f	3.1 12/17/92
      subroutine supertopvals(dat, ndat, istart, ntop,fracdif,delta,idx,
     &    next, bottom, status)
      
c     finds the ntop highest values in dat(ndat) places pointers to these
c     values in idx(ntop). the ordering of the ntop values can be recovered
c     from bottom (which points via the corresponding element in idx to the
c     lowest value of the top ntop) the next highest value pointer is pointed
c     to by next. the elements may thus be printed in ascending order by the
c     following code fragment.
c     modified version of topvals pah 23-may-1991 - will only return
c     'peaks' if they are the fractional difference if greater then 
c     fracdif points apart.
c     
c     
c     91/09/17 pah corrected an initialization bug. on entry the list is made 
c     to consist of entirely the first element.
c     and again 91/11/24 - also stop adjacent bins ever having peaks found
      
c     sample fortran to recover the list of values
c     
c     iptr=bottom      
c     do i=1,ntop
c     write(*,*) dat(idx(iptr))
c     iptr=next(iptr)
c     enddo 
      
      
      integer ndat,ntop,status,istart
      integer idx(ntop),next(ntop),bottom,i,iptr,ilast,delta
      integer itemp, idxprev
      real dat(ndat),minval, fracdif
      
c     test if ndat is smaller than ntop panic!

         status = 0
         if (ndat.lt.ntop) then
            status = -1
            return
         end if

c  fill next and initialise bottom and next(ntop) make 
c  the first point the highest
         do i=1,ntop
            next(i)=i+1 
            idx(i)=istart
         end do
         next(ntop) = 0
         bottom = 1
         idxprev = istart
         minval = dat(istart)

c start main loop
         do i=istart+1, ndat
            if (dat(i).gt.minval) then

c first check that the point is further from previously found point
               if(real(i-idxprev)/i.gt.fracdif
     &             .and.i-idxprev.gt.delta) then

c put the previous point into the list
                   ilast = bottom
                   iptr  = next(bottom)
                   do while(dat(idxprev).gt.dat(idx(max(1,iptr))).and.
     :                                                (iptr.ne.0))
                       ilast = iptr
                       iptr = next(iptr)
                   end do
                   idx(bottom)=idxprev
                   if(ilast.ne.bottom)then
                     itemp = next(bottom)
                     next(bottom)=iptr
                     next(ilast)=bottom
                     bottom=itemp
                   endif
                   minval = dat(idx(bottom))
                   idxprev=i
               elseif(dat(i).gt.dat(idxprev)) then

c make the current i then maximum found in this nres group
                   idxprev=i
               endif
            end if
         end do

c there might still be 1 idxprev left to be put into the list
         if(dat(idxprev).gt.minval) then
            ilast = bottom
            iptr  = next(bottom)
            do while(dat(idxprev).gt.dat(idx(max(1,iptr))).and.
     :                                            (iptr.ne.0))
                ilast = iptr
                iptr = next(iptr)
            end do
            idx(bottom)=idxprev
            if(ilast.ne.bottom)then
              itemp = next(bottom)
              next(bottom)=iptr
              next(ilast)=bottom
              bottom=itemp
            endif
         endif

         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
