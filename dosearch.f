c=============================================================================
      subroutine dosearch(llog,dump,rspc,oldwhite,pmzap,mmzap,
     &                    recon,prdh,sfile,pmax)
c=============================================================================
c
c     Does the search for the pulsars in the Frequency domain.
c
c     - form power spectrum
c     - zaps any interfering bins in spectrum
c     - normalises resulting spectrum to zero mean and unit rms
c     - harmonic summing for search for strongest peaks 
c     - Profile reconstruction added. .prd (only) file contains additional
c       SNR from reconstructed profiles. R.E. 2/11/06
c     - got rid of .frq output
c
c=============================================================================
      implicit none
      include 'seek.inc'
      include 'csamp.inc'
      real junk,minsamp,maxsamp
      integer llog
      logical dump,rspc,pmzap,mmzap,prdh,filex
      integer oldwhite
      character*80 sfile
c
c     Local variables
c      
      integer h,i,j,k,npf,fold,fb,lun,nc,ncal,istat,n
      real snrbest,rms,sumsq,snrc,thresh,fnyq,spcsnr
      real saverms
      integer indx(npts/8),snum(npts/8)
      integer top,nsm,nf1,c(nfolds),cmin,cmax
      parameter(top=1024,nsm=1024*8)
c      parameter(top=4096,nsm=4096*8)
      real rmea(nsm),rrms(nsm),fastmea(npts),sres,flo,fhi,fhr,smax
c      real ralphmea(npts)   ! to  check what effect the running mean has on snr
      real*8 freq,pmax,fbest,ratio,pc(top)
      real*8 freqff
      logical harm1,harm2,harm3,harm4,harm5,harm6,recon
      real*8 pcand(nfolds,top),ptmp
      real scand(nfolds,top),sc(top),spc(512),
     &     temp_snrecon(nfolds,top),prof_peak,cpower
      real snfact
c      parameter(snfact=1.1284) ! = sqrt(4/pi)
      parameter(snfact=0.5) ! = sqrt(4/pi)
      integer icand(nfolds,top)
      integer nav,ntot,fbin,blo,bhi,nharm
      character*80 pfile,ffile,hsum,hsums
      integer kk,lhsums,length

      fold=0
      if (rspc) then
         call readspec(sfile,fold,samp,refdm,refac,tsamp,npf)
         nf1=real(tsamp)*2*npf/real(pmax)
      else
         nf1=real(tsamp)*ntim/real(pmax)
         write(llog,*) 'Forming amplitude spectrum. (Pmax=',pmax,' s!)'
         call formspec(npf,nf1)
         sres=real(freq(tsamp,npf,1,2))-real(freq(tsamp,npf,1,1))
         write(llog,*) 'Raw spectral resolution: ',sres*1000,' mHz)'
      endif

      fnyq=0.5/real(tsamp)
      write(llog,*) 'Nyquist frequency:',fnyq,' Hz'

      pfile="zerodm.spc"
      inquire(file=pfile,exist=filex)
      if (refdm.gt.0.0.and.filex) then
         call readspec(pfile,fold,zerodm,junk,junk,tsamp,npf)
         maxsamp=0.0
         do i=1,npf
            samp(i)=samp(i)-zerodm(i)
            maxsamp=max(maxsamp,samp(i))
c            if (samp(i).lt.0.0) samp(i)=0.0
         enddo
c         minsamp=-1.0*maxsamp
c         do i=1,npf
c            if (samp(i).lt.minsamp) samp(i)=minsamp
c         enddo
c         write(*,*) maxsamp,minsamp
      endif

      if (dumpraw) then
	  fold=0
	  if (sfile.ne.' ') then
	    pfile=sfile
	  else
            write(pfile,'(a,i1,a)') filename(1:lst)//'_',fold,'.spc'
	  endif
          call writespec(llog,pfile,fold,samp,refdm,refac,tsamp,npf)
	  stop
      else if (dump) then
c         write(pfile,'(a4,i1,a4)') 'fold',0,'.spc'
         write(pfile,'(a,i2.2,a4)') filename(1:lst)//'_',fold,'.spc'
         call writespec(llog,pfile,fold,samp,refdm,refac,tsamp,npf)
      endif 
c
c     Spectral mask
c
      if (maskfile(1).ne.' ') then
        if (nmasks.gt.1) then
           write(llog,*) 'Reading spectral mask...'
        else
           write(llog,*) 'Reading spectral masks...'
        endif
        call glun(lun)
	istat=0
	open(unit=lun,file=maskfile(1),status='old',iostat=istat)
        if (istat.ne.0) then
           write(*,*) 'WARNING - mask file not found...'
        else
           j=0
           do while(istat.eq.0)
              read(lun,*,iostat=istat) i
              if (istat.eq.0) then
                 samp(i)=0.0
                 j=j+1
              endif
           enddo
           write(llog,*) 'Masked',j,' spectral bins from fold 1'
         endif
         close(unit=lun)
      endif
c
c     Zap birdies before doing any spectral manipulation
c        
      if (zapfile.ne.' ') call zapit(llog,1,zapfile,samp,npf,tsamp)

      write(llog,*) 'Whitening spectrum...'

      nav=max(128,npf/nsm)
      if (oldwhite.eq.1) then

         write(llog,*) 'Calculating spectral mean/rms every',
     &   nav,' bins...',real(freq(tsamp,npf,1,nav+1))
     &   -real(freq(tsamp,npf,1,1)),' Hz'

         call getrmea(samp,npf,nav,rmea,ncal)
         call getrrms(samp,npf,rmea,nav,rrms)

      else if (oldwhite.eq.0) then

         write(llog,*) 'Calculating AGL mean and rms every',
     &   nav,' bins...',real(freq(tsamp,npf,1,nav+1))
     &   -real(freq(tsamp,npf,1,1)),' Hz'

         call getmeanrms(samp,npf,nav,rmea,ncal,rrms)

      else if (oldwhite.eq.2) then

         write(llog,*) 'Calculating spectral median every',
     &   nav,' bins...',real(freq(tsamp,npf,1,nav+1))
     &   -real(freq(tsamp,npf,1,1)),' Hz'

         call getrmed(samp,npf,nf1,nav,rmea,ncal)

      else
         stop 'Do not know how to whiten spectrum...'

      endif

      
      call rotate_time(series,npf,nf1)

      if (fbrute.gt.0.0) 
     &   write(llog,*) 'Brutal zapping below',fbrute,' Hz!'
      n=0
      sumsq=0.0
      j=0
      do i=1,npf
         h=min(ncal,i/nav+1)

         if (samp(i).ne.0.0) then

            if (oldwhite.le.1) then ! new and old mean -> same normalization
c
c             Subtracts running mean and scale it so that the rms=1
c

               if (rrms(h).ne.0.0) samp(i)=(samp(i)-rmea(h))/rrms(h) 
               if (rmea(h).eq.0.0) samp(i)=0.0

               if (rrms(h).eq.0.0) then
                  series(2*i-1) = 0
                  series(2*i) = 0
               else
                  series(2*i-1)=series(2*i-1)/rrms(h)
                  series(2*i)=series(2*i)/rrms(h) 
               endif               
               
            else
c
c              New method (default) is to divide by running median
c              to minimize biases and then subtract unity from the result
c
	       if (rmea(h).ne.0.0) samp(i)=(samp(i)/rmea(h))-1.0
	       if (rmea(h).eq.0.0) samp(i)=0.0

               if (rmea(h).eq.0.0) then
                  series(2*i-1) = 0
                  series(2*i) = 0
               else
                  series(2*i-1)=series(2*i-1)/rmea(h)
                  series(2*i)=series(2*i)/rmea(h)
               endif

            endif

         endif

         if (mod(i,1024).eq.0.and.samp(i).lt.3.0) then
            n=n+1
            sumsq=sumsq+samp(i)
         endif

c
c this line is the -b option which brutally zaps all RFI and psrs < fbrute Hz
c changed for 47tuc analysis Mar 17, 2003
	 if (fbrute.gt.0.0.and.freq(tsamp,npf,1,i).lt.fbrute.and.
     &      samp(i).gt.10.0) samp(i)=1.0
      enddo

      sumsq=sumsq/real(n)
      rms = 0
      do i=1, npf
         if (mod(i,1024).eq.0.and.samp(i).lt.3.0) then
             rms = rms + (samp(i)-sumsq)**2
         endif

c         write(88,*) i, samp(i), series(2*i-1), series(2*i)

      enddo

      rms=sqrt(rms/real(n))

      write(llog,*) 'Resulting spectral RMS:',rms

c     save rms for later recall
      saverms = rms

c
c     original spectrum + 4 harmonic sums (Lyne/Ashworth code)
c      
      write(llog,759) 'Harmonic sums are: '
      lhsums=0
      do fold=1,nfolds
         write(llog,762) foldvals(fold)
         write(hsum,'(i4)') foldvals(fold)
         hsums=hsums(1:lhsums)//hsum(1:4)
         lhsums=lhsums+4
      enddo
 759  format(1x,a,$)
 762  format(1x,i4,2x,$)

      write(llog,*)


      write(llog,*) 'Doing harmonic summing...'
      if (nfolds.eq.5.and.foldvals(1).eq.1.and.foldvals(2).eq.2
     &    .and.foldvals(3).eq.4.and.foldvals(4).eq.8.and.
     &    foldvals(5).eq.16) then
         write(llog,*) 'Lyne-Ashworth harmonic summing'
         call oldsumhrm(samp,npf,nf1) ! Lyne-Ashworth code
      else
         write(llog,*) 'DB\'s slow-but-simple harmonic summing routine'
         call sumhrm(samp,npf,nf1)    ! David Barnes code
      endif
c      
c     mask out folds 2-nfolds in addition to fold 1 (done above) if 
c     mask files are present (drl - 28/04/05)
c
      if (nmasks.eq.nfolds) then
         do fold=2,nfolds
           j=0
           call glun(lun)
           istat=0
           open(unit=lun,file=maskfile(fold),status='old',iostat=istat)
           do while(istat.eq.0)
              read(lun,*,iostat=istat) i
              if (istat.eq.0) then
                 power(fold,i)=0.0
                 j=j+1
              endif
           enddo
           write(llog,*) 'Masked',j,' spectral bins from fold',fold
           close(lun)
         enddo
      endif

      fb=0
      fbest=0.0
      snrbest=0.0
      ntot=0
c
c     Search for candidates over all harmonic folds
c
      write(llog,*) 'Doing harmonic searching...'
      do fold=1,nfolds
         rms = saverms * sqrt(real(foldvals(fold)))

c
c        for Parkes Data - call zapping algorithms if selected
c
        if (mmzap.or.pmzap) then
           do i=1,npf
              samp(i)=power(fold,i)
           enddo
           if (pmzap) 
c     &     call zap_pmbrd(samp,npf,nf1,tsamp*1000.0,rms,0.0,fold)
     &     call zap_pmbrd(samp,series,npf,nf1,tsamp*1000.0,rms,0.0,fold)
           if (mmzap) 
c     &     call zap_mmbrd(samp,npf,nf1,tsamp*1000.0,rms,0.0,fold)
     &     call zap_mmbrd(samp,series,npf,nf1,tsamp*1000.0,rms,0.0,fold)
           do i=1,npf
              power(fold,i)=samp(i)
           enddo
        endif
        c(fold)=0
        thresh = 5.0 
 5      do i=1,npf
c          ptmp=1.0/freq(tsamp,npf,fold,i)
          ptmp=1.0/freqff(tsamp,npf,foldvals(fold),i)
          if (ptmp.gt.pmax) power(fold,i)=0.0 ! Zap P>Pmax signals

          snrc=power(fold,i)/rms 

          if (snrc.gt.thresh.and.c(fold).lt.top) then
             c(fold)=c(fold)+1
             samp(c(fold))=power(fold,i)
             snum(c(fold))=i
c          else if (snrc.gt.thresh) then
c             write(*,*) 'WARNING: not enough candidates saved!!'
          endif
        enddo
        if (c(fold).lt.2) then
           thresh=thresh-1.0
           goto 5
        endif
c
c       Sort amplitude spectrum in s/n order using the dreaded indexx
c       routine from numerical recipes
c
        call indexxf77(c(fold),samp,indx)
c     
        j=0
        do i=c(fold),1,-1
          j=j+1
c          pcand(fold,j)=1000.0/freq(tsamp,npf,fold,snum(indx(i)))
          icand(fold,j) = snum(indx(i))
          pcand(fold,j)=1000.0/freqff(tsamp,npf,foldvals(fold),
     &         snum(indx(i)))
          snrc=samp(indx(i))/rms   + real(foldvals(fold))/rms
          scand(fold,j)=snrc  
          if (snrc.gt.snrbest) then
            snrbest=snrc
            fbest=1000.0/pcand(fold,j)
            fb=foldvals(fold)
          endif
        enddo

        do i=1,c(fold)
           do j=1,c(fold)
              if (j.ne.i.and.scand(fold,j).gt.0.0.and.
     &           scand(fold,i).gt.0.0) then
                 ratio=pcand(fold,i)/pcand(fold,j)
                 if (ratio.lt.1.0) ratio=1.0/ratio
                 ratio=ratio-int(ratio)
                 harm1=ratio.gt.0.999.or.ratio.lt.0.001
                 harm3=ratio.gt.0.333.and.ratio.lt.0.334
                 harm5=ratio.gt.0.499.and.ratio.lt.0.501
                 harm6=ratio.gt.0.666.and.ratio.lt.0.667
                 if(harm1.or.harm3.or.harm5.or.harm6)scand(fold,j)=0.0
c                 harm1=ratio.gt.0.999.or.ratio.lt.0.001
c                 harm2=ratio.gt.0.249.and.ratio.lt.0.251
c                 harm3=ratio.gt.0.333.and.ratio.lt.0.334
c                 harm4=ratio.gt.0.499.and.ratio.lt.0.501
c                 harm5=ratio.gt.0.666.and.ratio.lt.0.667
c                 harm6=ratio.gt.0.749.and.ratio.lt.0.751
c                 if(harm1.or.harm2.or.harm3.or.harm4.or.harm5.or.harm6)
c     &              scand(fold,j)=0.0
              endif
           enddo
        enddo

        j=0
        do i=1,c(fold)
           if (scand(fold,i).gt.thresh) then  
              j=j+1
              scand(fold,j)=scand(fold,i)
              pcand(fold,j)=pcand(fold,i)
           endif
        enddo
        c(fold)=j

        do i=1,c(fold)
           ntot=ntot+1
           if (ntot.le.top) then
             sc(ntot)=scand(fold,i)
             pc(ntot)=pcand(fold,i)
           else
c             write(*,*) 'WARNING: master candidate list full!'
           endif
        enddo
        
        if (dump) then
          do i=1,npf
            samp(i)=power(fold,i)
          enddo
c         write(pfile,'(a,i1,a)') filename(1:lst)//'_',fold,'.spc'
c         write(pfile,'(a,a)') filename(1:lst),'.spc'
c         write(pfile,'(a4,i1,a4)') 'fold',fold,'.spc'
c         call writespec(llog,pfile,fold,samp,refdm,refac,tsamp,npf)
          write(pfile,'(a,i2.2,a4)') filename(1:lst)//'_',fold,'.spc'
          call writespec(llog,pfile,fold,samp,refdm,refac,tsamp,npf)
        endif
        if (recon) then
c      
c         Do profile reconstruction. R.E.     
c 
          nharm=2**(fold-1)  
          do i=1, top        
c             temp_snrecon(fold,i)=sqrt(cpower(npf,foldvals(fold),
c     &       pcand(fold,i)))/saverms
             call recon_prof(llog,npf,fold,pcand(fold,i),i,prof_peak)
             temp_snrecon(fold,i)=prof_peak/sqrt(real(nharm))*snfact
          enddo
       endif
      enddo

      cmin=top
      do i=1,nfolds
         cmin=min(cmin,c(i))
      enddo
      
      cmax=0
      do i=1,nfolds
         cmax=max(cmax,c(i))
      enddo
      
      call glun(lun)
      pfile=filename(1:lst)//'.prd'
      open(unit=lun,file=pfile,status='unknown',access=facc)
      if (facc.ne.'append'.and.prdh) then
         write(lun,'(a)') '##BEGIN HEADER##'
         write(lun,'(a,a)') 'SOURCEID = ',source_name
         write(lun,'(a,f8.1,a)') 'FREF = ',fref,' MHz'
         write(lun,'(a,f12.4)') 'TSTART = ', mjdstart
         write(lun,'(a,a)') 'TELESCOPE = ', telname
         write(lun,'(a,a)') 'RAJ = ', ra
         write(lun,'(a,a)') 'DECJ = ', dec
         write(lun,'(a,f8.1,a)') 'TSAMP = ', tsamp*1.0e6, ' us'
         write(lun,'(a)')   'PROGRAM = SEEK'
         write(lun,'(a,a)') 'VERSION = ', version(length(version)-2:)
         write(lun,'(a,a)') 'HARM_FOLDS = ', hsums(1:lhsums)
         if (recon) then
            write(lun,'(a,a)') 'COLS = SNR_SPEC SNR_RECON PERIOD'
         else
            write(lun,'(a,a)') 'COLS = SNR_SPEC PERIOD'
         endif
         write(lun,'(a)')   '##END HEADER##'
      endif
      write(lun,*) 'DM:',refdm,' AC:',refac,' AD:',refad
      do i=1,cmax
         do j=1,nfolds
            if (i.gt.c(j)) then
              pcand(j,i)=1.0
              scand(j,i)=0.0
              temp_snrecon(j,i)=0.0
              icand(j,i)=0
            endif
         enddo
         do kk=1,nfolds
            if (.not.recon) then
               write(lun,71) scand(kk,i),pcand(kk,i)
            else
               write(lun,72) scand(kk,i),temp_snrecon(kk,i),pcand(kk,i)
            endif
c            write(lun,81) scand(kk,i),pcand(kk,i),icand(kk,i)
         enddo
         write(lun,*)
c         write(lun,1) scand(1,i),pcand(1,i),scand(2,i),pcand(2,i),
c     &                scand(3,i),pcand(3,i),scand(4,i),pcand(4,i),
c     &                scand(5,i),pcand(5,i)
      enddo
      close(unit=lun)

 1    format(5(f7.1,1x,f13.8,1x))
 2    format(f7.1,1x,f9.4,3x,f7.3,1x,f5.1)
c 71   format(f7.1,1x,f13.8,1a,$)
 71   format(f8.1,f14.8,$)
 72   format(f8.1,f8.1,f14.8,$)
 81   format(f7.1,1x,f13.8,1x,i8,1x,$)
      
      write(llog,*) 'Best suspect:',1000.0/fbest,' ms'
      write(llog,'(x,a,x,f5.1)') 'S/N:',snrbest
      write(llog,*) 'Found peak at:',fbest,' Hz'
      write(llog,*) 'Number of harmonics:',fb

      ffile=filename(1:lst)//'.top'
      open(unit=lun,file=ffile,status='unknown',access=facc)
      write(lun,*) 1000.0/fbest,snrbest,refdm
      close(unit=lun)

      end
c=============================================================================
