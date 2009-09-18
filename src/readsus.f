c==============================================================================
      subroutine readsus(filename, fold, mc, nc, mt, nt,
     &                   par, snr, trid, fld, dm, ac, ad)
c==============================================================================
c
c     Reads in a list of suspects from a FIND output file. The signals
c     are either frequencies or periods. This is however irrelevant as
c     far as this routine is concerned.
c
c     filename - ch - name of the suspect file
c     fold     - i4 - harmonic fold number to read in (1-5 inclusive)
c     mc       - i4 - maximum number of suspects allowed
c     nc       - i4 - number of candidates read in
c     mt       - i4 - maximum number of trial loops allowed
c     nt       - i4 - number of trial loops read in
c     par(mc)  - r4 - array containing the suspects for this fold
c     snr(mc)  - r4 - array containing the suspect signal-to-noise ratios
c     trid(mc) - i4 - array containing the trial (loop) index for each suspect
c     fld(mc)  - i4 - array containing the fold for each suspect
c     dm(mt)   - r4 - array containing the dm value for each list
c     ac(mt)   - r4 - array containing the ac value for each list      
c     ad(mt)   - r4 - array containing the ad value for each list      
c
c     (dunc@mpifr-bonn.mpg.de - November 1997)
c      
c     Modification history:
c      
c     98/07/01 -> dunc@mpifr-bonn.mpg.de (read in ac array)
c     98/07/13 -> dunc@mpifr-bonn.mpg.de (read in ad array)
c     98/07/15 -> dunc@mpifr-bonn.mpg.de (changed calling sequence)
c     98/11/24 -> dunc@naic.edu          (trid now correct for multiple calls)
c      
c==============================================================================
c
      implicit none
      character*(*) filename
      integer fold, nc, nt, mc, mt
      integer trid(mc),fld(mc)
      real*8 par(mc)
      real snr(mc),dm(mt),ac(mt),ad(mt)
      integer lun,ltmp,istat,i,idx,jj
      real s(5)
      real*8 p(5)
      character*400 line

      idx=0
      call glun(lun)
      open(unit=lun,file=filename,status='old',err=999)
      do while(nc.lt.mc.and.nt.lt.mt)
         read(lun,'(a)',end=998) line
         if (index(line,'DM').gt.0) then
            nt=nt+1
            idx=idx+1
            ltmp=index(line,'DM:')+3
            read(line(ltmp:),*) dm(nt)
            ltmp=index(line,'AC:')+3
	    read(line(ltmp:),*) ac(nt)
            ltmp=index(line,'AD:')+3
            ad(nt)=0.0
            if (ltmp.gt.3) read(line(ltmp:),*) ad(nt)
         else 
            nc=nc+1
            fld(nc)=fold
            if (index(filename,'.top').gt.0) then
              read(line,*) s(1),p(1),s(2),p(2)
              par(nc)=p(1)
              snr(nc)=s(1)
            else
c	      do i=1,5
c		p(i)=0.0
c		s(i)=0.0
c	      enddo
c              read(line,*,iostat=istat)
c     &        s(1),p(1),s(2),p(2),s(3),p(3),s(4),p(4),s(5),p(5)
c              par(nc)=p(fold)
c              snr(nc)=s(fold)
c               do jj=1,fold
c                  read(line,*,iostat=istat) s(1),p(1)
c               enddo
               read(line(1+(fold-1)*22:),*,iostat=istat) s(1),p(1)
               par(nc)=p(1)
               snr(nc)=s(1)
            endif
            if (par(nc).eq.1.0.or.snr(nc).eq.0.0) then
               nc=nc-1
            else
               trid(nc)=idx
            endif
         endif
      enddo
c 1    format(5(f7.1,1x,f13.8,1x))
c 1    format(5(f5.1,1x,f13.8,1x))
c 1     format(5(f8.1,f14.8))
 998  close(unit=lun)
      return
 999  write(*,'(a)') 'Error handling file: ',filename
      stop
      end 
c      
c==============================================================================
