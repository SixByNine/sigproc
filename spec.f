      program spec
      implicit none
      include 'vers.inc'
      integer mp,mc
      parameter(mp=2**23,mc=20)
      real fbin(mp),samp(mp),smin,smax,lo,hi
      real snr,spcsnr,dm,ac,sum,sumsq,sigma,mean
      real*8 tsamp,freq,fmark
      real histmax,nmax,ppsr,fpsr,deltaf,fhar,samp2(mp),samp3(mp)
      integer narg,iargc,fold,npf,i,j,nbins,nh,ih,np,rbin,nadd
      logical stats,dump,filex,power
      character*80 comline,title,sfile,pgdev
      sfile=' '
      pgdev='?'
      fold=1
      narg=iargc()
      lo=0.0
      hi=-1.0
      stats=.false.
      nbins=25
      ppsr=0.0
      fpsr=0.0
      rbin=0
      dump=.false.
      power=.false.
      fmark=0.0
      if (narg.gt.0) then
         do i=1,narg
           call getarg(i,comline)
           inquire(file=comline,exist=filex)
           if (filex) sfile=comline
           if (index(comline,'-f').gt.0) read(comline(3:),*) fold
           if (index(comline,'-l').gt.0) read(comline(3:),*) lo
           if (index(comline,'-h').gt.0) read(comline(3:),*) hi
           if (index(comline,'-s').gt.0) stats=.true.
           if (index(comline,'-n').gt.0) read(comline(3:),*) nbins
           if (index(comline,'-p').gt.0) read(comline(3:),*) ppsr
           if (index(comline,'-r').gt.0) read(comline(3:),*) rbin
           if (index(comline,'-d').gt.0) dump=.true.
           if (index(comline,'-d').gt.0) pgdev='/null'
           if (index(comline,'-P').gt.0) power=.true.
           if (index(comline,'-F').gt.0) read(comline(3:),*) fmark
         enddo
      else
         write(*,*)
         write(*,*) 'SPEC: ',version
         write(*,*) 'Displays spectrum files produces by FIND'
         write(*,*)
         write(*,*)'usage: spec <SPECTRUM_FILE> -{options}'
         write(*,*)
         write(*,*)'The spectrum file is produced by running find -s'
         write(*,*)'This program is usually run on fold1.spc'
         write(*,*)
         write(*,*)'Available options are:'
         write(*,*)
         write(*,*)'-s: run stats mode'
         write(*,*)'-d: dump ASCII version of spectrum to file "dump"'
         write(*,*)'-P: plot powers (def=amplitudes)'
         write(*,*)
         write(*,*)'-l[f_hz]: lowest frequency to consider (def=0)'
         write(*,*)'-h[f_hz]: highest frequency to consider (def=Nyqst)'
         write(*,*)'-f[fold]: give fold number instead of filename'
         write(*,*)'-n[nbin]: number of bins (histograms) in stats mode'
         write(*,*)'-F[f_hz]: mark frequency and harmonics on plot'
         write(*,*)
         write(*,*)'Options given as: spec fold1.spc -l10'
         write(*,*) 
         stop
      endif
      
      call readspec(sfile,fold,samp,dm,ac,tsamp,npf)
      

      if (fold.eq.0) then
         fold=1
         title='Raw Spectrum'
      else
         title=' '
      endif
      
      do i=1,npf
         fbin(i)=real(freq(tsamp,npf,fold,i))
         if (power) samp(i)=samp(i)*samp(i)
      enddo

      if (rbin.eq.0) then
         rbin=npf
         nadd=1
      else
         rbin=2**rbin
	 nadd=rbin
      endif
c      if (rbin.gt.npf) rbin=npf
c      nadd=npf/rbin
      if (nadd.gt.1) call rebin2(fbin,samp,npf,nadd)
      nh=1
      np=1
      if (ppsr.gt.0.0) then
         fpsr=1000.0/ppsr
         deltaf=fpsr*0.05
         nh=int(fbin(npf)/fpsr)
         if (nh.gt.50) nh=50
         np=2
      endif

      call pgbegin(0,pgdev,1,np)
      call pgscf(2)
      call pgsch(1.5)
      call pgvport(0.15,0.85,0.15,0.85)
      if (title.eq.' ') write(title,'(''Fold # '',i1)') fold

      do ih=1,nh
      call pgadvance
      fhar=real(ih)*fpsr
      
      if (stats) then
        call minmax(samp,npf,smin,smax)
        if (lo.lt.0.0) lo=smin
        if (hi.lt.0.0) hi=smax
      else
        if (lo.lt.0.0) lo=fbin(1)
        if (hi.lt.0.0) hi=fbin(npf)
        if (fpsr.gt.0.0) then
           lo=fhar-deltaf
           hi=fhar+deltaf
        endif
        j=0
        do i=1,npf
           if (fbin(i).ge.lo.and.fbin(i).le.hi) then
              j=j+1
              fbin(j)=fbin(i)
              samp2(j)=samp(i)
              samp3(j)=samp3(j)+samp2(j)
           endif
        enddo
        call minmax(samp2,j,smin,smax)
c	smax=5000.0
        snr=spcsnr(samp3,j)
      endif


      if (stats) then
         sum=0.0
         sumsq=0.0
        do i=1,npf
           sum=sum+samp(i)
           sumsq=sumsq+samp(i)*samp(i)
        enddo
        mean=sum/float(npf)
        sigma=sqrt(sumsq/float(npf)-mean*mean)
        write(*,*) "Mean:",mean," Sigma:",sigma
        nmax=histmax(npf,samp,lo,hi,nbins)
        call pgwindow(lo,hi,0.0,nmax) 
        call pgbox('bnst',0.0,0,'bcnst',0.0,0)
        call pghist(npf,samp,lo,hi,nbins,1)
        call pglabel('Amplitude','Number of occurences',title)
      else
        call pgwindow(lo,hi,smin-abs(smin*0.1),smax*1.1)
        call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
c        call pgbox('',0.0,0,'',0.0,0)
        call pgline(j,fbin,samp2)
        if (dump) then
           open(unit=55,file='dump',status='unknown')
           do i=1,j
              if (fbin(i).ge.lo.and.fbin(i).le.hi)
     &             write(55,*)fbin(i),samp2(i)
           enddo
           close(unit=55)
        endif
        if (power) then
           call pglabel('Frequency (Hz)','Power',title)
        else
           call pglabel('Frequency (Hz)','Amplitude',title)
        endif
        if (fmark.gt.0.0) then
           call pgsch(2.5)
           call pgsci(7)
           do i=1,256
              call pgpoint(1,real(i)*real(fmark),smax/2.0,31)
           enddo
           call pgsch(1.5)
           call pgsci(1)
        endif
        if (fpsr.gt.0.0) then
          call pgsch(2.5)
          call pgpoint(1,fhar,smax,31)
          call pgsch(1.5)
          call pgadvance
          call minmax(samp3,j,smin,smax)
          call pgwindow(lo,hi,smin-abs(smin*0.1),smax*1.1)
          call pgbox('bcst',0.0,0,'bcnst',0.0,0)
          call pgline(j,fbin,samp3)
        endif
      endif

      enddo
      call pgend
      end
      
      real function spcsnr(dat,n)
      implicit none
      integer n
      real dat(n)
      integer i,j,nb
      real rms,peak,sumsq
      peak=-1.0e32
      sumsq=0.0
      nb=n/10
      j=0
      do i=1,n
         if (i.le.nb.or.i.ge.n-nb) then
            sumsq=sumsq+dat(i)*dat(i)
            j=j+1
         else
            peak=max(peak,dat(i))
         endif
      enddo
      rms=sqrt(sumsq/real(j))
      spcsnr=0.0
      if (peak.gt.0.0) spcsnr=peak/rms
      end

      real function spcsnr0(dat,n)
      implicit none
      integer n
      real dat(n)
      integer i
      real rms,peak,sumsq
      peak=-1.0e32
      sumsq=0.0
      do i=1,n
         peak=max(peak,dat(i))
         if (i.le.5.or.i.ge.n-5) sumsq=sumsq+dat(i)*dat(i)
      enddo
      rms=sqrt(sumsq/10.0)
      spcsnr0=0.0
      if (peak.gt.0.0) spcsnr0=peak/rms
      end
      
         
      subroutine rebin2(fbin,samp,npf,nadd)
      implicit none
      real fbin(*), samp(*)
      integer npf,nadd
      integer i,j,k,l
      real sum

      j=0
      k=0
      l=0
      do i=1,npf
         j=j+1
         if (samp(i).ne.0.0) then
           sum=sum+samp(i)
           l=l+1
         endif
         if (j.eq.nadd) then
           k=k+1
           samp(k)=sum/real(l)
           sum=0.0
           fbin(k)=(fbin(i)+fbin(i-nadd+1))/2.0
           j=0
           l=0
         endif
      enddo
      npf=k
      end
         
