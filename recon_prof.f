c================================================================
      subroutine recon_prof(llog,npf,fold,pcand,candno,peak)
c================================================================
c
C performs profile reconstruction from the complex spectrum
C does this by doing a 64point fft for all of the harmonics
C currently using hermitian style
c
c================================================================
      
      implicit none
      include 'seek.inc' 
      include 'csamp.inc'
      real peak,profile(256),df0,amp
      real*8 pcand
      integer llog, slun
      integer npf,nfft,ifun,i,candno,fold,nharm,nharmmax
      character*12 sfilename
            
c      integer nharm, nharmmax
      
      nharm=2**(fold-1)



c 
c calculate the position of the fundamental
c

c      fbest=(1/0.45508644271)

      df0=2.0*real(npf)*tsamp/(pcand/1000.0)
c      df0=2*fbest

cc      nharmmax=min(15,npf/int(df0)) 
      nharmmax=min(16,npf/int(df0))  !!!!!!!!!!!!!!!!!!!!!!!

c      write(llog,*)'period in ms', (1000.0/fbest)
c      write(llog,*)'number of spectral bins', npf
c      write(llog,*)'tsamp', tsamp
c      write(llog,*)'candidate period s', (pcand/1000.0)
c      write(llog,*)'fundamental df0', df0
c      write(llog,*)df0
c      write(llog,*)'harmonic fold', fold
c      write(llog,*)'nharm', nharm
c      write(llog,*)'nharmmax', nharmmax
C fill profile with the data
c      nfft=2*nharm
C for interpolation this should be 64 -PAH

c      nfft=32
      nfft=64  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C zero out the rest
      do i=1,2*nfft
         profile(i)=0.0
      enddo
c      write(llog,*)series(100000)
C put in the harmonics leaving dc 0 - up to 15
C if only doing 32 point transforms then cannot have the 16th harmonic in
C the nyquist position only do up to 15
c
c     include fundamental in the array
c
c      open(78,file='cdat_prof.tmp',status='unknown',access='append')

      do i=1,min(nharm,nharmmax)

c      do i=1,4
       ifun=nint(df0*i)*2+1 
c        ifun=nint(df0)*i
c         write(llog,*)'ifun', ifun
c        write(*,*) 'Harmonic ',i,' Power ',
c     &  sqrt(cdat(ifun)**2+cdat(ifun+1)**2)
c        write(llog,*)series(ifun),series(ifun+1)
        profile(i*2+1)=series(ifun)
        profile(i*2+2)=series(ifun+1)
        profile((nfft-i)*2+1)=series(ifun)
        profile((nfft-i)*2+2)=-series(ifun+1)
      enddo

c================================
c      Test here

c      do i=1,min(nharm,nharmmax)
c        ifun=nint(df0*i)*2+1
c        write(*,*) 'Harmonic ',i,' Power ',
c     &  sqrt(cdat(ifun)**2+cdat(ifun+1)**2)
c        write(llog,*)series(ifun),series(ifun+1)
c        profile(i*2+1)=samp(ifun)
c        profile(i*2+2)=samp(ifun+1)
c        profile((nfft-i)*2+1)=samp(ifun)
c        profile((nfft-i)*2+2)=-samp(ifun+1)
c      enddo  

c================================


c.      do i = 1,16
c.        write(78,*) i*2+1,profile(i*2+1),i*2+2,profile(i*2+2),
c.     &  sqrt(profile(i*2+1)**2+profile(i*2+2)**2)
c.        write(*,*) i*2+1,profile(i*2+1),i*2+2,profile(i*2+2),
c.     &  sqrt(profile(i*2+1)**2+profile(i*2+2)**2)
c.      end do


C reconstruct profile
         call sglfft(profile(1),profile(2),nfft,nfft,nfft,-2)
cc        call four1(profile,nfft,-1)

C find peak 
c.         write(llog,*)' ----------------------------------------------'
c.         write(llog,*)' Period   Nharm'
c.         write(llog,*)pcand,nharm
c.         write(llog,*)' Real     Imag      Amp'
         peak=0
         do i=1,nfft*2,2
c            amp=sqrt(profile(i)*profile(i)+profile(i+1)*profile(i+1))
c            write(llog,*)i,amp
c.           write(*,*)i/2,amp
             peak=max(peak,profile(i))
         enddo
c        write(llog,*)'***', amp,peak
c        write(llog,*)
         
c         call glun(slun)
c         write(*,*)"Dumping reconstructed profile"
c         write(sfilename,"(i2.2,a,i5.5)")fold,"prof",candno
c         write(*,*)"Dumping to  ",sfilename
c         open(unit=slun,file=sfilename,status='unknown')
cc        do i=1,nfft*2,2
cc         do i=1,nfft,2
c         do i=1,nfft*2
cc           amp=sqrt(profile(i)*profile(i)+profile(i+1)
cc     &           *profile(i+1))
cc            write(slun,*)i,amp 
c            write(slun,*)profile(i) 
c         enddo   
c         close(slun)


      end

