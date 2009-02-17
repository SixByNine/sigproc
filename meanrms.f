c ******************************************************
      subroutine meanrms ( d, n, mean, rms, nrej )
c ******************************************************
c
c Returns the mean and rms of the real data array d(n)
c with iterative spike rejection.
c AGL Feb07
c
      implicit none
      integer n,nrej,nrejold,i
      real*4 d(*),s,ss,sum,sumsq,mean,rms
    

c      WRITE(60,*)d(1)
c      WRITE(*,*) n

c
      if ( n.le.1 ) then
          rms = 0.
          mean=0.
          nrej=0
          return
      endif
      nrej=0
      nrejold=0
      sum=0.
      sumsq=0.
      do i=1,n
         s=d(i)
         sum=sum+s
         sumsq=sumsq+s*s
      enddo

c Store sums
      s=sum
      ss=sumsq
 10   rms=sqrt(max(0.0,(sumsq-sum*sum/(n-nrej))/(n-nrej-1.0)))
      mean=sum/real(n-nrej)
c      write(*,*)n,nrej,sum,sumsq,mean,rms
      sum=s
      sumsq=ss

c Remove points from the sums which are .gt. 3*rms from mean 
      nrej=0
      do i=1,n
         if(abs(d(i)-mean) .gt. 3.*rms) then
            sum=sum-d(i)
            sumsq=sumsq-d(i)*d(i)
            nrej=nrej+1
         endif
      enddo
c Exit if no more points rejected.
      if(nrej.eq.nrejold) return
      nrejold=nrej
      goto 10
c
      end

C
C ******************************************************************
      REAL FUNCTION GRAN ( A, S, INIT )
C ******************************************************************
C
C GENERATES GAUSSIAN RANDOM NOISE MEAN A, RMS S
C INIT SHOULD BE A LARGE ODD INTEGER
C
C
C THIS ROUTINE IS INSTALLATION DEPENDENT
C
C VAX-11 FORTRAN VERSION
C
      B = 0
      DO 10 J=1,12
         B = B+RAN(INIT)
   10 CONTINUE
      GRAN = (B-6.0)*S+A
      RETURN
C
C END OF REAL FUNCTION GRAN
C
      END
c==============================================================================
      function ran(idum)
c==============================================================================
c
c     Donald Knuth's portable random number generator, nabbed from p199
c     of numerical recipes. Gives same random numbers on both Sun's &
c     Hp's given the same initial seed. Pass idum down as negative to 
c     reshuffle the random numbers.
c
c     Last Change 93/06/13 DRL @ JB.
c
c      function ran3(idum)
c         implicit real*4(m)
c         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran=mj*fac

      end




