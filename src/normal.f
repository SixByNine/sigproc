c---------------------------------------------
      real function normal(seed, mu, sigma)
c---------------------------------------------
c
c     Gaussian (normal) distribution function
c
      implicit none
      real mu, sigma, rnd1, rnd2
      integer*4 seed
      real donran
      rnd1 = donran(seed)
      rnd2 = donran(seed)
      if (rnd1 .eq. 0.0) rnd1 = donran(seed)
      normal = ((sigma * ((- (2.0 * log(rnd1))) ** 0.5)) 
     &               * cos((2.0 *  3.1415926) * rnd2)) + mu
      end
c==============================================================================
      function donran(idum)
c==============================================================================
c
c     Donald Knuth's portable random number generator, nabbed from p199
c     of numerical recipes. Gives same random numbers on both Sun's &
c     Hp's given the same initial seed. Pass idum down as negative to 
c     reshuffle the random numbers.
c
c     Last Change 98/02/25 DRL @ MPIfR -> added SAVE
c
c      function ran3(idum)
c         implicit real*4(m)
c         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      dimension ma(55)
      data iff /0/
      save
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
      donran=mj*fac

      end

