C **************************************************************
      SUBROUTINE SUMHRM (SP,NF,NF1)
C **************************************************************
      include 'seek.inc'
      INTEGER NF, NF1
      REAL SP(NF)
      INTEGER ifold
      REAL fval
      INTEGER NFS
      REAL X, XDIV, XDEL
      INTEGER IDX, LSTIDX


c SPH(npts)
      do ifold=1,nfolds
         fval = foldvals(ifold)
         if (fval.eq.1) then
            do n=1,nf
               power(1,n) = sp(n)
            enddo
         else
            NFS = MAX(1,MIN(foldvals(ifold)*NF1-foldvals(ifold)/2,NF))
            XDIV = 1./fval
            do 42 N=NFS,NF,1
               power(ifold,N)=0.
               LSTIDX = -1
               do 44 I=1,foldvals(ifold),1
                  X = N * I * XDIV + 0.5
                  IDX = X
c                  IF (IDX.EQ.1) THEN
c                     write(*,*) 'N=',N,', I=',I,', fval=',
c     +                    foldvals(ifold),', IDX=',IDX
c                  ENDIF
                  IF (IDX.GT.1.AND.IDX.NE.LSTIDX) THEN
                     power(ifold,N) = power(ifold,N) + SP(IDX)
                  ENDIF
                  LSTIDX = IDX
 44            continue
 42         continue   
         endif
      enddo
c     
      RETURN
C
C END OF SUBROUTINE SUMHRM
C
      END
C **************************************************************
      SUBROUTINE OLDSUMHRM (SP,NF,NF1)
C **************************************************************
      include 'seek.inc'
      INTEGER NF, NF1
      REAL SP(NF), SPH(npts)

      NF2 = MAX(1,MIN(2*NF1-1,NF))
      NF4 = MAX(1,MIN(4*NF1-2,NF))
      NF8 = MAX(1,MIN(8*NF1-4,NF))
      NF16= MAX(1,MIN(16*NF1-8,NF))
      do n=1,nf
         power(1,n)=sp(n)
      enddo
      DO N=1,NF
         SPH(N)=0.0
      ENDDO

      K = (NF2+1)/2
      DO 20 N=NF2,NF-1,2
         SPH(N)=SP(N)+SP(K)
         SPH(1+N)=SP(1+N)+SP(K)
         K=K+1
   20 CONTINUE
      do n=1,nf
         power(2,n)=sph(n)
      enddo
      
      KA = (NF4+2)/4
      KB = (3*NF4+2)/4
      DO 40 N=NF4,NF-3,4
         SPH(N)=SPH(N)+SP(KA)+SP(KB)
         SPH(1+N)=SPH(1+N)+SP(KA)+SP(KB)
         SPH(2+N)=SPH(2+N)+SP(KA)+SP(1+KB)
         SPH(3+N)=SPH(3+N)+SP(KA)+SP(2+KB)
         KA = KA + 1
         KB = KB + 3
   40 CONTINUE
      do n=1,nf
         power(3,n)=sph(n)
      enddo

      JA = (NF8+4)/8
      JB = (3*NF8+4)/8
      JC = (5*NF8+4)/8
      JD = (7*NF8+4)/8
      DO 80 N=NF8,NF-7,8
         SPH(  N)=SPH(  N)+SP(JA)+SP(  JB)+SP(  JC)+SP(  JD)
         SPH(1+N)=SPH(1+N)+SP(JA)+SP(  JB)+SP(  JC)+SP(  JD)
         SPH(2+N)=SPH(2+N)+SP(JA)+SP(  JB)+SP(1+JC)+SP(1+JD)
         SPH(3+N)=SPH(3+N)+SP(JA)+SP(1+JB)+SP(1+JC)+SP(2+JD)
         SPH(4+N)=SPH(4+N)+SP(JA)+SP(1+JB)+SP(2+JC)+SP(3+JD)
         SPH(5+N)=SPH(5+N)+SP(JA)+SP(1+JB)+SP(3+JC)+SP(4+JD)
         SPH(6+N)=SPH(6+N)+SP(JA)+SP(2+JB)+SP(3+JC)+SP(5+JD)
         SPH(7+N)=SPH(7+N)+SP(JA)+SP(2+JB)+SP(4+JC)+SP(6+JD)
         JA = JA + 1
         JB = JB + 3
         JC = JC + 5
         JD = JD + 7
   80 CONTINUE
      do n=1,nf
         power(4,n)=sph(n)
      enddo

      LA = (   NF16+8)/16
      LB = (3*NF16+8)/16
      LC = (5*NF16+8)/16
      LD = (7*NF16+8)/16
      LE = (9*NF16+8)/16
      LF = (11*NF16+8)/16
      LG = (13*NF16+8)/16
      LH = (15*NF16+8)/16
      DO 160 N=NF16,NF-15,16
         SPH(  N)=SPH(  N)+SP(LA)+SP(  LB)+SP(  LC)+SP(  LD)+SP(  LE)+
     +    SP(   LF)+SP(   LG)+SP(   LH)
         SPH(1+N)=SPH(1+N)+SP(LA)+SP(  LB)+SP(  LC)+SP(  LD)+SP(  LE)+
     +    SP(   LF)+SP(   LG)+SP(   LH)
         SPH(2+N)=SPH(2+N)+SP(LA)+SP(  LB)+SP(  LC)+SP(  LD)+SP(1+LE)+
     +    SP( 1+LF)+SP( 1+LG)+SP( 1+LH)
         SPH(3+N)=SPH(3+N)+SP(LA)+SP(  LB)+SP(  LC)+SP(1+LD)+SP(1+LE)+
     +    SP( 2+LF)+SP( 2+LG)+SP( 2+LH)
         SPH(4+N)=SPH(4+N)+SP(LA)+SP(  LB)+SP(1+LC)+SP(1+LD)+SP(2+LE)+
     +    SP( 2+LF)+SP( 3+LG)+SP( 3+LH)
         SPH(5+N)=SPH(5+N)+SP(LA)+SP(  LB)+SP(1+LC)+SP(2+LD)+SP(2+LE)+
     +    SP( 3+LF)+SP( 4+LG)+SP( 4+LH)
         SPH(6+N)=SPH(6+N)+SP(LA)+SP(1+LB)+SP(1+LC)+SP(2+LD)+SP(3+LE)+
     +    SP( 4+LF)+SP( 4+LG)+SP( 5+LH)
         SPH(7+N)=SPH(7+N)+SP(LA)+SP(1+LB)+SP(2+LC)+SP(3+LD)+SP(3+LE)+
     +    SP( 4+LF)+SP( 5+LG)+SP( 6+LH)
         SPH(8+N)=SPH(8+N)+SP(LA)+SP(1+LB)+SP(2+LC)+SP(3+LD)+SP(4+LE)+
     +    SP( 5+LF)+SP( 6+LG)+SP( 7+LH)
         SPH(9+N)=SPH(9+N)+SP(LA)+SP(1+LB)+SP(2+LC)+SP(3+LD)+SP(5+LE)+
     +    SP( 6+LF)+SP( 7+LG)+SP( 8+LH)
         SPH(10+N)=SPH(10+N)+SP(LA)+SP(1+LB)+SP(3+LC)+SP(4+LD)+SP(5+LE)+
     +    SP( 6+LF)+SP( 8+LG)+SP( 9+LH)
         SPH(11+N)=SPH(11+N)+SP(LA)+SP(2+LB)+SP(3+LC)+SP(4+LD)+SP(6+LE)+
     +    SP( 7+LF)+SP( 8+LG)+SP(10+LH)
         SPH(12+N)=SPH(12+N)+SP(LA)+SP(2+LB)+SP(3+LC)+SP(5+LD)+SP(6+LE)+
     +    SP( 8+LF)+SP( 9+LG)+SP(11+LH)
         SPH(13+N)=SPH(13+N)+SP(LA)+SP(2+LB)+SP(4+LC)+SP(5+LD)+SP(7+LE)+
     +    SP( 8+LF)+SP(10+LG)+SP(12+LH)
         SPH(14+N)=SPH(14+N)+SP(LA)+SP(2+LB)+SP(4+LC)+SP(6+LD)+SP(7+LE)+
     +    SP( 9+LF)+SP(11+LG)+SP(13+LH)
         SPH(15+N)=SPH(15+N)+SP(LA)+SP(2+LB)+SP(4+LC)+SP(6+LD)+SP(8+LE)+
     +    SP(10+LF)+SP(12+LG)+SP(14+LH)
         LA = LA + 1
         LB = LB + 3
         LC = LC + 5
         LD = LD + 7
         LE = LE + 9
         LF = LF + 11
         LG = LG + 13
         LH = LH + 15
  160 CONTINUE
      do n=1,nf
         power(5,n)=sph(n)
      enddo

      RETURN
C
C END OF SUBROUTINE SUMHRM
C
      END
