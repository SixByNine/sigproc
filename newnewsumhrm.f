C **************************************************************
      SUBROUTINE NEWNEWSUMHRM (SP,NF,NF1,P1,P2,P4,P8,P16)
C **************************************************************
      IMPLICIT NONE
      INTEGER NF, NF1
      REAL SP(NF)
      INTEGER ifold
      REAL fval
      INTEGER NFS
      REAL X, XDIV
      INTEGER IDX, LSTIDX
      INTEGER I, N, NFOLDS
      INTEGER FOLDVALS(5)
      REAL P1(*),P2(*),P4(*),P8(*),P16(*)

      NFOLDS=5
      FOLDVALS(1)=1
      FOLDVALS(2)=2
      FOLDVALS(3)=4
      FOLDVALS(4)=8
      FOLDVALS(5)=16

c SPH(npts)
      do ifold=1,nfolds
         fval = foldvals(ifold)
         if (fval.eq.1) then
            do n=1,nf
               p1(n) = sp(n)
            enddo
         endif
         if (fval.eq.2) then
            NFS = MAX(1,MIN(foldvals(ifold)*NF1-foldvals(ifold)/2,NF))
            XDIV = 1./fval
            do N=NFS,NF
               p2(N)=0.
               LSTIDX = -1
               do I=1,foldvals(ifold)
                  X = N * I * XDIV + 0.5
                  IDX = X
                  IF (IDX.GT.1.AND.IDX.NE.LSTIDX) THEN
                     p2(N) = p2(N) + SP(IDX)
                  ENDIF
                  LSTIDX = IDX
               enddo
            enddo
         endif
         if (fval.eq.4) then
            NFS = MAX(1,MIN(foldvals(ifold)*NF1-foldvals(ifold)/2,NF))
            XDIV = 1./fval
            do N=NFS,NF
               p4(N)=0.
               LSTIDX = -1
               do I=1,foldvals(ifold)
                  X = N * I * XDIV + 0.5
                  IDX = X
                  IF (IDX.GT.1.AND.IDX.NE.LSTIDX) THEN
                     p4(N) = p4(N) + SP(IDX)
                  ENDIF
                  LSTIDX = IDX
               enddo
            enddo
         endif
         if (fval.eq.8) then
            NFS = MAX(1,MIN(foldvals(ifold)*NF1-foldvals(ifold)/2,NF))
            XDIV = 1./fval
            do N=NFS,NF
               p8(N)=0.
               LSTIDX = -1
               do I=1,foldvals(ifold)
                  X = N * I * XDIV + 0.5
                  IDX = X
                  IF (IDX.GT.1.AND.IDX.NE.LSTIDX) THEN
                     p8(N) = p8(N) + SP(IDX)
                  ENDIF
                  LSTIDX = IDX
               enddo
            enddo
         endif
         if (fval.eq.16) then
            NFS = MAX(1,MIN(foldvals(ifold)*NF1-foldvals(ifold)/2,NF))
            XDIV = 1./fval
            do N=NFS,NF
               p16(N)=0.
               LSTIDX = -1
               do I=1,foldvals(ifold)
                  X = N * I * XDIV + 0.5
                  IDX = X
                  IF (IDX.GT.1.AND.IDX.NE.LSTIDX) THEN
                     p16(N) = p16(N) + SP(IDX)
                  ENDIF
                  LSTIDX = IDX
               enddo
            enddo
         endif
      enddo
c     
      RETURN
C
C END OF SUBROUTINE NEWNEWSUMHRM
C
      END
