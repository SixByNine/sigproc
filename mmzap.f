cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zap_mmbrd(dat,dat2,nf,nf1,tsmp_data,sprms,mean,ifold)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     zap interference from Methanol Multibeam survey data.
c     adapted from zap_birdies.f to knock out just harmonics of
c     the sampler cards and the mains hum. In the language of
c     zap_birdies, these are treated as two "special sequences"
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      real perzapped
      integer zapped

      integer nf, nf1,ifold
      real dat(nf), sprms, mean
      real dat2(2*nf)
      double precision tsmp_data

      integer nb, nh, n1, n2, n1store, n2store, nt, tempn1
      real snhrm, sig4, fnyq, fnyq2, fnyq4
      integer j
      integer maxspseq, nspecial

      parameter (maxspseq=2500,nspecial=2)
      integer ii, jj, k,  nspseqf(nspecial)
      real freq, spseqf(nspecial,maxspseq), spseqw(nspecial,maxspseq)
      real binwidth

      fnyq=500./tsmp_data
      fnyq2=2.*fnyq
      fnyq4=4.*fnyq
      binwidth = fnyq/nf
c      
c     Special sequence 1: 'Sample counter': Harmonics of 8192 ms
c
      k=1
      do ii=1,15
        do jj=1,10
           freq=abs(2d0**(ii-1)/8.192d0)*(1.+2*(jj-1))
           if(freq.le.fnyq) then 
              spseqf(1,k)=freq
              spseqw(1,k)=binwidth
              k = k + 1
           endif
        enddo
      enddo
      nspseqf(1) = k-1
c
c     Special sequence 2: f = 50.0*i
c
      k = 1
      do ii=1,fnyq/50. - 1
         freq=ii*50.
         spseqf(2,k)=freq
         spseqw(2,k)=0.2
         k = k + 1
      enddo
      nspseqf(2) = k-1


      zapped=0
      sig4 = mean + 4.0 * sprms

      do nb=1,nspecial

         do nh = 1, nspseqf(nb)
            n1=nint(0.001*(spseqf(nb,nh) - spseqw(nb,nh))
     &           * tsmp_data*nf*2)
            n2=nint(0.001*(spseqf(nb,nh) + spseqw(nb,nh))
     &           * tsmp_data*nf*2)
         
            n1store=n1
            n2store=n2
         
cDCS remove multiples of 2*Nyquist frequency from the filter range and,
cDCS if the frequency is still greater than the Nyquist frequency, reflect
cDCS about the Nyquist frequency
            n1=max(mod(n1,2*nf),1) ! Avoid n1=0 case 
            n2=max(mod(n2,2*nf),1) ! Avoid n2=0 case 
            if (n1.gt.nf)n1=2*nf-n1
            if (n2.gt.nf)n2=2*nf-n2
c     DCS make sure that n2 is greater than n1 as later loops rely on this
            if (n1.gt.n2) then
               tempn1=n1
               n1=n2
               n2=tempn1
            endif
cDCS check to see if n1 and n2 are derived from a different number of
cDCS reflections.  If so then make the filter as big as possible by setting
cDCS n1 or n2 to the nearest extreme frequency.
            if ((n1store/nf).ne.(n2store/nf))then
c
c              check which frequency extreme is straddled
c              note that this is a check for an odd number of reflections
c
               if (n2store/nf.ne.2*(n2store/(2*nf)))then
                  n2=nf
               else
                  n1=1
               endif
            endif
c     
c  work out the s/n in each harmonic and only 
c  zap it if above threshhold
c
            snhrm=0.0
            do j=n1,n2
               snhrm=max(snhrm,dat(j))
            end do
            if (snhrm.gt.sig4) then
               do j=n1,n2
                  if ((dat(j).gt.sig4).or.(nh.le.2)) then
                     dat(j)=0.0
                     dat2(2*j) = 0.0           !!!!!
                     dat2(2*j+1) = 0.0           !!!!!
                     zapped = zapped+1
                  endif
               enddo
            endif
            nt=nt+n2-n1
         end do
      enddo
c
c     Delete Nyq freq. 
c      
      dat(nf-1)=0.
      dat(nf)=0.
      perzapped = (float(zapped)/float(nf))*100.0
      write(*,'('' Zap:'',i4,f8.5,f6.2)')ifold,tsmp_data,perzapped

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

