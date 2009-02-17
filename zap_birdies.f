C @(#)zap_birdies.f	3.2 18 Feb 1994
c------------------------------------------------------------
      subroutine zap_pmbrd(dat,dat2,nf,nf1,tsmp_data,sprms,mean,ifold)
c------------------------------------------------------------
c
c     Zaps the birdies from the spectrum dat(nf) if above
c     four sigma (sig4) threshold DRL 92/06/17 (hijacked from
c     kill_birdies.f)
c
c JFB 971010  Mods cMJFB to fix quantisation problems for low freq filters
c             Pass around freqs rather bins to stop quantisation spreading
c             to higher harmonics
cDCS 971013   routine now zaps reflected harmonics for any number
c             of reflections
c VMK 980618  added special case handling of the 135 Hz birdie and
c             keeps track of the number of bins zapped.
cDCS increased range of intermodulation zapping
c VMK 980723  Added internal zapping of 4.096 s birdie and harmonics,
c             and changed internal zapping of 135 Hz intermodulation
c             freqs so does not assume 135 Hz is in std filter list.
c IHS 981028  Commented out writes to fort.61 and fort.135 - hope that's OK
c AGL 000120  Restructure to allow multiple special sequences
c GBH 011019  Made changes to the sequences:
c              - added wide sidelobes to 135 Hz sequence
c              - slight changes to widths for most sequences.
c              - widths now do not increase with harmonic number
c              - added 0.14704*ii sequences (removed from bhbird)

      implicit none

      include 'birdies.inc'
      include 'find_summary.inc'

      integer nf, nf1,ifold
      real dat(nf), sprms, mean
      real dat2(2*nf)
      double precision tsmp_data

      integer nb, nh, n1, n2, n1store, n2store, nt, tempn1
      real snhrm, sig4, fnyq, fnyq2, fnyq4
      integer j
      integer maxspseq, nspecial
c      parameter (maxspseq=1166, nspecial=5)
cGBH Changed maxspseq=2500 as sequence 2 now has 2392 filters
      parameter (maxspseq=2500,nspecial=5)
      integer ii, jj, k,  maxnh, nspseqf(nspecial)
      real freq, spseqf(nspecial,maxspseq), spseqw(nspecial,maxspseq)
      real binwidth,f1,f2

cAGL  the special sequences are those which are essentially endemic to
cAGL  the whole of the survey database.  Many contain harmonics which 
cAGL  may be individually small but add up to something significant when
cAGL  folded.  Others are not straight harmonic sequences.
cDCS  maxspseq = max no of features in any special sequence = 41+75*15 = 1166

      fnyq=500./tsmp_data
      fnyq2=2.*fnyq
      fnyq4=4.*fnyq

cGBH: Frequencies for the 135 and 99.9Hz sequences
      f1=134.99776
      f2=99.99592 

cGBH: Width of a bin

      binwidth = fnyq/nf
c     Special sequence 1: 'Sample counter': Harmonics of 8192 ms
      k=1
      do ii=1,15
        do jj=1,10
           freq=abs(2d0**(ii-1)/8.192d0)*(1.+2*(jj-1))
           if(freq.le.fnyq) then 
              spseqf(1,k)=freq
c              spseqw(1,k)=0.0005
              spseqw(1,k)=binwidth

              k = k + 1
           endif
        enddo
      enddo
      nspseqf(1) = k-1
c Special sequence 2: f = 135*i + 100*j
      k=1
      do ii=0,160
         if(ii.eq.0) then
            do jj=0,120
c               freq = abs(ii*134.99779 + jj*99.99570)
               freq = abs(jj*f2)
               freq = mod(freq,fnyq2)

               if (freq.gt.fnyq) freq=fnyq2-freq
c     GBH: Don't zap below 2Hz
               if (freq.gt.2.0) then
                  spseqf(2,k)=freq
c     spseqw(2,k)=0.015
                  spseqw(2,k)=0.05
                  k = k + 1
               endif
            enddo
         else
            if (ii.le.50) then
               do jj=-20,20
                  freq = abs(ii*f1 + jj*f2)
                  freq = mod(freq,fnyq2)
                  if (freq.gt.fnyq) freq=fnyq2-freq
                  if (freq.gt.1.00) then
                     spseqf(2,k)=freq
                     spseqw(2,k)=0.020
                     k = k + 1                     
                  endif
               enddo
            else
               freq = abs(ii*f1)
               freq = mod(freq,fnyq2)
               if (freq.gt.fnyq) freq = fnyq2-freq
               spseqf(2,k)=freq
               spseqw(2,k)=0.020
               k = k + 1                     
               
            endif
         endif
      enddo
C 135Hz and 0.0308 sidelobes (covered with one 0.2 Hz filter)

      do ii=1,fnyq4/f1
         freq = abs(ii*f1)
         freq = mod(freq,fnyq2)

         if (freq.gt.fnyq) freq=fnyq2-freq
         spseqf(2,k)=freq
         spseqw(2,k)=0.2
         k = k + 1      
      enddo

c 135 Hz, wide sidelobes of f = 135ii
      do ii=1,fnyq2/f1
         freq = abs(ii*f1)
         freq = mod(freq,fnyq2)
         if (freq.gt.fnyq) freq=fnyq2-freq

         spseqf(2,k)=freq-3.90
         spseqw(2,k)=0.15
         k = k + 1      

         spseqf(2,k)=freq+3.90
         spseqw(2,k)=0.15
         k = k + 1      
      enddo
      nspseqf(2) = k-1
            
c     Special sequence 3: f = 50.0*i
      k = 1
      do ii=1,fnyq/50. - 1
         freq=ii*50.
         spseqf(3,k)=freq
         spseqw(3,k)=0.1
         k = k + 1
      enddo
      nspseqf(3) = k-1

c     Special sequence 4: f = 200 + 1.0*i
      k=1
      do ii=-20,20
         freq=200.+ii*1.0
         spseqf(4,k)=freq
C         spseqw(4,k)=0.1
         spseqw(4,k)=0.13
         k = k + 1
      enddo
      nspseqf(4) = k-1
      
cGBH     Special sequence 5: f=0.14704*ii
      k=1
      do ii=1,100
         freq = 0.14704*ii
         spseqf(4,k)=freq
         spseqw(4,k)=max(0.0012,binwidth)
         k = k + 1
      enddo
      nspseqf(5) = k-1

      zapped=0

      sig4 = mean + 4.0 * sprms
c      write(61,*)nf,nf1,tsmp_data,sprms
      
      do nb=1,nbrd + nspecial

cDCS if a standard filter is not valid for this beam then nba and nbb
cDCS will be set to zero
         if (nb.le.nbrd) then
            if(bffa(nb).eq.0..and.bffb(nb).eq.0.)go to 55
            if(ifold.eq.1)snbrd(nb)=0.
            maxnh = nint(abs(float(nhbrd(nb))))
         else
            maxnh = nspseqf(nb-nbrd)
         endif

cAGL Now the zapping loop.
         do nh = 1, maxnh
            if (nb.gt.nbrd)  then
               n1=nint(0.001*(spseqf(nb-nbrd,nh) - spseqw(nb-nbrd,nh))
     &              * tsmp_data*nf*2)
               n2=nint(0.001*(spseqf(nb-nbrd,nh) + spseqw(nb-nbrd,nh))
     &              * tsmp_data*nf*2)
            else
cMJFB the 0.5 is the replacements for the delta that was in find_find_hum
               n1=nint(nh*0.001*bffa(nb)*tsmp_data*nf*2 - 0.5)
               n2=nint(nh*0.001*bffb(nb)*tsmp_data*nf*2 + 0.5)
            endif
            
            n1store=n1
            n2store=n2


cVMK  Only do zapping of aliases of birdies past the Nyquist for the two
cVMK  special cases to avoid zapping of huge parts of the spectrum at low freq. 
            if ((n1.gt.nf).and.(nb.le.nbrd)) go to 55

c            if(nb.le.nbrd)then
c            write(61,610)nb,nh,nh*(bffb(nb)+bffa(nb))/2.,
c     +         nh*(bffb(nb)-bffa(nb)),n1,n2,n2-n1+1,1000.*(n2-n1+1)/nf
c            else
c            write(61,610)nb,nh,spseqf(nb-nbrd,nh),spseqw(nb-nbrd,nh)
c     +         ,n1,n2,n2-n1+1,1000.*(n2-n1+1)/nf
c            endif
c 610        format(2i6,f10.3,f10.5,3i10,f8.3)

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
              ! check which frequency extreme is straddled
              ! note that this is a check for an odd number of reflections
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

cvmk Make sure not to count intermodulation frequencies here.
cIHS I think the point here is that we always want to zap 2^n harmonics
cIHS since they will harmonically sum to trouble, and to let every other 
cIHS sequence be zapped only if it needs to be.

cIHS Make 2^n the first sequence so if the 135 sequence picks up 2^n things,
cIHS they won't be counted as zapped multiple times

            if ((ifold.eq.1).and.(nb.le.nbrd))
     +           snbrd(nb)=max(snbrd(nb),snhrm)

            if ((snhrm.gt.sig4).or.(nb.eq.nbrd+1)) then
               do j=n1,n2
                  if ( (dat(j).gt.sig4).or.(nh.le.2).or.
     +                 (nb.eq.nbrd+1))then
                     dat(j)=0.0 
                     dat2(2*j) = 0.0           !!!!!
                     dat2(2*j+1) = 0.0           !!!!!
                     zapped = zapped+1
C     write(61,*) '    zapping j=',j
                  endif
               enddo
            endif
            nt=nt+n2-n1
         end do
 55      continue
      end do
cvmk   endif

c     Delete Nyq freq. 
      
      dat(nf-1)=0.
      dat(nf)=0.

c     Turn birdie amps into snrs
      if(ifold.eq.1) then
         do nb=1,nbrd
            snbrd(nb)=snbrd(nb)/sprms
         enddo
      endif

      perzapped = (float(zapped)/float(nf))*100.0
c      write(*,'('' Zap:'',i4,f8.5,f6.2)')ifold,tsmp_data,perzapped

      end


