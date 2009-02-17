c=============================================================================
      real*8 function freq(tsamp,npf,fold,k)
c=============================================================================
c
c     Returns the fluctuation frequency (Hz) of bin k in the Fourier
c     spectrum having npf points. Tsamp is the sampling interval
c     of the corresponding time domain data (seconds), whilst fold
c     refers to 1 plus the number of harmonic sums that have produced
c     the present spectrum. e.g fold=1 -> refers to the raw spectrum
c     fold=5 refers to 4 harmonic sums (16 harmonics).
c
      real*8 tsamp
      integer npf,fold,k
      freq=real(k)/(2.0*tsamp*npf*2.0**(fold-1))
      end
c=============================================================================
      integer function fbin(tsamp,npf,fold,freq)
c=============================================================================
c
c     Just the mathematical inverse of freq. Returns the bin number for
c     a given frequency.
c      
      real*8 tsamp
      integer npf,fold
      real freq
      fbin=real(freq)*(2.0*real(tsamp)*npf*2.0**(fold-1))
      end
c=============================================================================
      real*8 function freqff(tsamp,npf,fold,k)
c=============================================================================
c
c     Returns the fluctuation frequency (Hz) of bin k in the Fourier
c     spectrum having npf points. Tsamp is the sampling interval
c     of the corresponding time domain data (seconds), whilst fold
c     refers to the number of harmonics summed, eg. 1,2,4,6,...
c
      real*8 tsamp
      integer npf,fold,k
      freqff=real(k)/(2.0*tsamp*npf*fold)
      end
c=============================================================================
