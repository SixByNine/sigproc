c==============================================================================
      subroutine getrrms(data,ndat,rmea,nrun,rrms)
c==============================================================================
c
c     Calculates the running rms of the data in the array data(1:ndat) after
c     subtracting off its running mean contained in the array rmea(1:nrun).
c      
c     This calculation AVOIDS rfi/psr spikes which bias the overall rms.
c      
c     data    - r4  - data array
c     ndat    - i4  - number of data points
c     rmea    - r4  - array to store running mean
c     nrun    - i4  - number of points to average
c      
c     rrms    - r4  - array to store running rms      
c
c     Created: November 1997 (dunc@mpifr-bonn.mpg.de)
c
c==============================================================================
c     
      implicit none
      real data(*),rmea(*),rrms(*)
      integer ndat,nrun
c
c     Local variables...
c      
      logical zap,spike
      integer i,j,k,l
      real sumsq,rms,mean,lastgood
c
c     Initialise...
c      
      j=0
      k=0
      l=0
      sumsq=0.0
      rms=0.0
      mean=rmea(1)
      lastgood=1.0
c
c     Main loop
c      
      do i=1,ndat
c
c       Check whether data point has been zapped i.e. it is exactly 0.0
c         
        zap=data(i).eq.0.0
c
c       Watch out for any remaining spikes >3 sigma (e.g. pulsars!)
c        
        spike=mean.ne.0.0.and.data(i).gt.mean+rms*3.0
c
c       Update running sum of squares if this not a zapped point or spike
c        
        if ((.not.zap).and.(.not.spike)) then
          l=l+1
          sumsq=sumsq+(data(i)-mean)*(data(i)-mean)
        endif
        j=j+1
c
c       Calculate rms if required sum reached...
c        
        if (j.eq.nrun) then
           k=k+1
           mean=rmea(k)
           if (l.gt.0) then
             rrms(k)=sqrt(sumsq/real(l))
             lastgood=rrms(k)
           else
             rrms(k)=lastgood
           endif
           rms=rrms(k)
           sumsq=0.0
           j=0
           l=0
        endif
      enddo
      end
c
c     That's all folks!
c      
c==============================================================================
