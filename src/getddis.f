c==============================================================================
      subroutine getddis(llog,filename,nchan,dmidx)
c==============================================================================
      implicit none
      include 'time.inc'
      integer dmidx, lun,nsrec,nchan,i,j,k,l,llog
      parameter(nsrec=4096)
      byte ddat(nsrec) 
      character filename*(*), uprow*3
c==============================================================================
      uprow = char(27)//'[A'
      write(llog,*) 'Reading DMIDX:',dmidx,'...'
      call glun(lun)
      open(lun,file=filename,status='old',access='direct',
     &     form='unformatted',recl=nsrec)
      i=dmidx
      j=0
      k=0
      l=0
      write(llog,*)
      do while(k.lt.npts)
        read(lun,rec=i,err=1) (ddat(j),j=1,nsrec)
        do j=1,nsrec
          k=k+1
          series(k)=ddat(j)
        enddo
        i=i+nchan
        l=l+1
        if (mod(l,8).eq.0) write(llog,*) uprow,'Read Rec.#:',l
      enddo
 1    close(lun)
      write(llog,*) uprow,'Read',l,' records.'
      ntim=k
      end
c==============================================================================


