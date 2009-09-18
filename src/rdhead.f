c==============================================================================
      subroutine rdhead(hdrfile,pspmfile,f0,chbw,tsamp,nsrec,nchan)
c==============================================================================
      implicit none
      character*(*) hdrfile,pspmfile
      real f0,chbw,tsamp
      integer nsrec,nchan,lun,idx
      character*80 line

      call glun(lun)
      open(unit=lun,file=hdrfile,status='old',err=999)
      read(lun,'(a)') line  ! filename
      idx=index(line,':')+1
      pspmfile=line(idx:)
      read(lun,'(a)') line  ! tsamp (us)
      idx=index(line,':')+1
      read(line(idx:),*) tsamp
      tsamp=tsamp*1.0e-3
      read(lun,'(a)') line  ! nsrec
      idx=index(line,':')+1
      read(line(idx:),*) nsrec
      read(lun,'(a)') line  ! fcent (MHz)
      idx=index(line,':')+1
      read(line(idx:),*) f0  
      read(lun,'(a)') line  ! channel band (kHz)
      idx=index(line,':')+1
      read(line(idx:),*) chbw
      chbw=chbw/1000.0
      read(lun,'(a)') line  !  nchans
      idx=index(line,':')+1
      read(line(idx:),*) nchan
      close(unit=lun)
      return
 999  stop 'Header file not found!'
      end
