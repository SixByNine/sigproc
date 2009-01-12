c=============================================================================
      subroutine readspec(sfile,fold,ssnr,dm,ac,tsamp,npf)
c=============================================================================
c
c     Reads in a power spectrum 
c 
c     99/07/12 - dunc@naic.edu -- added dm and ac to header info
c      
c=============================================================================
      implicit none
      character*80 sfile
      real dm, ac
      integer fold,npf
      real*8 tsamp
      real ssnr(*)
      integer i,lun
      character*80 fname
C      byte bdat(2**23)
C      real scale,offset
      logical filex

      if (sfile.ne.' ') then
         fname=sfile
      else
         write(fname,'(''fold'',i1,''.spc'')')fold
      endif

      write(*,*) 'Spectrum file: ',fname
      inquire(file=fname,exist=filex)
      if (.not.filex) stop 'File does not exist!'
      call glun(lun)
      open(unit=lun,file=fname,status='unknown',form='unformatted')
      read(lun) dm, ac ! new 
      read(lun) tsamp,npf,fold
      read(lun) (ssnr(i),i=1,npf)
C      read(lun) scale,offset
C      read(lun) (bdat(i),i=1,npf)
C      call reatfbin(ssnr,npf,bdat,scale,offset,-1)
      close(unit=lun)
      end
      
