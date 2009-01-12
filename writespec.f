c=============================================================================
      subroutine writespec(llg,pfile,fold,ssnr,dm,ac,tsamp,npf)
c=============================================================================
c
c     Writes out a power spectrum 
c
c     99/07/12 - dunc@naic.edu -- added dm and ac to header info
c      
c=============================================================================
      implicit none
      real dm, ac
      integer fold,npf
      real*8 tsamp
      real ssnr(*)
      integer i,lun,llg
      character*9 fname
      character*(*) pfile
C      byte bdat(2**23)
C      real scale,offset

      write(fname,'(''fold'',i1,''.spc'')')fold

      call glun(lun)
      open(unit=lun,file=pfile,status='unknown',form='unformatted')
c      open(unit=lun,file=fname,status='unknown',form='unformatted')
      write(lun) dm, ac ! new 
      write(lun) tsamp,npf,fold
      write(lun) (ssnr(i),i=1,npf)
C      call reatfbin(ssnr,npf,bdat,scale,offset,1)
C      write(lun) scale,offset
C      write(lun) (bdat(i),i=1,npf)
      close(unit=lun)
c      write(llg,*) 'Dumped fold to file: ',fname
      write(llg,'(a,a)') ' Dumped fold to file: ',
     &    pfile(1:index(pfile,' ')-1)
      end
