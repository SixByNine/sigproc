c=============================================================================
      subroutine srd(infile,header,ntim,mpts,tsamp,refdm,refac,series)
      implicit none
      character*80 header,infile
      integer ntim,mpts
      real tsamp,refdm,refac,series(mpts)
      integer lun,i
      call glun(lun)
      open(lun,file=infile,status='old',form='unformatted',err=999)
      read(lun) header
      read(lun) ntim,tsamp,refdm,refac
      if (ntim.gt.mpts) then
	write(*,*) 'WARNING - too many points in data file!'
	write(*,*) '** reading in max possible:',mpts/1024,' kpts'
	ntim=mpts
      endif
      read(lun) (series(i),i=1,ntim)
      close(unit=lun)
      return
 999  stop 'Error opening input file!'
      end
c=============================================================================
      subroutine swr(oufile,header,ntim,mpts,tsamp,refdm,refac,series)
      implicit none
      character*80 header,oufile
      integer ntim,mpts
      real tsamp,refdm,refac,series(mpts)
      integer lun,i
      call glun(lun)
      open(unit=lun,file=oufile,status='unknown',
     &     form='unformatted')
      write(lun) header
      write(lun) ntim,tsamp,refdm,refac
      write(lun) (series(i),i=1,ntim)
      close(unit=lun)
      end
c==============================================================================
      subroutine ope(infile,lun,header,ntim,mpts,tsamp,refdm,refac)
      implicit none
      character*80 header,infile
      integer ntim,mpts
      real tsamp,refdm,refac
      integer lun
      call glun(lun)
      open(lun,file=infile,status='old',form='unformatted',err=999)
      read(lun) header
      read(lun) ntim,tsamp,refdm,refac
      if (ntim.gt.mpts) then
	write(*,*) 'WARNING - too many points in data file!'
	write(*,*) '** will read in max possible:',mpts/1024,' kpts'
	ntim=mpts
      endif
      return
 999  stop 'Error opening input file!'
      end
c=============================================================================
      integer function glp(lun,gulp,ntim,gmean)
      implicit none
      integer lun,ntim,i
      real gulp(*),gmean
      gmean=0.0
      glp=0
      read(lun,end=999) (gulp(i),i=1,ntim)
      do i=1,ntim
         gmean=gmean+gulp(i)
      enddo
      gmean=gmean/real(ntim)
      return
 999  glp=-1 
      end
c=============================================================================
      
