c==============================================================================
      subroutine zapit(llog,fold,zapfile,spec,npts,tsamp)
c==============================================================================
c
c     Zeros portions of the amplitude spectrum that are regularly occupied
c     by RFI signals: the mains power line etc. The frequency ranges are read
c     from a free-format ASCII file containing, per line, the minimum
c     and maximum frequency to zap in the spectrum.
c
c     llog    - i4  - logical unit number for all but warning messages.
c     fold    - i4  - harmonic fold index
c     zapfile - ch  - name of the ASCII file containing zap-freq ranges.
c     spec    - r4  - array containing the amplitude spectrum.
c     npts    - i4  - number of points in the amplitude spectrum.
c     tsamp   - r8  - sampling time in the time domain (s).
c
c     Created: November 1997 (dunc@mpifr-bonn.mpg.de)
c      
c     Modification history:
c      
c     98/04/10 (drl@jb.man.ac.uk) works out the bin numbers to zap (faster!)
c     98/04/28 (dunc@mpifr-bonn.mpg.de) just return if zapfile is empty.
c     98/04/29 (dunc@mpifr-bonn.mpg.de) range check on bin numbers.
c     98/04/29 (dunc@mpifr-bonn.mpg.de) fold number passed down.
c     02/02/26 (drl@jb.man.ac.uk) widths of filters in proportion to frequency
c
c     To do:
c
c     Something more sensible... e.g. what the GW people are proposing!
c      
c==============================================================================
c      
      implicit none
      character*(*) zapfile
      real spec(*)
      real*8 tsamp
      integer npts,llog,fold
c      
c     local variables...
c
      character*80 line
      integer nb
      real blo,bhi,bf,df
      integer i,j,lun,nlo,nhi,fbin,nh,nz
      logical first,static
      data first/.true./
      save
c
c     Initialize
c      
      nz=0
      nb=0
c
c     Get a free logical unit number
c
      call glun(lun)
c
c     Attempt to open the file and read in the birdies
c      
      open(unit=lun,file=zapfile,status='old',err=3)
      do while(.true.)
         read(lun,'(a)',end=1) line
	 if (line.eq.' ') return      ! return if file empty
         if (line(1:1).ne.'#') then
            read(line,*,err=4) blo,bhi,nh
	    write(llog,*)
     &      'Filter:',blo,' ->',bhi,' Hz.',nh,' harmonic(s)'
            nb=nb+1
            bf=0.5*(bhi+blo)
            df=0.5*(bhi-blo)
            if (nh.lt.0) then
		static=.true.
		nh=-1*nh
            else
		static=.false.
            endif
            do i=1,nh
	       if (static) then
                  blo=bf*real(i)-df
                  bhi=bf*real(i)+df
                else
                  blo=(bf-df)*real(i)
                  bhi=(bf+df)*real(i)
	        endif
c
c              Convert frequencies to bin numbers
c
               nlo=fbin(tsamp,npts,fold,blo)
               nhi=fbin(tsamp,npts,fold,bhi)
c
c              Zero spectrum within these ranges provided
c              that they are sensible (i.e. within npts)
c
                do j=nlo,nhi
                   if (spec(j).ne.0.0.and.j.le.npts.and.j.ge.1) then
                     nz=nz+1
                     spec(j)=0.0
                   endif
                enddo
            enddo
         endif
      enddo
 1    close(unit=lun)
      write(llog,'(i5,a,f7.3,a)') nz,' spectral bins blown away (',
     & 100.0*real(nz)/real(npts),'% of spectrum)'
 2    format(i3,a,a,a)
      return
c
c     issue warning message (first time only) if the file was not opened
c      
 3    if (first) write(*,*) 'WARNING - birdie file not found...'
      first=.false.
      return
c
c     stop here if an error occured when reading a filter
c
 4    write(*,'('' Error reading filters from: '',a)')zapfile
      write(*,*) 'Please supply filters in the following way:'
      write(*,*) '#flo (Hz) fhi (Hz) nharm (integer)'
      write(*,*) '49.9 50.1 1'
      write(*,*)
      write(*,*) 'The # mark can be used to comment the file'
      write(*,*) 'Happy zapping!'
      stop
      end
c      
c=============================================================================
