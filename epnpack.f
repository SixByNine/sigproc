c==============================================================================
      subroutine hex2b10(hexstr,b10no)
c==============================================================================

      implicit none

      integer b10no
      character*4 hexstr
      character*1 hex(16)
      integer i,j,dig
      data hex/'0','1','2','3','4','5','6','7','8','9',
     &         'A','B','C','D','E','F'/

      b10no=0
      dig=0

      do i=4,1,-1
         do j=1,16
            if (hexstr(i:i).eq.hex(j)) goto 5
         enddo
 5       continue
         b10no=b10no+int(16.0**float(dig)*float(j-1))
         dig=dig+1
      enddo

      end
c==============================================================================
      subroutine b102hex(b10no,hexstr)
c==============================================================================
c
c	Converts a base 10 number passed down as the integer "b10no"
c	to a hexadeximal string returned as the character*4 "hexstr"
c       N.B. maximum integer for 4 character hex string is 65535
c
      implicit none

      integer b10no
      character*4 hexstr
      character*1 hex(16)
      integer i,j,no,num,dig
      data hex/'0','1','2','3','4','5','6','7','8','9',
     &         'A','B','C','D','E','F'/

      no=b10no
      dig=0

      do i=4,1,-1
         do j=15,0,-1
            num=int(float(j)*16.0**float(i-1))
            if (num.le.no) goto 5
         enddo
 5       continue
         dig=dig+1
         hexstr(dig:dig)=hex(j+1)
         no=no-num
      enddo

      end
c==============================================================================
      subroutine rawtfint(raw,nbin,intdata,scale,offset,rms,dirn,maxint)
c==============================================================================
c
c     Routine to convert "raw" (i.e. floating point) to/from integer
c     format. In this case, the integer may vary between 0 and maxint.
c     Raw data is passed down from site software and stored as integers
c     in the EPN files.
c
c     The direction of the conversion is controlled by "dirn":
c
c     dirn = +1 : raw -> int
c     dirn = -1 : int -> raw
c
c==============================================================================
c    
      implicit none
c
c     Passed down variables
c
      real raw(*)
      real*8 scale, offset, rms
      integer nbin, intdata(*), dirn, maxint
c
c     Local variables
c
      real dmin, dmax, sumsq, pmax, sum, mean, ndiv
      integer i, ibmax, nshift
      if (dirn.eq.1) then
c
c       Scale the data:
c
c       Find minimum & maximum values of the data
c
        dmin = raw(1)
        dmax = raw(1)
        do i=1,nbin
          dmax=max(raw(i),dmax)
          dmin=min(raw(i),dmin)
        enddo
        offset=dmin
        scale=float(maxint)/(dmax-dmin)
c
c       Now do the scaling
c
        do i=1,nbin
          intdata(i)=int((raw(i)-offset)*scale)
        enddo
c
c-------------------------------------------------------------
c       If no rms has been supplied by the user, have a go at 
c       calculating it from the wings of the profile after 
c       shifting it so that its peak lies at bin number nbin/2
c-------------------------------------------------------------
c
        if (rms.eq.0.0) then
          pmax=-1.0e32
          ibmax=0
          do i=1,nbin
             if (raw(i).gt.pmax) then
                ibmax=i
                pmax=raw(i)
             endif
          enddo
          nshift=nbin/2-ibmax
          call sprof(raw,nbin,nshift)
          sumsq=0.0
          sum=0.0
          do i=1,nbin/15
            sumsq=sumsq+raw(i)*raw(i)
            sumsq=sumsq+raw(nbin-i+1)*raw(nbin-i+1)
            sum=sum+raw(i)
            sum=sum+raw(nbin-i+1)
          enddo
	  ndiv=max(1.0,float(2*nbin/15))
          rms=sqrt(sumsq/ndiv)
          mean=sum/ndiv
        endif
      else if (dirn.eq.-1) then
c
c       Prepare to de-scale the data:
c
        do i=1,nbin
          raw(i)=real(offset+dble(intdata(i))/scale)
        enddo
      else
c
c       Silly option was passed down...
c
        write(*,*) 'Invalid value for dirn passed to rawtfint : ',dirn
        stop
      endif
c
c     Job done!
c
      end
c==============================================================================
      integer function nepnrec(filename)
c==============================================================================
c
c     Function to find out the number of records in an EPN file.
c
c     This can in principle be done with an INQUIRE statement but
c     it is system dependent. The following approach uses a simple
c     binary search algorithm to find the last record which is far
c     quicker than reading in every record in very large EPN files.
c
c     if the file doesn't exist nepnrec is returned as 0
c     nepnrec = -1 if an error occurred whilst reading the file
c
c     Created 96/12/03 DRL@MPIfR
c     Modified October 1997 to return -1 on error
c
c
      implicit none
      character*80 filename

      character lin1*12,recrd*80
      integer irec,istat,counter,recln,hirec,lorec,orec,lun
      logical filex
c
c     Check to see whether the EPN file exists...
c     return nepnrec=0 if it doesn't
c
      nepnrec=0
      inquire(file=filename,exist=filex)
      if (.not.filex) return 
c
c     It exists! Open it and find out the true record length..
c
      call glun(lun)
      nepnrec=-1 ! This will be the value of nepnrec if an error occurs
      irec=-1
      open(unit=lun,file=filename,status='old',access='direct',
     &     form='unformatted',recl=12,err=999)
      read(lun,rec=1,err=999)lin1
      read(lin1,'(8x,i4)',err=999)counter
      close(unit=lun)
      recln=counter*80
      if (recln.le.0) return ! Error
c
c     Find out crude lower and upper bounds for the last record
c
      open(unit=lun,file=filename,status='old',access='direct',
     &     form='unformatted',recl=recln,err=999)
      irec=1
      istat=0
      do while(.true.)
         read(lun,rec=irec,iostat=istat) recrd
         if (istat.ne.0) goto 998
         irec=irec*2
      enddo
 998  hirec=irec
      lorec=irec/2
c
c     Now search for the last record within these bounds iteratively
c
      istat=0
      orec=0
      do while(orec.ne.irec)
         orec=irec
         irec=lorec+(hirec-lorec)/2
         read(lun,rec=irec,iostat=istat) recrd
         if (istat.eq.0) then
            lorec=irec
         else
            hirec=irec
         endif
      enddo
 999  continue
      close(unit=lun)
      nepnrec=irec
      end
c============================================================================= 
      subroutine dattim(date,hh,mm,ss)
c============================================================================= 
c
c     Returns the time obtained by truncating the MJD passed down
c     as date
c
      implicit none 
      double precision date,dummy,ss
      integer hh,mm


      dummy=abs(date-int(date))

      dummy=dummy*24.0
      hh=int(dummy)

      dummy=abs(dummy-int(dummy))

      dummy=dummy*60.0
      mm=int(dummy)

      dummy=abs(dummy-int(dummy))

      dummy=dummy*60.0
      ss=dummy

      end
