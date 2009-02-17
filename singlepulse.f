      subroutine singlepulse(llog,append,spthresh,ncandsmax,nsmax)

      implicit none
      include 'seek.inc'
C
C Variables for single-pulse search
      logical append
      integer done_pulse, pulse, maxloopsize
      parameter (maxloopsize=10485760)
      integer loopsize,lsd1,llog,npulses
      real realdata(maxloopsize)
      character best_pulses*132, hist_pulses*132, scrdisk*80

c  initialize parameters for single-pulse search and write headerfile

      best_pulses=filename(1:lst)//'.pls'
      hist_pulses=filename(1:lst)//'.hst'

      dmidx=0
      if (append) dmidx=1
      loopsize = ntim2/64
      scrdisk='./'
      lsd1=2
      if (dmidx.eq.0) 
     .write(llog,*) 'Searching for single pulses... (thresh=',spthresh,
     .		     ', ncandsmax=',ncandsmax,', nsmax=',nsmax,')'

c       Search for pulses in this DM file and write out to file best.tmp.

	ntim2=ntim2-10./tsamp
      done_pulse = pulse(realdata, spthresh, nsmax,
     .     dmidx, ntim2, series, loopsize, scrdisk, lsd1)

c       Find best smoothing index for each pulse and append pulses to file
c       best_pulses.

      call best_smoothing(ncandsmax, dmidx, refdm, ntim2,
     .     best_pulses, hist_pulses, spthresh, nsmax, scrdisk,
     .     npulses,tsamp)
      write(llog,*) npulses, ' pulses found...'


      return
      end
 	subroutine best_smoothing(ncandsmax, dmidx, dm, length,
     .			           best_pulses, hist_pulses, thresh,
     .			 	   nsmax, scrdsk1, npulses, tsamp)
      implicit none

c reads output file from pulse.c that is 
c sorted according to time bin number.
c For a given time bin, it finds the smoothing
c that gives the largest s/n.
 
c program iterates nitmax times since it only looks at
c pairs of events in each iteration
 
c jmc 9 May 1995

      integer maxcands
      parameter (maxcands = 2**21)

      integer ndm, ns, timebin, ncandsmax, dmidx, npulses
      real snr, mean, rms, dm
      real*8 tsamp, time, width
      integer ndm1, ns1, timebin1, npulse, n_to_print
      real snr1, mean1, rms1, thresh

      integer nsmax

      integer length, power, lun
      real rlength

      integer best_dm(maxcands), best_ns(maxcands),
     .	      best_time(maxcands), best_snr(maxcands)
      integer itolerance, io_error, isamp
      character*80 inline, inline1
      character*132 scrdsk1
      data inline/' '/, inline1/' '/

      character *(*) best_pulses, hist_pulses

      integer ind
      logical latest

      integer lnblnk			!fn to find line length
      integer lunin, lunout

      integer nit, nitmax
      data nitmax/4/
      integer l0,l1,l2,l3,l4,l5,i
      npulses = 0


c      lun=index(scrdsk1,' ')-1        ! MBy: for accomplishing modif done in
      lun=index(scrdsk1,' ')-2         ! pulse.c of the string scrdsk1
   
      call system('sort -n -k 3 '//scrdsk1(1:lun)//'/'//'best_tmp >
     &		 '//scrdsk1(1:lun)//'/'//'best_tmp.sorted')

      call glun(l0)
      open(l0, file=scrdsk1(1:lun)//'/'//'best_tmp.sorted')

      call glun(l1)
      open(l1, file=scrdsk1(1:lun)//'/'//'best_tmp1')

      call glun(l2)
      open(l2, file=scrdsk1(1:lun)//'/'//'best_tmp2')

      do nit = 1, nitmax

       if(nit .eq. 1) then		!first time, use orig data file
          lunin = l0
          lunout = l1
       else
          if(mod(nit,2) .eq. 0) then    !other iterations: flip flop 
             lunin = l1
             lunout = l2
          else
             lunin = l2
             lunout = l1
          endif
       endif

       rewind(lunin)
       rewind(lunout)
       read(lunin, "(a80)",end=1000) inline
       ind = lnblnk(inline)-1

    1 continue
      read(inline, *)  
     .   ndm, ns, timebin, snr, mean, rms
      read(lunin, "(a80)",end=998) inline1
      latest = .false.
      read(inline1, *) 
     .   ndm1, ns1, timebin1, snr1, mean1, rms1
      itolerance=(2.**ns+2.**ns1)/2.+1
      if(abs(timebin1-timebin).le.itolerance) then 
           latest = .true.		!the better of the pair is written 
           if(snr1 .gt. snr) then 
             ind = lnblnk(inline1)
             write(lunout,*) inline1(1:ind)
           else
             ind = lnblnk(inline)
             write(lunout,*) inline(1:ind)
           endif
           read(lunin, "(a80)",end=999) inline
           latest = .false.
      else
           ind = lnblnk(inline)
           write(lunout,*) inline(1:ind)
      inline = inline1
      endif
      go to 1

c always save last line

  998 continue
      ind = lnblnk(inline)
      write(lunout,*) inline(1:ind)

  999 continue 

      enddo

      lunin = lunout
      rewind(lunin)
    2 continue
      read(lunin,"(a80)", end=9999) inline
      ind = lnblnk(inline)

      read (inline(1:ind),*) ndm1, ns1, timebin1, snr1, mean1, rms1
      if (timebin1.lt.length-10000) then
      npulses = npulses + 1
      best_dm(npulses) = ndm1
      best_ns(npulses) = ns1
      best_time(npulses) = timebin1
      best_snr(npulses) = snr1*100
      end if

      go to 2
 9999 continue
       
 1000	close(l0)
	close(l1)
	close(l2)

	if (npulses.gt.0) then

	call sort3(npulses,best_snr,best_time,best_ns)
        if (dmidx.eq.1) then
	   call glun(l3)
           open(l3, file=best_pulses,access='append')
        else
	   call glun(l3)
           open(l3, file=best_pulses)
        end if

	if (dmidx.eq.1) then
	   call glun(l4)
	   open (l4, file=hist_pulses,access='append')
	else
	   call glun(l4)
	   open (l4, file=hist_pulses)
	end if
	
	if (ncandsmax.gt.npulses) then 
 	   n_to_print = npulses
	else
	   n_to_print = ncandsmax
	end if

        rlength = length
        power = log(rlength)/log(2.) + 0.5

	isamp = tsamp/1.e-6
	if (dmidx.eq.0) 
     .	write (l3,100) isamp,thresh,nsmax,ncandsmax
	do npulse = npulses,npulses-n_to_print+1,-1
	   write (l3,200) dm, best_ns(npulse),best_time(npulse),best_snr(npulse)/100.,power
	end do

	write (l4,300) dm, npulses
	close (l3)
	close (l4)

 100	format('tsamp=',i4,' us, thresh= ',f4.1,' sigma, nsmax=',i3,', ncandsmax=',i8)
 200	format(f9.3,3x,i2,3x,i10,3x,f7.2,3x,i2)
 300	format(f9.3,3x,i8)
	else 
	   n_to_print = 0
	end if
	call system('rm best_tmp best_tmp1 best_tmp2 best_tmp.sorted')

  	 end
      SUBROUTINE sort3(n,arr,brr,crr)
      INTEGER n,M,NSTACK
      INTEGER arr(n),brr(n),crr(n) !changed fron REAL (MBy 01/11/01)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,c,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
	  c=crr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
	    crr(i+1)=crr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
	  crr(i+1)=c
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
	temp=crr(k)
	crr(k)=crr(l+1)
	crr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
	  temp=crr(l+1)
	  crr(l+1)=crr(ir)
	  crr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
	  temp=crr(l)
	  crr(l)=crr(ir)
	  crr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
	  temp=crr(l+1)
	  crr(l+1)=crr(l)
	  crr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
	c=crr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
	temp=crr(i)
	crr(i)=crr(j)
	crr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
	crr(l)=crr(j)
	crr(j)=c
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0!5,.
