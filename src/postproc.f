c==============================================================================
	program postproc
c==============================================================================

	parameter (maxdm = 2000, maxpulses=100000)

	integer npulses_orig(maxdm), npulses(maxdm)
  	real rms(maxdm), dm(maxdm), snr
	integer index, ns, length(maxdm), ncandsmax, nsmax
	integer iterate
	real thresh
	character command_line*80
	
	call getarg(1,command_line)
	if (command_line.eq.'') then
	   write(*,*) 'usage: postproc sample_time_in_seconds'
	   stop
	endif
	read(command_line,*) samp_int

	open (17,file='pulse_stats')
	open (27,file='pulse_best')
	open (37,file='dm.hist')

	read (17,*) ncandsmax
	read (17,*) thresh
	read (17,*) nsmax
	read (17,*) iterate

	print *, 'thresh: ', thresh
	print *, 'ncandsmax: ', ncandsmax
	print *, 'iterate: ', iterate 
	print *, 'nsmax: ', nsmax
	i = 1

 10	read (17,*,end=20) dm(i), npulses_orig(i), rms(i), length(i)
	if (npulses_orig(i).lt.ncandsmax) then
	   npulses(i) = npulses_orig(i)
	else
	   npulses(i) = ncandsmax
	end if
	if (i.gt.1) then
	   write (37,*) i, 0, npulses_orig(i),
     .			  npulses_orig(i),length(i),rms(i)
	   write (37,*) i+1, 0,  npulses_orig(i),
     .			 npulses_orig(i),length(i),rms(i)
	else
           write (37,*) i, 0, npulses_orig(i),
     .			  npulses_orig(i),length(i),rms(i)
           write (37,*) i+1, 0,  npulses_orig(i),
     .			 npulses_orig(i),length(i),rms(i)
	   length1 = length(i)
	endif
	i = i + 1
	goto 10

 20	close (17)
	close (37)

	ndm = i - 1
	

	open (17,file='dmlist')
	open (37,file='best')
	open (47,file='best.high')

	write (17,*) ndm
	timestart= 0.0

	npulse = 0
	do i = 1, ndm
	   write (17,100) dm(i)
	   do j = 1, npulses(i)
	      read (27,*) ns, index, snr
	      time = index*samp_int
	      write (37,200) i-1, ns, index, snr, rms(i)
	      if (time.ge.timestart.and.snr.ge.5)
     .		 write (47,*) i-1, ns, time - timestart, snr
	      npulse = npulse + 1
	   end do
	end do

 100	format(f10.4)
 200	format(i3,2x,i1,2x,i8,2x,f6.2,2x,2x,f7.5)

	close (17)
	close (27)
	close (37)
	close (47)

	end
c==============================================================================
