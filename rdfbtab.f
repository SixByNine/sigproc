	subroutine rdfbtab(skyfreq,maxchans,nchans)
	implicit none
c
c	Reads in the details of the filterbank contained in
c       the ASCII file fb.tab
c
	integer maxchans,nchans
	real skyfreq(maxchans)

	integer lun,i

	call glun(lun)
	nchans=maxchans
	open(unit=lun,file='fb.tab',status='old',err=999)
	do i=1,maxchans
	  read(lun,*,end=1) skyfreq(i)
        enddo
 1	close(unit=lun)
	if (i.lt.maxchans) nchans=i-1
	write(*,*) 'The filterbank has',nchans,' channels...'
	do i=1,nchans
	  write(*,*) 'Chan #',i,' Freq:',skyfreq(i),' MHz'
        enddo
	return
c
c	Couldn't find filterbank table... Tell the user what to do...
c	
 999    write(*,*) 'You have not made a filterbank table: "fb.tab"'
	write(*,*) 'This should be a free-format ASCII listing'
	write(*,*) 'For example:'
	write(*,*) '1424.0'
	write(*,*) '1420.0'
	write(*,*) '1416.0'
	write(*,*) '1412.0'
	write(*,*) 'would be the set up for a four channel filterbank.'
	write(*,*) 'For the PSE in Effelsberg, all numbers should be'
	write(*,*) 'identical and equal to the centre frequency used.'
	write(*,*)
	write(*,*) 'All numbers are sky frequencies in MHz. Good luck!'
	stop
        end 
