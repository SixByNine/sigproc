c==============================================================================
	subroutine baseline(llog)
c==============================================================================
c
c       Subtracts a SLOPING baseline from the time series contained in the
c       array series(ntim). Writes out the mean which is then subtracted
c       The time series is then normalized so that the rms is unity.
c
c	22/03/98 (dunc@mpifr-bonn.mpg.de)
c
c==============================================================================
	implicit none
        include 'time.inc'
	integer llog,i,j,nrun,nrav
	parameter(nrav=32)
	real mean,sums,sumd,rssq(nrav),rsum(nrav),var,mea,msq,rms
	real x(nrav),y(nrav),z(nrav),e(nrav),slope,inter,eslo,eint
	
	sumd=0.0
	do i=1,nrav
	   rsum(i)=0.0
	   rssq(i)=0.0
	enddo
	nrun=ntim/nrav
	
	do i=1,ntim
	  j=min(nrav,i/nrun+1)
	  rsum(j)=rsum(j)+series(i)
	  rssq(j)=rssq(j)+series(i)*series(i)
	  sumd=sumd+series(i)
        enddo
	mean=sumd/real(ntim)

	do i=1,nrav
	   mea=rsum(i)/real(nrun)
	   msq=rssq(i)/real(nrun)
	   var=sqrt(msq-mea*mea)
	   j=nrun/2+(i-1)*nrun
	   x(i)=i
	   y(i)=mea
	   e(i)=var/sqrt(real(nrun))
	enddo
	call slfit(x,y,z,nrav,e,.false.,inter,slope,eint,eslo)
	write(llog,*) 'Subtracting mean from function:',inter,slope

	do i=1,ntim
	  mean=inter+slope*real(i)*real(nrav)/real(ntim)
	  series(i)=series(i)-mean
	enddo
	
	sums=0.0
	do i=1,ntim
	  sums=sums+series(i)*series(i)
        enddo
	rms=sqrt(sums/real(ntim))

	do i=1,ntim
          series(i)=series(i)/rms
        enddo

        write(llog,*) 'Time series now has zero mean and unit rms...'
	end
