c==============================================================================
	real function histmax(ndata, data, datmin, datmax, nbin)
c==============================================================================
c
c	This function returns the maximum binheight in the array "data"
c	having "ndata" items for binning over the range "datmin" -> "datmax"
c       "nbin" signifies the number of bins to be used.
c
c       Last modified Mon Apr 19 1999 (dunc@naic.edu) -> Took old histmax
c       Binning routine is now called seperately by histmax.
c
c==============================================================================
c
	implicit none
	integer ndata
	real data (ndata), datmin, datmax
	integer nbin, ibin
	integer num(256)
c
c       Get the bin values - now packaged into seperate routine...
c
	call histval(ndata,data,datmin,datmax,nbin,num)
c
c       Find the maximum value
c
	histmax = 0
	do ibin=1,nbin
	  histmax = max(histmax,real(num(ibin)))
	end do
c
c       Increase it by 10% for plot scaling...
c
	histmax = histmax * 1.1
c
c       Job Done!
c
	end
c==============================================================================
	subroutine histval(ndata, data, datmin, datmax, nbin, binval)
c==============================================================================
c
c	This subroutine returns the binned values for a data set passed down
c       in the real*4 array "data(ndata)" [i.e. ndata items]. The binning
c       is carried out for "nbin" bins defined between "datmin" and "datmax".
c       The binned values are returned in the integer*4 array "binval(ndata)"
c
c       Created Mon Apr 19 1999 (dunc@naic.edu) -> Adapted from old histmax
c
c==============================================================================
c
	implicit none
c
c       Variables passed down
c
	integer ndata
	real data (ndata), datmin, datmax
c
c       Returned
c
	integer binval(256)
c
c       Local variables...
c
	integer nbin, ibin, i
c
c       First, Zero the bin values
c
	do ibin=1,nbin
	  binval(ibin) = 0
	end do
c
c       Now bin the data...
c
	do i =1,ndata
	  ibin = (data(i)-datmin)/(datmax-datmin)*nbin+1
	  if (ibin.ge.1.and.ibin.le.nbin) binval(ibin)=binval(ibin)+1
	end do
c
c       Job Done!
c
	end
c
c==============================================================================
