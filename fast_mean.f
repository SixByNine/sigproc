      subroutine fast_mean(val,num,nav,meanval,ncal)
      
      implicit none

      integer num, nav, ncal, i1, i2, i, nonzero, navmin, j
      integer nonzero0
      real val(num), meanval(num), suboff, subon, sum, sum0

      navmin = 10 ! minimum number of sample to compute average of

      ncal = 0
      
c     --- this treats the middle array, apart from edge effects --
      sum = 0.0
      nonzero = 0.0
      do i=nav/2, nav+nav/2
         sum = sum + val(i)
         if (val(i).ne.0.0) nonzero=nonzero+1
      enddo
      sum0 = sum
      nonzero0 = nonzero
      
      i1 = nav
      i2 = num-nav/2
      do i=i1, i2
         suboff  = val(i-nav/2)
         subon   = val(i+nav/2)
         sum = sum - suboff + subon
         if (suboff.ne.0.0) nonzero=nonzero-1
         if (subon .ne.0.0) nonzero=nonzero+1
         meanval(i) = sum/float(nonzero)
         ncal = ncal + 1
      enddo


c     --- this deals with last part of the array --
      do i=i2+1, num-navmin
         suboff = val(i-nav/2)
         sum = sum - suboff
         if (suboff.ne.0.0) nonzero=nonzero-1
         meanval(i) = sum/float(nonzero)
         ncal = ncal + 1
      enddo
c     make sure that it doesn't get too noisy at the end:
      do i=num-navmin+1, num
         meanval(i)=meanval(num-navmin) 
         ncal = ncal + 1
      enddo
 

c     ----- this deals with first part:

c      find first nonzero bin while copying the zero to the mean array:
      i=1
 100  continue
        if (val(i).ne.0.0) goto 200
        meanval(i) = 0.0
        i=i+1
        if (i.lt.i1-1) goto 100
 200  continue

c     - fill the last bit until we reached the already treated array
      i2 = i1-1
      i1 = i
      sum=sum0
      nonzero = nonzero0
      do i=i2, i1+navmin, -1
         suboff = val(i+nav/2)
         sum = sum - suboff
         if (suboff.ne.0.0) nonzero=nonzero-1
         meanval(i) = sum/float(nonzero)
         ncal = ncal + 1
      enddo
c     make sure that it doesn't get too noisy at the end:
      do i=i1, i1+navmin-1
         meanval(i)=meanval(i1+navmin) 
         ncal = ncal + 1
      enddo

c      do i=i1, i2
c         meanval(i) = meanval(i2+1)
c         ncal = ncal + 1
c      enddo

      end







