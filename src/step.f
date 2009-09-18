c==============================================================================
	program step
c==============================================================================
c       Produces DM/AC files for up to two loops (dunc@mpifr-bonn.mpg.de)
c==============================================================================
	real x,x1,x2,xs,y,y1,y2,ys
	character*80 tmp
	integer narg,iargc
	logical inner
	narg=iargc()
	if (narg.lt.3) stop 'Specify min max step'
	call getarg(1,tmp)
	read(tmp,*) x1
	call getarg(2,tmp)
	read(tmp,*) x2
	call getarg(3,tmp)
	read(tmp,*) xs
	call getarg(4,tmp)
	inner=tmp.ne.' '
	if (inner) then
	   read(tmp,*) y1
	   call getarg(5,tmp)
	   read(tmp,*) y2
	   call getarg(6,tmp)
	   read(tmp,*) ys
	endif
	x=x1
	do while (x.lt.x2+xs) 
	  if (inner) then
             y=y1	
	     do while (y.lt.y2+ys)
	     write(*,'(2(f8.3,3x))') x,y
	     y=y+ys
	     enddo
	  else
	     write(*,'(f8.3)') x
	  endif
	  x=x+xs
	enddo
	end
c==============================================================================
