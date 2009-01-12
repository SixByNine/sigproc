	subroutine quikgray(dat,nxd,nyd,nx,ny)

c  Coarse grey-scale plot using PG plot routines
c  Assumed that viewport and window defined outside routine
        integer nxd,nyd
        integer nsym
        parameter(nsym=10)
	integer*4 ksym(nsym)
	real*4 dat(nxd,nyd),tr(6),xw,yw
	data ksym/1,20,21,2,3,14,15,17,16,18/

	s=0.
	ss=0.
	smin=1.e30
	smax=-smin
	do j=1,ny
	  do i=1,nx
	    aa=dat(i,j)
	    s=s+aa
	  enddo
	enddo
	s=s/(nx*ny) 
	do j=1,ny
	  do i=1,nx
            aa=dat(i,j)-s
	    ss=ss+aa**2
	    smax=max(smax,dat(i,j))
	    smin=min(smin,dat(i,j))
	  enddo
	enddo
	rms=sqrt(ss/(nx*ny))
        
	tr(1)=0.0 
	tr(2)=1.0/real(nx)
	tr(3)=0.0
	tr(4)=0.0
	tr(5)=0.0
	tr(6)=1.0
	xmin=xw(tr,0.5,0.5)
	ymin=yw(tr,0.5,0.5)
	xmax=xw(tr,real(nx)+0.5,real(ny)+0.5)
	ymax=yw(tr,real(nx)+0.5,real(ny)+0.5)
	call pgwindow(xmin,xmax,ymin,ymax)
        call pggray(dat,nxd,nyd,1,nx,1,ny,smax,s+0.5*rms,tr)

	return

	do j=1,ny
	  do i=1,nx
	    k=min(int((dat(i,j)-s)/rms),nsym)
	    if(k.gt.0)then
	      x=i
	      y=j
	      call pgpoint(1,x,y,ksym(k))
	    endif
	  enddo
	enddo

	end


	real function xw(tr,i,j)
	real tr(6),i,j
	xw=tr(1)+tr(2)*i+tr(3)*j
	end

	real function yw(tr,i,j)
	real tr(6),i,j
	yw=tr(4)+tr(5)*i+tr(6)*j
	end

