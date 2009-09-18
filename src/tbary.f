      subroutine tbary(ra,dec,topo,obs,mjd)
c
c     quick utility to return the barycentric arrival time
c     by calling tempo given RA and DEC and a topocentric time
c
      implicit none
      character*80 ra,dec,topo,obs
      double precision mjd,pfphi,pfsec,orphi,freq,wgt,etoa,prephi,tmp
      
      if (obs.eq." ") stop "usage: tssb RA DEC MJD SITE"

      open(10,file="tssb.par",status="unknown")
      write(10,'(''PSR 0000+00'')')
      write(10,'(''RAJ '',a)') ra(1:index(ra," ")-1)
      write(10,'(''DECJ '',a)') dec(1:index(dec," ")-1)
      write(10,'(''F0 1.0'')')
      write(10,'(''DM 0.0'')')
      write(10,'(''PEPOCH '',a)') topo(1:index(topo," ")-1)
      close(10)

      read(topo,*) mjd
      open(10,file="tssb.tim",status="unknown")
      write(10,'(a1,7x,a7,1x,f8.3,1x,f19.13,1x,i4,3x,f6.2)')
     &       obs(1:1),"0000+00",9999.0,mjd,1024,10.0
      close(10)

      call system("tempo tssb.tim > /dev/null")
      open(10,file="resid2.tmp",form="unformatted")
      read(10) mjd,pfphi,pfsec,orphi,freq,wgt,etoa,prephi,tmp
      close(10)
      write(*,'(f19.13)') mjd

      call system("rm -f tssb.tim tssb.par tempo.lis resid2.tmp")
      call system("rm -f 0000+00.par")

      end
