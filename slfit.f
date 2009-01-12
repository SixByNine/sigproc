c==============================================================================
      subroutine slfit(x, y, z, ndata, sig, weight, a, b, siga, sigb)
c==============================================================================
c
c     Fits the data passed down in arrays x(ndata) y(ndata) to the
c     straight line y=a+bx by minimising chi**2. Standard deviations
c     in y are passed down by sig(ndata) and can be weighted into
c     the fit if the logical switch weight is on. a and b together
c     with their uncertainties siga and sigb are returned.
c     Adapted from the fit routine given in numerical recipes to
c     be used in the psrflux software.
c     DRL 93/06/01 @ JB
c
      logical weight, useweight
      dimension x(ndata),y(ndata),z(ndata),sig(ndata)

      sx=0.
      sy=0.
      st2=0.
      b=0.
c
c     for upper limits opt for an unweighted fit
c     and use only three quarters the flux values
c
      useweight=weight
      do i=1,ndata
         z(i)=y(i)
         if (sig(i).eq.-999.0) then
           useweight=.false.
           z(i)=0.75*y(i)
         endif
      enddo

      if(useweight) then
        ss=0.
        do i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+z(i)*wt
        enddo
      else
        do i=1,ndata
          sx=sx+x(i)
          sy=sy+z(i)
        enddo
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(useweight) then
        do i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          b=b+t*z(i)/sig(i)
          st2=st2+t*t
        enddo
      else
        do i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*z(i)
        enddo
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.0
      if (.not.useweight) then
         do i=1,ndata
           chi2=chi2+(z(i)-a-b*x(i))**2
         enddo
         q=1
         sigdat=0.0
c
c        no uncertainty for an unweighted two point fit.
c
         if (ndata.gt.2)  sigdat=sqrt(chi2/(ndata-2)) 
         siga=siga*sigdat
         sigb=sigb*sigdat
      endif

      end
