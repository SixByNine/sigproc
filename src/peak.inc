	logical wepn
	common /wrepn/ wepn

	real snrmin
	common /snrth/ snrmin

	logical pout,xout
	common /psout/ pout,xout

	integer nsave
	parameter(nsave=1024)
	real smoothed
	common /smser/ smoothed(nsave)

	integer pulidx,nfound,nloop
	real    pulsnr
	common /pulse/ nloop,nfound,pulidx(nsave),pulsnr(nsave)

	integer nbins
	parameter(nbins=128)

	logical pzero
