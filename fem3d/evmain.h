
c------------------------------------------------
c finite element information for each element
c------------------------------------------------
c
c please change this info also in main

	integer evdim
	parameter ( evdim = 16 )

	!real ev(evdim,neldim)
	double precision ev(evdim,neldim)
	common /ev/ev

c------------------------------------------------

