c
c $Id: basin.h,v 1.2 2007-03-20 13:19:42 georg Exp $

c------------------------------------------------------- basin
	integer nkn,nel
	integer ipv(nkndim)
	integer ipev(neldim)
	integer iarv(neldim)
	integer nen3v(3,neldim)
	real xgv(nkndim), ygv(nkndim)
	real hev(neldim), hkv(nkndim)
	common /nkon/ nkn,nel
	common /ipv/ipv
	common /ipev/ipev
	common /iarv/iarv
	common /nen3v/nen3v
	common /xgv/xgv, /ygv/ygv
	common /hev/hev, /hkv/hkv
c-------------------------------------------------------

