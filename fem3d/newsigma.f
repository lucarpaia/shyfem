c
c $Id: newsigma.f,v 1.3 2009-03-24 17:49:19 georg Exp $
c
c routines for handling initialization from sigma levels
c
c contents :
c
c revision log :
c
c 07.11.2008    ggu     program written from scratch
c 06.12.2008    ggu     bsigma set from STR (nsigma)
c 24.03.2009    ggu     use nsigma to check if initialized
c
c*****************************************************************

	subroutine handle_sigma_init

c initializes u/v/z from sigma level data

	implicit none

	include 'param.h'
	include 'sigma.h'

	character*(60) sigma_l,sigma_d,sigma_u,sigma_v,sigma_z,sigma_b
	logical bsigma
        integer nsig
        real getpar

        nsig = nint(getpar('nsigma'))

	bsigma = .false.		!if true read in sigma data
	bsigma = .true.			!if true read in sigma data
        bsigma = nsig .gt. 0

	nsigma = 0			!number of sigma levels read (init)

	if( .not. bsigma ) return	!no initialization with sigma data

	sigma_l = 'RosSigmas.txt'
	sigma_d = 'PtsBathyAll.txt'
	sigma_b = 'PtsBathyBoundary.txt'

	sigma_u = 'Uvelocities_init.txt'
	sigma_v = 'Vvelocities_init.txt'
	sigma_z = 'elevation_init.txt'

	call set_sigma_levels(sigma_l)
	call set_sigma_depths(sigma_d)
	call init_vel_from_sigma(sigma_u,sigma_v,sigma_z)

	end

c*****************************************************************

	subroutine set_sigma_levels(file)

c read in sigma levels

	implicit none

	character*(*) file

	include 'param.h'
	include 'sigma.h'

	integer i

	nsigma = 0

	open(1,file=file)
	read(1,*) nsigma
	nsigma = nsigma - 1		!number of layers, not levels

	if( nsigma .gt. nsidim ) goto 98

	do i=0,nsigma
	  read(1,*) siglev(i)
	end do

	close(1)

	if( siglev(0) .ne. 0. ) goto 99
	if( siglev(nsigma) .ne. 1. ) goto 99

	write(6,*) nsigma,' sigma layers read from file ',file

	return
   98	continue
	write(6,*) 'too many sigma levels'
	write(6,*) 'nsigma = ',nsigma,'    nsidim = ',nsidim
	stop 'error stop set_sigma_levels: nsidim'
   99	continue
	write(6,*) 'problems in sigma levels...'
	write(6,*) nsigma
	write(6,*) (siglev(i),i=0,nsigma)
	stop 'error stop set_sigma_levels: sigma'
	end

c*****************************************************************

	subroutine set_sigma_depths(file)

c read in depth file for sigma layers

	implicit none

	include 'param.h'
	include 'sigma.h'

	character*(*) file

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer k,l
	integer kint,kext
	integer nknaux
	real x,y
	real hfem,hfd

	do k=1,nkn
	  sigdep(k) = -1.
	end do

	open(1,file=file)
	read(1,*) nknaux

	if( nknaux .ne. nkn ) goto 98

	do k=1,nkn
	  !read(1,*) kint,kext,l,x,y,hfem,hfd
	  read(1,*) kext,x,y,hfem,hfd
	  sigdep(k) = hfd
	end do

	close(1)

	write(6,*) 'sigma layer depth read from file ',file

	return
   98	continue
	write(6,*) 'not compatible number of total nodes'
	write(6,*) 'nknaux = ',nknaux,'    nkn = ',nkn
	stop 'error stop set_sigma_depths: nkn'
	end

c*****************************************************************

	subroutine init_vel_from_sigma(fileu,filev,filez)

c initializes u/v/z from sigma level data

	implicit none

	include 'param.h'
	include 'sigma.h'

	character*(*) fileu,filev,filez

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v
        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1)
        common /znv/znv
        real zenv(3,1)
        common /zenv/zenv
        real ulnv(nlvdim,1),vlnv(nlvdim,1)
        common /ulnv/ulnv, /vlnv/vlnv

	real sigu(nsidim,nkndim)
	real sigv(nsidim,nkndim)
	real sigz(nkndim)

	real femu(nlvdim)
	real femv(nlvdim)

	integer k,l,lmax,ie,ii,ifact
	integer kint,kext
	integer nknaux
	integer itime
	real x,y,h

	integer ipext

	if( nsigma .gt. nsidim ) goto 99

c----------------------------------------------------------
c read in water levels
c----------------------------------------------------------

	write(6,*) 'start reading water levels from file ',filez
	open(1,file=filez)
	read(1,*) itime,nknaux
	if( nknaux .ne. nkn ) goto 98
	do k=1,nkn
	  !read(1,*) kint,kext,lmax,x,y,sigz(k)
	  read(1,*) kext,x,y,h,lmax,sigz(k)
	  if( ipext(k) .ne. kext ) goto 97
	end do
	close(1)

c----------------------------------------------------------
c read in velocities on nodes and sigma levels
c----------------------------------------------------------

	write(6,*) 'start reading u velocities from file ',fileu
	open(1,file=fileu)
	read(1,*) itime,nknaux
	if( nknaux .ne. nkn ) goto 98
	do k=1,nkn
	  !read(1,*) kint,kext,lmax,x,y,(sigu(l,k),l=1,lmax)
	  read(1,*) kext,x,y,h,lmax,(sigu(l,k),l=1,lmax)
	  if( ipext(k) .ne. kext ) goto 97
	end do
	close(1)

	write(6,*) 'start reading v velocities from file ',filev
	open(1,file=filev)
	read(1,*) itime,nknaux
	if( nknaux .ne. nkn ) goto 98
	do k=1,nkn
	  !read(1,*) kint,kext,lmax,x,y,(sigv(l,k),l=1,lmax)
	  read(1,*) kext,x,y,h,lmax,(sigv(l,k),l=1,lmax)
	  if( ipext(k) .ne. kext ) goto 97
	end do
	close(1)

	write(6,*) 'finished reading initial fields'

c----------------------------------------------------------
c initialize water levels and adjust depth
c----------------------------------------------------------

	do k=1,nkn
	  znv(k) = sigz(k)
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    zenv(ii,ie) = znv(k)
	  end do
	end do

c----------------------------------------------------------
c loop over elements and do interpolation
c----------------------------------------------------------

	ifact = 0		!do not compute fact in interpolation routine

	do ie=1,nel

	  lmax = ilhv(ie)

	  do l=1,lmax
	    femu(l) = 0.
	    femv(l) = 0.
	  end do

	  do ii=1,3
	    k = nen3v(ii,ie)
	    call add_sigma_to_element(ie,k,ifact,nsidim
     +                          ,sigu,sigv,femu,femv)
	  end do

	  do l=1,lmax
	    ulnv(l,ie) = femu(l) / 3.
	    vlnv(l,ie) = femv(l) / 3.
	  end do

	end do

c----------------------------------------------------------
c set all other variables
c----------------------------------------------------------

	call make_new_depth	!makes layer thickness
	call vtot		!makes utlnv, vtlnv
	call init_uv		!here we need zenv, utlnv, vtlnv, hdenv

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   97	continue
	write(6,*) 'incompatible nodes'
	write(6,*) k,ipext(k),kext
	stop 'error stop init_vel_from_sigma: kext'
   98	continue
	write(6,*) 'not compatible number of total nodes'
	write(6,*) 'nknaux = ',nknaux,'    nkn = ',nkn
	stop 'error stop init_vel_from_sigma: nkn'
   99	continue
	write(6,*) 'too many sigma levels'
	write(6,*) 'nsigma = ',nsigma,'    nsidim = ',nsidim
	stop 'error stop init_vel_from_sigma: nsidim'
	end

c*****************************************************************

	subroutine add_sigma_to_element(ie,k,ifact,ndim,uval,vval,u,v)

	implicit none

	include 'param.h'
	include 'sigma.h'

	integer ie,k,ifact,ndim
	real uval(ndim,nkndim)
	real vval(ndim,nkndim)
	real u(1),v(1)

        integer ilhv(1)
        common /ilhv/ilhv
        real znv(1),hev(1)
        common /znv/znv, /hev/hev
        real hlv(1)
	common /hlv/hlv

	real hlfem(0:nlvdim+1)
	real femuval(nlvdim)
	real femvval(nlvdim)

	real hsig(0:nsidim+1)
	real siguval(nsidim+1)
	real sigvval(nsidim+1)

	integer l,lmax
	real hfem,hfd,zfem
	real fact

	if( nsigma .le. 0 ) goto 98

	hfem = hev(ie)
	lmax = ilhv(ie)

	hfd = sigdep(k)
	zfem = znv(k)		!this might not be set at the boundary

	hlfem(0) = -zfem
	do l=1,lmax-1
	  hlfem(l) = hlv(l)
	end do
	hlfem(lmax) = hfem

	fact = 1.
	if( ifact .ne. 0 ) fact = (hfd+zfem) / (hfem+zfem)

	call make_sigma_layer_bottom(nsigma,siglev,hfem,zfem,hsig)

	do l=1,nsigma
	  siguval(l) = uval(l,k)
	  sigvval(l) = vval(l,k)
	end do

	call intp_vert(nsigma,hsig,siguval,lmax,hlfem,femuval)
	call intp_vert(nsigma,hsig,sigvval,lmax,hlfem,femvval)

	do l=1,lmax
	  u(l) = u(l) + fact * femuval(l)
	  v(l) = v(l) + fact * femvval(l)
	end do

	return
   98	continue
	stop 'error stop add_sigma_to_element: sigma layers not init'
	end

c*****************************************************************

	subroutine make_sigma_layer_thickness(nsigma,siglev,h,z,hlayer)

c compute layer thickness for depth h and water level z

	implicit none

	integer nsigma
	real siglev(0:nsigma)
	real h,z
	real hlayer(nsigma)

	integer l

	do l=1,nsigma
	  hlayer(l) = (h+z) * (siglev(l) - siglev(l-1))
	end do

	end

c*****************************************************************

	subroutine make_sigma_layer_bottom(nsigma,siglev,h,z,hbottom)

c compute bottom of layers for depth h and water level z

	implicit none

	integer nsigma
	real siglev(0:nsigma)
	real h,z
	real hbottom(0:nsigma)

	integer l

	do l=0,nsigma
	  hbottom(l) = - z - (h+z) * siglev(l)
	end do

	end

c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************
c*****************************************************************

	subroutine intp_vert(nl1,zb1,var1,nl2,zb2,var2)

c vertical interpolation of variables from one grid to another

	implicit none

	integer nl1		!number of levels of first grid
	real zb1(0:nl1+1)	!depth of bottom of vertical boxes
	real var1(nl1+1)	!value of variable at center of box
	integer nl2		!number of levels of second grid
	real zb2(0:nl2)		!depth of bottom of vertical boxes
	real var2(nl2)		!value of variable at center of box

c values are interpolated from first grid to second grid
c z levels refer to bottom of each grid (zb1(nl1) is total depth of column)
c zb(0) is surface ... normally 0
c variables are considered at center of box
c values are considered to be constant for every box
c
c ATTENTION: arrays zb1, var1 MUST be dimensioned one element bigger
c	then the available data. The last element is altered by this
c	subroutine, but should be of no concern for the calling program
c
c output is var1, all other variables are input values

	logical bmiss
	integer l,j,ltop1,lbot1
	real ztop2,ztop1,zbot2,zbot1
	real ztop,zbot
	real vint1,vint2,fact
	real val
	logical bcons

	zb1(nl1+1) = zb2(nl2)
	var1(nl1+1) = var1(nl1)
	bcons = .true.			!conserve total quantity
	bcons = .false.			!do not conserve total quantity

	ltop1 = 0
	vint2 = 0.

	do l=1,nl2
	  ztop2 = zb2(l-1)
	  zbot2 = zb2(l)

	  do while( ltop1 .lt. nl1 .and. zb1(ltop1+1) .le. ztop2 )
	    !write(6,*) 'adjourning top depth: ',ltop1,zb1(ltop1+1),ztop2
	    ltop1 = ltop1 + 1
	  end do

	  bmiss = .false.
	  do lbot1=ltop1+1,nl1+1
	    if( zb1(lbot1) .ge. zbot2 ) goto 1
	  end do
	  lbot1 = nl1
	  bmiss = .true.
    1	  continue

	  ztop1 = zb1(ltop1)
	  zbot1 = zb1(lbot1)

	  if( ztop1 .gt. ztop2 .or. zbot1 .lt. zbot2 ) goto 99
	  if( bmiss ) goto 98

	  !write(6,*) l,ltop1,lbot1,ztop2,zbot2,zb1(lbot1)

	  val = 0.
	  do j=ltop1+1,lbot1
	      ztop = max(zb1(j-1),ztop2)
	      zbot = min(zb1(j),zbot2)
	      val = val + var1(j) * ( zbot - ztop )
	  end do

	  vint2 = vint2 + val			!integrated value
	  var2(l) = val / (zbot2-ztop2)
	  ltop1 = lbot1 - 1

	end do

	zb1(nl1+1) = 0.
	var1(nl1+1) = 0.

	if( bcons ) then	!must conserve total content of scalar
	  vint1 = 0.
	  do l=1,nl1
	    ztop = zb1(l-1)
	    zbot = zb1(l)
	    vint1 = vint1 + var1(l) * ( zbot - ztop )
	  end do

	  if( vint1 .eq. vint2 ) return

	  fact = vint1 / vint2
	  do l=1,nl2
	    var2(l) = fact * var2(l)
	  end do
	end if

	return
   98	continue
	stop 'error stop intp_vert: missing value not possible'
   99	continue
	write(6,*) ztop1,ztop2
	write(6,*) zbot1,zbot2
	stop 'error stop intp_vert: interval in (1) must include (2)'
	end

c*****************************************************************

