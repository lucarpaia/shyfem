c
c $Id: laplap.f,v 1.9 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 20.08.2003	ggu	new laplacian interpolation
c 02.09.2003	ggu	some comments, write to .dat file
c 30.10.2003	ggu	subroutine prepare_bc included in this file
c 04.03.2004	ggu	writes also number of variables (1)
c 11.03.2009	ggu	bug fix -> declare hev() here
c
c notes :
c
c please prepare file like this:
c
c----------------- start
c
c k1	val1
c k2	val2
c ...
c kn	valn
c----------------- end
c
c first line of file must be empty !!!
c
c run memory and set the basin ( memory -b venlag62 )
c run laplap with input file ( laplap < input.dat )
c
c****************************************************************

        program laplap

c laplacian interpolation

	implicit none

	include 'param.h'

	integer matdim
	parameter (matdim = nkndim*100)

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

        real xgv(nkndim), ygv(nkndim)
        real hm3v(3,neldim)
        common /xgv/xgv, /ygv/ygv
        common /hm3v/hm3v

        integer nen3v(3,neldim)
        integer ipev(neldim), ipv(nkndim)
        integer iarv(neldim)
        common /nen3v/nen3v
        common /ipev/ipev, /ipv/ipv
        common /iarv/iarv

        real hev(neldim)
        common /hev/hev

	include 'evmain.h'

	real rmat(matdim)
	real zv(nkndim)
	real rzv(nkndim)

	integer k,ie
        integer ilev,ivar
	real flag
	real zmin,zmax

	integer iapini

	flag = 1.23456e+23

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn = ',nkn,'  nel = ',nel
        write(6,*) ' mbw = ',mbw,'  ngr = ',ngr
        write(6,*)
        write(6,*) ' dcor = ',dcor,'  dirn = ',dirn
        write(6,*)

c-----------------------------------------------------------------
c check dimension
c-----------------------------------------------------------------

	if( nkn*(mbw+1) .gt. matdim ) then
	  write(6,*) nkn*(mbw+1),matdim
	  stop 'error stop laplap: matdim too small'
	end if

c-----------------------------------------------------------------
c set up ev
c-----------------------------------------------------------------

	call set_ev
	call check_ev

	do ie=1,nel
	  hev(ie) = 1.
	end do

c-----------------------------------------------------------------
c read BC and interpolate
c-----------------------------------------------------------------

	call prepare_bc(' ',rzv,nkn,flag)

	call lapint(rmat,zv,rzv,flag)

c-----------------------------------------------------------------
c min/max of interpolated values
c-----------------------------------------------------------------

	call mima(zv,nkn,zmin,zmax)
	write(6,*) 'min/max: ',zmin,zmax

c-----------------------------------------------------------------
c write to NOS file laplace.nos
c-----------------------------------------------------------------

	call wrnos2d('laplace','laplace interpolation',zv)

c-----------------------------------------------------------------
c write to DAT file laplace.dat
c-----------------------------------------------------------------

        ilev = 0
        ivar = 1

	open(1,file='laplace.dat',status='unknown',form='unformatted')
	write(1) nkn,ilev,ivar
	write(1) (zv(k),k=1,nkn)
	close(1)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c****************************************************************

	subroutine prepare_bc(file,rzv,nkn,flag)

c reads boundary conditions from file and sets up array
c
c file must be made like this:
c
c	k1, val1
c	k2, val2
c	...
c	kn, valn

	implicit none

	character*(*) file
	real rzv(1)
	integer nkn
	real flag

	integer iunit
	integer k,kn
	real val

	integer ipint

	if( file .ne. ' ' ) then
	  open(1,file=file,status='old',form='formatted',err=97)
	  iunit = 1
	else
	  iunit = 5
	end if

	do k=1,nkn
	  rzv(k) = flag
	end do

	write(6,*) '...reading boundary conditions from unit :',iunit
	write(6,*) '   format: k  val'

    1	continue
	  read(iunit,*,end=2) k,val
	  kn = ipint(k)
	  if( kn .le. 0 ) goto 99
	  if( kn .gt. nkn ) goto 98
	  rzv(kn) = val
	  write(6,*) k,kn,val
	  goto 1
    2	continue

	if( iunit .ne. 5 ) close(iunit)

	return
   97	continue
	write(6,*) file
	stop 'error stop prepare_bc: cannot open file'
   98	continue
	write(6,*) k,kn,nkn,val
	stop 'error stop prepare_bc: error in internal node number'
   99	continue
	write(6,*) k,kn,nkn,val
	stop 'error stop prepare_bc: no such node'
	end

c******************************************************************

