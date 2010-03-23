c
c $Id: lagrange_inout.f,v 1.2 2009-03-11 16:25:59 georg Exp $
c
c handles insertion and output of particles
c
c revision log :
c
c 05.02.2009    ggu	copied from lagrange_cont.f and integrated from others
c
c*******************************************************************

	subroutine insert_particle(ie,x,y)

c inserts particle at position x,y - if element unknown set ie=0

	implicit none

	include 'param.h'
	include 'lagrange.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer ie
	real x,y

	nbdy = nbdy + 1

	if( nbdy .gt. nbdydim ) goto 99

        xst(nbdy)=x
        yst(nbdy)=y
        est(nbdy)=ie

        x_body(nbdy)=x
        y_body(nbdy)=y
        ie_body(nbdy)=ie

        tin(nbdy)=it

	zst(nbdy) = 0.5
	z_body(nbdy) = 0.5

	return
   99	continue
	write(6,*) 'nbdydim = ',nbdydim
	write(6,*) 'Cannot insert more than nbdydim particles'
	write(6,*) 'Please change nbdydim in param.h and recompile'
	stop 'error stop insert_particle: nbdy'
	end

c*******************************************************************

	subroutine copy_particle(ifrom,ito)

c copies particle from ifrom to ito

	implicit none

	include 'param.h'
	include 'lagrange.h'

	integer ifrom,ito

	if( ifrom .eq. ito ) return

        xst(ito) = xst(ifrom)
        yst(ito) = yst(ifrom)
        est(ito) = est(ifrom)

        x_body(ito)  = x_body(ifrom)
        y_body(ito)  = y_body(ifrom)
        ie_body(ito) = ie_body(ifrom)

        tin(ito) = tin(ifrom)
	
	end

c*******************************************************************

	subroutine delete_particle(ip)

c deletes particle ip

	implicit none

	include 'param.h'
	include 'lagrange.h'

	integer ip

	if( ie_body(ip) .gt. 0 ) ie_body(ip) = -ie_body(ip)

	end

c*******************************************************************

	subroutine ntot_particles(ntot)

c returns total number of particles

	implicit none

	include 'param.h'
	include 'lagrange.h'

	integer ntot

	ntot = nbdy

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine lgr_output(iu,it)

c outputs particles to file
c
c if body has exited the element number is negativ (last element)
c once it has been written to output with negative ie, it is set to 0

	implicit none

        include 'param.h'
        include 'lagrange.h'

	integer iu,it

	integer nvers,mtype
	parameter ( nvers = 3 , mtype = 367265 )

	integer nn,nout,ie,i
	real x,y,z
	real getpar

	integer icall
	save icall
	data icall / 0 /

c----------------------------------------------------------------
c initialize params
c----------------------------------------------------------------

	if( icall .eq. 0 ) then
	  write(iu) mtype,nvers
	  write(6,*) 'lgr file initialized...'
	end if

	icall = icall + 1

c----------------------------------------------------------------
c count active particles
c----------------------------------------------------------------

	nn = 0
	nout = 0
	do i=1,nbdy
	  ie = ie_body(i)
	  if( ie .ne. 0 ) nn = nn + 1
	  if( ie .lt. 0 ) nout = nout + 1
	end do
	  
c----------------------------------------------------------------
c write to file
c----------------------------------------------------------------

	write(iu) it,nbdy,nn,nout

	do i=1,nbdy
          x = x_body(i)
          y = y_body(i)
          z = z_body(i)
	  ie = ie_body(i)
	  if( ie .ne. 0 ) then
	    write(iu) i,x,y,z,ie,xst(i),yst(i),zst(i),est(i),tin(i)
	  end if
	  if( ie .lt. 0 ) ie_body(i) = 0		!flag as out
	end do

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	end

c*******************************************************************

        subroutine lgr_output_concentrations

c outputs particles as density (concentration) to NOS file

        implicit none

        include 'param.h'
        include 'lagrange.h'
        include 'ev.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer nen3v(3,1)
	common /nen3v/nen3v

	integer ie,ii,k
	integer ic,i
	integer nvar,ivar,nlvdi
	real area_el,z

        integer ecount(neldim)
        integer kcount(nkndim)
        real area(nkndim)
        real density(nkndim)

	integer ifem_open_file

        integer iunit
        save iunit
        data iunit / 0 /

c---------------------------------------------------------
c initialize
c---------------------------------------------------------

	do ie=1,nel
	  ecount(ie) = 0
	end do
	do k=1,nkn
	  kcount(k) = 0
	  area(k) = 0.
	end do

c---------------------------------------------------------
c count particles in elements
c---------------------------------------------------------

	do i=1,nbdy
          z = z_body(i)			! z == 1 => bottom
	  ie = ie_body(i)
	  if( ie .ne. 0 ) then
	    ecount(ie) = ecount(ie) + 1
	  end if
	end do

c---------------------------------------------------------
c copy to nodes
c---------------------------------------------------------

	do ie=1,nel
	  ic = ecount(ie)
	  area_el = ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    kcount(k) = kcount(k) + ic
	    area(k) = area(k) + area_el
	  end do
	end do

c---------------------------------------------------------
c compute density
c---------------------------------------------------------

	do k=1,nkn
	  density(k) = kcount(k) / area(k)
	end do

c---------------------------------------------------------
c write to file
c---------------------------------------------------------

	nvar = 1
	ivar = 131
	nlvdi = 1
	call conwrite(iunit,'lgc',nvar,ivar,nlvdi,density)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

        end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine print_elements(text)

c writes element numbers of particles to terminal

	implicit none

	include 'param.h'
	include 'lagrange.h'

	character*(*) text

	integer i,ii

	write(6,*) text
	write(6,*) 'particles: ',nbdy
	do i=1,nbdy,10
	  write(6,*) (ie_body(ii),ii=i,min(i+9,nbdy))
	end do

	end

c*******************************************************************
	
	subroutine compress_particles

c deletes particles not in system and compresses array

	implicit none

	include 'param.h'
	include 'lagrange.h'

        integer itanf,itend,idt,nits,niter,it
        common /femtim/ itanf,itend,idt,nits,niter,it

	integer ifree,icopy,i
	integer nbefore

	logical is_free,is_particle
	is_free(i) = ie_body(i) .le. 0
	is_particle(i) = ie_body(i) .gt. 0

	ifree = 0
	icopy = nbdy+1
	nbefore = nbdy

	do while( ifree .lt. icopy )

	  do i=ifree+1,icopy-1
	    if( is_free(i) ) goto 1
	  end do
    1	  continue
	  ifree = i
	  if( ifree .eq. icopy ) goto 3

	  do i=icopy-1,ifree+1,-1
	    if( is_particle(i) ) goto 2
	  end do
    2	  continue
	  icopy = i
	  if( ifree .eq. icopy ) goto 3

	  call copy_particle(icopy,ifree)
	  call delete_particle(icopy)

	end do
    3	continue

	nbdy = ifree - 1

c--------------------------------------------------------------
c check particle structure
c--------------------------------------------------------------

	write(6,*) 'compress_particles: ',nbefore,nbdy,nbefore-nbdy

	do i=1,nbdy
	  if( is_free(i) ) then
	    write(6,*) 'particle i = ',i,' is deleted'
	    stop 'error stop compress_particles: structure'
	  end if
	end do

	do i=nbdy+1,nbefore
	  if( is_particle(i) ) then
	    write(6,*) 'particle i = ',i,' is not deleted'
	    stop 'error stop compress_particles: structure'
	  end if
	end do

	end

c*******************************************************************

	subroutine test_compress

c tests compress routine

	implicit none

	integer nadd,nloop,i,ndel
	integer iloop,ntot,ip,ie
	real ggrand

	nadd = 3
	nloop = 20
	ie = 0

	do iloop=1,nloop

	do i=1,nadd
	  ie = ie + 1
	  call insert_particle(ie,1.,1.)
	end do

	call ntot_particles(ntot)

	ndel = nint(ggrand(99)*ntot)

	do i=1,ndel
	  ip = nint(ggrand(99)*ntot)
	  call delete_particle(ip)
	end do

	call print_elements('before compress')
	call compress_particles
	call print_elements('after compress')

	end do

	end

c*******************************************************************
c	call test_compress
c	end
c*******************************************************************

