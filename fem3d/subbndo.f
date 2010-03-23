c
c $Id: subbndo.f,v 1.18 2010-03-11 15:36:39 georg Exp $
c
c routines for open boundary conditions
c
c contents :
c
c subroutine bndo_init
c       sets up bndo data structure
c subroutine bndinsert(ib,area,kn)
c       inserts node kn with weight area into list (internal routine)
c
c subroutine bndo_info(iunit)
c       writes info on open boundary nodes to terminal
c
c function is_zeta_bound(k)
c	checks if node k is a zeta boundary
c
c subroutine bndo_setbc(it,what,nlvdim,cv,rbc,uprv,vprv)
c	sets open boundary condition for level boundaries
c subroutine bndo_impbc(it,what,nlvdim,cv,rbc)
c       imposes boundary conditions on open boundary
c subroutine bndo_adjbc(it,nlvdim,cv,uprv,vprv)
c       adjusts boundary conditions on open boundary (values on bnd already set)
c
c subroutine bndo_radiat(it,rzv)
c	imposes radiation condition for levels
c
c notes :
c
c	integer kbcdim
c	integer kopdim
c	parameter ( kbcdim = 100 )	!total number of open boundary nodes
c	parameter ( kopdim = 10 )	!maximum number of nodes close to OB
c
c       integer nbndo                   !total number of OB nodes
c
c	real xynorm(2,kbcdim)		!normal direction for OB node
c
c	integer iopbnd(nkndim)		!if >0 pointer into array irv
c					!if <0 internal boundary (= -ibc)
c
c	integer ibcnod(kbcdim)		!number of boundary
c	integer kbcnod(kbcdim)		!number of boundary node
c	integer itynod(kbcdim)          !type of boundary
c
c	integer nopnod(kbcdim)		!number of internal nodes close to OB
c	integer nopnodes(kopdim,kbcdim)	!nodes close to OB
c
c	real wopnodes(kopdim,kbcdim)	!weights of nodes close to OB
c
c	iopbnd(k) = 0            no open BC
c	iopbnd(k) > 0            external open BC (ibtyp=1,2)
c	iopbnd(k) < 0            internal open BC (ibtyp=3)
c
c the initialization routines should be called only after the
c ... arrays kantv and ieltv have been setup
c
c revision log :
c
c 15.01.2001    ggu     written from scratch
c 03.12.2001    ggu     LEVMX - look out for missing level of near node
c 05.12.2001    ggu     NTBC - BUG -> has not been set before
c 27.03.2003    ggu     in bndo_adjbc use ambient value (bamb)
c 13.03.2004    ggu     in bndo_adjbc only for level BC (LEVELBC)
c 05.10.2004    ggu     new routine bndo_radiat, ibtyp=31
c 31.05.2007    ggu     reset BC for flux to old type (DEBHELP)
c 23.08.2007    ggu     use iopbnd as indicator for ext/int boundary nodes
c 08.04.2008    ggu     file cleaned (new bndo_setbc, bndo_impbc)
c 17.04.2008    ggu     calls to infobnd deleted (subst by get_bnd_ipar)
c 03.09.2008    ggu     new routine bndo_info_file()
c 06.11.2008    ggu     better error handling
c 12.11.2009    ggu     new array itynod and is_zeta_bound()
c
c***********************************************************************

	subroutine bndo_init

c sets up bndo data structure

	implicit none

	include 'subbndo.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw                         

	integer iopbnd(1)
	common /iopbnd/iopbnd
	integer kantv(2,1)
	common /kantv/kantv
	integer nen3v(1)
	common /nen3v/nen3v
        real xgv(1),ygv(1)
        common /xgv/xgv, /ygv/ygv

	logical bexternal
	logical berror
	integer k,nodes,itype
	integer i,ibc
	integer inext,ilast,knext,klast
	integer ie,n,ibase
	integer ii,iii,ib,in,kn,nb,j
	real area
	real dx,dy

	integer nkbnds,itybnd,kbnds,ipext
	real areaele

c----------------------------------------------------------
c set up array iopbnd
c----------------------------------------------------------

	do k=1,nkn
	  iopbnd(k) = 0
	end do

	nbndo = 0
        ndebug = 0              !unit number for debug (in common block)

	do ibc = 1,nbc
	  nodes = nkbnds(ibc)
	  itype = itybnd(ibc)

	  bexternal = ( itype .ge. 1 .and. itype .le. 2 
     +      		    .or. itype .ge. 31 .and. itype .le. 39 )

	  do i=1,nodes
	    k = kbnds(ibc,i)
	    if( bexternal ) then
	      nbndo = nbndo + 1
	      if( nbndo .gt. kbcdim ) goto 99
	      iopbnd(k) = nbndo
	      ibcnod(nbndo) = ibc
	      kbcnod(nbndo) = k
	      itynod(nbndo) = itype
	    else
	      iopbnd(k) = -ibc
	    end if
	  end do

	end do

c----------------------------------------------------------
c set up normal direction
c----------------------------------------------------------

	do i=1,nbndo
	  k = kbcnod(i)
	  ibc = ibcnod(i)

	  knext = kantv(1,k)
	  klast = kantv(2,k)

	  inext = iopbnd(knext)
	  ilast = iopbnd(klast)

c	  -------------------------------
c	  internal consistency check
c	  -------------------------------

	  if( iopbnd(k) .ne. i ) then
	    stop 'internal error bndo (0)'
	  end if
	  if( inext .gt. 0 .and. kbcnod(inext) .ne. knext ) then
	    stop 'internal error bndo (1)'
	  end if
	  if( ilast .gt. 0 .and. kbcnod(ilast) .ne. klast ) then
	    stop 'internal error bndo (2)'
	  end if

c	  -------------------------------
c	  adjacent boundary nodes must be of same boundary
c	  -------------------------------

	  if( inext .gt. 0 .and. ibcnod(inext) .ne. ibc ) then
	    goto 98
	  end if
	  if( ilast .gt. 0 .and. ibcnod(ilast) .ne. ibc ) then
	    goto 98
	  end if

c	  -------------------------------
c	  get normal direction
c	  -------------------------------

	  if( inext .gt. 0 .and. ilast .gt. 0 ) then	!inner node in OB
	    dx = xgv(knext) - xgv(klast)
	    dy = ygv(knext) - ygv(klast)
	  else if( inext .gt. 0 ) then			!first node in OB
	    dx = xgv(knext) - xgv(k)
	    dy = ygv(knext) - ygv(k)
	  else if( ilast .gt. 0 ) then			!last node in OB
	    dx = xgv(k) - xgv(klast)
	    dy = ygv(k) - ygv(klast)
	  else
	    write(6,*) 'One node open boundary not permitted'
	    write(6,*) i,k,ipext(k)
	    stop 'error stop bndo'
	  end if

	  xynorm(1,i) = -dy			!x-component
	  xynorm(2,i) = dx			!y-component
	end do
  
c----------------------------------------------------------
c set up weights and node list
c----------------------------------------------------------

	do i=1,nbndo
	  nopnod(i) = 0
	end do

	do ie=1,nel

	  call elebase(ie,n,ibase)
	  area = areaele(ie)

	  do ii=1,n
	    k = nen3v(ibase+ii)
	    ib = iopbnd(k)
	    if( ib .gt. 0 ) then		!insert inner nodes
	      do iii=1,n
		kn = nen3v(ibase+iii)
	        in = iopbnd(kn)
		if( in .le. 0 ) then		!only inner nodes
		  call bndinsert(ib,area,kn)
		end if
	      end do
	    end if
	  end do

	end do

c----------------------------------------------------------
c scale weights to unit
c----------------------------------------------------------

	berror = .false.

	do i=1,nbndo

	  nb = nopnod(i)
	  if( nb .le. 0 ) then
	    k = kbcnod(i)
	    ibc = ibcnod(i)
	    !write(6,*) i,k,ipext(k)
	    write(6,*) '*** No inner nodes for boundary node'
	    write(6,*) '      boundary = ',ibc,'   node = ',ipext(k)
	    berror = .true.
	  end if

	  area = 0.
	  do j=1,nb
	    area = area + wopnodes(j,i)
	  end do
	  do j=1,nb
	    if( area .gt. 0. ) then
	      wopnodes(j,i) = wopnodes(j,i) / area
	    end if
	  end do
	  
	end do

	if( berror ) stop 'error stop bndo'

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	return
   98	continue
	write(6,*) 'different boundary : ',ibc
	stop 'error stop bndo'
   99	continue
	write(6,*) 'dimension error kbcdim : ',kbcdim
	stop 'error stop bndo'
	end

c***********************************************************************

	subroutine bndinsert(ib,area,kn)

c inserts node kn with weight area into list (internal routine)

	implicit none

	integer ib		!nodal index
	real area		!weight
	integer kn		!node number to insert

	include 'subbndo.h'

	integer nb,j,k

	nb = nopnod(ib)

	do j=1,nb
	  k = nopnodes(j,ib)
	  if( k .eq. kn ) then		!already there -> add weight
	    wopnodes(j,ib) = wopnodes(j,ib) + area
	    return			!done
	  end if
	end do

	nb = nb + 1
	if( nb .gt. kopdim ) then
	  write(6,*) 'Too much inner neighbors: ',nb
	  write(6,*) 'Please raise kopdim in subbndo.h'
	  stop 'error stop bndo'
	end if

	nopnodes(nb,ib) = kn
	wopnodes(nb,ib) = area
	nopnod(ib) = nb

	end

c***********************************************************************

	subroutine bndo_info_file(file)

c writes bndo info to file

	implicit none

	character*(*) file

	if( file .eq. ' ' ) return

	open(1,file=file,status='unknown',form='formatted')
	call bndo_info(1)
	close(1)

	end

c***********************************************************************

	subroutine bndo_info(iunit)

c writes info on open boundary nodes to terminal

	implicit none

        integer iunit

	include 'subbndo.h'

        integer iopbnd(1)
        common /iopbnd/iopbnd

	integer ilhkv(1)
	common /ilhkv/ilhkv

	integer i,k,ibc,nb,j,iu
	integer itybnd,ipext

        iu = iunit
        if( iu .le. 0 ) iu = 6

	write(iu,*) '--------------------------------'
	write(iu,*) 'Information on open boundary nodes'
	write(iu,*) 'Total number of open boundary nodes: ',nbndo

	do i=1,nbndo

	  k = kbcnod(i)
	  ibc = ibcnod(i)
	  nb = nopnod(i)

	  if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

	  write(iu,*) '-------------------------------- bndo_info'
	  write(iu,*) i,k,ipext(k),ibc,itybnd(ibc),ilhkv(k)
	  write(iu,*) nb
	  write(iu,*) (ipext(nopnodes(j,i)),j=1,nb)
	  write(iu,*) (wopnodes(j,i),j=1,nb)
	  write(iu,*) (ilhkv(nopnodes(j,i)),j=1,nb)
	  write(iu,*) (xynorm(j,i),j=1,2)
	end do

	write(iu,*) '--------------------------------'

	end

c***********************************************************************

	function is_zeta_bound(k)

c checks if node k is a zeta boundary

	implicit none

	logical is_zeta_bound
	integer k

	include 'subbndo.h'

        integer iopbnd(1)
        common /iopbnd/iopbnd

	integer ip

	is_zeta_bound = .false.
	ip = iopbnd(k)

	if( ip .gt. 0 ) then
	  if( itynod(ip) .eq. 1 ) then
	    is_zeta_bound = .true.
	  end if
	end if
	
	end

c***********************************************************************

        subroutine bndo_setbc(it,what,nlvdim,cv,rbc,uprv,vprv)

c sets open boundary condition for level boundaries
c
c simply calls bndo_impbc() and bndo_adjbc()

        implicit none

        integer it
        character*(*) what      !conz/temp/salt or else
        integer nlvdim
        real cv(nlvdim,1)
        real rbc(nlvdim,1)	!boundary condition (3D)
	real uprv(nlvdim,1)
	real vprv(nlvdim,1)

c----------------------------------------------------------
c simply imposes whatever is in rbc
c----------------------------------------------------------

        call bndo_impbc(it,what,nlvdim,cv,rbc)

c----------------------------------------------------------
c adjusts for ambient value, no gradient or outgoing flow
c----------------------------------------------------------

	call bndo_adjbc(it,what,nlvdim,cv,uprv,vprv)

c----------------------------------------------------------
c end of routine
c----------------------------------------------------------

	end

c***********************************************************************

        subroutine bndo_impbc(it,what,nlvdi,cv,rbc)

c imposes boundary conditions on open boundary

        implicit none

        integer it
        character*(*) what      !conz/temp/salt or else
        integer nlvdi
        real cv(nlvdi,1)
        real rbc(nlvdi,1)	!boundary condition (3D)

        include 'subbndo.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer iopbnd(1)
        common /iopbnd/iopbnd
        integer ilhkv(1)
        common /ilhkv/ilhkv

        logical bgrad0
        logical bdebug
        integer i,j,k,l
        integer ibc,ibcold
        integer nb,nlev
        integer ibtyp
        real value

        integer ifemopa

        bdebug = .true.
        bdebug = .false.

        do i=1,nbndo

          k = kbcnod(i)
          nb = nopnod(i)
          ibc = ibcnod(i)

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    ibcold = ibc
	  end if

          if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

	  if( ibtyp .eq. 1 ) then
            nlev = ilhkv(k)

            do l=1,nlev
              cv(l,k) = rbc(l,k)
            end do
	  end if

        end do

        end

c***********************************************************************

	subroutine bndo_adjbc(it,what,nlvdi,cv,uprv,vprv)

c adjusts boundary conditions on open boundary (values on bnd already set)
c
c adjusts for ambient value, no gradient or outgoing flow

	implicit none

	integer it
        character*(*) what	!conz/temp/salt or else
	integer nlvdi
	real cv(nlvdi,1)
	real uprv(nlvdi,1)
	real vprv(nlvdi,1)

	include 'subbndo.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw                         

	integer iopbnd(1)
	common /iopbnd/iopbnd
	integer ilhkv(1)
	common /ilhkv/ilhkv

	logical bgrad0
	logical blevel
	logical bdebug
	logical bdggu
	logical bout,bamb
	logical binside
	integer i,j,k,l
	integer ibtyp,igrad0
	integer ibc,ibcold
	integer nb,nlev,ko
        integer ntbc,nlevko
	real dx,dy
	real scal,bc,weight,tweight
	real value

	integer ipext

	integer ifemopa

	bdebug = .true.
	bdebug = .false.
	ibcold = 0

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  write(ndebug,*) 'bndo_adjbc ........... ',what,it
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)

	  if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	    call get_bnd_ipar(ibc,'igrad0',igrad0)
	    bgrad0 = igrad0 .gt. 0
	    blevel = ibtyp .eq. 1
	    ibcold = ibc
	  end if

          if( blevel ) then       !LEVELBC !DEBHELP

	  dx = xynorm(1,i)
	  dy = xynorm(2,i)
	  nlev = ilhkv(k)

	  do l=1,nlev
	    scal = dx * uprv(l,k) + dy * vprv(l,k)
	    bout = scal .le. 0.				!outgoing flow
	    bamb = cv(l,k) .le. -990.			!make ambient value
	    binside = bgrad0 .or. bout .or. bamb
	    if( binside ) then				!take from inside
	      bc = 0.
              ntbc = 0
              tweight = 0.
	      do j=1,nb
		ko = nopnodes(j,i)
                nlevko = ilhkv(ko)
                if( l .le. nlevko ) then        !LEVMX  !only if level exists
                  ntbc = ntbc + 1
		  weight = wopnodes(j,i)
                  tweight = tweight + weight
		  bc = bc + weight * cv(l,ko)
                end if
	      end do
	    else				!impose boundary value
              ntbc = nb                         !NTBC - BUG -> has not been set
	      tweight = 1.			!prob useless
	      bc = cv(l,k)
	    end if

            if( ntbc .gt. 0 ) then              !LEVMX  !at least 1 node found
	      value = bc / tweight
	      if( cv(l,k) .le. -5555. ) then	!differential value
		value = value - 10000. - cv(l,k)
	      end if
	      cv(l,k) = value
            else if( l .gt. 1 ) then            !take from above
              cv(l,k) = cv(l-1,k)
            else
	      write(6,*) i,k,ipext(k),nb,ntbc,ibc,binside
	      write(6,*) nlev,ntbc
              stop 'error stop bndo_adjbc: internal error (5)'
            end if

	  end do

          end if

	end do

	end

c***********************************************************************

	subroutine bndo_radiat(it,rzv)

c imposes radiation condition for levels

	implicit none

	integer it
	real rzv(1)

	include 'subbndo.h'

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw                         

        real znv(1)
        common /znv/znv

	integer iopbnd(1)
	common /iopbnd/iopbnd

	logical bdebug
	integer i,j,k,l
	integer ibc,ibcold
	integer nb,nlev,ko
        integer ntbc,nlevko
	integer ibtyp
	real bc,weight,tweight

	integer ipext

	integer ifemopa

	bdebug = .true.
	bdebug = .false.
	ibcold = 0

	if( bdebug ) then
	  if( ndebug .eq. 0 ) then
	    ndebug = ifemopa('bndo_adjbc (91)','.bndo','form','unknown')
            call bndo_info(ndebug)
	  end if
	  write(ndebug,*) 'bndo_adjbc ........... ','radiat',it
	end if

	do i=1,nbndo

	  k = kbcnod(i)
	  nb = nopnod(i)
	  ibc = ibcnod(i)

	  if( ibc .ne. ibcold ) then
	    call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	  end if
          !write(78,*) i,k,nb,ibc,ibtyp

	  if( iopbnd(k) .ne. i ) stop 'internal error bndo: (11)'

          if( ibtyp .eq. 31 ) then       !radiation condition only

	    bc = 0.
            ntbc = 0
            tweight = 0.
	    do j=1,nb
		ko = nopnodes(j,i)
                ntbc = ntbc + 1
		weight = wopnodes(j,i)
                tweight = tweight + weight
		bc = bc + weight * znv(ko)
	    end do

            if( ntbc .gt. 0 ) then         !LEVMX  !at least 1 node found
	      rzv(k) = bc / tweight
            else
              stop 'error stop bndo_radiat: internal error (5)'
            end if

          end if

	end do

	end

c***********************************************************************

