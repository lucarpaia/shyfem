c
c $Id: adj5el.f,v 1.5 2009-01-26 13:27:24 georg Exp $
c
c description :
c
c 5 grade routines
c
c contents :
c
c subroutine elim5(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c			eliminates low grades
c subroutine elim55(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c			eliminates 5-5 connections
c
c***********************************************************

	subroutine elim5(nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)

c eliminates low grades

	implicit none

        integer nkn,nel,ngrdim
        integer ngrade(1)
        integer nbound(1)
        integer ngri(2*ngrdim,1)
        integer nen3v(3,1)

        integer k,n

        write(6,*) 'eliminating grades for grade 5... '

        do k=1,nkn
          n = ngrade(k)
          if( n .eq. 5 .and. nbound(k) .eq. 0 ) then
            call elim55(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)
c	    call nodeinfo(1290)
	    call chkgrd
          end if
        end do

	end

c***********************************************************

	subroutine elim55(k,nkn,nel,ngrdim,ngrade,nbound,ngri,nen3v)

c eliminates 5-5 connections

	implicit none

	integer k
        integer nkn,nel,ngrdim
        integer ngrade(1)
        integer nbound(1)
        integer ngri(2*ngrdim,1)
        integer nen3v(3,1)

	logical bdebug
        integer n,i,nc,nmax,ii
	integer ie1,ie2
	integer ip1,ip2
	integer np,nt,nn
	integer nval,ip
	integer nga(0:30)
	integer ngr(0:30)
	integer nba(0:30)

        real xgv(1), ygv(1)
        common /xgv/xgv, /ygv/ygv

	integer ifindel

	if( k .gt. nkn ) return

	bdebug = .true.
	bdebug = .false.
	if( k .eq. 1138 ) bdebug = .true.

	if( bdebug ) then
	  write(6,*) '==============================================='
	  write(6,*) 'debug of new node: ',k
	end if

c make circular list
c
c nga 	node numbers around k
c ngr	grades of node numbers around k
c nba	boundary  flag for nodes around k

        n = ngrade(k)
	nga(0) = ngri(n,k)
	do i=1,n
	  nga(i) = ngri(i,k)
	end do
	nga(n+1) = ngri(1,k)

	do i=0,n+1
	  ngr(i) = ngrade(nga(i))
	  nba(i) = 0
	  if( nbound(nga(i)) .ne. 0 ) then
	    ngr(i) = 6	!FIXME
	    nba(i) = 1
	  end if
	end do

c check if exchange is possible

	nc = 0		!how many of this nmax value
	nmax = 0	!maximum sum of grades -> must be at least 3
	ip = 0		!pointer to node in list that has been chosen
	do i=1,n
	  np = ngr(i-1)
	  nt = ngr(i)
	  nn = ngr(i+1)

	  nval = np + nn - n - nt

	  if( nval .gt. nmax ) then
	    nc = 1
	    ip = i
	    nmax = nval
	  else if( nval .eq. nmax ) then
	    nc = nc + 1
	  end if
	end do

	if( nmax .lt. 3 ) return

	write(6,*) k,n,nmax,nc,ip

c nc gives number of occurences of this value of nmax ...
c ip is the pointer to the node to be exchanged
c
c we decide to take the first choice
c
c k is eliminated, nga(ip) is retained (to account for boundary node)

	if( bdebug ) then
	    write(6,*) 'exchanging with node ... ',nga(ip)
	    write(6,'(7i10)') (nga(i),i=0,n+1)
	    write(6,'(7i10)') (ngr(i),i=0,n+1)
	    write(6,'(7i10)') (nba(i),i=0,n+1)
	    call plosno(k)
	    call plosno(nga(ip))
	end if

c find elements that have to be deleted

	ie1 = ifindel(k,nga(ip),nga(ip+1))
	ie2 = ifindel(k,nga(ip-1),nga(ip))

	if( ie1 .eq. 0 .or. ie2 .eq. 0 ) then
	  stop 'error stop elim55: internal error (2)'
	end if

	if( bdebug ) then
	  write(6,*) 'elements to be deleted... ',ie1,ie2
	  write(6,*) ie1,k,nga(ip),nga(ip+1)
	  write(6,*) (nen3v(ii,ie1),ii=1,3)
	  write(6,*) ie2,k,nga(ip-1),nga(ip)
	  write(6,*) (nen3v(ii,ie2),ii=1,3)
	  call plosel2(ie1,ie2,nkn,nel,ngrdim,nen3v,ngrade,ngri)
	end if

c delete elements

	if( ie1 .gt. ie2 ) then		!to avoid bug
	  call delele(ie1,nkn,nel,ngrdim,ngrade,ngri)
	  call delele(ie2,nkn,nel,ngrdim,ngrade,ngri)
	else
	  call delele(ie2,nkn,nel,ngrdim,ngrade,ngri)
	  call delele(ie1,nkn,nel,ngrdim,ngrade,ngri)
	end if

	if( bdebug ) then
	  write(6,*) 'grade index befor manipulation:'
	  call prgr(k,ngrdim,ngrade,ngri)
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	  call prgr(nga(ip-1),ngrdim,ngrade,ngri)
	  call prgr(nga(ip+1),ngrdim,ngrade,ngri)
	end if

c new coordinates for node

	if( nba(ip) .le. 0 ) then	!no boundary node
	  xgv(k) = 0.5 * ( xgv(k) + xgv(nga(ip)) )
	  ygv(k) = 0.5 * ( ygv(k) + ygv(nga(ip)) )
	end if

c substitute all occurrences of k with nga(ip)

	call subnod(k,nga(ip),nkn,nel,ngrdim,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'after substitution...'
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	  call prgr(nga(ip-1),ngrdim,ngrade,ngri)
	  call prgr(nga(ip+1),ngrdim,ngrade,ngri)
	end if

c adjourn grade (delete) for nodes ip-1, ip+1

	call delgr(nga(ip-1),nga(ip),ngrdim,ngrade,ngri)
	call delgr(nga(ip+1),nga(ip),ngrdim,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'after deleting ip-1,ip+1...'
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	  call prgr(nga(ip-1),ngrdim,ngrade,ngri)
	  call prgr(nga(ip+1),ngrdim,ngrade,ngri)
	end if

c adjourn grade index for ip and delete node k finally

	call delgr(nga(ip),nga(ip),ngrdim,ngrade,ngri)
	call delnod(k,nkn,nel,ngrdim,ngrade,ngri)
	call subval(n+2,nga(0),nkn+1,k)	!if nkn is in nga

	if( bdebug ) then
	  write(6,*) 'after deleting ip...'
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	end if

	ip1 = mod(ip+2,n)
	ip2 = mod(ip+3,n)
	call insgr(nga(ip),nga(ip+1),nga(ip1),ngrdim,ngrade,ngri)
	call insgr(nga(ip),nga(ip1),nga(ip2),ngrdim,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'grade index after manipulation:'
	  call prgr(nga(ip),ngrdim,ngrade,ngri)
	  call prgr(nga(ip-1),ngrdim,ngrade,ngri)
	  call prgr(nga(ip+1),ngrdim,ngrade,ngri)
	end if

	if( bdebug ) then
	  call plosno(nga(ip))
	  call plosno(nga(ip-1))
	  call plosno(nga(ip+1))
	end if

	if( bdebug ) then
	  write(6,*) 'end of debug of node ',k
	  write(6,*) '==============================================='
	end if

	end

