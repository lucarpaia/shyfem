
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2019  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! compares two debug files or writes info on one

! revision log :
!
! 03.06.2020	ggu	newly written based on check_debug.f
! 10.11.2021    ggu     avoid warning for stack size
! 30.03.2022    ggu     ntime was not initialized
! 02.04.2022    ggu     new routine i_info()

!**************************************************************************

	program check_shympi_debug

	implicit none

	integer nc,ierr

	nc = command_argument_count()

	if( nc == 1 ) then
	  call read_file
	else if( nc == 2 ) then
	  call compare_files(ierr)
	else
	  write(6,*) 'nc = ',nc
	  stop 'error stop check_shympi_debug: need one or two files'
	end if

	if( ierr > 0 ) then
	  if( ierr == 99 ) ierr = 100	!terrible hack - FIXME
	  call exit(ierr)
	else
	  call exit(99)
	end if

	end

!**************************************************************************

	subroutine read_file

	implicit none

	integer nc
	integer ntime,nrec
	integer nh,nv,nt
	integer ios
	double precision dtime
	character*60 name_one,text

	nc = command_argument_count()
	if( nc .lt. 1 ) then
	  write(6,*) 'Usage: check_debug file'
	  stop 'error stop check_debug: no files given'
	end if

	call get_command_argument(1,name_one)

	open(1,file=name_one,status='old',form='unformatted',iostat=ios)

	if( ios /= 0 ) stop 'error opening file'

	write(6,*) 'file 1: ',trim(name_one)

	ntime = 0

	do while(.true.)

	  read(1,end=9) dtime
	  write(6,*) 'time = ',dtime
	  ntime = ntime + 1

	  write(6,*) '       irec          nh          nv' //
     +			'          nt name'

	  nrec = 0
	  do while(.true.)
	    read(1,end=9) nh,nv,nt
	    if( nh == 0 ) exit
	    nrec = nrec + 1
	    read(1,end=9) text
	    read(1)
	    write(6,*) nrec,nh,nv,nt,trim(text)
	  end do
	end do

    9	continue
	write(6,*) 'total time records read: ',ntime

	end

!**************************************************************************

	subroutine compare_files(idiff_end)

c checks two files written with check_debug from ht

	implicit none

	INTERFACE 
	  subroutine alloc_int(n,ival,niv,iv)
	  integer n,niv
	  integer ival(n)
	  integer, allocatable :: iv(:)
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine r_info(nh,nv,rval1,rval2,ipv,ipev,text)
	  integer nh,nv
	  real rval1(nh*nv)
	  real rval2(nh*nv)
	  integer ipv(:),ipev(:)
	  character*(*) text
	  END subroutine
	END INTERFACE

	INTERFACE 
	  subroutine i_info(nh,nv,ival1,ival2,ipv,ipev,text,iunit)
	  integer nh,nv
	  integer ival1(nh*nv)
	  integer ival2(nh*nv)
	  integer ipv(:),ipev(:)
	  character*(*) text
	  integer, optional :: iunit
	  END subroutine
	END INTERFACE

	integer idiff_end

	integer ndim
	parameter (ndim=2000000)

	character*60 name_one,name_two
	character*80 text1,text2,text
	logical bcheck,bstop,bverbose
	integer nt1,nt2,nt
	integer nh1,nh2,nh
	integer nv1,nv2,nv
	integer nrec,ntot,ntime
	integer i,idiff,idiff_tot
	integer nc
	integer nipv,nipev
	integer ios
	integer, allocatable :: ipv(:),ipev(:)
	double precision dtime,dtime1,dtime2

	integer, allocatable :: ival1(:),ival2(:)
	real, allocatable :: rval1(:),rval2(:)

	allocate(ival1(ndim))
	allocate(ival2(ndim))
	allocate(rval1(ndim))
	allocate(rval2(ndim))

	bstop = .false.			!stop on error
	bstop = .true.			!stop on error
	bverbose = .false.		!write info on all records
	bverbose = .true.		!write info on all records
	bcheck = .false.		!check for differences
	bcheck = .true.			!check for differences

	nc = command_argument_count()
	if( nc .ne. 2 ) then
	  write(6,*) 'Usage: check_debug file1 file2'
	  stop 'error stop check_debug: no files given'
	end if

	call get_command_argument(1,name_one)
	call get_command_argument(2,name_two)

	open(1,file=name_one,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error opening file 1'
	open(2,file=name_two,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error opening file 2'

	write(6,*) 'file 1: ',trim(name_one)
	write(6,*) 'file 2: ',trim(name_two)

	idiff_tot = 0
	ntime = 0
	nipv = 0
	nipev = 0

	do while(.true.)

	  read(1,end=9) dtime1
	  read(2,end=9) dtime2
	  if( dtime1 .ne. dtime2 ) goto 99
	  dtime = dtime1
	  write(6,*) 'time = ',dtime
	  ntime = ntime + 1

	  if( bverbose ) then
	    write(6,*) '       irec          nh          nv' //
     +			'          nt        diff name'
	  end if

	  nrec = 0
	  do while(.true.)
	    read(1) nh1,nv1,nt1
	    read(2) nh2,nv2,nt2
	    if( nh1 .ne. nh2 ) goto 98
	    if( nv1 .ne. nv2 ) goto 98
	    if( nt1 .ne. nt2 ) goto 98
	    nh = nh1
	    nv = nv1
	    nt = nt1
	    ntot = nh*nv

	    if( nt .eq. 0 ) exit
	    nrec = nrec + 1

	    read(1,end=9) text1
	    read(2,end=9) text2
	    text = text1
	    if( nt .gt. ndim ) goto 97
	    if( text1 .ne. text2 ) goto 96

	    idiff = 0

	    if( .not. bcheck ) then
	      read(1)
	      read(2)
	      if( bverbose ) write(6,*) nrec,nh,nv,nt,idiff,trim(text)
	      cycle
	    end if

	    if( nt == 1 ) then			!integer
	      read(1) (ival1(i),i=1,ntot)
	      read(2) (ival2(i),i=1,ntot)
	      if( bcheck ) then
	        call check_ival(dtime,nrec,nh,nv,ival1,ival2,idiff)
	      end if
	      if( text == 'ipv' ) call alloc_int(nh,ival1,nipv,ipv)
	      if( text == 'ipev' ) call alloc_int(nh,ival1,nipev,ipev)
	      if( idiff > 0 .and. bverbose ) then
	        call i_info(nh,nv,ival1,ival2,ipv,ipev,text)
	      end if
	    else if( nt == 2 ) then		!real
	      read(1) (rval1(i),i=1,ntot)
	      read(2) (rval2(i),i=1,ntot)
	      if( bcheck ) then
	        call check_rval(dtime,nrec,nh,nv,rval1,rval2,idiff)
	      end if
	      if( idiff > 0 .and. bverbose ) then
	        call r_info(nh,nv,rval1,rval2,ipv,ipev,text)
	      end if
	    else
	      write(6,*) 'cannot handle nt = ',nt
	      stop 'error stop: nt'
	    end if

	    if( idiff > 0 .or. bverbose ) then
	      write(6,*) nrec,nh,nv,nt,idiff,trim(text)
	    end if
	    idiff_tot = idiff_tot + idiff

	  end do

	  if( bstop .and. idiff_tot > 0 ) exit
	  if( bverbose ) write(6,*) 'nrecs checked: ',nrec
	end do

    9	continue

	close(1)
	close(2)

	write(6,*) 'total time records read: ',ntime
	write(6,*) 'total differences found: ',idiff_tot

	idiff_end = idiff_tot

	return
   99	continue
	write(6,*) dtime1,dtime2
	stop 'error stop check_debug: time mismatch'
   98	continue
	write(6,*) nh1,nh2,nv1,nv2,nt1,nt2
	stop 'error stop check_debug: size or type mismatch'
   97	continue
	write(6,*) dtime,nrec,ntot,ndim
	stop 'error stop check_debug: dimension'
   96	continue
	write(6,*) trim(text1),' - ',trim(text2)
	stop 'error stop check_debug: text mismatch'
	end

c*******************************************************************

	subroutine check_rval(dtime,nrec,nh,nv,val1,val2,idiff)

	implicit none

	double precision dtime
	integer nrec
	integer nh,nv,idiff
	real val1(nh*nv)
	real val2(nh*nv)

	integer i,k,l,ntot

	idiff = 0
	ntot = nh*nv

	do i=1,ntot
	  if( val1(i) .ne. val2(i) ) then
	    k = 1 + (i-1)/nv
	    l = 1 + mod(i-1,nv)
	    write(77,*) dtime,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

c*******************************************************************

	subroutine check_ival(dtime,nrec,nh,nv,val1,val2,idiff)

	implicit none

	double precision dtime
	integer nrec
	integer nh,nv,idiff
	integer val1(nh*nv)
	integer val2(nh*nv)

	integer i,k,l,ntot

	idiff = 0
	ntot = nh*nv

	do i=1,ntot
	  if( val1(i) .ne. val2(i) ) then
	    k = 1 + (i-1)/nv
	    l = 1 + mod(i-1,nv)
	    write(77,*) dtime,nrec,k,l,val1(i),val2(i)
	    idiff = idiff + 1
	  end if
	end do

	end

c*******************************************************************

	subroutine r_info(nh,nv,rval1,rval2,ipv,ipev,text)

	implicit none

	integer nh,nv
	real rval1(nh*nv)
	real rval2(nh*nv)
	integer ipv(:),ipev(:)
	character*(*) text

	integer i,ih,iv
	integer ipvv(nh)

	if( nh == size(ipv) ) then
	  ipvv = ipv
	else if( nh == size(ipev) ) then
	  ipvv = ipev
	else
	  write(6,*) nh,size(ipv),size(ipev)
	  stop 'error stop r_info: unknown nh'
	end if

	write(6,*) 'differences found reading ',trim(text)
	write(6,*) '         irec       k       l    kext' //
     +			'              val1              val2'
!                        123456789012345678901234567890123456

	do i=1,nh*nv
	  if( rval1(i) /= rval2(i) ) then
	    iv = 1 + mod((i-1),nv)
	    ih = 1 + (i-1)/nv
	    write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),rval1(i),rval2(i)
	  end if
	end do

 1000	format(a,4i8,2f18.6)
	end

c*******************************************************************

	subroutine i_info(nh,nv,ival1,ival2,ipv,ipev,text,iunit)

	implicit none

	integer nh,nv
	integer ival1(nh*nv)
	integer ival2(nh*nv)
	integer ipv(:),ipev(:)
	character*(*) text
	integer, optional :: iunit

	integer i,ih,iv,iu
	integer ipvv(nh)

	iu = 0
	if( present(iunit) ) iu = iunit

	if( nh == size(ipv) ) then
	  ipvv = ipv
	else if( nh == size(ipev) ) then
	  ipvv = ipev
	else
	  write(6,*) nh,size(ipv),size(ipev)
	  stop 'error stop i_info: unknown nh'
	end if

	write(6,*) 'differences found reading ',trim(text)
	write(6,*) '         irec       k       l    kext' //
     +			'              val1              val2'
	write(iu,*) '         irec       k       l    kext' //
     +			'              val1              val2'
!                        123456789012345678901234567890123456

	do i=1,nh*nv
	  iv = 1 + mod((i-1),nv)
	  ih = 1 + (i-1)/nv
	  if( iu > 0 ) then
	    write(iu,1000) 'diff: ',i,ih,iv,ipvv(ih),ival1(i),ival2(i)
	  else if( ival1(i) /= ival2(i) ) then
	    write(6,1000) 'diff: ',i,ih,iv,ipvv(ih),ival1(i),ival2(i)
	  end if
	end do

 1000	format(a,4i8,2i18)
 2000	format(4i8,2i18)
	end

c*******************************************************************

	subroutine alloc_int(n,ival,niv,iv)

	implicit none

	integer n,niv
	integer ival(n)
	integer, allocatable :: iv(:)

	allocate(iv(n))
	iv = ival
	niv = n

	end

c*******************************************************************

