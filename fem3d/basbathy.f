c
c $Id: basbathy.f,v 1.7 2010-02-26 17:35:06 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 08.09.2003	ggu	mode 5 -> write depth values from elements
c 23.09.2004    ggu     interpolq() changed for bathy interpolation
c 02.10.2004    ggu     interpole() for exponential interpolation
c 12.05.2005    ggu     pass hmin to interpolation functions
c 06.04.2009    ggu     read param.h
c 24.04.2009	ggu	new call to rdgrd()
c 21.05.2009	ggu	restructured to allow for nodal interpolation
c 16.12.2010	ggu	bug fix in transfer_depth()
c 02.12.2011	ggu	introduction of nminimum - hardcoded for now
c 16.03.2012	ggu	autoregression introduced (make_auto_corr,interpola)
c 16.03.2012	ggu	default value for umfact set to 3, new mode = 3
c 01.06.2012	ggu	some more changes
c 13.06.2013	ggu	copy_depth() renamed to transfer_depth()
c 13.02.2014	ggu	new data written, can read also bas file
c 05.03.2014	ggu	subroutines copied to other routine
c
c****************************************************************

        program basbathy

c performs bathymetry interpolation in basin
c
c takes care of lat/lon coordinates

	implicit none

	include 'param.h'

	integer ndim
	parameter(ndim=1200000)
	real xp(ndim)
	real yp(ndim)
	real dp(ndim)
	real ap(ndim)

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
        real hkv(nkndim)
        common /hkv/hkv

        real raux(neldim)
        integer iaux(neldim)
        integer ipaux(nkndim)

	include 'evmain.h'

        character*40 bfile,gfile,nfile
        character*60 line
	integer node,nit
	integer mode,np,n,i,nt
        integer ner,nco,nknh,nelh,nli
	integer nlidim,nlndim
	integer ike,idepth
	integer nminimum
	integer isphe
	real ufact,umfact
	real f(5)
	logical bstop,bbasin
	integer iscanf

        real xt(neldim)
        real yt(neldim)
        real at(neldim)
        real ht(neldim)

c-----------------------------------------------------------------
c what to do
c-----------------------------------------------------------------

	nminimum = 5	!for dwejra
	nminimum = 1	!minimum number of points to be used for interpolation

        write(6,*)
        write(6,*) 'I need the name of the basin file '
        write(6,*) '(the file can be in GRD or BAS format)'
        write(6,*) '(please include extension - default is GRD)'
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') gfile
        if( gfile .eq. ' ' ) stop
	write(6,*) 'grid is read from file : ', gfile
        write(6,*)

        write(6,*)
        write(6,*) 'I need the name of the bathymetry data file '
        write(6,*) '(the file must be in GRD format)'
        write(6,*)
	write(6,*) 'Enter file name: '
	read(5,'(a)') bfile
        if( bfile .eq. ' ' ) stop
	write(6,*) 'Bathymetry is read from file : ', bfile
        write(6,*)

        write(6,*)
        write(6,*) 'Two different algorithms are available:'
        write(6,*) '  1   exponential interpolation (default)'
        write(6,*) '  2   uniform interpolation on squares'
        write(6,*) '  3   exponential interpolation (autocorrelation)'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') mode
	if( mode .lt. 1 ) mode = 1
	if( mode .gt. 3 ) mode = 1
	write(6,*) 'Mode is : ', mode

        write(6,*)
        write(6,*) 'If there are some values missing you can:'
        write(6,*) '  1   interpolate on missing depth values (default)'
        write(6,*) '  2   interpolate on all elements/nodes'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') idepth
	if( idepth .ne. 2 ) idepth = 1
	write(6,*) 'Choice is : ', idepth

	ike = 1
	ufact = 1.
	umfact = 2.	!old default
	umfact = 3.

	if( mode .eq. 1 .or. mode .eq. 3 ) then

        write(6,*)
        write(6,*) 'For the exponential algorithm you can:'
        write(6,*) '  1   interpolate on elements (default)'
        write(6,*) '  2   interpolate on nodes'
        write(6,*)
	write(6,*) 'Enter choice: '
	read(5,'(i10)') ike
	if( ike .ne. 2 ) ike = 1
	write(6,*) 'Choice is : ', ike

        write(6,*)
	write(6,*) 'Enter parameters for expontential interpolation:'
        write(6,*)
	write(6,*) 'The std deviation is about the size of the elements'
	write(6,*) 'With ufact you can ultimately correct it (default=1)'
	write(6,*) 'The maximum radius is 3 times the standard deviation'
	write(6,*) 'With umfact you can correct it (default=3)'
        write(6,*)
	write(6,*) 'Enter params ufact and umfact (<CR> for default): '
        read(5,'(a)') line
        n = iscanf(line,f,2)
	if( n .lt. 0 .or. n .gt. 2 ) goto 95
        if( n .gt. 0 ) ufact = f(1)
        if( n .gt. 1 ) umfact = f(2)
        write(6,*) 'ufact,umfact :',ufact,umfact
        write(6,*)

	end if

c-----------------------------------------------------------------
c read in bathymetry file
c-----------------------------------------------------------------

	write(6,*) 'reading bathymetry file : ',bfile
	np = ndim
	call readgrd(bfile,np,xp,yp,dp)

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

	call check_basin_name(gfile,bbasin)

	if( bbasin ) then

	  write(6,*) 'reading basin as bas file...'
	  call read_basin(gfile,nkndim,neldim)

	else

	  write(6,*) 'reading basin as grd file...'

          ner = 6
          bstop = .false.

          nlidim = 0
          nlndim = 0
          call rdgrd(
     +                   gfile
     +                  ,bstop
     +                  ,nco,nkn,nel,nli
     +                  ,nkndim,neldim,nlidim,nlndim
     +                  ,ipv,ipev,iaux
     +                  ,iaux,iarv,iaux
     +                  ,hkv,hev,raux
     +                  ,xgv,ygv
     +                  ,nen3v
     +                  ,iaux,iaux
     +                  )

          if( bstop ) stop 'error stop rdgrd'

          call ex2in(nkn,3*nel,nlidim,ipv,ipaux,nen3v,iaux,bstop)
          if( bstop ) stop 'error stop ex2in'

	end if

c-----------------------------------------------------------------
c handling of depth and coordinates
c-----------------------------------------------------------------

	call check_spheric_ev			!sets lat/lon flag
	call get_coords_ev(isphe)
	call set_dist(isphe)

	call set_depth_i(idepth,nknh,nelh)

c-----------------------------------------------------------------
c general info
c-----------------------------------------------------------------

        write(6,*)
        write(6,*) ' nkn  = ',nkn, '  nel  = ',nel
        write(6,*) ' nknh = ',nknh,'  nelh = ',nelh
        write(6,*)

c-----------------------------------------------------------------
c node_test
c-----------------------------------------------------------------

	call node_test
	call set_ev

        if( ike .eq. 1 ) then                           !elementwise
          call prepare_on_elem(nt,xt,yt,at,ht,ufact)
        else                                            !nodewise
          call prepare_on_node(nt,xt,yt,at,ht,ufact)
        end if

c-----------------------------------------------------------------
c interpolate
c-----------------------------------------------------------------

	if( mode .eq. 1 ) then
	  call interpole(np,xp,yp,dp,nt,xt,yt,at,ht,umfact,nminimum)
        else if( mode .eq. 2 ) then
	  call interpolq(np,xp,yp,dp)
	else if( mode .eq. 3 ) then
	  call make_auto_corr(np,xp,yp,dp,ap,ufact)
	  call interpola(np,xp,yp,dp,ap,nt,xt,yt,at,ht)
        else
          write(6,*) 'wrong choice for mode : ',mode
          stop 'error stop'
	end if

	call transfer_depth(ike,ht)	!copy to nodes/elements

c-----------------------------------------------------------------
c write
c-----------------------------------------------------------------

	nfile = 'basbathy.grd'
	open(1,file=nfile,status='unknown',form='formatted')
	call wrgrd(1,ike)
	close(1)
        write(6,*) 'The new file has been written to ',nfile

	call write_data('basbathy.dat',nkn,hkv)
	!call write_xy('basbathy.xyz',nkn,ipv,xgv,ygv)

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	stop
   95	continue
	write(6,*) n,(f(i),i=1,n)
	write(6,*) line
	stop 'error stop basbathy: error in parameters'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

