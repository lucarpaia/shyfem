c
c $Id: subnetcdf_admin.f,v 1.1 2008-04-11 16:05:37 georg Exp $
c
c NetCDF file administration routines
c Writes standardized NetCDF files for Hydro-models.
c Link to Hydro_netcdfs_fem.f  
c output are written to runnam.nc
c
c contents :
c
c subroutine wrnetcdf
c subroutine set_ele
c subroutine set_bnd
c subroutine setvelocity
c
c revision log :
c
c 07.11.2007	ccf	start writing NetCDF output just for u,v and zeta
c
c********************************************************

	subroutine wrnetcdf

c writes and administers NetCDF file

	implicit none

	include 'param.h'

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer itanf,itend,idt,nits,niter,it
	common /femtim/ itanf,itend,idt,nits,niter,it
        integer nlvdi,nlv
        common /level/ nlvdi,nlv

        real xgv(nkndim),ygv(nkndim)
        common /xgv/xgv,/ygv/ygv
        real hkv(nkndim)
        common /hkv/hkv
        real hlv(nkndim)
	common /hlv/hlv
        real znv(nkndim)
        common /znv/znv
        character*80 descrp
        common /descrp/ descrp

	integer itmout
	integer iround
	real getpar

	integer idtout,itout
	integer icall
        integer date
        integer year,month,dayy,time

	integer netcdf
        character*80 nam,bas,dir
        character*80 ncfile
	integer ncid
	integer nb			!number of boundary segments
	integer ele(3,neldim)
	integer bnd(4,nkndim)
        real day			!time in days
	integer ibasedate(4)		!iyear, imonth, iday, ihour of base date (time = 0)
        character globalstr(9)*40
	real u(nkn,nlvdim),v(nkn,nlvdim)
	integer iaux(nkndim)

        data globalstr/
     &  'Triangular','Z','SHYFEM'
     & ,'','',''
     & ,'ISMAR-CNR','original','c.ferrarin@ismar.cnr.it'/

	data icall /0/

	save ncfile,ncid
	save ibasedate,day
	save nb,bnd,ele
	save idtout,itout
	save icall

c       ----------------------------------------------------------
c       Initialization
c       ----------------------------------------------------------

        if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  idtout=iround(getpar('idtout'))
	  itmout=iround(getpar('itmout'))

    	  if( itmout .le. itanf ) itmout = itanf !$$ITMOUT

 	  netcdf=iround(getpar('netcdf'))
          if( netcdf .le. 0 ) icall = -1
	  if( idtout .le. 0 ) icall = -1
	  if( icall .eq. -1 ) return
		
	  itout = itmout

c         ----------------------------------------------------------
c         Get NetCDF output file name (*.nc)
c         ----------------------------------------------------------

          call getfnm('datdir',dir)
          call getfnm('runnam',nam)
          call mkname(dir,nam,'.nc',ncfile)

          call getfnm('basnam',bas)
	  globalstr(4) = descrp		!title
	  globalstr(5) = nam		!name of simulation
	  globalstr(6) = bas		!name of basin

c         ----------------------------------------------------------
c         Set initial time and date from str file
c         ----------------------------------------------------------

          day = 0.
	  date = getpar('date')
	  call unpackdate(date,year,month,dayy)
	  time = getpar('time')

          ibasedate(1) = year + 2000
          ibasedate(2) = month
          ibasedate(3) = dayy
          ibasedate(4) = time

c         ----------------------------------------------------------
c         Set ele and bnd arrays
c         ----------------------------------------------------------

	  call set_ele(nel,ele)
	  call set_bnd(nkn,iaux,nb,bnd)

c         ----------------------------------------------------------
c         Initialize NetCDF output file
c         ----------------------------------------------------------

          call write_netcdf_Hydro_fem(ncfile,ncid,1, 
     & 	    globalstr,nel,nkn,nlv,nb,ele,bnd,day,ibasedate,
     &	    xgv,ygv,hlv,hkv,1.,1.,1.,-1.,-1.,-1.,-1.,-1.)

	  icall =  1

	end if

	if( it .lt. itout ) goto 123

c       ----------------------------------------------------------
c       Write NetCDF every output time step
c       ----------------------------------------------------------

        day = it / 86400.

	call setvelocity(nkn,u,v)

        call write_netcdf_Hydro_fem(ncfile,ncid,2, 
     & 	  globalstr,nel,nkn,nlv,nb,ele,bnd,day,ibasedate,
     &	  xgv,ygv,hlv,hkv,znv,u,v,-1.,-1.,-1.,-1.,-1.)

	itout=itout+idtout

 123	continue

c       ----------------------------------------------------------
c       Close NetCDF file if it = itend
c       ----------------------------------------------------------

	if( it .eq. itend ) then
          call write_netcdf_Hydro_fem(ncfile,ncid,3, 
     & 	    globalstr,nel,nkn,nlv,nb,ele,bnd,day,ibasedate,
     &	    xgv,ygv,hlv,hkv,znv,u,v,-1.,-1.,-1.,-1.,-1.)
	    return
	end if

	end

c********************************************************

	subroutine set_ele(nel,ele)

c return array ele containing node number for each element

	implicit none

        integer nel     	!number of element
	integer ie,ii
	integer nen3v(3,1)
        common /nen3v/nen3v

	integer ele(3,1)	!Connectivity of triangular elements

	do ie = 1,nel
	  do ii = 1,3
	    ele(ii,ie) = nen3v(ii,ie)
	  end do
	end do

	end

c********************************************************

	subroutine set_bnd(nkn,iaux,nb,bnd)

c return total number of boundary segment and array bnd
c bnd(4,nb)  Indices of nodes making up boundary segments
c            bnd(1),bnd(2) nodes of a boundary segment (water on right)
c            bnd(3)  Island number of this segment
c            bnd(4)  Land segment=0,  Water segment=1

	implicit none

	integer nkn
	integer nb			!total number of boundary nodes
	integer bnd(4,1)

	integer kantv(2,1)
	common /kantv/kantv
        integer inodv(1)
        common /inodv/inodv

	integer iaux(1)
	integer isl
	integer k,knext,kstart,k1

	isl = 0
	nb = 0
	do k=1,nkn
	  iaux(k) = 0
	end do

	do k=1,nkn
	  if( kantv(2,k) .ne. 0 .and. iaux(k) .eq. 0 ) then
		isl = isl + 1
		kstart = k
		k1 = k
		do while( kantv(2,k1) .ne. kstart )
		  nb = nb + 1
		  iaux(k1) = isl
		  knext = kantv(2,k1)
		  bnd(1,nb) = k1
		  bnd(2,nb) = knext
		  bnd(3,nb) = isl
		  bnd(4,nb) = 0
	          if(inodv(k1).gt.0. .and. inodv(knext).gt.0) bnd(4,nb)=1
		  k1 = knext
		end do
		nb = nb + 1
		iaux(knext) = isl
		bnd(1,nb) = knext
		bnd(2,nb) = kstart
		bnd(3,nb) = isl
		bnd(4,nb) = 0
	        if(inodv(knext).gt.0. .and. inodv(kstart).gt.0) bnd(4,nb)=1
	  end if
	end do

	write(*,*)'total number of boundary segments ', nb

	end

c********************************************************

	subroutine setvelocity(nkn,u,v)

	implicit none

	include 'param.h'

        real uprv(nlvdim,1),vprv(nlvdim,1)
	real wprv(0:nlvdim,1)
        common /uprv/uprv, /vprv/vprv, /wprv/wprv
	integer k,l,nkn

	real u(1,nlvdim),v(1,nlvdim)

	do k = 1,nkn
	  do l = 1,nlvdim
	     u(k,l) = uprv(l,k)
	     v(k,l) = vprv(l,k)
	  end do
	end do

	end
