c
c $Id: baschk.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	completely restructured
c 04.06.1999	ggu	new statistics are computed
c 06.04.2005    ggu     read param.h
c
c****************************************************************

        program baschk

c writes information on basin about nodes and elements

	implicit none

	include 'param.h'

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

	include 'evmain.h'

	logical bnode,belem
	integer iapini

c-----------------------------------------------------------------
c read in basin
c-----------------------------------------------------------------

        if( iapini(1,nkndim,neldim,0) .le. 0 ) stop

c-----------------------------------------------------------------
c set up element info
c-----------------------------------------------------------------

	call sp110a

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
c specific info
c-----------------------------------------------------------------

	call bascheck

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	end

c*******************************************************************

	subroutine bascheck

c writes statistics on basin

	implicit none

	include 'param.h'

        character*80 descrp
        common /descrp/ descrp

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real grav,fcor,dcor,dirn,rowass,roluft
        common /pkonst/ grav,fcor,dcor,dirn,rowass,roluft

	include 'evmain.h'

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

	integer ie,ii,k,i
	integer imin,imax
	real area,amin,amax
	real x(3),y(3)
	real xmin,xmax,ymin,ymax
	real dxmax,dymax
	real h,w,eps
	integer iang,ic

	integer count(nkndim)
	integer icount(20)

	real areatr
	integer ipext

	eps = 1.e-5

c-----------------------------------------------------------------
c area code
c-----------------------------------------------------------------

	do k=1,nkn
	  count(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    count(k) = count(k) + 1
	  end do
	end do

	do k=1,nkn
	  if( count(k) .le. 1 ) then
	    write(6,*) 'low count for node ',ipext(k),' : ',count(k)
	  end if
	end do

c-----------------------------------------------------------------
c angle
c-----------------------------------------------------------------

	do i=1,20
	  icount(i) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    w = ev(10+ii,ie)
	    iang = (w-90.)/10. + 1. - eps
	    if( iang .gt. 0 ) then
	      icount(iang) = icount(iang) + 1
c	      write(6,*) k,ii,w,iang
	    end if
	  end do
	end do

	do i=1,20
	  ic = icount(i)
	  if( ic .gt. 0 ) then
	    write(6,*) 'angle > ',80+10*i,' : ',ic
	  end if
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    w = ev(10+ii,ie)
	    if( w .gt. 120. ) then
	      write(6,*) 'big angle ',ipext(k),' : ',w
	    end if
	  end do
	end do

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	write(6,*)

	end

c*******************************************************************

        function areatr(ie)

c determination of area of element
c
c ie            number of element (internal)
c areatr        element area (return value)

	real areatr
	integer ie

	include 'param.h'

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer nen3v(3,neldim)
        common /nen3v/nen3v

	real aj
	integer ii,i1,i2,k1,k2

        aj=0.
        do ii=1,3
          i1=mod(ii,3)+1
          i2=mod(i1,3)+1
          k1=nen3v(i1,ie)
          k2=nen3v(i2,ie)
          aj=aj+xgv(k1)*ygv(k2)-xgv(k2)*ygv(k1)
        end do

        areatr = aj / 2.

        end

c*******************************************************************

