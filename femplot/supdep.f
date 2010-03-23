c
c $Id: supdep.f,v 1.5 2009-11-18 17:16:00 georg Exp $
c
c routines for averaging depth
c
c contents :
c
c subroutine mkhv(hv,auxv,nkn,nel)	makes hv (nodewise depth)
c subroutine mkhev(hev,nel)		makes hev (elementwise depth)
c subroutine mkht(hetv,href)		makes hetv (elem depth of actual layer)
c subroutine mkht3(nlvdim,het3v,href)	makes het3v (3D depth structure)
c
c revision log :
c
c 26.05.2000    ggu     routines written from scratch
c 17.09.2008    ggu     routine mkht changed for layer = -1
c 13.10.2009    ggu     new routine mkht3
c
c******************************************************************

	subroutine mkhv(hv,auxv,nkn,nel)

c makes hv (nodewise depth)

	implicit none

c arguments
	real hv(1)
	real auxv(1)
	integer nkn,nel
c common
	integer nen3v(3,1)
	real hm3v(3,1)
	common /nen3v/nen3v, /hm3v/hm3v
c local
	integer ie,ii,k,kn

        do k=1,nkn
          hv(k) = 0.
        end do

        do ie=1,nel
          do ii=1,3
            kn=nen3v(ii,ie)
            hv(kn)=hv(kn)+hm3v(ii,ie)
	    auxv(kn)=auxv(kn)+1.
          end do
        end do

        do k=1,nkn
          hv(k) = hv(k) / auxv(k)
        end do

	return
	end

c******************************************************************

	subroutine mkhev(hev,nel)

c makes hev (elementwise depth)

	implicit none

c arguments
	real hev(1)
	integer nel
c common
	real hm3v(3,1)
	common /hm3v/hm3v
c local
	integer ie,ii
	real h

        do ie=1,nel
	  h=0.
          do ii=1,3
            h=h+hm3v(ii,ie)
          end do
	  hev(ie) = h / 3.
        end do

	return
	end

c******************************************************************

	subroutine mkht(hetv,href)

c makes hetv (elementwise depth of actual layer)
c
c uses level to decide what to do

	implicit none

c arguments
	real hetv(1)
	real href
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhv(1)
        common /ilhv/ilhv
	real zenv(3,1)
	common /zenv/zenv
        real vev(1)
        common /vev/vev
c local
	logical bdebug
	integer ie,ii
	integer level,lmax
	real z
c functions
	integer getlev

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.

	level = getlev()

c-------------------------------------------------------------------
c compute water level variation -> store in vev
c-------------------------------------------------------------------

	if( level .le. 1 ) then		!we could need water level variation
          do ie=1,nel
            z = 0.
            do ii=1,3
              z = z + zenv(ii,ie)
            end do
            z = z / 3.
            vev(ie) = z - href
          end do
	end if

c-------------------------------------------------------------------
c handle different kind of levels
c-------------------------------------------------------------------

	if( level .eq. 0 ) then			! barotropic integrated

          do ie=1,nel
	    z = vev(ie)
            hetv(ie) = hev(ie) + z
          end do

	else if( level .eq. 1 ) then		! surface layer

          do ie=1,nel
            z = vev(ie)
	    if( ilhv(ie) .eq. 1 ) then
              hetv(ie) = hev(ie) + z
	    else
              hetv(ie) = hlv(1) + z
	    end if
	  end do

	else if( level .eq. -1 ) then		! bottom layer

          do ie=1,nel
	    lmax = ilhv(ie)
	    if( lmax .eq. 1 ) then
              hetv(ie) = hev(ie) + z
	    else
              hetv(ie) = hev(ie) - hlv(lmax-1)
	    end if
	  end do

	else					! inner layer

          do ie=1,nel
	    if( ilhv(ie) .eq. level ) then
              hetv(ie) = hev(ie) - hlv(level-1)
	    else
              hetv(ie) = hlv(level) - hlv(level-1)
	    end if
	  end do

	end if

c-------------------------------------------------------------------
c debug output
c-------------------------------------------------------------------

	if( bdebug ) then
	  ie = 1644
	  write(6,*) 'debugging mkht...'
	  write(6,*) ilhv(ie),hetv(ie),hev(ie)
	  do ii=1,ilhv(ie)
	    write(6,*) hlv(ii)
	  end do
	end if

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c******************************************************************

	subroutine mkht3(nlvdim,het3v,href)

c makes het3v (3D depth structure)

	implicit none

c arguments
	integer nlvdim
	real het3v(nlvdim,1)
	real href
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real hlv(1), hev(1)
        common /hlv/hlv, /hev/hev
        integer ilhv(1)
        common /ilhv/ilhv
	real zenv(3,1)
	common /zenv/zenv
        real vev(1)
        common /vev/vev
c local
	logical bdebug
	integer ie,ii
	integer l,lmax
	real z

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

        bdebug = .true.
        bdebug = .false.

c-------------------------------------------------------------------
c compute water level variation -> store in vev
c-------------------------------------------------------------------

        do ie=1,nel
          z = 0.
          do ii=1,3
            z = z + zenv(ii,ie)
          end do
          z = z / 3.
          vev(ie) = z - href
        end do

c-------------------------------------------------------------------
c compute layer thickness
c-------------------------------------------------------------------

	do ie=1,nel
	  z = vev(ie)
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( lmax .eq. 1 ) then			!only one layer
              het3v(l,ie) = hev(ie) + z
	    else if( l .eq. 1 ) then			!surface layer
              het3v(l,ie) = hlv(1) + z
	    else if( l .eq. lmax ) then			!bottom layer
              het3v(l,ie) = hev(ie) - hlv(lmax-1)
	    else					!inner layer
              het3v(l,ie) = hlv(l) - hlv(l-1)
	    end if
	  end do
	end do

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	end

c******************************************************************

