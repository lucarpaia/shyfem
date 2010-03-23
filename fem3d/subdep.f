c
c $Id: subdep.f,v 1.10 2008-07-16 15:41:39 georg Exp $
c
c depth utility routines
c
c contents :
c
c function igtdep(k,f)                  get depth for node
c function igtdpa(mode,h)               gets unique depth for all nodes
c subroutine huniqu(hev,hkv)            make depth unique for every node
c subroutine makehev(hev)		makes hev (elementwise depth)
c subroutine makehkv(hkv,haux)		makes hkv (nodewise depth)
c subroutine depadj(hmin,hmax,href)	adjusts depth to ref/min/max values
c
c revision log :
c
c 29.06.1997	ggu	depth routines in one file
c 06.11.1998	ggu	new huniqu to compute hev and hkv
c 19.10.1999	ggu	new routine makehv from subutl
c 25.03.2002	ggu	new routines makehkv (before makehv) and makehev
c 28.11.2005	ggu	makehkv changed (uses real aux value, area weight)
c 24.02.2006	ggu	bug in makehkv -> haux was integer
c 18.10.2006	ccf	bug in makehkv -> no area multiplication
c
c********************************************************************

	function igtdep(k,f,ndim)

c gets depth given a node number
c if different values of depth are associated
c with a node, it returns all these values
c
c k             node number (internal)
c f             vector in which the depth values are
c               ...stored at return
c ndim		dimension of f
c igtdep        number of different values found
c
c revised 02.02.94 by ggu	$$nmax - check error condtion nmax
c revised 29.06.97 by ggu	$$ndim - dimension of f is passed

	implicit none

	integer igtdep
	integer k,ndim
	real f(1)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

	integer iact,ie,i,ii

	iact=0
	do ie=1,nel
	  do ii=1,3
	    if(nen3v(ii,ie).eq.k) then
		do i=1,iact
		    if(f(i).eq.hm3v(ii,ie)) goto 1
		end do

		iact=iact+1     		!new depth
		if(iact.gt.ndim) goto 99	!$$nmax !$$ndim
		f(iact)=hm3v(ii,ie)

    1           continue       			!old depth
	    end if
	  end do
	end do

	igtdep=iact

	return
   99	continue
	stop 'error stop igtdep : nmax'		!$$nmax
	end

c********************************************************************

	function igtdpa(mode,h)

c gets unique depth for all nodes
c
c mode          switch
c               1       deepest value is returned
c               -1      most shallow value
c h             vector in which the depth values are
c               ...stored (return value)
c igtdpa        return status
c               1       unique depth
c               0       no unique depth
c               -1      error

	implicit none

	integer igtdpa
	integer mode
	real h(1)

	real high
	parameter(high=1.e+30)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

	logical buniq
	integer ie,ii,i,k
	real hh,hhh,hflag

	buniq=.true.

	if(mode.eq.1) then
		hflag=-high
	else if(mode.eq.-1) then
		hflag=high
	else
		write(6,*) 'Value for mode not allowed :',mode
		igtdpa=-1
		return
	end if

	do i=1,nkn
	   h(i)=hflag
	end do

	do ie=1,nel
	 do ii=1,3
	   k=nen3v(ii,ie)
	   hh=hm3v(ii,ie)
	   hhh=h(k)
	   if(mode.eq.1) then
		if(hh.gt.h(k)) h(k)=hh
	   else
		if(hh.lt.h(k)) h(k)=hh
	   end if
	   if(hhh.ne.hflag.and.hhh.ne.hh) buniq=.false.
	 end do
	end do

	do i=1,nkn
	   if(h(i).eq.hflag) then
		write(6,*) 'igtdpa : Nodes without depth'
		igtdpa=-1
		return
	   end if
	end do

	if(buniq) then
		igtdpa=1
	else
		igtdpa=0
	end if

	end

c********************************************************************

	subroutine huniqu(hev,hkv)

c make depth unique for every node (changes hm3v)
c nodal values are the highest (deepest) value
c
c hev		element averaged depth values
c hkv            array with unique depth values

	implicit none

	real hev(1)
	real hkv(1)

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

	integer ie,ii,k
	logical bstop
	real h,flag,hm

	integer ipext

c flag nodal values

	flag = -999.

	do k=1,nkn
	  hkv(k) = flag
	end do

c create element averaged depth values and assign to nodal values

	do ie=1,nel
	  hm = 0.
	  do ii=1,3
	    hm = hm + hm3v(ii,ie)
	  end do

	  hev(ie) = hm / 3.

	  do ii=1,3
	    hm3v(ii,ie) = hev(ie)
	  end do

	  h = hev(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 ) write(6,*) 'huniqu: ',ie,ii,k,hev(ie)
	    if( h .gt. hkv(k) ) hkv(k) = h
	  end do
	end do

c check if all depth values are available

	bstop = .false.

	do k=1,nkn
	  if( hkv(k) .eq. flag ) then
		write(6,*) 'No depth for node ',ipext(k)
		bstop = .true.
	  end if
	end do

	if( bstop ) stop 'error stop huniqu'

	end

c********************************************************************

        subroutine makehev(hev)

c makes hev (elementwise depth)

        implicit none

c arguments
        real hev(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real hm3v(3,1)
        common /hm3v/hm3v
c local
        integer ie,ii
	real hm

        do ie=1,nel
	  hm = 0.
          do ii=1,3
	    hm = hm + hm3v(ii,ie)
          end do
	  hev(ie) = hm / 3.
        end do

        end

c********************************************************************

        subroutine makehkv(hkv,haux)

c makes hkv (nodewise depth)

        implicit none

c arguments
        real hkv(1)
        real haux(1)   !aux array -> bug - was integer
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer nen3v(3,1)
        common /nen3v/nen3v
        real hm3v(3,1)
        common /hm3v/hm3v
	include 'ev.h'
c local
        integer ie,ii,k,kn
	real area

        do k=1,nkn
          hkv(k) = 0.
          haux(k) = 0.
        end do

        do ie=1,nel
	  area = ev(10,ie)
          do ii=1,3
            kn=nen3v(ii,ie)
            hkv(kn)=hkv(kn)+hm3v(ii,ie)*area	!ccf
            haux(kn)=haux(kn)+area
          end do
        end do

        do k=1,nkn
          hkv(k) = hkv(k) / haux(k)
        end do

        end

c********************************************************************

	subroutine depadj(hmin,hmax,href)

c adjusts depth to reference and min/max values

	implicit none

	real hmin,hmax,href

	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

	integer iaux,ie,ii
	real hmed

c adjust depth to constant in element %%%%%%%%%%%%%%%%%%%%%%

        do ie=1,nel
          hmed=0.
          do ii=1,3
            hmed=hmed+hm3v(ii,ie)
          end do
          hmed=hmed/3.
          do ii=1,3
            hm3v(ii,ie)=hmed
          end do
        end do

c adjust depth to minimum depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	iaux=0
	do ie=1,nel
	 do ii=1,3
	  if(hm3v(ii,ie).lt.hmin) then
	    hm3v(ii,ie)=hmin
	    iaux=iaux+1
	  end if
	 end do
	end do

	if(iaux.gt.0) then
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '*********** hmin = ',     hmin  ,' ************'
	  write(6,*) '*********** changes = ',  iaux  ,' ************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	end if

c adjust depth to maximum depth %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	iaux=0
	do ie=1,nel
	 do ii=1,3
	  if(hm3v(ii,ie).gt.hmax) then
	    hm3v(ii,ie)=hmax
	    iaux=iaux+1
	  end if
	 end do
	end do

	if(iaux.gt.0) then
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '*********** hmax = ',     hmax  ,' ************'
	  write(6,*) '*********** changes = ',  iaux  ,' ************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	  write(6,*) '***********************************************'
	end if

c adjust depth to reference level %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do ie=1,nel
	 do ii=1,3
	  hm3v(ii,ie)=hm3v(ii,ie)-href
	 end do
	end do

	end

c********************************************************************

