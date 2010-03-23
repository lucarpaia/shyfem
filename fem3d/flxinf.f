c
c $Id: flxinf.f,v 1.10 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c
c notes :
c
c ppp		total flux through section
c ppppm(1,)	positive flux
c ppppm(2,)	negative flux (number is positive)
c
c	=> ppp = ppppm(1,) - ppppm(2,)
c
c*****************************************************************

	program flxinf

c reads flx files

	implicit none

	integer nfxdim
	parameter( nfxdim = 3000 )

	integer nin
	integer idfile,nvers,nsect,kfluxm
	integer i,it
	integer kflux(nfxdim)
	real ppp(nfxdim)
	real ppppm(2,nfxdim)

	integer iapini,ideffi

c---------------------------------------------------------------

        if(iapini(2,1,1,0).eq.0) then
                stop 'error stop: iapini'
        end if

        nin=ideffi('datdir','runnam','.flx','unform','old')
        if(nin.le.0) then
	  stop 'error stop: cannot open file'
	end if

	read(nin) idfile,nvers
	read(nin) nsect,kfluxm

	if( nsect .gt. nfxdim .or. kfluxm .gt. nfxdim ) then
	  write(6,*) idfile,nvers
	  write(6,*) nsect,kfluxm
	  write(6,*) nfxdim
	  stop 'error stop: nfxdim'
	end if

	read(nin) (kflux(i),i=1,kfluxm)

	write(6,*) idfile,nvers,nsect,kfluxm

	do while(.true.)

	  if( nvers .eq. 1 ) then
	    read(nin,end=2) it,nsect,(ppp(i),i=1,nsect)
	  else
	    read(nin,end=2) it,nsect,(ppp(i),i=1,nsect)
     +				,(ppppm(1,i),ppppm(2,i),i=1,nsect)
	  end if

c	  write(6,*) it,nsect
c	  write(6,*) (ppp(i),i=1,nsect)

	  write(6,'(i10,10i7)') it,(nint(ppp(i)),i=1,nsect)
	  write(67,'(i10,20f10.2)') it,(ppp(i),i=1,nsect)

	  !write(68,'(i10,20f10.2)') it,(ppppm(1,i),i=1,nsect)
	  !write(68,'(i10,20f10.2)') it,(ppppm(2,i),i=1,nsect)

          write(69,'(i10,20f10.2)') it,(ppppm(1,i)+ppppm(2,i)
     +						,i=1,nsect)

	end do

    2	continue

	close(nin)

	end
