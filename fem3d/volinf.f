c
c $Id: volinf.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 06.04.1999	ggu	cosmetic changes
c
c*****************************************************************

	program volinf

c reads flx files

	implicit none

	integer nfxdim
	parameter( nfxdim = 5000 )

	integer nin
	integer idfile,nvers,nvols,kvolm,ivolm
	integer i,it
	integer kvol(nfxdim)
	integer ivol(nfxdim)
	real ppp(nfxdim)

	integer iapini,ideffi

c---------------------------------------------------------------

        if(iapini(2,1,1,0).eq.0) then
                stop 'error stop: iapini'
        end if

        nin=ideffi('datdir','runnam','.vol','unform','old')
        if(nin.le.0) then
	  stop 'error stop: cannot open file'
	end if

	read(nin) idfile,nvers
	read(nin) nvols,kvolm

	if( nvols .gt. nfxdim .or. kvolm .gt. nfxdim ) then
	  write(6,*) 'nvols,kvolm: ',nvols,kvolm
	  stop 'error stop: nfxdim'
	end if

	read(nin) (kvol(i),i=1,kvolm)

	read(nin) ivolm

	if( ivolm .gt. nfxdim ) then
	  write(6,*) 'ivolm: ',ivolm
	  stop 'error stop: nfxdim'
	end if

	read(nin) (ivol(i),i=1,ivolm)

	write(6,*) idfile,nvers,nvols,kvolm,ivolm

	do while(.true.)
	  read(nin,end=2) it,nvols,(ppp(i),i=1,nvols)
c	  write(6,*) it,nvols
c	  write(6,*) (ppp(i),i=1,nvols)
	  write(6 ,'(i8,10e12.4)') it,(ppp(i),i=1,nvols)
	  write(67,'(i8,10e12.4)') it,(ppp(i),i=1,nvols)
	end do

    2	continue

	close(nin)

	end
