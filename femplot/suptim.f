c
c $Id: suptim.f,v 1.8 2008-12-09 11:45:14 georg Exp $
c
c revision log :
c
c 12.02.1999  ggu     adapted to auto mode
c 27.05.2005  ggu     increase nrec always in oktime (even when it is the same)
c 13.11.2008  ggu     in oktime() increase irec only for new time
c 06.12.2008  ggu     in oktime() set itact to actual time
c
c******************************************************

        subroutine timeset(itanf,itend,idtout)

c set time limits

        implicit none

c arguments
        integer itanf,itend,idtout
c parameters
	integer ihigh
	parameter(ihigh=1000000000)
c common
        integer itmin,itmax,itfreq,nrec,idto,itact
	common /timlim/ itmin,itmax,itfreq,nrec,idto,itact
	save /timlim/

        itmin = itanf
        itmax = itend
	idto  = idtout
        itfreq = 1
	nrec = 0
	itact = itanf - 1

	if( idtout .le. 0 ) then
	  itmin = -ihigh
	  itmax =  ihigh
	end if

        end

c******************************************************

        subroutine timeask

c ask for time limits

        implicit none

c common
        integer itmin,itmax,itfreq,nrec,idto,itact
	common /timlim/ itmin,itmax,itfreq,nrec,idto,itact
c local
        character*80 line
        integer ianz
	integer iauto
	integer itanf,itend
        real f(10)
c functions
        integer iscan,iround
	real getpar

	iauto = nint(getpar('iauto'))

	if( idto .gt. 0 ) then
          write(6,*) 'itanf,itend,idtout : ',itmin,itmax,idto
	end if

	if( iauto .eq. 0 ) then
          write(6,*) 'Enter time limits : (itmin,itmax,itfreq)'
          read(5,'(a)') line

          ianz = iscan(line,1,f)

          if( ianz .le. 0 .or. ianz .gt. 3 ) then
            return
          else if( ianz .eq. 1 ) then
            itfreq = iround(f(1))
          else if( ianz .eq. 2 ) then
            itmin = iround(f(1))
            itmax = iround(f(2))
          else if( ianz .eq. 3 ) then
            itmin = iround(f(1))
            itmax = iround(f(2))
            itfreq = iround(f(3))
          end if
	else
	  itanf = nint(getpar('itanf'))
	  itend = nint(getpar('itend'))

	  if( itanf .ne. -1 ) itmin = itanf
	  if( itend .ne. -1 ) itmax = itend

	  itfreq = nint(getpar('nout'))
	end if

        if( itfreq .le. 0 ) itfreq = 1

	write(6,*) 'Using time parameters : ',itmin,itmax,itfreq
	write(6,*)

        end

c******************************************************

        function oktime(it)

c is time ok?

        implicit none

c arguments
	logical oktime
        integer it
c common
        integer itmin,itmax,itfreq,nrec,idto,itact
	common /timlim/ itmin,itmax,itfreq,nrec,idto,itact

	integer itold
	save itold
	integer icall
	save icall
	data icall /0/

	if( icall .eq. 0 .or. it .ne. itold ) then !increase only for new time
	  nrec = nrec + 1
	  itold = it
	end if

	icall = icall + 1
        itact = it

c        write(6,*) 'oktime: ',it,nrec,itmin,itmax,itfreq

	if( it .lt. itmin .or. it .gt. itmax ) then
	  oktime = .false.
	else if( mod(nrec,itfreq) .ne. 0 ) then
	  oktime = .false.
	else
	  oktime = .true.
	end if

	end

c******************************************************

        function gettime()

c return current time

        implicit none

c arguments
	integer gettime
c common
        integer itmin,itmax,itfreq,nrec,idto,itact
	common /timlim/ itmin,itmax,itfreq,nrec,idto,itact

	gettime = itact

	end

c******************************************************

