c
c $Id: splitflx.f,v 1.3 2009-04-07 10:43:57 georg Exp $
c
c revision log :
c
c 09.04.1999	ggu	restructured from readext
c 28.09.1999	ggu	reads now all data and then writes it
c
c***************************************************************

	program splitflx

c This routine reads an FLX file and writes the data to single files

	implicit none

	integer nfxdim			!total number of sections
	parameter (nfxdim=20)
	integer noddim			!total number of nodes in sections
	parameter (noddim=1000)
	integer datdim			!total number of data records
	parameter (datdim=400000)

	character*80 line,file
	integer nvers,kfluxm,idfile,nsect
	integer nrec,it,i,nin,kn,in,nout

	integer kflux(noddim)
	real ppp(nfxdim)

	integer itime(datdim)
	real pdata(datdim,nfxdim)

	integer iapini, ideffi, ifileo

c---------------------------------------------------------------
c open files
c---------------------------------------------------------------

        if(iapini(2,1,1,0).eq.0) then
                stop 'error stop : iapini'
        end if

        call defmak('datdir','runnam','.flx',file)
	nin = ifileo(1,file,'unform','old')
        if(nin.le.0) goto 97

c---------------------------------------------------------------
c read and write header information
c---------------------------------------------------------------

        read(nin) idfile,nvers
        read(nin) nsect,kfluxm

        write(6,*) 'file id:                  ',idfile
        write(6,*) 'file version:             ',nvers
        write(6,*) 'total number of sections: ',nsect
        write(6,*) 'nodes in sections:        ',kfluxm

        if( nsect .gt. nfxdim .or. kfluxm .gt. noddim ) then
          stop 'error stop: nfxdim'
        end if
                                                   
        read(nin) (kflux(i),i=1,kfluxm)

	write(6,*) '...reading data'

c---------------------------------------------------------------
c loop over data
c---------------------------------------------------------------

	nrec = 0

        do while( .true. )

          read(nin,end=100) it,nsect,(ppp(i),i=1,nsect)

          nrec = nrec + 1

	  if( nrec .gt. datdim ) then
	    write(6,*) 'Cannot read more than ',datdim,' data records'
	    nrec = nrec - 1
	    goto 100
	  end if

	  if(mod(nrec,100).eq.0) write(6,*) nrec,' data records read'

          itime(nrec) = it
          do i=1,nsect
            pdata(nrec,i) = ppp(i)
          end do

          call check_nan(it,nsect,ppp)
        end do

  100	continue

c---------------------------------------------------------------
c end of data
c---------------------------------------------------------------

	write(6,*)
	write(6,*) 'Total number of data records read : ',nrec
	write(6,*) 'Last time value read: ',itime(nrec)
	write(6,*)

c---------------------------------------------------------------
c writing files
c---------------------------------------------------------------

	do i=1,nsect
	  call wrts(nrec,itime,pdata(1,i),'p',i)
	end do

c---------------------------------------------------------------
c end of routine
c---------------------------------------------------------------

	stop
   89	continue
	stop 'error stop: datdim'
   97	continue
	stop 'Cannot open FLX file'
	end

c*******************************************************************

	subroutine wrts(n,it,data,name,number)

c writes data to file name.number

	implicit none

	integer n
	integer it(1)
	real data(1)
	character*(*) name
	integer number

	integer i
	integer in,nout
	character*80 numlin,file
	integer ialfa

	nout = 1
	in = ialfa(float(number),numlin,-1,-1)
	file = name // '.' // numlin(1:in)

	open(nout,file=file,status='unknown',form='formatted')

	do i=1,n
	  write(nout,*) it(i),data(i)
	end do

	close(nout)

	end

c*******************************************************************

        subroutine check_nan(it,n,p)

        implicit none

        integer it,n
        real p(1)

        integer i

        if( p(1) .gt. 0. ) then
        else if( p(1) .le. 0. ) then
        else
          write(6,*) 'nan found: ',it,n,(p(i),i=1,n)
        end if

        if( it .eq. 216000 ) then
          write(6,*) 'error check...'
          write(6,*) it,n,(p(i),i=1,n)
        end if

        end

c*******************************************************************

