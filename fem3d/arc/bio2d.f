
c***********************************************
c
c State variables used:
c
c nh3		71	1
c no3		72	2
c opo4		73	3
c phyto		74	4
c cbod		75	5
c do		76	6
c on		77	7
c op		78	8
c zoo		79	9
c***********************************************

	subroutine bio2d(it,idt)

c eco-model cosimo

	implicit none

	include 'param.h'

	integer it	!time in seconds
	integer idt	!time step in seconds

	integer ndim,ntot
	parameter( ndim = 9 , ntot = 1 )

	integer narr
	parameter( narr = 100 )

	real bioarr(narr,nbcdim)	!array containing boundary state
	real bioaux(narr)		!aux array for boundaries

	real e(nkndim,ndim)	!state vector
	real ee(3,neldim,ndim)	!state vector for element
	real eb(nkndim,ndim)	!boundary vector of state vectors

	real eload(3,neldim,ndim)	!vector of loadings

	save e,ee,eb,eload

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,higi
        common /mkonst/ eps1,eps2,pi,flag,high,higi

	character*80 bio2dn(1)
	common /bio2dn/bio2dn

	real soev(3,neldim), snv(nkndim)
	real toev(3,neldim), tnv(nkndim)
        common /soev/soev, /snv/snv   !$$ST
        common /toev/toev, /tnv/tnv   !$$ST

        real uov(1),vov(1)
        common /uov/uov, /vov/vov

	real v1v(1)
	common /v1v/v1v
	real ev(13,1)
	common /ev/ev

	real zeov(3,1)
	common /zeov/zeov
	real zenv(3,1)
	common /zenv/zenv
	real rsv(1)				!??
	common /rsv/rsv
	real hm3v(3,1)
	common /hm3v/hm3v

	integer k,i,l,lmax
	integer ibio
	real t,s,dt
	real eaux(ndim)
	real elaux(ndim)
	real einit(ndim)
	real ebound(ndim)
	real ebriver(ndim)
	integer icall,iunit
	integer itt,idtt,j
	real lux
	real dtt,dttday
	real rkpar,rvpar,azpar,adpar,aapar
	real area,vol
	integer istot,isact
	real oxysat
	real getpar
	integer iround
	integer ieint,ipint

	integer iub,itmcon,idtcon
	integer ie,ii
	real d
	real cbod,nh3,krear,sod
	real vel
	real windspeed,tempair
	real tday,t0,tsec
	real stp
	integer iespecial,inspecial
	save iespecial,inspecial

	integer idox			!index of DO state variable
	parameter ( idox = 6 )

	save rkpar,azpar,istot
	save iub,itmcon,idtcon

	save einit,ebound,icall
c------------------------------------------------------------------
c run1=wasp365 (temperature.con) and wasp365-2 (temperature.in)
c run2=new boundaries and initial conditions wasp365-4.* 
c------------------------------------------------------------------
c				   initial conditions  [mg/l]			
c run2:
        data einit /0.05, 0.4, 0.01, 0.05, 2.,   11.,0.2,0.01,0.015/  
c                  nh3,  no2, opo4, phyto, cbod,do, on, op    zoo   /
c
c       initial boundary conditions at the inlets (bocche) [mg/l]
c run2: insert ebound here as initial condition
	data ebound /0.05,  0.2, 0.01,  0.03,1.47,10.13,0.139,0.01,0.029/
c                 nh3,  no2, opo4,  phyto,cbod, do,   on,  op    zoo   /
c				    rivers		 [mg/l]
c run2 21 sett 2000:
c	data ebriver/4.02, 1.006, 0.5, 0.0,   11.2, 9.37, 0.0, 0.0, 0.0/ 
c                nh3, no2, opo4, phyto, cbod, do,   on,  op   zoo  /
c run7 24 luglio  2001, dati drain globali:
	data ebriver/0.45, 2.46, 0.068, 0.0,    5.6, 9.37, 0.82, 0.144, 0.0/ 
c                nh3, no2, opo4, phyto, cbod, do,   on,  op   zoo  /
	data icall /0/
c------------------------------------------------------------------

c-------------------------------------------------------------------
c initialization
c-------------------------------------------------------------------

	if( icall .le. -1 ) return

	if( icall .eq. 0 ) then
	  ibio = iround(getpar('ibio'))
	  if( ibio .le. 0 ) icall = -1
	  if( icall .le. -1 ) return
	  icall = 1

c	  initialize state variables with einit

	   do k=1,nkn		!loop on nodes
	       do i=1,ndim
	         e(k,i) = einit(i)
	       end do
	       if( einit(idox) .le. 0. ) then	!if DO < 0 -> saturation value
                 t = tnv(k)
                 s = snv(k)
	         e(k,idox) = oxysat(t,s)
	       end if
          end do

c	  distribute nodal value to vertices of element

	  do i=1,ndim
	    call cp2elem(e(1,i),ee(1,1,i))
	  end do
c
c	  set loadings in the interal areas
c
	  call setload(eload)
c
c	  set boundary conditions for all state variables
c
	  do i=1,ndim
	   do k=1,nkn
	    if( rsv(k) .eq. flag ) then
	      eb(k,i) = flag
	    else
	      eb(k,i) = ebound(i)
	    end if
	   end do
	  end do

	  call bnds_init(bio2dn,2,ndim,narr,bioarr,ebound)	!new_section

c	  initialize eutro

	  call eutroini

c	  parameters for transport/diffusion resolution

          rkpar=getpar('chpar')
          call getaz(azpar)
          istot=iround(getpar('istot'))         !$$istot

c	  meteorological forcing (HACK -> FIXME)

	  tempair = 22.
	  windspeed = 3.			!FIXME
	  call wmeteo(tempair,windspeed)

c	  initialize output 

	  iub = 55
          itmcon = iround(getpar('itmcon'))
          idtcon = iround(getpar('idtcon'))

          call confop(iub,itmcon,idtcon,ndim,'bio')

	  write(6,*) 'bio2d model initialized...'

c	  special element and node

	  iespecial = 7805
	  inspecial = 4333
	  iespecial = ieint(iespecial)
	  inspecial = ipint(inspecial)
	end if

c-------------------------------------------------------------------
c normal call
c-------------------------------------------------------------------

c	-------------------------------------------------------------------
c	time management
c	-------------------------------------------------------------------

	t0 = 0.
	dt = idt
	idtt = idt/ntot
	dtt = dt / ntot
	dttday = dtt / 86400.

	tsec = it
	tday = it / 86400. + t0		!time in days, FEM 0 is day t0

c	leggo la temperatura quando non voglio usare quella calcolata dal
c	modello idrodinamico
c	devo anche commentare toev=t... qui sotto
c
c       call rdtemp(tday,stp)           !read sea temperature

c	-------------------------------------------------------------------
c	boundary conditions	!new_section
c	-------------------------------------------------------------------

	call bnds_set(tsec,narr,bioarr,bioaux)
	call bnds_set_global(narr,bioarr,nkndim,1,ndim,eb,e)

	do i=1,ndim
          call n2ebar(e(1,i),ee(1,1,i))
	end do

c	call binlet(tday,ndim,ebound)	!read inlet concentrations

c	-------------------------------------------------------------------
c	loop on elements for biological reactor
c	-------------------------------------------------------------------

	do ie=1,nel		!loop on elements

	  area = 4. * ev(10,ie)

	  do ii=1,3

c	    t=stp
	    t = toev(ii,ie)			!temperature
	    s = soev(ii,ie)			!salinity
	    d = hm3v(ii,ie) + zeov(ii,ie)	!depth
	    vol = area * d			!volume of box

	    vel = sqrt( (uov(ie)/d)**2 + (vov(ie)/d)**2 )

	    do i=1,ndim
	      eaux(i) = ee(ii,ie,i)
	      elaux(i) = eload(ii,ie,i)
	    end do

	    if( ii .eq. 1 ) then
	      if( ie .eq. 1 ) then
	      end if
	    end if	
	    if( ii .eq. 2 ) then
	      if( ie .eq. 1 ) then

	      end if
	    end if

	    do j=1,ntot
	      itt = it - idt + j*idtt
	      call eutro0d(tday,dttday,vol,d,vel,t,s,eaux,elaux)
	    end do

	    do i=1,ndim
	      ee(ii,ie,i) = eaux(i)
	    end do

	  end do
	end do

c	-------------------------------------------------------------------
c	advection and diffusion
c	-------------------------------------------------------------------

	do i=1,ndim

          do isact=1,istot                !$$istot
	    call conz2d(e(1,i),ee(1,1,i),v1v,dt,rkpar,azpar,istot,isact)
	    call conzbc(e(1,i),ee(1,1,i),v1v,eb(1,i),flag,azpar)
          end do

	end do

c	-------------------------------------------------------------------
c	boundary (rivers)	!new_section
c	-------------------------------------------------------------------

c	call briver (tday,ndim,ebriver) !inserisce una conc. variabile 
c	call biobnd0(ndim,e,ebriver,dt)
c	do i=1,ndim
c          call n2ebar(e(1,i),ee(1,1,i))
c	end do

c	-------------------------------------------------------------------
c	write of results (file BIO)
c	-------------------------------------------------------------------

	do i=1,ndim
          call confil(iub,itmcon,idtcon,70+i,1,e(1,i))
	end do

c	-------------------------------------------------------------------
c	end of routine
c	-------------------------------------------------------------------

	end

c*************************************************************

	subroutine writee(iunit,it,k,l,e,t,s,nlvdim,nkndim,ndim)

c formatted write for debug

	implicit none

	integer iunit
	integer it,k,l
	integer nlvdim,nkndim,ndim
	real e(nlvdim,nkndim,ndim)
	real t(nlvdim,nkndim)
	real s(nlvdim,nkndim)

	integer i

	write(iunit,'(i10,11f12.4)') it,
     +			(e(l,k,i),i=1,ndim),
     +			t(l,k),
     +			s(l,k)

	end

c*************************************************************

	subroutine biobnd0(ndim,e,ebound,dt)

c handles open boundary conditions
c
c since this is called AFTER the T/D step, the first time step
c is missed -> nothing is injected in the first time step
c there may be also some inconsistencies, since the water volume
c used is the one of the old time step, but the actual
c water volume injected is determined only afterwards

        implicit none

c parameter
        include 'param.h'
c arguments
	integer ndim
        real e(nkndim,ndim)
        real ebound(ndim)
	real dt
c local
        integer ibc,nbc,ibtyp,levmax,n,j,kn,i
	real rw,vol
c functions
	integer nbnds,itybnd,levbnd,kbnds,nkbnds
	real zvbnds

	nbc = nbnds()

        do ibc=1,nbc

	  ibtyp = itybnd(ibc)
          levmax = levbnd(ibc)

	  if( ibtyp .eq. 3 ) then
            rw = zvbnds(ibc)
	    vol = rw * dt

	    n = nkbnds(ibc)

            do j=1,n

             kn = kbnds(ibc,j)

	     do i=1,ndim
               call volno0(kn,levmax,1,e(1,i),vol,ebound(i))
	     end do

	    end do
	  end if

	end do

	end

c*************************************************************

        subroutine biobnd(conz,rbound,rsv)

        implicit none

c parameter
        include 'param.h'
c arguments
        real conz(nlvdim,1)
        real rbound
	real rsv(1)
c common
        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        real eps1,eps2,pi,flag,high,hihi
        integer nlvdi,nlv
        integer ilhkv(1)
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /mkonst/ eps1,eps2,pi,flag,high,hihi
        common /level/ nlvdi,nlv
        common /ilhkv/ilhkv
c local
        integer k,l,lmax

        if(nlvdim.ne.nlvdi)stop'error stop : level dimension in biobnd'

        do k=1,nkn
          if( rsv(k) .ne. flag ) then
            lmax=ilhkv(k)
            do l=1,lmax
                conz(l,k) = rbound
            end do
          end if
        end do

        end

c*************************************************************

	subroutine cp2elem(e,ee)

	implicit none

	real e(1)
	real ee(3,1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer nen3v(3,1)
        common /nen3v/nen3v

	integer ie,ii,k

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ee(ii,ie) = e(k)
	  end do
	end do

	end

c*************************************************************

	subroutine setload(eload)

c sets up eload which is loading for specified areas
c
c the computed loadings in eload are in [g/(m**3 day)] == [mg/(l day)]
c the specified loadings in areaload are in [kg/day]
c
c variables to be specified:
c
c nimmis        total number of areas for which loading is specified
c nodes         total number of nodes used to identify all areas
c karee         node numbers that specify the areas of loading
c iaree         area numbers for the nodes [1-nimmis]
c areaload      total loadings [kg/day] for areas
c
c the node numbers in karee are external node numbers

	implicit none

        include 'param.h'
	integer ndim
	parameter(ndim=9)

	real eload(3,neldim,ndim)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real ev(13,1)
	common /ev/ev
	integer nen3v(3,1)
	common /nen3v/nen3v
	real hm3v(3,1)
	common /hm3v/hm3v

	integer nimmis
	parameter (nimmis=5)
	real volaux(nimmis)
	real areaload(nimmis,ndim)

	save volaux,areaload

	integer aree(nkndim)

	integer nodes
	parameter(nodes=152)
	integer karee(nodes)
	integer iaree(nodes)
	save karee,iaree

	logical berror
	integer k,ie,ii,ia,i
	integer itype
	real area
	real litri,kgs
	real fact,rlfact
	real hdep,vol,load

	real getpar

c loading is kg/day
c
c	1=venezia-giudecca, 2=Murano-S. Erasmo 3=Burano
c	4=Fusina, 5=Campalto 
c
c	loading for areas [kg/day]
c 	
c	21 sett 2000

c	donata inserire cavallino

c carichi ripartiti secondo % franco

c 	data areaload /
c     + 1055.,70.,58.,792.,348. 			!nh3
c     + ,264.,17.6,14.61,198,87			!no3
c     +	,142.4,9.5,6.33,183,42			!opo4
c     + ,0.,0.,0.,0.,0.				!phyto
c     +	,7907.,518.,353.,2193.,915.		!cbod
c     +	,0.,0.,0.,0.,0.				!do
c     +	,0.,0.,0.,0.,0.				!on
c     +	,0.,0.,0.,0.,0.				!op
c     +	,0.,0.,0.,0.,0./			!zoo

c carichi ripartiti assegnando 50% alla forma ridotta dell'azoto e 50% alla
c forma ossidata

 	data areaload /
     +  659.,44.,37.,495.,217. 			!nh3
     + , 659.,44.,37.,495.,217. 		 	!no3
     +	,125.4,8.36,5.57,183.,61.		!opo4
     + ,0.,0.,0.,0.,0.				!phyto
     +	,7907.,518.,353.,2193.,915.		!cbod
     +	,0.,0.,0.,0.,0.				!do
     +	,0.,0.,0.,0.,0.				!on
     +	,0.,0.,0.,0.,0.				!op
     +	,0.,0.,0.,0.,0./			!zoo

	data karee /
     : 2112,2113,2863,2858,2114,2508,2507,2506,2505,2504,
     : 2503,2502,2501,2500,2499,2598,2497,2496,2495,2494,
     : 2493,2492,2491,2490,2489,2488,2487,2486,2485,2484,
     : 2447,2555,4235,2448,2449,2133,2132,2131,2130,2129,
     : 2128,2125,2124,2123,2122,4211,2109,2379,2108,2111,
     : 2112,2139,2138,2137,2136,2135,2134,2450,2451,2452,
     : 2453,2454,2455,2144,2143,2142,2141,2140,
     :
     : 3001,3105,2999,2998,2997,3096,2995,2994,3007,3084,
     : 3006,3005,3004,3003,3002,3070,2989,2988,3013,3014,
     : 3015,3016,3017,3018,3019,3020,3301,3300,3299,4357,
     : 3296,3295,3294,4358,3293,3366,2537,2536,2535,2534,
     : 2533,
     
     : 3041,3040,3039,3038,3192,3188,3029,3030,4210,3309,
     : 3308,3307,3306,3305,3304,3023,3024,3025,3026,3032,
     : 3044,3043,3042,3047,3048,3049,3050,3806,3805,3844,
     : 3804,3319,3318,3317,3316,3315,3314,3313,3312,3311,
     : 3310,
     :
     : 2172,		! fusina 2 maggio 2000 
     :
     : 2954		!campalto 2 maggio 2000
     : 
     :  /	!node in area
c
	data iaree /
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,
     : 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     : 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     : 2,
     : 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     : 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     : 3,
     : 4,		!4,4,
     : 5
     :  /	!type of area
c
	litri = 1000.	!litri in m*3
	kgs  = 10.e+6	!mg in kg

c	rlfact = getpar('flbio')
	rlfact = 1.

	fact = rlfact*kgs/litri         ! [kg/m**3] -> [mg/l]
     
c extern to intern

	call n2int(nodes,karee,berror)

	if( berror) stop 'error stop: loading'

c intialize

	do i=1,nimmis
	  volaux(i) = 0.
	end do

	do k=1,nkn
	  aree(k) = 0
	end do

	do i=1,nodes
	  k = karee(i)
	  itype = iaree(i)
	  aree(k) = itype
	end do

c compute total volume for all areas given -> store in volaux

	do ie=1,nel
	  area = 4. * ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hdep = hm3v(ii,ie)
	    ia = aree(k)
	    if( ia .gt. nimmis ) stop 'error stop ia'
	    if( ia .gt. 0 ) then
		volaux(ia) = volaux(ia) + area * hdep
	    end if
	  end do
	end do

c compute and set loading in eload [g/(m**3 day)] == [mg/(l day]

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ia = aree(k)
	    vol = volaux(ia)
	    do i=1,ndim
	      load = fact * areaload(ia,i) / vol
	      if( ia .le. 0 ) load = 0.
	      eload(ii,ie,i) = load
	    end do
	  end do
	end do

	end

c*************************************************************
c*************************************************************
c*************************************************************

	subroutine binlet(t,neb,ebound)

c gets boundary conditions from file

	implicit none

	real t		!time [day]
	integer neb	!total number of boundary conditions
	real ebound(1)	!boundary values

	integer nstat,nyear
	parameter( nstat = 9 , nyear = 365 )

	real ebfile(nstat,nyear)
	save ebfile

	integer n,i
	real tdummy
	integer it2n

	integer icall,nold
	save icall,nold
	data icall,nold / 0 , 0 /

	if( neb .ne. nstat ) stop 'error stop binlet: neb <> nstat'

	if( icall .eq. 0 ) then
	  open(2,file='input/ebound.dat',status='old',form='formatted')
	  do n=1,nyear
	    read(2,*) tdummy,(ebfile(i,n),i=1,nstat)
	    if( nint(tdummy) .ne. n ) stop 'error stop binlet: tdummy'
	  end do
	  close(2)
	  write(6,*) 'eutro boundary file read: ebound.dat'
	  icall = 1
	end if

	n = it2n(t)		!day in year

	do i=1,neb
	  ebound(i) = ebfile(i,n)
	end do

	if( n .ne. nold ) then
	  nold = n
	  write(6,*) 'binlet: New boundary values read for day ',n,t
	  write(6,*) (ebound(i),i=1,neb)
	end if

	end

c*************************************************************

	subroutine briver (t,neb,ebriver)


c gets river boundary conditions from file

	implicit none

	real t		!time [day]
	integer neb	!total number of boundary conditions
	real ebriver(1)	!boundary values

	integer nstat,nyear
	parameter( nstat = 9 , nyear = 365 )

	real ebrfile(nstat,nyear)
	save ebrfile

	integer n,i
	real tdummy
	integer it2n

	integer icall,nold
	save icall,nold
	data icall,nold / 0 , 0 /

	if( neb .ne. nstat ) stop 'error stop briver: neb <> nstat'

	if( icall .eq. 0 ) then
	  open(2,file='input/ebriver.dat',status='old',form='formatted')
	  do n=1,nyear
	    read(2,*) tdummy,(ebrfile(i,n),i=1,nstat)
	    if( nint(tdummy) .ne. n ) stop 'error stop briver: tdummy'
	  end do
	  close(2)
	  write(6,*) 'eutro boundary file read: ebriver.dat'
	  icall = 1
	end if

	n = it2n(t)		!day in year

	do i=1,neb
	  ebriver(i) = ebrfile(i,n)
	end do

	if( n .ne. nold ) then
	  nold = n
	  write(6,*) 'briver: New boundary values read for day ',n,t
	  write(6,*) (ebriver(i),i=1,neb)
	end if

	end

c*************************************************************


