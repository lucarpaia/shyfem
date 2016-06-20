c$Id: seeding.f, v 1.0 2014-03-12 14.53 donata
c
c loading routines
c
c revision log :
c
c 15.06.2004	ggu     written from scratch
c 28.06.2004	dmk     loading factor was wrong: kgs  = 1.e+6
c 28.06.2004	dmk     areaload loops now over state variables
c 14.03.2008	ggu     new routine set_surface_load
c 12.03.2014	dmk	seed of shellfish
c 17.06.2016	dmk	last changes integrated
c
c notes :
c
c needed new routine init_load (sets to 0)
c in setseed_new add to load matrix
c
c*************************************************************

	subroutine setseed_new(eseed)

c sets up eseed which is loading for specified areas
c
c
c variables to be specified:
c
c nareas        total number of areas for which loading is specified
c nodes         total number of nodes used to identify all areas
c karee         node numbers that specify the areas of loading
c iaree         area numbers for the nodes [1-nareas]
c areaseed      total loadings [kg/day] for areas
c
c  the node numbers in karee are external node numbers

	implicit none

        include 'param.h'
	integer nshstate
	parameter(nshstate=3)

	real eseed(nlvdim,nkndim,nshstate)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

        integer ilhkv(1)
        common /ilhkv/ilhkv

	integer nareas
	parameter (nareas=1)
	real volume(nareas)

	real areaseed(nshstate,nareas)
	save areaseed

	integer aree(nkndim)
        save aree

	integer k,l,lmax,i,ia
	real litri,kgs
	real fact,rlfact
	real seed,vol

	real getpar

c        data areaseed/ shellfarm, shellsize, controllo
        data areaseed/
     +   25., 2., 3.
     +  /


        integer nnodes1
        parameter(nnodes1=180)
        integer nodes1(nnodes1)
        save nodes1
        data nodes1 /
     :	31,33,35,36,39,40,42,45,48,51,89,92,110,176,253,296,301,347,348,
     :	350,353,371,394,403,407,409,414,440,441,443,444,450,474,504,509,
     :	525,529,531,536,543,545,554,559,568,570,580,582,586,589,590,595,
     :	598,601,610,614,615,617,622,625,627,630,641,649,653,654,655,657,
     :	658,662,663,675,680,681,708,717,720,780,790,791,808,820,824,851,
     :	855,857,859,888,890,891,900,902,913,926,927,933,934,957,958,991,
     :	996,1006,1008,1009,1012,1015,1018,1023,1025,1043,1053,1071,1076,
     :	1114,1149,1153,1157,1160,1190,1202,1216,1219,1232,1248,1252,1271,
     :	1273,1292,1295,1318,1335,1338,1344,1351,1354,1357,1358,1360,1362,
     :	1365,1366,1368,1369,1370,1371,1373,1374,1377,1378,1379,1382,1384,
     :	1385,1387,1388,1391,1394,1396,1397,1399,1401,1403,1404,1405,1406,
     :	1407,1413,1420,1423,1426,1429,1431,1433,1434,1437,1439,1440,1442,
     :	1444,1477,1481/

        integer icall
        save icall
        data icall / 0 /

c---------------------------------------------------------
c initialization
c---------------------------------------------------------

        if( icall .eq. 0 ) then
          icall = 1

          call seed_init_area(nkn,aree)

          call seed_add_area(1,nnodes1,nodes1,aree)
        end if

     
c---------------------------------------------------------
c compute volumes
c---------------------------------------------------------

        call seed_make_volume(nareas,volume,aree)

c---------------------------------------------------------
c---------------------------------------------------------


        do k=1,nkn
	  ia = aree(k)
          if( ia .gt. 0 ) then
c	    vol = volume(ia)
            lmax = ilhkv(k)
	    do i=1,nshstate
              do l=1,lmax

	        eseed(l,k,i) =   areaseed(i,ia) 
              end do
	    end do
          end if
	end do
c	write(86,*) eseed

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c*************************************************************

c*************************************************************

        subroutine seed_init_area(nkn,aree)

c initializes area array

        implicit none

        integer nkn
        integer aree(1)

        integer k

        do k=1,nkn
          aree(k) = 0
        end do

        end

c*************************************************************

        subroutine seed_add_area(iarea,n,nodes,aree)

c initializes area array

        implicit none

        integer iarea           !number of area
        integer n               !total number of nodes
        integer nodes(1)        !nodes that make up area
        integer aree(1)         !areas of each node

        integer i,k
        logical berror

	call n2int(n,nodes,berror)
	if( berror ) stop 'error stop: load_init_area'

        do i=1,n
	  k = nodes(i)
	  aree(k) = iarea
        end do

        end

c*************************************************************

        subroutine seed_make_volume(nareas,volume,aree)

c makes total volume of areas

        implicit none

        integer nareas
        real volume(1)
        integer aree(1)

        integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
        integer ilhkv(1)
        common /ilhkv/ilhkv

        integer mode,i,k,ia,lmax,l

        real volnode

        mode = +1               !new time level

	do i=1,nareas
	  volume(i) = 0.
	end do

        do k=1,nkn
	  ia = aree(k)
	  if( ia .gt. nareas ) stop 'error stop ia'
	  if( ia .gt. 0 ) then
            lmax = ilhkv(k)
            do l=1,lmax
              volume(ia) = volume(ia) + volnode(l,k,mode)
            end do
          end if
        end do

        end

c*************************************************************


