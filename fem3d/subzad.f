
!--------------------------------------------------------------------------
!
!    Copyright (C) 2011-2017,2019  Georg Umgiesser
!    Copyright (C) 2022  Luca Arpaia
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c routines for handling z-adaptive layers
c
c revision log :
c
c 28.03.2022    lrp     introduce z-adaptive layers
c
c notes:
c this file is used also in:      
c      
c	compute_zadaptive_info (subele.f)
c	get_zadapt_info (newexpl.f)
c	set_zadapt_info (subele.f)
c	get_zadaptivelayer_thickness (subele.f)
c	init_zadaptation (shyfem.f)
c	zadaptation (shyfem.f)
c	copy_zadaptation (shyfem.f)
c	get_zadaptive_weights (new3di.f,newcon_omp.f)
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module zadapt
!==================================================================

	implicit none

        integer, save, allocatable :: nadapt_com(:,:) !number of adaptive layers
        real   , save, allocatable :: hadapt_com(:,:) !closing depth of adaptive layers
	integer, save, allocatable :: iskremap(:,:) !remap yes/no
	integer, save, allocatable :: iseremap(:,:) !remap yes/no
	integer, save, allocatable :: iskbathywall(:,:) !is layer on a resolved bathymetry wall yes/no
	real   , save :: nzadapt_com		      !minimum number of adaptive layers
	real, save :: rzmov_com		      !parameter for moving surface layers
	real, save :: rztop_com		      !parameter for removing surface layers

!==================================================================
        end module zadapt
!==================================================================

c******************************************************************

        subroutine get_rzmov_info(rzmov)

        use zadapt

        implicit none

        real rzmov

        rzmov = rzmov_com

        end

c******************************************************************

        subroutine set_rzmov_info(rzmov)

        use zadapt

        implicit none

        real rzmov

        rzmov_com = rzmov

        end

c******************************************************************

        subroutine get_rztop_info(rztop)

        use zadapt

        implicit none

        real rztop

        rztop = rztop_com

        end

c******************************************************************

        subroutine set_rztop_info(rztop)

        use zadapt

        implicit none

        real rztop

        rztop_com = rztop

        end

c******************************************************************

	subroutine compute_rzpar_info(nlv,nzadapt,hlv,rzmov,rztop)

	use zadapt

	implicit none

        integer nlv             !total number of layers
	integer nzadapt		!number of surface moving layers
        real hlv(nlv)           !layer structure
	real rzmov 		!parameter for moving surface layers (return)
	real rztop		!parameter for removing surface layers (return)

	real maxz
	
	maxz = 0.0				!estimate of max water level

	if (nzadapt .le. 1) then		!z-layers
	  rzmov = 0.125
	  rztop = 0.200
	else if (nzadapt .ge. nlv) then		!z-star
	  rzmov =  10000.
	  rztop = -10000.
	else 					!z + z-star
          rzmov = (maxz+hlv(nzadapt-1)) / (hlv(nzadapt)-hlv(nzadapt-1))
	  rztop = -10000.
	end if

	end 

c******************************************************************

        subroutine get_nzadapt_info(nzadapt)

        use zadapt

        implicit none

        integer nzadapt

        nzadapt = nzadapt_com

        end

c******************************************************************

        subroutine set_nzadapt_info(nzadapt)

        use zadapt

        implicit none

        integer nzadapt

        nzadapt_com = nzadapt

        end

c******************************************************************

        subroutine get_zadapt_info(ie,nadapt,hadapt)

        use zadapt

        implicit none

	integer ie
        integer nadapt(4)
        real hadapt(4)

	!call check_sigma_initialized  !lrp do some check here

        nadapt = nadapt_com(:,ie)
        hadapt = hadapt_com(:,ie)

        end

c******************************************************************

        subroutine set_zadapt_info(ie,nadapt,hadapt)

        use zadapt

        implicit none

        integer ie
        integer nadapt(4)
        real hadapt(4)

        nadapt_com(:,ie) = nadapt
        hadapt_com(:,ie) = hadapt

        end

c******************************************************************

        subroutine compute_zadaptive_info(ie,nlv,lmin,lmax,hlv,z,htot,
     +					  nadapt,ladapt,hadapt,hdl)

c returns adaptive layers info. Info is computed by node of element:
c number of adaptive layers       by node (3) + by ele (1)
c lowest index of adaptive layer  by node (3) + by ele (1)
c closing depth of adaptive layer by node (3) + by ele (1)
c coefficients of adaptive layers 

        implicit none

        integer ie              !element index
        integer nlv             !total number of layers
        integer lmin(3)         !top layer index	
	integer lmax		!bottom layer index
	real hlv(nlv)           !layer structure
        real z(3)               !water level
	real htot(3)		!water depth	
        integer nadapt(4)       !total number of adaptive layers (return)
	integer ladapt(4)       !lowest index of adaptive layers (return)
        real hadapt(4)          !closing depth of adaptive layers(return)
	real hdl(nlv,3)         !coefficient (return)

	integer l,ii,levmax,lmine,nsigma,nlev
	real getpar,r,hsigma,htop,hbot
	real den,check
        logical bsigma

	nadapt = 0         !no adaptation -> all to zero
	ladapt = 0	   !no adaptation -> all to zero
	hadapt = 0	   !no adaptation -> all to zero	
	hdl = 0.           !no adaptation -> all to zero

	call get_rzmov_info(r)
        call get_sigma_info(nlev,nsigma,hsigma)
        bsigma = nsigma .gt. 0

	lmine = maxval(lmin)

        if( .not.bsigma ) then

c---------------------------------------------------------
c loop over nodes: adaptation is node-driven
c---------------------------------------------------------

	do ii=1,3
	  levmax = 0       !no adapation -> all to zero
	  do l=lmin(ii),lmax-1 !-1 to skip bottom layer

c---------------------------------------------------------
c lowest index of adaptive deforming layers
c---------------------------------------------------------	  
	  
            if(z(ii).le.(-hlv(l)+r*(hlv(l)-hlv(l-1)))) then 
              levmax = l+1

c---------------------------------------------------------
c no deformation
c--------------------------------------------------------- 

	    else  
	      exit 
	    end if  
          end do

c---------------------------------------------------------
c compute nadapt, ladapt, hadapt 
c---------------------------------------------------------  	  

	  hadapt(ii) = hlv(levmax) 
	  nadapt(ii) = max(0,levmax-lmin(ii)+1) !+1 (min 2 adaptive layer)
	  ladapt(ii) = levmax

c---------------------------------------------------------
c compute hdl: different strategy tested
c--------------------------------------------------------- 

    	  if (nadapt(ii).gt.0) then	  
c	  den = (nsigma(ii)-1.)+r		!freezed
          den = hadapt(ii)-hlv(lmin(ii)-1) !zstar
	  if (ladapt(ii).eq.lmax) den = htot(ii)-hlv(lmin(ii)-1)
          do l=lmin(ii),ladapt(ii)
c	    hdl(l,ii) = - 1. / den		!freezed		
c           hdl(l,ii) = - 1. / nsigma(ii)       !constant
	    htop = hlv(l-1)
	    hbot = hlv(l)
	    if (l.eq.lmax) hbot = htot(ii)
            hdl(l,ii) = (htop-hbot)/den   	!zstar
	  end do
c         hdl(lmin(ii),ii) = - r / den          !freezed
	  check = 0.
	  do l=lmin(ii),ladapt(ii)
            check = check + hdl(l,ii)
          end do		    
	  if (abs(check+1.).gt.1e-5) then
	    write(6,*) 'error computing layer thickness'
	    write(6,*) 'you are using z-adaptive levels'
            write(6,*) 'with local layer deformation but'
	    write(6,*) 'the weights does not sum to -1 ', check
	    write(6,*) 'lmin,ladapt,lmax ', lmin(ii),ladapt(ii),lmax
            write(6,*) 'hadapt,htot ', hadapt(ii),htot(ii)	    
	    stop 'error stop in compute_zadaptive_info'	    
	  end if
	  end if

c---------------------------------------------------------
c compute element info: number of adaptive layers in ele
c---------------------------------------------------------

	  !flag adaptive layers in element with non-conformal edge:
	  !free-surface must span all layers greater then lmin
	  if (ladapt(ii).gt.lmine) then
	    ladapt(4) = max(ladapt(ii),ladapt(4))
	    nadapt(4) = ladapt(4)-lmine+1
	    hadapt(4) = max(hadapt(ii),hadapt(4))	  
	  end if
	end do

	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************

        subroutine get_zadaptivelayer_thickness(l,lmin,hzad,levz,hldv,
     +                          		hdl,hlv,z,hnod)

c returns nodal layer thickness for a z-adaptive layer

       	implicit none

	integer l            !layer index
        integer lmin         !top layer index	
	real hzad	     !closing depth of adaptive layer 
	integer levz	     !lowest index of adaptive layers
	real hldv	     !coefficient of adaptive layer
        real hdl             !z-layer thickness
	real hlv	     !layer structure
	real z		     !water level
	real hnod	     !nodal layer thickness computed (return)

c---------------------------------------------------------
c adaptive z-layer found
c---------------------------------------------------------

        if (l.le.levz) then
	  hnod = - hzad * hldv

c---------------------------------------------------------
c fixed z-layer found
c---------------------------------------------------------

	else
          hnod = hdl
 	  if( l .eq. lmin ) then
            hnod = hnod + hlv + z
          end if
        end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------	

	end

c*****************************************************************

        subroutine get_zadaptive_weights(ie,jlevelmin,weight,mode)

c returns weights to split transport across sublevels for
c non-conformal element. Non-conformal element are those:
c
c       x-------------x----------------x
c       !     U_1     !                !      
c       x-------------x                ! 
c       !             !       U_2      !
c       !             !                !
c       !     U_2     !                !
c       x-------------x----------------x
c       
c ele       i-1/2            i+1/2
c type    conformal      non-conformal
c jlhv        1               2
c jlhev 1             1                2
c
c Weights are needed for advection schemes on non-conformal elements:
c element is splitted connecting vertices of non-conformal edges 
c to half edge of conformal edges:
c
c 	x-------------x----------------x
c 	!     U_1     !     U_2*w_1    !      
c 	x-------------x--------        ! 
c 	!             !        --------!
c 	!             !                !
c 	!     U_2     !      U_2*w_2   !
c 	x-------------x----------------x
c
c TODO hdknv still global ... pass as argument 

	use basin
	use levels
        use mod_layer_thickness

	implicit none

	integer :: ie		!element index
	integer :: jlevelmin 	!minimum top-layer nodal index (return)
	double precision,dimension(nlvdi) :: weight !weights (return)
	integer mode		!mode

	integer jlevel
	integer jlevelk(3)
	integer lmiss,k,ii,n,numOfLev
	double precision w
	real depnode
	logical isenonconf

c---------------------------------------------------------
c mode
c---------------------------------------------------------    

	if( mode .gt. 0 ) then
	  jlevel = jlhv(ie)
	  jlevelk = jlhev(:,ie)
	else if( mode .lt. 0 ) then
          jlevel = jlhov(ie)
          jlevelk = jlheov(:,ie)		
	else
	  write(6,*) 'mode = 0 not implemented'
	  stop 'error stop'
	end if

	n = 3
        jlevelmin = huge(1)
        do ii=1,3
          jlevelmin=min(jlevelmin,jlevelk(ii))
        end do

c---------------------------------------------------------
c is element non-conformal?
c---------------------------------------------------------	

        isenonconf = .false.
        if (jlevel.gt.jlevelmin) isenonconf = .true.	

	if (isenonconf) then

        w = 0.

c---------------------------------------------------------
c loop over sub-layers
c---------------------------------------------------------

        do lmiss=jlevel,jlevelmin,-1

          weight(lmiss) = 0.

c---------------------------------------------------------
c loop over node -> compute element sub-layer depth
c---------------------------------------------------------

	  do ii=1,n

	    k=nen3v(ii,ie)

	    if (lmiss.gt.jlevelk(ii)) then
              weight(lmiss) = weight(lmiss) + depnode(lmiss,k,mode)
            else
              numOfLev = real(jlevelk(ii)-jlevelmin+1)
	      weight(lmiss) = weight(lmiss) + 
     + 				depnode(jlevelk(ii),k,mode)/numOfLev
            end if

	  end do

          w = w + weight(lmiss)

	end do

c---------------------------------------------------------
c compute sub-layer weight
c---------------------------------------------------------

        do lmiss=jlevel,jlevelmin,-1
          if (w.gt.0) then		  
            weight(lmiss) = weight(lmiss)/w
	  end if 
        end do

	end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end	    

c*****************************************************************

	subroutine set_jlhv

c sets elemental top layer index: jlhv and jlhvo 
c only needs zenv and hlv
c jlhv(ie)	new layer index of the upper existing layer 
c jlhvo(ie)	old layer index of the upper existing layer
c rlhv(ie)      lower layer index where insertion/removal
c               operations occurred: needed for remap operations
c e.g. jlhv(ie)=2 means first layer has been removed for element ie

        use zadapt
	use levels
	use basin
	use mod_hydro

	implicit none

c local
	logical bsigma
	integer ie,ii,l,lmax,nlev,nsigma
	real hsigma
	real zmin
        real getpar,r

	lmax=0
	call get_rztop_info(r)

	call get_sigma_info(nlev,nsigma,hsigma)
	bsigma = nsigma .gt. 0

	!rlhv = 0		! (moved to copy_zadaptation)
        iseremap = 1            ! initialize 1: no config change	

c---------------------------------------------------------
c loop over element
c---------------------------------------------------------

	do ie=1,nel

c---------------------------------------------------------
c set jlhv
c---------------------------------------------------------

	  zmin = 1000000.
	  do ii=1,3
	    zmin = min(zmin, zenv(ii,ie))
	  end do

	  if( bsigma ) then	! sigma levels always start from 1
	    l = 1
	  else
	    do l=1,nlv          ! a threshold is used to avoid very small
				! layers
	      if((-hlv(l)+r*(hlv(l)-hlv(l-1))).le.zmin) exit
              !if((-hlv(l)+hzoff).le.zmin) exit
	    end do
	  end if

	  jlhv(ie)=min(l,ilhv(ie)) !safety min: jlhv>=ilhv
	  lmax = max(lmax,jlhv(ie))

c---------------------------------------------------------
c set rlhv
c---------------------------------------------------------

	  if (jlhv(ie).gt.jlhov(ie)) then	!top layer removed	
            rlhv(ie) = jlhv(ie) 
	  else if (jlhv(ie).lt.jlhov(ie)) then	!top layer inserted
	    rlhv(ie) = jlhov(ie)
          end if
	  iseremap(1:rlhv(ie),ie) = 0

	end do

c	nlv = lmax
c	write(6,*) 'finished setting jlhv and nlv'
c	write(6,*) 'nlv,lmax: ',nlv,lmax
c	write(6,'(5g14.6)') (hlv(l),l=1,nlv)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end	

c*****************************************************************

        subroutine set_jlhkv
		
c set nodal index for top layer jlhkv,jlhev and jlhkov,jlheov array - 
c set nodal index for remap operations
c only needs zenv
c a double index is needed (by node and elem) for wet/dry interfaces:
c top layer depends on the free surface and at wet/dry interfaces 
c the free-surface is discontinous -> element storing is also required
c jlhkv(k)     new layer index of the upper existing node (node)
c jlhev(ii,ie) new layer index of the upper existing node (elem)
c jlhkov(k)    old layer index of the upper existing node   
c jlhkev(ii,ie)old layer index of the upper existing node
c rlhkv(k)     lowest layer index where insertion/removal
c              operations occurred: needed for remap operations
c 
c e.g. jlhkv(k)=jlhkv(ii,ie)=2 means first layer has been removed 
c      for node k <-> (ii,ie) but only if k is not wet/dry
c

        use zadapt
        use levels
        use basin
        use shympi
	use mod_hydro
	use mod_layer_thickness

        implicit none

	logical bsigma
        integer ie,ii,k,l,nlev,nsigma
	integer lminnv,lminov,ladaptov
        real hsigma
	real getpar,r
        integer nadapt(4)
        real hadapt(4)	

	call get_rztop_info(r)
        call get_sigma_info(nlev,nsigma,hsigma)
        bsigma = nsigma .gt. 0

        !rlhkv = 0		  !(moved to copy_zadaptation)
	iskremap = 1		  !initialize 1: no config change

        jlhkv=huge(1)

c---------------------------------------------------------
c loop over element
c---------------------------------------------------------

        do ie=1,nel

	call get_zadapt_info(ie,nadapt,hadapt)

	do ii=1,3
	    k=nen3v(ii,ie)

c---------------------------------------------------------
c set jlhev, jlhkv
c---------------------------------------------------------

            if( bsigma ) then     !sigma levels always start from 1
              l = 1
            else
              do l=1,nlv
                if((-hlv(l)+r*(hlv(l)-hlv(l-1))).lt.zenv(ii,ie)) exit
                !if((-hlv(l)+hzoff).lt.zenv(ii,ie)) exit		
              end do
            end if

	    jlhev(ii,ie)=min(l,jlhv(ie)) !safety min: jlhev>=jlhv, jlhkv>=ilhkv

	    if (jlhev(ii,ie).lt.jlhkv(k)) then 
	      jlhkv(k)=min(jlhev(ii,ie),ilhkv(k))
	    end if

c---------------------------------------------------------
c set rlhkv
c---------------------------------------------------------  	      
	      
            if (jlhev(ii,ie).gt.jlheov(ii,ie)) then       !top layer removed 
              rlhkv(k) = max(rlhkv(k),jlhev(ii,ie))	  !safety max for w/d nodes 
            else if (jlhev(ii,ie).lt.jlheov(ii,ie)) then  !top layer inserted
              rlhkv(k) = max(rlhkv(k),jlheov(ii,ie))
            end if  
            iskremap(1:rlhkv(k),k) = 0

c---------------------------------------------------------
c update nadapt
c--------------------------------------------------------- 

	    lminnv = jlhev(ii,ie)
	    lminov = jlheov(ii,ie)
	    if (lminnv.gt.lminov) then			  !update only with removal
	      ladaptov = nadapt(ii) + lminov - 1	
	      nadapt(ii) = ladaptov - lminnv
	    end if
          end do
	 
          lminnv = jlhv(ie)
          lminov = jlhov(ie)	  
          if (lminnv.gt.lminov) then                      !update only with removal	  
	    ladaptov = nadapt(4) + lminov - 1
	    nadapt(4) = ladaptov - lminnv
	  end if		  
	  call set_zadapt_info(ie,nadapt,hadapt)

  	end do

        if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange jlhkv - min')
          call shympi_exchange_2d_nodes_min(jlhkv)
        else
          call shympi_exchange_2d_node(jlhkv)
        end if

        if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange ilhkv - max')
          call shympi_exchange_2d_nodes_max(rlhkv)
        else
          call shympi_exchange_2d_node(rlhkv)
        end if

        !if( shympi_partition_on_elements() ) then !FIXME: parallel implementation
        !  call shympi_exchange_3d_nodes_min(iskremap)
        !else
        !  call shympi_exchange_3d_node(iskremap)
        !end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

        end

c*****************************************************************

        subroutine set_jwlhkv
		
c set nodal index for top layer jwlhkv array - 
c only needs zenv
c at wet/dry nodes (boundary nodes) jlhkv is computed only "on wet side".
c This is because many shyfem subroutines work with non-negative layers
c and, at wet/dry nodes with z-adaptation, it may happen that zero-depth 
c layers separate layers on "wet side" from layers on "dry side": 
c l=jwlhkv,ilhkv 
c means that, at boundary nodes, you loop only layers on "wet side"

        use levels
        use basin
        use shympi
	use mod_geom_dynamic

        implicit none

        integer ie,ii,k
	logical isein,iseout,iskin

	isein(ie)  = iwegv(ie).eq. 0
        iseout(ie) = iwegv(ie).gt. 0
	iskin(k)   = inodv(k) .ne.-1

        jwlhkv=huge(1)

c---------------------------------------------------------
c loop over element
c---------------------------------------------------------

        do ie=1,nel

	  do ii=1,3
	    k=nen3v(ii,ie)

c---------------------------------------------------------
c set jwlhkv
c---------------------------------------------------------

            if (isein(ie).or.                            !boundary nodes:
     +          (iseout(ie).and.iskin(k))) then          !jwlhkv computed on wet side 
              if (jlhev(ii,ie).lt.jwlhkv(k)) then
                jwlhkv(k)=min(jlhev(ii,ie),ilhkv(k))
              end if
            end if

          end do
  	end do

        if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange jlhkv - min')
          call shympi_exchange_2d_nodes_min(jwlhkv)
        else
          call shympi_exchange_2d_node(jwlhkv)
        end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

        end

c*****************************************************************

        subroutine check_bathywall

c returns a flag for layer belonging to resolved bathymetry vertical 
c slope in case of z-layer modelling:
c				
c				  flag	
c       -------x------ 
c			layer 1	=  0
c 	.......x------
c	       !        layer 2 =  1
c              x------
c	       !        layer 3 =  0 
c	       x......				
c       
c	ilhv=1   ilhv=3
c

        use zadapt
        use levels
        use basin

	implicit none

	integer :: ie,ii,k,lmin

        iskbathywall = 1

        do ie=1,nel
          do ii=1,3
            k = nen3v(ii,ie)
            iskbathywall(ilhv(ie),k) = 0 !flag for bottom layers
          end do
        end do

        do k=1,nkn
          lmin = ilmkv(k)
          iskbathywall(1:lmin,k) = 0	!flag for upper layers
        end do

        end

c*****************************************************************

        logical function layer_exist(l,k)

c returns a boolean to skip non existing layers in dry areas	
c in case of z-layer modelling:
c                               
c                                  h  
c        
c       _______         layer 1 =  >0
c       .......x------
c              !        layer 2 =  =0
c              x------
c              !______  layer 3 =  >0
c              x......                          
c       
c

        !use mod_geom_dynamic
        use mod_layer_thickness 
        use zadapt

	implicit none

	integer,  intent(in) :: l
	integer, intent(in) :: k
	!logical :: iskout

        !iskout(k)   = inodv(k) .lt. 0

	layer_exist = .true.
	if (iskbathywall(l,k)==1 .and. hdknv(l,k).eq.0) then
	  layer_exist = .false.
	end if

	end function

c*****************************************************************  

	subroutine modify_ccoeff(l,k,nlvddi,aa,ad,cn,cdiag,clow,chigh)

c returns modified coefficient for the tracer system in non existing 
c layers in dry areas, in case of z-layer modelling. These coefficient
c correspond to update the tracerv value in non-existing layer copying 
c the value from the upwind side.

        use mod_hydro_vel

	implicit none

        integer, intent(in) :: l
        integer, intent(in) :: k
        integer, intent(in) :: nlvddi
        double precision, intent(in) :: aa
        double precision, intent(in) :: ad
        double precision,dimension(nlvddi),intent(inout) :: cn
        double precision,dimension(nlvddi),intent(inout) :: cdiag
        double precision,dimension(nlvddi),intent(inout) :: clow
        double precision,dimension(nlvddi),intent(inout) :: chigh

	integer :: ll
	double precision :: wbool
        logical :: layer_exist

c Searching the upwind layer

	if (wlnv(l,k) .gt. 0) then
	  wbool = 1.
	  do ll=l+1,nlvddi
            if (layer_exist(ll,k)) then
              exit
            end if
	  end do
	else
	  wbool = 0.
          do ll=l-1,1,-1 
            if (layer_exist(ll,k)) then
              exit
            end if
          end do
	end if

c New coeff: c^{n+1}_{l} = c^{n+1}_{upwind}

	if (aa .eq. 0. .and. ad .eq. 0.) then 
          cn(l) = cn(ll)
          cdiag(l) = cdiag(ll)
	else
	  cn(l) = 0.
          chigh(l) = -cdiag(ll) * wbool
          cdiag(l) =  cdiag(ll) 
          clow(l)  = -cdiag(ll) * (1.-wbool)
	end if

	end 

c*****************************************************************	

        subroutine copy_zadaptation

c copies top layer index to old time step
c copies areas to old time step
c initialize remap index

        use levels
	use mod_area

        implicit none

        jlhov  = jlhv
        jlhkov = jlhkv
	jlheov = jlhev
	rlhv   = 0 	!initialize 0: no config change
	rlhkv  = 0 	!initialize 0: no config change
	areakov = areakv

        end

c*****************************************************************

        subroutine remap_hydro

c conservative remap of the transports when a top
c layer is removed/inserted	

	use levels	
        use mod_layer_thickness
        use mod_hydro
	use basin, only : nel

        implicit none

	integer ie,l,lminnv,lminov
	real htot  			 !total depth of inserted layers
	real utswap,vtswap		 !swap variable for transports

        do ie=1,nel

	  lminnv = jlhv(ie)
          lminov = jlhov(ie)

	  if (lminnv > lminov) then 	 !top layer element removed

	    do l=lminnv-1,lminov,-1	 !loop over removed layers 
	      utlnv(lminnv,ie) = utlnv(lminnv,ie) + utlnv(l,ie)
              vtlnv(lminnv,ie) = vtlnv(lminnv,ie) + vtlnv(l,ie)	      
      	      utlnv(l,ie) = 0.0		 !cleaning removed layer	   
	      vtlnv(l,ie) = 0.0          !
	    end do

	  else if (lminnv < lminov) then !top layer inserted

	    htot = 0.0
	    do l=lminov,lminnv,-1 
              htot = htot + hdenv(l,ie)   
            end do
	    utswap = utlnv(lminov,ie)
	    vtswap = vtlnv(lminov,ie)
	    do l=lminov,lminnv,-1 	 !loop over inserted layers 
	      utlnv(l,ie) = utswap * hdenv(l,ie)/htot
	      vtlnv(l,ie) = vtswap * hdenv(l,ie)/htot
            end do	      

c	  else				 !no insertion/removal happened 	
	  endif    

	end do

        end

c*****************************************************************

	subroutine remap_new_depth

c shell (helper) for remapdepth

	use mod_layer_thickness
	use mod_area
	use mod_hydro
	use levels, only : nlvdi,nlv

	implicit none

	call remapdepth(nlvdi,hdknv,hdenv,zenv,areakv)

	end

c*****************************************************************

	subroutine remapdepth(levdim,hdkn,hden,zenv,area)

c remap depth array for nodes

	use zadapt
	use mod_depth
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	integer levdim
	real hdkn(levdim,nkn)	  !depth at node, new time level
	real hden(levdim,nel)	  !depth at element, new time level
	real zenv(3,nel)    	  !water level at new time level
        real area(levdim,nkn)

        logical bdebug
	integer k,l,ie,ii
	integer lmax,lmin,lremap,lmink,lremapk,n
	real hfirst,hlast,h,htot,z,zmed,hm
	real hmin

	real areael,areafv
	real areaele

        logical :: islay,layer_exist

        bdebug = .false.
	hmin = -99999.
	hmin = 0.

c----------------------------------------------------------------
c initialize and copy
c----------------------------------------------------------------

	!cleaning layer to be remapped
	hdkn = hdkn * real(iskremap)
        hden = hden * real(iseremap)

c----------------------------------------------------------------
c compute volumes at node
c----------------------------------------------------------------

	do ie=1,nel

	  !call elebase(ie,n,ibase)
	  n = 3
	  areael = 12 * ev(10,ie)
	  areafv = areael / n

	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  lremap = rlhv(ie)
	  hm = hev(ie)
	  zmed = 0.

c	  -------------------------------------------------------
c	  nodal values
c	  -------------------------------------------------------

	  do ii=1,n

	    k = nen3v(ii,ie)
	    lmink = jlhev(ii,ie)
            lremapk = min(rlhkv(k),lmax)
	    htot = hm3v(ii,ie)
	    z = zenv(ii,ie)
	    zmed = zmed + z

	    do l=lmink,lremapk
	      hdkn(l,k) = hdkn(l,k) + areafv * hldv(l)
	      if ( l .eq. lmink) then
	        hdkn(lmink,k) = hdkn(lmink,k) + areafv *(hlv(lmink-1) + z)
	      end if
	    end do
            if( lremapk .eq. lmax ) then
              hlast = htot - hlv(lmax-1)
              !if( hlast .lt. 0. ) goto 77  !lrp: remove this line	      
	      hdkn(lmax,k) = hdkn(lmax,k) + areafv * (hlast-hldv(lremapk))
            end if

	    !do l=1,lmax
	    !  if( hdkn(l,k) <= 0. ) then
	    !    write(6,*) 'no volume in node ',k
	    !    write(6,*) 'depth: ',h
	    !    write(6,*) 'nlv,lmax: ',nlv,lmax
	    !    write(6,*) 'nsigma: ',nsigma
	    !    write(6,*) 'hlv: ',hlv
	    !    write(6,*) 'hldv: ',hldv
	    !    write(6,*) 'hdkn: ',hdkn(:,k)
	    !    stop 'error stop remapdepth: no volume in node'
	    !  end if
	    !end do

	  end do

!	  in hdkn is volume of finite volume around k
!	  in sigma layers hdkn already contains nodal depth

c	  -------------------------------------------------------
c	  element values
c	  -------------------------------------------------------

	  k = 0
	  ii = 0
	  zmed = zmed / n
	  htot = hev(ie)

	  do l=lmin,lremap
	    hden(l,ie) = hldv(l)
	    if ( l .eq. lmin) then
	      hden(lmin,ie) = hden(lmin,ie) + hlv(lmin-1) + zmed    
	    end if
	  end do
	  if( lremap .eq. lmax ) then
  	    hlast = htot - hlv(lmax-1)
	    if( hlast .lt. 0. ) goto 77
	    hden(lmax,ie) = hden(lmax,ie) + (hlast-hldv(lremap))
	  end if

          do l=lmin,lremap
            h = hden(l,ie)
            if( h <= hmin ) then
              write(6,*) 'error computing layer thickness'
              write(6,*) 'no layer depth in element: ',ie,l,lmax
              write(6,*) 'depth: ',h,htot,zmed
              write(6,*) 'additional information available in fort.666'
              call check_set_unit(666)
              call check_elem(ie)
              call check_nodes_in_elem(ie)
              stop 'error stop remapdepth: no layer in element'
            end if
          end do

	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange hdkn')
          call shympi_exchange_and_sum_3d_nodes(hdkn)
	end if

c----------------------------------------------------------------
c compute depth at nodes
c----------------------------------------------------------------

	do k=1,nkn_inner
	  lmin = jlhkv(k)
          lremapk = rlhkv(k)
	  do l=lmin,lremapk
	    areafv = area(l,k)
	    if( areafv .gt. 0. ) then
	      hdkn(l,k) = hdkn(l,k) / areafv
	    end if
	  end do
          do l=lmin,lremapk
            h = hdkn(l,k)
            islay = layer_exist(l,k)
            if( (h <= hmin) .and. islay ) then
              write(6,*) 'error computing layer thickness'
              write(6,*) 'no layer depth in node: ',k,l,lmax
              write(6,*) 'depth: ',h
              write(6,*) 'nlv,lmax: ',nlv,lmax
              write(6,*) 'hlv: ',hlv
              write(6,*) 'hldv: ',hldv
              write(6,*) 'additional information available in fort.666'
              call check_set_unit(666)
              call check_node(k)
              call check_elems_around_node(k)
              stop 'error stop remapdepth: no layer in node'
            end if
          end do
	end do

c----------------------------------------------------------------
c echange nodal values
c----------------------------------------------------------------

	if( shympi_partition_on_nodes() ) then
	  !call shympi_comment('exchanging hdkn')
	  call shympi_exchange_3d_node(hdkn)
	  !call shympi_barrier
	end if

c----------------------------------------------------------------
c end of routine
c----------------------------------------------------------------

	return
   77	continue
	write(6,*) 'last layer negative: ',hlast
	write(6,*) 'levdim,lmax: ',levdim,lmax
	write(6,*) 'ie,ii,k: ',ie,ii,k
	write(6,*) 'htot,hm: ',htot,hm
	write(6,*) 'hm3v  ',(hm3v(ii,ie),ii=1,3)
	write(6,*) 'hlv   ',(hlv(l),l=1,levdim)
	write(6,*) 'hldv  ',(hldv(l),l=1,levdim)
	stop 'error stop remapdepth: last layer negative'
	end

c*****************************************************************	

	subroutine remap_hydro_vertical

c remap (here it means only recompute) vertical velocities
c
c wlnv (dvol)   aux array for volume difference
c vv            aux array for area

	use mod_bound_geom
	use mod_geom_dynamic
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
        use mod_layer_thickness
	use evgeom
	use levels
	use basin
	use shympi
	use zadapt

	implicit none

	logical debug
	integer k,ie,ii,kk,l,lmin,lmiss,lmink,lremapk
	integer jlevel,ilevel,jlevelmin
        integer ibc,ibtyp
	real aj,wbot,wdiv,ff,atop,abot,wfold
	real b,c
	real am,az,azt,azpar,ampar
	real ffn,ffo
	real volo,voln,dt,dvdt
	real, allocatable :: vf(:,:)
	real, allocatable :: va(:,:)
        double precision,dimension(nlvdi) :: wein,weio
c statement functions

	logical is_zeta_bound
	real volnode

	!logical isein
        !isein(ie) = iwegv(ie).eq.0

c 2d -> nothing to be done

	wlnv(0:(nlvdi-1),:) = wlnv(0:(nlvdi-1),:) * real(iskremap)
	if( nlvdi == 1 ) return

c initialize

	call getazam(azpar,ampar)
	az=azpar
	am=ampar
	azt = 1. - az
	call get_timestep(dt)

	allocate(vf(nlvdi,nkn),va(nlvdi,nkn))
	vf = 0.
	va = 0.

c compute difference of velocities for each layer
c
c f(ii) > 0 ==> flux into node ii
c aj * ff -> [m**3/s]     ( ff -> [m/s]   aj -> [m**2]    b,c -> [1/m] )

	do ie=1,nel

	 !if( isein(ie) ) then		!FIXME	
	  aj=4.*ev(10,ie)		!area of triangle / 3
	  ilevel = ilhv(ie)
	  weio = 0.; wein = 0.
	  call get_zadaptive_weights(ie,jlevelmin,wein,+1)
          call get_zadaptive_weights(ie,jlevelmin,weio,-1)

	  do ii=1,3

  	    kk=nen3v(ii,ie)
            b = ev(ii+3,ie)
            c = ev(ii+6,ie)	  
	    lremapk = min(rlhkv(kk),ilevel)

	    !configuration may be changed between
	    !time^{n+1} and time^{n} due to element insertion/removal:
	    !we treat separately time level t^{n+1} and t^{n}	    
	    lmink = jlhev(ii,ie)	! t^{n+1}
	    jlevel = jlhv(ie)
	    do l=lmink,lremapk
              !element with non conformal edge
	      if (l.le.jlevel .and. jlevel.gt.lmink) then
		ffn = (utlnv(jlevel,ie)*b + vtlnv(jlevel,ie)*c)*wein(l)
	      !general case for inferior layers
	      !or conformal edge
	      else
	        ffn = utlnv(l,ie)*b + vtlnv(l,ie)*c
	      end if
              ff = ffn * az
              vf(l,kk) = vf(l,kk) + 3. * aj * ff
              va(l,kk) = va(l,kk) + aj
	    end do

	    lmink = jlheov(ii,ie)       ! t^{n}  
	    jlevel = jlhov(ie)  
            do l=lmink,lremapk
              !element with non conformal edge
              if (l.le.jlevel .and. jlevel.gt.lmink) then
                ffo = (utlov(jlevel,ie)*b + vtlov(jlevel,ie)*c)*weio(l)
              !general case for inferior layers
              !or conformal edge
              else
                ffo = utlov(l,ie)*b + vtlov(l,ie)*c
              end if
              ff = ffo * azt
              vf(l,kk) = vf(l,kk) + 3. * aj * ff
              !va(l,kk) = va(l,kk) + aj
            end do
	      
	  end do
	 !end if		
	end do

	if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange vf,va')
          call shympi_exchange_and_sum_3d_nodes(vf)
          call shympi_exchange_and_sum_3d_nodes(va)
	end if

c from vel difference get absolute velocity (w_bottom = 0)
c	-> wlnv(nlv,k) is already in place !
c	-> wlnv(nlv,k) = 0 + wlnv(nlv,k)
c w of bottom of last layer must be 0 ! -> shift everything up
c wlnv(nlv,k) is always 0
c
c dividing wlnv [m**3/s] by area [vv] gives vertical velocity
c
c in va(l,k) is the area of the upper interface: a(l) = a_i(l-1)
c =>  w(l-1) = flux(l-1) / a_i(l-1)  =>  w(l-1) = flux(l-1) / a(l)

	do k=1,nkn
	  lremapk = rlhkv(k)
	  lmin = jlhkv(k)
	  debug = k .eq. 0
	  abot = 0.
	  do l=lremapk,lmin,-1
	    atop = va(l,k)
            voln = volnode(l,k,+1)
            volo = volnode(l,k,-1)
	    dvdt = (voln-volo)/dt
	    wdiv = vf(l,k)
	    !configuration may be changed:
	    !removed layers are remapped on the new grid
	    if (l.eq.lmin .and. lmin.gt.jlhkov(k)) then	    
	      do lmiss=lmin-1,jlhkov(k),-1
	        dvdt = dvdt - volnode(lmiss,k,-1)/dt
	        wdiv = wdiv + vf(lmiss,k)
	      end do
	    end if	    
	    !wfold = azt * (atop*wlov(l-1,k)-abot*wlov(l,k))
	    !wlnv(l-1,k) = wlnv(l,k) + (wdiv-dvdt+wfold)/az
	    wlnv(l-1,k) = wlnv(l,k) 
	    if( atop .gt. 0. ) then  
	      wlnv(l-1,k) = wlnv(l-1,k) + (wdiv - dvdt)/atop
	    end if  
	    abot = atop
	    if( debug ) write(6,*) k,l,wdiv,wlnv(l,k),wlnv(l-1,k)
	    if (l.eq.lmin) then
              wlnv(lmin-1,k) = 0.   ! ensure no flux across surface - is very small
	    end if  
	  end do
	end do

c set w to zero at open boundary nodes (new 14.08.1998)
c
c FIXME	-> only for ibtyp = 1,2 !!!!

	do k=1,nkn
          !if( is_external_boundary(k) ) then	!bug fix 10.03.2010
          if( is_zeta_bound(k) ) then
            lremapk = rlhkv(k)		  
	    wlnv(0:lremapk,k) = 0.
          end if
	end do

	deallocate(vf,va)

	if( shympi_partition_on_nodes() ) then
	  !call shympi_comment('exchanging wlnv')
          call shympi_exchange_3d0_node(wlnv)
	  !call shympi_barrier
	end if

	end

c*****************************************************************	  

	subroutine remap_ttov

c transforms transports to velocities

	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use zadapt

	implicit none

	if( .not. mod_layer_thickness_is_initialized() ) then
	  write(6,*) 'layer thickness is not initialized'
	  stop 'error stop ttov: no layerthickness'
	end if

	ulnv = ulnv * real(iseremap)
	vlnv = vlnv * real(iseremap)
	where( hdenv > 0. .and. iseremap .eq. 0)
	  ulnv = utlnv / hdenv
	  vlnv = vtlnv / hdenv
	end where

	end

c*****************************************************************

	subroutine remap_uvtopr

c transforms velocities to nodal values

	use mod_geom_dynamic
	use mod_hydro_print
	use mod_hydro_vel
	use evgeom
	use levels
	use basin
	use shympi
	use zadapt

	implicit none

	integer ie,l,lmiss,k,ii
	integer lmax,lmin,lremapk,lmink
	real aj
	!real vv(nlvdi,nkn)
	real, allocatable :: vv(:,:)

	allocate(vv(nlvdi,nkn))
	uprv = uprv * real(iskremap)
	vprv = vprv * real(iskremap)
	vv   = 0.

c baroclinic part

	do ie=1,nel
	  if ( iwegv(ie) /= 0 ) cycle
          lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  aj=ev(10,ie)
          do ii=1,3
	    k=nen3v(ii,ie)
	    lremapk = min(rlhkv(k),lmax)
	    lmink = jlhev(ii,ie)
	    do l=lmink,lremapk
	      vv(l,k)=vv(l,k)+aj
	      !element with non-conformal edge
	      if (l.lt.lmin) then 
	        uprv(l,k)=uprv(l,k)+aj*ulnv(lmin,ie)
		vprv(l,k)=vprv(l,k)+aj*vlnv(lmin,ie)		
	      !standard element	
              else
                uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
                vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)      		      
	      end if
	    end do
	  end do
	end do

        !call shympi_comment('shympi_elem: exchange uprv, vprv')
        call shympi_exchange_and_sum_3d_nodes(uprv)
        call shympi_exchange_and_sum_3d_nodes(vprv)
        call shympi_exchange_and_sum_3d_nodes(vv)

	where ( vv > 0. ) 	!this automatically select only node 
	  uprv = uprv / vv      !to be remapped
	  vprv = vprv / vv
	end where

	call shympi_exchange_3d_node(uprv)
	call shympi_exchange_3d_node(vprv)

c vertical velocities -> we compute average over one layer

	do l=1,nlv		!FIXME not-optmized
	  wprv(l,:)=0.5*(wlnv(l,:)+wlnv(l-1,:))
	end do

	deallocate(vv)

	end

c*****************************************************************

	subroutine remap_velocities

c computes horizontal velocities from zenv, utlnv, vtlnv, hdenv

	implicit none

	call remap_ttov			!velocities ulnv/vlnv
	call remap_uvtopr		!nodal values uprv/vprv/wprv

	end

c*****************************************************************	

        subroutine scal_remap_removal(nlvdi,cn,hn,lminknv,lminkov)

c conservative remap of scalar vars in case of layer removal
c \tilde{hc}_{l}=hc_{l}+hc_{l-1} 
c -> \tilde{c}_{l}=(hc_{l}+hc_{l-1})/\tilde{h}

	implicit none

        integer,intent(in) :: nlvdi		
        real,dimension(nlvdi),intent(inout) :: cn
        real,dimension(nlvdi),intent(in) :: hn
        integer,intent(in) :: lminknv
	integer,intent(in) :: lminkov
        integer l
	real ctot,htot

	htot = 0.
	ctot = 0.
        do l=lminknv,lminkov,-1      	!loop over removed layers
	  ctot = ctot + cn(l) * hn(l)
	  htot = htot + hn(l)
	  cn(l) = 0.0
	end do
	cn(lminknv) = ctot/htot	

	end	

c*****************************************************************

        subroutine scal_remap_insertion(nlvdi,cn,lminknv,lminkov)

c conservative remap of scalar vars in case of layer insertion
c \tilde{hc}_{l-1}=hc_{l}*\tilde{h}/h -> \tilde{c}_{l-1} = c_{l}

	implicit none

        integer,intent(in) :: nlvdi		
        real,dimension(nlvdi),intent(inout) :: cn
        integer,intent(in) :: lminknv
        integer,intent(in) :: lminkov
        integer l	

        do l=lminkov-1,lminknv,-1	!loop over inserted layers
          cn(l) = cn(lminkov)
	end do

        end


c*****************************************************************

        subroutine remap_scalar

c conservative remap of scalar vars for
c configuration changes due top layer insertion/removal.
c On the new configuation we compute tracers on nodes

        use levels
        use mod_layer_thickness
        use mod_ts
        use basin !, only : nel,nkn,nen3v,hm3v
        use shympi

        implicit none

	integer k,lminknv,lminkov
	integer itemp,isalt
	real getpar

        itemp=nint(getpar('itemp'))
        isalt=nint(getpar('isalt'))

c----------------------------------------------------------------
c loop over nodes
c----------------------------------------------------------------

        do k=1,nkn

          lminknv = jlhkv(k)
          lminkov = jlhkov(k)

c----------------------------------------------------------------
c top layer node removed
c----------------------------------------------------------------

          if (lminknv > lminkov) then

	    !remap temp
            if( itemp .gt. 0 ) then
	      call scal_remap_removal(nlvdi,tempv(:,k),
     +				      hdknv(:,k),lminknv,lminkov)
      	    end if
            if( isalt .gt. 0 ) then
              call scal_remap_removal(nlvdi,saltv(:,k),
     +                                hdknv(:,k),lminknv,lminkov)
            end if

c----------------------------------------------------------------
c top layer node inserted
c----------------------------------------------------------------

	  else if (lminknv < lminkov) then

            !remap temp
            if( itemp .gt. 0 ) then
	      call scal_remap_insertion(nlvdi,tempv(:,k),
     +					lminknv,lminkov)
	    end if
            if( isalt .gt. 0 ) then
              call scal_remap_insertion(nlvdi,saltv(:,k),
     +                                  lminknv,lminkov)
            end if

c	  else
	  end if

	end do

	end

c******************************************************************

        subroutine init_zadaptation

        use levels, only : nlvdi,hlv
        use basin, only : nkn,nel		
	use zadapt

	implicit none		

	integer lmin,nzadapt
	real getpar,rzmov,rztop,testz

        allocate(nadapt_com(4,nel))
        allocate(hadapt_com(4,nel))	
	allocate(iseremap(nlvdi,nel))
	allocate(iskremap(nlvdi,nkn))
	allocate(iskbathywall(nlvdi,nkn))

        nadapt_com  = 0
	hadapt_com  = 0.	
	iseremap = 1
	iskremap = 1
	call check_bathywall

	nzadapt = nint(getpar('nzadapt'))
	call set_nzadapt_info(nzadapt)
        call compute_rzpar_info(nlvdi,nzadapt,hlv,rzmov,rztop)
        call set_rzmov_info(rzmov)
	call set_rztop_info(rztop)

        write(6,'(a)') ' Initializing z-layers parameters ...'
        write(6,*) ' nzadapt,rzmov,rztop: ', nzadapt,rzmov,rztop

        if (nzadapt .le. 1) then                !z-layers
          write(6,*) ' z-layers'
        else if (nzadapt .ge. nlvdi) then         !z-star
          write(6,*) ' z-star layers'
        else                                    !z + z-star
          write(6,*) ' z-star + z-layers:'
c          write(6,*) ' z (water level), (nzadapt) number of'
c          write(6,*) ' surface layers moving with z-star:'
c          lmin = 1
c          testz = -0.0
c          call compute_nadapt_info(testz,hlv,nlv,lmin,nzadapt)
c          write(6,*) ' z nzadapt: ', testz,nzadapt
c
c          testz = -0.5
c          call compute_nadapt_info(testz,hlv,nlv,lmin,nzadapt)
c          write(6,*) ' z nzadapt: ', testz,nzadapt
c
c          testz = -1.0
c          call compute_nadapt_info(testz,hlv,nlv,lmin,nzadapt)
c          write(6,*) ' z nzadapt: ', testz,nzadapt
c
c          testz = -1.5
c          call compute_nadapt_info(testz,hlv,nlv,lmin,nzadapt)
c          write(6,*) ' z nzadapt: ', testz,nzadapt
        end if


	end

c*****************************************************************

        subroutine zadaptation

c adaptation routines for z-layer coordinate driven free-surface 
c movement:
c 1/ change the vertical grid set_jlhv, set_jlhkv,set_area
c 2/ remap variable           remap_new_depth,remap_hydro,remap_scalar
c			      remap_velocities,remap_hydro_vertical		
c important: set_area must be always called immediatly after set_jlhkv

	use basin

	implicit none

c----------------------------------------------------------------
c change vertical grid
c----------------------------------------------------------------

	call set_jlhv                !set top layer index (elemental)
        call set_jlhkv               !set top layer index (nodal)
	call set_area		     !re-compute areas with new index

c----------------------------------------------------------------
c remap all variables
c----------------------------------------------------------------
	
	call remap_scalar            !remap ts vars (need old depth)	
	call remap_new_depth	     !remap depth vars: hdenv,hdknv
	call remap_hydro             !reamp hydro vars (need new depth)
	call remap_hydro_vertical    !remap vertical velocity 
        call mass_conserve           !check mass balance
	call remap_velocities        !re-compute velocities (elements and nodes)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************
c	program sigma_main
c	call sigma_test
c	end
c******************************************************************

