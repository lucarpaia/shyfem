
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
c	compute_zadaptive_info (subele.f)
c	get_zadaptivelayer_thickness (subele.f)
c	zadaptation (shyfem.f)
c	get_zadaptive_weights (new3di.f)
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module zadapt
!==================================================================

	implicit none

	real, save :: rmin_gridmov = 0.125 !threshold grid movement
	real, save :: rmin_gridtop = 0.200 !threshold grid topology

!==================================================================
        end module zadapt
!==================================================================

c******************************************************************

        subroutine compute_zadaptive_info(ie,nlv,lmin,lmax,hlv,z,
     +					  nsigma,lsigma,hsigma,hdl)

c returns adaptive layers info. Info is computed by node of element:
c number of adaptive layers by node (3) + by ele (1)
c lowest index of adaptive layer by node (3) + by ele (1)
c closing depth of adaptive layer
c coefficients of adaptive layers

	use zadapt

        implicit none

        integer ie              !element index
        integer nlv             !total number of layers
        integer lmin(3)         !top layer index	
	integer lmax		!bottom layer index
	real hlv(nlv)           !layer structure
        real z(3)               !water level	
        integer nsigma(4)       !total number of adaptive layers (return)
	integer lsigma(4)       !lowest index of adaptive layers (return)
        real hsigma(3)          !closing depth of adaptive layers(return)
	real hdl(nlv,3)         !coefficient (return)

	integer l,ii,levmax,lmine
	real r,q,den

	nsigma = 0         !no adaptation -> all to zero
	lsigma = 0	   !no adaptation -> all to zero
	hdl = 0.           !no adaptation -> all to zero

	r = rmin_gridmov
	q = 1./r

	lmine = maxval(lmin)

c---------------------------------------------------------
c loop over nodes: adaptation is node-driven
c---------------------------------------------------------

	do ii=1,3
	  levmax = 0       !no adapation -> all to zero
	  do l=lmin(ii),lmax-1 !-1 to skip bottom layer

c---------------------------------------------------------
c node movement for layer needed to avoid negative depth
c---------------------------------------------------------	  
	  
            if(z(ii).le.(-hlv(l)+r*(hlv(l)-hlv(l-1)))) then 
              levmax = l+1

c---------------------------------------------------------
c no movement needed 
c--------------------------------------------------------- 

	    else  
	      exit 
	    end if  
          end do

c---------------------------------------------------------
c compute nsigma, lsigma, hsigma 
c---------------------------------------------------------  	  

	  hsigma(ii) = hlv(levmax) 
	  nsigma(ii) = max(0,levmax-lmin(ii)+1) !+1 (min 2 adaptive layer)
	  lsigma(ii) = levmax

c---------------------------------------------------------
c compute hdl: we freeze the configuration
c--------------------------------------------------------- 

    	  if (nsigma(ii).gt.0) then	  
	  den = (nsigma(ii)-1.)*q+1.
	  hdl(lmin(ii),ii) = - 1. / den
c	  hdl(lmin(ii),ii) = - 1. / nsigma(ii)	 !constant hdl
          do l=lmin(ii)+1,lsigma(ii)
c	    hdl(l,ii) = - 1. / nsigma(ii)	 !constant hdl
	    hdl(l,ii) = - q / den
	  end do
	  end if

c---------------------------------------------------------
c compute element info: number of adaptive layers in ele
c---------------------------------------------------------

	  !flag adaptive layers in element with non-conformal edge:
	  !free-surface must span all layers greater then lmin
	  if (lsigma(ii).gt.lmine) then
            nsigma(4) = max(nsigma(ii),nsigma(4))
	    lsigma(4) = max(lsigma(ii),lsigma(4))
	  end if
	end do

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

c******************************************************************

        subroutine get_zadaptive_weights(ie,weightn,weighto)

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
	double precision,dimension(nlvdi) :: weightn !weights at new time-level(return)
        double precision,dimension(nlvdi) :: weighto !weights at old time-level(return)

	integer jlevel,jlevelmin,jlevelk
	integer lmiss,k,ii,n,numOfLev
	double precision wn,wo
	logical isenonconf

	n = 3
	jlevel = jlhv(ie)
        jlevelmin = huge(1)
        do ii=1,3
          jlevelmin=min(jlevelmin,jlhev(ii,ie))
        end do

c---------------------------------------------------------
c is element non-conformal?
c---------------------------------------------------------	

        isenonconf = .false.
        if (jlevel.gt.jlevelmin) isenonconf = .true.	

	if (isenonconf) then

        wn = 0.
        wo = 0.	

c---------------------------------------------------------
c loop over sub-layers
c---------------------------------------------------------

        do lmiss=jlevel,jlevelmin,-1

          weightn(lmiss) = 0.
	  weighto(lmiss) = 0.

c---------------------------------------------------------
c loop over node -> compute element sub-layer depth
c---------------------------------------------------------

	  do ii=1,n

	    k=nen3v(ii,ie)
            jlevelk = jlhev(ii,ie)

	    if (lmiss.gt.jlevelk) then
              weightn(lmiss) = weightn(lmiss) + hdknv(lmiss,k)
              weighto(lmiss) = weighto(lmiss) + hdkov(lmiss,k) 
            else
              numOfLev = real(jlevelk-jlevelmin+1)
	      weightn(lmiss) = weightn(lmiss) + hdknv(jlevelk,k)/numOfLev
	      weighto(lmiss) = weighto(lmiss) + hdkov(jlevelk,k)/numOfLev
            end if

	  end do

          wn = wn + weightn(lmiss)
          wo = wo + weighto(lmiss)

	end do

c---------------------------------------------------------
c compute sub-layer weight
c---------------------------------------------------------

        do lmiss=jlevel,jlevelmin,-1

	  !if (wn.gt.0) then 
	    weightn(lmiss) = weightn(lmiss)/wn
	  !end if
          if (wo.gt.0) then		  
            weighto(lmiss) = weighto(lmiss)/wo
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
c jlhv(ie)	layer index of the upper existing layer updated 
c jlhvo(ie)	layer index of the upper existing layer before update	
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

	!real getpar,hzoff
        real zmin,r

	lmax=0
        !hzoff=getpar('hzoff')
	r = rmin_gridtop

	call get_sigma_info(nlev,nsigma,hsigma)
	bsigma = nsigma .gt. 0

	jlhov=jlhv    		! swap for saving old index
				! first timestep: jlhov=1 even if some
				! layer are removed but it is unused

c---------------------------------------------------------
c loop over element
c---------------------------------------------------------

	do ie=1,nel

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
		
c set nodal top layer index jlhkv,jlhev and jlhkov,jlheov array - 
c only needs zenv
c a double index is needed (by node and elem) for wet/dry interfaces:
c top layer depends on the free surface and at wet/dry interfaces 
c the free-surface is discontinous -> element storing is also required
c jlhkv(k)     layer index of the upper existing node updated (node)
c jlhev(ii,ie) layer index of the upper existing node updated (elem)
c jlhkov(k)    layer index of the upper existing node before update   
c jlhkev(ii,ie)layer index of the upper existing node before update
c 
c e.g. jlhkv(k)=jlhkv(ii,ie)=2 means first layer has been removed 
c      for node k <-> (ii,ie) but only if k is not wet/dry

        use zadapt
        use levels
        use basin
        use shympi
	use mod_hydro
	use mod_layer_thickness

        implicit none

	logical bsigma
        integer ie,ii,k,l,nlev,nsigma
        real hsigma
	!logical isein,iseout,iskin
	!real getpar,hzoff 	  !lrp
	real r 

	r = rmin_gridtop
        !hzoff=getpar('hzoff')

        call get_sigma_info(nlev,nsigma,hsigma)
        bsigma = nsigma .gt. 0

        jlhkov=jlhkv              !swap for saving old index
	jlheov=jlhev			

        do k=1,nkn
          jlhkv(k)=huge(1)
        end do

c---------------------------------------------------------
c loop over element
c---------------------------------------------------------

        do ie=1,nel
          do ii=1,3
	    k=nen3v(ii,ie)
            if( bsigma ) then     !sigma levels always start from 1
              l = 1
            else
              do l=1,nlv
                if((-hlv(l)+r*(hlv(l)-hlv(l-1))).lt.zenv(ii,ie)) exit
                !if((-hlv(l)+hzoff).lt.zenv(ii,ie)) exit		
              end do
            end if
	    !safety min: jlhev>=jlhv, jlhkv>=ilhkv
	    jlhev(ii,ie)=min(l,jlhv(ie))
	    if (jlhev(ii,ie).lt.jlhkv(k)) then 
	      jlhkv(k)=min(jlhev(ii,ie),ilhkv(k))
	    end if
          end do
        end do

        if( shympi_partition_on_elements() ) then
          !call shympi_comment('shympi_elem: exchange ilhkv - max')
          call shympi_exchange_2d_nodes_max(jlhkv)
        else
          call shympi_exchange_2d_node(jlhkv)
        end if

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

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
	    do l=lminov,lminnv,-1 	 !loop over inserted layers 
	      utlnv(l,ie) = utlnv(lminov,ie) * hdenv(l,ie)/htot
	      vtlnv(l,ie) = vtlnv(lminov,ie) * hdenv(l,ie)/htot 
            end do	      

c	  else				 !no insertion/removal happened 	
	  endif    

	end do

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

c*****************************************************************

        subroutine zadaptation

c adaptation routines for z-layer coordinate driven free-surface 
c movement:
c 1/ change the vertical grid set_jlhv, set_jlhkv,set_area
c 2/ remap variable           make_new_depth,remap_hydro,remap_scalar
c important: set_area must be always called immediatly after set_jlhkv

	use basin

	implicit none

        real dzeta(nkn)

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
	call make_new_depth	     !remap depth vars: hdenv,hdknv
	call remap_hydro             !reamp hydro vars (need new depth)
	call hydro_vertical(dzeta)   !remap vertical velocity 
	call compute_velocities      !re-compute velocities (elements and nodes)

c---------------------------------------------------------
c end of routine
c---------------------------------------------------------

	end

c******************************************************************
c	program sigma_main
c	call sigma_test
c	end
c******************************************************************

