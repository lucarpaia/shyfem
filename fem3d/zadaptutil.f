
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

c routines for handling z-adaptive layers in output programs
c
c revision log :
c
c 28.03.2022    lrp     introduce z-adaptive layers
c
c notes:
c this file is used also in:      
c      
c	compute_zadapt_info (nosutil.f,ousutil.f,shyutil.f,shyelab_util.f,
c			     shyelab_nodes.f,shyelab_average.f)
c
c******************************************************************
c******************************************************************
c******************************************************************

!==================================================================
        module zadaptutil
!==================================================================

        implicit none

        real, save :: rzmov_com                    !parameter for moving surface layers
	real, save :: rztop_com		      !parameter for removing surface layers

!==================================================================
        end module zadaptutil
!==================================================================

c******************************************************************

        subroutine get_rzpar_info(rzmov,rztop)

        use zadaptutil

        implicit none

        real rzmov,rztop

        rzmov = rzmov_com
	rztop = rztop_com

        end

c******************************************************************

        subroutine set_rzpar_info(rzmov,rztop)

        use zadaptutil

        implicit none

        real rzmov,rztop

        rzmov_com = rzmov
	rztop_com = rztop

        end

c******************************************************************

        subroutine compute_rzpar_info2(nlv,nzadapt,hlv,rzmov,rztop)

        use zadaptutil

        implicit none

        integer nlv             !total number of layers
        integer nzadapt         !number of surface moving layers
        real hlv(nlv)           !layer structure
        real rzmov              !parameter for moving surface layers (return)
	real rztop		!parameter for removing surface layers (return)

        real maxz

        maxz = 0.0                              !estimate of max water level

        if (nzadapt .le. 0) then                !z-layers
          rzmov = 0.125
	  rztop = 0.200
        else if (nzadapt .ge. nlv) then         !z-star
          rzmov =  10000.
	  rztop = -10000.
        else                                    !z + z-star
          rzmov = (maxz+hlv(nzadapt-1)) / (hlv(nzadapt)-hlv(nzadapt-1))
	  rztop = -10000.
        end if

        end

c******************************************************************

        subroutine init_rzmov_info(nlv,nzadapt,hlv)

        implicit none

        integer nlv
        integer nzadapt
        real hlv(nlv)

        real rzmov,rztop


        call compute_rzpar_info2(nlv,nzadapt,hlv,rzmov,rztop)
        call set_rzpar_info(rzmov,rztop)

        end

c******************************************************************

	subroutine compute_zadapt_info(z,hlv,nsig,lmax,lmin,nzad,hzad)

c returns relevant info for z-adaptive layers

	implicit none

	real z			!water level
	real hlv(lmax)		!layer structure
	integer nsig		!total number of sigma layers
	integer lmax            !bottom layer  index
	integer lmin            !surface layer  index (return)
	integer nzad		!number of z-adaptive layers (return)
	real hzad		!closing depth of z-adaptive layers (return)

	integer l,levmax
	real rgridtop,rgridmov

	call get_rzpar_info(rgridmov,rgridtop)

        lmin = 1
        levmax = 0       !no adapation -> all to zero	

	if( nsig .eq. 0 ) then   

c---------------------------------------------------------
c surface layer index
c---------------------------------------------------------  

      	  do l=1,lmax              !a threshold is used
            if((-hlv(l)+rgridtop*(hlv(l)-hlv(l-1))).le.z) exit
          end do
          lmin=min(l,lmax)         !safety min: jlhv>=ilhv    

c---------------------------------------------------------
c lowest index of adaptive deforming layers
c--------------------------------------------------------- 

          do l=lmin,lmax-1 !-1 to skip bottom layer
            if(z.le.(-hlv(l)+rgridmov*(hlv(l)-hlv(l-1)))) then
              levmax = l+1
            else
              exit
            end if
          end do

	end if 

c---------------------------------------------------------
c compute nzad, hzad
c--------------------------------------------------------- 

        hzad = hlv(levmax)
        nzad = max(0,levmax-lmin+1) !+1 (min 2 adaptive layer)

	end

c******************************************************************
c	program sigma_main
c	call sigma_test
c	end
c******************************************************************

