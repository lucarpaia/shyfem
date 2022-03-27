
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010,2015,2019  Georg Umgiesser
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

! revision log :
!

        real rmin_gridmov
        common /rmin_gridmov/rmin_gridmov

        real rmin_gridtop
        common /rmin_gridtop/rmin_gridtop

        save /rmin_gridmov/,/rmin_gridtop/

        integer iskremap(nlvdim,nkndim)
        common /iskremap/iskremap
	save /iskremap/

        integer iseremap(nlvdim,neldim)
        common /iseremap/iseremap
        save /iseremap/
