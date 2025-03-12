!  Copyright 2013-2018, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
!  Guillermo Ramirez, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
! 
!  This file is part of DAM320.
! 
!  DAM320 is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
! 
!  DAM320 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
! 
!  You should have received a copy of the GNU General Public License
!  along with DAM320.  If not, see <http://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------
!
! Program for computing a density isosurface and its normals from the representation of the molecular density performed with
! DAM320
!
!
! Version of October 2019
!

! #define DBLPRCGRID    ! Uncomment this line  if double precision grid is wanted
!===============================================================================================
!                 MODULE DAMISODEN320_D
!=============================================================================================== 
MODULE DAMISODEN320_D
    USE DAM320_D
    USE DAM320_CONST_D
    IMPLICIT NONE
    logical :: lbinary, lexist
    character(300) :: gridname, outisodenname, outrootname, straux
    integer(KINT) :: i, i1, i2, ia, iaux, ierr, icube, iopt(5), iuni, ix, iy, iz, j, k, knt, kntgrid
    integer(KINT) :: kntind, kntvert, l, npoints, nx, ny, nz
    integer(KINT), allocatable ::  indices(:), interpolmat(:,:), tritable(:,:)
    real(KREAL), parameter :: angstromtobohr = 1.889725989d0
    real(KREAL) :: a, aux, b, bux, c, contourval, cux, d, den, denrep, dendrvx, disthresq, drvxtot, dendrvy
    real(KREAL) :: denmaxerr, drvytot, dendrvz, drvztot, errabs, s, surftot, surftrian, volume, voltetrahed, volvoxel
    real(KREAL) :: x, xini, xinterp, xfin, y, yini, yinterp, yfin, z, zini, zinterp, zfin
    real(KREAL), allocatable :: gradient(:,:), grid(:), fvoxel(:), vertices(:,:)
    real(KREAL), allocatable :: xvoxel(:), yvoxel(:), zvoxel(:)
    real(KREAL) :: xyz(3), xyztetr(3,0:3)
    real(KREAL4), allocatable :: grid4(:)
END MODULE
!
!                 END OF MODULE DAMISODEN320_D
!...............................................................................................

  program DAMISODEN320
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, fA_orig => fA
    USE DAMDEN320_D
    USE DAMISODEN320_D
    implicit none
    integer(KINT) :: jaux
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL) :: geomthr, xmax, xmin, ymax, ymin, zmax, zmin
    
#ifdef DBLPRCGRID
    real(KREAL) :: v4, x4, xini4, xfin4, y4, yini4, yfin4, z4, zini4, zfin4, drvxtot4, drvytot4, drvztot4
    integer(KINT), parameter :: i0 = 0
#else
    real(KREAL4) :: v4, x4, xini4, xfin4, y4, yini4, yfin4, z4, zini4, zfin4, drvxtot4, drvytot4, drvztot4
    integer(KINT), parameter :: i0 = 3
#endif

    namelist / options / contourval, filename, gridname, geomthr, iswindows, langstrom, lbinary, lmaxrep, &
        longoutput, lvalence, umbrlargo
    tiempo = dtime(tarray)
!    Namelist default values
    contourval = 1.d-3      ! value of density for isosurface
    filename = ""           ! root file name for output surface files (*.srf, *.sgh)
    geomthr = 1.d-5         ! Geometry threshold: two points at a distance lower than geomthr are considered to be the coincident
    gridname = ""           ! name of file with density grid
    langstrom = .true.      ! If .false. original grid distances in bohr
    lbinary = .true.        ! If true writes a file *.isoden with the surface in binary form, otherwise writes a file *.isoden_txt text mode
    lmaxrep = 5             ! highest "l" in the expansion of the density and potential
    longoutput = .false.    ! If true a more detailed output is given
    lvalence = .false.      ! If .true. only valence electrons are considered
    umbrlargo = 1.d-9       ! Threshold for determining the short-range radius
    iswindows = .false.     ! true if running on a MS-windows system
!    End of namelist defaults

    allocate(interpolmat(2,12), tritable(16,256), stat = ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Memory error when allocating interpolmat and tritable. Stop')
    endif 
    
    interpolmat = reshape((/1,2,3,1,5,6,7,5,1,2,3,4,2,3,4,4,6,7,8,8,5,6,7,8/),(/12,2/))
    tritable = reshape((/ &
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,9,4,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,10,0,0,0,0,0,0,0,0,0,0,0,0,0,                &
       2,9,4,10,9,2,0,0,0,0,0,0,0,0,0,0,2,3,11,0,0,0,0,0,0,0,0,0,0,0,0,0,1,9,4,2,3,11,0,0,0,0,0,0,0,0,0,0,              &
       10,3,11,1,3,10,0,0,0,0,0,0,0,0,0,0,3,9,4,3,11,9,11,10,9,0,0,0,0,0,0,0,4,12,3,0,0,0,0,0,0,0,0,0,0,0,0,0,          &
       1,12,3,9,12,1,0,0,0,0,0,0,0,0,0,0,2,10,1,3,4,12,0,0,0,0,0,0,0,0,0,0,2,12,3,2,10,12,10,9,12,0,0,0,0,0,0,0,        &
       4,11,2,12,11,4,0,0,0,0,0,0,0,0,0,0,1,11,2,1,9,11,9,12,11,0,0,0,0,0,0,0,4,10,1,4,12,10,12,11,10,0,0,0,0,0,0,0,    &
       10,9,11,11,9,12,0,0,0,0,0,0,0,0,0,0,5,8,9,0,0,0,0,0,0,0,0,0,0,0,0,0,5,4,1,8,4,5,0,0,0,0,0,0,0,0,0,0,             &
       1,2,10,9,5,8,0,0,0,0,0,0,0,0,0,0,5,2,10,5,8,2,8,4,2,0,0,0,0,0,0,0,2,3,11,9,5,8,0,0,0,0,0,0,0,0,0,0,              &
       4,5,8,4,1,5,2,3,11,0,0,0,0,0,0,0,10,3,11,10,1,3,9,5,8,0,0,0,0,0,0,0,3,11,10,3,10,8,3,8,4,8,10,5,0,0,0,0,         &
       9,5,8,4,12,3,0,0,0,0,0,0,0,0,0,0,12,5,8,12,3,5,3,1,5,0,0,0,0,0,0,0,10,1,2,9,5,8,3,4,12,0,0,0,0,0,0,0,            &
       5,8,12,10,5,12,10,12,3,10,3,2,0,0,0,0,4,11,2,4,12,11,8,9,5,0,0,0,0,0,0,0,2,12,11,2,5,12,2,1,5,8,12,5,0,0,0,0,    &
       5,8,9,10,1,12,10,12,11,12,1,4,0,0,0,0,5,8,12,5,12,10,10,12,11,0,0,0,0,0,0,0,10,6,5,0,0,0,0,0,0,0,0,0,0,0,0,0,    &
       10,6,5,1,9,4,0,0,0,0,0,0,0,0,0,0,1,6,5,2,6,1,0,0,0,0,0,0,0,0,0,0,9,6,5,9,4,6,4,2,6,0,0,0,0,0,0,0,                &
       2,3,11,10,6,5,0,0,0,0,0,0,0,0,0,0,4,1,9,2,3,11,5,10,6,0,0,0,0,0,0,0,6,3,11,6,5,3,5,1,3,0,0,0,0,0,0,0,            &
       3,11,6,4,3,6,4,6,5,4,5,9,0,0,0,0,10,6,5,3,4,12,0,0,0,0,0,0,0,0,0,0,1,12,3,1,9,12,5,10,6,0,0,0,0,0,0,0,           &
       1,6,5,1,2,6,3,4,12,0,0,0,0,0,0,0,3,2,6,3,6,9,3,9,12,5,9,6,0,0,0,0,11,4,12,11,2,4,10,6,5,0,0,0,0,0,0,0,           &
       5,10,6,1,9,2,9,11,2,9,12,11,0,0,0,0,6,5,1,6,1,12,6,12,11,12,1,4,0,0,0,0,6,5,9,6,9,11,11,9,12,0,0,0,0,0,0,0,      &
       10,8,9,6,8,10,0,0,0,0,0,0,0,0,0,0,10,4,1,10,6,4,6,8,4,0,0,0,0,0,0,0,1,8,9,1,2,8,2,6,8,0,0,0,0,0,0,0,             &
       2,6,4,4,6,8,0,0,0,0,0,0,0,0,0,0,10,8,9,10,6,8,11,2,3,0,0,0,0,0,0,0,11,2,3,10,6,1,6,4,1,6,8,4,0,0,0,0,            &
       9,1,3,9,3,6,9,6,8,11,6,3,0,0,0,0,3,11,6,3,6,4,4,6,8,0,0,0,0,0,0,0,8,10,6,8,9,10,4,12,3,0,0,0,0,0,0,0,            &
       10,6,8,10,8,3,10,3,1,3,8,12,0,0,0,0,3,4,12,1,2,9,2,8,9,2,6,8,0,0,0,0,12,3,2,12,2,8,8,2,6,0,0,0,0,0,0,0,          &
       10,6,9,9,6,8,11,2,4,11,4,12,0,0,0,0,6,8,1,6,1,10,8,12,1,2,1,11,12,11,1,0,12,11,1,12,1,4,11,6,1,9,1,8,6,8,1,0,    &
       12,11,6,8,12,6,0,0,0,0,0,0,0,0,0,0,11,7,6,0,0,0,0,0,0,0,0,0,0,0,0,0,1,9,4,6,11,7,0,0,0,0,0,0,0,0,0,0,            &
       10,1,2,6,11,7,0,0,0,0,0,0,0,0,0,0,2,9,4,2,10,9,6,11,7,0,0,0,0,0,0,0,2,7,6,3,7,2,0,0,0,0,0,0,0,0,0,0,             &
       2,7,6,2,3,7,4,1,9,0,0,0,0,0,0,0,10,7,6,10,1,7,1,3,7,0,0,0,0,0,0,0,6,10,9,6,9,3,6,3,7,4,3,9,0,0,0,0,              &
       3,4,12,11,7,6,0,0,0,0,0,0,0,0,0,0,12,1,9,12,3,1,11,7,6,0,0,0,0,0,0,0,1,2,10,3,4,12,6,11,7,0,0,0,0,0,0,0,         &
       6,11,7,2,10,3,10,12,3,10,9,12,0,0,0,0,7,4,12,7,6,4,6,2,4,0,0,0,0,0,0,0,1,9,12,1,12,6,1,6,2,6,12,7,0,0,0,0,       &
       4,12,7,1,4,7,1,7,6,1,6,10,0,0,0,0,7,6,10,7,10,12,12,10,9,0,0,0,0,0,0,0,6,11,7,5,8,9,0,0,0,0,0,0,0,0,0,0,         &
       5,4,1,5,8,4,7,6,11,0,0,0,0,0,0,0,2,10,1,6,11,7,9,5,8,0,0,0,0,0,0,0,11,7,6,2,10,8,2,8,4,8,10,5,0,0,0,0,           &
       7,2,3,7,6,2,5,8,9,0,0,0,0,0,0,0,2,3,6,6,3,7,4,1,5,4,5,8,0,0,0,0,9,5,8,10,1,6,1,7,6,1,3,7,0,0,0,0,                &
       8,4,10,8,10,5,4,3,10,6,10,7,3,7,10,0,4,12,3,8,9,5,11,7,6,0,0,0,0,0,0,0,6,11,7,5,8,3,5,3,1,3,8,12,0,0,0,0,        &
       1,2,10,5,8,9,3,4,12,6,11,7,0,0,0,0,10,3,2,10,12,3,10,5,12,8,12,5,6,11,7,0,9,5,8,4,12,6,4,6,2,6,12,7,0,0,0,0,     &
       6,2,12,6,12,7,2,1,12,8,12,5,1,5,12,0,1,6,10,1,7,6,1,4,7,12,7,4,9,5,8,0,7,6,10,7,10,12,5,8,10,8,12,10,0,0,0,0,    &
       11,5,10,7,5,11,0,0,0,0,0,0,0,0,0,0,5,11,7,5,10,11,1,9,4,0,0,0,0,0,0,0,11,1,2,11,7,1,7,5,1,0,0,0,0,0,0,0,         &
       9,4,2,9,2,7,9,7,5,7,2,11,0,0,0,0,2,5,10,2,3,5,3,7,5,0,0,0,0,0,0,0,4,1,9,2,3,10,3,5,10,3,7,5,0,0,0,0,             &
       1,3,5,5,3,7,0,0,0,0,0,0,0,0,0,0,9,4,3,9,3,5,5,3,7,0,0,0,0,0,0,0,11,5,10,11,7,5,12,3,4,0,0,0,0,0,0,0,             &
       1,9,3,3,9,12,5,10,11,5,11,7,0,0,0,0,4,12,3,1,2,7,1,7,5,7,2,11,0,0,0,0,7,5,2,7,2,11,5,9,2,3,2,12,9,12,2,0,        &
       10,7,5,10,4,7,10,2,4,12,7,4,0,0,0,0,9,12,2,9,2,1,12,7,2,10,2,5,7,5,2,0,4,12,7,4,7,1,1,7,5,0,0,0,0,0,0,0,         &
       7,5,9,12,7,9,0,0,0,0,0,0,0,0,0,0,8,11,7,8,9,11,9,10,11,0,0,0,0,0,0,0,1,8,4,1,11,8,1,10,11,7,8,11,0,0,0,0,        &
       11,7,8,2,11,8,2,8,9,2,9,1,0,0,0,0,11,7,8,11,8,2,2,8,4,0,0,0,0,0,0,0,2,3,7,2,7,9,2,9,10,9,7,8,0,0,0,0,            &
       3,7,10,3,10,2,7,8,10,1,10,4,8,4,10,0,8,9,1,8,1,7,7,1,3,0,0,0,0,0,0,0,8,4,3,7,8,3,0,0,0,0,0,0,0,0,0,0,            &
       3,4,12,11,7,9,11,9,10,9,7,8,0,0,0,0,3,1,8,3,8,12,1,10,8,7,8,11,10,11,8,0,2,9,1,2,8,9,2,11,8,7,8,11,3,4,12,0,     &
       12,3,2,12,2,8,11,7,2,7,8,2,0,0,0,0,9,10,7,9,7,8,10,2,7,12,7,4,2,4,7,0,1,10,2,12,7,8,0,0,0,0,0,0,0,0,0,0,         &
       8,9,1,8,1,7,4,12,1,12,7,1,0,0,0,0,8,12,7,0,0,0,0,0,0,0,0,0,0,0,0,0,8,7,12,0,0,0,0,0,0,0,0,0,0,0,0,0,             &
       4,1,9,12,8,7,0,0,0,0,0,0,0,0,0,0,1,2,10,12,8,7,0,0,0,0,0,0,0,0,0,0,9,2,10,9,4,2,12,8,7,0,0,0,0,0,0,0,            &
       11,2,3,7,12,8,0,0,0,0,0,0,0,0,0,0,2,3,11,4,1,9,7,12,8,0,0,0,0,0,0,0,3,10,1,3,11,10,7,12,8,0,0,0,0,0,0,0,         &
       7,12,8,3,11,4,11,9,4,11,10,9,0,0,0,0,8,3,4,7,3,8,0,0,0,0,0,0,0,0,0,0,8,1,9,8,7,1,7,3,1,0,0,0,0,0,0,0,            &
       3,8,7,3,4,8,1,2,10,0,0,0,0,0,0,0,2,7,3,2,9,7,2,10,9,9,8,7,0,0,0,0,11,8,7,11,2,8,2,4,8,0,0,0,0,0,0,0,             &
       11,8,7,2,8,11,2,9,8,2,1,9,0,0,0,0,1,4,8,1,8,11,1,11,10,7,11,8,0,0,0,0,8,7,11,8,11,9,9,11,10,0,0,0,0,0,0,0,       &
       7,9,5,12,9,7,0,0,0,0,0,0,0,0,0,0,4,7,12,4,1,7,1,5,7,0,0,0,0,0,0,0,9,7,12,9,5,7,10,1,2,0,0,0,0,0,0,0,             &
       10,5,7,10,7,4,10,4,2,12,4,7,0,0,0,0,7,9,5,7,12,9,3,11,2,0,0,0,0,0,0,0,2,3,11,4,1,12,1,7,12,1,5,7,0,0,0,0,        &
       5,12,9,5,7,12,1,3,10,3,11,10,0,0,0,0,11,10,4,11,4,3,10,5,4,12,4,7,5,7,4,0,9,3,4,9,5,3,5,7,3,0,0,0,0,0,0,0,       &
       1,5,3,5,7,3,0,0,0,0,0,0,0,0,0,0,2,10,1,3,4,5,3,5,7,5,4,9,0,0,0,0,2,10,5,2,5,3,3,5,7,0,0,0,0,0,0,0,               &
       9,2,4,9,7,2,9,5,7,7,11,2,0,0,0,0,11,2,1,11,1,7,7,1,5,0,0,0,0,0,0,0,5,7,4,5,4,9,7,11,4,1,4,10,11,10,4,0,          &
       11,10,5,7,11,5,0,0,0,0,0,0,0,0,0,0,5,10,6,8,7,12,0,0,0,0,0,0,0,0,0,0,1,9,4,5,10,6,12,8,7,0,0,0,0,0,0,0,          &
       6,1,2,6,5,1,8,7,12,0,0,0,0,0,0,0,12,8,7,9,4,5,4,6,5,4,2,6,0,0,0,0,10,6,5,11,2,3,8,7,12,0,0,0,0,0,0,0,            &
       7,12,8,2,3,11,1,9,4,5,10,6,0,0,0,0,8,7,12,6,5,11,5,3,11,5,1,3,0,0,0,0,4,5,9,4,6,5,4,3,6,11,6,3,12,8,7,0,         &
       8,3,4,8,7,3,6,5,10,0,0,0,0,0,0,0,10,6,5,1,9,7,1,7,3,7,9,8,0,0,0,0,4,7,3,4,8,7,2,6,1,6,5,1,0,0,0,0,               &
       7,3,9,7,9,8,3,2,9,5,9,6,2,6,9,0,10,6,5,11,2,7,2,8,7,2,4,8,0,0,0,0,2,7,11,2,8,7,2,1,8,9,8,1,10,6,5,0,             &
       5,1,11,5,11,6,1,4,11,7,11,8,4,8,11,0,8,7,11,8,11,9,6,5,11,5,9,11,0,0,0,0,7,10,6,7,12,10,12,9,10,0,0,0,0,0,0,0,   &
       4,7,12,1,7,4,1,6,7,1,10,6,0,0,0,0,1,12,9,1,6,12,1,2,6,6,7,12,0,0,0,0,7,12,4,7,4,6,6,4,2,0,0,0,0,0,0,0,           &
       2,3,11,10,6,12,10,12,9,12,6,7,0,0,0,0,1,12,4,1,7,12,1,10,7,6,7,10,2,3,11,0,12,9,6,12,6,7,9,1,6,11,6,3,1,3,6,0,   &
       7,12,4,7,4,6,3,11,4,11,6,4,0,0,0,0,6,9,10,6,3,9,6,7,3,4,9,3,0,0,0,0,10,6,7,10,7,1,1,7,3,0,0,0,0,0,0,0,           &
       2,6,9,2,9,1,6,7,9,4,9,3,7,3,9,0,2,6,7,3,2,7,0,0,0,0,0,0,0,0,0,0,2,4,7,2,7,11,4,9,7,6,7,10,9,10,7,0,              &
       11,2,1,11,1,7,10,6,1,6,7,1,0,0,0,0,1,4,9,6,7,11,0,0,0,0,0,0,0,0,0,0,11,6,7,0,0,0,0,0,0,0,0,0,0,0,0,0,            &
       12,6,11,8,6,12,0,0,0,0,0,0,0,0,0,0,12,6,11,12,8,6,9,4,1,0,0,0,0,0,0,0,6,12,8,6,11,12,2,10,1,0,0,0,0,0,0,0,       &
       11,8,6,11,12,8,10,9,2,9,4,2,0,0,0,0,12,2,3,12,8,2,8,6,2,0,0,0,0,0,0,0,1,9,4,2,3,8,2,8,6,8,3,12,0,0,0,0,          &
       10,8,6,10,3,8,10,1,3,3,12,8,0,0,0,0,8,6,3,8,3,12,6,10,3,4,3,9,10,9,3,0,3,6,11,3,4,6,4,8,6,0,0,0,0,0,0,0,         &
       9,3,1,9,6,3,9,8,6,11,3,6,0,0,0,0,10,1,2,6,11,4,6,4,8,4,11,3,0,0,0,0,10,9,3,10,3,2,9,8,3,11,3,6,8,6,3,0,          &
       2,4,6,4,8,6,0,0,0,0,0,0,0,0,0,0,1,9,8,1,8,2,2,8,6,0,0,0,0,0,0,0,10,1,4,10,4,6,6,4,8,0,0,0,0,0,0,0,               &
       10,9,8,6,10,8,0,0,0,0,0,0,0,0,0,0,6,9,5,6,11,9,11,12,9,0,0,0,0,0,0,0,6,1,5,6,12,1,6,11,12,12,4,1,0,0,0,0,        &
       1,2,10,9,5,11,9,11,12,11,5,6,0,0,0,0,11,12,5,11,5,6,12,4,5,10,5,2,4,2,5,0,3,6,2,3,9,6,3,12,9,5,6,9,0,0,0,0,      &
       1,5,12,1,12,4,5,6,12,3,12,2,6,2,12,0,1,3,6,1,6,10,3,12,6,5,6,9,12,9,6,0,10,5,6,3,12,4,0,0,0,0,0,0,0,0,0,0,       &
       3,6,11,4,6,3,4,5,6,4,9,5,0,0,0,0,6,11,3,6,3,5,5,3,1,0,0,0,0,0,0,0,4,11,3,4,6,11,4,9,6,5,6,9,1,2,10,0,            &
       6,11,3,6,3,5,2,10,3,10,5,3,0,0,0,0,9,5,6,9,6,4,4,6,2,0,0,0,0,0,0,0,1,5,6,2,1,6,0,0,0,0,0,0,0,0,0,0,              &
       9,5,6,9,6,4,10,1,6,1,4,6,0,0,0,0,10,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,5,12,8,5,10,12,10,11,12,0,0,0,0,0,0,0,         &
       1,9,4,5,10,8,10,12,8,10,11,12,0,0,0,0,2,11,12,2,12,5,2,5,1,8,5,12,0,0,0,0,4,2,5,4,5,9,2,11,5,8,5,12,11,12,5,0,   &
       5,12,8,10,12,5,10,3,12,10,2,3,0,0,0,0,10,8,5,10,12,8,10,2,12,3,12,2,1,9,4,0,12,8,5,12,5,3,3,5,1,0,0,0,0,0,0,0,   &
       12,8,5,12,5,3,9,4,5,4,3,5,0,0,0,0,3,10,11,3,8,10,3,4,8,8,5,10,0,0,0,0,10,11,8,10,8,5,11,3,8,9,8,1,3,1,8,0,       &
       4,8,11,4,11,3,8,5,11,2,11,1,5,1,11,0,2,11,3,9,8,5,0,0,0,0,0,0,0,0,0,0,5,10,2,5,2,8,8,2,4,0,0,0,0,0,0,0,          &
       5,10,2,5,2,8,1,9,2,9,8,2,0,0,0,0,5,1,4,8,5,4,0,0,0,0,0,0,0,0,0,0,5,9,8,0,0,0,0,0,0,0,0,0,0,0,0,0,                &
       10,11,9,11,12,9,0,0,0,0,0,0,0,0,0,0,4,1,10,4,10,12,12,10,11,0,0,0,0,0,0,0,1,2,11,1,11,9,9,11,12,0,0,0,0,0,0,0,   &
       4,2,11,12,4,11,0,0,0,0,0,0,0,0,0,0,2,3,12,2,12,10,10,12,9,0,0,0,0,0,0,0,4,1,10,4,10,12,2,3,10,3,12,10,0,0,0,0,   &
       1,3,12,9,1,12,0,0,0,0,0,0,0,0,0,0,4,3,12,0,0,0,0,0,0,0,0,0,0,0,0,0,3,4,9,3,9,11,11,9,10,0,0,0,0,0,0,0,           &
       10,11,3,1,10,3,0,0,0,0,0,0,0,0,0,0,3,4,9,3,9,11,1,2,9,2,11,9,0,0,0,0,2,11,3,0,0,0,0,0,0,0,0,0,0,0,0,0,           &
       2,4,9,10,2,9,0,0,0,0,0,0,0,0,0,0,1,10,2,0,0,0,0,0,0,0,0,0,0,0,0,0,1,4,9,0,0,0,0,0,0,0,0,0,0,0,0,0,               &
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0    /),(/16,256/))

    read(5,OPTIONS)
    read(5,*) projectname
    write(6,"(1x,'project name : ',a,/,1x,'==============')") trim(projectname)
    
    if (iswindows) then
        dirsep = "\\"
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        if (i .eq. 0) then    ! This is intended for MinGW, whose directory separator in windows is also /
                dirsep = "/"
                i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        endif
    else
        dirsep = "/"
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
    end if

    if (lvalence) then
        write(6,"(//'Using valence electrons only',//)")
    endif

    if (len_trim(gridname) .lt. 7) then        ! Checks name of grid for isosurface
        gridname = trim(projectname)//"-d.plt"
    else
        if (gridname(len_trim(gridname)-5:len_trim(gridname)) .ne. '-d.plt') then
            gridname = trim(gridname)//"-d.plt"
        endif
        gridname = projectname(1:i)//trim(gridname)
    endif
    
    if (len_trim(filename).ne.0) then
        filename = projectname(1:i)//trim(filename)
    endif
    
#ifdef DBLPRCGRID
    write(6,"(/'Computation in double precision',/)")
#endif

    write(6,"(/'Grid name for isosurface = ', a)") trim(gridname)
        
    inquire(file=trim(gridname), exist=lgrid, iostat=ierr)
    if (ierr .ne. 0) then
        call error (ierr, 'Grid file named '//trim(gridname)//' not found. Stop')
    endif
    iuni = 9
#if _WIN32
    open (unit=iuni, file=trim(gridname), form='binary', carriagecontrol='NONE', iostat=ierr)
#elif __INTEL_COMPILER
    open (unit=iuni, file=trim(gridname), form='binary', carriagecontrol='NONE', iostat=ierr)
#else
    open (unit=iuni, file=trim(gridname), form='unformatted', access='stream', iostat=ierr)
#endif
    if (ierr .ne. 0) then
        call error (ierr, 'Error opening grid file named '//trim(gridname)//'. Stop')
    endif
    read(iuni, iostat=ierr) iopt(:)
    if (ierr .ne. 0) then
        call error (ierr, 'Error reading iopt from file '//trim(gridname)//'. Stop')
    endif
    if (i0 .eq. 0 .and. iopt(1) .ne. 0) then
        call error (ierr, 'Program compiled in double precision and file '//trim(gridname)//' in single precision. Stop')
    else if(i0 .eq. 3 .and. iopt(1) .eq. 0) then
        call error (ierr, 'Program compiled in single precision and file '//trim(gridname)//' in double precision. Stop')
    endif
    write(6,"(/'iopt = ', 2(i5))") iopt(1:2)
    nz = iopt(3)
    ny = iopt(4)
    nx = iopt(5)
    npoints = nx * ny * nz
    write(6,"(/'Number of points in grid: nx = ', i4, ' ny = ', i4, ' nz = ', i4)") nx, ny, nz
    
        
    if (iopt(1) .eq. 3)  then   ! single precision data
        read(iuni, iostat=ierr) zini4, zfin4, yini4, yfin4, xini4, xfin4
        if (ierr .ne. 0) then
            call error (ierr, 'Error reading grid dimensions from file '//trim(gridname)//'. Stop')
        endif
        xini = xini4 ; yini = yini4 ; zini = zini4
        xfin = xfin4 ; yfin = yfin4 ; zfin = zfin4
    else if (iopt(1) .eq. 0)  then   ! double precision data
        read(iuni, iostat=ierr) zini, zfin, yini, yfin, xini, xfin
        if (ierr .ne. 0) then
            call error (ierr, 'Error reading grid dimensions from file '//trim(gridname)//'. Stop')
        endif
    else
        call error(ierr,'Wrong grid precision. Stop')
    endif
    
    if (langstrom) then
        xini = xini * angstromtobohr
        xfin = xfin * angstromtobohr
        yini = yini * angstromtobohr
        yfin = yfin * angstromtobohr
        zini = zini * angstromtobohr
        zfin = zfin * angstromtobohr
    endif
    
            
    write(6,"(/'Grid dimensions (bohr):',/3x,'xini = ', e17.10, ' xfin = ', e17.10, &
        /3x,'yini = ', e17.10, ' yfin = ', e17.10, &
        /3x,'zini = ', e17.10, ' zfin = ', e17.10, /)") xini, xfin, yini, yfin, zini, zfin
        
    dltx = (xfin - xini) / (nx - 1)
    dlty = (yfin - yini) / (ny - 1)
    dltz = (zfin - zini) / (nz - 1)
    volvoxel = dltx * dlty * dltz
    disthresq = 4.d0 * (dltx*dltx+dlty*dlty+dltz*dltz)
    if (min(dltx, dlty, dltz) .lt. 1.d-5) then
        write(6,"('Lowest dimension step size = ', e12.5)") min(dltx, dlty, dltz)
        call error (ierr, 'Wrong grid dimensions. Stop')
    endif
    
    write(6,"('dltx = ', e17.10, ' dlty = ', e17.10, ' dltz = ', e17.10)") dltx, dlty, dltz

    if (npoints .le. 0) then
        call error (ierr, 'Wrong number of grid points. Stop')
    endif
!     Read grid points
    if (iopt(1) .eq. 3)  then   ! single precision grid
        allocate (grid(npoints), grid4(npoints), stat = ierr)
        if (ierr .ne. 0) then
            call error(ierr,'Memory error when allocating grid, grid4. Stop')
        endif
        read(iuni, iostat=ierr) grid4
        if (ierr .ne. 0) then
            call error (ierr, 'Error reading grid points from file '//trim(gridname)//'. Stop')
        endif
        grid = grid4
        deallocate(grid4)
    else if (iopt(1) .eq. 0)  then   ! double precision grid
        allocate (grid(npoints), stat = ierr)
        if (ierr .ne. 0) then
            call error(ierr,'Memory error when allocating grid. Stop')
        endif
        read(iuni, iostat=ierr) grid
        if (ierr .ne. 0) then
            call error (ierr, 'Error reading grid points from file '//trim(gridname)//'. Stop')
        endif
    else
        call error(ierr,'Wrong grid precision. Stop')
    endif
    close(iuni)
!     End of grid points read
    allocate (indices(max(30000,npoints)), vertices(3,max(30000,npoints)), stat = ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Memory error when allocating indices, vertices. Stop')
    endif
    
    allocate (fvoxel(8), xvoxel(8), yvoxel(8), zvoxel(8), stat = ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Memory error when allocating fvoxel, xvoxel, yvoxel, zvoxel. Stop')
    endif
        
    call consta     !   Computes auxiliary constants
    
!     Compute surface triangles
    volume = cero
    surftot = cero
    knt = 0
    kntvert = 0
    kntind = 0
    z = zini
    do iz = 1, nz-1
        y = yini
        do iy = 1, ny-1
            x = xini
            kntgrid = (iy-1)*nx + (iz-1)*nx*ny 
            do ix = 1, nx-1
                knt = knt + 1
                kntgrid = kntgrid + 1             
                xvoxel(1) = x
                yvoxel(1) = y
                zvoxel(1) = z+dltz
                fvoxel(1) = grid(kntgrid+nx*ny)
                
                xvoxel(2) = x+dltx
                yvoxel(2) = y
                zvoxel(2) = z+dltz
                fvoxel(2) = grid(kntgrid+nx*ny+1)
                
                xvoxel(3) = x+dltx
                yvoxel(3) = y+dlty
                zvoxel(3) = z+dltz
                fvoxel(3) = grid(kntgrid+nx*ny+nx+1)
                
                xvoxel(4) = x
                yvoxel(4) = y+dlty
                zvoxel(4) = z+dltz
                fvoxel(4) = grid(kntgrid+nx*ny+nx)
                
                xvoxel(5) = x
                yvoxel(5) = y
                zvoxel(5) = z
                fvoxel(5) = grid(kntgrid)
                
                xvoxel(6) = x+dltx
                yvoxel(6) = y
                zvoxel(6) = z
                fvoxel(6) = grid(kntgrid+1)
                
                xvoxel(7) = x+dltx
                yvoxel(7) = y+dlty
                zvoxel(7) = z
                fvoxel(7) = grid(kntgrid+nx+1)
                
                xvoxel(8) = x
                yvoxel(8) = y+dlty
                zvoxel(8) = z
                fvoxel(8) = grid(kntgrid+nx)               
                icube = 1
                if (fvoxel(1).lt.contourval) icube = icube + 1
                if (fvoxel(2).lt.contourval) icube = icube + 2
                if (fvoxel(3).lt.contourval) icube = icube + 4
                if (fvoxel(4).lt.contourval) icube = icube + 8
                if (fvoxel(5).lt.contourval) icube = icube + 16
                if (fvoxel(6).lt.contourval) icube = icube + 32
                if (fvoxel(7).lt.contourval) icube = icube + 64
                if (fvoxel(8).lt.contourval) icube = icube + 128
                if (icube .eq. 1) then
                    volume = volume + volvoxel
                endif
                do k = 1, 16
                    iaux = tritable(k,icube)
                    if (iaux .le. 0) exit
                    i1 = interpolmat(iaux,1)
                    i2 = interpolmat(iaux,2)
                    xinterp = xvoxel(i1) + (contourval - fvoxel(i1)) * (xvoxel(i2)-xvoxel(i1)) / (fvoxel(i2)-fvoxel(i1))
                    yinterp = yvoxel(i1) + (contourval - fvoxel(i1)) * (yvoxel(i2)-yvoxel(i1)) / (fvoxel(i2)-fvoxel(i1))
                    zinterp = zvoxel(i1) + (contourval - fvoxel(i1)) * (zvoxel(i2)-zvoxel(i1)) / (fvoxel(i2)-fvoxel(i1))
                    lexist = .false.
                    do l = 1, kntvert
                        if (  (xinterp-vertices(1,l))*(xinterp-vertices(1,l)) &
                            + (yinterp-vertices(2,l))*(yinterp-vertices(2,l)) &
                            + (zinterp-vertices(3,l))*(zinterp-vertices(3,l)) .lt. geomthr*geomthr) then
                            kntind = kntind + 1
                            indices(kntind) = l
                            lexist = .true.
                            exit
                        endif
                    enddo
                    if (.not. lexist) then
                        kntvert = kntvert+1
                        vertices(1,kntvert) = xinterp
                        vertices(2,kntvert) = yinterp
                        vertices(3,kntvert) = zinterp
                        kntind = kntind + 1
                        indices(kntind) = kntvert
                    endif
                    xyztetr(:,mod((k-1),3)+1) = vertices(:,indices(kntind))
                    if (mod(k,3) == 0) then
                        xyztetr(:,0) = (/ xvoxel(i1), yvoxel(i1), zvoxel(i1) /)
                        a = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,1),xyztetr(:,2)-xyztetr(:,1)))
                        b = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,3),xyztetr(:,2)-xyztetr(:,3)))
                        c = sqrt(dot_product(xyztetr(:,1)-xyztetr(:,3),xyztetr(:,1)-xyztetr(:,3)))
                        s = umed * (a+b+c)                          ! Triangle semiperimeter
                        surftrian = sqrt(s*(s-a)*(s-b)*(s-c))     ! Heron's formula for the triangle surface
                        surftot = surftot + surftrian
                        xyz = (xyztetr(:,1)+xyztetr(:,2)+xyztetr(:,3)) * ri(3)
                        voltetrahed = surftrian * sqrt(dot_product(xyz-xyztetr(:,0),xyz-xyztetr(:,0))) * ri(3)
                        volume = volume + voltetrahed
                    endif
                enddo             
                x = x + dltx
            enddo
            y = y + dlty
        enddo
        z = z + dltz
    enddo
 
    write(6,"('Potential from expansion of the density: lmaxrep = ', i3)") lmaxrep

    call readdamqtisoden        !    Reads file .damqt  (generated by DAM2016)

    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zlma, zlmadx, zlmady, zlmadz. Stop')
    
    allocate (gradient(3,kntvert), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating gradient. Stop')
    
!     Writes file with the vertices and indices of MED isosurface 
    outrootname = ""
    if (len_trim(filename).ne.0) then
        outrootname = trim(filename)
    else
        outrootname = gridname(1:len_trim(gridname)-4)
    endif
    write(6,"(/3x,i9,' triangles generated for contour ',e12.5)") kntind/3, contourval 
    write(6,"(3x,'Number of vertices = ',i9)") kntvert
    write(6,"(3x,'Number of indices = ',i9)") kntind
    write(straux,"(e9.2)") contourval
    
    iuni = 10
    if (lbinary) then
        outisodenname = trim(outrootname)//"_"//trim(adjustl(straux))//".isoden"
#if _WIN32    
        open (unit=iuni+1, file=trim(outisodenname), form='binary', carriagecontrol='NONE', iostat=ierr)
        if (ierr .ne. 0) call error(ierr,'Error when opening file '//trim(outisodenname)//'. Stop')
#elif __INTEL_COMPILER
        open (unit=iuni+1, file=trim(outisodenname), form='binary', carriagecontrol='NONE', iostat=ierr)
        if (ierr .ne. 0) call error(ierr,'Error when opening file '//trim(outisodenname)//'. Stop')
#else
        open (unit=iuni+1, file=trim(outisodenname), form='unformatted', access='stream', iostat=ierr)
        if (ierr .ne. 0) call error(ierr,'Error when opening file '//trim(outisodenname)//'. Stop')
#endif
    else
        outisodenname = trim(outrootname)//"_"//trim(adjustl(straux))//".isoden_txt"
        open (unit=iuni+1, file=trim(outisodenname), iostat=ierr)
        if (ierr .ne. 0) call error(ierr,'Error when opening file '//trim(outisodenname)//'. Stop') 
    endif
    write(6,"('outisodenname = ', a)") trim(outisodenname)
    if (lbinary) then
        write(iuni+1) i0, 0, nx, ny, nz
        xini4 = xini ; xfin4 = xfin ; yini4 = yini ; yfin4 = yfin ; zini4 = zini ; zfin4 = zfin
        write(iuni+1) xini4, xfin4, yini4, yfin4, zini4, zfin4
        write(iuni+1) kntvert, kntind
    else
        write(iuni+1,"(5(i6))") i0, 0, nx, ny, nz
        write(iuni+1,"(6(1x,e14.7))") xini, xfin, yini, yfin, zini, zfin
        write(iuni+1,"(i10)") kntvert, kntind
    endif
!     Computes MED and MED gradient on the vertices
    errabs = cero
    do i = 1, kntvert
        x = vertices(1,i)
        y = vertices(2,i)
        z = vertices(3,i)
        den = cero
        drvxtot = cero
        drvytot = cero
        drvztot = cero
        do ia = 1, ncen
            call densrepr(ia, x, y, z, denrep, dendrvx, dendrvy, dendrvz)
            den = den + denrep
            drvxtot = drvxtot + dendrvx
            drvytot = drvytot + dendrvy
            drvztot = drvztot + dendrvz
        enddo
        if (abs(contourval-den) .gt. errabs) then
            errabs = abs(contourval-den)
            denmaxerr = den
        endif
!         Normalizes the gradient
        aux = sqrt(drvxtot*drvxtot+drvytot*drvytot+drvztot*drvztot)
        if (aux .lt. 1.d-12) then
            drvxtot = cero ; drvytot = cero ; drvztot = uno
        else
            drvxtot = drvxtot / aux ; drvytot = drvytot / aux ; drvztot = drvztot / aux
        endif
        gradient(:,i) = (/ drvxtot, drvytot, drvztot /)      
    enddo
    
    write(6,"(/'For density contour ', e12.5,':', /5x,'Molecular volume = ', e12.5, /5x,'Molecular surface = ', e12.5)") &
        contourval, volume, surftot
    write(6,"(/'Highest error in density, absolute = ', e10.3, ' relative = ', e10.3, ' actual value = ', e10.3,/)") &
        errabs, errabs / contourval, denmaxerr
    
!         Writes the vertex data to file
    xmax = -1.d99
    xmin =  1.d99
    ymax = -1.d99
    ymin =  1.d99
    zmax = -1.d99
    zmin =  1.d99
    do i = 1, kntvert
        x = vertices(1,i) ; y = vertices(2,i) ; z = vertices(3,i)
        xmax = max(xmax,x)
        xmin = min(xmin,x)
        ymax = max(ymax,y)
        ymin = min(ymin,y)
        zmax = max(zmax,z)
        zmin = min(zmin,z)
        drvxtot = gradient(1,i) ; drvytot = gradient(2,i) ; drvztot = gradient(3,i)  
        if (lbinary) then
            x4 = x; y4 = y; z4 = z; drvxtot4 = drvxtot; drvytot4 = drvytot; drvztot4 = drvztot;
            write(iuni+1) x4, y4, z4, drvxtot4, drvytot4, drvztot4
        else
            write(iuni+1,"(7(1x,e13.6))") x, y, z, drvxtot, drvytot, drvztot
        endif
    enddo

!     Writes triangles indices to file
    if (lbinary) then
        write(iuni+1) indices(1:kntind)
    else
        write(iuni+1,"(21(1x,i7))") indices(1:kntind)
    endif    

    
    tiempo = dtime(tarray)
    write(6,"(/1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tiempo
    stop
    end
    
!**********************************************************************
!    subroutine consta
!
!    Computes and stores auxiliary constants
!        re(i) = dfloat(i)
!        ri(i) = uno / dfloat(i)
!        fact(i) = dfloat(i!)
!        facti(i) = uno / dfloat(i!)
!        facts(i) = dfloat((i+1/2)!)
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDEN320_D
    implicit none
    integer(KINT) :: i
    real(KREAL) :: aux, ss, sd
!    auxiliary parameters and functions
    pi = acos(-uno)
    raizpi = sqrt(pi)
    re(0) = cero
    ri(0) = 1.d300
    dosl1(0) = uno
    dosl1i(0) = uno
    do i = 1, mxreal
        re(i) = re(i-1) + uno        ! dfloat(i)
        re(-i) = -re(i)
        ri(i) = uno / re(i)           ! uno / dfloat(i)
        ri(-i) = -ri(i)
        dosl1(i) = re(i) + re(i) + uno    ! dfloat(iì)
        dosl1(-i) = -re(i) - re(i) + uno
        dosl1i(i) = uno / dosl1(i)        ! dfloat( 1/(i+i+1) )
        dosl1i(-i) = uno / dosl1(-i)
    enddo
    fact(0) = uno
    facti(0) = uno
    facts(-1) = raizpi
    facts(0) = facts(-1) * umed
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)           !  i!
        facts(i) = facts(i-1) * re(i+i+1) * umed    ! (i+1/2)!
        facti(i) = uno / fact(i)             !  uno / i!
    enddo
    root(0) = cero
    do i = 1, mxroot
        root(i) = sqrt(re(i))        !  sqrt(i)
        rooti(i) = uno / root(i)     !  uno / sqrt(i)
    enddo
    return
    end
!
!    ***************************************************************
!
  subroutine readdamqtisoden
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, ncenbas, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, r, ra, ral, suml, summ
    real(KREAL) :: step, stepmed, t
    real(KREAL) :: tcheb(0:mxlenpol-1)
    inquire(file=trim(projectname)//"_2016.damqt", size=nsize, iostat=ierr)
    if (ierr .ne. 0) call error(ierr,'Error when inquiring file '//trim(projectname)//"_2016.damqt")
    if (nsize .eq. -1) call error(1,'Size of file '//trim(projectname)//"_2016.damqt cannot be determined")
    if (longoutput) write(6,"('Size of file ', a, ' = ', i12)") trim(projectname)//"_2016.damqt", nsize
#if _WIN32
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.damqt")
#elif __INTEL_COMPILER
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.damqt")
#else
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='unformatted', action = 'read', access='stream', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.damqt")
#endif
    if (longoutput) write(6,"('Open file ', a)") trim(projectname)//"_2016.damqt"
    read(10) ncen, nbas, ncaps
    nsize = nsize - sizeof(ncen) - sizeof(nbas) - sizeof(ncaps)
    write(6,"('ncen = ', i4, ' nbas = ', i6, ' nshells = ', i5)") ncen, nbas, ncaps

!    Allocates memory for geometry

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating atmnam. Stop')

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nzn. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rcen. Stop')

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zn. Stop')

    write(6,"(/24x,'GEOMETRY (BOHR)')")
    write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge')")
    do ia = 1, ncen
        read(10) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia)
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
                nzn(ia) = 0
        else
                nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
        write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5)") &
                ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia)
    enddo
    nsize = nsize - sizeof(rcen) - sizeof(zn)

!    Basis set
    read(10) lsto    ! .true. means STO basis, .false. means GTO basis
    if (lsto) then
        if (lexact) call error(1,'Exact potential not prepared for STO yet. Stop')
!        Allocates memory for the basis set

        allocate(ll(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ll. Stop')

        allocate(lmaxc(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating lmaxc. Stop')

        allocate(nf(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating nf. Stop')

        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')

        allocate(ngfin(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngfin. Stop')

        allocate(nn(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating nn. Stop')

        allocate(rnor(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating rnor. Stop')

        allocate(xx(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating xx. Stop')

        write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfin(ia)
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfin(ia))
            if (longoutput) write(6,"(t5,'center ', i4,/t12,'n', t16, 'l', t25,'exp', t35, 'ind func')") ia
            lmaxc(ia) = 0
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                rnor(i) = sqrt((dos * xx(i))**(2*nn(i)+1) / fact(2*nn(i)))
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (longoutput) write(6,"(t11,i2,t15,i2,t20,e12.5,t36,i4)") nn(i), ll(i), xx(i), nf(i)
            enddo
        enddo
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)

!        Allocates memory for the basis set

        allocate(cfcontr(nprimitot), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating cfcontr. Stop')

        allocate(ipntprim(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ipntprim. Stop')

        allocate(ll(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ll. Stop')

        allocate(lmaxc(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating lmaxc. Stop')

        allocate(ncontr(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ncontr. Stop')

        allocate(nf(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating nf. Stop')

        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')

        allocate(ngfin(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngfin. Stop')

        allocate(nprimit(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating nprimit. Stop')

        allocate(rnor(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating rnor. Stop')

        allocate(xxg(nprimitot), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating xxg. Stop')

        lmaxbase = 0
        knt = 0
        icarga = 0
        indnf = 1
        indng = 1
        ncenbas = 0
        do ia = 1, ncen
            lmaxc(ia) = 0
            read(10) ncontr(ia)
            nsize = nsize - sizeof(ncontr(ia))
            if (ncontr(ia) .le. 0) then
                ngini(ia) = -1
                ngfin(ia) = -1
                cycle
            endif
            ncenbas = ncenbas + 1
            ngini(ia) = indng
            ngfin(ia) = indng + ncontr(ia) - 1
            indng = indng + ncontr(ia)
            do j = 1, ncontr(ia)
                knt = knt + 1
                read(10) nprimit(knt), ll(knt)
                nsize = nsize - sizeof(nprimit(knt)) - sizeof(ll(knt))
                if (ll(knt) .gt. lmaxc(ia)) lmaxc(ia) = ll(knt)
                nf(knt) = indnf
                indnf = indnf + 2*ll(knt) + 1
                ipntprim(knt) = icarga+1
                read(10) xxg(icarga+1:icarga+nprimit(knt))
                nsize = nsize - sizeof(xxg(icarga+1:icarga+nprimit(knt)))
                read(10) cfcontr(icarga+1:icarga+nprimit(knt))
                nsize = nsize - sizeof(cfcontr(icarga+1:icarga+nprimit(knt)))
!                computes and stores the radial normalization factor
                aux = cero
                bux = ll(knt) + 1.5d0
                do k1 = 1, nprimit(knt)
                    do k2 = 1, k1-1
                        aux=aux + dos*cfcontr(icarga+k1)*cfcontr(icarga+k2)/(xxg(icarga+k1)+xxg(icarga+k2))**bux
                    enddo
                    aux = aux + cfcontr(icarga+k1) * cfcontr(icarga+k1) / (dos*xxg(icarga+k1))**bux
                enddo
                rnor(knt) = sqrt( dos / (facts(ll(knt))*aux) )
                icarga = icarga+nprimit(knt)    ! actualizes the index for loading primitives exponents and contraction coefficients
            enddo
        enddo
        write(6,"(/t22,'GTO Basis set',/t22,13('-'))")
    if (longoutput) then
            icarga = 0
            knt = 0
            do ia = 1, ncen
                if (ncontr(ia) .le. 0) cycle
                write(6,"(/1x,'atom no.',1x,i4,'(',a2,')')") ia, atmnam(ia)
                write(6,"(1x,'number of contractions = ',i4)") ncontr(ia)
                do j = 1, ncontr(ia)
                    knt = knt + 1
                    write(6,"(/1x,'contraction no. ',i4,' ; l = ',i2)") j,  ll(knt)
                    write(6,"('exponents: ', 8(1x,e12.5))") (xxg(icarga+k), k = 1, nprimit(knt))
                    write(6,"('coefficients: ', 8(1x,e12.5))") (cfcontr(icarga+k),k=1,nprimit(knt))
                    icarga = icarga+nprimit(knt)
                enddo
            enddo
        endif
        write(6,"('Number of basis functions = ', i4)") nbas
    endif

!    Data of density representation
    read(10) lmaxexp
    nsize = nsize - sizeof(lmaxexp)
    if (lmaxrep .gt. lmaxexp) then
        write(6,"('lmaxrep = ', i3, ' greater than lmaxexp ', i3)") lmaxrep, lmaxexp
        write(6,"('takes lmaxrep = ',i3)") lmaxexp
        lmaxrep = lmaxexp
    endif
    lmtop = (lmaxexp+1)*(lmaxexp+1)

    if (longoutput) write(6,"('lmaxexp = ', i2, ' nintervaj = ', i2)") lmaxexp, nintervaj

    allocate(icfposd(lmtop*nintervaj+1,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating icfposd. Stop')
    if (longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    nsize = nsize - sizeof(icfposd(:,1)) * ncenbas

    allocate(xajustd(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating xajustd. Stop')
    if (longoutput) write(6,"('Estimated highest size of xajustd   = ', i15, ' bytes')") size(xajustd)
    nsize = nsize - sizeof(xajustd(:,1)) * ncenbas

    allocate(cfajust(nsize/8), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating cfajust. Stop')
    if (longoutput) write(6,"('Size of cfajust   = ', i15, ' bytes')") size(cfajust)

    if (longoutput) write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
    icfposd = 0
    xajustd = cero
    cfajust = cero
    k = 0
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) cycle
        read(10) icfposd(1:lmtop*nintervaj+1,ia)
        if (k .gt. 0) then
            icfposd(1:lmtop*nintervaj+1,ia) = icfposd(1:lmtop*nintervaj+1,ia) + icfposd(lmtop*nintervaj+1,k) - 1
        endif
        k = ia
        read(10) xajustd(1:nintervaj,ia)        ! Exponents
        if (longoutput) write(6,"('fitting exponents: ',/, 8(1x,e17.10))")  xajustd(1:nintervaj,ia)
!     fitting coeficients
        read(10) cfajust(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)
    enddo

!    Generates an auxiliary index array for determining the interval to which a given r belongs
    step = rinterv(1)
    fct = uno / step
    lenindintrv = int(rinterv(nintervaj) * fct + udec)

    allocate(indintrv(lenindintrv), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating indintrv. Stop')
    if (longoutput) write(6,"('Size of indintrv   = ', i15, ' bytes')") size(indintrv)

    r = cero
    interv = 1
    do i = 1, lenindintrv-1
        r = r + step
        if (r .gt. (rinterv(interv))) interv = interv + 1
        indintrv(i) = interv
    enddo
    indintrv(lenindintrv) = interv
    
!    Determines the long-range radii and the highest l in the expansion for each interval

    allocate(rlargo(ncen), lcorto(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rlargo and lcorto. Stop')
    if (longoutput) write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)

    allocate(umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating umedpow in leedamqtden. Stop')
    umedpow(0) = uno                            !
    do i = 1, lmaxexp                            !
            umedpow(i) = umedpow(i-1) * umed            ! 1 / 2^i
    enddo
    write(6,"('Long-range threshold = ',e12.5)") umbrlargo
    nsamples = 4
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) then
                rlargo(ia) = cero
                cycle
        endif
        rlargo(ia) = rinterv(nintervaj)
        lcorto(1:nintervaj,ia) = 0
        do interv = 1, nintervaj
            dltsample = udec * (rinterv(interv) - rinterv(interv-1))
            do i = 0, nsamples-1    ! samples over nsamples points in each interval to decide the highest l
                ra = rinterv(interv-1) + dltsample + (rinterv(interv) - rinterv(interv-1) - dos * dltsample) &
                    * ri(nsamples-1) * i
                aux = exp(-xajustd(interv,ia)*ra)
                t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
                dost = t + t
                tcheb(0) = uno    ! Chebyshev T  polynomials
                tcheb(1) = t
                do j = 2, mxlenpol-1
                        tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
                enddo
                suml = cero
                lm = 0
                ral = uno
                do l = 0, lmaxrep
                    summ = cero
                    do m = -l, l
                        lm = lm + 1
                        icflm = icfposd((interv-1)*lmtop+lm,ia)
                        if(icflm .ge. icfposd((interv-1)*lmtop+lm+1,ia)) cycle
                        flm = cero
                        do j = 0, icfposd((interv-1)*lmtop+lm+1,ia)-icflm-1
                                flm   = flm + cfajust(icflm+j) * tcheb(j)
                        enddo
                        summ = summ + aux * abs(flm) * fact(l+abs(m)) * umedpow(abs(m)) * facti(l-abs(m)) &
                                    * facti(abs(m))
                    enddo
                    summ = summ * cuatro * pi * dosl1i(l) * ral
                    ral = ral * ra
                    if (summ .gt. umbrlargo .and. l .gt. lcorto(interv,ia)) lcorto(interv,ia) = l
                    suml = suml + summ
                enddo
                if (suml .gt. umbrlargo) rlargo(ia) = rinterv(interv)
            enddo
        enddo
        if (longoutput) then
            write(6,"('Long-range radius for center ',i8,' (',a2,') = ', e12.5, ' lcorto = ', 30(i3))") &
                    ia, atmnms(nzn(ia)), rlargo(ia), lcorto(1:nintervaj,ia)
        else
            write(6,"('Long-range radius for center ',i8,' (',a2,') = ', e12.5)") &
                    ia, atmnms(nzn(ia)), rlargo(ia)
        endif
    enddo
    deallocate(umedpow)
    return
    end
        
!   ***************************************************************

  subroutine densrepr(ia, x, y, z, denrep, dendrvx, dendrvy, dendrvz)
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, icflm, interv, j, jshft, kntlm, l, m
    real(KREAL) :: aux, bux, cux, denrep, dendrvx, dendrvy, dendrvz, dost, dux, drvflm
    real(KREAL) :: eux, flm, fux
    real(KREAL) :: r3inv, ra, ra2, rainv, rj2, sgn, t, umt2i, x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL) :: tcheb(0:mxlenpol-1), ucheb(0:mxlenpol-1), drvtcheb(0:mxlenpol-1), drv2tcheb(0:mxlenpol-1)
!    Contribution of atomic fragment ia to density  in point (x,y,z)
    denrep = cero
    dendrvx = cero
    dendrvy = cero
    dendrvz = cero
    if (ngini(ia) .le. 0) then
        return
    endif
    xa = x - rcen(1,ia)
    ya = y - rcen(2,ia)
    za = z - rcen(3,ia)
    ra2 = xa*xa+ya*ya+za*za
    if (ra2 .gt. rlargo(ia)*rlargo(ia)) then
        return
    endif
    ra = sqrt(ra2)
    if (ra .lt. rlargo(ia)) then
        interv = indintrv(int(fct*ra)+1)
    else
        interv = nintervaj
    endif
    rainv = uno / ra
    xadivra = xa * rainv
    yadivra = ya * rainv
    zadivra = za * rainv

    zlma(1) = uno        ! Regular spherical harmonics of r-R(ia)
    zlma(2) = ya
    zlma(3) = za
    zlma(4) = xa
    do l = 1, lcorto(interv,ia)
        zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))        ! zlm(l+1,l+1,ia)
        zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))        ! zlm(l+1,-(l+1),ia)
        zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)                ! zlm(l+1,l,ia)
        zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)                    ! zlm(l+1,-l,ia)
        do m = 0, l-1
            zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))    ! zlm(l+1,m,ia)
            zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
        enddo
    enddo
    call derivzlm(lcorto(interv,ia), idimzlm, zlma, zlmadx, zlmady, zlmadz)
    aux = exp(-xajustd(interv,ia)*ra)
    bux = xajustd(interv,ia)
    cux = dos / (rinterv(interv)-rinterv(interv-1))
    t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
    dost = t + t
    tcheb(0) = uno        ! Chebyshev T and U polynomials (U polynomial are used for derivatives of T polynomials
    tcheb(1) = t        ! and first and second derivatives of T polynomials according to:
    ucheb(0) = uno        !            D[T_n(t),t] = n * U_(n-1)(t)
    ucheb(1) = dost    !    and:        D[T_n(t),{t,2}] = (n/(t^2-1)) * ( (n-1)*t*U_(n-1) - n * U_(n-2) )
    drvtcheb(0) = cero
    drvtcheb(1) = uno
    sgn = uno
    do j = 2, mxlenpol-1
        tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
        ucheb(j) = dost * ucheb(j-1) - ucheb(j-2)
        drvtcheb(j) = re(j) * ucheb(j-1)
    enddo
    kntlm = 0
    do l = 0, lcorto(interv,ia)    !     Computes density terms 0 <= l <= lcorto(interv,ia)
        do m = -l, l
            kntlm = kntlm + 1
            icflm = icfposd(kntlm+(interv-1)*lmtop,ia)
            if(icflm .lt. icfposd(kntlm+(interv-1)*lmtop+1,ia)) then
                flm = cero
                drvflm = cero
                do j = 0, icfposd(kntlm+(interv-1)*lmtop+1,ia)-icflm-1
                    flm = flm + cfajust(j+icflm) * tcheb(j)
                    drvflm = drvflm + cfajust(j+icflm) * drvtcheb(j)
                enddo
                flm = aux * flm
                drvflm = cux * aux * drvflm        ! Converts the derivative with respect to t to derivative respect to r
                drvflm = -bux * flm + drvflm                            ! D[flm,r]
                denrep = denrep + flm * zlma(kntlm)
                dendrvx = dendrvx + xadivra * drvflm * zlma(kntlm) + flm * zlmadx(kntlm)
                dendrvy = dendrvy + yadivra * drvflm * zlma(kntlm) + flm * zlmady(kntlm)
                dendrvz = dendrvz + zadivra * drvflm * zlma(kntlm) + flm * zlmadz(kntlm)
            endif
        enddo
    enddo
    return
    end

!   ***************************************************************

  subroutine derivzlm(lmax, idimzlm, zlma, zlmadx, zlmady, zlmadz)
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm)
!    Derivatives of the regular harmonics with respecto to the Cartesian coordinates
    zlmadx(1) = cero    ! Derivatives of the S spherical harmonics
    zlmady(1) = cero
    zlmadz(1) = cero
    zlmadx(2) = cero    ! Derivatives of the P spherical harmonics
    zlmadx(3) = cero
    zlmadx(4) = zlma(1)
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    zlmadx(5) = re(6) * zlma(2)    ! Derivatives of the D spherical harmonics
    zlmadx(6) = cero
    zlmadx(7) = -zlma(4)
    zlmadx(8) = re(3) * zlma(3)
    zlmadx(9) = re(6) * zlma(4)
    zlmady(5) = re(6) * zlma(4)
    zlmady(6) = re(3) * zlma(3)
    zlmady(7) = -zlma(2)
    zlmady(8) = cero
    zlmady(9) = -re(6) * zlma(2)
    zlmadz(5) = cero
    zlmadz(6) = re(3) * zlma(2)
    zlmadz(7) = re(2) * zlma(3)
    zlmadz(8) = re(3) * zlma(4)
    zlmadz(9) = cero
    zlmadx(10) = re(15) * zlma(5)        ! Derivatives of the F spherical harmonics
    zlmadx(11) = re(10) * zlma(6)
    zlmadx(12) = -umed * zlma(5)
    zlmadx(13) = -zlma(8)
    zlmadx(14) = re(6) * zlma(7) - umed * zlma(9)
    zlmadx(15) = re(10) * zlma(8)
    zlmadx(16) = re(15) * zlma(9)
    zlmady(10) = re(15) * zlma(9)
    zlmady(11) = re(10) * zlma(8)
    zlmady(12) = re(6) * zlma(7) + umed * zlma(9)
    zlmady(13) = -zlma(6)
    zlmady(14) = -umed * zlma(5)
    zlmady(15) = -re(10) * zlma(6)
    zlmady(16) = -re(15) * zlma(5)
    zlmadz(10) = cero
    zlmadz(11) = re(5) * zlma(5)
    zlmadz(12) = re(4) * zlma(6)
    zlmadz(13) = re(3) * zlma(7)
    zlmadz(14) = re(4) * zlma(8)
    zlmadz(15) = re(5) * zlma(9)
    zlmadz(16) = cero
    do l = 4, lmax        ! Derivatives of the remaining spherical harmonics
        zlmadx(l*(l+1)+1) = - zlma((l-1)*l+2)
        zlmady(l*(l+1)+1) = - zlma((l-1)*l)
        zlmadz(l*(l+1)+1) = re(l) * zlma((l-1)*l+1)
        zlmadx(l*(l+1)+2) = umed * (re(l+1) * re(l) * zlma((l-1)*l+1) - zlma((l-1)*l+3))
        zlmadx(l*(l+1)) = umed * (- zlma((l-1)*l-1))
        zlmady(l*(l+1)+2) = -umed * (zlma((l-1)*l-1))
        zlmady(l*(l+1)) = umed * (re(l+1) * re(l) * zlma((l-1)*l+1) + zlma((l-1)*l+3))
        zlmadz(l*(l+1)+2) = re(l+1) * zlma((l-1)*l+2)
        zlmadz(l*(l+1)) = re(l+1) * zlma((l-1)*l)
        do m = 2, l-2
            zlmadx(l*(l+1)+m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l+m) - zlma((l-1)*l+m+2))
            zlmadx(l*(l+1)-m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l-m+2) - zlma((l-1)*l-m))
            zlmady(l*(l+1)+m+1) = -umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l-m+2) + zlma((l-1)*l-m))
            zlmady(l*(l+1)-m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l+m) + zlma((l-1)*l+m+2))
            zlmadz(l*(l+1)+m+1) = re(l+m) * zlma((l-1)*l+m+1)
            zlmadz(l*(l+1)-m+1) = re(l+m) * zlma((l-1)*l-m+1)
        enddo
        zlmadx(l*(l+2)) = re(l+l-1) * re(l-1) * zlma(l*l-1)
        zlmady(l*(l+2)) = -re(l+l-1) * re(l-1) * zlma(l*(l-2)+3)
        zlmadz(l*(l+2)) = re(l+l-1) * zlma(l*l)
        zlmadx(l*l+2) = re(l+l-1) * re(l-1) * zlma(l*(l-2)+3)
        zlmady(l*l+2) = re(l+l-1) * re(l-1) * zlma(l*l-1)
        zlmadz(l*l+2) = re(l+l-1) * zlma(l*(l-2)+2)
        zlmadx((l+1)*(l+1)) = re(l) * re(l+l-1) * zlma(l*l)
        zlmady((l+1)*(l+1)) = -re(l) * re(l+l-1) * zlma((l-2)*l+2)
        zlmadz((l+1)*(l+1)) = cero
        zlmadx(l*l+1) = re(l) * re(l+l-1) * zlma((l-2)*l+2)
        zlmady(l*l+1) = re(l) * re(l+l-1) * zlma(l*l)
        zlmadz(l*l+1) = cero
    enddo
    return
    end
    
!   ***************************************************************

  subroutine dzlm2y(lmax, idimzlm, zlma, zlmady, zlmadz)
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmady(idimzlm), zlmadz(idimzlm)
    zlmady(1) = cero
    zlmadz(1) = cero
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    zlmady(5) = re(6) * zlma(4)
    zlmady(6) = re(3) * zlma(3)
    zlmady(7) = -zlma(2)
    zlmady(8) = cero
    zlmady(9) = -re(6) * zlma(2)
    zlmadz(5) = cero
    zlmadz(6) = re(3) * zlma(2)
    zlmadz(7) = re(2) * zlma(3)
    zlmadz(8) = re(3) * zlma(4)
    zlmadz(9) = cero
    zlmady(10) = re(15) * zlma(9)
    zlmady(11) = re(10) * zlma(8)
    zlmady(12) = re(6) * zlma(7) + umed * zlma(9)
    zlmady(13) = -zlma(6)
    zlmady(14) = -umed * zlma(5)
    zlmady(15) = -re(10) * zlma(6)
    zlmady(16) = -re(15) * zlma(5)
    zlmadz(10) = cero
    zlmadz(11) = re(5) * zlma(5)
    zlmadz(12) = re(4) * zlma(6)
    zlmadz(13) = re(3) * zlma(7)
    zlmadz(14) = re(4) * zlma(8)
    zlmadz(15) = re(5) * zlma(9)
    zlmadz(16) = cero
    do l = 4, lmax        ! Derivatives of the remaining spherical harmonics
        zlmady(l*(l+1)+1) = - zlma((l-1)*l)
        zlmadz(l*(l+1)+1) = re(l) * zlma((l-1)*l+1)
        zlmady(l*(l+1)+2) = -umed * (zlma((l-1)*l-1))
        zlmady(l*(l+1)) = umed * (re(l+1) * re(l) * zlma((l-1)*l+1) + zlma((l-1)*l+3))
        zlmadz(l*(l+1)+2) = re(l+1) * zlma((l-1)*l+2)
        zlmadz(l*(l+1)) = re(l+1) * zlma((l-1)*l)
        do m = 2, l-2
            zlmady(l*(l+1)+m+1) = -umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l-m+2) + zlma((l-1)*l-m))
            zlmady(l*(l+1)-m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l+m) + zlma((l-1)*l+m+2))
            zlmadz(l*(l+1)+m+1) = re(l+m) * zlma((l-1)*l+m+1)
            zlmadz(l*(l+1)-m+1) = re(l+m) * zlma((l-1)*l-m+1)
        enddo
        zlmady(l*(l+2)) = -re(l+l-1) * re(l-1) * zlma(l*(l-2)+3)
        zlmadz(l*(l+2)) = re(l+l-1) * zlma(l*l)
        zlmady(l*l+2) = re(l+l-1) * re(l-1) * zlma(l*l-1)
        zlmadz(l*l+2) = re(l+l-1) * zlma(l*(l-2)+2)
        zlmady((l+1)*(l+1)) = -re(l) * re(l+l-1) * zlma((l-2)*l+2)
        zlmadz((l+1)*(l+1)) = cero
        zlmady(l*l+1) = re(l) * re(l+l-1) * zlma(l*l)
        zlmadz(l*l+1) = cero
    enddo
    return
    end
        
!   ***************************************************************

  subroutine dzlm2z(lmax, idimzlm, zlma, zlmadz)
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmadz(idimzlm)
    zlmadz(1) = cero
    do l = 1, lmax        ! Derivatives of the remaining spherical harmonics
        zlmadz(l*(l+1)+1) = re(l) * zlma((l-1)*l+1)
        zlmadz((l+1)*(l+1)) = cero
        zlmadz(l*l+1) = cero
        do m = 1, l-1
            zlmadz(l*(l+1)+m+1) = re(l+m) * zlma((l-1)*l+m+1)
            zlmadz(l*(l+1)-m+1) = re(l+m) * zlma((l-1)*l-m+1)
        enddo
    enddo
    return
    end
!
!    -------------------------------------------------------------------------------------------------------
!
  subroutine error(ierr, msg)
    USE DAM320_D
    implicit none
    integer(KINT) :: ierr
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    stop
    end
