!  Copyright 2013-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
! Program for tabulation of the electrostatic potential on an isodensity surface from the
! representation of the molecular density performed with DAM320
!
!
! Version of November 2019
!

! #define DBLPRCGRID    ! Uncomment this line  if double precision grid is wanted
!===============================================================================================
!                 MODULE DAMSGHOLE320_D
!=============================================================================================== 
MODULE DAMSGHOLE320_D
    USE DAM320_D
    USE DAM320_CONST_D
    IMPLICIT NONE
    logical :: lbinary, lcolor, lexist, lfound, lsghole
    character(300) :: gridname, indfile, outsgholename, outcolorname, outrootname, straux, vertfile
    integer(KINT) :: i, i1, i2, ia, iaux, ierr, icube, iopt(5), iuni, ix, iy, iz, j, k, knt, kntgrid
    integer(KINT) :: kntgroupmax, kntgroupmin, kntind, kntmax, kntmin, kntvert, l, maxhist, minhist, npoints, nx, ny, nz
    integer(KINT), allocatable ::  indices(:), interpolmat(:,:), kntlocal(:), localaux(:,:), localmax(:,:), localmin(:,:)
    integer(KINT), allocatable ::  nind(:), nvert(:), numlocalmax(:), numlocalmin(:), tritable(:,:)
    real(KREAL), parameter :: angstromtobohr = 1.889725989d0
    real(KREAL) :: a, aux, b, blue, bux, c, contourval, cux, d, den, denrep, dendrvx, disthresq, dlthist, dlthistinv, drvxtot
    real(KREAL) :: dendrvy, drvytot, dendrvz, drvztot, errabs, fpos, fneg, green, s, red, surfneg, surfpos, surftot, surftrian
    real(KREAL) :: thrslocal, thresmax, thresmin, topcolor, ve, vetot, vmax, vmaxabs, vmin, vminabs, vn, vntot, volume
    real(KREAL) :: voltetrahed, volvoxel, vtotia, x, xini, xinterp, xfin, y, yini, yinterp, yfin, z, zini, zinterp, zfin
    real(KREAL), allocatable :: gradient(:), grid(:), fvoxel(:), vertices(:,:), vertices2(:), vtot(:)
    real(KREAL), allocatable :: histogram(:), histpartial(:,:), xhist(:), xvoxel(:), yvoxel(:), zvoxel(:)
    real(KREAL) :: centroid(3), xyz(3), xyztetr(3,0:3)
    real(KREAL4), allocatable :: grid4(:)
    real(KREAL) :: areatot, apostot, anegtot, fdevtot, fmed, fmedneg, fmedpos, fvarneg, fvarpos, fvartot
#ifdef DBLPRCGRID
    real(KREAL) :: aux4, blue4, bux4, v4, drvxtot4, drvytot4, drvztot4, green4, red4
    real(KREAL) :: x4, xini4, xfin4, y4, yini4, yfin4, z4, zini4, zfin4 
    real(KREAL), allocatable :: gradient4(:), vertices4(:,:), vtot4(:)
    integer(KINT), parameter :: i0 = 0
#else
    real(KREAL4) :: aux4, blue4, bux4, v4, drvxtot4, drvytot4, drvztot4, green4, red4
    real(KREAL4) :: x4, xini4, xfin4, y4, yini4, yfin4, z4, zini4, zfin4
    real(KREAL4), allocatable :: gradient4(:), vertices4(:,:), vtot4(:)
    integer(KINT), parameter :: i0 = 3
#endif
END MODULE
!
!                 END OF MODULE DAMSGHOLE320_D
!...............................................................................................
  program DAMSGHOLE320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMPOT320_D
    USE DAMSGHOLE320_D
    USE PARALELO
    implicit none
    integer(KINT) :: jaux
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    logical :: lnamelist(8), ltimeprocs
    integer(KINT) :: inamelist(1)
    real(KREAL) :: xmax, xmin, ymax, ymin, zmax, zmin
    real(KREAL) :: rnamelist(5)

    
    namelist / options / contourval, filename, gridname, geomthr, iswindows, langstrom, lbinary, lexact, lmaxrep, &
        longoutput, lsghole, lvalence, thrslocal, topcolor, umbrlargo
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    abort = 0
    abortroot = 0
    tiempo = dtime(tarray)
!    Namelist default values
    contourval = 1.d-3      ! value of density for isosurface
    dlthist = 1.d-3         ! step size for histogram
    filename = ""           ! root file name for output surface files (*.srf, *.sgh)
    geomthr = 1.d-5         ! Geometry threshold: two points at a distance lower than geomthr are considered to be the coincident
    gridname = ""           ! name of file with density grid
    langstrom = .true.      ! If false original grid distances in bohr
    lbinary = .true.        ! If true writes a file *.srf with the surface in binary form, otherwise writes a file *.srf_txt text mode
    lcolor = .false.        ! If true generates a file .srf with vertices positions, normals and colors
    lexact  = .false.		! if .true. "exact" potential is tabulated
    lmaxrep = 5             ! highest "l" in the expansion of the density and potential
    longoutput = .false.    ! If true a more detailed output is given
    lsghole = .true.        ! If true generates a file *.sgh with the vertices positions and normals and mesp values on a med isosurface 
    thrslocal = 0.9d0       ! Threshold for search of local maxima and minima
    topcolor = 0.05d0       ! parameter for color assignment (blue: -topcolor, red: topcolor)
    umbrlargo = 1.d-9       ! Threshold for determining the short-range radius
    iswindows = .false.     ! true if running on a MS-windows system
!    End of namelist defaults

    if (myrank .eq. 0) write(6,"(/'Number of processors = ', i3,/)") nprocs
    ltimeprocs = .false.
    if (myrank .eq. 0) then
        allocate(timeprocs(2*nprocs), stat = ierr)
        if (ierr .eq. 0) then
            ltimeprocs = .true.
            timeprocs = 0.
        else
            write(6,"('WARNING: Memory error when allocating timeprocs, ierr =  ',i5)") ierr
            ltimeprocs = .false.
        endif
        lnamelist(1) = ltimeprocs
    endif
    CALL MPI_BCAST(lnamelist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (myrank .ne. 0) then
        ltimeprocs = lnamelist(1)
    endif
       
    allocate(interpolmat(2,12), tritable(16,256), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating interpolmat and tritable in processor ',i3)") myrank
        abort = 1
    endif

    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
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

    if (myrank .eq. 0) then   
        read(5,OPTIONS)
        read(5,*) projectname
        write(6,"(1x,'project name : ',a,/,1x,'==============')") projectname

        if((thrslocal .le. 0.2d0) .or. (thrslocal .gt. 0.99d0)) then
            write(6,"('Threshold for local extrema must be in the range (0.2,0.99]. Sets it to 0.9' )")
            thrslocal = 0.9d0
        endif
            
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

        ldengz = .false.	! Checks whether the .den file is gzipped or not
        if (lexact) then
            inquire(file=trim(projectname)//".den.gz", exist=ldengz, iostat=ierr)
            if (ierr .eq. 0 .and. ldengz) then
                    call system ("gunzip "//trim(projectname)//".den.gz")
            endif
        endif
            
#ifdef DBLPRCGRID
        write(6,"(/'Computation in double precision',/)")
#endif

        write(6,"(/'Grid name for isosurface = ', a)") trim(gridname)
        if (.not. lexact) then
            write(6,"(/'Potential from expansion of the density: lmaxrep = ', i3)") lmaxrep
        else
            write(6,"(/'Potential computed from density matrix and basis set')")
        endif
        
        lnamelist(1) = langstrom
        lnamelist(2) = lbinary
        lnamelist(3) = lcolor
        lnamelist(4) = longoutput
        lnamelist(5) = lsghole
        lnamelist(6) = iswindows
        lnamelist(7) = lvalence
        lnamelist(8) = lexact

        inamelist(1) = lmaxrep

        rnamelist(1) = contourval
        rnamelist(2) = geomthr 
        rnamelist(3) = topcolor
        rnamelist(4) = umbrlargo
        rnamelist(5) = thrslocal
    endif
    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(gridname,len(gridname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lnamelist,8,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,5,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (myrank .ne. 0) then
        langstrom  = lnamelist(1) 
        lbinary    = lnamelist(2) 
        lcolor     = lnamelist(3) 
        longoutput = lnamelist(4) 
        lsghole    = lnamelist(5) 
        iswindows  = lnamelist(6) 
        lvalence   = lnamelist(7)
        lexact     = lnamelist(8)

        lmaxrep    = inamelist(1)

        contourval = rnamelist(1)
        geomthr    = rnamelist(2) 
        topcolor   = rnamelist(3) 
        umbrlargo  = rnamelist(4) 
        thrslocal  = rnamelist(5)
    endif
        
    inquire(file=trim(gridname), exist=lgrid, iostat=ierr)
    if (myrank .eq. 0 .and. ierr .ne. 0) then
        write(6,"('Grid file named ',a,' not found in processor ',i3)") trim(gridname), myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
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
        write(6,"('Error opening grid file named ',a,' in processor ',i3)") trim(gridname), myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    read(iuni, iostat=ierr) iopt(:)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5, ' reading iopt from file ',a,' in processor ',i3)") ierr, trim(gridname), myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    if (myrank .eq. 0) then
        if (i0 .eq. 0 .and. iopt(1) .ne. 0) then
            write(6,"('Error ' i5, '. Program compiled in double precision and file ',a,' in single precision ')") &
                ierr, trim(gridname)
            abortroot = 1
        else if(i0 .eq. 3 .and. iopt(1) .eq. 0) then
            write(6,"('Error ' i5, '. Program compiled in single precision and file ',a,' in double precision ')") &
                ierr, trim(gridname)
            abortroot = 1
        endif
    endif
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    if (myrank .eq. 0) write(6,"(/'iopt = ', 2(i5))") iopt(1:2)
    nz = iopt(3)
    ny = iopt(4)
    nx = iopt(5)
    npoints = nx * ny * nz
    if (myrank .eq. 0) then
        if (npoints .le. 0) then
            write(6,"('Wrong number of grid points. Stop')")
            abortroot = 1
        endif
    endif
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    if (myrank .eq. 0) write(6,"(/'Number of points in grid: nx = ', i4, ' ny = ', i4, ' nz = ', i4)") nx, ny, nz
       
    if (iopt(1) .eq. 3)  then   ! single precision data
        read(iuni, iostat=ierr) zini4, zfin4, yini4, yfin4, xini4, xfin4
        if (ierr .ne. 0) then
            if (myrank .eq. 0) write(6,"('Error ' i5, ' reading grid dimensions from file ',a)") ierr, trim(gridname)
            abortroot = 1
        endif
        xini = xini4 ; yini = yini4 ; zini = zini4
        xfin = xfin4 ; yfin = yfin4 ; zfin = zfin4
    else if (iopt(1) .eq. 0)  then   ! double precision data
        read(iuni, iostat=ierr) zini, zfin, yini, yfin, xini, xfin
        if (ierr .ne. 0) then
            if (myrank .eq. 0) write(6,"('Error ' i5, ' reading grid dimensions from file ',a)") ierr, trim(gridname)
            abortroot = 1
        endif
    else
        if (myrank .eq. 0) write(6,"('Wrong grid precision in file ',a)") gridname
        abortroot = 1
    endif
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    if (langstrom) then
        xini = xini * angstromtobohr
        xfin = xfin * angstromtobohr
        yini = yini * angstromtobohr
        yfin = yfin * angstromtobohr
        zini = zini * angstromtobohr
        zfin = zfin * angstromtobohr
    endif
    
            
    if (myrank .eq. 0) write(6,"(/'Grid dimensions (bohr):',/3x,'xini = ', e17.10, ' xfin = ', e17.10, &
        /3x,'yini = ', e17.10, ' yfin = ', e17.10, &
        /3x,'zini = ', e17.10, ' zfin = ', e17.10, /)") xini, xfin, yini, yfin, zini, zfin
        
    dltx = (xfin - xini) / (nx - 1)
    dlty = (yfin - yini) / (ny - 1)
    dltz = (zfin - zini) / (nz - 1)
    volvoxel = dltx * dlty * dltz
    disthresq = 4.d0 * (dltx*dltx+dlty*dlty+dltz*dltz)
    if (myrank .eq. 0) write(6,*) 'disthresq = ', disthresq
    if (myrank .eq. 0) then
        if (min(dltx, dlty, dltz) .lt. 1.d-5) then
            write(6,"('Lowest dimension step size = ', e12.5)") min(dltx, dlty, dltz)
            write(6,"('Wrong grid dimensions. Stop')")
            abortroot = 1
        endif
    endif
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    
    if (myrank .eq. 0) write(6,"('dltx = ', e17.10, ' dlty = ', e17.10, ' dltz = ', e17.10)") dltx, dlty, dltz
    
    
!    Allocates buffers for para_range
    allocate(istav(0:nprocs), iendv(0:nprocs-1), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating istav and iendv in processor ',i3)") myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    call para_range(nz)
    
!     Read grid points
    if (iopt(1) .eq. 3)  then   ! single precision grid
        allocate (grid(npoints), grid4(npoints), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ' i5, ' allocating grid, grid4 in processor ', i3)") ierr, myrank
            abort= 1
        else
            read(iuni, iostat=ierr) grid4
            if (ierr .ne. 0) then
                write(6,"('Error ' i5, ' reading grid points from file ',a,' in processor ', i3)") ierr, trim(gridname), myrank
                abort= 1
            else
                grid = grid4
            endif
            deallocate(grid4)
        endif
    else if (iopt(1) .eq. 0)  then   ! double precision grid
        allocate (grid(npoints), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ' i5, ' allocating grid in processor ', i3)") ierr, myrank
            abort= 1
        else   
            read(iuni, iostat=ierr) grid
            if (ierr .ne. 0) then
                write(6,"('Error ' i5, ' reading grid points from file ',a,' in processor ', i3)") ierr, trim(gridname), myrank
                abort= 1
            endif
        endif
    else
        write(6,"('Wrong grid precision. Stop')")
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    close(iuni)
!     End of grid points read
    allocate (indices(max(30000,npoints)/nprocs), vertices(3,max(30000,npoints)/nprocs), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ' i5, ' allocating indices, vertices in processor ', i3)") ierr, myrank
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    allocate (fvoxel(8), xvoxel(8), yvoxel(8), zvoxel(8), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ' i5, ' allocating fvoxel, xvoxel, yvoxel, zvoxel in processor ', i3)") ierr, myrank
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
        
    call consta     !   Computes auxiliary constants
    
!     Compute triangles surface
    volume = cero
    surftot = cero
    knt = 0
    kntvert = 0
    kntind = 0
    z = zini + (istav(myrank)-1) * dltz
    do iz = istav(myrank), iendv(myrank)
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
  
    call readdamqtsghole        !    Reads file .damqt  (generated by DAM2016)

    if (lexact) then
        allocate(rl(-2*mxl:2*mxl,-2*mxl:2*mxl,0:2*mxl), dl(-2*mxl:2*mxl,-2*mxl:2*mxl,0:2*mxl), stat = ierr)
        if (ierr .ne. 0) then
                write(6,"('Memory error when allocating rl and dl in processor ',i3)") myrank
                abort = 1
        endif
    endif
    
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ' i5, ' allocating zlma, zlmadx, zlmady, zlmadz in processor ', i3)") ierr, myrank
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    allocate (ra2l1((lmaxrep+1)**2), ra2l1inv((lmaxrep+1)**2), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ' i5, ' allocating ra2l1 and ra2l1inv in processor ', i3)") ierr, myrank
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    allocate (gradient(3*kntvert), vtot(kntvert), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ' i5, ' allocating gradient and vtot in processor ', i3)") ierr, myrank
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

!     writes file with the values of MESP over MED isosurface 
    outrootname = ""
    if (len_trim(filename).ne.0) then
        outrootname = trim(filename)
    else
        outrootname = gridname(1:len_trim(gridname)-4)
    endif

    write(straux,"(i3)") myrank
    iuni = 10
    vertfile = trim(outrootname)//"_"//trim(adjustl(straux))//".vrttmp"
    indfile = trim(outrootname)//"_"//trim(adjustl(straux))//".indtmp"
#if _WIN32
    open (unit=iuni, file=trim(vertfile), form='binary', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(vertfile)
        abort= 1
    endif
    open (unit=iuni+1, file=trim(indfile), form='binary', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(indfile)
        abort= 1
    endif
#elif __INTEL_COMPILER
    open (unit=iuni, file=trim(vertfile), form='binary', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(vertfile)
        abort= 1
    endif
    open (unit=iuni+1, file=trim(indfile), form='binary', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(indfile)
        abort= 1
    endif
#else
    open (unit=iuni, file=trim(vertfile), form='unformatted', access='stream', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(vertfile)
        abort= 1
    endif
    open (unit=iuni+1, file=trim(indfile), form='unformatted', access='stream', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(indfile)
        abort= 1
    endif
#endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    aux4 = surftot
    bux4 = volume
    write(iuni) aux4, bux4, kntvert
    write(iuni+1) kntind
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
!     Computes MED, MED gradient and MESP on the vertices
    errabs = cero
    gradient = cero
    vtot = cero
    knt = 0
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
            if (lexact) cycle   ! If exact potential, mesp is computed outside this loop
            call mesp(ia, x, y, z, vn, ve, vtotia)
            vtot(i) = vtot(i) + vtotia
        enddo
        errabs = max(errabs, abs(contourval-den))
        if (lexact) then
            call GTOexactmesp(x, y, z, vn, ve, vtot(i))
        endif
        
!         Normalizes the gradient
        aux = sqrt(drvxtot*drvxtot+drvytot*drvytot+drvztot*drvztot)
        if (aux .lt. 1.d-12) then
            drvxtot = cero ; drvytot = cero ; drvztot = uno
        else
            drvxtot = drvxtot / aux ; drvytot = drvytot / aux ; drvztot = drvztot / aux
        endif
        gradient(3*(i-1)+1:3*i) = (/ drvxtot, drvytot, drvztot /)   
    enddo 
    aux4 = errabs
    write(iuni) aux4
!         writes the vertex data to vertices temporal file in unit iuni
    do i = 1, kntvert   
        x4 = vertices(1,i); y4 = vertices(2,i); z4 = vertices(3,i)
        drvxtot4 = gradient(3*(i-1)+1);  drvytot4 = gradient(3*(i-1)+2);  drvztot4 = gradient(3*i)
        v4 = vtot(i)
        write(iuni) x4, y4, z4, drvxtot4, drvytot4, drvztot4, v4
    enddo

!     writes triangles indices to temporal file in unit iuni+1
    write(iuni+1) indices(1:kntind)

    close(iuni)
    close(iuni+1)
    deallocate (gradient, indices, vertices, vtot)
    if (myrank .ne. 0) tiempo = dtime(tarray)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! Necessary for the processors to end writing temporal files
!     Gathers files generated by processors
    if (myrank .eq. 0) then
        if (ldengz) then	! restores density file back to its original gzipped status
            call system ("gzip "//trim(projectname)//".den")
        endif
        kntvert = 0
        kntind = 0
        errabs = cero
        vmax = -1.d-90
        vmaxabs = -1.d-90
        vmin = 1.d90
        vminabs = 1.d90
        surftot = cero
        volume = cero
!     reads temporal files and writes file with the values of MESP over MED isosurface 
        outrootname = ""
        if (len_trim(filename).ne.0) then
            outrootname = trim(filename)
        else
            outrootname = gridname(1:len_trim(gridname)-4)
        endif
!         Opens temporal files for reading
        allocate (nind(0:nprocs), nvert(0:nprocs), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ' i5, ' allocating nind and nvert in processor ', i3)") ierr, myrank
            abort= 1
        endif
        if (abort .eq. 0) then
            iuni = 10
            do i = 0, nprocs-1
                write(straux,"(i3)") i
                vertfile = trim(outrootname)//"_"//trim(adjustl(straux))//".vrttmp"
                indfile = trim(outrootname)//"_"//trim(adjustl(straux))//".indtmp"
#if _WIN32
                open (unit=iuni+2*i, file=trim(vertfile), form='binary', carriagecontrol='NONE', iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' opening file ',a)") ierr, trim(vertfile)
                    abort= 1
                    exit
                endif
                if (abort .eq. 0) then
                    open (unit=iuni+2*1+1, file=trim(indfile), form='binary', carriagecontrol='NONE', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(indfile)
                        abort= 1
                        exit
                    endif
                endif
#elif __INTEL_COMPILER
                open (unit=iuni+2*i, file=trim(vertfile), form='binary', carriagecontrol='NONE', iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' opening file ',a)") ierr, trim(vertfile)
                    abort= 1
                    exit
                endif
                if (abort .eq. 0) then
                    open (unit=iuni+2*i+1, file=trim(indfile), form='binary', carriagecontrol='NONE', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(indfile)
                        abort= 1
                        exit
                    endif
                endif
#else
                open (unit=iuni+2*i, file=trim(vertfile), form='unformatted', access='stream', iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' opening file ',a)") ierr, trim(vertfile)
                    abort= 1
                    exit
                endif
                if (abort .eq. 0) then
                    open (unit=iuni+2*i+1, file=trim(indfile), form='unformatted', access='stream', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(indfile)
                        abort= 1
                        exit
                    endif
                endif
#endif
                if (abort .eq. 0) then
                    read(iuni+2*i) aux4, bux4, i1
                    surftot = surftot + aux4
                    volume = volume + bux4
                    nvert(i) = i1
                    kntvert = kntvert + i1
                    read(iuni+2*i) aux4
                    errabs = max(errabs,aux4)
                    read(iuni+2*i+1) i2
                    nind(i) = i2
                    kntind = kntind + i2
                endif
            enddo
        endif
        if (abort .eq. 0) then
            allocate (indices(kntind), gradient4(3*kntvert), vertices4(3,kntvert), vtot4(kntvert), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Error ' i5, ' allocating indices, gradient4, vertices4, and vtot4 in processor ', i3)")ierr, myrank
                abort= 1
            endif
        endif
        if (abort .eq. 0) then
            i1 = 0
            i2 = 0
            knt = 0
            do i = 0, nprocs-1
                do j = 1, nvert(i)
                    read(iuni+2*i) vertices4(:,i1+1), gradient4(3*i1+1:3*i1+3), vtot4(i1+1)
                    vmax = max(vmax,vtot4(i1+1))
                    vmin = min(vmin,vtot4(i1+1))
                    i1 = i1 + 1
                enddo
                read(iuni+2*i+1) indices(i2+1:i2+nind(i))
                indices(i2+1:i2+nind(i)) = indices(i2+1:i2+nind(i)) + knt
                i2 = i2 + nind(i)
                knt = knt + nvert(i)
            enddo
            
            outrootname = ""
            if (len_trim(filename).ne.0) then
                outrootname = trim(filename)
            else
                outrootname = gridname(1:len_trim(gridname)-4)
            endif
            write(straux,"(e9.2)") contourval
            if (lexact) then
                straux = trim(straux)//"_exact"
            endif

            if (lbinary) then
                if (lcolor) outcolorname = trim(outrootname)//"_"//trim(adjustl(straux))//".srf"
                if (lsghole) outsgholename = trim(outrootname)//"_"//trim(adjustl(straux))//".sgh"
#if _WIN32
                if (lcolor) then
                    open (unit=iuni+2*nprocs, file=trim(outcolorname), form='binary', carriagecontrol='NONE', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outcolorname)
                        abort= 1
                    endif
                endif
            
                if (abort .eq. 0 .and. lsghole) then
                    open (unit=iuni+2*nprocs+1, file=outsgholename, form='binary', carriagecontrol='NONE', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outsgholename)
                        abort= 1
                    endif
                endif
#elif __INTEL_COMPILER
                if (lcolor) then
                    open (unit=iuni+2*nprocs, file=trim(outcolorname), form='binary', carriagecontrol='NONE', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outcolorname)
                        abort= 1
                    endif
                endif
            
                if (abort .eq. 0 .and. lsghole) then
                    open (unit=iuni+2*nprocs+1, file=outsgholename, form='binary', carriagecontrol='NONE', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outsgholename)
                        abort= 1
                    endif
                endif
#else
                if (lcolor) then
                    open (unit=iuni+2*nprocs, file=trim(outcolorname), form='unformatted', access='stream', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outcolorname)
                        abort= 1
                    endif
                endif
            
                if (abort .eq. 0 .and. lsghole) then
                    open (unit=iuni+2*nprocs+1, file=outsgholename, form='unformatted', access='stream', iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' reading file ',a)") ierr, trim(outsgholename)
                        abort= 1
                    endif
                endif
#endif
            else
                if (lcolor) then
                    outcolorname = trim(outrootname)//"_"//trim(adjustl(straux))//".srf_txt"
                    open (unit=iuni+2*nprocs, file=trim(outcolorname), iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' reading file ',a)") ierr, trim(outcolorname)
                        abort= 1
                    endif
                endif
                if (abort .eq. 0 .and. lsghole) then
                    outsgholename = trim(outrootname)//"_"//trim(adjustl(straux))//".sgh_txt"
                    open (unit=iuni+2*nprocs+1, file=outsgholename, iostat=ierr)
                    if (ierr .ne. 0) then
                        write(6,"('Error ', i5,' reading file ',a)") ierr, trim(outsgholename)
                        abort= 1
                    endif
                endif
            endif
        endif
        
        if (abort .eq. 0) then
            if (lbinary) then
                if (lcolor) then
                    write(iuni+2*nprocs) i0, 0, nx, ny, nz
                    xini4 = xini ; xfin4 = xfin ; yini4 = yini ; yfin4 = yfin ; zini4 = zini ; zfin4 = zfin
                    write(iuni+2*nprocs) xini4, xfin4, yini4, yfin4, zini4, zfin4
                    write(iuni+2*nprocs) kntvert, kntind
                endif
                if (lsghole) then
                    write(iuni+2*nprocs+1) i0, 0, nx, ny, nz
                    xini4 = xini ; xfin4 = xfin ; yini4 = yini ; yfin4 = yfin ; zini4 = zini ; zfin4 = zfin
                    write(iuni+2*nprocs+1) xini4, xfin4, yini4, yfin4, zini4, zfin4
                    write(iuni+2*nprocs+1) kntvert, kntind
                endif
            else
                if (lcolor) then
                    write(iuni+2*nprocs,"(5(i6))") i0, 0, nx, ny, nz
                    write(iuni+2*nprocs,"(6(1x,e14.7))") xini, xfin, yini, yfin, zini, zfin
                    write(iuni+2*nprocs,"(i10)") kntvert, kntind
                endif
                if (lsghole) then
                    write(iuni+2*nprocs+1,"(5(i6))") i0, 0, nx, ny, nz
                    write(iuni+2*nprocs+1,"(6(1x,e14.7))") xini, xfin, yini, yfin, zini, zfin
                    write(iuni+2*nprocs+1,"(i10)") kntvert, kntind
                endif
            endif
!         writes the vertices data to file
            xmax = -1.d99
            xmin =  1.d99
            ymax = -1.d99
            ymin =  1.d99
            zmax = -1.d99
            zmin =  1.d99
            do i = 1, kntvert
                x = vertices4(1,i) ; y = vertices4(2,i) ; z = vertices4(3,i)
                xmax = max(xmax,x)
                xmin = min(xmin,x)
                ymax = max(ymax,y)
                ymin = min(ymin,y)
                zmax = max(zmax,z)
                zmin = min(zmin,z)
                drvxtot = gradient4(3*(i-1)+1) ; drvytot = gradient4(3*(i-1)+2) ; drvztot = gradient4(3*i)     
                if (lbinary) then
                    x4 = x; y4 = y; z4 = z; drvxtot4 = drvxtot; drvytot4 = drvytot; drvztot4 = drvztot;
                    if (lcolor) then
                        aux = vtot4(i) 
                        call color_generator(aux, topcolor, red, green, blue)
                        red4 = red; green4 = green; blue4 = blue
                        write(iuni+2*nprocs) x4, y4, z4, drvxtot4, drvytot4, drvztot4, red4, green4, blue4
                    endif
                    if (lsghole) then
                        write(iuni+2*nprocs+1) x4, y4, z4, drvxtot4, drvytot4, drvztot4, vtot4(i)
                    endif
                else
                    if (lcolor) then
                        call color_generator(vtot(i), topcolor, red, green, blue)
                        write(iuni+2*nprocs,"(9(1x,e13.6))") x, y, z, drvxtot, drvytot, drvztot, red, green, blue
                    endif
                    if (lsghole) write(iuni+2*nprocs+1,"(7(1x,e13.6))") x, y, z, drvxtot, drvytot, drvztot, vtot4(i)
                endif
            enddo
        endif
        
!     Writes triangles indices to file
        if (abort .eq. 0) then
            if (lbinary) then
                if (lcolor) write(iuni+2*nprocs) indices(1:kntind)
                if (lsghole) write(iuni+2*nprocs+1) indices(1:kntind)
            else
                if (lcolor) write(iuni+2*nprocs,"(21(1x,i7))") indices(1:kntind)
                if (lsghole) write(iuni+2*nprocs+1,"(21(1x,i7))") indices(1:kntind)
            endif
            write(6,"(/3x,i9,' triangles generated for contour ',e12.5)") kntind/3, contourval 
            write(6,"(3x,'Number of vertices = ',i9)") kntvert
            write(straux,"(e9.2)") contourval
            write(6,"(/'For density contour ', e12.5,':', &
                &/5x,'Molecular volume / bohr^3 = ', e12.5, 5x,' Molecular volume / A^3 = ', e12.5, &
                &/5x,'Molecular surface / bohr^2 = ', e12.5, 5x,' Molecular surface / A^2 = ', e12.5)") &
                contourval, volume, volume*0.148184534296, surftot, surftot*0.280028297329d0
            write(6,"(/'Highest error in density, absolute = ', e10.3, ' relative = ', e10.3,/)") errabs, errabs / contourval 
            write(6,"('Highest value of potential = ', e12.5)") vmax
            write(6,"('Lowest value of potential  = ', e12.5)") vmin
    
!     Prepares data for histograms

            dlthistinv = uno / dlthist
            minhist = int(vmin * dlthistinv)-1
            maxhist = int(vmax  * dlthistinv)+2
            allocate (xhist(minhist:maxhist), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Error ', i5,' when allocating xhist')") ierr
                abort = 1
            endif
        endif
        if (abort .eq. 0) then
            xhist(minhist) = minhist * dlthist
            do i = minhist+1, maxhist
                xhist(i) = xhist(i-1) + dlthist
            enddo            
             
!     Search regions with local maxima and minima

            thresmax = thrslocal * vmax
            thresmin = thrslocal * vmin
!            topcolor = max(topcolor, thresmax, abs(thresmin))
            call localextrema
        endif
        if (abort .eq. 0) then
            allocate (histpartial(minhist:maxhist,kntgroupmax+kntgroupmin), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Error ', i5,' when allocating histpartial')") ierr
                abort = 1
            endif
        endif
        if (abort .eq. 0) then
            histpartial = 0

            write(6,"(//30('*'),2x,'FINAL RESULTS ON LOCAL EXTREMA',2x,30('*'),/)")

            write(6,"(/'Number of local maxima (higher than mesp = ', e9.2, ') found = ', i3)") thresmax, kntgroupmax
            
!     Writes number of regions with local maxima to file
            if (lbinary) then
                if (lcolor) write(iuni+2*nprocs) kntgroupmax
                if (lsghole) write(iuni+2*nprocs+1) kntgroupmax
            else
                if (lcolor) write(iuni+2*nprocs,"(i7)") kntgroupmax
                if (lsghole) write(iuni+2*nprocs,"(i7)") kntgroupmax
            endif
            if (kntgroupmax .gt. 0) then
                write(6,"(/,'Positions and mesp values of the local maxima in the mesh, and area for local surface with mesp > '&
                    &,e12.5,/111('-'))") thresmax
                do jaux = 1, kntgroupmax
                    surfpos = cero
                    call makehistmax(jaux)
                    vmax = 0.d0
                    i1 = 0
                    do i = 1, numlocalmax(jaux)
                        do k = 0, 2
                            if (vtot4(indices(localmax(i,jaux)+k)) .gt. vmax) then
                                i1 = localmax(i,jaux)+k
                                vmax = vtot4(indices(localmax(i,jaux)+k))
                            endif
                        enddo
                    enddo
                    vmaxabs = max(vmaxabs,vmax)
                    if (i1 .lt. 1) cycle
                    write(6,"('x = ', e12.5, ' y = ', e12.5, ' z = ', e12.5, ' mesp = ', e12.5, ' area / bohr^2 = ', &
                        e12.5, ' area / A^2 = ', e12.5)") vertices4(:,indices(i1)), vmax, surfpos, surfpos * 0.280028297329d0
!     Writes position of local maximum and value to file
                    if (lbinary) then
                        v4 = vmax
                        if (lcolor) write(iuni+2*nprocs) vertices4(:,indices(i1)), v4
                        if (lsghole) write(iuni+2*nprocs+1) vertices4(:,indices(i1)), v4
                    else
                        if (lcolor) write(iuni+2*nprocs,"(4(e12.5,1x))") vertices4(:,indices(i1)), vmax
                        if (lsghole) write(iuni+2*nprocs+1,"(4(e12.5,1x))") vertices4(:,indices(i1)), vmax
                    endif
                enddo
            endif
                  
!     Searches local minima in the mesh and computes local histogram

            write(6,"(/'Number of local minima (lower than mesp = ', e9.2, ') found = ', i3)") thresmin, kntgroupmin
!write(6,*) 'kntgroupmin = ', kntgroupmin
!     Writes number of regions with local minima to file
            if (lbinary) then
                if (lcolor) write(iuni+2*nprocs) kntgroupmin
                if (lsghole) write(iuni+2*nprocs+1) kntgroupmin
            else
                if (lcolor) write(iuni+2*nprocs,"(i7)") kntgroupmin
                if (lsghole) write(iuni+2*nprocs,"(i7)") kntgroupmin
            endif
            if (kntgroupmin .gt. 0) then
                write(6,"(/,'Positions and mesp values of the local minima in the mesh, and area for local surface with mesp < '&
                    &,e12.5,/111('-'))") thresmin
                do jaux = 1, kntgroupmin
                    surfneg = cero
                    call makehistmin(jaux)
                    vmin = 0.d0
                    i1 = 0
                    do i = 1, numlocalmin(jaux)
                        do k = 0, 2
                            if (vtot4(indices(localmin(i,jaux)+k)) .lt. vmin) then
                                i1 = localmin(i,jaux)+k
                                vmin = vtot4(indices(localmin(i,jaux)+k))
                            endif
                        enddo
                    enddo
                    vminabs = min(vminabs,vmin)
                    if (i1 .lt. 1) cycle
                    write(6,"('x = ', e12.5, ' y = ', e12.5, ' z = ', e12.5, ' mesp = ', e12.5, ' area / bohr^2 = ', &
                        e12.5, ' area / A^2 = ', e12.5)") vertices4(:,indices(i1)), vmin, surfneg, surfneg * 0.280028297329d0
!     Writes position of local minimum and value to file
                    if (lbinary) then
                        v4 = vmin
                        if (lcolor) write(iuni+2*nprocs) vertices4(:,indices(i1)), v4
                        if (lsghole) write(iuni+2*nprocs+1) vertices4(:,indices(i1)), v4
                    else
                        if (lcolor) write(iuni+2*nprocs,"(4(e12.5,1x))") vertices4(:,indices(i1)), vmin
                        if (lsghole) write(iuni+2*nprocs+1,"(4(e12.5,1x))") vertices4(:,indices(i1)), vmin
                    endif
                enddo
            endif
            do i = 0, nprocs
                if (lcolor) close(iuni+2*i)
                if (lsghole) close(iuni+2*i+1)
            enddo
            
            surfpos = cero
            surfneg = cero
    
!     Computes global histogram
    
            call makehistogram
    
!         Writes the mesp total and partial histograms  to file

            outrootname = ""
            if (len_trim(filename).ne.0) then
                outrootname = trim(filename)
            else
                outrootname = gridname(1:len_trim(gridname)-4)
            endif
            write(straux,"(e9.2)") contourval
            if (lexact) then
                straux = trim(straux)//"_exact"
            endif

            outcolorname = trim(outrootname)//"_"//trim(adjustl(straux))//".hst"
            open (unit=iuni+2, file=trim(outcolorname), iostat=ierr)
            if (ierr .ne. 0) then
                write(6,"('Error ', i5,' when opening file ', a)") ierr, trim(outcolorname)
                abort = 1
            endif
        endif
        if (abort .eq. 0) then
            write(iuni+2,"(i10)") maxhist-minhist+1
            write(iuni+2,"(10(1x,e12.5))") xhist
            write(iuni+2,"(10(1x,e12.5))") histogram
            write(iuni+2,"(i10)") kntgroupmax+kntgroupmin
            do i = 1, kntgroupmax+kntgroupmin
                write(iuni+2,"(10(1x,e12.5))") histpartial(:,i)
            enddo
            close(iuni+2)
            write(6,"(/5x,'Total surface with mesp higher than ', e12.5, ' = ', e12.5, 3x,'(',e9.2,'%)')") &
                thresmax, surfpos, 100.d0 * surfpos / surftot
            write(6,"(5x,'Total surface with mesp lower than  ', e12.5, ' = ', e12.5, 3x,'(',e9.2,'%)')") &
                thresmin, surfneg, 100.d0 * surfneg / surftot
        endif
            
!     Computes the mean and variances of positive, negative and total MESP

        call statMESP
           
        write(6,"(//30('*'),2x,'MESP statistics',2x,30('*'),/)")
        write(6,"(/5x,'Total surface with positive MESP = ', e12.5)") apostot
        write(6,"( 5x,'Total surface with negative MESP = ', e12.5)") anegtot
        write(6,"(/5x,'MESP average = ', e12.5)") fmed
        write(6,"( 5x,'Positive MESP average = ', e12.5)") fmedpos
        write(6,"( 5x,'Negative MESP average = ', e12.5)")  fmedneg
        write(6,"(/5x,'MESP variance = ', e12.5)") fvartot
        write(6,"( 5x,'Positive MESP variance = ', e12.5)") fvarpos
        write(6,"( 5x,'Negative MESP variance = ', e12.5)") fvarneg
        write(6,"(/5x,'MESP average deviation = ', e12.5)") fdevtot
        write(6,"(/5x,'MESP nu parameter = ', e12.5)") fvarpos * fvarneg / (fvartot**2)
        
        outrootname = ""
        if (len_trim(filename).ne.0) then
            outrootname = trim(filename)
        else
            outrootname = gridname(1:len_trim(gridname)-4)
        endif
        vertfile = trim(outrootname)//"_*.vrttmp"
        indfile = trim(outrootname)//"_*.indtmp"
        call system("rm -f "//trim(vertfile))
        call system("rm -f "//trim(indfile))
        outrootname = ""
        if (len_trim(filename).ne.0) then
            outrootname = trim(filename)
        else
            outrootname = gridname(1:len_trim(gridname)-6)
        endif
        write(straux,"(e9.2)") contourval
        if (lexact) then
            straux = trim(straux)//"_exact"
        endif
        outcolorname = trim(outrootname)//"_"//trim(adjustl(straux))//"_SGMESP_summary.txt"
        open (unit=69, file=trim(outcolorname), iostat=ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ', i5,' when opening file ', a)") ierr, trim(outcolorname)
            abort = 1
        endif
        if (abort .eq. 0) then
            write(69,*) surftot
            write(69,*) volume
            write(69,*) vmaxabs
            write(69,*) vminabs
            write(69,*) kntgroupmax
            write(69,*) kntgroupmin
            write(69,*) xmin
            write(69,*) xmax
            write(69,*) ymin
            write(69,*) ymax
            write(69,*) zmin
            write(69,*) zmax
            write(69,*) apostot
            write(69,*) anegtot
            write(69,*) fmed
            write(69,*) fmedpos
            write(69,*) fmedneg
            write(69,*) fvartot
            write(69,*) fvarpos
            write(69,*) fvarneg
            write(69,*) fdevtot
            write(69,*) fvarpos * fvarneg / (fvartot**2)
            close(69)
        endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (myrank .eq. 0) tiempo = dtime(tarray)
    if (ltimeprocs) then
        CALL MPI_GATHER(tarray, 2, mpi_real4, timeprocs, 2, mpi_real4, 0, MPI_COMM_WORLD, ierr)
        if (ierr .eq. 0 .and. myrank .eq. 0) then
            write(6,"(/30x,'TIMING (in seconds)',/)")
            do i = 0, nprocs-1
                write(6,"(1x,'Processor ', i2, ' (user, system, total):',5x,'(',e12.5,',',e12.5,',',e12.5')')") &
                    i, timeprocs(2*i+1), timeprocs(2*i+2), timeprocs(2*i+1)+timeprocs(2*i+2)
            enddo
            write(6,"(' ')")
        endif
    endif
    
    call MPI_FINALIZE(ierr)
    stop
    end
!
!    -------------------------------------------------------------------------------------------------------
!
  subroutine para_range(nz)
    USE DAM320_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, naux, nz
    naux = nz / nprocs
    istav(0) = 1
    iendv(0) = nz - naux * (nprocs-1)-1
    do i = 1, nprocs-1
        istav(i) = iendv(i-1) + 1
        iendv(i) = istav(i) + naux - 1
    enddo
    return
    end
!
!**********************************************************************
!
  subroutine localextrema
    USE DAMSGHOLE320_D, vold => vtot, vtot => vtot4, verticesold => vertices, vertices => vertices4
    USE PARALELO
    implicit none   
    integer(KINT) :: ii, jj, jaux, kntcoincid, maxaux, minaux
    real(KREAL), allocatable :: vmaxarray(:), vminarray(:)
    real(KREAL) :: vauxmax, vauxmin
    kntmax = 0
    kntmin = 0
    do i = 1, kntind, 3
        if (max(vtot(indices(i)), vtot(indices(i+1)), vtot(indices(i+2))) .gt. thresmax) kntmax = kntmax + 1
        if (min(vtot(indices(i)), vtot(indices(i+1)), vtot(indices(i+2))) .lt. thresmin) kntmin = kntmin + 1
    enddo
    allocate (localmax(max(1,kntmax),2), localmin(max(1,kntmin),2), numlocalmax(max(1,kntmax)), &
        numlocalmin(max(1,kntmin)), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ', i5,' when allocating localmax, localmin, numlocalmax and numlocalmin')") ierr
        abort = 1
        return
    endif
    kntmax = 0
    kntmin = 0
    kntgroupmax = 1
    kntgroupmin = 1
    numlocalmax = 0    
    numlocalmin = 0
    do i = 1, kntind, 3
        if (max(vtot(indices(i)), vtot(indices(i+1)), vtot(indices(i+2))) .gt. thresmax) then
            kntmax = kntmax + 1
            centroid = ri(3) * ( vertices(:,indices(i)) + vertices(:,indices(i+1)) + vertices(:,indices(i+2)) ) 
            if (kntmax .eq. 1) then
                localmax(kntmax,1) = i
                localmax(kntmax,2) = 1
                numlocalmax(1) = 1
            else
                lfound = .false.
                do j = 1, kntmax-1
                    xyz = ri(3) * ( vertices(:,indices(localmax(j,1))) + vertices(:,indices(localmax(j,1)+1)) &
                        + vertices(:,indices(localmax(j,1)+2)) )
                    if ( dot_product(centroid-xyz,centroid-xyz) .le. disthresq ) then
                        localmax(kntmax,1) = i
                        localmax(kntmax,2) = localmax(j,2)
                        numlocalmax(localmax(j,2)) = numlocalmax(localmax(j,2)) + 1
                        lfound = .true.
                        exit
                    endif
                enddo
                if (.not. lfound) then
                    kntgroupmax = kntgroupmax + 1
                    localmax(kntmax,1) = i
                    localmax(kntmax,2) = kntgroupmax
                    numlocalmax(kntgroupmax) = 1
                endif
            endif
        else if (min(vtot(indices(i)), vtot(indices(i+1)), vtot(indices(i+2))) .lt. thresmin) then
            kntmin = kntmin + 1
            centroid = ri(3) * ( vertices(:,indices(i)) + vertices(:,indices(i+1)) + vertices(:,indices(i+2)) ) 
            if (kntmin .eq. 1) then
                localmin(kntmin,1) = i
                localmin(kntmin,2) = 1
                numlocalmin(1) = 1
            else
                lfound = .false.
                do j = 1, kntmin-1
                    xyz = ri(3) * ( vertices(:,indices(localmin(j,1))) + vertices(:,indices(localmin(j,1)+1)) &
                        + vertices(:,indices(localmin(j,1)+2)) )
                    if ( dot_product(centroid-xyz,centroid-xyz) .le. disthresq ) then
                        localmin(kntmin,1) = i
                        localmin(kntmin,2) = localmin(j,2)
                        numlocalmin(localmin(j,2)) = numlocalmin(localmin(j,2)) + 1
                        lfound = .true.
                        exit
                    endif
                enddo
                if (.not. lfound) then
                    kntgroupmin = kntgroupmin + 1
                    localmin(kntmin,1) = i
                    localmin(kntmin,2) = kntgroupmin
                    numlocalmin(kntgroupmin) = 1
                endif
            endif
        endif
    enddo
    if (kntmax .eq. 0) kntgroupmax = 0
    if (kntmin .eq. 0) kntgroupmin = 0

!         Scan regions of local maxima to gather overlapping regions

    if (kntmax .gt. 0) then
        write(6,"(/'Number of starting regions with local maxima = ', i3)") kntgroupmax
        call flush(6)
        do i = 1, kntgroupmax
            write(6,"('Number triangles in region ', i3 '  = ', i7)") i, numlocalmax(i)
        enddo
    else
        write(6,"(/'No regions with local maxima')")
    endif

    if (kntmax .gt. 0) then
        allocate (localaux(kntmax,kntgroupmax), kntlocal(kntgroupmax), vmaxarray(kntgroupmax), &
                            stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ', i5,' when allocating localaux, kntlocal and vmaxarray')") ierr
            abort = 1
            return
        endif
        kntlocal = 0
        do i = 1, kntmax
            kntlocal(localmax(i,2)) = kntlocal(localmax(i,2))+1
            localaux(kntlocal(localmax(i,2)),localmax(i,2)) = localmax(i,1)
        enddo


        write(6,"(/,'Positions and mesp values of the local maxima in the starting regions for local surface with mesp > ',&
            &e12.5)") thresmax
        vmaxarray = cero
        do jaux = 1, kntgroupmax
            vmax = 0.d0
            i1 = 0
            do i = 1, kntlocal(jaux)
                do k = 0, 2
                    if (vtot(indices(localaux(i,jaux)+k)) .gt. vmax) then
                        i1 = localaux(i,jaux)+k
                        vmax = vtot(indices(localaux(i,jaux)+k))
                    endif
                enddo
            enddo
            vmaxarray(jaux) = vmax
            if (i1 .lt. 1) cycle
            write(6,"('x = ', e12.5, ' y = ', e12.5, ' z = ', e12.5, ' mesp = ', e12.5, e12.5)") vertices(:,indices(i1)), vmax
        enddo

        do j = 2, kntgroupmax
            dok1: do k = 1, numlocalmax(j)
                centroid = ri(3) * ( vertices(:,indices(localaux(k,j))) + vertices(:,indices(localaux(k,j)+1)) &
                    + vertices(:,indices(localaux(k,j)+2)) )
                do i = 1, j-1
                    do l = 1, numlocalmax(i)
                        xyz = ri(3) * ( vertices(:,indices(localaux(l,i))) + vertices(:,indices(localaux(l,i)+1)) &
                            + vertices(:,indices(localaux(l,i)+2)) )
                        if ( dot_product(centroid-xyz,centroid-xyz) .le. disthresq ) then
                            kntcoincid = 0
                            doii1: do ii = 0, 2
                                do jj = 0, 2
                                    if (indices(localaux(k,j)+ii) .eq. indices(localaux(l,i)+jj) &
                                        .and. (vtot(indices(localaux(k,j)+ii)) .ge. vmaxarray(i) &
                                        .or. vtot(indices(localaux(k,j)+ii)) .ge. vmaxarray(j)) ) then
                                            kntcoincid = kntcoincid + 1
                                            exit doii1
                                    endif
                                enddo
                            enddo doii1
                            if (kntcoincid .gt. 0) then
                                do i1 = 1, numlocalmax(j)
                                    localaux(numlocalmax(i)+i1,i) = localaux(i1,j)
                                enddo
                                numlocalmax(i) = numlocalmax(i) + numlocalmax(j)
                                numlocalmax(j) = 0
                                exit dok1
                            endif
                        endif
                    enddo
                enddo
            enddo dok1
        enddo
        knt = 0
        maxaux = 0
        do i = 1, kntgroupmax
            if (numlocalmax(i) .eq. 0) cycle
            maxaux = max(numlocalmax(i),maxaux)
            knt = knt + 1
        enddo
        deallocate(localmax)
        allocate (localmax(maxaux,knt), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ', i5,' when allocating localmax')") ierr
            abort = 1
            return
        endif

        localmax = 0
        knt = 0
        do i = 1, kntgroupmax
            if (numlocalmax(i) .eq. 0) cycle
            knt = knt + 1
            localmax(1:numlocalmax(knt),knt) = localaux(1:numlocalmax(i),i)
            numlocalmax(knt) = numlocalmax(i)
        enddo
        kntgroupmax = knt

        write(6,"(/'Number of final regions with local maxima = ', i3)") kntgroupmax
        do i = 1, kntgroupmax
            write(6,"('Number triangles in region ', i3 ' corrected = ', i7)") i, numlocalmax(i)
        enddo

        deallocate(localaux, kntlocal)
    endif
    
!         Scan regions of local minima to gather overlapping regions

    if (kntmin .gt. 0) then
        write(6,"(/111('-'),/'Number of starting regions with local minima = ', i3)") kntgroupmin
        do i = 1, kntgroupmin
            write(6,"('Number triangles in region ', i3 '  = ', i7)") i, numlocalmin(i)
        enddo
    else
        write(6,"(/'No regions with local minima')")
    endif

    if (kntmin .gt. 0) then
        allocate (localaux(kntmin,kntgroupmin), kntlocal(kntgroupmin), vminarray(kntgroupmin), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating localaux, kntlocal and vminarray. Stop')
        kntlocal = 0
        do i = 1, kntmin
            kntlocal(localmin(i,2)) = kntlocal(localmin(i,2))+1
            localaux(kntlocal(localmin(i,2)),localmin(i,2)) = localmin(i,1)
        enddo

        write(6,"(/,'Positions and mesp values of the local minima in the starting regions for local surface with mesp < ',&
            &e12.5)") thresmin
        vminarray = cero
        do jaux = 1, kntgroupmin
            vmin = 0.d0
            i1 = 0
            do i = 1, kntlocal(jaux)
                do k = 0, 2
                    if (vtot(indices(localaux(i,jaux)+k)) .lt. vmin) then
                        i1 = localaux(i,jaux)+k
                        vmin = vtot(indices(localaux(i,jaux)+k))
                    endif
                enddo
            enddo
            vminarray(jaux) = vmin
            if (i1 .lt. 1) cycle
            write(6,"('x = ', e12.5, ' y = ', e12.5, ' z = ', e12.5, ' mesp = ', e12.5, e12.5)") vertices(:,indices(i1)), vmin
        enddo

        do j = 2, kntgroupmin
            dok2: do k = 1, numlocalmin(j)
                centroid = ri(3) * ( vertices(:,indices(localaux(k,j))) + vertices(:,indices(localaux(k,j)+1)) &
                    + vertices(:,indices(localaux(k,j)+2)) )
                do i = 1, j-1
                    do l = 1, numlocalmin(i)
                        xyz = ri(3) * ( vertices(:,indices(localaux(l,i))) + vertices(:,indices(localaux(l,i)+1)) &
                            + vertices(:,indices(localaux(l,i)+2)) )
                        if ( dot_product(centroid-xyz,centroid-xyz) .le. disthresq ) then
                            kntcoincid = 0
                            doii2: do ii = 0, 2
                                do jj = 0, 2
                                    if (indices(localaux(k,j)+ii) .eq. indices(localaux(l,i)+jj) &
                                        .and. (vtot(indices(localaux(k,j)+ii)) .le. vminarray(i) &
                                        .or. vtot(indices(localaux(k,j)+ii)) .le. vminarray(j)) ) then
                                            kntcoincid = kntcoincid + 1
                                            exit doii2
                                    endif
                                enddo
                            enddo doii2
                            if (kntcoincid .gt. 0) then
                                do i1 = 1, numlocalmin(j)
                                    localaux(numlocalmin(i)+i1,i) = localaux(i1,j)
                                enddo
                                numlocalmin(i) = numlocalmin(i) + numlocalmin(j)
                                numlocalmin(j) = 0
                                exit dok2
                            endif
                        endif
                    enddo
                enddo
            enddo dok2
        enddo
        knt = 0
        minaux = 0
        do i = 1, kntgroupmin
            if (numlocalmin(i) .eq. 0) cycle
            minaux = max(numlocalmin(i),minaux)
            knt = knt + 1
        enddo
        deallocate(localmin)
        allocate (localmin(minaux,knt), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Error ', i5,' when allocating localmin')") ierr
            abort = 1
            return
        endif

        localmin = 0
        knt = 0
        do i = 1, kntgroupmin
            if (numlocalmin(i) .eq. 0) cycle
            knt = knt + 1
            numlocalmin(knt) = numlocalmin(i)
            localmin(1:numlocalmin(knt),knt) = localaux(1:numlocalmin(i),i)   
        enddo
        kntgroupmin = knt
        write(6,"(/'Number of final regions with local minima = ', i3)") kntgroupmin
        do i = 1, kntgroupmin
            write(6,"('Number triangles in region ', i3 ' corrected = ', i7)") i, numlocalmin(i)
        enddo
    endif
    return
    end
!
!**********************************************************************
!
  subroutine makehistogram
    USE DAMSGHOLE320_D, vold => vtot, vtot => vtot4, verticesold => vertices, vertices => vertices4
    implicit none    
    allocate (histogram(minhist:maxhist), stat = ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Memory error when allocating histogram. Stop')
    endif
    histogram = cero
    do i = 1, kntind, 3
        aux = max(vtot(indices(i)), vtot(indices(i+1)), vtot(indices(i+2)))             ! highest mesp in triangle vertices
        bux = min(vtot(indices(i)), vtot(indices(i+1)), vtot(indices(i+2)))             ! lowest mesp in triangle vertices
        cux = vtot(indices(i)) + vtot(indices(i+1)) + vtot(indices(i+2)) - aux - bux    ! intermediate mesp in triangle vertices
        xyztetr(:,mod(i-1,3)+1) = vertices(:,indices(i))
        xyztetr(:,mod(i-1,3)+2) = vertices(:,indices(i+1))
        xyztetr(:,mod(i-1,3)+3) = vertices(:,indices(i+2))
        a = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,1),xyztetr(:,2)-xyztetr(:,1)))
        b = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,3),xyztetr(:,2)-xyztetr(:,3)))
        c = sqrt(dot_product(xyztetr(:,1)-xyztetr(:,3),xyztetr(:,1)-xyztetr(:,3)))
        s = umed * (a+b+c)                          ! Triangle semiperimeter
        surftrian = sqrt(s*(s-a)*(s-b)*(s-c))     ! Heron's formula for the triangle surface
        histogram(int(aux*dlthistinv)) = histogram(int(aux*dlthistinv)) + ri(3) * surftrian
        histogram(int(bux*dlthistinv)) = histogram(int(bux*dlthistinv)) + ri(3) * surftrian
        histogram(int(cux*dlthistinv)) = histogram(int(cux*dlthistinv)) + ri(3) * surftrian
        if (aux .lt. thresmax .and. bux .gt. thresmin) cycle
        fpos = cero
        fneg = cero
        if (bux .ge. thresmax) then
            fpos = uno
        else if (aux .le. thresmin) then
            fneg = uno
        else 
            if( aux .gt. thresmax) then
                if (cux .gt. thresmax) then
                    fpos = dos * ri(3)
                else if (cux .lt. thresmax) then
                    fpos = uno * ri(3)
                else
                    fpos = umed
                endif
            endif
            if( bux .lt. thresmin) then
                if (cux .lt. thresmin) then
                    fneg = dos * ri(3)
                else if (cux .gt. thresmin) then
                    fneg = uno * ri(3)
                else
                    fneg = umed
                endif
            endif
        endif 
        surfpos = surfpos + fpos * surftrian
        surfneg = surfneg + fneg * surftrian
    enddo
    return
    end
!
!**********************************************************************
!
  subroutine makehistmax(jaux)
    USE DAMSGHOLE320_D, vold => vtot, vtot => vtot4, verticesold => vertices, vertices => vertices4
    implicit none  
    integer(KINT) :: jaux

    j = jaux
    do i = 1, numlocalmax(j)
        aux = max(vtot(indices(localmax(i,j))), vtot(indices(localmax(i,j)+1)), vtot(indices(localmax(i,j)+2)))
        bux = min(vtot(indices(localmax(i,j))), vtot(indices(localmax(i,j)+1)), vtot(indices(localmax(i,j)+2)))
        cux = vtot(indices(localmax(i,j))) + vtot(indices(localmax(i,j)+1)) + vtot(indices(localmax(i,j)+2)) - aux - bux
        xyztetr(:,1) = vertices(:,indices(localmax(i,j)))
        xyztetr(:,2) = vertices(:,indices(localmax(i,j)+1))
        xyztetr(:,3) = vertices(:,indices(localmax(i,j)+2))
        a = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,1),xyztetr(:,2)-xyztetr(:,1)))
        b = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,3),xyztetr(:,2)-xyztetr(:,3)))
        c = sqrt(dot_product(xyztetr(:,1)-xyztetr(:,3),xyztetr(:,1)-xyztetr(:,3)))
        s = umed * (a+b+c)                          ! Triangle semiperimeter
        surftrian = sqrt(s*(s-a)*(s-b)*(s-c))     ! Heron's formula for the triangle surface
        histpartial(int(aux*dlthistinv),j) = histpartial(int(aux*dlthistinv),j) + ri(3) * surftrian
        histpartial(int(bux*dlthistinv),j) = histpartial(int(bux*dlthistinv),j) + ri(3) * surftrian
        histpartial(int(cux*dlthistinv),j) = histpartial(int(cux*dlthistinv),j) + ri(3) * surftrian
        if (aux .lt. thresmax) cycle
        fpos = cero
        if (bux .ge. thresmax) then
            fpos = uno
        else 
            if (cux .gt. thresmax) then
                fpos = dos * ri(3)
            else if (cux .lt. thresmax) then
                fpos = uno * ri(3)
            else
                fpos = umed
            endif
        endif 
        surfpos = surfpos + fpos * surftrian
    enddo
    return
    end
!
!**********************************************************************
!
  subroutine makehistmin(jaux)
    USE DAMSGHOLE320_D, vold => vtot, vtot => vtot4, verticesold => vertices, vertices => vertices4
    implicit none    
    integer(KINT) :: jaux
    
    j = jaux
    do i = 1, numlocalmin(j)
        bux = min(vtot(indices(localmin(i,j))), vtot(indices(localmin(i,j)+1)), vtot(indices(localmin(i,j)+2)))
        aux = max(vtot(indices(localmin(i,j))), vtot(indices(localmin(i,j)+1)), vtot(indices(localmin(i,j)+2)))
        cux = vtot(indices(localmin(i,j))) + vtot(indices(localmin(i,j)+1)) + vtot(indices(localmin(i,j)+2)) - aux - bux
        xyztetr(:,1) = vertices(:,indices(localmin(i,j)))
        xyztetr(:,2) = vertices(:,indices(localmin(i,j)+1))
        xyztetr(:,3) = vertices(:,indices(localmin(i,j)+2))
        a = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,1),xyztetr(:,2)-xyztetr(:,1)))
        b = sqrt(dot_product(xyztetr(:,2)-xyztetr(:,3),xyztetr(:,2)-xyztetr(:,3)))
        c = sqrt(dot_product(xyztetr(:,1)-xyztetr(:,3),xyztetr(:,1)-xyztetr(:,3)))
        s = umed * (a+b+c)                          ! Triangle semiperimeter
        surftrian = sqrt(s*(s-a)*(s-b)*(s-c))     ! Heron's formula for the triangle surface
        histpartial(int(aux*dlthistinv),j+kntgroupmax) = histpartial(int(aux*dlthistinv),j+kntgroupmax) + ri(3) * surftrian
        histpartial(int(bux*dlthistinv),j+kntgroupmax) = histpartial(int(bux*dlthistinv),j+kntgroupmax) + ri(3) * surftrian
        histpartial(int(cux*dlthistinv),j+kntgroupmax) = histpartial(int(cux*dlthistinv),j+kntgroupmax) + ri(3) * surftrian
        if (bux .gt. thresmin) cycle
        fneg = cero
        if (aux .le. thresmin) then
            fneg = uno
        else 
            if (cux .lt. thresmin) then
                fneg = dos * ri(3)
            else if (cux .gt. thresmin) then
                fneg = uno * ri(3)
            else
                fneg = umed
            endif
        endif 
        surfneg = surfneg + fneg * surftrian
    enddo
    return
    end
!
!**********************************************************************
!
  subroutine color_generator(v, topcolor, red, green, blue)
    USE DAM320_D, only: KREAL
    implicit none
    real(KREAL) :: blue, green, red, topcolor, vmax, vmin, v 
    if (v .ge. topcolor) then
        red = 1.d0
        green = 0.d0
        blue = 0.d0
    else if(v .lt. -topcolor) then
        red = 0.d0
        green = 0.d0
        blue = 1.d0
    else if(v .ge. 0) then
        red = v / topcolor
        green = sqrt(topcolor*topcolor - v*v) / topcolor
        blue = 0
    else
        red = 0
        green = sqrt(topcolor*topcolor - v*v) / topcolor
        blue = -v / topcolor
    endif
    return
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
    USE DAMPOT320_D
    implicit none
    integer(KINT) :: i, ierr, k1, k12, l, l1, l1l1, l2, l2l2, lm
    integer(KINT) :: m, m1, m1a, m2, m2a, ma, mb, md, ms
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
    if (lexact) then
        allocate(ang((mxl+1)*(mxl+2)/2), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating ang. Stop')
!    ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
        ang(1) = umed / raizpi
        lm = 1
        do l = 1, mxl
            lm = lm + 1
            ang(lm) = ang(1) * sqrt(re(2*l+1))
            aux = ang(lm) * raiz2
            do m = 1, l
                lm = lm + 1
                aux = aux / sqrt(re(l-m+1)*re(l+m))
                ang(lm) = aux
            enddo
        enddo
        mxind = (mxldst+1)*(mxldst+2)/2
        allocate(ind(0:mxind), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating ind. Stop')
        ind(0) = 0
        do i = 1, mxind
            ind(i) = ind(i-1) + i         !  i*(i+1)/2
        enddo
        mxbin = max(mxldst,mxlenpol)
        allocate(bin((mxbin+1)*(mxbin+2)/2), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating bin. Stop')
        lm = 0
        do l = 0, mxbin
            do m = 0, l
                lm = lm + 1
                bin(lm) = fact(l) * facti(m) * facti(l-m)
            end do
        end do
!     Tabulates the coefficients for the decomposition of products
!     of two functions depending on phi (sin (m*phi), cos (m*phi))
!     into functions of the same type
        mxemes = mxldst
        allocate(ssv(-mxemes:mxemes,-mxemes:mxemes), sdv(-mxemes:mxemes,-mxemes:mxemes), &
            msv(-mxemes:mxemes,-mxemes:mxemes), mdv(-mxemes:mxemes,-mxemes:mxemes), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating msv, mdv, ssv and sdv. Stop')
        do m2 = -mxemes, mxemes
            do m1 = -mxemes, mxemes
                call emes ( m1, m2, ms, md, ss, sd )
                msv(m1,m2) = ms
                mdv(m1,m2) = md
                ssv(m1,m2) = ss
                sdv(m1,m2) = sd
            enddo
        enddo
!    Coefficients for the decomposition of products of regular spherical harmonics into
!    regular spherical harmonics
        mxlcof = mxldst*(mxldst+3)/2
        mxkcof = mxlcof*(mxlcof+3)/2
        allocate(app(0:2*mxl+1,0:mxkcof), bpp(0:2*mxl+1,0:mxkcof), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating app and bpp. Stop')
        if (longoutput) write(6,"('Size of app   = ', i15, ' bytes')") size(app)

        call acof
        call bcof
!    Tabulates some auxiliary indices for locating the previous coefficients
        allocate(indk12((mxl+1)**2,(mxl+1)**2), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating indk12. Stop')
        do l2 = 0,mxl
            do l1 = 0,mxl
                do m2 = -l2, l2
                    do m1 = -l1, l1
                        l1l1 = ind(l1)
                        l2l2 = ind(l2)
                        m1a = abs(m1)
                        m2a = abs(m2)
                        if ( l1.eq.l2 ) then
                            k1 = l1l1 + max(m1a,m2a)
                            k12 = ind(k1) + l1l1 + min(m1a,m2a)
                        elseif (l1.gt.l2) then
                            k1 = l1l1 + m1a
                            k12 = ind(k1) + l2l2 + m2a
                        else
                            k1 = l2l2 + m2a
                            k12 = ind(k1) + l1l1 + m1a
                        endif
                        indk12(l1*(l1+1)+m1+1,l2*(l2+1)+m2+1) = k12
                    end do
                end do
            end do
        end do
    endif
    return
    end
!
!    ***************************************************************
!
  subroutine readdamqtsghole
    USE MPI
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, nbasis, ncenbas, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, fr1, fr2l2, pi4d2l1, r, ra, ral1inv, rainv, ral
    real(KREAL) :: rinta, rintb, rlarex, step, stepmed, suml, suml1, suml2, summ, summ1, summ2, t
    real(KREAL) :: tcheb(0:mxlenpol-1)
    inquire(file=trim(projectname)//"_2016.damqt", size=nsize, iostat=ierr)
    if (ierr .ne. 0) call error(ierr,'Error when inquiring file '//trim(projectname)//"_2016.damqt")
    if (nsize .eq. -1) call error(1,'Size of file '//trim(projectname)//"_2016.damqt cannot be determined")
    if (longoutput) write(6,"('Size of file ', a, ' = ', i12)") trim(projectname)//"_2016.damqt", nsize
#if _WIN32
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when opening file '//trim(projectname)//"_2016.damqt"
        abort = 1
        return
    endif 
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when opening file '//trim(projectname)//"_2016.dmqtv"
        abort = 1
        return
    endif
#elif __INTEL_COMPILER
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when opening file '//trim(projectname)//"_2016.damqt"
        abort = 1
        return
    endif
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when opening file '//trim(projectname)//"_2016.dmqtv"
        abort = 1
        return
    endif
#else
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='unformatted', action = 'read', access='stream', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when opening file '//trim(projectname)//"_2016.damqt"
        abort = 1
        return
    endif
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='unformatted', action = 'read', access='stream', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when opening file '//trim(projectname)//"_2016.dmqtv"
        abort = 1
        return
    endif
#endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Opens files ', a, ' and ', a)") trim(projectname)//"_2016.damqt", &      
        trim(projectname)//"_2016.dmqtv"
    read(10) ncen, nbas, ncaps
    nsize = nsize - sizeof(ncen) - sizeof(nbas) - sizeof(ncaps)
    if (myrank .eq. 0) write(6,"('ncen = ', i4, ' nbas = ', i6, ' ncaps = ', i5)") ncen, nbas, ncaps

!    Allocates memory for geometry

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating atmnam'
        abort = 1
        return
    endif

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating nzn'
        abort = 1
        return
    endif

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating rcen'
        abort = 1
        return
    endif

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating zn'
        abort = 1
        return
    endif

    if (myrank .eq. 0) then
        write(6,"(/24x,'GEOMETRY (BOHR)')")
        write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge')")
    endif
    do ia = 1, ncen
        read(10) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia)
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
                nzn(ia) = 0
        else
                nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
        if (myrank .eq. 0) write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5)") &
                ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia)
    enddo
    nsize = nsize - sizeof(rcen) - sizeof(zn)

!    Basis set
    read(10) lsto    ! .true. means STO basis, .false. means GTO basis
    if (lsto) then
        if (lexact) then
            if (myrank .eq. 0) write(6,"('Exact potential not prepared for STO yet. Stop')")
            abort = 1
            return
        endif
!        Allocates memory for the basis set

        allocate(ll(ncaps), lmaxc(ncen), nf(ncaps), ngini(ncen), ngfin(ncen), nn(ncaps), rnor(ncaps), xx(ncaps), stat = ierr)
        if (ierr .ne. 0) then
            write(6,*) 'Memory error in proc ', myrank, ' when allocating ll, lmaxc, nf, ngini, ngfin, nn, rnor, xx'
            abort = 1
            return
        endif

        if (myrank .eq. 0) write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfin(ia)
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfin(ia))
            if (myrank .eq. 0 .and. longoutput) write(6,"(t5,'center ', i4,/t12,'n', t16, 'l', t25,'exp', t35, 'ind func')") ia
            lmaxc(ia) = 0
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                rnor(i) = sqrt((dos * xx(i))**(2*nn(i)+1) / fact(2*nn(i)))
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (myrank .eq. 0 .and. longoutput) write(6,"(t11,i2,t15,i2,t20,e12.5,t36,i4)") nn(i), ll(i), xx(i), nf(i)
            enddo
        enddo
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)

!        Allocates memory for the basis set

        allocate(cfcontr(nprimitot), ipntprim(ncaps), ll(ncaps), lmaxc(ncen), ncontr(ncen), nf(ncaps), ngini(ncen), &
            ngfin(ncen), nprimit(ncaps), rnor(ncaps), xxg(nprimitot), stat = ierr)
        if (ierr .ne. 0) then
            write(6,*) 'Memory error in proc ', myrank, ' when allocating cfcontr, ipntprim, ll, lmaxc, ncontr, &
                &nf, ngini, ngfin, nprimit, rnor, xxg'
            abort = 1
            return
        endif

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
        if (myrank .eq. 0) write(6,"(/t22,'GTO Basis set',/t22,13('-'))")
        if (myrank .eq. 0 .and. longoutput) then
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
        if (myrank .eq. 0) write(6,"('Number of basis functions = ', i4)") nbas
    endif

!    Data of density representation
    read(10) lmaxexp
    nsize = nsize - sizeof(lmaxexp)
    if (lmaxrep .gt. lmaxexp) then
        if (myrank .eq. 0) then
            write(6,"('lmaxrep = ', i3, ' greater than lmaxexp ', i3)") lmaxrep, lmaxexp
            write(6,"('takes lmaxrep = ',i3)") lmaxexp
        endif
        lmaxrep = lmaxexp
    endif
    lmtop = (lmaxexp+1)*(lmaxexp+1)

    if (myrank .eq. 0 .and. longoutput) write(6,"('lmaxexp = ', i2, ' nintervaj = ', i2)") lmaxexp, nintervaj

    allocate(icfposd(lmtop*nintervaj+1,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating icfposd'
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    
    nsize = nsize - sizeof(icfposd(:,1)) * ncenbas

    allocate(xajustd(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating xajustd'
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('Estimated highest size of xajustd   = ', i15, ' bytes')") size(xajustd)
    
    nsize = nsize - sizeof(xajustd(:,1)) * ncenbas

    allocate(cfajust(nsize/8), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating cfajust'
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) then
        write(6,"('Size of cfajust   = ', i15, ' bytes')") size(cfajust)
        write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
    endif
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
        if (myrank .eq. 0 .and. longoutput) write(6,"('fitting exponents: ',/, 8(1x,e17.10))")  xajustd(1:nintervaj,ia)
!     fitting coeficients
        read(10) cfajust(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)
    enddo

!    Generates an auxiliary index array for determining the interval to which a given r belongs
    step = rinterv(1)
    fct = uno / step
    lenindintrv = int(rinterv(nintervaj) * fct + udec)

    allocate(indintrv(lenindintrv), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating indintrv'
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of indintrv   = ', i15, ' bytes')") size(indintrv)

    r = cero
    interv = 1
    do i = 1, lenindintrv-1
        r = r + step
        if (r .gt. (rinterv(interv))) interv = interv + 1
        indintrv(i) = interv
    enddo
    indintrv(lenindintrv) = interv

!    Reads auxiliary integrals from file .dmqtv

    allocate(cfrint1(icfposd(lmtop*nintervaj+1,ncen)-1), cfrint2l2(icfposd(lmtop*nintervaj+1,ncen)-1), &
        QGacum(nintervaj*lmtop,ncen), Qgpart(nintervaj*lmtop), qppart(nintervaj*lmtop), qpacum(nintervaj*lmtop,ncen), &
        rlargo(ncen), rmultip(lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating cfrint1, cfrint2l2, QGacum, Qgpart, qppart&
            &qpacum, rlargo, rmultip'
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of cfrint1   = ', i15, ' bytes')") size(cfrint1)
    cfrint1 = cero
    cfrint2l2 = cero
    QGacum = cero
    Qgpart = cero
    qppart = cero
    qpacum = cero
    rlargo = cero
    rmultip = cero
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) cycle
!    multipolar moments
        read(11) rmultip(1:lmtop,ia)
!    Reads the integrals:
!        Qg(la,ma,i;ia) = 
!            Integrate[ r**(2*la+2) * fradtr[la,ma,r], {r,l_(i-1),l_i}]
!        qp(la,ma,i;ia) = 
!            Integrate[ r * fradtr[la,ma,r], {r,l_(i-1),l_i}]
!    and computes from them and stores the integrals:
!        QGacum(la,ma,i;ia) = 
!            Integrate[ r**(2*la+2) * fradtr[la,ma,r], {r,0,l_i}]
!        qpacum(la,ma,i;ia) = 
!            Integrate[ r * fradtr[la,ma,r], {r,l_(i-1),Infinity}]
!    where  ia  labels the center (atom),
        read(11) Qgpart(1:nintervaj*lmtop)
        read(11) qppart(1:nintervaj*lmtop)
        QGacum(1:lmtop,ia) = Qgpart(1:lmtop)
        do i = lmtop+1, nintervaj*lmtop, lmtop
            QGacum(i:i+lmtop-1,ia) = QGacum(i-lmtop:i-1,ia) + Qgpart(i:i+lmtop-1)
        enddo
        qpacum((nintervaj-1)*lmtop+1:nintervaj*lmtop,ia) = cero
        do i = (nintervaj-2)*lmtop+1, 1, -lmtop
            qpacum(i:i+lmtop-1,ia) = qpacum(i+lmtop:i+2*lmtop-1,ia) + qppart(i+lmtop:i+2*lmtop-1)
        enddo
!    Reads the fitting coeficients of auxiliary integrals for electrostatic potential and field
        read(11) cfrint1(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)    ! Expansion coefficients of auxiliary integrals rint1
        read(11) cfrint2l2(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)    ! Expansion coefficients of auxiliary integrals rint2l2
    enddo
    close(10)
    close(11)

!    Determines the long-range radii and the highest l in the expansion for each interval
    allocate(lcorto(nintervaj,ncen), llargo(0:mxlargo,ncen), Qllargo(0:lmaxrep), stat = ierr )
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating lcorto, llargo and Qllargo'
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) then
        write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)
        write(6,"('Size of llargo   = ', i15, ' bytes')") size(llargo)
        write(6,"('Size of Qllargo   = ', i15, ' bytes')") size(Qllargo)
    endif
    
!    long-range radii

    allocate(umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) then
        write(6,*) 'Memory error in proc ', myrank, ' when allocating umedpow'
        abort = 1
        return
    endif

    umedpow(0) = uno                            !
    do i = 1, lmaxexp                            !
        umedpow(i) = umedpow(i-1) * umed            ! 1 / 2^i
    enddo
    if (myrank .eq. 0 .and. .not. lexact) then
        write(6,"('Long-range threshold = ', e12.5)") umbrlargo
    endif
    nsamples = 4
    do ia = 1, ncen
        rlargo(ia) = cero
        if (ngini(ia) .le. 0) cycle
        kntlm = 0
        do l = 0, lmaxrep
            summ1 = cero
            do m = -l, l
                kntlm = kntlm + 1
                summ1 = summ1 + abs(rmultip(kntlm,ia)) * fact(l+abs(m)) * umedpow(abs(m)) &
                                        * facti(l-abs(m)) * facti(abs(m))
            enddo
            Qllargo(l) = summ1
        enddo
        rlargo(ia) = rinterv(nintervaj)
        lcorto(1:nintervaj,ia) = 0
        do interv = 1, nintervaj
            dltsample = udec * (rinterv(interv) - rinterv(interv-1))
            do i = 0, nsamples-1    ! samples over nsamples points in each interval to decide the highest l
                ra = rinterv(interv-1) + dltsample + (rinterv(interv) - rinterv(interv-1) - dos * dltsample) &
                        * ri(nsamples-1) * re(i)
                t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
                dost = t + t
                tcheb(0) = uno    ! Chebyshev T  polynomials
                tcheb(1) = t
                do j = 2, mxlenpol-1
                        tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
                enddo
                lm = 0
                rainv = uno / ra
                ral = uno
                ral1inv = rainv
                suml1 = cero
                suml2 = cero
                do l = 0, lmaxrep
                    pi4d2l1 = cuatro * pi * dosl1i(l)
                    summ2 = cero
                    do m = -l, l
                        lm = lm + 1
                        if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo) cycle
                        icflm = icfposd((interv-1)*lmtop+lm,ia)
                        rinta = cero
                        rintb = cero
                        do j = 0, icfposd((interv-1)*lmtop+lm+1,ia)-icflm-1
                                rinta = rinta + cfrint2l2(icflm+j) * tcheb(j)
                                rintb = rintb + cfrint1(icflm+j) * tcheb(j)
                        enddo
                        rinta = rinta * (ra-rinterv(interv-1))
                        if (interv .gt. 1) rinta = QGacum((interv-2)*lmtop+lm,ia) + rinta
                        rintb = qpacum((interv-1)*lmtop+lm,ia) + rintb * (rinterv(interv)-ra)
                        summ2 = summ2 + abs(pi4d2l1 * ( ral1inv * rinta + ral * rintb)) * fact(l+abs(m)) &
                                * umedpow(abs(m)) * facti(l-abs(m)) * facti(abs(m))
                    enddo
                    suml1 = suml1 + Qllargo(l) * ral1inv
                    suml2 = suml2 + summ2
                    if (abs(summ2) .gt. umbrlargo .and. l .gt. lcorto(interv,ia)) lcorto(interv,ia) = l
                    if (abs(suml2-suml1) .gt. umbrlargo) rlargo(ia) = rinterv(interv)
                    ral = ral * ra
                    ral1inv = ral1inv * rainv
                enddo
            enddo
        enddo
        if (myrank .eq. 0) then 
            if (longoutput) then
                write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5, ' lcorto = ', 30(i3))") &
                        ia, atmnms(nzn(ia)), rlargo(ia), lcorto(1:nintervaj,ia)
            else
                write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5)") &
                        ia, atmnms(nzn(ia)), rlargo(ia)
            endif
        endif
        llargo(0,ia) = lcorto(1,ia)
        do i = 1, mxlargo
            ra = re(i)
            rainv = uno / ra
            ral1inv = rainv
            kntlm = 0
            do l = 0, lmaxrep
                summ1 = cero
                do m = -l, l
                    kntlm = kntlm + 1
                    summ1 = summ1 + abs(rmultip(kntlm,ia)) * fact(l+abs(m)) * umedpow(abs(m)) &
                                            * facti(l-abs(m)) * facti(abs(m))
                enddo
                Qllargo(l) = summ1 * ral1inv
                ral1inv = ral1inv * rainv
            enddo
            llargo(i,ia) = 0
            suml1 = cero
            do l = lmaxrep, 0, -1
                suml1 = suml1 + Qllargo(l)
                if (suml1 .gt. umbrlargo) then
                    llargo(i,ia) = l
                    exit
                endif
            enddo
        enddo
        if (myrank .eq. 0 .and. longoutput) write(6,"('llargo: ', 51(i3))") llargo(0:mxlargo,ia)
    enddo
    deallocate (Qgpart, qppart)
    if (lexact .and. .not. lsto) then

!		Allocates the array containing the density matrix

        allocate(dmat(nbas,nbas), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating dmat in processor ',i3)") myrank
            abort = 1
            return
        endif
        if (myrank .eq. 0 .and. longoutput) write(6,"('Estimated highest size of dmat   = ', i15, ' bytes')") size(dmat)

        lgradient = .false.	! The gradient is not computed for the exact density
        open(16,file=trim(projectname)//".den",form='formatted', iostat=ierr)
        if (ierr .ne. 0) then
            write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//".den", myrank
            abort = 1
        else
            read(16,*, iostat = ierr) nbasis, ((dmat(i,j), j=1,i), i=1,nbasis)
            do j = 1, nbasis
                do i = 1, j-1
                    dmat(i,j) = dmat(j,i)
                enddo
            enddo
        endif
        if ( ierr .ne. 0 .or. nbas .ne. nbasis ) then
            write(6,"('ERROR reading density matrix in processor ',i3)") myrank
            write(6,"('Check whether the density matrix correspond to this basis set.')")
            abort = 1
            return
        endif
        close(16)
    endif
    deallocate(umedpow)
    return
    end
        
!   ***************************************************************

  subroutine densrepr(ia, x, y, z, denrep, dendrvx, dendrvy, dendrvz)
    USE DAM320_D
    USE DAMPOT320_D
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

!
!   **************************************************************************************
!     Calculates the electrostatic potential from the density expansion at point (x,y,z)
!
   subroutine mesp(ia, x, y, z, vnucl, vel, vtot)
    USE DAM320_D
    USE DAMSGHOLE320_D, only: kntvert, vertices
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    implicit none
    integer(KINT) :: i, ia, ierr, interv, icflm, j, jshft, kntlm, l, lm, lmtopot, ltop, m
    real(KREAL) :: aux, bux, dost, drvx, drvy, drvz
    real(KREAL) :: flm, fux
    real(KREAL) :: pi4exp, ra, ra2, rainv, ra2inv, rinta, rintb, sgn, t, tp, vnucl, vel, vlm, vqlm, vtot
    real(KREAL) :: x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL), parameter :: vtope = 1.d10        ! To prevent infinity, if the point coincides with a nucleus, loads the value of parameter vtope
    real(KREAL) :: tcheb(0:mxlenpol-1), pi4d2l1((lmaxexp+1)*(lmaxexp+1)), d2l1((lmaxexp+1)*(lmaxexp+1))
    real(KREAL), allocatable :: zn(:)
    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating zn in multmolec. Stop')
    if (lvalence) then
        do i = 1, ncen
            zn(i) = atom_core(zn_orig(i))
        enddo
    else
        zn = zn_orig
    endif
    pi4exp = cuatro * pi
    xa = x - rcen(1,ia)
    ya = y - rcen(2,ia)
    za = z - rcen(3,ia)
    ra2 = xa*xa+ya*ya+za*za
    vnucl = cero
    vel = cero
    vtot = cero
    if (ra2 .lt. geomthr*geomthr) then
        vnucl = vtope
        vtot = vtope
        if (ngini(ia) .le. 0) return
        icflm = icfposd(1,ia)
        sgn = uno
        rintb = cero
        do i = 0, icfposd(2,ia)-icflm-1
            rintb = rintb + cfrint1(icflm+i) * sgn
            sgn = -sgn
        enddo
        vel = - cuatro * pi * (qpacum(1,ia) + rintb * rinterv(1))
        return
    endif
    ra = sqrt(ra2)
    rainv = uno / ra
    ra2inv = rainv * rainv
    if (ngini(ia) .le. 0) then
        vnucl = zn(ia) * rainv
        vtot = vnucl
        return
    endif
    if (ra .lt. rlargo(ia)) then
        interv = indintrv(int(fct*ra)+1)
        ltop = lcorto(interv,ia)
    else
        ltop = llargo(min(int(ra),mxlargo),ia)
    endif
    lmtopot = (ltop+1)*(ltop+1)
    vnucl = zn(ia) * rainv
    xadivra = xa * rainv
    yadivra = ya * rainv
    zadivra = za * rainv
    zlma(1) = uno        ! Regular spherical harmonics of r-R(ia)
    zlma(2) = ya
    zlma(3) = za
    zlma(4) = xa
    ra2l1(1) = ra
    ra2l1inv(1) = rainv
    pi4d2l1(1) = cuatro * pi
    d2l1(1) = uno
    lm = 1
    do l = 1, ltop
        zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))        ! zlm(l+1,l+1,ia)
        zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))        ! zlm(l+1,-(l+1),ia)
        zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)                ! zlm(l+1,l,ia)
        zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)                    ! zlm(l+1,-l,ia)
        do m = 0, l-1
            zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))    ! zlm(l+1,m,ia)
            zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
        enddo
        aux = ra2l1(lm) * ra2
        bux = ra2l1inv(lm) * ra2inv
        do m = -l, l
            lm = lm + 1
            ra2l1(lm) = aux
            ra2l1inv(lm) = bux
            pi4d2l1(lm) = cuatro * pi * dosl1i(l)    !   4 * pi / (2l+1)
            d2l1(lm) = dosl1(l)
        enddo
    enddo
!   Derivatives of the potential are not computed
    if (ra .ge. rlargo(ia)) then  ! The point is in the long-range region (ra >= rlargo(ia) )
        kntlargo = kntlargo + 1
        vel = - dot_product( rmultip(1:lmtopot,ia),zlma(1:lmtopot)*ra2l1inv(1:lmtopot) )
    else        !     The point is in the short-range region (ra < rlargo(ia) )
        kntcorto = kntcorto + 1
        t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
        dost = t + t
        tcheb(0) = uno    ! Chebyshev T  polynomials
        tcheb(1) = t
        do j = 2, mxlenpol-1
            tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
        enddo
        do lm = 1, lmtopot
            if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo) cycle
            icflm = icfposd((interv-1)*lmtop+lm,ia)
            rinta = cero
            rintb = cero
            do i = 0, icfposd((interv-1)*lmtop+lm+1,ia)-icflm-1
                rinta = rinta + cfrint2l2(icflm+i) * tcheb(i)
                rintb = rintb + cfrint1(icflm+i) * tcheb(i)
            enddo
            rinta = rinta * (ra-rinterv(interv-1))
            if (interv .gt. 1) rinta = QGacum((interv-2)*lmtop+lm,ia) + rinta
            rintb = qpacum((interv-1)*lmtop+lm,ia) + rintb * (rinterv(interv)-ra)
            vel = vel - pi4d2l1(lm) * ( rinta + ra2l1(lm) * rintb ) * zlma(lm) * ra2l1inv(lm)
        enddo
    endif  ! End of test over long/short-range
    vtot = vnucl + vel
    return
    end
    
!
!   **************************************************************************************
!     Calculates the electrostatic potential from the density expansion at point (x,y,z)
!
  subroutine statMESP
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAMSGHOLE320_D, fpos_ori => fpos, fneg_ori => fneg
    implicit none
    real(KREAL) :: Ac(3), Bc(3), Cc(3)
    real(KREAL) :: fA, fB, fC, area, apos, aneg, fdev, fneg, fnegsq, fnegtot, fnegsqtot, fpos, fpossq, fpostot, fpossqtot 

    areatot = cero
    apostot = cero
    anegtot = cero
    fpostot = cero
    fnegtot = cero
    fpossqtot = cero
    fnegsqtot = cero
    do i = 1, kntind, 3
        if (indices(i) .eq. indices(i+1) .or. indices(i) .eq. indices(i+2) .or. indices(i+1) .eq. indices(i+2) ) cycle
        Ac(:) = vertices4(:,indices(i))
        Bc(:) = vertices4(:,indices(i+1))
        Cc(:) = vertices4(:,indices(i+2))
        fA = vtot4(indices(i)) 
        fB = vtot4(indices(i+1)) 
        fC = vtot4(indices(i+2))
        call sub_statMESP(Ac, Bc, Cc, fA, fB, fC, area, apos, aneg, fpos, fneg, fpossq, fnegsq)
        areatot = areatot + area
        apostot = apostot + apos
        anegtot = anegtot + aneg
        fpostot = fpostot + fpos
        fnegtot = fnegtot + fneg
        fpossqtot = fpossqtot + fpossq
        fnegsqtot = fnegsqtot + fnegsq
    enddo
    fmed = (fpostot + fnegtot) / areatot
    if (apostot .gt. cero) then
        fmedpos = fpostot / apostot
        fvarpos = fpossqtot / apostot - fmedpos**2
    else 
        fmedpos = cero
        fvarpos = cero
    endif
    if (anegtot .gt. cero) then
        fmedneg = fnegtot / anegtot
        fvarneg = fnegsqtot / anegtot - fmedneg**2
    else 
        fmedneg = cero
        fvarneg = cero
    endif
    fvartot = fvarpos + fvarneg
    
    fdevtot = cero
    do i = 1, kntind, 3
        if (indices(i) .eq. indices(i+1) .or. indices(i) .eq. indices(i+2) .or. indices(i+1) .eq. indices(i+2) ) cycle
        Ac(:) = vertices4(:,indices(i))
        Bc(:) = vertices4(:,indices(i+1))
        Cc(:) = vertices4(:,indices(i+2))
        fA = vtot4(indices(i)) 
        fB = vtot4(indices(i+1)) 
        fC = vtot4(indices(i+2))
        call sub_MESP_dev(Ac, Bc, Cc, fA, fB, fC, fmed, fdev)
        fdevtot = fdevtot + fdev
    enddo
    fdevtot = fdevtot / areatot

    return
    end
  
!   ------------------------------------------------------------------------------------------
 
!  subroutine sub_statMESP: 
!       receives:   three points, Ac, Bc, Cc, defining the vertices of a triangle and the values of a function in the points, fA, fB, fC
!       returns:    the total area, and the areas for positive values of f (apos) and negative values of f (aneg),
!                   the integrals of positive f (fpos) and negative f (fneg) in their respective surfaces
!                   the integrals of the squares of positive f (fpossq) and negative f (fnegsq) in their respective surfaces
!
!       Computations are made with the assumption that f is linearly fitted to the three points (a plane)
!  
  subroutine sub_statMESP(Ac, Bc, Cc, fA, fB, fC, area, apos, aneg, fpos, fneg, fpossq, fnegsq)
    implicit none
    logical :: leqsgn
    integer(4) :: i, j
    real(8) :: Ac(3), Bc(3), Cc(3), Paux(3)
    real(8) :: a, aAMNB, aMNC, aneg, apos, area, aux, b, c
    real(8) :: fA, fAMNB, fB, fC, fMNC, fneg, fnegsq, fpos, fpossq, fsqAMNB, fsqMNC 
    real(8) :: u2, u3, v3
    
    leqsgn = .true.
    if (fA*fB .lt. 0.d0) then
        if (fB*fC .lt. 0.d0) then
            aux = fB ; fB = fC ; fC = aux
            Paux = Bc ; Bc = Cc ; Cc = Paux
            leqsgn = .false.
        else
            aux = fA ; fA = fC ; fC = aux
            Paux = Ac ; Ac = Cc ; Cc = Paux
            leqsgn = .false.
        endif
    else if (fB*fC .lt. 0.d0) then
        leqsgn = .false.
    endif
    
    u2 = sqrt(dot_product(Bc-Ac,Bc-Ac))
    u3 = dot_product(Cc-Ac,Bc-Ac) / u2
    v3 = sqrt(dot_product(Cc-Ac-u3*(Bc-Ac)/u2,Cc-Ac-u3*(Bc-Ac)/u2))
    
    area = u2 * v3 * 0.5d0
    if (leqsgn) then
        aAMNB = 0.d0
        fAMNB = 0.d0
        fsqAMNB = 0.d0
        aMNC = 0.5d0 * u2 * v3
        fMNC = area * (fA + fB + fC) / 3.d0
        fsqMNC = (fA**2 + fB**2 + fC**2 + fA*fB + fA*fC + fB*fC) * u2 * v3 / 12.d0
    else
        aMNC = fC * fC * u2 * v3 / (2.d0 *(fA-fC) * (fB-fC))
        fMNC = aMNC * (fA-fC) * (fB-fC) * fC / (3.d0 * (fC-fA) * (fC-fB) )
        fsqMNC = fC**4 * u2 * v3 / (12.d0 * (fC-fA) * (fC-fB) )
        aAMNB = (fA*fC + fB*fC - fA*fB) * u2 * v3 / (2.d0 * (fA-fC) * (fC-fB))
        fAMNB = aAMNB * ((fA+fB)**2 * fC - fA * fB * (fA+fB+fC)) / (3.d0 * (fA*fC + fB*fC - fA*fB) )
        fsqAMNB = (-fA*fB * (fA**2 + fA*fB + fB**2) + (fA+fB) * (fA**2+fB**2) * fC) * u2 * v3 &
            / (12.d0 * (fA-fC) * (fC-fB))
    endif
    
    if (fC .lt. 0) then
        apos = aAMNB
        aneg = aMNC
        fpos = fAMNB
        fneg = fMNC
        fpossq = fsqAMNB
        fnegsq = fsqMNC
    else
        aneg = aAMNB
        apos = aMNC
        fneg = fAMNB
        fpos = fMNC
        fnegsq = fsqAMNB
        fpossq = fsqMNC
    endif
    return
    end    
      
!   ------------------------------------------------------------------------------------------
 
!  subroutine sub_MESP_dev: 
!       receives:   three points, Ac, Bc, Cc, defining the vertices of a triangle and the values of a function in the points, fA, fB, fC
!       returns:    |f - fmed|
!
!       Computations are made with the assumption that f is linearly fitted to the three points (a plane)
!  
  subroutine sub_MESP_dev(Ac, Bc, Cc, fA, fB, fC, fmed, fdev)
    implicit none
    logical :: leqsgn
    integer(4) :: i, j
    real(8) :: Ac(3), Bc(3), Cc(3), Paux(3)
    real(8) :: a, aAMNB, aMNC, aneg, apos, area, aux, b, c
    real(8) :: fA, fAMNB, fB, fC, fdev, fmed, fMNC, fsqAMNB, fsqMNC 
    real(8) :: u2, u3, v3
    
    leqsgn = .true.
    if (fA*fB .lt. 0.d0) then
        if (fB*fC .lt. 0.d0) then
            aux = fB ; fB = fC ; fC = aux
            Paux = Bc ; Bc = Cc ; Cc = Paux
            leqsgn = .false.
        else
            aux = fA ; fA = fC ; fC = aux
            Paux = Ac ; Ac = Cc ; Cc = Paux
            leqsgn = .false.
        endif
    else if (fB*fC .lt. 0.d0) then
        leqsgn = .false.
    endif
    
    u2 = sqrt(dot_product(Bc-Ac,Bc-Ac))
    u3 = dot_product(Cc-Ac,Bc-Ac) / u2
    v3 = sqrt(dot_product(Cc-Ac-u3*(Bc-Ac)/u2,Cc-Ac-u3*(Bc-Ac)/u2))
    
    area = u2 * v3 * 0.5d0
    if (leqsgn) then
        aAMNB = 0.d0
        fAMNB = 0.d0
        aMNC = 0.5d0 * u2 * v3
        fMNC = area * (fA + fB + fC) / 3.d0
    else
        aMNC = fC * fC * u2 * v3 / (2.d0 *(fA-fC) * (fB-fC))
        fMNC = aMNC * (fA-fC) * (fB-fC) * fC / (3.d0 * (fC-fA) * (fC-fB) )
        aAMNB = (fA*fC + fB*fC - fA*fB) * u2 * v3 / (2.d0 * (fA-fC) * (fC-fB))
        fAMNB = aAMNB * ((fA+fB)**2 * fC - fA * fB * (fA+fB+fC)) / (3.d0 * (fA*fC + fB*fC - fA*fB) )
    endif
    
    fdev = abs(fAMNB - fmed * aAMNB) + abs(fMNC - fmed * aMNC)
    
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
!   ***************************************************************
!     Calculates the electrostatic potential from the exact density at point (x,y,z)
!
   subroutine GTOexactmesp(x, y, z, vnucl, vel, vtot)
    USE MPI
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, i1, i1p, i2, i2p, ia, ib, ierr, indf1, k, k1, k12, k2, kmax, km1, km2
    integer(KINT) :: l, l1, l12, l2, l1l1, l2l2, lk, lm, lmax, ltop, m, m1, m1a, m2, m2a, ms, msa, md, mda
    real(KREAL) :: angnorm, aux, bux, c1, c2, cosal, cosbet, cosga, cux, expbux, funcF0, potaux, potia, potiab
    real(KREAL) :: rabinv, rra, rra2, rrab, rrab2, rrai, rrp, rrp2, rrpi, rn, rn1
    real(KREAL) :: umur, ur, saux, sinal, sinbet, singa, suma, ss, sd
    real(KREAL) :: vel, vnucl, vtot, x, xinv, xp, xx0, xxa, xxab, xxp, xy
    real(KREAL) :: y, yy0, yya, yyab, yyp, z, zz, zz0, zza, zzab, zzp
    real(KREAL) :: fv(0:20), gammaG(mxl), Qg(0:mxl), qp(0:mxl), rpw(0:4*mxl+2), urpw(0:mxl), umurpw(0:mxl)
    real(KREAL) :: vaux((mxl+1)**2,(mxl+1)**2)
    real(KREAL) :: roaux(-mxl:mxl,-mxl:mxl)
    real(KREAL), parameter :: vtope = 1.d10		! To prevent infinity, if the point coincides with a nucleus, loads the value of parameter vtope
    integer(KINT), parameter :: mxlcofpot = mxl*(mxl+3)/2, mxkcofpot = mxlcofpot*(mxlcofpot+3)/2
    real(KREAL), allocatable :: zn(:)
    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating zn in multmolec. Stop')
    if (lvalence) then
        do i = 1, ncen
            zn(i) = atom_core(zn_orig(i))
        enddo
    else
        zn = zn_orig
    endif

    vel = cero
    vnucl = cero
    do ia = 1, ncen    ! Loop over the first center
        ltop = 2*lmaxc(ia)
        xxa = x - rcen(1,ia)
        yya = y - rcen(2,ia)
        zza = z - rcen(3,ia)
        rra2 = xxa*xxa+yya*yya+zza*zza
        if (rra2 .lt. geomthr*geomthr) then
                vnucl = vtope
        endif
        call armonicos(mxldst, xxa, yya, zza, zlma)	! Regular spherical harmonics
        rra = sqrt(rra2)
!       rpw(i) = 1 / rrp**i
        if (rra .lt. 1.e-20) then
            rpw = uno               ! To prevent indeterminations below
        else
            rpw(0) = uno
            rrai = uno / rra
            do k = 1, 4*lmaxc(ia)+2
                rpw(k) = rpw(k-1) * rrai
            end do
        endif
        potia = cero
        do i1 = ngini(ia), ngfin(ia)
            l1 =  ll(i1)
            do i2 = ngini(ia), ngfin(ia)
                l2 =  ll(i2)
                kmax = (l1+l2)/2
                do k = 0, kmax
                    Qg(k) = cero
                    qp(k) = cero
                end do
                rn = rnor(i2) * rnor(i1)
                potaux = cero
                do i1p = ipntprim(i1), ipntprim(i1)+nprimit(i1)-1
                    do i2p = ipntprim(i2), ipntprim(i2)+nprimit(i2)-1
                        bux = (xxg(i1p)+xxg(i2p))*rra2
                        xinv = uno / (xxg(i1p)+xxg(i2p))
                        expbux = exp(-bux)
                        gammaG(1) = expbux * xinv	! gammaG(i) = Gamma[i,xcc*r] / xcc**i
                        cux = expbux * rra2
                        do i = 1, kmax
                            gammaG(i+1) = xinv * (re(i) * gammaG(i) + cux)
                            cux = cux * rra2
                        end do
!                           fv(i) = rra**(2i) * Gamma[i+1/2,0,(xx(i1p)+xx(i2p))*rra] / (xx(i1p)+xx(i2p))**i
                        aux = rra*rra
                        if (rra .lt. 1.e-20) then
                            fv = cero
                        else
                            call freq20(aux,bux,fv)
                        endif
                        do k = 0, kmax
                            Qg(k) = Qg(k) + cfcontr(i1p) * cfcontr(i2p) * fv(l1+l2-k+1)
                            qp(k) = qp(k) + cfcontr(i1p) * cfcontr(i2p) * gammaG(k+1)
                        end do
                    end do
                end do
                do m1 = -l1, l1
                    do m2 = -l2, l2
                        ms = msv(m1,m2)
                        md = mdv(m1,m2)
                        ss = ssv(m1,m2)
                        sd = sdv(m1,m2)
                        msa = abs(ms)
                        mda = abs(md)
                        angnorm = ang(ind(l1)+abs(m1)+1) * ang(ind(l2)+abs(m2)+1)
                        k12 = indk12(l1*(l1+1)+m1+1,l2*(l2+1)+m2+1)
                        if (abs(ss) .ge. 1.d-2) then
                            do k = 0, (l1+l2-msa)/2
                                lk = l1+l2-2*k
                                saux = ss*app(lk,k12) * zlma(lk*(lk+1)+ms+1) * angnorm * ri(2*lk+1) &
                                                * (Qg(k) * rpw(2*lk) + qp(k) ) * dmat(nf(i1)+l1+m1,nf(i2)+l2+m2)
                                potaux = potaux + saux
                            end do
                        end if
                        if (abs(sd) .ge. 1.d-2) then
                            do k = 0, (l1+l2-mda)/2
                                lk = l1+l2-2*k
                                saux = sd*bpp(lk,k12) * zlma(lk*(lk+1)+md+1) * angnorm * ri(2*lk+1) &
                                        * (Qg(k) * rpw(2*lk) + qp(k) ) * dmat(nf(i1)+l1+m1,nf(i2)+l2+m2)
                                potaux = potaux + saux
                            end do
                        end if
                    end do
                end do
                potia = potia + rn * potaux
            end do
        end do
        vel = vel - dos * pi * potia
        if (vnucl .lt. vtope) vnucl = vnucl + zn(ia) / rra
        do ib = 1, ncen      ! Loop over the second center
            if (ia .eq. ib) cycle
            xxab = rcen(1,ib) - rcen(1,ia)
            yyab = rcen(2,ib) - rcen(2,ia)
            zzab = rcen(3,ib) - rcen(3,ia)
            rrab2 = xxab*xxab + yyab*yyab + zzab*zzab
            xy = dsqrt(xxab*xxab + yyab*yyab)
            rrab = dsqrt(rrab2)
            if (rrab .lt. 1.d-10) then
                    write(6,"('Error: centers ', i3, ' and ', i3, ' coincide.')") ia, ib
                    exit
            end if
            if (xy .gt. 1.d-10) then
                    sinal = yyab / xy
                    cosal = xxab / xy
            else
                    sinal = cero
                    cosal = uno
            end if
            rabinv = uno / rrab
            sinbet = xy * rabinv
            cosbet = zzab * rabinv
            singa = cero
            cosga = uno
            lmax = max(1,lmaxc(ia) + lmaxc(ib))
            call rotar (lmax, cosal, sinal, cosbet, sinbet, cosga, singa)
!			coordinates of point  ir  with respect to point A in the original (molecular) axis system
            xx0 = x - rcen(1,ia)
            yy0 = y - rcen(2,ia)
            zz0 = z - rcen(3,ia)
!			coordinates of point  ir  with respect to point A in the lined-up axis system
            yyp = rl(-1,-1,1) * yy0 + rl(0,-1,1) * zz0 + rl(1,-1,1) * xx0
            zz = rl(-1,0,1) * yy0 + rl(0,0,1) * zz0 + rl(1,0,1) * xx0
            xxp = rl(-1,1,1) * yy0 + rl(0,1,1) * zz0 + rl(1,1,1) * xx0
            potiab = cero
            do i1 = ngini(ia), ngfin(ia)
                l1 =  ll(i1)
                rn1 = rnor(i1)
                do i2 = ngini(ib), ngfin(ib)
                    l2 =  ll(i2)
                    l12 = l1 + l2
                    kmax = l12 / 2
                    rn = rnor(i2) * rn1
! 						Reads the block of density matrix and rotates it to the AB alligned system. Loads the result in matrix  roblk.
! 						Angular normalization factors are introduced at the end of loading process.
                    do m1 = -l1, l1
                        indf1 = nf(i1)+l1+m1
                        do m2 = -l2, l2
                                roblk(m1,m2) = dmat(indf1,nf(i2)+l2+m2)
                        end do
                    end do
                    do m1 = -l1, l1		! Rotation on center B
                        do m2 = -l2, l2
                            suma = cero
                            do k = -l2, l2
                                suma = suma + roblk(m1,k) * rl(k,m2,l2)
                            end do
                            roaux(m1,m2) = suma
                        end do
                    end do
                    do m1 = -l1, l1		! Rotation on center A and introduction of the angular normalization
                        do m2 = -l2, l2
                            suma = cero
                            do k = -l1, l1
                                suma = suma + roaux(k,m2) * rl(k,m1,l1)
                            end do
                            roblk(m1,m2) = suma * ang(ind(l1)+abs(m1)+1) * ang(ind(l2)+abs(m2)+1) * rn
                        end do
                    end do
!	Loops over the primitives associated to the contraction
                    do i1p = ipntprim(i1), ipntprim(i1)+nprimit(i1)-1
                        do i2p = ipntprim(i2), ipntprim(i2)+nprimit(i2)-1
                            xp = xxg(i1p) + xxg(i2p)
                            ur = rrab * xxg(i2p) / xp		! ur = u * rrab = xx(ip2) * rrab / xp
                            umur = ur - rrab			! umur = (u-1) * rrab = -xx(ip1) * rrab / xp
                            urpw(0) = uno			! urpw(i) = (u * rrab)**i
                            do i = 1, l1
                                urpw(i) = urpw(i-1) * ur
                            end do
                            umurpw(0) = uno			! umurpw(i) = ((u-1) * rrab)**i
                            do i = 1, l2
                                umurpw(i) = umurpw(i-1) * umur
                            end do
                            zzp = zz - ur	! Z coordinate of point  ir  with respect to point P in the lined-up axis system
                            rrp2 = xxp*xxp + yyp*yyp + zzp*zzp
                            rrp = sqrt(rrp2)
!                           rpw(i) = 1 / rrp**i
                            if (rrp .lt. 1.e-20) then
                                rpw = uno               ! To prevent indeterminations below
                            else
                                rpw(0) = uno
                                rrpi = uno / rrp
                                do k = 1, 2*l12+2
                                    rpw(k) = rpw(k-1) * rrpi
                                end do
                            endif
                            call armonicos(l12, xxp, yyp, zzp, zlma)	! Regular spherical harmonics
                            bux = xp*rrp2
                            xinv = uno / xp
                            expbux = exp(-bux)
                            gammaG(1) = expbux * xinv	! gammaG(i) = Gamma[i,xp*rrp] / xp**i
                            cux = expbux * rrp2
                            do i = 1, kmax
                                gammaG(i+1) = xinv * (re(i) * gammaG(i) + cux)
                                cux = cux * rrp2
                            end do
                            aux = rrp*rrp
!                           fv(i) = rra**(2i) * Gamma[i+1/2,0,(xx(i1p)+xx(i2p))*rra] / (xx(i1p)+xx(i2p))**i
                            if (rrp .lt. 1.e-20) then
                                fv = cero
                            else
                                call freq20(aux,bux,fv)
                            endif

                            do m2 = -l2, l2
                            do k2 = abs(m2), l2
                                    km2 = k2*(k2+1)+m2+1
                                    do m1 = -l1, l1
                                        do k1 = abs(m1), l1
                                            km1 = k1*(k1+1)+m1+1
                                            vaux(km1,km2) = cero
                                            ms = msv(m1,m2)
                                            md = mdv(m1,m2)
                                            msa = abs(ms)
                                            mda = abs(md)
                                            ss = ssv(m1,m2)
                                            sd = sdv(m1,m2)
                                            k12 = indk12(k1*(k1+1)+m1+1,k2*(k2+1)+m2+1)
                                            if (abs(ss) .ge. 1.d-2) then
                                                do k = 0, (k1+k2-msa)/2
                                                    lk = k1+k2-2*k
                                                    vaux(km1,km2) = vaux(km1,km2) + ss * app(lk,k12) &
                                                            * zlma(lk*(lk+1)+ms+1) * ri(2*lk+1) *(fv(k1+k2-k+1) &
                                                            * rpw(2*lk) + gammaG(k+1) )
                                                end do
                                            end if
                                            if (abs(sd) .ge. 1.d-2) then
                                                do k = 0, (k1+k2-mda)/2
                                                    lk = k1+k2-2*k
                                                    vaux(km1,km2) = vaux(km1,km2) + sd * bpp(lk,k12) &
                                                            * zlma(lk*(lk+1)+md+1) * ri(2*lk+1) * (fv(k1+k2-k+1) &
                                                            * rpw(2*lk) + gammaG(k+1) )
                                                end do
                                            end if
                                        end do	! End of Do over m1
                                    end do	! End of Do over k1
                                end do	! End of Do over m2
                            end do	! End of Do over k2
                            potaux = cero
                            do m1 = -l1, l1
                                m1a = abs(m1)
                                do k1 = m1a, l1
                                    c1 = bin(ind(l1+m1a)+k1+m1a+1) * urpw(l1-k1)
                                    km1 = k1*(k1+1)+m1+1
                                    do m2 = -l2, l2
                                        m2a = abs(m2)
                                        do k2 = m2a, l2
                                            c2 = bin(ind(l2+m2a)+k2+m2a+1) * umurpw(l2-k2)
                                            km2 = k2*(k2+1)+m2+1
                                            potaux = potaux + c1 * c2 * vaux(km1, km2) * roblk(m1,m2)
                                        end do
                                    end do
                                end do
                            end do
                            potiab = potiab + potaux * cfcontr(i1p) * cfcontr(i2p) &
                                    * exp(-xxg(i1p) * xxg(i2p) * rrab2 / xp)
                        end do		! End of Do over the primitives of the contraction on B
                    end do		! End of Do over the primitives of the contraction on A
                end do		! End of Do over the contractions on B
            end do		 ! End of Do over the contractions on A
            vel = vel - dos * pi * potiab
        end do  ! End of loop over second center
    end do  ! End of loop over first center
    if (vnucl .lt. vtope) then
        vtot = vnucl + vel
    else
        vtot = vtope
    endif
    return
    end

!   ***************************************************************

  subroutine armonicos(lmax, xxa, yya, zza, zlma)
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: l, lmax, m
    real(KREAL) :: rra2, xxa, yya, zza
    real(KREAL) :: zlma((2*mxl+1)**2)
!	Tabulation of regular spherical harmonics associated to the position vector of point (x,y,z) relative
!	to the positions of the different centers, (XA,YA,ZA):
!
!     	zlm(l,m,ia) = zlm(l,m,x-XA,y-YA,z-ZA)
!
!     	where  ia   numerates the centers (nuclei), the coordinates of center ia being  (XA,YA,ZA)
!
!     	the indices (l,m) are contracted into a single one:	lm = l*(l+1)+m      lm = 0, 1, 2, ... (lmaxexp+1)**2-1
    rra2 = xxa*xxa+yya*yya+zza*zza
    zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
    zlma(2) = yya
    zlma(3) = zza
    zlma(4) = xxa
    do l = 1, lmax-1
        zlma((l+1)*(l+3)+1)=re(l+l+1)*(xxa*zlma(l*(l+2)+1)-yya * zlma(l*l+1))		! zlm(l+1,l+1)
        zlma((l+1)*(l+1)+1)=re(l+l+1)*(yya*zlma(l*(l+2)+1)+xxa*zlma(l*l+1))		! zlm(l+1,-(l+1))
        zlma((l+2)*(l+2)-1)=re(l+l+1)*zza* zlma(l*(l+2)+1)				! zlm(l+1,l)
        zlma(l*(l+2)+3) = re(l+l+1) * zza * zlma(l*l+1)					! zlm(l+1,-l)
        do m = 0, l-1
            zlma((l+1)*(l+2)+m+1)=ri(l-m+1)*(re(l+l+1)*zza*zlma(l*(l+1)+m+1)-(l+m)*rra2*zlma((l-1)*l+m+1))	! zlm(l+1,m)
            zlma((l+1)*(l+2)-m+1)=ri(l-m+1)*(re(l+l+1)*zza*zlma(l*(l+1)-m+1)-(l+m)*rra2*zlma((l-1)*l-m+1))	! zlm(l+1,-m)
        end do
    end do
    return
    end
! This file has been generated with the notebook:
!	<</home/rafa/math/notebooks/pargamma25_D_3.nb>>
! and modified to introduce a factor r^n
!	Functions fv(n) = r^n * Integrate[Exp[-x*t]* t**(n-1/2),{t,0,1}]
!		= r^n * ( Gamma(n+1/2) - Gamma(n+1/2,x) ) * x**(Gamma(-n-1/2)
!               = r^n * gamma(n+1/2,x) * x**(Gamma(-n-1/2)
!	 with 0 <= n <= 20
!  ********************************************************

  subroutine freq20(r, x, fv)
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: i, iz
    real(KREAL) :: ex, fv, r, x, y, z
    dimension fv(0:20)
    z = abs(x)
    if (z .gt. 200.d0) then
       y = 1.d0 / sqrt(z)
       fv(0) = sqrt(pi) * y
       do i = 1, 20
         fv(i) = fv(i-1) * (i+i-1.d0) * 0.5d0 * y * y
       enddo
       go to 3000
    endif
    iz = z
    ex = exp(-x)
    if (iz .lt. 15) then
       if (z .lt. uno) then
 !   Interval  0 <= z <= 1   ( eps = 1.233D-23 )
 !   polynomial approximation of F20(z)
       fv(20) = (4.8780487804878049D-2+z            &
          *(-4.6511627906976744D-2+z*(2.2222222222222222D-2+z*(-7.0921985815602829D-3+z    &
          *(1.7006802721088306D-3+z*(-3.2679738562078113D-4+z*(5.2410901466593768D-5+z     &
          *(-7.2150072107054325D-6+z*(8.7023111885020894D-7+z*(-9.3414605558779118D-8+z    &
          *(9.0351213814743503D-9+z*(-7.9521759542980011D-10+z*(6.4150755129047891D-11+z   &
          *(-4.7343391649764174D-12+z*(3.0436531511475638D-13+z*(-1.3202921548315112D-14))))))))))))))))
       go to 2000
       endif
       go to (5,10,15) iz/5+1
 5      continue
 !      Interval  1 <= z < 5   ( eps = 1.330D-18 )
 !      polynomial approximation of Exp[x] * F20(z)
       fv(20) = ex * (4.8780487804887270D-2+z   &
          *(2.2688598978478804D-3+z*(1.0083821782313802D-4+z*(4.2909877732148515D-6+z      &
          *(1.7514258408881220D-7+z*(6.8681734495339967D-9+z*(2.5926134770018764D-10+z     &
          *(9.3953223491830761D-12+z*(3.3867973398004913D-13+z*(9.6616173655492100D-15+z   &
          *(5.7513908577916034D-16+z*(-6.2419044563656976D-18+z*(1.1998040254144381D-18)))))))))))))
       go to 2000
 10      continue
 !      Interval  5 <= z < 10   ( eps = 2.197D-19 )
 !      polynomial approximation of Exp[x] * F20(z)
       fv(20) = ex * (4.8780489794684994D-2+z   &
          *(2.2688557223858922D-3+z*(1.0084228213483975D-4+z*(4.2885543942049550D-6+z      &
          *(1.7614425798911456D-7+z*(6.5681000211449926D-9+z*(3.2677371327339277D-10+z     &
          *(-2.2026297716826657D-12+z*(1.8688737085070686D-12+z*(-1.4483765011674838D-13+z &
          *(1.2354396010917641D-14+z*(-6.6752016236177459D-16+z*(2.7402216839405727D-17+z  &
          *(-6.7707030795245797D-19+z*(9.1143746166955684D-21)))))))))))))))
       go to 2000
 15      continue
 !      Interval  10 <= z < 15   ( eps = 1.851D-18 )
 !      polynomial approximation of Exp[x] * F20(z)
       fv(20) = ex * (4.8841694978314296D-2+z   &
          *(2.1952612298984881D-3+z*(1.4197437389228347D-4+z*(-9.8786540053342240D-6+z     &
          *(3.5368648978184568D-6+z*(-5.7450024066913360D-7+z*(7.5878158972779751D-8+z     &
          *(-7.5109257899927147D-9+z*(5.7546785517165086D-10+z*(-3.3687901762469887D-11+z  &
          *(1.4918741683111878D-12+z*(-4.8451102282376191D-14+z*(1.0969631911429229D-15+z  &
          *(-1.5550170752478661D-17+z*(1.0620901118859267D-19)))))))))))))))
       go to 2000
    else
       go to (25,35,45,55,65,75,85,95,105) (iz-15)/10+1
 !      If z > 105: asymptotical expression
       fv(20) = 5.406242982335075D17 / sqrt(z)**41
       go to 2000
 25      continue
 !      Interval  15 <= z < 25   ( eps = 6.625D-15 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (3.0503304845140581D-3+y       &
          *(1.1486927303924751D-1+y*(1.3891510674628955D2+y*(-9.5722910256661872D3+y       &
          *(5.7578266413379715D5+y*(-2.3718074066048430D7+y*(7.5220959099166695D8+y        &
          *(-1.8095278252442547D10+y*(3.3759474357408302D11+y*(-4.8252519908180661D12+y    &
          *(5.2571946506189288D13+y*(-4.2249439839665580D14+y*(2.4081147377852433D15+y     &
          *(-8.7280423804452341D15+y*(1.6439985245835246D16)))))))))))))))
       go to 2000
 35      continue
 !      Interval  25 <= z < 35   ( eps = 3.089D-18 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (5.8747691846974363D-5+y       &
          *(9.7764465514125199D-1+y*(2.3413516830415170D1+y*(-5.6335933979497473D1+y       &
          *(3.6483404782967801D4+y*(-1.4583609380850435D6+y*(6.1127998541254078D7+y        &
          *(-1.6602879996934203D9+y*(3.5553799778962732D10+y*(-5.3096505923262290D11+y     &
          *(5.6870782720347124D12+y*(-3.7315170959956712D13+y*(1.3243661404335882D14)))))))))))))
       go to 2000
 45      continue
 !      Interval  35 <= z < 45   ( eps = 1.547D-19 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (1.5989754165300373D-5+y       &
          *(9.9312556185710981D-1+y*(2.0838150065672113D1+y*(2.0523027473807872D2+y        &
          *(1.8285896271688041D4+y*(-5.3545710391705998D5+y*(2.5719413843168101D7+y        &
          *(-6.1308064271060475D8+y*(1.1700760809463463D10+y*(-1.2321136681771141D11+y     &
          *(7.6107340878398950D11)))))))))))
       go to 2000
 55      continue
 !      Interval  45 <= z < 55   ( eps = 4.160D-20 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (1.7437254480506169D-5+y       &
          *(9.9235105114453346D-1+y*(2.0982789209186459D1+y*(1.9436556162147273D2+y        &
          *(1.8182243066041124D4+y*(-4.5064318513701511D5+y*(1.8402159511206095D7+y        &
          *(-2.8453651344225732D8+y*(3.1367646260747387D9)))))))))
       go to 2000
 65      continue
 !      Interval  55 <= z < 65   ( eps = 3.182D-20 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (4.4813294293896279D-5+y       &
          *(9.8176159665610593D-1+y*(2.2650698969190298D1+y*(6.2837162193241420D1+y        &
          *(2.2790603420436922D4+y*(-4.1574900536061811D5+y*(9.5985867522416429D6)))))))
       go to 2000
 75      continue
 !      Interval  65 <= z < 75   ( eps = 5.328D-20 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (2.0651610448975019D-4+y       &
          *(9.3080337375748514D-1+y*(2.8584764207936303D1+y*(-2.1182416824959043D2+y       &
          *(2.2555785714611651D4)))))
       go to 2000
 85      continue
 !      Interval  75 <= z < 85   ( eps = 1.660D-19 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (1.4546785444263069D-3+y       &
          *(6.7563946245566736D-1+y*(4.1886052644826494D1)))
       go to 2000
 95      continue
 !      Interval  85 <= z < 95   ( eps = 6.308D-21 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (-3.7912464139565155D-3+y      &
          *(1.6136443019529395D0))
       go to 2000
 105     continue
 !      Interval  95 <= z < Infinity   ( eps = 1.692D-24 )
 !   polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
       y = uno / z
       fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (-2.9220486969845953D-3+y      &
          *(1.5312675193464893D0))
    endif
 2000 continue
    fv(19) = (ex + z * fv(20) ) * 5.128205128205128D-2
    fv(18) = (ex + z * fv(19) ) * 5.405405405405405D-2
    fv(17) = (ex + z * fv(18) ) * 5.714285714285714D-2
    fv(16) = (ex + z * fv(17) ) * 6.060606060606061D-2
    fv(15) = (ex + z * fv(16) ) * 6.451612903225806D-2
    fv(14) = (ex + z * fv(15) ) * 6.896551724137931D-2
    fv(13) = (ex + z * fv(14) ) * 7.407407407407407D-2
    fv(12) = (ex + z * fv(13) ) * 8.000000000000000D-2
    fv(11) = (ex + z * fv(12) ) * 8.695652173913043D-2
    fv(10) = (ex + z * fv(11) ) * 9.523809523809524D-2
    fv(9)  = (ex + z * fv(10) ) * 1.052631578947368D-1
    fv(8)  = (ex + z * fv(9) )  * 1.176470588235294D-1
    fv(7)  = (ex + z * fv(8) )  * 1.333333333333333D-1
    fv(6)  = (ex + z * fv(7) )  * 1.538461538461538D-1
    fv(5)  = (ex + z * fv(6) )  * 1.818181818181818D-1
    fv(4)  = (ex + z * fv(5) )  * 2.222222222222222D-1
    fv(3)  = (ex + z * fv(4) )  * 2.857142857142857D-1
    fv(2)  = (ex + z * fv(3) )  * 4.000000000000000D-1
    fv(1)  = (ex + z * fv(2) )  * 6.666666666666667D-1
    fv(0)  = (ex + z * fv(1) )  * 2.000000000000000D0
 3000 continue
    z = r
    do i = 1, 20
       fv(i) = fv(i) * z
       z = z * r
    enddo
    return
    end
!
!   ******************************************************************
!
  subroutine emes ( m1, m2, ms, md, ss, sd )
    USE DAM320_D
    USE DAM320_CONST_D, ONLY: uno, cero, umed
    implicit none
    integer(KINT) :: m1, m1a, m2, m2a, ms, md
    real(KREAL) :: s1, s2, s12, ss, sd
    s1 = sign(1,m1)
    s2 = sign(1,m2)
    s12 = s1 * s2
    m1a = abs(m1)
    m2a = abs(m2)
    ms = s12 * ( m1a + m2a )
    md = s12 * abs( m1a - m2a )
    if ( ms.eq.md ) then
        ss = uno
        sd = cero
        return
    endif
    if ( m1.lt.0 .and. m2.lt.0 ) then
        ss = -umed
    else
        ss = umed
    endif
    if ( s12.gt.cero ) then
        sd = umed
    elseif ( md.eq.0 ) then
        sd = cero
    elseif ( sign(1,m1a-m2a) .eq. s1 ) then
        sd = - umed
    else
        sd = umed
    endif
    return
    end
!**********************************************************************
!
!   subroutine rotar
!
!	this subroutine yields the rotation matrices rl(m',m;l) of reals spherical harmonics
!	receives the trigonometric functions of Euler angles defining the rotation
!
!**********************************************************************
  subroutine rotar(lmax, cosal, sinal, cosbet, sinbet, cosga, singa)
    USE DAM320_D
    USE DAM320_DATA_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: l, lmax
    real(KREAL) :: cosag, cosal, cosamg, cosbet, cosga, sinag, sinal, sinamg, singa, sinbet, tgbet2
!	Initial matrices d0, r0, d1 and r1
    rl(:,:,:) = cero
    dl(:,:,:) = cero
    dl(0,0,0)  = uno
    rl(0,0,0)  = uno
    if(lmax.eq.0) return
    dl(1,1,1)  = (uno + cosbet) * umed
    dl(1,0,1)  =-sinbet/raiz2
    dl(1,-1,1) = (uno - cosbet) * umed
    dl(0,1,1)  =-dl(1,0,1)
    dl(0,0,1)  = dl(1,1,1)-dl(1,-1,1)
    dl(0,-1,1) = dl(1,0,1)
    dl(-1,1,1) = dl(1,-1,1)
    dl(-1,0,1) = dl(0,1,1)
    dl(-1,-1,1)= dl(1,1,1)
    cosag  = cosal * cosga - sinal * singa
    cosamg = cosal * cosga + sinal * singa
    sinag  = sinal * cosga + cosal * singa
    sinamg = sinal * cosga - cosal * singa
    rl(0,0,1)  = dl(0,0,1)
    rl(1,0,1)  = raiz2 * dl(0,1,1) * cosal
    rl(-1,0,1) = raiz2 * dl(0,1,1) * sinal
    rl(0,1,1)  = raiz2 * dl(1,0,1) * cosga
    rl(0,-1,1) =-raiz2 * dl(1,0,1) * singa
    rl(1,1,1)  = dl(1,1,1) * cosag - dl(1,-1,1) * cosamg
    rl(1,-1,1) =-dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
    rl(-1,1,1) = dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
    rl(-1,-1,1)= dl(1,1,1) * cosag + dl(1,-1,1) * cosamg
!	the remaining matrices are calculated using symmetry and recurrence relations by means of the subroutine dlmn.
    if ( abs(sinbet) .lt. 1.d-14 ) then
        tgbet2 = cero
    elseif ( abs(sinbet) .lt. 1.d-10 ) then
        tgbet2 = cero
        write(6,"('WARNING in ROTAR: sinbet = ', e17.10, ' takes  0')") sinbet
    else
        tgbet2 = ( uno - cosbet ) / sinbet
    endif
    do l = 2, lmax
        call dlmn(l, sinal, cosal, cosbet, tgbet2, singa, cosga)
    enddo
    return
    end
!**********************************************************************
!
!   subroutine dlmn
!
!   this subroutine generates the matrices dl(m',m;l) for a fixed value
!   of the orbital quantum number l, and it needs the dl(l-2;m',m) and
!   dl(l-1;m',m) matrices. this subroutine uses symmetry and recurrence
!   relations. the matrices dl(m',m;l) are the rotation matrices for
!   complex spherical harmonics
!
!**********************************************************************
  subroutine dlmn(l, sinal, cosal, cosbet, tgbet2, singa, cosga)
    USE DAM320_D
    USE DAM320_DATA_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: iinf, isup, l, m, mp
    real(KREAL) :: al, al1, ali, aux, cosag, cosagm, cosal, cosaux, cosbet, cosga, cosmal, cosmga, cux, d1, d2
    real(KREAL) :: sgn, sinag, sinagm, sinal, singa, sinmal, sinmga, tal1, tgbet2
    iinf=1-l
    isup=-iinf
!	computation of the dl(m',m;l) matrix, mp is m' and m is m.
!	first row by recurrence: see equations 19 and 20 of reference (6)
    dl(l,l,l) = dl(isup,isup,l-1) * (uno + cosbet) * umed
    dl(l,-l,l) = dl(isup,-isup,l-1) * (uno - cosbet) * umed
    do m = isup, iinf, -1
         dl(l,m,l) = -tgbet2 * root(l+m+1) * rooti(l-m) * dl(l,m+1,l)
    enddo
!	the rows of the upper quarter triangle of the dl(m',m;l) matrix see equation 21 of reference (6)
    al = l
    al1 = al - uno
    tal1 = al + al1
    ali = uno / al1
    cosaux = cosbet * al * al1
    do mp = l-1, 0, -1
        aux = rooti(l+mp) * rooti(l-mp) * ali
        cux = root(l+mp-1) * root(l-mp-1) * al
        do m = isup, iinf, -1
            dl(mp,m,l) = aux * rooti(l+m) * rooti(l-m) * (tal1 * (cosaux - re(m) * re(mp)) * dl(mp,m,l-1) &
                    - root(l+m-1) * root(l-m-1) * cux * dl(mp,m,l-2) )
        enddo
        iinf=iinf+1
        isup=isup-1
    enddo
!	the remaining elements of the dl(m',m;l) matrix are calculated using the corresponding symmetry relations:
!		reflection ---> ((-1)**(m-m')) dl(m,m';l) = dl(m',m;l), m'<=m
!		inversion ---> ((-1)**(m-m')) dl(-m',-m;l) = dl(m',m;l)
!	reflection
    sgn = uno
    iinf = -l
    isup = l-1
    do m = l, 1, -1
        do mp = iinf, isup
            dl(mp,m,l) = sgn * dl(m,mp,l)
            sgn = -sgn
        enddo
        iinf=iinf+1
        isup=isup-1
    enddo
!	inversion
    iinf=-l
    isup=iinf
    do m = l-1, -l, -1
        sgn = -uno
        do mp = isup, iinf,- 1
            dl(mp,m,l) = sgn * dl(-mp,-m,l)
            sgn = -sgn
        enddo
        isup=isup+1
    enddo
!	computation of the rotation matrices rl(m',m;l) for real spherical harmonics using the matrices dl(m',m;l)
!	for complex spherical harmonics: see equations 10 to 18 of reference (6)
    rl(0,0,l) = dl(0,0,l)
    cosmal = cosal
    sinmal = sinal
    sgn = - uno
    do mp = 1, l
        cosmga = cosga
        sinmga = singa
        aux = raiz2 * dl(0,mp,l)
        rl(mp,0,l) = aux * cosmal
        rl(-mp,0,l)= aux * sinmal
        do m = 1, l
            aux = raiz2 * dl(m,0,l)
            rl(0,m,l) = aux * cosmga
            rl(0,-m,l)=-aux * sinmga
            d1 = dl(-mp,-m,l)
            d2 = sgn * dl(mp,-m,l)
            cosag = cosmal * cosmga - sinmal * sinmga
            cosagm= cosmal * cosmga + sinmal * sinmga
            sinag = sinmal * cosmga + cosmal * sinmga
            sinagm= sinmal * cosmga - cosmal * sinmga
            rl(mp,m,l)  = d1 * cosag + d2 * cosagm
            rl(mp,-m,l) =-d1 * sinag + d2 * sinagm
            rl(-mp,m,l) = d1 * sinag + d2 * sinagm
            rl(-mp,-m,l)= d1 * cosag - d2 * cosagm
            aux    = cosmga * cosga - sinmga * singa
            sinmga = sinmga * cosga + cosmga * singa
            cosmga = aux
        enddo
        sgn = - sgn
        aux    = cosmal * cosal - sinmal * sinal
        sinmal = sinmal * cosal + cosmal * sinal
        cosmal = aux
    enddo
    return
    end
!
!   *******************************************************************
!
  subroutine acof
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: k1, k2, k20, k200, kk, kk0, kk00, l, lp, m, m1, mp, n
    real(KREAL) :: aux, bux
    app = cero
!
!   starting elements app(00,lm)(n) = delta(l,n)
!
    k1 = 0
    do l = 0 , mxl
        do m = 0 , l
            kk = ind(k1)
            app(l,kk) = uno
            k1 = k1 + 1
        enddo
    enddo
!
!   elements app(lm,m'm')(n)
!
    do mp = 1 , mxl
        k2 = ind(mp) + mp
        k20 = ind(mp-1) + mp-1
        do l = mp , mxl
            if ( l.eq.mp ) then
                m1 = mp
            else
                m1 = 0
            endif
            do m = m1 , l
                k1 = ind(l) + m
                kk = ind(k1) + k2
                kk0 = ind(k1) + k20
                do n = l-mp , l+mp , 2
                    if ( n.ge.m+mp) then
                        app(n,kk) = (2*mp-1) * ( app(n-1,kk0) * ri(n+n-1) - app(n+1,kk0) * ri(n+n+3) )
                    endif
                enddo
            enddo
        enddo
    enddo
!
!   elements app(lm,l'm')(n)
!
    do mp = 0 , mxl
        k200 = 0
        do lp = mp+1 , mxl
            k2 = ind(lp) + mp
            k20 = ind(lp-1) + mp
            if ( lp.gt.mp+1 ) k200 = ind(lp-2) + mp
            do l = lp , mxl
                if ( l.eq.lp ) then
                    m1 = mp
                else
                    m1 = 0
                endif
                do m = m1 , l
                    k1 = ind(l) + m
                    kk = ind(k1) + k2
                    kk0 = ind(k1) + k20
                    kk00 = ind(k1) + k200
                    do n = l-lp , l+lp , 2
                        if ( n.ge.m+mp) then
                            aux = app(n+1,kk0) * re(n+m+mp+1) * dosl1i(n+1)
                            if ( n.gt.m+mp ) aux = aux + app(n-1,kk0) * re(n-m-mp) * dosl1i(n-1)
                            aux = aux * dosl1(lp-1)
                            if ( lp.gt.mp+1 ) aux = aux - re(lp+mp-1) * app(n,kk00)
                            app(n,kk) = aux * ri(lp-mp)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    return
    end
!
!   *******************************************************************
!
  subroutine bcof
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: k1, k2, k20, k200, kk, kk0, kk00, l, lp, m, m1, mmp, mp, n
    real(KREAL) :: aux, bux, t1, t2
    bpp = cero
!
!   starting elements bpp(lm,00)(n) = delta(l,n)
!
    k1 = 0
    do l = 0 , mxl
        do m = 0 , l
            kk = ind(k1)
            bpp(l,kk) = uno
            k1 = k1 + 1
        enddo
    enddo
!
!   elements bpp(lm,m'm')(n)
!
    do mp = 1 , mxl
        k2 = ind(mp) + mp
        k20 = ind(mp-1) + mp-1
        do l = mp , mxl
            if ( l.eq.mp ) then
                m1 = mp
            else
                m1 = 0
            endif
            do m = m1 , l
                k1 = ind(l) + m
                kk = ind(k1) + k2
                kk0 = ind(k1) + k20
                do n = l-mp , l+mp , 2
                    if ( mp.gt.m ) then
                        t1 = uno
                        t2 = uno
                    else
                        t1 = -re(n-(m-mp+1)) * re(n-(m-mp+1)+1)
                        t2 = -re(n+(m-mp+1)) * re(n+(m-mp+1)+1)
                    endif
                    if ( n.ge.abs(m-mp)) then
                        if (n.eq.0) then
                            bux=cero
                        else
                            bux=t1*bpp(n-1,kk0) * dosl1i(n-1)
                        endif
                        bpp(n,kk) = dosl1(mp-1) * ( bux - t2 * bpp(n+1,kk0) * dosl1i(n+1) )
                    endif
                enddo
            enddo
        enddo
    enddo
!
!   elements bpp(lm,l'm')(n)
!
    do mp = 0 , mxl
        k200 = 0
        do lp = mp+1 , mxl
            k2 = ind(lp) + mp
            k20 = ind(lp-1) + mp
            if ( lp.gt.mp+1 ) k200 = ind(lp-2) + mp
            do l = lp , mxl
                if ( l.eq.lp ) then
                     m1 = mp
                else
                    m1 = 0
                endif
                do m = m1 , l
                    k1 = ind(l) + m
                    kk = ind(k1) + k2
                    kk0 = ind(k1) + k20
                    kk00 = ind(k1) + k200
                    do n = l-lp , l+lp , 2
                        mmp = abs(m-mp)
                        if ( n.ge.mmp) then
                            aux = bpp(n+1,kk0) * re(n+mmp+1) * dosl1i(n+1)
                            if ( n.gt.mmp ) aux = aux + bpp(n-1,kk0) * re(n-mmp) * dosl1i(n-1)
                            aux = aux * dosl1(lp-1)
                            if ( lp.gt.mp+1 ) aux = aux - re(lp+mp-1) * bpp(n,kk00)
                            bpp(n,kk) = aux * ri(lp-mp)
                        endif
                    enddo
                enddo
            enddo
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
    integer(KINT) :: ierr, ierr2
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    call MPI_FINALIZE(ierr2)
    stop
    end
