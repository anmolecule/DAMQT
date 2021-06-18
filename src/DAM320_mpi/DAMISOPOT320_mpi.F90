!  Copyright 2013-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
! Program for computing the molecular electrostatic potential (mesp) isosurface and its normals 
! from the representation of the molecular density performed with DAM320
!
!
! Version of October 2019
!

! #define DBLPRCGRID    ! Uncomment this line  if double precision grid is wanted
!===============================================================================================
!                 MODULE DAM2017_CONST_D
!=============================================================================================== 
MODULE DAMISOPOT320_D
    USE DAM320_D
    USE DAM320_CONST_D
    IMPLICIT NONE
    logical :: lbinary, lexist
    character(300) :: gridname, indfile, outisopotname, outrootname, straux, vertfile
    integer(KINT) :: i, i1, i2, ia, iaux, ierr, icube, iopt(5), iuni, ix, iy, iz, j, k, knt, kntgrid
    integer(KINT) :: kntind, kntvert, l, npoints, nx, ny, nz
    integer(KINT), allocatable ::  indices(:), interpolmat(:,:), nind(:), nvert(:), tritable(:,:)
    real(KREAL), parameter :: angstromtobohr = 1.889725989d0
    real(KREAL) :: a, aux, b, bux, c, contourval, cux, d, drvx, disthresq, drvxtot
    real(KREAL) :: drvy, drvytot, drvz, drvztot, errabs, s, surftot, surftrian
    real(KREAL) :: vaux, vel, vnucl, volume, voltetrahed, volvoxel, vtot
    real(KREAL) :: x, xini, xinterp, xfin, y, yini, yinterp, yfin, z, zini, zinterp, zfin
    real(KREAL), allocatable :: gradient(:), grid(:), fvoxel(:), vertices(:,:), vertices2(:)
    real(KREAL), allocatable :: xvoxel(:), yvoxel(:), zvoxel(:)
    real(KREAL) :: xyz(3), xyztetr(3,0:3)
    real(KREAL4), allocatable :: grid4(:)
#ifdef DBLPRCGRID
    real(KREAL) :: aux4, bux4, drvxtot4, drvytot4, drvztot4
    real(KREAL) :: x4, xini4, xfin4, y4, yini4, yfin4, z4, zini4, zfin4 
    real(KREAL), allocatable :: gradient4(:), vertices4(:,:)
    integer(KINT), parameter :: i0 = 0
#else
    real(KREAL4) :: aux4, bux4, drvxtot4, drvytot4, drvztot4
    real(KREAL4) :: x4, xini4, xfin4, y4, yini4, yfin4, z4, zini4, zfin4
    real(KREAL4), allocatable :: gradient4(:), vertices4(:,:)
    integer(KINT), parameter :: i0 = 3
#endif
END MODULE
!
!                 END OF MODULE DAMISOPOT320_D
!...............................................................................................
  program DAMISOPOT320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMPOT320_D
    USE DAMISOPOT320_D
    USE PARALELO
    implicit none
    integer(KINT) :: jaux
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    logical :: lnamelist(5), ltimeprocs
    integer(KINT) :: inamelist(1)
    real(KREAL) :: xmax, xmin, ymax, ymin, zmax, zmin
    real(KREAL) :: rnamelist(3)

    
    namelist / options / contourval, filename, gridname, geomthr, iswindows, langstrom, lbinary, lmaxrep, &
        longoutput, lvalence, umbrlargo
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    abort = 0
    abortroot = 0
    tiempo = dtime(tarray)
!    Namelist default values
    contourval = 1.d-3      ! value of mesp for isosurface
    filename = ""           ! root file name for output surface files (*.srf, *.sgh)
    geomthr = 1.d-5         ! Geometry threshold: two points at a distance lower than geomthr are considered to be the coincident
    gridname = ""           ! name of file with mesp grid
    langstrom = .true.      ! If false original grid distances in bohr
    lbinary = .true.        ! If true writes a file *.srf with the surface in binary form, otherwise writes a file *.srf_txt text mode
    lmaxrep = 5             ! highest "l" in the expansion of the density and potential
    longoutput = .false.    ! If true a more detailed output is given
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
            gridname = trim(projectname)//"-v.plt"
        else
            if (gridname(len_trim(gridname)-5:len_trim(gridname)) .ne. '-v.plt') then
                gridname = trim(gridname)//"-v.plt"
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
        
        lnamelist(1) = langstrom
        lnamelist(2) = lbinary
        lnamelist(3) = longoutput
        lnamelist(4) = iswindows
        lnamelist(5) = lvalence

        inamelist(1) = lmaxrep

        rnamelist(1) = contourval
        rnamelist(2) = geomthr 
        rnamelist(3) = umbrlargo
    endif
    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(gridname,len(gridname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lnamelist,5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (myrank .ne. 0) then
        langstrom  = lnamelist(1) 
        lbinary    = lnamelist(2) 
        longoutput = lnamelist(3) 
        iswindows  = lnamelist(4) 
        lvalence   = lnamelist(5) 

        lmaxrep    = inamelist(1)

        contourval = rnamelist(1)
        geomthr    = rnamelist(2)  
        umbrlargo  = rnamelist(3) 
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

    if (myrank .eq. 0) write(6,"('Potential from expansion of the density: lmaxrep = ', i3)") lmaxrep
  
    call readdamqtisopot        !    Reads file .damqt  (generated by DAM2016)
    
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
    
    allocate (gradient(3*kntvert), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error ' i5, ' allocating gradient in processor ', i3)") ierr, myrank
        abort= 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
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
!     Computes MESP and MESP gradient on the vertices
    errabs = cero
    gradient = cero
    knt = 0
    do i = 1, kntvert
        x = vertices(1,i)
        y = vertices(2,i)
        z = vertices(3,i)
        drvxtot = cero
        drvytot = cero
        drvztot = cero
        vnucl = cero
        vel = cero
        vtot = cero
        do ia = 1, ncen
            call mespdrv(ia, x, y, z, vnucl, vel, vaux, drvx, drvy, drvz)
            vtot = vtot + vaux
            drvxtot = drvxtot + drvx
            drvytot = drvytot + drvy
            drvztot = drvztot + drvz
        enddo
        errabs = max(errabs, abs(contourval-vtot))
        
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
        write(iuni) x4, y4, z4, drvxtot4, drvytot4, drvztot4
    enddo

!     writes triangles indices to temporal file in unit iuni+1
    write(iuni+1) indices(1:kntind)

    close(iuni)
    close(iuni+1)
    deallocate (gradient, indices, vertices)
    if (myrank .ne. 0) tiempo = dtime(tarray)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! Necessary for the processors to end writing temporal files
!     Gathers files generated by processors
    if (myrank .eq. 0) then
        kntvert = 0
        kntind = 0
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
            allocate (indices(kntind), gradient4(3*kntvert), vertices4(3,kntvert), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Error ' i5, ' allocating indices, gradient4 and vertices4 in processor ', i3)")ierr, myrank
                abort= 1
            endif
        endif
        if (abort .eq. 0) then
            i1 = 0
            i2 = 0
            knt = 0
            do i = 0, nprocs-1
                do j = 1, nvert(i)
                    read(iuni+2*i) vertices4(:,i1+1), gradient4(3*i1+1:3*i1+3)
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

            if (lbinary) then
                outisopotname = trim(outrootname)//"_"//trim(adjustl(straux))//".isopot"
#if _WIN32         
                open (unit=iuni+2*nprocs+1, file=outisopotname, form='binary', carriagecontrol='NONE', iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outisopotname)
                    abort= 1
                endif
#elif __INTEL_COMPILER
                open (unit=iuni+2*nprocs+1, file=outisopotname, form='binary', carriagecontrol='NONE', iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' opening file ',a)") ierr, trim(outisopotname)
                    abort= 1
                endif
#else
                open (unit=iuni+2*nprocs+1, file=outisopotname, form='unformatted', access='stream', iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' reading file ',a)") ierr, trim(outisopotname)
                    abort= 1
                endif
#endif
            else
                outisopotname = trim(outrootname)//"_"//trim(adjustl(straux))//".isopot_txt"
                open (unit=iuni+2*nprocs+1, file=outisopotname, iostat=ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error ', i5,' reading file ',a)") ierr, trim(outisopotname)
                    abort= 1
                endif
            endif
        endif
        
        if (abort .eq. 0) then
            if (lbinary) then
                write(iuni+2*nprocs+1) i0, 0, nx, ny, nz
                xini4 = xini ; xfin4 = xfin ; yini4 = yini ; yfin4 = yfin ; zini4 = zini ; zfin4 = zfin
                write(iuni+2*nprocs+1) xini4, xfin4, yini4, yfin4, zini4, zfin4
                write(iuni+2*nprocs+1) kntvert, kntind
            else
                write(iuni+2*nprocs+1,"(5(i6))") i0, 0, nx, ny, nz
                write(iuni+2*nprocs+1,"(6(1x,e14.7))") xini, xfin, yini, yfin, zini, zfin
                write(iuni+2*nprocs+1,"(i10)") kntvert, kntind
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
                    write(iuni+2*nprocs+1) x4, y4, z4, drvxtot4, drvytot4, drvztot4
                else
                    write(iuni+2*nprocs+1,"(7(1x,e13.6))") x, y, z, drvxtot, drvytot, drvztot
                endif
            enddo
        endif
        
!     Writes triangles indices to file
        if (abort .eq. 0) then
            if (lbinary) then
                write(iuni+2*nprocs+1) indices(1:kntind)
            else
                write(iuni+2*nprocs+1,"(21(1x,i7))") indices(1:kntind)
            endif
            write(6,"(/3x,i9,' triangles generated for contour ',e12.5)") kntind/3, contourval 
            write(6,"(3x,'Number of vertices = ',i9)") kntvert
            write(straux,"(e9.2)") contourval
            write(6,"(/'For mesp contour ', e12.5,':', &
                &/5x,'Molecular volume / bohr^3 = ', e12.5, 5x,' Molecular volume / A^3 = ', e12.5, &
                &/5x,'Molecular surface / bohr^2 = ', e12.5, 5x,' Molecular surface / A^2 = ', e12.5)") &
                contourval, volume, volume*0.148184534296, surftot, surftot*0.280028297329d0
            write(6,"(/'Highest error in potential, absolute = ', e10.3, ' relative = ', e10.3,/)") errabs, errabs / contourval 
        endif
        vertfile = trim(outrootname)//"_*.vrttmp"
        indfile = trim(outrootname)//"_*.indtmp"
        call system("rm -f "//trim(vertfile))
        call system("rm -f "//trim(indfile))
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
  subroutine readdamqtisopot
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, ncenbas, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, fr1, fr2l2, pi4d2l1, r, ra, ral1inv, rainv, ral
    real(KREAL) :: rinta, rintb, rlarex, step, stepmed, suml, suml1, suml2, summ, summ1, summ2, t
    real(KREAL) :: tcheb(0:mxlenpol-1)
    inquire(file=trim(projectname)//"_2016.damqt", size=nsize, iostat=ierr)
    if (ierr .ne. 0) call error(ierr,'Error when inquiring file '//trim(projectname)//"_2016.damqt")
    if (nsize .eq. -1) call error(1,'Size of file '//trim(projectname)//"_2016.damqt cannot be determined")
    if (longoutput) write(6,"('Size of file ', a, ' = ', i12)") trim(projectname)//"_2016.damqt", nsize
#if _WIN32
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.damqt")
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.dmqtv")
#elif __INTEL_COMPILER
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.damqt")
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.dmqtv")
#else
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='unformatted', action = 'read', access='stream', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.damqt")
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='unformatted', action = 'read', access='stream', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//"_2016.dmqtv")
#endif
    if (longoutput) write(6,"('Opens files ', a, ' and ', a)") trim(projectname)//"_2016.damqt", trim(projectname)//"_2016.dmqtv"
    read(10) ncen, nbas, ncaps
    nsize = nsize - sizeof(ncen) - sizeof(nbas) - sizeof(ncaps)
    write(6,"('ncen = ', i4, ' nbas = ', i6, ' ncaps = ', i5)") ncen, nbas, ncaps

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

!    Reads auxiliary integrals from file .dmqtv

    allocate(cfrint1(icfposd(lmtop*nintervaj+1,ncen)-1), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating cfrint1. Stop')
    if (longoutput) write(6,"('Size of cfrint1   = ', i15, ' bytes')") size(cfrint1)

    allocate(cfrint2l2(icfposd(lmtop*nintervaj+1,ncen)-1), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating cfrint2l2. Stop')
    if (longoutput) write(6,"('Size of cfrint2l2 = ', i15, ' bytes')") size(cfrint2l2)

    allocate(QGacum(nintervaj*lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating QGacum. Stop')
    if (longoutput) write(6,"('Size of QGacum    = ', i15, ' bytes')") size(QGacum)

    allocate(Qgpart(nintervaj*lmtop), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating Qgpart. Stop')
    if (longoutput) write(6,"('Size of Qgpart    = ', i15, ' bytes')") size(Qgpart)

    allocate(qppart(nintervaj*lmtop), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating qppart. Stop')
    if (longoutput) write(6,"('Size of qppart    = ', i15, ' bytes')") size(qppart)

    allocate(qpacum(nintervaj*lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating qpacum. Stop')
    if (longoutput) write(6,"('Size of qpacum    = ', i15, ' bytes')") size(qpacum)

    allocate(rlargo(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rlargo. Stop')

    allocate(rmultip(lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rmultip. Stop')
    if (longoutput) write(6,"('Size of rmultip    = ', i15, ' bytes')") size(rmultip)

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
    if (ierr .ne. 0) call error(1,'Memory error when allocating lcorto, llargo and Qllargo. Stop')
    if (longoutput) write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)
    if (longoutput) write(6,"('Size of llargo   = ', i15, ' bytes')") size(llargo)
    if (longoutput) write(6,"('Size of Qllargo   = ', i15, ' bytes')") size(Qllargo)
!    long-range radii

    allocate(umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating umedpow in readdamqtisopot. Stop')
    umedpow(0) = uno                            !
    do i = 1, lmaxexp                            !
            umedpow(i) = umedpow(i-1) * umed            ! 1 / 2^i
    enddo
    write(6,"('Long-range threshold = ', e12.5)") umbrlargo
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
        if (longoutput) then
            write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5, ' lcorto = ', 30(i3))") &
                    ia, atmnms(nzn(ia)), rlargo(ia), lcorto(1:nintervaj,ia)
        else
            write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5)") &
                    ia, atmnms(nzn(ia)), rlargo(ia)
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
        if (longoutput) write(6,"('llargo: ', 51(i3))") llargo(0:mxlargo,ia)
    enddo
    deallocate (Qgpart, qppart)
    deallocate(umedpow)
    return
    end
!
!   ***************************************************************
!     Calculates the electrostatic potential and its gradient from the represented density at point (x,y,z)
!
   subroutine mespdrv(ia, x, y, z, vnucl, vel, vtot, drvx, drvy, drvz)
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    implicit none
    integer(KINT) :: i, ia, ierr, interv, icflm, j, jshft, kntlm, l, lm, lmtopot, ltop, m
    real(KREAL) :: aux, bux, dost, drvx, drvy, drvz, drvvlm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, flm, fux
    real(KREAL) :: pi4exp, ra, ra2, rainv, ra2inv, rinta, rintb, sgn, t, tp, vnucl, vel, vlm, vqlm, vtot
    real(KREAL) :: x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL), parameter :: vtope = 1.d10  ! To prevent infinity, if the point coincides with a nucleus
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
    drvx = cero
    drvy = cero
    drvz = cero
    dxx = cero
    dxy = cero
    dxz = cero
    dyy = cero
    dyz = cero
    dzz = cero
    if (ra2 .lt. geomthr*geomthr) then
        vnucl = vtope
        vtot = vtope
        if (lgradient) then
            drvx = -vtope
            drvy = -vtope
            drvz = -vtope
        endif
        if (ngini(ia) .eq. 0) return
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
    if (ngini(ia) .eq. 0) then
        vnucl = zn(ia) * rainv
        vtot = vnucl
        if (lgradient) then
            drvx = -vnucl * xa * ra2inv
            drvy = -vnucl * ya * ra2inv
            drvz = -vnucl * za * ra2inv
        endif
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
    call derivzlm(ltop, idimzlm, zlma, zlmadx, zlmady, zlmadz)
    if (ra .ge. rlargo(ia)) then  ! The point is in the long-range region (ra >= rlargo(ia) )
        kntlargo = kntlargo + 1
        aux = zn(ia)    ! aux is the nuclear charge for l = 0 and zero otherwise
        do lm = 1, lmtopot
            vel = vel - rmultip(lm,ia) * zlma(lm) * ra2l1inv(lm)
            vlm = (aux - rmultip(lm,ia) ) * ra2l1inv(lm)
            vtot = vtot + vlm * zlma(lm)
            drvvlm = -d2l1(lm) * vlm * rainv
            drvx = drvx + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
            drvy = drvy + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
            drvz = drvz + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
            aux = cero
        enddo
    else        !     The point is in the short-range region (ra < rlargo(ia) )
        kntcorto = kntcorto + 1
        t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
        dost = t + t
        tcheb(0) = uno    ! Chebyshev T  polynomials
        tcheb(1) = t
        do j = 2, mxlenpol-1
                tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
        enddo
        aux = zn(ia)    ! aux is the nuclear charge for l = 0 and zero otherwise
        pi4exp = pi4exp * exp(-xajustd(interv,ia)*ra)
        do lm = 1, lmtopot
            if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo2) cycle
            icflm = icfposd((interv-1)*lmtop+lm,ia)
            flm = cero
            rinta = cero
            rintb = cero
            do i = 0, icfposd((interv-1)*lmtop+lm+1,ia)-icflm-1
                    rinta = rinta + cfrint2l2(icflm+i) * tcheb(i)
                    rintb = rintb + cfrint1(icflm+i) * tcheb(i)
                    flm   = flm + cfajust(icflm+i) * tcheb(i)
            enddo
            rinta = rinta * (ra-rinterv(interv-1))
            if (interv .gt. 1) rinta = QGacum((interv-2)*lmtop+lm,ia) + rinta
            rintb = qpacum((interv-1)*lmtop+lm,ia) + rintb * (rinterv(interv)-ra)
            vel = vel - pi4d2l1(lm) * ( rinta + ra2l1(lm) * rintb) * zlma(lm) * ra2l1inv(lm)
            vlm = (aux - pi4d2l1(lm) * ( rinta + ra2l1(lm) * rintb)) * ra2l1inv(lm)
            vqlm = (aux - pi4d2l1(lm) * rinta) * ra2l1inv(lm)
            vtot = vtot + vlm * zlma(lm)
            drvvlm = -d2l1(lm) * vqlm * rainv     ! d2l1(lm) = (l+l+1)
            drvx = drvx + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
            drvy = drvy + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
            drvz = drvz + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
            aux = cero
        enddo
    endif  ! End of test over long/short-range
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
    integer(KINT) :: ierr, ierr2
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    call MPI_FINALIZE(ierr2)
    stop
    end
