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
! Program for computing density gradient lines from the representation of the molecular density performed with
! DAM320. Version parallelized with MPI only 3D lines.
!
!
! Version of September 2018
!
  program DAMDENGRAD320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    USE icosahedron
    USE PARALELO
    implicit none
    character(1) :: cpslabel
    integer(KINT) :: i, ia, icnt, ierr, ioerr, ioplines3D, ipnt, itheta, j, knt, kntinf, kntsup
    integer(KINT) :: naux, nlincomptot, nlineastot, ncntplane, nlinpernuc, npntot, npuntot, numpnt
    real(KREAL) :: aux, R, u0, ucnt, v0, vcnt, usalto
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz
    real(KREAL) :: dxdu, dxdv, dydu, dydv, dzdu, dzdv, ex, ey, ez, gu, gv, x0, xcnt, vmod, y0, ycnt, z0, zcnt
    real(KREAL) :: amat(2,3), bmat(2,3), eigval(2), eigvec(2,2), hxyz(3,3), huv(2,2)
    real(KREAL), allocatable :: rcntplane(:,:)
    real(KREAL), parameter :: angs2bohr = 1.889725989d0
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    logical :: lnamelist(1), ltimeprocs
    integer(KINT) :: inamelist(3)
    real(KREAL) :: rnamelist(9)
    character(256) :: filelines
    logical :: lnucleo, lexist
    character*4 :: strbux
    namelist / options / dlt0, filename, filelines, icntlines, ioplines3D, iswindows, lextralines, lmaxrep, longoutput &
            , lplot2d, nlines, nlinpernuc, numpnt, planeA, planeB, planeC &
            , thresh, umbrlargo, usalto, uvratio, rlines &
            , vmod, xinf, xsup, yinf, ysup, zinf, zsup, uinf, usup, vinf, vsup
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    abort = 0
    abortroot = 0
    if (myrank .eq. 0) write(6,"('number of processors = ', i3)") nprocs
    tiempo = dtime(tarray)
!     Namelist default values
    ioplines3D = 1      ! Set of lines per nucleus for 3D plots based on icosahedron vertices, C2 axes and C3 axes or combinations 
                        ! of them:
                        ! 1: vertices (12 points);   2: C3 axes (20 points);   3: C2 axes (30 points)
                        ! 4: vertices + c3 ;   5: vertices + C2;   6: C2 + C3;   7: vertices + C2 + C3 
    longoutput = .false.! If .true. a more detailed output is given
    iswindows = .false.          ! .true. if running on a MS-windows system
    lmaxrep = 5         ! highest value of  l  in the expansion of the density for computing the density gradient
    umbrlargo = 1.d-8   ! Threshold for determining the short-range radius
    usalto = 1.d-3      ! Threshold for convergence in jump
    numpnt = 2000       ! Maximum number of points in each density gradient line
    lextralines = .false. ! If .true. reads extra lines from a file
    lplot2d = .false.   ! Not used. Kept for compatibility of input files
    nlinpernuc = 16     ! Kept for compatibility of input files
    filename = ""       ! root file name for .dengr files
    filelines = ""      ! File with starting points for density gradient lines
    xinf = cero
    xsup = cero
    yinf = cero
    ysup = cero
    zinf = cero
    zsup = cero
    dlt0 = 1.d-2
    planeA = cero       ! Not used. Kept for compatibility of input files
    planeB = cero       ! Not used. Kept for compatibility of input files
    planeC = uno        ! Not used. Kept for compatibility of input files
    uinf = cero         ! Not used. Kept for compatibility of input files
    usup = uno          ! Not used. Kept for compatibility of input files
    uvratio = uno       ! Not used. Kept for compatibility of input files
    vinf = cero         ! Not used. Kept for compatibility of input files
    vsup = uno          ! Not used. Kept for compatibility of input files
    nlines = 0          ! Number of starting points for extra lines in namelist
    icntlines = 0       ! Array with indices for starting nuclei of extra lines (0 for lines starting in a point other than 
                        ! nuclei)
    rlines = cero       ! Array with coordinates (in au) of points defining lines
!     End of namelist defaults
    ltimeprocs = .false.
    if (myrank .eq. 0) then
        read(5,OPTIONS)
        read(5,*) projectname

        if (xsup-xinf .le. cero) then
            abort = 1
        endif
        if (ysup-yinf .le. cero) then
            abort = 1
        endif
        if (zsup-zinf .le. cero) then
            abort = 1
        endif
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Error in grid dimensions. Stop')
    endif
    if (myrank .eq. 0) then
        rnamelist = (/ xinf, xsup, yinf, ysup, zinf, zsup, dlt0, umbrlargo, usalto /)

        ioplines3D = max(0,min(7,ioplines3D))    ! Forces ioplines3D to be in the range [0,7]
        inamelist = (/ lmaxrep, ioplines3D, numpnt /)
        write(6,"(1x,'project name : ',a,/,1x,'==============')") projectname
        if (iswindows) then
            dirsep = "\\"
            i = index(projectname,dirsep,.true.)     ! Checks position of last directory name separator
            if (i .eq. 0) then     ! This is intended for MinGW, whose directory separator in windows is also /
                    dirsep = "/"
                    i = index(projectname,dirsep,.true.)     ! Checks position of last directory name separator
            endif
        else
            dirsep = "/"
            i = index(projectname,dirsep,.true.)     ! Checks position of last directory name separator
        end if
        if (len_trim(filename).eq.0) then
                filename = projectname
        else
                filename = projectname(1:i)//trim(filename)
        endif
        allocate(timeprocs(2*nprocs), stat = ierr)
        if (ierr .eq. 0) then
            ltimeprocs = .true.
            timeprocs = 0.
        else
            write(6,"('WARNING: Memory error when allocating timeprocs, ierr =  ',i5)") ierr
            ltimeprocs = .false.
        endif
        lnamelist = (/ ltimeprocs /)
    endif
    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lnamelist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (myrank .ne. 0) then
        lmaxrep = inamelist(1); ioplines3D = inamelist(2); numpnt = inamelist(3); 
        ltimeprocs = lnamelist(1); 
        xinf = rnamelist(1); xsup = rnamelist(2); 
        yinf = rnamelist(3); ysup = rnamelist(4); 
        zinf = rnamelist(5); zsup = rnamelist(6); dlt0 = rnamelist(7)
        umbrlargo = rnamelist(8); usalto = rnamelist(9)
    endif
    
    call consta          !     Computes auxiliary constants
    call leedamqtden     !     Reads file .damqt  (generated by DAM2016)
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
            call error(1,'Stop')
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (lmaxrep .gt. lmaxexp) then
            write(6,"('lmaxrep = ', i3, ' greater than lmaxexp ', i3)") lmaxrep, lmaxexp
            write(6,"('takes lmaxrep = ',i3)") lmaxexp
            lmaxrep = lmaxexp
    endif
    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma, zlmadx, zlmady, zlmadz in processor ',i3)") myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    call seekdenmax     ! Checks if density maxima are in nuclei, otherwise seeks the actual positions of the maxima
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    if (myrank .eq. 0) write(6,"('Density gradient lines from expansion of the density: lmaxrep = ', i3)") lmaxrep
    if (xsup-xinf .le. cero) then
        xsup = dos * (maxval(rcen(1,1:ncen)) + uno)
        xinf = dos * (minval(rcen(1,1:ncen)) - uno)
    endif
    if (myrank .eq. 0) write(6,"('xinf = ', e17.10, ' xsup = ', e17.10)") xinf, xsup
    if (ysup-yinf .le. cero) then
        ysup = dos * (maxval(rcen(1,1:ncen)) + uno)
        yinf = dos * (minval(rcen(1,1:ncen)) - uno)
    endif
    if (myrank .eq. 0) write(6,"('yinf = ', e17.10, ' ysup = ', e17.10)") yinf, ysup
    if (zsup-zinf .le. cero) then
        zsup = dos * (maxval(rcen(1,1:ncen)) + uno)
        zinf = dos * (minval(rcen(1,1:ncen)) - uno)
    endif
    if (myrank .eq. 0) write(6,"('zinf = ', e17.10, ' zsup = ', e17.10)") zinf, zsup
    
!     Opens file for density gradient lines    
    
    write(strbux,'(i4.2)') myrank
    open (unit=12, file=trim(filename)//".dengr_"//trim(adjustl(strbux)),status='unknown',form='formatted',iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error in processor ',i3, ' when opening file ', a, '.Stop')") myrank, &
            trim(filename)//".dengr_"//trim(adjustl(strbux))
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
    !    density gradient 3D lines

    npntot = 0
    nlineas = 0
    nlincomp = 0
    
    nlinpernuc = 0
    if (ioplines3D .eq. 1 .or. ioplines3D .eq. 4 .or. ioplines3D .eq. 5 .or. ioplines3D .eq. 7) nlinpernuc = nlinpernuc + idimvert
    if (ioplines3D .eq. 2 .or. ioplines3D .eq. 4 .or. ioplines3D .eq. 6 .or. ioplines3D .eq. 7) nlinpernuc = nlinpernuc + idimc3
    if (ioplines3D .eq. 3 .or. ioplines3D .eq. 5 .or. ioplines3D .eq. 6 .or. ioplines3D .eq. 7) nlinpernuc = nlinpernuc + idimc2
    naux = ncen * nlinpernuc / nprocs

    knt = 0
    kntinf = myrank * naux
    if (myrank .eq. nprocs-1) then
        kntsup = ncen * nlinpernuc
    else
        kntsup = (myrank+1) * naux - 1
    endif
    r = dlt0
    do ia = 1, ncen
        xcnt = rdenmax(1,ia)
        ycnt = rdenmax(2,ia)
        zcnt = rdenmax(3,ia)
        if (ioplines3D .eq. 1 .or. ioplines3D .eq. 4 .or. ioplines3D .eq. 5 .or. ioplines3D .eq. 7) then
            if (knt+idimvert .ge. kntinf .and. knt .le. kntsup) then
                do i = 1, idimvert
                    knt = knt + 1
                    if (knt .lt. kntinf .or. knt .gt. kntsup) cycle 
                    x0 = xcnt + r * vertices(i,1)
                    y0 = ycnt + r * vertices(i,2)
                    z0 = zcnt + r * vertices(i,3)
                    nlineas = nlineas + 1
                    write(12,"(3(1x,e12.5))") xcnt, ycnt, zcnt
                    call metgrad(x0, y0, z0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            else
                knt = knt + idimvert
            endif
        endif
        if (ioplines3D .eq. 2 .or. ioplines3D .eq. 4 .or. ioplines3D .eq. 6 .or. ioplines3D .eq. 7) then
            if (knt+idimc3 .ge. kntinf .and. knt .le. kntsup) then
                do i = 1, idimc3
                    knt = knt + 1
                    if (knt .lt. kntinf .or. knt .gt. kntsup) cycle
                    x0 = xcnt + r * c3axes(i,1)
                    y0 = ycnt + r * c3axes(i,2)
                    z0 = zcnt + r * c3axes(i,3)
                    nlineas = nlineas + 1
                    write(12,"(3(1x,e12.5))") xcnt, ycnt, zcnt
                    call metgrad(x0, y0, z0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            else
                knt = knt + idimc3
            endif
        endif
        if (ioplines3D .eq. 3 .or. ioplines3D .eq. 5 .or. ioplines3D .eq. 6 .or. ioplines3D .eq. 7) then
            if (knt+idimc2 .ge. kntinf .and. knt .le. kntsup) then
                do i = 1, idimc2
                    knt = knt + 1
                    if (knt .lt. kntinf .or. knt .gt. kntsup) cycle
                    x0 = xcnt + r * c2axes(i,1)
                    y0 = ycnt + r * c2axes(i,2)
                    z0 = zcnt + r * c2axes(i,3)
                    nlineas = nlineas + 1
                    write(12,"(3(1x,e12.5))") xcnt, ycnt, zcnt
                    call metgrad(x0, y0, z0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            else
                knt = knt + idimc2
            endif
        endif
    enddo
    tiempo = dtime(tarray)
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
    if (myrank .eq. 0 .and. lextralines) then
!       reads the starting points coordinates for extra lines from namelist
        if (nlines .gt. 0) write(6,"('Reads extra lines data from namelist')")
        do i = 1, nlines
            icnt = icntlines(i)
            x0 = rlines(1,i)
            y0 = rlines(2,i)
            z0 = rlines(3,i)
            if (icnt .lt. 1 .or. icntlines(i) .gt. ncen) then
                xcnt = x0
                ycnt = y0
                zcnt = z0
                write(12,'(3(1x,e12.5))') xcnt, ycnt, zcnt
            else
                vmod = sqrt(x0*x0+y0*y0+z0*z0)
                if (vmod .eq. 0.d0) cycle
                vmod = dlt0 / vmod
                xcnt = rdenmax(1,icnt)
                ycnt = rdenmax(2,icnt)
                zcnt = rdenmax(3,icnt)
                write(12,'(3(1x,e12.5))') xcnt, ycnt, zcnt
                x0 = rdenmax(1,icnt) + x0 * vmod
                y0 = rdenmax(2,icnt) + y0 * vmod
                z0 = rdenmax(3,icnt) + z0 * vmod
            endif
            nlineas = nlineas + 1
            call metgrad(x0, y0, z0, numpnt, npntot)
            if (lcompl) nlincomp = nlincomp + 1
        enddo
!         reads the starting points coordinates for extra lines from an external file:
!              icnt: number of center from which the line departs
!                    if (0 < icnt <= ncen) the line departs from center (x(icnt),y(icnt),z(icnt))
!                                       in the direction of (x(icnt),y(icnt),z(icnt)) + (x0,y0,z0)
!                    if (icnt < 1 .or. icnt > ncen) the line departs from point x0, y0, z0
        lexist = .false.
        inquire(file=TRIM(filelines), exist=lexist, iostat=ierr)
        if (ierr .eq. 0 .and. lexist) then
            open(7, file=filelines, status='unknown', form='formatted', iostat=ioerr)
            if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(filelines)//'.Stop')
            write(6,"('Reads extra lines data from file ',a)") trim(filelines)
            ioerr = 0
            do while (ioerr .eq. 0)
                read(7,*,iostat=ioerr) icnt, x0, y0, z0
                if (ioerr .eq. -1) exit ! End of file
                if (ioerr .gt. 0) call error(1,'Error when reading file '//trim(filelines)//'.Stop')
                if (icnt .lt. 1 .or. icnt .gt. ncen) then
                    xcnt = x0
                    ycnt = y0
                    zcnt = z0
                    write(12,'(3(1x,e12.5))') xcnt, ycnt, zcnt
                else
                    vmod = sqrt(x0*x0+y0*y0+z0*z0)
                    if (vmod .eq. 0.d0) cycle
                    vmod = dlt0 / vmod
                    xcnt = rdenmax(1,icnt)
                    ycnt = rdenmax(2,icnt)
                    zcnt = rdenmax(3,icnt)
                    write(12,'(3(1x,e12.5))') xcnt, ycnt, zcnt
                    x0 = rdenmax(1,icnt) + x0 * vmod
                    y0 = rdenmax(2,icnt) + y0 * vmod
                    z0 = rdenmax(3,icnt) + z0 * vmod
                endif
                nlineas = nlineas + 1
                call metgrad(x0, y0, z0, numpnt, npntot)
                if (lcompl) nlincomp = nlincomp + 1
            enddo
        endif
        tiempo = dtime(tarray)
        write(6,"(1x,'Timing in seconds of individual points tabulation in proc 0 (user, system, total):', &
                5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tarray(1)+tarray(2)
    endif
    close(12)
    CALL MPI_REDUCE(nlineas,nlineastot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(npntot,npuntot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(nlincomp,nlincomptot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    if (myrank .eq. 0) then
        call system("cat "//trim(filename)//".dengr_?? > "//trim(filename)//".dengr")
        call system("rm -f "//trim(filename)//".dengr_??")
        write(6,"(/7x,'STATISTICS',/)")
        write(6,"('Number of computed points = ', i9)") npuntot
        write(6,"('Total number of attempted lines = ', i9)") nlineastot
        write(6,"('Total number of completed lines = ', i9)") nlincomptot     
    endif
    call MPI_FINALIZE(ierr)
    stop
    end

!   ***************************************************************

  subroutine metgrad(x0, y0, z0, numpnt, npntot)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    implicit none
    integer(KINT) :: i, ik, irebound, j, k, numpnt, npntot
    real(KREAL) :: aux, den, dlt, dltjmp, dspmodi, ex, exfor, ey, eyfor, ez, ezfor
    real(KREAL) :: gxnext, gynext, gznext, xnext, ynext, znext, exnext, eynext, eznext
    real(KREAL) :: gmod, gx, gxfor, gy, gyfor, gz, gzfor, prdesc, r13sq
    real(KREAL) :: x, x0, xant, xdsp, xfor, y, y0, yant, ydsp, yfor
    real(KREAL) :: z, z0, zant, zdsp, zfor
    real(KREAL), parameter :: ud3 = 0.333333333333333d0
    x = x0
    y = y0
    z = z0
!    Minus density gradient in point (x,y,z)
    call densgrad(x, y, z, den, ex, ey, ez)
    gx = ex
    gy = ey
    gz = ez
    write(12,"(3(1x,e12.5))") x, y, z
!     remaining points along the density gradient line
    lcompl = .true.
    irebound = 0
    dlt = dlt0
    do j = 1, numpnt
!          determines jump 
        gmod = sqrt(gx*gx + gy*gy + gz*gz)
        if (gmod .lt. 1.d-10 .or. dlt .lt. 1.d-2*dlt0) exit  ! If critical point or too short jump exits
        dltjmp = dlt / gmod
        xfor = x + dltjmp * gx
        yfor = y + dltjmp * gy
        zfor = z + dltjmp * gz
    
!         Minus density gradient in point (xfor,yfor,zfor)
        call densgrad(xfor, yfor, zfor, den, exfor, eyfor, ezfor)
        gxfor = exfor
        gyfor = eyfor
        gzfor = ezfor
        prdesc = (gxfor*gx + gyfor*gy + gzfor*gz) / sqrt(((gx*gx)+(gy*gy)+(gz*gz)) &
            * ((gxfor*gxfor)+(gyfor*gyfor)+(gzfor*gzfor)))
        if(xfor .lt. xinf .or. xfor .gt. xsup .or. yfor .lt. yinf .or. yfor .gt. ysup .or. &
            zfor .lt. zinf .or. zfor .gt. zsup) exit
        if (prdesc .lt. cero) then 
            irebound = irebound + 1
            if (irebound .gt. 2) then
                exit ! Three consecutive rebounds: Exits for rebound of the line
            else   ! Checks for a strong change of direction
                gmod = sqrt(gxfor*gxfor + gyfor*gyfor + gzfor*gzfor)
                if (gmod .lt. 1.d-10 .or. dlt .lt. 1.d-2*dlt0) exit  ! If critical point or too short jump exits
                dltjmp = dlt / gmod
                xnext = xfor + dltjmp * gxfor
                ynext = yfor + dltjmp * gyfor
                znext = zfor + dltjmp * gzfor
!               Minus density gradient in point (xnext,ynext,znext)
                call densgrad(xnext, ynext, znext, den, exnext, eynext, eznext)
                gxnext = exnext
                gynext = eynext
                gznext = eznext
                prdesc = (gxfor*gxnext + gyfor*gynext + gzfor*gznext) / sqrt(((gxnext*gxnext)+(gynext*gynext)+(gznext*gznext)) &
                    * ((gxfor*gxfor)+(gyfor*gyfor)+(gzfor*gzfor)))
                if (prdesc .gt. cero) then
                    irebound = 0
                    write(12,"(3(2x,e12.5))") xfor, yfor, zfor
                    xfor = xnext
                    yfor = ynext
                    zfor = znext
                    gxfor = gxnext
                    gyfor = gynext
                    gzfor = gznext
                    write(12,"(3(2x,e12.5))") xfor, yfor, zfor
                    dlt = dlt0
                else     ! Halves the step size and tries again with the new point
                    dlt = dlt * umed
                    cycle
                endif
            endif
        else
            irebound = 0
            write(12,"(3(2x,e12.5))") xfor, yfor, zfor
        endif  
!          Updates the point
        x = xfor
        y = yfor
        z = zfor
        gx = gxfor
        gy = gyfor
        gz = gzfor
    enddo         ! End of Do over j
    if (j .gt. numpnt) then
        write(6,"('Highest number of steps in line reached')")
        lcompl = .false.
    endif
    npntot = npntot + j
    write(12,*) ' '
    return
    end
!
!    ***************************************************************
!
  subroutine leedamqtden
    USE MPI
    USE DAM320_D
    USE DAMDENGRAD320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, ncenbas, ncfaj, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, r, ra, ral, step, suml, summ, t
    real(KREAL) :: tcheb(0:mxlenpol-1)
    inquire(file=trim(projectname)//"_2016.damqt", size=nsize, iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error when inquiring file ', a, ' in processor ',i3)") trim(projectname)//"_2016.damqt", myrank
        abort = 1
        return
    endif
    if (nsize .eq. -1) then
        write(6,"('Size of file ', a, ' cannot be determined in processor ',i3)") trim(projectname)//"_2016.damqt", myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of file ', a, ' = ', i12)") trim(projectname)//".damqt", nsize
#if _WIN32
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#elif __INTEL_COMPILER
     open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#else
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='unformatted', action = 'read', access='stream', iostat=ierr)
#endif
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//"_2016.damqt", myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Opens file ', a)") trim(projectname)//"_2016.damqt"
    read(10) ncen, nbas, ncaps
    nsize = nsize - sizeof(ncen) - sizeof(nbas) - sizeof(ncaps)
    if (myrank .eq. 0) write(6,"('ncen = ', i4, ' nbas = ', i6, ' ncaps = ', i5)") ncen, nbas, ncaps
     
!    Allocates memory for geometry

    allocate(atmnam(ncen), nzn(ncen), rcen(3,ncen), zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam, nzn, rcen and zn in processor ',i3)") myrank
        abort = 1
        return
    endif

    if (myrank .eq. 0) write(6,"(/24x,'GEOMETRY (BOHR)')")
    if (myrank .eq. 0) write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge')")
!    Geometry and nuclear charges
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

    read(10) lsto  ! .true. means STO basis, .false. means GTO basis
! write(6,"('sizeof(lsto) = ', i3)") sizeof(lsto)
    nsize = nsize - sizeof(lsto)
    if (lsto) then
     
!         Allocates memory for the basis set

        allocate(ll(ncaps), lmaxc(ncen), nf(ncaps), ngini(ncen), ngfin(ncen), nn(ncaps), rlargo(ncen), rnor(ncaps), &
                xx(ncaps), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating ll, lmaxc, nf, ngini, ngfin, nn, rlargo, rnor and xx in processor ',i3)") &
            myrank
            abort = 1
            return
        endif
          
        if (myrank .eq. 0) write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        if (myrank .eq. 0 .and. longoutput) write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp')")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfin(ia)
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfin(ia))
            lmaxc(ia) = 0
            rlargo(ia) = cero
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            if (myrank .eq. 0 .and. longoutput) write(6,"(t5,'center ', i4,/t12,'n', t16, 'l', t25,'exp', t35, 'ind func')") ia
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (myrank .eq. 0 .and. longoutput) write(6,"(t11,i2,t15,i2,t20,e12.5,t36,i4)") nn(i), ll(i), xx(i), nf(i)
            enddo
        enddo
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)
          
!         Allocates memory for the basis set

        allocate(cfcontr(nprimitot), ipntprim(ncaps), ll(ncaps), lmaxc(ncen), ncontr(ncen), nf(ncaps), ngini(ncen), &
                ngfin(ncen), nprimit(ncaps), rlargo(ncen), rnor(ncaps), xxg(nprimitot), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating cfcontr, ipntprim, ll, lmaxc, ncontr, nf, ngini, ngfin, &
                    &nprimit, rlargo, rnor and xxg in processor ',i3)") myrank
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
            rlargo(ia) = cero
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
!                   computes and stores the radial normalization factor
                aux = cero
                bux = ll(knt) + 1.5d0
                do k1 = 1, nprimit(knt)
                    do k2 = 1, k1-1
                        aux=aux + dos*cfcontr(icarga+k1)*cfcontr(icarga+k2)/(xxg(icarga+k1)+xxg(icarga+k2))**bux
                    enddo
                    aux = aux + cfcontr(icarga+k1) * cfcontr(icarga+k1) / (dos*xxg(icarga+k1))**bux
                enddo
                icarga = icarga+nprimit(knt)  ! actualizes the index for loading primitives exponents and contraction coefficients
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
        if (myrank .eq. 0) write(6,"('Total number of primitives = ', i8)") nprimitot
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
    
    if (myrank .eq. 0 .and. longoutput) write(6,"('lmaxexp = ', i2, ' nintervaj = ', i2)") lmaxexp, nintervaj
    
    allocate(icfposd(lmtop*nintervaj+1,ncen), xajustd(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating icfpos and xajust in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    nsize = nsize - sizeof(icfposd(:,1)) * ncenbas
    if (myrank .eq. 0 .and. longoutput) write(6,"('Estimated highest size of xajust   = ', i15, ' bytes')") size(xajustd)
    nsize = nsize - sizeof(xajustd(:,1)) * ncenbas
    
    allocate(cfajust(nsize/8), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfajust in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of cfajust   = ', i15, ' bytes')") size(cfajust)

    if (myrank .eq. 0 .and. longoutput) write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
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
        write(6,"('Memory error when allocating indintrv in processor ',i3)") myrank
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
    close(10)
     
!    Determines the long-range radii and the highest l in the expansion for each interval

    allocate(lcorto(nintervaj,ncen), umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lcorto and umedpow in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)
    umedpow(0) = uno                                  !
    do i = 1, lmaxexp                                 !
        umedpow(i) = umedpow(i-1) * umed             ! 1 / 2^i
    enddo
    if (myrank .eq. 0) write(6,"('Long-range threshold = ',e12.5)") umbrlargo
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
            do i = 0, nsamples-1     ! samples over nsamples points in each interval to decide the highest l 
                ra = rinterv(interv-1) + dltsample + (rinterv(interv) - rinterv(interv-1) - dos * dltsample) &
                        * ri(nsamples-1) * i
                aux = exp(-xajustd(interv,ia)*ra)
                t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
                dost = t + t
                tcheb(0) = uno ! Chebyshev T  polynomials
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
        if (myrank .eq. 0) then
            if (longoutput) then
                write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5, ' lcorto = ', 30(i3))") &
                    ia, atmnms(nzn(ia)), rlargo(ia), lcorto(1:nintervaj,ia)
            else
                write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5)") &
                    ia, atmnms(nzn(ia)), rlargo(ia)
            endif
        endif
    enddo
    deallocate(umedpow)
    return
    end
!**********************************************************************
!    subroutine consta
!
!     Computes and stores auxiliary constants
!          re(i) = dfloat(i)
!          ri(i) = 1.d0 / dfloat(i)
!          fact(i) = dfloat(i!)
!          facti(i) = 1.d0 / dfloat(i!)
!          facts(i) = dfloat((i+1/2)!)
!          ind(i) = i*(i+1)/2
!         ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) 
!               / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    implicit none
    integer(KINT) :: i
!     auxiliary parameters and functions
    pi = acos(-uno)
    raizpi = sqrt(pi)
    re(0) = cero
    ri(0) = 1.d300
    dosl1(0) = uno
    dosl1i(0) = uno
    do i = 1, mxreal
        re(i) = re(i-1) + uno        ! dfloat(i)
        re(-i) = -re(i)
        ri(i) = uno / re(i)            ! uno / dfloat(i)
        ri(-i) = -ri(i)
        dosl1(i) = re(i) + re(i) + uno     ! dfloat(iì)
        dosl1(-i) = -re(i) - re(i) + uno
        dosl1i(i) = uno / dosl1(i)          ! dfloat( 1/(i+i+1) )
        dosl1i(-i) = -dosl1i(i)
    enddo
    fact(0) = uno
    facti(0) = uno
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)             !  i!
        facti(i) = uno / fact(i)               !  uno / i!
    enddo
    return
    end
     
!   ***************************************************************

  subroutine densgrad(x, y, z, den, dendrvx, dendrvy, dendrvz)
    USE DAM320_D
    USE DAMDENGRAD320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, icflm, interv, j, jshft, kntlm, l, lmax, m
    real(KREAL) :: aux, bux, cux, den, dendrvx, dendrvy, dendrvz, dost, dux, drvflm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, eux, flm, fux
    real(KREAL) :: r3inv, ra, ra2, rainv, rj2, sgn, t, umt2i, x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL) :: tcheb(0:mxlenpol-1), ucheb(0:mxlenpol-1), drvtcheb(0:mxlenpol-1)
    den = cero
    dendrvx = cero
    dendrvy = cero
    dendrvz = cero
    do ia = 1, ncen
!       Contribution of atomic fragment ia to density  in point (x,y,z)
    if (ngini(ia) .le. 0) cycle
    xa = x - rcen(1,ia)
    ya = y - rcen(2,ia)
    za = z - rcen(3,ia)
    ra2 = xa*xa+ya*ya+za*za
    if (ra2 .gt. rlargo(ia)*rlargo(ia)) cycle
    ra = sqrt(ra2)
    if (ra .lt. rlargo(ia)) then
        interv = indintrv(int(fct*ra)+1)
    else
        interv = nintervaj
    endif
    lmax = lcorto(interv,ia)
    if (ra .eq. cero) then
        xadivra = uno
        xadivra = uno
        zadivra = uno
        lmax = 0
    else
        rainv = uno / ra
        xadivra = xa * rainv
        yadivra = ya * rainv
        zadivra = za * rainv
    endif

    zlma(1) = uno       ! Regular spherical harmonics of r-R(ia)
    zlma(2) = ya
    zlma(3) = za
    zlma(4) = xa
    do l = 1, lmax
        zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))      ! zlm(l+1,l+1,ia)
        zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))      ! zlm(l+1,-(l+1),ia)
        zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)                  ! zlm(l+1,l,ia)
        zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)                         ! zlm(l+1,-l,ia)
        do m = 0, l-1
            zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1)) ! zlm(l+1,m,ia)
            zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1)) ! zlm(l+1,-m,ia)
        enddo
    enddo
    call derivzlm(lmax, idimzlm, zlma, zlmadx, zlmady, zlmadz)
    aux = exp(-xajustd(interv,ia)*ra)
    bux = xajustd(interv,ia)
    cux = dos / (rinterv(interv)-rinterv(interv-1))
    t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
    dost = t + t
    tcheb(0) = uno      ! Chebyshev T and U polynomials (U polynomial are used for derivatives of T polynomials
    tcheb(1) = t        ! and first and second derivatives of T polynomials according to:
    ucheb(0) = uno      !              D[T_n(t),t] = n * U_(n-1)(t)
    ucheb(1) = dost     !    and:      D[T_n(t),{t,2}] = (n/(t^2-1)) * ( (n-1)*t*U_(n-1) - n * U_(n-2) )
    drvtcheb(0) = cero
    drvtcheb(1) = uno
    sgn = uno
    do j = 2, mxlenpol-1
        tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
        ucheb(j) = dost * ucheb(j-1) - ucheb(j-2)
        drvtcheb(j) = re(j) * ucheb(j-1)
    enddo
    kntlm = 0
        do l = 0, lmax  !     Computes density terms 0 <= l <= lcorto(interv,ia)
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
                        drvflm = cux * aux * drvflm -bux * flm      ! Converts the derivative with respect to t to derivative 
                                                                    ! respect to r  D[flm,r]
                        den = den + flm * zlma(kntlm)
                        dendrvx = dendrvx + xadivra * drvflm * zlma(kntlm) + flm * zlmadx(kntlm)
                        dendrvy = dendrvy + yadivra * drvflm * zlma(kntlm) + flm * zlmady(kntlm)
                        dendrvz = dendrvz + zadivra * drvflm * zlma(kntlm) + flm * zlmadz(kntlm)
                    endif
            enddo
        enddo
    enddo
!      Changes the sign to gradient
    dendrvx = -dendrvx
    dendrvy = -dendrvy
    dendrvz = -dendrvz
    return
    end
!
!   ***************************************************************************************
!     Calculates the density derivatives from the represented density at point (x,y,z)
!
  subroutine densdrvs(x, y, z, dendrvx, dendrvy, dendrvz, dxx, dxy, dxz, dyy, dyz, dzz)
    USE DAM320_D
    USE DAMDENGRAD320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, icflm, interv, j, jshft, kntlm, l, m
    real(KREAL) :: aux, bux, cux, dendrvx, dendrvy, dendrvz, dost, dux, drvflm, drv2flm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, eux, flm, fux
    real(KREAL) :: r3inv, ra, ra2, rainv, rj2, sgn, t, umt2i, x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL) :: tcheb(0:mxlenpol-1), ucheb(0:mxlenpol-1), drvtcheb(0:mxlenpol-1), drv2tcheb(0:mxlenpol-1)
    dendrvx = cero
    dendrvy = cero
    dendrvz = cero
    dxx = cero
    dxy = cero
    dxz = cero
    dyy = cero
    dyz = cero
    dzz = cero
!    Computes the derivatives of the electron density
    do ia = 1, ncen
        xa = x - rcen(1,ia)
        ya = y - rcen(2,ia)
        za = z - rcen(3,ia)
        ra2 = xa*xa+ya*ya+za*za
        if (ra2 .gt. rlargo(ia)*rlargo(ia)) cycle
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
                ! zlm(l+1,m,ia)
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1)) 
                ! zlm(l+1,-m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))   
            enddo
        enddo
        call derivzlm(lcorto(interv,ia), idimzlm, zlma, zlmadx, zlmady, zlmadz)
        call derivzlm(lcorto(interv,ia), idimzlm, zlmadx, zlmadxx, zlmadxy, zlmadxz)
        call dzlm2y(lcorto(interv,ia), idimzlm, zlmady, zlmadyy, zlmadyz)
        call dzlm2z(lcorto(interv,ia), idimzlm, zlmadz, zlmadzz)
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
        drv2tcheb(0) = cero
        drv2tcheb(1) = cero
        if (uno-abs(t) .gt. 1.d-7) then
            umt2i = uno / (t*t - uno)
            do j = 2, mxlenpol-1
                drv2tcheb(j) = re(j) * umt2i * (re(j-1) * t * ucheb(j-1) - re(j) * ucheb(j-2))
            enddo
        else    ! For values of t very close to 1 or -1, takes a linear approximation (Taylor series) of D[T[j,t],{t,2}]
            do j = 2, mxlenpol-1
                rj2 = re(j) * re(j)
                drv2tcheb(j) = sgn * ri(3) * (rj2 * (rj2-uno) + ri(5) * rj2 * (re(4)+rj2*(-re(5)+rj2)) * (abs(t)-uno) ) 
                if (t .lt. cero) sgn = - sgn
            enddo
        endif
        kntlm = 0
        do l = 0, lcorto(interv,ia)    !     Computes density 
            do m = -l, l
                kntlm = kntlm + 1
                icflm = icfposd(kntlm+(interv-1)*lmtop,ia)
                if(icflm .lt. icfposd(kntlm+(interv-1)*lmtop+1,ia)) then
                    flm = cero
                    drvflm = cero
                    drv2flm = cero
                    do j = 0, icfposd(kntlm+(interv-1)*lmtop+1,ia)-icflm-1
                        flm = flm + cfajust(j+icflm) * tcheb(j)
                        drvflm = drvflm + cfajust(j+icflm) * drvtcheb(j)
                        drv2flm = drv2flm + cfajust(j+icflm) * drv2tcheb(j)
                    enddo
                    flm = aux * flm
                    drvflm = cux * aux * drvflm        ! Converts the derivative with respect to t to derivative respect to r
                    drv2flm = cux * cux * aux * drv2flm    ! Same for second derivative
                    drv2flm = bux * (bux * flm - (drvflm+drvflm)) + drv2flm    ! D[flm,{r,2}]
                    drvflm = -bux * flm + drvflm                            ! D[flm,r]
                    dendrvx = dendrvx + xadivra * drvflm * zlma(kntlm) + flm * zlmadx(kntlm)
                    dendrvy = dendrvy + yadivra * drvflm * zlma(kntlm) + flm * zlmady(kntlm)
                    dendrvz = dendrvz + zadivra * drvflm * zlma(kntlm) + flm * zlmadz(kntlm)
                    r3inv = rainv * rainv * rainv
                    fux = r3inv * (ra * drv2flm - drvflm)
                    dxx = dxx + (xa*xa * fux + drvflm * rainv) * zlma(kntlm) &
                        + dos * xadivra * drvflm * zlmadx(kntlm) + flm * zlmadxx(kntlm)
                    dyy = dyy + (ya*ya * fux + drvflm * rainv) * zlma(kntlm) &
                        + dos * yadivra * drvflm * zlmady(kntlm) + flm * zlmadyy(kntlm)
                    dzz = dzz + (za*za * fux + drvflm * rainv) * zlma(kntlm) &
                        + dos * zadivra * drvflm * zlmadz(kntlm) + flm * zlmadzz(kntlm)
                    dxy = dxy + xa*ya * fux * zlma(kntlm) + xadivra * drvflm * zlmady(kntlm) &
                        + yadivra * drvflm * zlmadx(kntlm) + flm * zlmadxy(kntlm)
                    dxz = dxz + xa*za * fux * zlma(kntlm) + xadivra * drvflm * zlmadz(kntlm) &
                        + zadivra * drvflm * zlmadx(kntlm) + flm * zlmadxz(kntlm)
                    dyz = dyz + ya*za * fux * zlma(kntlm) + yadivra * drvflm * zlmadz(kntlm) &
                        + zadivra * drvflm * zlmady(kntlm) + flm * zlmadyz(kntlm)
                endif
            enddo
        enddo
    enddo
!    Changes the sign of the derivatives of the density
    dendrvx = -dendrvx
    dendrvy = -dendrvy
    dendrvz = -dendrvz
    dxx = -dxx
    dxy = -dxy
    dxz = -dxz
    dyy = -dyy
    dyz = -dyz
    dzz = -dzz
    return
    end

!   ***************************************************************

  subroutine seekdenmax     ! Fine tunes the density maxima in the neighborhood of hydrogens
    USE MPI
    USE DAM320_D
    USE DAMDENGRAD320_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: icnt, ierr
    real(KREAL) :: x0, xf, y0, yf, z0, zf
    allocate(rdenmax(3,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rdenmax in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0) write(6,"(/'Positions of density local maxima')")
    do icnt = 1, ncen
        if (nzn(icnt) .eq. 1) then
            x0 = rcen(1,icnt)
            y0 = rcen(2,icnt)
            z0 = rcen(3,icnt)
            call findmax(x0, y0, z0, xf, yf, zf)
            rdenmax(1,icnt) = xf
            rdenmax(2,icnt) = yf
            rdenmax(3,icnt) = zf
        else
            rdenmax(1,icnt) = rcen(1,icnt)
            rdenmax(2,icnt) = rcen(2,icnt)
            rdenmax(3,icnt) = rcen(3,icnt)
        endif
        if (myrank .eq. 0) write(6,"('Center ', i5, 2x,' R: ',3(1x,f12.7))") icnt, rdenmax(1:3,icnt)
    enddo
    if (myrank .eq. 0) write(6,"(/)")
    return
    end

!   ***************************************************************

  subroutine findmax(x0, y0, z0, x, y, z)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    implicit none
    integer(KINT) :: i, ik, irebound, j, k, numpnt, npntot
    real(KREAL) :: aux, den, denfor, dennext, denmax, dlt, dltjmp, dspmodi, ex, exfor, ey, eyfor, ez, ezfor
    real(KREAL) :: gxnext, gynext, gznext, xnext, ynext, znext, exnext, eynext, eznext
    real(KREAL) :: gmod, gx, gxfor, gy, gyfor, gz, gzfor, prdesc, r13sq
    real(KREAL) :: x, x0, xant, xdsp, xfor, y, y0, yant, ydsp, yfor
    real(KREAL) :: z, z0, zant, zdsp, zfor
    real(KREAL), parameter :: ud3 = 0.333333333333333d0
    x = x0
    y = y0
    z = z0
!    Minus density gradient in point (x,y,z)
    call densgrad(x, y, z, den, ex, ey, ez)
    denmax = den
    gx = ex
    gy = ey
    gz = ez
!     remaining points along the density gradient line
    irebound = 0
    dlt = -1.d-3
    numpnt = 100
    do j = 1, numpnt
!          determines jump
        gmod = sqrt(gx*gx + gy*gy + gz*gz)
        if (gmod .lt. 1.d-10 .or. abs(dlt) .lt. 1.d-6) exit  ! If critical point or too short jump exits
        dltjmp = dlt / gmod
        xfor = x + dltjmp * gx
        yfor = y + dltjmp * gy
        zfor = z + dltjmp * gz
!         Minus density gradient in point (xfor,yfor,zfor)
        call densgrad(xfor, yfor, zfor, denfor, exfor, eyfor, ezfor)
        gxfor = exfor
        gyfor = eyfor
        gzfor = ezfor
        prdesc = (gxfor*gx + gyfor*gy + gzfor*gz) / sqrt(((gx*gx)+(gy*gy)+(gz*gz)) &
            * ((gxfor*gxfor)+(gyfor*gyfor)+(gzfor*gzfor)))
        if (prdesc .lt. cero) then
            irebound = irebound + 1
            if (irebound .gt. 2) then
                exit ! Three consecutive rebounds: Exits for rebound of the line
            else   ! Checks for a strong change of direction
                gmod = sqrt(gxfor*gxfor + gyfor*gyfor + gzfor*gzfor)
                if (gmod .lt. 1.d-10 .or. abs(dlt) .lt. 1.d-6) exit  ! If critical point or too short jump exits
                dltjmp = dlt / gmod
                xnext = xfor + dltjmp * gxfor
                ynext = yfor + dltjmp * gyfor
                znext = zfor + dltjmp * gzfor
!               Minus density gradient in point (xnext,ynext,znext)
                call densgrad(xnext, ynext, znext, dennext, exnext, eynext, eznext)
                gxnext = exnext
                gynext = eynext
                gznext = eznext
                prdesc = (gxfor*gxnext + gyfor*gynext + gzfor*gznext) / sqrt(((gxnext*gxnext)+(gynext*gynext)+(gznext*gznext)) &
                    * ((gxfor*gxfor)+(gyfor*gyfor)+(gzfor*gzfor)))
                if (prdesc .gt. cero) then
                    irebound = 0
                    xfor = xnext
                    yfor = ynext
                    zfor = znext
                    gxfor = gxnext
                    gyfor = gynext
                    gzfor = gznext
                    dlt = dlt0
                else     ! Halves the step size and tries again with the new point
                    dlt = dlt * umed
                    cycle
                endif
            endif
        else
            irebound = 0
        endif
!          Updates the point provided the density has increased
        if (denfor .lt. denmax) then
            exit
        else
            denmax = denfor
        endif
        x = xfor
        y = yfor
        z = zfor
        gx = gxfor
        gy = gyfor
        gz = gzfor
    enddo         ! End of Do over j
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
     zlmadx(5) = re(6) * zlma(2)   ! Derivatives of the D spherical harmonics
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
     zlmadx(10) = re(15) * zlma(5)      ! Derivatives of the F spherical harmonics
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
     do l = 4, lmax      ! Derivatives of the remaining spherical harmonics
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
    do l = 4, lmax       ! Derivatives of the remaining spherical harmonics
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
    do l = 1, lmax       ! Derivatives of the remaining spherical harmonics
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
!     -------------------------------------------------------------------------------------------------------
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
