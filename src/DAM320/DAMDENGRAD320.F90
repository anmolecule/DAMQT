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
! Program for computing density gradient lines from the representation of the molecular density performed with
! DAM320
!
! IMPORTANT!!!
! The 2D grid calculation is intended only for symmetry planes. The symmetry plane (passing through the origin)
!   must be supplied in terms of the three parameters (A, B, C) according to the plane definition:
!
!           A x + B y + C z = 0
!
! These parameters are supplied in the namelist variables  planeA, planeB, planeC
!
! The following cases are considered (controlled by planecase):
!
!                                           (A, B, C)
!                                           ---------
!       planecase = 1       plane: XY       (0, 0, C)    (z = 0)
!       planecase = 2       plane: XZ       (0, B, 0)    (y = 0)
!       planecase = 3       plane: YZ       (A, 0, 0)    (x = 0)
!       planecase = 4                       (A, B, 0)
!       planecase = 5                       (A, 0, C)
!       planecase = 6                       (0, B, C)
!       planecase = 7                       (A, B, C)
!
! Version of September 2018
!
  program DAMDENGRAD320
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    USE icosahedron
    implicit none
    character(1) :: cpslabel
    integer(KINT) :: i, ia, icnt, ierr, ioerr, ioplines3D, ipnt, itheta, j, k, ncps, ncntplane, nlinpernuc, npntot, numpnt
    integer(KINT) :: naux, nbux
    real(KREAL) :: aux, R, u0, ucnt, v0, vcnt, usalto
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz
    real(KREAL) :: dxdu, dxdv, dydu, dydv, dzdu, dzdv, ex, ey, ez, gu, gv, x0, xcnt, vmod, y0, ycnt, z0, zcnt
    real(KREAL) :: amat(2,3), bmat(2,3), eigval(2), eigvec(2,2), hxyz(3,3), huv(2,2)
    real(KREAL), allocatable :: rcntplane(:,:)
    real(KREAL), parameter :: angs2bohr = 1.889725989d0
    real(KREAL4) :: tarray(2), tiempo, dtime
    character(256) :: filelines
    logical :: lnucleo, lexist
    namelist / options / dlt0, filename, filelines, icntlines, ioplines3D, iswindows, lextralines, lmaxrep, longoutput &
            , lplot2d, nlines, nlinpernuc, numpnt, planeA, planeB, planeC &
            , thresh, umbrlargo, usalto, uvratio, rlines &
            , vmod, xinf, xsup, yinf, ysup, zinf, zsup, uinf, usup, vinf, vsup
!     Namelist default values
    ioplines3D = 1      ! Set of lines per nucleus for 3D plots based on icosahedron vertices, C2 axes and C3 axes or combinations of them:
                        ! 1: vertices (12 points);   2: C3 axes (20 points);   3: C2 axes (30 points)
                        ! 4: vertices + c3 ;   5: vertices + C2;   6: C2 + C3;   7: vertices + C2 + C3 
    longoutput = .false.! If .true. a more detailed output is given
    iswindows = .false.          ! .true. if running on a MS-windows system
    lmaxrep = 5          ! highest value of  l  in the expansion of the density for computing the density gradient
    umbrlargo = 1.d-8     ! Threshold for determining the short-range radius
    usalto = 1.d-3          ! Threshold for convergence in jump
    numpnt = 2000          ! Maximum number of points in each density gradient line
    lextralines = .false. ! If .true. reads extra lines from a file
    lplot2d = .false.       ! If .true. 2D plot (lines projection over a 2D surface)
    nlinpernuc = 16         ! Number of lines per nucleus for 2D plots
    filename = ""               ! root file name for .dengr files
    filelines = ""          ! File with starting points for density gradient lines
    xinf = cero
    xsup = cero
    yinf = cero
    ysup = cero
    zinf = cero
    zsup = cero
    dlt0 = 1.d-2
    planeA = cero        ! Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    planeB = cero
    planeC = uno
    uinf = cero
    usup = uno
    uvratio = uno
    vinf = cero
    vsup = uno
    nlines = 0          ! Number of starting points for extra lines in namelist
    icntlines = 0       ! Array with indices for starting nuclei of extra lines (0 for lines starting in a point other than nuclei)
    rlines = cero       ! Array with coordinates (in au) of points defining lines
!     End of namelist defaults
    tiempo = dtime(tarray)

    read(5,OPTIONS)
    read(5,*) projectname

    ioplines3D = max(0,min(7,ioplines3D))    ! Forces ioplines3D to be in the range [0,7]

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
    call consta          !     Computes auxiliary constants
    call leedamqtden     !     Reads file .damqt  (generated by DAM2016)
    if (lmaxrep .gt. lmaxexp) then
            write(6,"('lmaxrep = ', i3, ' greater than lmaxexp ', i3)") lmaxrep, lmaxexp
            write(6,"('takes lmaxrep = ',i3)") lmaxexp
            lmaxrep = lmaxexp
    endif
    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zlma. Stop')
    allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating zlmadx, zlmady, zlmadz. Stop')

    call seekdenmax     ! Checks if density maxima are in nuclei, otherwise seeks the actual positions of the maxima

!   Checks that the plane is a plane of symmetry and if it is, identifies the case
    planecase = 1
    if (lplot2d) then
        call checkplanecase
        write(6,"('Plane: A = ', e17.10, ' B = ', e17.10, ' C = ', e17.10)") planeA, planeB, planeC
        write(6,"('Case = ', i3)") planecase
        write(6,"('wu = ', 3(2x,e17.10))") wu
        write(6,"('wv = ', 3(2x,e17.10))") wv
    endif

    write(6,"('Density gradient lines from expansion of the density: lmaxrep = ', i3)") lmaxrep
    if (lplot2d) then
        if (usup-uinf .le. cero) then
            call error(100,'uinf must be lower than usup')
        endif
        write(6,"('uinf = ', e17.10, ' usup = ', e17.10)") uinf, usup
        if (vsup-vinf .le. cero) then
            call error(100,'vinf must be lower than vsup')
        endif
        write(6,"('vinf = ', e17.10, ' vsup = ', e17.10)") vinf, vsup
    else
        if (xsup-xinf .le. cero) then
            call error(100,'xinf must be lower than xsup')
        endif
        write(6,"('xinf = ', e17.10, ' xsup = ', e17.10)") xinf, xsup
        if (ysup-yinf .le. cero) then
            call error(100,'yinf must be lower than ysup')
        endif
        write(6,"('yinf = ', e17.10, ' ysup = ', e17.10)") yinf, ysup
        if (zsup-zinf .le. cero) then
            call error(100,'zinf must be lower than zsup')
        endif
        write(6,"('zinf = ', e17.10, ' zsup = ', e17.10)") zinf, zsup
    endif
!     Opens file for density gradient lines
    if (lplot2d) then
        open(12,file=trim(filename)//".dengr2D",status='unknown',form='formatted',iostat=ioerr)
        if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(filename)//".dengr2D"//'.Stop')
    else
        open(12,file=trim(filename)//".dengr",status='unknown',form='formatted',iostat=ioerr)
        if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(filename)//".dengr"//'.Stop')
    endif

!    density gradient lines
    npntot = 0
    nlineas = 0
    nlincomp = 0
    if (.not. lplot2d) then     ! ********** 3D lines **********
        r = dlt0
        do ia = 1, ncen
            xcnt = rdenmax(1,ia)
            ycnt = rdenmax(2,ia)
            zcnt = rdenmax(3,ia)
            if (ioplines3D .eq. 1 .or. ioplines3D .eq. 4 .or. ioplines3D .eq. 5 .or. ioplines3D .eq. 7) then
                do i = 1, idimvert
                    x0 = xcnt + r * vertices(i,1)
                    y0 = ycnt + r * vertices(i,2)
                    z0 = zcnt + r * vertices(i,3)
                    nlineas = nlineas + 1
                    write(12,"(3(1x,e12.5))") xcnt, ycnt, zcnt
                    call metgrad(x0, y0, z0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            endif
            if (ioplines3D .eq. 2 .or. ioplines3D .eq. 4 .or. ioplines3D .eq. 6 .or. ioplines3D .eq. 7) then
                do i = 1, idimc3
                    x0 = xcnt + r * c3axes(i,1)
                    y0 = ycnt + r * c3axes(i,2)
                    z0 = zcnt + r * c3axes(i,3)
                    nlineas = nlineas + 1
                    write(12,"(3(1x,e12.5))") xcnt, ycnt, zcnt
                    call metgrad(x0, y0, z0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            endif
            if (ioplines3D .eq. 3 .or. ioplines3D .eq. 5 .or. ioplines3D .eq. 6 .or. ioplines3D .eq. 7) then
                do i = 1, idimc2
                    x0 = xcnt + r * c2axes(i,1)
                    y0 = ycnt + r * c2axes(i,2)
                    z0 = zcnt + r * c2axes(i,3)
                    nlineas = nlineas + 1
                    write(12,"(3(1x,e12.5))") xcnt, ycnt, zcnt
                    call metgrad(x0, y0, z0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            endif
        enddo
        if (lextralines) then
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
        endif
        write(6,"(/7x,'STATISTICS',/)")
        write(6,"('Number of computed points = ', i9)") npntot
        write(6,"('Total number of attempted lines = ', i9)") nlineas
        write(6,"('Total number of completed lines = ', i9)") nlincomp
    else        ! ********** 2D lines **********
        allocate(rcntplane(6,ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating rcntplane. Stop')
        if (nlinpernuc .lt. 1) then ! If no lines per nucleus, writes the corners coordinates
            write(12,"(2(2x,e12.5),/2(2x,e12.5),/)") uinf, vinf, uinf, vinf
            write(12,"(2(2x,e12.5),/2(2x,e12.5),/)") uinf, vinf, uinf, vinf
            write(12,"(2(2x,e12.5),/2(2x,e12.5),/)") uinf, vsup, uinf, vsup
            write(12,"(2(2x,e12.5),/2(2x,e12.5),/)") usup, vinf, usup, vinf
            write(12,"(2(2x,e12.5),/2(2x,e12.5),/)") usup, vsup, usup, vsup
        endif
        ncntplane = 0
        r = dlt0
        do ia = 1, ncen
            xcnt = rcen(1,ia)
            ycnt = rcen(2,ia)
            zcnt = rcen(3,ia)
            if (abs(planeA * xcnt + planeB * ycnt + planeC * zcnt) .gt. 1.e-6) cycle
            ncntplane = ncntplane + 1
            if (planecase .eq. 1) then
                ucnt = xcnt
                vcnt = ycnt
            else if (planecase .eq. 2) then
                ucnt = xcnt
                vcnt = zcnt
            else if (planecase .eq. 3) then
                ucnt = ycnt
                vcnt = zcnt
            else if (planecase .eq. 4) then
                ucnt = xcnt / wu(1)
                vcnt = zcnt
            else if (planecase .eq. 5) then
                ucnt = ycnt
                vcnt = xcnt / wv(1)
            else if (planecase .eq. 6) then
                ucnt = -xcnt
                vcnt = ycnt / wv(2)
            else
                ucnt = wu(1) * xcnt + wu(2) * ycnt
                vcnt = (wv(1) * xcnt + wv(2) * ycnt) / (uno - wv(3)*wv(3))
            endif
            rcntplane(1,ncntplane) = ucnt
            rcntplane(2,ncntplane) = vcnt
            rcntplane(3,ncntplane) = zn(ia)
            rcntplane(4,ncntplane) = xcnt
            rcntplane(5,ncntplane) = ycnt
            rcntplane(6,ncntplane) = zcnt
            xcnt = rdenmax(1,ia)
            ycnt = rdenmax(2,ia)
            zcnt = rdenmax(3,ia)
            if (planecase .eq. 1) then
                ucnt = xcnt
                vcnt = ycnt
            else if (planecase .eq. 2) then
                ucnt = xcnt
                vcnt = zcnt
            else if (planecase .eq. 3) then
                ucnt = ycnt
                vcnt = zcnt
            else if (planecase .eq. 4) then
                ucnt = xcnt / wu(1)
                vcnt = zcnt
            else if (planecase .eq. 5) then
                ucnt = ycnt
                vcnt = xcnt / wv(1)
            else if (planecase .eq. 6) then
                ucnt = -xcnt
                vcnt = ycnt / wv(2)
            else
                ucnt = wu(1) * xcnt + wu(2) * ycnt
                vcnt = (wv(1) * xcnt + wv(2) * ycnt) / (uno - wv(3)*wv(3))
            endif

            do itheta = 0, nlinpernuc-1
                u0 = ucnt + r * cos(itheta * dos * pi / nlinpernuc)
                v0 = vcnt + uvratio * r * sin(itheta * dos * pi / nlinpernuc)
                nlineas = nlineas + 1
                write(12,"(2(2x,e12.5))") ucnt, vcnt
                call metgrad2d(u0, v0, numpnt, npntot)
                if (lcompl) nlincomp = nlincomp + 1
            enddo
        enddo
        if (lextralines) then
!       reads the starting points coordinates for extra lines from namelist
            if (nlines .gt. 0) write(6,"('Reads extra lines data from namelist')")
            do i = 1, nlines
                icnt = icntlines(i)
                dltu = rlines(1,i)
                dltv = rlines(2,i)
                if (icnt .lt. 1 .or. icnt .gt. ncen) then
                    u0 = dltu
                    v0 = dltv
                else
                    xcnt = rdenmax(1,icnt)
                    ycnt = rdenmax(2,icnt)
                    zcnt = rdenmax(3,icnt)
                    if (abs(planeA * xcnt + planeB * ycnt + planeC * zcnt) .gt. 1.e-6) cycle
                    if (planecase .eq. 1) then
                        ucnt = xcnt
                        vcnt = ycnt
                    else if (planecase .eq. 2) then
                        ucnt = xcnt
                        vcnt = zcnt
                    else if (planecase .eq. 3) then
                        ucnt = ycnt
                        vcnt = zcnt
                    else if (planecase .eq. 4) then
                        ucnt = xcnt / wu(1)
                        vcnt = zcnt
                    else if (planecase .eq. 5) then
                        ucnt = ycnt
                        vcnt = xcnt / wv(1)
                    else if (planecase .eq. 6) then
                        ucnt = -xcnt
                        vcnt = ycnt / wv(2)
                    else
                        ucnt = wu(1) * xcnt + wu(2) * ycnt
                        vcnt = (wv(1) * xcnt + wv(2) * ycnt) / (uno - wv(3)*wv(3))
                    endif
                    write(12,"(2(2x,e12.5))") ucnt, vcnt
                    vmod = sqrt(dltu*dltu+dltv*dltv)
                    if (vmod .eq. 0.d0) cycle
                    vmod = dlt0 / vmod
                    u0 = ucnt + vmod * dltu
                    v0 = vcnt + vmod * dltv
                endif
                nlineas = nlineas + 1
                call metgrad2d(u0, v0, numpnt, npntot)
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
                    read(7,*,iostat=ioerr) icnt, dltu, dltv
                    if (ioerr .eq. -1) exit ! End of file
                    if (ioerr .gt. 0) call error(1,'Error when reading file '//trim(filelines)//'.Stop')
                    if (icnt .lt. 1 .or. icnt .gt. ncen) then
                        u0 = dltu
                        v0 = dltv
                    else
                        xcnt = rdenmax(1,icnt)
                        ycnt = rdenmax(2,icnt)
                        zcnt = rdenmax(3,icnt)
                        if (abs(planeA * xcnt + planeB * ycnt + planeC * zcnt) .gt. 1.e-6) cycle
                        if (planecase .eq. 1) then
                            ucnt = xcnt
                            vcnt = ycnt
                        else if (planecase .eq. 2) then
                            ucnt = xcnt
                            vcnt = zcnt
                        else if (planecase .eq. 3) then
                            ucnt = ycnt
                            vcnt = zcnt
                        else if (planecase .eq. 4) then
                            ucnt = xcnt / wu(1)
                            vcnt = zcnt
                        else if (planecase .eq. 5) then
                            ucnt = ycnt
                            vcnt = xcnt / wv(1)
                        else if (planecase .eq. 6) then
                            ucnt = -xcnt
                            vcnt = ycnt / wv(2)
                        else
                            ucnt = wu(1) * xcnt + wu(2) * ycnt
                            vcnt = (wv(1) * xcnt + wv(2) * ycnt) / (uno - wv(3)*wv(3))
                        endif
                        write(12,"(2(2x,e12.5))") ucnt, vcnt
                        vmod = sqrt(dltu*dltu+dltv*dltv)
                        if (vmod .eq. 0.d0) cycle
                        vmod = dlt0 / vmod
                        u0 = ucnt + vmod * dltu
                        v0 = vcnt + vmod * dltv
                    endif
                    nlineas = nlineas + 1
                    call metgrad2d(u0, v0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
            endif
        endif
        naux = nlineas
        nbux = nlincomp
        nlineas = 0
        nlincomp = 0
!        2D Basins
        write(12,"(2x,i6)") ncntplane
        do ia = 1, ncntplane
            write(12,"(6(2x,e12.5))") rcntplane(1:6,ia)
        enddo
        write(12,"(3(2x,e12.5))") planeA, planeB, planeC
        close(12)
        inquire(file=trim(projectname)//"-cps-d.xyz", exist=lcpsd, iostat=ierr)
        if (lcpsd) then 
            allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                    zlmadzz(idimzlm), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating zlmadxx, zlmadxy,... zlmadzz. Stop')
            allocate (ra2l1((lmaxrep+1)**2), ra2l1inv((lmaxrep+1)**2), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating ra2l1 and ra2l1inv. Stop')
            open(17, file=trim(projectname)//"-cps-d.xyz", status='unknown', form='formatted', iostat=ierr)
            if (ierr .eq. 0) then
                open(12,file=trim(filename)//"-d.basin2D",status='unknown',form='formatted',iostat=ioerr)
                if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(filename)//".cam2D"//'.Stop')
                read(17,*) ncps
                do ipnt = 1, ncps

                    read(17,*) cpslabel, x0, y0, z0
                    x0 = x0 * angs2bohr
                    y0 = y0 * angs2bohr
                    z0 = z0 * angs2bohr
                    if (cpslabel .ne. 'z' .or. abs(planeA * x0 + planeB * y0 + planeC * z0) .gt. 1.e-6) cycle
                    if (planecase .eq. 1) then
                        ucnt = x0
                        vcnt = y0
                    else if (planecase .eq. 2) then
                        ucnt = x0
                        vcnt = z0
                    else if (planecase .eq. 3) then
                        ucnt = y0
                        vcnt = z0
                    else if (planecase .eq. 4) then
                        ucnt = x0 / wu(1)
                        vcnt = z0
                    else if (planecase .eq. 5) then
                        ucnt = y0
                        vcnt = x0 / wv(1)
                    else if (planecase .eq. 6) then
                        ucnt = -x0
                        vcnt = y0 / wv(2)
                    else
                        ucnt = wu(1) * x0 + wu(2) * y0
                        vcnt = (wv(1) * x0 + wv(2) * y0) / (uno - wv(3)*wv(3))
                    endif
                    dxdu = wu(1)
                    dxdv = wv(1)
                    dydu = wu(2)
                    dydv = wv(2)
                    dzdu = wu(3)
                    dzdv = wv(3)
                    call densdrvs(x0, y0, z0, ex, ey, ez, dxx, dxy, dxz, dyy, dyz, dzz)
                    gu = ex * dxdu + ey * dydu + ez * dzdu
                    gv = ex * dxdv + ey * dydv + ez * dzdv
                    amat(1,1) = dxdu
                    amat(1,2) = dydu
                    amat(1,3) = dzdu
                    amat(2,1) = dxdv
                    amat(2,2) = dydv
                    amat(2,3) = dzdv
                    hxyz(1,1) = dxx
                    hxyz(1,2) = dxy
                    hxyz(1,3) = dxz
                    hxyz(2,1) = dxy
                    hxyz(2,2) = dyy
                    hxyz(2,3) = dyz
                    hxyz(3,1) = dxz
                    hxyz(3,2) = dyz
                    hxyz(3,3) = dzz
                    do i = 1, 2
                        do j = 1, 3
                            aux = 0.d0
                            do k = 1, 3
                                aux = aux + amat(i,k) * hxyz(k,j)
                            enddo
                            bmat(i,j) = aux
                        enddo
                    enddo
                    do i = 1, 2
                        do j = 1, 2
                            aux = 0.d0
                            do k = 1, 3
                                aux = aux + bmat(i,k) * amat(j,k)
                            enddo
                            huv(i,j) = aux
                        enddo
                    enddo
                    aux = sqrt(huv(1,1)*huv(1,1) + cuatro * huv(1,2)*huv(1,2) - dos * huv(1,1)*huv(2,2) + huv(2,2)*huv(2,2))
                    eigval(1) = umed * (huv(1,1) + huv(2,2) - aux)
                    eigval(2) = umed * (huv(1,1) + huv(2,2) + aux)
                    if (eigval(1)* eigval(2) .gt. cero) then
                        write(6,"('Eigenvalues have both the same sign')")
                        cycle
                    endif
                    eigvec(1,1) = huv(1,1) - huv(2,2) - aux
                    eigvec(1,2) = dos * huv(1,2)
                    eigvec(2,1) = huv(1,1) - huv(2,2) + aux
                    eigvec(2,2) = dos * huv(1,2)
                    aux = sqrt(eigvec(1,1)*eigvec(1,1)+eigvec(1,2)*eigvec(1,2))
                    eigvec(1,1) = eigvec(1,1) / aux
                    eigvec(1,2) = eigvec(1,2) / aux
                    aux = sqrt(eigvec(2,1)*eigvec(2,1)+eigvec(2,2)*eigvec(2,2))
                    eigvec(2,1) = eigvec(2,1) / aux
                    eigvec(2,2) = eigvec(2,2) / aux
                    if (eigval(1) .gt. cero) then
                        dltu = eigvec(1,1) * dlt0
                        dltv = eigvec(1,2) * dlt0
                    else
                        dltu = eigvec(2,1) * dlt0
                        dltv = eigvec(2,2) * dlt0
                    endif

                    u0 = ucnt + dltu
                    v0 = vcnt + dltv
                    write(12,"(2(2x,e12.5))") ucnt, vcnt
                    call metgrad2d(u0, v0, numpnt, npntot)

                    u0 = ucnt - dltu
                    v0 = vcnt - dltv
                    write(12,"(2(2x,e12.5))") ucnt, vcnt
                    nlineas = nlineas + 1
                    call metgrad2d(u0, v0, numpnt, npntot)
                    if (lcompl) nlincomp = nlincomp + 1
                enddo
                close(12)
            endif
        endif
        write(6,"(/7x,'STATISTICS',/)")
        write(6,"('Number of computed points = ', i9)") npntot
        write(6,"('Total number of attempted lines = ', i9)") naux
        write(6,"('Total number of completed lines = ', i9)") nbux
        write(6,"('Total number of attempted lines in basins = ', i9)") nlineas
        write(6,"('Total number of completed lines in basins = ', i9)") nlincomp
    endif

    tiempo = dtime(tarray)
    write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tiempo
    stop
    end


!   ***************************************************************

  subroutine checkplanecase
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    implicit none
    logical :: found
    integer(KINT) :: ia, ib
    real(KREAL) :: aux, ABC(3), rdif(3), rtransf(3)
    ABC(1) = planeA
    ABC(2) = planeB
    ABC(3) = planeC
    aux = dos / dot_product(ABC,ABC)
    do ia = 1, ncen
        rtransf(1:3) = rcen(1:3,ia) - aux * dot_product(ABC,rcen(1:3,ia)) * ABC
        found = .false.
        do ib = 1, ncen
            rdif(1:3) = rcen(1:3,ib) - rtransf
            if (dot_product(rdif,rdif) .lt. 1.e-8) then
                found = .true.
                exit
            endif
        enddo
        if (.not. found) then
            write(6,"('The plane: A = ', e12.5, ' B = ', e12.5, ' C = ', e12.5, ' is not a symmetry plane.')") &
                planeA, planeB, planeC
            write(6,"('Reflection of center ', 3(1x,e17.10), ' gives ', 3(1x,e17.10), ' which is not a center')") &
                rcen(1:3,ia), rtransf
            call error(2,'Stop')
        endif
    enddo

    if (abs(planeA) .lt. 1.e-7) then
        if (abs(planeB) .lt. 1.e-7) then
            if (abs(planeC) .lt. 1.e-7) then
                write(6,"('Error in plane parameters: A = ', e12.5, ' B = ', e12.5, ' C = ', e12.5)") planeA, planeB, planeC
                stop
            else
                planecase = 1       ! A == 0, B == 0, C != 0
                wu = (/ uno, cero, cero /)
                wv = (/ cero, uno, cero /)
            endif
        else
            if (abs(planeC) .lt. 1.e-7) then
                planecase = 2       ! A == 0, B != 0, C == 0
                wu = (/ uno, cero, cero /)
                wv = (/ cero, cero, uno /)
            else
                planecase = 6       ! A == 0, B != 0, C != 0
                wu = (/ -uno, cero, cero /)
                wv = (/ cero, -ABC(3), ABC(2) /) / sqrt(ABC(2)*ABC(2)+ABC(3)*ABC(3))
            endif
        endif
    else
        if (abs(planeB) .lt. 1.e-7) then
            if (abs(planeC) .lt. 1.e-7) then
                planecase = 3       ! A != 0, B == 0, C == 0
                wu = (/ cero, uno, cero /)
                wv = (/ cero, cero, uno /)
            else
                planecase = 5       ! A != 0, B == 0, C != 0   
                wu = (/ cero, uno, cero /)
                wv = (/ -ABC(3), cero, ABC(1) /) / sqrt(ABC(1)*ABC(1)+ABC(3)*ABC(3))
            endif
        else
            if (abs(planeC) .lt. 1.e-7) then
                planecase = 4       ! A != 0, B != 0, C == 0
                wu = (/ -ABC(2), ABC(1), cero /) / sqrt(ABC(1)*ABC(1)+ABC(2)*ABC(2))
                wv = (/ cero, cero, uno /)
            else
                planecase = 7       ! A != 0, b != 0, C != 0
                wu = (/ -ABC(2), ABC(1), cero /) / sqrt(ABC(1)*ABC(1)+ABC(2)*ABC(2))
                wv = (/ -ABC(1)*ABC(3), -ABC(2)*ABC(3), ABC(1)*ABC(1)+ABC(2)*ABC(2) /) &
                    / sqrt((ABC(1)*ABC(1)+ABC(2)*ABC(2))*(ABC(1)*ABC(1)+ABC(2)*ABC(2)+ABC(3)*ABC(3)))
            endif
        endif
    endif
    return
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


!   ***************************************************************

  subroutine metgrad2d(u0, v0, numpnt, npntot)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDENGRAD320_D
    implicit none
    integer(KINT) :: i, ik, irebound, j, k, numpnt, npntot
    real(KREAL) :: aux, den, dlt, dltjmp, dspmodi, dxdu, dxdv, dydu, dydv, dzdu, dzdv, ex, exfor, ey, eyfor, ez, ezfor
    real(KREAL) :: gxnext, gynext, gznext, gunext, gvnext, xnext, ynext, znext, exnext, eynext, eznext, unext, vnext
    real(KREAL) :: gmod, gu, gufor, gv, gvfor, prdesc, r13sq
    real(KREAL) :: u, u0, uant, udsp, ufor, v, v0, vant, vdsp, vfor
    real(KREAL) :: x, x0, xant, xfor, y, y0, yant, yfor, z, z0, zant, zfor
    real(KREAL), parameter :: ud3 = 0.333333333333333d0

    dxdu = wu(1)
    dxdv = wv(1)
    dydu = wu(2)
    dydv = wv(2)
    dzdu = wu(3)
    dzdv = wv(3)
! write(6,"('diff_x_u = ', e12.5, ' diff_x_v = ', e12.5)") dxdu, dxdv
! write(6,"('diff_y_u = ', e12.5, ' diff_y_v = ', e12.5)") dydu, dydv
! write(6,"('diff_z_u = ', e12.5, ' diff_z_v = ', e12.5)") dzdu, dzdv
    u = u0
    v = v0
    x = u * wu(1) + v * wv(1)
    y = u * wu(2) + v * wv(2)
    z = u * wu(3) + v * wv(3)
! write(6,"(//'u0 = ', e12.5, ' v0 = ', e12.5, ' x = ', e12.5, ' y = ', e12.5, ' z = ', e12.5)") u, v, x, y, z

!   Minus density gradient in point (x,y,z)
    call densgrad(x, y, z, den, ex, ey, ez)

! write(6,"('ex = ', e12.5, ' ey = ', e12.5, ' ez = ', e12.5)") ex, ey, ez
    gu = ex * dxdu + ey * dydu + ez * dzdu
    gv = ex * dxdv + ey * dydv + ez * dzdv
! write(6,"('gu = ', e12.5, ' gv = ', e12.5)") gu, gv
    write(12,"(2(2x,e12.5))") u0, v0
!     remaining points along the density gradient line
    lcompl = .true.
    irebound = 0
    dlt = dlt0
    do j = 1, numpnt
!          determines jump
        gmod = sqrt(gu*gu + gv*gv)
        if (gmod .lt. 1.d-10 .or. dlt .lt. 1.d-2*dlt0) exit  ! If critical point or too short jump exits
        dltjmp = dlt / gmod
        ufor = u + dltjmp * gu
        vfor = v + dltjmp * gv
        xfor = ufor * wu(1) + vfor * wv(1)
        yfor = ufor * wu(2) + vfor * wv(2)
        zfor = ufor * wu(3) + vfor * wv(3)
! write(6,"('ufor = ', e12.5, ' vfor = ', e12.5,' xfor = ', e12.5,' yfor = ', e12.5,' zfor = ', e12.5)") ufor, vfor, xfor, yfor, zfor
!       Minus density gradient in point (xfor,yfor,zfor)
        call densgrad(xfor, yfor, zfor, den, exfor, eyfor, ezfor)
! write(6,"('exfor = ', e12.5, ' eyfor = ', e12.5, ' ezfor = ', e12.5)") exfor, eyfor, ezfor
        gufor = exfor * dxdu + eyfor * dydu + ezfor * dzdu
        gvfor = exfor * dxdv + eyfor * dydv + ezfor * dzdv
! write(6,"('gufor = ', e12.5, ' gvfor = ', e12.5)") gufor, gvfor
        prdesc = (gufor*gu + gvfor*gv) / sqrt(((gu*gu)+(gv*gv)) * ((gufor*gufor)+(gvfor*gvfor)))
        if(ufor .lt. uinf .or. ufor .gt. usup .or. vfor .lt. vinf .or. vfor .gt. vsup) exit
        if (prdesc .lt. cero) then
            irebound = irebound + 1
            if (irebound .gt. 2) then
                exit ! Three consecutive rebounds: Exits for rebound of the line
            else   ! Checks for a strong change of direction
                gmod = sqrt(gufor*gufor + gvfor*gvfor)
                if (gmod .lt. 1.d-10 .or. dlt .lt. 1.d-2*dlt0) exit  ! If critical point or too short jump exits
                dltjmp = dlt / gmod
                unext = ufor + dltjmp * gufor
                vnext = vfor + dltjmp * gvfor
                xnext = unext * wu(1) + vnext * wv(1)
                ynext = unext * wu(2) + vnext * wv(2)
                znext = unext * wu(3) + vnext * wv(3)
                call densgrad(xnext, ynext, znext, den, exnext, eynext, eznext)
                gunext = exnext * dxdu + eynext * dydu + eznext * dzdu
                gvnext = exnext * dxdv + eynext * dydv + eznext * dzdv
                prdesc = (gunext*gufor + gvnext*gvfor) / sqrt(((gufor*gufor)+(gvfor*gvfor)) * ((gunext*gunext)+(gvnext*gvnext)))
                if (prdesc .gt. cero) then
                    irebound = 0
                    write(12,"(2(2x,e12.5))") ufor, vfor
!write(6,"('ufor, vfor = ', 2(2x,e12.5))") ufor, vfor
                    ufor = unext
                    vfor = vnext
                    gufor = gunext
                    gvfor = gvnext
                    write(12,"(2(2x,e12.5))") ufor, vfor
!write(6,"('ufor, vfor = ', 2(2x,e12.5))") ufor, vfor
                    dlt = dlt0
                else     ! Halves the step size and tries again with the new point
!write(6,"('Halves the step size and tries again with the new point')")
                    dlt = dlt * umed
                    cycle
                endif
            endif
        else
            irebound = 0
            write(12,"(2(2x,e12.5))") ufor, vfor
        endif
!          Updates the point
        u = ufor
        v = vfor
        gu = gufor
        gv = gvfor
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
    USE DAM320_D
    USE DAMDENGRAD320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, ncenbas, ncfaj, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, r, ra, ral, step, suml, summ, t
    real(KREAL) :: tcheb(0:mxlenpol-1) 
    inquire(file=trim(projectname)//"_2016.damqt", size=nsize, iostat=ierr)
    if (ierr .ne. 0) call error(ierr,'Error when inquiring file '//trim(projectname)//"_2016.damqt")
    if (nsize .eq. -1) call error(1,'Size of file '//trim(projectname)//"_2016.damqt cannot be determined")
    if (longoutput) write(6,"('Size of file ', a, ' = ', i12)") trim(projectname)//".damqt", nsize
#if _WIN32
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#elif __INTEL_COMPILER
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#else
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='unformatted', action = 'read', access='stream', iostat=ierr)
#endif
    if (ierr .ne. 0) call error(ierr,'Error when opening file '//trim(projectname)//"_2016.damqt")
    if (longoutput) write(6,"('Opens file ', a)") trim(projectname)//"_2016.damqt"
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
!    Geometry and nuclear charges
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

    read(10) lsto  ! .true. means STO basis, .false. means GTO basis
! write(6,"('sizeof(lsto) = ', i3)") sizeof(lsto)
    nsize = nsize - sizeof(lsto)
    if (lsto) then
    
!         Allocates memory for the basis set
        
        allocate(ll(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating ll. Stop')
        
        allocate(lmaxc(ncen), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating lmaxc. Stop')
        
        allocate(nf(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating nf. Stop')
        
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating ngini. Stop')
        
        allocate(ngfin(ncen), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating ngfin. Stop')
        
        allocate(nn(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating nn. Stop')
        
        allocate(rlargo(ncen), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating rlargo. Stop')
        
        allocate(xx(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating xx. Stop')
        
        write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        if (longoutput) write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp')")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfin(ia)
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfin(ia))
            lmaxc(ia) = 0
            rlargo(ia) = cero
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            if (longoutput) write(6,"(t5,'center ', i4,/t12,'n', t16, 'l', t25,'exp', t35, 'ind func')") ia
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (longoutput) write(6,"(t11,i2,t15,i2,t20,e12.5,t36,i4)") nn(i), ll(i), xx(i), nf(i)
            enddo
        enddo
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)
        
!         Allocates memory for the basis set
        
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
        
        allocate(rlargo(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating rlargo. Stop')
        
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
! write(6,"('icarga = ', i12)") icarga
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
        write(6,"('Total number of primitives = ', i8)") nprimitot
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
! write(6,"('size(icfposd) = ', i12)") sizeof(icfposd)
    nsize = nsize - sizeof(icfposd(:,1)) * ncenbas
        
    allocate(xajustd(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating xajustd. Stop')
    if (longoutput) write(6,"('Estimated highest size of xajust   = ', i15, ' bytes')") size(xajustd) 
    nsize = nsize - sizeof(xajustd(:,1)) * ncenbas

! write(6,"('Tamaño estimado de cfajust = ', i12)") nsize/8
! call flush(6)
     
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
    close(10)
     
!    Determines the long-range radii and the highest l in the expansion for each interval

    allocate(lcorto(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating lcorto. Stop')
    if (longoutput) write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)

    allocate(umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating umedpow in leedamqtden. Stop')
    umedpow(0) = uno                                  !
    do i = 1, lmaxexp                                 !
        umedpow(i) = umedpow(i-1) * umed             ! 1 / 2^i
    enddo
    write(6,"('Long-range threshold = ',e12.5)") umbrlargo
    nsamples = 4
    do ia = 1, ncen      ! Do over centers
! write(6,"('---- centro = ', i3)") ia
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
! write(6,"('sample = ', i3, ' aux = ', e17.10)") i, aux
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
        if (longoutput) then
            write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5, ' lcorto = ', 30(i3))") &
                ia, atmnms(nzn(ia)), rlargo(ia), lcorto(1:nintervaj,ia)
        else 
            write(6,"('Long-range radius for center ',i4,' (',a2,') = ', e12.5)") &
                ia, atmnms(nzn(ia)), rlargo(ia)
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
                        drvflm = cux * aux * drvflm -bux * flm        ! Converts the derivative with respect to t to derivative respect to r
                                                                        ! D[flm,r]
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
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))    ! zlm(l+1,m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
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
     USE DAM320_D
     USE DAMDENGRAD320_D
     USE DAM320_DATA_D
     implicit none
     integer(KINT) :: icnt, ierr
     real(KREAL) :: x0, xf, y0, yf, z0, zf
     allocate(rdenmax(3,ncen), stat = ierr)
     if (ierr .ne. 0) call error(ierr,'Memory error when allocating denmax. Stop')
     write(6,"(/'Positions of density local maxima')")
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
        write(6,"('Center ', i5, 2x,' R: ',3(1x,f12.7))") icnt, rdenmax(1:3,icnt)
     enddo
     write(6,"(/)")
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
!write(6,*) 'den ini = ', den
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
