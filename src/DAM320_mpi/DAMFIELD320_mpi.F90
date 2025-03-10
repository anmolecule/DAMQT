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
! Program for computing electic field lines from the representation of the molecular density performed with
! DAM320. Version parallelized with MPI only 3D lines.
!
! Version of April 2019
!
  program DAMFIELD320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMFIELD320_D
    USE ICOSAHEDRON
    USE icosahedron
    USE PARALELO
    implicit none
    character(1) :: cpslabel
    integer(KINT) :: i, ia, icnt, ierr, ioerr, ioplines3D, ipnt, itheta, j, knt, kntcortotot, kntinf, kntlargotot, kntsup 
    integer(KINT) :: naux, nlincomptot, nlineastot, ncntplane, nlinpernuc, npntot, npuntot, numpnt
    real(KREAL) :: aux, R, u0, ucnt, v0, vcnt, usalto
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz
    real(KREAL) :: dxdu, dxdv, dydu, dydv, dzdu, dzdv, ex, ey, ez, gu, gv, x0, xcnt, vmod, y0, ycnt, z0, zcnt
    real(KREAL) :: amat(2,3), bmat(2,3), eigval(2), eigvec(2,2), hxyz(3,3), huv(2,2)
    real(KREAL), allocatable :: rcntplane(:,:)
    real(KREAL), parameter :: angs2bohr = 1.889725989d0
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    logical :: lnamelist(2), ltimeprocs
    integer(KINT) :: inamelist(3)
    real(KREAL) :: rnamelist(9)
    character(256) :: filelines
    logical :: lnucleo, lexist
    character*4 :: strbux
    namelist / options / basintol, dlt0, filename, filelines, icntlines, ioplines3D, iswindows, largo, lextralines &
            , lmaxrep, longoutput, lplot2d, lvalence, nlines, nlinpernuc, numpnt, planeA, planeB, planeC &
            , thresh, umbrlargo, usalto, uvratio, rlines &
            , vmod, xinf, xsup, yinf, ysup, zinf, zsup, uinf, usup, vinf, vsup
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    abort = 0
    abortroot = 0
    if (myrank .eq. 0) write(6,"('number of processors = ', i3)") nprocs
    tiempo = dtime(tarray)
!	Namelist default values
    basintol = 1.e-6    ! Not used. Kept for compatibility of input files
    ioplines3D = 1      ! Set of lines per nucleus for 3D plots based on icosahedron vertices, C2 axes and C3 axes or combinations of them:
                        ! 1: vertices (12 points);   2: C3 axes (20 points);   3: C2 axes (30 points)
                        ! 4: vertices + c3 ;   5: vertices + C2;   6: C2 + C3;   7: vertices + C2 + C3
    longoutput = .false.! If .true. a more detailed output is given
    lvalence = .false.  ! If .true. only valence electrons are considered
    iswindows = .false. ! .true. if running on a MS-windows system
    lmaxrep = 5         ! highest value of  l  in the expansion of the density for computing the field
    umbrlargo = 1.d-8   ! Threshold for determining the short-range radius
    usalto = 1.d-3      ! Threshold for convergence in jump
    numpnt = 2000       ! Maximum number of points in each field line
    lextralines = .false.	! If .true. reads extra lines from a file
    largo = .false.     ! If .true. computes the long-range field
    lplot2d = .false.   ! Not used. Kept for compatibility of input files
    nlinpernuc = 16     ! Kept for compatibility of input files
    filename = ""       ! root file name for .cam files
    filelines = ""      ! File with starting points for field lines
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
    icntlines = 0       ! Array with indices for starting nuclei of extra lines (0 for lines starting in a point other than nuclei)
    rlines = cero       ! Array with coordinates (in au) of points defining lines
!	End of namelist defaults
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
        
        ioplines3D = max(0,min(7,ioplines3D))    ! Forces ioplines3D to be in the range [1,7]
        inamelist = (/ lmaxrep, ioplines3D, numpnt /)
        write(6,"(1x,'project name : ',a,/,1x,'==============')") projectname
        if (iswindows) then
                dirsep = "\\"
                i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
                if (i .eq. 0) then	! This is intended for MinGW, whose directory separator in windows is also /
                        dirsep = "/"
                        i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
                endif
        else
                dirsep = "/"
                i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
        end if

        if (lvalence) then
            write(6,"(//'Using valence electrons only',//)")
        endif

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
        lnamelist = (/ ltimeprocs, lvalence /)
    endif
    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lnamelist,2,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (myrank .ne. 0) then
        lmaxrep = inamelist(1); ioplines3D = inamelist(2); numpnt = inamelist(3); 
        ltimeprocs = lnamelist(1); 
        lvalence = lnamelist(2);
        xinf = rnamelist(1); xsup = rnamelist(2); 
        yinf = rnamelist(3); ysup = rnamelist(4); 
        zinf = rnamelist(5); zsup = rnamelist(6); dlt0 = rnamelist(7)
        umbrlargo = rnamelist(8); usalto = rnamelist(9)
    endif
        
    call consta		!	Computes auxiliary constants
    call leedamqtfield	!	Reads file .damqt  (generated by DAM2016)
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
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma in processor ',i3)") myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    if (myrank .eq. 0) write(6,"('Electric field from expansion of the density: lmaxrep = ', i3)") lmaxrep

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
    
!	Opens file for field lines
    
    write(strbux,'(i4.2)') myrank
    open (unit=12, file=trim(filename)//".cam_"//trim(adjustl(strbux)),status='unknown',form='formatted',iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Error in processor ',i3, ' when opening file ', a, '.Stop')") myrank, &
            trim(filename)//".cam_"//trim(adjustl(strbux))
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    
!    3D Force lines

    npntot = 0
    nlineas = 0
    nlincomp = 0
    kntlargo = 0
    kntcorto = 0
    
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
    if (myrank .eq. 0 .and. longoutput) write(6,*) 'naux = ', naux, ' myrank = ', myrank, ' kntinf = ', kntinf, &
        ' kntsup = ', kntsup

    r = dlt0
    do ia = 1, ncen
        xcnt = rcen(1,ia)
        ycnt = rcen(2,ia)
        zcnt = rcen(3,ia)
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
                xcnt = rcen(1,icnt)
                ycnt = rcen(2,icnt)
                zcnt = rcen(3,icnt)
                write(12,'(3(1x,e12.5))') xcnt, ycnt, zcnt
                x0 = rcen(1,icnt) + x0 * vmod
                y0 = rcen(2,icnt) + y0 * vmod
                z0 = rcen(3,icnt) + z0 * vmod
            endif
            nlineas = nlineas + 1
            call metgrad(x0, y0, z0, numpnt, npntot)
            if (lcompl) nlincomp = nlincomp + 1
        enddo
!		reads the starting points coordinates for extra lines from an external file:
!			icnt: number of center from which the line departs
!				 if (0 < icnt <= ncen) the line departs from center (x(icnt),y(icnt),z(icnt))
!								in the direction of (x(icnt),y(icnt),z(icnt)) + (x0,y0,z0)
!				 if (icnt < 1 .or. icnt > ncen) the line departs from point x0, y0, z0
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
                    xcnt = rcen(1,icnt)
                    ycnt = rcen(2,icnt)
                    zcnt = rcen(3,icnt)
                    write(12,'(3(1x,e12.5))') xcnt, ycnt, zcnt
                    x0 = rcen(1,icnt) + x0 * vmod
                    y0 = rcen(2,icnt) + y0 * vmod
                    z0 = rcen(3,icnt) + z0 * vmod
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
    CALL MPI_REDUCE(kntlargo,kntlargotot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(kntcorto,kntcortotot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (myrank .eq. 0) then
        call system("cat "//trim(filename)//".cam_?? > "//trim(filename)//".cam")
        call system("rm -f "//trim(filename)//".cam_??")
        write(6,"(/7x,'STATISTICS',/)")
        write(6,"('Number of computed points = ', i9)") npuntot
        write(6,"('Number of long-range contributions = ', i9)") kntlargotot
        write(6,"('Number of short-range contributions = ', i9)") kntcortotot
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
    USE DAMFIELD320_D
    implicit none
    integer(KINT) :: i, ik, irebound, j, k, numpnt, npntot
    real(KREAL) :: aux, dlt, dltjmp, dspmodi, ex, exfor, ey, eyfor, ez, ezfor
    real(KREAL) :: gxnext, gynext, gznext, xnext, ynext, znext, exnext, eynext, eznext
    real(KREAL) :: gmod, gx, gxfor, gy, gyfor, gz, gzfor, prdesc, r13sq
    real(KREAL) :: x, x0, xdsp, xfor, y, y0, ydsp, yfor
    real(KREAL) :: z, z0, zdsp, zfor
    real(KREAL), parameter :: ud3 = 0.333333333333333d0
    logical lnucleo
    x = x0
    y = y0
    z = z0
!     Electric field in point (x,y,z)
    if (largo) then
        call campolargo(x, y, z, ex, ey, ez, lnucleo)
    else
        call campo(x, y, z, ex, ey, ez, lnucleo)
    endif
    if (lnucleo) return
    gx = ex
    gy = ey
    gz = ez
    write(12,"(3(1x,e12.5))") x, y, z
!     remaining points along the field line
    lcompl = .true.
    irebound = 0
    dlt = dlt0
    do j = 1, numpnt
!		determines jump 
        gmod = sqrt(gx*gx + gy*gy + gz*gz)
        if (gmod .lt. 1.d-10 .or. dlt .lt. 1.d-2*dlt0) exit  ! If critical point or too short jump exits
        dltjmp = dlt / gmod
        xfor = x + dltjmp * gx
        yfor = y + dltjmp * gy
        zfor = z + dltjmp * gz
!		Electric field in point (xfor,yfor,zfor)
        if (largo) then
            call campolargo(xfor, yfor, zfor, exfor, eyfor, ezfor, lnucleo)
        else
            call campo(xfor, yfor, zfor, exfor, eyfor, ezfor, lnucleo)
        endif
        if (lnucleo) exit
        gxfor = exfor
        gyfor = eyfor
        gzfor = ezfor
        prdesc = (gxfor*gx + gyfor*gy + gzfor*gz) / sqrt(((gx*gx)+(gy*gy)+(gz*gz)) &
            * ((gxfor*gxfor)+(gyfor*gyfor)+(gzfor*gzfor)))
        if(xfor .lt. xinf .or. xfor .gt. xsup .or. yfor .lt. yinf .or. yfor .gt. ysup .or. &
                zfor .lt. zinf .or. zfor .gt. zsup .or. (largo .and. lnucleo)) exit
        if (prdesc .lt. cero) then
            irebound = irebound + 1
            if (irebound .gt. 2) then
                exit ! Three consecutive rebounds: Exits for rebound of the line
            else   ! Checks for a strong change of direction
                gmod = sqrt(gxfor*gxfor + gyfor*gyfor + gzfor*gzfor)
                dltjmp = dlt / gmod
                xnext = xfor + dltjmp * gxfor
                ynext = yfor + dltjmp * gyfor
                znext = zfor + dltjmp * gzfor
!           	Electric field in point (xnext,ynext,znext)
                if (largo) then
                    call campolargo(xnext, ynext, znext, exnext, eynext, eznext, lnucleo)
                else
                    call campo(xnext, ynext, znext, exnext, eynext, eznext, lnucleo)
                endif
                gxnext = exnext
                gynext = eynext
                gznext = eznext
                prdesc = (gxfor*gxnext + gyfor*gynext + gzfor*gznext) / sqrt(((gxnext*gxnext)+(gynext*gy)+(gznext*gznext)) &
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
!		Updates the point
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
    write(12,"(' ')")
    return
    end

!   ***************************************************************

  subroutine campo(x, y, z, ex, ey, ez, lnucleo)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    USE DAMFIELD320_D
    implicit none
    logical :: lnucleo
    integer(KINT) :: i, ia, icflm, ierr, interv, isgm, isgmm1, isgmp1, j, l, lm, lm1l, lp1lp2, ltop, m, mm
    real(KREAL) :: c1, c2, dost, Ex, Exext, Exnucl, Ey, Eyext, Eynucl, Ez, Ezext, Eznucl
    real(KREAL) :: pi4d2l1, ra, ra2, ra2inv, ra3inv, ra2l3inv, rinta, rintb, t, x, xa, y, ya, z, za
    real(KREAL) :: tcheb(0:mxlenpol-1)
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
!    Electric field in point (x,y,z)
    Exnucl = cero
    Eynucl = cero
    Eznucl = cero
    Exext = cero
    Eyext = cero
    Ezext = cero
    lnucleo = .false.
    do ia = 1, ncen
        xa = x - rcen(1,ia)
        ya = y - rcen(2,ia)
        za = z - rcen(3,ia)
        ra2 = xa*xa + ya*ya + za*za
        ra = sqrt(ra2)
        if (ra .lt. 1.d-10) then
            write(6,"('The point ',3(1x,F10.5),' coincides with center ',i4)") x, y, z, ia
            lnucleo = .true.
            return
        endif
        ra3inv = uno / (ra*ra2)
        Exnucl = Exnucl + zn(ia) * xa * ra3inv
        Eynucl = Eynucl + zn(ia) * ya * ra3inv
        Eznucl = Eznucl + zn(ia) * za * ra3inv
        if (ngini(ia) .le. 0) cycle
        ra2inv = uno / ra2
        if (ra .lt. rlargo(ia)) then
            interv = indintrv(int(fct*ra)+1)
            ltop = lcorto(interv,ia)
        else
            ltop = llargo(min(int(ra),mxlargo),ia)
        endif
        zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
        zlma(2) = ya
        zlma(3) = za
        zlma(4) = xa
        do l = 1, ltop
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))		! zlma(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))		! zlma(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)				! zlma(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)					! zlma(l+1,-l,ia)
            do m = 0, l-1
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlma(l+1,m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlma(l+1,-m,ia)
            enddo
        enddo

        if (ra .gt. rlargo(ia)) then   ! Long-range
            lm = 0
            ra2l3inv = uno / ra
            do l = 0, ltop
                ra2l3inv = ra2l3inv  * ra2inv	! 1.d0 / (ra**(2*l+3))
                do m = -l, l
                    lm = lm + 1
                    if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo) cycle
                    mm = iabs(m)
                    if (m .ge. 0) then
                            isgm = 1
                    else
                            isgm = -1
                    endif
!      Ex
                    c1 = zlma( (l+1)*(l+2)+isgm*(mm+1)+1)
                    if (m .ne. 0 .and. m .ne. -1) then
                            c1 = c1 - re(l+1-mm) * re(l+2-mm) * zlma( (l+1)*(l+2)+isgm*(mm-1)+1 )
                    else
                            if (m .eq. 0) then
                                    c1 = c1 + c1
                            endif
                    endif
                    c1 = c1 * ra2l3inv
                    Ex = -umed * rmultip(lm,ia) * c1
                    Exext = Exext + Ex
!      Ey
                    c1 = zlma( (l+1)*(l+2)-isgm*(mm+1)+1 )
                    if (m .ne. 0 .and. m .ne. 1) then
                            c1 = c1 + re(l+1-mm) * re(l+2-mm) * zlma( (l+1)*(l+2)-isgm*(mm-1)+1 )
                    else
                            if (m .eq. 0) then
                                    c1 = c1 + c1
                            endif
                    endif
                    c1 = c1 * ra2l3inv
                    Ey = -isgm * umed * rmultip(lm,ia) * c1
                    Eyext = Eyext + Ey
!      Ez
                    Ez = -rmultip(lm,ia)  * re(l+1-mm) * zlma( (l+1)*(l+2)+m+1 ) * ra2l3inv
                    Ezext = Ezext + Ez
                enddo
            enddo
            kntlargo = kntlargo + 1
        else		! Short-range
!     locates the fitting interval for ra
            interv = indintrv(int(fct*ra)+1)
            t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
            dost = t + t
            tcheb(0) = uno	! Chebyshev T  polynomials
            tcheb(1) = t
            do j = 2, mxlenpol-1
                tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
            enddo
!     Integrals:
!		rinta = Integrate[ r**(2*la+2) * fradtr[la,ma,r], {r,l_(i-1),r0}]
!		rintb = Integrate[ r * fradtr[la,ma,r], {r,r0,l_i}]
!		computed from fitting polynomials	
            lm = 0
            ra2l3inv = uno / ra
            do l = 0, ltop
                pi4d2l1 = cuatro * pi * dosl1i(l)
                ra2l3inv = ra2l3inv  * ra2inv	! 1.d0 / (ra**(2*l+3))
                lp1lp2 = (l+1)*(l+2)
                lm1l = (l-1)*l
                do m = -l, l
                    lm = lm + 1
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
                    mm = iabs(m)
                    if (m .ge. 0) then
                        isgm = 1
                    else
                        isgm = -1
                    endif
                    isgmm1 = isgm*(mm-1)
                    isgmp1 = isgm*(mm+1)
!      Ex
                    c1 = zlma( lp1lp2+isgmp1+1 )
                    if (l .gt. mm+1) then
                        c2 = zlma(lm1l + isgmp1 + 1)
                    else
                        c2 = cero
                    endif
                    if (m .ne. 0 .and. m .ne. -1) then
                        c1 = c1 - re(l+1-mm) * re(l+2-mm) * zlma( lp1lp2+isgmm1+1 )
                        c2 = c2 - re(l+mm) * re(l+mm-1) * zlma( lm1l+isgmm1+1 )
                    else
                        if (m .eq. 0) then
                                c1 = c1 + c1
                                c2 = c2 + c2
                        endif
                    endif
                    Exext = Exext -umed * pi4d2l1 * (c1 * rinta * ra2l3inv + c2 * rintb)

!      Ey
                    c1 = zlma( lp1lp2-isgmp1+1 )
                    if (l .gt. mm+1) then
                        c2 = zlma(lm1l-isgmp1+1 )
                    else
                        c2 = cero
                    endif
                    if (m .ne. 0 .and. m .ne. 1) then
                        c1 = c1 + re(l+1-mm) * re(l+2-mm) * zlma( lp1lp2-isgmm1+1 )
                        c2 = c2 + re(l+mm) * re(l+mm-1) * zlma( lm1l-isgmm1+1 )
                    else
                        if (m .eq. 0) then
                            c1 = c1 + c1
                            c2 = c2 + c2
                        endif
                    endif
                    Eyext = Eyext - isgm * umed * pi4d2l1 * (c1 * rinta * ra2l3inv + c2 * rintb)
!      Ez	
                    c1 = re(l+1-mm) * zlma( lp1lp2+m+1 )
                    if (l .gt. mm) then
                        c2 = re(l+mm) * zlma( lm1l+m+1 )
                    else
                        c2 = cero
                    endif
                    Ezext = Ezext - pi4d2l1 * (c1 * rinta * ra2l3inv - c2 * rintb)
                enddo     ! End of Do over m
            enddo     ! End of Do over l
            kntcorto = kntcorto + 1
        endif
    enddo     ! Fin del Do en ia
    ex = Exext + Exnucl
    ey = Eyext + Eynucl
    ez = Ezext + Eznucl
    return
    end

!   ***************************************************************

  subroutine campolargo(x, y, z, ex, ey, ez, lnucleo)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    USE DAMFIELD320_D
    implicit none
    logical lnucleo
    integer(KINT) :: i, ia, ierr, interv, isgm, l, lm, ltop, m, mm
    real(KREAL) :: c1, c2, Ex, Exext, Exnucl, Ey, Eyext, Eynucl, Ez, Ezext, Eznucl
    real(KREAL) :: ra, ra2, ra2inv, ra3inv, ra2l3inv, x, xa, y, ya, z, za
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
!    Long-range Electric field in point (x,y,z)
    Exnucl = cero
    Eynucl = cero
    Eznucl = cero
    Exext = cero
    Eyext = cero
    Ezext = cero
    lnucleo = .false.
    do ia = 1, ncen
        xa = x - rcen(1,ia)
        ya = y - rcen(2,ia)
        za = z - rcen(3,ia)
        ra2 = xa*xa + ya*ya + za*za
        ra = sqrt(ra2)
        if (ra .lt. 1.d-1) then
            write(6,"('The point ',3(1x,F10.5),' is too close to center ',i4, ' for long-range computation')") x, y, z, ia
            lnucleo = .true.
            return
        endif
        ra3inv = uno / (ra * ra2)
        Exnucl = Exnucl + zn(ia) * xa * ra3inv
        Eynucl = Eynucl + zn(ia) * ya * ra3inv
        Eznucl = Eznucl + zn(ia) * za * ra3inv
        if (ngini(ia) .le. 0) cycle
        ra2inv = uno / ra2
        if (ra .lt. rlargo(ia)) then
            interv = indintrv(int(fct*ra)+1)
            ltop = lcorto(interv,ia)
        else
             ltop = llargo(min(int(ra),mxlargo),ia)
        endif
        zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
        zlma(2) = ya
        zlma(3) = za
        zlma(4) = xa
        do l = 1, ltop
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))		! zlma(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))		! zlma(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)				! zlma(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)					! zlma(l+1,-l,ia)
            do m = 0, l-1
                    zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlma(l+1,m,ia)
                    zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlma(l+1,-m,ia)
            enddo
        enddo
        lm = 0
        ra2l3inv = uno / ra
        do l = 0, ltop
            ra2l3inv = ra2l3inv  * ra2inv	! 1.d0 / (ra**(2*l+3))
            do m = -l, l
                lm = lm + 1
                mm = iabs(m)
                if (m .ge. 0) then
                    isgm = 1
                else
                    isgm = -1
                endif
!      Ex
                c1 = zlma( (l+1)*(l+2)+isgm*(mm+1)+1 )
                if (m .ne. 0 .and. m .ne. -1) then
                    c1 = c1 - re(l+1-mm) * re(l+2-mm) * zlma( (l+1)*(l+2) + isgm*(mm-1) + 1 )
                else
                    if (m .eq. 0) then
                          c1 = c1 + c1
                    endif
                endif
                c1 = c1 * ra2l3inv
                Ex = -umed * rmultip(lm,ia) * c1
                Exext = Exext + Ex
                c1 = zlma( (l+1)*(l+2)-isgm*(mm+1)+1 )
                if (m .ne. 0 .and. m .ne. 1) then
                    c1 = c1 + re(l+1-mm) * re(l+2-mm) * zlma( (l+1)*(l+2)-isgm*(mm-1)+1 )
                else
                    if (m .eq. 0) then
                        c1 = c1 + c1
                    endif
                endif
                c1 = c1 * ra2l3inv
                Ey = -isgm * umed * rmultip(lm,ia) * c1
                Eyext = Eyext + Ey
!      Ez
                c1 = re(l+1-mm) * zlma( (l+1)*(l+2)+m )
                if (l .gt. mm) then
                    c2 = re(l+mm) * zlma( (l-1)*l+m )
                else
                    c2 = cero
                endif
                c1 = c1 * ra2l3inv
                Ez = -rmultip(lm,ia)  * c1
                Ezext = Ezext + Ez
            enddo
        enddo
    enddo     ! End of DO over ia
    ex = Exext + Exnucl
    ey = Eyext + Eynucl
    ez = Ezext + Eznucl
    return
    end
!
!	***************************************************************
!
  subroutine leedamqtfield
    USE MPI
    USE DAM320_D
    USE DAMFIELD320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, ierr2, indnf, indng, interv, j, k, k1, k2, knt, kntlm, l, lenindintrv
    integer(KINT) :: lldummy, lm, m, ncenbas, ncflm, ndummy, nfdummy, nginidummy, ngfindummy, nndummy, nsamples
    real(KREAL) :: aux, bux, dltsample, dost, dummy, fr1, fr2l2, pi4d2l1, r, ra, ral1inv, rainv, ral
    real(KREAL) :: rinta, rintb, rlarex, step, stepmed, suml, suml1, suml2, summ, summ1, summ2, t, xxdummy
    real(KREAL), allocatable :: Qg(:), qp(:)
    real(KREAL) :: tcheb(0:mxlenpol-1)
    nsize = -1
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
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr2)
#elif __INTEL_COMPILER
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr2)
#else
    open (unit=10, file=trim(projectname)//"_2016.damqt", form='unformatted', action = 'read', access='stream', iostat=ierr)
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='unformatted', action = 'read', access='stream', iostat=ierr2)
#endif
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//"_2016.damqt", myrank
        abort = 1
        return
    endif
    if (ierr2 .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//"_2016.dmqtv", myrank
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('Opens files ', a, ' and ', a)") trim(projectname)//"_2016.damqt",&
        &trim(projectname)//"_2016.dmqtv"
    read(10) ncen, nbas, ncaps
    nsize = nsize - sizeof(ncen) - sizeof(nbas) - sizeof(ncaps)
    if (myrank .eq. 0) write(6,"('ncen = ', i4, ' nbas = ', i6, ' nshells = ', i5)") ncen, nbas, ncaps

!	Allocates memory for geometry and basis set

    allocate(atmnam(ncen), ncontr(ncen), nzn(ncen), rcen(3,ncen), zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
            write(6,"('Memory error when allocating atmnam, ncontr, nzn, rcen and zn in processor ',i3)") myrank
            abort = 1
            return
        endif

!	Geometry and nuclear charges

    if (myrank .eq. 0) write(6,"(/24x,'GEOMETRY (BOHR)')")
    if (myrank .eq. 0) write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge')")
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

!	Basis set
    i = 0
    read(10) lsto	! .true. means STO basis, .false. means GTO basis
    if (lsto) then
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating ngini in processor ',i3)") &
            myrank
            abort = 1
            return
        endif
!		Reads basis set data to dummies
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfindummy
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfindummy)
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfindummy
                i = i + 1
                read(10) nfdummy, nndummy, lldummy, xxdummy
                nsize = nsize - sizeof(nfdummy) - sizeof(nndummy) - sizeof(lldummy) - sizeof(xxdummy)
            enddo
        enddo
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating ngini in processor ',i3)") &
            myrank
            abort = 1
            return
        endif
!		Reads basis set data to dummies
        indng = 1
        ncenbas = 0
        do ia = 1, ncen
            read(10) ncontr(ia)
            nsize = nsize - sizeof(ncontr(ia))
            if (ncontr(ia) .le. 0) then
                ngini(ia) = -1
                cycle
            endif
            ncenbas = ncenbas + 1
            ngini(ia) = indng
            indng = indng + ncontr(ia)
            do j = 1, ncontr(ia)
                read(10) nndummy, lldummy
                nsize = nsize - sizeof(nndummy) - sizeof(lldummy)
                read(10) (xxdummy, k = 1, nndummy)
                read(10) (xxdummy, k = 1, nndummy)
                nsize = nsize - 2*nndummy*sizeof(xxdummy)
            enddo
        enddo
    endif

!	Data of density representation
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

    allocate(icfposd(lmtop*nintervaj+1,ncen), xajustd(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating icfpos and xajust in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    nsize = nsize - sizeof(icfposd(:,1)) * ncenbas

    if (myrank .eq. 0 .and. longoutput) write(6,"('Estimated highest size of xajustd   = ', i15, ' bytes')") size(xajustd)
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
        if (.not. lsto .and. ncontr(ia) .le. 0) cycle
        read(10) icfposd(1:lmtop*nintervaj+1,ia)
        nsize = nsize - sizeof(icfposd)
        if (k .gt. 0) then
            icfposd(1:lmtop*nintervaj+1,ia) = icfposd(1:lmtop*nintervaj+1,ia) + icfposd(lmtop*nintervaj+1,k) - 1 
        endif
        k = ia
        read(10) xajustd(1:nintervaj,ia)	! Reads xajust
!     fitting coeficients
        read(10) cfajust(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1) ! Reads cfajust
    enddo

!	Generates an auxiliary index array for determining the interval to which a given r belongs
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

!	Reads auxiliary integrals from file .dmqtv

    allocate(cfrint1(icfposd(lmtop*nintervaj+1,ncen)-1), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfrint1 in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of cfrint1   = ', i15, ' bytes')") size(cfrint1)

    allocate(cfrint2l2(icfposd(lmtop*nintervaj+1,ncen)-1), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfrint2l2 in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of cfrint2l2 = ', i15, ' bytes')") size(cfrint2l2)

    allocate(QGacum(nintervaj*lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating QGacum in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of QGacum    = ', i15, ' bytes')") size(QGacum)

    allocate(Qgpart(nintervaj*lmtop), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating Qgpart in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of Qgpart    = ', i15, ' bytes')") size(Qgpart)

    allocate(qppart(nintervaj*lmtop), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating qppart in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of qppart    = ', i15, ' bytes')") size(qppart)

    allocate(qpacum(nintervaj*lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating qpacum in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of qpacum    = ', i15, ' bytes')") size(qpacum)

    allocate(rlargo(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rlargo in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(rmultip(lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rmultip. Stop')
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rmultip in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of rmultip    = ', i15, ' bytes')") size(rmultip)
    QGacum = cero
    Qgpart = cero
    qppart = cero
    qpacum = cero
    rlargo = cero
    rmultip = cero
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) cycle
!	multipolar moments
            read(11) rmultip(1:lmtop,ia)
!	Reads the integrals:
!		Qg(la,ma,i;ia) = 
!			Integrate[ r**(2*la+2) * fradtr[la,ma,r], {r,l_(i-1),l_i}]
!		qp(la,ma,i;ia) = 
!			Integrate[ r * fradtr[la,ma,r], {r,l_(i-1),l_i}]
!	and computes from them and stores the integrals:
!		QGacum(la,ma,i;ia) = 
!			Integrate[ r**(2*la+2) * fradtr[la,ma,r], {r,0,l_i}]
!		qpacum(la,ma,i;ia) = 
!			Integrate[ r * fradtr[la,ma,r], {r,l_(i-1),Infinity}]
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
        read(11) cfrint1(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)	! Expansion coefficients of auxiliary integrals rint1
        read(11) cfrint2l2(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)	! Expansion coefficients of auxiliary integrals rint2l2
    enddo
    close(10)
    close(11)

!	Determines the long-range radii and the highest l in the expansion for each interval
    allocate(lcorto(nintervaj,ncen), llargo(0:mxlargo,ncen), Qllargo(0:lmaxrep), stat = ierr )
    if (ierr .ne. 0) call error(1,'Memory error when allocating lcorto, llargo and Qllargo. Stop')
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of llargo   = ', i15, ' bytes')") size(llargo)
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of Qllargo   = ', i15, ' bytes')") size(Qllargo)
!	long-range radii
    allocate(umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating umedpow in processor ',i3)") myrank
        abort = 1
        return
    endif
    umedpow(0) = uno							!
    do i = 1, lmaxexp							!
            umedpow(i) = umedpow(i-1) * umed			! 1 / 2^i
    enddo
    if (myrank .eq. 0) write(6,"('Long-range threshold = ', e12.5)") umbrlargo
    nsamples = 4
    do ia = 1, ncen
        rlargo(ia) = cero
        if (ngini(ia) .le. 0) cycle
        kntlm = 0
        do l = 0, lmaxrep
            summ1 = cero
            do m = -l, l
                kntlm = kntlm + 1
                summ1 = summ1 + abs(rmultip(kntlm,ia)) * fact(l+abs(m)) * umedpow(abs(m)) * facti(l-abs(m)) * facti(abs(m))
            enddo
            Qllargo(l) = summ1
        enddo
        rlargo(ia) = rinterv(nintervaj)
        lcorto(1:nintervaj,ia) = 0
        do interv = 1, nintervaj
            dltsample = udec * (rinterv(interv) - rinterv(interv-1))
            do i = 0, nsamples-1	! samples over nsamples points in each interval to decide the highest l
                ra = rinterv(interv-1) + dltsample + (rinterv(interv) - rinterv(interv-1) - dos * dltsample) &
                        * ri(nsamples-1) * re(i)
                t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
                dost = t + t
                tcheb(0) = uno	! Chebyshev T  polynomials
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
    deallocate (Qgpart, qppart, umedpow)
    return
    end
!**********************************************************************
!    subroutine consta
!
!	Computes and stores auxiliary constants
!		re(i) = dfloat(i)
!		ri(i) = 1.d0 / dfloat(i)
!		fact(i) = dfloat(i!)
!		facti(i) = 1.d0 / dfloat(i!)
!		facts(i) = dfloat((i+1/2)!)
!		ind(i) = i*(i+1)/2
!    	ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) 
!			/ (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMFIELD320_D
    implicit none
    integer(KINT) :: i
!	auxiliary parameters and functions
    pi = acos(-uno)
    raizpi = sqrt(pi)
    re(0) = cero
    ri(0) = 1.d300
    dosl1(0) = uno
    dosl1i(0) = uno
    do i = 1, mxreal
        re(i) = re(i-1) + uno        ! dfloat(i)
        re(-i) = -re(i)
        ri(i) = uno / re(i)       	! uno / dfloat(i)
        ri(-i) = -ri(i)
        dosl1(i) = re(i) + re(i) + uno	! dfloat(iì)
        dosl1(-i) = -re(i) - re(i) + uno
        dosl1i(i) = uno / dosl1(i)		! dfloat( 1/(i+i+1) )
        dosl1i(-i) = -dosl1i(i)
    enddo
    fact(0) = uno
    facti(0) = uno
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)   		!  i!
        facti(i) = uno / fact(i)     		!  uno / i!
    enddo
    return
    end
!
!   ***************************************************************
!     Calculates the electric field and its derivatives from the represented density at point (x,y,z)
!
  subroutine fielddrvs(x, y, z, ex, ey, ez, dxx, dxy, dxz, dyy, dyz, dzz, lnucleo)
    USE DAM320_D
    USE DAMFIELD320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    implicit none
    logical :: lnucleo
    integer(KINT) :: i, ia, ierr, interv, icflm, j, jshft, kntlm, l, lm, lmtopot, ltop, m
    real(KREAL) :: aux, bux, dost, ex, ey, ez, drvvlm, drv2vlm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, flm, fux
    real(KREAL) :: pi4exp, ra, ra2, rainv, ra2inv, rinta, rintb, sgn, t, tp, vlm, vqlm
    real(KREAL) :: x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL), parameter :: vtope = 1.d10		! To prevent infinity, if the point coincides with a nucleus, loads the value of parameter vtope
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
    ex = cero
    ey = cero
    ez = cero
    dxx = cero
    dxy = cero
    dxz = cero
    dyy = cero
    dyz = cero
    dzz = cero
    lnucleo = .false.
!    Computes the derivatives of the electrostatic potential
    do ia = 1, ncen
        xa = x - rcen(1,ia)
        ya = y - rcen(2,ia)
        za = z - rcen(3,ia)
        ra2 = xa*xa+ya*ya+za*za
        if (ra2 .lt. 1.d-20) then
            lnucleo = .true.
            return
        endif
        ra = sqrt(ra2)
        rainv = uno / ra
        ra2inv = rainv * rainv
        if (ngini(ia) .le. 0) then
            ex = ex + zn(ia) * rainv * xa * ra2inv
            ey = ey + zn(ia) * rainv * ya * ra2inv
            ez = ez + zn(ia) * rainv * za * ra2inv
            dxx =  dxx  - (re(3) * ex * xa + zn(ia) * rainv) * ra2inv
            dxy =  dxy  - re(3) * ex * ya * ra2inv
            dxz =  dxz  - re(3) * ex * za * ra2inv
            dyy =  dyy  - (re(3) * ey * ya + zn(ia) * rainv) * ra2inv
            dyz =  dyz  - re(3) * ey * za * ra2inv
            dzz =  dzz  - (re(3) * ez * za + zn(ia) * rainv) * ra2inv
            cycle
        endif
        if (ra .lt. rlargo(ia)) then
            interv = indintrv(int(fct*ra)+1)
            ltop = lcorto(interv,ia)
        else
            ltop = llargo(min(int(ra),mxlargo),ia)
        endif
        lmtopot = (ltop+1)*(ltop+1)
        xadivra = xa * rainv
        yadivra = ya * rainv
        zadivra = za * rainv
        zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
        zlma(2) = ya
        zlma(3) = za
        zlma(4) = xa
        ra2l1(1) = ra
        ra2l1inv(1) = rainv
        pi4d2l1(1) = cuatro * pi
        d2l1(1) = uno
        lm = 1
        do l = 1, ltop
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))		! zlm(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))		! zlm(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)				! zlm(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)					! zlm(l+1,-l,ia)
            do m = 0, l-1
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlm(l+1,m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlm(l+1,-m,ia)
            enddo
            aux = ra2l1(lm) * ra2
            bux = ra2l1inv(lm) * ra2inv
            do m = -l, l
                lm = lm + 1
                ra2l1(lm) = aux
                ra2l1inv(lm) = bux
                pi4d2l1(lm) = cuatro * pi * dosl1i(l)	!   4 * pi / (2l+1)
                d2l1(lm) = dosl1(l)
            enddo
        enddo
        call derivzlm(ltop, idimzlm, zlma, zlmadx, zlmady, zlmadz)
        call derivzlm(ltop, idimzlm, zlmadx, zlmadxx, zlmadxy, zlmadxz)
        call dzlm2y(ltop, idimzlm, zlmady, zlmadyy, zlmadyz)
        call dzlm2z(ltop, idimzlm, zlmadz, zlmadzz)
        if (ra .ge. rlargo(ia)) then  ! The point is in the long-range region (ra >= rlargo(ia) )
            aux = zn(ia)	! aux is the nuclear charge for l = 0 and zero otherwise
            do lm = 1, lmtopot
                vlm = (aux - rmultip(lm,ia) ) * ra2l1inv(lm)
                drvvlm = -d2l1(lm) * vlm * rainv
                drv2vlm = -(d2l1(lm)+uno) * drvvlm * rainv
                ex = ex + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
                ey = ey + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
                ez = ez + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
                fux = rainv * rainv * rainv * (ra * drv2vlm - drvvlm)
                dxx = dxx + (xa*xa * fux + drvvlm * rainv) * zlma(lm) &
                        + dos * xadivra * drvvlm * zlmadx(lm) + vlm * zlmadxx(lm)
                dyy = dyy + (ya*ya * fux + drvvlm * rainv) * zlma(lm) &
                        + dos * yadivra * drvvlm * zlmady(lm) + vlm * zlmadyy(lm)
                dzz = dzz + (za*za * fux + drvvlm * rainv) * zlma(lm) &
                        + dos * zadivra * drvvlm * zlmadz(lm) + vlm * zlmadzz(lm)
                dxy = dxy + xa*ya * fux * zlma(lm) + xadivra * drvvlm * zlmady(lm) &
                        + yadivra * drvvlm * zlmadx(lm) + vlm * zlmadxy(lm)
                dxz = dxz + xa*za * fux * zlma(lm) + xadivra * drvvlm * zlmadz(lm) &
                        + zadivra * drvvlm * zlmadx(lm) + vlm * zlmadxz(lm)
                dyz = dyz + ya*za * fux * zlma(lm) + yadivra * drvvlm * zlmadz(lm) &
                                + zadivra * drvvlm * zlmady(lm) + vlm * zlmadyz(lm)
                aux = cero
            enddo
        else		!     The point is in the short-range region (ra < rlargo(ia) )
            t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
            dost = t + t
            tcheb(0) = uno	! Chebyshev T  polynomials
            tcheb(1) = t
            do j = 2, mxlenpol-1
                tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
            enddo
            aux = zn(ia)	! aux is the nuclear charge for l = 0 and zero otherwise
            pi4exp = cuatro * pi * exp(-xajustd(interv,ia)*ra)
            do lm = 1, lmtopot
                if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo) cycle
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
                vlm = (aux - pi4d2l1(lm) * ( rinta + ra2l1(lm) * rintb)) * ra2l1inv(lm)
                vqlm = (aux - pi4d2l1(lm) * rinta) * ra2l1inv(lm)
                drvvlm = -d2l1(lm) * vqlm * rainv 	! d2l1(lm) = (l+l+1)
                drv2vlm = -(d2l1(lm)+uno) * drvvlm * rainv + pi4exp * flm
                ex = ex + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
                ey = ey + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
                ez = ez + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
                fux = rainv * rainv * rainv * (ra * drv2vlm - drvvlm)
                dxx = dxx + (xa*xa * fux + drvvlm * rainv) * zlma(lm) &
                        + dos * xadivra * drvvlm * zlmadx(lm) + vlm * zlmadxx(lm)
                dyy = dyy + (ya*ya * fux + drvvlm * rainv) * zlma(lm) &
                        + dos * yadivra * drvvlm * zlmady(lm) + vlm * zlmadyy(lm)
                dzz = dzz + (za*za * fux + drvvlm * rainv) * zlma(lm) &
                        + dos * zadivra * drvvlm * zlmadz(lm) + vlm * zlmadzz(lm)
                dxy = dxy + xa*ya * fux * zlma(lm) + xadivra * drvvlm * zlmady(lm) &
                        + yadivra * drvvlm * zlmadx(lm) + vlm * zlmadxy(lm)
                dxz = dxz + xa*za * fux * zlma(lm) + xadivra * drvvlm * zlmadz(lm) &
                        + zadivra * drvvlm * zlmadx(lm) + vlm * zlmadxz(lm)
                dyz = dyz + ya*za * fux * zlma(lm) + yadivra * drvvlm * zlmadz(lm) &
                        + zadivra * drvvlm * zlmady(lm) + vlm * zlmadyz(lm)
                aux = cero
            enddo
        endif  ! End of test over long/short-range
    enddo
!    Changes the sign of the derivatives of the potential to get the electric field and its derivatives
    ex = -ex
    ey = -ey
    ez = -ez
    dxx = -dxx
    dxy = -dxy
    dxz = -dxz
    dyy = -dyy
    dyz = -dyz
    dzz = -dzz
    return
    end

!   ***************************************************************

  subroutine derivzlm(lmax, idimzlm, zlma, zlmadx, zlmady, zlmadz)
    USE DAM320_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm)
!	Derivatives of the regular harmonics with respecto to the Cartesian coordinates
    zlmadx(1) = cero	! Derivatives of the S spherical harmonics
    zlmady(1) = cero
    zlmadz(1) = cero
    zlmadx(2) = cero	! Derivatives of the P spherical harmonics
    zlmadx(3) = cero
    zlmadx(4) = zlma(1)
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    zlmadx(5) = re(6) * zlma(2)	! Derivatives of the D spherical harmonics
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
    zlmadx(10) = re(15) * zlma(5)		! Derivatives of the F spherical harmonics
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
    do l = 4, lmax		! Derivatives of the remaining spherical harmonics
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
    do l = 4, lmax		! Derivatives of the remaining spherical harmonics
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
    do l = 1, lmax		! Derivatives of the remaining spherical harmonics
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
!	-------------------------------------------------------------------------------------------------------
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
