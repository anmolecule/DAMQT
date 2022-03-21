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
! Program for tabulation of the electronic density from the representation of the molecular density performed with
! DAM320
!
! Parallel version with MPI
!
! Version of September 2018
!

! #define DBLPRCGRID	! Uncomment this line  if double precision grid is wanted
  program DAMDEN320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDEN320_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, ierr, ios, ipoint, j, lmaxctot, nbasis
    real(KREAL) :: aux, den, denlp, denrep, dendrvx, dendrvy, dendrvz, denlplc, dxden, dxx, dxxden
    real(KREAL) :: dxy, dxyden, dxz, dxzden, dyden, dyy, dyyden, dyz, dyzden, dzden, dzz, dzzden
    real(KREAL) :: x, xmax, xmin, xyzmax, xyzmin, y, ymax, ymin, z, zmax, zmin
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    logical :: lnamelist(12), ltimeprocs
    integer(KINT) :: inamelist(2)
    real(KREAL) :: rnamelist(10)
    namelist / options / dltu, dltv, dltx, dlty, dltz, filename, iatomsel, iswindows, langstrom, laplacian, latomics &
            , latomsel, lboundsx, lboundsy, lboundsz, ldeform, ldensacc, lderiv2, lexact, lgradient, lgrid, lgrid2d, lmaxrep &
            , lminrep, lmolec, longoutput, lpoints, nsel, numrtab, planeA, planeB, planeC, planecase &
            , rtab, uinf, umbrlargo, usup, vinf, vsup, xboundinf, xboundsup, x_func_uv &
            , xinf, xsup, yboundinf, yboundsup, y_func_uv, yinf, ysup, zboundinf, zboundsup, z_func_uv, zinf, zsup
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    abort = 0
    abortroot = 0
    if (myrank .eq. 0) write(6,"('number of processors = ', i3)") nprocs
    tiempo = dtime(tarray)
!	Defaults for the NAMELIST OPTIONS
!		Some variables included in the namelist are kept for compatibility with the input of the scalar version of the program
!		but they have no effect in this parallel code. Therefore, they are neither assigned nor used.
    langstrom = .true.		! If .false. distances of grids (.plt files) in bohr
    longoutput = .false.	! If .true. a more detailed output is given
    ldeform = .false.           ! If .true. density deformations are computed (lminrep = 1)
    lexact  = .false.		! if .true. "exact" density is tabulated
    lminrep = 0			! lowest "l" in the expansion of the density. It does not hold if "lexact .eq. true"
    lmaxrep = 10			! highest "l" in the expansion of the density. It does not hold if "lexact .eq. true"
                                            ! the expansion of the density takes lminrep <= l <= lmaxrep
    lmolec = .true.		! If .true. tabulation of the whole molecular density
    latomics  = .false.		! If .true. tabulation of atomic densities stored in files named:
                                            !    projectname-cxx-d.plt (gOpenMol)
                                            ! where xx refers to the center corresponding to the tabulation
    latomsel = .false.		! If .true. indices of the centers for atomic tabulations will be user-supplied
                                            ! The indices of the selected atoms must be supplied in vector "iatomsel".
                                            ! Maximum number of centers that can be selected "mxsel" (parameter).
    lgradient = .false.		! If .true. gradient components of the density tabulated in files
                                            !    projectname-d-dx.pltd, projectname-d-dy.pltd, projectname-d-dz.pltd, for molecular density
    laplacian = .false.		! If .true. laplacian tabulated in file projectname-dlplc.plt
    lderiv2 = .false.		! If .true. second derivatives of the density computed and, if lgrid = .true., tabulated in files
                                            !	projectname-d-dxx.pltd, projectname-d-dxy.pltd, projectname-d-dxz.pltd, projectname-d-dyz.pltd, etc
    iatomsel = 0			! To cover the possibility of an input file with "latomsel = .true."
    iatomsel(1) = 1		! but without assignment of "iatomsel".
    ldensacc  = .false.		! If .true. computes and tabulates the accumulated density of the selected atoms in file projectname-frg-d.plt
    lboundsx = .false.		! If .true. the evaluation of deformation charges is constrained to a range in x
    xboundinf = cero		! Lower limit of that range
    xboundsup = cero		! Upper limit of that range
    lboundsy = .false.		! If .true. the evaluation of deformation charges is constrained to a range in y
    yboundinf = cero		! Lower limit of that range
    yboundsup = cero		! Upper limit of that range
    lboundsz = .false.		! If .true. the evaluation of deformation charges is constrained to a range in z
    zboundinf = cero		! Lower limit of that range
    zboundsup = cero		! Upper limit of that range
    lgrid = .true.			! If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    lpoints = .false.		! If .true. computes in selected points and prints out the results in the standard output.
                                            ! Points must be given in cartesian coordinates.
                                            ! If lgrid .eq. .true., these coordinates must be placed after the grid data
    filename = ""			! root file name for .plt nad .pltd files
    umbrlargo = 1.d-8		! Long-range threshold
    iswindows = .false.		! .true. if running on a MS-windows system
    numrtab = 0			! Number of tabulation points supplied in namelist
    rtab = cero			! Tabulation points supplied in namelist
    xinf = cero
    xsup = cero
    dltx = uno
    yinf = cero
    ysup = cero
    dlty = uno
    zinf = cero
    zsup = cero
    dltz = uno
    uinf = cero
    usup = uno
    dltu = uno
    vinf = cero
    vsup = uno
    dltv = uno

!       The following namelist variables are not used in this program. They are included for DAMQT compatibility.
    planeA = cero        ! Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    planeB = cero
    planeC = uno
    planecase = 1
!	End of namelist defaults

    ltimeprocs = .false.
    if (myrank .eq. 0) then
        read(5,OPTIONS)
        read(5,*) projectname
        if (lderiv2) lgradient = .true.
        inamelist = (/ lmaxrep, lminrep /)
        rnamelist = (/ xinf, xsup, dltx, yinf, ysup, dlty, zinf, zsup, dltz, umbrlargo /)
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
        if (len_trim(filename).eq.0) then
            filename = projectname
        else
            filename = projectname(1:i)//trim(filename)
        endif
        lgbsgz = .false.	! Checks whether the .ggbs or .sgbs file is gzipped or not
        inquire(file=trim(projectname)//".ggbs.gz", exist=lgbsgz, iostat=ierr)
        if (ierr .eq. 0 .and. lgbsgz) then
            call system ("gunzip "//trim(projectname)//".ggbs.gz")
        else
            inquire(file=trim(projectname)//".sgbs.gz", exist=lgbsgz, iostat=ierr)
            if (ierr .eq. 0 .and. lgbsgz) then
                    call system ("gunzip "//trim(projectname)//".sgbs.gz")
            endif
        endif
        ldengz = .false.	! Checks whether the .den file is gzipped or not
        if (lexact) then
            inquire(file=trim(projectname)//".den.gz", exist=ldengz, iostat=ierr)
            if (ierr .eq. 0 .and. ldengz) then
                    call system ("gunzip "//trim(projectname)//".den.gz")
            endif
        endif
        allocate(timeprocs(2*nprocs), stat = ierr)
        if (ierr .eq. 0) then
            ltimeprocs = .true.
            timeprocs = 0.
        else
            write(6,"('WARNING: Memory error when allocating timeprocs, ierr =  ',i5)") ierr
            ltimeprocs = .false.
        endif
        lnamelist = (/ langstrom, lgrid, lgradient, laplacian, lderiv2, lexact, ldensacc, latomics, latomsel, &
                    lmolec, ldeform, ltimeprocs /)
    endif
    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lnamelist,12,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,10,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    if (myrank .ne. 0) then
        lmaxrep = inamelist(1); lminrep = inamelist(2)
        langstrom = lnamelist(1); lgrid = lnamelist(2);  lgradient = lnamelist(3);  laplacian = lnamelist(4);
        lderiv2 = lnamelist(5); lexact = lnamelist(6); ldensacc = lnamelist(7); latomics = lnamelist(8);
        latomsel = lnamelist(9); lmolec = lnamelist(10); ldeform = lnamelist(11); ltimeprocs = lnamelist(12);
        xinf = rnamelist(1); xsup = rnamelist(2); dltx = rnamelist(3)
        yinf = rnamelist(4); ysup = rnamelist(5); dlty = rnamelist(6)
        zinf = rnamelist(7); zsup = rnamelist(8); dltz = rnamelist(9)
        umbrlargo = rnamelist(10)
    endif

    if (lexact) ldeform = .false. ! Deformations computation is allowed only for represented density

    if (ldeform) lminrep = 1

    call consta		! Computes auxiliary constants

    call leedamqtden		! Reads file .damqt  (generated by DAM2016)
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
            call error(1,'Stop')
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (myrank .eq. 0) then
        if (lgbsgz) then
            inquire(file=trim(projectname)//".ggbs", exist=lgbsgz, iostat=ierr)
            if (ierr .eq. 0 .and. lgbsgz) then
                call system ("gzip "//trim(projectname)//".ggbs")
            else
                inquire(file=trim(projectname)//".sgbs", exist=lgbsgz, iostat=ierr)
                if (ierr .eq. 0 .and. lgbsgz) then
                    call system ("gzip "//trim(projectname)//".sgbs")
                endif
            endif
        endif
    endif

#ifdef DBLPRCGRID
    if (myrank .eq. 0) then
        write(6,"('Grid generated in double precision')")
        write(6,"('WARNING! This grid will not be compatible with gOpenMol')")
    endif
#endif

    allocate(lselat(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lselat in processor ',i3)") myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    do ia = 1, ncen
        lselat(ia) = .false.
    enddo

    if (myrank .eq. 0) then
!	selected atomic fragments for tabulation	
        if (latomsel) then
            do i = 1, mxsel
                if (iatomsel(i) .le. ncen .and. iatomsel(i) .gt. 0) lselat(iatomsel(i)) = .true.
            enddo
            nsel = 0
            do ia = 1, ncen
                if (lselat(ia)) then
                    nsel = nsel+1
                    iatomsel(nsel) = ia
                endif
            enddo
        elseif (latomics ) then
            nsel = ncen
            latomsel = .true.
            if (ncen .gt. mxsel) then
                if (myrank .eq. 0) then
                    write(6,"('WARNING: Number of selected atoms larger than allowed for tabulation of the atomic densities')")
                    write(6,"('Only the first ',i3, ' atomic densities will be tabulated')") mxsel
                endif
                nsel = mxsel
            endif
            do i = 1, nsel
                lselat(i) = .true.
                iatomsel(i) = i
            enddo
        endif

        if (latomics) write(6,"('Indices of atoms selected for tabulation = ', 20(i6))") (iatomsel(i), i = 1, nsel)
        if (ldensacc) write(6,"('Indices of atoms selected for density accumulation = ', 20(i6))") (iatomsel(i), i = 1, nsel)
    endif
    CALL MPI_BCAST(nsel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lselat,ncen,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(iatomsel,nsel,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (lminrep .gt. lmaxexp .and. .not. lexact) then
        if (myrank .eq. 0) then
            write(6,"('lminrep = ', i3, ' greater than lmaxexp ', i3)") lminrep, lmaxexp
            write(6,"('takes lminrep = 0')")
        endif
        lminrep = 0
    endif
    if (myrank .eq. 0 .and. .not. lexact) write(6,"('Expansion of the density: lminrep = ', i3, ' : lmaxrep = ', i3)") &
        lminrep, lmaxrep
    if (lexact) then

        lgradient = .false.	! The gradient is not computed for the exact density
        lmolec = .true.
!         Seeks for density file
        ldengz = .false.	! Checks whether the .den file is gzipped or not
        inquire(file=trim(projectname)//".den.gz", exist=ldengz, iostat=ierr)
        if (ierr .eq. 0 .and. ldengz) then
            call system ("gunzip "//trim(projectname)//".den.gz")
        endif
        lden = .false.
        ldensprsbin = .false.
        inquire(file=trim(projectname)//".den", exist=lden, iostat=ierr)
        if (.not.(ierr .eq. 0 .and. lden)) then
            inquire(file=trim(projectname)//".densprsbin", exist=ldensprsbin, iostat=ierr)
        endif
        if (.not. (lden .or. ldensprsbin)) then
            call error(1,'No density matrix available. Stop')
        endif
        if (lden) then
            call readden
            if (ldengz) then	! restores files back to their original gzipped status
                call system ("gzip "//trim(projectname)//".den")
            endif
            if (lgbsgz) then
                call system ("gzip "//trim(projectname)//".ggbs")
            endif
        else if(ldensprsbin) then
            call readdensprsbin
        else
            write(6,"('No density matrix available ')")
            abort = 1
        endif
        CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if (abortroot .gt. 0) then
            call error(1,'Stop')
        endif
        close(16)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (myrank .eq. 0) then
            if (ldengz) then	! restores files back to their original gzipped status
                call system ("gzip "//trim(projectname)//".den")
            endif
        endif

! 		Allocates auxiliary arrays for exact density
        allocate(f(nbas), faux(nbas), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating f and faux in processor ',i3)") myrank
            abort = 1
        endif
    endif


    if (lgrid) then
        if ((xsup-xinf) .eq. cero .or. (ysup-yinf) .eq. cero .or. (zsup-zinf) .eq. cero) then	!	Default grid
            xmin = minval(rcen(1,1:ncen))
            xmax = maxval(rcen(1,1:ncen))
            ymin = minval(rcen(2,1:ncen))
            ymax = maxval(rcen(2,1:ncen))
            zmin = minval(rcen(3,1:ncen))
            zmax = maxval(rcen(3,1:ncen))
            xyzmin = min(xmin, ymin, zmin)
            xyzmax = max(xmax, ymax, zmax)
            if ((xmin-xmax) .eq. cero) then
                xmin = min(xyzmin,-dos)
                xmax = max(xyzmax,dos)
            endif
            if ((ymin-ymax) .eq. cero) then
                ymin = min(xyzmin,-dos)
                ymax = max(xyzmax,dos)
            endif
            if ((zmin-zmax) .eq. cero) then
                zmin = min(xyzmin,-dos)
                zmax = max(xyzmax,dos)
            endif
            if ((xsup-xinf) .eq. cero) then
                xinf = umed * (re(5) * xmin - re(4) * xmax)
                xsup = umed * (re(5) * xmax - re(4) * xmin)
                dltx = (xsup-xinf)/49
            endif
            if ((ysup-yinf) .eq. cero) then
                yinf = umed * (re(5) * ymin - re(4) * ymax)
                ysup = umed * (re(5) * ymax - re(4) * ymin)
                dlty = (ysup-yinf)/49
            endif
            if ((zsup-zinf) .eq. cero) then
                zinf = umed * (re(5) * zmin - re(4) * zmax)
                zsup = umed * (re(5) * zmax - re(4) * zmin)
                dltz = (zsup-zinf)/49
            endif
        endif
        if (xinf .gt. xsup) then
            aux = xsup
            xsup = xinf
            xinf = aux
            dltx = abs(dltx)
        endif
        if (yinf .gt. ysup) then
            aux = ysup
            ysup = yinf
            yinf = aux
            dlty = abs(dlty)
        endif
        if (zinf .gt. zsup) then
            aux = zsup
            zsup = zinf
            zinf = aux
            dltz = abs(dltz)
        endif
!		Computes the density on the grid points
        if (lexact) then
            call gridexac
        else
            call gridrep

        endif
        tiempo = dtime(tarray)
        CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!        write(6,"(1x,'Timing in seconds of processor ', i2, ' (user, system, total):',5x,'(',e12.5,',',e12.5,',',e12.5')')") &
!                myrank, tarray(1), tarray(2), tarray(1)+tarray(2)
        if (abortroot .gt. 0) then
            call error(1,'Stop')
        endif
    endif
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
!	Tabulates specific points if required
    if (lpoints .and. myrank .eq. 0) then
        allocate(lnegia(ncen), stat = ierr)
        if (ierr .eq. 0) then
            if (lexact) then
                write(6,"(//12x,'Density from the basis set and density matrix' &
                        ,/12x, 45('='), //9x,'X',17x,'Y',17x,'Z',20x,'den')")
            else
                if (lgradient) then
                    if (laplacian) then
                        write(6,"(//12x,'Density expansion from l = ',i2, ' to  l = ',i2,/10x, 45('='),  &
                                //9x,'X',17x,'Y',17x,'Z',20x,'den',t89,'der x', t112, 'der y', t135, 'der z',  &
                                t156,'laplacian')") lminrep, lmaxrep
                    else
                        write(6,"(//12x,'Density expansion from l = ',i2, ' to  l = ',i2,/10x, 45('='), &
                                //9x,'X',17x,'Y',17x,'Z',20x,'den',t89,'der x', t112, 'der y', t135, 'der z')") &
                                lminrep, lmaxrep
                    endif
                else
                    if (laplacian) then
                        write(6,"(//12x,'Density expansion from l = ',i2, ' to  l = ',i2,/10x, 45('='),  &
                                //9x,'X',17x,'Y',17x,'Z',t67,'den', t89,'laplacian')") lminrep, lmaxrep
                    else
                        write(6,"(//12x,'Density expansion from l = ',i2, ' to  l = ',i2,/10x, 45('='),  &
                                //9x,'X',17x,'Y',17x,'Z',t67,'den')") lminrep, lmaxrep
                    endif
                endif
            endif
            latomics = .false.
            if (lderiv2) then
                if (ldeform) then
                    open(9,file=trim(projectname)//"-deform-d.der2",form='formatted', iostat=ierr)
                else
                    open(9,file=trim(projectname)//"-d.der2",form='formatted', iostat=ierr)
                endif
                if (ierr .ne. 0) then
                    write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//"-d.der2", myrank
                    lderiv2 = .false.
                endif
            endif
            if (numrtab .gt. mxrtab) then
                write(6,"('Number of tabulation points numrtab = ', i4, ' exceeds the highest allowable = ', i4 &
                        / ' sets numrtab = ', i4)") numrtab, mxrtab, mxrtab
                numrtab = mxrtab
            endif
            ipoint = 1
            ierr = 0
            if (ipoint .le. numrtab) then
                x = rtab(1,ipoint)
                y = rtab(2,ipoint)
                z = rtab(3,ipoint)
                ipoint = ipoint + 1
            else
                read(5,*,iostat=ierr) x, y, z
            endif
            if (.not. lexact) then
                idimzlm = (lmaxexp+2)**2
                allocate(zlma(idimzlm), stat = ierr)
                if (ierr .ne. 0) then
                    write(6,"('Memory error when allocating zlma in processor ',i3)") myrank
                    abortroot = 1
                endif
                if (lgradient) then
                    allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
                    if (ierr .ne. 0) then
                            write(6,"('Memory error when allocating zlmadx, zlmady, zlmadz in processor ',i3)") myrank
                            lgradient = .false.
                    endif
                endif
                if (lderiv2) then
                    allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                            zlmadzz(idimzlm), stat = ierr)
                    if (ierr .ne. 0) then
                            write(6,"('Memory error when allocating zlmadxx, zlmadxy,... zlmadzz in processor ',i3)") myrank
                            lderiv2 = .false.
                    endif
                endif
            else
                if (.not. allocated(zlma)) then
                    lmaxctot = 0
                    do ia = 1, ncen
                        lmaxctot = max(lmaxc(ia), lmaxctot)
                    enddo
                    allocate(zlma((lmaxctot+2)**2), stat = ierr)
                    if (ierr .ne. 0) call error(1,'Memory error when allocating zlma. Stop')
                endif
            endif
            if (abortroot .gt. 0) then
                call error(1,'Stop')
            endif
            do while (ierr .eq. 0)
                if (lexact ) then       ! Tabula la densidad exacta
                    if (lsto) then
                        call densexactaSTO(x, y, z, den)
                    else
                        call densexactaGTO(x, y, z, den)
                    endif
                else
                    den = cero
                    dxden = cero
                    dyden = cero
                    dzden = cero
                    dxxden = cero
                    dxyden = cero
                    dxzden = cero
                    dyyden = cero
                    dyzden = cero
                    dzzden = cero
                    denlp = cero
                    do ia = 1, ncen
                        call densrepr(ia, x, y, z, denrep, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz)
                        den = den + denrep
                        if (lgradient) then
                            dxden = dxden + dendrvx
                            dyden = dyden + dendrvy
                            dzden = dzden + dendrvz
                        endif
                        if (lderiv2) then
                            dxxden = dxxden + dxx
                            dxyden = dxyden + dxy
                            dxzden = dxzden + dxz
                            dyyden = dyyden + dyy
                            dyzden = dyzden + dyz
                            dzzden = dzzden + dzz
                        endif
                        if (laplacian) denlp = denlp + denlplc
                    enddo
                endif
                if (lgradient) then
                    if (lderiv2) then
                        write(9,"(/'Second derivatives of density at point ',3(1x,e17.10))") x, y, z
                        write(9,"(6(1x,e22.15))") dxxden, dxyden, dxzden, dyyden, dyzden, dzzden
                        if (laplacian) then
                            write(9,"('laplacian comprobation ',2(1x,e22.15))") denlp, dxxden + dyyden + dzzden
                        endif
                    endif
                    if (laplacian) then
                        write(6,"(3(1x,e17.10), 3x, d22.15, 4(1x,e22.15),/)") x, y, z, den, dxden, dyden, dzden, denlp
                    else
                        write(6,"(3(1x,e17.10), 3x, d22.15, 3(1x,e22.15),/)") x, y, z, den, dxden, dyden, dzden
                    endif
                else
                    if (laplacian) then
                        write(6,"(3(1x,e17.10), 3x, d22.15,1x, e22.15/)") x, y, z, den, denlp
                    else
                        write(6,"(3(1x,e17.10), 3x, d22.15,/)") x, y, z, den
                    endif
                endif
                if (ipoint .le. numrtab) then
                    x = rtab(1,ipoint)
                    y = rtab(2,ipoint)
                    z = rtab(3,ipoint)
                    ipoint = ipoint + 1
                else
                    read(5,*,iostat=ierr) x, y, z
                endif
            enddo
            deallocate(lnegia)
            if (allocated(zlma)) deallocate(zlma)
            if (allocated(zlmadx)) deallocate(zlmadx, zlmady, zlmadz)
            if (allocated(zlmadxx)) deallocate(zlmadxx, zlmadxy, zlmadxz, zlmadyy, zlmadyz, zlmadzz)
        endif
        tiempo = dtime(tarray)
        write(6,"(1x,'Timing in seconds of individual points tabulation in proc 0 (user, system, total):', &
                5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tarray(1)+tarray(2)
    endif
    call MPI_FINALIZE(ierr)
    stop
    end
!
!	***************************************************************
!
  subroutine readden
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, icarga, ierr, ios, j, k, k1, k2, knt, nbasis
    real(KREAL) :: aux, bux

!	Allocates the array containing the density matrix
    
    allocate(dmat(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating dmat.')")
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Estimated highest size of dmat   = ', i15, ' bytes')") size(dmat)
    
    if (lsto) then
        open(16,file=trim(projectname)//".den",form='unformatted', iostat=ierr)
        if (ierr .ne. 0) call error(1,'Cannot open file '//trim(projectname)//'.den. Stop')
        read(16,iostat = ios) nbasis, ((dmat(i,j), i=1,nbasis), j=1,nbasis)
    else
        open(16,file=trim(projectname)//".den",form='formatted', iostat=ierr)
        if (ierr .ne. 0) call error(1,'Cannot open file '//trim(projectname)//'.den. Stop')
        read(16,*, iostat = ios) nbasis, ((dmat(i,j),j=1,i),i=1,nbasis)
    endif
    if ( ios .ne. 0 .or. nbas .ne. nbasis  .and. myrank .eq. 0) then
        write(6,"('nbas = ', i5,' nbasis = ', i5)") nbas, nbasis
        write(6,"('ERROR reading density matrix. Check whether the density matrix correspond to this basis set.')")
        abort = 1
        return
    endif
    close(16)
    if (.not. lsto) then
        do i = 2, nbasis
            do j = 1, i-1
                dmat(j,i) = dmat(i,j)
            enddo
        enddo
    endif

    return
    end
!
!	***************************************************************
!
  subroutine readdensprsbin
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ierr, ios

!	Reads the density matrix in lower triangle form:  read ((dmat(i,j), j = 1, i), i = 1, nbas)

#if _WIN32
    open (unit=16, file=trim(projectname)//".densprsbin", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//".densprsbin")
#elif __INTEL_COMPILER
    open (unit=16, file=trim(projectname)//".densprsbin", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//".densprsbin")
#else
    open (unit=16, file=trim(projectname)//".densprsbin", form='unformatted', action = 'read', access='stream', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Memory error when opening file '//trim(projectname)//".densprsbin")
#endif

!	Allocates the array containing the density matrix
    allocate(dvec(nbas*(nbas+1)/2), ivec(nbas*(nbas+1)/2), jvec(nbas*(nbas+1)/2), dmat(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error ', i6,' when allocating dvec, ivec, jvec and dmat. Stop')") ierr
        abort = 1
        return
    endif

    do i=1,nbas*(nbas+1)/2
        read(16,end = 100, err = 100, iostat = ios) ivec(i), jvec(i), dvec(i)
    enddo

100 continue
    if (ios .eq. -1) then
        numdvec = i-1
        if (longoutput) write(6,"('Number of elements of density matrix read = ', i16)") numdvec
    else if ( ios .ne. 0 ) then
        write(6,"('ERROR ', i6, ' reading density matrix. Stop')") ios
        abort = 1
        return
    endif
    close(16)
    
    if ( maxval(ivec) .gt. nbas .or. maxval(jvec) .gt. nbas ) then
        write(6,"('Index on density matrix higher than number of basis functions.',& 
            &/'Check whether the density matrix correspond to this basis set.')")
        write(6,*) 'maxval(ivec) = ', maxval(ivec)
        write(6,*) 'maxval(jvec) = ', maxval(jvec)
        write(6,*) 'nbas = ', nbas
        abort = 1
        return
    endif
    
    if (myrank .eq. 0) then
        write(6,"(/'Sparse density matrix read from binary file', a,/)") trim(projectname)//".densprsbin"
    endif
        
    dmat = cero
    do i = 1, numdvec
        dmat(ivec(i), jvec(i)) = dvec(i)
        dmat(jvec(i), ivec(i)) = dvec(i)
    enddo
    deallocate(dvec)
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
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, l, lm, m
    real(KREAL) :: aux
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
        dosl1i(-i) = uno / dosl1(-i)
    enddo
    fact(0) = uno
    facts(-1) = raizpi
    facts(0) = facts(-1) * umed
    facti(0) = uno
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)   		!  i!
        facts(i) = facts(i-1) * re(i+i+1) * umed	! (i+1/2)!
        facti(i) = uno / fact(i)     		!  uno / i!
    enddo									!
    return
    end
!
!	***************************************************************
!
  subroutine leedamqtden
    USE MPI
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, ncenbas, ncfaj, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, r, ra, ral, rlarex, step, suml, summ, t
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
    if (myrank .eq. 0) write(6,"('ncen = ', i8, ' nbas = ', i8, ' ncaps = ', i8)") ncen, nbas, ncaps

!	Allocates memory for geometry

    allocate(atmnam(ncen), nzn(ncen), rcen(3,ncen), zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam, nzn, rcen and zn in processor ',i3)") myrank
        abort = 1
        return
    endif

    if (myrank .eq. 0) write(6,"(/24x,'GEOMETRY (BOHR)')")
    if (myrank .eq. 0) write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge')")
!	Geometry and nuclear charges
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
    read(10) lsto	! .true. means STO basis, .false. means GTO basis
    nsize = nsize - sizeof(lsto)
    if (lsto) then

!		Allocates memory for the basis set

        allocate(ll(ncaps), lmaxc(ncen), nf(ncaps), ngini(ncen), ngfin(ncen), nn(ncaps), rlargo(ncen), rnor(ncaps), &
                xx(ncaps), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating ll, lmaxc, nf, ngini, ngfin, nn, rlargo, rnor and xx in processor ',i3)") &
            myrank
            abort = 1
            return
        endif

        if (myrank .eq. 0) write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfin(ia)
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfin(ia))
            if (myrank .eq. 0 .and. longoutput) write(6,"(t5,'center ', i8,/t16,'n', t20, 'l', t29,'exp', t39, 'ind func')") ia
            lmaxc(ia) = 0
            rlargo(ia) = cero
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                rnor(i) = sqrt((dos * xx(i))**(2*nn(i)+1) / fact(2*nn(i)))
                rlarex = (15.d0 + 2.5d0 * re(nn(i)-1) + log(rnor(i))) / xx(i)
                rlargo(ia) = max(rlargo(ia), rlarex)		! intended to accelerate calculations if lexact .eq. .true.
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (myrank .eq. 0 .and. longoutput)  write(6,"(t15,i2,t19,i2,t24,e12.5,t40,i4)") nn(i), ll(i), xx(i), nf(i)
            enddo
        enddo
        if (myrank .eq. 0) write(6,"('Number of basis functions = ', i4)") nbas
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)
!		Allocates memory for the basis set

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
            lmaxc(ia) = 0
            rlargo(ia) = cero
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
!				computes and stores the radial normalization factor
                aux = cero
                bux = ll(knt) + 1.5d0
                do k1 = 1, nprimit(knt)
                    do k2 = 1, k1-1
                        aux=aux + dos*cfcontr(icarga+k1)*cfcontr(icarga+k2)/(xxg(icarga+k1)+xxg(icarga+k2))**bux
                    enddo
                    aux = aux + cfcontr(icarga+k1) * cfcontr(icarga+k1) / (dos*xxg(icarga+k1))**bux
                enddo
                rnor(knt) = sqrt( dos / (facts(ll(knt))*aux) )
                rlarex = sqrt((15.d0 + 2.5d0 * re(ll(knt)) + log(rnor(knt))) / xxg(icarga+1))
                rlargo(ia) = max(rlargo(ia), rlarex)		! intended to accelerate calculations if lexact .eq. .true.
                icarga = icarga+nprimit(knt)	! actualizes the index for loading primitives exponents and contraction coefficients
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
        if (myrank .eq. 0) write(6,"('number of total primitives = ', i8)") nprimitot
    endif

    if (lexact) then
        if (myrank .eq. 0) then
            do ia = 1, ncen
                write(6,"('long range radius of basis set for atom ', i4, ': ', e12.5)") ia, rlargo(ia)
            enddo
        endif
        return
    endif

!	Data of density representation
    read(10) lmaxexp
    nsize = nsize - sizeof(lmaxexp)

    if (lmaxrep .gt. lmaxexp .and. .not. lexact) then
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
    nsize = nsize - (sizeof(icfposd(:,1)) - sizeof(xajustd(:,1))) * ncenbas

    allocate(cfajust(nsize/8), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfajust in processor ',i3)") myrank
        abort = 1
        return
    endif
    k = 0
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) cycle
        read(10) icfposd(1:lmtop*nintervaj+1,ia)
        if (k .gt. 0) icfposd(1:lmtop*nintervaj+1,ia) = icfposd(1:lmtop*nintervaj+1,ia) + icfposd(lmtop*nintervaj+1,k) - 1
        k = ia
        read(10) xajustd(1:nintervaj,ia)		! Exponents
        if (longoutput) write(6,"('fitting exponents: ',/, 8(1x,e17.10))")  xajustd(1:nintervaj,ia)
!     fitting coeficients
        read(10) cfajust(icfposd(1,ia):icfposd(lmtop*nintervaj+1,ia)-1)
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
    close(10)

!	Determines the long-range radii and the highest l in the expansion for each interval

    allocate(lcorto(nintervaj,ncen), umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lcorto and umedpow in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)
    umedpow(0) = uno							!
    do i = 1, lmaxexp							!
        umedpow(i) = umedpow(i-1) * umed			! 1 / 2^i
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
            do i = 0, nsamples-1	! samples over nsamples points in each interval to decide the highest l
                ra = rinterv(interv-1) + dltsample + (rinterv(interv) - rinterv(interv-1) - dos * dltsample) &
                        * ri(nsamples-1) * i
                aux = exp(-xajustd(interv,ia)*ra)
                t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
                dost = t + t
                tcheb(0) = uno	! Chebyshev T  polynomials
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
!	
!   ***************************************************************
!
   subroutine gridexac
    USE MPI
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ierr, iii, iuni, ix, iy, iz, knt, mpireal, nx, nxyz, nxyzrank, ny, nz
    real(KREAL) :: b2a, den, rx, ry, rz, x, y, z
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: array(:), arrayrank(:)
    mpireal = MPI_REAL8
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: array(:), arrayrank(:)
    mpireal = MPI_REAL4
#endif

    rx = (xsup - xinf) / dltx + umed
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + umed
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + umed
    nz = int(rz) + 1

    if (myrank .eq. 0) then
        write(6,"('GRID (inf,sup,dlt,npnts)')")
        write(6,"(3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
        write(6,"(3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
        write(6,"(3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz
    endif
    if (langstrom) then		! Converts grid distances to angstrom
        b2a = 0.5291772d0
        x1 = xinf * b2a
        y1 = yinf * b2a
        z1 = zinf * b2a
        x2 = (xinf+(nx-1)*dltx) * b2a
        y2 = (yinf+(ny-1)*dlty) * b2a
        z2 = (zinf+(nz-1)*dltz) * b2a
    else
        x1 = xinf
        y1 = yinf
        z1 = zinf
        x2 = (xinf+(nx-1)*dltx)
        y2 = (yinf+(ny-1)*dlty)
        z2 = (zinf+(nz-1)*dltz)
    endif
!	Determines the grid points for tabulation in each processor:
!		istav(myrank): starting index iz assigned to processor myrank
!		iendv(myrank): ending index iz assigned to processor myrank

    allocate(istav(0:nprocs-1), iendv(0:nprocs-1), ilenv(0:nprocs-1), idispv(0:nprocs-1), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating istav, iendv, ilenv and idispv in gridexac in processor ',i3)") myrank
        abort = 1
        return
    endif

    call para_range(nz)
    idispv(0) = 0
    do i = 0, nprocs-2
            ilenv(i) = (iendv(i)-istav(i)+1) * nx * ny
            idispv(i+1) = idispv(i) + ilenv(i)
    enddo
    ilenv(nprocs-1) = (iendv(nprocs-1)-istav(nprocs-1)+1) * nx * ny
!	Allocates arrays
    nxyz = nx*ny*nz
    if (myrank .eq. 0) then
        allocate(array(nxyz), stat = ierr)
        if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraysp in gridexact in processor ',i3)") myrank
                abort = 1
                return
        endif
        iuni = 1
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"_exact-d.plt")
        array = cero	! Array initialization
    endif
    nxyzrank = ilenv(myrank)

    allocate( arrayrank(nxyzrank), lnegia(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating arrayrank and lnegia in gridexac in processor ',i3)") myrank
        abort = 1
        return
    endif
    idimzlm = (mxl+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma in gridexac in processor ',i3)") myrank
        abort = 1
        return
    endif
! 	arrayranksp = cero
    knt = 0
    do iz =  1, nz
        if (mod(iz-1,nprocs) .ne. myrank) cycle
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix = 1, nx
                knt = knt + 1
                x = xinf + (ix - 1) * dltx
                if (lsto) then
                        call densexactaSTO(x, y, z, den)
                else
                        call densexactaGTO(x, y, z, den)
                endif
                arrayrank(knt) = den
            enddo
        enddo
    enddo
    CALL MPI_GATHERV(arrayrank, nxyzrank, mpireal, array, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0 .and. myrank .eq. 0) then
        write(6,"('Error in MPI_GATHERV for arraysp,   ierr = ',i7)") ierr
        abort = 1
        return
    endif
    if (myrank .eq. 0) then
        do iz =  1, nz
            knt = idispv(mod(iz-1,nprocs)) + nx * ny * ((iz-1)/nprocs)
            do iy = 1, ny
                do ix =1, nx
                    knt = knt + 1
                    arrayrank(ix) = array(knt)
                enddo
                call linea( 1, nx , arrayrank )
            enddo
        enddo
        close(1)	! 	Closes the grid file
        write(6,"('Total number of tabulated points = ', i12)") nxyz
    endif
    if (myrank .eq. 0) deallocate(array)
    deallocate(arrayrank, zlma)
    return
    end
!	
!   ***************************************************************
!
   subroutine gridrep
    USE MPI
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, ierr, isel, iuni, ix, iy, iz, knt, mpireal, nx, ny, nz, nxyz, nxyzrank
    real(KREAL) :: b2a, denrep, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz, dV, rx, ry, rz, x, y, z
    character*4 :: strbux
    character*256 :: straux
    character*7 :: strcux
    real(KREAL), allocatable :: arrayaccrank(:), arrayaccdxrank(:), arrayaccdyrank(:), arrayaccdzrank(:)
    real(KREAL), allocatable :: arrayatrank(:,:), arrayatdxrank(:,:), arrayatdyrank(:,:), arrayatdzrank(:,:)
    real(KREAL), allocatable :: arraydxrank(:), arraydyrank(:), arraydzrank(:)
    real(KREAL), allocatable :: arraydxxrank(:), arraydxyrank(:), arraydxzrank(:), arraydyyrank(:)
    real(KREAL), allocatable :: arraydyzrank(:), arraydzzrank(:), arraylplrank(:), arrayrank(:)
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: arrayranksp(:)
    real(KREAL), allocatable :: array(:), arrayacc(:), arrayaccdx(:), arrayaccdy(:), arrayaccdz(:)
    real(KREAL), allocatable :: arrayat(:,:), arrayatdx(:,:), arrayatdy(:,:), arrayatdz(:,:)
    real(KREAL), allocatable :: arraydx(:), arraydy(:), arraydz(:)
    real(KREAL), allocatable :: arraydxx(:), arraydxy(:), arraydxz(:), arraydyy(:), arraydyz(:), arraydzz(:), arraylpl(:)
    mpireal = MPI_REAL8
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: arrayranksp(:)
    real(KREAL4), allocatable :: array(:), arrayacc(:), arrayaccdx(:), arrayaccdy(:), arrayaccdz(:)
    real(KREAL4), allocatable :: arrayat(:,:), arrayatdx(:,:), arrayatdy(:,:), arrayatdz(:,:)
    real(KREAL4), allocatable :: arraydx(:), arraydy(:), arraydz(:)
    real(KREAL4), allocatable :: arraydxx(:), arraydxy(:), arraydxz(:), arraydyy(:), arraydyz(:), arraydzz(:), arraylpl(:)
    mpireal = MPI_REAL4
#endif

    rx = (xsup - xinf) / dltx + umed
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + umed
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + umed
    nz = int(rz) + 1
    if (myrank .eq. 0)  then
        write(6,"('GRID (inf,sup,dlt,npnts)')")
        write(6,"(3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
        write(6,"(3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
        write(6,"(3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz
    endif
    if (langstrom) then		! Converts grid distances to angstrom
        b2a = 0.5291772d0
        x1 = xinf * b2a
        y1 = yinf * b2a
        z1 = zinf * b2a
        x2 = (xinf+(nx-1)*dltx) * b2a
        y2 = (yinf+(ny-1)*dlty) * b2a
        z2 = (zinf+(nz-1)*dltz) * b2a
    else
        x1 = xinf
        y1 = yinf
        z1 = zinf
        x2 = (xinf+(nx-1)*dltx)
        y2 = (yinf+(ny-1)*dlty)
        z2 = (zinf+(nz-1)*dltz)
    endif
!	Determines the grid points for tabulation assigned to each processor:
!		istav(myrank): starting index iz assigned to processor myrank
!		iendv(myrank): ending index iz assigned to processor myrank

    allocate(istav(0:nprocs-1), iendv(0:nprocs-1), ilenv(0:nprocs-1), idispv(0:nprocs-1), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating istav, iendv, ilenv and idispv in gridrep in processor ',i3)") myrank
        abort = 1
        return
    endif

    call para_range(nz)
    idispv(0) = 0
    do i = 0, nprocs-2
        ilenv(i) = (iendv(i)-istav(i)+1) * nx * ny
        idispv(i+1) = idispv(i) + ilenv(i)
    enddo
    ilenv(nprocs-1) = (iendv(nprocs-1)-istav(nprocs-1)+1) * nx * ny
    if (myrank .eq. 0 .and. longoutput)  then
        write(6,"(/,'Grid points assignment')")
        do i = 0, nprocs-1
            write(6,"('In proc ', i3,' istart = ', i3, ' iend = ', i3, ' no points = ', i12, ' disp = ', i12)") &
                    i, istav(i), iendv(i), ilenv(i), idispv(i)
        enddo
    endif

!	Allocates arrays and opens files for molecular density tabulation

    nxyz = nx*ny*nz
    if (myrank .eq. 0) then
    if (ierr .ne. 0) call error(1,'Memory error when allocating arraysp in gridrep. Stop')
        if (ldeform) then
            strcux="-deform"
        else
            strcux=""
        endif
        if (lmolec) then
            allocate(array(nxyz), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating array in gridrep in processor ',i3)") myrank
                abort = 1
                return
            endif
            straux = trim(filename)//trim(strcux)
            iuni = 21
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d.plt")
            array = cero	! Array initialization
            if (lgradient) then
                allocate(arraydx(nxyz), arraydy(nxyz), arraydz(nxyz), stat = ierr)
                if (ierr .ne. 0) then
                    write(6,"('Memory error when allocating arraydx, arraydy and arraydz in gridrep in processor ',i3)") myrank
                    abort = 1
                    return
                endif
                iuni = 23
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dx.pltd")
                iuni = 24
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dy.pltd")
                iuni = 25
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dz.pltd")
                arraydx = cero	! Array initialization
                arraydy = cero	! Array initialization
                arraydz = cero	! Array initialization
            endif
            if (lderiv2) then
                allocate(arraydxx(nxyz), arraydxy(nxyz), arraydxz(nxyz), arraydyy(nxyz), arraydyz(nxyz) &
                        , arraydzz(nxyz), stat = ierr)
                if (ierr .ne. 0) then
                        write(6,"('Memory error when allocating arraydxx ... arraydzz in gridrep in processor ',i3)") myrank
                        abort = 1
                        return
                endif
                iuni = 28
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dxx.pltd")
                iuni = 29
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dxy.pltd")
                iuni = 30
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dxz.pltd")
                iuni = 31
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dyy.pltd")
                iuni = 32
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dyz.pltd")
                iuni = 33
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dzz.pltd")
                arraydxx = cero	! Array initialization
                arraydxy = cero	! Array initialization
                arraydxz = cero	! Array initialization
                arraydyy = cero	! Array initialization
                arraydyz = cero	! Array initialization
                arraydzz = cero	! Array initialization
            endif
            if (laplacian) then
                allocate(arraylpl(nxyz), stat = ierr)
                if (ierr .ne. 0) then
                    write(6,"('Memory error when allocating arraylpl in gridrep in processor ',i3)") myrank
                    abort = 1
                    return
                endif
                iuni = 27
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-lplc.plt")
                arraylpl = cero	! Array initialization
            endif
        endif
        if (ldensacc) then
            allocate(arrayacc(nxyz), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arrayacc in gridpot in processor ',i3)") myrank
                abort = 1
                return
            endif
            iuni = 22
            straux = trim(filename)//"-fgr"//trim(strcux)
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d.plt")
            arrayacc = cero	! Array initialization
            if (lgradient) then
                allocate(arrayaccdx(nxyz), arrayaccdy(nxyz), arrayaccdz(nxyz), stat = ierr)
                if (ierr .ne. 0) then
                    write(6,"('Memory error when allocating arrayaccdx, &
                            &arrayaccdy and arrayaccdz  in gridpot in processor ',i3)") myrank
                    abort = 1
                    return
                endif
                iuni = 34
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dx.pltd")
                iuni = 35
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dy.pltd")
                iuni = 36
                call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dz.pltd")
            endif
        endif
        if (latomics) then
            allocate(arrayat(nxyz,nsel), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arrayat in gridrep in processor ',i3)") myrank
                abort = 1
                return
            endif
            if (lgradient) then
                allocate(arrayatdx(nxyz,nsel), arrayatdy(nxyz,nsel), arrayatdz(nxyz,nsel), stat = ierr)
                if (ierr .ne. 0) then
                    write(6,"('Memory error when allocating arrayatdx, arrayatdy, arrayatdz  in gridrep in processor ',i3)") &
                            myrank
                    abort = 1
                    return
                endif
                iuni = 40
                do i = 1, nsel
                    iuni = iuni + 1
                    write(strbux,'(i4.1)') iatomsel(i)
                    straux = trim(filename)//"-"//trim(atmnms(nzn(iatomsel(i))))//trim(adjustl(strbux))//trim(strcux)
                    call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d.plt")
                    if (lgradient) then
                        iuni = iuni + 1
                        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dx.pltd")
                        iuni = iuni + 1
                        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dy.pltd")
                        iuni = iuni + 1
                        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dz.pltd")
                    endif
                enddo
            endif
        endif
    endif
    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma in gridrep in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (lgradient) then
        allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating zlmadx, zlmady, zlmadz in gridrep in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif
    if (lderiv2) then
        allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                zlmadzz(idimzlm), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating zlmadxx, zlmadxy,... zlmadzz in gridrep in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif
    nxyzrank = ilenv(myrank)
    allocate(arrayranksp(nxyzrank), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating arrayranksp in gridrep in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (lmolec) then
        allocate(arrayrank(nxyzrank), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arrayrank in gridrep in processor ',i3)") myrank
            abort = 1
            return
        endif
        arrayrank = cero	! Array initialization
        if (lgradient) then
            allocate(arraydxrank(nxyzrank), arraydyrank(nxyzrank), arraydzrank(nxyzrank), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraydxrank, arraydyrank and arraydzrank in gridrep in processor ' &
                        ,i3)") myrank
                abort = 1
                return
            endif
            arraydxrank = cero	! Array initialization
            arraydyrank = cero	! Array initialization
            arraydzrank = cero	! Array initialization
        endif
        if (lderiv2) then
            allocate(arraydxxrank(nxyzrank), arraydxyrank(nxyzrank), arraydxzrank(nxyzrank), arraydyyrank(nxyzrank) &
                    , arraydyzrank(nxyzrank), arraydzzrank(nxyzrank), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraydxxrank ... arraydzzrank in gridexact in processor ',i3)") myrank
                abort = 1
                return
            endif
            arraydxxrank = cero	! Array initialization
            arraydxyrank = cero	! Array initialization
            arraydxzrank = cero	! Array initialization
            arraydyyrank = cero	! Array initialization
            arraydyzrank = cero	! Array initialization
            arraydzzrank = cero	! Array initialization
        endif
        if (laplacian) then
            allocate(arraylplrank(nxyzrank), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraylplrank in gridexact in processor ',i3)") myrank
                abort = 1
                return
            endif
            arraylplrank = cero	! Array initialization
        endif
    endif
    if (ldensacc) then
        allocate(arrayaccrank(nxyzrank), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arrayaccrank in gridrep in processor ',i3)") myrank
            abort = 1
            return
        endif
        arrayaccrank = cero	! Array initialization
        if (lgradient) then
            allocate(arrayaccdxrank(nxyzrank), arrayaccdyrank(nxyzrank), arrayaccdzrank(nxyzrank), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arrayaccdxrank, arrayaccdyrank, arrayaccdzrank &
                        &in gridrep in processor ',i3)") myrank
                abort = 1
                return
            endif
            arrayaccdxrank = cero	! Array initialization
            arrayaccdyrank = cero	! Array initialization
            arrayaccdzrank = cero	! Array initialization
        endif
    endif
    if (latomics) then
        allocate(arrayatrank(nxyzrank,nsel), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arrayatrank in gridrep in processor ',i3)") myrank
            abort = 1
            return
        endif
        arrayatrank = cero	! Array initialization
        if (lgradient) then
            allocate(arrayatdxrank(nxyzrank,nsel), arrayatdyrank(nxyzrank,nsel), arrayatdzrank(nxyzrank,nsel), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arrayatdxrank, arrayatdyrank, arrayatdzrank  in gridrep &
                        &in processor ',i3)") myrank
                abort = 1
                return
            endif
            arrayatdxrank = cero	! Array initialization
            arrayatdyrank = cero	! Array initialization
            arrayatdzrank = cero	! Array initialization
        endif
    endif
! 	Grid tabulation
    knt = 0
    do iz =  1, nz
        if (mod(iz-1,nprocs) .ne. myrank) cycle
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix =1, nx
                knt = knt + 1
                x = xinf + (ix - 1) * dltx
                isel = 0
                do ia = 1, ncen
                    call densrepr(ia, x, y, z, denrep, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz)
                    if (lmolec) then
                        arrayrank(knt) = arrayrank(knt) + denrep
                        if (lgradient) then
                            arraydxrank(knt) = arraydxrank(knt) + dendrvx
                            arraydyrank(knt) = arraydyrank(knt) + dendrvy
                            arraydzrank(knt) = arraydzrank(knt) + dendrvz
                        endif
                        if (lderiv2) then
                            arraydxxrank(knt) = arraydxxrank(knt) + dxx
                            arraydxyrank(knt) = arraydxyrank(knt) + dxy
                            arraydxzrank(knt) = arraydxzrank(knt) + dxz
                            arraydyyrank(knt) = arraydyyrank(knt) + dyy
                            arraydyzrank(knt) = arraydyzrank(knt) + dyz
                            arraydzzrank(knt) = arraydzzrank(knt) + dzz
                        endif
                        if (laplacian) arraylplrank(knt) = arraylplrank(knt) + denlplc
                    endif
                    if (ldensacc .and. lselat(ia)) then
                        arrayaccrank(knt) = arrayaccrank(knt) + denrep
                        if (lgradient) then
                            arrayaccdxrank(knt) = arrayaccdxrank(knt) + dendrvx
                            arrayaccdyrank(knt) = arrayaccdyrank(knt) + dendrvy
                            arrayaccdzrank(knt) = arrayaccdzrank(knt) + dendrvz
                        endif
                    endif
                    if (latomics .and. lselat(ia)) then
                        isel = isel+1
                        arrayatrank(knt,isel) = denrep
                        if (lgradient) then
                            arrayatdxrank(knt,isel) = dendrvx
                            arrayatdyrank(knt,isel) = dendrvy
                            arrayatdzrank(knt,isel) = dendrvz
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo

    if (lmolec) then
        arrayranksp = arrayrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, array, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error in MPI_GATHERV for array,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        if (lgradient) then
            arrayranksp = arraydxrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydx, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydx,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydyrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydy,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydzrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydz,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
        endif
        if (lderiv2) then
            arrayranksp = arraydxxrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydxx, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydxx,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydxyrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydxy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydxy,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydxzrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydxz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydxz,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydyyrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydyy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydyy,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydyzrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydyz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydyz,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arraydzzrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydzz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraydzz,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
        endif
        if (laplacian) then
            arrayranksp = arraylplrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraylpl, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arraylpl,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
        endif
    endif
    if (ldensacc) then
        arrayranksp = arrayaccrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayacc, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error in MPI_GATHERV for arrayacc,   ierr = ',i7)") ierr
            abort = 1
            return
! 			write(6,"('array = ', 8(1x,e17.10))") (array(iii), iii = 1, nx)
        endif
        if (lgradient) then
            arrayranksp = arrayaccdxrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayaccdx, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arrayaccdx,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arrayaccdyrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayaccdy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arrayaccdy,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            arrayranksp = arrayaccdzrank
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayaccdz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arrayaccdz,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
        endif
    endif
    if (latomics) then
        do i = 1, nsel
            arrayranksp = arrayatrank(:,i)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayat(1,i), ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error in MPI_GATHERV for arrayat,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (lgradient) then
                arrayranksp = arrayatdxrank(:,i)
                CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayatdx(1,i), ilenv, &
                        idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
                if (ierr .ne. 0 .and. myrank .eq. 0) then
                    write(6,"('Error in MPI_GATHERV for arrayatdx,   ierr = ',i7)") ierr
                    abort = 1
                    return
                endif
                arrayranksp = arrayatdyrank(:,i)
                CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayatdy(1,i), ilenv, &
                        idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
                if (ierr .ne. 0 .and. myrank .eq. 0) then
                    write(6,"('Error in MPI_GATHERV for arrayatdy,   ierr = ',i7)") ierr
                    abort = 1
                    return
                endif
                arrayranksp = arrayatdzrank(:,i)
                CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arrayatdz(1,i), ilenv, &
                        idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
                if (ierr .ne. 0 .and. myrank .eq. 0) then
                    write(6,"('Error in MPI_GATHERV for arrayatdz,   ierr = ',i7)") ierr
                    abort = 1
                    return
                endif
            endif
        enddo
    endif

    if (myrank .eq. 0) then
        do iz =  1, nz
            knt = idispv(mod(iz-1,nprocs)) + nx * ny * ((iz-1)/nprocs)
            do iy = 1, ny
                do ix =1, nx
                    knt = knt + 1
                    if (lmolec) then
                        arrayrank(ix) = array(knt)
                        if (lgradient) then
                            arraydxrank(ix) = arraydx(knt)
                            arraydyrank(ix) = arraydy(knt)
                            arraydzrank(ix) = arraydz(knt)
                        endif
                        if (lderiv2) then
                            arraydxxrank(ix) = arraydxx(knt)
                            arraydxyrank(ix) = arraydxy(knt)
                            arraydxzrank(ix) = arraydxz(knt)
                            arraydyyrank(ix) = arraydyy(knt)
                            arraydyzrank(ix) = arraydyz(knt)
                            arraydzzrank(ix) = arraydzz(knt)
                        endif
                        if (laplacian ) arraylplrank(ix) = arraylpl(knt)
                    endif
                    if (ldensacc) then
                        arrayaccrank(ix) = arrayacc(knt)
                        if (lgradient) then
                            arrayaccdxrank(ix) = arrayaccdx(knt)
                            arrayaccdyrank(ix) = arrayaccdy(knt)
                            arrayaccdzrank(ix) = arrayaccdz(knt)
                        endif
                    endif
                    if (latomics ) then
                        do i = 1, nsel
                            arrayatrank(ix,i) = arrayat(knt,i)
                            if (lgradient) then
                                arrayatdxrank(ix,i) = arrayatdx(knt,i)
                                arrayatdyrank(ix,i) = arrayatdy(knt,i)
                                arrayatdzrank(ix,i) = arrayatdz(knt,i)
                            endif
                        enddo
                    endif
                enddo
                if (lmolec) then
                    arrayranksp(1:nx) = arrayrank(1:nx)
                    call linea( 21, nx , arrayranksp )
                    if (lgradient) then
                        arrayranksp(1:nx) = arraydxrank(1:nx)
                        call linea(23, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydyrank(1:nx)
                        call linea(24, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydzrank(1:nx)
                        call linea(25, nx , arrayranksp )
                    endif
                    if (lderiv2) then
                        arrayranksp(1:nx) = arraydxxrank(1:nx)
                        call linea(28, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydxyrank(1:nx)
                        call linea(29, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydxzrank(1:nx)
                        call linea(30, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydyyrank(1:nx)
                        call linea(31, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydyzrank(1:nx)
                        call linea(32, nx , arrayranksp )
                        arrayranksp(1:nx) = arraydzzrank(1:nx)
                        call linea(33, nx , arrayranksp )
                    endif
                    if (laplacian ) then
                        arrayranksp(1:nx) = arraylplrank(1:nx)
                        call linea(27, nx , arrayranksp )
                    endif
                endif
                if (ldensacc) then
                    arrayranksp(1:nx) = arrayaccrank(1:nx)
                    call linea( 22, nx , arrayranksp )
                    if (lgradient) then
                        arrayranksp(1:nx) = arrayaccdxrank(1:nx)
                        call linea(34, nx , arrayranksp )
                        arrayranksp(1:nx) = arrayaccdyrank(1:nx)
                        call linea(35, nx , arrayranksp )
                        arrayranksp(1:nx) = arrayaccdzrank(1:nx)
                        call linea(36, nx , arrayranksp )
                    endif
                endif
                if (latomics ) then
                    iuni = 40
                    do i = 1, nsel
                        iuni = iuni+1
                        arrayranksp(1:nx) = arrayatrank(1:nx,i)
                        call linea( iuni, nx , arrayranksp )
                        if (lgradient) then
                            iuni = iuni+1
                            arrayranksp(1:nx) = arrayatdxrank(1:nx,i)
                            call linea( iuni, nx , arrayranksp )
                            iuni = iuni+1
                            arrayranksp(1:nx) = arrayatdyrank(1:nx,i)
                            call linea( iuni, nx , arrayranksp )
                            iuni = iuni+1
                            arrayranksp(1:nx) = arrayatdzrank(1:nx,i)
                            call linea( iuni, nx , arrayranksp )
                        endif
                    enddo
                endif
            enddo
        enddo
! 	Deallocates arrays and closes the grid files
        if (lmolec) then
            close(21)
            deallocate(array)
            if (lgradient) then
                close(23)
                close(24)
                close(25)
                deallocate(arraydx, arraydy,  arraydz)
            endif
            if (lderiv2) then
                close(28)
                close(29)
                close(30)
                close(31)
                close(32)
                close(33)
                deallocate(arraydxx, arraydxy, arraydxz, arraydyy, arraydyz, arraydzz)
            endif
            if (laplacian ) then
                close(27)
                deallocate(arraylpl)
            endif
        endif
        if (ldensacc) then
            close(22)
            deallocate(arrayacc)
            if (lgradient) then
                close(34)
                close(35)
                close(36)
                deallocate(arrayaccdx, arrayaccdy,  arrayaccdz)
            endif
        endif
        if (latomics) then
            deallocate(arrayat)
            if (lgradient) deallocate(arrayatdx, arrayatdy, arrayatdz)
            iuni = 40
            do i = 1, nsel
                iuni = iuni + 1
                close(iuni)
            enddo
        endif
        write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    endif
    if (lmolec) then
        deallocate(arrayrank)
        if (lgradient) then
            deallocate(arraydxrank, arraydyrank,  arraydzrank)
        endif
        if (lderiv2) deallocate(arraydxxrank, arraydxyrank, arraydxzrank, arraydyyrank, arraydyzrank, arraydzzrank)
        if (laplacian ) deallocate(arraylplrank)
    endif
    if (ldensacc) then
        deallocate(arrayaccrank)
        if (lgradient) deallocate(arrayaccdxrank, arrayaccdyrank,  arrayaccdzrank)
    endif
    if (latomics) then
        deallocate(arrayatrank)
        if (lgradient) deallocate(arrayatdxrank, arrayatdyrank, arrayatdzrank)
    endif
    deallocate(arrayranksp)
    deallocate(zlma)
    if (allocated(zlmadx)) deallocate(zlmadx, zlmady, zlmadz)
    if (allocated(zlmadxx)) deallocate(zlmadxx, zlmadxy, zlmadxz, zlmadyy, zlmadyz, zlmadzz)
    return
    end
!	
!   ***************************************************************
!
!   Density tabulation in point (x,y,z) from density matrix and basis set	
   subroutine densexactaSTO(x, y, z, denex)
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, i1, ia, ib, ierr, j, knt, kntb, l, lm, m
    real(KREAL) :: aux, denex, frad, ra, ra2, x, xxa, y, yya, z, zza
    real(KREAL) :: rapow(0:mxn)
    logical :: ldenex
    allocate(ang((mxl+1)*(mxl+2)/2), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ang in densexactaSTO in processor ',i3)") myrank
        abort = 1
        return
    endif
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
    mxind = (mxl+1)*(mxl+2)/2
    allocate(ind(0:mxind), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ind in densexactaSTO in processor ',i3)") myrank
        abort = 1
        return
    endif
    ind(0) = 0
    do i = 1, mxind
        ind(i) = ind(i-1) + i         !  i*(i+1)/2
    enddo
    knt = 0
    denex = cero
    ldenex = .false.
    do ia = 1, ncen
        lnegia(ia) = .true.
        if (ngini(ia) .le. 0) cycle	! If center without basis set, skips to next atom
        xxa = x - rcen(1,ia)
        yya = y - rcen(2,ia)
        zza = z - rcen(3,ia)
        ra2 = xxa*xxa+yya*yya+zza*zza
        if (ra2 .gt. rlargo(ia)*rlargo(ia)) cycle	! Contributions of atom ia are negligible. Skips to next atom
        ldenex = .true.
        lnegia(ia) = .false.
        ra = sqrt(ra2)
        zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
        zlma(2) = yya
        zlma(3) = zza
        zlma(4) = xxa
        do l = 1, lmaxc(ia)
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xxa * zlma(l*(l+2)+1) - yya * zlma(l*l+1))		! zlm(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (yya * zlma(l*(l+2)+1) + xxa * zlma(l*l+1))		! zlm(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * zza* zlma(l*(l+2)+1)				! zlm(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * zza * zlma(l*l+1)					! zlm(l+1,-l,ia)
            do m = 0, l-1
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlm(l+1,m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlm(l+1,-m,ia)
            enddo
        enddo
        rapow(0) = uno
        do i = 1, mxn
            rapow(i) = rapow(i-1) * ra
        enddo
        do i1 = ngini(ia), ngfin(ia)
            frad = rnor(i1) * exp(-xx(i1)*ra) * rapow(nn(i1)-ll(i1)-1)
            do m = -ll(i1), ll(i1)
                knt = knt + 1
                f(knt) = frad * ang(ind(ll(i1))+abs(m)+1) * zlma(ll(i1)*(ll(i1)+1)+m+1)
            enddo
        enddo
    enddo
    deallocate (ang, ind)
    if (.not. ldenex) return
    knt = 0
    faux = cero
    do ia = 1, ncen
        if (lnegia(ia)) cycle		! Contributions of atom ia are negligible. Skips to next atom
        kntb = 0
        do ib = 1, ncen
            if (lnegia(ib)) cycle	! Contributions of atom ib are negligible. Skips to next atom
            do j = nf(ngini(ib)), nf(ngfin(ib))+2*ll(ngfin(ib))
                faux(kntb+j-nf(ngini(ib))+1) = faux(kntb+j-nf(ngini(ib))+1) + &
                        dot_product(f(knt+1:knt+nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1), &
                                        dmat(nf(ngini(ia)):nf(ngfin(ia))+2*ll(ngfin(ia)), j))
            enddo
            kntb = kntb -nf(ngini(ib))+ nf(ngfin(ib))+2*ll(ngfin(ib)) + 1
        enddo
        knt = knt + nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1
    enddo
    denex = max(dot_product(faux(1:knt),f(1:knt)),cero)	! Total density is positive defined
    return
    end
!	
!   ***************************************************************
!
!   Density tabulation in point (x,y,z) from density matrix and basis set	
   subroutine densexactaGTO(x, y, z, denex)
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, i1, i1p, ia, ib, ierr, j, knt, kntb, l, lm, m
    real(KREAL) :: aux, denex, frad, ra, ra2, x, xxa, y, yya, z, zza
    logical :: ldenex
    allocate(ang((mxl+1)*(mxl+2)/2), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ang in densexactaGTO in processor ',i3)") myrank
        abort = 1
        return
    endif
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
    mxind = (mxl+1)*(mxl+2)/2
    allocate(ind(0:mxind), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ind in densexactaGTO in processor ',i3)") myrank
        abort = 1
        return
    endif
    ind(0) = 0
    do i = 1, mxind
        ind(i) = ind(i-1) + i         !  i*(i+1)/2
    enddo
    knt = 0
    denex = cero
    ldenex = .false.
    do ia = 1, ncen
        lnegia(ia) = .true.
        if (ngini(ia) .le. 0) cycle	! If center without basis set, skips to next atom
        xxa = x - rcen(1,ia)
        yya = y - rcen(2,ia)
        zza = z - rcen(3,ia)
        ra2 = xxa*xxa+yya*yya+zza*zza
        if (ra2 .gt. rlargo(ia)*rlargo(ia)) cycle	! Contributions of atom ia are negligible. Skips to next atom
        ldenex = .true.
        lnegia(ia) = .false.
        ra = sqrt(ra2)
!		Computes the regular spherical harmonics of r-R(ia)
        zlma(1) = uno
        zlma(2) = yya
        zlma(3) = zza
        zlma(4) = xxa
        do l = 1, lmaxc(ia)
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xxa * zlma(l*(l+2)+1) - yya * zlma(l*l+1))		! zlm(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (yya * zlma(l*(l+2)+1) + xxa * zlma(l*l+1))		! zlm(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * zza* zlma(l*(l+2)+1)				! zlm(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * zza * zlma(l*l+1)					! zlm(l+1,-l,ia)
            do m = 0, l-1
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlm(l+1,m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlm(l+1,-m,ia)
            enddo
        enddo
        do i1 = ngini(ia), ngfin(ia)
            frad = cero
            do i1p = ipntprim(i1), ipntprim(i1)+nprimit(i1)-1
                frad = frad + cfcontr(i1p) * exp(-xxg(i1p)*ra2)
            enddo
            frad = rnor(i1) * frad
            do m = -ll(i1), ll(i1)
                knt = knt + 1
                f(knt) = frad * ang(ind(ll(i1))+abs(m)+1) * zlma(ll(i1)*(ll(i1)+1)+m+1)
            enddo
        enddo
    enddo
    deallocate (ang, ind)
    if (.not. ldenex) return
    knt = 0
    faux = cero
    do ia = 1, ncen
        if (lnegia(ia)) cycle		! Contributions of atom ia are negligible. Skips to next atom
        kntb = 0
        do ib = 1, ncen
            if (lnegia(ib)) cycle	! Contributions of atom ib are negligible. Skips to next atom
            do j = nf(ngini(ib)), nf(ngfin(ib))+2*ll(ngfin(ib))
                faux(kntb+j-nf(ngini(ib))+1) = faux(kntb+j-nf(ngini(ib))+1) + &
                        dot_product(f(knt+1:knt+nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1), &
                                        dmat(nf(ngini(ia)):nf(ngfin(ia))+2*ll(ngfin(ia)), j))
            enddo
            kntb = kntb -nf(ngini(ib))+ nf(ngfin(ib))+2*ll(ngfin(ib)) + 1
        enddo
        knt = knt + nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1
    enddo
    denex = max(dot_product(faux(1:knt),f(1:knt)),cero)	! Total density is positive defined
    return
    end
	
!   ***************************************************************

  subroutine densrepr(ia, x, y, z, denrep, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz)
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, icflm, interv, j, jshft, kntlm, l, m
    real(KREAL) :: aux, bux, cux, denrep, denlplc, dendrvx, dendrvy, dendrvz, dost, dux, drvflm, drv2flm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, eux, flm, fux
    real(KREAL) :: r3inv, ra, ra2, rainv, rj2, sgn, t, umt2i, x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL) :: tcheb(0:mxlenpol-1), ucheb(0:mxlenpol-1), drvtcheb(0:mxlenpol-1), drv2tcheb(0:mxlenpol-1)
!	Contribution of atomic fragment ia to density  in point (x,y,z)
    if (ngini(ia) .le. 0) then
        denrep = cero
        if (lgradient) then
            dendrvx = cero
            dendrvy = cero
            dendrvz = cero
        endif
        if (lderiv2) then
            dxx = cero
            dxy = cero
            dxz = cero
            dyy = cero
            dyz = cero
            dzz = cero
        endif
        if (laplacian) denlplc = cero
        return
    endif
    xa = x - rcen(1,ia)
    ya = y - rcen(2,ia)
    za = z - rcen(3,ia)
    ra2 = xa*xa+ya*ya+za*za
    if (ra2 .gt. rlargo(ia)*rlargo(ia)) then
        denrep = cero
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
    zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
    zlma(2) = ya
    zlma(3) = za
    zlma(4) = xa
    do l = 1, lcorto(interv,ia)
        zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))		! zlm(l+1,l+1,ia)
        zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))		! zlm(l+1,-(l+1),ia)
        zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)				! zlm(l+1,l,ia)
        zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)					! zlm(l+1,-l,ia)
        do m = 0, l-1
            zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlm(l+1,m,ia)
            zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlm(l+1,-m,ia)
        enddo
    enddo
    if (lgradient .or. laplacian) then
        if (lgradient) then
            dendrvx = cero
            dendrvy = cero
            dendrvz = cero
            call derivzlm(lcorto(interv,ia), idimzlm, zlma, zlmadx, zlmady, zlmadz)
        endif
        if (lderiv2) then
            dxx = cero
            dxy = cero
            dxz = cero
            dyy = cero
            dyz = cero
            dzz = cero
            call derivzlm(lcorto(interv,ia), idimzlm, zlmadx, zlmadxx, zlmadxy, zlmadxz)
            call dzlm2y(lcorto(interv,ia), idimzlm, zlmady, zlmadyy, zlmadyz)
            call dzlm2z(lcorto(interv,ia), idimzlm, zlmadz, zlmadzz)
        endif
        aux = exp(-xajustd(interv,ia)*ra)
        bux = xajustd(interv,ia)
        cux = dos / (rinterv(interv)-rinterv(interv-1))
        t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
        dost = t + t
        umt2i = uno / (t*t - uno)
        tcheb(0) = uno		! Chebyshev T and U polynomials (U polynomial are used for derivatives of T polynomials
        tcheb(1) = t		! and first and second derivatives of T polynomials according to:
        ucheb(0) = uno		!			D[T_n(t),t] = n * U_(n-1)(t)
        ucheb(1) = dost	!	and:		D[T_n(t),{t,2}] = (n/(t^2-1)) * ( (n-1)*t*U_(n-1) - n * U_(n-2) )
        drvtcheb(0) = cero
        drvtcheb(1) = uno
        drv2tcheb(0) = cero
        drv2tcheb(1) = cero
        sgn = uno
        do j = 2, mxlenpol-1
            tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
            ucheb(j) = dost * ucheb(j-1) - ucheb(j-2)
            drvtcheb(j) = re(j) * ucheb(j-1)
            if (uno-abs(t) .gt. 1.d-7) then
                drv2tcheb(j) = re(j) * umt2i * (re(j-1) * t * ucheb(j-1) - re(j) * ucheb(j-2))
            else	! For values of t very close to 1 or -1, takes a linear approximation (Taylor series) of D[T[j,t],{t,2}]
                rj2 = re(j) * re(j)
                drv2tcheb(j) = sgn * ri(3) * (rj2 * (rj2-uno) + ri(5) * rj2 * (re(4)+rj2*(-re(5)+rj2)) * (abs(t)-uno) )
                if (t .lt. cero) sgn = - sgn
            endif
        enddo

        denrep = cero
        denlplc = cero
        kntlm = lminrep*lminrep
        do l = lminrep, lcorto(interv,ia)	!     Computes density terms lminrep <= l <= lcorto(interv,ia)
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
                    drvflm = cux * aux * drvflm		! Converts the derivative with respect to t to derivative respect to r
                    drv2flm = cux * cux * aux * drv2flm	! Same for second derivative
                    drv2flm = bux * (bux * flm - (drvflm+drvflm)) + drv2flm	! D[flm,{r,2}]
                    drvflm = -bux * flm + drvflm							! D[flm,r]
                    denrep = denrep + flm * zlma(kntlm)
                    if (lgradient) then
                        dendrvx = dendrvx + xadivra * drvflm * zlma(kntlm) + flm * zlmadx(kntlm)
                        dendrvy = dendrvy + yadivra * drvflm * zlma(kntlm) + flm * zlmady(kntlm)
                        dendrvz = dendrvz + zadivra * drvflm * zlma(kntlm) + flm * zlmadz(kntlm)
                    endif
                    if (lderiv2) then
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
                    if (laplacian) denlplc = denlplc + zlma(kntlm) * (re(l+l+2) * rainv * drvflm + drv2flm)
                endif
            enddo
        enddo
    else
        aux = exp(-xajustd(interv,ia)*ra)
        t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
        dost = t + t
        tcheb(0) = uno	! Chebyshev T  polynomials
        tcheb(1) = t
        do j = 2, mxlenpol-1
            tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
        enddo
        denrep = cero
        kntlm = lminrep*lminrep
        do l = lminrep, lcorto(interv,ia)	!     Computes density terms lminrep <= l <= lcorto(interv,ia)
            do m = -l, l
                kntlm = kntlm + 1
                icflm = icfposd(kntlm+(interv-1)*lmtop,ia)
                if(icflm .lt. icfposd(kntlm+(interv-1)*lmtop+1,ia)) then
                    flm = cero
                    do j = 0, icfposd(kntlm+(interv-1)*lmtop+1,ia)-icflm-1
                        flm = flm + cfajust(j+icflm) * tcheb(j)
                    enddo
                    denrep = denrep + aux * flm * zlma(l*(l+1)+m+1)
                endif
            enddo
        enddo
    endif

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
!	***************************************************************
!	Subroutine cabecera: writes head for .plt files (binary)
!
  subroutine cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, s)
    USE DAM320_D
    implicit none
    integer(KINT) :: i
    integer(KINT) :: nx,ny,nz,iuni,ns,iaux(0:2)
    character*(*) :: s
#ifdef DBLPRCGRID
    real(KREAL) :: x1,x2,y1,y2,z1,z2,v(0:5)
    integer(KINT), parameter :: i0 = 0
#else
    real(KREAL4) :: x1,x2,y1,y2,z1,z2,v(0:5)
    integer(KINT), parameter :: i0 = 3
#endif
!	If the compiler is other than INTEL's, uses the OPEN
!	sentence for stream files according to Fortran 2003 standard
#if _WIN32
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#elif __INTEL_COMPILER
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#else
    open (unit=iuni, file=s, form='unformatted', access='stream')
#endif
    iaux(0) = i0	! 3 for single precision grid and gOpenMol compatibility, 0 for double precision grid (no gOpenMol compatibility)
    iaux(1) = 2
    write(iuni) iaux(0), iaux(1)
    iaux(0) = nz
    iaux(1) = ny
    iaux(2) = nx
    write(iuni) iaux(0), iaux(1), iaux(2)
    v(0) = z1
    v(1) = z2
    v(2) = y1
    v(3) = y2
    v(4) = x1
    v(5) = x2
    write(iuni) (v(i), i = 0, 5)
    END
!
!	Subroutine linea: writes lines of .plt files (binary)
!
  SUBROUTINE linea(iuni,n,v)
    USE DAM320_D
    implicit none
    integer(KINT) :: i, iuni, n
#ifdef DBLPRCGRID
    real(KREAL) :: v(0:n-1)
#else
    real(KREAL4) :: v(0:n-1)
#endif
    write(iuni) (v(i), i = 0, n-1)
    END
!
!	-------------------------------------------------------------------------------------------------------
!
  subroutine para_range(nz)
    USE MPI
    USE DAM320_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, n1, n2, nz
    integer(KINT) :: naux(0:nprocs)
    n1 = nz / nprocs
    n2 = mod(nz,nprocs)
    naux(0) = 1
    do i = 1, n2
            naux(i) = naux(i-1) + n1 + 1
    enddo
    do i = n2+1, nprocs
            naux(i) = naux(i-1) + n1
    enddo
    do i = 0, nprocs-1
            istav(i) = naux(i)
            iendv(i) = naux(i+1)-1
    enddo
    return
    end
!
!	-------------------------------------------------------------------------------------------------------
!
  subroutine error(ierr, msg)
    USE MPI
    USE DAM320_D
    implicit none
    integer(KINT) :: ierr, ierr2
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    call MPI_FINALIZE(ierr2)
    stop
    end
