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
!
! Version of September 2018
!

! #define DBLPRCGRID    ! Uncomment this line  if double precision grid is wanted
  program DAMDEN320
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDEN320_D
    implicit none
    integer(KINT) :: i, ia, ierr, ios, ipoint, j, lmaxctot, nbasis
    real(KREAL) :: aux, den, denlp, denrep, dendrvx, dendrvy, dendrvz, denlplc, dxden, dxx, dxxden
    real(KREAL) :: dxy, dxyden, dxz, dxzden, dyden, dyy, dyyden, dyz, dyzden, dzden, dzz, dzzden
    real(KREAL) :: x, xmax, xmin, xyzmax, xyzmin, y, ymax, ymin, z, zmax, zmin
    real(KREAL4) :: tarray(2), tiempo, dtime    
    namelist / options / dltu, dltv, dltx, dlty, dltz, filename, iatomsel, iswindows, langstrom, laplacian, latomics &
            , latomsel, lboundsx , lboundsy, lboundsz, ldeform, ldensacc, lderiv2, lexact, lgradient, lgrid, lgrid2d &
            , lmaxrep, lminrep, lmolec, longoutput, lpoints, nsel, numrtab, planeA, planeB, planeC, planecase &
            , rtab, umbrlargo, xinf, xsup, y_func_uv, yinf, ysup, z_func_uv, zinf, zsup &
            , xboundinf, xboundsup, yboundinf, uinf, usup, vinf, vsup, x_func_uv, yboundsup, zboundinf, zboundsup
    tiempo = dtime(tarray)

!    Namelist default values
    longoutput = .false.    ! If .true. a more detailed output is given
    langstrom = .true.        ! If .false. distances of grids (.plt files) in bohr
    ldeform = .false.           ! If .true. density deformations are computed (lminrep = 1)
    lexact  = .false.        ! if .true. "exact" density is tabulated
    lminrep = 0            ! lowest "l" in the expansion of the density. It does not hold if "lexact .eq. true"
    lmaxrep = 10                ! highest "l" in the expansion of the density. It does not hold if "lexact .eq. true"
                                            ! the expansion of the density takes lminrep <= l <= lmaxrep
    lmolec = .true.        ! If .true. tabulation of the whole molecular density
    latomics  = .false.        ! If .true. tabulation of atomic densities stored in files named:
                                            !    projectname-cxx-d.plt (gOpenMol)
                                            ! where xx refers to the center corresponding to the tabulation
    latomsel = .false.        ! If .true. indices of the centers for atomic tabulations will be user-supplied
                                            ! The indices of the selected atoms must be supplied in vector "iatomsel".
                                            ! Maximum number of centers that can be selected "mxsel" (parameter).
    nsel = 0                ! Number of centers for atomic tabulations
    lgradient = .false.        ! If .true. gradient components of the density computed and, if lgrid = .true., tabulated in files
                                            !    projectname-d-dx.pltd, projectname-d-dy.pltd, projectname-d-dz.pltd
    laplacian = .false.        ! If .true. laplacian computed tabulated in file projectname-dlplc.plt if lgrid = .true.
    lderiv2 = .false.        ! If .true. second derivatives of the density computed and, if lgrid = .true., tabulated in files
                                            !    projectname-d-dxx.pltd, projectname-d-dxy.pltd, projectname-d-dxz.pltd, projectname-d-dyz.pltd, etc
    iatomsel = 0            ! To cover the possibility of an input file with "latomsel = .true."
    iatomsel(1) = 1        ! but without assignment of "iatomsel".
    ldensacc  = .false.        ! If .true. computes and tabulates the accumulated density of the selected atoms in file projectname-frg-d.plt
    lboundsx = .false.        ! If .true. the evaluation of deformation charges is constrained to a range in x
    xboundinf = cero        ! Lower limit of that range
    xboundsup = cero        ! Upper limit of that range
    lboundsy = .false.        ! If .true. the evaluation of deformation charges is constrained to a range in y
    yboundinf = cero        ! Lower limit of that range
    yboundsup = cero        ! Upper limit of that range
    lboundsz = .false.        ! If .true. the evaluation of deformation charges is constrained to a range in z
    zboundinf = cero        ! Lower limit of that range
    zboundsup = cero        ! Upper limit of that range
    lgrid = .true.            ! If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    lgrid2d = .false.        ! If .true. computes a 2D grid. (x,y,z) are given in terms of (u,v)
    lpoints = .false.        ! If .true. computes in selected points and prints out the results in the standard output.
                                            ! Points must be given in cartesian coordinates.
                                            ! If lgrid .eq. .true., these coordinates must be placed after the grid data
    numrtab = 0            ! Number of tabulation points supplied in namelist
    rtab = cero            ! Tabulation points supplied in namelist
    umbrlargo = 1.d-8        ! Long-range threshold
    filename = ""            ! root file name for .plt and .pltd files
    iswindows = .false.        ! .true. if running on a MS-windows system
    xinf = cero
    xsup = cero
    dltx = uno
    yinf = cero
    ysup = cero
    dlty = uno
    zinf = cero
    zsup = cero
    dltz = uno

    x_func_uv = 'u'        ! x = u for 2D grids:  default: plane XY
    y_func_uv = 'v'        ! y = v for 2D grids
    z_func_uv = '0'        ! z = 0 for 2D grids
    uinf = cero
    usup = uno
    dltu = uno
    vinf = cero
    vsup = uno
    dltv = uno

    planeA = cero        ! Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    planeB = cero
    planeC = uno
    planecase = 1

!    End of namelist defaults

    planesuffix = ""
    read(5,OPTIONS)
    read(5,*) projectname
    call subplanesuffix(planecase,planesuffix)
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
    if (len_trim(filename).eq.0) then
        filename = projectname
    else
        filename = projectname(1:i)//trim(filename)
    endif
    if (lderiv2) lgradient = .true.

    call consta        ! Computes auxiliary constants

    ldengz = .false.    ! Checks whether the .den file is gzipped or not
    if (lexact) then
        ldeform = .false. ! Deformations computation is allowed only for represented density
        inquire(file=trim(projectname)//".den.gz", exist=ldengz, iostat=ierr)
        if (ierr .eq. 0 .and. ldengz) then
            call system ("gunzip "//trim(projectname)//".den.gz")
        endif
    endif

    call leedamqtden        ! Reads file .damqt  (generated by DAM2016)

    if (ldeform) lminrep = 1
    
#ifdef DBLPRCGRID
    write(6,"('Grid generated in double precision')")
    write(6,"('WARNING! This grid will not be compatible with gOpenMol')")
#endif

    allocate(lselat(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating lselat in main. Stop')
!     selected atomic fragments for tabulation    
    do ia = 1, ncen
        lselat(ia) = .false.
    enddo
    if (latomsel) then
        do i = 1, mxsel
            if (iatomsel(i) .le. ncen .and. iatomsel(i) .gt. 0) lselat(iatomsel(i)) = .true.
        enddo
        nsel = 0
        do ia = 1, ncen
            if (ngini(ia) .le. 0) cycle
            if (lselat(ia)) then
                    nsel = nsel+1
                    iatomsel(nsel) = ia
            endif
        enddo
    elseif (latomics ) then
        nsel = ncen
        latomsel = .true.
        if (ncen .gt. mxsel) then
            write(6,"('WARNING: Number of selected atoms larger than allowed for tabulation of the atomic densities')")
            write(6,"('Only the first ',i3, ' atomic densities will be tabulated')") mxsel
            nsel = mxsel
        endif
        do i = 1, nsel
            lselat(i) = .true.
            iatomsel(i) = i
        enddo
    endif

    if (lminrep .gt. lmaxexp .and. .not. lexact) then
        write(6,"('lminrep = ', i3, ' greater than lmaxexp ', i3)") lminrep, lmaxexp
        write(6,"('takes lminrep = 0')")
        lminrep = 0
    endif
    if (lexact) then
        lgradient = .false.    ! The gradient is not computed for the exact density
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
            call error(1,'No density matrix available. Stop')
        endif
!        Allocates auxiliary arrays for exact density

        allocate(f(nbas), faux(nbas), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating f and faux. Stop')

        if (ldengz) then    ! restores files back to their original gzipped status
            call system ("gzip "//trim(projectname)//".den")
        endif
        write(6,"('Density computed from basis set and density matrix')")
    else
        write(6,"('Expansion of the density: lminrep = ', i3, ' : lmaxrep = ', i3)") lminrep, lmaxrep
    endif
    if (lgrid) then
        if ((xsup-xinf) .eq. cero .or. (ysup-yinf) .eq. cero .or. (zsup-zinf) .eq. cero) then    !    Default grid
            xmin = minval(rcen(1,1:ncen))
            xmax = maxval(rcen(1,1:ncen))
            ymin = minval(rcen(2,1:ncen))
            ymax = maxval(rcen(2,1:ncen))
            zmin = minval(rcen(3,1:ncen))
            zmax = maxval(rcen(3,1:ncen))
            xyzmin = min(xmin, ymin, zmin)
            xyzmax = max(xmax, ymax, zmax)
            if ((xmin-xmax) .eq. cero) then
                xmin = min(xyzmin,-uno)
                xmax = max(xyzmax,uno)
            endif
            if ((ymin-ymax) .eq. cero) then
                ymin = min(xyzmin,-uno)
                ymax = max(xyzmax,uno)
            endif
            if ((zmin-zmax) .eq. cero) then
                zmin = min(xyzmin,-uno)
                zmax = max(xyzmax,uno)
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
!    Computes the density on the grid points
        if (.not. lgrid2d) then
            if (lexact) then
                call gridexac
            else
                call gridrep
            endif
        else
            if (lexact) then
                call gridexac2d
            else
                call gridrep2d
            endif
        endif
    endif
!    Tabulates specific points if required
    if (lpoints) then
        if (lexact) then
            allocate(lnegia(ncen), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating lnegia. Stop')
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
            open(9,file=trim(projectname)//"-d.der2",form='formatted')
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
            if (ierr .ne. 0) call error(1,'Memory error when allocating zlma. Stop')
            if (lgradient) then
                allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
                if (ierr .ne. 0) call error(1,'Memory error when allocating zlmadx, zlmady, zlmadz. Stop')
            endif
            if (lderiv2) then
                allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                        zlmadzz(idimzlm), stat = ierr)
                if (ierr .ne. 0) call error(1,'Memory error when allocating zlmadxx, zlmadxy,... zlmadzz. Stop')
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
            if (lexact) then
                 write(6,"(3(1x,e17.10), 3x, d22.15)") x, y, z, den
            elseif (lgradient) then
                if (lderiv2) then
                    write(9,"(/'Second derivatives of density at point ',3(1x,e17.10))") x, y, z
                    write(9,"(6(1x,e22.15))") dxxden, dxyden, dxzden, dyyden, dyzden, dzzden
                    if (laplacian) then
                        write(9,"('laplacian comprobation ',2(1x,e22.15))") denlp, dxxden + dyyden + dzzden
                    endif
                endif
                if (laplacian) then
                    write(6,"(3(1x,e17.10), 3x, d22.15, 4(1x,e22.15))") x, y, z, den, dxden, dyden, dzden, denlp
                else
!aux = uno / sqrt(dxden**2+dyden**2+dzden**2)
!dxden = dxden * aux
!dyden = dyden * aux
!dzden = dzden * aux
                    write(6,"(3(1x,e17.10), 3x, d22.15, 3(1x,e22.15))") x, y, z, den, dxden, dyden, dzden
                endif
            else
                if (laplacian) then
                    write(6,"(3(1x,e17.10), 3x, d22.15,1x, e22.15/)") x, y, z, den, denlp
                else
                    write(6,"(3(1x,e17.10), 3x, d22.15)") x, y, z, den
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

        if (allocated(lnegia)) deallocate(lnegia)
        if (allocated(zlma)) deallocate(zlma)
        if (allocated(zlmadx)) deallocate(zlmadx, zlmady, zlmadz)
        if (allocated(zlmadxx)) deallocate(zlmadxx, zlmadxy, zlmadxz, zlmadyy, zlmadyz, zlmadzz)
    endif
    tiempo = dtime(tarray)
    write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
            tarray(1), tarray(2), tarray(1)+tarray(2)
    stop
    end
!
!	***************************************************************
!
  subroutine readden
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, ierr, indnf, indng, ios, j, k, k1, k2, knt, nbasis, ncapserr
    real(KREAL) :: aux, bux

!	Allocates the array containing the density matrix

    allocate(dmat(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Memory error when allocating dmat. Stop')
    endif
    if (longoutput) write(6,"('Estimated highest size of dmat   = ', i15, ' bytes')") size(dmat)
    
    
    if (lsto) then
        open(16,file=trim(projectname)//".den",form='unformatted', iostat=ierr)
        if (ierr .ne. 0) call error(1,'Cannot open file '//trim(projectname)//'.den. Stop')
        read(16,iostat = ios) nbasis, ((dmat(i,j), i=1,nbasis), j=1,nbasis)
    else
        open(16,file=trim(projectname)//".den",form='formatted', iostat=ierr)
        if (ierr .ne. 0) call error(1,'Cannot open file '//trim(projectname)//'.den. Stop')
        read(16,*, iostat = ios) nbasis, ((dmat(i,j),j=1,i),i=1,nbasis)
    endif
    if ( ios .ne. 0 .or. nbas .ne. nbasis ) then
        write(6,"('ERROR reading density matrix')")
        write(6,"('Check whether the density matrix correspond to this basis set.')")
        call error(1,' Stop')
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
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
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
        call error(ierr,'Memory error when allocating dvec, ivec, jvec and dmat. Stop')
    endif

    do i=1,nbas*(nbas+1)/2
        read(16,end = 100, err = 100, iostat = ios) ivec(i), jvec(i), dvec(i)
    enddo

100 continue
    if (ios .eq. -1) then
        numdvec = i-1
        if (longoutput) write(6,"('Number of elements of density matrix read = ', i16)") numdvec
    else if ( ios .ne. 0 ) then
        call error(ios,'Error reading file '//trim(projectname)//".densprsbin"//'. Stop')
    endif
    close(16)
    
    if ( maxval(ivec) .gt. nbas .or. maxval(jvec) .gt. nbas ) then
        write(6,"('Index on density matrix higher than number of basis functions.',& 
            &/'Check whether the density matrix correspond to this basis set.')")
        write(6,*) 'maxval(ivec) = ', maxval(ivec)
        write(6,*) 'maxval(jvec) = ', maxval(jvec)
        write(6,*) 'nbas = ', nbas
        call error(1,'Stop.')
    endif
    
    write(6,"(/'Sparse density matrix read from binary file', a,/)") trim(projectname)//".densprsbin"
        
    dmat = cero
    do i = 1, numdvec
        dmat(ivec(i), jvec(i)) = dvec(i)
        dmat(jvec(i), ivec(i)) = dvec(i)
    enddo
    deallocate(dvec, ivec, jvec)
    return
    end
    
!**********************************************************************
!    subroutine consta
!
!    Computes and stores auxiliary constants
!        re(i) = dfloat(i)
!        ri(i) = 1.d0 / dfloat(i)
!        fact(i) = dfloat(i!)
!        facti(i) = 1.d0 / dfloat(i!)
!        facts(i) = dfloat((i+1/2)!)
!        ind(i) = i*(i+1)/2
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i
    real(KREAL) :: aux
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
        fact(i) = fact(i-1) * re(i)               !  i!
        facts(i) = facts(i-1) * re(i+i+1) * umed    ! (i+1/2)!
        facti(i) = uno / fact(i)                    !  uno / i!
    enddo                                !
    return
    end
!
!    ***************************************************************
!
  subroutine leedamqtden
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, icflm, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lm, m, ncenbas, ncfaj, ncflm, nsamples, nsize
    real(KREAL) :: aux, bux, dltsample, dost, flm, r, ra, ral, rlarex, step, suml, summ, t
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
    write(6,"('ncen = ', i8, ' nbas = ', i8, ' nshells = ', i8)") ncen, nbas, ncaps

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

    read(10) lsto    ! .true. means STO basis, .false. means GTO basis

    nsize = nsize - sizeof(lsto)
    if (lsto) then

!        Allocates memory for the basis set

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

        allocate(rnor(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating rnor. Stop')

        allocate(xx(ncaps), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating xx. Stop')

        write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        if (longoutput) write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp',t60,'rnor')")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfin(ia)
            nsize = nsize - sizeof(ngini(ia)) - sizeof(ngfin(ia))
            lmaxc(ia) = 0
            rlargo(ia) = cero
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            if (longoutput) write(6,"(t5,'center ', i8,/t16,'n', t20, 'l', t29,'exp', t39, 'ind func')") ia
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                rnor(i) = sqrt((dos * xx(i))**(2*nn(i)+1) / fact(2*nn(i)))
                rlarex = (15.d0 + 2.5d0 * re(nn(i)-1) + log(rnor(i))) / xx(i)
                rlargo(ia) = max(rlargo(ia), rlarex)        ! intended to accelerate calculations if lexact .eq. .true.
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (longoutput) write(6,"(t15,i2,t19,i2,t24,e12.5,t40,i4)") nn(i), ll(i), xx(i), nf(i)
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

        allocate(rlargo(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating rlargo. Stop')

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
                rlarex = sqrt((15.d0 + 2.5d0 * re(ll(knt)) + log(rnor(knt))) / xxg(icarga+1))
                rlargo(ia) = max(rlargo(ia), rlarex)        ! intended to accelerate calculations if lexact .eq. .true.
                icarga = icarga+nprimit(knt)    ! actualizes the index for loading primitives exponents and contraction coefficients
            enddo
        enddo

        write(6,"(/t22,'GTO Basis set',/t22,13('-'))")
        if (longoutput) then
            icarga = 0
            knt = 0       
            do ia = 1, ncen
                if (ncontr(ia) .le. 0) cycle
                write(6,"(/1x,'atom no.',1x,i8,'(',a2,')')") ia, atmnam(ia)
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
        write(6,"('Number of basis functions = ', i8)") nbas
        write(6,"('Total number of primitives = ', i8)") nprimitot
    endif
    if (lexact) then
        do ia = 1, ncen
                write(6,"('long range radius of basis set for atom ', i8, ': ', e12.5)") ia, rlargo(ia)
        enddo
        return    ! The density representation is not necessary, returns
    endif

!    Data of density representation
    read(10) lmaxexp

    nsize = nsize - sizeof(lmaxexp)
    if (lmaxrep .gt. lmaxexp .and. .not. lexact) then
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
    if (longoutput) write(6,"('Estimated highest size of xajust   = ', i15, ' bytes')") size(xajustd)
    nsize = nsize - sizeof(xajustd(:,1)) * ncenbas

    allocate(cfajust(nsize/8), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating cfajust. Stop')
    if (longoutput) write(6,"('Size of cfajust   = ', i15, ' bytes')") size(cfajust)

    if (longoutput) write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
    k = 0
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) cycle
        read(10) icfposd(1:lmtop*nintervaj+1,ia)
        if (k .gt. 0) icfposd(1:lmtop*nintervaj+1,ia) = icfposd(1:lmtop*nintervaj+1,ia) + icfposd(lmtop*nintervaj+1,k) - 1
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
!
!    ***************************************************************
!    Subroutine cabecera: writes head for .plt files (binary)
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
!    If the compiler is other than INTEL's, uses the OPEN
!    sentence for stream files according to Fortran 2003 standard
#if _WIN32
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#elif __INTEL_COMPILER
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#else
    open (unit=iuni, file=s, form='unformatted', access='stream')
#endif
    iaux(0) = i0    ! 3 for single precision grid and gOpenMol compatibility, 0 for double precision grid (no gOpenMol compatibility)
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
!    ***************************************************************
!    Subroutine cabecera: writes head for .plt files (binary)
!
  subroutine cabecera2d(nu, nv, u1, u2, v1, v2, iuni, s)
    USE DAM320_D
    implicit none
    integer(KINT) :: i
    integer(KINT) :: nu,nv,iuni,ns,iaux(0:1)
    character*(*) :: s
#ifdef DBLPRCGRID
    real(KREAL) :: u1,u2,v1,v2,v(0:3)
#else
    real(KREAL4) :: u1,u2,v1,v2,v(0:3)
#endif
!    If the compiler is other than INTEL's, uses the OPEN
!    sentence for stream files according to Fortran 2003 standard
#if _WIN32
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#elif __INTEL_COMPILER
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#else
    open (unit=iuni, file=s, form='unformatted', access='stream')
#endif
    iaux(0) = nu
    iaux(1) = nv
    write(iuni) iaux(0), iaux(1)
    v(0) = u1
    v(1) = u2
    v(2) = v1
    v(3) = v2
    write(iuni) (v(i), i = 0, 3)
    END
!
!    Subroutine linea: writes lines of .plt files (binary)
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
!   ***************************************************************
!
   subroutine gridexac
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: ierr, iuni, ix, iy, iz, nx, ny, nz
    real(KREAL) :: b2a, den, dendrvx, dendrvy, dendrvz, denlplc, rx, ry, rz, x, y, z
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: array(:)
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: array(:)
#endif
    rx = (xsup - xinf) / dltx + umed
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + umed
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + umed
    nz = int(rz) + 1

    allocate(array(nx), lnegia(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating array and arraysp in gridexact. Stop')

    write(6,"('3D GRID (inf,sup,dlt,npnts)')")
    write(6,"('x: ',3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
    write(6,"('y: ',3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
    write(6,"('z: ',3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz
    if (langstrom) then        ! Converts grid distances to angstrom
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
    iuni = 1
    call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"_exact-d.plt")
    idimzlm = (mxl+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zlma in gridexact. Stop')
    do iz = 1, nz
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix = 1, nx
                x = xinf + (ix - 1) * dltx
                if (lsto) then
                    call densexactaSTO(x, y, z, den)
                else
                    call densexactaGTO(x, y, z, den)
                endif
                array(ix) = den
            enddo
            call linea( 1, nx , array )
        enddo
    enddo
!     Closes the grid file
    deallocate(array, lnegia)
    close(1)
    write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    deallocate(zlma)
    return
    end
!    
!   ***************************************************************
!
   subroutine gridrep
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, ierr, isel, iuni, ix, iy, iz, nx, ny, nz
    real(KREAL) :: b2a, den, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz, dV, rx, ry, rz, x, y, z
    character*4 :: strbux
    character*256 :: straux
    character*7 :: strcux

    real(KREAL) :: qdefpos(0:mxsel), qdefneg(0:mxsel)
    real(KREAL), allocatable :: array(:), arrayacc(:), arrayaccdx(:), arrayaccdy(:), arrayaccdz(:)
    real(KREAL), allocatable :: arrayat(:,:), arrayatdx(:,:), arrayatdy(:,:), arrayatdz(:,:)
    real(KREAL), allocatable :: arraydx(:), arraydy(:), arraydz(:), arraydxx(:), arraydxy(:), arraydxz(:)
    real(KREAL), allocatable :: arraydyy(:), arraydyz(:), arraydzz(:), arraylpl(:)
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: arraysp(:)
#endif

    rx = (xsup - xinf) / dltx + umed
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + umed
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + umed
    nz = int(rz) + 1

    if (.not. (lmolec .or. latomsel .or. latomics)) then
            write(6,"('WARNING: No tabulation option for grid has been chosen.')")
            return
    endif

    write(6,"('3D GRID (inf,sup,dlt,npnts)')")
    write(6,"('x: ',3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
    write(6,"('y: ',3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
    write(6,"('z: ',3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz
    if (langstrom) then        ! Converts grid distances to angstrom
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

    allocate(arraysp(nx), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating arraysp in gridrep. Stop')
    if (ldeform) then
        strcux="-deform"
    else
        strcux=""
    endif
    if (lmolec) then
        allocate(array(nx), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating array in gridrep. Stop')
        straux = trim(filename)//trim(strcux)
        iuni = 21
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d.plt")
        if (lgradient) then
            allocate(arraydx(nx), arraydy(nx), arraydz(nx), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraydx, arraydy and arraydz in gridrep. Stop')
            iuni = 23
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dx.pltd")
            iuni = 24
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dy.pltd")
            iuni = 25
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dz.pltd")
        endif
        if (lderiv2) then
            allocate(arraydxx(nx), arraydxy(nx), arraydxz(nx), arraydyy(nx), arraydyz(nx), arraydzz(nx), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraydx, arraydy and arraydz in gridrep. Stop')
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
        endif
        if (laplacian) then
            allocate(arraylpl(nx), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraylpl in gridrep. Stop')
            iuni = 27
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-lplc.plt")
        endif
    endif
    if (ldensacc) then
        allocate(arrayacc(nx), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating arrayacc in gridrep. Stop')
        write(6,"('Indices of atoms selected for density accumulation = ', 20(i6))") (iatomsel(i), i = 1, nsel)
        straux = trim(filename)//trim(strcux)//"-fgr"
        iuni = 22
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d.plt")
        if (lgradient) then
            allocate(arrayaccdx(nx), arrayaccdy(nx), arrayaccdz(nx), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayaccdx, arrayaccdy and arrayaccdz in gridrep. Stop')
            iuni = 34
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dx.pltd")
            iuni = 35
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dy.pltd")
            iuni = 36
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(straux)//"-d-dz.pltd")
        endif
    endif
    if (latomics) then
        write(6,"('Indices of atoms selected for tabulation = ', 20(i6))") (iatomsel(i), i = 1, nsel)
        allocate(arrayat(nx,nsel), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating arrayat in gridrep. Stop')
        if (lgradient) then
            allocate(arrayatdx(nx,nsel), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayatdx in gridrep. Stop')
            allocate(arrayatdy(nx,nsel), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayatdy in gridrep. Stop')
            allocate(arrayatdz(nx,nsel), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayatdz in gridrep. Stop')
        endif
        iuni = 40
        do i = 1, nsel
            iuni = iuni + 1
            write(strbux,'(i4.1)') iatomsel(i)
            straux = trim(filename)//trim(strcux)//"-"//trim(atmnms(nzn(iatomsel(i))))//trim(adjustl(strbux))
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-"//trim(straux)//"-d.plt")
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
    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error in gridrep when allocating zlma. Stop')
    if (lgradient) then
        allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error in densrepr when allocating zlmadx, zlmady, zlmadz. Stop')
    endif
    if (lderiv2) then
        allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                zlmadzz(idimzlm), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error in densrepr when allocating zlmadxx, zlmadxy,... zlmadzz. Stop')
    endif
!    Initializes deformation charges
    qdefpos = cero
    qdefneg = cero
    do iz = 1, nz
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix = 1, nx
                x = xinf + (ix - 1) * dltx
                if (lmolec) then
                    array(ix) = cero
                    if (lgradient) then
                        arraydx(ix) = cero
                        arraydy(ix) = cero
                        arraydz(ix) = cero
                    endif
                    if (lderiv2) then
                        arraydxx(ix) = cero
                        arraydxy(ix) = cero
                        arraydxz(ix) = cero
                        arraydyy(ix) = cero
                        arraydyz(ix) = cero
                        arraydzz(ix) = cero
                    endif
                    if (laplacian) arraylpl(ix) = cero
                endif
                if (ldensacc) then
                    arrayacc(ix) = cero
                    if (lgradient) then
                        arrayaccdx(ix) = cero
                        arrayaccdy(ix) = cero
                        arrayaccdz(ix) = cero
                    endif
                endif
                if (latomics) then
                    arrayat(ix,1:nsel) = cero
                    if (lgradient) then
                        arrayatdx(ix,1:nsel) = cero
                        arrayatdy(ix,1:nsel) = cero
                        arrayatdz(ix,1:nsel) = cero
                    endif
                endif

                isel = 0
                do ia = 1, ncen
                    if (nintervaj .le. 0) cycle
                    call densrepr(ia, x, y, z, den, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz)
                    if (lmolec) then
                        array(ix) = array(ix) + den
                        if (lgradient) then
                            arraydx(ix) = arraydx(ix) + dendrvx
                            arraydy(ix) = arraydy(ix) + dendrvy
                            arraydz(ix) = arraydz(ix) + dendrvz
                        endif
                        if (lderiv2) then
                            arraydxx(ix) = arraydxx(ix) + dxx
                            arraydxy(ix) = arraydxy(ix) + dxy
                            arraydxz(ix) = arraydxz(ix) + dxz
                            arraydyy(ix) = arraydyy(ix) + dyy
                            arraydyz(ix) = arraydyz(ix) + dyz
                            arraydzz(ix) = arraydzz(ix) + dzz
                        endif
                        if (laplacian) arraylpl(ix) = arraylpl(ix) + denlplc
                    endif
                    if (ldensacc .and. lselat(ia)) then
                        arrayacc(ix) = arrayacc(ix) + den
                        if (lgradient) then
                            arrayaccdx(ix) = arrayaccdx(ix) + dendrvx
                            arrayaccdy(ix) = arrayaccdy(ix) + dendrvy
                            arrayaccdz(ix) = arrayaccdz(ix) + dendrvz
                        endif
                    endif
                    if (latomics .and. lselat(ia)) then
                        isel = isel+1
                        arrayat(ix,isel) = den
                        if (lgradient) then
                            arrayatdx(ix,isel) = dendrvx
                            arrayatdy(ix,isel) = dendrvy
                            arrayatdz(ix,isel) = dendrvz
                        endif
                    endif
                    if(    (.not. lboundsx .or. (x .gt. xboundinf .and. x .lt. xboundsup)) .and. &
                            (.not. lboundsy .or. (y .gt. yboundinf .and. y .lt. yboundsup)) .and. &
                            (.not. lboundsz .or. (z .gt. zboundinf .and. z .lt. yboundsup)) ) then
                        if (ldensacc) then
                            if(lselat(ia)) then
                                if (den .gt. cero) then
                                    qdefpos(0) = qdefpos(0) + den
                                else
                                    qdefneg(0) = qdefneg(0) + den
                                endif
                            endif
                        else
                            if (den .gt. cero) then
                                qdefpos(0) = qdefpos(0) + den
                            else
                                qdefneg(0) = qdefneg(0) + den
                            endif
                        endif
                    endif
                    if (latomics .and. lselat(ia)) then
                        if (den .gt. cero) then
                            qdefpos(isel) = qdefpos(isel) + den
                        else
                            qdefneg(isel) = qdefneg(isel) + den
                        endif
                    endif
                enddo
            enddo
            if (lmolec) then
                arraysp = array
                call linea(21, nx , arraysp )
                if (lgradient) then
                    arraysp = arraydx
                    call linea(23, nx , arraysp )
                    arraysp = arraydy
                    call linea(24, nx , arraysp )
                    arraysp = arraydz
                    call linea(25, nx , arraysp )
                endif
                if (lderiv2) then
                    arraysp = arraydxx
                    call linea(28, nx , arraysp )
                    arraysp = arraydxy
                    call linea(29, nx , arraysp )
                    arraysp = arraydxz
                    call linea(30, nx , arraysp )
                    arraysp = arraydyy
                    call linea(31, nx , arraysp )
                    arraysp = arraydyz
                    call linea(32, nx , arraysp )
                    arraysp = arraydzz
                    call linea(33, nx , arraysp )
                endif
                if (laplacian ) then
                    arraysp = arraylpl
                    call linea(27, nx , arraysp )
                endif
            endif
            if (ldensacc ) then
                arraysp = arrayacc
                call linea(22, nx , arraysp )
                if (lgradient) then
                    arraysp = arraydx
                    call linea(34, nx , arraysp )
                    arraysp = arraydy
                    call linea(35, nx , arraysp )
                    arraysp = arraydz
                    call linea(36, nx , arraysp )
                endif
            endif
            if (latomics ) then
                iuni = 40
                do i = 1, nsel
                    iuni = iuni+1
                    arraysp = arrayat(:,i)
                    call linea( iuni, nx , arraysp )
                    if (lgradient) then
                        iuni = iuni+1
                        arraysp = arrayatdx(:,i)
                        call linea( iuni, nx , arraysp )
                        iuni = iuni+1
                        arraysp = arrayatdy(:,i)
                        call linea( iuni, nx , arraysp )
                        iuni = iuni+1
                        arraysp = arrayatdz(:,i)
                        call linea( iuni, nx , arraysp )
                    endif
                enddo
            endif
        enddo
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(arraysp)
    if (lmolec) then
        deallocate(array)
        close(21)
        if (lgradient) then
            deallocate(arraydx, arraydy, arraydz)
            close(23)
            close(24)
            close(25)
        endif
        if (lderiv2) then
            deallocate(arraydxx, arraydxy, arraydxz, arraydyy, arraydyz, arraydzz)
            close(28)
            close(29)
            close(30)
            close(31)
            close(32)
            close(33)
        endif
        if (laplacian ) then
            deallocate(arraylpl)
            close(27)
        endif
    endif
    if (ldensacc ) then
        deallocate(arrayacc)
        close(22)
        if (lgradient) then
            deallocate(arrayaccdx, arrayaccdy, arrayaccdz)
            close(34)
            close(35)
            close(36)
        endif
    endif

    if (latomics ) then
        deallocate(arrayat, arrayatdx, arrayatdy, arrayatdz)
        iuni = 40
        do i = 1, nsel
            iuni = i + 1
            close(iuni)
            if (lgradient) then
                iuni = i + 1
                close(iuni)
                iuni = i + 1
                close(iuni)
                iuni = i + 1
                close(iuni)
            endif
        enddo
    endif
    write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    if (lminrep .eq. 1) then
        write(6,"('Deformation charges accumulated in the region:')")
        if (lboundsx) then
            write(6,"(7x, f10.5, ' < x < ', f10.5)") xboundinf, xboundsup
        else
            write(6,"(7x, f10.5, ' < x < ', f10.5)") xinf, xsup
        endif
        if (lboundsy) then
            write(6,"(7x, f10.5, ' < y < ', f10.5)") yboundinf, yboundsup
        else
            write(6,"(7x, f10.5, ' < y < ', f10.5)") yinf, ysup
        endif
        if (lboundsz) then
            write(6,"(7x, f10.5, ' < z < ', f10.5)") zboundinf, zboundsup
        else
            write(6,"(7x, f10.5, ' < z < ', f10.5)") zinf, zsup
        endif
        dV = dltx*dlty*dltz
        write(6,"(/2x,'Accumulated negative charge = ', f12.5, 7x,'Accumulated positive charge = ', f12.5)") &
            -qdefpos(0) * dV, -qdefneg(0) * dV
        if (latomics) then
            write(6,"(/9x,'Atom      Negative charge     Positive charge'/)")
            do i = 1, nsel
                if (lselat(ia)) write(6,"(7x, a2, i5, 8x, f12.5, 12x, f12.5)") atmnam(iatomsel(i)), iatomsel(i) &
                        , -qdefpos(i)*dV , -qdefneg(i)*dV
            enddo
        endif
    endif
    deallocate(zlma)
    if (allocated(zlmadx)) deallocate(zlmadx, zlmady, zlmadz)
    if (allocated(zlmadxx)) deallocate(zlmadxx, zlmadxy, zlmadxz, zlmadyy, zlmadyz, zlmadzz)
    return
    end
!    
!   ***************************************************************
!
   subroutine gridexac2d
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    use strings
    use evaluate
    implicit none
    integer(KINT) :: iuni, iu, iv, nu, nv
    real(KREAL) :: b2a, den, dendrvx, dendrvy, dendrvz, denlplc, ru, rv, u, v, x, y, z
#ifdef DBLPRCGRID
    real(KREAL) :: u1, u2, v1, v2
    real(KREAL), allocatable :: array(:)
#else
    real(KREAL4) :: u1, u2, v1, v2
    real(KREAL4), allocatable :: array(:)
#endif

    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1

    allocate(array(nu), lnegia(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating array, arraysp and lnegia in gridexact. Stop')

    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1
    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv
    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        u1 = uinf * b2a
        v1 = vinf * b2a
        u2 = (uinf+(nu-1)*dltu) * b2a
        v2 = (vinf+(nv-1)*dltv) * b2a
    else
        u1 = uinf
        v1 = vinf
        u2 = (uinf+(nu-1)*dltu)
        v2 = (vinf+(nv-1)*dltv)
    endif
    iuni = 1
    call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//"_exact-d.cnt")
    idimzlm = (mxl+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating zlma in gridexac2d. Stop')
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)
            call evalexpr(y_func_uv,y)
            call evalexpr(z_func_uv,z)
            if (lsto) then
                call densexactaSTO(x, y, z, den)
            else
                call densexactaGTO(x, y, z, den)
            endif
            array(iu) = den
        enddo
        call linea( 1, nu , array )
    enddo
!     Closes the grid file
    close(1)
    write(6,"('Total number of tabulated points = ', i12)") nu*nv
    deallocate(zlma, lnegia)
    return
    end
!    
!   ***************************************************************
!
   subroutine gridrep2d
    USE DAM320_D
    USE DAMDEN320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    use strings
    use evaluate
    implicit none
    integer(KINT) :: i, ia, isel, iuni, iu, iv, nu, nv
    real(KREAL) :: b2a, den, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz, ru, rv, u, v, x, y, z
    character*4 :: strbux
    character*256 :: straux
    character*7 :: strcux
    real(KREAL), allocatable :: array(:), arrayacc(:), arrayaccdx(:), arrayaccdy(:), arrayaccdz(:)
    real(KREAL), allocatable :: arrayat(:,:), arrayatdx(:,:), arrayatdy(:,:), arrayatdz(:,:)
    real(KREAL), allocatable :: arraydx(:), arraydy(:), arraydz(:), arraydxx(:), arraydxy(:), arraydxz(:)
    real(KREAL), allocatable :: arraydyy(:), arraydyz(:), arraydzz(:), arraylpl(:)
#ifdef DBLPRCGRID
    real(KREAL) :: u1, u2, v1, v2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: u1, u2, v1, v2
    real(KREAL4), allocatable :: arraysp(:)
#endif

    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1
    if (.not. (lmolec .or. latomsel .or. latomics)) then
            write(6,"('WARNING: No tabulation option for grid has been chosen.')")
            return
    endif

    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1

    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv

    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        u1 = uinf * b2a
        v1 = vinf * b2a
        u2 = (uinf+(nu-1)*dltu) * b2a
        v2 = (vinf+(nv-1)*dltv) * b2a
    else
        u1 = uinf
        v1 = vinf
        u2 = (uinf+(nu-1)*dltu)
        v2 = (vinf+(nv-1)*dltv)
    endif

    allocate(arraysp(nu), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating arraysp in gridrep. Stop')
    if (ldeform) then
        strcux="-deform"
    else
        strcux=""
    endif
    if (lmolec) then
        allocate(array(nu), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating array in gridrep. Stop')
        straux = trim(filename)//trim(strcux)//trim(planesuffix)
        iuni = 21
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d.cnt")
        if (lgradient) then
            allocate(arraydx(nu), arraydy(nu), arraydz(nu), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraydx, arraydy and arraydz in gridrep. Stop')
            iuni = 23
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dx.cnt")
            iuni = 24
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dy.cnt")
            iuni = 25
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dz.cnt")
        endif
        if (lderiv2) then
            allocate(arraydxx(nu), arraydxy(nu), arraydxz(nu), arraydyy(nu), arraydyz(nu), arraydzz(nu), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraydx, arraydy and arraydz in gridrep. Stop')
            iuni = 28
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dxx.cnt")
            iuni = 29
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dxy.cnt")
            iuni = 30
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dxz.cnt")
            iuni = 31
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dyy.cnt")
            iuni = 32
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dyz.cnt")
            iuni = 33
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dzz.cnt")
        endif
        if (laplacian) then
            allocate(arraylpl(nu), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraylpl in gridrep. Stop')
            iuni = 27
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-lplc.cnt")
        endif
    endif
    if (ldensacc) then
        allocate(arrayacc(nu), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating arrayacc in gridrep. Stop')
        write(6,"('Indices of atoms selected for density accumulation = ', 20(i6))") (iatomsel(i), i = 1, nsel)
        iuni = 22
        straux = trim(filename)//trim(strcux)//"-fgr"//trim(planesuffix)
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d.cnt")
        if (lgradient) then
            allocate(arrayaccdx(nu), arrayaccdy(nu), arrayaccdz(nu), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayaccdx, arrayaccdy and arrayaccdz in gridrep. Stop')
            iuni = 34
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dx.cnt")
            iuni = 35
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dy.cnt")
            iuni = 36
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dz.cnt")
        endif
    endif
    if (latomics) then
        write(6,"('Indices of atoms selected for tabulation = ', 20(i6))") (iatomsel(i), i = 1, nsel)
        allocate(arrayat(nu,nsel), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating arrayat in gridrep. Stop')
        if (lgradient) then
            allocate(arrayatdx(nu,nsel), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayatdx in gridrep. Stop')
            allocate(arrayatdy(nu,nsel), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayatdy in gridrep. Stop')
            allocate(arrayatdz(nu,nsel), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arrayatdz in gridrep. Stop')
        endif
        iuni = 40
        do i = 1, nsel
            iuni = iuni + 1
            write(strbux,'(i4.1)') iatomsel(i)
            straux = trim(filename)//trim(strcux)//"-"//trim(atmnms(nzn(iatomsel(i))))//trim(adjustl(strbux))//trim(planesuffix)
            call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d.cnt")
            if (lgradient) then
                iuni = iuni + 1
                call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dx.cnt")
                iuni = iuni + 1
                call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dy.cnt")
                iuni = iuni + 1
                call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(straux)//"-d-dz.cnt")
            endif
        enddo
    endif
    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error in gridrep when allocating zlma. Stop')
    if (lgradient) then
        allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error in densrepr when allocating zlmadx, zlmady, zlmadz. Stop')
    endif
    if (lderiv2) then
        allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                zlmadzz(idimzlm), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error in densrepr when allocating zlmadxx, zlmadxy,... zlmadzz. Stop')
    endif
!    Grid tabulation
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)
            call evalexpr(y_func_uv,y)
            call evalexpr(z_func_uv,z)
            if (lmolec) then
                array(iu) = cero
                if (lgradient) then
                    arraydx(iu) = cero
                    arraydy(iu) = cero
                    arraydz(iu) = cero
                endif
                if (lderiv2) then
                    arraydxx(iu) = cero
                    arraydxy(iu) = cero
                    arraydxz(iu) = cero
                    arraydyy(iu) = cero
                    arraydyz(iu) = cero
                    arraydzz(iu) = cero
                endif
                if (laplacian) arraylpl(iu) = cero
            endif
            if (ldensacc) then
                arrayacc(iu) = cero
                if (lgradient) then
                    arrayaccdx(iu) = cero
                    arrayaccdy(iu) = cero
                    arrayaccdz(iu) = cero
                endif
            endif
            if (latomics) then
                arrayat(iu,1:nsel) = cero
                if (lgradient) then
                    arrayatdx(iu,1:nsel) = cero
                    arrayatdy(iu,1:nsel) = cero
                    arrayatdz(iu,1:nsel) = cero
                endif
            endif

            isel = 0
            do ia = 1, ncen
                if (nintervaj .le. 0) cycle
                call densrepr(ia, x, y, z, den, dendrvx, dendrvy, dendrvz, denlplc, dxx, dxy, dxz, dyy, dyz, dzz)
                if (lmolec) then
                    array(iu) = array(iu) + den
                    if (lgradient) then
                        arraydx(iu) = arraydx(iu) + dendrvx
                        arraydy(iu) = arraydy(iu) + dendrvy
                        arraydz(iu) = arraydz(iu) + dendrvz
                    endif
                    if (lderiv2) then
                        arraydxx(iu) = arraydxx(iu) + dxx
                        arraydxy(iu) = arraydxy(iu) + dxy
                        arraydxz(iu) = arraydxz(iu) + dxz
                        arraydyy(iu) = arraydyy(iu) + dyy
                        arraydyz(iu) = arraydyz(iu) + dyz
                        arraydzz(iu) = arraydzz(iu) + dzz
                    endif
                    if (laplacian) arraylpl(iu) = arraylpl(iu) + denlplc
                endif
                if (ldensacc .and. lselat(ia)) then
                    arrayacc(iu) = arrayacc(iu) + den
                    if (lgradient) then
                        arrayaccdx(iu) = arrayaccdx(iu) + dendrvx
                        arrayaccdy(iu) = arrayaccdy(iu) + dendrvy
                        arrayaccdz(iu) = arrayaccdz(iu) + dendrvz
                    endif
                endif
                if (latomics .and. lselat(ia)) then
                    isel = isel+1
                    arrayat(iu,isel) = den
                    if (lgradient) then
                        arrayatdx(iu,isel) = dendrvx
                        arrayatdy(iu,isel) = dendrvy
                        arrayatdz(iu,isel) = dendrvz
                    endif
                endif
            enddo
        enddo
        if (lmolec) then
            arraysp = array
            call linea(21, nu , arraysp )
            if (lgradient) then
                arraysp = arraydx
                call linea(23, nu , arraysp )
                arraysp = arraydy
                call linea(24, nu , arraysp )
                arraysp = arraydz
                call linea(25, nu , arraysp )
            endif
            if (lderiv2) then
                arraysp = arraydxx
                call linea(28, nu , arraysp )
                arraysp = arraydxy
                call linea(29, nu , arraysp )
                arraysp = arraydxz
                call linea(30, nu , arraysp )
                arraysp = arraydyy
                call linea(31, nu , arraysp )
                arraysp = arraydyz
                call linea(32, nu , arraysp )
                arraysp = arraydzz
                call linea(33, nu , arraysp )
            endif
            if (laplacian ) then
                arraysp = arraylpl
                call linea(27, nu , arraysp )
            endif
        endif
        if (ldensacc ) then
            arraysp = arrayacc
            call linea(22, nu , arraysp )
            if (lgradient) then
                arraysp = arrayaccdx
                call linea(34, nu , arraysp )
                arraysp = arrayaccdy
                call linea(35, nu , arraysp )
                arraysp = arrayaccdz
                call linea(36, nu , arraysp )
            endif
        endif
        if (latomics ) then
            iuni = 40
            do i = 1, nsel
                iuni = iuni+1
                arraysp = arrayat(:,i)
                call linea( iuni, nu , arraysp )
                if (lgradient) then
                    iuni = iuni+1
                    arraysp = arrayatdx(:,i)
                    call linea( iuni, nu , arraysp )
                    iuni = iuni+1
                    arraysp = arrayatdy(:,i)
                    call linea( iuni, nu , arraysp )
                    iuni = iuni+1
                    arraysp = arrayatdz(:,i)
                    call linea( iuni, nu , arraysp )
                endif
            enddo
        endif
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(arraysp)
    if (lmolec) then
        deallocate(array)
        close(21)
        if (lgradient) then
            deallocate(arraydx, arraydy, arraydz)
            close(23)
            close(24)
            close(25)
        endif
        if (lderiv2) then
            deallocate(arraydxx, arraydxy, arraydxz, arraydyy, arraydyz, arraydzz)
            close(28)
            close(29)
            close(30)
            close(31)
            close(32)
            close(33)
        endif
        if (laplacian ) then
            deallocate(arraylpl)
            close(27)
        endif
    endif
    if (ldensacc ) then
        deallocate(arrayacc)
        close(22)
        if (lgradient) then
            deallocate(arrayaccdx, arrayaccdy, arrayaccdz)
            close(34)
            close(35)
            close(36)
        endif
    endif

    if (latomics ) then
        deallocate(arrayat)
        iuni = 40
        do i = 1, nsel
            iuni = i + 1
            close(iuni)
            if (lgradient) then
                deallocate(arrayatdx, arrayatdy, arrayatdz)
                iuni = i + 1
                close(iuni)
                iuni = i + 1
                close(iuni)
                iuni = i + 1
                close(iuni)
            endif
        enddo
    endif
    write(6,"('Total number of tabulated points = ', i12)") nu*nv
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
    implicit none
    integer(KINT) :: i, i1, ia, ib, ierr, j, knt, kntb, l, lm, m
    real(KREAL) :: aux, denex, frad, ra, ra2, x, xxa, y, yya, z, zza
    real(KREAL) :: rapow(0:mxn)
    logical :: ldenex
    allocate(ang((mxl+1)*(mxl+2)/2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ang in densexactaSTO. Stop')
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
    mxind = (mxl+1)*(mxl+2)/2
    allocate(ind(0:mxind), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ind in densexactaSTO. Stop')
    ind(0) = 0
    do i = 1, mxind
            ind(i) = ind(i-1) + i         !  i*(i+1)/2
    enddo
    knt = 0
    denex = cero
    ldenex = .false.
    do ia = 1, ncen
        lnegia(ia) = .true.
        if (ngini(ia) .le. 0) cycle    ! If center without basis set, skips to next atom
        xxa = x - rcen(1,ia)
        yya = y - rcen(2,ia)
        zza = z - rcen(3,ia)
        ra2 = xxa*xxa+yya*yya+zza*zza
        if (ra2 .gt. rlargo(ia)*rlargo(ia)) cycle    ! Contributions of atom ia are negligible. Skips to next atom
        ldenex = .true.
        lnegia(ia) = .false.
        ra = sqrt(ra2)
        zlma(1) = uno        ! Regular spherical harmonics of r-R(ia)
        zlma(2) = yya
        zlma(3) = zza
        zlma(4) = xxa
        do l = 1, lmaxc(ia)
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xxa * zlma(l*(l+2)+1) - yya * zlma(l*l+1))        ! zlm(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (yya * zlma(l*(l+2)+1) + xxa * zlma(l*l+1))        ! zlm(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * zza* zlma(l*(l+2)+1)                ! zlm(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * zza * zlma(l*l+1)                    ! zlm(l+1,-l,ia)
            do m = 0, l-1
                    zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))    ! zlm(l+1,m,ia)
                    zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
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
        if (lnegia(ia)) cycle        ! Contributions of atom ia are negligible. Skips to next atom
        kntb = 0
        do ib = 1, ncen
            if (lnegia(ib)) cycle    ! Contributions of atom ib are negligible. Skips to next atom
            do j = nf(ngini(ib)), nf(ngfin(ib))+2*ll(ngfin(ib))
                faux(kntb+j-nf(ngini(ib))+1) = faux(kntb+j-nf(ngini(ib))+1) + &
                        dot_product(f(knt+1:knt+nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1), &
                                        dmat(nf(ngini(ia)):nf(ngfin(ia))+2*ll(ngfin(ia)), j))
            enddo
            kntb = kntb -nf(ngini(ib))+ nf(ngfin(ib))+2*ll(ngfin(ib)) + 1
        enddo
        knt = knt + nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1
    enddo
    denex = max(dot_product(faux(1:knt),f(1:knt)),cero)    ! Total density is positive defined
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
    implicit none
    integer(KINT) :: i, i1, i1p, ia, ib, ierr, j, knt, kntb, l, lm, m
    real(KREAL) :: aux, denex, frad, ra, ra2, x, xxa, y, yya, z, zza
    logical :: ldenex
    allocate(ang((mxl+1)*(mxl+2)/2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ang in densexactaGTO. Stop')
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
    mxind = (mxl+1)*(mxl+2)/2
    allocate(ind(0:mxind), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ind in densexactaGTO. Stop')
    ind(0) = 0
    do i = 1, mxind
        ind(i) = ind(i-1) + i         !  i*(i+1)/2
    enddo
    knt = 0
    denex = cero
    ldenex = .false.
    do ia = 1, ncen
        lnegia(ia) = .true.
        if (ngini(ia) .le. 0) cycle    ! If center without basis set, skips to next atom
        xxa = x - rcen(1,ia)
        yya = y - rcen(2,ia)
        zza = z - rcen(3,ia)
        ra2 = xxa*xxa+yya*yya+zza*zza
        if (ra2 .gt. rlargo(ia)*rlargo(ia)) cycle    ! Contributions of atom ia are negligible. Skips to next atom
        ldenex = .true.
        lnegia(ia) = .false.
        ra = sqrt(ra2)
!        Computes the regular spherical harmonics of r-R(ia)
        zlma(1) = uno
        zlma(2) = yya
        zlma(3) = zza
        zlma(4) = xxa
        do l = 1, lmaxc(ia)
            zlma((l+1)*(l+3)+1) = dosl1(l) * (xxa * zlma(l*(l+2)+1) - yya * zlma(l*l+1))        ! zlm(l+1,l+1,ia)
            zlma((l+1)*(l+1)+1) = dosl1(l) * (yya * zlma(l*(l+2)+1) + xxa * zlma(l*l+1))        ! zlm(l+1,-(l+1),ia)
            zlma((l+2)*(l+2)-1) = dosl1(l) * zza* zlma(l*(l+2)+1)                ! zlm(l+1,l,ia)
            zlma(l*(l+2)+3) = dosl1(l) * zza * zlma(l*l+1)                    ! zlm(l+1,-l,ia)
            do m = 0, l-1
                zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))    ! zlm(l+1,m,ia)
                zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*zza*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
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
        if (lnegia(ia)) cycle        ! Contributions of atom ia are negligible. Skips to next atom
        kntb = 0
        do ib = 1, ncen
            if (lnegia(ib)) cycle    ! Contributions of atom ib are negligible. Skips to next atom
            do j = nf(ngini(ib)), nf(ngfin(ib))+2*ll(ngfin(ib))
                faux(kntb+j-nf(ngini(ib))+1) = faux(kntb+j-nf(ngini(ib))+1) + &
                        dot_product(f(knt+1:knt+nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1), &
                                        dmat(nf(ngini(ia)):nf(ngfin(ia))+2*ll(ngfin(ia)), j))
            enddo
            kntb = kntb -nf(ngini(ib))+ nf(ngfin(ib))+2*ll(ngfin(ib)) + 1
        enddo
        knt = knt + nf(ngfin(ia))-nf(ngini(ia))+2*ll(ngfin(ia))+1
    enddo
    denex = max(dot_product(faux(1:knt),f(1:knt)),cero)    ! Total density is positive defined
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
!    Contribution of atomic fragment ia to density  in point (x,y,z)
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

        denrep = cero
        denlplc = cero
        kntlm = lminrep*lminrep
        do l = lminrep, lcorto(interv,ia)    !     Computes density terms lminrep <= l <= lcorto(interv,ia)
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
        tcheb(0) = uno    ! Chebyshev T  polynomials
        tcheb(1) = t
        do j = 2, mxlenpol-1
            tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
        enddo
        denrep = cero
        kntlm = lminrep*lminrep
        do l = lminrep, lcorto(interv,ia)    !     Computes density terms lminrep <= l <= lcorto(interv,ia)
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
  subroutine subplanesuffix(pc, ps)
    USE DAM320_D
    implicit none
    integer(KINT) :: pc
    character(4) :: ps
    character(4), parameter :: suffixes(7) = (/"_XY0", "_X0Z", "_0YZ", "_AB0", "_A0C", "_0BC", "_ABC" /)
    if (pc .ge. 1 .and. pc .le. 7) then
        ps = suffixes(pc)
    else
        ps = ""
    endif
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
