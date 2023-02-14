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
! Program for tabulation of the electrostatic potential from the representation of the molecular density performed with
! DAM320
!
! Parallel version with MPI
!
! Version of April 2019
!

! #define DBLPRCGRID	! Uncomment this line  if double precision grid is wanted
  program DAMPOT320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMPOT320_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, ierr, ipoint, j, nbasis
    real(KREAL) :: aux, drvx, drvxtot, drvy, drvytot, drvz, drvztot
    real(KREAL) :: dxx, dxxtot, dxy, dxytot, dxz, dxztot, dyy,dyytot, dyz, dyztot, dzz, dzztot
    real(KREAL) :: ve, vetot, vn, vntot, vtot, vtotia, x, xmax, xmin, xyzmax, xyzmin, y, ymax, ymin, z, zmax, zmin
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    logical :: lnamelist(8), ltimeprocs
    integer(KINT) :: inamelist(1)
    real(KREAL) :: rnamelist(10)
    namelist / options / dltu, dltv, dltx, dlty, dltz, iswindows, filename, geomthr, langstrom, largo, lderiv2, lexact &
            , lgradient, lgrid, lgrid2d, lmaxrep, longoutput, lpoints, lvalence, numrtab &
            , planeA, planeB, planeC, planecase &
            , rtab, uinf, umbrlargo, usup, vinf, vsup &
            , x_func_uv, xinf, xsup, y_func_uv, yinf, ysup, z_func_uv, zinf, zsup
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
    geomthr = 1.d-10		! Geometry threshold: two points at a distance lower than geomthr are considered to be the coincident
    longoutput = .false.	! If .true. a more detailed output is given
    langstrom = .true.		! If .false. distances in bohr
    largo = .false.		! if .true. long-range potential
    lexact  = .false.		! if .true. "exact" potential is tabulated
    lgrid = .true.			! If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    lmaxrep = 10			! highest "l" in the expansion of the potential
    lgradient = .false.		! If .true. gradient components of the density computed and, if lgrid = .true., tabulated in files
                                            !    projectname-v-dx.pltd, projectname-v-dy.pltd, projectname-v-dz.pltd
    lderiv2 = .false.		! If .true. second derivatives of the density computed and, if lgrid = .true., tabulated in files
                                            !	projectname-v-dxx.pltd, projectname-v-dxy.pltd, projectname-v-dxz.pltd, projectname-v-dyz.pltd, etc
    lpoints = .false.		! If .true. computes in selected points and prints out the results in the standard output.
                                            ! Points must be given in cartesian coordinates.
                                            ! If lgrid .eq. .true., these coordinates must be placed after the grid data
    lvalence = .false.      ! If .true. only valence electrons are considered
    filename = ""			! root file name for .plt files
    umbrlargo = 1.d-9		! Threshold for determining the short-range radius
    numrtab = 0			! Number of tabulation points supplied in namelist
    rtab = cero			! Tabulation points supplied in namelist
    iswindows = .false.		! .true. if running on a MS-windows system
    xinf = cero
    xsup = cero
    dltx = uno
    yinf = cero
    ysup = cero
    dlty = uno
    zinf = cero
    zsup = cero
    dltz = uno
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
        if (lderiv2) lgradient = .true.
        read(5,*) projectname
        inamelist = (/ lmaxrep /)
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

        if (lvalence) then
            write(6,"(//'Using valence electrons only',//)")
        endif

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
        lnamelist = (/ langstrom, lgrid, largo, lgradient, lderiv2, lexact, ltimeprocs, lvalence /)
    endif
    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(lnamelist,8,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,10,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    if (myrank .ne. 0) then
        lmaxrep = inamelist(1)
        langstrom = lnamelist(1); lgrid = lnamelist(2); largo = lnamelist(3); lgradient = lnamelist(4);
        lderiv2 = lnamelist(5); lexact = lnamelist(6); ltimeprocs = lnamelist(7); lvalence = lnamelist(8)
        xinf = rnamelist(1); xsup = rnamelist(2); dltx = rnamelist(3)
        yinf = rnamelist(4); ysup = rnamelist(5); dlty = rnamelist(6)
        zinf = rnamelist(7); zsup = rnamelist(8); dltz = rnamelist(9)
        umbrlargo = rnamelist(10)
    endif
    if (myrank .eq. 0) then
        if (.not. lexact) then
            write(6,"('Potential from expansion of the density: lmaxrep = ', i3)") lmaxrep
        else
            write(6,"('Potential computed from density matrix and basis set')")
        endif
    endif
    call consta		!	Computes auxiliary constants

    call leedamqtpot		!	Reads file _2016.damqt and _2016.dmqtv (generated by DAM2016)
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
            call error(1,'Stop')
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (myrank .eq. 0) then
!       restores files back to their original gzipped status
        if (ldengz) then
            call system ("gzip "//trim(projectname)//".den")
        endif
        if (lgbsgz) then
            inquire(file=trim(projectname)//".ggbs", exist=lgbsgz, iostat=ierr)
            if (ierr .eq. 0 .and. lgbsgz) then
                call system ("gzip "//trim(projectname)//".ggbs")
            else
                    call system ("gzip "//trim(projectname)//".sgbs")
            endif
        endif
    endif

#ifdef DBLPRCGRID
    if (myrank .eq. 0) then
        write(6,"('Grid generated in double precision')")
        write(6,"('WARNING! This grid will not be compatible with gOpenMol')")
    endif
#endif

    allocate (ra2l1((lmaxrep+1)**2), ra2l1inv((lmaxrep+1)**2), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ra2l1 and ra2l1inv in processor ',i3)") myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
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
                xinf = umed * (re(3) * xmin - xmax)
                xsup = umed * (re(3) * xmax - xmin)
                dltx = (xsup-xinf)/49
            endif
            if ((ysup-yinf) .eq. cero) then
                yinf = umed * (re(3) * ymin - ymax)
                ysup = umed * (re(3) * ymax - ymin)
                dlty = (ysup-yinf)/49
            endif
            if ((zsup-zinf) .eq. cero) then
                zinf = umed * (re(3) * zmin - zmax)
                zsup = umed * (re(3) * zmax - zmin)
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
!	Computes the electrostatic potential on the grid points
        if (lexact) then
            call exactmespgrid
        else
            call mespgrid
        endif
        CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        tiempo = dtime(tarray)
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
        if (largo .and. .not. lexact) then       ! Long-range potential
            if (lgradient) then
                write(6,"(//12x,'Long-range electrostatic potential from density expansion up to  l = ',i2,/10x, 95('='),  &
                //9x,'X',17x,'Y',17x,'Z',t65,'Vnuc',t88,'Vel',t111,'Vtot',t133,'der x', t156, 'der y', t180, 'der z')") &
                lmaxrep
            else
                write(6,"(//12x,'Long-range electrostatic potential from density expansion up to  l = ',i2 &
                ,/10x, 71('='), //9x,'X',17x,'Y',17x,'Z',19x,'Vnuc',19x,'Vel',19x,'Vtot')") lmaxrep
            endif
        else
            if (.not. lexact) then
                if (lgradient) then
                    write(6,"(//12x,'Electrostatic potential from density expansion up to  l = ',i2,/10x, 95('='),  &
                    //9x,'X',17x,'Y',17x,'Z',t65,'Vnuc',t88,'Vel',t111,'Vtot',t133,'der x', t156, 'der y', t180, 'der z')") &
                    lmaxrep
                else
                    write(6,"(//12x,'Electrostatic potential from density expansion up to  l = ',i2,  &
                    /10x, 71('='), //9x,'X',17x,'Y',17x,'Z',19x,'Vnuc',19x,'Vel',19x,'Vtot')") lmaxrep
                endif
            else
                write(6,"(//12x,'Electrostatic potential from density and basis set'  &
                /10x, 71('='), //9x,'X',17x,'Y',17x,'Z',19x,'Vnuc',19x,'Vel',19x,'Vtot')")
            endif
        endif
        if (lderiv2 .and. .not. lexact) then
            open(9,file=trim(projectname)//"-v.der2",form='formatted', iostat=ierr)
            if (ierr .ne. 0) then
                write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//"-v.der2", myrank
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
            allocate(zlma((mxldst+1)**2), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating zlma in processor ',i3)") myrank
                abortroot = 1
            endif
        endif
        if (abortroot .gt. 0) then
            call error(1,'Stop')
        endif
        do while (ierr .eq. 0)
            if (lexact) then
                call GTOexactmesp(x, y, z, vntot, vetot, vtot)
                write(6,"(3(1x,e17.10), 3(1x,e22.15))") x, y, z, vntot, vetot, vtot
            else
                vntot = cero
                vetot = cero
                vtot = cero
                drvxtot = cero
                drvytot = cero
                drvztot = cero
                dxxtot = cero
                dxytot = cero
                dxztot = cero
                dyytot = cero
                dyztot = cero
                dzztot = cero
                do ia = 1, ncen
                    if (largo) then
                        call longmesp(ia, x, y, z, vn, ve, vtotia, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz)
                    else
                        call mesp(ia, x, y, z, vn, ve, vtotia, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz)
                    endif
                    if (lgradient) then
                        drvxtot = drvxtot + drvx
                        drvytot = drvytot + drvy
                        drvztot = drvztot + drvz
                    endif
                    if (lderiv2) then
                        dxxtot = dxxtot + dxx
                        dxytot = dxytot + dxy
                        dxztot = dxztot + dxz
                        dyytot = dyytot + dyy
                        dyztot = dyztot + dyz
                        dzztot = dzztot + dzz
                    endif
                    vntot = vntot + vn
                    vetot = vetot + ve
                    vtot = vtot + vtotia
                enddo
                if (lgradient) then
                    if (lderiv2) then
                        write(9,"(/'Second derivatives of potential at point ',3(1x,e17.10))") x, y, z
                        write(9,"(6(1x,e22.15))") dxxtot, dxytot, dxztot, dyytot, dyztot, dzztot
                    endif
                    write(6,"(3(1x,e17.10), 3(1x,e22.15), 3(1x,e22.15))") x, y, z, vntot, vetot, vtot, &
                            drvxtot, drvytot, drvztot
                else
                    write(6,"(3(1x,e17.10), 3(1x,e22.15))") x, y, z, vntot, vetot, vtot
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
        tiempo = dtime(tarray)
        write(6,"(1x,'Timing in seconds of individual points tabulation in proc 0 (user, system, total):', &
                        5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tarray(1)+tarray(2)
        if (allocated(zlma)) deallocate(zlma)
        if (allocated(zlmadx)) deallocate(zlmadx, zlmady, zlmadz)
        if (allocated(zlmadxx)) deallocate(zlmadxx, zlmadxy, zlmadxz, zlmadyy, zlmadyz, zlmadzz)
    endif
    call MPI_FINALIZE(ierr)
    stop
    end
!**********************************************************************
!    subroutine consta
!
!	Computes and stores auxiliary constants
!		re(i) = dfloat(i)
!		ri(i) = uno / dfloat(i)
!		fact(i) = dfloat(i!)
!		facti(i) = uno / dfloat(i!)
!		facts(i) = dfloat((i+1/2)!)
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMPOT320_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ierr, k1, k12, l, l1, l1l1, l2, l2l2, lm
    integer(KINT) :: m, m1, m1a, m2, m2a, ma, mb, md, ms
    real(KREAL) :: aux, ss, sd
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
    facti(0) = uno
    facts(-1) = raizpi
    facts(0) = facts(-1) * umed
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)   		!  i!
        facts(i) = facts(i-1) * re(i+i+1) * umed	! (i+1/2)!
        facti(i) = uno / fact(i)     		!  uno / i!
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
!	***************************************************************
!
  subroutine leedamqtpot
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
    if (myrank .eq. 0 .and. longoutput) write(6,"('Size of file ', a, ' = ', i12)") trim(projectname)//"_2016.damqt", nsize
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
#if _WIN32
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#elif __INTEL_COMPILER
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#else
    open (unit=11, file=trim(projectname)//"_2016.dmqtv", form='unformatted', action = 'read', access='stream', iostat=ierr)
#endif
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(projectname)//"_2016.dmqtv", myrank
        abort = 1
        return
    endif
    if (myrank .eq. 0 .and. longoutput) write(6,"('Opens files ', a, ' and ', a)") trim(projectname)//"_2016.damqt", &
            trim(projectname)//"_2016.dmqtv"
    read(10) ncen, nbas, ncaps
    if (myrank .eq. 0) write(6,"('ncen = ', i4, ' nbas = ', i6, ' ncaps = ', i5)") ncen, nbas, ncaps

!	Allocates memory for geometry

    allocate(atmnam(ncen), nzn(ncen), rcen(3,ncen), zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam, nzn, rcen and zn in processor ',i3)") myrank
        abort = 1
        return
    endif

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
    read(10) lsto	! .true. means STO basis, .false. means GTO basis
    if (lsto) then
        if (lexact) then
            if (myrank .eq. 0) write(6,"('Exact potential not prepared for STO yet. Stop')")
            abort = 1
            return
        endif

!		Allocates memory for the basis set

        allocate(ll(ncaps), lmaxc(ncen), nf(ncaps), ngini(ncen), ngfin(ncen), nn(ncaps), rnor(ncaps), &
                xx(ncaps), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating ll, lmaxc, nf, ngini, ngfin, nn, rnor and xx in processor ',i3)") &
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
            if (myrank.eq.0 .and. longoutput)  write(6,"(t5,'center ', i4,/t12,'n', t16, 'l', t25,'exp', t35, 'ind func')") ia
            lmaxc(ia) = 0
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(10) nf(i), nn(i), ll(i), xx(i)
                nsize = nsize - sizeof(nf(i)) - sizeof(nn(i)) - sizeof(ll(i)) - sizeof(xx(i))
                rnor(i) = sqrt((dos * xx(i))**(2*nn(i)+1) / fact(2*nn(i)))
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
                if (myrank .eq. 0 .and. longoutput)  write(6,"(t11,i2,t15,i2,t20,e12.5,t36,i4)") nn(i), ll(i), xx(i), nf(i)
            enddo
        enddo
    else
        read(10) nprimitot
        nsize = nsize - sizeof(nprimitot)

!		Allocates memory for the basis set

        allocate(cfcontr(nprimitot), ipntprim(ncaps), ll(ncaps), lmaxc(ncen), ncontr(ncen), nf(ncaps), ngini(ncen), &
                ngfin(ncen), nprimit(ncaps), rnor(ncaps), xxg(nprimitot), stat = ierr)
        if (ierr .ne. 0) then
                write(6,"('Memory error when allocating cfcontr, ipntprim, ll, lmaxc, ncontr, nf, ngini, ngfin, &
                        &nprimit, rnor and xxg in processor ',i3)") myrank
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
                icarga = icarga+nprimit(knt)	! actualizes the index for loading primitives exponents and contraction coefficients
            enddo
        enddo
        if (myrank .eq. 0) write(6,"(/t22,'GTO Basis set',/t22,13('-'))")
        if (myrank .eq. 0 .and. longoutput) then
            icarga = 0
            knt = 0
            do ia = 1, ncen
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

!	Data of density representation
    read(10) lmaxexp
    nsize = nsize - sizeof(lmaxexp)
    if (lmaxrep .gt. lmaxexp) then
            if (myrank .eq. 0) write(6,"('lmaxrep = ', i3, ' greater than lmaxexp ', i3)") lmaxrep, lmaxexp
            if (myrank .eq. 0) write(6,"('takes lmaxrep = ',i3)") lmaxexp
            lmaxrep = lmaxexp
    endif
    lmtop = (lmaxexp+1)*(lmaxexp+1)

    if (myrank .eq. 0 .and. longoutput) write(6,"('lmaxexp = ', i2, ' nintervaj = ', i2)") lmaxexp, nintervaj

    allocate(icfposd(lmtop*nintervaj+1,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating icfposd in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    nsize = nsize - sizeof(icfposd(:,1)) * ncenbas

    allocate(xajustd(nintervaj,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating xajustd in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (longoutput) write(6,"('Estimated highest size of xajustd   = ', i15, ' bytes')") size(xajustd)
    nsize = nsize - sizeof(xajustd(:,1)) * ncenbas

    allocate(cfajust(nsize/8), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfajust in processor ',i3)") myrank
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
    icfposd = 0
    xajustd = cero
    cfajust = cero
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
    if (longoutput) write(6,"('Size of indintrv   = ', i15, ' bytes')") size(indintrv)

    r = cero
    interv = 1
    do i = 1, lenindintrv-1
        r = r + step
        if (r .gt. (rinterv(interv))) interv = interv + 1
        indintrv(i) = interv
    enddo
    indintrv(lenindintrv) = interv

!	Reads auxiliary integrals from file .dmqtv

    allocate(cfrint1(icfposd(lmtop*nintervaj+1,ncen)-1), cfrint2l2(icfposd(lmtop*nintervaj+1,ncen)-1), &
            QGacum(nintervaj*lmtop,ncen), Qgpart(nintervaj*lmtop), qpacum(nintervaj*lmtop,ncen), qppart(nintervaj*lmtop), &
            rlargo(ncen), rmultip(lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfrint1, cfrint2l2, QGacum, Qgpart, &
        &qpacum, qppart, rlargo and rmultip in processor ',i3)") myrank
        abort = 1
        return
    endif
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

    if (myrank .eq. 0 .and. longoutput) then
        write(6,"('Size of icfpos   = ', i15, ' bytes')") size(icfpos)
        write(6,"('Size of cfajust   = ', i15, ' bytes')") size(cfajust)
        write(6,"('Size of xajust   = ', i15, ' bytes')") size(xajust)
        write(6,"('Size of rmultip   = ', i15, ' bytes')") size(rmultip)
        write(6,"('Size of cfrint1   = ', i15, ' bytes')") size(cfrint1)
        write(6,"('Size of cfrint2l2 = ', i15, ' bytes')") size(cfrint2l2)
        write(6,"('Size of Qgpart    = ', i15, ' bytes')") size(Qgpart)
        write(6,"('Size of Qgacum    = ', i15, ' bytes')") size(Qgacum)
        write(6,"('Size of qppart    = ', i15, ' bytes')") size(qppart)
        write(6,"('Size of qpacum    = ', i15, ' bytes')") size(qpacum)
    endif

!	Determines the long-range radii and the highest l in the expansion for each interval
    allocate(lcorto(nintervaj,ncen), llargo(0:mxlargo,ncen), Qllargo(0:lmaxrep), stat = ierr )
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lcorto, llargo, auxll,  and Qllargo in processor ',i3)") myrank
        abort = 1
        return
    endif
!	long-range radii
    if (myrank .eq. 0 .and. .not. lexact) then
        write(6,"('Long-range threshold = ', e12.5)") umbrlargo
        write(6,"('Highest l in expansion of the potential = ', i2)") lmaxrep
    endif
    if (myrank .eq. 0 .and. longoutput) then
        write(6,"('Size of lcorto   = ', i15, ' bytes')") size(lcorto)
        write(6,"('Size of Qllargo   = ', i15, ' bytes')") size(Qllargo)
    endif
    allocate(umedpow(0:lmaxexp), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating umedpow in processor ',i3)") myrank
        abort = 1
        return
    endif
    umbrlargo2 = umbrlargo*umbrlargo
    umedpow(0) = uno							!
    do i = 1, lmaxexp							!
        umedpow(i) = umedpow(i-1) * umed			! 1 / 2^i
    enddo
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
                    pi4d2l1 = cuatro * pi * dosl1i(l)	!   4 * pi / (2l+1)
                    summ2 = cero
                    do m = -l, l
                        lm = lm + 1
                        if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo2) cycle
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
                    if (abs(summ2) .gt. umbrlargo2 .and. l .gt. lcorto(interv,ia)) lcorto(interv,ia) = l
                    if (abs(suml2-suml1) .gt. umbrlargo2) rlargo(ia) = rinterv(interv)
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
                if (suml1 .gt. umbrlargo2) then
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
!
!   ***************************************************************
!
   subroutine exactmespgrid
    USE MPI
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ierr, iuni, ix, iy, iz, knt, mpireal, nx, nxyz, nxyzrank, ny, nz
    real(KREAL) :: b2a, vetot, vntot, vtot, rx, ry, rz, x, y, z
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
        write(6,"('Memory error when allocating istav, iendv, ilenv and idispv in exactmespgrid in processor ',i3)") myrank
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
            write(6,"('In proc ', i3,' istart = ', i3, ' iend = ', i3, ' no points = ', i12, ' offset = ', i12)") &
                    i, istav(i), iendv(i), ilenv(i), idispv(i)
        enddo
    endif
!	Opens file for electrostatic potential tabulation
    nxyz = nx*ny*nz
    if (myrank .eq. 0) then
        allocate(array(nxyz), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating array in mespgrid in processor ',i3)") myrank
            abort = 1
            return
        endif
        iuni = 21
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"_exact-v.plt")
        array = cero
    endif
    nxyzrank = ilenv(myrank)
    allocate(arrayrank(nxyzrank), zlma((mxldst+1)**2), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating arrayrank and zlma in mespgrid in processor ',i3)") myrank
        abort = 1
        return
    endif
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
                    if (ierr .ne. 0) then
                        write(6,"('Exact potential for STO not yet available.')")
                        abort = 1
                        return
                    endif
                else
                    call GTOexactmesp(x, y, z, vntot, vetot, vtot)
                endif
                arrayrank(knt) = vtot
            enddo
        enddo
    enddo
    CALL MPI_GATHERV(arrayrank, nxyzrank, mpireal, array, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0 .and. myrank .eq. 0) then
        write(6,"('Error en MPI_GATHERV for array,   ierr = ',i7)") ierr
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
                call linea( 21, nx , arrayrank )
            enddo
        enddo
        close(21)
    endif
    if (myrank .eq. 0) deallocate(array)
    deallocate(arrayrank)
    deallocate(zlma)
    if (myrank .eq. 0) write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    return
    end
!	
!   ***************************************************************
!
   subroutine mespgrid
    USE MPI
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer status(MPI_STATUS_SIZE)
    integer(KINT) :: i, ia, ierr, iii, isel, iuni, ix, iy, iz, knt, mpireal, nx, nxyz, nxyzrank, ny, nz
    real(KREAL) :: b2a, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz, dV, vel, vnuc, vtot, rx, ry, rz, x, y, z
    real(KREAL), allocatable :: arrayrank(:), arraydxrank(:), arraydyrank(:), arraydzrank(:)
    real(KREAL), allocatable :: arraydxxrank(:), arraydxyrank(:), arraydxzrank(:), arraydyyrank(:)
    real(KREAL), allocatable :: arraydyzrank(:), arraydzzrank(:), arraylplrank(:)
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: arrayranksp(:)
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
    real(KREAL), allocatable :: arraydxx(:), arraydxy(:), arraydxz(:), arraydyy(:), arraydyz(:), arraydzz(:)
    mpireal = MPI_REAL8
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: arrayranksp(:)
    real(KREAL4), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
    real(KREAL4), allocatable :: arraydxx(:), arraydxy(:), arraydxz(:), arraydyy(:), arraydyz(:), arraydzz(:)
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
        write(6,"('Memory error when allocating istav, iendv, ilenv and idispv in mespgrid in processor ',i3)") myrank
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
!	Opens file for electrostatic potential tabulation
    nxyz = nx*ny*nz
    if (myrank .eq. 0) then
        allocate(array(nxyz), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating array in mespgrid in processor ',i3)") myrank
            abort = 1
            return
        endif
        iuni = 21
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v.plt")
        array = cero
!	Opens files for electrostatic potential gradient tabulation
        if (lgradient) then
            allocate(arraydx(nxyz), arraydy(nxyz), arraydz(nxyz), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraydx, arraydy and arraydz in mespgrid in processor ',i3)") myrank
                abort = 1
                return
            endif
            iuni = 23
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dx.pltd")
            iuni = 24
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dy.pltd")
            iuni = 25
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dz.pltd")
            arraydx = cero	! Array initialization
            arraydy = cero	! Array initialization
            arraydz = cero	! Array initialization
        endif
!	Opens files for electrostatic potential second derivatives tabulation
        if (lderiv2) then
            allocate(arraydxx(nxyz), arraydxy(nxyz), arraydxz(nxyz), arraydyy(nxyz), arraydyz(nxyz) &
                    , arraydzz(nxyz), stat = ierr)
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraydxx ... arraydzz in mespgrid in processor ',i3)") myrank
                abort = 1
                return
            endif
            iuni = 28
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dxx.pltd")
            iuni = 29
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dxy.pltd")
            iuni = 30
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dxz.pltd")
            iuni = 31
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dyy.pltd")
            iuni = 32
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dyz.pltd")
            iuni = 33
            call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-v-dzz.pltd")
            arraydxx = cero	! Array initialization
            arraydxy = cero	! Array initialization
            arraydxz = cero	! Array initialization
            arraydyy = cero	! Array initialization
            arraydyz = cero	! Array initialization
            arraydzz = cero	! Array initialization
        endif
    endif
    nxyzrank = ilenv(myrank)
    allocate(arrayrank(nxyzrank), arrayranksp(nxyzrank), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating arrayrank and arrayranksp in mespgrid in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (lgradient) then
        allocate(arraydxrank(nxyzrank), arraydyrank(nxyzrank), arraydzrank(nxyzrank), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arraydxrank, arraydyrank and arraydzrank in mespgrid in processor ' &
                    ,i3)") myrank
            abort = 1
            return
        endif
    endif
    if (lderiv2) then
        allocate(arraydxxrank(nxyzrank), arraydxyrank(nxyzrank), arraydxzrank(nxyzrank), arraydyyrank(nxyzrank) &
                , arraydyzrank(nxyzrank), arraydzzrank(nxyzrank), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arraydxxrank ... arraydzzrank in mespgrid in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif
    idimzlm = (lmaxexp+2)**2
    allocate(zlma(idimzlm), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma in mespgrid in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (lgradient) then
        allocate(zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating zlmadx, zlmady, zlmadz in mespgrid in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif
    if (lderiv2) then
        allocate(zlmadxx(idimzlm), zlmadxy(idimzlm), zlmadxz(idimzlm), zlmadyy(idimzlm), zlmadyz(idimzlm), &
                zlmadzz(idimzlm), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating zlmadxx, zlmadxy,... zlmadzz in mespgrid in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif
!	Grid tabulation
    kntlargo = 0
    kntcorto = 0
    knt = 0
    do iz =  1, nz
        if (mod(iz-1,nprocs) .ne. myrank) cycle
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix =1, nx
                knt = knt + 1
                x = xinf + (ix - 1) * dltx
                arrayrank(knt) = cero
                if (lgradient) then
                    arraydxrank(knt) = cero
                    arraydyrank(knt) = cero
                    arraydzrank(knt) = cero
                endif
                if (lderiv2) then
                    arraydxxrank(knt) = cero
                    arraydxyrank(knt) = cero
                    arraydxzrank(knt) = cero
                    arraydyyrank(knt) = cero
                    arraydyzrank(knt) = cero
                    arraydzzrank(knt) = cero
                endif
                do ia = 1, ncen
                    if (largo) then
                        call longmesp(ia, x, y, z, vnuc, vel, vtot, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz)
                    else
                        call mesp(ia, x, y, z, vnuc, vel, vtot, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz)
                    endif
                    arrayrank(knt) = arrayrank(knt) + vtot
                    if (lgradient) then
                        arraydxrank(knt) = arraydxrank(knt) + drvx
                        arraydyrank(knt) = arraydyrank(knt) + drvy
                        arraydzrank(knt) = arraydzrank(knt) + drvz
                    endif
                    if (lderiv2) then
                        arraydxxrank(knt) = arraydxxrank(knt) + dxx
                        arraydxyrank(knt) = arraydxyrank(knt) + dxy
                        arraydxzrank(knt) = arraydxzrank(knt) + dxz
                        arraydyyrank(knt) = arraydyyrank(knt) + dyy
                        arraydyzrank(knt) = arraydyzrank(knt) + dyz
                        arraydzzrank(knt) = arraydzzrank(knt) + dzz
                    endif
                enddo
            enddo
        enddo
    enddo
    arrayranksp = arrayrank
    CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, array, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
    if (ierr .ne. 0 .and. myrank .eq. 0) then
        write(6,"('Error en MPI_GATHERV for array,   ierr = ',i7)") ierr
        abort = 1
        return
    endif
    if (lgradient) then
        arrayranksp = arraydxrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydx, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydx,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        arrayranksp = arraydyrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydy,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        arrayranksp = arraydzrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydz,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
    endif
    if (lderiv2) then
        arrayranksp = arraydxxrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydxx, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
                    write(6,"('Error en MPI_GATHERV for arraydxx,   ierr = ',i7)") ierr
                    abort = 1
                    return
        endif
        arrayranksp = arraydxyrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydxy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydxy,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        arrayranksp = arraydxzrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydxz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydxz,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        arrayranksp = arraydyyrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydyy, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydyy,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        arrayranksp = arraydyzrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydyz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydyz,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        arrayranksp = arraydzzrank
        CALL MPI_GATHERV(arrayranksp, nxyzrank, mpireal, arraydzz, ilenv, idispv, mpireal, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for arraydzz,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
    endif
    if (re(nx)*re(ny)*re(nz) .le. (2.d9 / ncen)) then
        CALL MPI_REDUCE(kntlargo, kntlargotot, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_REDUCE for kntlargotot,   ierr = ',i7)") ierr
            write(6,"('kntlargotot = ', i12)") kntlargotot
            abort = 1
            return
        endif
        CALL MPI_REDUCE(kntcorto, kntcortotot, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_REDUCE for kntcortotot,   ierr = ',i7)") ierr
            write(6,"('kntcortotot = ', i12)") kntcortotot
            abort = 1
            return
        endif
    else
        if (myrank .eq. 0) then
            kntlargotot = kntlargo
            kntcortotot = kntcorto
            do i = 1, nprocs-1
                CALL MPI_RECV(kntlargo, 2, MPI_INTEGER8, i, 0, MPI_COMM_WORLD, status, ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error en MPI_RECV for kntlargo, ierr = ',i7)") ierr
                    write(6,"('kntlargo = ', i12)") kntlargo
                    abort = 1
                    return
                endif
                kntlargotot = kntlargotot + kntlargo
                CALL MPI_RECV(kntcorto, 2, MPI_INTEGER8, i, 1, MPI_COMM_WORLD, status, ierr)
                if (ierr .ne. 0) then
                    write(6,"('Error en MPI_RECV for kntcorto, ierr = ',i7)") ierr
                    write(6,"('kntcorto = ', i12)") kntcorto
                    abort = 1
                    return
                endif
                kntcortotot = kntcortotot + kntcorto
            enddo
        else
            CALL MPI_SEND(kntlargo, 2, MPI_INTEGER8, 0, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) then
                write(6,"('Error en MPI_SEND for kntlargo, ierr = ',i7)") ierr
                write(6,"('kntlargo = ', i12)") kntlargo
                abort = 1
                return
            endif
            CALL MPI_SEND(kntcorto, 2, MPI_INTEGER8, 0, 1, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) then
                write(6,"('Error en MPI_SEND for kntcorto, ierr = ',i7)") ierr
                write(6,"('kntcorto = ', i12)") kntcorto
                abort = 1
                return
            endif
        endif
    endif
    if (myrank .eq. 0) then
        do iz =  1, nz
            knt = idispv(mod(iz-1,nprocs)) + nx * ny * ((iz-1)/nprocs)
            do iy = 1, ny
                do ix =1, nx
                    knt = knt + 1
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
                enddo
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
            enddo
        enddo
! 	Deallocates arrays and closes the grid files
        close(21)
        if (lgradient) then
            close(23)
            close(24)
            close(25)
            deallocate(arraydx, arraydy,  arraydz, arraydxrank, arraydyrank,  arraydzrank)
        endif
        if (lderiv2) then
            close(28)
            close(29)
            close(30)
            close(31)
            close(32)
            close(33)
            deallocate(arraydxx, arraydxy, arraydxz, arraydyy, arraydyz, arraydzz)
            deallocate(arraydxxrank, arraydxyrank, arraydxzrank, arraydyyrank, arraydyzrank, arraydzzrank)
        endif
    endif
    deallocate(arrayrank, arrayranksp)
    if (myrank .eq. 0) deallocate(array)
    write(6,"('In proc ', i3, ' Number of long-range contributions  = ', i15)") myrank, kntlargo
    write(6,"('In proc ', i3, ' Number of short-range contributions = ', i15)") myrank, kntcorto
    if (myrank .eq. 0) write(6,"('Total number of long-range contributions  = ',6x, i15)") kntlargotot
    if (myrank .eq. 0) write(6,"('Total number of short-range contributions = ',6x, i15)") kntcortotot
    if (myrank .eq. 0) write(6,"('Total number of tabulated points = ', i15)") nx*ny*nz
    deallocate(zlma)
    if (allocated(zlmadx)) deallocate(zlmadx, zlmady, zlmadz)
    if (allocated(zlmadxx)) deallocate(zlmadxx, zlmadxy, zlmadxz, zlmadyy, zlmadyz, zlmadzz)
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
    real(KREAL) :: angnorm, bux, c1, c2, cosal, cosbet, cosga, cux, expbux, funcF0, potaux, potia, potiab
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
!                           fv(i) = rra**(2i) * Gamma[i+1/2,0,(xx(i1p)+xx(i2p))*rra**2] / (xx(i1p)+xx(i2p))**i
                        if (rra .lt. 1.e-20) then
                            fv = cero
                        else
                            call freq20(rra2,bux,fv)
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
!                           fv(i) = rrp**(2i) * Gamma[i+1/2,0,(xx(i1p)+xx(i2p))*rrp**2] / (xx(i1p)+xx(i2p))**i
                            if (rrp .lt. 1.e-20) then
                                fv = cero
                            else
                                call freq20(rrp2,bux,fv)
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
!
!   **************************************************************************************
!     Calculates the electrostatic potential from the represented density at point (x,y,z)
!
   subroutine mesp(ia, x, y, z, vnucl, vel, vtot, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz)
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    implicit none
    integer(KINT) :: i, ia, ierr, interv, icflm, j, jshft, kntlm, l, lm, lmtopot, ltop, m
    real(KREAL) :: aux, bux, dost, drvx, drvy, drvz, drvvlm, drv2vlm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, flm, fux
    real(KREAL) :: pi4exp, ra, ra2, rainv, ra2inv, rinta, rintb, sgn, t, tp, vnucl, vel, vlm, vqlm, vtot
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
            if (lderiv2) then
                dxx = vtope
                dxy = vtope
                dxz = vtope
                dyy = vtope
                dyz = vtope
                dzz = vtope
            endif
        endif
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
        if (lgradient) then
            drvx = -vnucl * xa * ra2inv
            drvy = -vnucl * ya * ra2inv
            drvz = -vnucl * za * ra2inv
            if (lderiv2) then
                dxx = -(re(3) * drvx * xa + vnucl) * ra2inv
                dxy = -re(3) * drvx * ya * ra2inv
                dxz = -re(3) * drvx * za * ra2inv
                dyy = -(re(3) * drvy * ya + vnucl) * ra2inv
                dyz = -re(3) * drvy * za * ra2inv
                dzz = -(re(3) * drvz * za + vnucl) * ra2inv
            endif
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

    if (lgradient) then		! Derivatives of the potential are also computed
        call derivzlm(ltop, idimzlm, zlma, zlmadx, zlmady, zlmadz)
        if (lderiv2) then
            call derivzlm(ltop, idimzlm, zlmadx, zlmadxx, zlmadxy, zlmadxz)
            call dzlm2y(ltop, idimzlm, zlmady, zlmadyy, zlmadyz)
            call dzlm2z(ltop, idimzlm, zlmadz, zlmadzz)
        endif
        if (ra .ge. rlargo(ia)) then  ! The point is in the long-range region (ra >= rlargo(ia) )
            kntlargo = kntlargo + 1
            aux = zn(ia)	! aux is the nuclear charge for l = 0 and zero otherwise
            do lm = 1, lmtopot
                vel = vel - rmultip(lm,ia) * zlma(lm) * ra2l1inv(lm)
                vlm = (aux - rmultip(lm,ia) ) * ra2l1inv(lm)
                vtot = vtot + vlm * zlma(lm)
                drvvlm = -d2l1(lm) * vlm * rainv
                drv2vlm = -(d2l1(lm)+uno) * drvvlm * rainv
                drvx = drvx + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
                drvy = drvy + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
                drvz = drvz + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
                if (lderiv2) then
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
                endif
                aux = cero
            enddo
        else		!     The point is in the short-range region (ra < rlargo(ia) )
            kntcorto = kntcorto + 1
            t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
            dost = t + t
            tcheb(0) = uno	! Chebyshev T  polynomials
            tcheb(1) = t
            do j = 2, mxlenpol-1
                    tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
            enddo
            aux = zn(ia)	! aux is the nuclear charge for l = 0 and zero otherwise
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
                drvvlm = -d2l1(lm) * vqlm * rainv 	! d2l1(lm) = (l+l+1)
                drv2vlm = -(d2l1(lm)+uno) * drvvlm * rainv + pi4exp * flm
                drvx = drvx + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
                drvy = drvy + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
                drvz = drvz + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
                if (lderiv2) then
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
                endif
                aux = cero
            enddo
        endif  ! End of test over long/short-range
    else		! Derivatives of the potential are not computed
        if (ra .ge. rlargo(ia)) then  ! The point is in the long-range region (ra >= rlargo(ia) )
            kntlargo = kntlargo + 1
            vel = - dot_product( rmultip(1:lmtopot,ia),zlma(1:lmtopot)*ra2l1inv(1:lmtopot) )
        else		!     The point is in the short-range region (ra < rlargo(ia) )
            kntcorto = kntcorto + 1
            t = dos * (ra - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
            dost = t + t
            tcheb(0) = uno	! Chebyshev T  polynomials
            tcheb(1) = t
            do j = 2, mxlenpol-1
                tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
            enddo
            do lm = 1, lmtopot
                if(abs(QGacum((nintervaj-1)*lmtop+lm,ia)) .lt. umbrlargo2) cycle
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
    endif
    return
    end
!
!   ********************************************************************
!     Calculates the long-range electrostatic potential on point (x,y,z)
!
   subroutine longmesp(ia, x, y, z, vnucl, vel, vtot, drvx, drvy, drvz, dxx, dxy, dxz, dyy, dyz, dzz)
    USE DAM320_D
    USE DAMPOT320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    implicit none
    integer(KINT) :: i, ia, ierr, kntlm, l, m, lm, lmtopot, ltop
    real(KREAL) :: aux, bux, drvx, drvy, drvz, drvvlm, drv2vlm
    real(KREAL) :: dxx, dxy, dxz, dyy, dyz, dzz, flm, fux
    real(KREAL) :: ra, ra2, rainv, ra2inv, rinta, rintb, vnucl, vel, vlm, vtot, x, xa, xadivra, y, ya, yadivra, z, za, zadivra
    real(KREAL) , parameter :: vtope = 1.d10	! To prevent infinity, if the point coincides with a nucleus, loads the value of parameter vtope
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
        return
    endif
    ra = sqrt(ra2)
    rainv = uno / ra
    ra2inv = rainv * rainv
    if (nintervaj .eq. 0) then
        vnucl = zn(ia) * rainv
        vtot = vnucl
        if (lgradient) then
            drvx = -vnucl * xa * ra2inv
            drvy = -vnucl * ya * ra2inv
            drvz = -vnucl * za * ra2inv
            if (lderiv2) then
                dxx = -(re(3) * drvx * xa + vnucl) * ra2inv
                dxy = -re(3) * drvx * ya * ra2inv
                dxz = -re(3) * drvx * za * ra2inv
                dyy = -(re(3) * drvy * ya + vnucl) * ra2inv
                dyz = -re(3) * drvy * za * ra2inv
                dzz = -(re(3) * drvz * za + vnucl) * ra2inv
            endif
        endif
        return
    endif
    ltop = llargo(min(int(ra),mxlargo),ia)
    lmtopot = (ltop+1)*(ltop+1)
    vnucl = zn(ia) * rainv
    xadivra = xa * rainv
    yadivra = ya * rainv
    zadivra = za * rainv
    zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
    zlma(2) = ya
    zlma(3) = za
    zlma(4) = xa
    ra2l1inv(1) = rainv
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
        bux = ra2l1inv(lm) * ra2inv
        do m = -l, l
            lm = lm + 1
            ra2l1inv(lm) = bux
        enddo
    enddo
    if (lgradient) then		! Derivatives of the potential are also computed
        call derivzlm(ltop, idimzlm, zlma, zlmadx, zlmady, zlmadz)
        if (lderiv2) then
            call derivzlm(ltop, idimzlm, zlmadx, zlmadxx, zlmadxy, zlmadxz)
            call dzlm2y(ltop, idimzlm, zlmady, zlmadyy, zlmadyz)
            call dzlm2z(ltop, idimzlm, zlmadz, zlmadzz)
        endif
    endif
    kntlargo = kntlargo + 1
    kntlm = 0
    aux = zn(ia)	! aux is the nuclear charge for l = 0 and zero otherwise
    do lm = 1, lmtopot
        vel = vel - rmultip(lm,ia) * zlma(lm) * ra2l1inv(lm)
        if (lgradient) then
            vlm = (aux - rmultip(lm,ia) ) * ra2l1inv(lm)
            drvvlm = -re(l+l+1) * vlm * rainv
            drv2vlm = -re(l+l+2) * drvvlm * rainv
            drvx = drvx + xadivra * drvvlm * zlma(lm) + vlm * zlmadx(lm)
            drvy = drvy + yadivra * drvvlm * zlma(lm) + vlm * zlmady(lm)
            drvz = drvz + zadivra * drvvlm * zlma(lm) + vlm * zlmadz(lm)
        endif
        if (lderiv2) then
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
        endif
        aux = cero
    enddo
    vtot = vnucl + vel
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
    if (lmax .eq. 0) return
    zlmadx(2) = cero	! Derivatives of the P spherical harmonics
    zlmadx(3) = cero
    zlmadx(4) = zlma(1)
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    if (lmax .eq. 1) return
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
    if (lmax .eq. 2) return
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
    if (lmax .eq. 0) return
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    if (lmax .eq. 1) return
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
    if (lmax .eq. 2) return
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
!	<</home/rafa/math/notebooks/pargamma25_D_2.nb>>
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
    if (iz .lt. 25) then
            go to (5,10,15,20,25) iz/5+1
5		continue
!	Interval  0 <= z <= 5   ( eps = 1.916D-17 )
!	polynomial approximation of F20(z)
            fv(20) = (4.8780487804878048D-2+z            &
                    *(-4.6511627906976616D-2+z*(2.2222222222219208D-2+z*(-7.0921985815316068D-3+z    &
                    *(1.7006802719618550D-3+z*(-3.2679738515298565D-4+z*(5.2410900461138484D-5+z     &
                    *(-7.2150056709346961D-6+z*(8.7022937944453839D-7+z*(-9.3413127864992509D-8+z    &
                    *(9.0341765026013756D-9+z*(-7.9477750647592279D-10+z*(6.4019156095112397D-11+z   &
                    *(-4.7232683275372959D-12+z*(3.1446647223867795D-13+z*(-1.7942609417785333D-14+z &
                    *(7.5594370023985866D-16+z*(-1.0859764932712809D-17+z*(-1.3577300442739411D-18+z &
                    *(1.1258133087359119D-19+z*(-3.1078282520470412D-21)))))))))))))))))))))
            go to 2000
10		continue
!		Interval  5 <= z <= 10   ( eps = 1.360D-16 )
!		polynomial approximation of Exp[x] * F20(z)
            fv(20) = ex * (4.8780523607920746D-2+z   &
                    *(2.2687949516701193D-3+z*(1.0089218810483327D-4+z*(4.2638000853699621D-6+z      &
                    *(1.8439794950170495D-7+z*(4.6223943411867869D-9+z*(6.5815410035004138D-10+z     &
                    *(-4.2969248943043529D-11+z*(5.3950760658262024D-12+z*(-3.4191372957330341D-13+z &
                    *(1.7368928116802409D-14+z*(-5.0438414248714736D-16+z*(8.1969351643649301D-18)))))))))))))
            go to 2000
15		continue
!		Interval  10 <= z <= 15   ( eps = 1.851D-18 )
!		polynomial approximation of Exp[x] * F20(z)
            fv(20) = ex * (4.8841694978314296D-2+z   &
                    *(2.1952612298984881D-3+z*(1.4197437389228347D-4+z*(-9.8786540053342240D-6+z     &
                    *(3.5368648978184568D-6+z*(-5.7450024066913360D-7+z*(7.5878158972779751D-8+z     &
                    *(-7.5109257899927147D-9+z*(5.7546785517165086D-10+z*(-3.3687901762469887D-11+z  &
                    *(1.4918741683111878D-12+z*(-4.8451102282376191D-14+z*(1.0969631911429229D-15+z  &
                    *(-1.5550170752478661D-17+z*(1.0620901118859267D-19)))))))))))))))
            go to 2000
20		continue
!		Interval  15 <= z <= 20   ( eps = 2.813D-20 )
!		polynomial approximation of Exp[x] * F20(z)
            fv(20) = ex * (9.6026515044316689D-2+z   &
                    *(-4.3166209794092771D-2+z*(2.0602239582622524D-2+z*(-5.7582054167261001D-3+z    &
                    *(1.1296518257164943D-3+z*(-1.6371999527316147D-4+z*(1.8161294813656418D-5+z     &
                    *(-1.5728398492941456D-6+z*(1.0752043107951713D-7+z*(-5.8231954957349261D-9+z    &
                    *(2.4915484933233244D-10+z*(-8.3383411682966011D-12+z*(2.1414249476219437D-13+z  &
                    *(-4.0839112977641997D-15+z*(5.4623096440148860D-17+z*(-4.5874914857129447D-19+z &
                    *(1.8295571542161169D-21)))))))))))))))))
            go to 2000
25		continue
!		Interval  20 <= z <= 25   ( eps = 4.030D-22 )
!		polynomial approximation of Exp[x] * F20(z)
            fv(20) = ex * (3.3224438260503371D1+z  &
                    *(-2.7603662264922665D1+z*(1.0855563769831582D1+z*(-2.6805503362623642D0+z       &
                    *(4.6582428660998887D-1+z*(-6.0503010428719655D-2+z*(6.0875879597865257D-3+z     &
                    *(-4.8521941163595912D-4+z*(3.1063234536669398D-5+z*(-1.6094515995913748D-6+z    &
                    *(6.7662136341545240D-8+z*(-2.3031265455170785D-9+z*(6.3019542936826522D-11+z    &
                    *(-1.3678219667355413D-12+z*(2.3040658855046203D-14+z*(-2.9078867616752635D-16+z &
                    *(2.5909392709711021D-18+z*(-1.4556134182691300D-20+z*(3.8853572116505235D-23)))))))))))))))))))
    else
            go to (35,45,55,65,75,85,95) (iz-25)/10+1
!		If z > 95: asymptotical expression
            fv(20) = 5.406242982335075D17 / sqrt(z)**41
            go to 2000
35   	continue
!		Interval  25 <= z <= 35   ( eps = 3.089D-18 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (5.8747691846974363D-5+y       &
                    *(9.7764465514125199D-1+y*(2.3413516830415170D1+y*(-5.6335933979497473D1+y       &
                    *(3.6483404782967801D4+y*(-1.4583609380850435D6+y*(6.1127998541254078D7+y        &
                    *(-1.6602879996934203D9+y*(3.5553799778962732D10+y*(-5.3096505923262290D11+y     &
                    *(5.6870782720347124D12+y*(-3.7315170959956712D13+y*(1.3243661404335882D14)))))))))))))
            go to 2000
45   	continue
!		Interval  35 <= z <= 45   ( eps = 1.547D-19 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (1.5989754165300373D-5+y       &
                    *(9.9312556185710981D-1+y*(2.0838150065672113D1+y*(2.0523027473807872D2+y        &
                    *(1.8285896271688041D4+y*(-5.3545710391705998D5+y*(2.5719413843168101D7+y        &
                    *(-6.1308064271060475D8+y*(1.1700760809463463D10+y*(-1.2321136681771141D11+y     &
                    *(7.6107340878398950D11)))))))))))
            go to 2000
55   	continue
!		Interval  45 <= z <= 55   ( eps = 4.160D-20 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (1.7437254480506169D-5+y       &
                    *(9.9235105114453346D-1+y*(2.0982789209186459D1+y*(1.9436556162147273D2+y        &
                    *(1.8182243066041124D4+y*(-4.5064318513701511D5+y*(1.8402159511206095D7+y        &
                    *(-2.8453651344225732D8+y*(3.1367646260747387D9)))))))))
            go to 2000
65   	continue
!		Interval  55 <= z <= 65   ( eps = 3.182D-20 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (4.4813294293896279D-5+y       &
                    *(9.8176159665610593D-1+y*(2.2650698969190298D1+y*(6.2837162193241420D1+y        &
                    *(2.2790603420436922D4+y*(-4.1574900536061811D5+y*(9.5985867522416429D6)))))))
            go to 2000
75   	continue
!		Interval  65 <= z <= 75   ( eps = 5.328D-20 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (2.0651610448975019D-4+y       &
                    *(9.3080337375748514D-1+y*(2.8584764207936303D1+y*(-2.1182416824959043D2+y       &
                    *(2.2555785714611651D4)))))
            go to 2000
85   	continue
!		Interval  75 <= z <= 85   ( eps = 1.660D-19 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (1.4546785444263069D-3+y       &
                    *(6.7563946245566736D-1+y*(4.1886052644826494D1)))
            go to 2000
95   	continue
!		Interval  85 <= z <= 95   ( eps = 6.308D-21 )
!	polynomial approximation of Exp[x] * (Fasympt20 - F20(z))
            y = uno / z
            fv(20) = 5.406242982335075D17 * sqrt(y)**41 - ex * (-3.7912464139565155D-3+y      &
                    *(1.6136443019529395D0))
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
!	***************************************************************
!	Subroutine cabecera: writes header for .plt files (binary)
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
    end
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

