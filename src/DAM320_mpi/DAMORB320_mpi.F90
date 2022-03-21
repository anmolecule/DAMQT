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
!	Program for generating grids for plotting orbitals
! 
! Parallel version with MPI
!
! Version of August 2014
!
  program DAMORB320_mpi
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE PARALELO
    implicit none
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL4), allocatable :: timeprocs(:)
    real(KREAL) :: aux, x, xmax, xmin, xyzmax, xyzmin, y, ymax, ymin, z, zmax, zmin
    integer(KINT) :: i, ierr, knt, norbs
    logical :: existe
    logical :: lsgbs, lsgbsden, lsgbsgz, lsgbsdengz
    logical :: lnamelist(6), ltimeprocs
    integer(KINT) :: inamelist(1)
    real(KREAL) :: rnamelist(9)
    namelist / options / dltu, dltv, dltx, dlty, dltz, filename, fileMOname, iorbsinp, iswindows, langstrom, lgradient,  &
            lgrid2d, lm2c, norbs, planeA, planeB, planeC, planecase, &
            uinf, usup, vinf, vsup, x_func_uv, xinf, xsup, y_func_uv, yinf, ysup, z_func_uv, zinf, zsup
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
    iswindows = .false.            ! .true. if running on a MS-windows system
    langstrom = .true.		! If .false. distances in bohr
    lgradient = .false.
    lm2c = .false.          ! If .true. read data from a calculation with m2c
    norbs = 1
    iorbsinp(1) = 1
    dltx = uno
    dlty = uno
    dltz = uno
    xinf = cero
    xsup = cero
    yinf = cero
    ysup = cero
    zinf = cero
    zsup = cero
    uinf = cero
    usup = cero
    dltu = uno
    vinf = cero
    vsup = cero
    dltv = uno
    filename = ""			! root file name for .plt and .pltd files
    fileMOname = ""		! file with Molecular orbitals coefficients
    lsgbsgz = .false.
    ldengz = .false.
    lsgbsdengz = .false.
!	End of namelist defaults

    ltimeprocs = .false.
    if (myrank .eq. 0) then
        read(5,OPTIONS)
        read(5,*) projectname
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
        existe = .false.
        if (len_trim(fileMOname) .eq.0) then
            inquire(file=projectname(1:i)//".orb", exist=existe)
            if (existe) then
                fileMOname = projectname//".orb"
            else
                inquire(file=projectname//".GAorba", exist=existe)
                if (existe) then
                    fileMOname = projectname//".GAorba"
                else
                    write(6,"('Cannot find a file with MO named ',a,'+ .orb or GAorba')") trim(projectname)
                    abort = 1
                endif
            endif
        else
            i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
            if (i .eq. 0) i = index(projectname,"/",.true.)	! In case of windows with MinGW directory separator is "/"
            fileMOname = projectname(1:i)//trim(fileMOname)
            inquire(file=fileMOname, exist=existe)
            if (.not. existe) then
                write(6,"('Cannot find a file with MO named ',a)") trim(fileMOname)
                abort = 1
            endif
        endif
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    if (myrank .eq. 0) then
!        Checks if files .sgbs or .sgbs.gz and .den or .den.gz exist. If yes, geometry, basis set and density will be read
!        with subroutine leedamm2c otherwise read them with subroutine leedatgen
        inquire(file=trim(projectname)//".sgbs", exist=lsgbs, iostat=ierr)
        if (ierr .ne. 0 .or. .not. lsgbs) then
            inquire(file=trim(projectname)//".sgbs.gz", exist=lsgbs, iostat=ierr)
            if (ierr .eq. 0 .and. lsgbs) then
                    call system ("gunzip "//trim(projectname)//".sgbs.gz")
                    lsgbsgz = .true.
            endif
        endif
        lden = .false.
        ldengz = .false.
        if (lsgbs) then
            inquire(file=trim(projectname)//".den", exist=lden, iostat=ierr)
            if (ierr .ne. 0 .or. .not. lden) then
                inquire(file=trim(projectname)//".den.gz", exist=lden, iostat=ierr)
                if (ierr .eq. 0 .and. lden) then
                    if(myrank .eq. 0) call system ("gunzip "//trim(projectname)//".den.gz")
                    ldengz = .true.
                endif
            endif
        endif
        lsgbsden = .false.
        lsgbsdengz = .false.
        inquire(file=trim(projectname)//".sgbsden", exist=lsgbsden, iostat=ierr)
        if (ierr .ne. 0 .or. .not. lsgbsden) then
            inquire(file=trim(projectname)//".sgbsden.gz", exist=lsgbsden, iostat=ierr)
            if (ierr .eq. 0 .and. lsgbsden) then
                call system ("gunzip "//trim(projectname)//".sgbsden.gz")
                lsgbsdengz = .true.
            endif
        endif
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    CALL MPI_BCAST(projectname,len(projectname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(fileMOname,len(fileMOname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    call consta
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    if (myrank .eq. 0) then
        if (fileMOname(len_trim(fileMOname)-3:len_trim(fileMOname)) == ".orb" .or. &
            fileMOname(len_trim(fileMOname)-6:len_trim(fileMOname)-1) == ".SLorb") then
            lsto = .true.
            existe = .true.
        elseif (fileMOname(len_trim(fileMOname)-6:len_trim(fileMOname)-1) == ".GAorb") then
            lsto = .false.
            existe = .true.
        else
            existe = .false.
            i = index(fileMOname,".",.true.)	! Checks position of last directory name separator
            write(6,"('MO filename extension:  ',a,'  not acceptable. Must be .orb, .SLorba, .SLorbb, .GAorba or .GAorbb ')") &
                    trim(fileMOname(i:len_trim(fileMOname)))
            abort = 1
        endif
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    if (myrank .eq. 0) then
        if (lsto) then
            if (lm2c) then
                call leedatm2c
            else
                call leedatSTOgen(lsgbsden)
            endif
        else
            call leedatgauss
        endif
        allocate(iorbs(norbs), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating iorbs. Stop')")
            abort = 1
        endif
        knt = 0
        do i = 1, norbs
            if (iorbsinp(i) .le. nbas) then	!	Checks that the orbital indices are in the range 1:nbas
                knt = knt + 1
                iorbs(knt) = iorbsinp(i)
            endif
        enddo
        norbs = knt
        inamelist = (/ norbs /)
        rnamelist = (/ xinf, xsup, dltx, yinf, ysup, dlty, zinf, zsup, dltz /)
        allocate(timeprocs(2*nprocs), stat = ierr)
        if (ierr .eq. 0) then
            ltimeprocs = .true.
            timeprocs = 0.
        else
            write(6,"('WARNING: Memory error when allocating timeprocs, ierr =  ',i5)") ierr
            ltimeprocs = .false.
        endif
        lnamelist = (/ langstrom, lgradient, lsto, ltimeprocs, lm2c, lsgbsden /)
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
            call error(1,'Stop')
    endif

    CALL MPI_BCAST(lnamelist,6,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(inamelist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(rnamelist,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    if (myrank .ne. 0) then
        norbs = inamelist(1)
        langstrom = lnamelist(1); lgradient = lnamelist(2); lsto = lnamelist(3); ltimeprocs = lnamelist(4);
        lm2c = lnamelist(5); lsgbsden = lnamelist(6)
        xinf = rnamelist(1); xsup = rnamelist(2); dltx = rnamelist(3)
        yinf = rnamelist(4); ysup = rnamelist(5); dlty = rnamelist(6)
        zinf = rnamelist(7); zsup = rnamelist(8); dltz = rnamelist(9)
        allocate(iorbs(norbs), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating iorbs in processor ',i3)") myrank
            abort = 1
        endif
        if (lsto) then
            if (lm2c) then
                call leedatm2c
            else
                call leedatSTOgen(lsgbsden)
            endif
        else
            call leedatgauss
        endif
    endif
    if (myrank .eq. 0) then
        if (ldengz) then    ! restores files back to their original gzipped status
            call system ("gzip "//trim(projectname)//".den")
        endif
        if (lsgbsgz) then
            call system ("gzip "//trim(projectname)//".sgbs")
        endif
        if (lsgbsdengz) then
            call system ("gzip "//trim(projectname)//".sgbsden")
        endif
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    CALL MPI_BCAST(iorbs,norbs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    allocate(corbs(nbas,norbs), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating corbs in processor ',i3)") myrank
        abort = 1
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif

    if ((xsup-xinf) .eq. cero .or. (ysup-yinf) .eq. cero .or. (zsup-zinf) .eq. cero) then	!	Default grid
        xmin = cero
        xmax = cero
        ymin = cero
        ymax = cero
        zmin = cero
        zmax = cero
        do i = 1, ncen
            if (rcen(1,i) .lt. xmin) xmin = rcen(1,i)
            if (rcen(1,i) .gt. xmax) xmax = rcen(1,i)
            if (rcen(2,i) .lt. ymin) ymin = rcen(2,i)
            if (rcen(2,i) .gt. ymax) ymax = rcen(2,i)
            if (rcen(3,i) .lt. zmin) zmin = rcen(3,i)
            if (rcen(3,i) .gt. zmax) zmax = rcen(3,i)
        enddo
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
            xinf = umed * (re(3) * xmin - xmax)
            xsup = umed * (re(3) * xmax - xmin)
            dltx = (xsup-xinf)/65
        endif
        if ((ysup-yinf) .eq. cero) then
            yinf = umed * (re(3) * ymin - ymax)
            ysup = umed * (re(3) * ymax - ymin)
            dlty = (ysup-yinf)/65
        endif
        if ((zsup-zinf) .eq. cero) then
            zinf = umed * (re(3) * zmin - zmax)
            zsup = umed * (re(3) * zmax - zmin)
            dltz = (zsup-zinf)/65
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

!	reads orbitals and generates grids with orbitals
    if (lsto) then
        if (lm2c) then
            call leeorbSMILES(norbs)
        else
            call leeorbgen(norbs)
        endif
        call pltorbSMILES(norbs)
    else
        call leeorbGTO(norbs)
        call pltorbGTO(norbs)
    endif
    CALL MPI_REDUCE(abort,abortroot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(abortroot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    if (abortroot .gt. 0) then
        call error(1,'Stop')
    endif
    tiempo = dtime(tarray)
!    write(6,"(1x,'Timing in seconds of processor ', i2, ' (user, system, total):',5x,'(',e12.5,',',e12.5,',',e12.5')')") &
!            myrank, tarray(1), tarray(2), tarray(1)+tarray(2)
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
!   *************************************************************
!
!  	Reads the coefficients of Slater MO  and returns them in c
!
  subroutine leeorbSMILES(norbs)
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE PARALELO
    implicit none
    real(KREAL), allocatable :: caux(:,:)
    integer(KINT) :: i, ierr, j, nbasis, norbs
    allocate(caux(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating caux in processor ',i3)") myrank
        abort = 1
        return
    endif

    open(13,file=trim(fileMOname),form='unformatted', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(fileMOname), myrank
        abort = 1
        return
    endif

    read(13) nbasis, ((caux(i,j),i=1,nbas),j=1,nbas)
    close(13)
    if (myrank .eq. 0) then
        if (nbas .ne. nbasis) then
            write(6,"('Wrong number of basis functions. Stop')")
            abort = 1
            return
        endif
    endif

    do i = 1, norbs
        corbs(1:nbas,i) = caux(1:nbas,iorbs(i))
    enddo
    return
    deallocate(caux)
    end
!
!   *************************************************************
!
!  	Reads the coefficients of Slater MO  and returns them in c
!
subroutine leeorbgen(norbs)
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE PARALELO
    implicit none
    real(KREAL), allocatable :: caux(:,:)
    integer(KINT) :: i, ierr, j, nbasis, nmos, norbs, nvoid

    open(13,file=trim(fileMOname),form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(fileMOname), myrank
        abort = 1
    endif
    read(13,*) nbasis, nvoid, nmos

    allocate(caux(nbasis,nmos), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating caux in processor ',i3)") myrank
        abort = 1
        return
    endif

    read(13,*) ((caux(i,j),i=1,nbasis),j=1,nmos)
    close(13)
    i = 0
    do j = 1, norbs
        if (iorbs(j) .le. nmos) then
            i = i + 1
            corbs(1:nbas,i) = caux(1:nbas,iorbs(j))
        endif
    enddo
    deallocate(caux)
    return
    end
!
!   *************************************************************
!
!  	Reads the coefficients of Slater MO  and returns them in c
!
  subroutine leeorbGTO(norbs)
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE PARALELO
    implicit none
    real(KREAL), allocatable :: caux(:,:)
    integer(KINT) :: i, ierr, j, nbasis, numorb, norbs, nvoid
    allocate(caux(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating caux in processor ',i3)") myrank
        abort = 1
        return
    endif

    open(13,file=trim(fileMOname),form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, ' in processor ',i3)") trim(fileMOname), myrank
        abort = 1
    endif

    read(13,*) nbas, nvoid, numorb
    do j = 1, numorb
        read(13,*) (caux(i,j),i=1,nbas)
    enddo
    close(13)
    do i = 1, norbs
        corbs(1:nbas,i) = caux(1:nbas,iorbs(i))
    enddo
    return
    deallocate(caux)
    end
!                                                                               
! *******************************************************************
!                                                                               
  subroutine pltorbSMILES(norbs)
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE PARALELO
    implicit none
    real(KREAL4), allocatable :: array(:), arrayranksp(:)
    real(KREAL), allocatable :: arrayrank(:,:), arraydxrank(:,:), arraydyrank(:,:), arraydzrank(:,:)
    real(KREAL) :: b2a, x, y, z
    integer(KINT) :: iaux, ierr, iorb, iuni, ix, iy, iz, norbs, nx, nxyz, nxyzrank, ny, nz
    character(2) faux
    real(KREAL4) :: x1a, x2a, y1a, y2a, z1a, z2a
    integer(KINT) :: i, knt
    allocate(orbv(norbs), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating orbv in processor ',i3)") myrank
        abort = 1
        return
    endif
    allocate(chi(nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating chi in processor ',i3)") myrank
        abort = 1
        return
    endif
    allocate(zlma((lmaxbase+1)*(lmaxbase+1)), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma in processor ',i3)") myrank
        abort = 1
        return
    endif

    if (lgradient) then
        allocate(orbvdx(nbas), orbvdy(nbas), orbvdz(nbas), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating orbvdx, orbvdy, orbvdz in processor ',i3)") myrank
            abort = 1
            return
        endif
        allocate(chidx(nbas), chidy(nbas), chidz(nbas), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating chidx, chidy, chidz in processor ',i3)") myrank
            abort = 1
            return
        endif
        allocate(zlmadx((lmaxbase+1)*(lmaxbase+1)), zlmady((lmaxbase+1)*(lmaxbase+1)), zlmadz((lmaxbase+1)*(lmaxbase+1)) &
                , stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating zlmadx, zlmady, zlmadz in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif

    b2a = 0.5291772d0
    nx = (xsup-xinf) / dltx + 1
    ny = (ysup-yinf) / dlty + 1
    nz = (zsup-zinf) / dltz + 1
    x1a = xinf * b2a
    x2a = xsup * b2a
    y1a = yinf * b2a
    y2a = ysup * b2a
    z1a = zinf * b2a
    z2a = zsup * b2a
!	Determines the grid points for tabulation assigned to each processor:
!		istav(myrank): starting index iz assigned to processor myrank
!		iendv(myrank): ending index iz assigned to processor myrank
    allocate(istav(0:nprocs-1), iendv(0:nprocs-1), ilenv(0:nprocs-1), idispv(0:nprocs-1), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating istav in gridpot in processor ',i3)") myrank
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
    if (myrank .eq. 0)  then
        write(6,"(/,'Grid points assignment')")
        do i = 0, nprocs-1
            write(6,"('In proc ', i3,' istart = ', i3, ' iend = ', i3, ' no points = ', i12, ' disp = ', i12)") &
                    i, istav(i), iendv(i), ilenv(i), idispv(i)
        enddo
    endif
    nxyz = nx*ny*nz
    if (myrank .eq. 0) then
        allocate(array(nxyz), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating array in gridpot in processor ',i3)") myrank
            abort = 1
            return
        endif
        iuni = 10
        do iorb = 1, norbs
            iuni = iuni + 1
            write(faux,'(i2.2)') iorbs(iorb)
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//".plt")
        enddo
!	Opens files for orbitals gradients tabulation
        if (lgradient) then
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraydx, arraydy and arraydz in gridpot in processor ',i3)") myrank
                abort = 1
                return
            endif
            do iorb = 1, norbs
                iuni = iuni + 1
                write(faux,'(i2.2)') iorbs(iorb)
                call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dx.pltd")
            enddo
            do iorb = 1, norbs
                iuni = iuni + 1
                write(faux,'(i2.2)') iorbs(iorb)
                call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dy.pltd")
            enddo
            do iorb = 1, norbs
                iuni = iuni + 1
                write(faux,'(i2.2)') iorbs(iorb)
                call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dz.pltd")
            enddo
        endif
    endif
    nxyzrank = ilenv(myrank)
    allocate(arrayrank(nxyzrank,norbs), arrayranksp(nxyzrank), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating arrayrank and arrayranksp in gridpot in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (lgradient) then
        allocate(arraydxrank(nxyzrank,norbs), arraydyrank(nxyzrank,norbs), arraydzrank(nxyzrank,norbs), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arraydxrank, arraydyrank and arraydzrank in gridpot in processor ' &
                    ,i3)") myrank
            abort = 1
            return
        endif
    endif
!
!  Generates grids
!
    knt = 0
    do iz = istav(myrank), iendv(myrank)
        z = zinf + (iz-1)*dltz
        do iy = 1 , ny
            y = yinf + (iy-1)*dlty
            do ix = 1 , nx
                x = xinf + (ix-1)*dltx
                knt = knt + 1
                call orbitalesSMILES( norbs, x , y , z )
                do iorb = 1, norbs
                    arrayrank(knt,iorb) = orbv(iorb)
                end do
                if (lgradient) then
                    do iorb = 1, norbs
                        arraydxrank(knt,iorb) = orbvdx(iorb)
                        arraydyrank(knt,iorb) = orbvdy(iorb)
                        arraydzrank(knt,iorb) = orbvdz(iorb)
                    end do
                endif
            end do
        end do
    end do
    iuni = 10
    do iorb = 1, norbs
        iuni = iuni + 1
        arrayranksp = arrayrank(1:nxyzrank,iorb)
        CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for array,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        if (myrank .eq. 0) write(iuni) array
    enddo
    if (lgradient) then
        do iorb = 1, norbs
            iuni = iuni + 1
            arrayranksp = arraydxrank(1:nxyzrank,iorb)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error en MPI_GATHERV for arraydxrank,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (myrank .eq. 0) write(iuni) array
        enddo
        do iorb = 1, norbs
            iuni = iuni + 1
            arrayranksp = arraydyrank(1:nxyzrank,iorb)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error en MPI_GATHERV for arraydyrank,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (myrank .eq. 0) write(iuni) array
        enddo
        do iorb = 1, norbs
            iuni = iuni + 1
            arrayranksp = arraydzrank(1:nxyzrank,iorb)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error en MPI_GATHERV for arraydzrank,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (myrank .eq. 0) write(iuni) array
        enddo
    endif
!
!  Close files
!
    if (myrank .eq. 0) then
        iuni = 10
        iaux = norbs
        if (lgradient) iaux = 4*iaux
        do iorb = 1 , iaux
            close(iuni+iorb)
        end do
    endif
    return
    end
!                                                                               
! *******************************************************************
!                                                                               
  subroutine pltorbGTO(norbs)
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE PARALELO
    implicit none
    real(KREAL4), allocatable :: array(:), arrayranksp(:)
    real(KREAL), allocatable :: arrayrank(:,:), arraydxrank(:,:), arraydyrank(:,:), arraydzrank(:,:)
    real(KREAL) :: b2a, x, y, z
    integer(KINT) :: iaux, ierr, iorb, iuni, ix, iy, iz, norbs, nx, nxyz, nxyzrank, ny, nz
    character(2) faux
    real(KREAL4) :: x1a, x2a, y1a, y2a, z1a, z2a
    integer(KINT) :: i, knt
    allocate(orbv(norbs), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating orbv in processor ',i3)") myrank
        abort = 1
        return
    endif
    allocate(chi(nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating chi in processor ',i3)") myrank
        abort = 1
        return
    endif
    allocate(zlma((lmaxbase+1)*(lmaxbase+1)), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zlma in processor ',i3)") myrank
        abort = 1
        return
    endif

    if (lgradient) then
        allocate(orbvdx(nbas), orbvdy(nbas), orbvdz(nbas), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating orbvdx, orbvdy, orbvdz in processor ',i3)") myrank
            abort = 1
            return
        endif
        allocate(chidx(nbas), chidy(nbas), chidz(nbas), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating chidx, chidy, chidz in processor ',i3)") myrank
            abort = 1
            return
        endif
        allocate(zlmadx((lmaxbase+1)*(lmaxbase+1)), zlmady((lmaxbase+1)*(lmaxbase+1)), zlmadz((lmaxbase+1)*(lmaxbase+1)) &
                , stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating zlmadx, zlmady, zlmadz in processor ',i3)") myrank
            abort = 1
            return
        endif
    endif

    b2a = 0.5291772d0
    nx = (xsup-xinf) / dltx + 1
    ny = (ysup-yinf) / dlty + 1
    nz = (zsup-zinf) / dltz + 1
    x1a = xinf * b2a
    x2a = xsup * b2a
    y1a = yinf * b2a
    y2a = ysup * b2a
    z1a = zinf * b2a
    z2a = zsup * b2a
!	Determines the grid points for tabulation assigned to each processor:
!		istav(myrank): starting index iz assigned to processor myrank
!		iendv(myrank): ending index iz assigned to processor myrank
    allocate(istav(0:nprocs-1), iendv(0:nprocs-1), ilenv(0:nprocs-1), idispv(0:nprocs-1), stat = ierr)
    if (ierr .ne. 0) then
            write(6,"('Memory error when allocating istav in gridpot in processor ',i3)") myrank
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
    if (myrank .eq. 0)  then
        write(6,"(/,'Grid points assignment')")
        do i = 0, nprocs-1
            write(6,"('In proc ', i3,' istart = ', i3, ' iend = ', i3, ' no points = ', i12, ' disp = ', i12)") &
                    i, istav(i), iendv(i), ilenv(i), idispv(i)
        enddo
    endif
    nxyz = nx*ny*nz
    if (myrank .eq. 0) then
        allocate(array(nxyz), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating array in gridpot in processor ',i3)") myrank
            abort = 1
            return
        endif
        iuni = 10
        do iorb = 1, norbs
            iuni = iuni + 1
            write(faux,'(i2.2)') iorbs(iorb)
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//".plt")
        enddo
!	Opens files for orbitals gradients tabulation
        if (lgradient) then
            if (ierr .ne. 0) then
                write(6,"('Memory error when allocating arraydx, arraydy and arraydz in gridpot in processor ',i3)") myrank
                abort = 1
                return
            endif
            do iorb = 1, norbs
                iuni = iuni + 1
                write(faux,'(i2.2)') iorbs(iorb)
                call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dx.pltd")
            enddo
            do iorb = 1, norbs
                iuni = iuni + 1
                write(faux,'(i2.2)') iorbs(iorb)
                call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dy.pltd")
            enddo
            do iorb = 1, norbs
                iuni = iuni + 1
                write(faux,'(i2.2)') iorbs(iorb)
                call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dz.pltd")
            enddo
        endif
    endif
    nxyzrank = ilenv(myrank)
    allocate(arrayrank(nxyzrank,norbs), arrayranksp(nxyzrank), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating arrayrank and arrayranksp in gridpot in processor ',i3)") myrank
        abort = 1
        return
    endif
    if (lgradient) then
        allocate(arraydxrank(nxyzrank,norbs), arraydyrank(nxyzrank,norbs), arraydzrank(nxyzrank,norbs), stat = ierr)
        if (ierr .ne. 0) then
            write(6,"('Memory error when allocating arraydxrank, arraydyrank and arraydzrank in gridpot in processor ' &
                    ,i3)") myrank
            abort = 1
            return
        endif
    endif
    if (abort .gt. 0) return
!
!  Generates grids
!
    knt = 0
    do iz = istav(myrank), iendv(myrank)
        z = zinf + (iz-1)*dltz
        do iy = 1 , ny
            y = yinf + (iy-1)*dlty
            do ix = 1 , nx
                x = xinf + (ix-1)*dltx
                knt = knt + 1
                call orbitalesGTO( norbs, x , y , z )
                do iorb = 1, norbs
                    arrayrank(knt,iorb) = orbv(iorb)
                end do
                if (lgradient) then
                    do iorb = 1, norbs
                        arraydxrank(knt,iorb) = orbvdx(iorb)
                        arraydyrank(knt,iorb) = orbvdy(iorb)
                        arraydzrank(knt,iorb) = orbvdz(iorb)
                    end do
                endif
            end do
        end do
    end do
    iuni = 10
    do iorb = 1, norbs
        iuni = iuni + 1
        arrayranksp = arrayrank(1:nxyzrank,iorb)
        CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
        if (ierr .ne. 0 .and. myrank .eq. 0) then
            write(6,"('Error en MPI_GATHERV for array,   ierr = ',i7)") ierr
            abort = 1
            return
        endif
        if (myrank .eq. 0) write(iuni) array
    enddo
    if (lgradient) then
        do iorb = 1, norbs
            iuni = iuni + 1
            arrayranksp = arraydxrank(1:nxyzrank,iorb)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error en MPI_GATHERV for arraydxrank,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (myrank .eq. 0) write(iuni) array
        enddo
        do iorb = 1, norbs
            iuni = iuni + 1
            arrayranksp = arraydyrank(1:nxyzrank,iorb)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error en MPI_GATHERV for arraydyrank,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (myrank .eq. 0) write(iuni) array
        enddo
        do iorb = 1, norbs
            iuni = iuni + 1
            arrayranksp = arraydzrank(1:nxyzrank,iorb)
            CALL MPI_GATHERV(arrayranksp, nxyzrank, MPI_REAL4, array, ilenv, idispv, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0 .and. myrank .eq. 0) then
                write(6,"('Error en MPI_GATHERV for arraydzrank,   ierr = ',i7)") ierr
                abort = 1
                return
            endif
            if (myrank .eq. 0) write(iuni) array
        enddo
    endif
!
!  Close files
!
    if (myrank .eq. 0) then
        iuni = 10
        iaux = norbs
        if (lgradient) iaux = 4*iaux
        do iorb = 1 , iaux
            close(iuni+iorb)
        end do
    endif
    return
    end
!
!   *******************************************************
!
  subroutine orbitalesSMILES(norbs, x , y , z)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL) :: aux, ra, ra2, rai, raix, raiy, raiz, ranl, ranlm1, x, xa, y, ya, z, za
    integer(KINT) :: iat, icapa, ix, knt, l, m, n, norbs
!
!  Computes basis functions in point (x,y,z)
!
    do iat = 1 , ncen
        if (ngini(iat) .le. 0) cycle
        xa = x - rcen(1,iat)
        ya = y - rcen(2,iat)
        za = z - rcen(3,iat)
        ra2= xa*xa + ya*ya + za*za
        ra = sqrt(ra2)
        if (lgradient) then
            if (ra .ne. cero) then
                rai = uno / ra
                raix = rai * x
                raiy = rai * y
                raiz = rai * z
            else
                rai = cero	! This prevents dividing by zero and anihilates null contributions to derivatives
                raix = uno
                raiy = uno
                raiz = uno
            endif
        endif
        call solido(xa, ya, za)
        do icapa = ngini(iat) , ngfin(iat)
            n = nn(icapa)
            l = ll(icapa)
            aux = rnor(icapa) * exp(-xx(icapa)*ra) * ra**(n-l-1)
            knt = l*l
            do m = -l , l
                knt = knt + 1
                ix = nf(icapa)+l+m
                chi(ix) = aux * zlma(knt) * alm(knt)
                if (lgradient) then
                    chidx(ix) = aux * alm(knt) * ((re(n-l-1)*rai-xx(icapa)) * raix * zlma(knt) + zlmadx(knt))
                    chidy(ix) = aux * alm(knt) * ((re(n-l-1)*rai-xx(icapa)) * raiy * zlma(knt) + zlmady(knt))
                    chidz(ix) = aux * alm(knt) * ((re(n-l-1)*rai-xx(icapa)) * raiz * zlma(knt) + zlmadz(knt))
                endif
            end do
        end do
    end do
!
!  calculo los orbitales en x,y,z
!
    orbv(1:norbs) = matmul( chi(1:nbas) , corbs(1:nbas,1:norbs) )
    if (lgradient) then
        orbvdx(1:norbs) = matmul( chidx(1:nbas) , corbs(1:nbas,1:norbs) )
        orbvdy(1:norbs) = matmul( chidy(1:nbas) , corbs(1:nbas,1:norbs) )
        orbvdz(1:norbs) = matmul( chidz(1:nbas) , corbs(1:nbas,1:norbs) )
    endif
    return
    end
!
!   *******************************************************
!
  subroutine orbitalesGTO(norbs, x , y , z)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    USE GAUSS
    implicit none
    real(KREAL) :: aux, frad, fradder, ra, ra2, x, xa, y, ya, z, za
    integer(KINT) :: i1, i1p, iat, ix, knt, l, m, n, norbs
!
!  Computes basis functions in point (x,y,z)
!
    do iat = 1 , ncen
        if (ngini(iat) .le. 0) cycle
        xa = x - rcen(1,iat)
        ya = y - rcen(2,iat)
        za = z - rcen(3,iat)
        ra2= xa*xa + ya*ya + za*za
        ra = sqrt(ra2)
        call solido(xa, ya, za)
        do i1 = ngini(iat), ngfin(iat)
            frad = cero
            fradder = cero
            do i1p = ipntprim(i1), ipntprim(i1)+nprimit(i1)-1
                aux = cfcontr(i1p) * exp(-xxg(i1p)*ra2)
                frad = frad + aux
                if (lgradient) fradder = fradder - dos * aux * xxg(i1p)
            enddo
            frad = rnor(i1) * frad
            if (lgradient) fradder = rnor(i1) * fradder
            do m = -ll(i1), ll(i1)
                chi(nf(i1)+ll(i1)+m) = frad * alm(ll(i1)*(ll(i1)+1)+m+1) * zlma(ll(i1)*(ll(i1)+1)+m+1)
                if (lgradient) then
                    chidx(nf(i1)+ll(i1)+m) = (fradder*x* zlma(ll(i1)*(ll(i1)+1)+m+1) + frad* zlmadx(ll(i1)*(ll(i1)+1)+m+1)) &
                            * alm(ll(i1)*(ll(i1)+1)+m+1)
                    chidy(nf(i1)+ll(i1)+m) = (fradder*y* zlma(ll(i1)*(ll(i1)+1)+m+1) + frad* zlmady(ll(i1)*(ll(i1)+1)+m+1)) &
                            * alm(ll(i1)*(ll(i1)+1)+m+1)
                    chidz(nf(i1)+ll(i1)+m) = (fradder*z* zlma(ll(i1)*(ll(i1)+1)+m+1) + frad* zlmadz(ll(i1)*(ll(i1)+1)+m+1)) &
                            * alm(ll(i1)*(ll(i1)+1)+m+1)
                endif
            enddo
        enddo
    end do
!
!  calculo los orbitales en x,y,z
!
    orbv(1:norbs) = matmul( chi(1:nbas) , corbs(1:nbas,1:norbs) )
    if (lgradient) then
        orbvdx(1:norbs) = matmul( chidx(1:nbas) , corbs(1:nbas,1:norbs) )
        orbvdy(1:norbs) = matmul( chidy(1:nbas) , corbs(1:nbas,1:norbs) )
        orbvdz(1:norbs) = matmul( chidz(1:nbas) , corbs(1:nbas,1:norbs) )
    endif
    return
    end
!
!	***************************************************************
!	Subroutine cabecera: writes head for .plt files (binary)
!
  subroutine cabecera(nx, ny, nz, xinf, xsup, yinf, ysup, zinf, zsup, iuni, s)
    USE DAM320_D
    implicit none
    integer(KINT) :: i
    integer(KINT) :: nx,ny,nz,iuni,ns,iaux(0:2)
    real(KREAL4) :: xinf,xsup,yinf,ysup,zinf,zsup,v(0:5)
    character*(*) :: s
!	If the compiler is other than INTEL's, uses the OPEN
!	sentence for stream files according to Fortran 2003 standard
#if _WIN32
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#elif __INTEL_COMPILER
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#else
    open (unit=iuni, file=s, form='unformatted', access='stream')
#endif
    write(iuni) 3_4, 2_4
    iaux(0) = nz
    iaux(1) = ny
    iaux(2) = nx
    write(iuni) iaux(0), iaux(1), iaux(2)
    v(0) = zinf
    v(1) = zsup
    v(2) = yinf
    v(3) = ysup
    v(4) = xinf
    v(5) = xsup
    write(iuni) (v(i), i = 0, 5)
    return
    END
!
!   ***************************************************************
!
  subroutine consta
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    real(KREAL) :: aux
    integer(KINT) :: i, ierr, l, m, ma
    pi = acos(-uno)
    raizpi = sqrt(pi)
    re(0) = cero
    ri(0) = 1.d300
    dosl1(0) = uno
    do i = 1, mxreal
        re(i) = re(i-1) + uno        ! dfloat(i)
        re(-i) = -re(i)
        ri(i) = uno / re(i)       	! uno / dfloat(i)
        ri(-i) = -ri(i)
        dosl1(i) = re(i) + re(i) + uno	! dfloat(i+i+1)
        dosl1(-i) = -re(i) - re(i) + uno
    enddo
    fact(0) = uno
    facts(-1) = raizpi
    facts(0) = facts(-1) * umed
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)   	!  i!
        facts(i) = facts(i-1) * re(i+i+1) * umed	! (i+1/2)!
    enddo
    allocate(alm((mxl+1)*(mxl+1)), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating alm in processor ',i3)") myrank
        abort = 1
        return
    endif
    do l = 0 , mxl
        do m = -l , l
            ma = abs(m)
            aux = dos * pi * fact(l+ma) / (fact(l-ma) * dosl1(l))
            if (m.eq.0) aux = aux + aux
            alm(l*(l+1)+m+1) = uno / sqrt(aux)
        end do
    end do
    return
    end
!
!   *****************************************************************
!
  subroutine solido(xa, ya, za)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL) :: ra2, xa, ya, za
    integer(KINT) :: l, m

    zlma(1) = uno		! Regular spherical harmonics of r-R(ia)
    if (lmaxbase .eq. 0) return
    zlma(2) = ya
    zlma(3) = za
    zlma(4) = xa
    ra2 = xa*xa + ya*ya + za*za
    do l = 1, lmaxbase-1
        zlma((l+1)*(l+3)+1) = dosl1(l) * (xa * zlma(l*(l+2)+1) - ya * zlma(l*l+1))		! zlm(l+1,l+1,ia)
        zlma((l+1)*(l+1)+1) = dosl1(l) * (ya * zlma(l*(l+2)+1) + xa * zlma(l*l+1))		! zlm(l+1,-(l+1),ia)
        zlma((l+2)*(l+2)-1) = dosl1(l) * za* zlma(l*(l+2)+1)				! zlm(l+1,l,ia)
        zlma(l*(l+2)+3) = dosl1(l) * za * zlma(l*l+1)					! zlm(l+1,-l,ia)
        do m = 0, l-1
            zlma((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)+m+1) - re(l+m)*ra2*zlma((l-1)*l+m+1))	! zlm(l+1,m,ia)
            zlma((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*za*zlma(l*(l+1)-m+1) - re(l+m)*ra2*zlma((l-1)*l-m+1))	! zlm(l+1,-m,ia)
        enddo
    enddo
    if (lgradient) call derivzlm(lmaxbase, (lmaxbase+1)*(lmaxbase+1), zlma, zlmadx, zlmady, zlmadz)
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
    if (lmax .eq. 3) return
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
!
!	***************************************************************
!
  subroutine leedatm2c
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, ierr, igrespec, indgr, iopsim, ios, ivclass, ivopsim, j, k, m2cmxcap, m2cmxfun, m2cmxcen
    integer(KINT) :: nbasis, nclassgr, numelem
    real(KREAL) :: repnuc, umbrznm2c, vchar, vchari
    character(6) :: grupo
    character(5) :: repirred
    logical lcmplxgr, lgrespec
    real(KREAL) :: rijk(3,3)
    logical :: falso
!	Reads the geometry and basis set (actually, it reads more data than required, because of the format of the *.sgbs file
!	generated by SMILES. These "extra" data are just skipped)
    open(15,file=trim(projectname)//".sgbs",form='unformatted', iostat=ierr)
    if (ierr .ne. 0) then
            write(6,"('Cannot open file ', a, '.sgbs in processor ',i3)") trim(projectname), myrank
            abort = 1
    endif
    rewind(15)
    read(15) umbrznm2c
    if (myrank .eq. 0 .and. longoutput)  write(6,"('umbrznm2c = ', e22.15)") umbrznm2c
    read(15) m2cmxcap, m2cmxfun, m2cmxcen
    if (myrank .eq. 0 .and. longoutput)  write(6,"('m2cmxcap, m2cmxfun, m2cmxcen = ', 3(1x,i4))")  m2cmxcap, m2cmxfun, m2cmxcen
    read(15) ncen
    if (myrank .eq. 0) write(6,"('ncen = ', i3)") ncen
    read(15) nbasis, ncaps
    read(15) repnuc
    if (myrank .eq. 0) write(6,"('nbasis, ncaps, repnuc = ', 2(1x,i4),3x,e22.15)") nbasis, ncaps, repnuc
    if (myrank .eq. 0) write(6,"('Number of basis functions = ', i4)") nbasis
    if (myrank .eq. 0) write(6,"('Number of function shells = ', i4)") ncaps

!	Allocates memory for geometry and basis set

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ll(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(lmaxc(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lmaxc in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(lsdisf(ncaps,ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lsdisf in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(nf(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating nf in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ngini(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ngini in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ngfin(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ngfin in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(nn(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating nn in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating nzn in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rcen in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(rnor(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rnor in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(xx(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating xx in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zn in processor ',i3)") myrank
        abort = 1
        return
    endif

    lmaxbase = 0
    nbas = nbasis
    i = 0
    do ia = 1, ncen
        lmaxc(ia) = 0
        read(15) ngini(ia), ngfin(ia)
        if (ngini(ia) .le. 0) cycle
        do k = ngini(ia), ngfin(ia)
            i = i + 1
            read(15) nf(i), nn(i), ll(i), xx(i), rnor(i)
            if (ll(i) .gt. lmaxbase) lmaxbase = ll(i)
            if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
        enddo
    enddo

    if (lmaxbase .gt. mxl) then
        if (myrank .eq. 0) write(6,"('Basis functions with not allowed values of  l. ')")
        if (myrank .eq. 0) write(6,"('Highest allowed value: ',i2 , ' Highest value in basis set: ', i2)") mxl, lmaxbase
        abort = 1
        return
    endif

    do i = 1, 3
        read(15) rijk(i,1), rijk(i,2), rijk(i,3)
    enddo
!	Reads the geometry and nuclear charges
    do ia = 1, ncen
        read(15) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia)
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
            nzn(ia) = 0
        else
            nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
    enddo
    close(15)

!	prints out the input data to standard output
    if (myrank .eq. 0) then
        write(6,"(/24x,'GEOMETRY (BOHR)')")
        write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge', t68, 'n. of shells')")
        do ia = 1, ncen
            if (ngini(ia) .gt. 0) then
                write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
                        ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia), ngfin(ia)-ngini(ia)+1
            else
                write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
                        ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia), 0
            endif
        enddo
        if (longoutput) then
            write(6,"(27x,'BASIS SET')")
            write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp',t60,'rnor')")
            do i = 1, ncaps
                write(6,"(t2, i3, t12, i2, t23, i3, t38, f12.7, t55, d17.10)") i, nn(i), ll(i), xx(i), rnor(i)
            enddo
        endif
            write(6,"('Number of basis functions = ', i4)") nbas
    endif
    return
    end
!
!    ***************************************************************
!
  subroutine leedatSTOgen(lsgbsden)
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, ib, ierr, indnf, iunit, j, k, nbasis
    integer(KINT), allocatable :: nshells(:)
    logical :: ldst, lsgbsden
    real(KREAL) :: rab
write(6,*) 'en leedatSTOgen, lsgbsden = ', lsgbsden
    if (lsgbsden) then
        iunit = 17
        open(iunit,file=trim(projectname)//".sgbsden",form='formatted', iostat=ierr)
        if (ierr .ne. 0) then
                iunit = 5
        endif
    else
        iunit = 5
    endif
!    Reads the number of centers
    read(iunit,*) ncen

!    Allocates memory for geometry and basis set
    ncaps = mxcap ! just for allocating

    allocate(atmnam(ncen), ll(ncaps), lmaxc(ncen), lsdisf(ncaps,ncaps), nf(ncaps), ngini(ncen), &
            ngfin(ncen), nn(ncaps), nzn(ncen), nshells(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam, ll, lmaxc, lsdisf, nf, ngini, ngfin, nn, nzn, nshells in processor ', &
                i3)") myrank
        abort = 1
        return
    endif

    allocate(rcen(3,ncen), rnor(ncaps), xx(ncaps), zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rcen, rnor, xx, zn in processor ',i3)") myrank
        abort = 1
        return
    endif

!    Reads geometry
    ncaps = 0
    do ia = 1, ncen
        read(iunit,*) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia), nshells(ia)
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
            nzn(ia) = 0
        else
            nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
        if (nshells(ia) .le. 0) cycle
        if (ncaps + nshells(ia) .gt. mxcap) then
            if (myrank .eq. 0) write(6,"('Error: maximum number of shells in basis set exceeded.')")
            abort = 1
            return
        endif
        ncaps = ncaps + nshells(ia)
    enddo
!    Reads basis set
    lmaxbase = 0
    nbas = 0
    ncaps = 0
    i = 0
    ngini(1) = 1
    indnf = 1
    do ia = 1, ncen
        if (nshells(ia) .gt. 0) then
            ncaps = ncaps + nshells(ia)
            ngfin(ia) = ngini(ia) + nshells(ia) - 1
            if (ia .lt. ncen) ngini(ia+1) = ngfin(ia) + 1
            lmaxc(ia) = 0
            do k = ngini(ia), ngfin(ia)
                i = i + 1
                read(iunit,*) nn(i), ll(i), xx(i)
                rnor(i) = sqrt( (dos*xx(i))**(2*nn(i)+1) / fact(2*nn(i)) )
                nf(i) = indnf
                indnf = indnf + 2*ll(i) + 1
                nbas = nbas + 2*ll(i) + 1
                if (ll(i) .gt. lmaxbase) lmaxbase = ll(i)
                if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
            enddo
        else
            ngini(ia) = -1
            ngfin(ia) = -1
        endif
    enddo
    if (lmaxbase .gt. mxl) then
        if (myrank .eq. 0) write(6,"('Basis functions with not allowed values of  l. ')")
        if (myrank .eq. 0) write(6,"('Highest allowed value: ', i2, ' Highest value in basis set: ', i2)") mxl, lmaxbase
        abort = 1
        return
    endif
    nbasis = nbas
    if (myrank .eq. 0) then
!    prints out the input data
        write(6,"(/24x,'GEOMETRY (BOHR)')")
        write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge', t68, 'n. of shells')")
        do i = 1, ncen
            write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
                    i, rcen(1,i), rcen(2,i), rcen(3,i) , zn(i), ngfin(i)-ngini(i)+1
        enddo
        if (longoutput) then
            write(6,"(27x,'BASIS SET')")
            write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp',t60,'rnor')")
            do i = 1, ncaps
                write(6,"(t2, i3, t12, i2, t23, i3, t38, f12.7, t55, d17.10)") i, nn(i), ll(i), xx(i), rnor(i)
            enddo
        endif
        write(6,"('Number of basis functions = ', i4)") nbas
    endif
    return
    end
!
!	***************************************************************
!
  subroutine leedatgauss
    USE MPI
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ia, icarga, ierr, indnf, indng, ios, j, k, k1, k2, knt, nbasis, ncapserr
    real(KREAL) :: aux, bux
    real(KREAL) :: xaux(mxprimit), cfaux(mxprimit)
!	Reads the number of centers
    open(15,file=trim(projectname)//".ggbs",form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, '.ggbs in processor ',i3)") trim(projectname), myrank
        abort = 1
    endif
    read(15,*) ncen

!	Allocates memory for geometry and basis set
    ncaps = mxcap ! just for allocating

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating atmnam in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(cfcontr0(mxcap*mxprimit), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating cfcontr0 in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ipntprim(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ipntprim in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(isort(mxprimit), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating isort in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ll(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ll in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(lmaxc(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating lmaxc in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ncontr(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ncontr in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(nf(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating nf in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ngini(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ngini in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(ngfin(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating ngfin in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(nprimit(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating nprimit in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating nzn in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rcen in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(rnor(ncaps), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating rnor in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(xxg0(mxcap*mxprimit), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating xxg0 in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Memory error when allocating zn in processor ',i3)") myrank
        abort = 1
        return
    endif

!	Reads the number of centers, geometry and basis set 
    do ia = 1, ncen
        read(15,*) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia)
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
            nzn(ia) = 0
        else
            nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
    enddo

!     Basis set::
!         for each center:
!              number of contractions
!              for each contraction:
!                   number of primitives
!                   l  quantum number
!                   exponents of the primitives
!                   contraction coefficients
!     Primitives are sorted in increasing order of exponents.
!     An array is created with pointers to the beginning of
!     exponents and coefficients of a contracted function.
!     Arrays  nf, ngini  and  ngfin are initialized.
!
!	Contraction coefficients correspond to the expansion in UNNORMALIZED primitives	

    lmaxbase = 0
    ncaps = 0
    icarga = 0
    indnf = 1
    indng = 1
    ncaps = 0
    ncapserr = 0
    lmaxbase = 0
    do ia = 1, ncen
        read(15,*) ncontr(ia)
        if (ncontr(ia) .le. 0) then
            ngini(ia) = -1
            ngfin(ia) = -1
            cycle
        endif
        ngini(ia) = indng
        ngfin(ia) = indng + ncontr(ia) - 1
        indng = indng + ncontr(ia)
        lmaxc(ia) = 0
        do j = 1, ncontr(ia)
            ncaps = ncaps + 1
            if (ncaps .gt. mxcap) then
                write(6,"('Error: maximum number of shells in basis set exceeded in processor ',i3)") myrank
                abort = 1
                return
            endif
            read(15,*) nprimit(ncaps), ll(ncaps)
            if (icarga+nprimit(ncaps) .gt. mxcap*mxprimit)  then
                write(6,"('Error: maximum number of total primitives exceeded. &
                        &\nChange parameter mxcap in DAMGLOBALxxxx.F90 and remake')")
                abort = 1
                return
            endif
            if (ll(ncaps) .gt. lmaxbase) lmaxbase = ll(ncaps)
            if (ll(ncaps) .gt. lmaxc(ia)) lmaxc(ia) = ll(ncaps)
            nf(ncaps) = indnf
            indnf = indnf + 2*ll(ncaps) + 1
            if(nprimit(ncaps) .gt. mxprimit) then
                write(6,"('Error: maximum number of primitives per contraction exceeded in processor ',i3)") myrank
                abort = 1
                return
            endif
            ipntprim(ncaps) = icarga+1
            read(15,*) (xaux(k), k = 1, nprimit(ncaps))
            read(15,*) (cfaux(k), k = 1, nprimit(ncaps))
!			sorts primitives
            call sort(nprimit(ncaps),xaux)
            if (abort .gt. 0) return
            do k = 1, nprimit(ncaps)
                xxg0(icarga+k) = xaux(k)
                cfcontr0(icarga+k) = cfaux(isort(k))
            enddo
! 			kntprim = kntprim + nprimit(ncaps)
!			computes and stores the radial normalization factor
            aux = cero
            bux = ll(ncaps) + 1.5d0
            do k1 = 1, nprimit(ncaps)
                do k2 = 1, k1-1
                    aux=aux + dos*cfcontr0(icarga+k1)*cfcontr0(icarga+k2)/(xaux(k1)+xaux(k2))**bux
                enddo
                aux = aux + cfcontr0(icarga+k1) * cfcontr0(icarga+k1) / (dos*xaux(k1))**bux
            enddo
            rnor(ncaps) = sqrt( dos / (facts(ll(ncaps))*aux) )
!			actualizes the index for loading
            icarga = icarga+nprimit(ncaps)
        enddo
    enddo

    nprimitot = icarga
    nbas = indnf-1

    allocate(xxg(nprimitot), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error when allocating xxg in processor ',i3)") myrank
        abort = 1
        return
    endif

    allocate(cfcontr(nprimitot), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error when allocating cfcontr in processor ',i3)") myrank
        abort = 1
        return
    endif

    xxg(1:nprimitot) = xxg0(1:nprimitot)
    cfcontr(1:nprimitot) = cfcontr0(1:nprimitot)
    deallocate(xxg0, cfcontr0)

    if (ncaps .gt. mxcap .and. myrank .eq. 0) then
        write(6,"('Number of function shells = ', i4, ' higher  than maximum allowed = ',i4)")  ncaps, mxcap
        write(6,"('Modify parameter  mxcap  in module DAM320_D of file DAM320_GLOBAL.F90 and recompile.')")
        abort = 1
        return
    endif
    if (lmaxbase .gt. mxl .and. myrank .eq. 0) then
        write(6,"('Basis functions with not allowed values of  l. ')")
        write(6,"('Highest allowed value: ', i2 , ' Highest value in basis set: ', i2)") mxl, lmaxbase
        abort = 1
        return
    endif

!	Allocates the array containing the density matrix

    allocate(dmat(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) then
        write(6,"('Error when allocating dmat in processor ',i3)") myrank
        abort = 1
        return
    endif

    if (myrank .eq. 0 .and. longoutput) write(6,"('Estimated highest size of dmat   = ', i15, ' bytes')") size(dmat)

!	Reads the density matrix in lower triangle form:  read ((dmat(i,j), j = 1, i), i = 1, nbasis)
    open(16,file=trim(projectname)//".den",form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        write(6,"('Cannot open file ', a, '.den in processor ',i3)") trim(projectname), myrank
        abort = 1
    endif
    read(16,*, iostat = ios) nbasis, ((dmat(i,j),j=1,i),i=1,nbasis)
    if ( ios .ne. 0 .or. nbas .ne. nbasis  .and. myrank .eq. 0) then
        write(6,"('nbas = ', i5,' nbasis = ', i5)") nbas, nbasis
        write(6,"('ERROR reading density matrix. Check whether the density matrix correspond to this basis set.')")
        abort = 1
        return
    endif

    close(16)
    do i = 2, nbasis
        do j = 1, i-1
            dmat(j,i) = dmat(i,j)
        enddo
    enddo

!	prints out the input data
    if (myrank .eq. 0) then
        write(6,"(/24x,'GEOMETRY (BOHR)')")
        write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge', t68, 'n. of shells')")
        do ia = 1, ncen
            write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
                    ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia), ngfin(ia)-ngini(ia)+1
        enddo
        write(6,"(27x,'GTO BASIS SET')")
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
        write(6,"('Number of primitive functions = ', i4)") nprimitot
    endif
    return
    end
	
!   ***************************************************************
!
!	Subroutine sort: sorts the primitives of a contraction in ascending exponents 

  subroutine sort(nprim, xsort)
    USE MPI
    USE GAUSS
    USE PARALELO
    implicit none
    integer(KINT) :: i, ierr, iux, j, nprim
    real(KREAL) :: aux
    integer(KINT), parameter :: mxsteps = 1000
    real(KREAL) :: xsort(*)
    logical :: lend
    do i = 1, nprim
        isort(i) = i
    enddo
    do i = 1, mxsteps
        lend = .true.
        do j = nprim, 2, -2
            if (xsort(j) .lt. xsort(j-1)) then
                aux = xsort(j)
                xsort(j) = xsort(j-1)
                xsort(j-1) = aux
                iux = isort(j)
                isort(j) = isort(j-1)
                isort(j-1) = iux
                lend = .false.
            endif
        enddo
        do j = nprim-1, 2, -2
            if (xsort(j) .lt. xsort(j-1)) then
                aux = xsort(j)
                xsort(j) = xsort(j-1)
                xsort(j-1) = aux
                iux = isort(j)
                isort(j) = isort(j-1)
                isort(j-1) = iux
                lend = .false.
            endif
        enddo
        if (lend) return
    enddo
    write(6,"('Error in subroutine sort. Highest number of steps (',i4,') exceeded in processor ',i3)") mxsteps, myrank
    abort = 1
    return
    end

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
    USE PARALELO
    implicit none
    integer(KINT) :: ierr, ierr2
    character(*) :: msg
    if (myrank .eq. 0) then
        write(6,"(a)") msg
        write(6,"('Error code = ', i4)") ierr
    endif
    call MPI_FINALIZE(ierr2)
    stop
    end
