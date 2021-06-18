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
!	Program for generating grids for plotting orbitals
!
! Version of September 2018
!
  program DAMORB320
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL4) :: tarray(2), tiempo, dtime
    real(KREAL) :: aux, x, xmax, xmin, xyzmax, xyzmin, y, ymax, ymin, z, zmax, zmin
    integer(KINT) :: i, ierr, knt, norbs
    logical :: existe
    namelist / options / dltu, dltv, dltx, dlty, dltz, filename, fileMOname, iorbsinp, iswindows, langstrom, lgradient,  &
            lgrid2d, lm2c, norbs, planeA, planeB, planeC, planecase, &
            uinf, usup, vinf, vsup, x_func_uv, xinf, xsup, y_func_uv, yinf, ysup, z_func_uv, zinf, zsup
    tiempo = dtime(tarray)
!	Namelist default values
    iswindows = .false.            ! .true. if running on a MS-windows system
    langstrom = .true.		! If .false. distances in bohr
    lm2c = .false.          ! If .true. read data from a calculation with m2c
    lgradient = .false.
    lgrid2d = .false.
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
    x_func_uv = 'u'		! x = u for 2D grids:  default: plane XY
    y_func_uv = 'v'		! y = v for 2D grids
    z_func_uv = '0'		! z = 0 for 2D grids

    planeA = cero        ! Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    planeB = cero
    planeC = uno
    planecase = 1

    filename = ""			! root file name for .plt and .pltd files
    fileMOname = ""		! file with Molecular orbitals coefficients
!	End of namelist defaults

    planesuffix = ""
    read(5,OPTIONS)
    read(5,*) projectname
    call subplanesuffix(planecase,planesuffix)
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
    if (len_trim(fileMOname) .eq. 0) then
        inquire(file=projectname(1:i)//".orb", exist=existe)
        if (existe) then
            fileMOname = projectname//".orb"
        else
            inquire(file=projectname//".GAorba", exist=existe)
            if (existe) then
                fileMOname = projectname//".GAorba"
            else
                call error(2, 'Cannot find a file with MO named '//projectname//'+ .orb or GAorba')
            endif
        endif
    else
        i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
        if (i .eq. 0) i = index(projectname,"/",.true.)	! In case of windows with MinGW directory separator is "/"
        fileMOname = projectname(1:i)//trim(fileMOname)
        inquire(file=fileMOname, exist=existe)
        if (.not. existe) then
            call error(2, 'Cannot find a file with MO named '//fileMOname)
        endif
    endif
    if (lgrid2d) lgradient = .false.

    call consta
    allocate(iorbs(norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating iorbs. Stop')

    if (fileMOname(len_trim(fileMOname)-3:len_trim(fileMOname)) == ".orb" .or. &
        fileMOname(len_trim(fileMOname)-6:len_trim(fileMOname)-1) == ".SLorb") then
        lsto = .true.
    elseif (fileMOname(len_trim(fileMOname)-6:len_trim(fileMOname)-1) == ".GAorb") then
        lsto = .false.
    else
        i = index(fileMOname,".",.true.)	! Checks position of last directory name separator
        write(6,"('MO filename extension:  ', a, '  not acceptable. Must be .orb, .SLorba, .SLorbb, .GAorba or .GAorbb ')") &
                trim(fileMOname(i:len_trim(fileMOname)))
        call error(2,'Stop')
    endif
    if (lsto) then
        if (lm2c) then
            call leedatm2c        ! Reads input from SMILES files
        else        ! Reads input from general input file
            call leedatSTOgen
        endif
    else
        call leedatgauss
    endif

    knt = 0
    do i = 1, norbs
        if (iorbsinp(i) .le. nbas) then
            knt = knt + 1
            iorbs(knt) = iorbsinp(i)
        endif
    enddo
    norbs = knt

!		Checks that the orbital indices are in the range 1:nbas

    allocate(corbs(nbas,norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating corbs. Stop')

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
        if (.not. lgrid2d) then
            call pltorbSMILES(norbs)
        else
            call pltorbSMILES2d(norbs)
        endif
    else
        call leeorbGTO(norbs)
        if (.not. lgrid2d) then
            call pltorbGTO(norbs)
        else
            call pltorbGTO2d(norbs)
        endif
    endif
    tiempo = dtime(tarray)
    write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
            tarray(1), tarray(2), tarray(1)+tarray(2)
    stop
    end
!
!   *************************************************************
!
!  	Reads the coefficients of Slater MO from a calculation with m2c and returns them in c
!
subroutine leeorbSMILES(norbs)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL), allocatable :: caux(:,:)
    integer(KINT) :: i, ierr, j, nbasis, norbs
    allocate(caux(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error in leeorbSMILES when allocating caux. Stop')
    open(13,file=trim(fileMOname),form='unformatted', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Cannot open file'//trim(fileMOname)//'. Stop')

    read(13) nbasis, ((caux(i,j),i=1,nbas),j=1,nbas)
    close(13)
    if (nbas.ne.nbasis) then
        call error(1,'Wrong number of basis functions. Stop')
    endif
    do i = 1, norbs
        corbs(1:nbas,i) = caux(1:nbas,iorbs(i))
    enddo
    deallocate(caux)
    return
    end
!
!   *************************************************************
!
!  	Reads the coefficients of Slater MO  and returns them in c
!
subroutine leeorbgen(norbs)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL), allocatable :: caux(:,:)
    integer(KINT) :: i, ierr, j, nbasis, nmos, norbs, nvoid
    open(13,file=trim(fileMOname),form='formatted', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Cannot open file'//trim(fileMOname)//'. Stop')
    read(13,*) nbasis, nvoid, nmos
    allocate(caux(nbasis,nmos), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error in leeorbgen when allocating caux. Stop')
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
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL), allocatable :: caux(:,:)
    integer(KINT) :: i, ierr, j, nbasis, numorb, norbs, nvoid
    allocate(caux(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating caux. Stop')
    open(13,file=trim(fileMOname),form='formatted', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Cannot open file'//trim(fileMOname)//'. Stop')

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
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL4), allocatable :: array(:,:), arraydx(:,:), arraydy(:,:), arraydz(:,:)
    real(KREAL) :: b2a, x, y, z
    integer(KINT) :: iaux, ierr, iorb, iuni, ix, iy, iz, norbs, nx, ny, nz
    character(2) faux
    real(KREAL4) :: x1a, x2a, y1a, y2a, z1a, z2a

    allocate(orbv(norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating orbv. Stop')
    allocate(chi(nbas), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating chi. Stop')
    allocate(zlma((lmaxbase+1)*(lmaxbase+1)), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rz. Stop')

    b2a = 0.5291772d0
    nx = (xsup-xinf) / dltx + 1
    ny = (ysup-yinf) / dlty + 1
    nz = (zsup-zinf) / dltz + 1
    allocate(array(nx,norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating array. Stop')
    if (lgradient) then
        allocate(orbvdx(nbas), orbvdy(nbas), orbvdz(nbas), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating orbvdx, orbvdy, orbvdz. Stop')
        allocate(chidx(nbas), chidy(nbas), chidz(nbas), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating chidx, chidy, chidz. Stop')
        allocate(arraydx(nx,norbs), arraydy(nx,norbs), arraydz(nx,norbs), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating arraydx, arraydy, arraydz. Stop')
        allocate(zlmadx((lmaxbase+1)*(lmaxbase+1)), zlmady((lmaxbase+1)*(lmaxbase+1)), zlmadz((lmaxbase+1)*(lmaxbase+1)) &
                , stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating zlmadx, zlmady, zlmadz. Stop')
    endif
    x1a = xinf * b2a
    x2a = xsup * b2a
    y1a = yinf * b2a
    y2a = ysup * b2a
    z1a = zinf * b2a
    z2a = zsup * b2a
!
!  Headers for orbital plt files
!
    iuni = 10
    do iorb = 1, norbs
        iuni = iuni + 1
        write(faux,'(i2.2)') iorbs(iorb)
        call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//".plt")
        if (lgradient) then
            iuni = iuni + 1
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dx.pltd")
            iuni = iuni + 1
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dy.pltd")
            iuni = iuni + 1
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dz.pltd")
        endif
    end do
!
!  Generates grids
!
    do iz = 1 , nz
        z = zinf + (iz-1)*dltz
        do iy = 1 , ny
            y = yinf + (iy-1)*dlty
            do ix = 1 , nx
                x = xinf + (ix-1)*dltx
                call orbitalesSMILES( norbs, x , y , z )
                do iorb = 1, norbs
                    array(ix,iorb) = orbv(iorb)
                end do
                if (lgradient) then
                    do iorb = 1, norbs
                        arraydx(ix,iorb) = orbvdx(iorb)
                        arraydy(ix,iorb) = orbvdy(iorb)
                        arraydz(ix,iorb) = orbvdz(iorb)
                    end do
                endif
            end do
            iuni = 10
            do iorb = 1, norbs
                iuni = iuni + 1
                write(iuni) array(1:nx,iorb)
                if (lgradient) then
                    write(iuni+1) arraydx(1:nx,iorb)
                    write(iuni+2) arraydy(1:nx,iorb)
                    write(iuni+3) arraydz(1:nx,iorb)
                    iuni = iuni+3
                endif
            end do
        end do
    end do
!
!  Close files
!
    iuni = 10
    iaux = norbs
    if (lgradient) iaux = 4*iaux
    do iorb = 1 , iaux
        close(iuni+iorb)
    end do
    return
    end
!                                                                               
! *******************************************************************
!                                                                               
  subroutine pltorbSMILES2d(norbs)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    use strings
    use evaluate
    implicit none
    real(KREAL4), allocatable :: array(:,:), arraydx(:,:), arraydy(:,:), arraydz(:,:)
    real(KREAL) :: b2a, ru, rv, u, v, x, y, z
    integer(KINT) :: iaux, iorb, iuni, iu, iv, norbs, nu, nv
    character(2) faux
    real(KREAL4) :: u1, u2, v1, v2
    allocate(orbv(norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating orbv. Stop')
    allocate(chi(nbas), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating chi. Stop')
    allocate(zlma((lmaxbase+1)*(lmaxbase+1)), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rz. Stop')

    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1
    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1

    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv
    if (langstrom) then		! Converts grid distances to angstrom
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

!	Opens file for electrostatic potential tabulation
    allocate(array(nu,norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating array in pltorbSMILES2d. Stop')
!
!  Headers for orbital plt files
!
    iuni = 10
    do iorb = 1, norbs
        iuni = iuni + 1
        write(faux,'(i2.2)') iorbs(iorb)
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(planesuffix)//"_MO_"//faux//".cnt")
    end do
!
!  Generates grids
!
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)	! Computes x(u,v)
            call evalexpr(y_func_uv,y)	! Computes y(u,v)
            call evalexpr(z_func_uv,z)	! Computes z(u,v)
            call orbitalesSMILES( norbs, x , y , z )
            do iorb = 1, norbs
                array(iu,iorb) = orbv(iorb)
            end do
            iuni = 10
            do iorb = 1, norbs
                iuni = iuni + 1
                write(iuni) array(1:nu,iorb)
            end do
        end do
    end do
!
!  Close files
!
    iuni = 10
    iaux = norbs
    do iorb = 1 , iaux
        close(iuni+iorb)
    end do
    return
    end
!                                                                               
! *******************************************************************
!                                                                               
  subroutine pltorbGTO(norbs)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    implicit none
    real(KREAL4), allocatable :: array(:,:), arraydx(:,:), arraydy(:,:), arraydz(:,:)
    real(KREAL) :: b2a, x, y, z
    integer(KINT) :: iaux, ierr, iorb, iuni, ix, iy, iz, norbs, nx, ny, nz
    character(2) faux
    real(KREAL4) :: x1a, x2a, y1a, y2a, z1a, z2a
    allocate(orbv(norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating orbv. Stop')
    allocate(chi(nbas), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating chi. Stop')
    allocate(zlma((lmaxbase+1)*(lmaxbase+1)), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zlma. Stop')

    b2a = 0.5291772d0
    nx = (xsup-xinf) / dltx + 1
    ny = (ysup-yinf) / dlty + 1
    nz = (zsup-zinf) / dltz + 1
    allocate(array(nx,norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating array. Stop')
    if (lgradient) then
            allocate(orbvdx(nbas), orbvdy(nbas), orbvdz(nbas), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating orbvdx, orbvdy, orbvdz. Stop')
            allocate(chidx(nbas), chidy(nbas), chidz(nbas), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating chidx, chidy, chidz. Stop')
            allocate(arraydx(nx,norbs), arraydy(nx,norbs), arraydz(nx,norbs), stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating arraydx, arraydy, arraydz. Stop')
            allocate(zlmadx((lmaxbase+1)*(lmaxbase+1)), zlmady((lmaxbase+1)*(lmaxbase+1)), zlmadz((lmaxbase+1)*(lmaxbase+1)) &
                    , stat = ierr)
            if (ierr .ne. 0) call error(1,'Memory error when allocating zlmadx, zlmady, zlmadz. Stop')
    endif
    x1a = xinf * b2a
    x2a = xsup * b2a
    y1a = yinf * b2a
    y2a = ysup * b2a
    z1a = zinf * b2a
    z2a = zsup * b2a
!
!  Headers for orbital plt files
!
    iuni = 10
    do iorb = 1, norbs
        iuni = iuni + 1
        write(faux,'(i2.2)') iorbs(iorb)
        call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//".plt")
        if (lgradient) then
            iuni = iuni + 1
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dx.pltd")
            iuni = iuni + 1
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dy.pltd")
            iuni = iuni + 1
            call cabecera(nx, ny, nz, x1a, x2a, y1a, y2a, z1a, z2a, iuni, trim(filename)//"_MO_"//faux//"-dz.pltd")
        endif
    end do
!
!  Generates grids
!
    do iz = 1 , nz
        z = zinf + (iz-1)*dltz
        do iy = 1 , ny
            y = yinf + (iy-1)*dlty
            do ix = 1 , nx
                x = xinf + (ix-1)*dltx
                call orbitalesGTO( norbs, x , y , z )
                do iorb = 1, norbs
                        array(ix,iorb) = orbv(iorb)
                end do
                if (lgradient) then
                    do iorb = 1, norbs
                        arraydx(ix,iorb) = orbvdx(iorb)
                        arraydy(ix,iorb) = orbvdy(iorb)
                        arraydz(ix,iorb) = orbvdz(iorb)
                    end do
                endif
            end do
            iuni = 10
            do iorb = 1, norbs
                iuni = iuni + 1
                write(iuni) array(1:nx,iorb)
                if (lgradient) then
                    write(iuni+1) arraydx(1:nx,iorb)
                    write(iuni+2) arraydy(1:nx,iorb)
                    write(iuni+3) arraydz(1:nx,iorb)
                    iuni = iuni+3
                endif
            end do
        end do
    end do
!
!  Close files
!
    iuni = 10
    iaux = norbs
    if (lgradient) iaux = 4*iaux
    do iorb = 1 , iaux
        close(iuni+iorb)
    end do
    return
    end
!                                                                               
! *******************************************************************
!                                                                               
  subroutine pltorbGTO2d(norbs)
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMORB320_D
    use strings
    use evaluate
    implicit none
    real(KREAL4), allocatable :: array(:,:)
    real(KREAL) :: b2a, ru, rv, u, v, x, y, z
    integer(KINT) :: iaux, iorb, iuni, iu, iv, norbs, nu, nv
    character(2) faux
    real(KREAL4) :: u1, u2, v1, v2
    allocate(orbv(norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating orbv. Stop')
    allocate(chi(nbas), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating chi. Stop')
    allocate(zlma((lmaxbase+1)*(lmaxbase+1)), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rz. Stop')

    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1
    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1

    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv
    if (langstrom) then		! Converts grid distances to angstrom
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

!	Opens file for electrostatic potential tabulation
    allocate(array(nu,norbs), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating array in pltorbSMILES2d. Stop')
!
!  Headers for orbital plt files
!
    iuni = 10
    do iorb = 1, norbs
        iuni = iuni + 1
        write(faux,'(i2.2)') iorbs(iorb)
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(planesuffix)//"_MO_"//faux//".cnt")
    end do
!
!  Generates grids
!
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)	! Computes x(u,v)
            call evalexpr(y_func_uv,y)	! Computes y(u,v)
            call evalexpr(z_func_uv,z)	! Computes z(u,v)
            call orbitalesGTO( norbs, x , y , z )
            do iorb = 1, norbs
                array(iu,iorb) = orbv(iorb)
            end do
        end do
        iuni = 10
        do iorb = 1, norbs
            iuni = iuni + 1
            write(iuni) array(1:nu,iorb)
        end do
    end do

!
!  Close files
!
    iuni = 10
    iaux = norbs
    do iorb = 1 , iaux
        close(iuni+iorb)
    end do
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
!	***************************************************************
!	Subroutine cabecera: writes head for .plt files (binary)
!
  subroutine cabecera2d(nu, nv, u1, u2, v1, v2, iuni, s)
    USE DAM320_D
    implicit none
    integer(KINT) :: i
    integer(KINT) :: nu,nv,iuni,ns,iaux(0:1)
    real(KREAL4) :: u1,u2,v1,v2,v(0:3)
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
!   ***************************************************************
!
   subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
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
    if (ierr .ne. 0) call error(1,'Memory error when allocating alm. Stop')
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
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, igrespec, ierr, indgr, iopsim, ios, ivclass, ivopsim, j, k, m2cmxcap, m2cmxfun, m2cmxcen
    integer(KINT) :: nbasis, nclassgr, numelem
    real(KREAL) :: repnuc, umbrznm2c, vchar, vchari
    character(6) grupo
    character(5) repirred
    logical lcmplxgr, lgrespec
    real(KREAL) :: rijk(3,3)
    logical :: falso
!	Reads the geometry and basis set (actually, it reads more data than required, because of the format of the *.sgbs file
!	generated by SMILES. These "extra" data are just skipped)
    open(15,file=trim(projectname)//".sgbs",form='unformatted', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Cannot open file'//trim(projectname)//".sgbs"//'. Stop')
    rewind(15)
    read(15) umbrznm2c
    if (longoutput) write(6,"('umbrznm2c = ', e22.15)") umbrznm2c
    read(15) m2cmxcap, m2cmxfun, m2cmxcen
    if (longoutput) write(6,"('m2cmxcap, m2cmxfun, m2cmxcen = ', 3(1x,i4))")  m2cmxcap, m2cmxfun, m2cmxcen
    read(15) ncen
    write(6,"('ncen = ', i3)") ncen
    read(15) nbas, ncaps
    read(15) repnuc
    write(6,"('Nuclear repulsion = ', e22.15)") repnuc
    write(6,"('Number of basis functions = ', i4)") nbas
    write(6,"('Number of function shells = ', i4)") ncaps

!	Allocates memory for geometry and basis set

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating atmnam. Stop')

    allocate(ll(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ll. Stop')

    allocate(lmaxc(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating lmaxc. Stop')

    allocate(lsdisf(ncaps,ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating lsdisf. Stop')

    allocate(nf(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nf. Stop')

    allocate(ngini(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')

    allocate(ngfin(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ngfin. Stop')

    allocate(nn(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nn. Stop')

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nzn. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rcen. Stop')

    allocate(rnor(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rnor. Stop')

    allocate(xx(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating xx. Stop')

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zn. Stop')

!	Reads the basis set
    lmaxbase = 0
    i = 0
    do ia = 1, ncen
        read(15) ngini(ia), ngfin(ia)
        lmaxc(ia) = 0
        if (ngini(ia) .le. 0) cycle
        do k = ngini(ia), ngfin(ia)
            i = i + 1
            read(15) nf(i), nn(i), ll(i), xx(i), rnor(i)
            if (ll(i) .gt. lmaxbase) lmaxbase = ll(i)
            if (ll(i) .gt. lmaxc(ia)) lmaxc(ia) = ll(i)
        enddo
    enddo
    if (lmaxbase .gt. mxl) then
        write(6,"('Basis functions with not allowed values of  l. ')")
        write(6,"('Highest allowed value: ', i2 , ' Highest value in basis set: ', i2)") mxl, lmaxbase
        call error(1,' Stop')
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
        write(6,"(27x,'STO BASIS SET')")
        write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp',t60,'rnor')")
        do i = 1, ncaps
            write(6,"(t2, i3, t12, i2, t23, i3, t38, f12.7, t55, d17.10)") i, nn(i), ll(i), xx(i), rnor(i)
        enddo
    endif
    write(6,"('Number of basis functions = ', i4)") nbas
    return
    end
!
!    ***************************************************************
!
  subroutine leedatSTOgen
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: i, ia, ib, ierr, indnf, iunit, j, k, nbasis, ngi, ngf
    integer(KINT), allocatable :: nshells(:)
    logical :: ldst, lsgbs, lsgbsden, lsgbsdengz
    real(KREAL) :: rab

!     Checks if file .sgbsden or .sgbsden.gz exist. If yes, read geometry, basis set and density from it. Otherwise,
!     read these data from standard input.
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

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating atmnam. Stop')

    allocate(ll(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ll. Stop')

    allocate(lmaxc(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating lmaxc. Stop')

    allocate(lsdisf(ncaps,ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating lsdisf. Stop')

    allocate(nf(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nf. Stop')

    allocate(ngini(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')

    allocate(ngfin(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ngfin. Stop')

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nzn. Stop')

    allocate(nshells(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nshells. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rcen. Stop')

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zn. Stop')

    allocate(xx(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating xx. Stop')

    allocate(rnor(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rnor. Stop')

    allocate(nn(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nn. Stop')

!    Reads geometry, nuclear charge and number of function shells per center

    ncaps = 0
    do ia = 1, ncen
        read(iunit,*) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia), nshells(ia)
        if (nshells(ia) .le. 0) cycle
        if (ncaps + nshells(ia) .gt. mxcap) call error(1,'Error: maximum number of shells in basis set exceeded. Stop')
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
            nzn(ia) = 0
        else
            nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
        ncaps = ncaps + nshells(ia)
    enddo

!    Reads the basis set

    lmaxbase = 0
    nbas = 0
    i = 0
    ngi = 1
    ngf = 0
    indnf = 1
    do ia = 1, ncen
        if (nshells(ia) .gt. 0) then
            ngini(ia) = ngi
            ngf = ngi + nshells(ia) - 1
            ngfin(ia) = ngf
            ngi = ngf + 1
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
        write(6,"('Basis functions with not allowed values of  l. ')")
        write(6,"('Highest allowed value: ', i2 , ' Highest value in basis set: ', i2)") mxl, lmaxbase
        call error(1,' Stop')
    endif
    nbasis = nbas

!    prints out the input data to standard output
    write(6,"(/24x,'GEOMETRY (BOHR)')")
    write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge', t68, 'n. of shells')")
    do i = 1, ncen
        if (ngini(i) .gt. 0) then
            write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
                    i, rcen(1,i), rcen(2,i), rcen(3,i) , zn(i), ngfin(i)-ngini(i)+1
        else
            write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
                    i, rcen(1,i), rcen(2,i), rcen(3,i) , zn(i), 0
        endif
    enddo
    if (longoutput) then
        write(6,"(27x,'STO BASIS SET')")
        write(6,"(/t1,' shell:',t13,'n',t25,'l',t43,'exp',t60,'rnor')")
        do i = 1, ncaps
            write(6,"(t2, i3, t12, i2, t23, i3, t38, f12.7, t55, d17.10)") i, nn(i), ll(i), xx(i), rnor(i)
        enddo
    endif
    write(6,"('Number of basis functions = ', i4)") nbas
    deallocate(nshells)
    return
    end
!
!	***************************************************************
!
  subroutine leedatgauss
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, ierr, indnf, indng, ios, j, k, k1, k2, knt, nbasis, ncapserr
    real(KREAL) :: aux, bux
    real(KREAL) :: xaux(mxprimit), cfaux(mxprimit)
!	Reads the number of centers
    open(15,file=trim(projectname)//".ggbs",form='formatted', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Cannot open file'//trim(projectname)//".ggbs"//'. Stop')
    read(15,*) ncen

!	Allocates memory for geometry and basis set
    ncaps = mxcap ! just for allocating

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating atmnam. Stop')

    allocate(cfcontr0(mxcap*mxprimit))
    if (.not. allocated(cfcontr0)) call error(1,'Memory error when allocating cfcontr0. Stop')

    allocate(ipntprim(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ipntprim. Stop')

    allocate(isort(mxprimit), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating isort. Stop')

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

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nzn. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rcen. Stop')

    allocate(rnor(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rnor. Stop')

    allocate(xxg0(mxcap*mxprimit))
    if (.not. allocated(xxg0)) call error(1,'Memory error when allocating xxg0. Stop')

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zn. Stop')

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
            if (ncaps .gt. mxcap)  call error(1,'Error: maximum number of shells in basis set exceeded. Stop')
            read(15,*) nprimit(ncaps), ll(ncaps)
            if (icarga+nprimit(ncaps) .gt. mxcap*mxprimit)  call error(1,'Error: maximum number of total primitives exceeded. &
                            &\nChange parameter mxcap in DAMGLOBALxxxx.F90 and remake. \nStop')
            if (ll(ncaps) .gt. lmaxbase) lmaxbase = ll(ncaps)
            if (ll(ncaps) .gt. lmaxc(ia)) lmaxc(ia) = ll(ncaps)
            nf(ncaps) = indnf
            indnf = indnf + 2*ll(ncaps) + 1
            if(nprimit(ncaps) .gt. mxprimit) call error(1,'Error: maximum number of primitives per contraction exceeded. Stop')
            ipntprim(ncaps) = icarga+1
            read(15,*) (xaux(k), k = 1, nprimit(ncaps))
            read(15,*) (cfaux(k), k = 1, nprimit(ncaps))
!			sorts primitives in increasing exponents
            call sort(nprimit(ncaps),xaux)
            do k = 1, nprimit(ncaps)
                xxg0(icarga+k) = xaux(k)
                cfcontr0(icarga+k) = cfaux(isort(k))
            enddo
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

    allocate(xxg(nprimitot))
    if (.not. allocated(xxg)) call error(1,'Memory error when allocating xxg. Stop')

    allocate(cfcontr(nprimitot))
    if (.not. allocated(cfcontr)) call error(1,'Memory error when allocating cfcontr. Stop')

    xxg(1:nprimitot) = xxg0(1:nprimitot)
    cfcontr(1:nprimitot) = cfcontr0(1:nprimitot)

    deallocate(xxg0, cfcontr0)

    if (ncaps .gt. mxcap) then
        write(6,"('Number of function shells = ', i4, ' higher  than maximum allowed = ',i4)")  ncaps, mxcap
        write(6,"('Modify parameter  mxcap  in module DAM320_D of file DAM320_GLOBAL.F90 and recompile.')")
        call error(1,' Stop')
    endif
    if (lmaxbase .gt. mxl) then
        write(6,"('Basis functions with not allowed values of  l. ')")
        write(6,"('Highest allowed value: ', i2 , ' Highest value in basis set: ', i2)") mxl, lmaxbase
        call error(1,' Stop')
    endif

!	prints out the input data
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
    return
    end
	
!   ***************************************************************
!
!	Subroutine sort: sorts the primitives of a contraction in ascending exponents 

   subroutine sort(nprim, xsort)
    USE GAUSS
    implicit none
    integer(KINT) :: i, iux, j, nprim
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
    write(6,"('Error in subroutine sort. Highest number of steps (',i4,') exceeded')") mxsteps
    call error(1,'Stop.')
    return
    end
!
!	-------------------------------------------------------------------------------------------------------
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
