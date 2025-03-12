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
! Program for computing oriented multipolar moments of atomic fragments from the representation of the molecular density performed with
! DAM320
!
!
! Version of September 2018
!
  program DAMMULTROT320
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT), parameter :: nprt = 5, mxncntab = 100
    integer(KINT) :: i, i1, i2, i3, ierr, j, k, l, lm, lm0, lmax, lmin, m, ncntab
    integer(KINT) :: icntab(mxncntab+3)
    character(60) :: filename
    real(4) :: tarray(2), tiempo, dtime
    real(KREAL) :: alfa, aux, beta, cosag, cosal, cosamg, cosbet, cosga, gamma, sinag, sinal, sinamg, singa, sinbet
    real(KREAL) :: sum1, sum2, sum3, r21(3), r23(3), rip(3), rjp(3), rkp(3), un(3), vaux(3), rmlt(nprt)
    character(150) :: string
    namelist / options / filename, i1, i2, i3, iswindows, lmin, lmax, ncntab, icntab, rip

!	Namelist default values
    iswindows = .false.		! .true. if running on a MS-windows system
    i1 = 1
    i2 = 2
    i3 = 3
    lmin = 0		! Lowest l for multipolar moments
    lmax = 0
    rip = cero
    ncntab = 0
    icntab(:) = 0
    filename = ""		! root file name
!	End of namelist defaults

    tiempo = dtime(tarray)

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
    if (ncntab .gt. mxncntab) then
        write(6,"('Number of centers for tabulation = ', i4, ' higher than maximum allowable = ', i4)") ncntab, mxncntab
        write(6,"('Takes the maximum allowed value')")
        write(6,"('To increase maximum, change parameter mxncntab in file DAMMULTROT320.F90 and recompile')")
        ncntab = mxncntab
    endif
    icntab(ncntab+1) = i1
    icntab(ncntab+2) = i2
    icntab(ncntab+3) = i3
    ncntab = ncntab+3
    call consta
    call leedamqtmult

    if (ncen .lt. 3) then
            call error(1,'Number of centers too small. Must be equal or higher than 3. Stop')
    endif
    if (min(i1,i2,i3) .le. 0) then
            call error(1,'Wrong indices. Indices must be greater than zero. Stop')
    endif
    if (max(i1,i2,i3) .gt. ncen) then
            call error(1,'Wrong indices. Indices must be less or equal than number of centers. Stop')
    endif
    if (min(abs(i1-i2),abs(i1-i3),abs(i2-i3)) .eq. 0) then
            call error(1,'Wrong indices. Indices must not coincide. Stop')
    endif
    lmax = min(lmax,lmaxexp)
    allocate(rl(-lmax:lmax,-lmax:lmax,0:lmax))
    if (.not. allocated(rl)) call error(1,'Memory error when allocating rl. Stop')
    if (longoutput) write(6,"('Size of rl   = ', i15, ' bytes')") size(rl)

    allocate(dl(-lmax:lmax,-lmax:lmax,0:lmax))
    if (.not. allocated(dl)) call error(1,'Memory error when allocating dl. Stop')
    if (longoutput) write(6,"('Size of dl   = ', i15, ' bytes')") size(dl)


!	Defines a new axis frame (right-handed) with the Z axis perpendicular to plane defined by the centers
!    The X' axis can be user-supplied in options by means of array rip; otherwise the Y axis is defined
!	as the bisector of the angle formed by centers of indices i1, i2, and i3
    r21(1:3) = rcen(1:3,i1) - rcen(1:3,i2)	! r(i1) - r(i2)
    r23(1:3) = rcen(1:3,i3) - rcen(1:3,i2)	! r(i3) - r(i2)

!	Unitary vector of Z' axis 
    rkp(1) = r23(2)*r21(3) - r23(3)*r21(2)
    rkp(2) = r23(3)*r21(1) - r23(1)*r21(3)
    rkp(3) = r23(1)*r21(2) - r23(2)*r21(1)
    if(dot_product(rkp,rkp) .lt. 1.d-10) then
            call error(1,'The centers a alligned and do not define a plane. Stop')
    endif

    rkp = rkp / sqrt(dot_product(rkp,rkp))

    if (dot_product(rip,rkp) .gt. 1.d-10) then
        write(6,"('Vector rip is not orthogonal to rkp. Tries to take an orthogonal component')")
        rip = rip - rkp * dot_product(rip,rkp)
    endif

    if (dot_product(rip,rip) .lt. 1.d-10) then
! 		Unitary vector of Y' axis
        rjp(1:3) = rcen(1:3,i3) + rcen(1:3,i1) - dos * rcen(1:3,i2)
        rjp = rjp / sqrt(dot_product(rjp,rjp))
! 		Unitary vector of X' axis 
        rip(1) = rjp(2)*rkp(3) - rjp(3)*rkp(2)
        rip(2) = rjp(3)*rkp(1) - rjp(1)*rkp(3)
        rip(3) = rjp(1)*rkp(2) - rjp(2)*rkp(1)
    else
! 		Unitary vector of X' axis
        rip = rip / sqrt(dot_product(rip,rip))
! 		Unitary vector of Y' axis
        rjp(1) = rkp(2)*rip(3) - rkp(3)*rip(2)
        rjp(2) = rkp(3)*rip(1) - rkp(1)*rip(3)
        rjp(3) = rkp(1)*rip(2) - rkp(2)*rip(1)
    endif

!	Unitary vector along the line of nodes
    un(1) = -rkp(2)
    un(2) = rkp(1)
    un(3) = cero

!	Trigonometric functions of Euler angles for the rotation from the original frame to the new frame
    if (dot_product(un,un) .lt. 1.d-10) then
        cosal = rip(1)
        sinal = rip(2)
        sinbet = cero
        cosbet = uno
        singa = cero
        cosga = uno
        alfa = acos(cosal)
        beta = cero
        gamma = cero
    else
        un = un / sqrt(dot_product(un,un))
        sinbet = -un(1)*rkp(2) + un(2)*rkp(1)
        cosbet = max(-uno,min(uno,rkp(3))) ! To prevent numerical issues
        cosal  = un(2)
        sinal  = -un(1)
        cosga  = max(-uno,min(uno,dot_product(rjp,un))) ! To prevent numerical issues
        vaux(1) = un(2)*rjp(3) - un(3)*rjp(2)
        vaux(2) = un(3)*rjp(1) - un(1)*rjp(3)
        vaux(3) = un(1)*rjp(2) - un(2)*rjp(1)
        singa  = max(-uno,min(uno,dot_product(vaux,rkp))) ! To prevent numerical issues
        alfa = acos(cosal)
        if (sinal .lt. 0) alfa = dos*pi-alfa
        beta = acos(cosbet)
        gamma = acos(cosga)
        if (singa .lt. 0) gamma = dos*pi-gamma
    endif
    write(6,"(/,'Euler angles in radians: alfa = ', e19.12, ' beta = ', e19.12, ' gamma = ', e19.12)") alfa, beta, gamma
    write(6,"('Euler angles in degrees: alfa = ', f12.2, 8x, 'beta = ', f12.2, 8x, 'gamma = ', f12.2)") &
            alfa * 180.d0 / pi, beta * 180.d0 / pi, gamma * 180.d0 / pi

!	Builds the rotation matrices
    call rotar(lmax, cosal, sinal, cosbet, sinbet, cosga, singa)

    write(6,"(//'The following multipolar moments are referred to normalized spherical harmonics and defined in the frame:', &
            /7x, 'Z'' axis perpendicular to the plane defined by centers: (', a,1x,i3,', ',a,1x,i3,', ', a,1x,i3,'):', &
            2x,' k'' = ', 3(1x,e17.10), &
            /7x, 'Y'' axis in the bisector of angle: (', a,1x,i3,', ',a,1x,i3,', ', a,1x,i3,'):',  &
            21x, 2x, ' j'' = ', 3(1x,e17.10), &
            /7x, 'X'' axis orthogonal to Y and Z (right-handed system): ', 29x, ' i'' = ', 3(1x,e17.10)) ") &
            atmnam(i1), i1, atmnam(i2), i2, atmnam(i3), i3, rkp, &
            atmnam(i1), i1, atmnam(i2), i2, atmnam(i3), i3, rjp, rip
    call sort(ncntab, icntab)
!	Rotates the multipolar moments (in terms of normalized spherical harmonics) and prints out the results
    allocate(ang((lmax+1)*(lmax+2)/2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ang in densexactaGTO. Stop')
!    ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
    ang(1) = umed / raizpi
    lm = 1
    do l = 1, lmax
        lm = lm + 1
        ang(lm) = ang(1) * sqrt(re(2*l+1))
        aux = ang(lm) * raiz2
        do m = 1, l
            lm = lm + 1
            aux = aux / sqrt(re(l-m+1)*re(l+m))
            ang(lm) = aux
        enddo
    enddo
    do l = lmin, lmax
        lm0 = l*l
        string = "(/,'Multipolar moments of order l = ', i3)"
        write(6,string) l
        do k = 1, ncntab/nprt+1
            if (nprt*(k-1) .ge. ncntab) exit
            string = "(/3x,'m'"
            do i1 = (k-1)*nprt+1, min(ncntab,nprt*k)
                string = trim(string) // ",11x,'atom: ',a,1x,i3,2x"
            enddo
            string = trim(string) // ")"
            write(6,string) (atmnam(icntab(i1)), icntab(i1), i1 = (k-1)*nprt+1, min(ncntab,nprt*k))
            do i = -l, l
                rmlt = cero
                do j = -l, l
                    i2 = 0
                    do i1 = (k-1)*nprt+1, min(ncntab,nprt*k)
                        i2 = i2 + 1
                        rmlt(i2) = rmlt(i2) + rl(j,i,l) * rmultip(j+l+lm0+1,icntab(i1)) / ang(l*(l+1)/2+abs(j)+1)
                    enddo
                enddo
                write(6,"(i4,2x,5(3x,e22.15))") i, rmlt(1:i2)
            enddo
        enddo
    enddo
    tiempo = dtime(tarray)
    write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tiempo
    stop
9998 write(6,"(1x,'Input data error')")
    stop
9999 continue
    write(6,"(1x,'Input data ended before completed')")
    stop
    end
	
!
!	***************************************************************
!
  subroutine leedamqtmult
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMDEN320_D, only: icfposd
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, ierr, indnf, indng, interv, j, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lldummy, lm, m, ncenbas, ncflm, ndummy, nfdummy, nginidummy, ngfindummy, nndummy, nsamples
    real(KREAL) :: aux, bux, dummy, xxdummy
    logical :: latomtype(103) = .false.
    integer(KINT), allocatable :: ncfaj(:)
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
    write(6,"('ncen = ', i3, ' nbas = ', i5, ' nshells = ', i3)") ncen, nbas, ncaps
    !	Geometry and nuclear charges

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

!	Basis set
    i = 0
    read(10) lsto	! .true. means STO basis, .false. means GTO basis
    if (lsto) then
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')
!		Reads basis set data to dummies
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfindummy
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfindummy
                i = i + 1
                read(10) ndummy, ndummy, ndummy, xxdummy
            enddo
        enddo
    else
        read(10) nprimitot
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')
!		Reads basis set data to dummies
        indng = 1
        ncenbas = 0
        do ia = 1, ncen
            read(10) ndummy
            if (ndummy .le. 0) then
                ngini(ia) = -1
                cycle
            endif
            ncenbas = ncenbas + 1
            ngini(ia) = indng
            indng = indng + ndummy
            do j = 1, ndummy
                read(10) nndummy, lldummy
                read(10) (xxdummy, k = 1, nndummy)
                read(10) (xxdummy, k = 1, nndummy)
            enddo
        enddo
    endif

!	Data of density representation
    read(10) lmaxexp
    lmtop = (lmaxexp+1)*(lmaxexp+1)
    if (longoutput) write(6,"('lmaxexp = ', i2, ' nintervaj = ', i2)") lmaxexp, nintervaj

    allocate(icfposd(lmtop*nintervaj+1,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating icfposd. Stop')
    if (longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    if (longoutput) write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
    icfposd = 0
    k = 0
    do ia = 1, ncen      ! Do over centers
        if (ngini(ia) .le. 0) cycle
        read(10) icfposd(1:lmtop*nintervaj+1,ia)
        if (k .gt. 0) then
            icfposd(1:lmtop*nintervaj+1,ia) = icfposd(1:lmtop*nintervaj+1,ia) + icfposd(lmtop*nintervaj+1,k) - 1
        endif
        k = ia
        read(10) (xxdummy, i = 1, nintervaj)	! Reads xajust to dummy
!     fitting coeficients
        read(10) (xxdummy, i = icfposd(1,ia), icfposd(lmtop*nintervaj+1,ia)-1) ! Reads cfajust to dummy
    enddo

    allocate(rmultip(lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rmultip. Stop')
    if (longoutput) write(6,"('Size of rmultip    = ', i15, ' bytes')") size(rmultip)

    do ia = 1, ncen      ! Do over centers
        read(11) rmultip(1:lmtop,ia)	! multipolar moments
        read(11) (xxdummy, i = 1,nintervaj*lmtop)	! QGpart
        read(11) (xxdummy, i = 1,nintervaj*lmtop)	! qppart
        read(11) (xxdummy, i = icfposd(1,ia), icfposd(lmtop*nintervaj+1,ia)-1)	! crint1
        read(11) (xxdummy, i = icfposd(1,ia), icfposd(lmtop*nintervaj+1,ia)-1)	! cfrint2l2
    enddo
    close(10)
    close(11)
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
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMFIELD320_D
    integer(KINT) :: i, ierr
! 	implicit double precision (a-h,o-z)
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
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)   		!  i!
    enddo
    root(0) = cero
    do i = 1, mxroot
        root(i) = sqrt(re(i))        !  sqrt(i)
        rooti(i) = uno / root(i)     !  uno / sqrt(i)
    enddo
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
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    implicit none
    integer(KINT) :: l, lmax, ltot
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
        call dlmn(l, ltot, sinal, cosal, cosbet, tgbet2, singa, cosga)
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
  subroutine dlmn(l, ltot, sinal, cosal, cosbet, tgbet2, singa, cosga)
    USE DAM320_D
    USE DAM320_DATA_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: iinf, isup, l, ltot, m, mp
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
!	Subroutine sort: a silly subroutine for sorting indices of centers for tabulation
!

   subroutine sort(n, ix)
    USE DAM320_D
    implicit none
    integer(KINT) :: i, iux, j, k, n
    real(KREAL) :: aux
    integer(KINT), parameter :: mxsteps = 1000
    integer(KINT) :: ix(*)
    logical :: lend
    do i = 1, mxsteps
        lend = .true.
        do j = n, 2, -2
            if (ix(j) .lt. ix(j-1)) then
                aux = ix(j)
                ix(j) = ix(j-1)
                ix(j-1) = aux
                lend = .false.
            endif
        enddo
        do j = n-1, 2, -2
            if (ix(j) .lt. ix(j-1)) then
                aux = ix(j)
                ix(j) = ix(j-1)
                ix(j-1) = aux
                lend = .false.
            endif
        enddo
        if (lend) then
!			suppresses repeated indices
            k = 1
            do j = 2, n
                if (ix(j) .ne. ix(k) .and. ix(j) .ne. 0) then
                    k = k+1
                    ix(k) = ix(j)
                endif
            enddo
            n = k
            return
        endif
    enddo
    write(6,"('Error in subroutine sort. Highest number of steps (',i4,') exceeded')") mxsteps
    call error(1,'Stop.')
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
