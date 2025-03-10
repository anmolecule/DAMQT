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
! Program for computing the Hellman-Feynman forces from the representation of the molecular density performed with
! DAM320
!
!
! Version of April 2019
!
  program DAMFORCES320
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMFIELD320_D
    USE DAMFORCES320_D
    implicit none
    integer(KINT) :: i, j, k, ncntab
!     implicit double precision (a-h,o-z)
    real(4) :: tarray(2), tiempo, dtime
    namelist / options / filename, iatomsel, iswindows, latomsel, lmaxrep, longoutput, lvalence, ncntab, umbrlargo
!    Namelist default values
    longoutput = .false.    ! If .true. a more detailed output is given
    latomsel = .false.      ! If .true. indices of the centers for atomic tabulations will be user-supplied
                            ! The indices of the selected atoms must be supplied in vector "iatomsel".
    lvalence = .false.      ! If .true. only valence electrons are considered
    lmaxrep = 10            ! highest "l" for the computation of forces
    ncntab = 1              ! Number of atoms selected for tabulation (dummy, only used by DAMQT GUI)
    iatomsel(1) = 1         ! To cover the possibility of an input file with "latomsel = .true."
    iatomsel(2:mxsel) = 0   ! but without assignment of "iatomsel".
    umbrlargo = 1.d-8       ! Long-range threshold
    iswindows = .false.     ! .true. if running on a MS-windows system
    filename = ""           ! root file name
!    End of namelist defaults

    tiempo = dtime(tarray)
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

    if (len_trim(filename).eq.0) then
        filename = projectname
    else
        filename = projectname(1:i)//trim(filename)
    endif

    call consta

    call leedamqtforces

    if (latomsel) then
        k = 0
        doi: do i = 1, mxsel
            if (iatomsel(i) .le. ncen .and. iatomsel(i) .gt. 0) then
                do j = 1, k
                    if (iatomsel(i) .eq. iatomsel(j)) cycle doi
                enddo
                k = k+1
                iatomsel(k) = iatomsel(i)
            endif
        enddo doi
        nsel = k
    endif

!    Field and forces on nuclei
    call fuerzas
    tiempo = dtime(tarray)
    write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") tarray(1), tarray(2), tiempo
    stop
    end
    
!
!    ***************************************************************
!
  subroutine leedamqtforces
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMFIELD320_D
    USE DAMFORCES320_D
    USE GAUSS
    integer(KINT) :: i, ia, icarga, ierr, indnf, indng, interv, j, jshft, k, k1, k2, knt, kntlm
    integer(KINT) :: l, lenindintrv, lldummy, lm, m, ncenbas, ncflm, ndummy, nfdummy, nginidummy, ngfindummy, nndummy, nsamples
    real(KREAL) :: aux, bux, dltsample, dost, dummy, fr1, fr2l2, pi4d2l1, r, ra, ral1inv, rainv, ral
    real(KREAL) :: rinta, rintb, rlarex, step, stepmed, suml, suml1, suml2, summ, summ1, summ2, t, xxdummy
    real(KREAL), allocatable :: Qg(:), qp(:)
    real(KREAL) :: tcheb(0:mxlenpol-1)
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
    if (longoutput) write(6,"('Opens file ', a)") trim(projectname)//"_2016.damqt"
    read(10) ncen, nbas, ncaps
    write(6,"('ncen = ', i3, ' nbas = ', i5, ' nshells = ', i3)") ncen, nbas, ncaps

!    Allocates memory for geometry and basis set

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating atmnam. Stop')

    allocate(ncontr(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating ncontr. Stop')

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating nzn. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rcen. Stop')

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating zn. Stop')

!    Geometry and nuclear charges

    knt = 0
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

!    Basis set
    i = 0
    read(10) lsto    ! .true. means STO basis, .false. means GTO basis
    if (lsto) then
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')
!        Reads basis set data to dummies
        write(6,"(/t22,'STO Basis set',/t22,13('-'))")
        i = 0
        ncenbas = 0
        do ia = 1, ncen
            read(10) ngini(ia), ngfindummy
            if (ngini(ia) .le. 0) cycle
            ncenbas = ncenbas + 1
            do k = ngini(ia), ngfindummy
                i = i + 1
                read(10) nfdummy, nndummy, lldummy, xxdummy
            enddo
        enddo
    else
        write(6,"(/t22,'GTO Basis set',/t22,13('-'))")
        read(10) nprimitot
        allocate(ngini(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating ngini. Stop')
!        Reads basis set data to dummies
        indng = 1
        ncenbas = 0
        do ia = 1, ncen
            read(10) ncontr(ia)
            if (ncontr(ia) .le. 0) then
                ngini(ia) = -1
                cycle
            endif
            ncenbas = ncenbas + 1
            ngini(ia) = indng
            indng = indng + ncontr(ia)
            do j = 1, ncontr(ia)
                read(10) nndummy, lldummy
                read(10) (xxdummy, k = 1, nndummy)
                read(10) (xxdummy, k = 1, nndummy)
            enddo
        enddo
    endif

!    Data of density representation
    read(10) lmaxexp
    lmtop = (lmaxexp+1)*(lmaxexp+1)
    if (lmaxrep .gt. lmaxexp) then
        write(6,"('lmaxrep = ', i3, ' greater than lmaxexp ', i3)") lmaxrep, lmaxexp
        write(6,"('takes lmaxrep = ',i3)") lmaxexp
        lmaxrep = lmaxexp
    else
        write(6,"('Highest l in the computation of forces = ', i3)") lmaxrep
    endif

    allocate(icfposd(lmtop*nintervaj+1,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating icfposd. Stop')
    if (longoutput) write(6,"('Size of icfposd   = ', i15, ' bytes')") size(icfposd)
    if (longoutput) write(6,"('radii of fitting intervals: ',/, 8(1x,e17.10))") rinterv
    icfposd = 0
!    xajustd = cero
!    cfajust = cero
    k = 0
    do ia = 1, ncen      ! Do over centers
        if (.not. lsto .and. ncontr(ia) .le. 0) cycle
        read(10) icfposd(1:lmtop*nintervaj+1,ia)
        if (k .gt. 0) icfposd(1:lmtop*nintervaj+1,ia) = icfposd(1:lmtop*nintervaj+1,ia) + icfposd(lmtop*nintervaj+1,k) - 1
        k = ia
        read(10) (xxdummy, i = 1, nintervaj)    ! Reads xajust to dummy
!     fitting coeficients
        read(10) (xxdummy, i = icfposd(1,ia), icfposd(lmtop*nintervaj+1,ia)-1) ! Reads cfajust to dummy
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

    allocate(qptotal(lmtop,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating qptotal. Stop')

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
        qptotal(1:lmtop,ia) =  qpacum(1:lmtop,ia) + qppart(1:lmtop)
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
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating umedpow in consta. Stop')
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
    deallocate (Qgpart, qppart, umedpow)
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
!        ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) 
!            / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
!
!**********************************************************************
  subroutine consta
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE DAMFIELD320_D
    integer(KINT) :: i, ierr
!     implicit double precision (a-h,o-z)
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
        dosl1i(-i) = -dosl1i(i)
    enddo
    fact(0) = uno
    facti(0) = uno
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)           !  i!
        facti(i) = uno / fact(i)             !  uno / i!
    enddo
    return
    end

!   ***************************************************************

  subroutine fuerzas
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D, zn_orig => zn
    USE DAMFIELD320_D
    USE DAMFORCES320_D
    implicit none
    integer(KINT) :: i, ia, ib, icflm, ierr, interv, isgm, isgmm1, isgmp1, j, k, knt, l, lm, lm1l, lp1lp2
    integer(KINT) :: ltop, m, mm, ndimt, nord
    real(KREAL), parameter :: umbrlambda = 1.e-13
    real(KREAL) :: aux, c1, c2, dost, Ex, Exext, Exint, Exnucl, Ey, Eyext, Eyint, Eynucl, Ez, Ezext, Ezint, Eznucl
    real(KREAL) :: fconfx, fconfy, fconfz, fnoconfx, fnoconfy, fnoconfz, fxexte, fxextn
    real(KREAL) :: fyexte, fyextn, fzexte, fzextn
    real(KREAL) :: pi4d2l1, rab, rab2, rab2l1, rab2l1inv, rab2i, rab3inv, rab2l3inv, rinta, rintb
    real(KREAL) :: t, tconfx, tconfy, tconfz, tnoconfx, tnoconfy, tnoconfz
    real(KREAL) :: x, xab, vext, vin, vnucl, vrad, y, yab, z, zab
    logical, allocatable :: lselaux(:)
    real(KREAL) :: zlm((lmaxexp+2)**2)
    real(KREAL) :: rk(3,3), tt(6,6), ttvec(6,6), ttinv(6,6), ttval(6)
    real(KREAL), allocatable :: Exatome(:), Eyatome(:), Ezatome(:), Exatomn(:), Eyatomn(:), Ezatomn(:)
    real(KREAL), allocatable :: fnoconf(:), fconf(:), fxext(:), fxint(:), fxtot(:), fyext(:), fyint(:), fytot(:)
    real(KREAL), allocatable :: fzext(:), fzint(:), fztot(:), tmat(:,:), taux(:,:), tproy(:,:)
    real(KREAL), allocatable :: xangst(:), yangst(:), zangst(:)
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

    allocate(lselaux(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating lselaux. Stop')

    allocate(Exatome(ncen), Eyatome(ncen), Ezatome(ncen), Exatomn(ncen), Eyatomn(ncen), Ezatomn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating Exatome, Eyatome, Ezatome, Exatomn, Eyatomn, Ezatomn. Stop')

    allocate(xangst(ncen), yangst(ncen), zangst(ncen))
    if (ierr .ne. 0) call error(1,'Memory error when allocating xangst, yxangst, zangst. Stop')

    allocate(fnoconf(3*ncen), fconf(3*ncen), fxext(ncen), fxint(ncen), fxtot(ncen), &
            fyext(ncen), fyint(ncen), fytot(ncen), fzext(ncen), fzint(ncen), fztot(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating fnoconf, fconf, fxext, fxint, ..., fztot. Stop')

    allocate(tmat(3*ncen,6), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating tmat. Stop')

    allocate(taux(6,3*ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating taux. Stop')

    allocate(tproy(3*ncen,3*ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating tproy. Stop')

    lmtop = (lmaxexp+1) * (lmaxexp+1)
!    Opens file for Hellman-Feynman forces
    open(14,file=trim(filename)//".forces",status='unknown',form='formatted',iostat=ierr)
    if (ierr .ne. 0) call error(1,'Error when opening file '//trim(filename)//".forces"//'.Stop')
    write(14,*) ncen

    lmtop = (lmaxexp+1)*(lmaxexp+1)
    tmat = cero
    tt = cero
    tt(1,1) = ncen
    tt(2,2) = ncen
    tt(3,3) = ncen
    if (latomsel) then
        do ia = 1, ncen
            lselaux(ia) = .false.
            do j = 1, nsel
                if (ia .eq. iatomsel(j)) then
                    lselaux(ia) = .true.
                    exit
                endif
            enddo
        enddo
    else
            do ia = 1, ncen
                    lselaux(ia) = .false.
            enddo
    endif
    WRITE(6,"('Electrostatic potential, field and forces on nuclei computed with an atomic expansion with l up to ', i3)") lmaxexp
    do ia = 1, ncen
            write(6,"(//3x,18('*'),'  center:',1X,i3,1x,A,3x,20('*'))") ia, atmnam(ia)
!        Potential and field of distribution centered at  ia  over this center  ia
        vnucl = cero
        Exnucl = cero
        Eynucl = cero
        Eznucl = cero
        vin = -cuatro * pi * qptotal(1,ia)
        Exint = 1.333333333333333d0 * pi * qptotal(4,ia)
        Eyint = 1.333333333333333d0 * pi * qptotal(2,ia)
        Ezint = 1.333333333333333d0 * pi * qptotal(3,ia)
!        Potential and field of distribution placed at centers other than  ia  over center  ia
        vext = cero
        Exext = cero
        Eyext = cero
        Ezext = cero
        Exatome(ia) = cero
        Eyatome(ia) = cero
        Ezatome(ia) = cero
        Exatomn(ia) = cero
        Eyatomn(ia) = cero
        Ezatomn(ia) = cero
        do ib = 1, ncen
            if (ia .eq. ib) cycle   ! If ia = ib skips to next
            xab = rcen(1,ia) - rcen(1,ib)
            yab = rcen(2,ia) - rcen(2,ib)
            zab = rcen(3,ia) - rcen(3,ib)
            rab2 = xab*xab + yab*yab + zab*zab
            rab = dsqrt(rab2)
            if (rab .lt. 1.d-10) then
                    write(6,"('The center ',i4,' coincides with center ',i4)") ib, ia
                    return
            endif
            rab2i = uno / rab2
            rab3inv = uno / (rab * rab2)
            Exnucl = Exnucl + zn(ib) * xab * rab3inv
            Eynucl = Eynucl + zn(ib) * yab * rab3inv
            Eznucl = Eznucl + zn(ib) * zab * rab3inv
            if (rab .lt. rlargo(ib)) then
                interv = indintrv(int(fct*rab)+1)
                ltop = lcorto(interv,ib)
            else
                ltop = llargo(min(int(rab),mxlargo),ib)
            endif
            vnucl = vnucl + zn(ib) / rab
            if (lselaux(ia)) then
                Exatomn(ib) = zn(ib) * xab * rab3inv
                Eyatomn(ib) = zn(ib) * yab * rab3inv
                Ezatomn(ib) = zn(ib) * zab * rab3inv
                Exatome(ib) = cero
                Eyatome(ib) = cero
                Ezatome(ib) = cero
            endif
            if (ngini(ib) .le. 0) cycle
            zlm(1) = uno        ! Regular spherical harmonics of R(ia)-R(ib)
            zlm(2) = yab
            zlm(3) = zab
            zlm(4) = xab
            do l = 1, ltop
                zlm((l+1)*(l+3)+1) = dosl1(l) * (xab * zlm(l*(l+2)+1) - yab * zlm(l*l+1))        ! zlm(l+1,l+1)
                zlm((l+1)*(l+1)+1) = dosl1(l) * (yab * zlm(l*(l+2)+1) + xab * zlm(l*l+1))        ! zlm(l+1,-(l+1))
                zlm((l+2)*(l+2)-1) = dosl1(l) * zab* zlm(l*(l+2)+1)                ! zlm(l+1,l)
                zlm(l*(l+2)+3) = dosl1(l) * zab * zlm(l*l+1)                    ! zlm(l+1,-l)
                do m = 0, l-1
                    zlm((l+1)*(l+2)+m+1) = ri(l-m+1) * (dosl1(l)*zab*zlm(l*(l+1)+m+1) - re(l+m)*rab2*zlm((l-1)*l+m+1))    ! zlm(l+1,m)
                    zlm((l+1)*(l+2)-m+1) = ri(l-m+1) * (dosl1(l)*zab*zlm(l*(l+1)-m+1) - re(l+m)*rab2*zlm((l-1)*l-m+1))    ! zlm(l+1,-m)
                enddo
            enddo
            if (rab .gt. rlargo(ib)) then    ! Long-range
                lm = 0
                rab2l1inv = uno / rab
                rab2l3inv = rab2l1inv * rab2i
                do l = 0, ltop
                    do m = -l, l
                        lm = lm + 1
                        if(abs(QGacum((nintervaj-1)*lmtop+lm,ib)) .lt. umbrlargo) cycle
                        mm = abs(m)
                        if (m .ge. 0) then
                            isgm = 1
                        else
                            isgm = -1
                        endif
!                         Potential
                        vext = vext - rmultip(lm,ib) * zlm(l*(l+1)+m+1) * rab2l1inv
!                         Field
!                        Ex
                        c1 = zlm( (l+1)*(l+2)+isgm*(mm+1)+1)
                        if (m .ne. 0 .and. m .ne. -1) then
                            c1 = c1 - re(l+1-mm) * re(l+2-mm) * zlm( (l+1)*(l+2)+isgm*(mm-1)+1 )
                        else
                            if (m .eq. 0) then
                                c1 = c1 + c1
                            endif
                        endif
                        c1 = c1 * rab2l3inv
                        Ex = -umed * rmultip(lm,ib) * c1
                        Exext = Exext + Ex
!                        Ey
                        c1 = zlm( (l+1)*(l+2)-isgm*(mm+1)+1 )
                        if (m .ne. 0 .and. m .ne. 1) then
                            c1 = c1 + re(l+1-mm) * re(l+2-mm) * zlm( (l+1)*(l+2)-isgm*(mm-1)+1 )
                        else
                            if (m .eq. 0) then
                                c1 = c1 + c1
                            endif
                        endif
                        c1 = c1 * rab2l3inv
                        Ey = -isgm * umed * rmultip(lm,ib) * c1
                        Eyext = Eyext + Ey
!                        Ez
                        Ez = -rmultip(lm,ib)  * re(l+1-mm) * zlm( (l+1)*(l+2)+m+1 ) * rab2l3inv
                        Ezext = Ezext + Ez
                        if (lselaux(ia)) then
                            Exatome(ib) = Exatome(ib) + Ex
                            Eyatome(ib) = Eyatome(ib) + Ey
                            Ezatome(ib) = Ezatome(ib) + Ez
                        endif
                    enddo
                    rab2l1inv = rab2l1inv * rab2i
                    rab2l3inv = rab2l3inv * rab2i
                enddo
                kntlargo = kntlargo + 1
            else        ! Short-range
!                Seeks the fitting interval for rab
                interv = indintrv(int(fct*rab)+1)
                t = dos * (rab - rinterv(interv-1))/(rinterv(interv)-rinterv(interv-1)) - uno
                dost = t + t
                tcheb(0) = uno    ! Chebyshev T  polynomials
                tcheb(1) = t
                do j = 2, mxlenpol-1
                    tcheb(j) = dost * tcheb(j-1) - tcheb(j-2)
                enddo
!                 Integrals:
!                    rinta = Integrate[ r**(2*la+2) * fradtr[la,ma,r], {r,l_(i-1),r0}]
!                    rintb = Integrate[ r * fradtr[la,ma,r], {r,r0,l_i}]
!                computed from fitting polynomials
                lm = 0
                rab2l1inv = uno / rab
                rab2l3inv = rab2l1inv * rab2i    ! 1.d0 / (ra**(2*l+3))
                rab2l1 = rab
                do l = 0, ltop
                    pi4d2l1 = cuatro * pi * dosl1i(l)
                    lp1lp2 = (l+1)*(l+2)
                    lm1l = (l-1)*l
                    do m = -l, l
                        lm = lm + 1
                        if(abs(QGacum((nintervaj-1)*lmtop+lm,ib)) .lt. umbrlargo) cycle
                        icflm = icfposd((interv-1)*lmtop+lm,ib)
                        rinta = cero
                        rintb = cero
                        do i = 0, icfposd((interv-1)*lmtop+lm+1,ib)-icflm-1
                            rinta = rinta + cfrint2l2(icflm+i) * tcheb(i)
                            rintb = rintb + cfrint1(icflm+i) * tcheb(i)
                        enddo
                        rinta = rinta * (rab-rinterv(interv-1))
                        if (interv .gt. 1) rinta = QGacum((interv-2)*lmtop+lm,ib) + rinta
                        rintb = qpacum((interv-1)*lmtop+lm,ib) + rintb * (rinterv(interv)-rab)
                        vrad = pi4d2l1 * ( rinta + rab2l1 * rintb)
!                        Potential
                        vext = vext - vrad * zlm(l*(l+1)+m+1) * rab2l1inv
!                        Field
                        mm = iabs(m)
                        if (m .ge. 0) then
                            isgm = 1
                        else
                            isgm = -1
                        endif
                        isgmm1 = isgm*(mm-1)
                        isgmp1 = isgm*(mm+1)
!    Ex
                        c1 = zlm( lp1lp2+isgmp1+1 )
                        if (l .gt. mm+1) then
                            c2 = zlm(lm1l + isgmp1 + 1)
                        else
                            c2 = cero
                        endif
                        if (m .ne. 0 .and. m .ne. -1) then
                            c1 = c1 - re(l+1-mm) * re(l+2-mm) * zlm( lp1lp2+isgmm1+1 )
                            c2 = c2 - re(l+mm) * re(l+mm-1) * zlm( lm1l+isgmm1+1 )
                        else
                            if (m .eq. 0) then
                                c1 = c1 + c1
                                c2 = c2 + c2
                            endif
                        endif
                        Ex = - umed * pi4d2l1 * (c1 * rinta * rab2l3inv + c2 * rintb)
                        Exext = Exext + Ex
!      Ey
                        c1 = zlm( lp1lp2-isgmp1+1 )
                        if (l .gt. mm+1) then
                            c2 = zlm(lm1l-isgmp1+1 )
                        else
                            c2 = cero
                        endif
                        if (m .ne. 0 .and. m .ne. 1) then
                            c1 = c1 + re(l+1-mm) * re(l+2-mm) * zlm( lp1lp2-isgmm1+1 )
                            c2 = c2 + re(l+mm) * re(l+mm-1) * zlm( lm1l-isgmm1+1 )
                        else
                            if (m .eq. 0) then
                                c1 = c1 + c1
                                c2 = c2 + c2
                            endif
                        endif
                        Ey = - isgm * umed * pi4d2l1 * (c1 * rinta * rab2l3inv + c2 * rintb)
                        Eyext = Eyext + Ey
!      Ez    
                        c1 = re(l+1-mm) * zlm( lp1lp2+m+1 )
                        if (l .gt. mm) then
                            c2 = re(l+mm) * zlm( lm1l+m+1 )
                        else
                            c2 = cero
                        endif
                        Ez = - pi4d2l1 * (c1 * rinta * rab2l3inv - c2 * rintb)
                        Ezext = Ezext + Ez
                        if (lselaux(ia)) then
                            Exatome(ib) = Exatome(ib) + Ex
                            Eyatome(ib) = Eyatome(ib) + Ey
                            Ezatome(ib) = Ezatome(ib) + Ez
                        endif
                    enddo
                    rab2l1inv = rab2l1inv * rab2i
                    rab2l3inv = rab2l3inv * rab2i
                    rab2l1 = rab2l1 * rab2
                enddo
            endif
        enddo     ! End of Do on ib
        write(6,"(/3x,'Electrostatic force (nuclei+electrons) on nucleus',1x,i3,/3x,57('-'))") ia
        fxint(ia) = zn(ia) * Exint
        fyint(ia) = zn(ia) * Eyint
        fzint(ia) = zn(ia) * Ezint
        fxextn = zn(ia) * Exnucl
        fyextn = zn(ia) * Eynucl
        fzextn = zn(ia) * Eznucl
        fxexte = zn(ia) * Exext
        fyexte = zn(ia) * Eyext
        fzexte = zn(ia) * Ezext
        fxext(ia) = zn(ia) * ( Exnucl + Exext )
        fyext(ia) = zn(ia) * ( Eynucl + Eyext )
        fzext(ia) = zn(ia) * ( Eznucl + Ezext )
        fxtot(ia) = zn(ia) * ( Exnucl + Exint + Exext )
        fytot(ia) = zn(ia) * ( Eynucl + Eyint + Eyext )
        fztot(ia) = zn(ia) * ( Eznucl + Ezint + Ezext )
        write(6,"(t28,'Fx',t41,'Fy',t54,'Fz')")
        write(6,"(/3x,'Internal force',t21,3f13.8)") fxint(ia),fyint(ia),fzint(ia)
        write(6,"(/3x,'External force')")
        write(6,"(6X,'Nuclei',t21,3f13.8)") fxextn,fyextn,fzextn
        write(6,"(6X,'Electron',t21,3f13.8)") fxexte,fyexte,fzexte
        write(6,"(6X,'Total',t21,3f13.8)") fxext(ia),fyext(ia),fzext(ia)
        write(6,"(/3x,'Total force',t21,3f13.8)") fxtot(ia),fytot(ia),fztot(ia)
        if (lselaux(ia)) then
            write(6,"(//3x,'Atomic components of external force',/3x, 57('-'))")
            do ib = 1, ncen
                if (ia .ne. ib) then
                    write(6,"(/3x,'Forces exerted by atom',1X,I3,1X,A,1X,'on nucleus',1X,I3,1X,A2)") ib, &
                                    atmnam(ib), ia, atmnam(ia)
                    write(6,"(t28,'Fx',t41,'Fy',t54,'Fz')")
                    write(6,"(6x,'Nuclei',t21,3(f13.8))") zn(ia) * Exatomn(ib), zn(ia) * Eyatomn(ib), zn(ia) * Ezatomn(ib)
                    write(6,"(6x,'Electron',t21,3(f13.8))") zn(ia) * Exatome(ib), zn(ia) * Eyatome(ib), zn(ia) * Ezatome(ib)
                    write(6,"(6x,'Total',t21,3(f13.8))") zn(ia) * (Exatomn(ib) + Exatome(ib)) &
                                    , zn(ia) * (Eyatomn(ib) + Eyatome(ib)), zn(ia) * (Ezatomn(ib) + Ezatome(ib))
                endif
            enddo
        endif
        write(6,"(//3x,'Electrostatic potential at center',1X,i3,/3x,57('-'))") ia
        write(6,"(/3x,'Electron')")
        write(6,"(6X,'Internal',t47,f13.8)") vin
        write(6,"(6X,'External',t47,f13.8)") vext
        write(6,"(6X,'Total',t47,f13.8)") vin + vext
        write(6,"(/3x,'Nuclei (excluding nucleus at center)',t47,f13.8)") vnucl
        write(6,"(/3x,'Total (excluding nucleus at center)',t47,f13.8)") vin + vext + vnucl
        write(6,"(//3x,'Electric field at center',1X,i3,/3x,57('-'))") ia
        write(6,"(t28,'Ex',t41,'Ey',t54,'Ez')")
        write(6,"(3x,'Electron')")
        write(6,"(6X,'Internal',t21,3(f13.8))") Exint, Eyint, Ezint
        write(6,"(6X,'External',t21,3(f13.8))") Exext, Eyext, Ezext
        write(6,"(6X,'Total',t21,3(f13.8))") Exint+Exext, Eyint+Eyext, Ezint+Ezext
        write(6,"(/3x,'Nuclei  (excluding nucleus at center)',/t21,3(f13.8))") Exnucl, Eynucl, Eznucl
        write(6,"(/3x,'Total field (excluding nucleus at center)',/t21,3(f13.8))") &
                Exint+Exext+Exnucl, Eyint+Eyext+Eynucl, Ezint+Ezext+Eznucl
        rk(1,1) = cero
        rk(1,2) = rcen(3,ia)
        rk(1,3) = -rcen(2,ia)
        rk(2,1) = -rcen(3,ia)
        rk(2,2) = cero
        rk(2,3) = rcen(1,ia)
        rk(3,1) = rcen(2,ia)
        rk(3,2) = -rcen(1,ia)
        rk(3,3) = cero
        tmat(3*(ia-1)+1,1) = uno
        tmat(3*(ia-1)+2,2) = uno
        tmat(3*(ia-1)+3,3) = uno
        do j = 1, 3
            do i = 1, 3
                tmat(3*(ia-1)+i,j+3) = rk(i,j)
                tt(i,j+3) = tt(i,j+3) + rk(i,j)
                tt(i+3,j) = tt(i+3,j) + rk(j,i)
                aux = 0.d0
                do k = 1, 3
                    aux = aux + rk(k,i) * rk(k,j)
                enddo
                tt(i+3,j+3) = tt(i+3,j+3) + aux
            enddo
        enddo
        xangst(ia) = rcen(1,ia) * 5.291772083d-1
        yangst(ia) = rcen(2,ia) * 5.291772083d-1
        zangst(ia) = rcen(3,ia) * 5.291772083d-1
    enddo     ! End of Do on ia
    do ia = 1, ncen
            write(14,"(6(2x,F12.5))") xangst(ia), yangst(ia), zangst(ia), fxext(ia), fyext(ia), fzext(ia)
    enddo
    do ia = 1, ncen
            write(14,"(6(2x,F12.5))") xangst(ia), yangst(ia), zangst(ia), fxint(ia), fyint(ia), fzint(ia)
    enddo
    do ia = 1, ncen
            write(14,"(6(2x,F12.5))") xangst(ia), yangst(ia), zangst(ia), fxtot(ia), fytot(ia), fztot(ia)
    enddo

!    Computes the conformacional and nonconformational components of the total forces (see work by Jaime Fernandez Rico)
    nord = -1        ! orders eigenvalues and eigenvectors in descending eigenvalues
    call jacobi( 6 , 6 , tt , ttvec , ttval , nord )
    ndimt = 0
    do i = 1, 6
        if (abs(ttval(i)) .lt. umbrlambda) exit
        ndimt = ndimt + 1
    enddo
    do i = 1, ndimt
        do j = 1, ndimt
            aux = 0.d0
            do k = 1, ndimt
                if (abs(ttval(k)) .ge. umbrlambda)  aux = aux + ttvec(i,k) * ttvec(j,k) / ttval(k)
            enddo
            ttinv(i,j) = aux
        enddo
    enddo
!    Projector  T (T+ T)^-1 T+
    do i = 1, ndimt
        do j = 1, 3*ncen
            aux = cero
            do k = 1, ndimt
                aux = aux + ttinv(i,k) * tmat(j,k)
            enddo
            taux(i,j) = aux
        enddo
    enddo
    do i = 1, 3*ncen
        do j = 1, 3*ncen
            aux = cero
            do k = 1, ndimt
                aux = aux + tmat(i,k) * taux(k,j)
            enddo
            tproy(i,j) = aux
        enddo
    enddo
    write(6,"(//3X,11('*'),1X,'FORCES DECOMPOSITION',1X,12('*'),/)")
!    Nonconformational forces
    do i = 1, 3*ncen
        aux = cero
        knt = 0
        do k = 1, 3*ncen, 3
            knt = knt + 1
            aux = aux + tproy(i,k) * fxtot(knt) + tproy(i,k+1) * fytot(knt) + tproy(i,k+2) * fztot(knt)
        enddo
        fnoconf(i) = aux
    enddo
    knt = 0
    fnoconfx = cero
    fnoconfy = cero
    fnoconfz = cero
    tnoconfx = cero
    tnoconfy = cero
    tnoconfz = cero
    knt = 0
    do i = 1, 3*ncen, 3
        knt = knt + 1
        fnoconfx = fnoconfx + fnoconf(i)
        fnoconfy = fnoconfy + fnoconf(i+1)
        fnoconfz = fnoconfz + fnoconf(i+2)
        tnoconfx = tnoconfx + rcen(2,knt) * fnoconf(i+2) - rcen(3,knt) * fnoconf(i+1)
        tnoconfy = tnoconfy + rcen(3,knt) * fnoconf(i) - rcen(1,knt) * fnoconf(i+2)
        tnoconfz = tnoconfz + rcen(1,knt) * fnoconf(i+1) - rcen(2,knt) * fnoconf(i)
        write(14,"(6(2x,F12.5))") xangst((i-1)/3+1), yangst((i-1)/3+1), zangst((i-1)/3+1), fnoconf(i), fnoconf(i+1), fnoconf(i+2)
    enddo
    write(6,"(/3X,'Non conformational forces (spurious)')")
    write(6,"(3X,57('-'))")
    write(6,"(3X,'Atom',t28,'Fx',t41,'Fy',t54,'Fz',/)")
    knt = 0
    do i = 1, 3*ncen, 3
            knt = knt + 1
            write(6,'(2x,i3,t21,3(f13.8))') knt, fnoconf(i), fnoconf(i+1), fnoconf(i+2)
    enddo
    write(6,"(/3X,'Total',t21,3(f13.8))") fnoconfx, fnoconfy, fnoconfz
    write(6,"(/3x,'Torque',t21,3(f13.8))") tnoconfx, tnoconfy, tnoconfz
!    Conformational forces
    fconfx = cero
    fconfy = cero
    fconfz = cero
    tconfx = cero
    tconfy = cero
    tconfz = cero
    knt = 0
    do i = 1, 3*ncen, 3
        knt = knt + 1
        fconf(i) = fxtot(knt) - fnoconf(i)
        fconf(i+1) = fytot(knt) - fnoconf(i+1)
        fconf(i+2) = fztot(knt) - fnoconf(i+2)
        fconfx = fconfx + fconf(i)
        fconfy = fconfy + fconf(i+1)
        fconfz = fconfz + fconf(i+2)
        tconfx = tconfx + rcen(2,knt) * fconf(i+2) - rcen(3,knt) * fconf(i+1)
        tconfy = tconfy + rcen(3,knt) * fconf(i) - rcen(1,knt) * fconf(i+2)
        tconfz = tconfz + rcen(1,knt) * fconf(i+1) - rcen(2,knt) * fconf(i)
        write(14,"(6(2x,F12.5))") xangst((i-1)/3+1), yangst((i-1)/3+1), zangst((i-1)/3+1), fconf(i), fconf(i+1), fconf(i+2)
    enddo
    WRITE(6,"(/3X,'Conformational forces')")
    WRITE(6,"(3X,57('-'))")
    write(6,"(3X,'Atom',t28,'Fx',t41,'Fy',t54,'Fz',/)")
    knt = 0
    do i = 1, 3*ncen, 3
        knt = knt + 1
        write(6,"(2x,i3,t21,3(f13.8))") knt, fconf(i), fconf(i+1), fconf(i+2)
    enddo
    write(6,"(/3X,'Total',t21,3(f13.8))") fconfx, fconfy, fconfz
    write(6,"(/3x,'Torque',t21,3(f13.8))") tconfx, tconfy, tconfz
    close(14)
    deallocate(Exatome, Eyatome, Ezatome, Exatomn, Eyatomn, Ezatomn, fnoconf, fconf, fxtot, fytot, fztot, tmat, taux, tproy)
    return
    end
!
!   *********************************************************
!
  subroutine jacobi( jdim , n , b , x , v , nord )
!
!  if nord = 1   does not order the eigenvalues and eigenvectors
!  if nord = -1  order the eigenvalues and eigenvectors in descending order of eigenvalues
!
    USE DAM320_D
    implicit none
!     implicit real*8 (a-h,o-z)
    real(KREAL), parameter :: tol = 1.d-15
    integer(KINT), parameter :: itmax = 10000
    real(KREAL) :: al, ap, aq, apq, aux, be, c, r, s, top, vv
    real(KREAL) :: b(jdim,jdim), x(jdim,jdim), v(jdim), a(jdim,jdim)
    integer(KINT) :: i, ip, it, iq, j, jdim, k, n, nord
    integer(KINT) :: iv(jdim)
    logical lsuccess, lordering
    if (n .eq. 1) then
        x(1,1) = 1.d0
        v(1) = b(1,1)
        return
    endif
    do i = 1,n
        do j = i,n
            a(i,j) = b(i,j)
            a(j,i) = b(j,i)
            x(i,j) = 0.d0
            x(j,i) = 0.d0
        enddo
        x(i,i) = 1.d0
    enddo
    do j = 1 , n
        top = 0.d0
        do i = 1 , n
            if ( j.ne.i ) then
                aux = abs(a(i,j))
                if (aux.gt.top) then
                    top = aux
                    iv(j) = i
                end if
            end if
        end do
        v(j) = top
    end do
    lsuccess = .false.
    do it = 1, itmax
        top = v(1)
        iq = 1
        do i = 2 , n
            if ( v(i) .gt. top ) then
                top = v(i)
                iq = i
            end if
        end do
        if ( top.lt.tol )  then ! Convergence reached: exits loop
            lsuccess = .true.
            exit
        endif
        ip = iv(iq)
        if ( ip.gt.iq) then
            i = ip
            ip = iq
            iq = i
        endif
        ap = a(ip,ip)
        aq = a(iq,iq)
        apq= a(ip,iq)
        al = 2.d0 * apq
        be = aq - ap
        a(ip,iq) = 0.d0
        a(iq,ip) = 0.d0
        if (abs(be) .gt. 1.d-10 .or. abs(apq) .gt. 1.d-10) then
            r = dsqrt ( al*al + be*be )
            c = dsqrt ( 0.5d0*( r + abs(be) ) / r )
            s = al / ( 2.d0 * r * c )
            if ( be.lt.0.d0 ) s = - s
            do j = 1, n
                aux    = x(j,ip)*c - x(j,iq)*s
                x(j,iq) = x(j,ip)*s + x(j,iq)*c
                x(j,ip) = aux
                if ( j.ne.ip .and. j.ne.iq ) then
                    aux    = a(j,ip)*c - a(j,iq)*s
                    a(j,iq) = a(j,ip)*s + a(j,iq)*c
                    a(j,ip) = aux
                    a(ip,j) = aux
                    a(iq,j) = a(j,iq)
                endif
            enddo
            a(ip,iq) = 0.d0
            a(iq,ip) = 0.d0
            aux = s * ( s*be - al*c )
            a(ip,ip) = ap + aux
            a(iq,iq) = aq - aux
        endif
        top = 0.d0
        do i = 1 , n
            if ( i.ne.ip ) then
                aux = abs( a(i,ip) )
                if ( aux .gt. top ) then
                    top = aux
                    iv(ip) = i
                end if
            end if
            v(ip) = top
        end do
        top = 0.d0
        do i = 1 , n
            if ( i.ne.iq ) then
                aux = abs( a(i,iq) )
                if ( aux .gt. top ) then
                    top = aux
                    iv(iq) = i
                end if
            end if
            v(iq) = top
        end do
    enddo
    if (.not. lsuccess) then
            write(6,"('Warning: maximum number of iterations in jacobi reached. Highest error = ', e12.5)") top
    endif
    do i = 1 , n
        v(i) = a(i,i)
    enddo
    if (nord.eq.1) return      ! Exits without ordering
    do i = 1,n
        do j = i+1,n
            if ( nord .eq. -1) then
                lordering = v(j).gt.v(i)
            else
                lordering = v(j).lt.v(i)
            endif
            if ( lordering ) then
                vv=v(j)
                v(j)=v(i)
                v(i)=vv
                do k = 1,n
                    vv=x(k,j)
                    x(k,j)=x(k,i)
                    x(k,i)=vv
                enddo
            endif
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
