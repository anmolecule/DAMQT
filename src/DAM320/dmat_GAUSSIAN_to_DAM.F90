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
!	Program for reordering the density matrix from GAUSSIAN order to DAM order
!	input:
! 		nshells: number of shells in the basis set (each shell corresponds to a set of basis functions differing only
!				in the "m" quantum number (i.e., l = 1 implies three P functions; l = 2, five D functions; and so forth
!		ll: l quantum number of each shell (in the same order as occurring in the basis set)
!		dmatGAUSS: lower triangle of the density matrix as given by GAUSSIAN
!	output:
!		dmatDAM:	lower triangle of density matrix in DAM order (canonical order of spherical hemonics)
!
!	August 2014. Rafael Lopez (rafael.lopez@uam.es)
!
program dmat_GAUSSIAN_to_DAM
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    real(KREAL), allocatable :: dmatGAUSS(:,:), dmatDAM(:,:), daux(:,:)
    integer(KINT) :: i, ierr, ii, irow, ishell, j, jcol, jj, jshell, nbasis
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

    call readggbs
write(6,*) 'nbas = ', nbas

    allocate(dmatGAUSS(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating dmatGAUSS. Stop')
write(6,*) 'allocated(dmatGAUSS) =  ', allocated(dmatGAUSS)
    allocate(dmatDAM(nbas,nbas), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating dmatDAM. Stop')

    open(16,file=trim(projectname)//".den",form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Cannot open file '//trim(projectname)//".den"//'. Stop')
    endif

    read(16,*,iostat = ierr) nbasis
    if ( ierr .ne. 0 .or. nbas .ne. nbasis ) then
        write(6,"('ERROR reading density matrix. Check whether the density matrix correspond to this basis set.')")
        call error(1,' Stop')
    endif
write(6,*) 'nbasis = ', nbasis
nbasis = 1

    read(16,*,iostat = ierr) ((dmatGAUSS(i,j), j = 1, i), i = 1, nbas)
    if (ierr .ne. 0) call error(ierr,'Error reading dmatGAUSS. Stop')

    irow = 1
    do ishell = 1, ncaps
        jcol = 1
        do jshell = 1, ishell
            if (ishell .eq. jshell .and. ll(ishell) .ne. 0) then
                do i = irow, irow + 2*ll(ishell)
                    do j = i+1, irow + 2*ll(ishell)
                        dmatGAUSS(i,j) = dmatGAUSS(j,i)
                    enddo
                enddo
            endif
            if (ll(ishell) .eq. 0) then
                if (ll(jshell) .eq. 0) then
                    dmatDAM(irow,jcol) = dmatGAUSS(irow,jcol)
                elseif (ll(jshell) .eq. 1) then
                    dmatDAM(irow,jcol  ) = dmatGAUSS(irow,jcol+1)
                    dmatDAM(irow,jcol+1) = dmatGAUSS(irow,jcol+2)
                    dmatDAM(irow,jcol+2) = dmatGAUSS(irow,jcol)
                else
                    j = 0
                    do jj = 2*ll(jshell), 0, -2
                        dmatDAM(irow,jcol+j) = dmatGAUSS(irow,jcol+jj)
                        j = j+1
                    enddo
                    do jj = 1, 2*ll(jshell), 2
                        dmatDAM(irow,jcol+j) = dmatGAUSS(irow,jcol+jj)
                        j = j+1
                    enddo
                endif
            elseif (ll(ishell) .eq. 1) then
                if (ll(jshell) .eq. 0) then
                    dmatDAM(irow,jcol)   = dmatGAUSS(irow+1,jcol)
                    dmatDAM(irow+1,jcol) = dmatGAUSS(irow+2,jcol)
                    dmatDAM(irow+2,jcol) = dmatGAUSS(irow  ,jcol)
                elseif (ll(jshell) .eq. 1) then
                    dmatDAM(irow,jcol)     = dmatGAUSS(irow+1,jcol+1)
                    dmatDAM(irow,jcol+1)   = dmatGAUSS(irow+1,jcol+2)
                    dmatDAM(irow,jcol+2)   = dmatGAUSS(irow+1,jcol)
                    dmatDAM(irow+1,jcol)   = dmatGAUSS(irow+2,jcol+1)
                    dmatDAM(irow+1,jcol+1) = dmatGAUSS(irow+2,jcol+2)
                    dmatDAM(irow+1,jcol+2) = dmatGAUSS(irow+2,jcol)
                    dmatDAM(irow+2,jcol)   = dmatGAUSS(irow,jcol+1)
                    dmatDAM(irow+2,jcol+1) = dmatGAUSS(irow,jcol+2)
                    dmatDAM(irow+2,jcol+2) = dmatGAUSS(irow,jcol)
                else
                    j = 0
                    do jj = 2*ll(jshell), 0, -2
                        dmatDAM(irow,jcol+j)   = dmatGAUSS(irow+1,jcol+jj)
                        dmatDAM(irow+1,jcol+j) = dmatGAUSS(irow+2,jcol+jj)
                        dmatDAM(irow+2,jcol+j) = dmatGAUSS(irow,jcol+jj)
                        j = j+1
                    enddo
                    do jj = 1, 2*ll(jshell), 2
                        dmatDAM(irow,jcol+j)   = dmatGAUSS(irow+1,jcol+jj)
                        dmatDAM(irow+1,jcol+j) = dmatGAUSS(irow+2,jcol+jj)
                        dmatDAM(irow+2,jcol+j) = dmatGAUSS(irow,jcol+jj)
                        j = j+1
                    enddo
                endif
            else
                if (ll(jshell) .eq. 0) then
                    i = 0
                    do ii = 2*ll(ishell), 0, -2
                        dmatDAM(irow+i,jcol) = dmatGAUSS(irow+ii,jcol)
                        i = i+1
                    enddo
                    do ii = 1, 2*ll(ishell), 2
                        dmatDAM(irow+i,jcol) = dmatGAUSS(irow+ii,jcol)
                        i = i+1
                    enddo
                elseif (ll(jshell) .eq. 1) then
                    i = 0
                    do ii = 2*ll(ishell), 0, -2
                        dmatDAM(irow+i,jcol)   = dmatGAUSS(irow+ii,jcol+1)
                        dmatDAM(irow+i,jcol+1) = dmatGAUSS(irow+ii,jcol+2)
                        dmatDAM(irow+i,jcol+2) = dmatGAUSS(irow+ii,jcol)
                        i = i+1
                    enddo
                    do ii = 1, 2*ll(ishell), 2
                        dmatDAM(irow+i,jcol)   = dmatGAUSS(irow+ii,jcol+1)
                        dmatDAM(irow+i,jcol+1) = dmatGAUSS(irow+ii,jcol+2)
                        dmatDAM(irow+i,jcol+2) = dmatGAUSS(irow+ii,jcol)
                        i = i+1
                    enddo
                else
                    i = 0
                    do ii = 2*ll(ishell), 0, -2
                        j = 0
                        do jj = 2*ll(jshell), 0, -2
                            dmatDAM(irow+i,jcol+j)   = dmatGAUSS(irow+ii,jcol+jj)
                            j = j+1
                        enddo
                        do jj = 1, 2*ll(jshell), 2
                            dmatDAM(irow+i,jcol+j)   = dmatGAUSS(irow+ii,jcol+jj)
                            j = j+1
                        enddo
                        i = i + 1
                    enddo
                    do ii = 1, 2*ll(ishell), 2
                        j = 0
                        do jj = 2*ll(jshell), 0, -2
                            dmatDAM(irow+i,jcol+j)   = dmatGAUSS(irow+ii,jcol+jj)
                            j = j+1
                        enddo
                        do jj = 1, 2*ll(jshell), 2
                            dmatDAM(irow+i,jcol+j)   = dmatGAUSS(irow+ii,jcol+jj)
                            j = j+1
                        enddo
                        i = i + 1
                    enddo
                endif
            endif
            jcol = jcol + 2 * ll(jshell) + 1
        enddo
        irow = irow + 2 * ll(ishell) + 1

    enddo
    close(16)
    open(16,file=trim(projectname)//".den_DAM",form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Cannot open file '//trim(projectname)//".den_DAM"//'. Stop')
    endif

    write(16,*) nbas
    do i = 1, nbas
        write(16,"(7(1x,e17.10))") (dmatDAM(i,j), j = 1, i)
    enddo
    close(16)

    stop
    end
!
!	***************************************************************
!
  subroutine readggbs
    USE DAM320_D
    USE DAM320_CONST_D
    USE DAM320_DATA_D
    USE GAUSS
    implicit none
    integer(KINT) :: i, ia, icarga, ierr, indnf, indng, ios, j, k, k1, k2, knt, ncapserr
    real(KREAL) :: aux, bux
    real(KREAL) :: xaux(mxprimit), cfaux(mxprimit)

    open(15,file=trim(projectname)//".ggbs",form='formatted', iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Cannot open file '//trim(projectname)//".ggbs"//'. Stop')
    endif

!	Reads the number of centers

    read(15,*) ncen

!	Allocates memory for geometry and basis set
    ncaps = mxcap ! just for allocating

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating atmnam. Stop')

!    allocate(cfcontr0(mxcap*mxprimit))
!    if (.not. allocated(cfcontr0)) call error(1,'Memory error when allocating cfcontr0. Stop')

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

!    allocate(xxg0(mxcap*mxprimit))
!    if (.not. allocated(xxg0)) call error(1,'Memory error when allocating xxg0. Stop')

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
            read(15,*) nprimit(ncaps), ll(ncaps)
            if (ncaps .gt. mxcap)  call error(1,'Error: maximum number of shells in basis set exceeded. Stop')
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
!!			sorts primitives in increasing exponents
!            call sort(nprimit(ncaps),xaux)
!            do k = 1, nprimit(ncaps)
!                xxg0(icarga+k) = xaux(k)
!                cfcontr0(icarga+k) = cfaux(isort(k))
!            enddo
!!			computes and stores the radial normalization factor
!            aux = cero
!            bux = ll(ncaps) + 1.5d0
!            do k1 = 1, nprimit(ncaps)
!                do k2 = 1, k1-1
!                        aux=aux + dos*cfcontr0(icarga+k1)*cfcontr0(icarga+k2)/(xaux(k1)+xaux(k2))**bux
!                enddo
!                aux = aux + cfcontr0(icarga+k1) * cfcontr0(icarga+k1) / (dos*xaux(k1))**bux
!            enddo
!            rnor(ncaps) = sqrt( dos / (facts(ll(ncaps))*aux) )
!!			actualizes the index for loading
            icarga = icarga+nprimit(ncaps)
        enddo
    enddo
    nprimitot = icarga
    nbas = indnf-1

!    allocate(xxg(nprimitot))
!    if (.not. allocated(xxg)) call error(1,'Memory error when allocating xxg. Stop')

!    allocate(cfcontr(nprimitot))
!    if (.not. allocated(cfcontr)) call error(1,'Memory error when allocating cfcontr. Stop')

!    xxg(1:nprimitot) = xxg0(1:nprimitot)
!    cfcontr(1:nprimitot) = cfcontr0(1:nprimitot)

!    deallocate(xxg0, cfcontr0)

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

!!	prints out the input data
!    write(6,"(24x,'GEOMETRY (BOHR)')")
!    write(6,"(/t1, ' no. of center:', t20, 'x', t32, 'y', t44, 'z', t56, 'charge', t68, 'n. of shells')")
!    do ia = 1, ncen
!        write(6,"(t4, i5, t13, f12.7, t25, f12.7, t37, f12.7, t51, f10.5, t73, i3)") &
!                ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia), ngfin(ia)-ngini(ia)+1
!    enddo
!    write(6,"(27x,'GTO BASIS SET')")
!    if (longoutput) then
!        icarga = 0
!        knt = 0
!        do ia = 1, ncen
!            if (ncontr(ia) .le. 0) cycle
!            write(6,"(/1x,'atom no.',1x,i4,'(',a2,')')") ia, atmnam(ia)
!            write(6,"(1x,'number of contractions = ',i4)") ncontr(ia)
!            do j = 1, ncontr(ia)
!                knt = knt + 1
!                write(6,"(/1x,'contraction no. ',i4,' ; l = ',i2)") j,  ll(knt)
!                write(6,"('exponents: ', 8(1x,e12.5))") (xxg(icarga+k), k = 1, nprimit(knt))
!                write(6,"('coefficients: ', 8(1x,e12.5))") (cfcontr(icarga+k),k=1,nprimit(knt))
!                icarga = icarga+nprimit(knt)
!            enddo
!        enddo
!    endif
    write(6,"('Number of basis functions = ', i4)") nbas
    write(6,"('Number of primitive functions = ', i4)") nprimitot
    close(15)
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
  subroutine error(ierr, msg)
    integer(4) :: ierr
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    stop
    end
