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
!	Program for subtracting grid files (*.plt)
!	Input:
!		fname1, fname2: file names
!
!	Length units conversions from grid files to output file:
!
!	If grid files are in Bohr and output is wanted in Angstrom, 
!		run the program with the option "b2a" (without quotation marks) in the command line, i.e.
!			subtractplt.exe b2a
!	If length units in grid files are in Angstrom, and output is wanted in Angstrom
!		run the program with the option "a2a in the command line, i.e.
!			subtractplt.exe a2a
!	If length units in grid files are in Angstrom, and output is wanted in Bohr
!		run the program with the option "a2b in the command line or without options (default), i.e.
!			subtractplt.exe a2b
!			or	
!			subtractplt.exe
!	If grid files are in Bohr and output is wanted in Bohr too,
!		run the program with the option "b2b in the command line, i.e.
!			subtractplt.exe b2b
!	Version of September 2018
!
  program subtractplt320
    USE DAM320_D
    implicit none
    integer(KINT), parameter :: lengthfilename = 128, lengthfilename3 = 256
    real(KREAL4), allocatable :: v1(:), v2(:), v3(:)
    real(KREAL), allocatable :: v1d(:), v2d(:), v3d(:)
    real(KREAL) :: dif, difmax, dltx, dlty, dltz, errel, fconv, sigma2, umb, vmax, vmin
    real(KREAL) :: xinf, xsup, yinf, ysup, zinf, zsup, val1, val2, x, y, z, xmax, ymax, zmax
    integer(KINT) :: i, iaux(0:2), ierr, iuni1, iuni2, iuni3, ncifras(0:12), nerrabs(0:12), ix, ixmax, iy, iymax, iz, izmax
    integer(KINT) :: k, k2, npoints, nx1, nx2, ny1, ny2, nz1, nz2
    character(4) :: arg
    character(lengthfilename) :: fname1, fname2
    character(lengthfilename3) :: fname3
    real(KREAL), parameter :: ratob = 1.88972613288564d0, rbtoa = 5.2917720859d-1
    logical :: ldouble1, ldouble2, ldouble3
    arg = ""
    do i = 1, iargc()
    call getarg(i, arg)
 end do
    if (trim(arg) .eq. "a2a" .or. trim(arg) .eq. "b2b") then
        fconv = 1.d0
    elseif (trim(arg) .eq. "b2a") then
        fconv = rbtoa
    else
        fconv = ratob
    endif
    write(6,"('Name of first .plt file: ')")
    read(5,*) fname1
    k = min(len_trim(fname1), lengthfilename-4 )
    if (fname1(k-3:k) .ne. ".plt") fname1(k+1:k+4) = ".plt"
    write(6,"('Name of second .plt file: ')")
    read(5,*) fname2
    k = min(len_trim(fname2), lengthfilename-4 )
    if (fname2(k-3:k) .ne. ".plt") fname2(k+1:k+4) = ".plt"
    iuni1 = 10
    iuni2 = 11
    k = len_trim(fname1)
    k2 = len_trim(fname2)
    fname3 = fname1(1:k-4)//"-"//fname2(1:k2-4)//".plt"
#if _WIN32
    open (unit=iuni1, file=fname1, form='binary', carriagecontrol='NONE',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname1(1:len_trim(fname1)))//'. Stop')
    endif
    open (unit=iuni2, file=fname2, form='binary', carriagecontrol='NONE',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname2(1:len_trim(fname2)))//'. Stop')
    endif
    open (unit=iuni3, file=fname3, form='binary', carriagecontrol='NONE',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname3(1:len_trim(fname3)))//'. Stop')
    endif
#elif __INTEL_COMPILER
    open (unit=iuni1, file=fname1, form='binary', carriagecontrol='NONE',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname1(1:len_trim(fname1)))//'. Stop')
    endif
    open (unit=iuni2, file=fname2, form='binary', carriagecontrol='NONE',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname2(1:len_trim(fname2)))//'. Stop')
    endif
    open (unit=iuni3, file=fname3, form='binary', carriagecontrol='NONE',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname3(1:len_trim(fname3)))//'. Stop')
    endif
#else
    open (unit=iuni1, file=fname1, form='unformatted', access='stream',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname1(1:len_trim(fname1)))//'. Stop')
    endif
    open (unit=iuni2, file=fname2, form='unformatted', access='stream',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname2(1:len_trim(fname2)))//'. Stop')
    endif
    open (unit=iuni3, file=fname3, form='unformatted', access='stream',iostat=ierr)
    if (ierr .ne. 0) then
        call error(ierr,'Error opening file  '//trim(fname3(1:len_trim(fname3)))//'. Stop')
    endif
#endif
    write(6,"('Subtracting files ',a, ' from ', a)") trim(fname2), trim(fname1)
    read(iuni1) iaux(0),iaux(1)
    if (iaux(0) .eq. 0) then
        ldouble1 = .true.
    else
        ldouble1 = .false.
    endif
    read(iuni1) iaux(0),iaux(1),iaux(2)
    nz1 = iaux(0)
    ny1 = iaux(1)
    nx1 = iaux(2)
    read(iuni2) iaux(0),iaux(1)
    if (iaux(0) .eq. 0) then
        ldouble2 = .true.
    else
        ldouble2 = .false.
    endif

    if (ldouble1 .neqv. ldouble2) then
        write(6,"('Data in files have different precision.')")
        if (ldouble1) then
            write(6,"('data in file ', a, ' are in double precision'&
                    &' whereas data in file ', a, ' are in single precision')") trim(fname1), trim(fname2)
        else
            write(6,"('data in file ', a, ' are in single precision'&
                    &' whereas data in file ', a, ' are in double precision')") trim(fname1), trim(fname2)
        endif
        ldouble3 = .false.
        write(iuni3) 3,iaux(1)
    else
        ldouble3 = ldouble1
        if (ldouble3) then
            write(iuni3) 0,iaux(1)
        else
            write(iuni3) 3,iaux(1)
        endif
    endif
    read(iuni2) iaux(0),iaux(1),iaux(2)
    nz2 = iaux(0)
    ny2 = iaux(1)
    nx2 = iaux(2)
    if (nx1 .ne. nx2 .or. ny1 .ne. ny2 .or. nz1 .ne. nz2) then
        write(6,"('The dimensions of the grids are different')")
        write(6,"('Dimensions of grid in ', a, ' : nx = ', i3, ' ny = ', i3, ' nz = ', i3)") trim(fname1), nx1, ny1, nz1
        write(6,"('Dimensions of grid in ', a, ' : nx = ', i3, ' ny = ', i3, ' nz = ', i3)") trim(fname2), nx2, ny2, nz2
        write(6,"('Subtraction is not possible')")
        stop
    endif
    write(iuni3) iaux(0),iaux(1),iaux(2)
    allocate(v1(nx1), v1d(nx1), stat = ierr)
    if (.not. allocated(v1)) stop 'Memory error when allocating v1'
    allocate(v2(nx1), v2d(nx1), stat = ierr)
    if (.not. allocated(v2)) stop 'Memory error when allocating v2'
    allocate(v3(nx1), v3d(nx1), stat = ierr)
    if (.not. allocated(v3)) stop 'Memory error when allocating v3'
    if (ldouble1) then
        read(iuni1) (v1d(i), i = 1, 6)
        v1(1:6) = v1d(1:6)
    else
        read(iuni1) (v1(i), i = 1, 6)
    endif
    if (ldouble2) then
        read(iuni2) (v2d(i), i = 1, 6)
        v2(1:6) = v2d(1:6)
    else
        read(iuni2) (v2(i), i = 1, 6)
    endif
    umb = 1.d-5
    if (abs(v1(1)-v2(1)) .gt. umb .or. abs(v1(2)-v2(2)) .gt. umb .or. abs(v1(3)-v2(3)) .gt. umb &
            .or. abs(v1(4)-v2(4)) .gt. umb .or. abs(v1(5)-v2(5)) .gt. umb .or. abs(v1(6)-v2(6)) .gt. umb) then
        write(6,"('The boundaries of the grids are different')")
        write(6,"('Boundaries of grid in ', a, ' : xinf = ', f12.5, ' xsup = ', f12.5, ' yinf = ', f12.5 &
                , ' ysup = ', f12.5, ' zinf = ', f12.5, ' zsup = ', f12.5)") trim(fname1), v1(5), v1(6), v1(3), v1(4), v1(1), v1(2)
        write(6,"('Boundaries of grid in ', a, ' : xinf = ', f12.5, ' xsup = ', f12.5, ' yinf = ', f12.5 &
                , ' ysup = ', f12.5, ' zinf = ', f12.5, ' zsup = ', f12.5)") trim(fname2), v2(5), v2(6), v2(3), v2(4), v2(1), v2(2)
        write(6,"('Subtraction is not possible')")
        stop
    endif
    if (ldouble3) then
        write(iuni3) (v1d(i), i = 1, 6)
    else
        write(iuni3) (v1(i), i = 1, 6)
    endif
    xinf = v1(5) * fconv
    xsup = v1(6) * fconv
    yinf = v1(3) * fconv
    ysup = v1(4) * fconv
    zinf = v1(1) * fconv
    zsup = v1(2) * fconv
    dltx = (xsup-xinf) / (nx1-1)
    dlty = (ysup-yinf) / (ny1-1)
    dltz = (zsup-zinf) / (nz1-1)
    difmax = 0.d0
    sigma2 = 0.d0
    ixmax = 0
    iymax = 0
    izmax = 0
    val1 = 0.d0
    val2 = 0.d0
    do i = 0, 8
        ncifras(i) = 0
    enddo
    do i = 0, 10
        nerrabs(i) = 0
    enddo
    do iz = 1, nz1
        z = zinf + dltz * (iz-1)
        do iy = 1, ny1
            y = yinf + dlty * (iy-1)
            if (ldouble1) then
                read(iuni1,end=9997,err=9996) (v1d(i), i = 1, nx1)
            else
                read(iuni1,end=9997,err=9996) (v1(i), i = 1, nx1)
                v1d(1:nx1) = v1(1:nx1)
            endif
            if (ldouble2) then
                    read(iuni2,end=9999,err=9998) (v2d(i), i = 1, nx1)
            else
                read(iuni2,end=9999,err=9998) (v2(i), i = 1, nx1)
                v2d(1:nx1) = v2(1:nx1)
            endif
            do ix = 1, nx1
                x = xinf + dltx * (ix-1)
                dif = abs(v1d(ix)-v2d(ix))
                v3d(ix) = v1d(ix)-v2d(ix)
                if (.not. ldouble3) v3(ix) = v3d(ix)
                if (dif .gt. difmax) then
                    difmax = dif
                    ixmax = ix
                    iymax = iy
                    izmax = iz
                    xmax = x
                    ymax = y
                    zmax = z
                    val1 = v1d(ix)
                    val2 = v2d(ix)
                endif
                sigma2 = sigma2 + dif*dif
            enddo
            if (ldouble3) then
                write(iuni3) (v3d(i), i = 1, nx1)
            else
                write(iuni3) (v3(i), i = 1, nx1)
            endif
        enddo
    enddo
    npoints = nx1 * ny1 * nz1
    write(6,"(/10x,'Statistics',/10x,10('='),/)")
    write(6,"('Number of points examined = ', i10)") npoints
    if (trim(arg) .eq. "a2a" .or. trim(arg) .eq. "b2a") then
        write(6,"('Length units: Angstrom')")
    else
        write(6,"('Length units: Bohr')")
    endif
    write(6,"('Highest absolute difference = ', e12.5)") difmax
    if (ixmax+iymax+izmax .gt. 0) then
        write(6,"('Point with highest absolute difference = ', 3(1x,i3), 3(1x,e17.10))") ixmax, iymax, izmax, xmax, ymax, zmax
        write(6,"('Values in point with highest absolute difference = ', 3(1x,e22.15))") val1, val2
    endif
    write(6,"('Variance = ', e17.10)") sigma2 / npoints
    write(6,"('Standard deviation = ', e17.10)") sqrt(sigma2 / npoints)
    close(iuni3)
    close(iuni2)
    close(iuni1)
    stop
9996	continue
    write(6,"('Error reading data of file ', a,'. Stop')") trim(fname1)
    stop
9997	continue
    write(6,"('End of data of file ', a,'. Stop')") trim(fname1)
    stop
9998	continue
    write(6,"('Error reading data of file ', a,'. Stop')") trim(fname2)
    stop
9999	continue
    write(6,"('End of data of file ', a,'. Stop')") trim(fname2)
    stop
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
