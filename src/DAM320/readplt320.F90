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
!	Program for reading binary files .plt  and printing the content as plain text
!	The program prompts for the name of the file to be read
!
! Version of September 2018
!
  program readplt320
    USE DAM320_D
    implicit none
    integer(KINT), parameter :: mxlin = 1000, lengthfilename = 128
    real(KREAL4) :: x1, x2, y1, y2, z1, z2, v(0:mxlin)
    real(KREAL) :: x1d, x2d, y1d, y2d, z1d, z2d, vd(0:mxlin)
    integer(KINT) :: aux(0:2), i, ierr, iuni, ix, iy, iz, k, nx, ny, nz
    character(lengthfilename) :: fnameinp, fnameout
    logical :: existe, ldouble

    write(6,"('Name of file to be read: ')")
    read(5,*) fnameinp
    fnameinp = trim(fnameinp)
    k = min(len_trim(fnameinp), lengthfilename-5 )
    if (fnameinp(k-3:k) .ne. ".plt" .and. fnameinp(k-4:k) .ne. ".pltd" ) fnameinp = trim(fnameinp)//".plt"
    existe = .false.
    inquire(file=fnameinp, exist=existe)
    if (.not. existe) then
        fnameinp = trim(fnameinp)//"d"
        inquire(file=fnameinp, exist=existe)
        if (.not. existe) then
            call error(1,'Files '//trim(fnameinp(1:len_trim(fnameinp)-4))//'+ .plt or .pltd do not exist. Stop')
        endif
    endif

! 	Opens input file (.plt or .pltd)
    iuni = 10
#if _WIN32
    open (unit=iuni, file=trim(fnameinp), form='binary', carriagecontrol='NONE',iostat=ierr)
#elif __INTEL_COMPILER
    open (unit=iuni, file=trim(fnameinp), form='binary', carriagecontrol='NONE',iostat=ierr)
#else
    open (unit=iuni, file=trim(fnameinp), form='unformatted', access='stream',iostat=ierr)
#endif
    if (ierr .ne. 0) call error(1,'Error when opening file '//trim(fnameinp)//'.Stop')
! 	Opens output file (.plt_txt)
    fnameout = trim(fnameinp)//"_txt"
    open (unit=7, file=trim(fnameout), form='formatted', status='unknown', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Error when opening file '//trim(fnameout)//'.Stop')
! 	Reads input data and writes output 
    ldouble = .false.
    read(iuni) aux(0),aux(1)
    if (aux(0) .eq. 0) ldouble = .true.
    write(7,"(2i10)") aux(0),aux(1)
    read(iuni) aux(0),aux(1),aux(2)
    write(7,"(3i10)") aux(0),aux(1), aux(2)
    nz = aux(0)
    ny = aux(1)
    nx = aux(2)

    if (ldouble) then
        read(iuni) (vd(i), i = 0, 5)
        write(7,"(5(1x,e22.15))") (vd(i), i = 0, 5)
        do iz = 1, nz
            do iy = 1, ny
                read(iuni,end=9999,err=9998) (vd(i), i = 0, nx-1)
                write(7,"(/6(1x,e22.15))") (vd(i), i = 0, nx-1)
            enddo
        enddo
    else
        read(iuni) (v(i), i = 0, 5)
        write(7,"(5(1x,e15.8))") (v(i), i = 0, 5)
        do iz = 1, nz
            do iy = 1, ny
                read(iuni,end=9999,err=9998) (v(i), i = 0, nx-1)
                write(7,"(/6(1x,e15.8))") (v(i), i = 0, nx-1)
            enddo
        enddo
    endif

    stop
9998	continue
    write(6,"('Error reading data. Stop')")
    stop
9999	continue
    write(6,"('End of data file. Stop')")
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
