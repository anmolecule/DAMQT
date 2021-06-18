!  Copyright 2013-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
!  Guillermo Ramirez
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
!	Program for comparing grid files (*.plt)
!	Input:
!		fname1, fname2: file names
!		ndecs: number of coincident decimal figures for listing: data with less than ndecs coincident decimal figures
!				are printed to file compare.txt (if ndecs < 0: no listing)
!
!	Length units conversions from grid files to output file:
!
!	If grid files are in Bohr and output is wanted in Angstrom, 
!		run the program with the option "b2a" (without quotation marks) in the command line, i.e.
!			compare.exe b2a
!	If length units in grid files are in Angstrom, and output is wanted in Angstrom
!		run the program with the option "a2a in the command line, i.e.
!			compare.exe a2a
!	If length units in grid files are in Angstrom, and output is wanted in Bohr
!		run the program with the option "a2b in the command line or without options (default), i.e.
!			compare.exe a2b
!			or	
!			compare.exe
!	If grid files are in Bohr and output is wanted in Bohr too,
!		run the program with the option "b2b in the command line, i.e.
!			compare.exe b2b
!
! Version of September 2018
!
  program compareplt320
	USE DAM320_D
	implicit none
    integer(KINT), parameter :: lengthfilename = 256
	real(KREAL4), allocatable :: v1(:), v2(:)
	real(KREAL), allocatable :: v1d(:), v2d(:)
	real(KREAL) :: dif, difmax, dltx, dlty, dltz, errel, fconv, sigma2, umb, vmax, vmin
	real(KREAL) :: xinf, xsup, yinf, ysup, zinf, zsup, val1, val2, x, y, z, xmax, ymax, zmax
	integer(KINT) :: i, iaux(0:2), ierr, iuni1, iuni2, ix, ixmax, iy, iymax, iz, izmax, j, k
	integer(KINT) :: ncifras(0:12), ndecs, nerrabs(0:12), npoints, nx1, nx2, ny1, ny2, nz1, nz2
	character(4) :: arg
	character(60) :: fname1, fname2
	real(KREAL), parameter :: ratob = 1.88972613288564d0, rbtoa = 5.2917720859d-1
    logical :: ldouble1, ldouble2, lexist
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
    if (fname1(k-3:k) .ne. ".plt" .and. fname1(k-4:k) .ne. ".pltd") then
        fname1(k+1:k+4) = ".plt"
        inquire(file=trim(fname1), exist=lexist, iostat=ierr)
        if (.not. lexist) then
            fname1(k+1:k+5) = ".pltd"
        endif
    endif
    inquire(file=trim(fname1), exist=lexist, iostat=ierr)
    if (.not. lexist) then
        write(6,"('No grid file (plt or pltd) ', a, ' available. Stop')") fname1(1:k)
        stop
    endif
	write(6,"('Name of second .plt file: ')")
	read(5,*) fname2
	k = min(len_trim(fname2), lengthfilename-4 )
    if (fname2(k-3:k) .ne. ".plt" .and. fname2(k-4:k) .ne. ".pltd") then
        fname2(k+1:k+4) = ".plt"
        inquire(file=trim(fname2), exist=lexist, iostat=ierr)
        if (.not. lexist) then
            fname2(k+1:k+5) = ".pltd"
        endif
    endif
    inquire(file=trim(fname2), exist=lexist, iostat=ierr)
    if (.not. lexist) then
        write(6,"('No grid file (plt or pltd) ', a, ' available. Stop')") fname2(1:k)
        stop
    endif
	iuni1 = 10
	iuni2 = 11
#if _WIN32
	open (unit=iuni1, file=fname1, form='binary', carriagecontrol='NONE',iostat=ierr)
	if (ierr .ne. 0) then
		call error(ierr,'Error opening file  '//trim(fname1(1:len_trim(fname1)))//'. Stop')
	endif
	open (unit=iuni2, file=fname2, form='binary', carriagecontrol='NONE',iostat=ierr)
	if (ierr .ne. 0) then
		call error(ierr,'Error opening file  '//trim(fname2(1:len_trim(fname2)))//'. Stop')
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
#else
	open (unit=iuni1, file=fname1, form='unformatted', access='stream',iostat=ierr)
	if (ierr .ne. 0) then
		call error(ierr,'Error opening file  '//trim(fname1(1:len_trim(fname1)))//'. Stop')
	endif
	open (unit=iuni2, file=fname2, form='unformatted', access='stream',iostat=ierr)
	if (ierr .ne. 0) then
		call error(ierr,'Error opening file  '//trim(fname2(1:len_trim(fname2)))//'. Stop')
	endif
#endif
	write(6,"('Choose the lower number of decimal figures for listing (if < 0: none): ')",advance='no')
	read(5,*) ndecs
	if (ndecs .ge. 0) then
		open(unit = 7, file='compare.txt', form='formatted', iostat=ierr)
	endif
	if (ierr .ne. 0) then
		call error(ierr,'Error opening file compare.txt. Stop')
	endif
	
	write(6,"('comparing the files ',a, ' and ', a)") trim(fname1), trim(fname2) 
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
	endif
	read(iuni2) iaux(0),iaux(1),iaux(2)
	nz2 = iaux(0)
	ny2 = iaux(1)
	nx2 = iaux(2)
	if (nx1 .ne. nx2 .or. ny1 .ne. ny2 .or. nz1 .ne. nz2) then
		write(6,"('The dimensions of the grids are different')")
		write(6,"('Dimensions of grid in ', a, ' : nx = ', i3, ' ny = ', i3, ' nz = ', i3)") trim(fname1), nx1, ny1, nz1
		write(6,"('Dimensions of grid in ', a, ' : nx = ', i3, ' ny = ', i3, ' nz = ', i3)") trim(fname2), nx2, ny2, nz2
		write(6,"('No comparison is possible')")
		stop
	endif
	allocate(v1(nx1), v1d(nx1), stat = ierr)
	if (ierr .ne. 0) stop 'Memory error when allocating v1 and v1d'
	allocate(v2(nx1), v2d(nx1), stat = ierr)
	if (ierr .ne. 0) stop 'Memory error when allocating v2 and v2d'

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
		write(6,"('No comparison is possible')")
		stop
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
	if (ndecs .ge. 0) then
		if (trim(arg) .eq. "a2a" .or. trim(arg) .eq. "b2a") then
			write(7,"('Length units: Angstrom')")
		else
			write(7,"('Length units: Bohr')")
		endif
		write(7,"('xinf = ', e17.10, ' xsup = ', e17.10)") xinf, xsup
		write(7,"('yinf = ', e17.10, ' ysup = ', e17.10)") yinf, ysup
		write(7,"('zinf = ', e17.10, ' zsup = ', e17.10)") zinf, zsup
		write(7,"('dltx = ', e17.10, ' dlty = ', e17.10, ' dltz = ', e17.10)") dltx, dlty, dltz
	endif
	difmax = 0.d0
	sigma2 = 0.d0
	ixmax = 0
	iymax = 0
	izmax = 0
	val1 = 0.d0
	val2 = 0.d0
	ncifras = 0
	nerrabs = 0
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
				vmin = min(abs(v1d(ix)),abs(v2d(ix)))
				vmax = max(abs(v1d(ix)),abs(v2d(ix)))
				if (vmax .eq. 0.d0) then
					errel = 0.d0
				else
					errel = 1.d0 - vmin/vmax
				endif
				if (dif .le. 1.d-12) then
					nerrabs(12) = nerrabs(12)+1
				else
					j = int(-log10(dif))
					if (j .lt. 0) then
						nerrabs(0) = nerrabs(0) + 1
					else
						nerrabs(j) = nerrabs(j) + 1
					endif
				endif
				if (ndecs .ge. 0 .and. dif .gt. 1.d-10 .and. j .le. ndecs) write(7,"('Point: ', 3(1x,i3), 3(1x,e17.10) &
					,' values = ', 2(1x,e20.13), ' abs. err. = ', e12.5, ' j = ', i3)") &
					ix, iy, iz, x, y, z, v1d(ix),  v2d(ix), dif, j
				if (errel .le. 1.d-12) then
					ncifras(12) = ncifras(12) + 1
				else
					j = int(-log10(errel))
					if (j .lt. 0) then
						ncifras(0) = ncifras(0) + 1
					else
						ncifras(j) = ncifras(j) + 1
					endif
				endif	
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
		enddo
	enddo
	npoints = nx1 * ny1 * nz1
	write(6,"(/10x,'Statistics',/10x,10('='),/)")
	write(7,"(/10x,10('='),/10x,'Statistics',/10x,10('='),/)")
	write(6,"('Number of points examined = ', i10)") npoints
	write(7,"('Number of points examined = ', i10)") npoints
	if (trim(arg) .eq. "a2a" .or. trim(arg) .eq. "b2a") then
		write(6,"('Length units: Angstrom')")
		write(7,"('Length units: Angstrom')")
	else
		write(6,"('Length units: Bohr')")
		write(7,"('Length units: Bohr')")
	endif
	write(6,"('Highest absolute error = ', e12.5)") difmax
	write(7,"('Highest absolute error = ', e12.5)") difmax
	if (ixmax+iymax+izmax .gt. 0) then
		write(6,"('Point with highest absolute error = ', 3(1x,i3), 3(1x,e17.10))") ixmax, iymax, izmax, xmax, ymax, zmax
		write(6,"('Values in point with highest absolute error = ', 3(1x,e22.15))") val1, val2
		write(7,"('Point with highest absolute error = ', 3(1x,i3), 3(1x,e17.10))") ixmax, iymax, izmax, xmax, ymax, zmax
		write(7,"('Values in point with highest absolute error = ', 3(1x,e22.15))") val1, val2
	endif
	write(6,"('Variance = ', e17.10)") sigma2 / npoints
	write(6,"('Standard deviation = ', e17.10)") sqrt(sigma2 / npoints)
	write(6,"('Accuracy: ',/ ('Number of coincident figures = ', i2, ' number of points = ', i10))") &
		(i, ncifras(i), i = 0, 12) 
	write(6,"(/'Precision: ',/ ('Number of coincident decimals = ', i2, ' number of points = ', i10))") &
		(i, nerrabs(i), i = 0, 12)
	write(7,"('Variance = ', e17.10)") sigma2 / npoints
	write(7,"('Standard deviation = ', e17.10)") sqrt(sigma2 / npoints)
	write(7,"('Accuracy: ',/ ('Number of coincident figures = ', i2, ' number of points = ', i10))") &
		(i, ncifras(i), i = 0, 12) 
	write(7,"(/'Precision: ',/ ('Number of coincident decimals = ', i2, ' number of points = ', i10))") &
		(i, nerrabs(i), i = 0, 12)
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
	character(*) :: msg
	write(6,"(a)") msg
	write(6,"('Error code = ', i4)") ierr
	stop
	end
