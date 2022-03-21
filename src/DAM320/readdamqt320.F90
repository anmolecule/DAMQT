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
!	Program for reading files _2016.damqt
!
!	Version of September 2018
!
  program readdamqt2016
    use DAM320_D
    implicit double precision (a-h,o-z)
    parameter (lengthfilename = 128)
    logical :: lSTO
    integer(KINT), allocatable :: icfpos(:), ncontr(:)
    real(KREAL), allocatable :: cfajust(:), cfcontr(:), xajust(:), xxg(:)
    character(lengthfilename) :: fnameinp, fnameout
    if (iargc() .gt. 0) then
        call getarg(1, fnameinp)
    else
        write(6,"('Name of file to be read: ')")
        read(5,*) fnameinp
    endif
! 	Opens input file (.damqt)
    fnameinp = trim(fnameinp)
    k = min(len_trim(fnameinp), lengthfilename-6 )
    if (fnameinp(k-10:k) .ne. "_2016.damqt") fnameinp = trim(fnameinp)//"_2016.damqt"
#if _WIN32
    open (unit=10, file=fnameinp, form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#elif __INTEL_COMPILER
    open (unit=10, file=fnameinp, form='binary', action = 'read', carriagecontrol='NONE', iostat=ierr)
#else
    open (unit=10, file=fnameinp, form='unformatted', action = 'read', access='stream', iostat=ierr)
#endif
    if (ierr .ne. 0) call error(ierr,'Error when opening file '//fnameinp)
! 	Opens output file (.damqt_txt)
    fnameout = trim(fnameinp)//"_txt"
    open (unit=7, file=trim(fnameout), form='formatted', status='unknown', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Error when opening file '//trim(fnameout)//'.Stop')

    read(10) ncen, nbas, ncaps
    write(7,"('ncen = ', i3, ' nbas = ', i5, ' ncaps = ', i3)") ncen, nbas, ncaps
    call flush(7)

!	Geometry and nuclear charges
    do i = 1, ncen
        read(10) xcen, ycen, zcen, zn
        write(7,"('centro ',i4, ' Zn = ', f8.1,4x,'coords: ', 3(1x,e17.10))") i, zn, xcen, ycen, zcen
        call flush(7)
    enddo

!	Basis set
    i = 0
    read(10) lSTO
write(7,*) 'lSTO = ', lSTO
    if (lSTO) then
        write(7,"(/t22,'STO Basis set',/t22,13('-'))")
        do j = 1, ncen
            read(10) ngini, ngfin
            write(7,"(t5,'centro ', i4,/t12,'n', t16, 'l', t25,'exp', t35, 'ind func')") j
write(7,*) 'ngini = ', ngini, ' ngfin = ', ngfin
            do k = ngini, ngfin
                i = i + 1
                read(10) nf, nn, ll, xx
                write(7,"(t11,i2,t15,i2,t20,e12.5,t36,i4)") nn, ll, xx, nf
                call flush(7)
            enddo
        enddo
    else
        write(7,"(/t22,'GTO Basis set',/t22,13('-'))")
        read(10) nprimitot
        write(7,"('Total number of primitives = ', i7)") nprimitot
        allocate(cfcontr(nprimitot), xxg(nprimitot), ncontr(ncen), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating cfcontr, xxg and ncontr. Stop')
        do i = 1, ncen
            read(10) ncontr(i)
            write(7,"(/'Number of contractions in center ', i4,': ', i3)") i, ncontr(i)
            do j = 1, ncontr(i)
                read(10) nprimit, l
                write(7,"('nprimit = ', i3, ' l = ', i3)") nprimit, l
                read(10) xxg(1:nprimit)
                write(7,"('   exponents    = ', 5(1x,e17.10))") xxg(1:nprimit)
                read(10) cfcontr(1:nprimit)
                write(7,"('   coefficients = ', 5(1x,e17.10))") cfcontr(1:nprimit)
            enddo
        enddo
        deallocate(cfcontr, xxg)
    endif
!	Data of density representation
    read(10) lmaxexp
    lmtop = (lmaxexp+1)*(lmaxexp+1)
    write(7,"('lmaxexp = ', i2)") lmaxexp
    allocate(icfpos(lmtop*nintervaj+1), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating icfpos. Stop')

    allocate(xajust(nintervaj), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating xajust. Stop')
    do ia = 1, ncen      ! Do over centers
        if (.not. lSTO .and. ncontr(ia) .lt. 1) cycle
        read(10) icfpos(1:lmtop*nintervaj+1)
        write(7,"('icfpos(',i6,') = ',/ 15(1x,i9))") ia, icfpos(1:lmtop*nintervaj+1)
        read(10) xajust(1:nintervaj)
        write(7,"('xajust = ', 8(1x,e17.10))") xajust(1:nintervaj)
        allocate(cfajust(icfpos(lmtop*nintervaj+1)-1), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error when allocating cfajust. Stop')
        read(10) cfajust(1:icfpos(lmtop*nintervaj+1)-1)
        write(7,"('cfajust = ', 8(1x,e17.10))") cfajust(1:icfpos(lmtop*nintervaj+1)-1)
        deallocate(cfajust)
    enddo
    close(10)
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
