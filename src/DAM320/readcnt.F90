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
!	Program for reading binary files .cnt  and printing the content as plain text
!	The program prompts for the name of the file to be read
!
! Version of September 2018
!
  program readcnt
    parameter (mxlin = 1000)
    real(4) u1, u2, v1, v2, v(0:mxlin)
    integer(4) iu, iv, aux(0:1), iuni
    character(60) :: fnameinp, fnameout
    write(6,"('Name of file to be read: ')")
    read(5,*) fnameinp
    fnameinp = trim(fnameinp)
    k = min(len_trim(fnameinp), 60-5 )
    if (fnameinp(k-3:k) .ne. ".cnt") fnameinp = trim(fnameinp)//".cnt"
! 	Opens input file (.cnt)
    iuni = 10
#if _WIN32
    open (unit=iuni, file=trim(fnameinp), form='binary', carriagecontrol='NONE',iostat=ioerr)
#elif __INTEL_COMPILER
    open (unit=iuni, file=trim(fnameinp), form='binary', carriagecontrol='NONE',iostat=ioerr)
#else
    open (unit=iuni, file=trim(fnameinp), form='unformatted', access='stream',iostat=ioerr)
#endif
    if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(fnameinp)//'.Stop')
! 	Opens output file .cnt_txt
    fnameout = trim(fnameinp)//"_txt"
    open (unit=7, file=trim(fnameout), form='formatted', status='unknown', iostat=ioerr)
    if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(fnameout)//'.Stop')

! 	Opens output file .cnt_gnu_txt
    fnameout = trim(fnameinp)//"_gnu_txt"
    open (unit=8, file=trim(fnameout), form='formatted', status='unknown', iostat=ioerr)
    if (ioerr .ne. 0) call error(1,'Error when opening file '//trim(fnameout)//'.Stop')


! 	Reads input data and writes output 
    read(iuni) aux(0),aux(1)
    write(7,"(3i10)") aux(0),aux(1)
    nu = aux(0)
    nv = aux(1)
    read(iuni) (v(i), i = 0, 3)
    write(7,"(3(1x,e15.8))") (v(i), i = 0, 3)
    u1 = v(0)
    u2 = v(1)
    v1 = v(2)
    v2 = v(3)
    do iv = 0, nv-1
        read(iuni,end=9999,err=9998) (v(iu), iu = 0, nu-1)
        write(7,"(//6(1x,e15.8))") (v(iu), iu = 0, nu-1)
        do iu = 0, nu-1
            write(8,"(4(1x,e15.8))") cero, v1 + (v2-v1) * iv / nv, u1 + (u2-u1) * iu / nu, v(iu)
        enddo
        write(8,*)
    enddo
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
