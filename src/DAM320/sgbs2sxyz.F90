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
!  Program for generating a .xyz file with geometry in bohr from a .sgbs file 
!
! Version of September 2018
!
  program sgbs2sxyz
    USE DAM320_D
    USE DAM320_DATA_D
    USE DAM320_CONST_D
    implicit none
    integer(KINT) :: i, idummy, idummy2, idummy3, ierr, j, k
    logical lsgbs
    real(KREAL) :: dummy, dummy2, dummy3, x, y, z, znuc
    read(5,*) projectname
    open(17,file=trim(projectname)//".sgbs2sxyz",form='formatted', iostat=ierr)
    if (ierr .ne. 0) call error(1,'Error when opening file '//trim(projectname)//".sgbs2sxyz"//'.Stop')
    inquire(file=trim(projectname)//".sgbs", exist=lsgbs, iostat=ierr)
    if (ierr .eq. 0 .and. lsgbs) then
        open(15,file=trim(projectname)//".sgbs",form='unformatted',iostat=ierr)
        if (ierr .ne. 0) call error(1,'Error when opening file '//trim(projectname)//".sgbs"//'.Stop')
        rewind(15)
        open(16,file=trim(projectname)//".sxyz",form='formatted',iostat=ierr)
        if (ierr .ne. 0) call error(1,'Error when opening file '//trim(projectname)//".sxyz"//'.Stop')
        read(15) dummy		!	reads umbrzn to dummy
        read(15) idummy, idummy2, idummy3	!	reads m2cmxcap, m2cmxfun, m2cmxcen to dummy
        read(15) ncen
        write(16,"(i7)") ncen
        read(15) idummy, idummy	!	reads nbas, ncaps to dummy
        read(15) dummy	!	reads repnuc to dummy
        do j = 1, ncen
            read(15) idummy, idummy2
            if (idummy .gt. 0) then
                do k = idummy, idummy2
                    read(15) idummy3, idummy3, idummy3, dummy, dummy  !	reads nf(i), nn(i), ll(i), xx(i), rnor(i) to dummy
                enddo
            endif
        enddo
        do i = 1, 3
            read(15) dummy, dummy, dummy
        enddo
        do i = 1, ncen
            read(15) x, y, z, znuc
            write(16,"(4(3x,e17.10))") x, y, z, znuc
        enddo
        write(17,"('Geometry written to file ', a)") trim(projectname)//".sxyz"
        close(17)
        close(16)
        close(15)
    else
        inquire(file=trim(projectname)//".sgbsden", exist=lsgbs, iostat=ierr)
        if (ierr .eq. 0 .and. lsgbs) then
            open(15,file=trim(projectname)//".sgbsden",form='formatted',iostat=ierr)
            if (ierr .ne. 0) call error(1,'Error when opening file '//trim(projectname)//".sgbsden"//'.Stop')
            rewind(15)
            open(16,file=trim(projectname)//".sxyz",form='formatted',iostat=ierr)
            if (ierr .ne. 0) call error(1,'Error when opening file '//trim(projectname)//".sxyz"//'.Stop')
            read(15,*) ncen
            write(16,"(i7)") ncen
            do i = 1, ncen
                read(15,*) x, y, z, znuc
                write(16,"(4(3x,e17.10))") x, y, z, znuc
            enddo
            write(17,"('Geometry written to file ', a)") trim(projectname)//".sxyz"
            close(17)
            close(16)
            close(15)
        else
            write(17,"('Files ',a , ' and ', a, ' do not exist')") trim(projectname)//".sgbs", trim(projectname)//".sgbsden"
            close(17)
        endif
    endif
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
