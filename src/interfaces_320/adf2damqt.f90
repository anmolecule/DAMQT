#include <config.h>
program adf2damqt

   use DAMQTinterfaceModule

   use Vartypes
   implicit none

   integer(KINT)           :: narg, i, iu, ls
   integer(KINT), external :: nargc

   character(D_SCM_MAXSTRINGLENGTH) :: name, option
   
   logical lspin, lorbitals

   call init01 ('Utils')
   call ppserial

!  input admits up to three optional arguments: 
!     one is the common root name for the files to be generated; must be the first one (default: ADF) 
!     the other two can be given in any order and are options:
!        SPIN  for spin density instead of electron density
!        NOORBITALS to discard generation of files with molecular orbitals
   narg = nargc()
   lspin = .false.    ! Default option: electron density 
   lorbitals = .true. ! Default option: generates files with molecular orbitals
   if (narg < 1) then
      name = 'ADF'
   else 
      call gtarg (1, name)
      write(option,"(a)") name
      call caseup(option)
      if(trim(option) == 'SPIN') then
         lspin = .true.
         name = 'ADF'
      else if (trim(option) == 'NOORBITALS') then
         lorbitals = .false.
         name = 'ADF'
      endif
      if (narg > 1) then
         do i = 2, narg
            call gtarg (i, option)
            call caseup(option)
            if (trim(option) == 'SPIN') then
               lspin = .true.
            else if (trim(option) == 'NOORBITALS') then
               lorbitals = .false.
            endif
         end do
      end if
   endif
   
   call dint21
   call rdtp21
   call DAM_interface(trim(name), lspin, lorbitals)

   call exit01 ('Utils', 'NORMAL TERMINATION')

end
