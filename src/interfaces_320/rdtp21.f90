subroutine rdtp21 

!  ==================================================================
!  purpose:  READ miscellaneous data from TAPE21, not read by dint21.
!            and set spinsuffix
!  ==================================================================

   use KF
   use DimensionsModule
   use SpinsuffixModule
   use GeometrydataModule

   use Vartypes
   implicit none

   integer(KINT) :: iu

!  ===========
!  open tape21
!  ===========

   call kfopfl (iu, 'TAPE21')

!  ----------------
!  section Geometry 
!  ----------------

   call kfopsc (iu, 'Geometry')
   call kfread (iu, 'nr of atoms', gGeometrydata%natoms)

   call kfclfl (iu)

   gSpinsuffix%sspin(1) = '_A'
   gSpinsuffix%sspin(2) = '_B'

!  =====================================================================
end subroutine rdtp21

