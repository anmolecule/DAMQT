!  Copyright 2011-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
! Program for the tabulation of molecular density expanded as a sum of atomic densities
!       expressed as a linear combination of Gaussian radial factors times unnormalized real spherical harmonics
!       or from its expansion in Canterakis-Zernike or Jacobi expansion  
!
! Files with Canterakis-Zernike (.zernike) or Jacobi (.jacobi) expansions must be available. These files
!       can be generated with the programs: DAMZernike-Jacobi_GTO, DAMZernike-Jacobi_STO,
!       DAMZernike-Jacobi_GTO_mpi or DAMZernike-Jacobi_STO_mpi
!
! Version of May 2018
!
! #define DBLPRCGRID    ! Uncomment this line  if double precision grid is wanted
MODULE DAMDENZJ_atdens
    IMPLICIT NONE
        integer(4), parameter :: KINT = 4, KREAL = 8, KREAL4 = 4, KINT8 = 8
! 		Atomic symbols
    character(2) :: atmnms(0:103) = (/							 &
            "q ", "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F "   &
            , "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K " &
            , "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu" &
            , "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y " &
            , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In" &
            , "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr" &
            , "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm" &
            , "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au" &
            , "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac" &
            , "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es" &
            , "Fm", "Md", "No", "Lw" /)
    character dirsep	! character for directory names separator: "/" for unix, "\\" for MS-windows
    character(2), allocatable :: atmnam(:)
    character(300) :: projectname
    logical :: lechelon, lexact, ljacobi, lmthfile, lrstarrel
    logical iswindows	! .true. if running on a MS-windows system
    logical :: longoutput, lgbsgz
    integer(KINT), parameter ::  mxcen = 500, mxl = 6, mxcap = 8000
    integer(KINT), parameter ::  mxldst = mxl+mxl
    integer(KINT), parameter :: mxk = 350, mxkextra = 100, mxl4 = 4*mxl
    integer(KINT), parameter :: npol = 10
    integer(KINT), parameter :: mxroot = 800, mxreal = 3000, mxfact = 150
    integer(KINT), allocatable :: ll(:), lmaxc(:), ngini(:), ngfin(:), nf(:), nzn(:), ind(:)
    integer(KINT) :: ldimaux, lexpansion, nquadpoints, kexpansion
    integer(KINT) :: lmaxbase, mxltot, mxind, ncaps, ncen
    integer(KINT) :: mxemes, mxkcof, mxlcof
    real(KREAL), parameter :: toldstorig = 1.d-12
    real(KREAL), parameter :: cero = 0.d0, uno = 1.d0, dos = 2.d0, cuatro = 4.d0, ocho = 8.d0, r16 = 16.d0, udec = 0.1d0
    real(KREAL), parameter :: umed = .5d0, pt25 = .25d0, euler = 5.7721566490153280D-1, raiz2 = 1.4142135623730950d0
    real(KREAL), parameter :: alfacutoff = 1.d-20, umbrzn = 1.d-6
    real(KREAL), allocatable  :: akgkl(:), cfgkl1(:), cfgkl2(:), cfgkl3(:), gkl(:)
    real(KREAL), allocatable :: omeganlm(:,:), flm(:,:,:), flmmamb(:,:), ftot(:,:), rquadaux(:), rquad01(:), rquadscal(:)
    real(KREAL), allocatable :: coefs(:), rcen(:,:), xxg(:), zn(:)
    real(KREAL), allocatable :: ang(:), bin(:), dl(:,:,:), rl(:,:,:), rlt(:)
    real(KREAL), allocatable :: vaux(:), weights(:), weights01(:)
    real(KREAL) :: rstar, thresmult, thresoverlap
    real(KREAL) :: pi, raizpi, pimed, pimedsqr   ! Sqrt[pi/2]
    real(KREAL) :: roblk(-mxl:mxl)
    real(KREAL) :: fact(0:mxfact), facti(0:mxfact), re(-mxreal:mxreal), ri(-mxreal:mxreal)
    real(KREAL) :: root(0:mxroot), rooti(mxroot)
    real(KREAL) :: facts(-1:mxfact), factsi(-1:mxfact)
!       Coefficients and m values for the decomposition of the product of two Phi_m(phi) functions (of the real spherical harmonics)
    integer(KINT), allocatable :: msv(:,:), mdv(:,:), indk12(:,:)
    real(KREAL), allocatable :: app(:,:), bpp(:,:), ssv(:,:), sdv(:,:)
!       Polynomials P_k^(L,M;L',M') of the shift-operators technique 
    real(KREAL), allocatable :: polP(:)
!       Pointers to the elements P_0^(L,-L; L',-L')
    integer(KINT), allocatable  :: ipntpolP(:,:)
!
    logical :: langstrom, lgradient, lgrid, lgrid2d, lindividk, lindividlk, lindividlkm, lindividl, lpoints
    integer(KINT) :: idimzlm, kmaxrep, lmaxrep, lminrep, nindices
    integer(KINT), parameter ::  mxrtab = 200
    real(KREAL) :: dltu, dltv, dltx, dlty, dltz, uinf, usup, vinf, vsup, xinf, xsup, yinf, ysup, zinf, zsup
    real(KREAL) :: rtab(3,mxrtab)
    real(KREAL), allocatable :: cnlm(:,:), dgkl(:), radfunction(:), radderiv(:), zlm(:), zlmdx(:), zlmdy(:), zlmdz(:)
    character(256) :: filename, fileZJname	! Names for .plt and .zernike or .jacobi files
    integer(KINT), allocatable :: indicesv(:)
    character(256) :: x_func_uv, y_func_uv, z_func_uv  ! Expressions of (x,y,z) in terms of (u,v) for 2D grids
END MODULE
!...............................................................................................
!===============================================================================================
!                 MODULE PRECISION  (due to Dr. George Benthien:   http://www.gbenthien.net/strings/index.html)
!===============================================================================================
module precision

! Real kinds

integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

! Integer kinds

integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

!Complex kinds

integer, parameter :: kc4 = kr4                            ! single precision complex
integer, parameter :: kc8 = kr8                            ! double precision complex

end module precision
!
!                 END OF MODULE PRECISION
!...............................................................................................
!===============================================================================================
!                 MODULE STRINGS  (due to Dr. George Benthien:   http://www.gbenthien.net/strings/index.html)
!===============================================================================================
module strings

use precision

private :: value_dr,value_sr,value_di,value_si
private :: write_dr,write_sr,write_di,write_si
private :: writeq_dr,writeq_sr,writeq_di,writeq_si

interface value  ! Generic operator for converting a number string to a 
                 ! number. Calling syntax is 'call value(numstring,number,ios)' 
                 ! where 'numstring' is a number string and 'number' is a 
                 ! real number or an integer (single or double precision).         
   module procedure value_dr
   module procedure value_sr
   module procedure value_di
   module procedure value_si
end interface

interface writenum  ! Generic  interface for writing a number to a string. The 
                    ! number is left justified in the string. The calling syntax
                    ! is 'call writenum(number,string,format)' where 'number' is
                    ! a real number or an integer, 'string' is a character string
                    ! containing the result, and 'format' is the format desired, 
                    ! e.g., 'e15.6' or 'i5'.
   module procedure write_dr
   module procedure write_sr
   module procedure write_di
   module procedure write_si
end interface

interface writeq  ! Generic interface equating a name to a numerical value. The
                  ! calling syntax is 'call writeq(unit,name,value,format)' where
                  ! unit is the integer output unit number, 'name' is the variable
                  ! name, 'value' is the real or integer value of the variable, 
                  ! and 'format' is the format of the value. The result written to
                  ! the output unit has the form <name> = <value>.
   module procedure writeq_dr
   module procedure writeq_sr
   module procedure writeq_di
   module procedure writeq_si
end interface


!**********************************************************************

contains

!**********************************************************************

subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

character(len=*) :: str,delims
character(len=len_trim(str)) :: strsav
character(len=*),dimension(:) :: args

strsav=str
call compact(str)
na=size(args)
do i=1,na
  args(i)=' '
end do  
nargs=0
lenstr=len_trim(str)
if(lenstr==0) return
k=0

do
   if(len_trim(str) == 0) exit
   nargs=nargs+1
   call split(str,delims,args(nargs))
   call removebksl(args(nargs))
end do   
str=strsav

end subroutine parse

!**********************************************************************

subroutine compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str)):: outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  
  select case(ich)
  
    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1
      
    case(33:)      ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0
      
  end select
  
end do

str=adjustl(outstr)

end subroutine compact

!**********************************************************************

subroutine removesp(str)

! Removes spaces, tabs, and control characters in string str

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)    
    case(0:32)  ! space, tab, or control character
         cycle       
    case(33:)  
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!**********************************************************************

subroutine value_dr(str,rnum,ios)

! Converts number string to a double precision real number

character(len=*)::str
real(kr8)::rnum
integer :: ios

ilen=len_trim(str)
ipos=scan(str,'Ee')
if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
   ios=3
   return
end if
read(str,*,iostat=ios) rnum

end subroutine value_dr

!**********************************************************************

subroutine value_sr(str,rnum,ios)

! Converts number string to a single precision real number

character(len=*)::str
real(kr4) :: rnum
real(kr8) :: rnumd 

call value_dr(str,rnumd,ios)
if( abs(rnumd) > huge(rnum) ) then
  ios=15
  return
end if
if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
rnum=rnumd

end subroutine value_sr

!**********************************************************************

subroutine value_di(str,inum,ios)

! Converts number string to a double precision integer value

character(len=*)::str
integer(ki8) :: inum
real(kr8) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki8)

end subroutine value_di

!**********************************************************************

subroutine value_si(str,inum,ios)

! Converts number string to a single precision integer value

character(len=*)::str
integer(ki4) :: inum
real(kr8) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki4)

end subroutine value_si

!**********************************************************************

subroutine shiftstr(str,n)
 
! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift 
! are replaced by spaces.

character(len=*):: str

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
return

end subroutine shiftstr

!**********************************************************************

subroutine insertstr(str,strins,loc)

! Inserts the string 'strins' into the string 'str' at position 'loc'. 
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are 
! removed prior to insertion

character(len=*):: str,strins
character(len=len(str))::tempstr

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!**********************************************************************

subroutine delsubstr(str,substr)

! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.

character(len=*):: str,substr

lensubstr=len_trim(substr)
ipos=index(str,substr)
if(ipos==0) return
if(ipos == 1) then
   str=str(lensubstr+1:)
else
   str=str(:ipos-1)//str(ipos+lensubstr:)
end if   
return

end subroutine delsubstr

!**********************************************************************

subroutine delall(str,substr)

! Deletes all occurrences of substring 'substr' from string 'str' and
! shifts characters left to fill holes.

character(len=*):: str,substr

lensubstr=len_trim(substr)
do
   ipos=index(str,substr)
   if(ipos == 0) exit
   if(ipos == 1) then
      str=str(lensubstr+1:)
   else
      str=str(:ipos-1)//str(ipos+lensubstr:)
   end if
end do   
return

end subroutine delall

!**********************************************************************

function uppercase(str) result(ucstr)

! convert string to upper case

character (len=*):: str
character (len=len_trim(str)):: ucstr

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')     
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do
return

end function uppercase

!**********************************************************************

function lowercase(str) result(lcstr)

! convert string to lower case

character (len=*):: str
character (len=len_trim(str)):: lcstr

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do
return

end function lowercase

!**********************************************************************

subroutine readline(nunitr,line,ios)

! Reads line from unit=nunitr, ignoring blank lines
! and deleting comments beginning with an exclamation point(!)

character (len=*):: line

do  
  read(nunitr,'(a)', iostat=ios) line      ! read input line
  if(ios /= 0) return
  line=adjustl(line)
  ipos=index(line,'!')
  if(ipos == 1) cycle
  if(ipos /= 0) line=line(:ipos-1)
  if(len_trim(line) /= 0) exit
end do
return

end subroutine readline

!**********************************************************************

subroutine match(str,ipos,imatch)

! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.

character(len=*) :: str
character :: delim1,delim2,ch

lenstr=len_trim(str)
delim1=str(ipos:ipos)
select case(delim1)
   case('(')
      idelim2=iachar(delim1)+1
      istart=ipos+1
      iend=lenstr
      inc=1
   case(')')
      idelim2=iachar(delim1)-1
      istart=ipos-1
      iend=1
      inc=-1
   case('[','{','<')
      idelim2=iachar(delim1)+2
      istart=ipos+1
      iend=lenstr
      inc=1
   case(']','}','>')
      idelim2=iachar(delim1)-2
      istart=ipos-1
      iend=1
      inc=-1
   case default
      write(*,*) delim1,' is not a valid delimiter'
      return
end select
if(istart < 1 .or. istart > lenstr) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if
delim2=achar(idelim2) ! matching delimiter

isum=1
do i=istart,iend,inc
   ch=str(i:i)
   if(ch /= delim1 .and. ch /= delim2) cycle
   if(ch == delim1) isum=isum+1
   if(ch == delim2) isum=isum-1
   if(isum == 0) exit
end do
if(isum /= 0) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if   
imatch=i

return

end subroutine match

!**********************************************************************

subroutine write_dr(rnum,str,fmt)

! Writes double precision real number rnum to string str using format fmt

real(kr8) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_dr

!***********************************************************************

subroutine write_sr(rnum,str,fmt)

! Writes single precision real number rnum to string str using format fmt

real(kr4) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_sr

!***********************************************************************

subroutine write_di(inum,str,fmt)

! Writes double precision integer inum to string str using format fmt

integer(ki8) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_di

!***********************************************************************

subroutine write_si(inum,str,fmt)

! Writes single precision integer inum to string str using format fmt

integer(ki4) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_si

!***********************************************************************

subroutine trimzero(str)

! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.

character(len=*) :: str
character :: ch
character(len=10) :: exp

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle          
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!**********************************************************************

subroutine writeq_dr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr8) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_dr

!**********************************************************************

subroutine writeq_sr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr4) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_sr

!**********************************************************************

subroutine writeq_di(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki8) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_di

!**********************************************************************

subroutine writeq_si(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki4) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_si

!**********************************************************************

function is_letter(ch) result(res)

! Returns .true. if ch is a letter and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('A':'Z','a':'z')
  res=.true.
case default
  res=.false.
end select
return

end function is_letter

!**********************************************************************

function is_digit(ch) result(res)

! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

character :: ch
logical :: res

select case(ch)
case('0':'9')
  res=.true.
case default
  res=.false.
end select
return

end function is_digit

!**********************************************************************

subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

character(len=*) :: str,delims,before
character,optional :: sep
logical :: pres
character :: ch,cha

pres=present(sep)
str=adjustl(str)
call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\') then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if
   ipos=index(delims,ch)         
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)
   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do
if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return

end subroutine split

!**********************************************************************

subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

character(len=*):: str
character(len=1):: ch
character(len=len_trim(str))::outstr

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive

do i=1,lenstr
  ch=str(i:i)
  if(ibsl == 1) then          ! backslash active
   k=k+1
   outstr(k:k)=ch
   ibsl=0
   cycle
  end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!**********************************************************************

end module strings  
!
!                 END OF MODULE STRINGS
!...............................................................................................
!===============================================================================================
!                 MODULE EVALUATE  (due to Dr. George Benthien:   http://www.gbenthien.net/strings/index.html)
!===============================================================================================
module evaluate

! The user can assign values to parameters that can be used in expressions with 
! the subroutine defparam. The calling syntax is
!
!     call defparam(symbol,value) or call defparam(symbol,expr)
!
! where symbol is the desired parameter name; value is a real, integer, or 
! complex variable (single or double precision); and expr is a string
! containing an expression to be evaluated. The value obtained by evaluating the
! expression expr is associated with the parameter symbol. Parameter names must
! begin with a letter (a-z, A-Z) and must not be longer than 24 characters.
! Parameter names are not case dependent.
!
! An expression can be evaluated with the subroutine evalexpr. The calling
! syntax is 
!
!          call evalexpr(expr,value)
!
! where expr is a string containing the expression to be evaluated; value is the
! result (single or double precision real, complex or integer). The 
! expression can contain the arithmetic operations +, -, *, /, or ^ as well as
! the functions sin, cos, tan, log, ln, abs, exp, sqrt, real, imag, conjg, and 
! ang (the function ang calculates the phase angle of its complex argument). The
! expression can also contain numerical values and previously defined parameters
! Grouping by nested levels of parentheses is also allowed. The parameters pi
! and i (imaginary unit) are predefined. Complex numbers can be entered as a+i*b 
! if the parameter i has not been redefined by the user. Complex numbers can also
! be entered using complex(a,b).
! Example expression: 
!
!          conjg(((cos(x) + sqrt(a+i*b))^2+complex(ln(1.6e-4),20))/2)
!
! An equation of the form <symbol> = <expression> can be evaluated using the 
! subroutine evaleqn. The calling syntax is
!
!          call evaleqn(eqn)       
!
! where eqn is a string containing the equation. The right-hand-side of the
! equation is evaluated and assigned to the symbol given by the left-hand-side.
!
! The value assigned to a symbol can be retrieved using the subroutine getparam.
! The calling syntax is
!
!          call getparam(sym,value)
!
! where sym is a symbol string; value is a numeric variable (any of the six 
! standard types).
!
! The symbols and their values in the symbol table can be listed using the
! subroutine listvar. The variable ierr is always available following a call
! to any of the above subroutines and is zero if there were no errors. The
! possible nonzero values for ierr are
!
! 1       Expression empty
! 2       Parentheses don't match
! 3       Number string does not correspond to a valid number
! 4       Undefined symbol
! 5       Less than two operands for binary operation
! 6       No operand for unary plus or minus operators
! 7       No argument(s) for function
! 8       Zero or negative real argument for logarithm
! 9       Negative real argument for square root
! 10      Division by zero
! 11      Improper symbol format
! 12      Missing operator
! 13      Undefined function
! 14      Argument of tangent function a multiple of pi/2
!
use precision
use strings

save
private
public :: valuep,evalexpr,defparam,evaleqn,getparam,listvar,ierr

type item                         
  character(len=24):: char
  character :: type
end type item

type param                        
  character (len=24):: symbol
  complex(kc8):: value
end type param

interface defparam                
  module procedure strdef       ! value given by expression      
  module procedure valdef_dc    ! Double precision complex value
  module procedure valdef_sc    ! Single precision complex value
  module procedure valdef_dr    ! Double precision real value
  module procedure valdef_sr    ! Single precision real value
  module procedure valdef_di    ! Double precision integer value
  module procedure valdef_si    ! Single precision integer value
end interface

interface evalexpr
  module procedure evalexpr_dc  ! Double precision complex result
  module procedure evalexpr_sc  ! Single precision complex result
  module procedure evalexpr_dr  ! Double precision real result
  module procedure evalexpr_sr  ! Single precision real result
  module procedure evalexpr_di  ! Double precision integer result
  module procedure evalexpr_si  ! Single precision integer result
end interface

interface getparam
  module procedure getparam_dc  ! Double precision complex result
  module procedure getparam_sc  ! Single precision complex result
  module procedure getparam_dr  ! Double precision real result
  module procedure getparam_sr  ! Single precision real result
  module procedure getparam_di  ! Double precision integer result
  module procedure getparam_si  ! Single precision integer result
end interface

integer,parameter :: numtok=100  ! Maximum number of tokens 
type(param) :: params(100)       ! Symbol table
integer :: nparams=0,itop,ibin
complex(kc8) :: valstack(numtok) ! Stack used in evaluation of expression
type(item):: opstack(numtok)     ! Operator stack used in conversion to postfix
integer :: ierr                  ! Error flag

!**********************************************************************

contains

!**********************************************************************


SUBROUTINE EVALEXPR_DC(expr,val)    ! Evaluate expression expr for
                                    ! val double precision complex

character (len=*),intent(in) :: expr
complex(kc8) :: val
character (len=len(expr)+1) :: tempstr
character :: cop
integer :: isp(numtok)          ! On stack priority of operators in opstack
integer :: lstr
complex(kc8) :: cval,oper1,oper2
real(kr8) :: valr,vali
type(item):: token(numtok)      ! List of tokens ( a token is an operator or 
                                ! operand) in postfix order
type(item) :: x,junk,tok

ierr=0
token(1:)%char=' '

if(nparams == 0) then                  ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if

if(len_trim(expr) == 0) then           ! Expression empty
  ierr=1
  write(*,*) 'Error: expression being evaluated is empty'
  return
end if

tempstr=adjustl(expr)
call removesp(tempstr)   ! Removes spaces, tabs, and control characters

! ****************************************************************************
! STEP 1:  Convert string to token array. Each token is either an operator or
!          an operand. Token array will be in postfix (reverse Polish) order.
!*****************************************************************************

ntok=0
ibin=0
itop=0
do
  lstr=len_trim(tempstr)
  call get_next_token(tempstr(1:lstr),tok,icp,insp)
  select case(tok%type)
  case('S')
    ntok=ntok+1
    token(ntok)=tok
  case('E')
    do
      if(itop < 1)exit
      call popop(x)        ! Output remaining operators on stack
      ntok=ntok+1
      token(ntok)=x
    end do
    ntok=ntok+1
    token(ntok)=tok
    exit
  case('R')  ! Token is right parenthesis
    do
      if(itop .le. 0 .or. opstack(itop)%type == 'L') exit  ! Output operators on stack down
      call popop(x)                       ! to left parenthesis
      ntok=ntok+1
      token(ntok)=x
    end do
    call popop(junk)                      ! Remove left parenthesis from stack
    if(itop .gt. 0 .and. opstack(itop)%type == 'F') then    ! Output function name if present
      call popop(x)
      ntok=ntok+1
      token(ntok)=x
    end if
  case('D')  ! Token is comma
    do
      if(itop .le. 0 .or. opstack(itop)%type == 'L') exit  ! Output operators on stack down
      call popop(x)                       ! to left parenthesis
      ntok=ntok+1
      token(ntok)=x
    end do
  case('U','B','L','F') ! Token is operator, left parenthesis or function name
    do
      if(itop .le. 0 .or. isp(itop) < icp) exit            ! Output operators on stack having
      call popop(x)                       ! an instack priority that is
      ntok=ntok+1                         ! greater than or equal to the
      token(ntok)=x                       ! priority of the incoming operator
    end do
    call pushop(tok)     ! Put incoming operator on stack
    isp(itop)=insp
  end select
end do

!write(6,"('token = ', 80(a, ' '))") token(1:ntok)%type
!write(6,"('token = ', 80(a, ' '))") token(1:ntok)%char(1:2)

isum=0                                 ! Error check for matching parentheses
do i=1,ntok
  if(token(i)%type == 'L' ) isum=isum+1
  if(token(i)%type == 'R' ) isum=isum-1
end do
if(isum /= 0) then
  ierr=2
  write(*,*) 'Error in the evaluation of the expression ',trim(expr)
  write(*,*) "Parentheses don't match"
  write(*,*)
  return
end if


!*****************************************************************************
! STEP 2: Evaluate token string in postfix order
!*****************************************************************************

itop=0
do i=1,ntok
  x=token(i)
  select case(x%type)
  case('E')  ! Token is end token
    if(itop>1) then                
      ierr=12
      write(*,*) 'Error: missing operator in expression ',trim(expr)
      write(*,*)
      return
    end if
    call popval(val)               ! Final result left on stack of values
    exit
  case('S')  ! Token is operand
    call valuep(x%char,cval)       ! Evaluate operand
    if(ierr/=0) return
    call pushval(cval)             ! Put value of operand on stack
  case('B')  ! Token is a binary operator
    if(itop < 2) then
      ierr=5
      write(*,*) 'Error in evaluation of expression ',trim(expr)
      write(*,*) 'Less than two operands for binary operator  '&
                 ,trim(x%char)
      write(*,*)
      return
    end if                         
    call popval(oper1)             ! Pull off top two values from stack
    call popval(oper2)
    select case(trim(x%char))      ! Perform operation on values
    case('^')
      cval=oper2**oper1
    case('*')
      cval=oper2*oper1
    case('/')
      if(oper1 == (0._kr8,0._kr8)) then
        ierr=10
        write(*,*) 'Error in expression ',trim(expr)
        write(*,*) 'Division by zero'
        write(*,*)
        return
      end if
      cval=oper2/oper1
    case('+')
      cval=oper2+oper1
    case('-')
      cval=oper2-oper1
    end select
    call pushval(cval)             ! Put result back on stack
  case('U')  ! Token is unary operator
    if(itop == 0) then
      ierr=6
      write(*,*) 'Error in expression ',trim(expr)
      write(*,*) 'No operand for unary operator ',trim(x%char)
      write(*,*)
      return
    else
      call popval(oper1)           ! Pull top value off stack
    end if
    select case(trim(x%char))      ! Operate on value
    case('+')
      cval=oper1
    case('-')
      cval=-oper1
    end select
    call pushval(cval)             ! Put result back on stack
  case('F')  ! Token is a function name
    if(itop == 0) then
      ierr=7
      write(*,*) 'Error in expression ',trim(expr)
      write(*,*) 'Missing argument(s) for function ',trim(x%char)
      write(*,*)
      return
    else  
      call popval(oper1)           ! Pull top value off stack
    end if 
    tempstr=uppercase(x%char)
    select case(trim(tempstr))      ! Evaluate function
    case('SIN')
      cval=sin(oper1)
    case('COS')
      cval=cos(oper1)
    case('TAN')
      oper2=cos(oper1)
      if(abs(oper2) == 0.0_kr8) then
        ierr=14
        write(*,*) 'Error: argument of tan function a multiple',&
        ' of pi/2 in expression ',trim(expr)
        write(*,*)
        return
      else 
        cval=sin(oper1)/oper2
      endif
    case('SQRT')
      if(real(oper1,kr8) < 0. .and. aimag(oper1)==0.) then
        ierr=9
        write(*,*) 'Warning: square root of negative real number',&
                   ' in expression ',trim(expr)
        write(*,*)
      end if
      cval=sqrt(oper1)
    case('ABS')
      cval=abs(oper1)
    case('LN')
      if(real(oper1,kr8) <= 0. .and. aimag(oper1)==0.) then
        ierr=8
        write(*,*) 'Error: negative real or zero argument for',&
                   ' natural logarithm in expression ',trim(expr)
        write(*,*)
        return
      end if
      cval=log(oper1)
    case('LOG')
      if(real(oper1,kr8) <= 0. .and. aimag(oper1)==0.) then
        ierr=8
        write(*,*) 'Error: negative real or zero argument for base',&
                   '10 logarithm in expression ',trim(expr)
        write(*,*)
        return
      end if
      cval=log(oper1)/2.30258509299405_kr8
    case('EXP')
      cval=exp(oper1)
    case('COMPLEX')
      if(itop == 0) then
        ierr=7
        write(*,*) 'Error in expression ',trim(expr)
        write(*,*) 'Missing argument(s) for function ',trim(x%char)
        write(*,*)
        return
      else  
        call popval(oper2)  ! Pull second argument off stack
      end if 
      valr=real(oper2,kr8)
      vali=real(oper1,kr8)
      cval=cmplx(valr,vali,kc8)
    case('CONJG')
      cval=conjg(oper1)
    case('ANG')
      cval=atan2(aimag(oper1),real(oper1,kr8))
    case('REAL')
      cval=real(oper1,kr8)
    case('IMAG')
      cval=aimag(oper1)
    case default ! Undefined function
      ierr=13
      write(*,*) 'Error: the function ',trim(x%char), ' is undefined',&
                 ' in the expression ',trim(expr)
      write(*,*)
      return
    end select
    call pushval(cval)    ! Put result back on stack
  end select
end do

end subroutine evalexpr_dc

!**********************************************************************

SUBROUTINE GET_NEXT_TOKEN(str,tok,icp,isp)

character(len=*) :: str
character :: cop,chtemp
type(item) :: tok
integer :: icp

lstr=len_trim(str)
if(lstr == 0) then
  tok%char='#'             ! Output end token
  tok%type='E'
  return
end if
ipos=scan(str,'+-*/^(),')  ! Look for an arithmetic operator 
                           ! + - * / ^ ( ) or ,
cop=str(ipos:ipos)                 
select case (ipos)              
case(0)    ! Operators not present
  ntok=ntok+1
  tok%char=str
  tok%type='S'
  str=''
  icp=0
  isp=0
case(1) 
  tok%char=cop
  select case(cop)
  case('+','-')
    if(ibin==0) then
      tok%type='U'
      icp=4
      isp=3
    else
      tok%type='B'
      icp=1
      isp=1
    end if
    ibin=0
  case('*','/')
    tok%type='B'
    icp=2
    isp=2
    ibin=0
  case('^')
    tok%type='B'
    icp=4
    isp=3
    ibin=0
  case('(')
    tok%type='L'
    icp=4
    isp=0
    ibin=0
  case(')')
    tok%type='R'
    icp=0
    isp=0
    ibin=1
  case(',')
    tok%type='D'
    icp=0
    isp=0
    ibin=0
  end select
  str=str(2:)
case(2:)
  select case(cop)
  case('(')
    tok%char=str(1:ipos-1)
    tok%type='F'
    icp=4
    isp=0
    ibin=0
    str=str(ipos:)
  case('+','-')
    chtemp=uppercase(str(ipos-1:ipos-1))
    if(is_letter(str(1:1)) .or. chtemp/='E') then
      tok%char=str(1:ipos-1)
      tok%type='S'
      icp=0
      isp=0
      ibin=1
      str=str(ipos:)
    else
      inext=scan(str(ipos+1:),'+-*/^(),')
      if(inext==0) then
        tok%char=str
        tok%type='S'
        icp=0
        isp=0
        ibin=0
        str=''
      else
        tok%char=str(1:ipos+inext-1)
        tok%type='S'
        icp=0
        isp=0
        ibin=1
        str=str(ipos+inext:)
      end if
    end if
  case default
    tok%char=str(1:ipos-1)
    tok%type='S'
    icp=0
    isp=0
    ibin=1
    str=str(ipos:)
  end select
end select

end subroutine get_next_token


!**********************************************************************

SUBROUTINE EVALEXPR_SC(expr,val)    ! Evaluate expression expr for
                                    ! val single precision complex
character(len=*) :: expr
complex(kc4) :: val
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
val=vald

end subroutine evalexpr_sc

!**********************************************************************

SUBROUTINE EVALEXPR_SR(expr,val)    ! Evaluate expression expr for
                                    ! val single precision real
character(len=*) :: expr
real(kr4) :: val
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
val=real(vald)

end subroutine evalexpr_sr

!**********************************************************************

SUBROUTINE EVALEXPR_DR(expr,val)    ! Evaluate expression expr for 
                                    ! val double precision real
character(len=*) :: expr
real(kr8) :: val
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
val=real(vald,kr8)

end subroutine evalexpr_dr

!**********************************************************************

SUBROUTINE EVALEXPR_SI(expr,ival)   ! Evaluate expression expr for 
                                    ! ival single precision integer
character(len=*) :: expr
integer(ki4) :: ival
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
ival=nint(real(vald,kr8),ki4)

end subroutine evalexpr_si

!**********************************************************************

SUBROUTINE EVALEXPR_DI(expr,ival)   ! Evaluate expression expr for 
                                    ! ival double precision integer
character(len=*) :: expr
integer(ki8) :: ival
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
ival=nint(real(vald,kr8),ki8)

end subroutine evalexpr_di


!**********************************************************************
SUBROUTINE VALDEF_DC(sym,val)    ! Associates sym with val in symbol table,
                                 ! val double precision complex
character(len=*) :: sym
character(len=len_trim(sym)) :: usym
complex(kc8) :: val

ierr=0
if(nparams == 0) then               ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if

! Assign val to sym if sym is already in symbol table
usym=uppercase(sym)
if(.not. is_letter(sym(1:1)) .or. len_trim(sym)>24) then
  ierr=11
  write(*,*) 'Error: symbol ',trim(sym),' has improper format'
  write(*,*)
  return
end if
do i=1,nparams
  if(trim(usym)==trim(params(i)%symbol)) then
    params(i)%value=val
    return
  end if
end do

nparams=nparams+1    ! Otherwise assign val to new symbol sym
params(nparams)%symbol=usym
params(nparams)%value=val

end subroutine valdef_dc


!**********************************************************************

SUBROUTINE VALDEF_SC(sym,val)     ! Associates sym with val in symbol table,
                                  ! val single precision complex
character(len=*) :: sym
complex(kc4) :: val
complex(kc8) :: vald

vald=val
call valdef_dc(sym,vald)

end subroutine valdef_sc


!**********************************************************************

SUBROUTINE VALDEF_DR(sym,val)    ! Associates sym with val in symbol table,
                                 ! val double precision real
character(len=*) :: sym
real(kr8) :: val
complex(kc8) :: vald

vald=cmplx(val,0.0_kr8,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_dr


!**********************************************************************

SUBROUTINE VALDEF_SR(sym,val)    ! Associates sym with val in symbol table,
                                 ! val single precision real
character(len=*) :: sym
real(kr4) :: val
complex(kc8) :: vald

vald=cmplx(val,0.0,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_sr


!**********************************************************************

SUBROUTINE VALDEF_DI(sym,ival)   ! Associates sym with ival in symbol table,
                                 ! ival double precision integer 
character(len=*) :: sym
integer(ki8) :: ival
complex(kc8) :: vald

vald=cmplx(real(ival,kr8),0.0_kr8,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_di


!**********************************************************************

SUBROUTINE VALDEF_SI(sym,ival)   ! Associates sym with ival in symbol table,
                                 ! ival single precision integer
character(len=*) :: sym
integer(ki4) :: ival
complex(kc8) :: vald

vald=cmplx(real(ival,kr8),0.0,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_si


!**********************************************************************

SUBROUTINE STRDEF(sym,expr)      ! Associates sym with the value of the
                                 ! expression expr

character(len=*) :: sym,expr
complex(kc8) :: val

if(nparams == 0) then            ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if

call evalexpr_dc(expr,val)       ! val is value of expression expr
if(ierr==0 .or. ierr==9) then
  call valdef_dc(sym,val)          ! Assign val to symbol sym
end if

end subroutine strdef


!**********************************************************************

SUBROUTINE VALUEP(xinchar,cval)  ! Finds double precision complex value 
                                 ! corresponding to number string xinchar 
                                 ! or value in symbol table corresponding 
                                 ! to symbol name xinchar.

character (len=*):: xinchar
complex(kc8) :: cval
real(kr8) :: rval

ierr=0

if(is_letter(xinchar(1:1))) then   ! xinchar is a symbol
  call getparam(xinchar,cval)
else                               ! xinchar is a number string
  call value(xinchar,rval,ios)     ! rval is the value of xinchar
  if(ios > 0) then
    ierr=3
    write(*,*) 'Error: number string ',trim(xinchar),' does not correspond to a valid number' 
    write(*,*)
  end if
  cval=cmplx(rval,0.0_kr8,kc8)
  return
end if

end subroutine valuep


!**********************************************************************


SUBROUTINE PUSHOP(op)  ! Puts an operator on operator stack

type(item):: op

itop=itop+1
if(itop > numtok) then
  write(*,*) 'Error: operator stack overflow in evaluation of expression'
  write(*,*)
  return
end if
opstack(itop)=op

end subroutine pushop

SUBROUTINE POPOP(op) ! Takes top operator of operator stack and assigns it to op

type(item):: op

op=opstack(itop)
itop=itop-1

end subroutine popop

SUBROUTINE PUSHVAL(val) ! Puts value on value stack

complex(kc8) :: val

itop=itop+1
if(itop > numtok) then
  write(*,*) 'Error: value stack overflow in evaluation of expression'
  write(*,*)
  return
end if
valstack(itop)=val

end subroutine pushval

SUBROUTINE POPVAL(val) ! Takes top value off value stack and assigns it to val

complex(kc8) :: val

val=valstack(itop)
itop=itop-1

end subroutine popval

!**********************************************************************

SUBROUTINE GETPARAM_DC(sym,var)  ! Find double precision complex value var
                                 ! corresponding to symbol sym

character(len=*) :: sym
character(len=len_trim(sym)) :: usym
complex(kc8) :: var

ierr=0
sym=adjustl(sym)
if(.not.is_letter(sym(1:1)) .or. len_trim(sym)>24) then
  ierr=11
  write(*,*) 'Error: symbol ',trim(sym),' has incorrect format'
  write(*,*)
  return
end if
ifind=0
usym=uppercase(sym)
do j=1,nparams
  if(trim(usym) == trim(params(j)%symbol)) then
    var=params(j)%value
    ifind=j
    exit
  end if
end do
if(ifind == 0) then          
  ierr=4
  write(*,*) 'Error: symbol ',trim(sym), ' not in symbol table'
  write(*,*) 
  return
end if

end subroutine getparam_dc

!**********************************************************************

SUBROUTINE GETPARAM_SC(sym,var)  ! Find single precision complex value var
                                 ! corresponding to symbol sym


character(len=*) :: sym
complex(kc4) :: var
complex(kc8) :: vard

call getparam_dc(sym,vard)
var=vard

end subroutine getparam_sc

!**********************************************************************

SUBROUTINE GETPARAM_DR(sym,var)  ! Find double precision real value var
                                 ! corresponding to symbol sym


character(len=*) :: sym
real(kr8) :: var
complex(kc8) :: vard

call getparam_dc(sym,vard)
var=real(vard,kr8)

end subroutine getparam_dr

!**********************************************************************

SUBROUTINE GETPARAM_SR(sym,var)  ! Find single precision real value var
                                 ! corresponding to symbol sym


character(len=*) :: sym
real(kr4) :: var
complex(kc8) :: vard

call getparam_dc(sym,vard)
var=real(vard)

end subroutine getparam_sr

!**********************************************************************

SUBROUTINE GETPARAM_DI(sym,ivar)  ! Find double precision integer value ivar
                                  ! corresponding to symbol sym


character(len=*) :: sym
integer(ki8) :: ivar
complex(kc8) :: vard

call getparam_dc(sym,vard)
ivar=nint(real(vard,kr8),ki8)

end subroutine getparam_di

!**********************************************************************

SUBROUTINE GETPARAM_SI(sym,ivar)  ! Find single precision integer value ivar
                                  ! corresponding to symbol sym


character(len=*) :: sym
integer(ki4) :: ivar
complex(kc8) :: vard

call getparam_dc(sym,vard)
ivar=nint(real(vard,kr8),ki4)

end subroutine getparam_si

!**********************************************************************

SUBROUTINE EVALEQN(eqn)  ! Evaluate an equation

character(len=*) :: eqn
character(len=len(eqn)) :: args(2)
complex(kc8) :: val

call parse(eqn,'=',args,nargs)   ! Seperate right- and left-hand-sides
call defparam(adjustl(args(1)),args(2)) ! Evaluate right-hand-side and
                                        ! assign to symbol on the
                                        ! left-hand-side.
end subroutine evaleqn

!**********************************************************************

SUBROUTINE LISTVAR      ! List all variables and their values

write(*,'(/a)') ' VARIABLE LIST:'
if(nparams == 0) then            ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if
do i=1,nparams
  write(*,*) trim(params(i)%symbol),' = ',params(i)%value
end do

end subroutine listvar

!**********************************************************************

end module evaluate
!
!                 END OF MODULE EVALUATE
!...............................................................................................



  program DAMDENZJ_atdens_320
    USE DAMDENZJ_atdens
    USE strings, only: lowercase
    implicit none
    logical :: existe, lindivid
    integer(KINT) :: i, ierr, ipoint, k, kmax, knt, ktop, l, ltop, m, nu, numrtab
    integer(KINT), allocatable :: indices(:)
    real(KREAL) :: aux, bux, den, dxden, dyden, dzden, x, y, z
    real(KREAL4) :: tarray(2), tiempo, dtime
    namelist / options / dltu, dltv, dltx, dlty, dltz, filename, fileZJname, indices, iswindows, kmaxrep, lechelon, lexact, &
            lgradient, lgrid, lgrid2d, lindividk, lindividlk, lindividlkm, lindividl, ljacobi, lmaxrep, lminrep, lpoints, &
            numrtab, rtab, uinf, usup, vinf, vsup, x_func_uv , xinf, xsup, y_func_uv, yinf, ysup, z_func_uv, zinf, zsup
    external zernike3DR, jacobiP
    tiempo = dtime(tarray)
!    Defaults for the NAMELIST OPTIONS
    filename = ""            ! root file name for .plt and .pltd files
    fileZJname = ""        ! .zernike or .jacobi files
    iswindows = .false.        ! .true. if running on a MS-windows system
    xinf = 0.d0
    xsup = 0.d0
    dltx = 1.d0
    yinf = 0.d0
    ysup = 0.d0
    dlty = 1.d0
    zinf = 0.d0
    zsup = 0.d0
    dltz = 1.d0
    x_func_uv = 'u'        ! x = u for 2D grids:  default: plane XY
    y_func_uv = 'v'        ! y = v for 2D grids
    z_func_uv = '0'        ! z = 0 for 2D grids
    uinf = 0.d0
    usup = 1.d0
    dltu = 1.d0
    vinf = 0.d0
    vsup = 1.d0
    dltv = 1.d0
    langstrom = .true.      ! If false distances in bohr
    lechelon = .false.      ! if true number of functions per l equal to max(lexpansion+1,nexpansion)-l
    lexact = .false.        ! If true computes the density fro its atomic expansion
    lindividk = .false.     ! If true projection functions with given values of k index and all l and m compatible indices are used
    lindividl = .false.     ! If true projection functions with given values of l index and all k and m compatible indices are used
    lindividlk = .false.    ! If true projection functions with given values of k, l indices and all m compatible indices are used
    lindividlkm = .false.   ! If true only projection functions with given values of k, l, m indices are used
    ljacobi = .false.       ! If .true. expansion in Jacobi functions
    lgrid = .true.          ! If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    lgradient = .false.     ! If .true. gradient components of the density computed and, if lgrid = .true., tabulated in files
                                            !    projectname-d-dx.pltd, projectname-d-dy.pltd, projectname-d-dz.pltd
    lgrid2d = .false.       ! If .true. computes a 2D grid. (x,y,z) are given in terms of (u,v)
    lpoints = .false.       ! If .true. computes in selected points and prints out the results in the standard output.
                                            ! Points must be given in cartesian coordinates.
                                            ! If lgrid .eq. .true., these coordinates must be placed after the grid data
    kmaxrep = 10            ! Highest k in expansion
    lminrep = 0             ! Lowest l in expansion
    lmaxrep = 10            ! Highest l in expansion
    numrtab = 0             ! Number of tabulation points supplied in namelist
!    End of Defaults for the NAMELIST OPTIONS

    allocate(indices(3000), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating indices. Stop')
    indices = -1000

    read(5,OPTIONS)    !    Reads the namelist OPTIONS
    read(5,*) projectname
    write(6,"(1x,'project name : ',a,/,1x,'==============')") projectname

    if (iswindows) then
        dirsep = "\\"
    else
        dirsep = "/"
    endif
    if (len_trim(filename) .eq.0) then
        filename = projectname
    else
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        if (i .eq. 0) i = index(projectname,"/",.true.)    ! In case of windows with MinGW directory separator is "/"
        filename = projectname(1:i)//trim(filename)
    endif
        
!     Tabulates the density from its atomic expansion
    
    if (lexact) then
        call density_exact
        if (allocated(zlm)) deallocate(zlm)
        if (allocated(zlmdx)) deallocate(zlmdx, zlmdy, zlmdz)
        write(6,*) '******* END OF DAMDENZJ_320 ******'
        tiempo = dtime(tarray)
        write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
                tarray(1), tarray(2), tarray(1)+tarray(2)
        stop
    endif
    
!     Tabulates the density from its Canterakis-Zernike or Jacobi expansion

    existe = .false.
    if (len_trim(fileZJname) .eq. 0) then
        inquire(file=trim(projectname)//".zernike", exist=existe)
        if (existe) then
            fileZJname = trim(projectname)//".zernike"
        else
            inquire(file=trim(projectname)//".jacobi", exist=existe)
            if (existe) then
                fileZJname = trim(projectname)//".jacobi"
            else
                call error(1,'Cannot find a file with one-center expansion named '//trim(projectname)//'+ .zernike or .jacobi')
            endif
        endif
    else
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        if (i .eq. 0) i = index(projectname,"/",.true.)    ! In case of windows with MinGW directory separator is "/"
        fileZJname = projectname(1:i)//trim(fileZJname)
        inquire(file=fileZJname, exist=existe)
        if (.not. existe) then
            call error(1, 'Cannot find a file with one-center expansion  named '//fileZJname)
        endif
    endif
    i = index(fileZJname,".",.true.)
    if (lowercase(fileZJname(i+1:len_trim(fileZJname))) .eq. 'zernike') then
        ljacobi = .false.
    else if (lowercase(fileZJname(i+1:len_trim(fileZJname))) .eq. 'jacobi') then
        ljacobi = .true.
    else
        call error(1, 'Extension of file with one-center expansion named '//fileZJname//' must be "zernike" or "jacobi"')
    endif

    lindivid = .false.
    if (lindividk .or. lindividl .or. lindividlk .or. lindividlkm) then
        nindices = 0
        do i = 1, 3000
            if (indices(i) .eq. -1000) cycle
            nindices = nindices+1
            indices(nindices) = indices(i)
        enddo
        if (nindices .gt. 0) then
            lindivid = .true.
            allocate(indicesv(nindices), stat = ierr)
            if (ierr .ne. 0) call error(ierr,'Memory error when allocating indicesv. Stop')
            do i = 1, nindices
                indicesv(i) = indices(i)
            enddo
        endif
    endif
    deallocate(indices)

    if (lindividlk .and. mod(nindices,2) .ne. 0) then
        call error(1,'Wrong number of indices for projection functions of selected l and k. Two indices are required for&
                & each function. Stop')
    endif

    if (lindividlkm .and. mod(nindices,3) .ne. 0) then
        call error(1,'Wrong number of indices for individual projection functions. Three indices are required for&
                & each function. Stop')
    endif

    write(6,"(/30x,'Using expansion in file: ', a)") trim(fileZJname)
    open(16,file=trim(fileZJname),form='formatted', iostat=ierr)
    if (ierr .ne. 0) call error(ierr,'Cannot open file '//trim(fileZJname)//'. Stop')
    read(16,*) rstar
    read(16,*) ltop, ktop
    if (kmaxrep .gt. ktop) then
        write(6,"(/'WARNING! highest value of k required: ', i5, ' greater than available.',&
                /' Sets it to the highest available = ', i5, / )") kmaxrep, ktop
        kmaxrep = ktop
    endif
    if (lmaxrep .gt. ltop) then
        write(6,"(/'WARNING! highest value of l required: ', i5, ' greater than available.',&
                /' Sets it to the highest available = ', i5, / )") lmaxrep, ltop
        lmaxrep = ltop
    endif
    write(6,"(/'****  rstar = ', e17.10,/)") rstar

    if (.not. lindivid) then
        write(6,"(/'lmaxrep = ', i3, ' kmaxrep = ', i3)") lmaxrep, kmaxrep
    elseif (lindividl) then
        write(6,"(/'Selected values of l: ',50(1x,i2))") indicesv(1:nindices)
        write(6,"('All values of k and m (compatible) taken'/)")
    elseif (lindividk) then
        write(6,"(/'Selected values of k: ',50(1x,i2))") indicesv(1:nindices)
        write(6,"('All values of l and m (compatible) taken'/)")
    elseif (lindividlk) then
        write(6,"(/'Selected values of (l,k): ',30(1x,'(',i2,',',i2,')'))") (indicesv(2*i-1),indicesv(2*i), i = 1, nindices/2)
        write(6,"('All values of m (compatible) taken'/)")
    elseif (lindividlkm) then
        write(6,"(/'Selected values of (l,k,m): ',20(1x,'(',i2,',',i2,',',i3,')'))") &
                (indicesv(3*i-2), indicesv(3*i-1),  indicesv(3*i), i = 1, nindices/3)
    endif

    allocate(omeganlm(0:ktop,(lmaxrep+1)*(lmaxrep+1)), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating omeganlm. Stop')
    if (lechelon) then
        write(6,"('Length of expansions taken in echelon form (kmax(l) = ', i3, '-l')") kmaxrep
    endif
    omeganlm = 0.d0
    do l = 0, lmaxrep
        do m = -l, l
            read(16,*) omeganlm(0:ktop,l*(l+1)+m+1)
        enddo
    enddo
    close(16)

    call consta
    allocate(gkl(-1:kmaxrep), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating gkl. Stop')
    allocate(radfunction(0:kmaxrep), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating radfunction. Stop')
    idimzlm = (max(4,lmaxrep)+2)**2
    allocate(zlm(idimzlm), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error in gridrep when allocating zlm. Stop')
    if (lgradient) then
        allocate(dgkl(-1:kmaxrep), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating dgkl. Stop')
        allocate(radderiv(0:kmaxrep), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating radderiv. Stop')
        allocate(zlmdx(idimzlm), zlmdy(idimzlm), zlmdz(idimzlm), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error in densrepr when allocating zlmdx, zlmdy, zlmdz. Stop')
    endif

#ifdef DBLPRCGRID
    write(6,"('Grid generated in double precision')")
    write(6,"('WARNING! This grid will not be compatible with gOpenMol')")
#endif
    if (lgrid) then
        if (xinf .eq. xsup) then
            xinf = -rstar
            xsup = rstar
        elseif (xinf .gt. xsup) then
            aux = xsup
            xsup = xinf
            xinf = aux
            dltx = abs(dltx)
        endif
        if (yinf .eq. ysup) then
            yinf = -rstar
            ysup = rstar
        elseif (yinf .gt. ysup) then
            aux = ysup
            ysup = yinf
            yinf = aux
            dlty = abs(dlty)
        endif
        if (zinf .eq. zsup) then
            zinf = -rstar
            zsup = rstar
        elseif (zinf .gt. zsup) then
            aux = zsup
            zsup = zinf
            zinf = aux
            dltz = abs(dltz)
        endif
!    Computes the density on the grid points
        if (.not. lindivid) then
            if (.not. lgrid2d) then
                call grid3D
            else
                call grid2D
            endif
        else
            if (.not. lgrid2d) then
                call grid3D_part
            else
                call grid2D_part
            endif
        endif
    endif
!    Tabulates specific points if required
    if (lpoints) then
        if (lgradient) then
            write(6,"(//12x,'Density expansion',/10x, 21('='), &
                    //9x,'X',17x,'Y',17x,'Z',20x,'den',t89,'der x', t112, 'der y', t135, 'der z')")
        else
            write(6,"(//12x,'Density expansion',/10x, 21('='), //9x,'X',17x,'Y',17x,'Z',t67,'den')")
        endif
        if (numrtab .gt. mxrtab) then
            write(6,"('Number of tabulation points numrtab = ', i4, ' exceeds the highest allowable = ', i4 &
                    / ' sets numrtab = ', i4)") numrtab, mxrtab, mxrtab
            numrtab = mxrtab
        endif
        ipoint = 1
        ierr = 0
        if (ipoint .le. numrtab) then
            x = rtab(1,ipoint)
            y = rtab(2,ipoint)
            z = rtab(3,ipoint)
            ipoint = ipoint + 1
        else
            read(5,*,iostat=ierr) x, y, z
        endif

        do while (ierr .eq. 0)
            den = 0.d0
            dxden = 0.d0
            dyden = 0.d0
            dzden = 0.d0
            if (.not. lindivid) then
                if (ljacobi) then
                    call compute_density(jacobiP, x, y, z, den, dxden, dyden, dzden)
                else
                    call compute_density(zernike3DR, x, y, z, den, dxden, dyden, dzden)
                endif
            else
                if (ljacobi) then
                    call compute_density_part(jacobiP, x, y, z, den, dxden, dyden, dzden)
                else
                    call compute_density_part(zernike3DR, x, y, z, den, dxden, dyden, dzden)
                endif
            endif

            if (lgradient) then
                            write(6,"(3(1x,e17.10), 3x, d22.15, 3(1x,e22.15))") x, y, z, den, dxden, dyden, dzden
            else
                    write(6,"(3(1x,e17.10), 3x, d22.15)") x, y, z, den
            endif

            if (ipoint .le. numrtab) then
                x = rtab(1,ipoint)
                y = rtab(2,ipoint)
                z = rtab(3,ipoint)
                ipoint = ipoint + 1
            else
                read(5,*,iostat=ierr) x, y, z
            endif
        enddo
    endif

    if (allocated(zlm)) deallocate(zlm)
    if (allocated(zlmdx)) deallocate(zlmdx, zlmdy, zlmdz)
    write(6,*) '******* END OF DAMDENZJ_320 ******'
    tiempo = dtime(tarray)
    write(6,"(1x,'Timing in seconds (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
            tarray(1), tarray(2), tarray(1)+tarray(2)
    stop
    end
!    
!   ***************************************************************
!
   subroutine compute_density(frad_fun, x, y, z, den, dxden, dyden, dzden)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: l, k, knt, m, nu
    real(KREAL) :: rvec(0:kmaxrep)
    real(KREAL) :: aux, bux, cux, den, dxden, dyden, dzden, dosl1, r, r2, rsinv, t, x, xa, y, ya, z, za
    interface
        subroutine frad_fun(lrad, trad)
            integer(4) :: lrad
            real(8) :: trad
        end subroutine frad_fun
    end interface
    den = 0.d0
    dxden = 0.d0
    dyden = 0.d0
    dzden = 0.d0
    if ((x*x + y*y + z*z) .gt. rstar*rstar) return
    rsinv = 1.d0 / rstar
    xa = x * rsinv
    ya = y * rsinv
    za = z * rsinv
    r2 = xa*xa + ya*ya + za*za
    r = sqrt(r2)
! write(6,"('r2 = ', e17.10, ' r = ', e17.10)") r2, r
    zlm(1) = 1.d0        ! Regular spherical harmonics of r/rstar = (xa,ya,za)
    zlm(2) = ya
    zlm(3) = za
    zlm(4) = xa
    do l = 1, lmaxrep
        dosl1 = dble(l+l+1)
        zlm((l+1)*(l+3)+1) = dosl1 * (xa * zlm(l*(l+2)+1) - ya * zlm(l*l+1))        ! zlm(l+1,l+1,ia)
        zlm((l+1)*(l+1)+1) = dosl1 * (ya * zlm(l*(l+2)+1) + xa * zlm(l*l+1))        ! zlm(l+1,-(l+1),ia)
        zlm((l+2)*(l+2)-1) = dosl1 * za * zlm(l*(l+2)+1)                ! zlm(l+1,l,ia)
        zlm(l*(l+2)+3) = dosl1 * za * zlm(l*l+1)                    ! zlm(l+1,-l,ia)
        do m = 0, l-1
            aux = 1.d0 / dble(l-m+1)
            bux = dble(l+m)
            zlm((l+1)*(l+2)+m+1) = aux * (dosl1*za*zlm(l*(l+1)+m+1) - bux*r2*zlm((l-1)*l+m+1))    ! zlm(l+1,m,ia)
            zlm((l+1)*(l+2)-m+1) = aux * (dosl1*za*zlm(l*(l+1)-m+1) - bux*r2*zlm((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
        enddo
    enddo
    if (lgradient) then
        call derivzlm(lmaxrep)
    endif
    if (ljacobi) then
        t = r
    else
        t = r2
    endif
    if (r .gt. 0.d0) then
        aux = 1.d0 / (r * rstar)
        bux = ya * aux
        cux = za * aux
        aux = xa * aux
    else
        aux = 1.d0
        bux = 1.d0
        cux = 1.d0
    endif
    den = 0.d0
    knt = lminrep*lminrep
    do l = lminrep, lmaxrep    ! Loop over Zernike's or Jacobi's functions
        call frad_fun(l, t)        ! Calls to zernike3DR or jacobiP depending on input
        do m = -l, l
            knt = knt + 1
            den = den + dot_product(omeganlm(0:kmaxrep,knt),radfunction) * ang(ind(l)+abs(m)+1) * zlm(knt)
            if (lgradient) then
                dxden = dxden +  ang(ind(l)+abs(m)+1) * dot_product(omeganlm(0:kmaxrep,knt),radfunction * zlmdx(knt) * rsinv &
                        + radderiv * zlm(knt) * aux )
                dyden = dyden +  ang(ind(l)+abs(m)+1) * dot_product(omeganlm(0:kmaxrep,knt),radfunction * zlmdy(knt) * rsinv &
                        + radderiv * zlm(knt) * bux )
                dzden = dzden +  ang(ind(l)+abs(m)+1) * dot_product(omeganlm(0:kmaxrep,knt),radfunction * zlmdz(knt) * rsinv &
                        + radderiv * zlm(knt) * cux )
            endif
        enddo
    enddo
    aux = 1.d0 / (rstar*sqrt(rstar))
    den = den * aux
    if (lgradient) then
        dxden = dxden * aux
        dyden = dyden * aux
        dzden = dzden * aux
    endif

! write(6,"(40(1H=))")
    return
    end
!
!    ***************************************************************
!
!     Computes the radial part of Zernike 3D functions divided by (r/r*)^l and its derivative with respecto to (r/r*)
! 
  subroutine zernike3DR(l, t2)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: k, kmax, l
    real(KREAL) :: t, t2
    gkl(-1) = 0.d0
    gkl(0) = 1.d0
    radfunction(0) = root(2*l+3)
    kmax = kmaxrep-1
    if (lechelon) kmax = kmax - l
    do k = 0, kmax
        gkl(k+1) = (cfgkl1((kmaxrep+1)*l+k+1) - cfgkl2((kmaxrep+1)*l+k+1) * t2) * gkl(k) &
                - cfgkl3((kmaxrep+1)*l+k+1) * gkl(k-1)
        radfunction(k+1) = akgkl(k+1) * root(2*l+4*(k+1)+3) * gkl(k+1)
    enddo
    if (lgradient) then
        t = sqrt(t2)
        dgkl(-1) = 0.d0
        dgkl(0) = 0.d0
        radderiv(0) = 0.d0
        do k = 0, kmax
            dgkl(k+1) = (cfgkl1((kmaxrep+1)*l+k+1) - cfgkl2((kmaxrep+1)*l+k+1) * t2) * dgkl(k) &
                    - cfgkl3((kmaxrep+1)*l+k+1) * dgkl(k-1) - 2.d0 * t * cfgkl2((kmaxrep+1)*l+k+1) * gkl(k)
            radderiv(k+1) = akgkl(k+1) * root(2*l+4*(k+1)+3) * dgkl(k+1)
        enddo
    endif
    return
    end
!
!**********************************************************************
! 
!     Computes the Jacobi polynomials P(0,2+2l)(2t-1)
!
  subroutine jacobiP(l, t)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: k, kmax, l
    real(KREAL) :: t
    radfunction(0) = 1.d0
    radfunction(1) = 4.d0 * t - 3.d0 + dble(l+l) * (t - 1.d0)
    kmax = kmaxrep-1
    if (lechelon) kmax = kmax - l
    do k = 1, kmax
        radfunction(k+1) = (dble(2*k+2*l+3) * (dble((k+l+1)*(k+l+2)) * (2.d0*t-1.d0) - dble((l+1)*(l+1))) &
                * radfunction(k) - dble(k*(k+l+2)*(k+2*l+2))*radfunction(k-1)) / dble((k+1)*(k+l+1)*(k+2*l+3))
    enddo
    if (lgradient) then
        radderiv(0) = 0.d0
        radderiv(1) = dble(l+l+4)
        do k = 1, kmax
            radderiv(k+1) = (dble(2*k+2*l+3) * (dble((k+l+1)*(k+l+2)) * (2.d0*t-1.d0) - dble((l+1)*(l+2))) &
                    * radderiv(k) - dble(k*(k+l+2)*(k+2*l+3))*radderiv(k-1)) / dble(k*(k+l+1)*(k+2*l+3))
        enddo
    endif
    do k = 0, kmax+1
        radfunction(k) = radfunction(k) * root(2*(k+l)+3)
    enddo
    if (lgradient) then
        do k = 0, kmax+1
            radderiv(k) = radderiv(k) * root(2*(k+l)+3)
        enddo
    endif
    return
    end
!
!**********************************************************************
! 
!     Computes the radial part of Jacobi P(0,2) polynomials
!
  subroutine jacobiP02(l, t)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: k, kmax, l
    real(KREAL) :: t
    radfunction(0) = 1.d0
    radfunction(1) = 4.d0 * t - 3.d0
    kmax = kmaxrep-1
    if (lechelon) kmax = kmax - l
    do k = 0, kmax
        radfunction(k+1) = (dble(2*k+3) * (dble((k+1)*(k+2)) * (2.d0*t-1.d0) - 1.d0) &
                * radfunction(k) - dble(k*(k+2)*(k+2))*radfunction(k-1)) / dble((k+1)*(k+1)*(k+3))
    enddo
    do k = 0, kmax+1
        radfunction(k) = radfunction(k) * root(2*k+3)
    enddo
    return
    end
!
!**********************************************************************
! 
  subroutine consta
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: i, ierr, k, knt, l, lm, m
    real(KREAL) :: aux, bux, sgn

    fact(0) = 1.d0
    facti(0) = 1.d0
    do i = 1, mxfact
        fact(i) = fact(i-1) * dble(i)               !  i!
        facti(i) = 1.d0 / fact(i)                    !  1.d0 / i!
    enddo

    allocate(ind(0:lmaxrep), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ind in consta. Stop')
    ind(0) = 0
    do i = 1, lmaxrep
        ind(i) = ind(i-1) + i         !  i*(i+1)/2
    enddo
    root(0) = 0.d0
    do i = 1, mxroot
        root(i) = sqrt(dble(i))        !  sqrt(i)
    enddo

    allocate(ang((lmaxrep+1)*(lmaxrep+2)/2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ang in consta. Stop')
!    ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
    ang(1) = 0.5d0 / sqrt(acos(-1.d0))
    lm = 1
    do l = 1, lmaxrep
        lm = lm + 1
        ang(lm) = ang(1) * sqrt(dble(2*l+1))
        aux = ang(lm) * root(2)
        do m = 1, l
            lm = lm + 1
            aux = aux / sqrt(dble(l-m+1)*dble(l+m))
            ang(lm) = aux
        enddo
    enddo

!     Computes and stores the coefficients for recursion of Zernike 3D functions
    m = (kmaxrep+1)*(lmaxrep+1)
    allocate(akgkl(0:kmaxrep), cfgkl1(m), cfgkl2(m), cfgkl3(m), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating akgkl, cfgkl1, cfgkl2 and cfgkl3. Stop')
    akgkl(0) = 1.d0        ! akgkl(k) = (-1)^k (k-1/2)! / (sqrt(pi) * k!)
    do k = 0, kmaxrep-1
            akgkl(k+1) = - dble(2*k+1) * akgkl(k) / dble(2*k+2)
    enddo
! write(6,"('akgkl = ', 8(1x,e17.10))") akgkl
    knt = 0
    do l = 0, lmaxrep
        do k = 0, kmaxrep
            knt = knt + 1
            cfgkl1(knt) = dble(4*k+2*l+3) * dble(4*k*(2*k+2*l+3)+4*l*(l+2)+3) &
                    / (dble(2*k+2*l+3) * dble(2*k+1) * dble(4*k+2*l+1))
            cfgkl2(knt) = dble(4*k+2*l+3) * dble(4*k+2*l+5) / (dble(2*k+2*l+3) * dble(2*k+1))
            cfgkl3(knt) = dble(4*k*k) * dble(4*k+2*l+5) * dble(2*k+2*l+1) &
                    / (dble(4*k+2*l+1) * dble(2*k+2*l+3) * dble(2*k-1) * dble(2*k+1))
        enddo
    enddo
! write(6,"('cfgkl1 = ', 8(1x,e17.10))") cfgkl1
! write(6,"('cfgkl2 = ', 8(1x,e17.10))") cfgkl2
! write(6,"('cfgkl3 = ', 8(1x,e17.10))") cfgkl3

    return
    end
!    
!   ***************************************************************
!
   subroutine grid3D
    USE DAMDENZJ_atdens
    implicit none
    character(8) :: tipo
    integer(KINT) :: i, ia, ierr, iuni, ix, iy, iz, nx, ny, nz
    real(KREAL) :: b2a, den, dxden, dyden, dzden, rx, ry, rz, x, y, z
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: arraysp(:)
#endif
    external zernike3DR, jacobiP
    rx = (xsup - xinf) / dltx + 0.5d0
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + 0.5d0
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + 0.5d0
    nz = int(rz) + 1

    write(6,"('3D GRID (inf,sup,dlt,npnts)')")
    write(6,"('x: ',3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
    write(6,"('y: ',3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
    write(6,"('z: ',3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz

    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        x1 = xinf * b2a
        y1 = yinf * b2a
        z1 = zinf * b2a
        x2 = (xinf+(nx-1)*dltx) * b2a
        y2 = (yinf+(ny-1)*dlty) * b2a
        z2 = (zinf+(nz-1)*dltz) * b2a
    else
        x1 = xinf
        y1 = yinf
        z1 = zinf
        x2 = (xinf+(nx-1)*dltx)
        y2 = (yinf+(ny-1)*dlty)
        z2 = (zinf+(nz-1)*dltz)
    endif

!    Opens file for density expansion tabulation
    allocate(array(nx), arraysp(nx), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating array and arraysp in grid3D. Stop')
    iuni = 21
    if (ljacobi) then
        tipo = "-jacobi"
    else
        tipo = "-zernike"
    endif
    call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d.plt")
!    Opens files for gradient tabulation
    if (lgradient) then
        allocate(arraydx(nx), arraydy(nx), arraydz(nx), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in gridpot. Stop')
        iuni = 23
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d-dx.pltd")
        iuni = 24
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d-dy.pltd")
        iuni = 25
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d-dz.pltd")
    endif
!    Grid tabulation
    do iz = 1, nz
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix = 1, nx
                x = xinf + (ix - 1) * dltx
                if (ljacobi) then
                    call compute_density(jacobiP, x, y, z, den, dxden, dyden, dzden)
                else
                    call compute_density(zernike3DR, x, y, z, den, dxden, dyden, dzden)
                endif
                array(ix) = den
                if (lgradient) then
                    arraydx(ix) = dxden
                    arraydy(ix) = dyden
                    arraydz(ix) = dzden
                endif
            enddo
            arraysp = array
            call linea(21, nx , arraysp )
            if (lgradient) then
                arraysp = arraydx
                call linea(23, nx , arraysp )
                arraysp = arraydy
                call linea(24, nx , arraysp )
                arraysp = arraydz
                call linea(25, nx , arraysp )
            endif
        enddo
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(array, arraysp)
    close(21)
    if (lgradient) then
        close(23)
        close(24)
        close(25)
        deallocate(arraydx, arraydy, arraydz)
    endif
    write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    return
    end
!    
!   ***************************************************************
!
   subroutine grid2D
    USE DAMDENZJ_atdens
    USE strings
    USE evaluate
    implicit none
    character(8) :: tipo
    integer(KINT) :: i, ia, iuni, iu, iv, nu, nv
    real(KREAL) :: b2a, den, dxden, dyden, dzden, ru, rv, u, v, x, y, z
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
#ifdef DBLPRCGRID
    real(KREAL) :: u1, u2, v1, v2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: u1, u2, v1, v2
    real(KREAL4), allocatable :: arraysp(:)
#endif
    external zernike3DR, jacobiP
    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1
    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1

    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv
    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        u1 = uinf * b2a
        v1 = vinf * b2a
        u2 = (uinf+(nu-1)*dltu) * b2a
        v2 = (vinf+(nv-1)*dltv) * b2a
    else
        u1 = uinf
        v1 = vinf
        u2 = (uinf+(nu-1)*dltu)
        v2 = (vinf+(nv-1)*dltv)
    endif

!    Opens file for electrostatic potential tabulation
    allocate(array(nu), arraysp(nu), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating array and arraysp in grid2D. Stop')
    iuni = 21
    if (ljacobi) then
        tipo = "-jacobi"
    else
        tipo = "-zernike"
    endif
    call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d.cnt")
!    Opens files for gradient tabulation
    if (lgradient) then
        allocate(arraydx(nu), arraydy(nu), arraydz(nu), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in gridpot2d. Stop')
        iuni = 23
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d-dx.cnt")
        iuni = 24
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d-dy.cnt")
        iuni = 25
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d-dz.cnt")
    endif

!    Grid tabulation
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)    ! Computes x(u,v)
            call evalexpr(y_func_uv,y)    ! Computes y(u,v)
            call evalexpr(z_func_uv,z)    ! Computes z(u,v)
            if (ljacobi) then
                call compute_density(jacobiP, x, y, z, den, dxden, dyden, dzden)
            else
                call compute_density(zernike3DR, x, y, z, den, dxden, dyden, dzden)
            endif
            array(iu) = den
            if (lgradient) then
                arraydx(iu) = dxden
                arraydy(iu) = dyden
                arraydz(iu) = dzden
            endif
        enddo
        arraysp = array
        call linea(21, nu , arraysp )
        if (lgradient) then
            arraysp = arraydx
            call linea(23, nu , arraysp )
            arraysp = arraydy
            call linea(24, nu , arraysp )
            arraysp = arraydz
            call linea(25, nu , arraysp )
        endif
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(arraysp)
    deallocate(array)
    close(21)
    if (lgradient) then
        deallocate(arraydx, arraydy, arraydz)
        close(23)
        close(24)
        close(25)
    endif
    write(6,"('Total number of tabulated points = ', i12)") nu*nv
    return
    end
!    
!   ***************************************************************
!
   subroutine grid3D_part
    USE DAMDENZJ_atdens
    implicit none
    character(8) :: tipo
    integer(KINT) :: i, ia, ierr, iuni, ix, iy, iz, nx, ny, nz
    real(KREAL) :: b2a, den, dxden, dyden, dzden, rx, ry, rz, x, y, z
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: arraysp(:)
#endif
    external zernike3DR, jacobiP
    rx = (xsup - xinf) / dltx + 0.5d0
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + 0.5d0
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + 0.5d0
    nz = int(rz) + 1

    write(6,"('3D GRID (inf,sup,dlt,npnts)')")
    write(6,"('x: ',3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
    write(6,"('y: ',3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
    write(6,"('z: ',3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz

    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        x1 = xinf * b2a
        y1 = yinf * b2a
        z1 = zinf * b2a
        x2 = (xinf+(nx-1)*dltx) * b2a
        y2 = (yinf+(ny-1)*dlty) * b2a
        z2 = (zinf+(nz-1)*dltz) * b2a
    else
        x1 = xinf
        y1 = yinf
        z1 = zinf
        x2 = (xinf+(nx-1)*dltx)
        y2 = (yinf+(ny-1)*dlty)
        z2 = (zinf+(nz-1)*dltz)
    endif

!    Opens file for density expansion tabulation
    allocate(array(nx), arraysp(nx), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating array and arraysp in grid3D. Stop')
    iuni = 21
    if (ljacobi) then
        tipo = "-jacobi"
    else
        tipo = "-zernike"
    endif
    call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d.plt")
!    Opens files for gradient tabulation
    if (lgradient) then
        allocate(arraydx(nx), arraydy(nx), arraydz(nx), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in gridpot. Stop')
        iuni = 23
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d-dx.pltd")
        iuni = 24
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d-dy.pltd")
        iuni = 25
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//trim(tipo)//"-d-dz.pltd")
    endif
!    Grid tabulation
    do iz = 1, nz
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix = 1, nx
                x = xinf + (ix - 1) * dltx
                if (ljacobi) then
                    call compute_density_part(jacobiP, x, y, z, den, dxden, dyden, dzden)
                else
                    call compute_density_part(zernike3DR, x, y, z, den, dxden, dyden, dzden)
                endif
                array(ix) = den
                if (lgradient) then
                    arraydx(ix) = dxden
                    arraydy(ix) = dyden
                    arraydz(ix) = dzden
                endif
            enddo
            arraysp = array
            call linea(21, nx , arraysp )
            if (lgradient) then
                arraysp = arraydx
                call linea(23, nx , arraysp )
                arraysp = arraydy
                call linea(24, nx , arraysp )
                arraysp = arraydz
                call linea(25, nx , arraysp )
            endif
        enddo
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(array, arraysp)
    close(21)
    if (lgradient) then
        close(23)
        close(24)
        close(25)
        deallocate(arraydx, arraydy, arraydz)
    endif
    write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    return
    end
!    
!   ***************************************************************
!
   subroutine grid2D_part
    USE DAMDENZJ_atdens
    USE strings
    USE evaluate
    implicit none
    character(8) :: tipo
    integer(KINT) :: i, ia, iuni, iu, iv, nu, nv
    real(KREAL) :: b2a, den, dxden, dyden, dzden, ru, rv, u, v, x, y, z
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
#ifdef DBLPRCGRID
    real(KREAL) :: u1, u2, v1, v2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: u1, u2, v1, v2
    real(KREAL4), allocatable :: arraysp(:)
#endif
    external zernike3DR, jacobiP
    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1
    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1

    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv
    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        u1 = uinf * b2a
        v1 = vinf * b2a
        u2 = (uinf+(nu-1)*dltu) * b2a
        v2 = (vinf+(nv-1)*dltv) * b2a
    else
        u1 = uinf
        v1 = vinf
        u2 = (uinf+(nu-1)*dltu)
        v2 = (vinf+(nv-1)*dltv)
    endif

!    Opens file for electrostatic potential tabulation
    allocate(array(nu), arraysp(nu), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating array and arraysp in grid2D. Stop')
    iuni = 21
    if (ljacobi) then
        tipo = "-jacobi"
    else
        tipo = "-zernike"
    endif
    call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d.cnt")
!    Opens files for gradient tabulation
    if (lgradient) then
        allocate(arraydx(nu), arraydy(nu), arraydz(nu), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in gridpot2d. Stop')
        iuni = 23
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d-dx.cnt")
        iuni = 24
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d-dy.cnt")
        iuni = 25
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-d-dz.cnt")
    endif

!    Grid tabulation
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)    ! Computes x(u,v)
            call evalexpr(y_func_uv,y)    ! Computes y(u,v)
            call evalexpr(z_func_uv,z)    ! Computes z(u,v)
            if (ljacobi) then
                call compute_density_part(jacobiP, x, y, z, den, dxden, dyden, dzden)
            else
                call compute_density_part(zernike3DR, x, y, z, den, dxden, dyden, dzden)
            endif
            array(iu) = den
            if (lgradient) then
                arraydx(iu) = dxden
                arraydy(iu) = dyden
                arraydz(iu) = dzden
            endif
        enddo
        arraysp = array
        call linea(21, nu , arraysp )
        if (lgradient) then
            arraysp = arraydx
            call linea(23, nu , arraysp )
            arraysp = arraydy
            call linea(24, nu , arraysp )
            arraysp = arraydz
            call linea(25, nu , arraysp )
        endif
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(arraysp)
    deallocate(array)
    close(21)
    if (lgradient) then
        deallocate(arraydx, arraydy, arraydz)
        close(23)
        close(24)
        close(25)
    endif
    write(6,"('Total number of tabulated points = ', i12)") nu*nv
    return
    end
!    
!   ***************************************************************
!
   subroutine compute_density_part(frad_fun, x, y, z, den, dxden, dyden, dzden)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: i, k, knt, l, laux, m, nu
    real(KREAL) :: rvec(0:kmaxrep)
    real(KREAL) :: aux, bux, cux, den, dxden, dyden, dzden, dosl1, r, r2, rsinv, t, x, xa, y, ya, z, za
    interface
        subroutine frad_fun(lrad, trad)
            integer(4) :: lrad
            real(8) :: trad
        end subroutine frad_fun
    end interface
    den = 0.d0
    dxden = 0.d0
    dyden = 0.d0
    dzden = 0.d0
    if ((x*x + y*y + z*z) .gt. rstar*rstar) return
    rsinv = 1.d0 / rstar
    xa = x * rsinv
    ya = y * rsinv
    za = z * rsinv
    r2 = xa*xa + ya*ya + za*za
    r = sqrt(r2)
! write(6,"('r2 = ', e17.10, ' r = ', e17.10)") r2, r
    zlm(1) = 1.d0        ! Regular spherical harmonics of r/rstar = (xa,ya,za)
    zlm(2) = ya
    zlm(3) = za
    zlm(4) = xa
    do l = 1, lmaxrep
        dosl1 = dble(l+l+1)
        zlm((l+1)*(l+3)+1) = dosl1 * (xa * zlm(l*(l+2)+1) - ya * zlm(l*l+1))        ! zlm(l+1,l+1,ia)
        zlm((l+1)*(l+1)+1) = dosl1 * (ya * zlm(l*(l+2)+1) + xa * zlm(l*l+1))        ! zlm(l+1,-(l+1),ia)
        zlm((l+2)*(l+2)-1) = dosl1 * za * zlm(l*(l+2)+1)                ! zlm(l+1,l,ia)
        zlm(l*(l+2)+3) = dosl1 * za * zlm(l*l+1)                    ! zlm(l+1,-l,ia)
        do m = 0, l-1
            aux = 1.d0 / dble(l-m+1)
            bux = dble(l+m)
            zlm((l+1)*(l+2)+m+1) = aux * (dosl1*za*zlm(l*(l+1)+m+1) - bux*r2*zlm((l-1)*l+m+1))    ! zlm(l+1,m,ia)
            zlm((l+1)*(l+2)-m+1) = aux * (dosl1*za*zlm(l*(l+1)-m+1) - bux*r2*zlm((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
        enddo
    enddo
    if (lgradient) then
        call derivzlm(lmaxrep)
    endif
    if (ljacobi) then
        t = r
    else
        t = r2
    endif
    if (r .gt. 0.d0) then
        aux = 1.d0 / (r * rstar)
        bux = ya * aux
        cux = za * aux
        aux = xa * aux
    else
        aux = 1.d0
        bux = 1.d0
        cux = 1.d0
    endif
    if (lindividk) then    !     Case of lindividk .eq. .true.: projection functions with given values of k index and all (l,m) are used
! write(6,"('Case of lindividk .eq. .true.: projection functions with given values of k index and all (l,m) are used')")
! write(6,"('k values: ', 50i3)") indicesv
        do l = lminrep, lmaxrep
            call frad_fun(l, t)
            knt = l*l
            do m = -l, l
                knt = knt + 1
                do i = 1, nindices
                    k = indicesv(i)
                    den = den + omeganlm(k,knt) * radfunction(k) * ang(ind(l)+abs(m)+1) * zlm(knt)
                    if (lgradient) then
                        dxden = dxden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdx(knt) * rsinv &
                                + radderiv(k) * zlm(knt) * aux )
                        dyden = dyden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdy(knt) * rsinv &
                                + radderiv(k) * zlm(knt) * bux )
                        dzden = dzden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdz(knt) * rsinv &
                                + radderiv(k) * zlm(knt) * cux )
                    endif
                enddo
            enddo
        enddo
    elseif (lindividl) then    !     Case of lindividl .eq. .true.: projection functions with given values of l, index and all (k,m) are used
! write(6,"('Case of lindividl .eq. .true.: projection functions with given values of l index and all (k,m) are used')")
! write(6,"('l values: ', 50i3)") indicesv
        do i = 1, nindices
            l = indicesv(i)
            if (l .lt. lminrep .or. l .gt. lmaxrep) cycle
            call frad_fun(l, t)
            knt = l*l
            do m = -l, l
                knt = knt + 1
                den = den + dot_product(omeganlm(0:kmaxrep,knt),radfunction) * ang(ind(l)+abs(m)+1) * zlm(knt)
                if (lgradient) then
                    dxden = dxden +  ang(ind(l)+abs(m)+1) * dot_product(omeganlm(0:kmaxrep,knt),radfunction * zlmdx(knt) * rsinv &
                            + radderiv * zlm(knt) * aux )
                    dyden = dyden +  ang(ind(l)+abs(m)+1) * dot_product(omeganlm(0:kmaxrep,knt),radfunction * zlmdy(knt) * rsinv &
                            + radderiv * zlm(knt) * bux )
                    dzden = dzden +  ang(ind(l)+abs(m)+1) * dot_product(omeganlm(0:kmaxrep,knt),radfunction * zlmdz(knt) * rsinv &
                            + radderiv * zlm(knt) * cux )
                endif
            enddo
        enddo
    elseif (lindividlk) then    !     Case of lindividl .eq. .true.: projection functions with given values of l, k indices and all m are used
! write(6,"('Case of lindividlk .eq. .true.: projection functions with given values of l, k indices and all m are used')")
! write(6,"('l,k values: ', 50i3)") indicesv
        if (mod(nindices,2) .ne. 0) then
            call error(1,'Wrong number of indices for projection functions of selected l and k. Two indices are required for&
                    & each function. Stop')
        endif
        laux = -1
        do i = 1, nindices, 2
            l = indicesv(i)
            if (l .lt. lminrep .or. l .gt. lmaxrep) cycle
            if (l .ne. laux) then
                call frad_fun(l, t)
                laux = l
            endif
            k = indicesv(i+1)
            knt = l*l
            do m = -l, l
                knt = knt + 1
                den = den + omeganlm(k,knt) * radfunction(k) * ang(ind(l)+abs(m)+1) * zlm(knt)
                if (lgradient) then
                    dxden = dxden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdx(knt) * rsinv &
                            + radderiv(k) * zlm(knt) * aux )
                    dyden = dyden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdy(knt) * rsinv &
                            + radderiv(k) * zlm(knt) * bux )
                    dzden = dzden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdz(knt) * rsinv &
                            + radderiv(k) * zlm(knt) * cux )
                endif
            enddo
        enddo
        return
    elseif (lindividlkm) then    !     Case of lindividlkm .eq. .true.: only projection functions with given values of k, l, m indices are used
! write(6,"('Case of lindividlkm .eq. .true.: projection functions with given values of k index and all (l,m) are used')")
! write(6,"('l,k,m values: ', 50i3)") indicesv
        if (mod(nindices,3) .ne. 0) then
                call error(1,'Wrong number of indices for individual projection functions. Three indices are required for&
                        & each function. Stop')
        endif
        laux = -1
        do i = 1, nindices, 3
            l = indicesv(i)
            if (l .lt. lminrep .or. l .gt. lmaxrep) cycle
            if (l .ne. laux) then
                call frad_fun(l, t)
                laux = l
            endif
            k = indicesv(i+1)
            knt = l*l
            do m = -l, l
                knt = knt + 1
                if (m .ne. indicesv(i+2)) cycle
                den = den + omeganlm(k,knt) * radfunction(k) * ang(ind(l)+abs(m)+1) * zlm(knt)
                if (lgradient) then
                    dxden = dxden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdx(knt) * rsinv &
                            + radderiv(k) * zlm(knt) * aux )
                    dyden = dyden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdy(knt) * rsinv &
                            + radderiv(k) * zlm(knt) * bux )
                    dzden = dzden + ang(ind(l)+abs(m)+1) * ( omeganlm(k,knt) * radfunction(k) * zlmdz(knt) * rsinv &
                            + radderiv(k) * zlm(knt) * cux )
                endif
            enddo
        enddo
        return
    endif
    aux = 1.d0 / (rstar*sqrt(rstar))
    den = den * aux
    if (lgradient) then
        dxden = dxden * aux
        dyden = dyden * aux
        dzden = dzden * aux
    endif

! write(6,"(40(1H=))")
    return
    end
!
!    ***************************************************************
!    Subroutine density_exact: tabulates density from its atomic expansion
!    
  subroutine density_exact
    USE DAMDENZJ_atdens
    implicit none
    integer :: ierr
    call leegeomatdensGTO
    allocate(zlm(max(16,(lmaxbase+1)*(lmaxbase+1))), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in density_exact. Stop')
    if (lgradient) then
        allocate(zlmdx(max(16,(lmaxbase+1)*(lmaxbase+1))), zlmdy(max(16,(lmaxbase+1)*(lmaxbase+1))), &
            zlmdz(max(16,(lmaxbase+1)*(lmaxbase+1))), stat = ierr)
        if (ierr .ne. 0) call error(1,'Memory error in densrepr when allocating zlmdx, zlmdy, zlmdz. Stop')
    endif
    if (.not. lgrid2d) then
        call grid_exac_3D
    else
        call grid_exac_2D
    endif
    return
    end
!
!    ***************************************************************
!
  subroutine leegeomatdensGTO
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT), parameter :: mxcoefs = 5*mxcap
    integer(KINT) :: i, ia, ierr, indnf, indng, ios, j, k, k1, k2, knt, nbas 
    integer(KINT), allocatable :: ngauss(:)
    real(KREAL) :: aux, bux, distmax, zntot
    real(KREAL) :: centcharge(3), cfaux(0:2*mxl)
!    Reads the number of centers
    open(15,file=trim(projectname)//".atomic_ggden",form='formatted', iostat = ierr)
    if (ierr .ne. 0) call error(ierr,'Error when opening file '//trim(projectname)//'.atomic_ggden. Stop')
    read(15,*) ncen
!    Allocates memory for geometry and density data
    ncaps = mxcap ! just for allocating

    allocate(atmnam(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating atmnam. Stop')

    allocate(coefs(mxcoefs), stat = ierr)
    if (.not. allocated(coefs)) call error(ierr,'Memory error when allocating coefs. Stop')

    allocate(ll(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ll. Stop')
    
    allocate(lmaxc(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating lmaxc. Stop')
    
    allocate(nf(ncaps), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating nf. Stop')
    
    allocate(ngauss(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ngauss. Stop')

    allocate(ngini(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ngini. Stop')

    allocate(ngfin(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ngfin. Stop')

    allocate(nzn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating nzn. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating rcen. Stop')

    allocate(xxg(mxcap), stat = ierr)
    if (.not. allocated(xxg)) call error(ierr,'Memory error when allocating xxg. Stop')

    allocate(zn(ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating zn. Stop')

    centcharge = 0.d0
    zntot = 0.d0
!    Reads the number of centers, geometry and basis set 
    do ia = 1, ncen
        read(15,*) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia), ngauss(ia)
        centcharge = centcharge + zn(ia) * rcen(:,ia)
        zntot = zntot + zn(ia)
        if (abs(zn(ia)-re(int(zn(ia) + umbrzn))) .gt. umbrzn) then
            nzn(ia) = 0
        else
            nzn(ia) = int(zn(ia) + umbrzn)
        endif
        atmnam(ia) = atmnms(nzn(ia))
    enddo
!     Transforms the coordinates to get the center of positive (nucelar) charges at the origin of coordinates
    centcharge = centcharge / zntot
    if (dot_product(centcharge,centcharge) .gt. 1.d-14) then
        write(6,"('Position of the center of positive charges: ',3(1x,e22.15))") centcharge
        write(6,"('Shifts the nuclear coordinates to put the center of positive charges at the origin')")
        do ia = 1, ncen
            rcen(1,ia) = rcen(1,ia) - centcharge(1)
            rcen(2,ia) = rcen(2,ia) - centcharge(2)
            rcen(3,ia) = rcen(3,ia) - centcharge(3)
        enddo
    endif
    if (lrstarrel) then
        distmax = 0.d0
        do ia = 1, ncen
            distmax = max(distmax,dot_product(rcen(:,ia),rcen(:,ia)))
        enddo
        rstar = rstar + sqrt(distmax)
    endif


!     Density data::
!         for each center:
!              number of Gaussians
!              for each Gaussian:
!                   l  quantum number
!                   exponent
!                   expansion coefficients of density for m = -l, -l+1, ...l
!     Arrays  nf, ngini  and  ngfin are initialized.
!
!    Contraction coefficients correspond to the expansion in UNNORMALIZED primitives

    indnf = 1
    indng = 1
    lmaxbase = 0

    lmaxc = 0
    ncaps = 0
    do ia = 1, ncen
        if (ngauss(ia) .le. 0) then
            ngini(ia) = -1
            ngfin(ia) = -1
            cycle
        endif
        ngini(ia) = indng
        ngfin(ia) = indng + ngauss(ia) - 1
        indng = indng + ngauss(ia)
        do j = 1, ngauss(ia)
            ncaps = ncaps + 1
            if (ncaps .gt. mxcap)  call error(1,'Error: maximum number of Gaussians for density expansion exceeded. Stop')
            read(15,*) ll(ncaps), xxg(ncaps), cfaux(0:2*ll(ncaps))
            if ((indnf + 2*ll(ncaps)+1) .gt. mxcoefs) &
                call error(1,'Error: maximum number of functions for density expansion exceeded. Stop')
            if (ll(ncaps) .gt. lmaxbase) lmaxbase = ll(ncaps)
            if (ll(ncaps) .gt. lmaxc(ia)) lmaxc(ia) = ll(ncaps)
            nf(ncaps) = indnf
            coefs(indnf:indnf+2*ll(ncaps)) = cfaux(0:2*ll(ncaps))
            indnf = indnf + 2*ll(ncaps) + 1
        enddo
    enddo
    return
    end
!
!    ***************************************************************
!
  subroutine grid_exac_3D
    USE DAMDENZJ_atdens
    implicit none
        character(8) :: tipo
    integer(KINT) :: i, ia, ierr, iuni, ix, iy, iz, nx, ny, nz
    real(KREAL) :: b2a, den, dxden, dyden, dzden, rx, ry, rz, x, y, z
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
#ifdef DBLPRCGRID
    real(KREAL) :: x1, x2, y1, y2, z1, z2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: x1, x2, y1, y2, z1, z2
    real(KREAL4), allocatable :: arraysp(:)
#endif
    
    rx = (xsup - xinf) / dltx + 0.5d0
    nx = int(rx) + 1
    ry = (ysup - yinf) / dlty + 0.5d0
    ny = int(ry) + 1
    rz = (zsup - zinf) / dltz + 0.5d0
    nz = int(rz) + 1

    write(6,"('3D GRID (inf,sup,dlt,npnts)')")
    write(6,"('x: ',3(2x,f12.5),2x,i4)") xinf, xsup, dltx, nx
    write(6,"('y: ',3(2x,f12.5),2x,i4)") yinf, ysup, dlty, ny
    write(6,"('z: ',3(2x,f12.5),2x,i4)") zinf, zsup, dltz, nz

    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        x1 = xinf * b2a
        y1 = yinf * b2a
        z1 = zinf * b2a
        x2 = (xinf+(nx-1)*dltx) * b2a
        y2 = (yinf+(ny-1)*dlty) * b2a
        z2 = (zinf+(nz-1)*dltz) * b2a
    else
        x1 = xinf
        y1 = yinf
        z1 = zinf
        x2 = (xinf+(nx-1)*dltx)
        y2 = (yinf+(ny-1)*dlty)
        z2 = (zinf+(nz-1)*dltz)
    endif

!    Opens file for density expansion tabulation
    allocate(array(nx), arraysp(nx), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating array and arraysp in grid3D. Stop')
    iuni = 21
    call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-atdens.plt")
!    Opens files for gradient tabulation
    if (lgradient) then
        allocate(arraydx(nx), arraydy(nx), arraydz(nx), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in gridpot. Stop')
        iuni = 23
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-atdens-dx.pltd")
        iuni = 24
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-atdens-dy.pltd")
        iuni = 25
        call cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, trim(filename)//"-atdens-dz.pltd")
    endif
!    Grid tabulation
    do iz = 1, nz
        z = zinf + (iz-1) * dltz
        do iy = 1, ny
            y = yinf + (iy - 1) * dlty
            do ix = 1, nx
                x = xinf + (ix - 1) * dltx
                call compute_density_exact(x, y, z, den, dxden, dyden, dzden)
                array(ix) = den
                if (lgradient) then
                    arraydx(ix) = dxden
                    arraydy(ix) = dyden
                    arraydz(ix) = dzden
                endif
            enddo
            arraysp = array
            call linea(21, nx , arraysp )
            if (lgradient) then
                arraysp = arraydx
                call linea(23, nx , arraysp )
                arraysp = arraydy
                call linea(24, nx , arraysp )
                arraysp = arraydz
                call linea(25, nx , arraysp )
            endif
        enddo
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(array, arraysp)
    close(21)
    if (lgradient) then
        close(23)
        close(24)
        close(25)
        deallocate(arraydx, arraydy, arraydz)
    endif
    write(6,"('Total number of tabulated points = ', i12)") nx*ny*nz
    return
    end
!    
!   ***************************************************************
!
   subroutine grid_exac_2D
    USE DAMDENZJ_atdens
    USE strings
    USE evaluate
    implicit none
    character(8) :: tipo
    integer(KINT) :: i, ia, iuni, iu, iv, nu, nv
    real(KREAL) :: b2a, den, dxden, dyden, dzden, ru, rv, u, v, x, y, z
    real(KREAL), allocatable :: array(:), arraydx(:), arraydy(:), arraydz(:)
#ifdef DBLPRCGRID
    real(KREAL) :: u1, u2, v1, v2
    real(KREAL), allocatable :: arraysp(:)
#else
    real(KREAL4) :: u1, u2, v1, v2
    real(KREAL4), allocatable :: arraysp(:)
#endif
    
    ru = (usup - uinf) / dltu + umed
    nu = int(ru) + 1
    rv = (vsup - vinf) / dltv + umed
    nv = int(rv) + 1

    write(6,"('2D GRID (inf,sup,dlt,npnts)')")
    write(6,"('u: ',3(2x,f12.5),2x,i4)") uinf, usup, dltu, nu
    write(6,"('v: ',3(2x,f12.5),2x,i4)") vinf, vsup, dltv, nv
    write(6,"(/'x(u,v) = ', a)") x_func_uv
    write(6,"('y(u,v) = ', a)") y_func_uv
    write(6,"('z(u,v) = ', a,/)") z_func_uv
    if (langstrom) then        ! Converts grid distances to angstrom
        b2a = 0.5291772d0
        u1 = uinf * b2a
        v1 = vinf * b2a
        u2 = (uinf+(nu-1)*dltu) * b2a
        v2 = (vinf+(nv-1)*dltv) * b2a
    else
        u1 = uinf
        v1 = vinf
        u2 = (uinf+(nu-1)*dltu)
        v2 = (vinf+(nv-1)*dltv)
    endif

!    Opens file for electrostatic potential tabulation
    allocate(array(nu), arraysp(nu), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating array and arraysp in grid2D. Stop')
    iuni = 21
    call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-atdens.cnt")
!    Opens files for gradient tabulation
    if (lgradient) then
        allocate(arraydx(nu), arraydy(nu), arraydz(nu), stat = ierr)
        if (ierr .ne. 0) call error(ierr,'Memory error when allocating arraydx, arraydy and arraydz in gridpot2d. Stop')
        iuni = 23
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-atdens-dx.cnt")
        iuni = 24
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-atdens-dy.cnt")
        iuni = 25
        call cabecera2d(nu, nv, u1, u2, v1, v2, iuni, trim(filename)//trim(tipo)//"-atdens-dz.cnt")
    endif

!    Grid tabulation
    do iv = 1, nv
        v = vinf + (iv - 1) * dltv
        do iu = 1, nu
            u = uinf + (iu - 1) * dltu
            call defparam('u',u)
            call defparam('v',v)
            call evalexpr(x_func_uv,x)    ! Computes x(u,v)
            call evalexpr(y_func_uv,y)    ! Computes y(u,v)
            call evalexpr(z_func_uv,z)    ! Computes z(u,v)
            call compute_density_exact(x, y, z, den, dxden, dyden, dzden)
            array(iu) = den
            if (lgradient) then
                arraydx(iu) = dxden
                arraydy(iu) = dyden
                arraydz(iu) = dzden
            endif
        enddo
        arraysp = array
        call linea(21, nu , arraysp )
        if (lgradient) then
            arraysp = arraydx
            call linea(23, nu , arraysp )
            arraysp = arraydy
            call linea(24, nu , arraysp )
            arraysp = arraydz
            call linea(25, nu , arraysp )
        endif
    enddo
!     Deallocates arrays and closes the grid files
    deallocate(arraysp)
    deallocate(array)
    close(21)
    if (lgradient) then
        deallocate(arraydx, arraydy, arraydz)
        close(23)
        close(24)
        close(25)
    endif
    write(6,"('Total number of tabulated points = ', i12)") nu*nv
    return
    end
!    
!   ***************************************************************
!
   subroutine compute_density_exact(x, y, z, den, dxden, dyden, dzden)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: ia, i1, l, la, k, knt, m, nfa
    real(KREAL) :: rvec(0:kmaxrep)
    real(KREAL) :: aux, bux, cux, den, dxden, dyden, dzden, dosl1, exa, r2, x, xa, y, ya, z, za
    den = 0.d0
    dxden = 0.d0
    dyden = 0.d0
    dzden = 0.d0
    do ia = 1, ncen
        xa = x - rcen(1,ia)
        ya = y - rcen(2,ia) 
        za = z - rcen(3,ia)
        r2 = xa*xa + ya*ya + za*za
! write(6,"('r2 = ', e17.10, ' r = ', e17.10)") r2, r
        zlm(1) = 1.d0        ! Regular spherical harmonics of r/rstar = (xa,ya,za)
        zlm(2) = ya
        zlm(3) = za
        zlm(4) = xa
        do l = 1, lmaxc(ia)-1
            dosl1 = dble(l+l+1)
            zlm((l+1)*(l+3)+1) = dosl1 * (xa * zlm(l*(l+2)+1) - ya * zlm(l*l+1))        ! zlm(l+1,l+1,ia)
            zlm((l+1)*(l+1)+1) = dosl1 * (ya * zlm(l*(l+2)+1) + xa * zlm(l*l+1))        ! zlm(l+1,-(l+1),ia)
            zlm((l+2)*(l+2)-1) = dosl1 * za * zlm(l*(l+2)+1)                ! zlm(l+1,l,ia)
            zlm(l*(l+2)+3) = dosl1 * za * zlm(l*l+1)                    ! zlm(l+1,-l,ia)
            do m = 0, l-1
                aux = 1.d0 / dble(l-m+1)
                bux = dble(l+m)
                zlm((l+1)*(l+2)+m+1) = aux * (dosl1*za*zlm(l*(l+1)+m+1) - bux*r2*zlm((l-1)*l+m+1))    ! zlm(l+1,m,ia)
                zlm((l+1)*(l+2)-m+1) = aux * (dosl1*za*zlm(l*(l+1)-m+1) - bux*r2*zlm((l-1)*l-m+1))    ! zlm(l+1,-m,ia)
            enddo
        enddo
        if (lgradient) then
            call derivzlm(lmaxc(ia))
        endif
        do i1 = ngini(ia), ngfin(ia)
            la = ll(i1)
            nfa = nf(i1)
            exa = xxg(i1)
            aux = exp(-exa*r2)
            do m = -la, la
                den = den + coefs(nfa+la+m) * aux * zlm(la*(la+1)+m+1)
            enddo
            if (lgradient) then
                do m = -la, la
                    dxden = dxden + coefs(nfa+la+m) * aux * (-dos * xa * exa * zlm(la*(la+1)+m+1) + zlmdx(la*(la+1)+m+1) )
                    dyden = dyden + coefs(nfa+la+m) * aux * (-dos * ya * exa * zlm(la*(la+1)+m+1) + zlmdy(la*(la+1)+m+1) )
                    dzden = dxden + coefs(nfa+la+m) * aux * (-dos * za * exa * zlm(la*(la+1)+m+1) + zlmdz(la*(la+1)+m+1) )
                enddo
            endif
        enddo
    enddo
! write(6,"(40(1H=))")
    return
    end
!
!    ***************************************************************
!    Subroutine cabecera: writes header for .plt files (binary)
!
  subroutine cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, s)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: i
    integer(KINT) :: nx,ny,nz,iuni,ns,iaux(0:2)
#ifdef DBLPRCGRID
    real(KREAL) :: x1,x2,y1,y2,z1,z2,v(0:5)
    integer(KINT), parameter :: i0 = 0
#else
    real(KREAL4) :: x1,x2,y1,y2,z1,z2,v(0:5)
    integer(KINT), parameter :: i0 = 3
#endif
    character*(*) :: s
!    If the compiler is other than INTEL's, uses the OPEN
!    sentence for stream files according to Fortran 2003 standard
#if _WIN32
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#elif __INTEL_COMPILER
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#else
    open (unit=iuni, file=s, form='unformatted', access='stream')
#endif
    iaux(0) = i0
    iaux(1) = 2
    write(iuni) iaux(0), iaux(1)
    iaux(0) = nz
    iaux(1) = ny
    iaux(2) = nx
    write(iuni) iaux(0), iaux(1), iaux(2)
    v(0) = z1
    v(1) = z2
    v(2) = y1
    v(3) = y2
    v(4) = x1
    v(5) = x2
    write(iuni) (v(i), i = 0, 5)
    end
!
!    ***************************************************************
!    Subroutine cabecera: writes head for .plt files (binary)
!
  subroutine cabecera2d(nu, nv, u1, u2, v1, v2, iuni, s)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: i
    integer(KINT) :: nu,nv,iuni,ns,iaux(0:1)
#ifdef DBLPRCGRID
    real(KREAL) :: u1,u2,v1,v2,v(0:3)
#else
    real(KREAL4) :: u1,u2,v1,v2,v(0:3)
#endif
    character*(*) :: s
!    If the compiler is other than INTEL's, uses the OPEN
!    sentence for stream files according to Fortran 2003 standard
#if _WIN32
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#elif __INTEL_COMPILER
    open (unit=iuni, file=s, form='binary', carriagecontrol='NONE')
#else
    open (unit=iuni, file=s, form='unformatted', access='stream')
#endif
    iaux(0) = nu
    iaux(1) = nv
    write(iuni) iaux(0), iaux(1)
    v(0) = u1
    v(1) = u2
    v(2) = v1
    v(3) = v2
    write(iuni) (v(i), i = 0, 3)
    END
!
!    Subroutine linea: writes lines of .plt files (binary)
!
  SUBROUTINE linea(iuni,n,v)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: i, iuni, n
#ifdef DBLPRCGRID
    real(KREAL) :: v(0:n-1)
#else
    real(KREAL4) :: v(0:n-1)
#endif
    write(iuni) (v(i), i = 0, n-1)
    END
    
!   ***************************************************************

  subroutine derivzlm(lmax)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: lmax, l, m
!    Derivatives of the regular harmonics with respecto to the Cartesian coordinates
    zlmdx(1) = 0.d0    ! Derivatives of the S spherical harmonics
    zlmdy(1) = 0.d0
    zlmdz(1) = 0.d0
    zlmdx(2) = 0.d0    ! Derivatives of the P spherical harmonics
    zlmdx(3) = 0.d0
    zlmdx(4) = zlm(1)
    zlmdy(2) = zlm(1)
    zlmdy(3) = 0.d0
    zlmdy(4) = 0.d0
    zlmdz(2) = 0.d0
    zlmdz(3) = zlm(1)
    zlmdz(4) = 0.d0
    zlmdx(5) = 6.d0 * zlm(2)    ! Derivatives of the D spherical harmonics
    zlmdx(6) = 0.d0
    zlmdx(7) = -zlm(4)
    zlmdx(8) = 3.d0 * zlm(3)
    zlmdx(9) = 6.d0 * zlm(4)
    zlmdy(5) = 6.d0 * zlm(4)
    zlmdy(6) = 3.d0 * zlm(3)
    zlmdy(7) = -zlm(2)
    zlmdy(8) = 0.d0
    zlmdy(9) = -(6) * zlm(2)
    zlmdz(5) = 0.d0
    zlmdz(6) = 3.d0 * zlm(2)
    zlmdz(7) = 2.d0 * zlm(3)
    zlmdz(8) = 3.d0 * zlm(4)
    zlmdz(9) = 0.d0
    zlmdx(10) = 15.d0 * zlm(5)        ! Derivatives of the F spherical harmonics
    zlmdx(11) = 10.d0 * zlm(6)
    zlmdx(12) = -0.5d0 * zlm(5)
    zlmdx(13) = -zlm(8)
    zlmdx(14) = 6.d0 * zlm(7) - 0.5d0 * zlm(9)
    zlmdx(15) = 10.d0 * zlm(8)
    zlmdx(16) = 15.d0 * zlm(9)
    zlmdy(10) = 15.d0 * zlm(9)
    zlmdy(11) = 10.d0 * zlm(8)
    zlmdy(12) = 6.d0 * zlm(7) + 0.5d0 * zlm(9)
    zlmdy(13) = -zlm(6)
    zlmdy(14) = -0.5d0 * zlm(5)
    zlmdy(15) = -10.d0 * zlm(6)
    zlmdy(16) = -15.d0 * zlm(5)
    zlmdz(10) = 0.d0
    zlmdz(11) = 5.d0 * zlm(5)
    zlmdz(12) = 4.d0 * zlm(6)
    zlmdz(13) = 3.d0 * zlm(7)
    zlmdz(14) = 4.d0 * zlm(8)
    zlmdz(15) = 5.d0 * zlm(9)
    zlmdz(16) = 0.d0
    do l = 4, lmax        ! Derivatives of the remaining spherical harmonics
        zlmdx(l*(l+1)+1) = - zlm((l-1)*l+2)
        zlmdy(l*(l+1)+1) = - zlm((l-1)*l)
        zlmdz(l*(l+1)+1) = dble(l) * zlm((l-1)*l+1)
        zlmdx(l*(l+1)+2) = 0.5d0 * (dble(l+1) * dble(l) * zlm((l-1)*l+1) - zlm((l-1)*l+3))
        zlmdx(l*(l+1)) = 0.5d0 * (- zlm((l-1)*l-1))
        zlmdy(l*(l+1)+2) = -0.5d0 * (zlm((l-1)*l-1))
        zlmdy(l*(l+1)) = 0.5d0 * (dble(l+1) * dble(l) * zlm((l-1)*l+1) + zlm((l-1)*l+3))
        zlmdz(l*(l+1)+2) = dble(l+1) * zlm((l-1)*l+2)
        zlmdz(l*(l+1)) = dble(l+1) * zlm((l-1)*l)
        do m = 2, l-2
            zlmdx(l*(l+1)+m+1) = 0.5d0 * (dble(l+m) * dble(l+m-1) * zlm((l-1)*l+m) - zlm((l-1)*l+m+2))
            zlmdx(l*(l+1)-m+1) = 0.5d0 * (dble(l+m) * dble(l+m-1) * zlm((l-1)*l-m+2) - zlm((l-1)*l-m))
            zlmdy(l*(l+1)+m+1) = -0.5d0 * (dble(l+m) * dble(l+m-1) * zlm((l-1)*l-m+2) + zlm((l-1)*l-m))
            zlmdy(l*(l+1)-m+1) = 0.5d0 * (dble(l+m) * dble(l+m-1) * zlm((l-1)*l+m) + zlm((l-1)*l+m+2))
            zlmdz(l*(l+1)+m+1) = dble(l+m) * zlm((l-1)*l+m+1)
            zlmdz(l*(l+1)-m+1) = dble(l+m) * zlm((l-1)*l-m+1)
        enddo
        zlmdx(l*(l+2)) = dble(l+l-1) * dble(l-1) * zlm(l*l-1)
        zlmdy(l*(l+2)) = -dble(l+l-1) * dble(l-1) * zlm(l*(l-2)+3)
        zlmdz(l*(l+2)) = dble(l+l-1) * zlm(l*l)
        zlmdx(l*l+2) = dble(l+l-1) * dble(l-1) * zlm(l*(l-2)+3)
        zlmdy(l*l+2) = dble(l+l-1) * dble(l-1) * zlm(l*l-1)
        zlmdz(l*l+2) = dble(l+l-1) * zlm(l*(l-2)+2)
        zlmdx((l+1)*(l+1)) = dble(l) * dble(l+l-1) * zlm(l*l)
        zlmdy((l+1)*(l+1)) = -dble(l) * dble(l+l-1) * zlm((l-2)*l+2)
        zlmdz((l+1)*(l+1)) = 0.d0
        zlmdx(l*l+1) = dble(l) * dble(l+l-1) * zlm((l-2)*l+2)
        zlmdy(l*l+1) = dble(l) * dble(l+l-1) * zlm(l*l)
        zlmdz(l*l+1) = 0.d0
    enddo
    return
    end
!
!    -------------------------------------------------------------------------------------------------------
!
  subroutine error(ierr, msg)
    USE DAMDENZJ_atdens
    implicit none
    integer(KINT) :: ierr
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    stop
    end
