!  Copyright 2013-2018, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
!	Modules for programs Modules for programs TDAMCOM16_mpi.F90, TDAMDEN16_mpi.F90, 
!           TDAMPOT16_mpi.F90, TDAMTOPO16_mpi.F90
!
! Version of April 2019
!> @file TDAMGLOBAL15_mpi.F90

!> @author Anmol Kumar  (modified by Rafael Lopez in September 2017)
!> @details
!> The driver file of topography code utilizes the module daminitial_t.
!> Variables for other modules have been defined herein.
!> @param [in] KINT Decides 4-byte integer. Moved from DAM320_T    
!> @param [in] KINT8 Decides 8-byte integer. Moved from DAM320_T    
!> @param [in] KREAL Decides double precision floating point number. Moved from DAM320_T     
!> @param [in] KREAL4 Decides single precision floating point number. Moved from DAM320_T     
!> @param atmnms atomic symbols of all element in periodic table
!> @param atmwts atomic mass of all element in periodic table
!> @param iout Stores number 151 for output log file.
!> @param cox,coy,coz x,y,z coordinate of Center of mass
!> @param drcutcp Value of gradient to decide if the path is near to a CP.
!> @param xmpi,ympi,zmpi Input coordinate for abstract point at which field value has 
!> @param mesp logical used for enabling mesp topography calculation 
!> @param med logical used for enabling med topography calculation 
!> @param topomanual logical used for finding field value at 
!! provided points.Used for checking bug in the code in case 
!! any change is made in potential and density evaluation.
!> @param topograph logical used for mapping critical points
!> @param gradpath logical used for gradient path evaluation
!> @param basin logical used for basin generation
!> @param wireframe Decides the type of display for basin. Other option is Solidsurf.
!> @param addguesspt logical which allows
!! the user to provide extra guess points. 
!! These guess points will be appended with the program generated guess points.
!> @param lpoints logical which decides the field evaluation on discrete points.
!! Moved from DAMDEN320_T and DAMPOT320_T module. Other such logicals are
!! lderiv2, lgradient, largo, lexact and langstrom.
!> @param cnvg Decides the convergence criterion for optimizer of guess points to
!! decide if it has found the critical point.
!> @param boxl Puts simple box constraint around each guess point, within which
!! the optimization will be performed around a guess point. 
!> @param lmaxi l value used for expansion in spherical harmonics in a given calculation. 
!> @param dirsep character for directory names separator: "/" for unix, "\\" for MS-windows
!> @param iswindows .true. if running on a MS-windows system
!> @param Projectname Name of project used for input files
!> @param filename Name of project used for output files
!> @param rcen Cartesian coordinates of atomic centers
!> @param nzn atomic number of atomic centers
!> @param ncen number of atomic centers
!> @param atwt atomic mass of centers
!> @param atmnam atomic symbols of centers

MODULE DAMINITIAL_T
    IMPLICIT NONE
    INTEGER(4), PARAMETER    :: KINT = 4, KREAL = 8, KREAL4 = 4, KINT8 = 8, mxguess = 200
    INTEGER(KINT), PARAMETER :: IOUT=6   ! Modified by Rafa
    REAL(KREAL), PARAMETER   :: ZERO=0.d0,P_I=3.14159265359D0
    REAL(KREAL), PARAMETER   :: ANGAU=1.889725989d0, AUANG=0.529177249d0
    INTEGER(KINT):: LMAXI, NCEN, ncntguess
    INTEGER(KINT), ALLOCATABLE :: NZN(:)
    REAL(KREAL) :: COX,COY,COZ,TMASS
    REAL(KREAL) :: FDISP
    REAL(KREAL) :: CNVG, BOXL, BOXG , BOXB, BOXT,STEPSZT, ANGLE, EXLN, drcutcp
    REAL(KREAL) :: rcntguess(3,mxguess)
    REAL(KREAL), ALLOCATABLE   :: RCEN(:,:)
    REAL(KREAL), ALLOCATABLE   :: ATWT(:)
    LOGICAL:: MESP, MED, TOPOGRAPH, GRADPATH, BASIN
    LOGICAL:: WIREFRAME, SOLIDSURF,ADDGUESS,EXDRAW, iswindows
    LOGICAL:: LDERIV2,LGRADIENT,LARGO,LEXACT
    CHARACTER      :: DIRSEP
    CHARACTER(300) :: PROJECTNAME,GUESSFILE
    CHARACTER(256) :: FILENAME
    CHARACTER(2), allocatable :: ATMNAM(:)
    CHARACTER(2)   :: ATMNMS(0:103) = (/                             &
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
    REAL(KREAL) :: ATMWTS(0:103) = (/            &
                  0.0, 1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.011 &
        , 14.00674,	15.9994, 18.9984032, 20.1797, 22.989768, 24.3050 &
        , 26.981539, 28.0855, 30.97362,	32.066,	35.4527, 39.948, 39.0983 &
        , 40.078, 44.955910, 47.88, 50.9415, 51.9961, 54.93085, 55.847 &
        , 58.93320, 58.69, 63.546, 65.39, 69.723, 72.61, 74.92159, 78.96 &
        , 79.904, 83.80, 85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.94 &
        , 98.9063, 101.07, 102.90550, 106.42, 107.8682, 112.411, 114.82 &
        , 118.710, 121.75, 127.60, 126.90447, 131.29, 132.90543, 137.327 &
        , 138.9055, 140.115, 140.90765, 144.24, 146.9151, 150.36, 151.965 &
        , 157.25, 158.92534, 162.50, 164.93032, 167.26, 168.93421, 173.04 &
        , 174.967, 178.49, 180.9479, 183.85, 186.207, 190.2, 192.22, 195.08 &
        , 196.96654, 200.59, 204.3833, 207.2, 208.98037, 208.9824, 209.9871 &
        , 222.0176, 223.0197, 226.0254, 227.0278, 232.0381, 231.0359, 238.0289 &
        , 237.0482, 244.0642, 243.0614, 247.0703, 247.0703, 251.0796, 252.0829 &
        , 257.0951, 258.0986, 259.1009, 260.1053 /)
END MODULE

!	
!===============================================================================================
!                 MODULE DAM320_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @note The following fundamental parameters set the dimensions:
!> @param mxcap highest number of shells allowable in basis set
!> @param mxcen highest number of centers
!> @param mxl   highest l value of STO or GTO basis set
!> @param mxn   highest n value of STO
!> @param mxlexp highest l of expansions in spherical harmonics
!> @param mxmult highest order of molecular multipoles
!> @param mxsigexp highest expansion order for the B integrals in frgsigma
!> @param mxrtab  highest number of tabulation points for density and potential
!> @note The following axiliary parameters derive from the previous ones
!> @param mxldst highest l of distributions
MODULE DAM320_T
    USE DAMINITIAL_T, ONLY: KINT, KREAL, KREAL4, KINT8, IOUT, dirsep, iswindows
    IMPLICIT NONE
    integer(KINT), parameter ::  mxcen = 500, mxl = 6, mxn = 9, mxcap = 8000
    integer(KINT), parameter ::  mxlexp = 30, mxmult = 4, mxrtab = 200
    integer(KINT), parameter ::  mxsigexp = 100
    integer(KINT), parameter ::  mxldst = mxl+mxl, mxltot = mxldst + mxlexp
    integer(KINT), parameter :: nintervaj = 29	! Number of intervals for piecewise representation of the radial factors
!		Radii of fitting intervals
    real(KREAL), parameter :: rinterv(0:nintervaj) = (/0.d0, 1.d-2, 3.d-2, 5.d-2, .7d-1, 1.d-1, 2.d-1, 3.d-1, 4.d-1, 5.d-1, &
            6.d-1, 7.d-1, 8.d-1, 9.d-1, 1.d0, 1.5d0, 2.d0, 2.5d0, 3.d0, 4.d0, &
            5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 12.d0, 14.d0, 16.d0, 20.d0/)
    integer(KINT), parameter :: mxlenpol = 20	! Highest number of fitting polynomials
!   integer(KINT), parameter :: nintervaj = 21
!   real(KREAL) :: rinterv(0:nintervaj) = (/0.d0,1.d-2, 3.d-2, 5.d-2, 1.d-1, 2.d-1, 4.d-1, 7.d-1, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, &
!           6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 12.d0, 14.d0, 16.d0, 20.d0/)
! 		integer(KINT), parameter :: mxlenpol = 10
    contains
    function atom_core(x)result(q)
       implicit none
       real(KREAL) :: x, q
       if (x .lt. 3.d0) then
           q = x
       else if (x .lt. 11.d0) then
           q = x-2.d0
       else if (x .lt. 19.d0) then
           q = x-10.d0
       else if (x .lt. 37.d0) then
           q = x-18.d0
       else if (x .lt. 55.d0) then
           q = x-36.d0
       else if (x .lt. 87.d0) then
           q = x-54.d0
       else
           q = x-86.d0
       endif
    end function atom_core
END MODULE
!
!                 END OF MODULE DAM320_T
!...............................................................................................
!===============================================================================================
!                 MODULE DAM320_DATA_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @param umbral threshold for neglecting radial factors
!> @param xx basis set STO exponents 
!> @param zn atomic number
!> @param rnor basis set STO radial normalization factors
!> @param lmaxexp length of expansion on l
!> @param nbas 	number of basis functions 
!> @param ncaps	number of shells
!> @param ioptaj fitting option
!> @param nn n quantum numbers of STO basis set 
!> @param ll l quantum numbers of STO basis set
!> @param nf absolute order index of first function in a shell
!> @param ngini	absolute order index of the first shell of a center
!> @param ngfin absolute order index of the last shell of a center
!> @param lastfnc index of last function of a given center
!> @param lmaxc	highest l quantum number in a center
!> @param lmaxbase	highest l quantum number in the basis set
!> @param lsdisf if .false., negligible charge distribution
!> @param longoutput it .true. more detailed output 
!> @param lsto	if .true. Slater basis set, otherwise, Gaussian basis set
!> @param dmat	full density matrix
!> @param xajust exponential (nonlinear) fitting parameters
!> @param cfajust linear fitting parameters
!> @param cfrint1 auxiliary integrals for potential and field 
!> @param cfrint2l2  auxiliary integrals for potential and field
!> @param nintervaj number of fitting intervals of fragments
!> @param nexpaj	number of intervals with exponential term
!> @param lenpol   number of terms in the polynomial
!> @param ipmax	number of terms in the expansion of the B integrals in frgsigma
!> @param ncheb	array with different numbers of tabulation points per center
!> @param fct		auxiliary factor for locating the interval where a point is placed
!> @param wthreshold threshold for contributions for the B integrals (see subroutine Avk)
MODULE DAM320_DATA_T
    USE DAMINITIAL_T
    USE DAM320_T
    IMPLICIT NONE
    integer(KINT), parameter :: mxtpw = 2*mxlexp+mxlenpol+3
    integer(KINT), parameter :: npntintr = (mxlenpol*3+1) / 2
    integer(KINT), parameter :: npntaj = npntintr * nintervaj
    real(KREAL) :: umbral
    integer(KINT) :: nintstd, nfitpar, lmaxexp, nbas, ncaps, ioptaj, lmaxbase, lmtop
    logical :: longoutput, lsto, lgbsgz, ldengz, lvalence
    real(KREAL) :: roblk(-mxl:mxl,-mxl:mxl)
    logical*1, allocatable :: lsdisf(:,:)
    real(KREAL), allocatable :: dmat(:,:), dmataux(:,:), rnor(:), xx(:), zn(:)
    integer(KINT), allocatable :: ll(:), lmaxc(:), nexpaj(:), ngini(:), ngfin(:), nn(:), nf(:)    !, nzn(:)
    real(KREAL), allocatable :: dl(:,:,:), expvec(:), expvinv(:), f0(:,:), fa(:,:), ftab(:,:), ha(:,:), r2pow(:,:), radio(:,:)
    real(KREAL), allocatable :: residuals(:), rl(:,:,:), rlt(:), rpntaj(:), rpow(:,:), xajust(:)
    real(KREAL), allocatable :: Qgpart(:), qppart(:), rintr1(:), rintr2l2(:), vaux1c(:,:), ymat1(:,:), ymat2(:,:)
    real(KREAL), allocatable :: dosxcheb(:), fkvec0(:), tchvec0(:), tchvec0a(:), tchvec1(:), xcheb(:)
    real(KREAL), allocatable, target, dimension(:) :: cfajust, cfrint1, cfrint2l2, cfajustrank, cfrint1rank, cfrint2l2rank
    real(KREAL), allocatable, target, dimension(:,:) :: rmultip, rmultipfr
    integer(KINT), allocatable, dimension(:) :: indintrv
    integer(KINT), allocatable, dimension(:) :: icfpos, icfposrank
    integer(KINT) :: ipmax, lencfparank, lenexpaj, natomtype, lenxajust
    real(KREAL) :: fct, wthreshold
! 		indices and arrays for multipoles calculation
    integer(KINT) :: ilow, iupp
!		auxiliary arrays for multipolar moments of STO distributions
    real(KREAL), allocatable :: bkmat(:,:), qlm2c(:), qlmdst(:), qlmasint(:), powu(:,:), pow1mu(:,:)
    real(KREAL) :: bkv1(-2*mxn+1:mxltot), bkv2(-2*mxn+1:mxltot), cina(0:mxn), cinb(0:mxn), powR(0:mxlexp), pow2(0:mxldst) &
            , powzamzb(0:2*mxn+mxlexp+1), powzamzbi(0:2*mxn+mxlexp+1), auxu(0:mxn+mxl), scomp(0:mxn+2*mxl)
!		auxiliary arrays for multipolar moments of CGTO distributions
    real(KREAL), allocatable :: besselint(:,:)
    integer(KINT) :: kntintB, kmltan, kmltquad
    real(KREAL), allocatable :: av(:), bv(:), sol(:)	! for overlap integrals with ellipsoidal coords
END MODULE
!
!                 END OF MODULE DAM320_DATA_T
!...............................................................................................
!===============================================================================================
!                 MODULE DAM320_CONST_T
!=============================================================================================== 
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @note Coefficients and m values for the decomposition of the product of two Phi_m(phi) functions (of the real spherical harmonics)
!> Coefficients of Chebyshev T polynomials
!> Coefficients for the decomposition of products of pairs of real spherical harmonics
MODULE DAM320_CONST_T
    USE DAMINITIAL_T
    USE DAM320_T
    IMPLICIT NONE
    integer(KINT), parameter :: mxroot = 80, mxreal = 1000, mxfact = 150, mxangl = mxlexp, mxbin = max(mxltot,mxlenpol)
    integer(KINT), parameter :: mxangm = mxangl , mxind = (mxltot+1)*(mxltot+2)/2
    integer(KINT), parameter :: mxemes = mxltot, mxtot = mxlexp, mxlcof = mxltot*(mxltot+3)/2, mxkcof = mxlcof*(mxlcof+3)/2
    real(KREAL), parameter :: cero = 0.d0, uno = 1.d0, dos = 2.d0, cuatro = 4.d0, ocho = 8.d0, r16 = 16.d0, udec = 0.1d0
    real(KREAL), parameter :: umed = .5d0, pt25 = .25d0, euler = 5.7721566490153280D-1, raiz2 = 1.4142135623730950d0
    real(KREAL), parameter :: umbrzn = 1.d-6
!		re(i) = double(i)
! 		ri(i) = double(1/i)
!		fact(i) = double(i!)
!		facti(i) = double(1/i!)
!		dosl1(i) = double(2l+1))
!		dosl1i(i) = 1 / double(2l+1))
    real(KREAL) :: re(-mxreal:mxreal), ri(-mxreal:mxreal), umedpow(0:mxlexp), fact(0:mxfact), facti(0:mxfact) &
            , facts(-1:mxfact), dosl1(-mxreal:mxreal), dosl1i(-mxreal:mxreal), bin((mxbin+1)*(mxbin+2)/2)
    real(KREAL) :: pi, raizpi, pimed
!		ind(i) = i*(i+1) / 2 
    integer(KINT) :: ind(0:mxind)
!		root(i) = double( sqrt(i) )
!		rooti(i) = double( 1 / sqrt(i) )
    real(KREAL) :: root(0:mxroot), rooti(mxroot)
!    	ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
    real(KREAL) :: ang((mxangl+1)*(mxangl+2)/2)
!		Coefficients and m values for the decomposition of the product of two Phi_m(phi) functions (of the real spherical harmonics)
    real(KREAL) :: ssv(-mxemes:mxemes,-mxemes:mxemes), sdv(-mxemes:mxemes,-mxemes:mxemes)
    integer(KINT) :: msv(-mxemes:mxemes,-mxemes:mxemes), mdv(-mxemes:mxemes,-mxemes:mxemes)
    integer(KINT) :: indk12((mxltot+1)**2,(mxltot+1)**2)
!		Coefficients of Chebyshev T polynomials
    real(KREAL) ::  chebTcf((mxlenpol/2+1)*((mxlenpol+1)/2+1))
!		Coefficients for the decomposition of products of pairs of real spherical harmonics
    real(KREAL), allocatable :: app(:,:), bpp(:,:), ccl1l2(:)
    integer(KINT), allocatable :: i1l1l2(:,:), i2l1l2(:,:), lml1l2(:), npl1l2(:), llm(:), mlm(:), lll1l2(:), mml1l2(:)
    integer(KINT) :: idimlml1l2
!		Hypergeometrics 3F2(-j,j+2,1/2;3/2,3/2;1)
    real(KREAL),parameter :: h3f2(0:30) =	(/ 1.0000000000000000D0, 3.3333333333333330D-1, 2.8888888888888880D-1   &
    , 1.8095238095238100D-1, 1.6698412698412700D-1, 1.2400192400192400D-1, 1.1727637441923160D-1, 9.4283494283494300D-2 &
    , 9.0343498186635500D-2, 7.6045990473235050D-2, 7.3461722941036200D-2, 6.3716724290152750D-2 &
    , 6.1892360883217920D-2, 5.4825975317485430D-2, 5.3469760871032360D-2, 4.8111771784334780D-2 &
    , 4.7064198756022750D-2, 4.2862219460053240D-2, 4.2028788279367650D-2, 3.8645297583347980D-2 &
    , 3.7966485503304700D-2, 3.5183653709814110D-2, 3.4620113210353590D-2, 3.2291083670560490D-2 &
    , 3.1815766854350310D-2, 2.9837935700827030D-2, 2.9431638990733510D-2, 2.7731158377428100D-2 &
    , 2.7379872214383080D-2, 2.5902238055824550D-2, 2.5595502719280240D-2 /)
! for overlap integrals with ellipsoidal coords
    integer(KINT), allocatable :: ipntap(:,:), ipntalfa(:,:)
    real(KREAL), allocatable :: ap(:), alfasol(:), alm(:)
    real(KREAL), parameter :: alfacutoff = 1.d-20
END MODULE
!
!                 END OF MODULE DAM320_CONST_T
!...............................................................................................
!===============================================================================================
!                 MODULE DAM320_LEGENDRE_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @note Abscisae and weights for Legendre quadrature rules with 25 and 35 points
!> Data generated with notebook legendre_quadrature.nb
MODULE DAM320_LEGENDRE_T
    USE DAM320_T
!		Abscisae and weights for Legendre quadrature rules with 25 and 35 points
    integer(KINT) :: nqleg, nqleg2
    integer(KINT), parameter :: mxleg = 35
    real(KREAL), pointer :: u(:), w(:), u2(:), w2(:)
!	Data generated with notebook legendre_quadrature.nb
    real(KREAL), target :: uleg25(25) = (/ 2.2215151047509510D-3,1.1668039270241240D-2,2.8512714385512830D-2  &
            ,5.2504001060862320D-2,8.3278685619583000D-2,1.2037036848132120D-1,1.6321681576326580D-1          &
            ,2.1116853487938850D-1,2.6349863427714250D-1,3.1941384709530610D-1,3.7806655813950580D-1          &
            ,4.3856765369464480D-1,5.0000000000000000D-1,5.6143234630535520D-1,6.2193344186049420D-1          &
            ,6.8058615290469390D-1,7.3650136572285750D-1,7.8883146512061100D-1,8.3678318423673400D-1          &
            ,8.7962963151867900D-1,9.1672131438041700D-1,9.4749599893913800D-1,9.7148728561448700D-1          &
            ,9.8833196072975900D-1,9.9777848489524900D-1/)
    real(KREAL), target :: wleg25(25) = (/ 5.6968992505131440D-3,1.3177493307516069D-2,2.0469578350653156D-2  &
            ,2.7452347987917596D-2,3.4019166906178459D-2,4.0070350167500509D-2,4.5514130991481825D-2          &
            ,5.0267974533525322D-2,5.4259812237131827D-2,5.7429129572855824D-2,5.9727881767892386D-2          &
            ,6.1121221495155021D-2,6.1588026863357726D-2,6.1121221495155021D-2,5.9727881767892386D-2          &
            ,5.7429129572855824D-2,5.4259812237131827D-2,5.0267974533525322D-2,4.5514130991481825D-2          &
            ,4.0070350167500509D-2,3.4019166906178459D-2,2.7452347987917596D-2,2.0469578350653156D-2          &
            ,1.3177493307516069D-2,5.6968992505131440D-3/)
    real(KREAL), target :: uleg35(35) = (/ 1.1467154501998510D-3,6.0321177780742500D-3,1.4781191980385080D-2  &
            ,2.7327425896086340D-2,4.3572869320341200D-2,6.3390437487388800D-2,8.6625050453887300D-2          &
            ,1.1309487385654370D-1,1.4259274922168560D-1,1.7488781766705480D-1,2.0972732762511770D-1          &
            ,2.4683861337925570D-1,2.8593122924109290D-1,3.2669922278459300D-1,3.6882352939535200D-1          &
            ,4.1197446941700500D-1,4.5581432836217000D-1,5.0000000000000000D-1,5.4418567163783000D-1          &
            ,5.8802553058299500D-1,6.3117647060464800D-1,6.7330077721540700D-1,7.1406877075890700D-1          &
            ,7.5316138662074400D-1,7.9027267237488200D-1,8.2511218233294500D-1,8.5740725077831400D-1          &
            ,8.8690512614345600D-1,9.1337494954611300D-1,9.3660956251261100D-1,9.5642713067965900D-1          &
            ,9.7267257410391400D-1,9.8521880801961500D-1,9.9396788222192600D-1,9.9885328454980000D-1/)
    real(KREAL), target :: wleg35(35) = (/ 2.9417167102215425D-3,6.8254141741807461D-3,1.0661489955741790D-2  &
            ,1.4414630054447127D-2,1.8055057931731690D-2,2.1554211163085109D-2,2.4884685200676765D-2          &
            ,2.8020408106185064D-2,3.0936835983040094D-2,3.3611142634543452D-2,3.6022397386280032D-2          &
            ,3.8151728577721027D-2,3.9982471121162131D-2,4.1500296864428294D-2,4.2693326696049563D-2          &
            ,4.3552223498591767D-2,4.4070265215137731D-2,4.4243397453552145D-2,4.4070265215137731D-2          &
            ,4.3552223498591767D-2,4.2693326696049563D-2,4.1500296864428294D-2,3.9982471121162131D-2          &
            ,3.8151728577721027D-2,3.6022397386280032D-2,3.3611142634543452D-2,3.0936835983040094D-2          &
            ,2.8020408106185064D-2,2.4884685200676765D-2,2.1554211163085109D-2,1.8055057931731690D-2          &
            ,1.4414630054447127D-2,1.0661489955741790D-2,6.8254141741807461D-3,2.9417167102215425D-3/)
END MODULE
!
!                 END OF MODULE DAM320_LEGENDRE_T
!...............................................................................................
!===============================================================================================
!                 MODULE DAMDEN320_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @details Contains variables used by TDAMDEN320_mpi.F90
MODULE DAMDEN320_T
    USE DAMINITIAL_T
    USE DAM320_T
    IMPLICIT NONE
    integer(KINT), parameter :: mxsel = 200
    logical :: lmolec, latomics, latomsel, ldensacc, lboundsx, lboundsy, lboundsz, laplacian
    logical, allocatable :: lselat(:), lnegia(:)
    real(KREAL) :: umbrlargo, rtab(3,mxrtab), xboundinf, xboundsup, yboundinf, yboundsup, zboundinf, zboundsup
    real(KREAL), allocatable :: xajustd(:,:), rlargo(:)
    integer(KINT) :: iatomsel(mxsel), lminrep, lmaxrep, nsel, numrtab
    integer(KINT), allocatable :: icfposd(:,:), lcorto(:,:)
    real(KREAL), allocatable :: f(:), faux(:)
END MODULE
!
!                 END OF MODULE DAMDEN320_T
!...............................................................................................
!===============================================================================================
!                 MODULE DAMPOT320_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @details Contains variables used by TDAMPOT320_mpi.F90
MODULE DAMPOT320_T
    USE DAMINITIAL_T
    USE DAM320_T
    IMPLICIT NONE
    integer(KINT) :: mxarray
    integer(KINT), parameter :: mxlargo = 50
    real(KREAL) :: geomthr, umbrlargo, rtab(3,mxrtab)
    real(KREAL), allocatable :: xajustd(:,:), rlargo(:)
    integer(KINT8) :: kntlargo, kntcorto, kntcortotot, kntlargotot
    integer(KINT) :: lmaxrep, numrtab
    integer(KINT), allocatable, dimension (:,:) :: icfposd(:,:)
    real(KREAL), allocatable :: Qgacum(:,:), qpacum(:,:), Qllargo(:), ra2l1(:), ra2l1inv(:)
    integer(KINT), allocatable :: lcorto(:,:), llargo(:,:)
END MODULE
!
!                 END OF MODULE DAMPOT320_T
!    MODULE STEPIT_DATA_T
!        USE DAMINITIAL_T, ONLY: KREAL,KINT
!        integer(kint), parameter :: maxvar=3, npnt=5
!        real(kreal) ::xsi(maxvar),xmin(maxvar),xmax(maxvar),deltax(maxvar),delmin(maxvar), &
!                      err(maxvar,maxvar),mukhavata(maxvar)
!        real(kreal) :: chisq,ratio,colin,ack,signif,huge
!        integer(kint) :: nv,ntrace,matrix,nfmax,nflat,mosque,ncomp,jvary
!        integer(kint) :: ncalls
!    END MODULE
!...............................................................................................
!> @author Anmol Kumar
!> @details The module contains variable which are used by subroutine gradpath
!> and basin. 
!> @param fdisp It defines the first displacement (in a.u.) away from a critical point (CP) to start the gradient path.
!> Small variation in its value may lead the gradient path to acquire a
!> completely different direction. The displacement is made in the positive and
!> neative direction of the three eigenvectors. fdisp = 0.04 is fixed in the code
!> @param step The norm of the stepsize used to follow the gradient path.
!> @param drcutnuc Value of gradient to decide if the path is near to a nucleus. drcutnuc=1000.0 
!> @param distcp Value of distance to decide if the gradient path is near to a CP. distcp=0.04 au 
!> @param distnuc Value of distance to decide if the gradient path is near to a nucleus. distnuc = 0.04 au  
!> @param ncp number of CPs.
!> @param xcp x coordinate of CPs. ycp and zcp stores y and z coordinate respectively.
!> @param xboxgs strting x coordinate of the box within which gradient path will be generated. Other variables are
!> XBOXGE,YBOXGS,YBOXGE,ZBOXGS,ZBOXGE
MODULE GRADPATH_T
    USE DAMINITIAL_T
    INTEGER NONE
!        REAL(KREAL), PARAMETER     :: FDISP=0.02, STEP=0.001
    REAL(KREAL), PARAMETER     :: STEP=0.001
    ! Box size in choosen for grad path generation. 6 au away from farthest
    ! coordinates
!        Modificado porRafa: drcutcp moved to module DAMINITIAL_T (to be in namelist OPTIONS)
!        REAL(KREAL), PARAMETER     :: DRCUTcp = 10e-4,DRCUTnuc = 1000.0
    REAL(KREAL), PARAMETER     :: DRCUTnuc = 1000.0
!        Fin de la modificacion de Rafa
    REAL(KREAL), PARAMETER     :: DISTcp = 0.04,DISTnuc = 0.04
    INTEGER(KINT)             :: flag
    INTEGER(KINT)               :: NCP
    INTEGER(KINT)               :: cnt,cnt1, cnt2
    INTEGER(KINT)               :: CPNUM,IXPLN,IYPLN,IZPLN
    INTEGER(KINT)               :: TNLRM, NRM
    INTEGER(KINT),ALLOCATABLE :: isym(:)
    INTEGER(KINT), ALLOCATABLE  :: ncp1(:),ncp2(:),ngp(:),ind(:),nrma(:)
    REAL(KREAL)                :: XBOXGS,XBOXGE,YBOXGS,YBOXGE,ZBOXGS,ZBOXGE
    REAL(KREAL)                :: xp,yp,zp
    REAL(KREAL)                :: drsq, XDIFF, YDIFF, ZDIFF, TEMP1
    REAL(KREAL)                :: xpp,ypp,zpp
    REAL(KREAL)                :: vdummy
    REAL(KREAL), DIMENSION(6)  :: XNEW, YNEW, ZNEW
    REAL(KREAL), DIMENSION(72) :: XDVD, YDVD, ZDVD
    REAL(KREAL), DIMENSION(50,72,200)::xbas,ybas,zbas
    INTEGER(KINT), DIMENSION(50,72)::npbas
    REAL(KREAL), ALLOCATABLE   :: XCP(:), YCP(:), ZCP(:)
    LOGICAL               :: THERE
    CHARACTER             :: CPFILE*20, SYMB*2, CPSYMB*2, NSYMB*2
    CHARACTER             :: dummy*80, symfound*2, symthis*2
    CHARACTER (len=120)   :: DRM
    CHARACTER(2),ALLOCATABLE ::symcp2(:)!,atsy(:)
    CHARACTER*1, ALLOCATABLE ::symcp(:),symcp1(:)!,cpsy(:)
!        INTEGER(KINT), ALLOCATABLE :: atin(:),cpin(:)
    INTEGER(KINT):: NUMENT
    REAL(KREAL), ALLOCATABLE    :: Eval(:,:)
    REAL(KREAL), ALLOCATABLE    :: Evec(:,:,:)
END MODULE

!...............................................................................................
!===============================================================================================
!                 MODULE PARALELO_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @details Used for defining variables used in MPI. 
MODULE PARALELO_T
    USE DAM320_T
    IMPLICIT NONE
    integer(KINT) :: nprocs, myrank, istart, iend
    integer(KINT), allocatable :: nbasesac(:), istav(:), iendv(:), ilenv(:), idispv(:)
    integer(KINT) :: abort, abortroot
    character(256) :: fname, fnamerank
    logical lwrtcab
END MODULE
!
!                 END OF MODULE PARALELO_T
!...............................................................................................
!===============================================================================================
!                 MODULE GAUSS_T
!===============================================================================================
!> @author Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema, Guillermo Ramirez
!> @note Polynomials P_k^(L,M;L',M') of the shift-operators technique
!> Pointers to the elements P_0^(L,-L; L',-L')
!> Hypergeometrics 2F1 and overlap integrals between generalized gaussians < g(n,L,M)| g(0,L',M') >
!> Multipolar moments of two-center Gaussian distributions with shift-operators 
!> reciprocals of semiinteger factorials
MODULE GAUSS_T
    USE DAM320_T
    USE DAM320_CONST_T
    IMPLICIT NONE
    logical :: lbeta
    integer(KINT), parameter :: mxprimit = 20
    integer(KINT) :: nprimitot, ncontrtot, nocalfa, nocbeta
    integer(KINT), allocatable :: isort(:), ncontr(:), nprimit(:), ipntprim(:)
    real(KREAL), allocatable :: cfcontr(:), cfcontr0(:), aorb(:,:), borb(:,:), xxg(:), xxg0(:)
!		Polynomials P_k^(L,M;L',M') of the shift-operators technique 
    real(KREAL), allocatable :: polP(:)
!		Pointers to the elements P_0^(L,-L; L',-L')
    integer(KINT), allocatable  :: ipntpolP(:,:)
!		Hypergeometrics 2F1 and overlap integrals between generalized gaussians < g(n,L,M)| g(0,L',M') >
    real(KREAL), allocatable :: h2f1(:,:,:), sint(:,:,:)
!		Multipolar moments of two-center Gaussian distributions with shift-operators 
    real(KREAL), allocatable :: qlm2calt(:)
!		reciprocals of semiinteger factorials
    real(KREAL) :: factsi(-1:mxfact), occalfa, occbeta
END MODULE
!
!                 END OF MODULE GAUSS_T
