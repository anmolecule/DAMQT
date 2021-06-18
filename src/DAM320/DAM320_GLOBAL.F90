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
!	Modules for programs DAMSTO320, DAMGTO320, DAMDEN320, DAMPOT320, DAMFIELD320, DAMFIELDPNT320
!		DAMFRAD320, DAMSTO320_mpi, DAMGTO320_mpi, DAMDEN320_mpi, DAMPOT320_mpi, DAMMULTROT320,
!		DAMORB320, DAMORB320_mpi
!
! Version of April 2019
!	
!===============================================================================================
!                 MODULE DAM320_D
!===============================================================================================
MODULE DAM320_D
     IMPLICIT NONE
     integer(4), parameter :: KINT = 4, KREAL = 8, KREAL4 = 4, KINT8 = 8
!    The following fundamental parameters set the dimensions:
!    mxcap: highest number of shells allowable in basis set
!    mxcen: highest number of centers
!    mxl:   highest l value of STO or GTO basis set
!    mxn:   highest n value of STO
!    mxlexp: highest l of expansions in spherical harmonics
!    mxmult: highest order of molecular multipoles
!    mxsigexp:	  highest expansion order for the B integrals in frgsigma
!    mxrtab:	  highest number of tabulation points for density and potential
     integer(KINT), parameter ::  mxcen = 10000, mxl = 6, mxn = 9, mxcap = 30000
     integer(KINT), parameter ::  mxmult = 4, mxrtab = 200
     integer(KINT), parameter ::  mxsigexp = 100
!    The following axiliary parameters derive from the previous ones
!    mxldst: highest l of distributions
     integer(KINT), parameter ::  mxldst = mxl+mxl
     logical iswindows	! .true. if running on a MS-windows system
     character dirsep	! character for directory names separator: "/" for unix, "\\" for MS-windows
     integer(KINT), parameter :: nintervaj = 29	! Number of intervals for piecewise representation of the radial factors
!    Radii of fitting intervals
     real(KREAL), parameter :: rinterv(0:nintervaj) = (/0.d0, 1.d-2, 3.d-2, 5.d-2, .7d-1, 1.d-1, 2.d-1, 3.d-1, 4.d-1, 5.d-1, &
             6.d-1, 7.d-1, 8.d-1, 9.d-1, 1.d0, 1.5d0, 2.d0, 2.5d0, 3.d0, 4.d0, &
             5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 12.d0, 14.d0, 16.d0, 20.d0/)
     integer(KINT), parameter :: mxlenpol = 20	! Highest number of fitting polynomials
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
!                 END OF MODULE DAM320_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAM320_DATA_D
!===============================================================================================
MODULE DAM320_DATA_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT), parameter :: npntintr = (mxlenpol*3+1) / 2
     integer(KINT), parameter :: npntaj = npntintr * nintervaj
!    umbral:	threshold for neglecting radial factors
!    umbralres:	threshold for truncating radial factors expansions
!    xx:		basis set STO exponents
!    rcen:	Cartesian coordinates of atomic centers
!    rnor:	basis set STO radial normalization factors
!    zn:		atomic number
!    lmaxexp:	length of expansion on l
!    ncen:	number of atomic centers
!    nbas:	number of basis functions
!    ncaps:	number of shells
!    ioptaj:	fitting option
!    nn:		n quantum numbers of STO basis set
!    ll:		l quantum numbers of STO basis set
!    nf:		absolute order index of first function in a shell
!    ngini:	absolute order index of the first shell of a center
!    ngfin:	absolute order index of the last shell of a center
!    lastfnc:	index of last function of a given center
!    lmaxc:	highest l quantum number in a center
!    lmaxbase:	highest l quantum number in the basis set
!    lmultmx	highest l of multipoles whose modules are stored
!    atmnam:	atomic symbols of centers
!    projectname:	project name for input and output files
!    lsdisf: if .false., negligible charge distribution
!    longoutput: it .true. more detailed output
!    lsto:	if .true. Slater basis set, otherwise, Gaussian basis set
!    dmat:	full density matrix
     real(KREAL) :: umbral, umbralres, thresoverlap
     integer(KINT) :: ioptaj, lmaxbase, lmaxexp, lmtop, lmultmx, mxltot, mxtpw, nbas, ncaps, ncen, nfitpar, nintstd, numdvec
     character(2), allocatable :: atmnam(:)
     character(300) :: projectname
     logical :: lgbsgz, lden, ldengz, ldensprsbin, lm2c, longoutput, lsto, lvalence, lzdo
     real(KREAL) :: roblk(-mxl:mxl,-mxl:mxl)
     logical*1, allocatable :: lsdisf(:,:)
     real(KREAL), allocatable :: dmat(:,:), dmataux(:,:), dvec(:), rcen(:,:), rnor(:), xx(:), zn(:)
     integer(KINT), allocatable :: iblock(:), ivec(:), jblock(:), jvec(:), kntblock(:,:), nblocks(:), ndstblocks(:)
     integer(KINT), allocatable :: ll(:), lmaxc(:), nexpaj(:), ngini(:), ngfin(:), nn(:), nf(:), nzn(:)
!    xajust: 	exponential (nonlinear) fitting parameters
!    cfajust:	linear fitting parameters
!    cfrint1:  auxiliary integrals for potential and field
!    cfrint2l2:  auxiliary integrals for potential and field
!    nintervaj: number of fitting intervals of fragments
!    nexpaj:	number of intervals with exponential term
!    lenpol:   number of terms in the polynomial
!    ipmax:	number of terms in the expansion of the B integrals in frgsigma
!    ncheb:	array with different numbers of tabulation points per center
     real(KREAL), allocatable :: dl(:,:,:), expvec(:), expvinv(:), f0(:,:), f0ab(:,:,:), fa(:,:), ftab(:,:), ha(:,:)
     real(KREAL), allocatable :: r2pow(:,:), radio(:,:), residuals(:), rl(:,:,:), rlt(:), rpntaj(:), rpow(:,:), xajust(:)
     real(KREAL), allocatable :: Qgpart(:), qppart(:), rintr1(:), rintr2l2(:), vaux1c(:,:), ymat1(:,:), ymat2(:,:)
     real(KREAL), allocatable :: dosxcheb(:), fkvec0(:), rmultipmod(:), tchvec0(:), tchvec0a(:), tchvec1(:), xcheb(:)
     real(KREAL), allocatable, target, dimension(:) :: cfajust, cfrint1, cfrint2l2, cfajustrank, cfrint1rank, cfrint2l2rank
     real(KREAL), allocatable, target, dimension(:,:) :: rmultip, rmultipfr
     integer(KINT), allocatable, dimension(:) :: indintrv
     integer(KINT), allocatable, dimension(:) :: icfpos, icfposrank
     integer(KINT) :: ipmax, lencfparank, lenexpaj, natomtype, lenxajust
!    fct:		auxiliary factor for locating the interval where a point is placed
!    wthreshold: threshold for contributions for the B integrals (see subroutine Avk)
     real(KREAL) :: fct, wthreshold
!    indices and arrays for multipoles calculation
     integer(KINT) :: ilow, iupp
!    auxiliary arrays for multipolar moments of STO distributions
     real(KREAL), allocatable :: bkmat(:,:), qlm2c(:), qlmdst(:), qlmasint(:), powu(:,:), pow1mu(:,:)
     real(KREAL) :: auxu(0:mxn+mxl), cina(0:mxn), cinb(0:mxn), pow2(0:mxldst), scomp(0:mxn+2*mxl)
!    auxiliary arrays for multipolar moments of CGTO distributions
     real(KREAL), allocatable :: besselint(:,:)
     integer(KINT) :: kntintB, kmltan, kmltquad
     real(KREAL), allocatable :: av(:), bv(:), sol(:)	! for overlap integrals with ellipsoidal coords
END MODULE
!
!                 END OF MODULE DAM320_DATA_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAM320_CONST_D
!=============================================================================================== 
MODULE DAM320_CONST_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT), parameter :: mxroot = 800, mxreal = 3000, mxfact = 150
     real(KREAL), parameter :: cero = 0.d0, uno = 1.d0, dos = 2.d0, cuatro = 4.d0, ocho = 8.d0, r16 = 16.d0, udec = 0.1d0
     real(KREAL), parameter :: umed = .5d0, pt25 = .25d0, euler = 5.7721566490153280D-1, raiz2 = 1.4142135623730950d0
     real(KREAL), parameter :: umbrzn = 1.d-6
     integer(KINT) :: idimlml1l2, mxbin, mxemes, mxind, mxkcof, mxlcof
!    re(i) = double(i)
!    ri(i) = double(1/i)
!    fact(i) = double(i!)
!    facti(i) = double(1/i!)
!    dosl1(i) = double(2l+1))
!    dosl1i(i) = 1 / double(2l+1))
     real(KREAL) :: re(-mxreal:mxreal), ri(-mxreal:mxreal), fact(0:mxfact), facti(0:mxfact) &
             , facts(-1:mxfact), dosl1(-mxreal:mxreal), dosl1i(-mxreal:mxreal)
     real(KREAL) :: pi, raizpi, pimed, pimedsqr   ! Sqrt[pi/2]
!    	ind(i) = i*(i+1) / 2
     integer(KINT), allocatable :: ind(:)
!    	root(i) = double( sqrt(i) )
!    	rooti(i) = double( 1 / sqrt(i) )
     real(KREAL) :: root(0:mxroot), rooti(mxroot)
!    ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
     real(KREAL), allocatable :: app(:,:), ang(:), bin(:), bpp(:,:), ccl1l2(:), ssv(:,:), sdv(:,:), umedpow(:)
!    	Coefficients and m values for the decomposition of the product of two Phi_m(phi) functions (of the real spherical harmonics)
     integer(KINT), allocatable :: msv(:,:), mdv(:,:), indk12(:,:)
!    	Coefficients of Chebyshev T polynomials
     real(KREAL) ::  chebTcf((mxlenpol/2+1)*((mxlenpol+1)/2+1))
!    	Coefficients for the decomposition of products of pairs of real spherical harmonics
     integer(KINT), allocatable :: i1l1l2(:,:), i2l1l2(:,:), lml1l2(:), npl1l2(:), llm(:), mlm(:), lll1l2(:), mml1l2(:)
!    	Atomic symbols
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
!    Hypergeometrics 3F2(-j,j+2,1/2;3/2,3/2;1)
     real(KREAL),parameter :: h3f2(0:30) =	(/ 1.0000000000000000D0, 3.3333333333333330D-1, 2.8888888888888880D-1 &
     , 1.8095238095238100D-1, 1.6698412698412700D-1, 1.2400192400192400D-1, 1.1727637441923160D-1, 9.4283494283494300D-2 &
     , 9.0343498186635500D-2, 7.6045990473235050D-2, 7.3461722941036200D-2, 6.3716724290152750D-2 &
     , 6.1892360883217920D-2, 5.4825975317485430D-2, 5.3469760871032360D-2, 4.8111771784334780D-2 &
     , 4.7064198756022750D-2, 4.2862219460053240D-2, 4.2028788279367650D-2, 3.8645297583347980D-2 &
     , 3.7966485503304700D-2, 3.5183653709814110D-2, 3.4620113210353590D-2, 3.2291083670560490D-2 &
     , 3.1815766854350310D-2, 2.9837935700827030D-2, 2.9431638990733510D-2, 2.7731158377428100D-2 &
     , 2.7379872214383080D-2, 2.5902238055824550D-2, 2.5595502719280240D-2 /)
!    for overlap integrals with ellipsoidal coords
     integer(KINT), allocatable :: ipntap(:,:), ipntalfa(:,:)
     real(KREAL), allocatable :: ap(:), alfasol(:), alm(:)
     real(KREAL), parameter :: alfacutoff = 1.d-20
END MODULE
!
!                 END OF MODULE DAM320_CONST_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAM320_LEGENDRE_D
!===============================================================================================
MODULE DAM320_LEGENDRE_D
     USE DAM320_D
!    Abscisae and weights for Legendre quadrature rules with 25 and 35 points
     integer(KINT) :: nqleg, nqleg2
     integer(KINT), parameter :: mxleg = 35
     real(KREAL), pointer :: u(:), w(:), u2(:), w2(:)
!    Data generated with notebook legendre_quadrature.nb
     real(KREAL), target :: uleg25(25) = (/ 2.2215151047509510D-3,1.1668039270241240D-2,2.8512714385512830D-2   &
             ,5.2504001060862320D-2,8.3278685619583000D-2,1.2037036848132120D-1,1.6321681576326580D-1          &
             ,2.1116853487938850D-1,2.6349863427714250D-1,3.1941384709530610D-1,3.7806655813950580D-1          &
             ,4.3856765369464480D-1,5.0000000000000000D-1,5.6143234630535520D-1,6.2193344186049420D-1          &
             ,6.8058615290469390D-1,7.3650136572285750D-1,7.8883146512061100D-1,8.3678318423673400D-1          &
             ,8.7962963151867900D-1,9.1672131438041700D-1,9.4749599893913800D-1,9.7148728561448700D-1          &
             ,9.8833196072975900D-1,9.9777848489524900D-1/)
     real(KREAL), target :: wleg25(25) = (/ 5.6968992505131440D-3,1.3177493307516069D-2,2.0469578350653156D-2   &
             ,2.7452347987917596D-2,3.4019166906178459D-2,4.0070350167500509D-2,4.5514130991481825D-2          &
             ,5.0267974533525322D-2,5.4259812237131827D-2,5.7429129572855824D-2,5.9727881767892386D-2          &
             ,6.1121221495155021D-2,6.1588026863357726D-2,6.1121221495155021D-2,5.9727881767892386D-2          &
             ,5.7429129572855824D-2,5.4259812237131827D-2,5.0267974533525322D-2,4.5514130991481825D-2          &
             ,4.0070350167500509D-2,3.4019166906178459D-2,2.7452347987917596D-2,2.0469578350653156D-2          &
             ,1.3177493307516069D-2,5.6968992505131440D-3/)
     real(KREAL), target :: uleg35(35) = (/ 1.1467154501998510D-3,6.0321177780742500D-3,1.4781191980385080D-2   &
             ,2.7327425896086340D-2,4.3572869320341200D-2,6.3390437487388800D-2,8.6625050453887300D-2          &
             ,1.1309487385654370D-1,1.4259274922168560D-1,1.7488781766705480D-1,2.0972732762511770D-1          &
             ,2.4683861337925570D-1,2.8593122924109290D-1,3.2669922278459300D-1,3.6882352939535200D-1          &
             ,4.1197446941700500D-1,4.5581432836217000D-1,5.0000000000000000D-1,5.4418567163783000D-1          &
             ,5.8802553058299500D-1,6.3117647060464800D-1,6.7330077721540700D-1,7.1406877075890700D-1          &
             ,7.5316138662074400D-1,7.9027267237488200D-1,8.2511218233294500D-1,8.5740725077831400D-1          &
             ,8.8690512614345600D-1,9.1337494954611300D-1,9.3660956251261100D-1,9.5642713067965900D-1          &
             ,9.7267257410391400D-1,9.8521880801961500D-1,9.9396788222192600D-1,9.9885328454980000D-1/)
     real(KREAL), target :: wleg35(35) = (/ 2.9417167102215425D-3,6.8254141741807461D-3,1.0661489955741790D-2   &
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
!                 END OF MODULE DAM320_LEGENDRE_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAMDEN320_D
!===============================================================================================
MODULE DAMDEN320_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT), parameter :: mxsel = 200
     logical :: langstrom, lexact, lmolec, latomics, latomsel, ldensacc, lboundsx, lboundsy, lboundsz
     logical :: lgradient, laplacian, ldeform, lderiv2, lgrid, lgrid2d, lpoints
     logical, allocatable :: lselat(:), lnegia(:)
     real(KREAL) :: uinf, usup, dltu, vinf, vsup, dltv
     real(KREAL) :: xinf, xsup, dltx, yinf, ysup, dlty, zinf, zsup, dltz, planeA, planeB, planeC
     real(KREAL) :: umbrlargo, rtab(3,mxrtab)
     real(KREAL) :: xboundinf, xboundsup, yboundinf, yboundsup, zboundinf, zboundsup
     real(KREAL), allocatable :: xajustd(:,:), rlargo(:)
     integer(KINT) :: iatomsel(mxsel), idimzlm, lminrep, lmaxrep, nsel, numrtab, planecase
     integer(KINT), allocatable :: icfposd(:,:), lcorto(:,:)
     real(KREAL), allocatable :: f(:), faux(:)
     real(KREAL), allocatable :: zlma(:), zlmadx(:), zlmady(:), zlmadz(:)
     real(KREAL), allocatable :: zlmadxx(:), zlmadxy(:), zlmadxz(:), zlmadyy(:), zlmadyz(:), zlmadzz(:)
     character(4) :: planesuffix
     character(230) :: filename	! Name for .plt files
     character(256) :: x_func_uv, y_func_uv, z_func_uv  ! Expressions of (x,y,z) in terms of (u,v) for 2D grids
END MODULE
!
!                 END OF MODULE DAMDEN320_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAMPOT320_D
!===============================================================================================
MODULE DAMPOT320_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT) :: mxarray
     integer(KINT), parameter :: mxlargo = 50
     logical :: langstrom, largo, lexact, lgrid, lgrid2d, lpoints, lgradient, lderiv2
     real(KREAL) :: uinf, usup, dltu, vinf, vsup, dltv
     real(KREAL) :: xinf, xsup, dltx, yinf, ysup, dlty, zinf, zsup, dltz, planeA, planeB, planeC
     real(KREAL) :: geomthr, umbrlargo, umbrlargo2, rtab(3,mxrtab)
     real(KREAL), allocatable :: xajustd(:,:), rlargo(:)
     integer(KINT8) :: kntlargo, kntcorto, kntcortotot, kntlargotot
     integer(KINT) :: idimzlm, lmaxrep, numrtab, planecase
     integer(KINT), allocatable :: icfposd(:,:)
     real(KREAL), allocatable :: Qgacum(:,:), qpacum(:,:), Qllargo(:), ra2l1(:), ra2l1inv(:)
     real(KREAL), allocatable :: zlma(:), zlmadx(:), zlmady(:), zlmadz(:)
     real(KREAL), allocatable :: zlmadxx(:), zlmadxy(:), zlmadxz(:), zlmadyy(:), zlmadyz(:), zlmadzz(:)
     integer(KINT), allocatable :: lcorto(:,:), llargo(:,:)
     character(4) :: planesuffix
     character(230) :: filename	! Name for .plt files
     character(256) :: x_func_uv, y_func_uv, z_func_uv  ! Expressions of (x,y,z) in terms of (u,v) for 2D grids
END MODULE
!
!                 END OF MODULE DAMPOT320_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAMFIELD320_D
!===============================================================================================
MODULE DAMFIELD320_D
     USE DAM320_D
     IMPLICIT NONE
     logical :: lextralines, largo, lcompl, larrows, lplot2d, lstepindiv
     integer(KINT), parameter :: mxlargo = 50
     integer(KINT) :: idimzlm, lmaxrep, kntlargo, kntcorto, nlineas, nlincomp, nlines, narrstep, nsize, planecase
     integer(KINT) :: icntlines(mxrtab)
     integer(KINT), allocatable :: icfposd(:,:)
     real(KREAL) :: uinf, usup, dltu, vinf, vsup, dltv
     real(KREAL) :: xinf, xsup, yinf, ysup, zinf, zsup, dlt0, thresh, rlongarr, rwidearr
     real(KREAL) :: umbrlargo, uvratio, planeA, planeB, planeC
     real(KREAL) :: rlines(3,mxrtab), wu(3), wv(3)
     real(KREAL), allocatable :: xajustd(:,:), rlargo(:)
     real(KREAL), allocatable :: Qgacum(:,:), qpacum(:,:), Qllargo(:), ra2l1(:), ra2l1inv(:)
     real(KREAL), allocatable :: zlma(:), zlmadx(:), zlmady(:), zlmadz(:)
     real(KREAL), allocatable :: zlmadxx(:), zlmadxy(:), zlmadxz(:), zlmadyy(:), zlmadyz(:), zlmadzz(:)
     integer(KINT), allocatable :: lcorto(:,:), llargo(:,:)
     character(256) :: filename	! Name for .cam files
END MODULE
!
!                 END OF MODULE DAMFIELD320_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAMDENGRAD320_D
!===============================================================================================
MODULE DAMDENGRAD320_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT), parameter :: mxlargo = 50
     logical :: lextralines,  lcompl, lcpsd, lplot2d
     real(KREAL) :: uinf, usup, dltu, vinf, vsup, dltv
     real(KREAL) :: xinf, xsup, yinf, ysup, zinf, zsup, dlt0, thresh, rlongarr, rwidearr
     real(KREAL) :: umbrlargo, uvratio, planeA, planeB, planeC
     real(KREAL) :: rlines(3,mxrtab), wu(3), wv(3)
     real(KREAL), allocatable :: rdenmax(:,:), ra2l1(:), ra2l1inv(:), rlargo(:), xajustd(:,:), zlma(:), zlmadx(:), zlmady(:)
     real(KREAL), allocatable :: zlmadz(:), zlmadxx(:), zlmadxy(:), zlmadxz(:), zlmadyy(:), zlmadyz(:), zlmadzz(:)
     integer(KINT) :: idimzlm, lmaxrep, kntlargo, kntcorto, nlineas, nlincomp, nlines, narrstep, planecase
     integer(KINT) :: icntlines(mxrtab)
     integer(KINT), allocatable :: icfposd(:,:)
     integer(KINT), allocatable :: lcorto(:,:)
     character(256) :: filename	! Name for .dengr files
END MODULE
!
!                 END OF MODULE DAMDENGRAD320_D
!...............................................................................................
!===============================================================================================
!                 MODULE ICOSAHEDRON
!===============================================================================================
MODULE ICOSAHEDRON
  USE DAM320_D
  IMPLICIT NONE
  integer(KINT), parameter :: idimvert = 12, idimc2 = 30, idimc3 = 20
  real(KREAL) :: vertices(idimvert,3) = transpose(reshape((/ -0.525731112119134d0, -0.85065080835204d0, 0.d0, &
      -0.525731112119134d0, 0.85065080835204d0, 0.d0, 0.d0, -0.525731112119134d0, -0.85065080835204d0, 0.d0, &
      -0.525731112119134d0, 0.85065080835204d0, 0.d0, 0.525731112119134d0, &
      -0.85065080835204d0, 0.d0, 0.525731112119134d0, 0.85065080835204d0, &
       0.525731112119134d0, -0.85065080835204d0, 0.d0, 0.525731112119134d0, &
       0.85065080835204d0, 0.d0, -0.85065080835204d0, 0.d0, -0.525731112119134d0, &
      -0.85065080835204d0, 0.d0, 0.525731112119134d0, 0.85065080835204d0, 0.d0, &
      -0.525731112119134d0, 0.85065080835204d0, 0.d0, 0.525731112119134d0 /), (/3,idimvert/)))
  real(KREAL) :: c2axes(idimc2,3) = transpose(reshape((/ -0.309016994374947d0, -0.809016994374947d0, -0.5d0, &
      -0.309016994374947d0, -0.809016994374947d0, 0.5d0, -0.309016994374947d0, 0.809016994374947d0, -0.5d0, &
      -0.309016994374947d0, 0.809016994374947d0, 0.5d0, 0.d0, 0.d0, -1.d0, 0.d0, 0.d0, 1.d0, &
       0.d0, -1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.309016994374947d0, -0.809016994374947d0, -0.5d0, &
       0.309016994374947d0, -0.809016994374947d0, 0.5d0, 0.309016994374947d0, &
       0.809016994374947d0, -0.5d0, 0.309016994374947d0, 0.809016994374947d0, 0.5d0, &
      -0.809016994374947d0, -0.5d0, -0.309016994374947d0, -0.809016994374947d0, &
      -0.5d0, 0.309016994374947d0, -0.809016994374947d0, 0.5d0, -0.309016994374947d0, &
      -0.809016994374947d0, 0.5d0, 0.309016994374947d0, -1.d0, 0.d0, 0.d0, -0.5d0, &
      -0.309016994374947d0, -0.809016994374947d0, -0.5d0, -0.309016994374947d0, &
       0.809016994374947d0, -0.5d0, 0.309016994374947d0, -0.809016994374947d0, -0.5d0, &
       0.309016994374947d0, 0.809016994374947d0, 0.5d0, -0.309016994374947d0, &
      -0.809016994374947d0, 0.5d0, -0.309016994374947d0, 0.809016994374947d0, 0.5d0, &
       0.309016994374947d0, -0.809016994374947d0, 0.5d0, 0.309016994374947d0, &
       0.809016994374947d0, 1.d0, 0.d0, 0.d0, 0.809016994374947d0, -0.5d0, &
      -0.309016994374947d0, 0.809016994374947d0, -0.5d0, 0.309016994374947d0, &
       0.809016994374947d0, 0.5d0, -0.309016994374947d0, 0.809016994374947d0, 0.5d0, &
       0.309016994374947d0 /), (/3,idimc2/)))
  real(KREAL) :: c3axes(idimc3,3) = transpose(reshape((/ 0.d0, -0.934172358962716d0, -0.35682208977309d0, 0.d0, &
      -0.934172358962716d0, 0.35682208977309d0, 0.d0, 0.934172358962716d0, -0.35682208977309d0, 0.d0, &
       0.934172358962716d0, 0.35682208977309d0, -0.934172358962716d0, &
      -0.35682208977309d0, 0.d0, -0.934172358962716d0, 0.35682208977309d0, 0.d0, &
      -0.577350269189626d0, -0.577350269189626d0, -0.577350269189626d0, &
      -0.577350269189626d0, -0.577350269189626d0, 0.577350269189626d0, &
      -0.577350269189626d0, 0.577350269189626d0, -0.577350269189626d0, &
      -0.577350269189626d0, 0.577350269189626d0, 0.577350269189626d0, &
      -0.35682208977309d0, 0.d0, -0.934172358962716d0, -0.35682208977309d0, 0.d0, &
       0.934172358962716d0, 0.35682208977309d0, 0.d0, -0.934172358962716d0, &
       0.35682208977309d0, 0.d0, 0.934172358962716d0, 0.577350269189626d0, &
      -0.577350269189626d0, -0.577350269189626d0, 0.577350269189626d0, &
      -0.577350269189626d0, 0.577350269189626d0, 0.577350269189626d0, &
       0.577350269189626d0, -0.577350269189626d0, 0.577350269189626d0, &
       0.577350269189626d0, 0.577350269189626d0, 0.934172358962716d0, &
      -0.35682208977309d0, 0.d0, 0.934172358962716d0, 0.35682208977309d0, 0.d0 /), (/3,idimc3/)))
END MODULE
!
!                 END OF MODULE ICOSAHEDRON
!...............................................................................................
!===============================================================================================
!                 MODULE DAMORB320_D
!===============================================================================================
MODULE DAMORB320_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT), parameter :: mxsel = 100
     integer(KINT) :: iorbsinp(mxsel)
     character(256) :: filename, fileMOname
     logical :: langstrom, lgradient, lgrid2d
     real(KREAL) :: uinf, usup, dltu, vinf, vsup, dltv
     real(KREAL) :: xinf, xsup, dltx, yinf, ysup, dlty, zinf, zsup, dltz, planeA, planeB, planeC
     real(KREAL), allocatable :: chi(:), chidx(:), chidy(:), chidz(:), corbs(:,:), orbv(:), orbvdx(:), orbvdy(:), orbvdz(:)
     real(KREAL), allocatable :: zlma(:), zlmadx(:), zlmady(:), zlmadz(:)
     integer(KINT) :: planecase
     integer(KINT), allocatable :: iorbs(:) ! For orbitals plotting
     character(256) :: x_func_uv, y_func_uv, z_func_uv  ! Expressions of (x,y,z) in terms of (u,v) for 2D grids
     character(4) :: planesuffix
END MODULE
!
!                 END OF MODULE DAMORB320_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAMFORCES320_D
!===============================================================================================
MODULE DAMFORCES320_D
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT), parameter :: mxsel = 200
     logical :: latomics, latomsel
     integer(KINT) :: iatomsel(mxsel), nsel
     real(KREAL), allocatable :: qptotal(:,:)
END MODULE
!
!                 END OF MODULE DAMFORCES320_D
!...............................................................................................
!===============================================================================================
!                 MODULE PARALELO
!===============================================================================================
MODULE PARALELO
     USE DAM320_D
     IMPLICIT NONE
     integer(KINT) :: nprocs, myrank, istart, iend
     integer(KINT), allocatable :: nbasesac(:), istav(:), iendv(:), ilenv(:), idispv(:)
     integer(KINT) :: abort, abortroot
     character(256) :: fname, fnamerank
     logical lwrtcab
END MODULE
!
!                 END OF MODULE PARALELO
!...............................................................................................
!===============================================================================================
!                 MODULE GAUSS
!===============================================================================================
MODULE GAUSS
     USE DAM320_D
     USE DAM320_CONST_D
     IMPLICIT NONE
     logical :: lbeta
     integer(KINT), parameter :: mxprimit = 20
     integer(KINT) :: nprimitot, ncontrtot, nocalfa, nocbeta
     integer(KINT), allocatable :: isort(:), ncontr(:), nprimit(:), ipntprim(:)
     real(KREAL), allocatable :: cfcontr(:), cfcontr0(:), xxg(:), xxg0(:)
!    	Polynomials P_k^(L,M;L',M') of the shift-operators technique
     real(KREAL), allocatable :: polP(:)
!    	Pointers to the elements P_0^(L,-L; L',-L')
     integer(KINT), allocatable  :: ipntpolP(:,:)
!    	Hypergeometrics 2F1 and overlap integrals between generalized gaussians < g(n,L,M)| g(0,L',M') >
     real(KREAL), allocatable :: h2f1(:,:,:), sint(:,:,:)
!    	Multipolar moments of two-center Gaussian distributions with shift-operators
     real(KREAL), allocatable :: qlm2calt(:)
!    	reciprocals of semiinteger factorials
     real(KREAL) :: factsi(-1:mxfact), occalfa, occbeta
END MODULE
!
!                 END OF MODULE GAUSS
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
if (ipos > 0) cop=str(ipos:ipos)
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
