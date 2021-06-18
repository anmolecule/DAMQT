!  Copyright 2011-2019, Rafael Lopez, Gabriel Aires Urquiza de Carvalho
! 
!  This file is part of Zernike_Jacobi_320 package.
!  Zernike_Jacobi_320 is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
! 
!  You should have received a copy of the GNU General Public License
!  along with Zernike_Jacobi_320.  If not, see <http://www.gnu.org/licenses/>.
!
!------------------------------------------------------------------------
! 
! Program for computation of Zernike 3D and Jacobi moments of a molecular density 
!       resultant from a sum of atomic densities each one expressed as a linear combination of GTO radial factors
!       times unnormalized real spherical harmonics. 
!
! Version of May 2018
!
MODULE DAMZernike_Jacobi_atdens_GTO
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
    logical :: lechelon, ljacobi, lmthfile, lrstarrel
    logical iswindows	! .true. if running on a MS-windows system
    logical :: longoutput, lgbsgz
    integer(KINT), parameter ::  mxcen = 500, mxl = 6, mxcap = 8000
    integer(KINT), parameter ::  mxmult = 4
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
    real(KREAL), allocatable :: radfunction(:,:), vaux(:), weights(:), weights01(:)
    real(KREAL), allocatable :: qlm(:), qlmnuc(:)
    real(KREAL), allocatable, target :: rmultip(:,:)
    real(KREAL) :: rstar, thresmult, thresoverlap
    real(KREAL) :: pi, raizpi, pimed, pimedsqr   ! Sqrt[pi/2]
    real(KREAL) :: roblk(-mxl:mxl)
    real(KREAL) :: fact(0:mxfact), facti(0:mxfact), re(-mxreal:mxreal), ri(-mxreal:mxreal)
    real(KREAL) :: root(0:mxroot), rooti(mxroot)
    real(KREAL) :: dosl1(-mxreal:mxreal), dosl1i(-mxreal:mxreal), facts(-1:mxfact), factsi(-1:mxfact)
!       Coefficients and m values for the decomposition of the product of two Phi_m(phi) functions (of the real spherical harmonics)
    integer(KINT), allocatable :: msv(:,:), mdv(:,:), indk12(:,:)
    real(KREAL), allocatable :: app(:,:), bpp(:,:), ssv(:,:), sdv(:,:)
!       Polynomials P_k^(L,M;L',M') of the shift-operators technique 
    real(KREAL), allocatable :: polP(:)
!       Pointers to the elements P_0^(L,-L; L',-L')
    integer(KINT), allocatable  :: ipntpolP(:,:)
END MODULE
!
!                 END OF MODULE DAMZernike_Jacobi_atdens_GTO
!...............................................................................................
  program DAMZernike_Jacobi_atdens_GTO_320
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, ia, ierr, j, k, l, la, m, ma
    real(KREAL) :: aux
    real(KREAL4) :: tarray(2), tiempo, dtime
    external :: omeganlm_jacobi, omeganlm_zernike, frad_jacobi, frad_zernike
    namelist / options / iswindows, lechelon, ljacobi, lmthfile, longoutput, lexpansion, lrstarrel, &
            nquadpoints, kexpansion, rstar, thresmult, thresoverlap
    tiempo = dtime(tarray)
!    Defaults for the NAMELIST OPTIONS
    kexpansion = 10           ! Highest k in expansion functions
    lechelon = .false.        ! if true number of functions per l equal to max(lexpansion+1,kexpansion)-l
    lexpansion = 20           ! Highest l in expansion functions
    ljacobi = .false.         ! if .true. uses Jacobi P(0,2+2l) polynomials as radial functions instead of Zernike 3D
    lmthfile = .false.        ! if .true. generates a projectname.mth file with data about the moments computation
    longoutput = .false.      ! If .true. a more detailed output is given
    lrstarrel = .false.       ! If .true. ball radius equal to distance of farthest atom plus rstar
    nquadpoints = 128         ! Number of quadrature points
    rstar = 10.d0             ! Ball radius for expansion
    thresmult = 1.d-10        ! Threshold for printing multipole moments
    iswindows = .false.       ! .true. if running on a MS-windows system
!    End of Defaults for the NAMELIST OPTIONS

    read(5,OPTIONS)    !    Reads the namelist OPTIONS
    read(5,*) projectname
    write(6,"(1x,'project name : ',a,/,1x,'==============')") projectname
    write(6,"('lexpansion = ', i3, ' kexpansion = ', i3,/)") lexpansion, kexpansion
!     nquadpoints = min(120,nquadpoints)
    write(6,"(/'****  rstar = ', e17.10, ' ****     Number of quadrature points = ', i5,/)") rstar, nquadpoints

    if (iswindows) then
        dirsep = "\\"
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        if (i .eq. 0) then    ! This is intended for MinGW, whose directory separator in windows is also /
                dirsep = "/"
                i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        endif
    else
        dirsep = "/"
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
    end if
    mxltot = mxldst + lexpansion

    call consta    !    Computes and stores several auxiliary constants and functions

!    Reads geometry and atomic densities expanded in GTO 
    lgbsgz = .false.
    inquire(file=trim(projectname)//".atomic_ggden.gz", exist=lgbsgz, iostat=ierr)
    if (ierr .eq. 0 .and. lgbsgz) then
        call system ("gunzip "//trim(projectname)//".atomic_ggden.gz")
    endif

    call leegeomatdensGTO

    if (lgbsgz) then
        call system ("gzip "//trim(projectname)//".atomic_ggden")
    endif

    write(6,"(/'****  rstar = ', e17.10, ' ****     Number of quadrature points = ', i5,/)") rstar, nquadpoints


!    Computation of Zernike 3D or Jacobi moments

    if (lechelon) kexpansion = max(lexpansion, kexpansion)
    if (ljacobi) then
        call expand(omeganlm_jacobi, frad_jacobi)
    else
        call expand(omeganlm_zernike, frad_zernike)
    endif

!    Computes the total molecular multipolar moments from the Zernike expansion

    tiempo = dtime(tarray)
    write(6,"(1x,//80('='),/'Timing in seconds for fingerprints (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
            tarray(1), tarray(2), tarray(1)+tarray(2)
    write(6,"(1x,'Elapsed time = ', e12.5)") tiempo
    stop
    end
!
!    ***************************************************************
!
  subroutine leegeomatdensGTO
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT), parameter :: mxcoefs = 5*mxcap
    integer(KINT) :: i, ia, ierr, indnf, indng, ios, j, k, k1, k2, knt, kntlm, l, la, m, nbas 
    integer(KINT), allocatable :: ngauss(:)
    real(KREAL) :: aux, bux, distmax, qe, qn, qt, zntot
    real(KREAL) :: centcharge(3), cfaux(0:2*mxl)
    logical :: existe
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
    
    allocate(qlm((mxmult+1)**2), qlmnuc((mxmult+1)**2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating qlm and qlmnuc. Stop')

    allocate(rcen(3,ncen), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating rcen. Stop')
    
    allocate(rmultip(25,ncen), stat = ierr)
    if (ierr .ne. 0) call error(1,'Memory error when allocating rmultip. Stop')

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
    
    
!    writes geometry in .xyz file for gOpenMol compatibility
!
    existe = .false.
    inquire(file=trim(projectname)//".xyz", exist=existe)
    if (.not. existe) then
        open(99,file=trim(projectname)//".xyz",form='formatted', iostat=ierr)
        if (ierr .ne. 0) then
            write(6,"('Cannot open file ', a)") trim(projectname)//".xyz"
        else
            write(99,*) ncen
            write(99,*)
            do i = 1 , ncen
                    write(99,"(a2,3f15.10)") atmnms( int(zn(i)) ), (rcen(j,i)*0.529177249d0,j=1,3)
            end do
            close(99)
        endif
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
    rmultip = cero
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
            la = ll(ncaps)
            if (la .lt. 5) then
                do m = -la, la
                    rmultip(la*(la+1)+m+1,ia) = rmultip(la*(la+1)+m+1,ia) + dos * pi * cfaux(la+m) * facts(la) &
                        / ((2*la+1) * sqrt(xxg(ncaps))**(la+la+3))
                enddo
            endif
        enddo
    enddo
    nbas = indnf-1

    if (lmaxbase .gt. mxl) then
        write(6,"('Basis functions with not allowed values of  l. ')")
        write(6,"('Highest allowed value: ', i2 , ' Highest value in basis set: ', i2)") mxl, lmaxbase
        call error(1,' Stop')
    endif

!    prints out the input data
    write(6,"(27x,'GEOMETRY')")
    write(6,"(/t1, ' no. of center:', t22, 'x', t34, 'y', t46, 'z', t58, 'charge', t70, 'n. of gaussians')")
    do ia = 1, ncen
        write(6,"(t6, i5, t15, f12.7, t27, f12.7, t39, f12.7, t53, f10.5, t75, i3)") &
                ia, rcen(1,ia), rcen(2,ia), rcen(3,ia) , zn(ia), ngfin(ia)-ngini(ia)+1
    enddo
    write(6,"(/27x,'GTO EXPANSION OF DENSITY')")
    if (longoutput) then
        knt = 0
        do ia = 1, ncen
            if (ngini(ia) .le. 0) cycle
            write(6,"(/1x,'atom no.',1x,i5,'(',a2,')')") ia, atmnam(ia)
            write(6,"(1x,'Expansion: ', t14, 'L', t18, 'exponent', t35, 'coefficients(-L:L) ')")
            do j = ngini(ia), ngfin(ia)
                write(6,"(t14,i1,t16,e12.5,t32,5(1x,e22.15))") ll(j), xxg(j), coefs(nf(j):nf(j)+2*ll(j))
            enddo
        enddo
    endif
    write(6,"(/'Number of functions for density expansion = ', i8)") nbas
    call totalchargeGTO
!	Computes the total molecular multipolar moments from the atomic moments directly computed
    call multmolec
    write(6,"(/12x, ' molecular multipole components (nuclei, electrons, total)', /1x, 100('='))")
    kntlm = 0
    do l = 0, mxmult
        write(6,"(/57x,'L = ', i2,/57x,6('='),/)") l
        do m = -l, l
            kntlm = kntlm + 1
            qn = qlmnuc(kntlm)
            if (abs(qn) .lt. 1.d-15) qn = cero
            qe = qlm(kntlm)
            if (abs(qe) .lt. 1.d-15) qe = cero
            qt = qlm(kntlm)+qlmnuc(kntlm)
            if (abs(qt) .lt. 1.d-15) qt = cero
            write(6,"(5x,'q(',i2,',',i3,'): ',8x,3(2x,d22.15))") l, m, qn, qe, qt
        enddo
    enddo
    deallocate (qlm, qlmnuc, rmultip)
    return
    end
!
!	********************************************************
!
!	Subroutine for computing molecular multipolar moments from the multipolar moments of the &
!	 fragments. Generated with trasladamultipolos2.nb with the procedure described in the &
!	 libro de electrostatica.
!
  subroutine multmolec
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    real(KREAL) :: rab2, xab, yab, zab
    integer(KINT) :: i, ia, l, m
    real(KREAL) :: zlm(0:mxmult,-mxmult:mxmult)
    do i = 1, (mxmult+1)**2
        qlm(i) = cero
        qlmnuc(i) = cero
    enddo
    do ia = 1, ncen
            xab = rcen(1,ia)
            yab = rcen(2,ia)
            zab = rcen(3,ia)
            rab2 = xab*xab+yab*yab+zab*zab
            zlm(0,0) = uno
            zlm(1,-1) = yab
            zlm(1,0) = zab
            zlm(1,1) = xab
            do l = 1, mxmult-1
                zlm(l+1,l+1) = dosl1(l)*(xab*zlm(l,l)-yab*zlm(l,-l))
                zlm(l+1,-(l+1)) = dosl1(l)*(yab*zlm(l,l)+xab*zlm(l,-l))
                zlm(l+1,l) = dosl1(l)*zab*zlm(l,l)
                zlm(l+1,-l) = dosl1(l)*zab*zlm(l,-l)
                do m = 0, l-1
                    zlm(l+1,m)=(dosl1(l)*zab*zlm(l,m)-re(l+m)*rab2*zlm(l-1,m))*ri(l-m+1)
                    zlm(l+1,-m)=(dosl1(l)*zab*zlm(l,-m)-re(l+m)*rab2*zlm(l-1,-m))*ri(l-m+1)
                enddo
            enddo
!	Contribution of the nuclear charge of center ia to the multipoles
            qlmnuc(1) = qlmnuc(1) + zn(ia)
            qlmnuc(2) = qlmnuc(2) + zlm(1,-1)*zn(ia)
            qlmnuc(3) = qlmnuc(3) + zlm(1,0)*zn(ia)
            qlmnuc(4) = qlmnuc(4) + zlm(1,1)*zn(ia)
            qlmnuc(5) = qlmnuc(5) + 0.0833333333333333D0*zlm(2,-2)*zn(ia)
            qlmnuc(6) = qlmnuc(6) + 0.3333333333333333D0*zlm(2,-1)*zn(ia)
            qlmnuc(7) = qlmnuc(7) + zlm(2,0)*zn(ia)
            qlmnuc(8) = qlmnuc(8) + 0.3333333333333333D0*zlm(2,1)*zn(ia)
            qlmnuc(9) = qlmnuc(9) + 0.0833333333333333D0*zlm(2,2)*zn(ia)
            qlmnuc(10) = qlmnuc(10) + 0.0027777777777778D0*zlm(3,-3)*zn(ia)
            qlmnuc(11) = qlmnuc(11) + 0.0166666666666667D0*zlm(3,-2)*zn(ia)
            qlmnuc(12) = qlmnuc(12) + 0.1666666666666667D0*zlm(3,-1)*zn(ia)
            qlmnuc(13) = qlmnuc(13) + zlm(3,0)*zn(ia)
            qlmnuc(14) = qlmnuc(14) + 0.1666666666666667D0*zlm(3,1)*zn(ia)
            qlmnuc(15) = qlmnuc(15) + 0.0166666666666667D0*zlm(3,2)*zn(ia)
            qlmnuc(16) = qlmnuc(16) + 0.0027777777777778D0*zlm(3,3)*zn(ia)
            qlmnuc(17) = qlmnuc(17) + 0.0000496031746032D0*zlm(4,-4)*zn(ia)
            qlmnuc(18) = qlmnuc(18) + 0.0003968253968254D0*zlm(4,-3)*zn(ia)
            qlmnuc(19) = qlmnuc(19) + 0.0055555555555556D0*zlm(4,-2)*zn(ia)
            qlmnuc(20) = qlmnuc(20) + 0.1000000000000000D0*zlm(4,-1)*zn(ia)
            qlmnuc(21) = qlmnuc(21) + zlm(4,0)*zn(ia)
            qlmnuc(22) = qlmnuc(22) + 0.1000000000000000D0*zlm(4,1)*zn(ia)
            qlmnuc(23) = qlmnuc(23) + 0.0055555555555556D0*zlm(4,2)*zn(ia)
            qlmnuc(24) = qlmnuc(24) + 0.0003968253968254D0*zlm(4,3)*zn(ia)
            qlmnuc(25) = qlmnuc(25) + 0.0000496031746032D0*zlm(4,4)*zn(ia)
            if (ngini(ia) .le. 0) cycle	! If center without associated basis set, cycles
!	Contribution of the electronic fragments to the multipoles
            qlm(1) = qlm(1) - rmultip(1,ia)
            qlm(2) = qlm(2) - (zlm(1,-1)*rmultip(1,ia)+rmultip(2,ia))
            qlm(3) = qlm(3) - (zlm(1,0)*rmultip(1,ia)+rmultip(3,ia))
            qlm(4) = qlm(4) - (zlm(1,1)*rmultip(1,ia)+rmultip(4,ia))
            qlm(5) = qlm(5) - (0.0833333333333333D0*zlm(2,-2)*rmultip(1,ia)+0.5000000000000000D0*zlm(1,1) &
                    *rmultip(2,ia)+0.5000000000000000D0*zlm(1,-1)*rmultip(4,ia)+rmultip(5,ia))
            qlm(6) = qlm(6) - (0.3333333333333333D0*zlm(2,-1)*rmultip(1,ia)+zlm(1,0)*rmultip(2,ia)+zlm(1,-1) &
                    *rmultip(3,ia)+rmultip(6,ia))
            qlm(7) = qlm(7) - (zlm(2,0)*rmultip(1,ia)-1.0000000000000000D0*zlm(1,-1)*rmultip(2,ia) &
                    +2.0000000000000000D0*zlm(1,0)*rmultip(3,ia)-1.0000000000000000D0*zlm(1,1)*rmultip(4,ia) &
                    +rmultip(7,ia))
            qlm(8) = qlm(8) - (0.3333333333333333D0*zlm(2,1)*rmultip(1,ia)+zlm(1,1)*rmultip(3,ia)+zlm(1,0) &
                    *rmultip(4,ia)+rmultip(8,ia))
            qlm(9) = qlm(9) - (0.0833333333333333D0*zlm(2,2)*rmultip(1,ia)-0.5000000000000000D0*zlm(1,-1) &
                    *rmultip(2,ia)+0.5000000000000000D0*zlm(1,1)*rmultip(4,ia)+rmultip(9,ia))
            qlm(10) = qlm(10) - (0.0027777777777778D0*zlm(3,-3)*rmultip(1,ia)+0.0416666666666667D0*zlm(2,2) &
                    *rmultip(2,ia)+0.0416666666666667D0*zlm(2,-2)*rmultip(4,ia)+0.5000000000000000D0*zlm(1,1) &
                    *rmultip(5,ia)+0.5000000000000000D0*zlm(1,-1)*rmultip(9,ia)+rmultip(10,ia))
            qlm(11) = qlm(11) - (0.0166666666666667D0*zlm(3,-2)*rmultip(1,ia)+0.1666666666666667D0*zlm(2,1) &
                    *rmultip(2,ia)+0.0833333333333333D0*zlm(2,-2)*rmultip(3,ia)+0.1666666666666667D0*zlm(2,-1) &
                    *rmultip(4,ia)+zlm(1,0)*rmultip(5,ia)+0.5000000000000000D0*zlm(1,1)*rmultip(6,ia) &
                    +0.5000000000000000D0*zlm(1,-1)*rmultip(8,ia)+rmultip(11,ia))
            qlm(12) = qlm(12) - (0.1666666666666667D0*zlm(3,-1)*rmultip(1,ia)+(zlm(2,0)+0.0833333333333333D0 &
                    *zlm(2,2))*rmultip(2,ia)+0.6666666666666666D0*zlm(2,-1)*rmultip(3,ia)-0.0833333333333333D0 &
                    *zlm(2,-2)*rmultip(4,ia)-1.0000000000000000D0*zlm(1,1)*rmultip(5,ia)+2.0000000000000000D0 &
                    *zlm(1,0)*rmultip(6,ia)+zlm(1,-1)*rmultip(7,ia)+zlm(1,-1)*rmultip(9,ia)+rmultip(12,ia))
            qlm(13) = qlm(13) - (zlm(3,0)*rmultip(1,ia)-1.0000000000000000D0*zlm(2,-1)*rmultip(2,ia) &
                    +3.0000000000000000D0*zlm(2,0)*rmultip(3,ia)-1.0000000000000000D0*zlm(2,1) &
                    *rmultip(4,ia)-3.0000000000000000D0*zlm(1,-1)*rmultip(6,ia)+3.0000000000000000D0*zlm(1,0) &
                    *rmultip(7,ia)-3.0000000000000000D0*zlm(1,1)*rmultip(8,ia)+rmultip(13,ia))
            qlm(14) = qlm(14) - (0.1666666666666667D0*zlm(3,1)*rmultip(1,ia)-0.0833333333333333D0*zlm(2,-2) &
                    *rmultip(2,ia)+0.6666666666666666D0*zlm(2,1)*rmultip(3,ia)+(zlm(2,0)-0.0833333333333333D0 &
                    *zlm(2,2))*rmultip(4,ia)-1.0000000000000000D0*zlm(1,-1)*rmultip(5,ia)+zlm(1,1)*rmultip(7,ia) &
                    +2.0000000000000000D0*zlm(1,0)*rmultip(8,ia)-1.0000000000000000D0*zlm(1,1)*rmultip(9,ia) &
                    +rmultip(14,ia))
            qlm(15) = qlm(15) - (0.0166666666666667D0*zlm(3,2)*rmultip(1,ia)-0.1666666666666667D0*zlm(2,-1) &
                    *rmultip(2,ia)+0.0833333333333333D0*zlm(2,2)*rmultip(3,ia)+0.1666666666666667D0*zlm(2,1) &
                    *rmultip(4,ia)-0.5000000000000000D0*zlm(1,-1)*rmultip(6,ia)+0.5000000000000000D0*zlm(1,1) &
                    *rmultip(8,ia)+zlm(1,0)*rmultip(9,ia)+rmultip(15,ia))
            qlm(16) = qlm(16) - (0.0027777777777778D0*zlm(3,3)*rmultip(1,ia)-0.0416666666666667D0*zlm(2,-2) &
                    *rmultip(2,ia)+0.0416666666666667D0*zlm(2,2)*rmultip(4,ia)-0.5000000000000000D0*zlm(1,-1) &
                    *rmultip(5,ia)+0.5000000000000000D0*zlm(1,1)*rmultip(9,ia)+rmultip(16,ia))
            qlm(17) = qlm(17) - (0.0000496031746032D0*zlm(4,-4)*rmultip(1,ia)+0.0013888888888889D0*zlm(3,3) &
                    *rmultip(2,ia)+0.0013888888888889D0*zlm(3,-3)*rmultip(4,ia)+0.0416666666666667D0*zlm(2,2) &
                    *rmultip(5,ia)+0.0416666666666667D0*zlm(2,-2)*rmultip(9,ia)+0.5000000000000000D0*zlm(1,1) &
                    *rmultip(10,ia)+0.5000000000000000D0*zlm(1,-1)*rmultip(16,ia)+rmultip(17,ia))
            qlm(18) = qlm(18) - (0.0003968253968254D0*zlm(4,-3)*rmultip(1,ia)+0.0083333333333333D0*zlm(3,2) &
                    *rmultip(2,ia)+0.0027777777777778D0*zlm(3,-3)*rmultip(3,ia)+0.0083333333333333D0*zlm(3,-2) &
                    *rmultip(4,ia)+0.1666666666666667D0*zlm(2,1)*rmultip(5,ia)+0.0416666666666667D0*zlm(2,2) &
                    *rmultip(6,ia)+0.0416666666666667D0*zlm(2,-2)*rmultip(8,ia)+0.1666666666666667D0*zlm(2,-1) &
                    *rmultip(9,ia)+zlm(1,0)*rmultip(10,ia)+0.5000000000000000D0*zlm(1,1)*rmultip(11,ia) &
                    +0.5000000000000000D0*zlm(1,-1)*rmultip(15,ia)+rmultip(18,ia))
            qlm(19) = qlm(19) - (0.0055555555555556D0*zlm(4,-2)*rmultip(1,ia)+(0.0833333333333333D0*zlm(3,1) &
                    +0.0027777777777778D0*zlm(3,3))*rmultip(2,ia)+0.0333333333333333D0*zlm(3,-2)*rmultip(3,ia) &
                    +(0.0833333333333333D0*zlm(3,-1)-0.0027777777777778D0*zlm(3,-3))*rmultip(4,ia)+zlm(2,0) &
                    *rmultip(5,ia)+0.3333333333333333D0*zlm(2,1)*rmultip(6,ia)+0.0833333333333333D0*zlm(2,-2) &
                    *rmultip(7,ia)+0.3333333333333333D0*zlm(2,-1)*rmultip(8,ia)-1.0000000000000000D0*zlm(1,1) &
                    *rmultip(10,ia)+2.0000000000000000D0*zlm(1,0)*rmultip(11,ia)+0.5000000000000000D0*zlm(1,1) &
                    *rmultip(12,ia)+0.5000000000000000D0*zlm(1,-1)*rmultip(14,ia)+zlm(1,-1)*rmultip(16,ia) &
                    +rmultip(19,ia))
            qlm(20) = qlm(20) - (0.1000000000000000D0*zlm(4,-1)*rmultip(1,ia)+(zlm(3,0)+0.0500000000000000D0 &
                    *zlm(3,2))*rmultip(2,ia)+0.5000000000000000D0*zlm(3,-1)*rmultip(3,ia)-0.0500000000000000D0 &
                    *zlm(3,-2)*rmultip(4,ia)-1.0000000000000000D0*zlm(2,1)*rmultip(5,ia)+(3.0000000000000000D0 &
                    *zlm(2,0)+0.2500000000000000D0*zlm(2,2))*rmultip(6,ia)+zlm(2,-1) &
                    *rmultip(7,ia)-0.2500000000000000D0*zlm(2,-2)*rmultip(8,ia)+zlm(2,-1) &
                    *rmultip(9,ia)-3.0000000000000000D0*zlm(1,1)*rmultip(11,ia)+3.0000000000000000D0*zlm(1,0) &
                    *rmultip(12,ia)+zlm(1,-1)*rmultip(13,ia)+3.0000000000000000D0*zlm(1,-1)*rmultip(15,ia) &
                    +rmultip(20,ia))
            qlm(21) = qlm(21) - (zlm(4,0)*rmultip(1,ia)-1.0000000000000000D0*zlm(3,-1)*rmultip(2,ia) &
                    +4.0000000000000000D0*zlm(3,0)*rmultip(3,ia)-1.0000000000000000D0*zlm(3,1)*rmultip(4,ia) &
                    +zlm(2,-2)*rmultip(5,ia)-4.0000000000000000D0*zlm(2,-1)*rmultip(6,ia)+6.0000000000000000D0 &
                    *zlm(2,0)*rmultip(7,ia)-4.0000000000000000D0*zlm(2,1)*rmultip(8,ia)+zlm(2,2) &
                    *rmultip(9,ia)-6.0000000000000000D0*zlm(1,-1)*rmultip(12,ia)+4.0000000000000000D0*zlm(1,0) &
                    *rmultip(13,ia)-6.0000000000000000D0*zlm(1,1)*rmultip(14,ia)+rmultip(21,ia))
            qlm(22) = qlm(22) - (0.1000000000000000D0*zlm(4,1)*rmultip(1,ia)-0.0500000000000000D0*zlm(3,-2) &
                    *rmultip(2,ia)+0.5000000000000000D0*zlm(3,1)*rmultip(3,ia)+(zlm(3,0)-0.0500000000000000D0 &
                    *zlm(3,2))*rmultip(4,ia)-1.0000000000000000D0*zlm(2,-1)*rmultip(5,ia)-0.2500000000000000D0 &
                    *zlm(2,-2)*rmultip(6,ia)+zlm(2,1)*rmultip(7,ia)+(3.0000000000000000D0 &
                    *zlm(2,0)-0.2500000000000000D0*zlm(2,2))*rmultip(8,ia)-1.0000000000000000D0*zlm(2,1) &
                    *rmultip(9,ia)-3.0000000000000000D0*zlm(1,-1)*rmultip(11,ia)+zlm(1,1)*rmultip(13,ia) &
                    +3.0000000000000000D0*zlm(1,0)*rmultip(14,ia)-3.0000000000000000D0*zlm(1,1)*rmultip(15,ia) &
                    +rmultip(22,ia))
            qlm(23) = qlm(23) - (0.0055555555555556D0*zlm(4,2)*rmultip(1,ia)+(-0.0833333333333333D0 &
                    *zlm(3,-1)-0.0027777777777778D0*zlm(3,-3))*rmultip(2,ia)+0.0333333333333333D0*zlm(3,2) &
                    *rmultip(3,ia)+(0.0833333333333333D0*zlm(3,1)-0.0027777777777778D0*zlm(3,3)) &
                    *rmultip(4,ia)-0.3333333333333333D0*zlm(2,-1)*rmultip(6,ia)+0.0833333333333333D0*zlm(2,2) &
                    *rmultip(7,ia)+0.3333333333333333D0*zlm(2,1)*rmultip(8,ia)+zlm(2,0) &
                    *rmultip(9,ia)-1.0000000000000000D0*zlm(1,-1)*rmultip(10,ia)-0.5000000000000000D0*zlm(1,-1) &
                    *rmultip(12,ia)+0.5000000000000000D0*zlm(1,1)*rmultip(14,ia)+2.0000000000000000D0*zlm(1,0) &
                    *rmultip(15,ia)-1.0000000000000000D0*zlm(1,1)*rmultip(16,ia)+rmultip(23,ia))
            qlm(24) = qlm(24) - (0.0003968253968254D0*zlm(4,3)*rmultip(1,ia)-0.0083333333333333D0*zlm(3,-2) &
                    *rmultip(2,ia)+0.0027777777777778D0*zlm(3,3)*rmultip(3,ia)+0.0083333333333333D0*zlm(3,2) &
                    *rmultip(4,ia)-0.1666666666666667D0*zlm(2,-1)*rmultip(5,ia)-0.0416666666666667D0*zlm(2,-2) &
                    *rmultip(6,ia)+0.0416666666666667D0*zlm(2,2)*rmultip(8,ia)+0.1666666666666667D0*zlm(2,1) &
                    *rmultip(9,ia)-0.5000000000000000D0*zlm(1,-1)*rmultip(11,ia)+0.5000000000000000D0*zlm(1,1) &
                    *rmultip(15,ia)+zlm(1,0)*rmultip(16,ia)+rmultip(24,ia))
            qlm(25) = qlm(25) - (0.0000496031746032D0*zlm(4,4)*rmultip(1,ia)-0.0013888888888889D0*zlm(3,-3) &
                    *rmultip(2,ia)+0.0013888888888889D0*zlm(3,3)*rmultip(4,ia)-0.0416666666666667D0*zlm(2,-2) &
                    *rmultip(5,ia)+0.0416666666666667D0*zlm(2,2)*rmultip(9,ia)-0.5000000000000000D0*zlm(1,-1) &
                    *rmultip(10,ia)+0.5000000000000000D0*zlm(1,1)*rmultip(16,ia)+rmultip(25,ia))
    enddo
    return
    end
!
!    ***************************************************************
!
!    Subroutine for expanding the density of a molecule in orthogonal polynomials times spherical harmonics 
!    centered at origin
! 
  subroutine expand(clmn_subr, frad_fun)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, i1, i12p, i1p, i1pini, i1pfin, i2, i2p, i2pini, i2pfin, ia, ib, ierr, ii, ir, irinf, irshft, irsup
    integer(KINT) :: ishift, j, k, klmi, kml, knt, knticf, kntlm
    integer(KINT) :: l, la, ldmrot, lk, lm, lma, lmax, lmin, m, ma
    integer(KINT) :: n, na, nfa
    real(KREAL) :: amodule, aux, bpmodule, bpx, bpy, bpz, bux, bx, bx2, bz, cosalfa, cosbeta, cosgamma, cux
    real(KREAL) :: den, dosx, exa, exb, factor, pmodule, prdesc
    real(KREAL) :: rab, rabinv, rdif, rn, rna, rnab, rnb, rp2, sinalfa, sinbeta, singamma, suma
    real(KREAL) :: x, x12inv, xa, xab, xb, xinv, xip, xjp, xkp, xy, ya, yab, yb, yip, yjp, ykp, za, zab, zb, zip, zjp, zkp
    logical :: lf0, lprint, lmultmod, lrotar
    real(KREAL) :: roaux(-mxl:mxl), qlm1c(0:mxldst), vinvariant(0:kexpansion)
    real(KREAL4) :: tarray(2), tiempo, tarraynw(2), tiemponw, dtime
    real(KREAL), allocatable :: fcent(:,:)
    interface
        subroutine clmn_subr
        end subroutine clmn_subr
    end interface
    interface
        subroutine frad_fun
        end subroutine frad_fun
    end interface
    tiempo = dtime(tarray)
    ldimaux = lexpansion + lmaxbase
    ldmrot = max(lmaxbase,lexpansion)
!    Allocates memory for arrays dl, fa, fcent, flm, flmmamb, ftot, rl, rlt
    allocate(dl(-ldmrot:ldmrot,-ldmrot:ldmrot,0:ldmrot), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating dl. Stop')
    allocate(fcent(nquadpoints,(lexpansion+1)**2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating fcent. Stop')
    allocate(flm(nquadpoints,-lmaxbase:lmaxbase,(lexpansion+1)**2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating flm. Stop')
    allocate(flmmamb((ldimaux+1)*(ldimaux+1),-lmaxbase:lmaxbase), stat = ierr)  ! allocates input data arrays
    if (ierr .ne. 0) call error(2,'Memory error when allocating flmmamb. Stop')
    allocate(ftot(nquadpoints,(lexpansion+1)**2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ftot. Stop')
    allocate(gkl(-1:kexpansion), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating gkl. Stop')
    allocate(radfunction(nquadpoints, 0:kexpansion), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating radfunction. Stop')
    allocate(rl(-ldmrot:ldmrot,-ldmrot:ldmrot,0:ldmrot), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating rl. Stop')
    allocate(rlt((ldmrot+1)*(2*ldmrot+1)*(2*ldmrot+3)/3), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating rlt. Stop')
    allocate(rquad01(nquadpoints), rquadscal(nquadpoints), rquadaux(nquadpoints), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating rquad01, rquadscal and rquadaux. Stop')
    allocate(weights01(nquadpoints), weights(nquadpoints), vaux(nquadpoints), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating weights01, weights and weightaux. Stop')
    

!     Chebyshev quadrature rule:     rquad01: abscissae in the interval [0,1]
!                             rquadscal: abscissae in the interval [0,r^*]
!                            weights01: weights in the interval [0,1] 
!                             weights: weights in the interval [0,r^*]
    do i = 1, nquadpoints
        rquad01(i) = umed * ( uno - cos( pi*(dble(i)-umed) / dble(nquadpoints) ) )
        weights01(i) = pi * sqrt(1.d0 - (2.d0*rquad01(i)-1.d0)*(2.d0*rquad01(i)-1.d0)) / (2.d0*nquadpoints)
    enddo
!
!    Legendre quadrature rules (available up to order 120)
!     call qgleg(nquadpoints, rquad01, weights01, 0.d0, 1.d0, 0.d0)

    rquadscal = rstar * rquad01

!     if (longoutput) then
!         write(6,"('quadrature abscissae in the interval [0,1]: ', 8(1x,e17.10))") rquad01
!         write(6,"('quadrature abscissae in the interval [0,r^*]: ', 8(1x,e17.10))") rquadscal
!         write(6,"('quadrature weights in the interval [0,1]: ', 8(1x,e17.10))") rquad01
!     endif

    ftot = cero
    do ia = 1, ncen    ! Do over centers (ia)
        if (ngini(ia) .le. 0) cycle    ! If center without associated basis set, cycles
        xa = rcen(1,ia)
        ya = rcen(2,ia)
        za = rcen(3,ia)

!     Computes Euler angles to orient the frame so that A center lies on the local z axis and center B lies on the local xz plane

        pmodule = sqrt(xa*xa + ya*ya)
        amodule = sqrt(xa*xa + ya*ya + za*za)
        if (pmodule .gt. 1.d-15) then
            cosalfa = xa / pmodule
            sinalfa = ya / pmodule
        else
            cosalfa = uno
            sinalfa = cero
        endif
        if (amodule .gt. 1.d-15) then
            cosbeta = za / amodule
            sinbeta = pmodule / amodule
        else
            cosbeta = uno
            sinbeta = cero
        endif
        cosgamma = uno
        singamma = cero
        lrotar = .false.
        if (((max(lmaxc(ia), lexpansion)) .gt. 0) .and. ( (cosalfa  .ne. uno) .or. (cosbeta  .ne. uno) ) ) then
            lrotar = .true.
            call rotar (ldmrot, cosalfa, sinalfa, cosbeta, sinbeta, cosgamma, singamma)

!    Stores the transpose of rl multiplied and divided by the angular normalization to give radial factors which multiply UNNORMALIZED spherical harmonics

            knt = 0
            do i = 0, ldmrot
                do k = -i, i
                    do j = -i, i
                        knt = knt + 1
                        rlt(knt) = rl(k,j,i) * ang(ind(i)+abs(k)+1) / ang(ind(i)+abs(j)+1)
                    enddo
                enddo
            enddo
        endif

!       Loop over the expansion functions
! write(6,*) ' lrotar = ', lrotar
        fcent = cero    ! array for accumulating radial factors corresponding to center A
        do i1 = ngini(ia), ngfin(ia)
            la = ll(i1)
            nfa = nf(i1)
!    Reads the pertinent block of expansion coefficients and rotates it to the aligned system. Loads the result in roblk.
!    Angular normalization factors are introduced at the end of the loading process.
            if (lrotar) then
                do i = -la, la
                    roaux(i) = coefs(i+la+nfa)
                enddo
!               Rotation on center A and introduction of the angular normalization
                do i = -la, la
!                     roblk(i) = ang(ind(la)+abs(i)+1) * dot_product(coefs(nfa:nfa+2*la), rl(-la:la,i,la))
                    roblk(i) = dot_product(coefs(nfa:nfa+2*la), rl(-la:la,i,la))
                enddo
            else
!               Introduction of the angular and radial normalization in case of no rotation
                do i = -la, la
!                     roblk(i) = ang(ind(la)+abs(i)+1) * coefs(i+la+nfa)
                    roblk(i) = coefs(i+la+nfa)
                enddo
            endif

!           Computes the radial factors associated to the translation of the distribution

            call fradAgauss(i1, amodule)

!                    Multiplies the radial factors by the density matrix and accumulates

            do i = -la, la
                aux = roblk(i)
                do lm = 1, (lexpansion+1) * (lexpansion+1)
                    fcent(1:nquadpoints,lm) = fcent(1:nquadpoints,lm) + aux * flm(1:nquadpoints,i,lm)
                enddo
            enddo
        enddo    ! End of Do over shells in ia
        
!    Rotates the radial factors back to the molecular axis system and accumulates
!
!    IMPORTANT !!!     The radial factors computed here are those multiplying UNNORMALIZED real spherical harmonics in the
!                expansion of the density
!
        if (lrotar) then
            lm = 0
            kml = 0
            do l = 0, lexpansion
                do m = -l, l
                    lm = lm + 1
                    lk = l*(l+1)+1
                    do k = -l, l
                        kml = kml + 1
                        ftot(1:nquadpoints,lm) = ftot(1:nquadpoints,lm) + fcent(1:nquadpoints,lk+k) * rlt(kml)
                    enddo
                enddo
            enddo
        else
            lm = 0
            do l = 0, lexpansion
                do m = -l, l
                    lm = lm + 1
                    ftot(1:nquadpoints,lm) = ftot(1:nquadpoints,lm) + fcent(1:nquadpoints,lm)
                enddo
            enddo
        endif  
    enddo    ! End of Do over centers (ia)
    tiempo = dtime(tarray)
    write(6,"(1x,//80('='),/'Timing in seconds for radial factors (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
            tarray(1), tarray(2), tarray(1)+tarray(2)
    write(6,"(1x,'Elapsed time = ', e12.5)") tiempo
    tiemponw = dtime(tarraynw)
    deallocate(dl, fcent, flm, flmmamb, rl, rlt)

!     Prints the valus of radial factors in quadrature points (unseal for testing)
!     do lm = 1, (lexpansion+1)*(lexpansion+1)
!         if (abs(ftot(1,lm)) .gt. 1.d-10) &
!             write(6,"('ftot(',i3,') = ', 8(1x,e17.10))") lm, ftot(1:nquadpoints,lm)
!     enddo

!     Computes the multipolar moments from the radial factors (unseal for testing)
!     lm = 0
!     weights = weights01 * rstar
!     do l = 0, lexpansion
!         weights = weights * rquadscal * rquadscal
!         do m = -l, l
!             lm = lm + 1
!             aux = 4.d0 * pi * dot_product(weights, ftot(1:nquadpoints,lm)) * ri(l+l+1)
!             if (abs(aux) .gt. 1.d-15) write(6,"('qlm(',i2,',',i3,') = ',e22.15)") l, m, aux
!         enddo
!     enddo

!     Computes Zernike 3D or Jacobi moments (integrals with normalized Zernike 3D or Jacobi functions)

    allocate(omeganlm(0:kexpansion,(lexpansion+1)*(lexpansion+1)), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating omeganlm. Stop')

    omeganlm = cero
    call clmn_subr

!     Computes molecular multipoles from Zernike 3D or Jacobi expansions

    call multipoles

    tiemponw = dtime(tarraynw)
    write(6,"(1x,//80('='),/'Timing in seconds for projection (user, system, total):',/5x,'(',e12.5,',',e12.5,',',e12.5')')") &
            tarraynw(1), tarraynw(2), tarraynw(1)+tarraynw(2)
    write(6,"(1x,'Elapsed time = ', e12.5)") tiemponw

!     Writes Zernike 3D or Jacobi moments to .zernike or .jacobi file

    if (ljacobi) then
        open(99,file=trim(projectname)//".jacobi",form='formatted', iostat=ierr)
        if (ierr .ne. 0) call error(ierr,'Cannot open file '//trim(projectname)//'.jacobi. Stop')
    else
        open(99,file=trim(projectname)//".zernike",form='formatted', iostat=ierr)
        if (ierr .ne. 0) call error(ierr,'Cannot open file '//trim(projectname)//'.zernike. Stop')
    endif

    write(99,*) rstar
    write(99,*) lexpansion, kexpansion
    do l = 0, lexpansion
        do m = -l, l
            write(99,"(8(1x,e22.15))") omeganlm(0:kexpansion,l*(l+1)+m+1)
        enddo
    enddo

!     Computes and prints fingerprints (rotationally invariant) associated to Zernike 3D or Jacobi expansions
    if (ljacobi) then
        write(6,"(/' Rotationally invariant fingerprints for Jacobi expansion: sqrt( sum_{m=-l}^l |c_kl^m|^2 )',/91(1H=),/)")
    else
        write(6,"(/' Rotationally invariant fingerprints for Zernike 3D expansion: sqrt( sum_{m=-l}^l |c_kl^m|^2 )',/95(1H=),/)")
    endif
    tiemponw = dtime(tarraynw)
    lm = 0
    n = kexpansion
    do l = 0, lexpansion
        vinvariant = 0.d0
        do m = -l, l
            lm = lm + 1
            vinvariant = vinvariant + omeganlm(0:n,lm)*omeganlm(0:n,lm)
        enddo
        write(6,"('l = ', i3, ': ', 8(1x,e17.10))") l, sqrt(vinvariant)
    enddo

!     Prints data to file projectname.mth for checking (for instance, with Mathematica) 

    if (lmthfile) then
        if (ljacobi) then
            open(17,file=trim(projectname)//"-Jacobi.mth",form='formatted', iostat=ierr)
        else
            open(17,file=trim(projectname)//"-Zernike.mth",form='formatted', iostat=ierr)
        endif
        if (ierr .ne. 0) call error(ierr,'Cannot open file '//trim(projectname)//'.mth. Stop')
        if (ljacobi) then
            i = 1
        else
            i = 0
        endif
        write(17,*) i, rstar, lexpansion, kexpansion, nquadpoints
        write(17,"(10(1x,e22.15))") weights01(1:nquadpoints) * rstar
        write(17,"(10(1x,e22.15))") rquadscal(1:nquadpoints)
        lm = 0
        do l = 0, lexpansion
            do m = -l, l
                lm = lm + 1
                write(17,"(10(1x,e22.15))") ftot(1:nquadpoints,lm) / ang(ind(l)+abs(m)+1)
            enddo
        enddo
        do l = 0, lexpansion
            do m = -l, l
                write(17,"(10(1x,e22.15))") omeganlm(0:kexpansion,l*(l+1)+m+1)
            enddo
        enddo
        close(17)
    endif
    return
    end
!
!    ***************************************************************
!
!     Computes the Zernike 3D Moments (expansion coefficients of radial factors flm in Zernike 3D functions )
!
  subroutine omeganlm_zernike
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: k, l, m, n
    weights =  sqrt(rstar) * weights01
    rquadaux = rstar * rquad01 * rquad01
    n = kexpansion
    do l = 0, lexpansion    ! Loop over Zernike's functions
        weights = weights * rquadaux
        call zernike3DR(l)
        if (lechelon) n = kexpansion - l
        do k = 0, n
            vaux = weights * radfunction(1:nquadpoints,k)
            do m = -l, l
                omeganlm(k,l*(l+1)+m+1) = dot_product(vaux,ftot(1:nquadpoints, l*(l+1)+m+1)) / ang(ind(l)+abs(m)+1)
            enddo
        enddo
    enddo
    return
    end
!
!    ***************************************************************
!
!     Computes the Zernike 3D expansions of radial factors in the quadrature points
!
  subroutine frad_zernike
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, k, l, m, n
    ftot = 0.d0
    rquadaux = 1.d0 / (rstar*sqrt(rstar))
    n = kexpansion
    do l = 0, lexpansion    ! Loop over Zernike's functions
        call zernike3DR(l)
        if (lechelon) n = kexpansion - l
        do k = 0, n
            vaux = rquadaux * radfunction(1:nquadpoints,k)
            do m = -l, l
                do i = 1, nquadpoints
                    ftot(i,l*(l+1)+m+1) = ftot(i,l*(l+1)+m+1) + omeganlm(k,l*(l+1)+m+1) * vaux(i) * ang(ind(l)+abs(m)+1)
                enddo
            enddo
        enddo
        rquadaux = rquadaux * rquad01
    enddo
    return
    end
!
!    ***************************************************************
!
!     Computes the radial part of Zernike 3D functions divided by (r/r*)^l in the quadrature points
! 
  subroutine zernike3DR(l)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, jshift, k, kshift, l, n
    real(KREAL) :: rqsq(nquadpoints)
    rqsq = rquad01 * rquad01
    jshift = (kexpansion+1)*l+1
    kshift = 2*l+7
    if (lechelon) then
        n = kexpansion - l
    else
        n = kexpansion
    endif
    do i = 1, nquadpoints
        gkl(-1) = 0.d0
        gkl(0) = 1.d0
        radfunction(i,0) = root(2*l+3)
        do k = 0, n-1
            gkl(k+1) = (cfgkl1(jshift+k) - cfgkl2(jshift+k) * rqsq(i)) * gkl(k) - cfgkl3(jshift+k) * gkl(k-1)
            radfunction(i,k+1) = akgkl(k+1) * root(kshift+4*k) * gkl(k+1)
        enddo
    enddo
    return
    end
!
!    ***************************************************************
!
!     Computes the Jacobi Moments (expansion coefficients of radial factors flm in Jacobi P(0,2+2l) polynomials)
!
  subroutine omeganlm_jacobi
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: k, l, m, n
    weights =  sqrt(rstar) * weights01
    rquadaux = rstar * rquad01 * rquad01
    n = kexpansion
    do l = 0, lexpansion    ! Loop over Zernike's functions
        weights = weights * rquadaux
        call jacobiP(l)
        if (lechelon) n = kexpansion - l
        do k = 0, n
            vaux = weights * radfunction(1:nquadpoints,k)
            do m = -l, l
                omeganlm(k,l*(l+1)+m+1) = dot_product(vaux,ftot(1:nquadpoints, l*(l+1)+m+1)) / ang(ind(l)+abs(m)+1)
            enddo
        enddo
    enddo
    return
    end
!
!    ***************************************************************
!
!     Computes the Jacobi ( P(0,2+2l) ) expansions of radial factors in the quadrature points
!
  subroutine frad_jacobi
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, k, l, m, n
    real(KREAL) :: aux
    ftot = 0.d0
    rquadaux = 1.d0 / (rstar*sqrt(rstar))
    n = kexpansion
    do l = 0, lexpansion    ! Loop over Zernike's functions
        call jacobiP(l)
        if (lechelon) n = kexpansion - l
        do k = 0, n
            vaux = rquadaux * radfunction(1:nquadpoints,k)
            do m = -l, l
                do i = 1, nquadpoints
                    ftot(i,l*(l+1)+m+1) = ftot(i,l*(l+1)+m+1) + omeganlm(k,l*(l+1)+m+1) * vaux(i) * ang(ind(l)+abs(m)+1)
                enddo
            enddo
        enddo
        rquadaux = rquadaux * rquad01
    enddo
    return
    end
!
!**********************************************************************
! 
!     Computes the Jacobi polynomials P(0,2+2l)(2t-1)
!
  subroutine jacobiP(l)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: k, l, n
    real(KREAL) :: t
    if (lechelon) then
        n = kexpansion - l
    else
        n = kexpansion
    endif
    radfunction(1:nquadpoints,0) = 1.d0
    radfunction(1:nquadpoints,1) = 4.d0 * rquad01 - 3.d0 + dble(l+l) * (rquad01 - 1.d0)
    do k = 1, n-1
        radfunction(1:nquadpoints,k+1) = (dble(2*k+2*l+3) * (dble((k+l+1)*(k+l+2)) * (2.d0*rquad01-1.d0) - dble((l+1)*(l+1))) &
                * radfunction(1:nquadpoints,k) - dble(k*(k+l+2)*(k+2*l+2))*radfunction(1:nquadpoints,k-1)) &
                / dble((k+1)*(k+l+1)*(k+2*l+3))
    enddo
    do k = 0, n
        radfunction(1:nquadpoints,k) = radfunction(1:nquadpoints,k) * root(2*(k+l)+3)
    enddo
    return
    end
!
!    ***************************************************************
!
! Computes the multipolar moments from the Zernike 3D expansion or Jacobi expansion ( P^(0,2+2l) ). 
! Only the first Zernike 3D or Jacobi polynomial (n = l) of each  (l,m)  
! contributes to the corresponding multipolar moment, the remaining ones integrate to zero.
! Multipolar moments defined here as the coefficients that multiply the unnormalized irregular harmonics in the long-range 
! expansion of the electrostatic potential, q(l,m), are displayed in the third column.
! Rotationally invariant multipole moments: q(l,m) * sqrt((1+delta(m,0)) *  (l+|m|)!/(l-|m|)! ) 
! are listed in fourth column.
! Modules of the rotationally invariant multipole moments are displayed in the fifth column for each l
!
  subroutine multipoles
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: l, lm, m, n
    real(KREAL) :: aux, bux, buxmod
    write(6,"(//'Non-vanishing multipole components computed from the expansion',/63(1H=),&
            &//t5,'l', t10, 'm', t23,'q(l,m)',t45,'q(l,m) * sqrt((1+delta(m,0)) (l+|m|)!/(l-|m|)!)',t98,'qlmod'/)")
    lm = 0
    do l = 0, lexpansion
        buxmod = 0.d0
        do m = -l, l
            lm = lm + 1
            aux = omeganlm(0,lm) * fact(l-abs(m)) * facti(l+abs(m)) * sqrt(rstar**(2*l+3) * ri(2*l+3)) &
                    / (ang(ind(l)+abs(m)+1))
            if (m .ne. 0) aux = aux + aux
            bux = aux * sqrt(fact(l+abs(m)) * facti(l-abs(m)))
            if (m .eq. 0) bux = raiz2 * bux
            buxmod = buxmod + bux * bux
            if (m .lt. l) then
                if (abs(aux) .gt. thresmult) write(6,"(t4,i2,t8,i3,t16,1x,e22.15,15x,e22.15)") l, m, aux, bux
            else
                write(6,"(t4,i2,t8,i3,t16,1x,e22.15,15x,e22.15,14x,e22.15)") l, m, aux, bux, sqrt(buxmod)
            endif
        enddo
    enddo
    return
    end

!**********************************************************************
!    subroutine consta
!
!    Computes and stores auxiliary constants
!        re(i) = dfloat(i)
!        r1(i) = 1.d0 / dfloat(i)
!        fact(i) = dfloat(i!)
!        facti(i) = 1.d0 / dfloat(i!)
!        ind(i) = i*(i+1)/2
!        root(i) = dfloat(sqrt(i))
!        rooti(i) = 1.d0 / dfloat(sqrt(i))
!        ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) 
!            / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
!
!**********************************************************************
  subroutine consta
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, ierr, ikt, ip, j, k, k1, k12, knt, kntlm, l, l1, l1l1, l2, l2l2, la, lb, lm, lp
    integer(KINT) :: m, m1, m1a, m2, m2a, ma, mb, md, ms, mxang, n, np, nu
    real(KREAL) :: aux, aux1, aux2, auxk, auxm, bux, cux, sd, sgn, ss
!    auxiliary parameters and functions
    pi = acos(-uno)
    raizpi = sqrt(pi)
    pimed = umed * pi
    pimedsqr = sqrt(pimed)   ! Sqrt[pi/2]
    re(0) = cero
    ri(0) = 1.d300
    dosl1(0) = uno
    dosl1i(0) = uno
    do i = 1, mxreal
        re(i) = re(i-1) + uno        ! dfloat(i)
        re(-i) = -re(i)
        ri(i) = uno / re(i)           ! uno / dfloat(i)
        ri(-i) = -ri(i)
        dosl1(i) = re(i) + re(i) + uno    ! dfloat(iì)
        dosl1(-i) = -re(i) - re(i) + uno
        dosl1i(i) = uno / dosl1(i)        ! dfloat( 1/(i+i+1) )
        dosl1i(-i) = uno / dosl1(-i)
    enddo
    fact(0) = uno
    facti(0) = uno
    facts(-1) = raizpi
    facts(0) = facts(-1) * umed
    factsi(-1) = uno / facts(-1)    ! factsi in module GAUSS
    factsi(0) = uno / facts(0)
    do i = 1, mxfact
        fact(i) = fact(i-1) * re(i)               !  i!
        facts(i) = facts(i-1) * re(i+i+1) * umed    ! (i+1/2)!
        facti(i) = uno / fact(i)                    !  uno / i!
        factsi(i) = uno / facts(i)                !  uno / (i+1/2)!
    enddo
    mxind = (mxltot+1)*(mxltot+2)/2
    allocate(ind(0:mxind), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ind in consta. Stop')
    ind(0) = 0
    do i = 1, mxind
        ind(i) = ind(i-1) + i         !  i*(i+1)/2
    enddo
    root(0) = cero
    do i = 1, mxroot
        root(i) = sqrt(re(i))        !  sqrt(i)
        rooti(i) = uno / root(i)     !  uno / sqrt(i)
    enddo

    mxang = max(lexpansion,mxl)
    allocate(ang((mxang+1)*(mxang+2)/2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ang in consta. Stop')
!    ang(l*(l+1)/2+m+1) = sqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
    ang(1) = umed / raizpi
    lm = 1
    do l = 1, mxang
        lm = lm + 1
        ang(lm) = ang(1) * sqrt(re(2*l+1))
        aux = ang(lm) * raiz2
        do m = 1, l
            lm = lm + 1
            aux = aux / sqrt(re(l-m+1)*re(l+m))
            ang(lm) = aux
        enddo
    enddo

!     Computes and stores the coefficients for recursion of Zernike 3D functions
    n = (kexpansion+1)*(lexpansion+1)
    allocate(akgkl(0:kexpansion), cfgkl1(n), cfgkl2(n), cfgkl3(n), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating akgkl, cfgkl1, cfgkl2 and cfgkl3. Stop')
    akgkl(0) = 1.d0        ! akgkl(k) = (-1)^k (k-1/2)! / (sqrt(pi) * k!)
    do k = 0, kexpansion-1
        akgkl(k+1) = - re(2*k+1) * ri(2*k+2) * akgkl(k)
    enddo
    knt = 0
    do l = 0, lexpansion
        do k = 0, kexpansion
            knt = knt + 1
            cfgkl1(knt) = re(4*k+2*l+3) * dble(4*k*(2*k+2*l+3)+4*l*(l+2)+3) &
                    * ri(2*k+2*l+3) * ri(2*k+1) * ri(4*k+2*l+1)
            cfgkl2(knt) = re(4*k+2*l+3) * re(4*k+2*l+5) * ri(2*k+2*l+3) * ri(2*k+1)
            cfgkl3(knt) = dble(4*k*k) * re(4*k+2*l+5) * re(2*k+2*l+1) &
                    * ri(4*k+2*l+1) * ri(2*k+2*l+3) * ri(2*k-1) * ri(2*k+1)
        enddo
    enddo

!     Tabulates the coefficients for the decomposition of products
!     of two functions depending on phi (sin (m*phi), cos (m*phi))
!     into functions of the same type
    mxemes = mxltot
    allocate(ssv(-mxemes:mxemes,-mxemes:mxemes), sdv(-mxemes:mxemes,-mxemes:mxemes), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating ssv and sdv in consta. Stop')
    allocate(msv(-mxemes:mxemes,-mxemes:mxemes), mdv(-mxemes:mxemes,-mxemes:mxemes), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating msv and mdv in consta. Stop')
    do m2 = -mxemes, mxemes
        do m1 = -mxemes, mxemes
            call emes ( m1, m2, ms, md, ss, sd )
            msv(m1,m2) = ms
            mdv(m1,m2) = md
            ssv(m1,m2) = ss
            sdv(m1,m2) = sd
        enddo
    enddo
!    Coefficients for the decomposition of products of regular spherical harmonics into
!    regular spherical harmonics
    mxlcof = mxltot*(mxltot+3)/2
    mxkcof = mxlcof*(mxlcof+3)/2
    allocate(app(0:2*mxltot+1,0:mxkcof), stat = ierr)
    if (.not. allocated(app)) call error(ierr,'Memory error when allocating app. Stop')
    if (longoutput) write(6,"('Size of app   = ', i15, ' bytes')") size(app)
    allocate(bpp(0:2*mxltot+1,0:mxkcof), stat = ierr)
    if (.not. allocated(bpp)) call error(ierr,'Memory error when allocating bpp. Stop')
    if (longoutput) write(6,"('Size of bpp   = ', i15, ' bytes')") size(bpp)
    call acof
    call bcof
!    Tabulates some auxiliary indices for locating the previous coefficients
    allocate(indk12((mxltot+1)**2,(mxltot+1)**2), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating indk12 in consta. Stop')
    do l2 = 0,mxltot
        do l1 = 0,mxltot
            do m2 = -l2, l2
                do m1 = -l1, l1
                    l1l1 = ind(l1)
                    l2l2 = ind(l2)
                    m1a = abs(m1)
                    m2a = abs(m2)
                    if ( l1.eq.l2 ) then
                        k1 = l1l1 + max(m1a,m2a)
                        k12 = ind(k1) + l1l1 + min(m1a,m2a)
                    elseif (l1.gt.l2) then
                        k1 = l1l1 + m1a
                        k12 = ind(k1) + l2l2 + m2a
                    else
                        k1 = l2l2 + m2a
                        k12 = ind(k1) + l1l1 + m1a
                    endif
                    indk12(l1*(l1+1)+m1+1,l2*(l2+1)+m2+1) = k12
                end do
            end do
        end do
    end do

!    Polynomials P_k^(L,M;L',M')(0,0,1) of the shift-operators technique in the alligned-axes system

    allocate(ipntpolP(0:mxl+lexpansion,0:mxl), stat = ierr)    ! Pointers to the elements P_0^(L,0;L',0)
    if (ierr .ne. 0) call error(ierr,'Error allocating ipntpolP. Stop')

    n = (mxl+1) * (mxl+2) * (mxl+3) * (mxl+4) / 24 + lexpansion * (mxl+1) * (mxl+2) * (mxl+3) / 6

    allocate(polP(n), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Error allocating polP. Stop')

    call subpolP001(mxl+lexpansion, mxl)
    return
    end
!
!   *******************************************************************
!
  subroutine acof
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: k1, k2, k20, k200, kk, kk0, kk00, l, lp, m, m1, mp, n
    real(KREAL) :: aux, bux
    app = cero
!
!   starting elements app(00,lm)(n) = delta(l,n)
!
    k1 = 0
    do l = 0 , mxltot
        do m = 0 , l
            kk = ind(k1)
            app(l,kk) = uno
            k1 = k1 + 1
        enddo
    enddo
!
!   elements app(lm,m'm')(n)
!
    do mp = 1 , mxltot
        k2 = ind(mp) + mp
        k20 = ind(mp-1) + mp-1
        do l = mp , mxltot
            if ( l.eq.mp ) then
                m1 = mp
            else
                m1 = 0
            endif
            do m = m1 , l
                k1 = ind(l) + m
                kk = ind(k1) + k2
                kk0 = ind(k1) + k20
                do n = l-mp , l+mp , 2
                    if ( n.ge.m+mp) then
                        app(n,kk) = (2*mp-1) * ( app(n-1,kk0) * ri(n+n-1) - app(n+1,kk0) * ri(n+n+3) )
                    endif
                enddo
            enddo
        enddo
    enddo
!
!   elements app(lm,l'm')(n)
!
    do mp = 0 , mxltot
        k200 = 0
        do lp = mp+1 , mxltot
            k2 = ind(lp) + mp
            k20 = ind(lp-1) + mp
            if ( lp.gt.mp+1 ) k200 = ind(lp-2) + mp
            do l = lp , mxltot
                if ( l.eq.lp ) then
                    m1 = mp
                else
                    m1 = 0
                endif
                do m = m1 , l
                    k1 = ind(l) + m
                    kk = ind(k1) + k2
                    kk0 = ind(k1) + k20
                    kk00 = ind(k1) + k200
                    do n = l-lp , l+lp , 2
                        if ( n.ge.m+mp) then
                            aux = app(n+1,kk0) * re(n+m+mp+1) * dosl1i(n+1)
                            if ( n.gt.m+mp ) aux = aux + app(n-1,kk0) * re(n-m-mp) * dosl1i(n-1)
                            aux = aux * dosl1(lp-1)
                            if ( lp.gt.mp+1 ) aux = aux - re(lp+mp-1) * app(n,kk00)
                            app(n,kk) = aux * ri(lp-mp)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    return
    end
!
!   *******************************************************************
!
  subroutine bcof
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: k1, k2, k20, k200, kk, kk0, kk00, l, lp, m, m1, mmp, mp, n
    real(KREAL) :: aux, bux, t1, t2
    bpp = cero
!
!   starting elements bpp(lm,00)(n) = delta(l,n)
!
    k1 = 0
    do l = 0 , mxltot
        do m = 0 , l
            kk = ind(k1)
            bpp(l,kk) = uno
            k1 = k1 + 1
        enddo
    enddo
!
!   elements bpp(lm,m'm')(n)
!
    do mp = 1 , mxltot
        k2 = ind(mp) + mp
        k20 = ind(mp-1) + mp-1
        do l = mp , mxltot
            if ( l.eq.mp ) then
                m1 = mp
            else
                m1 = 0
            endif
            do m = m1 , l
                k1 = ind(l) + m
                kk = ind(k1) + k2
                kk0 = ind(k1) + k20
                do n = l-mp , l+mp , 2
                    if ( mp.gt.m ) then
                        t1 = uno
                        t2 = uno
                    else
                        t1 = -re(n-(m-mp+1)) * re(n-(m-mp+1)+1)
                        t2 = -re(n+(m-mp+1)) * re(n+(m-mp+1)+1)
                    endif
                    if ( n.ge.abs(m-mp)) then
                        if (n.eq.0) then
                            bux=cero
                        else
                            bux=t1*bpp(n-1,kk0) * dosl1i(n-1)
                        endif
                        bpp(n,kk) = dosl1(mp-1) * ( bux - t2 * bpp(n+1,kk0) * dosl1i(n+1) )
                    endif
                enddo
            enddo
        enddo
    enddo
!
!   elements bpp(lm,l'm')(n)
!
    do mp = 0 , mxltot
        k200 = 0
        do lp = mp+1 , mxltot
            k2 = ind(lp) + mp
            k20 = ind(lp-1) + mp
            if ( lp.gt.mp+1 ) k200 = ind(lp-2) + mp
            do l = lp , mxltot
                if ( l.eq.lp ) then
                    m1 = mp
                else
                    m1 = 0
                endif
                do m = m1 , l
                    k1 = ind(l) + m
                    kk = ind(k1) + k2
                    kk0 = ind(k1) + k20
                    kk00 = ind(k1) + k200
                    do n = l-lp , l+lp , 2
                        mmp = abs(m-mp)
                        if ( n.ge.mmp) then
                            aux = bpp(n+1,kk0) * re(n+mmp+1) * dosl1i(n+1)
                            if ( n.gt.mmp ) aux = aux + bpp(n-1,kk0) * re(n-mmp) * dosl1i(n-1)
                            aux = aux * dosl1(lp-1)
                            if ( lp.gt.mp+1 ) aux = aux - re(lp+mp-1) * bpp(n,kk00)
                            bpp(n,kk) = aux * ri(lp-mp)
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    return
    end
!
!   ******************************************************************
!
  subroutine emes ( m1, m2, ms, md, ss, sd )
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: m1, m1a, m2, m2a, ms, md
    real(KREAL) :: s1, s2, s12, ss, sd
    s1 = sign(1,m1)
    s2 = sign(1,m2)
    s12 = s1 * s2
    m1a = iabs(m1)
    m2a = iabs(m2)
    ms = s12 * ( m1a + m2a )
    md = s12 * iabs( m1a - m2a )
    if ( ms.eq.md ) then
        ss = uno
        sd = cero
        return
    endif
    if ( m1.lt.0 .and. m2.lt.0 ) then
        ss = -umed
    else
        ss = umed
    endif
    if ( s12.gt.cero ) then
        sd = umed
    elseif ( md.eq.0 ) then
        sd = cero
    elseif ( sign(1,m1a-m2a) .eq. s1 ) then
        sd = - umed
    else
        sd = umed
    endif
    return
    end
    
!**********************************************************************
!
!   subroutine rotar
!
!    this subroutine yields the rotation matrices rl(m',m;l) of reals spherical harmonics
!    receives the trigonometric functions of Euler angles defining the rotation
!
!**********************************************************************
  subroutine rotar(lmax, cosal, sinal, cosbet, sinbet, cosga, singa)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: l, lmax
    real(KREAL) :: cosag, cosal, cosamg, cosbet, cosga, sinag, sinal, sinamg, singa, sinbet, tgbet2
!    Initial matrices d0, r0, d1 and r1
    rl(:,:,:) = cero
    dl(:,:,:) = cero
    dl(0,0,0)  = uno
    rl(0,0,0)  = uno
    if(lmax.eq.0) return
    dl(1,1,1)  = (uno + cosbet) * umed
    dl(1,0,1)  =-sinbet/raiz2
    dl(1,-1,1) = (uno - cosbet) * umed
    dl(0,1,1)  =-dl(1,0,1)
    dl(0,0,1)  = dl(1,1,1)-dl(1,-1,1)
    dl(0,-1,1) = dl(1,0,1)
    dl(-1,1,1) = dl(1,-1,1)
    dl(-1,0,1) = dl(0,1,1)
    dl(-1,-1,1)= dl(1,1,1)
    cosag  = cosal * cosga - sinal * singa
    cosamg = cosal * cosga + sinal * singa
    sinag  = sinal * cosga + cosal * singa
    sinamg = sinal * cosga - cosal * singa
    rl(0,0,1)  = dl(0,0,1)
    rl(1,0,1)  = raiz2 * dl(0,1,1) * cosal
    rl(-1,0,1) = raiz2 * dl(0,1,1) * sinal
    rl(0,1,1)  = raiz2 * dl(1,0,1) * cosga
    rl(0,-1,1) =-raiz2 * dl(1,0,1) * singa
    rl(1,1,1)  = dl(1,1,1) * cosag - dl(1,-1,1) * cosamg
    rl(1,-1,1) =-dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
    rl(-1,1,1) = dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
    rl(-1,-1,1)= dl(1,1,1) * cosag + dl(1,-1,1) * cosamg
!    the remaining matrices are calculated using symmetry and recurrence relations by means of the subroutine dlmn.
    if ( abs(sinbet) .lt. 1.d-14 ) then
        tgbet2 = cero
    elseif ( abs(sinbet) .lt. 1.d-10 ) then
        tgbet2 = cero
        write(6,"('WARNING in ROTAR: sinbet = ', e17.10, ' takes  0')") sinbet
    else
        tgbet2 = ( uno - cosbet ) / sinbet
    endif
    do l = 2, lmax
        call dlmn(l, sinal, cosal, cosbet, tgbet2, singa, cosga)
    enddo
    return
    end
!**********************************************************************
!
!   subroutine dlmn
!
!   this subroutine generates the matrices dl(m',m;l) for a fixed value
!   of the orbital quantum number l, and it needs the dl(l-2;m',m) and 
!   dl(l-1;m',m) matrices. this subroutine uses symmetry and recurrence
!   relations. the matrices dl(m',m;l) are the rotation matrices for   
!   complex spherical harmonics
!
!**********************************************************************
  subroutine dlmn(l, sinal, cosal, cosbet, tgbet2, singa, cosga)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: iinf, isup, l, m, mp
    real(KREAL) :: al, al1, ali, aux, cosag, cosagm, cosal, cosaux, cosbet, cosga, cosmal, cosmga, cux, d1, d2
    real(KREAL) :: sgn, sinag, sinagm, sinal, singa, sinmal, sinmga, tal1, tgbet2
    iinf=1-l
    isup=-iinf
!    computation of the dl(m',m;l) matrix, mp is m' and m is m.
!    first row by recurrence: see equations 19 and 20 of reference (6)
    dl(l,l,l) = dl(isup,isup,l-1) * (uno + cosbet) * umed
    dl(l,-l,l) = dl(isup,-isup,l-1) * (uno - cosbet) * umed
    do m = isup, iinf, -1
            dl(l,m,l) = -tgbet2 * root(l+m+1) * rooti(l-m) * dl(l,m+1,l)
    enddo
!    the rows of the upper quarter triangle of the dl(m',m;l) matrix see equation 21 of reference (6)
    al = l
    al1 = al - uno
    tal1 = al + al1
    ali = uno / al1
    cosaux = cosbet * al * al1
    do mp = l-1, 0, -1
        aux = rooti(l+mp) * rooti(l-mp) * ali
        cux = root(l+mp-1) * root(l-mp-1) * al
        do m = isup, iinf, -1
            dl(mp,m,l) = aux * rooti(l+m) * rooti(l-m) * (tal1 * (cosaux - re(m) * re(mp)) * dl(mp,m,l-1) &
                    - root(l+m-1) * root(l-m-1) * cux * dl(mp,m,l-2) )
        enddo
        iinf=iinf+1
        isup=isup-1
    enddo
!    the remaining elements of the dl(m',m;l) matrix are calculated using the corresponding symmetry relations:
!        reflection ---> ((-1)**(m-m')) dl(m,m';l) = dl(m',m;l), m'<=m
!        inversion ---> ((-1)**(m-m')) dl(-m',-m;l) = dl(m',m;l)
!    reflection
    sgn = uno
    iinf = -l
    isup = l-1
    do m = l, 1, -1
        do mp = iinf, isup
            dl(mp,m,l) = sgn * dl(m,mp,l)
            sgn = -sgn
        enddo
        iinf=iinf+1
        isup=isup-1
    enddo
!    inversion
    iinf=-l
    isup=iinf
    do m = l-1, -l, -1
        sgn = -uno
        do mp = isup, iinf,- 1
            dl(mp,m,l) = sgn * dl(-mp,-m,l)
            sgn = -sgn
        enddo
        isup=isup+1
    enddo
!    computation of the rotation matrices rl(m',m;l) for real spherical harmonics using the matrices dl(m',m;l) 
!    for complex spherical harmonics: see equations 10 to 18 of reference (6)
    rl(0,0,l) = dl(0,0,l)
    cosmal = cosal
    sinmal = sinal
    sgn = - uno
    do mp = 1, l
        cosmga = cosga
        sinmga = singa
        aux = raiz2 * dl(0,mp,l)
        rl(mp,0,l) = aux * cosmal
        rl(-mp,0,l)= aux * sinmal
        do m = 1, l
            aux = raiz2 * dl(m,0,l)
            rl(0,m,l) = aux * cosmga
            rl(0,-m,l)=-aux * sinmga
            d1 = dl(-mp,-m,l)
            d2 = sgn * dl(mp,-m,l)
            cosag = cosmal * cosmga - sinmal * sinmga
            cosagm= cosmal * cosmga + sinmal * sinmga
            sinag = sinmal * cosmga + cosmal * sinmga
            sinagm= sinmal * cosmga - cosmal * sinmga
            rl(mp,m,l)  = d1 * cosag + d2 * cosagm
            rl(mp,-m,l) =-d1 * sinag + d2 * sinagm
            rl(-mp,m,l) = d1 * sinag + d2 * sinagm
            rl(-mp,-m,l)= d1 * cosag - d2 * cosagm
            aux    = cosmga * cosga - sinmga * singa
            sinmga = sinmga * cosga + cosmga * singa
            cosmga = aux
        enddo
        sgn = - sgn
        aux    = cosmal * cosal - sinmal * sinal
        sinmal = sinmal * cosal + cosmal * sinal
        cosmal = aux
    enddo
    return
    end

!! ************************************************************************
! 
!       Subroutines for computing radial factors of two center distributions of CGTO
!
  subroutine fradAgauss(i1, za)
    USE DAMZernike_Jacobi_atdens_GTO, lmxc => lexpansion, rflm => rquadscal
    implicit none

    integer(kint) :: i, i1, i2, ij, indmx, indmxk, ipoint, j, ja, jb, k1, kmax, la, lb, lm
    integer(kint) :: lmmaxc, lmmaxccab, lmxcab, ktop, m, mp, na, nb, ngsa, ngsb, nna, nnb
    real(kreal), parameter :: umbrfmax = 1.d-15
    real(kreal) :: zlmij((ldimaux+1)*(ldimaux+2)/2)
    real(kreal) :: argia, argib, argka, argkb, argsqa, argsqb, aux, aux2, bux, cosb, cux, cux2, dux, r
    real(kreal) :: rac, xxa, za

    la = ll(i1)
    xxa = xxg(i1)
    lmxcab = lmxc+la
    lmmaxc = (lmxc+1)*(lmxc+1)
    flm = 0.d0
    rac = abs(za)
    lmmaxccab = (lmxcab+1)*(lmxcab+2)/2
    call armonicosij(lmxcab, cero, za, zlmij)

    do ipoint = 1, nquadpoints
        r = rflm(ipoint)
        call flmgauss(la, rac, xxa, ngsa, lmxcab, zlmij, r)
        do m = -la, la
            do lm = 1, lmmaxc
                flm(ipoint,m,lm) =  flmmamb(lm,m)
            enddo
        enddo
    enddo
    return
    end subroutine fradAgauss
!  
! ************************************************************************
! 
  subroutine flmgauss(la, rac, xxa, nga, lmaxini, zlmij, r)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, ij, ip, itop, ja, jb, l, la, lb, ldim, lm, lmax, lmaxaux, lmaxbux, lmaxini, lmp, lmxcab, lmmax
    integer(KINT) :: m, nga, ngb
    real(KREAL) :: bsi(0:45), flmbas((ldimaux+1)*(ldimaux+2)/2), flmnnp((4+lmaxini)**2,mxl+1)
    real(KREAL) :: zlmij((ldimaux+1)*(ldimaux+2)/2)
    real(KREAL) :: a2, arg1, arg2, aux, bux, dlt, dosxxa, dosza, r, r2, ra, rac, rtop, xxa

    lmxcab = lmaxini
    ldim = 4+lmaxini
    r2 = r * r
    lmmax = (lmxcab+1)*(lmxcab+2)/2
    arg1 = 0.5d0 * xxa * (r+rac) * (r+rac)
    arg2 = 0.5d0 * xxa * (r-rac) * (r-rac)
    if (min(arg1,arg2) .gt. 100.d0) then
        flmmamb = cero
        return
    endif
    if (lmxcab .gt. 45) then
        write(6,"('Error in subroutine flmGTO, lmxcab = ', i3, ' higher than maximum allowed')")
        write(6,"('Higher value allowed due to Bessel functions parametrization = 45')")
        call error(1,'Stop')
    elseif (lmxcab .gt. 30) then
        itop = 44
        rtop = 91.d0
        call bibk91med(arg1,arg2,bsi(44),bsi(45))
    elseif (lmxcab .gt. 20) then
        itop = 29
        rtop = 61.d0
        call bibk61med(arg1,arg2,bsi(29),bsi(30))
    elseif(lmxcab .gt. 10) then
        itop = 19
        rtop = 41.d0
        call bibk41med(arg1,arg2,bsi(19),bsi(20))
    else
        itop = 9
        rtop = 21.d0
        call bibk21med(arg1,arg2,bsi(9),bsi(10))
    endif
    a2 = (arg1-arg2)*(arg1-arg2)
    do i = 1, itop
        bsi(itop-i) = (rtop-re(i+i))*bsi(itop+1-i) + a2*bsi(itop+2-i)
    enddo
    dosxxa = xxa + xxa
    aux = uno
    lm = 0
    do l = 0, lmxcab
        bux = aux * bsi(l) * re(l+l+1)
        lm = lm + 1
        flmbas(lm) = bux * zlmij(lm)
        bux = (bux + bux) * ri(l) * ri(l+1)
        do m = 1, l
            lm = lm + 1
            flmbas(lm) = bux * zlmij(lm)
            bux = bux * ri(l-m) * ri(l+m+1)
        enddo
        aux = aux * dosxxa
    enddo

    lmax = lmxcab
    if (la .eq. 0) then
        lm = 0
        lmp = 0
        do l = 0, lmax
            do m = -l, -1
                lm = lm + 1
                flmmamb(lm,0) = 0.d0
            enddo
            do m = 0, l
                lm = lm + 1
                lmp = lmp + 1
                flmmamb(lm,0) = flmbas(lmp)
            enddo
        enddo
        return
    endif
!     if .not. (la .eq. 0)  recursion on la
    lm = 0
    do l = 0, lmax
        do m = 0, l
            lm = lm + 1
            flmnnp(lm,1) = flmbas(lm)
        enddo
    enddo

!    recursion on  n  for the subsequent recursion on l
!    notice that recursion on n  runs in steps of two 2
!    whereas the storage index runs in steps of one, thus:
!    n = n0 + 2  (i-1)
    r2 = r * r
    lmaxaux = lmax
    ra = rac
    aux = r2 + ra*ra
    dosza = ra+ra
    do i = 2, la+1
        lmaxaux = lmaxaux - 1
        flmnnp(1,i) = aux * flmnnp(1,i-1)  - dosza * r2 * ri(3) * flmnnp(2,i-1)
        lm = 1
        do l = 1, lmaxaux
            do m = 0, l-1
                lm = lm + 1
                flmnnp(lm,i) = aux * flmnnp(lm,i-1)  - dosza * ( re(l-m) * ri(l+l-1) * flmnnp(ind(l-1)+m+1,i-1) &
                        + r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+m+1,i-1) )
            enddo
            lm = lm + 1
            flmnnp(lm,i) = aux * flmnnp(lm,i-1)  - dosza * r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+l+1,i-1)
        enddo
    enddo
    call recurflm(ldim, lmax, la, r2, rac, flmnnp)
    return
    end subroutine flmgauss
! 
! ***************************************************************
! 
!     Subrutine for recursion of the radial factors of the distribution
!     Indices (l,m) are contracted to a single one:
!         lm = l2 + l + m + 1
! 
  subroutine recurflm(ldim, lmaxtot, la, r2, za, fin)
    USE DAMZernike_Jacobi_atdens_GTO, fvoid => flm
    implicit none
    logical :: lcux
    integer(KINT) :: ierr, knt, l, la, lb, ldim, lg, lgg, lm, lmax, lmaxaux, lmaxbux, lmaxtot, lmpos
    integer(KINT) :: m, ma, mb, mm, n, nmax, np, npmax
    real(KREAL) :: fin(ldim*ldim,mxl+1), flmaux(ldim*ldim)
    real(KREAL), allocatable :: flm(:,:,:)
    real(KREAL) :: aux, aux2, bux, cux, r2, s1, s1m, s2, s2m, s3, s3m, s4, s4m, s5, s5m
    real(KREAL) :: s6, s6m, s7, s7m, s8, s8m, sxnn, sxnp, sxpn, sxpp, xb, za, zb, umdltm0, updltm1, umdltm1

    lmax = lmaxtot
!     nmax = max(1,la)
!     npmax = max(1,lb)
    nmax = la+1
    npmax = lb+1
    allocate(flm(ldim*ldim,la+1,-(la+1):(la+1)), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating flm. Stop')

!     loads the content of  fin  into  flm(lm,n,0)
    lmaxaux = lmax
    do n = 1, nmax
        lm = 0
        lmpos = 0
        do l = 0, lmaxaux
            do m = -l, -1
                lm = lm + 1
                flm(lm,n,0) = 0.d0
            enddo
            do m = 0, l
                lm = lm + 1
                lmpos = lmpos + 1
                flm(lm,n,0) = fin(lmpos,n)
            enddo
        enddo
        lmaxaux = lmaxaux -1
    enddo

!  =========================
!      recursion on L    
!  =========================

    if (la .ge. 1) then
        do lg = 1, la

!           recursion of factors corresponding to  n = 1

            aux = real(lg+lg-1) * 0.5d0
            aux2 = aux+aux
            lmaxaux = lmax-lg
            n = 1
!              elements with M = lg  and  M = -lg
            do l  = 0, lmaxaux
                do m  = 0, l
                    s1 = 0.d0
                    s2 = 0.d0
                    s3 = 0.d0
                    s4 = 0.d0
                    s5 = 0.d0
                    s6 = 0.d0
                    s7 = 0.d0
                    s8 = 0.d0
                    s1m = 0.d0
                    s2m = 0.d0
                    s3m = 0.d0
                    s4m = 0.d0
                    s5m = 0.d0
                    s6m = 0.d0
                    s7m = 0.d0
                    s8m = 0.d0
                    umdltm0 = 1.d0               ! umdltm0 = 1 - delta(m,0)
                    updltm1 = 1.d0               ! updltm1 = 1 + delta(m,1)
                    umdltm1 = 1.d0               ! umdltm1 = 1 - delta(m,1)
                    if (m .eq. 0) umdltm0 = 0.d0
                    if (m .eq. 1) then
                        umdltm1 = 0.d0
                        updltm1 = 2.d0
                    endif
                    if (m .gt. 0) then
                        s1 = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,lg-1)
                        s1m = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,lg-1)
                    endif
                    if (lg .gt. 1) then
                        s5 = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,-lg+1)
                        s5m = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,-lg+1)
                        s7 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,-lg+1)
                        s7m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,-lg+1)
                    endif
                    s3 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,lg-1)
                    s3m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,lg-1)
                    if (l .gt. 0) then
                        if (m .gt. 0) then
                            s2 = ri(l+l-1) * flm((l-1)*l+m,n,lg-1)
                            if (lg .gt. 1) s6m = ri(l+l-1) * flm((l-1)*l+m,n,-lg+1)
                        endif
                        if (m .gt. 1) then
                            s2m = ri(l+l-1) * flm((l-1)*l-m+2,n,lg-1)
                            if (lg .gt. 1) s6 = ri(l+l-1) * flm((l-1)*l-m+2,n,-lg+1)
                        endif
                        if (m .lt. l-1) then
                            s4 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,lg-1)
                            s4m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,lg-1)
                            if (lg .gt. 1) then
                                s8 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,-lg+1)
                                s8m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,-lg+1)
                            endif
                        endif
                    endif
                    if (m .gt. 0) then
                        flm(l*l+l-m+1,n,lg) = aux * (umdltm1*umdltm0*(-s1m+s2m) + s3m - s4m &
                                - updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                        flm(l*l+l-m+1,n,-lg) = aux *(updltm1*umdltm0*(-s1+s2) - s3 + s4 + umdltm1*umdltm0*(-s5+s6) &
                                + s7 - s8)
                    endif
                    flm(l*l+l+m+1,n,lg) = aux *(updltm1*umdltm0*(-s1+s2) + s3 - s4 - umdltm1*umdltm0*(s5-s6) &
                            - s7 + s8 )
                    flm(l*l+l+m+1,n,-lg) =aux*(umdltm1*umdltm0*(s1m-s2m) +s3m -s4m + updltm1*umdltm0*(-s5m+s6m) &
                            + s7m - s8m)
                enddo     ! End of Do on m
            enddo     ! End of Do on l

!             elements with   -lg < M < lg
            do mm = -lg+1, lg-1
                bux = ri(lg-abs(mm))
                if (lg-1-abs(mm) .le. 0) then
                    cux = 0.d0
                    lcux = .false.
                else
                    cux = bux * re(lg+abs(mm)-1)
                    lcux = .true.
                endif
                do l  = 0, lmaxaux
                    do m  = 0, l
                        s1 = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l+m+3,n,mm)
                        s1m = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l-m+3,n,mm)
                        s2 = 0.d0
                        s2m = 0.d0
                        if (m .lt. l) then
                            s2 = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l+m,n,mm)
                            s2m = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l-m,n,mm)
                        endif
                        if (m .gt. 0) then
                            flmaux(l*l+l-m+1) = bux * aux2 * (s1m+s2m - za * flm(l*l+l-m+1,n,mm) )
                            if (lcux) flmaux(l*l+l-m+1) = flmaux(l*l+l-m+1) - cux * flm(l*l+l-m+1,n+1,mm)
                        endif
                        flmaux(l*l+l+m+1) = bux * aux2 * (s1+s2 - za * flm(l*l+l+m+1,n,mm))
                        if (lcux) flmaux(l*l+l+m+1) = flmaux(l*l+l+m+1) - cux * flm(l*l+l+m+1,n+1,mm)
                    enddo     ! End of Do on m
                enddo     ! End of Do on l
                do l  = 0, lmaxaux
                    flm(l*l+l+1,n,mm) = flmaux(l*l+l+1)
                    do m  = 1, l
                        flm(l*l+l+m+1,n,mm) = flmaux(l*l+l+m+1)
                        flm(l*l+l-m+1,n,mm) = flmaux(l*l+l-m+1)
                    enddo     ! End of Do on m
                enddo     ! End of Do on l
            enddo     ! End of Do on mm

            if (lg .eq. la) exit  ! exits do on lg

!           recursion of factors with  n > 1

            if (lg .gt. 1) then
                lgg = lg
                do n  = 2, min(la,lg,la+1-lg)
                    lgg = lgg - 1
                    aux = real(lgg+lgg-1) * 0.5d0
                    aux2 = aux + aux
                    lmaxaux = lmax-lg-n+1
!                    elements with M = lgg  y  M = -lgg
                    do l  = 0, lmaxaux
                        do m  = 0, l
                            s1 = 0.d0
                            s2 = 0.d0
                            s3 = 0.d0
                            s4 = 0.d0
                            s5 = 0.d0
                            s6 = 0.d0
                            s7 = 0.d0
                            s8 = 0.d0
                            s1m = 0.d0
                            s2m = 0.d0
                            s3m = 0.d0
                            s4m = 0.d0
                            s5m = 0.d0
                            s6m = 0.d0
                            s7m = 0.d0
                            s8m = 0.d0
                            umdltm0 = 1.d0               ! umdltm0 = 1 - delta(m,0)
                            updltm1 = 1.d0               ! updltm1 = 1 + delta(m,1)
                            umdltm1 = 1.d0               ! umdltm1 = 1 - delta(m,1)
                            if (m .eq. 0) umdltm0 = 0.d0
                            if (m .eq. 1) then
                                umdltm1 = 0.d0
                                updltm1 = 2.d0
                            endif
                            if (m .gt. 0) then
                                s1 = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,lgg-1)
                                s1m = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,lgg-1)
                            endif
                            if (lgg .gt. 1) then
                                s5 = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,-lgg+1)
                                s5m = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,-lgg+1)
                                s7 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,-lgg+1)
                                s7m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,-lgg+1)
                            endif
                            s3 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,lgg-1)
                            s3m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,lgg-1)
                            if (l .gt. 0) then
                                if (m .gt. 0) then
                                    s2 = ri(l+l-1) * flm((l-1)*l+m,n,lgg-1)
                                    if (lgg .gt. 1) then
                                        s6m = ri(l+l-1) * flm((l-1)*l+m,n,-lgg+1)
                                    endif
                                endif
                                if (m .gt. 1) then
                                    s2m = ri(l+l-1) * flm((l-1)*l-m+2,n,lgg-1)
                                    if (lgg .gt. 1) then
                                        s6 = ri(l+l-1) * flm((l-1)*l-m+2,n,-lgg+1)
                                    endif
                                endif
                                if (m .lt. l-1) then
                                    s4 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,lgg-1)
                                    s4m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,lgg-1)
                                    if (lgg .gt. 1) then
                                        s8 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,-lgg+1)
                                        s8m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,-lgg+1)
                                    endif
                                endif
                            endif
                            if (m .gt. 0) then
                                flm(l*l+l-m+1,n,lgg) = aux  * (umdltm1*umdltm0*(-s1m+s2m) + s3m - s4m &
                                        - updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                                flm(l*l+l-m+1,n,-lgg) = aux  * (updltm1*umdltm0*(-s1+s2) - s3 + s4 &
                                        + umdltm1*umdltm0*(-s5+s6) + s7 - s8 )
                            endif
                            flm(l*l+l+m+1,n,lgg) =aux*(updltm1*umdltm0*(-s1+s2) + s3 - s4 &
                                    - umdltm1*umdltm0*(s5-s6) - s7 + s8)
                            flm(l*l+l+m+1,n,-lgg)=aux*(umdltm1*umdltm0*(s1m-s2m) + s3m - s4m &
                                    + updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                        enddo     ! End of Do on m
                    enddo     ! End of Do on l

!                    elements with   -lgg < M < lgg
                    do mm = -lgg+1, lgg-1
                        bux = ri(lgg-abs(mm))
                        if (lgg-1-abs(mm) .le. 0) then
                            cux = 0.d0
                            lcux = .false.
                        else
                            cux = bux * re(lgg+abs(mm)-1)
                            lcux = .true.
                        endif
                        do l  = 0, lmaxaux
                            do m  = 0, l
                                s1 = re(l+m+1) * ri(l+l+3) * r2  * flm((l+3)*l+m+3,n,mm)
                                s1m = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l-m+3,n,mm)
                                s2 = 0.d0
                                s2m = 0.d0
                                if (m .lt. l) then
                                    s2 = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l+m,n,mm)
                                    s2m = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l-m,n,mm)
                                endif
                                if (m .gt. 0) then
                                    flmaux(l*l+l-m+1) = bux * aux2 * (s1m+s2m - za * flm(l*l+l-m+1,n,mm) )
                                    if (lcux) flmaux(l*l+l-m+1) = flmaux(l*l+l-m+1)-cux * flm(l*l+l-m+1,n+1,mm)
                                endif
                                flmaux(l*l+l+m+1) = bux * aux2 * (s1+s2 - za * flm(l*l+l+m+1,n,mm) )
                                if (lcux) flmaux(l*l+l+m+1) = flmaux(l*l+l+m+1) - cux * flm(l*l+l+m+1,n+1,mm)
                            enddo     ! End of Do on m
                        enddo     ! End of Do on l
                        do l  = 0, lmaxaux
                            flm(l*l+l+1,n,mm) = flmaux(l*l+l+1)
                            do m  = 1, l
                                flm(l*l+l+m+1,n,mm) = flmaux(l*l+l+m+1)
                                flm(l*l+l-m+1,n,mm) = flmaux(l*l+l-m+1)
                            enddo     ! End of Do on m
                        enddo     ! End of Do on l
                    enddo     ! End of Do on mm
                enddo     ! End of Do on n
            endif     ! End of  if (lg .gt. 1)
        enddo     ! End of Do on lg
    endif  ! End of if (la .ge. 1)

    lmax = lmax - la
    do ma = -la, la
        knt = 0
        do l = 0, lmax
            do m = -l, l
                knt = knt + 1
                flmmamb(knt,ma) = flm(knt,1,ma)
            enddo
        enddo
    enddo
    return
    end subroutine recurflm
!
!   ***************************************************************
! 
!      Subroutine for tabulating regular spherical harmonics of (x,y=0,z), 
! 
!      zlm(l,m) = zlm(l,m,x,y=0,z)
!
!      Only stored for m .ge. 0  since the remaining ones are null because y = 0
! 
!      indices (l,m) contracted to a single one:
!          lm = l(l+1)+m      lm = 0, 1, 2, ... (lexpansion+1)2-1
!
  subroutine armonicosij(lmax, x, z, zlma)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: l, lmax, m
    real(KREAL) :: zlma((lmax+1)*(lmax+2)/2)
    real(KREAL) :: rra, rra2, x, xxa, z, zza

    xxa = x
    zza = z
    rra2 = xxa*xxa+zza*zza
    rra = sqrt(rra2)
    zlma(1) = 1.d0
    if (lmax .eq. 0) return
    zlma(2) = zza
    zlma(3) = xxa
    do l = 1, lmax-1
        zlma(ind(l+2)) = re(l+l+1) * xxa*zlma(ind(l+1))  ! element  zlma(l+1,l+1) = re(l+l+1)*(xxa*zlma(l,l)-yya*zlma(l,-l))
        zlma(ind(l+2)-1) = re(l+l+1) * zza * zlma(ind(l+1))  ! element    zlma(l+1,l)=re(l+l+1)*zza*zlma(l,l)
        do m = 0, l-1   !elements    zlma(l+1,m)=(re(l+l+1)*zza*zlma(l,m) - (l+m)*rra2*zlma(l-1,m)) * ri(l-m+1)
            zlma(ind(l+1)+m+1) = (re(l+l+1) * zza * zlma(ind(l)+m+1) - re(l+m) * rra2 * zlma(ind(l-1)+m+1) ) / re(l-m+1)
        enddo
    enddo
    return
    end subroutine armonicosij
!     
!     ---------------------------------------------------------------------------------
!
!    Polynomials P_k^(LML'M')(0,0,1) of the shift-operators technique in the axis alligned system:
!
  subroutine subpolP001(la, lb)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: ierr, k, k12, kntpol, l, la, lb, lk, llp, lm, lmp, lp, m, md, mp, ms, naux
    real(KREAL) :: aux, bux, sd, ss, suma, xa, ya, za

    kntpol = 0
    do l = 0, la
        do lp = 0, min(l, lb)
            ipntpolP(l,lp) = kntpol+1
            k12 = indk12(l*(l+1)+1,lp*(lp+1)+1)
            ss = ssv(0,0)
            aux = uno
            do k = 0, lp
                suma = cero
                do lk = k, min(l,lp)
                    suma = suma + fact(lk) * facts(l+lp-lk) * ss * app(l+lp-lk-lk, k12) * facti(lk-k) * factsi(l+lp-lk-k)
                enddo
                kntpol = kntpol+1
                polP(kntpol) = aux * suma
                aux = (aux+aux) * ri(k+1)
            enddo
            bux = dos
            do mp = 1, min(l, lb)
                k12 = indk12(l*(l+1)+mp+1,lp*(lp+1)+mp+1)
                sd = sdv(mp,mp)
                aux = bux
                do k = mp, lp
                    suma = cero
                    do lk = k, min(l,lp)
                        suma = suma + fact(lk) * facts(l+lp-lk) * sd * bpp(l+lp-lk-lk, k12) &
                                 * facti(lk-k) * factsi(l+lp-lk-k)
                    enddo
                    kntpol = kntpol+1
                    polP(kntpol) = aux * suma
                    aux = (aux+aux) * ri(k+1)
                enddo
                bux = (bux+bux) * ri(mp+1)
            enddo
        enddo
    enddo
    return
    end
!
!    -------------------------------------------------------------------------------------------------------
!
  subroutine totalchargeGTO
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: i, ia 
    real(KINT) :: charge
    charge = cero
    do ia = 1, ncen
        do i = ngini(ia), ngfin(ia)
            if (ll(i) .ne. 0) cycle
            charge = charge + coefs(nf(i)) / sqrt(xxg(i)**3)
        enddo
    enddo
    charge = charge * pi * raizpi
    write(6,"(/'Total electron charge = ', e22.15)") -charge
    return
    end
!
!    -------------------------------------------------------------------------------------------------------
!
  subroutine error(ierr, msg)
    USE DAMZernike_Jacobi_atdens_GTO
    implicit none
    integer(KINT) :: ierr
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    stop
    end
