!===============================================================================================
!                 MODULE Zernike_Jacobi_D
!===============================================================================================
    MODULE Zernike_Jacobi_MPI_D
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
        logical :: longoutput, lgbsgz, ldengz, lvalence, lzdo
        integer(KINT), parameter ::  mxcen = 500, mxl = 6, mxcap = 8000
        integer(KINT), parameter ::  mxldst = mxl+mxl
        integer(KINT), parameter :: mxk = 350, mxkextra = 100, mxl4 = 4*mxl, mxlckplm = 22, mxlckplmextra = 26
        integer(KINT), parameter :: npol = 10
        integer(KINT), parameter :: mxroot = 800, mxreal = 3000, mxfact = 150
        integer(KINT), allocatable :: ll(:), lmaxc(:), ngini(:), ngfin(:), nf(:), nzn(:), ind(:)
        integer(KINT), allocatable :: ipntap(:,:), ipntalfa(:,:)
        integer(KINT) :: jmax, kmax0, ldimaux, lexpansion, nquadpoints, kexpansion
        integer(KINT) :: lmaxbase, mxltot, mxind, nbas, ncaps, ncen
        integer(KINT) :: mxemes, mxkcof, mxlcof
        real(KREAL), parameter :: toldstorig = 1.d-12
        real(KREAL), parameter :: cero = 0.d0, uno = 1.d0, dos = 2.d0, cuatro = 4.d0, ocho = 8.d0, r16 = 16.d0, udec = 0.1d0
        real(KREAL), parameter :: umed = .5d0, pt25 = .25d0, euler = 5.7721566490153280D-1, raiz2 = 1.4142135623730950d0
        real(KREAL), parameter :: alfacutoff = 1.d-20, umbrzn = 1.d-6
        real(KREAL), allocatable  :: cfbk1(:), cfbk0(:), ckplm(:,:), ckplmextra(:,:), cffk21(:), cffk22(:), cffk23(:)
        real(KREAL), allocatable  :: akgkl(:), cfgkl1(:), cfgkl2(:), cfgkl3(:), gkl(:)
        real(KREAL), allocatable :: omeganlm(:,:), flm(:,:,:,:), flmmamb(:,:,:), rquadaux(:), rquad01(:), rquadscal(:)
        real(KREAL), allocatable :: dmat(:,:), rcen(:,:), rnor(:), xx(:), zn(:)
        real(KREAL), allocatable :: ang(:), bin(:), dl(:,:,:), rl(:,:,:), rlt(:)
        real(KREAL), allocatable :: alfasol(:), ap(:), av(:), bv(:)
        real(KREAL), allocatable :: vaux(:), weights(:), weights01(:)
        real(KREAL), allocatable :: cfa(:), cfb(:)
        real(KREAL), allocatable :: xia(:), xib(:)
        real(KREAL) :: rstar, thresmult, thresoverlap
        real(KREAL) :: pi, raizpi, pimed, pimedsqr   ! Sqrt[pi/2]
        real(KREAL) :: roblk(-mxl:mxl,-mxl:mxl)
        real(KREAL) :: fact(0:mxfact), facti(0:mxfact), re(-mxreal:mxreal), ri(-mxreal:mxreal)
        real(KREAL) :: root(0:mxroot), rooti(mxroot)
    END MODULE
!
!                 END OF MODULE Zernike_Jacobi_D
!...............................................................................................
!===============================================================================================
!                 MODULE Zernike_Jacobi_D
!===============================================================================================
    MODULE Zernike_Jacobi_GTO_MPI_D
        USE Zernike_Jacobi_MPI_D
        IMPLICIT NONE
        integer(KINT), parameter :: mxgauss = 30
        real(KREAL), allocatable :: sol(:)
        real(KREAL) :: dosl1(-mxreal:mxreal), dosl1i(-mxreal:mxreal), facts(-1:mxfact), factsi(-1:mxfact)
        logical :: lbeta
        integer(KINT), parameter :: mxprimit = 20
        integer(KINT) :: nprimitot, ncontrtot, nocalfa, nocbeta
        integer(KINT), allocatable :: isort(:), ncontr(:), nprimit(:), ipntprim(:)
        real(KREAL), allocatable :: cfcontr(:), cfcontr0(:), xxg(:), xxg0(:)
!		Coefficients and m values for the decomposition of the product of two Phi_m(phi) functions (of the real spherical harmonics)
        integer(KINT), allocatable :: msv(:,:), mdv(:,:), indk12(:,:)
        real(KREAL), allocatable :: app(:,:), bpp(:,:), ssv(:,:), sdv(:,:)
!		Polynomials P_k^(L,M;L',M') of the shift-operators technique 
        real(KREAL), allocatable :: polP(:)
!		Pointers to the elements P_0^(L,-L; L',-L')
        integer(KINT), allocatable  :: ipntpolP(:,:)
        real(KREAL), allocatable :: ftot(:,:), ftotacum(:,:), radderiv(:,:), radfunction(:,:)
    END MODULE
!
!                 END OF MODULE Zernike_Jacobi_GTO_MPI_D
!...............................................................................................
!===============================================================================================
!                 MODULE Zernike_Jacobi_STO_MPI_D
!===============================================================================================
    MODULE Zernike_Jacobi_STO_MPI_D
        USE Zernike_Jacobi_MPI_D
        IMPLICIT NONE
        logical*1, allocatable :: lsdisf(:,:)
        integer(KINT), parameter ::  mxn = 9, mxgauss = 30
        integer(KINT), allocatable :: nn(:)
        integer(KINT) :: nmaxbase, nmlmaxbase
        real(KREAL), allocatable :: sol(:), ftot(:,:), ftotacum(:,:), radderiv(:,:), radfunction(:,:)
    END MODULE
!
!                 END OF MODULE Zernike_Jacobi_STO_MPI_D
!...............................................................................................
!===============================================================================================
!                 MODULE DAMDENZERNIKEMPI320_D
!===============================================================================================
    MODULE DAMDENZERNIKEMPI320_D
        USE Zernike_Jacobi_MPI_D
        IMPLICIT NONE
        logical :: langstrom, lgradient, lgrid, lgrid2d, lindividk, lindividlk, lindividlkm, lindividl, lpoints
        integer(KINT) :: idimzlm, kmaxrep, lmaxrep, lminrep, nindices
        integer(KINT), allocatable :: indicesv(:)
        integer(KINT), parameter ::  mxrtab = 200
        real(KREAL) :: dltu, dltv, dltx, dlty, dltz, uinf, usup, vinf, vsup, xinf, xsup, yinf, ysup, zinf, zsup
        real(KREAL) :: rtab(3,mxrtab)
        real(KREAL), allocatable :: cnlm(:,:), dgkl(:), radfunction(:), radderiv(:), zlm(:), zlmdx(:), zlmdy(:), zlmdz(:)
        character(256) :: filename, fileZJname	! Names for .plt and .zernike or .jacobi files
        character(256) :: x_func_uv, y_func_uv, z_func_uv  ! Expressions of (x,y,z) in terms of (u,v) for 2D grids
    END MODULE
!
!                 END OF MODULE DAMZERNIKEMPI2017_D
!...............................................................................................
!===============================================================================================
!                 MODULE PARALELO
!===============================================================================================
    MODULE PARALELO
        USE Zernike_Jacobi_MPI_D, ONLY: KINT
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
