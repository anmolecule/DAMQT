#include <config.h>
module DAMQTinterfaceModule
   use Vartypes
   use DimensionsModule
   use PhyconModule
   use KF
   use ConstModule
   use TypescModule
   use BastypModule
   use GeometryDataModule
   use ModelDataModule, only: gModel
   use SpinsuffixModule, only: gSpinsuffix
   
!  Interface for DAMQT package
!  Builds two files (ADF.sgbs and ADF.den) which can be read by DAM2015 program for the analysis of the density.
!  File ADF.sgbs contains geometry and basis set.
!  File ADF.den contains density matrix
!  For details on DAMQT, see 
!  Rafael López, Jaime Fernández Rico, Guillermo Ramírez, Ignacio Ema, David Zorrilla
!  DAMQT 2.0: A new version of the DAMQT package for the analysis of electron density in molecules
!  http://dx.doi.org/10.1016/j.cpc.2015.02.027
! 
!  Current version: April 2015
!  Authors: Rafael lopez (rafael.lopez@uam.es) and Ignacio Ema (nacho.ema@uam.es)
! 
   integer (KINT), parameter :: mxl = 3
   integer (KINT) :: nbasnew, nbasDAM
   integer (KINT), allocatable :: ll(:), lmaxc(:), nfun(:), nn(:), ngini(:), ngfin(:), nzn(:)
   integer (KINT), allocatable :: llnew(:), nfnew(:), nnnew(:), ngininew(:), ngfinnew(:), ineglect(:)
   real(KREAL), allocatable :: ADFcharge(:), denmat(:,:), denmatxyz(:,:), rcen(:,:), xx(:), xxnew(:), zn(:)
   real(KREAL), allocatable   :: pkmat(:), pmat(:)
   real(KREAL) :: xarm((mxl+1)**2), fact(0:20)
   
   contains
   
   subroutine DAM_interface(name, lspin, lorbitals)
      
      implicit none
      character(*), intent(in) :: name
      logical, intent(in) :: lspin, lorbitals
      integer (KINT) :: i, ia, ib, ierr, ii, indnf, inuc, iset, ispin, iu, iu15, iu79, ityp, j, ja, jb, kk, n, ntyp
      integer (KINT) :: lmaxbase, nbasDAMcart, nbasis, ncaps, ncentros
      integer (KINT) :: kont, kontcapas, kontfbas, id, jd, l, ld, m, m1, m2, md, nd, nfd
      real(KREAL) :: ccc, rnor, xxd, xnorp2, xnorp11, xnorp3, xnorp21, xnorp111, xnorp4, xnorp31, xnorp22, xnorp211     
      logical :: lvalid, lbohr
      character(1000) :: units
      fact(0) = 1.0_KREAL
      do i = 1, 20
         fact(i) = fact(i-1) * i
      end do
      ncentros = 0
      ncaps = 0
      do ityp = 1, gDims%ntyp   !	Determines the number of centers and shells
         do inuc = gTypesc%nqptr(ityp), gTypesc%nqptr(ityp+1) - 1
            ncentros = ncentros + 1
            do iset = gTypesc%nbaspt(ityp), gTypesc%nbaspt(ityp+1) - 1
               ncaps = ncaps + 1
            end do
         end do
      end do
      allocate(nn(ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: nn', ierr)
      allocate(ll(ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: ll', ierr)
      allocate(nfun(ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: nf', ierr)
      allocate(xx(ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: xx', ierr)
      allocate(rcen(3,gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: rcen', ierr)
      allocate(zn(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: zn', ierr)
      allocate(lmaxc(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: lmaxc', ierr)
      allocate(ngini(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: ngini', ierr)
      allocate(ngfin(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: ngfin', ierr)
      allocate(nzn(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: nzn', ierr)
      allocate(nnnew(2*ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: nnew', ierr)
      allocate(llnew(2*ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: llnew', ierr)
      allocate(nfnew(2*ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: nfnew', ierr)
      allocate(xxnew(2*ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: xxnew', ierr)
      allocate(ngininew(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: ngininew', ierr)
      allocate(ngfinnew(gGeometrydata%natoms), stat = ierr)
      call chckmem (0, 'DAM_interface: ngfinnew', ierr)
      allocate(ineglect(2*ncaps), stat = ierr)
      call chckmem (0, 'DAM_interface: ineglect', ierr)
      n = 0
      i = 0
      ngini(1) = 1
      indnf = 1
      lmaxbase = 0
      nbasDAM = 0
      nbasDAMcart = 0
      do ityp = 1, gDims%ntyp
         do inuc = gTypesc%nqptr(ityp), gTypesc%nqptr(ityp+1) - 1
            n = n + 1
            lmaxc(n) = 0
            do iset = gTypesc%nbaspt(ityp), gTypesc%nbaspt(ityp+1) - 1
               i = i + 1
               ll(i)  = gBastyp%lqbas(iset)
               nn(i)  = gBastyp%nqbas(iset)
               xx(i)  = gBastyp%alfbas(iset)
               nfun(i) = indnf
               indnf = indnf + 2*ll(i) + 1
               nbasDAMcart = nbasDAMcart + (ll(i) + 1) * (ll(i) + 2) / 2
               if (ll(i) .gt. lmaxbase) lmaxbase = ll(i)
               if (ll(i) .gt. lmaxc(n)) lmaxc(n) = ll(i)
            end do
            ngfin(n) = i
            if (n .lt. gGeometrydata%natoms) ngini(n+1) = ngfin(n) + 1
         end do
      end do
      nbasDAM = indnf-1
      nbasis = indnf-1
      !  open section 'Geometry'
      call kfopfl (iu, 'TAPE21')
      call kfopsc (iu, 'General')
      call kfread(iu, 'Input', units)
      call caseup(units)
      i = index(units,"BOHR")
      lbohr = .false.
      if (i .ne. 0) lbohr = .true.
         
      if (lbohr) write(iuout,"('Original geometry in bohr')") 
      call kfopsc (iu, 'Geometry')   !  open section 'Geometry'
      call kfread(iu, 'ntyp', ntyp)
      call kfread(iu, 'xyz', rcen(1:3,1:gGeometrydata%natoms))
      allocate(ADFcharge(ntyp), stat = ierr)
      call chckmem (0, 'DAM_interface: ADFcharge', ierr)
      call kfread(iu, 'charge', ADFcharge)
      !  close section geometry
      call kfclsc (iu)
      call kfclfl (iu)
      i = 0
      do ityp = 1, gDims%ntyp
         do inuc = gTypesc%nqptr(ityp), gTypesc%nqptr(ityp+1) - 1
            i = i + 1
            zn(i) = ADFcharge(ityp)
         end do
      end do
      deallocate(ADFcharge, stat = ierr)
      call chckmem (1, 'DAM_interface: ADFcharge', ierr)
      do i = 1, gGeometrydata%natoms
         nzn(n) = int(zn(n)+0.1d0)
      end do
      kontcapas = 0
      kontfbas = 1
      lvalid = .true.
      dojd: do jd = 1, gGeometrydata%natoms
         ngininew(jd) = kontcapas+1
         do id = ngini(jd), ngfin(jd)
            kontcapas = kontcapas + 1
            nd  = nn(id)
            ld  = ll(id)
            xxd = xx(id)
            nfd = nfun(id)
            if (ld .le. 1) then
               nnnew(kontcapas) = nd
               llnew(kontcapas) = ld
               xxnew(kontcapas) = xxd
               nfnew(kontcapas) = kontfbas
               ineglect(kontcapas) = 0
               kontfbas = kontfbas + 2*ld+1
            else if(ld.eq.2) then
               nnnew(kontcapas) = nd
               llnew(kontcapas) = ld
               xxnew(kontcapas) = xxd
               nfnew(kontcapas) = kontfbas
               ineglect(kontcapas) = 0
               kontfbas = kontfbas + 2*ld+1
               kontcapas = kontcapas + 1
               nnnew(kontcapas) = nd
               llnew(kontcapas) = 0
               xxnew(kontcapas) = xxd
               nfnew(kontcapas) = kontfbas
               kontfbas = kontfbas + 1
               ineglect(kontcapas) = 1
            else if(ld.eq.3) then
               nnnew(kontcapas) = nd
               llnew(kontcapas) = ld
               xxnew(kontcapas) = xxd
               nfnew(kontcapas) = kontfbas
               ineglect(kontcapas) = 0
               kontfbas = kontfbas + 2*ld+1
               kontcapas = kontcapas + 1
               nnnew(kontcapas) = nd
               llnew(kontcapas) = 1
               xxnew(kontcapas) = xxd
               nfnew(kontcapas) = kontfbas
               ineglect(kontcapas) = 1
               kontfbas = kontfbas + 3
            else
               write(iuout,"('DAM interface prepared only for l .le. 3')")
               lvalid = .false.
               exit dojd
            end if
         end do
         ngfinnew(jd) = kontcapas
      end do dojd
      if (gModel%nspin .eq. 1 .and. lspin) then
         write(iuout,"(/15('-'), 3x,'DAMQT WARNING: Null spin density, loads full electron density in ADF.den', 3x,15('-')/)")
      end if
      if (lvalid) then
         nbasnew = kontfbas-1
         allocate(denmat(nbasnew,nbasnew), stat = ierr)
         call chckmem (0, 'DAM_interface: denmat', ierr)
         allocate(denmatxyz(nbasnew,nbasnew), stat = ierr)
         call chckmem (0, 'DAM_interface: denmatxyz', ierr)
         allocate (pmat(gDims%naosx), STAT = ierr)
         call chckmem (0, 'DAM_interface: pmat', ierr)
         call kfopfl (iu15, 'TAPE15')
         call kfread (iu15, 'Pmat%Pmat'//gSpinsuffix%sspin(1), pmat)
         if (gModel%nspin .eq. 2) then
            allocate (pkmat(gDims%naosx), STAT = ierr)
            call chckmem (0, 'DAM_interface: pkmat', ierr)
            call kfread (iu15, 'Pmat%Pmat'//gSpinsuffix%sspin(2), pkmat)
            if (lspin) then
               pmat = pmat - pkmat
            else
               pmat = pmat + pkmat
            end if
         end if
         call kfclfl (iu15)
         kont = 0
         do m1 = 1, nbasnew
            do m2 = 1, nbasnew
               denmat(m1,m2) = 0.0_KREAL
               if (m1.ge.m2) then
                  kont = kont+1
                  denmatxyz(m1,m2) = pmat(kont)
                  denmatxyz(m2,m1) = pmat(kont)
               end if
            end do
         end do
         
         do l = 0, mxl
            xarm(l*(l+1)+1) = (l+l+1.0_KREAL) / (2.0_KREAL * gConstants%pi)
            do m = 1, l
               xarm(l*(l+1)+m+1) = xarm(l*(l+1)+m) / ((l+m)*(l-m+1.0_KREAL))
               xarm(l*(l+1)-m+1) = xarm(l*(l+1)+m+1)
            end do
         end do
         do l = 0, mxl
            xarm(l*(l+1)+1) = xarm(l*(l+1)+1) * 0.5_KREAL
         end do
         xarm = sqrt(xarm)
!        Angular norm for Cartesians
         xnorp2 = 0.5_KREAL*sqrt(5.0_KREAL / gConstants%pi) !angular norm of x*x, y*y o z*z
         xnorp11 = 0.5_KREAL*sqrt(15.0_KREAL / gConstants%pi) !angular norm of x*y, x*z o y*z
         xnorp3 = 0.5_KREAL*sqrt(7.0_KREAL / gConstants%pi) !angular norm of x*x*x, ...
         xnorp21 = 0.5_KREAL*sqrt(35.0_KREAL / gConstants%pi) !angular norm of x*x*y, ...
         xnorp111 = 0.5_KREAL*sqrt(105.0_KREAL / gConstants%pi) !angular norm of x*y*z, ...
         xnorp4 = 1.5_KREAL*sqrt(1.0_KREAL / gConstants%pi) !angular norm of x*x*x*x, ...
         xnorp31 = 3.0_KREAL*xnorp3 !angular norm of x*x*x*y, ...
         xnorp22 = xnorp111 !angular norm of x*x*y*y, ...
         xnorp211 = 3.0_KREAL**xnorp21 !angular norm of x*x*y*z, ...
!        Transforms from Cartesian to spherical over the first index
         id = 0
         do ii = 1, ngfin(gGeometrydata%natoms)
            id = id+1
            nd  = nnnew(id)
            ld  = llnew(id)
            nfd = nfnew(id)
            if (ld.eq.0) then
               do m = 1, nbasnew
                  denmat(nfd,m) = denmatxyz(nfd,m)
               end do
            end if
            if (ld.eq.1) then ! order in adf: x,y,z. order in sphericals: y,z,x.
               do m = 1, nbasnew
                  denmat(nfd+2,m) = denmatxyz(nfd  ,m)
                  denmat(nfd  ,m) = denmatxyz(nfd+1,m)
                  denmat(nfd+1,m) = denmatxyz(nfd+2,m)
               end do
            end if
            if (ld.eq.2) then ! order in adf: x*x,x*y,x*z,y*y,y*z,z*z. order in sphericals: 3d(-2,1,0,1,2),3s.
               do m = 1, nbasnew
                  ccc = denmatxyz(nfd+0,m)*xnorp2  ! coefficient x*x (without angular norm) =
                  denmat(nfd+5,m) = denmat(nfd+5,m) + ccc/(xarm(1)*3.0_KREAL)    ! + 1/3 3s
                  denmat(nfd+2,m) = denmat(nfd+2,m) - ccc/(xarm(7)*3.0_KREAL)    ! - 1/3 3d0
                  denmat(nfd+4,m) = denmat(nfd+4,m) + ccc/(xarm(9)*6.0_KREAL)    ! + 1/6 3d2
                  ccc = denmatxyz(nfd+1,m)*xnorp11 ! coefficient x*y (without angular norm) =
                  denmat(nfd+0,m) = denmat(nfd+0,m) + ccc/(xarm(5)*6.0_KREAL)    ! + 1/6 3d-2
                  ccc = denmatxyz(nfd+2,m)*xnorp11 ! coefficient x*z (without angular norm) =
                  denmat(nfd+3,m) = denmat(nfd+3,m) + ccc/(xarm(8)*3.0_KREAL)    ! + 1/3 3d1
                  ccc = denmatxyz(nfd+3,m)*xnorp2  ! coefficient y*y (without angular norm) =
                  denmat(nfd+5,m) = denmat(nfd+5,m) + ccc/(xarm(1)*3.0_KREAL)    ! + 1/3 3s
                  denmat(nfd+2,m) = denmat(nfd+2,m) - ccc/(xarm(7)*3.0_KREAL)    ! - 1/3 3d0
                  denmat(nfd+4,m) = denmat(nfd+4,m) - ccc/(xarm(9)*6.0_KREAL)    ! - 1/6 3d2
                  ccc = denmatxyz(nfd+4,m)*xnorp11 ! coefficient y*z (without angular norm) =
                  denmat(nfd+1,m) = denmat(nfd+1,m) + ccc/(xarm(6)*3.0_KREAL)    ! + 1/3 3d-1
                  ccc = denmatxyz(nfd+5,m)*xnorp2  ! coefficient z*z (without angular norm) =
                  denmat(nfd+5,m) = denmat(nfd+5,m) + ccc/(xarm(1)*3.0_KREAL)    ! + 1/3 3s
                  denmat(nfd+2,m) = denmat(nfd+2,m) + 2.0_KREAL*ccc/(3.0_KREAL*xarm(7))! + 2/3 3d0
               end do
               id = id+1 ! skips the 3s
            end if
            if (ld.eq.3) then ! order in adf: x*x*x,x*x*y,x*x*z,x*y*y,x*y*z,x*z*z,y*y*y,y*y*z,y*z*z,z*z*z.
                              ! order in sphericals: 4f(-3,-2,1,0,1,2,3),4p(-1,0,1).
               do m = 1, nbasnew
                  ccc = denmatxyz(nfd+0,m)*xnorp3  ! coefficient x*x*x (without angular norm) =
                  denmat(nfd+9,m) = denmat(nfd+9,m) + 3.0_KREAL*ccc/(xarm(4)*5.0_KREAL)  ! + 3/5  4p1
                  denmat(nfd+4,m) = denmat(nfd+4,m) - ccc/(xarm(14)*10.0_KREAL)       ! - 1/10 4f1
                  denmat(nfd+6,m) = denmat(nfd+6,m) + ccc/(xarm(16)*60.0_KREAL)    ! + 1/60 4f3
                  ccc = denmatxyz(nfd+1,m)*xnorp21 ! coefficient x*x*y (without angular norm) =
                  denmat(nfd+7,m) = denmat(nfd+7,m) + ccc/(xarm(2)*5.0_KREAL)       ! + 1/5  4p-1
                  denmat(nfd+0,m) = denmat(nfd+0,m) + ccc/(xarm(10)*60.0_KREAL)    ! + 1/60 4f-3
                  denmat(nfd+2,m) = denmat(nfd+2,m) - ccc/(xarm(12)*30.0_KREAL)    ! - 1/30 4f-1
                  ccc = denmatxyz(nfd+2,m)*xnorp21 ! coefficient x*x*z (without angular norm) =
                  denmat(nfd+8,m) = denmat(nfd+8,m) + ccc/(xarm(3)*5.0_KREAL)       ! + 1/5  4p0
                  denmat(nfd+3,m) = denmat(nfd+3,m) - ccc/(xarm(13)*5.0_KREAL)      ! - 1/5  4f0
                  denmat(nfd+5,m) = denmat(nfd+5,m) + ccc/(xarm(15)*30.0_KREAL)    ! + 1/30 4f2
                  ccc = denmatxyz(nfd+3,m)*xnorp21 ! coefficient x*y*y (without angular norm) =
                  denmat(nfd+9,m) = denmat(nfd+9,m) + ccc/(xarm(4)*5.0_KREAL)       ! + 1/5  4p1
                  denmat(nfd+4,m) = denmat(nfd+4,m) - ccc/(xarm(14)*30.0_KREAL)    ! - 1/30 4f1
                  denmat(nfd+6,m) = denmat(nfd+6,m) - ccc/(xarm(16)*60.0_KREAL)    ! - 1/60 4f3
                  ccc = denmatxyz(nfd+4,m)*xnorp111! coefficient x*y*z (without angular norm) =
                  denmat(nfd+1,m) = denmat(nfd+1,m) + ccc/(xarm(11)*30.0_KREAL)    ! + 1/30 4f-2
                  ccc = denmatxyz(nfd+5,m)*xnorp21 ! coefficient x*z*z (without angular norm) =
                  denmat(nfd+9,m) = denmat(nfd+9,m) + ccc/(xarm(4)*5.0_KREAL)       ! + 1/5  4p1
                  denmat(nfd+4,m) = denmat(nfd+4,m) + 2.0_KREAL*ccc/(xarm(14)*15.0_KREAL) ! + 2/15 4f1
                  ccc = denmatxyz(nfd+6,m)*xnorp3  ! coefficient y*y*y (without angular norm) =
                  denmat(nfd+7,m) = denmat(nfd+7,m) + 3.0_KREAL*ccc/(xarm(2)*5.0_KREAL)  ! + 3/5  4p-1
                  denmat(nfd+0,m) = denmat(nfd+0,m) - ccc/(xarm(10)*60.0_KREAL)    ! - 1/60 4f-3
                  denmat(nfd+2,m) = denmat(nfd+2,m) - ccc/(xarm(12)*10.0_KREAL)       ! - 1/10 4f-1
                  ccc = denmatxyz(nfd+7,m)*xnorp21 ! coefficient y*y*z (without angular norm) =
                  denmat(nfd+8,m) = denmat(nfd+8,m) + ccc/(xarm(3)*5.0_KREAL)       ! + 1/5  4p0
                  denmat(nfd+3,m) = denmat(nfd+3,m) - ccc/(xarm(13)*5.0_KREAL)      ! - 1/5  4f0
                  denmat(nfd+5,m) = denmat(nfd+5,m) - ccc/(xarm(15)*30.0_KREAL)    ! - 1/30 4f2
                  ccc = denmatxyz(nfd+8,m)*xnorp21 ! coefficient y*z*z (without angular norm) =
                  denmat(nfd+7,m) = denmat(nfd+7,m) + ccc/(xarm(2)*5.0_KREAL)       ! + 1/5  4p-1
                  denmat(nfd+2,m) = denmat(nfd+2,m) + 2.0_KREAL*ccc/(xarm(12)*15.0_KREAL) ! + 2/15 4f-1
                  ccc = denmatxyz(nfd+9,m)*xnorp3  ! coefficient z*z*z (without angular norm) =
                  denmat(nfd+8,m) = denmat(nfd+8,m) + 3.0_KREAL*ccc/(xarm(3)*5.0_KREAL)  ! + 3/5  4p0
                  denmat(nfd+3,m) = denmat(nfd+3,m) + 2.0_KREAL*ccc/(xarm(13)*5.0_KREAL)  ! + 2/5  4f0
               end do
               id = id+1 ! skips the 4p
            end if
         end do
!        Transforms from Cartesian to spherical over the second index
         denmatxyz = 0.0_KREAL
         id = 0
         do ii = 1, ngfin(gGeometrydata%natoms)
            id = id+1
            nd  = nnnew(id)
            ld  = llnew(id)
            nfd = nfnew(id)
            If (ld.eq.0) then
               do m = 1, nbasnew
                  denmatxyz(m,nfd) = denmat(m,nfd)
               end do
            end if
            if (ld.eq.1) then
               do m = 1, nbasnew
                  denmatxyz(m,nfd+2) = denmat(m,nfd  )
                  denmatxyz(m,nfd  ) = denmat(m,nfd+1)
                  denmatxyz(m,nfd+1) = denmat(m,nfd+2)
               end do
            end if
            If (ld.eq.2) then
               do m = 1, nbasnew
                  ccc = denmat(m,nfd+0)*xnorp2
                  denmatxyz(m,nfd+5) = denmatxyz(m,nfd+5) + ccc/(xarm(1)*3.0_KREAL)
                  denmatxyz(m,nfd+2) = denmatxyz(m,nfd+2) - ccc/(xarm(7)*3.0_KREAL)
                  denmatxyz(m,nfd+4) = denmatxyz(m,nfd+4) + ccc/(xarm(9)*6.0_KREAL)
                  ccc = denmat(m,nfd+1)*xnorp11
                  denmatxyz(m,nfd+0) = denmatxyz(m,nfd+0) + ccc/(xarm(5)*6.0_KREAL)
                  ccc = denmat(m,nfd+2)*xnorp11
                  denmatxyz(m,nfd+3) = denmatxyz(m,nfd+3) + ccc/(xarm(8)*3.0_KREAL)
                  ccc = denmat(m,nfd+3)*xnorp2
                  denmatxyz(m,nfd+5) = denmatxyz(m,nfd+5) + ccc/(xarm(1)*3.0_KREAL)
                  denmatxyz(m,nfd+2) = denmatxyz(m,nfd+2) - ccc/(xarm(7)*3.0_KREAL)
                  denmatxyz(m,nfd+4) = denmatxyz(m,nfd+4) - ccc/(xarm(9)*6.0_KREAL)
                  ccc = denmat(m,nfd+4)*xnorp11
                  denmatxyz(m,nfd+1) = denmatxyz(m,nfd+1) + ccc/(xarm(6)*3.0_KREAL)
                  ccc = denmat(m,nfd+5)*xnorp2
                  denmatxyz(m,nfd+5) = denmatxyz(m,nfd+5) + ccc/(xarm(1)*3.0_KREAL)    
                  denmatxyz(m,nfd+2) = denmatxyz(m,nfd+2) + 2.0_KREAL*ccc/(xarm(7)*3.0_KREAL)
               end do
               id = id+1 ! skips the 3s
            end if
            if (ld.eq.3) then
               do m = 1, nbasnew
                  ccc = denmat(m,nfd+0)*xnorp3  ! coefficient x*x*x (without angular norm) =
                  denmatxyz(m,nfd+9) = denmatxyz(m,nfd+9) + 3.0_KREAL*ccc/(xarm(4)*5.0_KREAL) ! + 3/5  4p1
                  denmatxyz(m,nfd+4) = denmatxyz(m,nfd+4) - ccc/(xarm(14)*10.0_KREAL)      ! - 1/10 4f1
                  denmatxyz(m,nfd+6) = denmatxyz(m,nfd+6) + ccc/(xarm(16)*60.0_KREAL)   ! + 1/60 4f3
                  ccc = denmat(m,nfd+1)*xnorp21 ! coefficient x*x*y (without angular norm) =
                  denmatxyz(m,nfd+7) = denmatxyz(m,nfd+7) + ccc/(xarm(2)*5.0_KREAL)       ! + 1/5  4p-1
                  denmatxyz(m,nfd+0) = denmatxyz(m,nfd+0) + ccc/(xarm(10)*60.0_KREAL)    ! + 1/60 4f-3
                  denmatxyz(m,nfd+2) = denmatxyz(m,nfd+2) - ccc/(xarm(12)*30.0_KREAL)    ! - 1/30 4f-1
                  ccc = denmat(m,nfd+2)*xnorp21 ! coefficient x*x*z (without angular norm) =
                  denmatxyz(m,nfd+8) = denmatxyz(m,nfd+8) + ccc/(xarm(3)*5.0_KREAL)       ! + 1/5  4p0
                  denmatxyz(m,nfd+3) = denmatxyz(m,nfd+3) - ccc/(xarm(13)*5.0_KREAL)      ! - 1/5  4f0
                  denmatxyz(m,nfd+5) = denmatxyz(m,nfd+5) + ccc/(xarm(15)*30.0_KREAL)    ! + 1/30 4f2
                  ccc = denmat(m,nfd+3)*xnorp21 ! coefficient x*y*y (without angular norm) =
                  denmatxyz(m,nfd+9) = denmatxyz(m,nfd+9) + ccc/(xarm(4)*5.0_KREAL)       ! + 1/5  4p1
                  denmatxyz(m,nfd+4) = denmatxyz(m,nfd+4) - ccc/(xarm(14)*30.0_KREAL)    ! - 1/30 4f1
                  denmatxyz(m,nfd+6) = denmatxyz(m,nfd+6) - ccc/(xarm(16)*60.0_KREAL)    ! - 1/60 4f3
                  ccc = denmat(m,nfd+4)*xnorp111! coefficient x*y*z (without angular norm) =
                  denmatxyz(m,nfd+1) = denmatxyz(m,nfd+1) + ccc/(xarm(11)*30.0_KREAL)    ! + 1/30 4f-2
                  ccc = denmat(m,nfd+5)*xnorp21 ! coefficient x*z*z (without angular norm) =
                  denmatxyz(m,nfd+9) = denmatxyz(m,nfd+9) + ccc/(xarm(4)*5.0_KREAL)       ! + 1/5  4p1
                  denmatxyz(m,nfd+4) = denmatxyz(m,nfd+4) + 2.0_KREAL*ccc/(xarm(14)*15.0_KREAL) ! + 2/15 4f1
                  ccc = denmat(m,nfd+6)*xnorp3  ! coefficient y*y*y (without angular norm) =
                  denmatxyz(m,nfd+7) = denmatxyz(m,nfd+7) + 3.0_KREAL*ccc/(xarm(2)*5.0_KREAL)  ! + 3/5  4p-1
                  denmatxyz(m,nfd+0) = denmatxyz(m,nfd+0) - ccc/(xarm(10)*60.0_KREAL)    ! - 1/60 4f-3
                  denmatxyz(m,nfd+2) = denmatxyz(m,nfd+2) - ccc/(xarm(12)*10.0_KREAL)       ! - 1/10 4f-1
                  ccc = denmat(m,nfd+7)*xnorp21 ! coefficient y*y*z (sin la norma angular) =
                  denmatxyz(m,nfd+8) = denmatxyz(m,nfd+8) + ccc/(xarm(3)*5.0_KREAL)       ! + 1/5  4p0
                  denmatxyz(m,nfd+3) = denmatxyz(m,nfd+3) - ccc/(xarm(13)*5.0_KREAL)      ! - 1/5  4f0
                  denmatxyz(m,nfd+5) = denmatxyz(m,nfd+5) - ccc/(xarm(15)*30.0_KREAL)    ! - 1/30 4f2
                  ccc = denmat(m,nfd+8)*xnorp21 ! coefficient y*y*z (without angular norm) =
                  denmatxyz(m,nfd+7) = denmatxyz(m,nfd+7) + ccc/(xarm(2)*5.0_KREAL)       ! + 1/5  4p-1
                  denmatxyz(m,nfd+2) = denmatxyz(m,nfd+2) + 2.0_KREAL*ccc/(xarm(12)*15.0_KREAL) ! + 2/15 4f-1
                  ccc = denmat(m,nfd+9)*xnorp3  ! coefficient z*z*z (without angular norm) =
                  denmatxyz(m,nfd+8) = denmatxyz(m,nfd+8) + 3.0_KREAL*ccc/(xarm(3)*5.0_KREAL)  ! + 3/5  4p0
                  denmatxyz(m,nfd+3) = denmatxyz(m,nfd+3) + 2.0_KREAL*ccc/(xarm(13)*5.0_KREAL)  ! + 2/5  4f0
               end do
               id = id+1 ! skips the 4p
            end if
         end do
         do m1 = 1, nbasnew
            do m2 = 1, nbasnew
               denmat(m2,m1) = denmatxyz(m2,m1)
            end do
         end do
 !       Removes the redundant Cartesian
         kk = 0
         do jd = 1, gGeometrydata%natoms
            do id = ngininew(jd), ngfinnew(jd)
               nd  = nnnew(id)
               ld  = llnew(id)
               nfd = nfnew(id)
               do md = -ld,ld
                  if (ineglect(id).eq.0) then
                     kk = kk + 1
                     do i = 1, nbasnew
                        denmatxyz(kk,i) = denmat(nfd+ld+md,i)
                     end do
                  end if
               end do
            end do
         end do
         nbasDAM = kk
         kk = 0
         do jd = 1, gGeometrydata%natoms
            do id = ngininew(jd), ngfinnew(jd)
               nd  = nnnew(id)
               ld  = llnew(id)
               nfd = nfnew(id)
               do md = -ld,ld
                  if (ineglect(id).eq.0) then
                     kk = kk + 1
                     do i = 1, nbasDAM
                        denmat(i,kk) = denmatxyz(i,nfd+ld+md)
                     end do
                 end if
               end do
            end do
         end do   
         iu79 = 79
         open(iu79,file=trim(name)//".sgbs",form='unformatted', iostat=ierr)
         if (ierr .eq. 0) then
            write(iu79) 1.d-10
            write(iu79) ncaps, nbasis, ncentros
            write(iu79) ncentros
            write(iu79) nbasis, ncaps
            write(iu79) 0.0d0
            i = 0
            do ia = 1, ncentros
               write(iu79) ngini(ia), ngfin(ia)
               if (ngini(ia) .le. 0) cycle
               do kk = ngini(ia), ngfin(ia)
                  i = i + 1
                  rnor = sqrt( (2.0d0*xx(i))**(2*nn(i)+1) / fact(2*nn(i)) )
                  write(iu79) nfun(i), nn(i), ll(i), xx(i), rnor
               end do
            end do
            write(iu79) 1.0d0, 0.0d0, 0.0d0
            write(iu79) 0.0d0, 1.0d0, 0.0d0
            write(iu79) 0.0d0, 0.0d0, 1.0d0
            if (lbohr) rcen = rcen * 0.529177249_KREAL  !  Damqt requires distances in angstrom
            do ia = 1, ncentros
               write(iu79) rcen(1,ia), rcen(2,ia), rcen(3,ia), zn(ia)
            end do
            write(iu79) 'voidgr'
            write(iu79) 0
            write(iu79) .false.
            write(iu79) .false.
            write(iu79) 1, 1
            write(iu79) 0
            write(iu79) 0.0d0
            write(iu79) 0
            write(iu79) (0, i = 1, ncentros)
            write(iu79) 'XXXXX'
            write(iu79) .false., (.true., i = 1, ncaps*ncaps-1)
            close(iu79)
            open(iu79,file=trim(name)//".den",form='unformatted', iostat=ierr)
            if (ierr .eq. 0) then
               write(iu79, iostat = ierr) nbasis, ((denmat(i,j), i=1,nbasis), j=1,nbasis)
               if (ierr .ne. 0) write(iuout,"('Error writing density matrix to file ADF.den. Err code = ', i4)") ierr
            end if
         end if   
         deallocate(denmat, stat = ierr)
         call chckmem (1, 'DAM_interface: denmat', ierr)
         deallocate(denmatxyz, stat = ierr)
         call chckmem (1, 'DAM_interface: denmatxyz', ierr)
         deallocate(pmat, stat = ierr)
         call chckmem (1, 'DAM_interface: pmat', ierr)
         if (allocated(pkmat)) then
            deallocate(pkmat, stat = ierr)
            call chckmem (1, 'DAM_interface: pkmat', ierr)
         end if
      end if
      if (lorbitals) call orbitals(name)
      
      deallocate(nn, stat = ierr)
      call chckmem (1, 'DAM_interface: nn', ierr)
      deallocate(ll, stat = ierr)
      call chckmem (1, 'DAM_interface: ll', ierr)
      deallocate(nfun, stat = ierr)
      call chckmem (1, 'DAM_interface: nf', ierr)
      deallocate(xx, stat = ierr)
      call chckmem (1, 'DAM_interface: xx', ierr)
      deallocate(rcen, stat = ierr)
      call chckmem (1, 'DAM_interface: rcen', ierr)
      deallocate(zn, stat = ierr)
      call chckmem (1, 'DAM_interface: zn', ierr)
      deallocate(lmaxc, stat = ierr)
      call chckmem (1, 'DAM_interface: lmaxc', ierr)
      deallocate(ngini, stat = ierr)
      call chckmem (1, 'DAM_interface: ngini', ierr)
      deallocate(ngfin, stat = ierr)
      call chckmem (1, 'DAM_interface: ngfin', ierr)
      deallocate(nzn, stat = ierr)
      call chckmem (1, 'DAM_interface: nzn', ierr)
      deallocate(nnnew, stat = ierr)
      call chckmem (1, 'DAM_interface: nnew', ierr)
      deallocate(llnew, stat = ierr)
      call chckmem (1, 'DAM_interface: llnew', ierr)
      deallocate(nfnew, stat = ierr)
      call chckmem (1, 'DAM_interface: nfnew', ierr)
      deallocate(xxnew, stat = ierr)
      call chckmem (1, 'DAM_interface: xxnew', ierr)
      deallocate(ngininew, stat = ierr)
      call chckmem (1, 'DAM_interface: ngininew', ierr)
      deallocate(ngfinnew, stat = ierr)
      call chckmem (1, 'DAM_interface: ngfinnew', ierr)
      deallocate(ineglect, stat = ierr)
      call chckmem (1, 'DAM_interface: ineglect', ierr)
      write(iuout,"(/15('-'), 3x,'Creates files ADF.sgbs and ADF.den for DAMQT package', 3x,15('-')/)")
      if (gModel%nspin .eq. 2 .and. lspin) &
         write(iuout,"(/15('-'), 3x, 'File ', a, '.den contains spin density matrix', 3x,15('-')/)") name
   end subroutine
      
   subroutine orbitals(name)
!
!    Reads cartesian MO orbitals coefficients of ADF and transforms them for DAMQT
!
!    Order in DAMQT: for each basis shell (z,n,l), spherical harmonics ordered as m = -l, -l+1, ... l
!
      use KF
      use DimensionsModule
      use ModelDataModule
      use SpinsuffixModule
      use CharsaModule
      use SpinorbitModule
      use ParlnnucModule
      use SymptrModule
      
      use ScfdataModule
      use SymptrcModule
      use SharedArraysModule  
      implicit none
      character(*), intent(in) :: name
      integer(KINT) :: i, iaux, ibas, id, ieig, ierr, ifrc, ii, iorb, ispin, isym, itest, iu21, iu79, ivec
      integer(KINT) :: j, jd, k, kk, ld, m, md, naos, nbas, nd, nfd, nmo, nsym, norb
      integer(KINT), allocatable:: iorder(:), npart(:)
      real(KREAL), allocatable :: eigf(:), eigval(:), MOxyz(:,:), MOsph(:,:)
      real(KREAL) :: ccc, vmax, vval, xnorp2, xnorp11, xnorp3, xnorp21, xnorp111, xnorp4, xnorp31, xnorp22, xnorp211
      character(D_SCM_MAXSTRINGLENGTH) :: filename
      call kfopfl (iu21, 'TAPE21')
      call kfopsc (iu21, 'Symmetry')
      naos = gDims%naos
      call kfread(iu21, 'nsym', nsym)
      allocate (eigf(naos*naos), stat = ierr)
      call chckmem (0, 'orbitalscart: eigf', ierr)
      allocate (eigval(naos), stat = ierr)
      call chckmem (0, 'orbitalscart: eigval', ierr)
      allocate (MOxyz(naos,naos), stat = ierr)
      call chckmem (0, 'orbitalscart: MOxyz', ierr)
      allocate (MOsph(naos,naos), stat = ierr)
      call chckmem (0, 'orbitalscart: MOsph', ierr)
      allocate (npart(naos), stat = ierr)
      call chckmem (0, 'orbitalscart: npart', ierr)
      allocate (iorder(naos*naos), stat = ierr)
      call chckmem (0, 'orbitalscart: iorder', ierr)
      xnorp2 = 0.5_KREAL*sqrt(5.0_KREAL / gConstants%pi) !angular norm of x*x, y*y o z*z
      xnorp11 = 0.5_KREAL*sqrt(15.0_KREAL / gConstants%pi) !angular norm of x*y, x*z o y*z
      xnorp3 = 0.5_KREAL*sqrt(7.0_KREAL / gConstants%pi) !angular norm of x*x*x, ...
      xnorp21 = 0.5_KREAL*sqrt(35.0_KREAL / gConstants%pi) !angular norm of x*x*y, ...
      xnorp111 = 0.5_KREAL*sqrt(105.0_KREAL / gConstants%pi) !angular norm of x*y*z, ...
      xnorp4 = 1.5_KREAL*sqrt(1.0_KREAL / gConstants%pi) !angular norm of x*x*x*x, ...
      xnorp31 = 3.0_KREAL*xnorp3 !angular norm of x*x*x*y, ...
      xnorp22 = xnorp111 !angular norm of x*x*y*y, ...
      xnorp211 = 3.0_KREAL**xnorp21 !angular norm of x*x*y*z, ...
      do ispin = 1, gModel%nspin
         MOxyz = 0.0_KREAL
         MOsph = 0.0_KREAL
         ieig   = 1
         ifrc   = 1
         norb = 0
         do isym = 1, gDims%nsym
            nbas   = gSymptrc%nfcn(isym)
            call kfopsc (iu21, gCharsa%bb(isym))
            call kfrdni (iu21, 'npart', npart, nbas, 1) ! Reads the indices of the Cartesian functions
            if (gModel%ioprel .ge. 10) then
               write(iuout,"('gModel%ioprel = ', i3, ' must be less than 10 to generate files with MO')")
               return
            else
               call kfread (iu21, 'nmo'//gSpinsuffix%sspin(ispin), nmo)
               norb = norb + nmo
               call kfrdnr (iu21, 'eps'//gSpinsuffix%sspin(ispin), eigval(ifrc), nmo, 1)
               call kfrdnr (iu21, 'Eigen-Bas'//                                                         &
                     gSpinsuffix%sspin(ispin), eigf(ieig), nbas*nmo, 1)
!              ---------------------------------------------------
!              try to standardize the phase to make nicer pictures
!              ---------------------------------------------------
               ivec = ieig
               do iorb = 1, nmo
                   vmax = 0.0_KREAL
                   itest = 0
                   do ibas = 0, nbas-1
                       vval = abs(eigf(ivec+ibas))
                       if (vval>vmax) then
                           itest = ibas
                           vmax = vval
                           if (vmax > 1.0e-10) exit
                       end if
                   end do
                   if (eigf(ivec+itest) < 0.0) then
                       eigf(ivec:ivec+nbas-1) = -eigf(ivec:ivec+nbas-1)
                   end if
                   do i = 1, nbas  ! Loads coefficients of Cartesian MO in the calculation Cartesian basis set
                      MOxyz(npart(i),ifrc+iorb-1) = eigf(ivec+i-1)
                   end do
                   ivec = ivec + nbas
               end do
!
               ifrc = ifrc + nmo
               ieig = ieig + nbas*nmo
             
            end if
         end do       
!        Sorts the MO in ascending MO energies
         do i = 1, norb
            iorder(i) = i
         end do
         do i = 2, norb
            do j = 1, i-1
               if (eigval(iorder(j)) .le. eigval(i)) cycle
               iaux = i
               do k = i, j+1, -1
                  iorder(k) = iorder(k-1)
               end do
               iorder(j) = iaux
               exit
            end do
         end do
         id = 0
         do ii = 1, ngfin(gGeometrydata%natoms)
            id = id+1
            nd  = nnnew(id)
            ld  = llnew(id)
            nfd = nfnew(id)
            if (ld.eq.0) then
               do m = 1, nbasnew
                  MOsph(nfd,m) = MOxyz(nfd,iorder(m))
               end do
            end if
            if (ld.eq.1) then ! order in adf: x,y,z. order in sphericals: y,z,x.
               do m = 1, nbasnew
                  MOsph(nfd+2,m) = MOxyz(nfd  ,iorder(m))
                  MOsph(nfd  ,m) = MOxyz(nfd+1,iorder(m))
                  MOsph(nfd+1,m) = MOxyz(nfd+2,iorder(m))
               end do
            end if
            if (ld.eq.2) then ! order in adf: x*x,x*y,x*z,y*y,y*z,z*z. order in sphericals: 3d(-2,1,0,1,2),3s.
               do m = 1, nbasnew
                  ccc = MOxyz(nfd+0,iorder(m))*xnorp2  ! coefficient x*x (without angular norm) =
                  MOsph(nfd+5,m) = MOsph(nfd+5,m) + ccc/(xarm(1)*3.0_KREAL)    ! + 1/3 3s
                  MOsph(nfd+2,m) = MOsph(nfd+2,m) - ccc/(xarm(7)*3.0_KREAL)    ! - 1/3 3d0
                  MOsph(nfd+4,m) = MOsph(nfd+4,m) + ccc/(xarm(9)*6.0_KREAL)    ! + 1/6 3d2
                  ccc = MOxyz(nfd+1,iorder(m))*xnorp11 ! coefficient x*y (without angular norm) =
                  MOsph(nfd+0,m) = MOsph(nfd+0,m) + ccc/(xarm(5)*6.0_KREAL)    ! + 1/6 3d-2
                  ccc = MOxyz(nfd+2,iorder(m))*xnorp11 ! coefficient x*z (without angular norm) =
                  MOsph(nfd+3,m) = MOsph(nfd+3,m) + ccc/(xarm(8)*3.0_KREAL)    ! + 1/3 3d1
                  ccc = MOxyz(nfd+3,iorder(m))*xnorp2  ! coefficient y*y (without angular norm) =
                  MOsph(nfd+5,m) = MOsph(nfd+5,m) + ccc/(xarm(1)*3.0_KREAL)    ! + 1/3 3s
                  MOsph(nfd+2,m) = MOsph(nfd+2,m) - ccc/(xarm(7)*3.0_KREAL)    ! - 1/3 3d0
                  MOsph(nfd+4,m) = MOsph(nfd+4,m) - ccc/(xarm(9)*6.0_KREAL)    ! - 1/6 3d2
                  ccc = MOxyz(nfd+4,iorder(m))*xnorp11 ! coefficient y*z (without angular norm) =
                  MOsph(nfd+1,m) = MOsph(nfd+1,m) + ccc/(xarm(6)*3.0_KREAL)    ! + 1/3 3d-1
                  ccc = MOxyz(nfd+5,iorder(m))*xnorp2  ! coefficient z*z (without angular norm) =
                  MOsph(nfd+5,m) = MOsph(nfd+5,m) + ccc/(xarm(1)*3.0_KREAL)    ! + 1/3 3s
                  MOsph(nfd+2,m) = MOsph(nfd+2,m) + 2.0_KREAL*ccc/(3.0_KREAL*xarm(7))! + 2/3 3d0
               end do
               id = id+1 ! skips the 3s
            end if
            if (ld.eq.3) then ! order in adf: x*x*x,x*x*y,x*x*z,x*y*y,x*y*z,x*z*z,y*y*y,y*y*z,y*z*z,z*z*z.
                              ! order in sphericals: 4f(-3,-2,1,0,1,2,3),4p(-1,0,1).
               do m = 1, nbasnew
                  ccc = MOxyz(nfd+0,iorder(m))*xnorp3  ! coefficient x*x*x (without angular norm) =
                  MOsph(nfd+9,m) = MOsph(nfd+9,m) + 3.0_KREAL*ccc/(xarm(4)*5.0_KREAL)  ! + 3/5  4p1
                  MOsph(nfd+4,m) = MOsph(nfd+4,m) - ccc/(xarm(14)*10.0_KREAL)       ! - 1/10 4f1
                  MOsph(nfd+6,m) = MOsph(nfd+6,m) + ccc/(xarm(16)*60.0_KREAL)    ! + 1/60 4f3
                  ccc = MOxyz(nfd+1,iorder(m))*xnorp21 ! coefficient x*x*y (without angular norm) =
                  MOsph(nfd+7,m) = MOsph(nfd+7,m) + ccc/(xarm(2)*5.0_KREAL)       ! + 1/5  4p-1
                  MOsph(nfd+0,m) = MOsph(nfd+0,m) + ccc/(xarm(10)*60.0_KREAL)    ! + 1/60 4f-3
                  MOsph(nfd+2,m) = MOsph(nfd+2,m) - ccc/(xarm(12)*30.0_KREAL)    ! - 1/30 4f-1
                  ccc = MOxyz(nfd+2,iorder(m))*xnorp21 ! coefficient x*x*z (without angular norm) =
                  MOsph(nfd+8,m) = MOsph(nfd+8,m) + ccc/(xarm(3)*5.0_KREAL)       ! + 1/5  4p0
                  MOsph(nfd+3,m) = MOsph(nfd+3,m) - ccc/(xarm(13)*5.0_KREAL)      ! - 1/5  4f0
                  MOsph(nfd+5,m) = MOsph(nfd+5,m) + ccc/(xarm(15)*30.0_KREAL)    ! + 1/30 4f2
                  ccc = MOxyz(nfd+3,iorder(m))*xnorp21 ! coefficient x*y*y (without angular norm) =
                  MOsph(nfd+9,m) = MOsph(nfd+9,m) + ccc/(xarm(4)*5.0_KREAL)       ! + 1/5  4p1
                  MOsph(nfd+4,m) = MOsph(nfd+4,m) - ccc/(xarm(14)*30.0_KREAL)    ! - 1/30 4f1
                  MOsph(nfd+6,m) = MOsph(nfd+6,m) - ccc/(xarm(16)*60.0_KREAL)    ! - 1/60 4f3
                  ccc = MOxyz(nfd+4,iorder(m))*xnorp111! coefficient x*y*z (without angular norm) =
                  MOsph(nfd+1,m) = MOsph(nfd+1,m) + ccc/(xarm(11)*30.0_KREAL)    ! + 1/30 4f-2
                  ccc = MOxyz(nfd+5,iorder(m))*xnorp21 ! coefficient x*z*z (without angular norm) =
                  MOsph(nfd+9,m) = MOsph(nfd+9,m) + ccc/(xarm(4)*5.0_KREAL)       ! + 1/5  4p1
                  MOsph(nfd+4,m) = MOsph(nfd+4,m) + 2.0_KREAL*ccc/(xarm(14)*15.0_KREAL) ! + 2/15 4f1
                  ccc = MOxyz(nfd+6,iorder(m))*xnorp3  ! coefficient y*y*y (without angular norm) =
                  MOsph(nfd+7,m) = MOsph(nfd+7,m) + 3.0_KREAL*ccc/(xarm(2)*5.0_KREAL)  ! + 3/5  4p-1
                  MOsph(nfd+0,m) = MOsph(nfd+0,m) - ccc/(xarm(10)*60.0_KREAL)    ! - 1/60 4f-3
                  MOsph(nfd+2,m) = MOsph(nfd+2,m) - ccc/(xarm(12)*10.0_KREAL)       ! - 1/10 4f-1
                  ccc = MOxyz(nfd+7,iorder(m))*xnorp21 ! coefficient y*y*z (without angular norm) =
                  MOsph(nfd+8,m) = MOsph(nfd+8,m) + ccc/(xarm(3)*5.0_KREAL)       ! + 1/5  4p0
                  MOsph(nfd+3,m) = MOsph(nfd+3,m) - ccc/(xarm(13)*5.0_KREAL)      ! - 1/5  4f0
                  MOsph(nfd+5,m) = MOsph(nfd+5,m) - ccc/(xarm(15)*30.0_KREAL)    ! - 1/30 4f2
                  ccc = MOxyz(nfd+8,iorder(m))*xnorp21 ! coefficient y*z*z (without angular norm) =
                  MOsph(nfd+7,m) = MOsph(nfd+7,m) + ccc/(xarm(2)*5.0_KREAL)       ! + 1/5  4p-1
                  MOsph(nfd+2,m) = MOsph(nfd+2,m) + 2.0_KREAL*ccc/(xarm(12)*15.0_KREAL) ! + 2/15 4f-1
                  ccc = MOxyz(nfd+9,iorder(m))*xnorp3  ! coefficient z*z*z (without angular norm) =
                  MOsph(nfd+8,m) = MOsph(nfd+8,m) + 3.0_KREAL*ccc/(xarm(3)*5.0_KREAL)  ! + 3/5  4p0
                  MOsph(nfd+3,m) = MOsph(nfd+3,m) + 2.0_KREAL*ccc/(xarm(13)*5.0_KREAL)  ! + 2/5  4f0
               end do
               id = id+1 ! skips the 4p
            end if
         end do
!        Removes the redundant Cartesian
         MOxyz = MOsph
         kk = 0
         do jd = 1, gGeometrydata%natoms
            do id = ngininew(jd), ngfinnew(jd)
               nd  = nnnew(id)
               ld  = llnew(id)
               nfd = nfnew(id)
               do md = -ld,ld
                  if (ineglect(id).eq.0) then
                     kk = kk + 1
                     do i = 1, nbasnew
                        MOsph(kk,i) = MOxyz(nfd+ld+md,i)
                     end do
                  end if
               end do
            end do
         end do 
         iu79 = 79
         if (ispin .eq. 1) then
            write(filename,"(a)") trim(name)//".SLorba"
         else
            write(filename,"(a)") trim(name)//".SLorbb"
         endif
         open(iu79,file=filename,form='unformatted', iostat=ierr)
         if (ierr .eq. 0) then
            write(iu79) norb, ((MOsph(i,j),i=1,norb),j=1,norb)
            close(iu79)
            write(iuout,"(/15('-'), 3x,'Creates file ', a, ' with orbitals for DAMQT package', 3x,15('-')/)") trim(filename)
         else
            write(iuout,"('Error writing file with orbitals')")
         end if
      end do
      deallocate(eigf, eigval, MOxyz, MOsph, iorder, npart, stat = ierr)
      call chckmem (1, 'orbitalscart: eigf', ierr)
      call chckmem (1, 'orbitalscart: eigval', ierr)     
      call chckmem (1, 'orbitalscart: MOxyz', ierr)
      call chckmem (1, 'orbitalscart: MOsph', ierr)
      call chckmem (1, 'orbitalscart: iorder', ierr)
      call chckmem (1, 'orbitalscart: npart', ierr)
   end subroutine
end module
