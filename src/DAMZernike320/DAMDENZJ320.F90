!  Copyright 2011-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
! Program for the tabulation of molecular density from its expansion in Canterakis-Zernike or Jacobi expansion  
!
! Requires file DAM320_GLOBAL.F90 which must be compiled first (contains the modules)
!
! Files with Canterakis-Zernike (.zernike) or Jacobi (.jacobi) expansions must be available. These files
!       can be generated with the programs: DAMZernike-Jacobi_GTO, DAMZernike-Jacobi_STO,
!       DAMZernike-Jacobi_GTO_mpi or DAMZernike-Jacobi_STO_mpi
!
! Version of September 2018
!
! #define DBLPRCGRID    ! Uncomment this line  if double precision grid is wanted
  program DAMDENZernike310
    USE DAMDENZERNIKE320_D
    USE strings, only: lowercase
    implicit none
    logical :: existe, lindivid
    integer(KINT) :: i, ierr, ipoint, k, kmax, knt, ktop, l, ltop, m, nu, numrtab
    integer(KINT), allocatable :: indices(:)
    real(KREAL) :: aux, bux, den, dxden, dyden, dzden, x, y, z
    real(KREAL4) :: tarray(2), tiempo, dtime
    namelist / options / dltu, dltv, dltx, dlty, dltz, filename, fileZJname, indices, iswindows, kmaxrep, lechelon, lgradient, &
            lgrid, lgrid2d, lindividk, lindividlk, lindividlkm, lindividl, ljacobi, lmaxrep, lminrep, lpoints, &
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
    langstrom = .true.        ! If .false. distances in bohr
    lechelon = .false.        ! if true number of functions per l equal to max(lexpansion+1,nexpansion)-l
    lindividk = .false.     ! If true projection functions with given values of k index and all l and m compatible indices are used
    lindividl = .false.     ! If true projection functions with given values of l index and all k and m compatible indices are used
    lindividlk = .false.     ! If true projection functions with given values of k, l indices and all m compatible indices are used
    lindividlkm = .false. ! If true only projection functions with given values of k, l, m indices are used
    ljacobi = .false.        ! If .true. expansion in Jacobi functions
    lgrid = .true.            ! If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    lgradient = .false.        ! If .true. gradient components of the density computed and, if lgrid = .true., tabulated in files
                                            !    projectname-d-dx.pltd, projectname-d-dy.pltd, projectname-d-dz.pltd
    lgrid2d = .false.        ! If .true. computes a 2D grid. (x,y,z) are given in terms of (u,v)
    lpoints = .false.        ! If .true. computes in selected points and prints out the results in the standard output.
                                            ! Points must be given in cartesian coordinates.
                                            ! If lgrid .eq. .true., these coordinates must be placed after the grid data
    kmaxrep = 10            ! Highest k in expansion
    lminrep = 0            ! Lowest l in expansion
    lmaxrep = 10            ! Highest l in expansion
    numrtab = 0            ! Number of tabulation points supplied in namelist
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
write(6,*) 'ljacobi = ', ljacobi
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
!    Subroutine cabecera: writes header for .plt files (binary)
!
  subroutine cabecera(nx, ny, nz, x1, x2, y1, y2, z1, z2, iuni, s)
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
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
    USE DAMDENZERNIKE320_D
    implicit none
    integer(KINT) :: ierr
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    stop
    end
