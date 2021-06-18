!  Copyright 2011-2019, Rafael Lopez
! 
!  This is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
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
! Program for computation of Zernike 3D and Jacobi moments of a molecular density of GTO
!
! Version of September 2018
! 
!       Subroutines for computing radial factors of two center distributions
   subroutine fradAB(na, la, xxa, nb, lb, xxb, za, xb, zb)
    USE Zernike_Jacobi_STO_MPI_D, lmxc => lexpansion, rflm => rquadscal
    implicit none

    logical napar, nbpar, lstng
    logical lneglig(mxgauss*mxgauss)
    integer(KINT), parameter :: k1conv = 17
    integer(KINT) :: i, ij, imxjmax, imxjmax1, indmx, indmxk, ipoint, j, ja, jb, jmax1, k1, kmax, la, lb, lm
    integer(KINT) :: lmmaxc, lmmaxccab, lmxcab, ktop, m, mp, na, nb, ngsa, ngsb, nna, nnb
    real(KREAL), parameter :: umbrfmax = 1.d-15
    real(KREAL) :: bia(-1:mxk), bia0(-1:mxk), bib(-1:mxk), bib0(-1:mxk), bka(-1:mxk), bka0(-1:mxk), bkb(-1:mxk), bkb0(-1:mxk)
    real(KREAL) :: bma(-3:mxk), bmb(-3:mxk), cijv(mxgauss*mxgauss)
    real(KREAL) :: plm(-1:mxk,0:ldimaux), zlm((ldimaux+1)*(ldimaux+2)/2), zlmij((ldimaux+1)*(ldimaux+2)/2,mxgauss*mxgauss)
    real(KREAL) :: argia, argib, argka, argkb, argsqa, argsqb, aux, aux2, bux, cosB, cux, cux2, dux, r, rabsq
    real(KREAL) :: rac, raizrac, raizrbc, rbc, rinf, rmax, rmin, sinB, xb, xxa, xxabni, xxasq, xxb, xxbsq, xiij, Xij, za, zb, Zij

    nna = na - la
    nnb = nb - lb

    if (nna-2*(nna/2) .eq. 0) then
        napar = .true.
    else
        napar = .false.
    endif

    if (nnb-2*(nnb/2) .eq. 0) then
        nbpar = .true.
    else
        nbpar = .false.
    endif

    flm = 0.d0

    lmxcab = lmxc+la+lb+(nna-1)/2+(nnb-1)/2
    lmmaxc = (lmxc+1)*(lmxc+1)
    rac = abs(za)
    rbc = sqrt(xb*xb+zb*zb)
    raizrac = sqrt(rac)
    raizrbc = sqrt(rbc)

!     Legendre functions

    if (rbc .gt. 1.d-20) then
        cosB = zb / rbc
        sinB = sqrt(1.d0 - cosB*cosB)
    else
        cosB = 1.d0
        sinB = 0.d0
    endif
    do m = 0, lmxcab
        do i = 0, m-2
            plm(i,m) = 0.d0
        enddo
        plm(m-1,m) = 0.d0
        plm(m,m) = (-0.5d0*sinB)**m *fact(m+m) / fact(m)
        plm(m+1,m) = re(m+m+1) * cosB * plm(m,m)
        do i = m+1, max(51,mxk-1)
            plm(i+1,m) = (re(i+i+1) * cosB * plm(i,m) - re(i+m) * plm(i-1,m) ) / re(i-m+1)
        enddo
    enddo

    rmax = max(rac,rbc)
    rmin = min(rac,rbc)
    lstng = .false.
    if (rmin/rmax .gt. .9d0) then  ! when rac and rbc have close values, the series converges slowly
        lstng = .true.              ! in this case it is preferable tu use STO-nG expansions of chi_A and chi_B
        rmin = .45d0 * (rmin+rmax)    ! rmin = .9 * (rmin+rmax) / 2
        rmax = .6d0 * (rmin+rmax)     ! rmax = 1.2 * (rmin+rmax) / 2

!      Reads data of the highest STO-nG expansions available for
!      functions at A and B

        ngsa = 30
        call STONG5(Na,La,NGsa,XiA,Cfa)
        ngsb = 30
        call STONG5(Nb,Lb,NGsb,Xib,Cfb)
!      Computes the regular spherical harmonics of theta_ij
        ij = 0
        lmmaxccab = (lmxcab+1)*(lmxcab+2)/2
        xxasq = xxa * xxa
        xxbsq = xxb * xxb
        xxabni = 1.d0 / (xxa**(nna-1) * xxb**(nnb-1))
        rabsq = xb*xb + (zb-za) * (zb-za)
        do ja = 1, ngsa
            do jb = 1, ngsb
                ij = ij + 1
                xiij = 1.d0 / (xiA(ja) * xxasq + xiB(jb) * xxbsq)
                aux = xiA(ja) * xiB(jb) * xxasq * xxbsq * rabsq * xiij
                if (aux .lt. 50.d0) then
                    cijv(ij) = cfA(ja) * cfB(jb) * xxabni * exp(-aux)
                    Xij = xb * xiB(jb) * xxbsq * xiij
                    Zij = za + (zb - za) * xiB(jb) * xxbsq * xiij
                    call armonicosij(lmxcab, Xij, Zij, zlm)
                    do lm = 1, lmmaxccab
                        if (abs(zlm(lm)) .gt. 1.d-100) then
                            zlmij(lm,ij) = zlm(lm)
                        else
                            zlmij(lm,ij) = 0.d0
                        endif
                    enddo
                    lneglig(ij) = .false.
                else
                    lneglig(ij) = .true.
                endif
            enddo
        enddo
    endif

!      Computes modified Bessel I functions 
!            bia(l+1/2,xxa*rac)   and    bib(l+1/2,xxb*rcb)

    bia0 = 0.d0
    bib0 = 0.d0
    bka0 = 0.d0
    bkb0 = 0.d0
    if (rac .lt. .4d0) then
        ktop = min(21,kmax0) + 2*lmxcab + 1
    else
        ktop = kmax0 + lmxcab + 1
    endif
    indmx = ktop
    argia = xxa * rac
    aux = 1.d0
    aux2 = 1.d0
    bux = .25d0 * argia * argia
    jmax1 = jmax + 1
    imxjmax = indmx+jmax
    imxjmax1 = indmx+jmax + 1
    do j = 1, jmax
        aux = 1.d0 + aux * bux * 2.d0 / (re(jmax1-j) * re(2*(imxjmax1-j)+1))
        aux2 = 1.d0 + aux2 * bux * 2.d0 / (re(jmax1-j) * re(2*(imxjmax-j)+1))
    enddo

    bux = exp(-argia)
    bia0(indmx) = aux * bux
    bia0(indmx-1) = aux2 * bux

    argsqa = argia * argia
    do i = indmx-1, 0, -1
        bia0(i-1) = argsqa * cfbk0(i) * bia0(i+1) + bia0(i)
    enddo

    if (rbc .lt. .4d0) then
        ktop = min(21,kmax0) + 2*lmxcab + 1
    else
        ktop = min(ktop + lmxcab,kmax0 + lmxcab + 1)
    endif

    indmx = ktop
    argib = xxb * rbc
    cux = 1.d0
    cux2 = 1.d0
    dux = .25d0 * argib * argib
    jmax1 = jmax + 1
    imxjmax = indmx+jmax
    imxjmax1 = indmx+jmax + 1
    do j = 1, jmax
        cux = 1.d0 + cux * dux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax1-j)+1) )
        cux2 = 1.d0 + cux2 * dux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax-j)+1) )
    enddo

    bux = exp(-argib)
    bib0(indmx) = cux * bux
    bib0(indmx-1) = cux2 * bux
    argsqb = argib * argib

    do i = indmx-1, 0, -1
        bib0(i-1) = argsqb * cfbk0(i) * bib0(i+1) + bib0(i)
    enddo

!      Computes modified Bessel K functions
!            bka(l+1/2,xxa*rac)   y    bkb(l+1/2,xxb*rcb)

    indmxk = ktop
    argka = xxa * rac
    bka0(-1) = 1.d0 / argka
    bka0(0) = 1.d0
    argkb = xxb * rbc
    bkb0(-1) = 1.d0 / argkb
    bkb0(0) = 1.d0
    argsqa = argka * argka
    argsqb = argkb * argkb

    bma(-3) = 0.d0
    bma(-2) = 0.d0
    bmb(-3) = 0.d0
    bmb(-2) = 0.d0
    do i = 0, indmxk-1
        bka0(i+1) = argsqa * cfbk0(i) * bka0(i-1) + cfbk1(i) * bka0(i)
        bkb0(i+1) = argsqb * cfbk0(i) * bkb0(i-1) + cfbk1(i) * bkb0(i)
    enddo
    do ipoint = 1, nquadpoints
        r = rflm(ipoint)
        if (r .le. rmin) then  ! ( 0 < r < min(rac,rbc)
            k1 = k1conv * log(10.d0)/(log(rmin/r)+log(rmax/r))
            kmax = min(indmx-lmxcab-1, kmax0, k1) + lmxcab
            if (rmin .lt. .4d0) then
                kmax = min(21 + lmxcab,kmax) + 1
            endif
            ktop = kmax + lmxcab + 1

            call subbmk1(xxa, xxb, ktop, rac, rbc, r, bka0, bia, bkb0, bib, bma, bmb)
            call subflm(nna, napar, la, xxa, Rac, nnb, nbpar, lb, xxb, Rbc , xb, zb, plm , bma, bmb, lmxcab, ktop, kmax, r)

            do mp = -lb, lb
                do m = -la, la
                    do lm = 1, lmmaxc
                        flm(ipoint,m,mp,lm) =  flmmamb(lm,m,mp)
                    enddo
                enddo
            enddo
        else if(r .gt. rmin .and. r .le. rmax) then ! ( min(rac,rbc) < r <= max(rac,rbc) )
            if (lstng) then ! If A and B distances to C are close to each other, uses STO-nG expansions
                call flmGTO(la, xxa, Rac, lb, xxb, xb, zb, Rbc, ngsa , ngsb, lmxcab, lneglig, zlmij, cijv, r)
                do mp = -lb, lb
                    do m = -la, la
                        do lm = 1, lmmaxc
                            flm(ipoint,m,mp,lm) =  flmmamb(lm,m,mp)
                        enddo
                    enddo
                enddo
            else if (rac .lt. rbc) then
                    k1 = k1conv*log(10.d0)/(log(r/rmin)+log(rmax/r))
                    kmax = min(kmax0, k1) + lmxcab
                    if (rmin .lt. .4d0) then
                        kmax = min(31 + lmxcab,kmax) + 1
                    endif
                    ktop = kmax + lmxcab + 1
                    call subbmk2(xxa, xxb, ktop, rac, rbc, r, bka, bia0, bkb0, bib, bma, bmb)
                    call subflm(nna, napar, la, xxa, Rac, nnb, nbpar, lb, xxb , Rbc , xb, zb, plm , bma, bmb, lmxcab, ktop, kmax, r)
                    do mp = -lb, lb
                        do m = -la, la
                            do lm = 1, lmmaxc
                                flm(ipoint,m,mp,lm) =  flmmamb(lm,m,mp)
                            enddo
                        enddo
                    enddo
            else
                k1 = k1conv*log(10.d0)/(log(r/rmin)+log(rmax/r))
                kmax = min(kmax0, k1) + lmxcab
                if (rmin .lt. .4d0) then
                    kmax = min(31 + lmxcab,kmax) + 1
                endif
                ktop = kmax + lmxcab + 1
                call subbmk2(xxb, xxa, ktop, rbc, rac, r, bkb, bib0, bka0, bia, bmb, bma)
                call subflm(nna, napar, la, xxa, Rac, nnb, nbpar, lb, xxb , Rbc , xb, zb, plm , bma, bmb, lmxcab, ktop, kmax, r)
                do mp = -lb, lb
                    do m = -la, la
                        do lm = 1, lmmaxc
                            flm(ipoint,m,mp,lm) =  flmmamb(lm,m,mp)
                        enddo
                    enddo
                enddo
            endif
        else ! ( max(rac,rbc) < r)
            k1 = k1conv*log(10.d0)/(log(r/rmin)+log(r/rmax))
            kmax = min(kmax0, k1) + lmxcab
            if (rmax .lt. .4d0) then
                kmax = min(31 + lmxcab,kmax) + 1
            endif
            ktop = kmax + lmxcab + 1
            call subbmk3(xxa, xxb, ktop, rac, rbc, r, bka, bia0, bkb, bib0, bma, bmb)
            call subflm(nna, napar, la, xxa, Rac, nnb, nbpar, lb, xxb, Rbc , xb, zb, plm , bma, bmb, lmxcab, ktop, kmax , r)
            do mp = -lb, lb
                do m = -la, la
                    do lm = 1, lmmaxc
                        flm(ipoint,m,mp,lm) =  flmmamb(lm,m,mp)
                    enddo
                enddo
            enddo
        endif
    enddo
    return
    end subroutine fradAB
! 
! ************************************************************************
! 
!     subroutine for computing Bessel I functions and M funcions for translating STO from A and B to C in case of
!           0 <= r < min(Ra,Rb)
! 
  subroutine subbmk1(xxa, xxb, ktop, rac, rbc, r, bka, bia, bkb, bib, bma, bmb)
    USE Zernike_Jacobi_STO_MPI_D
    implicit none

    integer(KINT) :: i, imxjmax, imxjmax1, indmx, j, jmax1, k, ktop
    real(KREAL) :: argia, argib, argsqa, argsqb, aux, aux2, auxa, auxb, bux, buxa, buxb, cux, cux2
    real(KREAL) :: dux, r, rac, rbc, xxa, xxb
    real(KREAL) :: bia(-1:mxk), bib(-1:mxk), bka(-1:mxk), bkb(-1:mxk), bma(-3:mxk), bmb(-3:mxk)

    if (r .le. 0) call error(1,'Error in subbmk1. Case not prepared: r .le. 0 ')

!     computes Bessel functions: to prevent overflow  and  underflow, they are redefined as:
! 
!              bi(n) = (n+1/2)! (2/z)^(n+1/2) Exp[-z]  BesselI[n+1/2,z]
!              bk(n) = Exp[z] (z/2)^(n+1/2) BesselK[n+1/2,z] / (n+1/2)!

!     functions I(n+1/2,z)(n+1/2)!: descending recursion

    bia = 0.d0
    bib = 0.d0
    argia = xxa * r
    argib = xxb * r
    indmx = ktop
    aux = 1.d0
    aux2 = 1.d0
    bux = .25d0 * argia * argia
    cux = 1.d0
    cux2 = 1.d0
    dux = .25d0 * argib * argib
    jmax1 = jmax + 1
    imxjmax = indmx+jmax
    imxjmax1 = indmx+jmax + 1
    do j = 1, jmax
        aux = 1.d0 + aux * bux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax1-j)+1) )
        aux2 = 1.d0 + aux2 * bux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax-j)+1) )
        cux = 1.d0 + cux * dux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax1-j)+1) )
        cux2 = 1.d0 + cux2 * dux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax-j)+1) )
    enddo
    bux = exp(-argia)
    bia(indmx) = aux * bux
    bia(indmx-1) = aux2 * bux
    bux = exp(-argib)
    bib(indmx) = cux * bux
    bib(indmx-1) = cux2 * bux
    argsqa = argia * argia
    argsqb = argib * argib

    do i = indmx-1, 0, -1
        bia(i-1) = argsqa * cfbk0(i) * bia(i+1) + bia(i)
        bib(i-1) = argsqb * cfbk0(i) * bib(i+1) + bib(i)
    enddo

    auxa = argia / (xxa * rac)
    buxa = exp(argia-xxa*rac) / sqrt(auxa)
    auxb = argib / (xxb * rbc)
    buxb = exp(argib-xxb*rbc) / sqrt(auxb)
    do k = -1, ktop
        bma(k) = buxa * bka(k) * bia(k)
        bmb(k) = buxb * bkb(k) * bib(k)
        buxa = buxa * auxa
        buxb = buxb * auxb
    enddo
    return
    end subroutine subbmk1
! 
! ************************************************************************
! 
!     subroutine for computing Bessel I functions and M funcions for translating STO from A and B to C in case of
!           min(Ra,Rb) <= r < max(Ra,Rb)
! 
  subroutine subbmk2(xxa, xxb, ktop, rac, rbc, r, bka, bia, bkb, bib, bma, bmb)    
    USE Zernike_Jacobi_STO_MPI_D
    implicit none
    integer(KINT) :: i, imxjmax, imxjmax1, indmx, j, jmax1, k, ksup, ktop
    real(KREAL) :: argib, arginvka, argka, argsqa, argsqb, auxa, auxb, bux, buxa, buxb, cux, cux2
    real(KREAL) :: dux, r, rac, rbc, xxa, xxb
    real(KREAL) :: bia(-1:mxk), bib(-1:mxk), bka(-1:mxk), bkb(-1:mxk), bma(-3:mxk), bmb(-3:mxk)

!      computes Bessel functions: to prevent
!      overflow  and  underflow, they are redefined as:
! 
!              bi(n) = (n+1/2)! (2/z)^(n+1/2) Exp[-z]  BesselI[n+1/2,z]
!              bk(n) = Exp[z] (z/2)^(n+1/2) BesselK[n+1/2,z] / (n+1/2)!

!      functions I(n+1/2,z)(n+1/2)!: descending recursion

    bka = 0.d0
    bib = 0.d0
    bma = 0.d0
    bmb = 0.d0
    argib = xxb * r
    indmx = ktop
    cux = 1.d0
    cux2 = 1.d0
    dux = .25d0 * argib * argib
    jmax1 = jmax + 1
    imxjmax = indmx+jmax
    imxjmax1 = indmx+jmax + 1
    do j = 1, jmax
        cux = 1.d0 + cux * dux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax1-j)+1) )
        cux2 = 1.d0 + cux2 * dux * 2.d0 / ( re(jmax1-j) * re(2*(imxjmax-j)+1) )
    enddo

    bux = exp(-argib)
    bib(indmx) = cux  * bux
    bib(indmx-1) = cux2 * bux
    argsqb = argib * argib
    do i = indmx-1, 0, -1
        bib(i-1) = argsqb * cfbk0(i) * bib(i+1) + bib(i)
    enddo

!     functions K(n+1/2,z)/(n+1/2)! : ascending recursion
!     and functions M(n+1/2,z) = I(n+1/2,z_<) * K(n+1/2,z_>)

    argka = xxa * r
    arginvka = 1.d0 / argka

    bka(-1) = 1.d0 / argka
    bka(0) = 1.d0
    argsqa = argka * argka

    do i = 0, indmx-1
        bka(i+1) = argsqa * cfbk0(i) * bka(i-1) + cfbk1(i) * bka(i)
    enddo

    auxa = (xxa * rac) / argka
    buxa = exp(xxa*rac-argka) / sqrt(auxa)
    auxb = argib / (xxb * rbc)
    buxb = exp(argib-xxb*rbc) / sqrt(auxb)
    ksup = min(ktop, nint(-(300.d0+xxa*rac-argka) / log(auxa)))
    do k = -1, ksup
        bma(k) = buxa * bka(k) * bia(k)
        buxa = buxa * auxa
    enddo
    ksup = min(ktop, nint(-(300.d0+argib-xxb*rbc) / log(auxb)))
    do k = -1, ksup
        bmb(k) = buxb * bkb(k) * bib(k)
        buxb = buxb * auxb
    enddo
    return
    end subroutine subbmk2
! 
! ************************************************************************
! 
!     subroutine for computing Bessel I functions and M funcions for translating STO from A and B to C in case of
!           r > max(Ra,Rb)
! 
  subroutine subbmk3(xxa, xxb, ktop, rac, rbc, r, bka, bia, bkb, bib, bma, bmb) 
    USE Zernike_Jacobi_STO_MPI_D
    implicit none
    integer(KINT) :: i, indmxk, k, ksup, ktop
    real(KREAL) :: argka, argkb, argsqa, argsqb, auxa, auxb, buxa, buxb
    real(KREAL) :: r, rac, rbc, xxa, xxb
    real(KREAL) :: bia(-1:mxk), bib(-1:mxk), bka(-1:mxk), bkb(-1:mxk), bma(-3:mxk), bmb(-3:mxk)

!     computes Bessel functions: to prevent
!     overflow  and  underflow, they are redefined as:
! 
!              bi(n) = (n+1/2)! (2/z)^(n+1/2) Exp[-z]  BesselI[n+1/2,z]
!              bk(n) = Exp[z] (z/2)^(n+1/2) BesselK[n+1/2,z] / (n+1/2)!
!
!     functions K(n+1/2,z)/(n+1/2)! : ascending recursion
!     and functions M(n+1/2,z) = I(n+1/2,z_<) * K(n+1/2,z_>)

    bka = 0.d0
    bkb = 0.d0
    bma = 0.d0
    bmb = 0.d0
    indmxk = ktop
    argka = xxa * r
    bka(-1) = 1.d0 / argka
    bka(0) = 1.d0
    argkb = xxb * r
    bkb(-1) = 1.d0 / argkb
    bkb(0) = 1.d0
    argsqa = argka * argka
    argsqb = argkb * argkb
    do i = 0, indmxk-1
        bka(i+1) = argsqa * cfbk0(i) * bka(i-1) + cfbk1(i) * bka(i)
        bkb(i+1) = argsqb * cfbk0(i) * bkb(i-1) + cfbk1(i) * bkb(i)
    enddo
    auxa = (xxa * rac) / argka
    buxa = exp(xxa*rac-argka) / sqrt(auxa)
    auxb = (xxb * rbc) / argkb
    buxb = exp(xxb*rbc-argkb) / sqrt(auxb)
    ksup = min(ktop, nint(-(300.d0+xxa*rac-argka) / log(auxa)))
    do k = -1, ksup
        bma(k) = buxa * bka(k) * bia(k)
        buxa = buxa * auxa
    enddo
    ksup = min(ktop, nint(-(300.d0+xxb*rbc-argkb) / log(auxb)))
    do k = -1, ksup
        bmb(k) = buxb * bkb(k) * bib(k)
        buxb = buxb * auxb
    enddo
    return
    end subroutine subbmk3
! 
! ************************************************************************
! 
!     subrutine   subflm   for computing the radial factors of two-center
!     basic distributions:
!          1s-1s, 1s-2s, 2s-1s o 2s-2s
!     the basic function is chosen according to the parity of
!     na and nb:
!               na        nb     basic distribution
!              odd       odd           1s-1s
!              odd       even          1s-2s
!              even      odd           2s-1s
!              even      even          2s-2s
! 
  subroutine subflm(na, napar, la, xxa, Ra, nb, nbpar, lb, xxb, Rb, xb, zb, plm, bma, bmb, lmaxini, ktop, kmax, r) 
    USE Zernike_Jacobi_STO_MPI_D
    implicit none

    logical napar, nbpar
    integer(KINT) :: i, ip, k, l, la, lb, ldim, lm, lmax, lmaxaux, lmaxbux, lmaxini, kmax, ktop, m, n, n0, na, nb, np, np0
    real(KREAL) :: bma(-3:mxk), bmb(-3:mxk), fka(-lmaxini:mxk), fkb(0:mxk), flmbas(-1:lmaxini+2,-1:lmaxini)
    real(KREAL) :: flmm1(-1:lmaxini), flmnnp((4+lmaxini)**2,mxl+1,mxl+1), plm(-1:mxk,0:ldimaux)
    real(KREAL) :: aux, c1a, c1b, c2a, c2b, ca, cb, dlt, dosza, doszb, fl, flm1, r, r2, ra, rb, rinv, xb, xxa, xxb, zb

    lmax = lmaxini
    ldim = 4+lmaxini

    do i = 1, lmax
        fka(-i) = 0.d0
    enddo
    if (Ra .gt. toldstorig) then
            if (napar) then   ! Case na even: computes fka for n = 2
                n0 = 2
                c1a = r * Ra
                c2a = (r*r+Ra*Ra)
                ca = 1.d0 / sqrt(r * Ra)
                do k = 0, ktop-1
                    fka(k) = ca * ( cffk21(k) * c1a * bma(k-1) + cffk22(k) * c2a * bma(k) + cffk23(k) * c1a * bma(k+1) )
                enddo
            else              ! Case na odd: computes fka for n = 1
                n0 = 1
                ca = xxa * sqrt(r * Ra)
                do k = 0, ktop-1
                    fka(k) =  ca * ( bma(k-1) - bma(k+1) )
                enddo
            endif
    else
            do k = 0, ktop-1
                fka(k) = 0.d0
            enddo
            if (napar) then
                n0 = 2
                fka(0) = r * exp(-xxa * r)
            else
                n0 = 1
                fka(0) = exp(-xxa * r)
            endif
    endif

    if (Rb .gt. toldstorig) then
            if (nbpar) then   ! Case nb even: computes fkb for n' = 2
                np0 = 2
                c1b = r * Rb
                c2b = (r*r+Rb*Rb)
                cb = 1.d0 / sqrt(r * Rb)
                do k = 0, ktop-1
                    fkb(k) = cb * ( cffk21(k) * c1b * bmb(k-1) + cffk22(k) * c2b * bmb(k) + cffk23(k) * c1b * bmb(k+1) )
                enddo
            else              ! Case nb odd: computes fkb for n' = 1
                np0 = 1
                cb = xxb * sqrt(r * Rb)
                do k = 0, ktop-1
                    fkb(k) = cb * ( bmb(k-1) - bmb(k+1) )
                enddo
            endif
    else
            do k = 0, ktop-1
                fkb(k) = 0.d0
            enddo
            if (nbpar) then
                np0 = 2
                fkb(0) = r * exp(-xxb * r)
            else
                np0 = 1
                fkb(0) = exp(-xxb * r)
            endif
    endif

!      Computes the radial factors of the basic distribution
!      Notice that array flmbas contains indices physically meaningless 
!      (l = -1). In these cases, loads zeroes.
!      This is so to prevent errors in some if sentences in the 
!      recursion on n'

    flmbas = 0.d0
    rinv = 1.d0 / r
    call subflmbas0(ldimaux,kmax,lmax,rinv,plm,fka,fkb,flmbas)
!     Recursion on n up to na
    r2 = r*r
    aux = r2 + Ra*Ra
    dosza = Ra+Ra
    do n = n0+2, na, 2
        lmax = lmax - 1
        do m = 0, lmax
            flm1 = 0.d0
            fl   = flmbas(m,m)
            do l = m, lmax
                flmbas(m,l) = aux * fl  - dosza * ( re(l-m) * ri(l+l-1) * flm1 + r2 * re(l+m+1) * ri(l+l+3) * flmbas(m,l+1) )
                flm1 = fl
                fl = flmbas(m,l+1)
            enddo
        enddo
    enddo

!     Recursion on n' up to nb

    aux = r2 + Rb*Rb
    doszb = Zb+Zb
    do np = np0+2, nb, 2
        do l = -1, lmax
                flmm1(l) = 0.d0
        enddo
        lmax = lmax - 1
        do m = 0, lmax
            flm1 = 0.d0
            fl   = flmbas(m,m)
            dlt = 1.d0
            if (m .eq. 1) dlt = 2.d0
            do l = m, lmax
                flmbas(m,l) = aux * fl  - doszb * ( re(l-m) * ri(l+l-1) * flm1 + r2 * re(l+m+1) * ri(l+l+3) * flmbas(m,l+1) ) &
                        - xb * ( dlt * (ri(l+l-1) * flmm1(l-1) - r2 * ri(l+l+3) * flmm1(l+1) ) - re(l-m) * re(l-m-1) * ri(l+l-1) &
                                        * flmbas(m+1,l-1) + re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flmbas(m+1,l+1) )
                flmm1(l-1) = flm1
                flm1 = fl
                fl = flmbas(m,l+1)
            enddo
            flmm1(lmax) = flm1
            flmm1(lmax+1) = fl
    enddo
    enddo
    if (la .eq. 0 .and. lb .eq. 0) then  ! loads radial factors and return
        lm = 0
        do l = 0, lmax
            do m = -l, -1
                lm = lm + 1
                flmmamb(lm,0,0) = 0.d0
            enddo
            do m = 0, l
                lm = lm + 1
                flmmamb(lm,0,0) = flmbas(m,l)
            enddo
        enddo
        return
    endif

!     If .not. (la .eq. 0 .and. lb .eq. 0)  recursion on la and lb
!     It starts with a recurrence on n and n' as necessary

! 	flmnnp = 0.d0
    lm = 0
    do l = 0, lmax
        do m = 0, l
            lm = lm + 1
            flmnnp(lm,1,1) = flmbas(m,l)
        enddo
    enddo

!     recursion on  n  for the subsequent recursion on L
!     notice that recursion on n  runs in steps of two 2
!     whereas the storage index runs in steps of one, thus:
!     n = n0 + 2  (i-1)

    lmaxaux = lmax

    aux = r2 + ra*ra
    dosza = ra+ra
    do i = 2, la+1
        lmaxaux = lmaxaux - 1
        flmnnp(1,i,1) = aux * flmnnp(1,i-1,1)  - dosza * r2 * ri(3) * flmnnp(2,i-1,1)
        lm = 1
        do l = 1, lmaxaux
            do m = 0, l-1
                lm = lm + 1
                flmnnp(lm,i,1) = aux * flmnnp(lm,i-1,1)  - dosza * ( re(l-m) * ri(l+l-1) * flmnnp(ind(l-1)+m+1,i-1,1) &
                        + r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+m+1,i-1,1) )
            enddo
            lm = lm + 1
            flmnnp(lm,i,1) = aux * flmnnp(lm,i-1,1)  - dosza * r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+l+1,i-1,1)
        enddo
    enddo

!     recursion on  np  for the subsequent recursion on Lp
!     notice that recursion on np  runs in steps of two 2
!     whereas the storage index runs in steps of one, thus:
!     np = np0 + 2  (ip-1)

    aux = r2 + rb*rb
    doszb = zb+zb
    lmaxbux = lmax+1
    do ip = 2, lb+1
        lmaxbux = lmaxbux - 1
        lmaxaux = lmaxbux
        do i = 1, la+1
            lmaxaux = lmaxaux-1
!           l = 0  m = 0
            flmnnp(1,i,ip) = aux * flmnnp(1,i,ip-1) - doszb * r2 * ri(3) * flmnnp(2,i,ip-1) &
                    - xb * 0.6666666666666667d0 * r2 * flmnnp(3,i,ip-1)
                if (lmaxaux .eq. 0) cycle
!           l = 1  m = 0
            flmnnp(2,i,ip) = aux * flmnnp(2,i,ip-1) - doszb * ( flmnnp(1,i,ip-1) + r2 * .4d0 * flmnnp(4,i,ip-1) ) &
                    - xb * 1.2d0 * r2 * flmnnp(5,i,ip-1)
!           l = 1  m = 1
            flmnnp(3,i,ip) = aux * flmnnp(3,i,ip-1) - doszb * r2 * .6d0 * flmnnp(5,i,ip-1) - xb * ( 2.d0 * (flmnnp(1,i,ip-1) &
                    - r2 * .2d0 * flmnnp(4,i,ip-1)) + 2.4d0 * r2 * flmnnp(6,i,ip-1) )
            lm = 3
            do l = 2, lmaxaux
                lm = lm + 1
!              m = 0
                flmnnp(lm,i,ip) = aux * flmnnp(lm,i,ip-1) - doszb * ( re(l) * ri(l+l-1) * flmnnp(ind(l-1)+1,i,ip-1) &
                + r2 * re(l+1) * ri(l+l+3) * flmnnp(ind(l+1)+1,i,ip-1) ) - xb * (- re(l) * re(l-1) * ri(l+l-1) &
                        * flmnnp(ind(l-1)+2,i,ip-1) + re(l+1) * re(l+2) * ri(l+l+3) * r2 * flmnnp(ind(l+1)+2,i,ip-1) )
                dlt = 2.d0
                do m = 1, l-1
                    lm = lm + 1
                    flmnnp(lm,i,ip) = aux * flmnnp(lm,i,ip-1) - doszb * ( re(l-m) * ri(l+l-1) * flmnnp(ind(l-1)+m+1,i,ip-1) &
                            + r2 * re(l+m+1) * ri(l+l+3)  * flmnnp(ind(l+1)+m+1,i,ip-1) ) &
                            - xb * ( dlt * (ri(l+l-1) * flmnnp(ind(l-1)+m,i,ip-1) - r2 * ri(l+l+3) * flmnnp(ind(l+1)+m,i,ip-1)) &
                                    - re(l-m) * re(l-m-1) * ri(l+l-1) * flmnnp(ind(l-1)+m+2,i,ip-1) &
                                    + re(l+m+1) * re(l+m+2) * ri(l+l+3) * r2 * flmnnp(ind(l+1)+m+2,i,ip-1) )
                    dlt = 1.d0
                enddo
                lm = lm + 1
                flmnnp(lm,i,ip) = aux * flmnnp(lm,i,ip-1) - doszb * (r2 * re(l+l+1) * ri(l+l+3) * flmnnp(ind(l+1)+l+1,i,ip-1))&
                    - xb * ( ri(l+l-1) * flmnnp(ind(l-1)+l,i,ip-1) - r2 * ri(l+l+3) * flmnnp(ind(l+1)+l,i,ip-1) &
                        + re(l+l+1) * re(l+l+2) * ri(l+l+3) * r2 * flmnnp(ind(l+1)+l+2,i,ip-1) )
            enddo   ! End of Do on l
        enddo   ! End of Do on i
    enddo   ! End of Do on ip
    call recurflm(ldim, lmax, la, lb, r2, Ra, xb, zb, flmnnp)
    return
    end subroutine subflm
! 
!*******************************************************
!
!     subroutine for computing radial factors of two-center distributions with l = l' = 0
!     from the radial factors of the translation of individual STOs
! 
  subroutine subflmbas0(ldimaux, kmax, lmax, rinv, plm, fka, fkb, flmbas)
    USE Zernike_Jacobi_STO_MPI_D, only: KINT, KREAL, mxk, mxlckplm, mxkextra, mxlckplmextra, ckplm, ckplmextra
    implicit none
    integer(KINT) :: ip2, k, kmax, knt1, knt2, l, ldimaux, lmax, m
    real(KREAL) :: aux, bux, rinv, rinvl, sgn
    real(KREAL) :: fka(-lmax:mxk), fkb(0:mxk), flmbas(-1:lmax+2,-1:lmax), plm(-1:mxk,0:ldimaux)
    knt1 = 0
    rinvl = 1.d0
    do l = 0, min(lmax,mxlckplm)
        sgn = 1.d0
        do m = 0, l
            aux = 0.d0
            do k = m, kmax
                knt2 = knt1
                bux = 0.d0
                do ip2 = 0, l+l, 2
                    knt2 = knt2 + 1
                    bux = bux + ckplm(k,knt2) * fka(l+k-ip2)
                enddo
                aux = aux + plm(k,m) * fkb(k) * bux
            enddo
            flmbas(m,l) = sgn * aux * rinvl
            sgn = -sgn
            knt1 = knt2
        enddo
        rinvl = rinvl * rinv
    enddo
    if (lmax .gt. mxlckplm) then
        knt1 = 0
        do l = mxlckplm+1, min(lmax,mxlckplmextra)
            sgn = 1.d0
            do m = 0, l
                aux = 0.d0
                do k = m, min(kmax,mxkextra)
                    knt2 = knt1
                    bux = 0.d0
                    do ip2 = 0, l+l, 2
                        knt2 = knt2 + 1
                        bux = bux + ckplmextra(k,knt2) * fka(l+k-ip2)
                    enddo
                    aux = aux + plm(k,m) * fkb(k) * bux
                enddo
                flmbas(m,l) = sgn * aux * rinvl
                sgn = -sgn
                knt1 = knt2
            enddo
            rinvl = rinvl * rinv
        enddo
    endif
    if (lmax .gt. mxlckplmextra) then
        write(6,"(//'WARNING !!! lmax = ', i3,' greater than maximum available (',i3,') in subroutine subflmbas0. \n&
        &Loads zeroes in the remaining elements of flmbas',//)") lmax, mxlckplmextra
        do l = mxlckplmextra+1, lmax
            flmbas(0:l,l) = 0.d0
        enddo
    endif
    return
    end
! 
! ***************************************************************
! 
!     Subrutine for recursion of the radial factors of the distribution
!     Indices (l,m) are contracted to a single one:
!         lm = l2 + l + m + 1
! 
  subroutine recurflm(ldim, lmaxtot, la, lb, r2, za, xb, zb, fin)
    USE Zernike_Jacobi_STO_MPI_D, fvoid => flm
    implicit none
    logical :: lcux
    integer(KINT) :: ierr, knt, l, la, lb, ldim, lg, lgg, lm, lmax, lmaxaux, lmaxbux, lmaxtot, lmpos
    integer(KINT) :: m, ma, mb, mm, n, nmax, np, npmax
    real(KREAL) :: fin(ldim*ldim,mxl+1,mxl+1), flmaux(ldim*ldim)
    real(KREAL), allocatable :: flm(:,:,:,:,:)
    real(KREAL) :: aux, aux2, bux, cux, r2, s1, s1m, s2, s2m, s3, s3m, s4, s4m, s5, s5m
    real(KREAL) :: s6, s6m, s7, s7m, s8, s8m, sxnn, sxnp, sxpn, sxpp, xb, za, zb, umdltm0, updltm1, umdltm1

    lmax = lmaxtot
! 	nmax = max(1,la)
! 	npmax = max(1,lb)
    nmax = la+1
    npmax = lb+1
    allocate(flm(ldim*ldim,la+1,lb+1,-(la+1):(la+1),-(lb+1):(lb+1)), stat = ierr)
    if (ierr .ne. 0) call error(ierr,'Memory error when allocating flm. Stop')

!     loads the content of  fin  into  flm(lm,n,np,0,0)
    lmaxbux = lmax
    do np = 1, npmax
        lmaxaux = lmaxbux
        do n = 1, nmax
            lm = 0
            lmpos = 0
            do l = 0, lmaxaux
                do m = -l, -1
                    lm = lm + 1
                    flm(lm,n,np,0,0) = 0.d0
                enddo
                do m = 0, l
                    lm = lm + 1
                    lmpos = lmpos + 1
                    flm(lm,n,np,0,0) = fin(lmpos,n,np)
                enddo
            enddo
            lmaxaux = lmaxaux -1
        enddo
        lmaxbux = lmaxbux - 1
    enddo

!  =========================
!      recursion on L    
!  =========================

    if (la .ge. 1) then
        do lg = 1, la

!           recursion of factors corresponding to  n = 1

            aux = real(lg+lg-1) * 0.5d0
            aux2 = aux+aux
            lmaxaux = lmax-lg+1
            n = 1
            do np = 1, npmax

!              elements with M = lg  and  M = -lg
                lmaxaux = lmaxaux-1
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
                            s1 = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,lg-1,0)
                            s1m = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,lg-1,0)
                        endif
                        if (lg .gt. 1) then
                            s5 = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,-lg+1,0)
                            s5m = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,-lg+1,0)
                            s7 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,-lg+1,0)
                            s7m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,-lg+1,0)
                        endif
                        s3 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,lg-1,0)
                        s3m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,lg-1,0)
                        if (l .gt. 0) then
                            if (m .gt. 0) then
                                s2 = ri(l+l-1) * flm((l-1)*l+m,n,np,lg-1,0)
                                if (lg .gt. 1) s6m = ri(l+l-1) * flm((l-1)*l+m,n,np,-lg+1,0)
                            endif
                            if (m .gt. 1) then
                                s2m = ri(l+l-1) * flm((l-1)*l-m+2,n,np,lg-1,0)
                                if (lg .gt. 1) s6 = ri(l+l-1) * flm((l-1)*l-m+2,n,np,-lg+1,0)
                            endif
                            if (m .lt. l-1) then
                                s4 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,lg-1,0)
                                s4m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,lg-1,0)
                                if (lg .gt. 1) then
                                    s8 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,-lg+1,0)
                                    s8m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,-lg+1,0)
                                endif
                            endif
                        endif
                        if (m .gt. 0) then
                            flm(l*l+l-m+1,n,np,lg,0) = aux * (umdltm1*umdltm0*(-s1m+s2m) + s3m - s4m &
                                    - updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                            flm(l*l+l-m+1,n,np,-lg,0) = aux *(updltm1*umdltm0*(-s1+s2) -s3 +s4 + umdltm1*umdltm0*(-s5+s6) &
                                    + s7 - s8)
                        endif
                        flm(l*l+l+m+1,n,np,lg,0) = aux *(updltm1*umdltm0*(-s1+s2) + s3 - s4 - umdltm1*umdltm0*(s5-s6) &
                                - s7 + s8 )
                        flm(l*l+l+m+1,n,np,-lg,0) =aux*(umdltm1*umdltm0*(s1m-s2m) +s3m -s4m + updltm1*umdltm0*(-s5m+s6m) &
                                + s7m - s8m)
                    enddo     ! End of Do on m
                enddo     ! End of Do on l

!     		elements with   -lg < M < lg
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
                            s1 = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l+m+3,n,np,mm,0)
                            s1m = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l-m+3,n,np,mm,0)
                            s2 = 0.d0
                            s2m = 0.d0
                            if (m .lt. l) then
                                s2 = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l+m,n,np,mm,0)
                                s2m = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l-m,n,np,mm,0)
                            endif
                            if (m .gt. 0) then
                                flmaux(l*l+l-m+1) = bux * aux2 * (s1m+s2m - za * flm(l*l+l-m+1,n,np,mm,0) )
                                if (lcux) flmaux(l*l+l-m+1) = flmaux(l*l+l-m+1) - cux * flm(l*l+l-m+1,n+1,np,mm,0)
                            endif
                            flmaux(l*l+l+m+1) = bux * aux2 * (s1+s2 - za * flm(l*l+l+m+1,n,np,mm,0))
                            if (lcux) flmaux(l*l+l+m+1) = flmaux(l*l+l+m+1) - cux * flm(l*l+l+m+1,n+1,np,mm,0)
                        enddo     ! End of Do on m
                    enddo     ! End of Do on l
                    do l  = 0, lmaxaux
                        flm(l*l+l+1,n,np,mm,0) = flmaux(l*l+l+1)
                        do m  = 1, l
                            flm(l*l+l+m+1,n,np,mm,0) = flmaux(l*l+l+m+1)
                            flm(l*l+l-m+1,n,np,mm,0) = flmaux(l*l+l-m+1)
                        enddo     ! End of Do on m
                    enddo     ! End of Do on l
                enddo     ! End of Do on mm
            enddo     ! End of Do on np

            if (lg .eq. la) exit  ! exits do on lg

!           recursion of factors with  n > 1

            if (lg .gt. 1) then
                lgg = lg
                do n  = 2, min(la,lg,la+1-lg)
                    lgg = lgg - 1
                    aux = real(lgg+lgg-1) * 0.5d0
                    aux2 = aux + aux
                    lmaxaux = lmax-lg-n+2
                    do np = 1, npmax

!					elements with M = lgg  y  M = -lgg
                        lmaxaux = lmaxaux-1
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
                                    s1 = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,lgg-1,0)
                                    s1m = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,lgg-1,0)
                                endif
                                if (lgg .gt. 1) then
                                    s5 = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,-lgg+1,0)
                                    s5m = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,-lgg+1,0)
                                    s7 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,-lgg+1,0)
                                    s7m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,-lgg+1,0)
                                endif
                                s3 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,lgg-1,0)
                                s3m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,lgg-1,0)
                                if (l .gt. 0) then
                                    if (m .gt. 0) then
                                        s2 = ri(l+l-1) * flm((l-1)*l+m,n,np,lgg-1,0)
                                        if (lgg .gt. 1) then
                                            s6m = ri(l+l-1) * flm((l-1)*l+m,n,np,-lgg+1,0)
                                        endif
                                    endif
                                    if (m .gt. 1) then
                                        s2m = ri(l+l-1) * flm((l-1)*l-m+2,n,np,lgg-1,0)
                                        if (lgg .gt. 1) then
                                            s6 = ri(l+l-1) * flm((l-1)*l-m+2,n,np,-lgg+1,0)
                                        endif
                                    endif
                                    if (m .lt. l-1) then
                                        s4 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,lgg-1,0)
                                        s4m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,lgg-1,0)
                                        if (lgg .gt. 1) then
                                            s8 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,-lgg+1,0)
                                            s8m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,-lgg+1,0)
                                        endif
                                    endif
                                endif
                                if (m .gt. 0) then
                                    flm(l*l+l-m+1,n,np,lgg,0) = aux  * (umdltm1*umdltm0*(-s1m+s2m) + s3m - s4m &
                                            - updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                                    flm(l*l+l-m+1,n,np,-lgg,0) = aux  * (updltm1*umdltm0*(-s1+s2) - s3 + s4 &
                                            + umdltm1*umdltm0*(-s5+s6) + s7 - s8 )
                                endif
                                flm(l*l+l+m+1,n,np,lgg,0) =aux*(updltm1*umdltm0*(-s1+s2) + s3 - s4 &
                                        - umdltm1*umdltm0*(s5-s6) - s7 + s8)
                                flm(l*l+l+m+1,n,np,-lgg,0)=aux*(umdltm1*umdltm0*(s1m-s2m) + s3m - s4m &
                                        + updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                            enddo     ! End of Do on m
                        enddo     ! End of Do on l

!					elements with   -lgg < M < lgg
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
                                    s1 = re(l+m+1) * ri(l+l+3) * r2  * flm((l+3)*l+m+3,n,np,mm,0)
                                    s1m = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l-m+3,n,np,mm,0)
                                    s2 = 0.d0
                                    s2m = 0.d0
                                    if (m .lt. l) then
                                        s2 = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l+m,n,np,mm,0)
                                        s2m = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l-m,n,np,mm,0)
                                    endif
                                    if (m .gt. 0) then
                                        flmaux(l*l+l-m+1) = bux * aux2 * (s1m+s2m - za * flm(l*l+l-m+1,n,np,mm,0) )
                                        if (lcux) flmaux(l*l+l-m+1) = flmaux(l*l+l-m+1)-cux * flm(l*l+l-m+1,n+1,np,mm,0)
                                    endif
                                    flmaux(l*l+l+m+1) = bux * aux2 * (s1+s2 - za * flm(l*l+l+m+1,n,np,mm,0) )
                                    if (lcux) flmaux(l*l+l+m+1) = flmaux(l*l+l+m+1) - cux * flm(l*l+l+m+1,n+1,np,mm,0)
                                enddo     ! End of Do on m
                            enddo     ! End of Do on l
                            do l  = 0, lmaxaux
                                flm(l*l+l+1,n,np,mm,0) = flmaux(l*l+l+1)
                                do m  = 1, l
                                    flm(l*l+l+m+1,n,np,mm,0) = flmaux(l*l+l+m+1)
                                    flm(l*l+l-m+1,n,np,mm,0) = flmaux(l*l+l-m+1)
                                enddo     ! End of Do on m
                            enddo     ! End of Do on l
                        enddo     ! End of Do on mm
                    enddo     ! End of Do on np
                enddo     ! End of Do on n
            endif     ! End of  if (lg .gt. 1)
        enddo     ! End of Do on lg
    endif  ! End of if (la .ge. 1)

    lmax = lmax - la

! =========================
!     recursion on L'   
! =========================

    if (lb .ge. 1) then
        do lg = 1, lb

!           recursion of factors with  np = 1

            aux = re(lg+lg-1) * 0.5d0
            aux2 = aux+aux
            lmaxaux = lmax-lg

            n = 1
            np = 1
!		elements with M = lg  y  M = -lg
            do ma = -la, la
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
                        sxpp = 0.d0      ! xb * f(l,m,n,L,M,n',L',L')
                        sxnp = 0.d0      ! xb * f(l,-m,n,L,M,n',L',L')
                        s1m = 0.d0
                        s2m = 0.d0
                        s3m = 0.d0
                        s4m = 0.d0
                        s5m = 0.d0
                        s6m = 0.d0
                        s7m = 0.d0
                        s8m = 0.d0
                        sxpn = 0.d0      ! xb * f(l,m,n,L,M,n',L',-L')
                        sxnn = 0.d0      ! xb * f(l,-m,n,L,M,n',L',-L')
                        umdltm0 = 1.d0               ! umdltm0 = 1 - delta(m,0)
                        updltm1 = 1.d0               ! updltm1 = 1 + delta(m,1)
                        umdltm1 = 1.d0               ! umdltm1 = 1 - delta(m,1)
                        if (m .eq. 0) umdltm0 = 0.d0
                        if (m .eq. 1) then
                            umdltm1 = 0.d0
                            updltm1 = 2.d0
                        endif
                        if (m .gt. 0) then
                            s1 = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,ma,lg-1)
                            s1m = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,ma,lg-1)
                        endif
                        if (lg .gt. 1) then
                            s5 = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,ma,-lg+1)
                            s5m = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,ma,-lg+1)
                            s7 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,ma,-lg+1)
                            s7m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,ma,-lg+1)
                        endif
                        s3 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,ma,lg-1)
                        s3m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,ma,lg-1)
                        if (l .gt. 0) then
                            if (m .gt. 0) then
                                s2 = ri(l+l-1) * flm((l-1)*l+m,n,np,ma,lg-1)
                                if (lg .gt. 1) s6m = ri(l+l-1) * flm((l-1)*l+m,n,np,ma,-lg+1)
                            endif
                            if (m .gt. 1) then
                                s2m = ri(l+l-1) * flm((l-1)*l-m+2,n,np,ma,lg-1)
                                if (lg .gt. 1) s6 = ri(l+l-1) * flm((l-1)*l-m+2,n,np,ma,-lg+1)
                            endif
                            if (m .lt. l-1) then
                                s4 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,ma,lg-1)
                                s4m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,ma,lg-1)
                                if (lg .gt. 1) then
                                    s8 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,ma,-lg+1)
                                    s8m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,ma,-lg+1)
                                endif
                            endif
                        endif
                        sxpp = 2.d0 * xb * flm(l*l+l+m+1,n,np,ma,lg-1)
                        sxnp = 2.d0 * xb * flm(l*l+l-m+1,n,np,ma,lg-1)
                        if (lg .gt. 1) then
                            sxpn = 2.d0 * xb * flm(l*l+l+m+1,n,np,ma,-lg+1)
                            sxnn = 2.d0 * xb * flm(l*l+l-m+1,n,np,ma,-lg+1)
                        endif
                        if (m .gt. 0) then
                            flm(l*l+l-m+1,n,np,ma,lg) = aux * (umdltm1*umdltm0*(-s1m+s2m) + s3m - s4m - sxnp &
                                    - updltm1*umdltm0*(-s5m+s6m) + s7m - s8m)
                            flm(l*l+l-m+1,n,np,ma,-lg) = aux * (updltm1*umdltm0*(-s1+s2) - s3 + s4 &
                                    + umdltm1*umdltm0*(-s5+s6) + s7 - s8 - sxnn )
                        endif
                        flm(l*l+l+m+1,n,np,ma,lg) = aux*(updltm1*umdltm0*(-s1+s2) + s3 - s4 - sxpp - umdltm1*umdltm0*(s5-s6) &
                            - s7 + s8 )
                        flm(l*l+l+m+1,n,np,ma,-lg)=aux*(umdltm1*umdltm0*(s1m-s2m) + s3m -s4m &
                            + updltm1 * umdltm0 * (-s5m+s6m)+s7m-s8m-sxpn)
                    enddo     ! End of Do on m
                enddo     ! End of Do on l

!			elements with   -lg < M < lg
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
                            s1 = re(l+m+1) * ri(l+l+3) * r2  * flm((l+3)*l+m+3,n,np,ma,mm)
                            s1m = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l-m+3,n,np,ma,mm)
                            s2 = 0.d0
                            s2m = 0.d0
                            if (m .lt. l) then
                                s2 = re(l-m) * ri(l+l-1)  * flm((l-1)*(l-1)+l+m,n,np,ma,mm)
                                s2m = re(l-m) * ri(l+l-1)  * flm((l-1)*(l-1)+l-m,n,np,ma,mm)
                            endif
                            if (m .gt. 0) then
                                flmaux(l*l+l-m+1) = bux * aux2 * (s1m+s2m - zb * flm(l*l+l-m+1,n,np,ma,mm) )
                                if (lcux) flmaux(l*l+l-m+1) = flmaux(l*l+l-m+1) - cux * flm(l*l+l-m+1,n,np+1,ma,mm)
                            endif
                            flmaux(l*l+l+m+1) = bux * aux2 * (s1+s2 - zb * flm(l*l+l+m+1,n,np,ma,mm) )
                            if (lcux) flmaux(l*l+l+m+1) = flmaux(l*l+l+m+1) - cux * flm(l*l+l+m+1,n,np+1,ma,mm)
                        enddo     ! End of Do on m
                    enddo     ! End of Do on l
                    do l  = 0, lmaxaux
                        flm(l*l+l+1,n,np,ma,mm) = flmaux(l*l+l+1)
                        do m  = 1, l
                            flm(l*l+l+m+1,n,np,ma,mm) = flmaux(l*l+l+m+1)
                            flm(l*l+l-m+1,n,np,ma,mm) = flmaux(l*l+l-m+1)
                        enddo     ! End of Do on m
                    enddo     ! End of Do on l
                enddo     ! End of Do on mm
            enddo     ! End of Do on ma

            if (lg .eq. lb) exit

!           recurre los factores con  np > 1

            if (lg .gt. 1) then
                lgg = lg
                do np  = 2, min(lb,lg,lb+1-lg)
                    lgg = lgg - 1
                    aux = re(lgg+lgg-1) * 0.5d0
                    aux2 = aux + aux
                    lmaxaux = lmax-lg+1-np
                    do ma = -la, la

!					elements with M = lgg  y  M = -lgg
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
                                sxpp = 0.d0      ! xb * f(l,m,n,L,M,n',L',L')
                                sxnp = 0.d0      ! xb * f(l,-m,n,L,M,n',L',L')
                                s1m = 0.d0
                                s2m = 0.d0
                                s3m = 0.d0
                                s4m = 0.d0
                                s5m = 0.d0
                                s6m = 0.d0
                                s7m = 0.d0
                                s8m = 0.d0
                                sxpn = 0.d0      ! xb * f(l,m,n,L,M,n',L',-L')
                                sxnn = 0.d0      ! xb * f(l,-m,n,L,M,n',L',-L')
                                umdltm0 = 1.d0               ! umdltm0 = 1 - delta(m,0)
                                updltm1 = 1.d0               ! updltm1 = 1 + delta(m,1)
                                umdltm1 = 1.d0               ! umdltm1 = 1 - delta(m,1)
                                if (m .eq. 0) umdltm0 = 0.d0
                                if (m .eq. 1) then
                                    umdltm1 = 0.d0
                                    updltm1 = 2.d0
                                endif
                                if (m .gt. 0) then
                                    s1 = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,ma,lgg-1)
                                    s1m = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,ma,lgg-1)
                                endif
                                if (lgg .gt. 1) then
                                    s5 = r2 * ri(l+l+3) * flm((l+3)*l-m+4,n,np,ma,-lgg+1)
                                    s5m = r2 * ri(l+l+3) * flm((l+3)*l+m+2,n,np,ma,-lgg+1)
                                    s7 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,ma,-lgg+1)
                                    s7m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,ma,-lgg+1)
                                endif
                                s3 = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l+m+4,n,np,ma,lgg-1)
                                s3m = re(l+m+1) * re(l+m+2) * r2 * ri(l+l+3) * flm((l+3)*l-m+2,n,np,ma,lgg-1)
                                if (l .gt. 0) then
                                    if (m .gt. 0) then
                                        s2 = ri(l+l-1) * flm((l-1)*l+m,n,np,ma,lgg-1)
                                        if (lgg .gt. 1) then
                                            s6m = ri(l+l-1) * flm((l-1)*l+m,n,np,ma,-lgg+1)
                                        endif
                                    endif
                                    if (m .gt. 1) then
                                        s2m = ri(l+l-1) * flm((l-1)*l-m+2,n,np,ma,lgg-1)
                                        if (lgg .gt. 1) then
                                            s6 = ri(l+l-1) * flm((l-1)*l-m+2,n,np,ma,-lgg+1)
                                        endif
                                    endif
                                    if (m .lt. l-1) then
                                        s4 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,ma,lgg-1)
                                        s4m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,ma,lgg-1)
                                        if (lgg .gt. 1) then
                                            s8 = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l-m,n,np,ma,-lgg+1)
                                            s8m = re(l-m) * re(l-m-1) * ri(l+l-1) * flm((l-1)*l+m+2,n,np,ma,-lgg+1)
                                        endif
                                    endif
                                endif
                                sxpp = 2.d0 * xb * flm(l*l+l+m+1,n,np,ma,lgg-1)
                                sxnp = 2.d0 * xb * flm(l*l+l-m+1,n,np,ma,lgg-1)
                                if (lgg .gt. 1) then
                                    sxpn = 2.d0 * xb * flm(l*l+l+m+1,n,np,ma,-lgg+1)
                                    sxnn = 2.d0 * xb * flm(l*l+l-m+1,n,np,ma,-lgg+1)
                                endif
                                if (m .gt. 0) then
                                    flm(l*l+l-m+1,n,np,ma,lgg) = aux   * (umdltm1 * umdltm0 * (-s1m+s2m) + s3m - s4m &
                                            - sxnp - updltm1*umdltm0*(-s5m+s6m) + s7m - s8m )
                                    flm(l*l+l-m+1,n,np,ma,-lgg) = aux  * (updltm1*umdltm0*(-s1+s2) - s3 + s4 &
                                            + umdltm1*umdltm0*(-s5+s6) + s7 - s8 - sxnn )
                                endif
                                flm(l*l+l+m+1,n,np,ma,lgg) = aux*(updltm1*umdltm0*(-s1+s2) + s3 - s4 - sxpp &
                                        - umdltm1*umdltm0*(s5-s6) - s7 + s8 )
                                flm(l*l+l+m+1,n,np,ma,-lgg) = aux *(umdltm1*umdltm0*(s1m-s2m) + s3m - s4m &
                                        + updltm1*umdltm0*(-s5m+s6m) + s7m - s8m - sxpn)
                            enddo     ! End of Do on m
                        enddo     ! End of Do on l

!					elements with   -lgg < M < lgg
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
                                    s1 = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l+m+3,n,np,ma,mm)
                                    s1m = re(l+m+1) * ri(l+l+3) * r2 * flm((l+3)*l-m+3,n,np,ma,mm)
                                    s2 = 0.d0
                                    s2m = 0.d0
                                    if (m .lt. l) then
                                        s2 = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l+m,n,np,ma,mm)
                                        s2m = re(l-m) * ri(l+l-1) * flm((l-1)*(l-1)+l-m,n,np,ma,mm)
                                    endif
                                    if (m .gt. 0) then
                                        flmaux(l*l+l-m+1) = bux * aux2 * (s1m+s2m - zb * flm(l*l+l-m+1,n,np,ma,mm) )
                                        if (lcux) flmaux(l*l+l-m+1) = flmaux(l*l+l-m+1)- cux*flm(l*l+l-m+1,n,np+1,ma,mm)
                                    endif
                                    flmaux(l*l+l+m+1) = bux * aux2 * (s1+s2 - zb * flm(l*l+l+m+1,n,np,ma,mm) )
                                    if (lcux) flmaux(l*l+l+m+1) = flmaux(l*l+l+m+1) - cux * flm(l*l+l+m+1,n,np+1,ma,mm)
                                enddo     ! End of Do on m
                            enddo     ! End of Do on l
                            do l  = 0, lmaxaux
                                flm(l*l+l+1,n,np,ma,mm) = flmaux(l*l+l+1)
                                do m  = 1, l
                                    flm(l*l+l+m+1,n,np,ma,mm) = flmaux(l*l+l+m+1)
                                    flm(l*l+l-m+1,n,np,ma,mm) = flmaux(l*l+l-m+1)
                                enddo     ! End of Do on m
                            enddo     ! End of Do on l
                        enddo     ! End of Do on mm
                    enddo     ! End of Do on ma
                enddo     ! End of Do on np
            endif     ! End of  if (lg .gt. 1)
        enddo     ! End of Do on lg
    endif

    lmax = lmax - lb
    do mb = -lb, lb
        do ma = -la, la
        knt = 0
            do l = 0, lmax
                do m = -l, l
                    knt = knt + 1
                    flmmamb(knt,ma,mb) = flm(knt,1,1,ma,mb)
                enddo
            enddo
        enddo
    enddo
    return
    end subroutine recurflm
! 
! ************************************************************************
! 
  subroutine flmGTO(la, xxa, rac, lb, xxb, xb, zb, rbc, nga, ngb, lmaxini, lneglig, zlmij, cijv, r)
    USE Zernike_Jacobi_STO_MPI_D
    implicit none
    logical lneglig(mxgauss*mxgauss)
    integer(KINT) :: i, ij, ip, itop, ja, jb, l, la, lb, ldim, lm, lmax, lmaxaux, lmaxbux, lmaxini, lmp, lmxcab, lmmax
    integer(KINT) :: m, nga, ngb
    real(KREAL) :: bsi(0:45), cijv(mxgauss*mxgauss), flmbas((ldimaux+1)*(ldimaux+2)/2), flmnnp((4+lmaxini)**2,mxl+1,mxl+1)
    real(KREAL) :: zlmij((ldimaux+1)*(ldimaux+2)/2,mxgauss*mxgauss)
    real(KREAL) :: a2, arg1, arg2, aux, bux, dlt, dosxiij, dosza, doszb, r, r2, ra, rac, rb, rbc, rgij
    real(KREAL) :: rgijsq, rmrijsq, rprijsq, rtop, xb, xij, xiij, xxa, xxasq, xxb, xxbsq, za, zb, zij

    lmxcab = lmaxini
    ldim = 4+lmaxini
    za = rac
    xxasq = xxa * xxa
    xxbsq = xxb * xxb
    ij = 0
    lmmax = (lmxcab+1)*(lmxcab+2)/2
    do lm = 1, lmmax
        flmbas(lm) = 0.d0
    enddo
    do ja = 1, nga
        do jb = 1, ngb
            ij = ij + 1
            if (lneglig(ij)) cycle
            xiij = xia(ja) * xxasq + xib(jb) * xxbsq
            xij = xb * xib(jb) * xxbsq / xiij
            zij = za + (zb - za) * xib(jb) * xxbsq / xiij
            rgijsq = xij * xij + zij * zij
            rgij = sqrt(rgijsq)
            rmrijsq = (r-rgij) * (r-rgij)
            rprijsq = (r+rgij) * (r+rgij)
            arg1 = 0.5d0 * xiij * rprijsq
            arg2 = 0.5d0 * xiij * rmrijsq
            if (min(arg1,arg2) .gt. 100.d0) cycle
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
            dosxiij = xiij + xiij
            aux = cijv(ij)
            lm = 0
            do l = 0, lmxcab
                bux = aux * bsi(l) * re(l+l+1)
                lm = lm + 1
                flmbas(lm) = flmbas(lm) + bux * zlmij(lm,ij)
                bux = (bux + bux) * ri(l) * ri(l+1)
                do m = 1, l
                    lm = lm + 1
                    flmbas(lm) = flmbas(lm) + bux * zlmij(lm,ij)
                    bux = bux * ri(l-m) * ri(l+m+1)
                enddo
                aux = aux * dosxiij
            enddo
        enddo  ! End of do on jb
    enddo  ! End of do on ja
    lmax = lmxcab
    if (la .eq. 0 .and. lb .eq. 0) then
        lm = 0
        lmp = 0
        do l = 0, lmax
            do m = -l, -1
                lm = lm + 1
                flmmamb(lm,0,0) = 0.d0
            enddo
            do m = 0, l
                lm = lm + 1
                lmp = lmp + 1
                flmmamb(lm,0,0) = flmbas(lmp)
            enddo
        enddo
        return
    endif
!     if .not. (la .eq. 0 .and. lb .eq. 0)  recursion on la and lb
    lm = 0
    do l = 0, lmax
        do m = 0, l
            lm = lm + 1
            flmnnp(lm,1,1) = flmbas(lm)
        enddo
    enddo
!      recursion on  n  for the subsequent recursion on l
!      notice that recursion on n  runs in steps of two 2
!      whereas the storage index runs in steps of one, thus:
!      n = n0 + 2  (i-1)
    r2 = r * r
    lmaxaux = lmax
    ra = rac
    aux = r2 + ra*ra
    dosza = ra+ra
    do i = 2, la+1
        lmaxaux = lmaxaux - 1
        flmnnp(1,i,1) = aux * flmnnp(1,i-1,1)  - dosza * r2 * ri(3) * flmnnp(2,i-1,1)
        lm = 1
        do l = 1, lmaxaux
            do m = 0, l-1
                lm = lm + 1
                flmnnp(lm,i,1) = aux * flmnnp(lm,i-1,1)  - dosza * ( re(l-m) * ri(l+l-1) * flmnnp(ind(l-1)+m+1,i-1,1) &
                        + r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+m+1,i-1,1) )
            enddo
            lm = lm + 1
            flmnnp(lm,i,1) = aux * flmnnp(lm,i-1,1)  - dosza * r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+l+1,i-1,1)
        enddo
    enddo

!      recursion on  np  for the subsequent recursion on lp
!      notice that recursion on np  runs in steps of two 2
!      whereas the storage index runs in steps of one, thus:
!      np = np0 + 2  (ip-1)
    rb = rbc
    aux = r2 + rb*rb
    doszb = zb+zb
    lmaxbux = lmax+1
    do ip = 2, lb+1
    lmaxbux = lmaxbux - 1
    lmaxaux = lmaxbux
    do i = 1, la+1
        lmaxaux = lmaxaux-1
!           l = 0  m = 0
        flmnnp(1,i,ip) = aux * flmnnp(1,i,ip-1) - doszb * r2 * ri(3) * flmnnp(2,i,ip-1) &
                - xb * 0.6666666666666667d0 * r2 * flmnnp(3,i,ip-1)

        if (lmaxaux .ne. 0) then
!              l = 1  m = 0
            flmnnp(2,i,ip) = aux * flmnnp(2,i,ip-1) - doszb * ( flmnnp(1,i,ip-1) + r2 * .4d0 * flmnnp(4,i,ip-1) ) &
                    - xb * 1.2d0 * r2 * flmnnp(5,i,ip-1)

!              l = 1  m = 1
            flmnnp(3,i,ip) = aux * flmnnp(3,i,ip-1) - doszb * r2 * .6d0 * flmnnp(5,i,ip-1) - xb * ( 2.d0 * (flmnnp(1,i,ip-1) &
                    - r2 * .2d0 * flmnnp(4,i,ip-1)) + 2.4d0 * r2 * flmnnp(6,i,ip-1) )
            lm = 3
            do l = 2, lmaxaux
                lm = lm + 1
!                 m = 0
                flmnnp(lm,i,ip) = aux * flmnnp(lm,i,ip-1) - doszb * ( re(l) * ri(l+l-1) * flmnnp(ind(l-1)+1,i,ip-1) &
                        + r2 * re(l+1) * ri(l+l+3) * flmnnp(ind(l+1)+1,i,ip-1) )  - xb * (- re(l) * re(l-1) * ri(l+l-1) &
                        * flmnnp(ind(l-1)+2,i,ip-1) + re(l+1) * re(l+2) * ri(l+l+3) * r2 * flmnnp(ind(l+1)+2,i,ip-1) )
                dlt = 2.d0
                do m = 1, l-1
                    lm = lm + 1
                    flmnnp(lm,i,ip) = aux * flmnnp(lm,i,ip-1) - doszb * ( re(l-m) * ri(l+l-1) &
                    * flmnnp(ind(l-1)+m+1,i,ip-1) + r2 * re(l+m+1) * ri(l+l+3) * flmnnp(ind(l+1)+m+1,i,ip-1) ) &
                    - xb * ( dlt * (ri(l+l-1) * flmnnp(ind(l-1)+m,i,ip-1) - r2 * ri(l+l+3) * flmnnp(ind(l+1)+m,i,ip-1) ) &
                    - re(l-m) * re(l-m-1) * ri(l+l-1) * flmnnp(ind(l-1)+m+2,i,ip-1) &
                    + re(l+m+1) * re(l+m+2) * ri(l+l+3) * r2 * flmnnp(ind(l+1)+m+2,i,ip-1) )
                    dlt = 1.d0
                enddo
                lm = lm + 1
                flmnnp(lm,i,ip) = aux * flmnnp(lm,i,ip-1) - doszb * ( r2 * re(l+l+1) * ri(l+l+3) &
                        * flmnnp(ind(l+1)+l+1,i,ip-1))- xb * ( ri(l+l-1) * flmnnp(ind(l-1)+l,i,ip-1) - r2 * ri(l+l+3)  &
                        * flmnnp(ind(l+1)+l,i,ip-1) + re(l+l+1) * re(l+l+2) * ri(l+l+3) * r2 * flmnnp(ind(l+1)+l+2,i,ip-1) )
            enddo   ! end of do on l
        endif   ! end of  if (lmaxaux .ne. 0)
    enddo   ! end of do on i
    enddo   ! end of do on ip
    call recurflm(ldim, lmax, la, lb, r2, rac, xb, zb, flmnnp)
    return
    end subroutine flmGTO
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
!          lm = l(l+1)+m      lm = 0, 1, 2, ... (lmaxexp+1)2-1
!
  subroutine armonicosij(lmax, x, z, zlma)
    USE Zernike_Jacobi_STO_MPI_D
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
