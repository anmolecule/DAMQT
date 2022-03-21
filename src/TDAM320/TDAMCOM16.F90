!  Copyright 2008-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
!> @file TDAMCOM15.F90
!> @author Rafa
!
!> @details Error-------------------------------
  subroutine error(ierr, msg)
    USE DAM320_T
    implicit none
    integer(KINT) :: ierr, ierr2
    character(*) :: msg
    write(6,"(a)") msg
    write(6,"('Error code = ', i4)") ierr
    stop
  end
!
!   ***************************************************************
!> @details Calculates the derivatives
  subroutine derivzlm(lmax, idimzlm, zlma, zlmadx, zlmady, zlmadz)
    USE DAM320_T
    USE DAM320_CONST_T
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmadx(idimzlm), zlmady(idimzlm), zlmadz(idimzlm)
!	Derivatives of the regular harmonics with respecto to the Cartesian coordinates
    zlmadx(1) = cero	! Derivatives of the S spherical harmonics
    zlmady(1) = cero
    zlmadz(1) = cero
    zlmadx(2) = cero	! Derivatives of the P spherical harmonics
    zlmadx(3) = cero
    zlmadx(4) = zlma(1)
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    zlmadx(5) = re(6) * zlma(2)	! Derivatives of the D spherical harmonics
    zlmadx(6) = cero
    zlmadx(7) = -zlma(4)
    zlmadx(8) = re(3) * zlma(3)
    zlmadx(9) = re(6) * zlma(4)
    zlmady(5) = re(6) * zlma(4)
    zlmady(6) = re(3) * zlma(3)
    zlmady(7) = -zlma(2)
    zlmady(8) = cero
    zlmady(9) = -re(6) * zlma(2)
    zlmadz(5) = cero
    zlmadz(6) = re(3) * zlma(2)
    zlmadz(7) = re(2) * zlma(3)
    zlmadz(8) = re(3) * zlma(4)
    zlmadz(9) = cero
    zlmadx(10) = re(15) * zlma(5)		! Derivatives of the F spherical harmonics
    zlmadx(11) = re(10) * zlma(6)
    zlmadx(12) = -umed * zlma(5)
    zlmadx(13) = -zlma(8)
    zlmadx(14) = re(6) * zlma(7) - umed * zlma(9)
    zlmadx(15) = re(10) * zlma(8)
    zlmadx(16) = re(15) * zlma(9)
    zlmady(10) = re(15) * zlma(9)
    zlmady(11) = re(10) * zlma(8)
    zlmady(12) = re(6) * zlma(7) + umed * zlma(9)
    zlmady(13) = -zlma(6)
    zlmady(14) = -umed * zlma(5)
    zlmady(15) = -re(10) * zlma(6)
    zlmady(16) = -re(15) * zlma(5)
    zlmadz(10) = cero
    zlmadz(11) = re(5) * zlma(5)
    zlmadz(12) = re(4) * zlma(6)
    zlmadz(13) = re(3) * zlma(7)
    zlmadz(14) = re(4) * zlma(8)
    zlmadz(15) = re(5) * zlma(9)
    zlmadz(16) = cero
    do l = 4, lmax		! Derivatives of the remaining spherical harmonics
        zlmadx(l*(l+1)+1) = - zlma((l-1)*l+2)
        zlmady(l*(l+1)+1) = - zlma((l-1)*l)
        zlmadz(l*(l+1)+1) = re(l) * zlma((l-1)*l+1)
        zlmadx(l*(l+1)+2) = umed * (re(l+1) * re(l) * zlma((l-1)*l+1) - zlma((l-1)*l+3))
        zlmadx(l*(l+1)) = umed * (- zlma((l-1)*l-1))
        zlmady(l*(l+1)+2) = -umed * (zlma((l-1)*l-1))
        zlmady(l*(l+1)) = umed * (re(l+1) * re(l) * zlma((l-1)*l+1) + zlma((l-1)*l+3))
        zlmadz(l*(l+1)+2) = re(l+1) * zlma((l-1)*l+2)
        zlmadz(l*(l+1)) = re(l+1) * zlma((l-1)*l)
        do m = 2, l-2
            zlmadx(l*(l+1)+m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l+m) - zlma((l-1)*l+m+2))
            zlmadx(l*(l+1)-m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l-m+2) - zlma((l-1)*l-m))
            zlmady(l*(l+1)+m+1) = -umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l-m+2) + zlma((l-1)*l-m))
            zlmady(l*(l+1)-m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l+m) + zlma((l-1)*l+m+2))
            zlmadz(l*(l+1)+m+1) = re(l+m) * zlma((l-1)*l+m+1)
            zlmadz(l*(l+1)-m+1) = re(l+m) * zlma((l-1)*l-m+1)
        enddo
        zlmadx(l*(l+2)) = re(l+l-1) * re(l-1) * zlma(l*l-1)
        zlmady(l*(l+2)) = -re(l+l-1) * re(l-1) * zlma(l*(l-2)+3)
        zlmadz(l*(l+2)) = re(l+l-1) * zlma(l*l)
        zlmadx(l*l+2) = re(l+l-1) * re(l-1) * zlma(l*(l-2)+3)
        zlmady(l*l+2) = re(l+l-1) * re(l-1) * zlma(l*l-1)
        zlmadz(l*l+2) = re(l+l-1) * zlma(l*(l-2)+2)
        zlmadx((l+1)*(l+1)) = re(l) * re(l+l-1) * zlma(l*l)
        zlmady((l+1)*(l+1)) = -re(l) * re(l+l-1) * zlma((l-2)*l+2)
        zlmadz((l+1)*(l+1)) = cero
        zlmadx(l*l+1) = re(l) * re(l+l-1) * zlma((l-2)*l+2)
        zlmady(l*l+1) = re(l) * re(l+l-1) * zlma(l*l)
        zlmadz(l*l+1) = cero
    enddo
    return
    end
	
!   ***************************************************************
!> @details deals with derivative in z direction
  subroutine dzlm2y(lmax, idimzlm, zlma, zlmady, zlmadz)
    USE DAM320_T
    USE DAM320_CONST_T
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmady(idimzlm), zlmadz(idimzlm)
    zlmady(1) = cero
    zlmadz(1) = cero
    zlmady(2) = zlma(1)
    zlmady(3) = cero
    zlmady(4) = cero
    zlmadz(2) = cero
    zlmadz(3) = zlma(1)
    zlmadz(4) = cero
    zlmady(5) = re(6) * zlma(4)
    zlmady(6) = re(3) * zlma(3)
    zlmady(7) = -zlma(2)
    zlmady(8) = cero
    zlmady(9) = -re(6) * zlma(2)
    zlmadz(5) = cero
    zlmadz(6) = re(3) * zlma(2)
    zlmadz(7) = re(2) * zlma(3)
    zlmadz(8) = re(3) * zlma(4)
    zlmadz(9) = cero
    zlmady(10) = re(15) * zlma(9)
    zlmady(11) = re(10) * zlma(8)
    zlmady(12) = re(6) * zlma(7) + umed * zlma(9)
    zlmady(13) = -zlma(6)
    zlmady(14) = -umed * zlma(5)
    zlmady(15) = -re(10) * zlma(6)
    zlmady(16) = -re(15) * zlma(5)
    zlmadz(10) = cero
    zlmadz(11) = re(5) * zlma(5)
    zlmadz(12) = re(4) * zlma(6)
    zlmadz(13) = re(3) * zlma(7)
    zlmadz(14) = re(4) * zlma(8)
    zlmadz(15) = re(5) * zlma(9)
    zlmadz(16) = cero
    do l = 4, lmax		! Derivatives of the remaining spherical harmonics
        zlmady(l*(l+1)+1) = - zlma((l-1)*l)
        zlmadz(l*(l+1)+1) = re(l) * zlma((l-1)*l+1)
        zlmady(l*(l+1)+2) = -umed * (zlma((l-1)*l-1))
        zlmady(l*(l+1)) = umed * (re(l+1) * re(l) * zlma((l-1)*l+1) + zlma((l-1)*l+3))
        zlmadz(l*(l+1)+2) = re(l+1) * zlma((l-1)*l+2)
        zlmadz(l*(l+1)) = re(l+1) * zlma((l-1)*l)
        do m = 2, l-2
            zlmady(l*(l+1)+m+1) = -umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l-m+2) + zlma((l-1)*l-m))
            zlmady(l*(l+1)-m+1) = umed * (re(l+m) * re(l+m-1) * zlma((l-1)*l+m) + zlma((l-1)*l+m+2))
            zlmadz(l*(l+1)+m+1) = re(l+m) * zlma((l-1)*l+m+1)
            zlmadz(l*(l+1)-m+1) = re(l+m) * zlma((l-1)*l-m+1)
        enddo
        zlmady(l*(l+2)) = -re(l+l-1) * re(l-1) * zlma(l*(l-2)+3)
        zlmadz(l*(l+2)) = re(l+l-1) * zlma(l*l)
        zlmady(l*l+2) = re(l+l-1) * re(l-1) * zlma(l*l-1)
        zlmadz(l*l+2) = re(l+l-1) * zlma(l*(l-2)+2)
        zlmady((l+1)*(l+1)) = -re(l) * re(l+l-1) * zlma((l-2)*l+2)
        zlmadz((l+1)*(l+1)) = cero
        zlmady(l*l+1) = re(l) * re(l+l-1) * zlma(l*l)
        zlmadz(l*l+1) = cero
    enddo
    return
    end
		
!   ***************************************************************
!> @details deals with derivative in z direction
  subroutine dzlm2z(lmax, idimzlm, zlma, zlmadz)
    USE DAM320_T
    USE DAM320_CONST_T
    implicit none
    integer(KINT) :: idimzlm, lmax, l, m
    real(KREAL) :: zlma(idimzlm), zlmadz(idimzlm)
    zlmadz(1) = cero
    do l = 1, lmax		! Derivatives of the remaining spherical harmonics
        zlmadz(l*(l+1)+1) = re(l) * zlma((l-1)*l+1)
        zlmadz((l+1)*(l+1)) = cero
        zlmadz(l*l+1) = cero
        do m = 1, l-1
            zlmadz(l*(l+1)+m+1) = re(l+m) * zlma((l-1)*l+m+1)
            zlmadz(l*(l+1)-m+1) = re(l+m) * zlma((l-1)*l-m+1)
        enddo
    enddo
    return
    end
!
!   ******************************************************************
!> @details Used by consta
  subroutine emes ( m1, m2, ms, md, ss, sd)
    USE DAM320_T
    USE DAM320_CONST_T, ONLY: uno, cero, umed
    implicit none
    integer(KINT) :: m1, m1a, m2, m2a, ms, md
    real(KREAL) :: s1, s2, s12, ss, sd
    s1 = sign(1,m1)
    s2 = sign(1,m2)
    s12 = s1 * s2
    m1a = abs(m1)
    m2a = abs(m2)
    ms = s12 * ( m1a + m2a )
    md = s12 * abs( m1a - m2a )
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
!> @details subroutine rotar
!
!>	this subroutine yields the rotation matrices rl(m',m;l) of reals spherical harmonics
!>	receives the trigonometric functions of Euler angles defining the rotation
!
!**********************************************************************
  subroutine rotar(lmax, cosal, sinal, cosbet, sinbet, cosga, singa)
    USE DAM320_T
    USE DAM320_DATA_T
    USE DAM320_CONST_T
    implicit none
    integer(KINT) :: l, lmax
    real(KREAL) :: cosag, cosal, cosamg, cosbet, cosga, sinag, sinal, sinamg, singa, sinbet, tgbet2
!	Initial matrices d0, r0, d1 and r1
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
!	the remaining matrices are calculated using symmetry and recurrence relations by means of the subroutine dlmn.
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
!> @details   subroutine dlmn
!
!>   this subroutine generates the matrices dl(m',m;l) for a fixed value
!>   of the orbital quantum number l, and it needs the dl(l-2;m',m) and 
!>   dl(l-1;m',m) matrices. this subroutine uses symmetry and recurrence
!>   relations. the matrices dl(m',m;l) are the rotation matrices for   
!>   complex spherical harmonics
!
!**********************************************************************
  subroutine dlmn(l, sinal, cosal, cosbet, tgbet2, singa, cosga)
    USE DAM320_T
    USE DAM320_DATA_T
    USE DAM320_CONST_T
    implicit none
    integer(KINT) :: iinf, isup, l, m, mp
    real(KREAL) :: al, al1, ali, aux, cosag, cosagm, cosal, cosaux, cosbet, cosga, cosmal, cosmga, cux, d1, d2
    real(KREAL) :: sgn, sinag, sinagm, sinal, singa, sinmal, sinmga, tal1, tgbet2
    iinf=1-l
    isup=-iinf
!	computation of the dl(m',m;l) matrix, mp is m' and m is m.
!	first row by recurrence: see equations 19 and 20 of reference (6)
    dl(l,l,l) = dl(isup,isup,l-1) * (uno + cosbet) * umed
    dl(l,-l,l) = dl(isup,-isup,l-1) * (uno - cosbet) * umed
    do m = isup, iinf, -1
        dl(l,m,l) = -tgbet2 * root(l+m+1) * rooti(l-m) * dl(l,m+1,l)
    enddo
!	the rows of the upper quarter triangle of the dl(m',m;l) matrix see equation 21 of reference (6)
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
!	the remaining elements of the dl(m',m;l) matrix are calculated using the corresponding symmetry relations:
!		reflection ---> ((-1)**(m-m')) dl(m,m';l) = dl(m',m;l), m'<=m
!		inversion ---> ((-1)**(m-m')) dl(-m',-m;l) = dl(m',m;l)
!	reflection
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
!	inversion
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
!	computation of the rotation matrices rl(m',m;l) for real spherical harmonics using the matrices dl(m',m;l) 
!	for complex spherical harmonics: see equations 10 to 18 of reference (6)
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

