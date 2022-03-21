!  Copyright 2016, Rafael Lopez
! 
!  This file is part of Zernike_Jacobi_2016 package.
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
!
!*******************************************************
!
   subroutine subckplm( kmax, kdim, lmax, ldim, ckplm)
    implicit none

!	coeficients ckplm for computing the radial factors of two-center STO distributions from
!	the radial factors of the translation of individual STO
!	This subroutine was first generated with notebook integrales3C-7.nb and then modified by hand
!	Prepared for 0 <= l <= 22
!
!	OJO! OJO! OJO! OJO! OJO! OJO! OJO! OJO!
!
!	The ckplm defined herein differ from those appearing in the article of J. Comp. Chem. (2005) 26, 846-855
!	(Translation of STO Charge Distributions) in a factor (-1)**m
!
!	OJO! OJO! OJO! OJO! OJO! OJO! OJO! OJO!

    integer(4) :: ip, k, kdim, kmax, knt, l, lmax, ldim, m
    real(8) :: ckplm(0:kdim,ldim), den(0:lmax)
    real(8) :: rk, r2k
    ckplm = 0.d0
    do k = 0, kmax
        rk = k
        r2k = rk+rk

!    ckplm para l = 0
        den(0) = 1.d0/((r2k+1.D0))
        ckplm(k,1) = 1.D0 * den(0)
!    ckplm para l = 1
        den(0) = 1.d0/((r2k+1.D0)*(r2k+3.D0))
        den(1) = 1.d0/((r2k-1.D0)*(r2k+1.D0))
        ckplm(k,2) = 3.D0*(rk+1.D0) * den(0)
        ckplm(k,3) = 3.D0*rk * den(1)
        ckplm(k,4) = -3.D0 * den(0)
        ckplm(k,5) = 3.D0 * den(1)
!    ckplm para l = 2
        den(0) = 1.d0/((r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0))
        den(1) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0))
        den(2) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0))
        ckplm(k,6) = 7.5D0*(rk+1.D0)*(rk+2.D0) * den(0)
        ckplm(k,7) = 5.D0*(rk+1.D0)*rk * den(1)
        ckplm(k,8) = 7.5D0*(rk-1.D0)*rk * den(2)
        ckplm(k,9) = -5.D0*(rk+2.D0) * den(0)
        ckplm(k,10) = 5.D0 * den(1)
        ckplm(k,11) = 5.D0*(rk-1.D0) * den(2)
        ckplm(k,12) = 1.25D0 * den(0)
        ckplm(k,13) = -2.5D0 * den(1)
        ckplm(k,14) = 1.25D0 * den(2)
!    ckplm para l = 3
        den(0) = 1.d0/((r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0))
        den(1) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0))
        den(2) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0))
        den(3) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k-5.D0))
        ckplm(k,15) = 17.5D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0) * den(0)
        ckplm(k,16) = 10.5D0*(rk+1.D0)*(rk+2.D0)*rk * den(1)
        ckplm(k,17) = 10.5D0*(rk-1.D0)*(rk+1.D0)*rk * den(2)
        ckplm(k,18) = 17.5D0*(rk-1.D0)*(rk-2.D0)*rk * den(3)
        ckplm(k,19) = -8.75D0*(rk+2.D0)*(rk+3.D0) * den(0)
        ckplm(k,20) = -1.75D0*(rk+2.D0)*(rk-5.D0) * den(1)
        ckplm(k,21) = 1.75D0*(rk-1.D0)*(rk+6.D0) * den(2)
        ckplm(k,22) = 8.75D0*(rk-1.D0)*(rk-2.D0) * den(3)
        ckplm(k,23) = 1.75D0*(rk+3.D0) * den(0)
        ckplm(k,24) = -1.75D0*(rk+4.D0) * den(1)
        ckplm(k,25) = -1.75D0*(rk-3.D0) * den(2)
        ckplm(k,26) = 1.75D0*(rk-2.D0) * den(3)
        ckplm(k,27) = -0.291666666666667D0 * den(0)
        ckplm(k,28) = 0.875D0 * den(1)
        ckplm(k,29) = -0.875D0 * den(2)
        ckplm(k,30) = 0.291666666666667D0 * den(3)
!    ckplm para l = 4
        den(0) = 1.d0/((r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0))
        den(2) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0))
        den(3) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0))
        den(4) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0))
        ckplm(k,31) = 39.375D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0) * den(0)
        ckplm(k,32) = 22.5D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*rk * den(1)
        ckplm(k,33) = 20.25D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*rk * den(2)
        ckplm(k,34) = 22.5D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*rk * den(3)
        ckplm(k,35) = 39.375D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*rk * den(4)
        ckplm(k,36) = -15.75D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0) * den(0)
        ckplm(k,37) = -2.25D0*(rk+2.D0)*(rk+3.D0)*(-7.D0+2.D0*rk) * den(1)
        ckplm(k,38) = 20.25D0*(rk-1.D0)*(rk+2.D0) * den(2)
        ckplm(k,39) = 2.25D0*(rk-1.D0)*(rk-2.D0)*(9.D0+2.D0*rk) * den(3)
        ckplm(k,40) = 15.75D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0) * den(4)
        ckplm(k,41) = 2.625D0*(rk+3.D0)*(rk+4.D0) * den(0)
        ckplm(k,42) = -1.5D0*(rk+3.D0)*(rk+7.D0) * den(1)
        ckplm(k,43) = -2.25D0*(rk*rk+rk-9.D0) * den(2)
        ckplm(k,44) = -1.5D0*(rk-2.D0)*(rk-6.D0) * den(3)
        ckplm(k,45) = 2.625D0*(rk-2.D0)*(rk-3.D0) * den(4)
        ckplm(k,46) = -0.375D0*(rk+4.D0) * den(0)
        ckplm(k,47) = 0.375D0*(9.D0+2.D0*rk) * den(1)
        ckplm(k,48) = -3.375D0 * den(2)
        ckplm(k,49) = -0.375D0*(-7.D0+2.D0*rk) * den(3)
        ckplm(k,50) = 0.375D0*(rk-3.D0) * den(4)
        ckplm(k,51) = 0.046875D0 * den(0)
        ckplm(k,52) = -0.1875D0 * den(1)
        ckplm(k,53) = 0.28125D0 * den(2)
        ckplm(k,54) = -0.1875D0 * den(3)
        ckplm(k,55) = 0.046875D0 * den(4)
!    ckplm para l = 5
        den(0) = 1.d0/((r2k+11.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0))
        den(3) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0))
        den(4) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0))
        den(5) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,56) = 86.625D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0) * den(0)
        ckplm(k,57) = 48.125D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*rk * den(1)
        ckplm(k,58) = 41.25D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*rk * den(2)
        ckplm(k,59) = 41.25D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*rk * den(3)
        ckplm(k,60) = 48.125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*rk * den(4)
        ckplm(k,61) = 86.625D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*rk * den(5)
        ckplm(k,62) = -28.875D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0) * den(0)
        ckplm(k,63) = -9.625D0*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0) * den(1)
        ckplm(k,64) = -2.75D0*(rk-14.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0) * den(2)
        ckplm(k,65) = 2.75D0*(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0) * den(3)
        ckplm(k,66) = 9.625D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk+4.D0) * den(4)
        ckplm(k,67) = 28.875D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0) * den(5)
        ckplm(k,68) = 4.125D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0) * den(0)
        ckplm(k,69) = -1.375D0*(rk+12.D0)*(rk+3.D0)*(rk+4.D0) * den(1)
        ckplm(k,70) = -2.75D0*(rk+3.D0)*(-13.D0+(rk+3.D0)*rk) * den(2)
        ckplm(k,71) = -2.75D0*(rk-2.D0)*(-15.D0+(rk-1.D0)*rk) * den(3)
        ckplm(k,72) = -1.375D0*(rk-11.D0)*(rk-2.D0)*(rk-3.D0) * den(4)
        ckplm(k,73) = 4.125D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0) * den(5)
        ckplm(k,74) = -0.515625D0*(rk+4.D0)*(rk+5.D0) * den(0)
        ckplm(k,75) = 0.0572916666666667D0*(rk+4.D0)*(81.D0+13.D0*rk) * den(1)
        ckplm(k,76) = 0.34375D0*(-48.D0+(rk-7.D0)*rk) * den(2)
        ckplm(k,77) = -0.34375D0*(-40.D0+(rk+9.D0)*rk) * den(3)
        ckplm(k,78) = -0.0572916666666667D0*(rk-3.D0)*(-68.D0+13.D0*rk) * den(4)
        ckplm(k,79) = 0.515625D0*(rk-3.D0)*(rk-4.D0) * den(5)
        ckplm(k,80) = 0.0572916666666667D0*(rk+5.D0) * den(0)
        ckplm(k,81) = -0.0572916666666667D0*(16.D0+3.D0*rk) * den(1)
        ckplm(k,82) = 0.114583333333333D0*(rk+11.D0) * den(2)
        ckplm(k,83) = 0.114583333333333D0*(rk-10.D0) * den(3)
        ckplm(k,84) = -0.0572916666666667D0*(-13.D0+3.D0*rk) * den(4)
        ckplm(k,85) = 0.0572916666666667D0*(rk-4.D0) * den(5)
        ckplm(k,86) = -0.00572916666666667D0 * den(0)
        ckplm(k,87) = 0.0286458333333333D0 * den(1)
        ckplm(k,88) = -0.0572916666666667D0 * den(2)
        ckplm(k,89) = 0.0572916666666667D0 * den(3)
        ckplm(k,90) = -0.0286458333333333D0 * den(4)
        ckplm(k,91) = 0.00572916666666667D0 * den(5)
!    ckplm para l = 6
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0))
        den(4) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0))
        den(5) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,92) = 187.6875D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0) * den(0)
        ckplm(k,93) = 102.375D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*rk * den(1)
        ckplm(k,94) = 85.3125D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*rk * den(2)
        ckplm(k,95) = 81.25D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*rk * den(3)
        ckplm(k,96) = 85.3125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*rk * den(4)
        ckplm(k,97) = 102.375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*rk * den(5)
        ckplm(k,98) = 187.6875D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*rk * den(6)
        ckplm(k,99) = -53.625D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0) * den(0)
        ckplm(k,100) = -4.875D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(-11.D0+4.D0*rk) * den(1)
        ckplm(k,101) = -8.125D0*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk-9.D0) * den(2)
        ckplm(k,102) = 81.25D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0) * den(3)
        ckplm(k,103) = 8.125D0*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0) * den(4)
        ckplm(k,104) = 4.875D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(15.D0+4.D0*rk) * den(5)
        ckplm(k,105) = 53.625D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0) * den(6)
        ckplm(k,106) = 6.703125D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0) * den(0)
        ckplm(k,107) = -1.21875D0*(rk+22.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0) * den(1)
        ckplm(k,108) = -0.203125D0*(rk+3.D0)*(rk+4.D0)*(-306.D0+rk*(91.D0+17.D0*rk)) * den(2)
        ckplm(k,109) = -4.0625D0*(rk-2.D0)*(rk+3.D0)*(rk-4.D0)*(rk+5.D0) * den(3)
        ckplm(k,110) = -0.203125D0*(rk-2.D0)*(rk-3.D0)*(-380.D0+rk*(-57.D0+17.D0*rk)) * den(4)
        ckplm(k,111) = -1.21875D0*(rk-21.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0) * den(5)
        ckplm(k,112) = 6.703125D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0) * den(6)
        ckplm(k,113) = -0.744791666666667D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0) * den(0)
        ckplm(k,114) = 0.203125D0*(rk+4.D0)*(rk+5.D0)*(33.D0+4.D0*rk) * den(1)
        ckplm(k,115) = 0.609375D0*(rk+4.D0)*(-43.D0+(rk-2.D0)*rk) * den(2)
        ckplm(k,116) = -2.03125D0*(-40.D0+3.D0*(rk+1.D0)*rk) * den(3)
        ckplm(k,117) = -0.609375D0*(rk-3.D0)*(-40.D0+(rk+4.D0)*rk) * den(4)
        ckplm(k,118) = -0.203125D0*(rk-3.D0)*(rk-4.D0)*(-29.D0+4.D0*rk) * den(5)
        ckplm(k,119) = 0.744791666666667D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0) * den(6)
        ckplm(k,120) = 0.0744791666666667D0*(rk+5.D0)*(rk+6.D0) * den(0)
        ckplm(k,121) = -0.0135416666666667D0*(rk+5.D0)*(88.D0+13.D0*rk) * den(1)
        ckplm(k,122) = 0.0338541666666667D0*(216.D0+(rk+47.D0)*rk) * den(2)
        ckplm(k,123) = 0.135416666666667D0*(rk*rk+rk-50.D0) * den(3)
        ckplm(k,124) = 0.0338541666666667D0*(170.D0+(rk-45.D0)*rk) * den(4)
        ckplm(k,125) = -0.0135416666666667D0*(rk-4.D0)*(-75.D0+13.D0*rk) * den(5)
        ckplm(k,126) = 0.0744791666666667D0*(rk-4.D0)*(rk-5.D0) * den(6)
        ckplm(k,127) = -0.00677083333333334D0*(rk+6.D0) * den(0)
        ckplm(k,128) = 0.00677083333333334D0*(25.D0+4.D0*rk) * den(1)
        ckplm(k,129) = -0.0338541666666667D0*(rk+9.D0) * den(2)
        ckplm(k,130) = 0.338541666666667D0 * den(3)
        ckplm(k,131) = 0.0338541666666667D0*(rk-8.D0) * den(4)
        ckplm(k,132) = -0.00677083333333334D0*(-21.D0+4.D0*rk) * den(5)
        ckplm(k,133) = 0.00677083333333334D0*(rk-5.D0) * den(6)
        ckplm(k,134) = 0.00056423611111111D0 * den(0)
        ckplm(k,135) = -0.00338541666666667D0 * den(1)
        ckplm(k,136) = 0.00846354166666666D0 * den(2)
        ckplm(k,137) = -0.0112847222222222D0 * den(3)
        ckplm(k,138) = 0.00846354166666666D0 * den(4)
        ckplm(k,139) = -0.00338541666666667D0 * den(5)
        ckplm(k,140) = 0.00056423611111111D0 * den(6)
!    ckplm para l = 7
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(3) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(4) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0))
        den(5) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        ckplm(k,141) = 402.1875D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0) * den(0)
        ckplm(k,142) = 216.5625D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*rk &
                * den(1)
        ckplm(k,143) = 177.1875D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*rk &
                * den(2)
        ckplm(k,144) = 164.0625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*rk &
                * den(3)
        ckplm(k,145) = 164.0625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*rk &
                * den(4)
        ckplm(k,146) = 177.1875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*rk &
                * den(5)
        ckplm(k,147) = 216.5625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*rk &
                * den(6)
        ckplm(k,148) = 402.1875D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*rk &
                * den(7)
        ckplm(k,149) = -100.546875D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) &
                * den(0)
        ckplm(k,150) = -7.734375D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(-13.D0+5.D0&
                *rk) * den(1)
        ckplm(k,151) = -6.328125D0*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(-22.D0+3.D0&
                *rk) * den(2)
        ckplm(k,152) = -5.859375D0*(rk-1.D0)*(rk-27.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0) &
                * den(3)
        ckplm(k,153) = 5.859375D0*(rk-1.D0)*(rk+28.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0) &
                * den(4)
        ckplm(k,154) = 6.328125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(25.D0+3.D0*rk) &
                * den(5)
        ckplm(k,155) = 7.734375D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(18.D0+5.D0*rk) &
                * den(6)
        ckplm(k,156) = 100.546875D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0) &
                * den(7)
        ckplm(k,157) = 11.171875D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) * den(0)
        ckplm(k,158) = -0.859375D0*(rk+3.D0)*(rk+4.D0)*(rk+52.D0)*(rk+5.D0)*(rk+6.D0) * den(1)
        ckplm(k,159) = -0.234375D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(-462.D0+rk*(157.D0+19.D0*rk)) &
                * den(2)
        ckplm(k,160) = -5.859375D0*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(-26.D0+(rk+3.D0)*rk) * den(3)
        ckplm(k,161) = -5.859375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(-28.D0+(rk-1.D0)*rk) * den(4)
        ckplm(k,162) = -0.234375D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(-600.D0+rk*(-119.D0+19.D0*rk)) &
                * den(5)
        ckplm(k,163) = -0.859375D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-51.D0)*(rk-5.D0) * den(6)
        ckplm(k,164) = 11.171875D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0) * den(7)
        ckplm(k,165) = -1.1171875D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) * den(0)
        ckplm(k,166) = 0.0859375D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(117.D0+11.D0*rk) * den(1)
        ckplm(k,167) = 0.0703125D0*(rk+4.D0)*(rk+5.D0)*(-594.D0+rk*(9.D0+13.D0*rk)) * den(2)
        ckplm(k,168) = 0.1171875D0*(rk+4.D0)*(1250.D0+3.D0*rk*(-73.D0+(rk-22.D0)*rk)) * den(3)
        ckplm(k,169) = -0.1171875D0*(rk-3.D0)*(-1400.D0+3.D0*(rk-1.D0)*(rk+26.D0)*rk) * den(4)
        ckplm(k,170) = -0.0703125D0*(rk-3.D0)*(rk-4.D0)*(-590.D0+rk*(17.D0+13.D0*rk)) * den(5)
        ckplm(k,171) = -0.0859375D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(-106.D0+11.D0*rk) * den(6)
        ckplm(k,172) = 1.1171875D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0) * den(7)
        ckplm(k,173) = 0.1015625D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) * den(0)
        ckplm(k,174) = -0.0078125D0*(rk+5.D0)*(rk+6.D0)*(208.D0+25.D0*rk) * den(1)
        ckplm(k,175) = -0.0234375D0*(rk+5.D0)*(-456.D0+(rk-65.D0)*rk) * den(2)
        ckplm(k,176) = 0.1171875D0*(-400.D0+rk*(-38.D0+(rk+13.D0)*rk)) * den(3)
        ckplm(k,177) = 0.1171875D0*(350.D0+rk*(-61.D0+(rk-10.D0)*rk)) * den(4)
        ckplm(k,178) = -0.0234375D0*(rk-4.D0)*(-390.D0+(rk+67.D0)*rk) * den(5)
        ckplm(k,179) = 0.0078125D0*(rk-4.D0)*(rk-5.D0)*(183.D0-25.D0*rk) * den(6)
        ckplm(k,180) = 0.1015625D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0) * den(7)
        ckplm(k,181) = -0.00846354166666666D0*(rk+6.D0)*(rk+7.D0) * den(0)
        ckplm(k,182) = 0.00065104166666667D0*(rk+6.D0)*(325.D0+43.D0*rk) * den(1)
        ckplm(k,183) = 0.001953125D0*(-1050.D0-rk*(239.D0+11.D0*rk)) * den(2)
        ckplm(k,184) = -0.0162760416666667D0*(-138.D0+(rk-11.D0)*rk) * den(3)
        ckplm(k,185) = 0.0162760416666667D0*(-126.D0+(rk+13.D0)*rk) * den(4)
        ckplm(k,186) = 0.001953125D0*(822.D0+rk*(-217.D0+11.D0*rk)) * den(5)
        ckplm(k,187) = 0.00065104166666667D0*(rk-5.D0)*(282.D0-43.D0*rk) * den(6)
        ckplm(k,188) = 0.00846354166666666D0*(rk-5.D0)*(rk-6.D0) * den(7)
        ckplm(k,189) = 0.00065104166666667D0*(rk+7.D0) * den(0)
        ckplm(k,190) = 0.00065104166666667D0*(-36.D0-5.D0*rk) * den(1)
        ckplm(k,191) = 0.005859375D0*(rk+9.D0) * den(2)
        ckplm(k,192) = -0.00325520833333333D0*(rk+22.D0) * den(3)
        ckplm(k,193) = -0.00325520833333333D0*(rk-21.D0) * den(4)
        ckplm(k,194) = 0.005859375D0*(rk-8.D0) * den(5)
        ckplm(k,195) = 0.00065104166666667D0*(31.D0-5.D0*rk) * den(6)
        ckplm(k,196) = 0.00065104166666667D0*(rk-6.D0) * den(7)
        ckplm(k,197) = -0.00004650297619048D0 * den(0)
        ckplm(k,198) = 0.00032552083333333D0 * den(1)
        ckplm(k,199) = -0.0009765625D0 * den(2)
        ckplm(k,200) = 0.00162760416666667D0 * den(3)
        ckplm(k,201) = -0.00162760416666667D0 * den(4)
        ckplm(k,202) = 0.0009765625D0 * den(5)
        ckplm(k,203) = -0.00032552083333333D0 * den(6)
        ckplm(k,204) = 0.00004650297619048D0 * den(7)
!    ckplm para l = 8
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+1.D0)*(r2k+3.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+3.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k-9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k-7.D0)*(r2k-9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,205) = 854.6484375D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0) * den(0)
        ckplm(k,206) = 455.8125D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*rk * den(1)
        ckplm(k,207) = 368.15625D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*rk * den(2)
        ckplm(k,208) = 334.6875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*rk * den(3)
        ckplm(k,209) = 325.390625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*rk * den(4)
        ckplm(k,210) = 334.6875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*rk * den(5)
        ckplm(k,211) = 368.15625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*rk * den(6)
        ckplm(k,212) = 455.8125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*rk * den(7)
        ckplm(k,213) = 854.6484375D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*rk * den(8)
        ckplm(k,214) = -189.921875D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0) * den(0)
        ckplm(k,215) = -37.984375D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(-5.D0+2.D0*rk) * den(1)
        ckplm(k,216) = -20.453125D0*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(-13.D0+2.D0*rk) * den(2)
        ckplm(k,217) = -9.296875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(-33.D0+2.D0*rk) * den(3)
        ckplm(k,218) = 325.390625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0) &
                * den(4)
        ckplm(k,219) = 9.296875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(35.D0+2.D0*rk) * den(5)
        ckplm(k,220) = 20.453125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(15.D0+2.D0*rk) * den(6)
        ckplm(k,221) = 37.984375D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(7.D0+2.D0*rk) * den(7)
        ckplm(k,222) = 189.921875D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0) * den(8)
        ckplm(k,223) = 18.9921875D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) &
                * den(0)
        ckplm(k,224) = -75.96875D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) * den(1)
        ckplm(k,225) = -2.921875D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(-65.D0+2.D0*(rk+12.D0)&
                *rk) * den(2)
        ckplm(k,226) = -0.53125D0*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(-528.D0+rk*(83.D0+16.D0&
                *rk)) * den(3)
        ckplm(k,227) = -9.296875D0*(rk*rk+rk-35.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0) * den(4)
        ckplm(k,228) = -0.53125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(-595.D0+rk*(-51.D0+16.D0&
                *rk)) * den(5)
        ckplm(k,229) = -2.921875D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(-87.D0+2.D0*(rk-10.D0)&
                *rk) * den(6)
        ckplm(k,230) = 75.96875D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0) * den(7)
        ckplm(k,231) = 18.9921875D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0) &
                * den(8)
        ckplm(k,232) = -1.7265625D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) * den(0)
        ckplm(k,233) = 0.575520833333333D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(27.D0+2.D0*rk) &
                * den(1)
        ckplm(k,234) = 0.1328125D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(-507.D0+rk*(29.D0+10.D0*rk)) &
                * den(2)
        ckplm(k,235) = 0.3984375D0*(rk+4.D0)*(rk+5.D0)*(642.D0+rk*(-153.D0+rk*(-23.D0+2.D0*rk))) &
                * den(3)
        ckplm(k,236) = -4.6484375D0*(rk-3.D0)*(rk+4.D0)*(-70.D0+3.D0*(rk+1.D0)*rk) * den(4)
        ckplm(k,237) = -0.3984375D0*(rk-3.D0)*(rk-4.D0)*(-770.D0+rk*(-101.D0+rk*(29.D0+2.D0*rk))) &
                * den(5)
        ckplm(k,238) = -0.1328125D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(-526.D0+rk*(-9.D0+10.D0*rk)) &
                * den(6)
        ckplm(k,239) = -0.575520833333333D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(-25.D0+2.D0*rk) &
                * den(7)
        ckplm(k,240) = 1.7265625D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0) * den(8)
        ckplm(k,241) = 0.143880208333333D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) * den(0)
        ckplm(k,242) = -0.230208333333333D0*(rk+10.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) * den(1)
        ckplm(k,243) = -0.0265625D0*(rk-26.D0)*(rk+5.D0)*(rk+6.D0)*(23.D0+3.D0*rk) * den(2)
        ckplm(k,244) = 0.0885416666666667D0*(rk+5.D0)*(-876.D0+rk*(-10.D0+(rk+27.D0)*rk)) * den(3)
        ckplm(k,245) = 0.154947916666667D0*(2100.D0+(rk*rk+rk-122.D0)*(rk+1.D0)*rk) * den(4)
        ckplm(k,246) = 0.0885416666666667D0*(rk-4.D0)*(840.D0+rk*(-61.D0+(rk-24.D0)*rk)) * den(5)
        ckplm(k,247) = -0.0265625D0*(rk+27.D0)*(rk-4.D0)*(rk-5.D0)*(-20.D0+3.D0*rk) * den(6)
        ckplm(k,248) = -0.230208333333333D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-9.D0) * den(7)
        ckplm(k,249) = 0.143880208333333D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0) * den(8)
        ckplm(k,250) = -0.0110677083333333D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) * den(0)
        ckplm(k,251) = 0.00221354166666667D0*(rk+6.D0)*(rk+7.D0)*(125.D0+14.D0*rk) * den(1)
        ckplm(k,252) = -0.006640625D0*(rk+6.D0)*(425.D0+rk*(73.D0+2.D0*rk)) * den(2)
        ckplm(k,253) = -0.0110677083333333D0*(-1590.D0+rk*(-281.D0+rk*(9.D0+2.D0*rk))) * den(3)
        ckplm(k,254) = 0.387369791666667D0*(rk-6.D0)*(rk+7.D0) * den(4)
        ckplm(k,255) = 0.0110677083333333D0*(1302.D0+rk*(-293.D0+rk*(-3.D0+2.D0*rk))) * den(5)
        ckplm(k,256) = 0.006640625D0*(rk-5.D0)*(354.D0+rk*(-69.D0+2.D0*rk)) * den(6)
        ckplm(k,257) = -0.00221354166666667D0*(rk-5.D0)*(rk-6.D0)*(-111.D0+14.D0*rk) * den(7)
        ckplm(k,258) = 0.0110677083333333D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0) * den(8)
        ckplm(k,259) = 0.0007905505952381D0*(rk+7.D0)*(rk+8.D0) * den(0)
        ckplm(k,260) = -0.00021081349206349D0*(rk+7.D0)*(135.D0+16.D0*rk) * den(1)
        ckplm(k,261) = 0.00221354166666667D0*(183.D0+2.D0*(rk+20.D0)*rk) * den(2)
        ckplm(k,262) = -0.06640625D0*(rk+8.D0) * den(3)
        ckplm(k,263) = -0.00368923611111111D0*(rk*rk+rk-147.D0) * den(4)
        ckplm(k,264) = 0.06640625D0*(rk-7.D0) * den(5)
        ckplm(k,265) = 0.00221354166666667D0*(145.D0+2.D0*(rk-18.D0)*rk) * den(6)
        ckplm(k,266) = -0.00021081349206349D0*(rk-6.D0)*(-119.D0+16.D0*rk) * den(7)
        ckplm(k,267) = 0.0007905505952381D0*(rk-6.D0)*(rk-7.D0) * den(8)
        ckplm(k,268) = -0.00005270337301587D0*(rk+8.D0) * den(0)
        ckplm(k,269) = 0.00005270337301587D0*(49.D0+6.D0*rk) * den(1)
        ckplm(k,270) = -0.00036892361111111D0*(19.D0+2.D0*rk) * den(2)
        ckplm(k,271) = 0.00036892361111111D0*(31.D0+2.D0*rk) * den(3)
        ckplm(k,272) = -0.0129123263888889D0 * den(4)
        ckplm(k,273) = -0.00036892361111111D0*(-29.D0+2.D0*rk) * den(5)
        ckplm(k,274) = 0.00036892361111111D0*(-17.D0+2.D0*rk) * den(6)
        ckplm(k,275) = -0.00005270337301587D0*(-43.D0+6.D0*rk) * den(7)
        ckplm(k,276) = 0.00005270337301587D0*(rk-7.D0) * den(8)
        ckplm(k,277) = 3.29396081349206D-6 * den(0)
        ckplm(k,278) = -0.00002635168650794D0 * den(1)
        ckplm(k,279) = 0.00009223090277778D0 * den(2)
        ckplm(k,280) = -0.00018446180555556D0 * den(3)
        ckplm(k,281) = 0.00023057725694444D0 * den(4)
        ckplm(k,282) = -0.00018446180555556D0 * den(5)
        ckplm(k,283) = 0.00009223090277778D0 * den(6)
        ckplm(k,284) = -0.00002635168650794D0 * den(7)
        ckplm(k,285) = 3.29396081349206D-6 * den(8)
!    ckplm para l = 9
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,286) = 1804.2578125D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,287) = 955.1953125D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*rk * den(1)
        ckplm(k,288) = 764.15625D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*rk * den(2)
        ckplm(k,289) = 685.78125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*rk * den(3)
        ckplm(k,290) = 654.609375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*rk * den(4)
        ckplm(k,291) = 654.609375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*rk * den(5)
        ckplm(k,292) = 685.78125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*rk * den(6)
        ckplm(k,293) = 764.15625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*rk * den(7)
        ckplm(k,294) = 955.1953125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*rk * den(8)
        ckplm(k,295) = 1804.2578125D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*rk * den(9)
        ckplm(k,296) = -360.8515625D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,297) = -21.2265625D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(-17.D0+7.D0*rk) * den(1)
        ckplm(k,298) = -84.90625D0*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0) * den(2)
        ckplm(k,299) = -45.71875D0*(rk-13.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0) * den(3)
        ckplm(k,300) = -14.546875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-44.D0)&
                *(rk+4.D0)*(rk+5.D0) * den(4)
        ckplm(k,301) = 14.546875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+45.D0)&
                *(rk-4.D0)*(rk+4.D0) * den(5)
        ckplm(k,302) = 45.71875D0*(rk+14.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0) * den(6)
        ckplm(k,303) = 84.90625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+7.D0) * den(7)
        ckplm(k,304) = 21.2265625D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(24.D0+7.D0*rk) * den(8)
        ckplm(k,305) = 360.8515625D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0) * den(9)
        ckplm(k,306) = 32.8046875D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,307) = 1.9296875D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk-68.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0) * den(1)
        ckplm(k,308) = -3.859375D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(-87.D0+2.D0&
                *(rk+17.D0)*rk) * den(2)
        ckplm(k,309) = -2.078125D0*(6.D0*rk*rk+46.D0*rk-247.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0) * den(3)
        ckplm(k,310) = -14.546875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(-43.D0+(rk+3.D0)*rk) * den(4)
        ckplm(k,311) = -14.546875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(-45.D0+(rk-1.D0)*rk) * den(5)
        ckplm(k,312) = -2.078125D0*(6.D0*rk*rk-34.D0*rk-287.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0) * den(6)
        ckplm(k,313) = -3.859375D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(-119.D0+2.D0&
                *(rk-15.D0)*rk) * den(7)
        ckplm(k,314) = 1.9296875D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk+69.D0)*(rk-6.D0)&
                *(rk-7.D0) * den(8)
        ckplm(k,315) = 32.8046875D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0) * den(9)
        ckplm(k,316) = -2.73372395833333D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,317) = 1.447265625D0*(rk+17.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) &
                * den(1)
        ckplm(k,318) = 1.9296875D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(-57.D0+(rk+5.D0)*rk) &
                * den(2)
        ckplm(k,319) = 0.0494791666666667D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(8931.D0+rk*(-2516.D0+rk&
                *(-204.D0+29.D0*rk))) * den(3)
        ckplm(k,320) = 0.51953125D0*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(1176.D0+rk*(-124.D0+(rk-39.D0)&
                *rk)) * den(4)
        ckplm(k,321) = -0.51953125D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(-1260.D0+(rk-1.D0)*(rk+43.D0)&
                *rk) * den(5)
        ckplm(k,322) = -0.0494791666666667D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(-11214.D0+rk*(-2021.D0+rk&
                *(291.D0+29.D0*rk))) * den(6)
        ckplm(k,323) = -1.9296875D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(-61.D0+(rk-3.D0)*rk) &
                * den(7)
        ckplm(k,324) = -1.447265625D0*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0) &
                * den(8)
        ckplm(k,325) = 2.73372395833333D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0) * den(9)
        ckplm(k,326) = 0.210286458333333D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,327) = -0.0123697916666667D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(272.D0+23.D0&
                *rk) * den(1)
        ckplm(k,328) = -0.1484375D0*(rk-18.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+9.D0) * den(2)
        ckplm(k,329) = 0.0494791666666667D0*(rk+5.D0)*(rk+6.D0)*(-2556.D0+rk*(116.D0+(rk+69.D0)*rk)) &
                * den(3)
        ckplm(k,330) = 0.173177083333333D0*(rk+5.D0)*(3444.D0+rk*(-490.D0+rk*(-145.D0+(rk+10.D0)&
                *rk))) * den(4)
        ckplm(k,331) = 0.173177083333333D0*(rk-4.D0)*(3780.D0+(rk-1.D0)*rk*(-174.D0+(rk-5.D0)*rk)) &
                * den(5)
        ckplm(k,332) = 0.0494791666666667D0*(rk-4.D0)*(rk-5.D0)*(2604.D0+rk*(-19.D0+(rk-66.D0)*rk)) &
                * den(6)
        ckplm(k,333) = -0.1484375D0*(rk+19.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-8.D0) * den(7)
        ckplm(k,334) = -0.0123697916666667D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(-249.D0+23.D0&
                *rk) * den(8)
        ckplm(k,335) = 0.210286458333333D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0) * den(9)
        ckplm(k,336) = -0.0150204613095238D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,337) = 0.00088355654761905D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(425.D0+41.D0*rk) * den(1)
        ckplm(k,338) = -0.00212053571428571D0*(rk+6.D0)*(rk+7.D0)*(1875.D0+rk*(247.D0+3.D0*rk)) &
                * den(2)
        ckplm(k,339) = -0.0247395833333333D0*(rk+15.D0)*(rk+6.D0)*(-73.D0+(rk-3.D0)*rk) * den(3)
        ckplm(k,340) = -0.0123697916666667D0*(11760.D0+rk*(824.D0+rk*(-379.D0+(rk-26.D0)*rk))) &
                * den(4)
        ckplm(k,341) = 0.0123697916666667D0*(10584.D0+rk*(-1500.D0+rk*(-295.D0+(rk+30.D0)*rk))) &
                * den(5)
        ckplm(k,342) = 0.0247395833333333D0*(rk-14.D0)*(rk-5.D0)*(-69.D0+(rk+5.D0)*rk) * den(6)
        ckplm(k,343) = 0.00212053571428571D0*(rk-5.D0)*(rk-6.D0)*(1631.D0+rk*(-241.D0+3.D0*rk)) &
                * den(7)
        ckplm(k,344) = -0.00088355654761905D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(-384.D0+41.D0*rk) &
                * den(8)
        ckplm(k,345) = 0.0150204613095238D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0) * den(9)
        ckplm(k,346) = 0.00100136408730159D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,347) = -0.00053013392857143D0*(rk+7.D0)*(rk+8.D0)*(68.D0+7.D0*rk) * den(1)
        ckplm(k,348) = 0.00035342261904762D0*(rk+7.D0)*(1509.D0+2.D0*rk*(133.D0+5.D0*rk)) * den(2)
        ckplm(k,349) = 0.00082465277777778D0*(-5622.D0+rk*(-1193.D0+2.D0*(rk-21.D0)*rk)) * den(3)
        ckplm(k,350) = -0.00247395833333333D0*(-1911.D0+rk*(-124.D0+(rk+24.D0)*rk)) * den(4)
        ckplm(k,351) = -0.00247395833333333D0*(1764.D0+rk*(-169.D0+(rk-21.D0)*rk)) * den(5)
        ckplm(k,352) = 0.00082465277777778D0*(4473.D0+rk*(-1103.D0+2.D0*(rk+24.D0)*rk)) * den(6)
        ckplm(k,353) = 0.00035342261904762D0*(rk-6.D0)*(1253.D0+2.D0*rk*(-123.D0+5.D0*rk)) * den(7)
        ckplm(k,354) = -0.00053013392857143D0*(rk-6.D0)*(rk-7.D0)*(-61.D0+7.D0*rk) * den(8)
        ckplm(k,355) = 0.00100136408730159D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0) * den(9)
        ckplm(k,356) = -0.00006258525545635D0*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,357) = 3.68148561507937D-6*(rk+8.D0)*(833.D0+89.D0*rk) * den(1)
        ckplm(k,358) = -0.00001472594246032D0*(4116.D0+rk*(841.D0+41.D0*rk)) * den(2)
        ckplm(k,359) = 0.00010308159722222D0*(916.D0+rk*(131.D0+3.D0*rk)) * den(3)
        ckplm(k,360) = 0.00036078559027778D0*(-304.D0+(rk-15.D0)*rk) * den(4)
        ckplm(k,361) = -0.00036078559027778D0*(-288.D0+(rk+17.D0)*rk) * den(5)
        ckplm(k,362) = -0.00010308159722222D0*(788.D0+rk*(-125.D0+3.D0*rk)) * den(6)
        ckplm(k,363) = 0.00001472594246032D0*(3316.D0+rk*(-759.D0+41.D0*rk)) * den(7)
        ckplm(k,364) = -3.68148561507937D-6*(rk-7.D0)*(-744.D0+89.D0*rk) * den(8)
        ckplm(k,365) = 0.00006258525545635D0*(rk-7.D0)*(rk-8.D0) * den(9)
        ckplm(k,366) = 3.68148561507937D-6*(rk+9.D0) * den(0)
        ckplm(k,367) = -3.68148561507937D-6*(64.D0+7.D0*rk) * den(1)
        ckplm(k,368) = 0.00001472594246032D0*(51.D0+5.D0*rk) * den(2)
        ckplm(k,369) = -0.00010308159722222D0*(rk+14.D0) * den(3)
        ckplm(k,370) = 0.00005154079861111D0*(rk+37.D0) * den(4)
        ckplm(k,371) = 0.00005154079861111D0*(rk-36.D0) * den(5)
        ckplm(k,372) = -0.00010308159722222D0*(rk-13.D0) * den(6)
        ckplm(k,373) = 0.00001472594246032D0*(-46.D0+5.D0*rk) * den(7)
        ckplm(k,374) = -3.68148561507937D-6*(-57.D0+7.D0*rk) * den(8)
        ckplm(k,375) = 3.68148561507937D-6*(rk-8.D0) * den(9)
        ckplm(k,376) = -2.0452697861552D-7 * den(0)
        ckplm(k,377) = 1.84074280753968D-6 * den(1)
        ckplm(k,378) = -7.36297123015873D-6 * den(2)
        ckplm(k,379) = 0.0000171802662037D0 * den(3)
        ckplm(k,380) = -0.00002577039930556D0 * den(4)
        ckplm(k,381) = 0.00002577039930556D0 * den(5)
        ckplm(k,382) = -0.0000171802662037D0 * den(6)
        ckplm(k,383) = 7.36297123015873D-6 * den(7)
        ckplm(k,384) = -1.84074280753968D-6 * den(8)
        ckplm(k,385) = 2.0452697861552D-7 * den(9)
!    ckplm para l = 10
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,386) = 3788.94140625D0*(rk+10.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,387) = 1994.1796875D0*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,388) = 1583.61328125D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*rk * den(2)
        ckplm(k,389) = 1407.65625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*rk * den(3)
        ckplm(k,390) = 1326.4453125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*rk * den(4)
        ckplm(k,391) = 1302.328125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*rk * den(5)
        ckplm(k,392) = 1326.4453125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*rk * den(6)
        ckplm(k,393) = 1407.65625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*rk * den(7)
        ckplm(k,394) = 1583.61328125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*rk * den(8)
        ckplm(k,395) = 1994.1796875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*rk * den(9)
        ckplm(k,396) = 3788.94140625D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,397) = -688.8984375D0*(rk+10.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,398) = -36.2578125D0*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-19.D0+8.D0*rk) * den(1)
        ckplm(k,399) = -57.5859375D0*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(-17.D0+3.D0*rk) * den(2)
        ckplm(k,400) = -25.59375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(-45.D0+4.D0*rk) * den(3)
        ckplm(k,401) = -48.234375D0*(rk-1.D0)*(rk-26.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0) * den(4)
        ckplm(k,402) = 1302.328125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0) * den(5)
        ckplm(k,403) = 48.234375D0*(rk-1.D0)*(rk+27.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0) * den(6)
        ckplm(k,404) = 25.59375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(49.D0+4.D0*rk) * den(7)
        ckplm(k,405) = 57.5859375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(20.D0+3.D0*rk) * den(8)
        ckplm(k,406) = 36.2578125D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(27.D0+8.D0*rk) * den(9)
        ckplm(k,407) = 688.8984375D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,408) = 57.408203125D0*(rk+10.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,409) = 6.04296875D0*(rk-38.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,410) = -0.533203125D0*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(-1122.D0+rk*(457.D0+19.D0*rk)) * den(2)
        ckplm(k,411) = -1.421875D0*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(-660.D0+rk*(137.D0+13.D0*rk)) * den(3)
        ckplm(k,412) = -1.33984375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(-884.D0+rk*(87.D0+17.D0*rk)) * den(4)
        ckplm(k,413) = -24.1171875D0*(rk*rk+rk-54.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0) * den(5)
        ckplm(k,414) = -1.33984375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(-954.D0+rk*(-53.D0+17.D0*rk)) * den(6)
        ckplm(k,415) = -1.421875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(-784.D0+rk*(-111.D0+13.D0*rk)) * den(7)
        ckplm(k,416) = -0.533203125D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(-1560.D0+rk*(-419.D0+19.D0*rk)) * den(8)
        ckplm(k,417) = 6.04296875D0*(rk-2.D0)*(rk+39.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0) * den(9)
        ckplm(k,418) = 57.408203125D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,419) = -4.416015625D0*(rk+10.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,420) = 0.232421875D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(171.D0+8.D0*rk) * den(1)
        ckplm(k,421) = 0.123046875D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-1479.D0+rk&
                *(164.D0+23.D0*rk)) * den(2)
        ckplm(k,422) = 0.0546875D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(13935.D0+rk*(-4367.D0+rk&
                *(-177.D0+44.D0*rk))) * den(3)
        ckplm(k,423) = 1.33984375D0*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(834.D0+rk&
                *(-124.D0+(rk-21.D0)*rk)) * den(4)
        ckplm(k,424) = -36.17578125D0*(rk*rk+rk-36.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0) &
                * den(5)
        ckplm(k,425) = -1.33984375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(-936.D0+rk&
                *(-79.D0+(rk+24.D0)*rk)) * den(6)
        ckplm(k,426) = -0.0546875D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(-18081.D0+rk&
                *(-3881.D0+rk*(309.D0+44.D0*rk))) * den(7)
        ckplm(k,427) = -0.123046875D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(-1620.D0+rk&
                *(-118.D0+23.D0*rk)) * den(8)
        ckplm(k,428) = -0.232421875D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(-163.D0+8.D0*rk) * den(9)
        ckplm(k,429) = 4.416015625D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(10)
        ckplm(k,430) = 0.3154296875D0*(rk+10.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,431) = -0.033203125D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(152.D0+11.D0&
                *rk) * den(1)
        ckplm(k,432) = -0.0029296875D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-12648.D0+rk&
                *(-355.D0+83.D0*rk)) * den(2)
        ckplm(k,433) = -0.0078125D0*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(-5268.D0+(rk-598.D0)&
                *rk) * den(3)
        ckplm(k,434) = 0.013671875D0*(rk+5.D0)*(rk+6.D0)*(77028.D0+rk*(-15406.D0+rk*(-2209.D0+rk&
                *(274.D0+13.D0*rk)))) * den(4)
        ckplm(k,435) = 0.24609375D0*(rk-4.D0)*(rk+5.D0)*(5292.D0+(rk*rk+rk-198.D0)*(rk+1.D0)*rk) &
                * den(5)
        ckplm(k,436) = 0.013671875D0*(rk-4.D0)*(rk-5.D0)*(89964.D0+rk*(10218.D0+rk*(-2953.D0+rk&
                *(-222.D0+13.D0*rk)))) * den(6)
        ckplm(k,437) = -0.0078125D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(-4669.D0+(rk+600.D0)&
                *rk) * den(7)
        ckplm(k,438) = -0.0029296875D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(-12210.D0+rk&
                *(521.D0+83.D0*rk)) * den(8)
        ckplm(k,439) = -0.033203125D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(-141.D0+11.D0*rk) * den(9)
        ckplm(k,440) = 0.3154296875D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(10)
        ckplm(k,441) = -0.0210286458333333D0*(rk+10.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,442) = 0.00553385416666667D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(95.D0+8.D0*rk) &
                * den(1)
        ckplm(k,443) = 0.0009765625D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-5865.D0+(rk-596.D0)*rk) * den(2)
        ckplm(k,444) = -0.00651041666666667D0*(rk+6.D0)*(rk+7.D0)*(-6393.D0+rk*(-373.D0+rk&
                *(81.D0+4.D0*rk))) * den(3)
        ckplm(k,445) = -0.0227864583333333D0*(rk+6.D0)*(10878.D0+rk*(-22.D0+rk*(-325.D0+(rk-8.D0)&
                *rk))) * den(4)
        ckplm(k,446) = 0.123046875D0*(10584.D0+5.D0*(rk*rk+rk-100.D0)*(rk+1.D0)*rk) * den(5)
        ckplm(k,447) = 0.0227864583333333D0*(rk-5.D0)*(10584.D0+rk*(-600.D0+rk*(-295.D0+(rk+12.D0)&
                *rk))) * den(6)
        ckplm(k,448) = 0.00651041666666667D0*(rk-5.D0)*(rk-6.D0)*(5943.D0+rk*(-523.D0+rk*(-69.D0+4.D0&
                *rk))) * den(7)
        ckplm(k,449) = 0.0009765625D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(5268.D0-(rk+598.D0)*rk) * den(8)
        ckplm(k,450) = -0.00553385416666667D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-87.D0+8.D0&
                *rk) * den(9)
        ckplm(k,451) = 0.0210286458333333D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(10)
        ckplm(k,452) = 0.00131429036458333D0*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,453) = -0.00013834635416667D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(342.D0+31.D0*rk) * den(1)
        ckplm(k,454) = 0.00003662109375D0*(rk+7.D0)*(rk+8.D0)*(19686.D0+rk*(2845.D0+79.D0*rk)) &
                * den(2)
        ckplm(k,455) = 0.00003255208333333D0*(rk+7.D0)*(-206880.D0+rk*(-32234.D0+rk*(-159.D0+83.D0&
                *rk))) * den(3)
        ckplm(k,456) = -0.00005696614583333D0*(-818496.D0+rk*(-118998.D0+rk*(9533.D0+rk&
                *(1542.D0+19.D0*rk)))) * den(4)
        ckplm(k,457) = -0.003076171875D0*(14112.D0+(rk*rk+rk-338.D0)*(rk+1.D0)*rk) * den(5)
        ckplm(k,458) = -0.00005696614583333D0*(-691488.D0+rk*(133514.D0+rk*(5021.D0+rk&
                *(-1466.D0+19.D0*rk)))) * den(6)
        ckplm(k,459) = 0.00003255208333333D0*(rk-6.D0)*(174888.D0+rk*(-31667.D0+rk*(408.D0+83.D0&
                *rk))) * den(7)
        ckplm(k,460) = 0.00003662109375D0*(rk-6.D0)*(rk-7.D0)*(16920.D0+rk*(-2687.D0+79.D0*rk)) &
                * den(8)
        ckplm(k,461) = -0.00013834635416667D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-311.D0+31.D0*rk) &
                * den(9)
        ckplm(k,462) = 0.00131429036458333D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,463) = -0.00007731119791667D0*(rk+10.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,464) = 4.06901041666667D-6*(rk+8.D0)*(rk+9.D0)*(931.D0+88.D0*rk) * den(1)
        ckplm(k,465) = -0.00001220703125D0*(rk+8.D0)*(6321.D0+rk*(1084.D0+43.D0*rk)) * den(2)
        ckplm(k,466) = 0.00001627604166667D0*(55860.D0+rk*(12503.D0+rk*(723.D0+4.D0*rk))) * den(3)
        ckplm(k,467) = 0.00039876302083333D0*(-2616.D0+rk*(-304.D0+(rk+9.D0)*rk)) * den(4)
        ckplm(k,468) = -0.0107666015625D0*(rk*rk+rk-96.D0) * den(5)
        ckplm(k,469) = -0.00039876302083333D0*(2304.D0+rk*(-319.D0+(rk-6.D0)*rk)) * den(6)
        ckplm(k,470) = 0.00001627604166667D0*(44076.D0+rk*(-11069.D0+(711.D0-4.D0*rk)*rk)) * den(7)
        ckplm(k,471) = 0.00001220703125D0*(rk-7.D0)*(5280.D0+rk*(-998.D0+43.D0*rk)) * den(8)
        ckplm(k,472) = 4.06901041666667D-6*(rk-7.D0)*(rk-8.D0)*(843.D0-88.D0*rk) * den(9)
        ckplm(k,473) = 0.00007731119791667D0*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,474) = 4.29506655092593D-6*(rk+10.D0)*(rk+9.D0) * den(0)
        ckplm(k,475) = -4.52112268518519D-7*(rk+9.D0)*(608.D0+59.D0*rk) * den(1)
        ckplm(k,476) = 2.03450520833333D-6*(3552.D0+rk*(673.D0+31.D0*rk)) * den(2)
        ckplm(k,477) = 5.42534722222222D-6*(-2430.D0-rk*(359.D0+11.D0*rk)) * den(3)
        ckplm(k,478) = -9.49435763888889D-6*(-1842.D0+(rk-149.D0)*rk) * den(4)
        ckplm(k,479) = 0.00005696614583333D0*(rk*rk+rk-324.D0) * den(5)
        ckplm(k,480) = -9.49435763888889D-6*(-1692.D0+(rk+151.D0)*rk) * den(6)
        ckplm(k,481) = 5.42534722222222D-6*(-2082.D0+(337.D0-11.D0*rk)*rk) * den(7)
        ckplm(k,482) = 2.03450520833333D-6*(2910.D0+rk*(-611.D0+31.D0*rk)) * den(8)
        ckplm(k,483) = 4.52112268518519D-7*(rk-8.D0)*(549.D0-59.D0*rk) * den(9)
        ckplm(k,484) = 4.29506655092593D-6*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,485) = 2.26056134259259D-7*(-10.D0-rk) * den(0)
        ckplm(k,486) = 2.26056134259259D-7*(81.D0+8.D0*rk) * den(1)
        ckplm(k,487) = 6.103515625D-6*(-11.D0-rk) * den(2)
        ckplm(k,488) = 2.71267361111111D-6*(55.D0+4.D0*rk) * den(3)
        ckplm(k,489) = -9.49435763888889D-6*(rk+24.D0) * den(4)
        ckplm(k,490) = 0.00025634765625D0 * den(5)
        ckplm(k,491) = 9.49435763888889D-6*(rk-23.D0) * den(6)
        ckplm(k,492) = 2.71267361111111D-6*(51.D0-4.D0*rk) * den(7)
        ckplm(k,493) = 6.103515625D-6*(rk-10.D0) * den(8)
        ckplm(k,494) = 2.26056134259259D-7*(73.D0-8.D0*rk) * den(9)
        ckplm(k,495) = 2.26056134259259D-7*(rk-9.D0) * den(10)
        ckplm(k,496) = 1.1302806712963D-8 * den(0)
        ckplm(k,497) = -1.1302806712963D-7 * den(1)
        ckplm(k,498) = 5.08626302083333D-7 * den(2)
        ckplm(k,499) = -1.35633680555556D-6 * den(3)
        ckplm(k,500) = 2.37358940972222D-6 * den(4)
        ckplm(k,501) = -2.84830729166667D-6 * den(5)
        ckplm(k,502) = 2.37358940972222D-6 * den(6)
        ckplm(k,503) = -1.35633680555556D-6 * den(7)
        ckplm(k,504) = 5.08626302083333D-7 * den(8)
        ckplm(k,505) = -1.1302806712963D-7 * den(9)
        ckplm(k,506) = 1.1302806712963D-8 * den(10)
!    ckplm para l = 11
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,507) = 7922.33203125D0*(rk+10.D0)*(rk+11.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,508) = 4149.79296875D0*(rk+10.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,509) = 3276.15234375D0*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,510) = 2890.72265625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*rk * den(3)
        ckplm(k,511) = 2698.0078125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*rk * den(4)
        ckplm(k,512) = 2614.9921875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*rk * den(5)
        ckplm(k,513) = 2614.9921875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*rk * den(6)
        ckplm(k,514) = 2698.0078125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*rk * den(7)
        ckplm(k,515) = 2890.72265625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*rk * den(8)
        ckplm(k,516) = 3276.15234375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*rk * den(9)
        ckplm(k,517) = 4149.79296875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,518) = 7922.33203125D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,519) = -1320.388671875D0*(rk+10.D0)*(rk+11.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,520) = -188.626953125D0*(rk+10.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-7.D0+3.D0*rk) * den(1)
        ckplm(k,521) = -49.638671875D0*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-38.D0+7.D0*rk) * den(2)
        ckplm(k,522) = -43.798828125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-51.D0+5.D0*rk) * den(3)
        ckplm(k,523) = -122.63671875D0*(rk-1.D0)*(rk-20.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0) * den(4)
        ckplm(k,524) = -39.62109375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk-65.D0)*(rk+6.D0) * den(5)
        ckplm(k,525) = 39.62109375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+66.D0) * den(6)
        ckplm(k,526) = 122.63671875D0*(rk-1.D0)*(rk+21.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0) * den(7)
        ckplm(k,527) = 43.798828125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(56.D0+5.D0*rk) * den(8)
        ckplm(k,528) = 49.638671875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(45.D0+7.D0*rk) * den(9)
        ckplm(k,529) = 188.626953125D0*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(10.D0+3.D0*rk) * den(10)
        ckplm(k,530) = 1320.388671875D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,531) = 101.568359375D0*(rk+10.D0)*(rk+11.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,532) = 14.509765625D0*(rk+10.D0)*(rk-28.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,533) = -0.763671875D0*(rk+37.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-38.D0+17.D0*rk) * den(2)
        ckplm(k,534) = -0.673828125D0*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(-2550.D0+rk*(571.D0+41.D0*rk)) * den(3)
        ckplm(k,535) = -1.88671875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(236.D0+19.D0*rk) * den(4)
        ckplm(k,536) = -39.62109375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(-64.D0+(rk+3.D0)*rk) * den(5)
        ckplm(k,537) = -39.62109375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(-66.D0+(rk-1.D0)*rk) * den(6)
        ckplm(k,538) = -1.88671875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(-217.D0+19.D0*rk) * den(7)
        ckplm(k,539) = -0.673828125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(-3080.D0+rk*(-489.D0+41.D0*rk)) * den(8)
        ckplm(k,540) = -0.763671875D0*(rk-2.D0)*(rk-36.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(55.D0+17.D0*rk) * den(9)
        ckplm(k,541) = 14.509765625D0*(rk+29.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,542) = 101.568359375D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,543) = -7.2548828125D0*(rk+10.D0)*(rk+11.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,544) = 2.41829427083333D0*(rk+10.D0)*(rk+27.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,545) = 0.3818359375D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-798.D0+rk*(103.D0+11.D0*rk)) * den(2)
        ckplm(k,546) = 0.1123046875D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(11730.D0+rk&
                *(-3959.D0+rk*(-54.D0+35.D0*rk))) * den(3)
        ckplm(k,547) = 0.0149739583333333D0*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(134760.D0+rk*(-24302.D0+rk*(-2517.D0+179.D0*rk))) * den(4)
        ckplm(k,548) = 0.943359375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(2646.D0+rk&
                *(-187.D0+(rk-60.D0)*rk)) * den(5)
        ckplm(k,549) = -0.943359375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(-2772.D0+(rk-1.D0)*(rk+64.D0)*rk) * den(6)
        ckplm(k,550) = -0.0149739583333333D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(-156366.D0+rk*(-18731.D0+rk*(3054.D0+179.D0*rk))) * den(7)
        ckplm(k,551) = -0.1123046875D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(-15600.D0+rk*(-3746.D0+rk*(159.D0+35.D0*rk))) * den(8)
        ckplm(k,552) = -0.3818359375D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(-890.D0+rk*(-81.D0+11.D0*rk)) * den(9)
        ckplm(k,553) = -2.41829427083333D0*(rk-26.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,554) = 7.2548828125D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,555) = 0.483658854166667D0*(rk+10.D0)*(rk+11.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,556) = -0.483658854166667D0*(rk+10.D0)*(rk+16.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,557) = -0.3818359375D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-152.D0+(rk-1.D0)*rk) * den(2)
        ckplm(k,558) = -0.00748697916666667D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(44880.D0+rk&
                *(-5278.D0+rk*(-855.D0+13.D0*rk))) * den(3)
        ckplm(k,559) = 0.0149739583333333D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(121980.D0+rk*(-29522.D0+rk&
                *(-2135.D0+rk*(446.D0+11.D0*rk)))) * den(4)
        ckplm(k,560) = 0.314453125D0*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(7812.D0+rk*(-742.D0+rk&
                *(-229.D0+(rk+10.D0)*rk))) * den(5)
        ckplm(k,561) = 0.314453125D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(8316.D0+(rk-1.D0)*rk&
                *(-258.D0+(rk-5.D0)*rk)) * den(6)
        ckplm(k,562) = 0.0149739583333333D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(148932.D0+rk*(23958.D0+rk&
                *(-3407.D0+rk*(-402.D0+11.D0*rk)))) * den(7)
        ckplm(k,563) = -0.00748697916666667D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(-49290.D0+rk&
                *(-3529.D0+rk*(894.D0+13.D0*rk))) * den(8)
        ckplm(k,564) = -0.3818359375D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(-150.D0+(rk+3.D0)*rk) * den(9)
        ckplm(k,565) = -0.483658854166667D0*(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,566) = 0.483658854166667D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,567) = -0.0302286783854167D0*(rk+10.D0)*(rk+11.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,568) = 0.00431838262648809D0*(rk+10.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(175.D0+13.D0*rk) * den(1)
        ckplm(k,569) = 0.00340924944196429D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2470.D0+rk&
                *(-193.D0+3.D0*rk)) * den(2)
        ckplm(k,570) = -0.00033424014136905D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-192270.D0+rk&
                *(-4075.D0+rk*(2466.D0+79.D0*rk))) * den(3)
        ckplm(k,571) = -0.00467936197916667D0*(rk+16.D0)*(rk+6.D0)*(rk+7.D0)*(5514.D0+rk*(-643.D0+rk&
                *(-102.D0+7.D0*rk))) * den(4)
        ckplm(k,572) = -0.0028076171875D0*(rk+6.D0)*(-860832.D0+5.D0*rk*(20602.D0+rk*(6099.D0+rk&
                *(-443.D0+(rk-51.D0)*rk)))) * den(5)
        ckplm(k,573) = 0.0028076171875D0*(rk-5.D0)*(931392.D0+5.D0*(rk-1.D0)*rk*(-7284.D0+rk&
                *(-172.D0+(rk+57.D0)*rk))) * den(6)
        ckplm(k,574) = 0.00467936197916667D0*(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(-6048.D0+rk*(-418.D0+rk&
                *(123.D0+7.D0*rk))) * den(7)
        ckplm(k,575) = 0.00033424014136905D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(185808.D0+rk*(-8770.D0+rk&
                *(-2229.D0+79.D0*rk))) * den(8)
        ckplm(k,576) = -0.00340924944196429D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-2274.D0+rk&
                *(199.D0+3.D0*rk)) * den(9)
        ckplm(k,577) = -0.00431838262648809D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-162.D0+13.D0*rk) * den(10)
        ckplm(k,578) = 0.0302286783854167D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(11)
        ckplm(k,579) = 0.00177815755208333D0*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,580) = -0.00008467416914683D0*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(756.D0+61.D0&
                *rk) * den(1)
        ckplm(k,581) = 0.00006684802827381D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(14934.D0+rk&
                *(1789.D0+35.D0*rk)) * den(2)
        ckplm(k,582) = 0.00006684802827381D0*(rk+7.D0)*(rk+8.D0)*(-147300.D0+rk*(-16748.D0+rk&
                *(405.D0+53.D0*rk))) * den(3)
        ckplm(k,583) = 0.00031195746527778D0*(rk+7.D0)*(239040.D0+rk*(19998.D0+rk&
                *(-3625.D0+(rk-294.D0)*rk))) * den(4)
        ckplm(k,584) = -0.0028076171875D0*(169344.D0+rk*(9360.D0+rk*(-5066.D0+rk*(-289.D0+(rk+26.D0)&
                *rk)))) * den(5)
        ckplm(k,585) = -0.0028076171875D0*(-155232.D0+rk*(18526.D0+rk*(4053.D0+rk*(-383.D0+(rk-21.D0)&
                *rk)))) * den(6)
        ckplm(k,586) = 0.00031195746527778D0*(rk-6.D0)*(215712.D0+rk*(-26362.D0+rk&
                *(-2737.D0+(rk+298.D0)*rk))) * den(7)
        ckplm(k,587) = 0.00006684802827381D0*(rk-6.D0)*(rk-7.D0)*(130200.D0+rk*(-17399.D0+rk&
                *(-246.D0+53.D0*rk))) * den(8)
        ckplm(k,588) = 0.00006684802827381D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(13180.D0+rk&
                *(-1719.D0+35.D0*rk)) * den(9)
        ckplm(k,589) = -0.00008467416914683D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-695.D0+61.D0&
                *rk) * den(10)
        ckplm(k,590) = 0.00177815755208333D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(11)
        ckplm(k,591) = -0.0000987865306713D0*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,592) = 0.00001411236152447D0*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(343.D0+29.D0*rk) * den(1)
        ckplm(k,593) = -7.42755869708995D-7*(rk+8.D0)*(rk+9.D0)*(135926.D0+rk*(19839.D0+643.D0*rk)) &
                * den(2)
        ckplm(k,594) = -0.0000334240141369D0*(rk+8.D0)*(-37730.D0+rk*(-6707.D0+rk*(-254.D0+3.D0&
                *rk))) * den(3)
        ckplm(k,595) = 0.00003119574652778D0*(-359520.D0+rk*(-67686.D0+rk*(-719.D0+rk*(354.D0+11.D0&
                *rk)))) * den(4)
        ckplm(k,596) = 0.00021837022569444D0*(50976.D0+rk*(2734.D0+rk*(-757.D0+(rk-34.D0)*rk))) &
                * den(5)
        ckplm(k,597) = -0.00021837022569444D0*(47520.D0+rk*(-4142.D0+rk*(-649.D0+(rk+38.D0)*rk))) &
                * den(6)
        ckplm(k,598) = -0.00003119574652778D0*(-292896.D0+rk*(65230.D0+rk*(-1715.D0+rk*(-310.D0+11.D0&
                *rk)))) * den(7)
        ckplm(k,599) = 0.0000334240141369D0*(rk-7.D0)*(31280.D0+rk*(-6190.D0+rk*(263.D0+3.D0*rk))) &
                * den(8)
        ckplm(k,600) = 7.42755869708995D-7*(rk-7.D0)*(rk-8.D0)*(116730.D0+rk*(-18553.D0+643.D0*rk)) &
                * den(9)
        ckplm(k,601) = -0.00001411236152447D0*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-314.D0+29.D0*rk) &
                * den(10)
        ckplm(k,602) = 0.0000987865306713D0*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,603) = 5.19929108796296D-6*(rk+10.D0)*(rk+11.D0)*(rk+9.D0) * den(0)
        ckplm(k,604) = -7.42755869708995D-7*(rk+10.D0)*(rk+9.D0)*(448.D0+39.D0*rk) * den(1)
        ckplm(k,605) = 7.42755869708995D-7*(rk+9.D0)*(12064.D0+rk*(1971.D0+77.D0*rk)) * den(2)
        ckplm(k,606) = -0.0000334240141369D0*(4160.D0+(rk+30.D0)*(rk+31.D0)*rk) * den(3)
        ckplm(k,607) = -0.00003119574652778D0*(-5790.D0+rk*(-837.D0+(rk-14.D0)*rk)) * den(4)
        ckplm(k,608) = 0.00003119574652778D0*(-6264.D0+rk*(-286.D0+(rk+39.D0)*rk)) * den(5)
        ckplm(k,609) = 0.00003119574652778D0*(5940.D0+rk*(-361.D0+(rk-36.D0)*rk)) * den(6)
        ckplm(k,610) = -0.00003119574652778D0*(4968.D0+rk*(-806.D0+(rk+17.D0)*rk)) * den(7)
        ckplm(k,611) = -0.0000334240141369D0*(-3290.D0+rk*(811.D0+(rk-58.D0)*rk)) * den(8)
        ckplm(k,612) = 7.42755869708995D-7*(rk-8.D0)*(10170.D0+rk*(-1817.D0+77.D0*rk)) * den(9)
        ckplm(k,613) = -7.42755869708995D-7*(rk-8.D0)*(rk-9.D0)*(-409.D0+39.D0*rk) * den(10)
        ckplm(k,614) = 5.19929108796296D-6*(rk-10.D0)*(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,615) = -2.59964554398148D-7*(rk+10.D0)*(rk+11.D0) * den(0)
        ckplm(k,616) = 1.23792644951499D-8*(rk+10.D0)*(1701.D0+151.D0*rk) * den(1)
        ckplm(k,617) = -1.85688967427249D-7*(3834.D0+rk*(673.D0+29.D0*rk)) * den(2)
        ckplm(k,618) = 5.57066902281746D-7*(2690.D0+rk*(393.D0+13.D0*rk)) * den(3)
        ckplm(k,619) = -2.59964554398148D-6*(870.D0+(rk+85.D0)*rk) * den(4)
        ckplm(k,620) = -4.67936197916667D-6*(-570.D0+(rk-19.D0)*rk) * den(5)
        ckplm(k,621) = 4.67936197916667D-6*(-550.D0+(rk+21.D0)*rk) * den(6)
        ckplm(k,622) = 2.59964554398148D-6*(786.D0+(rk-83.D0)*rk) * den(7)
        ckplm(k,623) = -5.57066902281746D-7*(2310.D0+rk*(-367.D0+13.D0*rk)) * den(8)
        ckplm(k,624) = 1.85688967427249D-7*(3190.D0+rk*(-615.D0+29.D0*rk)) * den(9)
        ckplm(k,625) = -1.23792644951499D-8*(rk-9.D0)*(-1550.D0+151.D0*rk) * den(10)
        ckplm(k,626) = 2.59964554398148D-7*(rk-10.D0)*(rk-9.D0) * den(11)
        ckplm(k,627) = 1.23792644951499D-8*(rk+11.D0) * den(0)
        ckplm(k,628) = -1.23792644951499D-8*(100.D0+9.D0*rk) * den(1)
        ckplm(k,629) = 6.18963224757496D-8*(83.D0+7.D0*rk) * den(2)
        ckplm(k,630) = -9.28444837136243D-7*(rk+14.D0) * den(3)
        ckplm(k,631) = 3.71377934854497D-7*(61.D0+3.D0*rk) * den(4)
        ckplm(k,632) = -5.19929108796296D-7*(rk+56.D0) * den(5)
        ckplm(k,633) = -5.19929108796296D-7*(rk-55.D0) * den(6)
        ckplm(k,634) = 3.71377934854497D-7*(-58.D0+3.D0*rk) * den(7)
        ckplm(k,635) = -9.28444837136243D-7*(rk-13.D0) * den(8)
        ckplm(k,636) = 6.18963224757496D-8*(-76.D0+7.D0*rk) * den(9)
        ckplm(k,637) = -1.23792644951499D-8*(-91.D0+9.D0*rk) * den(10)
        ckplm(k,638) = 1.23792644951499D-8*(rk-10.D0) * den(11)
        ckplm(k,639) = -5.62693840688632D-10 * den(0)
        ckplm(k,640) = 6.18963224757496D-9 * den(1)
        ckplm(k,641) = -3.09481612378748D-8 * den(2)
        ckplm(k,642) = 9.28444837136244D-8 * den(3)
        ckplm(k,643) = -1.85688967427249D-7 * den(4)
        ckplm(k,644) = 2.59964554398148D-7 * den(5)
        ckplm(k,645) = -2.59964554398148D-7 * den(6)
        ckplm(k,646) = 1.85688967427249D-7 * den(7)
        ckplm(k,647) = -9.28444837136244D-8 * den(8)
        ckplm(k,648) = 3.09481612378748D-8 * den(9)
        ckplm(k,649) = -6.18963224757496D-9 * den(10)
        ckplm(k,650) = 5.62693840688632D-10 * den(11)
!    ckplm para l = 12
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,651) = 16504.8583984375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+1.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,652) = 8611.23046875D0*(rk+10.D0)*(rk+11.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,653) = 6765.966796875D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,654) = 5935.05859375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,655) = 5498.6572265625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*rk * den(4)
        ckplm(k,656) = 5278.7109375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*rk * den(5)
        ckplm(k,657) = 5211.03515625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*rk * den(6)
        ckplm(k,658) = 5278.7109375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*rk * den(7)
        ckplm(k,659) = 5498.6572265625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*rk * den(8)
        ckplm(k,660) = 5935.05859375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*rk * den(9)
        ckplm(k,661) = 6765.966796875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,662) = 8611.23046875D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,663) = 16504.8583984375D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,664) = -2539.208984375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,665) = -110.400390625D0*(rk+10.D0)*(rk+11.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-23.D0+10.D0*rk) * den(1)
        ckplm(k,666) = -173.486328125D0*(rk+10.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-21.D0+4.D0*rk) * den(2)
        ckplm(k,667) = -228.271484375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-19.D0+2.D0*rk) * den(3)
        ckplm(k,668) = -281.982421875D0*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) * den(4)
        ckplm(k,669) = -67.67578125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(-75.D0+2.D0*rk) * den(5)
        ckplm(k,670) = 5211.03515625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0) * den(6)
        ckplm(k,671) = 67.67578125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(77.D0+2.D0*rk) * den(7)
        ckplm(k,672) = 281.982421875D0*(rk+18.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0) * den(8)
        ckplm(k,673) = 228.271484375D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(21.D0+2.D0*rk) * den(9)
        ckplm(k,674) = 173.486328125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(25.D0+4.D0*rk) * den(10)
        ckplm(k,675) = 110.400390625D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(33.D0+10.D0*rk) * den(11)
        ckplm(k,676) = 2539.208984375D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,677) = 181.3720703125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,678) = 31.54296875D0*(rk+10.D0)*(rk+11.D0)*(rk-23.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,679) = -15.771484375D0*(rk+10.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-123.D0+(rk+53.D0)*rk) * den(2)
        ckplm(k,680) = -41.50390625D0*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-76.D0+(rk+18.D0)*rk) * den(3)
        ckplm(k,681) = -1.8310546875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(-2278.D0+rk*(309.D0+31.D0*rk)) * den(4)
        ckplm(k,682) = -1.7578125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(-2775.D0+rk*(188.D0+37.D0*rk)) * den(5)
        ckplm(k,683) = -67.67578125D0*(rk*rk+rk-77.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0) * den(6)
        ckplm(k,684) = -1.7578125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(-2926.D0+rk*(-114.D0+37.D0*rk)) * den(7)
        ckplm(k,685) = -1.8310546875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(-2556.D0+rk*(-247.D0+31.D0*rk)) * den(8)
        ckplm(k,686) = -41.50390625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(-93.D0+(rk-16.D0)*rk) * den(9)
        ckplm(k,687) = -15.771484375D0*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-175.D0+(rk-51.D0)*rk) * den(10)
        ckplm(k,688) = 31.54296875D0*(rk-10.D0)*(rk+24.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,689) = 181.3720703125D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,690) = -12.0914713541667D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,691) = 1.5771484375D0*(rk+10.D0)*(rk+11.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(69.D0+2.D0*rk) * den(1)
        ckplm(k,692) = 1.5771484375D0*(rk+10.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-327.D0+rk*(47.D0+4.D0*rk)) * den(2)
        ckplm(k,693) = 0.138346354166667D0*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(16530.D0+rk*(-5893.D0+rk*(33.D0+46.D0*rk))) * den(3)
        ckplm(k,694) = 0.0732421875D0*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(49470.D0+rk*(-10117.D0+rk*(-660.D0+67.D0*rk))) * den(4)
        ckplm(k,695) = 1.318359375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(3554.D0+rk*(-363.D0+rk*(-65.D0+2.D0*rk))) * den(5)
        ckplm(k,696) = -33.837890625D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)&
                *(-154.D0+3.D0*(rk+1.D0)*rk) * den(6)
        ckplm(k,697) = -1.318359375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(-3850.D0+rk*(-227.D0+rk*(71.D0+2.D0*rk))) * den(7)
        ckplm(k,698) = -0.0732421875D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(-58860.D0+rk*(-8596.D0+rk*(861.D0+67.D0*rk))) * den(8)
        ckplm(k,699) = -0.138346354166667D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(-22410.D0+rk*(-5821.D0+rk*(105.D0+46.D0*rk))) * den(9)
        ckplm(k,700) = -1.5771484375D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-370.D0+rk*(-39.D0+4.D0*rk)) * den(10)
        ckplm(k,701) = -1.5771484375D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-67.D0+2.D0*rk) * den(11)
        ckplm(k,702) = 12.0914713541667D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,703) = 0.755716959635417D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,704) = -0.131429036458333D0*(rk+10.D0)*(rk+11.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(92.D0+5.D0*rk) * den(1)
        ckplm(k,705) = -0.1971435546875D0*(rk+10.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-468.D0+rk*(5.D0+3.D0*rk)) * den(2)
        ckplm(k,706) = -0.0345865885416667D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(15960.D0+rk*(-2266.D0+rk*(-249.D0+7.D0*rk))) * den(3)
        ckplm(k,707) = 0.00050862630208333D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(6.17916D6+rk&
                *(-1.692866D6+rk*(-50819.D0+rk*(21686.D0+239.D0*rk)))) * den(4)
        ckplm(k,708) = 0.01220703125D0*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(368496.D0+rk&
                *(-50502.D0+rk*(-8155.D0+rk*(618.D0+31.D0*rk)))) * den(5)
        ckplm(k,709) = 0.469970703125D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(11088.D0+(rk&
                *rk+rk-290.D0)*(rk+1.D0)*rk) * den(6)
        ckplm(k,710) = 0.01220703125D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(410256.D0+rk&
                *(32462.D0+rk*(-9823.D0+rk*(-494.D0+31.D0*rk)))) * den(7)
        ckplm(k,711) = 0.00050862630208333D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(7.79976D6+rk&
                *(1.527126D6+rk*(-114443.D0+rk*(-20730.D0+239.D0*rk)))) * den(8)
        ckplm(k,712) = -0.0345865885416667D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(-17970.D0+rk*(-1747.D0+rk*(270.D0+7.D0*rk))) * den(9)
        ckplm(k,713) = -0.1971435546875D0*(3.D0*rk*rk+rk-470.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,714) = -0.131429036458333D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-87.D0+5.D0*rk) * den(11)
        ckplm(k,715) = 0.755716959635417D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,716) = -0.0444539388020833D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,717) = 0.00193277994791667D0*(rk+10.D0)*(rk+11.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(575.D0+38.D0*rk) * den(1)
        ckplm(k,718) = 0.00579833984375D0*(rk+10.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2175.D0+rk*(-129.D0+4.D0*rk)) * den(2)
        ckplm(k,719) = -0.00254313151041667D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-39330.D0+rk&
                *(305.D0+rk*(483.D0+10.D0*rk))) * den(3)
        ckplm(k,720) = -0.00254313151041667D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(267690.D0+rk&
                *(-25025.D0+rk*(-5723.D0+rk*(161.D0+17.D0*rk)))) * den(4)
        ckplm(k,721) = -0.0030517578125D0*(rk+6.D0)*(rk+7.D0)*(-1.414368D6+5.D0*rk*(48750.D0+rk&
                *(6995.D0+rk*(-908.D0+rk*(-47.D0+2.D0*rk))))) * den(5)
        ckplm(k,722) = 0.2349853515625D0*(rk-5.D0)*(rk+6.D0)*(22176.D0+5.D0*(rk*rk+rk-146.D0)&
                *(rk+1.D0)*rk) * den(6)
        ckplm(k,723) = 0.0030517578125D0*(rk-5.D0)*(rk-6.D0)*(1.618848D6+5.D0*rk*(32234.D0+rk&
                *(-9417.D0+rk*(-700.D0+rk*(57.D0+2.D0*rk))))) * den(7)
        ckplm(k,724) = 0.00254313151041667D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(286848.D0+rk*(13164.D0+rk&
                *(-6104.D0+rk*(-93.D0+17.D0*rk)))) * den(8)
        ckplm(k,725) = 0.00254313151041667D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(39162.D0+rk&
                *(-631.D0+rk*(-453.D0+10.D0*rk))) * den(9)
        ckplm(k,726) = -0.00579833984375D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-2042.D0+rk*(137.D0+4.D0*rk)) * den(10)
        ckplm(k,727) = -0.00193277994791667D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-537.D0+38.D0*rk) * den(11)
        ckplm(k,728) = 0.0444539388020833D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,729) = 0.00246966326678241D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,730) = -0.00128851996527778D0*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(69.D0+5.D0*rk) * den(1)
        ckplm(k,731) = 0.00009203714037698D0*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(15351.D0+rk&
                *(1535.D0+19.D0*rk)) * den(2)
        ckplm(k,732) = 0.0000161468667328D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-900600.D0+rk&
                *(-72992.D0+rk*(4152.D0+269.D0*rk))) * den(3)
        ckplm(k,733) = 0.00025431315104167D0*(rk+7.D0)*(rk+8.D0)*(463500.D0+rk*(17292.D0+rk&
                *(-7357.D0+rk*(-342.D0+7.D0*rk)))) * den(4)
        ckplm(k,734) = -0.00203450520833333D0*(rk+7.D0)*(406944.D0+rk*(-3246.D0+rk*(-11237.D0+rk&
                *(-157.D0+(rk+59.D0)*rk)))) * den(5)
        ckplm(k,735) = -0.00372992621527778D0*(-1.397088D6+(rk+1.D0)*rk*(55578.D0+(rk*rk+rk-575.D0)&
                *(rk+1.D0)*rk)) * den(6)
        ckplm(k,736) = -0.00203450520833333D0*(rk-6.D0)*(-399168.D0+rk*(18526.D0+rk*(10422.D0+rk&
                *(-383.D0+(rk-54.D0)*rk)))) * den(7)
        ckplm(k,737) = 0.00025431315104167D0*(rk-6.D0)*(rk-7.D0)*(439200.D0+rk*(-30952.D0+rk&
                *(-6289.D0+rk*(370.D0+7.D0*rk)))) * den(8)
        ckplm(k,738) = 0.0000161468667328D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(823725.D0+rk*(-80489.D0+rk&
                *(-3345.D0+269.D0*rk))) * den(9)
        ckplm(k,739) = 0.00009203714037698D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(13835.D0+rk&
                *(-1497.D0+19.D0*rk)) * den(10)
        ckplm(k,740) = -0.00128851996527778D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-64.D0+5.D0*rk) * den(11)
        ckplm(k,741) = 0.00246966326678241D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(12)
        ckplm(k,742) = -0.00012998227719907D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,743) = 5.65140335648148D-6*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)*(1127.D0+86.D0&
                *rk) * den(1)
        ckplm(k,744) = -8.07343336640212D-7*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(167727.D0+rk&
                *(21053.D0+556.D0*rk)) * den(2)
        ckplm(k,745) = -4.03671668320106D-6*(rk+8.D0)*(rk+9.D0)*(-438550.D0+rk*(-62359.D0+rk&
                *(-1301.D0+58.D0*rk))) * den(3)
        ckplm(k,746) = 0.00025431315104167D0*(rk+8.D0)*(-66850.D0+rk*(-9057.D0+rk*(249.D0+(rk+57.D0)&
                *rk))) * den(4)
        ckplm(k,747) = 0.00016954210069444D0*(780192.D0+rk*(95766.D0+rk*(-11125.D0+(2.D0*rk&
                *rk+rk-1484.D0)*rk))) * den(5)
        ckplm(k,748) = -0.0130547417534722D0*(9504.D0+(rk*rk+rk-218.D0)*(rk+1.D0)*rk) * den(6)
        ckplm(k,749) = -0.00016954210069444D0*(-674784.D0+rk*(113570.D0+rk*(6687.D0+rk*(-1468.D0+rk&
                *(9.D0+2.D0*rk))))) * den(7)
        ckplm(k,750) = -0.00025431315104167D0*(rk-7.D0)*(-57600.D0+rk*(9388.D0+rk*(84.D0+(rk-53.D0)&
                *rk))) * den(8)
        ckplm(k,751) = 4.03671668320106D-6*(rk-7.D0)*(rk-8.D0)*(377550.D0+rk*(-59583.D0+rk&
                *(1475.D0+58.D0*rk))) * den(9)
        ckplm(k,752) = 8.07343336640212D-7*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(147230.D0+rk&
                *(-19941.D0+556.D0*rk)) * den(10)
        ckplm(k,753) = -5.65140335648148D-6*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1041.D0+86.D0&
                *rk) * den(11)
        ckplm(k,754) = 0.00012998227719907D0*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(12)
        ckplm(k,755) = 6.4991138599537D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0) * den(0)
        ckplm(k,756) = -1.1302806712963D-6*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(368.D0+29.D0*rk) * den(1)
        ckplm(k,757) = 8.07343336640212D-8*(rk+10.D0)*(rk+9.D0)*(141456.D0+rk*(20159.D0+673.D0*rk)) &
                * den(2)
        ckplm(k,758) = -8.07343336640212D-7*(rk+9.D0)*(230240.D0+rk*(42758.D0+rk*(2167.D0+19.D0&
                *rk))) * den(3)
        ckplm(k,759) = 0.00001271565755208D0*(165280.D0+rk*(34866.D0+rk*(1663.D0-rk*(46.D0+3.D0&
                *rk)))) * den(4)
        ckplm(k,760) = 6.78168402777778D-6*(-333720.D0+rk*(-34842.D0+rk*(2099.D0+(rk+222.D0)*rk))) &
                * den(5)
        ckplm(k,761) = 0.00003729926215278D0*(59400.D0+(rk*rk+rk-722.D0)*(rk+1.D0)*rk) * den(6)
        ckplm(k,762) = 6.78168402777778D-6*(-297000.D0+rk*(38378.D0+rk*(1439.D0+(rk-218.D0)*rk))) &
                * den(7)
        ckplm(k,763) = 0.00001271565755208D0*(132120.D0+rk*(-31414.D0+rk*(1783.D0+(34.D0-3.D0*rk)&
                *rk))) * den(8)
        ckplm(k,764) = -8.07343336640212D-7*(rk-8.D0)*(-189630.D0+rk*(38481.D0+rk*(-2110.D0+19.D0&
                *rk))) * den(9)
        ckplm(k,765) = 8.07343336640212D-8*(rk-8.D0)*(rk-9.D0)*(121970.D0+rk*(-18813.D0+673.D0*rk)) &
                * den(10)
        ckplm(k,766) = 1.1302806712963D-6*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(339.D0-29.D0*rk) * den(11)
        ckplm(k,767) = 6.4991138599537D-6*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,768) = -3.09481612378748D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0) * den(0)
        ckplm(k,769) = 4.03671668320106D-8*(rk+10.D0)*(rk+11.D0)*(621.D0+50.D0*rk) * den(1)
        ckplm(k,770) = -4.03671668320106D-8*(rk+10.D0)*(21411.D0+rk*(3305.D0+124.D0*rk)) * den(2)
        ckplm(k,771) = 6.72786113866843D-8*(255690.D0+rk*(55633.D0+rk*(3747.D0+74.D0*rk))) * den(3)
        ckplm(k,772) = 6.05507502480159D-7*(-41730.D0+rk*(-6631.D0+(rk-240.D0)*rk)) * den(4)
        ckplm(k,773) = 2.17982700892857D-6*(13810.D0+rk*(1131.D0-rk*(31.D0+2.D0*rk))) * den(5)
        ckplm(k,774) = 0.00005594889322917D0*(-550.D0+3.D0*(rk+1.D0)*rk) * den(6)
        ckplm(k,775) = 2.17982700892857D-6*(12650.D0+rk*(-1187.D0+rk*(-25.D0+2.D0*rk))) * den(7)
        ckplm(k,776) = -6.05507502480159D-7*(35340.D0+rk*(-6148.D0+(rk+243.D0)*rk)) * den(8)
        ckplm(k,777) = -6.72786113866843D-8*(-203730.D0+rk*(48361.D0+rk*(-3525.D0+74.D0*rk))) * den(9)
        ckplm(k,778) = 4.03671668320106D-8*(rk-9.D0)*(18230.D0+rk*(-3057.D0+124.D0*rk)) * den(10)
        ckplm(k,779) = 4.03671668320106D-8*(rk-10.D0)*(rk-9.D0)*(571.D0-50.D0*rk) * den(11)
        ckplm(k,780) = 3.09481612378748D-7*(rk-10.D0)*(rk-11.D0)*(rk-9.D0) * den(12)
        ckplm(k,781) = 1.40673460172158D-8*(rk+11.D0)*(rk+12.D0) * den(0)
        ckplm(k,782) = -2.44649495951579D-9*(rk+11.D0)*(575.D0+47.D0*rk) * den(1)
        ckplm(k,783) = 1.34557222773369D-8*(4425.D0+rk*(721.D0+29.D0*rk)) * den(2)
        ckplm(k,784) = -6.72786113866843D-7*(212.D0+(rk+30.D0)*rk) * den(3)
        ckplm(k,785) = 5.04589585400132D-7*(482.D0+(rk+51.D0)*rk) * den(4)
        ckplm(k,786) = 1.61468667328042D-7*(-1983.D0+(rk-112.D0)*rk) * den(5)
        ckplm(k,787) = 5.65140335648148D-7*(-rk*rk-rk+605.D0) * den(6)
        ckplm(k,788) = 1.61468667328042D-7*(-1870.D0+(rk+114.D0)*rk) * den(7)
        ckplm(k,789) = 5.04589585400132D-7*(432.D0+(rk-49.D0)*rk) * den(8)
        ckplm(k,790) = -6.72786113866843D-7*(183.D0+(rk-28.D0)*rk) * den(9)
        ckplm(k,791) = 1.34557222773369D-8*(3733.D0+rk*(-663.D0+29.D0*rk)) * den(10)
        ckplm(k,792) = 2.44649495951579D-9*(rk-10.D0)*(528.D0-47.D0*rk) * den(11)
        ckplm(k,793) = 1.40673460172158D-8*(rk-10.D0)*(rk-11.D0) * den(12)
        ckplm(k,794) = 6.11623739878948D-10*(-12.D0-rk) * den(0)
        ckplm(k,795) = 6.11623739878948D-10*(121.D0+10.D0*rk) * den(1)
        ckplm(k,796) = 6.72786113866843D-9*(-51.D0-4.D0*rk) * den(2)
        ckplm(k,797) = 3.36393056933422D-8*(29.D0+2.D0*rk) * den(3)
        ckplm(k,798) = -1.00917917080027D-7*(rk+19.D0) * den(4)
        ckplm(k,799) = 4.03671668320106D-8*(69.D0+2.D0*rk) * den(5)
        ckplm(k,800) = -3.10827184606481D-6 * den(6)
        ckplm(k,801) = 4.03671668320106D-8*(67.D0-2.D0*rk) * den(7)
        ckplm(k,802) = 1.00917917080027D-7*(rk-18.D0) * den(8)
        ckplm(k,803) = -3.36393056933422D-8*(-27.D0+2.D0*rk) * den(9)
        ckplm(k,804) = 6.72786113866843D-9*(-47.D0+4.D0*rk) * den(10)
        ckplm(k,805) = 6.11623739878948D-10*(111.D0-10.D0*rk) * den(11)
        ckplm(k,806) = 6.11623739878948D-10*(rk-11.D0) * den(12)
        ckplm(k,807) = 2.54843224949562D-11 * den(0)
        ckplm(k,808) = -3.05811869939474D-10 * den(1)
        ckplm(k,809) = 1.68196528466711D-9 * den(2)
        ckplm(k,810) = -5.60655094889036D-9 * den(3)
        ckplm(k,811) = 1.26147396350033D-8 * den(4)
        ckplm(k,812) = -2.01835834160053D-8 * den(5)
        ckplm(k,813) = 2.35475139853395D-8 * den(6)
        ckplm(k,814) = -2.01835834160053D-8 * den(7)
        ckplm(k,815) = 1.26147396350033D-8 * den(8)
        ckplm(k,816) = -5.60655094889036D-9 * den(9)
        ckplm(k,817) = 1.68196528466711D-9 * den(10)
        ckplm(k,818) = -3.05811869939474D-10 * den(11)
        ckplm(k,819) = 2.54843224949562D-11 * den(12)
!    ckplm para l = 13
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-1.D0)*(r2k+1.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,820) = 34279.3212890625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+1.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,821) = 17825.2470703125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+1.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,822) = 13950.193359375D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,823) = 12178.740234375D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,824) = 11217.2607421875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,825) = 10689.3896484375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*rk * den(5)
        ckplm(k,826) = 10451.84765625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*rk * den(6)
        ckplm(k,827) = 10451.84765625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*rk * den(7)
        ckplm(k,828) = 10689.3896484375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*rk * den(8)
        ckplm(k,829) = 11217.2607421875D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*rk * den(9)
        ckplm(k,830) = 12178.740234375D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,831) = 13950.193359375D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,832) = 17825.2470703125D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,833) = 34279.3212890625D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,834) = -4897.0458984375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,835) = -195.8818359375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-25.D0+11.D0*rk) * den(1)
        ckplm(k,836) = -153.298828125D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-46.D0+9.D0*rk) * den(2)
        ckplm(k,837) = -936.826171875D0*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(3)
        ckplm(k,838) = -123.2666015625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-76.D0+5.D0*rk) * den(4)
        ckplm(k,839) = -117.4658203125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-85.D0+3.D0*rk) * den(5)
        ckplm(k,840) = -114.85546875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk-90.D0) * den(6)
        ckplm(k,841) = 114.85546875D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+91.D0) * den(7)
        ckplm(k,842) = 117.4658203125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(88.D0+3.D0*rk) * den(8)
        ckplm(k,843) = 123.2666015625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(81.D0+5.D0*rk) * den(9)
        ckplm(k,844) = 936.826171875D0*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,845) = 153.298828125D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(55.D0+9.D0*rk) * den(11)
        ckplm(k,846) = 195.8818359375D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(36.D0+11.D0*rk) * den(12)
        ckplm(k,847) = 4897.0458984375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,848) = 326.4697265625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,849) = 65.2939453125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-20.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,850) = -17.033203125D0*(rk+10.D0)*(rk+11.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-207.D0+(rk+91.D0)*rk) * den(2)
        ckplm(k,851) = -62.455078125D0*(rk+10.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-93.D0+(rk+23.D0)*rk) * den(3)
        ckplm(k,852) = -8.2177734375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-950.D0+rk*(141.D0+11.D0*rk)) * den(4)
        ckplm(k,853) = -2.6103515625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-3570.D0+rk*(299.D0+41.D0*rk)) * den(5)
        ckplm(k,854) = -114.85546875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(-89.D0+(rk+3.D0)*rk) * den(6)
        ckplm(k,855) = -114.85546875D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(-91.D0+(rk-1.D0)*rk) * den(7)
        ckplm(k,856) = -2.6103515625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(-3828.D0+rk*(-217.D0+41.D0*rk)) * den(8)
        ckplm(k,857) = -8.2177734375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-1080.D0+rk*(-119.D0+11.D0*rk)) * den(9)
        ckplm(k,858) = -62.455078125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-115.D0+(rk-21.D0)*rk) * den(10)
        ckplm(k,859) = -17.033203125D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-297.D0+(rk-89.D0)*rk) * den(11)
        ckplm(k,860) = 65.2939453125D0*(rk-10.D0)*(rk-11.D0)*(rk+21.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,861) = 326.4697265625D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,862) = -20.4043579101563D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,863) = 4.08087158203125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+45.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,864) = 3.1937255859375D0*(rk+10.D0)*(rk+11.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-276.D0+rk*(43.D0+3.D0*rk)) * den(2)
        ckplm(k,865) = 0.3548583984375D0*(rk+10.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(11244.D0+rk*(-4181.D0+rk*(84.D0+29.D0*rk))) * den(3)
        ckplm(k,866) = 0.04669189453125D0*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(139080.D0+rk*(-31094.D0+rk*(-1239.D0+185.D0*rk))) * den(4)
        ckplm(k,867) = 0.04449462890625D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(195330.D0+rk*(-24651.D0+rk*(-2846.D0+127.D0*rk))) * den(5)
        ckplm(k,868) = 0.652587890625D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(15488.D0+3.D0*rk*(-262.D0+(rk-85.D0)*rk)) * den(6)
        ckplm(k,869) = -0.652587890625D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(-16016.D0+3.D0*(rk-1.D0)*(rk+89.D0)*rk) * den(7)
        ckplm(k,870) = -0.04449462890625D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(-217008.D0+rk*(-18578.D0+rk*(3227.D0+127.D0*rk))) * den(8)
        ckplm(k,871) = -0.04669189453125D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(-168750.D0+rk*(-28061.D0+rk*(1794.D0+185.D0*rk))) * den(9)
        ckplm(k,872) = -0.3548583984375D0*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-15480.D0+rk*(-4262.D0+rk*(3.D0+29.D0*rk))) * den(10)
        ckplm(k,873) = -3.1937255859375D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-316.D0+rk*(-37.D0+3.D0*rk)) * den(11)
        ckplm(k,874) = -4.08087158203125D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-44.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,875) = 20.4043579101563D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,876) = 1.20025634765625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,877) = -0.04801025390625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(400.D0+19.D0*rk) * den(1)
        ckplm(k,878) = -0.0125244140625D0*(rk+10.D0)*(rk+11.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-11868.D0+rk*(295.D0+73.D0*rk)) * den(2)
        ckplm(k,879) = -0.0208740234375D0*(rk+10.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(43752.D0+rk*(-7070.D0+rk*(-549.D0+23.D0*rk))) * den(3)
        ckplm(k,880) = 0.00054931640625D0*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(9.81996D6+rk*(-2.936746D6+rk*(-4543.D0+rk*(31966.D0+43.D0*rk)))) * den(4)
        ckplm(k,881) = 0.01483154296875D0*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(546120.D0+rk*(-92358.D0+rk*(-8665.D0+rk*(994.D0+29.D0*rk)))) * den(5)
        ckplm(k,882) = 0.652587890625D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(15312.D0+rk*(-1042.D0+rk*(-329.D0+(rk+10.D0)*rk))) * den(6)
        ckplm(k,883) = 0.652587890625D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(16016.D0+(rk-1.D0)*rk*(-358.D0+(rk-5.D0)*rk)) * den(7)
        ckplm(k,884) = 0.01483154296875D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(628848.D0+rk*(72162.D0+rk*(-11473.D0+rk*(-878.D0+29.D0*rk)))) * den(8)
        ckplm(k,885) = 0.00054931640625D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(1.272024D7+rk*(2.831934D6+rk*(-100183.D0+rk*(-31794.D0+43.D0*rk)))) * den(9)
        ckplm(k,886) = -0.0208740234375D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-50250.D0+rk*(-5903.D0+rk*(618.D0+23.D0*rk))) * den(10)
        ckplm(k,887) = -0.0125244140625D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-12090.D0+rk*(-149.D0+73.D0*rk)) * den(11)
        ckplm(k,888) = -0.04801025390625D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-381.D0+19.D0*rk) * den(12)
        ckplm(k,889) = 1.20025634765625D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,890) = -0.066680908203125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,891) = 0.002667236328125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(625.D0+37.D0*rk) * den(1)
        ckplm(k,892) = 0.00069580078125D0*(rk+10.D0)*(rk+11.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-27600.D0+rk*(-1211.D0+61.D0*rk)) * den(2)
        ckplm(k,893) = -0.00115966796875D0*(rk+10.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-135480.D0+rk*(4193.D0+rk*(1548.D0+19.D0*rk))) * den(3)
        ckplm(k,894) = -0.000152587890625D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(7.33476D6+rk&
                *(-911570.D0+rk*(-125783.D0+rk*(6434.D0+359.D0*rk)))) * den(4)
        ckplm(k,895) = -0.001373291015625D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-5.493096D6+rk&
                *(1.16637D6+rk*(83741.D0+rk*(-18415.D0+rk*(-413.D0+37.D0*rk))))) * den(5)
        ckplm(k,896) = -0.0040283203125D0*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(-2.452032D6+5.D0*rk&
                *(41952.D0+rk*(12874.D0+rk*(-643.D0+(rk-76.D0)*rk)))) * den(6)
        ckplm(k,897) = 0.0040283203125D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(2.594592D6+5.D0*(rk-1.D0)*rk&
                *(-14584.D0+rk*(-247.D0+(rk+82.D0)*rk))) * den(7)
        ckplm(k,898) = 0.001373291015625D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(6.55776D6+rk*(945480.D0+rk&
                *(-136138.D0+rk*(-16393.D0+rk*(598.D0+37.D0*rk))))) * den(8)
        ckplm(k,899) = 0.000152587890625D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(8.114472D6+rk&
                *(642138.D0+rk*(-142931.D0+rk*(-4998.D0+359.D0*rk)))) * den(9)
        ckplm(k,900) = 0.00115966796875D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(138144.D0+rk*(1154.D0+rk*(-1491.D0+19.D0*rk))) * den(10)
        ckplm(k,901) = -0.00069580078125D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-26328.D0+rk*(1333.D0+61.D0*rk)) * den(11)
        ckplm(k,902) = -0.002667236328125D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-588.D0+37.D0*rk) * den(12)
        ckplm(k,903) = 0.066680908203125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,904) = 0.003509521484375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,905) = -0.000140380859375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(900.D0+59.D0*rk) * den(1)
        ckplm(k,906) = 0.00010986328125D0*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(18561.D0+rk*(1555.D0+9.D0*rk)) * den(2)
        ckplm(k,907) = 0.00006103515625D0*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-355854.D0+rk&
                *(-19567.D0+rk*(1977.D0+85.D0*rk))) * den(3)
        ckplm(k,908) = 0.000030517578125D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(6.0849D6+rk*(7524.D0+rk&
                *(-93023.D0+rk*(-2334.D0+113.D0*rk)))) * den(4)
        ckplm(k,909) = -0.000823974609375D0*(rk+7.D0)*(rk+8.D0)*(1.70424D6+rk*(-95052.D0+rk&
                *(-40552.D0+rk*(701.D0+(rk+202.D0)*rk)))) * den(5)
        ckplm(k,910) = -0.0040283203125D0*(rk+7.D0)*(-2.42352D6+rk*(250194.D0+rk*(74515.D0+rk&
                *(-5271.D0+rk*(-644.D0+(rk+21.D0)*rk))))) * den(6)
        ckplm(k,911) = -0.0040283203125D0*(rk-6.D0)*(-2.594592D6+(rk-1.D0)*rk*(88026.D0+rk&
                *(1757.D0+rk*(-748.D0+(rk-14.D0)*rk)))) * den(7)
        ckplm(k,912) = -0.000823974609375D0*(rk-6.D0)*(rk-7.D0)*(-1.75824D6+rk*(-12648.D0+rk&
                *(41453.D0+rk*(-97.D0+(rk-197.D0)*rk)))) * den(8)
        ckplm(k,913) = 0.000030517578125D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(5.9868D6+rk*(-186116.D0+rk&
                *(-85343.D0+rk*(2786.D0+113.D0*rk)))) * den(9)
        ckplm(k,914) = 0.00006103515625D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(334395.D0+rk&
                *(-23266.D0+rk*(-1722.D0+85.D0*rk))) * den(10)
        ckplm(k,915) = 0.00010986328125D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(17015.D0+rk*(-1537.D0+9.D0*rk)) * den(11)
        ckplm(k,916) = -0.000140380859375D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-841.D0+59.D0*rk) * den(12)
        ckplm(k,917) = 0.003509521484375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,918) = -0.00017547607421875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,919) = 0.00003509521484375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(245.D0+17.D0*rk) * den(1)
        ckplm(k,920) = -9.1552734375D-6*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)*(20286.D0+rk&
                *(2207.D0+47.D0*rk)) * den(2)
        ckplm(k,921) = -4.35965401785714D-7*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(-5.776806D6+rk&
                *(-657931.D0+rk*(-3516.D0+829.D0*rk))) * den(3)
        ckplm(k,922) = 7.62939453125D-6*(rk+8.D0)*(rk+9.D0)*(-3.36588D6+rk*(-318502.D0+rk&
                *(22385.D0+rk*(2266.D0+19.D0*rk)))) * den(4)
        ckplm(k,923) = 0.00001373291015625D0*(rk+8.D0)*(1.585332D7+rk*(1.067286D6+rk*(-263915.D0+rk&
                *(-18515.D0+rk*(515.D0+29.D0*rk))))) * den(5)
        ckplm(k,924) = 0.000201416015625D0*(-7.98336D6+rk*(-363588.D0+rk*(218464.D0+rk*(10485.D0+rk&
                *(-1385.D0+(rk-57.D0)*rk))))) * den(6)
        ckplm(k,925) = -0.000201416015625D0*(-7.41312D6+rk*(763812.D0+rk*(179284.D0+rk*(-15435.D0+rk&
                *(-1085.D0+(rk+63.D0)*rk))))) * den(7)
        ckplm(k,926) = -0.00001373291015625D0*(rk-7.D0)*(-1.454112D7+rk*(1.537656D6+rk*(205570.D0+rk&
                *(-20285.D0+rk*(-370.D0+29.D0*rk))))) * den(8)
        ckplm(k,927) = 7.62939453125D-6*(rk-7.D0)*(rk-8.D0)*(3.02724D6-rk*(356550.D0+rk*(15701.D0+rk&
                *(-2190.D0+19.D0*rk)))) * den(9)
        ckplm(k,928) = 4.35965401785714D-7*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(5.12322D6+rk*(-648412.D0+rk&
                *(6003.D0+829.D0*rk))) * den(10)
        ckplm(k,929) = 9.1552734375D-6*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(18126.D0+rk&
                *(-2113.D0+47.D0*rk)) * den(11)
        ckplm(k,930) = -0.00003509521484375D0*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-228.D0+17.D0*rk) * den(12)
        ckplm(k,931) = 0.00017547607421875D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(13)
        ckplm(k,932) = 8.35600353422619D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,933) = -1.67120070684524D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)*(320.D0+23.D0&
                *rk) * den(1)
        ckplm(k,934) = 1.30789620535714D-6*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(11408.D0+rk&
                *(1431.D0+41.D0*rk)) * den(2)
        ckplm(k,935) = -1.45321800595238D-7*(rk+10.D0)*(rk+9.D0)*(1.736352D6+rk*(270302.D0+rk&
                *(10347.D0+7.D0*rk))) * den(3)
        ckplm(k,936) = -3.63304501488095D-7*(rk+9.D0)*(-8.38416D6+rk*(-1.39037D6+rk*(-32711.D0+rk&
                *(3326.D0+107.D0*rk)))) * den(4)
        ckplm(k,937) = -6.53948102678571D-7*(4.421952D7+rk*(7.366056D6+rk*(-93910.D0+rk*(-58375.D0+rk&
                *(-1790.D0+19.D0*rk))))) * den(5)
        ckplm(k,938) = 0.00002877371651786D0*(rk-9.D0)*(-109560.D0+rk*(-17138.D0+rk&
                *(-193.D0+(rk+50.D0)*rk))) * den(6)
        ckplm(k,939) = 0.00002877371651786D0*(rk+10.D0)*(-92664.D0+rk*(16606.D0+rk&
                *(-337.D0+(rk-46.D0)*rk))) * den(7)
        ckplm(k,940) = -6.53948102678571D-7*(-3.681612D7+rk*(7.386006D6+rk*(-70285.D0+rk&
                *(-51025.D0+rk*(1885.D0+19.D0*rk))))) * den(8)
        ckplm(k,941) = 3.63304501488095D-7*(rk-8.D0)*(7.02972D6+rk*(-1.315398D6+rk&
                *(42047.D0+(2898.D0-107.D0*rk)*rk))) * den(9)
        ckplm(k,942) = 1.45321800595238D-7*(rk-8.D0)*(rk-9.D0)*(1.47639D6+rk&
                *(-249629.D0+(10326.D0-7.D0*rk)*rk)) * den(10)
        ckplm(k,943) = 1.30789620535714D-6*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(10018.D0+rk&
                *(-1349.D0+41.D0*rk)) * den(11)
        ckplm(k,944) = -1.67120070684524D-6*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)*(-297.D0+23.D0&
                *rk) * den(12)
        ckplm(k,945) = 8.35600353422619D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(13)
        ckplm(k,946) = -3.79818342464827D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0) * den(0)
        ckplm(k,947) = 1.51927336985931D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(2025.D0+149.D0*rk) &
                * den(1)
        ckplm(k,948) = -1.18899655032468D-8*(rk+10.D0)*(rk+11.D0)*(90666.D0+rk*(12425.D0+409.D0*rk)) &
                * den(2)
        ckplm(k,949) = 7.26609002976191D-8*(rk+10.D0)*(307206.D0+rk*(57199.D0+rk*(3180.D0+47.D0&
                *rk))) * den(3)
        ckplm(k,950) = 3.63304501488095D-8*(-8.71236D6+rk*(-1.926762D6+rk*(-125041.D0+rk&
                *(-1338.D0+61.D0*rk)))) * den(4)
        ckplm(k,951) = -2.94276646205357D-6*(-126840.D0+rk*(-17370.D0+rk*(-37.D0+(rk+54.D0)*rk))) &
                * den(5)
        ckplm(k,952) = -7.84737723214286D-7*(496100.D0+rk*(20300.D0+3.D0*rk*(-1343.D0+(rk-42.D0)&
                *rk))) * den(6)
        ckplm(k,953) = 7.84737723214286D-7*(471900.D0+rk*(-27968.D0+3.D0*rk*(-1211.D0+(rk+46.D0)&
                *rk))) * den(7)
        ckplm(k,954) = 2.94276646205357D-6*(-109560.D0+rk*(17138.D0+rk*(-193.D0+(rk-50.D0)*rk))) &
                * den(8)
        ckplm(k,955) = 3.63304501488095D-8*(6.90924D6-rk*(1.680938D6+rk*(-120661.D0+rk*(1582.D0+61.D0&
                *rk)))) * den(9)
        ckplm(k,956) = 7.26609002976191D-8*(rk-9.D0)*(253140.D0+rk*(-50980.D0+(3039.D0-47.D0*rk)&
                *rk)) * den(10)
        ckplm(k,957) = 1.18899655032468D-8*(rk-10.D0)*(rk-9.D0)*(78650.D0+rk*(-11607.D0+409.D0*rk)) &
                * den(11)
        ckplm(k,958) = -1.51927336985931D-8*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(-1876.D0+149.D0*rk) &
                * den(12)
        ckplm(k,959) = 3.79818342464827D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0) * den(13)
        ckplm(k,960) = 1.65138409767316D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0) * den(0)
        ckplm(k,961) = -6.60553639069264D-10*(rk+11.D0)*(rk+12.D0)*(2500.D0+187.D0*rk) * den(1)
        ckplm(k,962) = 1.18899655032468D-8*(rk+11.D0)*(5975.D0+rk*(869.D0+31.D0*rk)) * den(2)
        ckplm(k,963) = 7.26609002976191D-8*(-24390.D0-rk*(5099.D0+rk*(339.D0+7.D0*rk))) * den(3)
        ckplm(k,964) = 1.81652250744048D-7*(16122.D0+rk*(2657.D0+(rk+120.D0)*rk)) * den(4)
        ckplm(k,965) = 3.26974051339286D-7*(rk+12.D0)*(-982.D0+(rk-21.D0)*rk) * den(5)
        ckplm(k,966) = -2.61579241071429D-7*(-16335.D0+rk*(-548.D0+(rk+58.D0)*rk)) * den(6)
        ckplm(k,967) = -2.61579241071429D-7*(15730.D0+rk*(-661.D0+(rk-55.D0)*rk)) * den(7)
        ckplm(k,968) = 3.26974051339286D-7*(rk-11.D0)*(-960.D0+(rk+23.D0)*rk) * den(8)
        ckplm(k,969) = 1.81652250744048D-7*(-13584.D0+rk*(2420.D0+(rk-117.D0)*rk)) * den(9)
        ckplm(k,970) = 7.26609002976191D-8*(19623.D0+rk*(-4442.D0+(318.D0-7.D0*rk)*rk)) * den(10)
        ckplm(k,971) = 1.18899655032468D-8*(rk-10.D0)*(5137.D0+rk*(-807.D0+31.D0*rk)) * den(11)
        ckplm(k,972) = 6.60553639069264D-10*(rk-10.D0)*(rk-11.D0)*(2313.D0-187.D0*rk) * den(12)
        ckplm(k,973) = 1.65138409767316D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0) * den(13)
        ckplm(k,974) = -6.88076707363817D-10*(rk+12.D0)*(rk+13.D0) * den(0)
        ckplm(k,975) = 2.75230682945527D-11*(rk+12.D0)*(3025.D0+229.D0*rk) * den(1)
        ckplm(k,976) = 1.65138409767316D-10*(-26136.D0-rk*(3965.D0+149.D0*rk)) * den(2)
        ckplm(k,977) = 3.02753751240079D-9*(3840.D0+rk*(521.D0+17.D0*rk)) * den(3)
        ckplm(k,978) = 1.5137687562004D-9*(-14676.D0-rk*(1597.D0+37.D0*rk)) * den(4)
        ckplm(k,979) = 1.36239188058036D-8*(2388.D0+(rk+169.D0)*rk) * den(5)
        ckplm(k,980) = 3.99634951636905D-8*(-960.D0+(rk-23.D0)*rk) * den(6)
        ckplm(k,981) = -3.99634951636905D-8*(-936.D0+(rk+25.D0)*rk) * den(7)
        ckplm(k,982) = 1.36239188058036D-8*(-2220.D0-(rk-167.D0)*rk) * den(8)
        ckplm(k,983) = 1.5137687562004D-9*(13116.D0+rk*(-1523.D0+37.D0*rk)) * den(9)
        ckplm(k,984) = 3.02753751240079D-9*(-3336.D0+(487.D0-17.D0*rk)*rk) * den(10)
        ckplm(k,985) = 1.65138409767316D-10*(22320.D0+rk*(-3667.D0+149.D0*rk)) * den(11)
        ckplm(k,986) = 2.75230682945527D-11*(rk-11.D0)*(2796.D0-229.D0*rk) * den(12)
        ckplm(k,987) = 6.88076707363817D-10*(rk-11.D0)*(rk-12.D0) * den(13)
        ckplm(k,988) = 2.75230682945527D-11*(rk+13.D0) * den(0)
        ckplm(k,989) = 2.75230682945527D-11*(-144.D0-11.D0*rk) * den(1)
        ckplm(k,990) = 4.95415229301948D-10*(41.D0+3.D0*rk) * den(2)
        ckplm(k,991) = 6.05507502480159D-10*(-106.D0-7.D0*rk) * den(3)
        ckplm(k,992) = 1.5137687562004D-9*(93.D0+5.D0*rk) * den(4)
        ckplm(k,993) = -8.17435128348214D-9*(rk+28.D0) * den(5)
        ckplm(k,994) = 3.63304501488095D-9*(rk+79.D0) * den(6)
        ckplm(k,995) = 3.63304501488095D-9*(rk-78.D0) * den(7)
        ckplm(k,996) = -8.17435128348214D-9*(rk-27.D0) * den(8)
        ckplm(k,997) = 1.5137687562004D-9*(-88.D0+5.D0*rk) * den(9)
        ckplm(k,998) = 6.05507502480159D-10*(99.D0-7.D0*rk) * den(10)
        ckplm(k,999) = 4.95415229301948D-10*(-38.D0+3.D0*rk) * den(11)
        ckplm(k,1000) = 2.75230682945527D-11*(133.D0-11.D0*rk) * den(12)
        ckplm(k,1001) = 2.75230682945527D-11*(rk-12.D0) * den(13)
        ckplm(k,1002) = -1.05857954979049D-12 * den(0)
        ckplm(k,1003) = 1.37615341472763D-11 * den(1)
        ckplm(k,1004) = -8.2569204883658D-11 * den(2)
        ckplm(k,1005) = 3.02753751240079D-10 * den(3)
        ckplm(k,1006) = -7.56884378100199D-10 * den(4)
        ckplm(k,1007) = 1.36239188058036D-9 * den(5)
        ckplm(k,1008) = -1.81652250744048D-9 * den(6)
        ckplm(k,1009) = 1.81652250744048D-9 * den(7)
        ckplm(k,1010) = -1.36239188058036D-9 * den(8)
        ckplm(k,1011) = 7.56884378100199D-10 * den(9)
        ckplm(k,1012) = -3.02753751240079D-10 * den(10)
        ckplm(k,1013) = 8.2569204883658D-11 * den(11)
        ckplm(k,1014) = -1.37615341472763D-11 * den(12)
        ckplm(k,1015) = 1.05857954979049D-12 * den(13)
!    ckplm para l = 14
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        ckplm(k,1016) = 71007.1655273437D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,1017) = 36818.5302734375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+1.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,1018) = 28718.4536132813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,1019) = 24972.568359375D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,1020) = 22891.5209960938D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,1021) = 21686.7041015625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,1022) = 21048.8598632813D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*rk * den(6)
        ckplm(k,1023) = 20848.39453125D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*rk * den(7)
        ckplm(k,1024) = 21048.8598632813D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*rk * den(8)
        ckplm(k,1025) = 21686.7041015625D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*rk * den(9)
        ckplm(k,1026) = 22891.5209960938D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,1027) = 24972.568359375D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,1028) = 28718.4536132813D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,1029) = 36818.5302734375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,1030) = 71007.1655273437D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,1031) = -9467.6220703125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1032) = -1051.9580078125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-9.D0+4.D0*rk) * den(1)
        ckplm(k,1033) = -2735.0908203125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(2)
        ckplm(k,1034) = -237.833984375D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-69.D0+8.D0*rk) * den(3)
        ckplm(k,1035) = -1308.0869140625D0*(rk+10.D0)*(rk-14.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(4)
        ckplm(k,1036) = -206.5400390625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-95.D0+4.D0*rk) * den(5)
        ckplm(k,1037) = -400.9306640625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-51.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) * den(6)
        ckplm(k,1038) = 20848.39453125D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0) * den(7)
        ckplm(k,1039) = 400.9306640625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+52.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0) * den(8)
        ckplm(k,1040) = 206.5400390625D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(99.D0+4.D0*rk) * den(9)
        ckplm(k,1041) = 1308.0869140625D0*(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,1042) = 237.833984375D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(77.D0+8.D0*rk) * den(11)
        ckplm(k,1043) = 2735.0908203125D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,1044) = 1051.9580078125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(13.D0+4.D0*rk) * den(13)
        ckplm(k,1045) = 9467.6220703125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1046) = 591.726379394531D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1047) = 131.494750976563D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-18.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,1048) = -13.1494750976563D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-490.D0+(rk+219.D0)*rk) * den(2)
        ckplm(k,1049) = -2.286865234375D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-4692.D0+rk*(1201.D0+41.D0*rk)) * den(3)
        ckplm(k,1050) = -6.28887939453125D0*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2324.D0+rk*(369.D0+23.D0*rk)) * den(4)
        ckplm(k,1051) = -1.9859619140625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-8930.D0+rk*(861.D0+89.D0*rk)) * den(5)
        ckplm(k,1052) = -1.92755126953125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-10302.D0+rk*(511.D0+101.D0*rk)) * den(6)
        ckplm(k,1053) = -200.46533203125D0*(rk*rk+rk-104.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0) * den(7)
        ckplm(k,1054) = -1.92755126953125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(-10712.D0+rk*(-309.D0+101.D0*rk)) * den(8)
        ckplm(k,1055) = -1.9859619140625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-9702.D0+rk*(-683.D0+89.D0*rk)) * den(9)
        ckplm(k,1056) = -6.28887939453125D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2670.D0+rk*(-323.D0+23.D0*rk)) * den(10)
        ckplm(k,1057) = -2.286865234375D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-5852.D0+rk*(-1119.D0+41.D0*rk)) * den(11)
        ckplm(k,1058) = -13.1494750976563D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-708.D0+(rk-217.D0)*rk) * den(12)
        ckplm(k,1059) = 131.494750976563D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk+19.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,1060) = 591.726379394531D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1061) = -34.8074340820313D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1062) = 1.28916422526042D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(243.D0+4.D0*rk) * den(1)
        ckplm(k,1063) = 0.77349853515625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1965.D0+rk*(326.D0+19.D0*rk)) * den(2)
        ckplm(k,1064) = 0.0672607421875D0*(rk+10.D0)*(rk+11.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(104052.D0+rk*(-40007.D0+rk*(1257.D0+248.D0*rk))) * den(3)
        ckplm(k,1065) = 0.123311360677083D0*(rk+10.D0)*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(94416.D0+rk*(-22564.D0+rk*(-489.D0+121.D0*rk))) * den(4)
        ckplm(k,1066) = 0.05841064453125D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(274170.D0+rk*(-39809.D0+rk*(-3129.D0+188.D0*rk))) * den(5)
        ckplm(k,1067) = 5.78265380859375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(3334.D0+rk*(-249.D0+(rk-46.D0)*rk)) * den(6)
        ckplm(k,1068) = -100.232666015625D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(-208.D0+3.D0*(rk+1.D0)*rk) * den(7)
        ckplm(k,1069) = -5.78265380859375D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(-3536.D0+rk*(-154.D0+(rk+49.D0)*rk)) * den(8)
        ckplm(k,1070) = -0.05841064453125D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-310662.D0+rk*(-32987.D0+rk*(3693.D0+188.D0*rk))) * den(9)
        ckplm(k,1071) = -0.123311360677083D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-116370.D0+rk*(-21223.D0+rk*(852.D0+121.D0*rk))) * den(10)
        ckplm(k,1072) = -0.0672607421875D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-145068.D0+rk*(-41777.D0+rk*(-513.D0+248.D0*rk))) * den(11)
        ckplm(k,1073) = -0.77349853515625D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2272.D0+rk*(-288.D0+19.D0*rk)) * den(12)
        ckplm(k,1074) = -1.28916422526042D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-239.D0+4.D0*rk) * den(13)
        ckplm(k,1075) = 34.8074340820313D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1076) = 1.93374633789063D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1077) = -1.28916422526042D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+24.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,1078) = -0.128916422526042D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1880.D0+rk*(69.D0+11.D0*rk)) * den(2)
        ckplm(k,1079) = -0.0672607421875D0*(rk+10.D0)*(rk+11.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(22632.D0+rk*(-4022.D0+rk*(-223.D0+13.D0*rk))) * den(3)
        ckplm(k,1080) = -0.00560506184895833D0*(rk+10.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-1.652616D6+rk*(527670.D0+rk*(-10075.D0+rk*(-4890.D0+31.D0*rk)))) * den(4)
        ckplm(k,1081) = 0.00059000651041667D0*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(2.448948D7+rk*(-4.759886D6+rk*(-256109.D0+rk*(44834.D0+761.D0*rk)))) * den(5)
        ckplm(k,1082) = 0.009735107421875D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(1.922184D6+rk*(-192318.D0+rk*(-33077.D0+rk*(1722.D0+89.D0*rk)))) * den(6)
        ckplm(k,1083) = 1.012451171875D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(20592.D0+(rk*rk+rk-398.D0)*(rk+1.D0)*rk) * den(7)
        ckplm(k,1084) = 0.009735107421875D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(2.079792D6+rk*(121354.D0+rk*(-37709.D0+rk*(-1366.D0+89.D0*rk)))) * den(8)
        ckplm(k,1085) = 0.00059000651041667D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(2.8949184D7+rk*(4.11621D6+rk*(-386045.D0+rk*(-41790.D0+761.D0*rk)))) * den(9)
        ckplm(k,1086) = -0.00560506184895833D0*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-2.18544D6+rk*(-533026.D0+rk*(4781.D0+rk*(5014.D0+31.D0*rk)))) * den(10)
        ckplm(k,1087) = -0.0672607421875D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-26418.D0+rk*(-3537.D0+rk*(262.D0+13.D0*rk))) * den(11)
        ckplm(k,1088) = -0.128916422526042D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1938.D0+rk*(-47.D0+11.D0*rk)) * den(12)
        ckplm(k,1089) = -1.28916422526042D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-23.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,1090) = 1.93374633789063D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1091) = -0.101776123046875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1092) = 0.0339253743489583D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(75.D0+4.D0*rk) * den(1)
        ckplm(k,1093) = 0.00135701497395833D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-21875.D0+rk*(-678.D0+53.D0*rk)) * den(2)
        ckplm(k,1094) = -0.00177001953125D0*(rk+10.D0)*(rk+11.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-140760.D0+rk*(7031.D0+rk*(1471.D0+8.D0*rk))) * den(3)
        ckplm(k,1095) = -0.00029500325520833D0*(rk+10.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(6.24204D6+rk*(-930204.D0+rk*(-82741.D0+rk*(6396.D0+229.D0*rk)))) * den(4)
        ckplm(k,1096) = -0.00147501627604167D0*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-8.831688D6+rk*(2.15245D6+rk*(61695.D0+rk*(-28646.D0+rk*(-111.D0+52.D0*rk))))) * den(5)
        ckplm(k,1097) = -0.048675537109375D0*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk-9.D0)&
                *(41448.D0+rk*(-602.D0+rk*(-895.D0+(rk-28.D0)*rk))) * den(6)
        ckplm(k,1098) = 0.5062255859375D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(41184.D0+5.D0*(rk&
                *rk+rk-200.D0)*(rk+1.D0)*rk) * den(7)
        ckplm(k,1099) = 0.048675537109375D0*(rk+10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(41184.D0+rk*(-1100.D0+rk*(-805.D0+(rk+32.D0)*rk))) * den(8)
        ckplm(k,1100) = 0.00147501627604167D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(1.089396D7+rk&
                *(1.943826D6+rk*(-146447.D0+rk*(-27682.D0+rk*(371.D0+52.D0*rk))))) * den(9)
        ckplm(k,1101) = 0.00029500325520833D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(7.083336D6+rk*(746450.D0+rk*(-100555.D0+rk*(-5480.D0+229.D0*rk)))) * den(10)
        ckplm(k,1102) = 0.00177001953125D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(146328.D0+rk*(4113.D0+rk*(-1447.D0+8.D0*rk))) * den(11)
        ckplm(k,1103) = -0.00135701497395833D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-21144.D0+rk*(784.D0+53.D0*rk)) * den(12)
        ckplm(k,1104) = -0.0339253743489583D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-71.D0+4.D0*rk) * den(13)
        ckplm(k,1105) = 0.101776123046875D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1106) = 0.00508880615234375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1107) = -0.00037694860387731D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(486.D0+29.D0*rk) * den(1)
        ckplm(k,1108) = -0.00011308458116319D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-26490.D0+(rk-1861.D0)*rk) * den(2)
        ckplm(k,1109) = 0.00001966688368056D0*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.668696D6+rk*(-56234.D0+rk*(9969.D0+311.D0*rk))) * den(3)
        ckplm(k,1110) = 1.63890697337963D-6*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.78711848D8+rk&
                *(-4.93713D6+rk*(-2.521345D6+rk*(-23550.D0+3337.D0*rk)))) * den(4)
        ckplm(k,1111) = 0.00002950032552083D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-7.957512D7+rk&
                *(7.411164D6+rk*(1.5577D6+rk*(-70075.D0+rk*(-7180.D0+31.D0*rk))))) * den(5)
        ckplm(k,1112) = -0.00016225179036458D0*(rk+7.D0)*(rk+8.D0)*(-1.0855944D8+rk*(1.6438668D7+rk&
                *(2.397632D6+rk*(-301335.D0+rk*(-16255.D0+rk*(1107.D0+23.D0*rk)))))) * den(6)
        ckplm(k,1113) = -0.00562472873263889D0*(rk-6.D0)*(rk+7.D0)*(-3.70656D6+(rk+1.D0)*rk&
                *(108552.D0+(rk*rk+rk-818.D0)*(rk+1.D0)*rk)) * den(7)
        ckplm(k,1114) = -0.00016225179036458D0*(rk-6.D0)*(rk-7.D0)*(-1.2231648D8+rk*(-1.0809816D7+rk&
                *(3.193382D6+rk*(225705.D0+rk*(-21445.D0+rk*(-969.D0+23.D0*rk)))))) * den(8)
        ckplm(k,1115) = 0.00002950032552083D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(8.536572D7+rk&
                *(4.114414D6+rk*(-1.724535D6+rk*(-41045.D0+rk*(7335.D0+31.D0*rk))))) * den(9)
        ckplm(k,1116) = 1.63890697337963D-6*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.8115452D8+rk&
                *(-21562.D0+rk*(-2.430673D6+rk*(36898.D0+3337.D0*rk)))) * den(10)
        ckplm(k,1117) = 0.00001966688368056D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(1.602804D6+rk*(-75239.D0+rk*(-9036.D0+311.D0*rk))) * den(11)
        ckplm(k,1118) = -0.00011308458116319D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-24628.D0+(rk+1863.D0)*rk) * den(12)
        ckplm(k,1119) = -0.00037694860387731D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-457.D0+29.D0*rk) * den(13)
        ckplm(k,1120) = 0.00508880615234375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1121) = -0.00024232410249256D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1122) = 0.00018847430193866D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(63.D0+4.D0*rk) * den(1)
        ckplm(k,1123) = -0.00003769486038773D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(6895.D0+rk*(654.D0+11.D0*rk)) * den(2)
        ckplm(k,1124) = -3.27781394675926D-6*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.108002D6+rk*(-100643.D0+rk*(993.D0+152.D0*rk))) * den(3)
        ckplm(k,1125) = 0.00001147234881366D0*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(-3.386868D6+rk&
                *(-210232.D0+rk*(27263.D0+(rk+1708.D0)*rk))) * den(4)
        ckplm(k,1126) = 0.00010325113932292D0*(rk+8.D0)*(rk+9.D0)*(3.41252D6+rk*(82566.D0+rk&
                *(-57285.D0+rk*(-2110.D0+rk*(145.D0+4.D0*rk))))) * den(5)
        ckplm(k,1127) = 0.00037858751085069D0*(rk+8.D0)*(-7.52004D6+rk*(80778.D0+rk*(190189.D0+rk&
                *(1320.D0+rk*(-1190.D0+(rk-18.D0)*rk))))) * den(6)
        ckplm(k,1128) = -0.00281236436631944D0*(-7.41312D6+7.D0*(rk+1.D0)*rk*(36372.D0+(rk&
                *rk+rk-368.D0)*(rk+1.D0)*rk)) * den(7)
        ckplm(k,1129) = -0.00037858751085069D0*(rk-7.D0)*(-7.41312D6+rk*(290976.D0+rk*(179284.D0+rk&
                *(-5880.D0+rk*(-1085.D0+(rk+24.D0)*rk))))) * den(8)
        ckplm(k,1130) = -0.00010325113932292D0*(rk-7.D0)*(rk-8.D0)*(-3.27492D6+rk*(190246.D0+rk&
                *(50125.D0+rk*(-2650.D0+rk*(-125.D0+4.D0*rk))))) * den(9)
        ckplm(k,1131) = -0.00001147234881366D0*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-3.15108D6+rk&
                *(259638.D0+rk*(22145.D0+(rk-1704.D0)*rk))) * den(10)
        ckplm(k,1132) = 3.27781394675926D-6*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.006518D6+rk&
                *(-102173.D0+rk*(-537.D0+152.D0*rk))) * den(11)
        ckplm(k,1133) = 0.00003769486038773D0*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(6252.D0+rk*(-632.D0+11.D0*rk)) * den(12)
        ckplm(k,1134) = -0.00018847430193866D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-59.D0+4.D0*rk) * den(13)
        ckplm(k,1135) = 0.00024232410249256D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1136) = 0.00001101473193148D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,1137) = -2.44771820699555D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(288.D0+19.D0*rk) * den(1)
        ckplm(k,1138) = 2.44771820699555D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)*(81440.D0+rk&
                *(9057.D0+223.D0*rk)) * den(2)
        ckplm(k,1139) = 4.25690122955748D-8*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(-8.171808D6+rk&
                *(-1.072522D6+rk*(-29553.D0+283.D0*rk))) * den(3)
        ckplm(k,1140) = -1.17064783812831D-7*(rk+10.D0)*(rk+9.D0)*(-3.7799328D7+rk*(-4.915082D6+rk&
                *(-5987.D0+rk*(15158.D0+311.D0*rk)))) * den(4)
        ckplm(k,1141) = -2.10716610863095D-6*(rk+9.D0)*(2.150464D7+rk*(2.536152D6+rk*(-145130.D0+rk&
                *(-22585.D0+rk*(-290.D0+13.D0*rk))))) * den(5)
        ckplm(k,1142) = 3.51194351438492D-7*(1.12280256D9+rk*(1.18967208D8+rk*(-1.7079646D7+rk&
                *(-1.942665D6+rk*(31805.D0+rk*(5337.D0+41.D0*rk)))))) * den(6)
        ckplm(k,1143) = 0.0000365242125496D0*(-1.019304D7+(rk+1.D0)*rk*(220188.D0+(rk*rk+rk-1196.D0)&
                *(rk+1.D0)*rk)) * den(7)
        ckplm(k,1144) = 3.51194351438492D-7*(9.8872488D8+rk*(-1.47197724D8+rk*(-1.1113576D7+rk&
                *(2.017335D6+rk*(5735.D0+rk*(-5091.D0+41.D0*rk)))))) * den(8)
        ckplm(k,1145) = -2.10716610863095D-6*(rk-8.D0)*(-1.884564D7+rk*(2.759882D6+rk*(79245.D0+rk&
                *(-21295.D0+rk*(355.D0+13.D0*rk))))) * den(9)
        ckplm(k,1146) = -1.17064783812831D-7*(rk-8.D0)*(rk-9.D0)*(-3.290508D7+rk*(4.858878D6+rk&
                *(-49595.D0+rk*(-13914.D0+311.D0*rk)))) * den(10)
        ckplm(k,1147) = 4.25690122955748D-8*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(7.129122D6+rk&
                *(-1.012567D6+rk*(30402.D0+283.D0*rk))) * den(11)
        ckplm(k,1148) = 2.44771820699555D-7*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)*(72606.D0+rk&
                *(-8611.D0+223.D0*rk)) * den(12)
        ckplm(k,1149) = -2.44771820699555D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-269.D0+19.D0*rk) * den(13)
        ckplm(k,1150) = 0.00001101473193148D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(14)
        ckplm(k,1151) = -4.78901388325216D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0) &
                * den(0)
        ckplm(k,1152) = 1.77370884564895D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(2187.D0+148.D0*rk) * den(1)
        ckplm(k,1153) = -1.06422530738937D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(129465.D0+rk&
                *(15874.D0+461.D0*rk)) * den(2)
        ckplm(k,1154) = 2.12845061477874D-8*(rk+10.D0)*(rk+11.D0)*(1.384506D6+rk*(222499.D0+rk&
                *(10191.D0+104.D0*rk))) * den(3)
        ckplm(k,1155) = 3.90215946042769D-8*(rk+10.D0)*(-1.1324124D7+rk*(-2.067D6+rk*(-96655.D0+rk&
                *(660.D0+79.D0*rk)))) * den(4)
        ckplm(k,1156) = -3.51194351438492D-7*(-1.451628D7+rk*(-2.818794D6+rk*(-93125.D0+rk&
                *(9410.D0+rk*(545.D0+4.D0*rk))))) * den(5)
        ckplm(k,1157) = -3.16074916294643D-6*(1.67948D6+rk*(157074.D0+rk*(-13645.D0+rk&
                *(-1315.D0+(rk+5.D0)*rk)))) * den(6)
        ckplm(k,1158) = 0.0000547863188244D0*(94380.D0+(rk+1.D0)*rk*(-1216.D0+3.D0*(rk+1.D0)*rk)) &
                * den(7)
        ckplm(k,1159) = 3.16074916294643D-6*(-1.51008D6+rk*(180404.D0+(rk*rk*rk-1325.D0*rk+9680.D0)&
                *rk)) * den(8)
        ckplm(k,1160) = 3.51194351438492D-7*(1.179948D7+rk*(-2.606474D6+rk*(118125.D0+rk*(7270.D0+rk&
                *(-525.D0+4.D0*rk))))) * den(9)
        ckplm(k,1161) = -3.90215946042769D-8*(rk-9.D0)*(-9.35436D6+rk*(1.872026D6+rk*(-98161.D0+rk&
                *(-344.D0+79.D0*rk)))) * den(10)
        ckplm(k,1162) = -2.12845061477874D-8*(rk-10.D0)*(rk-9.D0)*(-1.172094D6+rk*(202429.D0+rk&
                *(-9879.D0+104.D0*rk))) * den(11)
        ckplm(k,1163) = 1.06422530738937D-8*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(114052.D0+rk&
                *(-14952.D0+461.D0*rk)) * den(12)
        ckplm(k,1164) = -1.77370884564895D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(-2039.D0+148.D0*rk) * den(13)
        ckplm(k,1165) = 4.78901388325216D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0) &
                * den(14)
        ckplm(k,1166) = 1.99542245135507D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0) * den(0)
        ckplm(k,1167) = -4.43427211412237D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(450.D0+31.D0*rk) &
                * den(1)
        ckplm(k,1168) = 8.86854422824475D-11*(rk+11.D0)*(rk+12.D0)*(981250.D0+rk*(128409.D0+4091.D0&
                *rk)) * den(2)
        ckplm(k,1169) = -1.77370884564895D-9*(rk+11.D0)*(1.26516D6+rk*(231206.D0+rk*(13161.D0+223.D0&
                *rk))) * den(3)
        ckplm(k,1170) = -4.87769932553461D-9*(-8.06316D6+rk*(-1.804186D6+rk*(-131449.D0+(rk-3086.D0)&
                *rk))) * den(4)
        ckplm(k,1171) = 2.92661959532077D-8*(-1.741896D6+rk*(-273550.D0+rk*(-8435.D0+rk*(310.D0+11.D0&
                *rk)))) * den(5)
        ckplm(k,1172) = -1.46330979766038D-8*(-3.90588D6+rk*(-305354.D0+rk*(15191.D0+(rk+1106.D0)&
                *rk))) * den(6)
        ckplm(k,1173) = -3.0436843791336D-7*(188760.D0+(rk*rk+rk-1322.D0)*(rk+1.D0)*rk) * den(7)
        ckplm(k,1174) = -1.46330979766038D-8*(-3.58644D6+rk*(332422.D0+rk*(11879.D0+(rk-1102.D0)&
                *rk))) * den(8)
        ckplm(k,1175) = 2.92661959532077D-8*(-1.47708D6+rk*(255794.D0+rk*(-9299.D0+rk*(-266.D0+11.D0&
                *rk)))) * den(9)
        ckplm(k,1176) = -4.87769932553461D-9*(-6.387336D6+rk*(1.55055D6+rk*(-122185.D0+(rk+3090.D0)&
                *rk))) * den(10)
        ckplm(k,1177) = -1.77370884564895D-9*(rk-10.D0)*(-1.046892D6+rk*(205553.D0+rk&
                *(-12492.D0+223.D0*rk))) * den(11)
        ckplm(k,1178) = 8.86854422824475D-11*(rk-10.D0)*(rk-11.D0)*(856932.D0+rk*(-120227.D0+4091.D0&
                *rk)) * den(12)
        ckplm(k,1179) = -4.43427211412237D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(-419.D0+31.D0*rk) &
                * den(13)
        ckplm(k,1180) = 1.99542245135507D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0) * den(14)
        ckplm(k,1181) = -7.98168980542027D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0) * den(0)
        ckplm(k,1182) = 8.86854422824475D-11*(rk+12.D0)*(rk+13.D0)*(1089.D0+76.D0*rk) * den(1)
        ckplm(k,1183) = -8.86854422824475D-11*(rk+12.D0)*(57233.D0+rk*(7842.D0+265.D0*rk)) * den(2)
        ckplm(k,1184) = 1.77370884564895D-10*(877008.D0+rk*(175067.D0+rk*(11283.D0+232.D0*rk))) &
                * den(3)
        ckplm(k,1185) = -9.75539865106923D-10*(294504.D0+rk*(48796.D0+rk*(2399.D0+31.D0*rk))) * den(4)
        ckplm(k,1186) = -2.92661959532077D-9*(-142524.D0+rk*(-16891.D0+rk*(-399.D0+4.D0*rk))) * den(5)
        ckplm(k,1187) = 3.21928155485284D-8*(-15636.D0+rk*(-949.D0+(rk+24.D0)*rk)) * den(6)
        ckplm(k,1188) = -1.67402640852348D-6*(rk*rk+rk-312.D0) * den(7)
        ckplm(k,1189) = -3.21928155485284D-8*(14664.D0+rk*(-994.D0+(rk-21.D0)*rk)) * den(8)
        ckplm(k,1190) = 2.92661959532077D-9*(126036.D0+rk*(-16081.D0+rk*(411.D0+4.D0*rk))) * den(9)
        ckplm(k,1191) = 9.75539865106923D-10*(-248076.D0+rk*(44091.D0+rk*(-2306.D0+31.D0*rk))) &
                * den(10)
        ckplm(k,1192) = -1.77370884564895D-10*(-712992.D0+rk*(153197.D0+rk*(-10587.D0+232.D0*rk))) &
                * den(11)
        ckplm(k,1193) = 8.86854422824475D-11*(rk-11.D0)*(49656.D0+rk*(-7312.D0+265.D0*rk)) * den(12)
        ckplm(k,1194) = -8.86854422824475D-11*(rk-11.D0)*(rk-12.D0)*(-1013.D0+76.D0*rk) * den(13)
        ckplm(k,1195) = 7.98168980542027D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0) * den(14)
        ckplm(k,1196) = 3.06988069439241D-11*(rk+13.D0)*(rk+14.D0) * den(0)
        ckplm(k,1197) = -2.27398569954994D-12*(rk+13.D0)*(1944.D0+137.D0*rk) * den(1)
        ckplm(k,1198) = 4.43427211412238D-11*(6216.D0+rk*(881.D0+31.D0*rk)) * den(2)
        ckplm(k,1199) = -1.77370884564895D-10*(4658.D0+rk*(603.D0+19.D0*rk)) * den(3)
        ckplm(k,1200) = 1.6258997751782D-10*(10794.D0+rk*(1175.D0+29.D0*rk)) * den(4)
        ckplm(k,1201) = -2.92661959532077D-9*(972.D0+(rk+77.D0)*rk) * den(5)
        ckplm(k,1202) = -1.46330979766038D-9*(-2524.D0+(rk-105.D0)*rk) * den(6)
        ckplm(k,1203) = 3.90215946042769D-9*(rk*rk+rk-1014.D0) * den(7)
        ckplm(k,1204) = -1.46330979766038D-9*(-2418.D0+(rk+107.D0)*rk) * den(8)
        ckplm(k,1205) = -2.92661959532077D-9*(896.D0+(rk-75.D0)*rk) * den(9)
        ckplm(k,1206) = 1.6258997751782D-10*(9648.D0+rk*(-1117.D0+29.D0*rk)) * den(10)
        ckplm(k,1207) = -1.77370884564895D-10*(4074.D0+rk*(-565.D0+19.D0*rk)) * den(11)
        ckplm(k,1208) = 4.43427211412238D-11*(5366.D0+rk*(-819.D0+31.D0*rk)) * den(12)
        ckplm(k,1209) = -2.27398569954994D-12*(rk-12.D0)*(-1807.D0+137.D0*rk) * den(13)
        ckplm(k,1210) = 3.06988069439241D-11*(rk-12.D0)*(rk-13.D0) * den(14)
        ckplm(k,1211) = -1.13699284977497D-12*(rk+14.D0) * den(0)
        ckplm(k,1212) = 1.13699284977497D-12*(169.D0+12.D0*rk) * den(1)
        ckplm(k,1213) = -1.47809070470746D-11*(73.D0+5.D0*rk) * den(2)
        ckplm(k,1214) = 2.95618140941492D-11*(127.D0+8.D0*rk) * den(3)
        ckplm(k,1215) = -1.6258997751782D-10*(56.D0+3.D0*rk) * den(4)
        ckplm(k,1216) = 1.6258997751782D-10*(101.D0+4.D0*rk) * den(5)
        ckplm(k,1217) = -4.87769932553461D-10*(rk+47.D0) * den(6)
        ckplm(k,1218) = 2.536403649278D-8 * den(7)
        ckplm(k,1219) = 4.87769932553461D-10*(rk-46.D0) * den(8)
        ckplm(k,1220) = -1.6258997751782D-10*(-97.D0+4.D0*rk) * den(9)
        ckplm(k,1221) = 1.6258997751782D-10*(-53.D0+3.D0*rk) * den(10)
        ckplm(k,1222) = -2.95618140941492D-11*(-119.D0+8.D0*rk) * den(11)
        ckplm(k,1223) = 1.47809070470746D-11*(-68.D0+5.D0*rk) * den(12)
        ckplm(k,1224) = -1.13699284977497D-12*(-157.D0+12.D0*rk) * den(13)
        ckplm(k,1225) = 1.13699284977497D-12*(rk-13.D0) * den(14)
        ckplm(k,1226) = 4.06068874919631D-14 * den(0)
        ckplm(k,1227) = -5.68496424887484D-13 * den(1)
        ckplm(k,1228) = 3.69522676176865D-12 * den(2)
        ckplm(k,1229) = -1.47809070470746D-11 * den(3)
        ckplm(k,1230) = 4.06474943794551D-11 * den(4)
        ckplm(k,1231) = -8.12949887589102D-11 * den(5)
        ckplm(k,1232) = 1.21942483138365D-10 * den(6)
        ckplm(k,1233) = -1.39362837872418D-10 * den(7)
        ckplm(k,1234) = 1.21942483138365D-10 * den(8)
        ckplm(k,1235) = -8.12949887589102D-11 * den(9)
        ckplm(k,1236) = 4.06474943794551D-11 * den(10)
        ckplm(k,1237) = -1.47809070470746D-11 * den(11)
        ckplm(k,1238) = 3.69522676176865D-12 * den(12)
        ckplm(k,1239) = -5.68496424887484D-13 * den(13)
        ckplm(k,1240) = 4.06068874919631D-14 * den(14)
!    ckplm para l = 15
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+3.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+3.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k-9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k-7.D0)*(r2k-9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-3.D0)*(r2k-5.D0)&
                *(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,1241) = 146748.142089844D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,1242) = 75904.2114257813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk &
                * den(1)
        ckplm(k,1243) = 59036.6088867188D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk &
                * den(2)
        ckplm(k,1244) = 51165.0610351563D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk &
                * den(3)
        ckplm(k,1245) = 46715.9252929688D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk &
                * den(4)
        ckplm(k,1246) = 44046.4438476563D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk &
                * den(5)
        ckplm(k,1247) = 42500.9545898438D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk &
                * den(6)
        ckplm(k,1248) = 41786.6528320313D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*rk &
                * den(7)
        ckplm(k,1249) = 41786.6528320313D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*rk &
                * den(8)
        ckplm(k,1250) = 42500.9545898438D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*rk &
                * den(9)
        ckplm(k,1251) = 44046.4438476563D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk &
                * den(10)
        ckplm(k,1252) = 46715.9252929688D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk &
                * den(11)
        ckplm(k,1253) = 51165.0610351563D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk &
                * den(12)
        ckplm(k,1254) = 59036.6088867188D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk &
                * den(13)
        ckplm(k,1255) = 75904.2114257813D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk &
                * den(14)
        ckplm(k,1256) = 146748.142089844D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk &
                * den(15)
        ckplm(k,1257) = -18343.5177612305D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,1258) = -632.535095214844D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-29.D0+13.D0&
                *rk) * den(1)
        ckplm(k,1259) = -491.971740722656D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-54.D0+11.D0&
                *rk) * den(2)
        ckplm(k,1260) = -1279.12652587891D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-25.D0+3.D0&
                *rk) * den(3)
        ckplm(k,1261) = -389.299377441406D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-92.D0+7.D0&
                *rk) * den(4)
        ckplm(k,1262) = -1835.26849365234D0*(rk+10.D0)*(rk-1.D0)*(rk-21.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(5)
        ckplm(k,1263) = -1062.52386474609D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-38.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(6)
        ckplm(k,1264) = -348.222106933594D0*(rk-119.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0) &
                * den(7)
        ckplm(k,1265) = 348.222106933594D0*(rk+120.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0) &
                * den(8)
        ckplm(k,1266) = 1062.52386474609D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+39.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0) &
                * den(9)
        ckplm(k,1267) = 1835.26849365234D0*(rk-1.D0)*(rk+22.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(10)
        ckplm(k,1268) = 389.299377441406D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(99.D0+7.D0&
                *rk) * den(11)
        ckplm(k,1269) = 1279.12652587891D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(28.D0+3.D0&
                *rk) * den(12)
        ckplm(k,1270) = 491.971740722656D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(65.D0+11.D0&
                *rk) * den(13)
        ckplm(k,1271) = 632.535095214844D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(42.D0+13.D0&
                *rk) * den(14)
        ckplm(k,1272) = 18343.5177612305D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(15)
        ckplm(k,1273) = 1079.03045654297D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1274) = 37.2079467773438D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-116.D0+7.D0*rk) &
                * den(1)
        ckplm(k,1275) = 4.13421630859375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(2862.D0+(rk-1297.D0)*rk) * den(2)
        ckplm(k,1276) = -10.7489624023438D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1850.D0+rk*(487.D0+13.D0*rk)) &
                * den(3)
        ckplm(k,1277) = -3.27142333984375D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-8372.D0+rk*(1401.D0+71.D0*rk)) &
                * den(4)
        ckplm(k,1278) = -15.4224243164063D0*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2184.D0+rk*(233.D0+19.D0*rk)) &
                * den(5)
        ckplm(k,1279) = -8.92877197265625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-4294.D0+rk*(267.D0+37.D0*rk)) &
                * den(6)
        ckplm(k,1280) = -348.222106933594D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(-118.D0+(rk+3.D0)*rk) * den(7)
        ckplm(k,1281) = -348.222106933594D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(-120.D0+(rk-1.D0)*rk) * den(8)
        ckplm(k,1282) = -8.92877197265625D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(-4524.D0+rk*(-193.D0+37.D0*rk)) &
                * den(9)
        ckplm(k,1283) = -15.4224243164063D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2398.D0+rk*(-195.D0+19.D0*rk)) &
                * den(10)
        ckplm(k,1284) = -3.27142333984375D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-9702.D0+rk*(-1259.D0+71.D0*rk)) &
                * den(11)
        ckplm(k,1285) = -10.7489624023438D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2324.D0+rk*(-461.D0+13.D0*rk)) &
                * den(12)
        ckplm(k,1286) = 4.13421630859375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(4160.D0+(rk+1299.D0)*rk) * den(13)
        ckplm(k,1287) = 37.2079467773438D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(123.D0+7.D0*rk) &
                * den(14)
        ckplm(k,1288) = 1079.03045654297D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1289) = -59.9461364746094D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1290) = 6.20132446289063D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+87.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,1291) = 2.06710815429688D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1278.D0+rk*(223.D0+11.D0*rk)) * den(2)
        ckplm(k,1292) = 0.137807210286458D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(89550.D0+rk*(-35381.D0+rk*(1434.D0+197.D0*rk))) &
                * den(3)
        ckplm(k,1293) = 0.125823974609375D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(166152.D0+rk*(-41822.D0+rk*(-333.D0+203.D0*rk))) &
                * den(4)
        ckplm(k,1294) = 0.197723388671875D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(148722.D0+rk*(-23881.D0+rk*(-1296.D0+103.D0*rk))) &
                * den(5)
        ckplm(k,1295) = 0.0381571451822917D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(953724.D0+rk*(-89264.D0+rk*(-10941.D0+341.D0*rk))) &
                * den(6)
        ckplm(k,1296) = 4.46438598632813D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(9126.D0+rk*(-349.D0+(rk-114.D0)*rk)) * den(7)
        ckplm(k,1297) = -4.46438598632813D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(-9360.D0+(rk+118.D0)*(rk-1.D0)*rk) * den(8)
        ckplm(k,1298) = -0.0381571451822917D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(-1.031706D6+rk*(-66359.D0+rk*(11964.D0+341.D0*rk))) &
                * den(9)
        ckplm(k,1299) = -0.197723388671875D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-171204.D0+rk*(-20980.D0+rk*(1605.D0+103.D0*rk))) &
                * den(10)
        ckplm(k,1300) = -0.125823974609375D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-207438.D0+rk*(-40547.D0+rk*(942.D0+203.D0*rk))) &
                * den(11)
        ckplm(k,1301) = -0.137807210286458D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-126168.D0+rk*(-37658.D0+rk*(-843.D0+197.D0*rk))) &
                * den(12)
        ckplm(k,1302) = -2.06710815429688D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1490.D0+rk*(-201.D0+11.D0*rk)) * den(13)
        ckplm(k,1303) = -6.20132446289063D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-86.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1304) = 59.9461364746094D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1305) = 3.15505981445313D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1306) = -0.108795166015625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(464.D0+17.D0*rk) * den(1)
        ckplm(k,1307) = -0.0362650553385417D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-11016.D0+rk*(515.D0+61.D0*rk)) * den(2)
        ckplm(k,1308) = -0.00725301106770833D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(352400.D0+rk*(-67378.D0+rk*(-2633.D0+211.D0*rk))) * den(3)
        ckplm(k,1309) = -0.006622314453125D0*(rk+10.D0)*(rk+11.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-2.407272D6+rk*(808734.D0+rk*(-28159.D0+rk*(-6386.D0+83.D0*rk)))) * den(4)
        ckplm(k,1310) = 0.00346883138020833D0*(rk+10.D0)*(rk+12.D0)*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(616182.D0+rk*(-183677.D0+rk*(11742.D0+113.D0*rk))) * den(5)
        ckplm(k,1311) = 0.0381571451822917D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(904872.D0+rk*(-113366.D0+rk*(-12089.D0+rk*(914.D0+29.D0*rk)))) * den(6)
        ckplm(k,1312) = 1.48812866210938D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(27144.D0+rk*(-1390.D0+rk*(-445.D0+(rk+10.D0)*rk))) * den(7)
        ckplm(k,1313) = 1.48812866210938D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(28080.D0+(rk-1.D0)*rk*(-474.D0+(rk-5.D0)*rk)) * den(8)
        ckplm(k,1314) = 0.0381571451822917D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(1.005264D6+rk*(86562.D0+rk*(-14657.D0+rk*(-798.D0+29.D0*rk)))) * den(9)
        ckplm(k,1315) = 0.00346883138020833D0*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-811488.D0+rk*(-206822.D0+rk*(-11403.D0+113.D0*rk))) * den(10)
        ckplm(k,1316) = -0.006622314453125D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-3.237696D6+rk*(-845562.D0+rk*(-8503.D0+rk*(6718.D0+83.D0*rk)))) * den(11)
        ckplm(k,1317) = -0.00725301106770833D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-416934.D0+rk*(-61479.D0+rk*(3266.D0+211.D0*rk))) * den(12)
        ckplm(k,1318) = -0.0362650553385417D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-11470.D0+rk*(-393.D0+61.D0*rk)) * den(13)
        ckplm(k,1319) = -0.108795166015625D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-447.D0+17.D0*rk) * den(14)
        ckplm(k,1320) = 3.15505981445313D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1321) = -0.157752990722656D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1322) = 0.0271987915039063D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(145.D0+7.D0*rk) * den(1)
        ckplm(k,1323) = 0.00906626383463541D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-5130.D0+rk*(-103.D0+13.D0*rk)) * den(2)
        ckplm(k,1324) = 0.00181325276692708D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(219950.D0+rk*(-14473.D0+(rk-2078.D0)*rk)) * den(3)
        ckplm(k,1325) = -0.00165557861328125D0*(rk+10.D0)*(rk+11.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(1.8354D6+rk*(-310794.D0+rk*(-17861.D0+rk*(2006.D0+49.D0*rk)))) * den(4)
        ckplm(k,1326) = -0.00001576741536458D0*(rk+10.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.4178906D9+rk*(3.81407274D8+rk*(-263425.D0+rk*(-4.264715D6+rk*(42865.D0+7001.D0*rk))))) &
                * den(5)
        ckplm(k,1327) = -0.00086720784505208D0*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-3.775992D7+rk*(5.934772D6+rk*(537000.D0+rk*(-71195.D0+rk*(-2160.D0+103.D0*rk))))) * den(6)
        ckplm(k,1328) = -0.0338211059570313D0*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(107640.D0+rk*(2866.D0+rk*(-1909.D0+(rk-94.D0)*rk))) * den(7)
        ckplm(k,1329) = 0.0338211059570313D0*(rk+12.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(102960.D0+rk*(-6398.D0+rk*(-1621.D0+(rk+98.D0)*rk))) * den(8)
        ckplm(k,1330) = 0.00086720784505208D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(4.308876D7+rk*(4.656342D6+rk*(-736595.D0+rk*(-61525.D0+rk*(2675.D0+103.D0*rk))))) * den(9)
        ckplm(k,1331) = 0.00001576741536458D0*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(1.79526072D9+rk*(3.69003524D8+rk*(-1.27179D7+rk*(-4.366165D6+rk*(-7860.D0+7001.D0*rk))))) &
                * den(10)
        ckplm(k,1332) = 0.00165557861328125D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(2.126376D6+rk*(269250.D0+rk*(-23585.D0+rk*(-1810.D0+49.D0*rk)))) * den(11)
        ckplm(k,1333) = -0.00181325276692708D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-232344.D0+rk*(-10314.D0+(rk+2081.D0)*rk)) * den(12)
        ckplm(k,1334) = -0.00906626383463541D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-5014.D0+rk*(129.D0+13.D0*rk)) * den(13)
        ckplm(k,1335) = -0.0271987915039063D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-138.D0+7.D0*rk) * den(14)
        ckplm(k,1336) = 0.157752990722656D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1337) = 0.00751204717726934D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1338) = -0.00077710832868304D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(348.D0+19.D0*rk) * den(1)
        ckplm(k,1339) = -0.00181325276692708D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-2466.D0+(rk-145.D0)*rk) * den(2)
        ckplm(k,1340) = 0.00012088351779514D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-414900.D0+rk*(-6632.D0+rk*(2523.D0+59.D0*rk))) * den(3)
        ckplm(k,1341) = 0.00001576741536458D0*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(2.9408904D7+rk*(-1.506978D6+rk*(-372797.D0+rk*(1362.D0+509.D0*rk)))) * den(4)
        ckplm(k,1342) = 0.00001576741536458D0*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.47044672D8+rk*(3.041748D7+rk*(3.80669D6+rk*(-282575.D0+rk*(-16058.D0+215.D0*rk))))) &
                * den(5)
        ckplm(k,1343) = -0.00005781385633681D0*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-5.3696808D8+rk&
                *(1.01608956D8+rk*(7.531568D6+rk*(-1.587075D6+rk*(-33895.D0+rk*(5199.D0+47.D0*rk)))))) * den(6)
        ckplm(k,1344) = -0.00676422119140625D0*(rk-6.D0)*(rk+7.D0)*(rk+8.D0)*(-5.86872D6+rk&
                *(454644.D0+rk*(139504.D0+rk*(-7185.D0+rk*(-905.D0+(rk+21.D0)*rk))))) * den(7)
        ckplm(k,1345) = -0.00676422119140625D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(-6.1776D6+(rk-1.D0)*rk&
                *(157800.D0+rk*(2366.D0+rk*(-1009.D0+(rk-14.D0)*rk)))) * den(8)
        ckplm(k,1346) = -0.00005781385633681D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(-6.2949744D8+rk&
                *(-8.1945888D7+rk*(1.2038138D7+rk*(1.400445D6+rk*(-59185.D0+rk*(-4917.D0+47.D0*rk)))))) * den(9)
        ckplm(k,1347) = 0.00001576741536458D0*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.7338916D8+rk&
                *(2.2021682D7+rk*(-4.555917D6+rk*(-216193.D0+rk*(17133.D0+215.D0*rk))))) * den(10)
        ckplm(k,1348) = 0.00001576741536458D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(3.0542232D7+rk*(759334.D0+rk*(-373829.D0+rk*(674.D0+509.D0*rk)))) * den(11)
        ckplm(k,1349) = 0.00012088351779514D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(405804.D0+rk*(-11501.D0+rk*(-2346.D0+59.D0*rk))) * den(12)
        ckplm(k,1350) = -0.00181325276692708D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-2320.D0+(rk+147.D0)*rk) * den(13)
        ckplm(k,1351) = -0.00077710832868304D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-329.D0+19.D0*rk) * den(14)
        ckplm(k,1352) = 0.00751204717726934D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1353) = -0.00034145668987588D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1354) = 0.00001177436861641D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(1421.D0+83.D0*rk) * den(1)
        ckplm(k,1355) = -9.15784225720749D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(40446.D0+rk*(3359.D0+43.D0*rk)) * den(2)
        ckplm(k,1356) = -1.8315684514415D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.89765D6+rk*(-207707.D0+rk*(5398.D0+359.D0*rk))) * den(3)
        ckplm(k,1357) = -5.57433876525673D-7*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.05954744D8+rk*(3.764138D6+rk*(-909343.D0+rk*(-38222.D0+283.D0*rk)))) * den(4)
        ckplm(k,1358) = 0.00005518595377604D0*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(1.0245192D7+rk&
                *(-105506.D0+rk*(-162119.D0+rk*(-2269.D0+rk*(439.D0+7.D0*rk))))) * den(5)
        ckplm(k,1359) = 0.00001839531792535D0*(rk+8.D0)*(rk+9.D0)*(-2.6655816D8+rk*(1.4594412D7+rk&
                *(5.830016D6+rk*(-144435.D0+rk*(-33925.D0+rk*(63.D0+29.D0*rk)))))) * den(6)
        ckplm(k,1360) = 0.00003416273328993D0*(rk+8.D0)*(1.15181352D9+7.D0*rk*(-1.4935092D7+rk&
                *(-4.481108D6+rk*(300409.D0+rk*(38080.D0+rk*(-1478.D0+(rk-92.D0)*rk)))))) * den(7)
        ckplm(k,1361) = -0.00003416273328993D0*(rk-7.D0)*(-1.2231648D9+7.D0*(rk-1.D0)*rk*(5.2308D6+rk&
                *(90324.D0+rk*(-44860.D0+rk*(-805.D0+(rk+100.D0)*rk))))) * den(8)
        ckplm(k,1362) = -0.00001839531792535D0*(rk-7.D0)*(rk-8.D0)*(-2.7521208D8+rk*(-2.636916D6+rk&
                *(6.059576D6+rk*(8685.D0+rk*(-33805.D0+rk*(111.D0+29.D0*rk)))))) * den(9)
        ckplm(k,1363) = -0.00005518595377604D0*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1.019128D7+rk&
                *(210204.D0+rk*(152748.D0+rk*(-3955.D0+rk*(-404.D0+7.D0*rk))))) * den(10)
        ckplm(k,1364) = 5.57433876525673D-7*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.01319768D8+rk&
                *(-5.467026D6+rk*(-792979.D0+rk*(39354.D0+283.D0*rk)))) * den(11)
        ckplm(k,1365) = 1.8315684514415D-6*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(2.684904D6+rk*(-217426.D0+rk*(-4321.D0+359.D0*rk))) * den(12)
        ckplm(k,1366) = 9.15784225720749D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(37130.D0+rk*(-3273.D0+43.D0*rk)) * den(13)
        ckplm(k,1367) = -0.00001177436861641D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-1338.D0+83.D0*rk) * den(14)
        ckplm(k,1368) = 0.00034145668987588D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1369) = 0.00001484594303808D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+9.D0) * den(0)
        ckplm(k,1370) = -5.1192907027868D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0)*(1856.D0+113.D0*rk) * den(1)
        ckplm(k,1371) = 3.98167054661195D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(68256.D0+rk*(6769.D0+143.D0*rk)) * den(2)
        ckplm(k,1372) = 7.96334109322391D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)*(-6.1136D6+rk&
                *(-678558.D0+rk*(-12163.D0+321.D0*rk))) * den(3)
        ckplm(k,1373) = -7.96334109322391D-8*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(-8.1168864D7+rk&
                *(-8.192398D6+rk*(173903.D0+rk*(29362.D0+397.D0*rk)))) * den(4)
        ckplm(k,1374) = -7.88370768229167D-6*(rk+10.D0)*(rk+9.D0)*(8.921472D6+rk*(706456.D0+rk&
                *(-80634.D0+rk*(-6929.D0+rk*(14.D0+5.D0*rk))))) * den(5)
        ckplm(k,1375) = -2.62790256076389D-6*(rk+9.D0)*(-2.5253184D8+rk*(-1.4095512D7+rk&
                *(4.265554D6+rk*(256095.D0+rk*(-13835.D0+(rk-783.D0)*rk))))) * den(6)
        ckplm(k,1376) = 0.00003416273328993D0*(-1.6308864D8+rk*(-6.308592D6+rk*(4.078972D6+rk&
                *(167704.D0+rk*(-27935.D0+rk*(-1073.D0+(rk+43.D0)*rk)))))) * den(7)
        ckplm(k,1377) = 0.00003416273328993D0*(1.528956D8+rk*(-1.38573D7+rk*(-3.419604D6+rk&
                *(267889.D0+rk*(21960.D0+rk*(-1310.D0+(rk-36.D0)*rk)))))) * den(8)
        ckplm(k,1378) = -2.62790256076389D-6*(rk-8.D0)*(-2.3443992D8+rk*(2.1806916D7+rk&
                *(3.422104D6+rk*(-303585.D0+rk*(-9905.D0+(rk+789.D0)*rk))))) * den(9)
        ckplm(k,1379) = -7.88370768229167D-6*(rk-8.D0)*(rk-9.D0)*(-8.14132D6+rk*(846906.D0+rk&
                *(59813.D0+rk*(-6935.D0+rk*(11.D0+5.D0*rk))))) * den(10)
        ckplm(k,1380) = -7.96334109322391D-8*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(-7.2831528D7+rk&
                *(8.453706D6+rk*(88199.D0+rk*(-27774.D0+397.D0*rk)))) * den(11)
        ckplm(k,1381) = 7.96334109322391D-8*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)*(5.447526D6+rk&
                *(-653269.D0+rk*(13126.D0+321.D0*rk))) * den(12)
        ckplm(k,1382) = 3.98167054661195D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(61630.D0+rk*(-6483.D0+143.D0*rk)) * den(13)
        ckplm(k,1383) = -5.1192907027868D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-1743.D0+113.D0*rk) * den(14)
        ckplm(k,1384) = 0.00001484594303808D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1385) = -6.18580959920071D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0) * den(0)
        ckplm(k,1386) = 6.3991133784835D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(783.D0+49.D0*rk) * den(1)
        ckplm(k,1387) = -7.11012597609277D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(253206.D0+rk*(27955.D0+719.D0*rk)) * den(2)
        ckplm(k,1388) = 4.74008398406185D-10*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(8.345835D7+rk&
                *(1.1650733D7+rk*(436038.D0+2479.D0*rk))) * den(3)
        ckplm(k,1389) = 9.95417636652988D-9*(rk+10.D0)*(rk+11.D0)*(-6.2374536D7+rk*(-9.433038D6+rk&
                *(-292187.D0+rk*(9402.D0+359.D0*rk)))) * den(4)
        ckplm(k,1390) = 4.69268314422123D-8*(rk+10.D0)*(1.63161432D8+rk*(2.4710514D7+rk*(152095.D0+rk&
                *(-108955.D0+(rk-3727.D0)*rk)))) * den(5)
        ckplm(k,1391) = -1.56422771474041D-8*(5.02255512D9+rk*(7.45691916D8+rk*(-2.4259792D7+rk&
                *(-7.124715D6+rk*(-161725.D0+rk*(8439.D0+197.D0*rk)))))) * den(6)
        ckplm(k,1392) = -1.83014642624628D-6*(-4.190472D7+rk*(-1.636356D6+rk*(655744.D0+rk&
                *(24705.D0+rk*(-2345.D0+(rk-69.D0)*rk))))) * den(7)
        ckplm(k,1393) = 1.83014642624628D-6*(-3.96396D7+rk*(2.8647D6+rk*(568264.D0+rk*(-33375.D0+rk&
                *(-1985.D0+(rk+75.D0)*rk))))) * den(8)
        ckplm(k,1394) = 1.56422771474041D-8*(4.25955816D9+rk*(-7.73525268D8+rk*(-3.937432D6+rk&
                *(6.397365D6+rk*(-200965.D0+rk*(-7257.D0+197.D0*rk)))))) * den(9)
        ckplm(k,1395) = -4.69268314422123D-8*(rk-9.D0)*(-1.3870824D8+rk*(2.4094372D7+rk&
                *(-456588.D0+rk*(-94037.D0+(rk+3732.D0)*rk)))) * den(10)
        ckplm(k,1396) = -9.95417636652988D-9*(rk-10.D0)*(rk-9.D0)*(-5.3242728D7+rk*(8.821894D6+rk&
                *(-318239.D0+rk*(-7966.D0+359.D0*rk)))) * den(11)
        ckplm(k,1397) = -4.74008398406185D-10*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(-7.2241176D7+rk&
                *(1.0786094D7+rk*(-428601.D0+2479.D0*rk))) * den(12)
        ckplm(k,1398) = 7.11012597609277D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)*(225970.D0+rk&
                *(-26517.D0+719.D0*rk)) * den(13)
        ckplm(k,1399) = -6.3991133784835D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(-734.D0+49.D0*rk) * den(14)
        ckplm(k,1400) = 6.18580959920071D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0) * den(15)
        ckplm(k,1401) = 2.47432383968028D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0) &
                * den(0)
        ckplm(k,1402) = -4.26607558565566D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(580.D0+37.D0*rk) * den(1)
        ckplm(k,1403) = 4.74008398406185D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(230310.D0+rk&
                *(27299.D0+781.D0*rk)) * den(2)
        ckplm(k,1404) = -4.74008398406185D-10*(rk+11.D0)*(rk+12.D0)*(6.10538D6+rk*(982824.D0+rk&
                *(48049.D0+657.D0*rk))) * den(3)
        ckplm(k,1405) = -3.31805878884329D-9*(rk+11.D0)*(-1.603644D7+rk*(-3.056242D6+rk&
                *(-178213.D0+rk*(-2462.D0+37.D0*rk)))) * den(4)
        ckplm(k,1406) = 3.12845542948082D-9*(-2.3690016D8+rk*(-4.9894008D7+rk*(-2.84335D6+rk&
                *(26005.D0+rk*(5470.D0+83.D0*rk))))) * den(5)
        ckplm(k,1407) = 1.56422771474041D-8*(5.283828D7+rk*(6.721766D6+rk*(-90575.D0+rk*(-32135.D0+rk&
                *(-625.D0+9.D0*rk))))) * den(6)
        ckplm(k,1408) = -2.03349602916253D-7*(4.15272D6+rk*(151684.D0+rk*(-38340.D0+rk&
                *(-1205.D0+(rk+60.D0)*rk)))) * den(7)
        ckplm(k,1409) = -2.03349602916253D-7*(-3.96396D6+rk*(224514.D0+rk*(34375.D0+rk&
                *(-1435.D0+(rk-55.D0)*rk)))) * den(8)
        ckplm(k,1410) = 1.56422771474041D-8*(-4.605744D7+rk*(6.809056D6+rk*(-1990.D0+rk*(-29545.D0+rk&
                *(670.D0+9.D0*rk))))) * den(9)
        ckplm(k,1411) = 3.12845542948082D-9*(1.8987012D8+rk*(-4.4150758D7+rk*(2.889375D6+rk&
                *(4955.D0+rk*(-5055.D0+83.D0*rk))))) * den(10)
        ckplm(k,1412) = -3.31805878884329D-9*(rk-10.D0)*(-1.3155912D7+rk*(2.70735D6+rk*(-170605.D0+rk&
                *(2610.D0+37.D0*rk)))) * den(11)
        ckplm(k,1413) = -4.74008398406185D-10*(rk-10.D0)*(rk-11.D0)*(-5.169948D6+rk*(888697.D0+rk&
                *(-46078.D0+657.D0*rk))) * den(12)
        ckplm(k,1414) = 4.74008398406185D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(203792.D0+rk&
                *(-25737.D0+781.D0*rk)) * den(13)
        ckplm(k,1415) = -4.26607558565566D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(-543.D0+37.D0*rk) * den(14)
        ckplm(k,1416) = 2.47432383968028D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0) &
                * den(15)
        ckplm(k,1417) = -9.51663015261648D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0) * den(0)
        ckplm(k,1418) = 3.28159660435051D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(3509.D0+227.D0*rk) &
                * den(1)
        ckplm(k,1419) = -3.64621844927834D-12*(rk+12.D0)*(rk+13.D0)*(1.679238D6+rk*(209225.D0+6397.D0&
                *rk)) * den(2)
        ckplm(k,1420) = 4.74008398406185D-11*(rk+12.D0)*(4.07165D6+rk*(721993.D0+rk*(40798.D0+719.D0&
                *rk))) * den(3)
        ckplm(k,1421) = -3.31805878884329D-10*(1.2420408D7+rk*(2.755154D6+rk*(209561.D0+rk&
                *(6034.D0+43.D0*rk)))) * den(4)
        ckplm(k,1422) = -1.56422771474041D-9*(-3.754296D6+rk*(-633282.D0+rk*(-29113.D0+rk&
                *(-42.D0+13.D0*rk)))) * den(5)
        ckplm(k,1423) = 1.72065048621445D-8*(-413352.D0+rk*(-42974.D0+rk*(119.D0+(rk+86.D0)*rk))) &
                * den(6)
        ckplm(k,1424) = 1.72065048621445D-8*(442104.D0+rk*(14090.D0+rk*(-2185.D0+(rk-50.D0)*rk))) &
                * den(7)
        ckplm(k,1425) = -1.72065048621445D-8*(425880.D0+rk*(-18306.D0+rk*(-2029.D0+(rk+54.D0)*rk))) &
                * den(8)
        ckplm(k,1426) = -1.72065048621445D-8*(-370344.D0+rk*(42958.D0+rk*(-133.D0+(rk-82.D0)*rk))) &
                * den(9)
        ckplm(k,1427) = 1.56422771474041D-9*(-3.150072D6+rk*(575234.D0+rk*(-28909.D0+rk*(94.D0+13.D0&
                *rk)))) * den(10)
        ckplm(k,1428) = 3.31805878884329D-10*(9.868824D6+rk*(-2.353962D6+rk*(191717.D0+rk&
                *(-5862.D0+43.D0*rk)))) * den(11)
        ckplm(k,1429) = -4.74008398406185D-11*(rk-11.D0)*(-3.389736D6+rk*(642554.D0+rk&
                *(-38641.D0+719.D0*rk))) * den(12)
        ckplm(k,1430) = 3.64621844927834D-12*(rk-11.D0)*(rk-12.D0)*(1.47641D6+rk*(-196431.D0+6397.D0&
                *rk)) * den(13)
        ckplm(k,1431) = -3.28159660435051D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(-3282.D0+227.D0*rk) &
                * den(14)
        ckplm(k,1432) = 9.51663015261648D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0) * den(15)
        ckplm(k,1433) = 3.5246778343024D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0) * den(0)
        ckplm(k,1434) = -3.64621844927834D-12*(rk+13.D0)*(rk+14.D0)*(1392.D0+91.D0*rk) * den(1)
        ckplm(k,1435) = 1.0938655347835D-11*(rk+13.D0)*(29272.D0+rk*(3783.D0+121.D0*rk)) * den(2)
        ckplm(k,1436) = -4.74008398406185D-11*(250000.D0+rk*(47518.D0+rk*(2943.D0+59.D0*rk))) * den(3)
        ckplm(k,1437) = 9.95417636652988D-10*(49.D0+3.D0*rk)*(498.D0+(rk+51.D0)*rk) * den(4)
        ckplm(k,1438) = -5.21409238246803D-10*(74424.D0+rk*(9398.D0+(rk+303.D0)*rk)) * den(5)
        ckplm(k,1439) = 1.73803079415601D-10*(293388.D0+rk*(23257.D0+(48.D0-13.D0*rk)*rk)) * den(6)
        ckplm(k,1440) = 1.56422771474041D-9*(-36504.D0+rk*(-934.D0+(rk+81.D0)*rk)) * den(7)
        ckplm(k,1441) = 1.56422771474041D-9*(35490.D0+rk*(-1093.D0+(rk-78.D0)*rk)) * den(8)
        ckplm(k,1442) = -1.73803079415601D-10*(270192.D0+rk*(-23122.D0+rk*(87.D0+13.D0*rk))) * den(9)
        ckplm(k,1443) = -5.21409238246803D-10*(-65328.D0+rk*(8795.D0+(rk-300.D0)*rk)) * den(10)
        ckplm(k,1444) = 9.95417636652988D-10*(-46.D0+3.D0*rk)*(448.D0+(rk-49.D0)*rk) * den(11)
        ckplm(k,1445) = -4.74008398406185D-11*(-205366.D0+rk*(41809.D0+rk*(-2766.D0+59.D0*rk))) &
                * den(12)
        ckplm(k,1446) = 1.0938655347835D-11*(rk-12.D0)*(25610.D0+rk*(-3541.D0+121.D0*rk)) * den(13)
        ckplm(k,1447) = -3.64621844927834D-12*(rk-12.D0)*(rk-13.D0)*(-1301.D0+91.D0*rk) * den(14)
        ckplm(k,1448) = 3.5246778343024D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0) * den(15)
        ckplm(k,1449) = -1.25881351225086D-12*(rk+14.D0)*(rk+15.D0) * den(0)
        ckplm(k,1450) = 4.34073624914089D-14*(rk+14.D0)*(4901.D0+323.D0*rk) * den(1)
        ckplm(k,1451) = -3.03851537439862D-13*(51714.D0+rk*(6871.D0+227.D0*rk)) * den(2)
        ckplm(k,1452) = 3.95006998671821D-12*(13150.D0+rk*(1621.D0+49.D0*rk)) * den(3)
        ckplm(k,1453) = -3.95006998671821D-12*(30814.D0+rk*(3303.D0+83.D0*rk)) * den(4)
        ckplm(k,1454) = 4.34507698539003D-11*(4998.D0+rk*(419.D0+7.D0*rk)) * den(5)
        ckplm(k,1455) = -4.34507698539003D-11*(7118.D0+(rk+381.D0)*rk) * den(6)
        ckplm(k,1456) = -2.42082860614587D-10*(-1498.D0+(rk-27.D0)*rk) * den(7)
        ckplm(k,1457) = 2.42082860614587D-10*(-1470.D0+(rk+29.D0)*rk) * den(8)
        ckplm(k,1458) = 4.34507698539003D-11*(6738.D0+(rk-379.D0)*rk) * den(9)
        ckplm(k,1459) = -4.34507698539003D-11*(4586.D0+rk*(-405.D0+7.D0*rk)) * den(10)
        ckplm(k,1460) = 3.95006998671821D-12*(27594.D0+rk*(-3137.D0+83.D0*rk)) * den(11)
        ckplm(k,1461) = -3.95006998671821D-12*(11578.D0+rk*(-1523.D0+49.D0*rk)) * den(12)
        ckplm(k,1462) = 3.03851537439862D-13*(45070.D0+rk*(-6417.D0+227.D0*rk)) * den(13)
        ckplm(k,1463) = -4.34073624914089D-14*(rk-13.D0)*(-4578.D0+323.D0*rk) * den(14)
        ckplm(k,1464) = 1.25881351225086D-12*(rk-13.D0)*(rk-14.D0) * den(15)
        ckplm(k,1465) = 4.34073624914089D-14*(rk+15.D0) * den(0)
        ckplm(k,1466) = -4.34073624914089D-14*(196.D0+13.D0*rk) * den(1)
        ckplm(k,1467) = 3.03851537439862D-13*(171.D0+11.D0*rk) * den(2)
        ckplm(k,1468) = -3.95006998671821D-12*(50.D0+3.D0*rk) * den(3)
        ckplm(k,1469) = 2.76504899070275D-11*(rk+19.D0) * den(4)
        ckplm(k,1470) = -4.34507698539003D-11*(rk+24.D0) * den(5)
        ckplm(k,1471) = 4.34507698539003D-11*(rk+37.D0) * den(6)
        ckplm(k,1472) = -1.86217585088144D-11*(rk+106.D0) * den(7)
        ckplm(k,1473) = -1.86217585088144D-11*(rk-105.D0) * den(8)
        ckplm(k,1474) = 4.34507698539003D-11*(rk-36.D0) * den(9)
        ckplm(k,1475) = -4.34507698539003D-11*(rk-23.D0) * den(10)
        ckplm(k,1476) = 2.76504899070275D-11*(rk-18.D0) * den(11)
        ckplm(k,1477) = -3.95006998671821D-12*(-47.D0+3.D0*rk) * den(12)
        ckplm(k,1478) = 3.03851537439862D-13*(-160.D0+11.D0*rk) * den(13)
        ckplm(k,1479) = -4.34073624914089D-14*(-183.D0+13.D0*rk) * den(14)
        ckplm(k,1480) = 4.34073624914089D-14*(rk-14.D0) * den(15)
        ckplm(k,1481) = -1.44691208304696D-15 * den(0)
        ckplm(k,1482) = 2.17036812457044D-14 * den(1)
        ckplm(k,1483) = -1.51925768719931D-13 * den(2)
        ckplm(k,1484) = 6.58344997786368D-13 * den(3)
        ckplm(k,1485) = -1.9750349933591D-12 * den(4)
        ckplm(k,1486) = 4.34507698539003D-12 * den(5)
        ckplm(k,1487) = -7.24179497565005D-12 * den(6)
        ckplm(k,1488) = 9.3108792544072D-12 * den(7)
        ckplm(k,1489) = -9.3108792544072D-12 * den(8)
        ckplm(k,1490) = 7.24179497565005D-12 * den(9)
        ckplm(k,1491) = -4.34507698539003D-12 * den(10)
        ckplm(k,1492) = 1.9750349933591D-12 * den(11)
        ckplm(k,1493) = -6.58344997786368D-13 * den(12)
        ckplm(k,1494) = 1.51925768719931D-13 * den(13)
        ckplm(k,1495) = -2.17036812457044D-14 * den(14)
        ckplm(k,1496) = 1.44691208304696D-15 * den(15)
!    ckplm para l = 16
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)*(r2k+3.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+3.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-3.D0)&
                *(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,1497) = 302668.043060303D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1498) = 156215.764160156D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(1)
        ckplm(k,1499) = 121201.885986328D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(2)
        ckplm(k,1500) = 104742.370605469D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(3)
        ckplm(k,1501) = 95315.5572509766D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(4)
        ckplm(k,1502) = 89513.7407226563D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(5)
        ckplm(k,1503) = 85961.6081542969D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(6)
        ckplm(k,1504) = 84022.6245117188D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*rk * den(7)
        ckplm(k,1505) = 83404.8110961914D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*rk * den(8)
        ckplm(k,1506) = 84022.6245117188D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*rk * den(9)
        ckplm(k,1507) = 85961.6081542969D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(10)
        ckplm(k,1508) = 89513.7407226563D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(11)
        ckplm(k,1509) = 95315.5572509766D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(12)
        ckplm(k,1510) = 104742.370605469D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(13)
        ckplm(k,1511) = 121201.885986328D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(14)
        ckplm(k,1512) = 156215.764160156D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(15)
        ckplm(k,1513) = 302668.043060303D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*rk * den(16)
        ckplm(k,1514) = -35608.005065918D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,1515) = -1148.64532470703D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-31.D0+14.D0*rk) * den(1)
        ckplm(k,1516) = -1782.38067626953D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-29.D0+6.D0*rk) * den(2)
        ckplm(k,1517) = -770.164489746094D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-81.D0+10.D0*rk) * den(3)
        ckplm(k,1518) = -2803.39874267578D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-25.D0+2.D0*rk) * den(4)
        ckplm(k,1519) = -658.189270019531D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-115.D0+6.D0*rk) * den(5)
        ckplm(k,1520) = -1264.14129638672D0*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-63.D0+2.D0*rk) * den(6)
        ckplm(k,1521) = -617.813415527344D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-133.D0+2.D0*rk) * den(7)
        ckplm(k,1522) = 83404.8110961914D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0) &
                * den(8)
        ckplm(k,1523) = 617.813415527344D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(135.D0+2.D0*rk) * den(9)
        ckplm(k,1524) = 1264.14129638672D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(65.D0+2.D0*rk) * den(10)
        ckplm(k,1525) = 658.189270019531D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(121.D0+6.D0*rk) * den(11)
        ckplm(k,1526) = 2803.39874267578D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(27.D0+2.D0*rk) * den(12)
        ckplm(k,1527) = 770.164489746094D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(91.D0+10.D0*rk) * den(13)
        ckplm(k,1528) = 1782.38067626953D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(35.D0+6.D0*rk) * den(14)
        ckplm(k,1529) = 1148.64532470703D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(45.D0+14.D0*rk) * den(15)
        ckplm(k,1530) = 35608.005065918D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(16)
        ckplm(k,1531) = 1978.22250366211D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,1532) = 255.254516601563D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-31.D0+2.D0&
                *rk) * den(1)
        ckplm(k,1533) = 13.2028198242188D0*(4.D0*rk*rk-758.D0*rk+1653.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(2)
        ckplm(k,1534) = -102.688598632813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-360.D0+rk*(97.D0+2.D0&
                *rk)) * den(3)
        ckplm(k,1535) = -186.893249511719D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-275.D0+2.D0*(rk+24.D0)&
                *rk) * den(4)
        ckplm(k,1536) = -48.7547607421875D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1311.D0+rk&
                *(151.D0+10.D0*rk)) * den(5)
        ckplm(k,1537) = -140.460144042969D0*(4.D0*rk*rk+38.D0*rk-525.D0)*(rk+10.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(6)
        ckplm(k,1538) = -27.4583740234375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2926.D0+rk&
                *(111.D0+22.D0*rk)) * den(7)
        ckplm(k,1539) = -617.813415527344D0*(rk*rk+rk-135.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0) * den(8)
        ckplm(k,1540) = -27.4583740234375D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(-3015.D0+rk&
                *(-67.D0+22.D0*rk)) * den(9)
        ckplm(k,1541) = -140.460144042969D0*(4.D0*rk*rk-30.D0*rk-559.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(10)
        ckplm(k,1542) = -48.7547607421875D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1452.D0+rk&
                *(-131.D0+10.D0*rk)) * den(11)
        ckplm(k,1543) = -186.893249511719D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-321.D0+2.D0*(rk-22.D0)&
                *rk) * den(12)
        ckplm(k,1544) = -102.688598632813D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-455.D0+rk*(-93.D0+2.D0&
                *rk)) * den(13)
        ckplm(k,1545) = 13.2028198242188D0*(4.D0*rk*rk+766.D0*rk+2415.D0)*(rk-10.D0)*(rk-11.D0)&
                *(rk-12.D0)*(rk-13.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(14)
        ckplm(k,1546) = 255.254516601563D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(33.D0+2.D0&
                *rk) * den(15)
        ckplm(k,1547) = 1978.22250366211D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(16)
        ckplm(k,1548) = -104.116973876953D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1549) = 3.35861206054688D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(279.D0+2.D0*rk) &
                * den(1)
        ckplm(k,1550) = 1.04232788085938D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-4437.D0+rk*(807.D0+34.D0*rk)) &
                * den(2)
        ckplm(k,1551) = 0.193023681640625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(113310.D0+rk*(-45791.D0+rk*(2199.D0+230.D0&
                *rk))) * den(3)
        ckplm(k,1552) = 0.702606201171875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(53550.D0+rk*(-14051.D0+rk*(39.D0+62.D0&
                *rk))) * den(4)
        ckplm(k,1553) = 0.274932861328125D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(196098.D0+rk*(-33979.D0+rk*(-1253.D0+134.D0&
                *rk))) * den(5)
        ckplm(k,1554) = 0.528045654296875D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(129234.D0+rk*(-14075.D0+rk*(-1221.D0+50.D0&
                *rk))) * den(6)
        ckplm(k,1555) = 6.86459350585938D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(11442.D0+rk*(-653.D0+rk*(-123.D0+2.D0&
                *rk))) * den(7)
        ckplm(k,1556) = -926.720123291016D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk-9.D0) * den(8)
        ckplm(k,1557) = -6.86459350585938D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(-11970.D0+rk*(-401.D0+rk*(129.D0+2.D0&
                *rk))) * den(9)
        ckplm(k,1558) = -0.528045654296875D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-142038.D0+rk*(-11483.D0+rk*(1371.D0+50.D0&
                *rk))) * den(10)
        ckplm(k,1559) = -0.274932861328125D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-228690.D0+rk*(-31071.D0+rk*(1655.D0+134.D0&
                *rk))) * den(11)
        ckplm(k,1560) = -0.702606201171875D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-67578.D0+rk*(-13943.D0+rk*(147.D0+62.D0&
                *rk))) * den(12)
        ckplm(k,1561) = -0.193023681640625D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-161070.D0+rk*(-49499.D0+rk&
                *(-1509.D0+230.D0*rk))) * den(13)
        ckplm(k,1562) = -1.04232788085938D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-5210.D0+rk*(-739.D0+34.D0*rk)) &
                * den(14)
        ckplm(k,1563) = -3.35861206054688D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-277.D0+2.D0*rk) &
                * den(15)
        ckplm(k,1564) = 104.116973876953D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1565) = 5.20584869384766D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1566) = -2.6868896484375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+31.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,1567) = -0.69488525390625D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-957.D0+rk*(53.D0+5.D0*rk)) * den(2)
        ckplm(k,1568) = -0.1544189453125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(27990.D0+rk*(-5672.D0+rk*(-150.D0+17.D0*rk))) * den(3)
        ckplm(k,1569) = -0.054046630859375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-509250.D0+rk*(178192.D0+rk*(-8417.D0+rk*(-1198.D0+23.D0&
                *rk)))) * den(4)
        ckplm(k,1570) = 0.0056396484375D0*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(8.05644D6+rk*(-1.865601D6+rk*(-13919.D0+rk*(13499.D0+31.D0&
                *rk)))) * den(5)
        ckplm(k,1571) = 0.01624755859375D0*(rk+10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(3.884832D6+rk*(-565965.D0+rk*(-38630.D0+rk*(4080.D0+83.D0&
                *rk)))) * den(6)
        ckplm(k,1572) = 2.1121826171875D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(36348.D0+rk*(-2776.D0+rk*(-496.D0+(rk+19.D0)*rk))) * den(7)
        ckplm(k,1573) = 2.37620544433594D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(35100.D0+(rk*rk+rk-522.D0)*(rk+1.D0)*rk) * den(8)
        ckplm(k,1574) = 2.1121826171875D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(38610.D0+rk*(1731.D0+rk*(-547.D0+(rk-15.D0)*rk))) * den(9)
        ckplm(k,1575) = 0.01624755859375D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(4.40817D6+rk*(476797.D0+rk*(-50372.D0+rk*(-3748.D0+83.D0&
                *rk)))) * den(10)
        ckplm(k,1576) = 0.0056396484375D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(9.894654D6+rk*(1.79739D6+rk*(-54230.D0+rk*(-13375.D0+31.D0&
                *rk)))) * den(11)
        ckplm(k,1577) = -0.054046630859375D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-694638.D0+rk*(-191340.D0+rk*(-4685.D0+rk*(1290.D0+23.D0&
                *rk)))) * den(12)
        ckplm(k,1578) = -0.1544189453125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-33495.D0+rk*(-5321.D0+rk*(201.D0+17.D0*rk))) &
                * den(13)
        ckplm(k,1579) = -0.69488525390625D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1005.D0+rk*(-43.D0+5.D0*rk)) * den(14)
        ckplm(k,1580) = -2.6868896484375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-30.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1581) = 5.20584869384766D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1582) = -0.247897556849888D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1583) = 0.00799669538225446D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(775.D0+34.D0*rk) * den(1)
        ckplm(k,1584) = 0.00413622174944196D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-17835.D0+rk*(-191.D0+46.D0*rk)) * den(2)
        ckplm(k,1585) = 0.0160853068033854D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(40050.D0+rk*(-3173.D0+rk*(-339.D0+2.D0*rk))) * den(3)
        ckplm(k,1586) = -0.00064341227213542D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(7.82565D6+rk*(-1.457951D6+rk*(-51344.D0+rk*(8699.D0+146.D0*rk)))) * den(4)
        ckplm(k,1587) = -0.00001678466796875D0*(rk+10.D0)*(rk+11.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-2.2798566D9+rk*(6.60671694D8+rk*(-1.4659925D7+rk*(-6.18104D6+rk*(135065.D0+9206.D0&
                *rk))))) * den(5)
        ckplm(k,1588) = -0.00067698160807292D0*(rk+10.D0)*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-8.61966D7+rk*(1.5740982D7+rk*(779975.D0+rk*(-164120.D0+rk*(-2255.D0+218.D0&
                *rk))))) * den(6)
        ckplm(k,1589) = -0.0440038045247396D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-1.70508D6+rk*(163358.D0+rk*(27585.D0+rk*(-1720.D0+rk*(-105.D0+2.D0*rk))))) * den(7)
        ckplm(k,1590) = 5.94051361083984D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(14040.D0+(rk*rk+rk-262.D0)*(rk+1.D0)*rk) * den(8)
        ckplm(k,1591) = 0.0440038045247396D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(1.83924D6+rk*(103458.D0+rk*(-32095.D0+rk*(-1280.D0+rk*(115.D0+2.D0*rk))))) * den(9)
        ckplm(k,1592) = 0.00067698160807292D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(1.0099596D8+rk*(1.3698782D7+rk*(-1.256625D6+rk*(-152920.D0+rk*(3345.D0+218.D0&
                *rk))))) * den(10)
        ckplm(k,1593) = 0.00001678466796875D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(2.94888132D9+rk*(6.70954194D8+rk*(-4.601525D6+rk*(-6.62924D6+rk*(-89035.D0+9206.D0&
                *rk))))) * den(11)
        ckplm(k,1594) = 0.00064341227213542D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(9.223704D6+rk*(1.32975D6+rk*(-76565.D0+rk*(-8115.D0+146.D0*rk)))) &
                * den(12)
        ckplm(k,1595) = -0.0160853068033854D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-42882.D0+rk*(-2489.D0+rk*(345.D0+2.D0*rk))) * den(13)
        ckplm(k,1596) = -0.00413622174944196D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-17598.D0+rk*(283.D0+46.D0*rk)) * den(14)
        ckplm(k,1597) = -0.00799669538225446D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-741.D0+34.D0*rk) * den(15)
        ckplm(k,1598) = 0.247897556849888D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1599) = 0.011268070765904D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1600) = -0.00145394461495536D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(279.D0+14.D0*rk) * den(1)
        ckplm(k,1601) = -0.00112806047712054D0*(4.D0*rk*rk-294.D0*rk-6003.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(2)
        ckplm(k,1602) = 0.00058492024739583D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-132480.D0+rk*(-134.D0+rk*(795.D0+14.D0*rk))) * den(3)
        ckplm(k,1603) = 0.00005849202473958D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(1.261575D7+rk*(-894729.D0+rk*(-140831.D0+2.D0*rk*(1083.D0+97.D0*rk)))) * den(4)
        ckplm(k,1604) = 9.1552734375D-6*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-7.0380828D8+rk*(1.03970142D8+rk*(8.077075D6+rk*(-898070.D0+rk*(-30925.D0+758.D0*rk))))) &
                * den(5)
        ckplm(k,1605) = -5.59488932291667D-6*(rk+10.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-9.63821628D9+rk*(2.117075166D9+rk*(6.4911151D7+rk*(-2.793036D7+rk*(-5795.D0+2.D0*rk&
                *(40497.D0+62.D0*rk)))))) * den(6)
        ckplm(k,1606) = -0.00072733561197917D0*(rk-6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.0079784D8+rk*(1.1628828D7+rk*(1.848136D6+rk*(-164289.D0+rk*(-9914.D0+rk*(453.D0+10.D0&
                *rk)))))) * den(7)
        ckplm(k,1607) = -0.00981903076171875D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)&
                *(-8.4942D6+(rk+1.D0)*rk*(190950.D0+(rk*rk+rk-1097.D0)*(rk+1.D0)*rk)) * den(8)
        ckplm(k,1608) = -0.00072733561197917D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(-1.104246D8+rk*(-7.48155D6+rk*(2.277139D6+rk*(120303.D0+rk*(-12029.D0+rk*(-393.D0+10.D0&
                *rk)))))) * den(9)
        ckplm(k,1609) = -5.59488932291667D-6*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.16625366D10+rk*(-1.90388919D9+rk*(1.47859381D8+rk*(2.709972D7+rk*(-408905.D0+2.D0*rk&
                *(-40125.D0+62.D0*rk)))))) * den(10)
        ckplm(k,1610) = 9.1552734375D-6*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(7.9883496D8+rk*(8.5249272D7+rk*(-1.0578155D7+rk*(-766790.D0+rk*(34715.D0+758.D0*rk))))) &
                * den(11)
        ckplm(k,1611) = 0.00005849202473958D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(1.3367676D7+rk*(607345.D0+rk*(-146165.D0+2.D0*rk*(-695.D0+97.D0*rk)))) * den(12)
        ckplm(k,1612) = 0.00058492024739583D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(131565.D0+rk*(-1682.D0+rk*(-753.D0+14.D0*rk))) * den(13)
        ckplm(k,1613) = -0.00112806047712054D0*(4.D0*rk*rk+302.D0*rk-5705.D0)*(rk-10.D0)*(rk-11.D0)&
                *(rk-12.D0)*(rk-13.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,1614) = -0.00145394461495536D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-265.D0+14.D0*rk) * den(15)
        ckplm(k,1615) = 0.011268070765904D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1616) = -0.0004899161202567D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1617) = 0.00001580374581473D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(1519.D0+82.D0*rk) * den(1)
        ckplm(k,1618) = -1.63487025669643D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(328251.D0+rk*(23939.D0+218.D0*rk)) * den(2)
        ckplm(k,1619) = -6.35782877604167D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-1.23543D6+rk*(-68591.D0+rk*(3063.D0+134.D0*rk))) * den(3)
        ckplm(k,1620) = -1.78019205729167D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(5.074545D7+rk*(683291.D0+rk*(-436996.D0+rk*(-12359.D0+214.D0*rk)))) * den(4)
        ckplm(k,1621) = 5.340576171875D-6*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)*(1.69292088D8+rk&
                *(-6.518358D6+rk*(-2.426155D6+rk*(10160.D0+rk*(6607.D0+58.D0*rk))))) * den(5)
        ckplm(k,1622) = 0.00001958211263021D0*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(-4.2396948D8+rk&
                *(3.8189106D7+rk*(7.628587D6+rk*(-414189.D0+rk*(-40109.D0+rk*(675.D0+34.D0*rk)))))) * den(6)
        ckplm(k,1623) = 0.00003636678059896D0*(rk+8.D0)*(rk+9.D0)*(1.96952184D9+7.D0*rk&
                *(-3.7997964D7+rk*(-5.656696D6+rk*(670043.D0+rk*(37535.D0+rk*(-3001.D0+rk*(-79.D0+2.D0&
                *rk))))))) * den(7)
        ckplm(k,1624) = -0.00490951538085938D0*(rk-7.D0)*(rk+8.D0)*(-1.69884D7+7.D0*(rk+1.D0)*rk&
                *(63900.D0+(rk*rk+rk-492.D0)*(rk+1.D0)*rk)) * den(8)
        ckplm(k,1625) = -0.00003636678059896D0*(rk-7.D0)*(rk-8.D0)*(-2.1915036D9+7.D0*rk&
                *(-2.48391D7+rk*(7.412832D6+rk*(491543.D0+rk*(-51285.D0+rk*(-2485.D0+rk*(93.D0+2.D0*rk))))))) &
                * den(9)
        ckplm(k,1626) = -0.00001958211263021D0*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4.5415656D8+rk&
                *(-2.1852972D7+rk*(8.62426D6+rk*(247683.D0+rk*(-42974.D0+rk*(-471.D0+34.D0*rk)))))) * den(10)
        ckplm(k,1627) = -5.340576171875D-6*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1.7338068D8+rk&
                *(-1.661706D6+rk*(2.417573D6+rk*(-15688.D0+rk*(-6317.D0+58.D0*rk))))) * den(11)
        ckplm(k,1628) = 1.78019205729167D-6*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(4.9637736D7+rk*(-1.51935D6+rk*(-398635.D0+rk*(13215.D0+214.D0*rk)))) * den(12)
        ckplm(k,1629) = 6.35782877604167D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(1.16391D6+rk*(-74315.D0+rk*(-2661.D0+134.D0*rk))) * den(13)
        ckplm(k,1630) = 1.63487025669643D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(304530.D0+rk*(-23503.D0+218.D0*rk)) * den(14)
        ckplm(k,1631) = -0.00001580374581473D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1437.D0+82.D0*rk) * den(15)
        ckplm(k,1632) = 0.0004899161202567D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1633) = 0.00002041317167736D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+9.D0) * den(0)
        ckplm(k,1634) = -0.00001053583054315D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+9.D0)*(124.D0+7.D0*rk) * den(1)
        ckplm(k,1635) = 1.63487025669643D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0)*(23084.D0+rk*(2051.D0+37.D0*rk)) * den(2)
        ckplm(k,1636) = 1.41285083912037D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(-489240.D0+rk*(-45950.D0+rk*(-411.D0+29.D0*rk))) * den(3)
        ckplm(k,1637) = -1.41285083912037D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)&
                *(-6.731382D8+rk*(-5.1704126D7+rk*(2.388331D6+rk*(207674.D0+1721.D0*rk)))) * den(4)
        ckplm(k,1638) = -5.08626302083333D-7*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(2.13835104D8+rk&
                *(1.0179528D7+rk*(-2.13739D6+rk*(-115165.D0+rk*(2026.D0+97.D0*rk))))) * den(5)
        ckplm(k,1639) = -9.32481553819445D-7*(rk+10.D0)*(rk+9.D0)*(-1.17455184D9+rk*(-1.8442392D7+rk&
                *(1.9643018D7+rk*(554361.D0+rk*(-71569.D0+rk*(-2097.D0+23.D0*rk)))))) * den(6)
        ckplm(k,1640) = 0.00002424452039931D0*(rk+9.D0)*(-4.1225184D8+rk*(4.944048D6+rk&
                *(9.559612D6+rk*(22504.D0+rk*(-63635.D0+rk*(-713.D0+(rk+103.D0)*rk)))))) * den(7)
        ckplm(k,1641) = 0.0000454584757487D0*(1.8347472D9+(rk+1.D0)*rk*(-5.54292D7+(rk+1.D0)*rk&
                *(536652.D0+(rk*rk+rk-1748.D0)*(rk+1.D0)*rk))) * den(8)
        ckplm(k,1642) = 0.00002424452039931D0*(rk-8.D0)*(4.077216D8+rk*(-1.38573D7+rk*(-9.118944D6+rk&
                *(267889.D0+rk*(58560.D0+rk*(-1310.D0+(rk-96.D0)*rk)))))) * den(9)
        ckplm(k,1643) = -9.32481553819445D-7*(rk-8.D0)*(rk-9.D0)*(-1.13709024D9+rk*(5.5789692D7+rk&
                *(1.7571836D7+rk*(-819207.D0+rk*(-60739.D0+rk*(2235.D0+23.D0*rk)))))) * den(10)
        ckplm(k,1644) = 5.08626302083333D-7*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(2.0163528D8+rk&
                *(-1.4101194D7+rk*(-1.780709D6+rk*(122299.D0+(1541.D0-97.D0*rk)*rk)))) * den(11)
        ckplm(k,1645) = 1.41285083912037D-8*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(6.19251696D8-rk*(5.586465D7+rk*(1.775635D6+rk*(-200790.D0+1721.D0*rk)))) * den(12)
        ckplm(k,1646) = 1.41285083912037D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(443730.D0+rk*(-45041.D0+rk*(498.D0+29.D0*rk))) * den(13)
        ckplm(k,1647) = 1.63487025669643D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(21070.D0+rk*(-1977.D0+37.D0*rk)) * den(14)
        ckplm(k,1648) = -0.00001053583054315D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-117.D0+7.D0*rk) * den(15)
        ckplm(k,1649) = 0.00002041317167736D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1650) = -8.16526867094494D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0) * den(0)
        ckplm(k,1651) = 2.63395763578869D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(2511.D0+146.D0*rk) * den(1)
        ckplm(k,1652) = -4.08717564174107D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(58725.D0+rk*(5869.D0+134.D0*rk)) * den(2)
        ckplm(k,1653) = 3.53212709780093D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(1.527741D7+rk*(1.861229D6+rk*(55995.D0+46.D0*rk))) * den(3)
        ckplm(k,1654) = 3.53212709780093D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(-2.4904341D8+rk&
                *(-3.1218351D7+rk*(-518236.D0+rk*(45867.D0+1090.D0*rk)))) * den(4)
        ckplm(k,1655) = 1.9073486328125D-8*(rk+10.D0)*(rk+11.D0)*(6.005286D8+rk*(7.0104666D7+rk&
                *(-1.443715D6+rk*(-359120.D0+rk*(-7385.D0+74.D0*rk))))) * den(5)
        ckplm(k,1656) = -2.33120388454861D-8*(rk+10.D0)*(5.42178252D9+rk*(5.62173594D8+rk&
                *(-4.6399481D7+rk*(-5.855265D6+rk*(-13265.D0+rk*(8751.D0+106.D0*rk)))))) * den(6)
        ckplm(k,1657) = -1.5152825249566D-6*(-8.0525016D8+rk*(-7.4967948D7+rk*(1.236492D7+rk&
                *(1.225403D6+rk*(-37425.D0+rk*(-4537.D0+rk*(-15.D0+2.D0*rk))))))) * den(7)
        ckplm(k,1658) = 0.00020456314086914D0*(-5.6628D6+(rk+1.D0)*rk*(114588.D0+(rk*rk+rk-668.D0)&
                *(rk+1.D0)*rk)) * den(8)
        ckplm(k,1659) = 1.5152825249566D-6*(7.191756D8+rk*(-9.589446D7+rk*(-8.509264D6+rk&
                *(1.330103D6+rk*(15035.D0+rk*(-4405.D0+rk*(29.D0+2.D0*rk))))))) * den(9)
        ckplm(k,1660) = 2.33120388454861D-8*(rk-9.D0)*(4.8190428D9+rk*(-6.3750294D8+rk&
                *(-2.8999196D7+rk*(5.716815D6+rk*(-55430.D0+rk*(-8115.D0+106.D0*rk)))))) * den(10)
        ckplm(k,1661) = -1.9073486328125D-8*(rk-10.D0)*(rk-9.D0)*(-5.2933188D8+rk*(7.1944646D7+rk&
                *(411405.D0+rk*(-328840.D0+rk*(7755.D0+74.D0*rk))))) * den(11)
        ckplm(k,1662) = 3.53212709780093D-9*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(2.18388072D8+rk&
                *(-3.0048638D7+rk*(649297.D0+(41507.D0-1090.D0*rk)*rk))) * den(12)
        ckplm(k,1663) = 3.53212709780093D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)*(1.347213D7+rk&
                *(-1.749377D6+(55857.D0-46.D0*rk)*rk)) * den(13)
        ckplm(k,1664) = 4.08717564174107D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(52990.D0+rk*(-5601.D0+134.D0*rk)) * den(14)
        ckplm(k,1665) = -2.63395763578869D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0)*(-2365.D0+146.D0*rk) * den(15)
        ckplm(k,1666) = 8.16526867094494D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-9.D0) * den(16)
        ckplm(k,1667) = 3.14048795036344D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0) * den(0)
        ckplm(k,1668) = -4.05224251659799D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(775.D0+46.D0*rk) * den(1)
        ckplm(k,1669) = 3.14398126287775D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(44515.D0+2.D0*rk*(2403.D0+62.D0*rk)) * den(2)
        ckplm(k,1670) = -3.88145834923179D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(9.7848D6+rk&
                *(1.39591D6+rk*(58629.D0+626.D0*rk))) * den(3)
        ckplm(k,1671) = -3.53212709780093D-9*(rk+11.D0)*(rk+12.D0)*(-2.057181D7+rk&
                *(-3.357493D6+(58.D0*rk*rk-578.D0*rk-153907.D0)*rk)) * den(4)
        ckplm(k,1672) = 1.27156575520833D-8*(rk+11.D0)*(-8.439228D7+rk*(-1.4602314D7+rk&
                *(-541025.D0+rk*(23890.D0+rk*(1535.D0+14.D0*rk))))) * den(5)
        ckplm(k,1673) = 3.33029126364087D-9*(3.8923038D9+rk*(6.9248973D8+rk*(1.2810257D7+rk&
                *(-3.453D6+(68.D0*rk*rk+30.D0*rk-173245.D0)*rk)))) * den(6)
        ckplm(k,1674) = -3.33029126364087D-8*(3.9752856D8+rk*(3.3471852D7+rk*(-3.659752D6+rk&
                *(-319245.D0+rk*(4670.D0+rk*(513.D0+2.D0*rk)))))) * den(7)
        ckplm(k,1675) = -2.49771844773066D-7*(-5.153148D7+(rk+1.D0)*rk*(673542.D0+(rk*rk+rk-2153.D0)&
                *(rk+1.D0)*rk)) * den(8)
        ckplm(k,1676) = -3.33029126364087D-8*(3.6072036D8+rk*(-3.9817494D7+rk*(-2.679097D6+rk&
                *(332835.D0+rk*(2135.D0+rk*(-501.D0+2.D0*rk)))))) * den(9)
        ckplm(k,1677) = 3.33029126364087D-9*(3.21590412D9+rk*(-6.57202938D8+rk*(2.2130507D7+rk&
                *(2.76108D6+(68.D0*rk*rk+378.D0*rk-172375.D0)*rk)))) * den(10)
        ckplm(k,1678) = 1.27156575520833D-8*(rk-10.D0)*(7.035336D7+rk*(-1.3454664D7+rk*(603625.D0+rk&
                *(17890.D0+rk*(-1465.D0+14.D0*rk))))) * den(11)
        ckplm(k,1679) = 3.53212709780093D-9*(rk-10.D0)*(rk-11.D0)*(1.7367588D7-rk*(3.051645D6+(58.D0&
                *rk*rk+810.D0*rk-151825.D0)*rk)) * den(12)
        ckplm(k,1680) = 3.88145834923179D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(8.446893D6+rk&
                *(-1.28053D6+(56751.D0-626.D0*rk)*rk)) * den(13)
        ckplm(k,1681) = 3.14398126287775D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(39833.D0+2.D0*rk*(-2279.D0+62.D0*rk)) * den(14)
        ckplm(k,1682) = -4.05224251659799D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(-729.D0+46.D0*rk) * den(15)
        ckplm(k,1683) = 3.14048795036344D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0) * den(16)
        ckplm(k,1684) = -1.16314368531979D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0) &
                * den(0)
        ckplm(k,1685) = 3.75207640425739D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(3751.D0+226.D0*rk) * den(1)
        ckplm(k,1686) = -1.94072917461589D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(389499.D0+rk&
                *(44383.D0+1234.D0*rk)) * den(2)
        ckplm(k,1687) = 1.94072917461589D-11*(rk+12.D0)*(rk+13.D0)*(1.257311D7+rk*(1.993801D6+rk&
                *(99215.D0+1494.D0*rk))) * den(3)
        ckplm(k,1688) = -1.76606354890046D-9*(rk+12.D0)*(3.07461D6+rk*(593977.D0+rk*(38016.D0+rk&
                *(835.D0+2.D0*rk)))) * den(4)
        ckplm(k,1689) = 3.53212709780093D-10*(2.5643772D8+rk*(5.6259882D7+rk*(3.948305D6+rk&
                *(66800.D0-rk*(2405.D0+62.D0*rk))))) * den(5)
        ckplm(k,1690) = 2.0351779944472D-9*(-5.34366D7+rk*(-8.031762D6+rk*(-174365.D0+rk*(18760.D0+rk&
                *(725.D0+2.D0*rk))))) * den(6)
        ckplm(k,1691) = 1.0175889972236D-8*(1.1547432D7+rk*(845118.D0+rk*(-59615.D0+rk*(-4280.D0+rk&
                *(23.D0+2.D0*rk))))) * den(7)
        ckplm(k,1692) = -1.37374514625186D-6*(85176.D0+(rk*rk+rk-678.D0)*(rk+1.D0)*rk) * den(8)
        ckplm(k,1693) = -1.0175889972236D-8*(-1.0647D7+rk*(951426.D0+rk*(46657.D0+rk*(-4352.D0+rk&
                *(-13.D0+2.D0*rk))))) * den(9)
        ckplm(k,1694) = -2.0351779944472D-9*(4.559724D7+rk*(-7.629642D6+rk*(226315.D0+rk*(15880.D0+rk&
                *(-715.D0+2.D0*rk))))) * den(10)
        ckplm(k,1695) = 3.53212709780093D-10*(2.04057D8+rk*(-4.8572982D7+rk*(3.734095D6+rk&
                *(-75800.D0+rk*(-2095.D0+62.D0*rk))))) * den(11)
        ckplm(k,1696) = 1.76606354890046D-9*(rk-11.D0)*(2.517816D6+rk*(-520442.D0+rk*(35523.D0+rk&
                *(-827.D0+2.D0*rk)))) * den(12)
        ckplm(k,1697) = 1.94072917461589D-11*(rk-11.D0)*(rk-12.D0)*(1.067703D7+rk&
                *(-1.799853D6+(94733.D0-1494.D0*rk)*rk)) * den(13)
        ckplm(k,1698) = 1.94072917461589D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(346350.D0+rk&
                *(-41915.D0+1234.D0*rk)) * den(14)
        ckplm(k,1699) = -3.75207640425739D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(-3525.D0+226.D0*rk) * den(15)
        ckplm(k,1700) = 1.16314368531979D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0) &
                * den(16)
        ckplm(k,1701) = 4.15408459042783D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0) * den(0)
        ckplm(k,1702) = -2.14404365957565D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(279.D0+17.D0*rk) &
                * den(1)
        ckplm(k,1703) = 3.32696429934153D-12*(rk+13.D0)*(rk+14.D0)*(114579.D0+rk*(13581.D0+397.D0&
                *rk)) * den(2)
        ckplm(k,1704) = -7.76291669846357D-11*(rk+13.D0)*(186170.D0+rk*(31832.D0+rk*(1758.D0+31.D0&
                *rk))) * den(3)
        ckplm(k,1705) = 5.04589585400132D-11*(7.35105D6+rk*(1.597088D6+rk*(122477.D0+rk&
                *(3798.D0+37.D0*rk)))) * den(4)
        ckplm(k,1706) = 6.05507502480159D-10*(-957768.D0+rk*(-167239.D0+rk*(-9073.D0+(rk-131.D0)&
                *rk))) * den(5)
        ckplm(k,1707) = -3.70032362626764D-10*(-2.05464D6+rk*(-249491.D0+rk*(-5090.D0+rk*(224.D0+5.D0&
                *rk)))) * den(6)
        ckplm(k,1708) = -1.05723532179075D-10*(8.219484D6+rk*(497880.D0+rk*(-21400.D0+(rk-1125.D0)&
                *rk))) * den(7)
        ckplm(k,1709) = 1.7840846055219D-9*(496860.D0+(rk*rk+rk-2186.D0)*(rk+1.D0)*rk) * den(8)
        ckplm(k,1710) = -1.05723532179075D-10*(7.70133D6+rk*(-537301.D0+rk*(-18019.D0+(rk+1129.D0)&
                *rk))) * den(9)
        ckplm(k,1711) = -3.70032362626764D-10*(-1.810458D6+rk*(238659.D0+rk*(-5732.D0+rk&
                *(-204.D0+5.D0*rk)))) * den(10)
        ckplm(k,1712) = 6.05507502480159D-10*(-799470.D0+rk*(149490.D0+rk*(-8674.D0+(rk+135.D0)&
                *rk))) * den(11)
        ckplm(k,1713) = 5.04589585400132D-11*(5.872678D6+rk*(-1.36338D6+rk*(111305.D0+rk&
                *(-3650.D0+37.D0*rk)))) * den(12)
        ckplm(k,1714) = 7.76291669846357D-11*(rk-12.D0)*(156065.D0+rk*(-28409.D0+(1665.D0-31.D0*rk)&
                *rk)) * den(13)
        ckplm(k,1715) = 3.32696429934153D-12*(rk-12.D0)*(rk-13.D0)*(101395.D0+rk*(-12787.D0+397.D0&
                *rk)) * den(14)
        ckplm(k,1716) = -2.14404365957565D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(-262.D0+17.D0*rk) &
                * den(15)
        ckplm(k,1717) = 4.15408459042783D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0) * den(16)
        ckplm(k,1718) = -1.43244296221649D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0) * den(0)
        ckplm(k,1719) = 4.62078374908546D-14*(rk+14.D0)*(rk+15.D0)*(5239.D0+322.D0*rk) * den(1)
        ckplm(k,1720) = -1.38623512472564D-13*(rk+14.D0)*(130299.D0+rk*(15911.D0+482.D0*rk)) * den(2)
        ckplm(k,1721) = 4.85182293653973D-12*(163930.D0+rk*(29643.D0+rk*(1757.D0+34.D0*rk))) * den(3)
        ckplm(k,1722) = 1.26147396350033D-11*(-142450.D0-rk*(22779.D0+rk*(1153.D0+18.D0*rk))) * den(4)
        ckplm(k,1723) = 1.26147396350033D-11*(249438.D0+rk*(32435.D0+rk*(1197.D0+10.D0*rk))) * den(5)
        ckplm(k,1724) = 4.62540453283455D-11*(-97230.D0+rk*(-8915.D0+rk*(-141.D0+2.D0*rk))) * den(6)
        ckplm(k,1725) = -8.59003698954987D-11*(-63294.D0+rk*(-2957.D0+rk*(69.D0+2.D0*rk))) * den(7)
        ckplm(k,1726) = 1.15965499358923D-8*(rk*rk+rk-490.D0) * den(8)
        ckplm(k,1727) = 8.59003698954987D-11*(60270.D0+rk*(-3089.D0+rk*(-63.D0+2.D0*rk))) * den(9)
        ckplm(k,1728) = -4.62540453283455D-11*(88458.D0+rk*(-8627.D0+rk*(147.D0+2.D0*rk))) * den(10)
        ckplm(k,1729) = 1.26147396350033D-11*(218190.D0+rk*(-30071.D0+(1167.D0-10.D0*rk)*rk)) &
                * den(11)
        ckplm(k,1730) = 1.26147396350033D-11*(-120806.D0+rk*(20527.D0+rk*(-1099.D0+18.D0*rk))) &
                * den(12)
        ckplm(k,1731) = 4.85182293653973D-12*(136010.D0+rk*(-26231.D0+(1655.D0-34.D0*rk)*rk)) &
                * den(13)
        ckplm(k,1732) = 1.38623512472564D-13*(rk-13.D0)*(114870.D0+rk*(-14947.D0+482.D0*rk)) * den(14)
        ckplm(k,1733) = 4.62078374908546D-14*(rk-13.D0)*(rk-14.D0)*(4917.D0-322.D0*rk) * den(15)
        ckplm(k,1734) = 1.43244296221649D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0) * den(16)
        ckplm(k,1735) = 4.77480987405498D-14*(rk+15.D0)*(rk+16.D0) * den(0)
        ckplm(k,1736) = -6.16104499878061D-15*(rk+15.D0)*(1519.D0+94.D0*rk) * den(1)
        ckplm(k,1737) = 4.62078374908546D-14*(68.D0*rk*rk+2186.D0*rk+17493.D0) * den(2)
        ckplm(k,1738) = 2.15636574957322D-13*(-13560.D0-rk*(1591.D0+46.D0*rk)) * den(3)
        ckplm(k,1739) = 9.81146416055813D-12*(765.D0+2.D0*(rk+40.D0)*rk) * den(4)
        ckplm(k,1740) = 1.17737569926698D-11*(-1245.D0-rk*(107.D0+2.D0*rk)) * den(5)
        ckplm(k,1741) = 3.0836030218897D-12*(4.D0*rk*rk+454.D0*rk+7395.D0) * den(6)
        ckplm(k,1742) = 4.40514717412814D-12*(-6618.D0+rk*(-211.D0+2.D0*rk)) * den(7)
        ckplm(k,1743) = -1.98231622835766D-11*(rk*rk+rk-1575.D0) * den(8)
        ckplm(k,1744) = 4.40514717412814D-12*(-6405.D0+rk*(215.D0+2.D0*rk)) * den(9)
        ckplm(k,1745) = 3.0836030218897D-12*(4.D0*rk*rk-446.D0*rk+6945.D0) * den(10)
        ckplm(k,1746) = 1.17737569926698D-11*(-1140.D0+(103.D0-2.D0*rk)*rk) * den(11)
        ckplm(k,1747) = 9.81146416055813D-12*(687.D0+2.D0*(rk-38.D0)*rk) * den(12)
        ckplm(k,1748) = 2.15636574957322D-13*(-12015.D0+(1499.D0-46.D0*rk)*rk) * den(13)
        ckplm(k,1749) = 4.62078374908546D-14*(68.D0*rk*rk-2050.D0*rk+15375.D0) * den(14)
        ckplm(k,1750) = 6.16104499878061D-15*(rk-14.D0)*(1425.D0-94.D0*rk) * den(15)
        ckplm(k,1751) = 4.77480987405498D-14*(rk-14.D0)*(rk-15.D0) * den(16)
        ckplm(k,1752) = 1.54026124969515D-15*(-16.D0-rk) * den(0)
        ckplm(k,1753) = 1.54026124969515D-15*(225.D0+14.D0*rk) * den(1)
        ckplm(k,1754) = 6.93117562362819D-14*(-33.D0-2.D0*rk) * den(2)
        ckplm(k,1755) = 2.69545718696652D-13*(35.D0+2.D0*rk) * den(3)
        ckplm(k,1756) = 7.00818868611295D-13*(-39.D0-2.D0*rk) * den(4)
        ckplm(k,1757) = 1.26147396350033D-12*(47.D0+2.D0*rk) * den(5)
        ckplm(k,1758) = -1.54180151094485D-12*(65.D0+2.D0*rk) * den(6)
        ckplm(k,1759) = 1.10128679353204D-12*(123.D0+2.D0*rk) * den(7)
        ckplm(k,1760) = -1.48673717126825D-10 * den(8)
        ckplm(k,1761) = -1.10128679353204D-12*(-121.D0+2.D0*rk) * den(9)
        ckplm(k,1762) = 1.54180151094485D-12*(-63.D0+2.D0*rk) * den(10)
        ckplm(k,1763) = 1.26147396350033D-12*(45.D0-2.D0*rk) * den(11)
        ckplm(k,1764) = 7.00818868611295D-13*(-37.D0+2.D0*rk) * den(12)
        ckplm(k,1765) = 2.69545718696652D-13*(33.D0-2.D0*rk) * den(13)
        ckplm(k,1766) = 6.93117562362819D-14*(-31.D0+2.D0*rk) * den(14)
        ckplm(k,1767) = 1.54026124969515D-15*(211.D0-14.D0*rk) * den(15)
        ckplm(k,1768) = 1.54026124969515D-15*(rk-15.D0) * den(16)
        ckplm(k,1769) = 4.81331640529735D-17 * den(0)
        ckplm(k,1770) = -7.70130624847577D-16 * den(1)
        ckplm(k,1771) = 5.77597968635683D-15 * den(2)
        ckplm(k,1772) = -2.69545718696652D-14 * den(3)
        ckplm(k,1773) = 8.76023585764119D-14 * den(4)
        ckplm(k,1774) = -2.10245660583388D-13 * den(5)
        ckplm(k,1775) = 3.85450377736212D-13 * den(6)
        ckplm(k,1776) = -5.50643396766017D-13 * den(7)
        ckplm(k,1777) = 6.1947382136177D-13 * den(8)
        ckplm(k,1778) = -5.50643396766017D-13 * den(9)
        ckplm(k,1779) = 3.85450377736212D-13 * den(10)
        ckplm(k,1780) = -2.10245660583388D-13 * den(11)
        ckplm(k,1781) = 8.76023585764119D-14 * den(12)
        ckplm(k,1782) = -2.69545718696652D-14 * den(13)
        ckplm(k,1783) = 5.77597968635683D-15 * den(14)
        ckplm(k,1784) = -7.70130624847577D-16 * den(15)
        ckplm(k,1785) = 4.81331640529735D-17 * den(16)
!    ckplm para l = 17
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)&
                *(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(17) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,1786) = 623140.088653565D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1787) = 321011.560821533D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,1788) = 248525.079345703D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,1789) = 214245.758056641D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,1790) = 194408.187866211D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,1791) = 181966.063842773D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,1792) = 174054.495849609D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(6)
        ckplm(k,1793) = 169318.319091797D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(7)
        ckplm(k,1794) = 167090.446472168D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(8)
        ckplm(k,1795) = 167090.446472168D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*rk * den(9)
        ckplm(k,1796) = 169318.319091797D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,1797) = 174054.495849609D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,1798) = 181966.063842773D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,1799) = 194408.187866211D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,1800) = 214245.758056641D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,1801) = 248525.079345703D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(15)
        ckplm(k,1802) = 321011.560821533D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(16)
        ckplm(k,1803) = 623140.088653565D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*rk * den(17)
        ckplm(k,1804) = -69237.7876281738D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1805) = -6294.34432983399D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-11.D0+5.D0*rk) * den(1)
        ckplm(k,1806) = -1624.34692382813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-62.D0+13.D0*rk) * den(2)
        ckplm(k,1807) = -1400.29907226563D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-87.D0+11.D0*rk) * den(3)
        ckplm(k,1808) = -11435.7757568359D0*(rk+10.D0)*(rk+11.D0)*(rk-12.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(4)
        ckplm(k,1809) = -1189.32067871094D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-125.D0+7.D0*rk) * den(5)
        ckplm(k,1810) = -1137.61108398438D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-138.D0+5.D0*rk) * den(6)
        ckplm(k,1811) = -3319.96704101563D0*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-49.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(7)
        ckplm(k,1812) = -1092.09442138672D0*(rk-152.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(8)
        ckplm(k,1813) = 1092.09442138672D0*(rk+153.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0) * den(9)
        ckplm(k,1814) = 3319.96704101563D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+50.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(10)
        ckplm(k,1815) = 1137.61108398438D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(143.D0+5.D0*rk) * den(11)
        ckplm(k,1816) = 1189.32067871094D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(132.D0+7.D0*rk) * den(12)
        ckplm(k,1817) = 11435.7757568359D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(13)
        ckplm(k,1818) = 1400.29907226563D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(98.D0+11.D0*rk) * den(14)
        ckplm(k,1819) = 1624.34692382813D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(75.D0+13.D0*rk) * den(15)
        ckplm(k,1820) = 6294.34432983399D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(16.D0+5.D0*rk) * den(16)
        ckplm(k,1821) = 69237.7876281738D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,1822) = 3644.09408569336D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,1823) = 331.281280517578D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-44.D0+3.D0*rk) * den(1)
        ckplm(k,1824) = 21.3729858398438D0*(8.D0*rk*rk-876.D0*rk+1891.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(2)
        ckplm(k,1825) = -18.4249877929688D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-3741.D0+4.D0&
                *rk*(257.D0+4.D0*rk)) * den(3)
        ckplm(k,1826) = -300.941467285156D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-321.D0+2.D0&
                *(rk+29.D0)*rk) * den(4)
        ckplm(k,1827) = -31.2979125976563D0*(26.D0*rk*rk+474.D0*rk-3875.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk+12.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(5)
        ckplm(k,1828) = -14.9685668945313D0*(64.D0*rk*rk+764.D0*rk-9453.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(6)
        ckplm(k,1829) = -43.6837768554688D0*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-3577.D0+4.D0&
                *rk*(43.D0+6.D0*rk)) * den(7)
        ckplm(k,1830) = -1092.09442138672D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-151.D0+(rk+3.D0)*rk) * den(8)
        ckplm(k,1831) = -1092.09442138672D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(-153.D0+(rk-1.D0)*rk) * den(9)
        ckplm(k,1832) = -43.6837768554688D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-3725.D0+4.D0&
                *rk*(-31.D0+6.D0*rk)) * den(10)
        ckplm(k,1833) = -14.9685668945313D0*(64.D0*rk*rk-636.D0*rk-10153.D0)*(rk-10.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(11)
        ckplm(k,1834) = -31.2979125976563D0*(26.D0*rk*rk-422.D0*rk-4323.D0)*(rk-10.D0)*(rk-11.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,1835) = -300.941467285156D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-377.D0+2.D0&
                *(rk-27.D0)*rk) * den(13)
        ckplm(k,1836) = -18.4249877929688D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4753.D0+4.D0&
                *rk*(-249.D0+4.D0*rk)) * den(14)
        ckplm(k,1837) = 21.3729858398438D0*(8.D0*rk*rk+892.D0*rk+2775.D0)*(rk-10.D0)*(rk-11.D0)&
                *(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(15)
        ckplm(k,1838) = 331.281280517578D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(47.D0+3.D0*rk) * den(16)
        ckplm(k,1839) = 3644.09408569336D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(17)
        ckplm(k,1840) = -182.204704284668D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,1841) = 5.52135467529297D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+297.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(1)
        ckplm(k,1842) = 2.13729858398438D0*(26.D0*rk*rk+718.D0*rk-3813.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(2)
        ckplm(k,1843) = 1.10549926757813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(35235.D0+(66.D0*rk*rk+788.D0&
                *rk-14513.D0)*rk) * den(3)
        ckplm(k,1844) = 2.00627644856771D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(33840.D0+rk*(-9187.D0+rk&
                *(105.D0+37.D0*rk))) * den(4)
        ckplm(k,1845) = 9.38937377929688D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(10525.D0+rk*(-1936.D0+rk&
                *(-46.D0+7.D0*rk))) * den(5)
        ckplm(k,1846) = 1.49685668945313D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(85146.D0+(34.D0*rk*rk-654.D0&
                *rk-10351.D0)*rk) * den(6)
        ckplm(k,1847) = 1.45612589518229D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(103047.D0+(22.D0*rk*rk-948.D0&
                *rk-7453.D0)*rk) * den(7)
        ckplm(k,1848) = 10.9209442138672D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(15000.D0+rk*(-448.D0+(rk-147.D0)&
                *rk)) * den(8)
        ckplm(k,1849) = -10.9209442138672D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(-15300.D0+(rk+151.D0)*(rk-1.D0)&
                *rk) * den(9)
        ckplm(k,1850) = -1.45612589518229D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-109530.D0+rk*(-5491.D0+2.D0*rk&
                *(507.D0+11.D0*rk))) * den(10)
        ckplm(k,1851) = -1.49685668945313D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-94809.D0+(34.D0*rk*rk+756.D0&
                *rk-8941.D0)*rk) * den(11)
        ckplm(k,1852) = -9.38937377929688D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-12408.D0+rk*(-1823.D0+rk&
                *(67.D0+7.D0*rk))) * den(12)
        ckplm(k,1853) = -2.00627644856771D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-43095.D0+rk*(-9286.D0+rk&
                *(6.D0+37.D0*rk))) * den(13)
        ckplm(k,1854) = -1.10549926757813D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-50470.D0+(66.D0*rk*rk-590.D0&
                *rk-15891.D0)*rk) * den(14)
        ckplm(k,1855) = -2.13729858398438D0*(26.D0*rk*rk-666.D0*rk-4505.D0)*(rk-10.D0)*(rk-11.D0)&
                *(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(15)
        ckplm(k,1856) = -5.52135467529297D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-296.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(16)
        ckplm(k,1857) = 182.204704284668D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(17)
        ckplm(k,1858) = 8.67641448974609D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1859) = -0.788764953613281D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(176.D0+5.D0*rk) &
                * den(1)
        ckplm(k,1860) = -0.61065673828125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1829.D0+rk*(115.D0+9.D0*rk)) &
                * den(2)
        ckplm(k,1861) = -0.03509521484375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(209670.D0+rk*(-44551.D0+rk*(-729.D0+127.D0&
                *rk))) * den(3)
        ckplm(k,1862) = -0.0409444173177083D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.16487D6+rk*(421416.D0+rk*(-24151.D0+rk&
                *(-2406.D0+61.D0*rk)))) * den(4)
        ckplm(k,1863) = -0.0212910970052083D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-3.77895D6+rk*(928348.D0+rk*(-6863.D0+rk&
                *(-5902.D0+17.D0*rk)))) * den(5)
        ckplm(k,1864) = 0.06109619140625D0*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.877904D6+rk*(-305199.D0+rk*(-12996.D0+rk&
                *(1966.D0+25.D0*rk)))) * den(6)
        ckplm(k,1865) = 2.91225179036458D0*(rk+10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(49464.D0+rk*(-4785.D0+rk*(-550.D0+(rk+30.D0)*rk))) &
                * den(7)
        ckplm(k,1866) = 3.64031473795573D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(44700.D0+rk*(-1786.D0+rk*(-577.D0+(rk+10.D0)*rk))) &
                * den(8)
        ckplm(k,1867) = 3.64031473795573D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(45900.D0+(rk-1.D0)*rk*(-606.D0+(rk-5.D0)*rk)) * den(9)
        ckplm(k,1868) = 2.91225179036458D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(53670.D0+rk*(3599.D0+rk*(-634.D0+(rk-26.D0)*rk))) &
                * den(10)
        ckplm(k,1869) = 0.06109619140625D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.168166D6+rk*(273409.D0+rk*(-18744.D0+rk&
                *(-1866.D0+25.D0*rk)))) * den(11)
        ckplm(k,1870) = -0.0212910970052083D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4.708242D6+rk*(-924300.D0+rk*(10945.D0+rk&
                *(5970.D0+17.D0*rk)))) * den(12)
        ckplm(k,1871) = -0.0409444173177083D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1.60797D6+rk*(-462256.D0+rk*(-16567.D0+rk&
                *(2650.D0+61.D0*rk)))) * den(13)
        ckplm(k,1872) = -0.03509521484375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-253365.D0+rk*(-42712.D0+rk*(1110.D0+127.D0&
                *rk))) * den(14)
        ckplm(k,1873) = -0.61065673828125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1935.D0+rk*(-97.D0+9.D0*rk)) &
                * den(15)
        ckplm(k,1874) = -0.788764953613281D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-171.D0+5.D0*rk) &
                * den(16)
        ckplm(k,1875) = 8.67641448974609D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,1876) = -0.394382476806641D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1877) = 0.394382476806641D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+25.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,1878) = 0.152664184570313D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-775.D0+2.D0*(rk-1.D0)*rk) * den(2)
        ckplm(k,1879) = 0.0438690185546875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(23925.D0+rk*(-2171.D0+2.D0*(rk-90.D0)*rk)) * den(3)
        ckplm(k,1880) = -0.102361043294271D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(81900.D0+rk*(-16437.D0+rk*(-304.D0+(rk+90.D0)*rk))) * den(4)
        ckplm(k,1881) = -0.00081888834635417D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-7.99554D7+rk*(2.4561091D7+rk*(-951075.D0+rk*(-191060.D0+rk&
                *(6285.D0+259.D0*rk))))) * den(5)
        ckplm(k,1882) = -0.0001068115234375D0*(rk+10.D0)*(rk+11.D0)*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-9.66483D8+rk*(1.96753734D8+rk*(4.287575D6+(2166.D0*rk*rk-3260.D0&
                *rk-1778565.D0)*rk))) * den(6)
        ckplm(k,1883) = -0.00509134928385417D0*(rk+10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-2.71557D7+rk*(3.293391D6+rk*(337750.D0+rk*(-30985.D0+2.D0*rk&
                *(-545.D0+17.D0*rk))))) * den(7)
        ckplm(k,1884) = -0.0636418660481771D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-2.53968D6+rk*(127264.D0+rk*(40464.D0+rk*(-1139.D0+(rk-138.D0)*rk)))) &
                * den(8)
        ckplm(k,1885) = 0.0636418660481771D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(2.62548D6+(rk-1.D0)*rk*(-43476.D0+rk*(-433.D0+(rk+144.D0)*rk))) * den(9)
        ckplm(k,1886) = 0.00509134928385417D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(3.008148D7+rk*(2.529466D6+rk*(-423825.D0+rk*(-26285.D0+2.D0*rk&
                *(630.D0+17.D0*rk))))) * den(10)
        ckplm(k,1887) = 0.0001068115234375D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(1.15717602D9+rk*(1.82866759D8+rk*(-9.58205D6+rk*(-1.743865D6+2.D0*rk&
                *(7045.D0+1083.D0*rk))))) * den(11)
        ckplm(k,1888) = 0.00081888834635417D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(1.0527048D8+rk*(2.5866216D7+rk*(342775.D0+rk*(-213610.D0+rk&
                *(-4990.D0+259.D0*rk))))) * den(12)
        ckplm(k,1889) = 0.102361043294271D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(97944.D0+rk*(15563.D0+rk*(-568.D0+(rk-86.D0)*rk))) * den(13)
        ckplm(k,1890) = -0.0438690185546875D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-25914.D0+rk*(-1805.D0+2.D0*(rk+93.D0)*rk)) * den(14)
        ckplm(k,1891) = -0.152664184570313D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-771.D0+2.D0*(rk+3.D0)*rk) * den(15)
        ckplm(k,1892) = -0.394382476806641D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-24.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,1893) = 0.394382476806641D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,1894) = 0.0171470642089844D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1895) = -0.00571568806966146D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(108.D0+5.D0*rk) * den(1)
        ckplm(k,1896) = -0.00110626220703125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-9393.D0+4.D0*rk*(-95.D0+2.D0*rk)) * den(2)
        ckplm(k,1897) = 0.00057220458984375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-211410.D0+rk*(2497.D0+4.D0*rk*(307.D0+4.D0*rk))) * den(3)
        ckplm(k,1898) = 0.00014834933810764D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(7.96257D6+rk*(-697299.D0+rk*(-76891.D0+2.D0*rk*(1032.D0+53.D0*rk)))) &
                * den(4)
        ckplm(k,1899) = 0.00005340576171875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-1.995126D8+rk*(3.3570674D7+rk*(1.563875D6+(226.D0*rk*rk-5360.D0*rk-264165.D0)&
                *rk))) * den(5)
        ckplm(k,1900) = 0.00001780192057292D0*(rk+10.D0)*(rk+11.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(5.21539956D9+rk*(-1.276017642D9+rk*(-1.080785D6+rk*(1.411587D7+rk*(-227255.D0+4.D0*rk&
                *(-9027.D0+40.D0*rk)))))) * den(6)
        ckplm(k,1901) = -0.00084855821397569D0*(rk+10.D0)*(rk-6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.5634476D8+rk*(2.2816422D7+rk*(2.054327D6+rk*(-282210.D0+(8.D0*rk*rk+708.D0*rk-8455.D0)&
                *rk)))) * den(7)
        ckplm(k,1902) = -0.0127283732096354D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.26126D7+rk*(760950.D0+rk*(238009.D0+rk*(-9363.D0+rk*(-1202.D0+(rk+21.D0)*rk))))) * den(8)
        ckplm(k,1903) = -0.0127283732096354D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(-1.31274D7+(rk-1.D0)*rk*(261750.D0+rk*(3059.D0+rk*(-1306.D0+(rk-14.D0)*rk)))) * den(9)
        ckplm(k,1904) = -0.00084855821397569D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.768338D8+rk*(-1.789845D7+rk*(2.843267D6+rk*(241470.D0+(8.D0*rk*rk-660.D0*rk-11875.D0)&
                *rk)))) * den(10)
        ckplm(k,1905) = 0.00001780192057292D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(6.47602956D9+rk*(1.230780942D9+rk*(-4.4428445D7+rk*(-1.466061D7+rk*(-44315.D0+4.D0*rk&
                *(9267.D0+40.D0*rk)))))) * den(11)
        ckplm(k,1906) = 0.00005340576171875D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(2.3126082D8+rk*(2.9672999D7+rk*(-2.32195D6+(226.D0*rk*rk+6490.D0*rk-240465.D0)&
                *rk))) * den(12)
        ckplm(k,1907) = 0.00014834933810764D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(8.58102D6+rk*(537749.D0+rk*(-82447.D0+2.D0*rk*(-820.D0+53.D0*rk)))) &
                * den(13)
        ckplm(k,1908) = 0.00057220458984375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(212695.D0+rk*(89.D0+4.D0*rk*(-295.D0+4.D0*rk))) * den(14)
        ckplm(k,1909) = -0.00110626220703125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-9005.D0+4.D0*rk*(99.D0+2.D0*rk)) * den(15)
        ckplm(k,1910) = -0.00571568806966146D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-103.D0+5.D0*rk) * den(16)
        ckplm(k,1911) = 0.0171470642089844D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,1912) = -0.00071446100870768D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,1913) = 0.00006495100079161D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0)*(539.D0+27.D0*rk) * den(1)
        ckplm(k,1914) = -0.00001676154859138D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(47089.D0+rk*(3021.D0+17.D0*rk)) * den(2)
        ckplm(k,1915) = -2.88992217092803D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-4.071165D6+rk*(-169762.D0+rk*(11652.D0+379.D0*rk))) * den(3)
        ckplm(k,1916) = -0.00002360103106258D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(5.89554D6+rk*(-30408.D0+rk*(-49199.D0+rk*(-882.D0+29.D0*rk)))) * den(4)
        ckplm(k,1917) = 9.44041242503157D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.5310407D9+rk*(-9.4912048D7+rk*(-1.9300275D7+rk*(393155.D0+rk*(51495.D0+173.D0*rk))))) &
                * den(5)
        ckplm(k,1918) = 0.00002076890733507D0*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-6.6986388D8+rk*(7.9739226D7+rk*(9.421063D6+rk*(-819885.D0+rk*(-43640.D0+rk*(1479.D0+37.D0&
                *rk)))))) * den(6)
        ckplm(k,1919) = 0.00003857082790799D0*(rk+10.D0)*(rk+8.D0)*(rk+9.D0)*(3.29971356D9+7.D0*rk&
                *(-8.0464986D7+rk*(-6.227749D6+rk*(1.218342D6+rk*(26435.D0+rk*(-4839.D0+rk*(-46.D0+3.D0&
                *rk))))))) * den(7)
        ckplm(k,1920) = 0.00004821353488498D0*(rk-7.D0)*(rk+8.D0)*(rk+9.D0)*(3.3070752D9+7.D0*rk&
                *(-3.33648D7+rk*(-1.0261124D7+rk*(522532.D0+rk*(68209.D0+rk*(-1973.D0+(rk-125.D0)*rk)))))) &
                * den(8)
        ckplm(k,1921) = -0.00004821353488498D0*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(-3.4656336D9+7.D0&
                *(rk-1.D0)*rk*(1.15569D7+rk*(155268.D0+rk*(-77233.D0+rk*(-1069.D0+(rk+133.D0)*rk))))) * den(9)
        ckplm(k,1922) = -0.00003857082790799D0*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-3.8110644D9+7.D0*rk&
                *(-6.44841D7+rk*(9.676528D6+rk*(1.065237D6+rk*(-49835.D0+rk*(-4500.D0+rk*(67.D0+3.D0*rk))))))) &
                * den(10)
        ckplm(k,1923) = -0.00002076890733507D0*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-7.3940724D8+rk*(-5.8619178D7+rk*(1.1604643D7+rk*(631275.D0+rk*(-50480.D0+rk*(-1257.D0+37.D0&
                *rk)))))) * den(11)
        ckplm(k,1924) = -9.44041242503157D-7*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.60631064D9+rk*(-5.5337148D7+rk*(2.01725D7+rk*(188905.D0+rk*(-50630.D0+173.D0*rk))))) &
                * den(12)
        ckplm(k,1925) = 0.00002360103106258D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(5.87766D6+rk*(-65228.D0+rk*(-46379.D0+rk*(998.D0+29.D0*rk)))) * den(13)
        ckplm(k,1926) = 2.88992217092803D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(3.89013D6+rk*(-191929.D0+rk*(-10515.D0+379.D0*rk))) * den(14)
        ckplm(k,1927) = 0.00001676154859138D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(44085.D0+rk*(-2987.D0+17.D0*rk)) * den(15)
        ckplm(k,1928) = -0.00006495100079161D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-512.D0+27.D0*rk) * den(16)
        ckplm(k,1929) = 0.00071446100870768D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,1930) = 0.00002857844034831D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+9.D0) * den(0)
        ckplm(k,1931) = -2.5980400316643D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+9.D0)*(704.D0+37.D0*rk) * den(1)
        ckplm(k,1932) = 6.70461943655303D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+9.D0)*(79484.D0+rk*(6351.D0+97.D0*rk)) * den(2)
        ckplm(k,1933) = 1.15596886837121D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0)*(-8.60604D6+rk*(-682786.D0+rk*(-159.D0+517.D0*rk))) * den(3)
        ckplm(k,1934) = -1.34863034643308D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(-1.0460268D8+rk*(-5.909382D6+rk*(458239.D0+rk*(26562.D0+101.D0*rk)))) * den(4)
        ckplm(k,1935) = -2.69726069286616D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)&
                *(6.23450352D9+rk*(1.33273856D8+rk*(-6.342705D7+rk*(-2.138965D6+rk*(82230.D0+2129.D0*rk))))) &
                * den(5)
        ckplm(k,1936) = -5.93397352430556D-7*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(-3.01056336D9+rk&
                *(5.1479208D7+rk*(4.6949594D7+rk*(228705.D0+rk*(-173425.D0+rk*(-2553.D0+71.D0*rk)))))) * den(6)
        ckplm(k,1937) = 7.71416558159722D-6*(rk+10.D0)*(rk+9.D0)*(-2.26058976D9+rk*(1.19110608D8+rk&
                *(4.5589412D7+rk*(-1.318796D6+rk*(-279535.D0+rk*(2347.D0+(rk+443.D0)*rk)))))) * den(7)
        ckplm(k,1938) = 0.00004821353488498D0*(rk+9.D0)*(3.284424D9+rk*(-2.65973328D8+rk&
                *(-8.0389252D7+rk*(5.058756D6+rk*(658929.D0+rk*(-27288.D0+rk*(-1902.D0+(rk+36.D0)*rk))))))) &
                * den(8)
        ckplm(k,1939) = 0.00004821353488498D0*(rk-8.D0)*(3.4656336D9+(rk-1.D0)*rk*(-9.2779056D7+rk&
                *(-1.410732D6+rk*(778672.D0+rk*(13023.D0+rk*(-2153.D0+(rk-27.D0)*rk)))))) * den(9)
        ckplm(k,1940) = 7.71416558159722D-6*(rk-8.D0)*(rk-9.D0)*(2.3330736D9+rk*(2.510262D7+rk&
                *(-4.7851744D7+rk*(-186011.D0+rk*(284660.D0+rk*(-290.D0+(rk-436.D0)*rk)))))) * den(10)
        ckplm(k,1941) = -5.93397352430556D-7*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(-3.01549248D9+rk&
                *(4.1053356D7+rk*(4.5249524D7+rk*(-895455.D0+rk*(-159595.D0+rk*(2979.D0+71.D0*rk)))))) * den(11)
        ckplm(k,1942) = -2.69726069286616D-8*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-6.04002168D9+rk*(2.53392786D8+rk*(5.6538065D7+rk*(-2.446595D6+rk*(-71585.D0+2129.D0*rk))))) &
                * den(12)
        ckplm(k,1943) = -1.34863034643308D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-9.826152D7+rk*(6.746578D6+rk*(379159.D0+rk*(-26158.D0+101.D0*rk)))) * den(13)
        ckplm(k,1944) = 1.15596886837121D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(7.92393D6+rk*(-680917.D0+rk*(1710.D0+517.D0*rk))) * den(14)
        ckplm(k,1945) = 6.70461943655303D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0)*(73230.D0+rk*(-6157.D0+97.D0*rk)) * den(15)
        ckplm(k,1946) = -2.5980400316643D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-8.D0)*(rk-9.D0)*(-667.D0+37.D0*rk) * den(16)
        ckplm(k,1947) = 0.00002857844034831D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,1948) = -1.0991707826272D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0) * den(0)
        ckplm(k,1949) = 3.33082055341577D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(2673.D0+145.D0*rk) * den(1)
        ckplm(k,1950) = -2.57869978328963D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(126387.D0+rk*(11485.D0+233.D0*rk)) * den(2)
        ckplm(k,1951) = 4.44603410912005D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(1.6807095D7+rk*(1.792648D6+(42132.D0-211.D0*rk)*rk)) * den(3)
        ckplm(k,1952) = 5.76337754885932D-10*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(-2.18234898D9+rk*(-2.26084332D8+rk*(-648803.D0+rk*(410622.D0+6893.D0*rk)))) * den(4)
        ckplm(k,1953) = 1.34863034643308D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(1.27062054D9+rk&
                *(1.11805572D8+rk*(-5.607565D6+rk*(-635135.D0+rk*(-6575.D0+203.D0*rk))))) * den(5)
        ckplm(k,1954) = -2.96698676215278D-7*(rk+10.D0)*(rk+11.D0)*(6.7641588D8+rk*(4.5442686D7+rk&
                *(-6.993295D6+rk*(-515745.D0+rk*(9110.D0+rk*(879.D0+5.D0*rk)))))) * den(6)
        ckplm(k,1955) = -9.88995587384259D-8*(rk+10.D0)*(-2.113526844D10+rk*(-1.001594898D9+rk&
                *(3.50351595D8+rk*(1.8143068D7+rk*(-1.421025D6+rk*(-73607.D0+rk*(870.D0+37.D0*rk))))))) * den(7)
        ckplm(k,1956) = -1.85436672634549D-6*(1.06007616D10+rk*(3.56119776D8+rk*(-2.4315786D8+rk&
                *(-8.709544D6+rk*(1.699669D6+rk*(59276.D0+rk*(-3602.D0+(rk-100.D0)*rk))))))) * den(8)
        ckplm(k,1957) = 1.85436672634549D-6*(1.00118304D10+rk*(-8.09825472D8+rk*(-2.07475876D8+rk&
                *(1.4846976D7+rk*(1.352829D6+rk*(-78732.D0+rk*(-2874.D0+(rk+108.D0)*rk))))))) * den(9)
        ckplm(k,1958) = 9.88995587384259D-8*(rk-9.D0)*(1.98028116D10+rk*(-1.64255778D9+rk&
                *(-2.88144584D8+rk*(2.3074993D7+rk*(1.041235D6+rk*(-78050.D0+rk*(-611.D0+37.D0*rk))))))) &
                * den(10)
        ckplm(k,1959) = 2.96698676215278D-7*(rk-10.D0)*(rk-9.D0)*(6.2450388D8+rk*(-5.7849966D7+rk&
                *(-5.400115D6+rk*(543495.D0+rk*(4790.D0+rk*(-849.D0+5.D0*rk)))))) * den(11)
        ckplm(k,1960) = -1.34863034643308D-8*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(-1.15383576D9+rk&
                *(1.21142612D8+rk*(3.74364D6+rk*(-606805.D0+rk*(7590.D0+203.D0*rk))))) * den(12)
        ckplm(k,1961) = -5.76337754885932D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(-1.95731718D9+rk*(2.23582432D8+rk*(-1.839311D6+rk*(-383050.D0+6893.D0*rk)))) * den(13)
        ckplm(k,1962) = 4.44603410912005D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(1.505679D7+rk*(-1.707751D6+rk*(42765.D0+211.D0*rk))) * den(14)
        ckplm(k,1963) = 2.57869978328963D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0)*(115135.D0+rk*(-11019.D0+233.D0*rk)) * den(15)
        ckplm(k,1964) = -3.33082055341577D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-9.D0)*(-2528.D0+145.D0*rk) * den(16)
        ckplm(k,1965) = 1.0991707826272D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-9.D0) * den(17)
        ckplm(k,1966) = 4.07100289861927D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0) * den(0)
        ckplm(k,1967) = -3.70091172601752D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(1100.D0+61.D0*rk) * den(1)
        ckplm(k,1968) = 2.38768498452743D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(766475.D0+4.D0*rk*(18933.D0+442.D0*rk)) * den(2)
        ckplm(k,1969) = -2.05834912459261D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(2.461665D7+rk*(3.127235D6+4.D0*rk*(28137.D0+220.D0*rk))) * den(3)
        ckplm(k,1970) = -1.44084438721483D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(-6.957615D7+rk&
                *(-9.757115D6+rk*(-339939.D0+2.D0*rk*(1384.D0+93.D0*rk)))) * den(4)
        ckplm(k,1971) = 7.49239081351712D-10*(rk+11.D0)*(rk+12.D0)*(-2.07440436D9+rk&
                *(-2.94546218D8+rk*(-5.561775D6+rk*(703705.D0+2.D0*rk*(13560.D0+59.D0*rk))))) * den(5)
        ckplm(k,1972) = 8.24162989486883D-9*(rk+11.D0)*(2.4315102D9+rk*(3.3518277D8+rk&
                *(-3.492907D6+rk*(-1.99695D6+rk*(-57805.D0+4.D0*rk*(255.D0+8.D0*rk)))))) * den(6)
        ckplm(k,1973) = 8.24162989486883D-9*(-2.70172188D10+rk*(-3.60752598D9+rk*(1.76635228D8+rk&
                *(3.7669067D7+rk*(509320.D0+rk*(-75175.D0+4.D0*rk*(-467.D0+2.D0*rk))))))) * den(7)
        ckplm(k,1974) = -2.06040747371721D-7*(-1.05271452D9+rk*(-3.6013842D7+rk*(1.6128503D7+rk&
                *(546094.D0+rk*(-66805.D0+rk*(-1973.D0+(rk+62.D0)*rk)))))) * den(8)
        ckplm(k,1975) = -2.06040747371721D-7*(1.00118304D9+rk*(-6.6375576D7+rk*(-1.411003D7+rk&
                *(792379.D0+rk*(56045.D0+(rk+28.D0)*(rk-83.D0)*rk))))) * den(9)
        ckplm(k,1976) = 8.24162989486883D-9*(2.327014404D10+rk*(-3.850191126D9+rk*(-6.7407509D7+rk&
                *(3.4917677D7+rk*(-856895.D0+rk*(-63799.D0+4.D0*rk*(481.D0+2.D0*rk))))))) * den(10)
        ckplm(k,1977) = 8.24162989486883D-9*(rk-10.D0)*(2.09477268D9+rk*(-3.36413862D8+rk&
                *(2.141393D6+rk*(1.75617D6+rk*(-62425.D0+4.D0*rk*(-207.D0+8.D0*rk)))))) * den(11)
        ckplm(k,1978) = 7.49239081351712D-10*(rk-10.D0)*(rk-11.D0)*(1.78609662D9+rk*(-2.81419443D8+rk&
                *(7.51135D6+rk*(596405.D0+2.D0*rk*(-13265.D0+59.D0*rk))))) * den(12)
        ckplm(k,1979) = -1.44084438721483D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(-6.0161556D7+rk&
                *(9.069677D6+rk*(-347127.D0+2.D0*rk*(-1012.D0+93.D0*rk)))) * den(13)
        ckplm(k,1980) = -2.05834912459261D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(-2.1601083D7+rk*(2.904779D6+4.D0*rk*(-27477.D0+220.D0*rk))) * den(14)
        ckplm(k,1981) = 2.38768498452743D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(692511.D0+4.D0*rk*(-18049.D0+442.D0*rk)) * den(15)
        ckplm(k,1982) = -3.70091172601752D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(-1039.D0+61.D0*rk) * den(16)
        ckplm(k,1983) = 4.07100289861927D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0) * den(17)
        ckplm(k,1984) = -1.45392960664974D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0) * den(0)
        ckplm(k,1985) = 1.3217541878634D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(1331.D0+75.D0*rk) * den(1)
        ckplm(k,1986) = -1.70548927466245D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(558899.D0+2.D0*rk*(29265.D0+743.D0*rk)) * den(2)
        ckplm(k,1987) = 2.94049874941802D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(1.06796415D8+(8542.D0&
                *rk*rk+669156.D0*rk+15228599.D0)*rk) * den(3)
        ckplm(k,1988) = 2.05834912459261D-11*(rk+12.D0)*(rk+13.D0)*(-3.5152194D8+rk*(-5.9474339D7+rk&
                *(-3.191716D6+rk*(-49666.D0+211.D0*rk)))) * den(4)
        ckplm(k,1989) = -5.3517077239408D-11*(rk+12.D0)*(-2.37451368D9+rk*(-4.42321529D8+rk&
                *(-2.3929155D7+rk*(-47840.D0+rk*(24645.D0+379.D0*rk))))) * den(5)
        ckplm(k,1990) = 5.88687849633488D-10*(-3.0593178D9+rk*(-6.04996974D8+rk*(-2.8364425D7+rk&
                *(1.1226D6+rk*(110195.D0+2.D0*(867.D0-5.D0*rk)*rk))))) * den(6)
        ckplm(k,1991) = 9.25080906566909D-10*(2.097459D9+rk*(2.46678858D8+rk*(-6.901193D6+rk&
                *(-1.48238D6+rk*(-22105.D0+2.D0*rk*(601.D0+9.D0*rk)))))) * den(7)
        ckplm(k,1992) = 1.15635113320864D-8*(-1.69329888D8+rk*(-5.543568D6+rk*(1.664632D6+rk&
                *(49797.D0+rk*(-3701.D0+(rk-81.D0)*rk))))) * den(8)
        ckplm(k,1993) = -1.15635113320864D-8*(-1.62175104D8+rk*(8.709048D6+rk*(1.49386D6+rk&
                *(-63771.D0+rk*(-3281.D0+(rk+87.D0)*rk))))) * den(9)
        ckplm(k,1994) = -9.25080906566909D-10*(1.84533804D9+rk*(-2.56128426D8+rk*(-2.598433D6+rk&
                *(1.3823D6+rk*(-27845.D0+2.D0*rk*(-547.D0+9.D0*rk)))))) * den(10)
        ckplm(k,1995) = 5.88687849633488D-10*(2.4836994D9+rk*(-5.45332374D8+rk*(3.1088545D7+rk&
                *(699360.D0+rk*(-101375.D0+2.D0*rk*(897.D0+5.D0*rk)))))) * den(11)
        ckplm(k,1996) = 5.3517077239408D-11*(rk-11.D0)*(1.9560492D9+rk*(-3.94703424D8+rk&
                *(2.3641555D7+rk*(-142630.D0+rk*(-22750.D0+379.D0*rk))))) * den(12)
        ckplm(k,1997) = 2.05834912459261D-11*(rk-11.D0)*(rk-12.D0)*(2.9518944D8-rk*(5.3240749D7+rk&
                *(-3.041452D6+rk*(50510.D0+211.D0*rk)))) * den(13)
        ckplm(k,1998) = 2.94049874941802D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(9.222843D7+(-8542.D0&
                *rk*rk+643530.D0*rk-13915913.D0)*rk) * den(14)
        ckplm(k,1999) = 1.70548927466245D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(501855.D0+2.D0*rk*(-27779.D0+743.D0*rk)) * den(15)
        ckplm(k,2000) = -1.3217541878634D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(-1256.D0+75.D0*rk) * den(16)
        ckplm(k,2001) = 1.45392960664974D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0) * den(17)
        ckplm(k,2002) = 5.01355036775773D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0) &
                * den(0)
        ckplm(k,2003) = -1.51925768719931D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(4752.D0+271.D0*rk) * den(1)
        ckplm(k,2004) = 1.17619949976721D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(394599.D0+rk&
                *(43111.D0+1157.D0*rk)) * den(2)
        ckplm(k,2005) = -1.76429924965081D-11*(rk+13.D0)*(rk+14.D0)*(1.01817D6+rk*(157527.D0+rk&
                *(7793.D0+121.D0*rk))) * den(3)
        ckplm(k,2006) = 2.05834912459261D-11*(rk+13.D0)*(2.321529D7+rk*(4.466584D6+rk*(296977.D0+rk&
                *(7626.D0+53.D0*rk)))) * den(4)
        ckplm(k,2007) = 1.60551231718224D-10*(-5.894364D7+rk*(-1.3142762D7+rk*(-1.017D6+rk&
                *(-28865.D0+rk*(-30.D0+7.D0*rk))))) * den(5)
        ckplm(k,2008) = 1.17737569926698D-9*(1.0416D7-rk*(-1.725461D6+rk*(-71715.D0+rk&
                *(1150.D0+(rk+105.D0)*rk)))) * den(6)
        ckplm(k,2009) = 5.60655094889036D-11*(-2.512692D8+rk*(-2.5166718D7+rk*(357365.D0+rk&
                *(78770.D0+(1000.D0-17.D0*rk)*rk)))) * den(7)
        ckplm(k,2010) = 1.05122830291694D-9*(1.4011452D7+rk*(415890.D0+rk*(-83075.D0+rk&
                *(-2023.D0+(rk+83.D0)*rk)))) * den(8)
        ckplm(k,2011) = 1.05122830291694D-9*(-1.3514592D7+rk*(575644.D0+rk*(76518.D0+rk&
                *(-2345.D0+(rk-78.D0)*rk)))) * den(9)
        ckplm(k,2012) = 5.60655094889036D-11*(2.2582287D8-rk*(2.5649223D7+rk*(127225.D0+rk&
                *(-74600.D0+rk*(1085.D0+17.D0*rk))))) * den(10)
        ckplm(k,2013) = 1.17737569926698D-9*(-8.7633D6-rk*(-1.578996D6+rk*(74545.D0+rk&
                *(740.D0+(rk-100.D0)*rk)))) * den(11)
        ckplm(k,2014) = 1.60551231718224D-10*(4.678905D7+rk*(-1.1195202D7+rk*(930655.D0+rk&
                *(-28675.D0+rk*(65.D0+7.D0*rk))))) * den(12)
        ckplm(k,2015) = 2.05834912459261D-11*(rk-12.D0)*(1.903811D7+rk*(-3.895296D6+rk*(274417.D0+rk&
                *(-7414.D0+53.D0*rk)))) * den(13)
        ckplm(k,2016) = 1.76429924965081D-11*(rk-12.D0)*(rk-13.D0)*(868315.D0+rk&
                *(-142304.D0+(7430.D0-121.D0*rk)*rk)) * den(14)
        ckplm(k,2017) = 1.17619949976721D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(352645.D0+rk&
                *(-40797.D0+1157.D0*rk)) * den(15)
        ckplm(k,2018) = -1.51925768719931D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(-4481.D0+271.D0*rk) * den(16)
        ckplm(k,2019) = 5.01355036775773D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0) &
                * den(17)
        ckplm(k,2020) = -1.67118345591924D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0) * den(0)
        ckplm(k,2021) = 1.51925768719931D-13*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(1859.D0+107.D0*rk) &
                * den(1)
        ckplm(k,2022) = -5.88099749883604D-14*(rk+14.D0)*(rk+15.D0)*(361491.D0+2.D0*rk&
                *(20387.D0+569.D0*rk)) * den(2)
        ckplm(k,2023) = 9.80166249806007D-14*(rk+14.D0)*(9.752145D6+rk*(1.602557D6+2.D0*rk&
                *(42894.D0+743.D0*rk))) * den(3)
        ckplm(k,2024) = 6.86116374864205D-13*(-4.244604D7-rk*(8.966249D6+rk*(680422.D0+rk&
                *(21586.D0+233.D0*rk)))) * den(4)
        ckplm(k,2025) = 1.78390257464693D-12*(2.791614D7+rk*(4.933391D6+rk*(289468.D0+rk&
                *(5914.D0+17.D0*rk)))) * den(5)
        ckplm(k,2026) = 5.88687849633488D-11*(-1.20057D6+rk*(-159519.D0+rk*(-5357.D0+2.D0*(rk+12.D0)&
                *rk))) * den(6)
        ckplm(k,2027) = -3.64425811677873D-11*(-2.37447D6+rk*(-193629.D0+rk*(1093.D0+2.D0*(rk+132.D0)&
                *rk))) * den(7)
        ckplm(k,2028) = -9.11064529194683D-11*(1.029D6+rk*(26128.D0+rk*(-3331.D0+(rk-58.D0)*rk))) &
                * den(8)
        ckplm(k,2029) = 9.11064529194683D-11*(999600.D0+rk*(-32612.D0+rk*(-3151.D0+(rk+62.D0)*rk))) &
                * den(9)
        ckplm(k,2030) = 3.64425811677873D-11*(-2.18001D6+rk*(195031.D0+rk*(313.D0+2.D0*(rk-128.D0)&
                *rk))) * den(10)
        ckplm(k,2031) = 5.88687849633488D-11*(1.04643D6+rk*(-148741.D0+rk*(5417.D0-2.D0*(rk-8.D0)&
                *rk))) * den(11)
        ckplm(k,2032) = 1.78390257464693D-12*(-2.326632D7+rk*(4.372129D6+rk&
                *(-271828.D0+(5846.D0-17.D0*rk)*rk))) * den(12)
        ckplm(k,2033) = 6.86116374864205D-13*(3.413886D7+rk*(-7.669231D6+rk*(617062.D0+rk&
                *(-20654.D0+233.D0*rk)))) * den(13)
        ckplm(k,2034) = 9.80166249806007D-14*(rk-13.D0)*(8.23389D6+rk*(-1.435439D6+2.D0&
                *(40665.D0-743.D0*rk)*rk)) * den(14)
        ckplm(k,2035) = 5.88099749883604D-14*(rk-13.D0)*(rk-14.D0)*(321855.D0+2.D0*rk&
                *(-19249.D0+569.D0*rk)) * den(15)
        ckplm(k,2036) = -1.51925768719931D-13*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(-1752.D0+107.D0*rk) &
                * den(16)
        ckplm(k,2037) = 1.67118345591924D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0) * den(17)
        ckplm(k,2038) = 5.39091437393304D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0) * den(0)
        ckplm(k,2039) = -4.90083124903003D-15*(rk+15.D0)*(rk+16.D0)*(2156.D0+125.D0*rk) * den(1)
        ckplm(k,2040) = 2.94049874941802D-14*(rk+15.D0)*(31311.D0+4.D0*rk*(905.D0+26.D0*rk)) * den(2)
        ckplm(k,2041) = 4.90083124903003D-14*(-971670.D0-rk*(167197.D0+4.D0*rk*(2367.D0+44.D0*rk))) &
                * den(3)
        ckplm(k,2042) = 2.40140731202472D-12*(49065.D0+(6.D0*rk*rk+380.D0*rk+7619.D0)*rk) * den(4)
        ckplm(k,2043) = -6.24365901126426D-12*(36060.D0+rk*(4729.D0+2.D0*(rk+93.D0)*rk)) * den(5)
        ckplm(k,2044) = 2.06040747371721D-10*(1695.D0+rk*(169.D0+4.D0*rk)) * den(6)
        ckplm(k,2045) = 1.40163773722259D-12*(-324990.D0+rk*(-20163.D0+4.D0*rk*(5.D0+2.D0*rk))) &
                * den(7)
        ckplm(k,2046) = 7.00818868611295D-12*(72975.D0-rk*(-1468.D0+(rk+108.D0)*rk)) * den(8)
        ckplm(k,2047) = 7.00818868611295D-12*(-71400.D0+rk*(1681.D0-(rk-105.D0)*rk)) * den(9)
        ckplm(k,2048) = 1.40163773722259D-12*(304815.D0+(8.D0*rk*rk+4.D0*rk-20179.D0)*rk) * den(10)
        ckplm(k,2049) = -2.06040747371721D-10*(1530.D0+rk*(-161.D0+4.D0*rk)) * den(11)
        ckplm(k,2050) = -6.24365901126426D-12*(-31515.D0+rk*(4363.D0+2.D0*(rk-90.D0)*rk)) * den(12)
        ckplm(k,2051) = 2.40140731202472D-12*(-41820.D0+(6.D0*rk*rk-362.D0*rk+6877.D0)*rk) * den(13)
        ckplm(k,2052) = 4.90083124903003D-14*(813765.D0+rk*(-148789.D0+4.D0*(2235.D0-44.D0*rk)*rk)) &
                * den(14)
        ckplm(k,2053) = 2.94049874941802D-14*(rk-14.D0)*(27795.D0+4.D0*rk*(-853.D0+26.D0*rk)) &
                * den(15)
        ckplm(k,2054) = 4.90083124903003D-15*(rk-14.D0)*(rk-15.D0)*(2031.D0-125.D0*rk) * den(16)
        ckplm(k,2055) = 5.39091437393304D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0) * den(17)
        ckplm(k,2056) = -1.68466074185407D-15*(rk+16.D0)*(rk+17.D0) * den(0)
        ckplm(k,2057) = 5.10503255107295D-17*(rk+16.D0)*(7425.D0+433.D0*rk) * den(1)
        ckplm(k,2058) = 1.22520781225751D-15*(-30900.D0-rk*(3643.D0+107.D0*rk)) * den(2)
        ckplm(k,2059) = 9.18905859193131D-14*(rk+20.D0)*(81.D0+5.D0*rk) * den(3)
        ckplm(k,2060) = -3.57352278575107D-14*(rk+24.D0)*(485.D0+29.D0*rk) * den(4)
        ckplm(k,2061) = 1.6724086637315D-13*(5288.D0+rk*(457.D0+9.D0*rk)) * den(5)
        ckplm(k,2062) = -1.22643302006977D-12*(rk+20.D0)*(rk+61.D0) * den(6)
        ckplm(k,2063) = 5.84015723842746D-14*(35580.D0+(rk+1489.D0)*rk) * den(7)
        ckplm(k,2064) = 1.09502948220515D-12*(-2208.D0+(rk-31.D0)*rk) * den(8)
        ckplm(k,2065) = 1.09502948220515D-12*(2176.D0-(rk+33.D0)*rk) * den(9)
        ckplm(k,2066) = 5.84015723842746D-14*(-34092.D0-(rk-1487.D0)*rk) * den(10)
        ckplm(k,2067) = 1.22643302006977D-12*(rk-19.D0)*(rk-60.D0) * den(11)
        ckplm(k,2068) = 1.6724086637315D-13*(-4840.D0+(439.D0-9.D0*rk)*rk) * den(12)
        ckplm(k,2069) = 3.57352278575107D-14*(rk-23.D0)*(-456.D0+29.D0*rk) * den(13)
        ckplm(k,2070) = -9.18905859193131D-14*(rk-19.D0)*(-76.D0+5.D0*rk) * den(14)
        ckplm(k,2071) = 1.22520781225751D-15*(27364.D0+rk*(-3429.D0+107.D0*rk)) * den(15)
        ckplm(k,2072) = 5.10503255107295D-17*(rk-15.D0)*(6992.D0-433.D0*rk) * den(16)
        ckplm(k,2073) = 1.68466074185407D-15*(rk-15.D0)*(rk-16.D0) * den(17)
        ckplm(k,2074) = 5.10503255107295D-17*(rk+17.D0) * den(0)
        ckplm(k,2075) = 5.10503255107295D-17*(-256.D0-15.D0*rk) * den(1)
        ckplm(k,2076) = 4.08402604085836D-16*(227.D0+13.D0*rk) * den(2)
        ckplm(k,2077) = 2.04201302042918D-15*(-202.D0-11.D0*rk) * den(3)
        ckplm(k,2078) = 7.14704557150213D-15*(181.D0+9.D0*rk) * den(4)
        ckplm(k,2079) = 1.85823184859056D-14*(-164.D0-7.D0*rk) * den(5)
        ckplm(k,2080) = 3.71646369718111D-14*(151.D0+5.D0*rk) * den(6)
        ckplm(k,2081) = 5.84015723842746D-14*(-142.D0-3.D0*rk) * den(7)
        ckplm(k,2082) = 7.30019654803432D-14*(rk+137.D0) * den(8)
        ckplm(k,2083) = 7.30019654803432D-14*(rk-136.D0) * den(9)
        ckplm(k,2084) = 5.84015723842746D-14*(139.D0-3.D0*rk) * den(10)
        ckplm(k,2085) = 3.71646369718111D-14*(-146.D0+5.D0*rk) * den(11)
        ckplm(k,2086) = 1.85823184859056D-14*(157.D0-7.D0*rk) * den(12)
        ckplm(k,2087) = 7.14704557150213D-15*(-172.D0+9.D0*rk) * den(13)
        ckplm(k,2088) = 2.04201302042918D-15*(191.D0-11.D0*rk) * den(14)
        ckplm(k,2089) = 4.08402604085836D-16*(-214.D0+13.D0*rk) * den(15)
        ckplm(k,2090) = 5.10503255107295D-17*(241.D0-15.D0*rk) * den(16)
        ckplm(k,2091) = 5.10503255107295D-17*(rk-16.D0) * den(17)
        ckplm(k,2092) = -1.50148016208028D-18 * den(0)
        ckplm(k,2093) = 2.55251627553648D-17 * den(1)
        ckplm(k,2094) = -2.04201302042918D-16 * den(2)
        ckplm(k,2095) = 1.02100651021459D-15 * den(3)
        ckplm(k,2096) = -3.57352278575107D-15 * den(4)
        ckplm(k,2097) = 9.29115924295277D-15 * den(5)
        ckplm(k,2098) = -1.85823184859056D-14 * den(6)
        ckplm(k,2099) = 2.92007861921373D-14 * den(7)
        ckplm(k,2100) = -3.65009827401716D-14 * den(8)
        ckplm(k,2101) = 3.65009827401716D-14 * den(9)
        ckplm(k,2102) = -2.92007861921373D-14 * den(10)
        ckplm(k,2103) = 1.85823184859056D-14 * den(11)
        ckplm(k,2104) = -9.29115924295277D-15 * den(12)
        ckplm(k,2105) = 3.57352278575107D-15 * den(13)
        ckplm(k,2106) = -1.02100651021459D-15 * den(14)
        ckplm(k,2107) = 2.04201302042918D-16 * den(15)
        ckplm(k,2108) = -2.55251627553648D-17 * den(16)
        ckplm(k,2109) = 1.50148016208028D-18 * den(17)
!    ckplm para l = 18
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)&
                *(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(17) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(18) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,2110) = 1.28089907112122D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2111) = 658748.093719482D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,2112) = 509032.617874146D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,2113) = 437877.520751953D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,2114) = 396354.652404785D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,2115) = 369931.008911133D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,2116) = 352667.561828613D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(6)
        ckplm(k,2117) = 341715.153076172D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(7)
        ckplm(k,2118) = 335613.09677124D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(8)
        ckplm(k,2119) = 333650.44708252D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk+9.D0)*rk * den(9)
        ckplm(k,2120) = 335613.09677124D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*rk * den(10)
        ckplm(k,2121) = 341715.153076172D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,2122) = 352667.561828613D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,2123) = 369931.008911133D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,2124) = 396354.652404785D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,2125) = 437877.520751953D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(15)
        ckplm(k,2126) = 509032.617874146D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(16)
        ckplm(k,2127) = 658748.093719482D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(17)
        ckplm(k,2128) = 1.28089907112122D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(18)
        ckplm(k,2129) = -134831.481170654D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2130) = -3852.32803344727D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-35.D0+16.D0*rk) * den(1)
        ckplm(k,2131) = -5953.59786987305D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-33.D0+7.D0*rk) * den(2)
        ckplm(k,2132) = -7682.06176757813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-31.D0+4.D0*rk) * den(3)
        ckplm(k,2133) = -4635.72692871094D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-58.D0+5.D0*rk) * den(4)
        ckplm(k,2134) = -2163.33923339844D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-135.D0+8.D0*rk) * den(5)
        ckplm(k,2135) = -12374.3004150391D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk-25.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(6)
        ckplm(k,2136) = -1998.33422851563D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-161.D0+4.D0*rk) * den(7)
        ckplm(k,2137) = -3925.29937744141D0*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-84.D0)*(rk+8.D0)*(rk+9.D0) * den(8)
        ckplm(k,2138) = 333650.44708252D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(9)
        ckplm(k,2139) = 3925.29937744141D0*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+85.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0) * den(10)
        ckplm(k,2140) = 1998.33422851563D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(165.D0+4.D0*rk) * den(11)
        ckplm(k,2141) = 12374.3004150391D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+26.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,2142) = 2163.33923339844D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(143.D0+8.D0*rk) * den(13)
        ckplm(k,2143) = 4635.72692871094D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(63.D0+5.D0*rk) * den(14)
        ckplm(k,2144) = 7682.06176757813D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(35.D0+4.D0*rk) * den(15)
        ckplm(k,2145) = 5953.59786987305D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(40.D0+7.D0*rk) * den(16)
        ckplm(k,2146) = 3852.32803344727D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(51.D0+16.D0*rk) * den(17)
        ckplm(k,2147) = 134831.481170654D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2148) = 6741.57405853272D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2149) = 1926.16401672363D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-14.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,2150) = 87.5529098510742D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk-78.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-11.D0+5.D0*rk) * den(2)
        ckplm(k,2151) = -45.1885986328125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+92.D0)&
                *(rk+9.D0)*(-31.D0+9.D0*rk) * den(3)
        ckplm(k,2152) = -13.6344909667969D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-13340.D0+rk*(2481.D0+71.D0*rk)) * den(4)
        ckplm(k,2153) = -12.7255249023438D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-18090.D0+rk*(2323.D0+107.D0*rk)) * den(5)
        ckplm(k,2154) = -181.975006103516D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1490.D0+rk*(131.D0+9.D0*rk)) * den(6)
        ckplm(k,2155) = -58.7745361328125D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-5152.D0+rk*(291.D0+31.D0*rk)) * den(7)
        ckplm(k,2156) = -11.5449981689453D0*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-28056.D0+rk*(841.D0+167.D0*rk)) * den(8)
        ckplm(k,2157) = -1962.6496887207D0*(rk*rk+rk-170.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(9)
        ckplm(k,2158) = -11.5449981689453D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(-28730.D0+rk*(-507.D0+167.D0*rk)) * den(10)
        ckplm(k,2159) = -58.7745361328125D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-5412.D0+rk*(-229.D0+31.D0*rk)) * den(11)
        ckplm(k,2160) = -181.975006103516D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1612.D0+rk*(-113.D0+9.D0*rk)) * den(12)
        ckplm(k,2161) = -12.7255249023438D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-20306.D0+rk*(-2109.D0+107.D0*rk)) * den(13)
        ckplm(k,2162) = -13.6344909667969D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-15750.D0+rk*(-2339.D0+71.D0*rk)) * den(14)
        ckplm(k,2163) = -45.1885986328125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-91.D0)&
                *(rk-9.D0)*(40.D0+9.D0*rk) * den(15)
        ckplm(k,2164) = 87.5529098510742D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk+79.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(16.D0+5.D0*rk) * den(16)
        ckplm(k,2165) = 1926.16401672363D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk+15.D0)*(rk-16.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(17)
        ckplm(k,2166) = 6741.57405853272D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2167) = -321.027336120606D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,2168) = 2889.24602508545D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(1)
        ckplm(k,2169) = 87.5529098510742D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-165.D0+(rk+32.D0)*rk) * den(2)
        ckplm(k,2170) = 3.76571655273438D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(18507.D0+2.D0*rk&
                *(-3874.D0+rk*(231.D0+16.D0*rk))) * den(3)
        ckplm(k,2171) = 1.94778442382813D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(63075.D0+rk&
                *(-17618.D0+rk*(327.D0+65.D0*rk))) * den(4)
        ckplm(k,2172) = 0.908966064453125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(199395.D0+rk&
                *(-38489.D0+rk*(-519.D0+128.D0*rk))) * den(5)
        ckplm(k,2173) = 8.66547648111979D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(27405.D0+rk&
                *(-3623.D0+rk*(-168.D0+11.D0*rk))) * den(6)
        ckplm(k,2174) = 2.09909057617188D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(135723.D0+2.D0*rk&
                *(-5764.D0+rk*(-531.D0+16.D0*rk))) * den(7)
        ckplm(k,2175) = 34.6349945068359D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(9186.D0+rk&
                *(-414.D0+(rk-79.D0)*rk)) * den(8)
        ckplm(k,2176) = -981.324844360352D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-340.D0+3.D0*(rk+1.D0)&
                *rk) * den(9)
        ckplm(k,2177) = -34.6349945068359D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(-9520.D0+rk&
                *(-253.D0+(rk+82.D0)*rk)) * den(10)
        ckplm(k,2178) = -2.09909057617188D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-146157.D0+2.D0*rk&
                *(-4654.D0+rk*(579.D0+16.D0*rk))) * den(11)
        ckplm(k,2179) = -8.66547648111979D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-30849.D0+rk&
                *(-3254.D0+rk*(201.D0+11.D0*rk))) * den(12)
        ckplm(k,2180) = -0.908966064453125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-237237.D0+rk&
                *(-37067.D0+rk*(903.D0+128.D0*rk))) * den(13)
        ckplm(k,2181) = -1.94778442382813D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-80955.D0+rk&
                *(-18077.D0+rk*(-132.D0+65.D0*rk))) * den(14)
        ckplm(k,2182) = -3.76571655273438D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-26685.D0+2.D0*rk&
                *(-4288.D0+rk*(-183.D0+16.D0*rk))) * den(15)
        ckplm(k,2183) = -87.5529098510742D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-196.D0+(rk-30.D0)*rk) * den(16)
        ckplm(k,2184) = 2889.24602508545D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(17)
        ckplm(k,2185) = 321.027336120606D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(18)
        ckplm(k,2186) = 14.5921516418457D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,2187) = -5.83686065673828D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+40.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(1)
        ckplm(k,2188) = -8.75529098510742D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+24.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) &
                * den(2)
        ckplm(k,2189) = -7.53143310546875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1674.D0+rk*(-370.D0+(rk-3.D0)&
                *rk)) * den(3)
        ckplm(k,2190) = -0.0649261474609375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.27803D6+rk*(475394.D0+rk*(-31213.D0+rk&
                *(-2294.D0+73.D0*rk)))) * den(4)
        ckplm(k,2191) = -0.302988688151042D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-470610.D0+rk*(121254.D0+rk*(-2321.D0+rk&
                *(-678.D0+5.D0*rk)))) * den(5)
        ckplm(k,2192) = 0.0787770589192708D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(2.63865D6+rk*(-466114.D0+rk*(-11191.D0+rk&
                *(2686.D0+19.D0*rk)))) * den(6)
        ckplm(k,2193) = 0.07633056640625D0*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(3.510444D6+rk*(-398620.D0+rk*(-30995.D0+rk&
                *(2270.D0+51.D0*rk)))) * den(7)
        ckplm(k,2194) = 0.17492421468099D0*(rk+10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.786356D6+rk*(-107658.D0+rk*(-19723.D0+rk&
                *(582.D0+31.D0*rk)))) * den(8)
        ckplm(k,2195) = 5.94742329915365D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(56100.D0+(rk*rk+rk-662.D0)*(rk+1.D0)*rk) &
                * den(9)
        ckplm(k,2196) = 0.17492421468099D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(1.87374D6+rk*(66590.D0+rk*(-21283.D0+rk&
                *(-458.D0+31.D0*rk)))) * den(10)
        ckplm(k,2197) = 0.07633056640625D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(3.87585D6+rk*(330024.D0+rk*(-37499.D0+rk&
                *(-2066.D0+51.D0*rk)))) * den(11)
        ckplm(k,2198) = 0.0787770589192708D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(3.090906D6+rk*(435750.D0+rk*(-19135.D0+rk&
                *(-2610.D0+19.D0*rk)))) * den(12)
        ckplm(k,2199) = -0.302988688151042D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-593502.D0+rk*(-123842.D0+rk*(-257.D0+rk&
                *(698.D0+5.D0*rk)))) * den(13)
        ckplm(k,2200) = -0.0649261474609375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1.78227D6+rk*(-530646.D0+rk*(-23893.D0+rk&
                *(2586.D0+73.D0*rk)))) * den(14)
        ckplm(k,2201) = -7.53143310546875D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2040.D0+rk*(-361.D0+(rk+6.D0)&
                *rk)) * den(15)
        ckplm(k,2202) = -8.75529098510742D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-23.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(16)
        ckplm(k,2203) = -5.83686065673828D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-39.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(17)
        ckplm(k,2204) = 14.5921516418457D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(18)
        ckplm(k,2205) = -0.634441375732422D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2206) = 0.0181268964494978D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(875.D0+32.D0*rk) &
                * den(1)
        ckplm(k,2207) = 0.0543806893484933D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-3525.D0+rk*(16.D0+9.D0*rk)) &
                * den(2)
        ckplm(k,2208) = 0.0116947719029018D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(147405.D0+2.D0*rk*(-7424.D0+rk&
                *(-489.D0+8.D0*rk))) * den(3)
        ckplm(k,2209) = -0.0141143798828125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(993975.D0+rk*(-211771.D0+rk*(-1141.D0+rk&
                *(1060.D0+7.D0*rk)))) * den(4)
        ckplm(k,2210) = -0.00094095865885417D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.191834D8+rk*(3.8372001D7+rk*(-1.99106D6+rk*(-245860.D0+rk&
                *(10865.D0+304.D0*rk))))) * den(5)
        ckplm(k,2211) = -0.00048929850260417D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-3.716664D8+rk*(8.2189631D7+rk*(121800.D0+rk*(-642835.D0+rk&
                *(4935.D0+719.D0*rk))))) * den(6)
        ckplm(k,2212) = -0.0190826416015625D0*(rk+10.D0)*(rk+11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.32027D7+rk*(1.878509D6+rk*(118825.D0+2.D0*rk*(-7845.D0+rk&
                *(-155.D0+8.D0*rk))))) * den(7)
        ckplm(k,2213) = -0.17492421468099D0*(rk+10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.75428D6+rk*(132534.D0+rk*(23260.D0+rk*(-1105.D0+(rk-70.D0)&
                *rk)))) * den(8)
        ckplm(k,2214) = 14.8685582478841D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(22440.D0+(rk*rk+rk-332.D0)*(rk+1.D0)*rk) * den(9)
        ckplm(k,2215) = 0.17492421468099D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(1.86252D6+rk*(82984.D0+rk*(-26145.D0+rk*(-815.D0+(rk+75.D0)&
                *rk)))) * den(10)
        ckplm(k,2216) = 0.0190826416015625D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(1.494702D7+rk*(1.595109D6+rk*(-163875.D0+2.D0*rk*(-7145.D0+rk&
                *(195.D0+8.D0*rk))))) * den(11)
        ckplm(k,2217) = 0.00048929850260417D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(4.5308718D8+rk*(8.0001381D7+rk*(-2.072725D6+rk*(-655385.D0+rk&
                *(-1340.D0+719.D0*rk))))) * den(12)
        ckplm(k,2218) = 0.00094095865885417D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.5929004D8+rk*(4.1574601D7+rk*(1.19133D6+rk*(-286280.D0+rk&
                *(-9345.D0+304.D0*rk))))) * den(13)
        ckplm(k,2219) = 0.0141143798828125D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.203552D6+rk*(206337.D0+rk*(-4279.D0+rk&
                *(-1032.D0+7.D0*rk)))) * den(14)
        ckplm(k,2220) = -0.0116947719029018D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-161259.D0+2.D0*rk*(-6422.D0+rk&
                *(513.D0+8.D0*rk))) * den(15)
        ckplm(k,2221) = -0.0543806893484933D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-3532.D0+rk*(2.D0+9.D0*rk)) &
                * den(16)
        ckplm(k,2222) = -0.0181268964494978D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-843.D0+32.D0*rk) &
                * den(17)
        ckplm(k,2223) = 0.634441375732422D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2224) = 0.0264350573221843D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2225) = -0.0135951723371233D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(70.D0+3.D0*rk) * den(1)
        ckplm(k,2226) = -0.00226586205618722D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-7122.D0+rk*(-235.D0+7.D0*rk)) * den(2)
        ckplm(k,2227) = 0.00194912865048363D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-97836.D0+rk*(2242.D0+rk*(543.D0+5.D0*rk))) * den(3)
        ckplm(k,2228) = 0.00035285949707031D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(5.39226D6+rk*(-549162.D0+rk*(-44221.D0+rk*(1722.D0+61.D0&
                *rk)))) * den(4)
        ckplm(k,2229) = 0.00023523966471354D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-7.500276D7+rk*(1.3925568D7+rk*(337832.D0+rk*(-98831.D0+rk*(-992.D0+83.D0&
                *rk))))) * den(5)
        ckplm(k,2230) = 1.56826443142361D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(1.01413053D11+rk*(-2.693311428D10+rk*(5.70012256D8+rk*(2.47665375D8+rk&
                *(-7.517285D6+rk*(-556035.D0+5569.D0*rk)))))) * den(6)
        ckplm(k,2231) = -0.00024464925130208D0*(rk+10.D0)*(rk+11.D0)*(rk-6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-9.6799248D8+rk*(1.65626856D8+rk*(8.054806D6+rk*(-1.779585D6+rk*(-19445.D0+rk&
                *(4029.D0+19.D0*rk)))))) * den(7)
        ckplm(k,2232) = -0.00504589080810547D0*(rk+10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-5.97168D7+rk*(5.42904D6+rk*(910598.D0+rk*(-60899.D0+rk*(-3929.D0+rk*(131.D0+3.D0&
                *rk)))))) * den(8)
        ckplm(k,2233) = -0.019062254163954D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-1.75032D7+(rk+1.D0)*rk*(311700.D0+(rk*rk+rk-1412.D0)*(rk+1.D0)*rk)) * den(9)
        ckplm(k,2234) = -0.00504589080810547D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(-6.41784D7+rk*(-3.4415D6+rk*(1.068456D6+rk*(43933.D0+rk*(-4539.D0+rk*(-113.D0+3.D0&
                *rk)))))) * den(10)
        ckplm(k,2235) = -0.00024464925130208D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-1.1238084D9+rk*(-1.442763D8+rk*(1.3236886D7+rk*(1.661895D6+rk*(-39305.D0+rk&
                *(-3915.D0+19.D0*rk)))))) * den(11)
        ckplm(k,2236) = 1.56826443142361D-6*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(1.2866155848D11+rk*(2.7302887116D10+rk*(-2.12443694D8+rk*(-2.72062785D8+rk&
                *(-4.653575D6+rk*(589449.D0+5569.D0*rk)))))) * den(12)
        ckplm(k,2237) = 0.00023523966471354D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(8.849274D7+rk*(1.2957794D7+rk*(-627543.D0+rk*(-94033.D0+rk*(1407.D0+83.D0&
                *rk))))) * den(13)
        ckplm(k,2238) = 0.00035285949707031D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(5.89554D6+rk*(455798.D0+rk*(-49021.D0+rk*(-1478.D0+61.D0&
                *rk)))) * den(14)
        ckplm(k,2239) = 0.00194912865048363D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(99540.D0+rk*(1171.D0+rk*(-528.D0+5.D0*rk))) * den(15)
        ckplm(k,2240) = -0.00226586205618722D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-6880.D0+rk*(249.D0+7.D0*rk)) * den(16)
        ckplm(k,2241) = -0.0135951723371233D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-67.D0+3.D0*rk) * den(17)
        ckplm(k,2242) = 0.0264350573221843D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2243) = -0.00105740229288737D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2244) = 0.00015105747041248D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+8.D0)*(rk+9.D0)*(343.D0+16.D0*rk) * den(1)
        ckplm(k,2245) = -0.00015105747041248D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0)*(7791.D0+(rk+440.D0)*rk) * den(2)
        ckplm(k,2246) = -0.00003898257300967D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(-457219.D0+rk*(-13599.D0+rk*(1411.D0+36.D0*rk))) * den(3)
        ckplm(k,2247) = -0.00002352396647135D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(9.17763D6+rk*(-193538.D0+rk*(-72563.D0+rk*(-688.D0+47.D0*rk)))) * den(4)
        ckplm(k,2248) = -0.00001097785101997D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-2.1084462D8+rk*(1.7258052D7+rk*(2.278535D6+rk*(-82610.D0+rk*(-5885.D0+8.D0&
                *rk))))) * den(5)
        ckplm(k,2249) = 0.00001097785101997D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.11118838D9+rk*(3.02338568D8+rk*(2.1683145D7+rk*(-2.84056D6+rk*(-84960.D0+rk*(5072.D0+75.D0&
                *rk)))))) * den(6)
        ckplm(k,2250) = 0.0002854241265191D0*(rk+10.D0)*(rk+11.D0)*(rk+8.D0)*(rk+9.D0)&
                *(7.7969268D8+rk*(-1.5593193D8+rk*(-5.285783D6+rk*(2.005261D6+rk*(-1730.D0+rk*(-6995.D0+rk&
                *(13.D0+4.D0*rk))))))) * den(7)
        ckplm(k,2251) = 0.00011213090684679D0*(rk+10.D0)*(rk-7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(2.63845296D9+7.D0*rk*(-4.0088232D7+rk*(-6.409128D6+rk*(562924.D0+rk*(35025.D0+rk&
                *(-1973.D0+(rk-57.D0)*rk)))))) * den(8)
        ckplm(k,2252) = -0.009531127081977D0*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-3.50064D7+7.D0&
                *(rk+1.D0)*rk*(104220.D0+(rk*rk+rk-632.D0)*(rk+1.D0)*rk)) * den(9)
        ckplm(k,2253) = -0.00011213090684679D0*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(-2.8705248D9+7.D0*rk*(-2.573082D7+rk*(7.868896D6+rk*(404269.D0+rk*(-44000.D0+rk&
                *(-1610.D0+(rk+64.D0)*rk)))))) * den(10)
        ckplm(k,2254) = -0.0002854241265191D0*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-9.2833884D8+rk*(-1.39372686D8+rk*(1.1241885D7+rk*(1.942111D6+rk*(-33300.D0+rk*(-6989.D0+rk&
                *(15.D0+4.D0*rk))))))) * den(11)
        ckplm(k,2255) = -0.00001097785101997D0*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-2.3890932D9+rk*(-2.50815348D8+rk*(2.964547D7+rk*(2.4515D6+rk*(-109195.D0+rk*(-4622.D0+75.D0&
                *rk)))))) * den(12)
        ckplm(k,2256) = 0.00001097785101997D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(2.2574742D8+rk*(1.2476732D7+rk*(-2.490975D6+rk*(-58990.D0+rk*(5925.D0+8.D0*rk))))) &
                * den(13)
        ckplm(k,2257) = 0.00002352396647135D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(9.29934D6+rk*(50664.D0+rk*(-70217.D0+rk*(876.D0+47.D0*rk)))) * den(14)
        ckplm(k,2258) = 0.00003898257300967D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(442245.D0+rk*(-16313.D0+rk*(-1303.D0+36.D0*rk))) * den(15)
        ckplm(k,2259) = 0.00015105747041248D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(7352.D0+(rk-438.D0)*rk) * den(16)
        ckplm(k,2260) = -0.00015105747041248D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-327.D0+16.D0*rk) * den(17)
        ckplm(k,2261) = 0.00105740229288737D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2262) = 0.00004066931895721D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+9.D0) * den(0)
        ckplm(k,2263) = -0.00001161980541634D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+9.D0)*(224.D0+11.D0*rk) * den(1)
        ckplm(k,2264) = 5.28172973470215D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+9.D0)*(144672.D0+rk*(10427.D0+133.D0*rk)) * den(2)
        ckplm(k,2265) = 2.72605405662046D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+9.D0)*(-5.319848D6+rk*(-355018.D0+rk*(3047.D0+307.D0*rk))) * den(3)
        ckplm(k,2266) = 8.22516310187209D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0)*(2.5656996D8+rk*(1.0043594D7+rk*(-1.241389D6+rk*(-52118.D0+25.D0*rk)))) * den(4)
        ckplm(k,2267) = -1.09668841358295D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(2.375784432D10+rk*(-1.7788272D7+rk*(-2.3497069D8+rk*(-4.441805D6+rk*(343990.D0+5777.D0&
                *rk))))) * den(5)
        ckplm(k,2268) = -7.84132215711806D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)&
                *(-3.69142032D9+rk*(1.63928272D8+rk*(5.1596702D7+rk*(-762515.D0+rk*(-185425.D0+rk&
                *(-797.D0+83.D0*rk)))))) * den(6)
        ckplm(k,2269) = -3.13652886284722D-6*(rk+10.D0)*(rk+11.D0)*(rk+25.D0)*(rk+9.D0)&
                *(3.80907072D8+rk*(-4.8018384D7+rk*(-4.42474D6+rk*(551220.D0+rk*(12593.D0+rk*(-1446.D0+5.D0&
                *rk)))))) * den(7)
        ckplm(k,2270) = 3.92066107855903D-7*(rk+10.D0)*(rk+9.D0)*(7.4081655648D11+rk&
                *(-9.0288947952D10+rk*(-1.3719116172D10+rk*(1.517669788D9+rk*(8.8252043D7+rk*(-7.432568D6+rk&
                *(-216298.D0+rk*(9292.D0+107.D0*rk)))))))) * den(8)
        ckplm(k,2271) = 0.0000666512383355D0*(rk-8.D0)*(rk+9.D0)*(5.0059152D9+(rk+1.D0)*rk&
                *(-1.1959848D8+(rk+1.D0)*rk*(910732.D0+(rk*rk+rk-2308.D0)*(rk+1.D0)*rk))) * den(9)
        ckplm(k,2272) = 3.92066107855903D-7*(rk-8.D0)*(rk-9.D0)*(8.159641776D11+rk*(5.868651528D10+rk&
                *(-1.7671724204D10+rk*(-1.094981124D9+rk*(1.21852683D8+rk*(5.94564D6+rk*(-278346.D0+rk&
                *(-8436.D0+107.D0*rk)))))))) * den(10)
        ckplm(k,2273) = -3.13652886284722D-6*(rk-10.D0)*(rk-24.D0)*(rk-8.D0)*(rk-9.D0)&
                *(4.2396354D8+rk*(3.7572876D7+rk*(-5.988307D6+rk*(-486288.D0+rk*(19898.D0+rk*(1476.D0+5.D0&
                *rk)))))) * den(11)
        ckplm(k,2274) = -7.84132215711806D-7*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-3.80317392D9+rk*(-5.918454D7+rk*(5.2780912D7+rk*(30445.D0+rk*(-180195.D0+rk*(1295.D0+83.D0&
                *rk)))))) * den(12)
        ckplm(k,2275) = -1.09668841358295D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-2.354544192D10+rk*(4.37480618D8+rk*(2.19639105D8+rk*(-5.759995D6+rk*(-315105.D0+5777.D0&
                *rk))))) * den(13)
        ckplm(k,2276) = 8.22516310187209D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(2.4533712D8+rk*(-1.2369918D7+rk*(-1.084885D6+rk*(52218.D0+25.D0*rk)))) * den(14)
        ckplm(k,2277) = 2.72605405662046D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0)*(4.96209D6+rk*(-360191.D0+rk*(-2126.D0+307.D0*rk))) * den(15)
        ckplm(k,2278) = 5.28172973470215D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-8.D0)*(rk-9.D0)*(134378.D0+rk*(-10161.D0+133.D0*rk)) * den(16)
        ckplm(k,2279) = -0.00001161980541634D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-8.D0)*(rk-9.D0)*(-213.D0+11.D0*rk) * den(17)
        ckplm(k,2280) = 0.00004066931895721D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2281) = -1.50627107248913D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0) * den(0)
        ckplm(k,2282) = 3.87326847211491D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(315.D0+16.D0*rk) * den(1)
        ckplm(k,2283) = -3.5211531564681D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(127809.D0+rk*(10600.D0+191.D0*rk)) * den(2)
        ckplm(k,2284) = -1.51447447590026D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(-6.930081D6+rk*(-648379.D0+rk*(-11325.D0+148.D0*rk))) * den(3)
        ckplm(k,2285) = 5.48344206791473D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(-3.3112287D8+rk*(-2.8162038D7+rk*(308581.D0+rk*(58848.D0+719.D0*rk)))) * den(4)
        ckplm(k,2286) = 1.64503262037442D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(1.56236346D9+rk*(9.9678624D7+rk*(-8.713951D6+rk*(-623342.D0+rk*(-839.D0+248.D0*rk))))) &
                * den(5)
        ckplm(k,2287) = -1.56826443142361D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(2.0170773D9+rk&
                *(7.436154D7+rk*(-2.2113941D7+rk*(-993660.D0+rk*(42520.D0+(rk+1920.D0)*rk))))) * den(6)
        ckplm(k,2288) = -9.40958658854167D-7*(rk+10.D0)*(rk+11.D0)*(-3.72687588D9+rk*(-3.7246806D7+rk&
                *(6.0555555D7+rk*(1.299481D6+rk*(-261000.D0+rk*(-6779.D0+rk*(225.D0+4.D0*rk))))))) * den(7)
        ckplm(k,2289) = -3.52859497070313D-6*(rk+10.D0)*(1.010017008D10+rk*(-1.25973432D8+rk&
                *(-2.15607864D8+rk*(145212.D0+rk*(1.454189D6+rk*(11132.D0+rk*(-3126.D0+(rk-32.D0)*rk))))))) &
                * den(8)
        ckplm(k,2290) = 0.00009997685750326D0*(3.3372768D9+(rk+1.D0)*rk*(-8.9980608D7+(rk+1.D0)*rk&
                *(826292.D0+(rk+1.D0)*rk*(-2920.D0+3.D0*(rk+1.D0)*rk)))) * den(9)
        ckplm(k,2291) = 3.52859497070313D-6*(rk-9.D0)*(1.00118304D10+rk*(-2.9993536D8+rk&
                *(-2.07475876D8+rk*(5.49888D6+rk*(1.352829D6+rk*(-29160.D0+rk*(-2874.D0+(rk+40.D0)*rk))))))) &
                * den(10)
        ckplm(k,2292) = 9.40958658854167D-7*(rk-10.D0)*(rk-9.D0)*(3.630627D9+rk*(-1.5345069D8+rk&
                *(-5.5162193D7+rk*(2.271331D6+rk*(223870.D0+rk*(-8045.D0+rk*(-197.D0+4.D0*rk))))))) * den(11)
        ckplm(k,2293) = 1.56826443142361D-7*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(1.92163608D9+rk&
                *(-1.15447956D8+rk*(-1.8897026D7+rk*(1.14456D6+rk*(32935.D0+(rk-1914.D0)*rk))))) * den(12)
        ckplm(k,2294) = -1.64503262037442D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(-1.45459314D9+rk*(1.15241096D8+rk*(6.851439D6+rk*(-617506.D0+rk*(2079.D0+248.D0*rk))))) &
                * den(13)
        ckplm(k,2295) = -5.48344206791473D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(-3.0271038D8+rk*(2.8605532D7+rk*(136351.D0+rk*(-55972.D0+719.D0*rk)))) * den(14)
        ckplm(k,2296) = 1.51447447590026D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0)*(6.293175D6+rk*(-625285.D0+rk*(11769.D0+148.D0*rk))) * den(15)
        ckplm(k,2297) = 3.5211531564681D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-9.D0)*(117400.D0+rk*(-10218.D0+191.D0*rk)) * den(16)
        ckplm(k,2298) = -3.87326847211491D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-9.D0)*(-299.D0+16.D0*rk) * den(17)
        ckplm(k,2299) = 1.50627107248913D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-9.D0) * den(18)
        ckplm(k,2300) = 5.37953954460404D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0) * den(0)
        ckplm(k,2301) = -2.15181581784162D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(250.D0+13.D0*rk) * den(1)
        ckplm(k,2302) = 9.78098099018916D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(249150.D0+rk*(22621.D0+479.D0*rk)) * den(2)
        ckplm(k,2303) = -2.52412412650043D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(2.71994D6+rk*(308854.D0+rk*(9469.D0+47.D0*rk))) * den(3)
        ckplm(k,2304) = -1.08798453728467D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(-1.2847899D9+rk*(-1.5505043D8+rk*(-3.855047D6+rk*(96206.D0+2951.D0*rk)))) * den(4)
        ckplm(k,2305) = -1.52317835219853D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(1.48229676D9+rk&
                *(1.71839208D8+rk*(247304.D0+rk*(-502025.D0+rk*(-12992.D0+5.D0*rk))))) * den(5)
        ckplm(k,2306) = 4.35629008728781D-9*(rk+11.D0)*(rk+12.D0)*(7.0550886D9+rk*(7.3845896D8+rk&
                *(-3.0418136D7+rk*(-4.894525D6+rk*(-68465.D0+rk*(3785.D0+61.D0*rk)))))) * den(6)
        ckplm(k,2307) = 2.48930862130732D-9*(rk+11.D0)*(-1.469967408D11+rk*(-1.354966488D10+rk&
                *(1.407450268D9+rk*(1.51648252D8+rk*(-1.370705D6+rk*(-348425.D0+rk*(-3983.D0+73.D0*rk))))))) &
                * den(7)
        ckplm(k,2308) = -1.55581788831708D-9*(rk+24.D0)*(-1.040935896D11+rk*(-4.29971922D9+rk&
                *(1.751289292D9+rk*(6.5001959D7+rk*(-8.60024D6+rk*(-233530.D0+rk*(11908.D0+71.D0*rk))))))) &
                * den(8)
        ckplm(k,2309) = -2.64489041013903D-7*(1.401656256D10+(rk+1.D0)*rk*(-2.65502304D8+(rk+1.D0)*rk&
                *(1.586308D6+(rk*rk+rk-3100.D0)*(rk+1.D0)*rk))) * den(9)
        ckplm(k,2310) = -1.55581788831708D-9*(rk-23.D0)*(9.811593792D10+rk*(-7.574129568D9+rk&
                *(-1.507194404D9+rk*(9.6831944D7+rk*(7.256455D6+rk*(-303487.D0+rk*(-11411.D0+71.D0*rk))))))) &
                * den(10)
        ckplm(k,2311) = 2.48930862130732D-9*(rk-10.D0)*(1.3219230024D11+rk*(-1.5905855556D10+rk&
                *(-9.47704254D8+rk*(1.53729037D8+rk*(-309120.D0+rk*(-322994.D0+rk*(4494.D0+73.D0*rk))))))) &
                * den(11)
        ckplm(k,2312) = 4.35629008728781D-9*(rk-10.D0)*(rk-11.D0)*(6.29103384D9+rk*(-7.84904076D8+rk&
                *(-1.6182286D7+rk*(4.584035D6+rk*(-86475.D0+rk*(-3419.D0+61.D0*rk)))))) * den(12)
        ckplm(k,2313) = -1.52317835219853D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(-1.311193884D9+rk&
                *(1.69890518D8+rk*(-1.675377D6+rk*(-450007.D0+rk*(13017.D0+5.D0*rk))))) * den(13)
        ckplm(k,2314) = -1.08798453728467D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(-1.133687772D9+rk*(1.47063522D8+rk*(-4.125959D6+rk*(-84402.D0+2951.D0*rk)))) * den(14)
        ckplm(k,2315) = -2.52412412650043D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(-2.420508D6+rk*(290057.D0+rk*(-9328.D0+47.D0*rk))) * den(15)
        ckplm(k,2316) = 9.78098099018916D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(227008.D0+rk*(-21663.D0+479.D0*rk)) * den(16)
        ckplm(k,2317) = -2.15181581784162D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(-237.D0+13.D0*rk) * den(17)
        ckplm(k,2318) = 5.37953954460404D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0) * den(18)
        ckplm(k,2319) = -1.85501363607036D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0) * den(0)
        ckplm(k,2320) = 3.71002727214072D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(605.D0+32.D0*rk) * den(1)
        ckplm(k,2321) = -3.37275206558247D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(363363.D0+rk*(35120.D0+817.D0*rk)) * den(2)
        ckplm(k,2322) = 2.17596907456934D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(1.8886285D7+2.D0*rk*(1.216332D6+rk*(47251.D0+504.D0*rk))) * den(3)
        ckplm(k,2323) = 2.17596907456934D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(-4.48274145D8+rk&
                *(-6.6698267D7+rk*(-2.976281D6+rk*(-27388.D0+491.D0*rk)))) * den(4)
        ckplm(k,2324) = -1.52317835219853D-10*(rk+12.D0)*(rk+13.D0)*(-1.17455184D9+rk&
                *(-1.86208689D8+rk*(-7.428224D6+rk*(121532.D0+rk*(10859.D0+112.D0*rk))))) * den(5)
        ckplm(k,2325) = -4.35629008728781D-9*(rk+12.D0)*(6.170934D8+rk*(9.9885965D7+rk*(2.580267D6+rk&
                *(-292945.D0+rk*(-15750.D0+rk*(-100.D0+3.D0*rk)))))) * den(6)
        ckplm(k,2326) = 6.2232715532683D-10*(5.52089538D10+rk*(9.015097446D9+rk*(4.8187405D7+rk&
                *(-5.5601126D7+rk*(-2.325995D6+2.D0*rk*(18032.D0+rk*(1295.D0+8.D0*rk))))))) * den(7)
        ckplm(k,2327) = 1.71139967714878D-8*(rk-12.D0)*(1.6886142D8+rk*(2.6961066D7+rk*(593164.D0+rk&
                *(-81495.D0+rk*(-3485.D0+(rk+9.D0)*rk))))) * den(8)
        ckplm(k,2328) = -1.45468972557647D-6*(-2.3167872D7+(rk+1.D0)*rk*(300312.D0+(rk*rk+rk-1100.D0)&
                *(rk+1.D0)*rk)) * den(9)
        ckplm(k,2329) = -1.71139967714878D-8*(rk+13.D0)*(1.4257152D8+rk*(-2.5544232D7+rk&
                *(816664.D0+rk*(67485.D0+rk*(-3515.D0+(rk-3.D0)*rk))))) * den(10)
        ckplm(k,2330) = -6.2232715532683D-10*(-4.62952854D10+rk*(8.76138813D9+rk*(-2.00712687D8+rk&
                *(-4.5987746D7+rk*(2.468025D6+2.D0*rk*(10430.D0+rk*(-1239.D0+8.D0*rk))))))) * den(11)
        ckplm(k,2331) = 4.35629008728781D-9*(rk-11.D0)*(5.20065D8+rk*(-9.3909078D7+rk*(3.365647D6+rk&
                *(231005.D0+rk*(-15205.D0+rk*(118.D0+3.D0*rk)))))) * den(12)
        ckplm(k,2332) = 1.52317835219853D-10*(rk-11.D0)*(rk-12.D0)*(9.9588216D8+rk*(-1.71030521D8+rk&
                *(7.728786D6+rk*(79216.D0+rk*(-10299.D0+112.D0*rk))))) * den(13)
        ckplm(k,2333) = -2.17596907456934D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(-3.8452428D8+rk&
                *(6.0829833D7+rk*(-2.891171D6+rk*(29352.D0+491.D0*rk)))) * den(14)
        ckplm(k,2334) = -2.17596907456934D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(-1.6547115D7+2.D0*rk*(1.123342D6+rk*(-45739.D0+504.D0*rk))) * den(15)
        ckplm(k,2335) = 3.37275206558247D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(329060.D0+rk*(-33486.D0+817.D0*rk)) * den(16)
        ckplm(k,2336) = -3.71002727214072D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(-573.D0+32.D0*rk) * den(17)
        ckplm(k,2337) = 1.85501363607036D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0) * den(18)
        ckplm(k,2338) = 6.18337878690119D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0) * den(0)
        ckplm(k,2339) = -1.59001168806031D-10*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(56.D0+3.D0*rk) * den(1)
        ckplm(k,2340) = 2.40910861827319D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(239448.D0+rk*(24217.D0+599.D0*rk)) * den(2)
        ckplm(k,2341) = -2.07235149958984D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(1.0968234D7+rk&
                *(1.54307D6+rk*(68613.D0+937.D0*rk))) * den(3)
        ckplm(k,2342) = 1.08798453728467D-11*(rk+13.D0)*(rk+14.D0)*(5.729931D7+rk*(9.814446D6+rk&
                *(566293.D0+rk*(11814.D0+47.D0*rk)))) * den(4)
        ckplm(k,2343) = 4.35193814913867D-12*(rk+13.D0)*(-2.96878932D9+rk*(-5.75136822D8+rk&
                *(-3.674345D7+rk*(-688735.D0+rk*(10250.D0+307.D0*rk))))) * den(5)
        ckplm(k,2344) = -1.03721192554472D-10*(-2.06177076D9+rk*(-4.34095566D8+rk*(-2.773384D7+rk&
                *(-117405.D0+rk*(44705.D0+rk*(1221.D0+5.D0*rk)))))) * den(6)
        ckplm(k,2345) = -1.24465431065366D-9*(1.9633068D8+rk*(2.788902D7+rk*(327874.D0+rk&
                *(-95325.D0+rk*(-3425.D0+(rk+15.D0)*rk))))) * den(7)
        ckplm(k,2346) = 6.66779094993032D-11*(3.86656452D9+rk*(2.63347896D8+rk*(-2.3472488D7+rk&
                *(-1.602595D6+rk*(21785.D0+rk*(1639.D0+3.D0*rk)))))) * den(8)
        ckplm(k,2347) = 1.25947162387573D-9*(-2.0271888D8+(rk+1.D0)*rk*(1.726932D6+(rk*rk+rk-3518.D0)&
                *(rk+1.D0)*rk)) * den(9)
        ckplm(k,2348) = 6.66779094993032D-11*(3.58136688D9+rk*(-3.05406124D8+rk*(-1.8550338D7+rk&
                *(1.673405D6+rk*(13635.D0+rk*(-1621.D0+3.D0*rk)))))) * den(10)
        ckplm(k,2349) = -1.24465431065366D-9*(1.6886142D8+rk*(-2.6961066D7+rk*(593164.D0+rk&
                *(81495.D0+rk*(-3485.D0+(rk-9.D0)*rk))))) * den(11)
        ckplm(k,2350) = -1.03721192554472D-10*(-1.65524814D9+rk*(3.79152846D8+rk*(-2.712553D7+rk&
                *(284115.D0+rk*(38675.D0+rk*(-1191.D0+5.D0*rk)))))) * den(12)
        ckplm(k,2351) = 4.35193814913867D-12*(rk-12.D0)*(2.42969727D9+rk*(-5.03755592D8+rk&
                *(3.4618815D7+rk*(-726665.D0+rk*(-8715.D0+307.D0*rk))))) * den(13)
        ckplm(k,2352) = 1.08798453728467D-11*(rk-12.D0)*(rk-13.D0)*(4.803939D7+rk*(-8.717114D6+rk&
                *(531133.D0+rk*(-11626.D0+47.D0*rk)))) * den(14)
        ckplm(k,2353) = -2.07235149958984D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(-9.49284D6+rk&
                *(1.408655D6+rk*(-65802.D0+937.D0*rk))) * den(15)
        ckplm(k,2354) = 2.40910861827319D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(215830.D0+rk*(-23019.D0+599.D0*rk)) * den(16)
        ckplm(k,2355) = -1.59001168806031D-10*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(-53.D0+3.D0*rk) * den(17)
        ckplm(k,2356) = 6.18337878690119D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0) * den(18)
        ckplm(k,2357) = -1.99463831835522D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0) &
                * den(0)
        ckplm(k,2358) = 2.84948331193603D-13*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(1183.D0+64.D0*rk) * den(1)
        ckplm(k,2359) = -7.77131812346191D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(329043.D0+rk&
                *(34432.D0+889.D0*rk)) * den(2)
        ckplm(k,2360) = 5.18087874897461D-13*(rk+14.D0)*(rk+15.D0)*(2.259699D6+2.D0*rk*(169594.D0+rk&
                *(8235.D0+128.D0*rk))) * den(3)
        ckplm(k,2361) = -3.62661512428223D-12*(rk+14.D0)*(1.0147605D7+rk*(1.923203D6+rk*(129119.D0+rk&
                *(3532.D0+31.D0*rk)))) * den(4)
        ckplm(k,2362) = -7.25323024856445D-13*(-1.18245582D9+rk*(-2.63763447D8+rk*(-2.137165D7+rk&
                *(-728960.D0+rk*(-7925.D0+32.D0*rk))))) * den(5)
        ckplm(k,2363) = 1.03721192554472D-10*(-1.157142D7+rk*(-2.029247D6+rk*(-108528.D0+rk&
                *(-929.D0+(rk+63.D0)*rk)))) * den(6)
        ckplm(k,2364) = -4.04512650962439D-9*(-363510.D0+rk*(-43529.D0+rk*(-617.D0+2.D0*(rk+37.D0)&
                *rk))) * den(7)
        ckplm(k,2365) = -9.6312535943438D-11*(rk+21.D0)*(801780.D0+rk*(8854.D0+rk*(-3244.D0+(rk-1.D0)&
                *rk))) * den(8)
        ckplm(k,2366) = 8.18656555519223D-9*(199920.D0+(rk*rk+rk-1052.D0)*(rk+1.D0)*rk) * den(9)
        ckplm(k,2367) = 9.6312535943438D-11*(rk-20.D0)*(789684.D0+rk*(-15335.D0+rk&
                *(-3235.D0+(rk+5.D0)*rk))) * den(10)
        ckplm(k,2368) = -4.04512650962439D-9*(-320670.D0+rk*(42081.D0+rk*(-827.D0+2.D0*(rk-33.D0)&
                *rk))) * den(11)
        ckplm(k,2369) = -1.03721192554472D-10*(9.64971D6+rk*(-1.815225D6+rk*(105373.D0+rk&
                *(-1171.D0+(rk-58.D0)*rk)))) * den(12)
        ckplm(k,2370) = 7.25323024856445D-13*(9.3934302D8+rk*(-2.23175167D8+rk*(1.923264D7+rk&
                *(-696940.D0+rk*(8085.D0+32.D0*rk))))) * den(13)
        ckplm(k,2371) = 3.62661512428223D-12*(rk-13.D0)*(8.35002D6+rk*(-1.675437D6+rk*(118709.D0+rk&
                *(-3408.D0+31.D0*rk)))) * den(14)
        ckplm(k,2372) = -5.18087874897461D-13*(rk-13.D0)*(rk-14.D0)*(-1.936725D6+2.D0*rk&
                *(153508.D0+rk*(-7851.D0+128.D0*rk))) * den(15)
        ckplm(k,2373) = 7.77131812346191D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(295500.D0+rk&
                *(-32654.D0+889.D0*rk)) * den(16)
        ckplm(k,2374) = -2.84948331193603D-13*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(-1119.D0+64.D0*rk) * den(17)
        ckplm(k,2375) = 1.99463831835522D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0) &
                * den(18)
        ckplm(k,2376) = 6.23324474486008D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0) * den(0)
        ckplm(k,2377) = -3.56185413992004D-15*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(3430.D0+187.D0*rk) &
                * den(1)
        ckplm(k,2378) = 4.85707382716369D-16*(rk+15.D0)*(rk+16.D0)*(2.208822D6+rk*(237305.D0+6323.D0&
                *rk)) * den(2)
        ckplm(k,2379) = -1.29521968724365D-14*(rk+15.D0)*(4.36296D6+rk*(687982.D0+rk*(35541.D0+599.D0&
                *rk))) * den(3)
        ckplm(k,2380) = 2.26663445267639D-14*(8.915592D7+rk*(1.8231202D7+rk*(1.354561D6+rk&
                *(42878.D0+479.D0*rk)))) * den(4)
        ckplm(k,2381) = -3.17328823374695D-13*(rk+24.D0)*(494370.D0+rk*(66557.D0+rk*(2526.D0+19.D0&
                *rk))) * den(5)
        ckplm(k,2382) = -4.53780217425813D-12*(-1.2732D6+rk*(-178210.D0+rk*(-7393.D0+(rk-62.D0)&
                *rk))) * den(6)
        ckplm(k,2383) = 7.77908944158537D-12*(rk+16.D0)*(-60930.D0+rk*(-2081.D0+(rk+50.D0)*rk)) &
                * den(7)
        ckplm(k,2384) = 2.31520519094803D-13*(3.775968D7+rk*(1.821654D6+rk*(-69413.D0+rk&
                *(-2742.D0+5.D0*rk)))) * den(8)
        ckplm(k,2385) = -7.8716976492233D-12*(1.1424D6+(rk*rk+rk-3362.D0)*(rk+1.D0)*rk) * den(9)
        ckplm(k,2386) = 2.31520519094803D-13*(3.587136D7+rk*(-1.952234D6+rk*(-61157.D0+rk&
                *(2762.D0+5.D0*rk)))) * den(10)
        ckplm(k,2387) = 7.77908944158537D-12*(rk-15.D0)*(58800.D0+rk*(-2178.D0+(rk-47.D0)*rk)) &
                * den(11)
        ckplm(k,2388) = -4.53780217425813D-12*(-1.10232D6+rk*(163614.D0+rk*(-7201.D0+(rk+66.D0)&
                *rk))) * den(12)
        ckplm(k,2389) = -3.17328823374695D-13*(rk-23.D0)*(-430320.D0+rk*(61562.D0+rk*(-2469.D0+19.D0&
                *rk))) * den(13)
        ckplm(k,2390) = 2.26663445267639D-14*(7.223688D7+rk*(-1.5648798D7+rk*(1.228801D6+rk&
                *(-40962.D0+479.D0*rk)))) * den(14)
        ckplm(k,2391) = -1.29521968724365D-14*(rk-14.D0)*(-3.70992D6+rk*(618697.D0+rk&
                *(-33744.D0+599.D0*rk))) * den(15)
        ckplm(k,2392) = 4.85707382716369D-16*(rk-14.D0)*(rk-15.D0)*(1.97784D6+rk*(-224659.D0+6323.D0&
                *rk)) * den(16)
        ckplm(k,2393) = -3.56185413992004D-15*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(-3243.D0+187.D0*rk) &
                * den(17)
        ckplm(k,2394) = 6.23324474486008D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0) * den(18)
        ckplm(k,2395) = -1.88886204389699D-15*(rk+16.D0)*(rk+17.D0)*(rk+18.D0) * den(0)
        ckplm(k,2396) = 4.85707382716369D-16*(rk+16.D0)*(rk+17.D0)*(875.D0+48.D0*rk) * den(1)
        ckplm(k,2397) = -1.61902460905457D-16*(rk+16.D0)*(264225.D0+rk*(28984.D0+791.D0*rk)) * den(2)
        ckplm(k,2398) = 2.15869947873942D-15*(1.19046D6+rk*(195091.D0+rk*(10551.D0+188.D0*rk))) &
                * den(3)
        ckplm(k,2399) = -1.1333172263382D-13*(61080.D0+rk*(9176.D0+rk*(447.D0+7.D0*rk))) * den(4)
        ckplm(k,2400) = 1.1333172263382D-13*(126792.D0+rk*(16549.D0+rk*(669.D0+8.D0*rk))) * den(5)
        ckplm(k,2401) = -1.96441652565287D-14*(1.22796D6+rk*(128423.D0+rk*(3642.D0+19.D0*rk))) &
                * den(6)
        ckplm(k,2402) = -1.17864991539172D-13*(-287580.D0+rk*(-20989.D0+rk*(-225.D0+4.D0*rk))) &
                * den(7)
        ckplm(k,2403) = 6.94561557284408D-13*(-58800.D0+rk*(-2178.D0+(rk+47.D0)*rk)) * den(8)
        ckplm(k,2404) = -1.96792441230582D-11*(-2176.D0+3.D0*(rk+1.D0)*rk) * den(9)
        ckplm(k,2405) = -6.94561557284408D-13*(56576.D0+rk*(-2269.D0+(rk-44.D0)*rk)) * den(10)
        ckplm(k,2406) = 1.17864991539172D-13*(266820.D0+rk*(-20527.D0+rk*(237.D0+4.D0*rk))) * den(11)
        ckplm(k,2407) = 1.96441652565287D-14*(-1.10316D6+rk*(121196.D0+rk*(-3585.D0+19.D0*rk))) &
                * den(12)
        ckplm(k,2408) = -1.1333172263382D-13*(-110904.D0+rk*(15235.D0+rk*(-645.D0+8.D0*rk))) * den(13)
        ckplm(k,2409) = 1.1333172263382D-13*(-52344.D0+rk*(8303.D0+rk*(-426.D0+7.D0*rk))) * den(14)
        ckplm(k,2410) = -2.15869947873942D-15*(-1.005732D6+rk*(174553.D0+rk*(-9987.D0+188.D0*rk))) &
                * den(15)
        ckplm(k,2411) = 1.61902460905457D-16*(rk-15.D0)*(236032.D0+rk*(-27402.D0+791.D0*rk)) * den(16)
        ckplm(k,2412) = -4.85707382716369D-16*(rk-15.D0)*(rk-16.D0)*(-827.D0+48.D0*rk) * den(17)
        ckplm(k,2413) = 1.88886204389699D-15*(rk-15.D0)*(rk-16.D0)*(rk-17.D0) * den(18)
        ckplm(k,2414) = 5.55547659969704D-17*(rk+17.D0)*(rk+18.D0) * den(0)
        ckplm(k,2415) = -3.17455805696973D-18*(rk+17.D0)*(4480.D0+247.D0*rk) * den(1)
        ckplm(k,2416) = 2.69837434842428D-17*(60288.D0+rk*(6725.D0+187.D0*rk)) * den(2)
        ckplm(k,2417) = -2.15869947873942D-15*(3202.D0+rk*(341.D0+9.D0*rk)) * den(3)
        ckplm(k,2418) = 3.77772408779399D-15*(5534.D0+rk*(543.D0+13.D0*rk)) * den(4)
        ckplm(k,2419) = -7.55544817558797D-15*(6372.D0+rk*(547.D0+11.D0*rk)) * den(5)
        ckplm(k,2420) = 9.82208262826436D-15*(8980.D0+rk*(623.D0+9.D0*rk)) * den(6)
        ckplm(k,2421) = -3.92883305130574D-14*(3374.D0+(rk+165.D0)*rk) * den(7)
        ckplm(k,2422) = -3.85867531824671D-14*(-4326.D0+(rk-109.D0)*rk) * den(8)
        ckplm(k,2423) = 7.71735063649343D-14*(rk*rk+rk-2312.D0) * den(9)
        ckplm(k,2424) = -3.85867531824671D-14*(-4216.D0+(rk+111.D0)*rk) * den(10)
        ckplm(k,2425) = -3.92883305130574D-14*(3210.D0+(rk-163.D0)*rk) * den(11)
        ckplm(k,2426) = 9.82208262826436D-15*(8366.D0+rk*(-605.D0+9.D0*rk)) * den(12)
        ckplm(k,2427) = -7.55544817558797D-15*(5836.D0+rk*(-525.D0+11.D0*rk)) * den(13)
        ckplm(k,2428) = 3.77772408779399D-15*(5004.D0+rk*(-517.D0+13.D0*rk)) * den(14)
        ckplm(k,2429) = -2.15869947873942D-15*(2870.D0+rk*(-323.D0+9.D0*rk)) * den(15)
        ckplm(k,2430) = 2.69837434842428D-17*(53750.D0+rk*(-6351.D0+187.D0*rk)) * den(16)
        ckplm(k,2431) = -3.17455805696973D-18*(rk-16.D0)*(-4233.D0+247.D0*rk) * den(17)
        ckplm(k,2432) = 5.55547659969704D-17*(rk-16.D0)*(rk-17.D0) * den(18)
        ckplm(k,2433) = -1.58727902848487D-18*(rk+18.D0) * den(0)
        ckplm(k,2434) = 1.58727902848487D-18*(289.D0+16.D0*rk) * den(1)
        ckplm(k,2435) = -2.69837434842428D-17*(129.D0+7.D0*rk) * den(2)
        ckplm(k,2436) = 2.15869947873942D-16*(77.D0+4.D0*rk) * den(3)
        ckplm(k,2437) = -5.39674869684855D-16*(104.D0+5.D0*rk) * den(4)
        ckplm(k,2438) = 7.55544817558797D-16*(189.D0+8.D0*rk) * den(5)
        ckplm(k,2439) = -9.82208262826436D-15*(rk+29.D0) * den(6)
        ckplm(k,2440) = 2.80630932236125D-15*(163.D0+4.D0*rk) * den(7)
        ckplm(k,2441) = -7.71735063649343D-15*(rk+78.D0) * den(8)
        ckplm(k,2442) = 6.55974804101941D-13 * den(9)
        ckplm(k,2443) = 7.71735063649343D-15*(rk-77.D0) * den(10)
        ckplm(k,2444) = -2.80630932236125D-15*(-159.D0+4.D0*rk) * den(11)
        ckplm(k,2445) = 9.82208262826436D-15*(rk-28.D0) * den(12)
        ckplm(k,2446) = -7.55544817558797D-16*(-181.D0+8.D0*rk) * den(13)
        ckplm(k,2447) = 5.39674869684855D-16*(-99.D0+5.D0*rk) * den(14)
        ckplm(k,2448) = -2.15869947873942D-16*(-73.D0+4.D0*rk) * den(15)
        ckplm(k,2449) = 2.69837434842428D-17*(-122.D0+7.D0*rk) * den(16)
        ckplm(k,2450) = -1.58727902848487D-18*(-273.D0+16.D0*rk) * den(17)
        ckplm(k,2451) = 1.58727902848487D-18*(rk-17.D0) * den(18)
        ckplm(k,2452) = 4.40910841245797D-20 * den(0)
        ckplm(k,2453) = -7.93639514242434D-19 * den(1)
        ckplm(k,2454) = 6.74593587106069D-18 * den(2)
        ckplm(k,2455) = -3.5978324645657D-17 * den(3)
        ckplm(k,2456) = 1.34918717421214D-16 * den(4)
        ckplm(k,2457) = -3.77772408779399D-16 * den(5)
        ckplm(k,2458) = 8.18506885688697D-16 * den(6)
        ckplm(k,2459) = -1.40315466118062D-15 * den(7)
        ckplm(k,2460) = 1.92933765912336D-15 * den(8)
        ckplm(k,2461) = -2.14370851013706D-15 * den(9)
        ckplm(k,2462) = 1.92933765912336D-15 * den(10)
        ckplm(k,2463) = -1.40315466118062D-15 * den(11)
        ckplm(k,2464) = 8.18506885688697D-16 * den(12)
        ckplm(k,2465) = -3.77772408779399D-16 * den(13)
        ckplm(k,2466) = 1.34918717421214D-16 * den(14)
        ckplm(k,2467) = -3.5978324645657D-17 * den(15)
        ckplm(k,2468) = 6.74593587106069D-18 * den(16)
        ckplm(k,2469) = -7.93639514242434D-19 * den(17)
        ckplm(k,2470) = 4.40910841245797D-20 * den(18)
!    ckplm para l = 19
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-3.D0)*(r2k+3.D0)&
                *(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(17) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(18) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(19) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,2471) = 2.62921388282776D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2472) = 1.35013685874939D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,2473) = 1.0415341481781D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,2474) = 894246.490859985D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,2475) = 807706.507873535D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,2476) = 752002.610778809D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,2477) = 714866.679382324D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(6)
        ckplm(k,2478) = 690356.964660645D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(7)
        ckplm(k,2479) = 675349.204559326D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(8)
        ckplm(k,2480) = 668202.652130127D0*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(9)
        ckplm(k,2481) = 668202.652130127D0*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(10)
        ckplm(k,2482) = 675349.204559326D0*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*rk * den(11)
        ckplm(k,2483) = 690356.964660645D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,2484) = 714866.679382324D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,2485) = 752002.610778809D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,2486) = 807706.507873535D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(15)
        ckplm(k,2487) = 894246.490859985D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(16)
        ckplm(k,2488) = 1.0415341481781D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(17)
        ckplm(k,2489) = 1.35013685874939D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(18)
        ckplm(k,2490) = 2.62921388282776D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(19)
        ckplm(k,2491) = -262921.388282776D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2492) = -7105.98346710205D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-37.D0+17.D0*rk) * den(1)
        ckplm(k,2493) = -27408.7933731079D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-14.D0+3.D0*rk) * den(2)
        ckplm(k,2494) = -4706.56047821045D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-99.D0+13.D0*rk) * den(3)
        ckplm(k,2495) = -4251.08688354492D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-124.D0+11.D0*rk) * den(4)
        ckplm(k,2496) = -3957.9084777832D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-145.D0+9.D0*rk) * den(5)
        ckplm(k,2497) = -3762.45620727539D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-162.D0+7.D0*rk) * den(6)
        ckplm(k,2498) = -18167.2885437012D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-35.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(7)
        ckplm(k,2499) = -3554.46949768066D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-184.D0+3.D0*rk) * den(8)
        ckplm(k,2500) = -3516.85606384277D0*(rk+10.D0)*(rk-189.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0) * den(9)
        ckplm(k,2501) = 3516.85606384277D0*(rk+190.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(10)
        ckplm(k,2502) = 3554.46949768066D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(187.D0+3.D0*rk) * den(11)
        ckplm(k,2503) = 18167.2885437012D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk+36.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0) * den(12)
        ckplm(k,2504) = 3762.45620727539D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(169.D0+7.D0*rk) * den(13)
        ckplm(k,2505) = 3957.9084777832D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(154.D0+9.D0*rk) * den(14)
        ckplm(k,2506) = 4251.08688354492D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(135.D0+11.D0*rk) * den(15)
        ckplm(k,2507) = 4706.56047821045D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(112.D0+13.D0*rk) * den(16)
        ckplm(k,2508) = 27408.7933731079D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(17.D0+3.D0*rk) * den(17)
        ckplm(k,2509) = 7105.98346710205D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(54.D0+17.D0*rk) * den(18)
        ckplm(k,2510) = 262921.388282776D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2511) = 12520.0661087036D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2512) = 338.380165100098D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-148.D0+11.D0*rk) * den(1)
        ckplm(k,2513) = 1015.14049530029D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(138.D0+(rk-65.D0)*rk) * den(2)
        ckplm(k,2514) = -522.951164245606D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-462.D0+(rk+131.D0)*rk) * den(3)
        ckplm(k,2515) = -67.4775695800781D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-5084.D0+rk*(969.D0+23.D0*rk)) * den(4)
        ckplm(k,2516) = -20.9413146972656D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-20880.D0+rk*(2791.D0+109.D0*rk)) * den(5)
        ckplm(k,2517) = -59.7215270996094D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-8694.D0+rk*(817.D0+47.D0*rk)) * den(6)
        ckplm(k,2518) = -288.369659423828D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-2030.D0+rk*(129.D0+11.D0*rk)) * den(7)
        ckplm(k,2519) = -18.8067169189453D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-33672.D0+rk*(1291.D0+181.D0*rk)) * den(8)
        ckplm(k,2520) = -3516.85606384277D0*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-188.D0+(rk+3.D0)*rk) * den(9)
        ckplm(k,2521) = -3516.85606384277D0*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(-190.D0+(rk-1.D0)*rk) * den(10)
        ckplm(k,2522) = -18.8067169189453D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(-34782.D0+rk*(-929.D0+181.D0*rk)) * den(11)
        ckplm(k,2523) = -288.369659423828D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-2148.D0+rk*(-107.D0+11.D0*rk)) * den(12)
        ckplm(k,2524) = -59.7215270996094D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-9464.D0+rk*(-723.D0+47.D0*rk)) * den(13)
        ckplm(k,2525) = -20.9413146972656D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-23562.D0+rk*(-2573.D0+109.D0*rk)) * den(14)
        ckplm(k,2526) = -67.4775695800781D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-6030.D0+rk*(-923.D0+23.D0*rk)) * den(15)
        ckplm(k,2527) = -522.951164245606D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-592.D0+(rk-129.D0)*rk) * den(16)
        ckplm(k,2528) = 1015.14049530029D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(204.D0+(rk+67.D0)*rk) * den(17)
        ckplm(k,2529) = 338.380165100098D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(159.D0+11.D0*rk) * den(18)
        ckplm(k,2530) = 12520.0661087036D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2531) = -569.093914031982D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2532) = -15.380916595459D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-333.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,2533) = 138.428249359131D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-186.D0+(rk+37.D0)*rk) * den(2)
        ckplm(k,2534) = 15.380916595459D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(8142.D0+rk&
                *(-3457.D0+rk*(222.D0+13.D0*rk))) * den(3)
        ckplm(k,2535) = 1.98463439941406D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(112344.D0+rk&
                *(-32146.D0+rk*(789.D0+109.D0*rk))) * den(4)
        ckplm(k,2536) = 1.84776306152344D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(180090.D0+rk&
                *(-36173.D0+rk*(-188.D0+111.D0*rk))) * den(5)
        ckplm(k,2537) = 0.159683227539063D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(2.767716D6+rk&
                *(-390976.D0+rk*(-13179.D0+1099.D0*rk))) * den(6)
        ckplm(k,2538) = 0.771041870117188D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(698670.D0+rk&
                *(-66749.D0+rk*(-4614.D0+173.D0*rk))) * den(7)
        ckplm(k,2539) = 0.150856018066406D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(4.063824D6+rk&
                *(-234234.D0+rk*(-30541.D0+551.D0*rk))) * den(8)
        ckplm(k,2540) = 9.40335845947266D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(69938.D0+3.D0&
                *rk*(-559.D0+(rk-184.D0)*rk)) * den(9)
        ckplm(k,2541) = -9.40335845947266D0*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(-71060.D0+3.D0*(rk+188.D0)*(rk-1.D0)*rk) * den(10)
        ckplm(k,2542) = -0.150856018066406D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(-4.266966D6+rk*(-171499.D0+rk*(32194.D0+551.D0*rk))) * den(11)
        ckplm(k,2543) = -0.771041870117188D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-760632.D0+rk&
                *(-57002.D0+rk*(5133.D0+173.D0*rk))) * den(12)
        ckplm(k,2544) = -0.159683227539063D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-3.144414D6+rk*(-361321.D0+rk*(16476.D0+1099.D0*rk))) * den(13)
        ckplm(k,2545) = -1.84776306152344D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-215964.D0+rk&
                *(-35464.D0+rk*(521.D0+111.D0*rk))) * den(14)
        ckplm(k,2546) = -1.98463439941406D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-145170.D0+rk&
                *(-33397.D0+rk*(-462.D0+109.D0*rk))) * den(15)
        ckplm(k,2547) = -15.380916595459D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-11808.D0+rk&
                *(-3862.D0+rk*(-183.D0+13.D0*rk))) * den(16)
        ckplm(k,2548) = -138.428249359131D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-222.D0+(rk-35.D0)*rk) * den(17)
        ckplm(k,2549) = 15.380916595459D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk+334.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(18)
        ckplm(k,2550) = 569.093914031982D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2551) = 24.7432136535645D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,2552) = -0.668735504150391D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(592.D0+13.D0*rk) * den(1)
        ckplm(k,2553) = -2.00620651245117D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1608.D0+rk*(121.D0+7.D0*rk)) * den(2)
        ckplm(k,2554) = -0.668735504150391D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(32496.D0+rk&
                *(-7426.D0+rk*(-9.D0+19.D0*rk))) * den(3)
        ckplm(k,2555) = -0.0862884521484375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.679766D6+rk*(639730.D0+rk&
                *(-46505.D0+rk*(-2590.D0+101.D0*rk)))) * den(4)
        ckplm(k,2556) = -0.0267791748046875D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-9.45081D6+rk*(2.532582D6+rk&
                *(-72667.D0+rk*(-12458.D0+143.D0*rk)))) * den(5)
        ckplm(k,2557) = 0.0069427490234375D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(5.4151902D7+rk*(-1.0216974D7+rk&
                *(-100661.D0+rk*(52746.D0+137.D0*rk)))) * den(6)
        ckplm(k,2558) = 0.0335235595703125D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.478211D7+rk*(-1.887254D6+rk&
                *(-100001.D0+rk*(9746.D0+149.D0*rk)))) * den(7)
        ckplm(k,2559) = 0.0502853393554688D0*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.1800476D7+rk*(-909198.D0+rk&
                *(-109697.D0+rk*(4562.D0+157.D0*rk)))) * den(8)
        ckplm(k,2560) = 9.40335845947266D0*(rk+10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(69564.D0+rk*(-2230.D0+rk&
                *(-725.D0+(rk+10.D0)*rk))) * den(9)
        ckplm(k,2561) = 9.40335845947266D0*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(71060.D0+(rk-1.D0)*rk&
                *(-754.D0+(rk-5.D0)*rk)) * den(10)
        ckplm(k,2562) = 0.0502853393554688D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(1.2595572D7+rk*(676746.D0+rk&
                *(-122441.D0+rk*(-3934.D0+157.D0*rk)))) * den(11)
        ckplm(k,2563) = 0.0335235595703125D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(1.6559766D7+rk*(1.65861D6+rk&
                *(-128345.D0+rk*(-9150.D0+149.D0*rk)))) * den(12)
        ckplm(k,2564) = 0.0069427490234375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(6.4215606D7+rk*(9.857962D6+rk&
                *(-258077.D0+rk*(-52198.D0+137.D0*rk)))) * den(13)
        ckplm(k,2565) = -0.0267791748046875D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1.2043458D7+rk*(-2.63997D6+rk&
                *(-34435.D0+rk*(13030.D0+143.D0*rk)))) * den(14)
        ckplm(k,2566) = -0.0862884521484375D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2.36331D6+rk*(-724566.D0+rk&
                *(-38129.D0+rk*(2994.D0+101.D0*rk)))) * den(15)
        ckplm(k,2567) = -0.668735504150391D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-39894.D0+rk&
                *(-7351.D0+rk*(66.D0+19.D0*rk))) * den(16)
        ckplm(k,2568) = -2.00620651245117D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1722.D0+rk&
                *(-107.D0+7.D0*rk)) * den(17)
        ckplm(k,2569) = -0.668735504150391D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-579.D0+13.D0*rk) * den(18)
        ckplm(k,2570) = 24.7432136535645D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(19)
        ckplm(k,2571) = -1.03096723556519D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,2572) = 0.0278639793395996D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(925.D0+31.D0*rk) * den(1)
        ckplm(k,2573) = 0.0167183876037598D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-18750.D0+rk&
                *(203.D0+47.D0*rk)) * den(2)
        ckplm(k,2574) = 0.0278639793395996D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(102390.D0+rk*(-11209.D0+rk&
                *(-594.D0+13.D0*rk))) * den(3)
        ckplm(k,2575) = -0.00359535217285156D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(6.56394D6+rk*(-1.468894D6+rk*(7649.D0+rk&
                *(6706.D0+19.D0*rk)))) * den(4)
        ckplm(k,2576) = -0.00334739685058594D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-5.74983D7+rk*(1.9242182D7+rk*(-1.205525D6+rk&
                *(-100155.D0+rk*(5685.D0+113.D0*rk))))) * den(5)
        ckplm(k,2577) = -0.0002892812093099D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.105083D9+rk*(2.60869992D8+rk*(-3.68948D6+rk&
                *(-1.760665D6+rk*(28040.D0+1813.D0*rk))))) * den(6)
        ckplm(k,2578) = -0.00027936299641927D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.6311813D9+rk*(2.60820002D8+rk*(9.578475D6+rk&
                *(-1.929445D6+rk*(-17355.D0+1823.D0*rk))))) * den(7)
        ckplm(k,2579) = -0.00209522247314453D0*(rk+10.D0)*(rk+11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.7408384D8+rk*(2.6461152D7+rk*(2.93275D6+rk&
                *(-200345.D0+rk*(-7690.D0+173.D0*rk))))) * den(8)
        ckplm(k,2580) = -0.130602200826009D0*(rk+10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-4.98168D6+rk*(200154.D0+rk*(64255.D0+rk&
                *(-1435.D0+(rk-175.D0)*rk)))) * den(9)
        ckplm(k,2581) = 0.130602200826009D0*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(5.11632D6+(rk-1.D0)*rk*(-68044.D0+rk&
                *(-544.D0+(rk+181.D0)*rk))) * den(10)
        ckplm(k,2582) = 0.00209522247314453D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(2.9741976D8+rk*(2.0026242D7+rk*(-3.485915D6+rk&
                *(-167855.D0+rk*(8555.D0+173.D0*rk))))) * den(11)
        ckplm(k,2583) = 0.00027936299641927D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(1.88051256D9+rk*(2.35953252D8+rk*(-1.524445D7+rk&
                *(-1.841795D6+rk*(26470.D0+1823.D0*rk))))) * den(12)
        ckplm(k,2584) = 0.0002892812093099D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.36785558D9+rk*(2.62863862D8+rk*(-1.742625D6+rk&
                *(-1.854695D6+rk*(-18975.D0+1813.D0*rk))))) * den(13)
        ckplm(k,2585) = 0.00334739685058594D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(7.784028D7+rk*(2.1330592D7+rk*(872080.D0+rk&
                *(-121765.D0+rk*(-5120.D0+113.D0*rk))))) * den(14)
        ckplm(k,2586) = 0.00359535217285156D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(8.033796D6+rk*(1.46415D6+rk*(-12355.D0+rk&
                *(-6630.D0+19.D0*rk)))) * den(15)
        ckplm(k,2587) = -0.0278639793395996D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-112992.D0+rk*(-9982.D0+rk&
                *(633.D0+13.D0*rk))) * den(16)
        ckplm(k,2588) = -0.0167183876037598D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-18906.D0+rk&
                *(-109.D0+47.D0*rk)) * den(17)
        ckplm(k,2589) = -0.0278639793395996D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-894.D0+31.D0*rk) * den(18)
        ckplm(k,2590) = 1.03096723556519D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) &
                * den(19)
        ckplm(k,2591) = 0.0412386894226074D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2592) = -0.00111455917358398D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1332.D0+53.D0*rk) &
                * den(1)
        ckplm(k,2593) = -0.00143300465175084D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-17682.D0+rk*(-467.D0+19.D0&
                *rk)) * den(2)
        ckplm(k,2594) = 0.00015922273908343D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.905324D6+rk*(62164.D0+rk&
                *(10011.D0+59.D0*rk))) * den(3)
        ckplm(k,2595) = 0.00014381408691406D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(2.1446172D7+rk*(-2.44935D6+rk*(-146195.D0+rk&
                *(7710.D0+203.D0*rk)))) * den(4)
        ckplm(k,2596) = 0.00013389587402344D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.1873888D8+rk*(4.3879548D7+rk*(319250.D0+rk*(-279155.D0+rk&
                *(-530.D0+227.D0*rk))))) * den(5)
        ckplm(k,2597) = 1.65303548177083D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(1.6439282748D11+rk*(-4.6576666488D10+rk*(1.780170928D9+rk&
                *(3.51855405D8+rk*(-1.5558275D7+rk*(-686697.D0+11047.D0*rk)))))) * den(6)
        ckplm(k,2598) = 0.00027936299641927D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(1.49951844D9+rk*(-2.88174504D8+rk*(-5.864936D6+rk*(2.675655D6+rk&
                *(-8285.D0+(rk-5451.D0)*rk))))) * den(7)
        ckplm(k,2599) = -0.00125713348388672D0*(rk+10.D0)*(rk+11.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-4.4201376D8+rk*(5.1327504D7+rk*(5.180442D6+rk*(-514045.D0+rk&
                *(-18295.D0+rk*(1021.D0+13.D0*rk)))))) * den(8)
        ckplm(k,2600) = -0.0261204401652018D0*(rk+10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-2.477376D7+rk*(1.197624D6+rk*(379534.D0+rk*(-11805.D0+rk&
                *(-1535.D0+(rk+21.D0)*rk))))) * den(9)
        ckplm(k,2601) = -0.0261204401652018D0*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(-2.55816D7+(rk-1.D0)*rk*(409380.D0+rk*(3836.D0+rk*(-1639.D0+(rk-14.D0)&
                *rk)))) * den(10)
        ckplm(k,2602) = -0.00125713348388672D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(-4.8766608D8+rk*(-3.9502692D7+rk*(6.602792D6+rk*(430915.D0+rk&
                *(-23205.D0+rk*(-943.D0+13.D0*rk)))))) * den(11)
        ckplm(k,2603) = 0.00027936299641927D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(1.77914952D9+rk*(2.68411788D8+rk*(-1.3887086D7+rk*(-2.654265D6+rk&
                *(18985.D0+(rk+5457.D0)*rk))))) * den(12)
        ckplm(k,2604) = 1.65303548177083D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(2.1238294896D11+rk*(4.9022708796D10+rk*(6.38287738D8+rk*(-4.07000595D8+rk&
                *(-1.1959085D7+rk*(752979.D0+11047.D0*rk)))))) * den(13)
        ckplm(k,2605) = 0.00013389587402344D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.6202078D8+rk*(4.2406838D7+rk*(-1.151265D6+rk*(-274765.D0+rk&
                *(1665.D0+227.D0*rk))))) * den(14)
        ckplm(k,2606) = 0.00014381408691406D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.374182D7+rk*(2.134642D6+rk*(-168107.D0+rk&
                *(-6898.D0+203.D0*rk)))) * den(15)
        ckplm(k,2607) = 0.00015922273908343D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.957536D6+rk*(42319.D0+rk*(-9834.D0+59.D0&
                *rk))) * den(16)
        ckplm(k,2608) = -0.00143300465175084D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-17196.D0+rk*(505.D0+19.D0&
                *rk)) * den(17)
        ckplm(k,2609) = -0.00111455917358398D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1279.D0+53.D0*rk) &
                * den(18)
        ckplm(k,2610) = 0.0412386894226074D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2611) = -0.00158610343933106D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2612) = 0.00004286766052246D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+8.D0)*(rk+9.D0)*(1813.D0+79.D0*rk) * den(1)
        ckplm(k,2613) = 0.00001837185450963D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+8.D0)*(rk+9.D0)*(-96726.D0+rk*(-4807.D0+5.D0*rk)) * den(2)
        ckplm(k,2614) = -6.12395150320871D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0)*(-4.456158D6+rk*(-85997.D0+rk*(14262.D0+293.D0&
                *rk))) * den(3)
        ckplm(k,2615) = -5.53131103515625D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(6.1014324D7+rk*(-2.127574D6+rk*(-449719.D0+rk*(-974.D0+307.D0&
                *rk)))) * den(4)
        ckplm(k,2616) = -5.14984130859375D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-7.2192078D8+rk*(7.1371038D7+rk*(6.510435D6+rk*(-347795.D0+rk&
                *(-16195.D0+97.D0*rk))))) * den(5)
        ckplm(k,2617) = 4.45048014322917D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-8.649941076D10+rk*(1.4158113816D10+rk*(5.78900108D8+rk*(-1.19625165D8+rk&
                *(-1.718395D6+rk*(204369.D0+1787.D0*rk)))))) * den(6)
        ckplm(k,2618) = 5.78562418619792D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+8.D0)*(rk+9.D0)&
                *(6.654308076D10+rk*(-1.4937639696D10+rk*(-5.8060388D7+rk*(1.61537611D8+rk*(-2.9843D6+rk&
                *(-490574.D0+rk*(5608.D0+259.D0*rk))))))) * den(7)
        ckplm(k,2619) = 8.67843627929688D-6*(rk+10.D0)*(rk+11.D0)*(rk-7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(6.194582784D10+rk*(-8.410666176D9+rk*(-7.64909052D8+rk*(1.03820168D8+rk*(3.115525D6+rk&
                *(-328869.D0+rk*(-4393.D0+157.D0*rk))))))) * den(8)
        ckplm(k,2620) = 0.00007727940877279D0*(rk+10.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)&
                *(8.32802256D9+7.D0*rk*(-6.7278852D7+rk*(-2.1041148D7+rk*(844099.D0+rk*(112350.D0+rk&
                *(-2528.D0+(rk-162.D0)*rk)))))) * den(9)
        ckplm(k,2621) = -0.00007727940877279D0*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(-8.6465808D9+7.D0*(rk-1.D0)*rk*(2.312532D7+rk*(248804.D0+rk*(-123890.D0+rk&
                *(-1365.D0+(rk+170.D0)*rk))))) * den(10)
        ckplm(k,2622) = -8.67843627929688D-6*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(-6.949120464D10+rk*(-6.583466556D9+rk*(1.054456908D9+rk*(8.8162733D7+rk*(-4.68848D6+rk&
                *(-299214.D0+rk*(5492.D0+157.D0*rk))))))) * den(11)
        ckplm(k,2623) = -5.78562418619792D-6*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-8.125863408D10+rk*(-1.4327453592D10+rk*(5.555946D8+rk*(1.68465976D8+rk*(456375.D0+rk&
                *(-518783.D0+rk*(-3795.D0+259.D0*rk))))))) * den(12)
        ckplm(k,2624) = -4.45048014322917D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-9.996092028D10+rk*(-1.2649322808D10+rk*(9.25448348D8+rk*(1.10743635D8+rk&
                *(-2.713435D6+rk*(-193647.D0+1787.D0*rk)))))) * den(13)
        ckplm(k,2625) = 5.14984130859375D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(7.8644988D8+rk*(5.7372048D7+rk*(-7.45568D6+rk*(-282045.D0+rk&
                *(16680.D0+97.D0*rk))))) * den(14)
        ckplm(k,2626) = 5.53131103515625D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(6.269346D7+rk*(1.232286D6+rk*(-444955.D0+rk*(2202.D0+307.D0&
                *rk)))) * den(15)
        ckplm(k,2627) = 6.12395150320871D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(4.356192D6+rk*(-113642.D0+rk*(-13383.D0+293.D0&
                *rk))) * den(16)
        ckplm(k,2628) = -0.00001837185450963D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-91914.D0+rk*(4817.D0+5.D0*rk)) * den(17)
        ckplm(k,2629) = -0.00004286766052246D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1734.D0+79.D0*rk) * den(18)
        ckplm(k,2630) = 0.00158610343933106D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2631) = 0.00005874457182708D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+9.D0) * den(0)
        ckplm(k,2632) = -1.58769113046152D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+9.D0)*(2368.D0+109.D0*rk) * den(1)
        ckplm(k,2633) = 2.04131716773624D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+9.D0)*(54432.D0+rk*(3547.D0+37.D0*rk)) * den(2)
        ckplm(k,2634) = 6.80439055912078D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+9.D0)*(-3.144384D6+rk*(-175246.D0+rk*(3141.D0+169.D0*rk))) * den(3)
        ckplm(k,2635) = 6.14590115017361D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+9.D0)*(5.1806952D7+rk*(1.251658D6+rk*(-261437.D0+rk*(-8182.D0+41.D0*rk)))) &
                * den(4)
        ckplm(k,2636) = 6.35782877604167D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0)*(-6.37211664D9+rk*(1.25293344D8+rk*(5.948107D7+rk*(409075.D0-rk*(91690.D0+1039.D0&
                *rk))))) * den(5)
        ckplm(k,2637) = -6.35782877604167D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(-7.352802576D10+rk*(4.962406416D9+rk*(8.92039998D8+rk*(-2.9665075D7+rk*(-3.052625D6+rk&
                *(14179.D0+1427.D0*rk)))))) * den(6)
        ckplm(k,2638) = -8.26517740885417D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)&
                *(6.114160416D10+rk*(-6.964680256D9+rk*(-7.9631086D8+rk*(7.4241968D7+rk*(3.722585D6+rk&
                *(-189289.D0+rk*(-5485.D0+57.D0*rk))))))) * den(7)
        ckplm(k,2639) = 4.13258870442708D-7*(rk+10.D0)*(rk+11.D0)*(rk+9.D0)*(1.25835274656D12+rk&
                *(-1.95666866064D11+rk*(-1.5826024596D10+rk*(2.847569476D9+rk*(6.5173789D7+rk*(-1.2367496D7+rk&
                *(-116294.D0+rk*(14164.D0+61.D0*rk)))))))) * den(8)
        ckplm(k,2640) = 0.00007727940877279D0*(rk+10.D0)*(rk-8.D0)*(rk+9.D0)*(8.28251424D9+rk&
                *(-5.36720976D8+rk*(-1.65602916D8+rk*(8.178004D6+rk*(1.094049D6+rk*(-34984.D0+rk&
                *(-2494.D0+(rk+36.D0)*rk))))))) * den(9)
        ckplm(k,2641) = 0.00007727940877279D0*(rk-8.D0)*(rk-9.D0)*(rk+9.D0)*(8.6465808D9+(rk-1.D0)*rk&
                *(-1.8551704D8+rk*(-2.256108D6+rk*(1.246944D6+rk*(16575.D0+rk*(-2745.D0+(rk-27.D0)*rk)))))) &
                * den(10)
        ckplm(k,2642) = 4.13258870442708D-7*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(1.43542342944D12+rk&
                *(1.55793844656D11+rk*(-2.3856055476D10+rk*(-2.466017564D9+rk*(1.24775389D8+rk*(1.1375704D7+rk&
                *(-213734.D0+rk*(-13676.D0+61.D0*rk)))))))) * den(11)
        ckplm(k,2643) = -8.26517740885417D-7*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-6.723963792D10+rk*(-5.165136108D9+rk*(9.94891836D8+rk*(5.7570433D7+rk*(-4.58476D6+rk&
                *(-155182.D0+rk*(5884.D0+57.D0*rk))))))) * den(12)
        ckplm(k,2644) = 6.35782877604167D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(7.757179248D10+rk*(3.101604028D9+rk*(-9.62599088D8+rk*(-1.7341325D7+rk&
                *(3.102115D6+(5617.D0-1427.D0*rk)*rk))))) * den(13)
        ckplm(k,2645) = 6.35782877604167D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(6.43842864D9+rk*(7.919994D6+rk*(-5.7714095D7+rk*(765445.D0+(86495.D0-1039.D0*rk)&
                *rk)))) * den(14)
        ckplm(k,2646) = 6.14590115017361D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0)*(5.030208D7+rk*(-1.749822D6+rk*(-236645.D0+rk*(8346.D0+41.D0*rk)))) &
                * den(15)
        ckplm(k,2647) = 6.80439055912078D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-8.D0)*(rk-9.D0)*(2.966166D6+rk*(-181021.D0+rk*(-2634.D0+169.D0*rk))) * den(16)
        ckplm(k,2648) = 2.04131716773624D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-8.D0)*(rk-9.D0)*(50922.D0+rk*(-3473.D0+37.D0*rk)) * den(17)
        ckplm(k,2649) = -1.58769113046152D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-8.D0)*(rk-9.D0)*(-2259.D0+109.D0*rk) * den(18)
        ckplm(k,2650) = 0.00005874457182708D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2651) = -2.09802042239558D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2652) = 5.67032546593399D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(2997.D0+143.D0*rk) * den(1)
        ckplm(k,2653) = -1.53098787580218D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(4122.D0+rk*(313.D0+5.D0*rk)) * den(2)
        ckplm(k,2654) = -1.54645239980018D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(-9.659826D6+rk*(-793379.D0+rk*(-9366.D0+251.D0*rk))) * den(3)
        ckplm(k,2655) = 1.99542245135507D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(-1.326233196D9+rk*(-9.157785D7+rk*(2.357375D6+rk*(213150.D0+1861.D0*rk)))) * den(4)
        ckplm(k,2656) = 1.85780710988231D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(2.083454946D10+rk*(8.99201034D8+rk*(-1.28995275D8+rk*(-6.368645D6+rk*(49515.D0+2951.D0&
                *rk))))) * den(5)
        ckplm(k,2657) = 1.58945719401042D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(-3.125886372D10+rk*(-3.53165568D8+rk*(3.40418188D8+rk*(8.543145D6+rk*(-744785.D0+rk&
                *(-20157.D0+97.D0*rk)))))) * den(6)
        ckplm(k,2658) = -2.95184907459077D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(-1.9658918916D11+rk&
                *(4.156820856D9+rk*(2.962205652D9+rk*(-5.684623D6+rk*(-1.261554D7+rk*(-126826.D0+rk&
                *(11928.D0+113.D0*rk))))))) * den(7)
        ckplm(k,2659) = -1.32833208356585D-7*(rk+10.D0)*(rk+11.D0)*(4.7329918272D11+rk&
                *(-2.3794635168D10+rk*(-8.843036652D9+rk*(2.74863372D8+rk*(5.4690013D7+rk*(-728952.D0+rk&
                *(-113078.D0+rk*(108.D0+37.D0*rk)))))))) * den(8)
        ckplm(k,2660) = -7.52721514020647D-7*(rk+10.D0)*(-8.4566594112D11+rk*(6.1815535008D10+rk&
                *(1.881082922D10+rk*(-1.10827788D9+rk*(-1.47525595D8+rk*(6.206739D6+rk*(461930.D0+3.D0*rk&
                *(-3690.D0+(rk-145.D0)*rk)))))))) * den(9)
        ckplm(k,2661) = 7.52721514020647D-7*(rk-9.D0)*(8.877156288D11+(rk-1.D0)*rk&
                *(-2.148733344D10+rk*(-2.91742408D8+rk*(1.7396714D8+rk*(2.693422D6+rk*(-534005.D0+3.D0*rk&
                *(-2339.D0+(rk+155.D0)*rk))))))) * den(10)
        ckplm(k,2662) = 1.32833208356585D-7*(rk-10.D0)*(rk-9.D0)*(4.8803122368D11+rk&
                *(5.505697632D9+rk*(-9.333894572D9+rk*(-5.1077068D7+rk*(5.6637413D7+rk*(50288.D0+rk&
                *(-112798.D0+rk*(188.D0+37.D0*rk)))))))) * den(11)
        ckplm(k,2663) = 2.95184907459077D-8*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(1.9779059664D11+rk&
                *(-1.734887064D9+rk*(-2.905011088D9+rk*(4.3274672D7+rk*(1.1806445D7+rk*(-196021.D0+rk&
                *(-11137.D0+113.D0*rk))))))) * den(12)
        ckplm(k,2664) = 1.58945719401042D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(3.057454764D10-rk*(1.005494736D9+rk*(3.10523068D8+rk*(-1.1318775D7+rk*(-642545.D0+rk&
                *(20739.D0+97.D0*rk)))))) * den(13)
        ckplm(k,2665) = -1.85780710988231D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(-1.981276836D10+rk*(1.137902344D9+rk*(1.0962176D8+rk*(-6.537195D6+rk*(-34760.D0+2951.D0&
                *rk))))) * den(14)
        ckplm(k,2666) = -1.99542245135507D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0)*(-1.23250926D9+rk*(9.5660594D7+rk*(1.729091D6+rk*(-205706.D0+1861.D0*rk)))) * den(15)
        ckplm(k,2667) = 1.54645239980018D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-9.D0)*(8.876064D6+rk*(-773894.D0+rk*(10119.D0+251.D0*rk))) * den(16)
        ckplm(k,2668) = 1.53098787580218D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-9.D0)*(3814.D0+rk*(-303.D0+5.D0*rk)) * den(17)
        ckplm(k,2669) = -5.67032546593399D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-9.D0)*(-2854.D0+143.D0*rk) * den(18)
        ckplm(k,2670) = 2.09802042239558D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-9.D0) * den(19)
        ckplm(k,2671) = 7.2345531806744D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2672) = -1.95528464342551D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(3700.D0+181.D0*rk) * den(1)
        ckplm(k,2673) = 3.51951235816592D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(93750.D0+rk*(7849.D0+151.D0*rk)) * den(2)
        ckplm(k,2674) = -5.33259448206958D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(1.769922D7+rk*(1.801652D6+rk*(46587.D0+91.D0*rk))) * den(3)
        ckplm(k,2675) = 6.88076707363817D-11*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(2.86294362D9+rk*(2.97268202D8+rk*(4.622693D6-rk*(260498.D0+5357.D0*rk)))) * den(4)
        ckplm(k,2676) = -1.03211506104573D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(3.19271328D9+rk*(2.99265564D8+rk*(-4.88965D6+rk*(-994909.D0+rk*(-17294.D0+109.D0*rk))))) &
                * den(5)
        ckplm(k,2677) = 7.94728597005208D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(5.90877D8+rk&
                *(4.535064D7+rk*(-3.54764D6+rk*(-326623.D0+rk*(-367.D0+rk*(307.D0+3.D0*rk)))))) * den(6)
        ckplm(k,2678) = 5.90369814918155D-9*(rk+11.D0)*(rk+12.D0)*(-1.003188732D11+rk&
                *(-5.80073932D9+rk*(1.100699152D9+rk*(7.0378378D7+rk*(-2.394245D6+rk*(-178325.D0+rk&
                *(-287.D0+47.D0*rk))))))) * den(7)
        ckplm(k,2679) = 1.47592453729539D-8*(rk+11.D0)*(4.577017536D11+rk*(1.876855392D10+rk&
                *(-7.331045772D9+rk*(-3.31196852D8+rk*(3.3203569D7+rk*(1.54504D6+rk*(-37478.D0+(rk-1628.D0)&
                *rk))))))) * den(8)
        ckplm(k,2680) = -2.50907171340216D-7*(2.803312512D11+rk*(8.31887424D9+rk*(-5.924967504D9+rk&
                *(-1.87563356D8+rk*(4.1111896D7+rk*(1.296449D6+rk*(-102536.D0+rk*(-2854.D0+(rk+64.D0)&
                *rk)))))))) * den(9)
        ckplm(k,2681) = -2.50907171340216D-7*(-2.6631468864D11+rk*(1.9448748576D10+rk*(5.1300469D9+rk&
                *(-3.3709912D8+rk*(-3.3195855D7+rk*(1.848273D6+rk*(80850.D0+rk*(-3330.D0+(rk-55.D0)*rk)))))))) &
                * den(10)
        ckplm(k,2682) = 1.47592453729539D-8*(rk-10.D0)*(4.3196497344D11+rk*(-3.2312179296D10+rk&
                *(-6.154212156D9+rk*(4.47868204D8+rk*(2.4973249D7+rk*(-1.735664D6+rk*(-26054.D0+(rk+1636.D0)&
                *rk))))))) * den(11)
        ckplm(k,2683) = 5.90369814918155D-9*(rk-10.D0)*(rk-11.D0)*(9.349002936D10+rk&
                *(-7.782315084D9+rk*(-8.76976506D8+rk*(7.8179493D7+rk*(1.50857D6+rk*(-175616.D0+rk&
                *(616.D0+47.D0*rk))))))) * den(12)
        ckplm(k,2684) = 7.94728597005208D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(5.42304672D8+rk&
                *(-5.1469036D7+rk*(-2.572998D6+rk*(322145.D0+rk*(-1857.D0+rk*(-289.D0+3.D0*rk)))))) * den(13)
        ckplm(k,2685) = -1.03211506104573D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(-2.889535572D9+rk*(3.06129858D8+rk*(2.009777D6+rk*(-924643.D0+rk*(17839.D0+109.D0*rk))))) &
                * den(14)
        ckplm(k,2686) = 6.88076707363817D-11*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(2.570553252D9+rk*(-2.8726275D8+rk*(5.372045D6+(239070.D0-5357.D0*rk)*rk))) * den(15)
        ckplm(k,2687) = -5.33259448206958D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(-1.5944064D7+rk*(1.708751D6+rk*(-46314.D0+91.D0*rk))) * den(16)
        ckplm(k,2688) = 3.51951235816592D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(86052.D0+rk*(-7547.D0+151.D0*rk)) * den(17)
        ckplm(k,2689) = -1.95528464342551D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(-3519.D0+181.D0*rk) * den(18)
        ckplm(k,2690) = 7.2345531806744D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0) * den(19)
        ckplm(k,2691) = -2.41151772689147D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2692) = 6.51761547808504D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(4477.D0+223.D0*rk) * den(1)
        ckplm(k,2693) = -1.95528464342551D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(82038.D0+rk*(7345.D0+157.D0*rk)) * den(2)
        ckplm(k,2694) = 2.96255249003866D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(1.8439674D7+rk*(2.153513D6+rk*(73938.D0+643.D0*rk))) * den(3)
        ckplm(k,2695) = 3.82264837424343D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(-3.480192804D9+rk*(-4.56702302D8+rk*(-1.6697615D7+rk*(-40558.D0+4259.D0*rk)))) * den(4)
        ckplm(k,2696) = -6.88076707363817D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(-3.679811586D10+rk&
                *(-4.966082106D9+rk*(-1.33819225D8+rk*(5.971145D6+rk*(276025.D0+1861.D0*rk))))) * den(5)
        ckplm(k,2697) = -5.88687849633488D-11*(rk+12.D0)*(rk+13.D0)*(6.789323772D10+rk&
                *(8.943976392D9+rk*(5.1887848D7+rk*(-3.2813325D7+rk*(-1.127255D6+rk*(3513.D0+307.D0*rk)))))) &
                * den(6)
        ckplm(k,2698) = 1.09327743503362D-10*(rk+12.D0)*(4.9754580876D11+rk*(6.2595122256D10+rk&
                *(-1.580620832D9+rk*(-4.39593329D8+rk*(-9.51692D6+rk*(492814.D0+rk*(15232.D0+19.D0*rk))))))) &
                * den(7)
        ckplm(k,2699) = 1.63991615255043D-10*(-3.98302736832D12+rk*(-4.82472137376D11+rk&
                *(3.0227853876D10+rk*(5.239599484D9+rk*(2.5251541D7+rk*(-1.3706504D7+rk*(-293846.D0+rk&
                *(5836.D0+109.D0*rk)))))))) * den(8)
        ckplm(k,2700) = 1.02221440175643D-8*(6.220573632D10+rk*(1.889993376D9+rk*(-9.22314396D8+rk&
                *(-2.8162004D7+rk*(4.122649D6+rk*(116104.D0+rk*(-5534.D0+(rk-116.D0)*rk))))))) * den(9)
        ckplm(k,2701) = -1.02221440175643D-8*(5.942559168D10+rk*(-3.634258464D9+rk*(-8.14334076D8+rk&
                *(4.3384996D7+rk*(3.463249D6+rk*(-146816.D0+rk*(-4694.D0+(rk+124.D0)*rk))))))) * den(10)
        ckplm(k,2702) = -1.63991615255043D-10*(-3.47552831808D12+rk*(5.27376782304D11+rk&
                *(1.4793102516D10+rk*(-5.007603356D9+rk*(8.9179741D7+rk*(1.1826976D7+rk*(-331646.D0+rk&
                *(-4964.D0+109.D0*rk)))))))) * den(11)
        ckplm(k,2703) = -1.09327743503362D-10*(rk-11.D0)*(-4.3379966448D11+rk*(6.4478024424D10+rk&
                *(3.23642424D8+rk*(-3.96901484D8+rk*(1.1753175D7+rk*(401821.D0+rk*(-15099.D0+19.D0*rk))))))) &
                * den(12)
        ckplm(k,2704) = 5.88687849633488D-11*(rk-11.D0)*(rk-12.D0)*(5.903283204D10+rk&
                *(-8.746285464D9+rk*(1.43533768D8+rk*(2.8275315D7+rk*(-1.140215D6+rk*(-1671.D0+307.D0*rk)))))) &
                * den(13)
        ckplm(k,2705) = 6.88076707363817D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(3.197154996D10+rk&
                *(-4.681625016D9+rk*(1.5009512D8+rk*(4.885655D6+rk*(-266720.D0+1861.D0*rk))))) * den(14)
        ckplm(k,2706) = 3.82264837424343D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(3.0401433D9-rk*(4.23445782D8+rk*(-1.6550387D7+rk*(57594.D0+4259.D0*rk)))) * den(15)
        ckplm(k,2707) = -2.96255249003866D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(-1.6359456D7+rk*(2.007566D6+rk*(-72009.D0+643.D0*rk))) * den(16)
        ckplm(k,2708) = 1.95528464342551D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(74850.D0+rk*(-7031.D0+157.D0*rk)) * den(17)
        ckplm(k,2709) = -6.51761547808504D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(-4254.D0+223.D0*rk) * den(18)
        ckplm(k,2710) = 2.41151772689147D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0) * den(19)
        ckplm(k,2711) = 7.77908944158537D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2712) = -2.10245660583388D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(5328.D0+269.D0*rk) * den(1)
        ckplm(k,2713) = 1.8922109452505D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(38616.D0+rk*(3629.D0+83.D0*rk)) * den(2)
        ckplm(k,2714) = -9.55662093560856D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(3.0603408D7+rk*(3.931298D6+rk*(157497.D0+1885.D0*rk))) * den(3)
        ckplm(k,2715) = 3.82264837424343D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(2.15448678D8+rk&
                *(3.2990958D7+rk*(1.649785D6+rk*(27006.D0+11.D0*rk)))) * den(4)
        ckplm(k,2716) = 2.06423012209145D-11*(rk+13.D0)*(rk+14.D0)*(-8.5951692D8+rk*(-1.45183482D8+rk&
                *(-7.56015D6+rk*(-67735.D0+rk*(4050.D0+67.D0*rk))))) * den(5)
        ckplm(k,2717) = 5.88687849633488D-11*(rk+13.D0)*(5.25724668D9+rk*(9.37624218D8+rk&
                *(4.4028904D7+rk*(-783585.D0+rk*(-100355.D0+(rk-1743.D0)*rk))))) * den(6)
        ckplm(k,2718) = -1.09327743503362D-10*(4.154118696D10+rk*(7.679577636D9+rk*(2.81268918D8+rk&
                *(-2.3151646D7+rk*(-1.741845D6+rk*(-17731.D0+rk*(777.D0+11.D0*rk))))))) * den(7)
        ckplm(k,2719) = -4.91974845765129D-10*(-9.72411804D9+rk*(-1.057979952D9+rk*(4.3277808D7+rk&
                *(7.245469D6+rk*(63840.D0+rk*(-9758.D0+(rk-168.D0)*rk)))))) * den(8)
        ckplm(k,2720) = 9.29285819778577D-10*(-5.15485152D9+rk*(-1.52230056D8+rk*(5.208172D7+rk&
                *(1.455424D6+rk*(-140825.D0+rk*(-3269.D0+(rk+85.D0)*rk)))))) * den(9)
        ckplm(k,2721) = 9.29285819778577D-10*(4.95213264D9+rk*(-2.51480772D8+rk*(-4.6904442D7+rk&
                *(1.984369D6+rk*(123240.D0+rk*(-3758.D0+(rk-78.D0)*rk)))))) * den(10)
        ckplm(k,2722) = -4.91974845765129D-10*(8.63003232D9+rk*(-1.123102296D9+rk*(-2.201948D7+rk&
                *(6.895924D6+rk*(-110075.D0+rk*(-8729.D0+(rk+175.D0)*rk)))))) * den(11)
        ckplm(k,2723) = -1.09327743503362D-10*(-3.416430654D10+rk*(7.054459002D9+rk*(-3.4046152D8+rk&
                *(-1.6376731D7+rk*(1.64192D6+rk*(-22162.D0+rk*(-700.D0+11.D0*rk))))))) * den(12)
        ckplm(k,2724) = 5.88687849633488D-11*(rk-12.D0)*(4.36433634D9+rk*(-8.47608354D8+rk&
                *(4.5794974D7+rk*(399615.D0+rk*(-91625.D0+(rk+1749.D0)*rk))))) * den(13)
        ckplm(k,2725) = 2.06423012209145D-11*(rk-12.D0)*(rk-13.D0)*(7.2182187D8+rk*(-1.30282252D8+rk&
                *(7.333315D6+rk*(-83265.D0+rk*(-3715.D0+67.D0*rk))))) * den(14)
        ckplm(k,2726) = 3.82264837424343D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(1.8408051D8+rk&
                *(-2.9772362D7+rk*(1.568833D6+rk*(-26962.D0+11.D0*rk)))) * den(15)
        ckplm(k,2727) = 9.55662093560856D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(2.6827722D7+rk*(-3.621959D6+(151842.D0-1885.D0*rk)*rk)) * den(16)
        ckplm(k,2728) = 1.8922109452505D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(35070.D0+rk*(-3463.D0+83.D0*rk)) * den(17)
        ckplm(k,2729) = 2.10245660583388D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(5059.D0-269.D0*rk) * den(18)
        ckplm(k,2730) = 7.77908944158537D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0) * den(19)
        ckplm(k,2731) = -2.43096545049543D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0) * den(0)
        ckplm(k,2732) = 6.57017689323089D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(6253.D0+319.D0*rk) * den(1)
        ckplm(k,2733) = -2.81579009709895D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(1.114386D6+rk*(108599.D0+2603.D0*rk)) * den(2)
        ckplm(k,2734) = 4.26634863196811D-15*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(3.42104334D8+rk&
                *(4.7111303D7+rk*(2.082222D6+29077.D0*rk))) * den(3)
        ckplm(k,2735) = -1.19457761695107D-13*(rk+14.D0)*(rk+15.D0)*(3.94575792D8+rk*(6.7418354D7+rk&
                *(4.014401D6+rk*(94162.D0+643.D0*rk)))) * den(4)
        ckplm(k,2736) = -2.15023971051193D-13*(rk+14.D0)*(-5.31826776D9+rk*(-1.048525746D9+rk&
                *(-7.2821975D7+rk*(-1.953305D6+rk*(-8425.D0+251.D0*rk))))) * den(5)
        ckplm(k,2737) = 2.62807075729236D-13*(-8.366497776D10+rk*(-1.8253245876D10+rk&
                *(-1.356514828D9+rk*(-3.0835755D7+rk*(678935.D0+rk*(35151.D0+293.D0*rk)))))) * den(6)
        ckplm(k,2738) = 4.44143957982408D-11*(6.0200784D8+rk*(9.6096924D7+rk*(3.367012D6+rk&
                *(-144591.D0+rk*(-9653.D0+(rk-93.D0)*rk))))) * den(7)
        ckplm(k,2739) = -6.66215936973612D-11*(4.4497152D8+rk*(4.2395436D7+rk*(-1.078316D6+rk&
                *(-171855.D0+rk*(-1685.D0+(rk+99.D0)*rk))))) * den(8)
        ckplm(k,2740) = -5.39317663264353D-11*(-5.6617344D8+rk*(-1.5552852D7+rk*(3.735124D6+rk&
                *(90225.D0+rk*(-5525.D0+(rk-93.D0)*rk))))) * den(9)
        ckplm(k,2741) = 5.39317663264353D-11*(-5.4698112D8+rk*(2.2730796D7+rk*(3.432244D6+rk&
                *(-111375.D0+rk*(-5045.D0+(rk+99.D0)*rk))))) * den(10)
        ckplm(k,2742) = 6.66215936973612D-11*(4.0166784D8+rk*(-4.4043732D7+rk*(-573836.D0+rk&
                *(164145.D0+rk*(-2165.D0+(rk-93.D0)*rk))))) * den(11)
        ckplm(k,2743) = -4.44143957982408D-11*(5.0941296D8+rk*(-8.8967268D7+rk*(3.743812D6+rk&
                *(106929.D0+rk*(-9173.D0+(rk+99.D0)*rk))))) * den(12)
        ckplm(k,2744) = 2.62807075729236D-13*(6.673676688D10-rk*(1.5635265228D10+rk&
                *(-1.260281068D9+rk*(3.3205845D7+rk*(507575.D0+rk*(-33393.D0+293.D0*rk)))))) * den(13)
        ckplm(k,2745) = 2.15023971051193D-13*(rk-13.D0)*(4.34061936D9+rk*(-9.08706756D8+rk&
                *(6.701512D7+rk*(-1.917095D6+rk*(9680.D0+251.D0*rk))))) * den(14)
        ckplm(k,2746) = 1.19457761695107D-13*(rk-13.D0)*(rk-14.D0)*(3.3107832D8+rk*(-5.9669466D7+rk&
                *(3.735773D6+rk*(-91590.D0+643.D0*rk)))) * den(15)
        ckplm(k,2747) = 4.26634863196811D-15*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(2.97046176D8+rk&
                *(-4.303409D7+(1.994991D6-29077.D0*rk)*rk)) * den(16)
        ckplm(k,2748) = 2.81579009709895D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(1.00839D6+rk*(-103393.D0+2603.D0*rk)) * den(17)
        ckplm(k,2749) = 6.57017689323089D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(5934.D0-319.D0*rk) * den(18)
        ckplm(k,2750) = 2.43096545049543D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0) * den(19)
        ckplm(k,2751) = 7.36656197119827D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0) &
                * den(0)
        ckplm(k,2752) = -1.99096269491845D-15*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(7252.D0+373.D0*rk) * den(1)
        ckplm(k,2753) = 8.53269726393622D-16*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(1.496166D6+rk&
                *(149945.D0+3719.D0*rk)) * den(2)
        ckplm(k,2754) = -1.42211621065604D-15*(rk+15.D0)*(rk+16.D0)*(4.8044892D7+rk*(6.974092D6+rk&
                *(329985.D0+5057.D0*rk))) * den(3)
        ckplm(k,2755) = 3.9819253898369D-14*(rk+15.D0)*(6.3070896D7+rk*(1.1693422D7+rk*(780055.D0+rk&
                *(21818.D0+209.D0*rk)))) * den(4)
        ckplm(k,2756) = 7.16746570170642D-14*(-9.5268096D8-rk*(2.10398016D8+rk*(1.735435D7+rk&
                *(642355.D0+rk*(9650.D0+29.D0*rk))))) * den(5)
        ckplm(k,2757) = 5.57469554577166D-14*(1.84952808D9+rk*(3.34747518D8+rk*(2.0335285D7+rk&
                *(401785.D0-rk*(2485.D0+103.D0*rk))))) * den(6)
        ckplm(k,2758) = 1.0353006013576D-13*(-1.30332528D9+rk*(-1.74264228D8+rk*(-5.401828D6+rk&
                *(112751.D0+rk*(6148.D0+37.D0*rk))))) * den(7)
        ckplm(k,2759) = 1.55295090203639D-13*(1.01378592D9+rk*(8.1937026D7+rk*(-1.197545D6+rk&
                *(-177695.D0+rk*(-1495.D0+29.D0*rk))))) * den(8)
        ckplm(k,2760) = -4.14859740972579D-12*(4.021248D7+rk*(983144.D0+rk*(-162510.D0+rk&
                *(-3145.D0+(rk+110.D0)*rk)))) * den(9)
        ckplm(k,2761) = -4.14859740972579D-12*(-3.907008D7+rk*(1.298294D6+rk*(152425.D0+rk&
                *(-3575.D0+(rk-105.D0)*rk)))) * den(10)
        ckplm(k,2762) = 1.55295090203639D-13*(-9.3082752D8+rk*(8.3805156D7+rk*(673720.D0+rk&
                *(-171425.D0+rk*(1640.D0+29.D0*rk))))) * den(11)
        ckplm(k,2763) = 1.0353006013576D-13*(1.13456952D9+rk*(-1.63146726D8+rk*(5.703563D6+rk&
                *(88529.D0+rk*(-5963.D0+37.D0*rk))))) * den(12)
        ckplm(k,2764) = 5.57469554577166D-14*(-1.53471168D9+rk*(2.95291728D8+rk*(-1.911605D7+rk&
                *(410695.D0+(1970.D0-103.D0*rk)*rk)))) * den(13)
        ckplm(k,2765) = 7.16746570170642D-14*(7.5900456D8+rk*(-1.77577926D8+rk*(1.5484895D7+rk&
                *(-604045.D0+(9505.D0-29.D0*rk)*rk)))) * den(14)
        ckplm(k,2766) = 3.9819253898369D-14*(rk-14.D0)*(5.213592D7+rk*(-1.019793D7+rk*(715855.D0+rk&
                *(-20982.D0+209.D0*rk)))) * den(15)
        ckplm(k,2767) = 1.42211621065604D-15*(rk-14.D0)*(rk-15.D0)*(4.1395728D7+rk&
                *(-6.329293D6+(314814.D0-5057.D0*rk)*rk)) * den(16)
        ckplm(k,2768) = 8.53269726393622D-16*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(1.34994D6+rk&
                *(-142507.D0+3719.D0*rk)) * den(17)
        ckplm(k,2769) = 1.99096269491845D-15*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(6879.D0-373.D0*rk) * den(18)
        ckplm(k,2770) = 7.36656197119827D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0) &
                * den(19)
        ckplm(k,2771) = -2.16663587388184D-15*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2772) = 5.85577263211309D-17*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(8325.D0+431.D0*rk) &
                * den(1)
        ckplm(k,2773) = -1.50577010540051D-17*(rk+16.D0)*(rk+17.D0)*(3.28125D6+rk*(336193.D0+8557.D0&
                *rk)) * den(2)
        ckplm(k,2774) = 1.42211621065604D-16*(rk+16.D0)*(2.120871D7+rk*(3.207541D6+rk&
                *(159546.D0+2603.D0*rk))) * den(3)
        ckplm(k,2775) = 3.9819253898369D-15*(-3.147048D7-rk*(6.213798D6+rk*(449207.D0+rk&
                *(13998.D0+157.D0*rk)))) * den(4)
        ckplm(k,2776) = 5.37559927627982D-13*(rk+18.D0)*(rk+24.D0)*(1091.D0+(rk+84.D0)*rk) * den(5)
        ckplm(k,2777) = 1.99096269491845D-14*(-2.10708D7+rk*(-3.033262D6+rk*(-140329.D0+(rk-2042.D0)&
                *rk))) * den(6)
        ckplm(k,2778) = -1.0353006013576D-14*(-5.700888D7+rk*(-6.118314D6+rk*(-153527.D0+rk&
                *(1314.D0+47.D0*rk)))) * den(7)
        ckplm(k,2779) = 2.32942635305459D-13*(rk+18.D0)*(-172800.D0+rk*(-1741.D0+(rk+180.D0)*rk)) &
                * den(8)
        ckplm(k,2780) = 1.22017570874288D-13*(6.4736D6+rk*(133682.D0+3.D0*rk*(-4829.D0+(rk-66.D0)&
                *rk))) * den(9)
        ckplm(k,2781) = -1.22017570874288D-13*(6.325632D6+rk*(-162050.D0+3.D0*rk*(-4625.D0+(rk+70.D0)&
                *rk))) * den(10)
        ckplm(k,2782) = -2.32942635305459D-13*(rk-17.D0)*(170880.D0+rk*(-2098.D0+(rk-177.D0)*rk)) &
                * den(11)
        ckplm(k,2783) = 1.0353006013576D-14*(-5.104536D7+rk*(5.807506D6+rk*(-157187.D0+rk&
                *(-1126.D0+47.D0*rk)))) * den(12)
        ckplm(k,2784) = 1.99096269491845D-14*(1.8175824D7-rk*(2.758734D6+rk*(-134197.D0+(rk+2046.D0)&
                *rk))) * den(13)
        ckplm(k,2785) = -5.37559927627982D-13*(rk-17.D0)*(rk-23.D0)*(1008.D0+(rk-82.D0)*rk) * den(14)
        ckplm(k,2786) = 3.9819253898369D-15*(2.5692048D7+rk*(-5.35675D6+rk*(408155.D0+rk&
                *(-13370.D0+157.D0*rk)))) * den(15)
        ckplm(k,2787) = 1.42211621065604D-16*(rk-15.D0)*(1.8158112D7+rk&
                *(-2.896258D6+(151737.D0-2603.D0*rk)*rk)) * den(16)
        ckplm(k,2788) = 1.50577010540051D-17*(rk-15.D0)*(rk-16.D0)*(2.953614D6+rk*(-319079.D0+8557.D0&
                *rk)) * den(17)
        ckplm(k,2789) = 5.85577263211309D-17*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(7894.D0-431.D0*rk) &
                * den(18)
        ckplm(k,2790) = 2.16663587388184D-15*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0) * den(19)
        ckplm(k,2791) = 6.19038821109098D-17*(rk+17.D0)*(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2792) = -1.67307789488946D-18*(rk+17.D0)*(rk+18.D0)*(9472.D0+493.D0*rk) * den(1)
        ckplm(k,2793) = 5.01923368466837D-18*(rk+17.D0)*(363648.D0+rk*(37927.D0+985.D0*rk)) * den(2)
        ckplm(k,2794) = 2.84423242131207D-17*(-4.433664D6-rk*(692806.D0+rk*(35799.D0+611.D0*rk))) &
                * den(3)
        ckplm(k,2795) = 1.13769296852483D-16*(3.232914D6+rk*(468871.D0+rk*(22194.D0+341.D0*rk))) &
                * den(4)
        ckplm(k,2796) = 1.02392367167235D-15*(-806472.D0-rk*(103858.D0+rk*(4227.D0+53.D0*rk))) &
                * den(5)
        ckplm(k,2797) = 7.96385077967381D-16*(1.877796D6+rk*(201359.D0+rk*(6276.D0+49.D0*rk))) &
                * den(6)
        ckplm(k,2798) = 1.47900085908228D-15*(-1.527576D6+rk*(-123182.D0+rk*(-2187.D0+5.D0*rk))) &
                * den(7)
        ckplm(k,2799) = -2.21850128862342D-15*(-1.315014D6+rk*(-65551.D0+rk*(186.D0+19.D0*rk))) &
                * den(8)
        ckplm(k,2800) = 2.44035141748576D-14*(-134096.D0+rk*(-2174.D0+(rk+139.D0)*rk)) * den(9)
        ckplm(k,2801) = 2.44035141748576D-14*(131784.D0+rk*(-2449.D0+(rk-136.D0)*rk)) * den(10)
        ckplm(k,2802) = -2.21850128862342D-15*(1.249296D6+rk*(-65866.D0+rk*(-129.D0+19.D0*rk))) &
                * den(11)
        ckplm(k,2803) = 1.47900085908228D-15*(1.406586D6+rk*(-118793.D0+rk*(2202.D0+5.D0*rk))) &
                * den(12)
        ckplm(k,2804) = 7.96385077967381D-16*(-1.682664D6+rk*(188954.D0+rk*(-6129.D0+49.D0*rk))) &
                * den(13)
        ckplm(k,2805) = 1.02392367167235D-15*(706788.D0+rk*(-95563.D0+(4068.D0-53.D0*rk)*rk)) &
                * den(14)
        ckplm(k,2806) = 1.13769296852483D-16*(-2.785896D6+rk*(425506.D0+rk*(-21171.D0+341.D0*rk))) &
                * den(15)
        ckplm(k,2807) = 2.84423242131207D-17*(3.776046D6+rk*(-623041.D0+(33966.D0-611.D0*rk)*rk)) &
                * den(16)
        ckplm(k,2808) = 5.01923368466837D-18*(rk-16.D0)*(326706.D0+rk*(-35957.D0+985.D0*rk)) * den(17)
        ckplm(k,2809) = 1.67307789488946D-18*(rk-16.D0)*(rk-17.D0)*(8979.D0-493.D0*rk) * den(18)
        ckplm(k,2810) = 6.19038821109098D-17*(rk-16.D0)*(rk-17.D0)*(rk-18.D0) * den(19)
        ckplm(k,2811) = -1.71955228085861D-18*(rk+18.D0)*(rk+19.D0) * den(0)
        ckplm(k,2812) = 4.64743859691515D-20*(rk+18.D0)*(10693.D0+559.D0*rk) * den(1)
        ckplm(k,2813) = 4.18269473722364D-19*(-154326.D0-rk*(16331.D0+431.D0*rk)) * den(2)
        ckplm(k,2814) = 2.37019368442673D-18*(124398.D0+rk*(12643.D0+319.D0*rk)) * den(3)
        ckplm(k,2815) = 9.48077473770691D-18*(-101298.D0-rk*(9595.D0+223.D0*rk)) * den(4)
        ckplm(k,2816) = 2.84423242131207D-17*(83826.D0+rk*(7091.D0+143.D0*rk)) * den(5)
        ckplm(k,2817) = 6.63654231639484D-17*(-70974.D0-rk*(5035.D0+79.D0*rk)) * den(6)
        ckplm(k,2818) = 1.2325007159019D-16*(61926.D0+rk*(3331.D0+31.D0*rk)) * den(7)
        ckplm(k,2819) = 1.84875107385285D-16*(-56058.D0+(rk-1883.D0)*rk) * den(8)
        ckplm(k,2820) = -3.84129389789425D-15*(-3114.D0+(rk-35.D0)*rk) * den(9)
        ckplm(k,2821) = 3.84129389789425D-15*(-3078.D0+(rk+37.D0)*rk) * den(10)
        ckplm(k,2822) = -1.84875107385285D-16*(-54174.D0+(rk+1885.D0)*rk) * den(11)
        ckplm(k,2823) = -1.2325007159019D-16*(58626.D0+rk*(-3269.D0+31.D0*rk)) * den(12)
        ckplm(k,2824) = 6.63654231639484D-17*(66018.D0+rk*(-4877.D0+79.D0*rk)) * den(13)
        ckplm(k,2825) = 2.84423242131207D-17*(-76878.D0+(6805.D0-143.D0*rk)*rk) * den(14)
        ckplm(k,2826) = 9.48077473770691D-18*(91926.D0+rk*(-9149.D0+223.D0*rk)) * den(15)
        ckplm(k,2827) = 2.37019368442673D-18*(-112074.D0+(12005.D0-319.D0*rk)*rk) * den(16)
        ckplm(k,2828) = 4.18269473722364D-19*(138426.D0+rk*(-15469.D0+431.D0*rk)) * den(17)
        ckplm(k,2829) = 4.64743859691515D-20*(rk-17.D0)*(10134.D0-559.D0*rk) * den(18)
        ckplm(k,2830) = 1.71955228085861D-18*(rk-17.D0)*(rk-18.D0) * den(19)
        ckplm(k,2831) = 4.64743859691515D-20*(rk+19.D0) * den(0)
        ckplm(k,2832) = 4.64743859691515D-20*(-324.D0-17.D0*rk) * den(1)
        ckplm(k,2833) = 1.25480842116709D-18*(97.D0+5.D0*rk) * den(2)
        ckplm(k,2834) = 2.37019368442673D-18*(-262.D0-13.D0*rk) * den(3)
        ckplm(k,2835) = 9.48077473770691D-18*(237.D0+11.D0*rk) * den(4)
        ckplm(k,2836) = 2.55980917918087D-16*(-24.D0-rk) * den(5)
        ckplm(k,2837) = 6.63654231639484D-17*(199.D0+7.D0*rk) * den(6)
        ckplm(k,2838) = -1.2325007159019D-16*(186.D0+5.D0*rk) * den(7)
        ckplm(k,2839) = 5.54625322155854D-16*(rk+59.D0) * den(8)
        ckplm(k,2840) = -2.25958464582015D-16*(rk+172.D0) * den(9)
        ckplm(k,2841) = -2.25958464582015D-16*(rk-171.D0) * den(10)
        ckplm(k,2842) = 5.54625322155854D-16*(rk-58.D0) * den(11)
        ckplm(k,2843) = -1.2325007159019D-16*(-181.D0+5.D0*rk) * den(12)
        ckplm(k,2844) = 6.63654231639484D-17*(-192.D0+7.D0*rk) * den(13)
        ckplm(k,2845) = 2.55980917918087D-16*(23.D0-rk) * den(14)
        ckplm(k,2846) = 9.48077473770691D-18*(-226.D0+11.D0*rk) * den(15)
        ckplm(k,2847) = 2.37019368442673D-18*(249.D0-13.D0*rk) * den(16)
        ckplm(k,2848) = 1.25480842116709D-18*(-92.D0+5.D0*rk) * den(17)
        ckplm(k,2849) = 4.64743859691515D-20*(307.D0-17.D0*rk) * den(18)
        ckplm(k,2850) = 4.64743859691515D-20*(rk-18.D0) * den(19)
        ckplm(k,2851) = -1.22301015708294D-21 * den(0)
        ckplm(k,2852) = 2.32371929845758D-20 * den(1)
        ckplm(k,2853) = -2.09134736861182D-19 * den(2)
        ckplm(k,2854) = 1.18509684221336D-18 * den(3)
        ckplm(k,2855) = -4.74038736885346D-18 * den(4)
        ckplm(k,2856) = 1.42211621065604D-17 * den(5)
        ckplm(k,2857) = -3.31827115819742D-17 * den(6)
        ckplm(k,2858) = 6.16250357950949D-17 * den(7)
        ckplm(k,2859) = -9.24375536926424D-17 * den(8)
        ckplm(k,2860) = 1.12979232291007D-16 * den(9)
        ckplm(k,2861) = -1.12979232291007D-16 * den(10)
        ckplm(k,2862) = 9.24375536926424D-17 * den(11)
        ckplm(k,2863) = -6.16250357950949D-17 * den(12)
        ckplm(k,2864) = 3.31827115819742D-17 * den(13)
        ckplm(k,2865) = -1.42211621065604D-17 * den(14)
        ckplm(k,2866) = 4.74038736885346D-18 * den(15)
        ckplm(k,2867) = -1.18509684221336D-18 * den(16)
        ckplm(k,2868) = 2.09134736861182D-19 * den(17)
        ckplm(k,2869) = -2.32371929845758D-20 * den(18)
        ckplm(k,2870) = 1.22301015708294D-21 * den(19)
!    ckplm para l = 20
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+41.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-3.D0)&
                *(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(17) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(18) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(19) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(20) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-39.D0)*(r2k-3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,2871) = 5.38988845979691D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+20.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2872) = 2.76404536399841D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,2873) = 2.12906196956635D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,2874) = 1.8249102596283D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,2875) = 1.64518424921036D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,2876) = 1.52842923797607D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,2877) = 1.44937255325317D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(6)
        ckplm(k,2878) = 1.39569208831787D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(7)
        ckplm(k,2879) = 1.36079978610992D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(8)
        ckplm(k,2880) = 1.34107805007935D6*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(9)
        ckplm(k,2881) = 1.33469196412659D6*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(10)
        ckplm(k,2882) = 1.34107805007935D6*(rk-10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(11)
        ckplm(k,2883) = 1.36079978610992D6*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*rk * den(12)
        ckplm(k,2884) = 1.39569208831787D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,2885) = 1.44937255325317D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,2886) = 1.52842923797607D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(15)
        ckplm(k,2887) = 1.64518424921036D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(16)
        ckplm(k,2888) = 1.8249102596283D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(17)
        ckplm(k,2889) = 2.12906196956635D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(18)
        ckplm(k,2890) = 2.76404536399841D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(19)
        ckplm(k,2891) = 5.38988845979691D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(20)
        ckplm(k,2892) = -513322.710456848D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2893) = -39486.3623428345D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-13.D0+6.D0*rk) * den(1)
        ckplm(k,2894) = -20276.7806625366D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-1.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-37.D0+8.D0*rk) * den(2)
        ckplm(k,2895) = -60830.3419876099D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-15.D0+2.D0*rk) * den(3)
        ckplm(k,2896) = -94010.5285263062D0*(rk+10.D0)*(rk-11.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(4)
        ckplm(k,2897) = -36391.1723327637D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-31.D0+2.D0*rk) * den(5)
        ckplm(k,2898) = -13803.5481262207D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-87.D0+4.D0*rk) * den(6)
        ckplm(k,2899) = -19938.458404541D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-63.D0+2.D0*rk) * den(7)
        ckplm(k,2900) = -25919.9959259033D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-50.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(8)
        ckplm(k,2901) = -6386.08595275879D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-207.D0+2.D0*rk) * den(9)
        ckplm(k,2902) = 1.33469196412659D6*(rk+10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(10)
        ckplm(k,2903) = 6386.08595275879D0*(rk-10.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(209.D0+2.D0*rk) * den(11)
        ckplm(k,2904) = 25919.9959259033D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+51.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0) * den(12)
        ckplm(k,2905) = 19938.458404541D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(65.D0+2.D0*rk) * den(13)
        ckplm(k,2906) = 13803.5481262207D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(91.D0+4.D0*rk) * den(14)
        ckplm(k,2907) = 36391.1723327637D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(33.D0+2.D0*rk) * den(15)
        ckplm(k,2908) = 94010.5285263062D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk+12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,2909) = 60830.3419876099D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(17.D0+2.D0*rk) * den(17)
        ckplm(k,2910) = 20276.7806625366D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(45.D0+8.D0*rk) * den(18)
        ckplm(k,2911) = 39486.3623428345D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(19.D0+6.D0*rk) * den(19)
        ckplm(k,2912) = 513322.710456848D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,2913) = 23332.8504753113D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2914) = 7179.33860778809D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-13.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,2915) = 97.0180892944336D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(2701.D0+rk*(-1281.D0+23.D0*rk)) * den(2)
        ckplm(k,2916) = -582.108535766602D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-780.D0+(rk+224.D0)*rk) * den(3)
        ckplm(k,2917) = -2473.96127700806D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-262.D0+(rk+51.D0)*rk) * den(4)
        ckplm(k,2918) = -3830.64971923828D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-217.D0+(rk+30.D0)*rk) * den(5)
        ckplm(k,2919) = -66.0456848144531D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-15051.D0+rk*(1493.D0+73.D0*rk)) * den(6)
        ckplm(k,2920) = -190.798645019531D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-5922.D0+rk*(412.D0+29.D0*rk)) * den(7)
        ckplm(k,2921) = -62.0095596313477D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-19900.D0+rk*(903.D0+97.D0*rk)) * den(8)
        ckplm(k,2922) = -61.1108703613281D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-21321.D0+rk*(518.D0+103.D0*rk)) * den(9)
        ckplm(k,2923) = -6386.08595275879D0*(rk*rk+rk-209.D0)*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(10)
        ckplm(k,2924) = -61.1108703613281D0*(rk-10.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(-21736.D0+rk*(-312.D0+103.D0*rk)) * den(11)
        ckplm(k,2925) = -62.0095596313477D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(-20706.D0+rk*(-709.D0+97.D0*rk)) * den(12)
        ckplm(k,2926) = -190.798645019531D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-6305.D0+rk*(-354.D0+29.D0*rk)) * den(13)
        ckplm(k,2927) = -66.0456848144531D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-16471.D0+rk*(-1347.D0+73.D0*rk)) * den(14)
        ckplm(k,2928) = -3830.64971923828D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-246.D0+(rk-28.D0)*rk) * den(15)
        ckplm(k,2929) = -2473.96127700806D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-312.D0+(rk-49.D0)*rk) * den(16)
        ckplm(k,2930) = -582.108535766602D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-1003.D0+(rk-222.D0)*rk) * den(17)
        ckplm(k,2931) = 97.0180892944336D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(4005.D0+rk*(1327.D0+23.D0*rk)) * den(18)
        ckplm(k,2932) = 7179.33860778809D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk+14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2933) = 23332.8504753113D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,2934) = -1014.47175979614D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2935) = -26.0120964050293D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-351.D0+2.D0*rk) * den(1)
        ckplm(k,2936) = 2.10908889770508D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-21867.D0+rk*(4447.D0+104.D0*rk)) * den(2)
        ckplm(k,2937) = 2.10908889770508D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(107130.D0+rk*(-46049.D0+rk*(3141.D0+158.D0*rk))) * den(3)
        ckplm(k,2938) = 11.9515037536621D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(33954.D0+rk*(-9919.D0+rk*(294.D0+31.D0*rk))) * den(4)
        ckplm(k,2939) = 13.8791656494141D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(44082.D0+rk*(-9155.D0+rk*(15.D0+26.D0*rk))) * den(5)
        ckplm(k,2940) = 0.478591918945313D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.717902D6+rk*(-256097.D0+rk*(-6093.D0+668.D0*rk))) * den(6)
        ckplm(k,2941) = 0.230433146158854D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(4.408614D6+rk*(-461011.D0+rk*(-24321.D0+1114.D0*rk))) * den(7)
        ckplm(k,2942) = 0.898689270019531D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.3071D6+rk*(-89152.D0+rk*(-8547.D0+199.D0*rk))) * den(8)
        ckplm(k,2943) = 45.8331527709961D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)&
                *(28018.D0+rk*(-1023.D0+rk*(-197.D0+2.D0*rk))) * den(9)
        ckplm(k,2944) = -3193.04297637939D0*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(-418.D0+3.D0*(rk+1.D0)*rk) * den(10)
        ckplm(k,2945) = -45.8331527709961D0*(rk-10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(-28842.D0+rk*(-623.D0+rk*(203.D0+2.D0*rk))) * den(11)
        ckplm(k,2946) = -0.898689270019531D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(-1.387506D6+rk*(-71461.D0+rk*(9144.D0+199.D0*rk))) * den(12)
        ckplm(k,2947) = -0.230433146158854D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-4.84419D6+rk*(-409027.D0+rk*(27663.D0+1114.D0*rk))) * den(13)
        ckplm(k,2948) = -0.478591918945313D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.967238D6+rk*(-241907.D0+rk*(8097.D0+668.D0*rk))) * den(14)
        ckplm(k,2949) = -13.8791656494141D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-53226.D0+rk*(-9107.D0+rk*(63.D0+26.D0*rk))) * den(15)
        ckplm(k,2950) = -11.9515037536621D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-44136.D0+rk*(-10414.D0+rk*(-201.D0+31.D0*rk))) * den(16)
        ckplm(k,2951) = -2.10908889770508D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-156162.D0+rk*(-51857.D0+rk*(-2667.D0+158.D0*rk))) * den(17)
        ckplm(k,2952) = -2.10908889770508D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-26210.D0+rk*(-4239.D0+104.D0*rk)) * den(18)
        ckplm(k,2953) = 26.0120964050293D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(353.D0+2.D0*rk) * den(19)
        ckplm(k,2954) = 1014.47175979614D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,2955) = 42.2696566581726D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2956) = -13.0060482025147D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+52.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,2957) = -0.52727222442627D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-10508.D0+rk*(845.D0+43.D0*rk)) * den(2)
        ckplm(k,2958) = -0.35151481628418D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(107160.D0+rk*(-25198.D0+rk*(117.D0+61.D0*rk))) * den(3)
        ckplm(k,2959) = -0.0878787040710449D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.892696D6+rk&
                *(1.124182D6+rk*(-88487.D0+rk*(-3778.D0+179.D0*rk)))) * den(4)
        ckplm(k,2960) = -0.136070251464844D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-3.307824D6+rk&
                *(916118.D0+rk*(-33553.D0+rk*(-3962.D0+61.D0*rk)))) * den(5)
        ckplm(k,2961) = -0.00234603881835938D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.89684944D8+rk&
                *(5.7655842D7+rk*(-77707.D0+rk*(-267018.D0+307.D0*rk)))) * den(6)
        ckplm(k,2962) = 0.00677744547526042D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.34710128D8+rk&
                *(-1.881813D7+rk*(-660815.D0+rk*(88170.D0+887.D0*rk)))) * den(7)
        ckplm(k,2963) = 0.00220266977945964D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(5.075508D8+rk&
                *(-4.6259806D7+rk*(-3.902689D6+rk*(213394.D0+5101.D0*rk)))) * den(8)
        ckplm(k,2964) = 0.149781545003255D0*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(8.449272D6+rk&
                *(-412314.D0+rk*(-76861.D0+rk*(1806.D0+97.D0*rk)))) * den(9)
        ckplm(k,2965) = 15.6521714528402D0*(rk+10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(85272.D0+(rk&
                *rk+rk-818.D0)*(rk+1.D0)*rk) * den(10)
        ckplm(k,2966) = 0.149781545003255D0*(rk-10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(8.783016D6+rk&
                *(253562.D0+rk*(-81697.D0+rk*(-1418.D0+97.D0*rk)))) * den(11)
        ckplm(k,2967) = 0.00220266977945964D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(5.49699624D8+rk&
                *(3.783465D7+rk*(-4.512265D6+rk*(-192990.D0+5101.D0*rk)))) * den(12)
        ckplm(k,2968) = 0.00677744547526042D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(1.5278016D8+rk&
                *(1.7235538D7+rk*(-920003.D0+rk*(-84622.D0+887.D0*rk)))) * den(13)
        ckplm(k,2969) = -0.00234603881835938D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-3.47151168D8+rk&
                *(-5.7008974D7+rk*(725189.D0+rk*(268246.D0+307.D0*rk)))) * den(14)
        ckplm(k,2970) = -0.136070251464844D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4.253472D6+rk&
                *(-971094.D0+rk*(-21301.D0+rk*(4206.D0+61.D0*rk)))) * den(15)
        ckplm(k,2971) = -0.0878787040710449D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4.101408D6+rk&
                *(-1.289106D6+rk*(-76079.D0+rk*(4494.D0+179.D0*rk)))) * den(16)
        ckplm(k,2972) = -0.35151481628418D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-132414.D0+rk*(-25249.D0+rk*(66.D0+61.D0*rk))) * den(17)
        ckplm(k,2973) = -0.52727222442627D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-11310.D0+rk*(-759.D0+43.D0*rk)) * den(18)
        ckplm(k,2974) = -13.0060482025147D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-4.D0)*(rk-51.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(19)
        ckplm(k,2975) = 42.2696566581726D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,2976) = -1.6907862663269D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,2977) = 0.650302410125732D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(65.D0+2.D0*rk) * den(1)
        ckplm(k,2978) = 0.052727222442627D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-9805.D0+rk*(161.D0+24.D0*rk)) * den(2)
        ckplm(k,2979) = 0.017575740814209D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(270570.D0+rk&
                *(-31709.D0+rk*(-1359.D0+38.D0*rk))) * den(3)
        ckplm(k,2980) = 0.017575740814209D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.27037D6+rk*(529409.D0+rk&
                *(-7399.D0+(rk-2201.D0)*rk))) * den(4)
        ckplm(k,2981) = -0.00136070251464844D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.4338286D8+rk*(8.4134842D7+rk&
                *(-6.026475D6+rk*(-348620.D0+rk*(25095.D0+358.D0*rk))))) * den(5)
        ckplm(k,2982) = -0.00023460388183594D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.39364666D9+rk*(5.95789566D8+rk&
                *(-1.5848125D7+rk*(-3.45641D6+rk*(81565.D0+3284.D0*rk))))) * den(6)
        ckplm(k,2983) = -0.00033887227376302D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2.42042094D9+rk*(4.23278994D8+rk&
                *(7.348295D6+rk*(-2.76962D6+rk*(-1895.D0+2426.D0*rk))))) * den(7)
        ckplm(k,2984) = -0.0110133488972982D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-9.6589128D7+rk*(1.10269D7+rk*(802440.D0+rk&
                *(-75293.D0+rk*(-1788.D0+61.D0*rk))))) * den(8)
        ckplm(k,2985) = -0.187226931254069D0*(rk+10.D0)*(rk+11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-6.661008D6+rk*(407262.D0+rk*(73391.D0+rk&
                *(-2756.D0+rk*(-179.D0+2.D0*rk))))) * den(9)
        ckplm(k,2986) = 7.82608572642009D0*(rk+10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(170544.D0+5.D0*(rk*rk+rk-410.D0)*(rk+1.D0)&
                *rk) * den(10)
        ckplm(k,2987) = 0.187226931254069D0*(rk-10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(6.992304D6+rk*(252938.D0+rk*(-80565.D0+rk&
                *(-2020.D0+rk*(189.D0+2.D0*rk))))) * den(11)
        ckplm(k,2988) = 0.0110133488972982D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(1.06740144D8+rk*(9.203598D6+rk&
                *(-1.016981D6+rk*(-67531.D0+rk*(2093.D0+61.D0*rk))))) * den(12)
        ckplm(k,2989) = 0.00033887227376302D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(2.83358634D9+rk*(4.00293254D8+rk&
                *(-1.5621525D7+rk*(-2.73778D6+rk*(14025.D0+2426.D0*rk))))) * den(13)
        ckplm(k,2990) = 0.00023460388183594D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(3.00174966D9+rk*(6.16806746D8+rk&
                *(5.022345D6+rk*(-3.74983D6+rk*(-65145.D0+3284.D0*rk))))) * den(14)
        ckplm(k,2991) = 0.00136070251464844D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(3.3317082D8+rk*(9.5043342D7+rk&
                *(4.833625D6+rk*(-445420.D0+rk*(-23305.D0+358.D0*rk))))) * den(15)
        ckplm(k,2992) = -0.017575740814209D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2.804976D6+rk*(-537600.D0+rk&
                *(-790.D0+(rk+2205.D0)*rk))) * den(16)
        ckplm(k,2993) = -0.017575740814209D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-300882.D0+rk&
                *(-28877.D0+rk*(1473.D0+38.D0*rk))) * den(17)
        ckplm(k,2994) = -0.052727222442627D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-9942.D0+rk*(-113.D0+24.D0*rk)) * den(18)
        ckplm(k,2995) = -0.650302410125732D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-63.D0+2.D0*rk) * den(19)
        ckplm(k,2996) = 1.6907862663269D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(20)
        ckplm(k,2997) = 0.0650302410125732D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,2998) = -0.0867069880167643D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+27.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(1)
        ckplm(k,2999) = -0.0035151481628418D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-11433.D0+rk&
                *(-235.D0+13.D0*rk)) * den(2)
        ckplm(k,3000) = 0.00703029632568359D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-69240.D0+rk&
                *(2852.D0+(rk+342.D0)*rk)) * den(3)
        ckplm(k,3001) = 0.0005858580271403D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(8.587332D6+rk*(-1.073664D6+rk&
                *(-47341.D0+rk*(3306.D0+67.D0*rk)))) * den(4)
        ckplm(k,3002) = 0.00272140502929688D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.7921844D7+rk*(3.827358D6+rk*(-23635.D0+rk&
                *(-21725.D0+rk*(109.D0+17.D0*rk))))) * den(5)
        ckplm(k,3003) = 0.00004692077636719D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(9.88773966D9+rk*(-2.952200706D9+rk*(1.53213137D8+rk&
                *(1.8020835D7+rk*(-1.06561D6+rk*(-30009.D0+713.D0*rk)))))) * den(6)
        ckplm(k,3004) = 3.47561306423611D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(2.1194552988D11+rk*(-4.4525503968D10+rk*(2.1735302D7+rk&
                *(3.5536545D8+rk*(-5.201545D6+rk*(-649902.D0+2843.D0*rk)))))) * den(7)
        ckplm(k,3005) = -0.00001694361368815D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-5.972733D10+rk*(8.19779568D9+rk*(5.00846258D8+rk&
                *(-7.2763575D7+rk*(-1.302865D6+rk*(132255.D0+887.D0*rk)))))) * den(8)
        ckplm(k,3006) = -0.00115216573079427D0*(rk+10.D0)*(rk+11.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-1.06658136D9+rk*(7.8434934D7+rk*(1.3643591D7+rk*(-715410.D0+rk&
                *(-48220.D0+rk*(1236.D0+29.D0*rk)))))) * den(9)
        ckplm(k,3007) = -0.0401337729560004D0*(rk+10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-3.325608D7+(rk+1.D0)*rk*(480882.D0+(rk*rk+rk-1763.D0)&
                *(rk+1.D0)*rk)) * den(10)
        ckplm(k,3008) = -0.00115216573079427D0*(rk-10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-1.13070672D9+rk*(-4.9200408D7+rk*(1.5488576D7+rk*(510750.D0+rk&
                *(-53965.D0+rk*(-1062.D0+29.D0*rk)))))) * den(11)
        ckplm(k,3009) = -0.00001694361368815D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(-6.735295008D10+rk*(-6.983679852D9+rk*(7.10010548D8+rk&
                *(6.6247305D7+rk*(-1.950835D6+rk*(-126933.D0+887.D0*rk)))))) * den(12)
        ckplm(k,3010) = 3.47561306423611D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(2.561328549D11+rk*(4.348533861D10+rk*(-1.069028653D9+rk&
                *(-3.6961575D8+rk*(-1.90939D6+rk*(666960.D0+2843.D0*rk)))))) * den(13)
        ckplm(k,3011) = 0.00004692077636719D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.297409778D10+rk*(3.200456358D9+rk*(9.3067757D7+rk&
                *(-2.1968925D7+rk*(-904870.D0+rk*(34287.D0+713.D0*rk)))))) * den(14)
        ckplm(k,3012) = 0.00272140502929688D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.175102D7+rk*(3.809102D6+rk*(-42024.D0+rk&
                *(-21991.D0+rk*(-24.D0+17.D0*rk))))) * den(15)
        ckplm(k,3013) = 0.0005858580271403D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(9.610416D6+rk*(969332.D0+rk*(-56857.D0+rk&
                *(-3038.D0+67.D0*rk)))) * den(16)
        ckplm(k,3014) = 0.00703029632568359D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(71751.D0+rk&
                *(2171.D0+(rk-339.D0)*rk)) * den(17)
        ckplm(k,3015) = -0.0035151481628418D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-11185.D0+rk&
                *(261.D0+13.D0*rk)) * den(18)
        ckplm(k,3016) = -0.0867069880167643D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-26.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(19)
        ckplm(k,3017) = 0.0650302410125732D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(20)
        ckplm(k,3018) = -0.00240852744491012D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3019) = 0.00240852744491012D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+8.D0)*(rk+9.D0)*(49.D0+2.D0*rk) &
                * den(1)
        ckplm(k,3020) = 0.00006509533634892D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+8.D0)*(rk+9.D0)*(-41699.D0+rk*(-1821.D0+8.D0&
                *rk)) * den(2)
        ckplm(k,3021) = -0.0000278980012924D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+8.D0)*(rk+9.D0)*(-1.51263D6+rk*(-15271.D0+rk&
                *(4899.D0+82.D0*rk))) * den(3)
        ckplm(k,3022) = -0.00019528600904677D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0)*(2.721474D6+rk*(-127667.D0+rk*(-18463.D0+rk&
                *(83.D0+13.D0*rk)))) * den(4)
        ckplm(k,3023) = -0.00002519819471571D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(-2.37883212D8+rk*(2.7022378D7+rk*(1.733049D6+rk*(-128156.D0+rk&
                *(-4149.D0+46.D0*rk))))) * den(5)
        ckplm(k,3024) = 8.68903266059028D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-7.352531298D10+rk*(1.3325493018D10+rk*(2.46868589D8+rk*(-1.00295925D8+rk&
                *(-258865.D0+rk*(161607.D0+716.D0*rk)))))) * den(6)
        ckplm(k,3025) = 6.0823228624132D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+8.D0)&
                *(rk+9.D0)*(1.0874537172D11+rk*(-2.6670041832D10+rk*(4.90917532D8+rk*(2.39710209D8+rk&
                *(-8.169455D6+rk*(-625503.D0+rk*(13183.D0+306.D0*rk))))))) * den(7)
        ckplm(k,3026) = 0.00003953509860569D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(2.4346998D10+rk*(-3.90533232D9+rk*(-1.92696412D8+rk*(4.1984702D7+rk*(399935.D0+rk&
                *(-119155.D0+rk*(-403.D0+53.D0*rk))))))) * den(8)
        ckplm(k,3027) = 0.00067209667629666D0*(rk+10.D0)*(rk+11.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk+9.D0)*(1.80155664D9+rk*(-1.54915596D8+rk*(-2.5969992D7+rk*(1.772273D6+rk*(117825.D0+rk&
                *(-4999.D0+rk*(-153.D0+2.D0*rk))))))) * den(9)
        ckplm(k,3028) = -0.0200668864780002D0*(rk+10.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(-6.651216D7+7.D0*(rk+1.D0)*rk*(160692.D0+(rk*rk+rk-788.D0)*(rk+1.D0)*rk)) * den(10)
        ckplm(k,3029) = -0.00067209667629666D0*(rk-10.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(-1.92885264D9+rk*(-9.8154156D7+rk*(3.0532208D7+rk*(1.254113D6+rk*(-140455.D0+rk&
                *(-4039.D0+rk*(167.D0+2.D0*rk))))))) * den(11)
        ckplm(k,3030) = -0.00003953509860569D0*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(-2.801816784D10+rk*(-3.396178116D9+rk*(3.15066516D8+rk*(3.9203327D7+rk&
                *(-987810.D0+rk*(-115624.D0+rk*(774.D0+53.D0*rk))))))) * den(12)
        ckplm(k,3031) = -6.0823228624132D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-1.356590898D11+rk*(-2.690327292D10+rk*(2.70783476D8+rk*(2.65880049D8+rk&
                *(4.854905D6+rk*(-698175.D0+rk*(-11041.D0+306.D0*rk))))))) * den(13)
        ckplm(k,3032) = -8.68903266059028D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-8.650406124D10+rk*(-1.2532707264D10+rk*(5.44597844D8+rk*(9.7658715D7+rk&
                *(-1.05616D6+rk*(-157311.D0+716.D0*rk)))))) * den(14)
        ckplm(k,3033) = 0.00002519819471571D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.6304858D8+rk*(2.3188638D7+rk*(-2.092163D6+rk*(-111100.D0+rk&
                *(4379.D0+46.D0*rk))))) * den(15)
        ckplm(k,3034) = 0.00019528600904677D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.830608D6+rk*(90544.D0+rk*(-18634.D0+rk&
                *(-31.D0+13.D0*rk)))) * den(16)
        ckplm(k,3035) = 0.0000278980012924D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.492542D6+rk*(-24823.D0+rk&
                *(-4653.D0+82.D0*rk))) * den(17)
        ckplm(k,3036) = -0.00006509533634892D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-39870.D0+rk*(1837.D0+8.D0&
                *rk)) * den(18)
        ckplm(k,3037) = -0.00240852744491012D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-47.D0+2.D0*rk) &
                * den(19)
        ckplm(k,3038) = 0.00240852744491012D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,3039) = 0.00008601883731822D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+9.D0) * den(0)
        ckplm(k,3040) = -0.00002646733455945D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+9.D0)*(208.D0+9.D0*rk) * den(1)
        ckplm(k,3041) = 3.57666683235836D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+9.D0)*(457616.D0+rk*(27009.D0+223.D0*rk)) &
                * den(2)
        ckplm(k,3042) = 2.14600009941502D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+9.D0)*(-1.48752D6+rk*(-68494.D0+rk*(1941.D0+73.D0*rk))) &
                * den(3)
        ckplm(k,3043) = 5.36500024853754D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+9.D0)*(9.0299808D7+rk*(1.000966D6+rk*(-459291.D0+rk&
                *(-10714.D0+111.D0*rk)))) * den(4)
        ckplm(k,3044) = -2.76903238634196D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+9.D0)*(2.285005536D9+rk*(-8.2296776D7+rk*(-1.9701162D7+rk*(71503.D0+rk&
                *(30942.D0+229.D0*rk))))) * den(5)
        ckplm(k,3045) = -4.7741937695551D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+9.D0)*(-1.57992098832D12+rk*(1.37852157912D11+rk*(1.6122703366D10+rk*(-8.47418355D8+rk&
                *(-5.1795905D7+rk*(702723.D0+24979.D0*rk)))))) * den(6)
        ckplm(k,3046) = -1.24129038008433D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(6.8337301536D11+rk*(-9.3957591216D10+rk*(-6.454234724D9+rk*(9.08492312D8+rk*(2.3932825D7+rk&
                *(-2.211769D6+rk*(-35021.D0+713.D0*rk))))))) * den(7)
        ckplm(k,3047) = -1.55161297510541D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+9.D0)&
                *(-5.8995197352D13+rk*(1.083140865816D13+rk*(4.06356836516D11+rk*(-1.34810795364D11+rk&
                *(5.7271683D7+rk*(5.1315348D8+rk*(-2.197146D6+rk*(-534036.D0+307.D0*rk)))))))) * den(8)
        ckplm(k,3048) = 1.05509682307168D-6*(rk+10.D0)*(rk+11.D0)*(rk-8.D0)*(rk+9.D0)&
                *(1.13065007328D12+rk*(-1.11361766832D11+rk*(-1.7960528148D10+rk*(1.528234148D9+rk&
                *(9.7611577D7+rk*(-6.035848D6+rk*(-194222.D0+rk*(5972.D0+73.D0*rk)))))))) * den(9)
        ckplm(k,3049) = 0.00011025761801099D0*(rk+10.D0)*(rk-8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.210521312D10+(rk+1.D0)*rk*(-2.34553008D8+(rk+1.D0)*rk*(1.443004D6+(rk*rk+rk-2932.D0)&
                *(rk+1.D0)*rk))) * den(10)
        ckplm(k,3050) = 1.05509682307168D-6*(rk-10.D0)*(rk-8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.22262652512D12+rk*(7.1275427088D10+rk*(-2.1902239348D10+rk*(-1.081518732D9+rk&
                *(1.24673577D8+rk*(4.749192D6+rk*(-233982.D0+rk*(-5388.D0+73.D0*rk)))))))) * den(11)
        ckplm(k,3051) = -1.55161297510541D-8*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-6.928589592288D13+rk*(-9.616608721872D12+rk*(8.05979584068D11+rk*(1.29883112828D11+rk&
                *(-2.522740157D9+rk*(-5.15104408D8+rk*(1.549702D6+rk*(536492.D0+307.D0*rk)))))))) * den(12)
        ckplm(k,3052) = -1.24129038008433D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-7.699939884D11+rk*(-7.843021986D10+rk*(9.014537308D9+rk*(7.91368697D8+rk*(-3.44414D7+rk&
                *(-1.98667D6+rk*(40012.D0+713.D0*rk))))))) * den(13)
        ckplm(k,3053) = -4.7741937695551D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-1.70085549816D12+rk*(-1.03275043476D11+rk*(1.8347530456D10+rk*(6.33707085D8+rk&
                *(-5.4934835D7+rk*(-552849.D0+24979.D0*rk)))))) * den(14)
        ckplm(k,3054) = -2.76903238634196D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-2.34756036D9+rk*(-4.2802566D7+rk*(1.9732309D7+rk*(-49975.D0+rk&
                *(-29797.D0+229.D0*rk))))) * den(15)
        ckplm(k,3055) = 5.36500024853754D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-8.D0)*(rk-9.D0)*(8.8850376D7+rk*(-1.886962D6+rk*(-426483.D0+rk*(11158.D0+111.D0&
                *rk)))) * den(16)
        ckplm(k,3056) = 2.14600009941502D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-8.D0)*(rk-9.D0)*(1.417158D6+rk*(-72157.D0+rk*(-1722.D0+73.D0*rk))) &
                * den(17)
        ckplm(k,3057) = 3.57666683235836D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-8.D0)*(rk-9.D0)*(430830.D0+rk*(-26563.D0+223.D0*rk)) &
                * den(18)
        ckplm(k,3058) = -0.00002646733455945D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-8.D0)*(rk-9.D0)*(-199.D0+9.D0*rk) * den(19)
        ckplm(k,3059) = 0.00008601883731822D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,3060) = -2.9661668040765D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3061) = 7.60555590788847D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(3159.D0+142.D0*rk) * den(1)
        ckplm(k,3062) = -6.166666952342D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(1.455543D6+rk*(101465.D0+1432.D0*rk)) * den(2)
        ckplm(k,3063) = -1.8500000857026D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(-1.163709D7+rk*(-838733.D0+rk*(-5343.D0+326.D0*rk))) * den(3)
        ckplm(k,3064) = 6.166666952342D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(-6.31428318D8+rk*(-3.4796859D7+rk*(1.485449D6+rk*(89211.D0+517.D0&
                *rk)))) * den(4)
        ckplm(k,3065) = 7.16129065433265D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(8.181249804D9+rk*(2.07336186D8+rk*(-5.2699255D7+rk*(-1.8287D6+rk*(33691.D0+974.D0&
                *rk))))) * den(5)
        ckplm(k,3066) = 2.14838719629979D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(-3.628954926D10+rk*(3.82909266D8+rk*(3.77209173D8+rk*(3.608995D6+rk*(-863305.D0+rk&
                *(-14401.D0+172.D0*rk)))))) * den(6)
        ckplm(k,3067) = -3.10322595021081D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(-3.0671183844D11+rk*(1.4557582224D10+rk*(4.11920922D9+rk*(-9.2275939D7+rk*(-1.6742355D7+rk&
                *(41741.D0+rk*(16275.D0+74.D0*rk))))))) * den(7)
        ckplm(k,3068) = -4.65483892531622D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(2.337250734D12+rk&
                *(-1.9124702784D11+rk*(-3.6293165364D10+rk*(2.210466566D9+rk*(1.98125263D8+rk*(-6.85762D6+rk&
                *(-386666.D0+rk*(4334.D0+127.D0*rk)))))))) * den(8)
        ckplm(k,3069) = -2.37396785191127D-6*(rk+10.D0)*(rk+11.D0)*(-4.9505881152D11+rk&
                *(5.4975487008D10+rk*(8.514145444D9+rk*(-8.77856604D8+rk*(-5.2900351D7+rk*(4.46873D6+rk&
                *(138866.D0+rk*(-7456.D0+rk*(-119.D0+2.D0*rk))))))))) * den(9)
        ckplm(k,3070) = 0.00016538642701649D0*(rk+10.D0)*(rk-9.D0)*(8.07014208D9+(rk+1.D0)*rk&
                *(-1.76358432D8+(rk+1.D0)*rk*(1.307084D6+(rk+1.D0)*rk*(-3700.D0+3.D0*(rk+1.D0)*rk)))) * den(10)
        ckplm(k,3071) = 2.37396785191127D-6*(rk-10.D0)*(rk-9.D0)*(5.4069951936D11+rk&
                *(3.5546686944D10+rk*(-1.0787862012D10+rk*(-6.24599348D8+rk*(7.2908633D7+rk*(3.485874D6+rk&
                *(-187558.D0+(rk-32.D0)*rk*(201.D0+2.D0*rk)))))))) * den(11)
        ckplm(k,3072) = 4.65483892531622D-8*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)*(2.49019872192D12+rk&
                *(1.12853737248D11+rk*(-4.1673124732D10+rk*(-1.357267212D9+rk*(2.26470573D8+rk*(4.453722D6+rk&
                *(-413448.D0+rk*(-3318.D0+127.D0*rk)))))))) * den(12)
        ckplm(k,3073) = 3.10322595021081D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(3.170747034D11+rk*(6.10941696D9+rk*(-4.295408068D9+rk*(-2.5212019D7+rk*(1.6709525D7+rk&
                *(-54355.D0+rk*(-15757.D0+74.D0*rk))))))) * den(13)
        ckplm(k,3074) = -2.14838719629979D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(-3.629970708D10+rk*(3.57301912D8+rk*(3.61348948D8+rk*(-6.914765D6+rk*(-788720.D0+rk&
                *(15433.D0+172.D0*rk)))))) * den(14)
        ckplm(k,3075) = -7.16129065433265D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0)*(-7.92307578D9+rk*(3.07118702D8+rk*(4.7020749D7+rk*(-1.953724D6+rk*(-28821.D0+974.D0&
                *rk))))) * den(15)
        ckplm(k,3076) = -6.166666952342D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-9.D0)*(-5.95234704D8+rk*(3.7502192D7+rk*(1.220918D6+rk*(-87143.D0+517.D0&
                *rk)))) * den(16)
        ckplm(k,3077) = 1.8500000857026D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-9.D0)*(1.0804026D7+rk*(-827069.D0+rk*(6321.D0+326.D0*rk))) * den(17)
        ckplm(k,3078) = 6.166666952342D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-9.D0)*(1.35551D6+rk*(-98601.D0+1432.D0*rk)) * den(18)
        ckplm(k,3079) = -7.60555590788847D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-9.D0)*(-3017.D0+142.D0*rk) * den(19)
        ckplm(k,3080) = 2.9661668040765D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-9.D0) * den(20)
        ckplm(k,3081) = 9.88722268025501D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3082) = -1.52111118157769D-7*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(65.D0+3.D0*rk) * den(1)
        ckplm(k,3083) = 2.05555565078067D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(220705.D0+rk*(17091.D0+299.D0*rk)) * den(2)
        ckplm(k,3084) = 3.7000001714052D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(-355480.D0+rk*(-32508.D0+(rk-698.D0)*rk)) * den(3)
        ckplm(k,3085) = -2.80303043288273D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(-1.00081806D9+rk*(-8.9237856D7+rk*(-588049.D0+rk*(96114.D0+1471.D0*rk)))) * den(4)
        ckplm(k,3086) = -2.89345076942733D-11*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(1.6716724146D11+rk*(1.2472057006D10+rk*(-4.46245875D8+rk*(-4.5712685D7+rk*(-488115.D0+7969.D0&
                *rk))))) * den(5)
        ckplm(k,3087) = 7.95698961592516D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(9.00881649D10+rk*(4.76920815D9+rk*(-6.23847779D8+rk*(-3.7715625D7+rk*(447790.D0+rk&
                *(40395.D0+229.D0*rk)))))) * den(6)
        ckplm(k,3088) = 2.06881730014054D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(-4.58844372D10+rk&
                *(-1.323543D9+rk*(5.21243212D8+rk*(1.9405358D7+rk*(-1.38152D6+rk*(-55615.D0+rk*(448.D0+17.D0&
                *rk))))))) * den(7)
        ckplm(k,3089) = 2.58602162517568D-9*(rk+11.D0)*(rk+12.D0)*(4.4439470712D12+rk&
                *(2.683881384D10+rk*(-6.9365407732D10+rk*(-1.128985872D9+rk*(3.24467409D8+rk*(7.1085D6+rk&
                *(-432138.D0+rk*(-8868.D0+61.D0*rk)))))))) * den(8)
        ckplm(k,3090) = -1.75849470511946D-7*(rk+11.D0)*(7.315255584D11+rk*(-9.17458992D9+rk&
                *(-1.4444019924D10+rk*(3.9002144D7+rk*(9.6362091D7+rk*(463869.D0+rk*(-239526.D0+rk&
                *(-2094.D0+(rk+159.D0)*rk)))))))) * den(9)
        ckplm(k,3091) = -3.34113993972698D-7*(-3.9947203296D12+(rk+1.D0)*rk*(9.724374288D10+(rk+1.D0)&
                *rk*(-8.44290804D8+(rk+1.D0)*rk*(3.086008D6+(rk*rk+rk-4165.D0)*(rk+1.D0)*rk)))) * den(10)
        ckplm(k,3092) = -1.75849470511946D-7*(rk-10.D0)*(-7.263127872D11+rk*(1.9448748576D10+rk&
                *(1.3991037D10+rk*(-3.3709912D8+rk*(-9.053415D7+rk*(1.848273D6+rk*(220500.D0+rk&
                *(-3330.D0+(rk-150.D0)*rk)))))))) * den(11)
        ckplm(k,3093) = 2.58602162517568D-9*(rk-10.D0)*(rk-11.D0)*(4.3491887712D12+rk&
                *(-1.60922874816D11+rk*(-6.4109024796D10+rk*(2.347441544D9+rk*(2.82757489D8+rk*(-9.511684D6+rk&
                *(-368354.D0+rk*(9356.D0+61.D0*rk)))))))) * den(12)
        ckplm(k,3094) = 2.06881730014054D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(4.406038182D10+rk&
                *(-2.302567914D9+rk*(-4.55300531D8+rk*(2.4366923D7+rk*(1.09732D6+rk*(-57946.D0+rk&
                *(-329.D0+17.D0*rk))))))) * den(13)
        ckplm(k,3095) = 7.95698961592516D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(8.473323222D10+rk*(-5.902166274D9+rk*(-5.08414679D8+rk*(3.9107415D7+rk*(249250.D0+rk&
                *(-39021.D0+229.D0*rk)))))) * den(14)
        ckplm(k,3096) = -2.89345076942733D-11*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(-1.5429415518D11+rk*(1.3229403006D10+rk*(3.121162D8+rk*(-4.3680535D7+rk*(527960.D0+7969.D0&
                *rk))))) * den(15)
        ckplm(k,3097) = -2.80303043288273D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(-9.12262896D8+rk*(8.77793D7+rk*(-867565.D0+rk*(-90230.D0+1471.D0*rk)))) * den(16)
        ckplm(k,3098) = 3.7000001714052D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(323671.D0+rk*(-31109.D0+(rk+701.D0)*rk)) * den(17)
        ckplm(k,3099) = 2.05555565078067D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(203913.D0+rk*(-16493.D0+299.D0*rk)) * den(18)
        ckplm(k,3100) = -1.52111118157769D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(-62.D0+3.D0*rk) * den(19)
        ckplm(k,3101) = 9.88722268025501D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3102) = -3.18942667105D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3103) = 2.45340513157693D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(1573.D0+74.D0*rk) * den(1)
        ckplm(k,3104) = -6.63082467993764D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(3.218963D6+rk*(267825.D0+5272.D0*rk)) * den(2)
        ckplm(k,3105) = 2.98387110597194D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(2.46961D6+rk*(262315.D0+rk*(7945.D0+54.D0*rk))) * den(3)
        ckplm(k,3106) = 9.04203365446041D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(-2.029665858D9+rk*(-2.35349799D8+rk*(-6.881567D6+rk*(37551.D0+2393.D0*rk)))) * den(4)
        ckplm(k,3107) = -7.23362692356833D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(-4.985418042D10+rk*(-5.716512014D9+rk*(-8.4369735D7+rk*(9.20506D6+rk*(286515.D0+1054.D0&
                *rk))))) * den(5)
        ckplm(k,3108) = -7.95698961592516D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(7.461686034D10+rk&
                *(7.916087634D9+rk*(-1.14268103D8+rk*(-3.3231405D7+rk*(-725165.D0+rk*(11391.D0+268.D0*rk)))))) &
                * den(6)
        ckplm(k,3109) = -1.03440865007027D-9*(rk+12.D0)*(rk+13.D0)*(-8.237488644D10+rk&
                *(-7.764465624D9+rk*(4.6412324D8+rk*(5.9276259D7+rk*(325535.D0+rk*(-82341.D0+rk*(-1435.D0+6.D0&
                *rk))))))) * den(7)
        ckplm(k,3110) = 5.17204325035136D-10*(rk+12.D0)*(-2.112087978D12+rk*(-1.7485756296D11+rk&
                *(2.1325547452D10+rk*(2.013630402D9+rk*(-3.9566499D7+rk*(-5.82414D6+rk*(-39942.D0+rk&
                *(3258.D0+29.D0*rk)))))))) * den(8)
        ckplm(k,3111) = 8.79247352559731D-9*(1.44115747776D12+rk*(1.07806965408D11+rk&
                *(-2.1128339676D10+rk*(-1.663206764D9+rk*(8.8247649D7+rk*(7.69377D6+rk*(-83454.D0+rk&
                *(-10656.D0+rk*(-39.D0+2.D0*rk))))))))) * den(9)
        ckplm(k,3112) = -1.83762696684984D-6*(6.60284352D9+(rk+1.D0)*rk*(-1.17234144D8+(rk+1.D0)*rk&
                *(700548.D0+(rk*rk+rk-1580.D0)*(rk+1.D0)*rk))) * den(10)
        ckplm(k,3113) = -8.79247352559731D-9*(-1.31396586048D12+rk*(1.44759929184D11+rk&
                *(1.5687200388D10+rk*(-1.937961188D9+rk*(-4.8896967D7+rk*(7.973154D6+rk*(10122.D0+rk&
                *(-10272.D0+rk*(57.D0+2.D0*rk))))))))) * den(11)
        ckplm(k,3114) = -5.17204325035136D-10*(rk-11.D0)*(-1.91795228352D12+rk*(2.11338359136D11+rk&
                *(1.5104831916D10+rk*(-2.114566244D9+rk*(-1.1156929D7+rk*(5.517694D6+rk*(-61936.D0+rk&
                *(-3026.D0+29.D0*rk)))))))) * den(12)
        ckplm(k,3115) = 1.03440865007027D-9*(rk-11.D0)*(rk-12.D0)*(7.42051674D10+rk*(-8.51658852D9+rk&
                *(-2.89049432D8+rk*(5.7179619D7+rk*(-715505.D0+rk*(-73605.D0+rk*(1477.D0+6.D0*rk))))))) &
                * den(13)
        ckplm(k,3116) = 7.95698961592516D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(6.661899972D10+rk&
                *(-8.047885632D9+rk*(-1.9034768D7+rk*(3.0222195D7+rk*(-778100.D0+rk*(-9783.D0+268.D0*rk)))))) &
                * den(14)
        ckplm(k,3117) = 7.23362692356833D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(4.423095774D10+rk*(-5.521298154D9+rk*(1.10276365D8+rk*(8.06954D6+rk*(-281245.D0+1054.D0&
                *rk))))) * den(15)
        ckplm(k,3118) = -9.04203365446041D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(-1.801232784D9+rk*(2.21483584D8+rk*(-6.979862D6+rk*(-27979.D0+2393.D0*rk)))) * den(16)
        ckplm(k,3119) = -2.98387110597194D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(-2.215186D6+rk*(246587.D0+rk*(-7783.D0+54.D0*rk))) * den(17)
        ckplm(k,3120) = 6.63082467993764D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(2.95641D6+rk*(-257281.D0+5272.D0*rk)) * den(18)
        ckplm(k,3121) = -2.45340513157693D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(-1499.D0+74.D0*rk) * den(19)
        ckplm(k,3122) = 3.18942667105D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3123) = 9.96695834703126D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3124) = -1.02225213815705D-11*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(1404.D0+67.D0*rk) * den(1)
        ckplm(k,3125) = 4.14426542496102D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(2.272836D6+rk*(199129.D0+4223.D0*rk)) * den(2)
        ckplm(k,3126) = -4.14426542496102D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(9.2238D6+rk*(1.085738D6+rk*(39249.D0+409.D0*rk))) * den(3)
        ckplm(k,3127) = -3.13959501890987D-14*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(-3.5084691576D10+rk*(-4.818715938D9+rk*(-2.08055771D8+rk*(-2.502426D6+11831.D0*rk)))) * den(4)
        ckplm(k,3128) = 9.04203365446041D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(-2.71349568D10+rk&
                *(-4.004891736D9+rk*(-1.6625357D8+rk*(314525.D0+rk*(128870.D0+1471.D0*rk))))) * den(5)
        ckplm(k,3129) = 1.49193555298597D-11*(rk+13.D0)*(rk+14.D0)*(2.999846016D10+rk&
                *(4.529708616D9+rk*(1.42934798D8+rk*(-7.345045D6+rk*(-454035.D0+rk*(-4691.D0+37.D0*rk)))))) &
                * den(6)
        ckplm(k,3130) = -1.43667868065315D-11*(rk+13.D0)*(4.8403547136D11+rk*(7.3010881296D10+rk&
                *(1.080602964D9+rk*(-2.92357352D8+rk*(-1.3330905D7+rk*(37429.D0+rk*(8421.D0+67.D0*rk))))))) &
                * den(7)
        ckplm(k,3131) = -5.38754505244933D-12*(-1.751899968D13+rk*(-2.63586849336D12+rk&
                *(1.6411965972D10+rk*(1.8354738892D10+rk*(6.21299651D8+rk*(-2.282084D7+rk*(-1.210762D6+rk&
                *(-4532.D0+179.D0*rk)))))))) * den(8)
        ckplm(k,3132) = 3.66353063566554D-10*(-2.5809009408D11+rk*(-1.7984438544D10+rk&
                *(2.596959644D9+rk*(1.87631556D8+rk*(-6.563151D6+rk*(-521976.D0+rk*(2226.D0+(rk+324.D0)&
                *rk))))))) * den(9)
        ckplm(k,3133) = 1.16011803462742D-9*(7.923412224D10+(rk+1.D0)*rk*(-1.005923088D9+(rk+1.D0)*rk&
                *(3.971244D6+(rk*rk+rk-5012.D0)*(rk+1.D0)*rk))) * den(10)
        ckplm(k,3134) = 3.66353063566554D-10*(-2.3770236672D11+rk*(2.2591831536D10+rk&
                *(1.999932444D9+rk*(-2.08631164D8+rk*(-3.931151D6+rk*(528584.D0+rk*(-14.D0+(rk-316.D0)&
                *rk))))))) * den(11)
        ckplm(k,3135) = -5.38754505244933D-12*(-1.488443104512D13+rk*(2.616220280016D12+rk&
                *(-3.4714305644D10+rk*(-1.5665378484D10+rk*(7.17413571D8+rk*(1.5661464D7+rk*(-1.174026D6+rk&
                *(5964.D0+179.D0*rk)))))))) * den(12)
        ckplm(k,3136) = -1.43667868065315D-11*(rk-12.D0)*(-4.123841904D11+rk*(7.002606402D10+rk&
                *(-1.877440208D9+rk*(-2.38825517D8+rk*(1.339408D7+rk*(-11690.D0+rk*(-7952.D0+67.D0*rk))))))) &
                * den(13)
        ckplm(k,3137) = 1.49193555298597D-11*(rk-12.D0)*(rk-13.D0)*(2.561858208D10+rk&
                *(-4.223596348D9+rk*(1.62293188D8+rk*(5.576555D6+rk*(-430025.D0+rk*(4913.D0+37.D0*rk)))))) &
                * den(14)
        ckplm(k,3138) = 9.04203365446041D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(2.329650576D10+rk&
                *(-3.671949146D9+rk*(1.66438635D8+rk*(-186245.D0+rk*(-121515.D0+1471.D0*rk))))) * den(15)
        ckplm(k,3139) = -3.13959501890987D-14*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(-3.0471517152D10+rk*(4.410158998D9+rk*(-2.00477507D8+rk*(2.54975D6+11831.D0*rk)))) * den(16)
        ckplm(k,3140) = -4.14426542496102D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(-8.176902D6+rk*(1.008467D6+rk*(-38022.D0+409.D0*rk))) * den(17)
        ckplm(k,3141) = 4.14426542496102D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(2.07793D6+rk*(-190683.D0+4223.D0*rk)) * den(18)
        ckplm(k,3142) = -1.02225213815705D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(-1337.D0+67.D0*rk) * den(19)
        ckplm(k,3143) = 9.96695834703126D-11*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3144) = -3.02029040819129D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3145) = 2.3233003139933D-13*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(2197.D0+106.D0*rk) * den(1)
        ckplm(k,3146) = -1.88375701134592D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(2.082249D6+rk*(189611.D0+4232.D0*rk)) * den(2)
        ckplm(k,3147) = 3.13959501890987D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(5.892861D7+rk*(7.473139D6+rk*(301473.D0+3782.D0*rk))) * den(3)
        ckplm(k,3148) = -3.13959501890987D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(1.948387818D9+rk&
                *(3.01327423D8+rk*(1.5931615D7+rk*(317129.D0+1535.D0*rk)))) * den(4)
        ckplm(k,3149) = -2.2605084136151D-13*(rk+14.D0)*(rk+15.D0)*(-6.79667976D9+rk&
                *(-1.188898918D9+rk*(-7.0528615D7+rk*(-1.40678D6+rk*(6115.D0+318.D0*rk))))) * den(5)
        ckplm(k,3150) = 7.53502804538368D-14*(rk+14.D0)*(-4.1096166696D11+rk*(-7.7771390046D10+rk&
                *(-4.668869513D9+rk*(-4.7578755D7+rk*(4.370485D6+rk*(123921.D0+628.D0*rk)))))) * den(6)
        ckplm(k,3151) = 3.26517881966626D-13*(1.59354297648D12+rk*(3.18934745028D11+rk&
                *(1.7955155652D10+rk*(-2.24722561D8+rk*(-4.9846965D7+rk*(-1.222753D6+rk*(1953.D0+206.D0&
                *rk))))))) * den(7)
        ckplm(k,3152) = -2.12236623278307D-12*(2.70910584D11+rk*(3.620786076D10+rk*(7.6092548D7+rk&
                *(-1.55193458D8+rk*(-4.850965D6+rk*(85645.D0+rk*(3857.D0+13.D0*rk))))))) * den(8)
        ckplm(k,3153) = -3.60802259573122D-11*(-1.654515072D10+rk*(-1.048131156D9+rk*(1.10653228D8+rk&
                *(7.099673D6+rk*(-156275.D0+rk*(-10759.D0+rk*(7.D0+2.D0*rk))))))) * den(9)
        ckplm(k,3154) = 7.54076722507825D-9*(rk-18.D0)*(rk+19.D0)*(228480.D0+(rk*rk+rk-1346.D0)&
                *(rk+1.D0)*rk) * den(10)
        ckplm(k,3155) = 3.60802259573122D-11*(1.539361152D10+rk*(-1.247567316D9+rk*(-8.8524212D7+rk&
                *(7.617113D6+rk*(102445.D0+rk*(-10759.D0+rk*(7.D0+2.D0*rk))))))) * den(11)
        ckplm(k,3156) = 2.12236623278307D-12*(-2.3492907648D11+rk*(3.5609904324D10+rk&
                *(-5.11768264D8+rk*(-1.35009833D8+rk*(5.22179D6+rk*(62776.D0+rk*(-3766.D0+13.D0*rk))))))) &
                * den(12)
        ckplm(k,3157) = -3.26517881966626D-13*(-1.2927394872D12+rk*(2.8254352986D11+rk&
                *(-1.8342494044D10+rk*(-3.7594081D7+rk*(4.3711115D7+rk*(-1.230145D6+rk*(-511.D0+206.D0&
                *rk))))))) * den(13)
        ckplm(k,3158) = -7.53502804538368D-14*(rk-13.D0)*(-3.3780732048D11+rk*(6.8593253388D10+rk&
                *(-4.501140128D9+rk*(6.3834045D7+rk*(3.7603D6+rk*(-120153.D0+628.D0*rk)))))) * den(14)
        ckplm(k,3159) = 2.2605084136151D-13*(rk-13.D0)*(rk-14.D0)*(5.67689688D9+rk*(-1.052084898D9+rk&
                *(6.6274765D7+rk*(-1.42806D6+rk*(-4525.D0+318.D0*rk))))) * den(15)
        ckplm(k,3160) = 3.13959501890987D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(1.662676416D9+rk&
                *(-2.7040944D8+rk*(1.4989438D7+rk*(-310989.D0+1535.D0*rk)))) * den(16)
        ckplm(k,3161) = -3.13959501890987D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(-5.1753162D7+rk*(6.881539D6+rk*(-290127.D0+3782.D0*rk))) * den(17)
        ckplm(k,3162) = 1.88375701134592D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(1.89687D6+rk*(-181147.D0+4232.D0*rk)) * den(18)
        ckplm(k,3163) = -2.3233003139933D-13*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(-2091.D0+106.D0*rk) * den(19)
        ckplm(k,3164) = 3.02029040819129D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3165) = 8.88320708291556D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(rk+20.D0) * den(0)
        ckplm(k,3166) = -2.73329448705094D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(637.D0+31.D0*rk) * den(1)
        ckplm(k,3167) = 1.10809235961525D-15*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(1.397823D6+rk*(131125.D0+3037.D0*rk)) * den(2)
        ckplm(k,3168) = -5.27663028388213D-16*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(1.5933036D8+rk&
                *(2.1373588D7+rk*(929358.D0+12977.D0*rk))) * den(3)
        ckplm(k,3169) = 1.56979750945493D-14*(rk+15.D0)*(rk+16.D0)*(2.01629148D8+rk*(3.4041584D7+rk&
                *(2.044985D6+rk*(50494.D0+409.D0*rk)))) * den(4)
        ckplm(k,3170) = 4.52101682723021D-13*(rk+15.D0)*(-1.9727064D8+rk*(-3.9022666D7+rk&
                *(-2.825915D6+rk*(-87885.D0+(rk-955.D0)*rk)))) * den(5)
        ckplm(k,3171) = -7.53502804538368D-14*(-2.628495576D10+rk*(-5.836235226D9+rk&
                *(-4.70642183D8+rk*(-1.5249405D7+rk*(-65090.D0+rk*(5751.D0+73.D0*rk)))))) * den(6)
        ckplm(k,3172) = 6.53035763933252D-13*(-3.93168816D9+rk*(-6.76033596D8+rk*(-3.3646326D7+rk&
                *(73010.D0+rk*(43925.D0+(rk+826.D0)*rk))))) * den(7)
        ckplm(k,3173) = 8.16294704916565D-14*(3.6775368D10+rk*(4.26487968D9+rk*(2.9453674D7+rk&
                *(-1.0512525D7+rk*(-277895.D0+rk*(1965.D0+61.D0*rk)))))) * den(8)
        ckplm(k,3174) = -3.26517881966626D-13*(9.89295552D9+rk*(5.54582886D8+rk*(-4.2012521D7+rk&
                *(-2.29119D6+rk*(29320.D0+(rk+1584.D0)*rk))))) * den(9)
        ckplm(k,3175) = -4.87444552364463D-12*(-6.6419136D8+(rk+1.D0)*rk*(3.894882D6+(rk&
                *rk+rk-5363.D0)*(rk+1.D0)*rk)) * den(10)
        ckplm(k,3176) = -3.26517881966626D-13*(9.29867904D9+rk*(-6.31624992D8+rk*(-3.4978856D7+rk&
                *(2.39265D6+rk*(21415.D0+(rk-1578.D0)*rk))))) * den(11)
        ckplm(k,3177) = 8.16294704916565D-14*(3.255017472D10+rk*(-4.175555796D9+rk*(5.9305144D7+rk&
                *(9.382515D6+rk*(-286805.D0+rk*(-1599.D0+61.D0*rk)))))) * den(12)
        ckplm(k,3178) = 6.53035763933252D-13*(-3.2893308D9+rk*(6.0869349D8+rk*(-3.3610051D7+rk&
                *(94450.D0+rk*(39810.D0+(rk-820.D0)*rk))))) * den(13)
        ckplm(k,3179) = -7.53502804538368D-14*(-2.090418408D10+rk*(4.940410398D9+rk*(-4.25340923D8+rk&
                *(1.4932995D7+rk*(-92750.D0+rk*(-5313.D0+73.D0*rk)))))) * den(14)
        ckplm(k,3180) = 4.52101682723021D-13*(rk-14.D0)*(1.6098696D8+rk*(-3.3630666D7+rk*(2.568D6+rk&
                *(-84055.D0+(rk+960.D0)*rk)))) * den(15)
        ckplm(k,3181) = 1.56979750945493D-14*(rk-14.D0)*(rk-15.D0)*(1.69582464D8+rk*(-3.010146D7+rk&
                *(1.895957D6+rk*(-48858.D0+409.D0*rk)))) * den(16)
        ckplm(k,3182) = -5.27663028388213D-16*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(-1.38873153D8+rk&
                *(1.9553803D7+rk*(-890427.D0+12977.D0*rk))) * den(17)
        ckplm(k,3183) = 1.10809235961525D-15*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(1.269735D6+rk*(-125051.D0+3037.D0*rk)) * den(18)
        ckplm(k,3184) = -2.73329448705094D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(-606.D0+31.D0*rk) * den(19)
        ckplm(k,3185) = 8.88320708291556D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(rk-19.D0) * den(20)
        ckplm(k,3186) = -2.5380591665473D-15*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0) &
                * den(0)
        ckplm(k,3187) = 3.25392200839398D-16*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(1755.D0+86.D0*rk) * den(1)
        ckplm(k,3188) = -2.63831514194106D-17*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(2.207235D6+rk&
                *(211957.D0+5048.D0*rk)) * den(2)
        ckplm(k,3189) = 2.63831514194106D-17*(rk+16.D0)*(rk+17.D0)*(1.3666263D8+rk*(1.9150283D7+rk&
                *(879273.D0+13174.D0*rk))) * den(3)
        ckplm(k,3190) = -1.49504524709994D-16*(rk+16.D0)*(1.02665457D9+rk*(1.85386797D8+rk&
                *(1.2169913D7+rk*(340707.D0+3373.D0*rk)))) * den(4)
        ckplm(k,3191) = 1.07643257791195D-15*(4.48579512D9+rk*(9.74036922D8+rk*(8.0427125D7+rk&
                *(3.08578D6+rk*(52495.D0+278.D0*rk))))) * den(5)
        ckplm(k,3192) = 1.61464886686793D-14*(-4.8644712D8+rk*(-8.9343662D7+rk*(-5.815675D6+rk&
                *(-149830.D0+rk*(-845.D0+12.D0*rk))))) * den(6)
        ckplm(k,3193) = -7.77423528491967D-15*(-1.41281496D9+rk*(-2.02580346D8+rk*(-8.443075D6+rk&
                *(-23740.D0+rk*(4315.D0+46.D0*rk))))) * den(7)
        ckplm(k,3194) = -5.83067646368975D-14*(2.32948512D8+rk*(2.27385D7+rk*(204700.D0+rk&
                *(-28993.D0+(rk-568.D0)*rk)))) * den(8)
        ckplm(k,3195) = 1.74920293910693D-13*(8.6598272D7+rk*(4.156182D6+rk*(-218409.D0+rk&
                *(-9476.D0+rk*(61.D0+2.D0*rk))))) * den(9)
        ckplm(k,3196) = -2.43722276182232D-12*(6.325632D6+5.D0*(rk+1.D0)*rk*(-4630.D0+3.D0*(rk+1.D0)&
                *rk)) * den(10)
        ckplm(k,3197) = -1.74920293910693D-13*(-8.2233216D7+rk*(4.564338D6+rk*(189635.D0+rk&
                *(-9700.D0+rk*(-51.D0+2.D0*rk))))) * den(11)
        ckplm(k,3198) = 5.83067646368975D-14*(-2.10443136D8+rk*(2.2244398D7+rk*(-288261.D0+rk&
                *(-26711.D0+(rk+573.D0)*rk)))) * den(12)
        ckplm(k,3199) = 7.77423528491967D-15*(1.21864968D9+rk*(-1.85782446D8+rk*(8.346425D6+rk&
                *(-40540.D0+rk*(-4085.D0+46.D0*rk))))) * den(13)
        ckplm(k,3200) = -1.61464886686793D-14*(4.0277016D8+rk*(-7.8158362D7+rk*(5.371375D6+rk&
                *(-146330.D0+rk*(905.D0+12.D0*rk))))) * den(14)
        ckplm(k,3201) = -1.07643257791195D-15*(-3.58915176D9+rk*(8.22231422D8+rk*(-7.1481975D7+rk&
                *(2.87858D6+rk*(-51105.D0+278.D0*rk))))) * den(15)
        ckplm(k,3202) = 1.49504524709994D-16*(rk-15.D0)*(8.53100352D8+rk*(-1.620556D8+rk&
                *(1.116803D7+rk*(-327215.D0+3373.D0*rk)))) * den(16)
        ckplm(k,3203) = -2.63831514194106D-17*(rk-15.D0)*(rk-16.D0)*(-1.18378446D8+rk*(1.7431259D7+rk&
                *(-839751.D0+13174.D0*rk))) * den(17)
        ckplm(k,3204) = 2.63831514194106D-17*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(2.000326D6+rk&
                *(-201861.D0+5048.D0*rk)) * den(18)
        ckplm(k,3205) = -3.25392200839398D-16*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(-1669.D0+86.D0*rk) * den(19)
        ckplm(k,3206) = 2.5380591665473D-15*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0) &
                * den(20)
        ckplm(k,3207) = 7.05016435152029D-17*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3208) = -2.16928133892932D-17*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(832.D0+41.D0*rk) &
                * den(1)
        ckplm(k,3209) = 2.9314612688234D-19*(rk+17.D0)*(rk+18.D0)*(7.134784D6+rk*(698145.D0+16991.D0&
                *rk)) * den(2)
        ckplm(k,3210) = -5.27663028388213D-18*(rk+17.D0)*(2.777792D7+rk*(4.029754D6+rk&
                *(192809.D0+3037.D0*rk))) * den(3)
        ckplm(k,3211) = 7.47522623549968D-18*(9.37838208D8+rk*(1.78558854D8+rk*(1.2513781D7+rk&
                *(380934.D0+4223.D0*rk)))) * den(4)
        ckplm(k,3212) = -1.19603619767995D-16*(1.28082024D8+rk*(2.1911882D7+rk*(1.344193D6+rk&
                *(34402.D0+299.D0*rk)))) * den(5)
        ckplm(k,3213) = 5.98018098839974D-17*(4.56534504D8+rk*(6.6524458D7+rk*(3.264617D6+rk&
                *(58478.D0+223.D0*rk)))) * den(6)
        ckplm(k,3214) = 1.55484705698393D-15*(-2.6451288D7+rk*(-3.03235D6+rk*(-98245.D0+rk&
                *(-410.D0+13.D0*rk)))) * den(7)
        ckplm(k,3215) = -5.83067646368975D-16*(-9.22178D7+rk*(-7.259978D6+rk*(-65167.D0+rk&
                *(3982.D0+43.D0*rk)))) * den(8)
        ckplm(k,3216) = -2.59141176163989D-16*(2.39833008D8+rk*(9.424394D6+rk*(-323699.D0+rk&
                *(-9926.D0+23.D0*rk)))) * den(9)
        ckplm(k,3217) = 2.70802529091368D-14*(2.372112D6+(rk*rk+rk-4898.D0)*(rk+1.D0)*rk) * den(10)
        ckplm(k,3218) = -2.59141176163989D-16*(2.30094864D8+rk*(-1.0041922D7+rk*(-293783.D0+rk&
                *(10018.D0+23.D0*rk)))) * den(11)
        ckplm(k,3219) = -5.83067646368975D-16*(-8.5026928D7+rk*(7.11787D6+rk*(-76855.D0+rk&
                *(-3810.D0+43.D0*rk)))) * den(12)
        ckplm(k,3220) = 1.55484705698393D-15*(-2.351676D7+rk*(2.837142D6+rk*(-96937.D0+rk&
                *(462.D0+13.D0*rk)))) * den(13)
        ckplm(k,3221) = 5.98018098839974D-17*(3.93216408D8+rk*(-6.0169766D7+rk*(3.090521D6+rk&
                *(-57586.D0+223.D0*rk)))) * den(14)
        ckplm(k,3222) = -1.19603619767995D-16*(1.07480232D8+rk*(-1.9325506D7+rk*(1.242781D6+rk&
                *(-33206.D0+299.D0*rk)))) * den(15)
        ckplm(k,3223) = 7.47522623549968D-18*(7.71416424D8+rk*(-1.54657202D8+rk*(1.1396317D7+rk&
                *(-364042.D0+4223.D0*rk)))) * den(16)
        ckplm(k,3224) = -5.27663028388213D-18*(rk-16.D0)*(-2.3937938D7+rk*(3.653247D6+rk&
                *(-183698.D0+3037.D0*rk))) * den(17)
        ckplm(k,3225) = 2.9314612688234D-19*(rk-16.D0)*(rk-17.D0)*(6.45363D6+rk*(-664163.D0+16991.D0&
                *rk)) * den(18)
        ckplm(k,3226) = -2.16928133892932D-17*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(-791.D0+41.D0*rk) &
                * den(19)
        ckplm(k,3227) = 7.05016435152029D-17*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3228) = -1.90544982473521D-18*(rk+18.D0)*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3229) = 1.4657306344117D-19*(rk+18.D0)*(rk+19.D0)*(3757.D0+186.D0*rk) * den(1)
        ckplm(k,3230) = -1.4657306344117D-19*(rk+18.D0)*(491011.D0+rk*(48789.D0+1208.D0*rk)) * den(2)
        ckplm(k,3231) = 1.31915757097053D-18*(4.29743D6+rk*(641171.D0+rk*(31681.D0+518.D0*rk))) &
                * den(3)
        ckplm(k,3232) = -7.47522623549968D-18*(2.383386D6+rk*(333339.D0+rk*(15286.D0+229.D0*rk))) &
                * den(4)
        ckplm(k,3233) = 2.99009049419987D-17*(1.440738D6+rk*(182015.D0+rk*(7365.D0+94.D0*rk))) &
                * den(5)
        ckplm(k,3234) = -2.99009049419987D-17*(2.805498D6+rk*(303847.D0+rk*(9963.D0+92.D0*rk))) &
                * den(6)
        ckplm(k,3235) = 3.88711764245983D-16*(349722.D0+rk*(30017.D0+rk*(667.D0+2.D0*rk))) * den(7)
        ckplm(k,3236) = 5.83067646368975D-16*(-322100.D0+rk*(-19144.D0+rk*(-139.D0+3.D0*rk))) * den(8)
        ckplm(k,3237) = -1.10134999869695D-15*(-204282.D0+rk*(-6143.D0+rk*(123.D0+2.D0*rk))) * den(9)
        ckplm(k,3238) = 2.30182149727663D-13*(rk*rk+rk-1026.D0) * den(10)
        ckplm(k,3239) = 1.10134999869695D-15*(198018.D0+rk*(-6383.D0+rk*(-117.D0+2.D0*rk))) * den(11)
        ckplm(k,3240) = -5.83067646368975D-16*(303098.D0+rk*(-18857.D0+rk*(148.D0+3.D0*rk))) * den(12)
        ckplm(k,3241) = -3.88711764245983D-16*(-320370.D0+rk*(28689.D0+rk*(-661.D0+2.D0*rk))) &
                * den(13)
        ckplm(k,3242) = 2.99009049419987D-17*(-2.511522D6+rk*(284197.D0+rk*(-9687.D0+92.D0*rk))) &
                * den(14)
        ckplm(k,3243) = -2.99009049419987D-17*(-1.265994D6+rk*(167567.D0+rk*(-7083.D0+94.D0*rk))) &
                * den(15)
        ckplm(k,3244) = 7.47522623549968D-18*(-2.065104D6+rk*(303454.D0+rk*(-14599.D0+229.D0*rk))) &
                * den(16)
        ckplm(k,3245) = -1.31915757097053D-18*(-3.687422D6+rk*(579363.D0+rk*(-30127.D0+518.D0*rk))) &
                * den(17)
        ckplm(k,3246) = 1.4657306344117D-19*(rk-17.D0)*(443430.D0+rk*(-46373.D0+1208.D0*rk)) * den(18)
        ckplm(k,3247) = -1.4657306344117D-19*(rk-17.D0)*(rk-18.D0)*(-3571.D0+186.D0*rk) * den(19)
        ckplm(k,3248) = 1.90544982473521D-18*(rk-17.D0)*(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3249) = 5.01434164404003D-20*(rk+19.D0)*(rk+20.D0) * den(0)
        ckplm(k,3250) = -5.14291450670773D-21*(rk+19.D0)*(3159.D0+157.D0*rk) * den(1)
        ckplm(k,3251) = 1.4657306344117D-19*(16227.D0+rk*(1633.D0+41.D0*rk)) * den(2)
        ckplm(k,3252) = -8.79438380647021D-19*(13220.D0+rk*(1284.D0+31.D0*rk)) * den(3)
        ckplm(k,3253) = 1.24587103924995D-18*(32586.D0+rk*(2977.D0+67.D0*rk)) * den(4)
        ckplm(k,3254) = -1.79405429651992D-16*(603.D0+(rk+50.D0)*rk) * den(5)
        ckplm(k,3255) = 2.69108144477988D-16*(853.D0+(rk+61.D0)*rk) * den(6)
        ckplm(k,3256) = -2.59141176163989D-16*(1542.D0+(rk+88.D0)*rk) * den(7)
        ckplm(k,3257) = 9.71779410614959D-17*(6000.D0+(rk+239.D0)*rk) * den(8)
        ckplm(k,3258) = 1.29570588081994D-16*(rk-151.D0)*(rk+37.D0) * den(9)
        ckplm(k,3259) = -2.37546078150323D-16*(rk*rk+rk-3249.D0) * den(10)
        ckplm(k,3260) = 1.29570588081994D-16*(rk+152.D0)*(rk-36.D0) * den(11)
        ckplm(k,3261) = 9.71779410614959D-17*(5762.D0+(rk-237.D0)*rk) * den(12)
        ckplm(k,3262) = -2.59141176163989D-16*(1455.D0+(rk-86.D0)*rk) * den(13)
        ckplm(k,3263) = 2.69108144477988D-16*(793.D0+(rk-59.D0)*rk) * den(14)
        ckplm(k,3264) = -1.79405429651992D-16*(554.D0+(rk-48.D0)*rk) * den(15)
        ckplm(k,3265) = 1.24587103924995D-18*(29676.D0+rk*(-2843.D0+67.D0*rk)) * den(16)
        ckplm(k,3266) = -8.79438380647021D-19*(11967.D0+rk*(-1222.D0+31.D0*rk)) * den(17)
        ckplm(k,3267) = 1.4657306344117D-19*(14635.D0+rk*(-1551.D0+41.D0*rk)) * den(18)
        ckplm(k,3268) = -5.14291450670773D-21*(rk-18.D0)*(-3002.D0+157.D0*rk) * den(19)
        ckplm(k,3269) = 5.01434164404003D-20*(rk-18.D0)*(rk-19.D0) * den(20)
        ckplm(k,3270) = -1.28572862667693D-21*(rk+20.D0) * den(0)
        ckplm(k,3271) = 1.28572862667693D-21*(361.D0+18.D0*rk) * den(1)
        ckplm(k,3272) = -2.44288439068617D-20*(163.D0+8.D0*rk) * den(2)
        ckplm(k,3273) = 7.32865317205851D-20*(295.D0+14.D0*rk) * den(3)
        ckplm(k,3274) = -1.24587103924995D-18*(67.D0+3.D0*rk) * den(4)
        ckplm(k,3275) = 4.98348415699979D-18*(49.D0+2.D0*rk) * den(5)
        ckplm(k,3276) = -4.98348415699979D-18*(113.D0+4.D0*rk) * den(6)
        ckplm(k,3277) = 4.98348415699979D-18*(211.D0+6.D0*rk) * den(7)
        ckplm(k,3278) = -3.23926470204986D-17*(rk+50.D0) * den(8)
        ckplm(k,3279) = 1.07975490068329D-17*(193.D0+2.D0*rk) * den(9)
        ckplm(k,3280) = -2.25668774242807D-15 * den(10)
        ckplm(k,3281) = -1.07975490068329D-17*(-191.D0+2.D0*rk) * den(11)
        ckplm(k,3282) = 3.23926470204986D-17*(rk-49.D0) * den(12)
        ckplm(k,3283) = -4.98348415699979D-18*(-205.D0+6.D0*rk) * den(13)
        ckplm(k,3284) = 4.98348415699979D-18*(-109.D0+4.D0*rk) * den(14)
        ckplm(k,3285) = -4.98348415699979D-18*(-47.D0+2.D0*rk) * den(15)
        ckplm(k,3286) = 1.24587103924995D-18*(-64.D0+3.D0*rk) * den(16)
        ckplm(k,3287) = -7.32865317205851D-20*(-281.D0+14.D0*rk) * den(17)
        ckplm(k,3288) = 2.44288439068617D-20*(-155.D0+8.D0*rk) * den(18)
        ckplm(k,3289) = -1.28572862667693D-21*(-343.D0+18.D0*rk) * den(19)
        ckplm(k,3290) = 1.28572862667693D-21*(rk-19.D0) * den(20)
        ckplm(k,3291) = 3.21432156669233D-23 * den(0)
        ckplm(k,3292) = -6.42864313338466D-22 * den(1)
        ckplm(k,3293) = 6.10721097671543D-21 * den(2)
        ckplm(k,3294) = -3.66432658602925D-20 * den(3)
        ckplm(k,3295) = 1.55733879906243D-19 * den(4)
        ckplm(k,3296) = -4.98348415699979D-19 * den(5)
        ckplm(k,3297) = 1.24587103924995D-18 * den(6)
        ckplm(k,3298) = -2.49174207849989D-18 * den(7)
        ckplm(k,3299) = 4.04908087756233D-18 * den(8)
        ckplm(k,3300) = -5.39877450341644D-18 * den(9)
        ckplm(k,3301) = 5.93865195375808D-18 * den(10)
        ckplm(k,3302) = -5.39877450341644D-18 * den(11)
        ckplm(k,3303) = 4.04908087756233D-18 * den(12)
        ckplm(k,3304) = -2.49174207849989D-18 * den(13)
        ckplm(k,3305) = 1.24587103924995D-18 * den(14)
        ckplm(k,3306) = -4.98348415699979D-19 * den(15)
        ckplm(k,3307) = 1.55733879906243D-19 * den(16)
        ckplm(k,3308) = -3.66432658602925D-20 * den(17)
        ckplm(k,3309) = 6.10721097671543D-21 * den(18)
        ckplm(k,3310) = -6.42864313338466D-22 * den(19)
        ckplm(k,3311) = 3.21432156669233D-23 * den(20)

!    ckplm para l = 21
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+41.D0)*(r2k+43.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+41.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k+21.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(17) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)*(r2k+9.D0))
        den(18) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0))
        den(19) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(20) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-39.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        den(21) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-39.D0)*(r2k-3.D0)*(r2k-41.D0)*(r2k-5.D0)*(r2k-7.D0)*(r2k-9.D0))
        ckplm(k,3312) = 1.10364382748222000D7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3313) = 5.65280984807968000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+20.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,3314) = 4.34831526775360000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-1.D0)*(rk+1.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,3315) = 3.72153108501434000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,3316) = 3.34937797651291000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,3317) = 3.10578685094833000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,3318) = 2.93880906326294000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(6)
        ckplm(k,3319) = 2.82299392776489000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(7)
        ckplm(k,3320) = 2.74457742977142000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(8)
        ckplm(k,3321) = 2.69578494213104000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(9)
        ckplm(k,3322) = 2.67234333393860000D6*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(10)
        ckplm(k,3323) = 2.67234333393860000D6*(rk-10.D0)*(rk+10.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(11)
        ckplm(k,3324) = 2.69578494213104000D6*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(12)
        ckplm(k,3325) = 2.74457742977142000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*rk * den(13)
        ckplm(k,3326) = 2.82299392776489000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,3327) = 2.93880906326294000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(15)
        ckplm(k,3328) = 3.10578685094833000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(16)
        ckplm(k,3329) = 3.34937797651291000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(17)
        ckplm(k,3330) = 3.72153108501434000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(18)
        ckplm(k,3331) = 4.34831526775360000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(19)
        ckplm(k,3332) = 5.65280984807968000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(20)
        ckplm(k,3333) = 1.10364382748222000D7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-20.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(21)
        ckplm(k,3334) = -1.00331257043839000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3335) = -24471.03830337520000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-41.00000000000000000D0+19.00000000000000000D0*rk) * den(1)
        ckplm(k,3336) = -18823.87561798100000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-1.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-78.00000000000000000D0+17.00000000000000000D0*rk) * den(2)
        ckplm(k,3337) = -48331.57253265380000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-37.00000000000000000D0+5.00000000000000000D0*rk) * den(3)
        ckplm(k,3338) = -14499.47175979610000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-140.00000000000000000D0+13.00000000000000000D0*rk) * den(4)
        ckplm(k,3339) = -147894.61194992100000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk-15.D0)*(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(5)
        ckplm(k,3340) = -38166.35147094730000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-62.00000000000000000D0+3.00000000000000000D0*rk) * den(6)
        ckplm(k,3341) = -85545.27053833010000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk-1.D0)*(rk-29.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(7)
        ckplm(k,3342) = -11881.28757476810000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-216.00000000000000000D0+5.00000000000000000D0*rk) * den(8)
        ckplm(k,3343) = -35010.19405364990000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-75.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0) * den(9)
        ckplm(k,3344) = -11568.58586120610000000D0*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-230.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(10)
        ckplm(k,3345) = 11568.58586120610000000D0*(rk-10.D0)*(rk+10.D0)*(rk-1.D0)*(rk+231.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(11)
        ckplm(k,3346) = 35010.19405364990000000D0*(rk-10.D0)*(rk-11.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+76.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(12)
        ckplm(k,3347) = 11881.28757476810000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(221.00000000000000000D0+5.00000000000000000D0*rk) * den(13)
        ckplm(k,3348) = 85545.27053833010000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk+30.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0) * den(14)
        ckplm(k,3349) = 38166.35147094730000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(65.00000000000000000D0+3.00000000000000000D0*rk) * den(15)
        ckplm(k,3350) = 147894.61194992100000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk+16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,3351) = 14499.47175979610000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(153.00000000000000000D0+13.00000000000000000D0*rk) * den(17)
        ckplm(k,3352) = 48331.57253265380000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(42.00000000000000000D0+5.00000000000000000D0*rk) * den(18)
        ckplm(k,3353) = 18823.87561798100000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(95.00000000000000000D0+17.00000000000000000D0*rk) * den(19)
        ckplm(k,3354) = 24471.03830337520000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(60.00000000000000000D0+19.00000000000000000D0*rk) * den(20)
        ckplm(k,3355) = 1.00331257043839000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-20.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3356) = 43622.28567123410000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3357) = 1063.95818710327000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-164.00000000000000000D0+13.00000000000000000D0*rk) * den(1)
        ckplm(k,3358) = 163.68587493896500000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(3003.00000000000000000D0+rk&
                *(-1433.00000000000000000D0+29.00000000000000000D0*rk)) * den(2)
        ckplm(k,3359) = -420.27454376220700000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2035.00000000000000000D0+(rk+591.D0)*rk) * den(3)
        ckplm(k,3360) = -126.08236312866200000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-9730.00000000000000000D0+rk&
                *(1929.00000000000000000D0+31.00000000000000000D0*rk)) * den(4)
        ckplm(k,3361) = -6430.20051956177000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+41.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(5)
        ckplm(k,3362) = -1659.40658569336000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1147.00000000000000000D0+rk&
                *(119.00000000000000000D0+5.00000000000000000D0*rk)) * den(6)
        ckplm(k,3363) = -743.87191772460900000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-2929.00000000000000000D0+rk&
                *(219.00000000000000000D0+13.00000000000000000D0*rk)) * den(7)
        ckplm(k,3364) = -103.31554412841800000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-23220.00000000000000000D0+rk&
                *(1193.00000000000000000D0+103.00000000000000000D0*rk)) * den(8)
        ckplm(k,3365) = -304.43647003173800000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-8400.00000000000000000D0+rk&
                *(263.00000000000000000D0+37.00000000000000000D0*rk)) * den(9)
        ckplm(k,3366) = -11568.58586120610000000D0*(rk+10.D0)*(rk+11.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-229.00000000000000000D0+(rk+3.D0)*rk) * den(10)
        ckplm(k,3367) = -11568.58586120610000000D0*(rk-10.D0)*(rk+10.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-231.00000000000000000D0+(rk-1.D0)*rk) * den(11)
        ckplm(k,3368) = -304.43647003173800000D0*(rk-10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-8626.00000000000000000D0+rk&
                *(-189.00000000000000000D0+37.00000000000000000D0*rk)) * den(12)
        ckplm(k,3369) = -103.31554412841800000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(-24310.00000000000000000D0+rk&
                *(-987.00000000000000000D0+103.00000000000000000D0*rk)) * den(13)
        ckplm(k,3370) = -743.87191772460900000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-3135.00000000000000000D0+rk&
                *(-193.00000000000000000D0+13.00000000000000000D0*rk)) * den(14)
        ckplm(k,3371) = -1659.40658569336000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1261.00000000000000000D0+rk&
                *(-109.00000000000000000D0+5.00000000000000000D0*rk)) * den(15)
        ckplm(k,3372) = -6430.20051956177000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-40.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0) * den(16)
        ckplm(k,3373) = -126.08236312866200000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-11628.00000000000000000D0+rk&
                *(-1867.00000000000000000D0+31.00000000000000000D0*rk)) * den(17)
        ckplm(k,3374) = -420.27454376220700000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-2625.00000000000000000D0+(rk-589.D0)*rk) * den(18)
        ckplm(k,3375) = 163.68587493896500000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(4465.00000000000000000D0+rk&
                *(1491.00000000000000000D0+29.00000000000000000D0*rk)) * den(19)
        ckplm(k,3376) = 1063.95818710327000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(177.00000000000000000D0+13.00000000000000000D0*rk) * den(20)
        ckplm(k,3377) = 43622.28567123410000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3378) = -1817.59523630142000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3379) = -132.99477338790900000D0*(rk+10.D0)*(rk+11.D0)*(rk-123.D0)*(rk+12.D0)&
                *(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(1)
        ckplm(k,3380) = 20.46073436737060000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-4056.00000000000000000D0+rk&
                *(841.00000000000000000D0+17.00000000000000000D0*rk)) * den(2)
        ckplm(k,3381) = 0.92165470123291000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(444000.00000000000000000D0+rk*(-192929.00000000000000000D0+rk&
                *(13836.00000000000000000D0+605.00000000000000000D0*rk))) * den(3)
        ckplm(k,3382) = 0.82948923110961900D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(892920.00000000000000000D0+rk*(-265586.00000000000000000D0+rk&
                *(9039.00000000000000000D0+767.00000000000000000D0*rk))) * den(4)
        ckplm(k,3383) = 126.91185235977200000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(8878.00000000000000000D0+rk*(-1897.00000000000000000D0+rk&
                *(14.00000000000000000D0+5.00000000000000000D0*rk))) * den(5)
        ckplm(k,3384) = 10.91714859008790000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(140120.00000000000000000D0+rk*(-21842.00000000000000000D0+rk&
                *(-345.00000000000000000D0+53.00000000000000000D0*rk))) * den(6)
        ckplm(k,3385) = 14.68168258666990000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(130268.00000000000000000D0+rk*(-14637.00000000000000000D0+rk&
                *(-592.00000000000000000D0+33.00000000000000000D0*rk))) * den(7)
        ckplm(k,3386) = 0.67970752716064500D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(3.29832000000000000D6+rk*(-254674.00000000000000000D0+rk&
                *(-18669.00000000000000000D0+535.00000000000000000D0*rk))) * den(8)
        ckplm(k,3387) = 0.66762383778890000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk+9.D0)*(3.73005000000000000D6+rk*(-175501.00000000000000000D0+rk&
                *(-23286.00000000000000000D0+337.00000000000000000D0*rk))) * den(9)
        ckplm(k,3388) = 76.10911750793460000D0*(rk+10.D0)*(rk+11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(34656.00000000000000000D0+rk*(-682.00000000000000000D0+(rk-225.D0)*rk)) * den(10)
        ckplm(k,3389) = -76.10911750793460000D0*(rk-10.D0)*(rk+10.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(-35112.00000000000000000D0+(rk-1.D0)*(rk+229.D0)*rk) * den(11)
        ckplm(k,3390) = -0.66762383778890000D0*(rk-10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(-3.88192800000000000D6+rk*(-127918.00000000000000000D0+rk&
                *(24297.00000000000000000D0+337.00000000000000000D0*rk))) * den(12)
        ckplm(k,3391) = -0.67970752716064500D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(-3.53379000000000000D6+rk*(-215731.00000000000000000D0+rk&
                *(20274.00000000000000000D0+535.00000000000000000D0*rk))) * den(13)
        ckplm(k,3392) = -14.68168258666990000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-144280.00000000000000000D0+rk*(-13354.00000000000000000D0+rk&
                *(691.00000000000000000D0+33.00000000000000000D0*rk))) * den(14)
        ckplm(k,3393) = -10.91714859008790000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-161564.00000000000000000D0+rk*(-20993.00000000000000000D0+rk&
                *(504.00000000000000000D0+53.00000000000000000D0*rk))) * den(15)
        ckplm(k,3394) = -126.91185235977200000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-10784.00000000000000000D0+(5.D0*rk*rk+rk-1910.D0)*rk) * den(16)
        ckplm(k,3395) = -0.82948923110961900D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-1.16677800000000000D6+rk*(-281363.00000000000000000D0+rk&
                *(-6738.00000000000000000D0+767.00000000000000000D0*rk))) * den(17)
        ckplm(k,3396) = -0.92165470123291000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-650160.00000000000000000D0+rk*(-218786.00000000000000000D0+rk&
                *(-12021.00000000000000000D0+605.00000000000000000D0*rk))) * den(18)
        ckplm(k,3397) = -20.46073436737060000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4880.00000000000000000D0+rk&
                *(-807.00000000000000000D0+17.00000000000000000D0*rk)) * den(19)
        ckplm(k,3398) = 132.99477338790900000D0*(rk-10.D0)*(rk-11.D0)*(rk+124.D0)*(rk-12.D0)&
                *(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(20)
        ckplm(k,3399) = 1817.59523630142000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3400) = 72.70380945205690000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3401) = -1.77326364517212000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(656.00000000000000000D0+11.00000000000000000D0*rk) * den(1)
        ckplm(k,3402) = -4.09214687347412000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-2340.00000000000000000D0+rk*(199.00000000000000000D0+9.00000000000000000D0*rk)) &
                * den(2)
        ckplm(k,3403) = -0.03686618804931640D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.78044000000000000D6+rk*(-429158.00000000000000000D0+rk&
                *(4155.00000000000000000D0+983.00000000000000000D0*rk))) * den(3)
        ckplm(k,3404) = -0.01843309402465820D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.42818800000000000D7+rk*(9.60376600000000000D6+rk*(-805943.00000000000000000D0+rk&
                *(-26386.00000000000000000D0+1523.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3405) = -0.18801755905151400D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-4.26708000000000000D6+rk*(1.21542600000000000D6+rk*(-52661.00000000000000000D0+rk&
                *(-4614.00000000000000000D0+89.00000000000000000D0*rk)))) * den(5)
        ckplm(k,3406) = -0.04852066040039060D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.53242720000000000D7+rk*(5.26855000000000000D6+rk*(-55105.00000000000000000D0+rk&
                *(-21910.00000000000000000D0+97.00000000000000000D0*rk)))) * den(6)
        ckplm(k,3407) = 0.10875320434570300D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.54303200000000000D7+rk*(-2.31536600000000000D6+rk*(-49811.00000000000000000D0+rk&
                *(9854.00000000000000000D0+59.00000000000000000D0*rk)))) * den(7)
        ckplm(k,3408) = 0.00302092234293620D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(6.93340560000000000D8+rk*(-7.15178460000000000D7+rk*(-4.31213300000000000D6+rk&
                *(302586.00000000000000000D0+5153.00000000000000000D0*rk)))) * den(8)
        ckplm(k,3409) = 0.66762383778890000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)&
                *(3.63192000000000000D6+rk*(-228326.00000000000000000D0+rk*(-28485.00000000000000000D0+rk&
                *(938.00000000000000000D0+33.00000000000000000D0*rk)))) * den(9)
        ckplm(k,3410) = 25.36970583597820000D0*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(103512.00000000000000000D0+rk*(-2722.00000000000000000D0+rk&
                *(-889.00000000000000000D0+(rk+10.D0)*rk))) * den(10)
        ckplm(k,3411) = 25.36970583597820000D0*(rk-10.D0)*(rk+10.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(105336.00000000000000000D0+(rk-1.D0)*rk*(-918.00000000000000000D0+(rk-5.D0)*rk)) * den(11)
        ckplm(k,3412) = 0.66762383778890000D0*(rk-10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(3.83085600000000000D6+rk*(168674.00000000000000000D0+rk*(-31101.00000000000000000D0+rk&
                *(-806.00000000000000000D0+33.00000000000000000D0*rk)))) * den(12)
        ckplm(k,3413) = 0.00302092234293620D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(7.60248840000000000D8+rk*(6.20064340000000000D7+rk*(-5.18897300000000000D6+rk&
                *(-281974.00000000000000000D0+5153.00000000000000000D0*rk)))) * den(13)
        ckplm(k,3414) = 0.10875320434570300D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(1.76860800000000000D7+rk*(2.18641800000000000D6+rk*(-79019.00000000000000000D0+rk&
                *(-9618.00000000000000000D0+59.00000000000000000D0*rk)))) * den(14)
        ckplm(k,3415) = -0.04852066040039060D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-3.06259200000000000D7+rk*(-5.31264200000000000D6+rk*(11207.00000000000000000D0+rk&
                *(22298.00000000000000000D0+97.00000000000000000D0*rk)))) * den(15)
        ckplm(k,3416) = -0.18801755905151400D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-5.53046400000000000D6+rk*(-1.30655000000000000D6+rk*(-38285.00000000000000000D0+rk&
                *(4970.00000000000000000D0+89.00000000000000000D0*rk)))) * den(16)
        ckplm(k,3417) = -0.01843309402465820D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-3.46636800000000000D7+rk*(-1.11304020000000000D7+rk*(-717647.00000000000000000D0+rk&
                *(32478.00000000000000000D0+1523.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3418) = -0.03686618804931640D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-2.21277000000000000D6+rk*(-434519.00000000000000000D0+rk&
                *(-1206.00000000000000000D0+983.00000000000000000D0*rk))) * den(18)
        ckplm(k,3419) = -4.09214687347412000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-2530.00000000000000000D0+rk*(-181.00000000000000000D0+9.00000000000000000D0*rk)) &
                * den(19)
        ckplm(k,3420) = -1.77326364517212000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-645.00000000000000000D0+11.00000000000000000D0*rk) * den(20)
        ckplm(k,3421) = 72.70380945205690000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3422) = -2.79630036354065000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3423) = 0.06820244789123540D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(1025.00000000000000000D0+29.00000000000000000D0*rk) * den(1)
        ckplm(k,3424) = 2.04607343673706000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-420.00000000000000000D0+(rk+9.D0)*rk) * den(2)
        ckplm(k,3425) = 0.09216547012329100D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(86580.00000000000000000D0+rk*(-10741.00000000000000000D0+rk&
                *(-372.00000000000000000D0+13.00000000000000000D0*rk))) * den(3)
        ckplm(k,3426) = 0.00921654701232910D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-7.35756000000000000D6+rk*(1.77674200000000000D6+rk*(-37967.00000000000000000D0+rk&
                *(-6718.00000000000000000D0+23.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3427) = -0.00552992820739746D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.03346280000000000D8+rk&
                *(3.67254660000000000D7+rk*(-2.91027500000000000D6+rk*(-117535.00000000000000000D0+rk&
                *(10835.00000000000000000D0+109.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3428) = -0.00142707824707031D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-6.91428960000000000D8+rk&
                *(1.79818792000000000D8+rk*(-6.61735000000000000D6+rk*(-891745.00000000000000000D0+rk&
                *(27820.00000000000000000D0+783.00000000000000000D0*rk))))) * den(6)
        ckplm(k,3429) = -0.00024604797363281D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-5.98197384000000000D9+rk&
                *(1.12326948400000000D9+rk*(2.78341500000000000D6+rk*(-6.49212500000000000D6+rk&
                *(40425.00000000000000000D0+5281.00000000000000000D0*rk))))) * den(7)
        ckplm(k,3430) = -0.00003417332967122D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-5.72483376000000000D10+rk&
                *(7.39405716000000000D9+rk*(3.48980278000000000D8+rk*(-4.54204810000000000D7+rk&
                *(-623098.00000000000000000D0+34381.00000000000000000D0*rk))))) * den(8)
        ckplm(k,3431) = -0.00755230585734050D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-3.12580632000000000D8+rk&
                *(2.46135700000000000D7+rk*(2.87415100000000000D6+rk*(-153165.00000000000000000D0+rk&
                *(-6223.00000000000000000D0+107.00000000000000000D0*rk))))) * den(9)
        ckplm(k,3432) = -0.05739752451578780D0*(rk+10.D0)*(rk+11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(-4.55507520000000000D7+5.00000000000000000D0*rk*(300112.00000000000000000D0+rk&
                *(97014.00000000000000000D0+rk*(-1763.00000000000000000D0+(rk-216.D0)*rk)))) * den(10)
        ckplm(k,3433) = 0.05739752451578780D0*(rk-10.D0)*(rk+10.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(4.65585120000000000D7+5.00000000000000000D0*(rk-1.D0)*rk*(-101664.00000000000000000D0+rk&
                *(-667.00000000000000000D0+(rk+222.D0)*rk))) * den(11)
        ckplm(k,3434) = 0.00755230585734050D0*(rk-10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(3.34173216000000000D8+rk&
                *(1.84312000000000000D7+rk*(-3.29523800000000000D6+rk*(-127203.00000000000000000D0+rk&
                *(6758.00000000000000000D0+107.00000000000000000D0*rk))))) * den(12)
        ckplm(k,3435) = 0.00003417332967122D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(6.42486514800000000D10+rk&
                *(6.56249945800000000D9+rk*(-4.81159323000000000D8+rk*(-4.25842790000000000D7+rk&
                *(795003.00000000000000000D0+34381.00000000000000000D0*rk))))) * den(13)
        ckplm(k,3436) = 0.00024604797363281D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(7.09593264000000000D9+rk&
                *(1.09809098400000000D9+rk*(-2.24495300000000000D7+rk*(-6.60101500000000000D6+rk&
                *(-14020.00000000000000000D0+5281.00000000000000000D0*rk))))) * den(14)
        ckplm(k,3437) = 0.00142707824707031D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(8.76946320000000000D8+rk&
                *(1.90270892000000000D8+rk*(3.78302500000000000D6+rk*(-995195.00000000000000000D0+rk&
                *(-23905.00000000000000000D0+783.00000000000000000D0*rk))))) * den(15)
        ckplm(k,3438) = 0.00552992820739746D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.42853760000000000D8+rk&
                *(4.21506160000000000D7+rk*(2.49375000000000000D6+rk*(-159785.00000000000000000D0+rk&
                *(-10290.00000000000000000D0+109.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3439) = -0.00921654701232910D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-9.16552800000000000D6+rk*(-1.83243000000000000D6+rk*(-17675.00000000000000000D0+rk&
                *(6810.00000000000000000D0+23.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3440) = -0.09216547012329100D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-96936.00000000000000000D0+rk*(-9958.00000000000000000D0+rk&
                *(411.00000000000000000D0+13.00000000000000000D0*rk))) * den(18)
        ckplm(k,3441) = -2.04607343673706000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-428.00000000000000000D0+(rk-7.D0)*rk) * den(19)
        ckplm(k,3442) = -0.06820244789123540D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-996.00000000000000000D0+29.00000000000000000D0*rk) * den(20)
        ckplm(k,3443) = 2.79630036354065000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3444) = 0.10356668013113500D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3445) = -0.00757804976569282D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(492.00000000000000000D0+17.00000000000000000D0*rk) * den(1)
        ckplm(k,3446) = -0.07578049765692820D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-849.00000000000000000D0+(rk-13.D0)*rk) * den(2)
        ckplm(k,3447) = 0.00068270718609845D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.15329000000000000D6+rk*(56279.00000000000000000D0+(rk+5325.D0)*rk)) * den(3)
        ckplm(k,3448) = 0.00307218233744303D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-383940.00000000000000000D0+rk*(-3180.00000000000000000D0+rk&
                *(1193.00000000000000000D0+17.00000000000000000D0*rk))) * den(4)
        ckplm(k,3449) = 0.00061443646748861D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.32660720000000000D8+rk&
                *(2.98358520000000000D7+rk*(-511540.00000000000000000D0+rk*(-150455.00000000000000000D0+rk&
                *(1690.00000000000000000D0+113.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3450) = 0.00005285474989149D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1.49912893800000000D10+rk*(-4.67443822200000000D9+rk&
                *(2.95193863000000000D8+rk*(2.25115950000000000D7+rk*(-1.73087000000000000D6+rk&
                *(-30813.00000000000000000D0+1087.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3451) = 0.00002733866373698D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(4.72080565800000000D10+rk*(-1.06441584780000000D10+rk&
                *(1.77036211000000000D8+rk*(7.25558730000000000D7+rk*(-1.81455200000000000D6+rk&
                *(-118659.00000000000000000D0+997.00000000000000000D0*rk)))))) * den(7)
        ckplm(k,3452) = -0.00002050399780273D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-8.90998884000000000D10+rk*(1.38300465600000000D10+rk&
                *(4.70535778000000000D8+rk*(-1.08328583000000000D8+rk*(-568569.00000000000000000D0+rk&
                *(179663.00000000000000000D0+431.00000000000000000D0*rk)))))) * den(8)
        ckplm(k,3453) = -0.00050348705715603D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-4.56435720000000000D9+rk*(4.32140760000000000D8+rk&
                *(4.69867420000000000D7+rk*(-3.57192300000000000D6+rk*(-141101.00000000000000000D0+rk&
                *(5763.00000000000000000D0+79.00000000000000000D0*rk)))))) * den(9)
        ckplm(k,3454) = -0.05739752451578780D0*(rk+10.D0)*(rk+11.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-4.53492000000000000D7+rk*(1.79663400000000000D6+rk&
                *(574735.00000000000000000D0+rk*(-14511.00000000000000000D0+rk&
                *(-1904.00000000000000000D0+(rk+21.D0)*rk))))) * den(10)
        ckplm(k,3455) = -0.05739752451578780D0*(rk-10.D0)*(rk+10.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-4.65585120000000000D7+(rk-1.D0)*rk&
                *(611346.00000000000000000D0+rk*(4697.00000000000000000D0+rk&
                *(-2008.00000000000000000D0+(rk-14.D0)*rk)))) * den(11)
        ckplm(k,3456) = -0.00050348705715603D0*(rk-10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-4.94608608000000000D9+rk*(-3.28044252000000000D8+rk&
                *(5.67994600000000000D7+rk*(2.95146900000000000D6+rk*(-168731.00000000000000000D0+rk&
                *(-5289.00000000000000000D0+79.00000000000000000D0*rk)))))) * den(12)
        ckplm(k,3457) = -0.00002050399780273D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(-1.02351818400000000D11+rk&
                *(-1.25671592600000000D10+rk*(7.90319948000000000D8+rk*(1.04266297000000000D8+rk&
                *(-1.46041900000000000D6+rk*(-177077.00000000000000000D0+431.00000000000000000D0*rk)))))) &
                * den(13)
        ckplm(k,3458) = 0.00002733866373698D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(5.79550005000000000D10+rk*(1.07739043500000000D10+rk&
                *(-5.03171750000000000D7+rk*(-7.86075510000000000D7+rk*(-1.20630200000000000D6+rk&
                *(124641.00000000000000000D0+997.00000000000000000D0*rk)))))) * den(14)
        ckplm(k,3459) = 0.00005285474989149D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.99367109000000000D10+rk*(5.19052827000000000D9+rk&
                *(2.17598293000000000D8+rk*(-2.91052050000000000D7+rk*(-1.56050000000000000D6+rk&
                *(37335.00000000000000000D0+1087.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3460) = 0.00061443646748861D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.62856080000000000D8+rk&
                *(3.04013720000000000D7+rk*(51165.00000000000000000D0+rk*(-156085.00000000000000000D0+rk&
                *(-1125.00000000000000000D0+113.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3461) = 0.00307218233744303D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(379584.00000000000000000D0+rk*(-5515.00000000000000000D0+rk&
                *(-1142.00000000000000000D0+17.00000000000000000D0*rk))) * den(17)
        ckplm(k,3462) = 0.00068270718609845D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(1.20424500000000000D6+rk*(45632.00000000000000000D0+(rk-5322.D0)*rk)) * den(18)
        ckplm(k,3463) = -0.07578049765692820D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-835.00000000000000000D0+(rk+15.D0)*rk) * den(19)
        ckplm(k,3464) = -0.00757804976569282D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-475.00000000000000000D0+17.00000000000000000D0*rk) * den(20)
        ckplm(k,3465) = 0.10356668013113500D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3466) = -0.00369881000468340D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+8.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,3467) = 0.00063150414714107D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+8.D0)*(rk+9.D0)&
                *(287.00000000000000000D0+11.00000000000000000D0*rk) * den(1)
        ckplm(k,3468) = 0.00126300829428214D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-3318.00000000000000000D0+(rk-127.D0)*rk) * den(2)
        ckplm(k,3469) = -0.00017067679652461D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-385910.00000000000000000D0+rk*(-743.00000000000000000D0+rk&
                *(1244.00000000000000000D0+17.00000000000000000D0*rk))) * den(3)
        ckplm(k,3470) = -0.00005120303895738D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+8.D0)*(rk+9.D0)*(1.64631600000000000D7+rk&
                *(-946978.00000000000000000D0+rk*(-101677.00000000000000000D0+rk&
                *(1102.00000000000000000D0+73.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3471) = -0.00005120303895738D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0)*(-1.89557928000000000D8+rk*(2.39748420000000000D7+rk&
                *(1.06727500000000000D6+rk*(-108325.00000000000000000D0+rk&
                *(-2467.00000000000000000D0+43.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3472) = 0.00001321368747287D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(-8.02768932000000000D9+rk*(1.57712555600000000D9+rk&
                *(2.14850400000000000D6+rk*(-1.04970610000000000D7+rk*(75489.00000000000000000D0+rk&
                *(15833.00000000000000000D0+15.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3473) = 6.50920565166171000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+8.D0)*(rk+9.D0)*(1.73826238345200000D13+rk*(-4.57318224901200000D12+rk&
                *(1.63668451912000000D11+rk*(3.36425044690000000D10+rk*(-1.66613573000000000D9+rk&
                *(-7.37195480000000000D7+rk*(2.40557800000000000D6+33571.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3474) = 5.69555494520399000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(2.99471800320000000D12+rk*(-5.42994205680000000D11+rk&
                *(-1.12827339240000000D10+rk*(5.04683435200000000D9+rk*(-2.92725510000000000D7+rk&
                *(-1.27562330000000000D7+rk*(56595.00000000000000000D0+5281.00000000000000000D0*rk))))))) &
                * den(8)
        ckplm(k,3475) = 0.00002517435285780D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk+9.D0)*(8.88565860000000000D10+rk*(-9.83334126000000000D9+rk&
                *(-9.89572276000000000D8+rk*(1.00635591000000000D8+rk*(3.59618000000000000D6+rk&
                *(-260190.00000000000000000D0+rk*(-4144.00000000000000000D0+99.00000000000000000D0*rk))))))) &
                * den(9)
        ckplm(k,3476) = 0.00095662540859646D0*(rk+10.D0)*(rk+11.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(2.70885888000000000D9+rk*(-1.25479992000000000D8+rk&
                *(-3.97151720000000000D7+rk*(1.28997400000000000D6+rk*(174055.00000000000000000D0+rk&
                *(-3143.00000000000000000D0+(rk-203.D0)*rk)))))) * den(10)
        ckplm(k,3477) = -0.00095662540859646D0*(rk-10.D0)*(rk+10.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(-2.79351072000000000D9+(rk-1.D0)*rk*(4.28904360000000000D7+rk&
                *(378036.00000000000000000D0+rk*(-188383.00000000000000000D0+rk&
                *(-1693.00000000000000000D0+(rk+211.D0)*rk))))) * den(11)
        ckplm(k,3478) = -0.00002517435285780D0*(rk-10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(-9.76035715200000000D10+rk*(-7.56795004800000000D9+rk&
                *(1.26736430800000000D9+rk*(8.37353160000000000D7+rk*(-4.83150500000000000D6+rk&
                *(-233247.00000000000000000D0+rk*(4837.00000000000000000D0+99.00000000000000000D0*rk))))))) &
                * den(12)
        ckplm(k,3479) = -5.69555494520399000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(-3.52136617560000000D12+rk*(-5.05235228340000000D11+rk&
                *(2.64705719320000000D10+rk*(5.03541516100000000D9+rk*(-3.51727040000000000D7+rk&
                *(-1.29849020000000000D7+rk*(-19628.00000000000000000D0+5281.00000000000000000D0*rk))))))) &
                * den(13)
        ckplm(k,3480) = -6.50920565166171000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-2.20842419868000000D13+rk*(-4.79330989272000000D12+rk&
                *(-5.35166982840000000D10+rk*(3.95229153340000000D10+rk*(1.26262930500000000D9+rk&
                *(-8.74480250000000000D7+rk*(-2.17058100000000000D6+33571.00000000000000000D0*rk))))))) &
                * den(14)
        ckplm(k,3481) = -0.00001321368747287D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-9.59210964000000000D9+rk*(-1.54111448400000000D9+rk&
                *(3.39345160000000000D7+rk*(1.06409870000000000D7+rk*(-3451.00000000000000000D0+rk&
                *(-15743.00000000000000000D0+15.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3482) = 0.00005120303895738D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(2.12359680000000000D8+rk*(2.15254000000000000D7+rk&
                *(-1.37701800000000000D6+rk*(-98027.00000000000000000D0+rk&
                *(2682.00000000000000000D0+43.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3483) = 0.00005120303895738D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.73074320000000000D7+rk&
                *(740610.00000000000000000D0+rk*(-104545.00000000000000000D0+rk&
                *(-810.00000000000000000D0+73.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3484) = 0.00017067679652461D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(383940.00000000000000000D0+rk&
                *(-3180.00000000000000000D0+rk*(-1193.00000000000000000D0+17.00000000000000000D0*rk))) * den(18)
        ckplm(k,3485) = -0.00126300829428214D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-3190.00000000000000000D0+(rk+129.D0)*rk) * den(19)
        ckplm(k,3486) = -0.00063150414714107D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-276.00000000000000000D0+11.00000000000000000D0*rk) * den(20)
        ckplm(k,3487) = 0.00369881000468340D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(21)
        ckplm(k,3488) = 0.00012754517257529D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+9.D0) * den(0)
        ckplm(k,3489) = -3.11085786768999000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+9.D0)&
                *(2624.00000000000000000D0+107.00000000000000000D0*rk) * den(1)
        ckplm(k,3490) = 6.22171573537998000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+9.D0)*(39216.00000000000000000D0+rk&
                *(2099.00000000000000000D0+13.00000000000000000D0*rk)) * den(2)
        ckplm(k,3491) = 8.40772396672970000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+9.D0)*(-5.72464000000000000D6+rk&
                *(-214390.00000000000000000D0+rk*(8713.00000000000000000D0+253.00000000000000000D0*rk))) &
                * den(3)
        ckplm(k,3492) = 2.52231719001891000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+9.D0)*(2.94648480000000000D8+rk&
                *(-124714.00000000000000000D0+rk*(-1.47693100000000000D6+rk&
                *(-25034.00000000000000000D0+439.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3493) = -2.52231719001891000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+9.D0)*(3.93720998400000000D9+rk&
                *(-1.98021480000000000D8+rk*(-3.08046500000000000D7+rk*(412775.00000000000000000D0+rk&
                *(48446.00000000000000000D0+205.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3494) = -6.50920565166171000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+9.D0)*(-1.87095926880000000D11+rk*(1.95255763600000000D10+rk&
                *(1.54885053800000000D9+rk*(-1.16304773000000000D8+rk*(-4.61146300000000000D6+rk&
                *(110197.00000000000000000D0+2309.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3495) = -6.50920565166171000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+9.D0)*(2.17684667808000000D12+rk*(-3.43310152848000000D11+rk&
                *(-1.30835347720000000D10+rk*(2.96787293600000000D9+rk*(2.95405250000000000D7+rk&
                *(-6.75255700000000000D6+rk*(-48913.00000000000000000D0+2189.00000000000000000D0*rk))))))) &
                * den(7)
        ckplm(k,3496) = -8.13650706457714000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+9.D0)&
                *(-1.95681869712000000D13+rk*(4.05906674328000000D12+rk*(3.03388875240000000D10+rk&
                *(-4.27083412840000000D10+rk*(7.99073107000000000D8+rk*(1.40686280000000000D8+rk&
                *(-2.27263400000000000D6+rk*(-132436.00000000000000000D0+643.00000000000000000D0*rk)))))))) &
                * den(8)
        ckplm(k,3497) = 3.59633612254309000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-8.D0)*(rk+9.D0)&
                *(6.05372040000000000D11+rk*(-7.67027762400000000D10+rk*(-7.09383768400000000D9+rk&
                *(9.30071156000000000D8+rk*(2.84286730000000000D7+rk*(-3.31492000000000000D6+rk&
                *(-46046.00000000000000000D0+rk*(3044.00000000000000000D0+17.00000000000000000D0*rk)))))))) &
                * den(9)
        ckplm(k,3498) = 0.00013666077265664D0*(rk+10.D0)*(rk+11.D0)*(rk-8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.88773603200000000D10+rk*(-1.00155076800000000D9+rk*(-3.13576180000000000D8+rk&
                *(1.25043240000000000D7+rk*(1.70412900000000000D6+rk*(-43512.00000000000000000D0+rk&
                *(-3150.00000000000000000D0+(rk+36.D0)*rk))))))) * den(10)
        ckplm(k,3499) = 0.00013666077265664D0*(rk-10.D0)*(rk+10.D0)*(rk-8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.95545750400000000D10+(rk-1.D0)*rk*(-3.43900368000000000D8+rk*(-3.42313200000000000D6+rk&
                *(1.89376000000000000D6+rk*(20511.00000000000000000D0+rk*(-3401.00000000000000000D0+(rk-27.D0)&
                *rk)))))) * den(11)
        ckplm(k,3500) = 3.59633612254309000D-6*(rk-10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(6.74082601920000000D11+rk*(5.98548792480000000D10+rk*(-9.68108405200000000D9+rk&
                *(-7.84233772000000000D8+rk*(4.42072330000000000D7+rk*(2.97567200000000000D6+rk&
                *(-66878.00000000000000000D0+rk*(-2908.00000000000000000D0+17.00000000000000000D0*rk)))))))) &
                * den(12)
        ckplm(k,3501) = -8.13650706457714000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-2.35535502384000000D13+rk*(-3.86778378696000000D12+rk*(1.61820196868000000D11+rk&
                *(4.44569895000000000D10+rk*(6.62324670000000000D7+rk*(-1.51504920000000000D8+rk&
                *(-1.32757800000000000D6+rk*(137580.00000000000000000D0+643.00000000000000000D0*rk)))))))) &
                * den(13)
        ckplm(k,3502) = -6.50920565166171000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-2.50414166520000000D12+rk*(-3.08391080580000000D11+rk*(2.17431645240000000D10+rk&
                *(2.78324014100000000D9+rk*(-6.24930000000000000D7+rk*(-6.41311000000000000D6+rk&
                *(64236.00000000000000000D0+2189.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3503) = -6.50920565166171000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-8.D0)*(rk-9.D0)*(-2.04961067280000000D11+rk*(-1.60979439480000000D10+rk&
                *(1.86902874400000000D9+rk*(9.68031310000000000D7+rk*(-5.12781300000000000D6+rk&
                *(-96343.00000000000000000D0+2309.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3504) = -2.52231719001891000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-8.D0)*(rk-9.D0)*(-4.10406228000000000D9+rk&
                *(-1.35366614000000000D8+rk*(3.17543490000000000D7+rk*(221041.00000000000000000D0+rk&
                *(-47421.00000000000000000D0+205.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3505) = 2.52231719001891000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-8.D0)*(rk-9.D0)*(2.93321736000000000D8+rk*(-2.75229000000000000D6+rk&
                *(-1.39919500000000000D6+rk*(26790.00000000000000000D0+439.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3506) = 8.40772396672970000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-8.D0)*(rk-9.D0)*(5.50179000000000000D6+rk&
                *(-231057.00000000000000000D0+rk*(-7954.00000000000000000D0+253.00000000000000000D0*rk))) &
                * den(18)
        ckplm(k,3507) = 6.22171573537998000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-8.D0)*(rk-9.D0)*(37130.00000000000000000D0+rk&
                *(-2073.00000000000000000D0+13.00000000000000000D0*rk)) * den(19)
        ckplm(k,3508) = -3.11085786768999000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-2517.00000000000000000D0+107.00000000000000000D0*rk) * den(20)
        ckplm(k,3509) = 0.00012754517257529D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-8.D0)*(rk-9.D0) * den(21)
        ckplm(k,3510) = -4.25150575250965000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) &
                * den(0)
        ckplm(k,3511) = 3.11085786768999000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(1107.00000000000000000D0+47.00000000000000000D0*rk) * den(1)
        ckplm(k,3512) = -2.39296759053076000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(54054.00000000000000000D0+rk*(3467.00000000000000000D0+43.00000000000000000D0*rk)) * den(2)
        ckplm(k,3513) = -2.15582665813582000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(-1.45604250000000000D8+rk&
                *(-9.19550900000000000D6+rk*(-9420.00000000000000000D0+4139.00000000000000000D0*rk))) * den(3)
        ckplm(k,3514) = 9.70121996161120000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(-5.96022840000000000D8+rk*(-2.55387900000000000D7+rk&
                *(1.62940900000000000D6+rk*(72474.00000000000000000D0+227.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3515) = 5.82073197696672000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(1.53216457200000000D10+rk*(1.50291234000000000D8+rk&
                *(-9.87999050000000000D7+rk*(-2.34518500000000000D6+rk&
                *(77705.00000000000000000D0+1471.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3516) = 5.00708127050901000D-10*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(-2.44534890060000000D12+rk*(7.19813160360000000D10+rk&
                *(2.35853651360000000D10+rk*(-8.80274850000000000D7+rk*(-5.43376150000000000D7+rk&
                *(-503391.00000000000000000D0+12839.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3517) = -9.76380847749256000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(-1.59022300680000000D11+rk*(1.11305265480000000D10+rk*(1.83890452800000000D9+rk&
                *(-7.80528850000000000D7+rk*(-6.97359000000000000D6+rk*(99812.00000000000000000D0+rk&
                *(6822.00000000000000000D0+5.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3518) = -7.32285635811942000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(2.53655922240000000D12+rk*(-2.75328570960000000D11+rk*(-3.09321869880000000D10+rk&
                *(2.95150450000000000D9+rk*(1.41792343000000000D8+rk*(-8.95020000000000000D6+rk&
                *(-259042.00000000000000000D0+rk*(6260.00000000000000000D0+87.00000000000000000D0*rk)))))))) &
                * den(8)
        ckplm(k,3519) = -1.38320620097811000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)&
                *(-1.53174385728000000D13+rk*(2.18714613768000000D12+rk*(1.84338274812000000D11+rk&
                *(-3.04836989840000000D10+rk*(-7.47057675000000000D8+rk*(1.37972877000000000D8+rk&
                *(1.32757800000000000D6+rk*(-209706.00000000000000000D0+rk&
                *(-1035.00000000000000000D0+53.00000000000000000D0*rk))))))))) * den(9)
        ckplm(k,3520) = -3.15371013823010000D-6*(rk+10.D0)*(rk+11.D0)*(rk-9.D0)&
                *(-8.14350700800000000D11+rk*(4.87133049600000000D10+rk*(1.50840210960000000D10+rk&
                *(-7.15624876000000000D8+rk*(-9.75520140000000000D7+rk*(3.25938900000000000D6+rk&
                *(248304.00000000000000000D0+rk*(-4674.00000000000000000D0+(rk-186.D0)*rk)))))))) * den(10)
        ckplm(k,3521) = 3.15371013823010000D-6*(rk-10.D0)*(rk+10.D0)*(rk-9.D0)&
                *(8.47364918400000000D11+(rk-1.D0)*rk*(-1.68033720960000000D10+rk*(-1.86564876000000000D8+rk&
                *(1.11377224000000000D8+rk*(1.40326900000000000D6+rk*(-278684.00000000000000000D0+rk&
                *(-2954.00000000000000000D0+(rk+196.D0)*rk))))))) * den(11)
        ckplm(k,3522) = 1.38320620097811000D-7*(rk-10.D0)*(rk-11.D0)*(rk-9.D0)&
                *(1.72906462310400000D13+rk*(1.73068716153600000D12+rk*(-2.69951583552000000D11+rk&
                *(-2.61495683720000000D10+rk*(1.40974780800000000D9+rk*(1.25668221000000000D8+rk&
                *(-2.76208800000000000D6+rk*(-199518.00000000000000000D0+rk&
                *(1512.00000000000000000D0+53.00000000000000000D0*rk))))))))) * den(12)
        ckplm(k,3523) = 7.32285635811942000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(2.77815457920000000D12+rk*(2.05220006480000000D11+rk*(-3.88504590840000000D10+rk&
                *(-2.30022819600000000D9+rk*(1.82444703000000000D8+rk*(7.26936000000000000D6+rk&
                *(-300426.00000000000000000D0+rk*(-5564.00000000000000000D0+87.00000000000000000D0*rk)))))))) &
                * den(13)
        ckplm(k,3524) = 9.76380847749256000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(1.68242936400000000D11+rk*(7.24691136000000000D9+rk*(-2.03032574800000000D9+rk&
                *(-4.92966700000000000D7+rk*(7.37049500000000000D6+rk*(58985.00000000000000000D0+rk&
                *(-6787.00000000000000000D0+5.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3525) = -5.00708127050901000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-9.D0)*(-2.49371064540000000D12+rk*(-2.47612597800000000D10+rk&
                *(2.35286483960000000D10+rk*(-1.24032285000000000D8+rk*(-5.16280750000000000D7+rk&
                *(580425.00000000000000000D0+12839.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3526) = -5.82073197696672000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-9.D0)*(-1.50749760000000000D10+rk*(3.40552024000000000D8+rk&
                *(9.13128300000000000D7+rk*(-2.64129500000000000D6+rk&
                *(-70350.00000000000000000D0+1471.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3527) = -9.70121996161120000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-9.D0)*(-5.68926888000000000D8+rk&
                *(2.85810940000000000D7+rk*(1.41334900000000000D6+rk&
                *(-71566.00000000000000000D0+227.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3528) = 2.15582665813582000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-9.D0)*(1.36422300000000000D8+rk&
                *(-9.16425200000000000D6+rk*(21837.00000000000000000D0+4139.00000000000000000D0*rk))) * den(18)
        ckplm(k,3529) = 2.39296759053076000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-9.D0)*(50630.00000000000000000D0+rk&
                *(-3381.00000000000000000D0+43.00000000000000000D0*rk)) * den(19)
        ckplm(k,3530) = -3.11085786768999000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-9.D0)&
                *(-1060.00000000000000000D0+47.00000000000000000D0*rk) * den(20)
        ckplm(k,3531) = 4.25150575250965000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-9.D0) * den(21)
        ckplm(k,3532) = 1.37145346855150000D-7*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3533) = -3.34500845988171000D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(4100.00000000000000000D0+179.00000000000000000D0*rk) * den(1)
        ckplm(k,3534) = 2.57308343067824000D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(245895.00000000000000000D0+rk&
                *(17659.00000000000000000D0+281.00000000000000000D0*rk)) * den(2)
        ckplm(k,3535) = 3.47713977118681000D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(-5.33884100000000000D7+rk&
                *(-4.39264100000000000D6+rk*(-76459.00000000000000000D0+425.00000000000000000D0*rk))) * den(3)
        ckplm(k,3536) = -3.12942579406813000D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(-1.28961518000000000D9+rk*(-9.83707680000000000D7+rk&
                *(228463.00000000000000000D0+rk*(122102.00000000000000000D0+1443.00000000000000000D0*rk)))) &
                * den(4)
        ckplm(k,3537) = -2.08628386271209000D-11*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(3.42462364560000000D11+rk*(1.98477719160000000D10+rk&
                *(-1.16965330000000000D9+rk*(-7.94847350000000000D7+rk&
                *(-379190.00000000000000000D0+17609.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3538) = 1.66902709016967000D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(6.56295254100000000D11+rk*(2.12596783900000000D10+rk*(-4.81321050700000000D9+rk&
                *(-1.96316975000000000D8+rk*(5.37467000000000000D6+rk&
                *(237145.00000000000000000D0+557.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3539) = 1.08486760861028000D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(-1.39415380200000000D11+rk*(-5.69118000000000000D8+rk*(1.55031299200000000D9+rk&
                *(2.74332130000000000D7+rk*(-4.38369500000000000D6+rk*(-103190.00000000000000000D0+rk&
                *(2143.00000000000000000D0+37.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3540) = 1.35608451076286000D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(1.42024654800000000D12+rk*(-3.34086084000000000D10+rk*(-2.05054027800000000D10+rk&
                *(1.40334208000000000D8+rk*(9.33860270000000000D7+rk*(490420.00000000000000000D0+rk&
                *(-129430.00000000000000000D0+rk*(-1268.00000000000000000D0+23.00000000000000000D0*rk)))))))) &
                * den(8)
        ckplm(k,3541) = -4.61068733659371000D-8*(rk+11.D0)*(rk+12.D0)*(4.96835015040000000D12+rk&
                *(-2.37504095520000000D11+rk*(-8.64176695040000000D10+rk*(2.76823396400000000D9+rk&
                *(5.27766876000000000D8+rk*(-8.87273100000000000D6+rk*(-1.24899600000000000D6+rk&
                *(5886.00000000000000000D0+(rk+824.D0)*rk)))))))) * den(9)
        ckplm(k,3542) = -3.50412237581122000D-7*(rk+11.D0)*(-7.29614208960000000D12+rk&
                *(4.86002990880000000D11+rk*(1.48804542876000000D11+rk*(-8.21193340000000000D9+rk&
                *(-1.11296564500000000D9+rk*(4.63998150000000000D7+rk*(3.61122300000000000D6+rk&
                *(-97350.00000000000000000D0+rk*(-4455.00000000000000000D0+(rk+55.D0)*rk))))))))) * den(10)
        ckplm(k,3543) = -3.50412237581122000D-7*(rk-10.D0)*(-7.62628426560000000D12+(rk-1.D0)*rk&
                *(1.68419653920000000D11+rk*(2.06501684400000000D9+rk*(-1.30643355600000000D9+rk&
                *(-1.87258610000000000D7+rk*(4.21842400000000000D6+rk*(54901.00000000000000000D0+rk&
                *(-4949.00000000000000000D0+(rk-44.D0)*rk)))))))) * den(11)
        ckplm(k,3544) = -4.61068733659371000D-8*(rk-10.D0)*(rk-11.D0)*(-5.11720372800000000D12+rk&
                *(-5.85119571840000000D10+rk*(9.14858783400000000D10+rk*(5.93579020000000000D8+rk&
                *(-5.53247135000000000D8+rk*(-1.30116700000000000D6+rk*(1.26721000000000000D6+rk&
                *(-670.00000000000000000D0+(rk-815.D0)*rk)))))))) * den(12)
        ckplm(k,3545) = 1.35608451076286000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)&
                *(1.43310218688000000D12+rk*(-7.65287529600000000D9+rk*(-2.03729076200000000D10+rk&
                *(2.25762768000000000D8+rk*(8.90384670000000000D7+rk*(-1.23908400000000000D6+rk&
                *(-119910.00000000000000000D0+rk*(1452.00000000000000000D0+23.00000000000000000D0*rk)))))))) &
                * den(13)
        ckplm(k,3546) = 1.08486760861028000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(1.37327660820000000D11+rk*(-3.57043811400000000D9+rk*(-1.44277445100000000D9+rk&
                *(4.38945280000000000D7+rk*(3.83689500000000000D6+rk*(-115271.00000000000000000D0+rk&
                *(-1884.00000000000000000D0+37.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3547) = 1.66902709016967000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(6.30423820260000000D11+rk*(-3.02768321820000000D10+rk*(-4.19437465700000000D9+rk&
                *(2.15455345000000000D8+rk*(4.19730000000000000D6+rk&
                *(-233803.00000000000000000D0+557.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3548) = -2.08628386271209000D-11*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(-3.21524027280000000D11+rk*(2.19502291160000000D10+rk&
                *(9.33650325000000000D8+rk*(-7.77918850000000000D7+rk&
                *(467235.00000000000000000D0+17609.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3549) = -3.12942579406813000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(-1.19113660800000000D9+rk*(9.84671600000000000D7+rk&
                *(-129185.00000000000000000D0+rk*(-116330.00000000000000000D0+1443.00000000000000000D0*rk)))) &
                * den(17)
        ckplm(k,3550) = 3.47713977118681000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(4.90726530000000000D7+rk&
                *(-4.23844800000000000D6+rk*(77734.00000000000000000D0+425.00000000000000000D0*rk))) * den(18)
        ckplm(k,3551) = 2.57308343067824000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(228517.00000000000000000D0+rk&
                *(-17097.00000000000000000D0+281.00000000000000000D0*rk)) * den(19)
        ckplm(k,3552) = 3.34500845988171000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(3921.00000000000000000D0-179.00000000000000000D0*rk) * den(20)
        ckplm(k,3553) = 1.37145346855150000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3554) = -4.28579208922344000D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3555) = 1.04531514371303000D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(4961.00000000000000000D0+221.00000000000000000D0*rk) * den(1)
        ckplm(k,3556) = -8.04088572086950000D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(358644.00000000000000000D0+rk&
                *(27809.00000000000000000D0+505.00000000000000000D0*rk)) * den(2)
        ckplm(k,3557) = 2.17321235699176000D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(4.63727660000000000D8+rk*(4.49115710000000000D7+rk&
                *(1.19566000000000000D6+5869.00000000000000000D0*rk))) * den(3)
        ckplm(k,3558) = 9.77945560646290000D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(-2.61760268000000000D9+rk*(-2.68427538000000000D8+rk&
                *(-6.03065100000000000D6+rk*(95178.00000000000000000D0+2771.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3559) = -5.92694279179570000D-14*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(-8.72674989516000000D12+rk*(-8.46683992266000000D11+rk*(-2.42576004500000000D9+rk&
                *(1.62831225500000000D9+rk*(3.66452450000000000D7+23011.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3560) = -5.21570965678021000D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(1.69312046112000000D12+rk*(1.42220749112000000D11+rk*(-5.17839163800000000D9+rk&
                *(-6.56420995000000000D8+rk*(-7.89337500000000000D6+rk&
                *(321203.00000000000000000D0+4413.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3561) = -3.39021127690714000D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(-3.91196139840000000D11+rk*(-2.63457621840000000D10+rk*(2.75640809800000000D9+rk&
                *(2.16297547000000000D8+rk*(-2.16302000000000000D6+rk*(-340526.00000000000000000D0+rk&
                *(-3158.00000000000000000D0+43.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3562) = 8.47552819226785000D-11*(rk+12.D0)*(rk+13.D0)*(-2.12032453248000000D13+rk&
                *(-1.06925245704000000D12+rk*(2.37632924364000000D11+rk*(1.33106107400000000D10+rk&
                *(-6.62444699000000000D8+rk*(-4.18028800000000000D7+rk*(240226.00000000000000000D0+rk&
                *(26780.00000000000000000D0+109.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,3563) = 1.44083979268553000D-9*(rk+12.D0)*(1.54691881344000000D13+rk&
                *(5.57389784400000000D11+rk*(-2.37633562108000000D11+rk*(-9.47189572800000000D9+rk&
                *(1.13640263100000000D9+rk*(4.74330150000000000D7+rk*(-1.73632200000000000D6+rk&
                *(-71142.00000000000000000D0+rk*(439.00000000000000000D0+15.00000000000000000D0*rk))))))))) &
                * den(9)
        ckplm(k,3564) = 1.09503824244101000D-8*(-2.32420091904000000D13+rk&
                *(-6.17473690560000000D11+rk*(4.54852347456000000D11+rk*(1.28802906200000000D10+rk&
                *(-3.08724274000000000D9+rk*(-8.82826350000000000D7+rk*(8.40699300000000000D6+rk&
                *(223530.00000000000000000D0+rk*(-7710.00000000000000000D0+(rk-155.D0)*rk))))))))) * den(10)
        ckplm(k,3565) = -1.09503824244101000D-8*(-2.21855542272000000D13+rk&
                *(1.47667877280000000D12+rk*(3.98692046016000000D11+rk*(-2.41865375400000000D10+rk&
                *(-2.52806818000000000D9+rk*(1.33618485000000000D8+rk*(6.63963300000000000D6+rk&
                *(-279510.00000000000000000D0+rk*(-6270.00000000000000000D0+(rk+165.D0)*rk))))))))) * den(11)
        ckplm(k,3566) = -1.44083979268553000D-9*(rk-11.D0)*(-1.46847239884800000D13+rk&
                *(9.99942692544000000D11+rk*(2.02898328384000000D11+rk*(-1.35109629560000000D10+rk&
                *(-8.75711536000000000D8+rk*(5.63342710000000000D7+rk*(1.22729600000000000D6+rk&
                *(-74114.00000000000000000D0+rk*(-304.00000000000000000D0+15.00000000000000000D0*rk))))))))) &
                * den(12)
        ckplm(k,3567) = -8.47552819226785000D-11*(rk-11.D0)*(rk-12.D0)*(-1.99102909824000000D13+rk&
                *(1.50214696392000000D12+rk*(1.94147496812000000D11+rk*(-1.55384874120000000D10+rk&
                *(-4.50756579000000000D8+rk*(4.26879600000000000D7+rk*(55818.00000000000000000D0+rk&
                *(-25908.00000000000000000D0+109.00000000000000000D0*rk)))))))) * den(13)
        ckplm(k,3568) = 3.39021127690714000D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(3.62312092800000000D11+rk*(-3.12027170400000000D10+rk*(-2.09789432400000000D9+rk&
                *(2.21609032000000000D8+rk*(509265.00000000000000000D0+rk*(-320675.00000000000000000D0+rk&
                *(3459.00000000000000000D0+43.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3569) = 5.21570965678021000D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(1.54636953120000000D12+rk*(-1.50641422440000000D11+rk*(-3.25963473800000000D9+rk&
                *(6.21723725000000000D8+rk*(-9.43319500000000000D6+rk&
                *(-294725.00000000000000000D0+4413.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3570) = 5.92694279179570000D-14*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(7.88408335296000000D12+rk*(-8.37094001336000000D11+rk*(7.09105545000000000D9+rk&
                *(1.48196138500000000D9+rk*(-3.65301900000000000D7+23011.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3571) = -9.77945560646290000D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(-2.35529820000000000D9+rk*(2.56091786000000000D8+rk&
                *(-6.29955900000000000D6+rk*(-84094.00000000000000000D0+2771.00000000000000000D0*rk)))) &
                * den(17)
        ckplm(k,3572) = -2.17321235699176000D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(-4.20005880000000000D8+rk*(4.25378580000000000D7+rk&
                *(-1.17805300000000000D6+5869.00000000000000000D0*rk))) * den(18)
        ckplm(k,3573) = 8.04088572086950000D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(331340.00000000000000000D0+rk&
                *(-26799.00000000000000000D0+505.00000000000000000D0*rk)) * den(19)
        ckplm(k,3574) = 1.04531514371303000D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(4740.00000000000000000D0-221.00000000000000000D0*rk) * den(20)
        ckplm(k,3575) = 4.28579208922344000D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3576) = 1.29872487552226000D-10*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3577) = -9.50286494284577000D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(1968.00000000000000000D0+89.00000000000000000D0&
                *rk) * den(1)
        ckplm(k,3578) = 1.46197922197627000D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(844116.00000000000000000D0+rk&
                *(69149.00000000000000000D0+1363.00000000000000000D0*rk)) * den(2)
        ckplm(k,3579) = -6.58549199088411000D-14*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(7.69900440000000000D8+rk*(8.32921540000000000D7+rk&
                *(2.71865100000000000D6+24407.00000000000000000D0*rk))) * den(3)
        ckplm(k,3580) = -2.96347139589785000D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(-5.03103132000000000D9+rk*(-6.21337146000000000D8+rk*(-2.30078710000000000D7+rk&
                *(-176130.00000000000000000D0+2587.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3581) = 5.92694279179570000D-14*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(-5.76850762080000000D11+rk*(-7.44566860320000000D10+rk*(-2.37318005000000000D9+rk&
                *(3.66368150000000000D7+rk*(2.46251000000000000D6+20317.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3582) = 1.58051807781219000D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(4.10110535808000000D12+rk*(5.22803978088000000D11+rk*(8.91759085400000000D9+rk&
                *(-1.14612802500000000D9+rk*(-4.60190950000000000D7+rk&
                *(-189903.00000000000000000D0+6161.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3583) = -6.16402050346752000D-12*(rk+13.D0)*(rk+14.D0)*(1.71561182976000000D12+rk&
                *(2.09079571536000000D11+rk*(-1.42291767600000000D9+rk*(-9.80107532000000000D8+rk&
                *(-2.77658550000000000D7+rk*(538939.00000000000000000D0+rk&
                *(23211.00000000000000000D0+97.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3584) = -2.31150768880032000D-11*(rk+13.D0)*(-6.59419246080000000D12+rk&
                *(-7.60649628240000000D11+rk*(2.91103029240000000D10+rk*(5.86514171600000000D9+rk&
                *(8.29083570000000000D7+rk*(-9.75884000000000000D6+rk*(-267974.00000000000000000D0+rk&
                *(1524.00000000000000000D0+53.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,3585) = -8.73236237991233000D-12*(2.25747744768000000D14+rk&
                *(2.50001454124800000D13+rk*(-1.85413543281600000D12+rk*(-2.74707766964000000D11+rk&
                *(8.24773404000000000D8+rk*(8.39206221000000000D8+rk*(1.40493360000000000D7+rk&
                *(-644406.00000000000000000D0+rk*(-13284.00000000000000000D0+29.00000000000000000D0&
                *rk))))))))) * den(9)
        ckplm(k,3586) = 9.95489311310005000D-10*(1.92803030784000000D12+rk*(5.26185983520000000D10+rk&
                *(-2.75011582680000000D10+rk*(-7.61570776000000000D8+rk*(1.27637223000000000D8+rk&
                *(3.36749700000000000D6+rk*(-210882.00000000000000000D0+rk&
                *(-4674.00000000000000000D0+(rk+87.D0)*rk)))))))) * den(10)
        ckplm(k,3587) = 9.95489311310005000D-10*(-1.84879618560000000D12+rk&
                *(1.04843723040000000D11+rk*(2.44873602480000000D10+rk*(-1.23439543600000000D9+rk&
                *(-1.07806062000000000D8+rk*(4.52988900000000000D6+rk*(175812.00000000000000000D0+rk&
                *(-5334.00000000000000000D0+(rk-78.D0)*rk)))))))) * den(11)
        ckplm(k,3588) = -8.73236237991233000D-12*(-1.99168171937280000D14+rk&
                *(2.78851012143840000D13+rk*(1.03323165414000000D12+rk*(-2.69917592960000000D11+rk&
                *(3.13889698500000000D9+rk*(7.42125237000000000D8+rk*(-1.81857900000000000D7+rk&
                *(-537090.00000000000000000D0+rk*(13545.00000000000000000D0+29.00000000000000000D0*rk))))))))) &
                * den(12)
        ckplm(k,3589) = -2.31150768880032000D-11*(rk-12.D0)*(-5.81020527360000000D12+rk&
                *(8.01653618480000000D11+rk*(1.21058661880000000D10+rk*(-5.44132974000000000D9+rk&
                *(1.27633317000000000D8+rk*(8.12196000000000000D6+rk*(-277158.00000000000000000D0+rk&
                *(-1100.00000000000000000D0+53.00000000000000000D0*rk)))))))) * den(13)
        ckplm(k,3590) = -6.16402050346752000D-12*(rk-12.D0)*(rk-13.D0)*(-1.50606116640000000D12+rk&
                *(2.09098703820000000D11+rk*(-1.34576652800000000D9+rk*(-8.64115547000000000D8+rk&
                *(3.01157800000000000D7+rk*(401710.00000000000000000D0+rk&
                *(-22532.00000000000000000D0+97.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3591) = 1.58051807781219000D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(3.58831927584000000D12+rk*(-5.01713502204000000D11+rk*(1.20818518040000000D10+rk&
                *(9.64073895000000000D8+rk*(-4.49771650000000000D7+rk&
                *(226869.00000000000000000D0+6161.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3592) = 5.92694279179570000D-14*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(5.04801450720000000D11+rk*(-6.96101639420000000D10+rk*(2.46851860500000000D9+rk&
                *(2.69899450000000000D7+rk*(-2.36092500000000000D6+20317.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3593) = -2.96347139589785000D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(-4.43252332800000000D9+rk*(5.75860142000000000D8+rk*(-2.24639590000000000D7+rk&
                *(186478.00000000000000000D0+2587.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3594) = -6.58549199088411000D-14*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(-6.89302530000000000D8+rk*(7.79280730000000000D7+rk&
                *(-2.64543000000000000D6+24407.00000000000000000D0*rk))) * den(18)
        ckplm(k,3595) = 1.46197922197627000D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(776330.00000000000000000D0+rk&
                *(-66423.00000000000000000D0+1363.00000000000000000D0*rk)) * den(19)
        ckplm(k,3596) = -9.50286494284577000D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(-1879.00000000000000000D0+89.00000000000000000D0&
                *rk) * den(20)
        ckplm(k,3597) = 1.29872487552226000D-10*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3598) = -3.81977904565369000D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3599) = 9.31653425769193000D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(6929.00000000000000000D0+317.00000000000000000D0*rk) * den(1)
        ckplm(k,3600) = -4.29993888816551000D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(1.16001600000000000D6+rk&
                *(98999.00000000000000000D0+2063.00000000000000000D0*rk)) * den(2)
        ckplm(k,3601) = 1.93690940908356000D-15*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(1.22984004000000000D9+rk*(1.44094999000000000D8+rk&
                *(5.31687600000000000D6+59837.00000000000000000D0*rk))) * den(3)
        ckplm(k,3602) = -9.68454704541781000D-16*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(8.30268675600000000D10+rk*(1.16606291300000000D10+rk*(5.47337243000000000D8+rk&
                *(9.08704600000000000D6+24181.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3603) = -2.96347139589785000D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(-7.02335134800000000D10+rk*(-1.09309092660000000D10+rk*(-5.50026525000000000D8+rk&
                *(-7.23930500000000000D6+rk*(140445.00000000000000000D0+2771.00000000000000000D0*rk))))) &
                * den(5)
        ckplm(k,3604) = 7.90259038906093000D-14*(rk+14.D0)*(rk+15.D0)*(-5.52815115840000000D11+rk&
                *(-9.08871145440000000D10+rk*(-4.29573866200000000D9+rk*(1.65430350000000000D7+rk&
                *(5.74839500000000000D6+rk*(109149.00000000000000000D0+227.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3605) = 1.02733675057792000D-12*(rk+14.D0)*(7.52666091840000000D11+rk&
                *(1.27357492224000000D11+rk*(4.95098266600000000D9+rk*(-2.24821163000000000D8+rk&
                *(-1.87737200000000000D7+rk*(-263774.00000000000000000D0+rk&
                *(4174.00000000000000000D0+73.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3606) = 4.28056979407467000D-13*(-2.76839192448000000D13+rk&
                *(-4.78322801784000000D12+rk*(-1.28990037276000000D11+rk*(1.88164915960000000D10+rk&
                *(1.17017908700000000D9+rk*(1.88776000000000000D6+rk*(-1.07367400000000000D6+rk&
                *(-16556.00000000000000000D0+23.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,3607) = -1.11294814645941000D-12*(-1.10453857920000000D13+rk&
                *(-1.11478549464000000D12+rk*(5.82954094280000000D10+rk*(8.24106239600000000D9+rk&
                *(2.06850910000000000D7+rk*(-1.48288000000000000D7+rk*(-230258.00000000000000000D0+rk&
                *(4724.00000000000000000D0+59.00000000000000000D0*rk)))))))) * den(9)
        ckplm(k,3608) = -4.22920295654577000D-11*(2.89587432960000000D11+rk*(7.76685676800000000D9+rk&
                *(-2.94536808400000000D9+rk*(-7.64258880000000000D7+rk*(8.95122900000000000D6+rk&
                *(206052.00000000000000000D0+rk*(-8106.00000000000000000D0+(rk-132.D0)*rk))))))) * den(10)
        ckplm(k,3609) = 4.22920295654577000D-11*(2.78960371200000000D11+rk&
                *(-1.33935883200000000D10+rk*(-2.66456235600000000D9+rk*(1.10012840000000000D8+rk&
                *(7.80406900000000000D6+rk*(-251860.00000000000000000D0+rk&
                *(-7154.00000000000000000D0+(rk+140.D0)*rk))))))) * den(11)
        ckplm(k,3610) = 1.11294814645941000D-12*(-9.88051067136000000D12+rk&
                *(1.20680859652800000D12+rk*(3.38410693640000000D10+rk*(-8.01480122800000000D9+rk&
                *(9.12140110000000000D7+rk*(1.33513520000000000D7+rk*(-261674.00000000000000000D0+rk&
                *(-4252.00000000000000000D0+59.00000000000000000D0*rk)))))))) * den(12)
        ckplm(k,3611) = -4.28056979407467000D-13*(-2.30473305216000000D13+rk&
                *(4.47346342008000000D12+rk*(-1.78453071932000000D11+rk*(-1.41755455800000000D10+rk&
                *(1.14521624700000000D9+rk*(-7.98084000000000000D6+rk*(-957138.00000000000000000D0+rk&
                *(16740.00000000000000000D0+23.00000000000000000D0*rk)))))))) * den(13)
        ckplm(k,3612) = -1.02733675057792000D-12*(rk-13.D0)*(-6.30465897600000000D11+rk&
                *(1.16854814880000000D11+rk*(-5.51550265200000000D9+rk*(-1.52444948000000000D8+rk&
                *(1.73947950000000000D7+rk*(-287285.00000000000000000D0+rk&
                *(-3663.00000000000000000D0+73.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3613) = -7.90259038906093000D-14*(rk-13.D0)*(rk-14.D0)*(-4.66234643520000000D11+rk&
                *(8.22684573120000000D10+rk*(-4.31196548200000000D9+rk*(5.36359500000000000D6+rk&
                *(5.20605500000000000D6+rk*(-107787.00000000000000000D0+227.00000000000000000D0*rk)))))) &
                * den(15)
        ckplm(k,3614) = 2.96347139589785000D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(5.98452537600000000D10+rk*(-9.85312205600000000D9+rk*(5.27493650000000000D8+rk&
                *(-7.77337500000000000D6+rk*(-126590.00000000000000000D0+2771.00000000000000000D0*rk))))) &
                * den(16)
        ckplm(k,3615) = 9.68454704541781000D-16*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(7.19045128080000000D10+rk*(-1.05931190580000000D10+rk*(5.20221191000000000D8+rk&
                *(-8.99032200000000000D6+24181.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3616) = -1.93690940908356000D-15*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(-1.09100208000000000D9+rk*(1.33640758000000000D8+rk&
                *(-5.13736500000000000D6+59837.00000000000000000D0*rk))) * den(18)
        ckplm(k,3617) = 4.29993888816551000D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(1.06308000000000000D6+rk&
                *(-94873.00000000000000000D0+2063.00000000000000000D0*rk)) * den(19)
        ckplm(k,3618) = -9.31653425769193000D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(-6612.00000000000000000D0+317.00000000000000000D0*rk) &
                * den(20)
        ckplm(k,3619) = 3.81977904565369000D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3620) = 1.09136544161534000D-13*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3621) = -1.86330685153839000D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(rk+20.D0)*(1148.00000000000000000D0+53.00000000000000000D0*rk) * den(1)
        ckplm(k,3622) = 4.29993888816551000D-14*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(44499.00000000000000000D0+rk*(3919.00000000000000000D0+85.00000000000000000D0&
                *rk)) * den(2)
        ckplm(k,3623) = -3.87381881816712000D-16*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(2.71476030000000000D8+rk*(3.37694830000000000D7+rk&
                *(1.35310500000000000D6+17237.00000000000000000D0*rk))) * den(3)
        ckplm(k,3624) = 1.93690940908356000D-16*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(2.08762146600000000D10+rk*(3.22146665600000000D9+rk*(1.74635567000000000D8+rk&
                *(3.79509400000000000D6+25363.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3625) = 5.92694279179570000D-15*(rk+15.D0)*(rk+16.D0)*(-1.98300160800000000D10+rk&
                *(-3.52674682800000000D9+rk*(-2.24222060000000000D8+rk*(-5.75788500000000000D6+rk&
                *(-37090.00000000000000000D0+363.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3626) = -1.58051807781219000D-14*(rk+15.D0)*(-1.72286336880000000D11+rk&
                *(-3.37673610180000000D10+rk*(-2.30566564300000000D9+rk*(-5.39386950000000000D7+rk&
                *(502550.00000000000000000D0+rk*(34833.00000000000000000D0+293.00000000000000000D0*rk)))))) &
                * den(6)
        ckplm(k,3627) = -1.46762392939703000D-13*(3.55767340320000000D11+rk&
                *(7.49145295320000000D10+rk*(5.18025569600000000D9+rk*(7.23514190000000000D7+rk&
                *(-6.31046500000000000D6+rk*(-254002.00000000000000000D0+rk&
                *(-2191.00000000000000000D0+11.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3628) = 8.56113958814934000D-14*(7.09518398400000000D11+rk*(1.08180897120000000D11+rk&
                *(3.01448107800000000D9+rk*(-2.47350397000000000D8+rk*(-1.39120800000000000D7+rk&
                *(-71890.00000000000000000D0+rk*(5082.00000000000000000D0+47.00000000000000000D0*rk))))))) &
                * den(8)
        ckplm(k,3629) = 8.56113958814934000D-14*(-7.67432332800000000D11+rk&
                *(-6.92010100800000000D10+rk*(2.48618255200000000D9+rk*(3.29513058000000000D8+rk&
                *(1.79714500000000000D6+rk*(-314685.00000000000000000D0+rk&
                *(-3857.00000000000000000D0+27.00000000000000000D0*rk))))))) * den(9)
        ckplm(k,3630) = -3.25323304349675000D-12*(-2.05899321600000000D10+rk&
                *(-5.23065810000000000D8+rk*(1.44760735000000000D8+rk*(3.37023400000000000D6+rk&
                *(-270095.00000000000000000D0+rk*(-5033.00000000000000000D0+(rk+112.D0)*rk)))))) * den(10)
        ckplm(k,3631) = -3.25323304349675000D-12*(1.99257408000000000D10+rk&
                *(-8.01422028000000000D8+rk*(-1.33081452000000000D8+rk*(4.39807900000000000D6+rk&
                *(243285.00000000000000000D0+rk*(-5684.00000000000000000D0+(rk-105.D0)*rk)))))) * den(11)
        ckplm(k,3632) = 8.56113958814934000D-14*(6.96072545280000000D11+rk&
                *(-7.31935746840000000D10+rk*(-1.51151467600000000D9+rk*(3.19255713000000000D8+rk&
                *(-3.31177000000000000D6+rk*(-290976.00000000000000000D0+rk&
                *(4046.00000000000000000D0+27.00000000000000000D0*rk))))))) * den(12)
        ckplm(k,3633) = 8.56113958814934000D-14*(-6.04585497600000000D11+rk&
                *(1.01465142480000000D11+rk*(-3.67385393200000000D9+rk*(-1.92520972000000000D8+rk&
                *(1.34780450000000000D7+rk*(-101395.00000000000000000D0+rk&
                *(-4753.00000000000000000D0+47.00000000000000000D0*rk))))))) * den(13)
        ckplm(k,3634) = -1.46762392939703000D-13*(-2.85954656400000000D11+rk&
                *(6.47950574700000000D10+rk*(-4.92784557300000000D9+rk*(9.50974640000000000D7+rk&
                *(5.07370500000000000D6+rk*(-240625.00000000000000000D0+rk&
                *(2268.00000000000000000D0+11.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,3635) = -1.58051807781219000D-14*(rk-14.D0)*(-1.40770234800000000D11+rk&
                *(2.93196836100000000D10+rk*(-2.14117819300000000D9+rk*(5.56064250000000000D7+rk&
                *(332780.00000000000000000D0+rk*(-33075.00000000000000000D0+293.00000000000000000D0*rk)))))) &
                * den(15)
        ckplm(k,3636) = 5.92694279179570000D-15*(rk-14.D0)*(rk-15.D0)*(1.65217708800000000D10+rk&
                *(-3.09542618800000000D9+rk*(2.07174575000000000D8+rk*(-5.60589500000000000D6+rk&
                *(38905.00000000000000000D0+363.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3637) = 1.93690940908356000D-16*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(1.78256138400000000D10+rk*(-2.88347935200000000D9+rk*(1.63402463000000000D8+rk&
                *(-3.69364200000000000D6+25363.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3638) = -3.87381881816712000D-16*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(-2.39042415000000000D8+rk*(3.11149840000000000D7+rk&
                *(-1.30139400000000000D6+17237.00000000000000000D0*rk))) * den(18)
        ckplm(k,3639) = 4.29993888816551000D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(40665.00000000000000000D0+rk*(-3749.00000000000000000D0+85.00000000000000000D0&
                *rk)) * den(19)
        ckplm(k,3640) = 1.86330685153839000D-14*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(rk-19.D0)*(1095.00000000000000000D0-53.00000000000000000D0*rk) * den(20)
        ckplm(k,3641) = 1.09136544161534000D-13*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3642) = -3.03157067115372000D-15*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3643) = 2.21822244230760000D-16*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(rk+20.D0)*(3075.00000000000000000D0+143.00000000000000000D0*rk) * den(1)
        ckplm(k,3644) = -1.70632495562123000D-16*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(409890.00000000000000000D0+rk*(37003.00000000000000000D0+827.00000000000000000D0*rk)) * den(2)
        ckplm(k,3645) = 7.68614844874429000D-18*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(5.70541110000000000D8+rk*(7.43351810000000000D7+rk&
                *(3.16028400000000000D6+43585.00000000000000000D0*rk))) * den(3)
        ckplm(k,3646) = -6.91753360386986000D-18*(rk+16.D0)*(rk+17.D0)*(2.75694022800000000D10+rk&
                *(4.57077103800000000D9+rk*(2.73304567000000000D8+rk&
                *(6.88051800000000000D6+59837.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3647) = 2.35196142531575000D-17*(rk+16.D0)*(2.63348924760000000D11+rk&
                *(5.17686774060000000D10+rk*(3.81677732500000000D9+rk*(1.27413365000000000D8+rk&
                *(1.77183500000000000D6+5869.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3648) = 6.27189713417534000D-17*(-2.51795224080000000D12+rk&
                *(-5.61973182420000000D11+rk*(-4.73703321640000000D10+rk*(-1.79358982500000000D9+rk&
                *(-2.48315350000000000D7+rk*(129165.00000000000000000D0+4139.00000000000000000D0*rk)))))) &
                * den(6)
        ckplm(k,3649) = -1.22301994116419000D-14*(-1.78553332800000000D10+rk&
                *(-3.21830374800000000D9+rk*(-1.91148592000000000D8+rk*(-2.88259500000000000D6+rk&
                *(99815.00000000000000000D0+rk*(3183.00000000000000000D0+17.00000000000000000D0*rk)))))) &
                * den(7)
        ckplm(k,3650) = -1.52877492645524000D-14*(1.76508288000000000D10+rk*(2.32195548000000000D9+rk&
                *(6.16397620000000000D7+rk*(-2.84227500000000000D6+rk*(-136295.00000000000000000D0+rk&
                *(-765.00000000000000000D0+13.00000000000000000D0*rk)))))) * den(8)
        ckplm(k,3651) = 1.69863880717249000D-15*(1.79143895808000000D11+rk*(1.40923766640000000D10+rk&
                *(-3.34346294000000000D8+rk*(-4.05126870000000000D7+rk*(-257483.00000000000000000D0+rk&
                *(17103.00000000000000000D0+121.00000000000000000D0*rk)))))) * den(9)
        ckplm(k,3652) = 3.87289648035327000D-14*(-8.19801907200000000D9+rk&
                *(-1.90399692000000000D8+5.00000000000000000D0*rk*(7.62361600000000000D6+rk&
                *(151221.00000000000000000D0+rk*(-7889.00000000000000000D0+(rk-105.D0)*rk))))) * den(10)
        ckplm(k,3653) = -3.87289648035327000D-14*(-7.97029632000000000D9+rk&
                *(2.64212412000000000D8+5.00000000000000000D0*rk*(7.12368400000000000D6+rk&
                *(-181707.00000000000000000D0+rk*(-7349.00000000000000000D0+(rk+111.D0)*rk))))) * den(11)
        ckplm(k,3654) = -1.69863880717249000D-15*(1.64757411072000000D11+rk&
                *(-1.46406459120000000D10+rk*(-2.14522346000000000D8+rk*(3.93141450000000000D7+rk&
                *(-341183.00000000000000000D0+rk*(-16377.00000000000000000D0+121.00000000000000000D0*rk)))))) &
                * den(12)
        ckplm(k,3655) = 1.52877492645524000D-14*(1.53932198400000000D10+rk*(-2.19069040800000000D9+rk&
                *(6.93566620000000000D7+rk*(2.30500500000000000D6+rk*(-132275.00000000000000000D0+rk&
                *(843.00000000000000000D0+13.00000000000000000D0*rk)))))) * den(13)
        ckplm(k,3656) = 1.22301994116419000D-14*(-1.48251988800000000D10+rk*(2.84503779600000000D9+rk&
                *(-1.81933492000000000D8+rk*(3.25036500000000000D6+rk*(84155.00000000000000000D0+rk&
                *(-3081.00000000000000000D0+17.00000000000000000D0*rk)))))) * den(14)
        ckplm(k,3657) = -6.27189713417534000D-17*(-2.00158075728000000D12+rk&
                *(4.72513340436000000D11+rk*(-4.21397814640000000D10+rk*(1.69305481500000000D9+rk&
                *(-2.54152750000000000D7+rk*(-104331.00000000000000000D0+4139.00000000000000000D0*rk)))))) &
                * den(15)
        ckplm(k,3658) = -2.35196142531575000D-17*(rk-15.D0)*(-2.15271377280000000D11+rk&
                *(4.45103048560000000D10+rk*(-3.44510955000000000D9+rk*(1.20384715000000000D8+rk&
                *(-1.74249000000000000D6+5869.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3659) = 6.91753360386986000D-18*(rk-15.D0)*(rk-16.D0)*(2.32651151280000000D10+rk&
                *(-4.04456411000000000D9+rk*(2.53022035000000000D8+rk&
                *(-6.64117000000000000D6+59837.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3660) = -7.68614844874429000D-18*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(-4.99322628000000000D8+rk*(6.81453680000000000D7+rk&
                *(-3.02952900000000000D6+43585.00000000000000000D0*rk))) * den(18)
        ckplm(k,3661) = 1.70632495562123000D-16*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(373714.00000000000000000D0+rk*(-35349.00000000000000000D0+827.00000000000000000D0*rk)) &
                * den(19)
        ckplm(k,3662) = 2.21822244230760000D-16*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(rk-19.D0)*(2932.00000000000000000D0-143.00000000000000000D0*rk) * den(20)
        ckplm(k,3663) = 3.03157067115372000D-15*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3664) = 8.19343424636142000D-17*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(rk+21.D0) * den(0)
        ckplm(k,3665) = -1.99839859667352000D-18*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(10496.00000000000000000D0+491.00000000000000000D0*rk) * den(1)
        ckplm(k,3666) = 1.53722968974886000D-18*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(1.58995200000000000D6+rk*(146411.00000000000000000D0+3349.00000000000000000D0*rk)) * den(2)
        ckplm(k,3667) = -1.53722968974886000D-18*(rk+17.D0)*(rk+18.D0)*(1.12841600000000000D8+rk&
                *(1.52517020000000000D7+rk*(677935.00000000000000000D0+9883.00000000000000000D0*rk))) * den(3)
        ckplm(k,3668) = 6.91753360386986000D-18*(rk+17.D0)*(1.22558912000000000D9+rk&
                *(2.14953530000000000D8+rk*(1.38035630000000000D7+rk&
                *(382138.00000000000000000D0+3809.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3669) = -7.83987141771918000D-18*(3.90590899200000000D10+rk*(8.29974962400000000D9+rk&
                *(6.78880070000000000D8+rk*(2.63488150000000000D7+rk&
                *(473230.00000000000000000D0+2981.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3670) = -6.27189713417534000D-17*(-8.53682580000000000D9+rk&
                *(-1.57272725800000000D9+rk*(-1.06206995000000000D8+rk*(-3.09822500000000000D6+rk&
                *(-32605.00000000000000000D0+3.00000000000000000D0*rk))))) * den(6)
        ckplm(k,3671) = 4.07673313721397000D-15*(-1.96147152000000000D8+rk*(-2.93867800000000000D7+rk&
                *(-1.42938000000000000D6+rk*(-19465.00000000000000000D0+rk&
                *(252.00000000000000000D0+5.00000000000000000D0*rk))))) * den(7)
        ckplm(k,3672) = -3.05754985291048000D-15*(-3.43054440000000000D8+rk&
                *(-3.77947780000000000D7+rk*(-887925.00000000000000000D0+rk*(21855.00000000000000000D0+rk&
                *(805.00000000000000000D0+3.00000000000000000D0*rk))))) * den(8)
        ckplm(k,3673) = -3.39727761434498000D-16*(3.62974752000000000D9+rk*(2.41556712000000000D8+rk&
                *(-3.56343400000000000D6+rk*(-378611.00000000000000000D0+rk&
                *(-2114.00000000000000000D0+47.00000000000000000D0*rk))))) * den(9)
        ckplm(k,3674) = 1.29096549345109000D-14*(1.02000816000000000D8+rk*(2.08297000000000000D6+rk&
                *(-293901.00000000000000000D0+rk*(-4619.00000000000000000D0+(rk+141.D0)*rk)))) * den(10)
        ckplm(k,3675) = 1.29096549345109000D-14*(-9.96287040000000000D7+rk*(2.65635600000000000D6+rk&
                *(279208.00000000000000000D0+rk*(-5173.00000000000000000D0+(rk-136.D0)*rk)))) * den(11)
        ckplm(k,3676) = -3.39727761434498000D-16*(-3.38500382400000000D9+rk*(2.47556438000000000D8+rk&
                *(2.44075500000000000D6+rk*(-369685.00000000000000000D0+rk&
                *(2349.00000000000000000D0+47.00000000000000000D0*rk))))) * den(12)
        ckplm(k,3677) = -3.05754985291048000D-15*(3.06168640000000000D8+rk*(-3.59565680000000000D7+rk&
                *(948690.00000000000000000D0+rk*(18665.00000000000000000D0+rk&
                *(-790.00000000000000000D0+3.00000000000000000D0*rk))))) * den(13)
        ckplm(k,3678) = 4.07673313721397000D-15*(1.68170040000000000D8+rk*(-2.65873980000000000D7+rk&
                *(1.36952300000000000D6+rk*(-20423.00000000000000000D0+rk&
                *(-227.00000000000000000D0+5.00000000000000000D0*rk))))) * den(14)
        ckplm(k,3679) = -6.27189713417534000D-17*(7.06723992000000000D9+rk*(-1.36947750800000000D9+rk&
                *(9.71079800000000000D7+rk*(-2.96777500000000000D6+rk&
                *(32620.00000000000000000D0+3.00000000000000000D0*rk))))) * den(15)
        ckplm(k,3680) = -7.83987141771918000D-18*(-3.14123418000000000D10+rk&
                *(7.01915791400000000D9+rk*(-6.02643195000000000D8+rk*(2.44857050000000000D7+rk&
                *(-458325.00000000000000000D0+2981.00000000000000000D0*rk))))) * den(16)
        ckplm(k,3681) = 6.91753360386986000D-18*(rk-16.D0)*(1.02406082400000000D9+rk&
                *(-1.88477582000000000D8+rk*(1.26800030000000000D7+rk&
                *(-366902.00000000000000000D0+3809.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3682) = -1.53722968974886000D-18*(rk-16.D0)*(rk-17.D0)*(-9.82579500000000000D7+rk&
                *(1.39254810000000000D7+rk*(-648286.00000000000000000D0+9883.00000000000000000D0*rk))) * den(18)
        ckplm(k,3683) = 1.53722968974886000D-18*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(1.44689000000000000D6+rk*(-139713.00000000000000000D0+3349.00000000000000000D0*rk)) * den(19)
        ckplm(k,3684) = 1.99839859667352000D-18*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(10005.00000000000000000D0-491.00000000000000000D0*rk) * den(20)
        ckplm(k,3685) = 8.19343424636142000D-17*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(rk-20.D0) * den(21)
        ckplm(k,3686) = -2.15616690693721000D-18*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3687) = 5.25894367545662000D-20*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(11849.00000000000000000D0+557.00000000000000000D0*rk) * den(1)
        ckplm(k,3688) = -8.09068257762557000D-21*(rk+18.D0)*(rk+19.D0)*(1.01213580000000000D7+rk&
                *(947237.00000000000000000D0+22069.00000000000000000D0*rk)) * den(2)
        ckplm(k,3689) = 7.68614844874429000D-19*(rk+18.D0)*(8.49949000000000000D6+rk&
                *(1.18337500000000000D6+rk*(54452.00000000000000000D0+827.00000000000000000D0*rk))) * den(3)
        ckplm(k,3690) = -6.91753360386986000D-19*(5.15703160000000000D8+rk*(9.46235620000000000D7+rk&
                *(6.41537300000000000D6+rk*(189922.00000000000000000D0+2063.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3691) = 3.91993570885959000D-18*(2.14086456000000000D8+rk*(3.58137420000000000D7+rk&
                *(2.17109900000000000D6+rk*(55878.00000000000000000D0+505.00000000000000000D0*rk)))) * den(5)
        ckplm(k,3692) = -3.13594856708767000D-17*(5.11918200000000000D7+rk*(7.46872000000000000D6+rk&
                *(377903.00000000000000000D0+rk*(7490.00000000000000000D0+43.00000000000000000D0*rk)))) * den(6)
        ckplm(k,3693) = -4.07673313721397000D-16*(-6.34179600000000000D6+rk&
                *(-757420.00000000000000000D0+rk*(-28135.00000000000000000D0+(rk-290.D0)*rk))) * den(7)
        ckplm(k,3694) = 1.52877492645524000D-15*(rk+27.D0)*(-87160.00000000000000000D0+rk&
                *(-4482.00000000000000000D0+(rk+19.D0)*rk)) * den(8)
        ckplm(k,3695) = -5.77537194438646000D-16*(-7.63380000000000000D6+rk&
                *(-411726.00000000000000000D0+rk*(3431.00000000000000000D0+(rk+294.D0)*rk))) * den(9)
        ckplm(k,3696) = -1.15507438887729000D-15*(4.17171600000000000D6+rk&
                *(71324.00000000000000000D0+rk*(-6727.00000000000000000D0+(rk-74.D0)*rk))) * den(10)
        ckplm(k,3697) = 1.15507438887729000D-15*(4.09374000000000000D6+rk&
                *(-84552.00000000000000000D0+rk*(-6499.00000000000000000D0+(rk+78.D0)*rk))) * den(11)
        ckplm(k,3698) = 5.77537194438646000D-16*(-7.21893600000000000D6+rk&
                *(417710.00000000000000000D0+rk*(2555.00000000000000000D0+(rk-290.D0)*rk))) * den(12)
        ckplm(k,3699) = -1.52877492645524000D-15*(rk-26.D0)*(82660.00000000000000000D0+rk&
                *(-4517.00000000000000000D0+(rk-16.D0)*rk)) * den(13)
        ckplm(k,3700) = 4.07673313721397000D-16*(-5.61222000000000000D6+rk&
                *(702024.00000000000000000D0+rk*(-27259.00000000000000000D0+(rk+294.D0)*rk))) * den(14)
        ckplm(k,3701) = 3.13594856708767000D-17*(4.40935560000000000D7+rk*(-6.73521200000000000D6+rk&
                *(355691.00000000000000000D0+rk*(-7318.00000000000000000D0+43.00000000000000000D0*rk)))) &
                * den(15)
        ckplm(k,3702) = -3.91993570885959000D-18*(1.80388440000000000D8+rk*(-3.16371580000000000D7+rk&
                *(2.00649500000000000D6+rk*(-53858.00000000000000000D0+505.00000000000000000D0*rk)))) * den(16)
        ckplm(k,3703) = 6.91753360386986000D-19*(4.27307112000000000D8+rk*(-8.23543300000000000D7+rk&
                *(5.85798500000000000D6+rk*(-181670.00000000000000000D0+2063.00000000000000000D0*rk)))) &
                * den(17)
        ckplm(k,3704) = -7.68614844874429000D-19*(rk-17.D0)*(-7.36974000000000000D6+rk&
                *(1.07695200000000000D6+rk*(-51971.00000000000000000D0+827.00000000000000000D0*rk))) * den(18)
        ckplm(k,3705) = 8.09068257762557000D-21*(rk-17.D0)*(rk-18.D0)*(9.19619000000000000D6+rk&
                *(-903099.00000000000000000D0+22069.00000000000000000D0*rk)) * den(19)
        ckplm(k,3706) = 5.25894367545662000D-20*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(11292.00000000000000000D0-557.00000000000000000D0*rk) * den(20)
        ckplm(k,3707) = 2.15616690693721000D-18*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3708) = 5.52863309471081000D-20*(rk+19.D0)*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3709) = -4.04534128881278000D-21*(rk+19.D0)*(rk+20.D0)&
                *(4428.00000000000000000D0+209.00000000000000000D0*rk) * den(1)
        ckplm(k,3710) = 8.09068257762557000D-21*(rk+19.D0)*(325917.00000000000000000D0+rk&
                *(30913.00000000000000000D0+731.00000000000000000D0*rk)) * den(2)
        ckplm(k,3711) = -2.56204948291476000D-19*(917190.00000000000000000D0+rk&
                *(130847.00000000000000000D0+rk*(6189.00000000000000000D0+97.00000000000000000D0*rk))) * den(3)
        ckplm(k,3712) = 2.30584453462329000D-19*(3.43581000000000000D6+rk&
                *(463207.00000000000000000D0+rk*(20544.00000000000000000D0+299.00000000000000000D0*rk))) &
                * den(4)
        ckplm(k,3713) = -1.17598071265788000D-17*(174744.00000000000000000D0+rk&
                *(21574.00000000000000000D0+rk*(861.00000000000000000D0+11.00000000000000000D0*rk))) * den(5)
        ckplm(k,3714) = 3.13594856708767000D-17*(136755.00000000000000000D0+rk&
                *(14806.00000000000000000D0+rk*(498.00000000000000000D0+5.00000000000000000D0*rk))) * den(6)
        ckplm(k,3715) = -9.40784570126301000D-17*(79038.00000000000000000D0+rk&
                *(7049.00000000000000000D0+(rk+177.D0)*rk)) * den(7)
        ckplm(k,3716) = -3.91993570885959000D-17*(-279840.00000000000000000D0+rk&
                *(-18613.00000000000000000D0+(rk-252.D0)*rk)) * den(8)
        ckplm(k,3717) = 1.13242587144833000D-17*(-1.23720000000000000D6+rk&
                *(-50612.00000000000000000D0+rk*(201.00000000000000000D0+11.00000000000000000D0*rk))) * den(9)
        ckplm(k,3718) = -6.79455522868995000D-17*(-230679.00000000000000000D0+rk&
                *(-3076.00000000000000000D0+(rk+174.D0)*rk)) * den(10)
        ckplm(k,3719) = -6.79455522868995000D-17*(227430.00000000000000000D0+rk&
                *(-3421.00000000000000000D0+(rk-171.D0)*rk)) * den(11)
        ckplm(k,3720) = 1.13242587144833000D-17*(1.18639800000000000D6+rk&
                *(-50981.00000000000000000D0+rk*(-168.00000000000000000D0+11.00000000000000000D0*rk))) * den(12)
        ckplm(k,3721) = -3.91993570885959000D-17*(261480.00000000000000000D0+rk&
                *(-18106.00000000000000000D0+(rk+255.D0)*rk)) * den(13)
        ckplm(k,3722) = -9.40784570126301000D-17*(-72165.00000000000000000D0+rk&
                *(6698.00000000000000000D0+(rk-174.D0)*rk)) * den(14)
        ckplm(k,3723) = 3.13594856708767000D-17*(-122442.00000000000000000D0+rk&
                *(13825.00000000000000000D0+rk*(-483.00000000000000000D0+5.00000000000000000D0*rk))) * den(15)
        ckplm(k,3724) = -1.17598071265788000D-17*(-154020.00000000000000000D0+rk&
                *(19885.00000000000000000D0+rk*(-828.00000000000000000D0+11.00000000000000000D0*rk))) * den(16)
        ckplm(k,3725) = 2.30584453462329000D-19*(-2.99284800000000000D6+rk&
                *(423016.00000000000000000D0+rk*(-19647.00000000000000000D0+299.00000000000000000D0*rk))) &
                * den(17)
        ckplm(k,3726) = -2.56204948291476000D-19*(-792435.00000000000000000D0+rk&
                *(118760.00000000000000000D0+rk*(-5898.00000000000000000D0+97.00000000000000000D0*rk))) &
                * den(18)
        ckplm(k,3727) = 8.09068257762557000D-21*(rk-18.D0)*(295735.00000000000000000D0+rk&
                *(-29451.00000000000000000D0+731.00000000000000000D0*rk)) * den(19)
        ckplm(k,3728) = 4.04534128881278000D-21*(rk-18.D0)*(rk-19.D0)&
                *(4219.00000000000000000D0-209.00000000000000000D0*rk) * den(20)
        ckplm(k,3729) = 5.52863309471081000D-20*(rk-18.D0)*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3730) = -1.38215827367770000D-21*(rk+20.D0)*(rk+21.D0) * den(0)
        ckplm(k,3731) = 3.37111774067732000D-23*(rk+20.D0)&
                *(14801.00000000000000000D0+701.00000000000000000D0*rk) * den(1)
        ckplm(k,3732) = -3.37111774067732000D-22*(242592.00000000000000000D0+rk&
                *(23269.00000000000000000D0+557.00000000000000000D0*rk)) * den(2)
        ckplm(k,3733) = 6.40512370728691000D-21*(66520.00000000000000000D0+rk&
                *(6183.00000000000000000D0+143.00000000000000000D0*rk)) * den(3)
        ckplm(k,3734) = -9.60768556093036000D-21*(165340.00000000000000000D0+rk&
                *(14565.00000000000000000D0+317.00000000000000000D0*rk)) * den(4)
        ckplm(k,3735) = 3.26661309071632000D-20*(138540.00000000000000000D0+rk&
                *(11221.00000000000000000D0+221.00000000000000000D0*rk)) * den(5)
        ckplm(k,3736) = -2.61329047257306000D-19*(39320.00000000000000000D0+rk&
                *(2807.00000000000000000D0+47.00000000000000000D0*rk)) * den(6)
        ckplm(k,3737) = 1.30664523628653000D-18*(14656.00000000000000000D0+rk&
                *(867.00000000000000000D0+11.00000000000000000D0*rk)) * den(7)
        ckplm(k,3738) = -3.26661309071632000D-19*(91620.00000000000000000D0+rk&
                *(4069.00000000000000000D0+29.00000000000000000D0*rk)) * den(8)
        ckplm(k,3739) = -1.41553233931041000D-18*(-28140.00000000000000000D0+(rk-775.D0)*rk) * den(9)
        ckplm(k,3740) = 1.07580457787591000D-17*(-4240.00000000000000000D0+(rk-39.D0)*rk) * den(10)
        ckplm(k,3741) = -1.07580457787591000D-17*(-4200.00000000000000000D0+(rk+41.D0)*rk) * den(11)
        ckplm(k,3742) = 1.41553233931041000D-18*(-27364.00000000000000000D0+(rk+777.D0)*rk) * den(12)
        ckplm(k,3743) = 3.26661309071632000D-19*(87580.00000000000000000D0+rk&
                *(-4011.00000000000000000D0+29.00000000000000000D0*rk)) * den(13)
        ckplm(k,3744) = -1.30664523628653000D-18*(13800.00000000000000000D0+rk&
                *(-845.00000000000000000D0+11.00000000000000000D0*rk)) * den(14)
        ckplm(k,3745) = 2.61329047257306000D-19*(36560.00000000000000000D0+rk&
                *(-2713.00000000000000000D0+47.00000000000000000D0*rk)) * den(15)
        ckplm(k,3746) = -3.26661309071632000D-20*(127540.00000000000000000D0+rk&
                *(-10779.00000000000000000D0+221.00000000000000000D0*rk)) * den(16)
        ckplm(k,3747) = 9.60768556093036000D-21*(151092.00000000000000000D0+rk&
                *(-13931.00000000000000000D0+317.00000000000000000D0*rk)) * den(17)
        ckplm(k,3748) = -6.40512370728691000D-21*(60480.00000000000000000D0+rk&
                *(-5897.00000000000000000D0+143.00000000000000000D0*rk)) * den(18)
        ckplm(k,3749) = 3.37111774067732000D-22*(219880.00000000000000000D0+rk&
                *(-22155.00000000000000000D0+557.00000000000000000D0*rk)) * den(19)
        ckplm(k,3750) = 3.37111774067732000D-23*(rk-19.D0)&
                *(14100.00000000000000000D0-701.00000000000000000D0*rk) * den(20)
        ckplm(k,3751) = 1.38215827367770000D-21*(rk-19.D0)*(rk-20.D0) * den(21)
        ckplm(k,3752) = 3.37111774067732000D-23*(rk+21.D0) * den(0)
        ckplm(k,3753) = -3.37111774067732000D-23*(400.00000000000000000D0+19.00000000000000000D0*rk) &
                * den(1)
        ckplm(k,3754) = 3.37111774067732000D-22*(363.00000000000000000D0+17.00000000000000000D0*rk) &
                * den(2)
        ckplm(k,3755) = -3.20256185364345000D-20*(rk+22.D0) * den(3)
        ckplm(k,3756) = 9.60768556093036000D-21*(301.00000000000000000D0+13.00000000000000000D0*rk) &
                * den(4)
        ckplm(k,3757) = -3.26661309071632000D-20*(276.00000000000000000D0+11.00000000000000000D0*rk) &
                * den(5)
        ckplm(k,3758) = 2.61329047257306000D-19*(85.00000000000000000D0+3.00000000000000000D0*rk) &
                * den(6)
        ckplm(k,3759) = -1.30664523628653000D-18*(rk+34.D0) * den(7)
        ckplm(k,3760) = 1.63330654535816000D-18*(rk+45.D0) * den(8)
        ckplm(k,3761) = -1.41553233931041000D-18*(rk+72.D0) * den(9)
        ckplm(k,3762) = 5.66212935724163000D-19*(rk+211.D0) * den(10)
        ckplm(k,3763) = 5.66212935724163000D-19*(rk-210.D0) * den(11)
        ckplm(k,3764) = -1.41553233931041000D-18*(rk-71.D0) * den(12)
        ckplm(k,3765) = 1.63330654535816000D-18*(rk-44.D0) * den(13)
        ckplm(k,3766) = -1.30664523628653000D-18*(rk-33.D0) * den(14)
        ckplm(k,3767) = 2.61329047257306000D-19*(-82.00000000000000000D0+3.00000000000000000D0*rk) &
                * den(15)
        ckplm(k,3768) = 3.26661309071632000D-20*(265.00000000000000000D0-11.00000000000000000D0*rk) &
                * den(16)
        ckplm(k,3769) = 9.60768556093036000D-21*(-288.00000000000000000D0+13.00000000000000000D0*rk) &
                * den(17)
        ckplm(k,3770) = -3.20256185364345000D-20*(rk-21.D0) * den(18)
        ckplm(k,3771) = 3.37111774067732000D-22*(-346.00000000000000000D0+17.00000000000000000D0*rk) &
                * den(19)
        ckplm(k,3772) = 3.37111774067732000D-23*(381.00000000000000000D0-19.00000000000000000D0*rk) &
                * den(20)
        ckplm(k,3773) = 3.37111774067732000D-23*(rk-20.D0) * den(21)
        ckplm(k,3774) = -8.02647081113648000D-25 * den(0)
        ckplm(k,3775) = 1.68555887033866000D-23 * den(1)
        ckplm(k,3776) = -1.68555887033866000D-22 * den(2)
        ckplm(k,3777) = 1.06752061788115000D-21 * den(3)
        ckplm(k,3778) = -4.80384278046518000D-21 * den(4)
        ckplm(k,3779) = 1.63330654535816000D-20 * den(5)
        ckplm(k,3780) = -4.35548412095510000D-20 * den(6)
        ckplm(k,3781) = 9.33318025918950000D-20 * den(7)
        ckplm(k,3782) = -1.63330654535816000D-19 * den(8)
        ckplm(k,3783) = 2.35922056551735000D-19 * den(9)
        ckplm(k,3784) = -2.83106467862081000D-19 * den(10)
        ckplm(k,3785) = 2.83106467862081000D-19 * den(11)
        ckplm(k,3786) = -2.35922056551735000D-19 * den(12)
        ckplm(k,3787) = 1.63330654535816000D-19 * den(13)
        ckplm(k,3788) = -9.33318025918950000D-20 * den(14)
        ckplm(k,3789) = 4.35548412095510000D-20 * den(15)
        ckplm(k,3790) = -1.63330654535816000D-20 * den(16)
        ckplm(k,3791) = 4.80384278046518000D-21 * den(17)
        ckplm(k,3792) = -1.06752061788115000D-21 * den(18)
        ckplm(k,3793) = 1.68555887033866000D-22 * den(19)
        ckplm(k,3794) = -1.68555887033866000D-23 * den(20)
        ckplm(k,3795) = 8.02647081113648000D-25 * den(21)
!    ckplm para l = 22
        den(0) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k+1.D0)&
                *(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+41.D0)*(r2k+43.D0)*(r2k+45.D0)*(r2k+5.D0)&
                *(r2k+7.D0)*(r2k+9.D0))
        den(1) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k+3.D0)*(r2k+41.D0)*(r2k+43.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(2) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k+41.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(3) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k+39.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(4) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k+37.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k+9.D0))
        den(5) = 1.d0/((r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)*(r2k+33.D0)&
                *(r2k+35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(6) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)*(r2k+19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)*(r2k+31.D0)&
                *(r2k+33.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(7) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k+15.D0)*(r2k+17.D0)&
                *(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)*(r2k+29.D0)&
                *(r2k+31.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(8) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)*(r2k+27.D0)&
                *(r2k+29.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(9) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)*(r2k+25.D0)&
                *(r2k+27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(10) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k+21.D0)*(r2k+23.D0)&
                *(r2k+25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(11) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k+21.D0)&
                *(r2k+23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(12) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k+21.D0)&
                *(r2k-23.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(13) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k+19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)&
                *(r2k-25.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(14) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k+17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)&
                *(r2k-27.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(15) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k+15.D0)&
                *(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)&
                *(r2k-29.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(16) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k+13.D0)*(r2k-15.D0)*(r2k-17.D0)&
                *(r2k-19.D0)*(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)&
                *(r2k-31.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(17) = 1.d0/((r2k-11.D0)*(r2k+11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)&
                *(r2k-1.D0)*(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)&
                *(r2k-33.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(18) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)*(r2k-9.D0)&
                *(r2k+9.D0))
        den(19) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)*(r2k+7.D0)&
                *(r2k-9.D0))
        den(20) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-39.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-5.D0)*(r2k+5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        den(21) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-39.D0)*(r2k-3.D0)*(r2k+3.D0)*(r2k-41.D0)*(r2k-5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        den(22) = 1.d0/((r2k-11.D0)*(r2k-13.D0)*(r2k-15.D0)*(r2k-17.D0)*(r2k-19.D0)*(r2k-1.D0)&
                *(r2k+1.D0)*(r2k-21.D0)*(r2k-23.D0)*(r2k-25.D0)*(r2k-27.D0)*(r2k-29.D0)*(r2k-31.D0)*(r2k-33.D0)&
                *(r2k-35.D0)*(r2k-37.D0)*(r2k-39.D0)*(r2k-3.D0)*(r2k-41.D0)*(r2k-43.D0)*(r2k-5.D0)*(r2k-7.D0)&
                *(r2k-9.D0))
        ckplm(k,3796) = 2.25745328348637000D7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+22.D0)*(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) &
                * den(0)
        ckplm(k,3797) = 1.15497609852791000D7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+1.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(1)
        ckplm(k,3798) = 8.87359685454369000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-1.D0)*(rk+1.D0)*(rk+20.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(2)
        ckplm(k,3799) = 7.58427081584930000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(3)
        ckplm(k,3800) = 6.81559471964836000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(4)
        ckplm(k,3801) = 6.30929339761734000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(5)
        ckplm(k,3802) = 5.95877709774971000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(6)
        ckplm(k,3803) = 5.71163887710571000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(7)
        ckplm(k,3804) = 5.53930494546890000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(8)
        ckplm(k,3805) = 5.42532747745514000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*rk * den(9)
        ckplm(k,3806) = 5.36022354772568000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(10)
        ckplm(k,3807) = 5.33903689338684000D6*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(11)
        ckplm(k,3808) = 5.36022354772568000D6*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(12)
        ckplm(k,3809) = 5.42532747745514000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)*(rk+1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*rk * den(13)
        ckplm(k,3810) = 5.53930494546890000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-1.D0)&
                *(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*rk * den(14)
        ckplm(k,3811) = 5.71163887710571000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(15)
        ckplm(k,3812) = 5.95877709774971000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(16)
        ckplm(k,3813) = 6.30929339761734000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(17)
        ckplm(k,3814) = 6.81559471964836000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(18)
        ckplm(k,3815) = 7.58427081584930000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(19)
        ckplm(k,3816) = 8.87359685454369000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk+1.D0)*(rk-2.D0)*(rk+2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(20)
        ckplm(k,3817) = 1.15497609852791000D7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk+1.D0)*(rk-20.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(21)
        ckplm(k,3818) = 2.25745328348637000D7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-20.D0)*(rk-21.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*rk * den(22)
        ckplm(k,3819) = -1.96300285520554000D6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3820) = -45651.22919082640000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-43.00000000000000000D0+20.00000000000000000D0*rk) * den(1)
        ckplm(k,3821) = -70147.01070785520000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-1.D0)*(rk+20.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-41.00000000000000000D0+9.00000000000000000D0*rk) * den(2)
        ckplm(k,3822) = -29977.35500335690000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-117.00000000000000000D0+16.00000000000000000D0*rk) * den(3)
        ckplm(k,3823) = -53878.21912765500000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-74.00000000000000000D0+7.00000000000000000D0*rk) * den(4)
        ckplm(k,3824) = -24937.91856765750000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-175.00000000000000000D0+12.00000000000000000D0*rk) * den(5)
        ckplm(k,3825) = -47104.95729446410000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-99.00000000000000000D0+5.00000000000000000D0*rk) * den(6)
        ckplm(k,3826) = -22575.64773559570000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-217.00000000000000000D0+8.00000000000000000D0*rk) * den(7)
        ckplm(k,3827) = -43788.97190093990000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-116.00000000000000000D0+3.00000000000000000D0*rk) * den(8)
        ckplm(k,3828) = -21443.98212432860000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-243.00000000000000000D0+4.00000000000000000D0*rk) * den(9)
        ckplm(k,3829) = -42373.30867767330000000D0*(rk+10.D0)*(rk+11.D0)*(rk-125.D0)*(rk+12.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(10)
        ckplm(k,3830) = 5.33903689338684000D6*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(11)
        ckplm(k,3831) = 42373.30867767330000000D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk+126.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(12)
        ckplm(k,3832) = 21443.98212432860000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-1.D0)&
                *(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(247.00000000000000000D0+4.00000000000000000D0*rk) * den(13)
        ckplm(k,3833) = 43788.97190093990000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(119.00000000000000000D0+3.00000000000000000D0*rk) * den(14)
        ckplm(k,3834) = 22575.64773559570000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(225.00000000000000000D0+8.00000000000000000D0*rk) * den(15)
        ckplm(k,3835) = 47104.95729446410000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(104.00000000000000000D0+5.00000000000000000D0*rk) * den(16)
        ckplm(k,3836) = 24937.91856765750000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(187.00000000000000000D0+12.00000000000000000D0*rk) * den(17)
        ckplm(k,3837) = 53878.21912765500000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)*(rk+3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(81.00000000000000000D0+7.00000000000000000D0*rk) * den(18)
        ckplm(k,3838) = 29977.35500335690000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-1.D0)*(rk-2.D0)*(rk+2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(133.00000000000000000D0+16.00000000000000000D0*rk) * den(19)
        ckplm(k,3839) = 70147.01070785520000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-2.D0)&
                *(rk+2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(50.00000000000000000D0+9.00000000000000000D0*rk) * den(20)
        ckplm(k,3840) = 45651.22919082640000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-20.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(63.00000000000000000D0+20.00000000000000000D0*rk) * den(21)
        ckplm(k,3841) = 1.96300285520554000D6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-1.D0)*(rk-20.D0)*(rk-21.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3842) = 81791.78563356400000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+22.D0)*(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3843) = 3804.26909923554000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-86.00000000000000000D0+7.00000000000000000D0*rk) * den(1)
        ckplm(k,3844) = 139.18057680130000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(6642.00000000000000000D0+rk&
                *(-3187.00000000000000000D0+71.00000000000000000D0*rk)) * den(2)
        ckplm(k,3845) = 356.87327384948700000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk-2.D0)*(rk+3.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(4524.00000000000000000D0+(rk-1327.D0)*rk) * den(3)
        ckplm(k,3846) = -320.70368528366100000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-7252.00000000000000000D0+rk&
                *(1461.00000000000000000D0+19.00000000000000000D0*rk)) * den(4)
        ckplm(k,3847) = -98.95999431610110000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-30450.00000000000000000D0+rk&
                *(4441.00000000000000000D0+109.00000000000000000D0*rk)) * den(5)
        ckplm(k,3848) = -841.15995168685900000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-4334.00000000000000000D0+rk&
                *(467.00000000000000000D0+17.00000000000000000D0*rk)) * den(6)
        ckplm(k,3849) = -5643.91193389893000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-744.00000000000000000D0+rk&
                *(59.00000000000000000D0+3.00000000000000000D0*rk)) * den(7)
        ckplm(k,3850) = -608.18016529083300000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-7656.00000000000000000D0+rk&
                *(433.00000000000000000D0+31.00000000000000000D0*rk)) * den(8)
        ckplm(k,3851) = -255.28550148010300000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-19602.00000000000000000D0+rk&
                *(731.00000000000000000D0+79.00000000000000000D0*rk)) * den(9)
        ckplm(k,3852) = -252.22207546234100000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-20750.00000000000000000D0+rk&
                *(417.00000000000000000D0+83.00000000000000000D0*rk)) * den(10)
        ckplm(k,3853) = -21186.65433883670000000D0*(rk*rk+rk-252.D0)*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(11)
        ckplm(k,3854) = -252.22207546234100000D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-21084.00000000000000000D0+rk&
                *(-251.00000000000000000D0+83.00000000000000000D0*rk)) * den(12)
        ckplm(k,3855) = -255.28550148010300000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-20254.00000000000000000D0+rk&
                *(-573.00000000000000000D0+79.00000000000000000D0*rk)) * den(13)
        ckplm(k,3856) = -608.18016529083300000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(-8058.00000000000000000D0+rk&
                *(-371.00000000000000000D0+31.00000000000000000D0*rk)) * den(14)
        ckplm(k,3857) = -5643.91193389893000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-800.00000000000000000D0+rk&
                *(-53.00000000000000000D0+3.00000000000000000D0*rk)) * den(15)
        ckplm(k,3858) = -841.15995168685900000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-4784.00000000000000000D0+rk&
                *(-433.00000000000000000D0+17.00000000000000000D0*rk)) * den(16)
        ckplm(k,3859) = -98.95999431610110000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-34782.00000000000000000D0+rk&
                *(-4223.00000000000000000D0+109.00000000000000000D0*rk)) * den(17)
        ckplm(k,3860) = -320.70368528366100000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)*(rk+4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-8694.00000000000000000D0+rk&
                *(-1423.00000000000000000D0+19.00000000000000000D0*rk)) * den(18)
        ckplm(k,3861) = 356.87327384948700000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-2.D0)*(rk-3.D0)*(rk+3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(5852.00000000000000000D0+(rk+1329.D0)*rk) &
                * den(19)
        ckplm(k,3862) = 139.18057680130000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-2.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(9900.00000000000000000D0+rk&
                *(3329.00000000000000000D0+71.00000000000000000D0*rk)) * den(20)
        ckplm(k,3863) = 3804.26909923554000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-2.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(93.00000000000000000D0+7.00000000000000000D0*rk) * den(21)
        ckplm(k,3864) = 81791.78563356400000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)&
                *(rk-2.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3865) = -3271.67142534256000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+22.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3866) = -76.08538198471070000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-387.00000000000000000D0+4.00000000000000000D0*rk) * den(1)
        ckplm(k,3867) = 16.70166921615600000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+4.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-8979.00000000000000000D0+rk&
                *(1894.00000000000000000D0+33.00000000000000000D0*rk)) * den(2)
        ckplm(k,3868) = 1.42749309539795000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(521040.00000000000000000D0+rk*(-228593.00000000000000000D0+rk&
                *(17103.00000000000000000D0+656.00000000000000000D0*rk))) * den(3)
        ckplm(k,3869) = 2.56562948226929000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk-3.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(528360.00000000000000000D0+rk*(-159652.00000000000000000D0+rk&
                *(6045.00000000000000000D0+427.00000000000000000D0*rk))) * den(4)
        ckplm(k,3870) = 5.93759965896606000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(350070.00000000000000000D0+rk*(-76659.00000000000000000D0+rk&
                *(941.00000000000000000D0+188.00000000000000000D0*rk))) * den(5)
        ckplm(k,3871) = 33.64639806747440000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(84634.00000000000000000D0+rk*(-13699.00000000000000000D0+rk&
                *(-126.00000000000000000D0+31.00000000000000000D0*rk))) * den(6)
        ckplm(k,3872) = 112.87823867797900000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(31868.00000000000000000D0+rk*(-3797.00000000000000000D0+rk&
                *(-117.00000000000000000D0+8.00000000000000000D0*rk))) * den(7)
        ckplm(k,3873) = 14.59632396698000000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(292320.00000000000000000D0+rk*(-24842.00000000000000000D0+rk&
                *(-1423.00000000000000000D0+49.00000000000000000D0*rk))) * den(8)
        ckplm(k,3874) = 1.02114200592041000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk+9.D0)*(4.70529000000000000D6+rk*(-263659.00000000000000000D0+rk&
                *(-25995.00000000000000000D0+484.00000000000000000D0*rk))) * den(9)
        ckplm(k,3875) = 252.22207546234100000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(20502.00000000000000000D0+rk*(-619.00000000000000000D0+(rk-120.D0)*rk)) &
                * den(10)
        ckplm(k,3876) = -31779.98150825500000000D0*(rk*rk+rk-168.D0)*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0) * den(11)
        ckplm(k,3877) = -252.22207546234100000D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(-21000.00000000000000000D0+rk*(-376.00000000000000000D0+(rk+123.D0)*rk)) &
                * den(12)
        ckplm(k,3878) = -1.02114200592041000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)&
                *(rk-9.D0)*(rk+9.D0)*(-4.94247000000000000D6+rk*(-210217.00000000000000000D0+rk&
                *(27447.00000000000000000D0+484.00000000000000000D0*rk))) * den(13)
        ckplm(k,3879) = -14.59632396698000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-3.D0)&
                *(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(-315690.00000000000000000D0+rk*(-21849.00000000000000000D0+rk&
                *(1570.00000000000000000D0+49.00000000000000000D0*rk))) * den(14)
        ckplm(k,3880) = -112.87823867797900000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(-35540.00000000000000000D0+rk*(-3539.00000000000000000D0+rk&
                *(141.00000000000000000D0+8.00000000000000000D0*rk))) * den(15)
        ckplm(k,3881) = -33.64639806747440000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-98176.00000000000000000D0+rk*(-13354.00000000000000000D0+rk&
                *(219.00000000000000000D0+31.00000000000000000D0*rk))) * den(16)
        ckplm(k,3882) = -5.93759965896606000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-427482.00000000000000000D0+rk*(-77977.00000000000000000D0+rk&
                *(-377.00000000000000000D0+188.00000000000000000D0*rk))) * den(17)
        ckplm(k,3883) = -2.56562948226929000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-3.D0)*(rk-4.D0)*(rk+4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-693630.00000000000000000D0+rk*(-170461.00000000000000000D0+rk&
                *(-4764.00000000000000000D0+427.00000000000000000D0*rk))) * den(18)
        ckplm(k,3884) = -1.42749309539795000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-766080.00000000000000000D0+rk*(-260831.00000000000000000D0+rk&
                *(-15135.00000000000000000D0+656.00000000000000000D0*rk))) * den(19)
        ckplm(k,3885) = -16.70166921615600000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-10840.00000000000000000D0+rk&
                *(-1828.00000000000000000D0+33.00000000000000000D0*rk)) * den(20)
        ckplm(k,3886) = 76.08538198471070000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-3.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(391.00000000000000000D0+4.00000000000000000D0*rk) * den(21)
        ckplm(k,3887) = 3271.67142534256000000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)&
                *(rk-3.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3888) = 125.83351635932900000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+22.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3889) = -5.85272169113159000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+5.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(344.00000000000000000D0+5.00000000000000000D0*rk) &
                * den(1)
        ckplm(k,3890) = -0.21412396430969200D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+5.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-77736.00000000000000000D0+rk&
                *(6935.00000000000000000D0+281.00000000000000000D0*rk)) * den(2)
        ckplm(k,3891) = -1.42749309539795000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(80520.00000000000000000D0+rk*(-19834.00000000000000000D0+rk&
                *(279.00000000000000000D0+43.00000000000000000D0*rk))) * den(3)
        ckplm(k,3892) = -0.06751656532287600D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.17171600000000000D7+rk*(4.70612200000000000D6+rk*(-416357.00000000000000000D0+rk&
                *(-10342.00000000000000000D0+737.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3893) = -0.10416841506958000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk-4.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.37617200000000000D7+rk*(4.01585400000000000D6+rk*(-197059.00000000000000000D0+rk&
                *(-13346.00000000000000000D0+311.00000000000000000D0*rk)))) * den(5)
        ckplm(k,3894) = -0.88543152809143100D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.51090400000000000D6+rk*(542254.00000000000000000D0+rk*(-9779.00000000000000000D0+rk&
                *(-2026.00000000000000000D0+15.00000000000000000D0*rk)))) * den(6)
        ckplm(k,3895) = 0.45699691772460900D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(6.74014400000000000D6+rk*(-1.07217000000000000D6+rk*(-11535.00000000000000000D0+rk&
                *(4150.00000000000000000D0+11.00000000000000000D0*rk)))) * den(7)
        ckplm(k,3896) = 0.00984907150268555D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(3.96866160000000000D8+rk*(-4.50438660000000000D7+rk*(-1.93578700000000000D6+rk&
                *(174790.00000000000000000D0+2111.00000000000000000D0*rk)))) * den(8)
        ckplm(k,3897) = 0.00413417816162109D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.11573936000000000D9+rk*(-8.35158660000000000D7+rk*(-7.44631900000000000D6+rk&
                *(318486.00000000000000000D0+7939.00000000000000000D0*rk)))) * den(9)
        ckplm(k,3898) = 0.51057100296020500D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.00064640000000000D7+rk*(-403618.00000000000000000D0+rk*(-76195.00000000000000000D0+rk&
                *(1462.00000000000000000D0+79.00000000000000000D0*rk)))) * den(10)
        ckplm(k,3899) = 42.88796424865720000D0*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(124488.00000000000000000D0+(rk*rk+rk-990.D0)*(rk+1.D0)*rk) * den(11)
        ckplm(k,3900) = 0.51057100296020500D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.03325040000000000D7+rk*(247158.00000000000000000D0+rk*(-80107.00000000000000000D0+rk&
                *(-1146.00000000000000000D0+79.00000000000000000D0*rk)))) * den(12)
        ckplm(k,3901) = 0.00413417816162109D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.19149836000000000D9+rk*(6.76995260000000000D7+rk*(-8.35414300000000000D6+rk&
                *(-286730.00000000000000000D0+7939.00000000000000000D0*rk)))) * den(13)
        ckplm(k,3902) = 0.00984907150268555D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-4.D0)&
                *(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(4.39801560000000000D8+rk*(4.06563660000000000D7+rk*(-2.44749100000000000D6+rk&
                *(-166346.00000000000000000D0+2111.00000000000000000D0*rk)))) * den(14)
        ckplm(k,3903) = 0.45699691772460900D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(7.79664000000000000D6+rk*(1.03669400000000000D6+rk*(-23919.00000000000000000D0+rk&
                *(-4106.00000000000000000D0+11.00000000000000000D0*rk)))) * den(15)
        ckplm(k,3904) = -0.88543152809143100D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-3.06089600000000000D6+rk*(-555674.00000000000000000D0+rk*(-3611.00000000000000000D0+rk&
                *(2086.00000000000000000D0+15.00000000000000000D0*rk)))) * den(16)
        ckplm(k,3905) = -0.10416841506958000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-4.D0)*(rk-5.D0)*(rk+5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.79609760000000000D7+rk*(-4.36869000000000000D6+rk*(-155155.00000000000000000D0+rk&
                *(14590.00000000000000000D0+311.00000000000000000D0*rk)))) * den(17)
        ckplm(k,3906) = -0.06751656532287600D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.68285600000000000D7+rk*(-5.50486200000000000D6+rk*(-380909.00000000000000000D0+rk&
                *(13290.00000000000000000D0+737.00000000000000000D0*rk)))) * den(18)
        ckplm(k,3907) = -1.42749309539795000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-100590.00000000000000000D0+rk*(-20263.00000000000000000D0+rk&
                *(-150.00000000000000000D0+43.00000000000000000D0*rk))) * den(19)
        ckplm(k,3908) = -0.21412396430969200D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-84390.00000000000000000D0+rk&
                *(-6373.00000000000000000D0+281.00000000000000000D0*rk)) * den(20)
        ckplm(k,3909) = -5.85272169113159000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-4.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-339.00000000000000000D0+5.00000000000000000D0*rk) &
                * den(21)
        ckplm(k,3910) = 125.83351635932900000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)&
                *(rk-4.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3911) = -4.66050060590108000D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0)&
                *(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3912) = 0.10838373502095500D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+6.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(1075.00000000000000000D0+28.00000000000000000D0*rk) * den(1)
        ckplm(k,3913) = 0.00793051719665527D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+6.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(-181425.00000000000000000D0+rk&
                *(4706.00000000000000000D0+419.00000000000000000D0*rk)) * den(2)
        ckplm(k,3914) = 0.13217528661092100D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(101940.00000000000000000D0+rk*(-13273.00000000000000000D0+rk&
                *(-369.00000000000000000D0+16.00000000000000000D0*rk))) * den(3)
        ckplm(k,3915) = 0.01250306765238440D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-9.26184000000000000D6+rk*(2.30502800000000000D6+rk*(-63841.00000000000000000D0+rk&
                *(-7916.00000000000000000D0+49.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3916) = -0.00064301490783691D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.53855828000000000D9+rk*(5.59903866000000000D8+rk*(-4.80348250000000000D7+rk&
                *(-1.31971000000000000D6+rk*(161185.00000000000000000D0+1084.00000000000000000D0*rk))))) &
                * den(5)
        ckplm(k,3917) = -0.00364375114440918D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk-5.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-4.76281080000000000D8+rk*(1.28530686000000000D8+rk*(-5.82582500000000000D6+rk&
                *(-540685.00000000000000000D0+rk*(21185.00000000000000000D0+439.00000000000000000D0*rk))))) &
                * den(6)
        ckplm(k,3918) = -0.00094032287597656D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-2.80399464000000000D9+rk*(5.57980348000000000D8+rk*(-5.23732500000000000D6+rk&
                *(-2.84333000000000000D6+rk*(35055.00000000000000000D0+2152.00000000000000000D0*rk))))) * den(7)
        ckplm(k,3919) = -0.00060796737670898D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-5.88816000000000000D9+rk*(8.36545080000000000D8+rk*(2.39147540000000000D7+rk&
                *(-4.61782900000000000D6+rk*(-28358.00000000000000000D0+3265.00000000000000000D0*rk))))) &
                * den(8)
        ckplm(k,3920) = -0.00344514846801758D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-1.28516580000000000D9+rk*(1.20458790000000000D8+rk*(9.59662700000000000D6+rk&
                *(-684302.00000000000000000D0+rk*(-18167.00000000000000000D0+452.00000000000000000D0*rk))))) &
                * den(9)
        ckplm(k,3921) = -0.17019033432006800D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(-2.96577840000000000D7+5.00000000000000000D0*rk*(299650.00000000000000000D0+rk&
                *(55035.00000000000000000D0+rk*(-1679.00000000000000000D0+(rk-111.D0)*rk)))) * den(10)
        ckplm(k,3922) = 21.44398212432860000D0*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(248976.00000000000000000D0+5.00000000000000000D0*(rk*rk+rk-496.D0)*(rk+1.D0)*rk) * den(11)
        ckplm(k,3923) = 0.17019033432006800D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(3.08730240000000000D7+5.00000000000000000D0*rk*(184992.00000000000000000D0+rk&
                *(-59396.00000000000000000D0+rk*(-1225.00000000000000000D0+(rk+116.D0)*rk)))) * den(12)
        ckplm(k,3924) = 0.00344514846801758D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)&
                *(1.39536228000000000D9+rk*(9.92875580000000000D7+rk*(-1.15360110000000000D7+rk&
                *(-607114.00000000000000000D0+rk*(20427.00000000000000000D0+452.00000000000000000D0*rk))))) &
                * den(13)
        ckplm(k,3925) = 0.00060796737670898D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)&
                *(6.69620412000000000D9+rk*(7.74991842000000000D8+rk*(-3.75654430000000000D7+rk&
                *(-4.47174700000000000D6+rk*(44683.00000000000000000D0+3265.00000000000000000D0*rk))))) &
                * den(14)
        ckplm(k,3926) = 0.00094032287597656D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(3.36433608000000000D9+rk*(5.59795548000000000D8+rk*(-3.48147500000000000D6+rk&
                *(-2.96203000000000000D6+rk*(-24295.00000000000000000D0+2152.00000000000000000D0*rk))))) &
                * den(15)
        ckplm(k,3927) = 0.00364375114440918D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-5.D0)*(rk-6.D0)*(rk+6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(6.10076160000000000D8+rk*(1.38477736000000000D8+rk*(4.08105000000000000D6+rk&
                *(-621035.00000000000000000D0+rk*(-18990.00000000000000000D0+439.00000000000000000D0*rk))))) &
                * den(16)
        ckplm(k,3928) = 0.00064301490783691D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(2.14501716000000000D9+rk*(6.51375066000000000D8+rk*(4.31194250000000000D7+rk&
                *(-1.95361000000000000D6+rk*(-155765.00000000000000000D0+1084.00000000000000000D0*rk))))) &
                * den(17)
        ckplm(k,3929) = -0.01250306765238440D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-1.16227440000000000D7+rk*(-2.40876600000000000D6+rk*(-39799.00000000000000000D0+rk&
                *(8112.00000000000000000D0+49.00000000000000000D0*rk)))) * den(18)
        ckplm(k,3930) = -0.13217528661092100D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-114828.00000000000000000D0+rk*(-12487.00000000000000000D0+rk&
                *(417.00000000000000000D0+16.00000000000000000D0*rk))) * den(19)
        ckplm(k,3931) = -0.00793051719665527D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-5.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-185712.00000000000000000D0+rk*(-3868.00000000000000000D0+419.00000000000000000D0&
                *rk)) * den(20)
        ckplm(k,3932) = -0.10838373502095500D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-5.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(-1047.00000000000000000D0+28.00000000000000000D0*rk) * den(21)
        ckplm(k,3933) = 4.66050060590108000D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)*(rk-5.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3934) = 0.16644645021075300D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0)&
                *(rk+7.D0)*(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3935) = -0.00774169535863967D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+7.D0)&
                *(rk+8.D0)*(rk+9.D0)*(774.00000000000000000D0+25.00000000000000000D0*rk) * den(1)
        ckplm(k,3936) = -0.00594788789749146D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+7.D0)*(rk+8.D0)&
                *(rk+9.D0)*(-17466.00000000000000000D0+rk*(-185.00000000000000000D0+21.00000000000000000D0&
                *rk)) * den(2)
        ckplm(k,3937) = -0.01321752866109210D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(97080.00000000000000000D0+rk*(-5398.00000000000000000D0+(rk-417.D0)*rk)) * den(3)
        ckplm(k,3938) = 0.00008930762608846D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)&
                *(1.52633880000000000D8+rk*(-2.18356380000000000D7+rk*(-480407.00000000000000000D0+rk&
                *(61998.00000000000000000D0+767.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3939) = 0.00096452236175537D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(-1.41780240000000000D8+rk&
                *(3.33070680000000000D7+rk*(-875780.00000000000000000D0+rk*(-148535.00000000000000000D0+rk&
                *(2540.00000000000000000D0+107.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3940) = 0.00005358457565308D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(2.52851781600000000D10+rk&
                *(-8.17753472400000000D9+rk*(5.93641880000000000D8+rk*(2.99711850000000000D7+rk&
                *(-3.01217500000000000D6+rk*(-30981.00000000000000000D0+1775.00000000000000000D0*rk)))))) &
                * den(6)
        ckplm(k,3941) = 3.95093645368304000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk-6.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(5.71124070720000000D11+rk&
                *(-1.36412785128000000D11+rk*(4.03603140200000000D9+rk*(7.87530075000000000D8+rk&
                *(-2.76291850000000000D7+rk*(-1.14386700000000000D6+14543.00000000000000000D0*rk)))))) * den(7)
        ckplm(k,3942) = 7.66345432826451000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk+8.D0)*(rk+9.D0)*(4.27713812640000000D12+rk&
                *(-7.30018379160000000D11+rk*(-1.03148055780000000D10+rk*(5.02743432100000000D9+rk&
                *(-2.41027010000000000D7+rk*(-7.59633700000000000D6+8303.00000000000000000D0*rk)))))) * den(8)
        ckplm(k,3943) = -0.00002026557922363D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(-2.09682154800000000D11+rk&
                *(2.36230810200000000D10+rk*(1.65637103600000000D9+rk*(-1.75688205000000000D8+rk&
                *(-3.99654700000000000D6+rk*(262185.00000000000000000D0+2111.00000000000000000D0*rk)))))) &
                * den(9)
        ckplm(k,3944) = -0.00250279903411865D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-1.99234152000000000D9+rk&
                *(1.21011444000000000D8+rk*(2.16022840000000000D7+rk*(-915087.00000000000000000D0+rk&
                *(-63587.00000000000000000D0+rk*(1299.00000000000000000D0+31.00000000000000000D0*rk)))))) &
                * den(10)
        ckplm(k,3945) = -0.09010076522827150D0*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-5.92562880000000000D7+(rk+1.D0)*rk&
                *(709728.00000000000000000D0+(rk*rk+rk-2150.D0)*(rk+1.D0)*rk)) * den(11)
        ckplm(k,3946) = -0.00250279903411865D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-2.09090044800000000D9+rk&
                *(-7.53222720000000000D7+rk*(2.39534980000000000D7+rk*(648369.00000000000000000D0+rk&
                *(-69617.00000000000000000D0+rk*(-1113.00000000000000000D0+31.00000000000000000D0*rk)))))) &
                * den(12)
        ckplm(k,3947) = -0.00002026557922363D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-2.31477433200000000D11+rk&
                *(-1.98005587800000000D10+rk*(2.15686618400000000D9+rk*(1.57122387000000000D8+rk&
                *(-5.27580700000000000D6+rk*(-249519.00000000000000000D0+2111.00000000000000000D0*rk)))))) &
                * den(13)
        ckplm(k,3948) = 7.66345432826451000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(4.99179776760000000D12+rk&
                *(6.94248085740000000D11+rk*(-2.54656368320000000D10+rk*(-5.04771569500000000D9+rk&
                *(1.40035290000000000D7+rk*(7.64615500000000000D6+8303.00000000000000000D0*rk)))))) * den(14)
        ckplm(k,3949) = 3.95093645368304000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-6.D0)*(rk-7.D0)*(rk+7.D0)*(rk-8.D0)*(rk-9.D0)*(7.10758886400000000D11+rk&
                *(1.42017547560000000D11+rk*(1.51932288200000000D9+rk*(-8.86317285000000000D8+rk&
                *(-2.16917050000000000D7+rk*(1.23112500000000000D6+14543.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,3950) = 0.00005358457565308D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(3.40234041600000000D10+rk&
                *(9.26302178400000000D9+rk*(4.85991710000000000D8+rk*(-4.16745750000000000D7+rk&
                *(-2.83064500000000000D6+rk*(41631.00000000000000000D0+1775.00000000000000000D0*rk)))))) &
                * den(16)
        ckplm(k,3951) = 0.00096452236175537D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.75812120000000000D8+rk&
                *(3.46033980000000000D7+rk*(416005.00000000000000000D0+rk*(-157625.00000000000000000D0+rk&
                *(-2005.00000000000000000D0+107.00000000000000000D0*rk))))) * den(17)
        ckplm(k,3952) = 0.00008930762608846D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(1.73927880000000000D8+rk*(2.06918980000000000D7+rk*(-661799.00000000000000000D0+rk&
                *(-58930.00000000000000000D0+767.00000000000000000D0*rk)))) * den(18)
        ckplm(k,3953) = -0.01321752866109210D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-102060.00000000000000000D0+rk*(-4561.00000000000000000D0+(rk+420.D0)*rk)) * den(19)
        ckplm(k,3954) = -0.00594788789749146D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-6.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-17260.00000000000000000D0+rk*(227.00000000000000000D0+21.00000000000000000D0*rk)) * den(20)
        ckplm(k,3955) = -0.00774169535863967D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-6.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0)*(-749.00000000000000000D0+25.00000000000000000D0*rk) * den(21)
        ckplm(k,3956) = 0.16644645021075300D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)*(rk-6.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3957) = -0.00573953276588803D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0)&
                *(rk+8.D0)*(rk+9.D0) * den(0)
        ckplm(k,3958) = 0.00013347750618344D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+8.D0)&
                *(rk+9.D0)*(2107.00000000000000000D0+76.00000000000000000D0*rk) * den(1)
        ckplm(k,3959) = 0.00006836652755737D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-95571.00000000000000000D0+rk*(-3194.00000000000000000D0+37.00000000000000000D0*rk)) * den(2)
        ckplm(k,3960) = -0.00022788842519124D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+8.D0)*(rk+9.D0)&
                *(-455070.00000000000000000D0+rk*(2447.00000000000000000D0+rk&
                *(1443.00000000000000000D0+16.00000000000000000D0*rk))) * den(3)
        ckplm(k,3961) = -0.00002155701319377D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+8.D0)*(rk+9.D0)*(6.24552600000000000D7+rk&
                *(-4.18081600000000000D6+rk*(-347533.00000000000000000D0+rk&
                *(5956.00000000000000000D0+253.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3962) = -5.54323196411133000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+8.D0)*(rk+9.D0)*(-2.84807124000000000D9+rk&
                *(3.92571498000000000D8+rk*(1.15719250000000000D7+rk*(-1.66985000000000000D6+rk&
                *(-26065.00000000000000000D0+692.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3963) = -1.84774398803711000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+8.D0)*(rk+9.D0)*(9.54584492400000000D10+rk&
                *(-2.00258264460000000D10+rk*(2.47032257000000000D8+rk*(1.17054900000000000D8+rk&
                *(-1.86073000000000000D6+rk*(-164454.00000000000000000D0+353.00000000000000000D0*rk)))))) &
                * den(6)
        ckplm(k,3964) = 6.81195940290179000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+8.D0)*(rk+9.D0)*(2.83459294162800000D13+rk*(-7.89634296961200000D12+rk&
                *(3.93311378572000000D11+rk*(4.65454500490000000D10+rk*(-3.09696684500000000D9+rk&
                *(-8.19684530000000000D7+rk*(4.05271300000000000D6+34936.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,3965) = 2.55448477608817000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk-7.D0)*(rk+8.D0)*(rk+9.D0)*(1.17463988376000000D13+rk*(-2.34099509904000000D12+rk&
                *(9.05373588000000000D8+rk*(1.86663861160000000D10+rk*(-3.47803813000000000D8+rk&
                *(-4.17168990000000000D7+rk*(481705.00000000000000000D0+16103.00000000000000000D0*rk))))))) &
                * den(8)
        ckplm(k,3966) = 0.00001013278961182D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk+9.D0)*(4.02425679600000000D11+rk*(-5.29751327400000000D10+rk&
                *(-3.20558599200000000D9+rk*(4.80970889000000000D8+rk*(8.00349900000000000D6+rk&
                *(-1.13035300000000000D6+rk*(-7707.00000000000000000D0+404.00000000000000000D0*rk))))))) &
                * den(9)
        ckplm(k,3967) = 0.00250279903411865D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(1.96815528000000000D9+rk*(-1.39731516000000000D8+rk&
                *(-2.42196920000000000D7+rk*(1.32712300000000000D6+rk*(92302.00000000000000000D0+rk&
                *(-3080.00000000000000000D0+(rk-98.D0)*rk)))))) * den(10)
        ckplm(k,3968) = -0.31535267829895000D0*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-1.69303680000000000D7+(rk+1.D0)*rk&
                *(237060.00000000000000000D0+(rk*rk+rk-960.D0)*(rk+1.D0)*rk)) * den(11)
        ckplm(k,3969) = -0.00250279903411865D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-2.08243526400000000D9+rk*(-8.76947760000000000D7+rk&
                *(2.76179400000000000D7+rk*(929110.00000000000000000D0+rk*(-106197.00000000000000000D0+rk&
                *(-2471.00000000000000000D0+(rk+105.D0)*rk)))))) * den(12)
        ckplm(k,3970) = -0.00001013278961182D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk+8.D0)*(rk-9.D0)*(rk+9.D0)*(-4.51723381200000000D11+rk*(-4.51586647800000000D10+rk&
                *(4.58929822400000000D9+rk*(4.37821643000000000D8+rk*(-1.35255190000000000D7+rk&
                *(-1.07562700000000000D6+rk*(10535.00000000000000000D0+404.00000000000000000D0*rk))))))) &
                * den(13)
        ckplm(k,3971) = -2.55448477608817000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk+8.D0)*(rk-9.D0)*(-1.40693273028000000D13+rk*(-2.28562683462000000D12+rk&
                *(5.67565512360000000D10+rk*(1.96313618830000000D10+rk*(1.32557348000000000D8+rk&
                *(-4.42689660000000000D7+rk*(-368984.00000000000000000D0+16103.00000000000000000D0*rk))))))) &
                * den(14)
        ckplm(k,3972) = 6.81195940290179000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(3.65860273338000000D13+rk*(8.53137542322000000D12+rk&
                *(2.35972968924000000D11+rk*(-5.80338013990000000D10+rk*(-2.62755664500000000D9+rk&
                *(1.05551075000000000D8+(3.80816100000000000D6-34936.00000000000000000D0*rk)*rk)))))) * den(15)
        ckplm(k,3973) = 1.84774398803711000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(1.15612557120000000D11+rk*(2.01621077280000000D10+rk&
                *(-1.13646988000000000D8+rk*(-1.22846220000000000D8+rk*(-1.03316500000000000D6+rk&
                *(166572.00000000000000000D0+353.00000000000000000D0*rk)))))) * den(16)
        ckplm(k,3974) = 5.54323196411133000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(3.22742772000000000D9+rk&
                *(3.64525818000000000D8+rk*(-1.64181650000000000D7+rk*(-1.55867000000000000D6+rk&
                *(29525.00000000000000000D0+692.00000000000000000D0*rk))))) * den(17)
        ckplm(k,3975) = 0.00002155701319377D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)*(6.62828400000000000D7+rk&
                *(3.46889400000000000D6+rk*(-363883.00000000000000000D0+rk&
                *(-4944.00000000000000000D0+253.00000000000000000D0*rk)))) * den(18)
        ckplm(k,3976) = 0.00022788842519124D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(456090.00000000000000000D0+rk*(-391.00000000000000000D0+rk&
                *(-1395.00000000000000000D0+16.00000000000000000D0*rk))) * den(19)
        ckplm(k,3977) = -0.00006836652755737D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-7.D0)*(rk-8.D0)*(rk-9.D0)&
                *(-92340.00000000000000000D0+rk*(3268.00000000000000000D0+37.00000000000000000D0*rk)) * den(20)
        ckplm(k,3978) = -0.00013347750618344D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-7.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-2031.00000000000000000D0+76.00000000000000000D0*rk) * den(21)
        ckplm(k,3979) = 0.00573953276588803D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)*(rk-7.D0)&
                *(rk-8.D0)*(rk-9.D0) * den(22)
        ckplm(k,3980) = 0.00019131775886293D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0)&
                *(rk+9.D0) * den(0)
        ckplm(k,3981) = -8.89850041222951000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+9.D0)*(1376.00000000000000000D0+53.00000000000000000D0*rk) * den(1)
        ckplm(k,3982) = 6.83665275573731000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+9.D0)&
                *(53792.00000000000000000D0+rk*(2613.00000000000000000D0+11.00000000000000000D0*rk)) * den(2)
        ckplm(k,3983) = 0.00001519256167942D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+9.D0)&
                *(-482400.00000000000000000D0+rk*(-14362.00000000000000000D0+rk&
                *(807.00000000000000000D0+19.00000000000000000D0*rk))) * den(3)
        ckplm(k,3984) = 1.02652443779839000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+9.D0)*(1.11950160000000000D9+rk&
                *(-1.19006660000000000D7+rk*(-5.44292300000000000D6+rk&
                *(-63034.00000000000000000D0+1823.00000000000000000D0*rk)))) * den(4)
        ckplm(k,3985) = -3.69548797607422000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+9.D0)*(4.24001088000000000D9+rk&
                *(-2.66569752000000000D8+rk*(-2.96201500000000000D7+rk*(682025.00000000000000000D0+rk&
                *(46210.00000000000000000D0+67.00000000000000000D0*rk))))) * den(5)
        ckplm(k,3986) = -1.84774398803711000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+9.D0)*(-1.06601616000000000D11+rk&
                *(1.27192442160000000D10+rk*(6.82136678000000000D8+rk*(-7.15724750000000000D7+rk&
                *(-1.85034500000000000D6+rk*(70379.00000000000000000D0+987.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,3987) = -1.36239188058036000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+9.D0)*(1.73261754288000000D13+rk*(-3.03615406433600000D12+rk&
                *(-4.78432335640000000D10+rk*(2.32328144720000000D10+rk*(-6.57555850000000000D7+rk&
                *(-4.88391190000000000D7+rk*(-28651.00000000000000000D0+15583.00000000000000000D0*rk))))))) &
                * den(7)
        ckplm(k,3988) = -8.51494925362723000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+9.D0)*(-3.22540091899200000D14+rk*(7.35049814968800000D13+rk&
                *(-1.09881887847600000D12+rk*(-6.45032607604000000D11+rk*(2.25645103670000000D10+rk&
                *(1.80412484000000000D9+rk*(-5.15366740000000000D7+rk&
                *(-1.52515600000000000D6+14543.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,3989) = 2.02655792236328000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk-8.D0)&
                *(rk+9.D0)*(1.93060516320000000D12+rk*(-2.90860667280000000D11+rk*(-1.47913290840000000D10+rk&
                *(3.08420118000000000D9+rk*(2.83851790000000000D7+rk*(-9.82220000000000000D6+rk&
                *(-17706.00000000000000000D0+rk*(8300.00000000000000000D0+11.00000000000000000D0*rk)))))))) &
                * den(9)
        ckplm(k,3990) = 0.00005005598068237D0*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(9.72085296000000000D10+rk*(-7.90216808000000000D9+rk*(-1.32844623600000000D9+rk&
                *(9.01588440000000000D7+rk*(6.12830700000000000D6+rk*(-293400.00000000000000000D0+rk&
                *(-10074.00000000000000000D0+rk*(236.00000000000000000D0+3.00000000000000000D0*rk)))))))) &
                * den(10)
        ckplm(k,3991) = 0.00020022392272949D0*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(2.66653296000000000D10+(rk+1.D0)*rk*(-4.27586544000000000D8+(rk+1.D0)*rk&
                *(2.17090800000000000D6+(rk*rk+rk-3620.D0)*(rk+1.D0)*rk))) * den(11)
        ckplm(k,3992) = 0.00005005598068237D0*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(1.03698504000000000D11+rk*(5.00071723200000000D9+rk*(-1.55937490800000000D9+rk&
                *(-6.29211880000000000D7+rk*(7.43614700000000000D6+rk*(228168.00000000000000000D0+rk&
                *(-11642.00000000000000000D0+rk*(-212.00000000000000000D0+3.00000000000000000D0*rk)))))))) &
                * den(12)
        ckplm(k,3993) = 2.02655792236328000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-8.D0)*(rk-9.D0)&
                *(rk+9.D0)*(2.20362848160000000D12+rk*(2.52187893040000000D11+rk*(-2.37758391320000000D10+rk&
                *(-2.87308246800000000D9+rk*(7.69408590000000000D7+rk*(9.54228000000000000D6+rk&
                *(-75498.00000000000000000D0+rk*(-8212.00000000000000000D0+11.00000000000000000D0*rk)))))))) &
                * den(13)
        ckplm(k,3994) = 8.51494925362723000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-8.D0)&
                *(rk-9.D0)*(3.96478149278400000D14-1.00000000000000000D0*rk*(-7.36865824413600000D13+rk&
                *(9.52884143508000000D11+rk*(7.16272862060000000D11+rk*(1.28252345270000000D10+rk&
                *(-2.08050220000000000D9+rk*(-4.04533780000000000D7+rk&
                *(1.64150000000000000D6+14543.00000000000000000D0*rk)))))))) * den(14)
        ckplm(k,3995) = 1.36239188058036000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-8.D0)*(rk-9.D0)*(2.02912364844000000D13+rk*(2.87075004606000000D12-1.00000000000000000D0&
                *rk*(1.17448576308000000D11+rk*(2.30085640470000000D10+rk*(-1.77464840000000000D8+rk&
                *(-4.83399700000000000D7+rk*(137732.00000000000000000D0+15583.00000000000000000D0*rk))))))) &
                * den(15)
        ckplm(k,3996) = -1.84774398803711000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-8.D0)*(rk-9.D0)*(-1.18569070800000000D11+rk&
                *(-1.11480007880000000D10+rk*(8.85063048000000000D8+rk*(6.34870450000000000D7+rk&
                *(-2.18743500000000000D6+rk*(-64457.00000000000000000D0+987.00000000000000000D0*rk)))))) &
                * den(16)
        ckplm(k,3997) = -3.69548797607422000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-8.D0)*(rk-9.D0)*(-4.47632460000000000D9+rk&
                *(-2.05467882000000000D8+rk*(3.13896350000000000D7+rk*(497855.00000000000000000D0+rk&
                *(-45875.00000000000000000D0+67.00000000000000000D0*rk))))) * den(17)
        ckplm(k,3998) = 1.02652443779839000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-8.D0)*(rk-9.D0)*(1.12602420000000000D9+rk&
                *(1.21121400000000000D6+rk*(-5.24288300000000000D6+rk&
                *(70326.00000000000000000D0+1823.00000000000000000D0*rk)))) * den(18)
        ckplm(k,3999) = 0.00001519256167942D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-8.D0)*(rk-9.D0)*(467250.00000000000000000D0+rk&
                *(-15919.00000000000000000D0+rk*(-750.00000000000000000D0+19.00000000000000000D0*rk))) * den(19)
        ckplm(k,4000) = 6.83665275573731000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-8.D0)*(rk-9.D0)&
                *(51190.00000000000000000D0+rk*(-2591.00000000000000000D0+11.00000000000000000D0*rk)) * den(20)
        ckplm(k,4001) = -8.89850041222951000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-8.D0)&
                *(rk-9.D0)*(-1323.00000000000000000D0+53.00000000000000000D0*rk) * den(21)
        ckplm(k,4002) = 0.00019131775886293D0*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)*(rk-8.D0)&
                *(rk-9.D0) * den(22)
        ckplm(k,4003) = -6.17154060848176000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+22.D0) * den(0)
        ckplm(k,4004) = 1.43524200197250000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(3483.00000000000000000D0+140.00000000000000000D0*rk) * den(1)
        ckplm(k,4005) = -3.15053122384208000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(598887.00000000000000000D0+rk*(35410.00000000000000000D0+383.00000000000000000D0*rk)) * den(2)
        ckplm(k,4006) = -3.50059024871342000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(-1.32178500000000000D7+rk&
                *(-729479.00000000000000000D0+rk*(3189.00000000000000000D0+368.00000000000000000D0*rk))) &
                * den(3)
        ckplm(k,4007) = 3.31136915418837000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(-2.61608130000000000D9+rk&
                *(-8.36751360000000000D7+rk*(7.78185100000000000D6+rk&
                *(268836.00000000000000000D0+149.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4008) = 7.66345432826451000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(1.78201522800000000D10+rk*(-6.86789820000000000D7+rk&
                *(-1.12129355000000000D8+rk*(-1.66761000000000000D6+rk&
                *(98015.00000000000000000D0+1332.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4009) = 8.51494925362723000D-10*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(-2.26286931420000000D12+rk*(1.03908297798000000D11+rk&
                *(1.97962515550000000D10+rk*(-3.16648020000000000D8+rk*(-4.49643500000000000D7+rk&
                *(-143538.00000000000000000D0+11635.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,4010) = 6.81195940290179000D-9*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(3.70716984900000000D12+rk*(-3.31768528236000000D11+rk*(-3.56017251000000000D10+rk&
                *(2.29120594100000000D9+rk*(1.23383655000000000D8+rk*(-3.38824900000000000D6+rk&
                *(-121155.00000000000000000D0+344.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4011) = -7.66345432826451000D-8*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(rk+14.D0)*(4.09969985040000000D12+rk*(-5.38914202560000000D11+rk*(-3.63023885880000000D10+rk&
                *(5.22626152800000000D9+rk*(1.25455951000000000D8+rk*(-1.48816800000000000D7+rk&
                *(-213362.00000000000000000D0+rk*(10392.00000000000000000D0+79.00000000000000000D0*rk)))))))) &
                * den(8)
        ckplm(k,4012) = -1.44754137311663000D-7*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(-2.59301883456000000D13+rk*(4.40058470184000000D12+rk*(1.81071847932000000D11+rk&
                *(-5.28403257320000000D10+rk*(-4.36238670000000000D7+rk*(2.09832252000000000D8+rk&
                *(-1.15210200000000000D6+rk*(-287628.00000000000000000D0+rk&
                *(837.00000000000000000D0+68.00000000000000000D0*rk))))))))) * den(9)
        ckplm(k,4013) = -7.15085438319615000D-6*(rk+10.D0)*(rk+11.D0)*(rk+12.D0)*(rk-9.D0)&
                *(-6.72135609600000000D11+rk*(6.15822040800000000D10+rk*(1.00291879560000000D10+rk&
                *(-8.18812576000000000D8+rk*(-5.35833690000000000D7+rk*(3.43956900000000000D6+rk&
                *(117894.00000000000000000D0+rk*(-4674.00000000000000000D0+(rk-81.D0)*rk)))))))) * den(10)
        ckplm(k,4014) = 0.00090100765228271D0*(rk-10.D0)*(rk+10.D0)*(rk+11.D0)*(rk-9.D0)&
                *(5.92562880000000000D9+(rk+1.D0)*rk*(-1.07118144000000000D8+(rk+1.D0)*rk&
                *(654708.00000000000000000D0+(rk*rk+rk-1520.D0)*(rk+1.D0)*rk))) * den(11)
        ckplm(k,4015) = 7.15085438319615000D-6*(rk-10.D0)*(rk+10.D0)*(rk-11.D0)*(rk-9.D0)&
                *(7.22926713600000000D11+rk*(3.92981823360000000D10+rk*(-1.21315940400000000D10+rk&
                *(-5.72600260000000000D8+rk*(6.88550100000000000D7+rk*(2.63871300000000000D6+rk&
                *(-148260.00000000000000000D0+rk*(-3990.00000000000000000D0+(rk+90.D0)*rk)))))))) * den(12)
        ckplm(k,4016) = 1.44754137311663000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-9.D0)&
                *(3.00971151936000000D13+rk*(3.88114857864000000D12+rk*(-3.37221539052000000D11+rk&
                *(-5.05545738440000000D10+rk*(1.09994965500000000D9+rk*(2.10666372000000000D8+rk&
                *(-879018.00000000000000000D0+rk*(-291876.00000000000000000D0+rk&
                *(-225.00000000000000000D0+68.00000000000000000D0*rk))))))))) * den(13)
        ckplm(k,4017) = 7.66345432826451000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-9.D0)&
                *(4.59722551680000000D12+rk*(4.51205520720000000D11+rk*(-5.10830371160000000D10+rk&
                *(-4.58024746000000000D9+rk*(1.96305731000000000D8+rk*(1.33877000000000000D7+rk&
                *(-283894.00000000000000000D0+rk*(-9760.00000000000000000D0+79.00000000000000000D0*rk)))))))) &
                * den(14)
        ckplm(k,4018) = 6.81195940290179000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-9.D0)*(4.00117209660000000D12+rk*(2.54201206740000000D11-1.00000000000000000D0*rk&
                *(4.17029830520000000D10+rk*(1.76622397100000000D9+rk*(-1.38495535000000000D8+rk&
                *(-2.65409500000000000D6+rk*(123563.00000000000000000D0+344.00000000000000000D0*rk))))))) &
                * den(15)
        ckplm(k,4019) = 8.51494925362723000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-9.D0)*(2.34670952160000000D12-1.00000000000000000D0*rk&
                *(-6.35449205280000000D10+rk*(2.04780194200000000D10+rk*(1.38458700000000000D8+rk&
                *(-4.40721350000000000D7+rk*(213348.00000000000000000D0+11635.00000000000000000D0*rk)))))) &
                * den(16)
        ckplm(k,4020) = -7.66345432826451000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-9.D0)*(-1.77784662000000000D10+rk&
                *(1.50191498000000000D8+rk*(1.06551755000000000D8+rk*(-2.04635000000000000D6+rk&
                *(-91355.00000000000000000D0+1332.00000000000000000D0*rk))))) * den(17)
        ckplm(k,4021) = 3.31136915418837000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-9.D0)*(2.52489300000000000D9-1.00000000000000000D0*rk&
                *(9.84329260000000000D7+rk*(6.97623700000000000D6+rk&
                *(-268240.00000000000000000D0+149.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4022) = 3.50059024871342000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-9.D0)*(1.24855500000000000D7+rk&
                *(-734753.00000000000000000D0+rk*(-2085.00000000000000000D0+368.00000000000000000D0*rk))) &
                * den(19)
        ckplm(k,4023) = 3.15053122384208000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-9.D0)&
                *(563860.00000000000000000D0+rk*(-34644.00000000000000000D0+383.00000000000000000D0*rk)) &
                * den(20)
        ckplm(k,4024) = 1.43524200197250000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-9.D0)&
                *(3343.00000000000000000D0-140.00000000000000000D0*rk) * den(21)
        ckplm(k,4025) = 6.17154060848176000D-6*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0)*(rk-9.D0) &
                * den(22)
        ckplm(k,4026) = 1.92860644015055000D-7*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4027) = -8.97026251232813000D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(2150.00000000000000000D0+89.00000000000000000D0*rk) * den(1)
        ckplm(k,4028) = 9.84541007450649000D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(908150.00000000000000000D0+rk*(60623.00000000000000000D0+877.00000000000000000D0*rk)) * den(2)
        ckplm(k,4029) = 8.41488040556110000D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(-3.15104400000000000D7+rk&
                *(-2.33461400000000000D6+rk*(-31713.00000000000000000D0+353.00000000000000000D0*rk))) * den(3)
        ckplm(k,4030) = -3.98001100263025000D-11*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(-1.47310186800000000D10+rk&
                *(-9.55589978000000000D8+rk*(1.10048470000000000D7+rk&
                *(1.32163400000000000D6+12137.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4031) = 4.09372560270540000D-12*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(-2.59517667528000000D12+rk*(-1.12441049508000000D11+rk&
                *(1.00774271000000000D10+rk*(4.97092505000000000D8-1.00000000000000000D0*rk&
                *(464180.00000000000000000D0+129317.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4032) = -2.04686280135270000D-12*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(rk+16.D0)*(-8.19671528412000000D13+rk*(-1.18745462118000000D12+rk&
                *(6.06554579304000000D11+rk*(1.60862807750000000D10+rk*(-8.05263065000000000D8+rk&
                *(-2.26265150000000000D7+18521.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,4033) = 4.25747462681362000D-10*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(rk+15.D0)*(-5.64197160960000000D12+rk*(9.77047765600000000D10+rk*(5.91845440120000000D10+rk&
                *(1.00249288000000000D8+rk*(-1.69130795000000000D8+rk*(-2.09847500000000000D6+rk&
                *(95623.00000000000000000D0+987.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4034) = 1.33046082087926000D-9*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(2.40099914880000000D13+rk*(-1.17246441240000000D12+rk*(-3.08193088380000000D11+rk&
                *(8.54221949200000000D9+rk*(1.32264839500000000D9+rk*(-1.16703200000000000D7+rk&
                *(-1.83173000000000000D6+rk*(-4292.00000000000000000D0+355.00000000000000000D0*rk)))))))) &
                * den(8)
        ckplm(k,4035) = 4.52356679098947000D-8*(rk+11.D0)*(rk+12.D0)*(rk+13.D0)&
                *(-8.84395814400000000D12+rk*(6.86950783200000000D11+rk*(1.28745543240000000D11+rk&
                *(-7.88027558800000000D9+rk*(-6.90905474000000000D8+rk*(2.73452270000000000D7+rk&
                *(1.52194000000000000D6+rk*(-27562.00000000000000000D0+rk&
                *(-986.00000000000000000D0+3.00000000000000000D0*rk))))))))) * den(9)
        ckplm(k,4036) = -1.71895538057600000D-8*(rk+11.D0)*(rk+12.D0)*(-2.76174934963200000D14+rk&
                *(2.81664798369600000D13+rk*(4.43821814639200000D12+rk*(-4.26929760100000000D11+rk&
                *(-2.65926068900000000D10+rk*(2.19879145500000000D9+rk*(7.20115410000000000D7+rk&
                *(-4.29465000000000000D6+rk*(-79060.00000000000000000D0+rk&
                *(2335.00000000000000000D0+17.00000000000000000D0*rk)))))))))) * den(10)
        ckplm(k,4037) = -4.81307506561279000D-7*(rk-10.D0)*(rk+11.D0)&
                *(-1.10927771136000000D13+(rk+1.D0)*rk*(2.23270447680000000D11+(rk+1.D0)*rk&
                *(-1.59737990400000000D9+(rk+1.D0)*rk*(4.78450800000000000D6+(rk*rk+rk-5240.D0)*(rk+1.D0)&
                *rk)))) * den(11)
        ckplm(k,4038) = -1.71895538057600000D-8*(rk-10.D0)*(rk-11.D0)*(-2.99504982067200000D14+rk&
                *(-1.81261571702400000D13+rk*(5.53863193459200000D12+rk*(3.00157340100000000D11+rk&
                *(-3.63619031400000000D10+rk*(-1.68125184500000000D9+rk*(9.96678410000000000D7+rk&
                *(3.58015000000000000D6+rk*(-99310.00000000000000000D0+rk&
                *(-2165.00000000000000000D0+17.00000000000000000D0*rk)))))))))) * den(12)
        ckplm(k,4039) = 4.52356679098947000D-8*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)&
                *(9.39499981056000000D12+rk*(4.08709901328000000D11+rk*(-1.47990865076000000D11+rk&
                *(-4.87454942400000000D9+rk*(8.03907237000000000D8+rk*(1.76903790000000000D7+rk&
                *(-1.68701400000000000D6+rk*(-19566.00000000000000000D0+rk&
                *(1013.00000000000000000D0+3.00000000000000000D0*rk))))))))) * den(13)
        ckplm(k,4040) = 1.33046082087926000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(2.48670530841600000D13+rk*(5.35789564848000000D11+rk*(-3.25794529164000000D11+rk&
                *(-3.17138721200000000D9+rk*(1.35369911500000000D9+rk*(789952.00000000000000000D0+rk&
                *(-1.79174600000000000D6+rk*(7132.00000000000000000D0+355.00000000000000000D0*rk)))))))) &
                * den(14)
        ckplm(k,4041) = 4.25747462681362000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(5.68075902912000000D12+rk*(-1.96980996240000000D10+rk*(-5.78914097460000000D10+rk&
                *(7.53909803000000000D8+rk*(1.57238620000000000D8+rk*(-2.65148600000000000D6+rk&
                *(-88714.00000000000000000D0+987.00000000000000000D0*rk))))))) * den(15)
        ckplm(k,4042) = 2.04686280135270000D-12*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(8.01900125395200000D13-1.00000000000000000D0*rk&
                *(2.34919712890400000D12+rk*(5.53690701554000000D11+rk*(-1.90806974650000000D10+rk&
                *(-6.91852675000000000D8+rk*(2.27376410000000000D7+18521.00000000000000000D0*rk)))))) * den(16)
        ckplm(k,4043) = 4.09372560270540000D-12*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(2.47315562604000000D12-1.00000000000000000D0*rk&
                *(1.31103416058000000D11+rk*(8.58465767500000000D9+rk*(-4.97656055000000000D8+rk&
                *(182405.00000000000000000D0+129317.00000000000000000D0*rk))))) * den(17)
        ckplm(k,4044) = -3.98001100263025000D-11*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(-1.37657333520000000D10+rk&
                *(9.73683318000000000D8+rk*(7.11276700000000000D6+rk&
                *(-1.27308600000000000D6+12137.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4045) = 8.41488040556110000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(2.92078920000000000D7+rk&
                *(-2.27012900000000000D6+rk*(32772.00000000000000000D0+353.00000000000000000D0*rk))) * den(19)
        ckplm(k,4046) = 9.84541007450649000D-10*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(848404.00000000000000000D0+rk*(-58869.00000000000000000D0+877.00000000000000000D0*rk)) &
                * den(20)
        ckplm(k,4047) = -8.97026251232813000D-9*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)&
                *(-2061.00000000000000000D0+89.00000000000000000D0*rk) * den(21)
        ckplm(k,4048) = 1.92860644015055000D-7*(rk-10.D0)*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4049) = -5.84426193985015000D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4050) = 1.49504375205469000D-9*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(473.00000000000000000D0+20.00000000000000000D0*rk) * den(1)
        ckplm(k,4051) = -1.09393445272294000D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(361251.00000000000000000D0+rk&
                *(26170.00000000000000000D0+439.00000000000000000D0*rk)) * den(2)
        ckplm(k,4052) = 2.80496013518703000D-11*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(4.97725800000000000D7+rk&
                *(4.40429300000000000D6+rk*(102477.00000000000000000D0+304.00000000000000000D0*rk))) * den(3)
        ckplm(k,4053) = 2.65334066842017000D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(-1.36013539200000000D10+rk*(-1.23340733200000000D9+rk&
                *(-1.98342430000000000D7+rk*(629692.00000000000000000D0+12403.00000000000000000D0*rk)))) &
                * den(4)
        ckplm(k,4054) = 6.82287600450900000D-13*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(rk+17.D0)*(1.09525719996000000D12+rk*(8.93068893060000000D10+rk&
                *(-8.23595665000000000D8+rk*(-1.94222350000000000D8+rk&
                *(-3.19053500000000000D6+9244.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4055) = -2.50172120165330000D-12*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(rk+16.D0)*(5.26427650908000000D12+rk*(3.40992354378000000D11+rk*(-2.13164770250000000D10+rk&
                *(-1.70598744000000000D9+rk*(-5.65361000000000000D6+rk&
                *(1.01146200000000000D6+9235.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,4056) = -2.60179004971943000D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(-7.90648153920000000D11+rk*(-3.48901416720000000D10+rk*(6.13515085000000000D9+rk&
                *(3.15804457000000000D8+rk*(-8.74239500000000000D6+rk*(-550673.00000000000000000D0+rk&
                *(-1535.00000000000000000D0+88.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4057) = 9.75671268644787000D-10*(rk+12.D0)*(rk+13.D0)*(rk+14.D0)&
                *(-2.99589010560000000D12+rk*(-6.82107832800000000D10+rk*(3.42412619880000000D10+rk&
                *(1.05746915200000000D9+rk*(-1.07740631000000000D8+rk*(-3.78428000000000000D6+rk&
                *(77362.00000000000000000D0+(rk+2728.D0)*rk))))))) * den(8)
        ckplm(k,4058) = 5.52880385565380000D-9*(rk+12.D0)*(rk+13.D0)*(6.94036788480000000D12+rk&
                *(2.24217122400000000D10+rk*(-1.03546649604000000D11+rk*(-1.27370695600000000D9+rk&
                *(5.02880889000000000D8+rk*(9.12483600000000000D6+rk*(-846006.00000000000000000D0+rk&
                *(-16524.00000000000000000D0+rk*(321.00000000000000000D0+4.00000000000000000D0*rk))))))))) &
                * den(9)
        ckplm(k,4059) = 2.10094546514844000D-8*(rk+12.D0)*(-2.23176110976000000D13+rk&
                *(2.76480589680000000D11+rk*(4.09524616476000000D11+rk*(-1.66705954000000000D9+rk&
                *(-2.66830523500000000D9+rk*(-6.97473000000000000D6+rk*(7.17072300000000000D6+rk&
                *(49440.00000000000000000D0+rk*(-6765.00000000000000000D0+(rk-50.D0)*rk))))))))) * den(10)
        ckplm(k,4060) = -2.40653753280640000D-7*(-2.21855542272000000D13+11.00000000000000000D0&
                *(rk+1.D0)*rk*(4.47478416000000000D10+(rk+1.D0)*rk*(-3.67138944000000000D8+(rk+1.D0)*rk&
                *(1.35250800000000000D6+(rk*rk+rk-2120.D0)*(rk+1.D0)*rk)))) * den(11)
        ckplm(k,4061) = -2.10094546514844000D-8*(rk-11.D0)*(-2.21855542272000000D13+rk&
                *(5.36974099200000000D11+rk*(3.98692046016000000D11+rk*(-8.79510456000000000D9+rk&
                *(-2.52806818000000000D9+rk*(4.85885400000000000D7+rk*(6.63963300000000000D6+rk&
                *(-101640.00000000000000000D0+rk*(-6270.00000000000000000D0+(rk+60.D0)*rk))))))))) * den(12)
        ckplm(k,4062) = -5.52880385565380000D-9*(rk-11.D0)*(rk-12.D0)*(-6.81616615680000000D12+rk&
                *(2.23732949040000000D11+rk*(9.68118260040000000D10+rk*(-3.17765801200000000D9+rk&
                *(-4.45166925000000000D8+rk*(1.38363960000000000D7+rk*(721686.00000000000000000D0+rk&
                *(-18948.00000000000000000D0+rk*(-285.00000000000000000D0+4.00000000000000000D0*rk))))))))) &
                * den(13)
        ckplm(k,4063) = -9.75671268644787000D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)&
                *(-2.89459941120000000D12+rk*(1.33109303760000000D11+rk*(3.04613567160000000D10+rk&
                *(-1.44913706000000000D9+rk*(-8.77542110000000000D7+rk*(4.19122000000000000D6+rk&
                *(58294.00000000000000000D0+(rk-2720.D0)*rk))))))) * den(14)
        ckplm(k,4064) = 2.60179004971943000D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(7.49946859200000000D11+rk*(-4.61808039600000000D10+rk*(-5.14076496600000000D9+rk&
                *(3.45301087000000000D8+rk*(6.01513500000000000D6+rk*(-539615.00000000000000000D0+rk&
                *(2151.00000000000000000D0+88.00000000000000000D0*rk))))))) * den(15)
        ckplm(k,4065) = 2.50172120165330000D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(4.90366700928000000D12+rk*(-3.78534962448000000D11+rk*(-1.62424124600000000D10+rk&
                *(1.67344308000000000D9+rk*(-1.05723950000000000D7+rk&
                *(-956052.00000000000000000D0+9235.00000000000000000D0*rk)))))) * den(16)
        ckplm(k,4066) = 6.82287600450900000D-13*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(1.00531773756000000D12-1.00000000000000000D0*rk&
                *(9.03842219460000000D10+rk*(2.60164265000000000D8+rk*(-1.81367770000000000D8+rk&
                *(3.23675500000000000D6+9244.00000000000000000D0*rk))))) * den(17)
        ckplm(k,4067) = 2.65334066842017000D-12*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(1.23883981200000000D10+rk*(-1.19189938200000000D9+rk&
                *(2.16489010000000000D7+(580080.00000000000000000D0-12403.00000000000000000D0*rk)*rk))) &
                * den(18)
        ckplm(k,4068) = -2.80496013518703000D-11*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(-4.54704600000000000D7+rk&
                *(4.20025100000000000D6+rk*(-101565.00000000000000000D0+304.00000000000000000D0*rk))) * den(19)
        ckplm(k,4069) = 1.09393445272294000D-10*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(335520.00000000000000000D0+rk&
                *(-25292.00000000000000000D0+439.00000000000000000D0*rk)) * den(20)
        ckplm(k,4070) = 1.49504375205469000D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)&
                *(453.00000000000000000D0-20.00000000000000000D0*rk) * den(21)
        ckplm(k,4071) = 5.84426193985015000D-9*(rk-11.D0)*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4072) = 1.71890057054416000D-10*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4073) = -7.99488637462401000D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(3096.00000000000000000D0+133.00000000000000000D0*rk) * den(1)
        ckplm(k,4074) = 8.77487528922147000D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(1.87058400000000000D6+rk&
                *(143651.00000000000000000D0+2637.00000000000000000D0*rk)) * den(2)
        ckplm(k,4075) = -1.49997868191820000D-13*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(4.54557480000000000D8+rk*(4.53133580000000000D7+rk&
                *(1.33514700000000000D6+10159.00000000000000000D0*rk))) * den(3)
        ckplm(k,4076) = 7.09449376582932000D-15*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(rk+18.D0)*(2.87735563080000000D11+rk*(3.20126614380000000D10+rk&
                *(1.00603228900000000D9+(3.04798200000000000D6-165229.00000000000000000D0*rk)*rk))) * den(4)
        ckplm(k,4077) = 1.09457903815652000D-14*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(rk+17.D0)*(-4.38741572976000000D12+rk*(-4.95085377936000000D11+rk*(-1.13293136500000000D10+rk&
                *(4.01613695000000000D8+rk*(1.58309500000000000D7+92341.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4078) = 1.03376909159227000D-14*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(9.10825269393600000D13+rk*(9.74954843505600000D12+rk*(2.80859978180000000D10+rk&
                *(-2.55314910450000000D10+rk*(-7.13105935000000000D8+rk&
                *(1.72086900000000000D6+129317.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,4079) = -2.15023971051193000D-12*(rk+13.D0)*(rk+14.D0)*(rk+15.D0)&
                *(7.45176278016000000D12+rk*(7.22202625296000000D11+rk*(-2.14450266360000000D10+rk&
                *(-3.74731155200000000D9+rk*(-5.79670050000000000D7+rk*(2.97682900000000000D6+rk&
                *(70521.00000000000000000D0+67.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4080) = -1.20950983716296000D-11*(rk+13.D0)*(rk+14.D0)*(-2.01026956032000000D13+rk&
                *(-1.71665589816000000D12+rk*(1.31034510036000000D11+rk*(1.41866010840000000D10+rk&
                *(-4.55297170000000000D7+rk*(-2.71935600000000000D7+rk*(-394186.00000000000000000D0+rk&
                *(8076.00000000000000000D0+107.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,4081) = -2.68779963813991000D-12*(rk+13.D0)*(1.24470299842560000D15+rk&
                *(9.33958389196800000D13+rk*(-1.28356499395680000D13+rk*(-1.08015682924400000D12+rk&
                *(3.15784582680000000D10+rk*(3.56837163900000000D9+rk*(4.33120800000000000D6+rk&
                *(-3.22398600000000000D6+rk*(-30708.00000000000000000D0+311.00000000000000000D0*rk))))))))) &
                * den(9)
        ckplm(k,4082) = 5.10681931246583000D-12*(8.24351807784960000D15+rk*(5.61271289713920000D14+rk&
                *(-1.16516030993136000D14+rk*(-8.31306817830000000D12+rk*(5.14660671080000000D11+rk&
                *(3.98820574950000000D10+rk*(-6.99418923000000000D8+rk*(-6.79228500000000000D7+rk&
                *(-6330.00000000000000000D0+rk*(28935.00000000000000000D0+109.00000000000000000D0*rk)))))))))) &
                * den(10)
        ckplm(k,4083) = 1.28691846674139000D-9*(-3.14295351552000000D13+(rk+1.D0)*rk&
                *(5.24218615200000000D11+(rk+1.D0)*rk*(-3.08976794400000000D9+(rk+1.D0)*rk&
                *(7.55870800000000000D6+(rk*rk+rk-6670.D0)*(rk+1.D0)*rk)))) * den(11)
        ckplm(k,4084) = 5.10681931246583000D-12*(7.57451797240320000D15+rk&
                *(-7.67508636132000000D14+rk*(-8.88967491246960000D13+rk*(9.96127643702000000D12+rk&
                *(3.07132333490000000D11+rk*(-4.26561640050000000D10+rk*(-2.26543863000000000D8+rk&
                *(6.68436300000000000D7+rk*(-261840.00000000000000000D0+rk&
                *(-27845.00000000000000000D0+109.00000000000000000D0*rk)))))))))) * den(12)
        ckplm(k,4085) = 2.68779963813991000D-12*(rk-12.D0)*(1.13957968400640000D15+rk&
                *(-1.15718148029520000D14+rk*(-9.44126061781200000D12+rk*(1.17098466382400000D12+rk&
                *(1.39122189570000000D10+rk*(-3.47643951900000000D9+rk*(2.60131620000000000D7+rk&
                *(2.96712600000000000D6-1.00000000000000000D0*rk&
                *(33507.00000000000000000D0+311.00000000000000000D0*rk))))))))) * den(13)
        ckplm(k,4086) = 1.20950983716296000D-11*(rk-12.D0)*(rk-13.D0)&
                *(1.82692105344000000D13-1.00000000000000000D0*rk*(1.93611654312000000D12+rk&
                *(8.84673846920000000D10+rk*(-1.41049447400000000D10+rk*(8.42501230000000000D7+rk&
                *(2.46648400000000000D7+rk*(-447722.00000000000000000D0+rk&
                *(-7220.00000000000000000D0+107.00000000000000000D0*rk)))))))) * den(14)
        ckplm(k,4087) = 2.15023971051193000D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)&
                *(6.71180156640000000D12-1.00000000000000000D0*rk*(7.54097073420000000D11+rk&
                *(1.05796058920000000D10+rk*(-3.48708331700000000D9+rk*(7.17956800000000000D7+rk&
                *(2.55511000000000000D6+rk*(-70052.00000000000000000D0+67.00000000000000000D0*rk))))))) &
                * den(15)
        ckplm(k,4088) = 1.03376909159227000D-14*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(8.13858812956800000D13+rk*(-9.61964221846800000D12+rk*(1.00386566408000000D11+rk&
                *(2.26644449550000000D10+rk*(-7.19770525000000000D8+rk&
                *(-944967.00000000000000000D0+129317.00000000000000000D0*rk)))))) * den(16)
        ckplm(k,4089) = 1.09457903815652000D-14*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(3.90404554056000000D12+rk*(-4.71284771646000000D11+rk*(1.24400924450000000D10+rk&
                *(3.39213305000000000D8+rk*(-1.53692450000000000D7+92341.00000000000000000D0*rk))))) * den(17)
        ckplm(k,4090) = 7.09449376582932000D-15*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(2.56725720720000000D11-1.00000000000000000D0*rk&
                *(3.00104017220000000D10+rk*(-9.95896969000000000D8+rk&
                *(3.70889800000000000D6+165229.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4091) = -1.49997868191820000D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(-4.10569110000000000D8+rk*(4.26735410000000000D7+rk&
                *(-1.30467000000000000D6+10159.00000000000000000D0*rk))) * den(19)
        ckplm(k,4092) = 8.77487528922147000D-13*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(1.72957000000000000D6+rk&
                *(-138377.00000000000000000D0+2637.00000000000000000D0*rk)) * den(20)
        ckplm(k,4093) = 7.99488637462401000D-12*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)&
                *(2963.00000000000000000D0-133.00000000000000000D0*rk) * den(21)
        ckplm(k,4094) = 1.71890057054416000D-10*(rk-12.D0)*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4095) = -4.91114448726903000D-12*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4096) = 1.14212662494629000D-13*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(7267.00000000000000000D0+316.00000000000000000D0&
                *rk) * den(1)
        ckplm(k,4097) = -5.84991685948098000D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(1.10171100000000000D6+rk&
                *(88354.00000000000000000D0+1723.00000000000000000D0*rk)) * den(2)
        ckplm(k,4098) = 1.49997868191820000D-14*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(rk+19.D0)*(2.07484680000000000D8+rk*(2.25235930000000000D7+rk&
                *(761337.00000000000000000D0+7664.00000000000000000D0*rk))) * den(3)
        ckplm(k,4099) = -2.02699821880838000D-16*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(rk+18.D0)*(5.27873007480000000D11+rk*(6.75056172680000000D10+rk*(2.80883612900000000D9+rk&
                *(3.77982520000000000D7+1231.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4100) = -3.64859679385508000D-16*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(-7.79480114868000000D12+rk*(-1.08149462845800000D12+rk*(-4.56325160750000000D10+rk&
                *(-2.78285690000000000D8+rk*(1.84940750000000000D7+239348.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4101) = -2.06753818318455000D-15*(rk+14.D0)*(rk+15.D0)*(rk+16.D0)&
                *(2.99099665378800000D13+rk*(4.27285126447800000D12+rk*(1.52094057209000000D11+rk&
                *(-3.31584756000000000D9+rk*(-2.71145230000000000D8+rk&
                *(-3.52747800000000000D6+5021.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,4102) = 3.07177101501704000D-14*(rk+14.D0)*(rk+15.D0)*(3.72663080899200000D13+rk&
                *(5.31456385987200000D12+rk*(1.20001191898000000D11+rk*(-1.34147070890000000D10+rk&
                *(-6.84968585000000000D8+rk*(-3.42244700000000000D6+rk&
                *(240247.00000000000000000D0+2344.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4103) = 1.15191413063139000D-13*(rk+14.D0)*(-1.60952401209600000D14+rk&
                *(-2.26275900400800000D13+rk*(-1.03003762572000000D11+rk*(1.09171010032000000D11+rk&
                *(4.16651260900000000D9+rk*(-6.76398800000000000D7+rk*(-4.97247800000000000D6+rk&
                *(-38552.00000000000000000D0+361.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,4104) = 2.68779963813991000D-13*(9.94683891801600000D14-1.00000000000000000D0*rk&
                *(-1.38520196704080000D14+rk*(2.26248226546800000D12+rk*(1.03525803024400000D12+rk&
                *(2.74932024570000000D10+rk*(-1.78883636400000000D9+rk*(-7.81115580000000000D7+rk&
                *(207636.00000000000000000D0+rk*(32433.00000000000000000D0+164.00000000000000000D0*rk))))))))) &
                * den(9)
        ckplm(k,4105) = -6.63886510620558000D-11*(4.01437257984000000D12+rk&
                *(2.57038876656000000D11+rk*(-4.06683770280000000D10+rk*(-2.69696454400000000D9+rk&
                *(1.17829257000000000D8+rk*(8.42916900000000000D6+rk*(-82614.00000000000000000D0+rk&
                *(-7842.00000000000000000D0+(rk-15.D0)*rk)))))))) * den(10)
        ckplm(k,4106) = 8.36497003381903000D-9*(3.09955968000000000D10+(rk+1.D0)*rk&
                *(-3.82673952000000000D8+(rk+1.D0)*rk*(1.57281200000000000D6+(rk*rk+rk-2400.D0)*(rk+1.D0)&
                *rk))) * den(11)
        ckplm(k,4107) = 6.63886510620558000D-11*(-3.71947161600000000D12+rk&
                *(3.29856006816000000D11+rk*(3.19558745280000000D10+rk*(-3.08261114800000000D9+rk&
                *(-7.47174960000000000D7+rk*(8.76113700000000000D6+rk*(28224.00000000000000000D0+rk&
                *(-7686.00000000000000000D0+(rk+24.D0)*rk)))))))) * den(12)
        ckplm(k,4108) = 2.68779963813991000D-13*(8.54907267110400000D14+rk&
                *(-1.40057834271120000D14+rk*(6.61619378388000000D11+rk*(9.08964552724000000D11+rk&
                *(-3.52606932930000000D10+rk*(-1.31760224400000000D9+rk*(7.86706620000000000D7+rk&
                *(-45924.00000000000000000D0+rk*(-30957.00000000000000000D0+164.00000000000000000D0&
                *rk))))))))) * den(13)
        ckplm(k,4109) = 1.15191413063139000D-13*(rk-13.D0)&
                *(1.38532756723200000D14-1.00000000000000000D0*rk*(2.21110441725600000D13+rk&
                *(-4.04915085684000000D11+rk*(-9.19266408200000000D10+rk*(4.43149942900000000D9+rk&
                *(3.86348200000000000D7+rk*(-4.69250600000000000D6+rk&
                *(41440.00000000000000000D0+361.00000000000000000D0*rk)))))))) * den(14)
        ckplm(k,4110) = 3.07177101501704000D-14*(rk-13.D0)*(rk-14.D0)*(3.20844788208000000D13+rk&
                *(-5.03703869184000000D12+rk*(1.56173280606000000D11+rk*(1.07137801190000000D10+rk&
                *(-6.64334685000000000D8+rk&
                *(4.81470500000000000D6+(223839.00000000000000000D0-2344.00000000000000000D0*rk)*rk)))))) &
                * den(15)
        ckplm(k,4111) = 2.06753818318455000D-15*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)&
                *(2.57922575654400000D13+rk*(-3.95978252078400000D12+rk*(1.60450078604000000D11+rk&
                *(2.26664184000000000D9+rk*(-2.53432525000000000D8+rk&
                *(3.55760400000000000D6+5021.00000000000000000D0*rk)))))) * den(16)
        ckplm(k,4112) = 3.64859679385508000D-16*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(6.75864249588000000D12+rk*(-9.91137232938000000D11+rk*(4.46890880350000000D10+rk&
                *(-3.49868510000000000D8+rk*(-1.72973350000000000D7+239348.00000000000000000D0*rk))))) * den(17)
        ckplm(k,4113) = 2.02699821880838000D-16*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(4.63138429320000000D11+rk*(-6.20013348420000000D10+rk*(2.69544875900000000D9+rk&
                *(-3.77933280000000000D7+1231.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4114) = -1.49997868191820000D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(-1.85714760000000000D8+rk*(2.10239110000000000D7+rk&
                *(-738345.00000000000000000D0+7664.00000000000000000D0*rk))) * den(19)
        ckplm(k,4115) = 5.84991685948098000D-14*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(1.01508000000000000D6+rk&
                *(-84908.00000000000000000D0+1723.00000000000000000D0*rk)) * den(20)
        ckplm(k,4116) = -1.14212662494629000D-13*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(-6951.00000000000000000D0+316.00000000000000000D0&
                *rk) * den(21)
        ckplm(k,4117) = 4.91114448726903000D-12*(rk-13.D0)*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4118) = 1.36420680201918000D-13*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4119) = -6.34514791636826000D-15*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(4214.00000000000000000D0+185.00000000000000000D0*rk) * den(1)
        ckplm(k,4120) = 1.62497690541138000D-15*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(rk+20.D0)*(1.47919800000000000D6+rk&
                *(122635.00000000000000000D0+2497.00000000000000000D0*rk)) * den(2)
        ckplm(k,4121) = -8.33321489954555000D-16*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(rk+19.D0)*(1.60316520000000000D8+rk*(1.85472020000000000D7+rk&
                *(686403.00000000000000000D0+7981.00000000000000000D0*rk))) * den(3)
        ckplm(k,4122) = 3.94138542546074000D-17*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(1.32782679960000000D11+rk*(1.87864888460000000D10+rk*(9.20156579000000000D8+rk&
                *(1.75004740000000000D7+92341.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4123) = 2.02699821880838000D-17*(rk+15.D0)*(rk+16.D0)*(rk+17.D0)&
                *(-7.71209294184000000D12+rk*(-1.23686388980400000D12+rk*(-6.89156529400000000D10+rk&
                *(-1.41833702500000000D9+rk*(-1.73846000000000000D6+165229.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4124) = -5.74316161995707000D-17*(rk+15.D0)*(rk+16.D0)*(-6.54285667543200000D13+rk&
                *(-1.13471489431320000D13+rk*(-6.50103440600000000D11+rk*(-9.38131396500000000D9+rk&
                *(3.33222115000000000D8+rk*(1.05500970000000000D7+60685.00000000000000000D0*rk)))))) * den(6)
        ckplm(k,4125) = -1.70653945278724000D-15*(rk+15.D0)*(4.42678192588800000D13+rk&
                *(8.07515269828800000D12+rk*(4.39654063100000000D11+rk*(-8.77193128000000000D8+rk&
                *(-7.77662095000000000D8+rk*(-2.01527830000000000D7+rk&
                *(-52885.00000000000000000D0+1823.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4126) = 3.19976147397608000D-15*(4.07864803852800000D14+rk*(7.73029837994400000D13+rk&
                *(3.74838470859600000D12+rk*(-1.22901402316000000D11+rk*(-1.49443259770000000D10+rk&
                *(-2.91267760000000000D8+rk*(5.01565400000000000D6+rk&
                *(176876.00000000000000000D0+767.00000000000000000D0*rk)))))))) * den(8)
        ckplm(k,4127) = 4.97740673729613000D-15*(-2.82989565734400000D14+rk&
                *(-3.55585571575200000D13+rk*(2.43952668972000000D11+rk*(1.75831126940000000D11+rk&
                *(4.58058719300000000D9+rk*(-1.66571440000000000D8+rk*(-6.42170200000000000D6+rk&
                *(-5980.00000000000000000D0+737.00000000000000000D0*rk)))))))) * den(9)
        ckplm(k,4128) = -4.72853640043132000D-14*(-3.05820269798400000D13+rk&
                *(-1.80511697937600000D12+rk*(2.16791616708000000D11+rk*(1.30277359240000000D10+rk&
                *(-3.91196309000000000D8+rk*(-2.49032000000000000D7+rk*(116158.00000000000000000D0+rk&
                *(10588.00000000000000000D0+19.00000000000000000D0*rk)))))))) * den(10)
        ckplm(k,4129) = -3.97197057636231000D-12*(3.58663334400000000D11+(rk+1.D0)*rk&
                *(-3.20568811200000000D9+(rk+1.D0)*rk*(8.79994800000000000D6+(rk*rk+rk-7580.D0)*(rk+1.D0)&
                *rk))) * den(11)
        ckplm(k,4130) = -4.72853640043132000D-14*(-2.85735123072000000D13+rk&
                *(2.19817735876800000D12+rk*(1.75611783636000000D11+rk*(-1.43415355160000000D10+rk&
                *(-2.65307189000000000D8+rk*(2.53788640000000000D7+rk*(42574.00000000000000000D0+rk&
                *(-10436.00000000000000000D0+19.00000000000000000D0*rk)))))))) * den(12)
        ckplm(k,4131) = 4.97740673729613000D-15*(-2.47358146291200000D14+rk&
                *(3.55380858381600000D13+rk*(-2.54487653604000000D11+rk*(-1.55971247236000000D11+rk&
                *(5.31737975300000000D9+rk*(1.28208080000000000D8+rk*(-6.35920600000000000D6+rk&
                *(11876.00000000000000000D0+737.00000000000000000D0*rk)))))))) * den(13)
        ckplm(k,4132) = 3.19976147397608000D-15*(3.34418457945600000D14+rk&
                *(-6.94958022784800000D13+rk*(4.03040717917200000D12+rk*(6.61309413800000000D10+rk&
                *(-1.34188893370000000D10+rk*(3.17690240000000000D8+rk*(3.79899800000000000D6+rk&
                *(-170740.00000000000000000D0+767.00000000000000000D0*rk)))))))) * den(14)
        ckplm(k,4133) = 1.70653945278724000D-15*(rk-14.D0)&
                *(3.66324402528000000D13-1.00000000000000000D0*rk*(7.19622320724000000D12+rk&
                *(-4.37820366186000000D11+rk*(2.03304892700000000D9+rk*(6.77755260000000000D8+rk&
                *(-1.97971900000000000D7+rk*(65646.00000000000000000D0+1823.00000000000000000D0*rk))))))) &
                * den(15)
        ckplm(k,4134) = 5.74316161995707000D-17*(rk-14.D0)*(rk-15.D0)&
                *(5.47218172051200000D13-1.00000000000000000D0*rk*(1.00763665059120000D13+rk&
                *(-6.20064756710000000D11+rk*(1.06099151550000000D10+rk*(2.81381905000000000D8+rk&
                *(-1.01859870000000000D7+60685.00000000000000000D0*rk)))))) * den(16)
        ckplm(k,4135) = 2.02699821880838000D-17*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)&
                *(6.54272827164000000D12+rk*(-1.10327981501400000D12+rk*(6.46727249150000000D10+rk&
                *(-1.40973089500000000D9+rk*(2.56460500000000000D6+165229.00000000000000000D0*rk))))) * den(17)
        ckplm(k,4136) = 3.94138542546074000D-17*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(1.14898939560000000D11+rk*(-1.69983077460000000D10+rk*(8.68209203000000000D8+rk&
                *(-1.71311100000000000D7+92341.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4137) = -8.33321489954555000D-16*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(-1.42447740000000000D8+rk*(1.71983390000000000D7+rk&
                *(-662460.00000000000000000D0+7981.00000000000000000D0*rk))) * den(19)
        ckplm(k,4138) = 1.62497690541138000D-15*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(rk-19.D0)*(1.35906000000000000D6+rk&
                *(-117641.00000000000000000D0+2497.00000000000000000D0*rk)) * den(20)
        ckplm(k,4139) = 6.34514791636826000D-15*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(4029.00000000000000000D0-185.00000000000000000D0*rk) &
                * den(21)
        ckplm(k,4140) = 1.36420680201918000D-13*(rk-14.D0)*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4141) = -3.68704541086264000D-15*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4142) = 8.57452421130846000D-17*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(rk+20.D0)*(rk+21.D0)*(9675.00000000000000000D0+428.00000000000000000D0*rk) * den(1)
        ckplm(k,4143) = -1.31754884222545000D-16*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(rk+20.D0)*(648825.00000000000000000D0+rk*(55214.00000000000000000D0+1161.00000000000000000D0&
                *rk)) * den(2)
        ckplm(k,4144) = 5.63055060780105000D-17*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(9.62718900000000000D7+rk*(1.16982890000000000D7+rk&
                *(461733.00000000000000000D0+5872.00000000000000000D0*rk))) * den(3)
        ckplm(k,4145) = -2.81527530390053000D-17*(rk+16.D0)*(rk+17.D0)*(rk+18.D0)&
                *(8.53368894000000000D9+rk*(1.30338770400000000D9+rk*(7.11569690000000000D7+rk&
                *(1.61050800000000000D6+12199.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4146) = 3.04049732821257000D-17*(rk+16.D0)*(rk+17.D0)*(2.63896715880000000D11+rk&
                *(4.71323511780000000D10+rk*(3.10712202500000000D9+rk*(8.96926700000000000D7+rk&
                *(974755.00000000000000000D0+772.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4147) = 5.74316161995707000D-17*(rk+16.D0)*(-3.70085362020000000D12+rk&
                *(-7.39031313930000000D11+rk*(-5.42685217610000000D10+rk*(-1.66581510000000000D9+rk&
                *(-1.16057900000000000D7+rk*(357510.00000000000000000D0+4711.00000000000000000D0*rk)))))) &
                * den(6)
        ckplm(k,4148) = 8.53269726393622000D-16*(5.42970948240000000D12+rk&
                *(1.17979341066000000D12-1.00000000000000000D0*rk*(-9.15902037120000000D10+rk&
                *(-2.54032736300000000D9+rk*(2.73817950000000000D7+rk*(2.72597500000000000D6+rk&
                *(40677.00000000000000000D0+88.00000000000000000D0*rk))))))) * den(7)
        ckplm(k,4149) = 4.79964221096412000D-14*(-1.18496649600000000D11+rk&
                *(-1.97298950400000000D10+rk*(-8.82111636000000000D8+rk*(1.46141800000000000D7+rk&
                *(1.90858500000000000D6+rk*(33145.00000000000000000D0-1.00000000000000000D0*rk&
                *(189.00000000000000000D0+5.00000000000000000D0*rk))))))) * den(8)
        ckplm(k,4150) = 1.24435168432403000D-14*(5.16820780800000000D11+rk*(5.75362023600000000D10+rk&
                *(-2.92785260000000000D7+rk*(-1.79699579000000000D8+rk*(-4.26771500000000000D6+rk&
                *(82135.00000000000000000D0+rk*(2521.00000000000000000D0+4.00000000000000000D0*rk))))))) &
                * den(9)
        ckplm(k,4151) = 4.72853640043132000D-14*(-1.43806917888000000D11+rk&
                *(-7.65884743200000000D9+rk*(6.90130266000000000D8+5.00000000000000000D0*rk&
                *(7.28322700000000000D6+rk*(-143934.00000000000000000D0+rk&
                *(-7700.00000000000000000D0+(rk+12.D0)*rk)))))) * den(10)
        ckplm(k,4152) = -5.95795586454346000D-12*(-1.13861376000000000D9+(rk+1.D0)*rk&
                *(7.14087600000000000D6+5.00000000000000000D0*(rk*rk+rk-2456.D0)*(rk+1.D0)*rk)) * den(11)
        ckplm(k,4153) = -4.72853640043132000D-14*(1.35495037440000000D11+rk&
                *(-8.92717370400000000D9+rk*(-5.76949636000000000D8+5.00000000000000000D0*rk&
                *(7.78175800000000000D6+rk*(105289.00000000000000000D0+rk*(-7751.00000000000000000D0+(rk-5.D0)&
                *rk)))))) * den(12)
        ckplm(k,4154) = 1.24435168432403000D-14*(4.59430652160000000D11+rk&
                *(-5.70731271120000000D10+rk*(4.83430302000000000D8+rk&
                *(1.61857649000000000D8-1.00000000000000000D0*rk*(4.64071500000000000D6+rk&
                *(67093.00000000000000000D0+rk*(-2493.00000000000000000D0+4.00000000000000000D0*rk))))))) &
                * den(13)
        ckplm(k,4155) = 4.79964221096412000D-14*(-9.96616051200000000D10+rk&
                *(1.79292967440000000D10+rk*(-9.14836846000000000D8+rk*(-7.31489500000000000D6+rk&
                *(1.74020000000000000D6+rk*(-34174.00000000000000000D0+rk&
                *(-154.00000000000000000D0+5.00000000000000000D0*rk))))))) * den(14)
        ckplm(k,4156) = 8.53269726393622000D-16*(4.33894125168000000D12+rk&
                *(-1.00433012607600000D12+rk*(8.38315822960000000D10+rk*(-2.62340525300000000D9+rk&
                *(-1.43589950000000000D7+rk*(2.48376100000000000D6+rk&
                *(-40061.00000000000000000D0+88.00000000000000000D0*rk))))))) * den(15)
        ckplm(k,4157) = 5.74316161995707000D-17*(rk-15.D0)*(3.01443697152000000D12+rk&
                *(-6.35443533264000000D11+rk*(4.93442156360000000D10+rk*(-1.61591106000000000D9+rk&
                *(1.33226750000000000D7+(329244.00000000000000000D0-4711.00000000000000000D0*rk)*rk))))) &
                * den(16)
        ckplm(k,4158) = 3.04049732821257000D-17*(rk-15.D0)*(rk-16.D0)*(2.19782768040000000D11+rk&
                *(-4.11832899780000000D10+rk*(2.84388482500000000D9+rk&
                *(-8.58013700000000000D7+(970895.00000000000000000D0-772.00000000000000000D0*rk)*rk)))) &
                * den(17)
        ckplm(k,4159) = 2.81527530390053000D-17*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)&
                *(7.29985989600000000D9+rk*(-1.16585649400000000D9+rk*(6.63986390000000000D7+rk&
                *(-1.56171200000000000D6+12199.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4160) = 5.63055060780105000D-17*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(8.50294620000000000D7+rk&
                *(-1.07924390000000000D7+(444117.00000000000000000D0-5872.00000000000000000D0*rk)*rk)) * den(19)
        ckplm(k,4161) = 1.31754884222545000D-16*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(rk-19.D0)*(594772.00000000000000000D0+rk*(-52892.00000000000000000D0+1161.00000000000000000D0&
                *rk)) * den(20)
        ckplm(k,4162) = 8.57452421130846000D-17*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(rk-19.D0)*(rk-20.D0)*(9247.00000000000000000D0-428.00000000000000000D0*rk) * den(21)
        ckplm(k,4163) = 3.68704541086264000D-15*(rk-15.D0)*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4164) = 9.70275108121746000D-17*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4165) = -4.51290747963603000D-18*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(rk+21.D0)*(5504.00000000000000000D0+245.00000000000000000D0*rk) * den(1)
        ckplm(k,4166) = 1.65106371206196000D-19*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(1.76175360000000000D7+rk*(1.53090500000000000D6+32999.00000000000000000D0*rk)) * den(2)
        ckplm(k,4167) = -8.46699339518955000D-20*(rk+17.D0)*(rk+18.D0)*(rk+19.D0)&
                *(2.46592320000000000D9+rk*(3.11500286000000000D8+rk&
                *(1.29000990000000000D7+174463.00000000000000000D0*rk))) * den(3)
        ckplm(k,4168) = 2.81527530390053000D-18*(rk+17.D0)*(rk+18.D0)*(3.69485760000000000D9+rk&
                *(5.99073266000000000D8+rk*(3.53551190000000000D7+rk&
                *(891634.00000000000000000D0+7981.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4169) = -1.44785587057741000D-18*(rk+17.D0)*(2.67006620160000000D11+rk&
                *(5.18176934160000000D10+rk*(3.83263829000000000D9+rk*(1.32282605000000000D8+rk&
                *(2.04433000000000000D6+10159.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4170) = 4.10225829996934000D-18*(2.74239440640000000D12+rk*(6.09882565848000000D11+rk&
                *(5.25450626300000000D10+rk*(2.15387680500000000D9+rk&
                *(3.98723750000000000D7+(186987.00000000000000000D0-1765.00000000000000000D0*rk)*rk))))) &
                * den(6)
        ckplm(k,4171) = 8.53269726393622000D-16*(-1.94417172000000000D10+rk&
                *(-3.60604119600000000D9+rk*(-2.36975144000000000D8+rk*(-5.78236500000000000D6+rk&
                *(11245.00000000000000000D0+rk*(2241.00000000000000000D0+19.00000000000000000D0*rk)))))) &
                * den(7)
        ckplm(k,4172) = 1.59988073698804000D-15*(1.35573026400000000D10+rk*(1.93984918800000000D9+rk&
                *(7.65466240000000000D7+rk*(-519105.00000000000000000D0+rk&
                *(-81545.00000000000000000D0+(rk-1083.D0)*rk))))) * den(8)
        ckplm(k,4173) = 3.55529052664009000D-16*(-7.22595254400000000D10-1.00000000000000000D0*rk&
                *(6.96783568800000000D9+rk*(2.36859820000000000D7+rk*(-1.28464050000000000D7+rk&
                *(-255425.00000000000000000D0+rk*(2277.00000000000000000D0+43.00000000000000000D0*rk)))))) &
                * den(9)
        ckplm(k,4174) = 1.77764526332005000D-16*(1.57384886976000000D11-1.00000000000000000D0*rk&
                *(-7.37590236000000000D9+rk*(4.88219710000000000D8+rk*(2.16680790000000000D7+rk&
                *(-262319.00000000000000000D0+(rk-10599.D0)*rk))))) * den(10)
        ckplm(k,4175) = 1.49322202118884000D-14*(-1.89294537600000000D9+(rk+1.D0)*rk&
                *(7.96906800000000000D6+(rk*rk+rk-7760.D0)*(rk+1.D0)*rk)) * den(11)
        ckplm(k,4176) = 1.77764526332005000D-16*(1.49542684704000000D11-1.00000000000000000D0*rk&
                *(8.28634126800000000D9+rk*(4.21747564000000000D8+rk*(-2.26113450000000000D7+rk&
                *(-209309.00000000000000000D0+(rk+10605.D0)*rk))))) * den(12)
        ckplm(k,4177) = 3.55529052664009000D-16*(-6.53279644800000000D10+rk*(6.88295733600000000D9+rk&
                *(-6.06705220000000000D7+rk*(-1.18027950000000000D7+rk&
                *(266165.00000000000000000D0+(2019.00000000000000000D0-43.00000000000000000D0*rk)*rk))))) &
                * den(13)
        ckplm(k,4178) = 1.59988073698804000D-15*(1.16944387200000000D10+rk*(-1.78551938400000000D9+rk&
                *(7.76255140000000000D7+rk*(203775.00000000000000000D0+rk&
                *(-76115.00000000000000000D0+(rk+1089.D0)*rk))))) * den(14)
        ckplm(k,4179) = 8.53269726393622000D-16*(-1.60668597600000000D10+rk*(3.14947189200000000D9+rk&
                *(-2.19582704000000000D8+rk*(5.80531500000000000D6+rk*(325.00000000000000000D0+rk&
                *(-2127.00000000000000000D0+19.00000000000000000D0*rk)))))) * den(15)
        ckplm(k,4180) = 4.10225829996934000D-18*(2.18294271000000000D12-1.00000000000000000D0*rk&
                *(5.11095527028000000D11+rk*(-4.63207701200000000D10+rk*(1.99629247500000000D9+rk&
                *(-3.89109650000000000D7+rk*(197577.00000000000000000D0+1765.00000000000000000D0*rk)))))) &
                * den(16)
        ckplm(k,4181) = 1.44785587057741000D-18*(rk-16.D0)*(2.18891316600000000D11+rk&
                *(-4.45411381260000000D10+rk*(3.44795486500000000D9+rk&
                *(-1.24206875000000000D8+(1.99353500000000000D6-10159.00000000000000000D0*rk)*rk)))) * den(17)
        ckplm(k,4182) = 2.81527530390053000D-18*(rk-16.D0)*(rk-17.D0)*(3.13025580000000000D9+rk&
                *(-5.31006006000000000D8+rk*(3.27281030000000000D7+rk&
                *(-859710.00000000000000000D0+7981.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4183) = 8.46699339518955000D-20*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)&
                *(2.16714855000000000D9+rk&
                *(-2.86223477000000000D8+(1.23767100000000000D7-174463.00000000000000000D0*rk)*rk)) * den(19)
        ckplm(k,4184) = 1.65106371206196000D-19*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(1.61196300000000000D7+rk*(-1.46490700000000000D6+32999.00000000000000000D0*rk)) * den(20)
        ckplm(k,4185) = 4.51290747963603000D-18*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(rk-20.D0)*(5259.00000000000000000D0-245.00000000000000000D0*rk) * den(21)
        ckplm(k,4186) = 9.70275108121746000D-17*(rk-16.D0)*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4187) = -2.48788489261986000D-18*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(rk+22.D0) * den(0)
        ckplm(k,4188) = 5.78577882004619000D-20*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(12427.00000000000000000D0+556.00000000000000000D0*rk) * den(1)
        ckplm(k,4189) = -1.27004900927843000D-20*(rk+18.D0)*(rk+19.D0)*(rk+20.D0)&
                *(7.47671900000000000D6+rk*(660866.00000000000000000D0+14527.00000000000000000D0*rk)) * den(2)
        ckplm(k,4190) = 1.41116556586493000D-20*(rk+18.D0)*(rk+19.D0)*(5.42958750000000000D8+rk&
                *(7.07682230000000000D7+rk*(3.04154700000000000D6+43024.00000000000000000D0*rk))) * den(3)
        ckplm(k,4191) = -9.38425101300175000D-19*(rk+18.D0)*(4.54568100000000000D8+rk&
                *(7.73072080000000000D7+rk*(4.83829900000000000D6+rk&
                *(131492.00000000000000000D0+1301.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4192) = 2.41309311762902000D-19*(7.27792052400000000D10+rk*(1.50866423340000000D10+rk&
                *(1.21416567500000000D9+rk*(4.69962500000000000D7+rk&
                *(862105.00000000000000000D0+5836.00000000000000000D0*rk))))) * den(5)
        ckplm(k,4193) = 4.10225829996934000D-18*(-8.02654380000000000D9-1.00000000000000000D0*rk&
                *(1.47119957400000000D9+rk*(1.01058645000000000D8+rk*(3.14321500000000000D6+rk&
                *(40875.00000000000000000D0+131.00000000000000000D0*rk))))) * den(6)
        ckplm(k,4194) = 3.28180663997547000D-17*(1.59855210000000000D9+rk*(2.45850364000000000D8+rk&
                *(1.30756150000000000D7+rk&
                *(256550.00000000000000000D0+(395.00000000000000000D0-24.00000000000000000D0*rk)*rk)))) * den(7)
        ckplm(k,4195) = 4.10225829996934000D-17*(-1.77950304000000000D9+rk*(-2.12483748000000000D8+rk&
                *(-7.02380000000000000D6+rk*(9205.00000000000000000D0+rk&
                *(3080.00000000000000000D0+23.00000000000000000D0*rk))))) * den(8)
        ckplm(k,4196) = 7.74871012216430000D-17*(1.16853948000000000D9+rk*(9.47320260000000000D7+rk&
                *(507755.00000000000000000D0+rk*(-92830.00000000000000000D0+rk&
                *(-1355.00000000000000000D0+4.00000000000000000D0*rk))))) * den(9)
        ckplm(k,4197) = -1.00733231588136000D-15*(1.00822968000000000D8+rk*(4.02829800000000000D6+rk&
                *(-189055.00000000000000000D0+rk*(-6607.00000000000000000D0+(rk+43.D0)*rk)))) * den(10)
        ckplm(k,4198) = 1.26923871801051000D-13*(818748.00000000000000000D0+(rk*rk+rk-2168.D0)&
                *(rk+1.D0)*rk) * den(11)
        ckplm(k,4199) = 1.00733231588136000D-15*(-9.66122640000000000D7+rk*(4.38642000000000000D6+rk&
                *(168986.00000000000000000D0+rk*(-6769.00000000000000000D0+(rk-38.D0)*rk)))) * den(12)
        ckplm(k,4200) = -7.74871012216430000D-17*(-1.07440668000000000D9+rk*(9.34434660000000000D7+rk&
                *(-778075.00000000000000000D0+rk*(-87370.00000000000000000D0+rk&
                *(1375.00000000000000000D0+4.00000000000000000D0*rk))))) * den(13)
        ckplm(k,4201) = 4.10225829996934000D-17*(-1.57404924000000000D9+rk*(1.98420738000000000D8+rk&
                *(-7.03316500000000000D6+rk&
                *(2885.00000000000000000D0+(2965.00000000000000000D0-23.00000000000000000D0*rk)*rk)))) * den(14)
        ckplm(k,4202) = 3.28180663997547000D-17*(1.36552122000000000D9+rk*(-2.20467084000000000D8+rk&
                *(1.23085750000000000D7+rk*(-254730.00000000000000000D0+rk&
                *(515.00000000000000000D0+24.00000000000000000D0*rk))))) * den(15)
        ckplm(k,4203) = 4.10225829996934000D-18*(-6.65330040000000000D9+rk*(1.27834908400000000D9+rk&
                *(-9.18729400000000000D7+rk*(2.98102500000000000D6+rk&
                *(-40220.00000000000000000D0+131.00000000000000000D0*rk))))) * den(16)
        ckplm(k,4204) = 2.41309311762902000D-19*(5.88605886000000000D10+rk&
                *(-1.27958804940000000D10+rk*(1.07829119500000000D9+rk&
                *(-4.36061900000000000D7+(832925.00000000000000000D0-5836.00000000000000000D0*rk)*rk)))) &
                * den(17)
        ckplm(k,4205) = 9.38425101300175000D-19*(rk-17.D0)*(3.81969000000000000D8+rk&
                *(-6.80198820000000000D7+rk*(4.45162900000000000D6+rk&
                *(-126288.00000000000000000D0+1301.00000000000000000D0*rk)))) * den(18)
        ckplm(k,4206) = 1.41116556586493000D-20*(rk-17.D0)*(rk-18.D0)*(4.75189050000000000D8+rk&
                *(-6.48142010000000000D7+(2.91247500000000000D6-43024.00000000000000000D0*rk)*rk)) * den(19)
        ckplm(k,4207) = 1.27004900927843000D-20*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)&
                *(6.83038000000000000D6+rk*(-631812.00000000000000000D0+14527.00000000000000000D0*rk)) * den(20)
        ckplm(k,4208) = 5.78577882004619000D-20*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)&
                *(11871.00000000000000000D0-556.00000000000000000D0*rk) * den(21)
        ckplm(k,4209) = 2.48788489261986000D-18*(rk-17.D0)*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)&
                *(rk-21.D0) * den(22)
        ckplm(k,4210) = 6.21971223154966000D-20*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4211) = -2.89288941002310000D-21*(rk+19.D0)*(rk+20.D0)*(rk+21.D0)&
                *(6966.00000000000000000D0+313.00000000000000000D0*rk) * den(1)
        ckplm(k,4212) = 3.17512252319608000D-22*(rk+19.D0)*(rk+20.D0)*(9.38957400000000000D6+rk&
                *(841711.00000000000000000D0+18797.00000000000000000D0*rk)) * den(2)
        ckplm(k,4213) = -7.05582782932463000D-22*(rk+19.D0)*(3.80860920000000000D8+rk&
                *(5.09307580000000000D7+rk*(2.25437700000000000D6+32999.00000000000000000D0*rk))) * den(3)
        ckplm(k,4214) = 2.34606275325044000D-20*(7.08818040000000000D8+rk*(1.25327766000000000D8+rk&
                *(8.21138300000000000D6+rk*(235794.00000000000000000D0+2497.00000000000000000D0*rk)))) * den(4)
        ckplm(k,4215) = 1.08589190293306000D-19*(-3.86384040000000000D8-1.00000000000000000D0*rk&
                *(6.30270020000000000D7+rk*(3.75356100000000000D6+rk&
                *(95998.00000000000000000D0+879.00000000000000000D0*rk)))) * den(5)
        ckplm(k,4216) = 1.02556457499233000D-19*(8.37623160000000000D8+rk*(1.21425806000000000D8+rk&
                *(6.22460300000000000D6+rk*(129994.00000000000000000D0+877.00000000000000000D0*rk)))) * den(6)
        ckplm(k,4217) = 1.64090331998773000D-18*(-8.99650800000000000D7-1.00000000000000000D0*rk&
                *(1.10056980000000000D7+rk*(442669.00000000000000000D0+rk&
                *(6102.00000000000000000D0+11.00000000000000000D0*rk)))) * den(7)
        ckplm(k,4218) = 9.23008117493101000D-18*(2.36653200000000000D7+rk*(2.26157400000000000D6+rk&
                *(58207.00000000000000000D0+(66.00000000000000000D0-7.00000000000000000D0*rk)*rk))) * den(8)
        ckplm(k,4219) = 2.27903238887185000D-19*(-1.24465596000000000D9+rk*(-8.13253220000000000D7+rk&
                *(-467201.00000000000000000D0+rk*(35282.00000000000000000D0+281.00000000000000000D0*rk)))) &
                * den(9)
        ckplm(k,4220) = 2.96274210553341000D-19*(1.10608956000000000D9+rk*(3.60663780000000000D7+rk&
                *(-1.12987100000000000D6+rk*(-27618.00000000000000000D0+71.00000000000000000D0*rk)))) * den(10)
        ckplm(k,4221) = 7.46611010594419000D-17*(-4.54860000000000000D6-1.00000000000000000D0*(rk&
                *rk+rk-6842.D0)*(rk+1.D0)*rk) * den(11)
        ckplm(k,4222) = 2.96274210553341000D-19*(1.06892100000000000D9+rk*(-3.82429820000000000D7+rk&
                *(-1.04659100000000000D6+rk*(27902.00000000000000000D0+71.00000000000000000D0*rk)))) * den(12)
        ckplm(k,4223) = 2.27903238887185000D-19*(-1.16383284000000000D9+rk*(8.02861980000000000D7+rk&
                *(-571361.00000000000000000D0+rk*(-34158.00000000000000000D0+281.00000000000000000D0*rk)))) &
                * den(13)
        ckplm(k,4224) = 9.23008117493101000D-18*(2.14618800000000000D7-1.00000000000000000D0*rk&
                *(2.14538600000000000D6+rk*(-57967.00000000000000000D0+rk&
                *(94.00000000000000000D0+7.00000000000000000D0*rk)))) * den(14)
        ckplm(k,4225) = 1.64090331998773000D-18*(-7.93959600000000000D7+rk*(1.01386220000000000D7+rk&
                *(-424429.00000000000000000D0+(6058.00000000000000000D0-11.00000000000000000D0*rk)*rk))) &
                * den(15)
        ckplm(k,4226) = 1.02556457499233000D-19*(7.22292840000000000D8+rk*(-1.09363074000000000D8+rk&
                *(5.83988300000000000D6+rk*(-126486.00000000000000000D0+877.00000000000000000D0*rk)))) * den(16)
        ckplm(k,4227) = 1.08589190293306000D-19*(-3.27015480000000000D8+rk*(5.58043580000000000D7+rk&
                *(-3.47084100000000000D6+(92482.00000000000000000D0-879.00000000000000000D0*rk)*rk))) * den(17)
        ckplm(k,4228) = 2.34606275325044000D-20*(5.91468360000000000D8+rk*(-1.09602394000000000D8+rk&
                *(7.51898300000000000D6+rk*(-225806.00000000000000000D0+2497.00000000000000000D0*rk)))) &
                * den(18)
        ckplm(k,4229) = 7.05582782932463000D-22*(rk-18.D0)*(3.32151540000000000D8+rk&
                *(-4.65210010000000000D7+(2.15538000000000000D6-32999.00000000000000000D0*rk)*rk)) * den(19)
        ckplm(k,4230) = 3.17512252319608000D-22*(rk-18.D0)*(rk-19.D0)*(8.56666000000000000D6+rk&
                *(-804117.00000000000000000D0+18797.00000000000000000D0*rk)) * den(20)
        ckplm(k,4231) = 2.89288941002310000D-21*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)&
                *(6653.00000000000000000D0-313.00000000000000000D0*rk) * den(21)
        ckplm(k,4232) = 6.21971223154966000D-20*(rk-18.D0)*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4233) = -1.51700298330479000D-21*(rk+20.D0)*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4234) = 3.52791391466231000D-23*(rk+20.D0)*(rk+21.D0)&
                *(15523.00000000000000000D0+700.00000000000000000D0*rk) * den(1)
        ckplm(k,4235) = -3.17512252319608000D-22*(rk+20.D0)*(284107.00000000000000000D0+rk&
                *(25770.00000000000000000D0+583.00000000000000000D0*rk)) * den(2)
        ckplm(k,4236) = 3.52791391466231000D-22*(2.56021200000000000D7+rk*(3.49744300000000000D6+rk&
                *(158547.00000000000000000D0+2384.00000000000000000D0*rk))) * den(3)
        ckplm(k,4237) = 2.34606275325044000D-20*(-1.38684000000000000D6-1.00000000000000000D0*rk&
                *(180236.00000000000000000D0+rk*(7725.00000000000000000D0+109.00000000000000000D0*rk))) * den(4)
        ckplm(k,4238) = 6.03273279407255000D-21*(1.49515800000000000D7+rk*(1.79936900000000000D6+rk&
                *(70461.00000000000000000D0+892.00000000000000000D0*rk))) * den(5)
        ckplm(k,4239) = 1.02556457499233000D-19*(-1.96086000000000000D6-1.00000000000000000D0*rk&
                *(210701.00000000000000000D0+rk*(7156.00000000000000000D0+75.00000000000000000D0*rk))) * den(6)
        ckplm(k,4240) = 8.20451659993867000D-19*(453520.00000000000000000D0+rk&
                *(41347.00000000000000000D0+rk*(1115.00000000000000000D0+8.00000000000000000D0*rk))) * den(7)
        ckplm(k,4241) = 1.02556457499233000D-18*(-570360.00000000000000000D0-1.00000000000000000D0*rk&
                *(40862.00000000000000000D0+(rk+723.D0)*rk)) * den(8)
        ckplm(k,4242) = 1.13951619443593000D-19*(6.97722000000000000D6+rk&
                *(344069.00000000000000000D0+(1665.00000000000000000D0-44.00000000000000000D0*rk)*rk)) * den(9)
        ckplm(k,4243) = 5.62921000051348000D-18*(-168060.00000000000000000D0+rk&
                *(-4183.00000000000000000D0+(rk+78.D0)*rk)) * den(10)
        ckplm(k,4244) = -7.09280460064698000D-16*(rk*rk+rk-1400.D0) * den(11)
        ckplm(k,4245) = -5.62921000051348000D-18*(163800.00000000000000000D0+rk&
                *(-4336.00000000000000000D0+(rk-75.D0)*rk)) * den(12)
        ckplm(k,4246) = 1.13951619443593000D-19*(6.63486000000000000D6+rk&
                *(-340607.00000000000000000D0+rk*(1797.00000000000000000D0+44.00000000000000000D0*rk))) &
                * den(13)
        ckplm(k,4247) = 1.02556457499233000D-18*(-530220.00000000000000000D0+rk&
                *(39419.00000000000000000D0+(rk-720.D0)*rk)) * den(14)
        ckplm(k,4248) = 8.20451659993867000D-19*(413280.00000000000000000D0+rk&
                *(-39141.00000000000000000D0+(1091.00000000000000000D0-8.00000000000000000D0*rk)*rk)) * den(15)
        ckplm(k,4249) = 1.02556457499233000D-19*(-1.75724000000000000D6+rk&
                *(196614.00000000000000000D0+rk*(-6931.00000000000000000D0+75.00000000000000000D0*rk))) &
                * den(16)
        ckplm(k,4250) = 6.03273279407255000D-21*(1.32217800000000000D7+rk&
                *(-1.66112300000000000D6+(67785.00000000000000000D0-892.00000000000000000D0*rk)*rk)) * den(17)
        ckplm(k,4251) = 2.34606275325044000D-20*(-1.21422000000000000D6+rk&
                *(165113.00000000000000000D0+rk*(-7398.00000000000000000D0+109.00000000000000000D0*rk))) &
                * den(18)
        ckplm(k,4252) = 3.52791391466231000D-22*(2.22608400000000000D7+rk&
                *(-3.18750100000000000D6+(151395.00000000000000000D0-2384.00000000000000000D0*rk)*rk)) * den(19)
        ckplm(k,4253) = 3.17512252319608000D-22*(rk-19.D0)*(258920.00000000000000000D0+rk&
                *(-24604.00000000000000000D0+583.00000000000000000D0*rk)) * den(20)
        ckplm(k,4254) = 3.52791391466231000D-23*(rk-19.D0)*(rk-20.D0)&
                *(14823.00000000000000000D0-700.00000000000000000D0*rk) * den(21)
        ckplm(k,4255) = 1.51700298330479000D-21*(rk-19.D0)*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4256) = 3.61191186501142000D-23*(rk+21.D0)*(rk+22.D0) * den(0)
        ckplm(k,4257) = -1.67995900698205000D-24*(rk+21.D0)&
                *(8600.00000000000000000D0+389.00000000000000000D0*rk) * den(1)
        ckplm(k,4258) = 1.76395695733116000D-23*(149400.00000000000000000D0+rk&
                *(13687.00000000000000000D0+313.00000000000000000D0*rk)) * den(2)
        ckplm(k,4259) = 5.87985652443719000D-22*(-24798.00000000000000000D0-1.00000000000000000D0*rk&
                *(2209.00000000000000000D0+49.00000000000000000D0*rk)) * den(3)
        ckplm(k,4260) = 2.79293184910766000D-21*(20706.00000000000000000D0+rk&
                *(1759.00000000000000000D0+37.00000000000000000D0*rk)) * den(4)
        ckplm(k,4261) = 2.01091093135752000D-21*(-87276.00000000000000000D0-1.00000000000000000D0*rk&
                *(6889.00000000000000000D0+133.00000000000000000D0*rk)) * den(5)
        ckplm(k,4262) = 5.69758097217963000D-21*(74580.00000000000000000D0+rk&
                *(5279.00000000000000000D0+89.00000000000000000D0*rk)) * den(6)
        ckplm(k,4263) = 1.30230422221249000D-20*(-64890.00000000000000000D0-1.00000000000000000D0*rk&
                *(3917.00000000000000000D0+53.00000000000000000D0*rk)) * den(7)
        ckplm(k,4264) = 1.22091020832421000D-19*(11550.00000000000000000D0+rk&
                *(551.00000000000000000D0+5.00000000000000000D0*rk)) * den(8)
        ckplm(k,4265) = 1.89919365739321000D-19*(-10560.00000000000000000D0-1.00000000000000000D0&
                *(rk+349.D0)*rk) * den(9)
        ckplm(k,4266) = 4.93790350922235000D-20&
                *(49776.00000000000000000D0+(839.00000000000000000D0-7.00000000000000000D0*rk)*rk) * den(10)
        ckplm(k,4267) = 5.92548421106682000D-19*(rk*rk+rk-4410.D0) * den(11)
        ckplm(k,4268) = 4.93790350922235000D-20*(48930.00000000000000000D0-1.00000000000000000D0*rk&
                *(853.00000000000000000D0+7.00000000000000000D0*rk)) * den(12)
        ckplm(k,4269) = 1.89919365739321000D-19*(-10212.00000000000000000D0-1.00000000000000000D0&
                *(rk-347.D0)*rk) * den(13)
        ckplm(k,4270) = 1.22091020832421000D-19*(11004.00000000000000000D0+rk&
                *(-541.00000000000000000D0+5.00000000000000000D0*rk)) * den(14)
        ckplm(k,4271) = 1.30230422221249000D-20&
                *(-61026.00000000000000000D0+(3811.00000000000000000D0-53.00000000000000000D0*rk)*rk) * den(15)
        ckplm(k,4272) = 5.69758097217963000D-21*(69390.00000000000000000D0+rk&
                *(-5101.00000000000000000D0+89.00000000000000000D0*rk)) * den(16)
        ckplm(k,4273) = 2.01091093135752000D-21&
                *(-80520.00000000000000000D0+(6623.00000000000000000D0-133.00000000000000000D0*rk)*rk) * den(17)
        ckplm(k,4274) = 2.79293184910766000D-21*(18984.00000000000000000D0+rk&
                *(-1685.00000000000000000D0+37.00000000000000000D0*rk)) * den(18)
        ckplm(k,4275) = 5.87985652443719000D-22&
                *(-22638.00000000000000000D0+(2111.00000000000000000D0-49.00000000000000000D0*rk)*rk) * den(19)
        ckplm(k,4276) = 1.76395695733116000D-23*(136026.00000000000000000D0+rk&
                *(-13061.00000000000000000D0+313.00000000000000000D0*rk)) * den(20)
        ckplm(k,4277) = 1.67995900698205000D-24*(rk-20.D0)&
                *(8211.00000000000000000D0-389.00000000000000000D0*rk) * den(21)
        ckplm(k,4278) = 3.61191186501142000D-23*(rk-20.D0)*(rk-21.D0) * den(22)
        ckplm(k,4279) = 8.39979503491027000D-25*(-22.00000000000000000D0-1.00000000000000000D0*rk) &
                * den(0)
        ckplm(k,4280) = 8.39979503491027000D-25*(441.00000000000000000D0+20.00000000000000000D0*rk) &
                * den(1)
        ckplm(k,4281) = 5.29187087199347000D-23*(-67.00000000000000000D0-3.00000000000000000D0*rk) &
                * den(2)
        ckplm(k,4282) = 5.87985652443719000D-23*(367.00000000000000000D0+16.00000000000000000D0*rk) &
                * den(3)
        ckplm(k,4283) = 3.91010458875073000D-21*(-24.00000000000000000D0-1.00000000000000000D0*rk) &
                * den(4)
        ckplm(k,4284) = 3.01636639703628000D-21*(103.00000000000000000D0+4.00000000000000000D0*rk) &
                * den(5)
        ckplm(k,4285) = 5.69758097217963000D-21*(-143.00000000000000000D0-5.00000000000000000D0*rk) &
                * den(6)
        ckplm(k,4286) = 6.51152111106244000D-21*(267.00000000000000000D0+8.00000000000000000D0*rk) &
                * den(7)
        ckplm(k,4287) = 7.32546124994524000D-20*(-42.00000000000000000D0-1.00000000000000000D0*rk) &
                * den(8)
        ckplm(k,4288) = 1.89919365739321000D-20*(241.00000000000000000D0+4.00000000000000000D0*rk) &
                * den(9)
        ckplm(k,4289) = 4.93790350922235000D-20*(-117.00000000000000000D0-1.00000000000000000D0*rk) &
                * den(10)
        ckplm(k,4290) = 6.22175842162016000D-18 * den(11)
        ckplm(k,4291) = 4.93790350922235000D-20*(rk-116.D0) * den(12)
        ckplm(k,4292) = 1.89919365739321000D-20*(237.00000000000000000D0-4.00000000000000000D0*rk) &
                * den(13)
        ckplm(k,4293) = 7.32546124994524000D-20*(rk-41.D0) * den(14)
        ckplm(k,4294) = 6.51152111106244000D-21*(259.00000000000000000D0-8.00000000000000000D0*rk) &
                * den(15)
        ckplm(k,4295) = 5.69758097217963000D-21*(-138.00000000000000000D0+5.00000000000000000D0*rk) &
                * den(16)
        ckplm(k,4296) = 3.01636639703628000D-21*(99.00000000000000000D0-4.00000000000000000D0*rk) &
                * den(17)
        ckplm(k,4297) = 3.91010458875073000D-21*(rk-23.D0) * den(18)
        ckplm(k,4298) = 5.87985652443719000D-23*(351.00000000000000000D0-16.00000000000000000D0*rk) &
                * den(19)
        ckplm(k,4299) = 5.29187087199347000D-23*(-64.00000000000000000D0+3.00000000000000000D0*rk) &
                * den(20)
        ckplm(k,4300) = 8.39979503491027000D-25*(421.00000000000000000D0-20.00000000000000000D0*rk) &
                * den(21)
        ckplm(k,4301) = 8.39979503491027000D-25*(rk-21.D0) * den(22)
        ckplm(k,4302) = 1.90904432611597000D-26 * den(0)
        ckplm(k,4303) = -4.19989751745513000D-25 * den(1)
        ckplm(k,4304) = 4.40989239332789000D-24 * den(2)
        ckplm(k,4305) = -2.93992826221859000D-23 * den(3)
        ckplm(k,4306) = 1.39646592455383000D-22 * den(4)
        ckplm(k,4307) = -5.02727732839380000D-22 * den(5)
        ckplm(k,4308) = 1.42439524304491000D-21 * den(6)
        ckplm(k,4309) = -3.25576055553122000D-21 * den(7)
        ckplm(k,4310) = 6.10455104162104000D-21 * den(8)
        ckplm(k,4311) = -9.49596828696606000D-21 * den(9)
        ckplm(k,4312) = 1.23447587730559000D-20 * den(10)
        ckplm(k,4313) = -1.34670095706064000D-20 * den(11)
        ckplm(k,4314) = 1.23447587730559000D-20 * den(12)
        ckplm(k,4315) = -9.49596828696606000D-21 * den(13)
        ckplm(k,4316) = 6.10455104162104000D-21 * den(14)
        ckplm(k,4317) = -3.25576055553122000D-21 * den(15)
        ckplm(k,4318) = 1.42439524304491000D-21 * den(16)
        ckplm(k,4319) = -5.02727732839380000D-22 * den(17)
        ckplm(k,4320) = 1.39646592455383000D-22 * den(18)
        ckplm(k,4321) = -2.93992826221859000D-23 * den(19)
        ckplm(k,4322) = 4.40989239332789000D-24 * den(20)
        ckplm(k,4323) = -4.19989751745513000D-25 * den(21)
        ckplm(k,4324) = 1.90904432611597000D-26 * den(22)
    enddo
    return
    end
