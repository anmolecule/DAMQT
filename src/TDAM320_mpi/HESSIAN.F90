!> @file HESSIAN.F90
!> @details Calculates the eigenvalues and eigenvectors of derivatives
      SUBROUTINE HESSIAN(IRANK,IRNK,D,V,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot)  
!***************************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NP = 5,NTMP=6*6)
      DIMENSION D(NP),V(NP,NP),H(NP,NP),E(NP,NP)
      REAL*8 :: VTOT, DRVXTOT, DRVYTOT, DRVZTOT, DXXTOT, DXYTOT, DXZTOT
      REAL*8 :: DYXTOT, DYYTOT, DYZTOT, DZXTOT, DZYTOT, DZZTOT

      DO 61 i=1,3
      DO 62 j=1,3
      H(i,j) = 0.0
62    CONTINUE
61    CONTINUE
!        HESSIAN MATRIX AT THE POINT x1,y1,z1
      H(1,1) = dxxtot ; H(1,2) = dxytot ; H(1,3) = dxztot
      H(2,1) = H(1,2) ; H(2,2) = dyytot ; H(2,3) = dyztot
      H(3,1) = H(1,3) ; H(3,2) = H(2,3) ; H(3,3) = dzztot 

      TRACE = H(1,1) + H(2,2) + H(3,3)

!      WRITE(6,211)
!      WRITE(6,*)
!      WRITE(6,*)H(1,1),'  ',H(1,2),'  ',H(1,3)
!      WRITE(6,*)
!      WRITE(6,*)H(2,1),'  ',H(2,2),'  ',H(2,3)
!      WRITE(6,*)'  '
!      WRITE(6,*)H(3,1),'  ',H(3,2),'  ',H(3,3)
!      WRITE(6,*)' '
!      WRITE(6,212)TRACE

      N=3
      DO 51 i=1,3
      DO 52 j=1,3
      E(i,j) = H(i,j)
52    CONTINUE
51    CONTINUE

      CALL JACOBI (E,N,NP,D,V,NROT)

!      WRITE(6,213)
!      WRITE(6,*)'  '
!      WRITE(6,*)(D(i),i=1,N)

      IPC = 0
      IRNK = 3
      DO i = 1,3
        IF(D(i).GT.0) IPC = IPC + 1
        IF(DABS(D(i)).LT.1.0D-8) IRNK = IRNK - 1
      ENDDO

      IRANK = 2*IPC - IRNK
 
      IF(2*IPC-IRNK.GT.0) THEN
!        WRITE(*,215)IRNK,2*IPC-IRNK
      ELSE
!        WRITE(*,214)IRNK,2*IPC-IRNK
      ENDIF
    
!      WRITE(6,*)'  '
!      WRITE(6,*)'  '

!  -------------------------------------------------------
17    FORMAT(7F20.7)
211   FORMAT(/,' HESSIAN MATRIX AT CP : ')
212   FORMAT(//,' TRACE OF THE HESSIAN MATRIX : ',F25.16)
213   FORMAT(//,' EIGEN VALUES : ')
214   FORMAT(//,' CHARACTERIZATION OF CP :     (',I1,',',I2,')')
215   FORMAT(//,' CHARACTERIZATION OF CP :     (',I1,',+',I1,')')

      ENDSUBROUTINE 
!****************************************************************************

!****************************************************************************
!> @details Hessian employes jacobi subroutine
         SUBROUTINE JACOBI (A,N,NP,D,V,NROT)
!****************************************************************************
! THIS SUBROUTINE COMPUTES ALL EIGEN VALUES OF A REAL SYMMETRIC MATRIX A
! OF SIZE N X N, STORED IN A Np X Np PHYSICAL ARRAY. D RETURNS THE EIGEN
! VALUES OF A IN ITS FIRST N ELEMENTS. V IS A MATRIX WITH THE SAME LOGICAL
! AND PHYSICAL DIMENSION AS A.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),B(NMAX),Z(NMAX),D(NP),V(NP,NP)

        do 12 ip = 1,n
        do 11 iq = 1,n
                v(ip,iq) = 0.0d0
11      continue
                v(ip,ip) = 1.0d0
12      continue
        do 13 ip = 1,n
                b(ip) = a(ip,ip)
                d(ip) = b(ip)
                z(ip) = 0.0d0
13      continue
        nrot = 0
        do 24 i = 1,50
                sm = 0.0d0
                do 15 ip = 1,n-1
                        do 14 iq = ip + 1,n
                                sm = sm + dabs(a(ip,iq))
14                      continue
15              continue
                if(dabs(sm).le.1.0d-20)return

                if(i.lt.4)then
                        thresh = 0.2d0*sm/n**2
                else
                        thresh = 0.0d0
                endif

                do 22 ip = 1,n-1
                        do 21 iq = ip+1,n
                                g = 100.0d0*dabs(a(ip,iq))

                                if((i.gt.4).and.(dabs(d(ip))+g.eq.dabs(d(ip)))&
                                & .and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
                                        a(ip,iq) = 0.0d0
                                else if(dabs(a(ip,iq)).gt.thresh)then
                                        h = d(iq) - d(ip)

                                        if(dabs(h)+g.eq.dabs(h))then
                                                t = a(ip,iq)/h
                                        else
                                                theta  = 0.5d0*h/a(ip,iq)
                                                t = 1.0d0/(dabs(theta) + dsqrt(1.0d0 + theta*theta))
                                                if(theta.lt.0.0d0) t= -t
                                        endif

                                        c = 1.0d0/dsqrt(1 + t*t)
                                        s = t*c
                                        tau = s/(1 + c)
                                        h = t*a(ip,iq)
                                        z(ip) = z(ip) - h
                                        z(iq) = z(iq) + h
                                        d(ip) = d(ip) - h
                                        d(iq) = d(iq) + h
                                        a(ip,iq) = 0.0d0
                                        do 16 j = 1,ip-1
                                                g = a(j,ip)
                                                h = a(j,iq)
                                                a(j,ip) = g - s*(h+g*tau)
                                                a(j,iq) = h+s*(g-h*tau)
16                                      continue
                                        do 17 j = ip+1,iq-1
                                                g = a(ip,j)
                                                h = a(j,iq)
                                                a(ip,j) = g - s*(h+g*tau)
                                                a(j,iq) = h + s*(g-h*tau)
17                                      continue
                                        do 18 j = iq+1,n
                                                g = a(ip,j)
                                                h = a(iq,j)
                                                a(ip,j) = g - s*(h+g*tau)
                                                a(iq,j) = h + s*(g-h*tau)
18                                      continue
                                        do 19 j = 1,n
                                                g = v(j,ip)
                                                h = v(j,iq)
                                                v(j,ip) = g - s*(h+g*tau)
                                                v(j,iq) = h + s*(g-h*tau)
19                                      continue
                                        nrot = nrot + 1
                                endif

21                      continue
22              continue
                do 23 ip = 1,n
                        b(ip) = b(ip) + z(ip)
                        d(ip) = b(ip)
                        z(ip) = 0.0d0
23              continue
24     continue
       return
       end
!****************************************************************************
!
!****************************************************************************
