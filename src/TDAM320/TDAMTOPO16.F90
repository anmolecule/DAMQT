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
!> @file TDAMTOPO15.F90

!> @authors Anmol Kumar, Sachin D. Yeole and Shridhar R. Gadre
!> @date 10-07-2015
!> @details The main program is the driver file for performing various facets of 
!> topography of molecular electron density(MED) and molecular electrostatic potential(MESP).
!> The key features of this program are:
!> \a a. Mapping critical points (CPs).
!> \a b. Determination of Euler characteristics for MED and MESP.
!> \a c. Gradient path determination for MED and MESP
!> \a d. Three dimensional basin generation for MED and MESP.
!>
!> @param [in] Projectname Base name provided for input files.
!> @param [in] Filename Base name provided for output files. 
!! In case filename is not defined, projectname will be used for output files as well.
!> @param [in] med Enables the topography for molecular electron density.
!> @param [in] mesp Enables the topography for molecular electrostatic potential.
!> @param [in] topomanual Enables MED and MESP value evaluation at user defined points.
!> @param [in] topograph Enables mapping of critical points. 
!> @param [in] lmaxi Value of l used for accurate determination of scalar field. Default is 15.
!> @param [in] boxl Size of the box in au around the guess point within which CP will be searched. 
!! User may vary the boxl value till Euler characteristic is statisfied. The recommended range is 0.8 to 3.0.
!! Default value is 1.2
!> @param [in] cnvg The threshold value of norm of square of gradient to decide if the optimizer has found the CP. 
!> @param [in] gradpath Enables atomic interaction line determination for MED and MESP. 
!! In case of MESP, shortest path connecting various CPs are also shown.
!> @param [in] boxg The value provided is used to create a box around molecule within which gradient path will be 
!! generated, especially in case of MESP.
!> @param [in] basin Enables three dimensional creation of basin of each atom based on MED and MESP.
!> @param [in] boxb The value provided is used to create a box around molecule within which basin will be 
!! generated, especially in case of MESP.
!> @param [in] atoms Enables the user to show separate basin for user-defined atom/s.
!> @param [in] nselbasin Number of atoms for which separate basin is to be shown.
!> @param [in] iselect The index of atoms (predefined according to .xyz file) for which separate basin is to be shown.
!> @param [in] addguess Enables user to provide additional guess points.
!! to find field values at user defined points.
 
    Program Topodriver
    USE DAMINITIAL_T
    USE DAM320_DATA_T, ONLY: lvalence
    implicit none
    integer(kint)::iden,ierr,i, j  ! Modified by Rafa
    real :: tarray(2), tiempo, dtime
    logical :: there
    namelist / options / med, mesp, topograph,gradpath, &
                & basin, lgradient, lderiv2, largo, lexact, lvalence, &
                & wireframe,solidsurf,addguess,exdraw,cnvg, &
                & boxl, boxg,boxb,boxt,stepszt,angle,exln,lmaxi, &
                & filename,guessfile,iswindows,fdisp, drcutcp, ncntguess, rcntguess

    tiempo = dtime(tarray)
    lvalence   = .false.
    med        = .false.
    mesp       = .true.
    topograph  = .true.
    gradpath   = .false.
    basin      = .false.
    lgradient  = .true.
    lderiv2    = .true.
    largo      = .false.		! IF .TRUE. LONG-RANGE POTENTIAL
    lexact     = .false.		! IF .TRUE. "EXACT" POTENTIAL IS TABULATED
    addguess   = .false.
    wireframe  = .true.
    solidsurf  = .true.
    exdraw     = .false.
    iswindows  = .false.
    cnvg       = 4.0e-12
    boxl       = 1.2
    boxg       = 6.0
    boxb       = 4.5
    boxt       = 5.0
    stepszt    = 0.2
    angle      = 10.0    ! modified by Rafa 2019
    exln       = 3.0
    lmaxi      = 15
    fdisp      = 0.02
    filename   = ""
    guessfile  = ""
    drcutcp    = 1.0e-5
    ncntguess = 0
    rcntguess = 0.d0
    read(5,options)     ! Read from standard input
    read(5,*) projectname
if (.not.med.and..not.mesp) then
		call error(1,"Field option not provided")
	end if
! Decides if the user has provided projectname with full path.
	if (iswindows) then
		dirsep = "\\"
		i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
		if (i .eq. 0) then	! This is intended for MinGW, whose directory separator in windows is also /
			dirsep = "/"
			i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
		endif
	else
		dirsep = "/"
		i = index(projectname,dirsep,.true.)	! Checks position of last directory name separator
	end if  
	if (len_trim(filename).eq.0) then
		filename	 = projectname
	else
		filename = projectname(1:i)//trim(filename)
	endif
! 	Code introdced by Rafa for files deletion in windows 
	if (iswindows) then
		filenamewin = filename
		do j = 1, len_trim(filenamewin)
			i = index(filenamewin,'/') ; if (i == 0) exit
			filenamewin = filenamewin(:i-1) // '\' // filenamewin(i+1:)
		end do
	endif
! ----------      End of change made by Rafa   ----------------------------  
	mesp=.not.med 
      
!       open(iout,file=trim(filename)//"-DAMTOPO320.out",status="replace")  ! Sealed by Rafa
      
	if (gradpath.or.basin) THEN  
		if (med) inquire(file=trim(filename)//"-cps-d.xyz",exist=there)
		if (mesp) inquire(file=trim(filename)//"-cps-v.xyz",exist=there)
		if (.not.there) topograph=.true.
	end if      
      
	if (basin) gradpath=.true.
!	if (basin) wireframe=.true. ! modified by Rafa 2019
!	if (basin) solidsurf=.true. ! modified by Rafa 2019
	  
	call titles    ! projectname (for input files) and filename (for output files) passed through module DAMINITIAL_T

!      call gdam
      
    if (med) call predamden
    if (mesp) call predampot

    call flush(iout)

    write(iout,"(a)")"Nuclear Charges and Cartesian Coordinates:"
    write(iout,"(a)")"-------------------------------------------------------------------------------"
    write(iout,"(t1,a,t33,a,t52,a,t71,a)") "ATOM","x", "y","z"
    write(iout,"(a)")"-------------------------------------------------------------------------------"
     
!	  open(99,file=trim(filename)//".xyz",form='formatted', iostat=ierr)

!      if (ierr .ne. 0) then
!		write(6,"('Cannot open file ', a)") trim(projectname)//".xyz"
!	  else
!	    write(99,*) ncen
!	    write(99,*)
	    do i = 1 , ncen
		    write(iout,"(a2,i3,t27,f12.6,t44,f12.6,t64,f12.6)") &
                           & atmnam(i),i,rcen(1,i),rcen(2,i),rcen(3,i)
!		    write(99,"(a2,3e15.10)") atmnam(i), &
!               & rcen(1,i)*0.529177249d0,rcen(2,i)*0.529177249d0,rcen(3,i)*0.529177249d0
	    end do
!	    close(99)
!	   endif
    write(iout,"(a)")"-------------------------------------------------------------------------------"
    write(iout,"(a)")"                                                                               "
      
      if (topograph) call topo
  
!======Center of mass====================

    cox=0; coy=0; coz=0 ; tmass=0
    do i=1,ncen
        tmass=tmass+atwt(i)
        cox=cox+atwt(i)*rcen(1,i)
        coy=coy+atwt(i)*rcen(2,i)
    coz=coz+atwt(i)*rcen(3,i)
    end do
    cox=cox/tmass; coy=coy/tmass; coz=coz/tmass

    if (gradpath) call gradpath_s

    if (basin) then
        if (angle.le.12.0) angle = 5.0
        if (angle.gt.12.0.and.angle.le.20.0) angle = 15.0
        if (angle.gt.20.0.and.angle.le.35.0) angle = 30.0
        if (angle.gt.35.0.and.angle.le.50.0) angle = 45.0
        if (angle.gt.50.0.and.angle.lt.90.0) angle = 90.0
        if (angle.ge.90.0) angle = 90.0
        call basin_s
    end if


    tiempo = dtime(tarray)
    write(iout,*)""
    write(iout,"(1x,'Timing in seconds of processor ',' (user, system, total):',/5x,'(', &
        &e12.5,',',e12.5,',',e12.5')')") &
        &tarray(1), tarray(2), tarray(1)+tarray(2)
    write(iout,*)""
    close(iden)
    close(iout)
    stop
    end

!******************************************
!> @details Creates guess points for MED and MESP topography.
  SUBROUTINE TOPO
!******************************************      
    USE DAMINITIAL_T
    IMPLICIT NONE
    CHARACTER::dummy*10
    INTEGER(KINT):: NPTS,NPNTRC,nn_iln,np_iln,ierr,ig,ngspt
    REAL(KREAL):: xmin,ymin,zmin,xmax,ymax,zmax,xgpt,ygpt,zgpt
    REAL(KREAL), ALLOCATABLE:: SUMDER(:,:,:),FVAL(:,:,:)  !SUMDERFR(:,:,:)
    REAL(KREAL) :: X1,Y1,Z1,den,vtot,drvxtot,drvytot,drvztot,dxxtot,&
                 & dxytot,dxztot,dyytot,dyztot,dzztot
    REAL(KREAL) :: Sr,bx,by,bz,cdx,cdy,cdz
    REAL(KREAL):: c1,c2,c3,c4,c5,c6,c7
    INTEGER(KINT):: i,j,k
    INTEGER(KINT):: nstepx,nstepy,nstepz
    INTEGER(KINT) :: NCOUNTER
    CHARACTER(4) :: IRANKD
    CHARACTER(2) :: gsym
    logical :: lexist     ! Introduced by Rafa
!------------------------------------------------------------------------------------------------------------------------
! In case one wants to provide certain points where the value of scalar
! field are to be evaluated
!------------------------------------------------------------------------------------------------------------------------
!      if (topomanual) then
!            write(6,*)"Enter the number of points on which you want &
!                    & to calculate ESP"
!            npts = nguesspoints
!            write(iout,*)"MED/MESP and gradient at user defined points"
!            write(iout,126)"X","Y","Z","Val","dx","dy","dz"
!126         format(t9,a,t26,a,t44,a,t63,a,t87,a,t110,a,t132,a)       
!  -------- End of  change made by Rafa -----------------------------------  

!        open(unit=125, file = trim(filename)//".pv",access="append", status="replace")  
!        do i = 1,npts 
!            if(med) then
!               CALL DAMDEN(den,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xpt(i),ypt(i),zpt(i))
!               write(125,420)xpt(i),ypt(i),zpt(i),vtot,drvxtot,drvytot,drvztot
!            elseif(mesp) then
!               CALL DAMPOT(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xpt(i),ypt(i),zpt(i))
!               write(125,420)xpt(i),ypt(i),zpt(i),vtot,drvxtot,drvytot,drvztot
!            endif
!        end do
!420     format(3(1x,e17.10),1x,e22.15,3(1x,e22.15))
!        close(125)
!        close(iout)
!        call system ("cat "//trim(filename)//".pv >> "//trim(filename)//"-DAMTOPO320.out")
!        call system ("rm -f "//trim(filename)//".pv")
!
!        return 
!      end if
!----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------
      
      if (med) write(iout,"(a,I3)") "Euler Characteristic = -1"
              
      open ( unit = 31, file = trim(filename)//".ts",access="append", status="replace") 

      xmin = 0 ; ymin = 0 ; zmin = 0
      xmax = 0 ; ymax = 0 ; zmax = 0
     
      xmin = rcen(1,1); ymin = rcen(2,1); zmin = rcen(3,1) 
      xmax = rcen(1,1); ymax = rcen(2,1); zmax = rcen(3,1) 
      
      do 12 I=1,ncen
            if(rcen(1,i).lt.xmin) xmin = rcen(1,i) 
            if(rcen(1,i).gt.xmax) xmax = rcen(1,i) 
            if(rcen(2,i).lt.ymin) ymin = rcen(2,i) 
            if(rcen(2,i).gt.ymax) ymax = rcen(2,i) 
            if(rcen(3,i).lt.zmin) zmin = rcen(3,i) 
            if(rcen(3,i).gt.zmax) zmax = rcen(3,i) 
12    continue
 
      xmin=xmin-boxt-stepszt;ymin=ymin-boxt-stepszt;zmin=zmin-boxt-stepszt
      xmax=xmax+boxt+stepszt;ymax=ymax+boxt+stepszt;zmax=zmax+boxt+stepszt
            ! Adding 5 au to maximum and subtracting 5 au from minimum
            ! ensures that the guess points are formed within these limits
            ! and it is unlikely that a CP will be found beyond 5 au
            ! distance.
      bx = (xmax) - (xmin) 
      by = (ymax) - (ymin) 
      bz = (zmax) - (zmin) 
          ! Bx, By and Bz defines the box to be used for guess point
          ! generation
            !++++++++++++++++++++++++++++++
      if (mesp) then
          Cdx = ((xmax+1.20d0*boxl) + (xmin-1.20d0*boxl))/2.0d0
          Cdy = ((ymax+1.20d0*boxl) + (ymin-1.20d0*boxl))/2.0d0
          Cdz = ((zmax+1.20d0*boxl) + (zmin-1.2d00*boxl))/2.0d0
          Sr = max((xmax+1.20d0*boxl-Cdx),(ymax+1.20d0*boxl-Cdy),(zmax+1.20d0*boxl-Cdz))
          ! Sr is the radius thus found for the sphere on which 
          ! Euler characteristic is to be determined. Cdx, Cdy and Cdz
          ! is the center of the sphere such that all the CPs will be 
          ! inside the sphere thus defined.
!          print*,Sr,cdx,cdy,cdz,xmax,xmin,ymax,ymin,zmax,zmin
          CALL findeuler(Sr)
      end if     
!          ! An arbitrary grid size of 0.2 au has been choosen for
!          ! creating guess points as it allows optimal time consumption
!          ! and guess point produced are very close to actual CP
      nstepx = bx/stepszt
      nstepy = by/stepszt 
      nstepz = bz/stepszt

      write(iout,"(a)") "Guess points of MESP are created using &
                  & gradient calculation on grid points"
      write(iout,"(a)")"Dimensions of the grid (in a.u.) in & 
                  & x, y, and z are"
      write(iout,"(3(f6.1,2x))") bx,by,bz
      write(iout,"(a)") "Number of steps in x, y, and z"
      write(iout,"(3(I10,2x))") nstepx,nstepy,nstepz
      write(iout,"(a,2x,I10)") "Total number of steps", nstepx*nstepy*nstepz
      write(iout,"(a,2x,f6.2)") "Stepsize of the grid", stepszt     
        
      allocate (sumder(nstepx,nstepy,nstepz))
      allocate (fval(nstepx,nstepy,nstepz))
          ! Storing the first partial derivatives in array sumder. It is
          ! later used to filter those grid values which are having
          ! smaller magnitude of derivatives than their neighbours.

          ! Following stores the magnitude of first partial derivatives 
          ! at all the grid points which fall inside the box, in an array called sumder.  
      lderiv2=.false.
      do 111 i=1,nstepx 
         do 112 j=1,nstepy 
            do 113 k= 1,nstepz 
                x1 = xmin + stepszt * (i - 1) 
                y1 = ymin + stepszt * (j - 1) 
                z1 = zmin + stepszt * (k - 1)
                if (med) then
                    CALL DAMDEN(den,drvxtot,drvytot,drvztot,dxxtot, &
                        & dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
                    fval(i,j,k) = dabs(den)
                elseif (mesp) then
                    CALL DAMPOT(vtot,drvxtot,drvytot,drvztot, &
                        & dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
                    fval(i,j,k) = dabs(vtot)
                end if
                sumder(i,j,k) = drvxtot**2 + drvytot**2 + drvztot**2
113         continue
112      continue
111   continue
      lderiv2=.true.
!C****************************************************************************
!           ! Following filters the derivatives which are of smaller
!           ! magnitude relative to their neighbouring guess points.
!           ! Stores such grid point coordinates in Transition (ts) file.
            ! These serve as guess points to MESP.   

      do 211 i=2,nstepx-1
           do 212 j=2,nstepy-1
               do 213 k= 2,nstepz-1
                   if (fval(i,j,k).le.0.00005) cycle 
                   c1=sumder(i,j,k) 
                   c2=sumder(i+1,j,k) 
                   c3=sumder(i-1,j,k)
                   c4=sumder(i,j+1,k) 
                   c5=sumder(i,j-1,k) 
                   c6=sumder(i,j,k+1)
                   c7=sumder(i,j,k-1) 
                   if(c1.lt.c2.and.c1.lt.c3.and.c1.lt.c4.and.c1.lt.c5.and.c1.lt.c6.and.c1.lt.c7) then 
                        x1 = xmin + stepszt * (i - 1)
                        y1 = ymin + stepszt * (j - 1)
                        z1 = zmin + stepszt * (k - 1)
                        write(31,"(3(f9.4,2x))")x1,y1,z1
                   endif
213             continue
212        continue
211   continue
      close(31)
      deallocate (sumder,fval)
      call guessmedcp     ! Creates extra guess points between nuclear positions.
      call aroundnucleus  ! Creates extra guess points around nuclear positions.
      call sortgcps       ! Sorts the CPs is ts file and stores in gcp file.
      if (iswindows) then
            lexist = .false.	! Checks whether the file .ts exists
            inquire(file=TRIM(filenamewin)//".ts", exist=lexist, iostat=ierr)
            if (lexist) call system("del /f /q  """//TRIM(filenamewin)//".ts""")
      else
            call system ("rm -f "//trim(filename)//".ts")
      endif

    if (addguess) then
        open(33,file=trim(filename)//".gcp",access="append",status="old",iostat=ierr)
        if (ierr.ne.0) then
            write(6,"('Error reading file ', a, ' with guess points. ierr = ', i5)") trim(adjustl(guessfile)), ierr
        else
            write(6,"(/10x,'Additional guess points',/)")
            do i = 1, min(ncntguess,mxguess)    ! Guess points given in namelist
                write(33,"(3(2x,f9.5))") rcntguess(1,i), rcntguess(2,i), rcntguess(3,i)
                write(6,"('extra point (in au) =', 3(1x,e22.15))") rcntguess(1,i), rcntguess(2,i), rcntguess(3,i)
            enddo
            if (len_trim(guessfile) > 0) then   ! Gess points given in file guessfile
                inquire(file=trim(adjustl(guessfile)), exist=lexist, iostat=ierr)
                if (.not. lexist) then
                    write(6,"(/'File ', a, ' with guess points does not exist',/)") trim(adjustl(guessfile))
                else
                    open(35,file=trim(adjustl(guessfile)),status="old",iostat=ierr)
                    if (ierr.ne.0) then
                        write(6,"('Cannot open file ', a, ' with guess points. ierr = ', i5)") trim(adjustl(guessfile)), ierr
                    else
                        write(6,"(/'Reads file', a,' with extra guess points',/)") trim(adjustl(guessfile))
                        do while (ierr .eq. 0)
                            read(35,*, iostat=ierr) xgpt, ygpt, zgpt
                            if (ierr .eq. 0) then
                                write(33,"(3(2x,f9.5))") xgpt,ygpt,zgpt
                                write(6,"('extra point (in au) =', 3(1x,e22.15))") xgpt,ygpt,zgpt
                            endif
                        end do
                    endif
                    close(35)
                endif
            endif
            close(33)
        endif
    endif
    ! The variables used in above subroutines guessmedcp and
    ! aroundnucleus are not used further in the program.
    call optdriver
    return
    end subroutine

!-------------------------------------------------------------------------------------
!> @details Euler characteristic is determined by analysing the number of
!> islands with negative and positive MESP on a large spherical grid around the
!> molecule.
!***************************************************

    subroutine findeuler (Sr)
    !***************************************************
    USE DAMINITIAL_T
    implicit none
    real(kreal), dimension(37,19)::x_s,y_s,z_s
    integer(kint),dimension(800)::nil1,nil2,nil3
    integer(kint),dimension(800)::zil1,zil2,zil3
    integer(kint),dimension(800)::pil1,pil2,pil3
    integer(kint)::i,j,itheta,jphai,k,l,m,ist  !,ifn
    REAL(KREAL) :: X1,Y1,Z1,den,vtot,drvxtot,drvytot,drvztot,dxxtot,&
                 & dxytot,dxztot,dyytot,dyztot,dzztot
    real(kreal):: Sr
    integer(kint)::nn_iln,np_iln,maxn_iln,maxp_iln

    i = 0; j = 0
    do itheta=0,360,10
        i=(itheta/10)+1
        do jphai=0,180,10
            j=(jphai/10)+1
            x_s(i,j) = Sr*sin(jphai*p_i/180.0d0)*cos(itheta*p_i/180.0d0) ! Storing the x y z coordinate relative to points on spherical grid
            y_s(i,j) = Sr*sin(jphai*p_i/180.0d0)*sin(itheta*p_i/180.0d0)
            z_s(i,j) = Sr*cos(jphai*p_i/180.0d0)
!            write(120,"(a,3(2x,f10.5))") "x",x_s(i,j),y_s(i,j),z_s(i,j)
        end do
    end do
      
    LGRADIENT   = .FALSE.
    LDERIV2     = .FALSE.


    k=0; l=0; m=0
    do i=1,37
        do j=1,19
            if (i.gt.1.and.j.eq.1) cycle ! to ensure that when phai is zero
                                         ! the point correposnding to different
                                         ! theta are considered only once
            if (i.gt.1.and.j.eq.19) cycle ! to ensure that when phai is 180 degrees
                                          ! the point corresponding to different
                                          ! theta are considered only once
            CALL DAMPOT(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot, &
                              &  dyytot,dyztot,dzztot,x_s(i,j),y_s(i,j),z_s(i,j))
            if (vtot.le.0.0d0) then ! if MESP value is negative
                k = k + 1              ! increase the value of index (k) by 1
                nil1(k)=i          ! Store the corresponding index values of coordinates in another array nil1 and nil2
                nil2(k)=j
                nil3(k)=0          ! nil3 is an array with signature value 0
!               write(121,"(a,4(2x,f10.5))") "x",x_s(i,j),y_s(i,j),z_s(i,j), vtot
!               write(121,"(4(I4,2x))") nil1(k),nil2(k),nil3(k),k
            elseif (vtot.gt.0.0d0) then ! if MESP value is positive
                m = m + 1              ! increase the value of index (m) by 1
                pil1(m)=i              ! Store the corresponding index values of coordinates in another array pil1 and pil2
                pil2(m)=j
                pil3(m)=0              ! pil3 is an array with signature value 0
!               write(122,"(a,4(2x,f10.5))") "z",x_s(i,j),y_s(i,j),z_s(i,j), vtot
!               write(122,"(4(I4,2x))") pil1(m),pil2(m),pil3(m),m
            end if
        end do
    end do
    close (121)
    close (122)

    LGRADIENT   = .TRUE.
    LDERIV2     = .TRUE.
!        write(121,*)
!        write(122,*)
    if (k.eq.631) then  ! In case the molecule is an anion, all the
                          ! grid points will possess negative esp value, thus no further
                          ! analysis is needed. 631 is maximum number of
                          ! grid points formed from above loop.
          nn_iln=1
          np_iln=0
          go to 631
    end if
    if (m.eq.631) then  ! In case molecule is cation similar method is acquired
          np_iln=1
          nn_iln=0
          go to 631
    end if


! Calculating number of negative islands
    nn_iln=0 ! Initializing number of negative island
    maxn_iln=10 ! Taking maximum number of distinct islands possible to be 10
    do 10 ist=1,maxn_iln
          if (nil3(ist).eq.0) then
              nn_iln=nn_iln+1
            do 20 i=ist,k      ! loop starts from index where it finds nil3 to be zero
               if (i.eq.ist) then  ! if it is first cycle 
                  do 30 j=i+1,k      ! Run next cycle
                        if (abs(nil1(i)-nil1(j)).eq.0.and.abs(nil2(i)-nil2(j)).eq.1 &
                        &  .or.abs(nil1(i)-nil1(j)).eq.1.and.abs(nil2(i)-nil2(j)).eq.0  &
                        &  .or.abs(nil1(i)-nil1(j)).eq.36.and.abs(nil2(i)-nil2(j)).eq.0) then
                         ! if there is a difference of 1 or 36 between the nil1
                         ! array elements when nil2 array elements are the same
                         ! or the other way around
                           nil3(i)=nn_iln   ! Change the nil3 signature array to nth number of negative island 
                           nil3(j)=nn_iln   ! for both the loop elements
!                           write(121,"(4(I4,2x))") nil1(i),nil2(i),nil3(i),i
!                           write(121,"(4(I4,2x))") nil1(j),nil2(j),nil3(j),j
                        end if
30                continue   
               elseif (i.gt.ist.and.nil3(i).eq.nn_iln) then 
                  ! if it is not the first cycle then go to those elements which have already been
                  ! designated with island number  
                  do 40 j=i+1,k
                        if (abs(nil1(i)-nil1(j)).eq.0.and.abs(nil2(i)-nil2(j)).eq.1 &
                        &  .or.abs(nil1(i)-nil1(j)).eq.1.and.abs(nil2(i)-nil2(j)).eq.0  &
                        &  .or.abs(nil1(i)-nil1(j)).eq.36.and.abs(nil2(i)-nil2(j)).eq.0) then
                           nil3(i)=nn_iln
                           nil3(j)=nn_iln
!                           write(121,"(4(I4,2x))") nil1(i),nil2(i),nil3(i),i
!                           write(121,"(4(I4,2x))") nil1(j),nil2(j),nil3(j),j
                        end if
40                continue   
               end if
20          continue
            ! At the end of this loop all the grid elements which are connected
            ! with each other have been alloted same island number
         end if
         ! Next cycle will increase the island number by 1   
10    continue
       
!144   continue    
      ! Calculating number of positive islands 
      np_iln=0 ! Initializing number of positive island
      maxp_iln=10 ! Taking maximum number of distinct islands possible to be 10
      do 50 ist=1,maxp_iln
!         do ifn=1,m
!            if (pil3(ifn).eq.0) then !if the nil3 array initially set to zero is zero
!                np_iln=np_iln + 1
!                exit
!            else if (ifn.eq.m.and.pil3(ifn).ne.0) then
!                go to 631
!            end if    
!         end do   
!         print*,ifn
          if (pil3(ist).eq.0) then
              np_iln=np_iln+1

            do 60 i = ist , m      ! loop starts from index where it finds pil3 to be zero
               if (i.eq.ist) then  ! if it is first cycle 
                  do 70 j=i+1,m      ! Run next loop
                        if (abs(pil1(i)-pil1(j)).eq.0.and.abs(pil2(i)-pil2(j)).eq.1 &
                        &  .or.abs(pil1(i)-pil1(j)).eq.1.and.abs(pil2(i)-pil2(j)).eq.0  &
                        &  .or.abs(pil1(i)-pil1(j)).eq.36.and.abs(pil2(i)-pil2(j)).eq.0) then
                         ! if there is a difference of 1 or 36 between the pil1
                         ! array elements when pil2 array elements are the same
                         ! or the other way around
                           pil3(i)=np_iln   ! Change the pil3 signature array to nth number of positive island 
                           pil3(j)=np_iln   ! for both the loop elements
!                              write(122,"(4(I4,2x))") pil1(i),pil2(i),pil3(i),i
!                              write(122,"(4(I4,2x))") pil1(j),pil2(j),pil3(j),j
                        end if
70                continue   
               elseif (i.gt.ist.and.pil3(i).eq.np_iln) then 
                  ! if it is not the first cycle then go to those elements which have already been
                  ! designated with island number  
                  do 80 j=i+1,m
                        if (abs(pil1(i)-pil1(j)).eq.0.and.abs(pil2(i)-pil2(j)).eq.1 &
                        &  .or.abs(pil1(i)-pil1(j)).eq.1.and.abs(pil2(i)-pil2(j)).eq.0  &
                        &  .or.abs(pil1(i)-pil1(j)).eq.36.and.abs(pil2(i)-pil2(j)).eq.0) then
                           pil3(i)=np_iln
                           pil3(j)=np_iln
!                              write(122,"(4(I4,2x))") pil1(i),pil2(i),pil3(i),i
!                              write(122,"(4(I4,2x))") pil1(j),pil2(j),pil3(j),j
                        end if
80                continue   
               end if
60          continue
            ! At the end of this loop all the grid elements which are connected
            ! with each other have been alloted same island number
         end if    
50    continue
631   continue
      write(iout,*)  
      write(iout,"(1x,a,I3)") "Number of negative islands = ", nn_iln
      write(iout,"(1x,a,I3)") "Number of positive islands = ", np_iln
      write(iout,"(1x,a,I3)") "Euler Characteristic = ", nn_iln - np_iln
      return
      end
!============================================================================================

!=====================================================================================
!> @details Creates guess points for MED and MESP topography by using mid points
!! of two atoms.
      SUBROUTINE GUESSMEDCP  
!*************************************      
      USE DAMINITIAL_T
      IMPLICIT NONE
      INTEGER(KINT):: ierr,i,j,dstart,dend
      REAL(KREAL)::dx,dy,dz,xs,ys,zs,dist
!
      dstart = 0; dend = 8 
        
! ----------      change made by Rafa          -----------------------------
      open (2, file=trim(filename)//".ts",access="append",iostat=ierr)
! ----------      End of change made by Rafa   ----------------------------  

      if (ierr.eq.0) then
          do 13 i=1,ncen
             do 14 j=i+1,ncen
             dx = rcen(1,i) - rcen(1,j) 
             dy = rcen(2,i) - rcen(2,j) 
             dz = rcen(3,i) - rcen(3,j) 
             dist = dsqrt(dx**2 + dy**2 + dz**2)
             if(dist.gt.dstart.and.dist.lt.dend) then
                xs = (rcen(1,i) +  rcen(1,j))/2   
                ys = (rcen(2,i) +  rcen(2,j))/2   
                zs = (rcen(3,i) +  rcen(3,j))/2  
                write(2,11)xs,ys,zs  
                xs = (3*rcen(1,i) +  rcen(1,j))/4   
                ys = (3*rcen(2,i) +  rcen(2,j))/4   
                zs = (3*rcen(3,i) +  rcen(3,j))/4   
                write(2,11)xs,ys,zs  
                xs = (3*rcen(1,j) +  rcen(1,i))/4   
                ys = (3*rcen(2,j) +  rcen(2,i))/4   
                zs = (3*rcen(3,j) +  rcen(3,i))/4   
                write(2,11)xs,ys,zs  
                xs = rcen(1,i) + (rcen(1,j) -  rcen(1,i))/5   
                ys = rcen(2,i) + (rcen(2,j) -  rcen(2,i))/5   
                zs = rcen(3,i) + (rcen(3,j) -  rcen(3,i))/5   
                write(2,11)xs,ys,zs  
                xs = rcen(1,i) + 2*(rcen(1,j) -  rcen(1,i))/5   
                ys = rcen(2,i) + 2*(rcen(2,j) -  rcen(2,i))/5   
                zs = rcen(3,i) + 2*(rcen(3,j) -  rcen(3,i))/5   
                write(2,11)xs,ys,zs  
                xs = rcen(1,i) + 3*(rcen(1,j) -  rcen(1,i))/5   
                ys = rcen(2,i) + 3*(rcen(2,j) -  rcen(2,i))/5   
                zs = rcen(3,i) + 3*(rcen(3,j) -  rcen(3,i))/5   
                write(2,11)xs,ys,zs  
                xs = rcen(1,i) + 4*(rcen(1,j) -  rcen(1,i))/5   
                ys = rcen(2,i) + 4*(rcen(2,j) -  rcen(2,i))/5   
                zs = rcen(3,i) + 4*(rcen(3,j) -  rcen(3,i))/5   
                write(2,11)xs,ys,zs  
             endif
14           continue  
13        continue  
      end if
      close(2)


11    format(3f9.4,2X)
      return 
      END SUBROUTINE

!============================================================      

!> @details Used mainly for creating guess points for MESP topography.
      SUBROUTINE AROUNDNUCLEUS 
!*****************************************     
      USE DAMINITIAL_T
      IMPLICIT NONE 
      INTEGER(KINT)::ierr,i,j
      REAL(KREAL):: dx,dy,dz,dist,h
      REAL(KREAL), ALLOCATABLE:: xg(:),yg(:),zg(:)
!        
      h=3.7
      open(32,status="SCRATCH")
      do 20 i=1,ncen
            write(32,11)rcen(1,i)+h,rcen(2,i),rcen(3,i)
            write(32,11)rcen(1,i)-h,rcen(2,i),rcen(3,i)
            write(32,11)rcen(1,i),rcen(2,i)+h,rcen(3,i)
            write(32,11)rcen(1,i),rcen(2,i)-h,rcen(3,i)
            write(32,11)rcen(1,i),rcen(2,i),rcen(3,i)+h
            write(32,11)rcen(1,i),rcen(2,i),rcen(3,i)-h
            write(32,11)rcen(1,i)+h,rcen(2,i)+h,rcen(3,i)+h  
            write(32,11)rcen(1,i)+h,rcen(2,i)+h,rcen(3,i)-h
            write(32,11)rcen(1,i)+h,rcen(2,i)-h,rcen(3,i)+h
            write(32,11)rcen(1,i)+h,rcen(2,i)-h,rcen(3,i)-h
            write(32,11)rcen(1,i)-h,rcen(2,i)+h,rcen(3,i)+h
            write(32,11)rcen(1,i)-h,rcen(2,i)+h,rcen(3,i)-h
            write(32,11)rcen(1,i)-h,rcen(2,i)-h,rcen(3,i)+h 
            write(32,11)rcen(1,i)-h,rcen(2,i)-h,rcen(3,i)-h 
20    continue
      rewind(32)
      
      allocate(xg(NCEN*14),yg(NCEN*14),zg(NCEN*14))
        
      open(2,file = trim(filename)//".ts",access="append",iostat=ierr)
      do 30 i=1,ncen*14
          read(32,*)xg(i),yg(i),zg(i)
             do 40 j=1,ncen  
                dx = xg(i) - rcen(1,j)
                dy = yg(i) - rcen(2,j)
                dz = zg(i) - rcen(3,j)
                dist = dsqrt(dx**2 + dy**2 + dz**2)
                if(dist.lt.1.0) goto 30
40           continue
          write(2,11)xg(i),yg(i),zg(i)  
30    continue
      deallocate (xg,yg,zg)
!
      close(2)
      close(32)  
!
11    format(3f9.4,2X)
      return
      end subroutine
!======================================
      
!*************************************************        
!SUBROUTINE FOR SORTING OUT THE IDENTICAL GUESS POINTS 
!> @details Guess points generated by different methods are sorted out based on
!! distance criterion.
        SUBROUTINE SORTGCPS  
        USE DAMINITIAL_T, ONLY: filename,KINT,KREAL,KINT8,KREAL4
        IMPLICIT NONE 
        CHARACTER (len=50):: dummy
        INTEGER(KINT)::icps,i,j,ierr
        REAL(KREAL):: dxi,dyi,dzi,dist
        REAL(KREAL), ALLOCATABLE::  xi(:),yi(:),zi(:)
        
        open(1,file=trim(filename)//".ts",status="old",iostat=ierr)
        
        icps=0
        do
            read(1,*,end=123,err=123)dummy 
            icps = icps + 1
        enddo
123     rewind(1)
        allocate (xi(icps),yi(icps),zi(icps))
        do i = 1,icps
            read(1,*,end=223,err=223)xi(i),yi(i),zi(i) 
        enddo
223     close(1)
        
        open(3,file=trim(filename)//".gcp",status="replace")

        do 312 i=1,icps
           do 313 j=1,i-1
              dxi = xi(i)-xi(j)
              dyi = yi(i)-yi(j)
              dzi = zi(i)-zi(j)
              dist= sqrt(dxi**2+dyi**2+dzi**2)
              if (dist.lt.0.05) go to 312
313        continue
           write(3,"(3(2x,f9.5))") xi(i),yi(i),zi(i)
312     continue
        deallocate (xi,yi,zi)
        close(3)
        return 
        END                  

!> @details Guess points are subjected to optimization using lbfgs optimization subroutine.             
      SUBROUTINE OPTDRIVER 
!**********************************************      
      USE DAMINITIAL_T, ONLY: ISWINDOWS, CNVG, BOXL, PROJECTNAME,NCEN, RCEN, IOUT, &
                              MED,MESP,FILENAME, filenamewin, KINT,KREAL,KREAL4,KINT8	! Modified by Rafa
      USE GRADPATH_T, ONLY: EVAL,EVEC,AUANG,ANGAU
      IMPLICIT NONE
      INTEGER,  PARAMETER           :: N = 3, M = 6, IPRINT = -2, ICMAX=1000 
      INTEGER,  PARAMETER           :: DP = KIND(1.0D0)
      REAL(DP), PARAMETER           :: FACTR  = 1.0D+3, PGTOL  = 1.0D-18
      INTEGER,  PARAMETER           :: NP=5   
      CHARACTER(LEN=60)             :: TASK, CSAVE
      LOGICAL                       :: LSAVE(4), lexist ! Modified by Rafa introducing lexist
      INTEGER                       :: ISAVE(44)
      REAL(DP)                      :: F
      REAL(DP)                      :: DSAVE(29)
      INTEGER, DIMENSION(21)        :: NPTSR
      INTEGER,  ALLOCATABLE         :: NBD(:), IWA(:)
      REAL(DP), ALLOCATABLE         :: X(:), L(:), U(:), G(:), WA(:)
      REAL(KREAL), ALLOCATABLE      :: xmpi(:),ympi(:),zmpi(:)
      REAL(KREAL), ALLOCATABLE      :: xscp(:),yscp(:),zscp(:)
      REAL(KREAL), ALLOCATABLE      :: VESP(:),vscp(:)
      CHARACTER(1), ALLOCATABLE     :: symb(:)
      CHARACTER                     :: symbc
      INTEGER(KINT), ALLOCATABLE    :: ntype(:) 
      REAL(KREAL)                   :: dxi,dyi,dzi,dist,dxin,dyin,dzin,distn
      REAL(KREAL), DIMENSION(NP)    :: D  
      REAL(KREAL), DIMENSION(NP,NP) :: V 
      CHARACTER(1), DIMENSION(4)    :: sym
      INTEGER(KINT), DIMENSION(4)   :: ntyp,nstrt
      REAL(DP)                      :: T1, T2
      INTEGER(KINT)                 :: I, IERR, J, K, KNT, ICOUNT, IE, JE, IS, IFILP	! Modified by Rafa introducing and KNT
      INTEGER(KINT)                 :: NPTS, IOS, ITR, NTCP, NIND
      INTEGER(KINT)                 :: IRANK,IRNK!,NGCP
      CHARACTER(4)                  :: IRANKD
      CHARACTER(50)                 :: DUMMY
      REAL(KREAL)                   :: X1,Y1,Z1
      REAL(KREAL)                   :: den,vtot,drvxtot,drvytot,drvztot,dxxtot
      REAL(KREAL)                   :: dxytot,dxztot,dyytot,dyztot,dzztot
      
      allocate ( nbd(n), x(n), l(n), u(n), g(n) )
      allocate ( iwa(3*n) )
      allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
        
      open (unit = 15,file=trim(filename)//".tcp",access="append",status="replace")
 
      ! All the guesscps which converge according to provided criteria are put in filename.tcp file.
      itr=0 
      do !itr= 1,5 ! Introducing iteration. All the guesscps not
                  ! converged in first go are placed in filename.gcp_?? file 
                  ! for reoptimization process. 
        itr = itr + 1 
      
        open(unit = 1, file = trim(filename)//".gcp",status='old')    ! File contains already created  guess points
 
        npts = 0  
        do 
            read(1,"(a50)",end=104) dummy
            npts = npts + 1  
        end do
104     continue
        ! NPTS contains number of guess points
        rewind(1)   
        nptsr(itr)=npts 
        
!         if (itr.gt.1.and.nptsr(itr).eq.nptsr(itr-1).or.itr.gt.20) then
!            close(1)
!            exit
!         end if  
        if (itr.gt.1) then
		if(nptsr(itr).eq.nptsr(itr-1) .or. itr.gt.20) then
           close(1)
           exit
          end if 
        end if 
        
        allocate (xmpi(npts),ympi(npts),zmpi(npts),stat=ierr)
        do i= 1,npts
            read(1,*,err=105)xmpi(i),ympi(i),zmpi(i)
        end do
105     continue
        close(1)

        if (i.lt.npts) write(iout,"(a,2x,I5)")"Error in guess point at", I                

        open(unit = 16, file = trim(filename)//".gcp-t")
        do i = 1,npts

            icount = 0 ; f=0 ; g(1)=0 ; g(2)=0 ; g(3)=0
            vtot=0 ; den=0
            drvxtot=0 ; drvytot=0 ; drvztot=0 ; dxxtot=0 
            dxytot=0 ; dxztot=0 ; dyytot=0 ; dyztot=0 ; dzztot=0 
            x(1) = xmpi(i)
            x(2) = ympi(i)
            x(3) = zmpi(i)
             
            do 10 j=1, n
                nbd(j) = 2 
                l(j)   = 0
                u(j)   = 0
                l(j)   = x(j)-boxl
                u(j)   = x(j)+boxl
10          continue
            task = 'START'
            do while   (task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START')
!               This is the call to the L-BFGS-B code.
                call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,& 
                    & csave,lsave,isave,dsave)
                if (task(1:2) .eq. 'FG') then 
!                   the minimization routine has returned to request the
!                   function f and gradient g values at the current x.
!                   Compute function value f for the sample problem.
!                   and gradient g for the sample problem.
                    x1 = x(1) ; y1 = x(2) ; z1 = x(3) 
    
                    if(med) then
                        CALL DAMDEN(vtot,drvxtot,drvytot,drvztot,dxxtot, &
                            &   dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
!                       write(6,420)x1,y1,z1,vtot,drvxtot,drvytot,drvztot
                    elseif(mesp) then
                        CALL DAMPOT(vtot,drvxtot,drvytot,drvztot,dxxtot, &
                            &   dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
!                       write(6,420)x1,y1,z1,vtot,drvxtot,drvytot,drvztot
                    end if   
                    f = drvxtot**2 + drvytot**2 + drvztot**2 
!       
                    g(1) = 2*(drvxtot*dxxtot + drvytot*dxytot + drvztot*dxztot) 
                    g(2) = 2*(drvxtot*dxytot + drvytot*dyytot + drvztot*dyztot) 
                    g(3) = 2*(drvxtot*dxztot + drvytot*dyztot + drvztot*dzztot) 
                    icount=icount+1
                    ! Check for maximum number of steps for the
                    ! convergence of a point.
                    if (icount.gt.icmax) then
                        write(iout,*)"ICOUNT EXCEEDED THE MAXIMUM LIMIT"
                        write(16,"(3(f9.4,2x))") (x(k),k=1,3)
                        exit
                    end if    
                elseif (task(1:5).eq.'NEW_X')  then  
                    if (dsave(2).le.cnvg) then
                        write(15,"(4(ES14.7,2x))") (x(k),k=1,3),vtot
                        exit
                    end if    
                elseif (task(1:4).eq.'CONV')  then  
                    if (dsave(2).gt.cnvg) then
                        write(16,"(3(f9.4,2x))") (x(k),k=1,3)
                    else
                        write(15,"(4(ES14.7,2x))") (x(k),k=1,3),vtot
                    endif
                elseif (task(1:4).eq.'ABNO')  then  
                    if (dsave(2).gt.cnvg) then
                        write(16,"(3(f9.4,2x))") (x(k),k=1,3)
                    end if
                elseif (task(1:4).eq.'STOP') then
                    if (dsave(2).gt.cnvg) then
                        write(16,"(3(f9.4,2x))") (x(k),k=1,3)
                    endif
                elseif (task(1:5).eq.'ERROR')  then  
                    call error(1,"error called from lbfgs subroutine")
                endif
            enddo
        enddo  ! End to loop of guess points
        continue        
        if (allocated(xmpi)) deallocate(xmpi)
        if (allocated(ympi)) deallocate(ympi)
        if (allocated(zmpi)) deallocate(zmpi)
        close(16) ! File contains those points which did not converge 
            ! Putting all the unconverged guess points in new file .gcp 
! ----------      change made by Rafa          -----------------------------      
! 		CALL SYSTEM("rm -f "//trim(filename)//".gcp")
! 	    CALL SYSTEM ( "cat "//trim(filename)//".gcp-t > "//trim(filename)//".gcp")
! 	    CALL SYSTEM("rm -f "//trim(filename)//".gcp-t")
       if (iswindows) then
	    CALL SYSTEM ( "move /y """//trim(filenamewin)//".gcp-t""  """//trim(filenamewin)//".gcp""")
       else
! 	    CALL SYSTEM ( "mv "//trim(filename)//".gcp-t  "//trim(filename)//".gcp")
	    CALL SYSTEM("rm -f "//trim(filename)//".gcp")
	    CALL SYSTEM ( "cat "//trim(filename)//".gcp-t > "//trim(filename)//".gcp")
	    CALL SYSTEM("rm -f "//trim(filename)//".gcp-t")
	  endif
! ----------      End of change made by Rafa   ----------------------------  
      end do  
      close(15) ! File contains all the converged guess points
      
!=============================================================      
           
       
      ! Sorting the guess points to list only the distinct CPs
      ! CPs which are 0.05 au away from each other are considered to be
      ! different  
        
	  open (17,file = trim(filename)//".tcp",status="old")

      npts = 0     ! NPTS now stores total number of unsorted CPs 
      do 
          read(17,"(a50)",end=106) dummy 
          npts = npts + 1  
      end do
106   continue 
      rewind(17) 
      allocate (xmpi(npts),ympi(npts),zmpi(npts),vesp(npts),stat=ierr)
      ! Defined new array for storing the sorted CPs
      ! Currently the array length is equal to total number of unsorted
      ! CPs 
      allocate (xscp(npts),yscp(npts),zscp(npts),vscp(npts),stat=ierr) 
      do i = 1, npts
          read(17,*)xmpi(i),ympi(i),zmpi(i),vesp(i)
      end do
      close(17)
        
!       CALL SYSTEM ("rm -f "//trim(filename)//".tcp")
        if (iswindows) then
            lexist = .false.	! Checks whether the file .tcp exists
            inquire(file=TRIM(filenamewin)//".tcp", exist=lexist, iostat=ierr)
            if (lexist) call system("del /f /q  """//TRIM(filenamewin)//".tcp""")
        else
	    CALL SYSTEM ( "rm -f "//trim(filename)//trim(filename)//".tcp")
        endif
       
      ntcp = 0      ! Stores the number of sorted CPs
      do i = 1, npts
        if (med.and.vesp(i).lt.0.0001) go to 107 
            ! checking if the CP is the nucleus itself
            do k=1,ncen
               dxin = (xmpi(i)-rcen(1,k))    
               dyin = (ympi(i)-rcen(2,k))    
               dzin = (zmpi(i)-rcen(3,k))    
               distn = sqrt(dxin**2+dyin**2+dzin**2)
               if (distn.lt.0.1) go to 107
            end do
                        
            do j = 1, i-1
                dxi= (xmpi(i)-xmpi(j))
                dyi= (ympi(i)-ympi(j)) 
                dzi= (zmpi(i)-zmpi(j))
                dist= sqrt(dxi**2+dyi**2+dzi**2)
                if (dist.lt.0.1) go to 107
            end do
            ntcp=ntcp+1
            xscp(ntcp)=xmpi(i) 
            yscp(ntcp)=ympi(i) 
            zscp(ntcp)=zmpi(i) 
            vscp(ntcp)=vesp(i)
107         continue
      end do
      
      ! Performing Hessian calculation for determining the nature of
      ! each CPs and assinging symbols accordingly.

      ntyp=0    
      
      open(18,file=trim(filename)//".tfcp")

      do i = 1, ntcp 
        if(med) then
            CALL DAMDEN(vscp(i),drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xscp(i),yscp(i),zscp(i))
        elseif(mesp) then
            CALL DAMPOT(vscp(i),drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xscp(i),yscp(i),zscp(i))
        endif
        CALL HESSIAN(IRANK,IRNK,D,V,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot)  
        if(irnk.eq.3) then
            if(irank.eq.3)  then
                ntyp(1)=ntyp(1)+1
                write(18,216)"x",xscp(i),yscp(i),zscp(i),vscp(i)
            else if(irank.eq.1) then
                ntyp(2)=ntyp(2)+1
                write(18,216)"y",xscp(i),yscp(i),zscp(i),vscp(i)
            else if(irank.eq.-1) then
                ntyp(3)=ntyp(3)+1
                write(18,216)"z",xscp(i),yscp(i),zscp(i),vscp(i)
            else if(irank.eq.-3) then
                ntyp(4)=ntyp(4)+1
                write(18,216)"m",xscp(i),yscp(i),zscp(i),vscp(i)
            end if    
            write(18,217)(D(J),J=1,3)
            write(18,217)((V(J,K),K=1,3),J=1,3)
        end if  
      end do
216   format(a,4(3x,e15.7))
217   format(3(e15.7,3X))
      close(18)   ! File contains sorted CPs whose nature have been determined
       
      deallocate (xmpi,ympi,zmpi,vesp)
      deallocate (xscp,yscp,zscp,vscp) 
      
      ! Sorting these CPs such that (3,+3) CPs i.e. "x" CPs appear
      ! first. Additionally the one with higher potential values appear
      ! on top of their respective category
      
      
      if (iswindows) then
		if (med) then
                    lexist = .false.	! Checks whether the file -cps-d.xyz exists
                    inquire(file=TRIM(filenamewin)//"-cps-d.xyz", exist=lexist, iostat=ierr)
                    if (lexist) call system("del /f /q  """//TRIM(filenamewin)//"-cps-d.xyz""")
                endif
		if (mesp) then
                    lexist = .false.	! Checks whether the file -cps-v.xyz exists
                    inquire(file=TRIM(filenamewin)//"-cps-v.xyz", exist=lexist, iostat=ierr)
                    if (lexist) call system("del /f /q  """//TRIM(filenamewin)//"-cps-v.xyz""")
                endif
		if (med) then
                    lexist = .false.	! Checks whether the file -cps-d.eigv exists
                    inquire(file=TRIM(filenamewin)//"-cps-d.eigv", exist=lexist, iostat=ierr)
                    if (lexist) call system("del /f /q  """//TRIM(filenamewin)//"-cps-d.eigv""")
                endif
		if (mesp) then
                    lexist = .false.	! Checks whether the file -cps-v.eigv exists
                    inquire(file=TRIM(filenamewin)//"-cps-v.eigv", exist=lexist, iostat=ierr)
                    if (lexist) call system("del /f /q  """//TRIM(filenamewin)//"-cps-v.eigv""")
                endif
       else
		if (med) call system("rm -f "//TRIM(filename)//"-cps-d.xyz")
		if (mesp) call system("rm -f "//TRIM(filename)//"-cps-v.xyz")
		if (med) call system("rm -f "//TRIM(filename)//"-cps-d.eigv")
		if (mesp) call system("rm -f "//TRIM(filename)//"-cps-v.eigv")
	  endif
	    
      allocate (symb(ntcp),xscp(ntcp),yscp(ntcp),zscp(ntcp),vscp(ntcp),stat=ierr) 
      allocate (eval(ntcp,3),evec(ntcp,3,3),stat=ierr)
      
      open (19,file=TRIM(filename)//".tfcp",iostat=ios)
        
      if (ios.ne.0) call error(1,"Error reading final cps")
      knt = 0	! Counter introduced by Rafa to cover the possibility that some CPs have zero eigenvalues and have been discarded
      do i = 1, ntcp
          read(19,*,err=109,end=109)symb(i),xscp(i),yscp(i),zscp(i),vscp(i)
          knt = knt + 1
          read(19,*,err=109)(eval(I,J),J=1,3)
          read(19,*,err=109)((evec(I,J,K),K=1,3),J=1,3)
      end do
109   continue
      close(19)
      ntcp = knt	! Introduced by Rafa to cover the possibility that some CPs have zero eigenvalues and have been discarded
      sym(1)="x"
      sym(2)="y"
      sym(3)="z"
      sym(4)="m"
           ! Sorting on the basis of symbol and listing in descending
           ! order of MESP value
      nstrt(1)=1
      nstrt(2)=ntyp(1)+1
      nstrt(3)=ntyp(1)+ntyp(2)+1   
      nstrt(4)=ntyp(1)+ntyp(2)+ntyp(3)+1
      do k=1,4
          if (ntyp(k).ne.0) then
              do i = nstrt(k),ntcp-1
                  do j = i+1, ntcp
                      if (symb(i).ne.sym(k).and.symb(j).eq.sym(k)) then
                          symbc=symb(i)
                          symb(i)=symb(j)
                          symb(j)=symbc
                          xscp(i)= xscp(i)+xscp(j)
                          xscp(j)= xscp(i)-xscp(j)
                          xscp(i)= xscp(i)-xscp(j)
                          yscp(i)= yscp(i)+yscp(j)
                          yscp(j)= yscp(i)-yscp(j)
                          yscp(i)= yscp(i)-yscp(j)
                          zscp(i)= zscp(i)+zscp(j)
                          zscp(j)= zscp(i)-zscp(j)
                          zscp(i)= zscp(i)-zscp(j)
                          vscp(i)= vscp(i)+vscp(j)
                          vscp(j)= vscp(i)-vscp(j)
                          vscp(i)= vscp(i)-vscp(j)
                          do ie = 1,3
                              Eval(i,ie) = Eval(i,ie) + Eval(j,ie)
                              Eval(j,ie) = Eval(i,ie) - Eval(j,ie)
                              Eval(i,ie) = Eval(i,ie) - Eval(j,ie)
                          end do
                          do ie =1,3
                              do je=1,3
                                  Evec(i,ie,je) = Evec(i,ie,je) + Evec(j,ie,je)
                                  Evec(j,ie,je) = Evec(i,ie,je) - Evec(j,ie,je)
                                  Evec(i,ie,je) = Evec(i,ie,je) - Evec(j,ie,je)
                              end do
                          end do        
                      end if
                  end do
              end do
          end if    
      end do
      do k =1,4
          if (ntyp(k).gt.1) then
              do i = nstrt(k),nstrt(k)+ntyp(k)-2
                  do j = i+1, nstrt(k)+ntyp(k)-1
                      if (vscp(j).gt.vscp(i)) then
                          symbc=symb(i)
                          symb(i)=symb(j)
                          symb(j)=symbc
                          xscp(i)= xscp(i)+xscp(j)
                          xscp(j)= xscp(i)-xscp(j)
                          xscp(i)= xscp(i)-xscp(j)
                          yscp(i)= yscp(i)+yscp(j)
                          yscp(j)= yscp(i)-yscp(j)
                          yscp(i)= yscp(i)-yscp(j)
                          zscp(i)= zscp(i)+zscp(j)
                          zscp(j)= zscp(i)-zscp(j)
                          zscp(i)= zscp(i)-zscp(j)
                          vscp(i)= vscp(i)+vscp(j)
                          vscp(j)= vscp(i)-vscp(j)
                          vscp(i)= vscp(i)-vscp(j)
                          do ie = 1,3
                              Eval(i,ie) = Eval(i,ie) + Eval(j,ie)
                              Eval(j,ie) = Eval(i,ie) - Eval(j,ie)
                              Eval(i,ie) = Eval(i,ie) - Eval(j,ie)
                          end do
                          do ie =1,3
                              do je=1,3
                                  Evec(i,ie,je) = Evec(i,ie,je) + Evec(j,ie,je)
                                  Evec(j,ie,je) = Evec(i,ie,je) - Evec(j,ie,je)
                                  Evec(i,ie,je) = Evec(i,ie,je) - Evec(j,ie,je)
                              end do
                          end do        
                      end if
                  end do
              end do
          end if
      end do
! 		Recomputes ntcp to prevent accounting more than once repeated CPs (modification introduced by Rafa)
		ntcp = 0
		do k =1,4
			ntcp = ntcp + ntyp(k)
		enddo        
      write(iout,*)
      write(iout,*)
      if (med) then
            write(iout,"(a)")"Critical Points and Scalar field values of Electron density:"
      else
            write(iout,"(a)")"Critical Points and Scalar field values of Electrostatic potential:"
      endif

      write(iout,"(a)")"-------------------------------------------------------------------------------" 
      write(iout,"(t1,a,t20,a,t38,a,t56,a,t72,a)")"Type","X","Y","Z", "Fval"
      write(iout,"(a)")"-------------------------------------------------------------------------------" 

      if (med) open(20, file=trim(filename)//"-cps-d.xyz")
      if (mesp) open(20, file=trim(filename)//"-cps-v.xyz")
      if (med) open(55, file=trim(filename)//"-cps-d.eigv")
      if (mesp) open(55, file=trim(filename)//"-cps-v.eigv")
      nind=0
      write(20,"(1x,I4)") ntcp
      write(20,*)
      write(55,"(1x,I4)") ntcp
      
      do i = 1,4
          do j = 1, ntyp(i) 
               nind = nind+1
               write(iout,"(A,1x,I4,4(3x,e15.7))") symb(nind),nind,xscp(nind),yscp(nind),zscp(nind),vscp(nind)
               write(20,"(A,4(3x,e15.7))") symb(nind),xscp(nind)*auang,yscp(nind)*auang,zscp(nind)*auang,vscp(nind)
               ifilp=55
               call draw_eigvec(ifilp,xscp(nind),yscp(nind),zscp(nind),nind)
            end do
      end do

        !In case there exists nonnuclear maxima(only in MED CPs), the
        !total number of (3,-3) CPs will be sum of non nuclear maxima
        !and nuclear centers.
        
      close(20) 
      close(55) 
        
      write(iout,"(a)")"-------------------------------------------------------------------------------" 
      write(iout,"(a)")"-------------------------------------------------------------------------------" 
      write(iout,*)
      do k=1,3
          write(iout,"(2(a,1x),I4)")sym(k),"=",ntyp(k)
      end do
      write(iout,"(a,1x,I4)") "m =", ncen+ntyp(4)
      write(iout,"(a,1x,I4)")"Calculated Euler Characteristic =",ntyp(1)-ntyp(2)+ntyp(3)-(ncen+ntyp(4))
      write(iout,*)
      write(iout,*)

      nind=0 
      do i = 1,4
          do j = 1, ntyp(i) 
               nind = nind+1
               write(iout,"(A,I5,t20,4(e15.7,3x))") symb(nind),nind,xscp(nind),yscp(nind),zscp(nind),vscp(nind)
               write(iout,*)""
               write(iout,289)"Eigenvalues:","Eval1","Eval2","Eval3"
               write(iout,"(t20,3(e15.7,3x))") (Eval(NIND,IE),IE=1,3)
               write(iout,*)""
               write(iout,"(a,t28,a,t46, a, t64, a)")   "Eigenvectors:","Evec1","Evec2","Evec3"
               do ie=1,3
                    write(iout,"(t20,3(e15.7,3x))") (Evec(NIND,IE,JE),JE=1,3)
               end do
               write(iout,"(A)") ""
               write(iout,"(A)") ""
          end do
      end do
289   format(t1,a,t28,a,t46,a,t64,a)
      deallocate (symb,xscp,yscp,zscp,vscp) 
      deallocate (eval,evec)
      if (iswindows) then
            lexist = .false.	! Checks whether the file .tfcp exists
            inquire(file=TRIM(filenamewin)//".tfcp", exist=lexist, iostat=ierr)
            if (lexist) call system("del /f /q  """//TRIM(filenamewin)//".tfcp""")
      else
	  call system("rm -f "//TRIM(filename)//".tfcp")
	 endif
      return 
      end

!----------------------------------------------------------
!> @details Arranges the eigenvalues and corresponding eigenvectors in
!! decreasing order of magnitude of eigenvalues.
!----------------------------------------------------------      
      SUBROUTINE DRAW_EIGVEC(IFILP,XOCP,YOCP,ZOCP,NIND)
      USE GRADPATH_T, ONLY: EVAL,EVEC,AUANG,ANGAU
      IMPLICIT NONE 
      REAL*8 :: EIVX,EIVY,EIVZ
      REAL*8 :: XOCP,YOCP,ZOCP
      INTEGER :: IQ, JQ, KQ, IFILP,NIND
      
      do iq = 1,2 
         do jq = iq+1, 3
            if (abs(eval(nind,jq))>abs(eval(nind,iq))) then
               eval(nind,iq) = eval(nind,iq) + eval(nind,jq)
               eval(nind,jq) = eval(nind,iq) - eval(nind,jq)      
               eval(nind,iq) = eval(nind,iq) - eval(nind,jq) 
               do kq = 1,3   
                  evec(nind,kq,iq) = evec(nind,kq,iq) + evec(nind,kq,jq)
                  evec(nind,kq,jq) = evec(nind,kq,iq) - evec(nind,kq,jq)  
                  evec(nind,kq,iq) = evec(nind,kq,iq) - evec(nind,kq,jq)  
               end do
            end if
         end do
      end do

      eivx = xocp*auang
      eivy = yocp*auang
      eivz = zocp*auang
      if (eval(nind,1).lt.0.0) WRITE(ifilp,9)eivx, eivy, eivz,evec(nind,1,1),evec(nind,2,1),evec(nind,3,1),"-1"
      if (eval(nind,1).gt.0.0) WRITE(ifilp,9)eivx, eivy, eivz,evec(nind,1,1),evec(nind,2,1),evec(nind,3,1)," 1"
      if (eval(nind,2).lt.0.0) WRITE(ifilp,9)eivx, eivy, eivz,evec(nind,1,2),evec(nind,2,2),evec(nind,3,2),"-1"
      if (eval(nind,2).gt.0.0) WRITE(ifilp,9)eivx, eivy, eivz,evec(nind,1,2),evec(nind,2,2),evec(nind,3,2)," 1"
      if (eval(nind,3).lt.0.0) WRITE(ifilp,9)eivx, eivy, eivz,evec(nind,1,3),evec(nind,2,3),evec(nind,3,3),"-1"
      if (eval(nind,3).gt.0.0) WRITE(ifilp,9)eivx, eivy, eivz,evec(nind,1,3),evec(nind,2,3),evec(nind,3,3)," 1"
9     format(6(f6.2,1x),a)                  

      return
      end subroutine
      

!============================================================================================

!> @details Creates gradient path for MED and MESP topography
!****************************************************      
      SUBROUTINE GRADPATH_S
!****************************************************      
      USE GRADPATH_T 
      IMPLICIT none
      INTEGER(KINT)               :: irank, irnk
      INTEGER(KINT)               :: i, j ,k, l 
      INTEGER(KINT)               :: ierr
      INTEGER(KINT), PARAMETER    :: ifila = 190, irm = 101, irmp = 182
      REAL(KREAL)                 :: x1,y1,z1,vtot,drvxtot,drvytot,drvztot
      REAL(KREAL)                 :: dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot
      REAL(KREAL), DIMENSION(5)   :: D
      REAL(KREAL), DIMENSION(5,5) :: V
      CHARACTER                   :: IRANKD*4,SI*2
      INTEGER(KINT)               :: xcptyp,ycptyp,zcptyp,mcptyp 
      REAL(KREAL)                 :: evx, evy, evz,evx1,evy1,evz1
      logical                     :: lexist     ! Introduced by Rafa
      
      if(med) open(21,file=trim(filename)//"-cps-d.xyz",status="old")
      if(mesp) open(21,file=trim(filename)//"-cps-v.xyz",status="old")
       
      ncp=0 !Initialize Number of CPs
      read(21,*) ncp
      read(21,*)

      allocate (symcp(ncp),xcp(ncp),ycp(ncp),zcp(ncp),stat=ierr)
      
      do i=1,ncp
        read(21,*,err=121,end=121) symcp(i), xcp(i), ycp(i), zcp(i),vdummy
        xcp(i)=xcp(i)*angau
        ycp(i)=ycp(i)*angau
        zcp(i)=zcp(i)*angau
      end do
121   rewind(21)
      close(21)

      open ( irmp, file=trim(filename)//".rmdat",status="replace")
      ! Separate temporary file to store the gradient path information.
      
      !Initializing Variables
      drsq=0.0; x1=0.0; y1=0.0; z1=0.0

      !Generate a box around the molecule wherein the gradient path will be
      !generated---
      
      xdiff=0  !Xdiff,Ydiff,Zdiff stores coordinates of farthest entity
      ydiff=0  !either CP or Nucleus, from center of mass
      zdiff=0

      do i=1,ncen   ! Deciding farthest Nucleus from COM
         temp1=abs(rcen(1,i)-cox)
         if (temp1>xdiff) xdiff=temp1
         temp1=abs(rcen(2,i)-coy)
         if (temp1>ydiff) ydiff=temp1
         temp1=abs(rcen(3,i)-coz)
         if (temp1>zdiff) zdiff=temp1
      end do
      
      do i=1,ncp   ! Deciding farthest CP from center of mass
         temp1=abs(xcp(i)-cox)
         if (temp1>xdiff) xdiff=temp1
         temp1=abs(ycp(i)-coy)
         if (temp1>ydiff) ydiff=temp1
         temp1=abs(zcp(i)-coz)
         if (temp1>zdiff) zdiff=temp1
      end do
         
      ixpln = int(abs(cox-xdiff))
      iypln = int(abs(coy-ydiff))
      izpln = int(abs(coz-zdiff))
      ! Insures 3 dimensional boundary for gradient path generation.
      if (ixpln.eq.0) xdiff = xdiff + 1.0 
      if (iypln.eq.0) ydiff = ydiff + 1.0
      if (izpln.eq.0) zdiff = zdiff + 1.0

      xboxgs = cox - (boxg*xdiff) !x_boxgs s START
      xboxge = cox + (boxg*xdiff) !x_boxge e END
      yboxgs = coy - (boxg*ydiff)
      yboxge = coy + (boxg*ydiff)
      zboxgs = coz - (boxg*zdiff)
      zboxge = coz + (boxg*zdiff)
      ! The above coordinates of the box decides the extent to which
      ! gradient path should evolve, especially in case of MESP.
      ! In case of MED gradient path are stopped by checking the value
      ! of electron density if less than 0.001 au
       
      do cpnum=1,ncp 
         lderiv2=.true.
         if (med) then
         call damden(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,   &
                        & xcp(cpnum),ycp(cpnum),zcp(cpnum))
         elseif (mesp) then
         call dampot(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,   &
                        & xcp(cpnum),ycp(cpnum),zcp(cpnum))
         end if
         call hessian(irank,irnk,d,v,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot)  
         ! Hessian diagonalizes the Hessian matrix formed by dxx,dxy.. and find
         ! eigenvalues as well as eigenvectors
          !===================================================
          xnew(1)=xcp(cpnum)+(fdisp*(v(1,1))) !First eigenvector Direction 
          ynew(1)=ycp(cpnum)+(fdisp*(v(2,1))) 
          znew(1)=zcp(cpnum)+(fdisp*(v(3,1)))
          xnew(2)=xcp(cpnum)+(fdisp*(v(1,2))) !Second eigenvector Direction
          ynew(2)=ycp(cpnum)+(fdisp*(v(2,2)))
          znew(2)=zcp(cpnum)+(fdisp*(v(3,2)))
          xnew(3)=xcp(cpnum)+(fdisp*(v(1,3))) !Third eigenvector Direction
          ynew(3)=ycp(cpnum)+(fdisp*(v(2,3)))
          znew(3)=zcp(cpnum)+(fdisp*(v(3,3)))
          !---------------------------------------------------
          xnew(4)=xcp(cpnum)-(fdisp*(v(1,1))) !Negative of First eigenvector Direction
          ynew(4)=ycp(cpnum)-(fdisp*(v(2,1)))
          znew(4)=zcp(cpnum)-(fdisp*(v(3,1)))
          xnew(5)=xcp(cpnum)-(fdisp*(v(1,2))) !Negative of Second eigenvector Direction
          ynew(5)=ycp(cpnum)-(fdisp*(v(2,2)))
          znew(5)=zcp(cpnum)-(fdisp*(v(3,2)))
          xnew(6)=xcp(cpnum)-(fdisp*(v(1,3))) !Negative of third eigenvector Direction
          ynew(6)=ycp(cpnum)-(fdisp*(v(2,3)))
          znew(6)=zcp(cpnum)-(fdisp*(v(3,3)))
          !===================================================
            ! Arranges the eigenvalues and corresponding eigenvectors in
            ! decreasing order of magnitude of eigenvalues. Writes in
            ! the .eigv file in appended mode. The file contains the
            ! number of CPs on top followed by the three eigenvector
            ! directions of each CP. 
          lderiv2=.false.    
          do j= 1,2  ! when j=1, the gradient path is followed in the direction of the gradient  
          do i= 1,6  ! i=1,3 Moves in the direction of eigenvector, while i=4,6 moves in the opposite direction of the eigenvector 
               flag=0
               cnt=0  
               cnt1=0
               xpp=0
               ypp=0
               zpp=0
               xp=0
               yp=0
               zp=0
               x1=xnew(i)
               y1=ynew(i)
               z1=znew(i)
               write(irmp,8)x1,y1,z1 
               do while (flag==0)
                  if (med) then 
                  call damden(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
                  elseif (mesp) then
                  call dampot(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
                  end if
                  !rationalise
                  drsq=1.0/((drvxtot**2)+(drvytot**2)+(drvztot**2)) !for normalizing gradient  
                  if (cnt1==0) xpp=x1 !storing previous point in xpp
                  if (cnt1==0) ypp=y1 !storing previous point in xpp
                  if (cnt1==0) zpp=z1 !storing previous point in xpp
                  cnt=cnt+1
                  cnt1=cnt1+1
                  if (j==1)then
                      x1=x1+(step*drvxtot*sqrt(drsq))  !step size is 0.001 au
                      y1=y1+(step*drvytot*sqrt(drsq))
                      z1=z1+(step*drvztot*sqrt(drsq))
                  else if (j==2) then  
                      X1=X1-(step*DRVXTOT*sqrt(drsq))
                      Y1=Y1-(step*DRVYTOT*sqrt(drsq))
                      Z1=Z1-(step*DRVZTOT*sqrt(drsq))
                  end if    
                  if ((cnt1*step)>=0.05) then ! Print points only when you have moved 0.05 au from previous point
                     write(irmp,8)x1,y1,z1 
                     cnt1=0
                  end if
                   
                  ! xpp,ypp,zpp,xp,yp and zp are also transferred
                  ! to terminate subroutine.
                  call terminate(vtot, drvxtot, drvytot, drvztot, x1, y1, z1, irmp)
               end do

          end do  
          end do

      end do ! End of loop over number of CPs, from where the gradient path originate.
       
8     format(3(f6.2,1x))                  
      
      close(irmp)  ! Closing the file containing the points on the gradient path, each at a distance of 0.05 au.
      
      
      open(102,file=trim(filename)//".rmdat",status="old")
        ! Alternatively, .rmdat file could be appended with points on
        ! gradient path. It could initially contain all the eigenvector
        ! directions of all the CPs.
! ----------      End of change made by Rafa   ---------------------------- 

      allocate (symcp1(ncp*12),ncp1(ncp*12),symcp2(ncp*12),stat=ierr)
      allocate (ncp2(ncp*12),ngp(ncp*12),nrma(0:ncp*12),ind(ncp*12),stat=ierr)
        ! Following part of the code reads the .rmdat, which at the end
        ! of each gradient path connecting two CPs has the information
        ! of the indices of the two CPs connected by gradient path.
        ! Shortest of the gradient path is chosen in case more than one
        ! gradient path is connecting same two CPs.  
      tnlrm=0   ! Stores total number of lines in .rmdat file
      nrm=0     ! Stores the number of lines which contain the
                  ! information of connection of CP (wriiten in the file from
                  ! terminate subroutine.)
      do 
          read(102,"(a10)",end=108)drm    ! Reads line as string
          tnlrm = tnlrm+1
          if ((index(drm,"x")).ne.0.or.(index(drm,"y")).ne.0.or.(index(drm,"z")).ne.0.or.(index(drm,"m")).ne.0) then
              nrm= nrm+1           
              nrma(nrm)=tnlrm   ! nrma array stores the lines where CP connection information is provided
              backspace(102)
              read(102,*)symcp1(nrm),ncp1(nrm),symcp2(nrm),ncp2(nrm),ngp(nrm)
          end if    
      end do
108   rewind(102)
      
      do i = 1,nrm
          ind(i)=1   ! An indicator ind of the size of nrm is set to one. 
      end do  


      do i = 1,nrm-1
          do j= i+1,nrm
              if (symcp1(i).eq.symcp1(j).and.ncp1(i).ne.ncp1(j)) exit 
              if (symcp1(i).ne.symcp1(j)) exit 
              if (symcp1(i).eq.symcp1(j).and.ncp1(i).eq.ncp1(j).and.symcp2(i).eq.symcp2(j).and.ncp2(i).eq.ncp2(j)) then
                  if (ngp(i).lt.ngp(j)) ind(j)=-1  
                  if (ngp(i).gt.ngp(j)) ind(i)=-1
                  ! Indicator ind is set to -1 when same two CPs are
                  ! connected by two or more gradient path. Longer one
                  ! is eliminated by setting the indicator to -1.    
              end if
              if (trim(symcp2(i)).eq."ob") ind(i) = -1        !Latest update by Anmol
              if (trim(symcp2(j)).eq."ob") ind(j) = -1        !Latest update by Anmol
          end do
      enddo              
      
      if (med) open (ifila, file = trim(filename)//"-d.gpdat", access="append",status="replace")
      if (mesp) open (ifila, file = trim(filename)//"-v.gpdat", access="append",status="replace")
      if (med.and.basin) open(12, file = trim(filename)//"-d.cnct", status="replace")
      if (mesp.and.basin) open(12, file = trim(filename)//"-v.cnct", status="replace")
        ! The above information from rmdat file about the shortest path
        ! is finally stored in gpdat file
        ! filename.fin stores the information of connectivity of a CP
        ! to any nucleus with its index. It helps in creating separate
        ! basin for separate atom.
        
!      if (basin) allocate(atsy(nrm),atin(nrm),cpsy(nrm),cpin(nrm),stat=ierr)
      j=1
      k=0
      nrma(0)=0
      do i=1,tnlrm
         read(102,"(a120)",end=119)drm
         if (i.ne.nrma(j).and.ind(j).eq.1) then
              write(ifila,"(a)")drm
         elseif (i.eq.nrma(j).and.ind(j).eq.1) then
              backspace(102)
              read(102,*)symcp1(j),ncp1(j),symcp2(j),ncp2(j),ngp(j)
              if (symcp2(j).ne."ob") then
                write(iout,"(2(a,1x,a,I3,1x))")"Gradient path connects",symcp1(j),ncp1(j),"to",symcp2(j),ncp2(j)
              end if 
              if (basin) then
                  if (trim(adjustl(symcp2(j))).eq."x".or.trim(adjustl(symcp2(j))).eq."y" &
                      & .or.trim(adjustl(symcp2(j))).eq."z".or.trim(adjustl(symcp2(j))).eq."m" &
                      & .or.trim(adjustl(symcp2(j))).eq."ob") then
                      continue
                  else    
                      write(12,"(2(a,2x,I4,2x))")symcp2(j),ncp2(j),symcp1(j),ncp1(j)
!                      k=k+1
!                      atsy(k)=symcp2(j)
!                      atin(k)=ncp2(j)
!                      cpsy(k)=symcp1(j)
!                      cpin(k)=ncp1(j)   
!                      nument=k
                  end if
              end if
              write(ifila,"(a)")""
              j= j+1
         end if     
         if (i.ne.nrma(j).and.ind(j).eq.-1) cycle
         if (i.eq.nrma(j).and.ind(j).eq.-1) j=j+1 
      END DO
119   close(102)  ! Closing rmdat file
      if (basin) close(12) ! Closing .fin file      
      if (iswindows) then
          lexist = .false.	! Checks whether the file .rmdat exists
            inquire(file=TRIM(filenamewin)//".rmdat", exist=lexist, iostat=ierr)
            if (lexist) call system("del /f /q  """//TRIM(filenamewin)//".rmdat""")
      else
            call system("rm -f "//trim(projectname)//".rmdat")
      endif
      deallocate (symcp1,ncp1,symcp2,ncp2,ngp,nrma,ind)
      deallocate(xcp, ycp,zcp,symcp)
      return
      end subroutine 
!---------------------------------------------------------------------------------

!> @details Decides the termination of gradient path based on gradient value,
!! distance from a CP and box size.
!*********************************************************************************     
      subroutine terminate(vtot, drvxtot, drvytot, drvztot, x1, y1,z1,irmp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE GRADPATH_T 
      IMPLICIT NONE 
      REAL(KREAL) :: vtot, drvxtot, drvytot, drvztot, x1, y1, z1
      INTEGER(KINT) :: I,J,K,L,IRMP

      if (med) then
        if (vtot.lt.0.0005) then  ! Takes the density value of 0.001 au to be the cutoff for the gradient path
            write(irmp,*) symcp(cpnum),cpnum,"ob",cpnum,cnt
            flag=1 
            return
        end if    
      elseif (mesp) then ! Uses the box size in case of mesp
        if ((x1>xboxge).or.(x1<xboxgs).or.(y1>yboxge).or.(y1<yboxgs).or.(z1>zboxge).or.(z1<zboxgs)) THEN
            write(irmp,*) symcp(cpnum),cpnum,"ob",cpnum,cnt
            flag=1 
            return
        end if
      end if

      if ((abs(drvxtot)<drcutcp).and.(abs(drvytot)<drcutcp).and.(abs(drvztot)<drcutcp)) then
         ! if all the derivatives are smaller than 10e-5
         
         do l=1, ncp ! To Check if the gradient path is in close proximity of a CP.
            if (l.ne.cpnum) then ! if current CP in the loop is not the
                                  ! same CP from which gradient path has originated then
               if ((abs(x1-xcp(l))<distcp).and.(abs(y1-ycp(l))<distcp).and.(abs(z1-zcp(l))<distcp)) then 
                  ! And the gradient path has reached within 0.05 au or less to another CP then
!                 write(iout,"(2(a,1x,a,1x,I3,1x),I7)") "Gradient path connects",symcp(cpnum),cpnum,"to",symcp(l),l, cnt   
                  write(irmp,*) symcp(cpnum),cpnum,symcp(l),l,cnt   
                  flag=1
                  return
               end if
            elseif (l.eq.cpnum) then
               if ((abs(x1-xcp(l))<distcp).and.(abs(y1-ycp(l))<distcp).and.(abs(z1-zcp(l))<distcp)) then 
                  write(irmp,*) symcp(cpnum),cpnum,"ob",cpnum,cnt
                  flag=1
                  return
               end if
            end if   
         end do
      end if   


      if ((abs(drvxtot)>drcutnuc).or.(abs(drvytot)>drcutnuc).or.(abs(drvztot)>drcutnuc)) then
         do l=1, ncen
            if ((abs(x1-rcen(1,l))<distnuc).and.(abs(y1-rcen(2,l))<distnuc).and.(abs(z1-rcen(3,l))<distnuc)) then 
               write(irmp,*) symcp(cpnum),cpnum,atmnam(l),l,cnt 
               flag=1
               return
            end if
         end do
      end if   
             
      !If the path finds a loop
              
      if((abs(xpp-x1).lt.1.0e-07).and.(abs(ypp-y1)<1.0e-07)  &
         &  .and.(abs(zpp-z1)<1.0e-07)) then
         write(irmp,*) symcp(cpnum),cpnum,"ob",cpnum,cnt
         flag=1
         return
      end if
      flag=0

      end subroutine

!******************************************************
!> @details Partitions the space of molecules in basin for MED and MESP.
      SUBROUTINE BASIN_S 
!****************************************************      
      USE GRADPATH_T 
      IMPLICIT NONE
      INTEGER(KINT)             :: IRANK, IRNK,RANDINP
      INTEGER(KINT)               :: I, J, K, L 
      INTEGER(KINT)               :: IERR,CSTEP,NZCP,NXYCP, IERRQT ! modified by Rafa 2019
      INTEGER(KINT), PARAMETER    :: IBILW = 186, IBILS = 187, IBILQT = 188 ! modified by Rafa 2019
      INTEGER(KINT)               :: iars, npar, noen, noenmax, ndots1, ndots1s, kl, ins, jns, kntvert, kntind ! modified by Rafa 2019
      INTEGER(KINT), DIMENSION(200):: ndotss
      INTEGER(KINT), DIMENSION(200):: ndots
      REAL(KREAL), PARAMETER     :: pival = 3.14159265359D0
      REAL(KREAL), PARAMETER     :: rad = pival/180.0d0
      REAL(KREAL)                :: phi
      REAL(KREAL)                :: X1,Y1,Z1,VTOT,DRVXTOT,DRVYTOT,DRVZTOT
      REAL(KREAL)                :: DXXTOT,DXYTOT,DXZTOT,DYYTOT,DYZTOT,DZZTOT
      REAL(KREAL)                :: xyz8(18) ! modified by Rafa 2019
      REAL(KREAL4)               :: x4, y4, z4, dx4, dy4, dz4 ! modified by Rafa 2019
      REAL(KREAL), DIMENSION(5)  :: D
      REAL(KREAL), DIMENSION(5,5):: V
      REAL(KREAL), DIMENSION(5)  :: VAL
      REAL(KREAL), DIMENSION(5,5):: VEC
      REAL(KREAL), DIMENSION(2,200):: xlines,ylines,zlines
      REAL(KREAL), DIMENSION(200)  :: xline1s,yline1s,zline1s
      REAL(KREAL), DIMENSION(2,200):: xline,yline,zline
      REAL(KREAL), DIMENSION(200)  :: xline1,yline1,zline1
      REAL(KREAL)                :: xter,yter,zter,lenth,lenths
      CHARACTER             :: IRANKD*4
      character(300) :: basinsname  ! Mofified by Rafa 2019
      LOGICAL               :: INTEGRATE
      REAL(KREAL), DIMENSION(3):: vnx,vny,vnz
      REAL(KREAL)           ::vecx,vecy,vecz,tang
      
      irank=0;irnk=0;cpnum=0;ixpln=0;iypln=0;izpln=0;tnlrm=0;nrm=0
      TEMP1=0.0;vdummy=0.0
      
      do i =1,5
          d(i)=0.0
      end do
      do i =1,5
          do j=1,5
              v(i,j)=0.0
          end do
      end do
                  
      if(med) open(25,file=trim(filename)//"-cps-d.xyz",status="old")
      if(mesp) open(25,file=trim(filename)//"-cps-v.xyz",status="old")

      ncp=0 !Initialize Number of CPs
      read(25,*) ncp
      read(25,*)

      allocate (symcp(ncp),xcp(ncp),ycp(ncp),zcp(ncp),stat=ierr)
      
      do i=1,ncp
        read(25,*,err=185,end=185) symcp(i), xcp(i), ycp(i), zcp(i),vdummy
        xcp(i)=xcp(i)*angau
        ycp(i)=ycp(i)*angau
        zcp(i)=zcp(i)*angau
      end do
185   rewind(25)
      close(25)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
      ierrqt = 1    ! if file basisname is not opened, prevents to load vertices and indices for DAMQT
      open (unit=IBILQT+1, file=trim(filename)//".vert",status="replace", iostat=ierrqt)
      if (ierrqt .ne. 0) then
          write(6,*) 'Error when opening file ', trim(filename)//".vert"
      endif
!        End of modification by Rafa 2019  ----------------------------------------------------------------------------------

      if (med.and.solidsurf) open(ibils,file=trim(filename)//"-d.bssol",status="replace")
      if (mesp.and.solidsurf) open(ibils,file=trim(filename)//"-v.bssol",status="replace")
      if (med.and.wireframe) open(ibilw,file=trim(filename)//"-d.bswir",status="replace")
      if (mesp.and.wireframe) open(ibilw,file=trim(filename)//"-v.bswir",status="replace")
      
      !Initializing Variables
      drsq=0.0; x1=0.0; y1=0.0; z1=0.0
      !Generate a box around the molecule wherein the gradient path will be
      !generated---
      
      xdiff=0  !Xdiff,Ydiff,Zdiff stores coordinates of farthest entity
      ydiff=0  !either CP or Nucleus, from center of mass
      zdiff=0

      do i=1,ncen   ! Deciding farthest Nucleus from COM
         temp1=abs(rcen(1,i)-cox)
         if (temp1>xdiff) xdiff=temp1
         temp1=abs(rcen(2,i)-coy)
         if (temp1>ydiff) ydiff=temp1
         temp1=abs(rcen(3,i)-coz)
         if (temp1>zdiff) zdiff=temp1
      end do
      
      do i=1,ncp   ! Deciding farthest CP from center of mass
         temp1=abs(xcp(i)-cox)
         if (temp1>xdiff) xdiff=temp1
         temp1=abs(ycp(i)-coy)
         if (temp1>ydiff) ydiff=temp1
         temp1=abs(zcp(i)-coz)
         if (temp1>zdiff) zdiff=temp1
      end do
         
      ixpln = int(abs(cox-xdiff))
      iypln = int(abs(coy-ydiff))
      izpln = int(abs(coz-zdiff))
      ! Insures 3 dimensional boundary for gradient path generation.
      if (ixpln.eq.0) xdiff = xdiff + 1.0 
      if (iypln.eq.0) ydiff = ydiff + 1.0
      if (izpln.eq.0) zdiff = zdiff + 1.0

      xboxgs = cox - (boxb*xdiff) !x_boxgs s START
      xboxge = cox + (boxb*xdiff) !x_boxge e END
      yboxgs = coy - (boxb*ydiff)
      yboxge = coy + (boxb*ydiff)
      zboxgs = coz - (boxb*zdiff)
      zboxge = coz + (boxb*zdiff)

      nzcp = 0
      nxycp = 0  

      ! Basin will be created usind the eigenvectors of (3,-1) CP ==> z CP only.
      do i = 1, ncp
         if (symcp(i).eq."z") nzcp = nzcp + 1 
!                Bug fixed by Rafa (September 2017)
!		if (symcp(i).ne."z") nxycp = nxycp + 1  Wrong test (m CP's come after z CP's and must not be accounted)
             if (symcp(i).ne."z" .and. symcp(i).ne."m") nxycp = nxycp + 1
      end do       
      kntvert = 0 ! modified by Rafa 2019
      do cpnum=nxycp+1, nxycp+nzcp 
         if (wireframe) write(ibilw,*)symcp(cpnum),cpnum
         if (solidsurf) write(ibils,*)symcp(cpnum),cpnum
         lderiv2=.true.
         if (med) then
         call damden(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,   &
                        & xcp(cpnum),ycp(cpnum),zcp(cpnum))
         elseif (mesp) then
         call dampot(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,   &
                        & xcp(cpnum),ycp(cpnum),zcp(cpnum))
         end if
         call hessian(irank,irnk,d,v,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot)  
         ! Hessian diagonalizes the Hessian matrix formed by dxx,dxy.. and find
         ! eigenvalues as well as eigenvectors
         do i = 1, 3
            val(i)= d(i)
         end do
          
         do i = 1, 3
            do j = 1, 3
                vec(i,j) = v(i,j)
            end do
         end do         
      
         do i = 1,2 
            do j = i+1, 3
!                if (abs(val(j))>abs(val(i))) then
                if (val(j)>val(i)) then
                    val(i)= val(i)+val(j)
                    val(j)= val(i)-val(j)      
                    val(i)= val(i)-val(j) 
                    do k = 1,3   
                        vec(k,i)=vec(k,i)+vec(k,j)
                        vec(k,j)=vec(k,i)-vec(k,j)  
                        vec(k,i)=vec(k,i)-vec(k,j)  
                    end do
                end if
            end do
         end do

         ! Depending on the angle defined by the user, total number of
         ! points to be generated around the plane defined by two
         ! eigenvectors is decided.   
         iars=360.0/angle
         npar=90.0/angle

         xdvd(1) = xcp(cpnum) + fdisp*vec(1,2)  
         ydvd(1) = ycp(cpnum) + fdisp*vec(2,2)
         zdvd(1) = zcp(cpnum) + fdisp*vec(3,2)
         xdvd(npar+1) = xcp(cpnum) + fdisp*vec(1,3)  
         ydvd(npar+1) = ycp(cpnum) + fdisp*vec(2,3)
         zdvd(npar+1) = zcp(cpnum) + fdisp*vec(3,3)
         xdvd(npar*2+1) = xcp(cpnum) - fdisp*vec(1,2)  
         ydvd(npar*2+1) = ycp(cpnum) - fdisp*vec(2,2)
         zdvd(npar*2+1) = zcp(cpnum) - fdisp*vec(3,2)
         xdvd(npar*3+1) = xcp(cpnum) - fdisp*vec(1,3)  
         ydvd(npar*3+1) = ycp(cpnum) - fdisp*vec(2,3)
         zdvd(npar*3+1) = zcp(cpnum) - fdisp*vec(3,3)
         
         ! The above array elements define the coordinates of points on
         ! the lying on the eigenvector, making a distance of fdisp away
         ! from CP.

         ! creating other points in the plane, based on the angle.
         phi = 0.0
         j= 0
         do i= 1,iars-(npar-1),npar
            phi=0.0
            j=i+npar
            if (i.eq.(iars-(npar-1))) j=1
            do k=i+1,i+(npar-1)
                phi=phi+angle 
                xdvd(k) = xcp(cpnum) + cos(phi*rad)*(xdvd(i)-xcp(cpnum)) + sin(phi*rad)*(xdvd(j)-xcp(cpnum))
                ydvd(k) = ycp(cpnum) + cos(phi*rad)*(ydvd(i)-ycp(cpnum)) + sin(phi*rad)*(ydvd(j)-ycp(cpnum))
                zdvd(k) = zcp(cpnum) + cos(phi*rad)*(zdvd(i)-zcp(cpnum)) + sin(phi*rad)*(zdvd(j)-zcp(cpnum))
            end do
         end do

          !===================================================
          lderiv2=.false.    
          xline=0.0
          yline=0.0
          zline=0.0
          ! since the two eigenvectors of z cp correspond to negative
          ! valued eigenvalue, (implies the gradient, points in the
          ! direction of CP through the corresponding eigenvectors.)
          ! only negative direction of gradient is used to create the basin.
          noenmax = 200   ! modified by Rafa 2019
          do i= 1,iars  ! The gradient path is originated in all directions of points created earlier
               flag=0
               cnt=0  
               cnt1=0
               xpp=0
               ypp=0
               zpp=0
               noen=0 
               if (i.eq.1) kl = 1 
               if (i.gt.1) kl = 2 
               X1=xdvd(i)
               Y1=ydvd(i)
               Z1=zdvd(i)
               do while (flag==0 .and. noen .lt. noenmax)   ! modified by Rafa 2019
                  if (med) then 
                  call damden(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
                  elseif (mesp) then
                  call dampot(vtot,drvxtot,drvytot,drvztot,dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,x1,y1,z1)
                  end if
                  !rationalise
                  drsq=1.0/((drvxtot**2)+(drvytot**2)+(drvztot**2)) !for normalizing gradient  
                  if (cnt1==0) xpp=x1 !storing previous point in xpp
                  if (cnt1==0) ypp=y1 !storing previous point in xpp
                  if (cnt1==0) zpp=z1 !storing previous point in xpp
                  cnt=cnt+1
                  cnt1=cnt1+1
                  X1 = X1 - (step*DRVXTOT*sqrt(drsq))
                  Y1 = Y1 - (step*DRVYTOT*sqrt(drsq))
                  Z1 = Z1 - (step*DRVZTOT*sqrt(drsq))
                  if ((cnt1*step)>=0.2) then ! Print points only when you have moved 0.2 au from previous point
                      noen=noen+1
                      xlines(kl,noen)=x1
                      ylines(kl,noen)=y1
                      zlines(kl,noen)=z1
                      xline(kl,noen)=x1
                      yline(kl,noen)=y1
                      zline(kl,noen)=z1
                      if (i.eq.1) then
                          xline1s(noen)=x1
                          yline1s(noen)=y1
                          zline1s(noen)=z1
                          xline1(noen)=x1
                          yline1(noen)=y1
                          zline1(noen)=z1
                          if (noen.eq.3) then  ! Noen equals 3 is chosen to ensure that the path has acquired a considerable bend
                            Vecx = x1-xcp(cpnum)
                            Vecy = y1-ycp(cpnum)
                            Vecz = z1-zcp(cpnum)
                            tang = acos((Vecx*vec(1,1)+Vecy*vec(2,1) &
                                & +Vecz*vec(3,1))/((Vecx**2+Vecy**2 &
                                & +Vecz**2)*(Vec(1,1)**2+Vec(2,1)**2 &
                                & +Vec(3,1)**2)))
                            tang = tang/rad
                          end if    
                      end if    
                      
                     cnt1=0
                  end if
                   
                  ! xpp,ypp,zpp,xp,yp and zp are also transferred
                  ! to terminate subroutine.
!               modified by Rafa 2019  ------------------------------------------------------------
                  if (i .eq. 1 .or. mesp) then
                        call basinstop(vtot, drvxtot, drvytot, drvztot, x1, y1, z1)
                        if (med .and. flag .ne. 0) noenmax = noen   ! modified by Rafa 2019
                  endif
!               end of modified by Rafa 2019  -----------------------------------------------------
               end do   
!=====================================================
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!=====================================================
!               IF (SOLIDSURF) THEN     ! modified by Rafa 2019
               if (i.eq.1) then
               
                    ndotss(kl)=noen
                    ndots1s=noen
!                elseif (i.gt.1.and.i.lt.iars) then
               elseif (i.gt.1.and.i.le.iars) then   ! modified by Rafa 2019
               
                    ndotss(kl)=noen

                    if (ndotss(1).eq.0.or.ndotss(2).eq.0) go to 95
                    call vertex(xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & xlines(2,1),ylines(2,1),zlines(2,1),&
                              & xlines(1,1),ylines(1,1),zlines(1,1),&
                              & vnx(1),vny(1),vnz(1),tang)
                    call vertex(xlines(2,1),ylines(2,1),zlines(2,1),&
                              & xlines(1,1),ylines(1,1),zlines(1,1),&
                              & xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & vnx(2),vny(2),vnz(2),tang) 
                    call vertex(xlines(1,1),ylines(1,1),zlines(1,1),&
                              & xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & xlines(2,1),ylines(2,1),zlines(2,1),&
                              & vnx(3),vny(3),vnz(3),tang) 
                    IF (SOLIDSURF) write(ibils,89) xcp(cpnum),ycp(cpnum),zcp(cpnum),&   ! modified by Rafa 2019
                              & xlines(2,1),ylines(2,1),zlines(2,1),&
                              & xlines(1,1),ylines(1,1),zlines(1,1),&
                              & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                              & vnx(3),vny(3),vnz(3)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
!        Writes vertices and indices to basisfile for DAMQT
                    if (ierrqt .eq. 0) then
                        write(ibilqt+1,*) xcp(cpnum),ycp(cpnum),zcp(cpnum), vnx(1), vny(1), vnz(1), &
                            xlines(2,1),ylines(2,1),zlines(2,1), vnx(2), vny(2), vnz(2), &
                            xlines(1,1),ylines(1,1),zlines(1,1), vnx(3), vny(3), vnz(3)
                        kntvert = kntvert + 3
                    endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------
                    
                    do ins = 1, min(ndotss(1),ndotss(2))-1

                        lenths = dsqrt((xlines(1,ins)-xlines(2,ins))**2+ &
                        & (ylines(1,ins)-ylines(2,ins))**2+(zlines(1,ins)-zlines(2,ins))**2)   
                        call vertex(xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),tang)
                        call vertex(xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & vnx(2),vny(2),vnz(2),tang)
                        call vertex(xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                        if (lenths.le.exln) then
                            IF (SOLIDSURF) write(ibils,89) xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                  & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                            if (ierrqt .eq. 0) then
                                write(ibilqt+1,*) xlines(1,ins),ylines(1,ins),zlines(1,ins), vnx(1), vny(1), vnz(1), &
                                    xlines(2,ins),ylines(2,ins),zlines(2,ins), vnx(2), vny(2), vnz(2), &
                                    xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1), vnx(3), vny(3), vnz(3)
                                kntvert = kntvert + 3
                            endif
                        endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                        call vertex(xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),tang)
                        call vertex(xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & vnx(2),vny(2),vnz(2),tang)
                        call vertex(xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                  & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                        if (lenths.le.exln) then
                            IF (SOLIDSURF) write(ibils,89) xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                  & xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                  & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                            if (ierrqt .eq. 0) then
                                write(ibilqt+1,*) xlines(2,ins),ylines(2,ins),zlines(2,ins), vnx(1), vny(1), vnz(1), &
                                    xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1), vnx(2), vny(2), vnz(2), &
                                    xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1), vnx(3), vny(3), vnz(3)
                                kntvert = kntvert + 3
                            endif
                      endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                    end do
                        
                    if (ndotss(1).eq.ndotss(2)) go to 96
                    
                    if (exdraw) then
                    if (ndotss(1).eq.min(ndotss(1),ndotss(2))) then 
                        do ins = min(ndotss(1),ndotss(2)),max(ndotss(1),ndotss(2))-1
                            lenths = dsqrt((xlines(1,ndotss(1))-xlines(2,ins))**2+ &
                            & (ylines(1,ndotss(1))-ylines(2,ins))**2+(zlines(1,ndotss(1))-zlines(2,ins))**2)   
                            call vertex(xlines(1,ndotss(1)),ylines(1,ndotss(1)),zlines(1,ndotss(1)),&
                                      & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                      & xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                      & vnx(1),vny(1),vnz(1),tang)
                            call vertex(xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                      & xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                      & xlines(1,ndotss(1)),ylines(1,ndotss(1)),zlines(1,ndotss(1)),&
                                      & vnx(2),vny(2),vnz(2),tang)
                            call vertex(xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                      & xlines(1,ndotss(1)),ylines(1,ndotss(1)),zlines(1,ndotss(1)),&
                                      & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                      & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                            if (lenths.le.exln) then
                                IF (SOLIDSURF) write(ibils,89) xlines(1,ndotss(1)),ylines(1,ndotss(1)),&
                                      & zlines(1,ndotss(1)),&
                                      & xlines(2,ins),ylines(2,ins),zlines(2,ins),&
                                      & xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1),&
                                      & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                      & vnx(3),vny(3),vnz(3)

!        Writes vertices and indices to basisfile for DAMQT
                                  if (ierrqt .eq. 0) then
                                      write(ibilqt+1,*) xlines(1,ndotss(1)),ylines(1,ndotss(1)), zlines(1,ndotss(1)), &
                                          vnx(1), vny(1), vnz(1), &
                                          xlines(2,ins),ylines(2,ins),zlines(2,ins), vnx(2), vny(2), vnz(2), &
                                          xlines(2,ins+1),ylines(2,ins+1),zlines(2,ins+1), vnx(3), vny(3), vnz(3)
                                      kntvert = kntvert + 3
                                  endif
                              endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                        end do        
                    elseif (ndotss(2).eq.min(ndotss(1),ndotss(2))) then 
                        do ins = min(ndotss(1),ndotss(2)),max(ndotss(1),ndotss(2))-1
                            lenths = dsqrt((xlines(2,ndotss(2))-xlines(1,ins))**2+ &
                            & (ylines(2,ndotss(2))-ylines(1,ins))**2+(zlines(2,ndotss(2))-zlines(1,ins))**2)   
                            call vertex(xlines(2,ndotss(2)),ylines(2,ndotss(2)),zlines(2,ndotss(2)),&
                                      & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & vnx(1),vny(1),vnz(1),tang)
                            call vertex(xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & xlines(2,ndotss(2)),ylines(2,ndotss(2)),zlines(2,ndotss(2)),&
                                      & vnx(2),vny(2),vnz(2),tang)
                            call vertex(xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & xlines(2,ndotss(2)),ylines(2,ndotss(2)),zlines(2,ndotss(2)),&
                                      & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                            if (lenths.le.exln) then
                                IF (SOLIDSURF) write(ibils,89) xlines(2,ndotss(2)),ylines(2,ndotss(2)),&
                                      & zlines(2,ndotss(2)),&
                                      & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                      & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                                  if (ierrqt .eq. 0) then
                                      write(ibilqt+1,*) xlines(2,ndotss(2)),ylines(2,ndotss(2)),zlines(2,ndotss(2)), &
                                          vnx(1), vny(1), vnz(1), &
                                          xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1), vnx(2), vny(2), vnz(2), &
                                          xlines(1,ins),ylines(1,ins),zlines(1,ins), vnx(3), vny(3), vnz(3)
                                      kntvert = kntvert + 3
                                  endif
                            endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                        end do 
                    end if
                    end if
96                  continue
                    
                    do ins = 1,ndotss(1)
                        xlines(1,ins)=0.0
                        ylines(1,ins)=0.0
                        zlines(1,ins)=0.0
                    end do    
                    do ins = 1,ndotss(2) 
                        xlines(1,ins)=xlines(2,ins)
                        ylines(1,ins)=ylines(2,ins)
                        zlines(1,ins)=zlines(2,ins)
                    end do
                    ndotss(1)=ndotss(2)
               end if 
               
               if (i.eq.iars) then
                    ndotss(kl)=noen
                    if (ndots1s.eq.0.or.ndotss(1).eq.0) go to 95
                    call vertex(xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & xline1s(1),yline1s(1),zline1s(1),&
                              & xlines(1,1),ylines(1,1),zlines(1,1),&
                              & vnx(1),vny(1),vnz(1),tang)
                    call vertex( xline1s(1),yline1s(1),zline1s(1),&
                              & xlines(1,1),ylines(1,1),zlines(1,1),&
                              & xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & vnx(2),vny(2),vnz(2),tang)
                    call vertex(xlines(1,1),ylines(1,1),zlines(1,1),&
                              & xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & xline1s(1),yline1s(1),zline1s(1),&
                              & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                    IF (SOLIDSURF) write(ibils,89) xcp(cpnum),ycp(cpnum),zcp(cpnum),&
                              & xline1s(1),yline1s(1),zline1s(1),&
                              & xlines(1,1),ylines(1,1),zlines(1,1),&
                              & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                              & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                    if (ierrqt .eq. 0) then
                        write(ibilqt+1,*) xcp(cpnum),ycp(cpnum),zcp(cpnum), vnx(1), vny(1), vnz(1), &
                            xline1s(1),yline1s(1),zline1s(1), vnx(2), vny(2), vnz(2), &
                            xlines(1,1),ylines(1,1),zlines(1,1), vnx(3), vny(3), vnz(3)
                        kntvert = kntvert + 3
                    endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------
                    
                    do ins = 1, min(ndots1s,ndotss(1))-1
                        lenths = dsqrt((xline1s(ins)-xlines(1,ins))**2 + &
                        & (yline1s(ins)-ylines(1,ins))**2 + (zline1s(ins)-zlines(1,ins))**2)
                        call vertex(xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),tang)
                        call vertex(xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & vnx(2),vny(2),vnz(2),tang)
                        call vertex(xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                        if (lenths.le.exln) then
                            IF (SOLIDSURF) write(ibils,89) xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                  & xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                  & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                                if (ierrqt .eq. 0) then
                                    write(ibilqt+1,*) xlines(1,ins),ylines(1,ins),zlines(1,ins), vnx(1), vny(1), vnz(1), &
                                        xline1s(ins),yline1s(ins),zline1s(ins), vnx(2), vny(2), vnz(2), &
                                        xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1), vnx(3), vny(3), vnz(3)
                                    kntvert = kntvert + 3
                                endif
                          endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                        call vertex(xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),tang)
                        call vertex(xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & vnx(2),vny(2),vnz(2),tang)
                        call vertex(xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                  & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                        if (lenths.le.exln) then
                            IF (SOLIDSURF) write(ibils,89) xline1s(ins),yline1s(ins),zline1s(ins),&
                                  & xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                  & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                  & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                  & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                              if (ierrqt .eq. 0) then
                                  write(ibilqt+1,*) xline1s(ins),yline1s(ins),zline1s(ins), vnx(1), vny(1), vnz(1), &
                                      xline1s(ins+1),yline1s(ins+1),zline1s(ins+1), vnx(2), vny(2), vnz(2), &
                                      xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1), vnx(3), vny(3), vnz(3)
                                  kntvert = kntvert + 3
                              endif
                        endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                    end do
                    if (ndots1s.eq.ndotss(1)) go to 95
                    
                    if (exdraw) then
                    if (ndotss(1).eq.min(ndots1s,ndotss(1))) then
                        do ins = min(ndots1s,ndotss(1)),max(ndots1s,ndotss(1))-1
                            lenths = dsqrt((xline1s(ins)-xlines(1,ndotss(1)))**2 + &
                            & (yline1s(ins)-ylines(1,ndotss(1)))**2 + (zline1s(ins)-zlines(1,ndotss(1)))**2)
                            call vertex(xlines(1,ndotss(1)),ylines(1,ndotss(1)),zlines(1,ndotss(1)),&
                                      & xline1s(ins), yline1s(ins), zline1s(ins),&
                                      & xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                      & vnx(1),vny(1),vnz(1),tang)
                            call vertex(xline1s(ins), yline1s(ins), zline1s(ins),&
                                      & xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                      & xlines(1,ndotss(1)),ylines(1,ndotss(1)),zlines(1,ndotss(1)),&
                                      & vnx(2),vny(2),vnz(2),tang)
                            call vertex(xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                      & xlines(1,ndotss(1)),ylines(1,ndotss(1)),zlines(1,ndotss(1)),&
                                      & xline1s(ins), yline1s(ins), zline1s(ins),&
                                      & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                            if (lenths.le.exln) then
                                IF (SOLIDSURF) write(ibils,89) xlines(1,ndotss(1)),ylines(1,ndotss(1)),&
                                      & zlines(1,ndotss(1)),&
                                      & xline1s(ins), yline1s(ins), zline1s(ins),&
                                      & xline1s(ins+1),yline1s(ins+1),zline1s(ins+1),&
                                      & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                      & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                                if (ierrqt .eq. 0) then
                                    write(ibilqt+1,*) xlines(1,ndotss(1)),ylines(1,ndotss(1)), zlines(1,ndotss(1)), &
                                        vnx(1), vny(1), vnz(1), &
                                        xline1s(ins), yline1s(ins), zline1s(ins), vnx(2), vny(2), vnz(2), &
                                        xline1s(ins+1),yline1s(ins+1),zline1s(ins+1), vnx(3), vny(3), vnz(3)
                                    kntvert = kntvert + 3
                                endif
                            endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                        end do
                    elseif (ndots1s.eq.min(ndots1s,ndotss(1))) then
                        do ins = min(ndots1s,ndotss(1)),max(ndots1s,ndotss(1))-1
                            lenths = dsqrt((xline1s(ndots1s)-xlines(1,ins))**2 + &
                            & (yline1s(ndots1s)-ylines(1,ins))**2 + (zline1s(ndots1s)-zlines(1,ins))**2)
                            call vertex(xline1s(ndots1s),yline1s(ndots1s),zline1s(ndots1s),&
                                      & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & vnx(1),vny(1),vnz(1),tang)
                            call vertex(xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & xline1s(ndots1s),yline1s(ndots1s),zline1s(ndots1s),&
                                      & vnx(2),vny(2),vnz(2),tang)
                            call vertex(xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & xline1s(ndots1s),yline1s(ndots1s),zline1s(ndots1s),&
                                      & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & vnx(3),vny(3),vnz(3),tang)

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
                            if (lenths.le.exln) then
                                IF (SOLIDSURF) write(ibils,89) xline1s(ndots1s),yline1s(ndots1s),&
                                      & zline1s(ndots1s),&
                                      & xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1),&
                                      & xlines(1,ins),ylines(1,ins),zlines(1,ins),&
                                      & vnx(1),vny(1),vnz(1),vnx(2),vny(2),vnz(2),&
                                      & vnx(3),vny(3),vnz(3)
!        Writes vertices and indices to basisfile for DAMQT
                              if (ierrqt .eq. 0) then
                                  write(ibilqt+1,*) xline1s(ndots1s),yline1s(ndots1s), zline1s(ndots1s), &
                                      vnx(1), vny(1), vnz(1), &
                                      xlines(1,ins+1),ylines(1,ins+1),zlines(1,ins+1), vnx(2), vny(2), vnz(2), &
                                      xlines(1,ins),ylines(1,ins),zlines(1,ins), vnx(3), vny(3), vnz(3)
                                  kntvert = kntvert + 3
                              endif
                           endif
!        End of modification by Rafa 2019  --------------------------------------------------------------------------------

                        end do
                    end if
                    end if
               end if
95             continue
               
!               ENDIF   !If loop for Solidsurf      Modified by Rafa 2019
                        
                        
               IF (WIREFRAME) THEN    
               
               if (i.eq.1) then
                    ndots(kl)=noen
                    ndots1=noen  ! Treating very first array differently
               elseif (i.gt.1.and.i.le.iars) then
                    ndots(kl)=noen
                    if (ndots(1).eq.0.or.ndots(2).eq.0) go to 985
                    write(ibilw,90)xcp(cpnum),ycp(cpnum),zcp(cpnum)
                    write(ibilw,90)xline(1,1),yline(1,1),zline(1,1) 
                    write(ibilw,*)""
                    write(ibilw,90)xcp(cpnum),ycp(cpnum),zcp(cpnum)
                    write(ibilw,90)xline(2,1),yline(2,1),zline(2,1)
                    write(ibilw,*)""
                  
                    !Writing the point on a line connecting two CPs
                    do jns=1,2
                    do ins = 1,ndots(jns)
                        write(ibilw,90)xline(jns,ins),yline(jns,ins),zline(jns,ins)
                    end do  
                    write(ibilw,*)""
                    end do
                    write(ibilw,*)""

                    ! In case two adjacent lines possess unequal number
                    ! of points, minimum among them is used to connect
                    ! them transversely.
                    do ins = 1, min(ndots(1),ndots(2))
                        lenth = dsqrt((xline(1,ins)-xline(2,ins))**2+ &
                        & (yline(1,ins)-yline(2,ins))**2+(zline(1,ins)-zline(2,ins))**2)   
                        if (lenth.le.exln) write(ibilw,90) xline(1,ins),yline(1,ins),zline(1,ins)
                        if (lenth.le.exln) write(ibilw,90) xline(2,ins),yline(2,ins),zline(2,ins)
                        if (lenth.le.exln) write(ibilw,*)""
                    end do
                    
                    ! In both the adjacent lines have same number of
                    ! points no need to write any extra connection    
                    if (ndots(1).eq.ndots(2)) go to 986
                    
                    ! In case above is not true then depending on
                    ! which is larger, all the extra points on larger
                    ! line are connected from the terminal point of
                    ! smaller line
                    if (exdraw) then    
                    if (ndots(1).eq.min(ndots(1),ndots(2))) then 
                        do ins = min(ndots(1),ndots(2)),max(ndots(1),ndots(2))   
                            lenth = dsqrt((xline(1,ndots(1))-xline(2,ins))**2+ &
                            & (yline(1,ndots(1))-yline(2,ins))**2+(zline(1,ndots(1))-zline(2,ins))**2)   
                            if (lenth.le.exln) write(ibilw,90) xline(1,ndots(1)),yline(1,ndots(1)),zline(1,ndots(1))
                            if (lenth.le.exln) write(ibilw,90) xline(2,ins),yline(2,ins),zline(2,ins)
                            if (lenth.le.exln) write(ibilw,*)""
                        end do        
                    elseif (ndots(2).eq.min(ndots(1),ndots(2))) then 
                        do ins = min(ndots(1),ndots(2)),max(ndots(1),ndots(2))   
                            lenth = dsqrt((xline(2,ndots(2))-xline(1,ins))**2+ &
                            & (yline(2,ndots(2))-yline(1,ins))**2+(zline(2,ndots(2))-zline(1,ins))**2)   
                            if (lenth.le.exln) write(ibilw,90) xline(2,ndots(2)),yline(2,ndots(2)),zline(2,ndots(2))
                            if (lenth.le.exln) write(ibilw,90) xline(1,ins),yline(1,ins),zline(1,ins)
                            if (lenth.le.exln) write(ibilw,*)""
                        end do 
                    end if
                    end if

986                 continue
                    ! Swapping the array elements     
                    do ins = 1,ndots(1)
                         xline(1,ins)=0.0
                         yline(1,ins)=0.0
                         zline(1,ins)=0.0
                    end do    
                    do ins = 1,ndots(2) 
                         xline(1,ins)=xline(2,ins)
                         yline(1,ins)=yline(2,ins)
                         zline(1,ins)=zline(2,ins)
                    end do
                    ndots(1)=ndots(2)
                end if

985             continue

                if (i.eq.iars) then
                    ndots(kl)=noen
                    if (ndots1.eq.0.or.ndots(1).eq.0) go to 987

                    do ins = 1, ndots(1)
                        write(ibilw,90) xline(1,ins),yline(1,ins),zline(1,ins)
                    end do
                    write(ibilw,*)""
                            
                    do ins = 1, min(ndots1,ndots(1))
                            lenth = dsqrt((xline1(ins)-xline(1,ins))**2 + &
                            & (yline1(ins)-yline(1,ins))**2 + (zline1(ins)-zline(1,ins))**2)   
                            if (lenth.le.exln) write(ibilw,90) xline1(ins),yline1(ins),zline1(ins)
                            if (lenth.le.exln) write(ibilw,90) xline(1,ins),yline(1,ins),zline(1,ins)
                            if (lenth.le.exln) write(ibilw,*)""
                    end do
                   
                    if (ndots1.eq.ndots(1)) go to 987
                    if (exdraw) then 
                    if (ndots1.eq.min(ndots1,ndots(1))) then 
                        do ins = min(ndots1,ndots(1)),max(ndots1,ndots(1))   
                            lenth = dsqrt((xline1(ndots1)-xline(1,ins))**2 + &
                            & (yline1(ndots1)-yline(1,ins))**2 + (zline1(ndots1)-zline(1,ins))**2)   
                            if (lenth.le.exln) write(ibilw,90) xline1(ndots1),yline1(ndots1),zline1(ndots1)
                            if (lenth.le.exln) write(ibilw,90) xline(1,ins),yline(1,ins),zline(1,ins)
                            if (lenth.le.exln) write(ibilw,*)""
                        end do        
                    elseif (ndots(1).eq.min(ndots1,ndots(1))) then 
                        do ins = min(ndots1,ndots(1)),max(ndots1,ndots(1))   
                            lenth = dsqrt((xline1(ins)-xline(1,ndots(1)))**2 + &
                            & (yline1(ins)-yline(1,ndots(1)))**2 + (zline1(ins)-zline(1,ndots(1)))**2)   
                            if (lenth.le.exln) write(ibilw,90) xline1(ins),yline1(ins),zline1(ins)
                            if (lenth.le.exln) write(ibilw,90) xline(1,ndots(1)),yline(1,ndots(1)),zline(1,ndots(1))
                            if (lenth.le.exln) write(ibilw,*)""
                        end do 
                    end if
                    end if   
                end if  
987             continue   
                END IF ! If loop for wireframe
!====================================================================
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!====================================================================          
          end do !Loop over number of points around a CP
      end do ! Loop over number of CPs

!        Modified by Rafa 2019  -------------------------------------------------------------------------------------------
    if (ierrqt .eq. 0) then
      if (med) basinsname = trim(filename)//"-d.basins"
      if (mesp) basinsname = trim(filename)//"-v.basins"
#if _WIN32
      open (unit=IBILQT, file=basinsname,status="replace", form='binary', carriagecontrol='NONE', iostat=ierrqt)
      if (ierrqt .ne. 0) write(6,*) 'Error when opening file ', basinsname
#elif __INTEL_COMPILER
      open (unit=IBILQT, file=basinsname,status="replace", form='binary', carriagecontrol='NONE', iostat=ierrqt)
      if (ierrqt .ne. 0) write(6,*) 'Error when opening file ', basinsname
#else
       open (unit=IBILQT, file=basinsname,status="replace", form='unformatted', access='stream', iostat=ierr)
       if (ierrqt .ne. 0) write(6,*) 'Error when opening file ', basinsname
#endif

    endif

    if (ierrqt .eq. 0) then
        rewind (IBILQT+1)
        write(IBILQT) kntvert
        do kntind = 1, kntvert / 3
            read(IBILQT+1,*) xyz8(1:18)
            do i = 1, 13, 6
                x4 = xyz8(i)
                y4 = xyz8(i+1)
                z4 = xyz8(i+2)
                dx4 = xyz8(i+3)
                dy4 = xyz8(i+4)
                dz4 = xyz8(i+5)
                write(IBILQT) x4, y4, z4, dx4, dy4, dz4
            enddo
        enddo
    endif
    close(IBILQT)
    close(IBILQT+1)
    CALL SYSTEM("rm -f "//trim(filename)//".vert")

!        End of modification by Rafa 2019  ----------------------------------------------------------------------------------

89    format(9(e11.3,1x),9(f7.3,1x))        
90    format(3(e11.3,1x))        
      close(ibils)
      close(ibilw)
      
      RETURN
      END SUBROUTINE 
!---------------------------------------------------------------------------------

!> @details Decides the box size for terminating the creation of basin
!*********************************************************************************     
      SUBROUTINE BASINSTOP(vtot, drvxtot, drvytot, drvztot, x1, y1,z1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USE GRADPATH_T
    IMPLICIT NONE
    REAL(KREAL) :: vtot, drvxtot, drvytot, drvztot, x1, y1, z1, distcpnew ! Modified by Rafa 2019
    INTEGER(KINT) :: I,J,K,L

    distcpnew = 0.01d0 * distcp ! Modified by Rafa 2019
      
    if (med) then
        if (vtot .lt. 1.d-3) then
            FLAG=1
            return
        end if
    elseif (mesp) then
        if ((X1.ge.XBOXGE).OR.(X1.le.XBOXGS).OR.(Y1.ge.YBOXGE).OR.(Y1.le.YBOXGS).OR.(Z1.ge.ZBOXGE).OR.(Z1.le.ZBOXGS)) then
            FLAG=1
            return
        end if
    end if


    if ((abs(drvxtot)<drcutcp).AND.(abs(drvytot)<drcutcp).AND.(abs(drvztot)<drcutcp)) then
    ! if all the derivatives are smaller than 10e-5

        do l=1, NCP ! To Check if the gradient path is in close proximity of a CP.
            if (l.ne.cpnum) then ! if current CP in the loop is not the
                                 ! same CP from which gradient path has originated then
               if ((abs(X1-XCP(l))<distcpnew).AND.(abs(Y1-YCP(l))<distcpnew).AND.(abs(Z1-ZCP(l))<distcpnew)) then  ! Modified by Rafa 2019
! And the gradient path has reached within 0.05 au or less to another CP then
!                 write(iout,"(2(a,1x,a,1x,I3,1x),I7)") "Gradient path connects",symcp(cpnum),cpnum,"to",symcp(l),l, cnt
!                  write(IRMP,*) symcp(cpnum),cpnum,symcp(l),l,cnt
                    FLAG=1
                    return
               end if
            elseif (l.eq.cpnum) then
                if ((abs(X1-XCP(l))<distcpnew).AND.(abs(Y1-YCP(l))<distcpnew).AND.(abs(Z1-ZCP(l))<distcpnew)) then  ! Modified by Rafa 2019
                    FLAG=1
!                   write(iout,*) "returned to same cp"
                    return
                end if
            end if
        end do
    end if


      if ((abs(drvxtot)>drcutnuc).OR.(abs(drvytot)>drcutnuc).OR.(abs(drvztot)>drcutnuc)) then
         do l=1, NCEN
            if ((abs(X1-RCEN(1,l))<DISTnuc).AND.(abs(Y1-RCEN(2,l))<DISTnuc).AND.(abs(Z1-RCEN(3,l))<DISTnuc)) then 
!              write(irmp,*) symcp(cpnum),cpnum,atmnms(int(NCN(l))),l,cnt 
               FLAG=1
               return
            end if
         end do
      end if   
             
      !If the path finds a loop
              
      if((ABS(xpp-x1).lt.1.0E-07).AND.(ABS(ypp-y1)<1.0E-07)  &
         &  .AND.(ABS(zpp-z1)<1.0E-07)) then
!         write(iout,*) "The path stuck in a loop at", xpp, ypp, zpp, &
!            &   "and", x1, y1, z1
         FLAG=1
         return
      end if
      FLAG=0

      END SUBROUTINE

!> @details Information Subroutine
!**************************************************      
      SUBROUTINE TITLES
!**************************************************      
      USE DAMINITIAL_T, ONLY: iout,GRADPATH,TOPOGRAPH,MED,MESP,BASIN, PROJECTNAME, FILENAME
      implicit none
      character(len=200):: cwd
! ---------------     End of change made by Rafa   -------------- 

      write(iout,"(a)")"****************************************"    
      write(iout,"(a)")"Topographer (Version 15.05, Academic)"    
      write(iout,"(a)")"Written by Anmol Kumar, Sachin D. Yeole and Shridhar R. Gadre"    
      write(iout,"(a)")"Topographer utilizies damqt2.0 subroutines for evaluation of scalar field."
      write(iout,"(a)")"(Molecular electron density and Molecular electrostatic Potential)"    
      write(iout,"(a)")"Copyright 2008-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,"
      write(iout,"(a)")"Guillermo Ramirez, David Zorrilla, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre"
      write(iout,"(a)")"Topographer utilizes the L-BFGS-B version lbfgs3.0"
      write(iout,"(a)")"optimization package http://users.iems.northwestern.edu/~nocedal/lbfgsb.html)"
      write(iout,"(a)")"for fast estimation of critical points."    
      write(iout,"(a)")""    
      write(iout,"(a)")""    
!       CALL GETCWD(CWD)	! Changed by Rafa
	cwd = ""		! Changed by Rafa
! write(6,*) 'cwd = ', cwd
! write(6,*) 'filename = ', filename
! ---------------     change made by Rafa   --------------------
      
      write(iout,"(5(a))")"Input file: ",TRIM(CWD),"/",TRIM(PROJECTNAME),".fchk"
      write(iout,"(5(a))")"Input file: ", TRIM(CWD),"/",TRIM(PROJECTNAME),".den"
      write(iout,"(5(a))")"Input file: ", TRIM(CWD),"/",TRIM(PROJECTNAME),".ggbs"
      IF (TOPOGRAPH.AND.MED) write(iout,"(5(a))")"CPs stored in: ",TRIM(CWD),"/",TRIM(FILENAME),"-cps-d.xyz"
      IF (TOPOGRAPH.AND.MESP) write(iout,"(5(a))")"CPs stored in: ",TRIM(CWD),"/",TRIM(FILENAME),"-cps-v.xyz"
      IF (GRADPATH) write(iout,"(5(a))")"Gradient path in: ",TRIM(CWD),"/",TRIM(FILENAME),".gpdat"
      IF (BASIN) write(iout,"(5(a))")"Basins in: ",TRIM(CWD),"/",TRIM(FILENAME),".bsdat"
      
! ---------------     End of change made by Rafa   --------------       
      
      write(iout,*)""    
      write(iout,"(a)")"Atomic units (a.u.) used throughout, &
                        & except where stated otherwise"
      write(iout,"(a)")"Symbols viz. x,y,z and m (small letters) are used to denote the"    
      write(iout,"(a)")"nature of critical points."    
      write(iout,"(a)")"x implies (3,+3) CP"    
      write(iout,"(a)")"y implies (3,+1) CP"    
      write(iout,"(a)")"z implies (3,-1) CP"    
      write(iout,"(a)")"m implies (3,-3) CP"    
      write(iout,*)""    
      write(iout,*)""    
      RETURN
      ENDSUBROUTINE
!===================================================
      subroutine vertex(px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,verxn,veryn,verzn,tang)
      IMPLICIT NONE
      REAL*8:: px1,py1,pz1,px2,py2,pz2,px3,py3,pz3,verxn,veryn,verzn,tang
      REAL*8:: vx1,vy1,vz1,vx2,vy2,vz2,verx,very,verz,lver
         vx1=px2-px1
         vy1=py2-py1
         vz1=pz2-pz1
         vx2=px3-px1             
         vy2=py3-py1
         vz2=pz3-pz1
         verx = vy1 * vz2 - vz1 * vy2
         very = vz1 * vx2 - vx1 * vz2
         verz = vx1 * vy2 - vy1 * vx2
         lver = sqrt(verx**2+very**2+verz**2)
         if (tang.lt.90) then
            verxn=verx/lver 
            veryn=very/lver
            verzn=verz/lver
         else
            verxn=-verx/lver 
            veryn=-very/lver
            verzn=-verz/lver
         end if
         return
      end subroutine
               
      
