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
!> @file TDAMOPTIMIZER.F90

!> @authors Anmol Kumar and Shridhar R. Gadre
!> @date 01-01-2021
!
!  Modified by Rafael López
!   DATE: 23-04-2021
!
    Program BUILDDRIVER
    USE DAMINITIAL_T
    USE DAMBUILD_T
    implicit none
    integer(kint)::iden,ierr,i,j,k,iw,jw,dum,nbatms,nb,flag
    real :: tarray(2), tiempo, dtime
    integer:: nmols,splt,rssize,tn,tmaxcharg
    real(kreal) :: Enj,tssize
    real(kreal) :: shiftx,shifty,shiftz,tcox,tcoy,tcoz,chosenx,choseny,chosenz
    character(2),allocatable::pisym(:),ts(:)
    real(kreal),allocatable ::tx(:),ty(:),tz(:),tq(:),tw(:),tv(:),pix(:),piy(:),piz(:)
    character(256)::exec
    character(256)::acmdout
    character(256)::mstring
    character(256):: preprocfile,templatefile,insertlocfile,qmsoftware,qmkeywords
    character(256):: mespimizerinit,mespimizervis,mespimizerfinal,basename,mespimizerbase
    character(256):: clustername, path
    character(256):: mespimizerauxframes, mespimizerkntfrms, mespimizerframes, mespimizercurrframe
    character(10)::d1,d2,d3,d4,extn
    integer::icmdout
    logical:: nocharge, lwriteqm
!       character, external:: ncase
    namelist / options / lgradient, lderiv2, largo, lexact, &
            & lmaxi,filename,preprocfile,templatefile,insertlocfile,&
            & iswindows,tssize,rssize,nocharge,qmsoftware,qmkeywords, &
            & clustername, path, lwriteqm

    tiempo = dtime(tarray)
    lgradient    = .true.
    lderiv2      = .true.
    largo        = .false.        ! IF .TRUE. LONG-RANGE POTENTIAL
    lexact       = .false.        ! IF .TRUE. "EXACT" POTENTIAL IS TABULATED
    iswindows    = .false.
    lmaxi        = 15
    path         = ""
    clustername  = "cluster"
    filename     = ""
    preprocfile  = ""
    templatefile = ""
    insertlocfile= ""
    tssize       = 0.5
    rssize       = 20
    nocharge     = .true.
    lwriteqm     = .false.
    qmsoftware   = "gaussian"
    qmkeywords   = "nproc=8,mem=2000Mb,lot=MP2,basis=6-31+G(d),charge=0,multi=1"

    read(5,options)     ! Read from standard input
    read(5,*) projectname

    ! Decides if the user has provided projectname with full path.
    if (iswindows) then
        dirsep = "\\"
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        if (i .eq. 0) then  ! This is intended for MinGW, whose directory separator in windows is also /
            dirsep = "/"
            i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
        endif
    else
        dirsep = "/"
        i = index(projectname,dirsep,.true.)    ! Checks position of last directory name separator
    end if
    if (len_trim(filename).eq.0) then
        filename = projectname
    else
        filename = projectname(1:i)//trim(filename)
    endif
    i = index(path,'/')
    if (i==0) then
        path = trim(adjustl(projectname(1:i)))
    endif
    if (iswindows) then
        do j = 1, len_trim(filename)
            i = index(filename,'/') ; if (i == 0) exit
            filename = filename(:i-1)//dirsep//filename(i+1:)
        end do
        do j = 1, len_trim(path)
            i = index(path,'/') ; if (i == 0) exit
            path = path(:i-1)//dirsep//path(i+1:)
        end do
    endif

    preprocfile   = trim(adjustl(preprocfile))
    templatefile  = trim(adjustl(templatefile))
    insertlocfile = trim(adjustl(insertlocfile))

    ! Figure out what files are provided
    if (len_trim(templatefile)==0 .and. len_trim(insertlocfile)==0 &
        .and. len_trim(preprocfile)==0) then
        write(6,*) "Provide preprocfile or templatefile with &
                    & insertlocfile"
        stop
    endif
    if (len(trim(preprocfile))/=0) then
        write(6,*) "Preprocfile is provided. Templatefile &
                    & and insertlocfile will be ignored"
        write(6,*) "Note: Preprocfile should not contain &
                    & parent molecule. "
        i = index(preprocfile,dirsep)
        if (i == 0) preprocfile = trim(path)//trim(preprocfile)
        splt = index(preprocfile,'.',back=.True.)
        basename=preprocfile(1:splt-1)
        extn=preprocfile(splt:)

    else
        if (len_trim(templatefile)/=0 .and. len_trim(insertlocfile)==0) then
            write(6,*) "Templatefile requires insertlocfile to &
                    &place the  template molecule. Exiting code."
            stop
        elseif (len_trim(templatefile)==0 .and. len_trim(insertlocfile)/=0) then
            write(6,*) "Insertlocfile requires templatefile. &
                      & Exiting code."
            stop
        elseif (len_trim(templatefile)/=0 .and. len_trim(insertlocfile)/=0) then
            i = index(templatefile,dirsep)
            if (i == 0) templatefile = trim(path)//trim(templatefile)
            i = index(insertlocfile,dirsep)
            if (i == 0) insertlocfile = trim(path)//trim(insertlocfile)
            splt = index(templatefile,'.',back=.True.)
            basename=templatefile(1:splt-1)
            extn=templatefile(splt:)
        endif
    endif
    if (dirsep /= '/') then ! Ensures coherency in dir separation character
        do j = 1, len_trim(basename)
            i = index(basename,'/') ; if (i == 0) exit
            basename = basename(:i-1) // '\' // basename(i+1:)
        end do
    endif
    mespimizerinit   = trim(adjustl(path))//trim(clustername)//".xyz_init"
    mespimizervis    = trim(adjustl(path))//trim(clustername)//".xyz_vis"
    mespimizerfinal  = trim(adjustl(path))//trim(clustername)//".xyz_final"
    mespimizerbase   = trim(adjustl(path))//trim(clustername)
    mespimizerauxframes = trim(adjustl(path))//trim(clustername)//".frames"
    mespimizerframes = trim(adjustl(path))//trim(clustername)//".xyz_frames"
    mespimizercurrframe  = trim(adjustl(path))//trim(clustername)//".xyz_curr_frame"
    mespimizerkntfrms = trim(adjustl(path))//trim(clustername)//".kntframes"

    write(6,*) 'preprocfile         = ', trim(preprocfile)
    write(6,*) 'templatefile        = ', trim(templatefile)
    write(6,*) 'insertlocfile       = ', trim(insertlocfile)
    write(6,*) 'mespimizerinit      = ', trim(mespimizerinit)
    write(6,*) 'mespimizervis       = ', trim(mespimizervis)
    write(6,*) 'mespimizerfinal     = ', trim(mespimizerfinal)
    write(6,*) 'mespimizerbase      = ', trim(mespimizerbase)
    write(6,*) 'mespimizerauxframes = ', trim(mespimizerauxframes)
    write(6,*) 'mespimizerframes    = ', trim(mespimizerframes)
    write(6,*) 'mespimizercurrframe = ', trim(mespimizercurrframe)
    write(6,*) 'mespimizerkntfrms   = ', trim(mespimizerkntfrms)

    call predampot

    ! Be it templatefile or be it preprocfile both are read same way.
    ! It could had be concise but focus is to make this work.

    if (len(trim(templatefile))/=0 .and. len(trim(preprocfile))==0) preprocfile =templatefile


    preprocfile=trim(adjustl(preprocfile))
    splt = index(preprocfile,'.',back=.True.)
    basename=preprocfile(1:splt-1)
    extn=preprocfile(splt:)
    ! If xyz file does not contain partial charges, obabel is used to get initial charges.
    if (nocharge) then
        if (iswindows) then
!            exec = '"c:' // dirsep // "Program Files" // dirsep // "OpenBabel-3.1.1" // dirsep // &
            exec = &
            & '"obabel.exe" -ixyz '//trim(adjustl(basename))//trim(adjustl(extn)) &
            & //" --partialcharge qtpie -omol2 -O " &
            & //trim(adjustl(basename))//".mol2"
            print*,exec
            call system(exec, icmdout)
            !write(6,*) 'llama a execute_command_line'
            !call execute_command_line(exec,cmdstat=icmdout,cmdmsg=acmdout)
        else
            exec = "obabel -ixyz "//trim(adjustl(basename))//trim(adjustl(extn)) &
            & //" --partialcharge qtpie -omol2 -O " &
            & //trim(adjustl(basename))//".mol2"
            print*,exec
            call execute_command_line(exec,cmdstat=icmdout,cmdmsg=acmdout)
        endif

        if (icmdout/=0) then
            print *,"Exiting Code: ", icmdout
            if (iswindows) then
                print *,"Exiting Code: ", icmdout
            else
                print *,"Exiting Code: ", icmdout, ' Exiting message: ', acmdout
            endif
            print*, "Maybe input xyz file does not contain charge &
                &and openbabel is not installed on your computer."
            print *,"If you believe you have charges in input xyz file, insert nocharge = .false. in options and rerun."
            stop
        else
            !Read mol2 file
            open(unit=298, file=trim(adjustl(basename))//".mol2",iostat=ierr)
            if (ierr.ne.0) call error(1,"error reading mol2 file created internally")
            nmols=0
            do
                read(298,*,err=299,end=299) mstring
                if (len(trim(adjustl(mstring)))>=16 .and. index(mstring,"@<TRIPOS>MOLECULE")==1) then
                    nmols=nmols+1
                    read(298,*,err=299,end=299) nb
                    read(298,*,err=299,end=299)
                    read(298,*,err=299,end=299)
                    read(298,*,err=299,end=299)
                    read(298,*,err=299,end=299)
                elseif (len(trim(adjustl(mstring)))>=12 .and. index(mstring,"@<TRIPOS>ATOM")==1) then
                    molecules(nmols)%natoms = nb
                    allocate (molecules(nmols)%bawt(nb),molecules(nmols)%bs(nb),molecules(nmols)%bx(nb),&
                    molecules(nmols)%by(nb),molecules(nmols)%bz(nb),molecules(nmols)%bq(nb),molecules(nmols)%bvdw(nb),stat=ierr)
                    do i = 1,nb
                        read(298,*,err=299,end=299) d1,molecules(nmols)%bs(i),molecules(nmols)%bx(i),molecules(nmols)%by(i),&
                        molecules(nmols)%bz(i),d2,d3,d4,molecules(nmols)%bq(i)
                        molecules(nmols)%bx(i)=molecules(nmols)%bx(i)*angau
                        molecules(nmols)%by(i)=molecules(nmols)%by(i)*angau
                        molecules(nmols)%bz(i)=molecules(nmols)%bz(i)*angau
                    end do
                endif
            enddo
299             close(298)
        endif

    else
        !Read xyz file as it contains charges
        open(unit=301, file = trim(adjustl(preprocfile)),iostat=ierr)
        if (ierr.ne.0) call error(1,"error reading xyz file")
        nmols=0
        do
            read(301,*,err=302,end=302) nb
            nmols=nmols+1
            molecules(nmols)%natoms = nb
            read(301,*,err=302,end=302)
            allocate (molecules(nmols)%bawt(nb),molecules(nmols)%bs(nb),molecules(nmols)%bx(nb),&
            molecules(nmols)%by(nb),molecules(nmols)%bz(nb),molecules(nmols)%bq(nb),molecules(nmols)%bvdw(nb),stat=ierr)
            do i = 1,nb
                read(301,*,err=302,end=302) molecules(nmols)%bs(i),molecules(nmols)%bx(i),molecules(nmols)%by(i),&
                molecules(nmols)%bz(i), molecules(nmols)%bq(i)
                molecules(nmols)%bx(i)=molecules(nmols)%bx(i)*angau
                molecules(nmols)%by(i)=molecules(nmols)%by(i)*angau
                molecules(nmols)%bz(i)=molecules(nmols)%bz(i)*angau
            end do
        enddo
302         close(301)
    endif
    ! Allocate an array for atomic mass and van der Waals radii of atoms.
    ! atmnms are atomic names defined in global module of size 103. If the sym of the
    ! atom in the new file matches with the declared atomic symbol, the
    ! local variable for at wt atwtn and at vdw atvdwn are updated.
    do i=1,nmols
        do j = 1,molecules(i)%natoms
            call ncase(molecules(i)%bs(j))   ! Case correction of Atomic symbol
            do k =1,103
                if (trim(adjustl(molecules(i)%bs(j))).ne.trim(adjustl(atmnms(k)))) cycle
                if (trim(adjustl(molecules(i)%bs(j))).eq.trim(adjustl(atmnms(k)))) then
                    molecules(i)%bawt(j) = atmwts(k)
                    molecules(i)%bvdw(j) = atmvdw(k)*0.8d0
                    exit
                end if
            end do
        enddo
    enddo

    ! If insertlocfile is given, the com of first molecule of template file is replicated at mentioned locations.
    if (trim(adjustl(preprocfile))==trim(adjustl(templatefile)) .and. len_trim(insertlocfile)/=0) then
        open(unit=304, file = trim(adjustl(insertlocfile)),iostat=ierr)
        if (ierr.ne.0) call error(1,"error reading xyz file")
        read(304,*,err=305,end=305) nb
        read(304,*,err=305,end=305)
        allocate (pisym(nb),pix(nb),piy(nb),piz(nb),stat=ierr)
        do i = 1,nb
            read(304,*,err=305,end=305) pisym(i),pix(i),piy(i),piz(i)
            call ncase(pisym(i))
            pix(i)=pix(i)*angau
            piy(i)=piy(i)*angau
            piz(i)=piz(i)*angau
        end do
305         close(304)
        tn = molecules(1)%natoms
        allocate (ts(tn),tx(tn),ty(tn),tz(tn),tq(tn),tw(tn),tv(tn),stat=ierr)
        ts = molecules(1)%bs
        tx = molecules(1)%bx ; ty = molecules(1)%by ; tz = molecules(1)%bz
        tq = molecules(1)%bq
        tw = molecules(1)%bawt
        tv = molecules(1)%bvdw
        
        ! Choose atom with maximum positive charge 
        tmaxcharg=maxloc(tq,dim=1)
        ! Get the coordinates of chosen atom
        chosenx=tx(tmaxcharg)
        choseny=ty(tmaxcharg)
        chosenz=tz(tmaxcharg)
        nmols = 0
        do i =1,nb
            ! shiftx/y/z is difference between insertloc point and chosen coordinate
            if (trim(adjustl(pisym(i)))=="X".and.nmols<50) then
                nmols = nmols + 1
                shiftx = pix(i) - chosenx 
                shifty = piy(i) - choseny 
                shiftz = piz(i) - chosenz 
!write(6,*) 'guest no ', i
!write(6,*) 'pi (bohr) = ', pix(i), piy(i), piz(i)
!write(6,*) 'shift (bohr) = ', shiftx, shifty, shiftz
                if (allocated(molecules(nmols)%bawt)) then
                   deallocate(molecules(nmols)%bawt,molecules(nmols)%bs,molecules(nmols)%bx,&
                   molecules(nmols)%by,molecules(nmols)%bz,molecules(nmols)%bq,molecules(nmols)%bvdw)
                endif
                allocate (molecules(nmols)%bawt(tn),molecules(nmols)%bs(tn),molecules(nmols)%bx(tn),&
                molecules(nmols)%by(tn),molecules(nmols)%bz(tn),molecules(nmols)%bq(tn),molecules(nmols)%bvdw(tn),stat=ierr)
                molecules(nmols)%natoms=tn
!write(6,*) 'coordinates (bohr) of guest no ', nmols
                ! Move the coordinates obtained from templatefile by shiftx/y/z amount to place molecule at insertloc points
                do j = 1,tn
                    molecules(nmols)%bs(j)  =ts(j)
                    molecules(nmols)%bx(j)  =tx(j) +shiftx
                    molecules(nmols)%by(j)  =ty(j) +shifty
                    molecules(nmols)%bz(j)  =tz(j) +shiftz
                    molecules(nmols)%bq(j)  =tq(j)
                    molecules(nmols)%bawt(j)=tw(j)
                    molecules(nmols)%bvdw(j)=tv(j)
!write(6,*) molecules(nmols)%bx(j), molecules(nmols)%by(j), molecules(nmols)%bz(j)
                enddo
            endif
        enddo
    endif

    !COM for host molecule is evaluated here.
    call comass(ncen,rcen(1,:),rcen(2,:),rcen(3,:),atwt,hcox,hcoy,hcoz)

    print *,"Number of guest molecules in input file ",nmols

    call removeclash(nmols)

    if (.not. iswindows) then
        call system ("rm -f "//mespimizerinit)
        call system ("rm -f "//mespimizervis)
        call system ("rm -f "//mespimizerfinal)
        call system ("rm -f "//mespimizerauxframes)
        call system ("rm -f "//mespimizerframes)
        call system ("rm -f "//mespimizercurrframe)
        call system ("rm -f "//mespimizerkntfrms)
    endif

    open(unit=88, file = trim(adjustl(mespimizerauxframes)),iostat=ierr)
    if (ierr .ne. 0) call error(1,"Error in creating new file for &
    geometry frames.")
    open(unit=89, file = trim(adjustl(mespimizercurrframe)),iostat=ierr)
    if (ierr .ne. 0) call error(1,"Error in creating new file for &
    geometry of current frame.")
    open(unit=90, file = trim(adjustl(mespimizerkntfrms)),iostat=ierr)
    if (ierr .ne. 0) call error(1,"Error in creating new file for &
    frames counter.")


    call writeclustergeom(mespimizerinit,nmols,0.0d0)

    kntframes = 0

    call writeframesgeom(nmols,0.0d0)

    call optimize(nmols,tssize,rssize,mespimizervis,mespimizerfinal)

!   Writes an input file for QM programs if desired
    if (lwriteqm) call writeqminput(qmsoftware,qmkeywords,mespimizerbase,nmols)

    write(90,*) kntframes
    if (iswindows) then
        open(unit=91, file = trim(adjustl(mespimizerframes)),iostat=ierr)
        if (ierr .ne. 0) call error(1,"Error in creating new file for &
            final geometry frames.")
        close(88)
        open(unit=88, file = trim(adjustl(mespimizerauxframes)),iostat=ierr)
        if (ierr .ne. 0) call error(1,"Error in creating new file for &
        geometry frames.")
        write(91,*) kntframes
        do while (1 .eq. 1)
            read(88,"(A80)",end=600,err=500) mstring
            write(91,*) trim(mstring)
            write(6,*) 'mstring = ', trim(mstring)
        end do
500     continue
        call error(1,"Error in writing file with final geometry frames.")
600     continue
        close(91)
    else
        call system("cat "//trim(adjustl(mespimizerkntfrms))//" "//trim(adjustl(mespimizerauxframes))&
                &//" > "//trim(adjustl(mespimizerframes)))
        if (.not. iswindows) then
            call system ("rm -f "//mespimizerkntfrms)
        endif
    endif

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

!!=================================
!!    Optimization control
!!=================================
        SUBROUTINE OPTIMIZE(nmols,tssize,rssize,mespimizervis,mespimizerfinal) 
        USE DAMINITIAL_T
        USE DAMBUILD_T
        IMPLICIT NONE
        CHARACTER(2), ALLOCATABLE:: bsl(:)
        REAL(KREAL), ALLOCATABLE:: bxl(:),byl(:),bzl(:),bql(:),bawtl(:),bvdwl(:)
        integer:: nmols
        REAL(KREAL) :: Enj,tssize 
        integer(kint)::ierr,i,j,k,nbatms,flag,icurr,rssize
        character(256):: mespimizervis,mespimizerfinal
        logical,external:: accept
        

        open(11,file=trim(adjustl(mespimizervis)))           
        do i = 1,nmols
            Enj=0.0d0 ! Initialize Electrostatic interaction energy 
            icurr=i
            nbatms= molecules(i)%natoms
            allocate (bsl(nbatms),bxl(nbatms),byl(nbatms),bzl(nbatms),bql(nbatms),bawtl(nbatms),bvdwl(nbatms),stat=ierr)
            flag=0
            do j=1,nbatms
		bsl(j)=molecules(i)%bs(j)
                bxl(j)=molecules(i)%bx(j)
                byl(j)=molecules(i)%by(j)
                bzl(j)=molecules(i)%bz(j)
                bql(j)=molecules(i)%bq(j)
                bawtl(j)=molecules(i)%bawt(j)
                bvdwl(j)=molecules(i)%bvdw(j)
            end do 
            if (accept(icurr,nbatms,bxl,byl,bzl,bvdwl)) &
                call ROTATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,bawtl,bvdwl,flag,rssize,enj)
            do j=1,nbatms
		molecules(i)%bs(j)=bsl(j)
                molecules(i)%bx(j)=bxl(j)
                molecules(i)%by(j)=byl(j)
                molecules(i)%bz(j)=bzl(j)
                molecules(i)%bq(j)=bql(j)
                molecules(i)%bawt(j)=bawtl(j)
                molecules(i)%bvdw(j)=bvdwl(j)
            end do 
            if (flag/=-1) then
                write(11,*) nbatms
                write(11,*) Enj  
                do j=1,nbatms           
                    write(11,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
				    molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                end do
                call writeframesgeom(nmols,Enj)
            endif

            flag = 0
            do while (flag/=-1) 
                do j=1,nbatms
					bsl(j)=molecules(i)%bs(j)
                    bxl(j)=molecules(i)%bx(j)
                    byl(j)=molecules(i)%by(j)
                    bzl(j)=molecules(i)%bz(j)
                    bql(j)=molecules(i)%bq(j)
                    bawtl(j)=molecules(i)%bawt(j)
                    bvdwl(j)=molecules(i)%bvdw(j)
                end do 
                if (accept(icurr,nbatms,bxl,byl,bzl,bvdwl)) &
                    call TRANSLATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,bawtl,bvdwl,flag,tssize,enj)
                do j=1,nbatms
                    molecules(i)%bs(j)=bsl(j)
                    molecules(i)%bx(j)=bxl(j)
                    molecules(i)%by(j)=byl(j)
                    molecules(i)%bz(j)=bzl(j)
                    molecules(i)%bq(j)=bql(j)
                    molecules(i)%bawt(j)=bawtl(j)
                    molecules(i)%bvdw(j)=bvdwl(j)
                end do 

                if (.not. (accept(icurr,nbatms,bxl,byl,bzl,bvdwl))) flag = -1
                if (flag/=-1) then
                    write(11,*) nbatms
                    write(11,*) Enj  
                    do j=1,nbatms           
                        write(11,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
					    molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                    end do
                    call writeframesgeom(nmols,Enj)
                endif
            end do

            flag=0
            do j=1,nbatms
				bsl(j)=molecules(i)%bs(j)
                bxl(j)=molecules(i)%bx(j)
                byl(j)=molecules(i)%by(j)
                bzl(j)=molecules(i)%bz(j)
                bql(j)=molecules(i)%bq(j)
                bawtl(j)=molecules(i)%bawt(j)
                bvdwl(j)=molecules(i)%bvdw(j)
            end do 
            call ROTATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,bawtl,bvdwl,flag,rssize,enj)
            do j=1,nbatms
				molecules(i)%bs(j)=bsl(j)
                molecules(i)%bx(j)=bxl(j)
                molecules(i)%by(j)=byl(j)
                molecules(i)%bz(j)=bzl(j)
                molecules(i)%bq(j)=bql(j)
                molecules(i)%bawt(j)=bawtl(j)
                molecules(i)%bvdw(j)=bvdwl(j)
            end do
            if (flag/=-1) then
                write(11,*) nbatms
                write(11,*) Enj  
                do j=1,nbatms           
                    write(11,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
				    molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                end do
                call writeframesgeom(nmols,Enj)
            endif

            deallocate(bsl,bxl,byl,bzl,bql,bawtl,bvdwl)
        enddo
!
        close(11)
        call writeclustergeom(mespimizerfinal,nmols,Enj)
        call writeframesgeom(nmols,Enj)
        close(88)
        close(89)

        end subroutine

        subroutine rotmat(ua,va,wa,ang,rt)
        USE DAMBUILD_T
        REAL(KREAL), DIMENSION(3,3) :: rt
        REAL(KREAL) :: ang,ua,va,wa
        rt(1,1) = ua**2+(1-ua**2)*cos(ang)
        rt(1,2) = ua*va*(1-cos(ang)) - wa*sin(ang)
        rt(1,3) = ua*wa*(1-cos(ang)) + va*sin(ang)
        rt(2,1) = ua*va*(1-cos(ang)) + wa*sin(ang)
        rt(2,2) = va**2+(1-va**2)*cos(ang)
        rt(2,3) = va*wa*(1-cos(ang)) - ua*sin(ang)
        rt(3,1) = ua*wa*(1-cos(ang)) - va*sin(ang)
        rt(3,2) = va*wa*(1-cos(ang)) + ua*sin(ang)
        rt(3,3) = wa**2+(1-wa**2)*cos(ang)
        return
        end subroutine rotmat
!!=================================
!!    Rotation
!!=================================
        SUBROUTINE ROTATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,batwtl,batvdwl,flag,as,enj)
        USE DAMBUILD_T
        IMPLICIT NONE 
        INTEGER::nbatms,flag,icurr
        REAL(KREAL):: enj
        REAL(KREAL), ALLOCATABLE :: enjr(:,:,:)
        REAL(KREAL), DIMENSION(nbatms):: btx,bty,btz
        REAL(KREAL), DIMENSION(nbatms):: rcom,rcox,rcoy,rcoz
        CHARACTER(2), DIMENSION(nbatms):: bsl
        REAL(KREAL), DIMENSION(nbatms):: bxl,byl,bzl,bql,batwtl,batvdwl
        REAL(KREAL), DIMENSION(3,3) :: rt
        REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
        REAL(KREAL) :: theta,phi,ang,ua,va,wa,eval
        REAL(KREAL) :: coxl,coyl,cozl
        INTEGER(KINT), ALLOCATABLE :: ex(:)
        INTEGER(KINT) :: i,j,k,l,a1,a2,a3,i1,j1,k1,as,dm,is,js,ierr,dum
        LOGICAL,external:: accept
!--------------------------------------------------
        call COMASS(nbatms,bxl,byl,bzl,batwtl,coxl,coyl,cozl)
        do i=1,nbatms
            rcox(i)=bxl(i)-coxl
            rcoy(i)=byl(i)-coyl
            rcoz(i)=bzl(i)-cozl
            rcom(i)=dsqrt(rcox(i)**2+rcoy(i)**2+rcoz(i)**2)
        end do
   
!        as=20  It is received as argument from namelist rssize
        a1=180/as+1   !Number of theta's 
        a2=360/as     !Number of phi's
        a3=360/as     !Number of rotation angle's
!       allocate (btx(nb),bty(nb),btz(nb),stat=ierr)
        allocate (enjr(a1,a2,a3),stat=ierr)
        enjr = 0.0
!       Choosing a rotation matrix of general form.
!       Rotating the axis of rotation using theta and phi
!       Rotating molecule around every axis using ang
        do i=0,180,as
            i1=i/as+1
            theta = real(i)*p_i/180.0d0
            do j=0,345,as
                j1=j/as+1
                if (i.eq.0.and.j.gt.0) exit  
                if (i.eq.180.and.j.gt.0) exit    
   
                phi = real(j)*p_i/180.0d0
                ua=sin(theta)*cos(phi)
                va=sin(theta)*sin(phi)
                wa=cos(theta)
                do k = 0,345,as
                    k1=k/as+1
                    ang = real(k)*p_i/180.0d0
                    call rotmat(ua,va,wa,ang,rt)
            
!                     write(11,*) nb
!                   write(11,*) 
                    !print *, nbatms
		    !print *, ""
                    do l=1,nbatms
                        btx(l) = coxl + (rcox(l)*rt(1,1) + rcoy(l)*rt(1,2) + rcoz(l)*rt(1,3))
                        bty(l) = coyl + (rcox(l)*rt(2,1) + rcoy(l)*rt(2,2) + rcoz(l)*rt(2,3))
                        btz(l) = cozl + (rcox(l)*rt(3,1) + rcoy(l)*rt(3,2) + rcoz(l)*rt(3,3))
                     !   print *, bsl(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
                    end do  
                    do l=1,nbatms
                        bxl(l)=btx(l) 
                        byl(l)=bty(l) 
                        bzl(l)=btz(l)
                    end do 
                    dum=2
                    if (accept(icurr,nbatms,bxl,byl,bzl,batvdwl).eqv..false.) then
                        enjr(i1,j1,k1)= +10000.0d0
                        cycle
                    end if    
                    do l=1,nbatms
                        CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
                        & dyytot,dyztot,dzztot,btx(l),bty(l),btz(l))
                        enjr(i1,j1,k1)=enjr(i1,j1,k1)+vtot*bql(l)
!                       print*,"vtot,bq(l)"
!                       print*,vtot,bq(l)
                        !write(6,*)"One move of rotation produces energy = ",enjr(i1,j1,k1)
                    end do
                end do
            end do
        end do

        dm=size(shape(enjr))
        allocate (ex(dm),stat=ierr)
        ex=minloc(enjr)
        eval=minval(enjr)
        if (eval<enj .and. (enj-eval).gt.1e-6) then
            enj=eval
!           write(6,*)ex(1),ex(2),ex(3)
            write(6,"(A,E12.5)") "Best energy of all rotations    = ",eval
            call flush(6)       ! Forces writing that is necessary to carry out amimation in the 3D viewer
            i= (ex(1)-1)*as
            theta = real(i)*p_i/180.0d0
            j=(ex(2)-1)*as
            phi = real(j)*p_i/180.0d0
            ua=sin(theta)*cos(phi)
            va=sin(theta)*sin(phi)
            wa=cos(theta)
            k=(ex(3)-1)*as
            ang = real(k)*p_i/180.0d0
            call rotmat(ua,va,wa,ang,rt)

            do l=1,nbatms
                bxl(l) = coxl + (rcox(l)*rt(1,1) + rcoy(l)*rt(1,2) + rcoz(l)*rt(1,3))
                byl(l) = coyl + (rcox(l)*rt(2,1) + rcoy(l)*rt(2,2) + rcoz(l)*rt(2,3))
                bzl(l) = cozl + (rcox(l)*rt(3,1) + rcoy(l)*rt(3,2) + rcoz(l)*rt(3,3))
            end do
        elseif ((enj-eval).le.1e-6) then
            flag=-1 
        end if 
        deallocate(ex)
        return
        end subroutine  

!!=================================
!!    Translation
!!=================================
        SUBROUTINE TRANSLATION(icurr,nbatms,bsl,bxl,byl,bzl,bql,batwtl,batvdwl,flag,tssize,enj)
        USE DAMBUILD_T
        IMPLICIT NONE 
        INTEGER::nbatms,flag,icurr
        REAL(KREAL):: enj
        REAL(KREAL), DIMENSION(3,3,3) :: enjt
        REAL(KREAL), DIMENSION(nbatms):: btx,bty,btz
        REAL(KREAL), DIMENSION(nbatms):: rcom,rcox,rcoy,rcoz
        CHARACTER(2), DIMENSION(nbatms):: bsl
        REAL(KREAL), DIMENSION(nbatms):: bxl,byl,bzl,bql,batwtl,batvdwl
        REAL(KREAL), DIMENSION(3,3) :: rt
        REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
        REAL(KREAL) :: totdrvx,totdrvy,totdrvz,norm
        REAL(KREAL) :: theta,phi,ang,ua,va,wa,eval,tssize
        REAL(KREAL) :: coxl,coyl,cozl
        INTEGER(KINT), ALLOCATABLE :: ex(:)
        INTEGER(KINT) :: i,j,k,l,a1,a2,a3,i1,j1,k1,as,dm,is,js,ierr,dum,np,nt
        LOGICAL,external:: accept
!       
        if (accept(icurr,nbatms,bxl,byl,bzl,batvdwl).eqv..false.) then
            flag=-1
        endif 
        totdrvx=0.0
        totdrvy=0.0
        totdrvz=0.0     
        do l=1,nbatms
            CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
            & dyytot,dyztot,dzztot,bxl(l),byl(l),bzl(l))
            totdrvx=totdrvx+drvx
            totdrvy=totdrvy+drvy
            totdrvz=totdrvz+drvz
        end do
        norm=max(abs(totdrvx),abs(totdrvy),abs(totdrvz))
        totdrvx=-totdrvx/norm;totdrvy=-totdrvy/norm;totdrvz=-totdrvz/norm

        !print*,"tssize"
        !print*,tssize
        !print*,"totdrvx,totdrvy,totdrvz,totdrv"
        !print*,totdrvx,totdrvy,totdrvz,sqrt(totdrvx**2+totdrvy**2+totdrvz**2)
        !print*,"tssize*totdrvx,tssize*totdrvy,tssize*totdrvz"
        !print*,tssize*totdrvx,tssize*totdrvy,tssize*totdrvz
        
        ! Move coordinates in the direction of gradient. Or check with negative gradient
        !print*,"Actual coordinate"
        do l = 1,nbatms
         !   print*,bxl(l),byl(l),bzl(l)
            bxl(l)= bxl(l)+tssize*totdrvx
            byl(l)= byl(l)+tssize*totdrvy
            bzl(l)= bzl(l)+tssize*totdrvz
!           write(11,*) bs(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
        end do
      !  print*,""
      !  print*,"Modified coordinate"
      !  do l = 1,nbatms
      !      print*,btx(l),bty(l),btz(l)
      !  end do
      !  print*,""
      !  totdrvx=0.0
      !  totdrvy=0.0
      !  totdrvz=0.0     
      !  do l=1,nbatms
      !      CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
      !      & dyytot,dyztot,dzztot,bxl(l),byl(l),bzl(l))
      !      totdrvx=totdrvx+drvx
      !      totdrvy=totdrvy+drvy
      !      totdrvz=totdrvz+drvz
      !  enddo
      !  norm=max(abs(totdrvx),abs(totdrvy),abs(totdrvz))
      !  totdrvx=totdrvx/norm;totdrvy=totdrvy/norm;totdrvz=totdrvz/norm

        !print*,"tssize"
        !print*,tssize
     !   print*,"totdrvx,totdrvy,totdrvz,totdrv"
     !   print*,totdrvx,totdrvy,totdrvz,sqrt(totdrvx**2+totdrvy**2+totdrvz**2)
        !print*,"tssize*totdrvx,tssize*totdrvy,tssize*totdrvz"
        !print*,tssize*totdrvx,tssize*totdrvy,tssize*totdrvz
        !=======================================================================
       ! enjt = 0.0d0
       ! do i= -1,1
       !     do j=-1,1
       !         do k=-1,1
       !             do l = 1,nbatms
       !                 btx(l)= bxl(l)+i*tssize
       !                 bty(l)= byl(l)+j*tssize
       !                 btz(l)= bzl(l)+k*tssize
!      !                 write(11,*) bs(l),btx(l)*auang,bty(l)*auang,btz(l)*auang
       !             end do

       !             if (accept(icurr,nbatms,btx,bty,btz,batvdwl).eqv..false.) then
       !                 enjt(i+2,j+2,k+2)= +10000.0d0
       !                 cycle
       !             end if    
       !             do l=1,nbatms
       !                 CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
       !                 & dyytot,dyztot,dzztot,btx(l),bty(l),btz(l))
       !               ! print*,drvx,drvy,drvz
       !                 enjt(i+2,j+2,k+2)=enjt(i+2,j+2,k+2)+vtot*bql(l)
       !             end do
       !         end do   
       !      end do
       ! end do
       ! dm=size(shape(enjt))
       ! allocate (ex(dm),stat=ierr)
       ! ex=minloc(enjt)
       ! eval=minval(enjt)
       ! if (eval<enj .and. (enj-eval).gt.1e-6) then
       !     ! print*, "Ok in translation"
       !     Enj=eval
       !     write(6,"(A,E12.5)") "Best energy of all translations = ",eval
       !     i=ex(1)-2
       !     j=ex(2)-2
       !     k=ex(3)-2
 
       !     print*,""
       !     print*,"Finalized coordinate"
       !     print*,""
       !     do l = 1,nbatms
       !         bxl(l)= bxl(l)+i*tssize
       !         byl(l)= byl(l)+j*tssize
       !         bzl(l)= bzl(l)+k*tssize
       !         print*,bxl(l),byl(l),bzl(l)
       !     end do
       ! elseif ((enj-eval).le.1e-6) then
       !     flag=-1
       ! end if
       ! deallocate(ex)
        return 
        end subroutine

!!=================================
!!    Center of Mass
!!=================================
        SUBROUTINE COMASS(ncatms,cx,cy,cz,catwt,coxl,coyl,cozl)
        USE DAMBUILD_T
        IMPLICIT NONE
        INTEGER:: i,ncatms
		REAL(KREAL):: coxl,coyl,cozl,tmass
        REAL(KREAL), DIMENSION(ncatms):: cx,cy,cz,catwt
       
        coxl=0.0d0; coyl=0.0d0; cozl=0.0d0; tmass=0.0d0
!write(6,*) 'template coordinates'
        do i=1,ncatms
            tmass=tmass+catwt(i)
!write(6,*) cx(i), cy(i), cz(i)
            coxl=coxl+catwt(i)*cx(i)
            coyl=coyl+catwt(i)*cy(i)
            cozl=cozl+catwt(i)*cz(i)
        end do
        coxl=coxl/tmass; coyl=coyl/tmass; cozl=cozl/tmass
        return
        END SUBROUTINE  
!    
!!==================================
!! Case correction of Atomic Symbol
!!==================================
        SUBROUTINE NCASE(CHSYM)
        IMPLICIT NONE
        CHARACTER (LEN=*) , INTENT(IN OUT) :: CHSYM
        INTEGER                            :: I,IC,NLEN

        NLEN = LEN(TRIM(ADJUSTL(CHSYM)))
        IF (NLEN.EQ.2) THEN
            IC = ICHAR(CHSYM(1:1))
            IF (IC .GE. IACHAR("a") .AND. IC .LE. IACHAR("z")) CHSYM(1:1) = CHAR(IC-32)
            IC = ICHAR(CHSYM(2:2))
            IF (IC .GE. IACHAR("A") .AND. IC .LE. IACHAR("Z")) CHSYM(2:2) = CHAR(IC+32)
        ELSE IF (NLEN.EQ.1) then
            IC = ICHAR(CHSYM(1:1))
            IF (IC .GE. IACHAR("a") .AND. IC .LE. IACHAR("z")) CHSYM(1:1) = CHAR(IC-32)
        END IF
        RETURN    
        END 
!
!!====================================
!! Van der Waals overlap check
!!====================================
        logical function accept(icurr,nbatms,bxl,byl,bzl,bvdwl)
        use dambuild_t
        implicit none
        integer(kint):: iu,ju,ku,ierr,dum,nbatms,icurr
        real(kreal), dimension(nbatms):: bxl,byl,bzl,bvdwl
        real(kreal), dimension(nbatms):: rcom,rcox,rcoy,rcoz
        real(kreal), dimension(ncen):: DIST
        
        accept=.true.
        do ju=1,ncen
            do ku=1,nbatms
                dist(ju) = sqrt((rcen(1,ju)-bxl(ku))**2+(rcen(2,ju)-byl(ku))**2+(rcen(3,ju)-bzl(ku))**2)
                if ((atvdw(ju)+bvdwl(ku)).gt.dist(ju)) then
                    accept=.false.
                    goto 98 
                end if   
            end do
98      end do       
        do iu = 1, icurr-1
            do ju=1,molecules(iu)%natoms
                do ku=1,nbatms
                    dist(ju) = sqrt((molecules(iu)%bx(ju)-bxl(ku))**2+&
                                    (molecules(iu)%by(ju)-byl(ku))**2+&
                                    (molecules(iu)%bz(ju)-bzl(ku))**2)
                    if ((molecules(iu)%bvdw(ju)+bvdwl(ku)).gt.dist(ju)) then
                        accept=.false.
                        goto 99 
                    end if   
                end do
99          end do
        end do

        return
        end 
!
!!==================================================================
!!     Write xyz file of cluster
!!==================================================================
        subroutine writeclustergeom(dfilename,nmols,ener)
        use daminitial_t
        use dambuild_t
        integer:: nmols,totalatoms,ierr
        real(kreal):: ener
        character(len=*)::dfilename

        open(unit=8, file = trim(adjustl(dfilename)),iostat=ierr)

        if (ierr.ne.0) call error(1,"Error in creating new file with &
        initial geometry of the cluster.")
        totalatoms=0
        do i =1,nmols
            totalatoms=totalatoms+ molecules(i)%natoms
        enddo   
        totalatoms=totalatoms+ncen
        write(8,*) totalatoms
        write(8,*) "Energy= ",ener
        call writecoordinates(8,nmols)
!        do i = 1,ncen
!            write(8,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
!        end do

!        do i=1,nmols
!            do j = 1,molecules(i)%natoms
!                write(8,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
!                molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
!            enddo
!        enddo
        close(8)
        end subroutine

!!==================================================================
!!     Write input for QM calculations
!!==================================================================
        subroutine writeqminput(software,keywords,ifilebase,nmols)
        implicit none
        integer:: nmols,ierr
        integer,parameter:: funit=8
        character(len=*)::ifilebase,keywords
        character(20)::software,nproc,mem,lot,basis
        character(2)::charge,multi
        ! Initialize the qm keywords in case someone passes null in keywords
        nproc="8"
        mem="1500Mb"
        lot="MP2"
        basis="6-31+G(d)"
        charge="0"
        multi="1"
        ! End Initialization

        call processkeyword(keywords,nproc,mem,lot,basis,charge,multi)

        if (software=="gaussian") open(unit=funit, file = trim(adjustl(ifilebase))//".com",iostat=ierr)
        if (ierr .ne. 0) then
            write(6,*) "error opening QM file: ", trim(adjustl(ifilebase))//".com"
            return
        endif
        call top(funit,software,ifilebase,nproc,mem,lot,basis,charge,multi)
        call writecoordinates(funit,nmols)
        call bottom(funit,software,lot,basis)
        close(funit)
        end subroutine

!!==================================================================
!!     Process the qmkeywords. It is rigid in its format, but does the job
!!==================================================================
        subroutine processkeyword(keywords,nproc,mem,lot,basis,charge,multi)
        implicit none
        integer::i,n,p
        integer,dimension(7):: comma
        character(len=*)::keywords
        character(40)::string
        character(20),intent(out)::nproc,mem,lot,basis
        character(2),intent(out)::charge,multi
        n=1
        comma(n)=1
        do i=1,len(keywords)
            if (keywords(i:i)==",") then
                n=n+1
                comma(n)=i
            endif
        enddo
        do i=1,6
            string=keywords(comma(i)+1:comma(i+1)-1)
            if ((i==1) .and. (string(:index(string,'=')-1)=="nproc")) then
                nproc=string(index(string,'=')+1:)
            elseif ((i==2) .and. string(:index(string,"=")-1)=="mem") then
                mem=string(index(string,"=")+1:)
            elseif ((i==3) .and. string(:index(string,"=")-1)=="lot") then
                lot=string(index(string,"=")+1:)
            elseif ((i==4) .and. string(:index(string,"=")-1)=="basis") then
                basis=string(index(string,"=")+1:)
            elseif ((i==5) .and. string(:index(string,"=")-1)=="charge") then
                charge=string(index(string,"=")+1:)
            elseif ((i==6) .and. string(:index(string,"=")-1)=="multi") then
                multi=string(index(string,"=")+1:)
        endif
        enddo
        end subroutine

!!==========================================================================================
!!      Compile strings that usually appear on top of geometry description in qm input file
!!==========================================================================================
        subroutine top(funit,software,ifilebase,nproc,mem,lot,basis,charge,multi)
        implicit none
        integer:: funit
        character(len=*)::software,ifilebase,nproc,mem,lot,basis,charge,multi
        if (software=="gaussian") then
            write(funit,"(a)") "%chk="//trim(ifilebase)//".chk"
            write(funit,"(a)") "%nprocShared="//trim(nproc)
            write(funit,"(a)") "%mem="//trim(mem)
            write(funit,"(a)") "#p "//trim(lot)//"/"//trim(basis)//" OPT  "// " NOSYMM "
            !"
            write(funit,"(a)") " "
            write(funit,"(a)") "Created by DAMQT"
            write(funit,"(a)") " "
            write(funit,"(a)") charge//" "//multi
        elseif (software=="gamess") then
            write(funit,"(a)") "$CONTRL SCFTYP="//lot//" CHARGE="//charge//" MULT="//multi//" RUNTYP=OPTIMIZE COORD=CART $END"
            write(funit,"(a)") "$SYSTEM TIMLIM=1 MWORDS=1 $END"
            write(funit,"(a)") "$STATPT OPTTOL=1.0E-5  $END"
            write(funit,"(a)") "$BASIS  GBASIS=STO NGAUSS=2 $END"
            write(funit,"(a)") "$DATA"
            write(funit,"(a)") "Created by DAMQT"
            write(funit,"(a)") "Cnv  2"
        elseif (software=="nwchem") then
            write(funit,"(a)") "DUMMY"
        endif
        end subroutine

!!================================================================================================
!!      Compile strings that usually appear at bottom after geometry description in qm input file
!!================================================================================================
        subroutine bottom(funit,software,lot,basis)
        implicit none
        integer:: funit
        character(len=*)::software,lot,basis
        if (software=="gaussian") then
            write(funit,"(a)") " "
        elseif (software=="gamess") then
            write(funit,"(a)") "$END"
        elseif (software=="nwchem") then
            write(funit,"(a)") "DUMMY"
        endif
        end subroutine

!!==========================
!!      Write coordinates
!!==========================
        subroutine writecoordinates(funit,nmols)
        use daminitial_t
        use dambuild_t
        integer:: nmols,funit
        do i = 1,ncen
            write(funit,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
        end do
        do i=1,nmols
            do j = 1,molecules(i)%natoms
                write(funit,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
                molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
            enddo
        enddo
        end subroutine

!
!!==================================================================
!!     Write frames geometry
!!==================================================================
        subroutine writeframesgeom(nmols,ener)
        use daminitial_t
        use dambuild_t
        integer:: nmols,totalatoms,ierr
        real(kreal):: ener
        rewind(89)
        totalatoms=0
        do i =1,nmols
            totalatoms=totalatoms+ molecules(i)%natoms
        enddo
        totalatoms=totalatoms+ncen
        write(88,*) totalatoms
        write(88,*) nmols+1, ncen, (molecules(i)%natoms, i = 1, nmols)
        write(88,*) "Energy= ",ener
        write(89,*) totalatoms
        write(89,*) "Energy= ",ener
        write(6,*) "Energy= ",ener
        do i = 1,ncen
            write(88,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
            write(89,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
            write(6,"(a2,3f15.10)") atmnms( int(zn(i))), (rcen(j,i)*0.529177249d0,j=1,3)
        end do
        do i=1,nmols
            do j = 1,molecules(i)%natoms
                write(88,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
                molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                write(89,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
                molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
                write(6,"(a2,3f15.10)") molecules(i)%bs(j),molecules(i)%bx(j)*auang,&
                molecules(i)%by(j)*auang,molecules(i)%bz(j)*auang
            enddo
        enddo
        kntframes = kntframes + 1
        call flush(89)
        call flush(6)
        end subroutine


!!==================================================================
!!     Separate molecules if they are too close to parent molecule
!!==================================================================
        subroutine removeclash(nmols)
        use daminitial_t
        use dambuild_t
        integer:: nmols,iserr
! modified by Anmol 3/6/2021
!        real(kreal):: shiftx,shifty,shiftz,dval
        real(kreal):: shiftx,shifty,shiftz,dval,norm
        integer:: tn
        real(kreal):: tcox,tcoy,tcoz
        real(kreal),allocatable ::tx(:),ty(:),tz(:),tw(:)
! end of modification
        integer(kint):: iu,ju,ku,dmdist
        real(kreal),allocatable:: DIST(:,:)
        integer,allocatable::dary(:)
        
        do iu = 1,nmols
            allocate(dist(ncen,molecules(iu)%natoms))
            do ju=1,ncen
                do ku=1,molecules(iu)%natoms
                    dist(ju,ku) = sqrt((rcen(1,ju)-molecules(iu)%bx(ku))**2+&
                                    (rcen(2,ju)-molecules(iu)%by(ku))**2+&
                                    (rcen(3,ju)-molecules(iu)%bz(ku))**2)
! modified by Anmol 3/6/2021
!                    if ((atvdw(ju)+molecules(iu)%bvdw(ku)).le.dist(ju,ku)) dist(ju,ku) = 0.0d0
                    if ((atvdw(ju)+molecules(iu)%bvdw(ku)).le.dist(ju,ku))  then
                        dist(ju,ku) = 0.0d0
                    else
                        dist(ju,ku) = atvdw(ju)+molecules(iu)%bvdw(ku)-dist(ju,ku)
                    endif
! end of modification
                end do
            enddo    
            dmdist=size(shape(dist))
            allocate (dary(dmdist),stat=ierr)
            dary=maxloc(dist)
            dval=maxval(dist)
            if (dval>0.001) then
                !print *,"===Begin======="    
                !print *,"size(shape(dist))=",dmdist
                !print *,"max location(dist)=",dary
                !print *,"Value at maxloc(dist)=",dval 
                !print *, dist
                tn = molecules(iu)%natoms
                allocate (tx(tn),ty(tn),tz(tn),tw(tn),stat=ierr)
                tx = molecules(iu)%bx ; ty = molecules(iu)%by ; tz = molecules(iu)%bz
                tw = molecules(iu)%bawt
                call comass(tn,tx,ty,tz,tw,tcox,tcoy,tcoz)
                shiftx = tcox-hcox
                shifty = tcoy-hcoy
                shiftz = tcoz-hcoz
                norm=max(abs(shiftx),abs(shifty),abs(shiftz))
                shiftx=shiftx/norm;shifty=shifty/norm;shiftz=shiftz/norm
                print*,"shiftx,shifty,shiftz"
                print*,shiftx,shifty,shiftz
                ! adding distance slightly larger than overlap
                scldval=dval*1.2 
                do ku=1,molecules(iu)%natoms
                    molecules(iu)%bx(ku)=molecules(iu)%bx(ku)+shiftx*scldval
                    molecules(iu)%by(ku)=molecules(iu)%by(ku)+shifty*scldval
                    molecules(iu)%bz(ku)=molecules(iu)%bz(ku)+shiftz*scldval
                enddo
                !print *,"===End======="    
            endif
            deallocate(dist)
        end do       

        end subroutine    
!
!
!!******************************************
!!> @details Creates guess points for MED and MESP topography.
!!      SUBROUTINE TUNE 
!!******************************************      
!!      USE DAMINITIAL_T
!!      USE DAMBUILD_T
!!      IMPLICIT NONE 
!!      REAL(KREAL), ALLOCATABLE:: enjt(:),rcom(:),rcox(:),rcoy(:),rcoz(:)
!!      REAL(KREAL), ALLOCATABLE:: btx(:),bty(:),btz(:)
!!      REAL(KREAL) :: vtot,drvx,drvy,drvz,dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot
!!      REAL(KREAL):: eval
!!      REAL(KREAL), DIMENSION(2):: frvcom,frvxcom,frvycom,frvzcom
!!      INTEGER(KINT):: i,j,k,l,n,ierr,as,np,nt,dm,is,js,nopt
!!      CHARACTER::dummy*10
!!      INTEGER(KINT):: ierr,ig,newpt
!!!      REAL(KREAL):: xmin,ymin,zmin,xmax,ymax,zmax
!!      CHARACTER(2), ALLOCATABLE:: cnsym(:)
!!      REAL(KREAL), ALLOCATABLE:: Vt(:), drvxt(:), drvyt(:), drvzt(:)
!!    REAL(KREAL), ALLOCATABLE:: xnpt(:),ynpt(:),znpt(:), atwtn(:),atvdwn(:),qnpt(:),dist(:)
!!!      REAL(KREAL) :: X1,Y1,Z1,den,vtot,drvxtot,drvytot,drvztot,
!!      REAL(KREAL) :: dxxtot, dxytot,dxztot,dyytot,dyztot,dzztot,enerji
!!!      REAL(KREAL) :: Sr,bx,by,bz,cdx,cdy,cdz
!!      REAL(KREAL):: cnox,cnoy,cnoz,tnmass
!!    REAL(KREAL),ALLOCATABLE:: frvxcom(:),frvycom(:),frvzcom(:)
!!    INTEGER(KINT):: i,j,k
!!      allocate (vt(newpt),drvxt(newpt),drvyt(newpt),drvzt(newpt),stat=ierr) 
!!      
!!    frvxcom(1)=1000.0
!!    frvycom(1)=1000.0
!!    frvzcom(1)=1000.0
!!    frvxcom(2)=0.0E+0
!!    frvycom(2)=0.0E+0
!!    frvzcom(2)=0.0E+0
!!
!!        allocate (frvxcom(10000),frvycom(10000),frvzcom(10000),stat=ierr)
!!        frvxcom=0.0E+0;frvycom=0.0E+0;frvzcom=0.0E+0;enerji=0.0E+0
!!        do i =1,10000
!!            do j=1,ncen
!!                do k=1,newpt
!!                    dist(j) = sqrt((rcen(1,j)-xnpt(k))**2+(rcen(2,j)-ynpt(k))**2+(rcen(3,j)-znpt(k))**2)
!!                    if ((atvdw(j)+atvdwn(k)).ge.dist(j)) go to 168
!!                end do
!!            end do
!!          do ig = 1,newpt
!!              CALL DAMPOT(vt(ig),drvxt(ig),drvyt(ig),drvzt(ig),&
!!              &           dxxtot,dxytot,dxztot,dyytot,dyztot,dzztot,xnpt(ig),ynpt(ig),znpt(ig))
!!          end do
!!
!!          do ig = 1,newpt
!!              frvxcom(i)=frvxcom(i)+drvxt(ig)*qnpt(ig)
!!              frvycom(i)=frvycom(i)+drvyt(ig)*qnpt(ig)
!!              frvzcom(i)=frvzcom(i)+drvzt(ig)*qnpt(ig)
!!                enerji=enerji+vt(ig)*qnpt(ig)
!!          end do
!!            write(7,*)frvxcom(i),frvycom(i),frvzcom(i),enerji
!!            if (i.eq.1) then
!!                do ig = 1,newpt
!!                    xnpt(ig)=(xnpt(ig)+frvxcom(i)*100)
!!                    ynpt(ig)=(ynpt(ig)+frvycom(i)*100)
!!                    znpt(ig)=(znpt(ig)+frvzcom(i)*100)
!!                end do
!!            end if
!!            if (abs(frvxcom(i)).le.abs(frvxcom(i-1)).and.i.gt.1.and.abs(frvxcom(i)).gt.10e-9) then
!!                do ig = 1,newpt
!!                    xnpt(ig)=(xnpt(ig)+frvxcom(i)*100)
!!                end do
!!            end if    
!!            if (abs(frvycom(i)).le.abs(frvycom(i-1)).and.i.gt.1.and.abs(frvycom(i)).gt.10e-9) then
!!                do ig = 1,newpt
!!                    ynpt(ig)=(ynpt(ig)+frvycom(i)*100)
!!                end do
!!            end if    
!!            if (abs(frvzcom(i)).le.abs(frvzcom(i-1)).and.i.gt.1.and.abs(frvzcom(i)).gt.10e-9) then
!!                do ig = 1,newpt
!!                    znpt(ig)=(znpt(ig)+frvzcom(i)*100)
!!                end do
!!            end if    
!!        end do
!!168     continue     
!!
!!  do while (abs(frvxcom(2)).le.abs(frvxcom(1)).or. &
!!         &  abs(frvycom(2)).le.abs(frvycom(1)).or. &
!!         &  abs(frvzcom(2)).le.abs(frvzcom(1)))
!!     do l=1,nb
!!          CALL DAMPOT(vtot,drvx,drvy,drvz,dxxtot,dxytot,dxztot, &
!!          & dyytot,dyztot,dzztot,bx(l),by(l),bz(l))
!!          frvxcom(2)=frvxcom(2)+drvx*bq(l)
!!          frvycom(2)=frvycom(2)+drvy*bq(l)
!!          frvzcom(2)=frvzcom(2)+drvz*bq(l)
!!          enj=enj+vtot*bq(l)
!!     end do
!!! frvcom(2)=sqrt(frvxcom(2)**2+frvycom(2)**2+frvzcom(2)**2)
!!
!!     if (.not.allocated(btx)) allocate (btx(nb),stat=ierr)
!!     if (.not.allocated(bty)) allocate (bty(nb),stat=ierr)
!!     if (.not.allocated(btz)) allocate (btz(nb),stat=ierr)
!!           if (abs(frvxcom(2)).le.abs(frvxcom(1)).and.abs(frvxcom(2)).gt.10e-9) then
!!              do i = 1,nb
!!                   btx(i)=(bx(i)+frvxcom(2)*100)
!!                end do
!!           end if    
!!           if (abs(frvycom(2)).le.abs(frvycom(1)).and.abs(frvycom(2)).gt.10e-9) then
!!                  do i = 1,nb
!!                   bty(i)=(by(i)+frvycom(2)*100)
!!                end do
!!           end if    
!!           if (abs(frvzcom(2)).le.abs(frvzcom(1)).and.abs(frvzcom(2)).gt.10e-9) then
!!              do i = 1,nb
!!                   btz(i)=(bz(i)+frvzcom(2)*100)
!!                end do
!!           end if    
!!     write(11,*) frvxcom(2),frvycom(2),frvzcom(2),frvcom(2),enj
!!
!!     frvxcom(1)=frvxcom(2)
!!     frvycom(1)=frvycom(2)
!!     frvzcom(1)=frvzcom(2)
!!
!!           do j=1,ncen
!!              do k=1,nb
!!                  dist(j) = sqrt((rcen(1,j)-btx(k))**2+(rcen(2,j)-bty(k))**2+(rcen(3,j)-btz(k))**2)
!!                  if ((atvdwn(j)+atvdwn(k)).ge.dist(j)) return
!!              end do
!!           end do
!!
!!     do i=1,nb
!!        bx(i)=btx(i)
!!        by(i)=bty(i)
!!          bz(i)=btz(i)
!!     end do
!!
!!    end do
!!      return 
!!      end subroutine
!!!

