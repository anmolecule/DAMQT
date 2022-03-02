# DAMQT

!  Copyright 2013-2022, Jaime Fernández Rico, Rafael López, Ignacio Ema,
!  Guillermo Ramírez, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
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

##  A QUICK GUIDE TO DAM320 PACKAGE

--------------------------------
                   
Authors: R. López(1), I. Ema(1), G. Ramírez(1), D. Zorrilla(2), 
A. Kumar(3), S. Yeole(4) and S. Gadre(5)

(1) Departamento de Química Física Aplicada,
Facultad de Ciencias,
Universidad Autónoma de Madrid, Spain.

(2) Departamento de Química Física,
Universidad de Cádiz, Spain.

(3) School of Pharmacy, 
University of Maryland, USA.

(4) Department of Chemistry, Bhusawal Arts, Science and 
P. O. Nahata Commerce College, Bhusawal, India

(5) Savitribai Phule Pune University, India.

Feb 2022

Contact: rafael.lopez@uam.es

## Introduction

------------

The DAM320 package consists of a set of programs designed for the partition
and representation of the molecular electron density into atomic
fragments as well as for the computation of several one-electron
functionals from this representation.

In the Deformed Atoms in Molecules method, the atomic fragments are
defined to retain as much as possible their sphericity. In the LCAO
framework, this is accomplished by assigning to each atom all the
one-center charge distributions centered on it as well as a part of its
two-center distributions. The resultant atomic fragments are expanded
in spherical harmonics times radial factors, and these radial factors
are piecewise fitted to analytical functions.

The package is prepared to work with densities computed with Slater (STO)
or Gaussian (GTO) basis sets. Geometry, basis set and the density matrix
in that basis set are required. The density matrix may be computed with any 
standard package for molecular structure calculations within the LCAO 
framework, at any calculation level (HF, CI, MC, functional density, etc).

Interfaces are included in the package for reading data from GAUSSIAN(TM)
*.fchk; MOLPRO *.out; MOPAC *.aux; NWCHEM *.nwcout; PSI4 *.fchk; 
TURBOMOLE *.coor, *.basis and *.mos; and MOLEKEL *.mkl files.

Details on the method and its implementation can be found in the
BIBLIOGRAPHY quoted at the end of this file.


## DAM320 package contents
------------------------

The package contains the source files and a sample deck of input and output files. 

The parallel versions have been prepared using MPI and tested with openmpi and mpich2  1.4.1p1

## Install Requirements

* Fortran90 compiler
* Fortran90 MPI compiler
* C++ compiler 
* python 2 or 3
* Qt library (5.9 or higher) 
* OpenGL (3.3 or above)

## Installation on Linux OS
0. ```cd /path/to/DAMQTparentDirectory```

1. If your computer does not have C++ compiler or Fortran compiler, then install the build-essential package which is a reference for all the packages needed to compile a Debian package. It generally includes the GCC/g++ compilers and libraries and some other utilities.

```
sudo apt-get install build-essential
```

2. If your computer does not have Qt5, install 

```
sudo apt-get install qt5-default
```

3. Check OpenGL version in your computer
  
```
glxinfo | grep "OpenGL version"
```

4. If you get an error in above command stating glxinfo not installed, please install mesa-utils

```
sudo apt-get install mesa-utils 
```

Check step 3 again to verify the version of OpenGL. OpenGL 3.3 or higher is required and most machines with a dedicated GPU from AMD or NVIDIA support OpenGL 3.3 and above as long as this GPU was released within the last 4-5 years. If the GPU is 5 or more years old, there is a possibility that it does not support OpenGL 3.3. In this case you will need to check the specifications for your GPU model on the manufacturer’s website. If you are limited by the version of OpenGL, DAMQT 3D viewer may not work properly. You can still make use of other functionalities.

5. Install freeglut

```
sudo apt-get install freeglut3 freeglut3-dev libglew-dev
```

6. If you want to use parallel version of package

```
sudo apt-get install openmpi-bin libopenmpi-dev
```

You will also need to add the path to bin and lib of openmpi binaries to your system path in .bashrc or .bash_profile file which is located in your home directory.

```
export PATH="path/to/mpicc/bin":$PATH
export LD_LIBRARY_PATH="path/to/mpicc/lib":$LD_LIBRARY_PATH
```

If you changed your .bashrc file, the changes will be in effect if you source it. 
```
source ~/.bashrc
```

7. For compilation of code you will require cmake or cmake-gui

```
sudo apt-get install cmake cmake-qt-gui
```

8. You may configure the package using cmake as:

```
mkdir bin
cmake -S src/ -DCMAKE_INSTALL_PREFIX=/usr/local/bin -B bin
```

9. If you use cmake-gui, follow the instructions there and make use of native gnu compilers. This should end with successful "configuration done" and "generation done".
   If configuration or generation was not successful, you need to specify the paths to libraries explicity during cmake and run cmake again till all libraries path are found.

10. After successful configuration...

```
cd bin
make -j4
sudo make install
```

11. Run the DAMQT GUI using

```
DAMQT320.exe
```

12. Various exectuables available in DAMQT can also be run from commandline by providing appropriate input file. Refer to manual for more details.

### Uninstall

To uninstall the package, type 

```
rake uninstall
```

To remove executables, type 

```
sudo make clean
```

## Samples deck
------------

A set of files is included to test the installation. The set includes input samples and the 
corresponding output for the several programs in the package.

Once you have installed the programs, run the samples and compare the output with that supplied
in the package to test that everything works.


## USER's GUIDE
------------

In all the programs, the input file consists of two parts: the first one
must always contain a NAMELIST named OPTIONS (although with different
contents for each program), in which optional variables and parameters
can be set. The namelist must appear necessarily in the top of the input
file and, optionally, it can be empty (default values will be assigned
in this case). The second part must appear after the namelist and
contains input data which are mandatory. The first record after the namelist
is always the projectname used to identify all the files in the project.

The following paragraphs briefly summarize the usage of each program in
the package.

DAMSTO320.F90  (DAMSTO320_mpi.F90)
==================================

Program for the DAM partition/expansion of a molecular density expressed in terms of Slater functions.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires information about the geometry, basis set and LCAO density matrix 
computed with a program for electron structure calculation with STO in the LCAO framework, expressed 
in the basis set.

If the calculation has been carried out with SMILES, the geometry and basis set and the density matrix
will be read from the "projectname.sgbs" and "projectname.den" generated by SMILES.
Otherwise, these data will be supplied in a file "projectname.sgbsden" or in the standard input, after 
the projectname, according to the following structure:

RECORD 1: ncen (number of centers, integer*4)

RECORDS 2 to ncen+1: rcen(1,j), rcen(2,j), rcen(3,j), zn(j), nshells (cartesian coodinates, real*8; 
     nuclear charge, real*8; and number of basis shells --i.e. sets of (n, l, exp) functions--,
     integer*4, for each center)
     
RECORDS ncen+2 and following: n, l, exp (n and l quantum numbers, integer*4; and exponent, real*8, of
    each basis set shell, in order)

remaining RECORDS: density matrix written in lower triangle form: ((dmat(i,j),j=1,i),i=1,nbasis)

Options in namelist and defaults:
    
    longoutput = .false.    If .true. a more detailed output is given (mainly for debugging)
    
    lvalence = .false.      If .true. only valence electrons are considered
    
    lzdo = .false.          If .true. ZDO approximation holds
                        
    lmaxexp = 10            Highest value of  l  in the expansion
    
    lm2c = .false.          If .true. read data from a calculation with m2c
                            
    lmultmx = 5             highest l of multipoles whose modules are stored
    
    ioptaj = 1              1: fits the total density  (default)
                            2: fits the one-center part of the density
                            3: fits the two-center part of the density
                            
    umbral = 1.d-12         Threshold for neglecting radial factors 
    
    umbralres = 1.d-12      Threshold for truncating radial factors expansions
    
    leg35 = .false.         if .true. uses a Legendre quadrature rule with 35 points for the B integrals
                            if .false.  uses a Legendre quadrature rule with 25 points for the B integrals  (default)
                            Currently, this is an obsolete option, because the integrals are computed by an alternative
                            procedure. It is kept only for testing purposes.
                        
    ipmax = 20              length of the expansion of the B integrals in subroutine frgsigma
    
    nqleg = 25              Default Legendre quadrature length for numerical integration of B integrals (Obsolescent)
    
    thresoverlap = 1.d-12   Threshold for distributions neglect
    
    u => uleg25             Default pointer initialization for Legendre quadrature abscissae (Obsolescent)
    
    w => wleg25             Default pointer initialization for Legendre quadrature weights (Obsolescent)
    
    wthreshold = 200        Threshold for B integrals (see Avk subroutine)
    
    iswindows = .false.     .true. if running on a MS-windows system

Output files
------------

DAMSTO320 writes to the standard output some information of input, the multipolar moments of the atomic fragments
and the molecular multipolar moments.

DAMSTO320 also generates two unformatted files: "projectname.damqt" and "projectname.dmqtv". 
File "projectname.damqt" contains the data of the density fit to be used by the remaining programs of the package. 
The content of this file can be read with the ancillary program "readdamqt320" also included in 
the package (readdamqt320.F90). This file requires as input data the projectname and writes the content of 
the file "projectname.damqt" as text to the standard output.
File "projectname.dmqtv" contains auxiliary integrals for computation of electrostatic potential and forces. 

A file "projectname.xyz" with the geometry in angstrom is generated for compatibility with gOpenMol and other programs.

A file "projectname.mltmod"  with the modules of the atomic multipolar moments is also created.

DAMGTO320.F90  (DAMGTO320_mpi.F90)
==================================

Program for the DAM partition/expansion of a molecular density expressed in terms of Gaussian functions.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires information about the geometry, basis set and LCAO density matrix 
computed with a program for electron structure calculation with GTO in the LCAO framework, expressed 
in the basis set. Geometry and basis set are read from the formatted file "projectname.ggbs" according to the
following structure:

RECORD 1: ncen (number of centers, integer*4)

RECORD 2 to ncen+1: rcen(1,j), rcen(2,j), rcen(3,j), zn(j) (cartesian coodinates and nuclear charge for each center, real*8)

For each center I:

    RECORD(I,1): number of contractions in center I

    For each contraction (I,J):
        
        RECORD(I,J,1): number of primitives and L quantum number in the contraction (integer*4)
        
        RECORD(I,J,2): primitive exponents
        
        RECORD(I,J,2): contraction coefficients
        
Density matrix is read from file "projectname.den" with the number of contracted basis functions followed by the
elements of density matrix written in lower triangle form::
        
RECORDS: nbasis,  ((dmat(i,j),j=1,i),i=1,nbasis)

Options in namelist and defaults:
    
    longoutput = .false.    If .true. a more detailed output is given (mainly for debugging)
    
    lvalence = .false.      If .true. only valence electrons are considered
    
    lzdo = .false.          If .true. ZDO approximation holds

    lmaxexp = 10            Highest value of  l  in the expansion
    
    lmultmx = 5             Highest l of multipoles whose modules are stored
    
    ioptaj = 1            1: fits the total density  (default)
                          2: fits the one-center part of the density
                          3: fits the two-center part of the density
                          
    umbral = 1.d-12         ! Threshold for neglecting radial factors
    
    mbralres = 1.d-12      ! Threshold for truncating radial factors expansions
    
    iswindows = .false.		! .true. if running on a MS-windows system
    

Output files
------------

DAMGTO320 writes to the standard output some information of input, the multipolar moments of the atomic fragments
and the molecular multipolar moments.

DAMGTO320 also generates two unformatted files: "projectname.damqt" and "projectname.dmqtv". 
File "projectname.damqt" contains the data of the density fit to be used by the remaining programs of the package. 
The content of this file can be read with the ancillary program "readdamqt320" also included in 
the package (readdamqt320.F90). This file requires as input data the projectname and writes the content of 
the file "projectname.damqt" as text to the standard output.
File "projectname.dmqtv" contains auxiliary integrals for computation of electrostatic potential and forces. 

A file projectname.xyz with the geometry in angstrom is generated for compatibility with gOpenMol and other programs.

A file "projectname.mltmod"  with the modules of the atomic multipolar moments is also created.


DAMDEN320.F90  (DAMDEN320_mpi.F90) 
==================================

Computes electron density from DAM partition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt generated y DAMSTO320 or DAMGTO320.

Options in namelist and defaults:
    
    longoutput = .false.     If .true. a more detailed output is given (mainly for debugging)

    langstrom = .true.       If .false. distances in bohr for .plt files
    
    ldeform = .false.        If .true. density deformations are computed (lminrep = 1)
    
    lexact  = .false.        If .true. "exact" density is tabulated from basis set and density matrix
    
    lminrep = 0              lowest "l" in the expansion of the density. It does not hold if "lexact .eq. true"
    
    lmaxrep = 10             highest "l" in the expansion of the density. It does not hold if "lexact .eq. true"
                             the expansion of the density takes lminrep <= l <= lmaxrep 
                        
    lmolec = .true.          If .true. tabulation of the whole molecular density
    
    latomics  = .false.      If .true. tabulation of atomic densities stored in files named:
                             projectname-cxx-d.plt (gOpenMol)
                             where xx refers to the center corresponding to the tabulation

    latomsel = .false.       If .true. indices of the centers for atomic tabulations will be user-supplied
                             The indices of the selected atoms must be supplied in vector "iatomsel".
                             Maximum number of centers that can be selected "mxsel" (parameter).
                         
    nsel = 0                 Number of centers for atomic tabulations

    lgradient = .false.      If .true. gradient components of the density computed
                        
    laplacian = .false.      If .true. Laplacian of density computed 
    
    lderiv2 = .false.        If .true. second derivatives of the density computed
    
    iatomsel(1) = 1          To cover the possibility of an input file with "latomsel = .true." 
                             but without assignment of "iatomsel".
    
    ldensacc  = .false.      If .true. computes and tabulates the accumulated density of the selected atoms in file 
                             projectname-frg-d.plt
    
    lboundsx = .false.       If .true. the evaluation of deformation charges is constrained to a range in x
    
    xboundinf = cero         Lower limit of that range
    
    xboundsup = cero         Upper limit of that range
    
    lboundsy = .false.       If .true. the evaluation of deformation charges is constrained to a range in y
    
    yboundinf = cero         Lower limit of that range
    
    yboundsup = cero         Upper limit of that range
    
    lboundsz = .false.       If .true. the evaluation of deformation charges is constrained to a range in z
    
    zboundinf = cero         Lower limit of that range
    
    zboundsup = cero         Upper limit of that range
    
    lgrid = .true.           If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    
    lgrid2d = .false.        If .true. computes a 2D grid. (x,y,z) are given in terms of (u,v)
    
    lpoints = .false.        If .true. computes in selected points and prints out the results in the standard output. 
                             Points must be given in cartesian coordinates. 
                             If lgrid .eq. .true., these coordinates must be placed after the grid data
    
    numrtab = 0              Number of tabulation points supplied in namelist
    
    rtab = cero              Tabulation points supplied in namelist
    
    umbrlargo = 1.d-8        Threshold for determining the short-range radius
    
    filename = ""            root file name for .plt and .pltd files
    
    iswindows = .false.      .true. if running on a MS-windows system
    
    xinf, xsup, dltx:        lower and upper bounds, and step in x coordinate for 3D grid
    
    yinf, ysup, dlty:        lower and upper bounds, and step in y coordinate for 3D grid

    zinf, zsup, dltz:        lower and upper bounds, and step in z coordinate for 3D grid
    
    x_func_uv = 'u'          x = u for 2D grids:  default: plane XY
    
    y_func_uv = 'v'          y = v for 2D grids
    
    z_func_uv = '0'          z = 0 for 2D grids
    
    uinf, usup, dltu:        lower and upper bounds, and step in u coordinate for 2D grid
    
    vinf, vsup, dltv:        lower and upper bounds, and step in v coordinate for 2D grid
    
    planeA = cero            Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    
    planeB = cero
    
    planeC = uno
    
    planecase = 1

    
Output files
------------

DAMDEN320 writes to the standard output some information of input and the values of density and its derivatives at
individual selected points.

If the 3D grid option is selected, it generates files with extensions ".plt" and ".pltd" containing the values of the several 
properties in the grid, according to the following naming convention:

projectname-d.plt:         density
projectname-d-dx.pltd:     derivative of the density with respect to x
projectname-d-dy.pltd:     derivative of the density with respect to y
projectname-d-dz.pltd:     derivative of the density with respect to z
projectname-d-dxx.pltd:    second derivative of the density with respect to x twice
projectname-d-dyy.pltd:    second derivative of the density with respect to y twice
projectname-d-dzz.pltd:    second derivative of the density with respect to z twice
projectname-d-dxy.pltd:    second derivative of the density with respect to x and y
projectname-d-dxz.pltd:    second derivative of the density with respect to x and z
projectname-d-dyz.pltd:    second derivative of the density with respect to y and z
projectname-d-lplc.plt:    Laplacian of density
projectname-cx-d.plt:      density of atomic fragment: c stands for atomic symbol and x for the center index. 

These files are compatible with gOpenMol.

If the 2D grid option is selected, it generates a file "projectname_???-d.cnt" with the values of density in the grid points.

If the selected points tabulation and second derivatives options are chosen (lpoints = .true. .and. lderiv2 = .true.) 
a text file "projectname-d.der2" with the second derivatives (dxx, dxy, dxz, dyy, dyz, dzz) is printed. 

DAMPOT320.F90  (DAMPOT320_mpi.F90)
==================================

Computes electrostatic potential from DAM partition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file "projectname.damqt" generated and DAMSTO320 or DAMGTO320.

Options in namelist and defaults:
    
    longoutput = .false.    If .true. a more detailed output is given (mainly for debugging)

    geomthr = 1.d-10        Geometry threshold: two points at a distance lower than geomthr are considered to be the coincident
    
    langstrom = .true.      If .false. distances in bohr in .plt files
    
    largo = .false.         If .true. long-range potential
    
    lexact  = .false.       If .true. "exact" potential is tabulated
    
    lgrid = .true.          If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    
    lmaxrep = 10            Highest "l" in the expansion of the potential
    
    lgradient = .false.     If .true. gradient components of the electrostatic potential computed
    
    lderiv2 = .false.       If .true. second derivatives of the electrostatic potential computed
    
    lgrid2d = .false.       If .true. computes a 2D grid. (x,y,z) are given in terms of (u,v)
   
    lpoints = .false.       If .true. computes in selected points and prints out the results in the standard output. 
                            Points must be given in cartesian coordinates. 
                            If lgrid .eq. .true., these coordinates must be placed after the grid data
                            
    lvalence = .false.      If .true. only valence electrons are considered
                        
    umbrlargo = 1.d-8       Threshold for determining the short-range radius
        
    numrtab = 0             Number of tabulation points supplied in namelist
    
    rtab = cero             Tabulation points supplied in namelist
    
    iswindows = .false.     .true. if running on a MS-windows system
    
    filename = ""           root file name for .plt and .pltd files
    
    xinf, xsup, dltx:       lower and upper bounds, and step in x coordinate for 3D grid
    
    yinf, ysup, dlty:       lower and upper bounds, and step in y coordinate for 3D grid

    zinf, zsup, dltz:       lower and upper bounds, and step in z coordinate for 3D grid
    
    x_func_uv = 'u'         x = u for 2D grids:  default: plane XY
    
    y_func_uv = 'v'         y = v for 2D grids
    
    z_func_uv = '0'         z = 0 for 2D grids
    
    uinf, usup, dltu:       lower and upper bounds, and step in u coordinate for 2D grid
    
    vinf, vsup, dltv:       lower and upper bounds, and step in v coordinate for 2D grid
        
    planeA = cero        ! Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    
    planeB = cero
    
    planeC = uno
    
    planecase = 1
    
Output files
------------

DAMPOT320 writes to the standard output some information of input and the values of electrostatic potential and its derivatives at
individual selected points.

If the 3D grid option is selected, it generates files with extension ".plt" and ".pltd" containing the values of the several 
properties in the grid, according to the following naming convention:

projectname-v.plt:         electrostatic potential
projectname-v-dx.pltd:     derivative of the electrostatic potential with respect to x
projectname-v-dy.pltd:     derivative of the electrostatic potential with respect to y
projectname-v-dz.pltd:     derivative of the electrostatic potential with respect to z
projectname-v-dxx.pltd:    second derivative of the electrostatic potential with respect to x twice
projectname-v-dyy.pltd:    second derivative of the electrostatic potential with respect to y twice
projectname-v-dzz.pltd:    second derivative of the electrostatic potential with respect to z twice
projectname-v-dxy.pltd:    second derivative of the electrostatic potential with respect to x and y
projectname-v-dxz.pltd:    second derivative of the electrostatic potential with respect to x and z
projectname-v-dyz.plt:     second derivative of the electrostatic potential with respect to y and z

These files are compatible with gOpenMol.

If the 2D grid option is selected, it generates a file "projectname_???-v.cnt" with the values of electrostatic potential in the grid points.

If the selected points tabulation and second derivatives options are chosen (lpoints = .true. .and. lderiv2 = .true.) 
a text file "projectname-v.der2" with the second derivatives (dxx, dxy, dxz, dyy, dyz, dzz) is printed.


DAMSGHOLE320.F90  (DAMSGHOLE320_mpi.F90)
========================================

Computes electrostatic potential on a density isosurface from DAM partition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt generated y DAM320 or G-DAM320.

Options in namelist and defaults:

    contourval = 1.d-3       value of density for isosurface
    
    dlthist = 1.d-3          step size for histogram
    
    filename = ""            root file name for output surface files (*.srf, *.sgh)
    
    geomthr = 1.d-5          Geometry threshold: two points at a distance lower than geomthr are considered to be the coincident
    
    gridname = ""            name of file with density grid
    
    langstrom = .true.       If .false. original grid distances in bohr
    
    lbinary = .true.         If true writes a file *.srf with the surface in binary form, otherwise writes a file *.srf_txt text mode
    
    lcolor = .false.         If true generates a file .srf with vertices positions, normals and colors
    
    lexact  = .false.		     If .true. "exact" potential is tabulated
    
    lmaxrep = 5              highest "l" in the expansion of the density and potential
    
    longoutput = .false.     If true a more detailed output is given
    
    lsghole = .true.         If true generates a file *.sgh with the vertices positions and normals and mesp values on a med isosurface 
    
    lvalence = .false.       If .true. only valence electrons are considered
    
    thrslocal = 0.9d0        Threshold for search of local maxima and minima
    
    topcolor = 0.05d0        parameter for color assignment (blue: -topcolor, red: topcolor)
    
    umbrlargo = 1.d-9        Threshold for determining the short-range radius
    
    iswindows = .false.      true if running on a MS-windows system

    
Output files
------------

DAMSGHOLE320 writes to the standard output some information of input and the number of 
electrostatic potential extrema on the density isosurface and the values and coordinates.

A file with extension ".sgh" is generated for sigma hole 3D display and another with extension 
".hst" with a histogram of MESP vs surface area on the density isosurface that can be visualized
with the 2D ploter.

A file "summary.txt" is written with a summary of statistics on MESP (see Appendix C in manual for details)

DAMFIELD320.F90  (DAMFIELD320_mpi.F90)
======================================

Computes electric field from DAM patition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file "projectname.damqt" generated y DAMSTO320 or DAMGTO320.

Optionally extra field lines can be read from an external file
whose records contain the following data:

RECORD I: icnt, x0, y0, z0    

where icnt is a center index (or 0) (integer*4), x0, y0, z0 are Cartesian coordinates of a point in the line (real*8) 

If (icnt .eq. 0 .or. icnt .gt. ncen), x0, y0, y0 is the starting point for the line. Otherwise, the line starts at the center icnt and its second 
point is x0, y0, y0. 

Options in namelist and defaults:

    basintol = 1.e-6       A point lays on the plane if its distance is lower than this
      
    ioplines3D = 1         Set of lines per nucleus for 3D plots based on icosahedron vertices, C2 axes and C3 axes or combinations of them:
                           1: vertices (12 points);   2: C3 axes (20 points);   3: C2 axes (30 points)
                           4: vertices + c3 ;   5: vertices + C2;   6: C2 + C3;   7: vertices + C2 + C3
    
    longoutput = .false.   If .true. a more detailed output is given (mainly for debugging)
    
    iswindows = .false.    .true. if running on a MS-windows system

    lmaxrep = 5            highest value of  l  in the expansion of the density for computing the field
    
    umbrlargo = 1.d-5      Threshold for determining the short-range radius
    
    usalto = 1.d-3         Threshold for convergence in stride
    
    numpnt = 2000          Maximum number of points in each field line
    
    lextralines = .false.	 If .true. reads extra lines from a file
    
    largo = .false.        If .true. computes the long-range field
    
    lplot2d = .false.      If .true. 2D plot (lines projection over a 2D surface). Not available in mpi version
        
    lvalence = .false.     If .true. only valence electrons are considered
    
    nlinpernuc = 16        Number of lines per nucleus for 2D plots
    
    filename = ""          root file name for .cam files
    
    filelines = ""         File with starting points for field lines
    
    xinf, xsup, dltx:      lower and upper bounds, and step in x coordinate for grid

    yinf, ysup, dlty:      lower and upper bounds, and step in y coordinate for grid

    zinf, zsup, dltz:      lower and upper bounds, and step in z coordinate for grid
    
    dlt0 = 1.d-2           stride length
    
    planeA = cero          Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    
    planeB = cero
    
    planeC = uno
    
    uinf = cero
    
    usup = uno
    
    uvratio = uno
    
    vinf = cero
    
    vsup = uno
    
    nlines = 0            Number of starting points for extra lines in namelist
    
    icntlines = 0         Array with indices for starting nuclei of extra lines (0 for lines starting in a point other than nuclei)
    
    rlines = cero         Array with coordinates (in au) of points defining lines

Output files
------------

DAMFIELD320 writes to the standard output some information of input and some statistics.

If 3D option is selected (lplot2d .eq. .false.), a file "projectname.cam" is generated containing the data of electric field lines. 
The lines are given as sets of coordinates x, y, z of the points in the line. 
Different lines are separated by blank lines.

If 2D option is selected (lplot2d .eq. .true.), a file "projectname.cam2D" is generated containing the data of electric field lines. 
The lines are given as sets of coordinates x, y of the points in the line. 
Different lines are separated by blank lines.

The "projectname.cam" and "projectname.cam2D" files can be used with gnuplot to plot the electic field lines. 
An example is included within the samples.

DAMFIELDPNT320.F90
====================

Computes electric field in slected points from DAM patition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project and
the Cartesian coordinates where the field will be computed.

The available options are described below. 

Besides these data, it also requires the file projectname.damqt generated y DAMSTO320 or DAMGTO320.

Options in namelist and defaults:
    
    longoutput = .false.     If .true. a more detailed output is given (mainly for debugging)
    
    lvalence = .false.      ! If .true. only valence electrons are considered

    lmaxrep = 5             highest value of  l  in the expansion of the density for computing the field
    
    umbrlargo = 1.d-5         Threshold for determining the short-range radius
    
    largo = .false.         If .true. computes the long-range field
    
    iswindows = .false.        .true. if running on a MS-windows system

Output files
------------

DAMFIELDPNT320 writes to the standard output some information of input, the values of the electric field components
at the tabulation points and some statistics.


DAMDENGRAD320.F90  (DAMDENGRAD320_mpi.F90)
==========================================

Computes density gradient from DAM partition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file "projectname.damqt" generated y DAMSTO320 or DAMGTO320.

Optionally extra field lines can be read from an external file
whose records contain the following data:

RECORD I: icnt, x0, y0, z0    

where icnt is a center index (or 0) (integer*4), x0, y0, z0 are Cartesian coordinates of a point in the line (real*8) 

If (icnt .eq. 0 .or. icnt .gt. ncen), x0, y0, y0 is the starting point for the line. Otherwise, the line starts at the center icnt and its second 
point is x0, y0, y0. 

Options in namelist and defaults:

    basintol = 1.e-6       A point lays on the plane if its distance is lower than this
      
    ioplines3D = 1         Set of lines per nucleus for 3D plots based on icosahedron vertices, C2 axes and C3 axes or combinations of them:
                           1: vertices (12 points);   2: C3 axes (20 points);   3: C2 axes (30 points)
                           4: vertices + c3 ;   5: vertices + C2;   6: C2 + C3;   7: vertices + C2 + C3
    
    longoutput = .false.   If .true. a more detailed output is given (mainly for debugging)
    
    iswindows = .false.    .true. if running on a MS-windows system

    lmaxrep = 5            highest value of  l  in the expansion of the density for computing the field
    
    umbrlargo = 1.d-5      Threshold for determining the short-range radius
    
    usalto = 1.d-3         Threshold for convergence in stride
    
    numpnt = 2000          Maximum number of points in each field line
    
    lextralines = .false.	 If .true. reads extra lines from a file
    
    lplot2d = .false.      If .true. 2D plot (lines projection over a 2D surface). Not available in mpi version
    
    nlinpernuc = 16        Number of lines per nucleus for 2D plots
    
    filename = ""          root file name for .cam files
    
    filelines = ""         File with starting points for field lines
    
    xinf, xsup, dltx:      lower and upper bounds, and step in x coordinate for grid

    yinf, ysup, dlty:      lower and upper bounds, and step in y coordinate for grid

    zinf, zsup, dltz:      lower and upper bounds, and step in z coordinate for grid
    
    dlt0 = 1.d-2           stride length
    
    planeA = cero          Default: plane for 2D plotting: XY:  A = 0, B = 0, C = 1  (z = 0)
    
    planeB = cero
    
    planeC = uno
    
    uinf = cero
    
    usup = uno
    
    uvratio = uno
    
    vinf = cero
    
    vsup = uno
    
    nlines = 0            Number of starting points for extra lines in namelist
    
    icntlines = 0         Array with indices for starting nuclei of extra lines (0 for lines starting in a point other than nuclei)
    
    rlines = cero         Array with coordinates (in au) of points defining lines

Output files
------------

DAMDENGRAD320 writes to the standard output some information of input and some statistics.

If 3D option is selected (lplot2d .eq. .false.), a file "projectname.dengr" is generated containing the data of density gradient lines. 
The lines are given as sets of coordinates x, y, z of the points in the line. 
Different lines are separated by blank lines.

If 2D option is selected (lplot2d .eq. .true.), a file "projectname.dengr2D" is generated containing the data of density gradient lines. 
The lines are given as sets of coordinates x, y of the points in the line. 
Different lines are separated by blank lines.

The "projectname.dengr" and "projectname.dengr2D" files can be used with gnuplot to plot the electic field lines. 
An example is included within the samples.


DAMFORCES320.F90
==================

Computes Hellmann-Feynman forces on nuclei from DAM patition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt generated y DAM320 or G-DAM320.

Options in namelist and defaults:
    
    longoutput = .false.    If .true. a more detailed output is given (mainly for debugging)
    
    latomsel = .false.      If .true. indices of the centers for atomic tabulations will be user-supplied
                            The indices of the selected atoms must be supplied in vector "iatomsel".
                            
    lvalence = .false.      If .true. only valence electrons are considered
                            
    lmaxrep = 10            highest "l" for the computation of forces
    
    ncntab = 1              Number of atoms selected for tabulation (dummy, only used by DAMQT GUI)
    
    iatomsel(1) = 1         To cover the possibility of an input file with "latomsel = .true."
    iatomsel(2:mxsel) = 0   but without assignment of "iatomsel". 
    
    umbrlargo = 1.d-8       Long-range threshold
                        
    iatomsel(1) = 1         If "latomsel = .true." selected atoms for detailed forces decomposition
    
    iswindows = .false.     .true. if running on a MS-windows system
    
    filename = ""           root file name for .fre, .fri, .frt, .fcf, .fnc files

Output files
------------

DAMFRAD320 writes to the standard output some information of input and the tabulation of the selected radial factor
and its first and second derivatives.

It also generates a file "projectname.forces" which contains the components of the Hellmann-Feynman forces over nuclei in the following order: 

internal forces on nuclei (forces due by the electron cloud of the atom to which the nucleus belongs)
external forces on nuclei (forces due to the electron cloud of the remaining atoms)
total forces (sum of the two previous ones)
conformational forces
non-conformational forces (spurious forces)



DAMFRAD320.F90
================

Tabulates radial factors of DAM patition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt generated y DAM320 or G-DAM320.

Options in namelist and defaults:
    
    longoutput = .false.    If .true. a more detailed output is given (mainly for debugging)
    
    lderiv = .false.        If .true. tabulates the first derivative of the radial factors
    
    lderiv2 = .false.       If .true. tabulates the second derivative of the radial factors

    rini = cero             Lower bound of interval for tabulation
    
    rfin = dos              Upper bound of interval for tabulation
    
    dltr = ri(10)           Step for tabulation
    
    ncntab = 1              Number of atoms selected for tabulation
    
    iatomsel(1) = 1         Atoms selected for tabulation
    
    iatomsel(2:mxsel) = 0    
    
    ltab = 0                l  value for tabulation
    
    mtab = 0                m  value for tabulation
    
    lrlist = .false.        If true, list of r values for tabulation supplied in rlist
    
    nlist = 0               Number of elements in rlist
    
    iswindows = .false.     .true. if running on a MS-windows system
    
    filename = ""           root file name

Output files
------------

DAMFRAD320 writes to the standard output some information of input and the tabulation of the selected radial factor
and its first and second derivatives.
It also writes a file ".frad" with the tabulation data.

DAMMULTROT320.F90
===================

Computes oriented multipolar moments of atomic fragments from DAM patition/expansion.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project, 
and the indices of the three centers whose multipolar moments will be rotated.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt generated y DAM320 or G-DAM320.

Options in namelist and defaults:

    iswindows = .false.   .true. if running on a MS-windows system
    
    i1 = 1                index of first center defining the axes orientation
    
    i2 = 2                index of second center defining the axes orientation
    
    i3 = 3                index of third center defining the axes orientation
    
    lmin = 0              lowest order of multipolar moments to be rotated
    
    lmin = 0              highest order of multipolar moments to be rotated
    
    ncntab =              number of centers for tabulations (apart from the centers used to define the frame)
    
    icntab(:) = 0         indices of centers to be tabulated
    
    rip = 0               vector defining the X' axis of the rotated frame (should be in the plane defined by the
                          three centers, otherwise its projection in that plane is taken)
                          
    filename = ""		      root file name

Output files
------------

DAMMULTROT320 writes to the standard output some information of input and the tabulation of the selected multiplar moments
of the three selected centers referred to normalized spherical harmonics in the frame with:
    Z' axis perpendicular to the centers
    Y' axis bisector of the angle formed by the three centers (second center in the angle vertex)
    X' axis orthogonal to Z' and Y' (right-handed system)

    
DAMZernike-Jacobi_STO.F90   (DAMZernike-Jacobi_STO_mpi.F90)
=====================================================================

Computes Zernike 3D and Jacobi moments of a STO molecular density.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project, 
and the indices of the three centers whose multipolar moments will be rotated.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt file.

Options in namelist and defaults:

    kexpansion = 10            Highest k in expansion functions
    
    kmax0 = 200                length of expansion of radial functions for STO translation
    
    jmax = 100                 length of the series used for comuting  the starting BesselI functions
    
    lechelon = .false.         if true number of functions per l equal to max(lexpansion+1,kexpansion)-l
    
    lexpansion = 20            Highest l in expansion functions
    
    ljacobi = .false.          if .true. uses Jacobi P(0,2+2l) polynomials as radial functions instead of Zernike 3D
    
    lmthfile = .false.         if .true. generates a projectname.mth file with data about the moments computation
    
    longoutput = .false.       If .true. a more detailed output is given
    
    lrstarrel = .false.        If .true. ball radius equal to distance of farthest atom plus rstar
    
    lvalence = .false.         If .true. only valence electrons are considered
    
    lzdo = .false.             If .true. ZDO approximation holds
    
    nquadpoints = 128          Number of quadrature points
    
    rstar = 10.d0              Ball radius for expansion
    
    thresmult = 1.d-10         Threshold for printing multipole moments
    
    thresoverlap = 1.d-12      Threshold for distributions neglect
    
    iswindows = .false.        .true. if running on a MS-windows system

Output files
------------

DAMZernike-Jacobi_GTO(STO) writes to the standard output some information of input and the one-center expansion
of GTO(STO) density into Canterkis-Zernike or Jacobi-Zernike functions. In particular, rotationally invariant fingerprints
are quoted.

It also writes a file ".zernike" (".jacobi") with the expansion of the density in Zernike (Jacobi) functions.
    
DAMZernike-Jacobi_GTO.F90   (DAMZernike-Jacobi_GTO_mpi.F90)
=====================================================================

Computes Zernike 3D and Jacobi moments of a GTO molecular density.

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project, 
and the indices of the three centers whose multipolar moments will be rotated.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt file.

Options in namelist and defaults:

    kexpansion = 10            Highest k in expansion functions
    
    lechelon = .false.         if true number of functions per l equal to max(lexpansion+1,kexpansion)-l
    
    lexpansion = 20            Highest l in expansion functions
    
    ljacobi = .false.          if .true. uses Jacobi P(0,2+2l) polynomials as radial functions instead of Canterakis-Zernike
    
    lmthfile = .false.         if .true. generates a projectname.mth file with data about the moments computation
    
    longoutput = .false.       If .true. a more detailed output is given
    
    lrstarrel = .false.        If .true. ball radius equal to distance of farthest atom plus rstar
    
    lvalence = .false.         If .true. only valence electrons are considered
    
    lzdo = .false.             If .true. ZDO approximation holds
    
    nquadpoints = 128          Number of quadrature points
    
    rstar = 10.d0              Ball radius for expansion
    
    thresmult = 1.d-10         Threshold for printing multipole moments
    
    thresoverlap = 1.d-12      Threshold for distributions neglect
    
    iswindows = .false.        .true. if running on a MS-windows system

Output files
------------

DAMZernike-Jacobi_GTO(STO) writes to the standard output some information of input and the one-center expansion
of GTO(STO) density into Canterkis-Zernike or Jacobi-Zernike functions. In particular, rotationally invariant fingerprints
are quoted.

It also writes a file ".zernike" (".jacobi") with the expansion of the density in Zernike (Jacobi) functions.

    
DAMDENZJ320.F90
===============

Tabulates density from expansion in Canterakis-Zernike or Jacobi functions

Input files
-----------

In the standard input, a namelist (OPTIONS) containing optional parameters will be read first, 
followed by the projectname wich will be used to identify all the files corresponding to the project, 
and the indices of the three centers whose multipolar moments will be rotated.
The available options are described below. 

Besides these data, it also requires the file projectname.damqt file.

Options in namelist and defaults:

    langstrom = .true.         If .false. distances in bohr
    
    lechelon = .false.         if true number of functions per l equal to max(lexpansion+1,nexpansion)-l
    
    lindividk = .false.        If true projection functions with given values of k index and all l and m compatible indices are used
    
    lindividl = .false.        If true projection functions with given values of l index and all k and m compatible indices are used
    
    lindividlk = .false.       If true projection functions with given values of k, l indices and all m compatible indices are used
    
    lindividlkm = .false.      If true only projection functions with given values of k, l, m indices are used
    
    ljacobi = .false.          If .true. expansion in Jacobi functions
    
    lgrid = .true.             If .true. computes and tabulates on a grid. The results are stored in an external file *.plt
    
    lgradient = .false.        If .true. gradient components of the density computed and, if lgrid = .true., tabulated in files
                                  projectname-d-dx.pltd, projectname-d-dy.pltd, projectname-d-dz.pltd
                              
    lgrid2d = .false.          If .true. computes a 2D grid. (x,y,z) are given in terms of (u,v)
    
    lpoints = .false.          If .true. computes in selected points and prints out the results in the standard output.
                               Points must be given in cartesian coordinates.
                               If lgrid .eq. .true., these coordinates must be placed after the grid data
                              
    kmaxrep = 10               Highest k in expansion
    
    lminrep = 0                Lowest l in expansion
    
    lmaxrep = 10               Highest l in expansion
    
    numrtab = 0                Number of tabulation points supplied in namelist
    
    filename = ""              root file name for .plt and .pltd files
    
    fileZJname = ""            .zernike or .jacobi files
    
    iswindows = .false.        .true. if running on a MS-windows system
    
    xinf, xsup, dltx:          lower and upper bounds, and step in x coordinate for grid

    yinf, ysup, dlty:          lower and upper bounds, and step in y coordinate for grid

    zinf, zsup, dltz:          lower and upper bounds, and step in z coordinate for grid
    
    x_func_uv = 'u'         x = u for 2D grids:  default: plane XY
    
    y_func_uv = 'v'         y = v for 2D grids
    
    z_func_uv = '0'         z = 0 for 2D grids
    
    uinf, usup, dltu:       lower and upper bounds, and step in u coordinate for 2D grid
    
    vinf, vsup, dltv:       lower and upper bounds, and step in v coordinate for 2D grid

Output files
------------

DAMDENZJ320 writes to the standard output some information of input and the tabulation of one-center expansion
of GTO(STO) density into Canterkis-Zernike or Jacobi functions. 

If the grid option is selected, generates files with extensions ".plt" and ".pltd" containing the values of the several 
properties in the grid, according to the following naming convention:

projectname-type-d.plt:         density
projectname-type-d-dx.pltd:     derivative of the density with respect to x
projectname-type-d-dy.pltd:     derivative of the density with respect to y
projectname-type-d-dz.pltd:     derivative of the density with respect to z

where type can be "zernike" or "jacobi"

readdamqt320.F90
==================

Reads the content of an unformatted file "projectname.damqt" and dumps it as text to the standard output.

To run it, just type

readdamqt320.exe

and you will be prompted for the name of the ".damqt" file to be read. Type it (with or without the ".damqt" extension) and
press enter.

To get the output into a file, use the  >  character. For instance:

readdamqt320.exe > projectname.damqt_txt


readcnt.F90
===========

Reads the content of a binary file with 2D grids "fname.cnt" and dumps it as text to a file "fname.cnt_txt", where 
fname stands for the corresponding file name. I creates also a "fname.cnt_gnu_txt" that can be loaded with gnuplot.


To run it, just type

readcnt.exe

readplt320.F90
================

Reads the content of a binary grid file "fname.plt" ("fname.pltd") and prints it to a text file named "fname.plt_txt"
("fname.pltd_txt")

To run it, just type

readplt320.exe

and you will be prompted for the name of the ".plt" file to be read. Type it (with or without the .plt extension) and
press enter.

Output goes to a file with the same name as the original ".plt" file appended with "txt"

subtractplt320.F90
====================

Subtracts the content of two ".plt" files and generates another ".plt" file with the result. 
It also gives some statistics of the comparison of the input files. 

To run, just type compareplt320.exe

Then, supply the names of the ".plt" files to be compared. 

Grids of two files must be equal. Otherwise, it will issue an error message.

Output goes to standard output.

compareplt320.F90
===================

Compares the content of two ".plt" files. Gives some statistics of the comparison. Intended for accuracy tests.

To run, just type 

compareplt320.exe

Then, supply the names of the ".plt" files to be compared. 

Grids of two files must be equal. Otherwise, it will issue an error message.

Output goes to standard output.

dmat_GAUSSIAN_to_DAM.F90
========================

Transforms the density matrix from GAUSSIAN order of the angular functions to DAM order.

input:
    nshells: integer with the number of shells in the basis set (each shell corresponds to a set of basis functions differing only
            in the "m" quantum number (i.e., l = 1 implies three P functions; l = 2, five D functions; and so forth
    lshell: array of dimension nshells with the l quantum number of each shell (in the same order as occurring in the basis set)
    dmatGAUSS: lower triangle of the density matrix as given by GAUSSIAN
output:
    dmatDAM:    lower triangle of density matrix in DAM order (canonical order of spherical hemonics)

INTERFACES
----------

GAUSS_interface.cpp
===================

Interface for generating DAM input files with geometry and basis set (.ggbs) 
and density matrix (.den) from  from Gaussian(C) .fchk file.

IMPORTANT: Calculations must be carried out with spherical basis sets (Not Cartesian). 
This is achieved by including in the Gaussian input file the option 5D,7F.

When called without arguments (by just typing GAUSS_interface) the program asks for
the name of the ".fchk" file to be processed (without the expansion .fchk) and the
SCF density matrix will be read.

When called with one argument, the argument must be the name of the .fchk file to be 
processed. The SCF density matrix will be read.

When called with two arguments, the first one must be the name of the ".fchk" file to be 
processed, and the second will be used to select the density matrix to be used, accrding
to the following choices:

0 :    checks for the density matrices available in the .fchk file. Prints the list of
    existing density files to outputfile (named GAUSS-'filename'.out, where 'filename'
    stands for the name of the .fchk file.
    
1 :    takes the SCF density matrix
2 :    takes the SCF spin-density matrix
3 :    takes the CI density matrix
4 :    takes the CI spin density matrix
5 :    takes the MP2 electron density
6 :    takes the MP2 spin density

MOLPRO_interface.cpp
====================

Interface for generating DAM input files with geometry and basis set (.ggbs) 
and density matrix (.den) from MOLPRO output *.out files. 

The interface must be run without arguments. It will ask for the name of the .out file. 
DAM files .ggbs and .den will be built with the same name as that of the .out file.

In order to generate a suitable input for DAM, basis set and density matrix must be 
written in MOLPRO's output file. The following code in MOLPRO's input file does the job:

...
gprint,basis
...
{matrop
load,d,den
print,d}
...

WARNING! Beware that the bug  4126 has been fixed in your MOLPRO code. 
Otherwise, patch MOLPRO's file src/argos/arinp.F file with 
the file "arinp.F.diff" included in the package. This bug may yield wrong results when
DAM320 is run with MOLPRO symmetry adapted functions.

MOLEKEL_interface.cpp
=====================

Interface for generating DAM input files with geometry and basis set (.ggbs) 
and density matrix (.den) from  from MOLEKEL .mkl file.

The interface must be run without arguments. It will ask for the name of the .mkl file. 
DAM files .ggbs and .den will be built with the same name as that of the .mkl file. 

TURBOMOLE_interface.cpp
=======================

Interface for generating DAM input files with geometry and basis set (.ggbs) 
and density matrix (.den) from  from TURBOMOLE basis, coords and mos (or alpha and beta)
files.

Molecular orbitals in file mos must be expressed in the original basis set, not in 
symmetry-adapted functions. 


Mopac_aux_interface.cpp
=======================

Interface for generating DAM input files with geometry and basis set (".ggbs") 
and density matrix (".den") from MOPAC ".aux" file.

NWChem_interface.cpp
====================

To make the interface accesible by just clicking on the outputfile, it is necessary 
to set the output file extension as ".nwcout".

For the interface to work, the mov2asc executable must be available in directory 
        $NWCHEM TOP/contrib/
where $NWCHEM TOP stands for the NWCHEM home directory.



BIBLIOGRAPHY
------------

(1) Analysis of the molecular density
    Fernández Rico, J.; López, R.; Ramírez, G. J Chem Phys
    1999, 110, 4213-4220.

(2) Density and Binding Forces in Diatomics
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G. J Chem
    Phys 2002, 116, 1788-1799.

(3) Analysis of the molecular density: STO densities
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G. J Chem
    Phys 2002, 117, 533-540.

(4) Density and binding forces: Rotational barrier of ethane
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G. J Chem
    Phys 2003, 119, 12251-12256.

(5) Analytical method for the representation of atoms-in-molecules
    densities
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.;
    Ludeña, E. V. J Comput Chem 2004, 25, 1355-1363.

(6) Electrostatic potentials and fields from density expansions of deformed
    atoms in molecules
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.
    J Comput Chem 2004, 25, 1347-1354.
    
(7) Chemical notions from the Electron density
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.
    J Chem Theory Comput 2005, 1, 1083-1095
    
(8) Deformed atoms in molecules: analytical representation of atomic densities
    for Gaussian type orbitals
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.
    J Mol Struct (Theochem) 2005, 727, 115-121
    
(9) Chemical forces in terms of the electron density
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.
    Theor Chem Account 2007, 118, 709-721
    
(10)DAMQT: A package for the analysis of electron density in molecules 
    López, R.; Fernández Rico, J.; Ramírez, G.; Ema, I., Zorrilla, D.
    Comput. Phys. Commun. 2009, 180, 1654–1660

(11)Translation of real solid spherical harmonics
    Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.
    International J Quantum Chem, 2013, 113, 1544–1548
    
(12)Improved partition-expansion of two-center distributions involving Slater functions
    López, R.; Ramírez, G.; Ema, I.; Fernández Rico, J.
    J Comput Chem, 2013, 34, 1800–1809
    
(13)Multipole moments from the partition-expansion method
    López, R.; Ramírez, G.; Fernández Rico, J.; Ema, I.
    Theoret. Chem. Acc., 2013, 132, 1406
    
(14)DAMQT 2.0: A new version of the DAMQT package for the analysis of electron
    density in molecules
    López, R.; Fernández Rico, J.; Ramírez, G.; Ema, I., Zorrilla, D.
    Comput. Phys. Commun. 2015, 192, 289-294
    
(15)DAMQT 2.1.0: A New Version of the DAMQT Package Enabled with the 
    Topographical Analysis of Electron Density and Electrostatic Potential in Molecules
    Kumar, A.; Yeole, S.; Gadre, S.; López, R.; Fernández Rico, J.; Ramírez, G.; Ema, I., Zorrilla, D.
    J Comput Chem, 2015, 36, 2350-2359
    
(16)Topology of molecular electron density and electrostatic potential with DAMQT
    López, R.; Fernández Rico, J.; Ramírez, G.; Ema, I., Zorrilla, D.; Kumar, A.; Yeole, S.; Gadre, S.
    Comput. Phys. Commun. 2017, 214, 2207-215
    
(17)Efficient Algorithm for Expanding Theoretical Electron Densities in Canterakis-
    Zernike Functions
    Urquiza-Carvalho, G.; Rocha, G.; López, R
    J Comput Chem, 2018, DOI 10.1002/jcc.25376
    
(18)Efficient Evaluation of Molecular Electrostatic Potential in Large Systems
    R. López, F. Martínez, I. Ema, J.M. García de la Vega, G. Ramírez
    Computation, 7, 64, 2019, doi: 10.3390/computation7040064


(19)Molecular fingerprints based on Jacobi expansions of electron densities
    Rafael López, Frank Martínez, José Manuel García de la Vega
    Theoretical Chemistry Accounts (2021) 140:18, doi: 10.1007/s00214-020-02708-7


CREDITS
-------

gOpenMol: Copyright 1997-2005 Leif Laaksonen (gopenmol@csc.fi)
http://www.csc.fi/gopenmol/index.phtml

gnuplot: Copyright(C) 1986 - 1993, 1998 - 2002
Thomas Williams, Colin Kelley and many others

Qt: is property of The Qt Company (Espoo, Finland) www.qt.io

Dr. George Benthien, Fortran character string utilities and math evaluation module
(https://gbenthien.net/strings/index.html)

Raghavendra Chandrashekara, implementation of marching cubes algorithm based on source code
provided by Paul Bourke and Cory Gene Bloyd.



