# Project details
# ---------------
project(DAMZernike320)

enable_language (C Fortran)

# Fortran compiler flags
# ======================
if (NOT CMAKE_Fortran_FLAGS_RELEASE)
	set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif (NOT CMAKE_Fortran_FLAGS_RELEASE)
if (NOT CMAKE_Fortran_FLAGS_DEBUG)
	set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
endif (NOT CMAKE_Fortran_FLAGS_DEBUG)

# Defines auxiliary containers with program names
# ===============================================
set (DAMZernike320_PROGRAMS "DAMZernike-Jacobi_GTO" "DAMZernike-Jacobi_STO" "DAMDENZJ320")

set (DAMZernike320_EXECUTABLES "" )

set (ckplm_file "Zernike-Jacobi_ckplm.F90")

SET_SOURCE_FILES_PROPERTIES(${ckplm_file} PROPERTIES COMPILE_FLAGS -O0)	# Zernike-Jacobi_ckplm.F90 compiled without optimization to prevent long delay in compilation


# Creates a library with modules
# ==============================
add_library(DAMZernike_mods   DAMZernike320_GLOBAL.F90   Zernike-Jacobi_ckplm.F90   Zernike-Jacobi_paribkb.F90  
	Zernike-Jacobi_quadrat.f  Zernike-Jacobi_stngexp.f  Zernike-Jacobi_subVABC.F90)

# Create binaries
# ===============
foreach (p ${DAMZernike320_PROGRAMS})
    if (UNIX)
        add_executable(${p}.exe "${p}.F90")
        list (APPEND DAMZernike320_EXECUTABLES ${p}.exe )
    else (UNIX)  # To prevent double .exe extension
        add_executable(${p} "${p}.F90")
        list (APPEND DAMZernike320_EXECUTABLES ${p} )
    endif (UNIX)
endforeach (p)


# Links binaries to modules
# =========================
foreach (p ${DAMZernike320_EXECUTABLES})
   target_link_libraries("${p}" DAMZernike_mods)
   message("Adds ${p}")
endforeach (p)

# install executables and scripts
# ===============================
install (TARGETS ${DAMZernike320_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
