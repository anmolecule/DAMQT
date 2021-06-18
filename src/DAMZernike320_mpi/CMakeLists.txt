# Project details
# ---------------
project(DAMZernike320_mpi)

enable_language (C Fortran)

find_package(MPI REQUIRED)   # This is redundant (it has been previously checked), but convenient to remember that MPI is necessary
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

# Fortran compiler flags
# ======================
if (NOT CMAKE_Fortran_FLAGS_RELEASE)
	set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif (NOT CMAKE_Fortran_FLAGS_RELEASE)
if (NOT CMAKE_Fortran_FLAGS_DEBUG)
	set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
endif (NOT CMAKE_Fortran_FLAGS_DEBUG)

# Defines auxiliary containers with program names and directory with ancillary files
# ==================================================================================
set (DAMZernike320_mpi_PROGRAMS "DAMZernike-Jacobi_STO_mpi" "DAMZernike-Jacobi_GTO_mpi" "DAMDENZJ320_mpi" )

set (DAMZernike320_mpi_EXECUTABLES "" )

set (ancillarydir "../DAMZernike320" )

# Prepares Zernike-Jacobi_ckplm.F90 to be compiled without optimization (-O0)
# =============================================================================

set (ckplm_file "${ancillarydir}/Zernike-Jacobi_ckplm.F90")

SET_SOURCE_FILES_PROPERTIES(${ckplm_file} PROPERTIES COMPILE_FLAGS -O0)	# Zernike-Jacobi_ckplm.F90 compiled without optimization to prevent long delay in compilation

# Ancillary files for libraries
# ===========================================================================

set (ancillaryfiles "${ancillarydir}/DAMZernike320_GLOBAL.F90" "${ancillarydir}/Zernike-Jacobi_ckplm.F90"  "${ancillarydir}/Zernike-Jacobi_paribkb.F90"
    "${ancillarydir}/Zernike-Jacobi_quadrat.f"  "${ancillarydir}/Zernike-Jacobi_stngexp.f")

# Creates a library with modules
# ==============================
#add_library(DAMZernikeGLOBAL_mpi_mods DAMZernike320_mpi_GLOBAL.F90 ${ancillaryfiles})
add_library(DAMZernikeGLOBAL_mpi_mods DAMZernike320_mpi_GLOBAL.F90 Zernike-Jacobi_subVABC_mpi.F90 ${ancillaryfiles})

# Create binaries
# ===============
foreach (p ${DAMZernike320_mpi_PROGRAMS})
   add_executable(${p}.exe "${p}.F90")
   list (APPEND DAMZernike320_mpi_EXECUTABLES ${p}.exe )
endforeach (p)

# Links binaries to modules
# =========================
foreach (p ${DAMZernike320_mpi_EXECUTABLES})
	target_link_libraries("${p}" DAMZernikeGLOBAL_mpi_mods ${MPE_LOG_LIBRARIES} ${MPI_Fortran_LIBRARIES})
	set_target_properties("${p}" PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPE_COMPILE_FLAGS}")
	set_target_properties("${p}" PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPE_LINK_FLAGS}")
	message("Adds ${p}")
endforeach (p)

# install executables and scripts
# ===============================
install (TARGETS ${DAMZernike320_mpi_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
# install (PROGRAMS ${SCRIPTS} DESTINATION "DAMQT320")
