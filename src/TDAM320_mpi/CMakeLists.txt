# Project details
# ---------------
project(TDAM320_mpi)

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

# Defines auxiliary containers with program names
# ===============================================
set (TDAM320_mpi_PROGRAMS "TDAMTOPO16_mpi")

set (TDAM320_mpi_EXECUTABLES "DAMTOPOGRAPHY_mpi.exe")

# Creates a library with modules
# ==============================
add_library(TDAMGLOBAL_mpi_mods TDAMGLOBAL16_mpi.F90 TDAMCOM16_mpi.F90
    TDAMPOT16_mpi.F90 TDAMDEN16_mpi.F90 lbfgsb.f linpack.f blas.f timer.f HESSIAN.F90)

# Create binaries
# ===============
add_executable (${TDAM320_mpi_EXECUTABLES} "${TDAM320_mpi_PROGRAMS}.F90")

# Links binaries to modules
# =========================
foreach (p ${TDAM320_mpi_EXECUTABLES})
    target_link_libraries("${p}" TDAMGLOBAL_mpi_mods ${MPE_LOG_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    set_target_properties("${p}" PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPE_COMPILE_FLAGS}")
    set_target_properties("${p}" PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPE_LINK_FLAGS}")
    message("Adds ${p}")
endforeach (p)

# install executables and scripts
# ===============================
install (TARGETS ${TDAM320_mpi_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
