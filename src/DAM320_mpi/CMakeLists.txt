# Project details
# ---------------
project(DAM320_mpi)

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
set (DAM320_mpi_PROGRAMS "DAMSTO320_mpi" "DAMGTO320_mpi" "DAMDEN320_mpi" "DAMDENGRAD320_mpi" "DAMFIELD320_mpi" 
"DAMISODEN320_mpi" "DAMISOPOT320_mpi" "DAMORB320_mpi" "DAMPOT320_mpi" "DAMPOT320_mpi_new" "DAMSGHOLE320_mpi")

set (DAM320_mpi_EXECUTABLES "" )

# Creates a library with modules
# ==============================
add_library(DAMGLOBAL_mpi_mods DAM320_GLOBAL_mpi.F90)

# Create binaries
# ===============
foreach (p ${DAM320_mpi_PROGRAMS})
	add_executable(${p}.exe "${p}.F90")
	list (APPEND DAM320_mpi_EXECUTABLES ${p}.exe )
endforeach (p)

# Links binaries to modules
# =========================
foreach (p ${DAM320_mpi_EXECUTABLES})
	target_link_libraries("${p}" DAMGLOBAL_mpi_mods ${MPE_LOG_LIBRARIES} ${MPI_Fortran_LIBRARIES})
   	set_target_properties("${p}" PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${MPE_COMPILE_FLAGS}")
	set_target_properties("${p}" PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPE_LINK_FLAGS}")
	message("Adds ${p}")
endforeach (p)

# install executables and scripts
# ===============================
install (TARGETS ${DAM320_mpi_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
