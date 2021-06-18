# Project details
# ---------------
project(src)

enable_language (Fortran)	# This is necessary to seek MPI FORTRAN

# add subfolders
# ==============
add_subdirectory(DAM320)
add_subdirectory(TDAM320)
add_subdirectory(DAMZernike320)
add_subdirectory(DAMGUI320)
add_subdirectory(interfaces_320)
find_package(MPI)

# Checks whether MPI exists to create mpi executables
# ===================================================

if (${MPI_FORTRAN_FOUND})
	message ("\nMESSAGE: MPI available. Creates parallel programs.")
	if(${USE_MPE})
		find_package(MPE)
		if (MPE_FOUND)
			message("Uses MPE for clog2 generation")
			set (MPI_Fortran_COMPILE_FLAGS "-lfmpich -llmpe -lmpe -lpmpich" ${MPI_Fortran_COMPILE_FLAGS})
		endif (MPE_FOUND)
	endif(${USE_MPE})
    add_subdirectory(DAM320_mpi)
    add_subdirectory(TDAM320_mpi)
    add_subdirectory(DAMZernike320_mpi)
else(${MPI_FORTRAN_FOUND})
	message ("WARNING: MPI not available. Skips parallel programs compilation.")
endif (${MPI_FORTRAN_FOUND})
