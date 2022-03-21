# Project details
# ---------------
project(TDAM320)

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
#set (TDAM320_PROGRAMS "TDAMTOPO16" "TDAMOPTIMIZER")
set (TDAM320_PROGRAMS "DAMTOPOGRAPHY" "DAMMESPIMIZER")

set (TDAM320_EXECUTABLES "" )


# Creates a library with modules
# ==============================
add_library(TDAMGLOBAL_mods TDAMGLOBAL16.F90 TDAMCOM16.F90 
	TDAMPOT16.F90 TDAMDEN16.F90 lbfgsb.f linpack.f blas.f timer.f HESSIAN.F90)


#    if (UNIX)
#        set (TDAM320_EXECUTABLES "DAMTOPOGRAPHY.exe" "TDAMOPTIMIZER.exe")
#    else (UNIX)  # To prevent double .exe extension
#        set (TDAM320_EXECUTABLES "DAMTOPOGRAPHY" "TDAMOPTIMIZER"")
#    endif (UNIX)

# Create binaries
# ===============
#add_executable (${TDAM320_EXECUTABLES} "${TDAM320_PROGRAMS}.F90")

foreach (p ${TDAM320_PROGRAMS})
    if (UNIX)
        add_executable(${p}.exe "${p}.F90")
        list (APPEND TDAM320_EXECUTABLES ${p}.exe )
    else (UNIX)  # To prevent double .exe extension
        add_executable(${p} "${p}.F90")
        list (APPEND TDAM320_EXECUTABLES ${p} )
    endif (UNIX)
endforeach (p)

# Links binaries to modules
# =========================
foreach (p ${TDAM320_EXECUTABLES})
   target_link_libraries("${p}" TDAMGLOBAL_mods)
   message("Adds ${p}")
endforeach (p)

# install executables and scripts
# ===============================
install (TARGETS ${TDAM320_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
