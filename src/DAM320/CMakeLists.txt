# Project details
# ---------------
project(DAM320)

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
set (DAM320_PROGRAMS "DAMSTO320" "DAMGTO320" "DAMDEN320" "DAMDENGRAD320" "DAMFIELD320" "DAMFIELDPNT320" 
"DAMFORCES320" "DAMFRAD320" "DAMISODEN320" "DAMISOPOT320" "DAMMULTROT320" "DAMORB320" "DAMPOT320" "DAMSGHOLE320")

set (DAM320_EXECUTABLES "" )

set (utilities_PROGRAMS "compareplt320" "dmat_GAUSSIAN_to_DAM" "readcnt" "readdamqt320" 
"readplt320" "sgbs2sxyz" "subtractplt320" )

set (utilities_EXECUTABLES "" )

# Creates a library with modules
# ==============================
add_library(DAMGLOBAL_mods DAM320_GLOBAL.F90)

# Create binaries
# ===============
foreach (p ${DAM320_PROGRAMS})
    if (UNIX)
	add_executable(${p}.exe "${p}.F90")
	list (APPEND DAM320_EXECUTABLES ${p}.exe )
    else (UNIX)  # To prevent double .exe extension
        add_executable(${p} "${p}.F90")
        list (APPEND DAM320_EXECUTABLES ${p} )
    endif (UNIX)
endforeach (p)

foreach (p ${utilities_PROGRAMS})
    if (UNIX)
	add_executable(${p}.exe "${p}.F90")
	list (APPEND utilities_EXECUTABLES ${p}.exe )
    else (UNIX)  # To prevent double .exe extension
        add_executable(${p} "${p}.F90")
        list (APPEND utilities_EXECUTABLES ${p} )
    endif (UNIX)
endforeach (p)

# Links binaries to modules
# =========================
foreach (p ${DAM320_EXECUTABLES})
	target_link_libraries("${p}" DAMGLOBAL_mods)
	message("Adds ${p}")
endforeach (p)

foreach (p ${utilities_EXECUTABLES})
	target_link_libraries("${p}" DAMGLOBAL_mods)
	message("Adds ${p}")
endforeach (p)

# install executables and scripts
# ===============================
install (TARGETS ${DAM320_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
install (TARGETS ${utilities_EXECUTABLES} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
