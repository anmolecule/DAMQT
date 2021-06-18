# Project details
# ---------------
project(DAM320_doc)

set (SCRIPTSDOC "manual/DAMQT_3.2.0_manual.pdf" "quick_start/quick_start.pdf" "README_DAM320.txt")

if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local/bin")
    install (PROGRAMS ${SCRIPTSDOC} DESTINATION "../doc")
    message("Installs documentation in ${CMAKE_INSTALL_PREFIX}/../doc")
else ()
    install (PROGRAMS ${SCRIPTSDOC} DESTINATION ".")
    message("Installs documentation in ${CMAKE_INSTALL_PREFIX}")
endif ()
