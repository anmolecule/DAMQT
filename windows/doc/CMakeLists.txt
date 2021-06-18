# Project details
# ---------------
project(DAM2017_windows_doc)

add_subdirectory(images)

set (SCRIPTSDOC "Making_MSWindows_installer.pdf")

install (PROGRAMS ${SCRIPTSDOC} DESTINATION "../doc")
