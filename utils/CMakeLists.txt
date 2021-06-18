# Project details
# ---------------
project(DAM2016_utils)

set (PACKAGE_GENERATING_SHELLS "prepare_tarball.sh" "prepare_tarball_without_figures.sh" )

install (PROGRAMS ${PACKAGE_GENERATING_SHELLS} DESTINATION ".")
