# CMake project file for DAMQT

cmake_minimum_required(VERSION 3.1.0)

# about this project
# ==================
project (DAMQT)

# preparation for distribution package generation with CPack
#   (For the moment, I have not been able to make it work properly:
#    symbolic links are not replaced by actual files and it
#    stores too much stuff)
# ==========================================================
# include(InstallRequiredSystemLibraries)
# set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DAMQT: a program for the analysis of molecular electron density and related properties")
# set(CPACK_PACKAGE_VENDOR "Freeware")
# SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/doc/README_DAM2016.txt")
# SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license/COPYING")
# set(CPACK_PACKAGE_VERSION_MAJOR "2")
# set(CPACK_PACKAGE_VERSION_MINOR "1")
# set(CPACK_PACKAGE_VERSION_PATCH "0")
# set(CPACK_SOURCE_GENERATOR "TBZ2")
# set(CPACK_SOURCE_PACKAGE_FILE_NAME 
# "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
# set(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
# set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
# if(WIN32 AND NOT UNIX)
#   # There is a bug in NSI that does not handle full unix paths properly. Make
#   # sure there is at least one set of four (4) backslashes.
#   set(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/utils\\\\damqt.png")
#   set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\DAMQT2016")
#   set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} DAMQT Project")
# #   set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.my-project-home-page.org")
# #   set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.my-personal-home-page.com")
#   set(CPACK_NSIS_CONTACT "rafael.lopez@uam.es")
#   set(CPACK_NSIS_MODIFY_PATH ON)
# else(WIN32 AND NOT UNIX)
#   set(CPACK_STRIP_FILES "bin/DAMQT2016")
#   set(CPACK_SOURCE_STRIP_FILES "")
# endif(WIN32 AND NOT UNIX)
# set(CPACK_PACKAGE_EXECUTABLES "DAMQT2016" "DAM_2016.exe")
# 
# file (GLOB_RECURSE cmakefiles RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeFiles/*"
#   RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/*/*/CMakeFiles/*" 
#  )
# foreach(cmakefile ${cmakefiles})
#    message ("fichero: " ${cmakefile})
# endforeach(cmakefile)
# set(CPACK_SOURCE_IGNORE_FILES
#   "/build/"
#   "/bin/"
#   "/debug/"
#   "/.bzr/"
#   "~$"
#   cmakefiles
#   "${CPACK_SOURCE_IGNORE_FILES}"
#   )
# include(CPack)
# add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
#
# End of package preparation with CPack
# =====================================


# make sure that the default is a RELEASE
# =======================================
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Set c++ compiler to standard 11
set (CMAKE_CXX_STANDARD 11)

# default installation
# ====================
get_filename_component (default_prefix "/usr/local/bin/" ABSOLUTE)
message ("CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
if (NOT CMAKE_INSTALL_PREFIX)
    set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
          "Choose the installation directory; by default it installs in the /usr/local/bin directory. Root privileges may be required"
          FORCE)
endif(NOT CMAKE_INSTALL_PREFIX)
# CACHE variable to control if compilation with MPE is wanted for log files generation to analyze MPI performance 
# ===============================================================================================================
set(USE_MPE OFF CACHE FILEPATH "Set to ON for compiling with MPE if available") 
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils/" ${CMAKE_MODULE_PATH})

# add subfolders
# ==============
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(license)
add_subdirectory(samples)
add_subdirectory(utils)
add_subdirectory(windows)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
