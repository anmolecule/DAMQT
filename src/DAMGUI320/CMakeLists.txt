# Project details
# ---------------
project(DAMGUI320)

add_subdirectory(include)

enable_language (C CXX)

# Instruct CMake to run moc automatically when needed.
# ===================================================

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(DAMGUI320_INCLUDE_DIR ./include)
include_directories(${DAMGUI320_INCLUDE_DIR})
message ("INCLUDE DIR FOR DAMGUI320 = " ${DAMGUI320_INCLUDE_DIR})

# here we define DAMGUI320_SRCS variable that contains a list of all .cpp files
# note that we don't need \ at the end of line
set( DAMGUI320_SRCS
   ./axes.cpp
   ./ballscyls.cpp
   ./CIsoSurface.cpp 
   ./ColorButton.cpp 
   ./contours.cpp 
   ./criticalpoints.cpp
   ./dialog.cpp  
   ./elements.cpp
   ./externals.cpp
   ./fieldlines.cpp
   ./forces.cpp
   ./geometryprocessor.cpp
   ./glmath.cpp 
   ./glWidget.cpp
   ./glWindow.cpp
   ./GlobalInfo.cpp
   ./grid.cpp
   ./IniFile.cpp
   ./isosurface.cpp
   ./lights.cpp
   ./main.cpp
   ./mainwindow.cpp 
   ./measures.cpp
   ./mespimizer.cpp
   ./molecule.cpp
   ./plotter.cpp
   ./recorder.cpp
   ./rotations.cpp
   ./screenshot.cpp
   ./Sheet.cpp 
   ./surface.cpp
   ./translations.cpp
   ./viewer2D.cpp
   ./viewport.cpp
   ./VVBuffer.cpp 
   ./widgetsubclasses.cpp
  )

# another list, this time it includes all header files that should be treated with moc (Qt's meta-object compiler)
set( DAMGUI320_MOC_HDRS
   ${DAMGUI320_INCLUDE_DIR}/axes.h
   ${DAMGUI320_INCLUDE_DIR}/ballscyls.h
   ${DAMGUI320_INCLUDE_DIR}/ColorButton.h
   ${DAMGUI320_INCLUDE_DIR}/criticalpoints.h
   ${DAMGUI320_INCLUDE_DIR}/dialog.h
   ${DAMGUI320_INCLUDE_DIR}/externals.h
   ${DAMGUI320_INCLUDE_DIR}/fieldlines.h
   ${DAMGUI320_INCLUDE_DIR}/forces.h
   ${DAMGUI320_INCLUDE_DIR}/glWidget.h
   ${DAMGUI320_INCLUDE_DIR}/glWindow.h
   ${DAMGUI320_INCLUDE_DIR}/grid.h
   ${DAMGUI320_INCLUDE_DIR}/isosurface.h
   ${DAMGUI320_INCLUDE_DIR}/lights.h
   ${DAMGUI320_INCLUDE_DIR}/mainwindow.h
   ${DAMGUI320_INCLUDE_DIR}/measures.h
   ${DAMGUI320_INCLUDE_DIR}/mespimizer.h
   ${DAMGUI320_INCLUDE_DIR}/molecule.h
   ${DAMGUI320_INCLUDE_DIR}/plotter.h
   ${DAMGUI320_INCLUDE_DIR}/recorder.h
   ${DAMGUI320_INCLUDE_DIR}/rotations.h
   ${DAMGUI320_INCLUDE_DIR}/screenshot.h
   ${DAMGUI320_INCLUDE_DIR}/Sheet.h
   ${DAMGUI320_INCLUDE_DIR}/surface.h
   ${DAMGUI320_INCLUDE_DIR}/translations.h
   ${DAMGUI320_INCLUDE_DIR}/viewer2D.h
   ${DAMGUI320_INCLUDE_DIR}/viewport.h
   ${DAMGUI320_INCLUDE_DIR}/widgetsubclasses.h
  )
  
# remaining headers
set( DAMGUI320_HDRS
   ${DAMGUI320_INCLUDE_DIR}/CIsoSurface.h
   ${DAMGUI320_INCLUDE_DIR}/contours.h
   ${DAMGUI320_INCLUDE_DIR}/elements.h
   ${DAMGUI320_INCLUDE_DIR}/geometryprocessor.h
   ${DAMGUI320_INCLUDE_DIR}/glmath.h
   ${DAMGUI320_INCLUDE_DIR}/GlobalInfo.h
   ${DAMGUI320_INCLUDE_DIR}/IniFile.h
   ${DAMGUI320_INCLUDE_DIR}/VertexNormalData.h
   ${DAMGUI320_INCLUDE_DIR}/VVBuffer.h
   )

# and finally a resource file
set( DAMGUI320_RCS
   ./damqt.qrc
   ./shaders.qrc
  )
  
# enable warnings
  add_definitions( -Wall )

# This GUI requires OpenGL library available
# ==========================================
find_package(OpenGL REQUIRED)

# This GUI requires Qt5
# =====================
# Seeks Qt5 libraries and sets all required variables

#find_package(Qt5 COMPONENTS Widgets Core Gui OpenGL PrintSupport REQUIRED)

find_package(Qt5Widgets)
# BLock for Qt5
# =============
if (Qt5Widgets_FOUND)
    set(QT_USE_QTOPENGL TRUE)
    message ("Builds GUI with Qt5")
    find_package(Qt5Core)
    find_package(Qt5Gui)
    find_package(Qt5OpenGL)
    find_package(Qt5PrintSupport)

    include_directories( ${OPENGL_INCLUDE_DIRS} )

    QT5_ADD_RESOURCES( DAMGUI320_RC_SRCS ${DAMGUI320_RCS} )
    QT5_WRAP_CPP( DAMGUI320_MOC_SRCS ${DAMGUI320_MOC_HDRS} )
    if (UNIX)
        set (DAMGUI320_EXECUTABLE "DAMQT320.exe" )
    else (UNIX)  # To prevent double .exe extension
        set (DAMGUI320_EXECUTABLE "DAMQT320" )
    endif (UNIX)
    add_executable(${DAMGUI320_EXECUTABLE} ${DAMGUI320_SRCS} ${DAMGUI320_MOC_SRCS} ${DAMGUI320_RC_SRCS})
    target_link_libraries(${DAMGUI320_EXECUTABLE} ${QT_LIBRARIES} ${OPENGL_LIBRARIES}
            Qt5::Widgets Qt5::OpenGL Qt5::PrintSupport )
else (Qt5Widgets_FOUND) # If Qt5 not found
        message (FATAL_ERROR "Qt5 not available")
endif (Qt5Widgets_FOUND)  # End of Block for Qt5Widgets

# install executables and scripts
# ===============================
install (TARGETS ${DAMGUI320_EXECUTABLE} RUNTIME DESTINATION ".")  # Default installation directory: /usr/local/bin/
