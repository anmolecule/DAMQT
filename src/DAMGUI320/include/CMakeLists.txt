# Project details
# ---------------
project(DAM311_include)

# another list, this time it includes all header files that should be treated with moc (Qt's meta-object compiler)
set( DAMGUI311_MOC_HDRS
    ./ColorButton.h
    ./criticalpoints.h
    ./dialog.h
    ./fieldlines.h
    ./forces.h
    ./glWidget.h
    ./glWindow.h
    ./grid.h
    ./mainwindow.h
    ./measures.h
    ./molecule.h
    ./plotter.h
    ./Sheet.h
    ./isosurface.h
    ./surface.h
    ./viewer2D.h
    ./widgetsubclasses.h
  )

# remaining headers
set( DAMGUI311_HDRS
    ./CIsoSurface.h
    ./contours.h
    ./elements.h
    ./geometryprocessor.h
    ./glmateriales.h
    ./glmath.h
    ./GlobalInfo.h
    ./IniFile.h
    ./VVBuffer.h
    ./VertexNormalData.h
   )
