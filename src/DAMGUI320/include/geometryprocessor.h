//  Copyright 2008-2021, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
//  Guillermo Ramirez, David Zorrilla, Anmol Kumar, Sachin D. Yeole, Shridhar R. Gadre
//
//  This file is part of DAMQT.
//
//  DAMQT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  DAMQT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with DAMQT.  If not, see <http://www.gnu.org/licenses/>.
//
//------------------------------------------------------------------------
//
//  File:   geometryprocessor.h
//  Description: class geometryprocessor contains the OpenGL engline for 3D plotting
//	(openGL 3.3 or higher required : using shaders)
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//
#ifndef GEOMETRYPROCESSOR_H
#define GEOMETRYPROCESSOR_H

#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_4_5_Core>

#include <math.h>
#include "molecule.h"

#define THRESANGLES 0.000001

class geomProcessor : protected QOpenGLFunctions
{
public:
    geomProcessor(QList<molecule*> *m);
    virtual ~geomProcessor();

    void applytranslucence(int);
    void drawCubeGeometry(QOpenGLShaderProgram *program);
    void drawDihedrals(QOpenGLShaderProgram *program, QVector <GLuint> *dihedralindices,
                       QVector <VertexNormalData> *dihedralvertices);
    void drawLabAxes(QOpenGLShaderProgram *program, QVector <GLuint> *axesindices,
                       QVector <VertexNormalData> *axesvertices);
    void drawStructure(QOpenGLShaderProgram *program, QList<molecule*> *m, int indmol);
    void loadbuffers(QList<molecule*> *m);
    void checkbuffers();
//    void setsoliddihedplanes(bool);

private:
    void initCubeGeometry();

//    bool soliddihedplanes;
    GLuint queryID;
    QOpenGLFunctions_4_5_Core *glFuncs;

    QVector<GLuint> *allindices;           // Indices of vertices in structure
    QVector<GLuint> *allindicesoffset;     // Offsets of indices in objects (cylinders and spheres of molecules) in structure
    QVector<GLfloat> *alllineswidth;       // Lines width (only for field lines)
    QVector<GLboolean>  *allsolidsurface;        // true: solid isosurface, false wired frame isosurface
    QVector<GLuint> *allsurfindices;           // Indices of surfaces in structure
    QVector<GLboolean>  *alltranslucence;        // true: solid isosurface, false wired frame isosurface
    QVector<GLuint>  *alltypeelement;         // 0: structure; 1: surface; 2: field lines;
    QVector<GLuint> *moleculesoffsetindex; // Position of indices in allindicesoffset corresponding to each molecule
    QVector<VertexNormalData> *allvertices;  // Vertices of triangles in structure

    QOpenGLBuffer arrayBuf;
    QOpenGLBuffer arrayBuf2;
    QOpenGLBuffer arrayBuf3;
    QOpenGLBuffer indexBuf;
    QOpenGLBuffer indexBuf2;
    QOpenGLBuffer indexBuf3;

    QList<molecule*> *mol;

};

#endif // GEOMETRYPROCESSOR_H
