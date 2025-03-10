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
//	File:   geometryprocessor.cpp
//  Description: implements geometryprocessor class
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: July 2018

#include "geometryprocessor.h"

#include <QVector2D>
#include <QVector3D>
#include <QMessageBox>
#include <QObject>
#include <QDebug>

geomProcessor::geomProcessor(QList<molecule*> *m)
    : arrayBuf(QOpenGLBuffer::VertexBuffer), arrayBuf2(QOpenGLBuffer::VertexBuffer),
      arrayBuf3(QOpenGLBuffer::VertexBuffer),
      indexBuf(QOpenGLBuffer::IndexBuffer), indexBuf2(QOpenGLBuffer::IndexBuffer),
      indexBuf3(QOpenGLBuffer::IndexBuffer)
{
    initializeOpenGLFunctions();
    glFuncs = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_4_5_Core>();
    glFuncs->glGenQueries(1, &queryID);

    // Generate 4 VBOs
    arrayBuf.create();
    indexBuf.create();
    arrayBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    indexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    arrayBuf2.create();
    indexBuf2.create();
    arrayBuf2.setUsagePattern(QOpenGLBuffer::StaticDraw);
    indexBuf2.setUsagePattern(QOpenGLBuffer::StaticDraw);
    arrayBuf3.create();
    indexBuf3.create();
    arrayBuf3.setUsagePattern(QOpenGLBuffer::StaticDraw);
    indexBuf3.setUsagePattern(QOpenGLBuffer::StaticDraw);
    allvertices = new QVector<VertexNormalData>();
    allindices = new QVector<GLuint>();
    allindicesoffset = new QVector<GLuint>();
    alllineswidth = new QVector<GLfloat>();
    alltypeelement = new QVector<GLuint>();
    allsolidsurface = new QVector<GLboolean>();
    allsurfindices = new QVector<GLuint>();
    alltranslucence = new QVector<GLboolean>();
    moleculesoffsetindex = new QVector<GLuint>();
    mol = m;

    loadbuffers(m);

}

geomProcessor::~geomProcessor()
{
    arrayBuf.destroy();
    indexBuf.destroy();
    arrayBuf2.destroy();
    indexBuf2.destroy();
    arrayBuf3.destroy();
    indexBuf3.destroy();
}

void geomProcessor::drawDihedrals(QOpenGLShaderProgram *program, QVector <GLuint> *dihedralindices,
                                  QVector <VertexNormalData> *dihedralvertices)
{
    // Tell OpenGL which VBOs to use
    arrayBuf2.bind();
    indexBuf2.bind();
    //  Transfer dihedral vertexes data to VBO arrayBuf2
    arrayBuf2.allocate(dihedralvertices->constData(), dihedralvertices->length() * sizeof(VertexNormalData));
    //  Transfer dihedral indexes data to VBO indexBuf2
    indexBuf2.allocate(dihedralindices->constData(), dihedralindices->length() * sizeof(GLuint));
    // Offset for position
    quintptr offset = 0;

    // Tell OpenGL programmable pipeline how to locate vertex position data
    int vertexLocation = 0;

    program->enableAttributeArray(vertexLocation);
    program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(VertexNormalData));

    // Offset for normals
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex normals data
    int normalsLocation = 1;
    program->enableAttributeArray(normalsLocation);
    program->setAttributeBuffer(normalsLocation, GL_FLOAT, offset, 3, sizeof(VertexNormalData));

    // Offset for colors
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex colors data
    int colorLocation = 2;
    program->enableAttributeArray(colorLocation);
    program->setAttributeBuffer(colorLocation, GL_FLOAT, offset, 4, sizeof(VertexNormalData));
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glDrawElements(GL_TRIANGLES, dihedralindices->length(), GL_UNSIGNED_INT, (void*)(nullpointer));
    indexBuf2.release();
    arrayBuf2.release();
}

void geomProcessor::drawLabAxes(QOpenGLShaderProgram *program, QVector <GLuint> *axesindices,
                                  QVector <VertexNormalData> *axesvertices)
{
    // Tell OpenGL which VBOs to use
    arrayBuf3.bind();
    indexBuf3.bind();
    //  Transfer dihedral vertexes data to VBO arrayBuf3
    arrayBuf3.allocate(axesvertices->constData(), axesvertices->length() * sizeof(VertexNormalData));
    //  Transfer dihedral indexes data to VBO indexBuf2
    indexBuf3.allocate(axesindices->constData(), axesindices->length() * sizeof(GLuint));
    // Offset for position
    quintptr offset = 0;

    // Tell OpenGL programmable pipeline how to locate vertex position data
    int vertexLocation = 0;

    program->enableAttributeArray(vertexLocation);
    program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(VertexNormalData));

    // Offset for normals
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex normals data
    int normalsLocation = 1;
    program->enableAttributeArray(normalsLocation);
    program->setAttributeBuffer(normalsLocation, GL_FLOAT, offset, 3, sizeof(VertexNormalData));

    // Offset for colors
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex colors data
    int colorLocation = 2;
    program->enableAttributeArray(colorLocation);
    program->setAttributeBuffer(colorLocation, GL_FLOAT, offset, 4, sizeof(VertexNormalData));
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glDrawElements(GL_TRIANGLES, axesindices->length(), GL_UNSIGNED_INT, (void*)(nullpointer));
    indexBuf3.release();
    arrayBuf3.release();
}

void geomProcessor::drawStructure(QOpenGLShaderProgram *program, QList<molecule*> *m, int indmol)
{
    GLuint visibleSamples = 0;
    // Tell OpenGL which VBOs to use
    arrayBuf.bind();
    indexBuf.bind();

    // Offset for position
    quintptr offset = 0;

    // Tell OpenGL programmable pipeline how to locate vertex position data
//    int vertexLocation = program->attributeLocation("a_position");
    int vertexLocation = 0;

    program->enableAttributeArray(vertexLocation);
    program->setAttributeBuffer(vertexLocation, GL_FLOAT, offset, 3, sizeof(VertexNormalData));

    // Offset for normals
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex normals data
//    int normalsLocation = program->attributeLocation("a_normals");
    int normalsLocation = 1;
    program->enableAttributeArray(normalsLocation);
    program->setAttributeBuffer(normalsLocation, GL_FLOAT, offset, 3, sizeof(VertexNormalData));

    // Offset for colors
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex colors data
//    int texcoordLocation = program->attributeLocation("a_colors");
    int colorLocation = 2;
    program->enableAttributeArray(colorLocation);
    program->setAttributeBuffer(colorLocation, GL_FLOAT, offset, 4, sizeof(VertexNormalData));

    if (moleculesoffsetindex->length() < indmol+2) // To prevent a nasty error when resizing some dialog windows
        return;
    int kntmax = 0;
    int kntmin = 0;
    for (int i = moleculesoffsetindex->at(indmol) ; i < moleculesoffsetindex->at(indmol+1); i++){
        if (allsolidsurface->at(i+1)){
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable(GL_BLEND);
        }
        else{
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        if (alltypeelement->at(i+1) == 0){  // Draws structures
            glEnable(GL_CULL_FACE);
            glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
        }
        else if (alltypeelement->at(i+1) == 1){  // Draws surfaces
            if (!alltranslucence->at(i+1)){
                glDisable(GL_CULL_FACE);
                glEnable(GL_DEPTH_TEST);
                glFrontFace(GL_CCW);
                glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                        (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
            }
            else{
                applytranslucence(i);
            }
        }
        else if (alltypeelement->at(i+1) == 2){  // Draws lines
            glLineWidth(alllineswidth->at(i));
            glDrawElements(GL_LINES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
        }
        else if (alltypeelement->at(i+1) == 3){  // Draws arrows
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable(GL_BLEND);
            glDisable(GL_CULL_FACE);
            glDrawElements(GL_TRIANGLE_FAN, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
            glEnable(GL_CULL_FACE);
        }
        if (alltypeelement->at(i+1) == 4){  // Draws structures
            glEnable(GL_CULL_FACE);
            glFuncs->glBeginQuery(GL_SAMPLES_PASSED, queryID);
            glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
            glFuncs->glEndQuery(GL_SAMPLES_PASSED);
            glFuncs->glGetQueryObjectuiv(queryID, GL_QUERY_RESULT, &visibleSamples);
            if (visibleSamples > 0){
//                std::cout << "visible maximum" << "\n";
                m->at(indmol)->surfaces->at(allsurfindices->at(i+1))->setextremhidden(0,kntmax++,false);
                glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                    (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
            }
            else {
//                std::cout << "hidden maximum" << "\n";
                m->at(indmol)->surfaces->at(allsurfindices->at(i+1))->setextremhidden(0,kntmax++,true);
            }
        }
        if (alltypeelement->at(i+1) == 5){  // Draws structures
            glEnable(GL_CULL_FACE);
            glFuncs->glBeginQuery(GL_SAMPLES_PASSED, queryID);
            glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
            glFuncs->glEndQuery(GL_SAMPLES_PASSED);
            glFuncs->glGetQueryObjectuiv(queryID, GL_QUERY_RESULT, &visibleSamples);
            if (visibleSamples > 0){
//                std::cout << "visible minimum" << "\n";
                m->at(indmol)->surfaces->at(allsurfindices->at(i+1))->setextremhidden(1,kntmin++,false);
                glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
                    (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
            }
            else {
//                std::cout << "hidden minimum" << "\n";
                m->at(indmol)->surfaces->at(allsurfindices->at(i+1))->setextremhidden(1,kntmin++,true);
            }
        }
    }
    indexBuf.release();
    arrayBuf.release();
}

void geomProcessor::loadbuffers(QList<molecule*> *m){
    allvertices->clear();
    allindices->clear();
    allindicesoffset->clear();
    alllineswidth->clear();
    alltypeelement->clear();
    allsolidsurface->clear();
    allsurfindices->clear();
    alltranslucence->clear();
    moleculesoffsetindex->clear();
    moleculesoffsetindex->append(0);
    int ioffset = 0;
    for (int i = 0 ; i < m->count() ; i++){     // Loop over molecules
        if ((m->at(i)->isvisible())){
//            Loads structures
            int joffset = allvertices->length();
            ioffset = allindices->count();
            for (int j = 0; j < m->at(i)->allvertices.count() ; j++){
                allvertices->append(m->at(i)->allvertices.at(j));
            }
            for (int j = 0; j < m->at(i)->allindices.count() ; j++){
                allindices->append(m->at(i)->allindices.at(j) + joffset);
            }
            for (int j = 0; j < m->at(i)->allindicesoffset.count() ; j++){
                allindicesoffset->append(m->at(i)->allindicesoffset.at(j) + ioffset);
                alltypeelement->append(0);
                allsurfindices->append(0);
                allsolidsurface->append(true);
                alltranslucence->append(true);
                alllineswidth->append(0.);
            }
//            Loads axes
            if (m->at(i)->isaxesvisible()){
                int joffset = allvertices->length();
                for (int l = 0 ; l < m->at(i)->allaxesvertices.length() ; l++){
                    allvertices->append(m->at(i)->allaxesvertices.at(l));
                }
                for (int l = 0 ; l < m->at(i)->allaxesindices.length() ; l++){
                    allindices->append(m->at(i)->allaxesindices.at(l) + joffset);
                }
                ioffset = allindices->count();
                allindicesoffset->append(ioffset);
                alltypeelement->append(0);
                allsurfindices->append(0);
                allsolidsurface->append(true);
                alltranslucence->append(true);
                alllineswidth->append(0.);
            }
//            Loads surfaces
            for (int j = 0 ; j < m->at(i)->surfaces->count() ; j++){
                if (m->at(i)->surfaces->at(j)->getvisible()){
                    int joffset = allvertices->length();
                    for (int l = 0 ; l < m->at(i)->surfaces->at(j)->allvertices.length() ; l++){
                        allvertices->append(m->at(i)->surfaces->at(j)->allvertices.at(l));
                    }
                    for (int l = 0 ; l < m->at(i)->surfaces->at(j)->allindices.length() ; l++){
                        allindices->append(m->at(i)->surfaces->at(j)->allindices.at(l) + joffset);
                    }
                    ioffset = allindices->count();
                    allindicesoffset->append(ioffset);
                    alltypeelement->append(1);
                    allsurfindices->append(j);
                    if (m->at(i)->surfaces->at(j)->getsolidsurf()){
                        allsolidsurface->append(true);
                    }
                    else{
                        allsolidsurface->append(false);
                    }
                    if (m->at(i)->surfaces->at(j)->gettranslucence() && m->at(i)->surfaces->at(j)->getopacity() < 0.99){
                        alltranslucence->append(true);
                    }
                    else{
                        alltranslucence->append(false);
                    }
                    alllineswidth->append(0.);
                    if (m->at(i)->surfaces->at(j)->getshowgridbounds()){
                        int joffset = allvertices->length();
                        ioffset = allindices->count();
                        for (int l = 0 ; l < m->at(i)->surfaces->at(j)->gridvertices.length() ; l++){
                            allvertices->append(m->at(i)->surfaces->at(j)->gridvertices.at(l));
                        }
                        for (int l = 0 ; l < m->at(i)->surfaces->at(j)->gridindices.length() ; l++){
                            allindices->append(m->at(i)->surfaces->at(j)->gridindices.at(l) + joffset);
                        }
                        for (int l = 1; l < m->at(i)->surfaces->at(j)->gridindicesoffset.count() ; l++){
                            allindicesoffset->append(m->at(i)->surfaces->at(j)->gridindicesoffset.at(l) + ioffset);
                            alltypeelement->append(2);
                            allsurfindices->append(j);
                            allsolidsurface->append(false);
                            alltranslucence->append(false);
                            alllineswidth->append(1.);
                        }
                    }
                    if (m->at(i)->surfaces->at(j)->getshowlocalmax()){
                        int joffset = allvertices->length();
                        ioffset = allindices->count();
                        for (int l = 0 ; l < m->at(i)->surfaces->at(j)->allverticesextrema[0].length() ; l++){
                            allvertices->append(m->at(i)->surfaces->at(j)->allverticesextrema[0].at(l));
                        }
                        for (int l = 0 ; l < m->at(i)->surfaces->at(j)->allindicesextrema[0].length() ; l++){
                            allindices->append(m->at(i)->surfaces->at(j)->allindicesextrema[0].at(l) + joffset);
                        }
                        for (int l = 1; l < m->at(i)->surfaces->at(j)->allindicesoffsetextrema[0].count() ; l++){
                            allindicesoffset->append(m->at(i)->surfaces->at(j)->allindicesoffsetextrema[0].at(l) + ioffset);
                            alltypeelement->append(4);
                            allsurfindices->append(j);
                            allsolidsurface->append(true);
                            alltranslucence->append(true);
                            alllineswidth->append(1.);
                        }
                    }
                    if (m->at(i)->surfaces->at(j)->getshowlocalmin()){
                        int joffset = allvertices->length();
                        ioffset = allindices->count();
                        for (int l = 0 ; l < m->at(i)->surfaces->at(j)->allverticesextrema[1].length() ; l++){
                            allvertices->append(m->at(i)->surfaces->at(j)->allverticesextrema[1].at(l));
                        }
                        for (int l = 0 ; l < m->at(i)->surfaces->at(j)->allindicesextrema[1].length() ; l++){
                            allindices->append(m->at(i)->surfaces->at(j)->allindicesextrema[1].at(l) + joffset);
                        }
                        for (int l = 1; l < m->at(i)->surfaces->at(j)->allindicesoffsetextrema[1].count() ; l++){
                            allindicesoffset->append(m->at(i)->surfaces->at(j)->allindicesoffsetextrema[1].at(l) + ioffset);
                            alltypeelement->append(5);
                            allsurfindices->append(j);
                            allsolidsurface->append(true);
                            alltranslucence->append(true);
                            alllineswidth->append(1.);
                        }
                    }
                }
            }
//            Loads grid isosurfaces
            for (int j = 0 ; j < m->at(i)->grids->count() ; j++){
                for (int k = 0 ; k < m->at(i)->grids->at(j)->surfaces->count() ; k++){
                    if (m->at(i)->grids->at(j)->surfaces->at(k)->isvisible()){
                        int joffset = allvertices->length();
                        for (int l = 0 ; l < m->at(i)->grids->at(j)->surfaces->at(k)->allvertices.length() ; l++){
                            allvertices->append(m->at(i)->grids->at(j)->surfaces->at(k)->allvertices.at(l));
                        }
                        for (int l = 0 ; l < m->at(i)->grids->at(j)->surfaces->at(k)->allindices.length() ; l++){
                            allindices->append(m->at(i)->grids->at(j)->surfaces->at(k)->allindices.at(l) + joffset);
                        }
                        ioffset = allindices->count();
                        allindicesoffset->append(ioffset);
                        alltypeelement->append(1);
                        allsurfindices->append(k);
                        if (m->at(i)->grids->at(j)->surfaces->at(k)->getsolidsurf()){
                            allsolidsurface->append(true);
                        }
                        else{
                            allsolidsurface->append(false);
                        }
                        if (m->at(i)->grids->at(j)->surfaces->at(k)->gettranslucence() && m->at(i)->grids->at(j)->surfaces->at(k)->getopacity() < 0.99){
                            alltranslucence->append(true);
                        }
                        else{
                            alltranslucence->append(false);
                        }
                        alllineswidth->append(0.);
                        if (m->at(i)->grids->at(j)->surfaces->at(k)->getshowgridbound()){
                            int joffset = allvertices->length();
                            ioffset = allindices->count();
                            for (int l = 0 ; l < m->at(i)->grids->at(j)->surfaces->at(k)->gridvertices.length() ; l++){
                                allvertices->append(m->at(i)->grids->at(j)->surfaces->at(k)->gridvertices.at(l));
                            }
                            for (int l = 0 ; l < m->at(i)->grids->at(j)->surfaces->at(k)->gridindices.length() ; l++){
                                allindices->append(m->at(i)->grids->at(j)->surfaces->at(k)->gridindices.at(l) + joffset);
                            }
                            for (int l = 1; l < m->at(i)->grids->at(j)->surfaces->at(k)->gridindicesoffset.count() ; l++){
                                allindicesoffset->append(m->at(i)->grids->at(j)->surfaces->at(k)->gridindicesoffset.at(l) + ioffset);
                                alltypeelement->append(2);
                                allsurfindices->append(k);
                                allsolidsurface->append(false);
                                alltranslucence->append(false);
                                alllineswidth->append(1.);
                            }
                        }
                    }
                }
            }
//            Loads Hellmann-Feynman forces
            if (m->at(i)->hfforces && m->at(i)->hfforces->isvisible()){
                int joffset = allvertices->length();
                for (int l = 0 ; l < m->at(i)->hfforces->allvertices.length() ; l++){
                    allvertices->append(m->at(i)->hfforces->allvertices.at(l));
                }
                for (int l = 0 ; l < m->at(i)->hfforces->allindices.length() ; l++){
                    allindices->append(m->at(i)->hfforces->allindices.at(l) + joffset);
                }
                ioffset = allindices->count();
                allindicesoffset->append(ioffset);
                alltypeelement->append(0);
                allsurfindices->append(0);
                allsolidsurface->append(true);
                alltranslucence->append(true);
                alllineswidth->append(0.);
            }
//            Loads field lines
            if (m->at(i)->flines && m->at(i)->flines->isvisible()){
                int joffset = allvertices->length();
                ioffset = allindices->count();
                for (int l = 0 ; l < m->at(i)->flines->allvertices.length() ; l++){
                    allvertices->append(m->at(i)->flines->allvertices.at(l));
                }
                for (int l = 0 ; l < m->at(i)->flines->allindices.length() ; l++){
                    allindices->append(m->at(i)->flines->allindices.at(l) + joffset);
                }
                for (int j = 1; j < m->at(i)->flines->allindicesoffset.count() ; j++){
                    allindicesoffset->append(m->at(i)->flines->allindicesoffset.at(j) + ioffset);
                    alltypeelement->append(2);
                    allsurfindices->append(0);
                    allsolidsurface->append(true);
                    alltranslucence->append(true);
                    alllineswidth->append(m->at(i)->flines->getlineswidth());
                }
                if (m->at(i)->flines->getshowarrows()){
                    int joffset = allvertices->length();
                    ioffset = allindices->count();
                    for (int l = 0 ; l < m->at(i)->flines->allarrowsvertices.length() ; l++){
                        allvertices->append(m->at(i)->flines->allarrowsvertices.at(l));
                    }
                    for (int l = 0 ; l < m->at(i)->flines->allarrowsindices.length() ; l++){
                        allindices->append(m->at(i)->flines->allarrowsindices.at(l) + joffset);
                    }
                    for (int j = 1; j < m->at(i)->flines->allarrowsindicesoffset.count() ; j++){
                        allindicesoffset->append(m->at(i)->flines->allarrowsindicesoffset.at(j) + ioffset);
                        alltypeelement->append(3);
                        allsurfindices->append(0);
                        allsolidsurface->append(true);
                        alltranslucence->append(true);
                        alllineswidth->append(m->at(i)->flines->getlineswidth());
                    }
                }
            }
//            Loads critical points
            if (m->at(i)->cps){
                for (int type = 0 ; type < MAX_CPS ; type++){
                    if (!m->at(i)->cps->isvisiblecps(type) || m->at(i)->cps->allvertices[type].length() < 1)
                        continue;
                    int joffset = allvertices->length();
                    ioffset = allindices->count();
                    for (int l = 0 ; l < m->at(i)->cps->allvertices[type].length() ; l++){
                        allvertices->append(m->at(i)->cps->allvertices[type].at(l));
                    }
                    for (int l = 0 ; l < m->at(i)->cps->allindices[type].length() ; l++){
                        allindices->append(m->at(i)->cps->allindices[type].at(l) + joffset);
                    }
                    for (int j = 1; j < m->at(i)->cps->allindicesoffset[type].count() ; j++){
                        allindicesoffset->append(m->at(i)->cps->allindicesoffset[type].at(j) + ioffset);
                        alltypeelement->append(0);
                        allsurfindices->append(0);
                        allsolidsurface->append(true);
                        alltranslucence->append(true);
                        alllineswidth->append(0.);
                    }
                }
                if (m->at(i)->cps->isvisiblecpseigen()){
                    for (int type = 0 ; type < MAX_CPS ; type++){
                        if (!m->at(i)->cps->isvisiblecps(type) || m->at(i)->cps->alleigvertices[type].length() < 1)
                            continue;
                        int joffset = allvertices->length();
                        ioffset = allindices->count();
                        for (int l = 0 ; l < m->at(i)->cps->alleigvertices[type].length() ; l++){
                            allvertices->append(m->at(i)->cps->alleigvertices[type].at(l));
                        }
                        for (int l = 0 ; l < m->at(i)->cps->alleigindices[type].length() ; l++){
                            allindices->append(m->at(i)->cps->alleigindices[type].at(l) + joffset);
                        }
                        for (int j = 1; j < m->at(i)->cps->alleigindicesoffset[type].count() ; j++){
                            allindicesoffset->append(m->at(i)->cps->alleigindicesoffset[type].at(j) + ioffset);
                            alltypeelement->append(0);
                            allsurfindices->append(0);
                            allsolidsurface->append(true);
                            alltranslucence->append(true);
                            alllineswidth->append(0.);
                        }
                    }
                }
            }
        }
        if (allindicesoffset->count() > 0){
            moleculesoffsetindex->append(allindicesoffset->count()-1);
        }
        else
            moleculesoffsetindex->append(0);
    }

//  Transfer vertex data to VBO arrayBuf

    arrayBuf.bind();
    arrayBuf.allocate(allvertices->constData(), allvertices->length() * sizeof(VertexNormalData));

//  Transfer index data to VBO indexBuf

    indexBuf.bind();
    indexBuf.allocate(allindices->constData(), allindices->length() * sizeof(GLuint));
    indexBuf.release();
    arrayBuf.release();
}

void geomProcessor::checkbuffers(){
    if (arrayBuf.size() < 0 || indexBuf.size() < 0){
        QMessageBox msgBox;
        msgBox.setText(QObject::tr("loadbuffers"));
        msgBox.setInformativeText(QObject::tr("Failed allocating buffers")
                + QObject::tr("\nIf display is not correctly updated after closing this warning, the following workaround may help: ")
                + QObject::tr("\nin surface edit window, push hide/show button")
                + QObject::tr("\nPush twice if required.")
                );
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
    }
}

//  Applies translucence correction by drawing the surface five times with different openGL settings
//
void geomProcessor::applytranslucence(int i){
    double f = 0.75; // Or some other factor
//    glDisable(GL_CULL_FACE);
//    glEnable(GL_DEPTH_TEST);
//    if (alltypeelement->at(i+1) == 1){  // This is to keep the front faces brighter thatn the back faces
//        glFrontFace(GL_CCW);     // Positive contour faces are counterclockwise
//    }
//    else{
//        glFrontFace(GL_CW);    // Negative contour faces are clockwise
//    }
//    glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
//            (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
//    //	First view of surfaces
//    glEnable(GL_BLEND);
//    glDisable(GL_CULL_FACE);
//    glDepthFunc(GL_LESS);
////    glBlendFunc(GL_ZERO, GL_ONE);
//    glBlendFunc(GL_ONE, GL_ZERO);
//    glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
//            (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
    //	Second view of surfaces: render with alpha = f*alpha
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glDepthFunc(GL_ALWAYS);
    glBlendFunc(f*GL_SRC_ALPHA, GL_ONE-f*GL_SRC_ALPHA);
    glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
            (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
    //	Third view of surfaces: render with alpha = (alpha-f*alpha)/(1.0-f*alpha)
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glDepthFunc(GL_LEQUAL);
    glBlendFunc((1-f)*GL_SRC_ALPHA/(GL_ONE-f*GL_SRC_ALPHA),
            GL_ONE-(1-f)*GL_SRC_ALPHA/(GL_ONE-f*GL_SRC_ALPHA));
    glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
            (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
    //	Fourth view of surfaces: render with render with alpha = f*alpha
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glDepthFunc(GL_ALWAYS);
    glBlendFunc(f*GL_SRC_ALPHA, GL_ONE-f*GL_SRC_ALPHA);
    glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
            (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
    //	Fifth view of surfaces: render with render with alpha = (alpha-f*alpha)/(1.0-f*alpha)
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LEQUAL);
    glBlendFunc((1-f)*GL_SRC_ALPHA/(GL_ONE-f*GL_SRC_ALPHA),
            GL_ONE-(1-f)*GL_SRC_ALPHA/(GL_ONE-f*GL_SRC_ALPHA));
    glDrawElements(GL_TRIANGLES, allindicesoffset->at(i+1)-allindicesoffset->at(i), GL_UNSIGNED_INT,
            (GLvoid *) (allindicesoffset->at(i) * sizeof(GLuint)));
}

