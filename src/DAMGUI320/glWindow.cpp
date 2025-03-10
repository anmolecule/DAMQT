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
//    Class for defining an OpenGL window for 3D display
//  
//    File:   glWindow.cpp
//
//    Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//    Last version: January 2019
//

#include <math.h>

#include <QApplication>
#include <QMessageBox>
#include <QObject>
#include <QMouseEvent>
#include <QtDebug>
#include <QDir>
#include <QWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>
#include <QShortcut>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>

#include "glWindow.h"

glWindow::glWindow(QList<molecule*> *mol, QWidget *parent) :
    QOpenGLWindow(),
    geometries(0),
    texture(0),
    angularSpeed(0)
{
    // Set OpenGL Version information
    // Note: This format must be set before show() is called.
    QSurfaceFormat format;
    format.setRenderableType(QSurfaceFormat::OpenGL);
    #if QT_VERSION >= 0x050900
        format.setProfile(QSurfaceFormat::CoreProfile); // To be uncommented for Qt 5.9 or higher
    #endif
    format.setVersion(3,3);
    format.setDepthBufferSize(24);
    format.setAlphaBufferSize(24);
    format.setStencilBufferSize(8);

    this->setFormat(format);

    allaxesindices = new QVector <GLuint>();
    allaxesindicesoffset = new QVector <GLuint>();
    allaxesvertices = new QVector <VertexNormalData>();
    Ambient_Color = QVector3D(0.5,0.5,0.5);
    anglecenters = new QVector<centerData>();
    angles = false;
    anglescolor = QColor(Qt::yellow);
    anglesfont = QFont("Noto Sans", 12, QFont::Bold);
    anglesprecision = 4;
    anglestype = 3; // Dot Line
    angleswidth = 2;
    angstrom = false;
    arcradius = 40;
    dihedralcenters = new QVector<centerData>();
    dihedralindices = new QVector <GLuint>();
    dihedralvertices = new QVector <VertexNormalData>();
    dihedrals = false;
    dihedralscolor = QColor(Qt::yellow);
    dihedralplanescolor = QColor(182, 182, 182, 128);
    dihedralsfont = QFont("Noto Sans", 12, QFont::Bold);
    distancecenters = new QVector<centerData>();
    distances = false;
    distancescolor = QColor(Qt::yellow);
    distancesfont = QFont("Noto Sans", 12, QFont::Bold);
    disthres = 1.2 * ANGSTROM_TO_BOHR;
    disttranspbkg = false;
    dihedralsprecision = 4;
    displayEPIC = false;
    distprecision = 2;
    distvshift = 0;
    drawarcs = true;
    drawlines = true;
    elem = new Elements();
    epicenergy = 0.;
    EPICposition = QPoint(0,0);
    EPICpressed = false;
    energycolor  = QColor(255, 172, 0, 255);
    energyfont = QFont("Helvetica", 18, QFont::Bold);
    energyprecision = 4;
    EPICstring = "";
    fbo = nullpointer;
    fov = PERSPECTIVE_ANGLE;
    frameknt = 0;
    frames.clear();
    hartree = false;
    imagequality = 20;
    linestype = 3; // Dot Line
    lineswidth = 2;
    linearattenuation = false;
    lightPosition = QVector3D(0,0,0);
    Light_Color = QVector3D(1,1,1);
    Light_Power = 100.f;
    molecules = mol;
    mvp_list = new QVector<QMatrix4x4>();
    mvp_list_backup = new QVector<QMatrix4x4>();
    positionaxeslabels = new QVector <QVector3D>();
    record = false;
    recordfilename = QString("Unnamed");
    removeframes = false;
    setimagefilename(QString(""));
    setmeasures(false);
    settransparentbg(false);
    showangles = false;
    showaxes = false;
    showdihedrals = false;
    showdistances = false;
    Specular_Color = QVector3D(0.5,0.5,0.5);
    Specular_Index = 5.;
    v_list = new QVector<QMatrix4x4>();
    world_rotationAxis = QVector3D(1,0,0);   
    zFar = ZFAR;
    zNear = ZNEAR;
    makeAxesCylinder(15,15);    // computes vertices of a cylinder and its indices
    makeAxesCone(15,15,2.);    // computes vertices of a cone and its indices
    resetlabaxes();
}

glWindow::~glWindow()
{
    // Make sure the context is current when deleting the texture and the buffers.
    makeCurrent();
    delete texture;
    texture = nullpointer;
    delete geometries;
    geometries = nullpointer;
    delete elem;
    elem = nullpointer;
    doneCurrent();
}

void glWindow::closeEvent(QCloseEvent *event){
    event->ignore();
}

void glWindow::keyPressEvent(QKeyEvent *e){
    if (e->key() == Qt::Key_Escape)
        return;
    else if (e->key() == Qt::Key_W){
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                        + QVector3D(0,0,molecules->at(i)->getstepwheel()));
            }
        }
        if (!activemols){
            world_translation += QVector3D(0,0,stepwheel);
            emit_update_worldtranslation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_S){
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                        - QVector3D(0,0,molecules->at(i)->getstepwheel()));
            }
        }
        if (!activemols){
            world_translation -= QVector3D(0,0,stepwheel);
            emit_update_worldtranslation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_A){
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                        - QVector3D(molecules->at(i)->getstepwheel(),0,0));
            }
        }
        if (!activemols){
            world_translation -= QVector3D(stepwheel,0,0);
            emit_update_worldtranslation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_D){
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                        + QVector3D(molecules->at(i)->getstepwheel(),0,0));
            }
        }
        if (!activemols){
            world_translation += QVector3D(stepwheel,0,0);
            emit_update_worldtranslation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_Q){
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                        + QVector3D(0,molecules->at(i)->getstepwheel(),0));
            }
        }
        if (!activemols){
            world_translation += QVector3D(0,stepwheel,0);
            emit_update_worldtranslation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_Z){
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                        - QVector3D(0,molecules->at(i)->getstepwheel(),0));
            }
        }
        if (!activemols){
            world_translation -= QVector3D(0,stepwheel,0);
            emit_update_worldtranslation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_R){
        QVector3D n = QVector3D(0.0, 1.0, 0.0).normalized();
        qreal acc = 1.;
        if (e->text() == QString('R')) acc = -acc;  // If capital letter, reverses rotation
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->setrotationAxis((angularSpeed * molecules->at(i)->getrotationAxis() + n * abs(acc)).normalized());
                molecules->at(i)->setrotation(QQuaternion::fromAxisAndAngle(molecules->at(i)->getrotationAxis(), acc)
                                             * molecules->at(i)->getrotation());
            }
        }
        if (!activemols){
            setworld_rotationAxis((angularSpeed * getrotationAxis() + n * abs(acc)).normalized());
            setworld_rotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), acc) * getrotation());
            emit_update_worldtrotation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_E){
        bool activemols = false;
        QVector3D n = QVector3D(1.0, 0.0, 0.0).normalized();
        qreal acc = 1.;
        if (e->text() == QString('E')) acc = -acc;  // If capital letter, reverses rotation
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->setrotationAxis((angularSpeed * molecules->at(i)->getrotationAxis() + n * abs(acc)).normalized());
                molecules->at(i)->setrotation(QQuaternion::fromAxisAndAngle(molecules->at(i)->getrotationAxis(), acc)
                                             * molecules->at(i)->getrotation());
            }
        }
        if (!activemols){
            setworld_rotationAxis((angularSpeed * getrotationAxis() + n * abs(acc)).normalized());
            setworld_rotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), acc) * getrotation());
            emit_update_worldtrotation();
        }
        update();
        return;
    }
    else if (e->key() == Qt::Key_F){
        QVector3D n = QVector3D(0.0, 0.0, 1.0).normalized();
        qreal acc = 1.;
        if (e->text() == QString('F')) acc = -acc;  // If capital letter, reverses rotation
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i ++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->setrotationAxis((angularSpeed * molecules->at(i)->getrotationAxis() + n * abs(acc)).normalized());
                molecules->at(i)->setrotation(QQuaternion::fromAxisAndAngle(molecules->at(i)->getrotationAxis(), acc)
                                             * molecules->at(i)->getrotation());
            }
        }
        if (!activemols){
            setworld_rotationAxis((angularSpeed * getrotationAxis() + n * abs(acc)).normalized());
            setworld_rotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), acc) * getrotation());
            emit_update_worldtrotation();
        }
        update();
        return;
    }
}

// Mouse functions and related functions
// =====================================


/* Defines mouse double click event */
void glWindow::mouseDoubleClickEvent(QMouseEvent *event)
{
    const Qt::KeyboardModifiers modifiers = event->modifiers();
    int idistance;
    int iangle;
    int x = event->x();
    int y = event->y();
    this->makeCurrent();
    if (modifiers.testFlag(Qt::ShiftModifier)){     //  Doubleclick + shift modifier
        if (measures){
            centerData center;
            center = searchatom(x,y);
            if (center.number < 0) center = searchcps(x,y);
            if (center.number < 0) center = searchextrema(x,y);
            if (center.molecule >= 0){
                if (angles){
                    if (anglecenters->length()%3 != 0){
                        for (int i = 0 ; i < anglecenters->length()%3 ; i++){
                            if (center.molecule == anglecenters->at(anglecenters->length()-i-1).molecule
                                    &&center.symbol == anglecenters->at(anglecenters->length()-i-1).symbol){
                                QMessageBox msgBox;
                                msgBox.setText(tr("Angles selection"));
                                msgBox.setInformativeText(tr("Same center chosen twice.")
                                       +"\nCurrent selection deleted.");
                                msgBox.setIcon(QMessageBox::Warning);
                                msgBox.exec();
                                anglecenters->removeLast();
                                if (anglecenters->length()%3 != 0)
                                    anglecenters->removeLast();
                                center.number = -1;       // To notice that this center is not allowed
                                center.molecule = -1;
                                center.type = -1;
                                center.znuc = -1;
                                center.symbol = "";
                                anglecenters->append(center);
                                emit select_angle(anglecenters,v_list);
                                return;
                            }
                        }
                    }
                    anglecenters->append(center);
                    emit select_angle(anglecenters,v_list);
                }
                else if (dihedrals){
                    for (int i = 0 ; i < dihedralcenters->length()%4 ; i++){
                        if (center.molecule == dihedralcenters->at(dihedralcenters->length()-i-1).molecule
                                && center.symbol == dihedralcenters->at(dihedralcenters->length()-i-1).symbol){
                            QMessageBox msgBox;
                            msgBox.setText(tr("Dihedral angles selection"));
                            msgBox.setInformativeText(tr("Same center chosen twice.")
                                   +"\nCurrent selection deleted.");
                            msgBox.setIcon(QMessageBox::Warning);
                            msgBox.exec();
                            center.number = -1;       // To notice that this center is not allowed
                            center.molecule = -1;
                            center.type = -1;
                            center.znuc = -1;
                            center.symbol = "";
                            dihedralcenters->append(center);
                            emit select_dihedral(dihedralcenters,v_list);
                            return;
                        }
                    }
                    dihedralcenters->append(center);
                    emit select_dihedral(dihedralcenters,v_list);
                }
                else if (distances){
                    if (distancecenters->length()%2 != 0){
                        if (center.molecule == distancecenters->last().molecule && center.symbol == distancecenters->last().symbol){
                            QMessageBox msgBox;
                            msgBox.setText(tr("Distances selection"));
                            msgBox.setInformativeText(tr("Same center chosen twice.")
                                   +"\nCurrent selection deleted.");
                            msgBox.setIcon(QMessageBox::Warning);
                            msgBox.exec();
//                            distancecenters->removeLast();
                            center.number = -1;       // To notice that this center is not allowed
                            center.molecule = -1;
                            center.type = -1;
                            center.znuc = -1;
                            center.symbol = "";
                            distancecenters->append(center);
                            emit select_distance(distancecenters,v_list);
                            return;
                        }
                    }
                    distancecenters->append(center);
                    emit select_distance(distancecenters,v_list);
                }
                update();
                return;
            }   
        }   
        if (selectPopUpWindow(x,y)){
            update();
            return;
        }
        return;
    }       //   End of doubleclick + Shift modifier
    if (modifiers.testFlag(Qt::ControlModifier)){   // Toggles molecule activation
        searchmolecule(x,y);
        update();
        return;
    }
    else if (selectatom(x,y)){
        update();
        return;
    }
    else if(selectcps(x,y)){
        update();
        return;
    }
    else if(selectmespextrema(x,y)){
        update();
        return;
    }
    else if (showdistances && (idistance = searchdistance(x,y)) >= 0){
        update();
        emit resetlastselectdist();
        emit update_distances(distancecenters,v_list);
        return;
    }
    else if (showangles && (iangle = searchangle(x,y)) >= 0){
        update();
        emit resetlastselectangles();
        emit update_angles(anglecenters,v_list);
        return;
    }
    else if (showdihedrals && (iangle = searchdihedral(x,y)) >= 0){
        update();
        emit resetlastselectdihedrals();
        emit update_dihedrals(dihedralcenters,v_list);
        return;
    }
}

void glWindow::mousePressEvent(QMouseEvent *e)
{
    // Save mouse press position
    mousePressPosition = QVector2D(e->localPos());
    if (displayEPIC){
        EPICpressed = searchEPICenergy(mousePressPosition);
    }
    mousePreviousPosRot = mousePressPosition;
    mousePreviousPosTrs = mousePressPosition;
}

void glWindow::wheelEvent(QWheelEvent *e)
{
    bool activemols = false;
    for (int i = 0 ; i < molecules->count() ; i ++){
        if (molecules->at(i)->isactive()){
            activemols = true;
            molecules->at(i)->settranslation(molecules->at(i)->gettranslation()
                    + QVector3D(0,0,e->angleDelta().y()*molecules->at(i)->getstepwheel()/120));
        }
    }
    if (!activemols){
        world_translation -= QVector3D(0,0,e->angleDelta().y()*stepwheel/120);
        emit_update_worldtranslation();
    }
    update();
    return;
}

void glWindow::mouseMoveEvent(QMouseEvent *e)
{
    if (EPICpressed){
        EPICposition += QPoint(e->localPos().x(),e->localPos().y()) - QPoint(mousePressPosition.x(),mousePressPosition.y());
        mousePressPosition = QVector2D(e->localPos());
        update();
        return;
    }
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    // Mouse release position - mouse press position
    QVector2D diff = QVector2D(e->localPos()) - mousePressPosition;
    QVector2D difftrs = QVector2D(e->localPos()) - mousePreviousPosTrs;
    QVector2D diffrot = (QVector2D(e->localPos()) - mousePreviousPosRot) ;
    diffrot = diffrot * diff.length() / std::max(0.00001f,diffrot.length());

    QMatrix4x4 v_matrix;
    v_matrix.rotate(getrotation());

    // Accelerate angular speed relative to the length of the mouse sweep
    qreal acc = diff.length() / 100.0;
    if (e->buttons() & Qt::LeftButton) {
        if(modifiers.testFlag(Qt::ShiftModifier)){
            // Translation in plane XY
            bool activemols = false;
            QVector3D translation = QVector3D(difftrs.x(), -difftrs.y(),0)/100;
            for (int i = 0 ; i < molecules->count() ; i ++){
                if (molecules->at(i)->isactive()){
                    activemols = true;
                    molecules->at(i)->settranslation(molecules->at(i)->gettranslation() + translation*v_matrix);
                }
            }
            if (!activemols){
                world_translation += translation;
                emit_update_worldtranslation();
            }
            mousePreviousPosTrs = QVector2D(e->localPos());
            update();
        }
        else if(diffrot.length() > 1.e-7){
            // Rotation axis is perpendicular to the mouse position difference in plane XY
            QVector3D n = QVector3D(diffrot.y(), diffrot.x(), 0.0).normalized();
            // Calculate new rotation axis as weighted sum
            bool activemols = false;
            for (int i = 0 ; i < molecules->count() ; i ++){
                if (molecules->at(i)->isactive()){
                    activemols = true;
                    molecules->at(i)->setrotationAxis((angularSpeed * molecules->at(i)->getrotationAxis()
                                        + n * acc * v_matrix).normalized());
//                    molecules->at(i)->setrotationAxis((angularSpeed * molecules->at(i)->getrotationAxis()
//                                        + v_matrix * n * acc).normalized());
                }
            }
            if (!activemols){
                setworld_rotationAxis((angularSpeed * getrotationAxis() + n * abs(acc)).normalized());
                for (int i = 0 ; i < molecules->count() ; i ++){
                    molecules->at(i)->setworld_rotation(world_rotation);
                }
            }
            mousePreviousPosRot = QVector2D(e->localPos());
            angularSpeed += acc;
        }
    }
    else if (e->buttons() & Qt::RightButton) {
        if(modifiers.testFlag(Qt::ShiftModifier)){
            // Zoom
            bool activemols = false;
            QVector3D translation = - QVector3D(0,0,difftrs.y())/100;
            for (int i = 0 ; i < molecules->count() ; i ++){
                if (molecules->at(i)->isactive()){
                    activemols = true;
                    molecules->at(i)->settranslation(molecules->at(i)->gettranslation() + translation*v_matrix);
                }
            }
            if (!activemols){
                world_translation += translation;
                emit_update_worldtranslation();
            }
            mousePreviousPosTrs = QVector2D(e->localPos());
            update();
        }
        else{
            // Rotation axis is Z axis
            float df = -diffrot.x();
            if (abs(df) < abs(diffrot.y())) df = -diffrot.y();
            QVector3D n = QVector3D(0.0, 0.0, df).normalized();
            if (n.length() > 1.e-7){
                // Calculate new rotation axis as weighted sum
                bool activemols = false;
                for (int i = 0 ; i < molecules->count() ; i ++){
                    if (molecules->at(i)->isactive()){
                        activemols = true;
                        molecules->at(i)->setrotationAxis((angularSpeed * molecules->at(i)->getrotationAxis()
                                        + n * acc * v_matrix).normalized());
                    }
                }
                if (!activemols){
                    setworld_rotationAxis((angularSpeed * getrotationAxis() + n * abs(acc)).normalized());
                    for (int i = 0 ; i < molecules->count() ; i ++){
                        molecules->at(i)->setworld_rotation(world_rotation);
                    }
                }
                mousePreviousPosRot = QVector2D(e->localPos());
                angularSpeed += acc;
            }
        }
    }
}

void glWindow::mouseReleaseEvent(QMouseEvent *e)
{
    // Mouse release position - mouse press position
    EPICpressed = false;
    return;

}

void glWindow::timerEvent(QTimerEvent *)
{
    // Decrease angular speed (friction)
//    angularSpeed *= 0.99;
    angularSpeed *= 0.8;

    // Stop rotation when speed goes below threshold
    if (angularSpeed < 0.01) {
        angularSpeed = 0.0;
    } else {
        // Update rotation
        bool activemols = false;
        for (int i = 0 ; i < molecules->count() ; i++){
            if (molecules->at(i)->isactive()){
                activemols = true;
                molecules->at(i)->setrotation(QQuaternion::fromAxisAndAngle(molecules->at(i)->getrotationAxis(), angularSpeed)
                                             * molecules->at(i)->getrotation());
            }
        }
        if (!activemols){
            setworld_rotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), angularSpeed) * getrotation());
            emit_update_worldtrotation();
        }
        // Request an update
        update();
    }
}

// Painting functions and related functions
// ========================================

void glWindow::initializeGL()
{
    initializeOpenGLFunctions();

    glClearColor(0, 0, 0, 1);

    initShaders();
//    initTextures();

    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);

    // Enable back face culling
    glEnable(GL_CULL_FACE);

    // Creates the geometry processing engine
    geometries = new geomProcessor(molecules);
    geometries->loadbuffers(molecules);

    // Use QBasicTimer because its faster than QTimer
    timer.start(12, this);

    printContextInformation();

}

void glWindow::initShaders()
{
    // Compile vertex shader
    if (!program.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/VShader.vert")){
        QMessageBox msgBox;
        msgBox.setText(tr("initShaders"));
        msgBox.setInformativeText(tr("Error opening vertex shader"));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        close();
    }

    // Compile fragment shader
    if (!program.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/FShader.frag")){
        QMessageBox msgBox;
        msgBox.setText(tr("initShaders"));
        msgBox.setInformativeText(tr("Error opening fragment shader"));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        close();
    }

    // Link shader pipeline
    if (!program.link()){
        QMessageBox msgBox;
        msgBox.setText(tr("initShaders"));
        msgBox.setInformativeText(tr("Error linking shader pipeline"));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        close();
    }

    // Bind shader pipeline for use
    if (!program.bind()){
        QMessageBox msgBox;
        msgBox.setText(tr("initShaders"));
        msgBox.setInformativeText(tr("Error binding shader pipeline"));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        close();
    }
}

void glWindow::initTextures()
{
    // Load cube.png image
//    texture = new QOpenGLTexture(QImage(QString(app_path).append("/cube.png")).mirrored());
//    texture = new QOpenGLTexture(QImage(":/cube.png").mirrored());

    // Set nearest filtering mode for texture minification
//    texture->setMinificationFilter(QOpenGLTexture::Nearest);

    // Set bilinear filtering mode for texture magnification
//    texture->setMagnificationFilter(QOpenGLTexture::Linear);

    // Wrap texture coordinates by repeating
    // f.ex. texture coordinate (1.1, 1.2) is same as (0.1, 0.2)
//    texture->setWrapMode(QOpenGLTexture::Repeat);
}

void glWindow::resizeGL(int w, int h)
{
    // Calculate aspect ratio
    qreal aspect = qreal(w) / qreal(h ? h : 1);

    // Reset projection
    projection.setToIdentity();

    // Set perspective projection
    projection.perspective(fov, aspect, zNear, zFar);
}

void glWindow::paintGL()
{
    QPainter painter(this);
    this->makeCurrent();

    mvp_list->clear();
    v_list->clear();

    if (!molecules)
        return;

//    emit setrotationbuttons();

    painter.beginNativePainting();

        // Clear color and depth buffer
        glClearColor(Bkg_Color.x(), Bkg_Color.y(), Bkg_Color.z(), 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glShadeModel( GL_SMOOTH );

        // Calculate aspect ratio
        qreal aspect = qreal(width()) / qreal(height() ? height() : 1);
        // Reset projection
        projection.setToIdentity();

        // Set perspective projection
        projection.perspective(fov, aspect, zNear, zFar);

        //  Model matrix (set to identity: all models are suppossed to be generated this way)
        QMatrix4x4 m_matrix;
        m_matrix.setToIdentity();

        QVector3D LightPosition_worldspace = lightPosition;

        program.bind();
        for (int i = 0 ; i < molecules->count() ; i++){
            //  View matrix
            QMatrix4x4 v_matrix;
            v_matrix.translate(world_translation);
            v_matrix.rotate(getrotation());
            v_matrix.translate(molecules->at(i)->gettranslation());
            v_matrix.rotate(molecules->at(i)->getrotation());

            v_list->append(v_matrix);

            molecules->at(i)->setrotationButtons();
            molecules->at(i)->settranslationButtons();

            // Model-view transformation
            QMatrix4x4 mv_matrix;
            mv_matrix = v_matrix * m_matrix;

            // Modelview-projection matrix
            QMatrix4x4 mvp_matrix;
            mvp_matrix = projection * mv_matrix;
            mvp_list->append(mvp_matrix);

            // Set model matrix
            program.setUniformValue("m_matrix", m_matrix);

            // Set view matrix
            program.setUniformValue("v_matrix", v_matrix);

            // Set model-view transformation
            program.setUniformValue("mv_matrix", mv_matrix);

            // Set modelview-projection matrix
            program.setUniformValue("mvp_matrix", mvp_matrix);

            // Set lights properties
            program.setUniformValue("LightPosition_worldspace", LightPosition_worldspace);
            program.setUniformValue("Ambient_Color", Ambient_Color);
            program.setUniformValue("Light_Color", Light_Color);
            program.setUniformValue("Light_Power", Light_Power);
            program.setUniformValue("Specular_Color", Specular_Color);
            program.setUniformValue("Specular_Index", Specular_Index);
            program.setUniformValue("Linear_Attenuation", linearattenuation);

            geometries->drawStructure(&program, molecules, i);
        }

        if (showdihedrals && measures){
            makedihedralplanes();
            if (dihedralindices->length() > 3){
                QMatrix4x4 v_matrix;
                v_matrix.setToIdentity();
                v_matrix.rotate(getrotation());
                QMatrix4x4 mv_matrix;
                mv_matrix.setToIdentity();
                QMatrix4x4 mvp_matrix;
                mvp_matrix = projection;
                program.setUniformValue("m_matrix", m_matrix);
                program.setUniformValue("v_matrix", v_matrix);
                program.setUniformValue("mv_matrix", mv_matrix);
                program.setUniformValue("mvp_matrix", mvp_matrix);

                geometries->drawDihedrals(&program, dihedralindices, dihedralvertices);
            }
        }
        if (showaxes){
            QMatrix4x4 v_matrix;
            v_matrix.translate(world_translation);
            v_matrix.rotate(getrotation());
            // Model-view transformation
            QMatrix4x4 mv_matrix;
            mv_matrix = v_matrix * m_matrix;
            // Modelview-projection matrix
            QMatrix4x4 mvp_matrix;
            mvp_matrix = projection * mv_matrix;
            mvp_list->append(mvp_matrix);
            program.setUniformValue("m_matrix", m_matrix);
            program.setUniformValue("v_matrix", v_matrix);
            program.setUniformValue("mv_matrix", mv_matrix);
            program.setUniformValue("mvp_matrix", mvp_matrix);

            geometries->drawLabAxes(&program, allaxesindices, allaxesvertices);
        }
        program.release();

        GLint vport[4];
        glGetIntegerv (GL_VIEWPORT, vport);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // This is necessary to keep the atom symbols on top even when wired surface is chosen
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

    painter.endNativePainting();

    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    drawlabels(&painter, viewport);
    if (measures){
        drawmeasures(&painter, viewport);
    }
    if (displayEPIC){
        drawEPIC(&painter, viewport);
    }

    painter.end();
    if (record && frameknt++ < numframes){
        this->makeCurrent();
#if QT_VERSION < QT_VERSION_CHECK(5,0,0)
        frames.append(QPixmap::grabWindow(this->winId()));
#else
        QScreen *screen = this->screen();
        frames.append(screen->grabWindow(this->winId(),0,0,-1,-1));
#endif
        if (frameknt >= numframes){
            record = false;
            emit endrecording();
            int ndigs = QString("%1").arg(numframes).length();
            QString andigs = QString("%1").arg(ndigs);
            for (int i = 0 ; i < frames.count() ; i++){
                QString fileName = recordfilename+QString("%1.PNG").arg(i+1, ndigs, 10, QChar('0'));
                frames.at(i).save(fileName,"PNG");
            }
            qDebug() << "QT_VERSION:"
                << (((QT_VERSION) >> 16) & 0xff)
                << (((QT_VERSION) >> 8) & 0xff)
                << ((QT_VERSION) & 0xff);
#if QT_VERSION < 0x050F00
            QString make_movie = recordcommand + recordfilename
                    + QString("%"+andigs+"d.") + "PNG "
                    + recordfilename + QString(".mp4");
            int ierr = QProcess::execute(make_movie);
#else
            QStringList arguments = recordcommand.split(" ");
            arguments.removeAll(QString(""));
            QString command = arguments[0];
            arguments.removeFirst();
            arguments << recordfilename + QString("%"+andigs+"d.") + "PNG"
                      << recordfilename + QString(".mp4");
            int ierr = QProcess::execute(command,arguments);
#endif
            if (ierr != 0){
                QMessageBox msgBox;
                msgBox.setText(tr("create_film"));
                msgBox.setInformativeText(tr("It could not make film"));
                msgBox.setIcon(QMessageBox::Warning);
                msgBox.exec();
                emit endmakingmovie();
            }
            else{
                emit endmakingmovie();
                if(removeframes){
                    QProgressDialog progress("Making movie...", "Cancel", 0, frames.count(), nullpointer);
                    progress.setWindowModality(Qt::WindowModal);
                    progress.setMinimumDuration(100);
                    progress.show();
                    progress.setLabelText("Deleting frames...");
                    for (int i = 1 ; i < frames.count()+1 ; i++){
                        progress.setValue(i);
                        QString fname = recordfilename + QString("%1.PNG").arg(i, ndigs, 10, QChar('0'));
                        QFile file(fname);
                        bool removed = file.remove();
                        if (!removed){
                            QMessageBox msgBox;
                            msgBox.setText(tr("create_film"));
                            msgBox.setInformativeText(tr("It could not remove files with frames"));
                            msgBox.setIcon(QMessageBox::Warning);
                            msgBox.exec();
                            break;
                        }
                        if (progress.wasCanceled())
                            break;
                    }
                    progress.setValue(frames.count());
                    progress.close();
                }
            }
        }
    }
}

void glWindow::paintGLbuff(QSize size)
{
    QOpenGLFramebufferObjectFormat fboFormat;
    fboFormat.setSamples(16);
    fboFormat.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
    fbo = new QOpenGLFramebufferObject(size,fboFormat);
    saveGLState();
    fbo->bind();

    QOpenGLPaintDevice device(size);
    QPainter painterbuff(&device);

    mvp_list->clear();

    if (!molecules)
        return;

    painterbuff.beginNativePainting();

        // Clear color and depth buffer
        if (transparentbg)
            glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        else
            glClearColor(Bkg_Color.x(), Bkg_Color.y(), Bkg_Color.z(), 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glShadeModel( GL_SMOOTH );

        // Calculate aspect ratio
        qreal aspect = qreal(size.width()) / qreal(size.height() ? size.height() : 1);
        // Reset projection
        projection.setToIdentity();

        // Set perspective projection
        projection.perspective(fov, aspect, zNear, zFar);

        //  Model matrix (set to identity: all models are suppossed to be generated this way)
        QMatrix4x4 m_matrix;
        m_matrix.setToIdentity();

        QVector3D LightPosition_worldspace = lightPosition;

        program.bind();

        for (int i = 0 ; i < molecules->count() ; i++){

            //  View matrix
            QMatrix4x4 v_matrix;
            v_matrix.translate(world_translation);
            v_matrix.rotate(getrotation());
            v_matrix.translate(molecules->at(i)->gettranslation());
            v_matrix.rotate(molecules->at(i)->getrotation());

            molecules->at(i)->setrotationButtons();
            molecules->at(i)->settranslationButtons();

            // Model-view transformation
            QMatrix4x4 mv_matrix;
            mv_matrix = v_matrix * m_matrix;

            // Modelview-projection matrix
            QMatrix4x4 mvp_matrix;
            mvp_matrix = projection * mv_matrix;
            mvp_list->append(mvp_matrix);

            // Set model matrix
            program.setUniformValue("m_matrix", m_matrix);

            // Set view matrix
            program.setUniformValue("v_matrix", v_matrix);

            // Set model-view transformation
            program.setUniformValue("mv_matrix", mv_matrix);

            // Set modelview-projection matrix
            program.setUniformValue("mvp_matrix", mvp_matrix);

            // Set lights properties
            program.setUniformValue("LightPosition_worldspace", LightPosition_worldspace);
            program.setUniformValue("Ambient_Color", Ambient_Color);
            program.setUniformValue("Light_Color", Light_Color);
            program.setUniformValue("Light_Power", Light_Power);
            program.setUniformValue("Specular_Color", Specular_Color);
            program.setUniformValue("Specular_Index", Specular_Index);
            program.setUniformValue("Linear_Attenuation", linearattenuation);

            geometries->drawStructure(&program, molecules, i);
        }

        if (showdihedrals && measures){
            makedihedralplanes();
            if (dihedralindices->length() > 3){
                QMatrix4x4 v_matrix;
                v_matrix.setToIdentity();
                v_matrix.rotate(getrotation());
                QMatrix4x4 mv_matrix;
                mv_matrix.setToIdentity();
                QMatrix4x4 mvp_matrix;
                mvp_matrix = projection;
                program.setUniformValue("m_matrix", m_matrix);
                program.setUniformValue("v_matrix", v_matrix);
                program.setUniformValue("mv_matrix", mv_matrix);
                program.setUniformValue("mvp_matrix", mvp_matrix);

                geometries->drawDihedrals(&program, dihedralindices, dihedralvertices);
            }
        }

        if (showaxes){
            QMatrix4x4 v_matrix;
            v_matrix.translate(world_translation);
            v_matrix.rotate(getrotation());
            // Model-view transformation
            QMatrix4x4 mv_matrix;
            mv_matrix = v_matrix * m_matrix;
            // Modelview-projection matrix
            QMatrix4x4 mvp_matrix;
            mvp_matrix = projection * mv_matrix;
            mvp_list->append(mvp_matrix);
            program.setUniformValue("m_matrix", m_matrix);
            program.setUniformValue("v_matrix", v_matrix);
            program.setUniformValue("mv_matrix", mv_matrix);
            program.setUniformValue("mvp_matrix", mvp_matrix);

            geometries->drawLabAxes(&program, allaxesindices, allaxesvertices);
        }

        program.release();

        GLint vport[4];

        glGetIntegerv (GL_VIEWPORT, vport);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // This is necessary to keep the atom symbols on top even when wired surface is chosen
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);

    painterbuff.endNativePainting();

    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    drawlabels(&painterbuff, viewport);
    if (measures){
        drawmeasures(&painterbuff, viewport);
    }
    if (displayEPIC){
        drawEPIC(&painterbuff, viewport);
    }
    painterbuff.end();
    fbo->release();
    restoreGLState();

    QImage tempImage(size,QImage::Format_ARGB32);
    tempImage = fbo->toImage();

    if (!imagefilename.isEmpty()){
        QString filter=QFileInfo(imagefilename).suffix();
        const char * tipo;
        if (filter.toUpper() == "PNG") tipo = "PNG";
        else if (filter.toUpper() == "JPG") tipo = "JPG";
        else if (filter.toUpper() == "BMP") tipo = "BMP";
        else if (filter.toUpper() == "JPEG") tipo = "JPEG";
        else if (filter.toUpper() == "PPM") tipo = "PPM";
        else if (filter.toUpper() == "XBM") tipo = "XBM";
        else if (filter.toUpper() == "XPM") tipo = "XPM";
        else if (filter.toUpper() == "TIFF") tipo = "TIFF";
        else tipo = nullpointer; // Default: tries to imagine the format
        if( !tempImage.save( imagefilename, tipo, imagequality) ){
            QMessageBox msgBox;
            msgBox.setText(tr("paintGLbuff"));
            msgBox.setInformativeText(tr("Failed to save the image file ")+imagefilename);
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
        }
    }
    delete fbo;
    fbo = nullpointer;
}

void glWindow::checkbuffers(){
    if (geometries)
        geometries->checkbuffers();
}

void glWindow::reloadbuffers(){
    if (geometries)
        geometries->loadbuffers(molecules);
}

// Functions for drawing scenario and related functions
// =======================================================

//  Function drawlabels: draws the atom and cp labels
//

void glWindow::drawlabels(QPainter *painter, QRect viewport){
    painter->setBackgroundMode(Qt::TransparentMode);
    for (int i = 0 ; i < molecules->count() ; i++){
        if (molecules->at(i)->isvisible()){
            if (molecules->at(i)->getdrawatomcoords() || molecules->at(i)->getdrawatomindices()
                    || molecules->at(i)->getdrawatomsymbols()){
                drawatomsvalues(painter, viewport, i);
            }
            if (molecules->at(i)->isaxeslabelsvisible()){
                drawaxeslabels(painter, viewport, i);
            }
            if (molecules->at(i)->cps && (molecules->at(i)->cps->getdrawcpscoords() ||
                        molecules->at(i)->cps->getdrawcpsindices() || molecules->at(i)->cps->getdrawcpssymbols() ||
                        molecules->at(i)->cps->getdrawcpsvalues() )){
                drawcplabels(painter, viewport, i);
            }
            if (molecules->at(i)->surfaces){
                drawextremalabels(painter, viewport, i);
            }
        }
    }
    if (axeslabelsvisible){
        drawlabaxeslabels(painter, viewport);
    }
}
//  End of function drawlabels
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function drawatomsvalues: draws the atom indices
//

void glWindow::drawatomsvalues(QPainter *painter, QRect viewport, int i){
    if (molecules->at(i)->gethideatoms() && molecules->at(i)->gethidebonds())
        return;
    QVector <QVector3D> xyz = molecules->at(i)->getxyz();
    QVector <int> charges = molecules->at(i)->getcharges();
    QFontMetrics fm(molecules->at(i)->getfont());
    painter->setPen(molecules->at(i)->getfontcolor());
    painter->setFont(molecules->at(i)->getfont());
    for (int j = 0; j < xyz.length() ; j++){
        if ((molecules->at(i)->getonlyatomactive() && !(molecules->at(i)->isatomactive(j))) ||
              (molecules->at(i)->gethidehydrogens() && (molecules->at(i)->getznuc(j) == 1))
                ) continue;
        QVector4D vaux = mvp_list->at(i) * QVector4D(xyz[j], 1.0f);
        vaux /= vaux.w();
        vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
        vaux.setX(vaux.x()*viewport.width()+viewport.x());
        vaux.setY(vaux.y()*viewport.height()+viewport.y());
        QVector3D win = vaux.toVector3D();
        if (win.z() >= 1.f)
            continue;
        QString string;
        if (molecules->at(i)->getdrawatomsymbols())
            string = elem->getSymbol(charges[j]);
        if (molecules->at(i)->getdrawatomindices())
            string.append(QString("%1").arg(j+1));
        QSize fmsize = fm.size( Qt::TextSingleLine, string );
        if (molecules->at(i)->getdrawatomcoords()){
            QVector <QVector3D> xyz = molecules->at(i)->getxyz();
            double funits = 1.;
            if (molecules->at(i)->getangstromcoor()){
                funits = BOHR_TO_ANGSTROM;
            }
            if (string.length() > 0){
                string.append("\n");
                string.append(QString("(%1,%2,%3)")
                    .arg( QString::number(funits*xyz.at(j).x(), 'g', molecules->at(i)->getcoordprecision()))
                    .arg( QString::number(funits*xyz.at(j).y(), 'g', molecules->at(i)->getcoordprecision()))
                    .arg( QString::number(funits*xyz.at(j).z(), 'g', molecules->at(i)->getcoordprecision()))
                        );
                fmsize = fm.size( Qt::TextWordWrap, string );
            }
            else{
                string.append(QString("(%1,%2,%3)")
                    .arg( QString::number(funits*xyz.at(j).x(), 'g', molecules->at(i)->getcoordprecision()))
                    .arg( QString::number(funits*xyz.at(j).y(), 'g', molecules->at(i)->getcoordprecision()))
                    .arg( QString::number(funits*xyz.at(j).z(), 'g', molecules->at(i)->getcoordprecision()))
                        );
                fmsize = fm.size( Qt::TextSingleLine, string );
            }
        }
        QRect rect = QRect(win.x()-0.5f*fmsize.width(),
                viewport.height()-(win.y()+0.5f*fmsize.height())-molecules->at(i)->getlabelsvshift(),
                fmsize.width(), fmsize.height()) ;
        painter->drawText(rect, Qt::AlignCenter, string);
    }
}

//  Function drawaxeslabels: draws the axes labels
//
void glWindow::drawaxeslabels(QPainter *painter, QRect viewport, int i){
    fontaxeslabels = molecules->at(i)->getfontaxeslabels();
    QFontMetrics fm(fontaxeslabels);
    painter->setFont(fontaxeslabels);
    QVector <QVector3D> labelspositions = molecules->at(i)->getpositionaxeslabels();
//    X label
    painter->setPen(Xaxis_color);
    QVector4D vaux = mvp_list->at(i) * QVector4D(labelspositions.at(0),1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    QVector3D win = vaux.toVector3D();
    QString string("x");
    QSize fmsize = fm.size( Qt::TextSingleLine, string );
    QRect rect = QRect(win.x()-0.5f*fmsize.width(),
            viewport.height()-(win.y()+0.5f*fmsize.height()),
            fmsize.width(), fmsize.height()) ;
    if (win.z() < 1.f)
        painter->drawText(rect, Qt::AlignCenter, string);
//    Y label
    painter->setPen(Yaxis_color);
    vaux = mvp_list->at(i) * QVector4D(labelspositions.at(1),1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    win = vaux.toVector3D();
    string = QString("y");
    fmsize = fm.size( Qt::TextSingleLine, string );
    rect = QRect(win.x()-0.5f*fmsize.width(),
            viewport.height()-(win.y()+0.5f*fmsize.height()),
            fmsize.width(), fmsize.height()) ;
    if (win.z() < 1.f)
        painter->drawText(rect, Qt::AlignCenter, string);
//    Z label
    painter->setPen(Zaxis_color);
    vaux = mvp_list->at(i) * QVector4D(labelspositions.at(2),1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    win = vaux.toVector3D();
    string = QString("z");
    fmsize = fm.size( Qt::TextSingleLine, string );
    rect = QRect(win.x()-0.5f*fmsize.width(),
            viewport.height()-(win.y()+0.5f*fmsize.height()),
            fmsize.width(), fmsize.height()) ;
    if (win.z() < 1.f)
        painter->drawText(rect, Qt::AlignCenter, string);
}

//  Function drawcplabels: draws the critical points labels
//

void glWindow::drawcplabels(QPainter *painter, QRect viewport, int i){
    QFontMetrics fm(molecules->at(i)->cps->getfont());
    painter->setPen(molecules->at(i)->cps->getfontcolor());
    painter->setFont(molecules->at(i)->cps->getfont());
    QStringList CPtypes;
    CPtypes << "x" << "y" << "z" << "m";
    for (int type = 0 ; type < MAX_CPS ; type++){
        if (!molecules->at(i)->cps->isvisiblecps(type))
            continue;
        for (int j = 0 ; j < molecules->at(i)->cps->cpsxyzval[type].length() ; j++){
            if (molecules->at(i)->cps->getonlycpsactive() && !molecules->at(i)->cps->getcpsactive(type,j))
                continue;
            QVector4D vaux = mvp_list->at(i) * QVector4D(QVector3D(molecules->at(i)->cps->cpsxyzval[type].at(j)),1);
            vaux /= vaux.w();
            vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
            vaux.setX(vaux.x()*viewport.width()+viewport.x());
            vaux.setY(vaux.y()*viewport.height()+viewport.y());
            QVector3D win = vaux.toVector3D();
            if (win.z() >= 1.f)
                continue;
            QString string;
            if (molecules->at(i)->cps->getdrawcpssymbols())
                string = CPtypes[type];
            if (molecules->at(i)->cps->getdrawcpsindices())
                string.append(QString("%1").arg(j+1));
            QSize fmsize = fm.size( Qt::TextSingleLine, string );
            if (molecules->at(i)->cps->getdrawcpsvalues()){
                if (string.length() > 0){
                    string.append("\n");
                    string.append(QString::number( molecules->at(i)->cps->cpsxyzval[type].at(j).w(), 'g',
                            molecules->at(i)->cps->getcpprecision()));
                    fmsize = fm.size( Qt::TextWordWrap, string );
                }
                else{
                    string = QString::number( molecules->at(i)->cps->cpsxyzval[type].at(j).w(), 'g',
                            molecules->at(i)->cps->getcpprecision() );
                    fmsize = fm.size( Qt::TextSingleLine, string );
                }
            }
            if (molecules->at(i)->cps->getdrawcpscoords()){
                QVector4D xyzcp = molecules->at(i)->cps->cpsxyzval[type].at(j);
                double funits = 1.;
                if (molecules->at(i)->getangstromcp()){
                    funits = BOHR_TO_ANGSTROM;
                }
                if (string.length() > 0){
                    string.append("\n");
                    string.append(QString("(%1,%2,%3)")
                        .arg( QString::number(funits*xyzcp.x(), 'g',
                            molecules->at(i)->cps->getcpcoordprecision()))
                        .arg( QString::number(funits*xyzcp.y(), 'g',
                            molecules->at(i)->cps->getcpcoordprecision()))
                        .arg( QString::number(funits*xyzcp.z(), 'g',
                            molecules->at(i)->cps->getcpcoordprecision()))
                            );
                    fmsize = fm.size( Qt::TextWordWrap, string );
                }
                else{
                    string.append(QString("(%1,%2,%3)")
                      .arg( QString::number(funits*xyzcp.x(), 'g',
                          molecules->at(i)->cps->getcpcoordprecision()))
                      .arg( QString::number(funits*xyzcp.y(), 'g',
                          molecules->at(i)->cps->getcpcoordprecision()))
                      .arg( QString::number(funits*xyzcp.z(), 'g',
                          molecules->at(i)->cps->getcpcoordprecision()))
                            );
                    fmsize = fm.size( Qt::TextSingleLine, string );
                }
            }
            QRect rect = QRect(win.x()-0.5f*fmsize.width(),
                    viewport.height()-(win.y()+0.5f*fmsize.height())-molecules->at(i)->cps->getcpvshift(),
                    fmsize.width(), fmsize.height()) ;
            painter->drawText(rect, Qt::AlignCenter, string);
        }
    }
}

//  Function drawextremalabels: draws the local extrema labels
//

void glWindow::drawextremalabels(QPainter *painter, QRect viewport, int i){
    for (int j = 0 ; j < molecules->at(i)->surfaces->count() ; j++){
        if (molecules->at(i)->surfaces->at(j)->getvisible() &&
                molecules->at(i)->surfaces->at(j)->getshowlocalmax() &&
                (molecules->at(i)->surfaces->at(j)->getshowextremasymbols() ||
                 molecules->at(i)->surfaces->at(j)->getshowextremaindices() ||
                 molecules->at(i)->surfaces->at(j)->getshowextremacoords() ||
                 molecules->at(i)->surfaces->at(j)->getshowextremavalues())){
            QFontMetrics fm(molecules->at(i)->surfaces->at(j)->getfont());
            painter->setPen(molecules->at(i)->surfaces->at(j)->getfontcolor());
            painter->setFont(molecules->at(i)->surfaces->at(j)->getfont());
            for (int k = 0 ; k < molecules->at(i)->surfaces->at(j)->localextrema[0].length() ; k++){
                if (molecules->at(i)->surfaces->at(j)->getextremhidden(0,k))
                    continue;
                if (molecules->at(i)->surfaces->at(j)->getonlyextremactive() &&
                        !molecules->at(i)->surfaces->at(j)->getextremactive(0,k))
                    continue;
                QVector4D vaux = mvp_list->at(i) *
                        QVector4D(QVector3D(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k)),1);
                vaux /= vaux.w();
                vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
                vaux.setX(vaux.x()*viewport.width()+viewport.x());
                vaux.setY(vaux.y()*viewport.height()+viewport.y());
                QVector3D win = vaux.toVector3D();
                if (win.z() >= 1.f)
                    continue;
                QString string = "";
                QSize fmsize = fm.size( Qt::TextSingleLine, string );
                if (molecules->at(i)->surfaces->at(j)->getshowextremasymbols()){
                    string.append("M");
                    fmsize = fm.size( Qt::TextSingleLine, string );
                }
                if (molecules->at(i)->surfaces->at(j)->getshowextremaindices()){
                    string.append(QString("%1").arg(k+1));
                    fmsize = fm.size( Qt::TextSingleLine, string );
                }
                if (molecules->at(i)->surfaces->at(j)->getshowextremavalues()){
                    if (string.length() > 0){
                        string.append("\n ");
                        string.append(QString::number( molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).w(), 'g',
                                molecules->at(i)->surfaces->at(j)->getvalueprecision())+" ");
                        fmsize = fm.size( Qt::TextWordWrap, string );
                    }
                    else{
                        string.append(QString::number( molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).w(), 'g',
                                molecules->at(i)->surfaces->at(j)->getvalueprecision()));
                        fmsize = fm.size( Qt::TextSingleLine, string );
                    }
                }
                if (molecules->at(i)->surfaces->at(j)->getshowextremacoords()){
                    if (string.length() > 0){
                        string.append("\n");
                        string.append(QString("(%1,%2,%3)")
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).x(), 'g',
                                molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).y(), 'g',
                                    molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).z(), 'g',
                                        molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                                );
                        fmsize = fm.size( Qt::TextWordWrap, string );
                    }
                    else{
                        string.append(QString("(%1,%2,%3)")
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).x(), 'g',
                                molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).y(), 'g',
                                    molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[0].at(k).z(), 'g',
                                        molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                                );
                        fmsize = fm.size( Qt::TextSingleLine, string );
                    }
                }
                QRect rect = QRect(win.x()-0.5f*fmsize.width(),
                        viewport.height()-(win.y()+0.5f*fmsize.height())-molecules->at(i)->surfaces->at(j)->getvshift(),
                        fmsize.width(), fmsize.height()) ;
                painter->drawText(rect, Qt::AlignCenter, string);
            }
        }
        if (molecules->at(i)->surfaces->at(j)->getvisible() &&
                molecules->at(i)->surfaces->at(j)->getshowlocalmin() &&
                (molecules->at(i)->surfaces->at(j)->getshowextremasymbols() ||
                 molecules->at(i)->surfaces->at(j)->getshowextremaindices() ||
                 molecules->at(i)->surfaces->at(j)->getshowextremacoords() ||
                 molecules->at(i)->surfaces->at(j)->getshowextremavalues())){
            QFontMetrics fm(molecules->at(i)->surfaces->at(j)->getfont());
            painter->setPen(molecules->at(i)->surfaces->at(j)->getfontcolor());
            painter->setFont(molecules->at(i)->surfaces->at(j)->getfont());
            for (int k = 0 ; k < molecules->at(i)->surfaces->at(j)->localextrema[1].length() ; k++){
                if (molecules->at(i)->surfaces->at(j)->getextremhidden(1,k))
                    continue;
                if (molecules->at(i)->surfaces->at(j)->getonlyextremactive() &&
                        !molecules->at(i)->surfaces->at(j)->getextremactive(1,k))
                    continue;
                QVector4D vaux = mvp_list->at(i) *
                        QVector4D(QVector3D(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k)),1);
                vaux /= vaux.w();
                vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
                vaux.setX(vaux.x()*viewport.width()+viewport.x());
                vaux.setY(vaux.y()*viewport.height()+viewport.y());
                QVector3D win = vaux.toVector3D();
                if (win.z() >= 1.f)
                    continue;
                QString string = "";
                QSize fmsize = fm.size( Qt::TextSingleLine, string );
                if (molecules->at(i)->surfaces->at(j)->getshowextremasymbols()){
                    string.append("m");
                    fmsize = fm.size( Qt::TextSingleLine, string );
                }
                if (molecules->at(i)->surfaces->at(j)->getshowextremaindices()){
                    string.append(QString("%1").arg(k+1));
                    fmsize = fm.size( Qt::TextSingleLine, string );
                }
                if (molecules->at(i)->surfaces->at(j)->getshowextremavalues()){
                    if (string.length() > 0){
                        string.append("\n");
                        string.append(QString::number( molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).w(), 'g',
                                molecules->at(i)->surfaces->at(j)->getvalueprecision())+" ");
                        fmsize = fm.size( Qt::TextWordWrap, string );
                    }
                    else{
                        string.append(QString::number( molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).w(), 'g',
                                molecules->at(i)->surfaces->at(j)->getvalueprecision()));
                        fmsize = fm.size( Qt::TextSingleLine, string );
                    }
                }
                if (molecules->at(i)->surfaces->at(j)->getshowextremacoords()){
                    if (string.length() > 0){
                        string.append("\n");
                        string.append(QString("(%1,%2,%3)")
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).x(), 'g',
                                molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).y(), 'g',
                                    molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).z(), 'g',
                                        molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                                );
                        fmsize = fm.size( Qt::TextWordWrap, string );
                    }
                    else{
                        string.append(QString("(%1,%2,%3)")
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).x(), 'g',
                                molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).y(), 'g',
                                    molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                            .arg( QString::number(molecules->at(i)->surfaces->at(j)->localextrema[1].at(k).z(), 'g',
                                        molecules->at(i)->surfaces->at(j)->getcoordprecision()))
                                );
                        fmsize = fm.size( Qt::TextSingleLine, string );
                    }
                }
                QRect rect = QRect(win.x()-0.5f*fmsize.width(),
                        viewport.height()-(win.y()+0.5f*fmsize.height())-molecules->at(i)->surfaces->at(j)->getvshift(),
                        fmsize.width(), fmsize.height()) ;
                painter->drawText(rect, Qt::AlignCenter, string);
            }
        }
    }
}


//  Function drawaxeslabels: draws the Laboratory axes labels
//
void glWindow::drawlabaxeslabels(QPainter *painter, QRect viewport){
    QMatrix4x4 v_matrix;
    v_matrix.translate(world_translation);
    v_matrix.rotate(getrotation());
    // Modelview-projection matrix
    QMatrix4x4 mvp_matrix;
    mvp_matrix = projection * v_matrix;
    mvp_list->append(mvp_matrix);


    QFontMetrics fm(fontlabaxeslabels);
    painter->setFont(fontlabaxeslabels);
//    positionaxeslabels
//    X label
    painter->setPen(Xaxis_color);
    QVector4D vaux = mvp_matrix * QVector4D(positionaxeslabels->at(0),1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    QVector3D win = vaux.toVector3D();
    QString string("x");
    QSize fmsize = fm.size( Qt::TextSingleLine, string );
    QRect rect = QRect(win.x()-0.5f*fmsize.width(),
            viewport.height()-(win.y()+0.5f*fmsize.height()),
            fmsize.width(), fmsize.height()) ;
    if (win.z() < 1.f)
        painter->drawText(rect, Qt::AlignCenter, string);
//    Y label
    painter->setPen(Yaxis_color);
    vaux = mvp_matrix * QVector4D(positionaxeslabels->at(1),1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    win = vaux.toVector3D();
    string = QString("y");
    fmsize = fm.size( Qt::TextSingleLine, string );
    rect = QRect(win.x()-0.5f*fmsize.width(),
            viewport.height()-(win.y()+0.5f*fmsize.height()),
            fmsize.width(), fmsize.height()) ;
    if (win.z() < 1.f)
        painter->drawText(rect, Qt::AlignCenter, string);
//    Z label
    painter->setPen(Zaxis_color);
    vaux = mvp_matrix * QVector4D(positionaxeslabels->at(2),1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    win = vaux.toVector3D();
    string = QString("z");
    fmsize = fm.size( Qt::TextSingleLine, string );
    rect = QRect(win.x()-0.5f*fmsize.width(),
            viewport.height()-(win.y()+0.5f*fmsize.height()),
            fmsize.width(), fmsize.height()) ;
    if (win.z() < 1.f)
        painter->drawText(rect, Qt::AlignCenter, string);
}

//  Function drawmeasures: driving function for drawing angles, dihedral angles or distances
//
void glWindow::drawmeasures(QPainter *painter, QRect viewport){
    painter->setBackgroundMode(Qt::TransparentMode);

    if (showdihedrals){
        drawdihedrals(painter, viewport);
    }
    if (showangles){
        drawangles(painter, viewport);
    }
    if (showdistances){
        drawdistances(painter, viewport);
    }
    emit update_angles(anglecenters,v_list);
    emit update_dihedrals(dihedralcenters,v_list);
    emit update_distances(distancecenters,v_list);
}
//  End of function drawmeasures
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function drawEPIC: function for drawing EPIC energy
//
void glWindow::drawEPIC(QPainter *painter, QRect viewport){
    if (EPICposition == QPoint(0,0)){
        EPICposition = QPoint(viewport.width()/20,viewport.height()/20);
    }
    painter->setBackgroundMode(Qt::TransparentMode);
    painter->save();
    QPen pen(energycolor);
    painter->setPen(pen);
    painter->setFont(energyfont);
    QFontMetrics fm(energyfont);
    if (hartree){
        EPICstring = "EPIC Energy = "+QString::number( epicenergy, 'e', energyprecision-1)+" E"+QChar(0x2095);
    }
    else{
        EPICstring = "EPIC Energy = "+QString::number( epicenergy*HARTREETOKCALMOL, 'f', energyprecision-1)+" kcal/mol";
    }
    QSize fmsize = fm.size( Qt::TextSingleLine, EPICstring );
    QRect rect = QRect(EPICposition,fmsize) ;
    painter->drawText(rect,EPICstring);
    painter->restore();
}
//  End of function drawEPIC
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function drawangles: draw angles
//
void glWindow::drawangles(QPainter *painter, QRect viewport){
    painter->save();
    QPen pen(anglescolor, angleswidth, Qt::PenStyle(anglestype), Qt::RoundCap, Qt::RoundJoin);
    painter->setPen(pen);
    painter->setFont(anglesfont);
    QFontMetrics fm(anglesfont);
    for (int i = 0 ; i < anglecenters->length()-2 ; i +=3){
        if (anglecenters->at(i).molecule >= molecules->length() || anglecenters->at(i+1).molecule >= molecules->length()
                    || anglecenters->at(i+2).molecule >= molecules->length())
            continue;
        if (!molecules->at(anglecenters->at(i).molecule)->isvisible()) continue;
        QVector3D A3D = QVector3D(v_list->at(anglecenters->at(i).molecule)*QVector4D(anglecenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(v_list->at(anglecenters->at(i+1).molecule)*QVector4D(anglecenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(v_list->at(anglecenters->at(i+2).molecule)*QVector4D(anglecenters->at(i+2).xyz,1));
        double angle = 180. * acos(QVector3D::dotProduct((A3D-B3D).normalized(),(C3D-B3D).normalized())) / M_PI;
        QVector3D wincrntA = worldTocanvas(mvp_list->at(anglecenters->at(i).molecule), viewport, anglecenters->at(i).xyz);
        QVector3D wincrntB = worldTocanvas(mvp_list->at(anglecenters->at(i+1).molecule), viewport, anglecenters->at(i+1).xyz);
        QVector3D wincrntC = worldTocanvas(mvp_list->at(anglecenters->at(i+2).molecule), viewport, anglecenters->at(i+2).xyz);

        float dst = QVector3D(A3D-B3D).length();
        if (dst > disthres * (elem->getCovalentRadius(anglecenters->at(i).znuc)
                +elem->getCovalentRadius(anglecenters->at(i+1).znuc))){
            painter->drawLine(wincrntA.x(), viewport.height()-wincrntA.y(), wincrntB.x(), viewport.height()-wincrntB.y());
        }
        dst = QVector3D(C3D-B3D).length();
        if (dst > disthres * (elem->getCovalentRadius(anglecenters->at(i+2).znuc)
                +elem->getCovalentRadius(anglecenters->at(i+1).znuc))){
            painter->drawLine(wincrntC.x(), viewport.height()-wincrntC.y(), wincrntB.x(), viewport.height()-wincrntB.y());
        }

        if ((wincrntA-wincrntB).length() < arcradius || (wincrntC-wincrntB).length() < arcradius) continue;
        QRectF const rect(wincrntB.x() - arcradius, viewport.height()-(wincrntB.y() + arcradius), arcradius * 2, arcradius * 2);
        double ang1 = atan2(wincrntA.y() - wincrntB.y(), wincrntA.x() - wincrntB.x());
        double ang2 = atan2(wincrntC.y() - wincrntB.y(), wincrntC.x() - wincrntB.x());
        double angdif = ang2-ang1;
        QVector2D vtext = QVector2D((wincrntA-wincrntB).normalized() + (wincrntC-wincrntB).normalized());
        if (vtext.length() < 1.e-5){
            angdif = std::abs(angdif);
            vtext = QVector2D(-(wincrntA-wincrntB).y(),(wincrntA-wincrntB).x()).normalized();
        }
        else{
            if (angdif > M_PI)
                angdif = angdif - 2*M_PI;
            else if(angdif < -M_PI)
                angdif = 2*M_PI + angdif;
        }
        ang1 = 16 * ang1 * 180.0 / M_PI;
        double angspan = 16 * angdif * 180.0 / M_PI;
        if (drawarcs)
            painter->drawArc(rect, ang1, angspan);
        QString string;
        string = QString::number( std::abs(angle), 'g', anglesprecision);
        string.append(QChar(0260));
        QSize fmsize = fm.size( Qt::TextSingleLine, string );
        float vtextx = vtext.x();
        if (vtext.y() > 0.)
            vtext = QVector2D(wincrntB) + vtext.normalized() * arcradius + angleswidth * QVector2D(0,1);
        else
            vtext = QVector2D(wincrntB) + vtext.normalized() * arcradius + angleswidth * QVector2D(0,-1)
                    + 0.5f * QVector2D(0.,-fmsize.height());
        if (vtextx > 0.)
            vtext = vtext + angleswidth * QVector2D(1,0);
        else
            vtext = vtext - QVector2D(angleswidth+fmsize.width(),0);
        pen.setStyle(Qt::SolidLine);
        painter->setPen(pen);
        painter->drawText(vtext.x(),viewport.height()-(vtext.y()),string);
        pen.setStyle(Qt::PenStyle(anglestype));
        painter->setPen(pen);
    }
    painter->restore();
}
//  End of function drawangles
//

//  Function drawdihedrals: draw dihedrals
//
void glWindow::drawdihedrals(QPainter *painter, QRect viewport){
    painter->save();
    QPen pen;
    pen.setColor(dihedralscolor);
    painter->setPen(pen);
    painter->setFont(dihedralsfont);
    QFontMetrics fm(dihedralsfont);
    for (int i = 0 ; i < dihedralcenters->length()-3 ; i +=4){
        if (dihedralcenters->at(i).molecule >= molecules->length() || dihedralcenters->at(i+1).molecule >= molecules->length()
            || dihedralcenters->at(i+2).molecule >= molecules->length() || dihedralcenters->at(i+3).molecule >= molecules->length())
            continue;
        if (!molecules->at(dihedralcenters->at(i).molecule)->isvisible()) continue;
        QVector3D A3D = QVector3D(v_list->at(dihedralcenters->at(i).molecule)*QVector4D(dihedralcenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(v_list->at(dihedralcenters->at(i+1).molecule)*QVector4D(dihedralcenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(v_list->at(dihedralcenters->at(i+2).molecule)*QVector4D(dihedralcenters->at(i+2).xyz,1));
        QVector3D D3D = QVector3D(v_list->at(dihedralcenters->at(i+3).molecule)*QVector4D(dihedralcenters->at(i+3).xyz,1));
        QVector3D wincrntA = worldTocanvas(mvp_list->at(dihedralcenters->at(i).molecule), viewport, dihedralcenters->at(i).xyz);
        QVector3D wincrntB = worldTocanvas(mvp_list->at(dihedralcenters->at(i+1).molecule), viewport, dihedralcenters->at(i+1).xyz);
        QVector3D wincrntC = worldTocanvas(mvp_list->at(dihedralcenters->at(i+2).molecule), viewport, dihedralcenters->at(i+2).xyz);
        QVector3D wincrntD = worldTocanvas(mvp_list->at(dihedralcenters->at(i+3).molecule), viewport, dihedralcenters->at(i+3).xyz);
        QVector3D winstr = 0.25 * (wincrntA+wincrntB+wincrntC+wincrntD);
        QVector3D ABC = QVector3D::crossProduct((C3D-B3D),(A3D-B3D));
        QVector3D BCD = QVector3D::crossProduct((D3D-C3D),(B3D-C3D));
        double aux = QVector3D::dotProduct(ABC,BCD)/(ABC.length()*BCD.length());
        double angle;
        if (qAbs(qAbs(aux)-1.) > 1.e-6)
            angle = 180. * std::acos(aux) / M_PI;
        else{
            if (aux > 0.)
                angle = 0.;
            else
                angle = M_PI;
        }
        QString string;
        string = QString::number( std::abs(angle), 'g', dihedralsprecision);
        string.append(QChar(0260));
        QSize fmsize = fm.size( Qt::TextSingleLine, string );
        painter->drawText(QRect(winstr.x()-0.5f*fmsize.width(),
            viewport.height()-(winstr.y()+0.5f*fmsize.height()), fmsize.width(), fmsize.height()),
            string);
    }
    painter->restore();
}
//  End of function drawdihedrals
//

//  Function drawdistances: draw distances
//
void glWindow::drawdistances(QPainter *painter, QRect viewport){
    painter->save();
    QPen pen(distancescolor, lineswidth, Qt::PenStyle(linestype), Qt::RoundCap, Qt::RoundJoin);
    QColor *bkgcolor = new QColor(getBkgColor());
    painter->setPen(pen);
    painter->setFont(distancesfont);
    QFontMetrics fm(distancesfont);
    QString units("");
    if (angstrom)
        units.append(QChar(0x212B));
    else
        units.append(QString("a")+QChar(0x2080));
    for (int i = 0 ; i < distancecenters->length()-1 ; i +=2){
        if (distancecenters->at(i).molecule >= molecules->length() ||
                distancecenters->at(i+1).molecule >= molecules->length())
            continue;
        if (!molecules->at(distancecenters->at(i).molecule)->isvisible()
                || !molecules->at(distancecenters->at(i+1).molecule)->isvisible())
            continue;
        float dst = QVector3D(v_list->at(distancecenters->at(i+1).molecule)*QVector4D(distancecenters->at(i+1).xyz,1)
                    -v_list->at(distancecenters->at(i).molecule)*QVector4D(distancecenters->at(i).xyz,1)).length();
        QVector3D wincrntA = worldTocanvas(mvp_list->at(distancecenters->at(i+1).molecule), viewport, distancecenters->at(i+1).xyz);
        QVector3D wincrntB = worldTocanvas(mvp_list->at(distancecenters->at(i).molecule), viewport, distancecenters->at(i).xyz);
        QVector3D winstr = 0.5 * (wincrntA+wincrntB);
        QString string;
        if (angstrom)
            string = QString::number( dst * BOHR_TO_ANGSTROM, 'g', distprecision);
        else
            string = QString::number( dst, 'g', distprecision);
        string.append(units);
        QSize fmsize = fm.size( Qt::TextSingleLine, string );
        QRect rect = QRect(winstr.x()-0.5f*fmsize.width(),
                viewport.height()-(winstr.y()+0.5f*fmsize.height())-distvshift,
                fmsize.width(), fmsize.height()) ;

        if (drawlines && dst > disthres * (elem->getCovalentRadius(distancecenters->at(i+1).znuc)
                   +elem->getCovalentRadius(distancecenters->at(i).znuc))){
            painter->drawLine(wincrntA.x(), viewport.height()-wincrntA.y(), wincrntB.x(), viewport.height()-wincrntB.y());
        }
        if (!disttranspbkg)
            painter->fillRect(rect,QBrush(*bkgcolor));
        pen.setStyle(Qt::SolidLine);
        painter->setPen(pen);
        painter->drawText(rect, Qt::AlignCenter, string);
        pen.setStyle(Qt::PenStyle(linestype));
        painter->setPen(pen);
    }
    delete bkgcolor;
    painter->restore();
}
//  End of function drawdistances
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


void glWindow::emit_update_angles(){
    emit update_angles(anglecenters,v_list);
}

void glWindow::emit_update_dihedrals(){
    emit update_dihedrals(dihedralcenters,v_list);
}

void glWindow::emit_update_distances(){
    emit update_distances(distancecenters,v_list);
}

void glWindow::emit_update_worldtrotation(){
    emit update_worldrotation(world_rotation);
}

void glWindow::emit_update_worldtranslation(){
    emit update_worldtranslation(world_translation);
}

//  Function worldTocanvas: transforms point world coordinates to canvas coordinates
//
QVector3D glWindow::worldTocanvas(QMatrix4x4 mvp, QRect viewport, QVector3D point){
    QVector4D vaux = mvp * QVector4D(point,1);
    vaux /= vaux.w();
    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
    vaux.setX(vaux.x()*viewport.width()+viewport.x());
    vaux.setY(vaux.y()*viewport.height()+viewport.y());
    QVector3D win = vaux.toVector3D();
    return win;
}


//  Function getimagequality: retrieves image quality (compression)
//
int glWindow::getimagequality(){
    return imagequality;
}

//  Function getframeknt: retrieves frame counter
int glWindow::getframeknt(){
    return frameknt;
}


//  Function getLightPower: retrieves power of light source
//
float glWindow::getLightPower(){
    return Light_Power;
}

//  Function getrotation: retrieves quaternion  associated to world rotation
//
QQuaternion glWindow::getrotation(){
    return world_rotation;
}
//  Function getrotationAxis: retrieves axis  associated to world rotation
//
QVector3D glWindow::getrotationAxis(){
    return world_rotationAxis;
}

//  Function getSpecularIndex: retrieves material reflection index
//
float glWindow::getSpecularIndex(){
    return Specular_Index;
}

//  Function getAmbientColor: retrieves color of ambient light
//
QColor glWindow::getAmbientColor(){
    return QColor(255*Ambient_Color.x(),255*Ambient_Color.y(),255*Ambient_Color.z());
}

//  Function getBkgColor: retrieves background color
//
QColor glWindow::getBkgColor(){
    return QColor(255*Bkg_Color.x(),255*Bkg_Color.y(),255*Bkg_Color.z());
}

//  Function getLightColor: retrieves color of light source
//
//
QColor glWindow::getLightColor(){
    return QColor(255*Light_Color.x(),255*Light_Color.y(),255*Light_Color.z());
}

//  Function getlinearattenuation: retrieves light attenuation with distance (true: linear, false: quadratic)
//
bool glWindow::getlinearattenuation(){
    return linearattenuation;
}

//  Function getSpecularColor: retrieves material reflection color
//
QColor glWindow::getSpecularColor(){
    return QColor(255*Specular_Color.x(),255*Specular_Color.y(),255*Specular_Color.z());
}

//  Function getimagefilename: retrieves image file name
//
QString glWindow::getimagefilename(){
    return imagefilename;
}

//  Function getWindowName: retrieves window name
//
QString glWindow::getWindowName(){
    return windowname;
}

//  Function searchtatom: searchs for atom after double click
//
centerData glWindow::searchatom(int x,int y){
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    centerData center;
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!(molecules->at(i)->isvisible())) continue;
        QVector <QVector3D> xyz = molecules->at(i)->getxyz();
        int r2 = molecules->at(i)->getfont().pointSizeF() * 0.5;
        r2 *= r2;
        for (int j = 0; j < xyz.length() ; j++){
            QVector4D vaux = mvp_list->at(i) * QVector4D(xyz[j], 1.0f);
            vaux /= vaux.w();
            vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
            vaux.setX(vaux.x()*viewport.width()+viewport.x());
            vaux.setY(vaux.y()*viewport.height()+viewport.y());
            if (((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y())) <= r2){
                center.molecule = i;
                center.number= j;
                center.type = 0;  // 0 means it is an atom
                center.znuc = molecules->at(i)->getcharges().at(j);
                center.symbol = elem->getSymbol(center.znuc)+QString("%1").arg(j+1);
                center.x = x;
                center.y = y;
                QMatrix4x4 v_matrix;
                v_matrix.translate(molecules->at(center.molecule)->gettranslation());
                v_matrix.rotate(molecules->at(center.molecule)->getrotation());
                center.xyz = xyz[j];
                return center;
            }
        }
    }
    return center;
}

//  Function searchcps: searchs for critical points after double click
//
centerData glWindow::searchcps(int x,int y){
    QStringList CPstypes;
    CPstypes << "x" << "y" << "z" << "m";
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    centerData center;
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!molecules->at(i)->cps || !molecules->at(i)->isvisible())
            continue;
        for (int type = 0 ; type < MAX_CPS ; type++){
            int r2 = molecules->at(i)->getfont().pointSizeF() * 0.5;
            r2 *= r2;
            for (int j = 0; j < molecules->at(i)->cps->cpsxyzval[type].length() ; j++){
                QVector4D vaux = mvp_list->at(i) * QVector4D(QVector3D(molecules->at(i)->cps->cpsxyzval[type].at(j)),1);
                vaux /= vaux.w();
                vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
                vaux.setX(vaux.x()*viewport.width()+viewport.x());
                vaux.setY(vaux.y()*viewport.height()+viewport.y());
                if (((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y())) <= r2){
                    center.molecule = i;
                    center.number= j;
                    center.type = 1;  // 0 means it is a critical point
                    center.znuc = -1;
                    center.symbol = CPstypes.at(type)+QString("%1").arg(j+1);
                    center.x = x;
                    center.y = y;
                    center.xyz = QVector3D(molecules->at(i)->cps->cpsxyzval[type].at(j));
                    return center;
                }
            }
        }
    }
    return center;
}

//  Function searchangle: searchs for an angle label after double click
//
int glWindow::searchangle(int x,int y){
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    centerData center;

    QFontMetrics fm(anglesfont);
    QSize fmsize = fm.size( Qt::TextSingleLine, "XXXXX" );
    for (int i = 0 ; i < anglecenters->length()-2 ; i +=3){
        if (anglecenters->at(i).molecule >= molecules->length() || anglecenters->at(i+1).molecule >= molecules->length()
                    || anglecenters->at(i+2).molecule >= molecules->length())
            continue;
        if (!molecules->at(anglecenters->at(i).molecule)->isvisible()) continue;
        QVector3D wincrntA = worldTocanvas(mvp_list->at(anglecenters->at(i).molecule), viewport, anglecenters->at(i).xyz);
        QVector3D wincrntB = worldTocanvas(mvp_list->at(anglecenters->at(i+1).molecule), viewport, anglecenters->at(i+1).xyz);
        QVector3D wincrntC = worldTocanvas(mvp_list->at(anglecenters->at(i+2).molecule), viewport, anglecenters->at(i+2).xyz);

        if ((wincrntA-wincrntB).length() < arcradius || (wincrntC-wincrntB).length() < arcradius) continue;
        QVector2D vtext = QVector2D((wincrntA-wincrntB).normalized() + (wincrntC-wincrntB).normalized());
        QRect rect;
        if (vtext.length() < 1.e-5){
            vtext = QVector2D(-(wincrntA-wincrntB).y(),(wincrntA-wincrntB).x()).normalized();
        }
        float vtextx = vtext.x();
        if (vtext.y() > 0.)
            vtext = QVector2D(wincrntB) + vtext.normalized() * arcradius + angleswidth * QVector2D(0,1);
        else
            vtext = QVector2D(wincrntB) + vtext.normalized() * arcradius + angleswidth * QVector2D(0,-1)
                    + 0.5f * QVector2D(0.,-fmsize.height());
        if (vtextx > 0.){
            vtext = vtext + angleswidth * QVector2D(1,0);
            rect = QRect(vtext.x(),
                    viewport.height()-(vtext.y()+0.5f*fmsize.height()),
                    fmsize.width(), fmsize.height());
        }
        else{
            vtext = vtext - QVector2D(angleswidth+fmsize.width(),0);
            rect = QRect(vtext.x(),
                    viewport.height()-(vtext.y()+0.5f*fmsize.height()),
                    fmsize.width(), fmsize.height());
        }
        if (rect.contains(x,y)){
            anglecenters->remove(i+2);
            anglecenters->remove(i+1);
            anglecenters->remove(i);
            return i;
        }
    }
    return -1;
}
//  End of function searchangle
//

//  Function searchdihedral: searchs for a dihedral angle label after double click
//
int glWindow::searchdihedral(int x,int y){
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);

    QFontMetrics fm(dihedralsfont);
    QSize fmsize = fm.size( Qt::TextSingleLine, "XXXXX" );
    for (int i = 0 ; i < dihedralcenters->length()-3 ; i +=4){
        if (dihedralcenters->at(i).molecule >= molecules->length() || dihedralcenters->at(i+1).molecule >= molecules->length()
            || dihedralcenters->at(i+2).molecule >= molecules->length() || dihedralcenters->at(i+3).molecule >= molecules->length())
            continue;
        if (!molecules->at(dihedralcenters->at(i).molecule)->isvisible()) continue;
        QVector3D wincrntA = worldTocanvas(mvp_list->at(dihedralcenters->at(i).molecule), viewport, dihedralcenters->at(i).xyz);
        QVector3D wincrntB = worldTocanvas(mvp_list->at(dihedralcenters->at(i+1).molecule), viewport, dihedralcenters->at(i+1).xyz);
        QVector3D wincrntC = worldTocanvas(mvp_list->at(dihedralcenters->at(i+2).molecule), viewport, dihedralcenters->at(i+2).xyz);
        QVector3D wincrntD = worldTocanvas(mvp_list->at(dihedralcenters->at(i+3).molecule), viewport, dihedralcenters->at(i+3).xyz);
        QVector3D winstr = 0.25 * (wincrntA+wincrntB+wincrntC+wincrntD);
        QRect rect = QRect(winstr.x()-0.5f*fmsize.width(),
            viewport.height()-(winstr.y()+0.5f*fmsize.height()), fmsize.width(), fmsize.height());
        if (rect.contains(x,y)){
            dihedralcenters->remove(i+3);
            dihedralcenters->remove(i+2);
            dihedralcenters->remove(i+1);
            dihedralcenters->remove(i);
            return i;
        }
    }
    return -1;
}
//  End of function searchdihedral
//


//  Function searchdistance: searchs for distance label after double click
//
int glWindow::searchdistance(int x,int y){
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    centerData center;

    QFontMetrics fm(distancesfont);

    for (int i = 0 ; i < distancecenters->length()-1 ; i +=2){
        if (!molecules->at(distancecenters->at(i).molecule)->isvisible())
            continue;
        QVector3D winstr = 0.5*(worldTocanvas(mvp_list->at(distancecenters->at(i).molecule), viewport, distancecenters->at(i).xyz)
                    + worldTocanvas(mvp_list->at(distancecenters->at(i+1).molecule), viewport, distancecenters->at(i+1).xyz));
        QSize fmsize = fm.size( Qt::TextSingleLine, "xxxx" );
        QRect rect = QRect(winstr.x()-0.5f*fmsize.width(),
                viewport.height()-(winstr.y()+0.5f*fmsize.height())-distvshift,
                fmsize.width(), fmsize.height());
        if (rect.contains(x,y)){
            distancecenters->remove(i+1);
            distancecenters->remove(i);
            return i;
        }
    }
    return -1;
}
//  End of function searchdistance
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function searchextrema: searchs for mesp extrema after double click
//
centerData glWindow::searchextrema(int x,int y){
    QStringList extrematypes;
    extrematypes << "M" << "m";
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    centerData center;
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!molecules->at(i)->surfaces || !molecules->at(i)->isvisible())
            continue;
        for (int k = 0 ; k < molecules->at(i)->surfaces->count() ; k++){
            for (int type = 0 ; type < 2 ; type++){
                int r2 = molecules->at(i)->surfaces->at(k)->getfont().pointSizeF() * 0.5;
                r2 *= r2;
                for (int j = 0; j < molecules->at(i)->surfaces->at(k)->localextrema[type].count() ; j++){
                    QVector4D vaux = mvp_list->at(i) * QVector4D(molecules->at(i)->surfaces->at(k)->localextrema[type].at(j).toVector3D(),1);
                    vaux /= vaux.w();
                    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
                    vaux.setX(vaux.x()*viewport.width()+viewport.x());
                    vaux.setY(vaux.y()*viewport.height()+viewport.y());
                    if (((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y())) <= r2){
                        center.molecule = i;
                        center.number= j;
                        center.type = 2;  // 2 means it is a local extremum
                        center.znuc = -1;
                        center.symbol = extrematypes.at(type)+QString("%1").arg(j+1);
                        center.x = x;
                        center.y = y;
                        center.xyz = molecules->at(i)->surfaces->at(k)->localextrema[type].at(j).toVector3D();
                        return center;
                    }
                }
            }
        }
    }
    return center;
}

//  Function getLigthPosition: returns position of light source
//
QVector3D glWindow::getLigthPosition(){
    return lightPosition;
}

//  Function resetlabaxes: initializes laboratory axes
//
void glWindow::resetlabaxes(){
    axesarrowsize = 4;
    axesarrowwidth = 10;
    axeslength = 10;
    axesthickness = 4;
    axeslabelsvisible = false;
    Xaxis_color = QColor(0,255,0);
    Yaxis_color = QColor(0,0,255);
    Zaxis_color = QColor(255,0,0);
    fontaxeslabels = QFont("Noto Sans", 20, QFont::Bold);
    fontlabaxeslabels = QFont("Noto Sans", 20, QFont::Bold);
    make_axes();
    showaxes = false;
    update();
}

//  Function searchmolecule: activates or deactivates molecule
//
void glWindow::searchmolecule(int x,int y){
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    int distsqmin = 999999999;
    int distsq;
    int indexmin = 0;
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!molecules->at(i)->isvisible()) continue;
        QVector4D vaux = mvp_list->at(i) * QVector4D(0.f,0.f,0.f, 1.0f);
        vaux /= vaux.w();
        vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
        vaux.setX(vaux.x()*viewport.width()+viewport.x());
        vaux.setY(vaux.y()*viewport.height()+viewport.y());
        distsq = ((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y()));
        if (distsq < distsqmin){
            distsqmin = distsq;
            indexmin = i;
        }
    }
    emitcheckactivate(indexmin);
    return;
}


//  Function selectatom: activates or deactivates atom for label display
//
bool glWindow::selectatom(int x,int y){
    int isel = -1;
    int jsel = -1;
    qreal zmax = -zFar;
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!(molecules->at(i)->isvisible()))
            continue;
        if (isel >= 0)
            break;
        QVector <QVector3D> xyz = molecules->at(i)->getxyz();
        int r2 = molecules->at(i)->getfont().pointSizeF() * 0.5;
        r2 *= r2;
        for (int j = 0; j < xyz.length() ; j++){
            QVector4D vaux = mvp_list->at(i) * QVector4D(xyz[j], 1.0f);
            vaux /= vaux.w();
            vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
            vaux.setX(vaux.x()*viewport.width()+viewport.x());
            vaux.setY(vaux.y()*viewport.height()+viewport.y());
            if (((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y())) <= r2){
                if (isel == -1 || jsel == -1 || vaux.z() < zmax){
                    isel = i;
                    jsel = j;
                    zmax = vaux.z();
                    break;
                }
            }
        }
        if (isel != -1 && jsel != -1)
            break;
    }
    if (isel >= 0 && jsel >= 0){
        molecules->at(isel)->setatomactive(jsel,!molecules->at(isel)->isatomactive(jsel));
        return true;
    }
    else
        return false;
}

//  Function selectcps: activates or deactivates critical point for label display
//
bool glWindow::selectcps(int x,int y){
    int isel = -1;
    int typesel = -1;
    int jsel = -1;
    qreal zmax = -zFar;
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!molecules->at(i)->cps || !(molecules->at(i)->isvisible()))
            continue;
        if (isel >= 0)
            break;
        for (int type = 0 ; type < MAX_CPS ; type++){
            if (isel >= 0)
                break;
            int r2 = molecules->at(i)->getfont().pointSizeF() * 0.5;
            r2 *= r2;
            for (int j = 0; j < molecules->at(i)->cps->cpsxyzval[type].length() ; j++){
                QVector4D vaux = mvp_list->at(i) * QVector4D(QVector3D(molecules->at(i)->cps->cpsxyzval[type].at(j)),1);
                vaux /= vaux.w();
                vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
                vaux.setX(vaux.x()*viewport.width()+viewport.x());
                vaux.setY(vaux.y()*viewport.height()+viewport.y());
                if (((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y())) <= r2){
                    if (isel == -1 || jsel == -1 || typesel == -1 || vaux.z() < zmax){
                        isel = i;
                        typesel = type;
                        jsel = j;
                        zmax = vaux.z();
                        break;
                    }
                }
            }
        }
    }
    if (isel >= 0 && jsel >= 0 && typesel >= 0){
        molecules->at(isel)->cps->setcpsactive(typesel,jsel,!molecules->at(isel)->cps->getcpsactive(typesel,jsel));
        return true;
    }
    else
        return false;
}

//  Function selectmespmaximum: activates or deactivates MESP local extrema for label display
//
bool glWindow::selectmespextrema(int x,int y){
    int isel = -1;
    int typesel = -1;
    int jsel = -1;
    int ksel = -1;
    qreal zmax = -zFar;
    GLint vport[4];
    glGetIntegerv (GL_VIEWPORT, vport);
    QRect viewport = QRect(vport[0],vport[1],vport[2],vport[3]);
    for (int i = 0 ; i < molecules->count() ; i++){
        if (!(molecules->at(i)->surfaces) || !(molecules->at(i)->isvisible()))
            continue;
        if (isel >= 0)
            break;
        for (int k = 0 ; k < molecules->at(i)->surfaces->count() ; k++){
            if (!molecules->at(i)->surfaces->at(k)->getvisible())
                continue;
            if (isel >= 0)
                break;
            for (int type = 0 ; type < 2 ; type++){
                if (isel >= 0)
                    break;
                int r2 = molecules->at(i)->surfaces->at(k)->getfont().pointSizeF() * 0.5;
                r2 *= r2;
                for (int j = 0; j < molecules->at(i)->surfaces->at(k)->localextrema[type].count() ; j++){
                    QVector4D vaux = mvp_list->at(i) * QVector4D(molecules->at(i)->surfaces->at(k)->localextrema[type].at(j).toVector3D(),1);
                    vaux /= vaux.w();
                    vaux = vaux * 0.5f + QVector4D(0.5f, 0.5f, 0.5f, 0.5f);
                    vaux.setX(vaux.x()*viewport.width()+viewport.x());
                    vaux.setY(vaux.y()*viewport.height()+viewport.y());
                    if (((x-vaux.x())*(x-vaux.x())+(viewport.height()-y-vaux.y())*(viewport.height()-y-vaux.y())) <= r2){
                        if (isel == -1 || jsel == -1 || ksel == -1 || typesel == -1 || vaux.z() < zmax){
                            isel = i;
                            typesel = type;
                            jsel = j;
                            ksel = k;
                            zmax = vaux.z();
                            break;
                        }
                    }
                }
            }
        }
    }
    if (isel >= 0 && jsel >= 0 && ksel >= 0 && typesel >= 0){
        molecules->at(isel)->surfaces->at(ksel)->setextremactive(typesel,jsel,
                    !molecules->at(isel)->surfaces->at(ksel)->getextremactive(typesel,jsel));
        return true;
    }
    return false;
}

//  Function selectPopUpWindow: dialog window for choosing atoms to display or hide
//
bool glWindow::selectPopUpWindow(int x, int y){
    if (molecules->isEmpty()) return false;
    DialogSPB *select = new DialogSPB(molecules);
    select->setAttribute( Qt::WA_DeleteOnClose );
    select->move(x,y);
    select->exec();
    if (select->imolsel < 0 || select->imolsel >= molecules->length())
        return false;
    if (select->iatomsel >= 0 && select->iatomsel < molecules->at(select->imolsel)->getnumatoms()){
        molecules->at(select->imolsel)->setatomactive(select->iatomsel,
                !molecules->at(select->imolsel)->isatomactive(select->iatomsel));
    }
    else if (molecules->at(select->imolsel)->cps && select->icpsel >= 0
             && select->icpsel < molecules->at(select->imolsel)->getnumcps(select->icptypesel) ){
        molecules->at(select->imolsel)->cps->setcpsactive(select->icptypesel,select->icpsel,
                !molecules->at(select->imolsel)->cps->iscpsactive(select->icptypesel,select->icpsel) );
    }
    return true;
}

//  Function resetangles: clear list of centers for angles
//
void glWindow::resetangles(){
    anglecenters->clear();
    if (showangles) update();
}

//  Function make_axes: makes frame axes
void glWindow::make_axes(){
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    QQuaternion q;
    QVector3D scale = SCALE * QVector3D(axesthickness,axesthickness,MOLSCALEHEIGHT*axeslength);

    VertexNormalData v;
    positionaxeslabels->clear();
    allaxesindices->clear();
    allaxesvertices->clear();
    allaxesindicesoffset->clear();
    allaxesindicesoffset->append(0);
//    X axis
    uint maxindex = 0;
    int kshift = 0;
    for (int k = 0 ; k < cylinderindices.length() ; k++){
        allaxesindices->append(cylinderindices.at(k) + kshift);
        maxindex = std::max(maxindex,allaxesindices->last()+1);
    }
    allaxesindicesoffset->append(allaxesindices->length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(1.,0.,0.));
    color = QVector4D(Xaxis_color.redF(), Xaxis_color.greenF(), Xaxis_color.blueF(), 1.);
    for (int k = 0 ; k < cylindervertices.length() ; k++){
        position = q.rotatedVector( scale * cylindervertices.at(k));
        normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices->append(v);
    }

//    Y axis
    kshift = maxindex;
    for (int k = 0 ; k < cylinderindices.length() ; k++){
        allaxesindices->append(cylinderindices.at(k) + kshift);
        maxindex = std::max(maxindex,allaxesindices->last()+1);
    }
    allaxesindicesoffset->append(allaxesindices->length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,1.,0.));
    color = QVector4D(Yaxis_color.redF(), Yaxis_color.greenF(), Yaxis_color.blueF(), 1.);
    for (int k = 0 ; k < cylindervertices.length() ; k++){
        position = q.rotatedVector( scale * cylindervertices.at(k));
        normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices->append(v);
    }
//    Z axis
    kshift = maxindex;
    for (int k = 0 ; k < cylinderindices.length() ; k++){
        allaxesindices->append(cylinderindices.at(k) + kshift);
        maxindex = std::max(maxindex,allaxesindices->last()+1);
    }
    allaxesindicesoffset->append(allaxesindices->length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,0.,1.));
    color = QVector4D(Zaxis_color.redF(), Zaxis_color.greenF(), Zaxis_color.blueF(), 1.);
    for (int k = 0 ; k < cylindervertices.length() ; k++){
        position = q.rotatedVector( scale * cylindervertices.at(k));
        normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices->append(v);
    }

//      Arrows
    scale = SCALE * QVector3D(axesarrowwidth,axesarrowwidth,MOLSCALEARROWSHEIGHT*axesarrowsize);
    QVector3D shiftorig = SCALE * QVector3D(0,0,MOLSCALEHEIGHT*axeslength);
//    X axis arrow
    kshift = maxindex;
    for (int k = 0 ; k < coneindices.length() ; k++){
        allaxesindices->append(coneindices.at(k) + kshift);
        maxindex = std::max(maxindex,allaxesindices->last()+1);
    }
    allaxesindicesoffset->append(allaxesindices->length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(1.,0.,0.));
    color = QVector4D(Xaxis_color.redF(), Xaxis_color.greenF(), Xaxis_color.blueF(), 1.);
    for (int k = 0 ; k < conevertices.length() ; k++){
        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig);
        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices->append(v);
    }
    QFontMetrics fm(fontlabaxeslabels);
    QSize fmsize = fm.size( Qt::TextSingleLine, "x" );
    int shift = (int)std::sqrt((double)(fmsize.width()*fmsize.width())+(double)(fmsize.height()*fmsize.height()));
    positionaxeslabels->append(q.rotatedVector( QVector3D(0.,0.,
                    SCALE * (shift + 2. * MOLSCALEARROWSHEIGHT*axesarrowsize) ) + shiftorig));
//    Y axis arrow
    kshift = maxindex;
    for (int k = 0 ; k < coneindices.length() ; k++){
        allaxesindices->append(coneindices.at(k) + kshift);
        maxindex = std::max(maxindex,allaxesindices->last()+1);
    }
    allaxesindicesoffset->append(allaxesindices->length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,1.,0.));
    color = QVector4D(Yaxis_color.redF(), Yaxis_color.greenF(), Yaxis_color.blueF(), 1.);
    for (int k = 0 ; k < conevertices.length() ; k++){
        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig);
        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices->append(v);
    }
    fmsize = fm.size( Qt::TextSingleLine, "x" );
    shift = (int)std::sqrt((double)(fmsize.width()*fmsize.width())+(double)(fmsize.height()*fmsize.height()));
    positionaxeslabels->append(q.rotatedVector( QVector3D(0.,0.,
                    SCALE * (shift + 2. * MOLSCALEARROWSHEIGHT*axesarrowsize) ) + shiftorig));
//    Z axis arrow
    kshift = maxindex;
    for (int k = 0 ; k < coneindices.length() ; k++){
        allaxesindices->append(coneindices.at(k) + kshift);
        maxindex = std::max(maxindex,allaxesindices->last()+1);
    }
    allaxesindicesoffset->append(allaxesindices->length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,0.,1.));
    color = QVector4D(Zaxis_color.redF(), Zaxis_color.greenF(), Zaxis_color.blueF(), 1.);
    for (int k = 0 ; k < conevertices.length() ; k++){
        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig);
        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices->append(v);
    }
    fmsize = fm.size( Qt::TextSingleLine, "x" );
    shift = (int)std::sqrt((double)(fmsize.width()*fmsize.width())+(double)(fmsize.height()*fmsize.height()));
    positionaxeslabels->append(q.rotatedVector( QVector3D(0.,0.,
                    SCALE * (shift + 2. * MOLSCALEARROWSHEIGHT*axesarrowsize) ) + shiftorig));
}

// Function makeAxesCone: generates vertices for drawing a cone with a given number of "slices" (longitude)
//      and "stacks" (latitude), a radius equal to 1 and a "height"
//      Base centered at (0,0,0) top at (0,0,height)
//
//  Counterclockwise triangles generated
//
void glWindow::makeAxesCone(int slices, int stacks, float height){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines

    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = height / (stacks-1);

    conevertices.clear();
    conenormals.clear();

    // vertices on the main body
    float aux = 1. / stacks;
    float bux = 1 / height;
    for (int i = 0; i < stacks-1; i++) {
        float rad = 1. - i * aux;
        float h = deltaLat * i;
        for (int j = 0; j < slices; j++) {
            conevertices << QVector3D(rad*cos(deltaLong * j), rad*sin(deltaLong * j), h);
            conenormals << QVector3D(cos(deltaLong * j), sin(deltaLong * j), bux).normalized();
        }
    }
    conevertices << QVector3D(0,0, height);
    conenormals << QVector3D(0,0,1.);

    // vertices on the base
    int indexbase = conevertices.length();
    conevertices << QVector3D(0,0,0);
    conenormals << QVector3D(0,0,-1.);
    for (int j = 0; j < slices; j++) {
        conevertices << QVector3D(cos(deltaLong * j), sin(deltaLong * j), 0);
        conenormals << QVector3D(0,0,-1.);
    }

    // Generates the indices
    coneindices.clear();
    // indices of the body
    for (int i = 1; i < stacks-1 ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            coneindices << (i - 1) * slices + j;
            coneindices << (i - 1) * slices + j + 1;
            coneindices << i * slices + j;
            // triangle two
            coneindices << (i - 1) * slices + j + 1;
            coneindices << i * slices + j + 1;
            coneindices << i * slices + j;

        }
        // triangle one
        coneindices << i * slices-1;
        coneindices << (i-1) * slices;
        coneindices << (i+1) * slices - 1;
        // triangle two
        coneindices << (i-1) * slices;
        coneindices << i * slices;
        coneindices << (i+1) * slices - 1;
    }
//    triangles of cone vertex
    for (int j = 0; j < slices-1; j++) {
        coneindices << (stacks - 2) * slices + j;
        coneindices << (stacks - 2) * slices + j + 1;
        coneindices << (stacks - 1) * slices;
    }
    coneindices << (stacks - 1) * slices - 1;
    coneindices << (stacks - 2) * slices;
    coneindices << (stacks - 1) * slices;
    // indices of the base
    for (int j = 0; j < slices-1; j++) {
        coneindices << indexbase;
        coneindices << j + 2 + indexbase;
        coneindices << j + 1 + indexbase;
    }
    coneindices << indexbase;
    coneindices << 1 + indexbase;
    coneindices << slices + indexbase;
}

//  End of function makeAxesCone
//  ----------------------------------------------------------------------------------------------------------------------------


// Function makeAxesCylinder: generates vertices for drawing the axes cylinder side surface (neither top nor base) with a given number of
//      "slices" (longitude) and "stacks" (latitude), a radius equal to 1 and a height equal to 1 to be reescaled afterwards
//      as appropriate. Base center at (0,0,0). Orientation: (0,0,1)
//
//  Counterclockwise triangles generated
//
void glWindow::makeAxesCylinder(int slices, int stacks){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines
    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = 1. / (stacks-1);

    cylindervertices.clear();

    // vertices on the main body
    for (int i = 0; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
            cylindervertices << QVector3D(cos(deltaLong * j), sin(deltaLong * j), deltaLat * i);
        }
    }
    // Generates the indices
    cylinderindices.clear();
    // add body (no. of element is (stacks - 2) * slices * 6
    for (int i = 1; i < stacks ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            cylinderindices << (i - 1) * slices + j;
            cylinderindices << i * slices + j + 1;
            cylinderindices << i * slices + j;

            // triangle two
            cylinderindices << (i - 1) * slices + j;
            cylinderindices << (i - 1) * slices + j + 1;
            cylinderindices << i * slices + j + 1;


        }
        // triangle one
        cylinderindices << i * slices - 1;
        cylinderindices << i * slices;
        cylinderindices << (i+1) * slices - 1;

        // triangle two
        cylinderindices << i * slices - 1;
        cylinderindices << (i - 1) * slices;
        cylinderindices << i * slices;

    }
}

//  Function makedihedralplanes: build planes for dihedral angles
//
void glWindow::makedihedralplanes(){
    dihedralindices->clear();
    dihedralvertices->clear();
    uint indices = 0;
    VertexNormalData v;
    v.color = QVector4D(dihedralplanescolor.redF(),dihedralplanescolor.greenF(),
                        dihedralplanescolor.blueF(),dihedralplanescolor.alphaF());
    for (int i = 0 ; i < dihedralcenters->length()-3 ; i +=4){
        if (!(molecules->at(dihedralcenters->at(i).molecule)->isvisible()))
            continue;
//        Plane ABC
        QVector3D A3D = QVector3D(v_list->at(dihedralcenters->at(i).molecule)*QVector4D(dihedralcenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(v_list->at(dihedralcenters->at(i+1).molecule)*QVector4D(dihedralcenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(v_list->at(dihedralcenters->at(i+2).molecule)*QVector4D(dihedralcenters->at(i+2).xyz,1));
        QVector3D D3D = QVector3D(v_list->at(dihedralcenters->at(i+3).molecule)*QVector4D(dihedralcenters->at(i+3).xyz,1));
        v.normal = QVector3D::crossProduct(C3D-B3D,A3D-B3D);
        v.position = A3D;
        dihedralvertices->append(v);
        v.position = B3D;
        dihedralvertices->append(v);
        v.position = C3D;
        dihedralvertices->append(v);
        v.position = A3D-B3D+C3D;
        dihedralvertices->append(v);
//        First triangle of plane ABC
        dihedralindices->append(indices);
        dihedralindices->append(indices+1);
        dihedralindices->append(indices+2);
//        Second triangle of plane ABC
        dihedralindices->append(indices);
        dihedralindices->append(indices+2);
        dihedralindices->append(indices+3);
        indices += 4;
//        Plane BCD
        v.normal = QVector3D::crossProduct(B3D-C3D,D3D-C3D);
        v.position = B3D;
        dihedralvertices->append(v);
        v.position = C3D;
        dihedralvertices->append(v);
        v.position = D3D;
        dihedralvertices->append(v);
        v.position = C3D-B3D+D3D;
        dihedralvertices->append(v);
//        First triangle of plane BCD
        dihedralindices->append(indices);
        dihedralindices->append(indices+1);
        dihedralindices->append(indices+2);

//        Second triangle of plane BCD
        dihedralindices->append(indices+1);
        dihedralindices->append(indices+2);
        dihedralindices->append(indices+3);
        indices += 4;
    }
}
//  Function resetdihedrals: clear list of centers for dihedral angles
//
void glWindow::resetdihedrals(){
    dihedralcenters->clear();
    dihedralvertices->clear();
    dihedralindices->clear();
    update();
}

//  Function resetdistances: clear list of centers for distances
//
void glWindow::resetdistances(){
    distancecenters->clear();
    if (showdistances) update();
}


//  Function resetEPICpressed: sets EPICpressed flag to false
//
void glWindow::resetEPICpressed(){
    EPICpressed = false;
}

//  Function resetmeasures: clear list of centers for angles
//
void glWindow::resetmeasures(){
    anglecenters->clear();
    dihedralcenters->clear();
    anglecenters->clear();
    angles = false;
    dihedrals = false;
    distances = false;
    angstrom = false;
    update();
}

//  Function searchEPICenergy: checks if the cursor is on EPIC energy string
//
bool glWindow::searchEPICenergy(QVector2D pos){
    QFontMetrics fm(energyfont);
    QSize fmsize = fm.size( Qt::TextSingleLine, EPICstring );
    QRect rect = QRect(EPICposition,fmsize) ;
    return rect.contains(pos.x(),pos.y());
}

//  Function setAmbientColor: sets color of ambient light
//
void glWindow::setAmbientColor(QColor a){
    Ambient_Color = QVector3D(a.redF(),a.greenF(),a.blueF());;
}

//  Function setangstrom: set angstroms for distance measures
//
void glWindow::setangstrom(bool a){
    angstrom = a;
    if (showdistances) update();
}

//  Function setangles: set allow/deny angles measures
//
void glWindow::setangles(bool a){
    angles = a;
    if (angles){
        dihedrals = false;
        distances = false;
        update();
    }
}

//  Function setanglescolor: sets color for angles font
//
void glWindow::setanglescolor(QColor a){
    anglescolor = a;
    if (showangles) update();
}

//  Function setanglesfont: set angles font
void glWindow::setanglesfont(QFont a){
    anglesfont = a;
    if (showangles) update();
}

//  Function setanglesprecision: set number of decimals for angles
//
void glWindow::setanglesprecision(int a){
    anglesprecision = a;
    if(showangles) update();
}

//  Function setanglestype: sets type for angles lines (0: NoLine, 1: SolidLine, 2: DashLine, 3: DotLine,
//          4: DashDotLine, 5: DashDotDotLine
//
void glWindow::setanglestype(int a){
    anglestype = a;
    if (showangles && drawarcs) update();
}

//  Function setangleswidth: sets width for angles lines
//
void glWindow::setangleswidth(int a){
    angleswidth = a;
    if (showangles && drawarcs) update();
}

//  Function setarcradius: sets size of angle arc
//
void glWindow::setarcradius(int a){
    arcradius = a;
    if (showangles && drawarcs) update();
}

//  Function setaxesarrowsize: sets axes arrows size
//
void glWindow::setaxesarrowsize(int a){
    axesarrowsize = a;
    make_axes();
    if (showaxes) update();
}

//  Function setaxesarrowsize: sets axes arrows size
//
void glWindow::setaxesarrowwidth(int a){
    axesarrowwidth = a;
    make_axes();
    if (showaxes) update();
}

//  Function setaxeslength: sets axes length
//
void glWindow::setaxeslength(int a){
    axeslength = a;
    make_axes();
    if (showaxes) update();
}

//  Function setaxesthickness: sets axes thickness
//
void glWindow::setaxesthickness(int a){
    axesthickness = a;
    make_axes();
    if (showaxes) update();
}

//  Function setaxesvisible: toggles axes displey
//
void glWindow::setaxesvisible(bool a){
    showaxes = a;
    update();
}

//  Function setaxeslabelsvisible: toggles axes labels displey
//
void glWindow::setaxeslabelsvisible(bool a){
    axeslabelsvisible = a;
    if (showaxes) update();
}

//  Function setBkgColor: sets background color
//
void glWindow::setBkgColor(QColor a){
    Bkg_Color = QVector3D(a.redF(),a.greenF(),a.blueF());
}

//  Function setdihedrals: set allow/deny dihedral angles measures
//
void glWindow::setdihedrals(bool a){
    dihedrals = a;
    if (dihedrals){
        angles = false;
        distances = false;
        update();
    }
}

//  Function setdihedralscolor: sets color for dihedral font
//
void glWindow::setdihedralscolor(QColor a){
    dihedralscolor = a;
    if (showdihedrals) update();
}

//  Function setdihedralplanescolor: sets color for dihedral font
//
void glWindow::setdihedralplanescolor(QColor a){
    dihedralplanescolor = a;
    if (showdihedrals) update();
}

//  Function setdihedralsfont: set dihedral angles font
void glWindow::setdihedralsfont(QFont a){
    dihedralsfont = a;
    if (showdihedrals) update();
}

//  Function setdihedralsprecision: set number of decimals for dihedral angles
//
void glWindow::setdihedralsprecision(int a){
    dihedralsprecision = a;
    if(showdihedrals) update();
}

//  Function setdisplayEPIC: set allow/deny EPIC energy display
//
void glWindow::setdisplayEPIC(bool a){
    displayEPIC = a;
    if (!displayEPIC) EPICpressed = false;
    update();
}

//  Function setdistances: set allow/deny distances measures
//
void glWindow::setdistances(bool a){
    distances = a;
    if (distances){
        angles = false;
        dihedrals = false;
        update();
    }
}

//  Function setdistancescolor: sets color for distances font
//
void glWindow::setdistancescolor(QColor a){
    distancescolor = a;
    if (showdistances) update();
}

//  Function setdistancesfont: set distances font
void glWindow::setdistancesfont(QFont a){
    distancesfont = a;
    if (showdistances) update();
}

//  Function setdistprecision: set number of decimals for distances
//
void glWindow::setdistprecision(int a){
    distprecision = a;
    if(showdistances) update();
}

//  Function setdistprecision: set number of decimals for distances
//
void glWindow::setdistvshift(int a){
    distvshift = a;
    if(showdistances) update();
}


//  Function setdihedrals: set allow/deny transparent background for distances lables
//
void glWindow::setdisttranspbkg(bool a){
    disttranspbkg = a;
    if(showdistances) update();
}

//  Function setdrawarcs: set allow/deny show arcs
//
void glWindow::setdrawarcs(bool a){
    drawarcs = a;
    if (showangles) update();
}

//  Function setdrawlines: set allow/deny show distances lines between nonbonded centers
//
void glWindow::setdrawlines(bool a){
    drawlines = a;
    if(showdistances) update();
}

//  Function setenergycolor: sets color for EPIC energy  display
//
void glWindow::setenergycolor(QColor a){
    energycolor = a;
    update();
}

//  Function setenergyfont: sets font for EPIC energy  display
//
void glWindow::setenergyfont(QFont a){
    energyfont = a;
    update();
}

//  Function setepicenergy: sets EPIC energy for display
//
void glWindow::setepicenergy(qreal a){
    epicenergy = a;
}

//  Function setenergyprecision: sets precision of EPIC energy for display
//
void glWindow::setenergyprecision(int a){
    energyprecision = a;
    update();
}

//  Function setfontaxeslabels: font for axes labels
void glWindow::setfontaxeslabels(QFont a){
    fontaxeslabels = a;
    if (showaxes) update();
}

//  Function setfontaxeslabels: font for axes labels
void glWindow::setfontlabaxeslabels(QFont a){
    fontlabaxeslabels = a;
    if (showaxes) update();
}

//  Function setfov: sets field of view (angle in degrees)
void glWindow::setfov(double a){
    fov = a;
}

//  Function setframeknt: sets frames counter
//
void glWindow::setframeknt(int a){
    frameknt = a;
}

//  Function sethartree: sets units for energy (true: hartree, false: kcal/mol)
//
void glWindow::sethartree(bool a){
    hartree = a;
    update();
}

//  Function setnumframes: sets number of frames to be recorded
//
void glWindow::setnumframes(int a){
    numframes = a;
}

//  Function setimagequality: sets image quality (compression)
//
void glWindow::setimagequality(int a){
    imagequality = a;
}

//  Function setimagefilename: sets image file name
//
void glWindow::setimagefilename(QString a){
    imagefilename = a;
}

//  Function setLightColor: sets color of light source
//
void glWindow::setLightColor(QColor a){
    Light_Color = QVector3D(a.redF(),a.greenF(),a.blueF());
}

//  Function setLightPower: sets power of light source
//
void glWindow::setLightPower(float a){
    Light_Power = a;
}

//  Function setlightsposition: sets light position
//
void glWindow::setlightsposition(QVector3D v){
    lightPosition = v;
    update();
}

//  Function setlinearattenuation: sets light attenuation type (true: linear, false: quadratic)
//
void glWindow::setlinearattenuation(bool a){
    linearattenuation = a;
}

//  Function setlineswidth: sets width for distances lines
//
void glWindow::setlineswidth(int a){
    lineswidth = a;
    if (showdistances && drawlines) update();
}

//  Function setlinestype: sets type for distances lines (0: NoLine, 1: SolidLine, 2: DashLine, 3: DotLine,
//          4: DashDotLine, 5: DashDotDotLine
//
void glWindow::setlinestype(int a){
    linestype = a;
    if (showdistances && drawlines) update();
}

//  Function setmeasures: set allow/deny geometry measures
//
void glWindow::setmeasures(bool a){
    measures = a;
    update();
}

//  Function setrecordcommand: command for recording
//
void glWindow::setrecordcommand(QString a){
    recordcommand = a;
}

//  Function setrecord: set allow/deny recording
//
void glWindow::setrecord(bool a){
    record = a;
    if (record)
        frames.clear();
}

//  Function setrecordfilename: sets file name for recording frames and film
//
void glWindow::setrecordfilename(QString a){
    recordfilename = a;
}

//  Function setremoveframes: set allow/deny removing frame files after movie making
//
void glWindow::setremoveframes(bool a){
    removeframes = a;
}

//  Function setshowangles: set allow/deny display angles
//
void glWindow::setshowangles(bool a){
    showangles = a;
    update();
}

//  Function setshowdihedrals: set allow/deny display dihedral angles
//
void glWindow::setshowdihedrals(bool a){
    showdihedrals = a;
    update();
}

//  Function setshowdistances: set allow/deny distances
//
void glWindow::setshowdistances(bool a){
    showdistances = a;
    update();
}

//  Function settransparentbg: sets whether background is transparent or not in capture
//
void glWindow::settransparentbg(bool a){
    transparentbg = a;
}

//  Function setSpecularColor: sets material specular color
//
void glWindow::setSpecularColor(QColor a){
    Specular_Color = QVector3D(a.redF(),a.greenF(),a.blueF());
}

//  Function setSpecularIndex: sets material specular index
//
void glWindow::setSpecularIndex(float a){
    Specular_Index = a;
}

//  Function setSpecularIndex: sets material specular index
//
void glWindow::setstepwheel(float a){
    stepwheel = a;
    update();
}

//  Function setWindowName: sets window name
//
void glWindow::setWindowName(QString name){
    windowname = name;
}

//  Function setworld_rotation: translation of the full system
void glWindow::setworld_rotation(QQuaternion a){
    world_rotation = a;
}

void glWindow::setworld_rotationAxis(QVector3D a){
    if (a.length() > 1.e-7)
        world_rotationAxis = a;
    else
        world_rotationAxis = QVector3D(1,0,0);
}

//  Function setworld_translation: translation of the full system
void glWindow::setworld_translation(QVector3D( a)){
    world_translation = a;
    update();
}

//  Function setzfar: set X axis color
void glWindow::setXaxis_color(QColor a){
    Xaxis_color = a;
    make_axes();
    if (showaxes) update();
}

//  Function setzfar: set Y axis color
void glWindow::setYaxis_color(QColor a){
    Yaxis_color = a;
    make_axes();
    if (showaxes) update();
}

//  Function setzfar: set Z axis color
void glWindow::setZaxis_color(QColor a){
    Zaxis_color = a;
    make_axes();
    if (showaxes) update();
}

//  Function setzfar: Farthest distance for objects to be displayed
void glWindow::setzfar(double a){
    zFar = a;
}

//  Function setznear: Closest distance for objects to be displayed
void glWindow::setznear(double a){
    zNear = a;
}

//**************************************************************************
//****************  GL STATE SAVING AND RESTORING    ***********************
//**************************************************************************

void glWindow::restoreGLState()
{
    mvp_list->clear();
    mvp_list = mvp_list_backup;
}

void glWindow::saveGLState()
{
    mvp_list_backup->clear();
    mvp_list_backup = mvp_list;
}


/*******************************************************************************
 * Private Helpers
 ******************************************************************************/

void glWindow::printContextInformation()
{
  QString glType;
  QString glVersion;
  QString glProfile;

  // Get Version Information
  glType = (context()->isOpenGLES()) ? "OpenGL ES" : "OpenGL";
  glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

  // Get Profile Information
#define CASE(c) case QSurfaceFormat::c: glProfile = #c; break
  switch (format().profile())
  {
    CASE(NoProfile);
    CASE(CoreProfile);
    CASE(CompatibilityProfile);
  }
#undef CASE

  // qPrintable() will print our QString w/o quotes around it.
  qDebug() << qPrintable(glType) << qPrintable(glVersion) << "(" << qPrintable(glProfile) << ")";
}
