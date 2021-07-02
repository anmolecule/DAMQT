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
//    Headers for glWindow.cpp
//  Description: glWindow implements a QOpenGLWindow for 3D display of molecules and
//  related properties using OpenGL 3.3 or higher
//  
//    File:   glWindow.h
//
//    Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//    Last version: July 2018
//

#ifndef GLWINDOW_H
#define GLWINDOW_H

#include <QBasicTimer>
#include <QCheckBox>
#include <QFileDialog>
#include <QGLFormat>
#include <QOpenGLFramebufferObject>
#include <QLabel>
#include <QMatrix4x4>
#include <QMessageBox>
#include <QOpenGLWindow>
#include <QOpenGLFunctions>
#include <QOpenGLPaintDevice>
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QPainter>
#include <QProcess>
#include <QPushButton>
#include <QProgressDialog>
#include <QScreen>
#include <QtGlobal>
#include <QTimer>
#include <QQuaternion>
#include <QVector2D>
#include <QWheelEvent>

#include <QOpenGLFunctions_4_1_Compatibility>

#include "GlobalInfo.h"

#include "geometryprocessor.h"
#include "molecule.h"
#include "elements.h"
#include "widgetsubclasses.h"
#include "dialog.h"

class geomProcessor;

#define ANGSTROM_TO_BOHR 1.889725989
//#define BOHR_TO_ANGSTROM 0.529177248882
#define MOLSCALEHEIGHT 20.
#define MOLSCALEARROWSHEIGHT 2.
#define PERSPECTIVE_ANGLE 45.    // Angle for Perspective (in degrees)
#define ZNEAR 0.001             // Znear for Perspective
#define ZFAR 200.                // Zfar  for Perspective

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class glWindow : public QOpenGLWindow, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    glWindow(QList<molecule*> *, QWidget *parent = 0);
    ~glWindow();

    bool getlinearattenuation();
    bool selectatom(int,int);
    bool selectcps(int,int);
    bool selectmespextrema(int,int);
    bool selectPopUpWindow(int,int);

    centerData searchatom(int,int);
    centerData searchcps(int,int);
    centerData searchextrema(int,int);

    float getLightPower();

    int getimagequality();
    int getframeknt();
    int searchdistance(int,int);
    float getSpecularIndex();

    QQuaternion getrotation();

    QString app_path;
    QString getimagefilename();
    QString getWindowName();

    QColor getAmbientColor();
    QColor getBkgColor();
    QColor getLightColor();
    QColor getSpecularColor();

    QVector3D getLigthPosition();
    QVector3D getrotationAxis();
    QVector3D worldTocanvas(QMatrix4x4 mvp, QRect viewport, QVector3D point);

    void checkbuffers();
    void drawangles(QPainter *painter, QRect viewport);
    void drawaxeslabels(QPainter *painter, QRect viewport, int i);
    void drawcplabels(QPainter *painter, QRect viewport, int i);
    void drawdihedrals(QPainter *painter, QRect viewport);
    void drawdistances(QPainter *painter, QRect viewport);
    void drawextremalabels(QPainter *painter, QRect viewport, int i);
    void drawlabaxeslabels(QPainter *painter, QRect viewport);
    void drawlabels(QPainter *painter, QRect viewport);
    void drawatomsvalues(QPainter *painter, QRect viewport, int i);
    void drawmeasures(QPainter *painter, QRect viewport);
    void printContextInformation();
    void reloadbuffers();
    void resetlabaxes();
    void searchmolecule(int,int);
    void setfontaxeslabels(QFont);
    void setfov(double);
    void setimagequality(int);
    void setimagefilename(QString);
    void setrecordfilename(QString);
    void setremoveframes(bool);
    void settransparentbg(bool);
    void setWindowName(QString);
    void setzfar(double);
    void setznear(double);

signals:
    void emitcheckactivate(int);
    void endmakingmovie();
    void endrecording();
    void remove_distance(int);
    void select_angle(QVector<centerData>*, QVector<QMatrix4x4>*);
    void select_dihedral(QVector<centerData>*, QVector<QMatrix4x4>*);
    void select_distance(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_angles(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_dihedrals(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_distances(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_worldrotation(QQuaternion);
    void update_worldtranslation(QVector3D);

public slots:
    void emit_update_angles();
    void emit_update_dihedrals();
    void emit_update_distances();
    void emit_update_worldtrotation();
    void emit_update_worldtranslation();
    void make_axes();
    void makeAxesCone(int slices, int stacks, float height);    // Creates a triangle strip for a cone of width 1 and given height
    void makeAxesCylinder(int slices, int stacks);              // Creates a triangle strip for cylinder of radius 1 and height 1
    void makedihedralplanes();
    void resetangles();
    void resetdihedrals();
    void resetdistances();
    void resetmeasures();
    void setAmbientColor(QColor);
    void setangles(bool);
    void setanglescolor(QColor);
    void setanglesfont(QFont);
    void setanglesprecision(int);
    void setanglestype(int);
    void setangleswidth(int);
    void setangstrom(bool);
    void setarcradius(int);
    void setaxesarrowsize(int);
    void setaxesarrowwidth(int);
    void setaxeslength(int);
    void setaxesthickness(int);
    void setaxesvisible(bool);
    void setaxeslabelsvisible(bool);
    void setBkgColor(QColor);
    void setdihedrals(bool);
    void setdihedralscolor(QColor);
    void setdihedralplanescolor(QColor);
    void setdihedralsfont(QFont);
    void setdihedralsprecision(int);
    void setdistances(bool);
    void setdistancescolor(QColor);
    void setdistancesfont(QFont);
    void setdistprecision(int);
    void setdisttranspbkg(bool);
    void setdistvshift(int);
    void setdrawarcs(bool);
    void setdrawlines(bool);
    void setfontlabaxeslabels(QFont);
    void setframeknt(int);
    void setlightsposition(QVector3D);
    void setLightColor(QColor);
    void setLightPower(float);
    void setlinearattenuation(bool);
    void setlinestype(int);
    void setlineswidth(int);
    void setmeasures(bool);
    void setnumframes(int);
    void setrecord(bool);
    void setrecordcommand(QString);
    void setshowangles(bool);
    void setshowdihedrals(bool);
    void setshowdistances(bool);
    void setSpecularColor(QColor);
    void setSpecularIndex(float);
    void setstepwheel(float);
    void setworld_rotation(QQuaternion a);
    void setworld_rotationAxis(QVector3D);
    void setworld_translation(QVector3D);
    void setXaxis_color(QColor);
    void setYaxis_color(QColor);
    void setZaxis_color(QColor);

public:
    void paintGLbuff(QSize size);

protected:
    virtual void closeEvent(QCloseEvent *);

protected:
    void keyPressEvent(QKeyEvent *e) Q_DECL_OVERRIDE;
    void mouseDoubleClickEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
    void mouseMoveEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
    void mousePressEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
    void mouseReleaseEvent(QMouseEvent *e) Q_DECL_OVERRIDE;
    void wheelEvent(QWheelEvent *e) Q_DECL_OVERRIDE;
    void timerEvent(QTimerEvent *e) Q_DECL_OVERRIDE;

    void initializeGL() Q_DECL_OVERRIDE;
    void resizeGL(int w, int h) Q_DECL_OVERRIDE;
    void paintGL() Q_DECL_OVERRIDE;

    void initShaders();
    void initTextures();

private:
    bool angles;
    bool angstrom;
    bool axesvisible;
    bool axeslabelsvisible;
    bool dihedrals;
    bool distances;
    bool disttranspbkg;
    bool drawarcs;
    bool drawlines;
    bool linearattenuation;
    bool measures;
    bool record;
    bool removeframes;
    bool showangles;
    bool showaxes;
    bool showdihedrals;
    bool showdistances;
    bool transparentbg;



    Elements *elem;

    float Light_Power;
    float Specular_Index;
    float stepwheel;

    int anglesprecision;
    int anglestype;
    int angleswidth;
    int arcradius;
    int axesarrowsize;
    int axesarrowwidth;
    int axeslength;
    int axesthickness;
    int dihedralsprecision;
    int distprecision;
    int distvshift;
    int frameknt;
    int imagequality;
    int linestype;
    int lineswidth;
    int numframes;
    int numIndices;

    void close_QDLaddmoleculeTowindow();
    void restoreGLState();
    void saveGLState();

    geomProcessor *geometries;

    QBasicTimer timer;

    QColor anglescolor;
    QColor dihedralscolor;
    QColor dihedralplanescolor;
    QColor distancescolor;
    QColor Xaxis_color;
    QColor Yaxis_color;
    QColor Zaxis_color;

    QFont anglesfont;
    QFont dihedralsfont;
    QFont distancesfont;
    QFont fontaxeslabels;
    QFont fontlabaxeslabels;

    QList<molecule*> *molecules;
    QList<QPixmap> frames;

    QMatrix4x4 projection;

    QOpenGLFramebufferObject *fbo;

    QOpenGLShaderProgram program;           // Shaders program
    QOpenGLTexture *texture;

    QPixmap pixmap;

    QQuaternion world_rotation;                   // Quaternion for rotation of the whole scene

    QString imagefilename;
    QString recordcommand;
    QString recordfilename;                 // Filename for recording frames and film
    QString windowname;                     // Window name

    QTimer *recordtimer;

    QVector<centerData> *anglecenters;              // Data of centers for measuring angles
    QVector<centerData> *dihedralcenters;           // Data of centers for measuring dihedral angles
    QVector<centerData> *distancecenters;           // Data of centers for measuring distances
    QVector <GLuint> *dihedralindices;               // Indices of vertices in dihedral planes
    QVector <VertexNormalData> *dihedralvertices;    // Vertices of triangles in structure (position, normal, color)
    QVector<GLuint> allindices;                     // Indices of vertices in structure
    QVector<VertexNormalData> allvertices;          // Vertices of triangles in structure
    QVector<GLuint> allindicesoffset;               // Offsets of indices in objects (cylinders and spheres) in structure
    QVector<QMatrix4x4> *mvp_list;
    QVector<QMatrix4x4> *mvp_list_backup;
    QVector<QMatrix4x4> *v_list;
    QVector <GLuint> *allaxesindices;           // Indices of axes vertices
    QVector <GLuint> *allaxesindicesoffset;     // Offsets of axes indices in objects
    QVector <VertexNormalData> *allaxesvertices;  // Vertices of axes triangles(position, normal, color)
    QVector <GLuint> coneindices;
    QVector <GLuint> cylinderindices;
    QVector <QVector3D> conevertices;
    QVector <QVector3D> conenormals;
    QVector <QVector3D> cylindervertices;
    QVector <QVector3D> *positionaxeslabels;

    QVector2D mousePressPosition;           // Cursor position when mouse buttons pressed
    QVector2D mousePreviousPosRot;          // Last cursor position when rotating
    QVector2D mousePreviousPosTrs;          // Last cursor position when translating

    QVector3D lightPosition;                // Light position in world space
    QVector3D world_rotationAxis;           // Rotation axis
    QVector3D world_translation;            // Translation vector

    qreal angularSpeed;
    qreal fov;                              // Field of view (angle in degrees)
    qreal disthres;                         // Distance threshold for bonding
    qreal zFar;                             // Farthest distance for objects to be displayed
    qreal zNear;                            // Nearest distance for objects to be displayed

    QVector3D Ambient_Color;                // Ambient light color for the shaders
    QVector3D Bkg_Color;                    // Background color
    QVector3D Light_Color;                  // Source light color for the shaders
    QVector3D Specular_Color;               // Reflection color for the shaders

};
    
#endif // GLWINDOW_H
    
