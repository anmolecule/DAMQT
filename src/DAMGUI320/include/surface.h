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
//  File:   surface.h
//
//      Last version: March 2016
//
#ifndef SURFACE_H
#define SURFACE_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QFileInfo>
#include <QGroupBox>
#include <QKeyEvent>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QPoint>
#include <QProgressBar>
#include <QRadioButton>
#include <QSlider>
#include <QSpinBox>
#include <QToolButton>
#include <QVector>
#include <QVector3D>
#include <QVector4D>

#include "ColorButton.h"
#include "widgetsubclasses.h"
#include "VertexNormalData.h"

#define pi           3.14159265358979323846    // pi
#define LOGBASE      1.08                      // base for logaritmics scale
#define SCALE 0.01

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

// Class surface
//
class surface : public QWidget
{
    friend class geomProcessor;
    Q_OBJECT
public:
    explicit surface(QWidget *parent = 0);
    ~surface();

    QColor getsurfacecolor();
    QColor getfontcolor();

    QFont getfont();

    bool getextremactive(int,int);
    bool getextremhidden(int,int);
    bool getvisible();
    bool getonlyextremactive();
    bool getshowcolorrule();
    bool getshowextremacoords();
    bool getshowextremaindices();
    bool getshowextremasymbols();
    bool getshowextremavalues();
    bool getshowgridbounds();
    bool getshowlocalmax();
    bool getshowlocalmin();
    bool gettranslucence();
    bool getsolidsurf();
    bool readbasins(QString);
    bool readbasinsnew(QString);
    bool readisosurfnew(QString);
    bool readsgh(QString);
    bool readsghnew(QString);
    bool readsrf(QString);
    bool readsrfnew(QString);

    int  getcoordprecision();
    int  getradius();
    int  getvalueprecision();
    int  getvshift();

    float getfabstop();
    float getopacity();
    float gettopcolor();

    QPoint getinitialposition();

    QString getfullname();
    QString getname();

    QVector3D gettranslation();
    QVector3D getrotationAxis();

    QVector<QVector4D> localextrema[2];        // Positions and mesp values of local extrema (0: maxima, 1: minima)

    QVector <GLuint> getallindices();
    QVector <VertexNormalData> getallvertices();
    QVector <GLuint> allindicesextrema[2];               // Indices of vertices in critical points
    QVector <GLuint> allindicesoffsetextrema[2];         // Offsets of indices in critical points
    QVector <VertexNormalData> allverticesextrema[2];    // Vertices of polygon critical points (position, normal, color)
    QVector <GLuint> allindices;              // indices of vertices in surface
    QVector <GLfloat> allvalues;              // values of function on surface triangles vertices
    QVector <VertexNormalData> allvertices;   // vertices of triangles in surface (position, normal, color)
    QVector <GLuint> gridindices;             // indices of vertices in grid boundaries
    QVector <GLuint> gridindicesoffset;       // offsets of indices in grid boundaries
    QVector <VertexNormalData> gridvertices;  // vertices of triangles in grid boundaries (position, normal, color)

    void surfacecolor_changed(QColor);
    void resetsurface();
    void setextremactive(int, int, bool);
    void setextremhidden(int, int, bool);
    void setinitialposition(QPoint);
    void setfullname(QString);
    void setname(QString);
    void setopacity(float);
    void setvisible(bool);

signals:
    void surfaceColor(QColor *);
    void fontColor(QColor *);
    void maximaColor(QColor *);
    void minimaColor(QColor *);
    void opendialog();
    void toggleshowbutton();    
    void updatedisplay();
    void updatelabelcolor(QString);
    void updateRightMenu();
    void updatetext(QString);

public slots:
    int getballradius();

    QColor getextremacolor(int);

    void changeextremafont();
    void extremaselectall();
    void extremaselectnone();
    void opacity_changed(int);
    void opacity_released();
    void setshowgridbounds(bool);
    void setshowlocalmax(bool);
    void setshowlocalmin(bool);
    void settranslucence(bool);
    void selectsurfacecolor();
    void selectfontcolor();
    void selectmaximacolor();
    void selectminimacolor();
    void setballradius(int);
    void setcoordprecision(int);
    void setonlyextremactive(bool);
    void setshowextremacoords(bool);
    void setshowextremaindices(bool);
    void setshowextremasymbols(bool);
    void setshowextremavalues(bool);
    void setsolidsurf(bool);
    void settopcolor(float);
    void settrianglecolors();
    void setvalueprecision(int);
    void setvshift(int);
    void toggleshowsurf();

private:
    float getscalevalueFloat(int, int, int, float, float);
    int getscalevalueInt(float, int, int, float, float);
    void createballs();
    void generategridbounds(float *);
    void makeSphere(int slices,int stacks);                 // Creates a triangle strip for a sphere of radius 1

    bool active;                        // If true, the surface is active for transformations
    bool onlyextremactive;              // If true atom labels displayed only for active critical points
    bool showcolorrule;                 // If true, a file .sgh is loaded: shows a rule for color selection in spectrogram
    bool solidsurf;                     // If true, surface display is solid, if false, surface display is wire frame
    bool showgridbounds;                 // If true displays grid boundaries
    bool showextremacoords;             // If true displays coordinates of local extrema
    bool showextremaindices;            // If true displays symbols of local extrema
    bool showextremasymbols;            // If true displays symbols of local extrema
    bool showextremavalues;             // If true displays values of local extrema
    bool showlocalmax;                  // If true displays local maxima
    bool showlocalmin;                  // If true displays local minima
    bool translucence;                  // If true, translucence correction is applied
    bool visible;                       // If true the surface is displayed

    int ballradius;
    int coordprecision;
    int nindices;
    int nlocalmax;
    int nlocalmin;
    int nvertices;
    int valueprecision;
    int vshift;

    uint maxindex[2];

    float fmax;
    float fmin;
    float fabstop;
    float opacity;                        // Opacity: 1 (opaque) 0 (transparent)
    float topcolor;

    QColor surfacecolor;
    QColor extremacolor[2];
    QColor fontcolor;

    QDoubleValidator *myDoubleValidator;

    QFont font;

    QList<QMetaObject::Connection> connections;

    QPoint initialposition;

    QString fullname;                     // Full name for surface including path
    QString name;                         // Name for surface

    QVector<bool> extremactive[2];         // extrema points active for label visualization
    QVector<bool> extremhidden[2];         // extrema points hidden for label visualization
    QVector<float> griddimensions;      // Original grid dimensions: xmin, xmax, ymin, ymax, zmin, zmax
    QVector <GLuint> sphereindices;
    QVector<int> gridnxyz;              // Original grid number of points (nx, ny, nz)
    QVector <QVector3D> spherevertices;
};

#endif // SURFACE_H
