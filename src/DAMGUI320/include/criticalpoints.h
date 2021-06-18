//  Copyright 2008-2018, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//  File:   criticalpoints.h
//  Description: class criticalpoints manages MED and MESP critical points for 3D plotting
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//
#ifndef CRITICALPOINTS_H
#define CRITICALPOINTS_H

#include <QApplication>
#include <QColorDialog>
#include <QComboBox>
#include <QDialog>
#include <QFile>
#include <QFileDialog>
#include <QFontDialog>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QObject>
#include <QOpenGLWidget>
#include <QQuaternion>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QVBoxLayout>
#include <QVector>
#include <QVector3D>
#include <QVector4D>

#include <QDebug>

#include "ColorButton.h"
#include "surface.h"

#define MAX_CPS 4		 // Number of types of critical points
#define ANGSTROM_TO_BOHR 1.889725989
#define THRESANGLES 0.000001
#define SCALE 0.01
#define CPSCALEHEIGHT 5.

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class criticalpoints : public QWidget
{
    Q_OBJECT
public:
    explicit criticalpoints(QString path, QWidget *parent = 0);
    
    ~criticalpoints();
    
//    bool create_QDLcriticalpoints();

    QColor getcolor(int);
    QColor geteigcolor(int);
    QColor getfontcolor();

    QFont getfont();

    QString getcpsfilename();
    QString getname();
    QString getpath();

    QVector <GLuint> alleigindices[MAX_CPS];            // Indices of vertices in eigenvectors of critical points
    QVector <GLuint> alleigindicesoffset[MAX_CPS];      // Offsets of indices in eigenvectors of critical points
    QVector <VertexNormalData> alleigvertices[MAX_CPS]; // Vertices of polygon eigenvectors of critical points (position, normal, color)
    QVector <GLuint> allindices[MAX_CPS];               // Indices of vertices in critical points
    QVector <GLuint> allindicesoffset[MAX_CPS];         // Offsets of indices in critical points
    QVector <VertexNormalData> allvertices[MAX_CPS];    // Vertices of polygon critical points (position, normal, color)
    QVector <QVector4D> cpsxyzval[MAX_CPS];    // Coordinates and values of critical points of types: x(3,+3), y(3,+1), z(3,-1), m(3,-3)
    QVector <QVector4D> cpseigval[MAX_CPS];    // Components of eigenvectors of hessian at critial points and sign of eigenvalue
    

    bool getcpsactive(int,int);
    bool getdrawcpscoords();
    bool getdrawcpsindices();
    bool getdrawcpssymbols();
    bool getdrawcpsvalues();
    bool geterroreigvec();
    bool getonlycpsactive();
    bool iscpsactive(int,int);
    bool isvisiblecps(int);
    bool isvisiblecpseigen();
    bool readcpsfile(QString);    // Function for reading files with critical points
    bool readcpseigen(QString);   // Function for reading files with eigenvectors of hessian at critical points

    int  getcpsactivelength(int);
    int  getcpcoordprecision();
    int  getcpprecision();
    int  getcpvshift();
    int  geteigarrowsize();
    int  geteigarrowwidth();
    int  geteiglength();
    int  geteigthickness();
    int  getradius();
    
    void createCPeigarrows();
    void createCPs();
    void createCPseigen();
    void setcolor(int, QColor);
    void setcoordprecision(int);
    void setcpprecision(int);
    void setcpsactive(int, int, bool);
    void seteigarrowsize(int);
    void seteigarrowwidth(int);
    void seteiglength(int);
    void seteigthickness(int);
    void setcpvshift(int);
    void setdrawcpscoords(bool);
    void setdrawcpsindices(bool);
    void setdrawcpssymbols(bool);
    void setdrawcpsvalues(bool);
    void seteigcolor(int, QColor);
    void setfont(QFont);
    void setfontcolor(QColor);
    void setname(QString);
    void setonlycpsactive(bool);
    void set_ProjectFolder(QString);
    void setpath(QString);
    void setradius(int);
    void setvisiblecps(int, bool);
    void setvisiblecpseigen(bool);

public slots:
    
private:  
    bool drawcpscoords;              // If true displays critical points coordinates
    bool drawcpsindices;            // If true displays critical points indices
    bool drawcpssymbols;            // If true displays critical points symbols
    bool drawcpsvalues;              // If true displays function value at critical points symbols
    bool erroreigvec;               // If true critical points eigenvectors cannot be displayed
    bool onlycpsactive;             // If true CP labels displayed only for active critical points
    bool visiblecps[MAX_CPS];       // If true the critical points are displayed
    bool visiblecpseigen;           // If true the eigenvectors of hessian at critical points are displayed

    int ballradius;
    int cpcoordprecision;
    int cpeigarrowsize;
    int cpeigarrowwidth;
    int cpeiglength;
    int cpeigthickness;
    int cpprecision;
    int cpvshift;

    uint maxindex[MAX_CPS];

    void makeCone(int slices, int stacks, float height);    // Creates a triangle strip for a cone of width 1 and given height
    void makeCylinder(int slices, int stacks);              // Creates a triangle strip for cylinder of radius 1 and height 1
    void makeSphere(int slices,int stacks);                 // Creates a triangle strip for a sphere of radius 1

    QColor cpscolor[MAX_CPS];
    QColor cpseigcolor[3];
    QColor fontcolor;

    QFont font;

    QLabel *LBLcpcoordprecision;
    QLabel *LBLcpprecision;
    QLabel *LBLcpselect;
    
    QList<QMetaObject::Connection> connections;     // Stores connections to release in destructor

    QString cpsfilename;
    QString name;
    QString path;              // Path to molecule home directory (that which contains the file with critical points)
    QString ProjectFolder;
    
//    QToolButton *BTNcps;
    
    QVector<bool> cpsactive[MAX_CPS];   // critical point active for label visualization
    QVector <GLuint> coneindices;
    QVector <GLuint> cylinderindices;
    QVector <GLuint> sphereindices;
    QVector <QVector3D> conevertices;
    QVector <QVector3D> conenormals;
    QVector <QVector3D> cylindervertices;
    QVector <QVector3D> spherevertices;

signals:
    void updatedisplay();
    void updateRightMenu();
    
private slots:

};

#endif // CRITICALPOINTS_H
