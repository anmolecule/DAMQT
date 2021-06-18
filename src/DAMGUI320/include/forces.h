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
//  File:   forces.h
//  Description: class forces manages MED and MESP critical points for 3D plotting
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: February 2019
//
#ifndef forces_H
#define forces_H

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

#define MAX_FORCES 5		 // Number of types of force types
#define ANGSTROM_TO_BOHR 1.889725989
#define SCALE 0.01
#define FORCESCALEHEIGHT 10.
#define FORCESCALEARROWSHEIGHT 1.

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class forces : public QWidget
{
    Q_OBJECT
public:
    explicit forces(QString path, QWidget *parent = 0);
    
    ~forces();

    QColor getfontcolor();

    QFont getfont();

    QString getname();
    QString getpath();

    QVector <GLuint> allindices;                // Indices of vertices in critical points
    QVector <GLuint> allindicesoffset;          // Offsets of indices in critical points
    QVector <VertexNormalData> allvertices;     // Vertices of polygon critical points (position, normal, color)
    QVector <qreal> forcesval[MAX_FORCES];                  // Forces moduli
    QVector <QVector3D> forcesorig[MAX_FORCES];             // Forces origins
    QVector <QVector3D> forcesxyz[MAX_FORCES];              // Components forces

    bool getvisibleforces(int);

    int  getforcesprecision();
    
    QColor getcolor(int);

    void setcolor(int, QColor);
    void setfont(QFont);
    void setfontcolor(QColor);
    void setname(QString);
    void set_ProjectFolder(QString);
    void setpath(QString);
    void setvisibleforces(int, bool);

public slots:

    bool isvisible();
    bool readforcesfile(QString);    // Function for reading files with HF forces
    int  getarrowlength();
    int  getarrowwidth();
    int  getforceslength();
    int  getforcesthickness();
    QString getforcesfilename();
    void createforces();
    void setarrowlength(int);
    void setarrowwidth(int);
    void setforceslength(int);
    void setforcesthickness(int);
    void setvisible(bool);
    
private:  
    bool errorforces;                   // If true forces cannot be displayed
    bool visible;
    bool visibleforces[MAX_FORCES];     // If true the pertaining forces are displayed

    int arrowlength;
    int arrowwidth;
    int numcenters;
    int forceslength;
    int forcesthickness;


    

    void makeCone(int slices, int stacks, float height);    // Creates a triangle strip for a cone of width 1 and given height
    void makeCylinder(int slices, int stacks);              // Creates a triangle strip for cylinder of radius 1 and height 1

    QColor forcecolors[MAX_FORCES];
    
    QList<QMetaObject::Connection> connections;     // Stores connections to release in destructor

    QString forcesfilename;
    QString name;
    QString path;              // Path to molecule home directory (that which contains the file with critical points)
    QString ProjectFolder;

    QVector <GLuint> coneindices;
    QVector <GLuint> cylinderindices;
    QVector <QVector3D> conevertices;
    QVector <QVector3D> conenormals;
    QVector <QVector3D> cylindervertices;

signals:
    void updatedisplay();
    
private slots:

};

#endif // forces_H
