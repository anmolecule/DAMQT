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
//  File:   fieldlines.h
//  Description: class fieldlines manages MED gradient and electrostatic foces lines for 3D plotting
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: April 2018
//
#ifndef FIELDLINES_H
#define FIELDLINES_H

#include <QApplication>
#include <QColorDialog>
#include <QDialog>
#include <QFile>
#include <QFileDialog>
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

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class editFieldDialog : public QDialog
{
    Q_OBJECT
public:
    explicit editFieldDialog(QWidget *parent = 0);
    ~editFieldDialog();
signals:
    void closed();
protected:
    void closeEvent(QCloseEvent *event);
    virtual void reject();
};

class fieldlines : public QWidget
{
    Q_OBJECT
public:
    explicit fieldlines(QString path, QWidget *parent = 0);
    
    ~fieldlines();

    QString getname();
    QString getpath();

    QVector <GLuint> allarrowsindices;             // Indices of vertices in arrows
    QVector <GLuint> allarrowsindicesoffset;       // Offsets of indices in arrows
    QVector <VertexNormalData> allarrowsvertices;  // Vertices of polygon in arrows (position, normal, color)
    QVector <GLuint> allindices;             // Indices of vertices in lines
    QVector <GLuint> allindicesoffset;       // Offsets of indices in lines
    QVector <VertexNormalData> allvertices;  // Vertices of polygon line (position, normal, color)

    bool getshowarrows();
    bool isvisible();

    int  getarrowsseparation();
    int  getarrowssize();
    int  getarrowswidth();
    int  getlineswidth();

    QColor getlinescolor();

    void createarrows();
    void setarrowsseparation(int);
    void setarrowssize(int);
    void setarrowswidth(int);
    void setlinescolor(QColor);
    void setlineswidth(int);
    void setname(QString);
    void setpath(QString);
    void set_ProjectFolder(QString);
    void setshowarrows(bool);
    void setvisible(bool);

public slots:
    bool readfieldlines(QString);

    QString getfieldfilename();

signals:
    void updatedisplay();
    void updateRightMenu();
    
private:
    bool showarrows;                    // If true the arrows are displayed
    bool visible;                       // If true the lines are displayed


    int arrowsseparation;
    int arrowssize;
    int arrowswidth;

    int lineswidth;


    void lines2vertexes();
    void makeArrow(VertexNormalData,VertexNormalData);
    void makeCircumference();

    editFieldDialog *QDLfieldlines;

    QColor linescolor;
    
    QList<QMetaObject::Connection> connections;     // Stores connections to release in destructor

    QSpinBox *SPBarrowssep;              // Arrows separation
    QSpinBox *SPBlineswidth;             // Lines width
    QSpinBox *SPBarrowssize;             // Arrows size
    QSpinBox *SPBarrowswidth;            // Arrows size

    QString fieldfilename;
    QString name;
    QString path;                       // Path to molecule home directory (that which contains the file with lines)
    QString ProjectFolder;

    QVector <QVector3D> circumference;
    
private slots:
};

#endif // FIELDLINES_H
