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
//  File:   grid.h
//  Description: class grid manages grids for 3D surface generation
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: April 2018
//
#ifndef GRID_H
#define GRID_H

#include <QObject>
#include <QProgressBar>
#include <QVector>
#include <QVector3D>
#include <QWidget>

#include "isosurface.h"
#include "CIsoSurface.h"
#include "VVBuffer.h"

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class grid : public QWidget
{
    Q_OBJECT
public:
    explicit grid(QWidget *parent = 0);

    ~grid();

    bool loadderiv(FILE *, VVBuffer *, int *, float *, int, QProgressBar *);
    bool loadderivnew(QFile *f, VVBuffer*v, int *iref, float *vref, int kntbar, QProgressBar *bar);
    bool loadnormals(FILE *f1, FILE *f2, FILE *f3, VVBuffer *v1, VVBuffer *v2, VVBuffer *v3, int *, float *, int, QProgressBar *);
    bool readplt(QString);
    bool readpltnew(QString);

    float getmaxcontourvalue();
    float getmincontourvalue();

    QList<isosurface*> *surfaces;

    QPoint getinitialposition();

    QString getfullname();
    QString getname();

    void setinitialposition(QPoint);
    void setmaxcontourvalue(float);
    void setmincontourvalue(float);
    void setfullname(QString);
    void setname(QString);
    void set_ProjectName(QString);
    void set_ProjectFolder(QString);

signals:
    void surfaceadded();
    void deletesurface();
    void surfacedeleted();

public slots:
    void addisosurf();
    void deletesurf(int);
    void generatesurf(int);
    void toggleshowsurf(int);

private: 
    bool compatderiv;

    CIsoSurface<float> *cisosurface;

    float maxcontourvalue;                // Highest contourvalue available
    float mincontourvalue;                // Lowest contourvalue available

    int nameindex;
    int nx;
    int ny;
    int nz;

    QList<QColor> surfcolors;

    QPoint initialposition;

    QVector3D sup_corner_max;
    QVector3D sup_corner_min;

    VVBuffer *fun;
    VVBuffer *dxfun;
    VVBuffer *dyfun;
    VVBuffer *dzfun;

    QString fullname;                     // Full name for grid including path
    QString name;                         // Name for grid
    QString ProjectFolder;
    QString ProjectName;
};

#endif // GRID_H
