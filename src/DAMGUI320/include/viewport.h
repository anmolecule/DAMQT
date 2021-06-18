
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
//	Header file of viewport class
//  Description: viewport defines a widget that manages all the elements required for
//  viewport display
//
//	File:   viewport.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef VIEWPORT_H
#define VIEWPORT_H

#include <QObject>
#include <QWidget>
#include <QGroupBox>
#include <QList>
#include <QPushButton>

#include "widgetsubclasses.h"


#define PERSPECTIVE_ANGLE 45.	// Angle for Perspective (in degrees)
#define ZNEAR 0.001             // Znear for Perspective
#define ZFAR 200.				// Zfar  for Perspective

class viewport : public QWidget
{
    Q_OBJECT
public:
    explicit viewport(QWidget *parent = 0);
    ~viewport();

    double getfov();
    double getzfar();
    double getznear();

    QGroupBox * getFRMviewport();

    void setfov(double);
    void setznear(double);
    void setzfar(double);

signals:
    void viewportchanged();

private slots:
    void create_viewport_layouts();
    void viewport_changed();
    void SPBzfar_changed();

private:
    DoubleSpinBox *SPBfov;
    DoubleSpinBox *SPBzfar;
    DoubleSpinBox *SPBznear;

    QGroupBox *FRMviewport;

    QPushButton *BTNviewportapply;

    QList<QMetaObject::Connection> connections;
};


#endif // VIEWPORT_H
