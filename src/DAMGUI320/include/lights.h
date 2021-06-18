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
//	Header file of lights class
//  Description: lights defines a widget that manages all the elements required for
//  lights display
//
//	File:   lights.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef LIGHTS_H
#define LIGHTS_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>

#include "ColorButton.h"

#include "widgetsubclasses.h"

class lights : public QWidget
{
    Q_OBJECT
public:
    explicit lights(QWidget *parent = 0);
    ~lights();

    bool getlinearattenuation();

    QGroupBox * getFRMlights();

    void setambientcolor(QColor);
    void setbkgcolor(QColor);
    void setlightpower(float);
    void setlightcolor(QColor);
    void setlightsposition(QVector3D);
    void setlinearattenuation(bool);
    void setspecularcolor(QColor);
    void setspecularindex(float);
public slots:

signals:
    void ambientcolorchanged(QColor);
    void attenuationchanged(bool);
    void bkgcolorchanged(QColor);
    void lightcolorchanged(QColor);
    void lightschanged(QVector3D);
    void lightpowerchanged(float);
    void updateandmoveToTop();
    void specularcolorchanged(QColor);
    void specularindexchanged(float);

private slots:
    void BTNambientcolor_clicked();
    void BTNbkgcolor_clicked();
    void BTNlightcolor_clicked();
    void BTNspecularcolor_clicked();
    void create_lights_layouts();
    void lights_changed();
    void lightpower_changed();
    void RBTattenuation_changed();
    void specularindex_changed();

private:
    bool linearattenuation;

    ColorButton *BTNambientcolor;
    ColorButton *BTNbkgcolor;             // Opens dialog for background color
    ColorButton *BTNlightcolor;
    ColorButton *BTNspecularcolor;

    DoubleSpinBox *SPBlights_x;
    DoubleSpinBox *SPBlights_y;
    DoubleSpinBox *SPBlights_z;
    DoubleSpinBox *SPBlightpower;
    DoubleSpinBox *SPBspecularindex;

    float lightpower;
    float specularindex;

    QColor ambientcolor;
    QColor bkgcolor;
    QColor lightcolor;
    QColor specularcolor;

    QGroupBox *FRMattenuation;
    QGroupBox *FRMlights;

    QList<QMetaObject::Connection> connections;

    QRadioButton *RBTlinear;                // Light attenuation linear with distance
    QRadioButton *RBTsquared;               // Light attenuation squared with distance

    QVector3D lightPosition;                // Light position in world space

};

#endif // LIGHTS_H
