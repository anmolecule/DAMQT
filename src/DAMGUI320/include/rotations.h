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
//	Header file of rotations class
//  Description: rotations defines a widget that manages all the elements required for
//  rotations of lab frame
//
//	File:   rotations.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QtCore/qmath.h>

#include "recorder.h"
#include "widgetsubclasses.h"

#define INTERVAL_INI 100
#define INTERVAL_SCALE 200

class rotations : public QWidget
{
    Q_OBJECT
public:
    explicit rotations(QWidget *parent = 0);
    ~rotations();

    recorder *recorddialog;

    bool getrotatex();
    bool getrotatey();
    bool getrotatez();
    bool getstartanimation();

    int getspeed();

    QGroupBox * getFRMrotation();

    QQuaternion getrotation();

    void setrotationButtons(QQuaternion);

signals:
    void animationclicked(bool);
    void chkrotatechanged(bool);
    void rotation_changed();
    void SLDspeed_changed();

private slots:
    void BTNanimation_clicked();
    void CHKrotate_changed();
    void create_rotation_layouts();

public slots:
    void reset_rotation();

private:
    bool rotatex;
    bool rotatey;
    bool rotatez;
    bool startanimation;                  // If true animates rotation

    DoubleSpinBox *SPBrot_angle;         // Rotation angle
    DoubleSpinBox *SPBrot_x;             // x component of rotation axis
    DoubleSpinBox *SPBrot_y;             // y component of rotation axis
    DoubleSpinBox *SPBrot_z;             // z component of rotation axis

    QCheckBox *CHKrotatex;
    QCheckBox *CHKrotatey;
    QCheckBox *CHKrotatez;

    QGroupBox *FRMangle;
    QGroupBox *FRManimate;
    QGroupBox *FRMrecord;
    QGroupBox *FRMrotation;
    QGroupBox *FRMrotationaxis;

    QList<QMetaObject::Connection> connections;

    QPushButton *BTNanimation;
    QPushButton *BTNapplyrot;
    QPushButton *BTNrecord;
    QPushButton *BTNresetrot;

    QSlider *SLDspeed;                    // Speed of rotation animation

};

#endif // ROTATIONS_H
