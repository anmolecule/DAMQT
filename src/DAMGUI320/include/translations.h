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
//	Header file of translations class
//  Description: translations defines a widget that manages all the elements required for
//  translations of lab frame
//
//	File:   translations.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef TRANSLATIONS_H
#define TRANSLATIONS_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QPushButton>
#include <QRadioButton>
#include <QSlider>
#include <QSpinBox>
#include <QtCore/qmath.h>

#include "widgetsubclasses.h"

#define ANGSTROM_TO_BOHR 1.889725989
#define BOHR_TO_ANGSTROM 0.529177249
#define Z_TRANS_INI -20.

class translations : public QWidget
{
    Q_OBJECT
public:
    explicit translations(QWidget *parent = 0);
    ~translations();

    QGroupBox * getFRMtranslation();
public slots:
    QVector3D gettranslation();
    void set_translation(QVector3D);
signals:
    void setworldtranslation(QVector3D);
    void stepwheel_changed(float);

private slots:
    void create_translation_layouts();
    void RBTbohr_changed();
    void SPBstepwheel_changed();
    void translation_changed();

public slots:
    void reset_translation();

private:
    DoubleSpinBox *SPBstepwheel;         // stride for translation with mouse wheel
    DoubleSpinBox *SPBtras_x;            // x component of translation vector
    DoubleSpinBox *SPBtras_y;            // y component of translation vector
    DoubleSpinBox *SPBtras_z;            // z component of translation vector

    QPushButton *BTNapplytrans;
    QPushButton *BTNresetrans;

    QGroupBox *FRMtranslation;
    QGroupBox *FRMtranslationunits;

    QList<QMetaObject::Connection> connections;

    QRadioButton *RBTangstrom;
    QRadioButton *RBTbohr;

    QVector3D translation;                  // Translation vector
};

#endif // TRANSLATIONS_H
