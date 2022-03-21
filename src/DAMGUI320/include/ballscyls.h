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
//	Header file of axes class
//  Description: axes defines a widget that manages all the elements required for
//  axes display
//
//	File:   axes.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef BALLSCYLS_H
#define BALLSCYLS_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QList>
#include <QPushButton>

#include "widgetsubclasses.h"

#define INIT_BOND_THRESHOLD 1.2

class ballscyls : public QWidget
{
    Q_OBJECT
public:
    explicit ballscyls(QWidget *parent = 0);
    ~ballscyls();

    double getballrad();
    double getcylrad();
    double getbondthres();

    QGroupBox * getFRMballcyl();

    void setballrad(double);
    void setcylrad(double);

signals:
    void ballcylchanged();
    void scaleradiichanged(bool);
private slots:
    void CHKscaleradii_changed();
    void create_balls_cyl_layouts();
private:
    DoubleSpinBox *SPBballrad;
    DoubleSpinBox *SPBcylrad;
    DoubleSpinBox *SPBbondthres;

    QCheckBox *CHKscaleradii;

    QGroupBox *FRMballcyl;

    QList<QMetaObject::Connection> connections;

    QPushButton *BTNballcylapply;
};

#endif // BALLSCYLS_H
