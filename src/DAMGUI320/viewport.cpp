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
//	Implementation of class viewport
//
//	File:   viewport.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QSignalMapper>

#include "viewport.h"

viewport::viewport(QWidget *parent) : QWidget()
{
    connections.clear();

    FRMviewport = new QGroupBox(tr("View Port"));
    FRMviewport->setVisible(false);

    SPBfov=new DoubleSpinBox();
    SPBfov->setRange(10.,170.);
    SPBfov->setDecimals(0);
    SPBfov->setSingleStep(5);
    SPBfov->setValue(PERSPECTIVE_ANGLE);
    SPBfov->setToolTip(tr("Field of view (angle in degrees)"));
    connections << connect(SPBfov,SIGNAL(valueChanged(double)),this,SLOT(viewport_changed()));

    SPBzfar=new DoubleSpinBox();
    SPBzfar->setRange(1.e-3,1.e3);
    SPBzfar->setDecimals(2);
    SPBzfar->setSingleStep(1.);
    SPBzfar->setValue(ZFAR);
    SPBzfar->setMinimumWidth(70);
    SPBzfar->setToolTip(tr("Farthest distance for objects to be displayed"));
    connections << connect(SPBzfar,SIGNAL(valueChanged(double)),this,SLOT(SPBzfar_changed()));

    SPBznear=new DoubleSpinBox();
    SPBznear->setRange(1.e-3,1.e3);
    SPBznear->setDecimals(3);
    SPBznear->setSingleStep(1.);
    SPBznear->setValue(ZNEAR);
    SPBznear->setMinimumWidth(70);
    SPBznear->setToolTip(tr("Nearest distance for objects to be displayed"));
    connections << connect(SPBznear,SIGNAL(valueChanged(double)),this,SLOT(viewport_changed()));

    BTNviewportapply = new QPushButton();
    BTNviewportapply->setText(tr("Apply"));
    connections << connect(BTNviewportapply,SIGNAL(clicked(bool)),this,SLOT(viewport_changed()));

    create_viewport_layouts();   // Viewport layouts
}


viewport::~viewport()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

double viewport::getfov(){
    return SPBfov->value();
}

double viewport::getzfar(){
    return SPBzfar->value();
}

double viewport::getznear(){
    return SPBznear->value();
}

QGroupBox * viewport::getFRMviewport(){
    return FRMviewport;
}

//          Viewport layouts

void viewport::create_viewport_layouts(){

    QGroupBox *FRMfov = new QGroupBox(tr("Field of view"));

    QGroupBox *FRMangle = new QGroupBox(tr("Angle"));
    FRMangle->setAlignment (Qt::AlignHCenter);

    QGroupBox *FRMclippingplanes = new QGroupBox(tr("Clipping planes"));
    FRMclippingplanes->setAlignment (Qt::AlignHCenter);

    QLabel *LBLznear = new QLabel(tr("Near:"));
    QLabel *LBLzfar = new QLabel(tr("Far:"));

    QHBoxLayout *layout1 = new QHBoxLayout(FRMangle);
    layout1->addWidget(SPBfov,Qt::AlignCenter);

    QHBoxLayout *layout2 = new QHBoxLayout(FRMclippingplanes);
    layout2->addWidget(LBLznear);
    layout2->addWidget(SPBznear);
    layout2->addWidget(LBLzfar);
    layout2->addWidget(SPBzfar);

    QHBoxLayout *layout3 = new QHBoxLayout(FRMfov);
    layout3->addWidget(FRMangle);
    layout3->addWidget(FRMclippingplanes);

    QVBoxLayout *layout4 = new QVBoxLayout(FRMviewport);
    layout4->addStretch();
    layout4->addWidget(FRMfov);
    layout4->addWidget(BTNviewportapply);
    layout4->addStretch();
}

void viewport::setfov(double a){
    SPBfov->setValue(a);
}

void viewport::setzfar(double a){
    SPBzfar->setValue(a);
}

void viewport::setznear(double a){
    SPBznear->setValue(a);
}

void viewport::viewport_changed(){
    emit viewportchanged();
}

void viewport::SPBzfar_changed(){
    if (SPBznear){
        if (SPBzfar->value() <= SPBznear->value()){
            SPBznear->setValue(SPBzfar->value() - 0.01);
        }
         SPBznear->setMaximum(SPBzfar->value());
    }
    viewport_changed();
}
