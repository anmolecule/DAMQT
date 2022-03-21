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
//	Implementation of class axes
//
//	File:   ballscyls.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QPoint>
#include <QVBoxLayout>
#include <QCloseEvent>
#include <QColorDialog>
#include <QSignalMapper>

#include "ballscyls.h"

ballscyls::ballscyls(QWidget *parent){
    connections.clear();

    FRMballcyl = new QGroupBox();
    FRMballcyl->setVisible(false);

    CHKscaleradii = new QCheckBox();
    CHKscaleradii->setChecked(true);
    CHKscaleradii->setText(tr("Scale balls with atoms radii"));

    connections << connect(CHKscaleradii, SIGNAL(stateChanged(int)), this, SLOT(CHKscaleradii_changed()));

    SPBballrad=new DoubleSpinBox();
    SPBballrad->setRange(0.,10.);
    SPBballrad->setDecimals(2);
    SPBballrad->setSingleStep(0.1);
    SPBballrad->setValue(0.2);
    connections << connect(SPBballrad,SIGNAL(valueChanged(double)),this,SIGNAL(ballcylchanged()));

    SPBcylrad=new DoubleSpinBox();
    SPBcylrad->setRange(0.,10.);
    SPBcylrad->setDecimals(2);
    SPBcylrad->setSingleStep(0.05);
    SPBcylrad->setValue(0.05);
    connections << connect(SPBcylrad,SIGNAL(valueChanged(double)),this,SIGNAL(ballcylchanged()));

    SPBbondthres=new DoubleSpinBox();
    SPBbondthres->setRange(0.,10.);
    SPBbondthres->setDecimals(2);
    SPBbondthres->setSingleStep(0.05);
    SPBbondthres->setValue(INIT_BOND_THRESHOLD);
    connections << connect(SPBbondthres,SIGNAL(valueChanged(double)),this,SIGNAL(ballcylchanged()));

    BTNballcylapply = new QPushButton();
    BTNballcylapply->setText(tr("Apply"));
    connections << connect(BTNballcylapply,SIGNAL(clicked(bool)),this,SIGNAL(ballcylchanged()));

    create_balls_cyl_layouts();
}

ballscyls::~ballscyls(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

double ballscyls::getballrad(){
    return SPBballrad->value();
}

double ballscyls::getbondthres(){
    return SPBbondthres->value();
}

double ballscyls::getcylrad(){
    return SPBcylrad->value();
}

QGroupBox * ballscyls::getFRMballcyl(){
    return FRMballcyl;
}

void ballscyls::CHKscaleradii_changed(){
    emit scaleradiichanged(CHKscaleradii->isChecked());
}

void ballscyls::setballrad(double rad){
    SPBballrad->setValue(rad);
}

void ballscyls::setcylrad(double rad){
    SPBcylrad->setValue(rad);
}


//          Balls and cylinders layouts

void ballscyls::create_balls_cyl_layouts(){

    QLabel *LBLballrad = new QLabel(tr("Ball"));
    QLabel *LBLcylrad = new QLabel(tr("Cylinder"));
    QLabel *LBLbondthres = new QLabel(tr("Bond threshold"));

    QHBoxLayout *layout1 = new QHBoxLayout();
    layout1->addStretch();
    layout1->addWidget(CHKscaleradii,Qt::AlignCenter);
    layout1->addStretch();

    QGridLayout *layout3a = new QGridLayout();
    layout3a->addWidget(LBLballrad,0,0,Qt::AlignCenter);
    layout3a->addWidget(LBLcylrad,0,1,Qt::AlignCenter);
    layout3a->addWidget(LBLbondthres,0,2,Qt::AlignCenter);
    layout3a->addWidget(SPBballrad,1,0);
    layout3a->addWidget(SPBcylrad,1,1);
    layout3a->addWidget(SPBbondthres,1,2);

    QVBoxLayout *layout3b = new QVBoxLayout(FRMballcyl);
    layout3b->addStretch();
    layout3b->addLayout(layout1);
    layout3b->addLayout(layout3a);
    layout3b->addWidget(BTNballcylapply);
    layout3b->addStretch();
}
