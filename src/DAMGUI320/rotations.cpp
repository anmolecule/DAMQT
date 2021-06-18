//  Copyright 2008-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//	Implementation of class rotations
//
//	File:   rotations.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QFontDialog>
#include <QLabel>
#include <QPoint>
#include <QVBoxLayout>
#include <QCloseEvent>
#include <QColorDialog>
#include <QSignalMapper>

#include "rotations.h"

rotations::rotations(QWidget *parent) : QWidget()
{
    connections.clear();

    recorddialog = new recorder();

    startanimation = false;

    FRMrotation = new QGroupBox(tr("Rotation"));
    FRMrotation->setVisible(false);

    FRMangle = new QGroupBox(tr("Angle"));
    FRManimate = new QGroupBox(tr("Animate rotation"));
    FRMrotationaxis = new QGroupBox(tr("Axis"));

    SPBrot_x=new DoubleSpinBox();
    SPBrot_x->setDecimals(3);
    SPBrot_x->setSingleStep(0.01);
    SPBrot_x->setRange(-1.,1.);
    SPBrot_x->setValue(0);
    SPBrot_x->setEnabled(true);
    connections << connect(SPBrot_x,SIGNAL(valueChanged(double)),this,SIGNAL(rotation_changed()));

    SPBrot_y=new DoubleSpinBox();
    SPBrot_y->setDecimals(3);
    SPBrot_y->setSingleStep(0.01);
    SPBrot_y->setRange(-1.,1.);
    SPBrot_y->setValue(0);
    SPBrot_y->setEnabled(true);
    connections << connect(SPBrot_y,SIGNAL(valueChanged(double)),this,SIGNAL(rotation_changed()));

    SPBrot_z=new DoubleSpinBox();
    SPBrot_z->setDecimals(3);
    SPBrot_z->setSingleStep(0.01);
    SPBrot_z->setRange(-1.,1.);
    SPBrot_z->setValue(0);
    SPBrot_z->setEnabled(true);
    connections << connect(SPBrot_z,SIGNAL(valueChanged(double)),this,SIGNAL(rotation_changed()));

    SPBrot_angle=new DoubleSpinBox();
    SPBrot_angle->setDecimals(1);
    SPBrot_angle->setSingleStep(5);
    SPBrot_angle->setRange(0,360);
    SPBrot_angle->setValue(0);
//    SPBrot_angle->setValue(360. * qAcos(rotation.scalar()) / pi);
    SPBrot_angle->setEnabled(true);
    connections << connect(SPBrot_angle,SIGNAL(valueChanged(double)),this,SIGNAL(rotation_changed()));

    BTNapplyrot = new QPushButton();
    BTNapplyrot->setText(tr("Apply"));
    connections << connect(BTNapplyrot,SIGNAL(clicked()),this,SIGNAL(rotation_changed()));

    BTNresetrot = new QPushButton();
    BTNresetrot->setText(tr("Reset"));
    connections << connect(BTNresetrot,SIGNAL(clicked()),this,SLOT(reset_rotation()));

    CHKrotatex = new QCheckBox(FRManimate);
    CHKrotatex->setText(tr("X axis"));
    CHKrotatex->setChecked(false);
    connections << connect(CHKrotatex, SIGNAL(stateChanged(int)), this, SLOT(CHKrotate_changed()));

    CHKrotatey = new QCheckBox(FRManimate);
    CHKrotatey->setText(tr("Y axis"));
    CHKrotatey->setChecked(false);
    connections << connect(CHKrotatey, SIGNAL(stateChanged(int)), this, SLOT(CHKrotate_changed()));

    CHKrotatez = new QCheckBox(FRManimate);
    CHKrotatez->setText(tr("Z axis"));
    CHKrotatez->setChecked(false);
    connections << connect(CHKrotatez, SIGNAL(stateChanged(int)), this, SLOT(CHKrotate_changed()));

    BTNanimation = new QPushButton(QIcon(":/images/empezar.png"),tr("Start"));
    if (startanimation)
        BTNanimation->setText(tr("Stop"));
    connections << connect(BTNanimation,SIGNAL(clicked()), this, SLOT(BTNanimation_clicked()));

    SLDspeed = new QSlider(Qt::Horizontal);
    SLDspeed->setRange(0,INTERVAL_SCALE);
    SLDspeed->setSingleStep(1);
    SLDspeed->setPageStep(10);
    SLDspeed->setTickPosition(QSlider::TicksBelow);
    SLDspeed->setValue(INTERVAL_INI);
    connections << connect(SLDspeed,SIGNAL(valueChanged(int)), this, SIGNAL(SLDspeed_changed()));

    FRMrecord = recorddialog->getFRMrecord();

    create_rotation_layouts();   // Rotation layouts
}

rotations::~rotations()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

bool rotations::getrotatex(){
    return CHKrotatex->isChecked();
}

bool rotations::getrotatey(){
    return CHKrotatey->isChecked();
}

bool rotations::getrotatez(){
    return CHKrotatez->isChecked();
}

bool rotations::getstartanimation(){
    return startanimation;
}

int rotations::getspeed(){
    return SLDspeed->value();
}

QGroupBox * rotations::getFRMrotation(){
    return FRMrotation;
}

QQuaternion rotations::getrotation(){
    return QQuaternion::fromAxisAndAngle(SPBrot_x->value(),SPBrot_y->value(),SPBrot_z->value(),SPBrot_angle->value());
}

//  Function BTNanimation_clicked: toggles animation
//
void rotations::BTNanimation_clicked(){
    if (!(CHKrotatex->isChecked() || CHKrotatey->isChecked() || CHKrotatez->isChecked())){
        startanimation = false;
        BTNanimation->setText(tr("Start"));
        emit animationclicked(startanimation);
        return;
    }
    startanimation = !startanimation;
    if (startanimation){
        BTNanimation->setText(tr("Stop"));
//        animate(true);
    }
    else{
        BTNanimation->setText(tr("Start"));
//        animate(false);
    }
    emit animationclicked(startanimation);
}
//  End of function BTNanimation_clicked

void rotations::CHKrotate_changed(){
    bool animate = (CHKrotatex->isChecked() || CHKrotatey->isChecked() || CHKrotatez->isChecked());
    emit chkrotatechanged(animate);
}

void rotations::create_rotation_layouts(){

    QLabel *LBLfast = new QLabel(tr("Fast"));
    QLabel *LBLrot_x = new QLabel("x:");
    QLabel *LBLrot_y = new QLabel("y:");
    QLabel *LBLrot_z = new QLabel("z:");
    QLabel *LBLslow = new QLabel(tr("Slow"));
    QLabel *LBLrot_w = new QLabel("w:");

    QGridLayout *layout4 = new QGridLayout(FRMrotationaxis);
    layout4->addWidget(LBLrot_x,0,0);
    layout4->addWidget(SPBrot_x,0,1);
    layout4->addWidget(LBLrot_y,0,2);
    layout4->addWidget(SPBrot_y,0,3);
    layout4->addWidget(LBLrot_z,0,4);
    layout4->addWidget(SPBrot_z,0,5);

    QHBoxLayout *layout5 = new QHBoxLayout(FRMangle);
    layout5->addWidget(LBLrot_w);
    layout5->addWidget(SPBrot_angle);

    QHBoxLayout *layout6 = new QHBoxLayout();
    layout6->addWidget(FRMrotationaxis);
    layout6->addWidget(FRMangle);

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addStretch();
    layout7->addWidget(BTNresetrot);
    layout7->addWidget(BTNapplyrot);
    layout7->addStretch();

    QGridLayout *layout7b = new QGridLayout();
    layout7b->addWidget(CHKrotatex,0,0,Qt::AlignCenter);
    layout7b->addWidget(CHKrotatey,0,1,Qt::AlignCenter);
    layout7b->addWidget(CHKrotatez,0,2,Qt::AlignCenter);
    layout7b->addWidget(BTNanimation,1,1,Qt::AlignCenter);

    QHBoxLayout *layout7c = new QHBoxLayout();
    layout7c->addWidget(LBLslow);
    layout7c->addWidget(SLDspeed);
    layout7c->addWidget(LBLfast);

    QVBoxLayout *layout7d = new QVBoxLayout(FRManimate);
    layout7d->addLayout(layout7b);
    layout7d->addLayout(layout7c);

    QVBoxLayout *layout8 = new QVBoxLayout(FRMrotation);
    layout8->addLayout(layout6);
    layout8->addLayout(layout7);
    layout8->addWidget(FRManimate);
    layout8->addWidget(FRMrecord);
}

void rotations::reset_rotation(){
    SPBrot_x->setValue(0);
    SPBrot_y->setValue(0);
    SPBrot_z->setValue(0);
    SPBrot_angle->setValue(0);
    emit rotation_changed();
}

void rotations::setrotationButtons(QQuaternion rotation){
    SPBrot_x->setValue(rotation.x());
    SPBrot_y->setValue(rotation.y());
    SPBrot_z->setValue(rotation.z());
    SPBrot_angle->setValue(360. * qAcos(rotation.scalar()) / M_PI);
}
