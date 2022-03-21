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
//	Implementation of class lights
//
//	File:   lights.cpp
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

#include "lights.h"

lights::lights(QWidget *parent) : QWidget()
{
    ambientcolor = QColor(128,128,128);
    lightcolor = QColor(230,230,230);
    specularcolor = QColor(128,128,128);

    lightpower = 20.f;
    specularindex = 5.;

    linearattenuation = true;

    connections.clear();

    FRMlights = new QGroupBox(tr("Lights properties"));
    FRMlights->setVisible(false);

    SPBlights_x=new DoubleSpinBox();
    SPBlights_x->setRange(-100.,1000.);
    SPBlights_x->setDecimals(1);
    SPBlights_x->setSingleStep(0.5);
    SPBlights_x->setValue(lightPosition.x());
    connections << connect(SPBlights_x,SIGNAL(valueChanged(double)),this,SLOT(lights_changed()));

    SPBlights_y=new DoubleSpinBox();
    SPBlights_y->setRange(-100.,1000.);
    SPBlights_y->setDecimals(1);
    SPBlights_y->setSingleStep(0.5);
    SPBlights_y->setValue(lightPosition.y());
    connections << connect(SPBlights_y,SIGNAL(valueChanged(double)),this,SLOT(lights_changed()));

    SPBlights_z=new DoubleSpinBox();
    SPBlights_z->setRange(-100.,1000.);
    SPBlights_z->setDecimals(1);
    SPBlights_z->setSingleStep(0.5);
    SPBlights_z->setValue(lightPosition.z());
    connections << connect(SPBlights_z,SIGNAL(valueChanged(double)),this,SLOT(lights_changed()));

    BTNlightcolor = new ColorButton();
    BTNlightcolor->setIcon(QIcon(":/images/colores.png"));
    BTNlightcolor->setText(tr("Color"));
    BTNlightcolor->setColor(&lightcolor);
    BTNlightcolor->setEnabled(true);
    connections << connect(BTNlightcolor,SIGNAL(clicked()),this,SLOT(BTNlightcolor_clicked()));

    BTNambientcolor = new ColorButton();
    BTNambientcolor->setIcon(QIcon(":/images/colores.png"));
    BTNambientcolor->setText(tr("Color"));
    BTNambientcolor->setColor(&ambientcolor);
    BTNambientcolor->setEnabled(true);
    connections << connect(BTNambientcolor,SIGNAL(clicked()),this,SLOT(BTNambientcolor_clicked()));

    BTNspecularcolor = new ColorButton();
    BTNspecularcolor->setIcon(QIcon(":/images/colores.png"));
    BTNspecularcolor->setText(tr("Color"));
    BTNspecularcolor->setColor(&specularcolor);
    BTNspecularcolor->setEnabled(true);
    connections << connect(BTNspecularcolor,SIGNAL(clicked()),this,SLOT(BTNspecularcolor_clicked()));

    BTNbkgcolor = new ColorButton();
    BTNbkgcolor->setIcon(QIcon(":/images/colores.png"));
    BTNbkgcolor->setText(tr("Color"));
    BTNbkgcolor->setColor(&bkgcolor);
    BTNbkgcolor->setEnabled(true);
    connections << connect(BTNbkgcolor,SIGNAL(clicked()),this,SLOT(BTNbkgcolor_clicked()));

    SPBlightpower=new DoubleSpinBox();
    SPBlightpower->setRange(0.,3000.);
    SPBlightpower->setDecimals(0);
    SPBlightpower->setSingleStep(10);
    SPBlightpower->setValue(lightpower);
    connections << connect(SPBlightpower,SIGNAL(valueChanged(double)),this,SLOT(lightpower_changed()));

    SPBspecularindex=new DoubleSpinBox();
    SPBspecularindex->setRange(0.,20.);
    SPBspecularindex->setDecimals(0);
    SPBspecularindex->setSingleStep(1);
    SPBspecularindex->setValue(specularindex);
    connections << connect(SPBspecularindex,SIGNAL(valueChanged(double)),this,SLOT(specularindex_changed()));

    FRMattenuation = new QGroupBox(tr("Light attenuation"));
    RBTlinear=new QRadioButton(tr("Linear"),FRMattenuation);
    RBTlinear->setChecked(true);
    setlinearattenuation(true);
    RBTsquared=new QRadioButton(tr("Squared"),FRMattenuation);
    RBTsquared->setChecked(false);
    connections << connect(RBTlinear,SIGNAL(toggled(bool)),this,SLOT(RBTattenuation_changed()));

    create_lights_layouts();     // Lights menu layouts
}

lights::~lights()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}


bool lights::getlinearattenuation(){
    return linearattenuation;
}

QGroupBox * lights::getFRMlights(){
    return FRMlights;
}

void lights::BTNambientcolor_clicked(){
    QColor currcolor = ambientcolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        ambientcolor.setRgb(newcolor.red(), newcolor.green(), newcolor.blue());
        BTNambientcolor->setColor(&newcolor);
        emit ambientcolorchanged(newcolor);
        emit updateandmoveToTop();
    }
}

void lights::BTNbkgcolor_clicked(){
    QColor currcolor = bkgcolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        bkgcolor.setRgb(newcolor.red(), newcolor.green(), newcolor.blue());
        BTNbkgcolor->setColor(&newcolor);
        emit bkgcolorchanged(newcolor);
        emit updateandmoveToTop();
    }
}

void lights::BTNlightcolor_clicked(){
    QColor currcolor = lightcolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        lightcolor.setRgb(newcolor.red(), newcolor.green(), newcolor.blue());
        BTNlightcolor->setColor(&newcolor);
        emit lightcolorchanged(newcolor);
        emit updateandmoveToTop();
    }
}

void lights::BTNspecularcolor_clicked(){
    QColor currcolor = specularcolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        specularcolor.setRgb(newcolor.red(), newcolor.green(), newcolor.blue());
        BTNspecularcolor->setColor(&newcolor);
        emit specularcolorchanged(newcolor);
        emit updateandmoveToTop();
    }
}

//          Lights menu layouts

void lights::create_lights_layouts(){

    QLabel *LBLlights_x = new QLabel("x:");
    QLabel *LBLlights_y = new QLabel("y:");
    QLabel *LBLlights_z = new QLabel("z:");

    QHBoxLayout *layout2a = new QHBoxLayout();
    layout2a->addStretch();
    layout2a->addWidget(LBLlights_x);
    layout2a->addWidget(SPBlights_x);
    layout2a->addWidget(LBLlights_y);
    layout2a->addWidget(SPBlights_y);
    layout2a->addWidget(LBLlights_z);
    layout2a->addWidget(SPBlights_z);
    layout2a->addStretch();

    QLabel *LBLlightcolor = new QLabel(tr("Source light color"));
    QLabel *LBLambientcolor = new QLabel(tr("Ambient light color"));
    QLabel *LBLspecularcolor = new QLabel(tr("Material specular color"));
    QLabel *LBLbkgcolor = new QLabel(tr("Background color"));
    QLabel *LBLlightpower = new QLabel(tr("Light Power"));
    QLabel *LBLspecularindex = new QLabel(tr("Specular index"));

    QGridLayout *layout2b = new QGridLayout();
    layout2b->addWidget(LBLlightcolor,0,0);
    layout2b->addWidget(BTNlightcolor,0,1);
    layout2b->addWidget(LBLambientcolor,1,0);
    layout2b->addWidget(BTNambientcolor,1,1);
    layout2b->addWidget(LBLspecularcolor,2,0);
    layout2b->addWidget(BTNspecularcolor,2,1);
    layout2b->addWidget(LBLbkgcolor,3,0);
    layout2b->addWidget(BTNbkgcolor,3,1);
    layout2b->addWidget(LBLlightpower,4,0);
    layout2b->addWidget(SPBlightpower,4,1);
    layout2b->addWidget(LBLspecularindex,5,0);
    layout2b->addWidget(SPBspecularindex,5,1);

    QHBoxLayout *layout2d = new QHBoxLayout(FRMattenuation);
    layout2d->addStretch();
    layout2d->addWidget(RBTlinear);
    layout2d->addWidget(RBTsquared);
    layout2d->addStretch();

    QVBoxLayout *layout2e = new QVBoxLayout(FRMlights);
    layout2e->addLayout(layout2a);
    layout2e->addLayout(layout2b);
    layout2e->addWidget(FRMattenuation);
    layout2e->addStretch();
}

void lights::lights_changed(){
    lightPosition.setX(SPBlights_x->value());
    lightPosition.setY(SPBlights_y->value());
    lightPosition.setZ(SPBlights_z->value());
    emit lightschanged(lightPosition);
    emit updateandmoveToTop();
}

void lights::lightpower_changed(){
    lightpower = SPBlightpower->value();
    emit lightpowerchanged(lightpower);
    emit updateandmoveToTop();
}

void lights::RBTattenuation_changed(){
    if (RBTlinear->isChecked()){
        setlinearattenuation(true);
        lightpower = 20.f;
    }
    else{
        setlinearattenuation(false);
        lightpower = 300.f;
    }
    BTNambientcolor->setColor(&ambientcolor);
    BTNspecularcolor->setColor(&specularcolor);
    BTNlightcolor->setColor(&lightcolor);
    SPBlightpower->setValue(lightpower);
    SPBlights_x->setValue(lightPosition.x());
    SPBlights_y->setValue(lightPosition.y());
    SPBlights_z->setValue(lightPosition.z());
    emit attenuationchanged(linearattenuation);
    emit lightpowerchanged(lightpower);
    emit updateandmoveToTop();
}

void lights::setambientcolor(QColor a){
    ambientcolor = a;
    BTNambientcolor->setColor(&ambientcolor);
}

void lights::setbkgcolor(QColor a){
    bkgcolor = a;
    BTNbkgcolor->setColor(&bkgcolor);
}

void lights::setlightcolor(QColor a){
    lightcolor = a;
    BTNlightcolor->setColor(&lightcolor);
}

void lights::setlightpower(float a){
    lightpower = a;
    SPBlightpower->setValue(a);
}

void lights::setlightsposition(QVector3D a){
    lightPosition = a;
    SPBlights_x->setValue(lightPosition.x());
    SPBlights_y->setValue(lightPosition.y());
    SPBlights_z->setValue(lightPosition.z());
}

void lights::setlinearattenuation(bool a){
    linearattenuation = a;
    RBTlinear->setChecked(a);
}

void lights::setspecularcolor(QColor a){
    specularcolor = a;
    BTNspecularcolor->setColor(&specularcolor);
}

void lights::setspecularindex(float a){
    specularindex = a;
    SPBspecularindex->setValue(a);
}

void lights::specularindex_changed(){
    specularindex = SPBspecularindex->value();
    emit specularindexchanged(specularindex);
    emit updateandmoveToTop();
}
