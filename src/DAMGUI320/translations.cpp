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
//	Implementation of class translations
//
//	File:   translations.cpp
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

#include "translations.h"

translations::translations(QWidget *parent) : QWidget()
{
    translation = QVector3D(0,0,Z_TRANS_INI);

    connections.clear();

    FRMtranslation = new QGroupBox(tr("Translation"));
    FRMtranslation->setVisible(false);

    FRMtranslationunits = new QGroupBox(tr("Units"));

    RBTangstrom = new QRadioButton(tr("angstrom"));
    RBTbohr = new QRadioButton(tr("bohr"));
    RBTbohr->setChecked(true);
    RBTangstrom->setChecked(false);
    connections << connect(RBTbohr, SIGNAL(toggled(bool)),this,SLOT(RBTbohr_changed()));

    SPBtras_x = new DoubleSpinBox();
    SPBtras_x->setDecimals(2);
    SPBtras_x->setSingleStep(0.1);
    SPBtras_x->setRange(-1000,1000);
    SPBtras_x->setValue(translation.x());
    SPBtras_x->setEnabled(true);
    connections << connect(SPBtras_x,SIGNAL(valueChanged(double)),this,SLOT(translation_changed()));

    SPBtras_y = new DoubleSpinBox();
    SPBtras_y->setDecimals(2);
    SPBtras_y->setSingleStep(0.1);
    SPBtras_y->setRange(-1000,1000);
    SPBtras_y->setValue(translation.y());
    SPBtras_y->setEnabled(true);
    connections << connect(SPBtras_y,SIGNAL(valueChanged(double)),this,SLOT(translation_changed()));

    SPBtras_z=new DoubleSpinBox();
    SPBtras_z->setDecimals(2);
    SPBtras_z->setSingleStep(0.1);
    SPBtras_z->setRange(-1000,1000);
    SPBtras_z->setValue(translation.z());
    SPBtras_z->setEnabled(true);
    connections << connect(SPBtras_z,SIGNAL(valueChanged(double)),this,SLOT(translation_changed()));

    BTNapplytrans = new QPushButton();
    BTNapplytrans->setText(tr("Apply"));
    connections << connect(BTNapplytrans,SIGNAL(clicked()),this,SLOT(translation_changed()));
    connections << connect(BTNapplytrans,SIGNAL(clicked()),this,SLOT(SPBstepwheel_changed()));

    BTNresetrans = new QPushButton();
    BTNresetrans->setText(tr("Reset"));
    connections << connect(BTNresetrans,SIGNAL(clicked()),this,SLOT(reset_translation()));

    SPBstepwheel = new DoubleSpinBox();
    SPBstepwheel->setDecimals(2);
    SPBstepwheel->setSingleStep(0.1);
    SPBstepwheel->setRange(0,100);
    SPBstepwheel->setValue(0.1);
    SPBstepwheel->setEnabled(true);
    connections << connect(SPBstepwheel,SIGNAL(valueChanged(double)),this,SLOT(SPBstepwheel_changed()));

    create_translation_layouts();   // Translation layouts
}

translations::~translations()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

QGroupBox * translations::getFRMtranslation(){
    return FRMtranslation;
}

void translations::create_translation_layouts(){
    QLabel *LBLtras_x = new QLabel("x:");
    QLabel *LBLtras_y = new QLabel("y:");
    QLabel *LBLtras_z = new QLabel("z:");
    QLabel *LBLstepwheel = new QLabel(tr("Stride for zooming with mouse wheel: "));

    QHBoxLayout *layout8 = new QHBoxLayout(FRMtranslationunits);
    layout8->addWidget(RBTbohr);
    layout8->addWidget(RBTangstrom);

    QHBoxLayout *layout9 = new QHBoxLayout();
    layout9->addStretch();
    layout9->addWidget(LBLtras_x);
    layout9->addWidget(SPBtras_x);
    layout9->addWidget(LBLtras_y);
    layout9->addWidget(SPBtras_y);
    layout9->addWidget(LBLtras_z);
    layout9->addWidget(SPBtras_z);
    layout9->addStretch();

    QHBoxLayout *layout10 = new QHBoxLayout();
    layout10->addStretch();
    layout10->addWidget(LBLstepwheel);
    layout10->addWidget(SPBstepwheel);
    layout10->addStretch();

    QHBoxLayout *layout11 = new QHBoxLayout();
    layout11->addStretch();
    layout11->addWidget(BTNresetrans);
    layout11->addWidget(BTNapplytrans);
    layout11->addStretch();

    QVBoxLayout *layout12 = new QVBoxLayout(FRMtranslation);
    layout12->addWidget(FRMtranslationunits);
    layout12->addLayout(layout9);
    layout12->addLayout(layout10);
    layout12->addLayout(layout11);
}


//  Function RBTbohr_changed: toggles between units
//
void translations::RBTbohr_changed(){
    if (RBTbohr->isChecked()){
        SPBtras_x->setValue(translation.x());
        SPBtras_y->setValue(translation.y());
        SPBtras_z->setValue(translation.z());
    }
    else{
        SPBtras_x->setValue(translation.x()*BOHR_TO_ANGSTROM);
        SPBtras_y->setValue(translation.y()*BOHR_TO_ANGSTROM);
        SPBtras_z->setValue(translation.z()*BOHR_TO_ANGSTROM);
    }
}

QVector3D translations::gettranslation(){
    return translation;
}

void translations::reset_translation(){
    SPBtras_x->setValue(0);
    SPBtras_y->setValue(0);
    if (RBTangstrom->isChecked()){
        SPBtras_z->setValue(Z_TRANS_INI*BOHR_TO_ANGSTROM);
        translation = QVector3D(0,0,Z_TRANS_INI*BOHR_TO_ANGSTROM);
    }
    else{
        SPBtras_z->setValue(Z_TRANS_INI);
        translation = QVector3D(0,0,Z_TRANS_INI);
    }
    translation_changed();
}

void translations::set_translation(QVector3D a){
    translation = a;
    if (RBTangstrom->isChecked()){
        SPBtras_x->setValue(translation.x()*BOHR_TO_ANGSTROM);
        SPBtras_y->setValue(translation.y()*BOHR_TO_ANGSTROM);
        SPBtras_z->setValue(translation.z()*BOHR_TO_ANGSTROM);
    }
    else{
        SPBtras_x->setValue(translation.x());
        SPBtras_y->setValue(translation.y());
        SPBtras_z->setValue(translation.z());
    }
}

//  Function SPBstepwheel_changed
//
void translations::SPBstepwheel_changed(){
    emit stepwheel_changed(SPBstepwheel->value());
}
//  End of function SPBstepwheel_changed

//  Function translation_changed
//
void translations::translation_changed(){
    if (RBTangstrom->isChecked()){
        emit setworldtranslation(QVector3D(SPBtras_x->value()*ANGSTROM_TO_BOHR,SPBtras_y->value()*ANGSTROM_TO_BOHR,
                                SPBtras_z->value()*ANGSTROM_TO_BOHR) );
    }
    else{
        emit setworldtranslation(QVector3D(SPBtras_x->value(),SPBtras_y->value(),SPBtras_z->value()) );
    }
//    emit setworldtranslation(translation);
}

//  End of function translation_changed
