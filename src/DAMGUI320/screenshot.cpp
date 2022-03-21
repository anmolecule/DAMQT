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
//	Implementation of class screenshot
//
//	File:   screenshot.cpp
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

#include "screenshot.h"

screenshot::screenshot(QWidget *parent) : QWidget()
{
    connections.clear();
    QDoubleValidator *myDoubleValidator = new QDoubleValidator(0);
    myDoubleValidator->setLocale(QLocale::English);
    FRMcapture = new QGroupBox(tr("Image capture"));
    FRMcapture->setVisible(false);

    CaptureFolder = "";
    FRMresolution = new QGroupBox(tr("Resolution"));


    FRMscale = new QGroupBox();
    FRMscale->setVisible(false);
    TXTscalesize = new QLineEdit(tr("1.0"));
    TXTscalesize->setValidator(myDoubleValidator);
    TXTscalesize->setAlignment(Qt::AlignRight);

    CHKtranspbg = new QCheckBox(tr("Transparent background"));

    BTNshot = new QPushButton(tr("Take picture"));
    connections << connect(BTNshot, SIGNAL(clicked()), this, SIGNAL(take_picture()));

    RBTscreendef = new QRadioButton();
    RBTscreendef->setChecked(true);
    RBTscreendef->setText(tr("Screen resolution"));
    RBTscaledef = new QRadioButton();
    RBTscaledef->setText(tr("Scaled resolution"));
    connections << connect(RBTscaledef, SIGNAL(toggled(bool)), this, SLOT(RBTscaledef_changed()));


    SPBimagequality = new QSpinBox();
    SPBimagequality->setMinimum(-1);
    SPBimagequality->setMaximum(100);
    SPBimagequality->setValue(20);
    SPBimagequality->setSingleStep(5);

    create_capture_layouts();    // Capture menu layouts
}

screenshot::~screenshot(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

bool screenshot::getscaledef(){
    return RBTscaledef->isChecked();
}

bool screenshot::gettransparentbkg(){
    return CHKtranspbg->isChecked();
}

int screenshot::getimagequality(){
    return SPBimagequality->value();
}

QGroupBox * screenshot::getFRMcapture(){
    return FRMcapture;
}

QString screenshot::getCaptureFolder(){
    return CaptureFolder;
}

QString screenshot::getscalesize(){
    return TXTscalesize->text();
}

void screenshot::capturexyz(QVector3D v){
    qDebug() << "captured position = " << v;
}

//          Capture layouts

void screenshot::create_capture_layouts(){
    QLabel *LBLscaledef= new QLabel(tr("Scale factor"));
    QLabel *LBLimagequality= new QLabel(tr("Image quality"));

    QHBoxLayout *layout1 = new QHBoxLayout(FRMscale);
    layout1->addWidget(LBLscaledef);
    layout1->addWidget(TXTscalesize);
    layout1->addStretch();

    QVBoxLayout *layout2 = new QVBoxLayout(FRMresolution);
    layout2->addStretch();
    layout2->addWidget(RBTscreendef);
    layout2->addWidget(RBTscaledef);
    layout2->addWidget(FRMscale);
    layout2->addStretch();

    QHBoxLayout *layout3 = new QHBoxLayout();
    layout3->addWidget(LBLimagequality);
    layout3->addWidget(SPBimagequality);
    layout3->addStretch();

    QVBoxLayout *layout4 = new QVBoxLayout(FRMcapture);
    layout4->addWidget(FRMresolution);
    layout4->addWidget(CHKtranspbg);
    layout4->addLayout(layout3);
    layout4->addWidget(BTNshot);
    layout4->addStretch();
}

void screenshot::RBTscaledef_changed(){
    if (RBTscaledef->isChecked())
        FRMscale->setVisible(true);
    else
        FRMscale->setVisible(false);
    FRMresolution->adjustSize();
    FRMcapture->adjustSize();
    emit scaledef_changed();
}
