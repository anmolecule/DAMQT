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
//	Implementation of class axes
//
//	File:   axes.cpp
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

#include "axes.h"

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

axes::axes(QWidget *parent) : QWidget()
{
    axes_visible = false;
    axeslabels_visible = false;

    connections.clear();

    Xaxis_color = QColor(0,255,0);
    Yaxis_color = QColor(0,0,255);
    Zaxis_color = QColor(255,0,0);

    FRMaxes = new QGroupBox(tr("Axes"));
    FRMaxes->setVisible(false);

    CHKshowaxes = new QCheckBox(tr("Show axes"));
    CHKshowaxes->setChecked(axes_visible);
    connections << connect(CHKshowaxes, SIGNAL(stateChanged(int)), this, SLOT(CHKshowaxes_changed(int)));

    CHKshowaxeslabels = new QCheckBox(tr("Show axes labels"));
    CHKshowaxeslabels->setChecked(axes_visible && axeslabels_visible);
    CHKshowaxeslabels->setEnabled(axes_visible);
    connections << connect(CHKshowaxeslabels, SIGNAL(stateChanged(int)), this, SLOT(CHKshowaxeslabels_changed(int)));

//        Fonts
    fontaxeslabels = QFont("Noto Sans", 20, QFont::Bold);
    BTNfontaxeslabels = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNfontaxeslabels, SIGNAL(clicked()), this, SLOT(BTNfontaxeslabels_clicked()));

//    Colors

    BTNXaxiscolor = new ColorButton();
    BTNXaxiscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNXaxiscolor->setText(tr("X axis"));
    BTNXaxiscolor->setColor(&Xaxis_color);
    connections << connect(BTNXaxiscolor,SIGNAL(clicked()),this,SLOT(BTNXaxiscolor_clicked()));

    BTNYaxiscolor = new ColorButton();
    BTNYaxiscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNYaxiscolor->setText(tr("Y axis"));
    BTNYaxiscolor->setColor(&Yaxis_color);
    connections << connect(BTNYaxiscolor,SIGNAL(clicked()),this,SLOT(BTNYaxiscolor_clicked()));

    BTNZaxiscolor = new ColorButton();
    BTNZaxiscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNZaxiscolor->setText(tr("Z axis"));
    BTNZaxiscolor->setColor(&Zaxis_color);
    connections << connect(BTNZaxiscolor,SIGNAL(clicked()),this,SLOT(BTNZaxiscolor_clicked()));

//         Thickness

    SPBaxesthickness = new QSpinBox();
    SPBaxesthickness->setMinimum(1);
    SPBaxesthickness->setMaximum(50);
    SPBaxesthickness->setMaximumWidth(50);
    SPBaxesthickness->setValue(4);
    connections << connect(SPBaxesthickness,SIGNAL(valueChanged(int)),this,SLOT(SPBaxesthickness_changed(int)));

//         Vectors length

    SPBaxeslength = new QSpinBox();
    SPBaxeslength->setMinimum(1);
    SPBaxeslength->setMaximum(100);
    SPBaxeslength->setMaximumWidth(50);
    SPBaxeslength->setSingleStep(1);
    SPBaxeslength->setValue(10);
    connections << connect(SPBaxeslength,SIGNAL(valueChanged(int)),this,SLOT(SPBaxeslength_changed(int)));

//         Arrowlength

    SPBaxesarrowsize = new QSpinBox();
    SPBaxesarrowsize->setMinimum(1);
    SPBaxesarrowsize->setMaximum(50);
    SPBaxesarrowsize->setMaximumWidth(50);
    SPBaxesarrowsize->setSingleStep(1);
    SPBaxesarrowsize->setValue(4);
    connections << connect(SPBaxesarrowsize,SIGNAL(valueChanged(int)),this,SLOT(SPBaxesarrowssize_changed(int)));


//         Arrowwidth

    SPBaxesarrowwidth = new QSpinBox();
    SPBaxesarrowwidth->setMinimum(1);
    SPBaxesarrowwidth->setMaximum(50);
    SPBaxesarrowwidth->setMaximumWidth(50);
    SPBaxesarrowwidth->setSingleStep(1);
    SPBaxesarrowwidth->setValue(10);
    connections << connect(SPBaxesarrowwidth,SIGNAL(valueChanged(int)),this,SLOT(SPBaxesarrowswidth_changed(int)));

    create_axes_layouts();
}

axes::~axes()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

QGroupBox * axes::getFRMaxes(){
    return FRMaxes;
}

void axes::create_axes_layouts(){

    QLabel *LBLaxesthickness = new QLabel(tr("Axes thickness"));
    QLabel *LBLaxeslength = new QLabel(tr("Axes length"));
    QLabel *LBLaxesarrowsize = new QLabel(tr("Arrows length"));
    QLabel *LBLaxesarrowwidth = new QLabel(tr("Arrows width"));

    QVBoxLayout *layout1 = new QVBoxLayout();
    layout1->addWidget(CHKshowaxes);
    layout1->addWidget(CHKshowaxeslabels);
    layout1->addStretch();

    QHBoxLayout *layout2 = new QHBoxLayout();
    layout2->addStretch();
    layout2->addWidget(BTNfontaxeslabels);
    layout2->addStretch();

    QHBoxLayout *layout3 = new QHBoxLayout();
    layout3->addStretch();
    layout3->addWidget(BTNXaxiscolor);
    layout3->addWidget(BTNYaxiscolor);
    layout3->addWidget(BTNZaxiscolor);
    layout3->addStretch();

    QGridLayout *layout4 = new QGridLayout();
    layout4->addWidget(LBLaxesthickness,0,0);
    layout4->addWidget(SPBaxesthickness,0,1);
    layout4->addWidget(LBLaxeslength,1,0);
    layout4->addWidget(SPBaxeslength,1,1);
    layout4->addWidget(LBLaxesarrowsize,2,0);
    layout4->addWidget(SPBaxesarrowsize,2,1);
    layout4->addWidget(LBLaxesarrowwidth,3,0);
    layout4->addWidget(SPBaxesarrowwidth,3,1);

    QVBoxLayout *layout = new QVBoxLayout(FRMaxes);
    layout->addLayout(layout1);
    layout->addLayout(layout2);
    layout->addLayout(layout3);
    layout->addLayout(layout4);
    layout->addStretch();
}

void axes::BTNfontaxeslabels_clicked(){
    fontaxeslabels = QFontDialog::getFont(nullpointer, fontaxeslabels);
    emit axesfont(fontaxeslabels);
    return;
}

void axes::BTNXaxiscolor_clicked(){
    QColor currcolor = Xaxis_color;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        Xaxis_color = newcolor;
        BTNXaxiscolor->setColor(&newcolor);
    }
    emit xaxiscolor(Xaxis_color);
    return;
}

void axes::BTNYaxiscolor_clicked(){
    QColor currcolor = Yaxis_color;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        Yaxis_color = newcolor;
        BTNYaxiscolor->setColor(&newcolor);
    }
    emit yaxiscolor(Yaxis_color);
    return;
}

void axes::BTNZaxiscolor_clicked(){
    QColor currcolor = Zaxis_color;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        Zaxis_color = newcolor;
        BTNZaxiscolor->setColor(&newcolor);
    }
    emit zaxiscolor(Zaxis_color);
    return;
}


void axes::CHKshowaxes_changed(int a){
    if (a == 0){
        axes_visible = false;
        CHKshowaxeslabels->setChecked(false);
    }
    else{
        axes_visible = true;
    }
    CHKshowaxeslabels->setEnabled(axes_visible);
    emit axesvisible(axes_visible);
}

void axes::CHKshowaxeslabels_changed(int a){
    if (a == 0){
        axeslabels_visible = false;
    }
    else{
        axeslabels_visible = true;
    }
    emit axeslabelsvisible(axeslabels_visible);
}

void axes::SPBaxesarrowssize_changed(int i){
    emit axesarrowsize(i);
}

void axes::SPBaxesarrowswidth_changed(int i){
    emit axesarrowwidth(i);
}

void axes::SPBaxeslength_changed(int i){
    emit axeslength(i);
}

void axes::SPBaxesthickness_changed(int i){
    emit axesthickness(i);
}
