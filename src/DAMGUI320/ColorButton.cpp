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
//
//  File:   ColorButton.cpp
//  Description: implements ColorButton class
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//

#include "ColorButton.h"

ColorButton::ColorButton(QWidget *parent) : QPushButton(parent) {
    bkgcolor = QColor(255,255,255);
    textcolor = QColor(0,0,0);
    setColor(&bkgcolor,&textcolor);
}
ColorButton::ColorButton(const QString & text, QWidget *parent) : QPushButton(text, parent) {
    bkgcolor = QColor(255,255,255);
    textcolor = QColor(0,0,0);
    setColor(&bkgcolor,&textcolor);
}
ColorButton::ColorButton(const QIcon & icon, const QString & text, QWidget *parent) : QPushButton(icon, text, parent) {
    bkgcolor = QColor(255,255,255);
    textcolor = QColor(0,0,0);
    setColor(&bkgcolor,&textcolor);
}

ColorButton::~ColorButton() {
}

QColor ColorButton::getbkgColor(){
    return bkgcolor;
}

QColor ColorButton::gettextColor(){
    return textcolor;
}

void ColorButton::setColor(QColor *color){
    QString strcolor = "background-color: rgb(";
    strcolor.append(QString("%1").arg(color->red()));
    strcolor.append(QString(",%1").arg(color->green()));
    strcolor.append(QString(",%1").arg(color->blue()));
    bkgcolor = *color;
	if ((color->red() > 210 || color->blue())&& color->green() > 210 ){
        strcolor.append("); color: rgb(0, 0, 0)");
        textcolor = QColor(0,0,0);
	}
	else{
        strcolor.append("); color: rgb(255, 255, 255)");
        textcolor = QColor(255, 255, 255);
	}
    setStyleSheet(strcolor);
    emit colorChanged(bkgcolor);
    emit textcolorChanged(textcolor);
}

void ColorButton::setColor(QColor *color, QColor *txtcolor){
    QString strcolor = "background-color: rgb(";
    strcolor.append(QString("%1").arg(color->red()));
    strcolor.append(QString(",%1").arg(color->green()));
    strcolor.append(QString(",%1").arg(color->blue()));
    strcolor.append("); color: rgb(").append(QString("%1,").arg(txtcolor->red()))
        .append(QString("%1,").arg(txtcolor->green())).append(QString("%1").arg(txtcolor->blue()))
		.append(")");
    bkgcolor = *color;
    textcolor = *txtcolor;
    setStyleSheet(strcolor);
    emit colorChanged(bkgcolor);
    emit textcolorChanged(textcolor);
}

