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
//  File:   ColorButton.h
//  Description: class ColorButton implements a QPushButton intended for color selection,
//  the color of the button is set to the color selected
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//

#ifndef COLORBUTTON_H
#define	COLORBUTTON_H

#include <QPushButton>
#include <QColor>
#include <QDebug>
//#include <QtGui>

class QColor;
class QPushButton;

class ColorButton : public QPushButton {
	Q_OBJECT
public:
	ColorButton(QWidget *parent = 0);
	ColorButton(const QString & text, QWidget * parent = 0 );
	ColorButton( const QIcon & icon, const QString & text, QWidget * parent = 0 );
	virtual ~ColorButton();
    QColor getbkgColor();
    QColor gettextColor();
public slots:
	void setColor(QColor *);
    void setColor(QColor *color, QColor *txtcolor);
signals:
    void colorChanged(QColor);
    void textcolorChanged(QColor);
private:
    QColor bkgcolor;
    QColor textcolor;
};

#endif	/* COLORBUTTON_H */

