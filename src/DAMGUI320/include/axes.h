//  Copyright 2008-2018, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
#ifndef AXES_H
#define AXES_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QPushButton>
#include <QSpinBox>

#include "ColorButton.h"

class axes : public QWidget
{
    Q_OBJECT
public:
    explicit axes(QWidget *parent = 0);
    ~axes();

    QGroupBox * getFRMaxes();

signals:
    void axesarrowsize(int);
    void axesarrowwidth(int);
    void axesfont(QFont);
    void axeslabelsvisible(bool);
    void axeslength(int);
    void axesthickness(int);
    void axesvisible(bool);
    void xaxiscolor(QColor);
    void yaxiscolor(QColor);
    void zaxiscolor(QColor);

private slots:
    void BTNfontaxeslabels_clicked();
    void BTNXaxiscolor_clicked();
    void BTNYaxiscolor_clicked();
    void BTNZaxiscolor_clicked();
    void CHKshowaxes_changed(int);
    void CHKshowaxeslabels_changed(int);
    void create_axes_layouts();
    void SPBaxesarrowssize_changed(int);
    void SPBaxesarrowswidth_changed(int);
    void SPBaxeslength_changed(int);
    void SPBaxesthickness_changed(int);

private:
    bool axes_visible;
    bool axeslabels_visible;

    ColorButton *BTNXaxiscolor;
    ColorButton *BTNYaxiscolor;
    ColorButton *BTNZaxiscolor;

    QCheckBox *CHKshowaxes;
    QCheckBox *CHKshowaxeslabels;

    QColor Xaxis_color;
    QColor Yaxis_color;
    QColor Zaxis_color;

    QFont fontaxeslabels;

    QGroupBox *FRMaxes;

    QList<QMetaObject::Connection> connections;

    QPushButton *BTNfontaxeslabels;

    QSpinBox *SPBaxesarrowsize;
    QSpinBox *SPBaxesarrowwidth;
    QSpinBox *SPBaxeslength;
    QSpinBox *SPBaxesthickness;
};

#endif // AXES_H
