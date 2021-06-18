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
//	Header file of screenshot class
//  Description: screenshot defines a widget that manages all the elements required for
//  capturing screenshots from 3D display
//
//	File:   screenshot.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef SCREENSHOT_H
#define SCREENSHOT_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QVector3D>

#include <QDebug>

class screenshot : public QWidget
{
    Q_OBJECT
public:
    explicit screenshot(QWidget *parent = 0);
    ~screenshot();

    bool getscaledef();
    bool gettransparentbkg();
    int getimagequality();
    QGroupBox * getFRMcapture();
    QString getCaptureFolder();
    QString getscalesize();

signals:
    void take_picture();
    void scaledef_changed();

private slots:
    void capturexyz(QVector3D);
    void create_capture_layouts();
    void RBTscaledef_changed();

private:

    QCheckBox *CHKtranspbg;

    QGroupBox *FRMcapture;
    QGroupBox *FRMresolution;
    QGroupBox *FRMscale;

    QLineEdit *TXTscalesize;

    QList<QMetaObject::Connection> connections;

    QPushButton *BTNshot;

    QRadioButton *RBTscaledef;              // Capture resolution scaled from screen resolution
    QRadioButton *RBTscreendef;             // Capture resolution equal to screen resolution

    QSpinBox *SPBimagequality;

    QString CaptureFolder;

};

#endif // SCREENSHOT_H
