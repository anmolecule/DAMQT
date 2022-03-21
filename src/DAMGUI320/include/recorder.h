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
//	Header file of recorder class
//  Description: recorder defines a widget that manages all the elements required for
//  recording animations
//
//	File:   recorder.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef RECORDER_H
#define RECORDER_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QList>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QString>
#include <QToolButton>

//#include "ColorButton.h"

#include "widgetsubclasses.h"

#define MAX_FRAMES 10000

class recorder : public QWidget
{
    Q_OBJECT
public:
    explicit recorder(QWidget *parent = 0);
    ~recorder();

    bool getisrecording();
    bool getrecord();
    bool getremoveframes();

    int getframeknt();
    int getnumframes();

    QGroupBox * getFRMrecord();

    QString getrecordcommand();
    QString getrecordfilename();

    void BTNstartrecording_update(bool);
    void endmakingmovie();
    void setmolpath(QString);
    void setProjectFolder(QString);
    void setrecordfilename(QString);
public slots:
    void setnumframes(int);
    void startmakingmovie();
signals:
    void recordfilenamechanged(QString);
    void startrecording();
    void stoprecording();

private slots:
    void BTNstartrecording_clicked();
    void create_record_layouts();
    void importRecordFile();

private:
    bool isrecording;
    bool record;

    int frameknt;

    QCheckBox *CHKremoveframes;

    QGroupBox *FRMrecord;

    QLabel *LBLmakingmovie;

    QLineEdit *TXTrecordcommand;
    QLineEdit *TXTrecordfile;

    QList<QMetaObject::Connection> connections;
//    QList<molecule*> *molecules;

    QPushButton *BTNstartrecording;

    QSpinBox *SPBnumframes;

    QString molpath;
    QString recordfilename;                 // Name for frames and film files
    QString ProjectFolder;

    QToolButton *BTNrecordfile;

};

#endif // RECORDER_H
