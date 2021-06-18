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
//	Header file of mespimizer class
//  Description: mespimizer defines a widget that manages all the elements required for
//  cluster optimization with mesp (mespimizer)
//
//	File:   mespimizer.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef MESPIMIZER_H
#define MESPIMIZER_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QGroupBox>
#include <QLabel>
#include <QList>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>

#include "GlobalInfo.h"

#include "molecule.h"

#define INTERVAL_INI 100
#define INTERVAL_SCALE 500

class mespimizer : public QWidget
{
    Q_OBJECT
public:
    explicit mespimizer(QWidget *parent = 0);
    ~mespimizer();

//    bool getlineinterpol();
    bool getmakingmovie();
    bool getoptimizeselect();
    bool getoptimizetemplate();
    bool getqm();
    bool getrecord();
    bool getremoveframes();

    int  getdelay();
    int  gethost();
    int  getinterpolpoints();
    int  gettemplate();

    QGroupBox * getFRMoptimizeCluster();
    QString getmespimizerpath();
    QString getqmkeywords();
    QString getqmsoft();
    QString getrecordcommand();
    QString getrecordfilename();

    void enableBTNmespimize(bool);
    void endmakingmovie();
    void setBTNreset(bool);
    void setCHKoptimizeselect(bool);
    void setCHKrecordoptim(bool);
    void setclustername(QString);
    void setframesfile(QString);
    void setguessfromcanvas(bool);
    void sethostname(QString);
    void setinterpolpoints(int);
    void setmespimizerpath(QString);
    void setmolecules(QList<molecule*> *);
    void setrecordfile(QString);
    void setspeed(int);
    void setSPBhostmax(int);
    void setSPBtemplate(int);
    void setSPBtemplatemax(int);
    void setTXTframesfile(QString);
    void startmakingmovie();

signals:
    void adjustQDL();
    void clusterfile_changed(QString);
    void deletecluster();
    void exec_mespimizer();
    void framesfile_changed(QString);
    void interpol_changed(int);
    void movetotop();
    void optimizecanvas_changed(bool);
    void optimizeselect_changed(bool);
    void recordfilenamechanged(QString);
    void recordoptim_changed(bool);
    void replay(QString);
    void reset(QString);
    void speed_changed(int);
    void startrecording();
    void stoptrecording();
    void template_changed(int);

private slots:
    void BTNmespimize_clicked();
    void BTNreplay_clicked();
    void BTNreset_clicked();
    void CHKoptimizeselect_changed();
    void CHKqm_changed();
    void CHKrecordoptim_changed();
    void CMBqmsoft_changed(int);
    void create_optimize_cluster_layouts();
    void framefiles_dialog();
    void importRecordDir();
    void RBToptimizecanvas_changed();
    void RBToptimizetemplate_changed();
    void SLDspeed_changed();
    void SPBhost_changed();
    void SPBinterpol_changed();
    void SPBtemplate_changed();
    void TXTclusterfile_changed();
    void TXTframesfile_changed();

private:
    bool create_insertlocfile();
    bool create_mespimizer_input();
    bool create_preprocfile();
    bool create_templatefile();

    bool guestfromcanvas;

    QCheckBox *CHKoptimizeselect;
    QCheckBox *CHKqm;
    QCheckBox *CHKrecordoptim;
    QCheckBox *CHKremoveframes;

    QComboBox *CMBqmsoft;

    QDoubleSpinBox *SPBtssize;

    QGroupBox *FRManimation;
    QGroupBox *FRMinterpol;
    QGroupBox *FRMoptimizeCluster;
    QGroupBox *FRMoptimizeopt;
    QGroupBox *FRMqm;
    QGroupBox *FRMrecord;
    QGroupBox *FRMtemplate;

    QLabel *LBLhostindex;
    QLabel *LBLhostname;
    QLabel *LBLinterpol;
    QLabel *LBLmakingmovie;
    QLabel *LBLqm;
    QLabel *LBLqmsoft;
    QLabel *LBLtemplateindex;
    QLabel *LBLtemplatename;

    QLineEdit *TXTframesfile;
    QLineEdit *TXTclusterfile;
    QLineEdit *TXTmespimizerpath;
    QLineEdit *TXTqmkeywords;
    QLineEdit *TXTrecordcommand;
    QLineEdit *TXTrecordir;
    QLineEdit *TXTrecordfile;

    QList<QMetaObject::Connection> connections;
    QList<molecule*> *molecules;

    QPushButton *BTNmespimize;
    QPushButton *BTNreplay;
    QPushButton *BTNreset;

    QRadioButton *RBTlininterpol;
    QRadioButton *RBToptimizecanvas;
    QRadioButton *RBToptimizetemplate;
    QRadioButton *RBTquatinterpol;

    QSlider *SLDspeed;                    // Speed of animation

    QSpinBox *SPBhost;
    QSpinBox *SPBinterpol;
    QSpinBox *SPBrssize;
    QSpinBox *SPBtemplate;

    QString qmsoft;

    QToolButton *BTNframefile;
    QToolButton *BTNrecordir;

};

#endif // MESPIMIZER_H
