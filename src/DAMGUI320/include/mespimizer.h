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
#include <QTextCodec>

#include "GlobalInfo.h"

#include "ColorButton.h"

#include "molecule.h"
#include "elements.h"
#include "externals.h"

#define INTERVAL_INI 100
#define INTERVAL_SCALE 500


class chargescanvasDialog : public QDialog
{
    Q_OBJECT

    public:
        chargescanvasDialog(QVector <double> *, molecule *, QWidget *parent = 0);
        ~chargescanvasDialog();
        int iatomsel;
        int icpsel;
        int icptypesel;
        int imolsel;

    private slots:
        void BTNaccept_pressed();
        void BTNcancel_pressed();
        void charge_changed(QString);
        void updatecharges(int);

    private:

        QList<QMetaObject::Connection> connections;

        QString currentcharge;

        QVector <double> charges;
        QVector <double> *chargespntr;

};




class mespimizer : public QWidget
{
    Q_OBJECT
public:
    explicit mespimizer(QList<molecule*> *, QWidget *parent = 0);
    ~mespimizer();

    bool getdisplayenergy();
//    bool getlineinterpol();
    bool getmakingmovie();
    bool getoptimizeselect();
    bool getoptimizetemplate();
    bool getrecord();
    bool getremoveframes();

    int  getdelay();
    int  gethost();
    int  getinterpolpoints();
    int  gettemplate();

    QGroupBox * getFRMoptimizeCluster();
    QString getmespimizerpath();
    QString getrecordcommand();
    QString getrecordfilename();

    void enableBTNmespimize(bool);
    void endmakingmovie();
    void replacechargescanvas(int, QVector <double> *);
    void setBTNreplay(bool);
    void setBTNreset(bool);
    void setCHKoptimizeselect(bool);
    void setCHKrecordoptim(bool);
    void setchargescanvas(QVector < QVector <double> *> *);
    void setclustername(QString);
    void setdisplayEPIC(bool);
    void setenergycolor(QColor);
    void setenergyfont(QFont);
    void setenergyprecision(int);
    void setframesfile(QString);
    void setguestfromcanvas(bool);
    void sethartree(bool);
    void sethostname(QString);
    void setinterpolpoints(int);
    void setmespimizerpath(QString);
    void setmolecules(QList<molecule*> *);
    void setobabelcharges(bool);
    void setobabelindex(int);
    void setoptimizecanvas(bool);
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
    void displayepic_changed(bool);
    void energyprecision_changed(int);
    void exec_mespimizer();
    void font_clicked(QFont);
    void fontcolor_clicked(QColor);
    void framesfile_changed(QString);
    void hartree_units(bool);
    void interpol_changed(int);
    void movetotop();
    void obabelcharges_changed(bool);
    void obabelindex_changed(int);
    void optimizecanvas_changed(bool);
    void optimizeselect_changed(bool);
    void qmrun(QString);
    void recordfilenamechanged(QString);
    void recordoptim_changed(bool);
    void replay(QString);
    void reset(QString);
    void speed_changed(int);
    void startrecording();
    void stoptrecording();
    void template_changed(int);

private slots:
    void BTNfont_clicked();
    void BTNfontcolor_clicked();
    void BTNmespimize_clicked();
    void BTNreplay_clicked();
    void BTNreset_clicked();
    void canvaschargesvisible(bool);
    void CHKdisplayepic_changed();
    void CHKoptimizeselect_changed();
    void CHKrecordoptim_changed();
    void CMBobcharges_changed(int);
    void create_optimize_cluster_layouts();
    void editcharges(int);
    void editchargestemplate();
    void external_package();
    void framefiles_dialog();
    void importcharges(int);
    void importchargestemplate();
    void importchargesfromfile(QString,int);
    void importmltmod(QString,int);
    void importRecordDir();
    void mespath_dialog();
    void qmcomputing(QString);
    void RBTcharges_changed();
    void RBThartree_changed();
    void RBToptimizetemplate_changed();
    void SLDspeed_changed();
    void SPBenergyprecision_changed();
    void SPBhost_changed();
    void SPBinterpol_changed();
    void SPBtemplate_changed();
    void TXTclusterfile_changed();
    void TXTframesfile_changed();

private:
    bool checkobabelinstall();
    bool create_insertlocfile();
    bool create_mespimizer_input();
    bool create_preprocfile();
    bool create_templatefile();

    bool guestfromcanvas;
    bool obabelcharges;
    bool openbabelinstalled;

    int obabelindex;

    Elements *elem;

    ColorButton *BTNfontcolor;

    QCheckBox *CHKdisplayepic;
    QCheckBox *CHKoptimizeselect;
    QCheckBox *CHKrecordoptim;
    QCheckBox *CHKremoveframes;

    QColor energycolor;

    QComboBox *CMBobcharges;

    QDoubleSpinBox *SPBtssize;

    QFont energyfont;

    QGroupBox *FRManimation;
    QGroupBox *FRMcharges;
    QGroupBox *FRMtemplatecharges;
    QGroupBox *FRMenergy;
    QGroupBox *FRMinterpol;
    QGroupBox *FRMoptimizeCluster;
    QGroupBox *FRMoptimizeopt;
    QGroupBox *FRMrecord;
    QGroupBox *FRMtemplate;

    QLabel *LBLcanvascharges;
    QLabel *LBLdamchargestemplate;
    QLabel *LBLenergyprecision;
    QLabel *LBLhostindex;
    QLabel *LBLhostname;
    QLabel *LBLinterpol;
    QLabel *LBLmakingmovie;
    QLabel *LBLobcharges;
    QLabel *LBLtemplateindex;
    QLabel *LBLtemplatename;

    QLineEdit *TXTframesfile;
    QLineEdit *TXTclusterfile;
    QLineEdit *TXTmespimizerpath;
    QLineEdit *TXTrecordcommand;
    QLineEdit *TXTrecordir;
    QLineEdit *TXTrecordfile;

    QList<QMetaObject::Connection> connections;
    QList<molecule*> *molecules;

    QPushButton *BTNeditchargestemplate;
    QPushButton *BTNimportchargestemplate;
    QPushButton *BTNfont;
    QPushButton *BTNmespimize;
    QPushButton *BTNqm;
    QPushButton *BTNreplay;
    QPushButton *BTNreset;

    QRadioButton *RBThartree;
    QRadioButton *RBTkcalmol;
//    QRadioButton *RBTlininterpol;
    QRadioButton *RBToptimizecanvas;
    QRadioButton *RBToptimizetemplate;
//    QRadioButton *RBTquatinterpol;
    QRadioButton *RBTobabelcharges;
    QRadioButton *RBTusercharges;

    QSlider *SLDspeed;                    // Speed of animation

    QSpinBox *SPBenergyprecision;
    QSpinBox *SPBhost;
    QSpinBox *SPBinterpol;
    QSpinBox *SPBiterations;
    QSpinBox *SPBrssize;
    QSpinBox *SPBtemplate;

    QStringList obcharges;

    QToolButton *BTNmespath;
    QToolButton *BTNframefile;
    QToolButton *BTNrecordir;

//    QVector <double> chargestemplate;
    QVector < QLabel *> LBLdamchargescanvas;
    QVector < QPushButton *> BTNimportchargescanvas;
    QVector < QPushButton *> BTNeditchargescanvas;

    QVector < QVector <double> *> *chargescanvas;

};


#endif // MESPIMIZER_H
