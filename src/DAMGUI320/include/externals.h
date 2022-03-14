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
//  File:   externals.h
//
//      Last version: July 2021
//
#ifndef EXTERNALS_H
#define EXTERNALS_H

#include <QButtonGroup>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QFileDialog>
#include <QFileInfo>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QListView>
#include <QMessageBox>
#include <QPushButton>
#include <QProcess>
#include <QRadioButton>
#include <QSpinBox>
#include <QStandardItemModel>
#include <QString>
#include <QTextEdit>
#include <QToolButton>

#include "elements.h"

class Externals : public QWidget
{
    Q_OBJECT
public:
    explicit Externals(QWidget *parent = 0);
    ~Externals();

    void setTXTextgeometry(QString);
    void setTXTextworkdir(QString);

protected:
//    virtual void reject();

signals:
    void closed();
    void computing(QString);
    void updatetextedit(QString);
    
    
private slots:
    void BTNjob_clicked();
    void BTNpreview_clicked();
    void BTNextreset_clicked();
    void BTNextsave_clicked();
    void BTNextsubmit_clicked();
    
    
    void CMBengine_changed();
    
    void externalinputfile_changed();
    
    void external_geometry();

    void formchkError(QProcess::ProcessError);
    void formchkOutput(int, QProcess::ExitStatus);
    void formchkStart();
    
    void make_Gamess_input();
    void make_Gamess_template();
    void make_Gaussian_input();
    void make_Gaussian_template();
    void make_Molpro_input();
    void make_Molpro_template();
    void make_Mopac_input();
    void make_Mopac_template();
    void make_NWChem_input();
    void make_NWChem_template();
    void make_Psi4_input();
    void make_Psi4_template();

    void RBTlocal_changed();

    void runformchk();

    void save_external_input();
    void make_sge_script();
    
    void submitError(QProcess::ProcessError);
    void submitOutput(int, QProcess::ExitStatus);
    void submitStart();
    
    void TXTextgeometry_changed();
        
private:
    void hideCMBlevel(QVector<int>);
    void resetCMBlevel();

    bool preview;
    
    QButtonGroup *QBGjobcommand;
    QButtonGroup *QBGrunmode;

    QCheckBox *CHKformchk;

    QComboBox *CMBbasis;
    QComboBox *CMBengine;
    QComboBox *CMBlevel;
    QComboBox *CMBlevel2;
    QComboBox *CMBtype;

    QDialog *QDLexternal;

    QGroupBox *FRMextproc;

    QLabel *LBLextproc;
    QLabel *LBLextpathremote;
    QLabel *LBLextworkdir;

    QLineEdit *TXTextcommand;
    QLineEdit *TXTextgeometry;
    QLineEdit *TXTextmem;
    QLineEdit *TXTextpathremote;
    QLineEdit *TXTexttime;
    QLineEdit *TXTextworkdir;
    QLineEdit *TXTkeywords;
    QLineEdit *TXTinputfile;
    QLineEdit *TXTtitle;

    QList<QMetaObject::Connection> connectionsext;

    QPushButton *BTNextsubmit;
    QPushButton *BTNjob;
    QPushButton *BTNpreview;

    QRadioButton *RBTlocal;
    QRadioButton *RBTremote;
    QRadioButton *RBTPBS;
    QRadioButton *RBTSGE;
    QRadioButton *RBTSLURM;

    QSpinBox *SPBextproc;

    QString extgeomfile;
    QString extgeompath;
    QString extInputFileName;
    QString extOutputSuffix;
    QString extOutputFileName;
    QString extJobscriptFileName;
    QString execname;

    QStringList extexecname;

    QSpinBox *SPBcharge;
    QSpinBox *SPBmult;

    QTextEdit *extextEdit;

};
    


#endif    /* EXTERNALS_H */
