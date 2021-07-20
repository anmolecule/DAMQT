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
#include <QComboBox>
#include <QDialog>
#include <QFileDialog>
#include <QFileInfo>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QProcess>
#include <QRadioButton>
#include <QSpinBox>
#include <QString>
#include <QTextEdit>
#include <QToolButton>


class Externals : public QWidget
{
    Q_OBJECT
public:
    explicit Externals(QWidget *parent = 0);
    ~Externals();

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
    
    
    void CMBengine_changed(int);
    
    void externalinputfile_changed();
    
    void external_geometry();
    
    void make_Gamess_input();
    void make_Gaussian_input();
    void make_Molpro_input();

    void RBTlocal_changed();

    void save_external_input();
    
    void submitError(QProcess::ProcessError);
    void submitOutput(int, QProcess::ExitStatus);
    void submitStart();
    
    void TXTextgeometry_changed();
        
private: 
    bool preview;

    int indexternal;
    
    QButtonGroup *QBGjobcommand;
    QButtonGroup *QBGrunmode;

    QComboBox *CMBbasis;
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
    QString extOutputFileName;


    QSpinBox *SPBcharge;
    QSpinBox *SPBmult;

    QTextEdit *extextEdit;
};
    


#endif    /* EXTERNALS_H */