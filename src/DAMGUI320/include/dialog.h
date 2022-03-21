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
//  File:   dialog.h
//
//      Last version: March 2016
//
#ifndef DIALOG_H
#define DIALOG_H

#include <QApplication>
#include <QDialog>
#include <QGroupBox>
#include <QLabel>
#include <QLayout>
#include <QObject>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>

#include "molecule.h"

class QErrorMessage;
class QGroupBox;
class QLabel;
class QPushButton;
class QRadioButton;
class QSpinBox;

class Dialog : public QDialog
{
	Q_OBJECT

	public:
		Dialog(int nitems = 4, QWidget *parent = 0);
		~Dialog();
        QSize sizeHint() const
        {
            return this->size();
        }
		bool *buttons;
        bool cancel;
		QList<QRadioButton*> RBToption;
		void RadioChange();
		int val;
		void Implement();

	private slots:
        void setcancel();
		void setVal();
		void setClose();		

	private:
		int numitems;
		QLabel *Label;
		QPushButton *BTNok;
        QPushButton *BTNcancel;
};


class DialogSPB : public QDialog
{
	Q_OBJECT

	public:
        DialogSPB(QList<molecule*> *molecules, QWidget *parent = 0);
		~DialogSPB();
		int iatomsel;
		int icpsel;
        int icptypesel;
        int imolsel;

	private slots:
        void Cancel();
        void CMBmolecules_changed(int);
		void setVal();
		void setClose();
        void RBTcps_changed();
		void RBTtipo_changed();
		
	private:
        bool rbtchecked;
        int icptype;
        int molindex;
		int natoms;
		int ncps;
        QComboBox *CMBmolecules;
        QGroupBox *FRMcps;
		QGroupBox *FRMoption;
		QLabel *LBLatomsel;
        QList<molecule*> *molecules;
        QList<QMetaObject::Connection> connections;     // Stores connections to release in destructor
		QRadioButton *RBTatoms;
        QRadioButton *RBTcpoints;
        QRadioButton *RBTcps[MAX_CPS];
		QSpinBox *SPBatomsel;
};

class DialogSPBinfo : public QDialog
{
	Q_OBJECT

	public:
		DialogSPBinfo(QString *title, QString *info, QString *text, QWidget *parent = 0);
		~DialogSPBinfo();
		
	public slots:	
		void setinfo(QString *info);
		void settext(QString *text);
		void settitle(QString *title);

	private slots:
		void setClose();
		
	private:
		QGroupBox *FRMtitle;
		QGroupBox *FRMtext;
		QLabel *LBLinfo;
		QLabel *LBLtext;
		QLabel *LBLtitle;		
};

#endif 
