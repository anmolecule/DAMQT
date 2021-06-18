//  Copyright 2008-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//  Several classes of dialogs
//
//  File:   dialog.cpp
//
//      Last version: December 2018
//

#include "dialog.h"

#include <QtDebug>

// //////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Class defining a dialog window with the user to decide among the
//  different density matrices available
//
// //////////////////////////////////////////////////////////////////////////////////////////////////////////

Dialog::Dialog(int nitems, QWidget *parent) : QDialog(parent)
{
	numitems = nitems;
	buttons = new bool[nitems];
	for (int i = 0 ; i < nitems ; i++){
		QRadioButton * radio_btn = new QRadioButton();
		RBToption << radio_btn; // append radio button to the list
	}
	QString path=QApplication::applicationDirPath();
	Label = new QLabel(tr("Select density matrix:"));
	BTNok = new QPushButton(tr("Accept"));
    BTNcancel = new QPushButton(tr("Cancel"));

	connect(BTNok, SIGNAL(clicked()), this, SLOT(setClose()));
    connect(BTNcancel, SIGNAL(clicked()), this, SLOT(setcancel()));
 }

//	Destructor
Dialog::~Dialog(){}
//	Value checked
void Dialog::setVal()
{
	for(int i=0;i<numitems;i++){
		if(RBToption[i]->isChecked()){
			val=i+1;
			break;
		}
	}
}
//	Buttons enabling	
void Dialog::RadioChange()
{
	for(int i=0;i<numitems;i++){
		RBToption[i]->setEnabled(buttons[i]);
		RBToption[i]->setVisible(buttons[i]);
	}
}
//	canceling
void Dialog::setcancel()
{
    cancel = true;
    this->close();
}
//	Closing
void Dialog::setClose()
{
    cancel = false;
    this->close();
}
//	Implement dialog
void Dialog::Implement()
{
	for(int i=0; i<numitems; i++){
		connect(RBToption[i], SIGNAL(toggled(bool)), this, SLOT(setVal()));
	}

	QHBoxLayout *layout1=new QHBoxLayout();
	layout1->addWidget(Label,0,Qt::AlignCenter);
	QVBoxLayout *layout2=new QVBoxLayout();
	for(int i=0; i<numitems; i++){
		layout2->addWidget(RBToption[i],0,Qt::AlignLeft);
	}
	QHBoxLayout *layout3=new QHBoxLayout();
	layout3->addWidget(BTNok,0,Qt::AlignCenter);
    layout3->addWidget(BTNcancel,0,Qt::AlignCenter);
	QVBoxLayout *layout0=new QVBoxLayout(this);
	layout0->addLayout(layout1,0);
	layout0->addLayout(layout2,0);
	layout0->addLayout(layout3,0);
		
	setWindowTitle(tr("Density matrix"));
}


// //////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class defining a dialog window with a label and a spinbox
// it is intended for choosing atoms or critical points to display or hide
//
// //////////////////////////////////////////////////////////////////////////////////////////////////////////

DialogSPB::DialogSPB(QList<molecule*> *mol, QWidget *parent) : QDialog(){
    molecules = mol;
    RBTcps[0] = new QRadioButton("(3,+3)");
    RBTcps[1] = new QRadioButton("(3,+1)");
    RBTcps[2] = new QRadioButton("(3,-1)");
    RBTcps[3] = new QRadioButton("(3,-3)");
    LBLatomsel = new QLabel(tr("Atom no."));
    SPBatomsel = new QSpinBox();
    connections.clear();

    molindex = 0;
    icptype = -1;
    QLabel *LBLmolecules = new QLabel(tr("Molecule"));
    CMBmolecules = new QComboBox();
    for (int i = 0 ; i < molecules->count() ; i++){
        CMBmolecules->addItem(molecules->at(i)->getname());
    }
    CMBmolecules->setCurrentIndex(molindex);
    connections << connect(CMBmolecules,SIGNAL(currentIndexChanged(int)),this,SLOT(CMBmolecules_changed(int)));
    natoms = molecules->at(molindex)->getnumatoms();

    FRMoption = new QGroupBox();
    RBTatoms = new QRadioButton(tr("Atoms"));   
    RBTatoms->setChecked(true);
    connect(RBTatoms, SIGNAL(toggled(bool)), this, SLOT(RBTtipo_changed()));
    RBTcpoints = new QRadioButton(tr("Critical points"));
    RBTcpoints->setChecked(false);
    connect(RBTcpoints, SIGNAL(toggled(bool)), this, SLOT(RBTtipo_changed()));

    FRMcps = new QGroupBox();
    FRMcps->setVisible(false);

    rbtchecked = false;
    if (molecules->at(molindex)->cps){
        for (int i = 0 ; i < MAX_CPS; i++){
            if (!rbtchecked && !molecules->at(molindex)->cps->cpsxyzval[i].isEmpty()){
                RBTcps[i]->setChecked(true);
                RBTcpoints->setVisible(true);
                icptype = i;
                rbtchecked = true;
            }
            else{
                RBTcps[i]->setChecked(false);
            }
            if (molecules->at(molindex)->cps->cpsxyzval[i].isEmpty()){
                RBTcps[i]->setVisible(false);
            }
            else{
                RBTcps[i]->setVisible(true);
            }
            connect(RBTcps[i], SIGNAL(toggled(bool)), this, SLOT(RBTcps_changed()));
        }
    }

    if (rbtchecked){
	   FRMoption->setVisible(true);
    }
    else{
	   FRMoption->setVisible(false);
       FRMcps->setVisible(false);
    }
    
	iatomsel = -1;
	icpsel = -1;
	this->setWindowTitle(tr("Show/hide"));
        setCursor(Qt::ArrowCursor);

	SPBatomsel->setMinimum(1);
	SPBatomsel->setMaximum(natoms);
	
	QPushButton *BTNaccept = new QPushButton(tr("Accept"));
	connect(BTNaccept, SIGNAL(clicked()), this, SLOT(setVal()));
	BTNaccept->setFocus();
	BTNaccept->setDefault(true);
	BTNaccept->setChecked(true);
	QPushButton *BTNcancel = new QPushButton(tr("Cancel"));
    connect(BTNcancel, SIGNAL(clicked()), this, SLOT(Cancel()));

    QVBoxLayout *layout1 = new QVBoxLayout();
    layout1->addWidget(LBLmolecules);
    layout1->addWidget(CMBmolecules);

    QHBoxLayout *layout2 = new QHBoxLayout(FRMoption);
    layout2->addWidget(RBTatoms);
    layout2->addWidget(RBTcpoints);

    QHBoxLayout *layout3 = new QHBoxLayout(FRMcps);
    layout3->addWidget(RBTcps[0]);
    layout3->addWidget(RBTcps[1]);
    layout3->addWidget(RBTcps[2]);
    layout3->addWidget(RBTcps[3]);
	
    QHBoxLayout *layout4 = new QHBoxLayout();
    layout4->addWidget(LBLatomsel);
    layout4->addWidget(SPBatomsel);

    QHBoxLayout *layout5 = new QHBoxLayout();
    layout5->addWidget(BTNcancel);
    layout5->addWidget(BTNaccept);

    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addLayout(layout1,Qt::AlignCenter);
    layout->addWidget(FRMoption);
    layout->addWidget(FRMcps);
    layout->addLayout(layout4,0);
    layout->addLayout(layout5,0);
    layout->addStretch();
}

//	Destructor
DialogSPB::~DialogSPB(){
    while(QApplication::overrideCursor() != 0)
            QApplication::restoreOverrideCursor();
}

//	Molecule index
void DialogSPB::CMBmolecules_changed(int i){
    molindex = i;
    natoms = molecules->at(molindex)->getnumatoms();
    LBLatomsel = new QLabel(tr("Atom no."));
    SPBatomsel->setMaximum(natoms);
    SPBatomsel->setValue(1);
    RBTatoms->setChecked(true);
    rbtchecked = false;
    if (molecules->at(molindex)->cps){
        for (int i = 0 ; i < MAX_CPS; i++){
            if (!rbtchecked && !molecules->at(molindex)->cps->cpsxyzval[i].isEmpty()){
                RBTcps[i]->setChecked(true);
                icptype = i;
                rbtchecked = true;
            }
            else{
                RBTcps[i]->setChecked(false);
            }
            if (molecules->at(molindex)->cps->cpsxyzval[i].isEmpty()){
                RBTcps[i]->setVisible(false);
            }
            else{
                RBTcps[i]->setVisible(true);
            }
        }
    }
    if (rbtchecked){
       FRMoption->setVisible(true);
    }
    else{
       FRMoption->setVisible(false);
       FRMcps->setVisible(false);
    }
    this->adjustSize();
}

//	Type of center: atom or critical point
void DialogSPB::RBTtipo_changed()
{
    if (RBTatoms->isChecked()){
       LBLatomsel->setText(tr("Atom no."));
       SPBatomsel->setMaximum(natoms);
       FRMcps->setVisible(false);
    }
    else{
       LBLatomsel->setText(tr("Critical point no."));
       SPBatomsel->setMaximum(molecules->at(molindex)->cps->cpsxyzval[icptype].count());
       if (rbtchecked) FRMcps->setVisible(true);
    }
    this->adjustSize();
}

//	Type of critical point
void DialogSPB::RBTcps_changed(){
    for (int i = 0 ; i < MAX_CPS; i++){
        if (RBTcps[i]->isChecked()){
            icptype = i;
            SPBatomsel->setMaximum(molecules->at(molindex)->cps->cpsxyzval[icptype].count());
            break;
        }
    }
}

//	Value checked
void DialogSPB::setVal()
{
    imolsel = molindex;
    if (RBTatoms->isChecked()){
        iatomsel = SPBatomsel->value()-1;
        icpsel = -1;
        icptypesel = -1;
        setClose();
    }
    else{
        icpsel = SPBatomsel->value()-1;
        icptypesel = icptype;
        iatomsel = -1;
        setClose();
    }
}

//	Cancel
void DialogSPB::Cancel()
{
    imolsel = -1;
    setClose();
}


//	Closing
void DialogSPB::setClose()
{
    this->close();
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class defining a dialog for displaying information
//
// //////////////////////////////////////////////////////////////////////////////////////////////////////////

DialogSPBinfo::DialogSPBinfo(QString *title, QString *info, QString *text, QWidget *parent) : QDialog(parent){
	FRMtitle = new QGroupBox();
	FRMtext = new QGroupBox();
	LBLtitle = new QLabel(*title);
	LBLtitle->setFont(QFont("Helvetica", 12, QFont::Bold));
	LBLtitle->setStyleSheet("QLabel { color : blue; }");
	if (LBLtitle->text().isEmpty()){ 
		FRMtitle->setHidden(true);
	}
	else{
		FRMtitle->setHidden(false);
	}
	LBLinfo = new QLabel(*info);
	LBLinfo->setFont(QFont("Helvetica", 12));
	if (LBLinfo->text().isEmpty()){ 
		LBLinfo->setHidden(true);
	}
	else{
		LBLinfo->setHidden(false);
	}
	LBLtext = new QLabel(*text);
	LBLtext->setFont(QFont("Helvetica", 12, QFont::Bold));
	if (LBLtext->text().isEmpty()){ 
		FRMtext->setHidden(true);
	}
	else{
		FRMtext->setHidden(false);
	}
	
	QPushButton *BTNclose = new QPushButton(tr("Close"));
	connect(BTNclose, SIGNAL(clicked()), this, SLOT(setClose()));
	
	
	QHBoxLayout *layout1=new QHBoxLayout(FRMtitle);
	layout1->addWidget(LBLtitle);
	layout1->setAlignment(Qt::AlignCenter);
	
	QHBoxLayout *layout2=new QHBoxLayout();
	layout2->addWidget(LBLinfo);
	layout2->setAlignment(Qt::AlignCenter);
	
	QHBoxLayout *layout3=new QHBoxLayout(FRMtext);
	layout3->addWidget(LBLtext);
	layout3->setAlignment(Qt::AlignCenter);
	QHBoxLayout *layout4=new QHBoxLayout();
	layout4->addWidget(BTNclose);
	layout4->setAlignment(Qt::AlignCenter);
	
	QVBoxLayout *layout10=new QVBoxLayout(this);
	layout10->addWidget(FRMtitle);
	layout10->addLayout(layout2);
	layout10->addWidget(FRMtext);
	layout10->addLayout(layout4);
}

//	Destructor
DialogSPBinfo::~DialogSPBinfo(){
    while(QApplication::overrideCursor() != 0)
            QApplication::restoreOverrideCursor();   
}
//	Closing
void DialogSPBinfo::setClose()
{
	this->close();
}

//	Set info section
void DialogSPBinfo::setinfo(QString *info)
{
	LBLinfo->setText(*info);
	if (LBLinfo->text().isEmpty()){ 
		LBLinfo->setHidden(true);
	}
	else{
		LBLinfo->setHidden(false);
	}
}

//	Set text section
void DialogSPBinfo::settext(QString *text)
{
	LBLtext->setText(*text);
	if (LBLtext->text().isEmpty()){ 
		FRMtext->setHidden(true);
	}
	else{
		FRMtext->setHidden(false);
	}
}

//	Set title section
void DialogSPBinfo::settitle(QString *title)
{
	LBLtitle->setText(*title);
	if (LBLtitle->text().isEmpty()){ 
		FRMtitle->setHidden(true);
	}
	else{
		FRMtitle->setHidden(false);
	}
}
