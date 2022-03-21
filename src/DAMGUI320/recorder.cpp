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
//	Implementation of class recorder
//
//	File:   recorder.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QFileDialog>
#include <QSignalMapper>

#include "recorder.h"

recorder::recorder(QWidget *parent) : QWidget()
{

    connections.clear();

    isrecording = false;

    FRMrecord = new QGroupBox(tr("Record animation"));

    TXTrecordcommand = new QLineEdit("ffmpeg -y -framerate 30 -r 25 -i ");

    CHKremoveframes = new QCheckBox(tr("Remove frame files at end"),FRMrecord);
    CHKremoveframes->setChecked(true);

    recordfilename = molpath+"/film";

    TXTrecordfile = new QLineEdit(recordfilename);

    BTNrecordfile = new QToolButton(FRMrecord);
    BTNrecordfile->setText(tr("..."));
    connections << connect(BTNrecordfile, SIGNAL(clicked()), this, SLOT(importRecordFile()));

    SPBnumframes = new QSpinBox(FRMrecord);
    SPBnumframes->setMinimum(1);
    SPBnumframes->setMaximum(MAX_FRAMES);
    SPBnumframes->setSingleStep(10);

    BTNstartrecording = new QPushButton(QIcon(":/images/empezar.png"),tr("Start"));
    connections << connect(BTNstartrecording, SIGNAL(clicked()), this, SLOT(BTNstartrecording_clicked()));

    LBLmakingmovie = new QLabel("Making movie");
    LBLmakingmovie->setStyleSheet("color : blue");
    LBLmakingmovie->setVisible(false);

    create_record_layouts();     // Record menu layouts

}

recorder::~recorder()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

bool recorder::getisrecording(){
    return isrecording;
}

bool recorder::getrecord(){
    return record;
}

bool recorder::getremoveframes(){
    return CHKremoveframes->isChecked();
}

int recorder::getframeknt(){
    return frameknt;
}

int recorder::getnumframes(){
    return SPBnumframes->value();
}


QGroupBox * recorder::getFRMrecord(){
    return FRMrecord;
}

QString recorder::getrecordcommand(){
    return TXTrecordcommand->text();
}

QString recorder::getrecordfilename(){
    return recordfilename;
}

void recorder::BTNstartrecording_update(bool a){
    if (a) {
        BTNstartrecording->setText("Stop");
        LBLmakingmovie->setText("Recording frames");
        LBLmakingmovie->setVisible(!(LBLmakingmovie->isVisible()));
        repaint();
    }
    else{
        BTNstartrecording->setText("Start");
        LBLmakingmovie->setVisible(false);
        repaint();
    }
}

void recorder::BTNstartrecording_clicked(){
    if (!isrecording){
        isrecording = true;
        if (TXTrecordfile->text() != "")
            recordfilename = TXTrecordfile->text();
        else
            TXTrecordfile->setText(recordfilename);
        frameknt = 0;
        record = true;
        emit startrecording();
    }
    else{
        isrecording = false;
        emit stoprecording();
    }
}

//          Recording layouts

void recorder::create_record_layouts(){
    QLabel *LBLframesnumber= new QLabel(tr("Number of frames"));

    QGroupBox *FRMcommand = new QGroupBox(tr("Command for converting frames to film"));

    QHBoxLayout *layout4 = new QHBoxLayout(FRMcommand);
    layout4->addWidget(TXTrecordcommand);

    QGroupBox *FRMfile = new QGroupBox(tr("Name for record file"));

    QHBoxLayout *layout5 = new QHBoxLayout(FRMfile);
    layout5->addWidget(TXTrecordfile);
    layout5->addWidget(BTNrecordfile);

    QHBoxLayout *layout6 = new QHBoxLayout();
    layout6->addWidget(LBLframesnumber);
    layout6->addWidget(SPBnumframes);
    layout6->addStretch();

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addWidget(CHKremoveframes,Qt::AlignLeft);

    QHBoxLayout *layout8 = new QHBoxLayout();
    layout8->addStretch();
    layout8->addWidget(BTNstartrecording);
    layout8->addStretch();

    QHBoxLayout *layout9 = new QHBoxLayout();
    layout9->addWidget(LBLmakingmovie);
    layout9->addStretch();

    QVBoxLayout *layout10 = new QVBoxLayout(FRMrecord);
    layout10->addWidget(FRMcommand);
    layout10->addWidget(FRMfile);
    layout10->addLayout(layout6);
    layout10->addLayout(layout7);
    layout10->addLayout(layout8);
    layout10->addLayout(layout9);
}


void recorder::endmakingmovie(){
    if (LBLmakingmovie){
        LBLmakingmovie->setVisible(false);
    }
}

void recorder::importRecordFile(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString recordfilename = filedialog.getOpenFileName(this,tr("Open file ..."),molpath,
        tr("All files")+" (*)");
    if (recordfilename.length()==0) return;
    TXTrecordfile->setText(QFileInfo(QFile(recordfilename)).absolutePath()+"/"
            +QFileInfo(QFile(recordfilename)).completeBaseName());
    emit recordfilenamechanged(recordfilename);
}


void recorder::setnumframes(int i){
    SPBnumframes->setValue(i);
}

void recorder::setrecordfilename(QString filename){
    recordfilename = filename;
    TXTrecordfile->setText(recordfilename);
}

//void recorder::setmolecules(QList<molecule*> * mols){
//    molecules = mols;
//}

void recorder::setmolpath(QString a){
    molpath = a;
    recordfilename = molpath+"/film";
    TXTrecordfile->setText(recordfilename);
}

void recorder::setProjectFolder(QString a){
    ProjectFolder = a;
}

void recorder::startmakingmovie(){
    BTNstartrecording_clicked();
}
