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
//	Implementation of class mespimizer
//
//	File:   mespimizer.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QFontDialog>
#include <QLabel>
#include <QPoint>
#include <QVBoxLayout>
#include <QCloseEvent>
#include <QColorDialog>
#include <QSignalMapper>

#include "mespimizer.h"

mespimizer::mespimizer(QWidget *parent) : QWidget()
{
    FRMoptimizeCluster = new QGroupBox(tr("Optimize cluster (EPIC)"));
    FRMoptimizeCluster->setVisible(false);

    FRMtemplate = new QGroupBox();
    FRMtemplate->setVisible(false);

    SPBhost = new QSpinBox();
    SPBhost->setMinimum(1);
    SPBhost->setMaximum(1);
    SPBhost->setMaximumWidth(50);
    SPBhost->setSingleStep(1);
    SPBhost->setValue(1);
    SPBhost->setEnabled(false);
    connections << connect(SPBhost,SIGNAL(valueChanged(int)),this,SLOT(SPBhost_changed()));

    SPBtemplate = new QSpinBox();
    SPBtemplate->setMinimum(1);
    SPBtemplate->setMaximum(1);
    SPBtemplate->setMaximumWidth(50);
    SPBtemplate->setSingleStep(1);
    SPBtemplate->setValue(1);
    SPBtemplate->setVisible(false);
    connections << connect(SPBtemplate,SIGNAL(valueChanged(int)),this,SLOT(SPBtemplate_changed()));

    CHKoptimizeselect = new QCheckBox();
    CHKoptimizeselect->setChecked(false);
    CHKoptimizeselect->setText(tr("Only selected CPs"));
    CHKoptimizeselect->setVisible(false);
    connections << connect(CHKoptimizeselect, SIGNAL(stateChanged(int)), this, SLOT(CHKoptimizeselect_changed()));

    FRMoptimizeopt = new QGroupBox();

    guestfromcanvas = true;
    RBToptimizecanvas = new QRadioButton(tr("Choose guests from canvas"),FRMoptimizeopt);
    RBToptimizecanvas->setCheckable(true);
    RBToptimizecanvas->setChecked(guestfromcanvas);
    connections << connect(RBToptimizecanvas,SIGNAL(toggled(bool)),this,SLOT(RBToptimizecanvas_changed()));

    RBToptimizetemplate = new QRadioButton(tr("Use template for guests"),FRMoptimizeopt);
    RBToptimizetemplate->setChecked(!guestfromcanvas);
    RBToptimizetemplate->setCheckable(true);
    connections << connect(RBToptimizetemplate,SIGNAL(toggled(bool)),this,SLOT(RBToptimizetemplate_changed()));

    SPBtssize = new QDoubleSpinBox();
    SPBtssize->setMinimum(0.1);
    SPBtssize->setMaximum(1.);
    SPBtssize->setMaximumWidth(55);
    SPBtssize->setSingleStep(0.1);
    SPBtssize->setValue(0.5);
    SPBtssize->setDecimals(1);
    SPBtssize->setEnabled(true);

    SPBrssize = new QSpinBox();
    SPBrssize->setMinimum(5);
    SPBrssize->setMaximum(30);
    SPBrssize->setMaximumWidth(55);
    SPBrssize->setSingleStep(5);
    SPBrssize->setValue(20);
    SPBrssize->setEnabled(true);

    SLDspeed = new QSlider(Qt::Horizontal);
    SLDspeed->setRange(0,INTERVAL_SCALE);
    SLDspeed->setSingleStep(1);
    SLDspeed->setPageStep(10);
    SLDspeed->setTickPosition(QSlider::TicksBelow);
    SLDspeed->setValue(INTERVAL_INI);
    connections << connect(SLDspeed,SIGNAL(sliderReleased()),this,SLOT(SLDspeed_changed()));

    // Buttons for choosing interpolation type are maitained here for possible debugging,
    // but actually they are not displayed nor used
//    FRMinterpol = new QGroupBox();

//    RBTlininterpol = new QRadioButton(tr("Linear interpolation"));
//    RBTlininterpol->setChecked(false);

//    RBTquatinterpol = new QRadioButton(tr("Quaternion interpolation"));
//    RBTquatinterpol->setChecked(true);
    // End of buttons for interpolation type

    SPBinterpol = new QSpinBox();
    SPBinterpol->setMinimum(1);
    SPBinterpol->setMaximum(50);
    SPBinterpol->setMaximumWidth(50);
    SPBinterpol->setSingleStep(1);
    SPBinterpol->setValue(10);
    connections << connect(SPBinterpol,SIGNAL(valueChanged(int)),this,SLOT(SPBinterpol_changed()));

    TXTmespimizerpath = new QLineEdit("");

    BTNmespath = new QToolButton();
    BTNmespath->setText(tr("..."));
    BTNmespath->setToolTip(tr("File with frames ..."));
    connections << connect(BTNmespath, SIGNAL(clicked()), this, SLOT(mespath_dialog()));

    TXTclusterfile = new QLineEdit("cluster");
    connections << connect(TXTclusterfile, SIGNAL(textChanged(const QString &)), this, SLOT(TXTclusterfile_changed()));


    BTNmespimize = new QPushButton(QIcon(":/images/exec.png"), tr("Exec"));
    BTNmespimize->setMaximumWidth(110);
    connections << connect(BTNmespimize,SIGNAL(clicked(bool)),this,SLOT(BTNmespimize_clicked()));

    FRManimation = new QGroupBox(tr("Animation"));

    TXTframesfile = new QLineEdit("");
    connections << connect(TXTframesfile, SIGNAL(textChanged(const QString &)), this, SLOT(TXTframesfile_changed()));

    BTNframefile = new QToolButton();
    BTNframefile->setText(tr("..."));
    BTNframefile->setToolTip(tr("File with frames ..."));
    connections << connect(BTNframefile, SIGNAL(clicked()), this, SLOT(framefiles_dialog()));


    CHKrecordoptim = new QCheckBox();
    CHKrecordoptim->setChecked(false);
    connections << connect(CHKrecordoptim, SIGNAL(stateChanged(int)), this, SLOT(CHKrecordoptim_changed()));

    FRMrecord = new QGroupBox();
    FRMrecord->setVisible(false);

    TXTrecordcommand = new QLineEdit("ffmpeg -y -framerate 30 -r 3 -i ");
    TXTrecordir = new QLineEdit();
    TXTrecordfile = new QLineEdit("film");

    BTNrecordir = new QToolButton(FRMrecord);
    BTNrecordir->setText(tr("..."));
    connections << connect(BTNrecordir, SIGNAL(clicked()), this, SLOT(importRecordDir()));

    CHKremoveframes = new QCheckBox(tr("Remove frame files at end"),FRMrecord);
    CHKremoveframes->setChecked(true);

    BTNreplay = new QPushButton(QIcon(":/images/replay.png"), tr("Replay"));
    BTNreplay->setMaximumWidth(110);
    BTNreplay->setEnabled(false);
    connections << connect(BTNreplay,SIGNAL(clicked(bool)),this,SLOT(BTNreplay_clicked()));

    BTNreset = new QPushButton(QIcon(""), tr("Reset"));
    BTNreset->setMaximumWidth(110);
    BTNreset->setEnabled(false);
    connections << connect(BTNreset,SIGNAL(clicked(bool)),this,SLOT(BTNreset_clicked()));

    LBLmakingmovie = new QLabel("Making movie");
    LBLmakingmovie->setStyleSheet("color : blue");
    LBLmakingmovie->setVisible(false);

    FRMenergy = new QGroupBox(tr("EPIC interaction energy"));

    CHKdisplayepic = new QCheckBox(tr("Display energy"));
    CHKdisplayepic->setChecked(true);
    connections << connect(CHKdisplayepic, SIGNAL(stateChanged(int)), this, SLOT(CHKdisplayepic_changed()));

    LBLenergyprecision = new QLabel(tr("Precision"));
    LBLenergyprecision->setVisible(true);

    SPBenergyprecision = new QSpinBox();
    SPBenergyprecision->setMinimum(1);
    SPBenergyprecision->setMaximum(8);
    SPBenergyprecision->setMaximumWidth(50);
    SPBenergyprecision->setSingleStep(1);
    SPBenergyprecision->setValue(4);
    SPBenergyprecision->setVisible(true);
    connections << connect(SPBenergyprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBenergyprecision_changed()));

    energyfont = QFont("Helvetica", 18, QFont::Bold);
    BTNfont = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    BTNfont->setVisible(true);
    connections << connect(BTNfont, SIGNAL(clicked()), this, SLOT(BTNfont_clicked()));

    energycolor = QColor(255, 172, 0, 255);
    BTNfontcolor = new ColorButton();
    BTNfontcolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNfontcolor->setText(tr("Color"));
    BTNfontcolor->setColor(&energycolor);
    BTNfontcolor->setEnabled(true);
    BTNfontcolor->setVisible(true);
    connections << connect(BTNfontcolor, SIGNAL(clicked()), this, SLOT(BTNfontcolor_clicked()));

    RBThartree = new QRadioButton(tr("Hartree"),FRMoptimizeopt);
    RBThartree->setCheckable(true);
    RBThartree->setChecked(false);
    connections << connect(RBThartree,SIGNAL(toggled(bool)),this,SLOT(RBThartree_changed()));

    RBTkcalmol = new QRadioButton(tr("kcal/mol"),FRMoptimizeopt);
    RBTkcalmol->setCheckable(true);
    RBTkcalmol->setChecked(true);
    connections << connect(RBTkcalmol,SIGNAL(toggled(bool)),this,SLOT(RBThartree_changed()));

    emit hartree_units(false);

    emit energyprecision_changed(SPBenergyprecision->value());

    emit displayepic_changed(CHKdisplayepic->isChecked());

    BTNqm = new QPushButton(tr("Quantum mechanics"));
    BTNqm->setToolTip(tr("Open template for quantum mechanics calculation"));
    connections << connect(BTNqm, SIGNAL(clicked()), this, SLOT(external_package()));

    create_optimize_cluster_layouts();  // Optimize cluster layouts
    adjustSize();
}

mespimizer::~mespimizer()
{
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

//  Function create_mespimizer_input: creates input file for cluster optimization
//

bool mespimizer::create_mespimizer_input(){
    QString mespimizerpath = getmespimizerpath();
    QString fileName = mespimizerpath + "/mespimizer.inp";
    QFile fileout(fileName);
    if (!fileout.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Unable to create input file for cluster optimization"
                "File %1 cannot be opened").arg(fileName)+QString(":\n%1.").arg(fileout.errorString()));
        emit movetotop();
        return false;
    }
    if (!fileout.isOpen()){
        fileout.open(QIODevice::Text | QIODevice::ReadWrite);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    buff.append(QString("$OPTIONS\n").toLatin1());
    if (getoptimizetemplate()){
        buff.append(QString("templatefile=\"cluster_template.xyz\"\n").toLatin1());
        buff.append(QString("insertlocfile=\"host-insertloc-cps.xyz\"\n").toLatin1());
    }
    else{
        buff.append(QString("preprocfile=\"cluster_geometry.xyz\"\n").toLatin1());
    }
    buff.append(QString("clustername=\"%1\"\n").arg(TXTclusterfile->text().trimmed()).toLatin1());
    buff.append(QString("nocharge=.true.\n").toLatin1());
    buff.append(QString("tssize=%1\n").arg(SPBtssize->value()).toLatin1());
    buff.append(QString("rssize=%1\n").arg(SPBrssize->value()).toLatin1());
    if (!mespimizerpath.isEmpty()){
        buff.append(QString("path=\"%1/\"\n").arg(mespimizerpath).toLatin1());
    }
    if (iswindows){
        buff.append(QString("iswindows=T\n").toLatin1());
    }
    buff.append(QString("$END\n").toLatin1());
    int host = gethost();
    buff.append(QString("\""+molecules->at(host-1)->getpath()+"/"+molecules->at(host-1)->getname()+"\"\n").toLatin1());
#if QT_VERSION < 0x050E00
    outfile << buff << endl;
#else
    outfile << buff << Qt::endl;
#endif
    fileout.close();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".xyz_init";
    QFile(fileName).remove();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".xyz_vis";
    QFile(fileName).remove();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".xyz_final";
    QFile(fileName).remove();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".frames";
    QFile(fileName).remove();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".xyz_frames";
    QFile(fileName).remove();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".curr_frame";
    QFile(fileName).remove();
    fileName = mespimizerpath + "/" + TXTclusterfile->text().trimmed() + ".kntframes";
    QFile(fileName).remove();
    emit movetotop();
    return true;
}

//  Function create_insertlocfile: creates insertlocfile file for cluster optimization
//

bool mespimizer::create_insertlocfile(){
    int host = gethost();
    QString cpsfileName = molecules->at(host-1)->getpath()+QString("/")+molecules->at(host-1)->getname()
            +QString("-cps-v.xyz");
    QFile fileinp(cpsfileName);
    if (!fileinp.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Unable to create insertloc file for cluster optimization\n"
            "File %1 cannot be opened").arg(cpsfileName)+QString(":\n%1.").arg(fileinp.errorString()));
        return false;
    }
    QString fileName = getmespimizerpath() + QString("/host-insertloc-cps.xyz");
    QFile fileout(fileName);
    if (!fileout.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Unable to create insertloc file for cluster optimization\n"
            "File %1 cannot be opened").arg(cpsfileName)+QString(":\n%1.").arg(fileout.errorString()));
        return false;
    }
    QTextStream in(&fileinp);
    QList<QString> inputlist;
    int knt = 0;
    while(!in.atEnd()) {
        QString line = in.readLine();
        if (line.contains("x")){
            if (getoptimizeselect()){
                int host = gethost();
                if (molecules->at(host-1)->cps->getcpsactive(0,knt)){
                    inputlist << line;
                }
            }
            else{
                inputlist << line;
            }
            knt++;
        }
    }
    fileinp.close();
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    buff.append(QString("%1\n\n").arg(knt).toLatin1());
    for (int i = 0 ; i < inputlist.length() ; i++){
        buff.append(QString("%1\n").arg(inputlist[i]).toLatin1());
    }
#if QT_VERSION < 0x050E00
    outfile << buff << endl;
#else
    outfile << buff << Qt::endl;
#endif
    fileout.close();
    emit movetotop();
    return true;
}

//  Function create_preprocfile: creates preprocfile file for cluster optimization
//

bool mespimizer::create_preprocfile(){
    if (molecules->length() < 2){
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Unable to create preproc file for cluster optimization\n"
            "Two molecules are necessary at least for clustering:\n"
            "First molecule will be the host"));
        emit movetotop();
        return false;
    }
    QString fileName = getmespimizerpath() + "/cluster_geometry.xyz";
    QFile fileout(fileName);
    if (!fileout.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Unable to create preproc file for cluster optimization\n"
            "File %1 cannot be opened").arg(fileName)+QString(":\n%1.").arg(fileout.errorString()));
        emit movetotop();
        return false;
    }
    if (!fileout.isOpen()){
        fileout.open(QIODevice::Text | QIODevice::ReadWrite);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    Elements *elem = new Elements();
//  molecules[0] is host and it is not included in preprocfile
    for (int i = 1 ; i < molecules->length() ; i++){
        buff.append(QString("%1\n").arg(molecules->at(i)->getnumatoms()).toLatin1());
        buff.append(QString("%1\n").arg(molecules->at(i)->getnumatoms()).toLatin1());

        QVector <int> znuc;
        QVector <QVector3D> xyz;
        znuc =  molecules->at(i)->getcharges();
        xyz = molecules->at(i)->getxyz();

        QMatrix4x4 m;
        m.rotate(molecules->at(0)->getrotation().conjugated());
        m.translate(-molecules->at(0)->gettranslation());
        m.translate(molecules->at(i)->gettranslation());
        m.rotate(molecules->at(i)->getrotation());
        for (int j = 0 ; j < molecules->at(i)->getnumatoms() ; j++){
            QVector3D p = m*(xyz[j]-QVector3D(0.,0.,molecules->at(i)->getz_trans_ini()));
            buff.append(QString("%1  %2  %3  %4\n").arg(elem->getSymbol(znuc[j])).arg(p[0]*BOHR_TO_ANGSTROM,2,'E',8)
                    .arg(p[1]*BOHR_TO_ANGSTROM,0,'E',8).arg(p[2]*BOHR_TO_ANGSTROM,0,'E',8));
        }
    }

#if QT_VERSION < 0x050E00
    outfile << buff << endl;
#else
    outfile << buff << Qt::endl;
#endif
    fileout.close();
    emit movetotop();
    delete elem;
    return true;
}

//  Function create_templatefile: creates template file for cluster optimization
//

bool mespimizer::create_templatefile(){
    if (gethost() == gettemplate()){
        QMessageBox msgBox;
        int host = gethost();
        msgBox.setInformativeText(QString(tr("Template and host molecules are the same: ")+molecules->at(host-1)->getname())
                   +QString(tr("\nDo you want to continue?")));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Yes);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::No){
            return false;
        }else if (ret == QMessageBox::Cancel){
            return false;
        }
    }
    QString fileName = getmespimizerpath() + "/cluster_template.xyz";
    QFile fileout(fileName);
    if (!fileout.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Unable to create template file for cluster optimization.\n"
            "File %1 cannot be opened").arg(fileName)+QString(":\n%1.").arg(fileout.errorString()));
        emit movetotop();
        return false;
    }
    if (!fileout.isOpen()){
        fileout.open(QIODevice::Text | QIODevice::ReadWrite);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    int temp = gettemplate();
    buff.append(QString("%1\n").arg(molecules->at(temp-1)->getnumatoms()).toLatin1());
    buff.append(QString("%1\n").arg(molecules->at(temp-1)->getnumatoms()).toLatin1());
    Elements *elem = new Elements();
    QVector <int> znuc;
    QVector <QVector3D> xyz;
    znuc =  molecules->at(temp-1)->getcharges();
    xyz = molecules->at(temp-1)->getxyz();
    for (int j = 0 ; j < molecules->at(temp-1)->getnumatoms() ; j++){
        QVector3D p = xyz[j];
        buff.append(QString("%1  %2  %3  %4\n").arg(elem->getSymbol(znuc[j])).arg(p[0]*BOHR_TO_ANGSTROM,2,'E',8)
                .arg(p[1]*BOHR_TO_ANGSTROM,0,'E',8).arg(p[2]*BOHR_TO_ANGSTROM,0,'E',8).toLatin1());
    }
#if QT_VERSION < 0x050E00
    outfile << buff << endl;
#else
    outfile << buff << Qt::endl;
#endif
    fileout.close();
    emit movetotop();
    delete elem;
    return true;
}

bool mespimizer::getdisplayenergy(){
    return CHKdisplayepic->isChecked();
}

//bool mespimizer::getlineinterpol(){
//    return RBTlininterpol->isChecked();
//}

bool mespimizer::getmakingmovie(){
    return LBLmakingmovie->isVisible();
}

bool mespimizer::getoptimizeselect(){
    return CHKoptimizeselect->isChecked();
}

bool mespimizer::getoptimizetemplate(){
    return RBToptimizetemplate->isChecked();
}

bool mespimizer::getrecord(){
    return CHKrecordoptim->isChecked();
}

bool mespimizer::getremoveframes(){
    return CHKremoveframes->isChecked();
}

int mespimizer::gethost(){
    return SPBhost->value();
}

int mespimizer::getdelay(){
    return INTERVAL_SCALE-SLDspeed->value()+1;
}

int mespimizer::getinterpolpoints(){
    return SPBinterpol->value();
}

int mespimizer::gettemplate(){
    return SPBtemplate->value();
}

QGroupBox * mespimizer::getFRMoptimizeCluster(){
    return FRMoptimizeCluster;
};

QString mespimizer::getmespimizerpath(){
    return TXTmespimizerpath->text();
}

QString mespimizer::getrecordcommand(){
    return TXTrecordcommand->text();
}

QString mespimizer::getrecordfilename(){
    if (TXTrecordir->text().isEmpty()){
        return TXTmespimizerpath->text()+"/"+TXTrecordfile->text();
    }
    else{
        return TXTrecordir->text().trimmed()+"/"+TXTrecordfile->text().trimmed();
    }
}

void mespimizer::BTNfont_clicked(){
    energyfont = QFontDialog::getFont(nullpointer, energyfont);
    emit font_clicked(energyfont);
}

void mespimizer::BTNfontcolor_clicked(){
    QColor col = QColorDialog::getColor(energycolor, this);
    if(col.isValid()) {
        BTNfontcolor->setColor(&col);
        emit fontcolor_clicked(col);
    }
}


//  Function BTNmespimize_clicked: executes cluster optimization
//

void mespimizer::BTNmespimize_clicked(){
    if (molecules->length() > 0 && molecules->last()->getiscluster())
        emit deletecluster();
    QString mespimizerpath =  TXTmespimizerpath->text();
    if (mespimizerpath.isEmpty()){
        mespimizerpath = molecules->at(SPBhost->value()-1)->getpath();
        TXTmespimizerpath->setText(mespimizerpath);
    }
    if (!QDir(mespimizerpath).exists()){
        if (!QDir().mkpath(mespimizerpath)){
            QMessageBox::warning(this, tr("MESPIMIZER"),tr("Could not create folder %1").arg(mespimizerpath));
            BTNmespimize->setEnabled(true);
            emit movetotop();
            return;
        }
    }
    if (RBToptimizetemplate->isChecked()){
        if (!(create_templatefile())){ 
            return;
        }
        if (!(create_insertlocfile())){
            return;
        }
    }
    else{
        if (!(create_preprocfile())){
            return;
        }
    }
    if (!(create_mespimizer_input())){
        return;
    }
    BTNmespimize->setEnabled(false);
    emit exec_mespimizer();
}

//  Function BTNreplay_clicked: replays previous cluster optimization
//

void mespimizer::BTNreplay_clicked(){
    emit replay(TXTframesfile->text().trimmed());
}

void mespimizer::BTNreset_clicked(){
    emit reset(TXTframesfile->text().trimmed());
}

void mespimizer::CHKdisplayepic_changed(){
    if (CHKdisplayepic->isChecked()){
        LBLenergyprecision->setVisible(true);
        SPBenergyprecision->setVisible(true);
        BTNfont->setVisible(true);
        BTNfontcolor->setVisible(true);
    }
    else{
        LBLenergyprecision->setVisible(false);
        SPBenergyprecision->setVisible(false);
        BTNfont->setVisible(false);
        BTNfontcolor->setVisible(false);
    }
    emit displayepic_changed(CHKdisplayepic->isChecked());
}

void mespimizer::CHKoptimizeselect_changed(){
    emit optimizeselect_changed(CHKoptimizeselect->isChecked());
}

void mespimizer::CHKrecordoptim_changed(){
    if (CHKrecordoptim->isChecked()){
        FRMrecord->setVisible(true);
    }
    else{
        FRMrecord->setVisible(false);
    }
    emit recordoptim_changed(CHKrecordoptim->isChecked());
    emit adjustQDL();
    emit movetotop();
}

//          Save optimize cluster layouts

void mespimizer::create_optimize_cluster_layouts(){

    LBLhostindex = new QLabel(tr("Host: "));
    LBLhostname = new QLabel();
    LBLinterpol = new QLabel(tr("Interpolation points: "));
    LBLtemplateindex = new QLabel(tr("Template: "));
    LBLtemplateindex->setVisible(false);
    LBLtemplatename = new QLabel();
    LBLtemplatename->setVisible(false);

    QLabel *LBLmespimizerpath = new QLabel(tr("Folder: "));
    QLabel *LBLclusterfile = new QLabel(tr("File: "));

    QHBoxLayout *layout1a = new QHBoxLayout();
    layout1a->addWidget(LBLmespimizerpath);
    layout1a->addWidget(TXTmespimizerpath);
    layout1a->addWidget(BTNmespath);

    QHBoxLayout *layout1b = new QHBoxLayout();
    layout1b->addWidget(LBLclusterfile,Qt::AlignLeft);
    layout1b->addWidget(TXTclusterfile,Qt::AlignRight);

    QVBoxLayout *layout1 = new QVBoxLayout();
    layout1->addLayout(layout1a);
    layout1->addLayout(layout1b);

    QGridLayout *layout3 = new QGridLayout();
    layout3->addWidget(LBLhostindex,0,0,Qt::AlignLeft);
    layout3->addWidget(LBLhostname,0,1,Qt::AlignRight);
    layout3->addWidget(SPBhost,0,2,Qt::AlignRight);

    QLabel *LBLtssize = new QLabel(tr("Translation stride (bohr): "));
    QHBoxLayout *layout4a = new QHBoxLayout();
    layout4a->addWidget(LBLtssize);
    layout4a->addWidget(SPBtssize);

    QLabel *LBLrssize = new QLabel(tr("Rotation stride (degrees): "));
    QHBoxLayout *layout4b = new QHBoxLayout();
    layout4b->addWidget(LBLrssize);
    layout4b->addWidget(SPBrssize);

    QVBoxLayout *layout5 = new QVBoxLayout(FRMoptimizeopt);
    layout5->addLayout(layout3);
    layout5->addWidget(RBToptimizecanvas);
    layout5->addWidget(RBToptimizetemplate);
    layout5->addLayout(layout4a);
    layout5->addLayout(layout4b);

    QGridLayout *layout6 = new QGridLayout(FRMtemplate);
    layout6->addWidget(CHKoptimizeselect,0,0,1,3,Qt::AlignLeft);
    layout6->addWidget(LBLtemplateindex,1,0,Qt::AlignLeft);
    layout6->addWidget(LBLtemplatename,1,1,Qt::AlignRight);
    layout6->addWidget(SPBtemplate,1,2,Qt::AlignRight);

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addStretch();
    layout7->addWidget(BTNmespimize,Qt::AlignCenter);
    layout7->addStretch();

    QLabel *LBLframefile = new QLabel(tr("File: "));
    QHBoxLayout *layout8 = new QHBoxLayout();
    layout8->addWidget(LBLframefile);
    layout8->addWidget(TXTframesfile);
    layout8->addWidget(BTNframefile);

    QLabel *LBLfast = new QLabel(tr("Fast"));
    QLabel *LBLslow = new QLabel(tr("Slow"));
    QHBoxLayout *layout9 = new QHBoxLayout();
    layout9->addWidget(LBLslow);
    layout9->addWidget(SLDspeed);
    layout9->addWidget(LBLfast);

    QHBoxLayout *layout11 = new QHBoxLayout();
    layout11->addWidget(LBLinterpol);
    layout11->addWidget(SPBinterpol);

    QLabel *LBLrecordoptim = new QLabel(tr("Record optimization"));
    QHBoxLayout *layout13 = new QHBoxLayout();
    layout13->addWidget(LBLrecordoptim);
    layout13->addWidget(CHKrecordoptim);

    QHBoxLayout *layout14 = new QHBoxLayout();
    layout14->addStretch();
    layout14->addWidget(BTNreset,Qt::AlignCenter);
    layout14->addWidget(BTNreplay,Qt::AlignCenter);
    layout14->addStretch();

    QGroupBox *FRMcommand = new QGroupBox(tr("Command for converting frames to film"));

    QHBoxLayout *layout15 = new QHBoxLayout(FRMcommand);
    layout15->addWidget(TXTrecordcommand);

    QGroupBox *FRMfile = new QGroupBox(tr("Record file"));
    QLabel *LBLrecordir = new QLabel(tr("Folder:"));
    QLabel *LBLrecordfile = new QLabel(tr("File:"));

    QHBoxLayout *layout16a = new QHBoxLayout();
    layout16a->addWidget(LBLrecordir);
    layout16a->addWidget(TXTrecordir);
    layout16a->addWidget(BTNrecordir);

    QHBoxLayout *layout16b = new QHBoxLayout();
    layout16b->addWidget(LBLrecordfile);
    layout16b->addWidget(TXTrecordfile);

    QVBoxLayout *layout16 = new QVBoxLayout(FRMfile);
    layout16->addLayout(layout16a);
    layout16->addLayout(layout16b);

    QHBoxLayout *layout18 = new QHBoxLayout();
    layout18->addWidget(CHKremoveframes,Qt::AlignLeft);

    QVBoxLayout *layout19 = new QVBoxLayout(FRMrecord);
    layout19->addWidget(FRMcommand);
    layout19->addWidget(FRMfile);
    layout19->addLayout(layout18);

    QVBoxLayout *layout20 = new QVBoxLayout(FRManimation);
    layout20->addLayout(layout8);
    layout20->addLayout(layout9);
    layout20->addLayout(layout11);
    layout20->addLayout(layout13);
    layout20->addWidget(FRMrecord);
    layout20->addLayout(layout14);

    QHBoxLayout *layout21 = new QHBoxLayout();
    layout21->addWidget(LBLmakingmovie);
    layout21->addStretch();

    QHBoxLayout *layout25 = new QHBoxLayout();
    layout25->addWidget(LBLenergyprecision);
    layout25->addWidget(SPBenergyprecision);

    QHBoxLayout *layout26 = new QHBoxLayout();
    layout26->addStretch();
    layout26->addWidget(BTNfont);
    layout26->addWidget(BTNfontcolor);
    layout26->addStretch();

    QHBoxLayout *layout27 = new QHBoxLayout();
    layout27->addWidget(RBTkcalmol);
    layout27->addWidget(RBThartree);

    QVBoxLayout *layout28 = new QVBoxLayout(FRMenergy);
    layout28->addWidget(CHKdisplayepic);
    layout28->addLayout(layout25);
    layout28->addLayout(layout26);
    layout28->addLayout(layout27);

    QHBoxLayout *layout29 = new QHBoxLayout();
    layout29->addStretch();
    layout29->addWidget(BTNqm);
    layout29->addStretch();

    QVBoxLayout *layout = new QVBoxLayout(FRMoptimizeCluster);
    layout->addStretch();
    layout->addLayout(layout1);
    layout->addWidget(FRMoptimizeopt);
    layout->addWidget(FRMtemplate);
    layout->addLayout(layout7);
    layout->addWidget(FRManimation);
    layout->addLayout(layout21);
    layout->addWidget(FRMenergy);
    layout->addLayout(layout29);
    layout->addStretch();
}

void mespimizer::enableBTNmespimize(bool a){
    BTNmespimize->setEnabled(a);
}

/**************************************************************************************************/
/********** FUNCTIONS FOR EXTERNAL PACKAGES DIALOG                      ***************************/
/**************************************************************************************************/

void mespimizer::external_package(){
    Externals *external = new Externals(this);
    connections << connect(external, SIGNAL(computing(QString)), this, SLOT(qmcomputing(QString)));
    connections << connect(external, SIGNAL(updatetextedit(QString)), this, SLOT(update_textedit(QString)));
    external->setTXTextgeometry(TXTmespimizerpath->text().trimmed()+"/"
        + TXTclusterfile->text().trimmed() +  ".xyz_final");
    external->setTXTextworkdir(TXTmespimizerpath->text().trimmed());
}

void mespimizer::framefiles_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(getmespimizerpath());
    QString fileName = filedialog.getOpenFileName(this,tr("Open file with frames"),getmespimizerpath(), tr("Allowed files")
            + "*.xyz_frames" + " (*.xyz_frames);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTframesfile->setText(fileName);
    TXTrecordfile->setText("film");
    BTNreplay->setEnabled(true);
    if (TXTmespimizerpath->text().isEmpty()){
        TXTmespimizerpath->setText(QFileInfo(fileName).canonicalPath());
    }
    emit reset(TXTframesfile->text().trimmed());
}

void mespimizer::mespath_dialog()
{
    QFileDialog dirdialog(this);
    dirdialog.setFileMode(QFileDialog::DirectoryOnly);
    dirdialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    dirdialog.setDirectory(getmespimizerpath());
    QString dirName = dirdialog.getExistingDirectory(this,tr("Open file with frames"));
    if (dirName.length()==0) return;
    TXTmespimizerpath->setText(dirName);

}

void mespimizer::importRecordDir(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString recordir = filedialog.getExistingDirectory(this, tr("Open Directory"),
                                                    QFileInfo(TXTrecordir->text()).absolutePath(),
                                                    QFileDialog::ShowDirsOnly
                                                    | QFileDialog::DontResolveSymlinks);
    TXTrecordir->setText(recordir);
    QString recordfilename = TXTrecordir->text().trimmed()+"/"+TXTrecordfile->text().trimmed();
    emit recordfilenamechanged(recordfilename);
}

void mespimizer::RBThartree_changed(){
    if (RBThartree->isChecked()){
        SPBenergyprecision->setValue(std::max(SPBenergyprecision->value()-1,2));
    }
    else{
        SPBenergyprecision->setValue(std::max(SPBenergyprecision->value(),4));
    }
    emit hartree_units(RBThartree->isChecked());
}

void mespimizer::RBToptimizecanvas_changed(){
    guestfromcanvas = RBToptimizecanvas->isChecked();
    emit optimizecanvas_changed(guestfromcanvas);
}

void mespimizer::RBToptimizetemplate_changed()
{
    if (!RBToptimizetemplate->isChecked()){
        LBLtemplateindex->setVisible(false);
        LBLtemplatename->setVisible(false);
        SPBtemplate->setVisible(false);
        FRMtemplate->setVisible(false);
        CHKoptimizeselect->setVisible(false);
        SPBhost->setValue(1);
        SPBhost->setEnabled(false);
    }
    else{
        if (molecules->length() == 0){
            QMessageBox::warning(this, tr("MESPIMIZER"),tr("At least one molecule must be loaded"));
            return;
        }
        LBLtemplateindex->setVisible(true);
        if (LBLtemplatename->text().isEmpty()){
            LBLtemplatename->setText(molecules->at(SPBtemplate->value()-1)->getname());
        }
        LBLtemplatename->setVisible(true);
        SPBhost->setMaximum(molecules->length());
        SPBtemplate->setMaximum(molecules->length());
        SPBtemplate->setVisible(true);
        FRMtemplate->setVisible(true);
        CHKoptimizeselect->setVisible(true);
        SPBhost->setEnabled(true);
    }
    emit adjustQDL();
    emit movetotop();
}

void mespimizer::endmakingmovie(){
    if (LBLmakingmovie){
        LBLmakingmovie->setVisible(false);
    }
    emit adjustQDL();
    emit movetotop();
}

void mespimizer::qmcomputing(QString a){
        emit qmrun(a);
}

void mespimizer::setCHKoptimizeselect(bool a){
    CHKoptimizeselect->setChecked(a);
}

void mespimizer::setBTNreplay(bool a){
    BTNreplay->setEnabled(a);
}

void mespimizer::setBTNreset(bool a){
    BTNreset->setEnabled(a);
}

void mespimizer::setCHKrecordoptim(bool a){
    CHKrecordoptim->setChecked(a);
}

void mespimizer::setclustername(QString a){
    TXTclusterfile->setText(a);
}

//  Function setdisplayEPIC: updates EPIC energy display flag
//
void mespimizer::setdisplayEPIC(bool a){
    CHKdisplayepic->setChecked(a);
}

//  Function setenergycolor: updates EPIC energy color
//
void mespimizer::setenergycolor(QColor a){
    energycolor = a;
    BTNfontcolor->setColor(&energycolor);
}

//  Function setenergyfont: updates EPIC energy font
//
void mespimizer::setenergyfont(QFont a){
    energyfont = a;
}

//  Function setenergyprecision: updates EPIC energy precision for display
//
void mespimizer::setenergyprecision(int a){
    SPBenergyprecision->setValue(a);
}

//  Function setenergyprecision: updates EPIC energy precision for display
//
void mespimizer::sethartree(bool a){
    RBThartree->setChecked(a);
    RBTkcalmol->setChecked(!a);
}

void mespimizer::setframesfile(QString a){
    TXTframesfile->setText(a);
}

void mespimizer::setguessfromcanvas(bool a){
    guestfromcanvas = a;
    RBToptimizecanvas->setChecked(a);
    update();
}

void mespimizer::sethostname(QString a){
    LBLhostname->setText(a);
}

void mespimizer::setinterpolpoints(int i){
    SPBinterpol->setValue(i);
}

void mespimizer::setmespimizerpath(QString a){
    TXTmespimizerpath->setText(a);
}

void mespimizer::setmolecules(QList<molecule*> *a){
    molecules = a;
}

void mespimizer::setspeed(int i){
    SLDspeed->setValue(i);
}

void mespimizer::setSPBhostmax(int i){
    SPBhost->setMaximum(i);
}

void mespimizer::setSPBtemplate(int i){
    SPBtemplate->setValue(i);
}

void mespimizer::setSPBtemplatemax(int i){
    SPBtemplate->setMaximum(i);
}

void mespimizer::SPBtemplate_changed(){
    if (molecules->length() > 0 ){
        LBLtemplatename->setText(molecules->at(SPBtemplate->value()-1)->getname());
    }
    else{
        LBLtemplatename->setText("");
    }
    int val = SPBtemplate->value();
    emit template_changed(val);
}

void mespimizer::setTXTframesfile(QString a){
    TXTframesfile->setText(a);
    if (!TXTframesfile->text().isEmpty()){
        BTNreplay->setEnabled(true);
        BTNreset->setEnabled(true);
    }
}

void mespimizer::setrecordfile(QString a){
    TXTrecordir->setText(QFileInfo(a).absolutePath());
    TXTrecordfile->setText(QFileInfo(a).fileName());
    if (TXTrecordfile->text().isEmpty()){
        CHKrecordoptim->setChecked(false);
    }
}

void mespimizer::SLDspeed_changed(){
    emit speed_changed(SLDspeed->value());
}

void mespimizer::SPBenergyprecision_changed(){
    emit energyprecision_changed(SPBenergyprecision->value());
}

void mespimizer::SPBhost_changed(){
    if (molecules->length() > 0 ){
        LBLhostname->setText(molecules->at(SPBhost->value()-1)->getname());
        TXTmespimizerpath->setText(molecules->at(SPBhost->value()-1)->getpath());
    }
    else{
        LBLhostname->setText("");
        TXTmespimizerpath->setText("");
    }
}

void mespimizer::SPBinterpol_changed(){
    emit interpol_changed(SPBinterpol->value());
}

void mespimizer::startmakingmovie(){
    LBLmakingmovie->setVisible(true);
    emit adjustQDL();
    emit movetotop();
}

void mespimizer::TXTclusterfile_changed(){
    emit clusterfile_changed(TXTclusterfile->text());
}

void mespimizer::TXTframesfile_changed(){
    if (QFileInfo(TXTframesfile->text()).exists()){
        BTNreset->setEnabled(true);
        BTNreplay->setEnabled(true);
    }
    else{
        BTNreset->setEnabled(false);
        BTNreplay->setEnabled(false);
    }
    emit framesfile_changed(TXTframesfile->text());
}

