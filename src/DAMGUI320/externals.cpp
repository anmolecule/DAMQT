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
//    Class defining the externals object. Dialog for launching external programs
//
//  File:   externals.cpp
//
//      Last version: July 2021
//
#include <QBoxLayout>
#include <QtDebug>

#include "externals.h"


Externals::Externals(QWidget *parent) : QWidget(parent)
{
    extOutputSuffix = "log";

    QDLexternal = new QDialog(this);
    QDLexternal->setMinimumWidth(400);
    QDLexternal->resize(900,200);
    QDLexternal->setWindowTitle(tr("External package"));
    QDLexternal->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::MinimumExpanding);

    extInputFileName = "";

    QLabel *LBLtitle = new QLabel(tr("Title"));
    TXTtitle = new QLineEdit(tr("Title"),QDLexternal);
    connectionsext << connect(TXTtitle, SIGNAL(returnPressed()), this, SLOT(externalinputfile_changed()));

    QLabel *LBLengine = new QLabel(tr("Engine choice"));
    CMBengine = new QComboBox(QDLexternal);
    CMBengine->addItem(tr("Gaussian"));
    CMBengine->addItem(tr("Gamess"));
    CMBengine->addItem(tr("Molpro"));
    CMBengine->addItem(tr("Mopac"));
    CMBengine->addItem(tr("NWChem"));
    CMBengine->addItem(tr("Psi4"));
    CMBengine->setCurrentIndex(0);
    connectionsext << connect(CMBengine, SIGNAL(currentIndexChanged(int)), this, SLOT(CMBengine_changed()));
    extexecname << "g09" << "gamess" << "molpro" << "mopac" << "nwchem" << "psi4";

    TXTextgeometry = new QLineEdit(QDLexternal);
    TXTextgeometry->setPlaceholderText(tr("Load file with molecule coordinates"));
    connectionsext << connect(TXTextgeometry, SIGNAL(textChanged(const QString &)), this, SLOT(TXTextgeometry_changed()));
    connectionsext << connect(TXTextgeometry, SIGNAL(returnPressed()), this, SLOT(externalinputfile_changed()));

    QToolButton *BTNgeometry = new QToolButton();
    BTNgeometry->setText(tr("..."));
    connectionsext << connect(BTNgeometry, SIGNAL(clicked()), this, SLOT(external_geometry()));

    QLabel *LBLtype = new QLabel(tr("Calculation type"));
    CMBtype = new QComboBox(QDLexternal);
    CMBtype->addItem(tr("Energy"));
    CMBtype->addItem(tr("Geometry optimization"));
    CMBtype->addItem(tr("Frequencies"));
    CMBtype->addItem(tr("Optimization+Frequencies"));
    CMBtype->addItem(tr("NMR"));
    connectionsext << connect(CMBtype, SIGNAL(currentIndexChanged(int)), this, SLOT(externalinputfile_changed()));

    CMBlevel = new QComboBox(QDLexternal);
    CMBlevel->setEditable(true);
    CMBlevel->addItem(tr("HF"));
    CMBlevel->addItem(tr("B3LYP"));
    CMBlevel->addItem("MP2");
    CMBlevel->addItem("CCSD");
    CMBlevel->addItem("BD");
    CMBlevel->addItem("CASSCF");
    CMBlevel->addItem(tr("PM6"));
    CMBlevel->addItem(tr("UFF"));
    connectionsext << connect(CMBlevel, SIGNAL(currentIndexChanged(int)), this, SLOT(externalinputfile_changed()));

    CMBlevel2 = new QComboBox(QDLexternal);
    CMBlevel2->addItem(tr("Default Spin"));
    CMBlevel2->addItem(tr("Restricted"));
    CMBlevel2->addItem(tr("Unrestricted"));
    CMBlevel2->addItem(tr("Open Restricted"));
    connectionsext << connect(CMBlevel2, SIGNAL(currentIndexChanged(int)), this, SLOT(externalinputfile_changed()));

    QLabel *LBLcharge = new QLabel(tr("Charge"));
    SPBcharge = new QSpinBox(QDLexternal);
    SPBcharge->setMaximum(30);
    SPBcharge->setMinimum(-10);
    SPBcharge->setValue(0);
    connectionsext << connect(SPBcharge,SIGNAL(valueChanged(int)),this,SLOT(externalinputfile_changed()));

    QLabel *LBLmult = new QLabel(tr("Mult"));
    SPBmult = new QSpinBox(QDLexternal);
    SPBmult->setMaximum(8);
    SPBmult->setMinimum(1);
    SPBmult->setValue(1);
    SPBmult->setToolTip("2S+1");
    connectionsext << connect(SPBmult,SIGNAL(valueChanged(int)),this,SLOT(externalinputfile_changed()));

    QLabel *LBLbasis = new QLabel(tr("Basis set"));
    CMBbasis = new QComboBox(QDLexternal);
    CMBbasis->setMinimumWidth(150);
    CMBbasis->setEditable(true);
    CMBbasis->addItem(tr("STO-3G"));
    CMBbasis->addItem(tr("3-21G"));
    CMBbasis->addItem(tr("6-31G"));
    CMBbasis->addItem(tr("6-311G"));
    CMBbasis->addItem(tr("cc-pVDZ"));
    CMBbasis->addItem(tr("cc-pVTZ"));
    CMBbasis->addItem(tr("cc-pVQZ"));
    CMBbasis->addItem(tr("cc-pV5Z"));
    connectionsext << connect(CMBbasis, SIGNAL(currentIndexChanged(int)), this, SLOT(externalinputfile_changed()));

    QLabel *LBLkeywords = new QLabel(tr("Keywords"));
    TXTkeywords = new QLineEdit(QDLexternal);
    TXTkeywords->setPlaceholderText("Additional Keywords");
    TXTkeywords->setText("5D 7F");
    connectionsext << connect(TXTkeywords, SIGNAL(textChanged(const QString &)), this, SLOT(externalinputfile_changed()));

    QBGrunmode = new QButtonGroup();

    RBTlocal = new QRadioButton(tr("Run Locally"),QDLexternal);
    RBTlocal->setChecked(true);
    RBTremote = new QRadioButton(tr("Run Remotely"),QDLexternal);
    RBTremote->setChecked(false);
    QBGrunmode->addButton(RBTlocal);
    QBGrunmode->addButton(RBTremote);

    BTNjob = new QPushButton(tr("Generate Job Script"),QDLexternal);
    BTNjob->setEnabled(false);
    BTNjob->setAutoDefault(false);
    connectionsext << connect(BTNjob, SIGNAL(clicked()), this, SLOT(BTNjob_clicked()));

    LBLextproc = new QLabel(tr("number of processors"),QDLexternal);
    SPBextproc = new QSpinBox(QDLexternal);
    SPBextproc->setMinimum(1);
    connectionsext << connect(SPBextproc,SIGNAL(valueChanged(int)),this,SLOT(externalinputfile_changed()));

    TXTextmem = new QLineEdit(QDLexternal);
    TXTextmem->setPlaceholderText(tr("Memory..."));
    connectionsext << connect(TXTextmem, SIGNAL(textChanged(const QString &)), this, SLOT(externalinputfile_changed()));

    TXTexttime = new QLineEdit(QDLexternal);
    TXTexttime->setPlaceholderText(tr("Time limit..."));

    QBGjobcommand = new QButtonGroup(QDLexternal);

    RBTPBS = new QRadioButton(tr("PBS"),QDLexternal);
    RBTPBS->setEnabled(false);
    RBTPBS->setChecked(true);

    RBTSGE = new QRadioButton(tr("SGE"),QDLexternal);
    RBTSGE->setEnabled(false);

    RBTSLURM = new QRadioButton(tr("SLURM"),QDLexternal);
    RBTSLURM->setEnabled(false);

    QBGjobcommand->addButton(RBTPBS);
    QBGjobcommand->addButton(RBTSGE);
    QBGjobcommand->addButton(RBTSLURM);

    connectionsext << connect(RBTlocal,SIGNAL(toggled(bool)),this,SLOT(RBTlocal_changed()));

    LBLextworkdir  = new QLabel(tr("Working Dir"),QDLexternal);
    TXTextworkdir = new QLineEdit(QDLexternal);
    TXTextworkdir->setPlaceholderText("Working directory...");
    connectionsext << connect(TXTextworkdir, SIGNAL(textChanged(const QString &)), this, SLOT(externalinputfile_changed()));

    LBLextpathremote = new QLabel(tr("Path"),QDLexternal);
    TXTextpathremote = new QLineEdit(QDLexternal);
    TXTextpathremote->setPlaceholderText("Path to remote executable...");
    LBLextpathremote->setEnabled(false);
    TXTextpathremote->setEnabled(false);

    preview = true;
    BTNpreview = new QPushButton(tr("Hide preview"),QDLexternal);
    BTNpreview->setStyleSheet("QPushButton {color: red;}");
    BTNpreview->setAutoDefault(false);
    connectionsext << connect(BTNpreview, SIGNAL(clicked()), this, SLOT(BTNpreview_clicked()));

    extextEdit = new QTextEdit(QDLexternal);
    extextEdit->setFocusPolicy(Qt::ClickFocus);

    QLabel *LBLextcommand = new QLabel(tr("Exec command:"));
    TXTextcommand = new QLineEdit("g09",QDLexternal);

    CHKformchk = new QCheckBox("Formchk",QDLexternal);
    CHKformchk->setChecked(true);
    CHKformchk->setVisible(true);

    QPushButton *BTNextreset = new QPushButton(tr("Reset"));
    BTNextreset->setAutoDefault(false);
    BTNextreset->setToolTip(tr("Delete all modifications introduced by hand in the editor"));
    connectionsext << connect(BTNextreset, SIGNAL(clicked()), this, SLOT(BTNextreset_clicked()));

    QPushButton *BTNextsave = new QPushButton(tr("Save"));
    BTNextsave->setAutoDefault(false);
    BTNextsave->setToolTip(tr("Saves the editor content to a file"));
    connectionsext << connect(BTNextsave, SIGNAL(clicked()), this, SLOT(BTNextsave_clicked()));

    BTNextsubmit = new QPushButton(tr("Submit"),QDLexternal);
    BTNextsubmit->setAutoDefault(false);
    BTNextsubmit->setToolTip(tr("Saves the editor content to an input file and submits the job"));
    connectionsext << connect(BTNextsubmit, SIGNAL(clicked()), this, SLOT(BTNextsubmit_clicked()));

    make_Gaussian_input();

    QHBoxLayout *layout1 = new QHBoxLayout();
    layout1->addWidget(LBLtitle);
    layout1->addWidget(TXTtitle);

    QHBoxLayout *layout2 = new QHBoxLayout();
    layout2->addWidget(LBLengine);
    layout2->addWidget(CMBengine);
    layout2->addWidget(TXTextgeometry);
    layout2->addWidget(BTNgeometry);

    QHBoxLayout *layout3 = new QHBoxLayout();
    layout3->addWidget(LBLtype);
    layout3->addWidget(CMBtype);
    layout3->addWidget(CMBlevel);
    layout3->addWidget(CMBlevel2);
    layout3->addWidget(LBLcharge);
    layout3->addWidget(SPBcharge);
    layout3->addWidget(LBLmult);
    layout3->addWidget(SPBmult);

    QHBoxLayout *layout4 = new QHBoxLayout();
    layout4->addWidget(LBLbasis);
    layout4->addWidget(CMBbasis);
    layout4->addWidget(LBLkeywords);
    layout4->addWidget(TXTkeywords);

    QHBoxLayout *layout5 = new QHBoxLayout();
    layout5->addWidget(RBTlocal);
    layout5->addWidget(RBTremote);

    QHBoxLayout *layout6 = new QHBoxLayout();
    layout6->addWidget(LBLextproc);
    layout6->addWidget(SPBextproc);

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addWidget(RBTPBS);
    layout7->addWidget(RBTSGE);
    layout7->addWidget(RBTSLURM);

    QGridLayout *layout10 = new QGridLayout();
    layout10->addLayout(layout5,0,0,1,4);
    layout10->addWidget(BTNjob,0,5,1,2);
    layout10->addLayout(layout6,1,0,1,2);
    layout10->addWidget(TXTextmem,1,2,1,1);
    layout10->addWidget(TXTexttime,1,3,1,1);
    layout10->addLayout(layout7,1,4,1,3);
    layout10->addWidget(LBLextworkdir,2,0,1,1);
    layout10->addWidget(TXTextworkdir,2,1,1,6);
    layout10->addWidget(LBLextpathremote,3,0,1,1);
    layout10->addWidget(TXTextpathremote,3,1,1,6);
    layout10->addWidget(LBLextcommand,4,0,1,1);
    layout10->addWidget(TXTextcommand,4,1,1,5);
    layout10->addWidget(CHKformchk,4,6,1,1);
    layout10->addWidget(BTNpreview,5,0,1,2);
    layout10->addWidget(BTNextreset,5,4,1,1);
    layout10->addWidget(BTNextsave,5,5,1,1);
    layout10->addWidget(BTNextsubmit,5,6,1,1);

    QVBoxLayout *layout = new QVBoxLayout(QDLexternal);
    layout->addLayout(layout1);
    layout->addLayout(layout2);
    layout->addLayout(layout3);
    layout->addLayout(layout4);
    layout->addLayout(layout10);
    layout->addWidget(extextEdit);

    QDLexternal->setTabOrder(TXTtitle,CMBengine);
    QDLexternal->setTabOrder(CMBengine,TXTextgeometry);
    QDLexternal->setTabOrder(TXTextgeometry,BTNgeometry);
    QDLexternal->setTabOrder(BTNgeometry,CMBtype);
    QDLexternal->setTabOrder(CMBtype,CMBlevel);
    QDLexternal->setTabOrder(CMBlevel,CMBlevel2);
    QDLexternal->setTabOrder(CMBlevel2,SPBcharge);
    QDLexternal->setTabOrder(SPBcharge,SPBmult);
    QDLexternal->setTabOrder(SPBmult,CMBbasis);
    QDLexternal->setTabOrder(CMBbasis,TXTkeywords);
    QDLexternal->setTabOrder(TXTkeywords,RBTlocal);
    QDLexternal->setTabOrder(RBTlocal,BTNjob);
    QDLexternal->setTabOrder(BTNjob,SPBextproc);
    QDLexternal->setTabOrder(SPBextproc,TXTextmem);
    QDLexternal->setTabOrder(TXTextmem,TXTexttime);
    QDLexternal->setTabOrder(TXTexttime,RBTPBS);
    QDLexternal->setTabOrder(RBTPBS,RBTSGE);
    QDLexternal->setTabOrder(RBTSGE,RBTSLURM);
    QDLexternal->setTabOrder(RBTSLURM,TXTextworkdir);
    QDLexternal->setTabOrder(TXTextworkdir,TXTextpathremote);
    QDLexternal->setTabOrder(TXTextpathremote,TXTextcommand);
    QDLexternal->setTabOrder(TXTextcommand,CHKformchk);
    QDLexternal->setTabOrder(CHKformchk,BTNpreview);
    QDLexternal->setTabOrder(BTNpreview,BTNextreset);
    QDLexternal->setTabOrder(BTNextreset,BTNextsave);
    QDLexternal->setTabOrder(BTNextsave,BTNextsubmit);
    TXTtitle->setFocus();

    QDLexternal->adjustSize();
    QDLexternal->show();
}

//	Destructor
Externals::~Externals(){
    for (int i = 0 ; i < connectionsext.size() ; i++){
        QObject::disconnect(connectionsext.at(i));
    }
    connectionsext.clear();
    delete QDLexternal;
}

void Externals::BTNjob_clicked(){
qDebug() << "BTNjob_clicked";
}

void Externals::BTNpreview_clicked(){
    preview = !preview;
    if (preview){
        BTNpreview->setText(tr("Hide preview"));
        BTNpreview->setStyleSheet("QPushButton {color: red;}");
        extextEdit->setVisible(true);
    }
    else{
        BTNpreview->setText(tr("Edit preview"));
        BTNpreview->setStyleSheet("QPushButton {color: black;}");
        extextEdit->setVisible(false);
    }
    QDLexternal->adjustSize();
    QDLexternal->update();
}

void Externals::BTNextreset_clicked(){
    if (CMBengine->currentText().toLower() == "gaussian"){
        CHKformchk->setVisible(true);
    }
    else{
        CHKformchk->setVisible(false);
    }
    switch (CMBengine->currentIndex()) {
        case 0:     // Gaussian
            make_Gaussian_input();
            break;
        case 1:     // Gamess
            make_Gamess_input();
            break;
        case 2:     // Molpro
            make_Molpro_input();
            break;
        case 3:     // Molpac
            make_Mopac_input();
            break;
        case 4:     // NWChem
            make_NWChem_input();
            break;
        case 5:     // Psi4
            make_Psi4_input();
            break;
    }
}

void Externals::BTNextsave_clicked(){
    save_external_input();
}

void Externals::BTNextsubmit_clicked(){
    QDir dir(TXTextworkdir->text().trimmed());
    if (!dir.exists()){
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Directory %1 does not exist")).arg(TXTextworkdir->text().trimmed())
                        +"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            if (!QDir().mkpath(TXTextworkdir->text().trimmed())){
                QMessageBox msgBox;
                msgBox.setText(tr("submitOutput"));
                msgBox.setInformativeText(QString(tr("Could not create directory %1\n")
                            .arg(TXTextworkdir->text().trimmed())));
                msgBox.setIcon(QMessageBox::Information);
                msgBox.exec();
                return;
            }
        }else{
            return;
        }
    }
    save_external_input();
    if (extInputFileName.isEmpty())
        return;
    QString strprocess;
    QStringList Parameters = TXTextcommand->text().trimmed().split(" ");
    QString processname = Parameters.at(0);
    Parameters.removeFirst();
    QString stdinput = extInputFileName;
    extOutputFileName = QFileInfo(extInputFileName).absolutePath()
        + "/" + QFileInfo(extInputFileName).completeBaseName() + "." + extOutputSuffix;
    QProcess *myProcess = new QProcess();
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(extOutputFileName,QIODevice::Truncate);
    myProcess->setStandardErrorFile(extOutputFileName,QIODevice::Append);
    connectionsext << connect(myProcess, SIGNAL(started()), this, SLOT(submitStart()));
    connectionsext << connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(submitError(QProcess::ProcessError)));
    connectionsext << connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(submitOutput(int, QProcess::ExitStatus)));
    myProcess->start(processname,Parameters);
//    QByteArray qbarray = myProcess->readAllStandardError();
}

void Externals::CMBengine_changed(){
    if (CMBengine->currentText().toLower() == "gaussian"){
        CHKformchk->setVisible(true);
    }
    else{
        CHKformchk->setVisible(false);
    }
    switch (CMBengine->currentIndex()) {
        case 0:     // Gaussian
            make_Gaussian_input();
            break;
        case 1:     // Gamess
            make_Gamess_input();
            break;
        case 2:     // Molpro
            make_Molpro_input();
            break;
        case 3:     // Molpac
            make_Mopac_input();
            break;
        case 4:     // NWChem
            make_NWChem_input();
            break;
        case 5:     // Psi4
            make_Psi4_input();
            break;
    }
}

void Externals::external_geometry(){
    QFileDialog filedialog;
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),"",
        tr("Import geometry from")+" (*.xyz *.xyz_init *.xyz_final);;"+
        tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTextgeometry->setText(fileName);
    QDLexternal->raise();
}

void Externals::externalinputfile_changed(){
    if (CMBengine->currentText().toLower() == "gaussian"){
        CHKformchk->setVisible(true);
    }
    else{
        CHKformchk->setVisible(false);
    }
    switch (CMBengine->currentIndex()) {
        case 0:     // Gaussian
            make_Gaussian_input();
            break;
        case 1:     // Gamess
            make_Gamess_input();
            break;
        case 2:     // Molpro
            make_Molpro_input();
            break;
        case 3:     // Molpac
            make_Mopac_input();
            break;
        case 4:     // NWChem
            make_NWChem_input();
            break;
        case 5:     // Psi4
            make_Psi4_input();
            break;
    }
    QDLexternal->raise();
}

void Externals::formchkError(QProcess::ProcessError error){
    emit computing(QString(""));
    QString message = QString("Error %1 ").arg(error)
                + QString(tr("Process failed to start program formchk\n"))
                + QString(tr("Check that it is installed in your system\n"))
                + QString(tr("If it is installed, add it to your $PATH and try again)"));
    QMessageBox::critical(this,QString("Error %1").arg(error),message);
    BTNextsubmit->setEnabled(true);
}

void Externals::formchkOutput(int exitCode, QProcess::ExitStatus status){
    if(status==QProcess::NormalExit){
        QMessageBox msgBox;
        msgBox.setText(tr("formchkOutput"));
        msgBox.setInformativeText(QString(tr("Computation ended\n")));
        msgBox.setIcon(QMessageBox::Information);
        msgBox.exec();

        emit computing(QString(""));
        emit updatetextedit(extOutputFileName);
    }else if (status==QProcess::CrashExit){
        QMessageBox msgBox;
        msgBox.setText(tr("formchkOutput"));
        msgBox.setInformativeText(QString(tr("Process crashed, exit code = %1").arg(exitCode).toLatin1()));
        msgBox.setIcon(QMessageBox::Information);
        msgBox.exec();
        emit computing(QString(""));
    }
    BTNextsubmit->setEnabled(true);
}

void Externals::formchkStart(){
    BTNextsubmit->setEnabled(false);
    emit computing(QString(tr("formchk launched...")));
}

void Externals::make_Gaussian_input(){
    extextEdit->clear();
    QStringList type = {"","opt","freq","opt freq","nmr=giao"};
    QStringList level2 = {"","r","u","ro"};
    QStringList mm = {"uff","amber","dreiding"};
    QStringList sme = {"pm6","pddg","am1","pm3","indo","cndo"};
    extOutputSuffix = "log";

    TXTextcommand->setText(extexecname[0]);
    QString filepath = TXTextworkdir->text().trimmed();
    if (filepath.isEmpty())
        return;
    QString filename = QFileInfo(TXTextgeometry->text()).baseName();
    if (filename.isEmpty())
        return;
    extInputFileName = filepath+"/"+filename+".com";
    TXTextcommand->setText(extexecname[0]+" "+extInputFileName);

    QFile geometryInput(TXTextgeometry->text().trimmed());
    if (!geometryInput.open(QFile::ReadOnly | QFile::Text)){
        return;
    }

    QByteArray buff;
    if (SPBextproc->value() > 1){
        buff.append(QString("\%nprocshared=%1\n").arg(SPBextproc->value()));
    }
    if (!TXTextmem->text().isEmpty()){
        buff.append(QString("\%mem=%1\n").arg(TXTextmem->text().trimmed()));
    }
    buff.append(QString("\%chk=%1.chk\n").arg(filepath+"/"+filename));
    buff.append(QString("#p %1  ").arg(type.at(CMBtype->currentIndex())));
    if (!(mm.contains(CMBlevel->currentText(),Qt::CaseInsensitive))){
        buff.append(QString("%1").arg(level2.at(CMBlevel2->currentIndex())));
    }
    buff.append(QString("%1").arg(CMBlevel->currentText().toLower()));
    if (mm.contains(CMBlevel->currentText(),Qt::CaseInsensitive) || sme.contains(CMBlevel->currentText(),Qt::CaseInsensitive)){
            buff.append(QString("\n\n"));   // Calculations without explicit basis set (MM and semiempirical)
    }
    else{
        buff.append(QString("/%1 %2\n\n").arg(CMBbasis->currentText()).arg(TXTkeywords->text().trimmed()));
    }
    buff.append(QString("%1\n\n").arg(TXTtitle->text().trimmed()));
    buff.append(QString("%1 %2\n").arg(SPBcharge->value()).arg(SPBmult->value()));

    QTextStream in(&geometryInput); // Buffer for reading from fileinput
    QString line = in.readLine();
#if QT_VERSION < 0x050E00
    QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
    QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
    int ncen = xyz.at(0).toInt();
    int kntcen = 0;
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xyz.count() == 4){
            buff.append(QString("%1    %2    %3    %4\n").arg(xyz[0]).arg(xyz[1]).arg(xyz[2]).arg(xyz[3]));
            kntcen++;
        }
    }
    if (ncen != kntcen){
        QMessageBox msgBox;
        msgBox.setText(tr("make_Gaussian_input"));
        msgBox.setInformativeText(tr("Wrong number of centers in file:\n")+
            TXTextgeometry->text().trimmed());
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    extextEdit->setText(buff);

}


void Externals::make_Gamess_input(){
    extextEdit->clear();
    extextEdit->setText("\n\n   To be prepared");
}

void Externals::make_Molpro_input(){
    extextEdit->clear();
    extextEdit->setText("\n\n   To be prepared");
    extOutputSuffix = "out";
}

void Externals::make_Mopac_input(){
    extextEdit->clear();
    extextEdit->setText("\n\n   To be prepared");
}

void Externals::make_NWChem_input(){
    extextEdit->clear();
    QStringList type = {"energy","optimize","freq","opt freq","nmr=giao"};
    QStringList level2 = {"","r","u","ro"};
    QStringList mult = {"","singlet","doublet","triplet","quartet","quintet","sextet","octet"};
    QStringList mm = {"scf","mp2","ccsd"};
    extOutputSuffix = "nwcout";

    TXTextcommand->setText(extexecname[4]);
    QString filepath = TXTextworkdir->text().trimmed();
    if (filepath.isEmpty())
        return;
    QString filename = QFileInfo(TXTextgeometry->text()).baseName();
    if (filename.isEmpty())
        return;
    extInputFileName = filepath+"/"+filename+".nwcinp";
    TXTextcommand->setText(extexecname[4]+" "+extInputFileName);

    QFile geometryInput(TXTextgeometry->text().trimmed());
    if (!geometryInput.open(QFile::ReadOnly | QFile::Text)){
        return;
    }

    QByteArray buff;
    buff.append(QString("start %1\n").arg(filename));
    buff.append(QString("title '%1'\n\n").arg(TXTtitle->text().trimmed()));
    buff.append(QString("scratch_dir %1\n\n").arg(TXTextworkdir->text().trimmed()));
    if (!TXTextmem->text().isEmpty()){
        buff.append(QString("memory total %1\n").arg(TXTextmem->text().trimmed()));
    }

    if (CMBlevel->currentIndex() < CMBlevel->count()-2){
        buff.append(QString("basis spherical\n * library %1\nend\n\n").arg(CMBbasis->currentText()));
    }
    if (level2.at(CMBlevel2->currentIndex()) != "") {
        buff.append(QString("scf \n%1hf\n%2\nend\n\n").arg(level2.at(CMBlevel2->currentIndex())).arg(mult.at(SPBmult->value())));
    }
    buff.append(QString("charge %1\n").arg(SPBcharge->value()));
    buff.append(QString("geometry units ang\n"));

    QTextStream in(&geometryInput); // Buffer for reading from fileinput

    QString line = in.readLine();
#if QT_VERSION < 0x050E00
    QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
    QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
    int ncen = xyz.at(0).toInt();
    int kntcen = 0;
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xyz.count() == 4){
            buff.append(QString("%1    %2    %3    %4\n").arg(xyz[0]).arg(xyz[1]).arg(xyz[2]).arg(xyz[3]));
            kntcen++;
        }
    }
    buff.append(QString("end\n"));
    // Please resolve QString error for level
    QString level;
    if (CMBlevel->currentText().toLower() == "hf") {
        level="scf";
    } else {
        level=CMBlevel->currentText().toLower();
    }
    if (mm.contains(level)){
        buff.append(QString("scf\n    vectors output %1\nend\n").arg(filepath+"/"+filename+".movecs"));
    }
    if (type.at(CMBtype->currentIndex()) == "energy") {  
        buff.append(QString("task %1\n").arg(QString("%1").arg(level)));
    } else if (type.at(CMBtype->currentIndex()) == "optimize") {
        buff.append(QString("task %1 optimize\n").arg(QString("%1").arg(level)));
    } else if (type.at(CMBtype->currentIndex()) == "freq") {
        buff.append(QString("task %1 freq\n").arg(QString("%1").arg(level)));
    } else if (type.at(CMBtype->currentIndex()) == "opt freq") {
        buff.append(QString("task %1 optimize\n").arg(QString("%1").arg(level)));
        buff.append(QString("task %1 freq\n").arg(level));
    }

    if (ncen != kntcen){
        QMessageBox msgBox;
        msgBox.setText(tr("make_NWChem_input"));
        msgBox.setInformativeText(tr("Wrong number of centers in file:\n")+
            TXTextgeometry->text().trimmed());
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    extextEdit->setText(buff);

}

void Externals::make_Psi4_input(){
    extextEdit->clear();
    QStringList type = {"energy","optimize","freq","opt freq","nmr=giao"};
    QStringList level2 = {"","r","u","ro"};
    extOutputSuffix = "log";

    TXTextcommand->setText(extexecname[5]);
    QString filepath = TXTextworkdir->text().trimmed();
    if (filepath.isEmpty())
        return;
    QString filename = QFileInfo(TXTextgeometry->text()).baseName();
    if (filename.isEmpty())
        return;
    extInputFileName = filepath+"/"+filename+".psi4inp";
    TXTextcommand->setText(extexecname[5]+" "+extInputFileName);

    QFile geometryInput(TXTextgeometry->text().trimmed());
    if (!geometryInput.open(QFile::ReadOnly | QFile::Text)){
        return;
    }

    QByteArray buff;
    buff.append(QString("#! %1\n\n").arg(TXTtitle->text().trimmed()));
    if (SPBextproc->value() > 1){
        buff.append(QString("set_num_threads(%1)\n").arg(SPBextproc->value()));
    }
    if (!TXTextmem->text().isEmpty()){
        buff.append(QString("memory %1\n").arg(TXTextmem->text().trimmed()));
    }

    if (CMBlevel->currentIndex() < CMBlevel->count()-2){
        buff.append(QString("set basis %1\n").arg(CMBbasis->currentText()));
//        buff.append(QString("set basis %1\n").arg(CMBbasis->currentText()).arg(TXTkeywords->text().trimmed()));
    }
    if (level2.at(CMBlevel2->currentIndex()) != "") {
        buff.append(QString("set reference %1%2\n").arg(level2.at(CMBlevel2->currentIndex())).arg(CMBlevel->currentText().toLower()));
    }
    buff.append(QString("molecule {\n"));
    buff.append(QString("%1 %2\n").arg(SPBcharge->value()).arg(SPBmult->value()));

    QTextStream in(&geometryInput); // Buffer for reading from fileinput

    QString line = in.readLine();
#if QT_VERSION < 0x050E00
    QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
    QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
    int ncen = xyz.at(0).toInt();
    int kntcen = 0;
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xyz.count() == 4){
            buff.append(QString("%1    %2    %3    %4\n").arg(xyz[0]).arg(xyz[1]).arg(xyz[2]).arg(xyz[3]));
            kntcen++;
        }
    }
    buff.append(QString("}\n"));
    if (type.at(CMBtype->currentIndex()) == "energy") {
         buff.append(QString("energy,wfn=energy('%1',return_wfn=True)\n").arg(CMBlevel->currentText().toLower()));
    } else if (type.at(CMBtype->currentIndex()) == "optimize") {
         buff.append(QString("energy,wfn=optimize('%1',return_wfn=True)\n").arg(CMBlevel->currentText().toLower()));
    } else if (type.at(CMBtype->currentIndex()) == "freq") {
         buff.append(QString("frequencies('%1',ref_gradient=wfn.gradient())\n").arg(CMBlevel->currentText().toLower()));
    } else if (type.at(CMBtype->currentIndex()) == "opt freq") {
         buff.append(QString("energy,wfn=optimize('%1',return_wfn=True)\n").arg(CMBlevel->currentText().toLower()));
         buff.append(QString("frequencies('%1',ref_gradient=wfn.gradient())\n").arg(CMBlevel->currentText().toLower()));
    } 
    buff.append(QString("fchk_writer = psi4.core.FCHKWriter(wfn)\n"));
    buff.append(QString("fchk_writer.write('%1.fchk')\n").arg(filepath+"/"+filename));

    if (ncen != kntcen){
        QMessageBox msgBox;
        msgBox.setText(tr("make_Psi4_input"));
        msgBox.setInformativeText(tr("Wrong number of centers in file:\n")+
            TXTextgeometry->text().trimmed());
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    extextEdit->setText(buff);

}


void Externals::RBTlocal_changed(){
    if (RBTlocal->isChecked()){
        RBTPBS->setEnabled(false);
        RBTSGE->setEnabled(false);
        RBTSLURM->setEnabled(false);
        LBLextpathremote->setEnabled(false);
        TXTextpathremote->setEnabled(false);
        BTNjob->setEnabled(false);
        QDLexternal->setTabOrder(TXTkeywords,RBTlocal);
        QDLexternal->setTabOrder(RBTlocal,BTNjob);
    }
    else{
        RBTPBS->setEnabled(true);
        RBTSGE->setEnabled(true);
        RBTSLURM->setEnabled(true);
        LBLextpathremote->setEnabled(true);
        TXTextpathremote->setEnabled(true);
        BTNjob->setEnabled(true);
        QDLexternal->setTabOrder(TXTkeywords,RBTremote);
        QDLexternal->setTabOrder(RBTremote,BTNjob);
    }
    QDLexternal->adjustSize();
    QDLexternal->update();

}

void Externals::runformchk(){
    QString strprocess;
    QStringList Parameters;
    Parameters << QFileInfo(extInputFileName).absolutePath()
                  + "/" + QFileInfo(extInputFileName).completeBaseName() + ".chk";
    QString processname("formchk");
    QProcess *myProcess = new QProcess();
    myProcess->setStandardOutputFile(extOutputFileName,QIODevice::Append);
    myProcess->setStandardErrorFile(extOutputFileName,QIODevice::Append);
    connectionsext << connect(myProcess, SIGNAL(started()), this, SLOT(formchkStart()));
    connectionsext << connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(formchkError(QProcess::ProcessError)));
    connectionsext << connect(myProcess, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(formchkOutput(int,QProcess::ExitStatus)));
    myProcess->start(processname,Parameters);
}

void Externals::save_external_input(){
    if (extInputFileName.isEmpty())
        return;
    QFile extInputFile(extInputFileName);
    if (extInputFile.isOpen()){
        extInputFile.close();
    }
    if (!extInputFile.open(QFile::WriteOnly | QFile::Text)){
        return;
    }
    extInputFile.write(extextEdit->toPlainText().toLatin1());
#if QT_VERSION < 0x050E00
    extInputFile.write("\n");
#else
    extInputFile.write("\n");
#endif
}

void Externals::setTXTextgeometry(QString a){
    TXTextgeometry->setText(a.trimmed());
}

void Externals::setTXTextworkdir(QString a){
    TXTextworkdir->setText(a.trimmed());
}

void Externals::submitError(QProcess::ProcessError error){
    emit computing(QString(""));
    QString message = QString("Error %1 ").arg(error)
                + QString(tr("Process failed to start program %1\n").arg(TXTextcommand->text().trimmed())
                + QString(tr("Check that it is installed in your system\n"))
                + QString(tr("If it is installed, add it to your $PATH and try again)")));
    QMessageBox::critical(this,QString("Error %1").arg(error),message);
    BTNextsubmit->setEnabled(true);
}

void Externals::submitOutput(int exitCode, QProcess::ExitStatus status){
    if(status==QProcess::NormalExit){
        QMessageBox msgBox;
        msgBox.setText(tr("submitOutput"));
        msgBox.setInformativeText(QString(tr("Computation ended\n")));
        msgBox.setIcon(QMessageBox::Information);
        msgBox.exec();
        if (CMBengine->currentText().toLower() == "gaussian" && CHKformchk->isChecked()){  // Case of Gaussian: runs formchk if required
            runformchk();
        }
        else if (CMBengine->currentText().toLower() == "psi4"){
            QString strprocess;
            QStringList Parameters;
            Parameters << QFileInfo(TXTextgeometry->text()).absolutePath() + "/timer.dat"
                    << TXTextworkdir->text() + "/timer.dat";
            QString processname("mv");
            QProcess *myProcess = new QProcess();
            myProcess->start(processname,Parameters);
        }
        extexecname.replace(CMBengine->currentIndex(),TXTextcommand->text().split(" ").at(0));
        emit computing(QString(""));
        emit updatetextedit(extOutputFileName);
    }else if (status==QProcess::CrashExit){
        QMessageBox msgBox;
        msgBox.setText(tr("submitOutput"));
        msgBox.setInformativeText(QString(tr("Process crashed, exit code = %1").arg(exitCode).toLatin1()));
        msgBox.setIcon(QMessageBox::Information);
        msgBox.exec();
        emit computing(QString(""));
    }
    BTNextsubmit->setEnabled(true);
}

void Externals::submitStart(){
    BTNextsubmit->setEnabled(false);
    emit computing(QString(tr("Computing...")));
}

void Externals::TXTextgeometry_changed(){
    extgeomfile = QFileInfo(TXTextgeometry->text().trimmed()).fileName();
    extgeompath = QFileInfo(TXTextgeometry->text().trimmed()).path();
    if (TXTextworkdir->text().isEmpty()){
        TXTextworkdir->setText(extgeompath);
    }
    externalinputfile_changed();
}

