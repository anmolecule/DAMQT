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
//    Class defining the main window. Manages the graphic interface,
//    executes the auxiliary programs and displays the results
//
//  File:   mainwindow.cpp
//
//      Last version: October 2018
//
#include <string.h>
#include <QDesktopServices>
#include <QUrl>
#include <QtGlobal>
#include <QApplication>
#include <QPrinter>
#include <QPrintDialog>
#include <QSignalMapper>

#include <QtDebug>
#include <QLineEdit>
#include <QBoxLayout>
#include <QtCore/qprocess.h>

#include "mainwindow.h"
#include "GlobalInfo.h"

#define ANGSTROMTOBOHR 1.88971616463


/* Sets initial values */
MainWindow::MainWindow(QWidget *parent)
{
    QPixmap pixmap;
    pixmap.load(":/images/splash_3.2_en.png");

    initpointers();
    activebeware = false;
    ldZJjacobi = false;
    lzdo = false;
    lvalence = false;
    plotsknt = 1;
    topindex = -1;
    widgetsknt = 1;
    maxnumprocessors = MAX_NUM_PROCESSORS;
    BTNraiseplotslist.clear();
    BTNshowplotsslist.clear();
    BTNraisewidgetslist.clear();
    BTNshowwidgetslist.clear();
    connections2D.clear();
    connections3D.clear();

    splash = new QSplashScreen(pixmap);
    splash->show();

    FRMlanguage = new QDialog();
    FRMlanguage->setWhatsThis(tr("Check a language and push Start to start DAMQT."));
    FRMlanguage->setWindowTitle(tr("Language"));
    FRMlanguage->setFont(QFont("Helvetica",15));
    FRMlanguage->setMinimumWidth(480);
    FRMlanguage->setAttribute(Qt::WA_DeleteOnClose);

    createLanguageMenu();

    LBLlanguage = new QLabel();
    LBLlanguage->setText(tr("Choose language and push Start"));
    BTNlangstart=new QPushButton(QIcon(":/images/empezar.png"), tr("Start"));
    BTNlangstart->setMinimumWidth(120);
    connect(BTNlangstart, SIGNAL(clicked()), this, SLOT(start()));

    QHBoxLayout *languageLBLLayout = new QHBoxLayout();
    languageLBLLayout->addWidget(LBLlanguage);
    languageLBLLayout->setAlignment(Qt::AlignCenter);

    QHBoxLayout *languageHLayout=new QHBoxLayout();
    languageHLayout->addWidget(languageMenu);
    languageHLayout->setAlignment(Qt::AlignCenter);

    QHBoxLayout *languageBTNLayout=new QHBoxLayout();
    languageBTNLayout->addWidget(BTNlangstart);
    languageBTNLayout->setAlignment(Qt::AlignRight);

    QVBoxLayout *languageVLayout = new QVBoxLayout(FRMlanguage);
    languageVLayout->addLayout(languageLBLLayout);
    languageVLayout->addSpacing(20);
    languageVLayout->addLayout(languageHLayout);
    languageVLayout->addSpacing(10);
    languageVLayout->addLayout(languageBTNLayout);
    languageVLayout->addSpacing(10);

    FRMlanguage->exec();

    QString path=QApplication::applicationDirPath();

    #if defined(Q_WS_WIN) || defined(Q_OS_WIN)
        System="windows";
        LanguagePath="";
        iswindows = true;    // To be used by fortran programs
        mpi = false;
    #else //Q_WS_X11, Q_WS_MAC
        System="linux";
        LanguagePath=path;
        iswindows = false;    // To be used by fortran programs
        QProcess process;
        process.start("sh");    
        process.write("command -v mpirun");
        process.closeWriteChannel();
        process.waitForFinished(-1); // will wait forever until finished
        QByteArray outprocess = process.readAll();
        process.close();
        mpi = false;
        if (outprocess.count() != 0){
            mpi = true;
            mpicommand = QString("mpirun");
        }
        else{
            process.start("sh");    
            process.write("command -v mpiexec");
            process.closeWriteChannel();
            process.waitForFinished(-1); // will wait forever until finished
            outprocess = process.readAll();                   
            if (outprocess.count() != 0){
                    mpi = true;
                    mpicommand = QString("mpiexec");
            }
            process.close();
        }
    #endif
    setMinimumSize(900,600);
    setWindowState(Qt::WindowMaximized);
    setWindowIcon(QIcon(":/images/icon.png"));
    setWindowTitle(tr("DAMQT"));

    TAWprincipal = new QTabWidget(this);
    textEdit = new QTextEdit();
    textEdit->setReadOnly(true);

    int tabIndex=TAWprincipal->addTab(textEdit,QIcon(":/images/document_text.png"),tr("Results"));
    TAWprincipal->setCurrentIndex(tabIndex);
    connect(TAWprincipal, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));

    setCentralWidget(TAWprincipal);    

    CreateActions();
    CreateMenus();
    CreateToolBars();
    CreateStatusBar();
    CreateLeftMenu();    // Creates the left menu
    CreateRightMenu();    // Creates the right menu

    readSettings();

    SetCurrentFile("",true,false);
    for(int i = 0; i < NFORTRANPROCS; ++i){
        BTNstop[i]->setEnabled(false);
    }
    executing = -1;
    mden = 0;
}

void MainWindow::finishsplash(){
    splash->finish(this);
}

/* Defines action for TAB change */
void MainWindow::tabChanged(int index)
{
    if(index==0){
        AccPrint->setEnabled(true);
        AccPdf->setEnabled(true);
    }else{
        AccPrint->setEnabled(false);
        AccPdf->setEnabled(false);
    }
}

/* Ends all events */
void MainWindow::closeEvent(QCloseEvent *event)
{
    if(QDLwidget3D && !end_Viewer3DDialog()){
       event->ignore();
       return;
    }
    if(QDLviewer2D && !end_Viewer2DDialog()){
        event->ignore();
        return;
    }
    if (mustSave()) {
        writeSettings();
        event->accept();
    } 
    else {
        event->ignore();
    }
    delete QDLviewer2D;
    QDLviewer2D = nullpointer;
    delete QDLwidget3D;
    QDLwidget3D = nullpointer;
}

/* Changes a string to a qstring */
QString MainWindow::toQString(string v)
{
    QString qv=QString(v.c_str());
    return qv;
}

/* Changes a qstring to a string */
string MainWindow::toString(QString qv)
{
    QByteArray ba = qv.toLatin1();
    //QByteArray ba = qv.toUtf8();
    const char *v=ba.data();
    return v;
}

/*****************************************************************************************************/
/******************************* ACTIONS  ************************************************************/
/*****************************************************************************************************/
/* Creates program actions (Language, New, Open, Save, Save as, Print, Print PDF, 2D Viewer, 3D Viewer, Exit, Help, About) */
void MainWindow::CreateActions()
{
//    New
    AccNew = new QAction(QIcon(":/images/Nuevo.png"),tr("&New project"), this);
    AccNew->setShortcut(tr("Ctrl+N"));
    AccNew->setStatusTip(tr("Opens a new project"));
    connect(AccNew, SIGNAL(triggered()), this, SLOT(newProject()));
//    Open
    AccOpen = new QAction(QIcon(":/images/Abrir.png"), tr("&Open project..."), this);
    AccOpen->setShortcut(tr("Ctrl+A"));
    AccOpen->setStatusTip(tr("Open project file"));
    connect(AccOpen, SIGNAL(triggered()), this, SLOT(openProject()));
//    Save
    AccSave = new QAction(QIcon(":/images/Guardar.png"),tr("&Save project"), this);
    AccSave->setShortcut(tr("Ctrl+S"));
    AccSave->setStatusTip(tr("Save project file"));
    connect(AccSave, SIGNAL(triggered()), this, SLOT(saveProject()));
//    Save as
    AccSaveAs = new QAction(tr("Save project &as..."), this);
    AccSaveAs->setStatusTip(tr("Saves project file as"));
    connect(AccSaveAs, SIGNAL(triggered()), this, SLOT(SaveProjectAs()));
//    Print
    AccPrint = new QAction(QIcon(":/images/printer.png"),tr("&Print"), this);
    AccPrint->setShortcut(tr("Ctrl+P"));
    AccPrint->setStatusTip(tr("Print output file"));
    connect(AccPrint, SIGNAL(triggered()), this, SLOT(PrintFile()));
//    Print to PDF file
    AccPdf = new QAction(QIcon(":/images/acrobat.png"), tr("&Create Pdf"), this);
    AccPdf->setShortcut(tr("Ctrl+D"));
    AccPdf->setStatusTip(tr("Print output file as Pdf"));
    connect(AccPdf, SIGNAL(triggered()), this, SLOT(PrintFilePdf()));
//    2D Viewer2D
    Acc2Dplot = new QAction(QIcon(":/images/plot2D_tiny.png"),tr("&2D Viewer"), this);
    Acc2Dplot->setStatusTip(tr("2D Viewer"));
    connect(Acc2Dplot, SIGNAL(triggered()), this, SLOT(addviewer2D()));
//    3D Viewer
    Acc3Dview = new QAction(QIcon(":/images/cube_molecule.png"),tr("&3D Viewer"), this);
    Acc3Dview->setStatusTip(tr("3D Viewer"));
    connect(Acc3Dview, SIGNAL(triggered()), this, SLOT(addglWidget()));
//    Recent files
    for (int i = 0; i < MAX_ARCHIVOS_RECIENTES; ++i) {
        if (!AccRecentFiles[i]){
            AccRecentFiles[i] = new QAction(this);
        }
        AccRecentFiles[i]->setVisible(false);
        connect(AccRecentFiles[i], SIGNAL(triggered()), this, SLOT(openRecentProjects()));
    }
//    Exit
    AccExit = new QAction(QIcon(":/images/Salir.png"),tr("&Exit"), this);
    AccExit->setShortcut(tr("Ctrl+Q"));
    AccExit->setStatusTip(tr("Quit"));
    connect(AccExit, SIGNAL(triggered()), this, SLOT(close()));
//    Help
    AccHelp = new QAction(QIcon(":/images/ayuda.png"),tr("&Help"), this);
    AccHelp->setStatusTip(tr("Program help"));
    connect(AccHelp, SIGNAL(triggered()), this, SLOT(Help()));
//    About
    AccAbout = new QAction(QIcon(":/images/icon.png"),tr("&About DAMQT"), this);
    AccAbout->setStatusTip(tr("About DAMQT"));
    connect(AccAbout, SIGNAL(triggered()), this, SLOT(about()));
//    About Qt
    AccAboutQt = new QAction(QIcon(":/images/qtlogo.png"),tr("About &Qt"), this);
    AccAboutQt->setStatusTip(tr("About QT Library"));
    connect(AccAboutQt, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

/* Action: New file */
void MainWindow::newProject()
{
    if (mustSave()) {
        disable_pages();
        textEdit->clear();
        loadDefault(0);
        SetCurrentFile("",true,false);
    }
    lzdo = false;
    lvalence = false;
    emit initproject();
}

/* Action: Open existing project *.damproj file */
void MainWindow::openProject()
{
    if (mustSave()) { 
        disable_pages();
        QFileDialog filedialog(this);
        filedialog.setDirectory(ProjectFolder);
        filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
        QString fileName = filedialog.getOpenFileName(this,tr("Open input file"),ProjectFolder,
                tr("Project files")+" *.damproj (*.damproj);;"+tr("All files")+" (*)");
        if (!fileName.isEmpty()){
            bool res=Open(fileName);
            if (res){
                defineRanges();
                QMessageBox::information(this,tr("DAMQT"),tr("Project %1 open").arg(ProjectName), QMessageBox::Ok, 0); 
                page_atdens->setEnabled(true);
                if (QFile::exists(QString(fileName+"_2016.damqt"))){
                    page_densgrad->setEnabled(true);
                    page_Efield->setEnabled(true);
                    page_frad->setEnabled(true);
                    page_HFforces->setEnabled(true);
                    page_MED->setEnabled(true);
                    page_MESP->setEnabled(true);
                    page_orimult->setEnabled(true);
                    page_SGhole->setEnabled(true);
                    page_TOPO->setEnabled(true);
                    page_ZJdens->setEnabled(true);
                    page_ZJtab->setEnabled(true);
                }
                page_MO->setEnabled(true);
            }
            else
                return;
        }
    }
    emit initproject();
}

/* Action: Open file ggbs or sgbs*/
bool MainWindow::Open(const QString &fileName)
{
    textEdit->clear();
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("Open"),tr("File %1 cannot be read").arg(fileName)
                    +QString(": \n%1").arg(file.errorString()));
            return false;
    }
// Checks whether _2016.damqt, .ggbs or .sgbs, .xyz and .den files exist
    QString path = Path(fileName);
    if (path.at(path.length()-1) != '/') path.append('/');
    QString damqtfilename = path + FileWithoutExt(fileName)+"_2016.damqt";
    QString ggbsfilename  = path + FileWithoutExt(fileName)+".ggbs";
    QString sgbsfilename  = path + FileWithoutExt(fileName)+".sgbs";
    QString sgbsgzfilename  = path + FileWithoutExt(fileName)+".sgbs.gz";
    QString sgbsdenfilename  = path + FileWithoutExt(fileName)+".sgbsden";
    QString sgbsdengzfilename  = path + FileWithoutExt(fileName)+".sgbsden.gz";
    QString xyzfilename   = path + FileWithoutExt(fileName)+".xyz";
    QString denfilename   = path + FileWithoutExt(fileName)+".den";
    QString dengzfilename   = path + FileWithoutExt(fileName)+".den.gz";
    QString densprsbinfilename   = path + FileWithoutExt(fileName)+".densprsbin";
    if (((QFile::exists(ggbsfilename) || QFile::exists(sgbsfilename) || QFile::exists(sgbsgzfilename)) &&
                    (QFile::exists(denfilename) || QFile::exists(dengzfilename) || QFile::exists(densprsbinfilename)))
            || QFile::exists(sgbsdenfilename) || QFile::exists(sgbsdengzfilename) ){
        QApplication::setOverrideCursor(Qt::WaitCursor);
        if (QFile::exists(sgbsfilename) || QFile::exists(sgbsgzfilename)
                || QFile::exists(sgbsdenfilename) || QFile::exists(sgbsdengzfilename) ){
            lslater = true;
            CHKpotexact->setVisible(false);
        }
        else{ 
            lslater = false;
            CHKpotexact->setVisible(true);
        }
        CHKpotexact->setChecked(false);
        TXTImport->setText(fileName);
        ImportFile = FileWithoutPath(TXTImport->text());
        ImportFolder = path;
        TXTProjectFolder->setText(path);
        TXTProjectName->setText(FileWithoutExt(fileName)); 
        if (lslater){
            QString sxyzfilename = path + FileWithoutExt(fileName)+".sxyz";
            if (!(QFile::exists(sxyzfilename))){
                execsgbs2sxyz(sxyzfilename);
            }
            set_natom(read_natom(sxyzfilename));
        }
        else{
            set_natom(read_natom(ggbsfilename));
        }
        if(Extension(fileName)!="sgbs" && Extension(fileName)!="ggbs" && Extension(FileWithoutExt(fileName))!="sgbs"
                && Extension(fileName)!="sgbsden" && Extension(FileWithoutExt(fileName))!="sgbsden"){
            loadDefault(0);
            readOptions(fileName);
        }
        QApplication::restoreOverrideCursor();

        page_atdens->setEnabled(true);
        if (QFile::exists(damqtfilename)){
            page_densgrad->setEnabled(true);
            page_Efield->setEnabled(true);
            page_frad->setEnabled(true);
            page_HFforces->setEnabled(true);
            page_MED->setEnabled(true);
            page_MESP->setEnabled(true);
            page_orimult->setEnabled(true);
            page_SGhole->setEnabled(true);
            page_TOPO->setEnabled(true);
            page_ZJdens->setEnabled(true);
            page_ZJtab->setEnabled(true);
        }else{
            page_densgrad->setEnabled(false);
            page_Efield->setEnabled(false);
            page_frad->setEnabled(false);
            page_HFforces->setEnabled(false);
            page_MED->setEnabled(false);
            page_MESP->setEnabled(false);
            page_orimult->setEnabled(false);
            page_SGhole->setEnabled(false);
            page_TOPO->setEnabled(false);
            page_ZJdens->setEnabled(false);
            page_ZJtab->setEnabled(false);
        }
        QStringList filenames(QString(FileWithoutExt(fileName) +".GAorb*"));
        filenames << QString(FileWithoutExt(fileName) +".SLorb*");
        filenames << QString(FileWithoutExt(fileName) +".orb*");
        QStringList files;
        files = QDir(ProjectFolder).entryList(QStringList(filenames),QDir::Files);
        if (!files.isEmpty()) 
            page_MO->setEnabled(true);
        else 
            page_MO->setEnabled(false);
        SetCurrentFile(fileName,true,false);
        statusBar()->showMessage(tr("File succesfully loaded"), 2000);
        emit initproject();
        QMessageBox::information(this,tr("DAMQT"),tr("File succesfully loaded"), QMessageBox::Ok, 0);
        defineRanges();
        return true;
    }
    else{
        QMessageBox::warning(this,tr("DAMQT"),tr("Files %1 and/or %2 not found").arg(FileWithoutPath(ggbsfilename))
                .arg(FileWithoutPath(denfilename)));
        page_atdens->setEnabled(false);
        page_densgrad->setEnabled(false);
        page_Efield->setEnabled(false);
        page_frad->setEnabled(false);
        page_HFforces->setEnabled(false);
        page_MED->setEnabled(false);
        page_MESP->setEnabled(false);
        page_MO->setEnabled(false);
        page_orimult->setEnabled(false);
        page_SGhole->setEnabled(false);
        page_TOPO->setEnabled(false);
        page_ZJdens->setEnabled(false);
        page_ZJtab->setEnabled(false);
        statusBar()->showMessage(tr("Error loading file"), 2000);
        return false;
    }

}

/* Action: Open a recent file */
void MainWindow::openRecentProjects()
{
    QAction *action=qobject_cast<QAction *>(sender());
    if (action){
        QString filezdo = ProjectFolder + "zdo";
        if (QFileInfo(filezdo).exists()){
            lzdo = true;
        }
        else{
            lzdo = false;
        }
        QString filevalence = ProjectFolder + "valence";
        if (QFileInfo(filevalence).exists()){
            lvalence = true;
        }
        else{
            lvalence = false;
        }
        bool res=Open(action->data().toString());
        QStringList filenames(QString(ProjectName +".GAorb*"));
        filenames << QString(ProjectName +".SLorb*");
        filenames << QString(ProjectName +".orb*");
        QStringList files;
        files = QDir(ProjectFolder).entryList(QStringList(filenames),QDir::Files);
        if ((QFile::exists(QString(ProjectFolder+ProjectName+".ggbs")) ||
             QFile::exists(QString(ProjectFolder+ProjectName+".ggbs.gz")) ||
             QFile::exists(QString(ProjectFolder+ProjectName+".sgbs")) ||
                QFile::exists(QString(ProjectFolder+ProjectName+".sgbs.gz")))
                && (QFile::exists(QString(ProjectFolder+ProjectName+".den")) ||
                    QFile::exists(QString(ProjectFolder+ProjectName+".den.gz") ))){
            page_atdens->setEnabled(true);
            if (res){
                QMessageBox::information(this,tr("DAMQT"),tr("Project %1 open").arg(ProjectName), QMessageBox::Ok, 0);
                if (QFile::exists(QString(ProjectFolder+ProjectName+"_2016.damqt"))){
                    page_densgrad->setEnabled(true);
                    page_Efield->setEnabled(true);
                    page_frad->setEnabled(true);
                    page_HFforces->setEnabled(true);
                    page_MED->setEnabled(true);
                    page_MESP->setEnabled(true);
                    page_orimult->setEnabled(true);
                    page_SGhole->setEnabled(true);
                    page_TOPO->setEnabled(true);
                    page_ZJdens->setEnabled(true);
                    page_ZJtab->setEnabled(true);
                    readOptions(QString(ProjectFolder+ProjectName+".damproj"));
                }
            }
            if (!files.isEmpty())
                page_MO->setEnabled(true);
            else
                page_MO->setEnabled(false);
        }
        else{
            QMessageBox::information(this,tr("DAMQT"),tr("Project %1 cannot be opened").arg(ProjectName), QMessageBox::Ok, 0);
            disable_pages();
        }
    }
    else
        return;
}

/* Action: Update recent files */
void MainWindow::UpdateRecentFiles()
{
    QMutableStringListIterator i(ArchivosRecientes);
    while (i.hasNext()) {
        if (!QFile::exists(i.next())){
            i.remove();
        }
    }
    for (int j = 0; j < MAX_ARCHIVOS_RECIENTES; ++j) {
        if (j < ArchivosRecientes.count()) {
            QString text = tr("&%1 %2").arg(j+1).arg(FileWithoutPath(ArchivosRecientes[j]));
            AccRecentFiles[j]->setText(text);
            AccRecentFiles[j]->setData(ArchivosRecientes[j]);
            AccRecentFiles[j]->setVisible(true);
        } 
        else {
            AccRecentFiles[j]->setVisible(false);
        }
    }
}

/* Action: Save file */
bool MainWindow::saveProject()
{    
    if (ProjectName.isEmpty()) {
        return SaveProjectAs();
    } 
    else {
        return Save(ArchivoActual);
    }
}

/* Action: Save file as */
bool MainWindow::SaveProjectAs()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getSaveFileName(this);
    if (fileName.isEmpty()){
        return false;
    }
    if (Extension(fileName) != "options"){
        fileName.append(".damproj");
    }
    return Save(fileName);
}

/* Action: Save project */
bool MainWindow::Save(const QString &fullfileName)
{
    QString filePath = Path(fullfileName)+"/";
    QString fileName = ProjectName+".damproj";
    QDir path(filePath);
    if (!path.exists()) {
        filePath = ProjectFolder;
        path = ProjectFolder;
    }
    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
            return false;
        }
    }

    QFile file(filePath+fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("File succesfully saved")
            .arg(fileName)+QString(":\n%2.").arg(file.errorString()));
        return false;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    saveOptions(filePath+fileName,0);
    QApplication::restoreOverrideCursor();
    statusBar()->showMessage(tr("File succesfully saved"), 2000);
    SetCurrentFile(filePath+fileName,true,false);
    return true;
}

/* Returns whether a file must be saved or not */
bool MainWindow::mustSave()
{
    if (changes) {
        QMessageBox msgBox;
        msgBox.setInformativeText(tr("Document has been modified")+"\n"+tr("Do you want to save changes?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes) return saveProject();
        else if (ret == QMessageBox::Cancel) return false;
    }
    return true;
}

/* Action: About */
void MainWindow::about()
{
    QMessageBox::about(this,tr("About DAMQT"), "<h2>"+tr("DAMQT 3.2")+"</h2>" "<p>"
        +tr("Copyright &copy; 2008-2019")+"</p>"
        "<p>"+tr("DAMQT is a program for the analysis of electron molecular density, "
        "electrostatic potential and field and Hellman-Feynman forces on nuclei.")
        +"</p>" "<p>"+tr("Developed in the Departamento de Quimica-Fisica Aplicada of "
        "the Universidad Autonoma de Madrid (Spain) in collaboration with the Departamento de "
        "Quimica-Fisica of the Universidad de Cadiz (Spain) and with the "
        "Indian Institute of Technology Kanpur (India).")+"</p>");
}

/* Action: Print */
void MainWindow::PrintFile()
{
    if (TAWprincipal->currentIndex()==0){
        if (!textEdit->document()->isEmpty()){
            QTextDocument *document = textEdit->document();
            QPrinter printer(QPrinter::HighResolution);
            QPrintDialog *dialog = new QPrintDialog(&printer,this);
            dialog->setWindowTitle(tr("Printing options"));
            if (dialog->exec() != QDialog::Accepted)return;
            printer.setFullPage(true);
            printer.setPageSize(QPrinter::A4);
            document->print(&printer);
        }
    }
}

/* Action: Print PDF */
void MainWindow::PrintFilePdf()
{
    if (TAWprincipal->currentIndex()==0){
        if (!textEdit->document()->isEmpty()){
            QTextDocument *document = textEdit->document();
            QPrinter printer(QPrinter::HighResolution);
            printer.setOutputFormat(QPrinter::PdfFormat);
            QFileDialog filedialog(this);
            filedialog.setDirectory(ProjectFolder);
            filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
            QString fileName = filedialog.getSaveFileName(this,tr("Save pdf..."),ProjectFolder,
                    tr("pdf files")+" (*.pdf);;"+tr("All files")+" (*)");
            if (fileName.isEmpty())return;
            if (QFile::exists(fileName)) {
                QMessageBox msgBox;
                msgBox.setInformativeText(tr("File %1 exists").arg(QDir::toNativeSeparators(fileName))
                        +"\n"+tr("Do you want to overwrite?"));
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
                msgBox.setDefaultButton(QMessageBox::Cancel);
                msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
                msgBox.setButtonText(QMessageBox::No, tr("No"));
                msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
                msgBox.setIcon(QMessageBox::Warning);
                int ret = msgBox.exec();
                if (!(ret==QMessageBox::Yes)) return;
            }
            printer.setOutputFileName(fileName);
            printer.setFullPage(true);
            printer.setPaperSize(QPrinter::A4);
            document->print(&printer);
        }
    }
}

/* Action: Help */
void MainWindow::Help()
{
    QString path;
    if (QApplication::applicationDirPath() == "/usr/local/bin")
        path=QApplication::applicationDirPath()+"/../doc/DAMQT_3.2.0_manual.pdf";
    else{
        path=QApplication::applicationDirPath()+"/DAMQT_3.2.0_manual.pdf";
        if (!QFileInfo::exists(path))
            path = QCoreApplication::applicationDirPath()+"/../../doc/manual"+"/DAMQT_3.2.0_manual.pdf";
    }
    bool r = QDesktopServices::openUrl(QUrl::fromLocalFile(path));
    if (!r)
        QMessageBox::warning(this,tr("DAMQT"),tr("Help file does not exist in ")+path,
                QMessageBox::Ok);
}

/* Action: viewer2D */
void MainWindow::menu_viewer2D()
{
    for (int i = 0 ; i < connections2D.size() ; i++){
        QObject::disconnect(connections2D.at(i));
    }
    connections2D.clear();
    plots = new QList<Viewer2D*>();
    QDLviewer2D = new ViewerDialog();
    QDLviewer2D->setMinimumSize(250,80);
    BTNnewplot = new QPushButton(tr("New 2D Plotter"));
    BTNnewplot->setToolTip(tr("Creates a new window for 2D plotting"));
    connections2D << connect(BTNnewplot, SIGNAL(clicked()), this, SLOT(addviewer2D()));

    QLabel *label2D = new QLabel(tr("2D plotters"));
    label2D->setStyleSheet("QLabel { color : blue; }");

    QVBoxLayout *layout2=new QVBoxLayout(QDLviewer2D);
    layout2->addWidget(label2D);
    layout2->addWidget(BTNnewplot);
    layout2->addStretch();
}

/* Action: viewer3D */
void MainWindow::menu_viewer3D()
{
    for (int i = 0 ; i < connections3D.size() ; i++){
        QObject::disconnect(connections3D.at(i));
    }
    connections3D.clear();
    widgets = new QList<glWidget*>();
    QDLwidget3D = new ViewerDialog();
    QDLwidget3D->setMinimumSize(250,80);
    BTNnewwidget = new QPushButton(tr("New 3D Window"));
    BTNnewwidget->setToolTip(tr("Creates a new window for 3D Viewer"));
    connections3D << connect(BTNnewwidget, SIGNAL(clicked()), this, SLOT(addglWidget()));

    QLabel *label3D = new QLabel(tr("3D viewers"));
    label3D->setStyleSheet("QLabel { color : red; }");

    QVBoxLayout *layout2 = new QVBoxLayout(QDLwidget3D);
    layout2->addWidget(label3D);
    layout2->addWidget(BTNnewwidget);
    layout2->addStretch();

}

/***************************************************************************/
/******************************* MENUS *************************************/
/***************************************************************************/
/* Creates the menus of the graphic interface*/
void MainWindow::CreateMenus()
{
    FileMenu = menuBar()->addMenu(tr("&File"));
    FileMenu->addAction(AccNew);
    FileMenu->addAction(AccOpen);
    FileMenu->addAction(AccSave);
    FileMenu->addAction(AccSaveAs);
    FileMenu->addSeparator();
    FileMenu->addAction(AccPrint);
    FileMenu->addAction(AccPdf);
    if(MAX_ARCHIVOS_RECIENTES>0){
        FileMenu->addSeparator();
        for (int i = 0; i < MAX_ARCHIVOS_RECIENTES; ++i){
            FileMenu->addAction(AccRecentFiles[i]);
        }
    }
    FileMenu->addSeparator();
    FileMenu->addAction(AccExit);
    UpdateRecentFiles();
    menuBar()->addSeparator();
    GraphicsMenu = menuBar()->addMenu(tr("&Graphics"));
    GraphicsMenu->addAction(Acc2Dplot);
    GraphicsMenu->addAction(Acc3Dview);
    menuBar()->addSeparator();
    HelpMenu = menuBar()->addMenu(tr("&Help"));
    HelpMenu->addAction(AccHelp);
    HelpMenu->addSeparator();
    HelpMenu->addAction(AccAbout);
    HelpMenu->addAction(AccAboutQt);
}

/***************************************************************************/
/******************************* TOOLBARS **********************************/
/***************************************************************************/
/* Toolbars */
void MainWindow::CreateToolBars()
{
    ToolBarFile = addToolBar(tr("Project folder"));
    ToolBarFile->addAction(AccNew);
    ToolBarFile->addAction(AccOpen);
    ToolBarFile->addAction(AccSave);
    ToolBarFile->addAction(AccPrint);
    ToolBarFile->addAction(AccPdf);
    ToolBarHelp=addToolBar(tr("Graphics"));
    ToolBarHelp->addAction(Acc2Dplot);
    ToolBarHelp->addAction(Acc3Dview);
    ToolBarHelp=addToolBar(tr("Help"));
    ToolBarHelp->addAction(AccHelp);
    ToolBarHelp->addAction(AccAbout);
    ToolBarHelp->addAction(AccExit);
}

/* Statusbars*/
void MainWindow::CreateStatusBar()
{
    statusBar()->showMessage(tr("DAMQT"));
}

/*******************************************************************************************************/
/******************************** OTHER ****************************************************************/
/*******************************************************************************************************/
/* Reads initial settings */
void MainWindow::readSettings()
{
    QSettings settings("DAMQT", "Densidades");
    QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
    QSize size = settings.value("size", QSize(400, 400)).toSize();
    resize(size);
    move(pos);
    ArchivosRecientes = settings.value("recentFiles").toStringList();
    UpdateRecentFiles();
}

/* Writes settings for next session */
void MainWindow::writeSettings()
{
    QSettings settings("DAMQT", "Densidades");
    settings.setValue("pos", pos());
    settings.setValue("size", size());
    QStringList NewArchRecent;
    for (int i = 0 ; i < min(2*MAX_ARCHIVOS_RECIENTES,ArchivosRecientes.size()) ; i++ ){
        NewArchRecent << ArchivosRecientes.at(i);
    }
    settings.setValue("recentFiles", NewArchRecent);
}

/* Determines the file currently open */
void MainWindow::SetCurrentFile(const QString &fileName,bool usar,bool modificado)
{
    changes = modificado; //indicates whether a change has been made or not 
    if (usar){ // indicates whether variable filename is used
        ArchivoActual = fileName;
    }
    QString MostrarNombre = tr("Unnamed");
    if(!ArchivoActual.isEmpty()){
        MostrarNombre=FileWithoutPath(ArchivoActual);
        ArchivosRecientes.removeAll(ArchivoActual);
        ArchivosRecientes.prepend(ArchivoActual);
        UpdateRecentFiles();
    } 
    setWindowTitle(tr("%1 - %2 [*]").arg(tr("DAMQT")).arg(MostrarNombre));
    if (!changes){
        setWindowModified(false);
    }
    else{
        setWindowModified(true);
    }
}


/***************************************************************************/
/*****************************  DOCK  WINDOWS ******************************/
/***************************************************************************/

/* Creates DockWindows */

//      CreateLeftMenu
//
void MainWindow::CreateLeftMenu()
{
    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);
    QDockWidget *dock = new QDockWidget(tr("Options"),this);
    dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    dock->setMaximumSize(QSize(500, 2000));

    toolBox = new mainmenu(dock);
    toolBox->resize(QSize(400, 2000));
    toolBox->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

//    page_project: PROJECT

    page_project = new QWidget();
    page_project->setEnabled(true);
    page_project_widgets();
    page_project_layouts();

    toolBox->addItem(page_project,QIcon(":/images/icon.png"), tr("Project"));

//    page_atdens: ATOMIC DENSITIES

    page_atdens = new QWidget();
    page_atdens->setEnabled(false);
    page_atdens_widgets();
    page_atdens_layouts();

    toolBox->addItem(page_atdens,QIcon(":/images/icon.png"),tr("Atomic densities"));

//    page_MED:    DENSITY

    page_MED = new QWidget();
    page_MED->setEnabled(false);
    page_MED_widgets();
    page_MED_layouts();

    toolBox->addItem(page_MED,QIcon(":/images/icon.png"),tr("Density"));

//    page_MESP: ELECTROSTATIC POTENTIAL

    page_MESP = new QWidget();
    page_MESP->setEnabled(false);
    page_MESP_widgets();
    page_MESP_layouts();

    toolBox->addItem(page_MESP,QIcon(":/images/icon.png"),tr("Electrostatic potential"));

//    page_MO: MOLECULAR ORBITALS

    page_MO = new QWidget();
    page_MO->setEnabled(false);
    page_MO_widgets();
    page_MO_layouts();

    toolBox->addItem(page_MO,QIcon(":/images/icon.png"),tr("Molecular orbitals"));

//    page_TOPO: MOLECULAR TOPOGRAPHY

    page_TOPO = new QWidget();
    page_TOPO->setEnabled(false);
    page_TOPO_widgets();
    page_TOPO_layouts();

    toolBox->addItem(page_TOPO,QIcon(":/images/icon.png"),tr("Molecular topography"));

//    page_SGhole: MESP sigma hole

    page_SGhole = new QWidget();
    page_SGhole->setEnabled(false);
    page_SGhole_widgets();
    page_SGhole_layouts();

    toolBox->addItem(page_SGhole,QIcon(":/images/icon.png"),tr("Mesp sigma hole"));

//    page_Efield: ELECTRIC FIELD

    page_Efield = new QWidget();
    page_Efield->setEnabled(false);
    page_Efield_widgets();
    page_Efield_layouts();

    toolBox->addItem(page_Efield,QIcon(":/images/icon.png"),tr("Electric field"));

//    page_densgrad: DENSITY GRADIENT

    page_densgrad = new QWidget();
    page_densgrad->setEnabled(false);
    page_densgrad_widgets();
    page_densgrad_layouts();

    toolBox->addItem(page_densgrad,QIcon(":/images/icon.png"),tr("Density gradient"));

//    page_HFforces: HELLMANN-FEYNMAN FORCES ON NUCLEI

    page_HFforces = new QWidget();
    page_HFforces->setEnabled(false);
    page_HFforces_widgets();
    page_HFforces_layouts();

    toolBox->addItem(page_HFforces,QIcon(":/images/icon.png"),tr("Hellmann-Feynman forces on nuclei"));

//    page_frad: RADIAL FACTORS

    page_frad = new QWidget();
    page_frad->setEnabled(false);
    page_frad_widgets();
    page_frad_layouts();

    toolBox->addItem(page_frad,QIcon(":/images/icon.png"),tr("Radial factors"));

//    page_orimult: ORIENTED MULTIPOLES

    page_orimult = new QWidget();
    page_orimult->setEnabled(false);
    page_orimult_widgets();
    page_orimult_layouts();

    toolBox->addItem(page_orimult,QIcon(":/images/icon.png"),tr("Oriented multipoles"));

//    page_ZJdens: ZERNIKE-JACOBI EXPANSIONS OF DENSITY

    page_ZJdens = new QWidget();
    page_ZJdens->setEnabled(false);
    page_ZJdens_widgets();
    page_ZJdens_layouts();

    toolBox->addItem(page_ZJdens,QIcon(":/images/icon.png"),tr("Zernike-Jacobi density expansion"));
    
//    page_ZJtab: TABULATION OF DENSITY FROM ZERNIKE-JACOBI EXPANSIONS

    page_ZJtab = new QWidget();
    page_ZJtab->setEnabled(false);
    page_ZJtab_widgets();
    page_ZJtab_layouts();

    toolBox->addItem(page_ZJtab,QIcon(":/images/icon.png"),tr("Zernike-Jacobi density tabulation"));

    dock->setFeatures(QDockWidget::NoDockWidgetFeatures);
    QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Preferred);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(0);
    sizePolicy1.setHeightForWidth(dock->sizePolicy().hasHeightForWidth());
    dock->setSizePolicy(sizePolicy1);
    dock->setWidget(toolBox);
    addDockWidget(Qt::LeftDockWidgetArea,dock);

//    LOAD DEFAULT OPTIONS

    loadDefault(0);

    disable_pages();
}
    
//      End of CreateLeftMenu
//      =====================

//      CreateRightMenu
//

void MainWindow::CreateRightMenu()
{    
    if (BTNnewplot){
        delete BTNnewplot;
        BTNnewplot = nullpointer;
    }
    if (FRMplots){
        delete FRMplots;
        FRMplots = nullpointer;
    }
    if (QDLviewer2D){
        delete QDLviewer2D;
        QDLviewer2D = nullpointer;
    }

    if (BTNnewwidget){
        delete BTNnewwidget;
        BTNnewwidget = nullpointer;
    }
    if (FRMviewers){
        delete FRMviewers;
        FRMviewers = nullpointer;
    }
    if (QDLwidget3D){
        delete QDLwidget3D;
        QDLwidget3D = nullpointer;
    }
    if ( BTNraiseviewers){
        delete  BTNraiseviewers;
         BTNraiseviewers = nullpointer;
    }
    if (dockright){
        delete dockright;
        dockright = nullpointer;
    }
    dockright = new QDockWidget(tr(""),this);
    dockright->setAllowedAreas(Qt::RightDockWidgetArea);
    dockright->resize(QSize(500, this->height()));
    dockright->setFeatures(QDockWidget::DockWidgetMovable);
    dockright->setFeatures(QDockWidget::DockWidgetFloatable);
    QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(1000);
    sizePolicy1.setHeightForWidth(dockright->sizePolicy().hasHeightForWidth());
    dockright->setSizePolicy(sizePolicy1);
    addDockWidget(Qt::RightDockWidgetArea,dockright);

    menu_viewer2D();
    menu_viewer3D();

    BTNraiseviewers = new QPushButton(tr("Raise viewers"));
    BTNraiseviewers->setToolTip(tr("Moves viewers to front"));
    BTNraiseviewers->setStyleSheet("QPushButton {background-color: darkGreen; color: white;}");
    connect(BTNraiseviewers, SIGNAL(clicked()), this, SLOT(updatewindowsoverlay()));

    QVBoxLayout *layout1 = new QVBoxLayout();
    layout1->addWidget(BTNraiseviewers);

    QWidget *btnwidget = new QWidget();
    btnwidget->setLayout(layout1);

    dockwidget = new QWidget();
    QVBoxLayout *layout = new QVBoxLayout(dockwidget);
    layout->addWidget(btnwidget);
    layout->addWidget(QDLviewer2D);
    layout->addWidget(QDLwidget3D);
    layout->addStretch();

    dockright->setWidget(dockwidget);
    updatewindowsoverlay();
}

//      End of CreateRightMenu
//      =========================

/***************************************************************************/
/*                 PAGES WIDGETS AND LAYOUTS                               */
/***************************************************************************/

//    page_project: PROJECT
//    =====================

void MainWindow::page_project_widgets()
{
//            Create project

    FRMproject=new QGroupBox(tr("Create Project"),page_project);
    FRMproject->setMaximumSize(QSize(400, 2000));

//                    Import File

    LBLImport = new QLabel(tr("Import data from")+":");
    TXTImport = new QLineEdit();
    connect(TXTImport, SIGNAL(textChanged(const QString &)), this, SLOT(TXTImport_changed()));

    BTNImport = new QToolButton();
    BTNImport->setText(tr("..."));
    connect(BTNImport, SIGNAL(clicked()), this, SLOT(importFile()));

//                    Project folder

    LBLProjectFolder = new QLabel(tr("Project folder")+":");
    TXTProjectFolder = new QLineEdit();
    connect(TXTProjectFolder, SIGNAL(textChanged(const QString &)), this, SLOT(TXTProjectFolder_changed(const QString )));
    TXTProjectFolder->setEnabled(true);

//                    Project Name

    LBLProjectName=new QLabel(tr("Project name")+":");
    TXTProjectName=new QLineEdit();
    connect(TXTProjectName, SIGNAL(textChanged(const QString &)), this, SLOT(TXTProjectName_changed(const QString &)));
    TXTProjectName->setEnabled(true);

//                    Exec

    BTNexecImport=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"));
    BTNexecImport->setEnabled(false);
    connect(BTNexecImport, SIGNAL(clicked()), this, SLOT(execImport()));

//                    MPI
    if (mpi){
        FRMmpioptions=new QGroupBox(tr("MPI options"),page_project);
        FRMmpioptions->setMaximumSize(QSize(400, 400));
        LBLmpicommand = new QLabel(tr("MPI command")+":");
        LBLmpiflags = new QLabel(tr("MPI flags")+":");
        TXTmpicommand = new QLineEdit();
        TXTmpicommand->setText(mpicommand);
        connect(TXTmpicommand, SIGNAL(textChanged(const QString &)), this, SLOT(TXTmpicommand_changed()));
        TXTmpiflags = new QLineEdit();
        connect(TXTmpiflags, SIGNAL(textChanged(const QString &)), this, SLOT(TXTmpiflags_changed()));
    }
}

void MainWindow::page_project_layouts()
{
    QHBoxLayout *Layout1=new QHBoxLayout();
    Layout1->addWidget(LBLImport);

    QHBoxLayout *Layout2=new QHBoxLayout();
    Layout2->addWidget(TXTImport);
    Layout2->addWidget(BTNImport);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addWidget(LBLProjectFolder);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(TXTProjectFolder);

    QHBoxLayout *Layout5=new QHBoxLayout();
    Layout5->addWidget(LBLProjectName);

    QHBoxLayout *Layout6=new QHBoxLayout();
    Layout6->addWidget(TXTProjectName);

    QHBoxLayout *Layout7=new QHBoxLayout();
    Layout7->addWidget(BTNexecImport,0,Qt::AlignRight);

    QVBoxLayout *Layout8=new QVBoxLayout(FRMproject);
    Layout8->addLayout(Layout1);
    Layout8->addLayout(Layout2);
    Layout8->addLayout(Layout3);
    Layout8->addLayout(Layout4);
    Layout8->addLayout(Layout5);
    Layout8->addLayout(Layout6);
    Layout8->addLayout(Layout7);

    QVBoxLayout *page_projectLayout = new QVBoxLayout(page_project);
    page_projectLayout->addWidget(FRMproject);

    if (mpi){
        QHBoxLayout *Layout9=new QHBoxLayout();
        Layout9->addWidget(LBLmpicommand);
        QHBoxLayout *Layout10=new QHBoxLayout();
        Layout10->addWidget(TXTmpicommand);
        QHBoxLayout *Layout11=new QHBoxLayout();
        Layout11->addWidget(LBLmpiflags);
        QHBoxLayout *Layout12=new QHBoxLayout();
        Layout12->addWidget(TXTmpiflags);
        QVBoxLayout *Layout13=new QVBoxLayout(FRMmpioptions);
        Layout13->addLayout(Layout9);
        Layout13->addLayout(Layout10);
        Layout13->addLayout(Layout11);
        Layout13->addLayout(Layout12);
        page_projectLayout->addWidget(FRMmpioptions);
    }

    page_projectLayout->addStretch();
}

//    page_atdens: ATOMIC DENSITIES
//    =============================

void MainWindow::page_atdens_widgets()
{
//            Highest l in expansion

    FRMatdenslmaxexp = new QGroupBox(tr("Highest l in expansion"),page_atdens);
    FRMatdenslmaxexp->setMaximumSize(QSize(400, 2000));

    SPBatdenslmaxexp=new QSpinBox(FRMatdenslmaxexp);
    SPBatdenslmaxexp->setRange(0, 25);
    SPBatdenslmaxexp->setValue(10);
    SPBatdenslmaxexp->setMaximumWidth(50);
    SPBatdenslmaxexp->setToolTip(tr("Size of multipolar expansion"));
    connect(SPBatdenslmaxexp, SIGNAL(valueChanged(int)), this, SLOT(SPBatdenslmaxexp_changed()));

    FRMatdenslmaxdisp = new QGroupBox(tr("Highest l to be displayed"),page_atdens);
    FRMatdenslmaxdisp->setMaximumSize(QSize(400, 2000));

    SPBatdenslmaxdisp=new QSpinBox(FRMatdenslmaxdisp);
    SPBatdenslmaxdisp->setRange(0, SPBatdenslmaxexp->value());
    SPBatdenslmaxdisp->setValue(5);
    SPBatdenslmaxdisp->setMaximumWidth(50);
    SPBatdenslmaxdisp->setToolTip(tr("Highest order of multipoles to be displayed in output file and whose modules will be tabulated in file .mltmod"));

//            Type of fitting

    FRMatdenstype = new QGroupBox(tr("Type of fitting"),page_atdens);
    FRMatdenstype->setMaximumSize(QSize(400, 2000));

    RBTatdensDtotal=new QRadioButton(tr("Total density"),FRMatdenstype);
    RBTatdensDtotal->setChecked(true);
    RBTatdensD1center=new QRadioButton(tr("One-center terms"),FRMatdenstype);
    RBTatdensD2center=new QRadioButton(tr("Two-center terms"),FRMatdenstype);
    connect(RBTatdensDtotal, SIGNAL(toggled (bool)), this, SLOT(TXTValidate_changed()));
    connect(RBTatdensD1center, SIGNAL(toggled (bool)), this, SLOT(TXTValidate_changed()));
    connect(RBTatdensD2center, SIGNAL(toggled (bool)), this, SLOT(TXTValidate_changed()));

    FRMfitthrs = new QGroupBox(tr("Thresholds"),page_atdens);
    FRMfitthrs->setMaximumSize(QSize(400, 2000));


    LBLfitthreshold = new QLabel("Fitting threshold: 10^",FRMfitthrs);
    LBLfitthreshold->setToolTip(tr("Threshold for truncating expansions of radial factors "));
    LBLfradthreshold = new QLabel("Cutoff threshold: 10^",FRMfitthrs);
    LBLfradthreshold->setToolTip(tr("Threshold for neglecting radial factors "));

    SPBfitthreshold = new QSpinBox(FRMfitthrs);
    SPBfitthreshold->setRange(-20, 0);
    SPBfitthreshold->setValue(-14);
    SPBfitthreshold->setMaximumWidth(60);
    connect(SPBfitthreshold, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

    SPBfradthreshold = new QSpinBox(FRMfitthrs);
    SPBfradthreshold->setRange(-20, 0);
    SPBfradthreshold->setValue(-14);
    SPBfradthreshold->setMaximumWidth(60);
    connect(SPBfradthreshold, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

//            Generate input only

    FRMatdensinput = new QGroupBox(tr("Input only"),page_atdens);
    FRMatdensinput->setMaximumSize(QSize(400, 2000));

    CHKatdensinput=new QCheckBox(tr("Generate input file only"),FRMatdensinput);
    CHKatdensinput->setChecked(false);
    CHKatdensinput->setEnabled(true);
    connect(CHKatdensinput, SIGNAL(stateChanged(int)), this, SLOT(CHKatdensinput_changed(int)));

//            MPI

    FRMatdensmpi = new QGroupBox(tr("Parallel computing"),page_atdens);
    FRMatdensmpi->setMaximumSize(QSize(400, 2000));

    CHKatdensmpi = new QCheckBox(tr("MPI"),FRMatdensmpi);
    LBLatdensmpi = new QLabel(tr("Number of processors"),FRMatdensmpi);

    SPBatdensmpi = new QSpinBox(FRMatdensmpi);
    SPBatdensmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBatdensmpi->setValue(1);
    SPBatdensmpi->setMaximumWidth(50);
    SPBatdensmpi->setToolTip(tr("Number of processors"));

    if (mpi){
        FRMatdensmpi->setVisible(true);
        FRMatdensmpi->setEnabled(true);
        CHKatdensmpi->setChecked(true);
        SPBatdensmpi->setEnabled(true);
}
    else{
        FRMatdensmpi->setHidden(true);
        CHKatdensmpi->setChecked(false);
        SPBatdensmpi->setEnabled(false);
    }
    connect(CHKatdensmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKatdensmpi_changed(int)));
    connect(SPBatdensmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBatdensmpi_changed(int)));

//            List file, Stop, Execute

    BTNexecDam=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_atdens);
    BTNexecDam->setToolTip(tr("Electron density analysis"));
    connect(BTNexecDam, SIGNAL(clicked()), this, SLOT(execDam()));

    if (!BTNstop[0])
        BTNstop[0]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_atdens);
    BTNstop[0]->setToolTip(tr("Kill the process"));
    connect(BTNstop[0], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[0])
        BTNtexto[0]=new QPushButton(QIcon(":/images/document_text.png"),"",page_atdens);
    BTNtexto[0]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[0], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_atdens_layouts()
{
    QHBoxLayout *Layout1=new QHBoxLayout(FRMatdenslmaxexp);
    Layout1->addWidget(SPBatdenslmaxexp);

    QHBoxLayout *Layout2=new QHBoxLayout(FRMatdenslmaxdisp);
    Layout2->addWidget(SPBatdenslmaxdisp);

    QVBoxLayout *Layout3=new QVBoxLayout(FRMatdenstype);
    Layout3->addWidget(RBTatdensDtotal);
    Layout3->addWidget(RBTatdensD1center);
    Layout3->addWidget(RBTatdensD2center);

    QVBoxLayout *Layout4 = new QVBoxLayout(FRMatdensinput);
    Layout4->addWidget(CHKatdensinput);

    QGridLayout *Layout5 = new QGridLayout(FRMfitthrs);
    Layout5->addWidget(LBLfradthreshold,0,0,Qt::AlignRight);
    Layout5->addWidget(SPBfradthreshold,0,1,Qt::AlignLeft);
    Layout5->addWidget(LBLfitthreshold,1,0,Qt::AlignRight);
    Layout5->addWidget(SPBfitthreshold,1,1,Qt::AlignLeft);

    QHBoxLayout *Layout7 = new QHBoxLayout(FRMatdensmpi);
    Layout7->addWidget(CHKatdensmpi);
    Layout7->addWidget(LBLatdensmpi);
    Layout7->addWidget(SPBatdensmpi);

    QHBoxLayout *Layout8=new QHBoxLayout();
    Layout8->addWidget(BTNtexto[0],0,Qt::AlignLeft);
    Layout8->addWidget(BTNstop[0],0,Qt::AlignRight);
    Layout8->addWidget(BTNexecDam,0,Qt::AlignRight);

    QVBoxLayout *page_atdensLayout = new QVBoxLayout(page_atdens);
    page_atdensLayout->addWidget(FRMatdenslmaxexp);
    page_atdensLayout->addWidget(FRMatdenslmaxdisp);
    page_atdensLayout->addWidget(FRMatdenstype);
    page_atdensLayout->addWidget(FRMfitthrs);
    page_atdensLayout->addWidget(FRMatdensinput);
    page_atdensLayout->addWidget(FRMatdensmpi);
    page_atdensLayout->addLayout(Layout8);
    page_atdensLayout->addStretch();
}

//    page_MED:    DENSITY
//    ====================

void MainWindow::page_MED_widgets()
{
//            Output file prefix
    FRMdensdamdenfilename = new QGroupBox(tr("Output files prefix"),page_MED);
    FRMdensdamdenfilename->setMaximumSize(QSize(400, 2000));
    TXTdensdamdenfilename = new QLineEdit(FRMdensdamdenfilename);
    connect(TXTdensdamdenfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

//            Density settings

    FRMdensexpansion = new QGroupBox(tr("Density settings"),page_MED);
    FRMdensexpansion->setMaximumSize(QSize(400, 2000));

//            Density type

    FRMdensdensity = new QGroupBox(tr("Data origin"),FRMdensexpansion);
    FRMdensdensity->setMaximumSize(QSize(400, 2000));

    RBTdensExact = new QRadioButton(tr("Original density"),FRMdensdensity);
    RBTdensExact->setToolTip(tr("Original density without fitting"));
    connect(RBTdensExact, SIGNAL(toggled (bool)), this, SLOT(RBTdenslexact_changed()));

    RBTdensRep1 = new QRadioButton(tr("Fitted density"),FRMdensdensity);
    RBTdensRep1->setToolTip(tr("Density expansion in atomic contributions"));
    RBTdensRep1->setChecked(true);
    connect(RBTdensRep1, SIGNAL(toggled (bool)), this, SLOT(RBTdenslexact_changed()));

//          Full density and deformations

    FRMdensdeformation = new QGroupBox(tr("Density options"),FRMdensexpansion);
    FRMdensdeformation->setMaximumSize(QSize(400, 2000));

    RBTdensfulldensity =new QRadioButton(tr("Full electron density"),FRMdensdeformation);
    RBTdensfulldensity->setToolTip(tr("Full electron density"));
    RBTdensfulldensity->setChecked(true);
    connect(RBTdensfulldensity, SIGNAL(toggled (bool)), this, SLOT(RBTdenstype_changed()));

    RBTdensdeform =new QRadioButton(tr("Density deformations"),FRMdensdeformation);
    RBTdensdeform->setToolTip(tr("Electron density deformations (atomic terms with l=0 excluded)"));
    connect(RBTdensdeform, SIGNAL(toggled (bool)), this, SLOT(RBTdenstype_changed()));

    RBTdenslrange =new QRadioButton(tr("Contributions to density"),FRMdensdeformation);
    RBTdenslrange->setToolTip(tr("Atomic terms of given range on l"));
    connect(RBTdenslrange, SIGNAL(toggled (bool)), this, SLOT(RBTdenstype_changed()));

//            Atomic terms

    FRMdenslrange = new QGroupBox(tr("Atomic terms"),FRMdensexpansion);
    FRMdenslrange->setMaximumSize(QSize(400, 2000));

    LBLdenslmaxexp=new QLabel(tr("Highest l"),FRMdenslrange);
    SPBdenslmaxexp=new QSpinBox(FRMdenslrange);
    SPBdenslmaxexp->setRange(0, 25);
    SPBdenslmaxexp->setValue(10);
    SPBdenslmaxexp->setMaximumWidth(50);
    SPBdenslmaxexp->setToolTip(tr("Highest l of atomic multipolar expansion"));
    connect(SPBdenslmaxexp, SIGNAL(valueChanged(int)), this, SLOT(SPBdenslmaxexp_changed()));

    LBLldensminexp=new QLabel(tr("Lowest l"),FRMdenslrange);
    SPBdenslminexp=new QSpinBox(FRMdenslrange);
    SPBdenslminexp->setRange(0, 25);
    SPBdenslminexp->setValue(0);
    SPBdenslminexp->setMaximumWidth(50);
    SPBdenslminexp->setEnabled(false);
    SPBdenslminexp->setToolTip(tr("Lowest l of atomic multipolar expansion"));
    connect(SPBdenslminexp, SIGNAL(valueChanged(int)), this, SLOT(SPBdenslminexp_changed()));


//            Derivatives

    FRMdensderivs = new QGroupBox(tr("Derivatives"),page_MED);
    FRMdensderivs->setMaximumSize(QSize(400, 2000));

    CHKdensgrad = new QCheckBox(tr("Gradient"),FRMdensderivs);
    CHKdensgrad->setToolTip(tr("Tabulate density gradient"));
    CHKdensgrad->setChecked(true);
    connect(CHKdensgrad, SIGNAL(stateChanged(int)), this, SLOT(CHKdensgrad_changed(int)));

    CHKdensder2 = new QCheckBox(tr("Second derivatives"),FRMdensderivs);
    CHKdensder2->setToolTip(tr("Tabulate second derivatives of density"));
    connect(CHKdensder2, SIGNAL(stateChanged(int)), this, SLOT(CHKdensder2_changed(int)));

    CHKdenslaplacian = new QCheckBox(tr("Laplacian"),FRMdensderivs);
    CHKdenslaplacian->setToolTip(tr("Tabulate Laplacian of density"));

//            Atomic Fragments

    FRMdensfragments = new QGroupBox(tr("Molecular fragments"),page_MED);
    FRMdensfragments->setMaximumSize(QSize(400, 2000));

    CHKdenslmolec=new QCheckBox(tr("Molecule"),FRMdensfragments);
    CHKdenslmolec->setChecked(true);

    CHKdenslatomics=new QCheckBox(tr("Atomic fragments"),FRMdensfragments);
    connect(CHKdenslatomics, SIGNAL(stateChanged(int)), this, SLOT(CHKfragments_changed(int)));

//            Fuctional group

    CHKdensldensacc=new QCheckBox(tr("Functional group"),FRMdensfragments);
    CHKdensldensacc->setEnabled(true);
    connect(CHKdensldensacc, SIGNAL(stateChanged(int)), this, SLOT(CHKfragments_changed(int)));

    FRMdensatoms = new QGroupBox(tr("Centers"),page_MED);
    FRMdensatoms->setMaximumSize(QSize(400, 2000));
    FRMdensatoms->setVisible(false);

    QRegExp densrx("[1-9][-,\\d]*");
    densvalidator = new QRegExpValidator(densrx, nullpointer);
    denslist = new QStringList();

    LBLdensatoms = new QLabel("1,3-5,10,13-17,...");
    TXTdensatoms = new QLineEdit(FRMdensatoms);
    TXTdensatoms->setValidator(densvalidator);
    TXTdensatoms->setEnabled(true);
    connect(TXTdensatoms, SIGNAL(textChanged(const QString &)), this, SLOT(TXTdensatoms_changed()));

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

//            Grid

    FRMdensgrid = new QGroupBox(tr("Grid"),page_MED);
    FRMdensgrid->setMaximumSize(QSize(400, 2000));
    CHKdensgrid=new QCheckBox(tr("Generate grid"),FRMdensgrid);
    CHKdensgrid->setChecked(true);
    connect(CHKdensgrid, SIGNAL(stateChanged(int)), this, SLOT(CHKdensgrid_changed()));

    FRMdensgridtype = new QGroupBox(tr("Grid type"));
    RBTdens2D=new QRadioButton(tr("2D grid"),FRMdensgridtype);
    RBTdens2D->setChecked(false);
    RBTdens3D=new QRadioButton(tr("3D grid"),FRMdensgridtype);
    RBTdens3D->setChecked(true);
    connect(RBTdens2D, SIGNAL(toggled (bool)), this, SLOT(RBTdens2D3D_changed()));

//          2D Grid

    FRMdensgrid2D = new QGroupBox(tr("2D grid"));
    FRMdensgrid2D->setMaximumSize(QSize(400, 2000));
    FRMdensgrid2D->setVisible(false);

    LBLdensu=new QLabel(tr("u"),FRMdensgrid2D);
    LBLdensv=new QLabel(tr("v"),FRMdensgrid2D);
    LBLdensinf2d=new QLabel(tr("Lowest"),FRMdensgrid2D);
    LBLdenssup2d=new QLabel(tr("Highest"),FRMdensgrid2D);

    TXTdensuinf=new QLineEdit(FRMdensgrid2D);
    TXTdensuinf->setText("-4.0");
    TXTdensuinf->setAlignment(Qt::AlignRight);
    TXTdensuinf->setValidator(myDoubleValidator);
    connect(TXTdensuinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdensusup=new QLineEdit(FRMdensgrid2D);
    TXTdensusup->setText("4.0");
    TXTdensusup->setAlignment(Qt::AlignRight);
    TXTdensusup->setValidator(myDoubleValidator);
    connect(TXTdensusup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdensvinf=new QLineEdit(FRMdensgrid2D);
    TXTdensvinf->setText("-4.0");
    TXTdensvinf->setAlignment(Qt::AlignRight);
    TXTdensvinf->setValidator(myDoubleValidator);
    connect(TXTdensvinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdensvsup=new QLineEdit(FRMdensgrid2D);
    TXTdensvsup->setText("4.0");
    TXTdensvsup->setAlignment(Qt::AlignRight);
    TXTdensvsup->setValidator(myDoubleValidator);
    connect(TXTdensvsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdenssurfpar = new QGroupBox(tr("Parametric equations"),FRMdensgrid2D);
    FRMdenssurfpar->setMaximumSize(QSize(400, 2000));
    FRMdenssurfpar->setVisible(false);

    LBLdensxformula2D=new QLabel(tr("x(u,v) = "),FRMdenssurfpar);
    LBLdensyformula2D=new QLabel(tr("y(u,v) = "),FRMdenssurfpar);
    LBLdenszformula2D=new QLabel(tr("z(u,v) = "),FRMdenssurfpar);

    TXTdensxformula2D=new QLineEdit(FRMdenssurfpar);
    TXTdensxformula2D->setText("u");
    TXTdensyformula2D=new QLineEdit(FRMdenssurfpar);
    TXTdensyformula2D->setText("v");
    TXTdenszformula2D=new QLineEdit(FRMdenssurfpar);
    TXTdenszformula2D->setText("0");

    FRMdenssurftype = new QGroupBox(tr("Surface type"),FRMdensgrid2D);
    FRMdenssurftype->setMaximumSize(QSize(400, 2000));

    RBTdensplane = new QRadioButton(tr("Plane"),FRMdenssurftype);
    RBTdensplane->setChecked(true);
    connect(RBTdensplane, SIGNAL(toggled (bool)), this, SLOT(RBTdensplane_changed()));

    RBTdensothersurf = new QRadioButton(tr("Parametric surface"),FRMdenssurftype);
    RBTdensothersurf->setToolTip(tr("Surface equation supplied in parametric form: x=x(u,v), y=y(u,v), z=z(u,v)"));
    RBTdensothersurf->setChecked(false);

    FRMdensplane2D = new QGroupBox(tr("2D Planes"),FRMdensgrid2D);
    FRMdensplane2D->setMaximumSize(QSize(400, 2000));
    FRMdensplane2D->setVisible(true);

    RBTdensplaneXY = new QRadioButton(tr("XY "),FRMdensgrid2D);
    RBTdensplaneXY->setChecked(true);
    connect(RBTdensplaneXY, SIGNAL(toggled (bool)), this, SLOT(RBTdens2Dplanes_changed()));

    RBTdensplaneXZ = new QRadioButton(tr("XZ "),FRMdensgrid2D);
    RBTdensplaneXZ->setChecked(false);
    connect(RBTdensplaneXZ, SIGNAL(toggled (bool)), this, SLOT(RBTdens2Dplanes_changed()));

    RBTdensplaneYZ = new QRadioButton(tr("YZ "),FRMdensgrid2D);
    RBTdensplaneYZ->setChecked(false);
    connect(RBTdensplaneYZ, SIGNAL(toggled (bool)), this, SLOT(RBTdens2Dplanes_changed()));

    RBTdensplaneABC = new QRadioButton(tr("Other "),FRMdensgrid2D);
    RBTdensplaneABC->setChecked(false);
    connect(RBTdensplaneABC, SIGNAL(toggled (bool)), this, SLOT(RBTdens2Dplanes_changed()));

    FRMdensplaneABC = new QGroupBox(tr("Plane parameters"),FRMdensgrid2D);
    FRMdensplaneABC->setMaximumSize(QSize(400, 2000));
    FRMdensplaneABC->setVisible(false);

    TXTdensplaneA=new QLineEdit(FRMdensplaneABC);
    TXTdensplaneA->setText("0.");
    TXTdensplaneA->setValidator(myDoubleValidator);
    connect(TXTdensplaneA, SIGNAL(editingFinished()), this, SLOT(RBTdens2Dplanes_changed()));

    TXTdensplaneB=new QLineEdit(FRMdensplaneABC);
    TXTdensplaneB->setText("0.");
    TXTdensplaneB->setValidator(myDoubleValidator);
    connect(TXTdensplaneB, SIGNAL(editingFinished()), this, SLOT(RBTdens2Dplanes_changed()));

    TXTdensplaneC=new QLineEdit(FRMdensplaneABC);
    TXTdensplaneC->setText("1.");
    TXTdensplaneC->setValidator(myDoubleValidator);
    connect(TXTdensplaneC, SIGNAL(editingFinished()), this, SLOT(RBTdens2Dplanes_changed()));

    densplanecase = 1;

    FRMdensresol2D = new QGroupBox(tr("Custom resolution"));
    FRMdensresol2D->setHidden(true);

    LBLdensuresol = new QLabel(tr("u"),FRMdensresol2D);
    LBLdensvresol = new QLabel(tr("v"),FRMdensresol2D);

    SPBdensures = new QSpinBox(FRMdensresol2D);
    SPBdensures->setMinimum(4);
    SPBdensures->setMaximum(2049);
    SPBdensures->setValue(129);
    SPBdensures->setSingleStep(10);

    SPBdensvres = new QSpinBox(FRMdensresol2D);
    SPBdensvres->setMinimum(4);
    SPBdensvres->setMaximum(2049);
    SPBdensvres->setValue(129);
    SPBdensvres->setSingleStep(10);

//          3D Grid

    FRMdensgrid3D = new QGroupBox(tr("3D grid"));
    FRMdensgrid3D->setMaximumSize(QSize(400, 2000));
    FRMdensgrid3D->setVisible(true);

    LBLdensx=new QLabel(tr("x"),FRMdensgrid3D);
    LBLdensy=new QLabel(tr("y"),FRMdensgrid3D);
    LBLdensz=new QLabel(tr("z"),FRMdensgrid3D);
    LBLdensinf3d=new QLabel(tr("Lowest"),FRMdensgrid3D);
    LBLdenssup3d=new QLabel(tr("Highest"),FRMdensgrid3D);

    TXTdensxinf=new QLineEdit(FRMdensgrid3D);
    TXTdensxinf->setText("-4.0");
    TXTdensxinf->setAlignment(Qt::AlignRight);
    TXTdensxinf->setValidator(myDoubleValidator);
    connect(TXTdensxinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdensxsup=new QLineEdit(FRMdensgrid3D);
    TXTdensxsup->setText("4.0");
    TXTdensxsup->setAlignment(Qt::AlignRight);
    TXTdensxsup->setValidator(myDoubleValidator);
    connect(TXTdensxsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdensyinf=new QLineEdit(FRMdensgrid3D);
    TXTdensyinf->setText("-4.0");
    TXTdensyinf->setAlignment(Qt::AlignRight);
    TXTdensyinf->setValidator(myDoubleValidator);
    connect(TXTdensyinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdensysup=new QLineEdit(FRMdensgrid3D);
    TXTdensysup->setText("4.0");
    TXTdensysup->setAlignment(Qt::AlignRight);
    TXTdensysup->setValidator(myDoubleValidator);
    connect(TXTdensysup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdenszinf=new QLineEdit(FRMdensgrid3D);
    TXTdenszinf->setText("-4.0");
    TXTdenszinf->setAlignment(Qt::AlignRight);
    TXTdenszinf->setValidator(myDoubleValidator);
    connect(TXTdenszinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdenszsup=new QLineEdit(FRMdensgrid3D);
    TXTdenszsup->setText("4.0");
    TXTdenszsup->setAlignment(Qt::AlignRight);
    TXTdenszsup->setValidator(myDoubleValidator);
    connect(TXTdenszsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdensresol3D = new QGroupBox(tr("Custom resolution"));
    FRMdensresol3D->setHidden(true);

    LBLdensxresol = new QLabel(tr("x"),FRMdensresol3D);
    LBLdensyresol = new QLabel(tr("y"),FRMdensresol3D);
    LBLdenszresol = new QLabel(tr("z"),FRMdensresol3D);

    SPBdensxres = new QSpinBox(FRMdensresol3D);
    SPBdensxres->setMinimum(4);
    SPBdensxres->setMaximum(2049);
    SPBdensxres->setValue(129);
    SPBdensxres->setSingleStep(10);

    SPBdensyres = new QSpinBox(FRMdensresol3D);
    SPBdensyres->setMinimum(4);
    SPBdensyres->setMaximum(2049);
    SPBdensyres->setValue(129);
    SPBdensyres->setSingleStep(10);

    SPBdenszres = new QSpinBox(FRMdensresol3D);
    SPBdenszres->setMinimum(4);
    SPBdenszres->setMaximum(2049);
    SPBdenszres->setValue(129);
    SPBdenszres->setSingleStep(10);

    FRMdensxyz = new QGroupBox(tr("Tabulation"),page_MED);
    FRMdensxyz->setMaximumSize(QSize(400, 2000));

    CHKdensxyz=new QCheckBox(tr("Tabulation points"),FRMdensxyz);
    CHKdensxyz->setChecked(false);
    connect(CHKdensxyz, SIGNAL(stateChanged(int)), this, SLOT(CHKdensxyz_changed()));

    Wtableden=new QWidget(FRMdensxyz);
    SHTxyz = new Sheet(0, 3, 0,true, Wtableden);
    QStringList QSLxyz;
    QSLxyz << "x" << "y" << "z";
    SHTxyz->setHeader(QSLxyz);
    Wtableden->setVisible(false);
    Wtableden->setEnabled(false);

//          Resolution

    FRMdensgridres = new QGroupBox(tr("Resolution"));

    RBTdensrlow=new QRadioButton(tr("Low"),FRMdensgridres);
    RBTdensrlow->setChecked(true);
    connect(RBTdensrlow, SIGNAL(toggled (bool)), this, SLOT(density_resolution_changed()));

    RBTdensrmedium=new QRadioButton(tr("Medium"),FRMdensgridres);
    connect(RBTdensrmedium, SIGNAL(toggled (bool)), this, SLOT(density_resolution_changed()));

    RBTdensrhigh=new QRadioButton(tr("High"),FRMdensgridres);
    connect(RBTdensrhigh, SIGNAL(toggled (bool)), this, SLOT(density_resolution_changed()));

    RBTdensrcustom=new QRadioButton(tr("Custom"),FRMdensgridres);
    connect(RBTdensrcustom, SIGNAL(toggled (bool)), this, SLOT(density_resolution_changed()));

    if (RBTdens3D->isChecked()){
        RBTdensrlow->setToolTip("65x65x65");
        RBTdensrmedium->setToolTip("129x129x129");
        RBTdensrhigh->setToolTip("257x257x257");
    }
    else{
        RBTdensrlow->setToolTip("129x129");
        RBTdensrmedium->setToolTip("257x257");
        RBTdensrhigh->setToolTip("513x513");
    }

//            Generate input only

    FRMdensinput = new QGroupBox(tr("Input only"),page_MED);
    FRMdensinput->setMaximumSize(QSize(400, 2000));

    CHKdensinput=new QCheckBox(tr("Generate input file only"),FRMdensinput);
    CHKdensinput->setChecked(false);
    CHKdensinput->setEnabled(true);
    connect(CHKdensinput, SIGNAL(stateChanged(int)), this, SLOT(CHKdensinput_changed(int)));

//            MPI

    FRMdensmpi = new QGroupBox(tr("Parallel computing"),page_MED);
    FRMdensmpi->setMaximumSize(QSize(400, 2000));

    CHKdensmpi = new QCheckBox(tr("MPI"),FRMdensmpi);
    LBLdensmpi = new QLabel(tr("Number of processors"),FRMdensmpi);

    SPBdensmpi = new QSpinBox(FRMdensmpi);
    SPBdensmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBdensmpi->setValue(1);
    SPBdensmpi->setMaximumWidth(50);
    SPBdensmpi->setToolTip(tr("Number of processors"));

    if (mpi){
        FRMdensmpi->setVisible(true);
        FRMdensmpi->setEnabled(true);
        CHKdensmpi->setChecked(true);
        SPBdensmpi->setEnabled(true);
    }
    else{
        FRMdensmpi->setHidden(true);
        CHKdensmpi->setChecked(false);
        SPBdensmpi->setEnabled(false);
    }
    connect(CHKdensmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKdensmpi_changed(int)));
    connect(SPBdensmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBdensmpi_changed(int)));
//            List file, Stop, Execute
    BTNexecDamden=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_MED);
    BTNexecDamden->setToolTip(tr("Generate grids"));
    connect(BTNexecDamden, SIGNAL(clicked()), this, SLOT(execDamden()));

    if (!BTNstop[1])
        BTNstop[1]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_MED);
    BTNstop[1]->setToolTip(tr("Kill the process"));
    connect(BTNstop[1], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[1])
        BTNtexto[1]=new QPushButton(QIcon(":/images/document_text.png"),"",page_MED);
    BTNtexto[1]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[1], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_MED_layouts()
{
    QLabel *LBLdensityplaneeq=new QLabel(tr("(A x + B y + C z = 0)"),FRMdensplaneABC);
    QLabel *LBLdensityplaneA=new QLabel(tr("A: "),FRMdensplaneABC);
    QLabel *LBLdensityplaneB=new QLabel(tr("B: "),FRMdensplaneABC);
    QLabel *LBLdensityplaneC=new QLabel(tr("C: "),FRMdensplaneABC);

    QHBoxLayout *Layout1=new QHBoxLayout(FRMdensdamdenfilename);
    Layout1->addWidget(TXTdensdamdenfilename);

    QVBoxLayout *Layout2=new QVBoxLayout(FRMdensdensity);
    Layout2->addWidget(RBTdensExact);
    Layout2->addWidget(RBTdensRep1);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addWidget(LBLdenslmaxexp);
    Layout3->addWidget(SPBdenslmaxexp);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(LBLldensminexp);
    Layout4->addWidget(SPBdenslminexp);

    QVBoxLayout *Layout5=new QVBoxLayout(FRMdenslrange);
    Layout5->addLayout(Layout3);
    Layout5->addLayout(Layout4);

    QVBoxLayout *Layout6=new QVBoxLayout(FRMdensdeformation);
    Layout6->addWidget(RBTdensfulldensity);
    Layout6->addWidget(RBTdensdeform);
    Layout6->addWidget(RBTdenslrange);
    Layout6->addWidget(FRMdenslrange);

    QVBoxLayout *Layout7=new QVBoxLayout(FRMdensexpansion);
    Layout7->addWidget(FRMdensdensity);
    Layout7->addWidget(FRMdensdeformation);

    QVBoxLayout *Layout8=new QVBoxLayout(FRMdensderivs);
    Layout8->addWidget(CHKdensgrad);
    Layout8->addWidget(CHKdensder2);
    Layout8->addWidget(CHKdenslaplacian);

    QVBoxLayout *Layout9 = new QVBoxLayout(FRMdensatoms);
    Layout9->addWidget(LBLdensatoms);
    Layout9->addWidget(TXTdensatoms);

    QVBoxLayout *Layout10=new QVBoxLayout(FRMdensfragments);
    Layout10->addWidget(CHKdenslmolec);
    Layout10->addWidget(CHKdenslatomics);
    Layout10->addWidget(CHKdensldensacc);
    Layout10->addWidget(FRMdensatoms);

    QHBoxLayout *Layout11=new QHBoxLayout();
    Layout11->addWidget(CHKdensgrid,0,Qt::AlignCenter);

    QHBoxLayout *Layout12=new QHBoxLayout(FRMdensgridtype);
    Layout12->addWidget(RBTdens2D);
    Layout12->addWidget(RBTdens3D);

    QGridLayout *Layout13=new QGridLayout(FRMdensgrid3D);
    Layout13->addWidget(LBLdensinf3d,1,1,Qt::AlignCenter);
    Layout13->addWidget(LBLdenssup3d,1,2,Qt::AlignCenter);
    Layout13->addWidget(LBLdensx,2,0);
    Layout13->addWidget(TXTdensxinf,2,1);
    Layout13->addWidget(TXTdensxsup,2,2);
    Layout13->addWidget(LBLdensy,3,0);
    Layout13->addWidget(TXTdensyinf,3,1);
    Layout13->addWidget(TXTdensysup,3,2);
    Layout13->addWidget(LBLdensz,4,0);
    Layout13->addWidget(TXTdenszinf,4,1);
    Layout13->addWidget(TXTdenszsup,4,2);

    QGridLayout *Layout14=new QGridLayout();
    Layout14->addWidget(LBLdensinf2d,1,1,Qt::AlignCenter);
    Layout14->addWidget(LBLdenssup2d,1,2,Qt::AlignCenter);
    Layout14->addWidget(LBLdensu,2,0);
    Layout14->addWidget(TXTdensuinf,2,1);
    Layout14->addWidget(TXTdensusup,2,2);
    Layout14->addWidget(LBLdensv,3,0);
    Layout14->addWidget(TXTdensvinf,3,1);
    Layout14->addWidget(TXTdensvsup,3,2);

    QVBoxLayout *Layout15=new QVBoxLayout(FRMdenssurftype);
    Layout15->addWidget(RBTdensplane);
    Layout15->addWidget(RBTdensothersurf);

    QGridLayout *Layout16=new QGridLayout(FRMdensplaneABC);
    Layout16->addWidget(LBLdensityplaneeq,0,0,1,2,Qt::AlignCenter);
    Layout16->addWidget(LBLdensityplaneA,1,0);
    Layout16->addWidget(TXTdensplaneA,1,1);
    Layout16->addWidget(LBLdensityplaneB,2,0);
    Layout16->addWidget(TXTdensplaneB,2,1);
    Layout16->addWidget(LBLdensityplaneC,3,0);
    Layout16->addWidget(TXTdensplaneC,3,1);

    QGridLayout *Layout17=new QGridLayout(FRMdensplane2D);
    Layout17->addWidget(RBTdensplaneXY,0,0,Qt::AlignLeft);
    Layout17->addWidget(RBTdensplaneXZ,0,1,Qt::AlignLeft);
    Layout17->addWidget(RBTdensplaneYZ,1,0,Qt::AlignLeft);
    Layout17->addWidget(RBTdensplaneABC,1,1,Qt::AlignLeft);
    Layout17->addWidget(FRMdensplaneABC,2,0,1,2,Qt::AlignCenter);

    QGridLayout *Layout18=new QGridLayout(FRMdenssurfpar);
    Layout18->addWidget(LBLdensxformula2D,1,0);
    Layout18->addWidget(TXTdensxformula2D,1,1);
    Layout18->addWidget(LBLdensyformula2D,2,0);
    Layout18->addWidget(TXTdensyformula2D,2,1);
    Layout18->addWidget(LBLdenszformula2D,3,0);
    Layout18->addWidget(TXTdenszformula2D,3,1);

    QVBoxLayout *Layout19=new QVBoxLayout(FRMdensgrid2D);
    Layout19->addLayout(Layout14,Qt::AlignCenter);
    Layout19->addWidget(FRMdenssurftype);
    Layout19->addWidget(FRMdensplane2D);
    Layout19->addWidget(FRMdenssurfpar);

    QVBoxLayout *Layout20=new QVBoxLayout();
    Layout20->addWidget(FRMdensgrid3D);
    Layout20->addWidget(FRMdensgrid2D);

    QHBoxLayout *Layout21=new QHBoxLayout();
    Layout21->addWidget(RBTdensrlow);
    Layout21->addWidget(RBTdensrmedium);

    QHBoxLayout *Layout22=new QHBoxLayout();
    Layout22->addWidget(RBTdensrhigh);
    Layout22->addWidget(RBTdensrcustom);

    QGridLayout *Layout23=new QGridLayout(FRMdensresol2D);
    Layout23->addWidget(LBLdensuresol,0,0,Qt::AlignCenter);
    Layout23->addWidget(LBLdensvresol,0,1,Qt::AlignCenter);
    Layout23->addWidget(SPBdensures,1,0);
    Layout23->addWidget(SPBdensvres,1,1);

    QGridLayout *Layout24=new QGridLayout(FRMdensresol3D);
    Layout24->addWidget(LBLdensxresol,0,0,Qt::AlignCenter);
    Layout24->addWidget(LBLdensyresol,0,1,Qt::AlignCenter);
    Layout24->addWidget(LBLdenszresol,0,2,Qt::AlignCenter);
    Layout24->addWidget(SPBdensxres,1,0);
    Layout24->addWidget(SPBdensyres,1,1);
    Layout24->addWidget(SPBdenszres,1,2);

    QVBoxLayout *Layout25=new QVBoxLayout(FRMdensgridres);
    Layout25->addLayout(Layout21);
    Layout25->addLayout(Layout22);
    Layout25->addWidget(FRMdensresol2D);
    Layout25->addWidget(FRMdensresol3D);

    QVBoxLayout *Layout26=new QVBoxLayout(FRMdensgrid);
    Layout26->addLayout(Layout11);
    Layout26->addWidget(FRMdensgridtype);
    Layout26->addLayout(Layout20);
    Layout26->addWidget(FRMdensgridres);

    QVBoxLayout *Layout27 = new QVBoxLayout(FRMdensxyz);
    Layout27->addWidget(CHKdensxyz,0,Qt::AlignCenter);
    Layout27->addWidget(Wtableden,0,Qt::AlignCenter);

    QVBoxLayout *Layout28 = new QVBoxLayout(FRMdensinput);
    Layout28->addWidget(CHKdensinput);

    QHBoxLayout *Layout29 = new QHBoxLayout(FRMdensmpi);
    Layout29->addWidget(CHKdensmpi);
    Layout29->addWidget(LBLdensmpi);
    Layout29->addWidget(SPBdensmpi);

    QHBoxLayout *Layout30=new QHBoxLayout();
    Layout30->addWidget(BTNtexto[1],0,Qt::AlignLeft);
    Layout30->addWidget(BTNstop[1],0,Qt::AlignRight);
    Layout30->addWidget(BTNexecDamden,0,Qt::AlignRight);

    QVBoxLayout *page_MEDLayout = new QVBoxLayout(page_MED);
    page_MEDLayout->addWidget(FRMdensdamdenfilename);
    page_MEDLayout->addWidget(FRMdensexpansion);
    page_MEDLayout->addWidget(FRMdensderivs);
    page_MEDLayout->addWidget(FRMdensfragments);
    page_MEDLayout->addWidget(FRMdensgrid);
    page_MEDLayout->addWidget(FRMdensxyz);
    page_MEDLayout->addWidget(FRMdensinput);
    page_MEDLayout->addWidget(FRMdensmpi);
    page_MEDLayout->addLayout(Layout30);
    page_MEDLayout->addStretch();
}

//    page_MESP: ELECTROSTATIC POTENTIAL
//    ==================================

void MainWindow::page_MESP_widgets()
{

    FRMpotgdampotfilename = new QGroupBox(tr("Output files prefix"),page_MESP);
    FRMpotgdampotfilename->setMaximumSize(QSize(400, 2000));

    TXTpotgdampotfilename = new QLineEdit(FRMpotgdampotfilename);
    connect(TXTpotgdampotfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    CHKpotexact = new QCheckBox(tr("Exact potential"),page_MESP);
    CHKpotexact->setChecked(false);
    CHKpotexact->setToolTip(tr("Electrostatic potential computed from density matrix and basis set"));
    connect(CHKpotexact, SIGNAL(stateChanged(int)), this, SLOT(CHKpotexact_changed(int)));

    FRMpotlmaxexp = new QGroupBox(tr("Highest l in expansion"),page_MESP);
    FRMpotlmaxexp->setMaximumSize(QSize(400, 2000));
    FRMpotlmaxexp->setStyleSheet("QGroupBox::title{padding 2 2000}");

    SPBpotlmaxexp=new QSpinBox(FRMpotlmaxexp);
    SPBpotlmaxexp->setRange(0, 25);
    SPBpotlmaxexp->setValue(10);
    SPBpotlmaxexp->setMaximumWidth(50);
    SPBpotlmaxexp->setToolTip(tr("Size of multipolar expansion"));
    connect(SPBpotlmaxexp, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

    FRMpotlong = new QGroupBox(tr("Long-range"),page_MESP);
    FRMpotlong->setMaximumSize(QSize(400, 2000));

    CHKpotlong = new QCheckBox(tr("Long-range only"),FRMpotlong);
    CHKpotlong->setChecked(false);
    connect(CHKpotlong, SIGNAL(stateChanged(int)), this, SLOT(CHKpotlong_changed()));

    LBLpotlongthreshold = new QLabel(tr("Long-range threshold:")+" 10^",FRMpotlong);

    SPBpotlongthreshold = new QSpinBox(FRMpotlong);
    SPBpotlongthreshold->setRange(-10, 0);
    SPBpotlongthreshold->setValue(-7);
    SPBpotlongthreshold->setMaximumWidth(60);
    connect(SPBpotlongthreshold, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

//            Derivatives

    FRMpotderivs = new QGroupBox(tr("Derivatives"),page_MESP);
    FRMpotderivs->setMaximumSize(QSize(400, 2000));

    CHKpotgrad = new QCheckBox(tr("Gradient"),FRMpotderivs);
    CHKpotgrad->setToolTip(tr("Tabulate electrostatic potential gradient"));
    CHKpotgrad->setChecked(true);
    connect(CHKpotgrad, SIGNAL(stateChanged(int)), this, SLOT(CHKpotgrad_changed(int)));

    CHKpotder2 = new QCheckBox(tr("Second derivatives"),FRMpotderivs);
    CHKpotder2->setToolTip(tr("Tabulate second derivatives of electrostatic potential"));
    connect(CHKpotder2, SIGNAL(stateChanged(int)), this, SLOT(CHKpotder2_changed(int)));

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

//            Grid

    FRMpotgrid = new QGroupBox(tr("Grid"),page_MESP);
    FRMpotgrid->setMaximumSize(QSize(400, 2000));

    CHKpotgrid=new QCheckBox(tr("Generate grid"),FRMpotgrid);
    CHKpotgrid->setChecked(true);
    connect(CHKpotgrid, SIGNAL(stateChanged(int)), this, SLOT(CHKpotgrid_changed()));

    FRMpotgridtype = new QGroupBox(tr("Grid type"));

    RBTpot2D=new QRadioButton(tr("2D grid"),FRMpotgridtype);
    RBTpot2D->setChecked(false);
    RBTpot3D=new QRadioButton(tr("3D grid"),FRMpotgridtype);
    RBTpot3D->setChecked(true);
    connect(RBTpot2D, SIGNAL(toggled (bool)), this, SLOT(RBTpot2D3D_changed()));

//          2D Grid

    FRMpotgrid2D = new QGroupBox(tr("2D grid"));
    FRMpotgrid2D->setMaximumSize(QSize(400, 2000));
    FRMpotgrid2D->setVisible(false);

    LBLpotu=new QLabel(tr("u"),FRMpotgrid2D);
    LBLpotv=new QLabel(tr("v"),FRMpotgrid2D);
    LBLpotinf2d=new QLabel(tr("Lowest"),FRMpotgrid2D);
    LBLpotsup2d=new QLabel(tr("Highest"),FRMpotgrid2D);

    TXTpotuinf=new QLineEdit(FRMpotgrid2D);
    TXTpotuinf->setText("-4.0");
    TXTpotuinf->setAlignment(Qt::AlignRight);
    TXTpotuinf->setValidator(myDoubleValidator);
    connect(TXTpotuinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotusup=new QLineEdit(FRMpotgrid2D);
    TXTpotusup->setText("4.0");
    TXTpotusup->setAlignment(Qt::AlignRight);
    TXTpotusup->setValidator(myDoubleValidator);
    connect(TXTpotusup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotvinf=new QLineEdit(FRMpotgrid2D);
    TXTpotvinf->setText("-4.0");
    TXTpotvinf->setAlignment(Qt::AlignRight);
    TXTpotvinf->setValidator(myDoubleValidator);
    connect(TXTpotvinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotvsup=new QLineEdit(FRMpotgrid2D);
    TXTpotvsup->setText("4.0");
    TXTpotvsup->setAlignment(Qt::AlignRight);
    TXTpotvsup->setValidator(myDoubleValidator);
    connect(TXTpotvsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMpotsurfpar = new QGroupBox(tr("Parametric equations"),FRMdensgrid2D);
    FRMpotsurfpar->setMaximumSize(QSize(400, 2000));
    FRMpotsurfpar->setVisible(false);

    LBLpotxformula2D=new QLabel(tr("x(u,v) = "),FRMpotsurfpar);
    LBLpotyformula2D=new QLabel(tr("y(u,v) = "),FRMpotsurfpar);
    LBLpotzformula2D=new QLabel(tr("z(u,v) = "),FRMpotsurfpar);

    TXTpotxformula2D=new QLineEdit(FRMpotsurfpar);
    TXTpotxformula2D->setText("u");

    TXTpotyformula2D=new QLineEdit(FRMpotsurfpar);
    TXTpotyformula2D->setText("v");

    TXTpotzformula2D=new QLineEdit(FRMpotsurfpar);
    TXTpotzformula2D->setText("0");

    FRMpotsurftype = new QGroupBox(tr("Surface type"),FRMpotgrid2D);
    FRMpotsurftype->setMaximumSize(QSize(400, 2000));

    RBTpotplane = new QRadioButton(tr("Plane"),FRMpotsurftype);
    RBTpotplane->setChecked(true);
    connect(RBTpotplane, SIGNAL(toggled (bool)), this, SLOT(RBTpotplane_changed()));

    RBTpotothersurf = new QRadioButton(tr("Parametric surface"),FRMpotsurftype);
    RBTpotothersurf->setToolTip(tr("Surface equation supplied in parametric form: x=x(u,v), y=y(u,v), z=z(u,v)"));
    RBTpotothersurf->setChecked(false);

    FRMpotplane2D = new QGroupBox(tr("2D Planes"),FRMpotgrid2D);
    FRMpotplane2D->setMaximumSize(QSize(400, 2000));
    FRMpotplane2D->setVisible(true);

    RBTpotplaneXY = new QRadioButton(tr("XY "),FRMpotgrid2D);
    RBTpotplaneXY->setChecked(true);
    connect(RBTpotplaneXY, SIGNAL(toggled (bool)), this, SLOT(RBTpot2Dplanes_changed()));

    RBTpotplaneXZ = new QRadioButton(tr("XZ "),FRMpotgrid2D);
    RBTpotplaneXZ->setChecked(false);
    connect(RBTpotplaneXZ, SIGNAL(toggled (bool)), this, SLOT(RBTpot2Dplanes_changed()));

    RBTpotplaneYZ = new QRadioButton(tr("YZ "),FRMpotgrid2D);
    RBTpotplaneYZ->setChecked(false);
    connect(RBTpotplaneYZ, SIGNAL(toggled (bool)), this, SLOT(RBTpot2Dplanes_changed()));

    RBTpotplaneABC = new QRadioButton(tr("Other "),FRMpotgrid2D);
    RBTpotplaneABC->setChecked(false);
    connect(RBTpotplaneABC, SIGNAL(toggled (bool)), this, SLOT(RBTpot2Dplanes_changed()));

    FRMpotplaneABC = new QGroupBox(tr("Plane parameters"),FRMpotgrid2D);
    FRMpotplaneABC->setMaximumSize(QSize(400, 2000));
    FRMpotplaneABC->setVisible(false);

    TXTpotplaneA=new QLineEdit(FRMpotplaneABC);
    TXTpotplaneA->setText("0.");
    TXTpotplaneA->setValidator(myDoubleValidator);
    connect(TXTpotplaneA, SIGNAL(editingFinished()), this, SLOT(RBTpot2Dplanes_changed()));

    TXTpotplaneB=new QLineEdit(FRMpotplaneABC);
    TXTpotplaneB->setText("0.");
    TXTpotplaneB->setValidator(myDoubleValidator);
    connect(TXTpotplaneB, SIGNAL(editingFinished()), this, SLOT(RBTpot2Dplanes_changed()));

    TXTpotplaneC=new QLineEdit(FRMpotplaneABC);
    TXTpotplaneC->setText("1.");
    TXTpotplaneC->setValidator(myDoubleValidator);
    connect(TXTpotplaneC, SIGNAL(editingFinished()), this, SLOT(RBTpot2Dplanes_changed()));

    potplanecase = 1;

    FRMpotresol2D = new QGroupBox(tr("Custom resolution"));
    FRMpotresol2D->setHidden(true);

    LBLpoturesol = new QLabel(tr("u"),FRMpotresol2D);
    LBLpotvresol = new QLabel(tr("v"),FRMpotresol2D);

    SPBpotures = new QSpinBox(FRMpotresol2D);
    SPBpotures->setMinimum(4);
    SPBpotures->setMaximum(2049);
    SPBpotures->setValue(129);
    SPBpotures->setSingleStep(10);

    SPBpotvres = new QSpinBox(FRMpotresol2D);
    SPBpotvres->setMinimum(4);
    SPBpotvres->setMaximum(2049);
    SPBpotvres->setValue(129);
    SPBpotvres->setSingleStep(10);

//          3D Grid

    FRMpotgrid3D = new QGroupBox(tr("3D grid"));
    FRMpotgrid3D->setMaximumSize(QSize(400, 2000));
    FRMpotgrid3D->setVisible(true);

    LBLpotx=new QLabel(tr("x"),FRMpotgrid3D);
    LBLpoty=new QLabel(tr("y"),FRMpotgrid3D);
    LBLpotz=new QLabel(tr("z"),FRMpotgrid3D);
    LBLpotinf=new QLabel(tr("Lowest"),FRMpotgrid3D);
    LBLpotsup=new QLabel(tr("Highest"),FRMpotgrid3D);

    TXTpotxinf=new QLineEdit(FRMpotgrid3D);
    TXTpotxinf->setText("-4.0");
    TXTpotxinf->setAlignment(Qt::AlignRight);
    TXTpotxinf->setValidator(myDoubleValidator);
    connect(TXTpotxinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotxsup=new QLineEdit(FRMpotgrid3D);
    TXTpotxsup->setText("4.0");
    TXTpotxsup->setAlignment(Qt::AlignRight);
    TXTpotxsup->setValidator(myDoubleValidator);
    connect(TXTpotxsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotyinf=new QLineEdit(FRMpotgrid3D);
    TXTpotyinf->setText("-4.0");
    TXTpotyinf->setAlignment(Qt::AlignRight);
    TXTpotyinf->setValidator(myDoubleValidator);
    connect(TXTpotyinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotysup=new QLineEdit(FRMpotgrid3D);
    TXTpotysup->setText("4.0");
    TXTpotysup->setAlignment(Qt::AlignRight);
    TXTpotysup->setValidator(myDoubleValidator);
    connect(TXTpotysup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotzinf=new QLineEdit(FRMpotgrid3D);
    TXTpotzinf->setText("-4.0");
    TXTpotzinf->setAlignment(Qt::AlignRight);
    TXTpotzinf->setValidator(myDoubleValidator);
    connect(TXTpotzinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTpotzsup=new QLineEdit(FRMpotgrid3D);
    TXTpotzsup->setText("4.0");
    TXTpotzsup->setAlignment(Qt::AlignRight);
    TXTpotzsup->setValidator(myDoubleValidator);
    connect(TXTpotzsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMpotresol3D = new QGroupBox(tr("Custom resolution"));
    FRMpotresol3D->setHidden(true);

    LBLpotxresol = new QLabel(tr("x"),FRMpotresol3D);
    LBLpotyresol = new QLabel(tr("y"),FRMpotresol3D);
    LBLpotzresol = new QLabel(tr("z"),FRMpotresol3D);

    SPBpotxres = new QSpinBox(FRMpotresol3D);
    SPBpotxres->setMinimum(4);
    SPBpotxres->setMaximum(2049);
    SPBpotxres->setValue(129);
    SPBpotxres->setSingleStep(10);

    SPBpotyres = new QSpinBox(FRMpotresol3D);
    SPBpotyres->setMinimum(4);
    SPBpotyres->setMaximum(2049);
    SPBpotyres->setValue(129);
    SPBpotyres->setSingleStep(10);

    SPBpotzres = new QSpinBox(FRMpotresol3D);
    SPBpotzres->setMinimum(4);
    SPBpotzres->setMaximum(2049);
    SPBpotzres->setValue(129);
    SPBpotzres->setSingleStep(10);

    FRMpotxyz = new QGroupBox(tr("Tabulation"),page_MESP);
    FRMpotxyz->setMaximumSize(QSize(400, 2000));

    CHKpotxyz=new QCheckBox(tr("Tabulation points"),FRMpotxyz);
    CHKpotxyz->setChecked(false);
    connect(CHKpotxyz, SIGNAL(stateChanged(int)), this, SLOT(CHKpotxyz_changed()));

    Wtablepot=new QWidget(FRMpotxyz);
    SHTpotxyz = new Sheet(0, 3, 0,true, Wtablepot);
    QStringList QSLpotxyz;
    QSLpotxyz << "x" << "y" << "z";
    SHTpotxyz->setHeader(QSLpotxyz);
    Wtablepot->setVisible(false);
    Wtablepot->setEnabled(false);

//          Resolution

    FRMpotgridres = new QGroupBox(tr("Resolution"));

    RBTpotrlow=new QRadioButton(tr("Low"),FRMpotgridres);
    RBTpotrlow->setChecked(true);
    connect(RBTpotrlow, SIGNAL(toggled (bool)), this, SLOT(potential_resolution_changed()));

    RBTpotrmedium=new QRadioButton(tr("Medium"),FRMpotgridres);
    connect(RBTpotrmedium, SIGNAL(toggled (bool)), this, SLOT(potential_resolution_changed()));

    RBTpotrhigh=new QRadioButton(tr("High"),FRMpotgridres);
    connect(RBTpotrhigh, SIGNAL(toggled (bool)), this, SLOT(potential_resolution_changed()));

    RBTpotrcustom=new QRadioButton(tr("Custom"),FRMpotgridres);
    connect(RBTpotrcustom, SIGNAL(toggled (bool)), this, SLOT(potential_resolution_changed()));

    if (RBTpot3D->isChecked()){
        RBTpotrlow->setToolTip("65x65x65");
        RBTpotrmedium->setToolTip("129x129x129");
        RBTpotrhigh->setToolTip("257x257x257");
    }
    else{
        RBTpotrlow->setToolTip("129x129");
        RBTpotrmedium->setToolTip("257x257");
        RBTpotrhigh->setToolTip("513x513");
    }

//          Generate input only

    FRMpotinput = new QGroupBox(tr("Input only"),page_MESP);
    FRMpotinput->setMaximumSize(QSize(400, 2000));

    CHKpotinput=new QCheckBox(tr("Generate input file only"),FRMpotinput);
    CHKpotinput->setChecked(false);
    CHKpotinput->setEnabled(true);
    connect(CHKpotinput, SIGNAL(stateChanged(int)), this, SLOT(CHKpotinput_changed(int)));

//            MPI

    FRMpotmpi = new QGroupBox(tr("Parallel computing"),page_MESP);
    FRMpotmpi->setMaximumSize(QSize(400, 2000));

    CHKpotmpi = new QCheckBox(tr("MPI"),FRMpotmpi);

    LBLpotmpi = new QLabel(tr("Number of processors"),FRMpotmpi);

    SPBpotmpi = new QSpinBox(FRMpotmpi);
    SPBpotmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBpotmpi->setValue(1);
    SPBpotmpi->setMaximumWidth(50);
    SPBpotmpi->setToolTip(tr("Number of processors"));

    if (mpi){
        CHKpotmpi->setChecked(true);
        SPBpotmpi->setEnabled(true);
    }
    else{
        FRMpotmpi->setHidden(true);
        CHKpotmpi->setChecked(false);
        SPBpotmpi->setEnabled(false);
    }
    connect(CHKpotmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKpotmpi_changed(int)));
    connect(SPBpotmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBpotmpi_changed(int)));

//            List file, Stop, Execute

    BTNexecDampot=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_MESP);
    BTNexecDampot->setToolTip(tr("Electrostatic potential calculation"));
    connect(BTNexecDampot, SIGNAL(clicked()), this, SLOT(execDampot()));

    if (!BTNstop[2])
        BTNstop[2]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_MESP);
    BTNstop[2]->setToolTip(tr("Kill the process"));
    connect(BTNstop[2], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[2])
        BTNtexto[2]=new QPushButton(QIcon(":/images/document_text.png"),"",page_MESP);
    BTNtexto[2]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[2], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_MESP_layouts()
{
    QLabel *LBLpotplaneeq=new QLabel(tr("(A x + B y + C z = 0)"),FRMpotplaneABC);
    QLabel *LBLpotplaneA=new QLabel(tr("A: "),FRMpotplaneABC);
    QLabel *LBLpotplaneB=new QLabel(tr("B: "),FRMpotplaneABC);
    QLabel *LBLpotplaneC=new QLabel(tr("C: "),FRMpotplaneABC);

    QHBoxLayout *Layout1=new QHBoxLayout(FRMpotgdampotfilename);
    Layout1->addWidget(TXTpotgdampotfilename);

    QHBoxLayout *Layout2=new QHBoxLayout(FRMpotlmaxexp);
    Layout2->addWidget(SPBpotlmaxexp);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addWidget(CHKpotlong);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(LBLpotlongthreshold);
    Layout4->addWidget(SPBpotlongthreshold);
    Layout4->setAlignment(Qt::AlignCenter);

    QVBoxLayout *Layout5=new QVBoxLayout(FRMpotlong);
    Layout5->addLayout(Layout3);
    Layout5->addLayout(Layout4);

    QVBoxLayout *Layout6=new QVBoxLayout(FRMpotderivs);
    Layout6->addWidget(CHKpotgrad);
    Layout6->addWidget(CHKpotder2);

    QHBoxLayout *Layout7=new QHBoxLayout();
    Layout7->addWidget(CHKpotgrid,0,Qt::AlignCenter);

    QHBoxLayout *Layout8=new QHBoxLayout(FRMpotgridtype);
    Layout8->addWidget(RBTpot2D);
    Layout8->addWidget(RBTpot3D);

    QGridLayout *Layout9=new QGridLayout(FRMpotgrid3D);
    Layout9->addWidget(LBLpotinf,1,1,Qt::AlignCenter);
    Layout9->addWidget(LBLpotsup,1,2,Qt::AlignCenter);
    Layout9->addWidget(LBLpotx,2,0);
    Layout9->addWidget(TXTpotxinf,2,1);
    Layout9->addWidget(TXTpotxsup,2,2);
    Layout9->addWidget(LBLpoty,3,0);
    Layout9->addWidget(TXTpotyinf,3,1);
    Layout9->addWidget(TXTpotysup,3,2);
    Layout9->addWidget(LBLpotz,4,0);
    Layout9->addWidget(TXTpotzinf,4,1);
    Layout9->addWidget(TXTpotzsup,4,2);

    QGridLayout *Layout10=new QGridLayout();
    Layout10->addWidget(LBLpotinf2d,1,1,Qt::AlignCenter);
    Layout10->addWidget(LBLpotsup2d,1,2,Qt::AlignCenter);
    Layout10->addWidget(LBLpotu,2,0);
    Layout10->addWidget(TXTpotuinf,2,1);
    Layout10->addWidget(TXTpotusup,2,2);
    Layout10->addWidget(LBLpotv,3,0);
    Layout10->addWidget(TXTpotvinf,3,1);
    Layout10->addWidget(TXTpotvsup,3,2);

    QVBoxLayout *Layout11=new QVBoxLayout(FRMpotsurftype);
    Layout11->addWidget(RBTpotplane);
    Layout11->addWidget(RBTpotothersurf);

    QGridLayout *Layout12=new QGridLayout(FRMpotplaneABC);
    Layout12->addWidget(LBLpotplaneeq,0,0,1,2,Qt::AlignCenter);
    Layout12->addWidget(LBLpotplaneA,1,0);
    Layout12->addWidget(TXTpotplaneA,1,1);
    Layout12->addWidget(LBLpotplaneB,2,0);
    Layout12->addWidget(TXTpotplaneB,2,1);
    Layout12->addWidget(LBLpotplaneC,3,0);
    Layout12->addWidget(TXTpotplaneC,3,1);

    QGridLayout *Layout13=new QGridLayout(FRMpotplane2D);
    Layout13->addWidget(RBTpotplaneXY,0,0,Qt::AlignLeft);
    Layout13->addWidget(RBTpotplaneXZ,0,1,Qt::AlignLeft);
    Layout13->addWidget(RBTpotplaneYZ,1,0,Qt::AlignLeft);
    Layout13->addWidget(RBTpotplaneABC,1,1,Qt::AlignLeft);
    Layout13->addWidget(FRMpotplaneABC,2,0,1,2,Qt::AlignCenter);

    QGridLayout *Layout14=new QGridLayout(FRMpotsurfpar);
    Layout14->addWidget(LBLpotxformula2D,1,0);
    Layout14->addWidget(TXTpotxformula2D,1,1);
    Layout14->addWidget(LBLpotyformula2D,2,0);
    Layout14->addWidget(TXTpotyformula2D,2,1);
    Layout14->addWidget(LBLpotzformula2D,3,0);
    Layout14->addWidget(TXTpotzformula2D,3,1);

    QVBoxLayout *Layout15=new QVBoxLayout(FRMpotgrid2D);
    Layout15->addLayout(Layout10,Qt::AlignCenter);
    Layout15->addWidget(FRMpotsurftype);
    Layout15->addWidget(FRMpotplane2D);
    Layout15->addWidget(FRMpotsurfpar);

    QVBoxLayout *Layout16=new QVBoxLayout();
    Layout16->addWidget(FRMpotgrid3D);
    Layout16->addWidget(FRMpotgrid2D);

    QHBoxLayout *Layout17=new QHBoxLayout();
    Layout17->addWidget(RBTpotrlow);
    Layout17->addWidget(RBTpotrmedium);

    QHBoxLayout *Layout18=new QHBoxLayout();
    Layout18->addWidget(RBTpotrhigh);
    Layout18->addWidget(RBTpotrcustom);

    QGridLayout *Layout19=new QGridLayout(FRMpotresol2D);
    Layout19->addWidget(LBLpoturesol,0,0,Qt::AlignCenter);
    Layout19->addWidget(LBLpotvresol,0,1,Qt::AlignCenter);
    Layout19->addWidget(SPBpotures,1,0);
    Layout19->addWidget(SPBpotvres,1,1);

    QGridLayout *Layout20=new QGridLayout(FRMpotresol3D);
    Layout20->addWidget(LBLpotxresol,0,0,Qt::AlignCenter);
    Layout20->addWidget(LBLpotyresol,0,1,Qt::AlignCenter);
    Layout20->addWidget(LBLpotzresol,0,2,Qt::AlignCenter);
    Layout20->addWidget(SPBpotxres,1,0);
    Layout20->addWidget(SPBpotyres,1,1);
    Layout20->addWidget(SPBpotzres,1,2);

    QVBoxLayout *Layout21=new QVBoxLayout(FRMpotgridres);
    Layout21->addLayout(Layout17);
    Layout21->addLayout(Layout18);
    Layout21->addWidget(FRMpotresol2D);
    Layout21->addWidget(FRMpotresol3D);

    QVBoxLayout *Layout22=new QVBoxLayout(FRMpotgrid);
    Layout22->addLayout(Layout7);
    Layout22->addWidget(FRMpotgridtype);
    Layout22->addLayout(Layout16);
    Layout22->addWidget(FRMpotgridres);

    QVBoxLayout *Layout23 = new QVBoxLayout(FRMpotxyz);
    Layout23->addWidget(CHKpotxyz,0,Qt::AlignCenter);
    Layout23->addWidget(Wtablepot,0,Qt::AlignCenter);

    QVBoxLayout *Layout24 = new QVBoxLayout(FRMpotinput);
    Layout24->addWidget(CHKpotinput);

    QHBoxLayout *Layout25 = new QHBoxLayout(FRMpotmpi);
    Layout25->addWidget(CHKpotmpi);
    Layout25->addWidget(LBLpotmpi);
    Layout25->addWidget(SPBpotmpi);

    QHBoxLayout *Layout26 = new QHBoxLayout();
    Layout26->addWidget(BTNtexto[2],0,Qt::AlignLeft);
    Layout26->addWidget(BTNstop[2],0,Qt::AlignRight);
    Layout26->addWidget(BTNexecDampot,0,Qt::AlignRight);

    QVBoxLayout *page_MESPLayout = new QVBoxLayout(page_MESP);
    page_MESPLayout->addWidget(FRMpotgdampotfilename);
    page_MESPLayout->addWidget(CHKpotexact,0,Qt::AlignCenter);
    page_MESPLayout->addWidget(FRMpotlmaxexp);
    page_MESPLayout->addWidget(FRMpotlong);
    page_MESPLayout->addWidget(FRMpotderivs);
    page_MESPLayout->addWidget(FRMpotxyz);
    page_MESPLayout->addWidget(FRMpotgrid);
    page_MESPLayout->addWidget(FRMpotinput);
    page_MESPLayout->addWidget(FRMpotmpi);
    page_MESPLayout->addLayout(Layout26);
    page_MESPLayout->addStretch();
}

//    page_TOPO: MOLECULAR TOPOGRAPHY
//    ===============================

void MainWindow::page_TOPO_widgets()
{
    FRMtopofilename = new QGroupBox(tr("Output files prefix"),page_TOPO);
    FRMtopofilename->setMaximumSize(QSize(400, 2000));
    TXTtopofilename = new QLineEdit(FRMtopofilename);
    connect(TXTtopofilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMtopotype = new QGroupBox(tr("Topography type"),page_TOPO);
    FRMtopotype->setMaximumSize(QSize(400, 2000));

    RBTtopodensity = new QRadioButton(tr("Molecular density"),FRMtopotype);
    RBTtopodensity->setChecked(true);
    connect(RBTtopodensity, SIGNAL(toggled (bool)), this, SLOT(RBTtopodensity_changed()));

    RBTtopopotential = new QRadioButton(tr("Molecular potential"),FRMtopotype);
    RBTtopopotential->setChecked(false);

    LBLtopolmaxi = new QLabel(tr("Highest l in expansion"),FRMtopotype);

    SPBtopolmaxi=new QSpinBox(FRMdenslrange);
    SPBtopolmaxi->setRange(0, 25);
    SPBtopolmaxi->setValue(10);
    SPBtopolmaxi->setMaximumWidth(50);
    SPBtopolmaxi->setToolTip(tr("Highest l for MED/MESP computation"));

    FRMtopograph = new QGroupBox(tr("Topography mapping"),page_TOPO);
    FRMtopograph->setMaximumSize(QSize(400, 2000));

    CHKtopograph = new QCheckBox(tr("Map critical points"),FRMtopograph);
    CHKtopograph->setChecked(true);
    connect(CHKtopograph, SIGNAL(toggled (bool)), this, SLOT(CHKtopograph_changed()));

    LBLtopoboxl = new QLabel(tr("Box margins size"));
    LBLtopoboxl->setVisible(true);

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

    TXTtopoboxl = new QLineEdit(FRMtopograph);
    TXTtopoboxl->setVisible(true);
    TXTtopoboxl->setText("1.0");
    TXTtopoboxl->setValidator(myDoubleValidator);
    TXTtopoboxl->setToolTip(tr("Box size around default guess points"));

    LBLtopocnvg = new QLabel(tr("Convergence threshold"));

    TXTtopocnvg = new QLineEdit(FRMtopotype);
    TXTtopocnvg->setVisible(true);
    TXTtopocnvg->setText("4.0e-15");
    TXTtopocnvg->setValidator(myDoubleValidator);
    TXTtopocnvg->setToolTip(tr("Convergence threshold. Recommended: 4.0E-12 for MESP and 4.0E-15 for MED"));

    CHKtopoaddguess = new QCheckBox(tr("Add guess points for CPs"),FRMtopograph);
    CHKtopoaddguess->setChecked(false);
    connect(CHKtopoaddguess, SIGNAL(stateChanged(int)), this, SLOT(CHKtopoaddguess_changed(int)));

    FRMtopoxyz = new QGroupBox(tr("Guess points for CPs"),FRMtopograph);
    FRMtopoxyz->setMaximumSize(QSize(400, 2000));

    CHKtopoxyz=new QCheckBox(tr("Add guess points to table"),FRMtopoxyz);
    CHKtopoxyz->setVisible(false);
    CHKtopoxyz->setChecked(false);
    connect(CHKtopoxyz, SIGNAL(stateChanged(int)), this, SLOT(CHKtopoxyz_changed()));

    Wtabledentopo=new QWidget(FRMtopoxyz);
    SHTtopoxyz = new Sheet(0, 3, 0,true, Wtabledentopo);
    QStringList QSLtopoxyz;
    QSLtopoxyz << "x" << "y" << "z";
    SHTtopoxyz->setHeader(QSLtopoxyz);
    Wtabledentopo->setVisible(false);
    Wtabledentopo->setEnabled(false);

    FRMtopoguessfilename = new QGroupBox(tr("Load guess points from file"),FRMtopoxyz);
    FRMtopoguessfilename->setMaximumSize(QSize(400, 2000));
    FRMtopoguessfilename->setVisible(false);

    TXTtopoguessfilename = new QLineEdit(FRMtopoguessfilename);
    TXTtopoguessfilename->setVisible(false);
    TXTtopoguessfilename->setToolTip(tr("Load file with CPs guess points"));
    connect(TXTtopoguessfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    BTNtopoguessfilename = new QToolButton();
    BTNtopoguessfilename->setText(tr("..."));
    BTNtopoguessfilename->setVisible(false);
    connect(BTNtopoguessfilename, SIGNAL(clicked()), this, SLOT(ImportTopguessfilename()));

    LBLtopoboxt = new QLabel(tr("Box size"));

    TXTtopoboxt = new QLineEdit(FRMtopograph);
    TXTtopoboxt->setText("2.0");
    TXTtopoboxt->setValidator(myDoubleValidator);
    LBLtopoboxt->setVisible(false);
    TXTtopoboxt->setVisible(false);
    TXTtopoboxt->setToolTip(tr("Box size around optional guess points"));

    LBLtopostepszt = new QLabel(tr("Step size"));
    TXTtopostepszt = new QLineEdit(FRMtopograph);
    TXTtopostepszt->setText("0.2");
    LBLtopostepszt->setVisible(false);
    TXTtopostepszt->setVisible(false);
    TXTtopostepszt->setValidator(myDoubleValidator);

//        Molecular graph

    CHKtopomolgraph = new QCheckBox(tr("Construct molecular graph"),FRMtopograph);
    CHKtopomolgraph->setChecked(false);
    connect(CHKtopomolgraph, SIGNAL(stateChanged(int)), this, SLOT(CHKtopomolgraph_changed(int)));

    LBLtopoboxg = new QLabel(tr("Box size"));
    LBLtopoboxg->setVisible(false);

    TXTtopoboxg = new QLineEdit(FRMtopograph);
    TXTtopoboxg->setVisible(false);
    TXTtopoboxg->setText("4.5");
    TXTtopoboxg->setValidator(myDoubleValidator);
    TXTtopoboxg->setToolTip(tr("Square box size computed as distance of farthest CP plus this margin size"));

    LBLtopoggradthr = new QLabel(tr("Gradient convergence"));
    LBLtopoggradthr->setVisible(false);

    TXTtopoggradthr = new QLineEdit(FRMtopograph);
    TXTtopoggradthr->setVisible(false);
    TXTtopoggradthr->setText("1.0e-5");
    TXTtopoggradthr->setMinimumWidth(60);
    TXTtopoggradthr->setValidator(myDoubleValidator);
    TXTtopoggradthr->setToolTip(tr("Gradient convergence threshold for finalizing lines or basins"));

//        Basin

    CHKtopobasin = new QCheckBox(tr("Construct 3D atomic basin"),FRMtopograph);
    CHKtopobasin->setChecked(false);
    connect(CHKtopobasin, SIGNAL(stateChanged(int)), this, SLOT(CHKtopobasin_changed(int)));

    LBLtopoboxb = new QLabel(tr("Box size"));
    LBLtopoboxb->setVisible(false);

    TXTtopoboxb = new QLineEdit(FRMtopograph);
    TXTtopoboxb->setVisible(false);
    TXTtopoboxb->setText("6.0");
    TXTtopoboxb->setValidator(myDoubleValidator);
    TXTtopoboxb->setToolTip(tr("Square box size computed as distance of farthest CP plus this margin size"));

    LBLtopofdisp = new QLabel(tr("Initial step"));
    LBLtopofdisp->setVisible(false);

    TXTtopofdisp = new QLineEdit(FRMtopograph);
    TXTtopofdisp->setVisible(false);
    TXTtopofdisp->setEnabled(true);
    TXTtopofdisp->setText("0.02");
    TXTtopofdisp->setValidator(myDoubleValidator);

    CHKtopoexdraw = new QCheckBox(tr("Extra connections"),FRMtopograph);
    CHKtopoexdraw->setChecked(false);
    CHKtopoexdraw->setVisible(false);
    connect(CHKtopoexdraw, SIGNAL(stateChanged(int)), this, SLOT(CHKtopoexdraw_changed(int)));

    LBLtopoexln = new QLabel(tr("Connection threshold"));
    LBLtopoexln->setVisible(false);

    TXTtopoexln = new QLineEdit(FRMtopograph);
    TXTtopoexln->setVisible(false);
    TXTtopoexln->setEnabled(false);
    TXTtopoexln->setText("3.0");
    TXTtopoexln->setValidator(myDoubleValidator);
    TXTtopoexln->setToolTip(tr("Used to determine the number of connecting lines in basin"));

//          Generate input only

    FRMtopoinput = new QGroupBox(tr("Input only"),page_TOPO);
    FRMtopoinput->setMaximumSize(QSize(400, 2000));

    CHKtopoinput=new QCheckBox(tr("Generate input file only"),FRMtopoinput);
    CHKtopoinput->setChecked(false);
    CHKtopoinput->setEnabled(true);
    connect(CHKtopoinput, SIGNAL(stateChanged(int)), this, SLOT(CHKtopoinput_changed(int)));

//            MPI

    FRMtopompi = new QGroupBox(tr("Parallel computing"),page_TOPO);
    FRMtopompi->setMaximumSize(QSize(400, 2000));

    CHKtopompi = new QCheckBox(tr("MPI"),FRMtopompi);

    LBLtopompi = new QLabel(tr("Number of processors"),FRMtopompi);

    SPBtopompi = new QSpinBox(FRMtopompi);
    SPBtopompi->setRange(1, MAX_NUM_PROCESSORS);
    SPBtopompi->setValue(1);
    SPBtopompi->setMaximumWidth(50);
    SPBtopompi->setToolTip(tr("Number of processors"));
    if (mpi){
        CHKtopompi->setChecked(true);
        SPBtopompi->setEnabled(true);
    }
    else{
        FRMtopompi->setHidden(true);
        CHKtopompi->setChecked(false);
        SPBtopompi->setEnabled(false);
    }
    connect(CHKtopompi, SIGNAL(stateChanged(int)), this, SLOT(CHKtopompi_changed(int)));
    connect(SPBtopompi, SIGNAL(valueChanged(int)), this, SLOT(SPBtopompi_changed(int)));

//            List file, Stop, Execute

    BTNexecDamtopo=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_TOPO);
    BTNexecDamtopo->setToolTip(tr("Molecular topography calculation"));
    connect(BTNexecDamtopo, SIGNAL(clicked()), this, SLOT(execDamTopography()));

    if (!BTNstop[8])
        BTNstop[8]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_TOPO);
    BTNstop[8]->setToolTip(tr("Kill the process"));
    connect(BTNstop[8], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[8])
        BTNtexto[8]=new QPushButton(QIcon(":/images/document_text.png"),"",page_TOPO);
    BTNtexto[8]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[8], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_TOPO_layouts()
{
    QHBoxLayout *Layout1=new QHBoxLayout(FRMtopofilename);
    Layout1->addWidget(TXTtopofilename);

    QHBoxLayout *Layout2=new QHBoxLayout();
    Layout2->addWidget(LBLtopolmaxi);
    Layout2->addWidget(SPBtopolmaxi);

    QVBoxLayout *Layout3=new QVBoxLayout(FRMtopotype);
    Layout3->addWidget(RBTtopodensity);
    Layout3->addWidget(RBTtopopotential);
    Layout3->addLayout(Layout2);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addSpacing(20);
    Layout4->addWidget(LBLtopoboxl);
    Layout4->addWidget(TXTtopoboxl);
    Layout4->setAlignment(Qt::AlignRight);

    QHBoxLayout *Layout5=new QHBoxLayout();
    Layout5->addSpacing(20);
    Layout5->addWidget(LBLtopocnvg);
    Layout5->addWidget(TXTtopocnvg);
    Layout5->setAlignment(Qt::AlignRight);

    QHBoxLayout *Layout6=new QHBoxLayout();
    Layout6->addSpacing(20);
    Layout6->addWidget(LBLtopoboxt);
    Layout6->addSpacing(5);
    Layout6->addWidget(TXTtopoboxt);
    Layout6->addSpacing(80);
    Layout6->setAlignment(Qt::AlignRight);

    QHBoxLayout *Layout7=new QHBoxLayout();
    Layout7->addSpacing(20);
    Layout7->addWidget(LBLtopostepszt);
    Layout7->addWidget(TXTtopostepszt);
    Layout7->addSpacing(80);
    Layout7->setAlignment(Qt::AlignRight);

    QVBoxLayout *Layout8 = new QVBoxLayout();
    Layout8->addWidget(CHKtopoxyz,0,Qt::AlignCenter);
    Layout8->addWidget(Wtabledentopo,0,Qt::AlignCenter);

    QHBoxLayout *Layout9=new QHBoxLayout();
    Layout9->addWidget(TXTtopoguessfilename);
    Layout9->addWidget(BTNtopoguessfilename);

    QVBoxLayout *Layout10=new QVBoxLayout(FRMtopoguessfilename);
    Layout10->addLayout(Layout9);

    QVBoxLayout *Layout11=new QVBoxLayout(FRMtopoxyz);
    Layout11->addWidget(CHKtopoaddguess);
    Layout11->addLayout(Layout8);
    Layout11->addWidget(FRMtopoguessfilename);

    QHBoxLayout *Layout12=new QHBoxLayout();
    Layout12->addSpacing(40);
    Layout12->addWidget(LBLtopoboxg);
    Layout12->addWidget(TXTtopoboxg,5,Qt::AlignRight);

    QHBoxLayout *Layout13=new QHBoxLayout();
    Layout13->addSpacing(40);
    Layout13->addWidget(LBLtopoggradthr);
    Layout13->addWidget(TXTtopoggradthr,5,Qt::AlignRight);

    QVBoxLayout *Layout14=new QVBoxLayout();
    Layout14->addWidget(CHKtopomolgraph);
    Layout14->addLayout(Layout12);
    Layout14->addLayout(Layout13);

    QHBoxLayout *Layout15=new QHBoxLayout();
    Layout15->addSpacing(40);
    Layout15->addWidget(LBLtopoboxb);
    Layout15->addWidget(TXTtopoboxb,5,Qt::AlignRight);

    QHBoxLayout *Layout16=new QHBoxLayout();
    Layout16->addWidget(LBLtopofdisp);
    Layout16->addWidget(TXTtopofdisp,5,Qt::AlignRight);
    Layout16->addSpacing(100);
    Layout16->addStretch();

    QHBoxLayout *Layout17=new QHBoxLayout();
    Layout17->addSpacing(40);
    Layout17->addWidget(CHKtopoexdraw,Qt::AlignRight);

    QHBoxLayout *Layout18=new QHBoxLayout();
    Layout18->addSpacing(40);
    Layout18->addWidget(LBLtopoexln);
    Layout18->addWidget(TXTtopoexln,5,Qt::AlignRight);

    QVBoxLayout *Layout19=new QVBoxLayout();
    Layout19->addLayout(Layout15);
    Layout19->addLayout(Layout17);
    Layout19->addLayout(Layout18);

    QVBoxLayout *Layout20=new QVBoxLayout();
    Layout20->addWidget(CHKtopobasin);
    Layout20->addLayout(Layout19);

    QVBoxLayout *Layout21=new QVBoxLayout(FRMtopograph);
    Layout21->addWidget(CHKtopograph);
    Layout21->addLayout(Layout4);
    Layout21->addLayout(Layout5);
    Layout21->addWidget(FRMtopoxyz);
    Layout21->addLayout(Layout6);
    Layout21->addLayout(Layout7);
    Layout21->addLayout(Layout14);
    Layout21->addLayout(Layout20);
    Layout21->addLayout(Layout16);

    QVBoxLayout *Layout22 = new QVBoxLayout(FRMtopoinput);
    Layout22->addWidget(CHKtopoinput);

    QHBoxLayout *Layout23 = new QHBoxLayout(FRMtopompi);
    Layout23->addWidget(CHKtopompi);
    Layout23->addWidget(LBLtopompi);
    Layout23->addWidget(SPBtopompi);

    QHBoxLayout *Layout24=new QHBoxLayout();
    Layout24->addWidget(BTNtexto[8],0,Qt::AlignLeft);
    Layout24->addWidget(BTNstop[8],0,Qt::AlignRight);
    Layout24->addWidget(BTNexecDamtopo,0,Qt::AlignRight);

    QVBoxLayout *page_TOPOLayout = new QVBoxLayout(page_TOPO);
    page_TOPOLayout->addWidget(FRMtopofilename);
    page_TOPOLayout->addWidget(FRMtopotype);
    page_TOPOLayout->addWidget(FRMtopograph);
    page_TOPOLayout->addWidget(FRMtopoinput);
    page_TOPOLayout->addWidget(FRMtopompi);
    page_TOPOLayout->addLayout(Layout24);
    page_TOPOLayout->addStretch();
}

//    page_SGhole: MESP sigma hole
//    ==============================

void MainWindow::page_SGhole_widgets()
{
    FRMImportSGholeden = new QGroupBox(tr("Import density grid from")+":",page_SGhole);
    TXTImportSGholeden = new QLineEdit(FRMImportSGholeden);
    connect(TXTImportSGholeden, SIGNAL(textChanged(const QString &)), this, SLOT(TXTImportSGholeden_changed()));

    FRMSGholefilename = new QGroupBox(tr("Output files prefix"),page_SGhole);
    FRMSGholefilename->setMaximumSize(QSize(400, 2000));

    TXTSGholefilename = new QLineEdit(FRMSGholefilename);
    connect(TXTSGholefilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));


    BTNImportSGholeden = new QToolButton();
    BTNImportSGholeden->setText(tr("..."));
    connect(BTNImportSGholeden, SIGNAL(clicked()), this, SLOT(importFileSGholeden()));

    LBLSGholecontour = new QLabel(tr("Density value"));
    LBLSGholecontour->setVisible(true);
    LBLSGholecontour->setEnabled(false);

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

    TXTSGholecontour = new QLineEdit();
    TXTSGholecontour->setVisible(true);
    TXTSGholecontour->setText("0.001");
    TXTSGholecontour->setValidator(myDoubleValidator);
    TXTSGholecontour->setToolTip(tr("Density value on isosurface (it defines molecular boundaries)"));
    TXTSGholecontour->setEnabled(false);

    LBLSGholelocalextrema = new QLabel(tr("Threshold for local extrema"));
    LBLSGholelocalextrema->setVisible(true);
    LBLSGholelocalextrema->setEnabled(false);
    LBLSGholelocalextrema->setToolTip(tr("Defines regions for local extrema search"));

    SPBSGholelocalextrema = new QSpinBox();
    SPBSGholelocalextrema->setRange(20, 99);
    SPBSGholelocalextrema->setValue(90);
    SPBSGholelocalextrema->setSingleStep(5);
    SPBSGholelocalextrema->setMaximumWidth(60);
    SPBSGholelocalextrema->setEnabled(false);

    LBLSGholelocalpower = new QLabel("x 10<sup>-2</sup>");


    LBLSGholegeomthreshold = new QLabel(tr("Geometry threshold:")+" 10^");
    LBLSGholegeomthreshold->setEnabled(false);
    SPBSGholegeomthreshold = new QSpinBox();
    SPBSGholegeomthreshold->setRange(-10, 0);
    SPBSGholegeomthreshold->setValue(-5);
    SPBSGholegeomthreshold->setMaximumWidth(60);
    SPBSGholegeomthreshold->setEnabled(false);

    LBLSGholelongthreshold = new QLabel(tr("Long-range threshold:")+" 10^");
    LBLSGholelongthreshold->setEnabled(false);
    SPBSGholelongthreshold = new QSpinBox();
    SPBSGholelongthreshold->setRange(-10, 0);
    SPBSGholelongthreshold->setValue(-9);
    SPBSGholelongthreshold->setMaximumWidth(60);
    SPBSGholelongthreshold->setEnabled(false);


    CHKSGhexactMESP=new QCheckBox(tr("Exact potential"));
    CHKSGhexactMESP->setChecked(false);
    CHKSGhexactMESP->setEnabled(false);
    connect(CHKSGhexactMESP, SIGNAL(stateChanged(int)), this, SLOT(CHKSGhexactMESP_changed(int)));

    FRMSGholelmaxexp = new QGroupBox(tr("MESP expansion"),page_SGhole);
    FRMSGholelmaxexp->setMaximumSize(QSize(400, 2000));
    FRMSGholelmaxexp->setEnabled(false);

    LBLSGholelmaxexp = new QLabel(tr("Highest l"));
    LBLSGholelmaxexp->setEnabled(false);
    SPBSGholelmaxexp = new QSpinBox(FRMSGholelmaxexp);
    SPBSGholelmaxexp->setRange(0, 25);
    SPBSGholelmaxexp->setValue(10);
    SPBSGholelmaxexp->setMaximumWidth(50);
    SPBSGholelmaxexp->setToolTip(tr("Highest value of l in expansion of MESP on density isosurface"));
    SPBSGholelmaxexp->setEnabled(false);

//                      Generate input only

    FRMSGholeinput = new QGroupBox(tr("Input only"),page_SGhole);
    FRMSGholeinput->setMaximumSize(QSize(400, 2000));
    FRMSGholeinput->setEnabled(false);

    CHKSGholeinput=new QCheckBox(tr("Generate input file only"),FRMSGholeinput);
    CHKSGholeinput->setChecked(false);
    CHKSGholeinput->setEnabled(false);
    connect(CHKSGholeinput, SIGNAL(stateChanged(int)), this, SLOT(CHKSGholeinput_changed(int)));

//            MPI

    FRMSGholempi = new QGroupBox(tr("Parallel computing"),page_SGhole);
    FRMSGholempi->setMaximumSize(QSize(400, 2000));

    CHKSGholempi = new QCheckBox(tr("MPI"),FRMSGholempi);
    CHKSGholempi->setEnabled(false);

    LBLSGholempi = new QLabel(tr("Number of processors"),FRMSGholempi);
    LBLSGholempi->setEnabled(false);

    SPBSGholempi = new QSpinBox(FRMSGholempi);
    SPBSGholempi->setRange(1, MAX_NUM_PROCESSORS);
    SPBSGholempi->setValue(1);
    SPBSGholempi->setMaximumWidth(60);
    SPBSGholempi->setToolTip(tr("Number of processors"));
    SPBSGholempi->setEnabled(false);
    FRMSGholempi->setEnabled(false);

    if (mpi){
        FRMSGholempi->setVisible(true);
        CHKSGholempi->setChecked(true);
        SPBSGholempi->setEnabled(true);
    }
    else{
        FRMSGholempi->setHidden(true);
        CHKSGholempi->setChecked(false);
        SPBSGholempi->setEnabled(false);
    }
    connect(CHKSGholempi, SIGNAL(stateChanged(int)), this, SLOT(CHKSGholempi_changed(int)));
    connect(SPBSGholempi, SIGNAL(valueChanged(int)), this, SLOT(SPBSGholempi_changed(int)));

    BTNexecDamSGhole=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_SGhole);
    BTNexecDamSGhole->setToolTip(tr("Computes MESP on density isosurface"));
    BTNexecDamSGhole->setEnabled(false);
    connect(BTNexecDamSGhole, SIGNAL(clicked()), this, SLOT(execDamSGhole()));

    if (!BTNstop[12])
        BTNstop[12]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_SGhole);
    BTNstop[12]->setToolTip(tr("Kill the process"));
    connect(BTNstop[12], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[12])
        BTNtexto[12]=new QPushButton(QIcon(":/images/document_text.png"),"",page_SGhole);
    BTNtexto[12]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[12], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_SGhole_layouts()
{
    QHBoxLayout *Layout1=new QHBoxLayout(FRMImportSGholeden);
    Layout1->addWidget(TXTImportSGholeden);
    Layout1->addWidget(BTNImportSGholeden);

    QHBoxLayout *Layout3=new QHBoxLayout(FRMSGholefilename);
    Layout3->addWidget(TXTSGholefilename);

    QHBoxLayout *Layout3b=new QHBoxLayout();
    Layout3b->addWidget(SPBSGholelocalextrema);
    Layout3b->addWidget(LBLSGholelocalpower);

    QGridLayout *Layout4 = new QGridLayout();
    Layout4->addWidget(LBLSGholecontour,0,0);
    Layout4->addWidget(TXTSGholecontour,0,1);
    Layout4->addWidget(LBLSGholegeomthreshold,1,0);
    Layout4->addWidget(SPBSGholegeomthreshold,1,1);
    Layout4->addWidget(LBLSGholelongthreshold,2,0);
    Layout4->addWidget(SPBSGholelongthreshold,2,1);
    Layout4->addWidget(LBLSGholelocalextrema,3,0);
    Layout4->addLayout(Layout3b,3,1);

    QHBoxLayout *Layout5=new QHBoxLayout(FRMSGholelmaxexp);
    Layout5->addWidget(LBLSGholelmaxexp);
    Layout5->addWidget(SPBSGholelmaxexp);

    QVBoxLayout *Layout6 = new QVBoxLayout(FRMSGholeinput);
    Layout6->addWidget(CHKSGholeinput);

    QHBoxLayout *Layout7 = new QHBoxLayout(FRMSGholempi);
    Layout7->addWidget(CHKSGholempi);
    Layout7->addWidget(LBLSGholempi);
    Layout7->addWidget(SPBSGholempi);

    QHBoxLayout *Layout8=new QHBoxLayout();
    Layout8->addWidget(BTNtexto[12],0,Qt::AlignLeft);
    Layout8->addWidget(BTNstop[12],0,Qt::AlignRight);
    Layout8->addWidget(BTNexecDamSGhole,0,Qt::AlignRight);

    QVBoxLayout *page_SGholeLayout = new QVBoxLayout(page_SGhole);
    page_SGholeLayout->addWidget(FRMImportSGholeden);
    page_SGholeLayout->addWidget(FRMSGholefilename);
    page_SGholeLayout->addLayout(Layout4);
    page_SGholeLayout->addWidget(CHKSGhexactMESP,Qt::AlignLeft);
    page_SGholeLayout->addWidget(FRMSGholelmaxexp);
    page_SGholeLayout->addWidget(FRMSGholeinput);
    page_SGholeLayout->addWidget(FRMSGholempi);
    page_SGholeLayout->addLayout(Layout8);
    page_SGholeLayout->addStretch();
}

//    page_HFforces: HELLMANN-FEYNMAN FORCES ON NUCLEI
//    ================================================

void MainWindow::page_HFforces_widgets()
{
    FRMHFgdamforcesfilename = new QGroupBox(tr("Output files prefix"),page_HFforces);
    FRMHFgdamforcesfilename->setMaximumSize(QSize(400, 2000));

    TXTHFgdamforcesfilename = new QLineEdit(FRMHFgdamforcesfilename);
    connect(TXTHFgdamforcesfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMHFdensfragments1 = new QGroupBox(tr("Atomic contributions"),page_HFforces);
    FRMHFdensfragments1->setMaximumSize(QSize(400, 2000));

    CHKHFlatomsel=new QCheckBox(tr("Atomic fragments"),FRMHFdensfragments1);
    connect(CHKHFlatomsel, SIGNAL(stateChanged(int)), this, SLOT(CHKHFlatomsel_changed()));

    FRMHFatomsforces = new QGroupBox(tr("Centers"),page_HFforces);
    FRMHFatomsforces->setMaximumSize(QSize(400, 2000));

    QRegExp forcesrx("[1-9][-,\\d]*");
    HFforcesvalidator = new QRegExpValidator(forcesrx, nullpointer);
    HFforceslist = new QStringList();
    TXTHFforcesatoms = new QLineEdit(FRMHFatomsforces);
    TXTHFforcesatoms->setEnabled(false);
    TXTHFforcesatoms->setValidator(HFforcesvalidator);
    TXTHFforcesatoms->setEnabled(true);
    connect(TXTHFforcesatoms, SIGNAL(textChanged(const QString &)), this, SLOT(TXTHFforcesatoms_changed()));

    BTNexecDamforces=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_HFforces);
    BTNexecDamforces->setToolTip(tr("Compute Hellmann-Feynman forces on nuclei"));
    connect(BTNexecDamforces, SIGNAL(clicked()), this, SLOT(execDamforces()));

    if (!BTNstop[3])
        BTNstop[3]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_HFforces);
    BTNstop[3]->setToolTip(tr("Kill the process"));
    connect(BTNstop[3], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[3])
        BTNtexto[3]=new QPushButton(QIcon(":/images/document_text.png"),"",page_HFforces);
    BTNtexto[3]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[3], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_HFforces_layouts()
{
    QLabel *LBLatomsforces = new QLabel("1,3-5,10,13-17,...");
    QHBoxLayout *Layout1=new QHBoxLayout(FRMHFgdamforcesfilename);
    Layout1->addWidget(TXTHFgdamforcesfilename);

    QVBoxLayout *Layout2 = new QVBoxLayout(FRMHFatomsforces);
    Layout2->addWidget(LBLatomsforces);
    Layout2->addWidget(TXTHFforcesatoms);

    QVBoxLayout *Layout3=new QVBoxLayout(FRMHFdensfragments1);
    Layout3->addWidget(CHKHFlatomsel);
    Layout3->addWidget(FRMHFatomsforces);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(BTNtexto[3],0,Qt::AlignLeft);
    Layout4->addWidget(BTNstop[3],0,Qt::AlignRight);
    Layout4->addWidget(BTNexecDamforces,0,Qt::AlignRight);

    QVBoxLayout *page_HFforcesLayout = new QVBoxLayout(page_HFforces);
    page_HFforcesLayout->addWidget(FRMHFgdamforcesfilename);
    page_HFforcesLayout->addWidget(FRMHFdensfragments1);
    page_HFforcesLayout->addLayout(Layout4);
    page_HFforcesLayout->addStretch();
}


//    page_Efield: ELECTRIC FIELD
//    ===========================

void MainWindow::page_Efield_widgets()
{
    FRMeffilename = new QGroupBox(tr("Output files prefix"),page_Efield);
    FRMeffilename->setMaximumSize(QSize(400, 2000));
    TXTeffilename = new QLineEdit(FRMeffilename);
    connect(TXTeffilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMeflmaxexp = new QGroupBox(tr("Highest l in expansion"),page_Efield);
    FRMeflmaxexp->setMaximumSize(QSize(400, 2000));

    SPBeflmaxexp=new QSpinBox(FRMeflmaxexp);
    SPBeflmaxexp->setRange(0, 25);
    SPBeflmaxexp->setValue(10);
    SPBeflmaxexp->setMaximumWidth(50);
    SPBeflmaxexp->setToolTip(tr("Size of multipolar expansion"));
    connect(SPBeflmaxexp, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

    FRMeflong = new QGroupBox(tr("Long-range"),page_Efield);
    FRMeflong->setMaximumSize(QSize(400, 2000));

    CHKeflong = new QCheckBox(tr("Long-range only"),FRMeflong);
    CHKeflong->setChecked(false);
    connect(CHKeflong, SIGNAL(stateChanged(int)), this, SLOT(CHKeflong_changed()));

    LBLeflongthreshold = new QLabel(tr("Long-range threshold:")+" 10^",FRMeflong);

    SPBeflongthreshold = new QSpinBox(FRMeflong);
    SPBeflongthreshold->setRange(-10, 0);
    SPBeflongthreshold->setValue(-9);
    SPBeflongthreshold->setMaximumWidth(60);
    connect(SPBeflongthreshold, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

    FRMeflines = new QGroupBox(tr("Lines"),page_Efield);
    FRMeflines->setMaximumSize(QSize(400, 2000));

    TXTefnumpnt = new QLineEdit(FRMeflines);
    TXTefnumpnt->setValidator(new QIntValidator(this));
    TXTefnumpnt->setText("1000");
    connect(TXTefnumpnt, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefnlinpernuc = new QLineEdit(FRMeflines);
    TXTefnlinpernuc->setValidator(new QIntValidator(this));
    TXTefnlinpernuc->setText("16");
    connect(TXTefnlinpernuc, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

    TXTefdlt0 = new QLineEdit(FRMeflines);
    TXTefdlt0->setValidator(myDoubleValidator);
    TXTefdlt0->setText("0.02");
    TXTefdlt0->setToolTip(tr("Allowed range between 0.1 and 0.001"));
    connect(TXTefdlt0, SIGNAL(textChanged(const QString &)), this, SLOT(TXTefdlt0_changed()));

    CMBefdirectionset = new QComboBox();
    CMBefdirectionset->addItem("0");
    CMBefdirectionset->addItem("1");
    CMBefdirectionset->addItem("2");
    CMBefdirectionset->addItem("3");
    CMBefdirectionset->addItem("4");
    CMBefdirectionset->addItem("5");
    CMBefdirectionset->addItem("6");
    CMBefdirectionset->addItem("7");
    CMBefdirectionset->setCurrentIndex(1);
    CMBefdirectionset->setMaximumWidth(50);
    CMBefdirectionset->setMinimumWidth(50);
    CMBefdirectionset->setToolTip(tr("Based on icosahedron vertices, C2 axes and C3 axes\n")+
            tr("1: vertices (12 points), 2: C3 axes (20 points), 3: C2 axes (30 points)")+
            tr("4: vertices + C3 axes, 5: vertices + C2 axes, 6: C2 axes + C3 axes, 7: vertices + C2 axes + C3 axes"));

    FRMefplottype = new QGroupBox(tr("Plot type"));

    RBTef2D = new QRadioButton(tr("2D grid"),FRMefplottype);
    RBTef2D->setChecked(false);
    connect(RBTef2D, SIGNAL(toggled (bool)), this, SLOT(RBTef2D3D_changed()));

    RBTef3D = new QRadioButton(tr("3D grid"),FRMefplottype);
    RBTef3D->setChecked(true);

    FRMefplot2D = new QGroupBox(tr("2D Plot"),FRMeflines);
    FRMefplot2D->setMaximumSize(QSize(400, 2000));
    FRMefplot2D->setVisible(false);

    TXTefuinf=new QLineEdit(FRMefplot2D);
    TXTefuinf->setText("-4.0");
    TXTefuinf->setAlignment(Qt::AlignRight);
    TXTefuinf->setValidator(myDoubleValidator);
    connect(TXTefuinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefusup=new QLineEdit(FRMefplot2D);
    TXTefusup->setText("4.0");
    TXTefusup->setAlignment(Qt::AlignRight);
    TXTefusup->setValidator(myDoubleValidator);
    connect(TXTefusup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefvinf=new QLineEdit(FRMefplot2D);
    TXTefvinf->setText("-4.0");
    TXTefvinf->setAlignment(Qt::AlignRight);
    TXTefvinf->setValidator(myDoubleValidator);
    connect(TXTefvinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefvsup=new QLineEdit(FRMefplot2D);
    TXTefvsup->setText("4.0");
    TXTefvsup->setAlignment(Qt::AlignRight);
    TXTefvsup->setValidator(myDoubleValidator);
    connect(TXTefvsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMefplane2D = new QGroupBox(tr("2D Planes"),FRMefplot2D);
    FRMefplane2D->setMaximumSize(QSize(400, 2000));
    FRMefplane2D->setVisible(false);

    RBTefplaneXY = new QRadioButton(tr("XY "),FRMefplane2D);
    RBTefplaneXY->setChecked(true);
    connect(RBTefplaneXY, SIGNAL(toggled (bool)), this, SLOT(RBTef2Dplanes_changed()));

    RBTefplaneXZ = new QRadioButton(tr("XZ "),FRMefplane2D);
    RBTefplaneXZ->setChecked(false);
    connect(RBTefplaneXZ, SIGNAL(toggled (bool)), this, SLOT(RBTef2Dplanes_changed()));

    RBTefplaneYZ = new QRadioButton(tr("YZ "),FRMefplane2D);
    RBTefplaneYZ->setChecked(false);
    connect(RBTefplaneYZ, SIGNAL(toggled (bool)), this, SLOT(RBTef2Dplanes_changed()));

    RBTefplaneABC = new QRadioButton(tr("Other "),FRMefplane2D);
    RBTefplaneABC->setChecked(false);
    connect(RBTefplaneABC, SIGNAL(toggled (bool)), this, SLOT(RBTef2Dplanes_changed()));

    FRMefplaneABC = new QGroupBox(tr("Plane parameters"),FRMefplane2D);
    FRMefplaneABC->setMaximumSize(QSize(400, 2000));
    FRMefplaneABC->setVisible(false);

    TXTefplaneA=new QLineEdit(FRMefplaneABC);
    TXTefplaneA->setText("0.");
    TXTefplaneA->setValidator(myDoubleValidator);
    connect(TXTefplaneA, SIGNAL(editingFinished()), this, SLOT(RBTef2Dplanes_changed()));

    TXTefplaneB=new QLineEdit(FRMefplaneABC);
    TXTefplaneB->setText("0.");
    TXTefplaneB->setValidator(myDoubleValidator);
    connect(TXTefplaneB, SIGNAL(editingFinished()), this, SLOT(RBTef2Dplanes_changed()));

    TXTefplaneC=new QLineEdit(FRMefplaneABC);
    TXTefplaneC->setText("1.");
    TXTefplaneC->setValidator(myDoubleValidator);
    connect(TXTefplaneC, SIGNAL(editingFinished()), this, SLOT(RBTef2Dplanes_changed()));

    TXTefuvratio=new QLineEdit(FRMefplaneABC);
    TXTefuvratio->setText("1.");
    TXTefuvratio->setValidator(myDoubleValidator);
    TXTefuvratio->setToolTip(tr("Change this value to modify the slopes distribution of starting lines around atoms"));

    FRMefplot3D = new QGroupBox(tr("3D Plot"),FRMeflines);
    FRMefplot3D->setMaximumSize(QSize(400, 2000));
    FRMefplot3D->setVisible(true);

    TXTefxinf=new QLineEdit(FRMefplot3D);
    TXTefxinf->setText("-4.0");
    TXTefxinf->setAlignment(Qt::AlignRight);
    TXTefxinf->setValidator(myDoubleValidator);
    connect(TXTefxinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefxsup=new QLineEdit(FRMefplot3D);
    TXTefxsup->setText("4.0");
    TXTefxsup->setAlignment(Qt::AlignRight);
    TXTefxsup->setValidator(myDoubleValidator);
    connect(TXTefxsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefyinf=new QLineEdit(FRMefplot3D);
    TXTefyinf->setText("-4.0");
    TXTefyinf->setAlignment(Qt::AlignRight);
    TXTefyinf->setValidator(myDoubleValidator);
    connect(TXTefyinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefysup=new QLineEdit(FRMefplot3D);
    TXTefysup->setText("4.0");
    TXTefysup->setAlignment(Qt::AlignRight);
    TXTefysup->setValidator(myDoubleValidator);
    connect(TXTefysup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefzinf=new QLineEdit(FRMefplot3D);
    TXTefzinf->setText("-4.0");
    TXTefzinf->setAlignment(Qt::AlignRight);
    TXTefzinf->setValidator(myDoubleValidator);
    connect(TXTefzinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTefzsup=new QLineEdit(FRMefplot3D);
    TXTefzsup->setText("4.0");
    TXTefzsup->setAlignment(Qt::AlignRight);
    TXTefzsup->setValidator(myDoubleValidator);
    connect(TXTefzsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMefextralines = new QGroupBox();
    FRMefextralines->setMaximumSize(QSize(400, 2000));

    CHKefextralines = new QCheckBox(tr("Extra lines"),FRMefextralines);
    CHKefextralines->setChecked(false);
    connect(CHKefextralines, SIGNAL(stateChanged(int)), this, SLOT(CHKefextralines_changed()));

    FRMefxyz = new QGroupBox();
    FRMefxyz->setMaximumSize(QSize(400, 2000));

    CHKefxyz=new QCheckBox(tr("Add starting points to table"),FRMefxyz);
    CHKefxyz->setVisible(false);
    CHKefxyz->setChecked(false);
    connect(CHKefxyz, SIGNAL(stateChanged(int)), this, SLOT(CHKefxyz_changed()));

    Wtabledenef=new QWidget(FRMefxyz);
    SHTefxyz = new Sheet(0, 4, 0,true, Wtabledenef);
    QStringList QSLefxyz;
    QSLefxyz << "cntr" << "x" << "y" << "z";
    SHTefxyz->setHeader(QSLefxyz);
    Wtabledenef->setVisible(false);
    Wtabledenef->setEnabled(false);

    FRMefuv = new QGroupBox();
    FRMefuv->setMaximumSize(QSize(400, 2000));
    Wtable2ef=new QWidget(FRMefuv);
    SHTefuv = new Sheet(0, 3, 0,true, Wtable2ef);
    QStringList QSLefuv;
    QSLefuv << "cntr" << "u" << "v" ;
    SHTefuv->setHeader(QSLefuv);
    Wtable2ef->setVisible(false);
    Wtable2ef->setEnabled(false);

    LBLeffilelines=new QLabel(tr("Read from file:"),FRMefextralines);
    LBLeffilelines->setEnabled(false);
    LBLeffilelines->setVisible(false);

    TXTeffilelines = new QLineEdit(FRMeflines);
    TXTeffilelines->setEnabled(false);
    TXTeffilelines->setVisible(false);
    connect(TXTeffilelines, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    BTNeffilelines = new QToolButton(FRMeflines);
    BTNeffilelines->setText(tr("..."));
    BTNeffilelines->setEnabled(false);
    BTNeffilelines->setVisible(false);
    connect(BTNeffilelines, SIGNAL(clicked()), this, SLOT(BTNeffilelines_clicked()));

//          Generate input only

    FRMefinput = new QGroupBox(tr("Input only"),page_Efield);
    FRMefinput->setMaximumSize(QSize(400, 2000));

    CHKefinput=new QCheckBox(tr("Generate input file only"),FRMefinput);
    CHKefinput->setChecked(false);
    CHKefinput->setEnabled(true);
    connect(CHKefinput, SIGNAL(stateChanged(int)), this, SLOT(CHKefinput_changed(int)));

//            MPI

    FRMefmpi = new QGroupBox(tr("Parallel computing"),page_Efield);
    FRMefmpi->setMaximumSize(QSize(400, 2000));

    CHKefmpi = new QCheckBox(tr("MPI"),FRMefmpi);

    LBLefmpi = new QLabel(tr("Number of processors"),FRMefmpi);

    SPBefmpi = new QSpinBox(FRMefmpi);
    SPBefmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBefmpi->setValue(1);
    SPBefmpi->setMaximumWidth(50);
    SPBefmpi->setToolTip(tr("Number of processors"));
    if (mpi){
        FRMefmpi->setVisible(true);
        FRMefmpi->setEnabled(true);
        CHKefmpi->setChecked(true);
        SPBefmpi->setEnabled(true);
}
    else{
        FRMefmpi->setHidden(true);
        CHKefmpi->setChecked(false);
        SPBefmpi->setEnabled(false);
    }
    connect(CHKefmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKefmpi_changed(int)));
    connect(SPBefmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBefmpi_changed(int)));

//            List file, Stop, Execute

    BTNexecDamfield=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_Efield);
    BTNexecDamfield->setToolTip(tr("Electric field calculation"));
    connect(BTNexecDamfield, SIGNAL(clicked()), this, SLOT(execDamfield()));

    if (!BTNstop[4])
        BTNstop[4]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_Efield);
    BTNstop[4]->setToolTip(tr("Kill the process"));
    connect(BTNstop[4], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[4])
        BTNtexto[4]=new QPushButton(QIcon(":/images/document_text.png"),"",page_Efield);
    BTNtexto[4]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[4], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_Efield_layouts()
{
    QLabel *LBLefnumpnt = new QLabel(tr("Highest number of points"),FRMeflines);
    QLabel *LBLefnlinpernuc = new QLabel(tr("Number of lines per nucleus"),FRMeflines);
    QLabel *LBLefdlt0 = new QLabel(tr("Stride length"),FRMeflines);
    QLabel *LBLefdirectionset = new QLabel(tr("Set of starting directions"));
    QLabel *LBLefu=new QLabel(tr("u"),FRMefplot2D);
    QLabel *LBLefv=new QLabel(tr("v"),FRMefplot2D);
    QLabel *LBLefinfa = new QLabel(tr("Lowest"));
    QLabel *LBLefsupa = new QLabel(tr("Highest"));
    QLabel *LBLefplaneeq=new QLabel(tr("(A x + B y + C z = 0)"),FRMefplaneABC);
    QLabel *LBLefplaneA=new QLabel(tr("A: "),FRMefplaneABC);
    QLabel *LBLefplaneB=new QLabel(tr("B: "),FRMefplaneABC);
    QLabel *LBLefplaneC=new QLabel(tr("C: "),FRMefplaneABC);
    QLabel *LBLefuvratiohead=new QLabel(tr("Starting slopes ratio"),FRMefplaneABC);
    QLabel *LBLefuvratio=new QLabel(tr("v/u: "),FRMefplaneABC);
    QLabel *LBLefx=new QLabel(tr("x"),FRMefplot3D);
    QLabel *LBLefy=new QLabel(tr("y"),FRMefplot3D);
    QLabel *LBLefz=new QLabel(tr("z"),FRMefplot3D);
    QLabel *LBLefinfb = new QLabel(tr("Lowest"));
    QLabel *LBLefsupb = new QLabel(tr("Highest"));

    QHBoxLayout *Layout1=new QHBoxLayout(FRMeffilename);
    Layout1->addWidget(TXTeffilename);

    QHBoxLayout *Layout2=new QHBoxLayout(FRMeflmaxexp);
    Layout2->addWidget(SPBeflmaxexp);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addWidget(CHKeflong);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(LBLeflongthreshold);
    Layout4->addWidget(SPBeflongthreshold);
    Layout4->setAlignment(Qt::AlignCenter);

    QVBoxLayout *Layout5=new QVBoxLayout(FRMeflong);
    Layout5->addLayout(Layout3);
    Layout5->addLayout(Layout4);

    QHBoxLayout *Layout6=new QHBoxLayout();
    Layout6->addWidget(CHKefextralines);
    Layout6->setAlignment(Qt::AlignCenter);

    QHBoxLayout *Layout7=new QHBoxLayout();
    Layout7->addWidget(LBLefnumpnt);
    Layout7->addWidget(TXTefnumpnt);

    QHBoxLayout *Layout8=new QHBoxLayout();
    Layout8->addWidget(LBLefnlinpernuc);
    Layout8->addWidget(TXTefnlinpernuc);

    QHBoxLayout *Layout9=new QHBoxLayout();
    Layout9->addWidget(LBLefdirectionset);
    Layout9->addWidget(CMBefdirectionset);

    QHBoxLayout *Layout10=new QHBoxLayout();
    Layout10->addWidget(LBLefdlt0);
    Layout10->addWidget(TXTefdlt0);

    QHBoxLayout *Layout11=new QHBoxLayout(FRMefplottype);
    Layout11->addWidget(RBTef2D);
    Layout11->addWidget(RBTef3D);

    QGridLayout *Layout12=new QGridLayout(FRMefplaneABC);
    Layout12->addWidget(LBLefplaneeq,0,0,1,2,Qt::AlignCenter);
    Layout12->addWidget(LBLefplaneA,1,0);
    Layout12->addWidget(TXTefplaneA,1,1);
    Layout12->addWidget(LBLefplaneB,2,0);
    Layout12->addWidget(TXTefplaneB,2,1);
    Layout12->addWidget(LBLefplaneC,3,0);
    Layout12->addWidget(TXTefplaneC,3,1);
    Layout12->addWidget(LBLefuvratiohead,4,0,1,2,Qt::AlignCenter);
    Layout12->addWidget(LBLefuvratio,5,0);
    Layout12->addWidget(TXTefuvratio,5,1);

    QGridLayout *Layout13=new QGridLayout(FRMefplane2D);
    Layout13->addWidget(RBTefplaneXY,0,0,Qt::AlignLeft);
    Layout13->addWidget(RBTefplaneXZ,0,1,Qt::AlignLeft);
    Layout13->addWidget(RBTefplaneYZ,1,0,Qt::AlignLeft);
    Layout13->addWidget(RBTefplaneABC,1,1,Qt::AlignLeft);
    Layout13->addWidget(FRMefplaneABC,2,0,1,2,Qt::AlignCenter);

    QGridLayout *Layout14=new QGridLayout();
    Layout14->addWidget(LBLefinfa,1,1,Qt::AlignCenter);
    Layout14->addWidget(LBLefsupa,1,2,Qt::AlignCenter);
    Layout14->addWidget(LBLefu,2,0);
    Layout14->addWidget(TXTefuinf,2,1);
    Layout14->addWidget(TXTefusup,2,2);
    Layout14->addWidget(LBLefv,3,0);
    Layout14->addWidget(TXTefvinf,3,1);
    Layout14->addWidget(TXTefvsup,3,2);
    Layout14->addWidget(FRMefplane2D,4,0,1,3,Qt::AlignCenter);

    QVBoxLayout *Layout15=new QVBoxLayout(FRMefplot2D);
    Layout15->addLayout(Layout8);
    Layout15->addLayout(Layout14);

    QGridLayout *Layout16=new QGridLayout();
    Layout16->addWidget(LBLefinfb,1,1,Qt::AlignCenter);
    Layout16->addWidget(LBLefsupb,1,2,Qt::AlignCenter);
    Layout16->addWidget(LBLefx,2,0);
    Layout16->addWidget(TXTefxinf,2,1);
    Layout16->addWidget(TXTefxsup,2,2);
    Layout16->addWidget(LBLefy,3,0);
    Layout16->addWidget(TXTefyinf,3,1);
    Layout16->addWidget(TXTefysup,3,2);
    Layout16->addWidget(LBLefz,4,0);
    Layout16->addWidget(TXTefzinf,4,1);
    Layout16->addWidget(TXTefzsup,4,2);

    QVBoxLayout *Layout17=new QVBoxLayout(FRMefplot3D);
    Layout17->addLayout(Layout9);
    Layout17->addLayout(Layout16);

    QVBoxLayout *Layout18 = new QVBoxLayout(FRMefextralines);
    Layout18->addWidget(CHKefxyz,0,Qt::AlignCenter);
    Layout18->addWidget(Wtabledenef,0,Qt::AlignCenter);
    Layout18->addWidget(Wtable2ef,0,Qt::AlignCenter);

    QHBoxLayout *Layout19=new QHBoxLayout();
    Layout19->addWidget(LBLeffilelines);

    QHBoxLayout *Layout20=new QHBoxLayout();
    Layout20->addWidget(TXTeffilelines);
    Layout20->addWidget(BTNeffilelines);

    QVBoxLayout *Layout21=new QVBoxLayout(FRMeflines);
    Layout21->addLayout(Layout7);
    Layout21->addLayout(Layout10);
    Layout21->addWidget(FRMefplottype);
    Layout21->addWidget(FRMefplot2D);
    Layout21->addWidget(FRMefplot3D);
    Layout21->addLayout(Layout6);
    Layout21->addWidget(FRMefextralines);
    Layout21->addLayout(Layout19);
    Layout21->addLayout(Layout20);

    QVBoxLayout *Layout22 = new QVBoxLayout(FRMefinput);
    Layout22->addWidget(CHKefinput);

    QHBoxLayout *Layout23 = new QHBoxLayout(FRMefmpi);
    Layout23->addWidget(CHKefmpi);
    Layout23->addWidget(LBLefmpi);
    Layout23->addWidget(SPBefmpi);

    QHBoxLayout *Layout24=new QHBoxLayout();
    Layout24->addWidget(BTNtexto[4],0,Qt::AlignLeft);
    Layout24->addWidget(BTNstop[4],0,Qt::AlignRight);
    Layout24->addWidget(BTNexecDamfield,0,Qt::AlignRight);

    QVBoxLayout *page_EfieldLayout = new QVBoxLayout(page_Efield);
    page_EfieldLayout->addWidget(FRMeffilename);
    page_EfieldLayout->addWidget(FRMeflmaxexp);
    page_EfieldLayout->addWidget(FRMeflong);
    page_EfieldLayout->addWidget(FRMeflines);
    page_EfieldLayout->addWidget(FRMefinput);
    page_EfieldLayout->addWidget(FRMefmpi);
    page_EfieldLayout->addLayout(Layout24);
    page_EfieldLayout->addStretch();
}

//    page_densgrad: DENSITY GRADIENT
//    ===============================

void MainWindow::page_densgrad_widgets()
{
    FRMdgfilename = new QGroupBox(tr("Output files prefix"),page_densgrad);
    FRMdgfilename->setMaximumSize(QSize(400, 2000));

    TXTdgfilename = new QLineEdit(FRMdgfilename);
    connect(TXTdgfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdglmaxexp = new QGroupBox(tr("Highest l in expansion"),page_densgrad);
    FRMdglmaxexp->setMaximumSize(QSize(400, 2000));

    SPBdglmaxexp=new QSpinBox(FRMdglmaxexp);
    SPBdglmaxexp->setRange(0, 25);
    SPBdglmaxexp->setValue(10);
    SPBdglmaxexp->setMaximumWidth(50);
    SPBdglmaxexp->setToolTip(tr("Size of multipolar expansion"));
    connect(SPBdglmaxexp, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

    FRMdglong = new QGroupBox(tr("Long-range"),page_densgrad);
    FRMdglong->setMaximumSize(QSize(400, 2000));

    LBLdglongthreshold = new QLabel(tr("Long-range threshold:")+" 10^",FRMdglong);

    SPBdglongthreshold = new QSpinBox(FRMdglong);
    SPBdglongthreshold->setRange(-10, 0);
    SPBdglongthreshold->setValue(-9);
    SPBdglongthreshold->setMaximumWidth(60);
    connect(SPBdglongthreshold, SIGNAL(valueChanged(int)), this, SLOT(TXTValidate_changed()));

    FRMdglines = new QGroupBox(tr("Lines"),page_densgrad);
    FRMdglines->setMaximumSize(QSize(400, 2000));

    CHKdgextralines = new QCheckBox(tr("Extra lines"));
    CHKdgextralines->setChecked(false);
    connect(CHKdgextralines, SIGNAL(stateChanged(int)), this, SLOT(CHKdgextralines_changed()));

    TXTdgnumpnt = new QLineEdit(FRMdglines);
    TXTdgnumpnt->setValidator(new QIntValidator(this));
    TXTdgnumpnt->setText("1000");
    connect(TXTdgnumpnt, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgnlinpernuc = new QLineEdit(FRMdglines);
    TXTdgnlinpernuc->setValidator(new QIntValidator(this));
    TXTdgnlinpernuc->setText("16");
    connect(TXTdgnlinpernuc, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

    TXTdgdlt0 = new QLineEdit(FRMdglines);
    TXTdgdlt0->setValidator(myDoubleValidator);
    TXTdgdlt0->setText("0.02");
    TXTdgdlt0->setToolTip(tr("Allowed range between 0.1 and 0.001"));
    connect(TXTdgdlt0, SIGNAL(textChanged(const QString &)), this, SLOT(TXTdgdlt0_changed()));

    CMBdgdirectionset = new QComboBox();
    CMBdgdirectionset->addItem("0");
    CMBdgdirectionset->addItem("1");
    CMBdgdirectionset->addItem("2");
    CMBdgdirectionset->addItem("3");
    CMBdgdirectionset->addItem("4");
    CMBdgdirectionset->addItem("5");
    CMBdgdirectionset->addItem("6");
    CMBdgdirectionset->addItem("7");
    CMBdgdirectionset->setCurrentIndex(1);
    CMBdgdirectionset->setMaximumWidth(50);
    CMBdgdirectionset->setMinimumWidth(40);
    CMBdgdirectionset->setToolTip(tr("Based on icosahedron vertices, C2 axes and C3 axes\n")+
            tr("1: vertices (12 points), 2: C3 axes (20 points), 3: C2 axes (30 points)")+
            tr("4: vertices + C3 axes, 5: vertices + C2 axes, 6: C2 axes + C3 axes, 7: vertices + C2 axes + C3 axes"));

    FRMdgplottype = new QGroupBox(tr("Plot type"));

    RBTdg2D = new QRadioButton(tr("2D grid"),FRMdgplottype);
    RBTdg2D->setChecked(false);
    connect(RBTdg2D, SIGNAL(toggled (bool)), this, SLOT(RBTdg2D3D_changed()));

    RBTdg3D = new QRadioButton(tr("3D grid"),FRMdgplottype);
    RBTdg3D->setChecked(true);

    FRMdgplot2D = new QGroupBox(tr("2D Plot"),FRMdglines);
    FRMdgplot2D->setMaximumSize(QSize(400, 2000));
    FRMdgplot2D->setVisible(false);

    TXTdguinf=new QLineEdit(FRMdgplot2D);
    TXTdguinf->setText("-4.0");
    TXTdguinf->setAlignment(Qt::AlignRight);
    TXTdguinf->setValidator(myDoubleValidator);
    connect(TXTdguinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgusup=new QLineEdit(FRMdgplot2D);
    TXTdgusup->setText("4.0");
    TXTdgusup->setAlignment(Qt::AlignRight);
    TXTdgusup->setValidator(myDoubleValidator);
    connect(TXTdgusup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgvinf=new QLineEdit(FRMdgplot2D);
    TXTdgvinf->setText("-4.0");
    TXTdgvinf->setAlignment(Qt::AlignRight);
    TXTdgvinf->setValidator(myDoubleValidator);
    connect(TXTdgvinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgvsup=new QLineEdit(FRMdgplot2D);
    TXTdgvsup->setText("4.0");
    TXTdgvsup->setAlignment(Qt::AlignRight);
    TXTdgvsup->setValidator(myDoubleValidator);
    connect(TXTdgvsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdgplane2D = new QGroupBox(tr("2D Planes"),FRMdgplot2D);
    FRMdgplane2D->setMaximumSize(QSize(400, 2000));
    FRMdgplane2D->setVisible(false);

    RBTdgplaneXY = new QRadioButton(tr("XY "),FRMdgplane2D);
    RBTdgplaneXY->setChecked(true);
    connect(RBTdgplaneXY, SIGNAL(toggled (bool)), this, SLOT(RBTdg2Dplanes_changed()));

    RBTdgplaneXZ = new QRadioButton(tr("XZ "),FRMdgplane2D);
    RBTdgplaneXZ->setChecked(false);
    connect(RBTdgplaneXZ, SIGNAL(toggled (bool)), this, SLOT(RBTdg2Dplanes_changed()));

    RBTdgplaneYZ = new QRadioButton(tr("YZ "),FRMdgplane2D);
    RBTdgplaneYZ->setChecked(false);
    connect(RBTdgplaneYZ, SIGNAL(toggled (bool)), this, SLOT(RBTdg2Dplanes_changed()));

    RBTdgplaneABC = new QRadioButton(tr("Other "),FRMdgplane2D);
    RBTdgplaneABC->setChecked(false);
    connect(RBTdgplaneABC, SIGNAL(toggled (bool)), this, SLOT(RBTdg2Dplanes_changed()));

    FRMdgplaneABC = new QGroupBox(tr("Plane parameters"),FRMdgplane2D);
    FRMdgplaneABC->setMaximumSize(QSize(400, 2000));
    FRMdgplaneABC->setVisible(false);

    TXTdgplaneA=new QLineEdit(FRMdgplaneABC);
    TXTdgplaneA->setText("0.");
    TXTdgplaneA->setValidator(myDoubleValidator);
    connect(TXTdgplaneA, SIGNAL(editingFinished()), this, SLOT(RBTdg2Dplanes_changed()));

    TXTdgplaneB=new QLineEdit(FRMdgplaneABC);
    TXTdgplaneB->setText("0.");
    TXTdgplaneB->setValidator(myDoubleValidator);
    connect(TXTdgplaneB, SIGNAL(editingFinished()), this, SLOT(RBTdg2Dplanes_changed()));

    TXTdgplaneC=new QLineEdit(FRMdgplaneABC);
    TXTdgplaneC->setText("1.");
    TXTdgplaneC->setValidator(myDoubleValidator);
    connect(TXTdgplaneC, SIGNAL(editingFinished()), this, SLOT(RBTdg2Dplanes_changed()));

    TXTdguvratio=new QLineEdit(FRMdgplaneABC);
    TXTdguvratio->setText("1.");
    TXTdguvratio->setValidator(myDoubleValidator);
    TXTdguvratio->setToolTip(tr("Change this value to modify the slopes distribution of starting lines around atoms"));

    FRMdgplot3D = new QGroupBox(tr("3D Plot"),FRMdglines);
    FRMdgplot3D->setMaximumSize(QSize(400, 2000));
    FRMdgplot3D->setVisible(true);

    TXTdgxinf=new QLineEdit(FRMdgplot3D);
    TXTdgxinf->setText("-4.0");
    TXTdgxinf->setAlignment(Qt::AlignRight);
    TXTdgxinf->setValidator(myDoubleValidator);
    connect(TXTdgxinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgxsup=new QLineEdit(FRMdgplot3D);
    TXTdgxsup->setText("4.0");
    TXTdgxsup->setAlignment(Qt::AlignRight);
    TXTdgxsup->setValidator(myDoubleValidator);
    connect(TXTdgxsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgyinf=new QLineEdit(FRMdgplot3D);
    TXTdgyinf->setText("-4.0");
    TXTdgyinf->setAlignment(Qt::AlignRight);
    TXTdgyinf->setValidator(myDoubleValidator);
    connect(TXTdgyinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgysup=new QLineEdit(FRMdgplot3D);
    TXTdgysup->setText("4.0");
    TXTdgysup->setAlignment(Qt::AlignRight);
    TXTdgysup->setValidator(myDoubleValidator);
    connect(TXTdgysup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgzinf=new QLineEdit(FRMdgplot3D);
    TXTdgzinf->setText("-4.0");
    TXTdgzinf->setAlignment(Qt::AlignRight);
    TXTdgzinf->setValidator(myDoubleValidator);
    connect(TXTdgzinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdgzsup=new QLineEdit(FRMdgplot3D);
    TXTdgzsup->setText("4.0");
    TXTdgzsup->setAlignment(Qt::AlignRight);
    TXTdgzsup->setValidator(myDoubleValidator);
    connect(TXTdgzsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdgextralines = new QGroupBox();
    FRMdgextralines->setMaximumSize(QSize(400, 2000));
    CHKdgextralines = new QCheckBox(tr("Extra lines"),FRMdgextralines);
    CHKdgextralines->setChecked(false);
    connect(CHKdgextralines, SIGNAL(stateChanged(int)), this, SLOT(CHKdgextralines_changed()));

    FRMdgxyz = new QGroupBox();
    FRMdgxyz->setMaximumSize(QSize(400, 2000));

    CHKdgxyz=new QCheckBox(tr("Add starting points to table"),FRMdgxyz);
    CHKdgxyz->setVisible(false);
    CHKdgxyz->setChecked(false);
    connect(CHKdgxyz, SIGNAL(stateChanged(int)), this, SLOT(CHKdgxyz_changed()));

    Wtabledendg=new QWidget(FRMdgxyz);
    SHTdgxyz = new Sheet(0, 4, 0,true, Wtabledendg);
    QStringList QSLdgxyz;
    QSLdgxyz << "cntr" << "x" << "y" << "z";
    SHTdgxyz->setHeader(QSLdgxyz);
    Wtabledendg->setVisible(false);
    Wtabledendg->setEnabled(false);

    FRMdguv = new QGroupBox();
    FRMdguv->setMaximumSize(QSize(400, 2000));

    Wtable2dg=new QWidget(FRMdguv);
    SHTdguv = new Sheet(0, 3, 0,true, Wtable2dg);
    QStringList QSLdguv;
    QSLdguv << "cntr" << "u" << "v" ;
    SHTdguv->setHeader(QSLdguv);
    Wtable2dg->setVisible(false);
    Wtable2dg->setEnabled(false);

    LBLdgfilelines=new QLabel(tr("Read from file:"),FRMdglines);
    LBLdgfilelines->setEnabled(false);
    LBLdgfilelines->setVisible(false);

    TXTdgfilelines = new QLineEdit(FRMdglines);
    TXTdgfilelines->setEnabled(false);
    TXTdgfilelines->setVisible(false);
    connect(TXTdgfilelines, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    BTNdgfilelines = new QToolButton(FRMdglines);
    BTNdgfilelines->setText(tr("..."));
    connect(BTNdgfilelines, SIGNAL(clicked()), this, SLOT(BTNdgfilelines_clicked()));

    BTNdgfilelines->setEnabled(false);
    BTNdgfilelines->setVisible(false);

//          Generate input only

    FRMdginput = new QGroupBox(tr("Input only"),page_densgrad);
    FRMdginput->setMaximumSize(QSize(400, 2000));

    CHKdginput=new QCheckBox(tr("Generate input file only"),FRMdginput);
    CHKdginput->setChecked(false);
    CHKdginput->setEnabled(true);
    connect(CHKdginput, SIGNAL(stateChanged(int)), this, SLOT(CHKdginput_changed(int)));

//            MPI

    FRMdgmpi = new QGroupBox(tr("Parallel computing"),page_densgrad);
    FRMdgmpi->setMaximumSize(QSize(400, 2000));

    CHKdgmpi = new QCheckBox(tr("MPI"),FRMdgmpi);

    LBLdgmpi = new QLabel(tr("Number of processors"),FRMdgmpi);

    SPBdgmpi = new QSpinBox(FRMdgmpi);
    SPBdgmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBdgmpi->setValue(1);
    SPBdgmpi->setMaximumWidth(50);
    SPBdgmpi->setToolTip(tr("Number of processors"));
    if (mpi){
        FRMdgmpi->setVisible(true);
        FRMdgmpi->setEnabled(true);
        CHKdgmpi->setChecked(true);
        SPBdgmpi->setEnabled(true);
    }
    else{
        FRMdgmpi->setHidden(true);
        CHKdgmpi->setChecked(false);
        SPBdgmpi->setEnabled(false);
    }
    connect(CHKdgmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKdgmpi_changed(int)));
    connect(SPBdgmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBdgmpi_changed(int)));

//            List file, Stop, Execute

    BTNexecDamdengrad=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_densgrad);
    BTNexecDamdengrad->setToolTip(tr("Electric field calculation"));
    connect(BTNexecDamdengrad, SIGNAL(clicked()), this, SLOT(execDamdengrad()));

    if (!BTNstop[11])
        BTNstop[11]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_densgrad);
    BTNstop[11]->setToolTip(tr("Kill the process"));
    connect(BTNstop[11], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[11])
        BTNtexto[11]=new QPushButton(QIcon(":/images/document_text.png"),"",page_densgrad);
    BTNtexto[11]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[11], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_densgrad_layouts()
{
    QLabel *LBLdgnumpnt = new QLabel(tr("Highest number of points"),FRMdglines);
    QLabel *LBLdgnlinpernuc = new QLabel(tr("Number of lines per nucleus"),FRMdglines);
    QLabel *LBLdgdlt0 = new QLabel(tr("Stride length"),FRMdglines);
    QLabel *LBLdgdirectionset = new QLabel(tr("Set of starting directions"));
    QLabel *LBLdgu=new QLabel(tr("u"),FRMdgplot2D);
    QLabel *LBLdgv=new QLabel(tr("v"),FRMdgplot2D);
    QLabel *LBLdginfa = new QLabel(tr("Lowest"));
    QLabel *LBLdgsupa = new QLabel(tr("Highest"));
    QLabel *LBLdgplaneeq=new QLabel(tr("(A x + B y + C z = 0)"),FRMdgplaneABC);
    QLabel *LBLdgplaneA=new QLabel(tr("A: "),FRMdgplaneABC);
    QLabel *LBLdgplaneB=new QLabel(tr("B: "),FRMdgplaneABC);
    QLabel *LBLdgplaneC=new QLabel(tr("C: "),FRMdgplaneABC);
    QLabel *LBLdguvratiohead=new QLabel(tr("Starting slopes ratio"),FRMdgplaneABC);
    QLabel *LBLdguvratio=new QLabel(tr("v/u: "),FRMdgplaneABC);
    QLabel *LBLdgx=new QLabel(tr("x"),FRMdgplot3D);
    QLabel *LBLdgy=new QLabel(tr("y"),FRMdgplot3D);
    QLabel *LBLdgz=new QLabel(tr("z"),FRMdgplot3D);
    QLabel *LBLdginfb = new QLabel(tr("Lowest"));
    QLabel *LBLdgsupb = new QLabel(tr("Highest"));

    QHBoxLayout *Layout1 = new QHBoxLayout(FRMdgfilename);
    Layout1->addWidget(TXTdgfilename);

    QHBoxLayout *Layout2 = new QHBoxLayout(FRMdglmaxexp);
    Layout2->addWidget(SPBdglmaxexp);

    QHBoxLayout *Layout3 = new QHBoxLayout();
    Layout3->addWidget(LBLdglongthreshold);
    Layout3->addWidget(SPBdglongthreshold);
    Layout3->setAlignment(Qt::AlignCenter);

    QVBoxLayout *Layout4 = new QVBoxLayout(FRMdglong);
    Layout4->addLayout(Layout3);

    QHBoxLayout *Layout5 = new QHBoxLayout();
    Layout5->addWidget(CHKdgextralines);
    Layout5->setAlignment(Qt::AlignCenter);

    QHBoxLayout *Layout6 = new QHBoxLayout();
    Layout6->addWidget(LBLdgnumpnt);
    Layout6->addWidget(TXTdgnumpnt);

    QHBoxLayout *Layout7 = new QHBoxLayout();
    Layout7->addWidget(LBLdgnlinpernuc);
    Layout7->addWidget(TXTdgnlinpernuc);

    QHBoxLayout *Layout8 = new QHBoxLayout();
    Layout8->addWidget(LBLdgdirectionset);
    Layout8->addWidget(CMBdgdirectionset);

    QHBoxLayout *Layout9 = new QHBoxLayout();
    Layout9->addWidget(LBLdgdlt0);
    Layout9->addWidget(TXTdgdlt0);

    QHBoxLayout *Layout10 = new QHBoxLayout(FRMdgplottype);
    Layout10->addWidget(RBTdg2D);
    Layout10->addWidget(RBTdg3D);

    QGridLayout *Layout11 = new QGridLayout(FRMdgplaneABC);
    Layout11->addWidget(LBLdgplaneeq,0,0,1,2,Qt::AlignCenter);
    Layout11->addWidget(LBLdgplaneA,1,0);
    Layout11->addWidget(TXTdgplaneA,1,1);
    Layout11->addWidget(LBLdgplaneB,2,0);
    Layout11->addWidget(TXTdgplaneB,2,1);
    Layout11->addWidget(LBLdgplaneC,3,0);
    Layout11->addWidget(TXTdgplaneC,3,1);
    Layout11->addWidget(LBLdguvratiohead,4,0,1,2,Qt::AlignCenter);
    Layout11->addWidget(LBLdguvratio,5,0);
    Layout11->addWidget(TXTdguvratio,5,1);

    QGridLayout *Layout12 = new QGridLayout(FRMdgplane2D);
    Layout12->addWidget(RBTdgplaneXY,0,0,Qt::AlignLeft);
    Layout12->addWidget(RBTdgplaneXZ,0,1,Qt::AlignLeft);
    Layout12->addWidget(RBTdgplaneYZ,1,0,Qt::AlignLeft);
    Layout12->addWidget(RBTdgplaneABC,1,1,Qt::AlignLeft);
    Layout12->addWidget(FRMdgplaneABC,2,0,1,2,Qt::AlignCenter);

    QGridLayout *Layout13 = new QGridLayout();
    Layout13->addWidget(LBLdginfa,1,1,Qt::AlignCenter);
    Layout13->addWidget(LBLdgsupa,1,2,Qt::AlignCenter);
    Layout13->addWidget(LBLdgu,2,0);
    Layout13->addWidget(TXTdguinf,2,1);
    Layout13->addWidget(TXTdgusup,2,2);
    Layout13->addWidget(LBLdgv,3,0);
    Layout13->addWidget(TXTdgvinf,3,1);
    Layout13->addWidget(TXTdgvsup,3,2);
    Layout13->addWidget(FRMdgplane2D,4,0,1,3,Qt::AlignCenter);

    QVBoxLayout *Layout14 = new QVBoxLayout(FRMdgplot2D);
    Layout14->addLayout(Layout7);
    Layout14->addLayout(Layout13);

    QGridLayout *Layout15 = new QGridLayout();
    Layout15->addWidget(LBLdginfb,1,1,Qt::AlignCenter);
    Layout15->addWidget(LBLdgsupb,1,2,Qt::AlignCenter);
    Layout15->addWidget(LBLdgx,2,0);
    Layout15->addWidget(TXTdgxinf,2,1);
    Layout15->addWidget(TXTdgxsup,2,2);
    Layout15->addWidget(LBLdgy,3,0);
    Layout15->addWidget(TXTdgyinf,3,1);
    Layout15->addWidget(TXTdgysup,3,2);
    Layout15->addWidget(LBLdgz,4,0);
    Layout15->addWidget(TXTdgzinf,4,1);
    Layout15->addWidget(TXTdgzsup,4,2);

    QVBoxLayout *Layout16 = new QVBoxLayout(FRMdgplot3D);
    Layout16->addLayout(Layout8);
    Layout16->addLayout(Layout15);

    QVBoxLayout *Layout17 = new QVBoxLayout(FRMdgextralines);
    Layout17->addWidget(CHKdgxyz,0,Qt::AlignCenter);
    Layout17->addWidget(Wtabledendg,0,Qt::AlignCenter);
    Layout17->addWidget(Wtable2dg,0,Qt::AlignCenter);

    QHBoxLayout *Layout18=new QHBoxLayout();
    Layout18->addWidget(LBLdgfilelines);

    QHBoxLayout *Layout19 = new QHBoxLayout();
    Layout19->addWidget(TXTdgfilelines);
    Layout19->addWidget(BTNdgfilelines);

    QVBoxLayout *Layout20 = new QVBoxLayout(FRMdglines);
    Layout20->addLayout(Layout6);
    Layout20->addLayout(Layout9);
    Layout20->addWidget(FRMdgplottype);
    Layout20->addWidget(FRMdgplot2D);
    Layout20->addWidget(FRMdgplot3D);
    Layout20->addLayout(Layout5);
    Layout20->addWidget(FRMdgextralines);
    Layout20->addLayout(Layout18);
    Layout20->addLayout(Layout19);

    QVBoxLayout *Layout21 = new QVBoxLayout(FRMdginput);
    Layout21->addWidget(CHKdginput);

    QHBoxLayout *Layout22 = new QHBoxLayout(FRMdgmpi);
    Layout22->addWidget(CHKdgmpi);
    Layout22->addWidget(LBLdgmpi);
    Layout22->addWidget(SPBdgmpi);

    QHBoxLayout *Layout23 = new QHBoxLayout();
    Layout23->addWidget(BTNtexto[11],0,Qt::AlignLeft);
    Layout23->addWidget(BTNstop[11],0,Qt::AlignRight);
    Layout23->addWidget(BTNexecDamdengrad,0,Qt::AlignRight);

    QVBoxLayout *page_densgradLayout = new QVBoxLayout(page_densgrad);
    page_densgradLayout->addWidget(FRMdgfilename);
    page_densgradLayout->addWidget(FRMdglmaxexp);
    page_densgradLayout->addWidget(FRMdglong);
    page_densgradLayout->addWidget(FRMdglines);
    page_densgradLayout->addWidget(FRMdginput);
    page_densgradLayout->addWidget(FRMdgmpi);
    page_densgradLayout->addLayout(Layout23);
    page_densgradLayout->addStretch();
}

//    page_frad: RADIAL FACTORS
//    =========================

void MainWindow::page_frad_widgets()
{
    FRMfraddamfilename = new QGroupBox(tr("Output files prefix"),page_frad);
    FRMfraddamfilename->setMaximumSize(QSize(400, 2000));
    TXTfraddamfilename = new QLineEdit(FRMfraddamfilename);
    connect(TXTfraddamfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

    FRMfradios = new QGroupBox(tr("Tabulation points"),page_frad);
    FRMfradios->setMaximumSize(QSize(400, 2000));
    LBLfradrini = new QLabel(tr("Initial"),FRMfradios);
    TXTfradrini = new QLineEdit(FRMfradios);
    TXTfradrini->setValidator(myDoubleValidator);
    connect(TXTfradrini, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    LBLfradrfin = new QLabel(tr("Final"),FRMfradios);
    TXTfradrfin = new QLineEdit(FRMfradios);
    TXTfradrfin->setValidator(myDoubleValidator);
    connect(TXTfradrfin, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    LBLfradltr = new QLabel(tr("Step"),FRMfradios);
    TXTfradltr = new QLineEdit(FRMfradios);
    TXTfradltr->setValidator(myDoubleValidator);
    connect(TXTfradltr, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    CHKfradextras = new QCheckBox(tr("Extra values"),FRMfradios);
    connect(CHKfradextras, SIGNAL(stateChanged(int)), this, SLOT(CHKfradextras_changed()));

    Wtable5 = new QWidget(FRMfradios);
    SHTfradrlist = new Sheet(5, 1, 0,true, Wtable5);
    QStringList QSLrlist;
    QSLrlist << "R";
    SHTfradrlist->setHeader(QSLrlist);
    Wtable5->setVisible(false);
    Wtable5->setMaximumWidth(160);

    FRMfradialfactors = new QGroupBox(tr("Radial factors"),page_frad);
    FRMfradialfactors->setMaximumSize(QSize(400, 2000));

    LBLfradltab = new QLabel(tr("l"),FRMfradialfactors);

    SPBfradltab = new QSpinBox(FRMfradialfactors);
    SPBfradltab->setRange(0, 1);
    SPBfradltab->setValue(0);
    SPBfradltab->setMaximumWidth(60);
    connect(SPBfradltab, SIGNAL(valueChanged(int)), this, SLOT(SPBfradltab_changed()));

    LBLfradmtab = new QLabel(tr("m"),FRMfradialfactors);
    SPBfradmtab=new QSpinBox(FRMfradialfactors);
    SPBfradmtab->setRange(0, 1);
    SPBfradmtab->setValue(0);
    SPBfradmtab->setMaximumWidth(60);
    SPBfradmtab->setMinimumWidth(60);
    connect(SPBfradmtab, SIGNAL(valueChanged(int)), this, SLOT(SPBfradmtab_changed()));

    BTNexecDamfrad=new QPushButton(QIcon(":/images/exec.png"),tr("Exec"),page_frad);
    BTNexecDamfrad->setToolTip(tr("Tabulate radial factors"));
    connect(BTNexecDamfrad, SIGNAL(clicked()), this, SLOT(execDamfrad()));

    if (!BTNstop[5])
        BTNstop[5]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_frad);
    BTNstop[5]->setToolTip(tr("Kill the process"));
    connect(BTNstop[5], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[5])
        BTNtexto[5]=new QPushButton(QIcon(":/images/document_text.png"),"",page_frad);
    BTNtexto[5]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[5], SIGNAL(clicked()), this, SLOT(importOUT()));

    FRMfradderivadas = new QGroupBox(tr("Derivatives"),page_frad);
    FRMfradderivadas->setMaximumSize(QSize(400, 2000));

    CHKfradderiv1 = new QCheckBox(tr("First derivatives"),FRMfradios);
    CHKfradderiv2 = new QCheckBox(tr("Second derivatives"),FRMfradios);

    FRMfradatoms = new QGroupBox(tr("Centers"),page_frad);
    FRMfradatoms->setMaximumSize(QSize(400, 2000));

    LBLfradatoms = new QLabel("1,3-5,10,13-17,...");
    TXTfradatoms = new QLineEdit(FRMfradatoms);
    QRegExp fradrx("[1-9][-,\\d]*");
    fradvalidator = new QRegExpValidator(fradrx, nullpointer);
    fradlist = new QStringList();
    TXTfradatoms->setValidator(fradvalidator);
    TXTfradatoms->setEnabled(true);
    connect(TXTfradatoms, SIGNAL(textChanged(const QString &)), this, SLOT(TXTfradatoms_changed()));
}

void MainWindow::page_frad_layouts()
{
    QHBoxLayout *Layout0=new QHBoxLayout(FRMfraddamfilename);
    Layout0->addWidget(TXTfraddamfilename);

    QHBoxLayout *Layout1=new QHBoxLayout();
    Layout1->addWidget(BTNtexto[5],0,Qt::AlignLeft);
    Layout1->addWidget(BTNstop[5],0,Qt::AlignRight);
    Layout1->addWidget(BTNexecDamfrad,0,Qt::AlignRight);

    QGridLayout *Layout2=new QGridLayout();
    Layout2->addWidget(LBLfradrini,1,0);
    Layout2->addWidget(LBLfradrfin,1,1);
    Layout2->addWidget(LBLfradltr,1,2);
    Layout2->addWidget(TXTfradrini,2,0);
    Layout2->addWidget(TXTfradrfin,2,1);
    Layout2->addWidget(TXTfradltr,2,2);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addStretch(10);
    Layout3->addWidget(Wtable5);
    Layout3->addStretch(10);

    QVBoxLayout *Layout4 = new QVBoxLayout(FRMfradios);
    Layout4->addLayout(Layout2);
    Layout4->addWidget(CHKfradextras);
    Layout4->addLayout(Layout3);


    QHBoxLayout *Layout5=new QHBoxLayout();
    Layout5->addStretch(10);
    Layout5->addWidget(LBLfradltab,Qt::AlignRight);
    Layout5->addWidget(SPBfradltab,Qt::AlignLeft);
    Layout5->addWidget(LBLfradmtab,Qt::AlignRight);
    Layout5->addWidget(SPBfradmtab,Qt::AlignLeft);
    Layout5->addStretch(10);

    QVBoxLayout *Layout6 = new QVBoxLayout(FRMfradatoms);
    Layout6->addWidget(LBLfradatoms);
    Layout6->addWidget(TXTfradatoms);

    QVBoxLayout *Layout7=new QVBoxLayout(FRMfradialfactors);
    Layout7->addLayout(Layout5);

    QVBoxLayout *Layout8=new QVBoxLayout();
    Layout8->addWidget(CHKfradderiv1);
    Layout8->addWidget(CHKfradderiv2);

    QVBoxLayout *Layout9=new QVBoxLayout(FRMfradderivadas);
    Layout9->addLayout(Layout8);

    QVBoxLayout *page_fradLayout = new QVBoxLayout(page_frad);
    page_fradLayout->addWidget(FRMfraddamfilename);
    page_fradLayout->addWidget(FRMfradios);
    page_fradLayout->addWidget(FRMfradialfactors);
    page_fradLayout->addWidget(FRMfradatoms);
    page_fradLayout->addWidget(FRMfradderivadas);
    page_fradLayout->addLayout(Layout1);
    page_fradLayout->addStretch();
}

//    page_orimult: ORIENTED MULTIPOLES
//    =================================

void MainWindow::page_orimult_widgets()
{
    FRMmrotorimultfilename = new QGroupBox(tr("Output files prefix"),page_orimult);
    FRMmrotorimultfilename->setMaximumSize(QSize(400, 2000));

    TXTmrotorimultfilename = new QLineEdit(FRMmrotorimultfilename);
    connect(TXTmrotorimultfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMmrotmultipoles = new QGroupBox(tr("Multipoles"),page_orimult);
    FRMmrotmultipoles->setMaximumSize(QSize(400, 2000));

    LBLmrotlmin = new QLabel(tr("lmin"),FRMmrotmultipoles);

    SPBmrotlmin = new QSpinBox(FRMmrotmultipoles);
    SPBmrotlmin->setRange(0, 25);
    SPBmrotlmin->setValue(0);
    SPBmrotlmin->setMaximumWidth(60);
    connect(SPBmrotlmin, SIGNAL(valueChanged(int)), this, SLOT(SPBmrotlmin_changed()));

    LBLmrotlmax = new QLabel(tr("lmax"),FRMmrotmultipoles);
    SPBmrotlmax=new QSpinBox(FRMmrotmultipoles);
    SPBmrotlmax->setRange(0, 25);
    SPBmrotlmax->setValue(0);
    SPBmrotlmax->setMaximumWidth(60);
    connect(SPBmrotlmax, SIGNAL(valueChanged(int)), this, SLOT(SPBmrotlmax_changed()));

    FRMmrotcenters = new QGroupBox(tr("Centers defining plane"),page_orimult);
    FRMmrotcenters->setMaximumSize(QSize(400, 2000));

    LBLmrotleft = new QLabel(tr("Left"),FRMmrotcenters);

    SPBmrotleft = new QSpinBox(FRMmrotcenters);
    SPBmrotleft->setRange(1, 1);
    SPBmrotleft->setValue(1);
    spbmrotleft = 0;
    SPBmrotleft->setMaximumWidth(60);
    SPBmrotleft->setEnabled(true);
    spbmrotleft = 0;    // Used to store the previous value of SPBmrotleft (see SLOT SPBmrotleft_changed)
    connect(SPBmrotleft, SIGNAL(valueChanged(int)), this, SLOT(SPBmrotleft_changed()));

    LBLmrotmiddle = new QLabel(tr("Middle"),FRMmrotcenters);
    SPBmrotmiddle = new QSpinBox(FRMmrotcenters);
    SPBmrotmiddle->setRange(1, 1);
    SPBmrotmiddle->setValue(1);
    spbmrotmiddle = 0;
    SPBmrotmiddle->setMaximumWidth(60);
    SPBmrotmiddle->setEnabled(true);
    spbmrotmiddle = 0;    // Used to store the previous value of SPBmrotmiddle (see SLOT SPBmrotmiddle_changed)
    connect(SPBmrotmiddle, SIGNAL(valueChanged(int)), this, SLOT(SPBmrotmiddle_changed()));

    LBLmrotright = new QLabel(tr("Right"),FRMmrotcenters);
    SPBmrotright = new QSpinBox(FRMmrotcenters);
    SPBmrotright->setRange(1, 1);
    SPBmrotright->setValue(1);
    spbmrotright = 0;
    SPBmrotright->setMaximumWidth(60);
    SPBmrotright->setEnabled(true);
    spbmrotright = 0;    // Used to store the previous value of SPBmrotright (see SLOT SPBmrotright_changed)
    connect(SPBmrotright, SIGNAL(valueChanged(int)), this, SLOT(SPBmrotright_changed()));

    BTNexecDammultrot=new QPushButton(QIcon(":/images/exec.png"),tr("Exec"),page_orimult);
    BTNexecDammultrot->setToolTip(tr("Orient multipoles in a frame with the Z axis "
            "orthogonal to the plane defined by the three centers"));
    connect(BTNexecDammultrot, SIGNAL(clicked()), this, SLOT(execDammultrot()));

    if (!BTNstop[6])
        BTNstop[6]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_orimult);
    BTNstop[6]->setToolTip(tr("Kill the process"));
    connect(BTNstop[6], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[6])
        BTNtexto[6]=new QPushButton(QIcon(":/images/document_text.png"),"",page_orimult);
    BTNtexto[6]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[6], SIGNAL(clicked()), this, SLOT(importOUT()));

    FRMmrotatoms = new QGroupBox(tr("Atomic fragments"),page_orimult);
    FRMmrotatoms->setMaximumSize(QSize(400, 2000));

    LBLmrotorimultatoms = new QLabel("1,3-5,10,13-17,...");

    TXTmrotorimultatoms = new QLineEdit(FRMmrotatoms);
    QRegExp orirx("[1-9][-,\\d]*");
    mrotorimultvalidator = new QRegExpValidator(orirx, nullpointer);
    mrotorimultlist = new QStringList();
    TXTmrotorimultatoms->setValidator(mrotorimultvalidator);
    TXTmrotorimultatoms->setEnabled(true);
    connect(TXTmrotorimultatoms, SIGNAL(textChanged(const QString &)), this, SLOT(TXTmrotorimultatoms_changed()));
}

void MainWindow::page_orimult_layouts()
{
    QHBoxLayout *Layout0=new QHBoxLayout(FRMmrotorimultfilename);
    Layout0->addWidget(TXTmrotorimultfilename);

    QGridLayout *Layout1=new QGridLayout(FRMmrotmultipoles);
    Layout1->addWidget(LBLmrotlmin,0,0);
    Layout1->addWidget(SPBmrotlmin,0,1);
    Layout1->addWidget(LBLmrotlmax,1,0);
    Layout1->addWidget(SPBmrotlmax,1,1);

    QGridLayout *Layout2=new QGridLayout(FRMmrotcenters);
    Layout2->addWidget(LBLmrotleft,0,0);
    Layout2->addWidget(SPBmrotleft,1,0);
    Layout2->addWidget(LBLmrotmiddle,0,1);
    Layout2->addWidget(SPBmrotmiddle,1,1);
    Layout2->addWidget(LBLmrotright,0,2);
    Layout2->addWidget(SPBmrotright,1,2);

    QVBoxLayout *Layout3=new QVBoxLayout(FRMmrotatoms);
    Layout3->addWidget(LBLmrotorimultatoms);
    Layout3->addWidget(TXTmrotorimultatoms);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(BTNtexto[6],0,Qt::AlignLeft);
    Layout4->addWidget(BTNstop[6],0,Qt::AlignRight);
    Layout4->addWidget(BTNexecDammultrot,0,Qt::AlignRight);

    QVBoxLayout *Layout5=new QVBoxLayout(page_orimult);
    Layout5->addWidget(FRMmrotorimultfilename);
    Layout5->addWidget(FRMmrotmultipoles);
    Layout5->addWidget(FRMmrotcenters);
    Layout5->addWidget(FRMmrotatoms);
    Layout5->addLayout(Layout4);
    Layout5->addStretch();
}

//    page_MO: MOLECULAR ORBITALS
//    ===========================

void MainWindow::page_MO_widgets()
{
//          Import File

    FRMMOImportfile = new QGroupBox(tr("Molecular orbitals"),page_MO);
    FRMMOImportfile->setMaximumSize(QSize(400, 2000));

    LBLMOImportFile = new QLabel(tr("Import data from")+":");

    TXTMOImportfile = new QLineEdit();

    BTNMOImportFile = new QToolButton();
    BTNMOImportFile->setText(tr("..."));
    connect(BTNMOImportFile, SIGNAL(clicked()), this, SLOT(ImportFileNameMO()));

    FRMMOfilename = new QGroupBox(tr("Output files prefix"),page_MO);
    FRMMOfilename->setMaximumSize(QSize(400, 2000));

    TXTMOfilename = new QLineEdit(FRMMOfilename);

    connect(TXTMOfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMMOchoose = new QGroupBox(tr("Molecular orbitals"),page_MO);
    FRMMOchoose->setMaximumSize(QSize(400, 2000));

    LBLMOchoose = new QLabel("1,3-5,10,13-17,...");

    TXTMOchoose = new QLineEdit(FRMMOchoose);
    QRegExp rx("[1-9][-,\\d]*");
    MOvalidator = new QRegExpValidator(rx, nullpointer);
    MOlist = new QStringList();
    TXTMOchoose->setValidator(MOvalidator);
    connect(TXTMOchoose, SIGNAL(textChanged(const QString &)), this, SLOT(TXTMOchoose_changed()));

//            Derivatives

    FRMMOderivs = new QGroupBox(tr("Derivatives"),page_MO);
    FRMMOderivs->setMaximumSize(QSize(400, 2000));

    CHKMOgrad = new QCheckBox(tr("Gradient"),FRMMOderivs);
    CHKMOgrad->setToolTip(tr("Tabulate MO gradient"));
    CHKMOgrad->setChecked(true);

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

//            Grid

    FRMMOgrid = new QGroupBox(tr("Grid"),page_MO);
    FRMMOgrid->setMaximumSize(QSize(400, 2000));

    CHKMOgrid=new QCheckBox(tr("Generate grid"),FRMMOgrid);
    CHKMOgrid->setChecked(true);
    connect(CHKMOgrid, SIGNAL(stateChanged(int)), this, SLOT(CHKMOgrid_changed()));

    FRMMOgridtype = new QGroupBox(tr("Grid type"));

//            2D Grid

    RBTMO2D=new QRadioButton(tr("2D grid"),FRMMOgridtype);
    RBTMO2D->setChecked(false);
    connect(RBTMO2D, SIGNAL(toggled (bool)), this, SLOT(RBTMO2D3D_changed()));

    FRMMOgrid2D = new QGroupBox(tr("2D grid"));
    FRMMOgrid2D->setMaximumSize(QSize(400, 2000));
    FRMMOgrid2D->setVisible(false);

    LBLMOu=new QLabel(tr("u"),FRMMOgrid2D);
    LBLMOv=new QLabel(tr("v"),FRMMOgrid2D);
    LBLMOinf2d3=new QLabel(tr("Lowest"),FRMMOgrid2D);
    LBLMOsup2d=new QLabel(tr("Highest"),FRMMOgrid2D);

    TXTMOuinf=new QLineEdit(FRMMOgrid2D);
    TXTMOuinf->setText("-4.0");
    TXTMOuinf->setAlignment(Qt::AlignRight);
    TXTMOuinf->setValidator(myDoubleValidator);
    connect(TXTMOuinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOusup=new QLineEdit(FRMMOgrid2D);
    TXTMOusup->setText("4.0");
    TXTMOusup->setAlignment(Qt::AlignRight);
    TXTMOusup->setValidator(myDoubleValidator);
    connect(TXTMOusup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOvinf=new QLineEdit(FRMMOgrid2D);
    TXTMOvinf->setText("-4.0");
    TXTMOvinf->setAlignment(Qt::AlignRight);
    TXTMOvinf->setValidator(myDoubleValidator);
    connect(TXTMOvinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOvsup=new QLineEdit(FRMMOgrid2D);
    TXTMOvsup->setText("4.0");
    TXTMOvsup->setAlignment(Qt::AlignRight);
    TXTMOvsup->setValidator(myDoubleValidator);
    connect(TXTMOvsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMMOsurfpar = new QGroupBox(tr("Parametric equations"),FRMdensgrid2D);
    FRMMOsurfpar->setMaximumSize(QSize(400, 2000));
    FRMMOsurfpar->setVisible(false);

    LBLMOxformula2D=new QLabel(tr("x(u,v) = "),FRMMOgrid2D);
    LBLMOyformula2D=new QLabel(tr("y(u,v) = "),FRMMOgrid2D);
    LBLMOzformula2D=new QLabel(tr("z(u,v) = "),FRMMOgrid2D);

    TXTMOxformula2D=new QLineEdit(FRMMOgrid2D);
    TXTMOxformula2D->setText("u");

    TXTMOyformula2D=new QLineEdit(FRMMOgrid2D);
    TXTMOyformula2D->setText("v");

    TXTMOzformula2D=new QLineEdit(FRMMOgrid2D);
    TXTMOzformula2D->setText("0");

    FRMMOsurftype = new QGroupBox(tr("Surface type"),FRMMOgrid2D);
    FRMMOsurftype->setMaximumSize(QSize(400, 2000));
    RBTMOplane = new QRadioButton(tr("Plane"),FRMMOsurftype);
    RBTMOplane->setChecked(true);
    connect(RBTMOplane, SIGNAL(toggled (bool)), this, SLOT(RBTMOplane_changed()));

    RBTMOothersurf = new QRadioButton(tr("Parametric surface"),FRMMOsurftype);
    RBTMOothersurf->setToolTip(tr("Surface equation supplied in parametric form: x=x(u,v), y=y(u,v), z=z(u,v)"));
    RBTMOothersurf->setChecked(false);

    FRMMOplane2D = new QGroupBox(tr("2D Planes"),FRMMOgrid2D);
    FRMMOplane2D->setMaximumSize(QSize(400, 2000));
    FRMMOplane2D->setVisible(true);

    RBTMOplaneXY = new QRadioButton(tr("XY "),FRMMOgrid2D);
    RBTMOplaneXY->setChecked(true);
    connect(RBTMOplaneXY, SIGNAL(toggled (bool)), this, SLOT(RBTMO2Dplanes_changed()));

    RBTMOplaneXZ = new QRadioButton(tr("XZ "),FRMMOgrid2D);
    RBTMOplaneXZ->setChecked(false);
    connect(RBTMOplaneXZ, SIGNAL(toggled (bool)), this, SLOT(RBTMO2Dplanes_changed()));

    RBTMOplaneYZ = new QRadioButton(tr("YZ "),FRMMOgrid2D);
    RBTMOplaneYZ->setChecked(false);
    connect(RBTMOplaneYZ, SIGNAL(toggled (bool)), this, SLOT(RBTMO2Dplanes_changed()));

    RBTMOplaneABC = new QRadioButton(tr("Other "),FRMMOgrid2D);
    RBTMOplaneABC->setChecked(false);
    connect(RBTMOplaneABC, SIGNAL(toggled (bool)), this, SLOT(RBTMO2Dplanes_changed()));

    FRMMOplaneABC = new QGroupBox(tr("Plane parameters"),FRMMOgrid2D);
    FRMMOplaneABC->setMaximumSize(QSize(400, 2000));
    FRMMOplaneABC->setVisible(false);

    TXTMOplaneA=new QLineEdit(FRMMOplaneABC);
    TXTMOplaneA->setText("0.");
    TXTMOplaneA->setValidator(myDoubleValidator);
    connect(TXTMOplaneA, SIGNAL(editingFinished()), this, SLOT(RBTMO2Dplanes_changed()));

    TXTMOplaneB=new QLineEdit(FRMMOplaneABC);
    TXTMOplaneB->setText("0.");
    TXTMOplaneB->setValidator(myDoubleValidator);
    connect(TXTMOplaneB, SIGNAL(editingFinished()), this, SLOT(RBTMO2Dplanes_changed()));

    TXTMOplaneC=new QLineEdit(FRMMOplaneABC);
    TXTMOplaneC->setText("1.");
    TXTMOplaneC->setValidator(myDoubleValidator);
    connect(TXTMOplaneC, SIGNAL(editingFinished()), this, SLOT(RBTMO2Dplanes_changed()));
    MOplanecase = 1;

    FRMMOresol2D = new QGroupBox(tr("Custom resolution"));
    FRMMOresol2D->setHidden(true);

    LBLMOuresol = new QLabel(tr("u"),FRMMOresol2D);
    LBLMOvresol = new QLabel(tr("v"),FRMMOresol2D);

    SPBMOures = new QSpinBox(FRMMOresol2D);
    SPBMOures->setMinimum(4);
    SPBMOures->setMaximum(2049);
    SPBMOures->setValue(129);
    SPBMOures->setSingleStep(10);

    SPBMOvres = new QSpinBox(FRMMOresol2D);
    SPBMOvres->setMinimum(4);
    SPBMOvres->setMaximum(2049);
    SPBMOvres->setValue(129);
    SPBMOvres->setSingleStep(10);

//            3D Grid

    RBTMO3D=new QRadioButton(tr("3D grid"),FRMMOgridtype);
    RBTMO3D->setChecked(true);

    FRMMOgrid3D = new QGroupBox(tr("3D grid"));
    FRMMOgrid3D->setMaximumSize(QSize(400, 2000));
    FRMMOgrid3D->setVisible(true);

    LBLMOx=new QLabel(tr("x"),FRMMOgrid3D);
    LBLMOy=new QLabel(tr("y"),FRMMOgrid3D);
    LBLMOz=new QLabel(tr("z"),FRMMOgrid3D);
    LBLMOinf=new QLabel(tr("Lowest"),FRMMOgrid3D);
    LBLMOsup=new QLabel(tr("Highest"),FRMMOgrid3D);   

    TXTMOxinf=new QLineEdit(FRMpotgrid3D);
    TXTMOxinf->setText("-4.0");
    TXTMOxinf->setAlignment(Qt::AlignRight);
    TXTMOxinf->setValidator(myDoubleValidator);
    connect(TXTMOxinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOxsup=new QLineEdit(FRMpotgrid3D);
    TXTMOxsup->setText("4.0");
    TXTMOxsup->setAlignment(Qt::AlignRight);
    TXTMOxsup->setValidator(myDoubleValidator);
    connect(TXTMOxsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOyinf=new QLineEdit(FRMpotgrid3D);
    TXTMOyinf->setText("-4.0");
    TXTMOyinf->setAlignment(Qt::AlignRight);
    TXTMOyinf->setValidator(myDoubleValidator);
    connect(TXTMOyinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOysup=new QLineEdit(FRMpotgrid3D);
    TXTMOysup->setText("4.0");
    TXTMOysup->setAlignment(Qt::AlignRight);
    TXTMOysup->setValidator(myDoubleValidator);
    connect(TXTMOysup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));
    TXTMOzinf=new QLineEdit(FRMpotgrid3D);
    TXTMOzinf->setText("-4.0");
    TXTMOzinf->setAlignment(Qt::AlignRight);
    TXTMOzinf->setValidator(myDoubleValidator);
    connect(TXTMOzinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTMOzsup=new QLineEdit(FRMpotgrid3D);
    TXTMOzsup->setText("4.0");
    TXTMOzsup->setAlignment(Qt::AlignRight);
    TXTMOzsup->setValidator(myDoubleValidator);
    connect(TXTMOzsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMMOresol3D = new QGroupBox(tr("Custom resolution"));
    FRMMOresol3D->setHidden(true);

    LBLMOxresol = new QLabel(tr("x"),FRMMOresol3D);
    LBLMOyresol = new QLabel(tr("y"),FRMMOresol3D);
    LBLMOzresol = new QLabel(tr("z"),FRMMOresol3D);

    SPBMOxres = new QSpinBox(FRMMOresol3D);
    SPBMOxres->setMinimum(4);
    SPBMOxres->setMaximum(2049);
    SPBMOxres->setValue(129);
    SPBMOxres->setSingleStep(10);

    SPBMOyres = new QSpinBox(FRMMOresol3D);
    SPBMOyres->setMinimum(4);
    SPBMOyres->setMaximum(2049);
    SPBMOyres->setValue(129);
    SPBMOyres->setSingleStep(10);

    SPBMOzres = new QSpinBox(FRMMOresol3D);
    SPBMOzres->setMinimum(4);
    SPBMOzres->setMaximum(2049);
    SPBMOzres->setValue(129);
    SPBMOzres->setSingleStep(10);

//          Resolution

    FRMMOgridres = new QGroupBox(tr("Resolution"));

    RBTMOrlow=new QRadioButton(tr("Low"),FRMMOgridres);
    RBTMOrlow->setChecked(true);
    connect(RBTMOrlow, SIGNAL(toggled (bool)), this, SLOT(MO_resolution_changed()));

    RBTMOrmedium=new QRadioButton(tr("Medium"),FRMMOgridres);
    connect(RBTMOrmedium, SIGNAL(toggled (bool)), this, SLOT(MO_resolution_changed()));

    RBTMOrhigh=new QRadioButton(tr("High"),FRMMOgridres);
    connect(RBTMOrhigh, SIGNAL(toggled (bool)), this, SLOT(MO_resolution_changed()));

    RBTMOrcustom=new QRadioButton(tr("Custom"),FRMMOgridres);
    connect(RBTMOrcustom, SIGNAL(toggled (bool)), this, SLOT(MO_resolution_changed()));

    if (RBTMO3D->isChecked()){
        RBTMOrlow->setToolTip("65x65x65");
        RBTMOrmedium->setToolTip("129x129x129");
        RBTMOrhigh->setToolTip("257x257x257");
}
    else{
        RBTMOrlow->setToolTip("129x129");
        RBTMOrmedium->setToolTip("257x257");
        RBTMOrhigh->setToolTip("513x513");
    }

//          Generate input only

    FRMMOinput = new QGroupBox(tr("Input only"),page_MO);
    FRMMOinput->setMaximumSize(QSize(400, 2000));

    CHKMOinput=new QCheckBox(tr("Generate input file only"),FRMMOinput);
    CHKMOinput->setChecked(false);
    CHKMOinput->setEnabled(true);
    connect(CHKMOinput, SIGNAL(stateChanged(int)), this, SLOT(CHKMOinput_changed(int)));

//            MPI

    FRMMOmpi = new QGroupBox(tr("Parallel computing"),page_MO);
    FRMMOmpi->setMaximumSize(QSize(400, 2000));

    CHKMOmpi = new QCheckBox(tr("MPI"),FRMMOmpi);

    LBLMOmpi = new QLabel(tr("Number of processors"),FRMMOmpi);

    SPBMOmpi = new QSpinBox(FRMMOmpi);
    SPBMOmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBMOmpi->setValue(1);
    SPBMOmpi->setMaximumWidth(60);
    SPBMOmpi->setToolTip(tr("Number of processors"));
    if (mpi){
        CHKMOmpi->setChecked(true);
        SPBMOmpi->setEnabled(true);
    }
    else{
        FRMMOmpi->setHidden(true);
        CHKMOmpi->setChecked(false);
        SPBMOmpi->setEnabled(false);
    }
    connect(CHKMOmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKMOmpi_changed(int)));
    connect(SPBMOmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBMOmpi_changed(int)));

//            List file, Stop, Execute

    BTNexecDamorb=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_MO);
    BTNexecDamorb->setToolTip(tr("Molecular orbitals calculation"));
    connect(BTNexecDamorb, SIGNAL(clicked()), this, SLOT(execDamorb()));

    if (!BTNstop[7])
        BTNstop[7]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_MO);
    BTNstop[7]->setToolTip(tr("Kill the process"));
    connect(BTNstop[7], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[7])
        BTNtexto[7]=new QPushButton(QIcon(":/images/document_text.png"),"",page_MO);
    BTNtexto[7]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[7], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_MO_layouts()
{
    QLabel *LBLMOplaneeq=new QLabel(tr("(A x + B y + C z = 0)"),FRMMOplaneABC);
    QLabel *LBLMOplaneA=new QLabel(tr("A: "),FRMMOplaneABC);
    QLabel *LBLMOplaneB=new QLabel(tr("B: "),FRMMOplaneABC);
    QLabel *LBLMOplaneC=new QLabel(tr("C: "),FRMMOplaneABC);

    QHBoxLayout *Layout0=new QHBoxLayout();
    Layout0->addWidget(LBLMOImportFile);

    QHBoxLayout *Layout1=new QHBoxLayout();
    Layout1->addWidget(TXTMOImportfile);
    Layout1->addWidget(BTNMOImportFile);

    QVBoxLayout *Layout2=new QVBoxLayout(FRMMOImportfile);
    Layout2->addLayout(Layout0);
    Layout2->addLayout(Layout1);

    QHBoxLayout *Layout3=new QHBoxLayout(FRMMOfilename);
    Layout3->addWidget(TXTMOfilename);

    QVBoxLayout *Layout4=new QVBoxLayout(FRMMOchoose);
    Layout4->addWidget(LBLMOchoose,Qt::AlignCenter);
    Layout4->addWidget(TXTMOchoose);

    QVBoxLayout *Layout5=new QVBoxLayout(FRMMOderivs);
    Layout5->addWidget(CHKMOgrad);


    QHBoxLayout *Layout6=new QHBoxLayout();
    Layout6->addWidget(CHKMOgrid,0,Qt::AlignCenter);

    QHBoxLayout *Layout7=new QHBoxLayout(FRMMOgridtype);
    Layout7->addWidget(RBTMO2D);
    Layout7->addWidget(RBTMO3D);

    QGridLayout *Layout8=new QGridLayout(FRMMOgrid3D);
    Layout8->addWidget(LBLMOinf,1,1,Qt::AlignCenter);
    Layout8->addWidget(LBLMOsup,1,2,Qt::AlignCenter);
    Layout8->addWidget(LBLMOx,2,0);
    Layout8->addWidget(TXTMOxinf,2,1);
    Layout8->addWidget(TXTMOxsup,2,2);
    Layout8->addWidget(LBLMOy,3,0);
    Layout8->addWidget(TXTMOyinf,3,1);
    Layout8->addWidget(TXTMOysup,3,2);
    Layout8->addWidget(LBLMOz,4,0);
    Layout8->addWidget(TXTMOzinf,4,1);
    Layout8->addWidget(TXTMOzsup,4,2);

    QGridLayout *Layout9=new QGridLayout();
    Layout9->addWidget(LBLMOinf2d3,1,1,Qt::AlignCenter);
    Layout9->addWidget(LBLMOsup2d,1,2,Qt::AlignCenter);
    Layout9->addWidget(LBLMOu,2,0);
    Layout9->addWidget(TXTMOuinf,2,1);
    Layout9->addWidget(TXTMOusup,2,2);
    Layout9->addWidget(LBLMOv,3,0);
    Layout9->addWidget(TXTMOvinf,3,1);
    Layout9->addWidget(TXTMOvsup,3,2);

    QVBoxLayout *Layout10=new QVBoxLayout(FRMMOsurftype);
    Layout10->addWidget(RBTMOplane);
    Layout10->addWidget(RBTMOothersurf);

    QGridLayout *Layout11=new QGridLayout(FRMMOplaneABC);
    Layout11->addWidget(LBLMOplaneeq,0,0,1,2,Qt::AlignCenter);
    Layout11->addWidget(LBLMOplaneA,1,0);
    Layout11->addWidget(TXTMOplaneA,1,1);
    Layout11->addWidget(LBLMOplaneB,2,0);
    Layout11->addWidget(TXTMOplaneB,2,1);
    Layout11->addWidget(LBLMOplaneC,3,0);
    Layout11->addWidget(TXTMOplaneC,3,1);

    QGridLayout *Layout12=new QGridLayout(FRMMOplane2D);
    Layout12->addWidget(RBTMOplaneXY,0,0,Qt::AlignLeft);
    Layout12->addWidget(RBTMOplaneXZ,0,1,Qt::AlignLeft);
    Layout12->addWidget(RBTMOplaneYZ,1,0,Qt::AlignLeft);
    Layout12->addWidget(RBTMOplaneABC,1,1,Qt::AlignLeft);
    Layout12->addWidget(FRMMOplaneABC,2,0,1,2,Qt::AlignCenter);

    QGridLayout *Layout13=new QGridLayout(FRMMOsurfpar);
    Layout13->addWidget(LBLMOxformula2D,1,0);
    Layout13->addWidget(TXTMOxformula2D,1,1);
    Layout13->addWidget(LBLMOyformula2D,2,0);
    Layout13->addWidget(TXTMOyformula2D,2,1);
    Layout13->addWidget(LBLMOzformula2D,3,0);
    Layout13->addWidget(TXTMOzformula2D,3,1);

    QVBoxLayout *Layout14=new QVBoxLayout(FRMMOgrid2D);
    Layout14->addLayout(Layout9,Qt::AlignCenter);
    Layout14->addWidget(FRMMOsurftype);
    Layout14->addWidget(FRMMOplane2D);
    Layout14->addWidget(FRMMOsurfpar);

    QVBoxLayout *Layout15=new QVBoxLayout();
    Layout15->addWidget(FRMMOgrid3D);
    Layout15->addWidget(FRMMOgrid2D);

    QHBoxLayout *Layout16=new QHBoxLayout();
    Layout16->addWidget(RBTMOrlow);
    Layout16->addWidget(RBTMOrmedium);

    QHBoxLayout *Layout17=new QHBoxLayout();
    Layout17->addWidget(RBTMOrhigh);
    Layout17->addWidget(RBTMOrcustom);

    QGridLayout *Layout18=new QGridLayout(FRMMOresol2D);
    Layout18->addWidget(LBLMOuresol,0,0,Qt::AlignCenter);
    Layout18->addWidget(LBLMOvresol,0,1,Qt::AlignCenter);
    Layout18->addWidget(SPBMOures,1,0);
    Layout18->addWidget(SPBMOvres,1,1);

    QGridLayout *Layout19=new QGridLayout(FRMMOresol3D);
    Layout19->addWidget(LBLMOxresol,0,0,Qt::AlignCenter);
    Layout19->addWidget(LBLMOyresol,0,1,Qt::AlignCenter);
    Layout19->addWidget(LBLMOzresol,0,2,Qt::AlignCenter);
    Layout19->addWidget(SPBMOxres,1,0);
    Layout19->addWidget(SPBMOyres,1,1);
    Layout19->addWidget(SPBMOzres,1,2);

    QVBoxLayout *Layout20=new QVBoxLayout(FRMMOgridres);
    Layout20->addLayout(Layout16);
    Layout20->addLayout(Layout17);
    Layout20->addWidget(FRMMOresol2D);
    Layout20->addWidget(FRMMOresol3D);

    QVBoxLayout *Layout21=new QVBoxLayout(FRMMOgrid);
    Layout21->addLayout(Layout6);
    Layout21->addWidget(FRMMOgridtype);
    Layout21->addLayout(Layout15);
    Layout21->addWidget(FRMMOgridres);

    QVBoxLayout *Layout22 = new QVBoxLayout(FRMMOinput);
    Layout22->addWidget(CHKMOinput);

    QHBoxLayout *Layout23 = new QHBoxLayout(FRMMOmpi);
    Layout23->addWidget(CHKMOmpi);
    Layout23->addWidget(LBLMOmpi);
    Layout23->addWidget(SPBMOmpi);

    QHBoxLayout *Layout24 = new QHBoxLayout();
    Layout24->addWidget(BTNtexto[7],0,Qt::AlignLeft);
    Layout24->addWidget(BTNstop[7],0,Qt::AlignRight);
    Layout24->addWidget(BTNexecDamorb,0,Qt::AlignRight);

    QVBoxLayout *page_MOLayout = new QVBoxLayout(page_MO);
    page_MOLayout->addWidget(FRMMOImportfile);
    page_MOLayout->addWidget(FRMMOfilename);
    page_MOLayout->addWidget(FRMMOchoose);
    page_MOLayout->addWidget(FRMMOderivs);
    page_MOLayout->addWidget(FRMMOgrid);
    page_MOLayout->addWidget(FRMMOinput);
    page_MOLayout->addWidget(FRMMOmpi);
    page_MOLayout->addLayout(Layout24);
    page_MOLayout->addStretch();
}
//    page_ZJdens: ZERNIKE-JACOBI EXPANSIONS OF DENSITY
//    =================================================

void MainWindow::page_ZJdens_widgets()
{
    FRMZJrstar = new QGroupBox(tr("Ball size (in bohr)"),page_ZJdens);
    FRMZJrstar->setMaximumSize(QSize(400, 2000));

    QDoubleValidator *ZJrstarValidator = new QDoubleValidator(nullpointer);
    ZJrstarValidator->setLocale(QLocale::English);
    ZJrstarValidator->setBottom(0.0);
    TXTZJrstar=new QLineEdit(FRMZJrstar);
    TXTZJrstar->setText("10.0");
    TXTZJrstar->setAlignment(Qt::AlignRight);
    TXTZJrstar->setValidator(ZJrstarValidator);
    connect(TXTZJrstar, SIGNAL(textChanged(const QString &)), this, SLOT(TXTZJrstar_changed()));
    rstar = TXTZJrstar->text().toDouble();

    FRMZJrstartype = new QGroupBox(tr(""),FRMZJrstar);
    FRMZJrstartype->setMaximumSize(QSize(350, 2000));

    RBTZJrstarabs=new QRadioButton(tr("Absolute"),FRMZJrstartype);
    RBTZJrstarabs->setChecked(true);
    RBTZJrstarabs->setToolTip(tr("Value of the ball radius"));
    connect(RBTZJrstarabs, SIGNAL(toggled (bool)), this, SLOT(TXTZJrstar_changed()));

    RBTZJrstarrel=new QRadioButton(tr("Relative"),FRMZJrstartype);
    RBTZJrstarrel->setToolTip(tr("Value added to the largest distance of nuclei to origin to define the ball radius"));
    connect(RBTZJrstarrel, SIGNAL(toggled (bool)), this, SLOT(TXTZJrstar_changed()));

//            Expansion length

    FRMZJlength = new QGroupBox(tr("Expansion length"),page_ZJdens);
    FRMZJlength->setMaximumSize(QSize(400, 2000));

//            Highest l in expansion

    LBLZJlmax=new QLabel(tr("Highest l"),FRMZJlength);

    SPBZJlmax=new QSpinBox(FRMZJlength);
    SPBZJlmax->setRange(0, 22);
    SPBZJlmax->setValue(10);
    SPBZJlmax->setMaximumWidth(60);
    SPBZJlmax->setToolTip(tr("Size of expansion"));
    connect(SPBZJlmax, SIGNAL(valueChanged(int)), this, SLOT(SPBZJlmax_changed()));

//            Highest k in expansion

    LBLZJkmax=new QLabel(tr("Highest k"),FRMZJlength);

    SPBZJkmax=new QSpinBox(FRMZJlength);
    SPBZJkmax->setRange(0, MAX_KEXPZJ);
    SPBZJkmax->setValue(10);
    SPBZJkmax->setMaximumWidth(60);
    SPBZJkmax->setToolTip(tr("Size of expansion"));
    connect(SPBZJkmax, SIGNAL(valueChanged(int)), this, SLOT(SPBZJkmax_changed()));

    RBTZJechelon=new QRadioButton(tr("Echelon type"),FRMZJlength);
    RBTZJechelon->setToolTip(tr("Number of functions per l equal to max(highest l+1, highest k)-l"));
    RBTZJechelon->setChecked(false);

//            Type of fitting

    FRMZJtype = new QGroupBox(tr("Type of fitting"),page_ZJdens);
    FRMZJtype->setMaximumSize(QSize(400, 2000));

    RBTZJZernike=new QRadioButton(tr("Zernike 3D"),FRMZJtype);
    RBTZJZernike->setChecked(true);
    connect(RBTZJZernike, SIGNAL(toggled (bool)), this, SLOT(TXTValidate_changed()));

    RBTZJacobi=new QRadioButton(tr("Jacobi"),FRMZJtype);
    connect(RBTZJacobi, SIGNAL(toggled (bool)), this, SLOT(TXTValidate_changed()));

//            Number of quadrature points

    FRMZJnquad = new QGroupBox(tr("Quadrature length"),page_ZJdens);
    FRMZJnquad->setMaximumSize(QSize(400, 2000));

    LBLZJnquad=new QLabel(tr("No. of sampling points"),FRMZJnquad);

    SPBZJnquad=new QSpinBox(FRMZJnquad);
    SPBZJnquad->setMinimum(128);
    SPBZJnquad->setMaximum(8192);
    SPBZJnquad->setValue(256);
    SPBZJnquad->setSingleStep(128);
    SPBZJnquad->setMaximumWidth(100);
    SPBZJnquad->setToolTip(tr("Number of points for sampling radial factors in projection"));

//            Thresholds

    FRMZJthreshold = new QGroupBox(tr("Thresholds"),page_ZJdens);
    FRMZJthreshold->setMaximumSize(QSize(400, 2000));

    LBLZJthrmult=new QLabel(tr("Multipole cutoff:")+" 10^",FRMZJthreshold);

    SPBZJthrmult=new QSpinBox(FRMZJthreshold);
    SPBZJthrmult->setRange(-15, 0);
    SPBZJthrmult->setValue(-10);
    SPBZJthrmult->setMaximumWidth(65);
    SPBZJthrmult->setToolTip(tr("Cutoff for multipole printing"));

    LBLZJthrdist=new QLabel(tr("Distributions cutoff:") + " 10^",FRMZJthreshold);

    SPBZJthrdist=new QSpinBox(FRMZJthreshold);
    SPBZJthrdist->setRange(-15, 0);
    SPBZJthrdist->setValue(-12);
    SPBZJthrdist->setMaximumWidth(65);
    SPBZJthrdist->setToolTip(tr("Overlap cutoff for neglecting distributions"));

//                      Generate input only

    FRMZJinput = new QGroupBox(tr("Input only"),page_ZJdens);
    FRMZJinput->setMaximumSize(QSize(400, 2000));

    CHKZJinput=new QCheckBox(tr("Generate input file only"),FRMZJinput);
    CHKZJinput->setChecked(false);
    CHKZJinput->setEnabled(true);
    connect(CHKZJinput, SIGNAL(stateChanged(int)), this, SLOT(CHKZJinput_changed(int)));

//            MPI

    FRMZJmpi = new QGroupBox(tr("Parallel computing"),page_ZJdens);
    FRMZJmpi->setMaximumSize(QSize(400, 2000));

    CHKZJmpi = new QCheckBox(tr("MPI"),FRMZJmpi);

    LBLZJmpi = new QLabel(tr("Number of processors"),FRMZJmpi);

    SPBZJmpi = new QSpinBox(FRMZJmpi);
    SPBZJmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBZJmpi->setValue(1);
    SPBZJmpi->setMaximumWidth(60);
    SPBZJmpi->setToolTip(tr("Number of processors"));
    if (mpi){
        FRMZJmpi->setVisible(true);
        FRMZJmpi->setEnabled(true);
        CHKZJmpi->setChecked(true);
        SPBZJmpi->setEnabled(true);
}
    else{
        FRMZJmpi->setHidden(true);
        CHKZJmpi->setChecked(false);
        SPBZJmpi->setEnabled(false);
    }
    connect(CHKZJmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKZJmpi_changed(int)));
    connect(SPBZJmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBZJmpi_changed(int)));

//            List file, Stop, Execute

    BTNexecDamZJ=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_ZJdens);
    BTNexecDamZJ->setToolTip(tr("Electron density analysis"));
    connect(BTNexecDamZJ, SIGNAL(clicked()), this, SLOT(execDamZJ()));

    if (!BTNstop[9])
        BTNstop[9]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_ZJdens);
    BTNstop[9]->setToolTip(tr("Kill the process"));
    connect(BTNstop[9], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[9])
        BTNtexto[9]=new QPushButton(QIcon(":/images/document_text.png"),"",page_ZJdens);
    BTNtexto[9]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[9], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_ZJdens_layouts()
{
    QVBoxLayout *Layout0=new QVBoxLayout();
    Layout0->addWidget(RBTZJrstarabs);
    Layout0->addWidget(RBTZJrstarrel);

    QHBoxLayout *Layout1=new QHBoxLayout(FRMZJrstartype);
    Layout1->addStretch();
    Layout1->addLayout(Layout0);
    Layout1->addStretch();

    QVBoxLayout *Layout2=new QVBoxLayout();
    Layout2->addWidget(TXTZJrstar);
    Layout2->addWidget(FRMZJrstartype);

    QHBoxLayout *Layout3=new QHBoxLayout(FRMZJrstar);
    Layout3->addStretch();
    Layout3->addLayout(Layout2);
    Layout3->addStretch();

    QGridLayout *Layout4=new QGridLayout();
    Layout4->addWidget(LBLZJlmax,0,0,Qt::AlignLeft);
    Layout4->addWidget(SPBZJlmax,0,1,Qt::AlignCenter);
    Layout4->addWidget(LBLZJkmax,1,0,Qt::AlignLeft);
    Layout4->addWidget(SPBZJkmax,1,1,Qt::AlignCenter);

    QHBoxLayout *Layout5=new QHBoxLayout();
    Layout5->addStretch();
    Layout5->addLayout(Layout4);
    Layout5->addStretch();

    QHBoxLayout *Layout6=new QHBoxLayout();
    Layout6->addStretch();
    Layout6->addWidget(RBTZJechelon);
    Layout6->addStretch();

    QVBoxLayout *Layout7=new QVBoxLayout(FRMZJlength);
    Layout7->addLayout(Layout5);
    Layout7->addLayout(Layout6);

    QVBoxLayout *Layout8=new QVBoxLayout(FRMZJtype);
    Layout8->addWidget(RBTZJZernike);
    Layout8->addWidget(RBTZJacobi);

    QHBoxLayout *Layout9=new QHBoxLayout(FRMZJnquad);
    Layout9->addWidget(LBLZJnquad);
    Layout9->addWidget(SPBZJnquad);

    QHBoxLayout *Layout10=new QHBoxLayout();
    Layout10->addStretch();
    Layout10->addWidget(LBLZJthrmult);
    Layout10->addWidget(SPBZJthrmult);

    QHBoxLayout *Layout11=new QHBoxLayout();
    Layout11->addStretch();
    Layout11->addWidget(LBLZJthrdist);
    Layout11->addWidget(SPBZJthrdist);

    QVBoxLayout *Layout12=new QVBoxLayout(FRMZJthreshold);
    Layout12->addLayout(Layout10);
    Layout12->addLayout(Layout11);

    QVBoxLayout *Layout13 = new QVBoxLayout(FRMZJinput);
    Layout13->addWidget(CHKZJinput);

    QHBoxLayout *Layout14 = new QHBoxLayout(FRMZJmpi);
    Layout14->addWidget(CHKZJmpi);
    Layout14->addWidget(LBLZJmpi);
    Layout14->addWidget(SPBZJmpi);

    QHBoxLayout *Layout15=new QHBoxLayout();
    Layout15->addWidget(BTNtexto[9],0,Qt::AlignLeft);
    Layout15->addWidget(BTNstop[9],0,Qt::AlignRight);
    Layout15->addWidget(BTNexecDamZJ,0,Qt::AlignRight);

    QVBoxLayout *page_ZJdensLayout = new QVBoxLayout(page_ZJdens);
    page_ZJdensLayout->addWidget(FRMZJrstar);
    page_ZJdensLayout->addWidget(FRMZJlength);
    page_ZJdensLayout->addWidget(FRMZJtype);
    page_ZJdensLayout->addWidget(FRMZJnquad);
    page_ZJdensLayout->addWidget(FRMZJthreshold);
    page_ZJdensLayout->addWidget(FRMZJinput);
    page_ZJdensLayout->addWidget(FRMZJmpi);
    page_ZJdensLayout->addLayout(Layout15);
    page_ZJdensLayout->addStretch();
}
//    page_ZJtab: TABULATION OF DENSITY FROM ZERNIKE-JACOBI EXPANSIONS
//    ================================================================

void MainWindow::page_ZJtab_widgets()
{
//            Import File

    FRMdZJImportfile = new QGroupBox(tr("Zernike or Jacobi expansion"),page_ZJtab);
    FRMdZJImportfile->setMaximumSize(QSize(400, 2000));

    LBLdZJImportFile = new QLabel(tr("Import data from")+":");

    TXTdZJImportfile = new QLineEdit();

    BTNdZJImportFile = new QToolButton();
    BTNdZJImportFile->setText(tr("..."));
    connect(BTNdZJImportFile, SIGNAL(clicked()), this, SLOT(ImportFileNameZJ()));

//            Output file prefix

    FRMdZJfilename = new QGroupBox(tr("Output files prefix"),page_ZJtab);
    FRMdZJfilename->setMaximumSize(QSize(400, 2000));

    TXTdZJfilename = new QLineEdit(FRMdZJfilename);
    connect(TXTdZJfilename, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

//            Density expansion terms

    FRMdZJexpansion = new QGroupBox(tr("Density expansion terms"),page_ZJtab);
    FRMdZJexpansion->setMaximumSize(QSize(400, 2000));

    LBLdZJchoose = new QLabel("");
    LBLdZJchoose->setVisible(false);

    RBTdZJchooseall=new QRadioButton(tr("Choose all functions"),FRMdZJexpansion);
    RBTdZJchooseall->setChecked(true);
    RBTdZJchooseall->setToolTip(tr("Choose all available functions for projection"));
    connect(RBTdZJchooseall, SIGNAL(clicked()), this, SLOT(ChooseZJ_changed()));

    RBTdZJchoosel=new QRadioButton(tr("Choose l indices"),FRMdZJexpansion);
    RBTdZJchoosel->setChecked(false);
    RBTdZJchoosel->setToolTip(tr("Choose l indices, all compatible k and m indices are taken"));
    connect(RBTdZJchoosel, SIGNAL(clicked()), this, SLOT(ChooseZJ_changed()));

    RBTdZJchoosek=new QRadioButton(tr("Choose k indices"),FRMdZJexpansion);
    RBTdZJchoosek->setChecked(false);
    RBTdZJchoosek->setToolTip(tr("Choose k indices, all compatible l and m indices are taken"));
    connect(RBTdZJchoosek, SIGNAL(clicked()), this, SLOT(ChooseZJ_changed()));

    RBTdZJchooselk=new QRadioButton(tr("Choose pairs of (l,k) indices"),FRMdZJexpansion);
    RBTdZJchooselk->setChecked(false);
    RBTdZJchooselk->setToolTip(tr("Choose pairs of (l,k) indices, all compatible m indices are taken"));
    connect(RBTdZJchooselk, SIGNAL(clicked()), this, SLOT(ChooseZJ_changed()));

    RBTdZJchooselkm=new QRadioButton(tr("Choose triads of (l,k,m) indices"),FRMdZJexpansion);
    RBTdZJchooselkm->setChecked(false);
    RBTdZJchooselkm->setToolTip(tr("Choose triads of (l,k,m) indices"));
    connect(RBTdZJchooselkm, SIGNAL(clicked()), this, SLOT(ChooseZJ_changed()));

    TXTdZJchoose = new QLineEdit(FRMdZJexpansion);
    TXTdZJchoose->setToolTip(tr("Set selected indices here, separated by commas for individual values, or by hyphens for ranges where applicable"));
    TXTdZJchoose->setVisible(false);
    QRegExp rxZF("[0-9][-,\\d]*");
    QValidator *ZJvalidator = new QRegExpValidator(rxZF, this);
    ZJlist = new QStringList();
    TXTdZJchoose->setValidator(ZJvalidator);
    connect(TXTdZJchoose, SIGNAL(textChanged(const QString &)), this, SLOT(TXTdZJchoose_changed()));

    FRMdZJexplength = new QGroupBox(tr(""),page_ZJtab);
    FRMdZJexplength->setMaximumSize(QSize(300, 400));
    FRMdZJexplength->setVisible(true);

    LBLdZJkmax=new QLabel(tr("Highest k"),FRMdZJexplength);

    SPBdZJkmax=new QSpinBox(FRMdZJexplength);
    SPBdZJkmax->setRange(0, SPBZJkmax->value());
    SPBdZJkmax->setValue(SPBZJkmax->value());
    SPBdZJkmax->setMaximumWidth(60);
    SPBdZJkmax->setToolTip(tr("Highest number of Zernike or Jacobi functions per l"));

    LBLdZJlmax=new QLabel(tr("Highest l"),FRMdZJexplength);

    SPBdZJlmax=new QSpinBox(FRMdZJexplength);
    SPBdZJlmax->setRange(0, SPBZJlmax->value());
    SPBdZJlmax->setValue(SPBZJlmax->value());
    SPBdZJlmax->setMaximumWidth(60);
    SPBdZJlmax->setToolTip(tr("Highest l of atomic multipolar expansion"));

    LBLdZJlmin=new QLabel(tr("Lowest l "),FRMdZJexplength);

    SPBdZJlmin=new QSpinBox(FRMdZJexplength);
    SPBdZJlmin->setRange(0, SPBZJlmax->value());
    SPBdZJlmin->setValue(0);
    SPBdZJlmin->setMaximumWidth(60);
    SPBdZJlmin->setToolTip(tr("Lowest l of atomic multipolar expansion"));

//            Derivatives

    FRMdZJderivs = new QGroupBox(tr("Derivatives"),page_ZJtab);
    FRMdZJderivs->setMaximumSize(QSize(400, 2000));

    CHKdZJgrad = new QCheckBox(tr("Gradient"),FRMdZJderivs);
    CHKdZJgrad->setToolTip(tr("Tabulate ZJ density gradient"));
    CHKdZJgrad->setChecked(true);

    QDoubleValidator *myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);

//            Grid

    FRMdZJgrid = new QGroupBox(tr("Grid"),page_ZJtab);
    FRMdZJgrid->setMaximumSize(QSize(400, 2000));

    CHKdZJgrid=new QCheckBox(tr("Generate grid"),FRMdZJgrid);
    CHKdZJgrid->setChecked(true);
    connect(CHKdZJgrid, SIGNAL(stateChanged(int)), this, SLOT(CHKdZJgrid_changed()));

    FRMdZJgridtype = new QGroupBox(tr("Grid type"));

//          2D Grid

    RBTdZJ2D=new QRadioButton(tr("2D grid"),FRMdZJgridtype);
    RBTdZJ2D->setChecked(false);
    connect(RBTdZJ2D, SIGNAL(toggled (bool)), this, SLOT(RBTdZJ2D3D_changed()));

    FRMdZJgrid2D = new QGroupBox(tr("2D grid"));
    FRMdZJgrid2D->setMaximumSize(QSize(400, 2000));
    FRMdZJgrid2D->setVisible(false);

    LBLdZJu=new QLabel(tr("u"),FRMdZJgrid2D);
    LBLdZJv=new QLabel(tr("v"),FRMdZJgrid2D);
    LBLdZJinf2d=new QLabel(tr("Lowest"),FRMdZJgrid2D);
    LBLdZJsup2d=new QLabel(tr("Highest"),FRMdZJgrid2D);

    TXTdZJuinf=new QLineEdit(FRMdZJgrid2D);
    TXTdZJuinf->setText("-10.0");
    TXTdZJuinf->setAlignment(Qt::AlignRight);
    TXTdZJuinf->setValidator(myDoubleValidator);
    connect(TXTdZJuinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJusup=new QLineEdit(FRMdZJgrid2D);
    TXTdZJusup->setText("10.0");
    TXTdZJusup->setAlignment(Qt::AlignRight);
    TXTdZJusup->setValidator(myDoubleValidator);
    connect(TXTdZJusup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJvinf=new QLineEdit(FRMdZJgrid2D);
    TXTdZJvinf->setText("-10.0");
    TXTdZJvinf->setAlignment(Qt::AlignRight);
    TXTdZJvinf->setValidator(myDoubleValidator);
    connect(TXTdZJvinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJvsup=new QLineEdit(FRMdZJgrid2D);
    TXTdZJvsup->setText("10.0");
    TXTdZJvsup->setAlignment(Qt::AlignRight);
    TXTdZJvsup->setValidator(myDoubleValidator);
    connect(TXTdZJvsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdZJsurfpar = new QGroupBox(tr("Parametric equations"),FRMdZJgrid2D);
    FRMdZJsurfpar->setMaximumSize(QSize(400, 2000));
    FRMdZJsurfpar->setVisible(false);

    LBLdZJxformula2D=new QLabel(tr("x(u,v) = "),FRMdZJsurfpar);
    LBLdZJyformula2D=new QLabel(tr("y(u,v) = "),FRMdZJsurfpar);
    LBLdZJzformula2D=new QLabel(tr("z(u,v) = "),FRMdZJsurfpar);

    TXTdZJxformula2D=new QLineEdit(FRMdZJsurfpar);
    TXTdZJxformula2D->setText("u");

    TXTdZJyformula2D=new QLineEdit(FRMdZJsurfpar);
    TXTdZJyformula2D->setText("v");

    TXTdZJzformula2D=new QLineEdit(FRMdZJsurfpar);
    TXTdZJzformula2D->setText("0");

    FRMdZJsurftype = new QGroupBox(tr("Surface type"),FRMdZJgrid2D);
    FRMdZJsurftype->setMaximumSize(QSize(400, 2000));

    RBTdZJplane = new QRadioButton(tr("Plane"),FRMdZJsurftype);
    RBTdZJplane->setChecked(true);
    connect(RBTdZJplane, SIGNAL(toggled (bool)), this, SLOT(RBTdZJplane_changed()));

    RBTdZJothersurf = new QRadioButton(tr("Parametric surface"),FRMdZJsurftype);
    RBTdZJothersurf->setToolTip(tr("Surface equation supplied in parametric form: x=x(u,v), y=y(u,v), z=z(u,v)"));
    RBTdZJothersurf->setChecked(false);

    FRMdZJplane2D = new QGroupBox(tr("2D Planes"),FRMdZJgrid2D);
    FRMdZJplane2D->setMaximumSize(QSize(400, 2000));
    FRMdZJplane2D->setVisible(true);

    RBTdZJplaneXY = new QRadioButton(tr("XY "),FRMdZJgrid2D);
    RBTdZJplaneXY->setChecked(true);
    connect(RBTdZJplaneXY, SIGNAL(toggled (bool)), this, SLOT(RBTdZJ2Dplanes_changed()));

    RBTdZJplaneXZ = new QRadioButton(tr("XZ "),FRMdZJgrid2D);
    RBTdZJplaneXZ->setChecked(false);
    connect(RBTdZJplaneXZ, SIGNAL(toggled (bool)), this, SLOT(RBTdZJ2Dplanes_changed()));

    RBTdZJplaneYZ = new QRadioButton(tr("YZ "),FRMdZJgrid2D);
    RBTdZJplaneYZ->setChecked(false);
    connect(RBTdZJplaneYZ, SIGNAL(toggled (bool)), this, SLOT(RBTdZJ2Dplanes_changed()));

    RBTdZJplaneABC = new QRadioButton(tr("Other "),FRMdZJgrid2D);
    RBTdZJplaneABC->setChecked(false);
    connect(RBTdZJplaneABC, SIGNAL(toggled (bool)), this, SLOT(RBTdZJ2Dplanes_changed()));

    FRMdZJplaneABC = new QGroupBox(tr("Plane parameters"),FRMdZJgrid2D);
    FRMdZJplaneABC->setMaximumSize(QSize(400, 2000));
    FRMdZJplaneABC->setVisible(false);

    TXTdZJplaneA=new QLineEdit(FRMdZJplaneABC);
    TXTdZJplaneA->setText("0.");
    TXTdZJplaneA->setValidator(myDoubleValidator);
    connect(TXTdZJplaneA, SIGNAL(editingFinished()), this, SLOT(RBTdZJ2Dplanes_changed()));

    TXTdZJplaneB=new QLineEdit(FRMdZJplaneABC);
    TXTdZJplaneB->setText("0.");
    TXTdZJplaneB->setValidator(myDoubleValidator);
    connect(TXTdZJplaneB, SIGNAL(editingFinished()), this, SLOT(RBTdZJ2Dplanes_changed()));

    TXTdZJplaneC=new QLineEdit(FRMdZJplaneABC);
    TXTdZJplaneC->setText("1.");
    TXTdZJplaneC->setValidator(myDoubleValidator);
    connect(TXTdZJplaneC, SIGNAL(editingFinished()), this, SLOT(RBTdZJ2Dplanes_changed()));

    dZJplanecase = 1;

    FRMdZJresol2D = new QGroupBox(tr("Custom resolution"));
    FRMdZJresol2D->setHidden(true);

    LBLdZJuresol = new QLabel(tr("u"),FRMdZJresol2D);
    LBLdZJvresol = new QLabel(tr("v"),FRMdZJresol2D);

    SPBdZJures = new QSpinBox(FRMdZJresol2D);
    SPBdZJures->setMinimum(4);
    SPBdZJures->setMaximum(2049);
    SPBdZJures->setValue(129);
    SPBdZJures->setSingleStep(10);

    SPBdZJvres = new QSpinBox(FRMdZJresol2D);
    SPBdZJvres->setMinimum(4);
    SPBdZJvres->setMaximum(2049);
    SPBdZJvres->setValue(129);
    SPBdZJvres->setSingleStep(10);

//          3D Grid

    RBTdZJ3D=new QRadioButton(tr("3D grid"),FRMdZJgridtype);
    RBTdZJ3D->setChecked(true);
    connect(RBTdZJ3D, SIGNAL(toggled (bool)), this, SLOT(RBTdZJ2D3D_changed()));

    FRMdZJgrid3D = new QGroupBox(tr("3D grid"));
    FRMdZJgrid3D->setMaximumSize(QSize(400, 2000));
    FRMdZJgrid3D->setVisible(true);

    LBLdZJx=new QLabel(tr("x"),FRMdZJgrid3D);
    LBLdZJy=new QLabel(tr("y"),FRMdZJgrid3D);
    LBLdZJz=new QLabel(tr("z"),FRMdZJgrid3D);
    LBLdZJinf=new QLabel(tr("Lowest"),FRMdZJgrid3D);
    LBLdZJsup=new QLabel(tr("Highest"),FRMdZJgrid3D);

    TXTdZJxinf=new QLineEdit(FRMdZJgrid3D);
    TXTdZJxinf->setText("-10.0");
    TXTdZJxinf->setAlignment(Qt::AlignRight);
    TXTdZJxinf->setValidator(myDoubleValidator);
    connect(TXTdZJxinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJxsup=new QLineEdit(FRMdZJgrid3D);
    TXTdZJxsup->setText("10.0");
    TXTdZJxsup->setAlignment(Qt::AlignRight);
    TXTdZJxsup->setValidator(myDoubleValidator);
    connect(TXTdZJxsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJyinf=new QLineEdit(FRMdZJgrid3D);
    TXTdZJyinf->setText("-10.0");
    TXTdZJyinf->setAlignment(Qt::AlignRight);
    TXTdZJyinf->setValidator(myDoubleValidator);
    connect(TXTdZJyinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJysup=new QLineEdit(FRMdZJgrid3D);
    TXTdZJysup->setText("10.0");
    TXTdZJysup->setAlignment(Qt::AlignRight);
    TXTdZJysup->setValidator(myDoubleValidator);
    connect(TXTdZJysup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJzinf=new QLineEdit(FRMdZJgrid3D);
    TXTdZJzinf->setText("-10.0");
    TXTdZJzinf->setAlignment(Qt::AlignRight);
    TXTdZJzinf->setValidator(myDoubleValidator);
    connect(TXTdZJzinf, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    TXTdZJzsup=new QLineEdit(FRMdZJgrid3D);
    TXTdZJzsup->setText("10.0");
    TXTdZJzsup->setAlignment(Qt::AlignRight);
    TXTdZJzsup->setValidator(myDoubleValidator);
    connect(TXTdZJzsup, SIGNAL(textChanged(const QString &)), this, SLOT(TXTValidate_changed()));

    FRMdZJresol3D = new QGroupBox(tr("Custom resolution"));
    FRMdZJresol3D->setHidden(true);

    LBLdZJxresol = new QLabel(tr("x"),FRMdZJresol3D);
    LBLdZJyresol = new QLabel(tr("y"),FRMdZJresol3D);
    LBLdZJzresol = new QLabel(tr("z"),FRMdZJresol3D);

    SPBdZJxres = new QSpinBox(FRMdZJresol3D);
    SPBdZJxres->setMinimum(4);
    SPBdZJxres->setMaximum(2049);
    SPBdZJxres->setValue(129);
    SPBdZJxres->setSingleStep(10);

    SPBdZJyres = new QSpinBox(FRMdZJresol3D);
    SPBdZJyres->setMinimum(4);
    SPBdZJyres->setMaximum(2049);
    SPBdZJyres->setValue(129);
    SPBdZJyres->setSingleStep(10);

    SPBdZJzres = new QSpinBox(FRMdZJresol3D);
    SPBdZJzres->setMinimum(4);
    SPBdZJzres->setMaximum(2049);
    SPBdZJzres->setValue(129);
    SPBdZJzres->setSingleStep(10);

    FRMdZJxyz = new QGroupBox(tr("Tabulation"),page_ZJtab);
    FRMdZJxyz->setMaximumSize(QSize(400, 2000));

    CHKdZJxyz=new QCheckBox(tr("Tabulation points"),FRMdZJxyz);
    CHKdZJxyz->setChecked(false);
    connect(CHKdZJxyz, SIGNAL(stateChanged(int)), this, SLOT(CHKdZJxyz_changed()));

    WtabledZJ=new QWidget(FRMdZJxyz);
    SHTdZJxyz = new Sheet(0, 3, 0,true, WtabledZJ);
    QStringList QSLdZJxyz;
    QSLdZJxyz << "x" << "y" << "z";
    SHTdZJxyz->setHeader(QSLdZJxyz);
    WtabledZJ->setVisible(false);
    WtabledZJ->setEnabled(false);

//          Resolution

    FRMdZJgridres = new QGroupBox(tr("Resolution"));

    RBTdZJrlow=new QRadioButton(tr("Low"),FRMdZJgridres);
    RBTdZJrlow->setChecked(true);
    connect(RBTdZJrlow, SIGNAL(toggled (bool)), this, SLOT(ZJ_resolution_changed()));

    RBTdZJrmedium=new QRadioButton(tr("Medium"),FRMdZJgridres);
    connect(RBTdZJrmedium, SIGNAL(toggled (bool)), this, SLOT(ZJ_resolution_changed()));

    RBTdZJrhigh=new QRadioButton(tr("High"),FRMdZJgridres);
    connect(RBTdZJrhigh, SIGNAL(toggled (bool)), this, SLOT(ZJ_resolution_changed()));

    RBTdZJrcustom=new QRadioButton(tr("Custom"),FRMdZJgridres);
    connect(RBTdZJrcustom, SIGNAL(toggled (bool)), this, SLOT(ZJ_resolution_changed()));

    if (RBTdZJ3D->isChecked()){
        RBTdZJrlow->setToolTip("65x65x65");
        RBTdZJrmedium->setToolTip("129x129x129");
        RBTdZJrhigh->setToolTip("257x257x257");
    }
    else{
        RBTdZJrlow->setToolTip("129x129");
        RBTdZJrmedium->setToolTip("257x257");
        RBTdZJrhigh->setToolTip("513x513");
    }

//                      Generate input only

    FRMdZJinput = new QGroupBox(tr("Input only"),page_ZJtab);
    FRMdZJinput->setMaximumSize(QSize(400, 2000));

    CHKdZJinput=new QCheckBox(tr("Generate input file only"),FRMdZJinput);
    CHKdZJinput->setChecked(false);
    CHKdZJinput->setEnabled(true);
    connect(CHKdZJinput, SIGNAL(stateChanged(int)), this, SLOT(CHKdZJinput_changed(int)));

//            MPI

    FRMdZJmpi = new QGroupBox(tr("Parallel computing"),page_ZJtab);
    FRMdZJmpi->setMaximumSize(QSize(400, 2000));

    CHKdZJmpi = new QCheckBox(tr("MPI"),FRMdZJmpi);

    LBLdZJmpi = new QLabel(tr("Number of processors"),FRMdZJmpi);

    SPBdZJmpi = new QSpinBox(FRMdZJmpi);
    SPBdZJmpi->setRange(1, MAX_NUM_PROCESSORS);
    SPBdZJmpi->setValue(1);
    SPBdZJmpi->setMaximumWidth(60);
    SPBdZJmpi->setToolTip(tr("Number of processors"));
    if (mpi){
        FRMdZJmpi->setVisible(true);
        FRMdZJmpi->setEnabled(true);
        CHKdZJmpi->setChecked(true);
        SPBdZJmpi->setEnabled(true);
    }
    else{
        FRMdZJmpi->setHidden(true);
        CHKdZJmpi->setChecked(false);
        SPBdZJmpi->setEnabled(false);
    }
    connect(CHKdZJmpi, SIGNAL(stateChanged(int)), this, SLOT(CHKdZJmpi_changed(int)));
    connect(SPBdZJmpi, SIGNAL(valueChanged(int)), this, SLOT(SPBdZJmpi_changed(int)));

//            List file, Stop, Execute
    BTNexecDamdZJ=new QPushButton(QIcon(":/images/exec.png"), tr("Exec"),page_ZJtab);
    BTNexecDamdZJ->setToolTip(tr("Generate grids"));
    connect(BTNexecDamdZJ, SIGNAL(clicked()), this, SLOT(execDamdenZJ()));

    if (!BTNstop[10])
        BTNstop[10]=new QPushButton(QIcon(":/images/stop.png"), tr("Stop"),page_ZJtab);
    BTNstop[10]->setToolTip(tr("Kill the process"));
    connect(BTNstop[10], SIGNAL(clicked()), this, SLOT(processStop()));

    if (!BTNtexto[10])
        BTNtexto[10]=new QPushButton(QIcon(":/images/document_text.png"),"",page_ZJtab);
    BTNtexto[10]->setToolTip(tr("Open output file ..."));
    connect(BTNtexto[10], SIGNAL(clicked()), this, SLOT(importOUT()));
}

void MainWindow::page_ZJtab_layouts()
{
    QLabel *LBLdZJplaneeq=new QLabel(tr("(A x + B y + C z = 0)"),FRMdZJplaneABC);
    QLabel *LBLdZJplaneA=new QLabel(tr("A: "),FRMdZJplaneABC);
    QLabel *LBLdZJplaneB=new QLabel(tr("B: "),FRMdZJplaneABC);
    QLabel *LBLdZJplaneC=new QLabel(tr("C: "),FRMdZJplaneABC);
    QHBoxLayout *Layout0=new QHBoxLayout();
    Layout0->addWidget(LBLdZJImportFile);

    QHBoxLayout *Layout1=new QHBoxLayout();
    Layout1->addWidget(TXTdZJImportfile);
    Layout1->addWidget(BTNdZJImportFile);

    QVBoxLayout *Layout2=new QVBoxLayout(FRMdZJImportfile);
    Layout2->addLayout(Layout0);
    Layout2->addLayout(Layout1);

    QHBoxLayout *Layout3=new QHBoxLayout(FRMdZJfilename);
    Layout3->addWidget(TXTdZJfilename);

    QGridLayout *Layout4=new QGridLayout(FRMdZJexplength);
    Layout4->addWidget(LBLdZJlmax,0,0,Qt::AlignLeft);
    Layout4->addWidget(SPBdZJlmax,0,1,Qt::AlignCenter);
    Layout4->addWidget(LBLdZJlmin,1,0,Qt::AlignLeft);
    Layout4->addWidget(SPBdZJlmin,1,1,Qt::AlignCenter);
    Layout4->addWidget(LBLdZJkmax,2,0,Qt::AlignLeft);
    Layout4->addWidget(SPBdZJkmax,2,1,Qt::AlignCenter);

    QVBoxLayout *Layout5=new QVBoxLayout();
    Layout5->addWidget(RBTdZJchooseall,Qt::AlignLeft);
    Layout5->addWidget(FRMdZJexplength,Qt::AlignCenter);
    Layout5->addWidget(RBTdZJchoosel,Qt::AlignLeft);
    Layout5->addWidget(RBTdZJchoosek,Qt::AlignLeft);
    Layout5->addWidget(RBTdZJchooselk,Qt::AlignLeft);
    Layout5->addWidget(RBTdZJchooselkm,Qt::AlignLeft);
    Layout5->addWidget(LBLdZJchoose,Qt::AlignLeft);
    Layout5->addWidget(TXTdZJchoose,Qt::AlignLeft);

    QHBoxLayout *Layout6=new QHBoxLayout(FRMdZJexpansion);
    Layout6->addStretch();
    Layout6->addLayout(Layout5);
    Layout6->addStretch();

    QVBoxLayout *Layout7=new QVBoxLayout(FRMdZJderivs);
    Layout7->addWidget(CHKdZJgrad);

    QHBoxLayout *Layout8=new QHBoxLayout();
    Layout8->addWidget(CHKdZJgrid,0,Qt::AlignCenter);

    QHBoxLayout *Layout9=new QHBoxLayout(FRMdZJgridtype);
    Layout9->addWidget(RBTdZJ2D);
    Layout9->addWidget(RBTdZJ3D);

    QGridLayout *Layout10=new QGridLayout(FRMdZJgrid3D);
    Layout10->addWidget(LBLdZJinf,1,1,Qt::AlignCenter);
    Layout10->addWidget(LBLdZJsup,1,2,Qt::AlignCenter);
    Layout10->addWidget(LBLdZJx,2,0);
    Layout10->addWidget(TXTdZJxinf,2,1);
    Layout10->addWidget(TXTdZJxsup,2,2);
    Layout10->addWidget(LBLdZJy,3,0);
    Layout10->addWidget(TXTdZJyinf,3,1);
    Layout10->addWidget(TXTdZJysup,3,2);
    Layout10->addWidget(LBLdZJz,4,0);
    Layout10->addWidget(TXTdZJzinf,4,1);
    Layout10->addWidget(TXTdZJzsup,4,2);

    QGridLayout *Layout11=new QGridLayout();
    Layout11->addWidget(LBLdZJinf2d,1,1,Qt::AlignCenter);
    Layout11->addWidget(LBLdZJsup2d,1,2,Qt::AlignCenter);
    Layout11->addWidget(LBLdZJu,2,0);
    Layout11->addWidget(TXTdZJuinf,2,1);
    Layout11->addWidget(TXTdZJusup,2,2);
    Layout11->addWidget(LBLdZJv,3,0);
    Layout11->addWidget(TXTdZJvinf,3,1);
    Layout11->addWidget(TXTdZJvsup,3,2);

    QVBoxLayout *Layout11b=new QVBoxLayout(FRMdZJsurftype);
    Layout11b->addWidget(RBTdZJplane);
    Layout11b->addWidget(RBTdZJothersurf);

    QGridLayout *Layout11c=new QGridLayout(FRMdZJplaneABC);
    Layout11c->addWidget(LBLdZJplaneeq,0,0,1,2,Qt::AlignCenter);
    Layout11c->addWidget(LBLdZJplaneA,1,0);
    Layout11c->addWidget(TXTdZJplaneA,1,1);
    Layout11c->addWidget(LBLdZJplaneB,2,0);
    Layout11c->addWidget(TXTdZJplaneB,2,1);
    Layout11c->addWidget(LBLdZJplaneC,3,0);
    Layout11c->addWidget(TXTdZJplaneC,3,1);

    QGridLayout *Layout11d=new QGridLayout(FRMdZJplane2D);
    Layout11d->addWidget(RBTdZJplaneXY,0,0,Qt::AlignLeft);
    Layout11d->addWidget(RBTdZJplaneXZ,0,1,Qt::AlignLeft);
    Layout11d->addWidget(RBTdZJplaneYZ,1,0,Qt::AlignLeft);
    Layout11d->addWidget(RBTdZJplaneABC,1,1,Qt::AlignLeft);
    Layout11d->addWidget(FRMdZJplaneABC,2,0,1,2,Qt::AlignCenter);

    QGridLayout *Layout12=new QGridLayout(FRMdZJsurfpar);
    Layout12->addWidget(LBLdZJxformula2D,1,0);
    Layout12->addWidget(TXTdZJxformula2D,1,1);
    Layout12->addWidget(LBLdZJyformula2D,2,0);
    Layout12->addWidget(TXTdZJyformula2D,2,1);
    Layout12->addWidget(LBLdZJzformula2D,3,0);
    Layout12->addWidget(TXTdZJzformula2D,3,1);

    QVBoxLayout *Layout13=new QVBoxLayout(FRMdZJgrid2D);
    Layout13->addLayout(Layout11,Qt::AlignCenter);
    Layout13->addWidget(FRMdZJsurftype);
    Layout13->addWidget(FRMdZJplane2D);
    Layout13->addWidget(FRMdZJsurfpar);


    QVBoxLayout *Layout14=new QVBoxLayout();
    Layout14->addWidget(FRMdZJgrid3D);
    Layout14->addWidget(FRMdZJgrid2D);

    QHBoxLayout *Layout15=new QHBoxLayout();
    Layout15->addWidget(RBTdZJrlow);
    Layout15->addWidget(RBTdZJrmedium);

    QHBoxLayout *Layout16=new QHBoxLayout();
    Layout16->addWidget(RBTdZJrhigh);
    Layout16->addWidget(RBTdZJrcustom);

    QGridLayout *Layout17=new QGridLayout(FRMdZJresol2D);
    Layout17->addWidget(LBLdZJuresol,0,0,Qt::AlignCenter);
    Layout17->addWidget(LBLdZJvresol,0,1,Qt::AlignCenter);
    Layout17->addWidget(SPBdZJures,1,0);
    Layout17->addWidget(SPBdZJvres,1,1);

    QGridLayout *Layout18=new QGridLayout(FRMdZJresol3D);
    Layout18->addWidget(LBLdZJxresol,0,0,Qt::AlignCenter);
    Layout18->addWidget(LBLdZJyresol,0,1,Qt::AlignCenter);
    Layout18->addWidget(LBLdZJzresol,0,2,Qt::AlignCenter);
    Layout18->addWidget(SPBdZJxres,1,0);
    Layout18->addWidget(SPBdZJyres,1,1);
    Layout18->addWidget(SPBdZJzres,1,2);

    QVBoxLayout *Layout19=new QVBoxLayout(FRMdZJgridres);
    Layout19->addLayout(Layout15);
    Layout19->addLayout(Layout16);
    Layout19->addWidget(FRMdZJresol2D);
    Layout19->addWidget(FRMdZJresol3D);

    QVBoxLayout *Layout20=new QVBoxLayout(FRMdZJgrid);
    Layout20->addLayout(Layout8);
    Layout20->addWidget(FRMdZJgridtype);
    Layout20->addLayout(Layout14);
    Layout20->addWidget(FRMdZJgridres);

    QVBoxLayout *Layout21 = new QVBoxLayout(FRMdZJxyz);
    Layout21->addWidget(CHKdZJxyz,0,Qt::AlignCenter);
    Layout21->addWidget(WtabledZJ,0,Qt::AlignCenter);

    QVBoxLayout *Layout22 = new QVBoxLayout(FRMdZJinput);
    Layout22->addWidget(CHKdZJinput);

    QHBoxLayout *Layout23 = new QHBoxLayout(FRMdZJmpi);
    Layout23->addWidget(CHKdZJmpi);
    Layout23->addWidget(LBLdZJmpi);
    Layout23->addWidget(SPBdZJmpi);

    QHBoxLayout *Layout24=new QHBoxLayout();
    Layout24->addWidget(BTNtexto[10],0,Qt::AlignLeft);
    Layout24->addWidget(BTNstop[10],0,Qt::AlignRight);
    Layout24->addWidget(BTNexecDamdZJ,0,Qt::AlignRight);

    QVBoxLayout *page_ZJtabLayout = new QVBoxLayout(page_ZJtab);
    page_ZJtabLayout->addWidget(FRMdZJImportfile);
    page_ZJtabLayout->addWidget(FRMdZJfilename);
    page_ZJtabLayout->addWidget(FRMdZJexpansion);
    page_ZJtabLayout->addWidget(FRMdZJderivs);
    page_ZJtabLayout->addWidget(FRMdZJgrid);
    page_ZJtabLayout->addWidget(FRMdZJxyz);
    page_ZJtabLayout->addWidget(FRMdZJinput);
    page_ZJtabLayout->addWidget(FRMdZJmpi);
    page_ZJtabLayout->addLayout(Layout24);
    page_ZJtabLayout->addStretch();
}

/***************************************************************************/
/*  page_project: PROJECT                                             */
/***************************************************************************/

void MainWindow::chooseLanguage(QAction *action)
{
    QString locale = action->data().toString();
    QString qmPath = ":/translations/";
    if (locale.isEmpty()) locale.append("en");
    appTranslator.load("DAMQT_"+locale,qmPath);
    appTranslator.setObjectName("DAMQT_"+locale);
    qApp->installTranslator(&appTranslator);
    qtTranslator.load("qt_"+locale,qmPath);
    qApp->installTranslator(&qtTranslator);
    FRMlanguage->setWhatsThis(tr("Check a language and push Start to start DAMQT."));
    FRMlanguage->setWindowTitle(tr("Language"));
    LBLlanguage->setText(tr("Choose language and push Start"));
    languageMenu->setWhatsThis(tr("Check a language and push Start to start DAMQT."));
    BTNlangstart->setText(tr("Start"));
    if (locale == QString("es")){
        QPixmap pixmap;
        pixmap.load(QString(":/images/splash_3.2_%1.png").arg(locale));
        splash->setPixmap(pixmap);
        splash->show();
    }
    else{
        QPixmap pixmap;
        pixmap.load(QString(":/images/splash_3.2_en.png"));
        splash->setPixmap(pixmap);
        splash->show();
    }
}

void MainWindow::createLanguageMenu(){
    languageMenu = new QMenu(this);
    languageMenu->setWhatsThis(tr("Check a language and push Start to start DAMQT."));
    languageActionGroup = new QActionGroup(this);
    connect(languageActionGroup,SIGNAL(triggered(QAction *)),this, SLOT(chooseLanguage(QAction *)));
    QDir qmDir = QDir(":/translations");
    QStringList fileNames = qmDir.entryList(QStringList("DAMQT_*.qm"));
//    qDebug() << "fileNames = " << fileNames;
    for (int i = 0 ; i < fileNames.size(); ++i){
        QString locale = fileNames[i];
//        qDebug() << "locale = " << locale;
        locale.remove(0,locale.indexOf('_')+1);
        locale.chop(3);
        QTranslator translator;
        translator.load(fileNames[i], qmDir.absolutePath());
        QString language = translator.translate("MainWindow", "English");
//        qDebug() << "language = " << language;
        QAction *action = new QAction(tr("&%1 %2").arg(i+1).arg(language),this);
        action->setCheckable(true);
        action->setData(locale);
        languageMenu->addAction(action);
        languageActionGroup->addAction(action);
        if (language == "English")
            action->setChecked(true);
    }
}

/* Imports file */
void MainWindow::execImport()
{
    QString DirNombreImport = TXTImport->text();
    if (DirNombreImport.size()==0){
        QMessageBox::warning(this, tr("DAMQT"),tr("Navigate to a directory with a suitable import file "));
        return;
    }
    QFile file(DirNombreImport);
    if (!QFile::exists(DirNombreImport)){
        QMessageBox::warning(this, tr("DAMQT"),
                     tr("File %1 not found.").arg(DirNombreImport));
        return;
    }
    if (ProjectFolder.size()==0){
        QMessageBox::warning(this, tr("DAMQT"),tr("Introduce project folder"));
        return;
    }
    ImportFolder = QFileInfo(DirNombreImport).path();
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    ImportFile = QFileInfo(DirNombreImport).fileName();
    QString suffix=Extension(DirNombreImport);
    bool isgzipped = false;
    if (suffix =="gz"){
        int ios = QProcess::execute("gunzip "+DirNombreImport);
        if (ios != 0){
            return;
        }
        DirNombreImport = DirNombreImport.remove(-3,3);
        suffix=Extension(DirNombreImport);
        isgzipped = true;
    }
    if (suffix=="ggbs"){
        lslater = false;
        execGgbsDen();
    }else if (suffix=="sgbs" || suffix =="sgbsden"){
        lslater = true;
        execSxyzDen();
    }else if (suffix=="fchk"){
        lslater = false;
        execFchk();
    }else if (suffix=="mos" || suffix=="coord" || suffix=="basis"){
        lslater = false;
        execTurbom();
    }else if (suffix=="mkl"){
        lslater = false;
        execMOLEKEL();
    }else if (suffix=="out" || suffix=="xml"){
        execMolpro();
        lslater = false;
    }else if (suffix=="aux"){
        lslater = true;
        execMopac();
    }else if (suffix=="nwcout"){
        lslater = false;
        execNWChem();
    }else if (suffix=="psiauxden"){
        lslater = false;
        execPsi4auxiliary();
    }
    if (lslater){
        CHKpotexact->setVisible(false);
        CHKSGhexactMESP->setVisible(false);
    }
    else{
        CHKpotexact->setVisible(true);
        CHKSGhexactMESP->setVisible(true);
    }
    CHKpotexact->setChecked(false);
    CHKSGhexactMESP->setChecked(false);

    if (isgzipped){
        QProcess::execute("gzip "+DirNombreImport);
    }
}

/* Import name of a file */
void MainWindow::importFile()
{
    lzdo = false;
    lvalence = false;
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ImportFolder,
        tr("Import data from")+" (*.ggbs *.sgbs *.sgbsden *.sgbsden.gz *.fchk *.coord *.basis *.mos *.mkl "
            "*.out *.xml *.aux *.nwcout *.psiauxden);;"+
        tr("Geometry and basis set files")+" (*.ggbs *.sgbs *.sgbsden *.sgbsden.gz *.coord *.basis);;"+
        tr("fchk files")+" (*.fchk);;"+
        tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTImport->setText(fileName);
    ImportFolder = Path(fileName);
    ImportFile = QFileInfo(fileName).completeBaseName();
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    TXTProjectFolder->setText(Path(fileName));
    ProjectFolder = TXTProjectFolder->text();
    if (ProjectFolder.mid(ProjectFolder.length()-1,1) != "/") ProjectFolder.append("/");
    QString filezdo = ProjectFolder + "zdo";
    if (QFileInfo(filezdo).exists()){
        lzdo = true;
    }
    else{
        lzdo = false;
    }
    QString filevalence = ProjectFolder + "valence";
    if (QFileInfo(filevalence).exists()){
        lvalence = true;
    }
    else{
        lvalence = false;
    }
    TXTProjectName->setText(FileWithoutExt(TXTImport->text()));
    ProjectName = TXTProjectName->text();
    if (Extension(fileName)=="ggbs" || Extension(fileName)=="sgbs" || Extension(fileName)=="sgbsden"){
        bool res=Open(fileName);
        if (res) 
            QMessageBox::information(this,tr("DAMQT"),tr("Project %1 open").arg(ProjectName), QMessageBox::Ok, 0);
        else{
            QMessageBox::warning(this,tr("DAMQT"),tr("Cannot open project %1").arg(ProjectName));
            return;
        }
    }
    BTNexecImport->setEnabled(true);
    page_atdens->setEnabled(true);
    if (QFile::exists(QString(ImportFolder+FileWithoutExt(fileName)+"_2016.damqt"))){
        page_MED->setEnabled(true);
        page_MESP->setEnabled(true);
        page_TOPO->setEnabled(true);
        page_HFforces->setEnabled(true);
        page_Efield->setEnabled(true);
        page_densgrad->setEnabled(true);
        page_frad->setEnabled(true);
        page_orimult->setEnabled(true);
        page_ZJdens->setEnabled(true);
        page_ZJtab->setEnabled(true);
    }
    else{
        page_MED->setEnabled(false);
        page_MESP->setEnabled(false);
        page_TOPO->setEnabled(false);
        page_HFforces->setEnabled(false);
        page_Efield->setEnabled(false);
        page_densgrad->setEnabled(false);
        page_frad->setEnabled(false);
        page_orimult->setEnabled(false);
        page_ZJdens->setEnabled(false);
        page_ZJtab->setEnabled(false);
    }
    QStringList filenames(QString(FileWithoutExt(fileName) +".GAorb*"));
    filenames << QString(FileWithoutExt(fileName) +".SLorb*");
    filenames << QString(FileWithoutExt(fileName) +".orb*");
    QStringList files;
    files = QDir(ImportFolder).entryList(filenames,QDir::Files);
    if (!files.isEmpty()) 
        page_MO->setEnabled(true);
    else
        page_MO->setEnabled(false);
    
    page_ZJdens->setEnabled(true);
}

/* Textbox for introducing the mpi command */
void MainWindow::TXTmpicommand_changed()
{
    mpicommand = TXTmpicommand->text();
}

/* Textbox for introducing the mpi flags */
void MainWindow::TXTmpiflags_changed()
{
    mpiflags = TXTmpiflags->text();
}

/* Textbox for introducing the file to import */
void MainWindow::TXTImport_changed()
{
    if(TXTImport->text().isEmpty()){
        BTNexecImport->setEnabled(false);
    }else{
        ImportFile = FileWithoutPath(TXTImport->text());
        ImportFolder = Path(TXTImport->text());
        if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
        BTNexecImport->setEnabled(true);
        TXTProjectFolder->setEnabled(true);
    }
    if (TXTProjectName->text().isEmpty() && !TXTImport->text().isEmpty()){
        TXTProjectName->setText(FileWithoutExt(TXTImport->text()));
        TXTProjectName->setEnabled(true);
    }
    page_atdens->setEnabled(false);
}

/* Textbox for input data file name */
void MainWindow::TXTProjectFolder_changed(const QString &cad)
{
    ProjectFolder = cad;
    if (cad.mid(cad.length()-1,1) != "/"){
            ProjectFolder.append("/");
    }
}

/* Textbox with project name */
void MainWindow::TXTProjectName_changed(const QString &cad)
{
    if(TXTProjectName->text().isEmpty()){
        BTNexecImport->setEnabled(false);
    }
    else{
        if(!TXTImport->text().isEmpty()) BTNexecImport->setEnabled(true);
    }
    ProjectName = cad;
    QString NombreArchivo = ProjectName;
    NombreArchivo.append(".damproj");
    SetCurrentFile(NombreArchivo,true,true);
    TXTdensdamdenfilename->setText(ProjectName);
    TXTHFgdamforcesfilename->setText(ProjectName);
    TXTeffilename->setText(ProjectName);
    TXTdgfilename->setText(ProjectName);
    TXTpotgdampotfilename->setText(ProjectName);
    TXTtopofilename->setText(ProjectName);
    TXTfraddamfilename->setText(ProjectName);
    TXTmrotorimultfilename->setText(ProjectName);
    TXTMOfilename->setText(ProjectName);
    TXTdZJfilename->setText(ProjectName);
}

/* Validates file name */
void MainWindow::TXTValidate_changed()
{
}

/**************************************************************************************************/
/********************** DIRECT IMPORT OF FILES WITH GEOMETRY; BASIS SET AND DENSITY  **************/
/**************************************************************************************************/

/* Gaussian data files  */
void MainWindow::execGgbsDen()
{
    QDir path(ProjectFolder);
    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
            return;
        }
    }
    if (ProjectFolder.at(ProjectFolder.length()-1) != '/')
        ProjectFolder.append('/');
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    if (ImportFolder.at(ImportFolder.length()-1) != '/')
        ImportFolder.append('/');
    QString DirNombreImport = ImportFolder+ImportFile;
    if (ProjectFolder != ImportFolder){
        bool ok = false;
        ok = QFile::copy(DirNombreImport,ProjectFolder+ImportFile);
        if (ok){
            QString denfilename = QFileInfo(ImportFile).completeBaseName()+".den";
            ok = QFile::copy(ImportFolder+denfilename,ProjectFolder+denfilename);
            if (!ok){
                denfilename = denfilename+ ".gz";
                ok = QFile::copy(ImportFolder+denfilename,ProjectFolder+denfilename);
                if (!ok){
                    QMessageBox::warning(this,tr("DAMQT"),tr("Cannot copy file %1(.gz) to folder %3")
                            .arg(denfilename).arg(ProjectFolder));
                    return;
                }
            }
        }
    }
    QString ext;
    QString filepath;
    ext="_2016.damqt";
    QFile file2(Path(DirNombreImport)+"/"+FileWithoutExt(DirNombreImport)+ext);
    if (file2.exists()) {
        filepath=ProjectFolder+ ProjectName+ ext;
        file2.copy(filepath);
    }
    defineRanges();

    page_atdens->setEnabled(true);
    if (file2.exists())
    {
        page_MED->setEnabled(true);
        page_MESP->setEnabled(true);
        page_TOPO->setEnabled(true);
        page_HFforces->setEnabled(true);
        page_Efield->setEnabled(true);
        page_densgrad->setEnabled(true);
        page_frad->setEnabled(true);
        page_orimult->setEnabled(true);
        page_ZJdens->setEnabled(true);
        page_ZJtab->setEnabled(true);
    }
    page_MO->setEnabled(true);
}

/* Slater data files  */
void MainWindow::execSxyzDen()
{
    QDir path(ProjectFolder);
    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
            return;
        }
    }
    if (ProjectFolder.at(ProjectFolder.length()-1) != '/')
        ProjectFolder.append('/');
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    if (ImportFolder.at(ImportFolder.length()-1) != '/')
        ImportFolder.append('/');
    QString DirNombreImport = ImportFolder+ImportFile;
    bool ok = false;
    if (ProjectFolder != ImportFolder){
        bool ok = false;
        ok = QFile::copy(DirNombreImport,ProjectFolder+ImportFile);
        if (ok && QFileInfo(ImportFile).suffix() == "sgbs"){
            QString denfilename = QFileInfo(ImportFile).completeBaseName()+".den";
            ok = QFile::copy(ImportFolder+denfilename,ProjectFolder+denfilename);
            if (!ok){
                denfilename = denfilename+ ".gz";
                ok = QFile::copy(ImportFolder+denfilename,ProjectFolder+denfilename);
                if (!ok){
                    QMessageBox::warning(this,tr("DAMQT"),tr("Cannot copy file %1(.gz) to folder %3")
                            .arg(denfilename).arg(ProjectFolder));
                    return;
                }
            }
        }
    }
    QString ext;
    QString filepath;
    ext="_2016.damqt";
    QFile file2(Path(DirNombreImport)+"/"+FileWithoutExt(DirNombreImport)+ext);
    if (file2.exists())
    {
        filepath=ProjectFolder+ ProjectName+ ext;
        file2.copy(filepath);
    }
    defineRanges();
    page_atdens->setEnabled(true);
    if (file2.exists())
    {
        page_MED->setEnabled(true);
        page_MESP->setEnabled(true);
        page_TOPO->setEnabled(true);
        page_HFforces->setEnabled(true);
        page_Efield->setEnabled(true);
        page_densgrad->setEnabled(true);
        page_frad->setEnabled(true);
        page_orimult->setEnabled(true);
        page_ZJdens->setEnabled(true);
        page_ZJtab->setEnabled(true);
    }
    page_MO->setEnabled(true);
}
    
/**************************************************************************************************/
/********************** INTERFACES TO OTHER PROGRAMS FOR DATA IMPORT  *****************************/
/**************************************************************************************************/

/* Reads data from a GAUSSIAN fchk file */
void MainWindow::execFchk()
{
    if (ImportFolder.at(ImportFolder.length()-1) != '/')
        ImportFolder.append('/');
    QString DirNombreImport = ImportFolder+ImportFile;
    QFile file(DirNombreImport);
    QDir path(ProjectFolder);
    QString filepath=ProjectFolder + FileWithoutPath(DirNombreImport); // + ".fchk";
    if (mden==0){
        if (!path.exists(ProjectFolder)) {
            QMessageBox msgBox;
            msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
            msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
            msgBox.setButtonText(QMessageBox::No, tr("No"));
            msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
            msgBox.setIcon(QMessageBox::Warning);
            int ret = msgBox.exec();
            if (ret == QMessageBox::Yes){
                createDir(ProjectFolder);
                statusBar()->showMessage(tr("Project succesfully created"), 2000);
            }else{
                return;
            }
        }
    }
    QStringList Parameters;
    QString q=QString("%1").arg(mden);
    Parameters << FileWithoutPath(filepath) << q << ImportFolder << ProjectFolder << ProjectName ;
    QString strprocess;
    QString processname = "GAUSS_interface.exe";
    QString execName = get_execName(processname, QString("interfaces_320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(create_damproj(int, QProcess::ExitStatus)));
    executing = 0;
    myProcess->start(strprocess,Parameters);
}

void MainWindow::create_damproj(int exitCode, QProcess::ExitStatus exitStatus){
    if(exitStatus==QProcess::NormalExit && mden == 0){
        QString damprojFile = ProjectFolder+ProjectName+".damproj";
        if (!QFile::exists(damprojFile))  {
            loadDefault(1);
            saveOptions(damprojFile,0);
        }
    }
}

/* Reads data from an MOLEKEL .mkl  file */
void MainWindow::execMOLEKEL()
{
    QDir path(ProjectFolder);

    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
            return;
        }
    }
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    QStringList Parameters;
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    Parameters << QFileInfo(ImportFile).completeBaseName() << ImportFolder  << ProjectFolder << ProjectName  ;
    QString strprocess;
    QString processname = "MOLEKEL_interface.exe";
    QString execName = get_execName(processname, QString("interfaces_320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 2;
    myProcess->start(strprocess,Parameters);
}

/* Reads data from a MOLPRO out file */
void MainWindow::execMolpro()
{
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    QString DirNombreImport = ImportFolder+ImportFile;
    QString suffix=Extension(DirNombreImport);
    QFile file(DirNombreImport);
    QDir path(ProjectFolder);
    QString filepath=ProjectFolder + FileWithoutPath(DirNombreImport);
    QStringList Parameters;
    if (suffix=="out"){
        QString stdoutput = ProjectFolder + ProjectName + "-MOLPRO_out_interface.out";;
        if (mden==0){
            if (!path.exists(ProjectFolder)) {
                QMessageBox msgBox;
                msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
                msgBox.setDefaultButton(QMessageBox::Cancel);
                msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
                msgBox.setButtonText(QMessageBox::No, tr("No"));
                msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
                msgBox.setIcon(QMessageBox::Warning);
                int ret = msgBox.exec();
                if (ret == QMessageBox::Yes){
                    createDir(ProjectFolder);
                    statusBar()->showMessage(tr("Project succesfully created"), 2000);
                }else{
                    return;
                }
            }
            QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
            existsinp(DirNombreArchivo,0,1,true);

            file.copy(DirNombreImport,filepath);
            int ierror=file.error();
            while (ierror!=0)
            {
                ierror=file.error();
            }
        }
        QString q=QString("%1").arg(mden);
        Parameters << ImportFile << q << ImportFolder << ProjectFolder << ProjectName ;
        QString strprocess;
        QString processname = "MOLPRO_out_interface.exe";
        QString execName = get_execName(processname, QString("interfaces_320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        executing = 3;
        myProcess = new QProcess(this);
        myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
        connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
        connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
        connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
        myProcess->start(strprocess,Parameters);
    }
    else if(suffix=="xml"){
        if (!path.exists(ProjectFolder)) {
            QMessageBox msgBox;
            msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
            msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
            msgBox.setButtonText(QMessageBox::No, tr("No"));
            msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
            msgBox.setIcon(QMessageBox::Warning);
            int ret = msgBox.exec();
            if (ret == QMessageBox::Yes){
                createDir(ProjectFolder);
                statusBar()->showMessage(tr("Project succesfully created"), 2000);
            }else{
                return;
            }
        }
        QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
        existsinp(DirNombreArchivo,0,1,true);
        file.copy(DirNombreImport,filepath);
        int ierror=file.error();
        while (ierror!=0)
        {
            ierror=file.error();
        }
        QString q=QString("%1").arg(mden);
        QString scriptFile = QCoreApplication::applicationDirPath() + "/MOLPRO_xml_interface.py" ;
        QStringList pythonCommandArguments = QStringList() << scriptFile << ImportFile << ImportFolder
                    << ProjectFolder << ProjectName ;
        QString strprocess;
        QString execName = get_python();
        if (execName.isEmpty())
            return;
        strprocess = execName;
        executing = 4;
        myProcess = new QProcess(this);
        connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
        connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
        connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
        myProcess->start(strprocess,pythonCommandArguments);
    }
}

/* Reads data from a MOPAC aux file */
void MainWindow::execMopac()
{
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    QString DirNombreImport = ImportFolder+ImportFile;
    QString suffix=Extension(DirNombreImport);
    QFile file(DirNombreImport);
    QDir path(ProjectFolder);
    QString filepath=ProjectFolder + FileWithoutPath(DirNombreImport);
    QStringList Parameters;
    QString strprocess;
    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
            return;
        }
    }
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    Parameters << QFileInfo(ImportFile).completeBaseName() << ImportFolder  << ProjectFolder << ProjectName  ;
    QString processname = "Mopac_aux_interface.exe";
    QString execName = get_execName(processname, QString("interfaces_320"));
    if (execName.isEmpty())
        return;

    strprocess = execName;
    executing = 5;
//    QString stdoutput = ProjectFolder + ProjectName + "-mopac.out";
    myProcess = new QProcess(this);
//    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    myProcess->start(strprocess,Parameters);
}

/* Reads data from a TURBOMOLE coord, mos and basis  files */
void MainWindow::execTurbom()
{
    QDir path(ProjectFolder);

    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
                return;
        }
    }
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    QStringList Parameters;
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    Parameters << QFileInfo(ImportFile).completeBaseName() << ImportFolder  << ProjectFolder << ProjectName  ;
    QString strprocess;
    QString processname = "TURBOMOLE_interface.exe";
    QString execName = get_execName(processname, QString("interfaces_320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 1;
    QByteArray qbarray = myProcess->readAllStandardError();
    myProcess->start(strprocess,Parameters);
}

/* Reads data from a NWChem output file (IMPORTANT! extension must be nwcout) */
void MainWindow::execNWChem()
{
    QDir path(ProjectFolder);

    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
                return;
        }
    }
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    QStringList Parameters;
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    Parameters << QFileInfo(ImportFile).completeBaseName() << ImportFolder  << ProjectFolder << ProjectName  ;
    QString strprocess;
    QString processname = "NWChem_interface.exe";
    QString execName = get_execName(processname, QString("interfaces_320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 7;
    QByteArray qbarray = myProcess->readAllStandardError();
    myProcess->start(strprocess,Parameters);
}


/* Reads data from Psi4 auxiliary density files (IMPORTANT! extension must be psiauxden) */
void MainWindow::execPsi4auxiliary()
{
    QDir path(ProjectFolder);

    if (!path.exists(ProjectFolder)) {
        QMessageBox msgBox;
        msgBox.setInformativeText(QString(tr("Project %1 not found")).arg(ProjectFolder)+"\n"+tr("Do you wish to create?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
        msgBox.setIcon(QMessageBox::Warning);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes){
            createDir(ProjectFolder);
            statusBar()->showMessage(tr("Project succesfully created"), 2000);
        }else{
                return;
        }
    }
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,0,1,true);
    QStringList Parameters;
    if (ImportFolder.at(ImportFolder.length()-1) != '/') ImportFolder.append('/');
    Parameters << QFileInfo(ImportFile).completeBaseName() << ImportFolder  << ProjectFolder << ProjectName  ;
    QString strprocess;
    QString processname = "Psi4_auxiliary_interface.exe";
    QString execName = get_execName(processname, QString("interfaces_320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    myProcess = new QProcess(this);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 8;
    QByteArray qbarray = myProcess->readAllStandardError();
    myProcess->start(strprocess,Parameters);
}

/**************************************************************************************************/
/********************************  CONTROL OF PROCESSES *******************************************/
/**************************************************************************************************/

//    Returns the error code of the external process currently running
void MainWindow::processError(QProcess::ProcessError error)
{
    QString strprocess;
    strprocess = Who_executing(executing);
//    if (executing == 4  || executing == 5 )
    if (executing == 4)
        strprocess = strprocess + ".py";  // Python interface for MOLPRO xml files
    else
        strprocess = strprocess + ".exe";

    if(error==QProcess::FailedToStart){
        QStringList strlistproc = strprocess.split("/");
        QString message;
        if (strlistproc.length() >= 2){
            message = QString("Error %1").arg(error)
                    + QString(tr("Process failed to start program %1\n").arg(strprocess))
                    + QString(tr("Check that the program is installed in any of the following directories: \n %1 \n %2 \n")
                                  .arg(QCoreApplication::applicationDirPath()).arg(strlistproc.at(strlistproc.length()-2)))
                    + QString(tr("or in any other directory in your $PATH"));
        }
        else{
            message = QString("Error %1").arg(error)
                    + QString(tr("Process failed to start program %1\n").arg(strprocess))
                    + QString(tr("Check that the program is installed in the directory: \n %1 \n")
                                  .arg(QCoreApplication::applicationDirPath()))
                    + QString(tr("or in any other directory in your $PATH"));
        }
        QMessageBox::critical(this,QString("Error %1").arg(error),message);
    }
    else{
        QMessageBox::critical(this,QString("Error %1").arg(error),tr("Error when running program %1").arg(strprocess));
    }
}

//    Slot to be run when a process ends
void MainWindow::processOutput(int exitCode, QProcess::ExitStatus exitStatus)
{
//    Processes numbering:  0: execFchk ; 1: execTurbom ; 2: execMOLEKEL ; 3: execMolpro (.out) ; 4: execMolpro (.xml)
//           5: execMopac ; 6: execsgbs2sxyz ; 7: execNWChem ; 8: execPsi4auxiliary ; 9: void ; 10: void ;
//          11: execDam ; 12: execDamden ; 13: execDampot ; 14: execDamforces ; 15: execDamfield
//          16: execDamfrad ; 17: execDammultrot ; 18: execDamorb ; 19: execDamTopography
//          20: execDamZJ ; 21: execDamdenZJ ; 22: execDamdengrad ; 23: execDamSGhole
    if (executing < 0 || executing > 23)
        return;
    if (QFileInfo(ProjectFolder+"zdo").exists()){
        lzdo = true;
    }
    else{
        lzdo = false;
    }
    if (QFileInfo(ProjectFolder+"valence").exists()){
        lvalence = true;
    }
    else{
        lvalence = false;
    }
    QString strprocess,str;
    str=Who_executing(executing);
    strprocess = str + ".exe";
    str = "-" + str ;
    if(exitStatus==QProcess::NormalExit){
        QString fileName=ProjectFolder+ProjectName+str;    // Default name of files
        if (executing == 12){
            if (!TXTdensdamdenfilename->text().isEmpty())
                fileName=ProjectFolder+TXTdensdamdenfilename->text()+str;
            if (RBTdens2D->isChecked() && RBTdensplane->isChecked())
                rename_density_cntfile();
        }
        if (executing == 13){
            if (!TXTpotgdampotfilename->text().isEmpty())
                fileName=ProjectFolder+TXTpotgdampotfilename->text()+str;
            if (RBTpot2D->isChecked() && RBTpotplane->isChecked())
                rename_pot_cntfile();
        }
        else if (executing == 14 && !TXTHFgdamforcesfilename->text().isEmpty())
            fileName=ProjectFolder+TXTHFgdamforcesfilename->text()+str;
        else if (executing == 15 && !TXTeffilename->text().isEmpty())
            fileName=ProjectFolder+TXTeffilename->text()+str;
        else if (executing == 16 && !TXTfraddamfilename->text().isEmpty())
            fileName=ProjectFolder+TXTfraddamfilename->text()+str;
        else if (executing == 17 && !TXTmrotorimultfilename->text().isEmpty())
            fileName=ProjectFolder+TXTmrotorimultfilename->text()+str;
        else if (executing == 18 && !TXTMOfilename->text().isEmpty())
            fileName=ProjectFolder+TXTMOfilename->text()+str;
        else if (executing == 19 && !TXTtopofilename->text().isEmpty()){
            fileName=ProjectFolder+TXTtopofilename->text()+str;
        }
        else if (executing == 21 && !TXTdZJfilename->text().isEmpty())
            fileName=ProjectFolder+TXTdZJfilename->text()+str;
        else if (executing == 22 && !TXTdgfilename->text().isEmpty())
            fileName=ProjectFolder+TXTdgfilename->text()+str;
        else if (executing == 23 && !TXTSGholefilename->text().isEmpty()){
            fileName=ProjectFolder+TXTSGholefilename->text()+str;
        }
        if ((CHKatdensmpi->isChecked()   && (executing == 11)) || 
                (RBTdens3D->isChecked() && CHKdensmpi->isChecked() && CHKdensgrid->isChecked() && (executing == 12)) ||
                (RBTpot3D->isChecked() && CHKpotmpi->isChecked() && CHKpotgrid->isChecked() && (executing == 13)) ||
                (RBTef3D->isChecked() && CHKefmpi->isChecked() && (executing == 15)) ||
                (RBTMO3D->isChecked() && CHKMOmpi->isChecked() && CHKMOgrid->isChecked() && (executing == 18)) ||
                (CHKtopompi->isChecked() && (executing == 19)) ||
                (CHKZJmpi->isChecked()   && (executing == 20)) ||
                (RBTdZJ3D->isChecked() && CHKdZJmpi->isChecked() && CHKdZJgrid->isChecked() && (executing == 21)) ||
                (RBTdg3D->isChecked() && CHKdgmpi->isChecked() && (executing == 22)) ||
                (CHKSGholempi->isChecked() && (executing == 23)) )
                fileName.append("_mpi");
        if (executing == 12){
            if(RBTdensExact->isChecked()){
                fileName.insert(fileName.indexOf("-DAMDEN"),"_exact");
            }
        }
        if (executing == 13){
            if(CHKpotexact->isChecked()){
                fileName.insert(fileName.indexOf("-DAMPOT"),"_exact");
            }
        }
        if (executing == 19){
            if(RBTtopodensity->isChecked())
                fileName.append("-d");
            else {
                fileName.append("-v");
            }
        }
        if (executing == 23){
            if(CHKSGhexactMESP->isChecked()){
                fileName.insert(fileName.indexOf("-DAMSGHOLE"),"_exact");
            }
        }
        fileName.append(".out");
        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("ProcessOutput"),tr("File %1 cannot be read")
                    .arg(fileName)+QString(":\n%1.").arg(file.errorString()));
        }else{
            if (executing == 0 || executing == 1 || executing == 2 || executing == 3 || executing == 4
                     || executing == 5 || executing == 6 || executing == 7){
                page_atdens->setEnabled(true);
            }else if (executing == 11){
                page_MED->setEnabled(true);
                page_MESP->setEnabled(true);
                page_TOPO->setEnabled(true);
                page_HFforces->setEnabled(true);
                page_Efield->setEnabled(true);
                page_densgrad->setEnabled(true);
                page_frad->setEnabled(true);
                page_orimult->setEnabled(true);
                page_ZJdens->setEnabled(true);
                page_ZJtab->setEnabled(true);
                page_SGhole->setEnabled(true);
            }
            page_MO->setEnabled(true);
            QTextStream in(&file);
            textEdit->setFont(QFont("Courier",10));
            textEdit->setPlainText(in.readAll());
        }
        if (executing == 5){
            lzdo = true;
            lvalence = true;
            QFile filezdo(ProjectFolder+"zdo");
            filezdo.open(QFile::WriteOnly | QFile::Text);
            filezdo.close();
            QFile filevalence(ProjectFolder+"valence");
            filevalence.open(QFile::WriteOnly | QFile::Text);
            filevalence.close();
        }
        statusBar()->showMessage(tr("End of calculation"));
    }else if (exitStatus==QProcess::CrashExit){
        statusBar()->showMessage(tr(QString("Process crashed, exit code = %1").arg(exitCode).toLatin1()));
    }
//    Special cases in which interface must be run twice because there is more than one density matrix available
    if (executing == 0){    // Executing GAUSS_interface: two-pass procedure
        QString fileName=ProjectFolder+ProjectName+str+".out";
        if (!GAUSS_two_pass_case(fileName)){
            mden = 0;
            BTNexecImport->setEnabled(true);
            return;
        };
        if(mden>0){
            QApplication::setOverrideCursor(QPixmap(":/images/wait.png"));
            execFchk();
            QApplication::restoreOverrideCursor();
            return;
        }else{
            defineRanges();
        }
        mden = 0;
    }
    if (executing == 3){    // Executing MOLPRO_out_interface: 
        QString fileName=ProjectFolder+ProjectName+str+".out";
        MOLPRO_two_pass_case(fileName);
        if(mden>0){
            QApplication::setOverrideCursor(QPixmap(":/images/wait.png"));
            execMolpro();
            QApplication::restoreOverrideCursor();
            return;
        }else{
            defineRanges();
        }
        mden = 0;
    }
//    End of special cases
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    SetCurrentFile(DirNombreArchivo,true,false);
    QApplication::restoreOverrideCursor();
    BTNexecImport->setEnabled(true);
    BTNexecDam->setEnabled(true);
    BTNexecDamden->setEnabled(true);
    BTNexecDampot->setEnabled(true);
    BTNexecDamforces->setEnabled(true);
    BTNexecDamfield->setEnabled(true);
    BTNexecDamdengrad->setEnabled(true);
    BTNexecDamfrad->setEnabled(true);
    BTNexecDammultrot->setEnabled(true);
    BTNexecDamorb->setEnabled(true);
    BTNexecDamtopo->setEnabled(true);
    BTNexecDamZJ->setEnabled(true);
    BTNexecDamdZJ->setEnabled(true);
    BTNexecDamSGhole->setEnabled(true);
    if (executing > 10) BTNstop[executing-11]->setEnabled(false);
    executing = -1;
}

//    Dialog for GAUSS special case
bool MainWindow::GAUSS_two_pass_case(QString fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("GAUSS_two_pass_case"),tr("File %1 cannot be read")
                .arg(fileName)+QString(":\n%1.").arg(file.errorString()));
    }
    QTextStream in(&file);
    QString line = in.readLine();
    QString line1 = in.readLine();
    QString line2 = in.readLine();
    QString line3 = in.readLine();
    QString line4 = in.readLine();
    QString line5 = in.readLine();
    QString line6 = in.readLine();
    QString line7 = in.readLine();
    file.close();

    dialog = new Dialog(8);
    (dialog->RBToption[0])->setText(tr("SCF Total density matrix"));
    (dialog->RBToption[0])->setChecked(true);
    dialog->val=1;
    (dialog->RBToption[1])->setText(tr("SCF Spin density matrix"));
    (dialog->RBToption[2])->setText(tr("CI Total density matrix"));
    (dialog->RBToption[3])->setText(tr("CI Spin density matrix"));
    (dialog->RBToption[4])->setText(tr("MP2 Total density matrix"));
    (dialog->RBToption[5])->setText(tr("MP2 Spin density matrix"));
    (dialog->RBToption[6])->setText(tr("CC Total density matrix"));
    (dialog->RBToption[7])->setText(tr("CC Spin density matrix"));
    dialog->Implement();
    if(line=="existSCF" || line1=="existSCF" || line2=="existSCF" || line3=="existSCF"
           || line4=="existSCF" || line5=="existSCF"){
        dialog->buttons[0]=true;
    }else{
        dialog->buttons[0]=false;
    }
    if(line=="existSpSCF" || line1=="existSpSCF" || line2=="existSpSCF" || line3=="existSpSCF"
           || line4=="existSpSCF" || line5=="existSpSCF"){
        dialog->buttons[1]=true;
    }else{
        dialog->buttons[1]=false;
    }
    if(line=="existCI" || line1=="existCI" || line2=="existCI" || line3=="existCI"
           || line4=="existCI" || line5=="existCI"){
        dialog->buttons[2]=true;
    }else{
        dialog->buttons[2]=false;
    }
    if(line=="existSpCI" || line1=="existSpCI" || line2=="existSpCI" || line3=="existSpCI"
           || line4=="existSpCI" || line5=="existSpCI"){
        dialog->buttons[3]=true;
    }else{
        dialog->buttons[3]=false;
    }
    if(line=="existCI" || line1=="existCI" || line2=="existCI" || line3=="existCI"
           || line4=="existCI" || line5=="existCI"){
        dialog->buttons[2]=true;
    }else{
        dialog->buttons[2]=false;
    }
    if(line=="existSpCI" || line1=="existSpCI" || line2=="existSpCI" || line3=="existSpCI"
           || line4=="existSpCI" || line5=="existSpCI"){
        dialog->buttons[3]=true;
    }else{
        dialog->buttons[3]=false;
    }
    if(line=="existMP2" || line1=="existMP2" || line2=="existMP2" || line3=="existMP2"
           || line4=="existMP2" || line5=="existMP2"){
        dialog->buttons[4]=true;
    }else{
        dialog->buttons[4]=false;
    }
    if(line=="existSpMP2" || line1=="existSpMP2" || line2=="existSpMP2" || line3=="existSpMP2"
           || line4=="existSpMP2" || line5=="existSpMP2"){
        dialog->buttons[5]=true;
    }else{
        dialog->buttons[5]=false;
    }
    if(line=="existCC" || line1=="existCC" || line2=="existCC" || line3=="existCC"
           || line4=="existCC" || line5=="existCC" || line6=="existCC" || line7=="existCC"){
        dialog->buttons[6]=true;
    }else{
        dialog->buttons[6]=false;
    }
    if(line=="existSpCC" || line1=="existSpCC" || line2=="existSpCC" || line3=="existSpCC"
           || line4=="existSpCC" || line5=="existSpCC" || line6=="existSpCC" || line7=="existSpCC"){
        dialog->buttons[7]=true;
    }else{
        dialog->buttons[7]=false;
    }
    dialog->RadioChange();
    QApplication::restoreOverrideCursor();
    if(dialog->buttons[0]==false && dialog->buttons[1]==false && dialog->buttons[2]==false && dialog->buttons[3]==false
           && dialog->buttons[4]==false && dialog->buttons[5]==false
           && dialog->buttons[6]==false && dialog->buttons[7]==false){
        mden = 0;
        return true;
    }else if(dialog->buttons[0]==true && dialog->buttons[1]==false && dialog->buttons[2]==false && dialog->buttons[3]==false
           && dialog->buttons[4]==false && dialog->buttons[5]==false
           && dialog->buttons[6]==false && dialog->buttons[7]==false  ){
        mden = 1;
        return true;
    }else{
        dialog->exec();
        if (dialog->cancel){
            return false;
        }
        QStringList list;
        list << "SCF" << "SP_SCF" << "CI" << "SP_CI" << "MP2" << "SP_MP2" << "CC" << "SP_CC";
        mden=dialog->val;
        QString text = TXTProjectName->text();
        for (int i = 0 ; i < list.length() ; i++){
            if (text.contains("-"+list.at(i))){
                text.replace("-"+list.at(i),"");
            }
        }
        TXTProjectName->setText(text+"-"+list.at(mden-1));
        return true;
    }
}

//    Dialog for MOLPRO special case
void MainWindow::MOLPRO_two_pass_case(QString fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MOLPRO_two_pass_case"),tr("File %1 cannot be read")
                .arg(fileName)+QString(":\n%1.").arg(file.errorString()));
    }
    QTextStream in(&file);
    QString line = in.readLine();
    int nentries = line.mid(0,2).toInt();

    if (nentries > 0){
        dialog = new Dialog(nentries);
        for (int i = 0 ; i < nentries ; ++i){
            dialog->buttons[i]=true;
            line = in.readLine();
            dialog->RBToption[i]->setText(line.mid(2));
        }
        dialog->RBToption[0]->setChecked(true);
        dialog->val=1;
        dialog->Implement();
        dialog->RadioChange();
        QApplication::restoreOverrideCursor();
        if (nentries == 1){
            mden = 1;
        }else{
            dialog->exec();
            mden=dialog->val;
        }
    }
    else
            mden = 0;
    file.close();
}

//    Starts an external process
void MainWindow::processStart()
{
    BTNexecImport->setEnabled(false);
    BTNexecDam->setEnabled(false);
    BTNexecDamden->setEnabled(false);
    BTNexecDampot->setEnabled(false);
    BTNexecDamforces->setEnabled(false);
    BTNexecDamfield->setEnabled(false);
    BTNexecDamdengrad->setEnabled(false);
    BTNexecDamfrad->setEnabled(false);
    BTNexecDammultrot->setEnabled(false);
    BTNexecDamorb->setEnabled(false);
    BTNexecDamtopo->setEnabled(false);
    BTNexecDamZJ->setEnabled(false);
    BTNexecDamdZJ->setEnabled(false);
    BTNexecDamSGhole->setEnabled(false);
    if (executing > 10) BTNstop[executing-11]->setEnabled(true);

    statusBar()->showMessage(tr("Computing..."));
}


//    Stops the external process currently running
void MainWindow::processStop()
{
    #if defined(Q_WS_WIN) || defined(Q_OS_WIN) || QT_VERSION < 0x050000
        myProcess->kill();
    #else
        if (myProcess->arguments().count() > 2){
            QString processname = myProcess->arguments().at(2);
            if (processname.contains("_mpi")){      // In case of mpi runs, kill all associated processes
//                processname = processname.split('.exe').first();
                QProcess getprocesses;
                QString pgrep;
                pgrep = QString("pgrep -f %1").arg(processname);
                qDebug() << "pgrep = " << pgrep;
                getprocesses.start(pgrep);
                getprocesses.waitForFinished();
                QByteArray procnumbers = getprocesses.readAllStandardOutput();
                QString procstring(procnumbers);
                QStringList procslist = procstring.split("\n");
                QMessageBox msgBox;
                msgBox.setInformativeText(QString(tr("Do you want to kill all processes named %1?").arg(processname)+"\n"
                       + tr("The following processes will be killed: "))+procstring.replace("\n"," "));
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
                msgBox.setDefaultButton(QMessageBox::Cancel);
                msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
                msgBox.setButtonText(QMessageBox::No, tr("No"));
                msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
                msgBox.setIcon(QMessageBox::Warning);
                int ret = msgBox.exec();
                if (ret == QMessageBox::Yes){
                    for (int i = 1; i < procslist.length()-1 ; i++){  // Kill the processes one by one (the last element in procslist is empty)
                        QProcess killprocess;
                        QString strprocess;
                        strprocess = QString("kill -9 %1").arg(procslist.at(i));
                        killprocess.start(strprocess);
                        killprocess.waitForFinished();
                    }
                    myProcess->kill();
                }
            }
        }
        else{
            myProcess->kill();
        }
    #endif
}

// Determines which external program is being executed
// Processes corresponding to Fortran runs start on index 11
QString MainWindow::Who_executing(int caso)
{
    QString str;

    if (caso==0){
        str="GAUSS_interface";
    }else if (caso == 1){
        str="TURBOMOLE_interface";
    }else if (caso == 2){
        str="MOLEKEL_interface";
    }else if (caso == 3){
        str="MOLPRO_out_interface"; // From out file
    }else if (caso == 4){
        str="MOLPRO_xml_interface"; // From xml file
    }else if (caso == 5){
        str="MOPAC_aux_interface";
    }else if (caso == 6){
        str="sgbs2sxyz";
    }else if (caso == 7){
        str="NWChem_interface";
    }else if (caso == 8){
        str="Psi4_auxiliary_interface";
    }else if (caso == 11){
        if (lslater) 
            str="DAMSTO320";
        else
            str="DAMGTO320";
    }else if (caso == 12){
        str="DAMDEN320";
    }else if (caso == 13){
        str="DAMPOT320";
    }else if (caso == 14){
        str="DAMFORCES320";
    }else if (caso == 15){
        str="DAMFIELD320";
    }else if (caso == 16){
        str="DAMFRAD320";
    }else if (caso == 17){
        str="DAMMULTROT320";
    }else if (caso == 18){
        str="DAMORB320";
    }else if (caso == 19){
        str="DAMTOPO320";
    }else if (caso == 20){
        if (RBTZJacobi->isChecked()){
            str += "Jacobi-DAMZJ320";
        }
        else{
            str += "Zernike-DAMZJ320";
        }
    }else if (caso == 21){
        if (ldZJjacobi){
            str += "jacobi-DAMDENZJ320";
        }
        else{
            str += "zernike-DAMDENZJ320";
        }
    }else if (caso == 22){
        str="DAMDENGRAD320";
    }else if (caso == 23){
        str="DAMSGHOLE320";
    }
    return str;
}


/**************************************************************************************************/
/**********************  FUNCTIONS FOR READING AND WRITING OPTIONS IN PROJECT FILE  ***************/
/**************************************************************************************************/


/*    LOADS DEFAULT VALUES in the widgets of the mainmenu        */
void MainWindow::loadDefault(int all)
{
    if (all==0){      
        loadDefault_project();     //    page_project: Project
    }
//    textEdit->clear();
    loadDefault_atdens();
    loadDefault_MED();
    loadDefault_MESP();
    loadDefault_SGhole();
    loadDefault_TOPO();
    loadDefault_HFforces();
    loadDefault_Efield();
    loadDefault_densgrad();
    loadDefault_frad();
    loadDefault_orimult();
    loadDefault_MO();
    loadDefault_ZJdens();
    loadDefault_ZJtab();

//      End of defaults loading

    SetCurrentFile("",true,false);
}

void MainWindow::loadDefault_project(){       //    page_project: Project
        TXTImport->setText("");
        TXTImport->setEnabled(true);
        BTNImport->setEnabled(true);

        TXTProjectFolder->setText("");
        TXTProjectFolder->setEnabled(false);

        TXTProjectName->setText("");
        TXTProjectName->setEnabled(true);

        BTNexecImport->setEnabled(false);
        set_natom(0);
}

void MainWindow::loadDefault_atdens(){   //    page_atdens: Atomic densities
    SPBatdenslmaxexp->setRange(0,MAX_LEXP);
    SPBatdenslmaxexp->setValue(10);
    RBTatdensDtotal->setChecked(true);
}

void MainWindow::loadDefault_MED(){   //    page_MED: Density
    if (TXTProjectName){
            TXTdensdamdenfilename->setText(TXTProjectName->text());
    }
    else{
            TXTdensdamdenfilename->setText("");
    }
    RBTdensfulldensity->setChecked(true);
    SPBdenslmaxexp->setRange(0,MAX_LEXP);
    SPBdenslmaxexp->setValue(10);
    SPBdenslminexp->setRange(0,MAX_LEXP);
    SPBdenslminexp->setValue(0);
    RBTdensRep1->setChecked(true);
    CHKdenslmolec->setChecked(true);
    CHKdenslatomics->setChecked(false);
    CHKdensldensacc->setChecked(false);
    TXTdensxinf->setText("-4.0");
    TXTdensyinf->setText("-4.0");
    TXTdenszinf->setText("-4.0");
    TXTdensuinf->setText("-4.0");
    TXTdensvinf->setText("-4.0");
    TXTdensxsup->setText("4.0");
    TXTdensysup->setText("4.0");
    TXTdenszsup->setText("4.0");
    TXTdensusup->setText("4.0");
    TXTdensvsup->setText("4.0");
    TXTdensxformula2D->setText("u");
    TXTdensyformula2D->setText("v");
    TXTdenszformula2D->setText("0");
    RBTdensrlow->setChecked(true);
    RBTdens3D->setChecked(true);
    RBTdensrlow->setToolTip("65x65x65");
    RBTdensrmedium->setToolTip("129x129x129");
    RBTdensrhigh->setToolTip("257x257x257");
    SHTxyz->clear();
    CHKdensgrid->setChecked(true);
    Wtableden->setEnabled(false);
}

void MainWindow::loadDefault_MESP(){   //    page_MESP: Electrostatic potential
    if (TXTProjectName){
        TXTpotgdampotfilename->setText(TXTProjectName->text());
    }
    else{
        TXTpotgdampotfilename->setText("");
    }
    SPBpotlmaxexp->setRange(0,25);
    SPBpotlmaxexp->setValue(10);
    CHKpotlong->setChecked(false);
    SPBpotlongthreshold->setRange(-10,0);
    SPBpotlongthreshold->setValue(-7);
    TXTpotxinf->setText("-4.0");
    TXTpotyinf->setText("-4.0");
    TXTpotzinf->setText("-4.0");
    TXTpotuinf->setText("-4.0");
    TXTpotvinf->setText("-4.0");
    TXTpotxsup->setText("4.0");
    TXTpotysup->setText("4.0");
    TXTpotzsup->setText("4.0");
    TXTpotusup->setText("4.0");
    TXTpotvsup->setText("4.0");
    RBTpotrlow->setChecked(true);
    RBTpot3D->setChecked(true);
    RBTpotrlow->setToolTip("65x65x65");
    RBTpotrmedium->setToolTip("129x129x129");
    RBTpotrhigh->setToolTip("257x257x257");
    SHTpotxyz->clear();
    CHKpotgrid->setChecked(true);
    Wtablepot->setEnabled(false);
}

void MainWindow::loadDefault_SGhole(){   //    page_SGhole: HMESP sigma hole
    TXTImportSGholeden->setText("");
    TXTSGholecontour->setText("0.001");
    TXTSGholefilename->setText("");
    SPBSGholelocalextrema->setValue(90);
    SPBSGholegeomthreshold->setValue(-5);
    SPBSGholelongthreshold->setValue(-9);
}
        
void MainWindow::loadDefault_TOPO(){   //    page_TOPO: Molecular topography
    if (TXTProjectName){
        TXTtopofilename->setText(TXTProjectName->text());
    }
    else{
        TXTtopofilename->setText("");
    }
    CHKtopobasin->setChecked(false);
    CHKtopomolgraph->setChecked(false);
    RBTtopodensity->setChecked(true);
    RBTtopopotential->setChecked(false);
    SHTtopoxyz->clear();
    CHKtopoxyz->setEnabled(false);
    CHKtopoxyz->setVisible(false);
    Wtabledentopo->setEnabled(false);
    Wtabledentopo->setVisible(false);
}

void MainWindow::loadDefault_HFforces(){   //    page_HFforces: Hellmann-Feynman forces
    if (TXTProjectName){
        TXTHFgdamforcesfilename->setText(TXTProjectName->text());
    }
    else{
        TXTHFgdamforcesfilename->setText("");
    }
    CHKHFlatomsel->setChecked(false);
}

void MainWindow::loadDefault_Efield(){   //    page_Efield: Electric field
    page_Efield->setEnabled(false);
    if (TXTProjectName){
        TXTeffilename->setText(TXTProjectName->text());
    }
    else{
        TXTeffilename->setText("");
    }
    SPBeflmaxexp->setRange(0,25);
    SPBeflmaxexp->setValue(10);
    CHKeflong->setChecked(false);
    SPBeflongthreshold->setRange(-10,0);
    SPBeflongthreshold->setValue(-9);
    CHKefextralines->setChecked(false);
    TXTefnumpnt->setText("1000");
    TXTefnlinpernuc->setText("16");
    TXTefdlt0->setText("0.02");
    TXTefuinf->setText("-3.0");
    TXTefvinf->setText("-3.0");
    TXTefxinf->setText("-3.0");
    TXTefyinf->setText("-3.0");
    TXTefzinf->setText("-3.0");
    TXTefusup->setText("3.0");
    TXTefvsup->setText("3.0");
    TXTefxsup->setText("3.0");
    TXTefysup->setText("3.0");
    TXTefzsup->setText("3.0");
    TXTeffilelines->setText("");
    LBLeffilelines->setEnabled(false);
    TXTeffilelines->setEnabled(false);
    BTNeffilelines->setEnabled(false);

    CHKefxyz->setEnabled(false);
    CHKefxyz->setVisible(false);
    SHTefxyz->clear();
    Wtabledenef->setEnabled(false);
    Wtabledenef->setVisible(false);
    SHTefuv->clear();
    Wtable2ef->setEnabled(false);
    Wtable2ef->setVisible(false);
}

void MainWindow::loadDefault_densgrad(){   //    page_densgrad: Density gradient
    if (TXTProjectName){
        TXTdgfilename->setText(TXTProjectName->text());
    }
    else{
        TXTdgfilename->setText("");
    }
    SPBdglmaxexp->setRange(0,25);
    SPBdglmaxexp->setValue(10);
    SPBdglongthreshold->setRange(-10,0);
    SPBdglongthreshold->setValue(-9);
    CHKdgextralines->setChecked(false);
    TXTdgnumpnt->setText("1000");
    TXTdgnlinpernuc->setText("16");
    TXTdgdlt0->setText("0.02");
    TXTdguinf->setText("-3.0");
    TXTdgvinf->setText("-3.0");
    TXTdgxinf->setText("-3.0");
    TXTdgyinf->setText("-3.0");
    TXTdgzinf->setText("-3.0");
    TXTdgusup->setText("3.0");
    TXTdgvsup->setText("3.0");
    TXTdgxsup->setText("3.0");
    TXTdgysup->setText("3.0");
    TXTdgzsup->setText("3.0");
    TXTdgfilelines->setText("");
    LBLdgfilelines->setEnabled(false);
    TXTdgfilelines->setEnabled(false);
    BTNdgfilelines->setEnabled(false);

    CHKdgxyz->setEnabled(false);
    CHKdgxyz->setVisible(false);
    SHTdgxyz->clear();
    Wtabledendg->setEnabled(false);
    Wtabledendg->setVisible(false);
    SHTdguv->clear();
    Wtable2dg->setEnabled(false);
    Wtable2dg->setVisible(false);
}

void MainWindow::loadDefault_frad(){   //    page_frad: Radial factors
    if (TXTProjectName){
        TXTfraddamfilename->setText(TXTProjectName->text());
    }
    else{
        TXTfraddamfilename->setText("");
    }
    TXTfradrini->setText("0.0");
    TXTfradrfin->setText("2.0");
    TXTfradltr->setText("0.05");
    CHKfradextras->setChecked(false);
    SHTfradrlist->clear();
    SPBfradltab->setRange(0,25);
    SPBfradltab->setValue(0);
    SPBfradmtab->setRange(0,0);
    SPBfradmtab->setValue(0);
    CHKfradderiv1->setChecked(false);
    CHKfradderiv2->setChecked(false);
}

void MainWindow::loadDefault_orimult(){   //    page_orimult: Oriented multipoles
    if (TXTProjectName){
        TXTmrotorimultfilename->setText(TXTProjectName->text());
    }
    else{
        TXTmrotorimultfilename->setText("");
    }
    SPBmrotlmin->setValue(0);
    SPBmrotlmin->setEnabled(true);
    SPBmrotlmax->setValue(0);
    SPBmrotlmax->setEnabled(true);
    SPBmrotleft->setValue(get_natom());
    spbmrotleft = get_natom();
    SPBmrotleft->setEnabled(true);
    SPBmrotmiddle->setValue(get_natom());
    spbmrotmiddle = get_natom();
    SPBmrotmiddle->setEnabled(true);
    SPBmrotright->setValue(get_natom());
    spbmrotright = get_natom();
    SPBmrotright->setEnabled(true);
}

void MainWindow::loadDefault_MO(){  //    page_MO: Molecular orbitals
    if (TXTProjectName){
            TXTMOfilename->setText(TXTProjectName->text());
    }
    else{
            TXTMOfilename->setText("");
    }
    TXTMOImportfile->setText("");
    TXTMOxinf->setText("-4.0");
    TXTMOyinf->setText("-4.0");
    TXTMOzinf->setText("-4.0");
    TXTMOuinf->setText("-4.0");
    TXTMOvinf->setText("-4.0");
    TXTMOxsup->setText("4.0");
    TXTMOysup->setText("4.0");
    TXTMOzsup->setText("4.0");
    TXTMOusup->setText("4.0");
    TXTMOvsup->setText("4.0");
    RBTMOrlow->setChecked(true);
    RBTMO3D->setChecked(true);
    RBTMOrlow->setToolTip("65x65x65");
    RBTMOrmedium->setToolTip("129x129x129");
    RBTMOrhigh->setToolTip("257x257x257");
}

void MainWindow::loadDefault_ZJdens(){  //    page_ZJdens: One-center Zernike 3D-jacobi expansion
    TXTZJrstar->setText("10.0");
    RBTZJrstarabs->setChecked(true);
    SPBZJlmax->setRange(0,MAX_LEXPZJ);
    SPBZJlmax->setValue(10);
    SPBZJkmax->setValue(10);
    RBTZJechelon->setChecked(false);
    RBTZJZernike->setChecked(true);
    RBTZJacobi->setChecked(false);
    SPBZJnquad->setMinimum(128);
    SPBZJnquad->setMaximum(8192);
    SPBZJnquad->setValue(256);
    SPBZJthrmult->setRange(-15,0);
    SPBZJthrmult->setValue(-10);
    SPBZJthrdist->setRange(-15,0);
    SPBZJthrdist->setValue(-12);
}

void MainWindow::loadDefault_ZJtab(){  //    page_ZJtab: Tabulation of density from Zernike 3D-jacobi expansion
    TXTdZJImportfile->setText("");
    if (TXTProjectName){
        TXTdZJfilename->setText(TXTProjectName->text());
    }
    else{
        TXTdZJfilename->setText("");
    }
    SPBdZJkmax->setValue(SPBZJkmax->value());
    SPBdZJlmax->setRange(0,SPBZJlmax->value());
    SPBdZJlmax->setValue(SPBZJlmax->value());
    SPBdZJlmin->setRange(0,SPBZJlmax->value());
    SPBdZJlmin->setValue(0);
    TXTdensxinf->setText("-10.0");
    TXTdensyinf->setText("-10.0");
    TXTdenszinf->setText("-10.0");
    TXTdensuinf->setText("-10.0");
    TXTdensvinf->setText("-10.0");
    TXTdensxsup->setText("10.0");
    TXTdensysup->setText("10.0");
    TXTdenszsup->setText("10.0");
    TXTdensusup->setText("10.0");
    TXTdensvsup->setText("10.0");
    TXTdZJxformula2D->setText("u");
    TXTdZJyformula2D->setText("v");
    TXTdZJzformula2D->setText("0");
    RBTdZJrlow->setChecked(true);
    RBTdZJ3D->setChecked(true);
    RBTdZJrlow->setToolTip("65x65x65");
    RBTdZJrmedium->setToolTip("129x129x129");
    RBTdZJrhigh->setToolTip("257x257x257");
    SHTdZJxyz->clear();
    CHKdZJgrid->setChecked(true);
    WtabledZJ->setEnabled(false);
}

/*    READS OPTIONS FROM A PROJECT FILE (*.damproj)                */

void MainWindow::readOptions(const QString &fullFileName)
{
    QFile files(fullFileName);

    if (!files.isOpen()){
        if (!files.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("readOptions"),tr("File %1 cannot be read")
                .arg(fullFileName)+QString(":\n%1.").arg(files.errorString()));
            return;
        }
    }
    string file = toString(fullFileName);

    string v = CIniFile::GetValue("ImportFolder","PROJECTSECT",file);
    ImportFolder = QString(v.c_str());
    if (ImportFolder.isEmpty()) ImportFolder = Path(fullFileName);

    read_page_project(file, ImportFolder);      //  page_project: Project

    read_page_atdens(file);                     //    page_atdens: Atomic densities

    read_page_MED(file);                        //    page_MED: Density

    read_page_MESP(file);                       //    page_MESP: Electrostatic potential

    read_page_TOPO(file);                       //    page_TOPO: Topography

    read_page_SGhole(file);                      //   page_SGhole: MESP sigma hole
        
    read_page_HFforces(file);                   //    page_HFforces: Hellmann-Feynman forces

    read_page_Efield(file);                     //    page_Efield: Electric field

    read_page_densgrad(file);                   //    page_densgrad: Density gradient

    read_page_frad(file);                       //    page_frad: Radial factors

    read_page_orimult(file);                    //    page_orimult: Oriented multipoles

    read_page_MO(file);                         //    page_MO: Molecular orbitals

    read_page_ZJdens(file);                     //    page_ZJdens: Zernike-Jacobi expansions

    read_page_ZJtab(file);                      //    page_MO: Zernike-Jacobi density
            
//      End of options read

    page_atdens->setEnabled(true);
    page_MED->setEnabled(true);
    page_MESP->setEnabled(true);
    page_TOPO->setEnabled(true);
    page_HFforces->setEnabled(true);
    page_Efield->setEnabled(true);
    page_densgrad->setEnabled(true);
    page_frad->setEnabled(true);
    page_orimult->setEnabled(true);
    page_MO->setEnabled(true);
    page_ZJdens->setEnabled(true);
    page_ZJtab->setEnabled(true);

    files.close();
    SetCurrentFile("",false,false);
}

void MainWindow::read_page_project(string file, QString ImportFolder)
{
//    Read options of page_project: Project

    string v;
    QString qv;

    read_text_to_TXT("ProjectFolder", "PROJECTSECT", TXTProjectFolder, file);
    read_text_to_TXT("ProjectName", "PROJECTSECT", TXTProjectName, file);
    if (TXTProjectFolder->text().isEmpty())
        TXTProjectFolder->setText(ImportFolder);
    if (TXTProjectName->text().isEmpty())
        TXTProjectName->setText(FileWithoutExt(QString(file.c_str())));
    v = CIniFile::GetValue("ImportFile","PROJECTSECT",file);
    qv = toQString(v.c_str());
    if (!qv.isEmpty())
        TXTImport->setText(ImportFolder.append('/')+qv);
    else
        TXTImport->setText(ImportFolder.append('/')+TXTProjectName->text());

    TXTProjectFolder->setEnabled(true);
    TXTProjectName->setEnabled(true);
//    Default names are project name. Will be overwritten below if an alternative name has been given
    TXTdensdamdenfilename->setText(TXTProjectName->text());
    TXTHFgdamforcesfilename->setText(TXTProjectName->text());
    TXTeffilename->setText(TXTProjectName->text());
    TXTpotgdampotfilename->setText(TXTProjectName->text());
    TXTtopofilename->setText(TXTProjectName->text());
    TXTfraddamfilename->setText(TXTProjectName->text());
    TXTmrotorimultfilename->setText(TXTProjectName->text());
    TXTMOfilename->setText(TXTProjectName->text());

//    QString sgbsfile = ProjectFolder+"/"+FileWithoutExt(TXTProjectName->text())+".sgbs";
//    QString sgbsgzfile = ProjectFolder+"/"+FileWithoutExt(TXTProjectName->text())+".sgbs.gz";
//    QString sgbsdenfile = ProjectFolder+"/"+FileWithoutExt(TXTProjectName->text())+".sgbsden";
//    QString sgbsdengzfile = ProjectFolder+"/"+FileWithoutExt(TXTProjectName->text())+".sgbsden.gz";
//    QString ggbsfile = ProjectFolder+"/"+FileWithoutExt(TXTProjectName->text())+".ggbs";

    QString sgbsfile =      ProjectFolder+"/"+TXTProjectName->text()+".sgbs";
    QString sgbsgzfile =    ProjectFolder+"/"+TXTProjectName->text()+".sgbs.gz";
    QString sgbsdenfile =   ProjectFolder+"/"+TXTProjectName->text()+".sgbsden";
    QString sgbsdengzfile = ProjectFolder+"/"+TXTProjectName->text()+".sgbsden.gz";
    QString ggbsfile =      ProjectFolder+"/"+TXTProjectName->text()+".ggbs";

    if (QFile::exists(sgbsfile) || QFile::exists(sgbsgzfile)
            || QFile::exists(sgbsdenfile) || QFile::exists(sgbsdengzfile) ){
        lslater = true;
        CHKpotexact->setVisible(false);
    }
    else{
        lslater = false;
        CHKpotexact->setVisible(true);
    }
    CHKpotexact->setChecked(false);
    if (lslater){
        QString fileaux = FileWithoutExt(TXTImport->text());
        if (Extension(fileaux) == "sgbs" || Extension(fileaux) == "sgbsden" ) fileaux = FileWithoutExt(fileaux);
        QString sxyzfilename = ProjectFolder+"/"+ProjectName+".sxyz";
        if (!(QFile::exists(sxyzfilename))){
            execsgbs2sxyz(sxyzfilename);
        }
        set_natom(read_natom(sxyzfilename));
    }
    else{
        set_natom(read_natom(ggbsfile));
    }
}

void MainWindow::read_page_atdens(string file)
{
//    Read options of page_atdens: Atomic densities
    string v;
    QString qv;
    bool ok;
    int i;

    if (lslater){
        read_SPB("lmaxexp","DAMSECT",SPBatdenslmaxexp,file);
        read_SPB("lmultmx","DAMSECT",SPBatdenslmaxdisp,file);
        v = CIniFile::GetValue("umbral","DAMSECT",file);
        qv = toQString(v.c_str());
        qv.remove(0,3);
        i = qv.toInt(&ok);
        if (ok) SPBfradthreshold->setValue(qv.toInt());
        v = CIniFile::GetValue("umbralres","DAMSECT",file);
        qv = toQString(v.c_str());
        qv.remove(0,3);
        i = qv.toInt(&ok);
        if (ok) SPBfitthreshold->setValue(qv.toInt());
        v = CIniFile::GetValue("ioptaj","DAMSECT",file);
    }
    else{
        read_SPB("lmaxexp","G-DAMSECT",SPBatdenslmaxexp,file);
        read_SPB("lmultmx","G-DAMSECT",SPBatdenslmaxdisp,file);
        v = CIniFile::GetValue("umbral","G-DAMSECT",file);
        qv = toQString(v.c_str());
        qv.remove(0,3);
        i = qv.toInt(&ok);
        if (ok) SPBfradthreshold->setValue(qv.toInt());
        v = CIniFile::GetValue("umbralres","G-DAMSECT",file);
        qv = toQString(v.c_str());
        qv.remove(0,3);
        i = qv.toInt(&ok);
        if (ok) SPBfitthreshold->setValue(qv.toInt());
        v = CIniFile::GetValue("ioptaj","G-DAMSECT",file);
    }
    qv = toQString(v.c_str());
    int iv=qv.toInt();
    if (iv==1){
        RBTatdensDtotal->setChecked(true);
    }else if (iv==2){
        RBTatdensD1center->setChecked(true);
    }else if (iv==3){
        RBTatdensD2center->setChecked(true);
    }


}

void MainWindow::read_page_MED(string file)
{
//    Read options of page_MED: Density
    string v;
    QString qv;
    int i;
    bool ok;

    read_RBT("lexact","DAMDENSECT",RBTdensExact,file);
    RBTdensRep1->setChecked(!RBTdensExact->isChecked());

    read_RBT("ldeform","DAMDENSECT",RBTdensdeform,file);

    read_SPB("lmaxrep","DAMDENSECT",SPBdenslmaxexp,file);
    read_SPB("lminrep","DAMDENSECT",SPBdenslminexp,file);

    if (SPBdenslminexp->value() > 1){
        RBTdenslrange->setChecked(!RBTdensdeform->isChecked());
    }
    else if (SPBdenslminexp->value() == 1){
        RBTdensdeform->setChecked(true);
    }
    else{
        RBTdensfulldensity->setChecked(!RBTdensdeform->isChecked());
    }

    read_CHK("lgradient","DAMDENSECT",CHKdensgrad,file);
    read_CHK("lderiv2","DAMDENSECT",CHKdensder2,file);
    read_CHK("laplacian","DAMDENSECT",CHKdenslaplacian,file);
    read_CHK("lmolec","DAMDENSECT",CHKdenslmolec,file);
    read_CHK("latomics","DAMDENSECT",CHKdenslatomics,file);
    read_CHK("ldensacc","DAMDENSECT",CHKdensldensacc,file);

    v = CIniFile::GetValue("nsel","DAMDENSECT",file);
    qv = toQString(v.c_str());
    int nsel = qv.toInt();
    if (nsel > 0){
        denslist->clear();
        for(int i = 0 ; i < nsel ; ++i){
            qv=QString("%1%2%3").arg("iatomsel(").arg(i+1).arg(")");
            v=toString(qv);
            v = CIniFile::GetValue(v,"DAMDENSECT",file);
            qv = toQString(v.c_str());
            if (v!=""){
                denslist->append(qv);
            }
        }
        TXTdensatoms->setText(denslist->join(","));
    }

    read_text_to_TXT("filename", "DAMDENSECT", TXTdensdamdenfilename, file);

    read_double_to_TXT("xinf","DAMDENSECT",TXTdensxinf,file);
    read_double_to_TXT("xsup","DAMDENSECT",TXTdensxsup,file);
    read_double_to_TXT("yinf","DAMDENSECT",TXTdensyinf,file);
    read_double_to_TXT("ysup","DAMDENSECT",TXTdensysup,file);
    read_double_to_TXT("zinf","DAMDENSECT",TXTdenszinf,file);
    read_double_to_TXT("zsup","DAMDENSECT",TXTdenszsup,file);

    read_double_to_TXT("uinf","DAMDENSECT",TXTdensuinf,file);
    read_double_to_TXT("usup","DAMDENSECT",TXTdensusup,file);
    read_double_to_TXT("vinf","DAMDENSECT",TXTdensvinf,file);
    read_double_to_TXT("vsup","DAMDENSECT",TXTdensvsup,file);

    read_plot_dimension("lgrid2d","DAMDENSECT",file,RBTdens2D,RBTdensrlow,RBTdensrmedium,RBTdensrhigh);
    if (RBTdens3D->isChecked()) read_resolution("dltx","dlty","dltz","DAMDENSECT",file,TXTdensxinf,TXTdensxsup,TXTdensyinf,
        TXTdensysup,TXTdenszinf,TXTdenszsup,RBTdens3D,RBTdensrlow,RBTdensrmedium,RBTdensrhigh,RBTdensrcustom);
    else read_resolution("dltu","dltv","dltv","DAMDENSECT",file,TXTdensuinf,TXTdensusup,TXTdensvinf,TXTdensvsup,TXTdensvinf,TXTdensvsup,
        RBTdens3D,RBTdensrlow,RBTdensrmedium,RBTdensrhigh,RBTdensrcustom);

    read_text_to_TXT("x_func_uv", "DAMDENSECT", TXTdensxformula2D, file);
    read_text_to_TXT("y_func_uv", "DAMDENSECT", TXTdensyformula2D, file);
    read_text_to_TXT("z_func_uv", "DAMDENSECT", TXTdenszformula2D, file);

    v = CIniFile::GetValue("planecase","DAMDENSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) densplanecase = i;
    if (densplanecase > 0 && densplanecase < 8)
        RBTdensplane->setChecked(true);
    else
        RBTdensplane->setChecked(false);

    read_double_to_TXT("planeA","DAMDENSECT",TXTdensplaneA,file);
    read_double_to_TXT("planeB","DAMDENSECT",TXTdensplaneB,file);
    read_double_to_TXT("planeC","DAMDENSECT",TXTdensplaneC,file);

    if (TXTdensplaneC->text() == "0."){
        if (TXTdensplaneB->text() == "0."){
            RBTdensplaneYZ->setChecked(true);
        }
        else if(TXTdensplaneA->text() == "0."){
            RBTdensplaneXZ->setChecked(true);
        }
    }
    else if (TXTdensplaneA->text() == "0." && TXTdensplaneB->text() == "0."){
        RBTdensplaneXY->setChecked(true);
    }

    read_CHK("lgrid","DAMDENSECT",CHKdensgrid,file);
    read_CHK("lpoints","DAMDENSECT",CHKdensxyz,file);

    int numrtabulation = QString(CIniFile::GetValue("numrtab","DAMDENSECT",file).c_str()).toInt();
    SHTxyz->clear();
    QString valor;
    for (int i = 0 ; i < numrtabulation ; ++i){
        SHTxyz->resizeRows(numrtabulation+1);
        SHTxyz->resizeRows(SHTxyz->tabla->rowCount()+1);
        SHTxyz->tabla->insertRow(i);
        SHTxyz->resizeRows(SHTxyz->tabla->rowCount());
        for (int j = 0 ; j < 3 ; ++j){
            valor = toQString(CIniFile::GetValue(toString(QString("%1%2%3%4%5").arg("rtab(").arg(j+1).arg(",").arg(i+1).arg(")"))
                ,"DAMDENSECT",file).c_str()) ;
            SHTxyz->setcellvalue(valor,i,j);
        }
    }
    if (CHKdensxyz->isChecked()){
        Wtableden->setVisible(true);
        Wtableden->setEnabled(true);
    }
    else{
        Wtableden->setVisible(false);
        Wtableden->setEnabled(false);
    }
}

void MainWindow::read_page_MESP(string file)
{
    activebeware = false;
//    Read options of page_MESP: Electrostatic potential
    string v;
    QString qv;
    int i;
    bool ok;

    read_text_to_TXT("filename", "DAMPOTSECT", TXTpotgdampotfilename, file);

    read_SPB("lmaxrep","DAMPOTSECT",SPBpotlmaxexp,file);
    read_CHK("largo","DAMPOTSECT",CHKpotlong,file);

    v = CIniFile::GetValue("umbrlargo","DAMPOTSECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBpotlongthreshold->setValue(qv.toInt());

    read_CHK("lexact","DAMPOTSECT",CHKpotexact,file);
    read_CHK("lgradient","DAMPOTSECT",CHKpotgrad,file);
    read_CHK("lderiv2","DAMPOTSECT",CHKpotder2,file);

    read_double_to_TXT("xinf","DAMPOTSECT",TXTpotxinf,file);
    read_double_to_TXT("xsup","DAMPOTSECT",TXTpotxsup,file);
    read_double_to_TXT("yinf","DAMPOTSECT",TXTpotyinf,file);
    read_double_to_TXT("ysup","DAMPOTSECT",TXTpotysup,file);
    read_double_to_TXT("zinf","DAMPOTSECT",TXTpotzinf,file);
    read_double_to_TXT("zsup","DAMPOTSECT",TXTpotzsup,file);

    read_double_to_TXT("uinf","DAMPOTSECT",TXTpotuinf,file);
    read_double_to_TXT("usup","DAMPOTSECT",TXTpotusup,file);
    read_double_to_TXT("vinf","DAMPOTSECT",TXTpotvinf,file);
    read_double_to_TXT("vsup","DAMPOTSECT",TXTpotvsup,file);

    read_plot_dimension("lgrid2d","DAMPOTSECT",file,RBTpot2D,RBTpotrlow,RBTpotrmedium,RBTpotrhigh);
    if (RBTpot3D->isChecked()) read_resolution("dltx","dlty","dltz","DAMPOTSECT",file,TXTpotxinf,TXTpotxsup,TXTpotyinf,TXTpotysup,TXTpotzinf,TXTpotzsup,
        RBTpot3D,RBTpotrlow,RBTpotrmedium,RBTpotrhigh,RBTpotrcustom);
    else read_resolution("dltu","dltv","dltv","DAMPOTSECT",file,TXTpotuinf,TXTpotusup,TXTpotvinf,TXTpotvsup,TXTpotvinf,TXTpotvsup,
        RBTpot3D,RBTpotrlow,RBTpotrmedium,RBTpotrhigh,RBTpotrcustom);

    read_text_to_TXT("x_func_uv", "DAMPOTSECT", TXTpotxformula2D, file);
    read_text_to_TXT("y_func_uv", "DAMPOTSECT", TXTpotyformula2D, file);
    read_text_to_TXT("z_func_uv", "DAMPOTSECT", TXTpotzformula2D, file);

    v = CIniFile::GetValue("planecase","DAMPOTSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) potplanecase = i;
    if (potplanecase > 0 && potplanecase < 8)
        RBTpotplane->setChecked(true);
    else
        RBTpotplane->setChecked(false);

    read_double_to_TXT("planeA","DAMPOTSECT",TXTpotplaneA,file);
    read_double_to_TXT("planeB","DAMPOTSECT",TXTpotplaneB,file);
    read_double_to_TXT("planeC","DAMPOTSECT",TXTpotplaneC,file);

    if (TXTpotplaneC->text() == "0."){
        if (TXTpotplaneB->text() == "0."){
            RBTpotplaneYZ->setChecked(true);
        }
        else if(TXTpotplaneA->text() == "0."){
            RBTpotplaneXZ->setChecked(true);
        }
    }
    else if (TXTpotplaneA->text() == "0." && TXTpotplaneB->text() == "0."){
        RBTpotplaneXY->setChecked(true);
    }

    read_CHK("lgrid","DAMPOTSECT",CHKpotgrid,file);
    read_CHK("lpoints","DAMPOTSECT",CHKpotxyz,file);

    int numrtabpot = QString(CIniFile::GetValue("numrtab","DAMPOTSECT",file).c_str()).toInt();
    SHTpotxyz->clear();
    QString valorpot;
    for (int i = 0 ; i < numrtabpot ; ++i){
        SHTpotxyz->resizeRows(numrtabpot+1);
        SHTpotxyz->resizeRows(SHTpotxyz->tabla->rowCount()+1);
        SHTpotxyz->tabla->insertRow(i);
        SHTpotxyz->resizeRows(SHTpotxyz->tabla->rowCount());
        for (int j = 0 ; j < 3 ; ++j){
            valorpot = toQString(CIniFile::GetValue(toString(QString("%1%2%3%4%5").arg("rtab(").arg(j+1).arg(",").arg(i+1).arg(")"))
                ,"DAMPOTSECT",file).c_str()) ;
            SHTpotxyz->setcellvalue(valorpot,i,j);
        }
    }

    if (CHKpotxyz->isChecked()){
        Wtablepot->setVisible(true);
        Wtablepot->setEnabled(true);
    }
    else{
        Wtablepot->setVisible(false);
        Wtablepot->setEnabled(false);
    }
    activebeware = true;
}

void MainWindow::read_page_TOPO(string file)
{
//    Read options of page_TOPO: Topography

    read_text_to_TXT("filename", "DAMTOPOSECT", TXTtopofilename, file);
    read_CHK("TOPOGRAPH","DAMTOPOSECT",CHKtopograph,file);
    CHKtopograph_changed();

    read_CHK("addguess","DAMTOPOSECT",CHKtopoaddguess,file);
    read_CHK("addguess","DAMTOPOSECT",CHKtopoxyz,file);

    int numtopoguess = QString(CIniFile::GetValue("ncntguess","DAMTOPOSECT",file).c_str()).toInt();
    SHTtopoxyz->clear();
    QString valorcnt;
    for (int i = 0 ; i < numtopoguess ; ++i){
        SHTtopoxyz->resizeRows(numtopoguess+1);
        SHTtopoxyz->resizeRows(SHTtopoxyz->tabla->rowCount()+1);
        SHTtopoxyz->tabla->insertRow(i);
        SHTtopoxyz->resizeRows(SHTtopoxyz->tabla->rowCount());
        for (int j = 0 ; j < 3 ; ++j){
            valorcnt = toQString(CIniFile::GetValue(toString(QString("%1%2%3%4%5").arg("rcntguess(").arg(j+1).arg(",").arg(i+1).arg(")"))
                ,"DAMTOPOSECT",file).c_str()) ;
            SHTtopoxyz->setcellvalue(valorcnt,i,j);
        }
    }

    if (CHKtopoaddguess->isChecked() && CHKtopoxyz->isChecked()){
        Wtabledentopo->setVisible(true);
        Wtabledentopo->setEnabled(true);
    }
    else{
        Wtabledentopo->setVisible(false);
        Wtabledentopo->setEnabled(false);
    }

    read_text_to_TXT("guessfile","DAMTOPOSECT",TXTtopoguessfilename,file);

    read_double_to_TXT("BOXL","DAMTOPOSECT",TXTtopoboxl,file);
    read_double_to_TXT("BOXT","DAMTOPOSECT",TXTtopoboxt,file);
    read_double_to_TXT("stepszt","DAMTOPOSECT",TXTtopostepszt,file);

    read_RBT("MED","DAMTOPOSECT",RBTtopodensity,file);
    read_RBT("MESP","DAMTOPOSECT",RBTtopopotential,file);

    read_double_to_TXT("CNVG","DAMTOPOSECT",TXTtopocnvg,file);
    read_SPB("lmaxi","DAMFIELDSECT",SPBtopolmaxi,file);

    read_CHK("graphpath","DAMTOPOSECT",CHKtopomolgraph,file);
    read_double_to_TXT("BOXG","DAMTOPOSECT",TXTtopoboxg,file);
    read_double_to_TXT("DRCUTCP","DAMTOPOSECT",TXTtopoggradthr,file);

    read_CHK("basin","DAMTOPOSECT",CHKtopobasin,file);
    read_CHK("exdraw","DAMTOPOSECT",CHKtopoexdraw,file);
    read_double_to_TXT("exln","DAMTOPOSECT",TXTtopoexln,file);
    read_double_to_TXT("BOXB","DAMTOPOSECT",TXTtopoboxb,file);
    read_double_to_TXT("FDISP","DAMTOPOSECT",TXTtopofdisp,file);
}

void MainWindow::read_page_SGhole(string file)
{
    activebeware = false;
//    Read options of page_SGhole: MESP sigma hole
    string v;
    QString qv;
    int i;
    bool ok;
    read_text_to_TXT("filename", "DAMSGHOLESECT", TXTSGholefilename, file);
    read_text_to_TXT("gridname", "DAMSGHOLESECT", TXTImportSGholeden, file);
    read_double_to_TXT("contourval","DAMSGHOLESECT",TXTSGholecontour,file);
//    read_double_to_TXT("topcolor","DAMSGHOLESECT",TXTSGholecolor,file);

    v = CIniFile::GetValue("geomthr","DAMSGHOLESECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBSGholegeomthreshold->setValue(qv.toInt());

    v = CIniFile::GetValue("thrslocal","DAMSGHOLESECT",file);
    qv = toQString(v.c_str());
    qv.remove(2,4);
    i = qv.toInt(&ok);
    if (ok) SPBSGholelocalextrema->setValue(qv.toInt());

    v = CIniFile::GetValue("umbrlargo","DAMSGHOLESECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBSGholelongthreshold->setValue(qv.toInt());

    read_SPB("lmaxrep","DAMSGHOLESECT",SPBSGholelmaxexp,file);
    read_CHK("lexact","DAMSGHOLESECT",CHKSGhexactMESP,file);
    activebeware = true;
}

void MainWindow::read_page_HFforces(string file)
{
//    Read options of page_HFforces: Hellmann-Feynman forces
    string v;
    QString qv;

    read_CHK("latomsel","DAMFORCESSECT",CHKHFlatomsel,file);

    v = CIniFile::GetValue("ncntab","DAMFORCESSECT",file);
    qv = toQString(v.c_str());
    int nforces = qv.toInt();
    if (nforces > 0){
        HFforceslist->clear();
        for(int i = 0 ; i < nforces ; ++i){
            qv=QString("%1%2%3").arg("iatomsel(").arg(i+1).arg(")");
            v=toString(qv);
            v = CIniFile::GetValue(v,"DAMFORCESSECT",file);
            qv = toQString(v.c_str());
            if (v!=""){
                HFforceslist->append(qv);
            }
        }
        TXTHFforcesatoms->setText(HFforceslist->join(","));
    }

    if (CHKHFlatomsel->isChecked()){
        TXTHFforcesatoms->setEnabled(true);
    }
    else{
        TXTHFforcesatoms->setEnabled(false);
    }

    read_text_to_TXT("filename", "DAMFORCESSECT", TXTHFgdamforcesfilename, file);
}

void MainWindow::read_page_Efield(string file)
{
//    Read options of page_Efield: Electric field
    string v;
    QString qv;
    int i;
    bool ok;

    read_text_to_TXT("filename", "DAMFIELDSECT", TXTeffilename, file);

    read_SPB("lmaxrep","DAMFIELDSECT",SPBeflmaxexp,file);

    read_CHK("largo","DAMFIELDSECT",CHKeflong,file);

    v = CIniFile::GetValue("umbrlargo","DAMFIELDSECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBeflongthreshold->setValue(qv.toInt());

    v = CIniFile::GetValue("numpnt","DAMFIELDSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) TXTefnumpnt->setText(qv);
    v = CIniFile::GetValue("nlinpernuc","DAMFIELDSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) TXTefnlinpernuc->setText(qv);
    v = CIniFile::GetValue("ioplines3D","DAMFIELDSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) CMBefdirectionset->setCurrentIndex(i);

    v = CIniFile::GetValue("dlt0","DAMFIELDSECT",file);
    qv = toQString(v.c_str());
    qv.toDouble(&ok);
    if (qv=="" || !ok) qv="0.02";
    TXTefdlt0->setText(qv);

    read_double_to_TXT("xinf","DAMFIELDSECT",TXTefxinf,file);
    read_double_to_TXT("yinf","DAMFIELDSECT",TXTefyinf,file);
    read_double_to_TXT("zinf","DAMFIELDSECT",TXTefzinf,file);
    read_double_to_TXT("xsup","DAMFIELDSECT",TXTefxsup,file);
    read_double_to_TXT("ysup","DAMFIELDSECT",TXTefysup,file);
    read_double_to_TXT("zsup","DAMFIELDSECT",TXTefzsup,file);

    read_double_to_TXT("uinf","DAMFIELDSECT",TXTefuinf,file);
    read_double_to_TXT("vinf","DAMFIELDSECT",TXTefvinf,file);
    read_double_to_TXT("usup","DAMFIELDSECT",TXTefusup,file);
    read_double_to_TXT("vsup","DAMFIELDSECT",TXTefvsup,file);

    read_double_to_TXT("planeA","DAMFIELDSECT",TXTefplaneA,file);
    read_double_to_TXT("planeB","DAMFIELDSECT",TXTefplaneB,file);
    read_double_to_TXT("planeC","DAMFIELDSECT",TXTefplaneC,file);
    read_double_to_TXT("uvratio","DAMFIELDSECT",TXTefuvratio,file);

    if (TXTefplaneC->text() == "0."){
        if (TXTefplaneB->text() == "0."){
            RBTefplaneYZ->setChecked(true);
        }
        else if(TXTefplaneA->text() == "0."){
            RBTefplaneXZ->setChecked(true);
        }
    }
    else if (TXTefplaneA->text() == "0." && TXTefplaneB->text() == "0."){
        RBTefplaneXY->setChecked(true);
    }

    read_CHK("lextralines","DAMFIELDSECT",CHKefextralines,file);

    read_text_to_TXT("filelines", "DAMFIELDSECT", TXTeffilelines, file);

    int nlines = QString(CIniFile::GetValue("nlines","DAMFIELDSECT",file).c_str()).toInt();
    SHTefxyz->clear();
    QString valorpot;
    for (int i = 0 ; i < nlines ; ++i){
        SHTefxyz->resizeRows(nlines+1);
        SHTefxyz->resizeRows(SHTefxyz->tabla->rowCount()+1);
        SHTefxyz->tabla->insertRow(i);
        SHTefxyz->resizeRows(SHTefxyz->tabla->rowCount());
        valorpot = toQString(CIniFile::GetValue(toString(QString("%1%2%3").arg("icntlines(").arg(i+1).arg(")"))
            ,"DAMFIELDSECT",file).c_str()) ;
        SHTefxyz->setcellvalue(valorpot,i,0);
        for (int j = 1 ; j < 4 ; ++j){
            valorpot = toQString(CIniFile::GetValue(toString(QString("%1%2%3%4%5").arg("rlines(").arg(j).arg(",").arg(i+1).arg(")"))
                ,"DAMFIELDSECT",file).c_str()) ;
            SHTefxyz->setcellvalue(valorpot,i,j);
        }
    }

    if (CHKefextralines->isChecked() && nlines > 0){
        Wtabledenef->setVisible(true);
        Wtabledenef->setEnabled(true);
        CHKefxyz->setChecked(true);
    }
    else{
        Wtabledenef->setVisible(false);
        Wtabledenef->setEnabled(false);
        CHKefxyz->setChecked(false);
    }
}

void MainWindow::read_page_densgrad(string file)
{
//    Read options of page_densgrad: Density gradient
    string v;
    QString qv;
    int i;
    bool ok;

    read_text_to_TXT("filename", "DAMDENGRADSECT", TXTdgfilename, file);

    read_SPB("lmaxrep","DAMDENGRADSECT",SPBdglmaxexp,file);

    v = CIniFile::GetValue("umbrlargo","DAMDENGRADSECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBdglongthreshold->setValue(qv.toInt());

    read_CHK("lextralines","DAMDENGRADSECT",CHKdgextralines,file);

    v = CIniFile::GetValue("numpnt","DAMDENGRADSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) TXTdgnumpnt->setText(qv);
    v = CIniFile::GetValue("nlinpernuc","DAMDENGRADSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) TXTdgnlinpernuc->setText(qv);
    v = CIniFile::GetValue("ioplines3D","DAMDENGRADSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) CMBdgdirectionset->setCurrentIndex(i);

    v = CIniFile::GetValue("dlt0","DAMDENGRADSECT",file);
    qv = toQString(v.c_str());
    qv.toDouble(&ok);
    if (qv=="" || !ok) qv="0.02";
    TXTdgdlt0->setText(qv);

    read_double_to_TXT("xinf","DAMDENGRADSECT",TXTdgxinf,file);
    read_double_to_TXT("yinf","DAMDENGRADSECT",TXTdgyinf,file);
    read_double_to_TXT("zinf","DAMDENGRADSECT",TXTdgzinf,file);
    read_double_to_TXT("xsup","DAMDENGRADSECT",TXTdgxsup,file);
    read_double_to_TXT("ysup","DAMDENGRADSECT",TXTdgysup,file);
    read_double_to_TXT("zsup","DAMDENGRADSECT",TXTdgzsup,file);

    read_double_to_TXT("uinf","DAMDENGRADSECT",TXTdguinf,file);
    read_double_to_TXT("vinf","DAMDENGRADSECT",TXTdgvinf,file);
    read_double_to_TXT("usup","DAMDENGRADSECT",TXTdgusup,file);
    read_double_to_TXT("vsup","DAMDENGRADSECT",TXTdgvsup,file);

    read_double_to_TXT("planeA","DAMDENGRADSECT",TXTdgplaneA,file);
    read_double_to_TXT("planeB","DAMDENGRADSECT",TXTdgplaneB,file);
    read_double_to_TXT("planeC","DAMDENGRADSECT",TXTdgplaneC,file);
    read_double_to_TXT("uvratio","DAMDENGRADSECT",TXTdguvratio,file);

    if (TXTdgplaneC->text() == "0."){
        if (TXTdgplaneB->text() == "0."){
            RBTdgplaneYZ->setChecked(true);
        }
        else if(TXTdgplaneA->text() == "0."){
            RBTdgplaneXZ->setChecked(true);
        }
    }
    else if (TXTdgplaneA->text() == "0." && TXTdgplaneB->text() == "0."){
        RBTdgplaneXY->setChecked(true);
    }

    read_CHK("lextralines","DAMDENGRADSECT",CHKdgextralines,file);

    read_text_to_TXT("filelines", "DAMDENGRADSECT", TXTdgfilelines, file);

    int nlines = QString(CIniFile::GetValue("nlines","DAMDENGRADSECT",file).c_str()).toInt();
    SHTdgxyz->clear();
    QString valorpot;
    for (int i = 0 ; i < nlines ; ++i){
        SHTdgxyz->resizeRows(nlines+1);
        SHTdgxyz->resizeRows(SHTdgxyz->tabla->rowCount()+1);
        SHTdgxyz->tabla->insertRow(i);
        SHTdgxyz->resizeRows(SHTdgxyz->tabla->rowCount());
        valorpot = toQString(CIniFile::GetValue(toString(QString("%1%2%3").arg("icntlines(").arg(i+1).arg(")"))
            ,"DAMDENGRADSECT",file).c_str()) ;
        SHTdgxyz->setcellvalue(valorpot,i,0);
        for (int j = 1 ; j < 4 ; ++j){
            valorpot = toQString(CIniFile::GetValue(toString(QString("%1%2%3%4%5").arg("rlines(").arg(j).arg(",").arg(i+1).arg(")"))
                ,"DAMDENGRADSECT",file).c_str()) ;
            SHTdgxyz->setcellvalue(valorpot,i,j);
        }
    }

    if (CHKdgextralines->isChecked() && nlines > 0){
        Wtabledendg->setVisible(true);
        Wtabledendg->setEnabled(true);
        CHKdgxyz->setChecked(true);
    }
    else{
        Wtabledendg->setVisible(false);
        Wtabledendg->setEnabled(false);
        CHKdgxyz->setChecked(false);
    }
}

void MainWindow::read_page_frad(string file)
{
//    Read options of page_frad: Radial factors
    string v;
    QString qv;
    int i;
    bool ok;
    int row;
    int col;

    read_text_to_TXT("filename", "DAMFRADSECT", TXTfraddamfilename, file);

    read_double_to_TXT("rini","DAMFRADSECT",TXTfradrini,file);
    read_double_to_TXT("rfin","DAMFRADSECT",TXTfradrfin,file);
    read_double_to_TXT("dltr","DAMFRADSECT",TXTfradltr,file);

    v = CIniFile::GetValue("iatom","DAMFRADSECT",file);
    v = CIniFile::GetValue("ltab","DAMFRADSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBfradltab->setValue(qv.toInt());
    v = CIniFile::GetValue("mtab","DAMFRADSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBfradmtab->setValue(qv.toInt());

    read_CHK("lrlist","DAMFRADSECT",CHKfradextras,file);
    if (CHKfradextras->isChecked()){
        SHTfradrlist->clear();
        row=-1;
        col=0;
        for (int i=1 ; i<=Sheet::max_sel ; ++i){
            qv=QString("%1%2%3").arg("rlist(").arg(i).arg(")");
            v=toString(qv);
            v = CIniFile::GetValue(v,"DAMFRADSECT",file);
            if (v!=""){
                qv = toQString(v.c_str());
                row++;
                SHTfradrlist->resizeRows(SHTfradrlist->tabla->rowCount()+1);
                SHTfradrlist->tabla->insertRow(row);
                SHTfradrlist->resizeRows(SHTfradrlist->tabla->rowCount());
                SHTfradrlist->setcellvalue(qv,row,col);
            }
        }
    }

    v = CIniFile::GetValue("ncntab","DAMFRADSECT",file);
    qv = toQString(v.c_str());
    int nfrad = qv.toInt();
    if (nfrad > 0){
        fradlist->clear();
        for(int i = 0 ; i < nfrad ; ++i){
            qv=QString("%1%2%3").arg("iatomsel(").arg(i+1).arg(")");
            v=toString(qv);
            v = CIniFile::GetValue(v,"DAMFRADSECT",file);
            qv = toQString(v.c_str());
            if (v!=""){
                fradlist->append(qv);
            }
        }
        TXTfradatoms->setText(fradlist->join(","));
    }

    read_CHK("lderiv","DAMFRADSECT",CHKfradderiv1,file);
    read_CHK("lderiv2","DAMFRADSECT",CHKfradderiv2,file);
}

void MainWindow::read_page_orimult(string file)
{
//    Read options of page_orimult: Oriented multipoles
    string v;
    QString qv;
    int i;
    bool ok;

    SPBmrotlmin->setRange(0,SPBatdenslmaxexp->value());

    v = CIniFile::GetValue("lmax","DAMMULTROTSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBmrotlmax->setValue(qv.toInt());

    v = CIniFile::GetValue("lmin","DAMMULTROTSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBmrotlmin->setValue(i);

    SPBmrotleft->setRange(1,get_natom());
    SPBmrotmiddle->setRange(1,get_natom());
    SPBmrotright->setRange(1,get_natom());

    v = CIniFile::GetValue("i1","DAMMULTROTSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBmrotleft->setValue(i);

    v = CIniFile::GetValue("i2","DAMMULTROTSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBmrotmiddle->setValue(i);

    v = CIniFile::GetValue("i3","DAMMULTROTSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) SPBmrotright->setValue(i);

    read_text_to_TXT("filename", "DAMMULTROTSECT", TXTmrotorimultfilename, file);
}

void MainWindow::read_page_MO(string file)
{
//    Read options of page_MO: Molecular orbitals
    string v;
    QString qv;
    int i;
    bool ok;

    read_text_to_TXT("filename", "DAMORBSECT", TXTMOfilename, file);
    read_text_to_TXT("fileMOname", "DAMORBSECT", TXTMOImportfile, file);

    read_CHK("lgradient","DAMORBSECT",CHKMOgrad,file);
    read_double_to_TXT("xinf","DAMORBSECT",TXTMOxinf,file);
    read_double_to_TXT("xsup","DAMORBSECT",TXTMOxsup,file);

    read_double_to_TXT("yinf","DAMORBSECT",TXTMOyinf,file);
    read_double_to_TXT("ysup","DAMORBSECT",TXTMOysup,file);
    read_double_to_TXT("zinf","DAMORBSECT",TXTMOzinf,file);
    read_double_to_TXT("zsup","DAMORBSECT",TXTMOzsup,file);

    read_double_to_TXT("uinf","DAMORBSECT",TXTMOuinf,file);
    read_double_to_TXT("usup","DAMORBSECT",TXTMOusup,file);
    read_double_to_TXT("vinf","DAMORBSECT",TXTMOvinf,file);
    read_double_to_TXT("vsup","DAMORBSECT",TXTMOvsup,file);

    read_plot_dimension("lgrid2d","DAMORBSECT",file,RBTMO2D,RBTMOrlow,RBTMOrmedium,RBTMOrhigh);
    if (RBTMO3D->isChecked()) read_resolution("dltx","dlty","dltz","DAMORBSECT",file,TXTMOxinf,TXTMOxsup,TXTMOyinf,TXTMOysup,
        TXTMOzinf,TXTMOzsup,RBTMO3D,RBTMOrlow,RBTMOrmedium,RBTMOrhigh,RBTMOrcustom);
    else read_resolution("dltu","dltv","dltv","DAMORBSECT",file,TXTMOuinf,TXTMOusup,TXTMOvinf,TXTMOvsup,TXTMOvinf,TXTMOvsup,
        RBTMO3D,RBTMOrlow,RBTMOrmedium,RBTMOrhigh,RBTMOrcustom);

    read_text_to_TXT("x_func_uv", "DAMORBSECT", TXTMOxformula2D, file);
    read_text_to_TXT("y_func_uv", "DAMORBSECT", TXTMOyformula2D, file);
    read_text_to_TXT("z_func_uv", "DAMORBSECT", TXTMOzformula2D, file);

    v = CIniFile::GetValue("planecase","DAMORBSECT",file);
    qv = toQString(v.c_str());
    i = qv.toInt(&ok);
    if (ok) MOplanecase = i;
    if (MOplanecase > 0 && MOplanecase < 8)
        RBTMOplane->setChecked(true);
    else
        RBTMOplane->setChecked(false);

    read_double_to_TXT("planeA","DAMORBSECT",TXTMOplaneA,file);
    read_double_to_TXT("planeB","DAMORBSECT",TXTMOplaneB,file);
    read_double_to_TXT("planeC","DAMORBSECT",TXTMOplaneC,file);

    if (TXTMOplaneC->text() == "0."){
        if (TXTMOplaneB->text() == "0."){
            RBTMOplaneYZ->setChecked(true);
        }
        else if(TXTMOplaneA->text() == "0."){
            RBTMOplaneXZ->setChecked(true);
        }
    }
    else if (TXTMOplaneA->text() == "0." && TXTMOplaneB->text() == "0."){
        RBTMOplaneXY->setChecked(true);
    }

    v = CIniFile::GetValue("norbs","DAMORBSECT",file);
    qv = toQString(v.c_str());
    int norb = qv.toInt();
    if (norb > 0){
        MOlist->clear();
        for(int i = 0 ; i < norb ; ++i){
            qv=QString("%1%2%3").arg("iorbsinp(").arg(i+1).arg(")");
            v=toString(qv);
            v = CIniFile::GetValue(v,"DAMORBSECT",file);
            qv = toQString(v.c_str());
            if (v!=""){
                MOlist->append(qv);
            }
        }
        TXTMOchoose->setText(MOlist->join(","));
    }
}

void MainWindow::read_page_ZJdens(string file)
{
//    Read options of page_ZJdens: One-center Zernike 3D-jacobi expansion
    string v;
    QString qv;
    int i;
    bool ok;

    read_SPB("lexpansion","DAMZJSECT",SPBZJlmax,file);
    read_SPB("kexpansion","DAMZJSECT",SPBZJkmax,file);
    read_RBT("ljacobi","DAMZJSECT",RBTZJacobi,file);
    if (RBTZJacobi->isChecked()) RBTZJZernike->setChecked(false);
    else RBTZJZernike->setChecked(true);

    read_RBT("lstarrel","DAMZJSECT",RBTZJrstarrel,file);
    read_double_to_TXT("rstar","DAMZJSECT",TXTZJrstar,file);
    RBTZJrstarabs->setChecked(true);

    read_RBT("lechelon","DAMZJSECT",RBTZJechelon,file);
    read_SPB("nquadpoints","DAMZJSECT",SPBZJnquad,file);

    v = CIniFile::GetValue("thresmult","DAMZJSECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBZJthrmult->setValue(qv.toInt());

    v = CIniFile::GetValue("thresoverlap","DAMZJSECT",file);
    qv = toQString(v.c_str());
    qv.remove(0,3);
    i = qv.toInt(&ok);
    if (ok) SPBZJthrdist->setValue(qv.toInt());
}

void MainWindow::read_page_ZJtab(string file)
{
//    Read options of page_ZJtab: Tabulation of density from Zernike 3D-Jacobi expansion
    string v;
    QString qv;
    int numrtabulation;
    QString valor;

    read_SPB("kmaxrep","DAMDENZJSECT",SPBdZJkmax,file);
    read_SPB("lmaxrep","DAMDENZJSECT",SPBdZJlmax,file);
    read_SPB("lminrep","DAMDENZJSECT",SPBdZJlmin,file);

    ZJlist->clear();
    TXTdZJchoose->setText("");

    read_CHK("lgradient","DAMDENZJSECT",CHKdZJgrad,file);

    read_text_to_TXT("filename", "DAMDENZJSECT", TXTdZJfilename, file);
    read_text_to_TXT("fileZJname", "DAMDENZJSECT", TXTdZJImportfile, file);

    read_double_to_TXT("xinf","DAMDENZJSECT",TXTdZJxinf,file);
    read_double_to_TXT("xsup","DAMDENZJSECT",TXTdZJxsup,file);
    read_double_to_TXT("yinf","DAMDENZJSECT",TXTdZJyinf,file);
    read_double_to_TXT("ysup","DAMDENZJSECT",TXTdZJysup,file);
    read_double_to_TXT("zinf","DAMDENZJSECT",TXTdZJzinf,file);
    read_double_to_TXT("zsup","DAMDENZJSECT",TXTdZJzsup,file);

    read_double_to_TXT("uinf","DAMDENZJSECT",TXTdZJuinf,file);
    read_double_to_TXT("usup","DAMDENZJSECT",TXTdZJusup,file);
    read_double_to_TXT("vinf","DAMDENZJSECT",TXTdZJvinf,file);
    read_double_to_TXT("vsup","DAMDENZJSECT",TXTdZJvsup,file);

    read_plot_dimension("lgrid2d","DAMDENZJSECT",file,RBTdZJ2D,RBTdZJrlow,RBTdZJrmedium,RBTdZJrhigh);
    if (RBTdZJ3D->isChecked()) read_resolution("dltx","dlty","dltz","DAMDENZJSECT",file,TXTdZJxinf,TXTdZJxsup,TXTdZJyinf,
        TXTdZJysup,TXTdZJzinf,TXTdZJzsup,RBTdZJ3D,RBTdZJrlow,RBTdZJrmedium,RBTdZJrhigh,RBTdZJrcustom);
    else read_resolution("dltu","dltv","dltv","DAMDENZJSECT",file,TXTdZJuinf,TXTdZJusup,TXTdZJvinf,TXTdZJvsup,TXTdZJvinf,TXTdZJvsup,
        RBTdZJ3D,RBTdZJrlow,RBTdZJrmedium,RBTdZJrhigh,RBTdZJrcustom);

    read_text_to_TXT("x_func_uv", "DAMDENZJSECT", TXTdZJxformula2D, file);
    read_text_to_TXT("y_func_uv", "DAMDENZJSECT", TXTdZJyformula2D, file);
    read_text_to_TXT("z_func_uv", "DAMDENZJSECT", TXTdZJzformula2D, file);

    numrtabulation = QString(CIniFile::GetValue("numrtab","DAMDENZJSECT",file).c_str()).toInt();
    SHTdZJxyz->clear();
    for (int i = 0 ; i < numrtabulation ; ++i){
        SHTdZJxyz->resizeRows(numrtabulation+1);
        SHTdZJxyz->resizeRows(SHTdZJxyz->tabla->rowCount()+1);
        SHTdZJxyz->tabla->insertRow(i);
        SHTdZJxyz->resizeRows(SHTdZJxyz->tabla->rowCount());
        for (int j = 0 ; j < 3 ; ++j){
            valor = toQString(CIniFile::GetValue(toString(QString("%1%2%3%4%5").arg("rtab(").arg(j+1).arg(",").arg(i+1).arg(")"))
                ,"DAMDENZJSECT",file).c_str()) ;
            SHTdZJxyz->setcellvalue(valor,i,j);
        }
    }

    read_CHK("lgrid","DAMDENZJSECT",CHKdZJgrid,file);
    read_CHK("lpoints","DAMDENZJSECT",CHKdZJxyz,file);
}

// Reads the content of a section in options .damproj file
QByteArray MainWindow::ReadSectionOptions(const char *SectionName, QFile *FileName)
{
    QByteArray buff, line;
    bool lreadend = false;    // True when reading section ended
    buff.append("&OPTIONS\n");
    while(!(*FileName).atEnd() && !lreadend){
        line = (*FileName).readLine(50);
        if( line.contains(SectionName) ) {
            while(!(*FileName).atEnd()){
                line = (*FileName).readLine(50);
                if( line.contains("[") ) {
                    lreadend = true;
                    break;
                }
                buff.append(line);
            }
        }
    }
    buff.append("&END\n");
    return buff;
}

//    Saves options in a project file (*.damproj)
void MainWindow::saveOptions(const QString &fullFileName,int clase)
{
    QString filezdo = ProjectFolder + "zdo";
    if (QFileInfo(filezdo).exists()){
        lzdo = true;
    }
    else{
        lzdo = false;
    }
    QString filevalence = ProjectFolder + "valence";
    if (QFileInfo(filevalence).exists()){
        lvalence = true;
    }
    else{
        lvalence = false;
    }
    string file = toString(fullFileName);
    bool *printwarns = new bool;

    QString *warns = new QString(tr("Warning: failed saving the following options") + ":\n");
    QFile files(fullFileName);
    if (!files.isOpen()){
        files.open(QFile::Append | QFile::WriteOnly);
    }
    *printwarns = false;

    if (clase==0){        // Project
        saveOptions_0(file, printwarns, warns);
    }

    if (clase==0 || clase==1){    // Project || G-DAM || DAM
        saveOptions_1(file, printwarns, warns);
    }

    if (clase==0 || clase==2){    // Project || DAMDEN
        saveOptions_2(file, printwarns, warns);
    }

    if (clase==0 || clase==3){    // Project || DAMPOT
        saveOptions_3(file, printwarns, warns);
    }

    if (clase==0 || clase==4){    // Project || DAMFORCES
        saveOptions_4(file, printwarns, warns);
    }

    if (clase==0 || clase==5){    // Project || DAMFIELD
        saveOptions_5(file, printwarns, warns);
    }

    if (clase==0 || clase==6){    // Project || DAMFRAD
        saveOptions_6(file, printwarns, warns);
    }

    if (clase==0 || clase==7){    // Project || DAMMULTROT
        saveOptions_7(file, printwarns, warns);
    }

    if (clase==0 || clase==8){    // Project || DAMORB
        saveOptions_8(file, printwarns, warns);
    }

    if (clase==0 || clase==9){    // Project || DAMTOPOGRAPHY
        saveOptions_9(file, printwarns, warns);
    }

    if (clase==0 || clase==10){    // Project || Zernike-Jacobi
        saveOptions_10(file, printwarns, warns);
    }

    if (clase==0 || clase==11){    // Project || DAMDENZJ
        saveOptions_11(file, printwarns, warns);
    }

    if (clase==0 || clase==12){    // Project || DAMDENGRAD
        saveOptions_12(file, printwarns, warns);
    }

    if (clase==0 || clase==13){    // Project || DAMSGHOLE
        saveOptions_13(file, printwarns, warns);
    }
    
    if (*printwarns){
        QMessageBox::warning(this, tr("Error saving options"),*warns);
    }

    files.close();
    SetCurrentFile(fullFileName,true,false);
}

//    Saves options only clase 0: Project
void MainWindow::saveOptions_0(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("PROJECTSECT",file);

    write_option("ImportFolder", "PROJECTSECT", Path(TXTImport->text()), file, printwarns, warns);
    write_option("ImportFile", "PROJECTSECT", FileWithoutPath(TXTImport->text()), file, printwarns, warns);
    write_option("ProjectFolder", "PROJECTSECT", TXTProjectFolder->text(), file, printwarns, warns);
    write_option("ProjectName", "PROJECTSECT", TXTProjectName->text(), file, printwarns, warns);
}

//    Saves options clase 0 || clase 1: Project || G-DAM || DAM
void MainWindow::saveOptions_1(string file, bool *printwarns, QString *warns)
{
    QString qv;
    CIniFile::DeleteSection("DAMSECT",file);
    CIniFile::DeleteSection("G-DAMSECT",file);

    if (iswindows){
        write_option("iswindows", "DAMSECT", QString("T"), file, printwarns, warns);
        write_option("iswindows", "G-DAMSECT", QString("T"), file, printwarns, warns);
    }else{
        write_option("iswindows", "DAMSECT", QString("F"), file, printwarns, warns);
        write_option("iswindows", "G-DAMSECT", QString("F"), file, printwarns, warns);
    }

    if (lzdo){
        if (lslater){
            write_option("lzdo", "DAMSECT", QString("T"), file, printwarns, warns);
        }
        else{
            write_option("lzdo", "G-DAMSECT", QString("T"), file, printwarns, warns);
        }
    }

    if (lvalence){
        if (lslater){
            write_option("lvalence", "DAMSECT", QString("T"), file, printwarns, warns);
        }
        else{
            write_option("lvalence", "G-DAMSECT", QString("T"), file, printwarns, warns);
        }
    }

    if (RBTatdensD1center->isChecked()) qv = QString("2");
    else if (RBTatdensD2center->isChecked())qv = QString("3");
    else qv = QString("1");

    if (lslater){
        write_option("lmaxexp", "DAMSECT", SPBatdenslmaxexp->text(), file, printwarns, warns);
        write_option("lmultmx", "DAMSECT", SPBatdenslmaxdisp->text(), file, printwarns, warns);
        write_option("ioptaj", "DAMSECT", qv, file, printwarns, warns);
        write_option("umbral", "DAMSECT", QString("1.d"+SPBfradthreshold->text()), file, printwarns, warns);
        write_option("umbralres", "DAMSECT", QString("1.d"+SPBfitthreshold->text()), file, printwarns, warns);
    }
    else{
        write_option("lmaxexp", "G-DAMSECT", SPBatdenslmaxexp->text(), file, printwarns, warns);
        write_option("lmultmx", "G-DAMSECT", SPBatdenslmaxdisp->text(), file, printwarns, warns);
        write_option("ioptaj", "G-DAMSECT", qv, file, printwarns, warns);
        write_option("umbral", "G-DAMSECT", QString("1.d"+SPBfradthreshold->text()), file, printwarns, warns);
        write_option("umbralres", "G-DAMSECT", QString("1.d"+SPBfitthreshold->text()), file, printwarns, warns);
    }

}

//    Saves options clase 0 || clase 2: Project || DAMDEN
void MainWindow::saveOptions_2(string file, bool *printwarns, QString *warns)
{

    CIniFile::DeleteSection("DAMDENSECT",file);

    if (iswindows) write_option("iswindows", "DAMDENSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMDENSECT", QString("F"), file, printwarns, warns);

    if (RBTdensExact->isChecked()) {
        write_option("lexact", "DAMDENSECT", QString("T"), file, printwarns, warns);
        write_option("ldeform", "DAMDENSECT", QString("F"), file, printwarns, warns);
        write_option("lgradient", "DAMDENSECT", QString("F"), file, printwarns, warns);
        write_option("lderiv2", "DAMDENSECT", QString("F"), file, printwarns, warns);
        write_option("lmolec", "DAMDENSECT", QString("F"), file, printwarns, warns);
        write_option("latomics", "DAMDENSECT", QString("F"), file, printwarns, warns);
        write_option("ldensacc", "DAMDENSECT", QString("F"), file, printwarns, warns);
        write_option("latomsel", "DAMDENSECT", QString("F"), file, printwarns, warns);
    }
    else {
        write_option("lexact", "DAMDENSECT", QString("F"), file, printwarns, warns);
        if (RBTdensdeform->isChecked()){
            write_option("ldeform", "DAMDENSECT", QString("T"), file, printwarns, warns);
        }
        else{
            write_option("ldeform", "DAMDENSECT", QString("F"), file, printwarns, warns);
        }
        write_option("lmaxrep", "DAMDENSECT", SPBdenslmaxexp->text(), file, printwarns, warns);
        if (RBTdensdeform->isChecked()){
            write_option("lminrep", "DAMDENSECT", QString("1"), file, printwarns, warns);
        }
        else{
            write_option("lminrep", "DAMDENSECT", SPBdenslminexp->text(), file, printwarns, warns);
        }


        if (CHKdensgrad->isChecked()) write_option("lgradient", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("lgradient", "DAMDENSECT", QString("F"), file, printwarns, warns);

        if (CHKdensder2->isChecked()) write_option("lderiv2", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("lderiv2", "DAMDENSECT", QString("F"), file, printwarns, warns);

        if (CHKdenslaplacian->isChecked()) write_option("laplacian", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("laplacian", "DAMDENSECT", QString("F"), file, printwarns, warns);

        if (CHKdenslmolec->isChecked()) write_option("lmolec", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("lmolec", "DAMDENSECT", QString("F"), file, printwarns, warns);

        if (CHKdenslatomics->isChecked())
            write_option("latomics", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("latomics", "DAMDENSECT", QString("F"), file, printwarns, warns);

        if (CHKdensldensacc->isChecked()) write_option("ldensacc", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("ldensacc", "DAMDENSECT", QString("F"), file, printwarns, warns);

        if (CHKdenslatomics->isChecked() || CHKdensldensacc->isChecked())
            write_option("latomsel", "DAMDENSECT", QString("T"), file, printwarns, warns);
        else write_option("latomsel", "DAMDENSECT", QString("F"), file, printwarns, warns);
    }

    if (denslist->count() < 1){
        write_option("nsel", "DAMDENSECT", QString("0"), file, printwarns, warns);
    }
    else{
        QList <int> list;
        for (int i = 0 ; i < denslist->count() ; ++i){
            list.append(denslist->at(i).toInt());
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        int ncntab = 0;
        for (int k=0 ; k < list.count() ; k++){
            if (list[k] <= 0 ) continue;
            ncntab++;
            write_option(QString("%1%2%3").arg("iatomsel(").arg(ncntab).arg(")").toStdString().c_str(),
                "DAMDENSECT", QString::number(list[k]), file, printwarns, warns);
        }
        write_option("nsel", "DAMDENSECT", QString("%1").arg(ncntab), file, printwarns, warns);
    }

    write_option("filename", "DAMDENSECT", TXTdensdamdenfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    write_intervals("xinf", "xsup", "dltx", "DAMDENSECT", file, TXTdensxinf, TXTdensxsup, printwarns, warns);
    write_intervals("yinf", "ysup", "dlty", "DAMDENSECT", file, TXTdensyinf, TXTdensysup, printwarns, warns);
    write_intervals("zinf", "zsup", "dltz", "DAMDENSECT", file, TXTdenszinf, TXTdenszsup, printwarns, warns);
    write_intervals("uinf", "usup", "dltu", "DAMDENSECT", file, TXTdensuinf, TXTdensusup, printwarns, warns);
    write_intervals("vinf", "vsup", "dltv", "DAMDENSECT", file, TXTdensvinf, TXTdensvsup, printwarns, warns);

    if (RBTdens2D->isChecked()) write_option("lgrid2d", "DAMDENSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid2d", "DAMDENSECT", QString("F"), file, printwarns, warns);
    write_option("x_func_uv", "DAMDENSECT", QString("\"")+TXTdensxformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("y_func_uv", "DAMDENSECT", QString("\"")+TXTdensyformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("z_func_uv", "DAMDENSECT", QString("\"")+TXTdenszformula2D->text()+QString("\""), file, printwarns, warns);

    write_option("lboundsx", "DAMDENSECT", QString("F"), file, printwarns, warns);
    write_option("xboundinf", "DAMDENSECT", QString("0.0"), file, printwarns, warns);
    write_option("xboundsup", "DAMDENSECT", QString("0.0"), file, printwarns, warns);
    write_option("lboundsy", "DAMDENSECT", QString("F"), file, printwarns, warns);
    write_option("yboundinf", "DAMDENSECT", QString("0.0"), file, printwarns, warns);
    write_option("yboundsup", "DAMDENSECT", QString("0.0"), file, printwarns, warns);
    write_option("lboundsz", "DAMDENSECT", QString("F"), file, printwarns, warns);
    write_option("zboundinf", "DAMDENSECT", QString("0.0"), file, printwarns, warns);
    write_option("zboundsup", "DAMDENSECT", QString("0.0"), file, printwarns, warns);

    write_option("planecase", "DAMDENSECT", QString("%1").arg(densplanecase).toStdString().c_str(), file, printwarns, warns);

    write_option("planeA", "DAMDENSECT", TXTdensplaneA->text(), file, printwarns, warns);
    write_option("planeB", "DAMDENSECT", TXTdensplaneB->text(), file, printwarns, warns);
    write_option("planeC", "DAMDENSECT", TXTdensplaneC->text(), file, printwarns, warns);

    write_option("numrtab", "DAMDENSECT", QString::number(SHTxyz->tabla->rowCount()), file, printwarns, warns);

    for(int i = 0 ; i < SHTxyz->tabla->rowCount() ; ++i){
        write_option(QString("%1%2%3").arg("rtab(1,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMDENSECT", SHTxyz->getcellvalue(i,0), file, printwarns, warns);
        write_option(QString("%1%2%3").arg("rtab(2,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMDENSECT", SHTxyz->getcellvalue(i,1), file, printwarns, warns);
        write_option(QString("%1%2%3").arg("rtab(3,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMDENSECT", SHTxyz->getcellvalue(i,2), file, printwarns, warns);
    }
    for(int i = SHTxyz->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(1,").arg(i+1).arg(")")),"DAMDENSECT", file);
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(2,").arg(i+1).arg(")")),"DAMDENSECT", file);
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(3,").arg(i+1).arg(")")),"DAMDENSECT", file);
    }

    if (CHKdensgrid->isChecked()) write_option("lgrid", "DAMDENSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid", "DAMDENSECT", QString("F"), file, printwarns, warns);

    if (CHKdensxyz->isChecked()) write_option("lpoints", "DAMDENSECT", QString("T"), file, printwarns, warns);
    else write_option("lpoints", "DAMDENSECT", QString("F"), file, printwarns, warns);

}

//    Saves options clase 0 || clase 3: Project || DAMPOT
void MainWindow::saveOptions_3(string file, bool *printwarns, QString *warns)
{

    CIniFile::DeleteSection("DAMPOTSECT",file);

    if (iswindows) write_option("iswindows", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMPOTSECT", QString("F"), file, printwarns, warns);

    if (lvalence){
        write_option("lvalence", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    }

    write_option("filename", "DAMPOTSECT", TXTpotgdampotfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    write_option("lmaxrep", "DAMPOTSECT", SPBpotlmaxexp->text(), file, printwarns, warns);

    if (CHKpotlong->isChecked()) write_option("largo", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("largo", "DAMPOTSECT", QString("F"), file, printwarns, warns);

    write_option("umbrlargo", "DAMPOTSECT", QString("1.d"+SPBpotlongthreshold->text()), file, printwarns, warns);

    if (CHKpotexact->isChecked()) write_option("lexact", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("lexact", "DAMPOTSECT", QString("F"), file, printwarns, warns);

    if (CHKpotgrad->isChecked()) write_option("lgradient", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("lgradient", "DAMPOTSECT", QString("F"), file, printwarns, warns);

    if (CHKpotder2->isChecked()) write_option("lderiv2", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("lderiv2", "DAMPOTSECT", QString("F"), file, printwarns, warns);


    write_intervals("xinf", "xsup", "dltx", "DAMPOTSECT", file, TXTpotxinf, TXTpotxsup, printwarns, warns);
    write_intervals("yinf", "ysup", "dlty", "DAMPOTSECT", file, TXTpotyinf, TXTpotysup, printwarns, warns);
    write_intervals("zinf", "zsup", "dltz", "DAMPOTSECT", file, TXTpotzinf, TXTpotzsup, printwarns, warns);
    write_intervals("uinf", "usup", "dltu", "DAMPOTSECT", file, TXTpotuinf, TXTpotusup, printwarns, warns);
    write_intervals("vinf", "vsup", "dltv", "DAMPOTSECT", file, TXTpotvinf, TXTpotvsup, printwarns, warns);

    if (RBTpot2D->isChecked()) write_option("lgrid2d", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid2d", "DAMPOTSECT", QString("F"), file, printwarns, warns);

    write_option("x_func_uv", "DAMPOTSECT", QString("\"")+TXTpotxformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("y_func_uv", "DAMPOTSECT", QString("\"")+TXTpotyformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("z_func_uv", "DAMPOTSECT", QString("\"")+TXTpotzformula2D->text()+QString("\""), file, printwarns, warns);


    write_option("planecase", "DAMPOTSECT", QString("%1").arg(potplanecase).toStdString().c_str(), file, printwarns, warns);

    write_option("planeA", "DAMPOTSECT", TXTpotplaneA->text(), file, printwarns, warns);
    write_option("planeB", "DAMPOTSECT", TXTpotplaneB->text(), file, printwarns, warns);
    write_option("planeC", "DAMPOTSECT", TXTpotplaneC->text(), file, printwarns, warns);

    write_option("numrtab", "DAMPOTSECT", QString::number(SHTpotxyz->tabla->rowCount()), file, printwarns, warns);

    for(int i = 0 ; i < SHTpotxyz->tabla->rowCount() ; ++i){
        write_option(QString("%1%2%3").arg("rtab(1,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMPOTSECT", SHTpotxyz->getcellvalue(i,0), file, printwarns, warns);
        write_option(QString("%1%2%3").arg("rtab(2,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMPOTSECT", SHTpotxyz->getcellvalue(i,1), file, printwarns, warns);
        write_option(QString("%1%2%3").arg("rtab(3,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMPOTSECT", SHTpotxyz->getcellvalue(i,2), file, printwarns, warns);
    }
    for(int i = SHTpotxyz->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(1,").arg(i+1).arg(")")),"DAMPOTSECT", file);
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(2,").arg(i+1).arg(")")),"DAMPOTSECT", file);
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(3,").arg(i+1).arg(")")),"DAMPOTSECT", file);
    }

    if (CHKpotgrid->isChecked()) write_option("lgrid", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid", "DAMPOTSECT", QString("F"), file, printwarns, warns);

    if (CHKpotxyz->isChecked()) write_option("lpoints", "DAMPOTSECT", QString("T"), file, printwarns, warns);
    else write_option("lpoints", "DAMPOTSECT", QString("F"), file, printwarns, warns);
}

//    Saves options clase 0 || clase 4: Project || DAMFORCES
void MainWindow::saveOptions_4(string file, bool *printwarns, QString *warns)
{

    CIniFile::DeleteSection("DAMFORCESSECT",file);

    if (iswindows) write_option("iswindows", "DAMFORCESSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMFORCESSECT", QString("F"), file, printwarns, warns);

    if (lvalence){
        write_option("lvalence", "DAMFORCESSECT", QString("T"), file, printwarns, warns);
    }

    if (CHKHFlatomsel->isChecked())
        write_option("latomsel", "DAMFORCESSECT", QString("T"), file, printwarns, warns);
    else write_option("latomsel", "DAMFORCESSECT", QString("F"), file, printwarns, warns);

    if (HFforceslist->count() < 1){
        write_option("ncntab", "DAMFORCESSECT", QString("0"), file, printwarns, warns);
    }
    else{
        QList <int> list;
        for (int i = 0 ; i < HFforceslist->count() ; ++i){
            list.append(HFforceslist->at(i).toInt());
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        int ncntab = 0;
        for (int k = 0 ; k < list.count() ; k++){
            if (list[k] <= 0 ) continue;
            ncntab++;
            write_option(QString("%1%2%3").arg("iatomsel(").arg(ncntab).arg(")").toStdString().c_str(),
                "DAMFORCESSECT", QString::number(list[k]), file, printwarns, warns);
        }
        write_option("ncntab", "DAMFORCESSECT", QString("%1").arg(ncntab), file, printwarns, warns);
    }
    write_option("filename", "DAMFORCESSECT", TXTHFgdamforcesfilename->text().remove("\"").prepend("\"").append("\""),
        file, printwarns, warns);
}

//    Saves options clase 0 || clase 5: Project || DAMFIELD
void MainWindow::saveOptions_5(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("DAMFIELDSECT",file);

    if (lvalence){
        write_option("lvalence", "DAMFIELDSECT", QString("T"), file, printwarns, warns);
    }

    if (iswindows) write_option("iswindows", "DAMFIELDSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMFIELDSECT", QString("F"), file, printwarns, warns);

    write_option("filename", "DAMFIELDSECT", TXTeffilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);
    write_option("lmaxrep", "DAMFIELDSECT", SPBeflmaxexp->text(), file, printwarns, warns);

    if (CHKeflong->isChecked()) write_option("largo", "DAMFIELDSECT", QString("T"), file, printwarns, warns);
    else write_option("largo", "DAMFIELDSECT", QString("F"), file, printwarns, warns);

    write_option("umbrlargo", "DAMFIELDSECT", QString("1.d"+SPBeflongthreshold->text()), file, printwarns, warns);

    write_option("numpnt", "DAMFIELDSECT", TXTefnumpnt->text(), file, printwarns, warns);
    write_option("nlinpernuc", "DAMFIELDSECT", TXTefnlinpernuc->text(), file, printwarns, warns);
    write_option("ioplines3D", "DAMFIELDSECT", QString("%1").arg(CMBefdirectionset->currentIndex()), file, printwarns, warns);

    if (!TXTefdlt0->text().isEmpty()) write_option("dlt0", "DAMFIELDSECT", TXTefdlt0->text(), file, printwarns, warns);
    else write_option("dlt0", "DAMFIELDSECT", "0.02", file, printwarns, warns);

    write_option("umbrlargo", "DAMFIELDSECT", QString("1.d"+SPBeflongthreshold->text()), file, printwarns, warns);

    if (TXTefuinf->text().toDouble() > TXTefusup->text().toDouble()){
        QString aux = TXTefusup->text();
        TXTefusup->setText(TXTefuinf->text());
        TXTefuinf->setText(aux);
    }
    write_option("uinf", "DAMFIELDSECT", TXTefuinf->text(), file, printwarns, warns);
    write_option("usup", "DAMFIELDSECT", TXTefusup->text(), file, printwarns, warns);

    if (TXTefvinf->text().toDouble() > TXTefvsup->text().toDouble()){
        QString aux = TXTefvsup->text();
        TXTefvsup->setText(TXTefvinf->text());
        TXTefvinf->setText(aux);
    }
    write_option("vinf", "DAMFIELDSECT", TXTefvinf->text(), file, printwarns, warns);
    write_option("vsup", "DAMFIELDSECT", TXTefvsup->text(), file, printwarns, warns);

    if (TXTefxinf->text().toDouble() > TXTefxsup->text().toDouble()){
        QString aux = TXTefxsup->text();
        TXTefxsup->setText(TXTefxinf->text());
        TXTefxinf->setText(aux);
    }
    write_option("xinf", "DAMFIELDSECT", TXTefxinf->text(), file, printwarns, warns);
    write_option("xsup", "DAMFIELDSECT", TXTefxsup->text(), file, printwarns, warns);

    if (TXTefyinf->text().toDouble() > TXTefysup->text().toDouble()){
        QString aux = TXTefysup->text();
        TXTefysup->setText(TXTefyinf->text());
        TXTefyinf->setText(aux);
    }
    write_option("yinf", "DAMFIELDSECT", TXTefyinf->text(), file, printwarns, warns);
    write_option("ysup", "DAMFIELDSECT", TXTefysup->text(), file, printwarns, warns);
    if (TXTefzinf->text().toDouble() > TXTefzsup->text().toDouble()){
        QString aux = TXTefzsup->text();
        TXTefzsup->setText(TXTefzinf->text());
        TXTefzinf->setText(aux);
    }
    write_option("zinf", "DAMFIELDSECT", TXTefzinf->text(), file, printwarns, warns);
    write_option("zsup", "DAMFIELDSECT", TXTefzsup->text(), file, printwarns, warns);

    if (RBTef2D->isChecked()) write_option("lplot2d", "DAMFIELDSECT", QString("T"), file, printwarns, warns);
    else write_option("lplot2d", "DAMFIELDSECT", QString("F"), file, printwarns, warns);

    write_option("planeA", "DAMFIELDSECT", TXTefplaneA->text(), file, printwarns, warns);
    write_option("planeB", "DAMFIELDSECT", TXTefplaneB->text(), file, printwarns, warns);
    write_option("planeC", "DAMFIELDSECT", TXTefplaneC->text(), file, printwarns, warns);
    write_option("uvratio", "DAMFIELDSECT", TXTefuvratio->text(), file, printwarns, warns);

    if (CHKefextralines->isChecked()) write_option("lextralines", "DAMFIELDSECT", QString("T"), file, printwarns, warns);
    else write_option("lextralines", "DAMFIELDSECT", QString("F"), file, printwarns, warns);

    if (CHKefxyz->isChecked()){
        if (RBTef3D->isChecked()){
            write_option("nlines", "DAMFIELDSECT", QString::number(SHTefxyz->tabla->rowCount()), file, printwarns, warns);

            for(int i = 0 ; i < SHTefxyz->tabla->rowCount() ; ++i){
                write_option(QString("%1%2%3").arg("icntlines(").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", QString("%1").arg((SHTefxyz->getcellvalue(i,0)).toInt()), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", SHTefxyz->getcellvalue(i,1), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", SHTefxyz->getcellvalue(i,2), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(3,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", SHTefxyz->getcellvalue(i,3), file, printwarns, warns);
            }
            for(int i = SHTefxyz->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")")),"DAMFIELDSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")")),"DAMFIELDSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(3,").arg(i+1).arg(")")),"DAMFIELDSECT", file);
            }
        }
        else{
            write_option("nlines", "DAMFIELDSECT", QString::number(SHTefuv->tabla->rowCount()), file, printwarns, warns);

            for(int i = 0 ; i < SHTefuv->tabla->rowCount() ; ++i){
                write_option(QString("%1%2%3").arg("icntlines(").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", QString("%1").arg((SHTefuv->getcellvalue(i,0)).toInt()), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", SHTefuv->getcellvalue(i,1), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMFIELDSECT", SHTefuv->getcellvalue(i,2), file, printwarns, warns);
            }
            for(int i = SHTefuv->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")")),"DAMFIELDSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")")),"DAMFIELDSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(3,").arg(i+1).arg(")")),"DAMFIELDSECT", file);
            }
        }
    }
    else{
        write_option("nlines", "DAMFIELDSECT","0", file, printwarns, warns);
    }

    write_option("filelines", "DAMFIELDSECT", TXTeffilelines->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);
}

//    Saves options clase 0 || clase 6: Project || DAMFRAD
void MainWindow::saveOptions_6(string file, bool *printwarns, QString *warns)
{
    bool ok;
    double dv;
    int i;
    CIniFile::DeleteSection("DAMFRADSECT",file);

    if (iswindows) write_option("iswindows", "DAMFRADSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMFRADSECT", QString("F"), file, printwarns, warns);

    write_option("rini", "DAMFRADSECT", TXTfradrini->text(), file, printwarns, warns);
    write_option("rfin", "DAMFRADSECT", TXTfradrfin->text(), file, printwarns, warns);
    write_option("dltr", "DAMFRADSECT", TXTfradltr->text(), file, printwarns, warns);
    write_option("ltab", "DAMFRADSECT", SPBfradltab->text(), file, printwarns, warns);
    write_option("mtab", "DAMFRADSECT", SPBfradmtab->text(), file, printwarns, warns);

    if (CHKfradextras->isChecked()) write_option("lrlist", "DAMFRADSECT", QString("T"), file, printwarns, warns);
    else write_option("lrlist", "DAMFRADSECT", QString("F"), file, printwarns, warns);

    if (CHKfradextras->isChecked()){
        QList <double> list;
        int knt = 0;
        for(int i = 0 ; i < SHTfradrlist->tabla->rowCount() ; ++i){
            dv=SHTfradrlist->getcellvalue(i,0).toDouble(&ok);
            if (!ok || dv < 0 || list.contains(dv)) continue;
            list << dv;
            knt++;
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        for(int i = 0 ; i < knt ; ++i){
            write_option(QString("%1%2%3").arg("rlist(").arg(i+1).arg(")").toStdString().c_str(),
                "DAMFRADSECT", QString::number(list[i]), file, printwarns, warns);
        }
        for(int i =knt+1 ; i < Sheet::max_sel ; ++i){
            CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlist(").arg(i).arg(")")),"DAMFRADSECT", file);
        }
        write_option("nlist", "DAMFRADSECT", QString("%1").arg(knt), file, printwarns, warns);

    }else{
        write_option("nlist", "DAMFRADSECT", QString("0"), file, printwarns, warns);
    }

    if (fradlist->count() < 1){
        write_option("ncntab", "DAMFRADSECT", QString("0"), file, printwarns, warns);
    }
    else{
        QList <int> list;
        for (i = 0 ; i < fradlist->count() ; ++i){
            list.append(fradlist->at(i).toInt());
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        int ncntab = 0;
        for (int k=0 ; k < list.count() ; k++){
            if (list[k] <= 0 ) continue;
            ncntab++;
            write_option(QString("%1%2%3").arg("iatomsel(").arg(ncntab).arg(")").toStdString().c_str(),
                "DAMFRADSECT", QString::number(list[k]), file, printwarns, warns);
        }
        write_option("ncntab", "DAMFRADSECT", QString("%1").arg(ncntab), file, printwarns, warns);
    }

    write_option("filename", "DAMFRADSECT", TXTfraddamfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);


    if (CHKfradderiv1->isChecked()) write_option("lderiv", "DAMFRADSECT", QString("T"), file, printwarns, warns);
    else write_option("lderiv", "DAMFRADSECT", QString("F"), file, printwarns, warns);

    if (CHKfradderiv2->isChecked()) write_option("lderiv2", "DAMFRADSECT", QString("T"), file, printwarns, warns);
    else write_option("lderiv2", "DAMFRADSECT", QString("F"), file, printwarns, warns);
}

//    Saves options clase 0 || clase 7: Project || DAMMULTROT
void MainWindow::saveOptions_7(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("DAMMULTROTSECT",file);

    if (iswindows) write_option("iswindows", "DAMMULTROTSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMMULTROTSECT", QString("F"), file, printwarns, warns);

    write_option("filename", "DAMMULTROTSECT", TXTmrotorimultfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    write_option("lmin", "DAMMULTROTSECT", SPBmrotlmin->text(), file, printwarns, warns);
    write_option("lmax", "DAMMULTROTSECT", SPBmrotlmax->text(), file, printwarns, warns);
    write_option("i1", "DAMMULTROTSECT", SPBmrotleft->text(), file, printwarns, warns);
    write_option("i2", "DAMMULTROTSECT", SPBmrotmiddle->text(), file, printwarns, warns);
    write_option("i3", "DAMMULTROTSECT", SPBmrotright->text(), file, printwarns, warns);

    if (mrotorimultlist->count() < 1){
        write_option("ncntab", "DAMMULTROTSECT", QString("0"), file, printwarns, warns);
    }
    else{
        QList <int> list;
        for (int i = 0 ; i < mrotorimultlist->count() ; ++i){
            list.append(mrotorimultlist->at(i).toInt());
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        int ncntab = 0;
        for (int k=0 ; k < list.count() ; k++){
            if (list[k] <= 0 ) continue;
            ncntab++;
            write_option(QString("%1%2%3").arg("icntab(").arg(ncntab).arg(")").toStdString().c_str(),
                "DAMMULTROTSECT", QString::number(list[k]), file, printwarns, warns);
        }
        write_option("ncntab", "DAMMULTROTSECT", QString("%1").arg(ncntab), file, printwarns, warns);
    }
}

//    Saves options clase 0 || clase 8: Project || DAMORB
void MainWindow::saveOptions_8(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("DAMORBSECT",file);

    if (iswindows) write_option("iswindows", "DAMORBSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMORBSECT", QString("F"), file, printwarns, warns);

    write_option("filename", "DAMORBSECT", TXTMOfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);
    write_option("fileMOname", "DAMORBSECT", TXTMOImportfile->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    if (CHKMOgrad->isChecked()) write_option("lgradient", "DAMORBSECT", QString("T"), file, printwarns, warns);
    else write_option("lgradient", "DAMORBSECT", QString("F"), file, printwarns, warns);

    write_intervals("xinf", "xsup", "dltx", "DAMORBSECT", file, TXTMOxinf, TXTMOxsup, printwarns, warns);
    write_intervals("yinf", "ysup", "dlty", "DAMORBSECT", file, TXTMOyinf, TXTMOysup, printwarns, warns);
    write_intervals("zinf", "zsup", "dltz", "DAMORBSECT", file, TXTMOzinf, TXTMOzsup, printwarns, warns);
    write_intervals("uinf", "usup", "dltu", "DAMORBSECT", file, TXTMOuinf, TXTMOusup, printwarns, warns);
    write_intervals("vinf", "vsup", "dltv", "DAMORBSECT", file, TXTMOvinf, TXTMOvsup, printwarns, warns);

    if (RBTMO2D->isChecked()) write_option("lgrid2d", "DAMORBSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid2d", "DAMORBSECT", QString("F"), file, printwarns, warns);

    write_option("x_func_uv", "DAMORBSECT", QString("\"")+TXTMOxformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("y_func_uv", "DAMORBSECT", QString("\"")+TXTMOyformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("z_func_uv", "DAMORBSECT", QString("\"")+TXTMOzformula2D->text()+QString("\""), file, printwarns, warns);

    write_option("planecase", "DAMORBSECT", QString("%1").arg(MOplanecase).toStdString().c_str(), file, printwarns, warns);

    write_option("planeA", "DAMORBSECT", TXTMOplaneA->text(), file, printwarns, warns);
    write_option("planeB", "DAMORBSECT", TXTMOplaneB->text(), file, printwarns, warns);
    write_option("planeC", "DAMORBSECT", TXTMOplaneC->text(), file, printwarns, warns);

    if (MOlist->count() < 1){
        write_option("norbs", "DAMORBSECT", QString("1"), file, printwarns, warns);
        write_option("iorbsinp(1)", "DAMORBSECT", QString("1"), file, printwarns, warns);
    }
    else{
        QList <int> list;
        for (int i = 0 ; i < MOlist->count() ; ++i){
            list.append(MOlist->at(i).toInt());
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        int norbs = 0;
        for (int k=0 ; k < list.count() ; k++){
            if (list[k] <= 0 ) continue;
            norbs++;
            write_option(QString("%1%2%3").arg("iorbsinp(").arg(norbs).arg(")").toStdString().c_str(),
                "DAMORBSECT", QString::number(list[k]), file, printwarns, warns);
        }
        write_option("norbs", "DAMORBSECT", QString("%1").arg(norbs), file, printwarns, warns);
    }
}

//    Saves options clase 0 || clase 9: Project || DAMTOPO
void MainWindow::saveOptions_9(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("DAMTOPOSECT",file);

    if (iswindows) write_option("iswindows", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    if (lvalence){
        write_option("lvalence", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    }

    if (CHKtopograph->isChecked()) write_option("TOPOGRAPH", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("TOPOGRAPH", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    write_option("filename", "DAMTOPOSECT", TXTtopofilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    if (CHKtopoaddguess->isChecked()) write_option("addguess", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("addguess", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    if (CHKtopoxyz->isChecked()){
        write_option("ncntguess", "DAMTOPOSECT", QString::number(SHTtopoxyz->tabla->rowCount()), file, printwarns, warns);

        for(int i = 0 ; i < SHTtopoxyz->tabla->rowCount() ; ++i){
            write_option(QString("%1%2%3").arg("rcntguess(1,").arg(i+1).arg(")").toStdString().c_str(),
                    "DAMTOPOSECT", QString(SHTtopoxyz->getcellvalue(i,0)), file, printwarns, warns);
            write_option(QString("%1%2%3").arg("rcntguess(2,").arg(i+1).arg(")").toStdString().c_str(),
                    "DAMTOPOSECT", SHTtopoxyz->getcellvalue(i,1), file, printwarns, warns);
            write_option(QString("%1%2%3").arg("rcntguess(3,").arg(i+1).arg(")").toStdString().c_str(),
                    "DAMTOPOSECT", SHTtopoxyz->getcellvalue(i,2), file, printwarns, warns);
        }
        for(int i = SHTtopoxyz->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
            CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rcntguess(1,").arg(i+1).arg(")")),"DAMTOPOSECT", file);
            CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rcntguess(2,").arg(i+1).arg(")")),"DAMTOPOSECT", file);
            CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rcntguess(3,").arg(i+1).arg(")")),"DAMTOPOSECT", file);
        }
    }
    else{
        write_option("ncntguess", "DAMTOPOSECT","0", file, printwarns, warns);
    }

    write_option("guessfile", "DAMTOPOSECT", TXTtopoguessfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    write_option("BOXT", "DAMTOPOSECT", TXTtopoboxt->text(), file, printwarns, warns);
    write_option("stepszt", "DAMTOPOSECT", TXTtopostepszt->text(), file, printwarns, warns);

    write_option("BOXL", "DAMTOPOSECT", TXTtopoboxl->text(), file, printwarns, warns);

    if (RBTtopodensity->isChecked()) write_option("MED", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("MED", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    if (RBTtopopotential->isChecked()) write_option("MESP", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("MESP", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    write_option("CNVG", "DAMTOPOSECT", TXTtopocnvg->text(), file, printwarns, warns);
    write_option("lmaxi", "DAMTOPOSECT", SPBtopolmaxi->text(), file, printwarns, warns);

    if (CHKtopomolgraph->isChecked()) write_option("gradpath", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("gradpath", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    write_option("BOXG", "DAMTOPOSECT", TXTtopoboxg->text(), file, printwarns, warns);
    write_option("DRCUTCP", "DAMTOPOSECT", TXTtopoggradthr->text(), file, printwarns, warns);

    if (CHKtopobasin->isChecked()) write_option("basin", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("basin", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    write_option("BOXB", "DAMTOPOSECT", TXTtopoboxb->text(), file, printwarns, warns);

    write_option("FDISP", "DAMTOPOSECT", TXTtopofdisp->text(), file, printwarns, warns);

    write_option("wireframe", "DAMTOPOSECT", QString("F"), file, printwarns, warns);
    write_option("solidsurf", "DAMTOPOSECT", QString("F"), file, printwarns, warns);

    if (CHKtopoexdraw->isChecked()) write_option("exdraw", "DAMTOPOSECT", QString("T"), file, printwarns, warns);
    else write_option("exdraw", "DAMTOPOSECT", QString("F"), file, printwarns, warns);
    write_option("exln", "DAMTOPOSECT", TXTtopoexln->text(), file, printwarns, warns);
}

//    Saves options clase 0 || clase 10: Project || DAMZJ
void MainWindow::saveOptions_10(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("DAMZJSECT",file);

    if (iswindows){
        write_option("iswindows", "DAMZJSECT", QString("T"), file, printwarns, warns);
    }else{
        write_option("iswindows", "DAMZJSECT", QString("F"), file, printwarns, warns);
    }

    write_option("rstar", "DAMZJSECT", QString("%1").arg(rstar), file, printwarns, warns);
    if (RBTZJrstarabs->isChecked())
        write_option("lrstarrel", "DAMZJSECT", QString("F"), file, printwarns, warns);
    else
        write_option("lrstarrel", "DAMZJSECT", QString("T"), file, printwarns, warns);

    write_option("lexpansion", "DAMZJSECT", SPBZJlmax->text(), file, printwarns, warns);
    write_option("kexpansion", "DAMZJSECT", SPBZJkmax->text(), file, printwarns, warns);

    if (RBTZJacobi->isChecked()) write_option("ljacobi", "DAMZJSECT", QString("T"), file, printwarns, warns);
    else write_option("ljacobi", "DAMZJSECT", QString("F"), file, printwarns, warns);

    write_option("nquadpoints", "DAMZJSECT", SPBZJnquad->text(), file, printwarns, warns);

    if (RBTZJechelon->isChecked()) write_option("lechelon", "DAMZJSECT", QString("T"), file, printwarns, warns);
    else write_option("lechelon", "DAMZJSECT", QString("F"), file, printwarns, warns);

    write_option("thresmult", "DAMZJSECT", QString("1.d"+SPBZJthrmult->text()), file, printwarns, warns);
    write_option("thresoverlap", "DAMZJSECT", QString("1.d"+SPBZJthrdist->text()), file, printwarns, warns);

    if (lzdo){
            write_option("lzdo", "DAMZJSECT", QString("T"), file, printwarns, warns);
    }

    if (lvalence){
            write_option("lvalence", "DAMZJSECT", QString("T"), file, printwarns, warns);
    }
}

//    Saves options clase 0 || clase 11: Project || DAMDENZJ
void MainWindow::saveOptions_11(string file, bool *printwarns, QString *warns)
{

    CIniFile::DeleteSection("DAMDENZJSECT",file);

    if (iswindows) write_option("iswindows", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMDENZJSECT", QString("F"), file, printwarns, warns);

//        ljacobi and lechelon are taken from Zernike-Jacobi module (DAMdZJ)
    if (RBTZJechelon->isChecked()) write_option("lechelon", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
    else write_option("lechelon", "DAMDENZJSECT", QString("F"), file, printwarns, warns);



    if (CHKdZJgrad->isChecked()) write_option("lgradient", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
    else write_option("lgradient", "DAMDENZJSECT", QString("F"), file, printwarns, warns);

    write_option("filename", "DAMDENZJSECT", TXTdZJfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);
    write_option("fileZJname", "DAMDENZJSECT", TXTdZJImportfile->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);

    if (RBTdZJchoosel->isChecked() || RBTdZJchoosek->isChecked()){    // If only l values or only k values are chosen, sorts them
        QList <int> list;
        for (int i = 0 ; i < ZJlist->count() ; ++i){
            list.append(ZJlist->at(i).toInt());
        }
//        qSort(list.begin(), list.end());
        std::sort(list.begin(), list.end());
        ZJlist->clear();
        for (int k=0 ; k < list.count() ; k++){
        ZJlist->append(QString::number(list[k]));
        }
    }
    if (!RBTdZJchooseall->isChecked() && (ZJlist->count() > 0)){
        QString str = ZJlist->join(",");
        write_option("indices", "DAMDENZJSECT", str, file, printwarns, warns);
        write_option("kmaxrep", "DAMDENZJSECT", QString("%1").arg(MAX_KEXPZJ), file, printwarns, warns);
        write_option("lmaxrep", "DAMDENZJSECT", QString("%1").arg(MAX_LEXPZJ), file, printwarns, warns);
        write_option("lminrep", "DAMDENZJSECT", QString("0"), file, printwarns, warns);
        if (RBTdZJchoosel->isChecked()){
        write_option("lindividl", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
        write_option("lindividk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlkm", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        }
        else if (RBTdZJchoosek->isChecked()){
        write_option("lindividl", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividk", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
        write_option("lindividlk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlkm", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        }
        else if (RBTdZJchooselk->isChecked()){
        write_option("lindividl", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlk", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
        write_option("lindividlkm", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        }
        else if (RBTdZJchooselkm->isChecked()){
        write_option("lindividl", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlkm", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
        }
    }
    else{
        write_option("lindividl", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlk", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("lindividlkm", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
        write_option("kmaxrep", "DAMDENZJSECT", SPBdZJkmax->text(), file, printwarns, warns);
        write_option("lmaxrep", "DAMDENZJSECT", SPBdZJlmax->text(), file, printwarns, warns);
        write_option("lminrep", "DAMDENZJSECT", SPBdZJlmin->text(), file, printwarns, warns);
    }

    write_intervals("xinf", "xsup", "dltx", "DAMDENZJSECT", file, TXTdZJxinf, TXTdZJxsup, printwarns, warns);
    write_intervals("yinf", "ysup", "dlty", "DAMDENZJSECT", file, TXTdZJyinf, TXTdZJysup, printwarns, warns);
    write_intervals("zinf", "zsup", "dltz", "DAMDENZJSECT", file, TXTdZJzinf, TXTdZJzsup, printwarns, warns);
    write_intervals("uinf", "usup", "dltu", "DAMDENZJSECT", file, TXTdZJuinf, TXTdZJusup, printwarns, warns);
    write_intervals("vinf", "vsup", "dltv", "DAMDENZJSECT", file, TXTdZJvinf, TXTdZJvsup, printwarns, warns);

    if (RBTdZJ2D->isChecked()) write_option("lgrid2d", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid2d", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
    write_option("x_func_uv", "DAMDENZJSECT", QString("\"")+TXTdZJxformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("y_func_uv", "DAMDENZJSECT", QString("\"")+TXTdZJyformula2D->text()+QString("\""), file, printwarns, warns);
    write_option("z_func_uv", "DAMDENZJSECT", QString("\"")+TXTdZJzformula2D->text()+QString("\""), file, printwarns, warns);

    write_option("numrtab", "DAMDENZJSECT", QString::number(SHTdZJxyz->tabla->rowCount()), file, printwarns, warns);

    for(int i = 0 ; i < SHTdZJxyz->tabla->rowCount() ; ++i){
        write_option(QString("%1%2%3").arg("rtab(1,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMDENZJSECT", SHTdZJxyz->getcellvalue(i,0), file, printwarns, warns);
        write_option(QString("%1%2%3").arg("rtab(2,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMDENZJSECT", SHTdZJxyz->getcellvalue(i,1), file, printwarns, warns);
        write_option(QString("%1%2%3").arg("rtab(3,").arg(i+1).arg(")").toStdString().c_str(),
                "DAMDENZJSECT", SHTdZJxyz->getcellvalue(i,2), file, printwarns, warns);
    }
    for(int i = SHTdZJxyz->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(1,").arg(i+1).arg(")")),"DAMDENZJSECT", file);
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(2,").arg(i+1).arg(")")),"DAMDENZJSECT", file);
        CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rtab(3,").arg(i+1).arg(")")),"DAMDENZJSECT", file);
    }

    if (CHKdZJgrid->isChecked()) write_option("lgrid", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
    else write_option("lgrid", "DAMDENZJSECT", QString("F"), file, printwarns, warns);

    if (CHKdZJxyz->isChecked()) write_option("lpoints", "DAMDENZJSECT", QString("T"), file, printwarns, warns);
    else write_option("lpoints", "DAMDENZJSECT", QString("F"), file, printwarns, warns);
}

//    Saves options clase 0 || clase 12: Project || DAMDENGRAD
void MainWindow::saveOptions_12(string file, bool *printwarns, QString *warns)
{
    CIniFile::DeleteSection("DAMDENGRADSECT",file);

    if (iswindows) write_option("iswindows", "DAMDENGRADSECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMDENGRADSECT", QString("F"), file, printwarns, warns);

    write_option("filename", "DAMDENGRADSECT", TXTdgfilename->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);
    write_option("lmaxrep", "DAMDENGRADSECT", SPBdglmaxexp->text(), file, printwarns, warns);

    write_option("umbrlargo", "DAMDENGRADSECT", QString("1.d"+SPBdglongthreshold->text()), file, printwarns, warns);

    if (CHKdgextralines->isChecked()) write_option("lextralines", "DAMDENGRADSECT", QString("T"), file, printwarns, warns);
    else write_option("lextralines", "DAMDENGRADSECT", QString("F"), file, printwarns, warns);

    write_option("numpnt", "DAMDENGRADSECT", TXTdgnumpnt->text(), file, printwarns, warns);
    write_option("nlinpernuc", "DAMDENGRADSECT", TXTdgnlinpernuc->text(), file, printwarns, warns);
    write_option("ioplines3D", "DAMDENGRADSECT", QString("%1").arg(CMBdgdirectionset->currentIndex()), file, printwarns, warns);

    if (!TXTdgdlt0->text().isEmpty()) write_option("dlt0", "DAMDENGRADSECT", TXTdgdlt0->text(), file, printwarns, warns);
    else write_option("dlt0", "DAMDENGRADSECT", "0.02", file, printwarns, warns);

    write_option("umbrlargo", "DAMDENGRADSECT", QString("1.d"+SPBdglongthreshold->text()), file, printwarns, warns);

    if (TXTdguinf->text().toDouble() > TXTdgusup->text().toDouble()){
        QString aux = TXTdgusup->text();
        TXTdgusup->setText(TXTdguinf->text());
        TXTdguinf->setText(aux);
    }
    write_option("uinf", "DAMDENGRADSECT", TXTdguinf->text(), file, printwarns, warns);
    write_option("usup", "DAMDENGRADSECT", TXTdgusup->text(), file, printwarns, warns);

    if (TXTdgvinf->text().toDouble() > TXTdgvsup->text().toDouble()){
        QString aux = TXTdgvsup->text();
        TXTdgvsup->setText(TXTdgvinf->text());
        TXTdgvinf->setText(aux);
    }
    write_option("vinf", "DAMDENGRADSECT", TXTdgvinf->text(), file, printwarns, warns);
    write_option("vsup", "DAMDENGRADSECT", TXTdgvsup->text(), file, printwarns, warns);

    if (TXTdgxinf->text().toDouble() > TXTdgxsup->text().toDouble()){
        QString aux = TXTdgxsup->text();
        TXTdgxsup->setText(TXTdgxinf->text());
        TXTdgxinf->setText(aux);
    }
    write_option("xinf", "DAMDENGRADSECT", TXTdgxinf->text(), file, printwarns, warns);
    write_option("xsup", "DAMDENGRADSECT", TXTdgxsup->text(), file, printwarns, warns);

    if (TXTdgyinf->text().toDouble() > TXTdgysup->text().toDouble()){
        QString aux = TXTdgysup->text();
        TXTdgysup->setText(TXTdgyinf->text());
        TXTdgyinf->setText(aux);
    }
    write_option("yinf", "DAMDENGRADSECT", TXTdgyinf->text(), file, printwarns, warns);
    write_option("ysup", "DAMDENGRADSECT", TXTdgysup->text(), file, printwarns, warns);
    if (TXTdgzinf->text().toDouble() > TXTdgzsup->text().toDouble()){
        QString aux = TXTdgzsup->text();
        TXTdgzsup->setText(TXTdgzinf->text());
        TXTdgzinf->setText(aux);
    }
    write_option("zinf", "DAMDENGRADSECT", TXTdgzinf->text(), file, printwarns, warns);
    write_option("zsup", "DAMDENGRADSECT", TXTdgzsup->text(), file, printwarns, warns);

    if (RBTdg2D->isChecked()) write_option("lplot2d", "DAMDENGRADSECT", QString("T"), file, printwarns, warns);
    else write_option("lplot2d", "DAMDENGRADSECT", QString("F"), file, printwarns, warns);

    write_option("planeA", "DAMDENGRADSECT", TXTdgplaneA->text(), file, printwarns, warns);
    write_option("planeB", "DAMDENGRADSECT", TXTdgplaneB->text(), file, printwarns, warns);
    write_option("planeC", "DAMDENGRADSECT", TXTdgplaneC->text(), file, printwarns, warns);
    write_option("uvratio", "DAMDENGRADSECT", TXTdguvratio->text(), file, printwarns, warns);

    if (CHKdgextralines->isChecked()) write_option("lextralines", "DAMDENGRADSECT", QString("T"), file, printwarns, warns);
    else write_option("lextralines", "DAMDENGRADSECT", QString("F"), file, printwarns, warns);

    if (CHKdgxyz->isChecked()){
        if (RBTdg3D->isChecked()){
            write_option("nlines", "DAMDENGRADSECT", QString::number(SHTdgxyz->tabla->rowCount()), file, printwarns, warns);

            for(int i = 0 ; i < SHTdgxyz->tabla->rowCount() ; ++i){
                write_option(QString("%1%2%3").arg("icntlines(").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", QString("%1").arg((SHTdgxyz->getcellvalue(i,0)).toInt()), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", SHTdgxyz->getcellvalue(i,1), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", SHTdgxyz->getcellvalue(i,2), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(3,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", SHTdgxyz->getcellvalue(i,3), file, printwarns, warns);
            }
            for(int i = SHTdgxyz->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")")),"DAMDENGRADSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")")),"DAMDENGRADSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(3,").arg(i+1).arg(")")),"DAMDENGRADSECT", file);
            }
        }
        else{
            write_option("nlines", "DAMDENGRADSECT", QString::number(SHTdguv->tabla->rowCount()), file, printwarns, warns);

            for(int i = 0 ; i < SHTdguv->tabla->rowCount() ; ++i){
                write_option(QString("%1%2%3").arg("icntlines(").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", QString("%1").arg((SHTdguv->getcellvalue(i,0)).toInt()), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", SHTdguv->getcellvalue(i,1), file, printwarns, warns);
                write_option(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")").toStdString().c_str(),
                        "DAMDENGRADSECT", SHTdguv->getcellvalue(i,2), file, printwarns, warns);
            }
            for(int i = SHTdguv->tabla->rowCount() ; i<Sheet::max_sel ; ++i){
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(1,").arg(i+1).arg(")")),"DAMDENGRADSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(2,").arg(i+1).arg(")")),"DAMDENGRADSECT", file);
                CIniFile::DeleteRecord(toString(QString("%1%2%3").arg("rlines(3,").arg(i+1).arg(")")),"DAMDENGRADSECT", file);
            }
        }
    }
    else{
        write_option("nlines", "DAMDENGRADSECT","0", file, printwarns, warns);
    }

    write_option("filelines", "DAMDENGRADSECT", TXTdgfilelines->text().remove("\"").prepend("\"").append("\""), file, printwarns, warns);
}

//    Saves options clase 0 || clase 13: Project || DAMSGHOLE
void MainWindow::saveOptions_13(string file, bool *printwarns, QString *warns)
{

    CIniFile::DeleteSection("DAMSGHOLESECT",file);

    if (iswindows) write_option("iswindows", "DAMSGHOLESECT", QString("T"), file, printwarns, warns);
    else write_option("iswindows", "DAMSGHOLESECT", QString("F"), file, printwarns, warns);

    if (lvalence){
        write_option("lvalence", "DAMSGHOLESECT", QString("T"), file, printwarns, warns);
    }

    write_option("filename", "DAMSGHOLESECT", TXTSGholefilename->text().remove("\"").prepend("\"").append("\""),
                 file, printwarns, warns);
    write_option("gridname", "DAMSGHOLESECT", TXTImportSGholeden->text().remove("\"").prepend("\"").append("\""),
                 file, printwarns, warns);
    write_option("contourval", "DAMSGHOLESECT", TXTSGholecontour->text(), file, printwarns, warns);
//    write_option("topcolor", "DAMSGHOLESECT", TXTSGholecolor->text(), file, printwarns, warns);
    write_option("geomthr", "DAMSGHOLESECT",  QString("1.d"+SPBSGholegeomthreshold->text()), file, printwarns, warns);
    write_option("thrslocal", "DAMSGHOLESECT",  QString(SPBSGholelocalextrema->text()+".d-2"), file, printwarns, warns);
    write_option("umbrlargo", "DAMSGHOLESECT",  QString("1.d"+SPBSGholelongthreshold->text()), file, printwarns, warns);
    write_option("lmaxrep", "DAMSGHOLESECT", SPBSGholelmaxexp->text(), file, printwarns, warns);
    if (CHKSGhexactMESP->isChecked()) write_option("lexact", "DAMSGHOLESECT", QString("T"), file, printwarns, warns);
    else write_option("lexact", "DAMSGHOLESECT", QString("F"), file, printwarns, warns);
}

/*******************************************************************************************************/
/******************************  PROGRAMS FOR DAM ANALYSIS *********************************************/
/*******************************************************************************************************/


/***********************************************************************/
/*  page_atdens: ATOMIC DENSITIES                                      */
/***********************************************************************/


void MainWindow::CHKatdensinput_changed(int state)     
{
    if (state == 0){
        CHKatdensmpi->setEnabled(true);
        FRMatdensmpi->setVisible(true);
        if (CHKatdensmpi->isChecked()){
            LBLatdensmpi->setEnabled(true);
            SPBatdensmpi->setEnabled(true);
        }
        else{
            LBLatdensmpi->setEnabled(false);
            SPBatdensmpi->setEnabled(false);
        }
    }
    else{
        FRMatdensmpi->setVisible(false);
        CHKatdensmpi->setEnabled(false);
        LBLatdensmpi->setEnabled(false);
        SPBatdensmpi->setEnabled(false);
    }
}

void MainWindow::CHKatdensmpi_changed(int state)     
{
    if (state != 0 && !CHKatdensinput->isChecked()){
        CHKatdensmpi->setChecked(true);
        LBLatdensmpi->setEnabled(true);
        SPBatdensmpi->setEnabled(true);
    }
    else{
        if (state == 0) 
            CHKatdensmpi->setChecked(false);
        else{
            CHKatdensmpi->setChecked(true);
            CHKatdensmpi->setEnabled(false);
        }
        LBLatdensmpi->setEnabled(false);
        SPBatdensmpi->setEnabled(false);
    }
}

//    Executes external program DAM (Partition of molecular density into atomic densities)
void MainWindow::execDam(){
    if (lslater){
        QString sxyzfilename = ProjectFolder+ProjectName+".sxyz";
        if (!(QFile::exists(sxyzfilename))){
            execsgbs2sxyz(sxyzfilename);
            if (!(QFile::exists(sxyzfilename))){
                QMessageBox::warning(this, tr("DAMQT"),tr("File %1 does not exist").arg(sxyzfilename));
                return;
            }
        }
        set_natom(read_natom(sxyzfilename));
        QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
        existsinp(DirNombreArchivo,1,1,false);
        inputdatafile("DAMSTO320.inp","DAMSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAM from .damproj file
        QString fileName=ProjectFolder+ProjectName+"-DAMSTO320.inp";
        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
            return;
        }
        if (CHKatdensinput->isChecked()){
            QTextStream in(&file);
            textEdit->setFont(QFont("Courier",10));
            textEdit->setPlainText(in.readAll());
            QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMSTO320.inp"));
            return;
        }
        file.close();
        QString stdinput = ProjectFolder + ProjectName + "-DAMSTO320.inp";
        QString stdoutput;
        QString strprocess;
        if (CHKatdensmpi->isChecked()){
            QString processname = "DAMSTO320_mpi.exe";
            QString execName = get_execName(processname, QString("DAM320_mpi"));
            if (execName.isEmpty())
                return;
            strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                        .arg(SPBatdensmpi->value()).arg(TXTmpiflags->text()).arg(execName);
            stdoutput = ProjectFolder + ProjectName + "-DAMSTO320_mpi.out";
        }
        else{
            QString processname = "DAMSTO320.exe";
            QString execName = get_execName(processname, QString("DAM320"));
            if (execName.isEmpty())
                return;
            strprocess = QString(execName);
            stdoutput = ProjectFolder + ProjectName + "-DAMSTO320.out";
        }
        myProcess = new QProcess(this);
        myProcess->setStandardInputFile(stdinput);
        myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
        myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
        connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
        connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
        connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
        executing = 11;
        myProcess->start(strprocess);
    }
    else{
        QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
        existsinp(DirNombreArchivo,1,1,false);
        inputdatafile("DAMGTO320.inp","G-DAMSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for GDAM from .damproj file
        QString fileName=ProjectFolder+ProjectName+"-DAMGTO320.inp";
        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
            return;
        }
        if (CHKatdensinput->isChecked()){
            QTextStream in(&file);
            textEdit->setFont(QFont("Courier",10));
            textEdit->setPlainText(in.readAll());
            QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMGTO320.inp"));
            return;
        }
        file.close();
        QString stdinput = ProjectFolder + ProjectName + "-DAMGTO320.inp";
        QString stdoutput;
        QString strprocess;

        if (CHKatdensmpi->isChecked()){
            QString processname = "DAMGTO320_mpi.exe";
            QString execName = get_execName(processname, QString("DAM320_mpi"));
            if (execName.isEmpty())
                return;
            strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                        .arg(SPBatdensmpi->value()).arg(TXTmpiflags->text()).arg(execName);
            stdoutput = ProjectFolder + ProjectName + "-DAMGTO320_mpi.out";
        }
        else{
            QString processname = "DAMGTO320.exe";
            QString execName = get_execName(processname, QString("DAM320"));
            if (execName.isEmpty())
                return;
            strprocess = QString(execName);
            stdoutput = ProjectFolder + ProjectName + "-DAMGTO320.out";
        }
        myProcess = new QProcess(this);
        myProcess->setStandardInputFile(stdinput);
        myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
        myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
        connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
        connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
        connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
        executing = 11;
        myProcess->start(strprocess);
    }
    defineRanges();
}

void MainWindow::SPBatdensmpi_changed(int nprocessors)
{
    SPBatdensmpi->setValue(nprocessors);
}

void MainWindow::SPBatdenslmaxexp_changed()
{
    int valor=SPBatdenslmaxexp->value();
    int valor1=SPBdenslmaxexp->value();
    int valor2=SPBfradltab->value();
    int valor3=SPBpotlmaxexp->value();
    int valor4=SPBeflmaxexp->value();

    if (valor < valor1) SPBdenslmaxexp->setValue(valor);
    SPBdenslmaxexp->setRange(0,valor);
    if (valor < valor2) SPBfradltab->setValue(valor);
    SPBfradltab->setRange(0,valor);
    if (valor < valor3) SPBpotlmaxexp->setValue(valor);
    SPBpotlmaxexp->setRange(0,valor);
    if (valor < valor4) SPBeflmaxexp->setValue(valor);
    SPBeflmaxexp->setRange(0,valor);

    SPBdenslminexp->setRange(0,valor);
    SPBfradltab->setRange(0,valor);

    SPBmrotlmin->setRange(0,SPBatdenslmaxexp->value());
    SPBmrotlmax->setRange(0,SPBatdenslmaxexp->value());
    
    SPBatdenslmaxdisp->setRange(0, valor);
    if (valor < SPBatdenslmaxdisp->value()) SPBatdenslmaxdisp->setValue(valor);
}


/***************************************************************************/
/*  page_densgrad: DENSITY GRADIENT                                        */
/***************************************************************************/


void MainWindow::BTNdgfilelines_clicked()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
                    tr("lin files")+" (*.lin2D *.lin3D);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTdgfilelines->setText(fileName);
}


void MainWindow::CHKdgextralines_changed()
{
    if (!CHKdgextralines->isChecked()){
        LBLdgfilelines->setEnabled(false);
        LBLdgfilelines->setVisible(false);
        TXTdgfilelines->setEnabled(false);
        TXTdgfilelines->setVisible(false);
        BTNdgfilelines->setEnabled(false);
        BTNdgfilelines->setVisible(false);
        CHKdgxyz->setVisible(false);
        CHKdgxyz->setEnabled(false);
    }else{
        LBLdgfilelines->setEnabled(true);
        LBLdgfilelines->setVisible(true);
        TXTdgfilelines->setEnabled(true);
        TXTdgfilelines->setVisible(true);
        BTNdgfilelines->setEnabled(true);
        BTNdgfilelines->setVisible(true);
        CHKdgxyz->setVisible(true);
        CHKdgxyz->setEnabled(true);
    }
}

void MainWindow::CHKdginput_changed(int state)
{
    if (state == 0 && RBTdg3D->isChecked()){
        CHKdgmpi->setEnabled(true);
        FRMdgmpi->setVisible(true);
        if (CHKdgmpi->isChecked()){
            LBLdgmpi->setEnabled(true);
            SPBdgmpi->setEnabled(true);
        }
        else{
            LBLdgmpi->setEnabled(false);
            SPBdgmpi->setEnabled(false);
        }
    }
    else{
        FRMdgmpi->setVisible(false);
        CHKdgmpi->setEnabled(false);
        LBLdgmpi->setEnabled(false);
        SPBdgmpi->setEnabled(false);
    }
}

void MainWindow::CHKdgmpi_changed(int state)
{
    if (state != 0 && RBTdg3D->isChecked() && !CHKdginput->isChecked()){
        CHKdgmpi->setChecked(true);
        LBLdgmpi->setEnabled(true);
        SPBdgmpi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKdgmpi->setChecked(false);
        else{
            CHKdgmpi->setChecked(true);
            CHKdgmpi->setEnabled(false);
        }
        LBLdgmpi->setEnabled(false);
        SPBdgmpi->setEnabled(false);
    }
}

void MainWindow::CHKdgxyz_changed()
{
    if(CHKdgxyz->isChecked()){
        if (RBTdg3D->isChecked()){
            Wtabledendg->setEnabled(true);
            Wtabledendg->setVisible(true);
            Wtable2dg->setEnabled(false);
            Wtable2dg->setVisible(false);
        }
        else{
            Wtabledendg->setEnabled(false);
            Wtabledendg->setVisible(false);
            Wtable2dg->setEnabled(true);
            Wtable2dg->setVisible(true);
        }
    }else{
        Wtabledendg->setEnabled(false);
        Wtabledendg->setVisible(false);
        Wtable2dg->setEnabled(false);
        Wtable2dg->setVisible(false);
    }
}

//    Executes external program DAMDENGRAD  (Computes density gradient lines from the atomic partition)
void MainWindow::execDamdengrad()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,12,1,false);
    QString fileNameDg = ProjectName;
    if (!TXTdgfilename->text().isEmpty())
            fileNameDg=TXTdgfilename->text();
    inputdatafile("DAMDENGRAD320.inp","DAMDENGRADSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMDENGRAD from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMDENGRAD320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKdginput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMDENGRAD320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMDENGRAD320.inp";
    QString stdoutput;
    QString strprocess;
    if (RBTdg3D->isChecked() && CHKdgmpi->isChecked()){
        QString processname = "DAMDENGRAD320_mpi.exe";
        QString execName = get_execName(processname, QString("DAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBdgmpi->value()).arg(TXTmpiflags->text()).arg(execName);
        stdoutput = ProjectFolder + fileNameDg + "-DAMDENGRAD320_mpi.out";
    }
    else{
        QString processname = "DAMDENGRAD320.exe";
        QString execName = get_execName(processname, QString("DAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameDg + "-DAMDENGRAD320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 22;
    myProcess->start(strprocess);
}


void MainWindow::RBTdg2D3D_changed(){
    if(RBTdg2D->isChecked()){
        FRMdgplot2D->setVisible(true);
        TXTdguinf->setEnabled(true);
        TXTdgusup->setEnabled(true);
        TXTdgvinf->setEnabled(true);
        TXTdgvsup->setEnabled(true);
        FRMdgplane2D->setVisible(true);
        FRMdgplot3D->setVisible(false);
        TXTdgxinf->setEnabled(false);
        TXTdgxsup->setEnabled(false);
        TXTdgyinf->setEnabled(false);
        TXTdgysup->setEnabled(false);
        TXTdgzinf->setEnabled(false);
        TXTdgzsup->setEnabled(false);
        if (mpi){
            FRMdgmpi->setVisible(false);
            CHKdgmpi->setEnabled(false);
            LBLdgmpi->setEnabled(false);
            SPBdgmpi->setEnabled(false);
        }
    }
    else{
        FRMdgplot2D->setVisible(false);
        TXTdguinf->setEnabled(false);
        TXTdgusup->setEnabled(false);
        TXTdgvinf->setEnabled(false);
        TXTdgvsup->setEnabled(false);
        FRMdgplane2D->setVisible(false);
        FRMdgplot3D->setVisible(true);
        TXTdgxinf->setEnabled(true);
        TXTdgxsup->setEnabled(true);
        TXTdgyinf->setEnabled(true);
        TXTdgysup->setEnabled(true);
        TXTdgzinf->setEnabled(true);
        TXTdgzsup->setEnabled(true);
        if (mpi){
            if (!CHKdginput->isChecked()){
                FRMdgmpi->setVisible(true);
                CHKdgmpi->setEnabled(true);
                if (CHKdgmpi->isChecked()){
                    LBLdgmpi->setEnabled(true);
                    SPBdgmpi->setEnabled(true);
                }
                else{
                    LBLdgmpi->setEnabled(false);
                    SPBdgmpi->setEnabled(false);
                }
            }
            else{
                FRMdgmpi->setVisible(false);
                CHKdgmpi->setEnabled(false);
                LBLdgmpi->setEnabled(false);
                SPBdgmpi->setEnabled(false);
            }
        }
    }
    CHKdgextralines->setChecked(false);
    CHKdgxyz->setChecked(false);
    CHKdgxyz->setEnabled(false);
    CHKdgxyz->setVisible(false);
}

void MainWindow::RBTdg2Dplanes_changed(){
    if (RBTdgplaneXY->isChecked()){
        TXTdgplaneA->setText("0.");
        TXTdgplaneB->setText("0.");
        TXTdgplaneC->setText("1.");
        FRMdgplaneABC->setVisible(false);
    }
    if (RBTdgplaneXZ->isChecked()){
        TXTdgplaneA->setText("0.");
        TXTdgplaneB->setText("1.");
        TXTdgplaneC->setText("0.");
        FRMdgplaneABC->setVisible(false);
    }
    if (RBTdgplaneYZ->isChecked()){
        TXTdgplaneA->setText("1.");
        TXTdgplaneB->setText("0.");
        TXTdgplaneC->setText("0.");
        FRMdgplaneABC->setVisible(false);
    }
    if (RBTdgplaneABC->isChecked()){
        FRMdgplaneABC->setVisible(true);
    }
}

void MainWindow::SPBdgmpi_changed(int nprocessors)
{
    SPBdgmpi->setValue(nprocessors);
}

void MainWindow::TXTdgdlt0_changed()
{
    double dnum=TXTdgdlt0->text().toDouble();
    if (dnum < 0.0) TXTdgdlt0->setText("0.001");
    if (dnum > 0.1) TXTdgdlt0->setText("0.1");
}


/***************************************************************************/
/*  page_Efield: ELECTRIC FIELD                                            */
/***************************************************************************/


void MainWindow::BTNeffilelines_clicked()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
                    tr("lin files")+" (*.lin2D *.lin3D);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTeffilelines->setText(fileName);
}

void MainWindow::CHKefinput_changed(int state)
{
    if (state == 0 && RBTef3D->isChecked()){
        FRMefmpi->setVisible(true);
        CHKefmpi->setEnabled(true);
        if (CHKefmpi->isChecked()){
            LBLefmpi->setEnabled(true);
            SPBefmpi->setEnabled(true);
        }
        else{
            LBLefmpi->setEnabled(false);
            SPBefmpi->setEnabled(false);
        }
    }
    else{
        FRMefmpi->setVisible(false);
        CHKefmpi->setEnabled(false);
        LBLefmpi->setEnabled(false);
        SPBefmpi->setEnabled(false);
    }
}

void MainWindow::CHKefmpi_changed(int state)
{
    if (state != 0 && RBTef3D->isChecked() && !CHKefinput->isChecked()){
        CHKefmpi->setChecked(true);
        LBLefmpi->setEnabled(true);
        SPBefmpi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKefmpi->setChecked(false);
        else{
            CHKefmpi->setChecked(true);
            CHKefmpi->setEnabled(false);
        }
        LBLefmpi->setEnabled(false);
        SPBefmpi->setEnabled(false);
    }
}

void MainWindow::CHKefextralines_changed()
{
    if (!CHKefextralines->isChecked()){
        LBLeffilelines->setEnabled(false);
        LBLeffilelines->setVisible(false);
        TXTeffilelines->setEnabled(false);
        TXTeffilelines->setVisible(false);
        BTNeffilelines->setEnabled(false);
        BTNeffilelines->setVisible(false);
        CHKefxyz->setVisible(false);
        CHKefxyz->setEnabled(false);
    }else{
        LBLeffilelines->setEnabled(true);
        LBLeffilelines->setVisible(true);
        TXTeffilelines->setEnabled(true);
        TXTeffilelines->setVisible(true);
        BTNeffilelines->setEnabled(true);
        BTNeffilelines->setVisible(true);
        CHKefxyz->setVisible(true);
        CHKefxyz->setEnabled(true);
    }
}

void MainWindow::CHKeflong_changed()
{
    if (CHKeflong->isChecked()){
        LBLeflongthreshold->setEnabled(false);
        SPBeflongthreshold->setEnabled(false);
    }else{
        LBLeflongthreshold->setEnabled(true);
        SPBeflongthreshold->setEnabled(true);
    }
}

void MainWindow::CHKefxyz_changed()
{
    if(CHKefxyz->isChecked()){
        if (RBTef3D->isChecked()){
            Wtabledenef->setEnabled(true);
            Wtabledenef->setVisible(true);
            Wtable2ef->setEnabled(false);
            Wtable2ef->setVisible(false);
        }
        else{
            Wtabledenef->setEnabled(false);
            Wtabledenef->setVisible(false);
            Wtable2ef->setEnabled(true);
            Wtable2ef->setVisible(true);
        }
    }else{
        Wtabledenef->setEnabled(false);
        Wtabledenef->setVisible(false);
        Wtable2ef->setEnabled(false);
        Wtable2ef->setVisible(false);
    }
}

//    Executes external program DAMFIELD  (Computes electric field lines from the atomic partition)
void MainWindow::execDamfield()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,5,1,false);
    QString fileNameEf = ProjectName;
    if (!TXTeffilename->text().isEmpty())
            fileNameEf=TXTeffilename->text();
    inputdatafile("DAMFIELD320.inp","DAMFIELDSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMFIELD from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMFIELD320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKefinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMFIELD320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMFIELD320.inp";
    QString stdoutput;
    QString strprocess;
    if (RBTef3D->isChecked() && CHKefmpi->isChecked()){
        QString processname = "DAMFIELD320_mpi.exe";
        QString execName = get_execName(processname, QString("DAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBefmpi->value()).arg(TXTmpiflags->text()).arg(execName);
        stdoutput = ProjectFolder + fileNameEf + "-DAMFIELD320_mpi.out";
    }
    else{
        QString processname = "DAMFIELD320.exe";
        QString execName = get_execName(processname, QString("DAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameEf + "-DAMFIELD320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 15;
    myProcess->start(strprocess);
}


void MainWindow::RBTef2D3D_changed(){
    if(RBTef2D->isChecked()){
        FRMefplot2D->setVisible(true);
        TXTefuinf->setEnabled(true);
        TXTefusup->setEnabled(true);
        TXTefvinf->setEnabled(true);
        TXTefvsup->setEnabled(true);
        FRMefplane2D->setVisible(true);
        FRMefplot3D->setVisible(false);
        TXTefxinf->setEnabled(false);
        TXTefxsup->setEnabled(false);
        TXTefyinf->setEnabled(false);
        TXTefysup->setEnabled(false);
        TXTefzinf->setEnabled(false);
        TXTefzsup->setEnabled(false);
        if (mpi){
            FRMefmpi->setVisible(false);
            CHKefmpi->setEnabled(false);
            LBLefmpi->setEnabled(false);
            SPBefmpi->setEnabled(false);
        }
    }
    else{
        FRMefplot2D->setVisible(false);
        TXTefuinf->setEnabled(false);
        TXTefusup->setEnabled(false);
        TXTefvinf->setEnabled(false);
        TXTefvsup->setEnabled(false);
        FRMefplane2D->setVisible(false);
        FRMefplot3D->setVisible(true);
        TXTefxinf->setEnabled(true);
        TXTefxsup->setEnabled(true);
        TXTefyinf->setEnabled(true);
        TXTefysup->setEnabled(true);
        TXTefzinf->setEnabled(true);
        TXTefzsup->setEnabled(true);
        if (mpi){
            if (!CHKefinput->isChecked()){
                FRMefmpi->setVisible(true);
                CHKefmpi->setEnabled(true);
                if (CHKefmpi->isChecked()){
                    LBLefmpi->setEnabled(true);
                    SPBefmpi->setEnabled(true);
                }
                else{
                    LBLefmpi->setEnabled(false);
                    SPBefmpi->setEnabled(false);
                }
            }
            else{
                FRMefmpi->setVisible(false);
                CHKefmpi->setEnabled(false);
                LBLefmpi->setEnabled(false);
                SPBefmpi->setEnabled(false);
            }
        }
    }
    CHKefextralines->setChecked(false);
    CHKefxyz->setChecked(false);
    CHKefxyz->setEnabled(false);
    CHKefxyz->setVisible(false);
}

void MainWindow::RBTef2Dplanes_changed(){
    if (RBTefplaneXY->isChecked()){
        TXTefplaneA->setText("0.");
        TXTefplaneB->setText("0.");
        TXTefplaneC->setText("1.");
        FRMefplaneABC->setVisible(false);
    }
    if (RBTefplaneXZ->isChecked()){
        TXTefplaneA->setText("0.");
        TXTefplaneB->setText("1.");
        TXTefplaneC->setText("0.");
        FRMefplaneABC->setVisible(false);
    }
    if (RBTefplaneYZ->isChecked()){
        TXTefplaneA->setText("1.");
        TXTefplaneB->setText("0.");
        TXTefplaneC->setText("0.");
        FRMefplaneABC->setVisible(false);
    }
    if (RBTefplaneABC->isChecked()){
        FRMefplaneABC->setVisible(true);
    }
}

void MainWindow::SPBefmpi_changed(int nprocessors)
{
    SPBefmpi->setValue(nprocessors);
}

void MainWindow::TXTefdlt0_changed()
{
    double dnum=TXTefdlt0->text().toDouble();
    if (dnum < 0.0) TXTefdlt0->setText("0.001");
    if (dnum > 0.1) TXTefdlt0->setText("0.1");
}


/****************************************************************************/
/*  page_frad: RADIAL FACTORS                                               */
/****************************************************************************/


void MainWindow::CHKfradextras_changed()
{
    if (CHKfradextras->QAbstractButton::isChecked()) {
        Wtable5->setVisible(true);
        SHTfradrlist->setVisible(true);
        SHTfradrlist->BTNadd_clicked();
    } else {
        Wtable5->setVisible(false);
        SHTfradrlist->setVisible(false);
        SHTfradrlist->clear();
    }
}

//    Executes external program DAMFRAD  (Computes radial factors of the atomic partition)
void MainWindow::execDamfrad()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,6,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTfraddamfilename->text().isEmpty())
            fileNameOut=TXTfraddamfilename->text();
    inputdatafile("DAMFRAD320.inp","DAMFRADSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMFRAD from .damproj file
    QString strprocess;
    QString processname = "DAMFRAD320.exe";
    QString execName = get_execName(processname, QString("DAM320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    QString stdinput = ProjectFolder + ProjectName + "-DAMFRAD320.inp";
    QString stdoutput = ProjectFolder + fileNameOut + "-DAMFRAD320.out";
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 16;    
    myProcess->start(strprocess);
}

void MainWindow::SPBfradltab_changed()
{
    SPBfradmtab->setRange(-SPBfradltab->value(),SPBfradltab->value());
}

void MainWindow::SPBfradmtab_changed()
{
    if (abs(SPBfradmtab->value())>SPBfradltab->value()){
        if(SPBfradmtab->value()>0)
            SPBfradmtab->setValue(SPBfradmtab->value()-1);
        else
            SPBfradmtab->setValue(SPBfradmtab->value()+1);
    }
}

void MainWindow::TXTfradatoms_changed(){
    QString string = QString(TXTfradatoms->text());
    QStringList stringlist1 = string.split(",");
    fradlist->clear();
    int ncen = get_natom();
    for (int i = 0 ; i < stringlist1.length() ; ++i){
        if (stringlist1.at(i).length() != 0){
            if (stringlist1.at(i).contains(QString("-"))){
                QStringList stringlist2 = stringlist1.at(i).split("-");
                if (stringlist2.at(0).length() != 0 && stringlist2.at(1).length() != 0 
                        && stringlist2.at(1).toInt() > stringlist2.at(0).toInt()){
                    for (int k = stringlist2.at(0).toInt() ; k <= stringlist2.at(1).toInt() ; k++){
                        if (k > 0 && k <= ncen) fradlist->append(QString("%1").arg(k));
                    }
                }
            }
            else{
                   int k = stringlist1.at(i).toInt();
                   if (k > 0 && k <= ncen) fradlist->append(stringlist1.at(i));
            }
        }
    }
    fradlist->sort();
    fradlist->removeDuplicates();
}


/***************************************************************************/
/*  page_HFforces: HELLMANN-FEYNMAN FORCES                                 */
/***************************************************************************/


void MainWindow::CHKHFlatomsel_changed()
{
    if (CHKHFlatomsel->QAbstractButton::isChecked()) {
        TXTHFforcesatoms->setEnabled(true);
    } else {
        TXTHFforcesatoms->setEnabled(false);
    }
}
//    Executes external program DAMFORCES  (Computes Hellmann-Feynman forces on nuclei from the atomic partition)
void MainWindow::execDamforces()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,4,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTHFgdamforcesfilename->text().isEmpty())
            fileNameOut=TXTHFgdamforcesfilename->text();
    inputdatafile("DAMFORCES320.inp","DAMFORCESSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMFORCES from .damproj file
    QString strprocess;
    QString processname = "DAMFORCES320.exe";
    QString execName = get_execName(processname, QString("DAM320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    QString stdinput = ProjectFolder + ProjectName + "-DAMFORCES320.inp";
    QString stdoutput = ProjectFolder + fileNameOut + "-DAMFORCES320.out";
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 14;
    myProcess->start(strprocess);
}

void MainWindow::TXTHFforcesatoms_changed(){
    QString string = QString(TXTHFforcesatoms->text());
    QStringList stringlist1 = string.split(",");
    HFforceslist->clear();
    int ncen = get_natom();
    for (int i = 0 ; i < stringlist1.length() ; ++i){
        if (stringlist1.at(i).length() != 0){
            if (stringlist1.at(i).contains(QString("-"))){
                QStringList stringlist2 = stringlist1.at(i).split("-");
                if (stringlist2.at(0).length() != 0 && stringlist2.at(1).length() != 0
                        && stringlist2.at(1).toInt() > stringlist2.at(0).toInt()){
                    for (int k = stringlist2.at(0).toInt() ; k <= stringlist2.at(1).toInt() ; k++){
                        if (k > 0 && k <= ncen) HFforceslist->append(QString("%1").arg(k));
                    }
                }
            }
            else{
                int k = stringlist1.at(i).toInt();
                if (k > 0 && k <= ncen) HFforceslist->append(stringlist1.at(i));
            }
        }
    }
    HFforceslist->sort();
    HFforceslist->removeDuplicates();
}


/***************************************************************************/
/*  page_MED: DENSITY                                                      */
/***************************************************************************/


void MainWindow::CHKdensder2_changed(int state)     // If second derivatives are checked forces gradient to be checked
{
    if (state == Qt::Checked) {
        CHKdensgrad->setChecked(true);
    }
}

void MainWindow::CHKdensgrad_changed(int state)     // If second derivatives option is checked unckecks it
{
    if (CHKdensder2->QAbstractButton::isChecked() && state == Qt::Unchecked) {
        CHKdensgrad->setChecked(true);
    }
}

void MainWindow::CHKdensinput_changed(int state)
{
    if (state == 0 && RBTdens3D->isChecked()){
        FRMdensmpi->setVisible(true);
        CHKdensmpi->setEnabled(true);
        if (CHKdensmpi->isChecked()){
            LBLdensmpi->setEnabled(true);
            SPBdensmpi->setEnabled(true);
        }
        else{
            LBLdensmpi->setEnabled(false);
            SPBdensmpi->setEnabled(false);
        }
    }
    else{
        FRMdensmpi->setVisible(false);
        CHKdensmpi->setEnabled(false);
        LBLdensmpi->setEnabled(false);
        SPBdensmpi->setEnabled(false);
    }
}

void MainWindow::CHKdensmpi_changed(int state)
{
    if (state != 0 && RBTdens3D->isChecked() && !CHKdensinput->isChecked()){
        CHKdensmpi->setChecked(true);
        LBLdensmpi->setEnabled(true);
        SPBdensmpi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKdensmpi->setChecked(false);
        else{
            CHKdensmpi->setChecked(true);
            CHKdensmpi->setEnabled(false);
        }
        LBLdensmpi->setEnabled(false);
        SPBdensmpi->setEnabled(false);
    }
}

void MainWindow::CHKfragments_changed(int)
{
    if ((CHKdenslatomics->isChecked() || CHKdensldensacc->isChecked())) {
        FRMdensatoms->setVisible(true);
    } else {
        FRMdensatoms->setVisible(false);
    }
}

void MainWindow::CHKdensgrid_changed(){
    if(CHKdensgrid->isChecked()){
        if(RBTdens2D->isChecked()){
            FRMdensgrid2D->setVisible(true);
            TXTdensuinf->setEnabled(true);
            TXTdensusup->setEnabled(true);
            TXTdensvinf->setEnabled(true);
            TXTdensvsup->setEnabled(true);
            FRMdensgrid3D->setVisible(false);
            TXTdensxinf->setEnabled(false);
            TXTdensxsup->setEnabled(false);
            TXTdensyinf->setEnabled(false);
            TXTdensysup->setEnabled(false);
            TXTdenszinf->setEnabled(false);
            TXTdenszsup->setEnabled(false);
            if (RBTdensrcustom->isChecked()){
                FRMdensresol2D->setHidden(false);
            }
            FRMdensresol3D->setHidden(true);
            if (mpi)
                FRMdensmpi->setVisible(false);
        }
        else{
            FRMdensgrid2D->setVisible(false);
            TXTdensuinf->setEnabled(false);
            TXTdensusup->setEnabled(false);
            TXTdensvinf->setEnabled(false);
            TXTdensvsup->setEnabled(false);
            FRMdensgrid3D->setVisible(true);
            TXTdensxinf->setEnabled(true);
            TXTdensxsup->setEnabled(true);
            TXTdensyinf->setEnabled(true);
            TXTdensysup->setEnabled(true);
            TXTdenszinf->setEnabled(true);
            TXTdenszsup->setEnabled(true);
            if (RBTdensrcustom->isChecked()){
                FRMdensresol3D->setHidden(false);
            }
            FRMdensresol2D->setHidden(true);
            if (mpi){
                if( !CHKdensinput->isChecked())
                    FRMdensmpi->setVisible(true);
                else
                    FRMdensmpi->setVisible(false);
            }
        }
        FRMdensgridtype->setVisible(true);
        FRMdensgridres->setVisible(true);
        RBTdens2D->setEnabled(true);
        RBTdens3D->setEnabled(true);
        RBTdensrcustom->setEnabled(true);
        RBTdensrlow->setEnabled(true);
        RBTdensrmedium->setEnabled(true);
        RBTdensrhigh->setEnabled(true);
    }
    else{
        if (mpi)
            FRMdensmpi->setVisible(false);
        FRMdensgridtype->setVisible(false);
        FRMdensgridres->setVisible(false);
        FRMdensgrid2D->setVisible(false);
        TXTdensuinf->setEnabled(false);
        TXTdensusup->setEnabled(false);
        TXTdensvinf->setEnabled(false);
        TXTdensvsup->setEnabled(false);
        FRMdensgrid3D->setVisible(false);
        TXTdensxinf->setEnabled(false);
        TXTdensxsup->setEnabled(false);
        TXTdensyinf->setEnabled(false);
        TXTdensysup->setEnabled(false);
        TXTdenszinf->setEnabled(false);
        TXTdenszsup->setEnabled(false);
        RBTdens2D->setEnabled(false);
        RBTdens3D->setEnabled(false);
        RBTdensrcustom->setEnabled(false);
        RBTdensrlow->setEnabled(false);
        RBTdensrmedium->setEnabled(false);
        RBTdensrhigh->setEnabled(false);
    }
    if (RBTdens3D->isChecked()){
        RBTdensrlow->setToolTip("65x65x65");
        RBTdensrmedium->setToolTip("129x129x129");
        RBTdensrhigh->setToolTip("257x257x257");
    }
    else{
        RBTdensrlow->setToolTip("129x129");
        RBTdensrmedium->setToolTip("257x257");
        RBTdensrhigh->setToolTip("513x513");
    }
}

void MainWindow::CHKdensxyz_changed()
{
    if(CHKdensxyz->isChecked()){
        Wtableden->setEnabled(true);
        Wtableden->setVisible(true);
    }else{
        Wtableden->setEnabled(false);
        Wtableden->setVisible(false);
    }
}

void MainWindow::density_resolution_changed(){
    if (RBTdensrcustom->isChecked())
    if (RBTdens3D->isChecked()){
        FRMdensresol3D->setHidden(false);
        FRMdensresol2D->setHidden(true);
    }
    else{
        FRMdensresol2D->setHidden(false);
        FRMdensresol3D->setHidden(true);
    }
    else{
    FRMdensresol2D->setHidden(true);
    FRMdensresol3D->setHidden(true);
    }
}

//    Executes external program DAMDEN  (Computes molecular density or deformations from the atomic partition)
void MainWindow::execDamden()
{

    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,2,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTdensdamdenfilename->text().isEmpty())
            fileNameOut=TXTdensdamdenfilename->text();

    inputdatafile("DAMDEN320.inp","DAMDENSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMDEN from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMDEN320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKdensinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMDEN320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMDEN320.inp";
    QString stdoutput;
    QString strprocess;

    if (RBTdens3D->isChecked() && CHKdensmpi->isChecked() && CHKdensgrid->isChecked()){
        QString processname = "DAMDEN320_mpi.exe";
        QString execName = get_execName(processname, QString("DAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                .arg(SPBdensmpi->value()).arg(TXTmpiflags->text()).arg(execName);
        stdoutput = ProjectFolder + fileNameOut;
        if (RBTdensExact->isChecked()){
            stdoutput += "_exact";
        }
        stdoutput += "-DAMDEN320_mpi.out";
    }
    else{
        QString processname = "DAMDEN320.exe";
        QString execName = get_execName(processname, QString("DAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameOut;
        if (RBTdensExact->isChecked()){
            stdoutput += "_exact";
        }
        stdoutput += "-DAMDEN320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 12;

    myProcess->start(strprocess);

    if (RBTdens2D->isChecked() && RBTdensplane->isChecked()){
        QString fileinpstr=ProjectFolder+ProjectName+".xyz";
        QFile fileinp(fileinpstr);
        if (fileinp.open(QFile::ReadOnly | QFile::Text)) {
            densplanecase = get_plane_case(TXTdensplaneA->text().toDouble(), TXTdensplaneB->text().toDouble(), TXTdensplaneC->text().toDouble());
            if (densplanecase > 0 && densplanecase < 8){
                QVector<QString> plane;
                plane << "XY0" << "X0Z" << "0YZ" << "AB0" << "A0C" << "0BC" << "ABC";
                QString fileoutstr = ProjectFolder+fileNameOut;
                if (RBTdensExact->isChecked()){
                    fileoutstr += "_exact";
                }
                fileoutstr = fileoutstr+"_"+plane.at(densplanecase-1)+"-d.plane";
                QFile fileout(fileoutstr);
                if(fileout.open(QFile::Text | QFile::WriteOnly)){
                    QTextStream in(&fileinp); // Buffer for reading from fileinput
                    QTextStream out(&fileout); // Buffer for writing to fileout
                    QString line;
                    Elements *Elem = new Elements();
                    while (!in.atEnd()){
                        line = in.readLine();
#if QT_VERSION < 0x050E00
                        QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
                        QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
                        if (xyz.count() == 4){
                            double x = xyz[1].toDouble() * ANGSTROMTOBOHR;
                            double y = xyz[2].toDouble() * ANGSTROMTOBOHR;
                            double z = xyz[3].toDouble() * ANGSTROMTOBOHR;
                            if (qAbs(QVector3D::dotProduct(QVector3D(x,y,z),
                                QVector3D(TXTdensplaneA->text().toDouble(), TXTdensplaneB->text().toDouble(), TXTdensplaneC->text().toDouble()))) > 1.e-8)
                                continue;
                            double u = QVector3D::dotProduct(QVector3D(x,y,z),wu) / wu.length();
                            double v = QVector3D::dotProduct(QVector3D(x,y,z),wv) / wv.length();
                            int znuc = Elem->getZsymbol(xyz[0]);
                            out << QString("%1  %2  %3  %4  %5  %6\n").arg(u).arg(v).arg(znuc).arg(x).arg(y).arg(z);
                        }
                    }
                    out << QString("%1 %2 %3").arg(TXTdensplaneA->text().toDouble()).arg(TXTdensplaneB->text().toDouble()).arg(TXTdensplaneC->text().toDouble());
                    fileout.close();
                }
            }
        }
    }
}

void MainWindow::RBTdens2D3D_changed()
{
    if (RBTdens2D->isChecked()){
        FRMdensgrid2D->setVisible(true);
        FRMdensgrid3D->setVisible(false);
        CHKdensgrad->setChecked(false);
        CHKdensder2->setChecked(false);
        CHKdenslaplacian->setChecked(false);
        if (RBTdensrcustom->isChecked()){
            FRMdensresol2D->setHidden(false);
        }
        FRMdensresol3D->setHidden(true);
        if (mpi){
            FRMdensmpi->setVisible(false);
            CHKdensmpi->setEnabled(false);
            LBLdensmpi->setEnabled(false);
            SPBdensmpi->setEnabled(false);
        }
    }
    else{
        FRMdensgrid3D->setVisible(true);
        FRMdensgrid2D->setVisible(false);
        if (RBTdensrcustom->isChecked()){
            FRMdensresol3D->setHidden(false);
        }
        FRMdensresol2D->setHidden(true);
        if (mpi){
            if (!CHKdensinput->isChecked()){
                FRMdensmpi->setVisible(true);
                CHKdensmpi->setEnabled(true);
                if (CHKdensmpi->isChecked()){
                    LBLdensmpi->setEnabled(true);
                    SPBdensmpi->setEnabled(true);
                }
                else{
                    LBLdensmpi->setEnabled(false);
                    SPBdensmpi->setEnabled(false);
                }
            }
            else{
                FRMdensmpi->setVisible(false);
                CHKdensmpi->setEnabled(false);
                LBLdensmpi->setEnabled(false);
                SPBdensmpi->setEnabled(false);
            }
        }
        CHKdensgrad->setChecked(true);
    }
    if (RBTdens3D->isChecked()){
        RBTdensrlow->setToolTip("65x65x65");
        RBTdensrmedium->setToolTip("129x129x129");
        RBTdensrhigh->setToolTip("257x257x257");
    }
    else{
        RBTdensrlow->setToolTip("129x129");
        RBTdensrmedium->setToolTip("257x257");
        RBTdensrhigh->setToolTip("513x513");
    }
    CHKdensgrid_changed();
}

void MainWindow::RBTdenstype_changed(){
    if (RBTdensfulldensity->QAbstractButton::isChecked()) {
        SPBdenslminexp->setValue(0);
        SPBdenslminexp->setEnabled(false);
    }
    else if (RBTdensdeform->QAbstractButton::isChecked()){
        SPBdenslminexp->setValue(1);
        SPBdenslminexp->setEnabled(false);
    }
    else if (RBTdenslrange->QAbstractButton::isChecked()){
        SPBdenslminexp->setEnabled(true);
    }
}

void MainWindow::RBTdensplane_changed()
{
    if (RBTdensplane->QAbstractButton::isChecked()) {
        FRMdenssurfpar->setVisible(false);
        FRMdensplane2D->setVisible(true);
        densplanecase = get_plane_case(TXTdensplaneA->text().toDouble(),TXTdensplaneB->text().toDouble(),TXTdensplaneC->text().toDouble());
    }else{
        FRMdenssurfpar->setVisible(true);
        FRMdensplane2D->setVisible(false);
    }
}


void MainWindow::RBTdens2Dplanes_changed(){
    if (RBTdensplaneXY->isChecked()){
        TXTdensxformula2D->setText("u");
        TXTdensyformula2D->setText("v");
        TXTdenszformula2D->setText("0");
        TXTdensplaneA->setText("0.");
        TXTdensplaneB->setText("0.");
        TXTdensplaneC->setText("1.");
        FRMdensplaneABC->setVisible(false);
    }
    else if (RBTdensplaneXZ->isChecked()){
        TXTdensxformula2D->setText("u");
        TXTdensyformula2D->setText("0");
        TXTdenszformula2D->setText("v");
        TXTdensplaneA->setText("0.");
        TXTdensplaneB->setText("1.");
        TXTdensplaneC->setText("0.");
        FRMdensplaneABC->setVisible(false);
    }
    else if (RBTdensplaneYZ->isChecked()){
        TXTdensxformula2D->setText("0");
        TXTdensyformula2D->setText("u");
        TXTdenszformula2D->setText("v");
        TXTdensplaneA->setText("1.");
        TXTdensplaneB->setText("0.");
        TXTdensplaneC->setText("0.");
        FRMdensplaneABC->setVisible(false);
    }
    else if (RBTdensplaneABC->isChecked()){
        TXTdensxformula2D->setText("u");
        TXTdensyformula2D->setText("v");
        TXTdenszformula2D->setText("-(("+TXTdensplaneA->text()+"*u)+("+TXTdensplaneB->text()+"*v))/("+TXTdensplaneC->text()+")");
        FRMdensplaneABC->setVisible(true);
    }
    densplanecase = get_plane_case(TXTdensplaneA->text().toDouble(),TXTdensplaneB->text().toDouble(),TXTdensplaneC->text().toDouble());
}

void MainWindow::RBTdenslexact_changed()
{
    if (RBTdensExact->QAbstractButton::isChecked()) {
        FRMdensdeformation->setVisible(false);
        FRMdensderivs->setVisible(false);
        FRMdensfragments->setVisible(false);
        RBTdensfulldensity->setChecked(true);
    }else{
        FRMdensdeformation->setVisible(true);
        FRMdensderivs->setVisible(true);
        FRMdensfragments->setVisible(true);
    }
}

void MainWindow::rename_density_cntfile(){
    QString aux;
    if (RBTdensfulldensity->isChecked()){
        if (RBTdensExact){
            aux = "_exact";
        }
        else{
            aux = "";
        }
    }
    else{
        aux = "_deform";
    }
    QFile filecnt(ProjectFolder + ProjectName + aux + "-d.cnt");
    if (filecnt.exists() && planesuffix(densplanecase) != ""){
        QFile fileold(ProjectFolder + ProjectName + aux + planesuffix(densplanecase) + "-d.cnt");
        if (fileold.exists())
            fileold.remove();
        filecnt.rename(ProjectFolder + ProjectName + aux + planesuffix(densplanecase) + "-d.cnt");
    }
}

void MainWindow::SPBdensmpi_changed(int nprocessors)
{
    SPBdensmpi->setValue(nprocessors);
}

void MainWindow::SPBdenslmaxexp_changed()
{
    if (SPBdenslmaxexp->value()<SPBdenslminexp->value())
        SPBdenslmaxexp->setValue(SPBdenslmaxexp->value()+1);
}

void MainWindow::SPBdenslminexp_changed()
{
    if (SPBdenslminexp->value()>SPBdenslmaxexp->value())
        SPBdenslminexp->setValue(SPBdenslminexp->value()-1);
}

void MainWindow::TXTdensatoms_changed(){
    QString string = QString(TXTdensatoms->text());
    QStringList stringlist1 = string.split(",");
    denslist->clear();
    int ncen = get_natom();
    for (int i = 0 ; i < stringlist1.length() ; ++i){
        if (stringlist1.at(i).length() != 0){
            if (stringlist1.at(i).contains(QString("-"))){
                QStringList stringlist2 = stringlist1.at(i).split("-");
                if (stringlist2.at(0).length() != 0 && stringlist2.at(1).length() != 0
                        && stringlist2.at(1).toInt() > stringlist2.at(0).toInt()){
                    for (int k = stringlist2.at(0).toInt() ; k <= stringlist2.at(1).toInt() ; k++){
                        if (k > 0 && k <= ncen) denslist->append(QString("%1").arg(k));
                    }
                }
            }
            else{
                int k = stringlist1.at(i).toInt();
                if (k > 0 && k <= ncen) denslist->append(stringlist1.at(i));
            }
        }
    }
    denslist->sort();
    denslist->removeDuplicates();
}


/***************************************************************************/
/*  page_MESP: ELECTROSTATIC POTENTIAL                                     */
/***************************************************************************/

void MainWindow::CHKpotlong_changed()
{
    if (CHKpotlong->isChecked()){
        LBLpotlongthreshold->setEnabled(false);
        SPBpotlongthreshold->setEnabled(false);
    }else{
        LBLpotlongthreshold->setEnabled(true);
        SPBpotlongthreshold->setEnabled(true);
    }
}

void MainWindow::CHKpotgrid_changed(){
    if(CHKpotgrid->isChecked()){
        if(RBTpot2D->isChecked()){
            FRMpotgrid2D->setVisible(true);
            TXTpotuinf->setEnabled(true);
            TXTpotusup->setEnabled(true);
            TXTpotvinf->setEnabled(true);
            TXTpotvsup->setEnabled(true);
            FRMpotgrid3D->setVisible(false);
            TXTpotxinf->setEnabled(false);
            TXTpotxsup->setEnabled(false);
            TXTpotyinf->setEnabled(false);
            TXTpotysup->setEnabled(false);
            TXTpotzinf->setEnabled(false);
            TXTpotzsup->setEnabled(false);
            if (RBTpotrcustom->isChecked()){
                FRMpotresol2D->setHidden(false);
            }
            FRMpotresol3D->setHidden(true);
            if (mpi)
                FRMpotmpi->setVisible(false);
        }
        else{
            FRMpotgrid2D->setVisible(false);
            TXTpotuinf->setEnabled(false);
            TXTpotusup->setEnabled(false);
            TXTpotvinf->setEnabled(false);
            TXTpotvsup->setEnabled(false);
            FRMpotgrid3D->setVisible(true);
            TXTpotxinf->setEnabled(true);
            TXTpotxsup->setEnabled(true);
            TXTpotyinf->setEnabled(true);
            TXTpotysup->setEnabled(true);
            TXTpotzinf->setEnabled(true);
            TXTpotzsup->setEnabled(true);
            if (RBTpotrcustom->isChecked()){
                FRMpotresol3D->setHidden(false);
            }
            FRMpotresol2D->setHidden(true);
            if (mpi){
                if( !CHKpotinput->isChecked())
                    FRMpotmpi->setVisible(true);
                else
                    FRMpotmpi->setVisible(false);
            }
        }
        FRMpotgridtype->setVisible(true);
        FRMpotgridres->setVisible(true);
        RBTpot2D->setEnabled(true);
        RBTpot3D->setEnabled(true);
        RBTpotrcustom->setEnabled(true);
        RBTpotrlow->setEnabled(true);
        RBTpotrmedium->setEnabled(true);
        RBTpotrhigh->setEnabled(true);
        if (RBTpot3D->isChecked()){
            RBTpotrlow->setToolTip("65x65x65");
            RBTpotrmedium->setToolTip("129x129x129");
            RBTpotrhigh->setToolTip("257x257x257");
        }
        else{
            RBTpotrlow->setToolTip("129x129");
            RBTpotrmedium->setToolTip("257x257");
            RBTpotrhigh->setToolTip("513x513");
        }
    }
    else{
        if (mpi)
            FRMpotmpi->setVisible(false);
        FRMpotgridtype->setVisible(false);
        FRMpotgridres->setVisible(false);
        FRMpotgrid2D->setVisible(false);
        TXTpotuinf->setEnabled(false);
        TXTpotusup->setEnabled(false);
        TXTpotvinf->setEnabled(false);
        TXTpotvsup->setEnabled(false);
        FRMpotgrid3D->setVisible(false);
        TXTpotxinf->setEnabled(false);
        TXTpotxsup->setEnabled(false);
        TXTpotyinf->setEnabled(false);
        TXTpotysup->setEnabled(false);
        TXTpotzinf->setEnabled(false);
        TXTpotzsup->setEnabled(false);
        RBTpot2D->setEnabled(false);
        RBTpot3D->setEnabled(false);
        RBTpotrcustom->setEnabled(false);
        RBTpotrlow->setEnabled(false);
        RBTpotrmedium->setEnabled(false);
        RBTpotrhigh->setEnabled(false);
    }
}

void MainWindow::CHKpotder2_changed(int state)     // If second derivatives are checked forces gradient to be checked
{
    if (state == Qt::Checked) {
        CHKpotgrad->setChecked(true);
    }
}

void MainWindow::CHKpotexact_changed(int state)    // If CHKpotexact is checked unchecks CHKpotgrad and CHKpotder2
{
    if (state == Qt::Checked) {
        CHKpotgrad->setChecked(false);
        CHKpotder2->setChecked(false);
        FRMpotlmaxexp->setVisible(false);
        FRMpotlong->setVisible(false);
        FRMpotderivs->setVisible(false);
        if (activebeware){
            QCheckBox *cb = new QCheckBox(tr("Do not display this message again"));
            QMessageBox msgbox;
            msgbox.setText(tr("Beware that MESP calculation from density matrix and basis set ")+
                    tr("is a lengthy process. \nThis option is intended only for testing purposes. ")+
                    tr("\nIt may take too much time in MESP tabulation of large systems on grids"));
            msgbox.setIcon(QMessageBox::Icon::Information);
            msgbox.addButton(QMessageBox::Ok);
            msgbox.setCheckBox(cb);
            QObject::connect(cb, &QCheckBox::stateChanged, [this](int state){
                if (static_cast<Qt::CheckState>(state) == Qt::CheckState::Checked) {
                    this->activebeware = false;
                }
            });
            msgbox.exec();
        }
    }
    else{
        FRMpotlmaxexp->setVisible(true);
        FRMpotlong->setVisible(true);
        FRMpotderivs->setVisible(true);
    }
}

void MainWindow::CHKpotgrad_changed(int state)     // If second derivatives option is checked unckecks it
{
    if (CHKpotder2->QAbstractButton::isChecked() && state == Qt::Unchecked) {
            CHKpotder2->setChecked(false);
    }
}

void MainWindow::CHKpotinput_changed(int state)
{
    if (state == 0 && RBTpot3D->isChecked()){
        CHKpotmpi->setEnabled(true);
        FRMpotmpi->setVisible(true);
        if (CHKpotmpi->isChecked()){
            LBLpotmpi->setEnabled(true);
            SPBpotmpi->setEnabled(true);
        }
        else{
            LBLpotmpi->setEnabled(false);
            SPBpotmpi->setEnabled(false);
        }
    }
    else{
        CHKpotmpi->setEnabled(false);
        FRMpotmpi->setVisible(false);
        LBLpotmpi->setEnabled(false);
        SPBpotmpi->setEnabled(false);
    }
}

void MainWindow::CHKpotmpi_changed(int state)
{
    if (state != 0 && RBTpot3D->isChecked() && !CHKpotinput->isChecked()){
        CHKpotmpi->setChecked(true);
        LBLpotmpi->setEnabled(true);
        SPBpotmpi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKpotmpi->setChecked(false);
        else{
            CHKpotmpi->setChecked(true);
            CHKpotmpi->setEnabled(false);
        }
        LBLpotmpi->setEnabled(false);
        SPBpotmpi->setEnabled(false);
    }
}

void MainWindow::CHKpotxyz_changed()
{
    if(CHKpotxyz->isChecked()){
        Wtablepot->setEnabled(true);
        Wtablepot->setVisible(true);
    }else{
        Wtablepot->setEnabled(false);
        Wtablepot->setVisible(false);
    }
}

//    Executes external program DAMPOT  (Computes molecular electrostatic potential from the atomic partition)
void MainWindow::execDampot()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,3,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTpotgdampotfilename->text().isEmpty())
            fileNameOut=TXTpotgdampotfilename->text();
    inputdatafile("DAMPOT320.inp","DAMPOTSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMPOT from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMPOT320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKpotinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMPOT320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMPOT320.inp";
    QString stdoutput;
    QString strprocess;
    if (RBTpot3D->isChecked() && CHKpotmpi->isChecked() && CHKpotgrid->isChecked()){
        QString processname;
        if (natom >= 10*SPBpotmpi->value())
            processname = "DAMPOT320_mpi_new.exe";
        else
            processname = "DAMPOT320_mpi.exe";
        QString execName = get_execName(processname, QString("DAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBpotmpi->value()).arg(TXTmpiflags->text()).arg(execName);
        stdoutput = ProjectFolder + fileNameOut;
        if (CHKpotexact->isChecked()){
            stdoutput += "_exact";
        }
        stdoutput += "-DAMPOT320_mpi.out";
    }
    else{
        QString processname = "DAMPOT320.exe";
        QString execName = get_execName(processname, QString("DAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameOut;
        if (CHKpotexact->isChecked()){
            stdoutput += "_exact";
        }
        stdoutput += "-DAMPOT320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 13;
    myProcess->start(strprocess);

    if (RBTpot2D->isChecked() && RBTpotplane->isChecked()){
        QString fileinpstr=ProjectFolder+ProjectName+".xyz";
        QFile fileinp(fileinpstr);
        if (fileinp.open(QFile::ReadOnly | QFile::Text)) {
            potplanecase = get_plane_case(TXTpotplaneA->text().toDouble(), TXTpotplaneB->text().toDouble(), TXTpotplaneC->text().toDouble());
            if (potplanecase > 0 && potplanecase < 8){
                QVector<QString> plane;
                plane << "XY0" << "X0Z" << "0YZ" << "AB0" << "A0C" << "0BC" << "ABC";
//                QString fileoutstr = ProjectFolder+ProjectName+"_"+plane.at(potplanecase-1)+"-v.plane";
                QString fileoutstr = ProjectFolder+fileNameOut;
                if (CHKpotexact->isChecked()){
                    fileoutstr += "_exact";
                }
                fileoutstr = fileoutstr+"_"+plane.at(potplanecase-1)+"-v.plane";
                QFile fileout(fileoutstr);
                if(fileout.open(QFile::Text | QFile::WriteOnly)){
                    QTextStream in(&fileinp); // Buffer for reading from fileinput
                    QTextStream out(&fileout); // Buffer for writing to fileout
                    QString line;
                    Elements *Elem = new Elements();
                    while (!in.atEnd()){
                        line = in.readLine();
#if QT_VERSION < 0x050E00
                        QStringList xyz = line.split(' ',QString::SkipEmptyParts);
#else
                        QStringList xyz = line.split(' ',Qt::SkipEmptyParts);
#endif
                        if (xyz.count() == 4){
                            double x = xyz[1].toDouble() * ANGSTROMTOBOHR;
                            double y = xyz[2].toDouble() * ANGSTROMTOBOHR;
                            double z = xyz[3].toDouble() * ANGSTROMTOBOHR;
                            if (qAbs(QVector3D::dotProduct(QVector3D(x,y,z),
                                QVector3D(TXTpotplaneA->text().toDouble(), TXTpotplaneB->text().toDouble(), TXTpotplaneC->text().toDouble()))) > 1.e-8)
                                continue;
                            double u = QVector3D::dotProduct(QVector3D(x,y,z),wu) / wu.length();
                            double v = QVector3D::dotProduct(QVector3D(x,y,z),wv) / wv.length();
                            int znuc = Elem->getZsymbol(xyz[0]);
                            out << QString("%1  %2  %3  %4  %5  %6\n").arg(u).arg(v).arg(znuc).arg(x).arg(y).arg(z);
                        }
                    }
                    out << QString("%1 %2 %3").arg(TXTpotplaneA->text().toDouble()).arg(TXTpotplaneB->text().toDouble()).arg(TXTpotplaneC->text().toDouble());
                    fileout.close();
                }
            }
        }
    }
}

void MainWindow::potential_resolution_changed(){
    if (RBTpotrcustom->isChecked())
    if (RBTpot3D->isChecked()){
        FRMpotresol3D->setHidden(false);
        FRMpotresol2D->setHidden(true);
    }
    else{
        FRMpotresol2D->setHidden(false);
        FRMpotresol3D->setHidden(true);
    }
    else{
    FRMpotresol2D->setHidden(true);
    FRMpotresol3D->setHidden(true);
    }
}

void MainWindow::RBTpot2D3D_changed()
{
    if (RBTpot2D->isChecked()){
        FRMpotgrid2D->setVisible(true);
        FRMpotgrid3D->setVisible(false);
        CHKpotgrad->setChecked(false);
        CHKpotder2->setChecked(false);
        if (RBTpotrcustom->isChecked()){
            FRMpotresol2D->setHidden(false);
        }
        FRMpotresol3D->setHidden(true);
        RBTpotrlow->setToolTip("129x129");
        RBTpotrmedium->setToolTip("257x257");
        RBTpotrhigh->setToolTip("513x513");
        if (mpi){
            FRMpotmpi->setVisible(false);
            CHKpotmpi->setEnabled(false);
            LBLpotmpi->setEnabled(false);
            SPBpotmpi->setEnabled(false);
        }
    }
    else{
        FRMpotgrid3D->setVisible(true);
        FRMpotgrid2D->setVisible(false);
        if (RBTpotrcustom->isChecked()){
            FRMpotresol3D->setHidden(false);
        }
        FRMpotresol2D->setHidden(true);
        CHKpotgrad->setChecked(true);
        RBTpotrlow->setToolTip("65x65x65");
        RBTpotrmedium->setToolTip("129x129x129");
        RBTpotrhigh->setToolTip("257x257x257");
        if (mpi){
            if (!CHKpotinput->isChecked()){
                FRMpotmpi->setVisible(true);
                CHKpotmpi->setEnabled(true);
                if (CHKpotmpi->isChecked()){
                    LBLpotmpi->setEnabled(true);
                    SPBpotmpi->setEnabled(true);
                }
                else{
                    LBLpotmpi->setEnabled(false);
                    SPBpotmpi->setEnabled(false);
                }
            }
            else{
                FRMpotmpi->setVisible(false);
                CHKpotmpi->setEnabled(false);
                LBLpotmpi->setEnabled(false);
                SPBpotmpi->setEnabled(false);
            }
        }
    }
    CHKpotgrid_changed();
}


void MainWindow::RBTpotplane_changed()
{
    if (RBTpotplane->QAbstractButton::isChecked()) {
        FRMpotsurfpar->setVisible(false);
        FRMpotplane2D->setVisible(true);
        potplanecase = get_plane_case(TXTpotplaneA->text().toDouble(),TXTpotplaneB->text().toDouble(),TXTpotplaneC->text().toDouble());
    }else{
        FRMpotsurfpar->setVisible(true);
        FRMpotplane2D->setVisible(false);
    }
}


void MainWindow::RBTpot2Dplanes_changed(){
    if (RBTpotplaneXY->isChecked()){
        TXTpotxformula2D->setText("u");
        TXTpotyformula2D->setText("v");
        TXTpotzformula2D->setText("0");
        TXTpotplaneA->setText("0.");
        TXTpotplaneB->setText("0.");
        TXTpotplaneC->setText("1.");
        FRMpotplaneABC->setVisible(false);
    }
    else if (RBTpotplaneXZ->isChecked()){
        TXTpotxformula2D->setText("u");
        TXTpotyformula2D->setText("0");
        TXTpotzformula2D->setText("v");
        TXTpotplaneA->setText("0.");
        TXTpotplaneB->setText("1.");
        TXTpotplaneC->setText("0.");
        FRMpotplaneABC->setVisible(false);
    }
    else if (RBTpotplaneYZ->isChecked()){
        TXTpotxformula2D->setText("0");
        TXTpotyformula2D->setText("u");
        TXTpotzformula2D->setText("v");
        TXTpotplaneA->setText("1.");
        TXTpotplaneB->setText("0.");
        TXTpotplaneC->setText("0.");
        FRMpotplaneABC->setVisible(false);
    }
    else if (RBTpotplaneABC->isChecked()){
        TXTpotxformula2D->setText("u");
        TXTpotyformula2D->setText("v");
        TXTpotzformula2D->setText("-(("+TXTpotplaneA->text()+"*u)+("+TXTpotplaneB->text()+"*v))/("+TXTpotplaneC->text()+")");
        FRMpotplaneABC->setVisible(true);
    }
    potplanecase = get_plane_case(TXTpotplaneA->text().toDouble(),TXTpotplaneB->text().toDouble(),TXTpotplaneC->text().toDouble());
}


void MainWindow::rename_pot_cntfile(){
    QString aux;
    if (CHKpotexact->isChecked()){
        aux = "_exact";
    }
    else{
        aux = "";
    }
    QFile filecnt(ProjectFolder + ProjectName + aux + "-v.cnt");
    if (filecnt.exists() && planesuffix(potplanecase) != ""){
        QFile fileold(ProjectFolder + ProjectName + aux + planesuffix(potplanecase) + "-v.cnt");
        if (fileold.exists())
            fileold.remove();
        filecnt.rename(ProjectFolder + ProjectName + aux + planesuffix(potplanecase) + "-v.cnt");
    }
}

void MainWindow::SPBpotmpi_changed(int nprocessors)
{
    SPBpotmpi->setValue(nprocessors);
}


/***************************************************************************/
/*  page_MO: MOLECULAR ORBITALS                                            */
/***************************************************************************/


void MainWindow::CHKMOgrid_changed(){
    if(CHKMOgrid->isChecked()){
        if(RBTMO2D->isChecked()){
            FRMMOgrid2D->setVisible(true);
            TXTMOuinf->setEnabled(true);
            TXTMOusup->setEnabled(true);
            TXTMOvinf->setEnabled(true);
            TXTMOvsup->setEnabled(true);
            FRMMOgrid3D->setVisible(false);
            TXTMOxinf->setEnabled(false);
            TXTMOxsup->setEnabled(false);
            TXTMOyinf->setEnabled(false);
            TXTMOysup->setEnabled(false);
            TXTMOzinf->setEnabled(false);
            TXTMOzsup->setEnabled(false);
            if (RBTMOrcustom->isChecked()){
                FRMMOresol2D->setHidden(false);
            }
            FRMMOresol3D->setHidden(true);
            if (mpi)
                FRMMOmpi->setVisible(false);
        }
        else{
            FRMMOgrid2D->setVisible(false);
            TXTMOuinf->setEnabled(false);
            TXTMOusup->setEnabled(false);
            TXTMOvinf->setEnabled(false);
            TXTMOvsup->setEnabled(false);
            FRMMOgrid3D->setVisible(true);
            TXTMOxinf->setEnabled(true);
            TXTMOxsup->setEnabled(true);
            TXTMOyinf->setEnabled(true);
            TXTMOysup->setEnabled(true);
            TXTMOzinf->setEnabled(true);
            TXTMOzsup->setEnabled(true);
            if (RBTdensrcustom->isChecked()){
                FRMMOresol3D->setHidden(false);
            }
            FRMMOresol2D->setHidden(true);
            if (mpi){
                if( !CHKMOinput->isChecked())
                    FRMMOmpi->setVisible(true);
                else
                    FRMMOmpi->setVisible(false);
            }
        }
        FRMMOgridtype->setVisible(true);
        FRMMOgridres->setVisible(true);
        RBTMO2D->setEnabled(true);
        RBTMO3D->setEnabled(true);
        RBTMOrcustom->setEnabled(true);
        RBTMOrlow->setEnabled(true);
        RBTMOrmedium->setEnabled(true);
        RBTMOrhigh->setEnabled(true);
    }
    else{
        if (mpi)
            FRMMOmpi->setVisible(false);
        FRMMOgridtype->setVisible(false);
        FRMMOgridres->setVisible(false);
        FRMMOgrid2D->setVisible(false);
        TXTMOuinf->setEnabled(false);
        TXTMOusup->setEnabled(false);
        TXTMOvinf->setEnabled(false);
        TXTMOvsup->setEnabled(false);
        FRMMOgrid3D->setVisible(false);
        TXTMOxinf->setEnabled(false);
        TXTMOxsup->setEnabled(false);
        TXTMOyinf->setEnabled(false);
        TXTMOysup->setEnabled(false);
        TXTMOzinf->setEnabled(false);
        TXTMOzsup->setEnabled(false);
        RBTMO2D->setEnabled(false);
        RBTMO3D->setEnabled(false);
        RBTMOrcustom->setEnabled(false);
        RBTMOrlow->setEnabled(false);
        RBTMOrmedium->setEnabled(false);
        RBTMOrhigh->setEnabled(false);
    }
    if (RBTMO3D->isChecked()){
        RBTMOrlow->setToolTip("65x65x65");
        RBTMOrmedium->setToolTip("129x129x129");
        RBTMOrhigh->setToolTip("257x257x257");
    }
    else{
        RBTMOrlow->setToolTip("129x129");
        RBTMOrmedium->setToolTip("257x257");
        RBTMOrhigh->setToolTip("513x513");
    }
}

void MainWindow::CHKMOinput_changed(int state)     
{
    if (state == 0 && RBTMO3D->isChecked()){
        CHKMOmpi->setEnabled(true);
        FRMMOmpi->setVisible(true);
        if (CHKMOmpi->isChecked()){
            LBLMOmpi->setEnabled(true);
            SPBMOmpi->setEnabled(true);
        }
        else{
            LBLMOmpi->setEnabled(false);
            SPBMOmpi->setEnabled(false);
        }
    }
    else{
        FRMMOmpi->setVisible(false);
        CHKMOmpi->setEnabled(false);
        LBLMOmpi->setEnabled(false);
        SPBMOmpi->setEnabled(false);
    }
}

void MainWindow::CHKMOmpi_changed(int state)     
{
    if (state != 0 && RBTMO3D->isChecked() && !CHKMOinput->isChecked()){
        CHKMOmpi->setChecked(true);
        LBLMOmpi->setEnabled(true);
        SPBMOmpi->setEnabled(true);
    }
    else{
        if (state == 0) 
            CHKMOmpi->setChecked(false);
        else{
            CHKMOmpi->setChecked(true);
            CHKMOmpi->setEnabled(false);
        }
        LBLMOmpi->setEnabled(false);
        SPBMOmpi->setEnabled(false);
    }
}

//    Executes external program DAMORB  (Computes molecular orbitals)
void MainWindow::execDamorb()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    if (TXTMOImportfile->text().isEmpty()){
        QMessageBox::warning(this, tr("DAMQT"),tr("Choose a file with molecular orbitals coefficients"));
        return;
    }
    existsinp(DirNombreArchivo,8,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTMOfilename->text().isEmpty())
            fileNameOut=TXTMOfilename->text();
    inputdatafile("DAMORB320.inp","DAMORBSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMORB from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMORB320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKMOinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMORB320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMORB320.inp";
    QString stdoutput;
    QString strprocess;
    if (RBTMO3D->isChecked() && CHKMOmpi->isChecked()){
        QString processname = "DAMORB320_mpi.exe";
        QString execName = get_execName(processname, QString("DAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBMOmpi->value()).arg(TXTmpiflags->text()).arg(execName);
        stdoutput = ProjectFolder + fileNameOut + "-DAMORB320_mpi.out";
    }
    else{
        QString processname = "DAMORB320.exe";
        QString execName = get_execName(processname, QString("DAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameOut + "-DAMORB320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 18;    
    myProcess->start(strprocess);
}

/* Import name of a file with molecular orbitals*/
void MainWindow::ImportFileNameMO()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
            tr("Import data from")+" (*.orb *.GAorba *.GAorbb *.SLorba *.SLorbb);;"+
            tr("All files")+" (*)");
    if (fileName.length()==0){
        return;
    }
    TXTMOImportfile->setText(FileWithoutPath(fileName));
}

void MainWindow::MO_resolution_changed(){
    if (RBTMOrcustom->isChecked())
    if (RBTMO3D->isChecked()){
        FRMMOresol3D->setHidden(false);
        FRMMOresol2D->setHidden(true);
    }
    else{
        FRMMOresol2D->setHidden(false);
        FRMMOresol3D->setHidden(true);
    }
    else{
    FRMMOresol2D->setHidden(true);
    FRMMOresol3D->setHidden(true);
    }
}

void MainWindow::RBTMO2D3D_changed()
{
    if (RBTMO2D->isChecked()){
        FRMMOgrid2D->setVisible(true);
        FRMMOgrid3D->setVisible(false);
        CHKMOgrad->setChecked(false);
        if (RBTMOrcustom->isChecked()){
            FRMMOresol2D->setHidden(false);
        }
        FRMMOresol3D->setHidden(true);
        if (mpi){
            FRMMOmpi->setVisible(false);
            CHKMOmpi->setEnabled(false);
            LBLMOmpi->setEnabled(false);
            SPBMOmpi->setEnabled(false);
        }
    }
    else{
        FRMMOgrid3D->setVisible(true);
        FRMMOgrid2D->setVisible(false);
        if (RBTMOrcustom->isChecked()){
            FRMMOresol3D->setHidden(false);
        }
        FRMMOresol2D->setHidden(true);
        if (mpi){
            if (!CHKMOinput->isChecked()){
                CHKMOmpi->setEnabled(true);
                FRMMOmpi->setVisible(true);
                if (CHKMOmpi->isChecked()){
                    LBLMOmpi->setEnabled(true);
                    SPBMOmpi->setEnabled(true);
                }
                else{
                    LBLMOmpi->setEnabled(false);
                    SPBMOmpi->setEnabled(false);
                }
            }
            else{
                FRMMOmpi->setVisible(false);
                CHKMOmpi->setEnabled(false);
                LBLMOmpi->setEnabled(false);
                SPBMOmpi->setEnabled(false);
            }
        }
        CHKMOgrad->setChecked(true);
    }
    if (RBTMO3D->isChecked()){
        RBTMOrlow->setToolTip("65x65x65");
        RBTMOrmedium->setToolTip("129x129x129");
        RBTMOrhigh->setToolTip("257x257x257");
    }
    else{
        RBTMOrlow->setToolTip("129x129");
        RBTMOrmedium->setToolTip("257x257");
        RBTMOrhigh->setToolTip("513x513");
    }
    CHKMOgrid_changed();
}


void MainWindow::SPBMOmpi_changed(int nprocessors)
{
    SPBMOmpi->setValue(nprocessors);
}

void MainWindow::TXTMOchoose_changed(){
    QString string = QString(TXTMOchoose->text());
    QStringList stringlist1 = string.split(",");
    MOlist->clear();
    for (int i = 0 ; i < stringlist1.length() ; ++i){
        if (stringlist1.at(i).length() != 0){
            if (stringlist1.at(i).contains(QString("-"))){
                QStringList stringlist2 = stringlist1.at(i).split("-");
                if (stringlist2.at(0).length() != 0 && stringlist2.at(1).length() != 0 
                        && stringlist2.at(1).toInt() > stringlist2.at(0).toInt()){
                    for (int k = stringlist2.at(0).toInt() ; k <= stringlist2.at(1).toInt() ; k++){
                        MOlist->append(QString("%1").arg(k));
                    }
                }
            }
            else{
                MOlist->append(stringlist1.at(i));
            }
        }
    }
    MOlist->sort();
    MOlist->removeDuplicates();
}

/****************************************************************************/
/*  page_orimult: ORIENTED MULTIPOLES                                       */
/****************************************************************************/

//    Executes external program DAMMULTROT  (Orients multipoles in a frame with the Z axis orthogonal to the plane defined by three selected atoms)
void MainWindow::execDammultrot()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,7,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTmrotorimultfilename->text().isEmpty())
            fileNameOut=TXTmrotorimultfilename->text();
    inputdatafile("DAMMULTROT320.inp","DAMMULTROTSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMMULTROT from .damproj file
    QString strprocess;
    QString processname = "DAMMULTROT320.exe";
    QString execName = get_execName(processname, QString("DAM320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    QString stdinput = ProjectFolder + ProjectName + "-DAMMULTROT320.inp";
    QString stdoutput = ProjectFolder + fileNameOut + "-DAMMULTROT320.out";
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 17;
    myProcess->start(strprocess);
}

//    Executes external program sgbs2sxyz  (extracts the number of centers and geometry from .sgbs file to file .sxyz)
void MainWindow::execsgbs2sxyz(QString sxyzfilename)
{
    QString fileoutstr = sxyzfilename;
    fileoutstr.replace(".sxyz",".tmpinp");
    QFile fileout(fileoutstr);
    if (!fileout.isOpen()){
        fileout.open(QFile::Text | QFile::WriteOnly);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    string vproj = "\"" + (QFileInfo(sxyzfilename).path() +"/"+ QFileInfo(sxyzfilename).completeBaseName()).toStdString() + "\"" ;
    outfile << vproj.c_str();
#if QT_VERSION < 0x050E00
    outfile << endl;
#else
    outfile << Qt::endl;
#endif
    fileout.close();
    QString processname = "sgbs2sxyz.exe";
    QString strprocess;
    QString execName = get_execName(processname, QString("DAM320"));
    if (execName.isEmpty())
        return;
    strprocess = QString(execName);
    QString stdinput = fileoutstr;
    QString stderror = fileoutstr.replace(".tmpinp",".err");
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stderror,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stderror,QIODevice::Append);
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    executing = 6;
    myProcess->start(strprocess);
    myProcess->waitForFinished(20000);
    fileout.remove() ;
    QString importfilename = ProjectFolder + QFileInfo(ImportFile).completeBaseName();
    if ((importfilename+".sxyz") != sxyzfilename){
        QFile(importfilename+".sxyz").rename(sxyzfilename);
        QFile(importfilename+".sgbs2sxyz").rename(ProjectFolder+ProjectName+".sgbs2sxyz");
    }
}

void MainWindow::RBTMOplane_changed()
{
    if (RBTMOplane->QAbstractButton::isChecked()) {
        FRMMOsurfpar->setVisible(false);
        FRMMOplane2D->setVisible(true);
        MOplanecase = get_plane_case(TXTMOplaneA->text().toDouble(),TXTMOplaneB->text().toDouble(),TXTMOplaneC->text().toDouble());
    }else{
        FRMMOsurfpar->setVisible(true);
        FRMMOplane2D->setVisible(false);
    }
}

void MainWindow::RBTMO2Dplanes_changed(){
    if (RBTMOplaneXY->isChecked()){
        TXTMOxformula2D->setText("u");
        TXTMOyformula2D->setText("v");
        TXTMOzformula2D->setText("0");
        TXTMOplaneA->setText("0.");
        TXTMOplaneB->setText("0.");
        TXTMOplaneC->setText("1.");
        FRMMOplaneABC->setVisible(false);
    }
    else if (RBTMOplaneXZ->isChecked()){
        TXTMOxformula2D->setText("u");
        TXTMOyformula2D->setText("0");
        TXTMOzformula2D->setText("v");
        TXTMOplaneA->setText("0.");
        TXTMOplaneB->setText("1.");
        TXTMOplaneC->setText("0.");
        FRMMOplaneABC->setVisible(false);
    }
    else if (RBTMOplaneYZ->isChecked()){
        TXTMOxformula2D->setText("0");
        TXTMOyformula2D->setText("u");
        TXTMOzformula2D->setText("v");
        TXTMOplaneA->setText("1.");
        TXTMOplaneB->setText("0.");
        TXTMOplaneC->setText("0.");
        FRMMOplaneABC->setVisible(false);
    }
    else if (RBTMOplaneABC->isChecked()){
        TXTMOxformula2D->setText("u");
        TXTMOyformula2D->setText("v");
        TXTMOzformula2D->setText("-(("+TXTMOplaneA->text()+"*u)+("+TXTMOplaneB->text()+"*v))/("+TXTMOplaneC->text()+")");
        FRMMOplaneABC->setVisible(true);
    }
    MOplanecase = get_plane_case(TXTMOplaneA->text().toDouble(),TXTMOplaneB->text().toDouble(),TXTMOplaneC->text().toDouble());
}

void MainWindow::SPBmrotleft_changed()
{
    if (SPBmrotleft->value() > spbmrotleft){    // The value of SPBmrotleft->value() is increasing
        if (SPBmrotleft->value() == SPBmrotmiddle->value()){
            if(SPBmrotright->value() != (SPBmrotmiddle->value()+1)){
                if (SPBmrotmiddle->value() < get_natom())
                    SPBmrotleft->setValue(SPBmrotmiddle->value()+1);
                else
                    SPBmrotleft->setValue(spbmrotleft);
            }
            else if(SPBmrotright->value() < get_natom())
                SPBmrotleft->setValue(SPBmrotright->value()+1);
            else
                SPBmrotleft->setValue(spbmrotleft);
        }
        else if(SPBmrotleft->value()==SPBmrotright->value()){
            if(SPBmrotmiddle->value() != SPBmrotright->value()+1){
                if (SPBmrotright->value() < get_natom())
                    SPBmrotleft->setValue(SPBmrotright->value()+1);
                else
                    SPBmrotleft->setValue(spbmrotleft);
            }
            else if(SPBmrotmiddle->value() < get_natom())
                SPBmrotleft->setValue(SPBmrotmiddle->value()+1);
            else
                SPBmrotleft->setValue(spbmrotleft);
        }
    }
    else{    // The value of SPBmrotleft->value() is decreasing
        if (SPBmrotleft->value() == SPBmrotmiddle->value()){
            if(SPBmrotright->value() != SPBmrotmiddle->value()-1){
                if (SPBmrotmiddle->value() > 1)
                    SPBmrotleft->setValue(SPBmrotmiddle->value()-1);
                else
                    SPBmrotleft->setValue(spbmrotleft);
            }
            else if(SPBmrotright->value() > 1)
                SPBmrotleft->setValue(SPBmrotright->value()-1);
            else
                SPBmrotleft->setValue(spbmrotleft);
        }
        else if(SPBmrotleft->value()==SPBmrotright->value()){
            if(SPBmrotmiddle->value() != SPBmrotright->value()-1){
                if (SPBmrotright->value() > 1)
                    SPBmrotleft->setValue(SPBmrotright->value()-1);
                else
                    SPBmrotleft->setValue(spbmrotleft);
            }
            else if(SPBmrotmiddle->value() > 1)
                SPBmrotleft->setValue(SPBmrotmiddle->value()-1);
            else
                SPBmrotleft->setValue(spbmrotleft);
        }
    }
    spbmrotleft = SPBmrotleft->value();
}

void MainWindow::SPBmrotlmax_changed()
{
    if (SPBmrotlmin->value() > SPBmrotlmax->value())
        SPBmrotlmax->setValue(SPBmrotlmin->value());
    SPBmrotlmin->setMaximum(SPBmrotlmax->value());
}

void MainWindow::SPBmrotlmin_changed()
{
    if (SPBmrotlmin->value() > SPBmrotlmax->value())
        SPBmrotlmin->setValue(SPBmrotlmax->value());
    SPBmrotlmax->setMinimum(SPBmrotlmin->value());
}

void MainWindow::SPBmrotmiddle_changed()
{
    if (SPBmrotmiddle->value() > spbmrotmiddle){    // The value of SPBmrotmiddle->value() is increasing
        if (SPBmrotmiddle->value() == SPBmrotleft->value()){
            if(SPBmrotright->value() != SPBmrotleft->value()+1){
                if (SPBmrotleft->value() < get_natom())
                    SPBmrotmiddle->setValue(SPBmrotleft->value()+1);
                else
                    SPBmrotmiddle->setValue(spbmrotmiddle);
            }
            else if(SPBmrotright->value() < get_natom())
                SPBmrotmiddle->setValue(SPBmrotright->value()+1);
            else
                SPBmrotmiddle->setValue(spbmrotmiddle);
        }
        else if(SPBmrotmiddle->value()==SPBmrotright->value()){
            if(SPBmrotleft->value() != SPBmrotright->value()+1){
                if (SPBmrotright->value() < get_natom())
                    SPBmrotmiddle->setValue(SPBmrotright->value()+1);
                else
                    SPBmrotmiddle->setValue(spbmrotmiddle);
            }
            else if(SPBmrotleft->value() < get_natom())
                SPBmrotmiddle->setValue(SPBmrotleft->value()+1);
            else
                SPBmrotmiddle->setValue(spbmrotmiddle);
        }
    }
    else{    // The value of SPBmrotmiddle->value() is decreasing
        if (SPBmrotmiddle->value() == SPBmrotleft->value()){
            if(SPBmrotright->value() != SPBmrotleft->value()-1){
                if (SPBmrotleft->value() > 1)
                    SPBmrotmiddle->setValue(SPBmrotleft->value()-1);
                else
                    SPBmrotmiddle->setValue(spbmrotmiddle);
            }
            else if(SPBmrotright->value() > 1)
                SPBmrotmiddle->setValue(SPBmrotright->value()-1);
            else
                SPBmrotmiddle->setValue(spbmrotmiddle);
        }
        else if(SPBmrotmiddle->value()==SPBmrotright->value()){
            if(SPBmrotleft->value() != SPBmrotright->value()-1){
                if (SPBmrotright->value() > 1)
                    SPBmrotmiddle->setValue(SPBmrotright->value()-1);
                else
                    SPBmrotmiddle->setValue(spbmrotmiddle);
            }
            else if(SPBmrotleft->value() > 1)
                SPBmrotmiddle->setValue(SPBmrotleft->value()-1);
            else
                SPBmrotmiddle->setValue(spbmrotmiddle);
        }
    }
    spbmrotmiddle = SPBmrotmiddle->value();
}

void MainWindow::SPBmrotright_changed()
{
    if (SPBmrotright->value() > spbmrotright){    // The value of SPBmrotright->value() is increasing
        if (SPBmrotright->value() == SPBmrotleft->value()){
            if(SPBmrotmiddle->value() != SPBmrotleft->value()+1){
                if (SPBmrotleft->value() < get_natom())
                    SPBmrotright->setValue(SPBmrotleft->value()+1);
                else
                    SPBmrotright->setValue(spbmrotright);
            }
            else if(SPBmrotmiddle->value() <get_natom())
                SPBmrotright->setValue(SPBmrotmiddle->value()+1);
            else
                SPBmrotright->setValue(spbmrotright);
        }
        else if(SPBmrotright->value()==SPBmrotmiddle->value()){
            if(SPBmrotleft->value() != SPBmrotmiddle->value()+1){
                if (SPBmrotmiddle->value() < get_natom())
                    SPBmrotright->setValue(SPBmrotmiddle->value()+1);
                else
                    SPBmrotright->setValue(spbmrotright);
            }
            else if(SPBmrotleft->value() < get_natom())
                SPBmrotright->setValue(SPBmrotleft->value()+1);
            else
                SPBmrotright->setValue(spbmrotright);
        }
    }
    else{    // The value of SPBmrotright->value() is decreasing
        if (SPBmrotright->value() == SPBmrotleft->value()){
            if(SPBmrotmiddle->value() != SPBmrotleft->value()-1){
                if (SPBmrotleft->value() > 1)
                    SPBmrotright->setValue(SPBmrotleft->value()-1);
                else
                    SPBmrotright->setValue(spbmrotright);
            }
            else if(SPBmrotmiddle->value() > 1)
                SPBmrotright->setValue(SPBmrotmiddle->value()-1);
            else
                SPBmrotright->setValue(spbmrotright);
        }
        else if(SPBmrotright->value()==SPBmrotmiddle->value()){
            if(SPBmrotleft->value() != SPBmrotmiddle->value()-1){
                if (SPBmrotmiddle->value() > 1)
                    SPBmrotright->setValue(SPBmrotmiddle->value()-1);
                else
                    SPBmrotright->setValue(spbmrotright);
            }
            else if(SPBmrotleft->value() > 1)
                SPBmrotright->setValue(SPBmrotleft->value()-1);
            else
                SPBmrotright->setValue(spbmrotright);
        }
    }
    spbmrotright = SPBmrotright->value();
}

void MainWindow::TXTmrotorimultatoms_changed()
{
    QString string = QString(TXTmrotorimultatoms->text());
    QStringList stringlist1 = string.split(",");
    mrotorimultlist->clear();
    for (int i = 0 ; i < stringlist1.length() ; ++i){
        if (stringlist1.at(i).length() != 0){
            if (stringlist1.at(i).contains(QString("-"))){
                QStringList stringlist2 = stringlist1.at(i).split("-");
                if (stringlist2.at(0).length() != 0 && stringlist2.at(1).length() != 0
                        && stringlist2.at(1).toInt() > stringlist2.at(0).toInt()){
                    for (int k = stringlist2.at(0).toInt() ; k <= stringlist2.at(1).toInt() ; k++){
                        mrotorimultlist->append(QString("%1").arg(k));
                    }
                }
            }
            else{
                mrotorimultlist->append(stringlist1.at(i));
            }
        }
    }
    mrotorimultlist->sort();
    mrotorimultlist->removeDuplicates();
}



/***************************************************************************/
/*               page_SGhole: MESP SIGMA HOLE GENERATION                   */
/***************************************************************************/

void MainWindow::CHKSGhexactMESP_changed(int state)
{
    if (state == 0){
        FRMSGholelmaxexp->setVisible(true);
    }
    else{
        FRMSGholelmaxexp->setVisible(false);
        if (activebeware){
            QCheckBox *cb = new QCheckBox(tr("Do not display this message again"));
            QMessageBox msgbox;
            msgbox.setText(tr("Beware that MESP calculation from density matrix and basis set ")+
                    tr("is a lengthy process. \nThis option is intended only for testing purposes. ")+
                    tr("\nIt may take too much time in MESP tabulation of large systems on grids"));
            msgbox.setIcon(QMessageBox::Icon::Information);
            msgbox.addButton(QMessageBox::Ok);
            msgbox.setCheckBox(cb);
            QObject::connect(cb, &QCheckBox::stateChanged, [this](int state){
                if (static_cast<Qt::CheckState>(state) == Qt::CheckState::Checked) {
                    this->activebeware = false;
                }
            });
            msgbox.exec();
        }
    }
}

void MainWindow::CHKSGholeinput_changed(int state)
{
    if (state == 0){
        CHKSGholempi->setEnabled(true);
        FRMSGholempi->setVisible(true);
        if (CHKSGholempi->isChecked()){
            LBLSGholempi->setEnabled(true);
            SPBSGholempi->setEnabled(true);
        }
        else{
            LBLSGholempi->setEnabled(false);
            SPBSGholempi->setEnabled(false);
        }
    }
    else{
        CHKSGholempi->setEnabled(false);
        FRMSGholempi->setVisible(false);
        LBLSGholempi->setEnabled(false);
        SPBSGholempi->setEnabled(false);
    }
}

void MainWindow::CHKSGholempi_changed(int state)
{
    if (state != 0 && !CHKSGholeinput->isChecked()){
        CHKSGholempi->setChecked(true);
        LBLSGholempi->setEnabled(true);
        SPBSGholempi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKSGholempi->setChecked(false);
        else{
            CHKSGholempi->setChecked(true);
            CHKSGholempi->setEnabled(false);
        }
        LBLSGholempi->setEnabled(false);
        SPBSGholempi->setEnabled(false);
    }
}

//    Executes external program DAMDENZJ  (Computes molecular density or deformations from the Zernike or Jacobi expansion)
void MainWindow::execDamSGhole()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,13,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTSGholefilename->text().isEmpty())
            fileNameOut=TXTSGholefilename->text();
    inputdatafile("DAMSGHOLE320.inp","DAMSGHOLESECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMSGHOLE from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMSGHOLE320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKSGholeinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMSGHOLE320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMSGHOLE320.inp";
    QString stdoutput;
    QString strprocess;
    if (CHKSGholempi->isChecked()){
        QString processname = "DAMSGHOLE320_mpi.exe";
        QString execName = get_execName(processname, QString("DAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBSGholempi->value()).arg(TXTmpiflags->text()).arg(execName);
        stdoutput = ProjectFolder + fileNameOut;
        if (CHKSGhexactMESP->isChecked()){
            stdoutput += "_exact";
        }
        stdoutput += "-DAMSGHOLE320_mpi.out";
    }
    else{
        QString processname = "DAMSGHOLE320.exe";
        QString execName = get_execName(processname, QString("DAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameOut;
        if (CHKSGhexactMESP->isChecked()){
            stdoutput += "_exact";
        }
        stdoutput += "-DAMSGHOLE320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 23;
    myProcess->start(strprocess);
}

/* Import name of a density grid file*/
void MainWindow::importFileSGholeden()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
            tr("Import data from")+" (*-d.plt);;"+
            tr("All files")+" (*)");
    if (fileName.length()==0){
        return;
    }
    TXTImportSGholeden->setText(FileWithoutPath(fileName));
}

void MainWindow::SPBSGholempi_changed(int nprocessors)
{
    SPBSGholempi->setValue(nprocessors);
}

void MainWindow::TXTImportSGholeden_changed(){
    if (TXTImportSGholeden->text().isEmpty()){
        LBLSGholecontour->setEnabled(false);
        TXTSGholecontour->setEnabled(false);
        LBLSGholelocalextrema->setEnabled(false);
        SPBSGholelocalextrema->setEnabled(false);
        LBLSGholelocalpower->setEnabled(false);
        LBLSGholegeomthreshold->setEnabled(false);
        SPBSGholegeomthreshold->setEnabled(false);
        LBLSGholelongthreshold->setEnabled(false);
        SPBSGholelongthreshold->setEnabled(false);
        CHKSGhexactMESP->setEnabled(false);
        FRMSGholelmaxexp->setEnabled(false);
        LBLSGholelmaxexp->setEnabled(false);
        SPBSGholelmaxexp->setEnabled(false);
        FRMSGholeinput->setEnabled(false);
        CHKSGholeinput->setEnabled(false);
        FRMSGholempi->setEnabled(false);
        CHKSGholempi->setEnabled(false);
        LBLSGholempi->setEnabled(false);
        SPBSGholempi->setEnabled(false);
        BTNexecDamSGhole->setEnabled(false);
    }
    else{
        LBLSGholecontour->setEnabled(true);
        TXTSGholecontour->setEnabled(true);
        LBLSGholelocalextrema->setEnabled(true);
        SPBSGholelocalextrema->setEnabled(true);
        LBLSGholelocalpower->setEnabled(true);
        LBLSGholegeomthreshold->setEnabled(true);
        SPBSGholegeomthreshold->setEnabled(true);
        LBLSGholelongthreshold->setEnabled(true);
        SPBSGholelongthreshold->setEnabled(true);
        CHKSGhexactMESP->setEnabled(true);
        FRMSGholelmaxexp->setEnabled(true);
        LBLSGholelmaxexp->setEnabled(true);
        SPBSGholelmaxexp->setEnabled(true);
        FRMSGholeinput->setEnabled(true);
        CHKSGholeinput->setEnabled(true);
        FRMSGholempi->setEnabled(true);
        CHKSGholempi->setEnabled(true);
        if (CHKSGholempi->isChecked()){
            LBLSGholempi->setEnabled(true);
            SPBSGholempi->setEnabled(true);
        }
        else{
            LBLSGholempi->setEnabled(false);
            SPBSGholempi->setEnabled(false);
        }
        BTNexecDamSGhole->setEnabled(true);
    }
}


/***************************************************************************/
/*  page_TOPO: TOPOGRAPHY                                                  */
/***************************************************************************/


void MainWindow::CHKtopobasin_changed(int state)
{
    if (state != 0){
        LBLtopoboxb->setVisible(true);
        TXTtopoboxb->setVisible(true);
        LBLtopofdisp->setVisible(true);
        TXTtopofdisp->setVisible(true);
        CHKtopomolgraph->setChecked(true);
        CHKtopoexdraw->setVisible(true);
        LBLtopoexln->setVisible(true);
        TXTtopoexln->setVisible(true);
    }
    else{
        LBLtopoboxb->setVisible(false);
        TXTtopoboxb->setVisible(false);
        if (!CHKtopomolgraph->isChecked()){
            LBLtopofdisp->setVisible(false);
            TXTtopofdisp->setVisible(false);
        }
        CHKtopoexdraw->setVisible(false);
        LBLtopoexln->setVisible(false);
        TXTtopoexln->setVisible(false);
    }
}

void MainWindow::CHKtopoexdraw_changed(int state)
{
    if (state != 0){
        TXTtopoexln->setEnabled(true);
    }
    else{
        TXTtopoexln->setEnabled(false);
    }
}

void MainWindow::CHKtopomolgraph_changed(int state)
{
    if (state != 0){
        LBLtopoboxg->setVisible(true);
        TXTtopoboxg->setVisible(true);
        LBLtopoggradthr->setVisible(true);
        TXTtopoggradthr->setVisible(true);
        LBLtopofdisp->setVisible(true);
        TXTtopofdisp->setVisible(true);
    }
    else{
        LBLtopoboxg->setVisible(false);
        TXTtopoboxg->setVisible(false);
        LBLtopoggradthr->setVisible(false);
        TXTtopoggradthr->setVisible(false);
        CHKtopobasin->setChecked(false);
        LBLtopofdisp->setVisible(false);
        TXTtopofdisp->setVisible(false);
    }
}

void MainWindow::CHKtopoaddguess_changed(int state)
{
    if (state != 0){
        FRMtopoguessfilename->setVisible(true);
        TXTtopoguessfilename->setVisible(true);
        BTNtopoguessfilename->setVisible(true);
        CHKtopoxyz->setEnabled(true);
        CHKtopoxyz->setVisible(true);
    }
    else{
        FRMtopoguessfilename->setVisible(false);
        TXTtopoguessfilename->setVisible(false);
        BTNtopoguessfilename->setVisible(false);
        CHKtopoxyz->setChecked(false);
        CHKtopoxyz->setEnabled(false);
        CHKtopoxyz->setVisible(false);
    }
}

void MainWindow::CHKtopograph_changed()
{
    if (CHKtopograph->QAbstractButton::isChecked()) {
        TXTtopoboxl->setVisible(true);
        LBLtopoboxl->setVisible(true);
        LBLtopocnvg->setVisible(true);
        TXTtopocnvg->setVisible(true);
        LBLtopoboxt->setVisible(true);
        TXTtopoboxt->setVisible(true);
        LBLtopostepszt->setVisible(true);
        TXTtopostepszt->setVisible(true);
    }else{
        TXTtopoboxl->setVisible(false);
        LBLtopoboxl->setVisible(false);
        LBLtopocnvg->setVisible(false);
        TXTtopocnvg->setVisible(false);
        FRMtopoguessfilename->setVisible(false);
        TXTtopoguessfilename->setVisible(false);
        BTNtopoguessfilename->setVisible(false);
        LBLtopoboxt->setVisible(false);
        TXTtopoboxt->setVisible(false);
        LBLtopostepszt->setVisible(false);
        TXTtopostepszt->setVisible(false);
    }
}

void MainWindow::CHKtopoinput_changed(int state)
{
    if (state == 0){
        CHKtopompi->setEnabled(true);
        if (CHKtopompi->isChecked()){
            LBLtopompi->setEnabled(true);
            SPBtopompi->setEnabled(true);
        }
        else{
            LBLtopompi->setEnabled(false);
            SPBtopompi->setEnabled(false);
        }
    }
    else{
        CHKtopompi->setEnabled(false);
        LBLtopompi->setEnabled(false);
        SPBtopompi->setEnabled(false);
    }
}


void MainWindow::CHKtopompi_changed(int state)
{
    if (state == 0){
        CHKtopompi->setChecked(false);
        LBLtopompi->setEnabled(false);
        SPBtopompi->setEnabled(false);
    }
    else if(!CHKtopoinput->isChecked()){
        CHKtopompi->setChecked(true);
        LBLtopompi->setEnabled(true);
        SPBtopompi->setEnabled(true);
    }
}


void MainWindow::CHKtopoxyz_changed()
{
    if(CHKtopoxyz->isChecked()){
        Wtabledentopo->setEnabled(true);
        Wtabledentopo->setVisible(true);
    }else{
        Wtabledentopo->setEnabled(false);
        Wtabledentopo->setVisible(false);
    }
}

//    Executes external program DAMTOPOGRAPHER  (Carries out topography analysis)
void MainWindow::execDamTopography()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,9,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTtopofilename->text().isEmpty())
            fileNameOut=TXTtopofilename->text();
    inputdatafile("DAMTOPO320.inp","DAMTOPOSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMTOPO from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMTOPO320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKtopoinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMTOPO320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMTOPO320.inp";
    QString stdoutput;
    QString strprocess;
    if (CHKtopompi->isChecked()){
        QString processname = "DAMTOPOGRAPHY_mpi.exe";
        QString execName = get_execName(processname, QString("TDAM320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBtopompi->value()).arg(TXTmpiflags->text()).arg(execName);
        if (RBTtopodensity->isChecked())
            stdoutput = ProjectFolder + fileNameOut + "-DAMTOPO320_mpi-d.out";
        else
            stdoutput = ProjectFolder + fileNameOut + "-DAMTOPO320_mpi-v.out";
    }
    else{
        QString processname = "DAMTOPOGRAPHY.exe";
        QString execName = get_execName(processname, QString("TDAM320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        if (RBTtopodensity->isChecked())
            stdoutput = ProjectFolder + fileNameOut + "-DAMTOPO320-d.out";
        else
            stdoutput = ProjectFolder + fileNameOut + "-DAMTOPO320-v.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 19;
    myProcess->start(strprocess);
}

/* Import name of a file with molecular orbitals*/
void MainWindow::ImportTopguessfilename()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
            tr("Import data from")+" (*.xyz);;"+
            tr("All files")+" (*)");
    if (fileName.length()==0){
        return;
    }
    TXTtopoguessfilename->setText(fileName);
}

void MainWindow::RBTtopodensity_changed()
{
    if (RBTtopodensity->QAbstractButton::isChecked()) {
        TXTtopocnvg->setText("4.0e-15");
        TXTtopoboxl->setText("1.0");
        TXTtopoboxt->setText("2.0");
    }else{
        TXTtopocnvg->setText("4.0e-12");
        TXTtopoboxl->setText("1.0");
        TXTtopoboxt->setText("5.0");
    }
    CHKtopoaddguess->setChecked(false);
    TXTtopoguessfilename->setText("");
}

void MainWindow::SPBtopompi_changed(int nprocessors)
{
    SPBtopompi->setValue(nprocessors);
}


/***************************************************************************/
/*  page_ZJdens: ZERNIKE 3D-JACOBI EXPANSIONS                              */
/***************************************************************************/


void MainWindow::CHKZJinput_changed(int state)
{
    if (state == 0){
        FRMZJmpi->setVisible(true);
        CHKZJmpi->setEnabled(true);
        if (CHKZJmpi->isChecked()){
            LBLZJmpi->setEnabled(true);
            SPBZJmpi->setEnabled(true);
        }
        else{
            LBLZJmpi->setEnabled(false);
            SPBZJmpi->setEnabled(false);
        }
    }
    else{
        FRMZJmpi->setVisible(false);
        CHKZJmpi->setEnabled(false);
        LBLZJmpi->setEnabled(false);
        SPBZJmpi->setEnabled(false);
    }
}

void MainWindow::CHKZJmpi_changed(int state)
{
    if (state != 0 && !CHKZJinput->isChecked()){
        CHKZJmpi->setChecked(true);
        LBLZJmpi->setEnabled(true);
        SPBZJmpi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKZJmpi->setChecked(false);
        else{
            CHKZJmpi->setChecked(true);
            CHKZJmpi->setEnabled(false);
        }
        LBLZJmpi->setEnabled(false);
        SPBZJmpi->setEnabled(false);
    }
}

//    Executes external program DAM (Partition of molecular density into atomic densities)
void MainWindow::execDamZJ(){
        defineRanges();
    if (lslater){
        QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
        existsinp(DirNombreArchivo,10,1,false);
        inputdatafile("DAMZJ320.inp","DAMZJSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAM from .damproj file
        QString fileName=ProjectFolder+ProjectName+"-DAMZJ320.inp";
        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
            return;
        }
        if (CHKZJinput->isChecked()){
            QTextStream in(&file);
            textEdit->setFont(QFont("Courier",10));
            textEdit->setPlainText(in.readAll());
            QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMZJ320.inp"));
            return;
        }
        file.close();
        QString stdinput = ProjectFolder + ProjectName + "-DAMZJ320.inp";
        QString stdoutput;
        QString strprocess;
        if (CHKZJmpi->isChecked()){
            QString processname = "DAMZernike-Jacobi_STO_mpi.exe";
            QString execName = get_execName(processname, QString("DAMZernike320_mpi"));
            if (execName.isEmpty())
                return;
            strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                        .arg(SPBZJmpi->value()).arg(TXTmpiflags->text()).arg(execName);
            stdoutput = ProjectFolder + ProjectName;
            if (RBTZJacobi->isChecked()){
                stdoutput +=  + "-Jacobi-DAMZJ320_mpi.out";
            }
            else{
                stdoutput +=  + "-Zernike-DAMZJ320_mpi.out";
            }
        }
        else{
            QString processname = "DAMZernike-Jacobi_STO.exe";
            QString execName = get_execName(processname, QString("DAMZernike320"));
            if (execName.isEmpty())
                return;
            strprocess = QString(execName);
            stdoutput = ProjectFolder + ProjectName;
            if (RBTZJacobi->isChecked()){
                stdoutput +=  + "-Jacobi-DAMZJ320.out";
            }
            else{
                stdoutput +=  + "-Zernike-DAMZJ320.out";
            }
        }
        myProcess = new QProcess(this);
        myProcess->setStandardInputFile(stdinput);
        myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
        myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
        connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
        connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
        connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
        executing = 20;
        myProcess->start(strprocess);
    }
    else{
        QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
        existsinp(DirNombreArchivo,10,1,false);
        inputdatafile("DAMZJ320.inp","DAMZJSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for GDAM from .damproj file
        QString fileName=ProjectFolder+ProjectName+"-DAMZJ320.inp";
        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
            return;
        }
        if (CHKZJinput->isChecked()){
            QTextStream in(&file);
            textEdit->setFont(QFont("Courier",10));
            textEdit->setPlainText(in.readAll());
            QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMZJ320.inp"));
            return;
        }
        file.close();
        QString stdinput = ProjectFolder + ProjectName + "-DAMZJ320.inp";
        QString stdoutput;
        QString strprocess;

        if (CHKZJmpi->isChecked()){
            QString processname = "DAMZernike-Jacobi_GTO_mpi.exe";
            QString execName = get_execName(processname, QString("DAMZernike320_mpi"));
            if (execName.isEmpty())
                return;
            strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                        .arg(SPBZJmpi->value()).arg(TXTmpiflags->text()).arg(execName);
            stdoutput = ProjectFolder + ProjectName;
            if (RBTZJacobi->isChecked()){
                stdoutput +=  + "-Jacobi-DAMZJ320_mpi.out";
            }
            else{
                stdoutput +=  + "-Zernike-DAMZJ320_mpi.out";
            }
        }
        else{
            QString processname = "DAMZernike-Jacobi_GTO.exe";
            QString execName = get_execName(processname, QString("DAMZernike320"));
            if (execName.isEmpty())
                return;
            strprocess = QString(execName);
            stdoutput = ProjectFolder + ProjectName;
            if (RBTZJacobi->isChecked()){
                stdoutput +=  + "-Jacobi-DAMZJ320.out";
            }
            else{
                stdoutput +=  + "-Zernike-DAMZJ320.out";
            }
        }
        myProcess = new QProcess(this);
        myProcess->setStandardInputFile(stdinput);
        myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
        myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
        connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
        connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
        connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
        executing = 20;
        myProcess->start(strprocess);
    }
}

void MainWindow::SPBZJkmax_changed()
{
    int valor=SPBZJkmax->value();
    int valor1=SPBdZJkmax->value();

    if (valor < valor1) SPBdZJkmax->setValue(valor);
    SPBdZJkmax->setRange(0,valor);
}

void MainWindow::SPBZJlmax_changed()
{
    int valor=SPBZJlmax->value();
    int valor1=SPBdZJlmax->value();

    if (valor < valor1) SPBdZJlmax->setValue(valor);
    SPBdZJlmax->setRange(0,valor);

    SPBdZJlmin->setRange(0,valor);
}

void MainWindow::SPBZJmpi_changed(int nprocessors)
{
    SPBZJmpi->setValue(nprocessors);
}


void MainWindow::TXTZJrstar_changed()
{
    rstar = TXTZJrstar->text().toDouble();
}


/***************************************************************************/
/*  page_ZJtab: DENSITY TABULATION FROM ZERNIKE 3D-JACOBI EXPANSIONS       */
/***************************************************************************/


void MainWindow::CHKdZJgrid_changed(){
    if(CHKdZJgrid->isChecked()){
        if(RBTdZJ2D->isChecked()){
            FRMdZJgrid2D->setVisible(true);
            TXTdZJuinf->setEnabled(true);
            TXTdZJusup->setEnabled(true);
            TXTdZJvinf->setEnabled(true);
            TXTdZJvsup->setEnabled(true);
            FRMdZJgrid3D->setVisible(false);
            TXTdZJxinf->setEnabled(false);
            TXTdZJxsup->setEnabled(false);
            TXTdZJyinf->setEnabled(false);
            TXTdZJysup->setEnabled(false);
            TXTdZJzinf->setEnabled(false);
            TXTdZJzsup->setEnabled(false);
            if (RBTdZJrcustom->isChecked()){
                FRMdZJresol2D->setHidden(false);
            }
            FRMdZJresol3D->setHidden(true);
            if (mpi)
                FRMdZJmpi->setVisible(false);
        }
        else{
            FRMdZJgrid2D->setVisible(false);
            TXTdZJuinf->setEnabled(false);
            TXTdZJusup->setEnabled(false);
            TXTdZJvinf->setEnabled(false);
            TXTdZJvsup->setEnabled(false);
            FRMdZJgrid3D->setVisible(true);
            TXTdZJxinf->setEnabled(true);
            TXTdZJxsup->setEnabled(true);
            TXTdZJyinf->setEnabled(true);
            TXTdZJysup->setEnabled(true);
            TXTdZJzinf->setEnabled(true);
            TXTdZJzsup->setEnabled(true);
            if (RBTdensrcustom->isChecked()){
                FRMdZJresol3D->setHidden(false);
            }
            FRMdZJresol2D->setHidden(true);
            if (mpi){
                if( !CHKdZJinput->isChecked())
                    FRMdZJmpi->setVisible(true);
                else
                    FRMdZJmpi->setVisible(false);
            }
        }
        FRMdZJgridtype->setVisible(true);
        FRMdZJgridres->setVisible(true);
        RBTdZJ2D->setEnabled(true);
        RBTdZJ3D->setEnabled(true);
        RBTdZJrcustom->setEnabled(true);
        RBTdZJrlow->setEnabled(true);
        RBTdZJrmedium->setEnabled(true);
        RBTdZJrhigh->setEnabled(true);
    }
    else{
        if (mpi)
            FRMdZJmpi->setVisible(false);
        FRMdZJgridtype->setVisible(false);
        FRMdZJgridres->setVisible(false);
        FRMdZJgrid2D->setVisible(false);
        TXTdZJuinf->setEnabled(false);
        TXTdZJusup->setEnabled(false);
        TXTdZJvinf->setEnabled(false);
        TXTdZJvsup->setEnabled(false);
        FRMdZJgrid3D->setVisible(false);
        TXTdZJxinf->setEnabled(false);
        TXTdZJxsup->setEnabled(false);
        TXTdZJyinf->setEnabled(false);
        TXTdZJysup->setEnabled(false);
        TXTdZJzinf->setEnabled(false);
        TXTdZJzsup->setEnabled(false);
        RBTdZJ2D->setEnabled(false);
        RBTdZJ3D->setEnabled(false);
        RBTdZJrcustom->setEnabled(false);
        RBTdZJrlow->setEnabled(false);
        RBTdZJrmedium->setEnabled(false);
        RBTdZJrhigh->setEnabled(false);
    }
    if (RBTdZJ3D->isChecked()){
        RBTdZJrlow->setToolTip("65x65x65");
        RBTdZJrmedium->setToolTip("129x129x129");
        RBTdZJrhigh->setToolTip("257x257x257");
    }
    else{
        RBTdZJrlow->setToolTip("129x129");
        RBTdZJrmedium->setToolTip("257x257");
        RBTdZJrhigh->setToolTip("513x513");
    }
}

void MainWindow::CHKdZJinput_changed(int state)
{
    if (state == 0 && RBTdZJ3D->isChecked()){
        FRMdZJmpi->setVisible(true);
        CHKdZJmpi->setEnabled(true);
        if (CHKdZJmpi->isChecked()){
            LBLdZJmpi->setEnabled(true);
            SPBdZJmpi->setEnabled(true);
        }
        else{
            LBLdZJmpi->setEnabled(false);
            SPBdZJmpi->setEnabled(false);
        }
    }
    else{
        FRMdZJmpi->setVisible(false);
        CHKdZJmpi->setEnabled(false);
        LBLdZJmpi->setEnabled(false);
        SPBdZJmpi->setEnabled(false);
    }
}

void MainWindow::CHKdZJmpi_changed(int state)
{
    if (state != 0 && RBTdZJ3D->isChecked() && !CHKdZJinput->isChecked()){
        CHKdZJmpi->setChecked(true);
        LBLdZJmpi->setEnabled(true);
        SPBdZJmpi->setEnabled(true);
    }
    else{
        if (state == 0)
            CHKdZJmpi->setChecked(false);
        else{
            CHKdZJmpi->setChecked(true);
            CHKdZJmpi->setEnabled(false);
        }
        LBLdZJmpi->setEnabled(false);
        SPBdZJmpi->setEnabled(false);
    }
}

void MainWindow::CHKdZJxyz_changed()
{
    if(CHKdZJxyz->isChecked()){
        WtabledZJ->setEnabled(true);
        WtabledZJ->setVisible(true);
    }else{
        WtabledZJ->setEnabled(false);
        WtabledZJ->setVisible(false);
    }
}

void MainWindow::ChooseZJ_changed()
{
    TXTdZJchoose->clear();
    if (RBTdZJchooseall->isChecked()){
    LBLdZJchoose->setText("");
    LBLdZJchoose->setVisible(false);
    TXTdZJchoose->setVisible(false);
    FRMdZJexplength->setVisible(true);
    return;
    }
    LBLdZJchoose->setVisible(true);
    TXTdZJchoose->setVisible(true);
    FRMdZJexplength->setVisible(false);
    if (RBTdZJchoosel->isChecked() || RBTdZJchoosek->isChecked()){
    QRegExp rxZF("[0-9][-,\\d]*");
    QValidator *ZJvalidator = new QRegExpValidator(rxZF, this);
    LBLdZJchoose->setText(tr("Select indices: 1,3-5,10,13-17,..."));
    TXTdZJchoose->setValidator(ZJvalidator);
    }
    else if (RBTdZJchooselk->isChecked()){
    QRegExp rxZF("^\\(?[0-9]*[\\(,\\d*\\)]*");
    QValidator *ZJvalidator = new QRegExpValidator(rxZF, this);
    LBLdZJchoose->setText(tr("Select pairs of indices: (0,0),(1,2),..."));
    TXTdZJchoose->setValidator(ZJvalidator);
    }
    else{
    QRegExp rxZF("^\\(?[0-9]*[\\(,-\\d*\\)]*");
    QValidator *ZJvalidator = new QRegExpValidator(rxZF, this);
    LBLdZJchoose->setText(tr("Select triads of indices:(0,0,0),(1,2,-1),..."));
    TXTdZJchoose->setValidator(ZJvalidator);
    }
}

//    Executes external program DAMDENZJ  (Computes molecular density or deformations from the Zernike or Jacobi expansion)
void MainWindow::execDamdenZJ()
{
    QString DirNombreArchivo = ProjectFolder+ProjectName+".damproj";
    existsinp(DirNombreArchivo,11,1,false);
    QString fileNameOut = ProjectName;
    if (!TXTdZJfilename->text().isEmpty())
            fileNameOut=TXTdZJfilename->text();
    inputdatafile("DAMDENZJ320.inp","DAMDENZJSECT",DirNombreArchivo, ProjectFolder, ProjectName); // Creates a file with input data for DAMDENZJ from .damproj file
    QString fileName=ProjectFolder+ProjectName+"-DAMDENZJ320.inp";
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Cannot create input file %1").arg(fileName));
        return;
    }
    if (CHKdZJinput->isChecked()){
        QTextStream in(&file);
        textEdit->setFont(QFont("Courier",10));
        textEdit->setPlainText(in.readAll());
        QMessageBox::information(this, tr("DAMQT"),tr("File %1 has been successfully created").arg(ProjectFolder+ProjectName+"-DAMDENZJ320.inp"));
        return;
    }
    file.close();
    QString stdinput = ProjectFolder + ProjectName + "-DAMDENZJ320.inp";
    QString stdoutput;
    QString strprocess;
    QString type = QFileInfo(TXTdZJImportfile->text()).suffix();
    if (type=="jacobi"){
        ldZJjacobi = true;
    }
    else{
        ldZJjacobi = false;
    }
    if (RBTdZJ3D->isChecked() && CHKdZJmpi->isChecked()){
        QString processname = "DAMDENZJ320_mpi.exe";
        QString execName = get_execName(processname, QString("DAMZernike320_mpi"));
        if (execName.isEmpty())
            return;
        strprocess = QString("%1 %2 %3 %4 %5").arg(TXTmpicommand->text()).arg("-np")
                    .arg(SPBdZJmpi->value()).arg(TXTmpiflags->text()).arg(execName);

        stdoutput = ProjectFolder + fileNameOut + "-" + type + "-DAMDENZJ320_mpi.out";
    }
    else{
        QString processname = "DAMDENZJ320.exe";
        QString execName = get_execName(processname, QString("DAMZernike320"));
        if (execName.isEmpty())
            return;
        strprocess = QString(execName);
        stdoutput = ProjectFolder + fileNameOut + "-" + type + "-DAMDENZJ320.out";
    }
    myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(processOutput(int, QProcess::ExitStatus)));
    executing = 21;
    myProcess->start(strprocess);
}

/* Import name of a file with Zernike or Jacobi expansion*/
void MainWindow::ImportFileNameZJ()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
            tr("Import data from")+" (*.zernike *.jacobi);;"+
            tr("All files")+" (*)");
    if (fileName.length()==0){
        return;
    }
    TXTdZJImportfile->setText(FileWithoutPath(fileName));
}

void MainWindow::RBTdZJ2D3D_changed()
{
    if (RBTdZJ2D->isChecked()){
        FRMdZJgrid2D->setVisible(true);
        FRMdZJgrid3D->setVisible(false);
        CHKdZJgrad->setChecked(false);
        if (RBTdZJrcustom->isChecked()){
            FRMdZJresol2D->setHidden(false);
        }
        FRMdZJresol3D->setHidden(true);
        if (mpi){
            FRMdZJmpi->setVisible(false);
            CHKdZJmpi->setEnabled(false);
            LBLdZJmpi->setEnabled(false);
            SPBdZJmpi->setEnabled(false);
        }
    }
    else{
        FRMdZJgrid3D->setVisible(true);
        FRMdZJgrid2D->setVisible(false);
        if (RBTdZJrcustom->isChecked()){
            FRMdZJresol3D->setHidden(false);
        }
        FRMdZJresol2D->setHidden(true);
        if (mpi){
            if (!CHKdZJinput->isChecked()){
                FRMdZJmpi->setVisible(true);
                CHKdZJmpi->setEnabled(true);
                if (CHKdZJmpi->isChecked()){
                    LBLdZJmpi->setEnabled(true);
                    SPBdZJmpi->setEnabled(true);
                }
                else{
                    LBLdZJmpi->setEnabled(false);
                    SPBdZJmpi->setEnabled(false);
                }
            }
            else{
                FRMdZJmpi->setVisible(false);
                CHKdZJmpi->setEnabled(false);
                LBLdZJmpi->setEnabled(false);
                SPBdZJmpi->setEnabled(false);
            }
        }
        CHKdZJgrad->setChecked(true);
    }
    if (RBTdZJ3D->isChecked()){
        RBTdZJrlow->setToolTip("65x65x65");
        RBTdZJrmedium->setToolTip("129x129x129");
        RBTdZJrhigh->setToolTip("257x257x257");
    }
    else{
        RBTdZJrlow->setToolTip("129x129");
        RBTdZJrmedium->setToolTip("257x257");
        RBTdZJrhigh->setToolTip("513x513");
    }
    CHKdZJgrid_changed();
}

void MainWindow::RBTdZJplane_changed()
{
    if (RBTdZJplane->QAbstractButton::isChecked()) {
        FRMdZJsurfpar->setVisible(false);
        FRMdZJplane2D->setVisible(true);
        dZJplanecase = get_plane_case(TXTdZJplaneA->text().toDouble(),TXTdZJplaneB->text().toDouble(),TXTdZJplaneC->text().toDouble());
    }else{
        FRMdZJsurfpar->setVisible(true);
        FRMdZJplane2D->setVisible(false);
    }
}


void MainWindow::RBTdZJ2Dplanes_changed(){
    if (RBTdZJplaneXY->isChecked()){
        TXTdZJxformula2D->setText("u");
        TXTdZJyformula2D->setText("v");
        TXTdZJzformula2D->setText("0");
        TXTdZJplaneA->setText("0.");
        TXTdZJplaneB->setText("0.");
        TXTdZJplaneC->setText("1.");
        FRMdZJplaneABC->setVisible(false);
    }
    else if (RBTdZJplaneXZ->isChecked()){
        TXTdZJxformula2D->setText("u");
        TXTdZJyformula2D->setText("0");
        TXTdZJzformula2D->setText("v");
        TXTdZJplaneA->setText("0.");
        TXTdZJplaneB->setText("1.");
        TXTdZJplaneC->setText("0.");
        FRMdZJplaneABC->setVisible(false);
    }
    else if (RBTdZJplaneYZ->isChecked()){
        TXTdZJxformula2D->setText("0");
        TXTdZJyformula2D->setText("u");
        TXTdZJzformula2D->setText("v");
        TXTdZJplaneA->setText("1.");
        TXTdZJplaneB->setText("0.");
        TXTdZJplaneC->setText("0.");
        FRMdZJplaneABC->setVisible(false);
    }
    else if (RBTdZJplaneABC->isChecked()){
        TXTdZJxformula2D->setText("u");
        TXTdZJyformula2D->setText("v");
        TXTdZJzformula2D->setText("-(("+TXTdZJplaneA->text()+"*u)+("+TXTdZJplaneB->text()+"*v))/("+TXTdZJplaneC->text()+")");
        FRMdZJplaneABC->setVisible(true);
    }
    dZJplanecase = get_plane_case(TXTdZJplaneA->text().toDouble(),TXTdZJplaneB->text().toDouble(),TXTdZJplaneC->text().toDouble());
}


void MainWindow::SPBdZJmpi_changed(int nprocessors)
{
    SPBdZJmpi->setValue(nprocessors);
}

void MainWindow::TXTdZJchoose_changed(){
    QString string = QString(TXTdZJchoose->text());
    string.remove("(",Qt::CaseInsensitive);
    string.remove(")",Qt::CaseInsensitive);
#if QT_VERSION < 0x050E00
    QStringList stringlist1 = string.split(' ',QString::SkipEmptyParts);
#else
    QStringList stringlist1 = string.split(' ',Qt::SkipEmptyParts);
#endif
    ZJlist->clear();
    for (int i = 0 ; i < stringlist1.length() ; ++i){
        if (stringlist1.at(i).length() != 0){
            if ((RBTdZJchoosel->isChecked() || RBTdZJchoosek->isChecked()) && stringlist1.at(i).contains(QString("-"))){
                QStringList stringlist2 = stringlist1.at(i).split("-");
                if (stringlist2.at(0).length() != 0 && stringlist2.at(1).length() != 0
                        && stringlist2.at(1).toInt() > stringlist2.at(0).toInt()){
                    for (int k = stringlist2.at(0).toInt() ; k <= stringlist2.at(1).toInt() ; k++){
                        ZJlist->append(QString("%1").arg(k));
                    }
                }
            }
            else{
                ZJlist->append(stringlist1.at(i));
            }
        }
    }
    if (RBTdZJchoosel->isChecked() || RBTdZJchoosek->isChecked()){
        ZJlist->sort();
        ZJlist->removeDuplicates();
    }
}

void MainWindow::ZJ_resolution_changed(){
    if (RBTdZJrcustom->isChecked())
    if (RBTdZJ3D->isChecked()){
        FRMdZJresol3D->setHidden(false);
        FRMdZJresol2D->setHidden(true);
    }
    else{
        FRMdZJresol2D->setHidden(false);
        FRMdZJresol3D->setHidden(true);
    }
    else{
    FRMdZJresol2D->setHidden(true);
    FRMdZJresol3D->setHidden(true);
    }
}


/**************************************************************************************************/
/********** FUNCTIONS FOR READING/WRITING OPTIONS FROM/TO .damproj FILE ***************************/
/**************************************************************************************************/


void MainWindow::read_CHK(const char * a, const char * b, QCheckBox * c, const string file){
    QString qv;
    string v = CIniFile::GetValue(a,b,file);
    qv = toQString(v.c_str());
    if (qv==QString("T") || qv==QString("t")){
        c->setChecked(true);
    }else{
        c->setChecked(false);
    }
}

void MainWindow::read_double_to_TXT(const char * a, const char * b, QLineEdit * c, const string file){
    bool ok;
    QString qv;
    string v = CIniFile::GetValue(a,b,file);
    qv = toQString(v.c_str());
    qv.toDouble(&ok);
    if (!qv.isEmpty() && ok) c->setText(qv);
}

void MainWindow::read_plot_dimension(const char * a, const char * b, const string file,
        QRadioButton * d2, QRadioButton * l, QRadioButton * m, QRadioButton * h){
    QString qv;
    string v = CIniFile::GetValue(a,b,file);
    qv = toQString(v.c_str());
    if (qv==QString("T") || qv==QString("t")){
        d2->setChecked(true);
    }else{
        d2->setChecked(false);
    }
    if (d2->isChecked()){
        l->setToolTip("129x129");
        m->setToolTip("257x257");
        h->setToolTip("513x513");   
        }
    else{
        l->setToolTip("65x65x65");
        m->setToolTip("129x129x129");
        h->setToolTip("257x257x257");
    }
}

void MainWindow::read_RBT(const char * a, const char * b, QRadioButton * c, const string file){
    QString qv;
    string v = CIniFile::GetValue(a,b,file);
    qv = toQString(v.c_str());
    if (qv==QString("T") || qv==QString("t")){
        c->setChecked(true);
    }else{
        c->setChecked(false);
    }
}

void MainWindow::read_resolution(const char * a1, const char * a2, const char * a3, const char * b, const string file, 
    const QLineEdit * xi, const QLineEdit * xf, const QLineEdit * yi, const QLineEdit * yf, const QLineEdit * zi, const QLineEdit * zf, 
    QRadioButton * d3, QRadioButton * l, QRadioButton * m, QRadioButton * h, QRadioButton * c){
    bool ok;
    double dltx, dlty, dltz;
    int depthx, depthy, depthz ;
    QString qv;
    string v;
    v = CIniFile::GetValue(a1,b,file);
    qv = toQString(v.c_str());
    dltx=qv.toDouble(&ok);
    v = CIniFile::GetValue(a2,b,file);
    qv = toQString(v.c_str());
    dlty=qv.toDouble(&ok);
    v = CIniFile::GetValue(a3,b,file);
    qv = toQString(v.c_str());
    dltz=qv.toDouble(&ok);
    l->setChecked(true); //2**6  low  (default)
    if (ok && dltx > 0. && dlty > 0. && dltz > 0.){
    depthx = round((xf->text().toDouble() - xi->text().toDouble()) / dltx);
    depthy = round((yf->text().toDouble() - yi->text().toDouble()) / dlty);
    depthz = round((zf->text().toDouble() - zi->text().toDouble()) / dltz);
    if (depthx == depthy && depthx == depthz){
        if (depthx == 512)
        h->setChecked(true);
        else if (depthx == 256){
        if(d3->isChecked())
            h->setChecked(true); //2**8  3D high
        else
            m->setChecked(true); //2**8  2D medium
        }
        else if (depthx == 129){
        if(d3->isChecked())
            m->setChecked(true); //2**8  2D medium
        }
        else{
        c->setChecked(true);
        if (strcmp(b,"DAMDENSECT") == 0){
            if (d3->isChecked()){
            SPBdensxres->setValue(depthx);
            SPBdensyres->setValue(depthy);
            SPBdenszres->setValue(depthz);
            }
            else{
            SPBdensures->setValue(depthx);
            SPBdensvres->setValue(depthy);
            }                
        }
        else if (strcmp(b,"DAMPOTSECT") == 0){
            if (d3->isChecked()){
            SPBpotxres->setValue(depthx);
            SPBpotyres->setValue(depthy);
            SPBpotzres->setValue(depthz);
            }
            else{
            SPBpotures->setValue(depthx);
            SPBpotvres->setValue(depthy);
            }                
        }
        else if (strcmp(b,"DAMORBSECT") == 0){
            if (d3->isChecked()){
            SPBMOxres->setValue(depthx);
            SPBMOyres->setValue(depthy);
            SPBMOzres->setValue(depthz);
            }
            else{
            SPBMOures->setValue(depthx);
            SPBMOvres->setValue(depthy);
            }                
        }
        else if (strcmp(b,"DAMDENZJSECT") == 0){
            if (d3->isChecked()){
            SPBdZJxres->setValue(depthx);
            SPBdZJyres->setValue(depthy);
            SPBdZJzres->setValue(depthz);
            }
            else{
            SPBdZJures->setValue(depthx);
            SPBdZJvres->setValue(depthy);
            }                
        }
        }
    }
    else{
        c->setChecked(true);
        if (strcmp(b,"DAMDENSECT") == 0){
        if (d3->isChecked()){
            SPBdensxres->setValue(depthx);
            SPBdensyres->setValue(depthy);
            SPBdenszres->setValue(depthz);
        }
        else{
            SPBdensures->setValue(depthx);
            SPBdensvres->setValue(depthy);
        }                
        }
        else if (strcmp(b,"DAMPOTSECT") == 0){
        if (d3->isChecked()){
            SPBpotxres->setValue(depthx);
            SPBpotyres->setValue(depthy);
            SPBpotzres->setValue(depthz);
        }
        else{
            SPBpotures->setValue(depthx);
            SPBpotvres->setValue(depthy);
        }                
        }
        else if (strcmp(b,"DAMORBSECT") == 0){
        if (d3->isChecked()){
            SPBMOxres->setValue(depthx);
            SPBMOyres->setValue(depthy);
            SPBMOzres->setValue(depthz);
        }
        else{
            SPBMOures->setValue(depthx);
            SPBMOvres->setValue(depthy);
        }                
        }
        else if (strcmp(b,"DAMDENZJSECT") == 0){
        if (d3->isChecked()){
            SPBdZJxres->setValue(depthx);
            SPBdZJyres->setValue(depthy);
            SPBdZJzres->setValue(depthz);
        }
        else{
            SPBdZJures->setValue(depthx);
            SPBdZJvres->setValue(depthy);
        }                
        }
    }
    };
}

void MainWindow::read_SPB(const char * a, const char * b, QSpinBox * c, const string file){
    bool ok;
    QString qv;
    string v = CIniFile::GetValue(a,b,file);
    qv = toQString(v.c_str());
    qv.toInt(&ok);
    if (ok) c->setValue(qv.toInt());
}

void MainWindow::read_text_to_TXT(const char * a, const char * b, QLineEdit * c, const string file){
    QString qv;
    string v = CIniFile::GetValue(a,b,file);
    qv = toQString(v.c_str());
    qv.remove("\"");   // To prevent duplicate quotation marks
    if (!qv.isEmpty()) c->setText(qv);
}

void MainWindow::write_intervals(const char * a, const char * b, const char * c, const char * d, const string file, 
        QLineEdit * xi, QLineEdit * xs, bool * p, QString * w){
    double ini, fin, dlt;
    ini = xi->text().toDouble();      
    fin = xs->text().toDouble();
    if (strcmp(d,"DAMDENSECT") == 0){
        if (fin > ini) dlt = set_delta(c,ini,fin);
        else if (fin < ini) dlt = set_delta(c,fin,ini);
        else dlt = 1.0;
    }
    else if (strcmp(d,"DAMPOTSECT") == 0){
        if (fin > ini) dlt = set_deltapot(c,ini,fin);
        else if (fin < ini) dlt = set_deltapot(c,fin,ini);
        else dlt = 1.0;
    }    
    else if (strcmp(d,"DAMORBSECT") == 0){
        if (fin > ini) dlt = set_deltaorb(c,ini,fin);
        else if (fin < ini) dlt = set_deltaorb(c,fin,ini);
        else dlt = 1.0;
    }
    else if (strcmp(d,"DAMDENZJSECT") == 0){
        if (fin > ini) dlt = set_deltaZJ(c,ini,fin);
        else if (fin < ini) dlt = set_deltaZJ(c,fin,ini);
        else dlt = 1.0;
    }
    else dlt = 1.0;
    if ( fin >= ini ){
        write_option(a, d, QString::number(ini), file, p, w);
        write_option(b, d, QString::number(fin), file, p, w);
    }
    else{
        write_option(a, d, QString::number(fin), file, p, w);
        write_option(b, d, QString::number(ini), file, p, w);
        xi->setText(QString().setNum(fin));
        xs->setText(QString().setNum(ini));
    }
    if (fin != ini) write_option(c, d, QString::number(dlt), file, p, w);
}

void MainWindow::write_option(const char * a, const char * b, const QString c, const string file, bool * p, QString * w){
    if (!CIniFile::SetValue(a, toString(c), b, file)){ 
        w->append(a).append(" in ").append(b).append("\n");
        *p = true;
    }
}


/**************************************************************************************************/
/***************************************  ANCILLARY FUNCTIONS *************************************/
/**************************************************************************************************/


//    Creates a folder. Returns whether it has succeded or not
bool MainWindow::createDir(QString &fullPathName)
{
    QDir dir(fullPathName);
    if(!dir.exists(fullPathName)){
        dir.mkpath(fullPathName);
        return true;
    }else{
        return false;
    }
}

//    Sets plot ranges
void MainWindow::defineRanges()
{
    QVector<double> x;
    QVector<double> y;
    QVector<double> z;
    QVector<int> ncarga;
    int nat;
    double min;
    double max;
    QString qv;
    readGeometry(nat,x,y,z,ncarga);
    if (nat == 0) return;
    set_natom(nat);
    dminmax(x,min,max);
    xmax = qRound(max);
    xmin = qRound(min);
    dminmax(y,min,max);
    ymax = qRound(max);
    ymin = qRound(min);
    dminmax(z,min,max);
    zmax = qRound(max);
    zmin = qRound(min);
    double xyztop;
    QVector<double> vaux;
    vaux << xmax << ymax << zmax << std::abs(xmin) << std::abs(ymin) << std::abs(zmin);
    dmax(vaux,xyztop);
    xmin = -xyztop;
    ymin = -xyztop;
    zmin = -xyztop;
    xmax =  xyztop;
    ymax =  xyztop;
    zmax =  xyztop;
    TXTdensxinf->setText(qv.setNum(xmin-5,'g',3));
    TXTdensxsup->setText(qv.setNum(xmax+5,'g',3));
    TXTdensuinf->setText(qv.setNum(xmin-5,'g',3));
    TXTdensusup->setText(qv.setNum(xmax+5,'g',3));
    TXTefxinf->setText(qv.setNum(2.0*(xmin-5),'g',3));
    TXTefxsup->setText(qv.setNum(2.0*(xmax+5),'g',3));
    TXTdgxinf->setText(qv.setNum(2.0*(xmin-5),'g',3));
    TXTdgxsup->setText(qv.setNum(2.0*(xmax+5),'g',3));
    TXTpotxinf->setText(qv.setNum(2.0*(xmin-5),'g',3));
    TXTpotxsup->setText(qv.setNum(2.0*(xmax+5),'g',3));
    TXTMOxinf->setText(qv.setNum(2.0*(xmin-5),'g',3));
    TXTMOxsup->setText(qv.setNum(2.0*(xmax+5),'g',3));
    TXTdensyinf->setText(qv.setNum(ymin-5,'g',3));
    TXTdensysup->setText(qv.setNum(ymax+5,'g',3));
    TXTefyinf->setText(qv.setNum(2.0*(ymin-5),'g',3));
    TXTefysup->setText(qv.setNum(2.0*(ymax+5),'g',3));
    TXTdgyinf->setText(qv.setNum(2.0*(ymin-5),'g',3));
    TXTdgysup->setText(qv.setNum(2.0*(ymax+5),'g',3));
    TXTpotyinf->setText(qv.setNum(2.0*(ymin-5),'g',3));
    TXTpotysup->setText(qv.setNum(2.0*(ymax+5),'g',3));
    TXTMOyinf->setText(qv.setNum(2.0*(ymin-5),'g',3));
    TXTMOysup->setText(qv.setNum(2.0*(ymax+5),'g',3));
    TXTdenszinf->setText(qv.setNum(zmin-5,'g',3));
    TXTdenszsup->setText(qv.setNum(zmax+5,'g',3));
    TXTdensvinf->setText(qv.setNum(zmin-5,'g',3));
    TXTdensvsup->setText(qv.setNum(zmax+5,'g',3));
    TXTefzinf->setText(qv.setNum(2.0*(zmin-5)));
    TXTefzsup->setText(qv.setNum(2.0*(zmax+5)));
    TXTdgzinf->setText(qv.setNum(2.0*(zmin-5)));
    TXTdgzsup->setText(qv.setNum(2.0*(zmax+5)));
    TXTpotzinf->setText(qv.setNum(2.0*(zmin-5),'g',3));
    TXTpotzsup->setText(qv.setNum(2.0*(zmax+5),'g',3));
    TXTMOzinf->setText(qv.setNum(2.0*(zmin-5),'g',3));
    TXTMOzsup->setText(qv.setNum(2.0*(zmax+5),'g',3));
    SPBmrotleft->setRange(0,get_natom());
    SPBmrotleft->setValue(0);
    spbmrotleft = 0;
    SPBmrotleft->setEnabled(true);
    SPBmrotmiddle->setRange(0,get_natom());
    SPBmrotmiddle->setValue(0);
    spbmrotmiddle = 0;
    SPBmrotmiddle->setEnabled(true);
    SPBmrotright->setRange(0,get_natom());
    SPBmrotright->setValue(0);
    spbmrotright = 0;
    SPBmrotright->setEnabled(true);
}


void MainWindow::disable_pages(){
    page_atdens->setEnabled(false);
    page_densgrad->setEnabled(false);
    page_Efield->setEnabled(false);
    page_frad->setEnabled(false);
    page_HFforces->setEnabled(false);
    page_MED->setEnabled(false);
    page_MESP->setEnabled(false);
    page_orimult->setEnabled(false);
    page_SGhole->setEnabled(false);
    page_TOPO->setEnabled(false);
    page_ZJdens->setEnabled(false);
    page_ZJtab->setEnabled(false);
    page_MO->setEnabled(false);
}

// Determines the highest value in an array
void MainWindow::dmax(QVector<double> &v,double &max)
{
    max=v[0];
    for (int i=0;i<v.size();i++){
        if (v[i]>max) max=v[i];
    }
}

// Determines the lowest value in an array
void MainWindow::dmin(QVector<double> &v,double &min)
{
    min=v[0];
    for (int i=0;i<v.size();i++){
        if (v[i]<min) min=v[i];
    }
}

// Determines the highest and the lowest values in an array
void MainWindow::dminmax(QVector<double> &v,double &min,double &max)
{
    min=v[0];
    max=v[0];
    for (int i=0;i<v.size();i++){
        if (v[i]<min) min=v[i];
        if (v[i]>max) max=v[i];
    }    
}

//    Checks whether a project file (.damproj) exist or not
bool MainWindow::existsinp(QString fullinputName, int tab, int def, bool pregunta)
{
    if (fullinputName.size()==0){
        QMessageBox::warning(this, tr("DAMQT"),tr("Introduce options file name")+" (*.damproj)");
        return false;
    }
    if (pregunta==true){
        if (!QFile::exists(fullinputName))  {
            loadDefault(def);
            saveOptions(fullinputName,tab);
        }else{
            QMessageBox msgBox;
            msgBox.setInformativeText(QString(tr("Options file %1 exists")).arg(fullinputName)
                        +"\n"+tr("Do you want to overwrite?"));
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
            msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
            msgBox.setButtonText(QMessageBox::No, tr("No"));
            msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
            msgBox.setIcon(QMessageBox::Warning);
            int ret = msgBox.exec();
            if (ret == QMessageBox::Yes){
                saveOptions(fullinputName,tab);
            }else if (ret == QMessageBox::No){
                readOptions(fullinputName);
            }else if (ret == QMessageBox::Cancel){
                return false;
            }
        }
    }else{
        saveOptions(fullinputName,tab);
    }
    return true;
}

//    gets content of variable natom
int MainWindow::get_natom()
{
    return MainWindow::natom;
}

int MainWindow::get_plane_case(double a, double b, double c){
    // Return plane case and vectors wu and wv
    //      Case 1: A == 0, B == 0, C != 0
    //      Case 2: A == 0, B != 0, C == 0
    //      Case 3: A != 0, B == 0, C == 0
    //      Case 4: A != 0, B != 0, C == 0
    //      Case 5: A != 0, B == 0, C != 0
    //      Case 6: A == 0, B != 0, C != 0
    //      Case 7: A != 0, b != 0, C != 0
    //      Case -1: Error
    if (std::abs(a) < 1.e-7){
        if (std::abs(b) < 1.e-7){
            if (std::abs(c) < 1.e-7){
                QMessageBox::warning(this, tr("DAMQT"),tr("Parameters A=%1, B=%2, C=%3 do not define a plane")
                    .arg(a).arg(b).arg(c));
                return -1;
            }
            else{       // A == 0, B == 0, C != 0
                wu = QVector3D(1.0, 0.0, 0.0);
                wv = QVector3D(0.0, 1.0, 0.0);
                return 1;
            }
        }
        else{
            if (std::abs(c) < 1.e-7){   // A == 0, B != 0, C == 0
                wu = QVector3D(1.0, 0.0, 0.0);
                wv = QVector3D(0.0, 0.0, 1.0);
                return 2;
            }
            else{   // A == 0, B != 0, C != 0
                wu = QVector3D(-1.0, 0.0, 0.0);
                wv = QVector3D(0.0, -c, b) / QVector3D(0.0,b,c).length();
                return 6;
            }
        }
    }
    else{
        if (std::abs(b) < 1.e-7){
            if (std::abs(c) < 1.e-7){   // A != 0, B == 0, C == 0
                wu = QVector3D(0.0, 1.0, 0.0);
                wv = QVector3D(0.0, 0.0, 1.0);
                return 3;
            }
            else{   // A != 0, B == 0, C != 0
                wu = QVector3D(0.0, 1.0, 0.0);
                wv = QVector3D(-c, 0.0, a) / QVector3D(a,0.0,c).length();
                return 5;
            }
        }
        else{
            if (std::abs(c) < 1.e-7){   // A != 0, B != 0, C == 0
                wu = QVector3D(-b, a, 0.0) / QVector3D(a,b,0.0).length();
                wv = QVector3D(0.0, 0.0, 1.0);
                return 4;
            }
            else{   // A != 0, b != 0, C != 0
                wu = QVector3D(-b, a, 0.0) / QVector3D(a,b,0.0).length();
                wv = QVector3D(-a*c, -b*c, a*a+b*b) / (QVector3D(a,b,0.0).length() * QVector3D(a,b,c).length());
                return 7;
            }
        }
    }
}


QString MainWindow::planesuffix(int planecase){
    switch (planecase){
    case 1:
        return QString("_XY0");
    case 2:
        return QString("_X0Z");
    case 3:
        return QString("_0YZ");
    case 4:
        return QString("_AB0");
    case 5:
        return QString("_A0C");
    case 6:
        return QString("_0BC");
    case 7:
        return QString("_ABC");
    default:
        return QString("");
    }
}

//    Imports an output data file (*.out)
void MainWindow::importOUT()
{
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString fileName = filedialog.getOpenFileName(this,tr("Open output file ..."),ProjectFolder,
        tr("Output files")+" (*.out *.log);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;

    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("importOUT"),tr("File %1 cannot be read")
            .arg(fileName)+QString(":\n%1.").arg(file.errorString()));
        return;
    }
    QTextStream in(&file);
    textEdit->setFont(QFont("Courier",10));
    textEdit->setPlainText(in.readAll());
}

//    Iitializes pointers
void MainWindow::initpointers()
{
    Acc2Dplot = nullpointer;
    Acc3Dview = nullpointer;
    AccAbout = nullpointer;
    AccAboutQt = nullpointer;
    AccExit = nullpointer;
    AccHelp = nullpointer;
    AccNew = nullpointer;
    AccOpen = nullpointer;
    AccPdf = nullpointer;
    AccPrint = nullpointer;
    AccSave = nullpointer;
    AccSaveAs = nullpointer;

    BTNdgfilelines = nullpointer;
    BTNdZJImportFile = nullpointer;
    BTNeffilelines = nullpointer;
    BTNexecDam = nullpointer;
    BTNexecDamden = nullpointer;
    BTNexecDamdengrad = nullpointer;
    BTNexecDamdZJ = nullpointer;
    BTNexecDamfield = nullpointer;
    BTNexecDamforces = nullpointer;
    BTNexecDamfrad = nullpointer;
    BTNexecDammultrot = nullpointer;
    BTNexecDamorb = nullpointer;
    BTNexecDampot = nullpointer;
    BTNexecDamSGhole = nullpointer;
    BTNexecDamtopo = nullpointer;
    BTNexecDamZJ = nullpointer;
    BTNexecImport = nullpointer;
    BTNImport = nullpointer;
    BTNImportSGholeden = nullpointer;
    BTNlangstart = nullpointer;
    BTNMOImportFile = nullpointer;
    BTNnewplot = nullpointer;
    BTNnewwidget = nullpointer;
    BTNraiseviewers = nullpointer;
    BTNtopoguessfilename = nullpointer;

    CHKatdensmpi = nullpointer;
    CHKatdensinput = nullpointer;
    CHKdensder2 = nullpointer;
    CHKdensgrad = nullpointer;
    CHKdensgrid = nullpointer;
    CHKdensinput = nullpointer;
    CHKdenslaplacian = nullpointer;
    CHKdenslatomics = nullpointer;
    CHKdensldensacc = nullpointer;
    CHKdenslmolec = nullpointer;
    CHKdensmpi = nullpointer;
    CHKdensxyz = nullpointer;
    CHKdgextralines = nullpointer;
    CHKdginput = nullpointer;
    CHKdgmpi = nullpointer;
    CHKdgxyz = nullpointer;
    CHKdZJgrad = nullpointer;
    CHKdZJgrid = nullpointer;
    CHKdZJinput = nullpointer;
    CHKdZJmpi = nullpointer;
    CHKdZJxyz = nullpointer;
    CHKefinput = nullpointer;
    CHKefmpi = nullpointer;
    CHKefextralines = nullpointer;
    CHKeflong = nullpointer;
    CHKefxyz = nullpointer;
    CHKfradextras = nullpointer;
    CHKfradderiv1 = nullpointer;
    CHKfradderiv2 = nullpointer;
    CHKHFlatomsel = nullpointer;
    CHKMOgrid = nullpointer;
    CHKMOgrad = nullpointer;
    CHKMOmpi = nullpointer;
    CHKMOinput = nullpointer;
    CHKpotgrid = nullpointer;
    CHKpotlong = nullpointer;
    CHKpotmpi = nullpointer;
    CHKpotder2 = nullpointer;
    CHKpotexact = nullpointer;
    CHKpotgrad = nullpointer;
    CHKpotinput = nullpointer;
    CHKpotxyz = nullpointer;
    CHKSGhexactMESP = nullpointer;
    CHKSGholeinput = nullpointer;
    CHKSGholempi = nullpointer;
    CHKtopobasin = nullpointer;
    CHKtopoexdraw = nullpointer;
    CHKtopomolgraph = nullpointer;
    CHKtopoinput = nullpointer;
    CHKtopompi = nullpointer;
    CHKtopoaddguess = nullpointer;
    CHKtopograph = nullpointer;
    CHKtopoxyz = nullpointer;
    CHKZJinput = nullpointer;
    CHKZJmpi = nullpointer;
    CMBdgdirectionset = nullpointer;
    CMBefdirectionset = nullpointer;

    denslist = nullpointer;
    densvalidator = nullpointer;
    dialog = nullpointer;
    dockright = nullpointer;
    dockwidget = nullpointer;

    FileMenu = nullpointer;
    fradlist = nullpointer;
    fradvalidator = nullpointer;
    FRM2Dwidgets = nullpointer;
    FRM3Dwidgets = nullpointer;
    FRMatdensmpi = nullpointer;
    FRMatdensinput = nullpointer;
    FRMatdenslmaxdisp = nullpointer;
    FRMatdenslmaxexp = nullpointer;
    FRMatdenstype = nullpointer;
    FRMdensatoms = nullpointer;
    FRMdensdeformation = nullpointer;
    FRMdensderivs = nullpointer;
    FRMdensexpansion = nullpointer;
    FRMdensdamdenfilename = nullpointer;
    FRMdensgrid = nullpointer;
    FRMdensgrid2D = nullpointer;
    FRMdensgrid3D = nullpointer;
    FRMdensgridres = nullpointer;
    FRMdensgridtype = nullpointer;
    FRMdensinput = nullpointer;
    FRMdensdensity = nullpointer;
    FRMdensfragments = nullpointer;
    FRMdenslrange = nullpointer;
    FRMdensmpi = nullpointer;
    FRMdensplane2D = nullpointer;
    FRMdensplaneABC = nullpointer;
    FRMdensresol2D = nullpointer;
    FRMdensresol3D = nullpointer;
    FRMdenssurfpar = nullpointer;
    FRMdenssurftype = nullpointer;
    FRMdensxyz = nullpointer;
    FRMdgextralines = nullpointer;
    FRMdgfilename = nullpointer;
    FRMdginput = nullpointer;
    FRMdglmaxexp = nullpointer;
    FRMdglines = nullpointer;
    FRMdglong = nullpointer;
    FRMdgmpi = nullpointer;
    FRMdgplane2D = nullpointer;
    FRMdgplaneABC = nullpointer;
    FRMdgplot2D = nullpointer;
    FRMdgplot3D = nullpointer;
    FRMdgplottype = nullpointer;
    FRMdguv = nullpointer;
    FRMdgxyz = nullpointer;
    FRMdZJderivs = nullpointer;
    FRMdZJexpansion = nullpointer;
    FRMdZJexplength = nullpointer;
    FRMdZJgrid2D = nullpointer;
    FRMdZJgrid3D = nullpointer;
    FRMdZJgridres = nullpointer;
    FRMdZJgrid = nullpointer;
    FRMdZJgridtype = nullpointer;
    FRMdZJfilename = nullpointer;
    FRMdZJImportfile = nullpointer;
    FRMdZJinput = nullpointer;
    FRMdZJmpi = nullpointer;
    FRMdZJplane2D = nullpointer;
    FRMdZJplaneABC = nullpointer;
    FRMdZJresol2D = nullpointer;
    FRMdZJresol3D = nullpointer;
    FRMdZJsurfpar = nullpointer;
    FRMdZJsurftype = nullpointer;
    FRMdZJxyz = nullpointer;
    FRMeffilename = nullpointer;
    FRMefinput = nullpointer;
    FRMeflines = nullpointer;
    FRMeflmaxexp = nullpointer;
    FRMeflong = nullpointer;
    FRMefmpi = nullpointer;
    FRMefextralines = nullpointer;
    FRMefplane2D = nullpointer;
    FRMefplaneABC = nullpointer;
    FRMefplot2D = nullpointer;
    FRMefplot3D = nullpointer;
    FRMefplottype = nullpointer;
    FRMefuv = nullpointer;
    FRMefxyz = nullpointer;
    FRMfradatoms = nullpointer;
    FRMfraddamfilename = nullpointer;
    FRMfradialfactors = nullpointer;
    FRMfradios = nullpointer;
    FRMfradderivadas = nullpointer;
    FRMHFatomsforces = nullpointer;
    FRMHFdensfragments1 = nullpointer;
    FRMHFgdamforcesfilename = nullpointer;
    FRMlanguage = nullpointer;
    FRMMOchoose = nullpointer;
    FRMMOgridtype = nullpointer;
    FRMMOImportfile = nullpointer;
    FRMMOfilename = nullpointer;
    FRMMOderivs = nullpointer;
    FRMMOsurfpar = nullpointer;
    FRMMOgrid = nullpointer;
    FRMMOgrid2D = nullpointer;
    FRMMOgrid3D = nullpointer;
    FRMMOgridres = nullpointer;
    FRMMOinput = nullpointer;
    FRMMOmpi = nullpointer;
    FRMMOplane2D = nullpointer;
    FRMMOplaneABC = nullpointer;
    FRMMOsurftype = nullpointer;
    FRMMOresol2D = nullpointer;
    FRMMOresol3D = nullpointer;
    FRMmpioptions = nullpointer;
    FRMmrotorimultfilename = nullpointer;
    FRMmrotmultipoles = nullpointer;
    FRMmrotcenters = nullpointer;
    FRMmrotatoms = nullpointer;
    FRMplots = nullpointer;
    FRMpotgdampotfilename = nullpointer;
    FRMpotgrid = nullpointer;
    FRMpotgrid2D = nullpointer;
    FRMpotgrid3D = nullpointer;
    FRMpotgridres = nullpointer;
    FRMpotgridtype = nullpointer;
    FRMpotlong = nullpointer;
    FRMpotlmaxexp = nullpointer;
    FRMpotderivs = nullpointer;
    FRMpotinput = nullpointer;
    FRMpotmpi = nullpointer;
    FRMpotplane2D = nullpointer;
    FRMpotplaneABC = nullpointer;
    FRMpotresol2D = nullpointer;
    FRMpotresol3D = nullpointer;
    FRMpotsurfpar = nullpointer;
    FRMpotsurftype = nullpointer;
    FRMpotxyz = nullpointer;
    FRMproject = nullpointer;
    FRMImportSGholeden = nullpointer;
    FRMSGholefilename = nullpointer;
    FRMSGholeinput = nullpointer;
    FRMSGholelmaxexp = nullpointer;
    FRMSGholempi = nullpointer;
    FRMtopoboxsize = nullpointer;
    FRMtopofilename = nullpointer;
    FRMtopoguessfilename = nullpointer;
    FRMtopograph = nullpointer;
    FRMtopoinput = nullpointer;
    FRMtopompi = nullpointer;
    FRMtopotype = nullpointer;
    FRMtopoxyz = nullpointer;
    FRMviewers = nullpointer;
    FRMZJinput = nullpointer;
    FRMZJlength = nullpointer;
    FRMZJmpi = nullpointer;
    FRMZJnquad = nullpointer;
    FRMZJrstar = nullpointer;
    FRMZJrstartype = nullpointer;
    FRMZJthreshold = nullpointer;
    FRMZJtype = nullpointer;

    GraphicsMenu = nullpointer;

    HelpMenu = nullpointer;
    HFforceslist = nullpointer;
    HFforcesvalidator = nullpointer;

    languageActionGroup = nullpointer;
    languageMenu = nullpointer;
    LBLatdensmpi = nullpointer;
    LBLdensatoms = nullpointer;
    LBLdensmpi = nullpointer;
    LBLdensuresol = nullpointer;
    LBLdensvresol = nullpointer;
    LBLdensxresol = nullpointer;
    LBLdensyresol = nullpointer;
    LBLdenszresol = nullpointer;
    LBLdensinf3d = nullpointer;
    LBLdensinf2d = nullpointer;
    LBLdenslmaxexp = nullpointer;
    LBLldensminexp = nullpointer;
    LBLdenssup3d = nullpointer;
    LBLdenssup2d = nullpointer;
    LBLdensu = nullpointer;
    LBLdensv = nullpointer;
    LBLdensx = nullpointer;
    LBLdensxformula2D = nullpointer;
    LBLdensy = nullpointer;
    LBLdensyformula2D = nullpointer;
    LBLdensz = nullpointer;
    LBLdenszformula2D = nullpointer;
    LBLdgfilelines = nullpointer;
    LBLdglongthreshold = nullpointer;
    LBLdgmpi = nullpointer;
    LBLdZJchoose = nullpointer;
    LBLdZJImportFile = nullpointer;
    LBLdZJinf = nullpointer;
    LBLdZJinf2d = nullpointer;
    LBLdZJkmax = nullpointer;
    LBLdZJlmax = nullpointer;
    LBLdZJlmin = nullpointer;
    LBLdZJmpi = nullpointer;
    LBLdZJsup = nullpointer;
    LBLdZJsup2d = nullpointer;
    LBLdZJuresol = nullpointer;
    LBLdZJu = nullpointer;
    LBLdZJvresol = nullpointer;
    LBLdZJv = nullpointer;
    LBLdZJxresol = nullpointer;
    LBLdZJx = nullpointer;
    LBLdZJxformula2D = nullpointer;
    LBLdZJyresol = nullpointer;
    LBLdZJy = nullpointer;
    LBLdZJyformula2D = nullpointer;
    LBLdZJzresol = nullpointer;
    LBLdZJz = nullpointer;
    LBLdZJzformula2D = nullpointer;
    LBLefmpi = nullpointer;
    LBLeffilelines = nullpointer;
    LBLeflongthreshold = nullpointer;
    LBLfitthreshold = nullpointer;
    LBLfradthreshold = nullpointer;
    LBLfradatoms = nullpointer;
    LBLfradltr = nullpointer;
    LBLfradltab = nullpointer;
    LBLfradmtab = nullpointer;
    LBLfradrfin = nullpointer;
    LBLfradrini = nullpointer;
    LBLImport = nullpointer;
    LBLlanguage = nullpointer;
    LBLMOchoose = nullpointer;
    LBLMOImportFile = nullpointer;
    LBLMOinf = nullpointer;
    LBLMOinf2d3 = nullpointer;
    LBLMOmpi = nullpointer;
    LBLMOuresol = nullpointer;
    LBLMOvresol = nullpointer;
    LBLMOxresol = nullpointer;
    LBLMOyresol = nullpointer;
    LBLMOzresol = nullpointer;
    LBLMOsup = nullpointer;
    LBLMOsup2d = nullpointer;
    LBLMOu = nullpointer;
    LBLMOv = nullpointer;
    LBLMOx = nullpointer;
    LBLMOxformula2D = nullpointer;
    LBLMOy = nullpointer;
    LBLMOyformula2D = nullpointer;
    LBLMOz = nullpointer;
    LBLMOzformula2D = nullpointer;
    LBLmpicommand = nullpointer;
    LBLmpiflags = nullpointer;
    LBLmrotlmin = nullpointer;
    LBLmrotlmax = nullpointer;
    LBLmrotleft = nullpointer;
    LBLmrotmiddle = nullpointer;
    LBLmrotright = nullpointer;
    LBLmrotorimultatoms = nullpointer;
    LBLpotinf = nullpointer;
    LBLpotinf2d = nullpointer;
    LBLpotmpi = nullpointer;
    LBLpoturesol = nullpointer;
    LBLpotvresol = nullpointer;
    LBLpotxresol = nullpointer;
    LBLpotyresol = nullpointer;
    LBLpotzresol = nullpointer;
    LBLpotres = nullpointer;
    LBLpotsup = nullpointer;
    LBLpotsup2d = nullpointer;
    LBLpotlongthreshold = nullpointer;
    LBLpotu = nullpointer;
    LBLpotv = nullpointer;
    LBLpotx = nullpointer;
    LBLpotxformula2D = nullpointer;
    LBLpoty = nullpointer;
    LBLpotyformula2D = nullpointer;
    LBLpotz = nullpointer;
    LBLpotzformula2D = nullpointer;
    LBLProjectFolder = nullpointer;
    LBLProjectName = nullpointer;
    LBLSGholelocalextrema = nullpointer;
    LBLSGholelocalpower = nullpointer;
    LBLSGholecontour = nullpointer;
    LBLSGholegeomthreshold = nullpointer;
    LBLSGholelmaxexp = nullpointer;
    LBLSGholelongthreshold = nullpointer;
    LBLSGholempi = nullpointer;
    LBLtopoexln = nullpointer;
    LBLtopofdisp = nullpointer;
    LBLtopompi = nullpointer;
    LBLtopoboxb = nullpointer;
    LBLtopoboxg = nullpointer;
    LBLtopoboxl = nullpointer;
    LBLtopoboxt = nullpointer;
    LBLtopoggradthr = nullpointer;
    LBLtopocnvg = nullpointer;
    LBLtopolmaxi = nullpointer;
    LBLtopostepszt = nullpointer;
    LBLZJkmax = nullpointer;
    LBLZJlmax = nullpointer;
    LBLZJmpi = nullpointer;
    LBLZJnquad = nullpointer;
    LBLZJthrdist = nullpointer;
    LBLZJthrmult = nullpointer;

    MOlist = nullpointer;
    MOvalidator = nullpointer;
    mrotorimultlist = nullpointer;
    mrotorimultvalidator = nullpointer;
    myProcess = nullpointer;

    page_project = nullpointer;
    page_atdens = nullpointer;
    page_densgrad = nullpointer;
    page_Efield = nullpointer;
    page_frad = nullpointer;
    page_HFforces = nullpointer;
    page_MED = nullpointer;
    page_MESP = nullpointer;
    page_MO = nullpointer;
    page_projectLayout = nullpointer;
    page_orimult = nullpointer;
    page_SGhole = nullpointer;
    page_TOPO = nullpointer;
    page_ZJdens = nullpointer;
    page_ZJtab = nullpointer;
    plots = nullpointer;

    QDLviewer2D = nullpointer;
    QDLwidget3D = nullpointer;
    QSLdZJxyz = nullpointer;

    RBTatdensD1center = nullpointer;
    RBTatdensD2center = nullpointer;
    RBTatdensDtotal = nullpointer;
    RBTdens2D = nullpointer;
    RBTdens3D = nullpointer;
    RBTdensdeform = nullpointer;
    RBTdensExact = nullpointer;
    RBTdenslrange = nullpointer;
    RBTdensplane = nullpointer;
    RBTdensplaneABC = nullpointer;
    RBTdensplaneXY = nullpointer;
    RBTdensplaneXZ = nullpointer;
    RBTdensplaneYZ = nullpointer;
    RBTdensRep1 = nullpointer;
    RBTdensothersurf = nullpointer;
    RBTdensfulldensity = nullpointer;
    RBTdensrcustom = nullpointer;
    RBTdensrhigh = nullpointer;
    RBTdensrlow = nullpointer;
    RBTdensrmedium = nullpointer;
    RBTdg2D = nullpointer;
    RBTdg3D = nullpointer;
    RBTdgplaneABC = nullpointer;
    RBTdgplaneXY = nullpointer;
    RBTdgplaneXZ = nullpointer;
    RBTdgplaneYZ = nullpointer;
    RBTdZJ2D = nullpointer;
    RBTdZJ3D = nullpointer;
    RBTdZJchooseall = nullpointer;
    RBTdZJchoosek = nullpointer;
    RBTdZJchoosel = nullpointer;
    RBTdZJchooselk = nullpointer;
    RBTdZJchooselkm = nullpointer;
    RBTdZJothersurf = nullpointer;
    RBTdZJplane = nullpointer;
    RBTdZJplaneABC = nullpointer;
    RBTdZJplaneXY = nullpointer;
    RBTdZJplaneXZ = nullpointer;
    RBTdZJplaneYZ = nullpointer;
    RBTdZJrlow = nullpointer;
    RBTdZJrmedium = nullpointer;
    RBTdZJrhigh = nullpointer;
    RBTdZJrcustom = nullpointer;
    RBTef2D = nullpointer;
    RBTef3D = nullpointer;
    RBTefplaneABC = nullpointer;
    RBTefplaneXY = nullpointer;
    RBTefplaneXZ = nullpointer;
    RBTefplaneYZ = nullpointer;
    RBTMO2D = nullpointer;
    RBTMO3D = nullpointer;
    RBTMOothersurf = nullpointer;
    RBTMOplane = nullpointer;
    RBTMOplaneABC = nullpointer;
    RBTMOplaneXY = nullpointer;
    RBTMOplaneXZ = nullpointer;
    RBTMOplaneYZ = nullpointer;
    RBTMOrcustom = nullpointer;
    RBTMOrhigh = nullpointer;
    RBTMOrlow = nullpointer;
    RBTMOrmedium = nullpointer;
    RBTpot2D = nullpointer;
    RBTpot3D = nullpointer;
    RBTpotothersurf = nullpointer;
    RBTpotplane = nullpointer;
    RBTpotplaneABC = nullpointer;
    RBTpotplaneXY = nullpointer;
    RBTpotplaneXZ = nullpointer;
    RBTpotplaneYZ = nullpointer;
    RBTpotrcustom = nullpointer;
    RBTpotrhigh = nullpointer;
    RBTpotrlow = nullpointer;
    RBTpotrmedium = nullpointer;
    RBTtopodensity = nullpointer;
    RBTtopopotential = nullpointer;
    RBTZJacobi = nullpointer;
    RBTZJZernike = nullpointer;
    RBTZJechelon = nullpointer;
    RBTZJrstarabs = nullpointer;
    RBTZJrstarrel = nullpointer;
    rightBox = nullpointer;
    rightBoxLayout = nullpointer;

    SHTdguv = nullpointer;
    SHTdgxyz = nullpointer;
    SHTdZJxyz = nullpointer;
    SHTefuv = nullpointer;
    SHTefxyz = nullpointer;
    SHTfradrlist = nullpointer;
    SHTtopoxyz = nullpointer;
    SHTxyz = nullpointer;
    SPBatdensmpi = nullpointer;
    SPBatdenslmaxexp = nullpointer;
    SPBatdenslmaxdisp = nullpointer;
    SPBdensmpi = nullpointer;
    SPBdensures = nullpointer;
    SPBdensvres = nullpointer;
    SPBdensxres = nullpointer;
    SPBdensyres = nullpointer;
    SPBdenszres = nullpointer;
    SPBdenslmaxexp = nullpointer;
    SPBdenslminexp = nullpointer;
    SPBdglmaxexp = nullpointer;
    SPBdglongthreshold = nullpointer;
    SPBdgmpi = nullpointer;
    SPBdZJkmax = nullpointer;
    SPBdZJlmax = nullpointer;
    SPBdZJlmin = nullpointer;
    SPBdZJmpi = nullpointer;
    SPBdZJures = nullpointer;
    SPBdZJvres = nullpointer;
    SPBdZJxres = nullpointer;
    SPBdZJyres = nullpointer;
    SPBdZJzres = nullpointer;
    SPBefmpi = nullpointer;
    SPBeflmaxexp = nullpointer;
    SPBeflongthreshold = nullpointer;
    SPBfitthreshold = nullpointer;
    SPBfradthreshold = nullpointer;
    SPBfradltab = nullpointer;
    SPBfradmtab = nullpointer;
    SPBpotlmaxexp = nullpointer;
    SPBMOmpi = nullpointer;
    SPBMOures = nullpointer;
    SPBMOvres = nullpointer;
    SPBMOxres = nullpointer;
    SPBMOyres = nullpointer;
    SPBMOzres = nullpointer;
    SPBmrotleft = nullpointer;
    SPBmrotlmax = nullpointer;
    SPBmrotlmin = nullpointer;
    SPBmrotmiddle = nullpointer;
    SPBmrotright = nullpointer;
    SPBpotmpi = nullpointer;
    SPBpotures = nullpointer;
    SPBpotvres = nullpointer;
    SPBpotxres = nullpointer;
    SPBpotyres = nullpointer;
    SPBpotzres = nullpointer;
    SPBpotlongthreshold = nullpointer;
    SPBSGholegeomthreshold = nullpointer;
    SPBSGholelmaxexp = nullpointer;
    SPBSGholelocalextrema = nullpointer;
    SPBSGholelongthreshold = nullpointer;
    SPBSGholempi = nullpointer;
    SPBtopolmaxi = nullpointer;
    SPBtopompi = nullpointer;
    SPBZJkmax = nullpointer;
    SPBZJlmax = nullpointer;
    SPBZJmpi = nullpointer;
    SPBZJnquad = nullpointer;
    SPBZJthrdist = nullpointer;
    SPBZJthrmult = nullpointer;
    splash = nullpointer;

    TAWprincipal = nullpointer;
    textEdit = nullpointer;
    ToolBarFile = nullpointer;
    ToolBarHelp = nullpointer;
    toolBox = nullpointer;
    TXTdensatoms = nullpointer;
    TXTdensdamdenfilename = nullpointer;
    TXTdensplaneA = nullpointer;
    TXTdensplaneB = nullpointer;
    TXTdensplaneC = nullpointer;
    TXTdensuinf = nullpointer;
    TXTdensusup = nullpointer;
    TXTdensvinf = nullpointer;
    TXTdensvsup = nullpointer;
    TXTdensxformula2D = nullpointer;
    TXTdensxinf = nullpointer;
    TXTdensxsup = nullpointer;
    TXTdensyformula2D = nullpointer;
    TXTdensyinf = nullpointer;
    TXTdensysup = nullpointer;
    TXTdenszformula2D = nullpointer;
    TXTdenszinf = nullpointer;
    TXTdenszsup = nullpointer;
    TXTdgatomsforces = nullpointer;
    TXTdgdlt0 = nullpointer;
    TXTdgfilelines = nullpointer;
    TXTdgfilename = nullpointer;
    TXTdgnlinpernuc = nullpointer;
    TXTdgnumpnt = nullpointer;
    TXTdgplaneA = nullpointer;
    TXTdgplaneB = nullpointer;
    TXTdgplaneC = nullpointer;
    TXTdguinf = nullpointer;
    TXTdgusup = nullpointer;
    TXTdguvratio = nullpointer;
    TXTdgvinf = nullpointer;
    TXTdgvsup = nullpointer;
    TXTdgxinf = nullpointer;
    TXTdgxsup = nullpointer;
    TXTdgyinf = nullpointer;
    TXTdgysup = nullpointer;
    TXTdgzinf = nullpointer;
    TXTdgzsup = nullpointer;
    TXTdZJchoose = nullpointer;
    TXTdZJfilename = nullpointer;
    TXTdZJImportfile = nullpointer;
    TXTdZJplaneA = nullpointer;
    TXTdZJplaneB = nullpointer;
    TXTdZJplaneC = nullpointer;
    TXTdZJuinf = nullpointer;
    TXTdZJusup = nullpointer;
    TXTdZJvinf = nullpointer;
    TXTdZJvsup = nullpointer;
    TXTdZJxformula2D = nullpointer;
    TXTdZJxinf = nullpointer;
    TXTdZJxsup = nullpointer;
    TXTdZJyformula2D = nullpointer;
    TXTdZJyinf = nullpointer;
    TXTdZJysup = nullpointer;
    TXTdZJzformula2D = nullpointer;
    TXTdZJzinf = nullpointer;
    TXTdZJzsup = nullpointer;
    TXTefdlt0 = nullpointer;
    TXTeffilelines = nullpointer;
    TXTeffilename = nullpointer;
    TXTefnlinpernuc = nullpointer;
    TXTefnumpnt = nullpointer;
    TXTefplaneA = nullpointer;
    TXTefplaneB = nullpointer;
    TXTefplaneC = nullpointer;
    TXTefuinf = nullpointer;
    TXTefusup = nullpointer;
    TXTefuvratio = nullpointer;
    TXTefvinf = nullpointer;
    TXTefvsup = nullpointer;
    TXTefxinf = nullpointer;
    TXTefxsup = nullpointer;
    TXTefyinf = nullpointer;
    TXTefysup = nullpointer;
    TXTefzinf = nullpointer;
    TXTefzsup = nullpointer;
    TXTfradatoms = nullpointer;
    TXTfraddamfilename = nullpointer;
    TXTfradltr = nullpointer;
    TXTfradrfin = nullpointer;
    TXTfradrini = nullpointer;
    TXTHFforcesatoms = nullpointer;
    TXTHFgdamforcesfilename = nullpointer;
    TXTImport = nullpointer;
    TXTMOchoose = nullpointer;
    TXTMOxformula2D = nullpointer;
    TXTMOyformula2D = nullpointer;
    TXTMOzformula2D = nullpointer;
    TXTMOImportfile = nullpointer;
    TXTMOfilename = nullpointer;
    TXTMOplaneA = nullpointer;
    TXTMOplaneB = nullpointer;
    TXTMOplaneC = nullpointer;
    TXTMOuinf = nullpointer;
    TXTMOusup = nullpointer;
    TXTMOvinf = nullpointer;
    TXTMOvsup = nullpointer;
    TXTMOxinf = nullpointer;
    TXTMOxsup = nullpointer;
    TXTMOyinf = nullpointer;
    TXTMOysup = nullpointer;
    TXTMOzinf = nullpointer;
    TXTMOzsup = nullpointer;
    TXTmpicommand = nullpointer;
    TXTmpiflags = nullpointer;
    TXTmrotorimultatoms = nullpointer;
    TXTmrotorimultfilename = nullpointer;
    TXTpotgdampotfilename = nullpointer;
    TXTpotplaneA = nullpointer;
    TXTpotplaneB = nullpointer;
    TXTpotplaneC = nullpointer;
    TXTpotuinf = nullpointer;
    TXTpotusup = nullpointer;
    TXTpotvinf = nullpointer;
    TXTpotvsup = nullpointer;
    TXTpotxformula2D = nullpointer;
    TXTpotxinf = nullpointer;
    TXTpotxsup = nullpointer;
    TXTpotyformula2D = nullpointer;
    TXTpotyinf = nullpointer;
    TXTpotysup = nullpointer;
    TXTpotzformula2D = nullpointer;
    TXTpotzinf = nullpointer;
    TXTpotzsup = nullpointer;
    TXTProjectFolder = nullpointer;
    TXTProjectName = nullpointer;
    TXTImportSGholeden = nullpointer;
    TXTSGholecontour = nullpointer;
    TXTSGholefilename = nullpointer;
    TXTtopoexln = nullpointer;
    TXTtopofdisp = nullpointer;
    TXTtopoboxb = nullpointer;
    TXTtopoboxg = nullpointer;
    TXTtopoboxl = nullpointer;
    TXTtopoboxt = nullpointer;
    TXTtopocnvg = nullpointer;
    TXTtopofilename = nullpointer;
    TXTtopoggradthr = nullpointer;
    TXTtopoguessfilename = nullpointer;
    TXTtopostepszt = nullpointer;
    TXTZJrstar = nullpointer;

    widgets = nullpointer;
    Wtable1 = nullpointer;
    Wtable2dg = nullpointer;
    Wtable2ef = nullpointer;
    Wtable5 = nullpointer;
    Wtable8 = nullpointer;
    Wtableden = nullpointer;
    Wtabledendg = nullpointer;
    Wtabledenef = nullpointer;
    WtabledZJ = nullpointer;
    Wtablepot = nullpointer;
    Wtabledentopo = nullpointer;

    ZJlist = nullpointer;

    for (int i = 0; i < MAX_ARCHIVOS_RECIENTES; ++i) {
        AccRecentFiles[i] = nullpointer;;
    }
    for (int i = 0; i < NFORTRANPROCS; ++i) {
        BTNstop[i] = nullpointer;
        BTNtexto[i] = nullpointer;
    }
}


//    Creates a suitable input file from options file (*.damproj)
void MainWindow::inputdatafile(const char *suffix, const char *section, 
       const QString &fullFileName, const QString projectDir, const QString projectName)
{
//    string fileinpstr = toString(fullFileName);
    QString fileoutstr = projectDir+projectName+"-"+suffix;
    string v;
    // Opens file .damproj for reading options
    QFile fileinp(fullFileName);
    if (!fileinp.isOpen()){
        fileinp.open(QFile::ReadOnly);
    }
//      Opens file to be filled with suitable input file pointed by suffix
    QFile fileout(fileoutstr);

    if (!fileout.isOpen()){
        fileout.open(QFile::Text | QFile::WriteOnly);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff = ReadSectionOptions(section, &fileinp);
    outfile << buff;
    string vproj = "\"" + projectDir.toStdString() + projectName.toStdString() + "\"" ;
    outfile << vproj.c_str();
#if QT_VERSION < 0x050E00
    outfile << endl;
#else
    outfile << Qt::endl;
#endif
    fileout.close();
}

//    Reads geometry
void MainWindow::readGeometry(int &nats,QVector<double> &x,QVector<double> &y,QVector<double> &z,QVector<int> &ncarga)
{
    rmax = 0.;
    QString suffix;
    if (lslater){
        suffix = ".sxyz";
        QString sxyzfilename = ProjectFolder + FileWithoutExt(ProjectName)+".sxyz";
        if (!(QFile::exists(sxyzfilename))){
            execsgbs2sxyz(sxyzfilename);
        }
        set_natom(read_natom(sxyzfilename));
    }
    else
        suffix = ".ggbs";
    QString filename=ProjectFolder+ProjectName+suffix;
    QFile file(filename);

    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("readGeometry"),tr("File %1 cannot be read")
                .arg(ProjectName+suffix)+
                QString(":\n%1.").arg(file.errorString()));
        nats = 0;
        return;
    }
    QTextStream in(&file);
    QString line = in.readLine();
    if (line.size()==0)
        nats = 0;
    else
        nats = line.toInt();
    x.resize(nats);
    y.resize(nats);
    z.resize(nats);
    ncarga.resize(nats);
    double q;
    double raux;
    double xC = 0.;
    double yC = 0.;
    double zC = 0.;
    double qC = 0.;
    for (int i=0 ; i<nats ; i++){
        QString line = in.readLine();
        if (line.isEmpty()) continue;
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        x[i] = fields.takeFirst().toDouble();
        y[i] = fields.takeFirst().toDouble();
        z[i] = fields.takeFirst().toDouble();
        q = fields.takeLast().toDouble();
        ncarga[i] = (int)q;
        xC += q * x[i];
        yC += q * y[i];
        zC += q * z[i];
        qC += q;
        raux = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
        xC = xC / qC;
        yC = yC / qC;
        zC = zC / qC;
        if (raux > rmax) rmax = raux;
    }
    for (int i = 0 ; i < nats ; i++){
        x[i] = x[i] - xC;
        y[i] = y[i] - yC;
        z[i] = z[i] - zC;
    }
    rmax = sqrt(rmax);
    if (x.size() == 0 || y.size() == 0 || z.size() == 0 ){
        QMessageBox::warning(this, tr("DAMQT"),tr("Error reading geometry in file %1")
                             .arg(ProjectName+suffix));
        return;
    }
}

//    Reads the number of atoms from file fileName
int MainWindow::read_natom(QString fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("read_natom"),tr("File %1 cannot be read").arg(fileName) +
                        QString(":\n%1.").arg(file.errorString()));
        return 0;
    }
    QTextStream in(&file);
    QString line = in.readLine();
    if (line.size()==0){
        return 0;
    }
    else{
        return line.toInt();
    }
}

// Computes dlt for a given density grid resolution
double MainWindow::set_delta(const char * c, double ini, double fin)
{
    double dlt = 1.0;
    if (RBTdensrlow->isChecked()){
        if (RBTdens3D->isChecked())
            dlt=(fin-ini)/64; //2**6+1 (low 3D)
        else
            dlt=(fin-ini)/128; //2**7+1 (low 2D)
    }
    else if (RBTdensrmedium->isChecked()){
        if (RBTdens3D->isChecked())
            dlt=(fin-ini)/128; //2**7+1 (medium 3D)
        else
            dlt=(fin-ini)/256; //2**8+1 (medium 2D)
    }
    else if (RBTdensrhigh->isChecked()){
        if (RBTdens3D->isChecked())
            dlt=(fin-ini)/256; //2**8+1 (high 3D)
        else
            dlt=(fin-ini)/512; //2**9+1 (high 2D)
    }
    else if(RBTdensrcustom->isChecked()){
        if (QString(c).compare(QString("dltx")) == 0)
        dlt=(fin-ini)/SPBdensxres->value();
        else if (QString(c).compare(QString("dlty")) == 0)
        dlt=(fin-ini)/SPBdensyres->value();
        else if (QString(c).compare(QString("dltz")) == 0)
        dlt=(fin-ini)/SPBdenszres->value();
        else if (QString(c).compare(QString("dltu")) == 0)
        dlt=(fin-ini)/SPBdensures->value();
        else if (QString(c).compare(QString("dltv")) == 0)
        dlt=(fin-ini)/SPBdensvres->value();
    }
    return dlt;
}

// Computes dlt for a given orbitals grid resolution
double MainWindow::set_deltaorb(const char * c, double ini, double fin)
{
    double dlt = 1.0;
    if (RBTMOrlow->isChecked()){
        if (RBTMO3D->isChecked()){
            dlt=(fin-ini)/64; //2**6+1 (low 3D)
        }
        else{
            dlt=(fin-ini)/128; //2**7+1 (low 2D)
        }
    }
    else if (RBTMOrmedium->isChecked()){
        if (RBTMO3D->isChecked())
            dlt=(fin-ini)/128; //2**7+1 (medium 3D)
        else
            dlt=(fin-ini)/256; //2**8+1 (medium 2D)
    }
    else if (RBTMOrhigh->isChecked()){
        if (RBTMO3D->isChecked())
            dlt=(fin-ini)/256; //2**8+1 (high 3D)
        else
            dlt=(fin-ini)/512; //2**9+1 (high 2D)
    }
    else if(RBTMOrcustom->isChecked()){
        if (QString(c).compare(QString("dltx")) == 0)
        dlt=(fin-ini)/SPBMOxres->value();
        else if (QString(c).compare(QString("dlty")) == 0)
        dlt=(fin-ini)/SPBMOyres->value();
        else if (QString(c).compare(QString("dltz")) == 0)
        dlt=(fin-ini)/SPBMOzres->value();
        else if (QString(c).compare(QString("dltu")) == 0)
        dlt=(fin-ini)/SPBMOures->value();
        else if (QString(c).compare(QString("dltv")) == 0)
        dlt=(fin-ini)/SPBMOvres->value();
    }
    return dlt;
}

// Computes dlt for a given potential grid resolution
double MainWindow::set_deltapot(const char * c, double ini, double fin)
{
    double dlt = 1.0;
    if (RBTpotrlow->isChecked()){
        if (RBTpot3D->isChecked()){
            dlt=(fin-ini)/64; //2**6+1 (low 3D)
        }
        else{
            dlt=(fin-ini)/128; //2**7+1 (low 2D)
        }
    }
    else if (RBTpotrmedium->isChecked()){
        if (RBTpot3D->isChecked())
            dlt=(fin-ini)/128; //2**7+1 (medium 3D)
        else
            dlt=(fin-ini)/256; //2**8+1 (medium 2D)
    }
    else if (RBTpotrhigh->isChecked()){
        if (RBTpot3D->isChecked())
            dlt=(fin-ini)/256; //2**8+1 (high 3D)
        else
            dlt=(fin-ini)/512; //2**9+1 (high 2D)
    }
    else if(RBTpotrcustom->isChecked()){
        if (QString(c).compare(QString("dltx")) == 0)
        dlt=(fin-ini)/SPBpotxres->value();
        else if (QString(c).compare(QString("dlty")) == 0)
        dlt=(fin-ini)/SPBpotyres->value();
        else if (QString(c).compare(QString("dltz")) == 0)
        dlt=(fin-ini)/SPBpotzres->value();
        else if (QString(c).compare(QString("dltu")) == 0)
        dlt=(fin-ini)/SPBpotures->value();
        else if (QString(c).compare(QString("dltv")) == 0)
        dlt=(fin-ini)/SPBpotvres->value();
    }
    return dlt;
}

// Computes dlt for a given Jacobi-Zernike density grid resolution
double MainWindow::set_deltaZJ(const char * c, double ini, double fin)
{
    double dlt = 1.0;
    if (RBTdZJrlow->isChecked()){
        if (RBTdZJ3D->isChecked())
            dlt=(fin-ini)/64; //2**6+1 (low 3D)
        else
            dlt=(fin-ini)/128; //2**7+1 (low 2D)
    }
    else if (RBTdZJrmedium->isChecked()){
        if (RBTdZJ3D->isChecked())
            dlt=(fin-ini)/128; //2**7+1 (medium 3D)
        else
            dlt=(fin-ini)/256; //2**8+1 (medium 2D)
    }
    else if (RBTdZJrhigh->isChecked()){
        if (RBTdZJ3D->isChecked())
            dlt=(fin-ini)/256; //2**8+1 (high 3D)
        else
            dlt=(fin-ini)/512; //2**9+1 (high 2D)
    }
    else if(RBTdZJrcustom->isChecked()){
        if (QString(c).compare(QString("dltx")) == 0)
        dlt=(fin-ini)/SPBdZJxres->value();
        else if (QString(c).compare(QString("dlty")) == 0)
        dlt=(fin-ini)/SPBdZJyres->value();
        else if (QString(c).compare(QString("dltz")) == 0)
        dlt=(fin-ini)/SPBdZJzres->value();
        else if (QString(c).compare(QString("dltu")) == 0)
        dlt=(fin-ini)/SPBdZJures->value();
        else if (QString(c).compare(QString("dltv")) == 0)
        dlt=(fin-ini)/SPBdZJvres->value();
    }
    return dlt;
}

//    sets variable natom
void MainWindow::set_natom(int i)
{
    MainWindow::natom = i;
}

void MainWindow::start(){
FRMlanguage->close();
}


/*******************************************************************************************************/
/******************************** FILE NAME, PATH AND EXTENSION HANDLING *******************************/
/*******************************************************************************************************/


/* Returns the file extension */
QString MainWindow::Extension(const QString &fullFileName)
{
    return QFileInfo(fullFileName).suffix();
}

/* Returns the file name without any extension */
QString MainWindow::FileWithoutExt(const QString &fullFileName)
{
    return QFileInfo(fullFileName).completeBaseName();
}

/* Returns the file name without path */
QString MainWindow::FileWithoutPath(const QString &fullFileName)
{

    return QFileInfo(fullFileName).fileName();
}

/* Returns the path of a file */
QString MainWindow::Path(const QString &fullFileName)
{
    return QFileInfo(fullFileName).path();
}

/*******************************************************************************************************/
/*********************   2D AND 3D ViewerS SLOTS AND FUNCTIONS     *************************************/
/*******************************************************************************************************/

bool MainWindow::end_Viewer2DDialog(){
    if (plots && !plots->isEmpty()){
        QMessageBox msgBox;
        msgBox.setText(tr("2D Viewer: Delete Confirmation"));
        msgBox.setInformativeText(tr("Ending application will delete open 2D viewers. Do you want to exit?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setIcon(QMessageBox::Question);
        int ret = msgBox.exec();
        if (ret == QMessageBox::No)
            return false;
        for (int i = 0 ; i < plots->length(); i++){
            delete plots->at(i);
        }
        plots->clear();  
    }
    update_dockright();
    return true;
}

bool MainWindow::end_Viewer3DDialog(){
    if (widgets && !widgets->isEmpty()){
        QMessageBox msgBox;
        msgBox.setText(tr("3D Viewer: Delete Confirmation"));
        msgBox.setInformativeText(tr("Ending application will delete open 3D viewers. Do you want to exit?"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Yes);
        msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
        msgBox.setButtonText(QMessageBox::No, tr("No"));
        msgBox.setIcon(QMessageBox::Question);
        int ret = msgBox.exec();
        if (ret == QMessageBox::No)
            return false;

        qDeleteAll(widgets->begin(), widgets->end());
        widgets->clear();
    }
    update_dockright();
    return true;
}

void MainWindow::addviewer2D(){
    QString *name = new QString(tr("%1").arg(plotsknt++));
    Viewer2D *viewer = new Viewer2D(name);
    viewer->set_ProjectFolder(ProjectFolder);
    viewer->set_ProjectName(ProjectName);
    connect(viewer, SIGNAL(moveToTop(int)), this,SLOT(moveviewertotop(int)), Qt::UniqueConnection);
    plots->append(viewer);
    plots->last()->set_position(QPoint(75,75)*(widgets->length()+plots->length()-1));
    this->moveviewertotop(plotsknt-1);
    topindex = plots->length()-1;
    update_dockright();
}

void MainWindow::addglWidget(){
    QString *name = new QString(tr("%1").arg(widgetsknt++));
    glWidget *widget = new glWidget(name,this);
    connect(widget, SIGNAL(moveToTop(int)), this,SLOT(movewidgettotop(int)), Qt::UniqueConnection);
    widgets->append(widget);
    widgets->last()->set_position(QPoint(75,75)*(widgets->length()+plots->length()-1));
    widgets->last()->set_ProjectFolder(ProjectFolder);
    widgets->last()->set_ProjectName(ProjectName);
    this->movewidgettotop(widgetsknt-1);
    topindex = widgets->length()+plots->length()-1;
    update_dockright();
}

void MainWindow::exit(){

    this->close();
}

void MainWindow::deleteplot(int i){
    int number = plots->at(i)->getviewernumber();
    QMessageBox msgBox;
    msgBox.setInformativeText(QString(tr("Do you want to remove Plot %1")).arg(number)+"?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);
    msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
    msgBox.setButtonText(QMessageBox::No, tr("No"));
    msgBox.setIcon(QMessageBox::Question);
    int ret = msgBox.exec();
    if (ret == QMessageBox::No){
        plots->at(i)->raise_mainwindow();
        return;
    }
    delete plots->at(i);
    plots->removeAt(i);
    if (i == topindex)
        topindex = -1;
    update_dockright();
}

void MainWindow::deletewidget(int i){
    int number = widgets->at(i)->getwindownumber();
    QMessageBox msgBox;
    msgBox.setInformativeText(QString(tr("Do you want to remove 3D viewer %1")).arg(number)+"?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);
    msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
    msgBox.setButtonText(QMessageBox::No, tr("No"));
    msgBox.setIcon(QMessageBox::Question);
    int ret = msgBox.exec();
    if (ret == QMessageBox::No){
//        this->movewidgettotop(0);
        widgets->at(i)->raise_mainwindow();
        return;
    }
    delete widgets->at(i);
    widgets->removeAt(i);
    if (i == topindex-plots->length())
        topindex = -1;
    update_dockright();
}

void MainWindow::moveviewertotop(int num){
    for (int i = 0 ; i < plots->length() ; i++){
        if (num == plots->at(i)->getviewernumber()){
            topindex = i;
            plots->at(topindex)->raise_mainwindow();
            return;
        }
    }
}

void MainWindow::movewidgettotop(int num){
    for (int i = 0 ; i < widgets->length() ; i++){
        if (num == widgets->at(i)->getwindownumber()){
            topindex = i + plots->length();
            widgets->at(i)->raise_mainwindow();
            return;
        }
    }
}

void MainWindow::raiseplot(int i){
    plots->at(i)->raise_mainwindow();
    topindex = i;
}

void MainWindow::raisewidget(int i){
    widgets->at(i)->raise_mainwindow();
    topindex = i + plots->length();
}

void MainWindow::showplot(int i){
    if (plots->at(i)->isvisible()){
        plots->at(i)->set_visible(false);
        BTNshowplotsslist.at(i)->setText(tr("Show"));
        QPoint position = plots->at(i)->get_position();
        QSize size = plots->at(i)->get_size();
        plots->at(i)->set_size(size);
        plots->at(i)->set_position(position);
        if (topindex < 0 )
            return;
        if (i != topindex){
            if (topindex < plots->length())
                plots->at(topindex)->raise_mainwindow();
            else if (topindex < plots->length() + widgets->length()) {
                widgets->at(topindex-plots->length())->raise_mainwindow();
            }
        }
    }
    else{
        plots->at(i)->set_visible(true);
        moveviewertotop(plots->at(i)->getviewernumber());
        BTNshowplotsslist.at(i)->setText(tr("Hide"));
    }
}

void MainWindow::showwidget(int i){
    if (widgets->at(i)->isvisible()){
        widgets->at(i)->set_visible(false);
        BTNshowwidgetslist.at(i)->setText(tr("Show"));
        QPoint position = widgets->at(i)->get_position();
        QSize size = widgets->at(i)->get_size();
        widgets->at(i)->set_size(size);
        widgets->at(i)->set_position(position);
        if (topindex < 0 )
            return;
        if (i != topindex-plots->length()){
            if (topindex < plots->length())
                plots->at(topindex)->raise_mainwindow();
            else if (topindex < plots->length() + widgets->length()) {
                widgets->at(topindex-plots->length())->raise_mainwindow();
            }
        }
    }
    else{
        widgets->at(i)->set_visible(true);
        movewidgettotop(widgets->at(i)->getwindownumber());
        BTNshowwidgetslist.at(i)->setText(tr("Hide"));
    }
}

void MainWindow::update_dockright(){
    if (BTNnewplot){
        delete BTNnewplot;
        BTNnewplot = nullpointer;
    }
    if (FRMplots){
        delete FRMplots;
        FRMplots = nullpointer;
    }
    if (QDLviewer2D){
        delete QDLviewer2D;
        QDLviewer2D = nullpointer;
    }

    if (BTNnewwidget){
        delete BTNnewwidget;
        BTNnewwidget = nullpointer;
    }
    if (FRMviewers){
        delete FRMviewers;
        FRMviewers = nullpointer;
    }
    if (QDLwidget3D){
        delete QDLwidget3D;
        QDLwidget3D = nullpointer;
    }

    if ( BTNraiseviewers){
        delete  BTNraiseviewers;
         BTNraiseviewers = nullpointer;
    }

    if (dockright){
        delete dockright;
        dockright = nullpointer;
    }
    dockright = new QDockWidget(tr(""),this);
    dockright->setAllowedAreas(Qt::RightDockWidgetArea);
    dockright->resize(QSize(500, this->height()));
    dockright->setFeatures(QDockWidget::DockWidgetMovable);
    dockright->setFeatures(QDockWidget::DockWidgetFloatable);
    QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(1000);
    sizePolicy1.setHeightForWidth(dockright->sizePolicy().hasHeightForWidth());
    dockright->setSizePolicy(sizePolicy1);
    addDockWidget(Qt::RightDockWidgetArea,dockright);

    update_menu_viewer2D();
    update_menu_viewer3D();

    BTNraiseviewers = new QPushButton(tr("Raise all viewers"));
    BTNraiseviewers->setToolTip(tr("Move viewers to front"));
    BTNraiseviewers->setStyleSheet("QPushButton {background-color: darkGreen; color: white;}");
    connect(BTNraiseviewers, SIGNAL(clicked()), this, SLOT(updatewindowsoverlay()));

    QVBoxLayout *layout1 = new QVBoxLayout();
    layout1->addWidget(BTNraiseviewers);

    QWidget *btnwidget = new QWidget();
    btnwidget->setLayout(layout1);

    dockwidget = new QWidget();
    QVBoxLayout *layout = new QVBoxLayout(dockwidget);
    layout->addWidget(btnwidget);
    layout->addWidget(QDLviewer2D);
    layout->addWidget(QDLwidget3D);
    layout->addStretch();

    dockright->setWidget(dockwidget);
    updatewindowsoverlay();
}


void MainWindow::update_menu_viewer2D()
{
    for (int i = 0 ; i < connections2D.size() ; i++){
        QObject::disconnect(connections2D.at(i));
    }
    connections2D.clear();
    if (BTNnewplot){
        delete BTNnewplot;
        BTNnewplot = nullpointer;
    }
    QDLviewer2D = new ViewerDialog();
    QDLviewer2D->setWindowTitle("2D Viewer: choose an option");
    QDLviewer2D->resize(300,80);
    FRMplots = new QGroupBox();
    FRMplots->setTitle("Available plots");
    if (plots->count() > 0)
        FRMplots->setVisible(true);
    else
        FRMplots->setVisible(false);
    QGridLayout *layout1 = new QGridLayout(FRMplots);
    QSignalMapper* deletesignalMapper = new QSignalMapper (this) ;
    QSignalMapper* raisesignalMapper = new QSignalMapper (this) ;
    QSignalMapper* showsignalMapper = new QSignalMapper (this) ;
    BTNshowplotsslist.clear();
    BTNraiseplotslist.clear();
    for (int i = 0 ; i < plots->count() ; i++){
        QLabel *LBLplot = new QLabel();
        LBLplot->setText(plots->at(i)->get_viewername());
        QPushButton *BTNdeleteplot = new QPushButton();
        BTNdeleteplot->setText(tr("Delete"));
        connections2D << connect(BTNdeleteplot, SIGNAL(clicked()), deletesignalMapper,SLOT(map()), Qt::UniqueConnection);
        deletesignalMapper -> setMapping(BTNdeleteplot,i);
        QPushButton *BTNraise = new QPushButton(tr("Raise"));
        BTNraiseplotslist.append(BTNraise);
        connections2D << connect(BTNraise, SIGNAL(clicked()), raisesignalMapper,SLOT(map()), Qt::UniqueConnection);
        raisesignalMapper -> setMapping(BTNraise,i);
        QPushButton *BTNshow = new QPushButton();
        if (plots->at(i)->isvisible())
            BTNshow->setText(tr("Hide"));
        else
            BTNshow->setText(tr("Show"));
        BTNshowplotsslist.append(BTNshow);
        connections2D << connect(plots->at(i), SIGNAL(hideplotter()), BTNshow,SIGNAL(clicked()), Qt::UniqueConnection);
        connections2D << connect(BTNshow, SIGNAL(clicked()), showsignalMapper,SLOT(map()), Qt::UniqueConnection);
        showsignalMapper -> setMapping(BTNshow,i);
        layout1->addWidget(LBLplot,i,0);
        layout1->addWidget(BTNraise,i,1);
        layout1->addWidget(BTNshow,i,2);
        layout1->addWidget(BTNdeleteplot,i,3);
    }
    connections2D << connect (deletesignalMapper, SIGNAL(mapped(int)), this,SLOT(deleteplot(int)), Qt::UniqueConnection) ;
    connections2D << connect (showsignalMapper, SIGNAL(mapped(int)), this, SLOT(showplot(int)), Qt::UniqueConnection) ;
    connections2D << connect (raisesignalMapper, SIGNAL(mapped(int)), this, SLOT(raiseplot(int)), Qt::UniqueConnection) ;

    BTNnewplot = new QPushButton(tr("New 2D Plotter"));
    BTNnewplot->setToolTip(tr("Creates a new window for 2D plotting"));
    connections2D << connect(BTNnewplot, SIGNAL(clicked()), this, SLOT(addviewer2D()));

    QLabel *label2D = new QLabel(tr("2D plotters"));
    label2D->setStyleSheet("QLabel { color : blue; }");

    QHBoxLayout *layout2 = new QHBoxLayout();
    layout2->addWidget(label2D,Qt::AlignCenter);

    QVBoxLayout *layout3=new QVBoxLayout(QDLviewer2D);
    layout3->addLayout(layout2);
    layout3->addWidget(FRMplots);
    layout3->addWidget(BTNnewplot);
    layout3->addStretch();
}

void MainWindow::update_menu_viewer3D()
{
    for (int i = 0 ; i < connections3D.size() ; i++){
        QObject::disconnect(connections3D.at(i));
    }
    connections3D.clear();
    if (BTNnewwidget){
        delete BTNnewwidget;
        BTNnewwidget = nullpointer;
    }
    QDLwidget3D = new ViewerDialog();
    QDLwidget3D->setMinimumSize(250,80);
    FRMviewers = new QGroupBox();
    FRMviewers->setTitle("Available 3D viewers");
    if (widgets->count() > 0)
        FRMviewers->setVisible(true);
    else
        FRMviewers->setVisible(false);
    QGridLayout *layout1 = new QGridLayout(FRMviewers);
    QSignalMapper* deletesignalMapper = new QSignalMapper (this) ;
    QSignalMapper* raisesignalMapper = new QSignalMapper (this) ;
    QSignalMapper* showsignalMapper = new QSignalMapper (this) ;
    BTNshowwidgetslist.clear();
    BTNraisewidgetslist.clear();
    for (int i = 0 ; i < widgets->count() ; i++){
        QLabel *LBLwidget = new QLabel();
        LBLwidget->setText(widgets->at(i)->getWindowName());
        QPushButton *BTNdeletewidget = new QPushButton();
        BTNdeletewidget->setText(tr("Delete"));
        connections3D << connect(BTNdeletewidget, SIGNAL(clicked()), deletesignalMapper,SLOT(map()), Qt::UniqueConnection);
        deletesignalMapper -> setMapping(BTNdeletewidget,i);
        QPushButton *BTNraise = new QPushButton(tr("Raise"));
        BTNraisewidgetslist.append(BTNraise);
        connections3D << connect(BTNraise, SIGNAL(clicked()), raisesignalMapper,SLOT(map()), Qt::UniqueConnection);
        raisesignalMapper -> setMapping(BTNraise,i);
        QPushButton *BTNshow = new QPushButton();
        if (widgets->at(i)->isvisible())
            BTNshow->setText(tr("Hide"));
        else
            BTNshow->setText(tr("Show"));
        BTNshowwidgetslist.append(BTNshow);
        connections3D << connect(widgets->at(i), SIGNAL(hideviewer()), BTNshow,SIGNAL(clicked()), Qt::UniqueConnection);
        connections3D << connect(BTNshow, SIGNAL(clicked()), showsignalMapper,SLOT(map()), Qt::UniqueConnection);
        showsignalMapper -> setMapping(BTNshow,i);
        layout1->addWidget(LBLwidget,i,0);
        layout1->addWidget(BTNraise,i,1);
        layout1->addWidget(BTNshow,i,2);
        layout1->addWidget(BTNdeletewidget,i,3);
    }
    connections3D << connect (deletesignalMapper, SIGNAL(mapped(int)), this,SLOT(deletewidget(int)), Qt::UniqueConnection) ;
    connections3D << connect (showsignalMapper, SIGNAL(mapped(int)), this, SLOT(showwidget(int)), Qt::UniqueConnection) ;
    connections3D << connect (raisesignalMapper, SIGNAL(mapped(int)), this, SLOT(raisewidget(int)), Qt::UniqueConnection) ;

    BTNnewwidget = new QPushButton(tr("New 3D Viewer"));
    BTNnewwidget->setToolTip(tr("Creates a new window for 3D display"));
    connections3D << connect(BTNnewwidget, SIGNAL(clicked()), this, SLOT(addglWidget()));

    QLabel *label3D = new QLabel(tr("3D viewers"));
    label3D->setStyleSheet("QLabel { color : red; }");

    QHBoxLayout *layout2 = new QHBoxLayout();
    layout2->addWidget(label3D,Qt::AlignCenter);

    QVBoxLayout *layout3=new QVBoxLayout(QDLwidget3D);
    layout3->addLayout(layout2);
    layout3->addWidget(FRMviewers);
    layout3->addWidget(BTNnewwidget);
    layout3->addStretch();
}

void MainWindow::updatewindowsoverlay(){
    for (int i = 0 ; i < plots->length() ; i++){
        plots->at(i)->raise_mainwindow();
    }
    for (int i = 0 ; i < widgets->length() ; i++){
        widgets->at(i)->raise_mainwindow();
    }
    if (topindex >= 0){
        if (topindex < plots->length()){
            plots->at(topindex)->raise_mainwindow();
        }
        else if (topindex < widgets->length()+plots->length()){
                widgets->at(topindex-plots->length())->raise_mainwindow();
        }
    }
}

QString MainWindow::get_execName(QString processname, QString subdir){
#if defined(Q_WS_WIN) || defined(Q_OS_WIN)
    QString execName = processname;
    if (!QFileInfo::exists(execName)){
        execName = QCoreApplication::applicationDirPath()+"/"+processname;
    }
#else
    QString execName = QCoreApplication::applicationDirPath()+"/"+processname;
    if (!QFileInfo::exists(execName)){
        execName = processname;
    }
#endif
    if (!QFileInfo::exists(execName))
        execName = QCoreApplication::applicationDirPath()+"/../"+subdir+"/"+processname;
    if (!QFileInfo::exists(execName)){
        QString message1, message2, message3;
        QString direc = QString(QCoreApplication::applicationDirPath());
        direc.truncate(direc.lastIndexOf(QChar('/')));
        message1 = QString(tr("Executable file %1 does not exist\n\n").arg(processname));
        message2 = QString(tr("Check that the program is installed in any of the following directories: \n\n %1 \n %2 \n\n")
                        .arg(QCoreApplication::applicationDirPath()+"/")
                        .arg(direc+"/"+subdir+"/"));
        message3 = QString(tr("or in any other directory available in your $PATH"));
        int messagelen = qMax(qMax(message1.length(),message2.length()),message3.length());
        QMessageBox msg;
        msg.setText(message1+message2+message3);
        msg.setIcon(QMessageBox::Critical);
        QSpacerItem* horizontalSpacer = new QSpacerItem(messagelen * 4, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
        QGridLayout* layout = (QGridLayout*)msg.layout();
        layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
        msg.exec();
        return QString("");
    }
    return execName;
}

QString MainWindow::get_python(){
    QString execName = "python";
    QProcess *process = new QProcess(this);
    process->start(execName+" -h");
    if (process->error() == QProcess::FailedToStart){
        execName = "python2";
        process->start(execName+" -h");
        if (process->error() == QProcess::FailedToStart){
            execName = "python3";
            process->start(execName+" -h");
            if (process->error() == QProcess::FailedToStart){
                execName = "";
            }
        }
    }
    if (execName.isEmpty())
        QMessageBox::warning(this, tr("get_python"),tr("Python not found"));
    return execName;
}


/*******************************************************************************************************/
/********************************  Class ViewerDialog  implementation  *******************************/
/*******************************************************************************************************/

ViewerDialog::ViewerDialog()
{

}

ViewerDialog::~ViewerDialog(){

}

void ViewerDialog::reject(){
    emit closed();
}

/*******************************************************************************************************/
/********************************  Class mainmenu  implementation  *******************************/
/*******************************************************************************************************/

void mainmenu::resizeEvent(QResizeEvent *event){
    size = event->size();
}
