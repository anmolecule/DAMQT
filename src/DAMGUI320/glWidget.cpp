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
//	Implementation of class glWidget
//
//	File:   glWidget.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QPoint>
#include <QVBoxLayout>
#include <QCloseEvent>
#include <QColorDialog>
#include <QSignalMapper>

#include "glWidget.h"
#include "measures.h"

glWidget::glWidget(QString *name, QWidget *parent) : QWidget(parent)
{
    initpointers();
    initQDLpointers();
    connections.clear();
    molecule_names.clear();
    msrs_connections.clear();
    windownumber = (*name).toInt();
    setAttribute(Qt::WA_DeleteOnClose);
    setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);  
    setWindowName(QString("Display no. %1").arg(*name));
    BTNshowlist.clear();

    ambientcolor = QColor(128,128,128);
    animating = false;
    angstrom = false;
    ballradius = 0.2;
    bkgcolor = QColor(0, 0, 0, 255);
    cluster_exists = false;
    clustername = "cluster";
    cylradius = 0.05;
    deltaAngles = 4.f;
    delay = 100;
    displayEPIC = true;
    disthressq = pow((INIT_BOND_THRESHOLD * ANGSTROM_TO_BOHR),2);
    dltinterval = (MAX_INTERVAL-MIN_INTERVAL)/float(INTERVAL_SCALE);
    elem = new Elements();
    energycolor = QColor(255, 172, 0, 255);
    energyfont = QFont("Helvetica", 18, QFont::Bold);
    energyprecision = 4;
    fov = PERSPECTIVE_ANGLE;
    framesfile = "";
    fullmolname = QString("");
    gdockrewidth = 1;
    guestfromcanvas = true;
    hartree = false;
    interpoints = 10;
    interval = MAX_INTERVAL - dltinterval * INTERVAL_INI;
    lightcolor = QColor(230,230,230);
    lightpower = 20.f;
    lightPosition = QVector3D(0,0,50);
    mespimizerfile = "";
    mespimirecordfilename = "";
    molecules = new QList<molecule*>();
    molname = QString("");
    molpath = QString("");
    ninterpol = NINTERPOL;
    optimvisible = false;
    onlyselcp = false;
    position = QPoint(0,0);
    ProjectFolder = "";
    recordfilename = "";
    recordoptim = false;
    scaleradii = true;
    size = QSize(1000,500);
    specularcolor = QColor(128,128,128);
    specularindex = 5.;
    speed = INTERVAL_INI;
    stepwheel = 0.1;
    templateindex = 1;
    xyz.clear();
    visible = true;
    zFar = ZFAR;
    zNear = ZNEAR;
    znuc.clear();

    create_window();

    mainWin = new MainWindow3DViewer();

    gdock = new MyDockWidget(mainWin);

    windowArea = new QScrollArea();
    windowArea->setWidget(QWidget::createWindowContainer(window, this));
    windowArea->setWidgetResizable(true);
    windowArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    windowArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    windowArea->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    windowArea->setMinimumSize(500, 300);

    gdock->setAllowedAreas(Qt::RightDockWidgetArea);
    gdock->resize(QSize(500, this->height()));
    gdock->setFeatures(QDockWidget::DockWidgetMovable);
    gdock->setFeatures(QDockWidget::DockWidgetFloatable);
    QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(1000);
    sizePolicy1.setHeightForWidth(gdock->sizePolicy().hasHeightForWidth());
    gdock->setSizePolicy(sizePolicy1);

    createWindowDialog(gdock);
    scrollArea = new myScrollArea();
    scrollArea->setWidget(QDLwindow);
    connections << connect(gdock,SIGNAL(resizedialog(QSize)),QDLwindow,SLOT(renewwidth(QSize)));
    gdock->setWidget(scrollArea);

//      mainWin settings
    mainWin->setCentralWidget(windowArea);
    mainWin->setWindowTitle(windowname);
    mainWin->addDockWidget(Qt::RightDockWidgetArea,gdock);
    mainWin->resize(size);
    mainWin->move(position);
    connections << connect(mainWin,SIGNAL(hide_viewer()),this,SIGNAL(hideviewer()));

    mainWin->show();

}

glWidget::~glWidget(){
    if (window){
        window->makeCurrent();
        delete window;
    }
    if (QDLwindow){
        delete QDLwindow;
        QDLwindow = nullpointer;
    }
    if (scrollArea){
        delete scrollArea;
        scrollArea = nullpointer;
    }

    if (windowArea){
        delete windowArea;
        windowArea = nullpointer;
    }
    if (gdock){
        delete gdock;
        gdock = nullpointer;
    }

    if (mainWin){
        delete mainWin;
        mainWin = nullpointer;
    }


    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
    if (molecules){
        qDeleteAll(molecules->begin(), molecules->end());
        molecules->clear();
        delete molecules;
        molecules = nullpointer;
    }
    if (msrs){
        for (int i = 0 ; i < msrs_connections.size() ; i++){
            QObject::disconnect(msrs_connections.at(i));
        }
        msrs_connections.clear();
        delete msrs;
        msrs = nullpointer;
    }
}

void glWidget::closeEvent(QCloseEvent *event){
    event->ignore();
    QMessageBox msgBox;
    msgBox.setText(tr("Delete Confirmation"));
    msgBox.setInformativeText(tr("Are you sure you want to close this window?\nViewer and ancillary windows will be deleted"));
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);
    msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
    msgBox.setButtonText(QMessageBox::No, tr("No"));
    msgBox.setIcon(QMessageBox::Question);
    int ret = msgBox.exec();
    if (ret == QMessageBox::No)
        return;
    this->deleteLater();
    event->accept();
}


//  -----------------------------------------------------    SLOTS    ------------------------------------------------
//
//          Dialog window
//
//  ------------------------------------------------------------------------------------------------------------------


void glWidget::createWindowDialog(QDockWidget *gdock){
    gdock->resize(QSize(500, this->height()));
    QDLwindow = new moleculeDialog(this);
    QDLwindow->setWindowTitle(getWindowName()+QString(tr(" Main Menu")));
    QDLwindow->setMinimumWidth(410);
    QDLwindow->resize(410,50);
    QDLwindow->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    connections << connect(QDLwindow,SIGNAL(closedialog()),this,SLOT(close()));
    connections << connect(mainWin,SIGNAL(hide_viewer()),this,SIGNAL(hideviewer()));
    connections << connect(gdock,SIGNAL(resizedialog(QSize)),QDLwindow,SLOT(renewwidth(QSize)));

    //        Add molecule button

    QPushButton *BTNaddmolecule = new QPushButton(tr("Add molecule"));
    connections << connect(BTNaddmolecule, SIGNAL(clicked()), this, SLOT(addmoleculeTowindow()));

    //        Molecules

    create_molecules_menu_and_layouts();    // Molecules menu and layouts

    //        Measure distances and angles
    BTNmeasures = new QPushButton(tr("Geometry measures"));
    connections << connect(BTNmeasures, SIGNAL(clicked()), this, SLOT(BTNmeasures_clicked()));

    //        Rotation
    create_rotation_menu();      // Rotation menu

    //        Translation
    create_translation_menu();      // Translation menu

    //        Axes
    create_axes_menu();             // Axes menu

    //        Capture

    create_capture_menu();          // Capture menu

    //        Lights

    create_lights_menu();        // Lights menu

    //        Balls and cylinders

    create_balls_cyl_menu();     // Balls and cylinders menu

    //        Viewport

    create_viewport_menu();      // Viewport menu


    //        Optimize cluster (MESPIMIZE)

    create_optimize_cluster_menu();     // Optimize cluster menu (MESPIMIZE)

    //        Save geometry

    create_save_geometry_menu();        // Save geometry menu
    create_save_geometry_layouts();     // Save geometry layouts

    //        Save/Retrieve settings

    create_save_retrieve_menu();        // Save/Retrieve settings menu
    create_save_retrieve_layouts();     // Save/Retrieve settings layouts

    //        Main layout

    QVBoxLayout *layout = new QVBoxLayout(QDLwindow);
    layout->addWidget(BTNaddmolecule);
    layout->addLayout(layoutmolecules);
    layout->addWidget(BTNmeasures);
    layout->addWidget(BTNrotation);
    layout->addWidget(FRMrotation);
    layout->addWidget(BTNtranslation);
    layout->addWidget(FRMtranslation);
    layout->addWidget(BTNaddaxes);
    layout->addWidget(FRMaxes);
    layout->addWidget(BTNcapture);
    layout->addWidget(FRMcapture);
    layout->addWidget(BTNlights);
    layout->addWidget(FRMlights);
    layout->addWidget(BTNballcyl);
    layout->addWidget(FRMballcyl);
    layout->addWidget(BTNviewport);
    layout->addWidget(FRMviewport);
    layout->addWidget(BTNoptimizeCluster);
    layout->addWidget(FRMoptimizeCluster);
    layout->addWidget(BTNsaveGeometry);
    layout->addWidget(FRMsaveGeometry);
    layout->addWidget(BTNsettings);
    layout->addWidget(FRMsettings);
    layout->addStretch();

    QDLwindow->show();
}

// Molecules menu and layouts

void glWidget::create_molecules_menu_and_layouts(){

    activateboxes.clear();

    layoutmolecules = new QGridLayout();
    QSignalMapper* deletesignalMapper = new QSignalMapper (this) ;
    QSignalMapper* showsignalMapper = new QSignalMapper (this) ;
    BTNshowlist.clear();
    for (int i = 0 ; i < molecules->count() ; i++){
        molecules->at(i)->setparentposition(this->pos());
        QLabel *LBLmol = new QLabel();
        LBLmol->setText(molecules->at(i)->getname());
        QCheckBox *CHKmol = new QCheckBox();
        CHKmol->setText(tr("active"));
        if (molecules->at(i)->isactive()){
            CHKmol->setChecked(true);
        }
        else{
            CHKmol->setChecked(false);
        }
        connections << connect(CHKmol, SIGNAL(stateChanged(int)), molecules->at(i), SLOT(toggleactive()));
        activateboxes << CHKmol;
        QPushButton *BTNmol = new QPushButton();
        BTNmol->setText(tr("Edit"));
        connections << connect(BTNmol, SIGNAL(clicked()), molecules->at(i), SLOT(editmolecule()));
        QPushButton *BTNdeletemol = new QPushButton();
        BTNdeletemol->setText(tr("Delete"));
        connections << connect(BTNdeletemol, SIGNAL(clicked()), deletesignalMapper, SLOT(map()), Qt::UniqueConnection);
        deletesignalMapper -> setMapping(BTNdeletemol,i);
        QPushButton *BTNshow = new QPushButton();
        if (molecules->at(i)->isvisible())
            BTNshow->setText(tr("Hide"));
        else
            BTNshow->setText(tr("Show"));
        BTNshowlist.append(BTNshow); 
        connections << connect(BTNshow, SIGNAL(clicked()), showsignalMapper, SLOT(map()), Qt::UniqueConnection);
        showsignalMapper -> setMapping(BTNshow,i);
        layoutmolecules->addWidget(LBLmol,i,0);
        layoutmolecules->addWidget(CHKmol,i,1);
        layoutmolecules->addWidget(BTNmol,i,2);
        layoutmolecules->addWidget(BTNshow,i,3);
        layoutmolecules->addWidget(BTNdeletemol,i,4);
    }
    connections << connect (deletesignalMapper, SIGNAL(mapped(int)), this, SLOT(deletemolecule(int)), Qt::UniqueConnection) ;
    connections << connect (showsignalMapper, SIGNAL(mapped(int)), this, SLOT(showmolecule(int)), Qt::UniqueConnection) ;
}

void glWidget::deletemolecule(int i){
    QMessageBox msgBox;
    msgBox.setInformativeText(tr("Do you want to remove ")+molecules->at(i)->getname()+"?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);
    msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
    msgBox.setButtonText(QMessageBox::No, tr("No"));
    msgBox.setIcon(QMessageBox::Question);
    int ret = msgBox.exec();
    if (ret == QMessageBox::No)
        return;
    delete molecules->at(i);
    molecules->removeAt(i);
    molecule_names.removeAt(i);
    window->resetmeasures();
    window->resetlabaxes();
    if (msrs){
        msrs->reset_all();
        msrs->set_molecules(molecule_names);
    }
    QDLwindow->mespimizerdialog->setmolecules(molecules);
    updateWindowDialog(gdock);
    QDLwindow->mespimizerdialog->setSPBhostmax(molecules->length());
    QDLwindow->mespimizerdialog->setSPBtemplatemax(molecules->length());
    if (molecules->length() < 2){
        templateindex = 1;
        QDLwindow->mespimizerdialog->setSPBtemplate(templateindex);
    }
    if (molecules->isEmpty()){
        QDLwindow->rotationsdialog->reset_rotation();
        QDLwindow->translationsdialog->reset_translation();
    }
    if (window && window->isVisible()){
        updatedisplay();
    }
}

void glWidget::showmolecule(int i){
    if (molecules->at(i)->isvisible()){
        molecules->at(i)->setvisible(false);
        BTNshowlist.at(i)->setText("Show");
    }
    else{
        molecules->at(i)->setvisible(true);
        BTNshowlist.at(i)->setText("Hide");
    }
    molecules->last()->emitupdatedisplay();
}

//          Rotations menu

void glWidget::create_rotation_menu(){
    FRMrotation = QDLwindow->rotationsdialog->getFRMrotation();
    FRMrotation->setVisible(false);
    QDLwindow->rotationsdialog->recorddialog->setmolpath(molpath);
    BTNrotation = new QPushButton(tr("Rotations"));
    connections << connect(BTNrotation, SIGNAL(clicked()), this, SLOT(BTNrotation_clicked()));
    connections << connect(QDLwindow->rotationsdialog,SIGNAL(rotation_changed()),this,SLOT(rotation_changed()));
    connections << connect(QDLwindow->rotationsdialog, SIGNAL(chkrotatechanged(bool)), this, SLOT(CHKrotate_changed(bool)));
    connections << connect(QDLwindow->rotationsdialog, SIGNAL(animationclicked(bool)), this, SLOT(toggleanimation(bool)));
    connections << connect(timer, SIGNAL(timeout()), this, SLOT(animaterotation()));
    connections << connect(QDLwindow->rotationsdialog,SIGNAL(SLDspeed_changed()), this, SLOT(resetinterval()));
    connections << connect(QDLwindow->rotationsdialog->recorddialog, SIGNAL(startrecording()), this, SLOT(startrecording()));
    connections << connect(QDLwindow->rotationsdialog->recorddialog, SIGNAL(stoprecording()), this, SLOT(stoprecording()));
    connections << connect(QDLwindow->rotationsdialog->recorddialog, SIGNAL(recordfilenamechanged(QString)), this, SLOT(set_recordfilename(QString)));
}


//          Translations menu

void glWidget::create_translation_menu(){

    FRMtranslation = QDLwindow->translationsdialog->getFRMtranslation();
    FRMtranslation->setVisible(false);

    BTNtranslation = new QPushButton(tr("Translations"));
    connections << connect(BTNtranslation, SIGNAL(clicked()), this, SLOT(BTNtranslation_clicked()));

    connections << connect(QDLwindow->translationsdialog,SIGNAL(setworldtranslation(QVector3D)),window,SLOT(setworld_translation(QVector3D)));
    connections << connect(QDLwindow->translationsdialog,SIGNAL(stepwheel_changed(float)),window,SLOT(setstepwheel(float)));
}

//          Axes menu

void glWidget::create_axes_menu(){

    FRMaxes = QDLwindow->axesdialog->getFRMaxes();
    BTNaddaxes = new QPushButton(tr("Axes"));

    connections << connect(BTNaddaxes, SIGNAL(clicked()), this, SLOT(BTNaddaxes_clicked()));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axesvisible(bool)), window, SLOT(setaxesvisible(bool)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axeslabelsvisible(bool)), window, SLOT(setaxeslabelsvisible(bool)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axesfont(QFont)), window, SLOT(setfontlabaxeslabels(QFont)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axesarrowsize(int)), window, SLOT(setaxesarrowsize(int)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axesarrowwidth(int)), window, SLOT(setaxesarrowwidth(int)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axeslength(int)), window, SLOT(setaxeslength(int)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(axesthickness(int)), window, SLOT(setaxesthickness(int)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(xaxiscolor(QColor)), window, SLOT(setXaxis_color(QColor)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(yaxiscolor(QColor)), window, SLOT(setYaxis_color(QColor)));
    connections << connect(QDLwindow->axesdialog, SIGNAL(zaxiscolor(QColor)), window, SLOT(setZaxis_color(QColor)));
}

//          Capture menu

void glWidget::create_capture_menu(){

    FRMcapture = QDLwindow->scrshotdialog->getFRMcapture();
    BTNcapture = new QPushButton(tr("Manage capture"));

    connections << connect(BTNcapture, SIGNAL(clicked()), this, SLOT(BTNcapture_clicked()));
    connections << connect(QDLwindow->scrshotdialog, SIGNAL(scaledef_changed()), this, SLOT(adjustQDLSize()));
    connections << connect(QDLwindow->scrshotdialog, SIGNAL(take_picture()), this, SLOT(capture()));
}

//          Lights menu

void glWidget::create_lights_menu(){

    FRMlights = QDLwindow->lightsdialog->getFRMlights();
    BTNlights = new QPushButton(tr("Manage lights"));

    QDLwindow->lightsdialog->setambientcolor(window->getAmbientColor());
    QDLwindow->lightsdialog->setbkgcolor(window->getBkgColor());
    QDLwindow->lightsdialog->setlightcolor(window->getLightColor());
    QDLwindow->lightsdialog->setlightpower(window->getLightPower());
    QDLwindow->lightsdialog->setlinearattenuation(window->getlinearattenuation());
    QDLwindow->lightsdialog->setspecularcolor(window->getSpecularColor());
    QDLwindow->lightsdialog->setspecularindex(window->getSpecularIndex());
    QDLwindow->lightsdialog->setlightsposition(window->getLigthPosition());

    connections << connect(BTNlights,SIGNAL(clicked(bool)),this,SLOT(BTNlights_clicked()));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(ambientcolorchanged(QColor)), window, SLOT(setAmbientColor(QColor)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(bkgcolorchanged(QColor)), window, SLOT(setBkgColor(QColor)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(lightcolorchanged(QColor)), window, SLOT(setLightColor(QColor)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(specularcolorchanged(QColor)), window, SLOT(setSpecularColor(QColor)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(lightschanged(QVector3D)), window, SLOT(setlightsposition(QVector3D)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(lightpowerchanged(float)), window, SLOT(setLightPower(float)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(attenuationchanged(bool)), window, SLOT(setlinearattenuation(bool)));
    connections << connect(QDLwindow->lightsdialog, SIGNAL(specularindexchanged(float)), window, SLOT(setSpecularIndex(float)));

    connections << connect(QDLwindow->lightsdialog, SIGNAL(updateandmoveToTop()), this, SLOT(emitmovetotop()));

}


//          Balls and cylinders menu

void glWidget::create_balls_cyl_menu(){
    FRMballcyl = QDLwindow->ballsandcylsdialog->getFRMballcyl();

    BTNballcyl = new QPushButton();
    BTNballcyl->setText(tr("Manage balls and sticks"));

    connections << connect(BTNballcyl,SIGNAL(clicked(bool)),this,SLOT(BTNballcyl_clicked()));
    connections << connect(QDLwindow->ballsandcylsdialog, SIGNAL(scaleradiichanged(bool)), this, SLOT(CHKscaleradii_changed(bool)));
    connections << connect(QDLwindow->ballsandcylsdialog,SIGNAL(ballcylchanged()),this,SLOT(ballcyl_changed()));
}


//          Viewport menu

void glWidget::create_viewport_menu(){

    FRMviewport = QDLwindow->viewportdialog->getFRMviewport();

    BTNviewport = new QPushButton();
    BTNviewport->setText(tr("Manage viewport"));
    connections << connect(BTNviewport,SIGNAL(clicked(bool)),this,SLOT(BTNviewport_clicked()));

    connections << connect(QDLwindow->viewportdialog,SIGNAL(viewportchanged()),this,SLOT(viewport_changed()));

}


//          Optimize cluster menu

void glWidget::create_optimize_cluster_menu(){

    FRMoptimizeCluster = QDLwindow->mespimizerdialog->getFRMoptimizeCluster();


    QDLwindow->mespimizerdialog->setmolecules(molecules);
    QDLwindow->mespimizerdialog->setTXTframesfile(mespimizerfile);
    QDLwindow->mespimizerdialog->setrecordfile(mespimirecordfilename);
    QDLwindow->mespimizerdialog->setguessfromcanvas(guestfromcanvas);
    QDLwindow->mespimizerdialog->setclustername(clustername);
    QDLwindow->mespimizerdialog->setframesfile(framesfile);
    QDLwindow->mespimizerdialog->setspeed(speed);
    QDLwindow->mespimizerdialog->setinterpolpoints(interpoints);
    QDLwindow->mespimizerdialog->setCHKrecordoptim(recordoptim);
    QDLwindow->mespimizerdialog->setSPBtemplatemax(molecules->length());
    QDLwindow->mespimizerdialog->setSPBtemplate(templateindex);
    QDLwindow->mespimizerdialog->setCHKoptimizeselect(onlyselcp);
    QDLwindow->mespimizerdialog->setdisplayEPIC(displayEPIC);
    QDLwindow->mespimizerdialog->setenergyprecision(energyprecision);
    QDLwindow->mespimizerdialog->setenergyfont(energyfont);
    QDLwindow->mespimizerdialog->setenergycolor(energycolor);
    QDLwindow->mespimizerdialog->sethartree(hartree);

    BTNoptimizeCluster = new QPushButton();
    BTNoptimizeCluster->setText(tr("Optimize cluster"));

    connections << connect(BTNoptimizeCluster,SIGNAL(clicked(bool)),this,SLOT(BTNoptimizeCluster_clicked()));
    connections << connect(QDLwindow->mespimizerdialog,SIGNAL(movetotop()),this,SLOT(emitmovetotop()));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(adjustQDL()), this, SLOT(adjustQDLSize()));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(clusterfile_changed(QString)), this, SLOT(clusterfile_changed(QString)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(framesfile_changed(QString)), this, SLOT(framesfile_changed(QString)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(deletecluster()), this, SLOT(deletecluster()));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(exec_mespimizer()), this, SLOT(exec_mespimizer()));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(optimizecanvas_changed(bool)), this, SLOT(optimizecanvas_changed(bool)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(qmrun(QString)), this, SLOT(qmrun(QString)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(replay(QString)), this, SLOT(replay_mespimization(QString)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(reset(QString)), this, SLOT(reset_mespimization(QString)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(recordfilenamechanged(QString)), this, SLOT(set_recordfilename(QString)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(speed_changed(int)), this, SLOT(speed_changed(int)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(interpol_changed(int)), this, SLOT(interpol_changed(int)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(recordoptim_changed(bool)), this, SLOT(recordoptim_changed(bool)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(template_changed(int)), this, SLOT(template_changed(int)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(optimizeselect_changed(bool)), this, SLOT(optimizeselect_changed(bool)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(displayepic_changed(bool)), this, SLOT(setdisplayEPIC(bool)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(energyprecision_changed(int)), this, SLOT(setenergyprecision(int)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(font_clicked(QFont)), this, SLOT(setenergyfont(QFont)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(fontcolor_clicked(QColor)), this, SLOT(setenergycolor(QColor)));
    connections << connect(QDLwindow->mespimizerdialog, SIGNAL(hartree_units(bool)), this, SLOT(sethartree(bool)));

    if (optimvisible)
        BTNoptimizeCluster_clicked();
}


//          Save geometry menu

void glWidget::create_save_geometry_menu(){
    FRMsaveGeometry= new QGroupBox(tr("Save geometry"));
    FRMsaveGeometry->setVisible(false);

    BTNsaveGeometry = new QPushButton();
    BTNsaveGeometry->setText(tr("Save geometry"));
    connections << connect(BTNsaveGeometry,SIGNAL(clicked(bool)),this,SLOT(BTNsaveGeometry_clicked()));
}

//          Save geometry layouts

void glWidget::create_save_geometry_layouts(){

}

//          Save/Retrieve menu

void glWidget::create_save_retrieve_menu(){
    FRMsettings= new QGroupBox(tr("Settings"));
    FRMsettings->setVisible(false);

    BTNsettings = new QPushButton();
    BTNsettings->setText(tr("Save/retrieve settings"));
    connections << connect(BTNsettings,SIGNAL(clicked(bool)),this,SLOT(BTNsettings_clicked()));

    BTNsaveSettings = new QPushButton();
    BTNsaveSettings->setText(tr("Save"));
    connections << connect(BTNsaveSettings, SIGNAL(clicked()), this, SLOT(BTNsaveSettings_clicked()));
    BTNretrieveSettings = new QPushButton();
    BTNretrieveSettings->setText(tr("Retrieve"));
    connections << connect(BTNretrieveSettings, SIGNAL(clicked()), this, SLOT(BTNretrieveSettings_clicked()));
    CHKsettingspos = new QCheckBox(tr("Position"));
    CHKsettingspos->setChecked(true);
    CHKsettingsrot = new QCheckBox(tr("Rotation"));
    CHKsettingsrot->setChecked(true);
    CHKsettingsbkg = new QCheckBox(tr("Background color"));
    CHKsettingsbkg->setChecked(true);
    CHKsettingslights = new QCheckBox(tr("Lights"));
    CHKsettingslights->setChecked(true);
    CHKsettingsballs = new QCheckBox(tr("Balls and sticks"));
    CHKsettingsballs->setChecked(true);
    CHKsettingsviewport = new QCheckBox(tr("Viewport"));
    CHKsettingsviewport->setChecked(true);
    CHKsettingsmolecs = new QCheckBox(tr("Molecules"));
    CHKsettingsmolecs->setChecked(true);
    connections << connect(CHKsettingsmolecs, SIGNAL(stateChanged(int)), this, SLOT(CHKsettingsmolecs_changed()));
    CHKsettingssurf = new QCheckBox(tr("Surfaces"));
    CHKsettingssurf->setChecked(true);
    CHKsettingsgrids = new QCheckBox(tr("Grids"));
    CHKsettingsgrids->setChecked(true);
    connections << connect(CHKsettingsgrids, SIGNAL(stateChanged(int)), this, SLOT(CHKsettingsgrids_changed()));
    CHKsettingsisosurf = new QCheckBox(tr("Isosurfaces"));
    CHKsettingsisosurf->setChecked(true);
}

//          Save/Retrieve layouts

void glWidget::create_save_retrieve_layouts(){
    QGridLayout *layout1 = new QGridLayout();
    layout1->addWidget(CHKsettingsbkg,0,0);
    layout1->addWidget(CHKsettingslights,0,1);
    layout1->addWidget(CHKsettingsballs,1,0);
    layout1->addWidget(CHKsettingsviewport,1,1);
    layout1->addWidget(CHKsettingsmolecs,2,0);
    layout1->addWidget(CHKsettingssurf,2,1);
    layout1->addWidget(CHKsettingspos,3,0);
    layout1->addWidget(CHKsettingsrot,3,1);
    layout1->addWidget(CHKsettingsgrids,4,0);
    layout1->addWidget(CHKsettingsisosurf,4,1);

    QHBoxLayout *layout2=new QHBoxLayout();
    layout2->addStretch();
    layout2->addWidget(BTNsaveSettings);
    layout2->addWidget(BTNretrieveSettings);
    layout2->addStretch();
    layout2->setAlignment(Qt::AlignCenter);

    QVBoxLayout *layout3=new QVBoxLayout(FRMsettings);
    layout3->addLayout(layout1);
    layout3->addLayout(layout2);
    layout3->setAlignment(Qt::AlignCenter);

}

//********************************************************************************
//          Display functions
//********************************************************************************


//          Create window for 3D images display

void glWidget::create_window(){
    if (window){
        return;
    }
    else{
        window = new glWindow(molecules, this);
        window->makeCurrent();  // This is necessary to prevent erratic behaviour when changing from one display to another
        QSurfaceFormat fmt;         //
        fmt.setSwapInterval(1);     //  This is to synchronize to the vertical refresh rate of the display
        window->setFormat(fmt);     //
        window->setTitle(getWindowName()+" 3D Display");
        window->setlinearattenuation(true);
        window->setworld_translation(QVector3D(0,0,Z_TRANS_INI));
        connections << connect(window,SIGNAL(endrecording()),this,SLOT(stoprecording()));
        connections << connect(window,SIGNAL(endmakingmovie()),this,SLOT(endmakingmovie()));
        connections << connect(window,SIGNAL(emitcheckactivate(int)),this,SLOT(checkactivate(int)));
        connections << connect(window,SIGNAL(update_worldrotation(QQuaternion)),this,SLOT(setrotation(QQuaternion)));
        connections << connect(window,SIGNAL(update_worldtranslation(QVector3D)),this,SLOT(settranslation(QVector3D)));
    }
    window->setlightsposition(lightPosition);
    window->setAmbientColor(ambientcolor);
    window->setSpecularColor(specularcolor);
    window->setBkgColor(bkgcolor);
    window->setLightColor(lightcolor);
    window->setLightPower(lightpower);
    window->setSpecularIndex(specularindex);
    window->setfov(fov);
    window->setznear(zNear);
    window->setzfar(zFar);
    window->setstepwheel(stepwheel);
}

void glWidget::initpointers(){
    arguments = nullpointer;
    msrs = nullpointer;
    myProcess = nullpointer;
    QDLwindow = nullpointer;

    scrollArea = nullpointer;


    strprocess = nullpointer;

    timer = new QTimer();

    window = nullpointer;
}

void glWidget::initQDLpointers(){
    BTNaddaxes = nullpointer;
    BTNballcyl = nullpointer;
    BTNcapture = nullpointer;
    BTNlights = nullpointer;
    BTNmeasures = nullpointer;
    BTNoptimizeCluster = nullpointer;
    BTNretrieveSettings = nullpointer;
    BTNrotation = nullpointer;
    BTNsaveGeometry = nullpointer;
    BTNsaveSettings = nullpointer;
    BTNsettings = nullpointer;
    BTNviewport = nullpointer;

    CHKsettingsballs = nullpointer;
    CHKsettingsbkg = nullpointer;
    CHKsettingsgrids = nullpointer;
    CHKsettingsisosurf = nullpointer;
    CHKsettingslights = nullpointer;
    CHKsettingsmolecs = nullpointer;
    CHKsettingspos = nullpointer;
    CHKsettingsrot = nullpointer;
    CHKsettingssurf = nullpointer;
    CHKsettingsviewport = nullpointer;

    FRMballcyl = nullpointer;
    FRMcapture = nullpointer;
    FRMlights = nullpointer;
    FRMmeasures = nullpointer;
    FRMoptimizeCluster = nullpointer;
    FRMrotation = nullpointer;
    FRMsaveGeometry = nullpointer;
    FRMsettings = nullpointer;
    FRMtranslation = nullpointer;
    FRMviewport = nullpointer;

    layoutmolecules = nullpointer;

    SPBbondthres = nullpointer;
}

bool glWidget::isvisible(){
    return visible;
}

QPoint glWidget::get_position(){
    if (mainWin)
        return QPoint(mainWin->x(),mainWin->y());
    else
        return QPoint(0,0);
}

QSize glWidget::get_size(){
    if (mainWin)
        return mainWin->size();
    else
        return QSize(0,0);
}




void glWidget::set_visible(bool a){
    visible = a;
    if (mainWin){
        if (visible){
            mainWin->show();
            gdock->setVisible(true);
        }
        else{
            gdock->setVisible(false);
            mainWin->hide();
        }
    }
}

void glWidget::set_position(QPoint a){
    if (mainWin){
        position = a;
        mainWin->move(position);
    }
}

void glWidget::set_size(QSize a){
    if (mainWin){
        size = a;
        mainWin->resize(size);
    }
}

//          Update display images

void glWidget::updatedisplay(){
    if (window){
        window->makeCurrent();
        window->reloadbuffers();
        window->update();
    }
}

//          Update display images

void glWidget::updateGL(){
    if (window && window->isVisible()){
        window->update();
        if (!QDLwindow->rotationsdialog->recorddialog->getisrecording() && !animating)
            emit moveToTop(windownumber);
    }
}

//********************************************************************************
//          Geometry measures functions
//********************************************************************************

//  Function BTNmeasures_clicked: open menu for geometry measures
//

void glWidget::BTNmeasures_clicked(){
    if (!msrs) create_measures();
    if (!msrs->QDLmeasures){
        msrs->create_QDLmeasures();
    }
    else if(!msrs->QDLmeasures->isVisible()){
        window->setmeasures(!msrs->getmeasurenone());
        window->setangles(msrs->getmeasureangles());
        window->setdihedrals(msrs->getmeasuredihedrals());
        window->setdistances(msrs->getmeasuredistances());
        msrs->QDLmeasures->show();
    }
    if (msrs->QDLmeasures->isVisible())
        msrs->QDLmeasures->raise();
}

//  End of function BTNmeasures_clicked

//  Function create_measures: creates widget for geometry measures and establishes connections
//
void glWidget::create_measures(){
    msrs = new measures();
    msrs->set_ProjectFolder(ProjectFolder);
    msrs->set_molecules(molecule_names);
    msrs_connections << connect(msrs,SIGNAL(activate_measure(bool)),window,SLOT(setmeasures(bool)));
    msrs_connections << connect(msrs,SIGNAL(angles_precision(int)),window,SLOT(setanglesprecision(int)));
    msrs_connections << connect(msrs,SIGNAL(angles_type(int)),window,SLOT(setanglestype(int)));
    msrs_connections << connect(msrs,SIGNAL(angles_width(int)),window,SLOT(setangleswidth(int)));
    msrs_connections << connect(msrs,SIGNAL(angstrom_changed(bool)),window,SLOT(setangstrom(bool)));
    msrs_connections << connect(msrs,SIGNAL(arc_radius(int)),window,SLOT(setarcradius(int)));
    msrs_connections << connect(msrs,SIGNAL(color_angles(QColor)),window,SLOT(setanglescolor(QColor)));
    msrs_connections << connect(msrs,SIGNAL(color_dihedrals(QColor)),window,SLOT(setdihedralscolor(QColor)));
    msrs_connections << connect(msrs,SIGNAL(color_distances(QColor)),window,SLOT(setdistancescolor(QColor)));
    msrs_connections << connect(msrs,SIGNAL(dihedrals_precision(int)),window,SLOT(setdihedralsprecision(int)));
    msrs_connections << connect(msrs,SIGNAL(dist_precision(int)),window,SLOT(setdistprecision(int)));
    msrs_connections << connect(msrs,SIGNAL(dist_vshift(int)),window,SLOT(setdistvshift(int)));
    msrs_connections << connect(msrs,SIGNAL(dist_transpbkg(bool)),window,SLOT(setdisttranspbkg(bool)));
    msrs_connections << connect(msrs,SIGNAL(draw_arcs(bool)),window,SLOT(setdrawarcs(bool)));
    msrs_connections << connect(msrs,SIGNAL(draw_lines(bool)),window,SLOT(setdrawlines(bool)));
    msrs_connections << connect(msrs,SIGNAL(emit_update_angles()),window,SLOT(emit_update_angles()));
    msrs_connections << connect(msrs,SIGNAL(emit_update_dihedrals()),window,SLOT(emit_update_dihedrals()));
    msrs_connections << connect(msrs,SIGNAL(emit_update_distances()),window,SLOT(emit_update_distances()));
    msrs_connections << connect(msrs,SIGNAL(font_angles(QFont)),window,SLOT(setanglesfont(QFont)));
    msrs_connections << connect(msrs,SIGNAL(font_dihedrals(QFont)),window,SLOT(setdihedralsfont(QFont)));
    msrs_connections << connect(msrs,SIGNAL(font_distances(QFont)),window,SLOT(setdistancesfont(QFont)));
    msrs_connections << connect(msrs,SIGNAL(lines_type(int)),window,SLOT(setlinestype(int)));
    msrs_connections << connect(msrs,SIGNAL(lines_width(int)),window,SLOT(setlineswidth(int)));
    msrs_connections << connect(msrs,SIGNAL(measure_angles(bool)),window,SLOT(setangles(bool)));
    msrs_connections << connect(msrs,SIGNAL(measure_dihedrals(bool)),window,SLOT(setdihedrals(bool)));
    msrs_connections << connect(msrs,SIGNAL(measure_distances(bool)),window,SLOT(setdistances(bool)));
    msrs_connections << connect(msrs,SIGNAL(reset_angles()),window,SLOT(resetangles()));
    msrs_connections << connect(msrs,SIGNAL(reset_dihedrals()),window,SLOT(resetdihedrals()));
    msrs_connections << connect(msrs,SIGNAL(reset_distances()),window,SLOT(resetdistances()));
    msrs_connections << connect(msrs,SIGNAL(show_angles(bool)),window,SLOT(setshowangles(bool)));
    msrs_connections << connect(msrs,SIGNAL(show_dihedrals(bool)),window,SLOT(setshowdihedrals(bool)));
    msrs_connections << connect(msrs,SIGNAL(show_distances(bool)),window,SLOT(setshowdistances(bool)));
    msrs_connections << connect(msrs,SIGNAL(dihedralplanes_color(QColor)),window,SLOT(setdihedralplanescolor(QColor)));
    msrs_connections << connect(msrs,SIGNAL(QDLmeasures_rejected()),this,SLOT(QDLmeasures_rejected()));
    msrs_connections << connect(msrs,SIGNAL(get_systemnames()),this,SLOT(send_systemnames()));

    msrs_connections << connect(window, SIGNAL(select_angle(QVector<centerData>*,QVector<QMatrix4x4>*)),
                                msrs, SLOT(LBLangles_add(QVector<centerData>*,QVector<QMatrix4x4>*)));
    msrs_connections << connect(window, SIGNAL(select_dihedral(QVector<centerData>*,QVector<QMatrix4x4>*)),
                                msrs, SLOT(LBLdihedrals_add(QVector<centerData>*,QVector<QMatrix4x4>*)));
    msrs_connections << connect(window, SIGNAL(select_distance(QVector<centerData>*,QVector<QMatrix4x4>*)),
                                msrs, SLOT(LBLdistances_add(QVector<centerData>*,QVector<QMatrix4x4>*)));
    msrs_connections << connect(window, SIGNAL(update_angles(QVector<centerData>*,QVector<QMatrix4x4>*)), msrs,
                                SLOT(update_angles(QVector<centerData>*,QVector<QMatrix4x4>*)));
    msrs_connections << connect(window, SIGNAL(update_dihedrals(QVector<centerData>*,QVector<QMatrix4x4>*)), msrs,
                                SLOT(update_dihedrals(QVector<centerData>*,QVector<QMatrix4x4>*)));
    msrs_connections << connect(window, SIGNAL(update_distances(QVector<centerData>*,QVector<QMatrix4x4>*)), msrs,
                                SLOT(update_distances(QVector<centerData>*,QVector<QMatrix4x4>*)));
    msrs_connections << connect(window, SIGNAL(remove_distance(QVector<centerData>*,int)), msrs,
                                SLOT(distance_remove(QVector<centerData>*,int)));
    msrs_connections << connect(window, SIGNAL(resetlastselectangles()), msrs,
                                SLOT(resetlastselectangles()));
    msrs_connections << connect(window, SIGNAL(resetlastselectdihedrals()), msrs,
                                SLOT(resetlastselectdihedrals()));
    msrs_connections << connect(window, SIGNAL(resetlastselectdist()), msrs,
                                SLOT(resetlastselectdist()));
}

//  End of function create_measures

void glWidget::QDLmeasures_rejected(){
    if (window)
        window->setmeasures(false);

}

void glWidget::send_systemnames(){
    if (!msrs)
        return;
    msrs->clear_systemnames();
    for (int i = 0 ; i < molecules->length() ; i++){
        msrs->append_systemnames(QFileInfo(molecules->at(i)->getname()).baseName());
    }
}

//********************************************************************************
//          Rotations function
//********************************************************************************

void glWidget::adjustQDLSize(){
    if (QDLwindow){
        QDLwindow->adjustSize();
        QDLwindow->renewwidth(QDLwindow->size());
        gdock->resize(gdock->width()+gdockrewidth,gdock->height());
        gdockrewidth = -gdockrewidth;
    }
}

void glWidget::animaterotation(){
    if (QDLwindow->rotationsdialog->getrotatex())
        setrotation(QQuaternion::fromAxisAndAngle(QVector3D(1.,0.,0.), deltaAngles) * getrotation());
    if (QDLwindow->rotationsdialog->getrotatey())
        setrotation(QQuaternion::fromAxisAndAngle(QVector3D(0.,1.,0.), deltaAngles) * getrotation());
    if (QDLwindow->rotationsdialog->getrotatez())
        setrotation(QQuaternion::fromAxisAndAngle(QVector3D(0.,0.,1.), deltaAngles) * getrotation());
    window->setworld_rotation(rotation);
    for (int i = 0 ; i < molecules->count() ; i ++){
        molecules->at(i)->setworld_rotation(rotation);
    }
    window->update();
}



//  Function BTNrotation_clicked: open menu for rotations
//
void glWidget::BTNrotation_clicked(){
    FRMrotation->setVisible(!FRMrotation->isVisible());
    if (FRMrotation->isVisible()){
        BTNrotation->setStyleSheet("QPushButton {color: red;}");
        BTNrotation->setText(tr("Hide rotations manager"));
    }
    else{
        BTNrotation->setStyleSheet("QPushButton {color: black;}");
        BTNrotation->setText(tr("Rotations"));
    }
    QDLwindow->update();
    QDLwindow->adjustSize();
    QDLwindow->renewwidth(gdock->size());
}
//  End of function BTNrotation_clicked


void glWidget::CHKrotate_changed(bool animate){
    if (!animate){
        timer->stop();
        return;
    }
//    window->setworld_rotation(rotation);
//    window->update();
    if (QDLwindow->rotationsdialog->getstartanimation())
        timer->start(interval);
}

void glWidget::resetinterval(){
    int speed = QDLwindow->rotationsdialog->getspeed();
    interval = MAX_INTERVAL - dltinterval * speed;
    deltaAngles = 1.f;
    if (speed > INTERVAL_SCALE/4)
        deltaAngles = 2.f;
    if (speed > INTERVAL_SCALE/2)
        deltaAngles = 4.f;
    if (speed > 3*INTERVAL_SCALE/4)
        deltaAngles = 8.f;
    if (QDLwindow->rotationsdialog->getstartanimation()){
        timer->stop();
        timer->start(interval);
    }
}


void glWidget::rotation_changed(){
    rotation = QDLwindow->rotationsdialog->getrotation();
    window->setworld_rotation(rotation);
    for (int i = 0 ; i < molecules->count() ; i ++){
        molecules->at(i)->setworld_rotation(rotation);
    }
    window->update();
}

QQuaternion glWidget::getrotation(){
    return rotation;
}

QVector3D glWidget::getrotationAxis(){
    return rotationAxis;
}

void glWidget::setrotation(QQuaternion a){
    rotation = a;
    QDLwindow->rotationsdialog->setrotationButtons(a);
}

void glWidget::setrotationAxis(QVector3D a){
    if (a.length() > 1.e-7)
        rotationAxis = a;
    else
        rotationAxis = QVector3D(1,0,0);
}

void glWidget::settranslation(QVector3D a){
    QDLwindow->translationsdialog->set_translation(a);
}

void glWidget::toggleanimation(bool startanimation){
    if (startanimation)
            timer->start(interval);
        else
            timer->stop();
}

//********************************************************************************
//          Translations function
//********************************************************************************

//  Function BTNtranslation_clicked: open menu for translations
//
void glWidget::BTNtranslation_clicked(){
    FRMtranslation->setVisible(!FRMtranslation->isVisible());
    if (FRMtranslation->isVisible()){
        BTNtranslation->setStyleSheet("QPushButton {color: red;}");
        BTNtranslation->setText(tr("Hide translations manager"));
    }
    else{
        BTNtranslation->setStyleSheet("QPushButton {color: black;}");
        BTNtranslation->setText(tr("Translations"));
    }
    QDLwindow->update();
    QDLwindow->adjustSize();
    QDLwindow->renewwidth(gdock->size());
}

//  End of function BTNtranslation_clicked

//********************************************************************************
//          Capture function
//********************************************************************************

void glWidget::capture(){
    if (!window || !window->isVisible()){
        QMessageBox msgBox;
        msgBox.setText(tr("capture"));
        msgBox.setInformativeText(tr("You must display the image to enable capture"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        emit moveToTop(windownumber);
        return;
    }
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    CaptureFolder = QDLwindow->scrshotdialog->getCaptureFolder();
    filedialog.setDirectory(CaptureFolder);
    QString fileName = filedialog.getSaveFileName(this,tr("Save image as ..."), CaptureFolder,
            (tr("Image files") + ": *.png, *.jpg, *.bmp, *.jpeg, *.ppm, *.xbm, *.xpm, *.tiff"
                + " (*.png *.jpg *.bmp *.jpeg *.ppm *.xbm *.xpm *.tiff);;"+tr("All files")+" (*)"));
    if (fileName.isEmpty()){
        QMessageBox msgBox;
        msgBox.setText(tr("capture"));
        msgBox.setInformativeText(tr("A file name must be supplied for image saving"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        emit moveToTop(windownumber);
        return;
    }
    CaptureFolder = QFileInfo(fileName).path();
    QPoint windowpos = QPoint(window->x(), window->y());
    QSize windowsize = QSize(window->size());
    window->makeCurrent();
    window->reloadbuffers();
    window->setPosition(windowpos);
    window->setlinearattenuation(true);
    window->resize(windowsize);
    window->setlightsposition(lightPosition);
    window->setAmbientColor(ambientcolor);
    window->setSpecularColor(specularcolor);
    window->setLightColor(lightcolor);
    window->setLightPower(lightpower);
    window->setSpecularIndex(specularindex);
    window->setimagequality(QDLwindow->scrshotdialog->getimagequality());
    window->settransparentbg(QDLwindow->scrshotdialog->gettransparentbkg());
    window->setimagefilename(fileName);
    double scalesize = QDLwindow->scrshotdialog->getscalesize().toDouble();
    if (QDLwindow->scrshotdialog->getscaledef() && scalesize > 0.)
        window->paintGLbuff(scalesize * window->size());
    else
        window->paintGLbuff(window->size());
    emit moveToTop(windownumber);
}

//********************************************************************************
//          Axes SLOTS
//********************************************************************************

//  Function BTNaxes_clicked: open menu for axes display
//
void glWidget::BTNaddaxes_clicked(){
    FRMaxes->setVisible(!FRMaxes->isVisible());
    if (FRMaxes->isVisible()){
        BTNaddaxes->setStyleSheet("QPushButton {color: red;}");
        BTNaddaxes->setText(tr("Hide axes manager"));
    }
    else{
        BTNaddaxes->setStyleSheet("QPushButton {color: black;}");
        BTNaddaxes->setText(tr("Axes"));
    }
    QDLwindow->update();
    QDLwindow->adjustSize();
    QDLwindow->renewwidth(gdock->size());
}



//********************************************************************************
//          Capture SLOTS
//********************************************************************************

void glWidget::BTNcapture_clicked(){
        FRMcapture->setVisible(!FRMcapture->isVisible());
        if (FRMcapture->isVisible()){
            BTNcapture->setStyleSheet("QPushButton {color: red;}");
            BTNcapture->setText(tr("Hide capture manager"));
        }
        else{
            BTNcapture->setStyleSheet("QPushButton {color: black;}");
            BTNcapture->setText(tr("Manage capture"));
        }
        QDLwindow->update();
        QDLwindow->adjustSize();
        QDLwindow->renewwidth(gdock->size());
}


//********************************************************************************
//          Record SLOTS
//********************************************************************************


void glWidget::endmakingmovie(){

    if (QDLwindow->mespimizerdialog->getmakingmovie()){
        QDLwindow->mespimizerdialog->endmakingmovie();
    }
    else{
        QDLwindow->rotationsdialog->recorddialog->endmakingmovie();
    }
    QDLwindow->repaint();
    timer->stop();
}

void glWidget::set_recordfilename(QString a){
    recordfilename = a;
}

void glWidget::startrecording(){
    window->setrecordcommand(QDLwindow->rotationsdialog->recorddialog->getrecordcommand());
    window->setframeknt(QDLwindow->rotationsdialog->recorddialog->getframeknt());
    window->setnumframes(QDLwindow->rotationsdialog->recorddialog->getnumframes());
    window->setrecord(QDLwindow->rotationsdialog->recorddialog->getrecord());
    window->setrecordfilename(QDLwindow->rotationsdialog->recorddialog->getrecordfilename());
    window->setremoveframes(QDLwindow->rotationsdialog->recorddialog->getremoveframes());
    window->raise();
    bool isrecording = (QDLwindow->rotationsdialog->getrotatex() ||
                        QDLwindow->rotationsdialog->getrotatey() ||
                        QDLwindow->rotationsdialog->getrotatez());
    for (int i = 0 ; i < molecules->count() ; i++){
        bool ok;
        ok = molecules->at(i)->startanimate();
        isrecording = (isrecording || ok);
    }
    QDLwindow->rotationsdialog->recorddialog->BTNstartrecording_update(isrecording);
    QDLwindow->adjustSize();
    timer->start(interval);
}

void glWidget::stoprecording(){
    QDLwindow->rotationsdialog->recorddialog->BTNstartrecording_update(false);
    QDLwindow->repaint();
    for (int i = 0 ; i < molecules->count() ; i++){
        molecules->at(i)->stopanimate();
    }
    if (window)
        window->setrecord(false);
    timer->stop();
}

//********************************************************************************
//          Lights SLOTS
//********************************************************************************

void glWidget::BTNlights_clicked(){
        FRMlights->setVisible(!FRMlights->isVisible());
        if (FRMlights->isVisible()){
            BTNlights->setStyleSheet("QPushButton {color: red;}");
            BTNlights->setText(tr("Hide lights manager"));
        }
        else{
            BTNlights->setStyleSheet("QPushButton {color: black;}");
            BTNlights->setText(tr("Manage lights"));
        }
        QDLwindow->update();
        QDLwindow->adjustSize();
        QDLwindow->renewwidth(gdock->size());
}



//********************************************************************************
//          Balls and cylinders SLOTS
//********************************************************************************

void glWidget::BTNballcyl_clicked(){
        FRMballcyl->setVisible(!FRMballcyl->isVisible());
        if (FRMballcyl->isVisible()){
            FRMballcyl->adjustSize();
            BTNballcyl->setStyleSheet("QPushButton {color: red;}");
            BTNballcyl->setText(tr("Hide balls and sticks manager"));
        }
        else{
            BTNballcyl->setStyleSheet("QPushButton {color: black;}");
            BTNballcyl->setText(tr("Manage balls and sticks"));
        }
        QDLwindow->update();
        QDLwindow->adjustSize();
        QDLwindow->renewwidth(gdock->size());
}

void glWidget::ballcyl_changed(){
    for (int i = 0 ; i < molecules->count() ; i++){
        ballradius = QDLwindow->ballsandcylsdialog->getballrad();
        cylradius = QDLwindow->ballsandcylsdialog->getcylrad();
        disthressq = pow((QDLwindow->ballsandcylsdialog->getbondthres()) * ANGSTROM_TO_BOHR,2);
        molecules->at(i)->setballradius(ballradius);
        molecules->at(i)->setcylradius(cylradius);
        molecules->at(i)->setdisthressq(disthressq);
        molecules->at(i)->makeStructureBondsandSticks();
    }
    updatedisplay();
}

void glWidget::CHKscaleradii_changed(bool scaleradii){
    for (int i = 0 ; i < molecules->count() ; i++){
        molecules->at(i)->setscaleradii(scaleradii);
        molecules->at(i)->makeStructureBondsandSticks();
    }
    updatedisplay();
}

//********************************************************************************
//          Viewport SLOTS
//********************************************************************************

void glWidget::BTNviewport_clicked(){
        FRMviewport->setVisible(!FRMviewport->isVisible());
        if (FRMviewport->isVisible()){
            BTNviewport->setStyleSheet("QPushButton {color: red;}");
            BTNviewport->setText(tr("Hide viewport manager"));
        }
        else{
            BTNviewport->setStyleSheet("QPushButton {color: black;}");
            BTNviewport->setText(tr("Manage viewport"));
        }
        QDLwindow->update();
        QDLwindow->adjustSize();
        QDLwindow->renewwidth(gdock->size());
}


void glWidget::updateviewport(){
    if (QDLwindow){
        QDLwindow->viewportdialog->setfov(fov);
        QDLwindow->viewportdialog->setznear(zNear);
        QDLwindow->viewportdialog->setzfar(zFar);
    }
    if (window){
        window->setfov(fov);
        window->setznear(zNear);
        window->setzfar(zFar);
        window->update();
        emit moveToTop(windownumber);
    }
}

void glWidget::viewport_changed(){
    fov = QDLwindow->viewportdialog->getfov();
    zNear = QDLwindow->viewportdialog->getznear();
    zFar = QDLwindow->viewportdialog->getzfar();
    if (window){
        window->setfov(fov);
        window->setznear(zNear);
        window->setzfar(zFar);
        window->update();
        emit moveToTop(windownumber);
    }
}


//********************************************************************************
//          Optimize cluster (MESPIMIZER) SLOTS
//********************************************************************************


//  Function addcluster: adds a cluster
//
bool glWidget::addcluster(QString filename){
    if (window)
        window->makeCurrent();
    molecules->append(new molecule());
    if (!loadclustergeometry(filename)){
        molecules->removeLast();
        return false;
    }
    else{
        molecules->last()->loadcharges(znuc);
        molecules->last()->loadxyz(xyz);
        molecules->last()->setpath(QDLwindow->mespimizerdialog->getmespimizerpath());
//        molecules->last()->setrotation(rotation.normalized());
        molecules->last()->setvisible(true);
//        molecules->last()->initatomactive(false);
        molecules->last()->set_ProjectFolder(ProjectFolder);
        molecules->last()->set_ProjectName(ProjectName);
        window->setdisplayEPIC(displayEPIC);
        connections << connect(molecules->last(), SIGNAL(updatedisplay()), this, SLOT(updatedisplay()));
        connections << connect(molecules->last(), SIGNAL(updateGL()), this, SLOT(updateGL()));
        connections << connect(molecules->last(), SIGNAL(animate(bool)), this, SLOT(animate(bool)));
        return true;
    }
}
//  End of function addcluster
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function addclusterTowindow: adds a cluster to the window
//
bool glWidget::addclusterTowindow(QString filename){
    if (addcluster(filename)){
        QString molname = QFileInfo(filename).completeBaseName();
        mespimizerpath = QDLwindow->mespimizerdialog->getmespimizerpath();
        molecules->last()->setname(molname);
        molecules->last()->setfullname(QFileInfo(filename).canonicalPath()
                    +"/"+QFileInfo(filename).completeBaseName());
        molecules->last()->setballradius(ballradius);
        molecules->last()->setcylradius(cylradius);
        molecules->last()->setdisthressq(disthressq);
        molecules->last()->setscaleradii(scaleradii);
        molecules->last()->makeStructureBondsandSticks();
        molecules->last()->setvisible(true);
        molecules->last()->set_ProjectFolder(mespimizerpath);
        molecules->last()->set_ProjectName(molname);
        molecules->last()->setiscluster(true);
        molecule_names.append(molname);
        if (msrs){
            msrs->set_molecules(molecule_names);
        }
        for (int i = 0 ; i < molecules->length()-1 ; i++){
            molecules->at(i)->setvisible(false);
            BTNshowlist.at(i)->setText("Show");
        }
        updatedisplay();
        if (QDLwindow){
            updateWindowDialog(gdock);
        }
        this->moveToTop(windownumber);
        return true;
    }
    else{
        return false;
    }
}
//  End of function addclusterTowindow
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function addemptyclusterTowindow: adds an empty cluster to the window
//
bool glWidget::addemptyclusterTowindow(QString filename){
    mespimizerpath = QDLwindow->mespimizerdialog->getmespimizerpath();
    molecules->append(new molecule());
    molecules->last()->setname(clustername);
    molecules->last()->setpath(mespimizerpath);
    molecules->last()->setfullname(mespimizerpath + "/" +clustername);
    molecules->last()->set_ProjectFolder(mespimizerpath);
    molecules->last()->set_ProjectName(clustername);
    molecules->last()->setiscluster(true);
    molecules->last()->setballradius(ballradius);
    molecules->last()->setcylradius(cylradius);
    molecules->last()->setdisthressq(disthressq);
    molecules->last()->setscaleradii(scaleradii);
    molecules->last()->setvisible(true);
    molecules->last()->setrotation(rotation.normalized());
    molecules->last()->initatomactive(false);
    window->setdisplayEPIC(displayEPIC);
    connections << connect(molecules->last(), SIGNAL(updatedisplay()), this, SLOT(updatedisplay()));
    connections << connect(molecules->last(), SIGNAL(updateGL()), this, SLOT(updateGL()));
    connections << connect(molecules->last(), SIGNAL(animate(bool)), this, SLOT(animate(bool)));

    molecule_names.append("cluster");
    if (msrs){
        msrs->set_molecules(molecule_names);
    }
    for (int i = 0 ; i < molecules->length()-1 ; i++){
        molecules->at(i)->setvisible(false);
        BTNshowlist.at(i)->setText("Show");
    }
    updatedisplay();
    mespimizerfile = QFileInfo(filename).fileName();
    mespimirecordfilename = QDLwindow->mespimizerdialog->getrecordfilename();
    QDLwindow->mespimizerdialog->setmolecules(molecules);
    updateWindowDialog(gdock);
//    BTNoptimizeCluster_clicked();
    this->moveToTop(windownumber);
    return true;
}
//  End of function addemptyclusterTowindow
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function BTNoptimizeCluster_clicked: open menu for cluster optimization
//

void glWidget::BTNoptimizeCluster_clicked(){
    int host = QDLwindow->mespimizerdialog->gethost();
    if (molecules->length() < 1){
        QDLwindow->mespimizerdialog->enableBTNmespimize(false);
        QDLwindow->mespimizerdialog->setmespimizerpath(ProjectFolder);
        QDLwindow->mespimizerdialog->sethostname("");
    }
    else{
        if (molecules->last()->getiscluster()){
            QDLwindow->mespimizerdialog->enableBTNmespimize(false);
        }
        else{
            QDLwindow->mespimizerdialog->enableBTNmespimize(true);
        }
        QDLwindow->mespimizerdialog->setmespimizerpath(molecules->at(host-1)->getpath());
        QDLwindow->mespimizerdialog->sethostname(molecules->at(host-1)->getname());
    }
    FRMoptimizeCluster->setVisible(!FRMoptimizeCluster->isVisible());
    if (FRMoptimizeCluster->isVisible()){
        BTNoptimizeCluster->setStyleSheet("QPushButton {color: red;}");
        BTNoptimizeCluster->setText(tr("Hide optimize cluster"));
        optimvisible = true;
    }
    else{
        BTNoptimizeCluster->setStyleSheet("QPushButton {color: black;}");
        BTNoptimizeCluster->setText(tr("Optimize cluster"));
        optimvisible = false;
    }
//    QDLwindow->mespimizerdialog->setguessfromcanvas(guestfromcanvas);
    QDLwindow->update();
    QDLwindow->adjustSize();
    QDLwindow->renewwidth(gdock->size());
}
//  End of function BTNoptimizeCluster_clicked
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function clusterfile_changed: updates cluster name
//
void glWidget::clusterfile_changed(QString a){
    clustername = a;
}
//  End of function clusterfile_changed
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function deletecluster: removes cluster from molecules list
//
void glWidget::deletecluster(){
    delete molecules->at(molecules->length()-1);
    molecules->removeLast();
    molecule_names.removeLast();
    window->resetmeasures();
    window->setdisplayEPIC(false);
    if (msrs){
        msrs->reset_all();
        msrs->set_molecules(molecule_names);
    }
    if (molecules->isEmpty()){
        QDLwindow->rotationsdialog->reset_rotation();
        QDLwindow->translationsdialog->reset_translation();
    }
    if (window && window->isVisible()){
        updatedisplay();
    }
}
//  End of function deletecluster
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function displayClusterGeometry: displays frame in animation
//
void glWidget::displayClusterGeometry(QString str){
    QString filename = mespimizerpath + "/"+clustername.trimmed()+".xyz_curr_frame";
    if (QFileInfo::exists(filename)){
        if (!cluster_exists){
            QString save_mespimizerpath = mespimizerpath;
            if (!addclusterTowindow(filename)){
                cluster_exists = false;
                return;
            }
            mespimizerpath = save_mespimizerpath;
            cluster_exists = true;
        }
        else{
            if (reloadxyz(filename)){
                molecules->last()->loadxyz(xyz);
                molecules->last()->loadcharges(znuc);
                molecules->last()->makeStructureBondsandSticks();
            }
        }
        updatedisplay();
    }
}
//  End of function displayClusterGeometry
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function exec_mespimizer: executes cluster optimization
//
void glWidget::exec_mespimizer(){
    cluster_exists = false;
    QString stdinput = QDLwindow->mespimizerdialog->getmespimizerpath() + "/mespimizer.inp";
    QString stdoutput = QDLwindow->mespimizerdialog->getmespimizerpath() + "/mespimizer.out";
    if (!strprocess) strprocess = new QString();
    if (!arguments) arguments = new QStringList();
    *strprocess = QString("DAMMESPIMIZER.exe");
    QString execName = get_execName(*strprocess, QString("TDAM320"));
    if (execName.isEmpty())
        return;
    arguments->clear();
    *arguments << " < " << stdinput << " > " << stdoutput;
    if (!myProcess) myProcess = new QProcess(this);
    myProcess->setStandardInputFile(stdinput);
    myProcess->setStandardOutputFile(stdoutput,QIODevice::Truncate);
    myProcess->setStandardErrorFile(stdoutput,QIODevice::Append);
    connect(myProcess, SIGNAL(started()), this, SLOT(processStart()));
    connect(myProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(processError(QProcess::ProcessError)));
    connect(myProcess, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(processOutput(int,QProcess::ExitStatus)));
//    qDebug() << "myProcess = " << execName << *arguments;
    myProcess->start(execName,*arguments);
    QFileSystemWatcher *watcher = new QFileSystemWatcher();
    connect(watcher, SIGNAL(fileChanged(QString)), this, SLOT(displayClusterGeometry(QString)));
    mespimizerpath = QDLwindow->mespimizerdialog->getmespimizerpath();
    watcher->addPath(QDLwindow->mespimizerdialog->getmespimizerpath() + "/mespimizer.out");
    emit moveToTop(windownumber);
}
//  End of function exec_mespimizer
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function framesfile_changed: updates frames file name
//
void glWidget::framesfile_changed(QString a){
    framesfile = a;
}
//  End of function clusterfile_changed
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//  Function get_execName: creates execution command
//
QString glWidget::get_execName(QString processname, QString subdir){
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
//  End of function get_execName
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function interpol_changed: updates number of interpolation points
//
void glWidget::interpol_changed(int i){
    interpoints = i;
}
//  End of function interpol_changed
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function linearinterpolation: linear interpolation of positions
//
void glWidget::linearinterpolation(){
    QVector<QVector3D > xyzant = molecules->last()->getxyz();
    double dlt = 1. / ninterpol;
    delay = QDLwindow->mespimizerdialog->getdelay();
//    qDebug() << "delay = " << delay;
    for (int i = 0; i < ninterpol+1 ; i++){
        QVector<QVector3D > xyzi;
        for (int j = 0 ; j < xyz.length(); j++){
            xyzi.append(xyzant[j] * (1 - dlt*i) + xyz[j] * dlt * i);
        }
        if (msrs)
            msrs->update_measure_centers(&xyzi);
        molecules->last()->loadxyz(xyzi);
        molecules->last()->makeStructureBondsandSticks();

        window->setframeknt(1);     // Dirty trick to guarantee that takes all frames
        updatedisplay();
        pause(delay);
    }
}
//  End of function linearinterpolation
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function loadclustergeometry: loads geometry and charges of centers
//
bool glWidget::loadclustergeometry(QString filename){
    if (!reloadxyz(filename)){
        return false;
    }
    else{
        return true;
    }
}
//  End of function loadclustergeometry
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function optimizecanvas_changed: sets guestfromcanvas
//
void glWidget::optimizecanvas_changed(bool a){
    guestfromcanvas = a;
}
//  End of function optimizecanvas_changed
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function optimizecanvas_changed: sets guestfromcanvas
//
void glWidget::optimizeselect_changed(bool a){
    onlyselcp = a;
}

//
// Function pause: pausing execution
void glWidget::pause( int millisecondsToWait )
{
    QTime dieTime = QTime::currentTime().addMSecs( millisecondsToWait );
    while( QTime::currentTime() < dieTime )
    {
        QCoreApplication::processEvents( QEventLoop::AllEvents, 100 );
    }
}
//  End of function pause
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function processError
void glWidget::processError(QProcess::ProcessError error)
{
    mainWin->statusBar()->setStyleSheet("color: red");
    mainWin->statusBar()->showMessage(tr("End of computation"));
    pause(1000);
    mainWin->statusBar()->clearMessage();
    mainWin->statusBar()->setStyleSheet("color: black");
    if(error==QProcess::FailedToStart){
        QMessageBox::critical(this,QString("Error %1").arg(error),QString(tr("Could not start running program %1")).arg(*strprocess));
    }
    else{
        QMessageBox::critical(this,QString("Error %1").arg(error),QString(tr("Error when running program %1")).arg(*strprocess)
                    + QString(arguments->join(" ")));
    }
    delete myProcess;
    myProcess = nullpointer;
    QDLwindow->mespimizerdialog->enableBTNmespimize(true);
    emit moveToTop(windownumber);
}
//  End of function processError
//  ---------------------------------------------------------------------------------------------------------------------------


//  Function processOutput
void glWidget::processOutput(int exitCode, QProcess::ExitStatus exitStatus)
{
    mainWin->statusBar()->setStyleSheet("color: red");
    mainWin->statusBar()->showMessage(tr("End of computation"));
    pause(1000);
    mainWin->statusBar()->clearMessage();
    mainWin->statusBar()->setStyleSheet("color: black");
    if(exitStatus==QProcess::NormalExit){
        if(QFileInfo(mespimizerpath+"/mespimizer.err").exists()){
            QFile file(mespimizerpath+"/mespimizer.err");
            if (!file.open(QFile::ReadOnly | QFile::Text)) {
                return;
            }
            QTextStream in(&file);
            QString line = in.readLine();
#if QT_VERSION < 0x050E00
            QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
            QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
            if (fields.count() < 2)
                return;
            int iter = fields.takeFirst().toInt();
            int molec = fields.takeFirst().toInt();
            QMessageBox::information(this, tr("MESPIMIZER"),
                QString(tr("Cluster optimization aborted, highest number of iterations: %1\nreached in molecule %2"))
                        .arg(iter).arg(molec));
        }
        else{
            QMessageBox::information(this, tr("MESPIMIZER"),tr("Cluster optimization finished"));
        }

    }else if (exitStatus==QProcess::CrashExit){
        QMessageBox::warning(this, tr("MESPIMIZER"),
                tr(QString("MESPIMIZER crashed, exit code = %1").arg(exitCode).toLatin1()));
    }
    delete myProcess;
    myProcess = nullpointer;
    QDLwindow->mespimizerdialog->enableBTNmespimize(true);
    emit moveToTop(windownumber);
}
//  End of function processOutput
//  ---------------------------------------------------------------------------------------------------------------------------


//  Function processStart:  Starts an external process
void glWidget::processStart()
{
    mainWin->statusBar()->setStyleSheet("color: blue");
    mainWin->statusBar()->showMessage(tr("Computing..."));
}
//  End of function processStart
//  ---------------------------------------------------------------------------------------------------------------------------


//  Function quaterninterpolation: Interpolation of rotations by means of Quaternions
//
void glWidget::quaterninterpolation(){
//    qDebug() << "Ener0 = " << Ener0;
//    qDebug() << "Ener1 = " << Ener1;
    QVector<QVector3D > xyzant = molecules->last()->getxyz();
    QVector<QVector3D > cmv;
    QVector<QVector3D > cmantv;
    for (int kmol = 0 ; kmol < indmols.length()-1 ; kmol++){
        QVector3D cmant = QVector3D(0,0,0);
        QVector3D cm = QVector3D(0,0,0);
        int ztot = 0;
        for (int i = indmols[kmol]; i < indmols[kmol+1] ; i++){
            cmant += znuc[i] * xyzant[i];
            cm += znuc[i] * xyz[i];
            ztot += znuc[i];
        }
        cmant /= ztot;
        cmantv << cmant;
        cm /= ztot;
        cmv << cm;
    }
    QVector<QVector3D > viv;
    QVector<QVector3D > vipv;
    QVector<QVector3D > vjv;
    QVector<QVector3D > vjpv;
    QVector<int> transformation;    // transformations are defined as:
                                    //      0: no change;
                                    //      1: only translation
                                    //      2: rotation with Z and Z' coincident
                                    //      3: rotation with Z and Z' opposed
                                    //      4: general rotation
    QVector<QQuaternion> qtotv;
    QQuaternion q1;
    for (int kmol = 0 ; kmol < indmols.length()-1 ; kmol++){
        int knt = 0;
        QVector3D vi = QVector3D(0,0,0);
        QVector3D vip=  QVector3D(0,0,0);
        for (int i = indmols[kmol] ; i < indmols[kmol+1] ; i++){
            if ((xyz.at(i)-xyzant.at(i)).lengthSquared() < 1.e-10)
                continue;
            vi =xyzant.at(i)-cmantv[kmol];
            if (vi.length() > 1.e-5){   // center defining X axis must not be at the charges center
                knt = i+1;
                vi.normalize();
                vip = xyz.at(i) - cmv[kmol];
                vip.normalize();
//                qDebug() << "i = " << i << "vi = " << vi << " vip " << vip;
                break;
            }
        }
        viv << vi;
        vipv << vip;
        if (vi.lengthSquared() < 1.e-10){   // Geometry of molecule does not change at all
            transformation << 0;
            qtotv << q1;
        }
        else{
            QVector3D vj = QVector3D(0,0,0);
            QVector3D vjp = QVector3D(0,0,0);
            for (int i = knt; i < indmols[kmol+1] ; i++){
                vj = xyzant.at(i)-cmantv[kmol];
                vj -= QVector3D::dotProduct(vj,vi) * vi;
                if (vj.length() > 1.e-5){
                    vj /= vj.length();
                    vj.normalize();
                    vjp = xyz.at(i) - cmv[kmol];
                    vjp -= QVector3D::dotProduct(vjp,vip) * vip;
        //            vjp /= vjp.length();
                    vjp.normalize();
//                    qDebug() << "i = " << i << "vj = " << vj << " vjp " << vjp;
                    break;
                }
            }
            if (vj.lengthSquared() < 1.e-10){   // The unlikely case that only one center is rotated
                vj = QVector3D(vi[1],-vi[0],0);
                vjp = QVector3D(vip[1],-vip[0],0);
                if (vj.lengthSquared() < 1.e-10){
                    vj = QVector3D(vi[2],0,0);
                    vjp = QVector3D(vip[2],0,0);
                }
            }
            vjv << vj;
            vjpv << vjp;
            if ((vi-vip).lengthSquared() < 1.e-10 && (vj-vjp).lengthSquared() < 1.e-10){  // Only translation is carried out
//                qDebug() << "only translation";
                transformation << 1;
                qtotv << q1;
            }
            else{
//                qDebug() << "hay rotacion: (vi-vip).lengthSquared() = " << (vi-vip).lengthSquared()
//                         << "  (vj-vjp).lengthSquared() = " << (vj-vjp).lengthSquared();
                QVector3D vk = QVector3D::crossProduct(vi,vj);
                QVector3D vkp = QVector3D::crossProduct(vip,vjp);
                if ((vk-vkp).lengthSquared() < 1.e-10){     // Z and Z' axes aligned
                    transformation << 2;
//                    qDebug() << "rotation Z and Zp aligned";
                    float alfa;
                    float cosalfa = QVector3D::dotProduct(vi,vip);
                    float sinalfa = QVector3D::dotProduct(QVector3D::crossProduct(vi,vip),vk);
                    if (sinalfa >= 0){
                        alfa = acos(cosalfa);
                    }
                    else{
                        alfa = -acos(cosalfa);
                    }
//                    qDebug() << "alfa = " << alfa << "cosalfa = " << cosalfa;
                    qtotv << QQuaternion::fromAxisAndAngle(vk,alfa);
                }
                else if((vk+vkp).lengthSquared() < 1.e-10){     // Z and Z' axes opposed
                    transformation << 3;
//                    qDebug() << "rotation Z and Zp opposed";
                    float alfa;
                    float cosalfa = QVector3D::dotProduct(vj,vjp);
                    float sinalfa = QVector3D::dotProduct(QVector3D::crossProduct(vj,vjp),vk);
                    if (sinalfa >= 0){
                        alfa = acos(cosalfa);
                    }
                    else{
                        alfa = -acos(cosalfa);
                    }
                    float beta = M_PI;
                    QQuaternion q2 = QQuaternion::fromAxisAndAngle(vk,alfa);
                    QQuaternion q3 = QQuaternion::fromAxisAndAngle(vjp,beta);
                    qtotv << q3*q2;
                }
                else{
                    transformation << 4;
                    QVector3D un = QVector3D::crossProduct(vk,vkp);
//                    qDebug() << "general rotation";
        //            un.normalize();
                    un /= un.length();
                    float alfa;
                    float cosalfa = QVector3D::dotProduct(vj,un);
                    float sinalfa = QVector3D::dotProduct(QVector3D::crossProduct(vj,un),vk);
                    if (sinalfa >= 0){
                        alfa = acos(cosalfa) * 180. / M_PI;
                    }
                    else{
                        alfa = -acos(cosalfa) * 180. / M_PI;
                    }
                    float beta = acos(QVector3D::dotProduct(vk,vkp)) * 180. / M_PI;
                    float gamma;
                    float cosgam =  QVector3D::dotProduct(vjp,un);
                    float singam = QVector3D::dotProduct(QVector3D::crossProduct(un,vjp),vkp);
                    if (singam >= 0){
                        gamma = acos(cosgam) * 180. / M_PI;
                    }
                    else{
                        gamma = -acos(cosgam) * 180. / M_PI;
                    }
                    QQuaternion q2 = QQuaternion::fromAxisAndAngle(vk,alfa);
                    QQuaternion q3 = QQuaternion::fromAxisAndAngle(un,beta);
                    QQuaternion q4 = QQuaternion::fromAxisAndAngle(vkp,gamma);
                    qtotv << q4*q3*q2;
//                    qDebug() << "alfa = " << alfa << "  beta = " << beta << " gamma = " << gamma;
                }
            }
        }
    }
    QQuaternion qinterp;
    QVector3D trs;
    double dlt = 1. / ninterpol;
    delay = QDLwindow->mespimizerdialog->getdelay();
//    qDebug() << "delay = " << delay;
    for (int i = 0; i < ninterpol ; i++){
        QVector<QVector3D > xyzi;
        for (int kmol = 0 ; kmol < indmols.length()-1 ; kmol++){
//            qDebug() << "mole: " << kmol << " case = " << transformation[kmol];
            switch (transformation[kmol]){
                case 0:         // No transformation
                    for (int j = indmols[kmol] ; j < indmols[kmol+1]; j++){
                        xyzi.append(xyz[j]);
                    }
                    break;
                case 1:         // Only translation
                    for (int j = indmols[kmol] ; j < indmols[kmol+1]; j++){
                        xyzi.append(xyzant[j] * (1 - dlt*i) + xyz[j] * dlt * i);
                    }
                    break;
                case 2:         // rotation Z and Z' axes aligned (and maybe translation)
                    qinterp = QQuaternion::nlerp(q1,qtotv[kmol],i*dlt);
                    trs = (cmv[kmol]-cmantv[kmol]) * i * dlt;
                    for (int j = indmols[kmol] ; j < indmols[kmol+1]; j++){
                        if ((xyz[j]-xyzant[j]).lengthSquared() < 1.e-10){
                            xyzi.append(xyz[j]);
                        }
                        else{
                            xyzi.append(qinterp.rotatedVector(xyzant[j]+trs));
                        }
                    }
                    break;
                case 3:         // rotation Z and Z' axes opposed (and maybe translation)
                    qinterp = QQuaternion::nlerp(q1,qtotv[kmol],i*dlt);
                    trs = (cmv[kmol]-cmantv[kmol]) * i * dlt;
                    for (int j = indmols[kmol] ; j < indmols[kmol+1]; j++){
                        if ((xyz[j]-xyzant[j]).lengthSquared() < 1.e-10){
                            xyzi.append(xyz[j]);
                        }
                        else{
                            xyzi.append(qinterp.rotatedVector(xyzant[j]+trs));
                        }
                    }
                    break;
                case 4:
                    qinterp = QQuaternion::nlerp(q1,qtotv[kmol],i*dlt);
                    trs = (cmv[kmol]-cmantv[kmol]) * i * dlt;
                    for (int j = indmols[kmol] ; j < indmols[kmol+1]; j++){
                        if ((xyz[j]-xyzant[j]).lengthSquared() < 1.e-10){
                            xyzi.append(xyz[j]);
                        }
                        else{
                            xyzi.append(qinterp.rotatedVector(xyzant[j]-cmantv[kmol])+cmantv[kmol]+trs);
                        }
                    }
                    break;
                default:
                    qDebug() << "Case not expected: " << transformation[kmol];
            }
        }
        molecules->last()->loadxyz(xyzi);
        molecules->last()->makeStructureBondsandSticks();
        Enerinterp = Ener0 * (1 - dlt*i) + Ener1 * dlt * i;
        window->setepicenergy(Enerinterp);
//        qDebug() << "Enerinterp = " << Enerinterp;
        if (msrs)
            msrs->update_measure_centers(&xyzi);
        window->setframeknt(1);     // Dirty trick to guarantee that takes all frames
        updatedisplay();
        pause(delay);
    }
    Enerinterp = Ener1;
//    qDebug() << "Enerinterp = " << Enerinterp;
    molecules->last()->loadxyz(xyz);
    molecules->last()->makeStructureBondsandSticks();
    window->setepicenergy(Enerinterp);
    if (msrs)
        msrs->update_measure_centers(&xyz);
    window->setframeknt(1);     // Dirty trick to guarantee that takes all frames
    updatedisplay();
    pause(delay);
}
//  End of function quaterninterpolation
//  ---------------------------------------------------------------------------------------------------------------------------

//
// Slots for setting status bar when external code is launched
void glWidget::qmrun(QString a){
    if (a.isEmpty()){
        mainWin->statusBar()->setStyleSheet("color: black");
    }
    else{
        mainWin->statusBar()->setStyleSheet("color: blue");
    }
    mainWin->statusBar()->showMessage(a);
}

//  End slots for setting status bar when external code is launched
//  ---------------------------------------------------------------------------------------------------------------------------


//  Function recordoptim_changed: updates recording option
//
void glWidget::recordoptim_changed(bool a){
    recordoptim = a;
}

//  Function reloadxyz: reads geometry and charges of centers from a xyz file
//
bool glWidget::reloadxyz(QString filename){
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        return false;
    }
    QTextStream in(&file);
    QString line = in.readLine();
//        qDebug() << "line = " << line;
    if (line.isEmpty()){
        return false;
    }
    znuc.clear();
    xyz.clear();
    int ncen = line.toInt();
    if (ncen <= 0){
        return false;
    }
    line = in.readLine();   // line with energy
//    qDebug() << "line = " << line;
#if QT_VERSION < 0x050E00
    QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
    QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
    if (fields.count() < 2)
        return false;
    fields.takeFirst();
    Enerinterp = fields.takeFirst().toFloat();
    window->setepicenergy(Enerinterp);
    while (!in.atEnd()){
        line = in.readLine();
//        qDebug() << "line = " << line;
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (fields.count() == 1 && znuc.length() == ncen){
            ncen += fields.takeFirst().toInt();
        }
        if (fields.count() != 4){
            continue;
        }
        QString ba = fields.takeFirst().toLatin1();
        znuc << elem->getZsymbol(ba);
        double x = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
        double y = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
        double z = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
        xyz << QVector3D(x, y, z);

    }
//    qDebug() << "ncen = " << ncen;
//    qDebug() << "znuc.length() = " << znuc.length();
    file.close();
    if (znuc.length() != ncen){
        return false;
    }
    return true;
}
//  End of function reloadxyz
//  ---------------------------------------------------------------------------------------------------------------------------


//  Function replay_mespimization: animates minimzation
//
void glWidget::replay_mespimization(QString filename){
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("File %1 cannot be opened")
                            .arg(filename)+QString(":\n%1.").arg(file.errorString()));
        return;
    }
    if (molecules->length() == 0 || !molecules->last()->getiscluster()){
        if (!addemptyclusterTowindow(filename)){
            return;
        }
    }
    QDLwindow->mespimizerdialog->setBTNreset(false);
    QDLwindow->mespimizerdialog->setBTNreplay(false);
    mainWin->statusBar()->setStyleSheet("color: blue");
    mainWin->statusBar()->showMessage(tr("Animating..."));
    mespimirecordfilename = QDLwindow->mespimizerdialog->getrecordfilename();
    window->setrecordfilename(mespimirecordfilename);
    window->setrecord(QDLwindow->mespimizerdialog->getrecord());
    window->setrecordcommand(QDLwindow->mespimizerdialog->getrecordcommand());
    window->setremoveframes(QDLwindow->mespimizerdialog->getremoveframes());
    window->setdisplayEPIC(QDLwindow->mespimizerdialog->getdisplayenergy());
    ninterpol = QDLwindow->mespimizerdialog->getinterpolpoints();
    QTextStream in(&file);
    QString line = in.readLine();   // first line of frames file
    int numcaptures = line.toInt();
    numframes = (numcaptures-1)*(ninterpol+1)+1;
//    qDebug() << "number of original captures = " << numcaptures;
//    qDebug() << "Total number of frames = " << numframes;
    window->setnumframes(numframes);
    kntframes = 0;
    window->setframeknt(kntframes);
    znuc.clear();
    xyz.clear();
    if (QDLwindow->mespimizerdialog->getrecord()){
        QDLwindow->mespimizerdialog->startmakingmovie();
    }
    line = in.readLine();
    int ncen = line.toInt();        // second line of frames file
    if (ncen <= 0){
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Error in file %1")
                            .arg(filename)+QString(":\n%1.").arg(tr("wrong format")));
        return;
    }
    bool interpol = false;
    int knt = 0;
    indmols.clear();
    knt = 0;
    line = in.readLine();           // third line of frames file
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
    indmols << 0;
    int itop = fields.takeFirst().toInt();
    for (int i = 0 ; i < itop ; i++){
        indmols << knt + fields.takeFirst().toInt();
        knt = indmols.last();
    }
    line = in.readLine();   // fourth line of frames file
//    qDebug() << "line = " << line;
#if QT_VERSION < 0x050E00
        fields = line.split(' ',QString::SkipEmptyParts);
#else
        fields = line.split(' ',Qt::SkipEmptyParts);
#endif
    fields.takeFirst();
    Ener0 = fields.takeFirst().toFloat();
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        fields = line.split(' ',QString::SkipEmptyParts);
#else
        fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (fields.count() == 1){       // beginning line of new frame
            if (!interpol){
                molecules->last()->loadxyz(xyz);
                molecules->last()->loadcharges(znuc);
                molecules->last()->makeStructureBondsandSticks();
//                molecules->last()->initatomactive(false);
                if (msrs)
                    msrs->update_measure_centers(&xyz);
                xyz.clear();
                interpol = true;
                window->setframeknt(1);     // Dirty trick to guarantee that takes all frames
                updatedisplay();
                delay = QDLwindow->mespimizerdialog->getdelay();
//                qDebug() << "delay = " << delay;
                pause(delay);
                line = in.readLine();
                line = in.readLine();
//                qDebug() << "line = " << line;
#if QT_VERSION < 0x050E00
        fields = line.split(' ',QString::SkipEmptyParts);
#else
        fields = line.split(' ',Qt::SkipEmptyParts);
#endif
                fields.takeFirst();
                Ener1 = fields.takeFirst().toFloat();
            }
            else{
//                if (QDLwindow->mespimizerdialog->getlineinterpol()){
//                    linearinterpolation();
//                }
//                else{
//                    quaterninterpolation();
//                }
                quaterninterpolation();     // Forced to use quaternions interpolation
                xyz.clear();
                line = in.readLine();       // skips lines other than coordinates
                Ener0 = Ener1;
                line = in.readLine();
//                qDebug() << "line = " << line;
#if QT_VERSION < 0x050E00
        fields = line.split(' ',QString::SkipEmptyParts);
#else
        fields = line.split(' ',Qt::SkipEmptyParts);
#endif
                fields.takeFirst();
                Ener1 = fields.takeFirst().toFloat();
            }
        }
        else{
            QString ba = fields.takeFirst().toLatin1();
            znuc << elem->getZsymbol(ba);
            double x = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
            double y = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
            double z = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
            xyz << QVector3D(x, y, z);
        }
    }
//    if (QDLwindow->mespimizerdialog->getlineinterpol()){
//        linearinterpolation();
//    }
//    else{
//        quaterninterpolation();
//    }
    quaterninterpolation();     // Forced to use quaternions interpolation
    window->setframeknt(numframes-1);     // Dirty trick to guarantee that takes all frames
    updatedisplay();
    xyz.clear();
    file.close();
    QDLwindow->mespimizerdialog->setBTNreplay(true);
    QDLwindow->mespimizerdialog->setBTNreset(true);
    mainWin->statusBar()->setStyleSheet("color: red");
    mainWin->statusBar()->showMessage(tr("End of animation"));
    pause(1000);
    mainWin->statusBar()->clearMessage();
    mainWin->statusBar()->setStyleSheet("color: black");
}
//  End of function replay_mespimization
//  ---------------------------------------------------------------------------------------------------------------------------

//  Function reset_mespimization: loads first frame of animation
//
void glWidget::reset_mespimization(QString filename){
    QString fileaux = QFileInfo(filename).canonicalPath()+"/"+QFileInfo(filename).completeBaseName();
    if (!molecules->isEmpty() && fileaux.trimmed() != molecules->last()->getfullname().trimmed()){
        if (molecules->last()->getiscluster()){
            deletecluster();
        }
        QString fileinit = filename;
        fileinit.replace("frames","init");
        addclusterTowindow(fileinit);
        return;
    }
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("File %1 cannot be opened")
                            .arg(filename)+QString(":\n%1.").arg(file.errorString()));
        return;
    }
    if (molecules->length() == 0 || !molecules->last()->getiscluster()){
        if (!addemptyclusterTowindow(filename)){
            return;
        }
    }
    QTextStream in(&file);
    QString line = in.readLine();   // first line of frames file
    in.readLine();  // second line of frames file
    znuc.clear();
    xyz.clear();
    int ncen;
    if ((ncen = line.toInt()) <= 0){
        QMessageBox::warning(this, tr("MESPIMIZER"),tr("Error in file %1")
                            .arg(filename)+QString(":\n%1.").arg(tr("wrong format")));
        return;
    }
    int knt = 0;
    indmols.clear();
    knt = 0;
    line = in.readLine();    // third line of frames file
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
    indmols << 0;
    int itop = fields.takeFirst().toInt();
    for (int i = 0 ; i < itop ; i++){
        indmols << knt + fields.takeFirst().toInt();
        knt = indmols.last();
    }
    line = in.readLine();   // fourth line of frames file
//    qDebug() << "line = " << line;
#if QT_VERSION < 0x050E00
    fields = line.split(' ',QString::SkipEmptyParts);
#else
    fields = line.split(' ',Qt::SkipEmptyParts);
#endif
    fields.takeFirst();
    Enerinterp = fields.takeFirst().toFloat();
    window->setepicenergy(Enerinterp);
    window->setdisplayEPIC(QDLwindow->mespimizerdialog->getdisplayenergy());
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        fields = line.split(' ',QString::SkipEmptyParts);
#else
        fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (fields.count() == 1){       // beginning line of frame
            molecules->last()->loadxyz(xyz);
            molecules->last()->loadcharges(znuc);
            molecules->last()->makeStructureBondsandSticks();
            if (msrs)
                msrs->update_measure_centers(&xyz);
            xyz.clear();
            file.close();
            updatedisplay();
            return;
        }
        else{
            QString ba = fields.takeFirst().toLatin1();
            znuc << elem->getZsymbol(ba);
            double x = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
            double y = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
            double z = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
            xyz << QVector3D(x, y, z);
        }
    }
//    molecules->last()->initatomactive(false);
}
//  End of function reset_mespimization
//  ---------------------------------------------------------------------------------------------------------------------------

//  Function setdisplayEPIC: updates EPIC energy display flag
//
void glWidget::setdisplayEPIC(bool a){
    displayEPIC = a;
    window->setdisplayEPIC(displayEPIC);
}

//  Function setenergycolor: updates EPIC energy color
//
void glWidget::setenergycolor(QColor a){
    energycolor = a;
    window->setenergycolor(energycolor);
}

//  Function setenergyfont: updates EPIC energy font
//
void glWidget::setenergyfont(QFont a){
    energyfont = a;
    window->setenergyfont(energyfont);
}

//  Function setenergyprecision: updates EPIC energy precision for display
//
void glWidget::setenergyprecision(int a){
    energyprecision = a;
    window->setenergyprecision(energyprecision);
}

//  Function setenergyprecision: updates EPIC energy units for display
//
void glWidget::sethartree(bool a){
    hartree = a;
    window->sethartree(hartree);
}

//  Function speed_changed: updates animation speed
//
void glWidget::speed_changed(int i){
    speed = i;
}

//  Function template_changed: updates molecule index for template
//
void glWidget::template_changed(int i){
    templateindex = i;
}
//********************************************************************************
//          Save geometry SLOTS
//********************************************************************************


void glWidget::BTNsaveGeometry_clicked(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName;
    if (molecules->count() == 1){
        fileName = filedialog.getSaveFileName(this,tr("Open file ..."),molpath,
            tr("Settings saving files") + " *.xyz" + " (*.xyz);;"+tr("All files")+" (*)");
        if (!fileName.contains(".xyz"))
            fileName.append(".xyz");
    }
    else if (molecules->count() == 0){
        QMessageBox::warning(this, tr("DAMQT"),tr("No molecule available to save geometry"));
        return;
    }
    else{       // Geometry of a cluster of molecules is to be saved
        fileName = filedialog.getSaveFileName(this,tr("Open file ..."),molpath,
            tr("Settings saving files") + " *-cluster.xyz" + " (*-cluster.xyz);;"+tr("All files")+" (*)");
        if (!fileName.contains("-cluster.xyz"))
            fileName.append("-cluster.xyz");
    }
    QFile fileout(fileName);
    if (!fileout.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be opened")
                            .arg(fileName)+QString(":\n%1.").arg(fileout.errorString()));
        emit moveToTop(windownumber);
        return;
    };
    if (!fileout.isOpen()){
        fileout.open(QIODevice::Text | QIODevice::ReadWrite);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    bool ishost = true;
    Elements *elem = new Elements();
    for (int i = 0 ; i < molecules->count() ; i++){
        QVector <int> znuc;
        QVector <QVector3D> xyz;
        znuc =  molecules->at(i)->getcharges();
        xyz = molecules->at(i)->getxyz();
        QMatrix4x4 m;

        m.rotate(molecules->at(0)->getrotation().conjugated());
        m.translate(-molecules->at(0)->gettranslation());
        m.translate(molecules->at(i)->gettranslation());
        m.rotate(molecules->at(i)->getrotation());

        buff.append(QString("      %1\n").arg(molecules->at(i)->getnumatoms()));
        if (molecules->count() == 1){
            buff.append(QString("      %1\n").arg(molecules->at(0)->getname()));
        }
        else{
            if (ishost){
                buff.append(QString("      Host\n"));
                ishost = false;
            }
            else{
                buff.append(QString("      Guest %1\n").arg(i));
            }
        }

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
    emit moveToTop(windownumber);
    delete elem;
}


//********************************************************************************
//          Save/Retrieve settings SLOTS
//********************************************************************************

void glWidget::BTNsettings_clicked(){
        FRMsettings->setVisible(!FRMsettings->isVisible());
        if (FRMsettings->isVisible()){
            BTNsettings->setStyleSheet("QPushButton {color: red;}");
            BTNsettings->setText(tr("Hide save/retrieve settings"));
        }
        else{
            BTNsettings->setStyleSheet("QPushButton {color: black;}");
            BTNsettings->setText(tr("Save/retrieve settings"));
        }
        QDLwindow->update();
        QDLwindow->adjustSize();
        QDLwindow->renewwidth(gdock->size());
}

void glWidget::BTNretrieveSettings_clicked(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),molpath, tr("Settings saving files")
                + " *.3Dsettings" + " (*.3Dsettings);;"+tr("All files")+" (*)");
    QFile filein(fileName);
    if (!filein.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                    .arg(fileName)+QString(":\n%1.").arg(filein.errorString()));
            emit moveToTop(windownumber);
            return;
    };
    CaptureFolder = QFileInfo(fileName).path();
    std::string file = QString(fileName).toStdString();
    QString qv;
    QStringList qlv;
    std::string v;
    if (CHKsettingsmolecs->isChecked()){
        int nummolecs = QString::fromStdString(CIniFile::GetValue("nummolecs","V3DMOLECULES",file)).toInt();
        for (int i = 0 ; i < nummolecs ; i++){
            std::string section = QString("V3DMOLECULE_%1").arg(i).toStdString();
            QString filename = QString::fromStdString(CIniFile::GetValue("name",section,file));
            if (QFileInfo(filename).suffix() == "xyz"){
                if (readxyz(filename))
                    retrievemolecule();
                else
                    return;
            }
            else if(QFileInfo(filename).suffix() == "ggbs"){
                if (readggbs(filename))
                    retrievemolecule();
                else
                    return;
            }
            QFont font;
            font.fromString(QString::fromStdString(CIniFile::GetValue("font",section,file)));
            molecules->last()->setfont(font);
            bool showsymbols = QString::fromStdString(CIniFile::GetValue("showsymbols",section,file)).toInt();
            molecules->last()->setdrawatomsymbols(showsymbols);
            bool showindices = QString::fromStdString(CIniFile::GetValue("showindices",section,file)).toInt();
            molecules->last()->setdrawatomindices(showindices);
            bool showselect = QString::fromStdString(CIniFile::GetValue("showselect",section,file)).toInt();
            molecules->last()->setonlyatomactive(showselect);
            if (showselect){
#if QT_VERSION < 0x050E00
                qlv = QString::fromStdString(CIniFile::GetValue("atomactive",section,file))
                        .split(",",QString::SkipEmptyParts);
#else
                qlv = QString::fromStdString(CIniFile::GetValue("atomactive",section,file))
                        .split(",",Qt::SkipEmptyParts);
#endif
                for (int i = 0 ; i < qlv.length() ; i++){
                    molecules->last()->setatomactive(i,qlv.at(i).toInt());
                }
            }
            int vshift = QString::fromStdString(CIniFile::GetValue("vshift",section,file)).toInt();
            molecules->last()->setlabelsvshift(vshift);
            if (CHKsettingspos->isChecked()){
#if QT_VERSION < 0x050E00
                qlv = QString::fromStdString(CIniFile::GetValue("position",section,file))
                        .split(",",QString::SkipEmptyParts);
#else
                qlv = QString::fromStdString(CIniFile::GetValue("position",section,file))
                        .split(",",Qt::SkipEmptyParts);
#endif
                if (qlv.length() == 3){
                    molecules->last()->settranslation(QVector3D(qlv[0].toFloat(),qlv[1].toFloat(),qlv[2].toFloat()));
                    molecules->last()->settranslationButtons();
                }
            }
            if (CHKsettingsrot->isChecked()){
#if QT_VERSION < 0x050E00
                qlv = QString::fromStdString(CIniFile::GetValue("rotation",section,file))
                        .split(",",QString::SkipEmptyParts);
#else
                qlv = QString::fromStdString(CIniFile::GetValue("rotation",section,file))
                        .split(",",Qt::SkipEmptyParts);
#endif
                if (qlv.length() == 4){
                    QQuaternion rotation = QQuaternion(qlv[0].toFloat(),qlv[1].toFloat(),qlv[2].toFloat(),qlv[3].toFloat());
                    molecules->last()->setrotation(rotation.normalized());
                    molecules->last()->setrotationButtons();
                }
            }
#if QT_VERSION < 0x050E00
            qlv = QString::fromStdString(CIniFile::GetValue("fontcolor",section,file))
                    .split(",",QString::SkipEmptyParts);
#else
            qlv = QString::fromStdString(CIniFile::GetValue("fontcolor",section,file))
                    .split(",",Qt::SkipEmptyParts);
#endif
            if (qlv.length() == 3){
                molecules->last()->setfontcolor(QColor(qlv.at(0).toInt(),qlv.at(1).toInt(),qlv.at(2).toInt()));
            }
            if (CHKsettingssurf->isChecked()){
                int nsurf = QString::fromStdString(CIniFile::GetValue("nsurf",section,file)).toInt();
                for (int j = 0 ; j < nsurf ; j++){
                    std::string surfname = QString("surfname_%1").arg(j).toStdString();
                    QString surfilename = QString::fromStdString(CIniFile::GetValue(surfname,section,file));
                    if (molecules->last()->retrievesurf(surfilename)){
                        connect (molecules->last()->surfaces->last(), SIGNAL(updatedisplay()), molecules->last(), SIGNAL(updatedisplay()));
                        if (molecules->last()->QDLeditMolecule_isVisible())
                            molecules->last()->updateeditMoleculeDialog();
                        std::string solid = QString("solid_%1").arg(j).toStdString();
                        bool surfsolid = QString::fromStdString(CIniFile::GetValue(solid,section,file)).toInt();
                        molecules->last()->surfaces->last()->setsolidsurf(surfsolid);
                        std::string trasluc = QString("trasluc_%1").arg(j).toStdString();
                        bool surftrasluc = QString::fromStdString(CIniFile::GetValue(trasluc,section,file)).toInt();
                        molecules->last()->surfaces->last()->settranslucence(surftrasluc);
                        std::string opacity = QString("opacity_%1").arg(j).toStdString();
                        float surfopacity =  QString::fromStdString(CIniFile::GetValue(opacity,section,file)).toFloat();
                        molecules->last()->surfaces->last()->setopacity(surfopacity);
                        std::string topcolor = QString("topcolor_%1").arg(j).toStdString();
                        float surftopcolor =  QString::fromStdString(CIniFile::GetValue(topcolor,section,file)).toFloat();
                        molecules->last()->surfaces->last()->settopcolor(surftopcolor);
                        molecules->last()->surfaces->last()->resetsurface();
                    }
                }
                molecules->last()->emitupdatedisplay();
            }
            if (CHKsettingsgrids->isChecked()){
                int ngrids = QString::fromStdString(CIniFile::GetValue("ngrids",section,file)).toInt();
                for (int j = 0 ; j < ngrids ; j++){
                    std::string grid = QString("gridname_%1").arg(j).toStdString();
                    QString gridname = QString::fromStdString(CIniFile::GetValue(grid,section,file));
                    if (molecules->last()->retrievegrid(gridname)){
                        std::string nisosurf = QString("nisosurf_%1").arg(j).toStdString();
                        int numisosurf = QString::fromStdString(CIniFile::GetValue(nisosurf,section,file)).toInt();
                        std::string maxcontour = QString("maxcontourval_%1").arg(j).toStdString();
                        float maxcontourval =  QString::fromStdString(CIniFile::GetValue(maxcontour,section,file)).toFloat();
                        std::string mincontour = QString("mincontourval_%1").arg(j).toStdString();
                        float mincontourval =  QString::fromStdString(CIniFile::GetValue(mincontour,section,file)).toFloat();
                        molecules->last()->grids->at(j)->setmaxcontourvalue(maxcontourval);
                        molecules->last()->grids->at(j)->setmincontourvalue(mincontourval);
                        molecules->last()->grids->at(j)->set_ProjectFolder(QFileInfo(gridname).absolutePath());
                        molecules->last()->grids->at(j)->set_ProjectName(QFileInfo(gridname).fileName().remove("-deform").remove("-d.plt").remove("v-plt"));
                        for (int k = 0 ; k < numisosurf ; k++){
                            std::string contour = QString("contourval_%1_%2").arg(j).arg(k).toStdString();
                            float contourval =  QString::fromStdString(CIniFile::GetValue(contour,section,file)).toFloat();
                            std::string solid = QString("isosurfsolid_%1_%2").arg(j).arg(k).toStdString();
                            bool solidsurf = QString::fromStdString(CIniFile::GetValue(solid,section,file)).toInt();
                            std::string trasluc = QString("isosurftrasluc_%1_%2").arg(j).arg(k).toStdString();
                            bool surftrasluc = QString::fromStdString(CIniFile::GetValue(trasluc,section,file)).toInt();
                            std::string opacity = QString("isosurfopacity_%1_%2").arg(j).arg(k).toStdString();
                            float surfopacity =  QString::fromStdString(CIniFile::GetValue(opacity,section,file)).toFloat();
                            std::string color = QString("isosurfcolor_%1_%2").arg(j).arg(k).toStdString();
                            QStringList colorlist = QString::fromStdString(CIniFile::GetValue(color,section,file)).split(",");
                            molecules->last()->grids->at(j)->addisosurf();
                            molecules->last()->grids->at(j)->surfaces->last()->setcontourvalue(contourval);
                            molecules->last()->grids->at(j)->surfaces->last()->setsolidsurf(solidsurf);
                            molecules->last()->grids->at(j)->surfaces->last()->settranslucence(surftrasluc);
                            molecules->last()->grids->at(j)->surfaces->last()->setopacity(surfopacity);
                            if (colorlist.length() == 3)
                                molecules->last()->grids->at(j)->surfaces->last()->setsurfcolor(QColor(colorlist[0].toInt(),
                                                colorlist[1].toInt(),colorlist[2].toInt()));
                            molecules->last()->grids->at(j)->generatesurf(k);
                        }
                    }
                }
            }
        }
    }
    if (CHKsettingslights->isChecked()){
#if QT_VERSION < 0x050E00
        qlv = QString::fromStdString(CIniFile::GetValue("ambientcolor","V3DLIGHTS",file))
                .split(",",QString::SkipEmptyParts);
#else
        qlv = QString::fromStdString(CIniFile::GetValue("ambientcolor","V3DLIGHTS",file))
                .split(",",Qt::SkipEmptyParts);
#endif
        if (qlv.length() == 3){
            ambientcolor = QColor(qlv.at(0).toInt(),qlv.at(1).toInt(),qlv.at(2).toInt());
            QDLwindow->lightsdialog->setambientcolor(ambientcolor);
            if (window)
                window->setAmbientColor(ambientcolor);
        }
#if QT_VERSION < 0x050E00
        qlv = QString::fromStdString(CIniFile::GetValue("lightcolor","V3DLIGHTS",file))
                .split(",",QString::SkipEmptyParts);
#else
        qlv = QString::fromStdString(CIniFile::GetValue("lightcolor","V3DLIGHTS",file))
                .split(",",Qt::SkipEmptyParts);
#endif
        if (qlv.length() == 3){
            lightcolor = QColor(qlv.at(0).toInt(),qlv.at(1).toInt(),qlv.at(2).toInt());
            QDLwindow->lightsdialog->setlightcolor(lightcolor);
            if (window)
                window->setLightColor(lightcolor);
        }
#if QT_VERSION < 0x050E00
        qlv = QString::fromStdString(CIniFile::GetValue("specularcolor","V3DLIGHTS",file))
                .split(",",QString::SkipEmptyParts);
#else
        qlv = QString::fromStdString(CIniFile::GetValue("specularcolor","V3DLIGHTS",file))
                .split(",",Qt::SkipEmptyParts);
#endif
        if (qlv.length() == 3){
            specularcolor = QColor(qlv.at(0).toInt(),qlv.at(1).toInt(),qlv.at(2).toInt());
            QDLwindow->lightsdialog->setspecularcolor(specularcolor);
            if (window)
                window->setSpecularColor(specularcolor);
        }
    }
    if (CHKsettingsbkg->isChecked()){
#if QT_VERSION < 0x050E00
        qlv = QString::fromStdString(CIniFile::GetValue("bkgcolor","V3DBKG",file))
                .split(",",QString::SkipEmptyParts);
#else
        qlv = QString::fromStdString(CIniFile::GetValue("bkgcolor","V3DBKG",file))
                .split(",",Qt::SkipEmptyParts);
#endif
        if (qlv.length() == 3){
            bkgcolor = QColor(qlv.at(0).toInt(),qlv.at(1).toInt(),qlv.at(2).toInt());
            QDLwindow->lightsdialog->setbkgcolor(bkgcolor);
            if (window)
                window->setBkgColor(bkgcolor);
        }
    }
    if (CHKsettingsballs->isChecked()){
        ballradius =  QString::fromStdString(CIniFile::GetValue("ballradius","V3DBALLS",file)).toFloat();
        cylradius =  QString::fromStdString(CIniFile::GetValue("cylradius","V3DBALLS",file)).toFloat();
        if (ballradius > 0)
            QDLwindow->ballsandcylsdialog->setballrad(ballradius);
        if (cylradius > 0)
            QDLwindow->ballsandcylsdialog->setcylrad(cylradius);
        ballcyl_changed();
    }
    if (CHKsettingsviewport->isChecked()){
        float aux;
        aux = QString::fromStdString(CIniFile::GetValue("fov","V3DVIEWPORT",file)).toFloat();
        if (aux > 0)
            fov =  aux;
        aux =  QString::fromStdString(CIniFile::GetValue("zNear","V3DVIEWPORT",file)).toFloat();
        if (aux > 0)
            zNear =  aux;
        aux =  QString::fromStdString(CIniFile::GetValue("zFar","V3DVIEWPORT",file)).toFloat();
        if (aux > zNear)
            zFar =  aux;
        updateviewport();
    }
    updatedisplay();
    emit moveToTop(windownumber);
}

void glWidget::BTNsaveSettings_clicked(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getSaveFileName(this,tr("Open file ..."),molpath,
        tr("Settings saving files") + " *.3Dsettings" + " (*.3Dsettings);;"+tr("All files")+" (*)");
    if (!fileName.contains(".3Dsettings"))
        fileName.append(".3Dsettings");
    QFile fileout(fileName);
    if (!fileout.open(QFile::ReadWrite | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be opened")
                            .arg(fileName)+QString(":\n%1.").arg(fileout.errorString()));
        emit moveToTop(windownumber);
        return;
    };
    if (!fileout.isOpen()){
        fileout.open(QIODevice::Text | QIODevice::ReadWrite);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    if (CHKsettingslights->isChecked()){
        buff.append("[V3DLIGHTS]\n");
        buff.append("ambientcolor="+QString("%1,%2,%3\n").arg(ambientcolor.red())
                    .arg(ambientcolor.green()).arg(ambientcolor.blue()));
        buff.append("lightcolor="+QString("%1,%2,%3\n").arg(lightcolor.red())
                    .arg(lightcolor.green()).arg(lightcolor.blue()));
        buff.append("specularcolor="+QString("%1,%2,%3\n").arg(specularcolor.red())
                    .arg(specularcolor.green()).arg(specularcolor.blue()));
    }
    if (CHKsettingsbkg->isChecked()){
        buff.append("[V3DBKG]\n");
        buff.append("bkgcolor="+QString("%1,%2,%3\n").arg(bkgcolor.red())
                    .arg(bkgcolor.green()).arg(bkgcolor.blue()));
    }
    if (CHKsettingsballs->isChecked()){
        buff.append("[V3DBALLS]\n");
        buff.append("ballradius="+QString("%1\n").arg(ballradius));
        buff.append("cylradius="+QString("%1\n").arg(cylradius));
    }
    if (CHKsettingsviewport->isChecked()){
        buff.append("[V3DVIEWPORT]\n");
        buff.append("fov="+QString("%1\n").arg(fov));
        buff.append("zFar="+QString("%1\n").arg(zFar));
        buff.append("zNear="+QString("%1\n").arg(zNear));
    }
    if (CHKsettingsmolecs->isChecked()){
        buff.append("[V3DMOLECULES]\n");
        int kont = 0;
        for (int i = 0 ; i < molecules->count() ; i++){
            if (molecules->at(i)->getiscluster())
                continue;
            kont++;
        }
        buff.append("nummolecs="+QString("%1\n").arg(kont));
        for (int i = 0 ; i < molecules->count() ; i++){
            if (molecules->at(i)->getiscluster())
                continue;
            buff.append(QString("[V3DMOLECULE_%1]\n").arg(i));
            buff.append("name="+molecules->at(i)->getfullname()+"\n");
            buff.append("font="+molecules->at(i)->getfont().toString()+"\n");
            buff.append("showsymbols="+QString("%1\n").arg(molecules->at(i)->getdrawatomsymbols()));
            buff.append("showindices="+QString("%1\n").arg(molecules->at(i)->getdrawatomindices()));
            buff.append("showselect="+QString("%1\n").arg(molecules->at(i)->getonlyatomactive()));
            if (molecules->at(i)->getonlyatomactive()){
                buff.append("atomactlen="+QString("%1\n").arg(molecules->at(i)->getnumatoms()));
                buff.append("atomactive=");
                for (int j = 0 ; j < molecules->at(i)->getnumatoms() ; j++ ){
                    buff.append(QString("%1,").arg(molecules->at(i)->isatomactive(j)));
                }
                buff.replace(buff.count()-1,1,"\n");
            }
            buff.append("vshift="+QString("%1\n").arg(molecules->at(i)->getlabelsvshift()));
            if (CHKsettingspos->isChecked()){
                QVector3D translation=molecules->at(i)->gettranslation();
                buff.append("position="+QString("%1,%2,%3\n").arg(translation.x()).arg(translation.y())
                            .arg(translation.z()));
            }
            if (CHKsettingsrot->isChecked()){
                QQuaternion rotation=molecules->at(i)->getrotation();
                buff.append("rotation="+QString("%1,%2,%3,%4\n").arg(rotation.scalar()).arg(rotation.x())
                            .arg(rotation.y()).arg(rotation.z()));
            }
            QColor fontcolor=molecules->at(i)->getfontcolor();
            buff.append("fontcolor="+QString("%1,%2,%3\n").arg(fontcolor.red())
                        .arg(fontcolor.green()).arg(fontcolor.blue()));
            if (CHKsettingssurf->isChecked()){
                buff.append("nsurf="+QString("%1\n").arg(molecules->at(i)->surfaces->count()));
                for (int j = 0 ; j < molecules->at(i)->surfaces->count() ; j++){
                    buff.append(QString("surfname_%1=").arg(j)+molecules->at(i)->surfaces->at(j)->getfullname()+"\n");
                    buff.append(QString("solid_%1=").arg(j)+QString("%1\n")
                                .arg(molecules->at(i)->surfaces->at(j)->getsolidsurf()));
                    buff.append(QString("opacity_%1=%2\n").arg(j)
                                .arg(molecules->at(i)->surfaces->at(j)->getopacity()));
                    buff.append(QString("trasluc_%1=").arg(j)+QString("%1\n")
                                .arg(molecules->at(i)->surfaces->at(j)->gettranslucence()));
                    buff.append(QString("topcolor_%1=%2\n").arg(j)
                                .arg(molecules->at(i)->surfaces->at(j)->gettopcolor()));
                }
            }
            if (CHKsettingsgrids->isChecked()){
                buff.append("ngrids="+QString("%1\n").arg(molecules->at(i)->grids->count()));
                for (int j = 0 ; j < molecules->at(i)->grids->count() ; j++){
                    buff.append(QString("gridname_%1=").arg(j)+molecules->at(i)->grids->at(j)->getfullname()+"\n");
                    if (CHKsettingsisosurf->isChecked()){
                        buff.append(QString("nisosurf_%1=").arg(j)+QString("%1\n")
                                    .arg(molecules->at(i)->grids->at(j)->surfaces->count()));
                        buff.append(QString("maxcontourval_%1=%2\n").arg(j)
                                .arg(molecules->at(i)->grids->at(j)->getmaxcontourvalue()));
                        buff.append(QString("mincontourval_%1=%2\n").arg(j)
                                .arg(molecules->at(i)->grids->at(j)->getmincontourvalue()));
                        for (int k = 0 ; k < molecules->at(i)->grids->at(j)->surfaces->count() ; k++){
                            buff.append(QString("contourval_%1_%2=%3\n").arg(j).arg(k)
                                    .arg(molecules->at(i)->grids->at(j)->surfaces->at(k)->getcontourvalue()));

                            buff.append(QString("isosurfsolid_%1_%2=%3\n").arg(j).arg(k)
                                    .arg(molecules->at(i)->grids->at(j)->surfaces->at(k)->getsolidsurf()));
                            buff.append(QString("isosurfopacity_%1_%2=%3\n").arg(j).arg(k)
                                    .arg(molecules->at(i)->grids->at(j)->surfaces->at(k)->getopacity()));
                            buff.append(QString("isosurftrasluc_%1_%2=%3\n").arg(j).arg(k)
                                    .arg(molecules->at(i)->grids->at(j)->surfaces->at(k)->gettranslucence()));
                            QColor surfcolor=molecules->at(i)->grids->at(j)->surfaces->at(k)->getsurfcolor();
                            buff.append(QString("isosurfcolor_%1_%2=%3,%4,%5\n").arg(j).arg(k)
                                        .arg(surfcolor.red()).arg(surfcolor.green()).arg(surfcolor.blue()));
                        }
                    }
                }
            }
        }
    }
#if QT_VERSION < 0x050E00
    outfile << buff << endl;
#else
    outfile << buff << Qt::endl;
#endif
    fileout.close();
    emit moveToTop(windownumber);
}

void glWidget::CHKsettingsgrids_changed(){
    if(CHKsettingsgrids->isChecked()){
        CHKsettingsisosurf->setEnabled(true);
    }
    else{
        CHKsettingsisosurf->setEnabled(false);
    }
}

void glWidget::CHKsettingsmolecs_changed(){
    if(CHKsettingsmolecs->isChecked()){
        CHKsettingspos->setEnabled(true);
        CHKsettingsrot->setEnabled(true);
        CHKsettingssurf->setEnabled(true);
        CHKsettingsgrids->setEnabled(true);
        if(CHKsettingsgrids->isChecked()){
            CHKsettingsisosurf->setEnabled(true);
        }
    }
    else{
        CHKsettingspos->setEnabled(false);
        CHKsettingsrot->setEnabled(false);
        CHKsettingssurf->setEnabled(false);
        CHKsettingsgrids->setEnabled(false);
        CHKsettingsisosurf->setEnabled(false);
    }
}

//********************************************************************************
//          Window functions
//********************************************************************************

int glWidget::getwindownumber(){
    return windownumber;
}

void glWidget::deleteWindowDialog(){
    delete QDLwindow;
    QDLwindow = nullpointer;
}

void glWidget::emitmovetotop(){
    window->update();
    emit moveToTop(windownumber);
}

void glWidget::updateWindowDialog(QDockWidget *gdock){
    bool QDLwindowexist = false;
    if (QDLwindow){
        QDLwindowexist = true;
        delete QDLwindow;
        QDLwindow = nullpointer;
        delete scrollArea;
        scrollArea = nullpointer;
    }
    if (QDLwindowexist){
        createWindowDialog(gdock);
        QDLwindow->renewwidth(gdock->size());
        scrollArea = new myScrollArea();
        scrollArea->setWidget(QDLwindow);
        scrollArea->updatesize(QDLwindow->size());
        gdock->setWidget(scrollArea);
    }
}

//********************************************************************************
//          Lights functions
//********************************************************************************

QColor glWidget::getbkgcolor(){
    return bkgcolor;
}


void glWidget::lower_mainwindow(){
    if (mainWin){
        mainWin->lower();
    }
}

void glWidget::raise_mainwindow(){
    if (mainWin){
        mainWin->raise();
    }
}

void glWidget::setbkgcolor(QColor a){
    bkgcolor = a;
}


void glWidget::set_ProjectFolder(QString name){
    ProjectFolder = name;
    CaptureFolder = ProjectFolder;
    recordfilename = ProjectFolder + "/film";
}

void glWidget::set_ProjectName(QString name){
    ProjectName = name;
}

//  Function setWindowName: sets window name
//
void glWidget::setWindowName(QString name){
    windowname = name;
}

//  Function getWindowName: retrieves window name
//
QString glWidget::getWindowName(){
    return windowname;
}

//********************************************************************************
//          Activate molecule function
//********************************************************************************

void glWidget::activatemolecules(){
    for (int i = 0 ; i < molecules->count() ; i++){
        if (activateboxes[i]->isChecked()){
            molecules->at(i)->setactive(true);
        }
        else{
            molecules->at(i)->setactive(false);
        }
    }
    window->update();
    emit moveToTop(windownumber);
}


void glWidget::checkactivate(int i){
    if (activateboxes.length() > i)
        activateboxes.at(i)->setChecked(! activateboxes.at(i)->isChecked());
}

//********************************************************************************
//          Loading data functions and related functions
//********************************************************************************

//  Function addmolecule: adds a new molecule
//
bool glWidget::addmolecule(){
    if (window)
        window->makeCurrent();
    molecules->append(new molecule());
    molecules->last()->set_ProjectFolder(ProjectFolder);
    molecules->last()->set_ProjectName(ProjectName);
    if (!loadmoleculargeometry()){
        molecules->removeLast();
        return false;
    }
    else{
        molecules->last()->loadcharges(znuc);
        molecules->last()->loadxyz(xyz);
        molecules->last()->setpath(molpath);
//        molecules->last()->setrotation(rotation.normalized());
        molecules->last()->setvisible(true);
        molecules->last()->initatomactive(false);
        molecules->last()->set_ProjectFolder(ProjectFolder);
        molecules->last()->set_ProjectName(ProjectName);
        if (window)
            window->resetlabaxes();
        connections << connect(molecules->last(), SIGNAL(updatedisplay()), this, SLOT(updatedisplay()));
        connections << connect(molecules->last(), SIGNAL(updateGL()), this, SLOT(updateGL()));
        connections << connect(molecules->last(), SIGNAL(animate(bool)), this, SLOT(animate(bool)));
        QDLwindow->mespimizerdialog->setSPBhostmax(molecules->length());
        QDLwindow->mespimizerdialog->setSPBtemplatemax(molecules->length());
        if (molecules->length() > 1 && QDLwindow->mespimizerdialog->gettemplate() == 1){
            templateindex = 2;
            QDLwindow->mespimizerdialog->setSPBtemplate(templateindex);
        }
        return true;
    }
}
//  End of function addmolecule
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function addmoleculeTowindow: adds a molecule to the window
//
void glWidget::addmoleculeTowindow(){
    if (!molecules->isEmpty() && molecules->last()->getiscluster()){
        QMessageBox msgBox;
        msgBox.setText(tr("addmoleculeTowindow"));
        msgBox.setInformativeText(tr("last system is a cluster.\n"
            "You must delete it before adding a new molecule."));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    if (addmolecule()){
        molecules->last()->setname(molname);
        molecules->last()->setfullname(fullmolname);
        molecules->last()->setballradius(ballradius);
        molecules->last()->setcylradius(cylradius);
        molecules->last()->setdisthressq(disthressq);
        molecules->last()->setscaleradii(scaleradii);
        molecules->last()->makeStructureBondsandSticks();
        molecules->last()->setvisible(true);
        molecules->last()->set_ProjectFolder(ProjectFolder);
        molecules->last()->set_ProjectName(ProjectName);
        molecule_names.append(molname);
        if (msrs){
            msrs->set_molecules(molecule_names);
        }
        updatedisplay();
        if (QDLwindow){
            updateWindowDialog(gdock);
        }
        this->moveToTop(windownumber);
    }
}
//  End of function addmoleculeTowindow
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------


void glWidget::animate(bool a){
    animating = a;
}

//  Function loadmoleculargeometry: loads geometry and charges of centers
//
bool glWidget::loadmoleculargeometry(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    if (!ProjectFolder.isEmpty())
        filedialog.setDirectory(ProjectFolder);
    QString filename = filedialog.getOpenFileName(this,tr("Open geometry file ..."),ProjectFolder, tr("Geometry and basis set files")
                + " *.ggbs, *.xyz " + " (*.ggbs *.xyz);;"+tr("All files")+" (*)");
    if (filename.length()==0){
        return false;
    }
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf=="ggbs"){
        if (!readggbs(filename)){
            QMessageBox msgBox;
            msgBox.setText(tr("loadmoleculargeometry"));
            msgBox.setInformativeText(tr("Error reading file %1").arg(filename));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return false;
        }
        else{
            return true;
        }
    }else if (suf=="xyz"){
        if (!readxyz(filename)){
            QMessageBox msgBox;
            msgBox.setText(tr("loadmoleculargeometry"));
            msgBox.setInformativeText(tr("Error reading file %1").arg(filename));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return false;
        }
        else{
            return true;
        }
    }else{
        QMessageBox msgBox;
        msgBox.setText(tr("loadmoleculargeometry"));
        msgBox.setInformativeText(tr("Extension of file %1 must be ggbs or xyz").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
}
//  End of function loadmoleculargeometry
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function readggbs: reads geometry and charges of centers from a ggbs file
//
bool glWidget::readggbs(QString filename){
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox msgBox;
        msgBox.setText(tr("readggbs"));
        msgBox.setInformativeText(tr("File %1 cannot be read").arg(filename)+QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        return false;
    }
    QTextStream in(&file);
    QString line = in.readLine();
    znuc.clear();
    xyz.clear();
    int ncen;
    if ((ncen = line.toInt()) <= 0){
        QMessageBox msgBox;
        msgBox.setText(tr("readggbs"));
        msgBox.setInformativeText(tr("First line of file %1 does not contain the number of centers").arg(filename)
                                  +QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        return false;
    }
    QVector3D xyzC = QVector3D(0.,0.,0.);
    double q;
    double qtot = 0.;
    while (!in.atEnd() && znuc.length() < ncen){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (fields.count() != 4) continue;
//        xyz << QVector3D(fields.takeFirst().toDouble(),fields.takeFirst().toDouble(),fields.takeFirst().toDouble());
        double x = fields.takeFirst().toDouble();
        double y = fields.takeFirst().toDouble();
        double z = fields.takeFirst().toDouble();
        xyz << QVector3D(x, y, z);
        q = fields.takeFirst().toDouble();
        znuc << (int) q ;
        xyzC += xyz.last() * q;
        qtot += q;
    }
    xyzC /= qtot;
    for (int i = 0 ; i < xyz.length() ; i++){
       xyz.replace(i,xyz.at(i)-xyzC);
    }
    QFileInfo fileInfo(filename);
    fullmolname = QString(fileInfo.absoluteFilePath());
    molname = QString(fileInfo.baseName());
    molpath = QString(fileInfo.absolutePath());
    if (znuc.length() != ncen){
        if ((ncen = line.toInt()) <= 0){
            QMessageBox msgBox;
            msgBox.setText(tr("readggbs"));
            msgBox.setInformativeText(tr("Wrong number of coordinates"));
            msgBox.setIcon(QMessageBox::Critical);
            msgBox.exec();
            return false;
        }
    }
    return true;
}
//  End of function readggbs
//  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//  Function readxyz: reads geometry and charges of centers from a xyz file
//
bool glWidget::readxyz(QString filename){
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox msgBox;
        msgBox.setText(tr("readxyz"));
        msgBox.setInformativeText(tr("File %1 cannot be read").arg(filename)+QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        return false;
    }
    QTextStream in(&file);
    QString line = in.readLine();
    znuc.clear();
    xyz.clear();
    int ncen;
    if ((ncen = line.toInt()) <= 0){
        QMessageBox msgBox;
        msgBox.setText(tr("readxyz"));
        msgBox.setInformativeText(tr("First line of file %1 does not contain the number of centers").arg(filename)
                                  +QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        return false;
    }
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (fields.count() == 1 && znuc.length() == ncen){
            ncen += fields.takeFirst().toInt();
        }
        if (fields.count() != 4){
            continue;
        }
        QString ba = fields.takeFirst().toLatin1();
        znuc << elem->getZsymbol(ba);
        double x = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
        double y = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
        double z = fields.takeFirst().toDouble()*ANGSTROM_TO_BOHR;
        xyz << QVector3D(x, y, z);

    }
    QFileInfo fileInfo(filename);
    fullmolname = QString(fileInfo.absoluteFilePath());
    molname = QString(fileInfo.baseName());
    molpath = QString(fileInfo.absolutePath());
    if (znuc.length() != ncen){
        QMessageBox msgBox;
        msgBox.setText(tr("readxyz"));
        msgBox.setInformativeText(tr("Wrong number of coordinates"));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        return false;
    }
    return true;
}
//  End of function readxyz
//  ---------------------------------------------------------------------------------------------------------------------------

//  Function retrievemolecule: retrieves a molecule
//
void glWidget::retrievemolecule(){
    molecules->append(new molecule);
    molecules->last()->loadcharges(znuc);
    molecules->last()->loadxyz(xyz);
    molecules->last()->setpath(molpath);
    molecules->last()->setrotation(rotation.normalized());
    molecules->last()->setvisible(true);
    molecules->last()->initatomactive(false);
    connections << connect(molecules->last(), SIGNAL(updateGL()), this, SLOT(updateGL()));
    connections << connect(molecules->last(), SIGNAL(updatedisplay()), this, SLOT(updatedisplay()));
    molecules->last()->setname(molname);
    molecules->last()->setfullname(fullmolname);
    molecules->last()->makeStructureBondsandSticks();
    if (window && window->isVisible()){
        molecules->last()->emitupdatedisplay();
    }
    if (QDLwindow){
        updateWindowDialog(gdock);
    }
}

//-------------------------------------------------------------------------------------------------
//
//          moleculeDialog SLOTS
//
//-------------------------------------------------------------------------------------------------

moleculeDialog::moleculeDialog(QWidget *parent) : QDialog()
{
    axesdialog = new axes(this);
    ballsandcylsdialog = new ballscyls(this);
    lightsdialog = new lights(this);
    mespimizerdialog = new mespimizer(this);
    rotationsdialog = new rotations(this);
    scrshotdialog = new screenshot(this);
    translationsdialog = new translations(this);
    viewportdialog = new viewport(this);
}

moleculeDialog::~moleculeDialog(){

}

void moleculeDialog::closeEvent(QCloseEvent *event){
    event->ignore();
    emit closedialog();
}

void moleculeDialog::renewwidth(QSize a){
    resize(a.width()-10,height());
}

//-------------------------------------------------------------------------------------------------
//
//    MainWindow3DViewer SLOTS
//
//-------------------------------------------------------------------------------------------------

MainWindow3DViewer::MainWindow3DViewer(QWidget *parent) : QMainWindow()
{
    statusBar();
}

MainWindow3DViewer::~MainWindow3DViewer(){

}

void MainWindow3DViewer::reject(){
    emit hide_viewer();
}

void MainWindow3DViewer::closeEvent(QCloseEvent *event){
    emit hide_viewer();
    event->accept();
}

//-------------------------------------------------------------------------------------------------
//
//    MyDockWidget SLOTS
//
//-------------------------------------------------------------------------------------------------

MyDockWidget::MyDockWidget(QWidget *parent) : QDockWidget()
{

}

MyDockWidget::~MyDockWidget(){
}


void MyDockWidget::resizeEvent(QResizeEvent *event){
    event->ignore();
    emit resizedialog(this->size());
    event->accept();
}


