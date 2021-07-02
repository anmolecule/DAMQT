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
//	Implementation of class molecule
//
//	File:   molecule.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#include <QtCore/qmath.h>
#include <QDialog>
#include <QCheckBox>
#include <QColorDialog>
#include <QFileDialog>
#include <QFontDialog>
#include <QGroupBox>
#include <QMessageBox>
#include <QPushButton>
#include <QPoint>
#include <QSignalMapper>
#include <QVBoxLayout>

#include "molecule.h"
#include <cmath>

molecule::molecule(QWidget *parent) : QWidget(parent)
{
    allaxesindices.clear();
    allaxesindicesoffset.clear();
    allaxesvertices.clear();
    allindices.clear();
    allindicesoffset.clear();
    allvertices.clear();
    atomactive.clear();
    coneindices.clear();
    conenormals.clear();
    conevertices.clear();
    connections.clear();
    cylinderindices.clear();
    cylindervertices.clear();
    positionaxeslabels.clear();
    xyz.clear();
    znuc.clear();
    setdrawatomcoords(false);
    setdrawatomindices(false);
    setdrawatomsymbols(false);
    setvisible(true);
    coordprecision = 2;
    darkenshift = QVector4D(0.3f,0.3f,0.3f,0.f);
    labelsvshift = 0;
    translation = QVector3D(0,0,0);
    ballradius = 0.2;
    cylradius = 0.05;
    disthressq = pow((INIT_BOND_THRESHOLD * ANGSTROM_TO_BOHR),2);
    deltaAngles = 4.f;
    elem = new Elements();
    font = QFont("Helvetica", 20, QFont::Bold);
    fontaxeslabels = QFont("Noto Sans", 20, QFont::Bold);
    fontcolor = QColor(255, 172, 0, 255);
    dltinterval = (MAX_INTERVAL-MIN_INTERVAL)/float(INTERVAL_SCALE);
    interval = MAX_INTERVAL - dltinterval * INTERVAL_INI;
    angstrom = false;
    angstromcoor = false;
    angstromcp = false;
    axes_visible = false;
    axeslabels_visible = false;
    for (int i = 0 ; i < MAX_CPS; i++){
        cpschecked[i] = true;
    }
    hideatoms = false;
    hidebonds = false;
    hidehydrogens = false;
    iscluster = false;
    startanimation = false;
    visible = true;
    onlyatomactive = false;
    rotatex = false;
    rotatey = false;
    rotatez = false;
    scaleradii = true;
    path = QString(".");
    ProjectFolder = "";
    grids = new QList<grid*>();
    parentposition = QPoint(200,200);
    scrollAreaposition = QPoint(200,200);
    surfaces = new QList<surface*>();
    BTNaddaxes = nullpointer;
    BTNaddcriticalpoints = nullpointer;
    BTNaddfieldlines = nullpointer;
    BTNaddforces = nullpointer;
    BTNaddgrid = nullpointer;
    BTNaddsurface = nullpointer;
    BTNanimation = nullpointer;
    BTNcpcolorfont = nullpointer;
    BTNcpeigcolor[0] = nullpointer;
    BTNcpeigcolor[1] = nullpointer;
    BTNcpeigcolor[2] = nullpointer;
    BTNcps = nullpointer;
    BTNcpselectall = nullpointer;
    BTNcpselectnone = nullpointer;
    BTNcplblfont = nullpointer;
    BTNfieldlines = nullpointer;
    BTNfont = nullpointer;
    BTNfontcolor = nullpointer;
    BTNfontaxeslabels = nullpointer;
    BTNhide = nullpointer;
    BTNrotation = nullpointer;
    BTNselectall = nullpointer;
    BTNselectnone = nullpointer;   
    BTNskeleton = nullpointer;
    BTNsymbols = nullpointer;
    BTNtranslation = nullpointer;
    BTNXaxiscolor = nullpointer;
    BTNYaxiscolor = nullpointer;
    BTNZaxiscolor = nullpointer;
    CHKactiveonly = nullpointer;
    CHKcpactiveonly = nullpointer;
    CHKcpcoords = nullpointer;
    CHKcpeigvec = nullpointer;
    CHKcpindices = nullpointer;
    for (int i = 0 ; i < MAX_CPS; i++){
        CHKcps[i] = nullpointer;
    }
    CHKcpsymbols = nullpointer;
    CHKcpvalues = nullpointer;
    CHKflines = nullpointer;
    CHKflinesarrows = nullpointer;
    for (int i = 0 ; i < MAX_FORCES; i++){
        CHKforces[i] = nullpointer;
    }
    CHKhideatoms = nullpointer;
    CHKhidebonds = nullpointer;
    CHKhidehydrogens = nullpointer;
    CHKrotatex = nullpointer;
    CHKrotatey = nullpointer;
    CHKrotatez = nullpointer;
    CHKshowaxes = nullpointer;
    CHKshowaxeslabels = nullpointer;
    CHKshowcoords = nullpointer;
    CHKshowindices = nullpointer;
    CHKshowsymbols = nullpointer;
    cps = nullpointer;
    flines = nullpointer;
    hfforces = nullpointer;
    FRMaxes = nullpointer;
    FRMcps = nullpointer;
    FRMcpeigvec = nullpointer;
    FRMcpunits = nullpointer;
    FRMcriticalpoints = nullpointer;
    FRMfield = nullpointer;
    FRMrotation = nullpointer;
    FRMsymbols = nullpointer;
    FRMtranslation = nullpointer;
    FRMtranslationunits = nullpointer;
    layoutgrids = nullpointer;
    layoutsurfs = nullpointer;
    LBLcoordprecision = nullpointer;
    LBLcpcoordprecision = nullpointer;
    LBLcpprecision = nullpointer;
    LBLcpselect = nullpointer;
    LBLloadinggrid = nullpointer;
    LBLselect = nullpointer;
    QDLeditMolecule = nullpointer;
    RBTangstrom = nullpointer;
    RBTangstromcoor = nullpointer;
    RBTangstromcp = nullpointer;
    RBTbohr = nullpointer;
    RBTbohrcoor = nullpointer;
    RBTbohrcp = nullpointer;
    scrollArea = nullpointer;
    SPBaxesarrowsize = nullpointer;
    SPBaxesarrowwidth = nullpointer;
    SPBaxeslength = nullpointer;
    SPBaxesthickness = nullpointer;
    SPBcoordprecision = nullpointer;
    SPBcpballradius = nullpointer;
    SPBcpcoordprecision = nullpointer;
    SPBcpeigarrowsize = nullpointer;
    SPBcpeigarrowwidth = nullpointer;
    SPBcpeigthickness = nullpointer;
    SPBcpeiglength = nullpointer;
    SPBcpprecision = nullpointer;
    SPBcpvshift = nullpointer;
    SPBflinesarrowssep = nullpointer;
    SPBflinesarrowssize = nullpointer;
    SPBflinesarrowswidth = nullpointer;
    SPBflineslinewidth = nullpointer;
    SPBlabelsvshift = nullpointer;
    SPBrot_angle = nullpointer;
    SPBrot_x = nullpointer;
    SPBrot_y = nullpointer;
    SPBrot_z = nullpointer;
    SPBtras_x = nullpointer;
    SPBtras_y = nullpointer;
    SPBtras_z = nullpointer;
    TXTcps = nullpointer;
    TXTfieldlines = nullpointer;
    TXTforces = nullpointer;
    timer = new QTimer();
    setactive(true);
    axesarrowssize = 9;
    axesarrowswidth = 4;
    axeslength = 10;
    axesthickness = 2;
    stepwheel = 0.1;
    name = "";
    Xaxis_color = QColor(0,255,0);
    Yaxis_color = QColor(0,0,255);
    Zaxis_color = QColor(255,0,0);
    makeAxesCylinder(15,15);    // computes vertices of a cylinder and its indices
    makeAxesCone(15,15,2.);    // computes vertices of a cone and its indices
}

molecule::~molecule(){
    if (scrollArea){
        delete scrollArea;
        scrollArea = nullpointer;
        QDLeditMolecule = nullpointer;
    }
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    if (grids){
        for (int i = grids->count()-1 ; i >= 0  ; i--){
            delete grids->at(i);
        }
        grids->clear();
        delete grids;
        grids = nullpointer;
    }
    if (surfaces){
        for (int i = surfaces->count()-1 ; i >= 0  ; i--){
            delete surfaces->at(i);
        }
        surfaces->clear();
        delete surfaces;
        surfaces = nullpointer;
    }
    if (hfforces){
        delete hfforces;
        hfforces = nullpointer;
    }
    if (flines){
        delete flines;
        flines = nullpointer;
    }
    if (cps){
        delete cps;
        cps = nullpointer;
    }
    delete timer;
    timer = nullpointer;

}

//  ------------------------------------------------------------------------------------------------------------------
//
//          Molecule editor
//
//  ------------------------------------------------------------------------------------------------------------------

void molecule::editmolecule(){
    if (QDLeditMolecule && QDLeditMolecule->isVisible()){
        if (scrollArea)
            scrollArea->raise();
        return;
    }
    if (QDLeditMolecule){
        QDLeditMolecule_delete();
    }
    if (scrollArea){
        scrollAreaposition = scrollArea->pos();
    }
    createeditMoleculeDialog();
    emitupdateRightMenu();
}

//      Molecule editor dialog
void molecule::createeditMoleculeDialog(){
    closeisosurfeditors();
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
    QDLeditMolecule = new editMoleculeDialog(this);
    QDLeditMolecule->setMinimumWidth(400);

    connections << connect(QDLeditMolecule,SIGNAL(closed()),this,SLOT(QDLeditMolecule_close()));
    connections << connect(QDLeditMolecule,SIGNAL(closed()),this,SLOT(emitupdateRightMenu()));

//        Molecular skeleton
    create_molecular_skeleton_widgets_and_layouts();

//        Symbols and indices
    create_symbols_indices_widgets_and_layouts();

//        Rotation
    create_rotation_widgets_and_layouts();

//        Translation
    create_translation_widgets_and_layouts();

//        Axes
    create_axes_widgets_and_layouts();

//        Hellmann-Feynman forces
    create_forces_widgets_and_layouts();

//        Field lines
    create_field_lines_widgets_and_layouts();

//        Critical points
    create_critical_points_widgets_and_layouts();

//        Surfaces
    create_surfaces_widgets_and_layouts();

//        Grids
    create_grids_widgets_and_layouts();

//        Hide/show
    if (!BTNhide)
        BTNhide = new QPushButton(QDLeditMolecule);
    BTNhide->setText(tr("Hide this menu"));
    connections << connect(BTNhide, SIGNAL(clicked()), this, SLOT(QDLeditMolecule_close()));

//        Layouts
    QVBoxLayout *layout = new QVBoxLayout(QDLeditMolecule);
    layout->addStretch();
    layout->addWidget(BTNhide);
    layout->addWidget(BTNskeleton);
    layout->addWidget(FRMskeleton);
    layout->addWidget(BTNsymbols);
    layout->addWidget(FRMsymbols);
    layout->addWidget(BTNrotation);
    layout->addWidget(FRMrotation);
    layout->addWidget(BTNtranslation);
    layout->addWidget(FRMtranslation);
    layout->addWidget(BTNaddaxes);
    layout->addWidget(FRMaxes);
    layout->addWidget(BTNaddforces);
    layout->addWidget(FRMforces);
    layout->addWidget(BTNaddfieldlines);
    layout->addWidget(FRMfield);
    layout->addWidget(BTNaddcriticalpoints);
    layout->addWidget(FRMcriticalpoints);
    layout->addWidget(BTNaddsurface);
    layout->addLayout(layoutsurfs);
    layout->addWidget(BTNaddgrid);
    layout->addLayout(layoutgrids);
//    layout->addWidget(LBLloadinggrid);
    layout->addStretch();

    if (scrollArea){
        scrollAreaposition = scrollArea->pos();
        delete scrollArea;
        scrollArea = nullpointer;
    }
    scrollArea = new myScrollArea();
    scrollArea->setWidget(QDLeditMolecule);
    scrollArea->setWindowTitle(name);
    scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    connections << connect(scrollArea,SIGNAL(closed()),this,SLOT(QDLeditMolecule_close()));
    connections << connect(scrollArea,SIGNAL(closed()),this,SLOT(emitupdateRightMenu()));
    if (!scrollArea->isVisible()){
        scrollArea->move(scrollAreaposition);
        scrollArea->updatesize(QDLeditMolecule->size());
        scrollArea->show();
    }
}
// End of molecule editor dialog


void molecule::create_molecular_skeleton_widgets_and_layouts(){

//        Show/hide molecular skeleton widgets
    BTNskeleton = new QPushButton(tr("Molecular skeleton"));
    connections << connect(BTNskeleton, SIGNAL(clicked()), this, SLOT(BTNskeleton_clicked()));

    FRMskeleton = new QGroupBox(tr("Atoms and bonds"));
    FRMskeleton->setVisible(false);

    CHKhideatoms = new QCheckBox(FRMskeleton);
    CHKhideatoms->setText(tr("Hide atoms"));
    CHKhideatoms->setChecked(false);
    connections << connect(CHKhideatoms,SIGNAL(stateChanged(int)), this, SLOT(CHKhideatoms_changed()));

    CHKhidebonds = new QCheckBox(FRMskeleton);
    CHKhidebonds->setText(tr("Hide bonds"));
    CHKhidebonds->setChecked(false);
    connections << connect(CHKhidebonds,SIGNAL(stateChanged(int)), this, SLOT(CHKhidebonds_changed()));

    CHKhidehydrogens = new QCheckBox(tr("Hide hydrogens"));
    CHKhidehydrogens->setChecked(false);
    connections << connect(CHKhidehydrogens, SIGNAL(stateChanged(int)), this, SLOT(CHKhidehydrogens_changed()));

//          Layouts

    QVBoxLayout *layout = new QVBoxLayout(FRMskeleton);
    layout->addWidget(CHKhideatoms);
    layout->addWidget(CHKhidebonds);
    layout->addWidget(CHKhidehydrogens);
}

void molecule::create_symbols_indices_widgets_and_layouts(){

//        Symbols and indices
    BTNsymbols = new QPushButton(tr("Labels"));
    connections << connect(BTNsymbols, SIGNAL(clicked()), this, SLOT(BTNsymbols_clicked()));

    FRMsymbols = new QGroupBox(tr("Symbols and indices"));
    FRMsymbols->setVisible(false);

    CHKshowsymbols = new QCheckBox(FRMsymbols);
    CHKshowsymbols->setText(tr("Atom symbols"));
    CHKshowsymbols->setChecked(drawatomsymbols);
    connections << connect(CHKshowsymbols,SIGNAL(stateChanged(int)), this, SLOT(CHKshowsymbols_changed(int)));

    CHKshowindices = new QCheckBox(FRMsymbols);
    CHKshowindices->setText(tr("Atom indices"));
    CHKshowindices->setChecked(drawatomindices);
    connections << connect(CHKshowindices,SIGNAL(stateChanged(int)), this, SLOT(CHKshowindices_changed(int)));


//         Atom coordinates

    CHKshowcoords = new QCheckBox(tr("Atom coordinates"));
    CHKshowcoords->setChecked(drawatomcoords);
    connections << connect(CHKshowcoords, SIGNAL(stateChanged(int)), this, SLOT(CHKshowcoords_changed(int)));

    LBLcoordprecision = new QLabel(tr("Precision"));
    LBLcoordprecision->setHidden(true);
    SPBcoordprecision = new QSpinBox();
    SPBcoordprecision->setMinimum(0);
    SPBcoordprecision->setMaximum(7);
    SPBcoordprecision->setMaximumWidth(50);
    SPBcoordprecision->setValue(coordprecision);
    SPBcoordprecision->setHidden(true);
    connections << connect(SPBcoordprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBcoordprecision_changed(int)));


    FRMcoorunits = new QGroupBox();
    RBTangstromcoor = new QRadioButton(tr("angstrom"),FRMcoorunits);
    RBTbohrcoor = new QRadioButton(tr("bohr"),FRMcoorunits);
    RBTbohrcoor->setChecked(!angstromcp);
    RBTangstromcoor->setVisible(CHKshowcoords->isChecked());
    RBTangstromcoor->setChecked(angstromcoor);
    RBTbohrcoor->setVisible(CHKshowcoords->isChecked());
    connections << connect(RBTbohrcoor, SIGNAL(toggled (bool)),this,SLOT(RBTbohrcoor_changed()));


    CHKactiveonly = new QCheckBox(tr("Only selected atoms"));
    CHKactiveonly->setChecked(onlyatomactive);
    CHKactiveonly->setVisible(drawatomsymbols || drawatomindices);
    connections << connect(CHKactiveonly, SIGNAL(stateChanged(int)), this, SLOT(CHKactiveonly_changed(int)));

    LBLselect = new QLabel(tr("Select")+":");
    LBLselect->setHidden(true);
    QColor *btnbkgcolor = new QColor(29,124,31);
    BTNselectall = new ColorButton(FRMsymbols);
    BTNselectall->setText(tr("All"));
    BTNselectall->setColor(btnbkgcolor);
    BTNselectall->setHidden(true);
    connections << connect(BTNselectall, SIGNAL(clicked()), this, SLOT(BTNselectall_clicked()));
    btnbkgcolor = new QColor(255,0,0);
    BTNselectnone = new ColorButton(FRMsymbols);
    BTNselectnone->setText(tr("None"));
    BTNselectnone->setColor(btnbkgcolor);
    BTNselectnone->setHidden(true);
    connections << connect(BTNselectnone, SIGNAL(clicked()), this, SLOT(BTNselectnone_clicked()));

//        Fonts
    BTNfont = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNfont, SIGNAL(clicked()), this, SLOT(BTNfont_clicked()));
    BTNfontcolor = new ColorButton(FRMsymbols);
    BTNfontcolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNfontcolor->setText(tr("Color"));
    BTNfontcolor->setColor(&fontcolor);
    BTNfontcolor->setEnabled(true);
    connections << connect(BTNfontcolor, SIGNAL(clicked()), this, SLOT(BTNfontcolor_clicked()));

    QLabel *LBLlabelsvshift = new QLabel(tr("Vertical shift"));
    SPBlabelsvshift = new QSpinBox();
    SPBlabelsvshift->setRange(-50,50);
    SPBlabelsvshift->setMaximumWidth(70);
    SPBlabelsvshift->setValue(labelsvshift);
    connections << connect(SPBlabelsvshift,SIGNAL(valueChanged(int)),this,SLOT(SPBlabelsvshift_changed(int)));


    QGridLayout *Layout1=new QGridLayout();
    Layout1->addWidget(CHKshowsymbols,0,0);
    Layout1->addWidget(CHKshowindices,1,0);
    Layout1->addWidget(CHKshowcoords,2,0,Qt::AlignLeft);
    Layout1->addWidget(LBLcoordprecision,2,1,Qt::AlignLeft);
    Layout1->addWidget(SPBcoordprecision,2,2,Qt::AlignLeft);
    Layout1->addWidget(RBTbohrcoor,3,1,Qt::AlignLeft);
    Layout1->addWidget(RBTangstromcoor,3,2,Qt::AlignLeft);


    QHBoxLayout *Layout2=new QHBoxLayout();
    Layout2->addWidget(CHKactiveonly);
    Layout2->setAlignment(Qt::AlignLeft);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addWidget(LBLselect);
    Layout3->addWidget(BTNselectall);
    Layout3->addWidget(BTNselectnone);

    QHBoxLayout *Layout4=new QHBoxLayout();
    Layout4->addWidget(LBLlabelsvshift);
    Layout4->addWidget(SPBlabelsvshift);

    QHBoxLayout *Layout5=new QHBoxLayout();
    Layout5->addStretch();
    Layout5->addWidget(BTNfont);
    Layout5->addWidget(BTNfontcolor);
    Layout5->addStretch();
    Layout5->setAlignment(Qt::AlignCenter);

    QVBoxLayout *layout = new QVBoxLayout(FRMsymbols);
    layout->addLayout(Layout1);
    layout->addLayout(Layout2);
    layout->addLayout(Layout3);
    layout->addLayout(Layout4);
    layout->addLayout(Layout5);

//    QHBoxLayout *layout1=new QHBoxLayout();
//    layout1->addWidget(LBLlabelsvshift);
//    layout1->addWidget(SPBlabelsvshift);
//    layout1->setAlignment(Qt::AlignLeft);

//    QHBoxLayout *layout2=new QHBoxLayout();
//    layout2->addStretch();
//    layout2->addWidget(LBLcoordprecision);
//    layout2->addWidget(SPBcoordprecision);
//    layout2->addStretch();

//    QHBoxLayout *layout3 = new QHBoxLayout();
//    layout3->addStretch();
//    layout3->addWidget(LBLselect);
//    layout3->addWidget(BTNselectall);
//    layout3->addWidget(BTNselectnone);
//    layout3->setAlignment(Qt::AlignLeft);
//    layout3->addStretch();

//    QHBoxLayout *layout4 = new QHBoxLayout();
//    layout4->addStretch();
//    layout4->addWidget(BTNfontcolor);
//    layout4->addWidget(BTNfont);
//    layout4->addStretch();

//    QVBoxLayout *layout = new QVBoxLayout(FRMsymbols);
//    layout->addWidget(CHKshowsymbols);
//    layout->addWidget(CHKshowindices);
//    layout->addWidget(CHKshowcoords);
//    layout->addLayout(layout2);
//    layout->addWidget(CHKactiveonly);
//    layout->addLayout(layout3);
//    layout->addLayout(layout1);
//    layout->addLayout(layout4);
}

void molecule::create_rotation_widgets_and_layouts(){

//        Rotations
    BTNrotation = new QPushButton(tr("Rotations"));
    connections << connect(BTNrotation, SIGNAL(clicked()), this, SLOT(BTNrotation_clicked()));

    FRMrotation = new QGroupBox(tr("Rotation"));
    FRMrotation->setVisible(false);
    QGroupBox *FRMaxis = new QGroupBox(tr("Axis"));
    QGroupBox *FRMangle = new QGroupBox(tr("Angle"));

    QLabel *LBLrot_x = new QLabel("x:");
    SPBrot_x=new DoubleSpinBox();
    SPBrot_x->setDecimals(3);
    SPBrot_x->setSingleStep(0.01);
    SPBrot_x->setRange(-1.,1.);
    SPBrot_x->setValue(rotation.x());
    SPBrot_x->setEnabled(true);
    connections << connect(SPBrot_x,SIGNAL(valueChanged(double)),this,SLOT(rotation_changed()));

    QLabel *LBLrot_y = new QLabel("y:");
    SPBrot_y=new DoubleSpinBox();
    SPBrot_y->setDecimals(3);
    SPBrot_y->setSingleStep(0.01);
    SPBrot_y->setRange(-1.,1.);
    SPBrot_y->setValue(rotation.y());
    SPBrot_y->setEnabled(true);
    connections << connect(SPBrot_y,SIGNAL(valueChanged(double)),this,SLOT(rotation_changed()));

    QLabel *LBLrot_z = new QLabel("z:");
    SPBrot_z=new DoubleSpinBox();
    SPBrot_z->setDecimals(3);
    SPBrot_z->setSingleStep(0.01);
    SPBrot_z->setRange(-1.,1.);
    SPBrot_z->setValue(rotation.z());
    SPBrot_z->setEnabled(true);
    connections << connect(SPBrot_z,SIGNAL(valueChanged(double)),this,SLOT(rotation_changed()));

    QLabel *LBLrot_w = new QLabel("w:");
    SPBrot_angle=new DoubleSpinBox();
    SPBrot_angle->setDecimals(1);
    SPBrot_angle->setSingleStep(5);
    SPBrot_angle->setRange(0,360);
    SPBrot_angle->setValue(360. * qAcos(rotation.scalar()) / pi);
    SPBrot_angle->setEnabled(true);
    connections << connect(SPBrot_angle,SIGNAL(valueChanged(double)),this,SLOT(rotation_changed()));

    QPushButton *BTNapplyrot = new QPushButton();
    BTNapplyrot->setText(tr("Apply"));
    connections << connections << connect(BTNapplyrot,SIGNAL(clicked()),this,SLOT(rotation_changed()));

    QPushButton *BTNresetrot = new QPushButton();
    BTNresetrot->setText(tr("Reset"));
    connections << connections << connect(BTNresetrot,SIGNAL(clicked()),this,SLOT(reset_rotation()));

    QGroupBox *FRManimate = new QGroupBox(tr("Animate rotation"));

    CHKrotatex = new QCheckBox(FRManimate);
    CHKrotatex->setText(tr("X axis"));
    CHKrotatex->setChecked(rotatex);
    connections << connect(CHKrotatex, SIGNAL(stateChanged(int)), this, SLOT(CHKrotate_changed()));

    CHKrotatey = new QCheckBox(FRManimate);
    CHKrotatey->setText(tr("Y axis"));
    CHKrotatey->setChecked(rotatey);
    connections << connect(CHKrotatey, SIGNAL(stateChanged(int)), this, SLOT(CHKrotate_changed()));

    CHKrotatez = new QCheckBox(FRManimate);
    CHKrotatez->setText(tr("Z axis"));
    CHKrotatez->setChecked(rotatez);
    connections << connect(CHKrotatez, SIGNAL(stateChanged(int)), this, SLOT(CHKrotate_changed()));

    BTNanimation = new QPushButton(QIcon(":/images/empezar.png"),tr("Start"));
    if (startanimation)
        BTNanimation->setText(tr("Stop"));
    connections << connect(BTNanimation,SIGNAL(clicked()), this, SLOT(BTNanimation_clicked()));

    connections << connect(timer, SIGNAL(timeout()), this, SLOT(animaterotation()));

    SLDspeed = new QSlider(Qt::Horizontal);
    SLDspeed->setRange(0,INTERVAL_SCALE);
    SLDspeed->setSingleStep(1);
    SLDspeed->setPageStep(10);
    SLDspeed->setTickPosition(QSlider::TicksBelow);
    SLDspeed->setValue(INTERVAL_INI);
    connections << connect(SLDspeed,SIGNAL(valueChanged(int)), this, SLOT(resetinterval()));

    QLabel *LBLslow = new QLabel(tr("Slow"));
    QLabel *LBLfast = new QLabel(tr("Fast"));

    QGridLayout *layout4 = new QGridLayout(FRMaxis);
    layout4->addWidget(LBLrot_x,0,0);
    layout4->addWidget(SPBrot_x,0,1);
    layout4->addWidget(LBLrot_y,0,2);
    layout4->addWidget(SPBrot_y,0,3);
    layout4->addWidget(LBLrot_z,0,4);
    layout4->addWidget(SPBrot_z,0,5);

    QHBoxLayout *layout5 = new QHBoxLayout(FRMangle);
    layout5->addWidget(LBLrot_w);
    layout5->addWidget(SPBrot_angle);

    QHBoxLayout *layout6 = new QHBoxLayout();
    layout6->addWidget(FRMaxis);
    layout6->addWidget(FRMangle);

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addStretch();
    layout7->addWidget(BTNresetrot);
    layout7->addWidget(BTNapplyrot);
    layout7->addStretch();

    QGridLayout *layout7b = new QGridLayout();
    layout7b->addWidget(CHKrotatex,0,0,Qt::AlignCenter);
    layout7b->addWidget(CHKrotatey,0,1,Qt::AlignCenter);
    layout7b->addWidget(CHKrotatez,0,2,Qt::AlignCenter);
    layout7b->addWidget(BTNanimation,1,1,Qt::AlignCenter);

    QHBoxLayout *layout7c = new QHBoxLayout();
    layout7c->addWidget(LBLslow);
    layout7c->addWidget(SLDspeed);
    layout7c->addWidget(LBLfast);

    QVBoxLayout *layout7d = new QVBoxLayout(FRManimate);
    layout7d->addLayout(layout7b);
    layout7d->addLayout(layout7c);

    QVBoxLayout *layout8 = new QVBoxLayout(FRMrotation);
    layout8->addLayout(layout6);
    layout8->addLayout(layout7);
    layout8->addWidget(FRManimate);
}

void molecule::create_translation_widgets_and_layouts(){

//        Translations
    BTNtranslation = new QPushButton(tr("Translations"));
    connections << connect(BTNtranslation, SIGNAL(clicked()), this, SLOT(BTNtranslation_clicked()));

    FRMtranslation = new QGroupBox(tr("Translation"));
    FRMtranslation->setVisible(false);

    FRMtranslationunits = new QGroupBox(tr("Units"));

    RBTangstrom = new QRadioButton(tr("angstrom"),FRMtranslationunits);
    RBTbohr = new QRadioButton(tr("bohr"),FRMtranslationunits);
    RBTbohr->setChecked(true);
    RBTangstrom->setChecked(false);
    angstrom = false;
    connections << connect(RBTbohr, SIGNAL(toggled (bool)),this,SLOT(RBTbohr_changed()));

    QLabel *LBLtras_x = new QLabel("x:");
    SPBtras_x=new DoubleSpinBox();
    SPBtras_x->setDecimals(2);
    SPBtras_x->setSingleStep(0.1);
    SPBtras_x->setRange(-1000,1000);
    SPBtras_x->setValue(translation.x());
    SPBtras_x->setEnabled(true);
    connections << connect(SPBtras_x,SIGNAL(valueChanged(double)),this,SLOT(translation_changed()));

    QLabel *LBLtras_y = new QLabel("y:");
    SPBtras_y=new DoubleSpinBox();
    SPBtras_y->setDecimals(2);
    SPBtras_y->setSingleStep(0.1);
    SPBtras_y->setRange(-1000,1000);
    SPBtras_y->setValue(translation.y());
    SPBtras_y->setEnabled(true);
    connections << connect(SPBtras_y,SIGNAL(valueChanged(double)),this,SLOT(translation_changed()));

    QLabel *LBLtras_z = new QLabel("z:");
    SPBtras_z=new DoubleSpinBox();
    SPBtras_z->setDecimals(2);
    SPBtras_z->setSingleStep(0.1);
    SPBtras_z->setRange(-1000,1000);
    SPBtras_z->setValue(translation.z());
    SPBtras_z->setEnabled(true);
    connections << connect(SPBtras_z,SIGNAL(valueChanged(double)),this,SLOT(translation_changed()));

    QPushButton *BTNapplytrans = new QPushButton();
    BTNapplytrans->setText(tr("Apply"));
    connections << connect(BTNapplytrans,SIGNAL(clicked()),this,SLOT(translation_changed()));

    QPushButton *BTNresetrans = new QPushButton();
    BTNresetrans->setText(tr("Reset"));
    connections << connect(BTNresetrans,SIGNAL(clicked()),this,SLOT(reset_translation()));

    QLabel *LBLstepwheel = new QLabel(tr("Stride for zooming with mouse wheel: "));
    SPBstepwheel=new DoubleSpinBox();
    SPBstepwheel->setDecimals(2);
    SPBstepwheel->setSingleStep(0.1);
    SPBstepwheel->setRange(0,100);
    SPBstepwheel->setValue(stepwheel);
    SPBstepwheel->setEnabled(true);
    connections << connect(SPBstepwheel,SIGNAL(valueChanged(double)),this,SLOT(SPBstepwheel_changed()));
    connections << connect(SPBstepwheel,SIGNAL(textChanged(QString)),this,SLOT(SPBstepwheel_changed()));
    connections << connect(BTNapplytrans,SIGNAL(clicked()),this,SLOT(SPBstepwheel_changed()));

    QHBoxLayout *layout8 = new QHBoxLayout(FRMtranslationunits);
    layout8->addWidget(RBTbohr);
    layout8->addWidget(RBTangstrom);

    QHBoxLayout *layout9 = new QHBoxLayout();
    layout9->addStretch();
    layout9->addWidget(LBLtras_x);
    layout9->addWidget(SPBtras_x);
    layout9->addWidget(LBLtras_y);
    layout9->addWidget(SPBtras_y);
    layout9->addWidget(LBLtras_z);
    layout9->addWidget(SPBtras_z);
    layout9->addStretch();

    QHBoxLayout *layout10 = new QHBoxLayout();
    layout10->addStretch();
    layout10->addWidget(LBLstepwheel);
    layout10->addWidget(SPBstepwheel);
    layout10->addStretch();

    QHBoxLayout *layout11 = new QHBoxLayout();
    layout11->addStretch();
    layout11->addWidget(BTNresetrans);
    layout11->addWidget(BTNapplytrans);
    layout11->addStretch();

    QVBoxLayout *layout12 = new QVBoxLayout(FRMtranslation);
    layout12->addWidget(FRMtranslationunits);
    layout12->addLayout(layout9);
    layout12->addLayout(layout10);
    layout12->addLayout(layout11);
}

void molecule::create_axes_widgets_and_layouts(){

    BTNaddaxes = new QPushButton(tr("Axes"));
    connections << connect(BTNaddaxes, SIGNAL(clicked()), this, SLOT(BTNaddaxes_clicked()));

    FRMaxes = new QGroupBox(tr("Axes"));
    FRMaxes->setVisible(false);

    CHKshowaxes = new QCheckBox(tr("Show axes"));
    CHKshowaxes->setChecked(axes_visible);
    connections << connect(CHKshowaxes, SIGNAL(stateChanged(int)), this, SLOT(CHKshowaxes_changed(int)));

    CHKshowaxeslabels = new QCheckBox(tr("Show axes labels"));
    CHKshowaxeslabels->setChecked(axes_visible && axeslabels_visible);
    CHKshowaxeslabels->setEnabled(axes_visible);
    connections << connect(CHKshowaxeslabels, SIGNAL(stateChanged(int)), this, SLOT(CHKshowaxeslabels_changed(int)));

//        Fonts
    BTNfontaxeslabels = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNfontaxeslabels, SIGNAL(clicked()), this, SLOT(BTNfontaxeslabels_clicked()));

//    Colors

    BTNXaxiscolor = new ColorButton();
    BTNXaxiscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNXaxiscolor->setText(tr("X axis"));
    BTNXaxiscolor->setColor(&Xaxis_color);
    connections << connect(BTNXaxiscolor,SIGNAL(clicked()),this,SLOT(BTNXaxiscolor_clicked()));

    BTNYaxiscolor = new ColorButton();
    BTNYaxiscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNYaxiscolor->setText(tr("Y axis"));
    BTNYaxiscolor->setColor(&Yaxis_color);
    connections << connect(BTNYaxiscolor,SIGNAL(clicked()),this,SLOT(BTNYaxiscolor_clicked()));

    BTNZaxiscolor = new ColorButton();
    BTNZaxiscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNZaxiscolor->setText(tr("Z axis"));
    BTNZaxiscolor->setColor(&Zaxis_color);
    connections << connect(BTNZaxiscolor,SIGNAL(clicked()),this,SLOT(BTNZaxiscolor_clicked()));

//         Thickness

    QLabel *LBLaxesthickness = new QLabel(tr("Axes thickness"));
    SPBaxesthickness = new QSpinBox();
    SPBaxesthickness->setMinimum(1);
    SPBaxesthickness->setMaximum(50);
    SPBaxesthickness->setMaximumWidth(50);
    SPBaxesthickness->setValue(4);
    connections << connect(SPBaxesthickness,SIGNAL(valueChanged(int)),this,SLOT(SPBaxesthickness_changed(int)));

//         Vectors length

    QLabel *LBLaxeslength = new QLabel(tr("Axes length"));
    SPBaxeslength = new QSpinBox();
    SPBaxeslength->setMinimum(1);
    SPBaxeslength->setMaximum(100);
    SPBaxeslength->setMaximumWidth(50);
    SPBaxeslength->setSingleStep(1);
    SPBaxeslength->setValue(10);
    connections << connect(SPBaxeslength,SIGNAL(valueChanged(int)),this,SLOT(SPBaxeslength_changed(int)));

//         Arrowlength

    QLabel *LBLaxesarrowsize = new QLabel(tr("Arrows length"));
    SPBaxesarrowsize = new QSpinBox();
    SPBaxesarrowsize->setMinimum(1);
    SPBaxesarrowsize->setMaximum(50);
    SPBaxesarrowsize->setMaximumWidth(50);
    SPBaxesarrowsize->setSingleStep(1);
    SPBaxesarrowsize->setValue(4);
    connections << connect(SPBaxesarrowsize,SIGNAL(valueChanged(int)),this,SLOT(SPBaxesarrowssize_changed(int)));


//         Arrowwidth

    QLabel *LBLaxesarrowwidth = new QLabel(tr("Arrows width"));
    SPBaxesarrowwidth = new QSpinBox();
    SPBaxesarrowwidth->setMinimum(1);
    SPBaxesarrowwidth->setMaximum(50);
    SPBaxesarrowwidth->setMaximumWidth(50);
    SPBaxesarrowwidth->setSingleStep(1);
    SPBaxesarrowwidth->setValue(10);
    connections << connect(SPBaxesarrowwidth,SIGNAL(valueChanged(int)),this,SLOT(SPBaxesarrowswidth_changed(int)));

    QVBoxLayout *layout1 = new QVBoxLayout();
    layout1->addWidget(CHKshowaxes);
    layout1->addWidget(CHKshowaxeslabels);
    layout1->addStretch();

    QHBoxLayout *layout2 = new QHBoxLayout();
    layout2->addStretch();
    layout2->addWidget(BTNfontaxeslabels);
    layout2->addStretch();

    QHBoxLayout *layout3 = new QHBoxLayout();
    layout3->addStretch();
    layout3->addWidget(BTNXaxiscolor);
    layout3->addWidget(BTNYaxiscolor);
    layout3->addWidget(BTNZaxiscolor);
    layout3->addStretch();

    QGridLayout *layout4 = new QGridLayout();
    layout4->addWidget(LBLaxesthickness,0,0);
    layout4->addWidget(SPBaxesthickness,0,1);
    layout4->addWidget(LBLaxeslength,1,0);
    layout4->addWidget(SPBaxeslength,1,1);
    layout4->addWidget(LBLaxesarrowsize,2,0);
    layout4->addWidget(SPBaxesarrowsize,2,1);
    layout4->addWidget(LBLaxesarrowwidth,3,0);
    layout4->addWidget(SPBaxesarrowwidth,3,1);

    QVBoxLayout *layout = new QVBoxLayout(FRMaxes);
    layout->addLayout(layout1);
    layout->addLayout(layout2);
    layout->addLayout(layout3);
    layout->addLayout(layout4);
    layout->addStretch();
}

void molecule::create_forces_widgets_and_layouts(){
//        Add Hellmann-Feynman forces
    BTNaddforces = new QPushButton(tr("Hellmann-Feynman forces"));
    connections << connect(BTNaddforces, SIGNAL(clicked()), this, SLOT(BTNaddforces_clicked()));
    addforces();

    FRMforces = new QGroupBox(tr("Hellmann-Feynman forces"));
    FRMforces->setVisible(false);

    TXTforces = new QLineEdit();
    TXTforces->setText(hfforces->getforcesfilename());
    connections << connect(TXTforces, SIGNAL(returnPressed()), this, SLOT(TXTforces_changed()));
    BTNforces = new QToolButton();
    BTNforces->setText(tr("..."));
    BTNforces->setToolTip(tr("Open file with Hellmann-Feynman forces ..."));
    connections << connect(BTNforces, SIGNAL(clicked()), this, SLOT(readforcefiles_dialog()));

    QLabel *LBLforcecolors0 = new QLabel(tr("External")+":");
    QLabel *LBLforcecolors1 = new QLabel(tr("Internal")+":");
    QLabel *LBLforcecolors2 = new QLabel(tr("Total")+":");
    QLabel *LBLforcecolors3 = new QLabel(tr("Non-conformational")+":");
    QLabel *LBLforcecolors4 = new QLabel(tr("Conformational")+":");

    for (int i=0 ; i<MAX_FORCES ; i++){
        CHKforces[i] = new QCheckBox();
        CHKforces[i]->setChecked(hfforces->getvisibleforces(i));
        connections << connect(CHKforces[i], SIGNAL(stateChanged(int)), this, SLOT(CHKforces_changed()));
        hfforces->setvisibleforces(i,true);
        BTNforcecolors[i] = new ColorButton();
        BTNforcecolors[i]->setIcon(QIcon(":/images/colores48.png"));
        BTNforcecolors[i]->setText(tr("Color"));
        QColor color = hfforces->getcolor(i);
        BTNforcecolors[i]->setColor(&color);
        connections << connect(BTNforcecolors[i],SIGNAL(clicked()),this,SLOT(BTNforcescolor_changed()));
    }

//         Vectors thickness

    QLabel *LBLforcesthickness = new QLabel(tr("Vectors thickness"));
    SPBforcesthickness = new QSpinBox();
    SPBforcesthickness->setMinimum(1);
    SPBforcesthickness->setMaximum(50);
    SPBforcesthickness->setMinimumWidth(60);
    SPBforcesthickness->setMaximumWidth(60);
    SPBforcesthickness->setValue(hfforces->getforcesthickness());
    connections << connect(SPBforcesthickness,SIGNAL(valueChanged(int)),this,SLOT(SPBforcesthickness_changed(int)));

//         Vectors length

    QLabel *LBLforceslength = new QLabel(tr("Vectors length"));
    SPBforceslength = new QSpinBox();
    SPBforceslength->setMinimum(1);
    SPBforceslength->setMaximum(200);
    SPBforceslength->setMinimumWidth(60);
    SPBforceslength->setMaximumWidth(60);
    SPBforceslength->setSingleStep(1);
    SPBforceslength->setValue(hfforces->getforceslength());
    connections << connect(SPBforceslength,SIGNAL(valueChanged(int)),this,SLOT(SPBforceslength_changed(int)));

//         Arrows length

    QLabel *LBLforcesarrowlength = new QLabel(tr("Arrows length"));
    SPBforcesarrowlength = new QSpinBox();
    SPBforcesarrowlength->setMinimum(1);
    SPBforcesarrowlength->setMaximum(50);
    SPBforcesarrowlength->setMinimumWidth(60);
    SPBforcesarrowlength->setMaximumWidth(60);
    SPBforcesarrowlength->setSingleStep(1);
    SPBforcesarrowlength->setValue(hfforces->getarrowlength());
    connections << connect(SPBforcesarrowlength,SIGNAL(valueChanged(int)),this,SLOT(SPBforcesarrowlength_changed(int)));

//         Arrows width

    QLabel *LBLforcesarrowwidth = new QLabel(tr("Arrows width"));
    SPBforcesarrowwidth = new QSpinBox();
    SPBforcesarrowwidth->setMinimum(1);
    SPBforcesarrowwidth->setMaximum(50);
    SPBforcesarrowwidth->setMinimumWidth(60);
    SPBforcesarrowwidth->setMaximumWidth(60);
    SPBforcesarrowwidth->setSingleStep(1);
    SPBforcesarrowwidth->setValue(hfforces->getarrowwidth());
    connections << connect(SPBforcesarrowwidth,SIGNAL(valueChanged(int)),this,SLOT(SPBforcesarrowwidth_changed(int)));

//    Layouts
    QLabel *LBLforces = new QLabel(tr("File")+":");

    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(LBLforces);
    Layout1->addWidget(TXTforces);
    Layout1->addWidget(BTNforces);

    QGridLayout *Layout2=new QGridLayout();
    Layout2->addWidget(LBLforcecolors0,0,0,Qt::AlignLeft);
    Layout2->addWidget(CHKforces[0],0,1,Qt::AlignCenter);
    Layout2->addWidget(BTNforcecolors[0],0,2,Qt::AlignRight);
    Layout2->addWidget(LBLforcecolors1,1,0,Qt::AlignLeft);
    Layout2->addWidget(CHKforces[1],1,1,Qt::AlignCenter);
    Layout2->addWidget(BTNforcecolors[1],1,2,Qt::AlignRight);
    Layout2->addWidget(LBLforcecolors2,2,0,Qt::AlignLeft);
    Layout2->addWidget(CHKforces[2],2,1,Qt::AlignCenter);
    Layout2->addWidget(BTNforcecolors[2],2,2,Qt::AlignRight);
    Layout2->addWidget(LBLforcecolors3,3,0,Qt::AlignLeft);
    Layout2->addWidget(CHKforces[3],3,1,Qt::AlignCenter);
    Layout2->addWidget(BTNforcecolors[3],3,2,Qt::AlignRight);
    Layout2->addWidget(LBLforcecolors4,4,0,Qt::AlignLeft);
    Layout2->addWidget(CHKforces[4],4,1,Qt::AlignCenter);
    Layout2->addWidget(BTNforcecolors[4],4,2,Qt::AlignRight);

//        Arrows layout

    QGridLayout *Layout15=new QGridLayout();
    Layout15->addWidget(LBLforcesthickness,0,0,Qt::AlignLeft);
    Layout15->addWidget(SPBforcesthickness,0,1,Qt::AlignRight);
    Layout15->addWidget(LBLforceslength,1,0,Qt::AlignLeft);
    Layout15->addWidget(SPBforceslength,1,1,Qt::AlignRight);
    Layout15->addWidget(LBLforcesarrowlength,2,0,Qt::AlignLeft);
    Layout15->addWidget(SPBforcesarrowlength,2,1,Qt::AlignRight);
    Layout15->addWidget(LBLforcesarrowwidth,3,0,Qt::AlignLeft);
    Layout15->addWidget(SPBforcesarrowwidth,3,1,Qt::AlignRight);


//        HF forces layout

    QVBoxLayout *Layout14=new QVBoxLayout(FRMforces);
    Layout14->addLayout(Layout1);
    Layout14->addLayout(Layout2);
    Layout14->addLayout(Layout15);
}

void molecule::create_field_lines_widgets_and_layouts(){
//        Add field lines
    BTNaddfieldlines = new QPushButton(tr("3D lines"));
    connections << connect(BTNaddfieldlines, SIGNAL(clicked()), this, SLOT(BTNaddfieldlines_clicked()));
    addfield();

    FRMfield = new QGroupBox(tr("3D lines"));
    FRMfield->setVisible(false);

    TXTfieldlines = new QLineEdit();
    TXTfieldlines->setText(flines->getfieldfilename());
    connections << connect(TXTfieldlines, SIGNAL(returnPressed()), this, SLOT(TXTfieldlines_changed()));

    QLabel *LBLfieldlines = new QLabel(tr("File")+":");
    BTNfieldlines = new QToolButton();
    BTNfieldlines->setText(tr("..."));
    BTNfieldlines->setToolTip(tr("Open geometry file ..."));
    connections << connect(BTNfieldlines, SIGNAL(clicked()), this, SLOT(readfieldlines_dialog()));

    CHKflines = new QCheckBox(tr("Show lines"));
    CHKflines->setChecked(flines->isvisible());
    connections << connect(CHKflines, SIGNAL(stateChanged(int)), this, SLOT(CHKflines_changed()));

    BTNflinescolor = new ColorButton(QIcon(":/images/colores48.png"),tr("Color"));
    QColor currcolor = flines->getlinescolor();
    BTNflinescolor->setColor(&currcolor);
    connections << connect(BTNflinescolor, SIGNAL(clicked()), this, SLOT(BTNflinescolor_clicked()));

    QLabel *LBLlinewidth = new QLabel(tr("Width"));
    SPBflineslinewidth = new QSpinBox();
    SPBflineslinewidth->setMinimum(1);
    SPBflineslinewidth->setMaximum(20);
    SPBflineslinewidth->setMaximumWidth(50);
    SPBflineslinewidth->setValue(flines->getlineswidth());
    connections << connect(SPBflineslinewidth,SIGNAL(valueChanged(int)),this,SLOT(SPBflineslinewidth_changed(int)));

    FRMflinesarrows = new QGroupBox();
    FRMflinesarrows->setTitle(tr("Arrows"));

    CHKflinesarrows = new QCheckBox(tr("Show arrows"));
    CHKflinesarrows->setChecked(flines->getshowarrows());
    connections << connect(CHKflinesarrows, SIGNAL(stateChanged(int)), this, SLOT(CHKflinesarrows_changed()));

    QLabel *LBLarrowssep = new QLabel(tr("Arrows spacing"));
    SPBflinesarrowssep = new QSpinBox();
    SPBflinesarrowssep->setMinimum(1);
    SPBflinesarrowssep->setMaximum(500);
    SPBflinesarrowssep->setMaximumWidth(50);
    SPBflinesarrowssep->setSingleStep(5);
    SPBflinesarrowssep->setValue(flines->getarrowsseparation());
    connections << connect(SPBflinesarrowssep,SIGNAL(valueChanged(int)),this,SLOT(SPBflinesarrowssep_changed(int)));

    QLabel *LBLarrowssize = new QLabel(tr("Arrows length"));
    SPBflinesarrowssize = new QSpinBox();
    SPBflinesarrowssize->setMinimum(1);
    SPBflinesarrowssize->setMaximum(25);
    SPBflinesarrowssize->setMaximumWidth(50);
    SPBflinesarrowssize->setSingleStep(1);
    SPBflinesarrowssize->setValue(flines->getarrowssize());
    flines->setarrowssize(SPBflinesarrowssize->value());
    connections << connect(SPBflinesarrowssize,SIGNAL(valueChanged(int)),this,SLOT(SPBflinesarrowssize_changed(int)));

    QLabel *LBLarrowswidth = new QLabel(tr("Arrows width"));
    SPBflinesarrowswidth = new QSpinBox();
    SPBflinesarrowswidth->setMinimum(1);
    SPBflinesarrowswidth->setMaximum(20);
    SPBflinesarrowswidth->setMaximumWidth(50);
    SPBflinesarrowswidth->setSingleStep(1);
    SPBflinesarrowswidth->setValue(flines->getarrowswidth());
    flines->setarrowswidth(SPBflinesarrowswidth->value());
    connections << connect(SPBflinesarrowswidth,SIGNAL(valueChanged(int)),this,SLOT(SPBflinesarrowswidth_changed(int)));

    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(LBLfieldlines);
    Layout1->addWidget(TXTfieldlines);
    Layout1->addWidget(BTNfieldlines);

    QHBoxLayout *Layout2=new QHBoxLayout();
    Layout2->addWidget(CHKflines);
    Layout2->setAlignment(Qt::AlignCenter);

    QHBoxLayout *Layout3=new QHBoxLayout();
    Layout3->addWidget(LBLlinewidth,Qt::AlignCenter);
    Layout3->addWidget(SPBflineslinewidth);
    Layout3->addWidget(BTNflinescolor);

    QHBoxLayout *Layout4 = new QHBoxLayout();
    Layout4->addWidget(CHKflinesarrows,Qt::AlignCenter);

    QGridLayout *Layout5 = new QGridLayout();
    Layout5->addWidget(LBLarrowssep,0,0,Qt::AlignLeft);
    Layout5->addWidget(SPBflinesarrowssep,0,1,Qt::AlignRight);
    Layout5->addWidget(LBLarrowssize,1,0,Qt::AlignLeft);
    Layout5->addWidget(SPBflinesarrowssize,1,1,Qt::AlignRight);
    Layout5->addWidget(LBLarrowswidth,2,0,Qt::AlignLeft);
    Layout5->addWidget(SPBflinesarrowswidth,2,1,Qt::AlignRight);

    QVBoxLayout *Layout6 = new  QVBoxLayout(FRMflinesarrows);
    Layout6->addLayout(Layout4);
    Layout6->addLayout(Layout5);

    QVBoxLayout *Layout = new  QVBoxLayout(FRMfield);
    Layout->addLayout(Layout1);
    Layout->addLayout(Layout2);
    Layout->addLayout(Layout3);
    Layout->addWidget(FRMflinesarrows);
}

void molecule::create_critical_points_widgets_and_layouts(){
//        Add critical points
    BTNaddcriticalpoints = new QPushButton(tr("Critical points"));
    connections << connect(BTNaddcriticalpoints, SIGNAL(clicked()), this, SLOT(BTNaddcriticalpoints_clicked()));
    addcriticalpoints();

    FRMcriticalpoints= new QGroupBox(tr("Critical points"));
    FRMcriticalpoints->setVisible(false);
    FRMcriticalpoints->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

    FRMcps= new QGroupBox();
    FRMcps->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

    QLabel *LBLcps = new QLabel(tr("File")+":");
    TXTcps = new QLineEdit();
    TXTcps->setText(cps->getcpsfilename());
    connections << connect(TXTcps, SIGNAL(returnPressed()), this, SLOT(TXTcps_changed()));

    BTNcps = new QToolButton();
    BTNcps->setText(tr("..."));
    BTNcps->setToolTip(tr("Open file with critical points ..."));
    connections << connect(BTNcps, SIGNAL(clicked()), this, SLOT(readcpsfiles_dialog()));

    QLabel *LBLcpscolor0 = new QLabel(tr("(3,+3) CP")+":");
    QLabel *LBLcpscolor1 = new QLabel(tr("(3,+1) CP")+":");
    QLabel *LBLcpscolor2 = new QLabel(tr("(3,-1) CP")+":");
    QLabel *LBLcpscolor3 = new QLabel(tr("(3,-3) CP")+":");
    for (int i=0 ; i<MAX_CPS ; i++){
        CHKcps[i] = new QCheckBox();
        CHKcps[i]->setChecked(cpschecked[i]);
        connections << connect(CHKcps[i], SIGNAL(stateChanged(int)), this, SLOT(CHKcps_changed()));
        BTNcpcolor[i] = new ColorButton();
        BTNcpcolor[i]->setIcon(QIcon(":/images/colores48.png"));
        BTNcpcolor[i]->setText(tr("Color"));
        QColor currcolor = cps->getcolor(i);
        BTNcpcolor[i]->setColor(&currcolor);
        connect(BTNcpcolor[i],SIGNAL(clicked()),this,SLOT(BTNcpcolor_change()));
    }
//         Ball radius

    QLabel *LBLcpsballradius = new QLabel(tr("Ball radius"));
    SPBcpballradius = new QSpinBox();
    SPBcpballradius->setMinimum(1);
    SPBcpballradius->setMaximum(100);
    SPBcpballradius->setMaximumWidth(50);
    SPBcpballradius->setValue(cps->getradius());
    connections << connect(SPBcpballradius,SIGNAL(valueChanged(int)),this,SLOT(SPBcpballradius_changed(int)));

//         CP coordinates

    CHKcpcoords = new QCheckBox(tr("CP coordinates"));
    CHKcpcoords->setChecked(cps->getdrawcpscoords());
    connections << connect(CHKcpcoords, SIGNAL(stateChanged(int)), this, SLOT(CHKcpcoords_changed(int)));

    FRMcpunits = new QGroupBox();
    RBTangstromcp = new QRadioButton(tr("angstrom"),FRMcpunits);
    RBTbohrcp = new QRadioButton(tr("bohr"),FRMcpunits);
    RBTbohrcp->setChecked(!angstromcp);
    RBTangstromcp->setVisible(CHKcpcoords->isChecked());
    RBTangstromcp->setChecked(angstromcp);
    RBTbohrcp->setVisible(CHKcpcoords->isChecked());
    connections << connect(RBTbohrcp, SIGNAL(toggled (bool)),this,SLOT(RBTbohrcp_changed()));

    LBLcpcoordprecision = new QLabel(tr("Precision"));
    LBLcpcoordprecision->setVisible(CHKcpcoords->isChecked());
    SPBcpcoordprecision = new QSpinBox();
    SPBcpcoordprecision->setMinimum(0);
    SPBcpcoordprecision->setMaximum(7);
    SPBcpcoordprecision->setMaximumWidth(50);
    SPBcpcoordprecision->setValue(cps->getcpcoordprecision());
    SPBcpcoordprecision->setVisible(CHKcpcoords->isChecked());
    connections << connect(SPBcpcoordprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBcpcoordprecision_changed(int)));

//         CP symbols, CP indices

    CHKcpsymbols = new QCheckBox(tr("CP symbols"));
    CHKcpsymbols->setChecked(cps->getdrawcpssymbols());
    connections << connect(CHKcpsymbols, SIGNAL(stateChanged(int)), this, SLOT(CHKcpsymbols_changed(int)));
    CHKcpindices = new QCheckBox(tr("CP indices"));
    CHKcpindices->setChecked(cps->getdrawcpsindices());
    connections << connect(CHKcpindices, SIGNAL(stateChanged(int)), this, SLOT(CHKcpindices_changed(int)));

//         CP values

    CHKcpvalues = new QCheckBox(tr("CP field values"));
    CHKcpvalues->setChecked(cps->getdrawcpsvalues());
    connections << connect(CHKcpvalues, SIGNAL(stateChanged(int)), this, SLOT(CHKcpvalues_changed(int)));

    LBLcpprecision = new QLabel(tr("Precision"));
    LBLcpprecision->setVisible(CHKcpvalues->isChecked());
    SPBcpprecision = new QSpinBox();
    SPBcpprecision->setMinimum(0);
    SPBcpprecision->setMaximum(7);
    SPBcpprecision->setMaximumWidth(50);
    SPBcpprecision->setValue(cps->getcpprecision());
    SPBcpprecision->setVisible(CHKcpvalues->isChecked());
    connections << connect(SPBcpprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBcpprecision_changed(int)));

    CHKcpactiveonly = new QCheckBox(tr("Only selected CPs"));
    CHKcpactiveonly->setChecked(cps->getonlycpsactive());
    if (CHKcpvalues->isChecked() || CHKcpcoords->isChecked() || CHKcpsymbols->isChecked() || CHKcpindices->isChecked()){
        CHKcpactiveonly->setVisible(true);
    }
    else {
        CHKcpactiveonly->setVisible(false);
    }
    connections << connect(CHKcpactiveonly, SIGNAL(stateChanged(int)), this, SLOT(CHKcpactiveonly_changed(int)));

//         CP selection

    LBLcpselect = new QLabel(tr("Select")+":");
    QColor *cpbtnbkgcolor = new QColor(29,124,31);
    BTNcpselectall = new ColorButton();
    BTNcpselectall->setText(tr("All"));
    BTNcpselectall->setColor(cpbtnbkgcolor);
    connections << connect(BTNcpselectall, SIGNAL(clicked()), this, SLOT(BTNcpselectall_clicked()));
    cpbtnbkgcolor = new QColor(255,0,0);
    BTNcpselectnone = new ColorButton();
    BTNcpselectnone->setText(tr("None"));
    BTNcpselectnone->setColor(cpbtnbkgcolor);
    connections << connect(BTNcpselectnone, SIGNAL(clicked()), this, SLOT(BTNcpselectnone_clicked()));
    if (CHKcpactiveonly->isChecked()){
        LBLcpselect->setVisible(true);
        BTNcpselectall->setVisible(true);
        BTNcpselectnone->setVisible(true);
    }
    else {
        LBLcpselect->setVisible(false);
        BTNcpselectall->setVisible(false);
        BTNcpselectnone->setVisible(false);
    }


    BTNcplblfont = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNcplblfont, SIGNAL(clicked()), this, SLOT(BTNcplblfont_clicked()));

    QLabel *LBLcpvshift = new QLabel(tr("Vertical shift"));
    SPBcpvshift = new QSpinBox();
    SPBcpvshift->setRange(-150,150);
    SPBcpvshift->setMaximumWidth(70);
    SPBcpvshift->setValue(cps->getcpvshift());
    connections << connect(SPBcpvshift,SIGNAL(valueChanged(int)),this,SLOT(SPBcpvshift_changed(int)));

    BTNcpcolorfont = new ColorButton();
    BTNcpcolorfont->setIcon(QIcon(":/images/fonts48.png"));
    BTNcpcolorfont->setText(tr("Color"));
    QColor fontcolor = cps->getfontcolor();
    BTNcpcolorfont->setColor(&fontcolor);
    BTNcpcolorfont->setVisible(true);
    connections << connect(BTNcpcolorfont, SIGNAL(clicked()), this, SLOT(BTNcpcolorfont_clicked()));


//        Hessian eigenvectors

    FRMcpeigvec= new QGroupBox();
    FRMcpeigvec->setTitle(tr("Hessian eigenvectors"));
    if (cps->geterroreigvec()){
        cps->setvisiblecpseigen(false);
        FRMcpeigvec->setVisible(false);
    }
    else
        FRMcpeigvec->setVisible(true);

    CHKcpeigvec = new QCheckBox(tr("Display hessian eigenvectors"));
    connections << connect(CHKcpeigvec, SIGNAL(stateChanged(int)), this, SLOT(CHKcpeigvec_changed(int)));

    QLabel *LBLcolorcpeigvec0 = new QLabel(tr("Highest abs. eigenvalue"));
    QLabel *LBLcolorcpeigvec1 = new QLabel(tr("Interm. abs. eigenvalue"));
    QLabel *LBLcolorcpeigvec2 = new QLabel(tr("Lowest abs. eigenvalue"));
    for (int i=0 ; i<3 ; ++i){
        BTNcpeigcolor[i] = new ColorButton();
        BTNcpeigcolor[i]->setIcon(QIcon(":/images/colores48.png"));
        BTNcpeigcolor[i]->setText(tr("Color"));
        QColor color = cps->geteigcolor(i);
        BTNcpeigcolor[i]->setColor(&color);
        connections << connect(BTNcpeigcolor[i],SIGNAL(clicked()),this,SLOT(BTNcpeigcolor_change()));
    }

//         Thickness

    QLabel *LBLcpeigthickness = new QLabel(tr("Vectors thickness"));
    SPBcpeigthickness = new QSpinBox();
    SPBcpeigthickness->setMinimum(1);
    SPBcpeigthickness->setMaximum(50);
    SPBcpeigthickness->setMaximumWidth(50);
    SPBcpeigthickness->setValue(cps->geteigthickness());
    connections << connect(SPBcpeigthickness,SIGNAL(valueChanged(int)),this,SLOT(SPBcpeigthickness_changed(int)));

//         Vectors length

    QLabel *LBLcpeiglength = new QLabel(tr("Vectors length"));
    SPBcpeiglength = new QSpinBox();
    SPBcpeiglength->setMinimum(1);
    SPBcpeiglength->setMaximum(100);
    SPBcpeiglength->setMaximumWidth(50);
    SPBcpeiglength->setSingleStep(1);
    SPBcpeiglength->setValue(cps->geteiglength());
    connections << connect(SPBcpeiglength,SIGNAL(valueChanged(int)),this,SLOT(SPBcpeiglength_change(int)));

//         Arrowlength

    QLabel *LBLcpeigarrowsize = new QLabel(tr("Arrows length"));
    SPBcpeigarrowsize = new QSpinBox();
    SPBcpeigarrowsize->setMinimum(1);
    SPBcpeigarrowsize->setMaximum(50);
    SPBcpeigarrowsize->setMaximumWidth(50);
    SPBcpeigarrowsize->setSingleStep(1);
    SPBcpeigarrowsize->setValue(cps->geteigarrowsize());
    connections << connect(SPBcpeigarrowsize,SIGNAL(valueChanged(int)),this,SLOT(SPBcpeigarrowsize_changed(int)));


//         Arrowwidth

    QLabel *LBLcpeigarrowwidth = new QLabel(tr("Arrows width"));
    SPBcpeigarrowwidth = new QSpinBox();
    SPBcpeigarrowwidth->setMinimum(1);
    SPBcpeigarrowwidth->setMaximum(50);
    SPBcpeigarrowwidth->setMaximumWidth(50);
    SPBcpeigarrowwidth->setSingleStep(1);
    SPBcpeigarrowwidth->setValue(10);
    connections << connect(SPBcpeigarrowwidth,SIGNAL(valueChanged(int)),this,SLOT(SPBcpeigarrowwidth_changed(int)));

//    Layouts

    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(LBLcps);
    Layout1->addWidget(TXTcps);
    Layout1->addWidget(BTNcps);

    QGridLayout *Layout2=new QGridLayout();
    Layout2->addWidget(LBLcpscolor0,0,0,Qt::AlignLeft);
    Layout2->addWidget(CHKcps[0],0,1,Qt::AlignCenter);
    Layout2->addWidget(BTNcpcolor[0],0,2,Qt::AlignRight);
    Layout2->addWidget(LBLcpscolor1,1,0,Qt::AlignLeft);
    Layout2->addWidget(CHKcps[1],1,1,Qt::AlignCenter);
    Layout2->addWidget(BTNcpcolor[1],1,2,Qt::AlignRight);
    Layout2->addWidget(LBLcpscolor2,2,0,Qt::AlignLeft);
    Layout2->addWidget(CHKcps[2],2,1,Qt::AlignCenter);
    Layout2->addWidget(BTNcpcolor[2],2,2,Qt::AlignRight);
    Layout2->addWidget(LBLcpscolor3,3,0,Qt::AlignLeft);
    Layout2->addWidget(CHKcps[3],3,1,Qt::AlignCenter);
    Layout2->addWidget(BTNcpcolor[3],3,2,Qt::AlignRight);

//             Ball radius

    QGridLayout *Layout3=new QGridLayout();
    Layout3->addWidget(LBLcpsballradius,0,0,Qt::AlignLeft);
    Layout3->addWidget(SPBcpballradius,0,1,Qt::AlignRight);

//             Show/hide CP symbols and CP indices layouts

    QGridLayout *Layout6=new QGridLayout();
    Layout6->addWidget(CHKcpsymbols,0,0);
    Layout6->addWidget(CHKcpindices,1,0);
    Layout6->addWidget(CHKcpcoords,2,0,Qt::AlignLeft);
    Layout6->addWidget(LBLcpcoordprecision,2,1,Qt::AlignLeft);
    Layout6->addWidget(SPBcpcoordprecision,2,2,Qt::AlignLeft);
    Layout6->addWidget(RBTbohrcp,3,1,Qt::AlignLeft);
    Layout6->addWidget(RBTangstromcp,3,2,Qt::AlignLeft);
    Layout6->addWidget(CHKcpvalues,4,0,Qt::AlignLeft);
    Layout6->addWidget(LBLcpprecision,4,1,Qt::AlignLeft);
    Layout6->addWidget(SPBcpprecision,4,2,Qt::AlignLeft);


    QHBoxLayout *Layout10=new QHBoxLayout();
    Layout10->addWidget(CHKcpactiveonly);
    Layout10->setAlignment(Qt::AlignLeft);

    QHBoxLayout *Layout11=new QHBoxLayout();
    Layout11->addWidget(LBLcpselect);
    Layout11->addWidget(BTNcpselectall);
    Layout11->addWidget(BTNcpselectnone);

    QHBoxLayout *Layout12=new QHBoxLayout();
    Layout12->addWidget(LBLcpvshift);
    Layout12->addWidget(SPBcpvshift);

    QHBoxLayout *Layout13=new QHBoxLayout();
    Layout13->addStretch();
    Layout13->addWidget(BTNcplblfont);
    Layout13->addWidget(BTNcpcolorfont);
    Layout13->addStretch();
    Layout13->setAlignment(Qt::AlignCenter);

//              CPs layout

    QVBoxLayout *Layout14=new QVBoxLayout(FRMcps);
    Layout14->addLayout(Layout2);
    Layout14->addLayout(Layout3);
    Layout14->addLayout(Layout6);
    Layout14->addLayout(Layout10);
    Layout14->addLayout(Layout11);
    Layout14->addLayout(Layout12);
    Layout14->addLayout(Layout13);
    Layout14->addStretch();


//      Hessian eigenvectors layout

//          Color buttons layout

    QGridLayout *Layout15=new QGridLayout();
    Layout15->addWidget(LBLcolorcpeigvec0,0,0,Qt::AlignLeft);
    Layout15->addWidget(BTNcpeigcolor[0],0,2,Qt::AlignRight);
    Layout15->addWidget(LBLcolorcpeigvec1,1,0,Qt::AlignLeft);
    Layout15->addWidget(BTNcpeigcolor[1],1,2,Qt::AlignRight);
    Layout15->addWidget(LBLcolorcpeigvec2,2,0,Qt::AlignLeft);
    Layout15->addWidget(BTNcpeigcolor[2],2,2,Qt::AlignRight);

//        Arrows layout

    QGridLayout *Layout16=new QGridLayout();
    Layout16->addWidget(LBLcpeigthickness,0,0,Qt::AlignLeft);
    Layout16->addWidget(SPBcpeigthickness,0,1,Qt::AlignRight);
    Layout16->addWidget(LBLcpeiglength,1,0,Qt::AlignLeft);
    Layout16->addWidget(SPBcpeiglength,1,1,Qt::AlignRight);
    Layout16->addWidget(LBLcpeigarrowsize,2,0,Qt::AlignLeft);
    Layout16->addWidget(SPBcpeigarrowsize,2,1,Qt::AlignRight);
    Layout16->addWidget(LBLcpeigarrowwidth,3,0,Qt::AlignLeft);
    Layout16->addWidget(SPBcpeigarrowwidth,3,1,Qt::AlignRight);

    QVBoxLayout *Layout17=new QVBoxLayout(FRMcpeigvec);
    Layout17->addWidget(CHKcpeigvec,Qt::AlignLeft);
    Layout17->addLayout(Layout15);
    Layout17->addLayout(Layout16);
    Layout17->addStretch();

//    Dialog layout

    QVBoxLayout *Layout = new  QVBoxLayout(FRMcriticalpoints);
    Layout->addLayout(Layout1);
    Layout->addWidget(FRMcps);
    Layout->addWidget(FRMcpeigvec);
    Layout->addStretch();
}

void molecule::create_surfaces_widgets_and_layouts(){
//        Add surface
    BTNaddsurface = new QPushButton(tr("Add surface"));
    connections << connect(BTNaddsurface, SIGNAL(clicked()), this, SLOT(addsurface()));

    layoutsurfs = new QVBoxLayout();
    if (surfaces->count() > 0){
        QLabel *LBLavailablesurfs = new QLabel();
        LBLavailablesurfs->setText(tr("<font color=\"black\">Loaded surfaces</font>"));
        layoutsurfs->addWidget(LBLavailablesurfs,Qt::AlignCenter);
        QSignalMapper* deletesurfsignalMapper = new QSignalMapper (this) ;
        QSignalMapper* showsurfsignalMapper = new QSignalMapper (this) ;
        for (int i = 0, knt = 0 ; i < surfaces->count() ; i++){
            togglingGroupBox *FRMsurfaceeditor = surfaces_editor(surfaces->at(i));
            FRMsurfaceeditor->setVisible(false);

            QLabel *LBLsurf = new QLabel();
            LBLsurf->setText(surfaces->at(i)->getname());
            editclosePushButton *BTNeditsurf = new editclosePushButton();
            connections << connect(BTNeditsurf, SIGNAL(clicked()), BTNeditsurf,
                            SLOT(toggletext()), Qt::UniqueConnection);
            connections << connect(BTNeditsurf, SIGNAL(clicked()), FRMsurfaceeditor,
                            SLOT(toggleVisible()), Qt::UniqueConnection);
            connections << connect(FRMsurfaceeditor, SIGNAL(isvisible(bool)), this,
                            SLOT(updateQDLeditMolecule(bool)), Qt::UniqueConnection);
            connections << connect(surfaces->at(i), SIGNAL(opendialog()), this,
                            SLOT(closeisosurfeditors()), Qt::UniqueConnection);
            QPushButton *BTNdeletesurf = new QPushButton();
            BTNdeletesurf->setText(tr("Delete"));
            connections << connect(BTNdeletesurf, SIGNAL(clicked()), deletesurfsignalMapper, SLOT(map()), Qt::UniqueConnection);
            deletesurfsignalMapper -> setMapping(BTNdeletesurf,i);
            showhidePushButton *BTNshowsurf = new showhidePushButton();
            if (surfaces->at(i)->getvisible())
                BTNshowsurf->inittext(1);
            else
                BTNshowsurf->inittext(0);
            connections << connect(BTNshowsurf, SIGNAL(clicked()), showsurfsignalMapper, SLOT(map()), Qt::UniqueConnection);
            showsurfsignalMapper -> setMapping(BTNshowsurf,i);
            connections << connect(BTNshowsurf, SIGNAL(clicked()), BTNshowsurf,
                            SLOT(toggletext()), Qt::UniqueConnection);

            QGridLayout *layout0 = new QGridLayout();
            layout0->addWidget(LBLsurf,knt+1,0,1,3);
            layout0->addWidget(BTNeditsurf,knt+1,3);
            layout0->addWidget(BTNshowsurf,knt+1,4);
            layout0->addWidget(BTNdeletesurf,knt+1,5);

            layoutsurfs->addLayout(layout0);


            layoutsurfs->addWidget(FRMsurfaceeditor);

            knt++;
        }
        connections << connect (deletesurfsignalMapper, SIGNAL(mapped(int)), this,
                        SLOT(deletesurf(int)), Qt::UniqueConnection) ;
        connections << connect (showsurfsignalMapper, SIGNAL(mapped(int)), this,
                        SLOT(toggleshowsurf(int)), Qt::UniqueConnection) ;
    }
}

togglingGroupBox *molecule::surfaces_editor(surface *surf){

    togglingGroupBox *FRMsurfaceeditor = new togglingGroupBox();
//        Surface type

    QGroupBox *FRMsurftype = new QGroupBox(tr("Surface type"));
    QRadioButton *RBTsolidsurf = new QRadioButton(tr("Solid surface"),FRMsurftype);
    RBTsolidsurf->setChecked(surf->getsolidsurf());
    QRadioButton *RBTwiresurf = new QRadioButton(tr("Wire frame"),FRMsurftype);
    RBTwiresurf->setChecked(!surf->getsolidsurf());
    connections << connect(RBTsolidsurf, SIGNAL(toggled (bool)), surf, SLOT(setsolidsurf(bool)));

//        Opacity
    labelSlider *SLDopacity = new labelSlider();
    SLDopacity->setlabeltext(QString(tr("Opacity:")));
    SLDopacity->setRange(0,100);
    SLDopacity->setSingleStep(1);
    SLDopacity->setPageStep(10);
    SLDopacity->setValue(100.f * surf->getopacity());
    connections << connect(SLDopacity, SIGNAL(valueChanged(int)), surf, SLOT(opacity_changed(int)));
    connections << connect(SLDopacity, SIGNAL(sliderReleased()), surf, SLOT(opacity_released()));


//        Translucency correction
    QCheckBox *CHKtranslucence = new QCheckBox(tr("Translucence correction"));
    CHKtranslucence->setChecked(surf->gettranslucence());
    connections << connect(CHKtranslucence, SIGNAL(toggled (bool)), surf, SLOT(settranslucence(bool)));

//        Basins color
    ColorButton *BTNsurfacecolor = new ColorButton();
    BTNsurfacecolor->setIcon(QIcon(":/images/colores48.png"));
    BTNsurfacecolor->setText(tr("Color"));
    QColor surfacecolor = surf->getsurfacecolor();
    BTNsurfacecolor->setColor(&surfacecolor);
    BTNsurfacecolor->setEnabled(true);
    connect(BTNsurfacecolor,SIGNAL(clicked()),surf,SLOT(selectsurfacecolor()));
    connect(surf,SIGNAL(surfaceColor(QColor *)),BTNsurfacecolor,SLOT(setColor(QColor *)));

//        Color boundaries

    QDoubleValidator *myDoubleValidator = new QDoubleValidator();
    myDoubleValidator->setLocale(QLocale::English);
    myDoubleValidator->setRange(0,surf->getfabstop(),6);

    QGroupBox *FRMcolorbounds = new QGroupBox(tr("Color boundaries"));
    FRMcolorbounds->setVisible(surf->getshowcolorrule());

    lineEditSlider *SLDcolorbounds = new lineEditSlider();
    SLDcolorbounds->setRange(0,100);
    SLDcolorbounds->setSingleStep(1);
    SLDcolorbounds->setPageStep(10);
    SLDcolorbounds->setlabeltext(tr((QString("Range: ")+QChar(0x00B1)).toUtf8()));
    SLDcolorbounds->setValidator(myDoubleValidator);
    SLDcolorbounds->setTXTtext(QString("%1").arg(surf->gettopcolor()));

    SLDcolorbounds->setValue(surf->getfabstop());
    connections << connect(SLDcolorbounds, SIGNAL(valueChanged(float)), surf, SLOT(settopcolor(float)));
    connections << connect(SLDcolorbounds, SIGNAL(valueChanged(float)), surf, SLOT(settrianglecolors()));
    connections << connect(SLDcolorbounds, SIGNAL(sliderReleased()), surf, SLOT(settrianglecolors()));
    connections << connect(SLDcolorbounds, SIGNAL(sliderReleased()), this, SLOT(emitupdatedisplay()));

//        Show local maxima and minima
    QLabel *LBLshowlocalmax = new QLabel(tr("Local maxima"));
    QCheckBox *CHKshowlocalmax = new QCheckBox();
    CHKshowlocalmax->setChecked(surf->getshowlocalmax());

    QLabel *LBLshowlocalmin = new QLabel(tr("Local minima"));
    QCheckBox *CHKshowlocalmin = new QCheckBox();
    CHKshowlocalmin->setChecked(surf->getshowlocalmin());

    combineSignalsOr *signalshowlocal = new combineSignalsOr();
    QGroupBox *FRMlocalextremaoptions = new QGroupBox();
    FRMlocalextremaoptions->setVisible(CHKshowlocalmax->isChecked()||CHKshowlocalmin->isChecked());

    connections << connect(CHKshowlocalmax, SIGNAL(toggled(bool)), signalshowlocal, SLOT(setSignal1(bool)));
    connections << connect(CHKshowlocalmax, SIGNAL(toggled(bool)), surf, SLOT(setshowlocalmax(bool)));

    connections << connect(CHKshowlocalmin, SIGNAL(toggled(bool)), signalshowlocal, SLOT(setSignal2(bool)));
    connections << connect(CHKshowlocalmin, SIGNAL(toggled(bool)), surf, SLOT(setshowlocalmin(bool)));

    connections << connect(signalshowlocal, SIGNAL(signalsCombined(bool)), FRMlocalextremaoptions, SLOT(setVisible(bool)));
    connections << connect(signalshowlocal, SIGNAL(signalsCombined(bool)), this, SLOT(resizeQDLeditMolecule()));


//        Show symbols of local extrema

    QCheckBox *CHKshowextremasymbols = new QCheckBox(tr("Show symbols"));
    CHKshowextremasymbols->setChecked(surf->getshowextremasymbols());
    connections << connect(CHKshowextremasymbols, SIGNAL(toggled(bool)), surf, SLOT(setshowextremasymbols(bool)));

//        Show indices of local extrema

    QCheckBox *CHKshowextremaindices = new QCheckBox(tr("Show indices"));
    CHKshowextremaindices->setChecked(surf->getshowextremaindices());
    connections << connect(CHKshowextremaindices, SIGNAL(toggled(bool)), surf, SLOT(setshowextremaindices(bool)));

//        Show values of local extrema

    QCheckBox *CHKshowextremavalues = new QCheckBox(tr("Show MESP values"));
    CHKshowextremavalues->setChecked(surf->getshowextremavalues());
    connections << connect(CHKshowextremavalues, SIGNAL(toggled(bool)), surf, SLOT(setshowextremavalues(bool)));

    QLabel *LBLvalueprecision = new QLabel(tr("Precision"));
    LBLvalueprecision->setVisible(CHKshowextremavalues->isChecked());
    QSpinBox *SPBvalueprecision = new QSpinBox();
    SPBvalueprecision->setMinimum(0);
    SPBvalueprecision->setMaximum(7);
    SPBvalueprecision->setMaximumWidth(50);
    SPBvalueprecision->setValue(surf->getvalueprecision());
    SPBvalueprecision->setVisible(CHKshowextremavalues->isChecked());
    connections << connect(SPBvalueprecision,SIGNAL(valueChanged(int)),surf,SLOT(setvalueprecision(int)));
    connections << connect(CHKshowextremavalues, SIGNAL(toggled(bool)), LBLvalueprecision, SLOT(setVisible(bool)));
    connections << connect(CHKshowextremavalues, SIGNAL(toggled(bool)), SPBvalueprecision, SLOT(setVisible(bool)));
    connections << connect(CHKshowextremavalues, SIGNAL(toggled(bool)), this, SLOT(resizeQDLeditMolecule()));

//        Show coordinates of local extrema

    QCheckBox *CHKshowextremacoords = new QCheckBox(tr("Show coordinates"));
    CHKshowextremacoords->setChecked(surf->getshowextremacoords());
    connections << connect(CHKshowextremacoords, SIGNAL(toggled(bool)), surf, SLOT(setshowextremacoords(bool)));

    LBLcoordprecision = new QLabel(tr("Precision"));
    LBLcoordprecision->setVisible(CHKshowextremacoords->isChecked());
    SPBcoordprecision = new QSpinBox();
    SPBcoordprecision->setMinimum(0);
    SPBcoordprecision->setMaximum(7);
    SPBcoordprecision->setMaximumWidth(50);
    SPBcoordprecision->setValue(surf->getcoordprecision());
    SPBcoordprecision->setVisible(CHKshowextremacoords->isChecked());
    connections << connect(SPBcoordprecision,SIGNAL(valueChanged(int)),surf,SLOT(setcoordprecision(int)));
    connections << connect(CHKshowextremacoords, SIGNAL(toggled(bool)), LBLcoordprecision, SLOT(setVisible(bool)));
    connections << connect(CHKshowextremacoords, SIGNAL(toggled(bool)), SPBcoordprecision, SLOT(setVisible(bool)));
    connections << connect(CHKshowextremacoords, SIGNAL(toggled(bool)), this, SLOT(resizeQDLeditMolecule()));

//            Only selected local extrema

    QCheckBox *CHKextremactiveonly = new QCheckBox(tr("Only selected extrema"));
    CHKextremactiveonly->setChecked(surf->getonlyextremactive());
    connections << connect(CHKextremactiveonly, SIGNAL(toggled(bool)), surf, SLOT(setonlyextremactive(bool)));

    QColor *cpbtnbkgcolor = new QColor(29,124,31);
    ColorButton *BTNextremaselectall = new ColorButton();
    BTNextremaselectall->setText(tr("All"));
    BTNextremaselectall->setColor(cpbtnbkgcolor);
    BTNextremaselectall->setHidden(true);
    connections << connect(BTNextremaselectall, SIGNAL(clicked()), surf, SLOT(extremaselectall()));
    cpbtnbkgcolor = new QColor(255,0,0);
    ColorButton *BTNextremaselectnone = new ColorButton();
    BTNextremaselectnone->setText(tr("None"));
    BTNextremaselectnone->setColor(cpbtnbkgcolor);
    BTNextremaselectnone->setHidden(true);
    connections << connect(BTNextremaselectnone, SIGNAL(clicked()), surf, SLOT(extremaselectnone()));
    connections << connect(CHKextremactiveonly, SIGNAL(toggled(bool)), BTNextremaselectall, SLOT(setVisible(bool)));
    connections << connect(CHKextremactiveonly, SIGNAL(toggled(bool)), BTNextremaselectnone, SLOT(setVisible(bool)));
    connections << connect(CHKextremactiveonly, SIGNAL(toggled(bool)), this, SLOT(resizeQDLeditMolecule()));

//         Ball radius

    QLabel *LBLballradius = new QLabel(tr("Ball radius"));
    QSpinBox *SPBballradius = new QSpinBox();
    SPBballradius->setMinimum(1);
    SPBballradius->setMaximum(100);
    SPBballradius->setMaximumWidth(50);
    SPBballradius->setValue(surf->getballradius());
    connections << connect(SPBballradius,SIGNAL(valueChanged(int)),surf,SLOT(setballradius(int)));

//         Extrema colors

    ColorButton *BTNmaximacolor = new ColorButton();
    BTNmaximacolor->setIcon(QIcon(":/images/colores48.png"));
    BTNmaximacolor->setText(tr("Color"));
    QColor maximacolor = surf->getextremacolor(0);
    BTNmaximacolor->setColor(&maximacolor);
    connect(BTNmaximacolor,SIGNAL(clicked()),surf,SLOT(selectmaximacolor()));
    connect(surf,SIGNAL(maximaColor(QColor *)),BTNmaximacolor,SLOT(setColor(QColor *)));

    ColorButton *BTNminimacolor = new ColorButton();
    BTNminimacolor->setIcon(QIcon(":/images/colores48.png"));
    BTNminimacolor->setText(tr("Color"));
    QColor minimacolor = surf->getextremacolor(1);
    BTNminimacolor->setColor(&minimacolor);
    connect(BTNminimacolor,SIGNAL(clicked()),surf,SLOT(selectminimacolor()));
    connect(surf,SIGNAL(minimaColor(QColor *)),BTNminimacolor,SLOT(setColor(QColor *)));

//        Vertical shift

    QLabel *LBLvshift = new QLabel(tr("Vertical shift"));
    QSpinBox *SPBvshift = new QSpinBox();
    SPBvshift->setRange(-150,150);
    SPBvshift->setMaximumWidth(70);
    SPBvshift->setValue(surf->getvshift());
    connections << connect(SPBvshift,SIGNAL(valueChanged(int)),surf,SLOT(setvshift(int)));

//         Extrema labels format

    QLabel *LBLformat = new QLabel(tr("Labels format")+":");

    QPushButton *BTNextremafont = new QPushButton();
    BTNextremafont = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNextremafont, SIGNAL(clicked()), surf, SLOT(changeextremafont()));
    ColorButton *BTNextremacolorfont = new ColorButton();
    BTNextremacolorfont->setIcon(QIcon(":/images/fonts48.png"));
    BTNextremacolorfont->setText(tr("Color"));
    BTNextremacolorfont->setColor(&fontcolor);
    BTNextremacolorfont->setEnabled(true);
    connections << connect(BTNextremacolorfont, SIGNAL(clicked()), surf, SLOT(selectfontcolor()));
    connect(surf,SIGNAL(fontColor(QColor *)),BTNextremacolorfont,SLOT(setColor(QColor *)));

//        Show grid
    QCheckBox *CHKshowgrid = new QCheckBox(tr("Show grid boundaries"));
    CHKshowgrid->setChecked(surf->getshowgridbounds());
    connections << connect(CHKshowgrid, SIGNAL(toggled(bool)), surf, SLOT(setshowgridbounds(bool)));

    QGroupBox *FRMlocalextrema = new QGroupBox();
    if (QFileInfo(surf->getname()).suffix() == "sgh" || QFileInfo(surf->getname()).suffix() == "srf"){
        FRMlocalextrema->setVisible(true);
        BTNsurfacecolor->setVisible(false);
    }
    else{
        FRMlocalextrema->setVisible(false);
        CHKshowgrid->setVisible(false);
        BTNsurfacecolor->setVisible(true);
    }

    QHBoxLayout *layout1 = new QHBoxLayout();
    layout1->addStretch();
    layout1->addWidget(RBTsolidsurf);
    layout1->addWidget(RBTwiresurf);
    layout1->addStretch();

    QHBoxLayout *layout4 = new QHBoxLayout();
    layout4->addStretch();
    layout4->addWidget(CHKtranslucence);
    layout4->addStretch();

    QVBoxLayout *layout5 = new QVBoxLayout(FRMsurftype);
    layout5->addLayout(layout1);
    layout5->addWidget(SLDopacity);
    layout5->addLayout(layout4);

    QVBoxLayout *layout6 = new QVBoxLayout(FRMcolorbounds);
    layout6->addWidget(SLDcolorbounds);
    layout6->addStretch();

    QHBoxLayout *layout7 = new QHBoxLayout();
    layout7->addWidget(CHKshowgrid);

    QGridLayout *layout8 = new QGridLayout();
    layout8->addWidget(LBLshowlocalmax,0,0,Qt::AlignLeft);
    layout8->addWidget(CHKshowlocalmax,0,1,Qt::AlignCenter);
    layout8->addWidget(BTNmaximacolor,0,2,Qt::AlignRight);
    layout8->addWidget(LBLshowlocalmin,1,0,Qt::AlignLeft);
    layout8->addWidget(CHKshowlocalmin,1,1,Qt::AlignCenter);
    layout8->addWidget(BTNminimacolor,1,2,Qt::AlignRight);

    QHBoxLayout *layout9 = new QHBoxLayout();
    layout9->addStretch();
    layout9->addWidget(LBLballradius);
    layout9->addWidget(SPBballradius);
    layout9->addStretch();

    QVBoxLayout *layout10 = new QVBoxLayout(FRMlocalextrema);
    layout10->addLayout(layout8);
    layout10->addLayout(layout9);
    layout10->addStretch();

    QGridLayout * layout12 = new QGridLayout();
    layout12->addWidget(CHKshowextremasymbols,0,0);
    layout12->addWidget(CHKshowextremaindices,1,0);
    layout12->addWidget(CHKshowextremavalues,2,0);
    layout12->addWidget(LBLvalueprecision,2,1,Qt::AlignRight);
    layout12->addWidget(SPBvalueprecision,2,2,Qt::AlignRight);
    layout12->addWidget(CHKshowextremacoords,3,0);
    layout12->addWidget(LBLcoordprecision,3,1,Qt::AlignRight);
    layout12->addWidget(SPBcoordprecision,3,2,Qt::AlignRight);
    layout12->addWidget(CHKextremactiveonly,4,0);
    layout12->addWidget(BTNextremaselectall,4,1);
    layout12->addWidget(BTNextremaselectnone,4,2);
    layout12->addWidget(LBLformat,5,0);
    layout12->addWidget(BTNextremafont,5,1);
    layout12->addWidget(BTNextremacolorfont,5,2);
    layout12->addWidget(LBLvshift,6,0);
    layout12->addWidget(SPBvshift,6,1,Qt::AlignLeft);

    QVBoxLayout * layout17 = new QVBoxLayout(FRMlocalextremaoptions);
    layout17->addLayout( layout12);
    layout17->addStretch();

    QHBoxLayout * layout18 = new QHBoxLayout();
    layout18->addStretch();
    layout18->addWidget(BTNsurfacecolor);
    layout18->addStretch();

    QVBoxLayout *layout = new QVBoxLayout(FRMsurfaceeditor);
    layout->addStretch();
    layout->addWidget(FRMsurftype);
    layout->addLayout( layout18);
    layout->addWidget(FRMcolorbounds);
    layout->addLayout(layout7);
    layout->addWidget(FRMlocalextrema);
    layout->addWidget(FRMlocalextremaoptions);

    return FRMsurfaceeditor;
}

// Function getscalevalueInt
//    Returns slider position ix corresponding to function value fx.
//      i0, i1: range of slider scale
//      f0, f1: function values at i0 and i1
int molecule::getscalevalueInt(float fx,int i0,int i1,float f0,float f1)
{
    if (fx > f1) fx = f1;
    if (fx < f0) fx = f0;
    return (int)(((fx-f0)*(i1-i0)/(f1-f0))+i0);
}

void molecule::create_grids_widgets_and_layouts(){
//        Add grid for isosurface
    int knt = 0 ; ;
    BTNaddgrid = new QPushButton(tr("Add grid for isosurfaces"));
    connections << connect(BTNaddgrid, SIGNAL(clicked()), this, SLOT(addgrid()));
    LBLloadinggrid = new QLabel("Loading grid");
    LBLloadinggrid->setStyleSheet("QLabel { color : red; }");
    LBLloadinggrid->setVisible(false);
    QSizePolicy sp_retain = LBLloadinggrid->sizePolicy();
    sp_retain.setRetainSizeWhenHidden(true);
    LBLloadinggrid->setSizePolicy(sp_retain);
    layoutgrids = new QGridLayout();
    if (grids->count() > 0){
        QLabel *LBLavailablegrids = new QLabel();
        LBLavailablegrids->setText(tr("<font color=\"black\">Loaded grids</font>"));
        QSignalMapper* deletegridsignalMapper = new QSignalMapper (this) ;
        layoutgrids->addWidget(LBLavailablegrids,0,0,1,5,Qt::AlignCenter);
        for (int i = 0 ; i < grids->count() ; i++){
            QLabel *LBLgrid = new QLabel();
            LBLgrid->setText(grids->at(i)->getname());
            LBLgrid->setStyleSheet("QLabel { color : black; }");
            QPushButton *BTNaddisosurf = new QPushButton();
            BTNaddisosurf->setText(tr("Add surface"));
            connections << connect(BTNaddisosurf, SIGNAL(clicked()), this, SLOT(closeisosurfeditors()), Qt::UniqueConnection);
            connections << connect(BTNaddisosurf, SIGNAL(clicked()), grids->at(i), SLOT(addisosurf()), Qt::UniqueConnection);
            QPushButton *BTNdeletegrid = new QPushButton();
            BTNdeletegrid->setText(tr("Delete"));
            connections << connect(BTNdeletegrid, SIGNAL(clicked()), deletegridsignalMapper, SLOT(map()), Qt::UniqueConnection);
            deletegridsignalMapper -> setMapping(BTNdeletegrid,i);
            layoutgrids->addWidget(LBLgrid,knt+1,0,1,2);
            layoutgrids->addWidget(BTNaddisosurf,knt+1,2);
            layoutgrids->addWidget(BTNdeletegrid,knt+1,3);
            knt++;
            QSignalMapper* delisosurfsignalMapper = new QSignalMapper (this) ;
            QSignalMapper* genisosurfsignalMapper = new QSignalMapper (this) ;
            QSignalMapper* showisosurfsignalMapper = new QSignalMapper (this) ;
            QSignalMapper* updateisosurfsignalMapper = new QSignalMapper (this) ;
            for (int j = 0 ; j < grids->at(i)->surfaces->count() ; j++){
                LBLavailablegrids->setText(tr("<font color=\"black\">Loaded grids</font> <font color=\"blue\">(Surfaces)</font>"));
                QLabel *LBLisosurf = new QLabel();
                LBLisosurf->setText(grids->at(i)->surfaces->at(j)->getname());
                LBLisosurf->setStyleSheet("color : blue");
                QString lblstring;
                int r = grids->at(i)->surfaces->at(j)->getsurfcolor().red();
                int g = grids->at(i)->surfaces->at(j)->getsurfcolor().green();
                int b = grids->at(i)->surfaces->at(j)->getsurfcolor().blue();
                if ((r > 210 || b > 0)&& g > 210 ){
                    lblstring = QString("QPushButton { color : rgb(0, 0, 0); background-color : rgb(%1,%2,%3); }").arg(r).arg(g).arg(b);
                }
                else{
                    lblstring = QString("QPushButton { color : rgb(255, 255, 255); background-color : rgb(%1,%2,%3); }").arg(r).arg(g).arg(b);
                }

                togglingGroupBox *FRMisosurface = new togglingGroupBox();
                FRMisosurface = grids->at(i)->surfaces->at(j)->editisosurface();
                FRMisosurface->setVisible(false);
                FRMisosurface->setAttribute(Qt::WA_DeleteOnClose);

                editclosePushButton *BTNeditisosurf = new editclosePushButton();
                BTNeditisosurf->setText(tr("Edit"));
                BTNeditisosurf->setStyleSheet(QString(lblstring));

                connections << connect(BTNeditisosurf, SIGNAL(clicked()), BTNeditisosurf,
                                SLOT(toggletext()), Qt::UniqueConnection);
                connections << connect(BTNeditisosurf, SIGNAL(isClose(bool)), FRMisosurface,
                                SLOT(setVisible(bool)), Qt::UniqueConnection);
                connections << connect(BTNeditisosurf, SIGNAL(isClose(bool)), this,
                                SLOT(updateQDLeditMolecule(bool)), Qt::UniqueConnection);

                QPushButton *BTNdeleteisosurf = new QPushButton();
                BTNdeleteisosurf->setText(tr("Delete"));
                BTNdeleteisosurf->setStyleSheet(QString(lblstring));
                connections << connect(BTNdeleteisosurf, SIGNAL(clicked()), delisosurfsignalMapper, SLOT(map()), Qt::UniqueConnection);
                delisosurfsignalMapper -> setMapping(BTNdeleteisosurf,j);

                showhidePushButton *BTNshowisosurf = new showhidePushButton();
                if (grids->at(i)->surfaces->at(j)->isvisible())
                    BTNshowisosurf->inittext(1);
                else
                    BTNshowisosurf->inittext(0);
                BTNshowisosurf->setStyleSheet(QString(lblstring));
                connections << connect(BTNshowisosurf, SIGNAL(clicked()), showisosurfsignalMapper, SLOT(map()), Qt::UniqueConnection);
                showisosurfsignalMapper -> setMapping(BTNshowisosurf,j);
                connections << connect(BTNshowisosurf, SIGNAL(clicked()), BTNshowisosurf,
                                SLOT(toggletext()), Qt::UniqueConnection);

                connections << connect(grids->at(i)->surfaces->at(j), SIGNAL(generatesurface()), genisosurfsignalMapper,
                                SLOT(map()), Qt::UniqueConnection);
                genisosurfsignalMapper->setMapping(grids->at(i)->surfaces->at(j),j);

                connections << connect(grids->at(i)->surfaces->at(j), SIGNAL(updatedisplay()), updateisosurfsignalMapper,
                                SLOT(map()), Qt::UniqueConnection);
                updateisosurfsignalMapper->setMapping(grids->at(i)->surfaces->at(j),j);

                connections << connect(grids->at(i)->surfaces->at(j), SIGNAL(updatelabelcolor(QString)), BTNeditisosurf,
                                SLOT(setStyleSheet(QString)), Qt::UniqueConnection);
                connections << connect(grids->at(i)->surfaces->at(j), SIGNAL(updatelabelcolor(QString)), BTNdeleteisosurf,
                                SLOT(setStyleSheet(QString)), Qt::UniqueConnection);
                connections << connect(grids->at(i)->surfaces->at(j), SIGNAL(updatelabelcolor(QString)), BTNshowisosurf,
                                SLOT(setStyleSheet(QString)), Qt::UniqueConnection);

                layoutgrids->addWidget(LBLisosurf,knt+1,0,1,2);
                layoutgrids->addWidget(BTNeditisosurf,knt+1,2);
                layoutgrids->addWidget(BTNshowisosurf,knt+1,3);
                layoutgrids->addWidget(BTNdeleteisosurf,knt+1,4);
                layoutgrids->addWidget(FRMisosurface,++knt+1,0,1,5);
                knt++;
            }
            connections << connect (delisosurfsignalMapper, SIGNAL(mapped(int)), this,
                            SLOT(closeisosurfeditors()), Qt::UniqueConnection) ;
            connections << connect (delisosurfsignalMapper, SIGNAL(mapped(int)), grids->at(i),
                            SLOT(deletesurf(int)), Qt::UniqueConnection) ;
            connections << connect (delisosurfsignalMapper, SIGNAL(mapped(int)), this, SLOT(emitupdatedisplay()), Qt::UniqueConnection) ;
            connections << connect (genisosurfsignalMapper, SIGNAL(mapped(int)), grids->at(i),
                            SLOT(generatesurf(int)), Qt::UniqueConnection) ;
            connections << connect (showisosurfsignalMapper, SIGNAL(mapped(int)), grids->at(i),
                            SLOT(toggleshowsurf(int)), Qt::UniqueConnection) ;
            connections << connect (updateisosurfsignalMapper, SIGNAL(mapped(int)), this,
                            SLOT(emitupdatedisplay()), Qt::UniqueConnection) ;
        }
        connections << connect (deletegridsignalMapper, SIGNAL(mapped(int)), this, SLOT(deletegrid(int)), Qt::UniqueConnection) ;
    }
    layoutgrids->addWidget(LBLloadinggrid,knt+1,0,1,5);
}


//  ------------------------------------------------------------------------------------------------------------------
//
//          Slots and functions for molecules editor
//
//  ------------------------------------------------------------------------------------------------------------------

//  Function addcriticalpoints: open menu for critical points
//
void molecule::addcriticalpoints(){
    if (!cps){
        cps = new criticalpoints(path);
        cps->set_ProjectFolder(ProjectFolder);
        connect(cps,SIGNAL(updatedisplay()),this,SIGNAL(updatedisplay()));
    }
}
//  End of function addcriticalpoints

//  Function addfield: open menu for field lines
//
void molecule::addfield(){
    if (!flines){
        flines = new fieldlines(path);
        flines->set_ProjectFolder(ProjectFolder);
        connect(flines,SIGNAL(updatedisplay()),this,SIGNAL(updatedisplay()));
    }
}
//  End of function addfield

//  Function addforces: open menu for Hellmann-Feyman forces
//
void molecule::addforces(){
    if (!hfforces){
        hfforces = new forces(path);
        hfforces->set_ProjectFolder(ProjectFolder);
        connect(hfforces,SIGNAL(updatedisplay()),this,SIGNAL(updatedisplay()));
    }
}
//  End of function addforces

//  Function addgrid: attaches a grid to the molecule for isosurfaces plotting
//
void molecule::addgrid(){
    closeisosurfeditors();
    grids->append(new grid());
    if (!loadgrid()){
        grids->removeLast();
        LBLloadinggrid->setVisible(false);
    }
    else{
        LBLloadinggrid->setVisible(false);
        grids->last()->set_ProjectFolder(ProjectFolder);
        grids->last()->set_ProjectName(ProjectName);
        if (QDLeditMolecule)
            updateeditMoleculeDialog();
        connect (grids->last(), SIGNAL(surfaceadded()), this, SLOT(updateeditMoleculeDialog())) ;
        connect (grids->last(), SIGNAL(surfacedeleted()), this, SLOT(updateeditMoleculeDialog())) ;
    }

}
//  End of function addgrid

//  Function addsurface: attaches a surface to the molecule for 3D plotting
//
void molecule::addsurface(){
    closeisosurfeditors();
    surfaces->append(new surface());
    if (!loadsurf()){
        delete surfaces->last();
        surfaces->removeLast();
    }
    else{
        connect (surfaces->last(), SIGNAL(updatedisplay()), this, SLOT(emitupdatedisplay()));
        if (QDLeditMolecule)
            updateeditMoleculeDialog();
    }
    emitupdatedisplay();
}


void molecule::animaterotation(){
    if (rotatex){
        setrotationAxis((getrotationAxis() + deltaAngles * QVector3D(1.,0.,0.)).normalized());
        setrotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), deltaAngles) * getrotation());
    }
    if (rotatey){
        setrotationAxis((getrotationAxis() + deltaAngles * QVector3D(0.,1.,0.)).normalized());
        setrotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), deltaAngles) * getrotation());
    }
    if (rotatez){
        setrotationAxis((getrotationAxis() + deltaAngles * QVector3D(0.,0.,1.)).normalized());
        setrotation(QQuaternion::fromAxisAndAngle(getrotationAxis(), deltaAngles) * getrotation());
    }
    emit updateGL();
}

//  Function BTNanimation_clicked: toggles animation
//
void molecule::BTNanimation_clicked(){
    if (!(CHKrotatex->isChecked() || CHKrotatey->isChecked() || CHKrotatez->isChecked())){
        rotatex = false;
        rotatey = false;
        rotatez = false;
        startanimation = false;
        timer->stop();
        BTNanimation->setText(tr("Start"));
        return;
    }
    startanimation = !startanimation;
    if (startanimation){
        BTNanimation->setText(tr("Stop"));
        emit animate(true);
    }
    else{
        BTNanimation->setText(tr("Start"));
        emit animate(false);
    }
    if (startanimation)
        timer->start(interval);
    else
        timer->stop();
}
//  End of function BTNanimation_clicked

//  Function BTNaxes_clicked: open menu for axes display
//
void molecule::BTNaddaxes_clicked(){
    FRMaxes->setVisible(!FRMaxes->isVisible());
    if (FRMaxes->isVisible()){
        BTNaddaxes->setStyleSheet("QPushButton {color: red;}");
        BTNaddaxes->setText(tr("Hide axes manager"));
    }
    else{
        BTNaddaxes->setStyleSheet("QPushButton {color: black;}");
        BTNaddaxes->setText(tr("Axes"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNaxes_clicked

//  Function BTNaddcriticalpoints_clicked: open menu for HF forces display
//
void molecule::BTNaddcriticalpoints_clicked(){
    if (!cps)
        addcriticalpoints();
    FRMcriticalpoints->setVisible(!FRMcriticalpoints->isVisible());
    if (FRMcriticalpoints->isVisible()){
        BTNaddcriticalpoints->setStyleSheet("QPushButton {color: red;}");
        BTNaddcriticalpoints->setText(tr("Hide critical points manager"));
    }
    else{
        BTNaddcriticalpoints->setStyleSheet("QPushButton {color: black;}");
        BTNaddcriticalpoints->setText(tr("Critical points"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNaddcriticalpoints_clicked

//  Function BTNaddfieldlines_clicked: open menu for field lines display
//
void molecule::BTNaddfieldlines_clicked(){
    if (!flines)
        addfield();
    FRMfield->setVisible(!FRMfield->isVisible());
    if (FRMfield->isVisible()){
        BTNaddfieldlines->setStyleSheet("QPushButton {color: red;}");
        BTNaddfieldlines->setText(tr("Hide 3D lines manager"));
    }
    else{
        BTNaddfieldlines->setStyleSheet("QPushButton {color: black;}");
        BTNaddfieldlines->setText(tr("3D lines"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNaddfieldlines_clicked

//  Function BTNaddforces_clicked: open menu for HF forces display
//
void molecule::BTNaddforces_clicked(){
    if (!hfforces)
        addforces();
    FRMforces->setVisible(!FRMforces->isVisible());
    if (FRMforces->isVisible()){
        BTNaddforces->setStyleSheet("QPushButton {color: red;}");
        BTNaddforces->setText(tr("Hide HF forces manager"));
    }
    else{
        BTNaddforces->setStyleSheet("QPushButton {color: black;}");
        BTNaddforces->setText(tr("Hellmann-Feynman forces"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNaddforces_clicked


void molecule::BTNcpcolor_change(){
    QPushButton *BTNp = (QPushButton *)sender();
    int tipo=-1;
    for(int i = 0; i < MAX_CPS; ++i){
        if(BTNcpcolor[i] == BTNp){
                tipo = i;
                break;
        }
    }
    if (tipo>-1){
        QColor colact = cps->getcolor(tipo);
        QColor col = QColorDialog::getColor(colact, this);
        if(col.isValid()) {
            cps->setcolor(tipo,col);
            BTNcpcolor[tipo]->setColor(&col);
            cps->createCPs();
            emit updatedisplay();
        }
    }
}

void molecule::BTNcpcolorfont_clicked(){
    QColor col = QColorDialog::getColor(fontcolor, this);
    if(col.isValid()) {
        cps->setfontcolor(col);
        BTNcpcolorfont->setColor(&col);
        emit updatedisplay();
    }
}

void molecule::BTNcpeigcolor_change(){
    QPushButton *BTNp = (QPushButton *)sender();
    int tipo=-1;
    for(int i = 0; i < 3; ++i){
        if(BTNcpeigcolor[i] == BTNp){
                tipo = i;
                break;
        }
    }
    if (tipo>-1){
        QColor colact = cps->geteigcolor(tipo);
        QColor col = QColorDialog::getColor(colact, this);
        if(col.isValid()) {
            cps->seteigcolor(tipo,col);
            BTNcpeigcolor[tipo]->setColor(&col);
            cps->createCPseigen();
            emit updatedisplay();
        }
    }
}

void molecule::BTNcplblfont_clicked(){
    font = QFontDialog::getFont(nullpointer, font);
    cps->setfont(font);
    emit updatedisplay();
}

void molecule::BTNcpselectall_clicked(){
    for (int i=0 ; i<MAX_CPS ; i++){
        for (int j = 0 ; j < cps->getcpsactivelength(i) ; j++){
            cps->setcpsactive(i, j, true);
        }
    }
    emit updatedisplay();
}

void molecule::BTNcpselectnone_clicked(){
    for (int i=0 ; i<MAX_CPS ; i++){
        for (int j = 0 ; j < cps->getcpsactivelength(i) ; j++){
            cps->setcpsactive(i, j, false);
        }
    }
    emit updatedisplay();
}

void molecule::BTNflinescolor_clicked(){
    QColor currcolor = flines->getlinescolor();
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        flines->setlinescolor(newcolor);
        int r = newcolor.red();
        int g = newcolor.green();
        int b = newcolor.blue();
        if ((r > 210 || b > 0)&& g > 210 ){
            QColor textcolor(0,0,0);
            BTNflinescolor->setColor(&newcolor, &textcolor);
        }
        else{
            QColor textcolor(255,255,255);
            BTNflinescolor->setColor(&newcolor, &textcolor);
        }
    }
    emit updatedisplay();
}

void molecule::BTNfontcolor_clicked(){
    QColor currcolor = fontcolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        setfontcolor(newcolor);
        BTNfontcolor->setColor(&newcolor);
    }
    if (visible)
        emit updatedisplay();
    return;
}

void molecule::BTNfont_clicked(){
    setfont(QFontDialog::getFont(nullpointer, font));
    if (visible)
        emit updatedisplay();
    return;
}

void molecule::BTNfontaxeslabels_clicked(){
    setfontaxeslabels(QFontDialog::getFont(nullpointer, fontaxeslabels));
    if (axeslabels_visible)
        emit updatedisplay();
    return;
}

//    Connects push button with CPs colors
void molecule::BTNforcescolor_changed()
{
    QPushButton *BTNp = (QPushButton *)sender();
    int tipo=-1;
    for(int i = 0; i < MAX_FORCES; ++i){
        if(BTNforcecolors[i] == BTNp){
                tipo = i;
                break;
        }
    }
    if (tipo>-1){
        QColor colact = hfforces->getcolor(tipo);
        QColor col = QColorDialog::getColor(colact, this);
        if(col.isValid()) {
            hfforces->setcolor(tipo,col);
            BTNforcecolors[tipo]->setColor(&col);
        }
    }
    hfforces->createforces();
}

//  Function BTNrotation_clicked: open menu for rotations
//
void molecule::BTNrotation_clicked(){
    FRMrotation->setVisible(!FRMrotation->isVisible());
    if (FRMrotation->isVisible()){
        BTNrotation->setStyleSheet("QPushButton {color: red;}");
        BTNrotation->setText(tr("Hide rotations manager"));
    }
    else{
        BTNrotation->setStyleSheet("QPushButton {color: black;}");
        BTNrotation->setText(tr("Rotations"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNrotation_clicked


void molecule::BTNselectall_clicked(){
    if (visible){
        setallatomactive(true);
        emit updatedisplay();
    }
    return;
}

void molecule::BTNselectnone_clicked(){
    if (visible)
        setallatomactive(false);
        emit updatedisplay();
    return;
}

//  Function BTNskeleton_clicked: open menu for molecule skeleton
//
void molecule::BTNskeleton_clicked(){
    FRMskeleton->setVisible(!FRMskeleton->isVisible());
    if (FRMskeleton->isVisible()){
        BTNskeleton->setStyleSheet("QPushButton {color: red;}");
        BTNskeleton->setText(tr("Hide molecular skeleton manager"));
    }
    else{
        BTNskeleton->setStyleSheet("QPushButton {color: black;}");
        BTNskeleton->setText(tr("Molecular skeleton"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNskeleton_clicked

//  Function BTNsymbols_clicked: open menu for symbols and indices
//
void molecule::BTNsymbols_clicked(){
    FRMsymbols->setVisible(!FRMsymbols->isVisible());
    if (FRMsymbols->isVisible()){
        BTNsymbols->setStyleSheet("QPushButton {color: red;}");
        BTNsymbols->setText(tr("Hide labels manager"));
    }
    else{
        BTNsymbols->setStyleSheet("QPushButton {color: black;}");
        BTNsymbols->setText(tr("Labels"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNsymbols_clicked

//  Function BTNtranslation_clicked: open menu for translations
//
void molecule::BTNtranslation_clicked(){
    FRMtranslation->setVisible(!FRMtranslation->isVisible());
    if (FRMtranslation->isVisible()){
        BTNtranslation->setStyleSheet("QPushButton {color: red;}");
        BTNtranslation->setText(tr("Hide translations manager"));
    }
    else{
        BTNtranslation->setStyleSheet("QPushButton {color: black;}");
        BTNtranslation->setText(tr("Translations"));
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    scrollArea->updatesize(QDLeditMolecule->size());
    scrollArea->raise();
}
//  End of function BTNtranslation_clicked


void molecule::BTNXaxiscolor_clicked(){
    QColor currcolor = Xaxis_color;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        Xaxis_color = newcolor;
        BTNXaxiscolor->setColor(&newcolor);
    }
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
    return;
}

void molecule::BTNYaxiscolor_clicked(){
    QColor currcolor = Yaxis_color;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        Yaxis_color = newcolor;
        BTNYaxiscolor->setColor(&newcolor);
    }
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
    return;
}

void molecule::BTNZaxiscolor_clicked(){
    QColor currcolor = Zaxis_color;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        Zaxis_color = newcolor;
        BTNZaxiscolor->setColor(&newcolor);
    }
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
    return;
}

void molecule::CHKactiveonly_changed(int a){
    if (a == 0){
        setonlyatomactive(false);
        LBLselect->setVisible(false);
        BTNselectall->setVisible(false);
        BTNselectnone->setVisible(false);
    }
    else{
        setonlyatomactive(true);
        LBLselect->setVisible(true);
        BTNselectall->setVisible(true);
        BTNselectnone->setVisible(true);
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    if (visible)
        emit updatedisplay();
}

void molecule::CHKcpactiveonly_changed(int a){
    if (a == 0){
        cps->setonlycpsactive(false);
        LBLcpselect->setVisible(false);
        BTNcpselectall->setVisible(false);
        BTNcpselectnone->setVisible(false);
    }
    else{
        cps->setonlycpsactive(true);
        LBLcpselect->setVisible(true);
        BTNcpselectall->setVisible(true);
        BTNcpselectnone->setVisible(true);
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    emit updatedisplay();
}

void molecule::CHKcpcoords_changed(int a){
    if (a == 0){
        LBLcpcoordprecision->setVisible(false);
        SPBcpcoordprecision->setVisible(false);
        cps->setdrawcpscoords(false);
        if (!cps->getdrawcpssymbols() && !cps->getdrawcpsindices() && !cps->getdrawcpsvalues())
            CHKcpactiveonly->setVisible(false);
        RBTangstromcp->setVisible(false);
        RBTbohrcp->setVisible(false);
    }
    else{
        LBLcpcoordprecision->setVisible(true);
        SPBcpcoordprecision->setVisible(true);
        cps->setdrawcpscoords(true);
        CHKcpactiveonly->setVisible(true);
        RBTangstromcp->setVisible(true);
        RBTbohrcp->setVisible(true);
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    emit updatedisplay();
}

void molecule::CHKcpeigvec_changed(int a){
    if (a == 0){
        cps->setvisiblecpseigen(false);
    }
    else{
        cps->setvisiblecpseigen(true);
    }
    emit updatedisplay();
}

void molecule::CHKcpindices_changed(int a){
    if (a == 0){
        cps->setdrawcpsindices(false);
        if (!cps->getdrawcpscoords() && !cps->getdrawcpssymbols() && !cps->getdrawcpsvalues())
            CHKcpactiveonly->setVisible(false);
    }
    else{
        cps->setdrawcpsindices(true);
        CHKcpactiveonly->setVisible(true);
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    emit updatedisplay();
}

void molecule::CHKcps_changed(){
    QCheckBox *CHKp = (QCheckBox *)sender();
    int tipo=-1;
    for(int i = 0; i < MAX_CPS; ++i){
        if(CHKcps[i] == CHKp){
                tipo = i;
                break;
        }
    }
    if(CHKcps[tipo]->isChecked()){
        cpschecked[tipo] = true;
        cps->setvisiblecps(tipo,true);
    }else{
        cpschecked[tipo] = false;
        cps->setvisiblecps(tipo,false);
    }
    emit updatedisplay();
}

void molecule::CHKcpsymbols_changed(int a){
    if (a == 0){
        cps->setdrawcpssymbols(false);
        if (!cps->getdrawcpscoords() && !cps->getdrawcpsindices() && !cps->getdrawcpsvalues())
            CHKcpactiveonly->setVisible(false);
    }
    else{
        cps->setdrawcpssymbols(true);
        CHKcpactiveonly->setVisible(true);
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    emit updatedisplay();
}

void molecule::CHKcpvalues_changed(int a){
    if (a == 0){
        LBLcpprecision->setVisible(false);
        SPBcpprecision->setVisible(false);
        cps->setdrawcpsvalues(false);
        if (!cps->getdrawcpscoords() && !cps->getdrawcpssymbols() && !cps->getdrawcpsindices())
            CHKcpactiveonly->setVisible(false);
    }
    else{
        LBLcpprecision->setVisible(true);
        SPBcpprecision->setVisible(true);
        cps->setdrawcpsvalues(true);
        CHKcpactiveonly->setVisible(true);
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    emit updatedisplay();
}

//	Connects Checkbox with show/hide arrows in electric field lines
void molecule::CHKflinesarrows_changed(){
    if(CHKflinesarrows->isChecked()){
        if (flines->allarrowsvertices.isEmpty())
            flines->createarrows();
        flines->setshowarrows(true);
    }
    else{
        flines->setshowarrows(false);
    }
    emit updatedisplay();
}

//	Connects Checkbox with show/hide electric field lines
void molecule::CHKflines_changed(){
    if(CHKflines->isChecked())
        flines->setvisible(true);
    else
        flines->setvisible(false);
    emit updatedisplay();
}

//    Connects Checkbox with show/hide critical points
void molecule::CHKforces_changed(){
    QCheckBox *CHKp = (QCheckBox *)sender();
    int tipo=-1;
    for(int i = 0; i < MAX_FORCES; ++i){
        if(CHKforces[i] == CHKp){
                tipo = i;
                break;
        }
    }
    if(CHKforces[tipo]->isChecked()){
        hfforces->setvisibleforces(tipo,true);
    }else{
        hfforces->setvisibleforces(tipo,false);
    }
    hfforces->setvisible(false);
    for (int i = 0 ; i < MAX_FORCES ; i++){
        if (CHKforces[i]->isChecked()){
            hfforces->setvisible(true);
            break;
        }
    }
    hfforces->createforces();
}

void molecule::CHKhideatoms_changed(){
    hideatoms = CHKhideatoms->isChecked();
    makeStructureBondsandSticks();
    emit updatedisplay();
}

void molecule::CHKhidebonds_changed(){
    hidebonds = CHKhidebonds->isChecked();
    makeStructureBondsandSticks();
    emit updatedisplay();
}

void molecule::CHKhidehydrogens_changed(){
    hidehydrogens = CHKhidehydrogens->isChecked();
    makeStructureBondsandSticks();
    emit updatedisplay();
}

void molecule::CHKrotate_changed(){
    rotatex = CHKrotatex->isChecked();
    rotatey = CHKrotatey->isChecked();
    rotatez = CHKrotatez->isChecked();
    if (!(rotatex || rotatey || rotatez)){
        startanimation = false;
        timer->stop();
        BTNanimation->setText(tr("Start"));
        return;
    }
}

void molecule::CHKshowaxes_changed(int a){
    if (a == 0){
        axes_visible = false;
        CHKshowaxeslabels->setChecked(0);
    }
    else{
        make_axes();
        axes_visible = true;
    }
    CHKshowaxeslabels->setEnabled(axes_visible);
    emit updatedisplay();
}

void molecule::CHKshowaxeslabels_changed(int a){
    if (a == 0){
        axeslabels_visible = false;
    }
    else{
        make_axes();
        axeslabels_visible = true;
    }
    emit updatedisplay();
}

void molecule::CHKshowcoords_changed(int a){
    if (a == 0){
        LBLcoordprecision->setVisible(false);
        SPBcoordprecision->setVisible(false);
        RBTangstromcoor->setVisible(false);
        RBTbohrcoor->setVisible(false);
        setdrawatomcoords(false);
        if (!getdrawatomindices() && !getdrawatomsymbols()){
            CHKactiveonly->setVisible(false);
            LBLselect->setVisible(false);
            BTNselectall->setVisible(false);
            BTNselectnone->setVisible(false);
        }
    }
    else{
        setdrawatomcoords(true);
        LBLcoordprecision->setVisible(true);
        SPBcoordprecision->setVisible(true);
        RBTangstromcoor->setVisible(true);
        RBTbohrcoor->setVisible(true);
        CHKactiveonly->setVisible(true);
        if (CHKactiveonly->isChecked()){
            LBLselect->setVisible(true);
            BTNselectall->setVisible(true);
            BTNselectnone->setVisible(true);
        }
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    if (visible)
        emit updatedisplay();
    return;
}

void molecule::CHKshowindices_changed(int a){
    if (a == 0){
        setdrawatomindices(false);
        if (!getdrawatomcoords() && !getdrawatomsymbols()){
            CHKactiveonly->setVisible(false);
            LBLselect->setVisible(false);
            BTNselectall->setVisible(false);
            BTNselectnone->setVisible(false);
        }
    }
    else{
        setdrawatomindices(true);
        CHKactiveonly->setVisible(true);
        if (CHKactiveonly->isChecked()){
            LBLselect->setVisible(true);
            BTNselectall->setVisible(true);
            BTNselectnone->setVisible(true);
        }
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    if (visible)
        emit updatedisplay();
    return;
}

void molecule::CHKshowsymbols_changed(int a){
    if (a == 0){
        setdrawatomsymbols(false);
        if (!getdrawatomcoords() && !getdrawatomindices()){
            CHKactiveonly->setVisible(false);
            LBLselect->setVisible(false);
            BTNselectall->setVisible(false);
            BTNselectnone->setVisible(false);
        }
    }
    else{
        setdrawatomsymbols(true);
        CHKactiveonly->setVisible(true);
        if (CHKactiveonly->isChecked()){
            LBLselect->setVisible(true);
            BTNselectall->setVisible(true);
            BTNselectnone->setVisible(true);
        }
    }
    QDLeditMolecule->update();
    QDLeditMolecule->adjustSize();
    if (visible)
        emit updatedisplay();
}

void molecule::closeisosurfeditors(){
    if (grids){
        for (int i = 0 ; i < grids->count() ; i++){
            for (int j = 0 ; j < grids->at(i)->surfaces->count() ; j++){
                grids->at(i)->surfaces->at(j)->closeeditor();
            }
        }
    }
}

void molecule::deletegrid(int i){
    closeisosurfeditors();
    if (grids->count() > i){
        delete grids->at(i);
        grids->removeAt(i);
        if (QDLeditMolecule)
            updateeditMoleculeDialog();
    }
}

void molecule::deletesurf(int i){
    if (surfaces->count() > i){
        closeisosurfeditors();
        delete surfaces->at(i);
        surfaces->removeAt(i);
        if (QDLeditMolecule)
            updateeditMoleculeDialog();
    }
}

//  Function loadgrid: loads a grid for isosurfaces plotting
//
bool molecule::loadgrid(){
    LBLloadinggrid->setVisible(true);
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString filename = filedialog.getOpenFileName(nullpointer,tr("Open 3D grid file ..."),path, tr("3D grid files")
                + " *.plt" + " (*.plt);;"+tr("All files")+" (*)");
    if (filename.length()==0){
        return false;
    }
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf=="plt"){
        if (grids->last()->readpltnew(filename)){
            grids->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            grids->last()->setname(QFileInfo(filename).completeBaseName());
            return true;
        }
        else
            return false;
    }else{
        QMessageBox msgBox;
        msgBox.setText(tr("Loadgrid"));
        msgBox.setInformativeText(tr("Extension of file %1 must be plt").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
}

//  Function loadsurf: loads a surface for 3D plotting
//
bool molecule::loadsurf(){
    QFileDialog filedialog(this);
    filedialog.setDirectory(ProjectFolder);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    QString filename = filedialog.getOpenFileName(this,tr("Open a surface file"),path, tr("Allowed files")
            + " (*.srf *.sgh *.basins *.isoden *.isopot);;  *.srf *.sgh *.basins *.isoden *.isopot;; "
            + tr("All files")+" (*)");
    if (filename.length()==0){
        return false;
    }
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf=="srf"){
        if (surfaces->last()->readsrfnew(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else
            return false;
    }
    else if (suf=="sgh"){
        if (surfaces->last()->readsgh(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else
            return false;
    }
    else if (suf=="isoden" || suf=="isopot"){
        if (surfaces->last()->readisosurfnew(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else
            return false;
    }
    else if (suf=="basins"){
        if (surfaces->last()->readbasinsnew(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else{
            surfaces->removeLast();
            return false;
        }
    }
    else{
        QMessageBox msgBox;
        msgBox.setText(tr("loadsurf"));
        msgBox.setInformativeText(tr("Extension of file %1 must be srf, sgh, isoden or isopot").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
}

void molecule::QDLeditMolecule_close(){
    if (QDLeditMolecule){
        for (int i = 0 ; i < grids->count() ; i++){
            for (int j = 0 ; j < grids->at(i)->surfaces->count() ; j++){
                grids->at(i)->surfaces->at(j)->closeeditor();
            }
        }
        if (scrollArea){
            scrollAreaposition = scrollArea->pos();
            scrollArea->hide();
        }
    }
    emitupdateRightMenu();
}

bool molecule::QDLeditMolecule_isVisible(){
    if (scrollArea && !scrollArea->isHidden())
        return true;
    else
        return false;
}

double molecule::getz_trans_ini(){
//    return Z_TRANS_INI;
    return 0;
}

void molecule::QDLeditMolecule_raise(){
    if (scrollArea)
        scrollArea->raise();
}

void molecule::QDLeditMolecule_delete(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    closeisosurfeditors();
    delete LBLloadinggrid;
    LBLloadinggrid = nullpointer;
    delete BTNfont;
    BTNfont = nullpointer;
    delete BTNfontcolor;
    BTNfontcolor = nullpointer;
    delete BTNselectall;
    BTNselectall = nullpointer;
    delete BTNselectnone;
    BTNselectnone = nullpointer;
    delete BTNhide;
    BTNhide = nullpointer;
    delete CHKactiveonly;
    CHKactiveonly = nullpointer;
    delete SPBrot_angle;
    SPBrot_angle = nullpointer;
    delete SPBrot_x;
    SPBrot_x = nullpointer;
    delete SPBrot_y;
    SPBrot_y = nullpointer;
    delete SPBrot_z;
    SPBrot_z = nullpointer;
    delete SPBtras_x;
    SPBtras_x = nullpointer;
    delete SPBtras_y;
    SPBtras_y = nullpointer;
    delete SPBtras_z;
    SPBtras_z = nullpointer;
    delete BTNanimation;
    BTNanimation = nullpointer;
    delete CHKrotatex;
    CHKrotatex = nullpointer;
    delete CHKrotatey;
    CHKrotatey = nullpointer;
    delete CHKrotatez;
    CHKrotatez = nullpointer;
    delete SLDspeed;
    SLDspeed = nullpointer;
    delete QDLeditMolecule;
    QDLeditMolecule = nullpointer;
    return;
}


//	Imports a file with critical points
void molecule::readcpsfiles_dialog()
{
    bool update = false;
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file with critical points"),path, tr("Allowed files")
            + " *cps-?.xyz " + " (*cps-?.xyz);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTcps->setText(QFileInfo(fileName).fileName());
    if (cps->readcpsfile(fileName)){
        cps->createCPs();
        update = true;
    }
    else {
        return;
    }
    if (cps->readcpseigen(fileName)){
        cps->createCPseigen();
    }
    if (cps->geterroreigvec())
        FRMcpeigvec->setVisible(false);
    else {
        update = true;
    }
    if (update){
        emitupdatedisplay();
    }
}

//	Imports a file with field lines
void molecule::readfieldlines_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file with field lines"),path, tr("Allowed files")
            + " (*.cam *.dengr *.gpdat);; *.cam;; *.dengr;; *gpdat;;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    ProjectFolder = QFileInfo(fileName).path();
    TXTfieldlines->setText(QFileInfo(fileName).fileName());
    if (flines->readfieldlines(fileName)){
        if (CHKflinesarrows)
            CHKflinesarrows->setChecked(false);
    }
    emitupdatedisplay();
}

//	Imports a file with HF forces
void molecule::readforcefiles_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file with HF forces"),path, tr("Allowed files")
            + " *.forces " + " (*.forces);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    TXTforces->setText(fileName.trimmed());
    if (hfforces->readforcesfile(TXTforces->text()))
        hfforces->createforces();
}

void molecule::resetinterval(){
    interval = MAX_INTERVAL - dltinterval * SLDspeed->value();
    deltaAngles = 1.f;
    if (SLDspeed->value() > INTERVAL_SCALE/4)
        deltaAngles = 2.f;
    if (SLDspeed->value() > INTERVAL_SCALE/2)
        deltaAngles = 4.f;
    if (SLDspeed->value() > 3*INTERVAL_SCALE/4)
        deltaAngles = 8.f;
    if (startanimation){
        timer->stop();
        timer->start(interval);
    }
}

void molecule::resizeQDLeditMolecule(){
    if (QDLeditMolecule){
        QDLeditMolecule->update();
        QDLeditMolecule->adjustSize();
        if (scrollArea){
            scrollArea->updatesize(QDLeditMolecule->size());
        }
    }
}

//  Function retrievegrid: retrieves a surface for 3D plotting
//
bool molecule::retrievegrid(QString filename){
    grids->append(new grid());
    grids->last()->set_ProjectFolder(ProjectFolder);
    grids->last()->set_ProjectName(ProjectName);
    if (QDLeditMolecule && QDLeditMolecule->isVisible()){
        grids->last()->setinitialposition(QPoint(QDLeditMolecule->x()+QDLeditMolecule->width(), QDLeditMolecule->y()));
    }
    else{
        grids->last()->setinitialposition(QPoint(400,10));
    }
    if (grids->last()->readpltnew(filename)){
        grids->last()->setfullname(QFileInfo(filename).absoluteFilePath());
        grids->last()->setname(QFileInfo(filename).fileName());
        if (QDLeditMolecule)
            updateeditMoleculeDialog();
        connect (grids->last(), SIGNAL(surfaceadded()), this, SLOT(updateeditMoleculeDialog())) ;
        connect (grids->last(), SIGNAL(surfacedeleted()), this, SLOT(updateeditMoleculeDialog())) ;
        return true;
    }
    else{
        grids->removeLast();
        return false;
    }
}

//  Function retrievesurf: retrieves a surface for 3D plotting
//
bool molecule::retrievesurf(QString filename){
    surfaces->append(new surface());
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf=="srf"){
        if (surfaces->last()->readsrfnew(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else{
            surfaces->removeLast();
            return false;
        }
    }
    else if (suf=="sgh"){
        if (surfaces->last()->readsgh(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else{
            surfaces->removeLast();
            return false;
        }
    }
    else if (suf=="basins"){
        if (surfaces->last()->readbasinsnew(filename)){
            surfaces->last()->setfullname(QFileInfo(filename).absoluteFilePath());
            surfaces->last()->setname(QFileInfo(filename).fileName());
            return true;
        }
        else{
            surfaces->removeLast();
            return false;
        }
    }
    else{
        surfaces->removeLast();
        return false;
    }
}

void molecule::reset_rotation(){
    SPBrot_x->setValue(0);
    SPBrot_y->setValue(0);
    SPBrot_z->setValue(0);
    SPBrot_angle->setValue(0);
    rotationAxis = QVector3D(0,0,0);
    rotation_changed();
}


void molecule::reset_translation(){
    SPBtras_x->setValue(0);
    SPBtras_y->setValue(0);
    SPBtras_z->setValue(0);
    translation_changed();
}

void molecule::rotation_changed(){
    rotation = QQuaternion::fromAxisAndAngle(SPBrot_x->value(),SPBrot_y->value(),SPBrot_z->value(),SPBrot_angle->value());
    emit updatedisplay();
}

bool molecule::startanimate(){
    if (!rotatex && !rotatey && !rotatez){
        startanimation = false;
        timer->stop();
        if (BTNanimation)
            BTNanimation->setText(tr("Start"));
        return false;
    }
    else{
        startanimation = true;
        resetinterval();
        timer->start();
        if (BTNanimation)
            BTNanimation->setText(tr("Stop"));
        return true;
    }
}

bool molecule::stopanimate(){
    if (!CHKrotatex || !CHKrotatey || !CHKrotatez){
        return false;
    }
    else{
        startanimation = false;
        timer->stop();
        BTNanimation->setText(tr("Start"));
        return true;
    }
}

void molecule::toggleshowsurf(int i){
    surfaces->at(i)->toggleshowsurf();
}

void molecule::translation_changed(){
    if (angstrom){
        translation = QVector3D(SPBtras_x->value()*BOHR_TO_ANGSTROM,SPBtras_y->value()*BOHR_TO_ANGSTROM,
                                SPBtras_z->value()*BOHR_TO_ANGSTROM);
    }
    else{
        translation = QVector3D(SPBtras_x->value(),SPBtras_y->value(),SPBtras_z->value());
    }
    emit updatedisplay();
}

void molecule::SPBstepwheel_changed(){
    setstepwheel(SPBstepwheel->value());
    stepwheel = SPBstepwheel->value();
    emit updatedisplay();
}

void molecule::TXTcps_changed(){
//    if (!TXTcps->text().isEmpty())
        cps->readcpsfile(TXTcps->text());
}

void molecule::TXTfieldlines_changed(){
//    if (!TXTfieldlines->text().isEmpty())
        flines->readfieldlines(TXTfieldlines->text());
}

void molecule::TXTforces_changed(){
//    if (!TXTforces->text().isEmpty())
        hfforces->readforcesfile(TXTforces->text());
}

void molecule::updateeditMoleculeDialog(){
    if (QDLeditMolecule){
        delete BTNhide;
        BTNhide = nullpointer;
        QDLeditMolecule_delete();
    }
    if (scrollArea){
        scrollAreaposition = scrollArea->pos();
        delete scrollArea;
        scrollArea = nullpointer;
    }
    createeditMoleculeDialog();
    emit updatedisplay();
}

void molecule::updateQDLeditMolecule(bool a){
    if (QDLeditMolecule){
        if (!a){
            for (int i = 0 ; i < grids->count() ; i++){
                for (int j = 0 ; j < grids->at(i)->surfaces->count() ; j++){
                    grids->at(i)->surfaces->at(j)->closeeditor();
                }
            }
            updateeditMoleculeDialog();
        }
    }
    resizeQDLeditMolecule();
}

// Functions for building structures and related functions
// =======================================================

void molecule::makeStructureBondsandSticks(){
    QVector <int> znuc;
    QVector <QVector3D> xyz;
    znuc =  getcharges();
    xyz = getxyz();
    allindices.clear();
    allvertices.clear();
    allindicesoffset.clear();
    allindicesoffset.append(0);
//    Make Balls
    if (!hideatoms){
        for (int i = 0 ; i < xyz.length() ; i++){
            if (hidehydrogens && znuc.at(i) == 1)
                continue;
            float scale;
            if (scaleradii){
                scale = elem->getAtomicRadius(znuc[i]);
                if (scale <= -999)
                    scale = 1.;
            }
            else
                scale = 1.;

            makeSphere(30,30,ballradius * scale,xyz[i],elem->getjmolColor(znuc[i]));
        }
    }
//    Make Sticks
    if (hidebonds)
        return;
    for (int i = 1 ; i < xyz.length() ; i++){
        if (hidehydrogens && znuc.at(i) == 1)
            continue;
        for (int j = 0 ; j < i ; j++){
            if (hidehydrogens && znuc.at(j) == 1)
                continue;
            QVector3D Rdiff = xyz[j]-xyz[i];
            qreal Rcov = elem->getCovalentRadius(znuc[i]) + elem->getCovalentRadius(znuc[j]);
            qreal Rdiffmod = QVector3D::dotProduct(Rdiff,Rdiff);
            if ( Rdiffmod < disthressq * Rcov * Rcov){
                makeCylinder(30, 20, cylradius, 0.5*sqrt(Rdiffmod), xyz[i], Rdiff, elem->getjmolColor(znuc[i]));
                makeCylinder(30, 20, cylradius, 0.5*sqrt(Rdiffmod), xyz[j], -Rdiff, elem->getjmolColor(znuc[j]));
            }
        }
    }
}


// Function makeSphere: generates vertices for drawing a sphere with a given number of "slices" (longitude)
//      and "stacks" (latitude) and a "radius"
//      optionally the "center" can be supplied (default: (0,0,0))
//
// Counterclockwise triangles generated
//
void molecule::makeSphere(int slices, int stacks, double radius){
    makeSphere(slices, stacks, radius, QVector3D(0,0,0), QColor(1,0,0));
}

void molecule::makeSphere(int slices, int stacks, double radius, QVector3D center){
//    int stacks = 40, slices = 40; // stacks = no. of Latitude lines,
//                                  // slices = no. of Longitude lines
    makeSphere(slices, stacks, radius, center, QColor(1,0,0));
}

void molecule::makeSphere(int slices, int stacks, double radius, QColor color){
//    int stacks = 40, slices = 40; // stacks = no. of Latitude lines,
//                                  // slices = no. of Longitude lines
    makeSphere(slices, stacks, radius, QVector3D(0,0,0), color);
}

void molecule::makeSphere(int slices, int stacks, double radius, QVector3D center, QColor color){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines
    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = PI / stacks;

    int numIndices;

    QVector4D vcolor = QVector4D(color.redF(), color.greenF(), color.blueF(), 1.);

    // Generate vertices coordinates, normal values, and texture coordinates

    int numVertices = slices * (stacks - 1) + 2;
    QVector3D *vertices = new QVector3D[numVertices];
    QVector3D *normals = new QVector3D [numVertices];
    QVector2D *textures = new QVector2D [numVertices];

    // North pole point
    normals[0] = QVector3D(0.,0.,1.);
    vertices[0] = QVector3D(0.,0.,radius) + center;
    textures[0] = QVector2D(0.5f,1.0f);

    // vertices on the main body
    int k = 1;
    for (int i = 1; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
            normals[k] = QVector3D(sin(deltaLat * i) * cos(deltaLong * j), sin(deltaLat * i) * sin(deltaLong * j),
                                cos(deltaLat * i));
            vertices[k] = radius * normals[k] + center;
            textures[k] = QVector2D((float) j / slices, 1 - (float) i / stacks);
            k++;
        }
    }

    // South pole point
    normals[k] = QVector3D(0,0,-1);
    vertices[k] = radius * normals[k] + center;
    textures[k] = QVector2D(0.5f,0.0f);

    // Generates the indices
    numIndices = (stacks - 1) * slices * 6; //why multiply by 6?
    GLuint *indices = new GLuint[numIndices];

    //add indices in North Pole region (no. of elements is slices * 3)
    k = 0;
    for (int j = 1; j<= slices-1; j++){
        indices[k++] = 0;
        indices[k++] = j;
        indices[k++] = j+1;
    }
    indices[k++] = 0;
    indices[k++] = slices;
    indices[k++] = 1;

    //add indices in South Pole Region (no. of element is slices * 3)
    int temp = numVertices  - 1;
    for (int j  = temp-1; j > temp - slices; j--){
        indices [k++] = temp;
        indices [k++] = j;
        indices [k++] = j - 1;
    }
    indices [k++] = temp;
    indices [k++] = temp-slices;
    indices [k++] = temp-1;

    // add body (no. of element is (stacks - 2) * slices * 6
    for (int i = 1; i < stacks - 1; i++) {
        for (int j = 1; j <= slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            indices[k++] = (i - 1) * slices + j;
            indices[k++] = i * slices + j;
            indices[k++] = i * slices + j + 1;
            // triangle two
            indices[k++] = (i - 1) * slices + j;
            indices[k++] = i * slices + j + 1;
            indices[k++] = (i - 1) * slices + j + 1;
        }
        // triangle one
        indices[k++] = i * slices;
        indices[k++] = (i+1) * slices;
        indices[k++] = i * slices + 1;
        // triangle two
        indices[k++] = i * slices;
        indices[k++] = i * slices + 1;
        indices[k++] = (i - 1) * slices + 1;
    }

    int j = allvertices.length();
    for (int i = 0 ; i < numIndices ; i++){
        allindices.append(indices[i]+j);
    }
    allindicesoffset.append(allindicesoffset.last()+numIndices);    // Prepares index offset for next object

    VertexNormalData Vertices;
    for (int i = 0 ; i < numVertices; i++){
        Vertices.position = vertices[i];
        Vertices.normal = normals[i];
        Vertices.color = vcolor;
        allvertices.append(Vertices);
    }
}
//  End of function makeSphere
//  ---------------------------------------------------------------------------------------------------------------------------

// Function makeCylinder: generates vertices for drawing the lateral surface of a cylinder with a given number of "slices" (longitude)
//      and "stacks" (latitude), a "radius" and a "height"
//      optionally the "center" of the base can be supplied (default: (0,0,0)) as well as the "axis" orientation (default (0,0,1))
//
// Counterclockwise triangles generated
//
void molecule::makeCylinder(int slices, int stacks, double radius, double height){
    makeCylinder(slices, stacks, radius, height, QVector3D(0,0,0), QVector3D(0,0,1), QColor(0,0,1));
}

void molecule::makeCylinder(int slices, int stacks, double radius, double height, QVector3D center){
    makeCylinder(slices, stacks, radius, height, center, QVector3D(0,0,1), QColor(0,0,1));
}

void molecule::makeCylinder(int slices, int stacks, double radius, double height, QVector3D center, QVector3D axis){
    makeCylinder(slices, stacks, radius, height, center, axis, QColor(0,0,1));
}

void molecule::makeCylinder(int slices, int stacks, double radius, double height, QVector3D center,
    QVector3D axis, QColor color){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines
    // radius = base radius
    // height = cylinder height
    // center = position of base center
    // axis = direction of cylinder axis
    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = height / (stacks-1);

    QVector3D vertex;
    QVector <QVector3D> cylindernormals;
    QQuaternion q;
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),axis);
    int k = allvertices.length();   // Starting value of indices

//    Generate normals
    for (int j = 0; j < slices; j++) {
        cylindernormals << q.rotatedVector(QVector3D(cos(deltaLong * j), sin(deltaLong * j), 0));
    }

    VertexNormalData Vertices;
    QVector4D vcolor = QVector4D(color.redF(), color.greenF(), color.blueF(), 1.);
//  Generates and stores vertices, normals and colors on the main body
    for (int i = 0; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
            vertex = q.rotatedVector(QVector3D(radius * cos(deltaLong * j), radius * sin(deltaLong * j), deltaLat * i))+center;
            Vertices.position = vertex;
            Vertices.normal = cylindernormals.at(j);
            Vertices.color = vcolor;
            allvertices.append(Vertices);
        }
    }
// Generates and stores the indices
// (no. of element is (stacks - 2) * slices * 6

    for (int i = 1; i < stacks ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            allindices.append((i - 1) * slices + j + k);
            allindices.append(i * slices + j + 1 + k);
            allindices.append(i * slices + j + k);

            // triangle two
            allindices.append((i - 1) * slices + j + k);
            allindices.append((i - 1) * slices + j + 1 + k);
            allindices.append(i * slices + j + 1 + k);


        }
        // triangle one
        allindices.append(i * slices - 1 + k);
        allindices.append(i * slices + k);
        allindices.append((i+1) * slices - 1 + k);

        // triangle two
        allindices.append(i * slices - 1 + k);
        allindices.append((i - 1) * slices + k);
        allindices.append(i * slices + k);
    }
    allindicesoffset.append(allindicesoffset.last()+(stacks - 1) * slices * 6);  // Prepares index offset for next object
}

//  End of function makeCylinder
//  ----------------------------------------------------------------------------------------------------------------------------


//  ------------------------------------------------------------------------------------------------------------------
//
//      General purpose functions
//
//  ------------------------------------------------------------------------------------------------------------------

bool molecule::isactive(){
    return active;
}

bool molecule::isaxeslabelsvisible(){
    return axeslabels_visible;
}

bool molecule::isaxesvisible(){
    return axes_visible;
}

bool molecule::isvisible(){
    return visible;
}

bool molecule::isatomactive(int i){
    return atomactive[i];
}

bool molecule::getangstromcoor(){
    return angstromcoor;
}

bool molecule::getangstromcp(){
    return angstromcp;
}

bool molecule::getdrawatomcoords(){
    return drawatomcoords;
}

bool molecule::getdrawatomindices(){
    return drawatomindices;
}

bool molecule::getdrawatomsymbols(){
    return drawatomsymbols;
}

bool molecule::gethideatoms(){
    return hideatoms;
}

bool molecule::gethidebonds(){
    return hidebonds;
}

bool molecule::gethidehydrogens(){
    return hidehydrogens;
}

bool molecule::getiscluster(){
    return iscluster;
}

int molecule::getcoordprecision(){
    return coordprecision;
}

int molecule::getlabelsvshift(){
    return labelsvshift;
}

int molecule::getnumcps(int i){
    if (cps)
        return cps->cpsxyzval[i].count();
    else {
        return 0;
    }
}

int molecule::getznuc(int i){
    return znuc.at(i);
}

int molecule::getnumatoms(){
    return znuc.count();
}

QColor molecule::getfontcolor(){
    return fontcolor;
}

QColor molecule::getXaxis_color(){
    return Xaxis_color;
}

QColor molecule::getYaxis_color(){
    return Yaxis_color;
}

QColor molecule::getZaxis_color(){
    return Zaxis_color;
}

QFont molecule::getfont(){
    return font;
}

QFont molecule::getfontaxeslabels(){
    return fontaxeslabels;
}

QPoint molecule::getparentposition(){
    return parentposition;
}

bool molecule::getonlyatomactive(){
    return onlyatomactive;
}

QQuaternion molecule::getrotation(){
    return rotation;
}

QVector3D molecule::getrotationAxis(){
    return rotationAxis;
}

QString molecule::getfullname(){
    return fullname;
}

QString molecule::getname(){
    return name;
}

QString molecule::getpath(){
    return path;
}

QVector<QVector3D> molecule::getpositionaxeslabels(){
    return positionaxeslabels;
}

QVector3D molecule::gettranslation(){
    return translation;
}

void molecule::darken(){
    VertexNormalData Vertices;
    for (int i = 0 ; i < allvertices.length() ; i++){
        Vertices = allvertices.at(i);
        Vertices.color -=  darkenshift;
        allvertices.replace(i,Vertices);
    }
    if (surfaces && !surfaces->isEmpty()){
        for (int i = 0 ; i < surfaces->length() ; i++){
            for (int j = 0 ; j < surfaces->at(i)->allvertices.length() ; j++){
                Vertices = surfaces->at(i)->allvertices.at(j);
                Vertices.color -=  darkenshift;
                surfaces->at(i)->allvertices.replace(j,Vertices);
            }
        }
    }
    if (grids && !grids->isEmpty()){
        for (int i = 0 ; i < grids->count() ; i++){
            if (grids->at(i)->surfaces && !grids->at(i)->surfaces->isEmpty()){
                for (int j = 0 ; j < grids->at(i)->surfaces->count() ; j++){
                    for (int k = 0 ; k < grids->at(i)->surfaces->at(j)->allvertices.length() ; k++){
                        Vertices = grids->at(i)->surfaces->at(j)->allvertices.at(k);
                        Vertices.color -=  darkenshift;
                        grids->at(i)->surfaces->at(j)->allvertices.replace(k,Vertices);
                    }
                }
            }
        }
    }
    emit updatedisplay();
}

void molecule::emitupdatedisplay(){
    emit updatedisplay();
}

void molecule::emitupdateGL(){
    emit updateGL();
}

void molecule::emitupdateRightMenu(){
    emit updateRightMenu();
}

void molecule::initatomactive(bool a){
    for (int i = 0 ; i < znuc.count() ; i++)
        atomactive << a;
}

void molecule::lighten(){
    VertexNormalData Vertices;
    for (int i = 0 ; i < allvertices.length() ; i++){
        Vertices = allvertices.at(i);
        Vertices.color +=  darkenshift;
        allvertices.replace(i,Vertices);
    }
    if (surfaces && !surfaces->isEmpty()){
        for (int i = 0 ; i < surfaces->length() ; i++){
            for (int j = 0 ; j < surfaces->at(i)->allvertices.length() ; j++){
                Vertices = surfaces->at(i)->allvertices.at(j);
                Vertices.color +=  darkenshift;
                surfaces->at(i)->allvertices.replace(j,Vertices);
            }
        }
    }
    if (grids && !grids->isEmpty()){
        for (int i = 0 ; i < grids->count() ; i++){
            if (grids->at(i)->surfaces && !grids->at(i)->surfaces->isEmpty()){
                for (int j = 0 ; j < grids->at(i)->surfaces->count() ; j++){
                    for (int k = 0 ; k < grids->at(i)->surfaces->at(j)->allvertices.length() ; k++){
                        Vertices = grids->at(i)->surfaces->at(j)->allvertices.at(k);
                        Vertices.color +=  darkenshift;
                        grids->at(i)->surfaces->at(j)->allvertices.replace(k,Vertices);
                    }
                }
            }
        }
    }

    emit updatedisplay();
}

void molecule::loadallindices(QVector<GLuint> allind){
    allindices = allind;
}

void molecule::loadallindicesoffset(QVector<GLuint> allindoff){
    allindicesoffset = allindoff;
}

void molecule::loadallvertices(QVector<VertexNormalData> v){
    allvertices = v;
}

void molecule::loadcharges(QVector <int> q){
    znuc = q;
}

void molecule::loadxyz(QVector <QVector3D> v){
    xyz = v;
}

void molecule::make_axes(){
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    QQuaternion q;
    QVector3D scale = SCALE * QVector3D(SPBaxesthickness->value(),SPBaxesthickness->value(),MOLSCALEHEIGHT*SPBaxeslength->value());

    VertexNormalData v;
    positionaxeslabels.clear();
    allaxesindices.clear();
    allaxesvertices.clear();
    allaxesindicesoffset.clear();
    allaxesindicesoffset.append(0);
//    X axis
    uint maxindex = 0;
    int kshift = 0;
    for (int k = 0 ; k < cylinderindices.length() ; k++){
        allaxesindices << cylinderindices.at(k) + kshift;
        maxindex = std::max(maxindex,allaxesindices.last()+1);
    }
    allaxesindicesoffset.append(allaxesindices.length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(1.,0.,0.));
    color = QVector4D(Xaxis_color.redF(), Xaxis_color.greenF(), Xaxis_color.blueF(), 1.);
    for (int k = 0 ; k < cylindervertices.length() ; k++){
        position = q.rotatedVector( scale * cylindervertices.at(k));
        normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices << v;
    }

//    Y axis
    kshift = maxindex;
    for (int k = 0 ; k < cylinderindices.length() ; k++){
        allaxesindices << cylinderindices.at(k) + kshift;
        maxindex = std::max(maxindex,allaxesindices.last()+1);
    }
    allaxesindicesoffset.append(allaxesindices.length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,1.,0.));
    color = QVector4D(Yaxis_color.redF(), Yaxis_color.greenF(), Yaxis_color.blueF(), 1.);
    for (int k = 0 ; k < cylindervertices.length() ; k++){
        position = q.rotatedVector( scale * cylindervertices.at(k));
        normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices << v;
    }
//    Z axis
    kshift = maxindex;
    for (int k = 0 ; k < cylinderindices.length() ; k++){
        allaxesindices << cylinderindices.at(k) + kshift;
        maxindex = std::max(maxindex,allaxesindices.last()+1);
    }
    allaxesindicesoffset.append(allaxesindices.length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,0.,1.));
    color = QVector4D(Zaxis_color.redF(), Zaxis_color.greenF(), Zaxis_color.blueF(), 1.);
    for (int k = 0 ; k < cylindervertices.length() ; k++){
        position = q.rotatedVector( scale * cylindervertices.at(k));
        normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices << v;
    }

//      Arrows
    scale = SCALE * QVector3D(SPBaxesarrowwidth->value(),SPBaxesarrowwidth->value(),MOLSCALEARROWSHEIGHT*SPBaxesarrowsize->value());
    QVector3D shiftorig = SCALE * QVector3D(0,0,MOLSCALEHEIGHT*SPBaxeslength->value());
//    X axis arrow
    kshift = maxindex;
    for (int k = 0 ; k < coneindices.length() ; k++){
        allaxesindices << coneindices.at(k) + kshift;
        maxindex = std::max(maxindex,allaxesindices.last()+1);
    }
    allaxesindicesoffset.append(allaxesindices.length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(1.,0.,0.));
    color = QVector4D(Xaxis_color.redF(), Xaxis_color.greenF(), Xaxis_color.blueF(), 1.);
    for (int k = 0 ; k < conevertices.length() ; k++){
        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig);
        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices << v;
    }
    QFontMetrics fm(getfontaxeslabels());
    QSize fmsize = fm.size( Qt::TextSingleLine, "x" );
    int shift = (int)std::sqrt((double)(fmsize.width()*fmsize.width())+(double)(fmsize.height()*fmsize.height()));
    positionaxeslabels << q.rotatedVector( QVector3D(0.,0.,
                    SCALE * (shift + 2. * MOLSCALEARROWSHEIGHT*SPBaxesarrowsize->value()) ) + shiftorig);
//    Y axis arrow
    kshift = maxindex;
    for (int k = 0 ; k < coneindices.length() ; k++){
        allaxesindices << coneindices.at(k) + kshift;
        maxindex = std::max(maxindex,allaxesindices.last()+1);
    }
    allaxesindicesoffset.append(allaxesindices.length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,1.,0.));
    color = QVector4D(Yaxis_color.redF(), Yaxis_color.greenF(), Yaxis_color.blueF(), 1.);
    for (int k = 0 ; k < conevertices.length() ; k++){
        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig);
        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices << v;
    }
    fmsize = fm.size( Qt::TextSingleLine, "x" );
    shift = (int)std::sqrt((double)(fmsize.width()*fmsize.width())+(double)(fmsize.height()*fmsize.height()));
    positionaxeslabels << q.rotatedVector( QVector3D(0.,0.,
                    SCALE * (shift + 2. * MOLSCALEARROWSHEIGHT*SPBaxesarrowsize->value()) ) + shiftorig);
//    Z axis arrow
    kshift = maxindex;
    for (int k = 0 ; k < coneindices.length() ; k++){
        allaxesindices << coneindices.at(k) + kshift;
        maxindex = std::max(maxindex,allaxesindices.last()+1);
    }
    allaxesindicesoffset.append(allaxesindices.length());
    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(0.,0.,1.));
    color = QVector4D(Zaxis_color.redF(), Zaxis_color.greenF(), Zaxis_color.blueF(), 1.);
    for (int k = 0 ; k < conevertices.length() ; k++){
        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig);
        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
        v.position.setX(position.x());
        v.position.setY(position.y());
        v.position.setZ(position.z());
        v.normal.setX(normal.x());
        v.normal.setY(normal.y());
        v.normal.setZ(normal.z());
        v.color = color;
        allaxesvertices << v;
    }
    fmsize = fm.size( Qt::TextSingleLine, "x" );
    shift = (int)std::sqrt((double)(fmsize.width()*fmsize.width())+(double)(fmsize.height()*fmsize.height()));
    positionaxeslabels << q.rotatedVector( QVector3D(0.,0.,
                    SCALE * (shift + 2. * MOLSCALEARROWSHEIGHT*SPBaxesarrowsize->value()) ) + shiftorig);
}

// Function makeAxesCone: generates vertices for drawing a cone with a given number of "slices" (longitude)
//      and "stacks" (latitude), a radius equal to 1 and a "height"
//      Base centered at (0,0,0) top at (0,0,height)
//
//  Counterclockwise triangles generated
//
void molecule::makeAxesCone(int slices, int stacks, float height){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines

    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = height / (stacks-1);

    conevertices.clear();
    conenormals.clear();

    // vertices on the main body
    float aux = 1. / stacks;
    float bux = 1 / height;
    for (int i = 0; i < stacks-1; i++) {
        float rad = 1. - i * aux;
        float h = deltaLat * i;
        for (int j = 0; j < slices; j++) {
            conevertices << QVector3D(rad*cos(deltaLong * j), rad*sin(deltaLong * j), h);
            conenormals << QVector3D(cos(deltaLong * j), sin(deltaLong * j), bux).normalized();
        }
    }
    conevertices << QVector3D(0,0, height);
    conenormals << QVector3D(0,0,1.);

    // vertices on the base
    int indexbase = conevertices.length();
    conevertices << QVector3D(0,0,0);
    conenormals << QVector3D(0,0,-1.);
    for (int j = 0; j < slices; j++) {
        conevertices << QVector3D(cos(deltaLong * j), sin(deltaLong * j), 0);
        conenormals << QVector3D(0,0,-1.);
    }

    // Generates the indices
    coneindices.clear();
    // indices of the body
    for (int i = 1; i < stacks-1 ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            coneindices << (i - 1) * slices + j;
            coneindices << (i - 1) * slices + j + 1;
            coneindices << i * slices + j;
            // triangle two
            coneindices << (i - 1) * slices + j + 1;
            coneindices << i * slices + j + 1;
            coneindices << i * slices + j;

        }
        // triangle one
        coneindices << i * slices-1;
        coneindices << (i-1) * slices;
        coneindices << (i+1) * slices - 1;
        // triangle two
        coneindices << (i-1) * slices;
        coneindices << i * slices;
        coneindices << (i+1) * slices - 1;
    }
//    triangles of cone vertex
    for (int j = 0; j < slices-1; j++) {
        coneindices << (stacks - 2) * slices + j;
        coneindices << (stacks - 2) * slices + j + 1;
        coneindices << (stacks - 1) * slices;
    }
    coneindices << (stacks - 1) * slices - 1;
    coneindices << (stacks - 2) * slices;
    coneindices << (stacks - 1) * slices;
    // indices of the base
    for (int j = 0; j < slices-1; j++) {
        coneindices << indexbase;
        coneindices << j + 2 + indexbase;
        coneindices << j + 1 + indexbase;
    }
    coneindices << indexbase;
    coneindices << 1 + indexbase;
    coneindices << slices + indexbase;
}

//  End of function makeAxesCone
//  ----------------------------------------------------------------------------------------------------------------------------

// Function makeAxesCylinder: generates vertices for drawing the axes cylinder side surface (neither top nor base) with a given number of
//      "slices" (longitude) and "stacks" (latitude), a radius equal to 1 and a height equal to 1 to be reescaled afterwards
//      as appropriate. Base center at (0,0,0). Orientation: (0,0,1)
//
//  Counterclockwise triangles generated
//
void molecule::makeAxesCylinder(int slices, int stacks){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines
    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = 1. / (stacks-1);

    cylindervertices.clear();

    // vertices on the main body
    for (int i = 0; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
            cylindervertices << QVector3D(cos(deltaLong * j), sin(deltaLong * j), deltaLat * i);
        }
    }
    // Generates the indices
    cylinderindices.clear();
    // add body (no. of element is (stacks - 2) * slices * 6
    for (int i = 1; i < stacks ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            cylinderindices << (i - 1) * slices + j;
            cylinderindices << i * slices + j + 1;
            cylinderindices << i * slices + j;

            // triangle two
            cylinderindices << (i - 1) * slices + j;
            cylinderindices << (i - 1) * slices + j + 1;
            cylinderindices << i * slices + j + 1;


        }
        // triangle one
        cylinderindices << i * slices - 1;
        cylinderindices << i * slices;
        cylinderindices << (i+1) * slices - 1;

        // triangle two
        cylinderindices << i * slices - 1;
        cylinderindices << (i - 1) * slices;
        cylinderindices << i * slices;

    }
}


void molecule::RBTbohr_changed(){
    if (RBTbohr->isChecked()){
        angstrom = false;
        SPBtras_x->setValue(translation.x());
        SPBtras_y->setValue(translation.y());
        SPBtras_z->setValue(translation.z());
    }
    else{
        angstrom = true;
        SPBtras_x->setValue(translation.x()*BOHR_TO_ANGSTROM);
        SPBtras_y->setValue(translation.y()*BOHR_TO_ANGSTROM);
        SPBtras_z->setValue(translation.z()*BOHR_TO_ANGSTROM);
    }
}

void molecule::RBTbohrcoor_changed(){
    if (RBTbohrcoor->isChecked()){
        angstromcoor = false;
    }
    else{
        angstromcoor = true;
    }
    emit updatedisplay();
}

void molecule::RBTbohrcp_changed(){
    if (RBTbohrcp->isChecked()){
        angstromcp = false;
    }
    else{
        angstromcp = true;
    }
    emit updatedisplay();
}

void molecule::setactive(bool a){
    active = a;
    if (active)
        lighten();
    else
        darken();
}

void molecule::setallatomactive(bool a){
    for(int i = 0 ; i < atomactive.count() ; i++)
        atomactive[i] = a;
}

void molecule::setatomactive(int i, bool a){
    atomactive[i] = a;
}

void molecule::setballradius(qreal a){
    ballradius = a;
};

void molecule::setcylradius(qreal a){
    cylradius = a;
};

void molecule::setdisthressq(qreal a){
    disthressq = a;
}

void molecule::setdrawatomcoords(bool a){
    drawatomcoords = a;
}

void molecule::setdrawatomindices(bool a){
    drawatomindices = a;
}

void molecule::setdrawatomsymbols(bool a){
    drawatomsymbols = a;
}

void molecule::setfontcolor(QColor a){
    fontcolor = a;
}

void molecule::setfont(QFont a){
    font = a;
}

void molecule::setfontaxeslabels(QFont a){
    fontaxeslabels = a;
}

void molecule::setfullname(QString a){
    fullname = a;
}

void molecule::setiscluster(bool a){
    iscluster = a;
}

void molecule::setlabelsvshift(int i){
    labelsvshift = i;
}

void molecule::setname(QString a){
    name = a;
}

void molecule::setonlyatomactive(bool a){
    onlyatomactive = a;
}

void molecule::setparentposition(QPoint a){
    parentposition = a;
}

void molecule::setpath(QString a){
    path = a;
}

void molecule::set_ProjectFolder(QString a){
    ProjectFolder = a;
}

void molecule::set_ProjectName(QString name){
    ProjectName = name;
}

void molecule::setrotation(QQuaternion a){
    rotation = a;
}

void molecule::setrotationAxis(QVector3D a){
    if (a.length() > 1.e-7)
        rotationAxis = a;
    else
        rotationAxis = QVector3D(1,0,0);
}

void molecule::setstepwheel(float a){
    stepwheel = a;
}

void molecule::settranslationButtons(){
    if (!SPBtras_x || !SPBtras_y || !SPBtras_z )
        return;
    if (!angstrom){
        SPBtras_x->setValue(translation.x());
        SPBtras_y->setValue(translation.y());
        SPBtras_z->setValue(translation.z());
    }
    else{
        SPBtras_x->setValue(translation.x()*BOHR_TO_ANGSTROM);
        SPBtras_y->setValue(translation.y()*BOHR_TO_ANGSTROM);
        SPBtras_z->setValue(translation.z()*BOHR_TO_ANGSTROM);
    }
}

void molecule::setrotationButtons(){
    if (!SPBrot_x || !SPBrot_y || !SPBrot_z || !SPBrot_angle  )
        return;
    SPBrot_x->setValue(rotation.x());
    SPBrot_y->setValue(rotation.y());
    SPBrot_z->setValue(rotation.z());
    SPBrot_angle->setValue(360. * qAcos(rotation.scalar())/pi);
}

void molecule::setscaleradii(bool a){
    scaleradii = a;
}

void molecule::settranslation(QVector3D a){
    translation = a;
}

void molecule::setvisible(bool a){
    visible = a;
}

void molecule::setworld_rotation(QQuaternion a){
    world_rotation = a;
}

void molecule::toggleactive(){
    setactive(!isactive());
}

void molecule::SPBaxesarrowssize_changed(int i){
    axesarrowssize = i;
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
}

void molecule::SPBaxesarrowswidth_changed(int i){
    axesarrowswidth = i;
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
}

void molecule::SPBaxeslength_changed(int i){
    axeslength = i;
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
}

void molecule::SPBaxesthickness_changed(int i){
    axesthickness = i;
    if (axes_visible){
        make_axes();
        emitupdatedisplay();
    }
}

// Changes precision for coordinates display
void molecule::SPBcoordprecision_changed(int i){
    coordprecision = i;
    emitupdatedisplay();
}

// Changes ball radius for critical
void molecule::SPBcpballradius_changed(int i){
    cps->setradius(i);
    emitupdatedisplay();
}

// Changes precision for critical points coordinates display
void molecule::SPBcpcoordprecision_changed(int i){
    cps->setcoordprecision(i);
    emitupdatedisplay();
}

// Changes width of critical points eigenvector arrows
void molecule::SPBcpeigarrowwidth_changed(int i){
    cps->seteigarrowwidth(i);
    emitupdatedisplay();
}

// Changes length critical points eigenvector arrows
void molecule::SPBcpeigarrowsize_changed(int i){
    cps->seteigarrowsize(i);
    emitupdatedisplay();
}

// Changes length of critical points eigenvectors
void molecule::SPBcpeiglength_change(int i){
    cps->seteiglength(i);
    emitupdatedisplay();
}

// Changes thickness of critical points eigenvectors
void molecule::SPBcpeigthickness_changed(int i){
    cps->seteigthickness(i);
    emitupdatedisplay();
}

// Changes precision for critical points
void molecule::SPBcpprecision_changed(int i){
    cps->setcpprecision(i);
    emitupdatedisplay();
}

// Changes vertical display for critical points labels
void molecule::SPBcpvshift_changed(int i){
    cps->setcpvshift(i);
    emitupdatedisplay();
}

void molecule::SPBflineslinewidth_changed(int i){
    flines->setlineswidth(i);
    emitupdatedisplay();
}

void molecule::SPBflinesarrowssep_changed(int i){
    flines->setarrowsseparation(i);
    flines->createarrows();
    emitupdatedisplay();
}

void molecule::SPBflinesarrowssize_changed(int i){
    flines->setarrowssize(i);
    flines->createarrows();
    emitupdatedisplay();
}

void molecule::SPBflinesarrowswidth_changed(int i){
    flines->setarrowswidth(i);
    flines->createarrows();
    emitupdatedisplay();
}

// Changes length of force arrows
void molecule::SPBforcesarrowlength_changed(int i){
    hfforces->setarrowlength(i);
    hfforces->createforces();
}

// Changes width of force arrows
void molecule::SPBforcesarrowwidth_changed(int i){
    hfforces->setarrowwidth(i);
    hfforces->createforces();
}

// Changes length of force vectors
void molecule::SPBforceslength_changed(int i){
    hfforces->setforceslength(i);
    hfforces->createforces();
}

// Changes thickness of force vectors
void molecule::SPBforcesthickness_changed(int i){
    hfforces->setforcesthickness(i);
    hfforces->createforces();
}

// Changes vertical display for critical points labels
void molecule::SPBlabelsvshift_changed(int i){
    labelsvshift = i;
    emitupdatedisplay();
}

QVector <GLuint> molecule::getallindices(){
    return allindices;
}

QVector <GLuint> molecule::getallindicesoffset(){
    return allindicesoffset;
}

QVector <VertexNormalData> molecule::getallvertices(){
    return allvertices;
}

QVector <int> molecule::getcharges(){
    return znuc;
}

QVector <QVector3D> molecule::getxyz(){
    return xyz;
}

float molecule::getstepwheel(){
    return stepwheel;
}


/*******************************************************************************************************/
/********************************  Class editMoleculeDialog  implementation  *******************************/
/*******************************************************************************************************/

editMoleculeDialog::editMoleculeDialog(QWidget *parent) : QDialog(parent)
{
    this->setAttribute(Qt::WA_DeleteOnClose);
}

editMoleculeDialog::~editMoleculeDialog(){

}

void editMoleculeDialog::reject(){
    emit closed();
}
void editMoleculeDialog::closeEvent(QCloseEvent *event){
    event->ignore();
    emit closed();
    event->accept();
}
