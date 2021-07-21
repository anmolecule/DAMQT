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
//  File:   measures.cpp
//  Description: implements measures class
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//
#include "measures.h"
#include <cmath>


measures::measures(QWidget *parent) : QWidget(parent)
{
    connections.clear();
    connectionsInfo.clear();
    anglecenters = nullpointer;
    anglestext = nullpointer;
    BTNanglescolor = nullpointer;
    BTNdihedralscolor = nullpointer;
    BTNdistancescolor = nullpointer;
    BTNfontangles = nullpointer;
    BTNfontdihedrals = nullpointer;
    BTNfontdistances = nullpointer;
    BTNresetangles = nullpointer;
    BTNresetdihedrals = nullpointer;
    BTNresetdistances = nullpointer;
    CHKshowangles = nullpointer;
    CHKshowdihedrals = nullpointer;
    CHKshowdistances = nullpointer;
    CHKshowdistwin = nullpointer;
    dihedralcenters = nullpointer;
    dihedralstext = nullpointer;
    distancecenters = nullpointer;
    distancestext = nullpointer;
    elem = new Elements();
    LBLalpha = nullpointer;
    LBLdistvshift = nullpointer;
    LBLlastangles = nullpointer;
    LBLlastdihedrals = nullpointer;
    LBLlastdistances = nullpointer;
    FRMmeasureangles = nullpointer;
    FRMmeasuredihedrals = nullpointer;
    FRMmeasuredistances = nullpointer;
    FRMmeasuregeom = nullpointer;
    measuresInfo = nullpointer;
    printfilename = nullpointer;
    QDLmeasures = nullpointer;
    RBTangstrom = nullpointer;
    RBTbohr = nullpointer;
    RBTmeasureangle = nullpointer;
    RBTmeasuredihed = nullpointer;
    RBTmeasuredist = nullpointer;
    RBTmeasurenone = nullpointer;
    SPBanglesprecision = nullpointer;
    SPBdihedralsprecision = nullpointer;
    SPBdistprecision = nullpointer;
    anglescolor = QColor(Qt::yellow);
    dihedralscolor = QColor(Qt::yellow);
    dihedralplanescolor = QColor(182, 182, 182, 128);
    distancescolor = QColor(Qt::yellow);
    anglesfont = QFont("Noto Sans", 12, QFont::Bold);
    dihedralsfont = QFont("Noto Sans", 12, QFont::Bold);
    distancesfont = QFont("Noto Sans", 12, QFont::Bold);
    measuresInfopos = 0.2*QPoint(this->x(),this->y());
    anglecenters = new QVector<centerData>();
    anglesprinttext.clear();
    dihedralcenters = new QVector<centerData>();
    dihedralsprinttext.clear();
    distancecenters = new QVector<centerData>();
    distancesprinttext.clear();
    lastselectangles.clear();
    lastselectdihedrals.clear();
    lastselectdist.clear();
    mtrsf = new QVector<QMatrix4x4>();
    selectangles.clear();
    selectdihedrals.clear();
    selectdist.clear();
    systemnames.clear();
}

//	Destructor
measures::~measures(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
    if (QDLmeasures){
        delete QDLmeasures;
        QDLmeasures = nullpointer;
    }
    if (measuresInfo){
        delete measuresInfo;
        measuresInfo = nullpointer;
    }
}

bool measures::create_QDLmeasures(){
    QDLmeasures = new editMeasuresDialog();
    QDLmeasures->setWindowTitle(tr("Geometry measures"));
    connections << connect(QDLmeasures, SIGNAL(rejected()), this, SIGNAL(QDLmeasures_rejected()));
    connections << connect(QDLmeasures, SIGNAL(closed()), this, SLOT(emitupdateRightMenu()));

    FRMmeasuregeom = new QGroupBox();
    FRMmeasuregeom->setTitle(tr("Type of measure"));
    FRMmeasuregeom->setMinimumWidth(355);

    RBTmeasurenone= new QRadioButton(tr("None"));
    RBTmeasurenone->setChecked(true);
    RBTmeasurenone->setMinimumWidth(355);
    connections << connect(RBTmeasurenone, SIGNAL(toggled(bool)), this, SLOT(RBTmeasurenone_changed()));
    RBTmeasuredist= new QRadioButton(tr("Distances"));
    RBTmeasuredist->setChecked(false);
    RBTmeasuredist->setToolTip(tr("Shift+double click on center to select"));
    connections << connect(RBTmeasuredist, SIGNAL(toggled(bool)), this, SLOT(RBTmeasuredist_changed()));
    RBTmeasureangle= new QRadioButton(tr("Angles"));
    RBTmeasureangle->setChecked(false);
    RBTmeasureangle->setToolTip(tr("Shift+double click on center to select"));
    connections << connect(RBTmeasureangle, SIGNAL(toggled(bool)), this, SLOT(RBTmeasureangle_changed()));
    RBTmeasuredihed= new QRadioButton(tr("Dihedral angles"));
    RBTmeasuredihed->setChecked(false);
    RBTmeasuredihed->setToolTip(tr("Shift+double click on center to select"));
    connections << connect(RBTmeasuredihed, SIGNAL(toggled(bool)), this, SLOT(RBTmeasuredihed_changed()));

//    Distances

    FRMmeasuredistances = new QGroupBox();
    FRMmeasuredistances->setTitle(tr("Distances"));
    FRMmeasuredistances->setVisible(false);
    FRMmeasuredistances->setToolTip(tr("Select centers with shift+double click. Pairs of centers required"));
    CHKshowdistwin = new QCheckBox(tr("Show/hide distances in a window"));
    CHKshowdistwin->setChecked(false);
    connections << connect(CHKshowdistwin, SIGNAL(stateChanged(int)), this, SLOT(CHKshowdistwin_changed()));

    CHKshowdistances = new QCheckBox(tr("Show/hide distances in the viewer"));
    CHKshowdistances->setChecked(true);
    connections << connect(CHKshowdistances, SIGNAL(stateChanged(int)), this, SLOT(CHKshowdistances_changed()));
    emit show_distances(CHKshowdistances->isChecked());

    LBLdistvshift = new QLabel(tr("Vertical shift"));
    LBLdistvshift->setVisible(false);
    SPBdistvshift = new QSpinBox();
    SPBdistvshift->setMinimum(-200);
    SPBdistvshift->setMaximum(200);
    SPBdistvshift->setMaximumWidth(50);
    SPBdistvshift->setValue(0);
    SPBdistvshift->setVisible(false);
    connections << connect(SPBdistvshift,SIGNAL(valueChanged(int)),this,SLOT(SPBdistvshift_changed(int)));

    RBTbohr= new QRadioButton(tr("Bohr"));
    RBTbohr->setChecked(true);
    RBTangstrom= new QRadioButton(tr("Angstrom"));
    RBTangstrom->setChecked(false);
    connections << connect(RBTangstrom, SIGNAL(toggled(bool)), this, SLOT(RBTangstrom_changed()));

    QLabel *LBLdistprecision = new QLabel(tr("Precision"));
    SPBdistprecision = new QSpinBox();
    SPBdistprecision->setMinimum(0);
    SPBdistprecision->setMaximum(7);
    SPBdistprecision->setMaximumWidth(50);
    SPBdistprecision->setValue(2);
    connections << connect(SPBdistprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBdistprecision_changed(int)));

    BTNfontdistances = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNfontdistances, SIGNAL(clicked()), this, SLOT(BTNfontdistances_clicked()));
    BTNdistancescolor = new ColorButton();
    BTNdistancescolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNdistancescolor->setText(tr("Color"));
    BTNdistancescolor->setColor(&distancescolor);
    BTNdistancescolor->setEnabled(true);
    connections << connect(BTNdistancescolor, SIGNAL(clicked()), this, SLOT(BTNdistancescolor_clicked()));
    BTNresetdistances = new QPushButton(QIcon(":/images/centrar.png"),tr("Reset"));
    connections << connect(BTNresetdistances, SIGNAL(clicked()), this, SLOT(BTNresetdistances_clicked()));

    CHKdrawlines = new QCheckBox(tr("Show/hide lines"));
    CHKdrawlines->setChecked(true);
    connections << connect(CHKdrawlines, SIGNAL(stateChanged(int)), this, SLOT(CHKdrawlines_changed()));

    QLabel *LBLlineswidth = new QLabel(tr("Line width"));
    SPBlineswidth = new QSpinBox();
    SPBlineswidth->setMinimum(0);
    SPBlineswidth->setMaximum(7);
    SPBlineswidth->setMaximumWidth(50);
    SPBlineswidth->setValue(2);
    connections << connect(SPBlineswidth,SIGNAL(valueChanged(int)),this,SLOT(SPBlineswidth_changed(int)));

    QLabel *LBLlinestype = new QLabel(tr("Line type"));
    CMBlinestype = new QComboBox();
    CMBlinestype->addItem(tr("No Line"));
    CMBlinestype->addItem(tr("Solid"));
    CMBlinestype->addItem(tr("Dash"));
    CMBlinestype->addItem(tr("Dot"));
    CMBlinestype->addItem(tr("Dash Dot"));
    CMBlinestype->addItem(tr("Dash Dot Dot"));
    CMBlinestype->setCurrentIndex(3);
    connections << connect(CMBlinestype,SIGNAL(currentIndexChanged(int)),this,SLOT(CMBlinestype_changed(int)));

    CHKdisttranspbkg = new QCheckBox(tr("Label transparent background"));
    CHKdisttranspbkg->setChecked(false);
    connections << connect(CHKdisttranspbkg, SIGNAL(stateChanged(int)), this, SLOT(CHKdisttranspbkg_changed()));

    LBLlastdistances = new QLabel();
    LBLlastdistances->setVisible(false);

    distancestext = new QString("  ");

//    Angles

    FRMmeasureangles = new QGroupBox();
    FRMmeasureangles->setTitle(tr("Angles"));
    FRMmeasureangles->setVisible(false);
    FRMmeasureangles->setToolTip(tr("Select centers with shift+double click. Three bonding centers required for an angle"));
    CHKshowangwin = new QCheckBox(tr("Show/hide angles in a window"));
    CHKshowangwin->setChecked(false);
    connections << connect(CHKshowangwin, SIGNAL(stateChanged(int)), this, SLOT(CHKshowangwin_changed()));
    CHKshowangles = new QCheckBox(tr("Show/hide angles in the viewer"));
    CHKshowangles->setChecked(true);
    connections << connect(CHKshowangles, SIGNAL(stateChanged(int)), this, SLOT(CHKshowangles_changed()));
    emit show_angles(CHKshowangles->isChecked());

    QLabel *LBLanglesprecision = new QLabel(tr("Precision"));
    SPBanglesprecision = new QSpinBox();
    SPBanglesprecision->setMinimum(0);
    SPBanglesprecision->setMaximum(7);
    SPBanglesprecision->setMaximumWidth(50);
    SPBanglesprecision->setValue(4);
    connections << connect(SPBanglesprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBanglesprecision_changed(int)));

    BTNfontangles = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNfontangles, SIGNAL(clicked()), this, SLOT(BTNfontangles_clicked()));
    BTNanglescolor = new ColorButton();
    BTNanglescolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNanglescolor->setText(tr("Color"));
    BTNanglescolor->setColor(&anglescolor);
    BTNanglescolor->setEnabled(true);
    connections << connect(BTNanglescolor, SIGNAL(clicked()), this, SLOT(BTNanglescolor_clicked()));

    BTNresetangles = new QPushButton(QIcon(":/images/centrar.png"),tr("Reset"));
    connections << connect(BTNresetangles, SIGNAL(clicked()), this, SLOT(BTNresetangles_clicked()));

    CHKdrawarcs = new QCheckBox(tr("Show/hide arcs"));
    CHKdrawarcs->setChecked(true);
    connections << connect(CHKdrawarcs, SIGNAL(stateChanged(int)), this, SLOT(CHKdrawarcs_changed()));

    QLabel *LBLangleswidth = new QLabel(tr("Line width"));
    SPBangleswidth = new QSpinBox();
    SPBangleswidth->setMinimum(0);
    SPBangleswidth->setMaximum(7);
    SPBangleswidth->setMaximumWidth(50);
    SPBangleswidth->setValue(4);
    connections << connect(SPBangleswidth,SIGNAL(valueChanged(int)),this,SLOT(SPBangleswidth_changed(int)));

    QLabel *LBLanglestype = new QLabel(tr("Line type"));
    CMBanglestype = new QComboBox();
    CMBanglestype->addItem(tr("No Line"));
    CMBanglestype->addItem(tr("Solid"));
    CMBanglestype->addItem(tr("Dash"));
    CMBanglestype->addItem(tr("Dot"));
    CMBanglestype->addItem(tr("Dash Dot"));
    CMBanglestype->addItem(tr("Dash Dot Dot"));
    CMBanglestype->setCurrentIndex(1);
    connections << connect(CMBanglestype,SIGNAL(currentIndexChanged(int)),this,SLOT(CMBanglestype_changed(int)));

    QLabel *LBLarcradius = new QLabel(tr("Arc radius"));
    SPBarcradius = new QSpinBox();
    SPBarcradius->setMinimum(0);
    SPBarcradius->setMaximum(200);
    SPBarcradius->setMaximumWidth(50);
    SPBarcradius->setValue(40);
    connections << connect(SPBarcradius,SIGNAL(valueChanged(int)),this,SLOT(SPBarcradius_changed(int)));

    LBLlastangles = new QLabel();
    LBLlastangles->setVisible(false);

    anglestext = new QString("  ");

//    Dihedral angles

    FRMmeasuredihedrals = new QGroupBox();
    FRMmeasuredihedrals->setTitle(tr("Dihedral angles"));
    FRMmeasuredihedrals->setVisible(false);
    FRMmeasuredihedrals->setToolTip(tr("Select centers with shift+double click. Four centers required:")
            +tr(" first two centers common to both dihedral planes"));
    CHKshowdihedwin = new QCheckBox(tr("Show/hide dihedral angles in a window"));
    CHKshowdihedwin->setChecked(false);
    connections << connect(CHKshowdihedwin, SIGNAL(stateChanged(int)), this, SLOT(CHKshowdihedwin_changed()));
    CHKshowdihedrals = new QCheckBox(tr("Show/hide dihedral angles in the viewer"));
    CHKshowdihedrals->setChecked(true);
    connections << connect(CHKshowdihedrals, SIGNAL(stateChanged(int)), this, SLOT(CHKshowdihedrals_changed()));
    emit show_dihedrals(CHKshowdihedrals->isChecked());

    QLabel *LBLdihedralsprecision = new QLabel(tr("Precision"));
    SPBdihedralsprecision = new QSpinBox();
    SPBdihedralsprecision->setMinimum(0);
    SPBdihedralsprecision->setMaximum(7);
    SPBdihedralsprecision->setMaximumWidth(50);
    SPBdihedralsprecision->setValue(4);
    connections << connect(SPBdihedralsprecision,SIGNAL(valueChanged(int)),this,SLOT(SPBdihedralsprecision_changed(int)));

    BTNfontdihedrals = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connections << connect(BTNfontdihedrals, SIGNAL(clicked()), this, SLOT(BTNfontdihedrals_clicked()));
    BTNdihedralscolor = new ColorButton();
    BTNdihedralscolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNdihedralscolor->setText(tr("Color"));
    BTNdihedralscolor->setColor(&dihedralscolor);
    BTNdihedralscolor->setEnabled(true);
    connections << connect(BTNdihedralscolor, SIGNAL(clicked()), this, SLOT(BTNdihedralscolor_clicked()));
    BTNresetdihedrals = new QPushButton(QIcon(":/images/centrar.png"),tr("Reset"));
    connections << connect(BTNresetdihedrals, SIGNAL(clicked()), this, SLOT(BTNresetdihedrals_clicked()));


//        Planes color and opacity
    QGroupBox *FRMsurfcolor = new QGroupBox(tr("Planes color and opacity"));

    BTNsurfcolor = new ColorButton(QIcon(":/images/colores48.png"),tr("Color"));
    BTNsurfcolor->setColor(&dihedralplanescolor);
    connections << connect(BTNsurfcolor, SIGNAL(clicked()), this, SLOT(BTNsurfcolor_clicked()));

//        Opacity
    QLabel *LBLopacity = new QLabel();
    LBLopacity->setText(tr("Opacity:"));
    LBLalpha = new QLabel();
    LBLalpha->setText("156");
    SLDopacity = new QSlider(Qt::Horizontal);
    SLDopacity->setRange(0,255);
    SLDopacity->setSingleStep(1);
    SLDopacity->setPageStep(10);
    SLDopacity->setTickPosition(QSlider::TicksBelow);
    SLDopacity->setValue(156);
    connections << connect(SLDopacity, SIGNAL(valueChanged(int)), LBLalpha, SLOT(setNum(int)));
    connections << connect(SLDopacity, SIGNAL(sliderReleased()), this, SLOT(SLDopacity_released()));

    LBLlastdihedrals = new QLabel();
    LBLlastdihedrals->setVisible(false);

    dihedralstext = new QString("  ");

//    Print report

    BTNprint = new QPushButton(tr("Print report"));
    QString style = QString("QPushButton { color : rgb(0, 0, 0); background-color : rgb(170,255,127); }");
    BTNprint->setStyleSheet(style);
    connections << connect(BTNprint, SIGNAL(clicked()), this, SLOT(BTNprint_clicked()));
    printfilename = new QString();

//    Layouts

    QVBoxLayout *layout1=new QVBoxLayout(FRMmeasuregeom);
    layout1->addWidget(RBTmeasurenone);
    layout1->addWidget(RBTmeasuredist);
    layout1->addWidget(RBTmeasureangle);
    layout1->addWidget(RBTmeasuredihed);

    QHBoxLayout *layout1b=new QHBoxLayout();
    layout1b->addStretch();
    layout1b->addWidget(LBLdistvshift);
    layout1b->addWidget(SPBdistvshift);
    layout1b->addStretch();

    QHBoxLayout *layout2=new QHBoxLayout();
    layout2->addStretch();
    layout2->addWidget(RBTbohr);
    layout2->addWidget(RBTangstrom);
    layout2->addStretch();

    QHBoxLayout *layout3=new QHBoxLayout();
    layout3->addStretch();
    layout3->addWidget(LBLdistprecision);
    layout3->addWidget(SPBdistprecision);
    layout3->addStretch();

    QHBoxLayout *layout4=new QHBoxLayout();
    layout4->addStretch();
    layout4->addWidget(BTNfontdistances);
    layout4->addWidget(BTNdistancescolor);
    layout4->addStretch();
    layout4->setAlignment(Qt::AlignCenter);

    QGridLayout *layout5=new QGridLayout();
    layout5->addWidget(CHKdrawlines,0,0,1,2);
    layout5->addWidget(LBLlineswidth,1,0);
    layout5->addWidget(SPBlineswidth,1,1);
    layout5->addWidget(LBLlinestype,2,0);
    layout5->addWidget(CMBlinestype,2,1);
    layout5->addWidget(CHKdisttranspbkg,3,0,1,2);

    QHBoxLayout *layout6=new QHBoxLayout();
    layout6->addWidget(BTNresetdistances);
    layout6->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout7=new QHBoxLayout();
    layout7->addWidget(LBLlastdistances);

    QVBoxLayout *layout8=new QVBoxLayout(FRMmeasuredistances);
    layout8->addWidget(CHKshowdistwin);
    layout8->addWidget(CHKshowdistances);
    layout8->addLayout(layout1b);
    layout8->addLayout(layout2);
    layout8->addLayout(layout3);
    layout8->addLayout(layout4);
    layout8->addLayout(layout5);
    layout8->addLayout(layout6);
    layout8->addLayout(layout7);

    QHBoxLayout *layout9=new QHBoxLayout();
    layout9->addStretch();
    layout9->addWidget(LBLanglesprecision);
    layout9->addWidget(SPBanglesprecision);
    layout9->addStretch();

    QHBoxLayout *layout10=new QHBoxLayout();
    layout10->addStretch();
    layout10->addWidget(BTNfontangles);
    layout10->addWidget(BTNanglescolor);
    layout10->addStretch();
    layout10->setAlignment(Qt::AlignCenter);

    QGridLayout *layout11=new QGridLayout();
    layout11->addWidget(CHKdrawarcs,0,0,1,2);
    layout11->addWidget(LBLangleswidth,1,0);
    layout11->addWidget(SPBangleswidth,1,1);
    layout11->addWidget(LBLanglestype,2,0);
    layout11->addWidget(CMBanglestype,2,1);
    layout11->addWidget(LBLarcradius,3,0);
    layout11->addWidget(SPBarcradius,3,1);

    QHBoxLayout *layout12=new QHBoxLayout();
    layout12->addWidget(BTNresetangles);
    layout12->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout13=new QHBoxLayout();
    layout13->addWidget(LBLlastangles);

    QVBoxLayout *layout14=new QVBoxLayout(FRMmeasureangles);
    layout14->addWidget(CHKshowangwin);
    layout14->addWidget(CHKshowangles);
    layout14->addLayout(layout9);
    layout14->addLayout(layout10);
    layout14->addLayout(layout11);
    layout14->addLayout(layout12);
    layout14->addLayout(layout13);

    QHBoxLayout *layout15 = new QHBoxLayout();
    layout15->addStretch();
    layout15->addWidget(LBLdihedralsprecision);
    layout15->addWidget(SPBdihedralsprecision);
    layout15->addStretch();

    QHBoxLayout *layout16 = new QHBoxLayout();
    layout16->addStretch();
    layout16->addWidget(BTNfontdihedrals);
    layout16->addWidget(BTNdihedralscolor);
    layout16->addStretch();
    layout16->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout17 = new QHBoxLayout();
    layout17->addWidget(BTNresetdihedrals);
    layout17->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout18 = new QHBoxLayout();
    layout18->addStretch();
    layout18->addWidget(BTNsurfcolor);
    layout18->addStretch();

    QHBoxLayout *layout20 = new QHBoxLayout();
    layout20->addWidget(LBLopacity,0,Qt::AlignLeft);
    layout20->addWidget(LBLalpha,0,Qt::AlignRight);

    QVBoxLayout *layout21 = new QVBoxLayout();
    layout21->addLayout(layout20);
    layout21->addWidget(SLDopacity);
    layout21->addStretch();

    QVBoxLayout *layout22 = new QVBoxLayout(FRMsurfcolor);
    layout22->addLayout(layout18);
    layout22->addLayout(layout21);

    QHBoxLayout *layout23 = new QHBoxLayout();
    layout23->addWidget(LBLlastdihedrals);

    QVBoxLayout *layout24 = new QVBoxLayout(FRMmeasuredihedrals);
    layout24->addWidget(CHKshowdihedwin);
    layout24->addWidget(CHKshowdihedrals);
    layout24->addLayout(layout15);
    layout24->addLayout(layout16);
    layout24->addLayout(layout17);
    layout24->addWidget(FRMsurfcolor);
    layout24->addLayout(layout23);

    QVBoxLayout *layout=new QVBoxLayout(QDLmeasures);
    layout->addWidget(FRMmeasuregeom);
    layout->addWidget(FRMmeasuredistances);
    layout->addWidget(FRMmeasureangles);
    layout->addWidget(FRMmeasuredihedrals);
    layout->addWidget(BTNprint);
    layout->addStretch();

    QDLmeasures->show();
    return true;
}

//********************************************************************************
//          FUNCTIONS AND SLOTS
//********************************************************************************


void measures::append_systemnames(QString a){
    systemnames.append(a);
}

void measures::BTNanglescolor_clicked(){
    QColor currcolor = anglescolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        anglescolor =  newcolor;
        BTNanglescolor->setColor(&newcolor);
    }
    emit color_angles(anglescolor);
}

void measures::BTNdihedralscolor_clicked(){
    QColor currcolor = dihedralscolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        dihedralscolor =  newcolor;
        BTNdihedralscolor->setColor(&newcolor);
    }
    emit color_dihedrals(dihedralscolor);
}

void measures::BTNdistancescolor_clicked(){
    QColor currcolor = distancescolor;
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        distancescolor =  newcolor;
        BTNdistancescolor->setColor(&newcolor);
    }
    emit color_distances(distancescolor);
}

void measures::BTNfontangles_clicked(){
    anglesfont = QFontDialog::getFont(0, anglesfont);
    emit font_angles(anglesfont);
}

void measures::BTNfontdihedrals_clicked(){
    dihedralsfont = QFontDialog::getFont(0, dihedralsfont);
    emit font_dihedrals(dihedralsfont);
}


void measures::BTNfontdistances_clicked(){
    distancesfont = QFontDialog::getFont(0, distancesfont);
    emit font_distances(distancesfont);
}

void measures::BTNprint_clicked(){
    emit get_systemnames();
    if (anglecenters->isEmpty() && dihedralcenters->isEmpty() && distancecenters->isEmpty()){
        QMessageBox msgBox;
        msgBox.setText(tr("Print measurements report"));
        msgBox.setInformativeText(tr("No distances, angles or dihedrals angles chosen for report"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getSaveFileName(this,tr("Save measurements report as ..."), ProjectFolder,
            (tr("Geometry report files") + ": *.measures *.html" + " (*.measures *.html);;"+tr("All files")+" (*)"));
    if (fileName.isEmpty()){
        return;
    }
    if (QFileInfo(fileName).suffix() != "measures" && QFileInfo(fileName).suffix() != "html"){
        fileName.append("_measures.html");
    }
    QFile fileout(fileName);
    if (!fileout.isOpen()){
        fileout.open(QFile::Text | QFile::WriteOnly);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QList<int> systems;
    if (!distancecenters->isEmpty()){
        for (int i = 0 ; i < distancecenters->length() ; i++){
            if (systems.indexOf(distancecenters->at(i).molecule) == -1)
                systems << distancecenters->at(i).molecule;
        }
    }
    if (!anglecenters->isEmpty()){
        for (int i = 0 ; i < anglecenters->length() ; i++){
            if (systems.indexOf(anglecenters->at(i).molecule) == -1)
                systems << anglecenters->at(i).molecule;
        }
    }
    if (!dihedralcenters->isEmpty()){
        for (int i = 0 ; i < dihedralcenters->length() ; i++){
            if (systems.indexOf(dihedralcenters->at(i).molecule) == -1)
                systems << dihedralcenters->at(i).molecule;
        }
    }
    std::sort(systems.begin(), systems.end());
    bool printheader;
    int knt;
    for (int j = 0 ; j < systems.length() ; j++){
        outfile << QString("\n\n<H2>") + QString(" ").repeated(10)
                   + QString(tr("   System %1: %2   ")).arg(j+1).arg(systemnames.at(j))
                    << QString("</H2>\n");
    }
    printheader = true;
    knt = 0;
    for (int j = 0 ; j < systems.length() ; j++){
        for (int i = 0 ; i < distancecenters->length() ; i +=2){
            if (distancecenters->at(i).molecule != systems.at(j)) continue;
            if (printheader){
                outfile << QString("<p><H3>") + QString("-").repeated(18) + QString(tr("  Distances  "))
                           + QString("-").repeated(18) + QString("</H3></p><p>");
                printheader = false;
            }
            outfile << distancesprinttext.at(i/2) << "   ";
            if (knt++ > 4){
                knt = 0;
                outfile << QString("</p><p>");
            }
        }
    }
    printheader = true;
    knt = 0;
    for (int j = 0 ; j < systems.length() ; j++){
        for (int i = 0 ; i < anglecenters->length() ; i +=3){
            if (anglecenters->at(i).molecule != systems.at(j)) continue;
            if (printheader){
                outfile << QString("<p><H3>") + QString("-").repeated(19) + QString(tr("  Angles  "))
                           + QString("-").repeated(20) + QString("</H3></p><p>");
                printheader = false;
            }
            outfile << anglesprinttext.at(i/3) << "   ";
            if (knt++ > 4){
                knt = 0;
                outfile << QString("</p><p>");
            }
        }
    }
    printheader = true;
    knt = 0;
    for (int j = 0 ; j < systems.length() ; j++){
        for (int i = 0 ; i < dihedralcenters->length() ; i +=4){
            if (dihedralcenters->at(i).molecule != systems.at(j)) continue;
            if (printheader){
                outfile << QString("<p><H3>") + QString("-").repeated(15) + QString(tr("  Dihedral angles  "))
                           + QString("-").repeated(15) + QString("</H3></p><p>");
                printheader = false;
            }
            outfile << dihedralsprinttext.at(i/4) << "   ";
            if (knt++ > 4){
                knt = 0;

            }
        }
        outfile << QString("</p><p>");
    }
    fileout.close();
}

void measures::BTNresetangles_clicked(){
    anglecenters->clear();
    lastselectangles.clear();
    selectangles.clear();
    anglestext->clear();
    anglestext->append("  ");
    anglesprinttext.clear();
    LBLlastangles->setVisible(false);
    CHKshowangwin->setChecked(false);
    update_measuresInfo();
    emit reset_angles();
}

void measures::BTNresetdihedrals_clicked(){
    dihedralcenters->clear();
    lastselectdihedrals.clear();
    selectdihedrals.clear();
    dihedralstext->clear();
    dihedralstext->append("  ");
    dihedralsprinttext.clear();
    LBLlastdihedrals->setVisible(false);
    CHKshowdihedwin->setChecked(false);
    update_measuresInfo();
    emit reset_dihedrals();
}

void measures::BTNresetdistances_clicked(){
    distancecenters->clear();
    lastselectdist.clear();
    selectdist.clear();
    distancestext->clear();
    distancestext->append("  ");
    distancesprinttext.clear();
    LBLlastdistances->setVisible(false);
    CHKshowdistwin->setChecked(false);
    update_measuresInfo();
    emit reset_distances();
}

void measures::BTNsurfcolor_clicked(){
    QColor currcolor = BTNsurfcolor->getbkgColor();
    QColor newcolor = QColorDialog::getColor(currcolor, this);
    if(newcolor.isValid()) {
        BTNsurfcolor->setColor(&newcolor);
    }
    newcolor.setAlpha(SLDopacity->value());
    emit dihedralplanes_color(newcolor);
}

void measures::CHKdrawarcs_changed(){
    emit draw_arcs(CHKdrawarcs->isChecked());
}

void measures::CHKdrawlines_changed(){
    emit draw_lines(CHKdrawlines->isChecked());
}

void measures::CHKdisttranspbkg_changed(){
    emit dist_transpbkg(CHKdisttranspbkg->isChecked());
}

void measures::CHKshowangles_changed(){
    emit show_angles(CHKshowangles->isChecked());
}

void measures::CHKshowdihedrals_changed(){
    emit show_dihedrals(CHKshowdihedrals->isChecked());
}

void measures::CHKshowdistances_changed(){
    SPBdistvshift->setVisible(CHKshowdistances->isChecked());
    LBLdistvshift->setVisible(CHKshowdistances->isChecked());
    emit show_distances(CHKshowdistances->isChecked());
}

void measures::CHKshowangwin_changed(){
    update_measuresInfo();
}

void measures::CHKshowdihedwin_changed(){
    update_measuresInfo();
}

void measures::CHKshowdistwin_changed(){
    update_measuresInfo();
}

void measures::clear_systemnames(){
    systemnames.clear();
}

void measures::close_measuresInfoWindow(){
    for (int i = 0 ; i < connectionsInfo.size() ; i++){
        QObject::disconnect(connectionsInfo.at(i));
    }
    connectionsInfo.clear();
    if (measuresInfo){
        measuresInfopos = QPoint(measuresInfo->x(),measuresInfo->y());
        delete measuresInfo;
        measuresInfo = nullpointer;
    }
    CHKshowangwin->setChecked(false);
    CHKshowdihedwin->setChecked(false);
    CHKshowdistwin->setChecked(false);
}

void measures::distance_remove(QVector<centerData> *dstcenters, int i){
    dstcenters->remove(i+1);
    dstcenters->remove(i);
    distancecenters = dstcenters;
    distancestext->clear();
    distancestext->append("  ");
    distancesprinttext.clear();
    for (int i = 0 ; i < distancecenters->length() ; i += 2){
        float dst = QVector3D(mtrsf->at(distancecenters->at(i+1).molecule)*distancecenters->at(i+1).xyz
                              -mtrsf->at(distancecenters->at(i).molecule)*distancecenters->at(i).xyz).length();
        distancestext->append("("+selectdist.at(i/2).join(",")+")"
                + QString("<sub>%1</sub>").arg(distancecenters->at(i).molecule + 1) + " = ");
        QString distvalue;
        if (RBTangstrom->isChecked())
            distvalue = QString::number( dst * BOHR_TO_ANGSTROM, 'g', SPBdistprecision->value())
                    + QChar(0x212B) + QString(" ; ");
        else
            distvalue = QString::number( dst, 'g', SPBdistprecision->value())
                    + QString("a") + QChar(0x2080) + QString(" ; ");
        distancestext->append(distvalue);
        distancesprinttext.append(QString("("+selectdist.at(i/2).join(",")+") = ") + distvalue);
        if (i > 0 && i%8 == 0){
            distancestext->chop(2);
            distancestext->append(" <br/> ");
        }

    }
    if (measuresInfo){
        update_measuresInfo();
    }
}

void measures::CMBanglestype_changed(int a){
    emit angles_type(a);
}

void measures::CMBlinestype_changed(int a){
    emit lines_type(a);
}

void measures::update_angles(QVector<centerData> *angcenters, QVector<QMatrix4x4> *transform_matrices){
    if (!CHKshowangwin->isChecked()) return;
    anglecenters = angcenters;
    mtrsf = transform_matrices;
    anglestext->clear();
    anglesprinttext.clear();
    for(int i = 0, k = 1 ; i < anglecenters->length() ; i += 3, k++){
        if (i >= anglecenters->length()-2) break;
        QVector3D A3D = QVector3D(mtrsf->at(anglecenters->at(i).molecule)*QVector4D(anglecenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(mtrsf->at(anglecenters->at(i+1).molecule)*QVector4D(anglecenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(mtrsf->at(anglecenters->at(i+2).molecule)*QVector4D(anglecenters->at(i+2).xyz,1));
        double angle = 180. * acos(QVector3D::dotProduct((A3D-B3D).normalized(),(C3D-B3D).normalized())) / M_PI;

        QString a = QString("(%1<sub>%2</sub>,%3<sub>%4</sub>,%5<sub>%6</sub>)")
                .arg(anglecenters->at(i).symbol).arg(anglecenters->at(i).molecule + 1)
                .arg(anglecenters->at(i+1).symbol).arg(anglecenters->at(i+1).molecule + 1)
                .arg(anglecenters->at(i+2).symbol).arg(anglecenters->at(i+2).molecule + 1);
        anglestext->append(a + " = " + QString::number(angle, 'g', SPBanglesprecision->value())
                              + QChar(0260) + QString(" ; "));
//        if (anglecenters->length()%10 == 0){
        if (k%3 == 0){
            anglestext->chop(2);
            anglestext->append(" <br/> ");
        }
        anglesprinttext.append(a + " = " + QString::number(angle, 'g', SPBanglesprecision->value()) + QChar(0260) + QString(" ; "));
    }
    if (measuresInfo) measuresInfo->update_angles(anglestext);
}

void measures::update_dihedrals(QVector<centerData> *dihcenters, QVector<QMatrix4x4> *transform_matrices){
    if (!CHKshowdihedwin->isChecked()) return;
    dihedralcenters = dihcenters;
    mtrsf = transform_matrices;
    dihedralstext->clear();
    dihedralsprinttext.clear();
    for(int i = 0, k = 1 ; i < dihedralcenters->length() ; i += 4, k++){
        if (i >= dihedralcenters->length()-3) break;
        QVector3D A3D = QVector3D(mtrsf->at(dihedralcenters->at(i).molecule)*QVector4D(dihedralcenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(mtrsf->at(dihedralcenters->at(i+1).molecule)*QVector4D(dihedralcenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(mtrsf->at(dihedralcenters->at(i+2).molecule)*QVector4D(dihedralcenters->at(i+2).xyz,1));
        QVector3D D3D = QVector3D(mtrsf->at(dihedralcenters->at(i+3).molecule)*QVector4D(dihedralcenters->at(i+3).xyz,1));
        QVector3D ABC = QVector3D::crossProduct((C3D-B3D),(A3D-B3D));
        QVector3D ABD = QVector3D::crossProduct((D3D-B3D),(A3D-B3D));
        double aux = QVector3D::dotProduct(ABC,ABD)/(ABC.length()*ABD.length());
        double angle;
        if (qAbs(qAbs(aux)-1.) > 1.e-6)
            angle = 180. * std::acos(aux) / M_PI;
        else{
            if (aux > 0.)
                angle = 0.;
            else
                angle = M_PI;
        }
        QString a = QString("(%1<sub>%2</sub>,%3<sub>%4</sub>,%5<sub>%6</sub>,%7<sub>%8</sub>)")
                .arg(dihedralcenters->at(i).symbol).arg(dihedralcenters->at(i).molecule + 1)
                .arg(dihedralcenters->at(i+1).symbol).arg(dihedralcenters->at(i+1).molecule + 1)
                .arg(dihedralcenters->at(i+2).symbol).arg(dihedralcenters->at(i+2).molecule + 1)
                .arg(dihedralcenters->at(i+3).symbol).arg(dihedralcenters->at(i+3).molecule + 1);
        dihedralstext->append(a + " = " + QString::number(angle, 'g', SPBdihedralsprecision->value())
                              + QChar(0260) + QString(" ; "));
//        if (dihedralcenters->length()%10 == 0){
        if (k%2 == 0){
            dihedralstext->chop(2);
            dihedralstext->append(" <br/> ");
        }
        dihedralsprinttext.append(a + " = " + QString::number(angle, 'g', SPBdihedralsprecision->value()) + QChar(0260) + QString(" ; "));     
    }
    if (measuresInfo) measuresInfo->update_dihedrals(dihedralstext);
}

void measures::update_distances(QVector<centerData> *dstcenters, QVector<QMatrix4x4> *transform_matrices){
//    if (!CHKshowdistwin->isChecked()) return;
    distancecenters = dstcenters;
    mtrsf = transform_matrices;
    distancestext->clear();
    distancesprinttext.clear();
    for(int i = 0, k = 1 ; i < distancecenters->length() ; i += 2, k++){
        if (i == distancecenters->length()-1) break;
        float dst = QVector3D( mtrsf->at(distancecenters->at(i+1).molecule)*QVector4D(distancecenters->at(i+1).xyz,1)
                        - mtrsf->at(distancecenters->at(i).molecule)*QVector4D(distancecenters->at(i).xyz,1) ).length();

        QString a = QString("(%1<sub>%2</sub>,%3<sub>%4</sub>)").arg(distancecenters->at(i).symbol)
                .arg(distancecenters->at(i).molecule + 1).arg(distancecenters->at(i+1).symbol)
                .arg(distancecenters->at(i+1).molecule + 1);
        distancestext->append(a + " = ");
        QString distvalue;
        if (RBTangstrom->isChecked())
            distvalue = QString::number( dst * BOHR_TO_ANGSTROM, 'g', SPBdistprecision->value())
                    + QChar(0x212B) + QString(" ; ");
        else
            distvalue = QString::number( dst, 'g', SPBdistprecision->value())
                    + QString("a") + QChar(0x2080) + QString(" ; ");
        distancestext->append(distvalue);
        if (k%4 == 0){
            distancestext->chop(2);
            distancestext->append(" <br/> ");
        }
        distancesprinttext.append(a + " = " + distvalue);
    }
    if (measuresInfo) measuresInfo->update_distances(distancestext);
}

void measures::emitupdateRightMenu(){
    emit updateRightMenu();
}

bool measures::getmeasureangles(){
    return RBTmeasureangle->isChecked();
}

bool measures::getmeasuredihedrals(){
    return RBTmeasuredihed->isChecked();
}

bool measures::getmeasuredistances(){
    return RBTmeasuredist->isChecked();
}

bool measures::getmeasurenone(){
    return RBTmeasurenone->isChecked();
}

bool measures::QDLmeasures_isVisible(){
    if (QDLmeasures)
        return QDLmeasures->isVisible();
    else
        return false;
}

void measures::LBLangles_add(QVector<centerData> * angcnts, QVector<QMatrix4x4> *transform_matrices){
    anglecenters = angcnts;
    mtrsf = transform_matrices;
    if (anglecenters->length() == 0) return;
    if (anglecenters->last().number < 0){
        for (int i = 0 ; i < qMin(lastselectangles.length(),anglecenters->length()) ; i++){
            anglecenters->removeLast();
        }
        lastselectangles.clear();
        LBLlastangles->setText(QString("Last selection: "));
        return;
    }
    QString a = anglecenters->last().symbol + QString("<sub>%1</sub>").arg(anglecenters->last().molecule + 1);
    lastselectangles.append(a);
    LBLlastangles->setText(QString("Last selection: ")+lastselectangles.join(", "));
    LBLlastangles->setStyleSheet("QLabel { color : red; }");
    LBLlastangles->setVisible(true);
    if (lastselectangles.length() == 3){
        int i = anglecenters->length() - 3;
        QVector3D A3D = QVector3D(mtrsf->at(anglecenters->at(i).molecule)*QVector4D(anglecenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(mtrsf->at(anglecenters->at(i+1).molecule)*QVector4D(anglecenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(mtrsf->at(anglecenters->at(i+2).molecule)*QVector4D(anglecenters->at(i+2).xyz,1));
        double angle = 180. * acos(QVector3D::dotProduct((A3D-B3D).normalized(),(C3D-B3D).normalized())) / M_PI;
        anglestext->append("("+lastselectangles.join(",")+")" + " = "
                + QString::number( std::abs(angle), 'g', SPBanglesprecision->value()) + QChar(0260) + QString(" ; "));
        anglesprinttext.append(QString("("+lastselectangles.join(",")+")" + " = "
                + QString::number( std::abs(angle), 'g', SPBanglesprecision->value()) + QChar(0260) + QString(" ; ")));
        if (anglecenters->length()%12 == 0){
            anglestext->chop(2);
            anglestext->append(" <br/> ");
        }
        selectangles.append(lastselectangles);
        lastselectangles.clear();
        if (measuresInfo)
            measuresInfo->update_angles(anglestext);
    }
}

void measures::LBLdihedrals_add(QVector<centerData> *dihcenters, QVector<QMatrix4x4> *transform_matrices){
    dihedralcenters = dihcenters;
    mtrsf = transform_matrices;
    if (dihedralcenters->length() == 0) return;
    if (dihedralcenters->last().number < 0){
        for (int i = 0 ; i <= qMin(lastselectdihedrals.length(),dihedralcenters->length()) ; i++){
            dihedralcenters->removeLast();
        }
        lastselectdihedrals.clear();
        LBLlastdihedrals->setText(QString("Last selection: "));
        return;
    }
    QString a = dihedralcenters->last().symbol + QString("<sub>%1</sub>").arg(dihedralcenters->last().molecule + 1);
    lastselectdihedrals.append(a);
    LBLlastdihedrals->setText(QString("Last selection: ")+lastselectdihedrals.join(", "));
    LBLlastdihedrals->setStyleSheet("QLabel { color : red; }");
    LBLlastdihedrals->setVisible(true);
    if (lastselectdihedrals.length() == 4){
        int i = dihedralcenters->length() - 4;
        QVector3D A3D = QVector3D(mtrsf->at(dihedralcenters->at(i).molecule)*QVector4D(dihedralcenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(mtrsf->at(dihedralcenters->at(i+1).molecule)*QVector4D(dihedralcenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(mtrsf->at(dihedralcenters->at(i+2).molecule)*QVector4D(dihedralcenters->at(i+2).xyz,1));
        QVector3D D3D = QVector3D(mtrsf->at(dihedralcenters->at(i+3).molecule)*QVector4D(dihedralcenters->at(i+3).xyz,1));
        QVector3D ABC = QVector3D::crossProduct((C3D-B3D),(A3D-B3D));
        QVector3D ABD = QVector3D::crossProduct((D3D-B3D),(A3D-B3D));
        double angle = 180. * std::acos(QVector3D::dotProduct(ABC,ABD)/(ABC.length()*ABD.length())) / M_PI;
        dihedralstext->append("("+lastselectdihedrals.join(",")+") = "
                + QString::number( std::abs(angle), 'g', SPBdihedralsprecision->value()) + QChar(0260) + QString(" ; "));
        dihedralsprinttext.append(QString("("+lastselectdihedrals.join(",")+")" + " = "
                + QString::number( std::abs(angle), 'g', SPBdihedralsprecision->value()) + QChar(0260) + QString(" ; ")));
        if (dihedralcenters->length()%12 == 0){
            dihedralstext->chop(2);
            dihedralstext->append(" <br/> ");
        }
        selectdihedrals.append(lastselectdihedrals);
        lastselectdihedrals.clear();
        if (measuresInfo)
            measuresInfo->update_dihedrals(dihedralstext);
    }
}

void measures::LBLdistances_add(QVector<centerData> *dstcenters, QVector<QMatrix4x4> *transform_matrices){
    distancecenters = dstcenters;
    mtrsf = transform_matrices;
    if (distancecenters->length() == 0) return;
    if (distancecenters->last().number < 0){
        for (int i = 0 ; i <= qMin(lastselectdist.length(),distancecenters->length()) ; i++){
            distancecenters->removeLast();
        }
        lastselectdist.clear();
        LBLlastdistances->setText(QString("Last selection: "));
        return;
    }
    QString a = distancecenters->last().symbol + QString("<sub>%1</sub>").arg(distancecenters->last().molecule + 1);
    if (lastselectdist.length() >= 2){
        lastselectdist.clear();
    }
    lastselectdist.append(a);
    LBLlastdistances->setText(QString("Last selection: ")+lastselectdist.join(", "));
    LBLlastdistances->setStyleSheet("QLabel { color : red; }");
    LBLlastdistances->setVisible(true);
    if (lastselectdist.length() == 2){
        int i = distancecenters->length() - 2;
        float dst = QVector3D( mtrsf->at(distancecenters->at(i+1).molecule)*QVector4D(distancecenters->at(i+1).xyz,1)
                        - mtrsf->at(distancecenters->at(i).molecule)*QVector4D(distancecenters->at(i).xyz,1) ).length();
        distancestext->append("("+lastselectdist.join(",")+")" + " = ");
        QString distvalue;
        if (RBTangstrom->isChecked())
            distvalue = QString::number( dst * BOHR_TO_ANGSTROM, 'g', SPBdistprecision->value())
                    + QChar(0x212B) + QString(" ; ");
        else
            distvalue = QString::number( dst, 'g', SPBdistprecision->value())
                    + QString("a") + QChar(0x2080) + QString(" ; ");
        distancestext->append(distvalue);
        if (distancecenters->length()%10 == 0){
            distancestext->chop(2);
            distancestext->append(" <br/> ");
        }
        distancesprinttext.append(QString("("+lastselectdist.join(",")+") = ") + distvalue);
        selectdist.append(lastselectdist);
        if (measuresInfo)
            measuresInfo->update_distances(distancestext);
    }
}

void measures::QDLmeasures_raise(){
    if (QDLmeasures)
        QDLmeasures->raise();
}

void measures::RBTangstrom_changed(){
    emit angstrom_changed(RBTangstrom->isChecked());
    distancestext->clear();
    distancestext->append("  ");
    distancesprinttext.clear();
    for (int i = 0 ; i < distancecenters->length() ; i += 2){
        if (distancecenters->length() <= i+1) break;  // if number of centers is odd, last center is skipped
        float dst = QVector3D( mtrsf->at(distancecenters->at(i+1).molecule)*QVector4D(distancecenters->at(i+1).xyz,1)
                        - mtrsf->at(distancecenters->at(i).molecule)*QVector4D(distancecenters->at(i).xyz,1) ).length();
        distancestext->append("("+selectdist.at(i/2).join(",")+")" + " = ");
        QString distvalue;
        if (RBTangstrom->isChecked())
            distvalue = QString::number( dst * BOHR_TO_ANGSTROM, 'g', SPBdistprecision->value())
                    + QChar(0x212B) + QString(" ; ");
        else
            distvalue = QString::number( dst, 'g', SPBdistprecision->value())
                    + QString("a") + QChar(0x2080) + QString(" ; ");
        distancestext->append(distvalue);
        distancesprinttext.append(QString("("+selectdist.at(i/2).join(",")+") = ") + distvalue);
        if (i > 0 && i%8 == 0){
            distancestext->chop(2);
            distancestext->append(" <br/> ");
        }
    }
    if (measuresInfo)
        measuresInfo->update_distances(distancestext);
}

void measures::RBTmeasureangle_changed(){
    if (RBTmeasureangle->isChecked()){
        FRMmeasureangles->setVisible(true);
        emit measure_angles(true);
    }
    else{
        // Prevents troubled caused by incomplete selection
        int knt = anglecenters->length()%3;
        if ( knt != 0){
            for (int k = 0 ; k < knt ; k++){
                anglecenters->removeLast();
            }
            lastselectangles.clear();
            LBLlastangles->setText(QString(""));
        }
        FRMmeasureangles->setVisible(false);
        QDLmeasures->adjustSize();
        emit measure_angles(false);
    }
}

void measures::RBTmeasuredihed_changed(){
    if (RBTmeasuredihed->isChecked()){
        FRMmeasuredihedrals->setVisible(true);
        emit measure_dihedrals(true);
    }
    else{
        // Prevents troubled caused by incomplete selection
        int knt = dihedralcenters->length()%4;
        if (knt != 0){
            for (int k = 0 ; k < knt ; k++){
                dihedralcenters->removeLast();
            }
            lastselectdihedrals.clear();
            LBLlastdihedrals->setText(QString(""));
        }
        FRMmeasuredihedrals->setVisible(false);
        QDLmeasures->adjustSize();
        emit measure_dihedrals(false);
    }
}

void measures::RBTmeasuredist_changed(){
    if (RBTmeasuredist->isChecked()){
        FRMmeasuredistances->setVisible(true);
        emit measure_distances(true);
    }
    else{
        // Prevents troubled caused by incomplete selection
        if (distancecenters->length()%2 == 1){
            distancecenters->removeLast();
            lastselectdist.clear();
            LBLlastdistances->setText(QString(""));
        }
        FRMmeasuredistances->setVisible(false);
        QDLmeasures->adjustSize();
        emit measure_distances(false);
    }
}

void measures::RBTmeasurenone_changed(){
    if (RBTmeasurenone->isChecked()){
        FRMmeasureangles->setVisible(false);
        FRMmeasuredihedrals->setVisible(false);
        FRMmeasuredistances->setVisible(false);
        QDLmeasures->adjustSize();
        emit activate_measure(false);
    }
    else{
        emit activate_measure(true);
    }
}


void measures::reset_all(){
    if (anglecenters) anglecenters->clear();
    lastselectangles.clear();
    selectangles.clear();
    if (anglestext){
        anglestext->clear();
        anglestext->append("  ");
    }
    anglesprinttext.clear();
    LBLlastangles->setVisible(false);
    CHKshowangwin->setChecked(false);
    if (dihedralcenters) dihedralcenters->clear();
    lastselectdihedrals.clear();
    selectdihedrals.clear();
    if (dihedralstext){
        dihedralstext->clear();
        dihedralstext->append("  ");
    }
    dihedralsprinttext.clear();
    LBLlastdihedrals->setVisible(false);
    CHKshowdihedwin->setChecked(false);
    if (distancecenters) distancecenters->clear();
    lastselectdist.clear();
    selectdist.clear();
    if (distancestext){
        distancestext->clear();
        distancestext->append("  ");
    }
    distancesprinttext.clear();
    LBLlastdistances->setVisible(false);
    RBTangstrom->setChecked(false);
    CHKshowdistwin->setChecked(false);
    RBTmeasurenone->setChecked(true);
    update_measuresInfo();
    emit reset_distances();
    emit reset_angles();
    emit reset_dihedrals();
}

void measures::resetlastselectangles(){
    lastselectangles.clear();
}

void measures::resetlastselectdihedrals(){
    lastselectdihedrals.clear();
}

void measures::resetlastselectdist(){
    lastselectdist.clear();
}

void measures::set_molecules(QStringList mol){
    molecules = mol;
}

void measures::set_ProjectFolder(QString a){
    ProjectFolder = a;
}

void measures::SLDopacity_released(){
    QColor newcolor = BTNsurfcolor->getbkgColor();
    newcolor.setAlpha(SLDopacity->value());
    emit dihedralplanes_color(newcolor);
}

void measures::SPBanglesprecision_changed(int a){
    anglestext->clear();
    anglestext->append("  ");
    anglesprinttext.clear();
    for (int i = 0 ; i < anglecenters->length() ; i += 3){
        QVector3D A3D = QVector3D(mtrsf->at(anglecenters->at(i).molecule)*QVector4D(anglecenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(mtrsf->at(anglecenters->at(i+1).molecule)*QVector4D(anglecenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(mtrsf->at(anglecenters->at(i+2).molecule)*QVector4D(anglecenters->at(i+2).xyz,1));
        float angle = 180. * acos(QVector3D::dotProduct((A3D-B3D).normalized(),(C3D-B3D).normalized())) / M_PI;
        anglestext->append("("+selectangles.at(i/3).join(",")+")"
                + QString("<sub>%1</sub>").arg(anglecenters->at(i).molecule + 1) + " = "
                + QString::number( std::abs(angle), 'g', SPBanglesprecision->value()) + QChar(0260) + QString(" ; "));
        anglesprinttext.append(QString("("+selectangles.at(i/3).join(",")+")" + " = "
                + QString::number( std::abs(angle), 'g', SPBanglesprecision->value()) + QChar(0260) + QString(" ; ")));
        if (i > 0 && i%9 == 0){
            anglestext->chop(2);
            anglestext->append(" <br/> ");
        }
    }
    if (measuresInfo)
        measuresInfo->update_angles(anglestext);
    emit angles_precision(a);
}

void measures::SPBangleswidth_changed(int a){
    emit angles_width(a);
}

void measures::SPBarcradius_changed(int a){
    emit arc_radius(a);
}

void measures::SPBdihedralsprecision_changed(int a){ 
    dihedralstext->clear();
    dihedralstext->append("  ");
    dihedralsprinttext.clear();
    for (int i = 0 ; i < dihedralcenters->length() ; i+= 4){
        QVector3D A3D = QVector3D(mtrsf->at(dihedralcenters->at(i).molecule)*QVector4D(dihedralcenters->at(i).xyz,1));
        QVector3D B3D = QVector3D(mtrsf->at(dihedralcenters->at(i+1).molecule)*QVector4D(dihedralcenters->at(i+1).xyz,1));
        QVector3D C3D = QVector3D(mtrsf->at(dihedralcenters->at(i+2).molecule)*QVector4D(dihedralcenters->at(i+2).xyz,1));
        QVector3D D3D = QVector3D(mtrsf->at(dihedralcenters->at(i+3).molecule)*QVector4D(dihedralcenters->at(i+3).xyz,1));
        QVector3D ABC = QVector3D::crossProduct((C3D-B3D),(A3D-B3D));
        QVector3D ABD = QVector3D::crossProduct((D3D-B3D),(A3D-B3D));
        double angle = 180. * std::acos(QVector3D::dotProduct(ABC,ABD)/(ABC.length()*ABD.length())) / M_PI;
        dihedralstext->append("("+selectdihedrals.at(i/4).join(",")+")"
                + QString("<sub>%1</sub>").arg(dihedralcenters->at(i).molecule + 1) + " = "
                + QString::number( std::abs(angle), 'g', SPBdihedralsprecision->value()) + QChar(0260) + QString(" ; "));
        dihedralsprinttext.append(QString("("+selectdihedrals.at(i/4).join(",")+")" + " = "
                + QString::number( std::abs(angle), 'g', SPBdihedralsprecision->value()) + QChar(0260) + QString(" ; ")));
        if (i > 0 && i%8 == 0){
            dihedralstext->chop(2);
            dihedralstext->append(" <br/> ");
        }
    }
    if (measuresInfo)
        measuresInfo->update_dihedrals(dihedralstext);
    emit dihedrals_precision(a);
}

void measures::SPBdistprecision_changed(int a){
    distancestext->clear();
    distancestext->append("  ");
    distancesprinttext.clear();
    for (int i = 0 ; i < distancecenters->length() ; i += 2){
        float dst = QVector3D( mtrsf->at(distancecenters->at(i+1).molecule)*QVector4D(distancecenters->at(i+1).xyz,1)
                        - mtrsf->at(distancecenters->at(i).molecule)*QVector4D(distancecenters->at(i).xyz,1) ).length();
        distancestext->append("("+selectdist.at(i/2).join(",")+")" + " = ");
        QString distvalue;
        if (RBTangstrom->isChecked())
            distvalue = QString::number( dst * BOHR_TO_ANGSTROM, 'g', SPBdistprecision->value())
                    + QChar(0x212B) + QString(" ; ");
        else
            distvalue = QString::number( dst, 'g', SPBdistprecision->value())
                    + QString("a") + QChar(0x2080) + QString(" ; ");
        distancestext->append(distvalue);
        distancesprinttext.append(QString("("+selectdist.at(i/2).join(",")+") = ") + distvalue);
        if (i > 0 && i%8 == 0){
            distancestext->chop(2);
            distancestext->append(" <br/> ");
        }

    }
//    if (measuresInfo){
//        update_measuresInfo();
//    }
    if (measuresInfo)
        measuresInfo->update_distances(distancestext);
    emit dist_precision(a);
}

void measures::SPBdistvshift_changed(int a){
    emit dist_vshift(a);
}

void measures::SPBlineswidth_changed(int a){
    emit lines_width(a);
}

void measures::update_measuresInfo(){
    if (measuresInfo){
        measuresInfopos = measuresInfo->pos();
        delete measuresInfo;
        measuresInfo = nullpointer;
    }
    if (CHKshowdistwin->isChecked() || CHKshowangwin->isChecked() || CHKshowdihedwin->isChecked()){
        measuresInfo = new measuresInfoWindow(molecules, CHKshowdistwin->isChecked(), distancestext,
                CHKshowangwin->isChecked(), anglestext, CHKshowdihedwin->isChecked(), dihedralstext,this);
        measuresInfo->setWindowTitle("Geometry measures info");
        measuresInfo->setAttribute( Qt::WA_DeleteOnClose );
        measuresInfo->move(measuresInfopos);
        connectionsInfo << connect(measuresInfo, SIGNAL(window_closed()), this, SLOT(close_measuresInfoWindow()));
        measuresInfo->show();
    }
}

void measures::update_measure_centers(QVector<QVector3D>* xyz){
    if (xyz->length() < 1)
        return;
    if (distancecenters->length()>0){
        for (int i = 0 ; i < distancecenters->length() ; i++){
            if (distancecenters->at(i).type > 0)
                continue;
            centerData cnt;
            cnt.type = distancecenters->at(i).type;
            cnt.molecule = distancecenters->at(i).molecule;
            cnt.number = distancecenters->at(i).number;
            cnt.znuc = distancecenters->at(i).znuc;
            cnt.x = distancecenters->at(i).x;
            cnt.y = distancecenters->at(i).y;
            cnt.symbol = distancecenters->at(i).symbol;
            cnt.xyz.setX(xyz->at(cnt.number).x());
            cnt.xyz.setY(xyz->at(cnt.number).y());
            cnt.xyz.setZ(xyz->at(cnt.number).z());
            cnt.mtrsf = distancecenters->at(i).mtrsf;
            distancecenters->replace(i,cnt);
//            if (i == 1) qDebug() << "i = " << i << " center = " << cnt.xyz;
        }
    }
    if (anglecenters->length()>0){
        for (int i = 0 ; i < anglecenters->length() ; i++){
            if (anglecenters->at(i).type > 0)
                continue;
            centerData cnt;
            cnt.type = anglecenters->at(i).type;
            cnt.molecule = anglecenters->at(i).molecule;
            cnt.number = anglecenters->at(i).number;
            cnt.znuc = anglecenters->at(i).znuc;
            cnt.x = anglecenters->at(i).x;
            cnt.y = anglecenters->at(i).y;
            cnt.symbol = anglecenters->at(i).symbol;
            cnt.xyz.setX(xyz->at(cnt.number).x());
            cnt.xyz.setY(xyz->at(cnt.number).y());
            cnt.xyz.setZ(xyz->at(cnt.number).z());
            cnt.mtrsf = anglecenters->at(i).mtrsf;
            anglecenters->replace(i,cnt);
        }
    }
    for (int i = 0 ; i < dihedralcenters->length() ; i++){
        if (dihedralcenters->at(i).type > 0)
            continue;
        centerData cnt;
        cnt.type = dihedralcenters->at(i).type;
        cnt.molecule = dihedralcenters->at(i).molecule;
        cnt.number = dihedralcenters->at(i).number;
        cnt.znuc = dihedralcenters->at(i).znuc;
        cnt.x = dihedralcenters->at(i).x;
        cnt.y = dihedralcenters->at(i).y;
        cnt.symbol = dihedralcenters->at(i).symbol;
        cnt.xyz.setX(xyz->at(cnt.number).x());
        cnt.xyz.setY(xyz->at(cnt.number).y());
        cnt.xyz.setZ(xyz->at(cnt.number).z());
        cnt.mtrsf = dihedralcenters->at(i).mtrsf;
        dihedralcenters->replace(i,cnt);
    }
}


/*******************************************************************************************************/
/********************************  Class editMeasuresDialog  implementation  ***************************/
/*******************************************************************************************************/

editMeasuresDialog::editMeasuresDialog(QWidget *parent) : QDialog(parent)
{

}

editMeasuresDialog::~editMeasuresDialog(){

}

void editMeasuresDialog::reject(){
    emit closed();
}
void editMeasuresDialog::closeEvent(QCloseEvent *event){
    event->ignore();
    this->setVisible(false);
    emit closed();
    event->accept();
}

//----------------------------------------------------
//----------------------------------------------------
//
//  Functions for measuresInfoWindow class
//
//----------------------------------------------------
//----------------------------------------------------


measuresInfoWindow::measuresInfoWindow(QStringList molecules, bool showdist, QString *distances, bool showang, QString *angles,
    bool showdihed, QString *dihedrals, QWidget *parent = 0) : QDialog(parent)
{
    QGroupBox *FRMmolecules = new QGroupBox(tr("Molecules"));
    FRMmolecules->setFont(QFont("Helvetica", 12, QFont::Bold));
    FRMmolecules->setStyleSheet("QLabel { color : red; }");

    QGroupBox *FRMdistances = new QGroupBox(tr("Distances"));
    FRMdistances->setFont(QFont("Helvetica", 12, QFont::Bold));
    FRMdistances->setStyleSheet("QLabel { color : blue; }");
    FRMdistances->setVisible(showdist);

    QGroupBox *FRMangles = new QGroupBox(tr("Angles"));
    FRMangles->setFont(QFont("Helvetica", 12, QFont::Bold));
    FRMangles->setStyleSheet("QLabel { color : blue; }");
    FRMangles->setVisible(showang);

    QGroupBox *FRMdihedrals = new QGroupBox(tr("Dihedral angles"));
    FRMdihedrals->setStyleSheet("QGroupBox {title: padding 2 2000 2 2;}");
    FRMdihedrals->setFont(QFont("Helvetica", 12, QFont::Bold));
    FRMdihedrals->setStyleSheet("QLabel { color : blue; }");
    FRMdihedrals->setVisible(showdihed);

    LBLmolecules = new QLabel();
    LBLdistances = new QLabel();
    LBLangles = new QLabel();
    LBLdihedrals = new QLabel();
    QString listmolecules;
    for (int i = 0 ; i < molecules.length() ; i++){
        listmolecules.append(QString("%1: %2\n").arg(i+1,2).arg(molecules.at(i)));
    }
#if QT_VERSION < 0x050A00
    listmolecules.chop(1);
    LBLmolecules->setText(listmolecules);
    if ((*distances).length() > 2){
        QString straux = *distances;
        straux.chop(2);
        LBLdistances->setText(straux);
    }
    else
        LBLdistances->setText(QString(""));
    if ((*angles).length() > 2){
        QString strbux = *angles;
        strbux.chop(2);
        LBLangles->setText(strbux);
    }
    else
        LBLangles->setText(QString(""));
    if ((*dihedrals).length() > 2){
        QString strcux = *dihedrals;
        strcux.chop(2);
        LBLdihedrals->setText(strcux);
    }
    else
        LBLdihedrals->setText(QString(""));
#else
    LBLmolecules->setText(listmolecules.chopped(1));
    if ((*distances).length() > 2)
        LBLdistances->setText(QString(*distances).chopped(2));
    else
        LBLdistances->setText(QString(""));
    if ((*angles).length() > 2)
        LBLangles->setText(QString(*angles).chopped(2));
    else
        LBLangles->setText(QString(""));
    if ((*dihedrals).length() > 2)
        LBLdihedrals->setText(QString(*dihedrals).chopped(2));
    else
        LBLdihedrals->setText(QString(""));
#endif
    LBLmolecules->setFont(QFont("Helvetica", 12, QFont::Bold));
    LBLdistances->setFont(QFont("Helvetica", 12, QFont::Bold)); 
    LBLangles->setFont(QFont("Helvetica", 12, QFont::Bold));
    LBLdihedrals->setFont(QFont("Helvetica", 12, QFont::Bold));

    QPushButton *BTNclose = new QPushButton(tr("Close"));
    connect(BTNclose, SIGNAL(clicked()), this, SLOT(measuresInfoWindow_close()));
    QHBoxLayout *layout0=new QHBoxLayout(FRMmolecules);
    layout0->addWidget(LBLmolecules);
    layout0->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout1=new QHBoxLayout(FRMdistances);
    layout1->addWidget(LBLdistances);
    layout1->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout2=new QHBoxLayout(FRMangles);
    layout2->addWidget(LBLangles);
    layout2->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout3=new QHBoxLayout(FRMdihedrals);
    layout3->addWidget(LBLdihedrals);
    layout3->setAlignment(Qt::AlignCenter);

    QHBoxLayout *layout4=new QHBoxLayout();
    layout4->addStretch();
    layout4->addWidget(BTNclose);
    layout4->addStretch();

    QVBoxLayout *layout=new QVBoxLayout(this);
    layout->addWidget(FRMmolecules);
    layout->addWidget(FRMdistances);
    layout->addWidget(FRMangles);
    layout->addWidget(FRMdihedrals);
    layout->addLayout(layout4);
}


measuresInfoWindow::~measuresInfoWindow(){

}


void measuresInfoWindow::closeEvent(QCloseEvent *event){
    event->ignore();
    emit window_closed();
    event->accept();
}

void measuresInfoWindow::measuresInfoWindow_close(){
    emit window_closed();
}

void measuresInfoWindow::update_angles(QString* str){
    if ((*str).length() > 2){
#if QT_VERSION < 0x050A00
        QString string = *str;
        string.chop(2);
        LBLangles->setText(string);
#else
        LBLangles->setText(QString(*str).chopped(2));
#endif
    }
    else{
        LBLangles->setText("");
    }
}

void measuresInfoWindow::update_distances(QString* str){
//    if ((*str).length() > 2)
//        LBLdistances->setText(QString(*str).chopped(2));
    if ((*str).length() > 2){
#if QT_VERSION < 0x050A00
        QString string = *str;
        string.chop(2);
        LBLdistances->setText(string);
#else
        LBLdistances->setText(QString(*str).chopped(2));
#endif
    }
    else{
        LBLdistances->setText("");
    }
}

void measuresInfoWindow::update_dihedrals(QString * str){
//    if ((*str).length() > 2)
//        LBLdihedrals->setText(QString(*str).chopped(2));
    if ((*str).length() > 2){
#if QT_VERSION < 0x050A00
        QString string = *str;
        string.chop(2);
        LBLdihedrals->setText(string);
#else
        LBLdihedrals->setText(QString(*str).chopped(2));
#endif
    }
    else{
        LBLdihedrals->setText("");
    }
}
