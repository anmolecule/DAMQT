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
//
//  File:   Viewer2D.cpp
//
//      Last version: September 2018
//

#include "viewer2D.h"

#include <QtDebug>
#include <cmath>
#include "IniFile.h"
#define ANGSTROMTOBOHR 1.889725989

class QByteArray;

Viewer2D::Viewer2D(QString *name, QWidget *parent): QWidget(parent) {
    myDoubleValidator = new QDoubleValidator(nullpointer);
    myDoubleValidator->setLocale(QLocale::English);
    loadedfiles.clear();
    nhistpart.clear();
    frdata.clear();
    frads_data.clear();
    sgdata.clear();
    sgdataori.clear();
    sghist_data.clear();
    sgx.clear();
    sgx.clear();
    initpointers();
    superimpose = false;
    basinscolor = new QColor(Qt::red);
    bkgcolor = new QColor(Qt::white);
    bondscolor = QColor(0,0,255);
    centerscolor = QColor(0,0,0);
    centerslabelcolor = QColor(0,0,255);
    fonttitlecolor = QColor(0,0,255);
    fontXlabelcolor = QColor(0,0,255);
    fontYlabelcolor = QColor(0,0,255);
    fonttitle = QFont("Liberation Mono", 16);
    fontXlabel = QFont("Liberation Mono", 16);
    fontYlabel = QFont("Liberation Mono", 16);
    distthr = 1.e-3;
    minx = 0.0;
    maxx = 0.0;
    miny = 0.0;
    maxy = 0.0;
    ncurves = 0;
    nfradcolors = 0;
    nsghist = 0;
    nsghistcolors = -1;
    planeA = 0.0;
    planeB = 0.0;
    planeC = 1.0;
    planecase = 1;
    plotcps= false;
    position = QPoint(0,0);
    appendsghstrect = true;
    printhistfile = QString("");
    retrieve = false;
    resetsghstrect = true;
    size = QSize(1000,500);
    CaptureFolder = QString("");
    ProjectFolder = QString("");
    viewername = QString("unnamed");
    visible = true;
    smoothfactor = 5;
    zero = 0.0;

    set_viewername(QString("Plot no. %1").arg(*name));
    viewernumber = (*name).toInt();

    mainWin = new MainWindow2DViewer();
        
    grafica = new Plotter(this);
    grafica->setxmin(minx);
    grafica->setxmax(maxx);
    grafica->setymin(miny);
    grafica->setymax(maxy);
    grafica->setfonttitle(fonttitle);

    grafica->Pal->setColor(QPalette::Background, bkgcolor->rgb());
    grafica->setPalette(*(grafica->Pal));
    grafica->setAutoFillBackground(true);
    grafica->setbackgroundcolor(bkgcolor->rgb());
    grafica->setfontcontour(QFont("Liberation Mono", 7));
    grafica->setfontXlabel(fontXlabel);
    grafica->setfontYlabel(fontYlabel);
    grafica->setfonttitlecolor(fonttitlecolor);
    grafica->setfontXlabelcolor(fontXlabelcolor);
    grafica->setfontYlabelcolor(fontYlabelcolor);
    grafica->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    connect(grafica, SIGNAL(deletefrad(int)), this, SLOT(deletefrad(int)));
    connect(grafica, SIGNAL(deletesghist(int)), this, SLOT(deletesghist(int)));

    graficaArea = new QScrollArea();
    graficaArea->setWidget(grafica);
    graficaArea->setWidgetResizable(true);
    graficaArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    graficaArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    graficaArea->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    graficaArea->setMinimumSize(500, 300);
    
    QDockWidget *gdock = new QDockWidget(mainWin);
        
    gdock->setAllowedAreas(Qt::RightDockWidgetArea);
    gdock->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
    gdock->setFeatures(QDockWidget::NoDockWidgetFeatures);
    gdock->setFeatures(QDockWidget::DockWidgetFloatable | QDockWidget::DockWidgetMovable);
    
    toolBox = new menuViewer2D(gdock);
    toolBox->setMinimumSize(400,120);
    gdock->setWidget(toolBox);
    
//      Contours (only available for contour plots)
        
    page_contours = new QWidget();
    page_contours->setEnabled(true);
    page_contours_widgets();
    page_contours_layouts();
    toolBox->addItem(page_contours,QIcon(":/images/opciones.png"),tr("Contour plots"));

//      Field lines

    page_field = new QWidget();
    page_field_widgets();
    page_field_layouts();
    toolBox->addItem(page_field,QIcon(":/images/opciones.png"),tr("Field lines"));

//      Sigma hole histogram (only available for sigma hole histogram plots)

    page_sghistogram = new QWidget();
    page_sghistogram_widgets();
    page_sghistogram_layouts();
    toolBox->addItem(page_sghistogram,QIcon(":/images/opciones.png"),tr("MESP sigma hole histogram"));

//      Radial factors (only available for radial factors plots)
        
    page_frad = new QWidget();
    page_frad_widgets();
    page_frad_layouts();
    toolBox->addItem(page_frad,QIcon(":/images/opciones.png"),tr("Radial factors"));

//      Critical points

    page_CPs = new QWidget();
    page_CPs->setEnabled(false);
    page_CPs_widgets();
    page_CPs_layouts();
    toolBox->addItem(page_CPs,QIcon(":/images/opciones.png"),tr("Critical points"));

//      Basins

    page_basins = new QWidget();
    page_basins->setEnabled(false);
    page_basins_widgets();
    page_basins_layouts();
    toolBox->addItem(page_basins,QIcon(":/images/opciones.png"),tr("Basins"));

//      Options
    page_options = new QWidget();
    page_options->setEnabled(true);
    page_options_widgets();
    page_options_layouts();
    toolBox->addItem(page_options,QIcon(":/images/opciones.png"),tr("Options"));

//      Image capture
    page_capture = new QWidget();
    page_capture->setEnabled(true);
    page_capture_widgets();
    page_capture_layouts();
    toolBox->addItem(page_capture,QIcon(":/images/opciones.png"),tr("Image capture"));

//         Save/Retrieve settings
    page_save = new QWidget();
    page_save->setEnabled(true);
    page_save_widgets();
    page_save_layouts();
    toolBox->addItem(page_save,QIcon(":/images/opciones.png"),tr("Save/retrieve settings"));

//      mainWin settings
    mainWin->setCentralWidget(graficaArea);
    mainWin->setWindowTitle(viewername);
    mainWin->addDockWidget(Qt::RightDockWidgetArea,gdock);
    mainWin->resize(size);
    mainWin->move(position);
    connect(mainWin,SIGNAL(hideplotter()),this,SIGNAL(hideplotter()));
    mainWin->show();
}

Viewer2D::~Viewer2D() {
    if (grafica){
        delete grafica;
        grafica = nullpointer;
    }
    if (mainWin){
        delete mainWin;
        mainWin = nullpointer;
    }
}

void Viewer2D::initpointers()
{
    basinscolor = nullpointer;
    bkgcolor = nullpointer;
    data = nullpointer;
    elines = nullpointer;
    grafica = nullpointer;
    gridcolor = nullpointer;
    levels = nullpointer;

    BTNactualizarcontour = nullpointer;
    BTNbkgcolor = nullpointer;
    BTNbasins = nullpointer;
    BTNbasinscolor = nullpointer;
    BTNbondscolor = nullpointer;
    BTNcaptura = nullpointer;
    BTNcenterscolor = nullpointer;
    BTNcenterslabelcolor = nullpointer;
    for (int i = 0 ; i < max_cps ; i++){
        BTNcolorcps[i] = nullpointer;
    }
    BTNcolortitlefont = nullpointer;
    BTNcolorXlabelfont = nullpointer;
    BTNcolorYlabelfont = nullpointer;
    BTNcontour = nullpointer;
    BTNcontourslabelcolor = nullpointer;
    BTNcps = nullpointer;
    BTNedittitle = nullpointer;
    BTNeditXlabel = nullpointer;
    BTNeditYlabel = nullpointer;
    BTNexeczoomreg = nullpointer;
    BTNfield = nullpointer;
    BTNfile = nullpointer;
    BTNfontcenterslabel = nullpointer;
    BTNfontcontourslabel = nullpointer;
    BTNfontcurvelabels = nullpointer;
    BTNfontsghistlabels = nullpointer;
    BTNfonttitle = nullpointer;
    BTNfontXlabel = nullpointer;
    BTNfontYlabel = nullpointer;
    BTNfrads = nullpointer;
    BTNgridcolor = nullpointer;
    BTNnegativecontourcolor = nullpointer;
    BTNpositivecontourcolor = nullpointer;
    BTNprinthist = nullpointer;
    BTNsghist = nullpointer;
    BTNretrieveSettings = nullpointer;
    BTNsaveSettings = nullpointer;
    BTNsghist = nullpointer;
    BTNsmooth = nullpointer;
    BTNzerocontourcolor = nullpointer;

    CHKarrows = nullpointer;
    CHKautomaticticks = nullpointer;
    CHKaspectratio = nullpointer;
    CHKbasins = nullpointer;
    CHKbonds = nullpointer;
    CHKcenters = nullpointer;
    CHKcenterslabel = nullpointer;
    CHKcontourslabel = nullpointer;
    for (int i = 0 ; i < max_cps ; i++){
        CHKcps[i] = nullpointer;
    }
    CHKcurvelabel = nullpointer;
    CHKgrid = nullpointer;
    CHKmulticolor = nullpointer;
    CHKprinthist = nullpointer;
    CHKsettingsfile = nullpointer;
    CHKsettingstitle = nullpointer;
    CHKsettingsXlabel = nullpointer;
    CHKsettingsYlabel = nullpointer;
    CHKsettingsZoom = nullpointer;
    CHKsettingsMargins = nullpointer;
    CHKsettingscntPen = nullpointer;
    CHKsettingselinesPen = nullpointer;
    CHKsettingsbkg = nullpointer;
    CHKshowparthist = nullpointer;
    CHKshowfield = nullpointer;
    CHKsmoothlines = nullpointer;
    CHKtranspbg = nullpointer;
    CHKtranspbgcapture = nullpointer;
    CHKsghistlabel = nullpointer;
    CHKXhide = nullpointer;
    CHKXscalebottom = nullpointer;
    CHKXscaletop = nullpointer;
    CHKYscaleleft = nullpointer;
    CHKYscaleright = nullpointer;
    CHKticks = nullpointer;
    CHKtitlehide = nullpointer;
    CHKtitle = nullpointer;
    CHKtitlebox = nullpointer;
    CHKXlabel = nullpointer;
    CHKXlabelbox = nullpointer;
    CHKYhide = nullpointer;
    CHKYlabel = nullpointer;
    CHKYlabelbox = nullpointer;
    CHKYhorizontal = nullpointer;
    CHKYlabelhorizontal = nullpointer;

    FRMarrows = nullpointer;
    FRMaspectratio = nullpointer;
    FRMbasins = nullpointer;
    FRMbasinspen = nullpointer;
    FRMbkgcolor = nullpointer;
    FRMbonds = nullpointer;
    FRMcaptura = nullpointer;
    FRMcenters = nullpointer;
    FRMcntpen = nullpointer;
    FRMcontours = nullpointer;
    FRMcontourscolor = nullpointer;
    FRMcontourslabel = nullpointer;
    FRMcontournegstyle = nullpointer;
    FRMcontourpostyle = nullpointer;
    FRMcontourzerostyle = nullpointer;
    FRMcontourstyle = nullpointer;
    FRMcps = nullpointer;
    FRMcurvespen = nullpointer;
    FRMelinespen = nullpointer;
    FRMfile = nullpointer;
    FRMfontcurvelabels = nullpointer;
    FRMfontsghistlabels = nullpointer;
    FRMgridcolor = nullpointer;
    FRMgridstyle = nullpointer;
    FRMmargins = nullpointer;
    FRMradfactortype = nullpointer;
    FRMscales = nullpointer;
    FRMsettings = nullpointer;
    FRMsghistpen = nullpointer;
    FRMsmoothlines = nullpointer;
    FRMticks = nullpointer;
    FRMtitle = nullpointer;
    FRMtitlegroup = nullpointer;
    FRMXlabel = nullpointer;
    FRMXlabelgroup = nullpointer;
    FRMYlabel = nullpointer;
    FRMYlabelgroup = nullpointer;
    FRMzoomregion = nullpointer;

    graficaArea = nullpointer;

    LBLfile = nullpointer;
    LBLBottomMargin = nullpointer;
    LBLcontourslabel = nullpointer;
    LBLinf = nullpointer;
    LBLLeftMargin = nullpointer;
    LBLnegativecontourcolor = nullpointer;
    LBLbasins = nullpointer;
    LBLbasinspenwidth = nullpointer;
    LBLcntpenwidth = nullpointer;
    LBLcontourfile = nullpointer;
    for (int i = 0 ; i < max_cps ; i++){
        LBLcps[i] = nullpointer;
    }
    LBLcpsballradius = nullpointer;
    LBLcpsfile = nullpointer;
    LBLcurvespenwidth = nullpointer;
    LBLelinespenwidth = nullpointer;
    LBLfield = nullpointer;
    LBLfrads = nullpointer;
    LBLpor = nullpointer;
    LBLpositivecontourcolor = nullpointer;
    LBLpostitle = nullpointer;
    LBLposXlabel = nullpointer;
    LBLposYlabel = nullpointer;
    LBLradfactor = nullpointer;
    LBLr2flm = nullpointer;
    LBLr2l2flm = nullpointer;
    LBLRightMargin = nullpointer;
    LBLscaledef = nullpointer;
    LBLsghist = nullpointer;
    LBLsghistpenwidth = nullpointer;
    LBLprthist = nullpointer;
    LBLsmoothfactor = nullpointer;
    LBLsup = nullpointer;
    LBLtolerance = nullpointer;
    LBLTopMargin = nullpointer;
    LBLx = nullpointer;
    LBLXcifras = nullpointer;
    LBLXticks = nullpointer;
    LBLxtitle = nullpointer;
    LBLxXlabel = nullpointer;
    LBLxYlabel = nullpointer;
    LBLy = nullpointer;
    LBLYcifras = nullpointer;
    LBLYticks = nullpointer;
    LBLytitle = nullpointer;
    LBLyXlabel = nullpointer;
    LBLyYlabel = nullpointer;
    LBLYlabelhorizontal = nullpointer;
    LBLzerocontourcolor = nullpointer;
    LBLzerotolerance = nullpointer;

    myDoubleValidator = nullpointer;

    page_contours = nullpointer;
    page_field = nullpointer;
    page_frad = nullpointer;
    page_CPs = nullpointer;
    page_basins = nullpointer;
    page_options = nullpointer;
    page_capture = nullpointer;
    page_save = nullpointer;
    page_sghistogram = nullpointer;

    RBTdashed = nullpointer;
    RBTnegsolid = nullpointer;
    RBTnegdashed = nullpointer;
    RBTnegdotted = nullpointer;
    RBTnegdashdot = nullpointer;
    RBTpossolid = nullpointer;
    RBTposdashed = nullpointer;
    RBTposdotted = nullpointer;
    RBTposdashdot = nullpointer;
    RBTradfactor = nullpointer;
    RBTr2flm = nullpointer;
    RBTr2l2flm = nullpointer;
    RBTscaledef = nullpointer;
    RBTscreendef = nullpointer;
    RBTsolid = nullpointer;
    RBTuserdef = nullpointer;
    RBTzerosolid = nullpointer;
    RBTzerodashed = nullpointer;
    RBTzerodotted = nullpointer;
    RBTzerodashdot = nullpointer;

    SHTclist = nullpointer;

    SPBarrowssep = nullpointer;
    SPBarrowssize = nullpointer;
    SPBarrowsskew = nullpointer;
    SPBarrowswidth = nullpointer;
    SPBbondswidth = nullpointer;
    SPBBottomMargin = nullpointer;
    SPBcpsballradius = nullpointer;
    SPBshftcenterlabels = nullpointer;
    SPBcenterradius = nullpointer;
    SPBbasinspenwidth = nullpointer;
    SPBcntpenwidth = nullpointer;
    SPBcurvespenwidth = nullpointer;
    SPBelinespenwidth = nullpointer;
    SPBimagequality = nullpointer;
    SPBLeftMargin = nullpointer;
    SPBRightMargin = nullpointer;
    SPBsghistpenwidth = nullpointer;
    SPBsmoothfactor = nullpointer;
    SPBtolerance = nullpointer;
    SPBTopMargin = nullpointer;
    SPBXcifras = nullpointer;
    SPBXticks = nullpointer;
    SPBYcifras = nullpointer;
    SPBYticks = nullpointer;
    SPBzerotolerance = nullpointer;

    TXTbondthreshold = nullpointer;
    TXTdistancethreshold = nullpointer;
    TXTbasins = nullpointer;
    TXTcontour = nullpointer;
    TXTcps = nullpointer;
    TXThsize = nullpointer;
    TXTfield = nullpointer;
    TXTfile = nullpointer;
    TXTfrads = nullpointer;
    TXTprinthist = nullpointer;
    TXTscalesize = nullpointer;
    TXTsghist = nullpointer;
    TXTtitle = nullpointer;
    TXTvsize = nullpointer;
    TXTxinf = nullpointer;
    TXTXlabel = nullpointer;
    TXTYlabel = nullpointer;
    TXTxsup = nullpointer;
    TXTxtitle = nullpointer;
    TXTXlabel = nullpointer;
    TXTxXlabel = nullpointer;
    TXTxYlabel = nullpointer;
    TXTyinf = nullpointer;
    TXTysup = nullpointer;
    TXTytitle = nullpointer;
    TXTYlabel = nullpointer;
    TXTyXlabel = nullpointer;
    TXTyYlabel = nullpointer;

    Wtablacont = nullpointer;
}

/***************************************************************************/
/*                 PAGES WIDGETS AND LAYOUTS                               */
/***************************************************************************/

//    page_contours: Contour Plots
//    ============================

void Viewer2D::page_contours_widgets()
{
    FRMcontours = new QGroupBox(page_contours);
    FRMcontours->setTitle(tr("Level values"));

    TXTcontour = new QLineEdit();
    BTNcontour = new QToolButton();
    BTNcontour->setText(tr("..."));
    LBLcontourfile = new QLabel(tr("File"));
    connect(BTNcontour, SIGNAL(clicked()), this, SLOT(importcntfile_dialog()));

    Wtablacont = new QWidget(FRMcontours);
    QStringList QSLclist;
    QSLclist << tr("Value") << tr("Value") << tr("Value") << tr("Value");
    SHTclist = new Sheet(4, 3, 0,true, Wtablacont);
    SHTclist->setHeader(QSLclist);
    Wtablacont->setVisible(true);

    BTNactualizarcontour = new QPushButton(QIcon(":/images/actualizar.png"),tr("Update"));
    connect(BTNactualizarcontour, SIGNAL(clicked()), this, SLOT(BTNactualizarcontour_click()));


//        Contours labels

    FRMcontourslabel = new QGroupBox(page_contours);
    FRMcontourslabel->setTitle(tr("Contour labels"));

    CHKcontourslabel = new QCheckBox(tr("Show/hide"));
    CHKcontourslabel->setChecked(true);
    CHKcontourslabel->setEnabled(true);
    connect(CHKcontourslabel, SIGNAL(stateChanged(int)), this, SLOT(CHKcontourslabel_changed()));

//            Contours label Font

    BTNfontcontourslabel = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));

    connect(BTNfontcontourslabel, SIGNAL(clicked()), this, SLOT(BTNfontcontourslabel_click()));

//            Contours label color

    contourslabelcolor = QColor(0,0,255);
    grafica->setcontourslabelcolor(contourslabelcolor);
    BTNcontourslabelcolor = new ColorButton();
    BTNcontourslabelcolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNcontourslabelcolor->setText(tr("Color"));
    BTNcontourslabelcolor->setColor(&contourslabelcolor);
    BTNcontourslabelcolor->setEnabled(true);
    connect(BTNcontourslabelcolor, SIGNAL(clicked()), this, SLOT(BTNcontourslabelcolor_clicked()));

    LBLtolerance = new QLabel(tr("Tolerance")+" (%): ");
    SPBtolerance = new QSpinBox();
    SPBtolerance->setMinimum(60);
    SPBtolerance->setMaximum(100);
    SPBtolerance->setMaximumWidth(50);
    SPBtolerance->setValue((int) 100*grafica->getcontourtolerance());
    connect(SPBtolerance,SIGNAL(valueChanged(int)),this,SLOT(SPBtolerance_changed(int)));

    LBLzerotolerance = new QLabel(tr("Tolerance for zero")+": 10^ ");
    SPBzerotolerance = new QSpinBox();
    SPBzerotolerance->setMinimum(-10);
    SPBzerotolerance->setMaximum(0);
    SPBzerotolerance->setMaximumWidth(50);
    SPBzerotolerance->setValue((int) grafica->getzerotolerancepow());
    connect(SPBzerotolerance,SIGNAL(valueChanged(int)),this,SLOT(SPBzerotolerance_changed(int)));

//        Contours colors

    FRMcontourscolor = new QGroupBox(page_contours);
    FRMcontourscolor->setTitle(tr("Colors"));

    CHKmulticolor = new QCheckBox(tr("Multicolor"));
    CHKmulticolor->setChecked(false);
    CHKmulticolor->setEnabled(true);
    connect(CHKmulticolor, SIGNAL(stateChanged(int)), this, SLOT(CHKmulticolor_changed()));

    LBLpositivecontourcolor = new QLabel(tr("Positive levels")+": ");
    BTNpositivecontourcolor = new ColorButton();
    BTNpositivecontourcolor->setIcon(QIcon(":/images/colores48.png"));
    BTNpositivecontourcolor->setText(tr("Color"));
    QColor positivecontourcolor = grafica->getpositivecontourcolor();
    BTNpositivecontourcolor->setColor(&(positivecontourcolor));
    BTNpositivecontourcolor->setEnabled(true);
    positivecontourcolor = grafica->getpositivecontourcolor();
    connect(BTNpositivecontourcolor, SIGNAL(clicked()), this, SLOT(BTNpositivecontourcolor_clicked()));

    LBLnegativecontourcolor = new QLabel(tr("Negative levels")+": ");
    BTNnegativecontourcolor = new ColorButton();
    BTNnegativecontourcolor->setIcon(QIcon(":/images/colores48.png"));
    BTNnegativecontourcolor->setText(tr("Color"));
    QColor negativecontourcolor = grafica->getnegativecontourcolor();
    BTNnegativecontourcolor->setColor(&(negativecontourcolor));
    BTNnegativecontourcolor->setEnabled(true);
    negativecontourcolor = grafica->getnegativecontourcolor();
    connect(BTNnegativecontourcolor, SIGNAL(clicked()), this, SLOT(BTNnegativecontourcolor_clicked()));

    LBLzerocontourcolor = new QLabel(tr("Zero level")+": ");
    BTNzerocontourcolor = new ColorButton();
    BTNzerocontourcolor->setIcon(QIcon(":/images/colores48.png"));
    BTNzerocontourcolor->setText(tr("Color"));
    QColor zerocontourcolor = grafica->getzerocontourcolor();
    BTNzerocontourcolor->setColor(&(zerocontourcolor));
    BTNzerocontourcolor->setEnabled(true);
    zerocontourcolor = grafica->getzerocontourcolor();
    connect(BTNzerocontourcolor, SIGNAL(clicked()), this, SLOT(BTNzerocontourcolor_clicked()));

    CHKtranspbg = new QCheckBox(tr("Transparent labels background"));
    CHKtranspbg->setChecked(false);
    CHKtranspbg->setEnabled(false);
    CHKtranspbg->setHidden(true);
    grafica->settranspbackground(false);
    connect(CHKtranspbg, SIGNAL(stateChanged(int)), this, SLOT(CHKtranspbg_changed()));

//         Pen width

    FRMcntpen = new QGroupBox(FRMcontours);
    FRMcntpen->setTitle(tr("Pen"));
//    FRMcntpen->setVisible(false);
    LBLcntpenwidth = new QLabel(tr("Thickness"),FRMcntpen);
    SPBcntpenwidth = new QSpinBox(FRMcntpen);
    SPBcntpenwidth->setRange(0, 8);
    SPBcntpenwidth->setValue(2);
    SPBcntpenwidth->setMaximumWidth(40);
    grafica->setcntpenwidth(2);
    connect(SPBcntpenwidth, SIGNAL(valueChanged(int)), this, SLOT(SPBcntpenwidth_changed()));

//        Contours styles

    FRMcontourstyle = new QGroupBox(page_contours);
    FRMcontourstyle->setTitle(tr("Styles"));
//    FRMcontourstyle->setVisible(false);

    FRMcontourpostyle = new QGroupBox(tr("Positive levels"));
    RBTpossolid=new QRadioButton(tr("Solid"),FRMcontourpostyle);
    RBTposdashed=new QRadioButton(tr("Dashed"),FRMcontourpostyle);
    RBTposdotted=new QRadioButton(tr("Dotted"),FRMcontourpostyle);
    RBTposdashdot=new QRadioButton(tr("Dashed dotted"),FRMcontourpostyle);
    RBTpossolid->setChecked(true);
    RBTposdashed->setChecked(false);
    RBTposdotted->setChecked(false);
    RBTposdashdot->setChecked(false);
    grafica->setpositivecontourstyle(0);
    connect(RBTpossolid, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTposdashed, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTposdotted, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTposdashdot, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));

    FRMcontournegstyle = new QGroupBox(tr("Negative levels"));
    RBTnegsolid=new QRadioButton(tr("Solid"),FRMcontournegstyle);
    RBTnegdashed=new QRadioButton(tr("Dashed"),FRMcontournegstyle);
    RBTnegdotted=new QRadioButton(tr("Dotted"),FRMcontournegstyle);
    RBTnegdashdot=new QRadioButton(tr("Dashed dotted"),FRMcontournegstyle);
    RBTnegsolid->setChecked(false);
    RBTnegdashed->setChecked(true);
    RBTnegdotted->setChecked(false);
    RBTnegdashdot->setChecked(false);
    grafica->setnegativecontourstyle(1);
    connect(RBTnegsolid, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTnegdashed, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTnegdotted, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTnegdashdot, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));

    FRMcontourzerostyle = new QGroupBox(tr("Zero level"));
    RBTzerosolid=new QRadioButton(tr("Solid"),FRMcontourzerostyle);
    RBTzerodashed=new QRadioButton(tr("Dashed"),FRMcontourzerostyle);
    RBTzerodotted=new QRadioButton(tr("Dotted"),FRMcontourzerostyle);
    RBTzerodashdot=new QRadioButton(tr("Dashed dotted"),FRMcontourzerostyle);
    RBTzerosolid->setChecked(false);
    RBTzerodashed->setChecked(false);
    RBTzerodotted->setChecked(true);
    RBTzerodashdot->setChecked(false);
    grafica->setzerocontourstyle(2);
    connect(RBTzerosolid, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTzerodashed, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTzerodotted, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
    connect(RBTzerodashdot, SIGNAL(toggled (bool)), this, SLOT(RBTcontpostyle_changed()));
}

void Viewer2D::page_contours_layouts()
{
//              Contours layouts

    QHBoxLayout *page1Layout0=new QHBoxLayout();
    page1Layout0->addWidget(LBLcontourfile);
    page1Layout0->addWidget(TXTcontour);
    page1Layout0->addWidget(BTNcontour);

    QHBoxLayout *page1Layout1=new QHBoxLayout();
    page1Layout1->addStretch();
    page1Layout1->addWidget(BTNactualizarcontour);
    page1Layout1->addStretch();

    QHBoxLayout *page1Layout2=new QHBoxLayout();
    page1Layout2->addStretch();
    page1Layout2->addWidget(Wtablacont);
    page1Layout2->addStretch();

    QVBoxLayout *page1Layout3  = new QVBoxLayout(FRMcontours);
    page1Layout3->addLayout(page1Layout2);
    page1Layout3->addLayout(page1Layout1);
    page1Layout3->addStretch();


//            Contours label layouts

    QHBoxLayout *page1Layout3b=new QHBoxLayout();
    page1Layout3b->addWidget(CHKcontourslabel);
    page1Layout3b->setAlignment(Qt::AlignCenter);

//             Contours label font layouts

    QHBoxLayout *page1Layout3c=new QHBoxLayout();
    page1Layout3c->addStretch();
    page1Layout3c->addWidget(BTNfontcontourslabel);
    page1Layout3c->addWidget(BTNcontourslabelcolor);
    page1Layout3c->addStretch();

    QHBoxLayout *page1Layout3d=new QHBoxLayout();
    page1Layout3d->addWidget(LBLtolerance);
    page1Layout3d->addWidget(SPBtolerance);
    page1Layout3d->setAlignment(Qt::AlignCenter);

    QHBoxLayout *page1Layout3e=new QHBoxLayout();
    page1Layout3e->addWidget(LBLzerotolerance);
    page1Layout3e->addWidget(SPBzerotolerance);
    page1Layout3e->setAlignment(Qt::AlignCenter);

    QHBoxLayout *page1Layout3g=new QHBoxLayout();
    page1Layout3g->addWidget(CHKtranspbg);
    page1Layout3g->setAlignment(Qt::AlignCenter);

    QVBoxLayout *page1Layout3f = new QVBoxLayout(FRMcontourslabel);
    page1Layout3f->addLayout(page1Layout3b);
    page1Layout3f->addLayout(page1Layout3c);
    page1Layout3f->addLayout(page1Layout3d);
    page1Layout3f->addLayout(page1Layout3e);
    page1Layout3f->addLayout(page1Layout3g);

//          Contours pen with layouts

    QHBoxLayout *page1Layout4a=new QHBoxLayout(FRMcntpen);
    page1Layout4a->addStretch();
    page1Layout4a->addWidget(LBLcntpenwidth);
    page1Layout4a->addWidget(SPBcntpenwidth);
    page1Layout4a->addStretch();

//             Contours color layouts

    QHBoxLayout *page1Layout4b=new QHBoxLayout();
    page1Layout4b->addWidget(CHKmulticolor);
    page1Layout4b->setAlignment(Qt::AlignCenter);

    QGridLayout *page1Layout5=new QGridLayout();
    page1Layout5->addWidget(LBLpositivecontourcolor,0,0);
    page1Layout5->addWidget(BTNpositivecontourcolor,0,1);
    page1Layout5->addWidget(LBLzerocontourcolor,1,0);
    page1Layout5->addWidget(BTNzerocontourcolor,1,1);
    page1Layout5->addWidget(LBLnegativecontourcolor,2,0);
    page1Layout5->addWidget(BTNnegativecontourcolor,2,1);
    page1Layout5->setAlignment(Qt::AlignCenter);

    QVBoxLayout *page1Layout6  = new QVBoxLayout(FRMcontourscolor);
    page1Layout6->addLayout(page1Layout4b);
    page1Layout6->addLayout(page1Layout5);

//             Contours style layouts

    QVBoxLayout *page1Layout7  = new QVBoxLayout(FRMcontourpostyle);
    page1Layout7->addWidget(RBTpossolid);
    page1Layout7->addWidget(RBTposdashed);
    page1Layout7->addWidget(RBTposdotted);
    page1Layout7->addWidget(RBTposdashdot);

    QVBoxLayout *page1Layout8  = new QVBoxLayout(FRMcontourzerostyle);
    page1Layout8->addWidget(RBTzerosolid);
    page1Layout8->addWidget(RBTzerodashed);
    page1Layout8->addWidget(RBTzerodotted);
    page1Layout8->addWidget(RBTzerodashdot);

    QVBoxLayout *page1Layout9  = new QVBoxLayout(FRMcontournegstyle);
    page1Layout9->addWidget(RBTnegsolid);
    page1Layout9->addWidget(RBTnegdashed);
    page1Layout9->addWidget(RBTnegdotted);
    page1Layout9->addWidget(RBTnegdashdot);

    QVBoxLayout *page1Layout10  = new QVBoxLayout(FRMcontourstyle);
    page1Layout10->addWidget(FRMcontourpostyle);
    page1Layout10->addWidget(FRMcontourzerostyle);
    page1Layout10->addWidget(FRMcontournegstyle);

    QVBoxLayout *page1Layout  = new QVBoxLayout(page_contours);
    page1Layout->addLayout(page1Layout0);
    page1Layout->addWidget(FRMcontours);
    page1Layout->addWidget(FRMcontourslabel);
    page1Layout->addWidget(FRMcntpen);
    page1Layout->addWidget(FRMcontourscolor);
    page1Layout->addWidget(FRMcontourstyle);
    page1Layout->addStretch();
}

//    page_field: Field Lines
//

void Viewer2D::page_field_widgets()
{

    TXTfield = new QLineEdit();
    BTNfield = new QToolButton();
    BTNfield->setText(tr("..."));
    LBLfield = new QLabel(tr("File"));
    connect(BTNfield, SIGNAL(clicked()), this, SLOT(importfieldfile_dialog()));

    CHKshowfield = new QCheckBox(tr("Show field lines"));
    CHKshowfield->setChecked(true);
    connect(CHKshowfield, SIGNAL(stateChanged(int)), this, SLOT(CHKshowfield_changed()));

//          Arrows
    FRMarrows = new QGroupBox(page_field);
    FRMarrows->setTitle(tr("Arrows"));
    CHKarrows = new QCheckBox(tr("Show arrows"));
    CHKarrows->setChecked(false);
    connect(CHKarrows, SIGNAL(stateChanged(int)), this, SLOT(CHKarrows_changed()));

//          Arrows spacing

    SPBarrowssep = new QSpinBox();
    SPBarrowssep->setMinimum(1);
    SPBarrowssep->setMaximum(1000);
    SPBarrowssep->setMaximumWidth(90);
    SPBarrowssep->setSingleStep(10);
    SPBarrowssep->setValue(150);
    connect(SPBarrowssep,SIGNAL(valueChanged(int)),this,SLOT(SPBarrowssep_changed(int)));

//          Arrows size

    SPBarrowssize = new QSpinBox();
    SPBarrowssize->setMinimum(1);
    SPBarrowssize->setMaximum(50);
    SPBarrowssize->setMaximumWidth(90);
    SPBarrowssize->setSingleStep(5);
    SPBarrowssize->setValue(15);
    connect(SPBarrowssize,SIGNAL(valueChanged(int)),this,SLOT(SPBarrowssize_changed(int)));

//          Arrows width

    SPBarrowswidth = new QSpinBox();
    SPBarrowswidth->setMinimum(1);
    SPBarrowswidth->setMaximum(30);
    SPBarrowswidth->setMaximumWidth(90);
    SPBarrowswidth->setSingleStep(1);
    SPBarrowswidth->setValue(6);
    connect(SPBarrowswidth,SIGNAL(valueChanged(int)),this,SLOT(SPBarrowswidth_changed(int)));

//          Arrows skew

    SPBarrowsskew = new QSpinBox();
    SPBarrowsskew->setMinimum(1);
    SPBarrowsskew->setMaximum(30);
    SPBarrowsskew->setMaximumWidth(90);
    SPBarrowsskew->setSingleStep(1);
    SPBarrowsskew->setValue(6);
    connect(SPBarrowsskew,SIGNAL(valueChanged(int)),this,SLOT(SPBarrowsskew_changed(int)));

//         Pen width

    FRMelinespen = new QGroupBox(page_field);
    FRMelinespen->setTitle(tr("Pen"));
    LBLelinespenwidth = new QLabel(tr("Thickness"),FRMelinespen);
    SPBelinespenwidth = new QSpinBox(FRMelinespen);
    SPBelinespenwidth->setRange(0, 8);
    SPBelinespenwidth->setValue(2);
    SPBelinespenwidth->setMaximumWidth(60);
    SPBelinespenwidth_changed();
    connect(SPBelinespenwidth, SIGNAL(valueChanged(int)), this, SLOT(SPBelinespenwidth_changed()));

}

void Viewer2D::page_field_layouts()
{
    QLabel *LBLarrowssep = new QLabel(tr("Arrows spacing"));
    QLabel *LBLarrowssize = new QLabel(tr("Arrow Length"));
    QLabel *LBLarrowswidth = new QLabel(tr("Arrow Width"));
    QLabel *LBLarrowsskew = new QLabel(tr("Arrows Skew"));
    QHBoxLayout *page2Layout13 = new QHBoxLayout();

    QHBoxLayout *page2Layout0=new QHBoxLayout();
    page2Layout0->addWidget(LBLfield);
    page2Layout0->addWidget(TXTfield);
    page2Layout0->addWidget(BTNfield);

    QHBoxLayout *page2Layout1=new QHBoxLayout();
    page2Layout1->addStretch();
    page2Layout1->addWidget(CHKshowfield);
    page2Layout1->addStretch();

    page2Layout13->addStretch();
    page2Layout13->addWidget(CHKarrows);
    page2Layout13->addStretch();

    QGridLayout *page2Layout14=new QGridLayout();
    page2Layout14->addWidget(LBLarrowssep,0,0,Qt::AlignLeft);
    page2Layout14->addWidget(SPBarrowssep,0,1,Qt::AlignRight);
    page2Layout14->addWidget(LBLarrowssize,1,0,Qt::AlignLeft);
    page2Layout14->addWidget(SPBarrowssize,1,1,Qt::AlignRight);
    page2Layout14->addWidget(LBLarrowswidth,2,0,Qt::AlignLeft);
    page2Layout14->addWidget(SPBarrowswidth,2,1,Qt::AlignRight);
    page2Layout14->addWidget(LBLarrowsskew,3,0,Qt::AlignLeft);
    page2Layout14->addWidget(SPBarrowsskew,3,1,Qt::AlignRight);


    QVBoxLayout *page2Layout15 = new QVBoxLayout(FRMarrows);
    page2Layout15->addLayout(page2Layout13);
    page2Layout15->addLayout(page2Layout14);

    QHBoxLayout *page2Layout_16=new QHBoxLayout(FRMelinespen);
    page2Layout_16->addStretch();
    page2Layout_16->addWidget(LBLelinespenwidth);
    page2Layout_16->addWidget(SPBelinespenwidth);
    page2Layout_16->addStretch();

    QVBoxLayout *page2Layout = new QVBoxLayout(page_field);
    page2Layout->addLayout(page2Layout0);
    page2Layout->addLayout(page2Layout1);
    page2Layout->addWidget(FRMarrows);
    page2Layout->addWidget(FRMelinespen);
    page2Layout->addStretch();
}

//    page_frad: Radial factors
//

void Viewer2D::page_sghistogram_widgets()
{

    TXTsghist = new QLineEdit();
    BTNsghist = new QToolButton();
    BTNsghist->setText(tr("..."));
    LBLsghist = new QLabel(tr("File"));
    connect(BTNsghist, SIGNAL(clicked()), this, SLOT(importhstfile_dialog()));

//             Sigma histogram labels

    CHKsghistlabel = new QCheckBox(tr("Show/hide"));
    CHKsghistlabel->setChecked(true);
    CHKsghistlabel->setEnabled(true);
    connect(CHKsghistlabel, SIGNAL(stateChanged(int)), this, SLOT(CHKsghistlabel_changed()));
    FRMfontsghistlabels = new QGroupBox(page_sghistogram);
    FRMfontsghistlabels->setTitle(tr("Label")+"s");
    BTNfontsghistlabels = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connect(BTNfontsghistlabels, SIGNAL(clicked()), this, SLOT(BTNfontsghistlabels_click()));

//              Pen width

    FRMsghistpen = new QGroupBox(page_sghistogram);
    FRMsghistpen->setTitle(tr("Pen"));
    LBLsghistpenwidth = new QLabel(tr("Thickness"),FRMsghistpen);
    SPBsghistpenwidth = new QSpinBox(FRMsghistpen);
    SPBsghistpenwidth->setRange(0, 8);
    SPBsghistpenwidth->setValue(2);
    SPBsghistpenwidth->setMaximumWidth(40);
    grafica->setsghistpenwidth(SPBsghistpenwidth->value());
    connect(SPBsghistpenwidth, SIGNAL(valueChanged(int)), this, SLOT(SPBsghistpenwidth_changed()));

//              Smoothing lines

    FRMsmoothlines = new QGroupBox(page_sghistogram);
    FRMsmoothlines->setTitle(tr("Smooth lines"));
    CHKsmoothlines = new QCheckBox(tr("Set smoothing"),FRMsmoothlines);
    CHKsmoothlines->setChecked(true);
    CHKsmoothlines->setEnabled(true);
    connect(CHKsmoothlines, SIGNAL(stateChanged(int)), this, SLOT(CHKsmoothlines_changed()));
    LBLsmoothfactor = new QLabel(tr("Smooth factor"),FRMsmoothlines);
    SPBsmoothfactor = new QSpinBox(FRMsmoothlines);
    SPBsmoothfactor->setRange(0, 9);
    SPBsmoothfactor->setValue(smoothfactor);
    SPBsmoothfactor->setMaximumWidth(40);
    connect(SPBsmoothfactor, SIGNAL(valueChanged(int)), this, SLOT(SPBsmoothfactor_changed()));
    BTNsmooth = new QPushButton();
    BTNsmooth->setText(tr("Apply"));
    connect(BTNsmooth, SIGNAL(clicked()), this, SLOT(BTNsmooth_click()));

//              Partial histograms

    CHKshowparthist = new QCheckBox(tr("Show partial histograms"));
    CHKshowparthist->setChecked(true);
    CHKshowparthist->setEnabled(true);
    connect(CHKshowparthist, SIGNAL(stateChanged(int)), this, SLOT(CHKshowparthist_changed()));

//              Print histograms data to file


    FRMprinthist = new QGroupBox(page_sghistogram);

    CHKprinthist = new QCheckBox(tr("Print histograms data to file"),FRMprinthist);
    CHKprinthist->setChecked(true);
    CHKprinthist->setEnabled(true);
    connect(CHKprinthist, SIGNAL(stateChanged(int)), this, SLOT(CHKprinthist_changed()));

    BTNprinthist = new QPushButton(FRMprinthist);
    BTNprinthist->setText(tr("Print"));
    connect(BTNprinthist, SIGNAL(clicked()), this, SLOT(BTNprinthist_click()));

    TXTprinthist = new QLineEdit(FRMprinthist);
    BTNprthist = new QToolButton(FRMprinthist);
    BTNprthist->setText(tr("..."));
    LBLprthist = new QLabel(tr("File"),FRMprinthist);
    connect(BTNprthist, SIGNAL(clicked()), this, SLOT(BTNprthist_click()));

}

void Viewer2D::page_sghistogram_layouts()
{

    QHBoxLayout *Layout0=new QHBoxLayout();
    Layout0->addWidget(LBLsghist);
    Layout0->addWidget(TXTsghist);
    Layout0->addWidget(BTNsghist);

    QVBoxLayout *Layout3  = new QVBoxLayout(FRMfontsghistlabels);
    Layout3->addWidget(CHKsghistlabel);
    Layout3->addWidget(BTNfontsghistlabels);
    Layout3->setAlignment(Qt::AlignCenter);

    QHBoxLayout *Layout4=new QHBoxLayout(FRMsghistpen);
    Layout4->addStretch();
    Layout4->addWidget(LBLsghistpenwidth);
    Layout4->addWidget(SPBsghistpenwidth);
    Layout4->addStretch();

    QHBoxLayout *Layout5=new QHBoxLayout();
    Layout5->addWidget(LBLsmoothfactor);
    Layout5->addWidget(SPBsmoothfactor);
    Layout5->addStretch();

    QHBoxLayout *Layout6=new QHBoxLayout();
    Layout6->addStretch();
    Layout6->addWidget(BTNsmooth);
    Layout6->addStretch();

    QVBoxLayout *Layout7  = new QVBoxLayout(FRMsmoothlines);
    Layout7->addWidget(CHKsmoothlines);
    Layout7->addLayout(Layout5);
    Layout7->addLayout(Layout6);

    QHBoxLayout *Layout8=new QHBoxLayout();
    Layout8->addWidget(LBLprthist);
    Layout8->addWidget(TXTprinthist);
    Layout8->addWidget(BTNprthist);

    QHBoxLayout *Layout9=new QHBoxLayout();
    Layout9->addStretch();
    Layout9->addWidget(BTNprinthist);
    Layout9->addStretch();

    QVBoxLayout *Layout10=new QVBoxLayout(FRMprinthist);
    Layout10->addWidget(CHKprinthist);
    Layout10->addLayout(Layout8);
    Layout10->addLayout(Layout9);
    Layout10->addStretch();

    QVBoxLayout *Layout  = new QVBoxLayout(page_sghistogram);
    Layout->addLayout(Layout0);
    Layout->addWidget(FRMfontsghistlabels);
    Layout->addWidget(FRMsghistpen);
    Layout->addWidget(FRMsmoothlines);
    Layout->addWidget(CHKshowparthist);
    Layout->addWidget(FRMprinthist);
    Layout->addStretch();
}

//    page_frad: Radial factors
//

void Viewer2D::page_frad_widgets()
{

//              Type of radial factor

    TXTfrads = new QLineEdit();
    BTNfrads = new QToolButton();
    BTNfrads->setText(tr("..."));
    LBLfrads = new QLabel(tr("File"));
    connect(BTNfrads, SIGNAL(clicked()), this, SLOT(importfradfile_dialog()));

    FRMradfactortype = new QGroupBox(page_frad);
    QString string = QString("f")+QChar(0x2097)+QChar(0x2098)+QString("(r)");
    RBTradfactor=new QRadioButton(FRMradfactortype);
    LBLradfactor = new QLabel("f<sub>lm</sub>(r)");
    LBLradfactor->setFont(QFont("Helvetica",15));
    RBTradfactor->setChecked(true);
    RBTr2flm=new QRadioButton(FRMradfactortype);
    LBLr2flm = new QLabel("r<sup>2l </sup>f<sub>lm</sub>(r)");
    LBLr2flm->setFont(QFont("Helvetica",15));
    RBTr2flm->setChecked(false);
    RBTr2l2flm=new QRadioButton(FRMradfactortype);
    LBLr2l2flm = new QLabel("r<sup>2l+2 </sup>f<sub>lm</sub>(r)");
    LBLr2l2flm->setFont(QFont("Helvetica",15));
    RBTr2l2flm->setChecked(false);
    connect(RBTradfactor, SIGNAL(toggled (bool)), this, SLOT(RBTfrad_changed()));
    connect(RBTr2flm, SIGNAL(toggled (bool)), this, SLOT(RBTfrad_changed()));
    connect(RBTr2l2flm, SIGNAL(toggled (bool)), this, SLOT(RBTfrad_changed()));

//             Radial factors labels

    CHKcurvelabel = new QCheckBox(tr("Show/hide"));
    CHKcurvelabel->setChecked(true);
    CHKcurvelabel->setEnabled(true);
    connect(CHKcurvelabel, SIGNAL(stateChanged(int)), this, SLOT(CHKcurvelabel_changed()));
    FRMfontcurvelabels = new QGroupBox(page_frad);
    FRMfontcurvelabels->setTitle(tr("Label")+"s");
    BTNfontcurvelabels = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connect(BTNfontcurvelabels, SIGNAL(clicked()), this, SLOT(BTNfontcurvelabels_click()));

//              Pen width

    FRMcurvespen = new QGroupBox(page_frad);
    FRMcurvespen->setTitle(tr("Pen"));
    LBLcurvespenwidth = new QLabel(tr("Thickness"),FRMcurvespen);
    SPBcurvespenwidth = new QSpinBox(FRMcurvespen);
    SPBcurvespenwidth->setRange(0, 8);
    SPBcurvespenwidth->setValue(2);
    SPBcurvespenwidth->setMaximumWidth(40);
    grafica->setcntpenwidth(2);
    connect(SPBcurvespenwidth, SIGNAL(valueChanged(int)), this, SLOT(SPBcurvespenwidth_changed()));
}

void Viewer2D::page_frad_layouts()
{
//              Layouts of type of radial factor

    QHBoxLayout *page3Layout0=new QHBoxLayout();
    page3Layout0->addWidget(LBLfrads);
    page3Layout0->addWidget(TXTfrads);
    page3Layout0->addWidget(BTNfrads);

    QGridLayout *page3Layout1  = new QGridLayout();
    page3Layout1->addWidget(LBLradfactor,0,1);
    page3Layout1->addWidget(RBTradfactor,0,0);
    page3Layout1->addWidget(LBLr2flm,1,1);
    page3Layout1->addWidget(RBTr2flm,1,0);
    page3Layout1->addWidget(LBLr2l2flm,2,1);
    page3Layout1->addWidget(RBTr2l2flm,2,0);

    QHBoxLayout *page3Layout2=new QHBoxLayout(FRMradfactortype);
    page3Layout2->addLayout(page3Layout1);
    page3Layout2->setAlignment(Qt::AlignLeft);

    QVBoxLayout *page3Layout3  = new QVBoxLayout(FRMfontcurvelabels);
    page3Layout3->addWidget(CHKcurvelabel);
    page3Layout3->addWidget(BTNfontcurvelabels);
    page3Layout3->setAlignment(Qt::AlignCenter);

    QHBoxLayout *page3Layout4=new QHBoxLayout(FRMcurvespen);
    page3Layout4->addStretch();
    page3Layout4->addWidget(LBLcurvespenwidth);
    page3Layout4->addWidget(SPBcurvespenwidth);
    page3Layout4->addStretch();

    QVBoxLayout *page3Layout  = new QVBoxLayout(page_frad);
    page3Layout->addLayout(page3Layout0);
    page3Layout->addWidget(FRMradfactortype);
    page3Layout->addWidget(FRMfontcurvelabels);
    page3Layout->addWidget(FRMcurvespen);
    page3Layout->addStretch();
}

//    page_CPs: Critical Points
//

void Viewer2D::page_CPs_widgets()
{
//           CPs

    FRMcps= new QGroupBox(page_CPs);
    FRMcps->setTitle(tr("Critical points"));
    TXTcps = new QLineEdit();
    BTNcps = new QToolButton();
    BTNcps->setText(tr("..."));
    LBLcpsfile = new QLabel(tr("File"));
    connect(BTNcps, SIGNAL(clicked()), this, SLOT(importcpsfile_dialog()));

    cpscolor << QColor(255,0,0) << QColor(0,255,0) << QColor(195,195,195) << QColor(255,128,0);
    grafica->setcpscolor(cpscolor);
    for (int i=0;i<max_cps;i++){
        LBLcps[i] = new QLabel();
        CHKcps[i] = new QCheckBox();
        CHKcps[i]->setChecked(true);
        connect(CHKcps[i], SIGNAL(stateChanged(int)), this, SLOT(CHKcps_changed()));
//        QRgb colorRgb = cpscolor[i].rgba();
        BTNcolorcps[i] = new ColorButton();
        BTNcolorcps[i]->setIcon(QIcon(":/images/colores48.png"));
        BTNcolorcps[i]->setText(tr("Color"));
        BTNcolorcps[i]->setColor(&cpscolor[i]);
        connect(BTNcolorcps[i],SIGNAL(clicked()),this,SLOT(BTNcolorcps_change()));
    }

    BTNcps->setToolTip(tr("Open critial points file ..."));
    LBLcps[0]->setText(tr("(3,+3) CP")+":");
    LBLcps[1]->setText(tr("(3,+1) CP")+":");
    LBLcps[2]->setText(tr("(3,-1) CP")+":");
    LBLcps[3]->setText(tr("(3,-3) CP")+":");

//         Ball radius

    LBLcpsballradius = new QLabel(tr("Ball radius"));
    SPBcpsballradius = new QSpinBox();
    SPBcpsballradius->setMinimum(1);
    SPBcpsballradius->setMaximum(20);
    SPBcpsballradius->setMaximumWidth(50);
    SPBcpsballradius->setValue(3);
    grafica->setcenterradius(SPBcpsballradius->value());
    connect(SPBcpsballradius,SIGNAL(valueChanged(int)),this,SLOT(SPBcpsballradius_changed(int)));

//          Distance to plane tolerance

    TXTdistancethreshold = new QLineEdit(tr("1.e-3"));
    TXTdistancethreshold->setMaximumWidth(50);
    TXTdistancethreshold->setValidator(myDoubleValidator);
    TXTdistancethreshold->setAlignment(Qt::AlignRight);
    connect(TXTdistancethreshold, SIGNAL(textChanged(const QString &)), this, SLOT(TXTdistancethreshold_changed()));
}

void Viewer2D::page_CPs_layouts()
{
    QLabel *LBLdistancethreshold = new QLabel(tr("Distance to plane tolerance")+":");
//             Critical points layouts

    QHBoxLayout *page4Layout1=new QHBoxLayout();
    page4Layout1->addWidget(LBLcpsfile);
    page4Layout1->addWidget(TXTcps);
    page4Layout1->addWidget(BTNcps);

    QGridLayout *page4Layout2=new QGridLayout();
    page4Layout2->addWidget(LBLcps[0],0,0,Qt::AlignLeft);
    page4Layout2->addWidget(CHKcps[0],0,1,Qt::AlignCenter);
    page4Layout2->addWidget(BTNcolorcps[0],0,2,Qt::AlignRight);
    page4Layout2->addWidget(LBLcps[1],1,0,Qt::AlignLeft);
    page4Layout2->addWidget(CHKcps[1],1,1,Qt::AlignCenter);
    page4Layout2->addWidget(BTNcolorcps[1],1,2,Qt::AlignRight);
    page4Layout2->addWidget(LBLcps[2],2,0,Qt::AlignLeft);
    page4Layout2->addWidget(CHKcps[2],2,1,Qt::AlignCenter);
    page4Layout2->addWidget(BTNcolorcps[2],2,2,Qt::AlignRight);
    page4Layout2->addWidget(LBLcps[3],3,0,Qt::AlignLeft);
    page4Layout2->addWidget(CHKcps[3],3,1,Qt::AlignCenter);
    page4Layout2->addWidget(BTNcolorcps[3],3,2,Qt::AlignRight);

//             Ball radius

    QGridLayout *page4Layout3=new QGridLayout();
    page4Layout3->addWidget(LBLcpsballradius,0,0,Qt::AlignLeft);
    page4Layout3->addWidget(SPBcpsballradius,0,1,Qt::AlignRight);
    page4Layout3->addWidget(LBLdistancethreshold,1,0,Qt::AlignLeft);
    page4Layout3->addWidget(TXTdistancethreshold,1,1,Qt::AlignRight);


//              CPs layout

    QVBoxLayout *page4Layout7=new QVBoxLayout(FRMcps);
    page4Layout7->addLayout(page4Layout2);
    page4Layout7->addLayout(page4Layout3);
    page4Layout7->addStretch();

//         Page 4 layout

    QVBoxLayout *page4Layout=new QVBoxLayout(page_CPs);
    page4Layout->addLayout(page4Layout1);
    page4Layout->addWidget(FRMcps);
}

//    page_basins: Basins
//

void Viewer2D::page_basins_widgets()
{
    TXTbasins = new QLineEdit();
    BTNbasins = new QToolButton();
    BTNbasins->setText(tr("..."));
    LBLbasins = new QLabel(tr("File"));

    connect(BTNbasins, SIGNAL(clicked()), this, SLOT(importbasinsfile_dialog()));

    CHKbasins = new QCheckBox(tr("Show basins"));
    CHKbasins->setChecked(false);
    connect(CHKbasins, SIGNAL(stateChanged(int)), this, SLOT(CHKbasins_changed()));

//         Pen width

    FRMbasinspen = new QGroupBox(page_basins);
    FRMbasinspen->setTitle(tr("Pen"));
    LBLbasinspenwidth = new QLabel(tr("Thickness"),FRMbasinspen);
    SPBbasinspenwidth = new QSpinBox(FRMbasinspen);
    SPBbasinspenwidth->setRange(0, 8);
    SPBbasinspenwidth->setValue(3);
    SPBbasinspenwidth->setMaximumWidth(40);
    grafica->setbasinspenwidth(3);
    connect(SPBbasinspenwidth, SIGNAL(valueChanged(int)), this, SLOT(SPBbasinspenwidth_changed()));

//         Basins frontier color

    BTNbasinscolor = new ColorButton();
    BTNbasinscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNbasinscolor->setText(tr("Color"));
    BTNbasinscolor->setColor(basinscolor);
    connect(BTNbasinscolor, SIGNAL(clicked()), this, SLOT(BTNbasinscolor_clicked()));
    BTNbasinscolor->setEnabled(true);
}

void Viewer2D::page_basins_layouts()
{
//             Basins layouts

    QHBoxLayout *page5Layout1=new QHBoxLayout();
    page5Layout1->addWidget(LBLbasins);
    page5Layout1->addWidget(TXTbasins);
    page5Layout1->addWidget(BTNbasins);

    QHBoxLayout *page5Layout2=new QHBoxLayout();
    page5Layout2->addStretch();
    page5Layout2->addWidget(CHKbasins);
    page5Layout2->addStretch();

    QHBoxLayout *page5Layout3=new QHBoxLayout(FRMbasinspen);
    page5Layout3->addStretch();
    page5Layout3->addWidget(LBLbasinspenwidth);
    page5Layout3->addWidget(SPBbasinspenwidth);
    page5Layout3->addWidget(BTNbasinscolor);
    page5Layout3->addStretch();

    QVBoxLayout *page5Layout=new QVBoxLayout(page_basins);
    page5Layout->addLayout(page5Layout1);
    page5Layout->addLayout(page5Layout2);
    page5Layout->addWidget(FRMbasinspen);
    page5Layout->addStretch();
}

//    page_options: Options
//

void Viewer2D::page_options_widgets()
{
//      Title
    title = tr("Title");
    grafica->title = title;
    QFontMetrics fm( grafica->getfonttitle() );
    QSize fmsize = fm.size( Qt::TextSingleLine, title );
    grafica->titlerect = QRect(0,0, fmsize.width(), fmsize.height() );
    grafica->titlerect.moveCenter(QPoint(grafica->width()/2,fmsize.height()));
    grafica->setTopMargin(max(50,2*fmsize.height()+20));

    FRMtitlegroup = new QGroupBox(page_options);
    FRMtitlegroup->setTitle(tr("Title"));
    CHKtitle = new QCheckBox(tr("Show/hide")+" "+tr("label"));
    CHKtitle->setChecked(true);
    CHKtitle->setEnabled(true);
        grafica->setshowtitle(true);
    connect(CHKtitle, SIGNAL(stateChanged(int)), this, SLOT(CHKtitle_changed()));

//      Title Font
    BTNedittitle = new QPushButton();
    BTNedittitle->setText(tr("Edit"));
    connect(BTNedittitle, SIGNAL(clicked()), this, SLOT(editar_title()));
    connect(grafica, SIGNAL(title_dialog()), this, SLOT(editar_title()));

//      Position title
    LBLpostitle = new QLabel(tr("Position")+": ");
    LBLxtitle = new QLabel("X");
    LBLytitle = new QLabel("Y");
    TXTxtitle = new QLineEdit(QString("%1").arg(grafica->titlerect.center().x()));
    TXTxtitle->setValidator(myDoubleValidator);
    TXTxtitle->setMaximumWidth(40);
    TXTxtitle->setAlignment(Qt::AlignCenter);
    TXTytitle = new QLineEdit(QString("%1").arg(grafica->titlerect.center().y()));
    TXTytitle->setValidator(myDoubleValidator);
    TXTytitle->setMaximumWidth(40);
    TXTytitle->setAlignment(Qt::AlignCenter);
    connect(TXTxtitle, SIGNAL(returnPressed()), this, SLOT(TXTxtitle_changed()));
    connect(TXTytitle, SIGNAL(returnPressed()), this, SLOT(TXTytitle_changed()));
    connect(grafica,SIGNAL(titleposChanged()),this, SLOT(TXTtitlepos_changed()));
    CHKtitlebox = new QCheckBox(tr("Show/hide")+" "+tr("Box"));
    CHKtitlebox->setChecked(false);
    CHKtitlebox->setEnabled(true);
    connect(CHKtitlebox, SIGNAL(stateChanged(int)), this, SLOT(CHKtitlebox_changed()));

//      X axis label
    Xlabel = tr("u axis");
    grafica->Xlabel = Xlabel;

    fm = QFontMetrics( grafica->getfontXlabel() );
    fmsize = fm.size( Qt::TextSingleLine, Xlabel );
    grafica->Xlabelrect = QRect(0,0, fmsize.width(), fmsize.height() );
    grafica->Xlabelrect.moveCenter(QPoint(grafica->width()/2,grafica->height()-fmsize.height()));
    grafica->setBottomMargin(max(50,2*fmsize.height()+20));
    FRMXlabelgroup = new QGroupBox(page_options);
    FRMXlabelgroup->setTitle(tr("X axis"));
    CHKXlabel = new QCheckBox(tr("Show/hide")+" "+tr("label"));
    CHKXlabel->setChecked(true);
    CHKXlabel->setEnabled(true);
    connect(CHKXlabel, SIGNAL(stateChanged(int)), this, SLOT(CHKXlabel_changed()));

//            X axis label Font
    BTNeditXlabel = new QPushButton();
    BTNeditXlabel->setText(tr("Edit"));
    connect(BTNeditXlabel, SIGNAL(clicked()), this, SLOT(editar_Xlabel()));
    connect(grafica, SIGNAL(Xlabel_dialog()), this, SLOT(editar_Xlabel()));

//        Position X axis label
    LBLposXlabel = new QLabel(tr("Position")+": ");
    LBLxXlabel = new QLabel("X");
    LBLyXlabel = new QLabel("Y");
    TXTxXlabel = new QLineEdit(QString("%1").arg(grafica->Xlabelrect.center().x()));
    TXTxXlabel->setValidator(myDoubleValidator);
    TXTxXlabel->setMaximumWidth(40);
    TXTxXlabel->setAlignment(Qt::AlignCenter);
    TXTyXlabel = new QLineEdit(QString("%1").arg(grafica->Xlabelrect.center().y()));
    TXTyXlabel->setValidator(myDoubleValidator);
    TXTyXlabel->setMaximumWidth(40);
    TXTyXlabel->setAlignment(Qt::AlignCenter);
    connect(TXTxXlabel, SIGNAL(returnPressed()), this, SLOT(TXTxXlabel_changed()));
    connect(TXTyXlabel, SIGNAL(returnPressed()), this, SLOT(TXTyXlabel_changed()));
    connect(grafica,SIGNAL(XlabelposChanged()),this, SLOT(TXTXlabelpos_changed()));
    CHKXlabelbox = new QCheckBox(tr("Show/hide")+" "+tr("Box"));
    CHKXlabelbox->setChecked(false);
    CHKXlabelbox->setEnabled(true);
    connect(CHKXlabelbox, SIGNAL(stateChanged(int)), this, SLOT(CHKXlabelbox_changed()));

//        Y axis label
    Ylabel = tr("v axis");
    grafica->Ylabel = Ylabel;

    fm = QFontMetrics( grafica->getfontYlabel() );
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    grafica->Ylabelrect = QRect(0,0, fmsize.height(), fmsize.width() );
    grafica->Ylabelrect.moveCenter(QPoint(fmsize.height(),grafica->height()/2));
    grafica->setYlabelvert(true);
    grafica->setLeftMargin(max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    grafica->setRightMargin(max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    FRMYlabelgroup = new QGroupBox(page_options);
    FRMYlabelgroup->setTitle(tr("Y axis"));
    CHKYlabel = new QCheckBox(tr("Show/hide")+" "+tr("label"));
    CHKYlabel->setChecked(true);
    CHKYlabel->setEnabled(true);
    grafica->setshowXlabel(true);
    connect(CHKYlabel, SIGNAL(stateChanged(int)), this, SLOT(CHKYlabel_changed()));

//            Y axis label Font
    BTNeditYlabel = new QPushButton();
    BTNeditYlabel->setText(tr("Edit"));
    connect(BTNeditYlabel, SIGNAL(clicked()), this, SLOT(editar_Ylabel()));
    connect(grafica, SIGNAL(Ylabel_dialog()), this, SLOT(editar_Ylabel()));

//        Position Y axis label
    LBLposYlabel = new QLabel(tr("Position")+": ");
    LBLxYlabel = new QLabel("X");
    LBLyYlabel = new QLabel("Y");
    TXTxYlabel = new QLineEdit(QString("%1").arg(grafica->Ylabelrect.center().x()));
    TXTxYlabel->setValidator(myDoubleValidator);
    TXTxYlabel->setMaximumWidth(40);
    TXTxYlabel->setAlignment(Qt::AlignCenter);
    TXTyYlabel = new QLineEdit(QString("%1").arg(grafica->Ylabelrect.center().y()));
    TXTyYlabel->setValidator(myDoubleValidator);
    TXTyYlabel->setMaximumWidth(40);
    TXTyYlabel->setAlignment(Qt::AlignCenter);
    connect(TXTxYlabel, SIGNAL(returnPressed()), this, SLOT(TXTxYlabel_changed()));
    connect(TXTyYlabel, SIGNAL(returnPressed()), this, SLOT(TXTyYlabel_changed()));
    connect(grafica,SIGNAL(YlabelposChanged()),this, SLOT(TXTYlabelpos_changed()));
    CHKYlabelhorizontal = new QCheckBox(tr("Horizontal"));
    LBLYlabelhorizontal = new QLabel(tr("Orientation")+":");
    CHKYlabelhorizontal->setChecked(false);
    CHKYlabelhorizontal->setEnabled(true);
    connect(CHKYlabelhorizontal, SIGNAL(stateChanged(int)), this, SLOT(CHKYlabelhorizontal_changed()));
    CHKYlabelbox = new QCheckBox(tr("Show/hide")+" "+tr("Box"));
    CHKYlabelbox->setChecked(false);
    CHKYlabelbox->setEnabled(true);
    grafica->setshowYlabel(true);
    connect(CHKYlabelbox, SIGNAL(stateChanged(int)), this, SLOT(CHKYlabelbox_changed()));

//        Aspect ratio
    FRMaspectratio = new QGroupBox(page_options);
    FRMaspectratio->setEnabled(true);
    FRMaspectratio->setTitle(tr("Aspect ratio"));
    CHKaspectratio = new QCheckBox(tr("Keep aspect ratio"));
    CHKaspectratio->setChecked(true);
    connect(CHKaspectratio, SIGNAL(stateChanged(int)), this, SLOT(CHKaspectratio_changed()));

//        Zoom region
    FRMzoomregion = new QGroupBox(page_options);
    FRMzoomregion->setTitle(tr("Zoom region"));
    LBLx=new QLabel(tr("x"),FRMzoomregion);
    LBLy=new QLabel(tr("y"),FRMzoomregion);
    LBLinf=new QLabel(tr("Lowest"),FRMzoomregion);
    LBLsup=new QLabel(tr("Highest"),FRMzoomregion);
    TXTxinf=new QLineEdit(FRMzoomregion);
    TXTxinf->setText(QString::number(grafica->getxmin()));
    TXTxinf->setAlignment(Qt::AlignRight);
    TXTxinf->setValidator(myDoubleValidator);
    TXTxsup=new QLineEdit(FRMzoomregion);
    TXTxsup->setText(QString::number(grafica->getxmax()));
    TXTxsup->setAlignment(Qt::AlignRight);
    TXTxsup->setValidator(myDoubleValidator);
    TXTyinf=new QLineEdit(FRMzoomregion);
    TXTyinf->setText(QString::number(grafica->getymin()));
    TXTyinf->setAlignment(Qt::AlignRight);
    TXTyinf->setValidator(myDoubleValidator);
    TXTysup=new QLineEdit(FRMzoomregion);
    TXTysup->setText(QString::number(grafica->getymax()));
    TXTysup->setAlignment(Qt::AlignRight);
    TXTysup->setValidator(myDoubleValidator);
    BTNexeczoomreg=new QPushButton(QIcon(":/images/ejecutar.png"), tr("Apply"),FRMzoomregion);
    BTNexeczoomreg->setToolTip(tr("Generate grids"));
    connect(BTNexeczoomreg, SIGNAL(clicked()), this, SLOT(execzoomreg()));
    connect(grafica,SIGNAL(zoomChanged()),this, SLOT(zoom_changed()));
    connect(grafica,SIGNAL(moveToTop()),this, SLOT(emit_moveToTopPlotter()));

//        Margins
    FRMmargins = new QGroupBox(page_options);
    FRMmargins->setTitle(tr("Margins"));

    LBLBottomMargin = new QLabel(tr("Bottom"));
    SPBBottomMargin = new QSpinBox();
    SPBBottomMargin->setMinimum(0);
    SPBBottomMargin->setMaximum(400);
    SPBBottomMargin->setMaximumWidth(90);
    SPBBottomMargin->setValue(int(5 + 2.5*QFontMetrics(grafica->getfontXlabel()).height()));
    SPBBottomMargin->setKeyboardTracking(false);  // To prevent updating with each key pressed
    connect(SPBBottomMargin,SIGNAL(valueChanged(int)),this,SLOT(SPBBottomMargin_changed(int)));

    LBLLeftMargin = new QLabel(tr("Left"));
    SPBLeftMargin = new QSpinBox();
    SPBLeftMargin->setMinimum(0);
    SPBLeftMargin->setMaximum(400);
    SPBLeftMargin->setMaximumWidth(90);
    SPBLeftMargin->setValue(grafica->getLeftMargin());
    SPBLeftMargin->setKeyboardTracking(false);  // To prevent updating with each key pressed
    connect(SPBLeftMargin,SIGNAL(valueChanged(int)),this,SLOT(SPBLeftMargin_changed(int)));
    connect(grafica,SIGNAL(Left_Margin_changed(int)),this,SLOT(SPBLeftMargin_changed(int)));

    LBLRightMargin = new QLabel(tr("Right"));
    SPBRightMargin = new QSpinBox();
    SPBRightMargin->setMinimum(0);
    SPBRightMargin->setMaximum(400);
    SPBRightMargin->setMaximumWidth(90);
    SPBRightMargin->setValue(grafica->getRightMargin());
    SPBRightMargin->setKeyboardTracking(false);  // To prevent updating with each key pressed
    connect(SPBRightMargin,SIGNAL(valueChanged(int)),this,SLOT(SPBRightMargin_changed(int)));

    LBLTopMargin = new QLabel(tr("Top"));
    SPBTopMargin = new QSpinBox();
    SPBTopMargin->setMinimum(0);
    SPBTopMargin->setMaximum(400);
    SPBTopMargin->setMaximumWidth(90);
    SPBTopMargin->setValue(grafica->getTopMargin());
    SPBTopMargin->setKeyboardTracking(false);
    connect(SPBTopMargin,SIGNAL(valueChanged(int)),this,SLOT(SPBTopMargin_changed(int)));

    connect(grafica,SIGNAL(MarginsChanged(int, int, int, int)),this,SLOT(set_SPBmargins(int, int, int, int)));


//         Background color

    FRMbkgcolor = new QGroupBox(tr("Background color"),page_options);
    BTNbkgcolor = new ColorButton();
    BTNbkgcolor->setIcon(QIcon(":/images/colores48.png"));
    BTNbkgcolor->setText(tr("Color"));
    BTNbkgcolor->setColor(bkgcolor);
    connect(BTNbkgcolor, SIGNAL(clicked()), this, SLOT(BTNbkgcolor_clicked()));
    BTNbkgcolor->setEnabled(true);

//         Grid
    FRMgridcolor = new QGroupBox(tr("Grid"),page_options);
    gridcolor = new QColor(0,0,255);
    BTNgridcolor = new ColorButton();
    BTNgridcolor->setIcon(QIcon(":/images/colores48.png"));
    BTNgridcolor->setText(tr("Color"));
    BTNgridcolor->setColor(gridcolor);
    grafica->setgridcolor(*gridcolor);
    connect(BTNgridcolor, SIGNAL(clicked()), this, SLOT(BTNgridcolor_clicked()));
    BTNgridcolor->setEnabled(true);
    CHKgrid = new QCheckBox(tr("Show grid"));
    CHKgrid->setChecked(true);
    CHKgrid->setEnabled(true);
    grafica->setshowgrid(true);
    connect(CHKgrid, SIGNAL(stateChanged(int)), this, SLOT(CHKgrid_changed()));

    FRMgridstyle = new QGroupBox();
    RBTsolid = new QRadioButton(tr("Solid"),FRMgridstyle);
    RBTsolid->setChecked(true);
    grafica->setsolidgrid(true);
    RBTdashed = new QRadioButton(tr("Dashed"),FRMgridstyle);
    connect(RBTsolid, SIGNAL(toggled (bool)), this, SLOT(RBTsolid_changed()));
    connect(RBTdashed, SIGNAL(toggled (bool)), this, SLOT(RBTsolid_changed()));

//         Scales
    FRMscales = new QGroupBox(tr("Scales"),page_options);
    CHKXscalebottom = new QCheckBox(tr("Bottom"));
    CHKXscalebottom->setChecked(true);
    CHKXscalebottom->setEnabled(true);
    grafica->setshowXscalebottom(true);
    connect(CHKXscalebottom, SIGNAL(stateChanged(int)), this, SLOT(CHKXscale_changed()));
    CHKXscaletop = new QCheckBox(tr("Top"));
    CHKXscaletop->setChecked(false);
    CHKXscaletop->setEnabled(true);
    grafica->setshowXscaletop(false);
    connect(CHKXscaletop, SIGNAL(stateChanged(int)), this, SLOT(CHKXscale_changed()));
    CHKYscaleleft = new QCheckBox(tr("Left"));
    CHKYscaleleft->setChecked(true);
    CHKYscaleleft->setEnabled(true);
    grafica->setshowYscaleleft(true);
    connect(CHKYscaleleft, SIGNAL(stateChanged(int)), this, SLOT(CHKYscale_changed()));
    CHKYscaleright = new QCheckBox(tr("Right"));
    CHKYscaleright->setChecked(false);
    CHKYscaleright->setEnabled(true);
    grafica->setshowYscaleright(false);
    connect(CHKYscaleright, SIGNAL(stateChanged(int)), this, SLOT(CHKYscale_changed()));
    CHKticks = new QCheckBox(tr("Show ticks"));
    CHKticks->setChecked(true);
    CHKticks->setEnabled(true);
    grafica->setshowXticksbottom(true);
    grafica->setshowXtickstop(false);
    grafica->setshowYticksleft(true);
    grafica->setshowYticksright(false);
    connect(CHKticks, SIGNAL(stateChanged(int)), this, SLOT(CHKXscale_changed()));
    connect(CHKticks, SIGNAL(stateChanged(int)), this, SLOT(CHKYscale_changed()));

//        Ticks
    FRMticks = new QGroupBox(page_options);
    FRMticks->setTitle(tr("Ticks labels"));

    CHKautomaticticks = new QCheckBox(tr("Automatic"));
    CHKautomaticticks->setChecked(true);
    CHKautomaticticks->setEnabled(true);
    connect(CHKautomaticticks, SIGNAL(stateChanged(int)), this, SLOT(CHKautomaticticks_changed()));

    LBLXticks = new QLabel(tr("X ticks"));
    LBLXticks->setEnabled(false);
    SPBXticks = new QSpinBox();
    SPBXticks->setMinimum(2);
    SPBXticks->setMaximum(20);
    SPBXticks->setMaximumWidth(50);
    SPBXticks->setValue(5);
    SPBXticks->setEnabled(false);
    connect(SPBXticks,SIGNAL(valueChanged(int)),grafica,SLOT(Xticks_changed(int)));

    LBLXcifras = new QLabel(tr("Decimal figures"));
    LBLXcifras->setEnabled(true);
    SPBXcifras = new QSpinBox();
    SPBXcifras->setMinimum(0);
    SPBXcifras->setMaximum(6);
    SPBXcifras->setMaximumWidth(50);
    SPBXcifras->setValue(1);
    SPBXcifras->setEnabled(true);
    connect(SPBXcifras,SIGNAL(valueChanged(int)),this,SLOT(SPBXcifras_changed(int)));

    LBLYticks = new QLabel(tr("Y ticks"));
    LBLYticks->setEnabled(false);
    SPBYticks = new QSpinBox();
    SPBYticks->setMinimum(0);
    SPBYticks->setMaximum(20);
    SPBYticks->setMaximumWidth(50);
    SPBYticks->setValue(5);
    SPBYticks->setEnabled(false);
    connect(SPBYticks,SIGNAL(valueChanged(int)),grafica,SLOT(Yticks_changed(int)));

    LBLYcifras = new QLabel(tr("Decimal figures"));
    LBLYcifras->setEnabled(true);
    SPBYcifras = new QSpinBox();
    SPBYcifras->setMinimum(0);
    SPBYcifras->setMaximum(6);
    SPBYcifras->setMaximumWidth(50);
    SPBYcifras->setValue(1);
    SPBYcifras->setEnabled(true);
    connect(SPBYcifras,SIGNAL(valueChanged(int)),this,SLOT(SPBYcifras_changed(int)));
    grafica->setXcifras(SPBXcifras->value());
    grafica->setYcifras(SPBYcifras->value());

//          Centers and bonds
    FRMcenters = new QGroupBox(page_options);
    FRMcenters->setTitle(tr("Centers"));
    FRMcenters->setEnabled(false);

//          Centers

    CHKcenters = new QCheckBox(tr("Show centers"));
    CHKcenters->setChecked(true);
    connect(CHKcenters, SIGNAL(stateChanged(int)), this, SLOT(CHKcenters_changed()));

//          Centers radius

    SPBcenterradius = new QSpinBox();
    SPBcenterradius->setRange(0, 100);
    SPBcenterradius->setValue(3);
    SPBcenterradius->setMaximumWidth(40);
    connect(SPBcenterradius, SIGNAL(valueChanged(int)), this, SLOT(SPBcenterradius_changed()));

//            Centers color


    BTNcenterscolor = new ColorButton();
    BTNcenterscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNcenterscolor->setText(tr("Color"));
    BTNcenterscolor->setColor(&centerscolor);
    BTNcenterscolor->setEnabled(true);
    connect(BTNcenterscolor, SIGNAL(clicked()), this, SLOT(BTNcenterscolor_clicked()));

//          Centers labels

    CHKcenterslabel = new QCheckBox(tr("Show centers labels"));
    CHKcenterslabel->setChecked(true);
    connect(CHKcenterslabel, SIGNAL(stateChanged(int)), this, SLOT(CHKcenterslabel_changed()));

//            Centers labels font

    BTNfontcenterslabel = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));

    connect(BTNfontcenterslabel, SIGNAL(clicked()), this, SLOT(BTNfontcenterslabel_click()));

//            Centers labels color

    BTNcenterslabelcolor = new ColorButton();
    BTNcenterslabelcolor->setIcon(QIcon(":/images/fonts48.png"));
    BTNcenterslabelcolor->setText(tr("Color"));
    BTNcenterslabelcolor->setColor(&centerslabelcolor);
    BTNcenterslabelcolor->setEnabled(true);
    connect(BTNcenterslabelcolor, SIGNAL(clicked()), this, SLOT(BTNcenterslabelcolor_clicked()));

    SPBshftcenterlabels = new QSpinBox();
    SPBshftcenterlabels->setRange(-100, 100);
    SPBshftcenterlabels->setValue(0);
    SPBshftcenterlabels->setMaximumWidth(40);
    connect(SPBshftcenterlabels, SIGNAL(valueChanged(int)), this, SLOT(SPBshftcenterlabels_changed(int)));


//          Bonds

    FRMbonds = new QGroupBox(page_options);
    FRMbonds->setTitle(tr("Bonds"));
    FRMbonds->setEnabled(false);

    CHKbonds = new QCheckBox(tr("Show bonds"));
    CHKbonds->setChecked(true);
    connect(CHKbonds, SIGNAL(stateChanged(int)), this, SLOT(CHKbonds_changed()));

//          Bonds thickness

    SPBbondswidth = new QSpinBox();
    SPBbondswidth->setRange(0, 8);
    SPBbondswidth->setValue(2);
    SPBbondswidth->setMaximumWidth(40);
    connect(SPBbondswidth, SIGNAL(valueChanged(int)), this, SLOT(SPBbondswidth_changed(int)));

//          Bonds lenght tolerance

    TXTbondthreshold = new QLineEdit(tr("1.2"));
    TXTbondthreshold->setMaximumWidth(50);
    TXTbondthreshold->setValidator(myDoubleValidator);
    TXTbondthreshold->setAlignment(Qt::AlignRight);
    connect(TXTbondthreshold, SIGNAL(textChanged(const QString &)), this, SLOT(TXTbondthreshold_changed()));

//            Bonds color
    BTNbondscolor = new ColorButton();
    BTNbondscolor->setIcon(QIcon(":/images/colores48.png"));
    BTNbondscolor->setText(tr("Color"));
    BTNbondscolor->setColor(&bondscolor);
    BTNbondscolor->setEnabled(true);
    connect(BTNbondscolor, SIGNAL(clicked()), this, SLOT(BTNbondscolor_clicked()));

}

void Viewer2D::page_options_layouts()
{
//      Title layouts
    QHBoxLayout *page6Layout1=new QHBoxLayout();
    page6Layout1->addWidget(CHKtitle);
    page6Layout1->setAlignment(Qt::AlignCenter);

//      Title font layouts
    QHBoxLayout *page6Layout2=new QHBoxLayout();
    page6Layout2->addStretch();
    page6Layout2->addWidget(BTNedittitle);
    page6Layout2->addStretch();

//      Position title layouts
    QHBoxLayout *page6Layout3=new QHBoxLayout();
    page6Layout3->addStretch();
    page6Layout3->addWidget(LBLpostitle);
    page6Layout3->addStretch();

    QHBoxLayout *page6Layout4=new QHBoxLayout();
    page6Layout4->addStretch();
    page6Layout4->addWidget(LBLxtitle);
    page6Layout4->addWidget(TXTxtitle);
    page6Layout4->addStretch();

    QHBoxLayout *page6Layout5=new QHBoxLayout();
    page6Layout5->addStretch();
    page6Layout5->addWidget(LBLytitle);
    page6Layout5->addWidget(TXTytitle);
    page6Layout5->addStretch();

    QHBoxLayout *page6Layout6=new QHBoxLayout();
    page6Layout6->addLayout(page6Layout3);
    page6Layout6->addLayout(page6Layout4);
    page6Layout6->addLayout(page6Layout5);

    QHBoxLayout *page6Layout_7=new QHBoxLayout();
    page6Layout_7->addStretch();
    page6Layout_7->addWidget(CHKtitlebox);
    page6Layout_7->addStretch();

    QVBoxLayout *page6Layout_8 = new QVBoxLayout(FRMtitlegroup);
    page6Layout_8->addLayout(page6Layout1);
    page6Layout_8->addLayout(page6Layout2);
    page6Layout_8->addLayout(page6Layout6);
    page6Layout_8->addLayout(page6Layout_7);

//      X axis label layouts
    QHBoxLayout *page6Layout_9=new QHBoxLayout();
    page6Layout_9->addWidget(CHKXlabel);
    page6Layout_9->setAlignment(Qt::AlignCenter);

//      X axis label font layouts
    QHBoxLayout *page6Layout_10=new QHBoxLayout();
    page6Layout_10->addStretch();
    page6Layout_10->addWidget(BTNeditXlabel);
    page6Layout_10->addStretch();

//      Position X axis label layouts
    QHBoxLayout *page6Layout_11=new QHBoxLayout();
    page6Layout_11->addStretch();
    page6Layout_11->addWidget(LBLposXlabel);
    page6Layout_11->addStretch();

    QHBoxLayout *page6Layout_12=new QHBoxLayout();
    page6Layout_12->addStretch();
    page6Layout_12->addWidget(LBLxXlabel);
    page6Layout_12->addWidget(TXTxXlabel);
    page6Layout_12->addStretch();

    QHBoxLayout *page6Layout_13=new QHBoxLayout();
    page6Layout_13->addStretch();
    page6Layout_13->addWidget(LBLyXlabel);
    page6Layout_13->addWidget(TXTyXlabel);
    page6Layout_13->addStretch();

    QHBoxLayout *page6Layout_14=new QHBoxLayout();
    page6Layout_14->addLayout(page6Layout_11);
    page6Layout_14->addLayout(page6Layout_12);
    page6Layout_14->addLayout(page6Layout_13);

    QHBoxLayout *page6Layout_15=new QHBoxLayout();
    page6Layout_15->addStretch();
    page6Layout_15->addWidget(CHKXlabelbox);
    page6Layout_15->addStretch();

    QVBoxLayout *page6Layout_16 = new QVBoxLayout(FRMXlabelgroup);
    page6Layout_16->addLayout(page6Layout_9);
    page6Layout_16->addLayout(page6Layout_10);
    page6Layout_16->addLayout(page6Layout_14);
    page6Layout_16->addLayout(page6Layout_15);

//      Y axis label layouts
    QHBoxLayout *page6Layout_17=new QHBoxLayout();
    page6Layout_17->addWidget(CHKYlabel);
    page6Layout_17->setAlignment(Qt::AlignCenter);

//      Y axis label font layouts
    QHBoxLayout *page6Layout_18=new QHBoxLayout();
    page6Layout_18->addStretch();
    page6Layout_18->addWidget(BTNeditYlabel);
    page6Layout_18->addStretch();

//      Position Y axis label layouts
    QHBoxLayout *page6Layout_19=new QHBoxLayout();
    page6Layout_19->addStretch();
    page6Layout_19->addWidget(LBLposYlabel);
    page6Layout_19->addStretch();

    QHBoxLayout *page6Layout_20=new QHBoxLayout();
    page6Layout_20->addStretch();
    page6Layout_20->addWidget(LBLxYlabel);
    page6Layout_20->addWidget(TXTxYlabel);
    page6Layout_20->addStretch();

    QHBoxLayout *page6Layout_21=new QHBoxLayout();
    page6Layout_21->addStretch();
    page6Layout_21->addWidget(LBLyYlabel);
    page6Layout_21->addWidget(TXTyYlabel);
    page6Layout_21->addStretch();

    QHBoxLayout *page6Layout_22=new QHBoxLayout();
    page6Layout_22->addLayout(page6Layout_19);
    page6Layout_22->addLayout(page6Layout_20);
    page6Layout_22->addLayout(page6Layout_21);

    QHBoxLayout *page6Layout_23=new QHBoxLayout();
    page6Layout_23->addStretch();
    page6Layout_23->addWidget(CHKYlabelbox);
    page6Layout_23->addStretch();

    QHBoxLayout *page6Layout_24=new QHBoxLayout();
    page6Layout_24->addStretch();
    page6Layout_24->addWidget(LBLYlabelhorizontal);
    page6Layout_24->addWidget(CHKYlabelhorizontal);
    page6Layout_24->addStretch();

    QVBoxLayout *page6Layout_25 = new QVBoxLayout(FRMYlabelgroup);
    page6Layout_25->addLayout(page6Layout_17);
    page6Layout_25->addLayout(page6Layout_18);
    page6Layout_25->addLayout(page6Layout_22);
    page6Layout_25->addLayout(page6Layout_23);
    page6Layout_25->addLayout(page6Layout_24);

//      Aspect ratio layouts
    QHBoxLayout *page6Layout_26=new QHBoxLayout();
    page6Layout_26->addWidget(CHKaspectratio);
    QVBoxLayout *page6Layout_27  = new QVBoxLayout(FRMaspectratio);
    page6Layout_27->addLayout(page6Layout_26);

//        Zoom region layouts
    QGridLayout *page6Layout_28=new QGridLayout();
    page6Layout_28->addWidget(LBLinf,0,1);
    page6Layout_28->addWidget(LBLsup,0,2);
    page6Layout_28->addWidget(LBLx,1,0);
    page6Layout_28->addWidget(TXTxinf,1,1);
    page6Layout_28->addWidget(TXTxsup,1,2);
    page6Layout_28->addWidget(LBLy,2,0);
    page6Layout_28->addWidget(TXTyinf,2,1);
    page6Layout_28->addWidget(TXTysup,2,2);

    QHBoxLayout *page6Layout_29=new QHBoxLayout();
    page6Layout_29->addStretch();
    page6Layout_29->addWidget(BTNexeczoomreg);
    page6Layout_29->addStretch();

    QVBoxLayout *page6Layout_30 = new QVBoxLayout(FRMzoomregion);
    page6Layout_30->addLayout(page6Layout_28);
    page6Layout_30->addLayout(page6Layout_29);

//      Margins layouts
    QHBoxLayout *page6Layout_31=new QHBoxLayout();
    page6Layout_31->addWidget(SPBBottomMargin);
    page6Layout_31->addWidget(LBLBottomMargin);
    page6Layout_31->addWidget(SPBTopMargin);
    page6Layout_31->addWidget(LBLTopMargin);

    QHBoxLayout *page6Layout_32=new QHBoxLayout();
    page6Layout_32->addWidget(SPBLeftMargin);
    page6Layout_32->addWidget(LBLLeftMargin);
    page6Layout_32->addWidget(SPBRightMargin);
    page6Layout_32->addWidget(LBLRightMargin);

    QVBoxLayout *page6Layout_33 = new QVBoxLayout(FRMmargins);
    page6Layout_33->addLayout(page6Layout_31);
    page6Layout_33->addLayout(page6Layout_32);

//      Background color layouts

    QHBoxLayout *page6Layout_35=new QHBoxLayout();
    page6Layout_35->addStretch();
    page6Layout_35->addWidget(BTNbkgcolor);
    page6Layout_35->addStretch();


    QVBoxLayout *page6Layout_36=new QVBoxLayout(FRMbkgcolor);
    page6Layout_36->addLayout(page6Layout_35);
    page6Layout_36->addStretch();

//      Grid layouts
    QHBoxLayout *page6Layout_37=new QHBoxLayout();
    page6Layout_37->addWidget(CHKgrid);
    page6Layout_37->setAlignment(Qt::AlignCenter);


    QVBoxLayout *page6Layout_38=new QVBoxLayout();
    page6Layout_38->addWidget(RBTsolid);
    page6Layout_38->addWidget(RBTdashed);
    page6Layout_38->setAlignment(Qt::AlignCenter);

    QHBoxLayout *page6Layout_39=new QHBoxLayout();
    page6Layout_39->addStretch();
    page6Layout_39->addWidget(BTNgridcolor);
    page6Layout_39->addLayout(page6Layout_38);
    page6Layout_39->addStretch();

    QVBoxLayout *page6Layout_40 = new QVBoxLayout(FRMgridcolor);
    page6Layout_40->addLayout(page6Layout_37);
    page6Layout_40->addLayout(page6Layout_39);

//      Scales layouts
    QGridLayout *page6Layout_41=new QGridLayout();
    page6Layout_41->addWidget(CHKXscalebottom, 0,0);
    page6Layout_41->addWidget(CHKXscaletop, 0,1);
    page6Layout_41->addWidget(CHKYscaleleft, 1,0);
    page6Layout_41->addWidget(CHKYscaleright, 1,1);

    QHBoxLayout *page6Layout_42=new QHBoxLayout();
    page6Layout_42->addWidget(CHKticks);
    page6Layout_42->setAlignment(Qt::AlignCenter);

    QVBoxLayout *page6Layout_43 = new QVBoxLayout(FRMscales);
    page6Layout_43->addLayout(page6Layout_41);
    page6Layout_43->addLayout(page6Layout_42);

//      Ticks layouts
    QHBoxLayout *page6Layout_44=new QHBoxLayout();
    page6Layout_44->addWidget(CHKautomaticticks);
    page6Layout_44->setAlignment(Qt::AlignCenter);

    QHBoxLayout *page6Layout_45=new QHBoxLayout();
    page6Layout_45->addWidget(SPBXticks);
    page6Layout_45->addWidget(LBLXticks);
    page6Layout_45->addWidget(SPBXcifras);
    page6Layout_45->addWidget(LBLXcifras);

    QHBoxLayout *page6Layout_46=new QHBoxLayout();
    page6Layout_46->addWidget(SPBYticks);
    page6Layout_46->addWidget(LBLYticks);
    page6Layout_46->addWidget(SPBYcifras);
    page6Layout_46->addWidget(LBLYcifras);

    QVBoxLayout *page6Layout_47 = new QVBoxLayout(FRMticks);
    page6Layout_47->addLayout(page6Layout_44);
    page6Layout_47->addLayout(page6Layout_45);
    page6Layout_47->addLayout(page6Layout_46);
    page6Layout_47->setAlignment(Qt::AlignCenter);

//      Centers layouts
    QHBoxLayout *page6Layout48 = new QHBoxLayout();
    page6Layout48->addStretch();
    page6Layout48->addWidget(CHKcenters);
    page6Layout48->addStretch();

    QLabel *LBLcenterradius = new QLabel(tr("Center thickness"));
    QHBoxLayout *page6Layout49 = new QHBoxLayout();
    page6Layout49->addStretch();
    page6Layout49->addWidget(LBLcenterradius);
    page6Layout49->addWidget(SPBcenterradius);
    page6Layout49->addWidget(BTNcenterscolor);
    page6Layout49->addStretch();

    QHBoxLayout *page6Layout50 = new QHBoxLayout();
    page6Layout50->addStretch();
    page6Layout50->addWidget(CHKcenterslabel);
    page6Layout50->addStretch();

    QHBoxLayout *page6Layout51 = new QHBoxLayout();
    page6Layout51->addStretch();
    page6Layout51->addWidget(BTNfontcenterslabel);
    page6Layout51->addWidget(BTNcenterslabelcolor);
    page6Layout51->addStretch();

    QLabel *LBLshftcenterlabels = new QLabel(tr("Labels vertical displacement"));
    QHBoxLayout *page6Layout52 = new QHBoxLayout();
    page6Layout52->addStretch();
    page6Layout52->addWidget(LBLshftcenterlabels);
    page6Layout52->addWidget(SPBshftcenterlabels);
    page6Layout52->addStretch();

    QVBoxLayout *page6Layout53 = new QVBoxLayout(FRMcenters);
    page6Layout53->addLayout(page6Layout48);
    page6Layout53->addLayout(page6Layout49);
    page6Layout53->addLayout(page6Layout50);
    page6Layout53->addLayout(page6Layout51);
    page6Layout53->addLayout(page6Layout52);

//      Bonds layouts
    QHBoxLayout *page6Layout54 = new QHBoxLayout();
    page6Layout54->addStretch();
    page6Layout54->addWidget(CHKbonds);
    page6Layout54->addStretch();

    QLabel *LBLbondthreshold = new QLabel(tr("Bonding tolerance")+":");
    QLabel *LBLbondswidth = new QLabel(tr("Bonds thickness"));
    QGridLayout *page6Layout55 = new QGridLayout();
    page6Layout55->addWidget(LBLbondthreshold,0,0);
    page6Layout55->addWidget(TXTbondthreshold,0,1);
    page6Layout55->addWidget(LBLbondswidth,1,0);
    page6Layout55->addWidget(SPBbondswidth,1,1);

    QHBoxLayout *page6Layout56 = new QHBoxLayout();
    page6Layout56->addStretch();
    page6Layout56->addWidget(BTNbondscolor);
    page6Layout56->addStretch();

    QVBoxLayout *page6Layout57 = new QVBoxLayout(FRMbonds);
    page6Layout57->addLayout(page6Layout54);
    page6Layout57->addLayout(page6Layout55);
    page6Layout57->addLayout(page6Layout56);

//      Options Layouts
    QVBoxLayout *page6Layout58 = new QVBoxLayout(page_options);
    page6Layout58->addWidget(FRMtitlegroup);
    page6Layout58->addWidget(FRMXlabelgroup);
    page6Layout58->addWidget(FRMYlabelgroup);
    page6Layout58->addWidget(FRMaspectratio);
    page6Layout58->addWidget(FRMzoomregion);
    page6Layout58->addWidget(FRMmargins);
    page6Layout58->addWidget(FRMbkgcolor);
    page6Layout58->addWidget(FRMgridcolor);
    page6Layout58->addWidget(FRMscales);
    page6Layout58->addWidget(FRMticks);
    page6Layout58->addWidget(FRMcenters);
    page6Layout58->addWidget(FRMbonds);
    page6Layout58->addStretch();
}

//    page_capture: Image Capture
//

void Viewer2D::page_capture_widgets()
{
    FRMcaptura= new QGroupBox(tr("Image capture"),page_capture);
    BTNcaptura = new QPushButton(QIcon(":/images/capturar.png"),tr("Capture"));
    connect(BTNcaptura, SIGNAL(clicked()), this, SLOT(BTNcapture_click()));

    CHKtranspbgcapture = new QCheckBox(tr("Transparent background"));
    CHKtranspbgcapture->setChecked(true);
    grafica->settranspbckgrcapture(true);
    connect(CHKtranspbgcapture, SIGNAL(stateChanged(int)), this, SLOT(CHKtranspbgcapture_changed()));

    SPBimagequality = new QSpinBox();
    SPBimagequality->setMinimum(-1);
    SPBimagequality->setMaximum(100);
    SPBimagequality->setValue(20);
    SPBimagequality->setSingleStep(5);

    TXTscalesize = new QLineEdit(tr("1.0"));
    TXTscalesize->setValidator(myDoubleValidator);
    TXTscalesize->setAlignment(Qt::AlignRight);
    TXTscalesize->setVisible(false);

    TXThsize = new QLineEdit(tr("1000"));
    TXThsize->setValidator(new QIntValidator(this));
    TXThsize->setAlignment(Qt::AlignRight);
    TXThsize->setVisible(false);
    TXTvsize = new QLineEdit(tr("1000"));
    TXTvsize->setValidator(new QIntValidator(this));
    TXTvsize->setAlignment(Qt::AlignRight);
    TXTvsize->setVisible(false);
    LBLpor = new QLabel(tr("x"));
    LBLpor->setVisible(false);

    LBLscaledef = new QLabel(tr("Scale factor"));
    LBLscaledef->setVisible(false);
    RBTscreendef = new QRadioButton();
    RBTscreendef->setChecked(true);
    RBTscreendef->setText(tr("Screen resolution"));
    connect(RBTscreendef, SIGNAL(toggled (bool)), this, SLOT(RBTscreendef_changed()));

    RBTscaledef = new QRadioButton();
    RBTscaledef->setChecked(false);
    RBTscaledef->setText(tr("Scaled resolution"));
    connect(RBTscaledef, SIGNAL(toggled (bool)), this, SLOT(RBTscreendef_changed()));

    RBTuserdef = new QRadioButton();
    RBTuserdef->setChecked(false);
    RBTuserdef->setText(tr("User defined"));
    connect(RBTuserdef, SIGNAL(toggled (bool)), this, SLOT(RBTscreendef_changed()));
}

void Viewer2D::page_capture_layouts()
{

    QLabel *LBLimagequality= new QLabel(tr("Image quality"));

    QGroupBox *FRMresolution= new QGroupBox(tr("Resolution"),FRMcaptura);

    QHBoxLayout *page7Layout1 = new QHBoxLayout();
    page7Layout1->addWidget(LBLimagequality);
    page7Layout1->addWidget(SPBimagequality);
    page7Layout1->addStretch();

    QVBoxLayout *page7Layout2 = new QVBoxLayout();
    page7Layout2->addStretch();
    page7Layout2->addWidget(RBTscreendef);
    page7Layout2->addWidget(RBTscaledef);
    page7Layout2->addWidget(RBTuserdef);
    page7Layout2->addStretch();

    QHBoxLayout *page7Layout3 = new QHBoxLayout();
    page7Layout3->addStretch();
    page7Layout3->addWidget(LBLscaledef);
    page7Layout3->addWidget(TXTscalesize);
    page7Layout3->addStretch();

    QHBoxLayout *page7Layout4 = new QHBoxLayout();
    page7Layout4->addStretch();
    page7Layout4->addWidget(TXThsize);
    page7Layout4->addWidget(LBLpor);
    page7Layout4->addWidget(TXTvsize);
    page7Layout4->addStretch();

    QVBoxLayout *page7Layout5=new QVBoxLayout(FRMresolution);
    page7Layout5->addLayout(page7Layout2);
    page7Layout5->addLayout(page7Layout3);
    page7Layout5->addLayout(page7Layout4);
    page7Layout5->addStretch();

    QVBoxLayout *page7Layout6=new QVBoxLayout(FRMcaptura);
    page7Layout6->addWidget(FRMresolution);
    page7Layout6->addWidget(CHKtranspbgcapture);
    page7Layout6->addWidget(BTNcaptura);
    page7Layout6->setAlignment(Qt::AlignCenter);
    page7Layout6->addStretch();

    QVBoxLayout *page7Layout = new QVBoxLayout(page_capture);
    page7Layout->addWidget(FRMcaptura);
    page7Layout->setAlignment(Qt::AlignCenter);
    page7Layout->addStretch();
}

//    page_save: Save/Retrieve settings
//

void Viewer2D::page_save_widgets()
{
    FRMsettings= new QGroupBox(tr("Settings"),page_save);
    BTNsaveSettings = new QPushButton();
    BTNsaveSettings->setText(tr("Save"));
    connect(BTNsaveSettings, SIGNAL(clicked()), this, SLOT(BTNsaveSettings_click()));
    BTNretrieveSettings = new QPushButton();
    BTNretrieveSettings->setText(tr("Retrieve"));
    connect(BTNretrieveSettings, SIGNAL(clicked()), this, SLOT(BTNretrieveSettings_click()));
    CHKsettingsfile = new QCheckBox(tr("Files"));
    CHKsettingsfile->setChecked(true);
    CHKsettingstitle = new QCheckBox(tr("Title"));
    CHKsettingstitle->setChecked(true);
    CHKsettingsXlabel = new QCheckBox(tr("X axis"));
    CHKsettingsXlabel->setChecked(true);
    CHKsettingsYlabel = new QCheckBox(tr("Y axis"));
    CHKsettingsYlabel->setChecked(true);
    CHKsettingsZoom = new QCheckBox(tr("Zoom"));
    CHKsettingsZoom->setChecked(true);
    CHKsettingsMargins = new QCheckBox(tr("Margins"));
    CHKsettingsMargins->setChecked(true);
    CHKsettingscntPen = new QCheckBox(tr("Curves or contours pen"));
    CHKsettingscntPen->setChecked(true);
    CHKsettingselinesPen = new QCheckBox(tr("Field lines pen"));
    CHKsettingselinesPen->setChecked(true);
    CHKsettingsbkg = new QCheckBox(tr("Background color"));
    CHKsettingsbkg->setChecked(true);
}

void Viewer2D::page_save_layouts()
{
    QGridLayout *page8Layout0=new QGridLayout();
    page8Layout0->addWidget(CHKsettingsfile, 0,0);
    page8Layout0->addWidget(CHKsettingstitle, 0,1);
    page8Layout0->addWidget(CHKsettingsXlabel, 1,0);
    page8Layout0->addWidget(CHKsettingsYlabel, 1,1);
    page8Layout0->addWidget(CHKsettingsZoom, 2,0);
    page8Layout0->addWidget(CHKsettingsMargins, 2,1);
    page8Layout0->addWidget(CHKsettingscntPen, 3,0);
    page8Layout0->addWidget(CHKsettingselinesPen, 3,1);
    page8Layout0->addWidget(CHKsettingsbkg, 4,0);

    QHBoxLayout *page8Layout1=new QHBoxLayout();
    page8Layout1->addStretch();
    page8Layout1->addWidget(BTNsaveSettings);
    page8Layout1->addWidget(BTNretrieveSettings);
    page8Layout1->addStretch();
    page8Layout1->setAlignment(Qt::AlignCenter);

    QVBoxLayout *page8Layout2=new QVBoxLayout(FRMsettings);
    page8Layout2->addLayout(page8Layout0);
    page8Layout2->addLayout(page8Layout1);
    page8Layout2->setAlignment(Qt::AlignCenter);

    QVBoxLayout *page8Layout3 = new QVBoxLayout(page_save);
    page8Layout3->addWidget(FRMsettings);
    page8Layout3->setAlignment(Qt::AlignCenter);
    page8Layout3->addStretch();
}

//***************************************************************************
//*********************  IMPORT FILES FUNCTIONS   ***************************
//***************************************************************************




//    Imports a file with field lines
void Viewer2D::importfieldfile_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder,
        tr("Import data from")+" *.cam2D *.dengr2D;;"+
        tr("Electric field files")+" (*.cam2D);;"+
        tr("Density gradient files")+" (*.dengr2D);;"+
        tr("All files")+" (*)");
    if (fileName.length()==0) return;
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    TXTfield->setText(fileName);
    if (readfieldfile(fileName)){
        plot_cam2D();
    }
}

//    Imports a file with a sigma hole histogram
void Viewer2D::importhstfile_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open sigma hole histrogram file ..."),ProjectFolder, tr("Radial factors files")
                + " *.hst " + " (*.hst);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    TXTsghist->setText(fileName);
    printhistfile = ProjectFolder+"/"+QFileInfo(TXTsghist->text()).completeBaseName()+".hst_txt";
    TXTprinthist->setText(printhistfile);
    if (readhstfile(fileName)){
        plot_sghist();
    }
}

//***************************************************************************
//**************  READ FILES AND PLOTTING FUNCTIONS    **********************
//***************************************************************************




//***************************************************************************
//***************************  CONTOUR FUNCTIONS   **************************
//***************************************************************************


bool  Viewer2D::checkplanesuffix(QString qstr){
    QStringList qstrlst = (QStringList() << "_XY0" << "_X0Z"
                << "_0YZ" << "_AB0" << "_A0C" << "_0BC" << "_ABC" );
    return qstrlst.contains(qstr);
}

void Viewer2D::BTNactualizarcontour_click(){
    int rows = SHTclist->tabla->rowCount();
    int cols = SHTclist->tabla->columnCount();
    if (levels != nullpointer)
        delete(levels);
    levels = new QList<double>;
    QStringList strlist;
    for (int i=0 ;i < rows ; ++i){
        for (int j = 0 ; j < cols ; ++j){
            if (!SHTclist->getcellvalue(i,j).isEmpty()){
                if (SHTclist->getcellvalue(i,j).toDouble() != 0.0)
                    strlist << SHTclist->getcellvalue(i,j);
                else
                    strlist << QString("%1").arg(zero);
            }
        }
    }
    strlist.removeDuplicates();
    for (int i=0 ;i < strlist.count() ; ++i){
        levels->append(strlist.at(i).toDouble());
    }
    for (int i=0 ;i < rows ; ++i){
        for (int j = 0 ; j < cols ; ++j){
            SHTclist->setcellvalue("",i,j);
        }
    }
//    qSort(levels->begin(),levels->end());
    std::sort(levels->begin(),levels->end());
    int knt = 0;
    for (int i = 0; knt < levels->count() && i < rows ; ++i){
        for (int j = 0 ; knt < levels->count() && j < cols ; ++j){
            if (levels->at(knt) != zero)
                SHTclist->setcellfloatvalue(QString::number(levels->at(knt)),i,j);
            else
                SHTclist->setcellfloatvalue(QString("0."),i,j);
            knt++;
        }
    }
    plot_cnt();
}

//    Button for choosing fonts for contour labels
void Viewer2D::BTNfontcontourslabel_click()
{
    bool OK;
    grafica->setfontcontour(QFontDialog::getFont(&OK, grafica->getfontcontour()));
    grafica->repaint();
    emit moveToTop(viewernumber);
}

//    Connects push button with colors for contour labels
void Viewer2D::BTNcontourslabelcolor_clicked()
{
    QColor col = QColorDialog::getColor(contourslabelcolor, this);
    if(col.isValid()) {
        contourslabelcolor = col;
        grafica->setcontourslabelcolor(contourslabelcolor);
        BTNcontourslabelcolor->setColor(&col);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

//    Connects push button with negative contour line color
void Viewer2D::BTNnegativecontourcolor_clicked()
{
    QColor col = QColorDialog::getColor(negativecontourcolor, this);
    if(col.isValid()) {
        negativecontourcolor = col;
        grafica->setnegativecontourcolor(negativecontourcolor);
        BTNnegativecontourcolor->setColor(&col);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

//    Connects push button with positive contour line color
void Viewer2D::BTNpositivecontourcolor_clicked()
{
    QColor col = QColorDialog::getColor(positivecontourcolor, this);
    if(col.isValid()) {
        positivecontourcolor = col;
        grafica->setpositivecontourcolor(positivecontourcolor);
        BTNpositivecontourcolor->setColor(&col);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}


//    Connects push button with negative contour line color
void Viewer2D::BTNzerocontourcolor_clicked()
{
    QColor col = QColorDialog::getColor(zerocontourcolor, this);
    if(col.isValid()) {
        zerocontourcolor = col;
        grafica->setzerocontourcolor(zerocontourcolor);
        BTNzerocontourcolor->setColor(&col);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKbasins_changed(){
    grafica->setplotbasins(CHKbasins->isChecked());
    if (grafica){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::CHKcontourslabel_changed(){
    if (CHKcontourslabel->isChecked()){
        grafica->setshowcontourlabels(true);
    }
    else{
        grafica->setshowcontourlabels(false);
    }
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKmulticolor_changed()
{
    if (CHKmulticolor->isChecked()){
        grafica->setmulticolor(true);
        LBLpositivecontourcolor->setVisible(false);
        LBLnegativecontourcolor->setVisible(false);
        LBLzerocontourcolor->setVisible(false);
        BTNpositivecontourcolor->setVisible(false);
        BTNnegativecontourcolor->setVisible(false);
        BTNzerocontourcolor->setVisible(false);
    }
    else{
        grafica->setmulticolor(false);
        LBLpositivecontourcolor->setVisible(true);
        LBLnegativecontourcolor->setVisible(true);
        LBLzerocontourcolor->setVisible(true);
        BTNpositivecontourcolor->setVisible(true);
        BTNnegativecontourcolor->setVisible(true);
        BTNzerocontourcolor->setVisible(true);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKaspectratio_changed(){
    grafica->setaspectratio(CHKaspectratio->isChecked());
    SPBLeftMargin->setDisabled(CHKaspectratio->isChecked());
    SPBRightMargin->setDisabled(CHKaspectratio->isChecked());
    SPBBottomMargin->setDisabled(CHKaspectratio->isChecked());
    SPBTopMargin->setDisabled(CHKaspectratio->isChecked());
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

//    Imports a file with CPs Cartesian coordinates
void Viewer2D::importcntfile_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open contour lines file ..."),ProjectFolder, tr("Contour lines files")
                + " .cnt " + " (*.cnt);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    TXTcontour->setText(fileName);
    if (readcntfile(fileName)){
        plot_cnt();
    }
}

//  Plots contour lines
//
void Viewer2D::plot_cnt()
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    grafica->contourMap.clear();
    ToMap *m=new ToMap;
    m->set_lower_bound(SPoint(0,0));
    m->set_upper_bound(SPoint(nu-1,nv-1));
    CContourMap *map=new CContourMap;
    map->n_levels = levels->count();
    map->levels = new double[map->n_levels];
    grafica->contourlevels->clear();
    grafica->setnlevels(map->n_levels);
    for (int i = 0 ; i < map->n_levels ; ++i){
        map->levels[i] = levels->at(i);
        if(map->levels[i] <= zmaxcnt && map->levels[i] >= zmincnt){
            grafica->contourlevels->append(map->levels[i]);
        }
    }
    int kntlines = 0;
    QVector<QPointF> *data = nullpointer;
    QVector<int> contours;
    if (grafica->contourlevels->count() > 0){
        grafica->setzero(zero);
        m->surf2d = grafica->surf2d;
        m->ncurves = 0;
        map->contour(m);
        map->consolidate();
        vector<CContourLevel*>::iterator it=map->contour_level->begin();
//      Counts the number of lines to be plotted


        while(it!=map->contour_level->end())
        {
           if(*it) {
               vector<CContour*>::iterator clit=(*it)->contour_lines->begin();
               while(clit!=(*it)->contour_lines->end())
               {
                    kntlines++;
                    clit++;
               }
           }
           it++;
        }
        data=new QVector<QPointF>[kntlines];

        it=map->contour_level->begin();
        int kntline = 0;
        int kntcontour = 0;
        grafica->curvecolors->clear();
        double uscale = (umaxcnt-umincnt) / (nu-1);
        double vscale = (vmaxcnt-vmincnt) / (nv-1);
        while(it!=map->contour_level->end())
        {
           if(*it) {
               vector<CContour*>::iterator clit=(*it)->contour_lines->begin();
               while(clit!=(*it)->contour_lines->end())
               {
                    vector<SVector>::iterator cit=(*clit)->contour->begin();
                    SPoint p=(*clit)->_start;
                    data[kntline].append(QPointF(umincnt+uscale*p.x,vmincnt+vscale*p.y));
                    while(cit!=(*clit)->contour->end())
                    {
                       p.x+=(*cit).dx;
                       p.y+=(*cit).dy;
                       data[kntline].append(QPointF(umincnt+uscale*p.x,
                               vmincnt+vscale*p.y));
                       cit++;
                    }
                    contours << kntcontour;
                    kntline++; clit++;
                }
               kntcontour++;
           }
           it++;
        }
    }
    grafica->setumincnt(umincnt);
    grafica->setumaxcnt(umaxcnt);
    grafica->setvmincnt(vmincnt);
    grafica->setvmaxcnt(vmaxcnt);
    if (superimpose){
        grafica->setumin(std::min(umincnt,grafica->getumin()));
        grafica->setumax(std::max(umaxcnt,grafica->getumax()));
        grafica->setvmin(std::min(vmincnt,grafica->getvmin()));
        grafica->setvmax(std::max(vmaxcnt,grafica->getvmax()));
        grafica->setwhoisabove(2);
    }
    else{
        grafica->setumin(umincnt);
        grafica->setumax(umaxcnt);
        grafica->setvmin(vmincnt);
        grafica->setvmax(vmaxcnt);
    }
    QFontMetrics fm( grafica->getfontYlabel() );
    QSize fmsize;
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    grafica->setLeftMargin(max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    grafica->setRightMargin(grafica->getLeftMargin());
    fm = QFontMetrics( grafica->getfonttitle() );
    fmsize = fm.size( Qt::TextSingleLine, title );
    grafica->setTopMargin(max(50,2*fmsize.height()+20));
    grafica->setBottomMargin(max(50,2*fm.size( Qt::TextSingleLine, "0" ).height()+20));
    grafica->Ylabelrect.moveCenter(QPoint(grafica->getLeftMargin() - 2*grafica->Ylabelrect.width() - 20, grafica->height()/2));
    if (grafica->getaspectratio()){
        double w = grafica->size().rwidth();
        double h = grafica->size().rheight();
        double deltax = (w - grafica->getLeftMargin() - grafica->getRightMargin())
            / (grafica->getumax()-grafica->getumin());
        double deltay = (h - grafica->getTopMargin() - grafica->getBottomMargin())
            / (grafica->getvmax()-grafica->getvmin());
        if (deltax > deltay){
            grafica->setLeftMargin((int) (0.5*(w - deltay *
                    (grafica->getumax() - grafica->getumin()))));
            grafica->setRightMargin(grafica->getLeftMargin());
        }
        else if (deltax < deltay){
            grafica->setTopMargin((int) (0.5*(h - deltax *
                    (grafica->getvmax() - grafica->getvmin()))));
            grafica->setBottomMargin(grafica->getTopMargin());
        }
        grafica->Ylabelrect.moveCenter(QPoint(grafica->getLeftMargin() - 2*grafica->Ylabelrect.width() - 20, grafica->height()/2));
    }
    PlotSettings settings;
    settings.minX = grafica->getumin();
    settings.maxX = grafica->getumax();
    settings.minY = grafica->getvmin();
    settings.maxY = grafica->getvmax();
    settings.BottomMargin = grafica->getBottomMargin();
    settings.TopMargin = grafica->getTopMargin();
    settings.LeftMargin = grafica->getLeftMargin();
    settings.RightMargin = grafica->getRightMargin();
    settings.set_numXticks(SPBXticks->value());
    settings.set_numYticks(SPBYticks->value());
    settings.adjust();
    grafica->setxmin(settings.minX);
    grafica->setxmax(settings.maxX);
    grafica->setymin(settings.minY);
    grafica->setymax(settings.maxY);
    minx = settings.minX;
    maxx = settings.maxX;
    miny = settings.minY;
    maxy = settings.maxY;
    TXTxinf->setText(QString::number(grafica->getxmin()));
    TXTxsup->setText(QString::number(grafica->getxmax()));
    TXTyinf->setText(QString::number(grafica->getymin()));
    TXTysup->setText(QString::number(grafica->getymax()));
    grafica->curvecolors->clear();
    grafica->curvecontnum->clear();
    grafica->setPlotSettings(settings);
    for (int i = 0, j = 0 ; i < kntlines ; i++ ){
        grafica->curvecontnum->append(contours[i]);
        grafica->setContourData(i, data[i]);
        if (i==0 || contours[i] != contours[i-1]){
            grafica->curvecolors->append(colorForIds[j%13]);
            j++;
        }
    }
    grafica->updatecontourlabels();
    grafica->setcenterlabelrect();
    grafica->refreshPixmap();
    grafica->updateallcenterslabelsrect();
//    Sets the plot type to contourplot plot, and optionally to field (superimpose), disables curve and sigma hole histogram plot
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    grafica->setcontourplot(true);
    if (!superimpose) grafica->setfieldplot(false);
    grafica->refreshPixmap();
    grafica->update();
    QApplication::restoreOverrideCursor();
    emit moveToTop(viewernumber);
}

//  Reads file with grid for countour lines plotting
//
bool Viewer2D::readcntfile(QString filename)
{
    QFile file(filename);
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf != "cnt"){
        QMessageBox::warning(this, tr("readcntfile"),tr("Invalid extension of file %1. Must be .cnt")
                .arg(filename));
        return false;
    }

    if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                    .arg(filename)+QString(":\n%1.").arg(file.errorString()));
            return false;
    }
    int  iaux[2];
    FILE *f;
    QByteArray ba = filename.toLatin1();
    const char *archivo = ba.data();
    f = fopen(archivo , "rb" );
    fread( iaux , sizeof(int) , 2 , f); //reads nu, nv (in this order)
    nu=-1;
    nv=-1;
    nu=iaux[0]; nv=iaux[1];
    if ( nu < 0 || nv < 0 ) {
            QMessageBox::warning(this, tr("Viewer")+"2D",tr("Error: wrong dimensions"));
            return false;
    }
    if (grafica->getfieldplot()){
        if (!retrieve){
            QMessageBox msgBox;
            msgBox.setInformativeText(tr("Do you want to suprimpose curves over existing ones?"));
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
            msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
            msgBox.setButtonText(QMessageBox::No, tr("No"));
            msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
            msgBox.setIcon(QMessageBox::Question);
            int reply = msgBox.exec();
            if (reply == QMessageBox::Yes) {
                superimpose = true;
            }
            else if (reply == QMessageBox::No) {
                superimpose = false;
                loadedfiles.clear();
            }
            else{
                return false;
            }
        }
        else{
            superimpose = true;
        }
    }
    loadedfiles.append(filename);
    page_contours->setEnabled(true);
    page_frad->setEnabled(false);
    page_sghistogram->setEnabled(false);
    page_CPs->setEnabled(true);
    page_basins->setEnabled(true);
    Wtablacont->setVisible(true);
    FRMbonds->setEnabled(true);
    FRMcenters->setEnabled(true);
    TXTcps->clear();

//    Initializes to false the curve and sigma hole histogram plot types
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    QApplication::setOverrideCursor(Qt::WaitCursor);
    float vaux[4];
    fread(vaux , sizeof(float) , 4 , f);
    umincnt = (double) vaux[0] * ANGSTROMTOBOHR;
    umaxcnt = (double) vaux[1] * ANGSTROMTOBOHR;
    vmincnt = (double) vaux[2] * ANGSTROMTOBOHR;
    vmaxcnt = (double) vaux[3] * ANGSTROMTOBOHR;
    CHKtranspbg->setEnabled(true);
    CHKtranspbg->setHidden(false);
    if (superimpose){
        grafica->setumin(std::min(umincnt,grafica->getumin()));
        grafica->setumax(std::max(umaxcnt,grafica->getumax()));
        grafica->setvmin(std::min(vmincnt,grafica->getvmin()));
        grafica->setvmax(std::max(vmaxcnt,grafica->getvmax()));
        grafica->setwhoisabove(2);
    }
    else{
        grafica->setumin(umincnt);
        grafica->setumax(umaxcnt);
        grafica->setvmin(vmincnt);
        grafica->setvmax(vmaxcnt);
    }
    grafica->setnu(nu);
    grafica->setnv(nv);
    grafica->setaspectratio(true);
    grafica->setmulticolor(false);
    grafica->setplotcps(false);
    grafica->curvelabels->clear();
    CHKbasins->setChecked(false);
    CHKmulticolor->setChecked(false);
    CHKaspectratio->setChecked(true);
    CHKaspectratio_changed();
    if (CHKcontourslabel->isChecked()){
            grafica->setshowcontourlabels(true);
    }
    else{
            grafica->setshowcontourlabels(false);
    }
    grafica->setmulticolor(false);
    if (!(grafica->surf2d == nullpointer))
        delete(grafica->surf2d);
    grafica->surf2d = new QVector<double>(nu*nv);
    zmincnt = FLT_MAX;
    zmaxcnt = FLT_MIN;
    float *buff=new float[nu];
    grafica->surf2d->begin();
    zero = 0.;
    int kntzeros = 0;
    for (int j = 0; j < nv; ++j){
            fread(buff, sizeof(float),nu, f);
            for (int k = 0; k < nu; k++) {
                    if ((double) buff[k] == 0.)
                        kntzeros++;
                    if (zmaxcnt < (double) buff[k]) zmaxcnt = (double) buff[k];
                    if (zmincnt > (double) buff[k]) zmincnt = (double) buff[k];
                    grafica->surf2d->replace(j*nv+k,(double) buff[k]);
            }
    }
    fclose(f);
    if (std::abs(zmaxcnt-zmincnt) < 1.e-20){
            QMessageBox::warning(this,tr("Viewer")+"2D",tr("Highest contour value equal to lowest = ")
                    +QString("%1").arg(zmaxcnt));
            return false;
    }
    if (kntzeros > nv*nu/10)
        zero = 1.e-20;
    existinp = true;
    if (levels != nullpointer)
        delete(levels);
    levels = new QList<double>;
    int logz;
    int ncont = 5;
    if (zmincnt < 0){
        logz = std::min<int>(2,std::floor(std::log10(-zmincnt)));
        for(int i = 0; i < ncont ; ++i){
            if(-std::pow(10.,logz-i) > zmaxcnt)
                break;
            levels->append(-std::pow(10.,logz-i));
        }
    }
    if (zmaxcnt > 0){
        if (zmincnt < 0) levels->append(zero);
        logz = std::min<int>(2,std::floor(std::log10(zmaxcnt)));
        for(int i = 0; i < ncont ; ++i){
            if(std::pow(10.,logz-i) < zmincnt)
                break;
            levels->append(std::pow(10.,logz-i));
        }
    }
//    qSort(levels->begin(),levels->end());
    std::sort(levels->begin(),levels->end());
    SHTclist->tabla->setRowCount(levels->count()/3+1);
    for (int i = 0; i < SHTclist->tabla->rowCount() ; ++i){
        for (int j = 0 ;  j < SHTclist->tabla->columnCount() ; ++j){
            SHTclist->setcellvalue("",i,j);
        }
    }
    int knt = 0;
    for (int i = 0; knt < levels->count() && i < SHTclist->tabla->rowCount() ; ++i){
        for (int j = 0 ; knt < levels->count() && j < SHTclist->tabla->columnCount() ; ++j){
            if (levels->at(knt) != zero)
                SHTclist->setcellfloatvalue(QString::number(levels->at(knt)),i,j);
            else
                SHTclist->setcellfloatvalue(QString("0."),i,j);
            knt++;
        }
    }
    SHTclist->resizeRows(SHTclist->tabla->rowCount());
    int ini = filename.length()-10;
    if (ini > -1){
        QString planesuffix(filename.mid(ini,4));
        if (checkplanesuffix(planesuffix)){
            QFileInfo fileinfo(file.fileName());
            QString filelblsname = fileinfo.absolutePath()+"/"+fileinfo.baseName()+".plane";
            QFile filelbls(filelblsname);
            if (!filelbls.open(QFile::ReadOnly | QFile::Text)) {
                    QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                            .arg(filelblsname)+QString(":\n%1.").arg(filelbls.errorString()));
                    grafica->setexistcntatlbls(false);
            }
            else{
                grafica->uvznuc->clear();
                grafica->rcen->clear();
                QTextStream in(&filelbls);
                QString line;
                while (!in.atEnd()){
                    line = in.readLine();
#if QT_VERSION < 0x050E00
                    QStringList xy = line.split(' ',QString::SkipEmptyParts);
#else
                    QStringList xy = line.split(' ',Qt::SkipEmptyParts);
#endif
                    if (xy.count() == 6){
                        grafica->uvznuc->append(QVector3D(xy[0].toFloat(),xy[1].toFloat(),xy[2].toFloat()));
                        grafica->rcen->append(QVector3D(xy[3].toDouble(),xy[4].toDouble(),xy[5].toDouble()));
                    }
                    else if (xy.count() == 3){
                        planeA = xy[0].toDouble();
                        planeB = xy[1].toDouble();
                        planeC = xy[2].toDouble();
                    }
                }
                file.close();
                planecase = get_plane_case(planeA, planeB, planeC);
                grafica->setplaneABC(QVector3D(planeA, planeB, planeC));
                grafica->setcenterslabel();
                grafica->setbonds();
                grafica->setexistcntatlbls(true);
            }
        }
    }
    QApplication::restoreOverrideCursor();
    return true;
}

void Viewer2D::RBTcontpostyle_changed()
{
    if(RBTpossolid->isChecked())
        grafica->setpositivecontourstyle(0);
    else if(RBTposdashed->isChecked())
        grafica->setpositivecontourstyle(1);
    else if(RBTposdotted->isChecked())
        grafica->setpositivecontourstyle(2);
    else if(RBTposdashdot->isChecked())
        grafica->setpositivecontourstyle(3);

    if(RBTnegsolid->isChecked())
        grafica->setnegativecontourstyle(0);
    else if(RBTnegdashed->isChecked())
        grafica->setnegativecontourstyle(1);
    else if(RBTnegdotted->isChecked())
        grafica->setnegativecontourstyle(2);
    else if(RBTnegdashdot->isChecked())
        grafica->setnegativecontourstyle(3);

    if(RBTzerosolid->isChecked())
        grafica->setzerocontourstyle(0);
    else if(RBTzerodashed->isChecked())
        grafica->setzerocontourstyle(1);
    else if(RBTzerodotted->isChecked())
        grafica->setzerocontourstyle(2);
    else if(RBTzerodashdot->isChecked())
        grafica->setzerocontourstyle(3);

    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBtolerance_changed(int tol){
    grafica->setcontourtolerance(double(tol)/100.);
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBzerotolerance_changed(int tol){
    grafica->setzerotolerancepow(tol);
    grafica->update();
    emit moveToTop(viewernumber);
}

//***************************************************************************
//************************  FIELD LINES FUNCTIONS   *************************
//***************************************************************************


void Viewer2D::CHKarrows_changed(){
    grafica->setdrawarrows(CHKarrows->isChecked());
    if (grafica){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::CHKshowfield_changed(){
    grafica->setfieldplot(CHKshowfield->isChecked());
    if (grafica){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

//  Plots field lines
//
void Viewer2D::plot_cam2D()
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QVector<QPointF> *data=new QVector<QPointF>[elines->count()];
    grafica->efieldMap.clear();
    grafica->plotarrows->clear();
    uminelines = 1.e20;
    umaxelines = -1.e20;
    vminelines = 1.e20;
    vmaxelines = -1.e20;
    for (int i = 0 ; i < elines->count() ; ++i){
        grafica->plotarrows->append(true);
        for (int j = 0 ; j < elines->at(i).count() ; ++j){
            float x = elines->at(i).at(j).x();
            float y = elines->at(i).at(j).y();
            data[i].append(QPointF(x,y));
            if (x < uminelines) uminelines = x;
            if (x > umaxelines) umaxelines = x;
            if (y < vminelines) vminelines = y;
            if (y > vmaxelines) vmaxelines = y;
        }
    }
    grafica->setaspectratio(true);
    grafica->setuminelines(uminelines);
    grafica->setumaxelines(umaxelines);
    grafica->setvminelines(vminelines);
    grafica->setvmaxelines(vmaxelines);
    if (superimpose){
        grafica->setumin(std::min(uminelines,grafica->getumin()));
        grafica->setumax(std::max(umaxelines,grafica->getumax()));
        grafica->setvmin(std::min(vminelines,grafica->getvmin()));
        grafica->setvmax(std::max(vmaxelines,grafica->getvmax()));
        grafica->setwhoisabove(3);
    }
    else{
        grafica->setumin(uminelines);
        grafica->setumax(umaxelines);
        grafica->setvmin(vminelines);
        grafica->setvmax(vmaxelines);
    }
    QFontMetrics fm( grafica->getfontYlabel() );
    QSize fmsize;
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    grafica->setLeftMargin(max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    grafica->setRightMargin(grafica->getLeftMargin());
    fm = QFontMetrics( grafica->getfonttitle() );
    fmsize = fm.size( Qt::TextSingleLine, title );
    grafica->setTopMargin(max(50,2*fmsize.height()+20));
    grafica->setBottomMargin(max(50,2*fm.size( Qt::TextSingleLine, "0" ).height()+20));
    grafica->Ylabelrect.moveCenter(QPoint(grafica->getLeftMargin() - 2*grafica->Ylabelrect.width() - 20, grafica->height()/2));
    if (grafica->getaspectratio()){
        double w = grafica->size().rwidth();
        double h = grafica->size().rheight();
        double deltax = (w - grafica->getLeftMargin() - grafica->getRightMargin())
            / (grafica->getumax()-grafica->getumin());
        double deltay = (h - grafica->getTopMargin() - grafica->getBottomMargin())
            / (grafica->getvmax()-grafica->getvmin());

        if (deltax > deltay){
            grafica->setLeftMargin((int) (0.5*(w - deltay *
                    (grafica->getumax() - grafica->getumin()))));
            grafica->setRightMargin(grafica->getLeftMargin());
        }
        else if (deltax < deltay){
            grafica->setTopMargin((int) (0.5*(h - deltax *
                    (grafica->getvmax() - grafica->getvmin()))));
            grafica->setBottomMargin(grafica->getTopMargin());
        }
        grafica->Ylabelrect.moveCenter(QPoint(grafica->getLeftMargin() - 2*grafica->Ylabelrect.width() - 20, grafica->height()/2));
    }
    PlotSettings settings;
    settings.minX = grafica->getumin();
    settings.maxX = grafica->getumax();
    settings.minY = grafica->getvmin();
    settings.maxY = grafica->getvmax();
    settings.BottomMargin = grafica->getBottomMargin();
    settings.TopMargin = grafica->getTopMargin();
    settings.LeftMargin = grafica->getLeftMargin();
    settings.RightMargin = grafica->getRightMargin();
    settings.set_numXticks(SPBXticks->value());
    settings.set_numYticks(SPBYticks->value());
    settings.adjust();
    grafica->setxmin(settings.minX);
    grafica->setxmax(settings.maxX);
    grafica->setymin(settings.minY);
    grafica->setymax(settings.maxY);
    minx = settings.minX;
    maxx = settings.maxX;
    miny = settings.minY;
    maxy = settings.maxY;
    TXTxinf->setText(QString::number(grafica->getxmin()));
    TXTxsup->setText(QString::number(grafica->getxmax()));
    TXTyinf->setText(QString::number(grafica->getymin()));
    TXTysup->setText(QString::number(grafica->getymax()));
    grafica->setPlotSettings(settings);
    grafica->efieldcolors->clear();
    grafica->efieldpenwidth->clear();
    for (int i = 0 ; i < elines->count() ; i++ ){
        grafica->setEfieldData(i, data[i]);
        grafica->efieldcolors->append(colorForIds[8]);
        grafica->efieldpenwidth->append(SPBelinespenwidth->value());
    }
    grafica->updatecontourlabels();
    grafica->setcenterlabelrect();
//        Sets the plot type to field plot, and optionally to contourplot (superimpose), disables curve and sigma hole histrogram plot
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    if (!superimpose) grafica->setcontourplot(false);
    grafica->setfieldplot(true);
    grafica->refreshPixmap();
    grafica->updateallcenterslabelsrect();
    grafica->update();
    QApplication::restoreOverrideCursor();
    emit moveToTop(viewernumber);
}


//  Reads file with field lines
//
bool Viewer2D::readfieldfile(QString filename)
{
    QFile file(filename);
    QString suf=QFileInfo(filename).suffix().toLower();
    if ((suf != "cam2d") && (suf != "dengr2d")){
        QMessageBox::warning(this, tr("readfieldfile"),tr("Invalid extension of file %1. Must be .cam2D or dengr2D")
                .arg(filename));
        return false;
    }
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                .arg(filename)+QString(":\n%1.").arg(file.errorString()));
        return false;
    }
    if (grafica->getcontourplot()){
        if (!retrieve){
            QMessageBox msgBox;
            msgBox.setInformativeText(tr("Do you want to suprimpose curves over existing ones?"));
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
            msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
            msgBox.setButtonText(QMessageBox::No, tr("No"));
            msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
            msgBox.setIcon(QMessageBox::Question);
            int reply = msgBox.exec();
            if (reply == QMessageBox::Yes) {
                superimpose = true;
            }
            else if (reply == QMessageBox::No){
                superimpose = false;
                loadedfiles.clear();
            }
            else{
                return false;
            }
        }
        else{
            superimpose = true;
        }
    }
    loadedfiles.append(filename);
    page_field->setEnabled(true);
    page_frad->setEnabled(false);
    page_sghistogram->setEnabled(false);
    Wtablacont->setVisible(false);
    page_CPs->setEnabled(true);
    page_basins->setEnabled(true);
    FRMbonds->setEnabled(true);
    FRMcenters->setEnabled(true);
    if (elines)
        elines->clear();
    else
        elines = new QVector<QVector<QPointF> >();
    TXTcps->clear();
//    Initializes to false the curve and sigma hole histogram plot types
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QString line;
    QVector<QPointF> curve;
    QVector <QVector3D> rcen;
    grafica->centerslabel->clear();

    bool readlines = true; // While readlines is true, read lines, when false reads coordinates of centers lying on the surface
    grafica->bonds->clear();
    grafica->uvznuc->clear();
    grafica->rcen->clear();
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xy = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xy = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xy.count() == 1){
            readlines = false;
        }
        if (readlines){
            if (xy.count() < 2){
                if (curve.count() > 1)
                    elines->append(curve);
                curve.clear();
            }
            else{
                curve.append(QPointF(xy[0].toFloat(),xy[1].toFloat()));
            }
        }
        else{
            if (xy.count() == 6){
                grafica->uvznuc->append(QVector3D(xy[0].toFloat(),xy[1].toFloat(),xy[2].toFloat()));
                grafica->rcen->append(QVector3D(xy[3].toDouble(),xy[4].toDouble(),xy[5].toDouble()));
            }
            else if (xy.count() == 3){
                planeA = xy[0].toDouble();
                planeB = xy[1].toDouble();
                planeC = xy[2].toDouble();
            }
        }
        planecase = get_plane_case(planeA, planeB, planeC);
        grafica->setplaneABC(QVector3D(planeA, planeB, planeC));
        grafica->setcenterslabel();
        grafica->setbonds();
        plotcps = false;
        CHKbasins->setChecked(false);
    }
    file.close();
    QApplication::restoreOverrideCursor();
    return true;
}

void Viewer2D::SPBarrowssep_changed(int a){
    if (grafica){
        grafica->setarrowsseparation(a);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::SPBarrowssize_changed(int a){
    if (grafica){
        grafica->setarrowssize(a);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::SPBarrowsskew_changed(int a){
    if (grafica){
        grafica->setarrowsskew(a);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::SPBarrowswidth_changed(int a){
    if (grafica){
        grafica->setarrowswidth(a);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}


//***************************************************************************
//***********************  MESP SIGMA HOLE FUNCTIONS   **********************
//***************************************************************************

//    Button for choosing fonts for plotting Xlabel
void Viewer2D::BTNfontsghistlabels_click()
{
    bool OK;
    grafica->setfontsghistlabels(QFontDialog::getFont(&OK, grafica->getfontsghistlabels()));
    QFontMetrics fm( grafica->getfontsghistlabels() );
    for (int i = 0 ; i < grafica->sghistlabels->count() ; i++ ){
    QFontMetrics fm( grafica->getfontsghistlabels() );
    QSize fmsize = fm.size( Qt::TextSingleLine, grafica->sghistlabels->at(i) );
        QPoint center = grafica->sghistlabelrect->at(i).center();
    (grafica->sghistlabelrect)->replace(i,QRect(0,0, fmsize.width(), fmsize.height() ));
    (grafica->sghistlabelrect->operator [](i)).moveCenter(center);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::BTNprinthist_click(){
    printhistfile = TXTprinthist->text();
    QString fileoutstr = printhistfile;
    string v;
//      Opens file to be filled with suitable input file pointed by suffix
    QFile fileout(fileoutstr);

    if (!fileout.isOpen()){
        fileout.open(QFile::Text | QFile::WriteOnly);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    buff.append(QString("%1 ").arg(sgx.count()));
    for (int k = 0 ; k < sgx.count() ; k++){
        buff.append(QString("\n%1 ").arg(sgx.at(k)->count()));
        for (int i = 0 ; i < sgx.at(k)->count() ; i++){
            if ((i)%20 == 0){
                buff.append(QString("\n"));
            }
            buff.append(QString("%1 ").arg(sgx.at(k)->at(i)));
        }
        buff.append(QString("\n"));
        for (int i = 0 ; i < sgx.at(k)->count() ; i++){
            if ((i)%20 == 0){
                buff.append(QString("\n"));
            }
            buff.append(QString("%1 ").arg(sgdata.at(k).at(i)));
        }
        buff.append(QString("\n\n"));
    }
    outfile << buff;
    fileout.close();
    QString listhist;
    if (loadedfiles.length() > 1){
        listhist = tr("Histograms:\n");
    }
    else{
        listhist = tr("Histogram:\n");
    }
    for (int i = 0; i < loadedfiles.length(); i++){
        listhist.append(QFileInfo(loadedfiles.at(i)).fileName()+"\n");
    }
    QMessageBox::information(this,tr("DAMQT"),listhist+
        tr("written fo file\n")+QFileInfo(fileoutstr).fileName(), QMessageBox::Ok, 0);
}


//  Connects push button with file dialog for printing sigma histogram
void Viewer2D::BTNprthist_click(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getSaveFileName(this,tr("File for printing sigma hole histrogram ..."),ProjectFolder, tr("")
                + " *.hst_txt " + " (*.hst_txt);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    printhistfile = fileName;
    TXTprinthist->setText(printhistfile);
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    if (!(QFileInfo(printhistfile).suffix() == ".hst_txt")){
        printhistfile = printhistfile+".hst_txt";
    }
}

void Viewer2D::BTNsmooth_click(){
    if (CHKsmoothlines->isChecked()){
        double sfactor = smoothfactor * 0.1;
        double norm = 1. / (1. + 2.*sfactor);
        QVector<QVector<double>> sgdataaux;
        sgdataaux = sgdata;
        sgdata.clear();
        for (int i = 0 ; i < sgdataaux.count() ; i++){
            QVector<double> sgyaux;
            for (int j = 0 ; j < sgdataaux.at(i).count() ; j++){
                if (j%sgx.at(i)->count() == 0){
                    sgyaux.append((sgdataaux.at(i).at(j)+sfactor*sgdataaux.at(i).at(j+1))/(1.+sfactor));
                }
                else if (j%sgx.at(i)->count() < sgx.at(i)->count()-1){
                    sgyaux.append( (sfactor*(sgdataaux.at(i).at(j-1)+sgdataaux.at(i).at(j+1))
                            + sgdataaux.at(i).at(j))*norm);
                }
                else{
                    sgyaux.append((sgdataaux.at(i).at(j)+sfactor*sgdataaux.at(i).at(j-1))/(1.+sfactor));
                }
            }
            sgdata.append(sgyaux);
        }
        resetsghstrect = false;
        appendsghstrect = false;
        plot_sghist();
    }
}

void Viewer2D::CHKprinthist_changed(){
    if (CHKprinthist->isChecked()){
        BTNprinthist->setVisible(true);
        TXTprinthist->setVisible(true);
        BTNprthist->setVisible(true);
        LBLprthist->setVisible(true);
    }
    else{
        BTNprinthist->setVisible(false);
        TXTprinthist->setVisible(false);
        BTNprthist->setVisible(false);
        LBLprthist->setVisible(false);
    }
}

void Viewer2D::CHKsghistlabel_changed(){
    if (CHKsghistlabel->isChecked()){
        grafica->setshowsghistlabels(true);       
    }
    else{
        grafica->setshowsghistlabels(false);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKshowparthist_changed(){
    if (CHKshowparthist->isChecked()){
        int n = grafica->sghisthide->length();
        grafica->sghisthide->clear();
        for (int i = 0 ; i < n ; i++){
            grafica->sghisthide->append(false);
        }
    }
    else{
        grafica->sghisthide->clear();
        for (int i = 0 ; i < nhistpart.count() ; i++){
            grafica->sghisthide->append(false);
            for (int j = 0 ; j < nhistpart.at(i) ; j++){
                grafica->sghisthide->append(true);
            }
        }
    }
    resetsghstrect = false;
    appendsghstrect = false;
    plot_sghist();
}

void Viewer2D::CHKsmoothlines_changed(){
    sgdata.clear();
    sgdata = sgdataori;
    resetsghstrect = false;
    appendsghstrect = false;
    plot_sghist();
}

void Viewer2D::deletesghist(int indfun){
    QMessageBox msgBox;
    msgBox.setInformativeText(tr("Do you want to delete curve %1?").arg(grafica->sghistlabels->at(indfun)));
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
    msgBox.setDefaultButton(QMessageBox::Cancel);
    msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
    msgBox.setButtonText(QMessageBox::No, tr("No"));
    msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
    msgBox.setIcon(QMessageBox::Question);
    int reply = msgBox.exec();
    if (!(reply == QMessageBox::Yes)) {
        return;
    }
    if (nsghist > indfun){
        int knt = 0;

        for (int i = 0 ; i < nhistpart.count() ; i++){
            if (indfun <= knt){
                sgdata.remove(i);
                sgdataori.remove(i);
                sgx.remove(i);
                grafica->sghistlabels->remove(indfun,nhistpart.at(i)+1);
                grafica->sghisthide->remove(indfun,nhistpart.at(i)+1);
                grafica->sghistcolors->remove(indfun,nhistpart.at(i)+1);
                grafica->sghistpenstyle->remove(indfun,nhistpart.at(i)+1);
                nsghist -= nhistpart.at(i)+1;
                nhistpart.remove(i);
                break;
            }
            knt += nhistpart.at(i)+1;
        }
        resetsghstrect = true;
        appendsghstrect = false;
        plot_sghist();
    }
}

//  Plots Sigma hole histogram
//
void Viewer2D::plot_sghist()
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    grafica->curveMap.clear();
//    Sets the plot type to curve plot, disables other types
    grafica->setfradplot(false);
    grafica->setsghistplot(true);
    grafica->setcontourplot(false);
    grafica->setfieldplot(false);
    grafica->setplotcps(false);
    plotcps = false;
    CHKbasins->setChecked(false);
    sghist_data.clear();

    maxx = -1.e20;
    maxy = -1.e20;
    minx = 1.e20;
    miny = 1.e20;
    QVector<QPointF> sghistaux;
    sghistaux.clear();
    for (int k = 0 ; k < sgdata.count() ; k++){
        for (int i = 0 ; i < sgdata.at(k).count() ; ++i){
            int ix = i%sgx.at(k)->count();
            if (ix == 0){
                sghistaux.append(QPointF(2.*sgx.at(k)->at(0)-sgx.at(k)->at(1),0.));
            }
            sghistaux.append(QPointF(sgx.at(k)->at(ix),sgdata.at(k).at(i)));
            minx = std::min<double>(minx,sghistaux.last().x());
            miny = std::min<double>(miny,sghistaux.last().y());
            maxx = std::max<double>(maxx,sghistaux.last().x());
            maxy = std::max<double>(maxy,sghistaux.last().y());
            if (ix == sgx.at(k)->count()-1){
                sghistaux.append(QPointF(2.*sgx.at(k)->last()-sgx.at(k)->at(sgx.at(k)->count()-2),0.));
                sghist_data.append(sghistaux);
                sghistaux.clear();
            }
        }
    }
    grafica->setmulticolor(true);

    CHKaspectratio->setChecked(false);
    QFontMetrics fm( grafica->getfontYlabel() );
    QSize fmsize;
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    grafica->setLeftMargin(max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    grafica->setRightMargin(grafica->getLeftMargin());
    fm = QFontMetrics( grafica->getfonttitle() );
    fmsize = fm.size( Qt::TextSingleLine, title );
    grafica->setTopMargin(max(50,2*fmsize.height()+20));
    grafica->setBottomMargin(max(50,2*fm.size( Qt::TextSingleLine, "0" ).height()+20));
    grafica->Ylabelrect.moveCenter(QPoint(fmsize.height(),grafica->height()/2));
    int LeftMargin = grafica->getLeftMargin();
    int RightMargin = grafica->getRightMargin();
    int TopMargin = grafica->getTopMargin();
    int BottomMargin = grafica->getBottomMargin();
    set_SPBmargins(BottomMargin, TopMargin, LeftMargin, RightMargin);
    emit SPBLeftMargin_changed(LeftMargin);
    PlotSettings settings;
    settings.minX = minx;
    settings.maxX = maxx;
    if (minx < maxx){
    settings.minX = minx;
    settings.maxX = maxx;
    }
    else if (minx > maxx){
    settings.minX = maxx;
    settings.maxX = minx;
    }
    else{
    settings.minX = minx;
    settings.maxX = minx+1.;
    }
    if (miny < maxy){
    settings.minY = miny;
    settings.maxY = maxy;
    }
    else if (miny > maxy){
    settings.minY = maxy;
    settings.maxY = miny;
    }
    else{
    settings.minY = miny;
    settings.maxY = miny+1.;
    }
    settings.set_numXticks(SPBXticks->value());
    settings.set_numYticks(SPBYticks->value());
    settings.adjust();
    grafica->setxmin(settings.minX);
    grafica->setxmax(settings.maxX);
    grafica->setymin(settings.minY);
    grafica->setymax(settings.maxY);
    minx = settings.minX;
    maxx = settings.maxX;
    miny = settings.minY;
    maxy = settings.maxY;
    TXTxinf->setText(QString::number(grafica->getxmin()));
    TXTxsup->setText(QString::number(grafica->getxmax()));
    TXTyinf->setText(QString::number(grafica->getymin()));
    TXTysup->setText(QString::number(grafica->getymax()));
    grafica->setPlotSettings(settings);

    if (resetsghstrect){
        int ilbl = 0;
        grafica->sghistlabelrect->clear();
        grafica->sghisthaslabel->clear();
        for (int i = 0 ; i < nhistpart.count() ; i++ ){
            QFontMetrics fm( grafica->getfontsghistlabels() );
            QSize fmsize = fm.size( Qt::TextSingleLine, grafica->sghistlabels->at(ilbl) );
            grafica->sghistlabelrect->append(QRect(0,0, fmsize.width(), fmsize.height() ));
            grafica->sghistlabelrect->last().moveCenter(QPoint(4*(grafica->width())/5,
                                            grafica->getTopMargin()+(i+1)*fmsize.height()));
            grafica->sghisthaslabel->append(true);
            ilbl++;
            for (int j = 0 ; j < nhistpart.at(i) ; j++){
                grafica->sghistlabelrect->append(QRect(0,0, 0,0));  // Just a dummy
                grafica->sghisthaslabel->append(false);
                ilbl++;
            }
        }
    }
    else if (appendsghstrect){
        int ilbl = 0;
        int imin = grafica->sghistlabelrect->count();
        for (int i = 0 ; i < nhistpart.count() ; i++ ){
            if (ilbl < imin){
                ilbl += nhistpart.at(i)+1;
                continue;
            }
            QFontMetrics fm( grafica->getfontsghistlabels() );
            QSize fmsize = fm.size( Qt::TextSingleLine, grafica->sghistlabels->at(ilbl) );
            grafica->sghistlabelrect->append(QRect(0,0, fmsize.width(), fmsize.height() ));
            grafica->sghistlabelrect->last().moveCenter(QPoint(4*(grafica->width())/5,
                                            grafica->getTopMargin()+(i+1)*fmsize.height()));
            grafica->sghisthaslabel->append(true);
            ilbl++;
            for (int j = 0 ; j < nhistpart.at(i) ; j++){
                grafica->sghistlabelrect->append(QRect(0,0, 0,0));  // Just a dummy
                grafica->sghisthaslabel->append(false);
                ilbl++;
            }
        }
    }

    for (int i = 0 ; i < sghist_data.count() ; i++ ){
        grafica->setCurveData(i, sghist_data.at(i));
    }

    grafica->refreshPixmap();
    grafica->update();
    QApplication::restoreOverrideCursor();
    emit moveToTop(viewernumber);
}

//  Reads file with radial factors (and derivatives if applicable)
//
bool Viewer2D::readhstfile(QString filename)
{
    QFile file(filename);
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf != "hst"){
        QMessageBox::warning(this, tr("readhstfile"),tr("Invalid extension of file %1. Must be .hst").arg(filename));
        return false;
    }
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                    .arg(filename)+QString(":\n%1.").arg(file.errorString()));
            return false;
    }
    QTextStream in(&file);
    QString line = in.readLine();
#if QT_VERSION < 0x050E00
    QStringList header = line.split(' ',QString::SkipEmptyParts);
#else
    QStringList header = line.split(' ',Qt::SkipEmptyParts);
#endif
    if (header.count() != 1){
        QMessageBox::warning(this, tr("DAMQT"), tr("Invalid %1 file").arg(filename));
        return false;
    }
    npoints = header[0].toInt();
    page_contours->setEnabled(false);
    page_field->setEnabled(false);
    page_frad->setEnabled(false);
    page_CPs->setEnabled(false);
    page_basins->setEnabled(false);
    page_sghistogram->setEnabled(true);
    Wtablacont->setVisible(false);
    FRMbonds->setEnabled(false);
    FRMcenters->setEnabled(false);
    TXTcps->clear();
//    Initializes to false the plot types, they will be appropriately set in the corresponding plot function
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    grafica->setfieldplot(false);
    grafica->setplotcps(false);
    plotcps = false;
    CHKbasins->setChecked(false);
    if (!retrieve){
        superimpose = false;
        if (!sgdataori.isEmpty()){
            QMessageBox msgBox;
            msgBox.setInformativeText(tr("Do you want to suprimpose curves over existing ones?"));
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
            msgBox.setDefaultButton(QMessageBox::Cancel);
            msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
            msgBox.setButtonText(QMessageBox::No, tr("No"));
            msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
            msgBox.setIcon(QMessageBox::Question);
            int reply = msgBox.exec();
            if (reply == QMessageBox::Yes) {
                superimpose = true;
            }
            else if (reply == QMessageBox::No) {
                superimpose = false;
                loadedfiles.clear();
            }
            else{
                return false;
            }
        }
    }
    else{
        superimpose = true;
    }
    loadedfiles.append(filename);
    if (sgdataori.isEmpty() || !superimpose){
        sgdata.clear();
        sgdataori.clear();
        sgx.clear();
        grafica->sghistlabels->clear();
        grafica->sghistcolors->clear();
        grafica->sghistpenstyle->clear();
        grafica->sghisthide->clear();
        grafica->sghisthaslabel->clear();
        grafica->sghistlabelrect->clear();
        nsghistcolors = -1;
        Xlabel = QString("MESP / au");
        grafica->Xlabel = Xlabel;
        QFontMetrics fmx( fontXlabel );
        QSize fmxsize = fmx.size( Qt::TextSingleLine, Xlabel );
        grafica->Xlabelrect.setSize(fmxsize);
        grafica->setBottomMargin(max(50,fmxsize.height()+10));
        grafica->refreshPixmap();
        if (TXTXlabel)
            TXTXlabel->setText(Xlabel);
        Ylabel = QString("A/a")+QChar(0x2080)+QChar(0xB2);
        QFontMetrics fmori( grafica->getfontYlabel() );
        grafica->Ylabel = Ylabel;
        QFontMetrics fmy( fontYlabel );
        QSize fmysize = fmy.size( Qt::TextSingleLine, Ylabel );
        grafica->Ylabelrect.setSize(fmysize);
        if (CHKYlabelhorizontal->isChecked()){
            grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmysize.width(), fmysize.height() );
        }
        else{
            grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmysize.height(), fmysize.width() );
        }
        grafica->setLeftMargin(max(50,fmysize.height()+10));
        grafica->setRightMargin(max(50,fmysize.height()+10));
        SPBLeftMargin->setValue(grafica->getLeftMargin()+fmy.maxWidth()-fmori.maxWidth());
        if (TXTYlabel)
            TXTYlabel->setText(Ylabel);
        nhistpart.clear();
        nsghist = 0;
        resetsghstrect = true;
        appendsghstrect = false;
    }
    else{
        resetsghstrect = false;
        appendsghstrect = true;
    }
    grafica->setshowcontourlabels(false); // Hides contour labels if they exist from a previous plot
    QVector<double> *sgxaux = new QVector<double>();
    for (int i = 0 ; i < npoints ; ){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList coords = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList coords = line.split(' ',Qt::SkipEmptyParts);
#endif
        for (int j = 0 ; j < coords.count() ; ++j){
            sgxaux->append(coords[j].toDouble());
        }
        i += coords.count();
    }
    sgx.append(sgxaux);
    nsghistcolors++;
    grafica->sghistcolors->append(colorForIds[(nsghistcolors)%13]);
    grafica->sghistpenstyle->append(0);
    QString hstlabel = QFileInfo(filename).fileName();
    hstlabel.chop(4);
    grafica->sghistlabels->append(hstlabel);
    grafica->sghisthide->append(false);
    QVector<double> sgyauxori;
    for (int i = 0 ; i < npoints ; ){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList coords = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList coords = line.split(' ',Qt::SkipEmptyParts);
#endif
        for (int j = 0 ; j < coords.count() ; ++j){
                sgyauxori.append(coords[j].toDouble());
        }
        i += coords.count();
    }
    line = in.readLine();
    if (!in.atEnd()){
#if QT_VERSION < 0x050E00
        QStringList histpart = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList histpart = line.split(' ',Qt::SkipEmptyParts);
#endif
        nhistpart.append(histpart[0].toInt());
        nsghist += nhistpart.last() + 1;
        for (int k = 0 ; k < nhistpart.last() ; k++){
            grafica->sghistcolors->append(colorForIds[(nsghistcolors)%13]);
            grafica->sghistpenstyle->append(2);
            grafica->sghistlabels->append("");
            grafica->sghisthide->append(false);
            for (int i = 0 ; i < npoints ; ){
                line = in.readLine();
#if QT_VERSION < 0x050E00
                QStringList coords = line.split(' ',QString::SkipEmptyParts);
#else
                QStringList coords = line.split(' ',Qt::SkipEmptyParts);
#endif
                for (int j = 0 ; j < coords.count() ; ++j){
                        sgyauxori.append(coords[j].toDouble());
                }
                i += coords.count();
            }
        }
    }
    else{
        nhistpart.append(0);
        nsghist += nhistpart.last() + 1;
    }
    sgdataori.append(sgyauxori);
    if (CHKsmoothlines){
        double sfactor = smoothfactor * 0.1;
        double norm = 1. / (1. + 2.*sfactor);
        sgdata.clear();
        for (int i = 0 ; i < sgdataori.count() ; i++){
            QVector<double> sgyaux;
            for (int j = 0 ; j < sgdataori.at(i).count() ; j++){
                if (j%sgx.at(i)->count() == 0){
                    sgyaux.append((sgdataori.at(i).at(j)+sfactor*sgdataori.at(i).at(j+1))/(1.+sfactor));
                }
                else if (j%sgx.at(i)->count() < sgx.count()-1){
                    sgyaux.append( (sfactor*(sgdataori.at(i).at(j-1)+sgdataori.at(i).at(j+1))
                            + sgdataori.at(i).at(j))*norm);
                }
                else{
                    sgyaux.append((sgdataori.at(i).at(j)+sfactor*sgdataori.at(i).at(j-1))/(1.+sfactor));
                }
            }
            sgdata.append(sgyaux);
        }
    }
    else{
        sgdata = sgdataori;
    }
    appendsghstrect = !resetsghstrect;
    return true;
}

void Viewer2D::SPBsmoothfactor_changed(){
    smoothfactor = SPBsmoothfactor->value();
    sgdata.clear();
    sgdata = sgdataori;
    BTNsmooth_click();
}

//***************************************************************************
//************************  RADIAL FACTORS FUNCTIONS   **********************
//***************************************************************************

//    Button for choosing fonts for plotting Xlabel
void Viewer2D::BTNfontcurvelabels_click()
{
    bool OK;
    grafica->setfontcurvelabels(QFontDialog::getFont(&OK, grafica->getfontcurvelabels()));
    QFontMetrics fm( grafica->getfontcurvelabels() );
    for (int i = 0 ; i < grafica->curvelabels->count() ; i++ ){
        QFontMetrics fm( grafica->getfontcurvelabels() );
        QSize fmsize = fm.size( Qt::TextSingleLine, grafica->curvelabels->at(i) );
        QPoint center = grafica->curvelabelrect->at(i).center();
        (grafica->curvelabelrect)->replace(i,QRect(0,0, fmsize.width(), fmsize.height() ));
        (grafica->curvelabelrect->operator [](i)).moveCenter(center);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}


void Viewer2D::CHKcurvelabel_changed(){
    if (CHKcurvelabel->isChecked()){
        grafica->setshowcurvelabels(true);
    }
    else{
        grafica->setshowcurvelabels(false);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::deletefrad(int indfun){
    QMessageBox msgBox;
    msgBox.setInformativeText(tr("Do you want to delete curve %1?").arg(grafica->curvelabels->at(indfun)));
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
    msgBox.setDefaultButton(QMessageBox::Cancel);
    msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
    msgBox.setButtonText(QMessageBox::No, tr("No"));
    msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
    msgBox.setIcon(QMessageBox::Question);
    int reply = msgBox.exec();
    if (!(reply == QMessageBox::Yes)) {
        return;
    }
    if (frdata.count() > indfun){
        frdata.remove(indfun);
        frx.removeAt(indfun);
        grafica->curvelabels->remove(indfun);
        grafica->curvehide->remove(indfun);
        grafica->curvecolors->remove(indfun);
        plot_frads();
    }
}

//    Imports a file with radial factors
void Viewer2D::importfradfile_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open radial factors file ..."),ProjectFolder, tr("Radial factors files")
                + " *.frad *.drvfrad *.drv2frad " + " (*.frad *.drvfrad *.drv2frad);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    TXTfrads->setText(fileName);
    if (readfradfile(fileName)){
        plot_frads();
    }
}


//  Plots radial factors
//
void Viewer2D::plot_frads()
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    grafica->curveMap.clear();
//    Sets the plot type to curve plot, disables other types
    grafica->setfradplot(true);
    grafica->setsghistplot(false);
    grafica->setcontourplot(false);
    grafica->setfieldplot(false);
    grafica->setplotcps(false);
    plotcps = false;
    CHKbasins->setChecked(false);
    frads_data.clear();
    if (grafica->curvelabelrect != nullpointer)
        grafica->curvelabelrect->clear();
    miny = 1.e20;
    maxy = -1.e20;
    if (RBTr2flm->isChecked()){
        double aux;
        for (int k = 0 ; k < frdata.count() ; k++){
            QVector<QPointF> fradaux;
            fradaux.clear();
            for (int i = 0 ; i < frdata.at(k).count() ; ++i){
                aux = std::pow(frx.at(k)->at(i),2*ltab) * frdata.at(k).at(i);
                fradaux.append(QPointF(frx.at(k)->at(i),aux));
                if (aux < miny) miny = aux;
                if (aux > maxy) maxy = aux;
            }
            frads_data.append(fradaux);
        }
    }
    else if (RBTr2l2flm->isChecked()){
        double aux;
        for (int k = 0 ; k < frdata.count() ; k++){
            QVector<QPointF> fradaux;
            fradaux.clear();
            for (int i = 0 ; i < npoints ; ++i){
                aux = std::pow(frx.at(k)->at(i),2*ltab+2) * frdata.at(k).at(i);
                fradaux.append(QPointF(frx.at(k)->at(i),aux));
                if (aux < miny) miny = aux;
                if (aux > maxy) maxy = aux;
            }
            frads_data.append(fradaux);
        }
    }
    else{
        for (int k = 0 ; k < frdata.count() ; k++){
            QVector<QPointF> fradaux;
            fradaux.clear();
            for (int i = 0 ; i < npoints ; ++i){
                fradaux.append(QPointF(frx.at(k)->at(i),frdata.at(k).at(i)));
                if (frdata.at(k).at(i) < miny) miny = frdata.at(k).at(i);
                if (frdata.at(k).at(i) > maxy) maxy = frdata.at(k).at(i);
            }
            frads_data.append(fradaux);
        }
    }
    grafica->setmulticolor(true);
    CHKaspectratio->setChecked(false);
    CHKaspectratio_changed();
    for (int i = 0 ; i < frdata.count() ; i++ ){
        QFontMetrics fm( grafica->getfontcurvelabels() );
        QSize fmsize = fm.size( Qt::TextSingleLine, grafica->curvelabels->at(i) );
        grafica->curvelabelrect->append(QRect(0,0, fmsize.width(), fmsize.height() ));
        grafica->curvelabelrect->last().moveCenter(QPoint(4*(grafica->width())/5,grafica->getTopMargin()+(i+1)*fmsize.height()));
        grafica->setCurveData(i, frads_data.at(i));
    }
    QFontMetrics fm( grafica->getfontYlabel() );
    QSize fmsize;
    fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    grafica->setLeftMargin(max(50,fmsize.height()+fm.size( Qt::TextSingleLine, "xx" ).width()+20));
    grafica->setRightMargin(grafica->getLeftMargin());
    fm = QFontMetrics( grafica->getfonttitle() );
    fmsize = fm.size( Qt::TextSingleLine, title );
    grafica->setTopMargin(max(50,2*fmsize.height()+20));
    grafica->setBottomMargin(max(50,2*fm.size( Qt::TextSingleLine, "0" ).height()+20));
    grafica->Ylabelrect.moveCenter(QPoint(fmsize.height(),grafica->height()/2));
    int LeftMargin = grafica->getLeftMargin();
    int RightMargin = grafica->getRightMargin();
    int TopMargin = grafica->getTopMargin();
    int BottomMargin = grafica->getBottomMargin();
    set_SPBmargins(BottomMargin, TopMargin, LeftMargin, RightMargin);
    emit SPBLeftMargin_changed(LeftMargin);
    PlotSettings settings;
    settings.minX = minx;
    settings.maxX = maxx;
    if (minx < maxx){
    settings.minX = minx;
    settings.maxX = maxx;
    }
    else if (minx > maxx){
    settings.minX = maxx;
    settings.maxX = minx;
    }
    else{
    settings.minX = minx;
    settings.maxX = minx+1.;
    }
    if (miny < maxy){
    settings.minY = miny;
    settings.maxY = maxy;
    }
    else if (miny > maxy){
    settings.minY = maxy;
    settings.maxY = miny;
    }
    else{
    settings.minY = miny;
    settings.maxY = miny+1.;
    }
    settings.set_numXticks(SPBXticks->value());
    settings.set_numYticks(SPBYticks->value());
    settings.adjust();
    grafica->setxmin(settings.minX);
    grafica->setxmax(settings.maxX);
    grafica->setymin(settings.minY);
    grafica->setymax(settings.maxY);
    minx = settings.minX;
    maxx = settings.maxX;
    miny = settings.minY;
    maxy = settings.maxY;
    TXTxinf->setText(QString::number(grafica->getxmin()));
    TXTxsup->setText(QString::number(grafica->getxmax()));
    TXTyinf->setText(QString::number(grafica->getymin()));
    TXTysup->setText(QString::number(grafica->getymax()));
    grafica->setPlotSettings(settings);
    grafica->refreshPixmap();
    grafica->update();
    QApplication::restoreOverrideCursor();
    emit moveToTop(viewernumber);
}

//  Reads file with radial factors (and derivatives if applicable)
//
bool Viewer2D::readfradfile(QString filename)
{
    QFile file(filename);
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf != "frad" && suf != "drvfrad" && suf != "drv2frad"){
        QMessageBox::warning(this, tr("readfradfile"),tr("Invalid extension of file %1. Must be .frad, .drvfrad or .drv2frad")
                .arg(filename));
        return false;
    }
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                    .arg(filename)+QString(":\n%1.").arg(file.errorString()));
            return false;
    }
    QTextStream in(&file);
    QString line = in.readLine();
#if QT_VERSION < 0x050E00
    QStringList header = line.split(' ',QString::SkipEmptyParts);
#else
    QStringList header = line.split(' ',Qt::SkipEmptyParts);
#endif
    if (header.count() < 4){
        QMessageBox::warning(this, tr("DAMQT"), tr("Invalid %1 file").arg(filename));
        return false;
    }
    ltab = header[0].toInt();
    mtab = header[1].toInt();
    ncurves = header[2].toInt();
    if (ncurves < 1){
        QMessageBox::warning(this, tr("DAMQT"), tr("No tabulated curves"));
        return false;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    npoints = header[3].toInt();
    if (npoints < 2){
        QMessageBox::warning(this, tr("DAMQT"), tr("Insufficient number of tabulation points"));
    QApplication::restoreOverrideCursor();
        return false;
    }

    page_contours->setEnabled(false);
    page_field->setEnabled(false);
    page_frad->setEnabled(true);
    page_CPs->setEnabled(false);
    page_basins->setEnabled(false);
    page_sghistogram->setEnabled(false);
    Wtablacont->setVisible(false);
    FRMbonds->setEnabled(false);
    FRMcenters->setEnabled(false);
    TXTcps->clear();
//    Initializes to false the plot types, they will be appropriately set in the corresponding plot function
    grafica->setcontourplot(false);
    grafica->setsghistplot(false);
    grafica->setfieldplot(false);
    grafica->setplotcps(false);
    grafica->setshowcontourlabels(false); // Hides contour labels if they exist from a previous plot
    plotcps = false;
    CHKbasins->setChecked(false);
    CHKtranspbg->setEnabled(false);
    CHKtranspbg->setHidden(true);
    grafica->curvecontnum->clear();
    if (!retrieve){
        superimpose = false;
        if (grafica->getfradplot()){
            if (!frdata.isEmpty()){
                QMessageBox msgBox;
                msgBox.setInformativeText(tr("Do you want to suprimpose curves over existing ones?"));
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox:: Cancel);
                msgBox.setDefaultButton(QMessageBox::Cancel);
                msgBox.setButtonText(QMessageBox::Yes, tr("Yes"));
                msgBox.setButtonText(QMessageBox::No, tr("No"));
                msgBox.setButtonText(QMessageBox::Cancel, tr("Cancel"));
                msgBox.setIcon(QMessageBox::Question);
                int reply = msgBox.exec();
                if (reply == QMessageBox::Yes) {
                    superimpose = true;
                }
                else if (reply == QMessageBox::No) {
                    superimpose = false;
                    nfradcolors = 0;
                    frdata.clear();
                    loadedfiles.clear();
                }
                else{
                    return false;
                }
            }
        }
    }
    else{
        superimpose = true;
    }
    loadedfiles.append(filename);
    if (frdata.isEmpty() || !superimpose){
        frdata.clear();
        frx.clear();
        grafica->curvelabels->clear();
        grafica->curvecolors->clear();
        grafica->curvehide->clear();
        grafica->curvelabelrect->clear();
        frads_data.clear();
        frdata.clear();
        Xlabel = QString("r / au");
        grafica->Xlabel = Xlabel;
        QFontMetrics fmx( fontXlabel );
        QSize fmxsize = fmx.size( Qt::TextSingleLine, Xlabel );
        grafica->Xlabelrect.setSize(fmxsize);
        grafica->setBottomMargin(max(50,fmxsize.height()+10));
        if (TXTXlabel)
            TXTXlabel->setText(Xlabel);
        maxx = -1.e3;
        maxy = -1.e3;
        minx = 1.e3;
        miny = 1.e3;
    }
    QVector<double> *frxaux = new QVector<double>();
    for (int i = 0 ; i < npoints ; ){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList coords = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList coords = line.split(' ',Qt::SkipEmptyParts);
#endif
        for (int j = 0 ; j < coords.count() ; ++j){
            frxaux->append(coords[j].toDouble());
        }
        i += coords.count();
    }
    maxx = std::max(maxx,frxaux->last());
    minx = std::min(minx,frxaux->first());
    QFontMetrics fm( grafica->getfontcurvelabels() );
    for (int k = 0 ; k < ncurves ; k++){
        frx.append(frxaux);
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList labels = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList labels = line.split(' ',Qt::SkipEmptyParts);
#endif
        grafica->curvelabels->append(QString(labels[1])+QString(labels[0])
                +"("+QFileInfo(filename).fileName()+")");
        QSize fmsize = fm.size( Qt::TextSingleLine, grafica->curvelabels->at(k) );
        grafica->curvelabelrect->append(QRect(0,0, fmsize.width(), fmsize.height() ));
        grafica->curvelabelrect->last().moveCenter(QPoint(4*(grafica->width())/5,
                    grafica->getTopMargin()+(k+1)*fmsize.height()));
        grafica->curvehide->append(false);
        QVector<double> fryauxori;
        for (int i = 0 ; i < npoints ; ){
            line = in.readLine();
#if QT_VERSION < 0x050E00
            QStringList coords = line.split(' ',QString::SkipEmptyParts);
#else
            QStringList coords = line.split(' ',Qt::SkipEmptyParts);
#endif
            for (int j = 0 ; j < coords.count() ; ++j){
                    fryauxori.append(coords[j].toDouble());
                    maxy = std::max<double>(maxy,coords[j].toDouble());
                    miny = std::min<double>(miny,coords[j].toDouble());
            }
            i += coords.count();
        }
        frdata.append(fryauxori);
        grafica->curvecolors->append(colorForIds[(nfradcolors++)%13]);
    }
    QApplication::restoreOverrideCursor();
    return true;
}

void Viewer2D::RBTfrad_changed(){
    plot_frads();
}


//***************************************************************************
//*****************  CENTERS AND BONDS FUNCTIONS   **************************
//***************************************************************************

//    Connects push button with colors for centers labels
void Viewer2D::BTNbondscolor_clicked()
{
    QColor col = QColorDialog::getColor(bondscolor, this);
    if(col.isValid()) {
        bondscolor = col;
        grafica->setbondscolor(bondscolor);
        BTNbondscolor->setColor(&col);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

//    Connects push button with colors for centers
void Viewer2D::BTNcenterscolor_clicked()
{
    QColor col = QColorDialog::getColor(centerscolor, this);
    if(col.isValid()) {
        centerscolor = col;
        grafica->setcenterscolor(centerscolor);
        BTNcenterscolor->setColor(&col);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

//    Connects push button with colors for centers labels
void Viewer2D::BTNcenterslabelcolor_clicked()
{
    QColor col = QColorDialog::getColor(centerslabelcolor, this);
    if(col.isValid()) {
        centerslabelcolor = col;
        grafica->setcenterslabelcolor(centerslabelcolor);
        BTNcenterslabelcolor->setColor(&col);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

//    Button for choosing fonts for centers labels
void Viewer2D::BTNfontcenterslabel_click()
{
    bool OK;
    grafica->setfontcenterslabel(QFontDialog::getFont(&OK, grafica->getfontcenterslabel()));
    grafica->updateallcenterslabelsrect();
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKbonds_changed(){
    grafica->setplotbonds(CHKbonds->isChecked());
    if (grafica){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::CHKcenters_changed(){
    grafica->setplotcenters(CHKcenters->isChecked());
    if (grafica){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::CHKcenterslabel_changed(){
    grafica->setshowcenterslabel(CHKcenterslabel->isChecked());
    if (grafica){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::SPBbondswidth_changed(int a){
    if (grafica){
        grafica->setbondswidth(a);
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::SPBcenterradius_changed(){
    if (grafica){
        grafica->setcenterradius(SPBcenterradius->value());
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::SPBshftcenterlabels_changed(int a){
    if (grafica){
        grafica->setcenterslabelshft(a);
        grafica->refreshPixmap();
        grafica->updateallcenterslabelsrect();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::TXTbondthreshold_changed(){
    if (grafica){
        grafica->setbondsthr(TXTbondthreshold->text().toDouble());
        grafica->setbonds();
        grafica->update();
        emit moveToTop(viewernumber);
    }
}


//***************************************************************************
//*****************  CRITICAL POINTS FUNCTIONS   ****************************
//***************************************************************************


//    Connects push button with CPs colors
void Viewer2D::BTNcolorcps_change()
{
    QPushButton *BTNp = (QPushButton *)sender();
    int tipo=-1;
    for(int i = 0; i < max_cps; ++i){
        if(BTNcolorcps[i] == BTNp){
                tipo = i;
                break;
        }
    }

    if (tipo>-1){
        QColor colact = QColor(cpscolor[tipo].red(),cpscolor[tipo].green(),cpscolor[tipo].blue());
        QColor col = QColorDialog::getColor(colact, this);
        if(col.isValid()) {
            cpscolor[tipo].setRgb(col.red(), col.green(), col.blue(), col.alpha());
            grafica->setcpscolor(cpscolor);
            BTNcolorcps[tipo]->setColor(&cpscolor[tipo]);
            grafica->update();
            emit moveToTop(viewernumber);
        }
    }
}

//    Connects Checkbox with show/hide critical points
void Viewer2D::CHKcps_changed(){
    QCheckBox *CHKp = (QCheckBox *)sender();
    int tipo=-1;
    for(int i = 0; i < max_cps; ++i){
        if(CHKcps[i] == CHKp){
                tipo = i;
                break;
        }
    }
    grafica->setshowcps(tipo,CHKcps[tipo]->isChecked());
    grafica->update();
    emit moveToTop(viewernumber);
}

//    Imports a file with CPs Cartesian coordinates
void Viewer2D::importcpsfile_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open critical points file ..."),ProjectFolder, tr("Critical points files")
                + " *cps-?.xyz " + " (*cps-?.xyz);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    TXTcps->setText(fileName);
    if (readcpsfile(fileName)){
        plotcps = true;
        plot_cps();
    }
}

//  Plots critical points
//
void Viewer2D::plot_cps()
{
    grafica->setplotcps(true);
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QVector <QVector3D> cpsuvval;
    QVector <QColor> cpscoloraux;
    QVector<int> cpstypeaux;
    for (int i = 0 ; i < cpsxyz.count() ; i++ ){
        if (abs(QVector3D::dotProduct(cpsxyz.at(i),QVector3D(planeA, planeB, planeC))) < distthr){
            QVector2D cpsuv = xyzTouv(cpsxyz.at(i));
            cpsuvval.append(QVector3D(cpsuv,cpsval.at(i)));
            cpscoloraux.append(cpscolor[cpstype.at(i)]);
            cpstypeaux.append(cpstype.at(i));
        }
    }
    for (int i = 0 ; i < 4 ; i++){
        grafica->setshowcps(i,true);
    }
    grafica->setcpstype(cpstypeaux);
    grafica->setcpsuvval(cpsuvval);
    grafica->setcpscolor(cpscoloraux);
    grafica->update();
    QApplication::restoreOverrideCursor();
    emit moveToTop(viewernumber);
}

//  Reads file with critical points
//
bool Viewer2D::readcpsfile(QString filename)
{
    QFile file(filename);
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf != "xyz" || !filename.contains("cps")){
        QMessageBox::warning(this, tr("readcpsfile"),tr("Invalid file %1. Must be *cps-?.xyz")
                .arg(filename));
        return false;
    }
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                .arg(filename)+QString(":\n%1.").arg(file.errorString()));
        return false;
    }
    loadedfiles.append(filename);
    page_field->setEnabled(true);
    page_frad->setEnabled(false);
    page_sghistogram->setEnabled(false);
    Wtablacont->setVisible(false);
//    Initializes to false the curve and sigma hole histogram plot types
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    QTextStream in(&file);
    QString line;
    QStringList CPstypes;
    CPstypes << "x" << "y" << "z" << "m";
    cpstype.clear();
    cpsxyz.clear();
    cpsval.clear();
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xy = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xy = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xy.count() == 5){
            cpstype.append(CPstypes.indexOf(xy[0]));
            cpsxyz.append(ANGSTROMTOBOHR * QVector3D(xy[1].toFloat(),xy[2].toFloat(),xy[3].toFloat()));
            cpsval.append(xy[4].toDouble());
        }
    }
    file.close();
    return true;
}

//    Sets the ball radius
void Viewer2D::SPBcpsballradius_changed(int radius)
{
    grafica->setcpsradius(radius);
    grafica->update();
    emit moveToTop(viewernumber);
}


void Viewer2D::TXTdistancethreshold_changed(){
    distthr = TXTdistancethreshold->text().toDouble();
    if (plotcps){
        plot_cps();
    }
}

QVector2D Viewer2D::xyzTouv(QVector3D a){
    if (planecase == 1){
        return QVector2D(a.x(),a.y());
    }
    if (planecase == 2){
        return QVector2D(a.x(),a.z());
    }
    if (planecase == 3){
        return QVector2D(a.y(),a.z());
    }
    if (planecase == 4){
        return QVector2D(a.x()/wu.x(),a.z());
    }
    if (planecase == 5){
        return QVector2D(a.y(),a.x()/wu.x());
    }
    if (planecase == 6){
        return QVector2D(-a.x(),a.y()/wu.y());
    }
    else{
        return QVector2D(wu.x()*a.x()+wu.y()*a.y(),(wv.x()*a.x()+wv.y()*a.y())/(1.-wv.z()*wv.z()));
    }
}

//***************************************************************************
//**************************  BASINS FUNCTIONS   ****************************
//***************************************************************************

//    Connects push button with font colors
void Viewer2D::BTNbasinscolor_clicked()
{
    QColor colact = QColor(basinscolor->red(),basinscolor->green(),basinscolor->blue());
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        basinscolor->setRgb(col.red(), col.green(), col.blue(), 1);
        BTNbasinscolor->setColor(basinscolor);
        grafica->setbasinscolor(col);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

//    Imports a file with basins
void Viewer2D::importbasinsfile_dialog()
{
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open 2D basins file ..."),ProjectFolder, tr("Basins files")
                + " *.basin2D " + " (*.basin2D);;"+tr("All files")+" (*)");
    if (fileName.length()==0){
        QMessageBox::warning(this, tr("DAMQT"), tr("If file with 2D basins is not available:")+
                tr("carry out Molecular topography analysis,\n and run 2D Electric field or density gradient generation again"));
        return;
    }
    ProjectFolder = QFileInfo(fileName).path();
    CaptureFolder = ProjectFolder;
    TXTbasins->setText(fileName);
    if (readbasinsfile(fileName)){
        CHKbasins->setChecked(true);
        plot_basins();
    }
}

//  Plots basins
//
void Viewer2D::plot_basins()
{
    CHKbasins->setChecked(true);
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QVector<QVector <QVector2D> > basinsuv;
    for (int i = 0 ; i < basinsxyz.count() ; i++){
        QVector <QVector2D> bsuvaux;
        for (int j = 0 ; j < basinsxyz.at(i).count() ; j++){
            bsuvaux.append(basinsxyz.at(i).at(j));
        }
        basinsuv.append(bsuvaux);
    }
    grafica->setbasinsuv(basinsuv);
    grafica->refreshPixmap();
    grafica->update();
    QApplication::restoreOverrideCursor();
    emit moveToTop(viewernumber);
}


//  Reads file with basins
//
bool Viewer2D::readbasinsfile(QString filename)
{
    QFile file(filename);
    QString suf=QFileInfo(filename).suffix().toLower();
    if (suf != "basin2d"){
        QMessageBox::warning(this, tr("readbasinsfile"),tr("Invalid extension of file %1. Must be .basin2D")
                .arg(filename));
        return false;
    }
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                .arg(filename)+QString(":\n%1.").arg(file.errorString()));
        return false;
    }
    loadedfiles.append(filename);
    page_field->setEnabled(true);
    page_frad->setEnabled(false);
    page_sghistogram->setEnabled(false);
    Wtablacont->setVisible(false);
//    Initializes to false the curve and sigma hole histogram plot types
    grafica->setfradplot(false);
    grafica->setsghistplot(false);
    QTextStream in(&file);
    QString line;
    basinsxyz.clear();
    QVector<QVector2D> bsxyz;
    bsxyz.clear();
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xy = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xy = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xy.count() == 2){
            bsxyz.append(QVector2D(xy[0].toFloat(),xy[1].toFloat()));
        }
        else if (xy.count() < 2 && bsxyz.count() > 0){
            basinsxyz.append(bsxyz);
            bsxyz.clear();
        }
    }
    file.close();
    return true;
}


//***************************************************************************
//***************************  TITLE FUNCTIONS   ****************************
//***************************************************************************

void Viewer2D::CHKtitle_changed(){
    if (CHKtitle->isChecked()){
        grafica->setshowtitle(true);
        grafica->update();
        emit moveToTop(viewernumber);
    }
    else{
        grafica->setshowtitle(false);
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::CHKtitlebox_changed(){
    if (CHKtitlebox->isChecked()){
        grafica->setshowtitlebox(true);
        grafica->update();
        emit moveToTop(viewernumber);
    }
    else{
        grafica->setshowtitlebox(false);
        grafica->update();
        emit moveToTop(viewernumber);
    }
}

void Viewer2D::TXTtitlepos_changed(){
    TXTxtitle->setText(QString("%1").arg(grafica->titlerect.center().x()));
    TXTytitle->setText(QString("%1").arg(grafica->titlerect.center().y()));
}

void Viewer2D::TXTxtitle_changed(){
    grafica->titlerect.moveCenter(QPoint(TXTxtitle->text().toInt(),grafica->titlerect.center().y()));
    emit CHKtitle_changed();
}

void Viewer2D::TXTytitle_changed(){
    grafica->titlerect.moveCenter(QPoint(grafica->titlerect.center().x(),TXTytitle->text().toInt()));
    emit CHKtitle_changed();
}

void Viewer2D::aceptar_title()
{
    title = TXTtitle->text();
    grafica->title = title;
    grafica->setfonttitle(fonttitle);
    QFontMetrics fm( grafica->getfonttitle() );
    grafica->titlerect.setSize(fm.size( Qt::TextSingleLine, title));
    grafica->setTopMargin(std::max<int>(grafica->getTopMargin(),fm.height()+20));
    grafica->setfonttitlecolor(fonttitlecolor);
    BTNcolortitlefont->setColor(&fonttitlecolor);
    CHKtitle->setChecked(!CHKtitlehide->isChecked());
    emit CHKtitle_changed();
    cerrar_title();
}

void Viewer2D::cerrar_title()
{
    if (FRMtitle != nullpointer){
        delete FRMtitle;
        FRMtitle = nullpointer;
    }
    fonttitlecolor = grafica->getfonttitlecolor();
}

void Viewer2D::editar_title()
{
    FRMtitle = new QDialog();
    FRMtitle->setWindowTitle(tr("Title"));
    FRMtitle->setFixedHeight(60);
    FRMtitle->setAttribute(Qt::WA_DeleteOnClose);

    TXTtitle = new QLineEdit();
    TXTtitle->setText(title);
        
    BTNfonttitle = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connect(BTNfonttitle, SIGNAL(clicked()), this, SLOT(BTNfonttitle_click()));
    
    grafica->setfonttitlecolor(fonttitlecolor);
    BTNcolortitlefont = new ColorButton();
    BTNcolortitlefont->setIcon(QIcon(":/images/fonts48.png"));
    BTNcolortitlefont->setText(tr("Color"));
    BTNcolortitlefont->setColor(&fonttitlecolor);
    BTNcolortitlefont->setEnabled(true);
    connect(BTNcolortitlefont, SIGNAL(clicked()), this, SLOT(BTNcolortitlefont_clicked()));
        
        
    QPushButton *BTNtitleOK = new QPushButton();
    BTNtitleOK->setText(tr("Accept"));
    connect(BTNtitleOK, SIGNAL(clicked()), this, SLOT(aceptar_title()));
    QPushButton *BTNtitleCancel = new QPushButton();
    BTNtitleCancel->setText(tr("Cancel"));
    connect(BTNtitleCancel, SIGNAL(clicked()), this, SLOT(cerrar_title()));
    CHKtitlehide = new QCheckBox(tr("Hide"));
    CHKtitlehide->setChecked(!CHKtitle->isChecked());
    CHKtitlehide->setEnabled(true); 
    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(TXTtitle);
    Layout1->addWidget(BTNfonttitle);
    Layout1->addWidget(BTNcolortitlefont);
    Layout1->addWidget(CHKtitlehide);
    Layout1->addWidget(BTNtitleOK);
    Layout1->addWidget(BTNtitleCancel);
    QVBoxLayout *Layout2 = new QVBoxLayout(FRMtitle);
    Layout2->addLayout(Layout1);
    FRMtitle->exec();
}

//    Button for choosing fonts for plot title
void Viewer2D::BTNfonttitle_click()
{ 
    bool OK;
    fonttitle = QFontDialog::getFont(&OK, grafica->getfonttitle());
}

//    Connects push button with font colors
void Viewer2D::BTNcolortitlefont_clicked()
{
    QColor col = QColorDialog::getColor(fonttitlecolor, this);
    if(col.isValid()) {
        fonttitlecolor = col;
                BTNcolortitlefont->setColor(&fonttitlecolor);
    }
}

//***************************************************************************
//********************  X AXIS LABEL FUNCTIONS    ***************************
//***************************************************************************

void Viewer2D::CHKXlabel_changed(){
    if (CHKXlabel->isChecked()){
        grafica->setshowXlabel(true);   
    }
    else{
        grafica->setshowXlabel(false);
    }
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKXlabelbox_changed(){
    if (CHKXlabelbox->isChecked()){
        grafica->setshowXlabelbox(true);
    }
    else{
        grafica->setshowXlabelbox(false);
    }
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::TXTXlabelpos_changed(){
    TXTxXlabel->setText(QString("%1").arg(grafica->Xlabelrect.center().x()));
    TXTyXlabel->setText(QString("%1").arg(grafica->Xlabelrect.center().y()));
}

void Viewer2D::TXTxXlabel_changed(){
    grafica->Xlabelrect.moveCenter(QPoint(TXTxXlabel->text().toInt(),grafica->Xlabelrect.center().y()));
    emit CHKXlabel_changed();
}

void Viewer2D::TXTyXlabel_changed(){
    grafica->Xlabelrect.moveCenter(QPoint(grafica->Xlabelrect.center().x(),TXTyXlabel->text().toInt()));
    emit CHKXlabel_changed();
}

void Viewer2D::aceptar_Xlabel()
{
    Xlabel = TXTXlabel->text();
    grafica->Xlabel = Xlabel;
    grafica->setfontXlabel(fontXlabel);
    QFontMetrics fm( grafica->getfontXlabel() );
    grafica->Xlabelrect.setSize(fm.size( Qt::TextSingleLine, Xlabel));
    QSize fmsize = fm.size( Qt::TextSingleLine, Xlabel );
    grafica->setBottomMargin(max(50,fmsize.height()+10));
        
    grafica->setfontXlabelcolor(fontXlabelcolor);
    BTNcolorXlabelfont->setColor(&fontXlabelcolor);
    CHKXlabel->setChecked(!CHKXhide->isChecked());
    grafica->refreshPixmap();

    emit CHKXlabel_changed();
    cerrar_Xlabel();
    grafica->updateplot();
    emit moveToTop(viewernumber);
}

void Viewer2D::cerrar_Xlabel()
{
    if (FRMXlabel != nullpointer){
        delete FRMXlabel;
        FRMXlabel = nullpointer;
    }
    fontXlabelcolor = grafica->getfontXlabelcolor();
}

void Viewer2D::editar_Xlabel()
{
    FRMXlabel = new QDialog();
    FRMXlabel->setWindowTitle(tr("Title"));
    FRMXlabel->setFixedHeight(60);
    FRMXlabel->setAttribute(Qt::WA_DeleteOnClose);

    TXTXlabel = new QLineEdit();
    TXTXlabel->setText(Xlabel);
        
    BTNfontXlabel = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connect(BTNfontXlabel, SIGNAL(clicked()), this, SLOT(BTNfontXlabel_click()));
    grafica->setfontXlabelcolor(fontXlabelcolor);
    BTNcolorXlabelfont = new ColorButton();
    BTNcolorXlabelfont->setIcon(QIcon(":/images/fonts48.png"));
    BTNcolorXlabelfont->setText(tr("Color"));
    BTNcolorXlabelfont->setColor(&fontXlabelcolor);
    BTNcolorXlabelfont->setEnabled(true);
    connect(BTNcolorXlabelfont, SIGNAL(clicked()), this, SLOT(BTNcolorXlabelfont_clicked()));
        
    QPushButton *BTNXlabelOK = new QPushButton();
    BTNXlabelOK->setText(tr("Accept"));
    connect(BTNXlabelOK, SIGNAL(clicked()), this, SLOT(aceptar_Xlabel()));
    QPushButton *BTNXlabelCancel = new QPushButton();
    BTNXlabelCancel->setText(tr("Cancel"));
    connect(BTNXlabelCancel, SIGNAL(clicked()), this, SLOT(cerrar_Xlabel()));
        CHKXhide = new QCheckBox(tr("Hide"));
    CHKXhide->setChecked(!CHKXlabel->isChecked());
    CHKXhide->setEnabled(true); 
        
    QHBoxLayout *Layout1 = new QHBoxLayout();
    Layout1->addWidget(TXTXlabel);
        Layout1->addWidget(BTNfontXlabel);
        Layout1->addWidget(BTNcolorXlabelfont);
        Layout1->addWidget(CHKXhide);
    Layout1->addWidget(BTNXlabelOK);
    Layout1->addWidget(BTNXlabelCancel);
    QVBoxLayout *Layout2 = new QVBoxLayout(FRMXlabel);
    Layout2->addLayout(Layout1);
    FRMXlabel->exec();
}

//    Button for choosing fonts for plotting Xlabel
void Viewer2D::BTNfontXlabel_click()
{ 
    bool OK;
    fontXlabel = QFontDialog::getFont(&OK, grafica->getfontXlabel());
}

//    Connects push button with font colors
void Viewer2D::BTNcolorXlabelfont_clicked()
{
    QColor col = QColorDialog::getColor(fontXlabelcolor, this);
    if(col.isValid()) {
        fontXlabelcolor = col;
        BTNcolorXlabelfont->setColor(&fontXlabelcolor);
    }
}

//***************************************************************************
//********************  Y AXIS LABEL FUNCTIONS    ***************************
//***************************************************************************

void Viewer2D::CHKYlabel_changed(){
    if (CHKYlabel->isChecked()){
        grafica->setshowYlabel(true);
    }
    else{
        grafica->setshowYlabel(false);
    }
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKYlabelbox_changed(){
    if (CHKYlabelbox->isChecked()){
        grafica->setshowYlabelbox(true);
    }
    else{
        grafica->setshowYlabelbox(false);
    }
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::TXTYlabelpos_changed(){
    TXTxYlabel->setText(QString("%1").arg(grafica->Ylabelrect.center().x()));
    TXTyYlabel->setText(QString("%1").arg(grafica->Ylabelrect.center().y()));
}

void Viewer2D::TXTxYlabel_changed(){
    grafica->Ylabelrect.moveCenter(QPoint(TXTxYlabel->text().toInt(),grafica->Ylabelrect.center().y()));
}

void Viewer2D::TXTyYlabel_changed(){
    grafica->Ylabelrect.moveCenter(QPoint(grafica->Ylabelrect.center().x(),TXTyYlabel->text().toInt()));
}

void Viewer2D::CHKYlabelhorizontal_changed(){
    QFontMetrics fm( grafica->getfontYlabel() );
    QSize fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    if (CHKYlabelhorizontal->isChecked()){
        grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmsize.width(), fmsize.height() );
        grafica->setYlabelvert(false);
    }
    else{
        grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmsize.height(), fmsize.width() );
        grafica->setYlabelvert(true);
    }
}

void Viewer2D::aceptar_Ylabel()
{
    Ylabel = TXTYlabel->text();
    grafica->Ylabel = Ylabel;
    QFontMetrics fmori( grafica->getfontYlabel() );
    grafica->setfontYlabel(fontYlabel);
    QFontMetrics fm( grafica->getfontYlabel() );
    grafica->Ylabelrect.setSize(fm.size( Qt::TextSingleLine, Ylabel));
    QSize fmsize = fm.size( Qt::TextSingleLine, Ylabel );
    grafica->setLeftMargin(max(grafica->getLeftMargin(),max(50,fmsize.height()+10)));
    grafica->setRightMargin(max(grafica->getRightMargin(),max(50,fmsize.height()+10)));
    if (CHKYlabelhorizontal->isChecked()){
        grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmsize.width(), fmsize.height() );
    }
    else{
        grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmsize.height(), fmsize.width() );
    }
    SPBLeftMargin->setValue(grafica->getLeftMargin()+fm.maxWidth()-fmori.maxWidth());
    CHKYlabel->setChecked(!CHKYhide->isChecked());
    CHKYlabelhorizontal->setChecked(CHKYhorizontal->isChecked());
    grafica->setfontYlabelcolor(fontYlabelcolor);
    BTNcolorYlabelfont->setColor(&fontYlabelcolor);
    grafica->refreshPixmap();
    
//    emit CHKYlabel_changed();
    cerrar_Ylabel();
    grafica->updateplot();
    emit moveToTop(viewernumber);
}

void Viewer2D::cerrar_Ylabel()
{
    delete FRMYlabel;
    FRMYlabel = nullpointer;
    fontYlabelcolor = grafica->getfontYlabelcolor();
}

void Viewer2D::editar_Ylabel()
{
//     lang = new IDIOMA();
    FRMYlabel = new QDialog();
    FRMYlabel->setWindowTitle(tr("Title"));
    FRMYlabel->setFixedHeight(60);
    FRMYlabel->setAttribute(Qt::WA_DeleteOnClose);

    TXTYlabel = new QLineEdit();
    TXTYlabel->setText(Ylabel);
        
        BTNfontYlabel = new QPushButton(QIcon(":/images/fonts48.png"),tr("Font"));
    connect(BTNfontYlabel, SIGNAL(clicked()), this, SLOT(BTNfontYlabel_click()));
    
    grafica->setfontYlabelcolor(fontYlabelcolor);
    BTNcolorYlabelfont = new ColorButton();
    BTNcolorYlabelfont->setIcon(QIcon(":/images/fonts48.png"));
    BTNcolorYlabelfont->setText(tr("Color"));
    BTNcolorYlabelfont->setColor(&fontYlabelcolor);
    BTNcolorYlabelfont->setEnabled(true);
    connect(BTNcolorYlabelfont, SIGNAL(clicked()), this, SLOT(BTNcolorYlabelfont_clicked()));
    CHKYhide = new QCheckBox(tr("Hide"));
    CHKYhide->setChecked(!CHKYlabel->isChecked());
    CHKYhide->setEnabled(true); 
    CHKYhorizontal = new QCheckBox(tr("Horizontal"));
    CHKYhorizontal->setChecked(CHKYlabelhorizontal->isChecked());
    CHKYhorizontal->setEnabled(true);      
    QPushButton *BTNYlabelOK = new QPushButton();
    BTNYlabelOK->setText(tr("Accept"));
    connect(BTNYlabelOK, SIGNAL(clicked()), this, SLOT(aceptar_Ylabel()));
    QPushButton *BTNYlabelCancel = new QPushButton();
    BTNYlabelCancel->setText(tr("Cancel"));
    connect(BTNYlabelCancel, SIGNAL(clicked()), this, SLOT(cerrar_Ylabel())); 
    QVBoxLayout *Layout1 = new QVBoxLayout();
    Layout1->addWidget(CHKYhide);
    Layout1->addWidget(CHKYhorizontal);
    QHBoxLayout *Layout2 = new QHBoxLayout();
    Layout2->addWidget(TXTYlabel);
    Layout2->addWidget(BTNfontYlabel);
    Layout2->addWidget(BTNcolorYlabelfont);
    Layout2->addLayout(Layout1);
    Layout2->addWidget(BTNYlabelOK);
    Layout2->addWidget(BTNYlabelCancel);
    QVBoxLayout *Layout3 = new QVBoxLayout(FRMYlabel);
    Layout3->addLayout(Layout2);
    FRMYlabel->exec();
}

//    Button for choosing fonts for plotting Ylabel
void Viewer2D::BTNfontYlabel_click()
{ 
    bool OK;
    fontYlabel = QFontDialog::getFont(&OK, grafica->getfontYlabel());
}

//    Connects push button with font colors
void Viewer2D::BTNcolorYlabelfont_clicked()
{
    QColor col = QColorDialog::getColor(fontYlabelcolor, this);
    if(col.isValid()) {
        fontYlabelcolor = col;
        BTNcolorYlabelfont->setColor(&fontYlabelcolor);
    }
}

//***************************************************************************
//********************* ZOOM REGION FUNCTIONS    ****************************
//***************************************************************************
void Viewer2D::execzoomreg(){
    grafica->setxmin(TXTxinf->text().toDouble());
    grafica->setxmax(TXTxsup->text().toDouble());
    grafica->setymin(TXTyinf->text().toDouble());
    grafica->setymax(TXTysup->text().toDouble());
    grafica->setzoomRegion();
}

void Viewer2D::zoom_changed(){
    TXTxinf->setText(QString::number(grafica->getxmin()));
    TXTxsup->setText(QString::number(grafica->getxmax()));
    TXTyinf->setText(QString::number(grafica->getymin()));
    TXTysup->setText(QString::number(grafica->getymax()));
}

//***************************************************************************
//****************************  PEN FUNCTIONS    ****************************
//***************************************************************************

void Viewer2D::SPBbasinspenwidth_changed(){
    grafica->setbasinspenwidth(SPBbasinspenwidth->value());
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBcntpenwidth_changed(){
    grafica->setcntpenwidth(SPBcntpenwidth->value());
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBcurvespenwidth_changed(){
    grafica->setcurvespenwidth(SPBcurvespenwidth->value());
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBelinespenwidth_changed(){
    grafica->efieldpenwidth->clear();
    if (!elines) return;
    for (int i = 0 ; i < elines->count() ; i++ ){
        grafica->efieldpenwidth->append(SPBelinespenwidth->value());
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBsghistpenwidth_changed(){
    grafica->setsghistpenwidth(SPBsghistpenwidth->value());
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

//***************************************************************************
//***************  BACKGROUND COLORS FUNCTIONS    ***************************
//***************************************************************************

//    Connects push button with background colors
void Viewer2D::BTNbkgcolor_clicked()
{
    QColor colact = QColor(bkgcolor->red(),bkgcolor->green(),bkgcolor->blue());
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        bkgcolor->setRgb(col.red(), col.green(), col.blue(), 1);
        BTNbkgcolor->setColor(bkgcolor);
        grafica->setbackgroundcolor(col);
        grafica->Pal->setColor(QPalette::Background, bkgcolor->rgb());
        grafica->setPalette(*(grafica->Pal));
        grafica->setAutoFillBackground(true);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKtranspbg_changed(){
    if (CHKtranspbg->isChecked())
        grafica->settranspbackground(true);
    else
        grafica->settranspbackground(false);
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

//***************************************************************************
//***************  GRID COLORS FUNCTIONS    ***************************
//***************************************************************************

//    Connects push button with font colors
void Viewer2D::BTNgridcolor_clicked()
{
    QColor colact = QColor(gridcolor->red(),gridcolor->green(),gridcolor->blue());
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        gridcolor->setRgb(col.red(), col.green(), col.blue(), 1);
        BTNgridcolor->setColor(gridcolor);
        grafica->setgridcolor(col);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKgrid_changed(){
    if (grafica->getshowgrid()){
    grafica->setshowgrid(false);
    }
    else{
        grafica->setshowgrid(true);
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::RBTsolid_changed(){
    if (RBTsolid->isChecked())
        grafica->setsolidgrid(true);
    else
        grafica->setsolidgrid(false);
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}
//***************************************************************************
//***********************  SCALES FUNCTIONS    ******************************
//***************************************************************************

//    Connect check button with scales or ticks changes 
void Viewer2D::CHKXscale_changed(){
    if (CHKXscalebottom->isChecked()){
        grafica->setshowXscalebottom(true);
        if (CHKticks->isChecked()) 
            grafica->setshowXticksbottom(true);
        else
            grafica->setshowXticksbottom(false);
    }
    else{
        grafica->setshowXscalebottom(false);
        grafica->setshowXticksbottom(false);    
    }
    if (CHKXscaletop->isChecked()){
        grafica->setshowXscaletop(true);
        if (CHKticks->isChecked()) 
            grafica->setshowXtickstop(true);
        else
            grafica->setshowXtickstop(false);
    }
    else{
        grafica->setshowXscaletop(false);
        grafica->setshowXtickstop(false);    
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::CHKYscale_changed(){
    if (CHKYscaleleft->isChecked()){
        grafica->setshowYscaleleft(true);
        if (CHKticks->isChecked()) 
            grafica->setshowYticksleft(true);
        else
            grafica->setshowYticksleft(false);
    }
    else{
        grafica->setshowYscaleleft(false);
        grafica->setshowYticksleft(false);    
    }
    if (CHKYscaleright->isChecked()){
        grafica->setshowYscaleright(true);
        if (CHKticks->isChecked()) 
            grafica->setshowYticksright(true);
        else
            grafica->setshowYticksright(false);
    }
    else{
        grafica->setshowYscaleright(false);
        grafica->setshowYticksright(false);    
    }
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}
//***************************************************************************
//***********************  MARGINS FUNCTIONS    ******************************
//***************************************************************************

void Viewer2D::SPBBottomMargin_changed(int bottom)
{
    grafica->setBottomMargin(bottom);
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBLeftMargin_changed(int left)
{
    grafica->setLeftMargin(left);
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBRightMargin_changed(int right)
{
    grafica->setRightMargin(right);
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBTopMargin_changed(int top)
{
    grafica->setTopMargin(top);
    grafica->refreshPixmap();
    grafica->update();
    emit moveToTop(viewernumber);
}

void Viewer2D::set_SPBmargins(int BottomMargin, int TopMargin, int LeftMargin, int RightMargin){
    SPBBottomMargin->setMaximum(grafica->height()/2);
    SPBTopMargin->setMaximum(grafica->height()/2);
    SPBLeftMargin->setMaximum(grafica->width()/2);
    SPBRightMargin->setMaximum(grafica->width()/2);
    SPBBottomMargin->setValue(BottomMargin);
    SPBLeftMargin->setValue(LeftMargin);
    SPBRightMargin->setValue(RightMargin);
    SPBTopMargin->setValue(TopMargin);
    
    
}
//***************************************************************************
//***********************  TICKS FUNCTIONS    ******************************
//***************************************************************************
void Viewer2D::CHKautomaticticks_changed()
{
    if (CHKautomaticticks->isChecked()){
        SPBXticks->setEnabled(false);
        SPBYticks->setEnabled(false);
        LBLXticks->setEnabled(false);
        LBLYticks->setEnabled(false);
        grafica->setXcifras(SPBXcifras->value());
        grafica->setYcifras(SPBYcifras->value());
        grafica->setautomaticticks(true);
        grafica->restore_automatic();
        grafica->refreshPixmap();
        grafica->update();
        TXTxinf->setText(QString("%1").arg(grafica->getxmin()));
        TXTxsup->setText(QString("%1").arg(grafica->getxmax()));
        TXTyinf->setText(QString("%1").arg(grafica->getymin()));
        TXTysup->setText(QString("%1").arg(grafica->getymax()));
    }
    else{
        SPBXticks->setEnabled(true);
        SPBYticks->setEnabled(true);
        LBLXticks->setEnabled(true);
        LBLYticks->setEnabled(true);
        grafica->setXcifras(SPBXcifras->value());
        grafica->setYcifras(SPBYcifras->value());
        grafica->setautomaticticks(false);
        grafica->Xticks_changed(SPBXticks->value());
        grafica->Yticks_changed(SPBYticks->value());
    } 
    this->repaint();
    emit moveToTop(viewernumber);
}

void Viewer2D::SPBXcifras_changed(int xcifras){
    grafica->setXcifras(xcifras);
    if (CHKautomaticticks->isChecked()){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
    else{
        grafica->Xticks_changed(SPBXticks->value());
    }
}

void Viewer2D::SPBYcifras_changed(int ycifras){
    grafica->setYcifras(ycifras);
    if (CHKautomaticticks->isChecked()){
        grafica->refreshPixmap();
        grafica->update();
        emit moveToTop(viewernumber);
    }
    else{
        grafica->Yticks_changed(SPBYticks->value());
    }
}
//***************************************************************************
//*********************  CAPTURE FUNCTIONS   ***************************
//***************************************************************************
//    Function for image capture
void Viewer2D::BTNcapture_click()
{
    if (!grafica)
        return;
    if (RBTscaledef->isChecked()){
        grafica->setuserdefresol(true);
        QSize sizenew = TXTscalesize->text().toDouble() * grafica->size();
        if (sizenew.width() > HIGHEST_RESOL || sizenew.height() > HIGHEST_RESOL){
            QMessageBox::warning(this, tr("DAMQT"),
                    tr("Size chosen exceeds the highest allowable (%1 x %2)\nTakes the highest compatible ").arg(HIGHEST_RESOL).arg(HIGHEST_RESOL));
            double scale = ((double) HIGHEST_RESOL) / ((double)(max(sizenew.width(),sizenew.height())));
            sizenew = (scale *  TXTscalesize->text().toDouble()) * grafica->size();

        }
        grafica->setimagesize(sizenew);
    }
    else if (RBTuserdef->isChecked()){
        grafica->setuserdefresol(true);
        QSize sizenew = QSize(TXThsize->text().toInt(),TXTvsize->text().toInt());
        if (sizenew.width() > HIGHEST_RESOL || sizenew.height() > HIGHEST_RESOL){
            QMessageBox::warning(this, tr("DAMQT"),
                    tr("Size chosen exceeds the highest allowable (%1 x %2)\nTakes the highest compatible ").arg(HIGHEST_RESOL).arg(HIGHEST_RESOL));
            double scale = ((double) HIGHEST_RESOL) / ((double)(max(sizenew.width(),sizenew.height())));
            sizenew = (scale *  TXTscalesize->text().toDouble()) * grafica->size();

        }
        grafica->setimagesize(sizenew);
    }
    else{
        grafica->setuserdefresol(false);
    }
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(CaptureFolder);
    QString fileName = filedialog.getSaveFileName(this,tr("Save image as ..."), CaptureFolder,
            (tr("Image files") + ": *.png, *.jpg, *.bmp, *.jpeg, *.ppm, *.xbm, *.xpm, *.tiff"
            + " (*.png *.jpg *.bmp *.jpeg *.ppm *.xbm *.xpm *.tiff);;"+tr("All files")+" (*)"));
    grafica->plotfile = fileName;
    grafica->setimagequality(SPBimagequality->value());
    if (!grafica->plotfile.isEmpty()){
        QString filter=QFileInfo(grafica->plotfile).suffix();
        if (filter.toUpper() == "PNG") grafica->tipo = "PNG";
        else if (filter.toUpper() == "JPG") grafica->tipo = "JPG";
        else if (filter.toUpper() == "BMP") grafica->tipo = "BMP";    
        else if (filter.toUpper() == "JPEG") grafica->tipo = "JPEG";    
        else if (filter.toUpper() == "PPM") grafica->tipo = "PPM";
        else if (filter.toUpper() == "XBM") grafica->tipo = "XBM";
        else if (filter.toUpper() == "XPM") grafica->tipo = "XPM";    
        else if (filter.toUpper() == "TIFF") grafica->tipo = "TIFF";
        else grafica->tipo = nullpointer; // Default: tries to imagine the format
        grafica->capture_image();
    }    
}

void Viewer2D::CHKtranspbgcapture_changed(){
    if (CHKtranspbgcapture->isChecked())
        grafica->settranspbckgrcapture(true);
    else
        grafica->settranspbckgrcapture(false);
}

void Viewer2D::RBTscreendef_changed()
{
    if (RBTscreendef->isChecked()){
        TXTscalesize->setVisible(false);
        LBLscaledef->setVisible(false);
        TXThsize->setVisible(false);
        TXTvsize->setVisible(false);
        LBLpor->setVisible(false);
    }
    else if (RBTscaledef->isChecked()){
        TXTscalesize->setVisible(true);
        LBLscaledef->setVisible(true);
        TXThsize->setVisible(false);
        TXTvsize->setVisible(false);
        LBLpor->setVisible(false);
    }
    else{
        TXTscalesize->setVisible(false);
        LBLscaledef->setVisible(false);
        TXThsize->setVisible(true);
        TXTvsize->setVisible(true);
        LBLpor->setVisible(true);
    }

}

//***************************************************************************
//***************  RETRIEVE/SAVE SETTINGS FUNCTIONS   ***********************
//***************************************************************************
void Viewer2D::BTNsaveSettings_click(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getSaveFileName(this,tr("Open file ..."),ProjectFolder, tr("Geometry and basis set files")
            + " *.2Dsettings" + " (*.2Dsettings);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    if (!fileName.contains(".2Dsettings"))
        fileName.append(".2Dsettings");
    QFile fileout(fileName);
    if (!fileout.open(QFile::WriteOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("DAMQT"),tr("File %1 cannot be read")
                                    .arg(fileName)+QString(":\n%1.").arg(fileout.errorString()));
            return;
    };
    if (!fileout.isOpen()){
            fileout.open(QFile::Text | QFile::WriteOnly);
    }
    QTextStream outfile(&fileout); // Buffer for writing to fileout
    QByteArray buff;
    if (CHKsettingsfile->isChecked()){
        buff.append("[V2DFILE]\n");
        QString files = loadedfiles.at(0);
        for (int i = 1 ; i < loadedfiles.length() ; i++){
            files.append(QString(",%1").arg(loadedfiles.at(i)));
        }
        buff.append("files="+files+"\n");
        if (grafica){
            if (!grafica->curvelabels->isEmpty()){
                QString labels = grafica->curvelabels->at(0);
                for (int i = 1 ; i < grafica->curvelabels->length() ; i++){
                    labels.append(QString(",%1").arg(grafica->curvelabels->at(i)));
                }
                buff.append("labels="+labels+"\n");
            }
            if (!grafica->sghistlabels->isEmpty()){
                QString labels = grafica->sghistlabels->at(0);
                for (int i = 1 ; i < grafica->sghistlabels->length() ; i++){
                    labels.append(QString(",%1").arg(grafica->sghistlabels->at(i)));
                }
                buff.append("labels="+labels+"\n");
            }
        }
    }
    if (CHKsettingstitle->isChecked()){
        buff.append("[V2DTITLE]\n");
        buff.append("chktitle="+QString("%1").arg(CHKtitle->isChecked())+"\n");
        buff.append("xtitle="+QString("%1").arg(TXTxtitle->text().toDouble())+"\n");
        buff.append("ytitle="+QString("%1").arg(TXTytitle->text().toDouble())+"\n");
        buff.append("chktitlebox="+QString("%1").arg(CHKtitlebox->isChecked())+"\n");
        buff.append("titlecolorR="+QString("%1").arg(fonttitlecolor.red())+"\n");
        buff.append("titlecolorG="+QString("%1").arg(fonttitlecolor.green())+"\n");
        buff.append("titlecolorB="+QString("%1").arg(fonttitlecolor.blue())+"\n");
        buff.append("titlefontfamily="+QString(fonttitle.family())+"\n");
        buff.append("titlefontsize="+QString("%1").arg(fonttitle.pointSize())+"\n");
        buff.append("titlefontweight="+QString("%1").arg(fonttitle.weight())+"\n");
        buff.append("titlefontitalic="+QString("%1").arg(fonttitle.italic())+"\n");
    }
    if (CHKsettingsXlabel->isChecked()){
        buff.append("[V2DXAXIS]\n");
        buff.append("chkXlabel="+QString("%1").arg(CHKXlabel->isChecked())+"\n");
        buff.append("xXlabel="+QString("%1").arg(TXTxXlabel->text().toDouble())+"\n");
        buff.append("yXlabel="+QString("%1").arg(TXTyXlabel->text().toDouble())+"\n");
        buff.append("chkXlabelbox="+QString("%1").arg(CHKXlabelbox->isChecked())+"\n");
        buff.append("XlabelcolorR="+QString("%1").arg(fontXlabelcolor.red())+"\n");
        buff.append("XlabelcolorG="+QString("%1").arg(fontXlabelcolor.green())+"\n");
        buff.append("XlabelcolorB="+QString("%1").arg(fontXlabelcolor.blue())+"\n");
        buff.append("Xlabelfontfamily="+QString(fontXlabel.family())+"\n");
        buff.append("Xlabelfontsize="+QString("%1").arg(fontXlabel.pointSize())+"\n");
        buff.append("Xlabelfontweight="+QString("%1").arg(fontXlabel.weight())+"\n");
        buff.append("Xlabelfontitalic="+QString("%1").arg(fontXlabel.italic())+"\n");
    }
    if (CHKsettingsYlabel->isChecked()){
        buff.append("[V2DYAXIS]\n");
        buff.append("chkYlabel="+QString("%1").arg(CHKYlabel->isChecked())+"\n");
        buff.append("xYlabel="+QString("%1").arg(TXTxYlabel->text().toDouble())+"\n");
        buff.append("yYlabel="+QString("%1").arg(TXTyYlabel->text().toDouble())+"\n");
        buff.append("chkYlabelbox="+QString("%1").arg(CHKYlabelbox->isChecked())+"\n");
        buff.append("chkYlabelhorizontal="+QString("%1").arg(CHKYlabelhorizontal->isChecked())+"\n");
        buff.append("YlabelcolorR="+QString("%1").arg(fontYlabelcolor.red())+"\n");
        buff.append("YlabelcolorG="+QString("%1").arg(fontYlabelcolor.green())+"\n");
        buff.append("YlabelcolorB="+QString("%1").arg(fontYlabelcolor.blue())+"\n");
        buff.append("Ylabelfontfamily="+QString(fontYlabel.family())+"\n");
        buff.append("Ylabelfontsize="+QString("%1").arg(fontYlabel.pointSize())+"\n");
        buff.append("Ylabelfontweight="+QString("%1").arg(fontYlabel.weight())+"\n");
        buff.append("Ylabelfontitalic="+QString("%1").arg(fontYlabel.italic())+"\n");
    }
    if (CHKsettingsZoom->isChecked()){
        buff.append("[V2DZOOM]\n");
        buff.append("xinfzoom="+QString("%1").arg(TXTxinf->text().toDouble())+"\n");
        buff.append("xsupzoom="+QString("%1").arg(TXTxsup->text().toDouble())+"\n");
        buff.append("yinfzoom="+QString("%1").arg(TXTyinf->text().toDouble())+"\n");
        buff.append("ysupzoom="+QString("%1").arg(TXTysup->text().toDouble())+"\n");
    }
    if (CHKsettingsMargins->isChecked()){
        buff.append("[V2DMARGINS]\n");
        buff.append("bottommargin="+QString("%1").arg(SPBBottomMargin->value())+"\n");
        buff.append("topmargin="+QString("%1").arg(SPBTopMargin->value())+"\n");
        buff.append("leftmargin="+QString("%1").arg(SPBLeftMargin->value())+"\n");
        buff.append("rightmargin="+QString("%1").arg(SPBRightMargin->value())+"\n");
    }
    if (CHKsettingscntPen->isChecked()){
        buff.append("[V2DPEN]\n");
        buff.append("cntpenwidth="+QString("%1").arg(SPBcntpenwidth->value())+"\n");
        buff.append("elinespenwidth="+QString("%1").arg(SPBelinespenwidth->value())+"\n");
    }
    if (CHKsettingsbkg->isChecked()){
        buff.append("[V2DBKG]\n");
        buff.append("bkgcolorR="+QString("%1").arg(bkgcolor->red())+"\n");
        buff.append("bkgcolorG="+QString("%1").arg(bkgcolor->green())+"\n");
        buff.append("bkgcolorB="+QString("%1").arg(bkgcolor->blue())+"\n");
    }
    outfile << buff << endl;    
    fileout.close();
}

void Viewer2D::BTNretrieveSettings_click(){
    QFileDialog filedialog(this);
    filedialog.setWindowFlags(Qt::WindowStaysOnTopHint);
    filedialog.setDirectory(ProjectFolder);
    QString fileName = filedialog.getOpenFileName(this,tr("Open file ..."),ProjectFolder, tr("Geometry and basis set files")
            + " *.2Dsettings" + " (*.2Dsettings);;"+tr("All files")+" (*)");
    if (fileName.length()==0) return;
    QFile filein(fileName);
    if (!filein.open(QFile::ReadOnly | QFile::Text)) {
            return;
    };
    QString suf = QFileInfo(fileName).suffix().toLower();
    if (!(suf == "2dsettings")) {
        QMessageBox::warning(this, tr("DAMQT"),tr("Extension of file %1 not allowed")
                .arg(fileName));
            return;
    };
    ProjectFolder = QFileInfo(fileName).path();
    retrieve = true;
    string file = QString(fileName).toStdString();
    QString qv;
    string v;
    int rojo, verde, azul;
    if (CHKsettingscntPen->isChecked()){
        v = CIniFile::GetValue("cntpenwidth","V2DPEN",file);
        if (!v.empty()){
            SPBcntpenwidth->setValue(QString(v.c_str()).toInt());
        }
    }

    if (CHKsettingstitle->isChecked()){
        v = CIniFile::GetValue("chktitle","V2DTITLE",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKtitle->setChecked(true);
            }
            else{
               CHKtitle->setChecked(false);
            }
        }
        v = CIniFile::GetValue("xtitle","V2DTITLE",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTxtitle->setText(qv);
        }
        v = CIniFile::GetValue("ytitle","V2DTITLE",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTytitle->setText(qv);
        }
        v = CIniFile::GetValue("chktitlebox","V2DTITLE",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKtitlebox->setChecked(true);
            }
            else{
               CHKtitlebox->setChecked(false);
            }
        }
        v = CIniFile::GetValue("titlecolorR","V2DTITLE",file);
        if (!v.empty())
            rojo = QString(v.c_str()).toInt();
        else
            rojo = fonttitlecolor.red();
        v = CIniFile::GetValue("titlecolorG","V2DTITLE",file);
        if (!v.empty())
            verde = QString(v.c_str()).toInt();
        else
            verde = fonttitlecolor.green();
        v = CIniFile::GetValue("titlecolorB","V2DTITLE",file);
        if (!v.empty())
            azul = QString(v.c_str()).toInt();
        else
            azul = fonttitlecolor.blue();
        fonttitlecolor.setRgb(rojo,verde,azul);
        grafica->setfonttitlecolor(fonttitlecolor);

        v = CIniFile::GetValue("titlefontfamily","V2DTITLE",file);
        if (!v.empty()){
            fonttitle.setFamily(QString(v.c_str()));
            v = CIniFile::GetValue("titlefontsize","V2DTITLE",file);
            if (!v.empty())
                fonttitle.setPointSize(QString(v.c_str()).toInt());
            v = CIniFile::GetValue("titlefontweight","V2DTITLE",file);
            if (!v.empty())
                fonttitle.setWeight(QString(v.c_str()).toInt());
            v = CIniFile::GetValue("titlefontitalic","V2DTITLE",file);
            if (!v.empty())
                fonttitle.setItalic(QString(v.c_str()).toInt());
            grafica->setfonttitle(fonttitle);
            QFontMetrics fm( grafica->getfonttitle() );
            grafica->titlerect.setSize(fm.size( Qt::TextSingleLine, title));
            grafica->setTopMargin(std::max<int>(grafica->getTopMargin(),fm.height()+20));
        }
    }
    if (CHKsettingsXlabel->isChecked()){
        v = CIniFile::GetValue("chkXlabel","V2DXAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKXlabel->setChecked(true);
            }
            else{
               CHKXlabel->setChecked(false);
            }
        }
        v = CIniFile::GetValue("xXlabel","V2DXAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTxXlabel->setText(qv);
        }
        v = CIniFile::GetValue("yXlabel","V2DXAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTyXlabel->setText(qv);
        }
        v = CIniFile::GetValue("chkXlabelbox","V2DXAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKXlabelbox->setChecked(true);
            }
            else{
               CHKXlabelbox->setChecked(false);
            }
        }
        v = CIniFile::GetValue("XlabelcolorR","V2DXAXIS",file);
        if (!v.empty())
            rojo = QString(v.c_str()).toInt();
        else
            rojo = fontXlabelcolor.red();
        v = CIniFile::GetValue("XlabelcolorG","V2DXAXIS",file);
        if (!v.empty())
            verde = QString(v.c_str()).toInt();
        else
            verde = fontXlabelcolor.green();
        v = CIniFile::GetValue("XlabelcolorB","V2DXAXIS",file);
        if (!v.empty())
            azul = QString(v.c_str()).toInt();
        else
            azul = fontXlabelcolor.blue();
        fontXlabelcolor.setRgb(rojo,verde,azul);
        grafica->setfontXlabelcolor(fontXlabelcolor);

        v = CIniFile::GetValue("Xlabelfontfamily","V2DXAXIS",file);
        if (!v.empty()){
            fontXlabel.setFamily(QString(v.c_str()));
            v = CIniFile::GetValue("Xlabelfontsize","V2DXAXIS",file);
            if (!v.empty())
                fontXlabel.setPointSize(QString(v.c_str()).toInt());
            v = CIniFile::GetValue("Xlabelfontweight","V2DXAXIS",file);
            if (!v.empty())
                fontXlabel.setWeight(QString(v.c_str()).toInt());
            v = CIniFile::GetValue("Xlabelfontitalic","V2DXAXIS",file);
            if (!v.empty())
                fontXlabel.setItalic(QString(v.c_str()).toInt());
            grafica->setfontXlabel(fontXlabel);
            QFontMetrics fm( grafica->getfontXlabel() );
            grafica->Xlabelrect.setSize(fm.size( Qt::TextSingleLine, Xlabel));
            QSize fmsize = fm.size( Qt::TextSingleLine, Ylabel );
            grafica->setBottomMargin(max(50,fmsize.height()+10));
        }
    }
    if (CHKsettingsYlabel->isChecked()){
        v = CIniFile::GetValue("chkYlabel","V2DYAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKYlabel->setChecked(true);
            }
            else{
               CHKYlabel->setChecked(false);
            }
        }
        v = CIniFile::GetValue("xYlabel","V2DYAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTxYlabel->setText(qv);
        }
        v = CIniFile::GetValue("yYlabel","V2DYAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTyYlabel->setText(qv);
        }
        v = CIniFile::GetValue("chkYlabelbox","V2DYAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKYlabelbox->setChecked(true);
            }
            else{
               CHKYlabelbox->setChecked(false);
            }
        }
        v = CIniFile::GetValue("chkYlabelhorizontal","V2DYAXIS",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            if (qv.toInt()==1){
               CHKYlabelhorizontal->setChecked(true);
            }
            else{
               CHKYlabelhorizontal->setChecked(false);
            }
        }
        v = CIniFile::GetValue("YlabelcolorR","V2DYAXIS",file);
        if (!v.empty())
            rojo = QString(v.c_str()).toInt();
        else
            rojo = fontYlabelcolor.red();
        v = CIniFile::GetValue("YlabelcolorG","V2DYAXIS",file);
        if (!v.empty())
            verde = QString(v.c_str()).toInt();
        else
            verde = fontYlabelcolor.green();
        v = CIniFile::GetValue("YlabelcolorB","V2DYAXIS",file);
        if (!v.empty())
            azul = QString(v.c_str()).toInt();
        else
            azul = fontYlabelcolor.blue();
        fontYlabelcolor.setRgb(rojo,verde,azul);
        grafica->setfontYlabelcolor(fontYlabelcolor);

        v = CIniFile::GetValue("Ylabelfontfamily","V2DYAXIS",file);
        if (!v.empty()){
            fontYlabel.setFamily(QString(v.c_str()));
            v = CIniFile::GetValue("Ylabelfontsize","V2DYAXIS",file);
            if (!v.empty())
                fontYlabel.setPointSize(QString(v.c_str()).toInt());
            v = CIniFile::GetValue("Ylabelfontweight","V2DYAXIS",file);
            if (!v.empty())
                fontYlabel.setWeight(QString(v.c_str()).toInt());
            v = CIniFile::GetValue("Ylabelfontitalic","V2DYAXIS",file);
            if (!v.empty())
                fontYlabel.setItalic(QString(v.c_str()).toInt());
            QFontMetrics fmori( grafica->getfontYlabel() );
            grafica->setfontYlabel(fontYlabel);
            QFontMetrics fm( grafica->getfontYlabel() );
            grafica->Ylabelrect.setSize(fm.size( Qt::TextSingleLine, Ylabel));
            QSize fmsize = fm.size( Qt::TextSingleLine, Ylabel );
//            grafica->setLeftMargin(max(50,fmsize.height()+10));
//            grafica->setRightMargin(max(50,fmsize.height()+10));
            if (CHKYlabelhorizontal->isChecked()){
                grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmsize.width(), fmsize.height() );
            }
            else{
                grafica->Ylabelrect = QRect(grafica->Ylabelrect.x(), grafica->Ylabelrect.y(), fmsize.height(), fmsize.width() );
            }
            SPBLeftMargin->setValue(grafica->getLeftMargin()+fm.maxWidth()-fmori.maxWidth());
        }
    }
    if (CHKsettingsZoom->isChecked()){
        v = CIniFile::GetValue("xinfzoom","V2DZOOM",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTxinf->setText(qv);
            grafica->setxmin(TXTxinf->text().toInt());
        }
        v = CIniFile::GetValue("xsupzoom","V2DZOOM",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTxsup->setText(qv);
            grafica->setxmax(TXTxsup->text().toInt());
        }
        v = CIniFile::GetValue("yinfzoom","V2DZOOM",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTyinf->setText(qv);
            grafica->setymin(TXTyinf->text().toInt());
        }
        v = CIniFile::GetValue("ysupzoom","V2DZOOM",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            TXTysup->setText(qv);
            grafica->setymax(TXTysup->text().toInt());
        }
        execzoomreg();
    }
    if (CHKsettingsMargins->isChecked()){
        v = CIniFile::GetValue("bottommargin","V2DMARGINS",file);
        if (!v.empty()){
            SPBBottomMargin->setValue(QString(v.c_str()).toInt());
        }
        v = CIniFile::GetValue("topmargin","V2DMARGINS",file);
        if (!v.empty()){
            SPBTopMargin->setValue(QString(v.c_str()).toInt());
        }
        v = CIniFile::GetValue("leftmargin","V2DMARGINS",file);
        if (!v.empty()){
            SPBLeftMargin->setValue(QString(v.c_str()).toInt());
        }
        v = CIniFile::GetValue("rightmargin","V2DMARGINS",file);
        if (!v.empty()){
            SPBRightMargin->setValue(QString(v.c_str()).toInt());
        }
        grafica->setLeftMargin(SPBLeftMargin->value());
        grafica->setRightMargin(SPBRightMargin->value());
        grafica->setTopMargin(SPBTopMargin->value());
        grafica->setBottomMargin(SPBBottomMargin->value());
    }
    if (CHKsettingselinesPen->isChecked()){
        v = CIniFile::GetValue("elinespenwidth","V2DPEN",file);
        if (!v.empty()){
            SPBelinespenwidth->setValue(QString(v.c_str()).toInt());
        }
    }
    if (CHKsettingsbkg->isChecked()){
        v = CIniFile::GetValue("bkgcolorR","V2DBKG",file);
        if (!v.empty())
            rojo = QString(v.c_str()).toInt();
        else
            rojo = bkgcolor->red();
        v = CIniFile::GetValue("bkgcolorG","V2DBKG",file);
        if (!v.empty())
            verde = QString(v.c_str()).toInt();
        else
            verde = bkgcolor->green();
        v = CIniFile::GetValue("bkgcolorB","V2DBKG",file);
        if (!v.empty())
            azul = QString(v.c_str()).toInt();
        else
            azul = bkgcolor->blue();
        bkgcolor->setRgb(rojo,verde,azul);
        BTNbkgcolor->setColor(bkgcolor);
        grafica->Pal->setColor(QPalette::Background, bkgcolor->rgb());
        grafica->setPalette(*(grafica->Pal));
        grafica->setAutoFillBackground(true);
        grafica->refreshPixmap();
    }
    if (CHKsettingsfile->isChecked()){
        loadedfiles.clear();
        v = CIniFile::GetValue("files","V2DFILE",file);
        QString suf;
        if (!v.empty()){
            qv=QString(v.c_str());
            QStringList filelist = qv.split(",");
            frdata.clear();
            frx.clear();
            sgdata.clear();
            sgdataori.clear();
            sgx.clear();
            for (int i = 0 ; i < filelist.length() ; i++){
                QString file = filelist.at(i);
                suf = QFileInfo(file).suffix().toLower();
                if (suf == "basin2d"){
                    readbasinsfile(file);
                    plot_basins();
                    ProjectFolder = QFileInfo(file).path();
                    TXTbasins->setText(file);
                }
                else if (suf == "cnt"){
                    readcntfile(file);
                    plot_cnt();
                    ProjectFolder = QFileInfo(file).path();
                    TXTcontour->setText(file);
                }
                else if (suf == "cam2d" || suf == "dengr2d"){
                    readfieldfile(file);
                    plot_cam2D();
                    ProjectFolder = QFileInfo(file).path();
                    TXTfield->setText(file);
                }
                else if (suf == "frad" || suf == "drvfrad" || suf == "drv2frad"){
                    readfradfile(file);
                    plot_frads();
                    ProjectFolder = QFileInfo(file).path();
                    TXTfrads->setText(file);
                }
                else if (suf == "hst"){
                    readhstfile(file);
                    plot_sghist();
                    ProjectFolder = QFileInfo(file).path();
                    TXTsghist->setText(file);
                }
                else if (suf == "xyz" && file.contains("cps-")){
                    readcpsfile(file);
                    plot_cps();
                    ProjectFolder = QFileInfo(file).path();
                    TXTcps->setText(file);
                }
                else{
                    QMessageBox::warning(this, tr("DAMQT"),tr("Extension of file %1 not allowed")
                            .arg(fileName));
                }
            }
        }
        v = CIniFile::GetValue("labels","V2DFILE",file);
        if (!v.empty()){
            qv=QString(v.c_str());
            QStringList labelslist = qv.split(",");
            if (suf == "hst"){
                grafica->sghistlabels->clear();
                for (int i = 0 ; i < labelslist.length() ; i++){
                    grafica->sghistlabels->append(labelslist.at(i));
                }
            }
            else if (suf == "frad" || suf == "drvfrad" || suf == "drv2frad"){
                grafica->curvelabels->clear();
                for (int i = 0 ; i < labelslist.length() ; i++){
                    grafica->curvelabels->append(labelslist.at(i));
                }
            }
            grafica->refreshPixmap();
        }
    }
    CaptureFolder = ProjectFolder;
    filein.close();
    retrieve = false;
}

//***************************************************************************
//***************************  ANCILLARY FUNCTIONS  *************************
//***************************************************************************
int Viewer2D::get_plane_case(double a, double b, double c){
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
                wv = QVector3D(0.0, -c, b) / QVector3D(0.0,-c,b).length();
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
                wv = QVector3D(-c, 0.0, a) / QVector3D(-c,0.0,a).length();
                return 5;
            }
        }
        else{
            if (std::abs(c) < 1.e-7){   // A != 0, B != 0, C == 0
                wu = QVector3D(-b, a, 0.0) / QVector3D(-b,a,0.0).length();
                wv = QVector3D(0.0, 0.0, 1.0);
                return 4;
            }
            else{   // A != 0, b != 0, C != 0
                wu = QVector3D(-b, a, 0.0) / QVector3D(-b,a,0.0).length();
                wv = QVector3D(-a*c, -b*c, a*a+b*b) / (QVector3D(a,b,0.0).length() * QVector3D(a,b,c).length());
                return 7;
            }
        }
    }
}

bool Viewer2D::isvisible(){
    return visible;
}

void Viewer2D::emit_moveToTopPlotter(){
    emit moveToTop(viewernumber);
}

QPoint Viewer2D::get_position(){
    if (mainWin)
        return QPoint(mainWin->x(),mainWin->y());
    else
        return QPoint(0,0);
}

QSize Viewer2D::get_size(){
    if (mainWin)
        return mainWin->size();
    else
        return QSize(0,0);
}

int Viewer2D::getviewernumber(){
    return viewernumber;
}

QString Viewer2D::get_viewername(){
    return viewername;
}

void Viewer2D::lower_mainwindow(){
    if (mainWin){
        mainWin->lower();
    }
}

void Viewer2D::raise_mainwindow(){
    if (mainWin){
        mainWin->raise();
    }
}

void Viewer2D::set_position(QPoint a){
    if (mainWin){
        position = a;
        mainWin->move(position);
    }
}

void Viewer2D::set_size(QSize a){
    if (mainWin){
        size = a;
        mainWin->resize(size);
    }
}

void Viewer2D::set_visible(bool a){
    visible = a;
    if (mainWin){
        if (visible)
            mainWin->show();
        else
            mainWin->hide();
    }
}

void Viewer2D::set_ProjectFolder(QString a){
    ProjectFolder = a;
    CaptureFolder = ProjectFolder;
}

void Viewer2D::set_ProjectName(QString a){
    ProjectName = a;
}

void Viewer2D::set_viewername(QString a){
    viewername = a;
}

void Viewer2D::showMainWindow(){
    visible = true;
    if (mainWin)
        mainWin->show();
}

//
//    MainWindow2DViewer SLOTS
//

MainWindow2DViewer::MainWindow2DViewer(QWidget *parent) : QMainWindow()
{
}

MainWindow2DViewer::~MainWindow2DViewer(){
}

void MainWindow2DViewer::reject(){
    emit hideplotter();
}

void MainWindow2DViewer::closeEvent(QCloseEvent *event){
    emit hideplotter();
    event->accept();
}
