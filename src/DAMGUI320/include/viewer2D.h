//  Copyright 2008-2016, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//  File:   Viewer2D.h
//
//      Last version: March 2016
//

#ifndef VIEWER2D_H
#define	VIEWER2D_H

#include <QApplication>
#include <QCloseEvent>
#include <QColorDialog>
#include <QDir>
#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDialog>
#include <QtGlobal>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLayout>
#include <QList>
#include <QMainWindow>
#include <QMessageBox>
#include <QObject>
#include <QRadioButton>
#include <QScrollArea>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QVector2D>
#include <QVector3D>
#include <QWidget>

#include <cstring>

#include "ColorButton.h"
#include "contours.h"
#include "float.h"
#include "plotter.h"
#include "Sheet.h"

class QByteArray;
class QCheckBox;
class QFile;
class QFont;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QRadioButton;
class QScrollArea;
class QToolBox;
class QToolButton;
class ContourLabel;
class ContourRect;

#define HIGHEST_RESOL 8192

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif


class MainWindow2DViewer : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow2DViewer(QWidget *parent = 0);
    ~MainWindow2DViewer();
signals:
    void hideplotter();
protected:
    virtual void reject();
    virtual void closeEvent(QCloseEvent *event);

};

class Viewer2D : public QWidget {
	Q_OBJECT
public:
    explicit Viewer2D(QString *name, QWidget *parent = 0);
//    Viewer2D(QWidget *parent = 0);
    virtual ~Viewer2D();
    QString CaptureFolder;
    QString ProjectFolder;
    QString ProjectName;
	bool existinp;
	int nu;
	int nv;
	QList<double> *levels;
    double umaxcnt;
    double umaxelines;
    double umincnt;
    double uminelines;
    double vmaxcnt;
    double vmaxelines;
    double vmincnt;
    double vminelines;
    double zmaxcnt;
    double zmincnt;
    int ltab;
    int mtab;
signals:
    void hideplotter();
    void moveToTop(int);
public slots:
    bool isvisible();
    QPoint get_position();
    QSize get_size();
    QString get_viewername();
    int getviewernumber();
    void lower_mainwindow();
    void raise_mainwindow();
    void set_position(QPoint);
    void set_ProjectFolder(QString);
    void set_ProjectName(QString);
    void set_size(QSize);
    void set_viewername(QString);
    void set_SPBmargins(int BottomMargin, int TopMargin, int LeftMargin, int RightMargin);
    void set_visible(bool);
    void showMainWindow();
private slots:
    void aceptar_title();
    void aceptar_Xlabel();
    void aceptar_Ylabel();
    void CHKautomaticticks_changed();
    void cerrar_title();
    void cerrar_Xlabel();
    void cerrar_Ylabel();
    void deletefrad(int);
    void deletesghist(int);
    void editar_title();
    void editar_Xlabel();
    void editar_Ylabel();
    void execzoomreg();
// 		Buttons functions
    void BTNactualizarcontour_click();
    void BTNbkgcolor_clicked();
    void BTNbasinscolor_clicked();
    void BTNbondscolor_clicked();
    void BTNcapture_click();
    void BTNcenterscolor_clicked();
    void BTNcenterslabelcolor_clicked();
    void BTNcolorcps_change();
    void BTNcolortitlefont_clicked();
    void BTNcolorXlabelfont_clicked();
    void BTNcolorYlabelfont_clicked();
    void BTNcontourslabelcolor_clicked();
    void BTNfontcenterslabel_click();
    void BTNfontcontourslabel_click();
    void BTNfontcurvelabels_click();
    void BTNfontsghistlabels_click();
    void BTNfonttitle_click();
    void BTNfontXlabel_click();
    void BTNfontYlabel_click();
    void BTNgridcolor_clicked();
    void BTNnegativecontourcolor_clicked();
    void BTNpositivecontourcolor_clicked();
    void BTNprinthist_click();
    void BTNprthist_click();
    void BTNretrieveSettings_click();
    void BTNsaveSettings_click();
    void BTNsmooth_click();
    void BTNzerocontourcolor_clicked();
// 		Check boxes functions
    void CHKarrows_changed();
    void CHKaspectratio_changed();
    void CHKbasins_changed();
    void CHKbonds_changed();
    void CHKcenters_changed();
    void CHKcenterslabel_changed();
    void CHKcontourslabel_changed();
    void CHKcps_changed();
    void CHKcurvelabel_changed();
    void CHKgrid_changed();
    void CHKmulticolor_changed();
    void CHKprinthist_changed();
    void CHKsghistlabel_changed();
    void CHKshowfield_changed();
    void CHKshowparthist_changed();
    void CHKsmoothlines_changed();
    void CHKtitle_changed();
    void CHKtitlebox_changed();
    void CHKtranspbg_changed();
    void CHKtranspbgcapture_changed();
    void CHKXscale_changed();
    void CHKYscale_changed();
    void CHKXlabel_changed();
    void CHKXlabelbox_changed();
    void CHKYlabel_changed();
    void CHKYlabelbox_changed();
    void CHKYlabelhorizontal_changed();
    void emit_moveToTopPlotter();
//      Layouts pages
    void page_basins_layouts();
    void page_capture_layouts();
    void page_contours_layouts();
    void page_CPs_layouts();
    void page_field_layouts();
    void page_frad_layouts();
    void page_options_layouts();
    void page_save_layouts();
    void page_sghistogram_layouts();
//      Widgets pages
    void page_basins_widgets();
    void page_capture_widgets();
    void page_contours_widgets();
    void page_CPs_widgets();
    void page_field_widgets();
    void page_frad_widgets();
    void page_options_widgets();
    void page_save_widgets();
    void page_sghistogram_widgets();
//		QLine Edit functions
    void TXTbondthreshold_changed();
    void TXTdistancethreshold_changed();
    void TXTtitlepos_changed();
    void TXTxtitle_changed();
    void TXTytitle_changed();
    void TXTXlabelpos_changed();
    void TXTxXlabel_changed();
    void TXTyXlabel_changed();
    void TXTYlabelpos_changed();
    void TXTxYlabel_changed();
    void TXTyYlabel_changed();
// 		Radio button functions
    void RBTcontpostyle_changed();
    void RBTfrad_changed();
    void RBTscreendef_changed();
    void RBTsolid_changed();
// 		Spin boxes functions
    void SPBarrowssep_changed(int);
    void SPBarrowssize_changed(int);
    void SPBarrowsskew_changed(int);
    void SPBarrowswidth_changed(int);
    void SPBbondswidth_changed(int);
    void SPBcenterradius_changed();
    void SPBbasinspenwidth_changed();
    void SPBcntpenwidth_changed();
    void SPBcpsballradius_changed(int radius);
    void SPBcurvespenwidth_changed();
    void SPBelinespenwidth_changed();
    void SPBBottomMargin_changed(int);
    void SPBLeftMargin_changed(int);
    void SPBRightMargin_changed(int);
    void SPBsghistpenwidth_changed();
    void SPBshftcenterlabels_changed(int);
    void SPBsmoothfactor_changed();
    void SPBtolerance_changed(int);
    void SPBTopMargin_changed(int);
    void SPBXcifras_changed(int);
    void SPBYcifras_changed(int);
    void SPBzerotolerance_changed(int);
// 		Import files functions
    void importbasinsfile_dialog();
    void importcntfile_dialog();
    void importcpsfile_dialog();
    void importhstfile_dialog();
    void importfieldfile_dialog();
    void importfradfile_dialog();
//              Zoom functions
    void zoom_changed();
	
private:
    static const int max_cps=4;

    bool appendsghstrect;
    bool plotcps;
    bool resetsghstrect;
    bool retrieve;
    bool superimpose;
    bool visible;

    int ncurves;
    int nfradcolors;
    int nsghist;
    int nsghistcolors;
    int npoints;
    int planecase;
    int viewernumber;
    int smoothfactor;

    double distthr;
    double minx;
    double maxx;
    double maxy;
    double miny;
    double planeA;
    double planeB;
    double planeC;
    double zero;

    MainWindow2DViewer *mainWin;

    Plotter *grafica;
    QVector<int> cpstype;
    QVector<int> nhistpart;
    QVector<QColor> cpscolor;
    QVector<double> *data;
    QVector<QVector<QPointF> > frads_data;
    QVector<QVector<QPointF> > sghist_data;
    QVector<QString> loadedfiles;
    QVector<QVector<double> > frdata;
    QVector<QVector<double> > sgdata;
    QVector<QVector<double> > sgdataori;
    QVector<QVector<double>* > frx;
    QVector<QVector<double>* > sgx;
    QVector<QVector<QPointF> > *elines;
    QVector <QVector<QVector2D> > basinsxyz;
    QVector <QVector3D> cpsxyz;
    QVector <double> cpsval;
    QVector3D wu;
    QVector3D wv;
//		ColorButton variables
    ColorButton *BTNbkgcolor;
    ColorButton *BTNbasinscolor;
    ColorButton *BTNbondscolor;
    ColorButton *BTNcenterscolor;
    ColorButton *BTNcenterslabelcolor;
    ColorButton *BTNcolorcps[max_cps];
    ColorButton *BTNcolortitlefont;
    ColorButton *BTNcolorXlabelfont;
    ColorButton *BTNcolorYlabelfont;
    ColorButton *BTNcontourslabelcolor;
    ColorButton *BTNgridcolor;
    ColorButton *BTNnegativecontourcolor;
    ColorButton *BTNpositivecontourcolor;
    ColorButton *BTNzerocontourcolor;
	
// 		QCheckbox variables
    QCheckBox *CHKarrows;
    QCheckBox *CHKautomaticticks;
    QCheckBox *CHKaspectratio;
    QCheckBox *CHKbasins;
    QCheckBox *CHKbonds;
    QCheckBox *CHKcenters;
    QCheckBox *CHKcenterslabel;
    QCheckBox *CHKcontourslabel;
    QCheckBox *CHKcps[max_cps];
    QCheckBox *CHKcurvelabel;
    QCheckBox *CHKgrid;
    QCheckBox *CHKmulticolor;
    QCheckBox *CHKprinthist;
    QCheckBox *CHKsettingsfile;
    QCheckBox *CHKsettingstitle;
    QCheckBox *CHKsettingsXlabel;
    QCheckBox *CHKsettingsYlabel;
    QCheckBox *CHKsettingsZoom;
    QCheckBox *CHKsettingsMargins;
    QCheckBox *CHKsettingscntPen;
    QCheckBox *CHKsettingselinesPen;
    QCheckBox *CHKsettingsbkg;
    QCheckBox *CHKshowparthist;
    QCheckBox *CHKshowfield;
    QCheckBox *CHKsmoothlines;
    QCheckBox *CHKtranspbg;
    QCheckBox *CHKtranspbgcapture;
    QCheckBox *CHKsghistlabel;
    QCheckBox *CHKXhide;
    QCheckBox *CHKXscalebottom;
    QCheckBox *CHKXscaletop;
    QCheckBox *CHKYscaleleft;
    QCheckBox *CHKYscaleright;
    QCheckBox *CHKticks;
    QCheckBox *CHKtitlehide;
    QCheckBox *CHKtitle;
    QCheckBox *CHKtitlebox;
    QCheckBox *CHKXlabel;
    QCheckBox *CHKXlabelbox;
    QCheckBox *CHKYhide;
    QCheckBox *CHKYlabel;
    QCheckBox *CHKYlabelbox;
    QCheckBox *CHKYhorizontal;
    QCheckBox *CHKYlabelhorizontal;
// 		QColor variables
    QColor *basinscolor;
    QColor *bkgcolor;
    QColor bondscolor;
    QColor centerscolor;
    QColor centerslabelcolor;
    QColor contourslabelcolor;
    QColor fonttitlecolor;
    QColor fontXlabelcolor;
    QColor fontYlabelcolor;
    QColor *gridcolor;
    QColor negativecontourcolor;
    QColor positivecontourcolor;
    QColor zerocontourcolor;
// 		QDialog variables
    QDialog *FRMtitle;
    QDialog *FRMXlabel;
    QDialog *FRMYlabel;
// 		QDoubleValidator variables
    QDoubleValidator *myDoubleValidator;
//		QFont variables
    QFont fonttitle;
    QFont fontXlabel;
    QFont fontYlabel;
// 		QGroup box variables
    QGroupBox *FRMarrows;
    QGroupBox *FRMaspectratio;
    QGroupBox *FRMbasins;
    QGroupBox *FRMbasinspen;
    QGroupBox *FRMbkgcolor;
    QGroupBox *FRMbonds;
    QGroupBox *FRMcaptura;    
    QGroupBox *FRMcenters;
    QGroupBox *FRMcntpen;
    QGroupBox *FRMcontours;
    QGroupBox *FRMcontourscolor;
    QGroupBox *FRMcontourslabel;
    QGroupBox *FRMcontournegstyle;
    QGroupBox *FRMcontourpostyle;
    QGroupBox *FRMcontourzerostyle;
    QGroupBox *FRMcontourstyle;
    QGroupBox *FRMcps;
    QGroupBox *FRMcurvespen;
    QGroupBox *FRMelinespen;
    QGroupBox *FRMfile;
    QGroupBox *FRMfontcurvelabels;
    QGroupBox *FRMfontsghistlabels;
    QGroupBox *FRMgridcolor;
    QGroupBox *FRMgridstyle;
    QGroupBox *FRMmargins;
    QGroupBox *FRMprinthist;
    QGroupBox *FRMradfactortype;
    QGroupBox *FRMscales;
    QGroupBox *FRMsettings;
    QGroupBox *FRMsghistpen;
    QGroupBox *FRMsmoothlines;
    QGroupBox *FRMticks;
    QGroupBox *FRMtitlegroup;
    QGroupBox *FRMXlabelgroup;
    QGroupBox *FRMYlabelgroup;
    QGroupBox *FRMzoomregion;
// 		QLabel variables 
    QLabel *LBLfile;
    QLabel *LBLBottomMargin;
    QLabel *LBLcontourslabel;
    QLabel *LBLinf;
    QLabel *LBLLeftMargin;
    QLabel *LBLnegativecontourcolor;
    QLabel *LBLbasins;
    QLabel *LBLbasinspenwidth;
    QLabel *LBLcntpenwidth;
    QLabel *LBLcontourfile;
    QLabel *LBLcps[max_cps];
    QLabel *LBLcpsballradius;
    QLabel *LBLcpsfile;
    QLabel *LBLcurvespenwidth;
    QLabel *LBLelinespenwidth;
    QLabel *LBLfield;
    QLabel *LBLfrads;
    QLabel *LBLpor;
    QLabel *LBLpositivecontourcolor;
    QLabel *LBLpostitle;
    QLabel *LBLposXlabel;
    QLabel *LBLposYlabel;
    QLabel *LBLprthist;
    QLabel *LBLradfactor;
    QLabel *LBLr2flm;
    QLabel *LBLr2l2flm;
    QLabel *LBLRightMargin;
    QLabel *LBLscaledef;
    QLabel *LBLsghist;
    QLabel *LBLsghistpenwidth;
    QLabel *LBLsmoothfactor;
    QLabel *LBLsup;
    QLabel *LBLtolerance;
    QLabel *LBLTopMargin;
    QLabel *LBLx;
    QLabel *LBLXcifras;
    QLabel *LBLXticks;
    QLabel *LBLxtitle;
    QLabel *LBLxXlabel;
    QLabel *LBLxYlabel;
    QLabel *LBLy;
    QLabel *LBLYcifras;
    QLabel *LBLYticks;
    QLabel *LBLytitle;
    QLabel *LBLyXlabel;
    QLabel *LBLyYlabel;
    QLabel *LBLYlabelhorizontal;
    QLabel *LBLzerocontourcolor;
    QLabel *LBLzerotolerance;
// 		QLineEdit variables
    QLineEdit *TXTbondthreshold;
    QLineEdit *TXTdistancethreshold;
    QLineEdit *TXTbasins;
    QLineEdit *TXTcontour;
    QLineEdit *TXTcps;
    QLineEdit *TXThsize;
    QLineEdit *TXTfield;
    QLineEdit *TXTfile;
    QLineEdit *TXTfrads;
    QLineEdit *TXTprinthist;
    QLineEdit *TXTscalesize;
    QLineEdit *TXTsghist;
    QLineEdit *TXTtitle;
    QLineEdit *TXTvsize;
    QLineEdit *TXTxinf;
    QLineEdit *TXTXlabel;
    QLineEdit *TXTYlabel;
    QLineEdit *TXTxsup;
    QLineEdit *TXTxtitle;
    QLineEdit *TXTxXlabel;
    QLineEdit *TXTxYlabel;
    QLineEdit *TXTyinf;
    QLineEdit *TXTysup;
    QLineEdit *TXTytitle;
    QLineEdit *TXTyXlabel;
    QLineEdit *TXTyYlabel;
// 		QPoint variables
    QPoint position;
// 		QPushButton variables
    QPushButton *BTNactualizarcontour;
    QPushButton *BTNcaptura;
    QPushButton *BTNedittitle;
    QPushButton *BTNeditXlabel;
    QPushButton *BTNeditYlabel;
    QPushButton *BTNexeczoomreg;
    QPushButton *BTNfontcenterslabel;
    QPushButton *BTNfontcontourslabel;
    QPushButton *BTNfontcurvelabels;
    QPushButton *BTNfontsghistlabels;
    QPushButton *BTNfonttitle;
    QPushButton *BTNfontXlabel;
    QPushButton *BTNfontYlabel;
    QPushButton *BTNprinthist;
    QPushButton *BTNretrieveSettings;
    QPushButton *BTNsaveSettings;
    QPushButton *BTNsmooth;
// 		QRadioButton variables
    QRadioButton *RBTdashed;
    QRadioButton *RBTnegsolid;
    QRadioButton *RBTnegdashed;
    QRadioButton *RBTnegdotted;
    QRadioButton *RBTnegdashdot;
    QRadioButton *RBTpossolid;
    QRadioButton *RBTposdashed;
    QRadioButton *RBTposdotted;
    QRadioButton *RBTposdashdot;
    QRadioButton *RBTradfactor;
    QRadioButton *RBTr2flm;
    QRadioButton *RBTr2l2flm;
    QRadioButton *RBTscaledef;
    QRadioButton *RBTscreendef;
    QRadioButton *RBTsolid;
    QRadioButton *RBTuserdef;
    QRadioButton *RBTzerosolid;
    QRadioButton *RBTzerodashed;
    QRadioButton *RBTzerodotted;
    QRadioButton *RBTzerodashdot;
// 		QScrollArea variables
    QScrollArea *graficaArea;
// 		QScrollArea variables
    QSize size;
// 		QSpin box variables
    QSpinBox *SPBarrowssep;
    QSpinBox *SPBarrowssize;
    QSpinBox *SPBarrowsskew;
    QSpinBox *SPBarrowswidth;
    QSpinBox *SPBbondswidth;
    QSpinBox *SPBBottomMargin;
    QSpinBox *SPBcpsballradius;
    QSpinBox *SPBshftcenterlabels;
    QSpinBox *SPBcenterradius;
    QSpinBox *SPBbasinspenwidth;
    QSpinBox *SPBcntpenwidth;
    QSpinBox *SPBcurvespenwidth;
    QSpinBox *SPBelinespenwidth;
    QSpinBox *SPBimagequality;
    QSpinBox *SPBLeftMargin;
    QSpinBox *SPBRightMargin;
    QSpinBox *SPBsghistpenwidth;
    QSpinBox *SPBsmoothfactor;
    QSpinBox *SPBtolerance;
    QSpinBox *SPBTopMargin;
    QSpinBox *SPBXcifras;
    QSpinBox *SPBXticks;
    QSpinBox *SPBYcifras;
    QSpinBox *SPBYticks;
    QSpinBox *SPBzerotolerance;
//		QString variables
    QString printhistfile;
    QString title;
    QString viewername;
    QString Xlabel;
    QString Ylabel;
// 		QToolBox variables
    QToolBox *toolBox;
// 		QToolButton variables
    QToolButton *BTNbasins;
    QToolButton *BTNcontour;
    QToolButton *BTNcps;
    QToolButton *BTNfield;
    QToolButton *BTNfile;
    QToolButton *BTNfrads;
    QToolButton *BTNprthist;
    QToolButton *BTNsghist;
// 		QWidget variables
    QWidget *page_contours;
    QWidget *page_field;
    QWidget *page_frad;
    QWidget *page_CPs;
    QWidget *page_basins;
    QWidget *page_options;
    QWidget *page_capture;
    QWidget *page_save;
    QWidget *page_sghistogram;
    QWidget *Wtablacont;
//              Table
    Sheet *SHTclist;
//        Bool functions
    bool readfieldfile(QString filename);
    bool checkplanesuffix(QString);
    bool readbasinsfile(QString filename);
    bool readcntfile(QString filename);
    bool readcpsfile(QString filename);
    bool readfradfile(QString filename);
    bool readhstfile(QString filename);
//        Int functions
    int get_plane_case(double, double, double);
//        QVector2D functions
    QVector2D xyzTouv(QVector3D);
// 		Void functions
    void initpointers();
    void plot_cam2D();
    void plot_basins();
    void plot_cps();
    void plot_frads();
    void plot_cnt();
    void plot_sghist();
    void setTitle();
};

#endif	/* VIEWER2D_H */

#ifndef MENUVIEWER2D_H
#define	MENUVIEWER2D_H

#include <QToolBox>

class QToolBox;

class menuViewer2D : public QToolBox
{
    Q_OBJECT

public:
    menuViewer2D(QWidget *parent) : QToolBox(parent){
	};
	QSize sizeHint() const
    {          
        return this->size();
    } 
};
#endif	/* MENUVIEWER2D_H */
