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
//  File:   mainwindow.h
//
//      Last version: October 2018
//
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QApplication>
#include <QButtonGroup>
#include <QComboBox>
#include <QDockWidget>
#include <QFileDialog>
#include <QList>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QProcess>
#include <QSettings>
#include <QSplashScreen>
#include <QStatusBar>
#include <QString>
#include <QTextEdit>
#include <QTextStream>
#include <QTableWidget>
#include <QtGlobal>
#include <QToolBar>
#include <QTranslator>
#include <QVector3D>
#include <stdlib.h>
#include <string>

#include <QtDebug>

#include "dialog.h"
#include "IniFile.h"
#include "math.h"
#include "viewer2D.h"
#include "Sheet.h"
#include "glWidget.h"
using namespace std;

#define MAX_LEXP 25
#define MAX_LEXPZJ 22
#define MAX_KEXPZJ 40
#define MAX_ARCHIVOS_RECIENTES  20
#define MAX_NUM_PROCESSORS 64
#define NFORTRANPROCS 13

class QAction;
class QActionGroup;
class QComboBox;
class QDir;
class QMenu;
class QTextEdit;
class QToolBox;
class QWidget;
class QLabel;
class QLineEdit;
class QGroupBox;
class QToolButton;
class QSpinBox;
class QRadioButton;
class QCheckBox;
class QTableWidget;
class QTableWidgetItem;
class QPushButton;
class QProcess;
class QTabWidget;
class QStringList;
class QPrinter;
class QPainter;
class Sheet;
//class Visor3D;
//class Visor2D;
class Surface2D;
class QDoubleValidator;

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class ViewerDialog : public QDialog
{
    Q_OBJECT
public:
    explicit ViewerDialog();
    ~ViewerDialog();

protected:
    virtual void reject();

signals:
    void closed();

};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    static void dmax(QVector<double> &v,double &max);
    static void dmin(QVector<double> &v,double &min);
    static void dminmax(QVector<double> &v,double &min,double &max);
    void finishsplash();
    
    
protected:
    void closeEvent(QCloseEvent *event);       
    
signals:
    void initproject();

public slots:

    void deleteplot(int); 
    void deletewidget(int);
    void external_package();
    void raiseplot(int);
    void raisewidget(int);
    void showplot(int);
    void showwidget(int);
    void exit();

private slots:      // Alphabetically sorted (including type in sort)
    bool end_Viewer2DDialog();
    bool end_Viewer3DDialog();
    bool GAUSS_two_pass_case(QString);
    bool saveProject();
    bool SaveProjectAs();

    void about();
    void addglWidget();
    void addviewer2D();
    
    void BTNdgfilelines_clicked();
    void BTNeffilelines_clicked();
    void BTNjob_clicked();
    void BTNpreview_clicked();

    void CHKatdensinput_changed(int state);
    void CHKatdensmpi_changed(int state);
    void CHKtopobasin_changed(int);
    void CHKdensder2_changed(int state);
    void CHKdensgrad_changed(int state);
    void CHKdensgrid_changed();
    void CHKdensinput_changed(int state);
    void CHKdensmpi_changed(int state);
    void CHKdensxyz_changed();
    void CHKdgextralines_changed();
    void CHKdginput_changed(int state);
    void CHKdgmpi_changed(int state);
    void CHKdgxyz_changed();
    void CHKdZJgrid_changed();
    void CHKdZJinput_changed(int state);
    void CHKdZJmpi_changed(int state);
    void CHKdZJxyz_changed();
    void CHKefextralines_changed();
    void CHKefinput_changed(int state);
    void CHKeflong_changed();
    void CHKefmpi_changed(int state);
    void CHKefxyz_changed();
    void CHKSGhexactMESP_changed(int);
    void CHKtopoexdraw_changed(int);
    void CHKfradextras_changed();
    void CHKfragments_changed(int);
    void CHKHFlatomsel_changed();
    void CHKpotlong_changed();
    void CHKtopomolgraph_changed(int);
    void CHKMOgrid_changed();
    void CHKMOinput_changed(int state);
    void CHKMOmpi_changed(int state);
    void CHKpotder2_changed(int state);
    void CHKpotexact_changed(int state);
    void CHKpotgrad_changed(int state);
    void CHKpotgrid_changed();
    void CHKpotinput_changed(int state);
    void CHKpotmpi_changed(int state);
    void CHKpotxyz_changed();
    void CHKSGholeinput_changed(int state);
    void CHKSGholempi_changed(int state);
    void CHKtopoaddguess_changed(int state);
    void CHKtopograph_changed();
    void CHKtopoinput_changed(int state);
    void CHKtopompi_changed(int state);
    void CHKtopoxyz_changed();
    void CHKZJinput_changed(int state);
    void CHKZJmpi_changed(int state);

    void chooseLanguage(QAction *action);
    void ChooseZJ_changed();

    void create_damproj(int exitCode, QProcess::ExitStatus exitStatus);
    void createLanguageMenu();

    void density_resolution_changed();
    void disable_pages();

    void execDam();
    void execDamden();
    void execDamdengrad();
    void execDamdenZJ();
    void execDamfield();
    void execDamforces();
    void execDamfrad();
    void execDamSGhole();
    void execDamZJ();
    void execFchk    ();
    void execGgbsDen();
    void execDammultrot();
    void execDamorb();
    void execDampot();
    void execDamTopography();
    void execImport();
    void execMolpro();
    void execMopac();
    void execNWChem();
    void execPsi4auxiliary();
    void execsgbs2sxyz(QString qstring);
    void execSxyzDen();
    void execTurbom();
    void execMOLEKEL();
    void external_geometry();

//    void gaussian_external();

    void Help();

    void importFile();
    void ImportFileNameMO();
    void ImportFileNameZJ();
    void importOUT();
    void importFileSGholeden();
    void ImportTopguessfilename();
//    void initviewers();

    void loadDefault_project();
    void loadDefault_atdens();
    void loadDefault_MED();
    void loadDefault_MESP();
    void loadDefault_SGhole();
    void loadDefault_TOPO();
    void loadDefault_HFforces();
    void loadDefault_Efield();
    void loadDefault_densgrad();
    void loadDefault_frad();
    void loadDefault_orimult();
    void loadDefault_MO();
    void loadDefault_ZJdens();
    void loadDefault_ZJtab();

    void menu_viewer2D();
    void menu_viewer3D();
    void moveviewertotop(int);
    void movewidgettotop(int);
    void MO_resolution_changed();
    void MOLPRO_two_pass_case(QString);

    void newProject();

    void openProject();
    void openRecentProjects();

    void page_project_layouts();
    void page_atdens_layouts();
    void page_MED_layouts();
    void page_MESP_layouts();
    void page_SGhole_layouts();
    void page_TOPO_layouts();
    void page_HFforces_layouts();
    void page_Efield_layouts();
    void page_densgrad_layouts();
    void page_frad_layouts();
    void page_orimult_layouts();
    void page_MO_layouts();
    void page_ZJdens_layouts();
    void page_ZJtab_layouts();

    void page_project_widgets();
    void page_atdens_widgets();
    void page_MED_widgets();
    void page_MESP_widgets();
    void page_SGhole_widgets();
    void page_TOPO_widgets();
    void page_HFforces_widgets();
    void page_Efield_widgets();
    void page_densgrad_widgets();
    void page_frad_widgets();
    void page_orimult_widgets();
    void page_MO_widgets();
    void page_ZJdens_widgets();
    void page_ZJtab_widgets();

    void potential_resolution_changed();

    void PrintFile();
    void PrintFilePdf();

    void processError(QProcess::ProcessError error);
    void processOutput(int exitCode, QProcess::ExitStatus exitStatus);
    void processStart();
    void processStop();

    void RBTdens2D3D_changed();
    void RBTdens2Dplanes_changed();
    void RBTdenslexact_changed();
    void RBTdensplane_changed();
    void RBTdenstype_changed();
    void RBTdg2D3D_changed();
    void RBTdg2Dplanes_changed();
    void RBTdZJ2Dplanes_changed();
    void RBTdZJplane_changed();
    void RBTef2D3D_changed();
    void RBTef2Dplanes_changed();
    void RBTlocal_changed();
    void RBTMO2D3D_changed();
    void RBTMO2Dplanes_changed();
    void RBTMOplane_changed();
    void RBTpot2D3D_changed();
    void RBTpot2Dplanes_changed();
    void RBTpotplane_changed();
    void RBTtopodensity_changed();
    void RBTdZJ2D3D_changed();

    void read_page_project(string file, QString ImportFolder);
    void read_page_atdens(string file);
    void read_page_MED(string file);
    void read_page_MESP(string file);
    void read_page_SGhole(string file);
    void read_page_TOPO(string file);
    void read_page_HFforces(string file);
    void read_page_Efield(string file);
    void read_page_densgrad(string file);
    void read_page_frad(string file);
    void read_page_orimult(string file);
    void read_page_MO(string file);
    void read_page_ZJdens(string file);
    void read_page_ZJtab(string file);

    void read_CHK(const char * a, const char * b, QCheckBox * c, const string file);
    void read_plot_dimension(const char * a, const char * b, const string file,
        QRadioButton * d2, QRadioButton * l, QRadioButton * m, QRadioButton * h);
    void read_double_to_TXT(const char * a, const char * b, QLineEdit * c, const string file);
    void read_text_to_TXT(const char * a, const char * b, QLineEdit * c, const string file);
    void read_RBT(const char * a, const char * b, QRadioButton * c, const string file);
    void read_SPB(const char * a, const char * b, QSpinBox * c, const string file);
    void read_resolution(const char * a1, const char * a2, const char * a3, const char * b,
        const string file, const QLineEdit * xi, const QLineEdit * xf,
        const QLineEdit * yi, const QLineEdit * yf, const QLineEdit * zi, const QLineEdit * zf,
        QRadioButton * d3, QRadioButton * l, QRadioButton * m, QRadioButton * h, QRadioButton * c);

    void saveOptions_0(string file, bool *printwarns, QString *warns);
    void saveOptions_1(string file, bool *printwarns, QString *warns);
    void saveOptions_2(string file, bool *printwarns, QString *warns);
    void saveOptions_3(string file, bool *printwarns, QString *warns);
    void saveOptions_4(string file, bool *printwarns, QString *warns);
    void saveOptions_5(string file, bool *printwarns, QString *warns);
    void saveOptions_6(string file, bool *printwarns, QString *warns);
    void saveOptions_7(string file, bool *printwarns, QString *warns);
    void saveOptions_8(string file, bool *printwarns, QString *warns);
    void saveOptions_9(string file, bool *printwarns, QString *warns);
    void saveOptions_10(string file, bool *printwarns, QString *warns);
    void saveOptions_11(string file, bool *printwarns, QString *warns);
    void saveOptions_12(string file, bool *printwarns, QString *warns);
    void saveOptions_13(string file, bool *printwarns, QString *warns);

    void SPBatdenslmaxexp_changed();
    void SPBatdensmpi_changed(int nprocessors);
    void SPBdenslmaxexp_changed();
    void SPBdenslminexp_changed();
    void SPBdensmpi_changed(int nprocessors);
    void SPBdgmpi_changed(int nprocessors);
    void SPBefmpi_changed(int nprocessors);
    void SPBfradltab_changed();
    void SPBfradmtab_changed();
    void SPBmrotleft_changed();
    void SPBmrotlmax_changed();
    void SPBmrotlmin_changed();
    void SPBmrotmiddle_changed();
    void SPBmrotright_changed();
    void SPBdZJmpi_changed(int nprocessors);
    void SPBMOmpi_changed(int nprocessors);
    void SPBpotmpi_changed(int nprocessors);
    void SPBSGholempi_changed(int nprocessors);
    void SPBtopompi_changed(int nprocessors);
    void SPBZJkmax_changed();
    void SPBZJlmax_changed();
    void SPBZJmpi_changed(int nprocessors);
        
    void start();
    
    void tabChanged(int index);
    void TXTdensatoms_changed();
    void TXTdgdlt0_changed();
    void TXTdZJchoose_changed();
    void TXTefdlt0_changed();
    void TXTextgeometry_changed();
    void TXTHFforcesatoms_changed();
    void TXTfradatoms_changed();
    void TXTImport_changed();
    void TXTMOchoose_changed();
    void TXTmpicommand_changed();
    void TXTmpiflags_changed();
    void TXTmrotorimultatoms_changed();
    void TXTProjectFolder_changed(const QString &cad);
    void TXTProjectName_changed(const QString &cad);
    void TXTImportSGholeden_changed();
    void TXTValidate_changed();
    void TXTZJrstar_changed();

    void update_dockright();
    void update_menu_viewer2D();
    void update_menu_viewer3D();
    void updatewindowsoverlay();
        
    void write_intervals(const char * a, const char * b, const char * c, const char * d, const string file, 
        QLineEdit * xi, QLineEdit * xs, bool * p, QString * w);
    void write_option(const char * a, const char * b, const QString c, const string file,
        bool * p, QString * w);

    void ZJ_resolution_changed();

    QString get_execName(QString, QString);
    QString get_python();
    
private:       // Grouped in blocks according to menu tabs. Block variables grouped in types. First block: general variables

    int mden;    // type of density matrix: 0: SCF density, 1: SCF spin density, 2: CI density, 3: CI spin 
    int densplanecase;
    int dZJplanecase;
    int MOplanecase;
    int plotsknt;
    int potplanecase;
    int topindex;
    int widgetsknt;

    double rmax;
    double rstar;
    double xmax;
    double xmin;
    double ymax;
    double ymin;
    double zmax;
    double zmin;
    
    QSplashScreen *splash;
    
    QTranslator appTranslator;
    QTranslator qtTranslator;
    
//    Actions
    
    QAction *Acc2Dplot;
    QAction *Acc3Dview;
    QAction *AccAbout;
    QAction *AccAboutQt;
    QAction *AccExit;
    QAction *AccExternal;
    QAction *AccHelp;
//    QAction *AccLanguage[10];
    QAction *AccNew;
    QAction *AccOpen;
    QAction *AccPdf;
    QAction *AccPrint;
    QAction *AccRecentFiles[MAX_ARCHIVOS_RECIENTES];
    QAction *AccSave;
    QAction *AccSaveAs;

    QActionGroup *languageActionGroup;

//    Dialogs

    ViewerDialog *QDLviewer2D;
    ViewerDialog *QDLwidget3D;

//    Labels

    QLabel *LBLlanguage;

//    Lists

    QList<QPushButton*> BTNraiseplotslist;             // Stores buttons for rising 2D plots
    QList<QPushButton*> BTNshowplotsslist;              // Stores buttons for hide/show 2D plots
    QList<QPushButton*> BTNraisewidgetslist;            // Stores buttons for rising 3D plots
    QList<QPushButton*> BTNshowwidgetslist;             // Stores buttons for hide/show 3D viewers
    QList<QMetaObject::Connection> connections2D;
    QList<QMetaObject::Connection> connections3D;
    QList<QMetaObject::Connection> connectionsext;
//    QList<QMetaObject::Connection> connectionsRightMenu;
    QList<glWidget*> *widgets;
    QList<Viewer2D*> *plots;

//    Menus

    QMenu *languageMenu;
    QMenu *FileMenu;
    QMenu *GraphicsMenu;
    QMenu *HelpMenu;

//   Push Buttons

    QPushButton *BTNlangstart;
    QPushButton *BTNexit2D;
//    QPushButton *BTNexit3D;
    QPushButton *BTNnewplot;
    QPushButton *BTNnewwidget;
        
//    Directories and files
    QString ImportFolder;    // path to the folder where the original input data of the project reside
    QString ImportFile;    // name of the file for data import (.fchk for Gaussian, .basis or .coor or .mos for TURBMOOLE, .mkl for MOLEKEL, .out for MOLPRO))
    QString ProjectFolder;    // path to the project folder 
    QString ProjectName;    // common name for all the files of the project
    QString DataFile; // Data input filenames for FORTRAN programs

//    Toolbars
    
    QToolBar *ToolBarFile;
    QToolBar *ToolBarHelp;
    
//    General purpose variables
    
    QString ArchivoActual;
    QStringList ArchivosRecientes;
    
    QPushButton *BTNstop[NFORTRANPROCS];
    QPushButton *BTNtexto[NFORTRANPROCS];
    
    bool changes;
    bool activebeware;
    bool lzdo;
    bool lvalence;

    Dialog *dialog;

    int executing;      // Number of process currently running
        
    QDialog *FRMlanguage;
    
    QString LanguagePath;
    bool lslater;        // true if slater calculation, false if gaussian calculation

    QProcess *myProcess;
    
    int natom;

    QPushButton *BTNraiseviewers;
    
    QWidget *dockwidget;
    QWidget *page_project;
    QWidget *page_atdens;
    QWidget *page_MED;
    QWidget *page_MESP;
    QWidget *page_SGhole;
    QWidget *page_TOPO;
    QWidget *page_HFforces;
    QWidget *page_Efield;
    QWidget *page_densgrad;
    QWidget *page_frad;
    QWidget *page_orimult;
    QWidget *page_MO;
    QWidget *page_ZJdens;
    QWidget *page_ZJtab;
    
    QString System;

    QTabWidget *TAWprincipal;
    QTextEdit *textEdit;

    QToolBox *toolBox;
    QWidget *rightBox;
    QVBoxLayout *rightBoxLayout;
    
//   page_project: Project
//   ---------------------

    QToolButton *BTNImport;
    QPushButton *BTNexecImport;
    
    QGroupBox *FRM2Dwidgets;
    QGroupBox *FRM3Dwidgets;
    QGroupBox *FRMmpioptions;
    QGroupBox *FRMplots;
    QGroupBox *FRMproject;
    QGroupBox *FRMviewers;

    QLabel *LBLImport;
    QLabel *LBLmpicommand;
    QLabel *LBLmpiflags;
    QLabel *LBLProjectFolder;
    QLabel *LBLProjectName;

    QLineEdit *TXTImport;
    QLineEdit *TXTmpicommand;
    QLineEdit *TXTmpiflags;
    QLineEdit *TXTProjectFolder;
    QLineEdit *TXTProjectName;    

    QVBoxLayout *page_projectLayout;
    
//   page_atdens: Atomic densities
//   -----------------------------
    
    QPushButton *BTNexecDam;
    
    QLabel *LBLatdensmpi;
    
    QCheckBox *CHKatdensmpi;
    QCheckBox *CHKatdensinput;
    
    QGroupBox *FRMatdensmpi;
    QGroupBox *FRMatdensinput;
    QGroupBox *FRMatdenslmaxdisp;
    QGroupBox *FRMatdenslmaxexp;
    QGroupBox *FRMatdenstype;
    QGroupBox *FRMfitthrs;

    QLabel *LBLfitthreshold;
    QLabel *LBLfradthreshold;
    
    QRadioButton *RBTatdensD1center;
    QRadioButton *RBTatdensD2center;
    QRadioButton *RBTatdensDtotal;
    
    QSpinBox *SPBatdensmpi;
    QSpinBox *SPBatdenslmaxexp;
    QSpinBox *SPBatdenslmaxdisp;
    QSpinBox *SPBfitthreshold;
    QSpinBox *SPBfradthreshold;

    QWidget *Wtable1;
    
//   page_MED: Density
//   -----------------
    
    QPushButton *BTNexecDamden;

    QCheckBox *CHKdensder2;
    QCheckBox *CHKdensgrad;
    QCheckBox *CHKdensgrid;
    QCheckBox *CHKdensinput;
    QCheckBox *CHKdenslaplacian;
    QCheckBox *CHKdenslatomics;
    QCheckBox *CHKdensldensacc;
    QCheckBox *CHKdenslmolec;
    QCheckBox *CHKdensmpi;
    QCheckBox *CHKdensxyz;
    
    QGroupBox *FRMdensatoms;
    QGroupBox *FRMdensdeformation;
    QGroupBox *FRMdensderivs;
    QGroupBox *FRMdensexpansion;
    QGroupBox *FRMdensdamdenfilename;
    QGroupBox *FRMdensgrid;
    QGroupBox *FRMdensgrid2D;
    QGroupBox *FRMdensgrid3D;
    QGroupBox *FRMdensgridres;
    QGroupBox *FRMdensgridtype;
    QGroupBox *FRMdensinput;
    QGroupBox *FRMdensdensity;
    QGroupBox *FRMdensfragments;
    QGroupBox *FRMdenslrange;
    QGroupBox *FRMdensmpi;
    QGroupBox *FRMdensplane2D;
    QGroupBox *FRMdensplaneABC;
    QGroupBox *FRMdensresol2D;
    QGroupBox *FRMdensresol3D;
    QGroupBox *FRMdenssurfpar;
    QGroupBox *FRMdenssurftype;
    QGroupBox *FRMdensxyz;
    
    QLabel *LBLdensatoms;
    QLabel *LBLdensmpi;
    QLabel *LBLdensuresol;
    QLabel *LBLdensvresol;
    QLabel *LBLdensxresol;
    QLabel *LBLdensyresol;
    QLabel *LBLdenszresol;
    QLabel *LBLdensinf3d;
    QLabel *LBLdensinf2d;
    QLabel *LBLdenslmaxexp;
    QLabel *LBLldensminexp;
    QLabel *LBLdenssup3d;
    QLabel *LBLdenssup2d;
    QLabel *LBLdensu;
    QLabel *LBLdensv;
    QLabel *LBLdensx;
    QLabel *LBLdensxformula2D;
    QLabel *LBLdensy;
    QLabel *LBLdensyformula2D;
    QLabel *LBLdensz;
    QLabel *LBLdenszformula2D;

    QRadioButton *RBTdens2D;
    QRadioButton *RBTdens3D;
    QRadioButton *RBTdensdeform;
    QRadioButton *RBTdensExact;
    QRadioButton *RBTdenslrange;
    QRadioButton *RBTdensplane;
    QRadioButton *RBTdensplaneABC;
    QRadioButton *RBTdensplaneXY;
    QRadioButton *RBTdensplaneXZ;
    QRadioButton *RBTdensplaneYZ;
    QRadioButton *RBTdensRep1;
    QRadioButton *RBTdensothersurf;
    QRadioButton *RBTdensfulldensity;
    QRadioButton *RBTdensrcustom;
    QRadioButton *RBTdensrhigh;
    QRadioButton *RBTdensrlow;
    QRadioButton *RBTdensrmedium;
    
    QSpinBox *SPBdensmpi;
    QSpinBox *SPBdensures;
    QSpinBox *SPBdensvres;
    QSpinBox *SPBdensxres;
    QSpinBox *SPBdensyres;
    QSpinBox *SPBdenszres;
    QSpinBox *SPBdenslmaxexp;
    QSpinBox *SPBdenslminexp;
    
    QLineEdit *TXTdensatoms;
    QLineEdit *TXTdensdamdenfilename;
    QLineEdit *TXTdensplaneA;
    QLineEdit *TXTdensplaneB;
    QLineEdit *TXTdensplaneC;
    QLineEdit *TXTdensuinf;
    QLineEdit *TXTdensusup;
    QLineEdit *TXTdensvinf;
    QLineEdit *TXTdensvsup;
    QLineEdit *TXTdensxformula2D;
    QLineEdit *TXTdensxinf;
    QLineEdit *TXTdensxsup;
    QLineEdit *TXTdensyformula2D;
    QLineEdit *TXTdensyinf;
    QLineEdit *TXTdensysup;
    QLineEdit *TXTdenszformula2D;
    QLineEdit *TXTdenszinf;
    QLineEdit *TXTdenszsup;

    QWidget *Wtableden;
    
    QStringList *denslist;
    
    QValidator *densvalidator;

    Sheet *SHTxyz;

//   page_MESP: Electrostatic potential
//   ----------------------------------
    
    QPushButton *BTNexecDampot;
    
    QCheckBox *CHKpotgrid;
    QCheckBox *CHKpotlong;
    QCheckBox *CHKpotmpi;
    QCheckBox *CHKpotder2;
    QCheckBox *CHKpotexact;
    QCheckBox *CHKpotgrad;
    QCheckBox *CHKpotinput;
    QCheckBox *CHKpotxyz;
    
    QGroupBox *FRMpotgdampotfilename;
    QGroupBox *FRMpotgrid;
    QGroupBox *FRMpotgrid2D;
    QGroupBox *FRMpotgrid3D;
    QGroupBox *FRMpotgridres;
    QGroupBox *FRMpotgridtype;
    QGroupBox *FRMpotlong;
    QGroupBox *FRMpotlmaxexp;
    QGroupBox *FRMpotderivs;
    QGroupBox *FRMpotinput;
    QGroupBox *FRMpotmpi;
    QGroupBox *FRMpotplane2D;
    QGroupBox *FRMpotplaneABC;
    QGroupBox *FRMpotresol2D;
    QGroupBox *FRMpotresol3D;
    QGroupBox *FRMpotsurfpar;
    QGroupBox *FRMpotsurftype;
    QGroupBox *FRMpotxyz;
    
    QLabel *LBLpotinf;
    QLabel *LBLpotinf2d;
    QLabel *LBLpotmpi;
    QLabel *LBLpoturesol;
    QLabel *LBLpotvresol;
    QLabel *LBLpotxresol;
    QLabel *LBLpotyresol;
    QLabel *LBLpotzresol;
    QLabel *LBLpotres;
    QLabel *LBLpotsup;
    QLabel *LBLpotsup2d;
    QLabel *LBLpotlongthreshold;
    QLabel *LBLpotu;
    QLabel *LBLpotv;
    QLabel *LBLpotx;
    QLabel *LBLpotxformula2D;
    QLabel *LBLpoty;
    QLabel *LBLpotyformula2D;
    QLabel *LBLpotz;
    QLabel *LBLpotzformula2D;
    
    QRadioButton *RBTpot2D;
    QRadioButton *RBTpot3D;
    QRadioButton *RBTpotothersurf;
    QRadioButton *RBTpotplane;
    QRadioButton *RBTpotplaneABC;
    QRadioButton *RBTpotplaneXY;
    QRadioButton *RBTpotplaneXZ;
    QRadioButton *RBTpotplaneYZ;
    QRadioButton *RBTpotrcustom;
    QRadioButton *RBTpotrhigh;
    QRadioButton *RBTpotrlow;
    QRadioButton *RBTpotrmedium;
    
    QSpinBox *SPBpotlmaxexp;
    QSpinBox *SPBpotmpi;
    QSpinBox *SPBpotures;
    QSpinBox *SPBpotvres;
    QSpinBox *SPBpotxres;
    QSpinBox *SPBpotyres;
    QSpinBox *SPBpotzres;
    QSpinBox *SPBpotlongthreshold;
    
    QLineEdit *TXTpotgdampotfilename;
    QLineEdit *TXTpotplaneA;
    QLineEdit *TXTpotplaneB;
    QLineEdit *TXTpotplaneC;
    QLineEdit *TXTpotuinf;
    QLineEdit *TXTpotusup;
    QLineEdit *TXTpotvinf;
    QLineEdit *TXTpotvsup;
    QLineEdit *TXTpotxformula2D;
    QLineEdit *TXTpotxinf;
    QLineEdit *TXTpotxsup;
    QLineEdit *TXTpotyformula2D;
    QLineEdit *TXTpotyinf;
    QLineEdit *TXTpotysup;
    QLineEdit *TXTpotzformula2D;
    QLineEdit *TXTpotzinf;
    QLineEdit *TXTpotzsup;

    QWidget *Wtablepot;

    Sheet *SHTpotxyz;

//   page_MO: Molecular orbitals
//   ---------------------------

    QPushButton *BTNexecDamorb;

    QToolButton *BTNMOImportFile;

    QCheckBox *CHKMOgrid;
    QCheckBox *CHKMOgrad;
    QCheckBox *CHKMOmpi;
    QCheckBox *CHKMOinput;

    QGroupBox *FRMMOchoose;
    QGroupBox *FRMMOgridtype;
    QGroupBox *FRMMOImportfile;
    QGroupBox *FRMMOfilename;
    QGroupBox *FRMMOderivs;
    QGroupBox *FRMMOsurfpar;
    QGroupBox *FRMMOgrid;
    QGroupBox *FRMMOgrid2D;
    QGroupBox *FRMMOgrid3D;
    QGroupBox *FRMMOgridres;
    QGroupBox *FRMMOinput;
    QGroupBox *FRMMOmpi;
    QGroupBox *FRMMOplane2D;
    QGroupBox *FRMMOplaneABC;
    QGroupBox *FRMMOsurftype;
    QGroupBox *FRMMOresol2D;
    QGroupBox *FRMMOresol3D;

    QLabel *LBLMOchoose;
    QLabel *LBLMOImportFile;
    QLabel *LBLMOinf;
    QLabel *LBLMOinf2d3;
    QLabel *LBLMOmpi;
    QLabel *LBLMOuresol;
    QLabel *LBLMOvresol;
    QLabel *LBLMOxresol;
    QLabel *LBLMOyresol;
    QLabel *LBLMOzresol;
    QLabel *LBLMOsup;
    QLabel *LBLMOsup2d;
    QLabel *LBLMOu;
    QLabel *LBLMOv;
    QLabel *LBLMOx;
    QLabel *LBLMOxformula2D;
    QLabel *LBLMOy;
    QLabel *LBLMOyformula2D;
    QLabel *LBLMOz;
    QLabel *LBLMOzformula2D;

    QLineEdit *TXTMOchoose;
    QLineEdit *TXTMOxformula2D;
    QLineEdit *TXTMOyformula2D;
    QLineEdit *TXTMOzformula2D;
    QLineEdit *TXTMOImportfile;
    QLineEdit *TXTMOfilename;
    QLineEdit *TXTMOplaneA;
    QLineEdit *TXTMOplaneB;
    QLineEdit *TXTMOplaneC;
    QLineEdit *TXTMOuinf;
    QLineEdit *TXTMOusup;
    QLineEdit *TXTMOvinf;
    QLineEdit *TXTMOvsup;
    QLineEdit *TXTMOxinf;
    QLineEdit *TXTMOxsup;
    QLineEdit *TXTMOyinf;
    QLineEdit *TXTMOysup;
    QLineEdit *TXTMOzinf;
    QLineEdit *TXTMOzsup;

    QRadioButton *RBTMO2D;
    QRadioButton *RBTMO3D;
    QRadioButton *RBTMOothersurf;
    QRadioButton *RBTMOplane;
    QRadioButton *RBTMOplaneABC;
    QRadioButton *RBTMOplaneXY;
    QRadioButton *RBTMOplaneXZ;
    QRadioButton *RBTMOplaneYZ;
    QRadioButton *RBTMOrcustom;
    QRadioButton *RBTMOrhigh;
    QRadioButton *RBTMOrlow;
    QRadioButton *RBTMOrmedium;

    QSpinBox *SPBMOmpi;
    QSpinBox *SPBMOures;
    QSpinBox *SPBMOvres;
    QSpinBox *SPBMOxres;
    QSpinBox *SPBMOyres;
    QSpinBox *SPBMOzres;

    QStringList *MOlist;

    QValidator *MOvalidator;

//   page_TOPO: Topography
//   ---------------------
    
    QPushButton *BTNexecDamtopo;
    QToolButton *BTNtopoguessfilename;
    
    QCheckBox *CHKtopobasin;
    QCheckBox *CHKtopoexdraw;
    QCheckBox *CHKtopomolgraph;
    QCheckBox *CHKtopoinput;
    QCheckBox *CHKtopompi;
    QCheckBox *CHKtopoaddguess;
    QCheckBox *CHKtopograph;
    QCheckBox *CHKtopoxyz;

    QGroupBox *FRMtopoboxsize;
    QGroupBox *FRMtopofilename;
    QGroupBox *FRMtopoguessfilename;
    QGroupBox *FRMtopograph;
    QGroupBox *FRMtopoinput;
    QGroupBox *FRMtopompi;
    QGroupBox *FRMtopotype;
    QGroupBox *FRMtopoxyz;
    
    QLabel *LBLtopoexln;
    QLabel *LBLtopofdisp;
    QLabel *LBLtopompi;
    QLabel *LBLtopoboxb;
    QLabel *LBLtopoboxg;
    QLabel *LBLtopoboxl;        
    QLabel *LBLtopoboxt;
    QLabel *LBLtopoggradthr;
    QLabel *LBLtopocnvg;
    QLabel *LBLtopolmaxi;
    QLabel *LBLtopostepszt;

    QLineEdit *TXTtopoexln;
    QLineEdit *TXTtopofdisp;
    QLineEdit *TXTtopoboxb;
    QLineEdit *TXTtopoboxg;
    QLineEdit *TXTtopoboxl;
    QLineEdit *TXTtopoboxt;        
    QLineEdit *TXTtopocnvg;
    QLineEdit *TXTtopofilename;
    QLineEdit *TXTtopoggradthr;
    QLineEdit *TXTtopoguessfilename;
    QLineEdit *TXTtopostepszt;

    QRadioButton *RBTtopodensity;
    QRadioButton *RBTtopopotential;
    
    QSpinBox *SPBtopolmaxi;
    QSpinBox *SPBtopompi;               
    
    QWidget *Wtabledentopo;
    QWidget *Wtable8;

    Sheet *SHTtopoxyz;

//   page_SGhole: MESP sigma hole
//   ----------------------------

    QPushButton *BTNexecDamSGhole;

    QCheckBox *CHKSGhexactMESP;
    QCheckBox *CHKSGholeinput;
    QCheckBox *CHKSGholempi;

    QGroupBox *FRMImportSGholeden;
    QGroupBox *FRMSGholefilename;
    QGroupBox *FRMSGholeinput;
    QGroupBox *FRMSGholelmaxexp;
    QGroupBox *FRMSGholempi;

    QLabel *LBLSGholelocalextrema;
    QLabel *LBLSGholelocalpower;
    QLabel *LBLSGholecontour;
    QLabel *LBLSGholegeomthreshold;
    QLabel *LBLSGholelmaxexp;
    QLabel *LBLSGholelongthreshold;
    QLabel *LBLSGholempi;

    QLineEdit *TXTImportSGholeden;
    QLineEdit *TXTSGholecontour;
    QLineEdit *TXTSGholefilename;

    QSpinBox *SPBSGholegeomthreshold;
    QSpinBox *SPBSGholelmaxexp;
    QSpinBox *SPBSGholelocalextrema;
    QSpinBox *SPBSGholelongthreshold;
    QSpinBox *SPBSGholempi;

    QToolButton *BTNImportSGholeden;

//   page_Efield: Electric field
//   ---------------------------
    
    QPushButton *BTNexecDamfield;
    QToolButton *BTNeffilelines;

    QCheckBox *CHKefinput;
    QCheckBox *CHKefmpi;
    QCheckBox *CHKefextralines;
    QCheckBox *CHKeflong;
    QCheckBox *CHKefxyz;

    QComboBox *CMBefdirectionset;

    QGroupBox *FRMeffilename;
    QGroupBox *FRMefinput;
    QGroupBox *FRMeflines;
    QGroupBox *FRMeflmaxexp;
    QGroupBox *FRMeflong;
    QGroupBox *FRMefmpi;
    QGroupBox *FRMefextralines;
    QGroupBox *FRMefplane2D;
    QGroupBox *FRMefplaneABC;
    QGroupBox *FRMefplot2D;
    QGroupBox *FRMefplot3D;
    QGroupBox *FRMefplottype;
    QGroupBox *FRMefuv;
    QGroupBox *FRMefxyz;

    QLabel *LBLefmpi;
    QLabel *LBLeffilelines;
    QLabel *LBLeflongthreshold;

    QRadioButton *RBTef2D;
    QRadioButton *RBTef3D;
    QRadioButton *RBTefplaneABC;
    QRadioButton *RBTefplaneXY;
    QRadioButton *RBTefplaneXZ;
    QRadioButton *RBTefplaneYZ;

    QSpinBox *SPBefmpi;
    QSpinBox *SPBeflmaxexp;
    QSpinBox *SPBeflongthreshold;

    QLineEdit *TXTefdlt0;
    QLineEdit *TXTeffilelines;
    QLineEdit *TXTeffilename;
    QLineEdit *TXTefnlinpernuc;
    QLineEdit *TXTefnumpnt;
    QLineEdit *TXTefplaneA;
    QLineEdit *TXTefplaneB;
    QLineEdit *TXTefplaneC;
    QLineEdit *TXTefuinf;
    QLineEdit *TXTefusup;
    QLineEdit *TXTefuvratio;
    QLineEdit *TXTefvinf;
    QLineEdit *TXTefvsup;
    QLineEdit *TXTefxinf;
    QLineEdit *TXTefxsup;
    QLineEdit *TXTefyinf;
    QLineEdit *TXTefysup;
    QLineEdit *TXTefzinf;
    QLineEdit *TXTefzsup;

    QWidget *Wtable2ef;
    QWidget *Wtabledenef;

    Sheet *SHTefuv;
    Sheet *SHTefxyz;


//   page_densgrad: Density gradient
//   -------------------------------

    QPushButton *BTNexecDamdengrad;
    QToolButton *BTNdgfilelines;

    QCheckBox *CHKdgextralines;
    QCheckBox *CHKdginput;
    QCheckBox *CHKdgmpi;
    QCheckBox *CHKdgxyz;

    QComboBox *CMBdgdirectionset;

    QGroupBox *FRMdgextralines;
    QGroupBox *FRMdgfilename;
    QGroupBox *FRMdginput;
    QGroupBox *FRMdglmaxexp;
    QGroupBox *FRMdglines;
    QGroupBox *FRMdglong;
    QGroupBox *FRMdgmpi;
    QGroupBox *FRMdgplane2D;
    QGroupBox *FRMdgplaneABC;
    QGroupBox *FRMdgplot2D;
    QGroupBox *FRMdgplot3D;
    QGroupBox *FRMdgplottype;
    QGroupBox *FRMdguv;
    QGroupBox *FRMdgxyz;

    QLabel *LBLdgfilelines;
    QLabel *LBLdglongthreshold;
    QLabel *LBLdgmpi;

    QRadioButton *RBTdg2D;
    QRadioButton *RBTdg3D;
    QRadioButton *RBTdgplaneABC;
    QRadioButton *RBTdgplaneXY;
    QRadioButton *RBTdgplaneXZ;
    QRadioButton *RBTdgplaneYZ;


    QSpinBox *SPBdglmaxexp;
    QSpinBox *SPBdglongthreshold;
    QSpinBox *SPBdgmpi;

    QLineEdit *TXTdgatomsforces;
    QLineEdit *TXTdgdlt0;
    QLineEdit *TXTdgfilelines;
    QLineEdit *TXTdgfilename;
    QLineEdit *TXTdgnlinpernuc;
    QLineEdit *TXTdgnumpnt;
    QLineEdit *TXTdgplaneA;
    QLineEdit *TXTdgplaneB;
    QLineEdit *TXTdgplaneC;
    QLineEdit *TXTdguinf;
    QLineEdit *TXTdgusup;
    QLineEdit *TXTdguvratio;
    QLineEdit *TXTdgvinf;
    QLineEdit *TXTdgvsup;
    QLineEdit *TXTdgxinf;
    QLineEdit *TXTdgxsup;
    QLineEdit *TXTdgyinf;
    QLineEdit *TXTdgysup;
    QLineEdit *TXTdgzinf;
    QLineEdit *TXTdgzsup;

    QWidget *Wtable2dg;
    QWidget *Wtabledendg;

    Sheet *SHTdguv;
    Sheet *SHTdgxyz;

//   page_HFforces: Hellmann-Feynman forces
//   --------------------------------------

    QPushButton *BTNexecDamforces;

    QCheckBox *CHKHFlatomsel;

    QGroupBox *FRMHFatomsforces;
    QGroupBox *FRMHFdensfragments1;
    QGroupBox *FRMHFgdamforcesfilename;

    QLineEdit *TXTHFforcesatoms;
    QLineEdit *TXTHFgdamforcesfilename;

    QStringList *HFforceslist;

    QValidator *HFforcesvalidator;

//   page_frad: Radial factors
//   -------------------------
    
    QPushButton *BTNexecDamfrad;
    
    QCheckBox *CHKfradextras;
    QCheckBox *CHKfradderiv1;
    QCheckBox *CHKfradderiv2;
    
    QGroupBox *FRMfradatoms;
    QGroupBox *FRMfraddamfilename;
    QGroupBox *FRMfradialfactors;
    QGroupBox *FRMfradios;
    QGroupBox *FRMfradderivadas;
    
    QLabel *LBLfradatoms;
    QLabel *LBLfradltr;
    QLabel *LBLfradltab;
    QLabel *LBLfradmtab;
    QLabel *LBLfradrfin;
    QLabel *LBLfradrini;

    Sheet *SHTfradrlist;

    QSpinBox *SPBfradltab;
    QSpinBox *SPBfradmtab;
    
    QLineEdit *TXTfradatoms;
    QLineEdit *TXTfraddamfilename;
    QLineEdit *TXTfradltr;
    QLineEdit *TXTfradrfin;
    QLineEdit *TXTfradrini;
    
    QStringList *fradlist;
    
    QValidator *fradvalidator;
    
    QWidget *Wtable5;
    
//   page_orimult: Oriented multipoles
//   ---------------------------------
    
    QGroupBox *FRMmrotorimultfilename;
    QGroupBox *FRMmrotmultipoles;
    QGroupBox *FRMmrotcenters;
    QGroupBox *FRMmrotatoms;
        
    QLabel *LBLmrotlmin;
    QLabel *LBLmrotlmax;
    QLabel *LBLmrotleft;
    QLabel *LBLmrotmiddle;
    QLabel *LBLmrotright;
    QLabel *LBLmrotorimultatoms;

    QLineEdit *TXTmrotorimultatoms;
    QLineEdit *TXTmrotorimultfilename;
                
    QPushButton *BTNexecDammultrot;

    int spbmrotleft;    // Used to store the previous value of SPBmrotleft (see SLOT SPBmrotleft_changed)
    int spbmrotmiddle;    // Used to store the previous value of SPBmrotmiddle (see SLOT SPBmrotmiddle_changed)
    int spbmrotright;    // Used to store the previous value of SPBmrotright (see SLOT SPBmrotright_changed)
        
    QSpinBox *SPBmrotleft;
    QSpinBox *SPBmrotlmax;
    QSpinBox *SPBmrotlmin;
    QSpinBox *SPBmrotmiddle;
    QSpinBox *SPBmrotright;
    
    QStringList *mrotorimultlist;
    
    QValidator *mrotorimultvalidator;

//   page_ZJdens: Zernike-Jacobi
//   ---------------------------
    
    QPushButton *BTNexecDamZJ;
    
    QCheckBox *CHKZJinput;
    QCheckBox *CHKZJmpi;
    
    QGroupBox *FRMZJinput;
    QGroupBox *FRMZJlength;
    QGroupBox *FRMZJmpi;
    QGroupBox *FRMZJnquad;
    QGroupBox *FRMZJrstar;
    QGroupBox *FRMZJrstartype;
    QGroupBox *FRMZJthreshold;
    QGroupBox *FRMZJtype;
    
    QLabel *LBLZJkmax;
    QLabel *LBLZJlmax;
    QLabel *LBLZJmpi;
    QLabel *LBLZJnquad;
    QLabel *LBLZJthrdist;
    QLabel *LBLZJthrmult;
    
    QLineEdit *TXTZJrstar;
        
    QRadioButton *RBTZJacobi;
    QRadioButton *RBTZJZernike;
    QRadioButton *RBTZJechelon;
    QRadioButton *RBTZJrstarabs;
    QRadioButton *RBTZJrstarrel;
    
    QSpinBox *SPBZJkmax;
    QSpinBox *SPBZJlmax;
    QSpinBox *SPBZJmpi;
    QSpinBox *SPBZJnquad;
    QSpinBox *SPBZJthrdist;
    QSpinBox *SPBZJthrmult;

//   page_ZJtab: One-center density from Zernike-Jacobi expansion
//   ------------------------------------------------------------

    bool ldZJjacobi;
    
    QPushButton *BTNexecDamdZJ;
    QToolButton *BTNdZJImportFile;
    
    QCheckBox *CHKdZJgrad;
    QCheckBox *CHKdZJgrid;
    QCheckBox *CHKdZJinput;
    QCheckBox *CHKdZJmpi;
    QCheckBox *CHKdZJxyz;

    QDockWidget *dockright;
    
    QGroupBox *FRMdZJderivs;
    QGroupBox *FRMdZJexpansion;
    QGroupBox *FRMdZJexplength;
    QGroupBox *FRMdZJgrid2D;
    QGroupBox *FRMdZJgrid3D;
    QGroupBox *FRMdZJgridres;
    QGroupBox *FRMdZJgrid;
    QGroupBox *FRMdZJgridtype;
    QGroupBox *FRMdZJfilename;
    QGroupBox *FRMdZJImportfile;
    QGroupBox *FRMdZJinput;
    QGroupBox *FRMdZJmpi;
    QGroupBox *FRMdZJplane2D;
    QGroupBox *FRMdZJplaneABC;
    QGroupBox *FRMdZJresol2D;
    QGroupBox *FRMdZJresol3D;
    QGroupBox *FRMdZJsurfpar;
    QGroupBox *FRMdZJsurftype;
    QGroupBox *FRMdZJxyz;
    
    QLabel *LBLdZJchoose;
    QLabel *LBLdZJImportFile;
    QLabel *LBLdZJinf;
    QLabel *LBLdZJinf2d;
    QLabel *LBLdZJkmax;
    QLabel *LBLdZJlmax;
    QLabel *LBLdZJlmin;
    QLabel *LBLdZJmpi;
    QLabel *LBLdZJsup;
    QLabel *LBLdZJsup2d;
    QLabel *LBLdZJuresol;
    QLabel *LBLdZJu;
    QLabel *LBLdZJvresol;
    QLabel *LBLdZJv;
    QLabel *LBLdZJxresol;
    QLabel *LBLdZJx;
    QLabel *LBLdZJxformula2D;
    QLabel *LBLdZJyresol;
    QLabel *LBLdZJy;
    QLabel *LBLdZJyformula2D;
    QLabel *LBLdZJzresol;
    QLabel *LBLdZJz;
    QLabel *LBLdZJzformula2D;
    
    QRadioButton *RBTdZJ2D;
    QRadioButton *RBTdZJ3D;
    QRadioButton *RBTdZJchooseall;
    QRadioButton *RBTdZJchoosek;
    QRadioButton *RBTdZJchoosel;
    QRadioButton *RBTdZJchooselk;
    QRadioButton *RBTdZJchooselkm;
    QRadioButton *RBTdZJothersurf;
    QRadioButton *RBTdZJplane;
    QRadioButton *RBTdZJplaneABC;
    QRadioButton *RBTdZJplaneXY;
    QRadioButton *RBTdZJplaneXZ;
    QRadioButton *RBTdZJplaneYZ;
    QRadioButton *RBTdZJrlow;
    QRadioButton *RBTdZJrmedium;
    QRadioButton *RBTdZJrhigh;
    QRadioButton *RBTdZJrcustom;
        
    QLineEdit *TXTdZJchoose;
    QLineEdit *TXTdZJfilename;
    QLineEdit *TXTdZJImportfile;
    QLineEdit *TXTdZJplaneA;
    QLineEdit *TXTdZJplaneB;
    QLineEdit *TXTdZJplaneC;
    QLineEdit *TXTdZJuinf;
    QLineEdit *TXTdZJusup;
    QLineEdit *TXTdZJvinf;
    QLineEdit *TXTdZJvsup;
    QLineEdit *TXTdZJxformula2D;
    QLineEdit *TXTdZJxinf;
    QLineEdit *TXTdZJxsup;
    QLineEdit *TXTdZJyformula2D;
    QLineEdit *TXTdZJyinf;
    QLineEdit *TXTdZJysup;
    QLineEdit *TXTdZJzformula2D;
    QLineEdit *TXTdZJzinf;
    QLineEdit *TXTdZJzsup;
    
    QSpinBox *SPBdZJkmax;
    QSpinBox *SPBdZJlmax;
    QSpinBox *SPBdZJlmin;
    QSpinBox *SPBdZJmpi;
    QSpinBox *SPBdZJures;
    QSpinBox *SPBdZJvres;
    QSpinBox *SPBdZJxres;
    QSpinBox *SPBdZJyres;
    QSpinBox *SPBdZJzres;
    
    QStringList *QSLdZJxyz;
    QStringList *ZJlist;
    
    QWidget *WtabledZJ;
    
    Sheet *SHTdZJxyz;

//  External packages
//  --------------------------------------------------------

    bool preview;
    
    QButtonGroup *QBGjobcommand;
    QButtonGroup *QBGrunmode;

    QDialog *QDLexternal;

    QGroupBox *FRMextproc;

    QLabel *LBLextproc;
    QLabel *LBLextpathremote;
    QLabel *LBLextworkdir;

    QLineEdit *TXTextgeometry = new QLineEdit();
    QLineEdit *TXTextpathremote;
    QLineEdit *TXTextworkdir;
    QLineEdit *TXTkeywords;

    QPushButton *BTNjob;
    QPushButton *BTNpreview;

    QRadioButton *RBTlocal;
    QRadioButton *RBTremote;
    QRadioButton *RBTPBS;
    QRadioButton *RBTSGE;
    QRadioButton *RBTSLURM;

    QSpinBox *SPBextproc;

    QString extgeomfile;
    QString extgeompath;

    QTextEdit *extextEdit;
    
//  Functions (alphabetically sorted including type in sort)
//  --------------------------------------------------------

        bool Open(const QString &fileName);
        bool compareIntegers(const QString& s1, const QString& s2);
        bool createDir(QString &fullPathName);
        bool mustSave();
        bool existsinp(QString fullinputName,int tab,int def, bool pregunta);
        bool Save(const QString &fileName);

        void initpointers();
        void UpdateRecentFiles();
        void CreateActions();
        void loadDefault(int all);
        void CreateLeftMenu();
        void CreateRightMenu();
        void CreateMenus();
        void CreateStatusBar();
        void CreateToolBars();
        void DAMDENdatafile(const QString &fullFileName, const QString projectDir, const QString projectName);
        void defineRanges();
        void SetCurrentFile(const QString &fileName,bool usar,bool modificado);
        void SetDir(const QString &carpeta,const QString &nombre);
        void setuvxyz();
        void GDAMdatafile(const QString &fullinputName, const QString projectDir, const QString projectName);
        void import(const QString &fileName);
        void inputdatafile(const char *suffix, const char *section, const QString &fullFileName,
            const QString projectDir, const QString projectName);
        void saveOptions(const QString &fullFileName,int clase);
        void readGeometry(int &natom,QVector<double> &x,QVector<double> &y,QVector<double> &z,QVector<int> &ncarga);
        void readOptions(const QString &fullFileName);
        void readSettings();
        void rename_density_cntfile();
        void rename_pot_cntfile();
        void set_natom(int);
        void writeSettings();

        int get_natom();
        int get_plane_case(double, double, double);
        int read_natom(QString fileName);

        double set_delta(const char * c, double ini, double fin);
        double set_deltaorb(const char * c, double ini, double fin);
        double set_deltapot(const char * c, double ini, double fin);
        double set_deltaZJ(const char * c, double ini, double fin);

        QByteArray ReadSectionOptions(const char *SectionName, QFile *FileName);

        QString FileWithoutExt(const QString &fullFileName);
        QString FileWithoutPath(const QString &fullFileName);
        QString Extension(const QString &fullFileName);
        QString Path(const QString &fullFileName);
        QString planesuffix(int);
        QString Who_executing(int caso);
        QString toQString(string v);

        QVector3D wu;
        QVector3D wv;

        string toString(QString qv);

};

#endif

/* 
 * Estracted from File:   mainmenu.h
 * Author: rafa
 *
 * Created on 13 de febrero de 2013, 8:51
 */

#ifndef MAINMENU_H
#define    MAINMENU_H

#include <QToolBox>

class QToolBox;

class mainmenu : public QToolBox
{
    Q_OBJECT

public:
    mainmenu(QWidget *parent) : QToolBox(parent){    
    };
    QSize sizeHint() const
    {          
        return size;
    } 
protected:
    void resizeEvent(QResizeEvent *event) override;
private:
    QSize size = QSize(400,2000);
};

#endif    /* MAINMENU_H */
