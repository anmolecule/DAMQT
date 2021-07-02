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
//  File:   measures.h
//  Description: class measures manages geometry measures (distances, angles
//  and dihedral angles)in 3D viewer
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//
#ifndef MEASURES_H
#define MEASURES_H

#include <QCheckBox>
#include <QColorDialog>
#include <QComboBox>
#include <QDialog>
#include <QFileDialog>
#include <QFontDialog>
#include <QtGlobal>
#include <QGroupBox>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QVBoxLayout>

#include <QDebug>

#include "ColorButton.h"
#include "elements.h"
#include "widgetsubclasses.h"

#define ANGSTROM_TO_BOHR 1.889725989
#define BOHR_TO_ANGSTROM 0.529177248882

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

//--------------------------------------------------
//
//          class measuresInfoWindow
//
//--------------------------------------------------

class measuresInfoWindow : public QDialog
{
    Q_OBJECT

public:
    measuresInfoWindow(QStringList molecules, bool showdist, QString *distances, bool showang, QString *angles,
            bool showdihed, QString *dihedrals, QWidget *parent);
    ~measuresInfoWindow();
    void update_angles(QString*);
    void update_dihedrals(QString*);
    void update_distances(QString*);


protected:
    virtual void closeEvent(QCloseEvent *);

signals:
    void window_closed();

public slots:


private slots:
    void measuresInfoWindow_close();

private:
    QLabel *LBLmolecules;
    QLabel *LBLdistances;
    QLabel *LBLangles;
    QLabel *LBLdihedrals;

};

//--------------------------------------------------
//
//          class editMeasuresDialog
//
//--------------------------------------------------

class editMeasuresDialog : public QDialog
{
    Q_OBJECT
public:
    explicit editMeasuresDialog(QWidget *parent = 0);
    ~editMeasuresDialog();
signals:
    void closed();
protected:
    void closeEvent(QCloseEvent *event);
    virtual void reject();
};

// Class measures
//


//--------------------------------------------------
//
//          class measures
//
//--------------------------------------------------

class measures : public QWidget
{
    Q_OBJECT
public:
    explicit measures(QWidget *parent = 0);

    ~measures();

    bool create_QDLmeasures();
    bool getmeasureangles();
    bool getmeasuredihedrals();
    bool getmeasuredistances();
    bool getmeasurenone();
    bool QDLmeasures_isVisible();

    void append_systemnames(QString);
    void clear_systemnames();
    void set_ProjectFolder(QString);

    editMeasuresDialog *QDLmeasures;

signals:

    void activate_measure(bool);
    void angles_precision(int);
    void angles_type(int);
    void angles_width(int);
    void angstrom_changed(bool);
    void arc_radius(int);
    void color_angles(QColor);
    void color_dihedrals(QColor);
    void color_distances(QColor);
    void dihedralplanes_color(QColor);
    void dihedrals_precision(int);
    void dist_precision(int);
    void dist_transpbkg(bool);
    void dist_vshift(int);
    void draw_arcs(bool);
    void draw_lines(bool);
    void emit_update_angles();
    void emit_update_dihedrals();
    void emit_update_distances();
    void font_angles(QFont);
    void font_dihedrals(QFont);
    void font_distances(QFont);
    void get_systemnames();
    void lines_type(int);
    void lines_width(int);
    void make_dihedral_planes();
    void measure_angles(bool);
    void measure_dihedrals(bool);
    void measure_distances(bool);
    void QDLmeasures_rejected(); 
    void reset_angles();
    void reset_dihedrals();
    void reset_distances();
    void show_angles(bool);
    void show_dihedrals(bool);
    void show_distances(bool);
    void surftype_solid(bool);
    void updateGL();
    void updateRightMenu();

public slots:
    void close_measuresInfoWindow();
    void distance_remove(int);
    void update_angles(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_dihedrals(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_distances(QVector<centerData>*, QVector<QMatrix4x4>*);
    void update_measure_centers(QVector<QVector3D>*);
    void LBLangles_add(QVector<centerData>*, QVector<QMatrix4x4>*);
    void LBLdihedrals_add(QVector<centerData>*, QVector<QMatrix4x4>*);
    void LBLdistances_add(QVector<centerData>*, QVector<QMatrix4x4>*);
    void QDLmeasures_raise();
    void reset_all();
    void set_molecules(QStringList);

private slots:
    void BTNanglescolor_clicked();
    void BTNdihedralscolor_clicked();
    void BTNdistancescolor_clicked();
    void BTNfontangles_clicked();
    void BTNfontdihedrals_clicked();
    void BTNfontdistances_clicked();
    void BTNprint_clicked();
    void BTNresetangles_clicked();
    void BTNresetdihedrals_clicked();
    void BTNresetdistances_clicked();
    void BTNsurfcolor_clicked();
    void CHKdrawarcs_changed();
    void CHKdrawlines_changed();
    void CHKdisttranspbkg_changed();
    void CHKshowangles_changed();
    void CHKshowdihedrals_changed();
    void CHKshowdistances_changed();
    void CHKshowangwin_changed();
    void CHKshowdihedwin_changed();
    void CHKshowdistwin_changed();
    void CMBanglestype_changed(int);
    void CMBlinestype_changed(int);
    void emitupdateRightMenu();
    void RBTangstrom_changed();
    void RBTmeasureangle_changed();
    void RBTmeasuredihed_changed();
    void RBTmeasuredist_changed();
    void RBTmeasurenone_changed();
    void SLDopacity_released();
    void SPBanglesprecision_changed(int);
    void SPBarcradius_changed(int);
    void SPBdihedralsprecision_changed(int);
    void SPBdistprecision_changed(int);
    void SPBangleswidth_changed(int);
    void SPBlineswidth_changed(int);
    void SPBdistvshift_changed(int);
    void update_measuresInfo();

private:
    
    ColorButton *BTNanglescolor;
    ColorButton *BTNdihedralscolor;
    ColorButton *BTNdistancescolor;
    ColorButton *BTNsurfcolor;

    Elements *elem;

    measuresInfoWindow *measuresInfo;
    
    QCheckBox *CHKdisttranspbkg;
    QCheckBox *CHKdrawlines;
    QCheckBox *CHKdrawarcs;
    QCheckBox *CHKshowangles;
    QCheckBox *CHKshowdihedrals;
    QCheckBox *CHKshowdistances;
    QCheckBox *CHKshowangwin;
    QCheckBox *CHKshowdihedwin;
    QCheckBox *CHKshowdistwin;

    QColor anglescolor;
    QColor dihedralscolor;
    QColor dihedralplanescolor;
    QColor distancescolor;

    QComboBox *CMBanglestype;
    QComboBox *CMBlinestype;

    QFont anglesfont;
    QFont dihedralsfont;
    QFont distancesfont;
    
    QGroupBox *FRMmeasureangles;
    QGroupBox *FRMmeasuredihedrals;
    QGroupBox *FRMmeasuredistances;
    QGroupBox *FRMmeasuregeom;

    QLabel *LBLalpha;
    QLabel *LBLdistvshift;
    QLabel *LBLlastangles;
    QLabel *LBLlastdihedrals;
    QLabel *LBLlastdistances;

    QList<QMetaObject::Connection> connections;     // Stores connections to release in destructor
    QList<QMetaObject::Connection> connectionsInfo;     // Stores connections to release in destructor

    QPoint measuresInfopos;
    
    QPushButton *BTNfontangles;
    QPushButton *BTNfontdihedrals;
    QPushButton *BTNfontdistances;
    QPushButton *BTNprint;
    QPushButton *BTNresetangles;
    QPushButton *BTNresetdihedrals;
    QPushButton *BTNresetdistances;
    
    QRadioButton *RBTangstrom;
    QRadioButton *RBTbohr;
    QRadioButton *RBTmeasureangle;
    QRadioButton *RBTmeasuredihed;
    QRadioButton *RBTmeasuredist;
    QRadioButton *RBTmeasurenone;
    QRadioButton *RBTsolidsurf;
    QRadioButton *RBTwiresurf;
    
    QSlider *SLDopacity;

    QSpinBox *SPBanglesprecision;
    QSpinBox *SPBangleswidth;
    QSpinBox *SPBarcradius;
    QSpinBox *SPBdihedralsprecision;
    QSpinBox *SPBdistprecision;
    QSpinBox *SPBlineswidth;
    QSpinBox *SPBdistvshift;

    QString *anglestext;
    QString *dihedralstext;
    QString *distancestext;
    QString *printfilename;
    QString ProjectFolder;    // path to the project folder

    QStringList anglesprinttext;
    QStringList dihedralsprinttext;
    QStringList distancesprinttext;
    QStringList lastselectangles;
    QStringList lastselectdihedrals;
    QStringList lastselectdist;
    QStringList molecules;
    QStringList systemnames;

    QVector<centerData> *anglecenters;
    QVector<centerData> *dihedralcenters;
    QVector<centerData> *distancecenters;

    QVector<QMatrix4x4> *mtrsf;

    QVector<QStringList> selectangles;
    QVector<QStringList> selectdihedrals;
    QVector<QStringList> selectdist;
};

#endif // MEASURES_H
