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
//	Header file of glWidget class
//  Description: glWidget defines a widget that manages all the elements required for
//  3D display of one or more molecules and properties related to them
//  such as MED, MESP, CPs, etc on a QOpenGLWindow
//
//	File:   glWidget.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <iostream>

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QDialog>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QFileSystemWatcher>
#include <QMainWindow>
#include <QPushButton>
#include <QQuaternion>
#include <QRadioButton>
#include <QStatusBar>
#include <QStringList>
#include <QTime>
#include <QtCore/qmath.h>

#include "ColorButton.h"

#include "IniFile.h"

#include "axes.h"
#include "ballscyls.h"
#include "glWindow.h"
#include "lights.h"
#include "measures.h"
#include "mespimizer.h"
#include "rotations.h"
#include "screenshot.h"
#include "surface.h"
#include "translations.h"
#include "viewport.h"
#include "elements.h"
#include "widgetsubclasses.h"

#define ANGSTROM_TO_BOHR 1.889725989
//#define BOHR_TO_ANGSTROM 0.529177249
//#define INIT_BOND_THRESHOLD 1.2
#define MAX_FRAMES 1000
#define MAX_INTERVAL 220.
#define MIN_INTERVAL 20.
#define NINTERPOL 10
#define INTERVAL_INI 100
#define INTERVAL_SCALE 200
//#define PERSPECTIVE_ANGLE 45.	// Angle for Perspective (in degrees)
//#define ZNEAR 0.001             // Znear for Perspective
//#define ZFAR 200.				// Zfar  for Perspective
#define Z_TRANS_INI -20.

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif


class MainWindow3DViewer : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow3DViewer(QWidget *parent = 0);
    ~MainWindow3DViewer();
signals:
    void hide_viewer();
protected:
    virtual void reject();
    virtual void closeEvent(QCloseEvent *event);
};


class MyDockWidget : public QDockWidget
{
    Q_OBJECT
public:
    explicit MyDockWidget(QWidget *parent = 0);
    ~MyDockWidget();
signals:
    void resizedialog(QSize);
protected:
    virtual void resizeEvent(QResizeEvent *event);

};

class moleculeDialog;

class glWidget : public QWidget
{
//    friend class moleculeDialog;
    Q_OBJECT
public:
    explicit glWidget(QString *name, QWidget *parent = 0);

    glWidget();
    ~glWidget();

    QColor getbkgcolor();

    measures *msrs;

    QList<molecule*> *molecules;

    bool addmolecule();
    bool addcluster(QString);
    bool loadmoleculargeometry();
    bool loadclustergeometry(QString);
    bool readggbs(QString);
    bool readxyz(QString);
    bool reloadxyz(QString filename);

    int getwindownumber();

    QQuaternion getrotation();

    QString getWindowName();

    QVector3D getrotationAxis();

    void create_axes_menu();
    void create_capture_menu();
    void create_lights_menu();
    void create_balls_cyl_menu();
    void create_molecules_menu_and_layouts();
    void create_optimize_cluster_menu();
    void create_save_geometry_menu();
    void create_save_geometry_layouts();
    void create_save_retrieve_menu();
    void create_save_retrieve_layouts();
    void create_translation_menu();
    void create_rotation_menu();
    void create_viewport_menu();

    void createWindowDialog(QDockWidget *);
    void deleteWindowDialog();

    void lower_mainwindow();
    void raise_mainwindow();
    void retrievemolecule();
    void savegeometry();
    void updateviewport();
    void setbkgcolor(QColor);  
    void set_ProjectFolder(QString);
    void set_ProjectName(QString);
    void setWindowName(QString);


protected:
    virtual void closeEvent(QCloseEvent *);

signals:
    void hideviewer();
    void moveToTop(int);
    void setangstromtranslation(bool);


public slots:
    bool isvisible();
    bool addclusterTowindow(QString);
    bool addemptyclusterTowindow(QString);
    QPoint get_position();
    QSize get_size();

    void activatemolecules();
    void addmoleculeTowindow();
    void animaterotation();
    void ballcyl_changed(); 
    void BTNcapture_clicked();
    void BTNmeasures_clicked();
    void BTNoptimizeCluster_clicked();
    void BTNballcyl_clicked();
    void BTNlights_clicked();
    void BTNretrieveSettings_clicked();
    void BTNrotation_clicked();
    void BTNsaveGeometry_clicked();
    void BTNsaveSettings_clicked();
    void BTNsettings_clicked();
    void BTNtranslation_clicked();
    void BTNviewport_clicked();
    void checkactivate(int);
    void CHKscaleradii_changed(bool);
    void CHKsettingsgrids_changed();
    void CHKsettingsmolecs_changed();
    void deletecluster();
    void deletemolecule(int);
    void emitmovetotop();
    void rotation_changed();
    void showmolecule(int);
    void set_position(QPoint);
    void setrotation(QQuaternion);
    void setrotationAxis(QVector3D);
    void set_size(QSize);
    void settranslation(QVector3D);
    void set_visible(bool);
    void updatedisplay();
    void updateGL();
    void updateWindowDialog(QDockWidget *);
    void viewport_changed();

private slots:
    QString get_execName(QString processname, QString subdir);
    void adjustQDLSize();
    void animate(bool);
    void BTNaddaxes_clicked();
    void capture();                 // Name for frames and film files
    void CHKrotate_changed(bool);
    void clusterfile_changed(QString);
    void displayClusterGeometry(QString);
    void endmakingmovie();
    void exec_mespimizer();
    void framesfile_changed(QString);
    void interpol_changed(int);
    void optimizecanvas_changed(bool);
    void optimizeselect_changed(bool);
    void processError(QProcess::ProcessError);
    void processOutput(int exitCode, QProcess::ExitStatus exitStatus);
    void processStart();
    void QDLmeasures_rejected();
    void recordoptim_changed(bool);
    void replay_mespimization(QString);
    void reset_mespimization(QString);
    void resetinterval();
    void send_systemnames();
    void setdisplayEPIC(bool);
    void setenergyfont(QFont);
    void setenergycolor(QColor);
    void setenergyprecision(int);
    void sethartree(bool);
    void set_recordfilename(QString);
    void speed_changed(int);
    void startrecording();
    void stoprecording();
    void template_changed(int);
    void toggleanimation(bool);

private: 
    void create_measures();
    void create_window();
    void initpointers();
    void initQDLpointers();
    void linearinterpolation();
    void quaterninterpolation();
    void resetQDLpointers();
    void pause(int);

    bool angstrom;
    bool animating;
    bool cluster_exists;
    bool displayEPIC;
    bool guestfromcanvas;
    bool hartree;
    bool recordoptim;
    bool onlyselcp;
    bool optimvisible;
    bool scaleradii;
    bool visible;

    Elements *elem;

    float deltaAngles;
    float dltinterval;
    float Ener0;
    float Ener1;
    float Enerinterp;
    float interval;
    float lightpower;
    float specularindex;
    float stepwheel;

    glWindow *window;

    int delay;
    int energyprecision;
    int gdockrewidth;
    int interpoints;
    int kntframes;
    int ninterpol;
    int numframes;
    int templateindex;
    int speed;
    int windownumber;

    DoubleSpinBox *SPBbondthres;

    MainWindow3DViewer *mainWin;

    myScrollArea *scrollArea;

    moleculeDialog *QDLwindow;

    QCheckBox *CHKsettingsballs;
    QCheckBox *CHKsettingsbkg;
    QCheckBox *CHKsettingsgrids;
    QCheckBox *CHKsettingsisosurf;
    QCheckBox *CHKsettingslights;
    QCheckBox *CHKsettingsmolecs;
    QCheckBox *CHKsettingspos;
    QCheckBox *CHKsettingsrot;
    QCheckBox *CHKsettingssurf;
    QCheckBox *CHKsettingsviewport;

    QColor ambientcolor;
    QColor bkgcolor;
    QColor dihedralplanescolor;
    QColor energycolor;
    QColor lightcolor;
    QColor specularcolor;

    MyDockWidget *gdock;

    QGridLayout *layoutmolecules;

    QFont energyfont;

    QGroupBox *FRMaxes;
    QGroupBox *FRMballcyl;
    QGroupBox *FRMcanvas;
    QGroupBox *FRMcapture;
    QGroupBox *FRMlights;
    QGroupBox *FRMmeasures;
    QGroupBox *FRMoptimizeCluster;
    QGroupBox *FRMrotation;
    QGroupBox *FRMsaveGeometry;
    QGroupBox *FRMsettings;
    QGroupBox *FRMtranslation;
    QGroupBox *FRMviewport;

    QList<QCheckBox*> activateboxes;
    QList<QMetaObject::Connection> connections;         // Stores connections to release in destructor
    QList<QMetaObject::Connection> msrs_connections;    // Stores measures connections to release in destructor
    QList<QPushButton*> BTNshowlist;                    // Stores buttons for hide/show molecules

    QPoint position;

    QProcess *myProcess;

    QPushButton *BTNaddaxes;
    QPushButton *BTNballcyl;
    QPushButton *BTNcapture;
    QPushButton *BTNlights;
    QPushButton *BTNmeasures;
    QPushButton *BTNoptimizeCluster;
    QPushButton *BTNretrieveSettings;
    QPushButton *BTNrotation;
    QPushButton *BTNsaveGeometry;
    QPushButton *BTNsaveSettings;
    QPushButton *BTNsettings;
    QPushButton *BTNtranslation;
    QPushButton *BTNviewport;

    QQuaternion rotation;                   // Quaternion for rotation

    QScrollArea *windowArea;

    QSize size;

    QString CaptureFolder;
    QString clustername;
    QString framesfile;
    QString fullmolname;                    // Auxiliar variable to read full molecule name (including path)
    QString molname;                        // Auxiliar variable to read molecule name
    QString molpath;                        // Auxiliar variable to read path to molecule home directory
    QString mespimizerpath;
    QString mespimizerfile;
    QString mespimirecordfilename;
    QString ProjectFolder;
    QString ProjectName;
    QString recordfilename;                 // Name for frames and film files
    QString *strprocess;
    QString windowname;                     // Window name

    QStringList *arguments;
    QStringList molecule_names;

    QTimer *timer;


    QVector<QVector3D > xyz;                // Cartesian coordinates of centers
    QVector<QVector3D > xyzinit;            // Initial cartesian coordinates of centers for cluster optimization
    QVector<int> indmols;
    QVector<int> znuc;                      // Atomic number of centers

    QVector3D lightPosition;                // Light position in world space
    QVector3D rotationAxis;                 // Rotation axis


    qreal ballradius;                       // Radius of atom spheres
    qreal cylradius;                        // Radius of bond cylinders
    qreal disthressq;                       // Threshold for bonding
    qreal fov;                              // Field of view (angle in degrees)
    qreal zFar;                             // Farthest distance for objects to be displayed
    qreal zNear;                            // Nearest distance for objects to be displayed

};

//-------------------------------------------------------------------------------------------------
//
//          Class moleculeDialog
//
//-------------------------------------------------------------------------------------------------

class moleculeDialog : public QDialog
{
    Q_OBJECT
public:
    explicit moleculeDialog(QWidget *parent = 0);
    ~moleculeDialog();

    axes *axesdialog;
    ballscyls *ballsandcylsdialog;
    lights *lightsdialog;
    mespimizer *mespimizerdialog;
    rotations *rotationsdialog;
    screenshot *scrshotdialog;
    translations *translationsdialog;
    viewport *viewportdialog;
signals:
    void closedialog();
    void resetQDLpointers();

protected:
    virtual void closeEvent(QCloseEvent *);

public slots:
    void renewwidth(QSize);

private:
    QList<QMetaObject::Connection> connections;
};


#endif // GLWIDGET_H
