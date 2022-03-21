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
//	Header file of molecule class
//  Description: molecule defines a molecule object with geometry, surfaces, critical points, etc
//
//	File:   molecule.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef MOLECULE_H
#define MOLECULE_H

#include <QApplication>
#include <QCheckBox>
#include <QDesktopWidget>
#include <QDoubleSpinBox>
#include <QFont>
#include <QLabel>
#include <QGridLayout>
#include <QObject>
#include <QOpenGLShaderProgram>
#include <QOpenGLWidget>
#include <QPushButton>
#include <QtGlobal>
#include <QTimer>
#include <QVector>
#include <QVector3D>
#include <QtCore/qmath.h>

#include "ColorButton.h"

#include "grid.h"
#include "criticalpoints.h"
#include "elements.h"
#include "fieldlines.h"
#include "forces.h"
#include "surface.h"
#include "widgetsubclasses.h"

#define pi           3.14159265358979323846  /* pi */
#define ANGSTROM_TO_BOHR 1.889725989
#define BOHR_TO_ANGSTROM 0.529177249
#define INIT_BOND_THRESHOLD 1.2
#define MAX_CPS 4		 // Number of types of critical points
#define MAX_FORCES 5		 // Number of types of force types
#define MAX_INTERVAL 220.
#define MIN_INTERVAL 20.
#define INTERVAL_INI 100
#define INTERVAL_SCALE 200
#define SCALE 0.01
#define MOLSCALEHEIGHT 20.
#define MOLSCALEARROWSHEIGHT 2.
#define Z_TRANS_INI -20.

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class editMoleculeDialog : public QDialog
{
    Q_OBJECT
public:
    explicit editMoleculeDialog(QWidget *parent = 0);
    ~editMoleculeDialog();
signals:
    void closed();
protected:
    void closeEvent(QCloseEvent *event);
    virtual void reject();
};

class molecule : public QWidget
{
    friend class geomProcessor;
    Q_OBJECT
public:
    explicit molecule(QWidget *parent = 0);

    ~molecule();

    bool getangstromcoor();
    bool getangstromcp();
    bool getdrawatomcoords();
    bool getdrawatomindices();
    bool getdrawatomsymbols();
    bool gethideatoms();
    bool gethidebonds();
    bool gethidehydrogens();
    bool getiscluster();
    bool getonlyatomactive();
    bool isactive();
    bool isatomactive(int);
    bool isaxeslabelsvisible();
    bool isaxesvisible();
    bool isvisible();
    bool loadgrid();
    bool loadsurf();
    bool retrievegrid(QString);
    bool retrievesurf(QString);

    double getz_trans_ini();

    float getstepwheel();

    int  getcoordprecision();
    int  getlabelsvshift();
    int  getnumatoms();
    int  getnumcps(int);
    int  getznuc(int);

    void setlabelsvshift(int);
    void set_ProjectFolder(QString);
    void set_ProjectName(QString);

    fieldlines *flines;
    forces *hfforces;
    criticalpoints *cps;

    QColor getfontcolor();
    QColor getXaxis_color();
    QColor getYaxis_color();
    QColor getZaxis_color();

    QFont getfont();
    QFont getfontaxeslabels();

    QList<grid*> *grids;
    QList<surface*> *surfaces;

    QPoint getparentposition();

    QQuaternion getrotation();

    QString getfullname();
    QString getname();
    QString getpath();

    QVector <QVector3D> getpositionaxeslabels();
    QVector3D gettranslation();
    QVector3D getrotationAxis();

    QVector <GLuint> getallindices();
    QVector <GLuint> getallindicesoffset();
    QVector <VertexNormalData> getallvertices();
    QVector <GLuint> allaxesindices;           // Indices of axes vertices
    QVector <GLuint> allaxesindicesoffset;     // Offsets of axes indices in objects
    QVector <VertexNormalData> allaxesvertices;  // Vertices of axes triangles(position, normal, color)
    QVector <GLuint> allindices;           // Indices of vertices in structure
    QVector <GLuint> allindicesoffset;     // Offsets of indices in objects (cylinders and spheres) in structure
    QVector <VertexNormalData> allvertices;  // Vertices of triangles in structure (position, normal, color)
    QVector <int> getcharges();
    QVector <QVector3D> getxyz();

    void createeditMoleculeDialog();
    void darken();
    void initatomactive(bool);
    void lighten();
    void loadallindices(QVector<GLuint>);
    void loadallindicesoffset(QVector<GLuint>);
    void loadallvertices(QVector<VertexNormalData>);
    void loadcharges(QVector <int>);
    void loadxyz(QVector <QVector3D>);
    void setactive(bool);
    void setallatomactive(bool);
    void setatomactive(int, bool);
    void setdrawatomcoords(bool);
    void setdrawatomindices(bool);
    void setdrawatomsymbols(bool);
    void setfont(QFont);
    void setfontaxeslabels(QFont);
    void setfontcolor(QColor);
    void setfullname(QString);
    void setiscluster(bool);
    void setname(QString);
    void setparentposition(QPoint);
    void setonlyatomactive(bool);
    void setpath(QString);
    void setrotation(QQuaternion);
    void setrotationAxis(QVector3D);
    void setrotationButtons();
    void setscaleradii(bool);
    void setstepwheel(float);
    void settranslationButtons();
    void settranslation(QVector3D);
    void setvisible(bool);
signals:
    void animate(bool);
    void updatedisplay();
    void updateGL();
    void updateRightMenu();

public slots:
//    bool QDLcriticalpoints_isVisible();
    bool QDLeditMolecule_isVisible();
    bool startanimate();
    bool stopanimate();

    togglingGroupBox *surfaces_editor(surface *);

    void addcriticalpoints();
    void addgrid();
    void addfield();
    void addforces();
    void addsurface();
    void animaterotation();
    void BTNaddaxes_clicked();
    void BTNaddcriticalpoints_clicked();
    void BTNaddfieldlines_clicked();
    void BTNaddforces_clicked();
    void BTNanimation_clicked();
    void BTNcpcolor_change();
    void BTNcpcolorfont_clicked();
    void BTNcpeigcolor_change();
    void BTNcplblfont_clicked();
    void BTNcpselectall_clicked();
    void BTNcpselectnone_clicked();
    void BTNflinescolor_clicked();
    void BTNfontcolor_clicked();
    void BTNfont_clicked();
    void BTNfontaxeslabels_clicked();
    void BTNforcescolor_changed();
    void BTNrotation_clicked();
    void BTNselectall_clicked();
    void BTNselectnone_clicked();
    void BTNtranslation_clicked();
    void BTNskeleton_clicked(); // Opens dialog for surface color
    void BTNsymbols_clicked();
    void BTNXaxiscolor_clicked();
    void BTNYaxiscolor_clicked();
    void BTNZaxiscolor_clicked();
    void CHKactiveonly_changed(int);
    void CHKcpactiveonly_changed(int);
    void CHKcpcoords_changed(int);
    void CHKcpeigvec_changed(int);
    void CHKcpindices_changed(int);
    void CHKcps_changed();
    void CHKcpsymbols_changed(int);
    void CHKcpvalues_changed(int);
    void CHKflines_changed();
    void CHKflinesarrows_changed();
    void CHKforces_changed();
    void CHKhideatoms_changed();
    void CHKhidebonds_changed();
    void CHKhidehydrogens_changed();
    void CHKrotate_changed();
    void CHKshowaxes_changed(int);
    void CHKshowaxeslabels_changed(int);
    void CHKshowcoords_changed(int);
    void CHKshowindices_changed(int);
    void CHKshowsymbols_changed(int);
    void closeisosurfeditors();
    void deletegrid(int);
    void deletesurf(int);
    void editmolecule();
    void emitupdatedisplay();
    void emitupdateGL();
    void emitupdateRightMenu();
    void makeStructureBondsandSticks();
    void QDLeditMolecule_close();
    void QDLeditMolecule_delete();
    void QDLeditMolecule_raise();
    void RBTbohr_changed();
    void RBTbohrcoor_changed();
    void RBTbohrcp_changed();
    void readcpsfiles_dialog();
    void readfieldlines_dialog();
    void readforcefiles_dialog();
    void resetinterval();
    void reset_rotation();
    void reset_translation();
    void resizeQDLeditMolecule();
    void rotation_changed();
    void setballradius(qreal);
    void setcylradius(qreal);
    void setdisthressq(qreal);
    void setworld_rotation(QQuaternion);
    void SPBaxesarrowssize_changed(int);
    void SPBaxesarrowswidth_changed(int);
    void SPBaxeslength_changed(int);
    void SPBaxesthickness_changed(int);
    void SPBcoordprecision_changed(int);
    void SPBcpballradius_changed(int);
    void SPBcpcoordprecision_changed(int);
    void SPBcpeigarrowsize_changed(int);
    void SPBcpeigarrowwidth_changed(int);
    void SPBcpeiglength_change(int);
    void SPBcpeigthickness_changed(int);
    void SPBcpprecision_changed(int);
    void SPBcpvshift_changed(int);
    void SPBflinesarrowssep_changed(int);
    void SPBflinesarrowssize_changed(int);
    void SPBflinesarrowswidth_changed(int);
    void SPBflineslinewidth_changed(int);
    void SPBforcesarrowlength_changed(int);
    void SPBforcesarrowwidth_changed(int);
    void SPBforceslength_changed(int);
    void SPBforcesthickness_changed(int);
    void SPBlabelsvshift_changed(int);
    void SPBstepwheel_changed();
    void toggleactive();
    void toggleshowsurf(int);
    void translation_changed();
    void TXTcps_changed();
    void TXTfieldlines_changed();
    void TXTforces_changed();
    void updateeditMoleculeDialog();
    void updateQDLeditMolecule(bool);

private slots:
    void create_axes_widgets_and_layouts();
    void create_critical_points_widgets_and_layouts();
    void create_field_lines_widgets_and_layouts();
    void create_forces_widgets_and_layouts();
    void create_grids_widgets_and_layouts();
    void create_molecular_skeleton_widgets_and_layouts();
    void create_rotation_widgets_and_layouts();
    void create_surfaces_widgets_and_layouts();
    void create_symbols_indices_widgets_and_layouts();
    void create_translation_widgets_and_layouts();
    void make_axes();
    void makeAxesCone(int slices, int stacks, float height);    // Creates a triangle strip for a cone of width 1 and given height
    void makeAxesCylinder(int slices, int stacks);              // Creates a triangle strip for cylinder of radius 1 and height 1
    void makeCylinder(int slices, int stacks, double radius, double height);
    void makeCylinder(int slices, int stacks, double radius, double height, QVector3D center);
    void makeCylinder(int slices, int stacks, double radius, double height, QVector3D center, QVector3D axis);
    void makeCylinder(int slices, int stacks, double radius, double height, QVector3D center, QVector3D axis, QColor color);
    void makeSphere(int slices, int stacks, double radius);
    void makeSphere(int slices, int stacks, double radius, QVector3D center);
    void makeSphere(int slices, int stacks, double radius, QColor color);
    void makeSphere(int slices, int stacks, double radius, QVector3D center, QColor color);


private:
    int getscalevalueInt(float, int, int, float, float);

    QVector <GLuint> coneindices;
    QVector <GLuint> cylinderindices;
    QVector <QVector3D> conevertices;
    QVector <QVector3D> conenormals;
    QVector <QVector3D> cylindervertices;

    bool active;                          // If true, the molecule is active for transformations
    bool angstrom;                        // If true translation distances displayed in angstrom
    bool angstromcoor;                    // If true atoms coordinates displayed in angstrom
    bool angstromcp;                      // If true CPs coordinates displayed in angstrom
    bool axes_visible;
    bool axeslabels_visible;
    bool cpschecked[MAX_CPS];
    bool drawatomcoords;                  // If true displays atom coordinates
    bool drawatomindices;                 // If true displays atom indices
    bool drawatomsymbols;                 // If true displays atom symbols
    bool hideatoms;                       // If true, does not display atoms
    bool hidebonds;                       // If true, does not display bonds
    bool hidehydrogens;                   // If true, does not display hydrogen atoms
    bool iscluster;
    bool startanimation;                  // If true animates rotation  
    bool onlyatomactive;                  // If true atom labels displayed only for active atoms
    bool rotatex;
    bool rotatey;
    bool rotatez;
    bool scaleradii;
    bool visible;                         // If true molecule is displayed

    Elements *elem;

    float deltaAngles;
    float dltinterval;
    float interval;
    float stepwheel;

    int axesarrowssize;
    int axesarrowswidth;
    int axeslength;
    int axesthickness;
    int coordprecision;
    int labelsvshift;

    ColorButton *BTNcpcolor[MAX_CPS];
    ColorButton *BTNcpcolorfont;
    ColorButton *BTNcpeigcolor[3];
    ColorButton *BTNcpselectall;
    ColorButton *BTNcpselectnone;
    ColorButton *BTNflinescolor;          // Opens dialog for lines color
    ColorButton *BTNfontcolor;            // Opens dialog for font color
    ColorButton *BTNforcecolors[MAX_FORCES];
    ColorButton *BTNselectall;            // Marks all centers as active for indices display
    ColorButton *BTNselectnone;           // Marks all centers as nonactive for indices display
    ColorButton *BTNXaxiscolor;           // Dialog for X axis color
    ColorButton *BTNYaxiscolor;           // Dialog for Y axis color
    ColorButton *BTNZaxiscolor;           // Dialog for Z axis color

    DoubleSpinBox *SPBrot_angle;         // Rotation angle
    DoubleSpinBox *SPBrot_x;             // x component of rotation axis
    DoubleSpinBox *SPBrot_y;             // y component of rotation axis
    DoubleSpinBox *SPBrot_z;             // z component of rotation axis
    DoubleSpinBox *SPBstepwheel;         // stride for translation with mouse wheel
    DoubleSpinBox *SPBtras_x;            // x component of translation vector
    DoubleSpinBox *SPBtras_y;            // y component of translation vector
    DoubleSpinBox *SPBtras_z;            // z component of translation vector

    editMoleculeDialog *QDLeditMolecule;

    myScrollArea *scrollArea;

    QCheckBox *CHKactiveonly;
    QCheckBox *CHKcpactiveonly;
    QCheckBox *CHKcpcoords;
    QCheckBox *CHKcpeigvec;
    QCheckBox *CHKcpindices;
    QCheckBox *CHKcps[MAX_CPS];
    QCheckBox *CHKcpsymbols;
    QCheckBox *CHKcpvalues;
    QCheckBox *CHKflines;
    QCheckBox *CHKflinesarrows;
    QCheckBox *CHKforces[MAX_FORCES];
    QCheckBox *CHKhideatoms;
    QCheckBox *CHKhidebonds;
    QCheckBox *CHKhidehydrogens;
    QCheckBox *CHKrotatex;
    QCheckBox *CHKrotatey;
    QCheckBox *CHKrotatez;
    QCheckBox *CHKshowaxes;
    QCheckBox *CHKshowaxeslabels;
    QCheckBox *CHKshowcoords;
    QCheckBox *CHKshowindices;
    QCheckBox *CHKshowsymbols;

    QColor fontcolor;
    QColor Xaxis_color;
    QColor Yaxis_color;
    QColor Zaxis_color;


    QFont font;
    QFont fontaxeslabels;

    QGridLayout *layoutgrids;

    QGroupBox *FRMaxes;
    QGroupBox *FRMcoorunits;
    QGroupBox *FRMcps;
    QGroupBox *FRMcpeigvec;
    QGroupBox *FRMcpunits;
    QGroupBox *FRMcriticalpoints;
    QGroupBox *FRMfield;
    QGroupBox *FRMflinesarrows;
    QGroupBox *FRMforces;
    QGroupBox *FRMrotation;
    QGroupBox *FRMskeleton;
    QGroupBox *FRMsymbols;
    QGroupBox *FRMtranslation;
    QGroupBox *FRMtranslationunits;

    QLabel *LBLcoordprecision;
    QLabel *LBLcpcoordprecision;
    QLabel *LBLcpprecision;
    QLabel *LBLcpselect;
    QLabel *LBLloadinggrid;
    QLabel *LBLselect;

    QLineEdit *TXTcps;
    QLineEdit *TXTfieldlines;
    QLineEdit *TXTforces;

    QList<QMetaObject::Connection> connections;

    QPoint parentposition;
    QPoint scrollAreaposition;

    QPushButton *BTNaddaxes;
    QPushButton *BTNaddcriticalpoints;
    QPushButton *BTNaddfieldlines;
    QPushButton *BTNaddforces;
    QPushButton *BTNaddgrid;
    QPushButton *BTNaddsurface;
    QPushButton *BTNanimation;
    QPushButton *BTNcplblfont;
    QPushButton *BTNfont;
    QPushButton *BTNfontaxeslabels;
    QPushButton *BTNhide;
    QPushButton *BTNrotation;
    QPushButton *BTNskeleton;
    QPushButton *BTNsymbols;
    QPushButton *BTNtranslation;

    QQuaternion rotation;                 // Quaternion for rotation
    QQuaternion world_rotation;

    QRadioButton *RBTangstrom;
    QRadioButton *RBTangstromcoor;
    QRadioButton *RBTangstromcp;
    QRadioButton *RBTbohr;
    QRadioButton *RBTbohrcoor;
    QRadioButton *RBTbohrcp;

    QSlider *SLDspeed;                    // Speed of rotation animation

    QSpinBox *SPBaxesarrowsize;
    QSpinBox *SPBaxesarrowwidth;
    QSpinBox *SPBaxeslength;
    QSpinBox *SPBaxesthickness;
    QSpinBox *SPBcoordprecision;
    QSpinBox *SPBcpballradius;
    QSpinBox *SPBcpcoordprecision;
    QSpinBox *SPBcpeigarrowsize;
    QSpinBox *SPBcpeigarrowwidth;
    QSpinBox *SPBcpeigthickness;
    QSpinBox *SPBcpeiglength;
    QSpinBox *SPBcpprecision;
    QSpinBox *SPBcpvshift;
    QSpinBox *SPBflinesarrowssep;             // Arrows separation
    QSpinBox *SPBflineslinewidth;             // Lines width
    QSpinBox *SPBflinesarrowssize;            // Arrows size
    QSpinBox *SPBflinesarrowswidth;            // Arrows size
    QSpinBox *SPBforcesarrowlength;
    QSpinBox *SPBforcesarrowwidth;
    QSpinBox *SPBforceslength;
    QSpinBox *SPBforcesthickness;
    QSpinBox *SPBlabelsvshift;

    QString fullname;                     // Full geometry file name including path
    QString name;                         // Name for window
    QString path;                         // Path to molecule home directory (that which contains the file with geometry)
    QString ProjectFolder;
    QString ProjectName;

    QTimer *timer;

    QToolButton *BTNcps;
    QToolButton *BTNfieldlines;
    QToolButton *BTNforces;

    QVBoxLayout *layoutsurfs;

    QVector <QVector3D> positionaxeslabels;
    QVector3D rotationAxis;               // Rotation axis
    QVector3D translation;                // Translation vector

    QVector4D darkenshift;

    QVector<bool> atomactive;   // Atom active for visualization
    QVector<int> znuc;          // Atomic number of centers
    QVector<QVector3D > xyz;    // Cartesian coordinates

    qreal ballradius;                       // Radius of atom spheres
    qreal cylradius;                        // Radius of bond cylinders
    qreal disthressq;                       // Threshold for bonding

};

#endif // MOLECULE_H
