#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include <QObject>
#include <QWidget>
#include <QCheckBox>
#include <QDialog>
#include <QDoubleSpinBox>
#include <QDoubleValidator>
#include <QGroupBox>
#include <QKeyEvent>
#include <QLabel>
#include <QLineEdit>
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QPoint>
#include <QProcess>
#include <QRadioButton>
#include <QSlider>
#include <QSpinBox>
#include <QToolButton>
#include <QVector>
#include <QVector3D>
#include <QVector4D>

#include "ColorButton.h"
#include "widgetsubclasses.h"
#include "VertexNormalData.h"

# define pi           3.14159265358979323846    // pi
# define LOGBASE      1.08                      // base for logaritmics scale

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class editIsoSurfaceDialog : public QDialog
{
    Q_OBJECT
public:
    explicit editIsoSurfaceDialog(QWidget *parent = 0);
    ~editIsoSurfaceDialog();
signals:
    void closed();
protected:
    void closeEvent(QCloseEvent *event);
    virtual void reject();
};

// Class surface
//
class isosurface : public QWidget
{
    friend class geomProcessor;
    Q_OBJECT
public:
    explicit isosurface(QWidget *parent = 0);
    ~isosurface();

    bool isvisible();
    bool getnormalgrad();
    bool getshowgridbound();
    bool gettranslucence();
    bool getsolidsurf();

    int getscalevalueInt(float, int, int, float, float, bool);

    float getcontourvalue();
    float getmaxcontourvalue();
    float getmincontourvalue();
    float getopacity();
    float getscalevalueFloat(int, int, int, float, float, bool);

    QColor getsurfcolor();

    QPoint getinitialposition();

    QString get_basename();
    QString get_execName(QString processname, QString subdir);
    QString getfullname();
    QString getname();

    QVector3D gettranslation();
    QVector3D getrotationAxis();

    QVector <GLuint> getallindices();
    QVector <VertexNormalData> getallvertices();
    QVector <GLuint> allindices;             // Indices of vertices in surface
    QVector <VertexNormalData> allvertices;  // Vertices of triangles in surface (position, normal, color)
    QVector <GLuint> gridindices;             // indices of vertices in grid boundaries
    QVector <GLuint> gridindicesoffset;       // Offsets of indices in grid boundaries
    QVector <VertexNormalData> gridvertices;  // vertices of triangles in grid boundaries (position, normal, color)

    void generategridbounds(float *);
    void setcompatderiv(bool);
    void setcontourvalue(float);
    void setinitialposition(QPoint);
    void setmaxcontourvalue(float);
    void setmincontourvalue(float);
    void setfullname(QString);
    void setname(QString);
    void setopacity(float);
    void set_ProjectFolder(QString);
    void set_ProjectName(QString);
    void setsolidsurf(bool);
    void setsurfcolor(QColor);
    void settranslucence(bool);
    void setvisible(bool);

signals:
    void generatesurface();
    void opendialog();
    void updatedisplay();
    void updatelabelcolor(QString);
    void updateRightMenu();
    void updatetext(QString);

public slots:
//    bool QDLeditSurface_isVisible();

    togglingGroupBox * editisosurface();

    void toggleshowsurf();
    void BTNexec_clicked();
    void BTNsurfcolor_clicked();
    void CHKmpi_changed(int state);
    void CHKnormalgrad_changed();
    void CHKshowgrid_changed();
    void CHKtranslucence_changed();
    void closeeditor();
    void emitupdateRightMenu();
    void processError(QProcess::ProcessError error);
    void processOutput(int exitCode, QProcess::ExitStatus exitStatus);
    void processStart();
    void processStop();
    void RBTscale_changed();
    void RBTsurftype_changed();
    void SLDcontourvalue_changed(int);
    void SLDcontourvalue_released();
    void SLDopacity_changed(int);
    void SLDopacity_released();
    void SPBmpi_changed(int nprocessors);
    void SPBsensitive_changed(int);
    void TXTcontourvalue_changed();

private:
    bool active;                        // If true, the surface is active for transformations
    bool compatderiv;
    bool editoropen;
    bool isdensity;
    bool normalgrad;                    // If ture, computes normals from inerpolated gradient
    bool solidsurf;                     // If true, surface display is solid, if false, surface display is wire frame
    bool showgridbound;                 // If true displays grid boundaries
    bool translucence;                  // If true, translucence correction is applied
    bool visible;                       // If true the surface is displayed

    int nindices;
    int nprocessors;
    int nvertices;

    QVector<float> griddimensions;      // Original grid dimensions: xmin, xmax, ymin, ymax, zmin, zmax
    QVector<int> gridnxyz;              // Original grid number of points (nx, ny, nz)

    ColorButton *BTNsurfcolor;            // Opens dialog for surface color

    float contourvalue;                   // Function value for isosurface
    float maxcontourvalue;                // Highest contourvalue available
    float mincontourvalue;                // Lowest contourvalue available
    float opacity;                        // Opacity: 1 (opaque) 0 (transparent)
    float scalemin;                       // Lowest value for isocontour scale
    float scalemax;                       // Highest value for isocontour scale

    int logdlt;

    LineEdit *TXTcontourvalue;

    QCheckBox *CHKmpi;
    QCheckBox *CHKnormalgrad;
    QCheckBox *CHKshowgrid;
    QCheckBox *CHKtranslucence;

    QColor surfcolor;

    QDoubleSpinBox *SPBopacity;           // Surface opacity/transparency

    QDoubleValidator *myDoubleValidator;

    QGroupBox *FRMhighquality;
    QGroupBox *FRMsurfcolor;
    QGroupBox *FRMsurftype;

    togglingGroupBox *FRMisosurface;

    QLabel *LBLalpha;
    QLabel *LBLcontourvalue;
    QLabel *LBLfilename;
    QLabel *LBLmpi;
    QLabel *LBLopacity;
    QLabel *LBLscale;
    QLabel *LBLsensitive;
    QLabel *LBLstatus;

    QLineEdit *TXTisosurffilename;

    QList<QMetaObject::Connection> connections;

    QPoint initialposition;

    QProcess *myProcess;

    QPushButton *BTNexec;
    QPushButton *BTNstop;

    QRadioButton *RBTscalelin;
    QRadioButton *RBTscalelog;
    QRadioButton *RBTsolidsurf;           // Solid surface
    QRadioButton *RBTwiresurf;

    QSlider *SLDcontourvalue;
    QSlider *SLDopacity;

    QSpinBox *SPBmpi;
    QSpinBox *SPBsensitive;

    QString basename;
    QString fullname;                     // Full name for surface including path
    QString name;                         // Name for surface
    QString processname;
    QString ProjectFolder;
    QString ProjectName;
};

#endif // ISOSURFACE_H
