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
//  File:   plotter.h
//
//      Last version: March 2016
//
#ifndef PLOTTER_H
#define PLOTTER_H

#include <QCheckBox>
#include <QColorDialog>
#include <QDialog>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLineEdit>
#include <QLabel>
#include <QMap>
#include <QMessageBox>
#include <QMouseEvent>
#include <QObject>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <QPixmap>
#include <QRadioButton>
#include <QSize>
#include <QSpinBox>
#include <QStyle>
#include <QStylePainter>
#include <QStyleOptionFocusRect>
#include <QVBoxLayout>
#include <QVector>
#include <QVector2D>
#include <QVector3D>
#include <QWidget>
#include <QToolButton>

#include "ColorButton.h"
#include "elements.h"
    
static const int penstyle[4] = {
    Qt::SolidLine, Qt::DashLine, Qt::DotLine, Qt::DashDotLine
};
static const QColor colorForIds[13] = {
    Qt::red, Qt::green, Qt::blue, Qt::cyan, Qt::magenta, Qt::yellow, Qt::black,
    Qt::darkRed, Qt::darkGreen, Qt::darkBlue, Qt::darkCyan, Qt::darkMagenta, 
    Qt::darkYellow
};

class QColor;
class QDialog;
class QDockWidget;
class QFont;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPoint;
class QRect;
class QSpinBox;
class QToolButton;
class QToolButton;
class QWindow;
class PlotSettings;
class ContourLabel;

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

class Plotter : public QWidget
{
    Q_OBJECT

public:
    Plotter(QWidget *parent = 0);

    bool   getaspectratio();
    bool   getautomaticticks();
    bool   getcontourplot();
    bool   getfieldplot();
    bool   getfradplot();
    bool   getmulticolor();
    bool   getsghistplot();
    bool   getshowcoords();
    bool   getshowgrid();
    bool   getshowcontourlabels();
    bool   getshowcurvelabels();
    bool   getshowsghistlabels();
    bool   getshowXscalebottom();
    bool   getshowXscaletop();
    bool   getshowYscaleleft();
    bool   getshowYscaleright();
    bool   getshowXticksbottom();
    bool   getshowXtickstop();
    bool   getshowYticksleft();
    bool   getshowYticksright();
    bool   getshowtitle();
    bool   getshowtitlebox();
    bool   getshowXlabel();
    bool   getshowXlabelbox();
    bool   getshowYlabel();
    bool   getshowYlabelbox();
    bool   getsolidgrid();
    bool   gettranspbackground();
    bool   gettranspbckgrcapture();
    bool   getYlabelvert();
    double getcontourtolerance();
    double getumin();
    double getumax();
    double getvmin();
    double getvmax();
    double getxmin();
    double getxmax();
    double getymin();
    double getymax();

    int    getbasinspenwidth();
    int    getBottomMargin();
    int    getcntpenwidth();
    int    getcurvespenwidth();
    int    getimagequality();
    int    getLeftMargin();
    int    getnegativecontourstyle();
    int    getnlevels();
    int    getnu();
    int    getnv();
    int    getpositivecontourstyle();
    int    getRightMargin();
    int    getTopMargin();
    int    getXcifras();
    int    getYcifras();
    int    getzerocontourstyle();
    int    getzerotolerancepow();

    QColor getbasinscolor();
    QColor getbackgroundcolor();
    QColor getcenterscolor();
    QColor getcenterslabelcolor();
    QColor getcontourslabelcolor();
    QColor getfonttitlecolor();
    QColor getfontXlabelcolor();
    QColor getfontYlabelcolor();
    QColor getgridcolor();
    QColor getnegativecontourcolor();
    QColor getpositivecontourcolor();
    QColor getzerocontourcolor();
    QFont  getfontcenterslabel();
    QFont  getfontcurvelabels();
    QFont  getfontcontour();
    QFont  getfontsghistlabels();
    QFont  getfonttitle();
    QFont  getfontXlabel();
    QFont  getfontYlabel();
    
    void setaspectratio(bool aspectratio);
    void setautomaticticks(bool);
    void setcontourplot(bool);
    void setContourData(int id, const QVector<QPointF> &data);
    void setCurveData(int id, const QVector<QPointF> &data);
    void setdrawarrows(bool);
    void setEfieldData(int id, const QVector<QPointF> &data);
    void setexistcntatlbls(bool);
    void setfieldplot(bool);
    void setfradplot(bool);
    void setinvertarrows(bool);
    void setmulticolor(bool);
    void setplaneABC(QVector3D);
    void setplotbasins(bool);
    void setplotcenters(bool);
    void setplotcps(bool);
    void setshowcenterslabel(bool);
    void setplotbonds(bool);
    void setPlotSettings(const PlotSettings &settings);
    void setsghistplot(bool);
    void setshowcoords(bool);
    void setshowgrid(bool);
    void setshowcontourlabels(bool);
    void setshowcurvelabels(bool);
    void setshowsghistlabels(bool);
    void setshowXscalebottom(bool);
    void setshowXscaletop(bool);
    void setshowYscaleleft(bool);
    void setshowYscaleright(bool);
    void setshowXticksbottom(bool);
    void setshowXtickstop(bool);
    void setshowYticksleft(bool);
    void setshowYticksright(bool);
    void setshowtitle(bool);
    void setshowtitlebox(bool);
    void setshowXlabel(bool);
    void setshowXlabelbox(bool);
    void setshowYlabel(bool);
    void setshowYlabelbox(bool);
    void setsolidgrid(bool);
    void settranspbackground(bool);
    void settranspbckgrcapture(bool);
    void setuserdefresol(bool);
    void setYlabelvert(bool);

    void setbondsthr(double);
    void setcontourtolerance(double);
    void setumin(double);
    void setumax(double);
    void setvmin(double);
    void setvmax(double);
    void setumincnt(double);
    void setumaxcnt(double);
    void setvmincnt(double);
    void setvmaxcnt(double);
    void setuminelines(double);
    void setumaxelines(double);
    void setvminelines(double);
    void setvmaxelines(double);
    void setxmin(double);
    void setxmax(double);
    void setymin(double);
    void setymax(double);
    void setzero(double);

    void setarrowsseparation(int);
    void setarrowssize(int);
    void setarrowsskew(int);
    void setarrowswidth(int);
    void setbondswidth(int);
    void setBottomMargin(int);
    void setcenterslabelshft(int);
    void setcenterradius(int);
    void setbasinspenwidth(int);
    void setcntpenwidth(int);
    void setcurvespenwidth(int);
    void setcpsradius(int);
    void setcpstype(QVector<int>);
    void setimagequality(int);
    void setLeftMargin(int);
    void setnegativecontourstyle(int);
    void setnlevels(int);
    void setnu(int);
    void setnv(int);
    void setpositivecontourstyle(int);
    void setRightMargin(int);
    void setsghistpenwidth(int);
    void setshowcps(int, bool);
    void setTopMargin(int);
    void setwhoisabove(int);
    void setXcifras(int);
    void setYcifras(int); 
    void setzerocontourstyle(int);
    void setzerotolerancepow(int);

    void setbackgroundcolor(QColor);
    void setbasinscolor(QColor);
    void setbondscolor(QColor);
    void setcenterscolor(QColor);
    void setcenterslabel();
    void setcenterslabelcolor(QColor);
    void setcontourslabelcolor(QColor);
    void setcpscolor(QVector<QColor>);
    void setfonttitlecolor(QColor);
    void setfontXlabelcolor(QColor);
    void setfontYlabelcolor(QColor);
    void setgridcolor(QColor);
    void setnegativecontourcolor(QColor);
    void setpositivecontourcolor(QColor);
    void setzerocontourcolor(QColor);

    void setfontcenterslabel(QFont) ;
    void setfontcurvelabels(QFont) ;
    void setfontcontour(QFont) ;
    void setfontsghistlabels(QFont) ;
    void setfonttitle(QFont) ;
    void setfontXlabel(QFont) ;
    void setfontYlabel(QFont) ;

    void setimagesize(QSize);

    void setbasinsuv( QVector<QVector <QVector2D> >) ;
    void setcpsuvval(QVector <QVector3D>);

    void setatomlabels();
    void setbonds();
    void setcenterlabelrect();
    void updateallcenterslabelsrect();
    void updatecenterlabelrect(int);

    void capture_image();
    void clearCurve(int id, QMap<int, QVector<QPointF> > map);

    const char *tipo;

    int elinespenwidth;
    QLineEdit *texto;
    
    QRect titlerect;
    QRect Xlabelrect;
    QRect Ylabelrect;
    QString coordslabel;
    QString plotfile;
//	Qpalette variables
    QPalette *Pal;
//	QSize variables
    QSize *fmtitlesize;
//	QSring variables
    QString title;
    QString Xlabel;
    QString Ylabel;
//	QMap variables
    QMap<int, QVector<QPointF> > curveMap;
    QMap<int, QVector<QPointF> > efieldMap;
    QMap<int, QVector<QPointF> > contourMap;
//	QVector variables
    QVector<int> cpstype;
    QVector<int> *curvecontnum;
//    QVector<int> *sghistcontnum;
    QVector<bool> *curvehide;
    QVector<bool> *sghisthide;
    QVector<bool> *sghisthaslabel;
    QVector<bool> *plotarrows;
    QVector<ContourLabel> contourlabels;
    QVector<double> *contourlevels;
    QVector<double> *surf2d;
    QVector<QColor> *curvecolors;
    QVector<QColor> *efieldcolors;
    QVector<QColor> *sghistcolors;
    QVector<int> *efieldpenwidth;
    QVector<int> *sghistpenstyle;
    QVector<QRect> *curvelabelrect;
    QVector<QRect> *sghistlabelrect;
    QVector<QString> *centerslabel;
    QVector<QString> *curvelabels;
    QVector<QString> *sghistlabels;
    QVector<QVector <QVector2D> > basinsuv;
    QVector <QVector3D> cpsuvval;
    QVector<QVector3D> *rcen;
    QVector<QVector3D> *uvznuc;  // (u,v,znuc)  with znuc: nuclear charge (double)
    QVector< QVector2D> *bonds;
    QVector3D planeABC;
       
public slots:
    void zoomIn();
    void zoomOut();
    void refreshPixmap();
    void Xticks_changed(int);
    void Yticks_changed(int);
    void restore_automatic();
    void aceptar_label();
    void aceptarcontour_label();
    void aceptar_elinescolor();
    void aceptar_sghistlabel();
    void cerrar_label();
    void cerrar_elinescolor();
    void BTNcolor_clicked();
    void BTNelinescolor_clicked();
    void CHKhide_changed();
    void CHKsghisthide_changed();
    void updateplot();
    void setzoomRegion();
    void updatecontourlabels();

private slots:
    void emit_deletefrad();
    void emit_deletesghist();

protected:
    void paintEvent(QPaintEvent *event);
    void resizeEvent(QResizeEvent *event);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);
	
signals:
    void deletefrad(int);
    void deletesghist(int);
    void moveToTop();
    void title_dialog();
    void titleposChanged();
    void Xlabel_dialog();
    void XlabelposChanged();
    void Ylabel_dialog();
    void YlabelposChanged();
    int  Left_Margin_changed(int);
    void zoomChanged();
    void MarginsChanged(int bot, int top, int left, int right);

private:
    static const int max_cps=4;

    int  elines(int, int);
    void contourlabels_dialog(int);
    void curvelabels_dialog(int);
    void drawBasins(QPainter *painter, double scalewidth);
    void drawCurves(QPainter *painter, double scalewidth, QMap<int, QVector<QPointF> > map, int caseplot);
    void drawGrid(QPainter *painter, bool display, double scalefont);
    void elines_dialog(int);
    void sghistlabels_dialog(int);
    void Title_changed(QWidget *parent = 0);
    void updateRubberBandRegion();
    void zoom_changed();
    double value(double x, double y);

    ColorButton *BTNcolor;
    ColorButton *BTNelinescolor;

    QCheckBox *CHKhide;
    QCheckBox *CHKplotarrows;
    QDockWidget *gdock;
    QGroupBox *FRMelinescolor;
    QLabel *LBLefpenwidth;
    QLineEdit *TXTtitle;
    QPixmap pixmap;
    QRadioButton *RBTallelines;
    QRadioButton *RBTatomelines;
    QRadioButton *RBTsingleeline;
    QRect rubberBandRect;
    QSpinBox *SPBelinespenwidth;
    QToolButton *zoomInButton;
    QToolButton *zoomOutButton;
    QVector<QColor> cpscolor;
    QVector<PlotSettings> zoomStack;
    QVector<QRect> *centerlabelrect;
    QVector<QVector2D> *centerlabelspos;

    bool aspectratio;
    bool automaticticks;
    bool drawarrows;
    bool contourplot;
    bool fieldplot;
    bool fradplot;
    bool drawrecttitle;
    bool drawrectXlabel;
    bool drawrectYlabel;
    bool existcntatlbls;
    bool invertarrows;
    bool movecenterlabel;
    bool movecurvelabel;
    bool movesghistlabel;
    bool movetitle;
    bool moveXlabel;
    bool moveYlabel;
    bool multicolor;
    bool painterinuse;
    bool plotbasins;
    bool plotbonds;
    bool plotcenters;
    bool plotcps;
    bool rubberBandIsShown;
    bool sghistplot;
    bool showcenterslabel;
    bool showcoords;
    bool showcps[max_cps];
    bool showgrid;
    bool showcontourlabels;
    bool showcurvelabels;
    bool showsghistlabels;
    bool showXscalebottom;
    bool showXscaletop;
    bool showYscaleleft;
    bool showYscaleright;
    bool showXticksbottom;
    bool showXtickstop;
    bool showYticksleft;
    bool showYticksright;
    bool showtitle;
    bool showtitlebox;
    bool showXlabel;
    bool showXlabelbox;
    bool showYlabel;
    bool showYlabelbox;
    bool solidgrid;
    bool transpbackground;
    bool transpbckgrcapture;
    bool updateYpos;
    bool userdefresol;
    bool Ylabelvert;
    double bondsthr;
    double contourtolerance;
    double umin;
    double umax;
    double vmin;
    double vmax;
    double umincnt;
    double umaxcnt;
    double vmincnt;
    double vmaxcnt;
    double uminelines;
    double umaxelines;
    double vminelines;
    double vmaxelines;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zero;

    int arrowsseparation;
    int arrowssize;
    int arrowsskew;
    int arrowswidth;
    int bondswidth;
    int centerslabelshft;
    int centerradius;
    int cpsradius;
    int curZoom;
    int plotcase;
    int dotx;
    int doty;
    int imagequality;
    int indcenter;
    int indcurve;
    int indeline;
    int nlevels;
    int nu;
    int nv;
    int Xcifras;
    int Ycifras;
    int negativecontourstyle;
    int positivecontourstyle;
    int zerocontourstyle;
    int basinspenwidth;
    int cntpenwidth;
    int curvespenwidth;
    int sghistpenwidth;
    int TopMargin;
    int BottomMargin;
    int LeftMargin;
    int RightMargin;
    int whoisabove;
    int zerotolerancepow;
    QColor backgroundcolor;
    QColor basinscolor;
    QColor bondscolor;
    QColor color;
    QColor elinescolor;
    QColor centerscolor;
    QColor centerslabelcolor;
    QColor contourslabelcolor;
    QColor fonttitlecolor;
    QColor fontXlabelcolor;
    QColor fontYlabelcolor;
    QColor gridcolor;
    QColor negativecontourcolor;
    QColor positivecontourcolor;
    QColor zerocontourcolor;
    QDialog *dialog;
    QDialog *dialogelines;
    QFont fontcenterslabel;
    QFont fontcurvelabels;
    QFont fontcontour;
    QFont fontsghistlabels;
    QFont fonttitle;
    QFont fontXlabel;
    QFont fontYlabel;
    QRect coordsrect;
    QSize canvassize;
    QSize imagesize;
	
};

class PlotSettings
{
public:
    PlotSettings();

    void scroll(int dx, int dy);
    void set_numXticks(int numx);
    void set_numYticks(int numy);
    void adjust();
    void adjustforced();
    double spanX() const { return maxX - minX; }
    double spanY() const { return maxY - minY; }
    void adjustAxisforced(double &min, double &max, int &numTicks);
    

    double minX;
    double maxX;
    int ncifrasx;
    int ncifrasy;
    int numXTicks;
    double minY;
    double maxY;
    int numYTicks;
    int TopMargin;
    int BottomMargin;
    int LeftMargin;
    int RightMargin;

private:
    static void adjustAxis(double &min, double &max, int &numTicks);
};

class ContourLabel
{
public:
    ContourLabel();
    ~ContourLabel();

    double getu();
    double getv();
    int getx();
    int gety();
    QRect getrect();
    QString getlabel();

    void setlabel(QString);
    void setrect(QRect);
    void setu(double);
    void setv(double);
    void setx(int);
    void sety(int);

private:
    QRect rect; 
    QString label;
    int x;
    int y;
    double u;
    double v;
};

#endif

