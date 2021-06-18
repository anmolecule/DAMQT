#include <QtCore/qmath.h>
#include <QDialog>
#include <QCheckBox>
#include <QColorDialog>
#include <QFontDialog>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>

#include "surface.h"
#include "CIsoSurface.h"
#include <math.h>

surface::surface(QWidget *parent) : QWidget(parent)
{
    allindices.clear();
    allvalues.clear();
    allvertices.clear();
    connections.clear();
    extremactive[0].clear();
    extremactive[1].clear();
    griddimensions.clear();
    gridindices.clear();
    gridindicesoffset.clear();
    gridvertices.clear();
    gridnxyz.clear();
    localextrema[0].clear();
    localextrema[1].clear();
    surfacecolor = QColor(255, 255, 0);
    font = QFont("Noto Sans", 14, QFont::Bold);
    fontcolor = QColor(255, 255, 0, 255);
    initialposition = QPoint(200,300);
    myDoubleValidator = nullpointer;;
    onlyextremactive = false;
    showcolorrule = false;
    showgridbounds = false;
    showextremacoords = false;
    showextremaindices = true;
    showextremasymbols = true;
    showextremavalues = false;
    showlocalmax = false;
    showlocalmin = false;
    solidsurf = true;

    ballradius = 10;
    coordprecision = 2;
    fabstop = 1.0;
    name = "";
    nlocalmax = 0;
    nlocalmin = 0;
    opacity = 1.0;
    topcolor = 0.9*fabstop;
    valueprecision = 2;
    vshift = 0;

    extremacolor[0] = QColor(245,114,15);
    extremacolor[1] = QColor(162,12,255);

    makeSphere(15,15);    // computes vertices of a sphere and its indices
    settranslucence(false);
    setvisible(true);

}

surface::~surface(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    allindicesextrema->clear();
    allindicesoffsetextrema->clear();
    allverticesextrema->clear();
    allindices.clear();
    allvalues.clear();
    allvertices.clear();
    gridindices.clear();
    gridindicesoffset.clear();
    gridvertices.clear();
}

//  ------------------------------------------------------------------------------------------------------------------
//
//          Functions and slots for surface editor
//
//  ------------------------------------------------------------------------------------------------------------------

void surface::surfacecolor_changed(QColor color)
{
    surfacecolor = color;
    VertexNormalData v;
    for (int i = 0 ; i < allvertices.length() ; i++){
        v = allvertices.at(i);
        v.color.setX(surfacecolor.redF());
        v.color.setY(surfacecolor.greenF());
        v.color.setZ(surfacecolor.blueF());
        allvertices.replace(i,v);
    }
}

//    Button for choosing fonts for atom labels
void surface::changeextremafont()
{
    font = QFontDialog::getFont(nullpointer, font);
    emit updatedisplay();
}

void surface::resetsurface(){
    float tc2i = 1.f / (topcolor*topcolor);
    VertexNormalData v;
    for (int i = 0 ; i < allvalues.count() ; i++){
        v = allvertices.at(i);
        if (allvalues.at(i) >= topcolor){
            v.color.setX(1.);
            v.color.setY(0.);
            v.color.setZ(0.);
        }
        else if(allvalues.at(i) <= -topcolor)
        {
            v.color.setX(0.);
            v.color.setY(0.);
            v.color.setZ(1.);
        }
        else if(allvalues.at(i) >= 0.f){
            v.color.setX(allvalues.at(i)/topcolor);
            v.color.setY(sqrtf(1.f-tc2i*allvalues.at(i)*allvalues.at(i)));
            v.color.setZ(0.);
        }
        else {
            v.color.setX(0.);
            v.color.setY(sqrtf(1.f-tc2i*allvalues.at(i)*allvalues.at(i)));
            v.color.setZ(-allvalues.at(i)/topcolor);
        }
        v.color.setW(opacity);
        allvertices.replace(i,v);
    }
}

//  ------------------------------------------------------------------------------------------------------------------
//
//      General purpose functions
//
//  ------------------------------------------------------------------------------------------------------------------


int surface::getballradius(){
    return ballradius;
}

// Changes ball radius
void surface::setballradius(int i){
    ballradius = i;
    createballs();
    emit updatedisplay();
}

void surface::extremaselectall(){
//    qDebug() << "entra en extremaselectall";
    for (int i=0 ; i<2 ; i++){
        for (int j = 0 ; j < extremactive[i].length() ; j++){
            setextremactive(i, j, true);
        }
    }
    emit updatedisplay();
}

void surface::extremaselectnone(){
//    qDebug() << "entra en extremaselectnone";
    for (int i=0 ; i<2 ; i++){
        for (int j = 0 ; j < extremactive[i].length() ; j++){
            setextremactive(i, j, false);
        }
    }
    emit updatedisplay();
}

void surface::generategridbounds(float *a){
    float xmin, xmax, ymin, ymax, zmin, zmax;
    griddimensions.clear();
    gridindices.clear();
    gridindicesoffset.clear();
    gridvertices.clear();
    xmin = a[0];
    xmax = a[1];
    ymin = a[2];
    ymax = a[3];
    zmin = a[4];
    zmax = a[5];
    QVector4D color = QVector4D(1.,0.5,0.,1.);
    VertexNormalData v;
    v.normal.setX(0.);
    v.normal.setY(0.);
    v.normal.setZ(1.);
    v.color.setX(color.x());
    v.color.setY(color.y());
    v.color.setZ(color.z());
    v.color.setW(color.w());
    v.position.setX(xmin);
    v.position.setY(ymin);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymin);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymax);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmin);
    v.position.setY(ymax);
    v.position.setZ(zmin);
    gridvertices.append(v);
    v.position.setX(xmin);
    v.position.setY(ymin);
    v.position.setZ(zmax);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymin);
    v.position.setZ(zmax);
    gridvertices.append(v);
    v.position.setX(xmax);
    v.position.setY(ymax);
    v.position.setZ(zmax);
    gridvertices.append(v);
    v.position.setX(xmin);
    v.position.setY(ymax);
    v.position.setZ(zmax);
    gridvertices.append(v);
    gridindicesoffset << 0 << 2 << 4 << 6 << 8 << 10 << 12 << 14 << 16 << 18 << 20 << 22 << 24;
    gridindices << 0 << 1 << 1 << 2 << 2 << 3 << 3 << 0 << 4 << 5 << 5 << 6 << 6 << 7 << 7 << 4
                << 0 << 4 << 1 << 5 << 2 << 6 << 3 << 7;
}

QPoint surface::getinitialposition(){
    return initialposition;
}

bool surface::getextremactive(int i, int j){
    return extremactive[i].at(j);
}

bool surface::getshowcolorrule(){
    return showcolorrule;
}

bool surface::getonlyextremactive(){
    return onlyextremactive;
}

bool surface::getshowextremacoords(){
    return showextremacoords;
}

bool surface::getshowgridbounds(){
    return showgridbounds;
}

bool surface::getshowlocalmax(){
    return showlocalmax;
}

bool surface::getshowlocalmin(){
    return showlocalmin;
}

bool surface::getshowextremaindices(){
    return showextremaindices;
}

bool surface::getshowextremasymbols(){
    return showextremasymbols;
}

bool surface::getshowextremavalues(){
    return showextremavalues;
}

bool surface::getsolidsurf(){
    return solidsurf;
}

bool surface::gettranslucence(){
    return translucence;
}

bool surface::getvisible(){
    return visible;
}

float surface::getopacity(){
    return opacity;
}

float surface::gettopcolor(){
    return topcolor;
}

int surface::getcoordprecision(){
    return coordprecision;
}

int surface::getvalueprecision(){
    return valueprecision;
}

int surface::getvshift(){
    return vshift;
}

QColor surface::getsurfacecolor(){
    return surfacecolor;
}

QColor surface::getfontcolor(){
    return fontcolor;
}

QFont surface::getfont(){
    return font;
}

QString surface::getfullname(){
    return fullname;
}

QString surface::getname(){
    return name;
}

QVector <GLuint> surface::getallindices(){
    return allindices;
}

QVector <VertexNormalData> surface::getallvertices(){
    return allvertices;
}


void surface::opacity_changed(int ix){
    setopacity(float(ix)/100.);
}

void surface::opacity_released(){
    VertexNormalData v;
    for (int i = 0 ; i < allvertices.count() ; i++){
        v.position = allvertices.at(i).position;
        v.normal = allvertices.at(i).normal;
        v.color = allvertices.at(i).color;
        v.color.setW(opacity);
        allvertices.replace(i,v);
    }
    emit updatedisplay();
}

bool surface::readbasins(QString filename){
    allindices.clear();
    allvalues.clear();
    allvertices.clear();
    int  iaux[5];
    float vaux[7];
    FILE *f;
    QByteArray ba = filename.toLatin1();

    const char *file = ba.data();

//    Open files with function and existderivs
    f   = fopen(file   , "rb" );

    fread( iaux , sizeof(int) , 1 , f);
    nvertices = iaux[0];
    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *bar;
    bar = new QProgressBar(win);
    bar->resize(300,25);
    bar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    bar->setMinimumWidth(300);
    bar->setMaximumWidth(300);
    bar->setMinimum(0);
    bar->setMaximum(nvertices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(bar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    for(int i = 0 ; i < nvertices ; i++){
        size_t nread = fread(vaux , sizeof(float) , 6 , f);
        if (nread != 6){
            QMessageBox msgBox;
            msgBox.setText(tr("readbasins"));
            msgBox.setInformativeText(tr("Error reading vertices from file %1").arg(filename));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return false;
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(surfacecolor.redF());
        v.color.setY(surfacecolor.greenF());
        v.color.setZ(surfacecolor.blueF());
        v.color.setW(opacity);
        allvertices.append(v);
        allvalues.append(topcolor);
    }
    if (allvertices.count() != nvertices){
        QMessageBox msgBox;
        msgBox.setText(tr("readbasins"));
        msgBox.setInformativeText(tr("Wrong number of vertices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    QVector3D r[3];
    for(int i = 0 ; i < nvertices ; i += 3){  // Loads indices
        r[0] = allvertices.at(i).position;
        r[1] = allvertices.at(i+1).position;
        r[2] = allvertices.at(i+2).position;
        allindices.append(i);   // Indices coming from file .basins start in 1 (Fortran convention)
        if (QVector3D::dotProduct(r[0]+r[1]+r[2],QVector3D::crossProduct(r[0]-r[1],r[2]-r[1])) >= 0){
            allindices.append(i+1);
            allindices.append(i+2);
        }
        else{
            allindices.append(i+2);
            allindices.append(i+1);
        }
    }
    showcolorrule = false;
    return true;
}

bool surface::readbasinsnew(QString filename){
    allindices.clear();
    allvalues.clear();
    allvertices.clear();
    int  iaux;
    float vaux[7];
//    Open files with function
    QFile inputfile(filename);
    if (!inputfile.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename;
        return false;
    }
    QDataStream data (&inputfile);

    QByteArray bar;
    bar = inputfile.read(sizeof(int));
    memcpy(&iaux, bar.constData(), sizeof(int));
    nvertices = iaux;
    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *prgbar;
    prgbar = new QProgressBar(win);
    prgbar->resize(300,25);
    prgbar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    prgbar->setMinimumWidth(300);
    prgbar->setMaximumWidth(300);
    prgbar->setMinimum(0);
    prgbar->setMaximum(nvertices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(prgbar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    for(int i = 0 ; i < nvertices ; i++){
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(float));
            memcpy(&vaux[i], bar.constData(), sizeof(float));
//            qDebug() << "vaux[" << i << "] = " << vaux[i];
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(surfacecolor.redF());
        v.color.setY(surfacecolor.greenF());
        v.color.setZ(surfacecolor.blueF());
        v.color.setW(opacity);
        allvertices.append(v);
        allvalues.append(topcolor);
    }
    if (allvertices.count() != nvertices){
        QMessageBox msgBox;
        msgBox.setText(tr("readbasins"));
        msgBox.setInformativeText(tr("Wrong number of vertices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    QVector3D r[3];
    for(int i = 0 ; i < nvertices ; i += 3){  // Loads indices
        r[0] = allvertices.at(i).position;
        r[1] = allvertices.at(i+1).position;
        r[2] = allvertices.at(i+2).position;
        allindices.append(i);   // Indices coming from file .basins start in 1 (Fortran convention)
        if (QVector3D::dotProduct(r[0]+r[1]+r[2],QVector3D::crossProduct(r[0]-r[1],r[2]-r[1])) >= 0){
            allindices.append(i+1);
            allindices.append(i+2);
        }
        else{
            allindices.append(i+2);
            allindices.append(i+1);
        }
    }
    showcolorrule = false;
    return true;
}

bool surface::readsgh(QString filename){
    allindices.clear();
    allvalues.clear();
    allvertices.clear();
    int  iaux[5];
    bool fdouble;
    float vaux[7];
    double dvaux[7];
    FILE *f;
    QByteArray ba = filename.toLatin1();

    const char *file = ba.data();

    f   = fopen(file   , "rb" );

//    Read file with function
    fread( iaux , sizeof(int) , 5 , f);  //iaux[0] (== 0 double, != 0 float); iaux[2] = not used; iaux[3-5] = nx, ny, nz
    if (iaux[0] == 0)
        fdouble = true;
    else
        fdouble = false;
    for (int i = 0 ; i < 3 ; i++){
        gridnxyz.append(iaux[i+2]);
    }
    if (fdouble){       // (d)vaux[0-5] = xini, xfin, yini, yfin, zini, zfin
        fread(dvaux , sizeof(double) , 6 , f);
        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append((float)dvaux[i]);
            vaux[i] = dvaux[i];
        }
    }
    else{
        fread(vaux , sizeof(float) , 6 , f);
        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append(vaux[i]);
        }
    }
    generategridbounds(vaux);
    fread( iaux , sizeof(int) , 2 , f);  //iaux[0] = number of vertices, iaux[1] = number of indices
    nvertices = iaux[0];
    nindices = iaux[1];
    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *bar;
    bar = new QProgressBar(win);
    bar->resize(300,25);
    bar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    bar->setMinimumWidth(300);
    bar->setMaximumWidth(300);
    bar->setMinimum(0);
    bar->setMaximum(nvertices+nindices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(bar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    fmax = -1.e20;
    fmin = 1.e20;
    for(int i = 0 ; i < nvertices ; i++){
        if (fdouble){
            size_t nread = fread(dvaux , sizeof(double) , 7 , f);
            if (nread != 7){
                QMessageBox msgBox;
                msgBox.setText(tr("readsgh"));
                msgBox.setInformativeText(tr("Error reading vertices from file %1").arg(filename));
                msgBox.setIcon(QMessageBox::Warning);
                msgBox.exec();
                return false;
            }
            for(int i = 0 ; i < 7 ; i++){
                vaux[i] = (float)dvaux[i];
            }
        }
        else{
            size_t nread = fread(vaux , sizeof(float) , 7 , f);
            if (nread != 7){
                QMessageBox msgBox;
                msgBox.setText(tr("readsgh"));
                msgBox.setInformativeText(tr("Error reading vertices from file %1").arg(filename));
                msgBox.setIcon(QMessageBox::Warning);
                msgBox.exec();
                return false;
            }
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(1.);
        v.color.setY(0.);
        v.color.setZ(0.);
        v.color.setW(opacity);
        allvertices.append(v);
        allvalues.append(vaux[6]);
        fmax = qMax(fmax,vaux[6]);
        fmin = qMin(fmin,vaux[6]);
    }
    topcolor = 0.9f*qMax(fmax,qAbs(fmin));
    settrianglecolors();
    for(int i = 0 ; i < nindices/3 ; i++){  // Reads indices in groups of three corresponding to vertices of a triangle
        size_t nread = fread(iaux , sizeof(int) , 3 , f);
        if (nread != 3){
            QMessageBox msgBox;
            msgBox.setText(tr("readsgh"));
            msgBox.setInformativeText(tr("Error reading allindices from file %1").arg(filename));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return false;
        }
        allindices.append(iaux[0]-1);   // Indices coming from file .srf start in 1 (Fortran convention)
        allindices.append(iaux[1]-1);
        allindices.append(iaux[2]-1);
    }
    if (allvertices.count() != nvertices || allindices.count() != nindices){
        QMessageBox msgBox;
        msgBox.setText(tr("readsgh"));
        msgBox.setInformativeText(tr("Wrong number of vertices and allindices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    showcolorrule = true;
    fabstop = 1.1*qMax(fmax,qAbs(fmin));
    topcolor = 0.9*fabstop;
    localextrema[0].clear();
    size_t nread = fread(iaux , sizeof(int) , 1 , f);
    if (!feof(f) && nread == 1 && iaux[0] > 0){
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                size_t nread = fread(dvaux , sizeof(double) , 4 , f);
                if (nread != 9){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsgh"));
                    msgBox.setInformativeText(tr("Error reading local maxima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
                for(int i = 0 ; i < 4 ; i++){
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                size_t nread = fread(vaux , sizeof(float) , 4 , f);
                if (nread != 4){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsgh"));
                    msgBox.setInformativeText(tr("Error reading local maxima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
            }
            localextrema[0].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[0].append(false);
        }
    }
    localextrema[1].clear();
    nread = fread(iaux , sizeof(int) , 1 , f);
    if (!feof(f) && nread == 1 && iaux[0] > 0){
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                size_t nread = fread(dvaux , sizeof(double) , 4 , f);
                if (nread != 9){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsgh"));
                    msgBox.setInformativeText(tr("Error reading local minima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
                for(int i = 0 ; i < 4 ; i++){
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                size_t nread = fread(vaux , sizeof(float) , 4 , f);
                if (nread != 4){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsgh"));
                    msgBox.setInformativeText(tr("Error reading local minima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
            }
            localextrema[1].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[1].append(false);
        }
    }
    createballs();
    return true;
}


bool surface::readsrf(QString filename){
    allindices.clear();
    allvertices.clear();
    int  iaux[5];
    bool fdouble;
    float vaux[9];
    double dvaux[9];
    FILE *f;
    QByteArray ba = filename.toLatin1();

    const char *file = ba.data();

//    Open files with function
    f   = fopen(file   , "rb" );

//    Read file with function
    fread( iaux , sizeof(int) , 5 , f);  //iaux[0] (== 0 double, != 0 float); iaux[2] = not used; iaux[3-5] = nx, ny, nz
    if (iaux[0] == 0)
        fdouble = true;
    else
        fdouble = false;
    for (int i = 0 ; i < 3 ; i++){
        gridnxyz.append(iaux[i+2]);
    }
    if (fdouble){       // (d)vaux[0-5] = xini, xfin, yini, yfin, zini, zfin
        fread(dvaux , sizeof(double) , 6 , f);
        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append((float)dvaux[i]);
        }
    }
    else{
        fread(vaux , sizeof(float) , 6 , f);
        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append(vaux[i]);
        }
    }
    generategridbounds(vaux);
    fread( iaux , sizeof(int) , 2 , f);  //iaux[0] = number of vertices, iaux[2] = number of indices
    nvertices = iaux[0];
    nindices = iaux[1];
    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *bar;
    bar = new QProgressBar(win);
    bar->resize(300,25);
    bar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    bar->setMinimumWidth(300);
    bar->setMaximumWidth(300);
    bar->setMinimum(0);
    bar->setMaximum(nvertices+nindices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(bar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    for(int i = 0 ; i < nvertices ; i++){
        if (fdouble){
            size_t nread = fread(dvaux , sizeof(double) , 9 , f);
            if (nread != 9){
                QMessageBox msgBox;
                msgBox.setText(tr("readsrf"));
                msgBox.setInformativeText(tr("Error reading vertices from file %1").arg(filename));
                msgBox.setIcon(QMessageBox::Warning);
                msgBox.exec();
                return false;
            }
            for(int i = 0 ; i < 9 ; i++){
                vaux[i] = (float)dvaux[i];
            }
        }
        else{
            size_t nread = fread(vaux , sizeof(float) , 9 , f);
            if (nread != 9){
                QMessageBox msgBox;
                msgBox.setText(tr("readsrf"));
                msgBox.setInformativeText(tr("Error reading vertices from file %1").arg(filename));
                msgBox.setIcon(QMessageBox::Warning);
                msgBox.exec();
                return false;
            }
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(vaux[6]);
        v.color.setY(vaux[7]);
        v.color.setZ(vaux[8]);
        v.color.setW(opacity);
        allvertices.append(v);
    }
    for(int i = 0 ; i < nindices/3 ; i++){  // Reads indices in groups of three corresponding to vertices of a triangle
        size_t nread = fread(iaux , sizeof(int) , 3 , f);
        if (nread != 3){
            QMessageBox msgBox;
            msgBox.setText(tr("readsrf"));
            msgBox.setInformativeText(tr("Error reading allindices from file %1").arg(filename));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            return false;
        }
        allindices.append(iaux[0]-1);   // Indices coming from file .srf start in 1 (Fortran convention)
        allindices.append(iaux[1]-1);
        allindices.append(iaux[2]-1);
    }
    if (allvertices.count() != nvertices || allindices.count() != nindices){
        QMessageBox msgBox;
        msgBox.setText(tr("readsrf"));
        msgBox.setInformativeText(tr("Wrong number of vertices and allindices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    showcolorrule = false;
    localextrema[0].clear();
    size_t nread = fread(iaux , sizeof(int) , 1 , f);
    if (!feof(f) && nread == 1 && iaux[0] > 0){
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                size_t nread = fread(dvaux , sizeof(double) , 4 , f);
                if (nread != 9){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsrf"));
                    msgBox.setInformativeText(tr("Error reading local maxima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
                for(int i = 0 ; i < 4 ; i++){
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                size_t nread = fread(vaux , sizeof(float) , 4 , f);
                if (nread != 4){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsrf"));
                    msgBox.setInformativeText(tr("Error reading local maxima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
            }
            localextrema[0].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[0].append(false);
        }
    }
    localextrema[1].clear();
    nread = fread(iaux , sizeof(int) , 1 , f);
    if (!feof(f) && nread == 1 && iaux[0] > 0){
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                size_t nread = fread(dvaux , sizeof(double) , 4 , f);
                if (nread != 9){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsrf"));
                    msgBox.setInformativeText(tr("Error reading local minima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
                for(int i = 0 ; i < 4 ; i++){
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                size_t nread = fread(vaux , sizeof(float) , 4 , f);
                if (nread != 4){
                    QMessageBox msgBox;
                    msgBox.setText(tr("readsrf"));
                    msgBox.setInformativeText(tr("Error reading local minima from file %1").arg(filename));
                    msgBox.setIcon(QMessageBox::Warning);
                    msgBox.exec();
                    break;
                }
            }
            localextrema[1].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[1].append(false);
        }
    }
    createballs();
    return true;
}

bool surface::readisosurfnew(QString filename){
    allindices.clear();
    allvertices.clear();
    int  iaux[5];
    bool fdouble;
    float vaux[7];
    double dvaux[7];
    QFile inputfile(filename);
    if (!inputfile.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename;
        return false;
    }
    QDataStream data (&inputfile);
//    qDebug() << "sizeof(int) = " << sizeof(int);
    ;

//      Read file with function
//      iaux[0] (== 0 double, != 0 float); iaux[2] = not used; iaux[3-5] = nx, ny, nz
    for (int i = 0 ; i < 5 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }

    if (iaux[0] == 0)
        fdouble = true;
    else
        fdouble = false;
    for (int i = 0 ; i < 3 ; i++){
        gridnxyz.append(iaux[i+2]);
    }

    if (fdouble){       // (d)vaux[0-5] = xini, xfin, yini, yfin, zini, zfin
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(double));
            memcpy(&dvaux[i], bar.constData(), sizeof(double));
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
        for (int i = 0 ; i < 6 ; i++){
            vaux[i] = dvaux[i];
            griddimensions.append(vaux[i]);
        }
    }
    else{
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(float));
            memcpy(&vaux[i], bar.constData(), sizeof(float));
//            qDebug() << "vaux[" << i << "] = " << vaux[i];
        }

        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append(vaux[i]);
        }
    }

    generategridbounds(vaux);
    for (int i = 0; i < 2 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }

    nvertices = iaux[0];
    nindices = iaux[1];

    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *prgbar;
    prgbar = new QProgressBar(win);
    prgbar->resize(300,25);
    prgbar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    prgbar->setMinimumWidth(300);
    prgbar->setMaximumWidth(300);
    prgbar->setMinimum(0);
    prgbar->setMaximum(nvertices+nindices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(prgbar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    fmax = -1.e20;
    fmin = 1.e20;
    for(int i = 0 ; i < nvertices ; i++){
        if (fdouble){
            for (int i = 0 ; i < 6 ; i++){
                QByteArray bar;
                bar = inputfile.read(sizeof(double));
                memcpy(&dvaux[i], bar.constData(), sizeof(double));
                vaux[i] = (float)dvaux[i];
            }
        }
        else{
            for (int i = 0 ; i < 6 ; i++){
                QByteArray bar;
                bar = inputfile.read(sizeof(float));
                memcpy(&vaux[i], bar.constData(), sizeof(float));
            }
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(surfacecolor.redF());
        v.color.setY(surfacecolor.greenF());
        v.color.setZ(surfacecolor.blueF());
        v.color.setW(opacity);
        allvertices.append(v);
        fmax = qMax(fmax,vaux[6]);
        fmin = qMin(fmin,vaux[6]);
    }
    for(int i = 0 ; i < nindices/3 ; i++){  // Reads indices in groups of three corresponding to vertices of a triangle
        for (int i = 0; i < 3 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(int));
            memcpy(&iaux[i], bar.constData(), sizeof(int));
        }
        allindices.append(iaux[0]-1);   // Indices coming from file .srf start in 1 (Fortran convention)
        allindices.append(iaux[1]-1);
        allindices.append(iaux[2]-1);
    }

    if (allvertices.count() != nvertices || allindices.count() != nindices){
        QMessageBox msgBox;
        msgBox.setText(tr("readsgh"));
        msgBox.setInformativeText(tr("Wrong number of vertices and allindices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    showcolorrule = false;

    return true;
}


bool surface::readsghnew(QString filename){
    allindices.clear();
    allvalues.clear();
    allvertices.clear();
    int  iaux[5];
    bool fdouble;
    float vaux[7];
    double dvaux[7];
    QFile inputfile(filename);
    if (!inputfile.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename;
        return false;
    }
    QDataStream data (&inputfile);
//    qDebug() << "sizeof(int) = " << sizeof(int);

//      Read file with function
//      iaux[0] (== 0 double, != 0 float); iaux[2] = not used; iaux[3-5] = nx, ny, nz
    for (int i = 0 ; i < 5 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }

    if (iaux[0] == 0)
        fdouble = true;
    else
        fdouble = false;
    for (int i = 0 ; i < 3 ; i++){
        gridnxyz.append(iaux[i+2]);
    }

    if (fdouble){       // (d)vaux[0-5] = xini, xfin, yini, yfin, zini, zfin
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(double));
            memcpy(&dvaux[i], bar.constData(), sizeof(double));
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
        for (int i = 0 ; i < 6 ; i++){
            vaux[i] = dvaux[i];
            griddimensions.append(vaux[i]);
        }
    }
    else{
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(float));
            memcpy(&vaux[i], bar.constData(), sizeof(float));
//            qDebug() << "vaux[" << i << "] = " << vaux[i];
        }

        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append(vaux[i]);
        }
    }

    generategridbounds(vaux);
    for (int i = 0; i < 2 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }

    nvertices = iaux[0];
    nindices = iaux[1];

    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *prgbar;
    prgbar = new QProgressBar(win);
    prgbar->resize(300,25);
    prgbar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    prgbar->setMinimumWidth(300);
    prgbar->setMaximumWidth(300);
    prgbar->setMinimum(0);
    prgbar->setMaximum(nvertices+nindices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(prgbar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    fmax = -1.e20;
    fmin = 1.e20;
    for(int i = 0 ; i < nvertices ; i++){
        if (fdouble){
            for (int i = 0 ; i < 7 ; i++){
                QByteArray bar;
                bar = inputfile.read(sizeof(double));
                memcpy(&dvaux[i], bar.constData(), sizeof(double));
                vaux[i] = (float)dvaux[i];
            }
        }
        else{
            for (int i = 0 ; i < 7 ; i++){
                QByteArray bar;
                bar = inputfile.read(sizeof(float));
                memcpy(&vaux[i], bar.constData(), sizeof(float));
            }
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(1.);
        v.color.setY(0.);
        v.color.setZ(0.);
        v.color.setW(opacity);
        allvertices.append(v);
        allvalues.append(vaux[6]);
        fmax = qMax(fmax,vaux[6]);
        fmin = qMin(fmin,vaux[6]);
    }
    topcolor = 0.9f*qMax(fmax,qAbs(fmin));
    settrianglecolors();
    for(int i = 0 ; i < nindices/3 ; i++){  // Reads indices in groups of three corresponding to vertices of a triangle
        for (int i = 0; i < 3 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(int));
            memcpy(&iaux[i], bar.constData(), sizeof(int));
        }
        allindices.append(iaux[0]-1);   // Indices coming from file .srf start in 1 (Fortran convention)
        allindices.append(iaux[1]-1);
        allindices.append(iaux[2]-1);
    }

    if (allvertices.count() != nvertices || allindices.count() != nindices){
        QMessageBox msgBox;
        msgBox.setText(tr("readsgh"));
        msgBox.setInformativeText(tr("Wrong number of vertices and allindices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    showcolorrule = true;
    fabstop = 1.1*qMax(fmax,qAbs(fmin));
    topcolor = 0.9*fabstop;
    localextrema[0].clear();

//    qDebug() << "allindices.count() = " << allindices.count();
//    qDebug() << "allvertices.count() = " << allvertices.count();

    QByteArray bar;
    bar = inputfile.read(sizeof(int));
    if (!inputfile.atEnd()){
        memcpy(&iaux[0], bar.constData(), sizeof(int));
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(double));
                    memcpy(&dvaux[i], bar.constData(), sizeof(double));
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(float));
                    memcpy(&vaux[i], bar.constData(), sizeof(float));
                }
            }
            localextrema[0].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[0].append(false);
        }
    }

//    qDebug() << "localextrema[0].count() = " << localextrema[0].count();

    localextrema[1].clear();

    bar = inputfile.read(sizeof(int));
    if (!inputfile.atEnd()){
        memcpy(&iaux[0], bar.constData(), sizeof(int));
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(double));
                    memcpy(&dvaux[i], bar.constData(), sizeof(double));
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(float));
                    memcpy(&vaux[i], bar.constData(), sizeof(float));
                }
            }
            localextrema[1].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[1].append(false);
        }
    }
//qDebug() << "localextrema[1].count() = " << localextrema[1].count();
    createballs();
    return true;
}

bool surface::readsrfnew(QString filename){
    allindices.clear();
    allvertices.clear();
    int  iaux[5];
    bool fdouble;
    float vaux[9];
    double dvaux[9];

    QFile inputfile(filename);
    if (!inputfile.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename;
        return false;
    }
    QDataStream data (&inputfile);

//      Read file with function
//      iaux[0] (== 0 double, != 0 float); iaux[2] = not used; iaux[3-5] = nx, ny, nz
    for (int i = 0 ; i < 5 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }
    if (iaux[0] == 0)
        fdouble = true;
    else
        fdouble = false;
    for (int i = 0 ; i < 3 ; i++){
        gridnxyz.append(iaux[i+2]);
    }
    if (fdouble){       // (d)vaux[0-5] = xini, xfin, yini, yfin, zini, zfin
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(double));
            memcpy(&dvaux[i], bar.constData(), sizeof(double));
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
        for (int i = 0 ; i < 6 ; i++){
            vaux[i] = dvaux[i];
            griddimensions.append((float)vaux[i]);
        }
    }
    else{
        for (int i = 0 ; i < 6 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(float));
            memcpy(&vaux[i], bar.constData(), sizeof(float));
//            qDebug() << "vaux[" << i << "] = " << vaux[i];
        }
        for (int i = 0 ; i < 6 ; i++){
            griddimensions.append(vaux[i]);
        }
    }
    generategridbounds(vaux);
//      iaux[0] = number of vertices, iaux[2] = number of indices
    for (int i = 0; i < 2 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }
    nvertices = iaux[0];
    nindices = iaux[1];
    QWidget *win;
    win = new QWidget(this);
    win->setAutoFillBackground(true);
    win->setFixedSize(320,50);
    win->setWindowTitle(tr("Loading files"));
    win->raise();
    QVBoxLayout* layout;
    layout = new QVBoxLayout(this);
    QLabel *label;
    label = new QLabel(tr("Loading files"), win);
    QProgressBar *prgbar;
    prgbar = new QProgressBar(win);
    prgbar->resize(300,25);
    prgbar->setOrientation(Qt::Horizontal);	//Orientation can be vertical too
    prgbar->setMinimumWidth(300);
    prgbar->setMaximumWidth(300);
    prgbar->setMinimum(0);
    prgbar->setMaximum(nvertices+nindices);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(prgbar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    VertexNormalData v;
    for(int i = 0 ; i < nvertices ; i++){
        if (fdouble){
            for (int i = 0 ; i < 9 ; i++){
                QByteArray bar;
                bar = inputfile.read(sizeof(double));
                memcpy(&dvaux[i], bar.constData(), sizeof(double));
                vaux[i] = (float)dvaux[i];
            }
        }
        else{
            for (int i = 0 ; i < 9 ; i++){
                QByteArray bar;
                bar = inputfile.read(sizeof(float));
                memcpy(&vaux[i], bar.constData(), sizeof(float));
            }
        }
        v.position.setX(vaux[0]);
        v.position.setY(vaux[1]);
        v.position.setZ(vaux[2]);
        v.normal.setX(vaux[3]);
        v.normal.setY(vaux[4]);
        v.normal.setZ(vaux[5]);
        v.color.setX(vaux[6]);
        v.color.setY(vaux[7]);
        v.color.setZ(vaux[8]);
        v.color.setW(opacity);
        allvertices.append(v);
    }
    for(int i = 0 ; i < nindices/3 ; i++){  // Reads indices in groups of three corresponding to vertices of a triangle
        for (int i = 0; i < 3 ; i++){
            QByteArray bar;
            bar = inputfile.read(sizeof(int));
            memcpy(&iaux[i], bar.constData(), sizeof(int));
        }
        allindices.append(iaux[0]-1);   // Indices coming from file .srf start in 1 (Fortran convention)
        allindices.append(iaux[1]-1);
        allindices.append(iaux[2]-1);
    }
    if (allvertices.count() != nvertices || allindices.count() != nindices){
        QMessageBox msgBox;
        msgBox.setText(tr("readsrf"));
        msgBox.setInformativeText(tr("Wrong number of vertices and allindices from file %1").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    showcolorrule = false;
    localextrema[0].clear();
    QByteArray bar;
    bar = inputfile.read(sizeof(int));
    if (!inputfile.atEnd()){
        memcpy(&iaux[0], bar.constData(), sizeof(int));
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(double));
                    memcpy(&dvaux[i], bar.constData(), sizeof(double));
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(float));
                    memcpy(&vaux[i], bar.constData(), sizeof(float));
                }
            }
            localextrema[0].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[0].append(false);
        }
    }
    localextrema[1].clear();

    bar = inputfile.read(sizeof(int));
    if (!inputfile.atEnd()){
        memcpy(&iaux[0], bar.constData(), sizeof(int));
        for (int i = 0 ; i < iaux[0] ; i++){
            if (fdouble){
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(double));
                    memcpy(&dvaux[i], bar.constData(), sizeof(double));
                    vaux[i] = (float)dvaux[i];
                }
            }
            else{
                for(int i = 0 ; i < 4 ; i++){
                    QByteArray bar;
                    bar = inputfile.read(sizeof(float));
                    memcpy(&vaux[i], bar.constData(), sizeof(float));
                }
            }
            localextrema[1].append(QVector4D(vaux[0],vaux[1],vaux[2],vaux[3]));
            extremactive[1].append(false);
        }
    }
    createballs();
    return true;
}

QColor surface::getextremacolor(int i){
    return extremacolor[i];
}

float surface::getfabstop()
{
    return fabstop;
}

// Function getscalevalueInt
//    Returns slider position ix corresponding to function value fx.
//      i0, i1: range of slider scale
//      f0, f1: function values at i0 and i1
int surface::getscalevalueInt(float fx,int i0,int i1,float f0,float f1)
{
    if (fx > f1) fx = f1;
    if (fx < f0) fx = f0;
    return (int)(((fx-f0)*(i1-i0)/(f1-f0))+i0);
}


// Function getscalevalueFloat
//    Returns function value fx corresponding to slider position ix.
//      i0, i1: range of slider scale
//      f0, f1: function values at i0 and i1
float surface::getscalevalueFloat(int ix,int i0,int i1,float f0,float f1)
{
    if (ix > i1) ix = i1;
    if (ix < i0) ix = i0;
    return (f1-f0) * (ix-i0) / (i1-i0) + f0;
}


//  Creates critical points
void surface::createballs(){
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    VertexNormalData v;
    for (int i = 0 ; i < 2 ; i++){
        allindicesextrema[i].clear();
        allverticesextrema[i].clear();
        allindicesoffsetextrema[i].clear();
        allindicesoffsetextrema[i].append(0);
        color = QVector4D(extremacolor[i].redF(), extremacolor[i].greenF(), extremacolor[i].blueF(), 1.);
        maxindex[i] = 0;
        for (int j = 0 ; j < localextrema[i].length() ; j++ ){
            int kshift = 0;
            if (allindicesextrema[i].length() > 0)
                kshift = maxindex[i];
            for (int k = 0 ; k < sphereindices.length() ; k++){
                allindicesextrema[i] << sphereindices.at(k) + kshift;
                maxindex[i] = std::max(maxindex[i],allindicesextrema[i].last()+1);
            }
            allindicesoffsetextrema[i].append(allindicesextrema[i].length());
            for (int k = 0 ; k < spherevertices.length() ; k++){
                position = SCALE * ballradius * spherevertices.at(k) + localextrema[i].at(j).toVector3D();
                normal = spherevertices.at(k);
                v.position.setX(position.x());
                v.position.setY(position.y());
                v.position.setZ(position.z());
                v.normal.setX(normal.x());
                v.normal.setY(normal.y());
                v.normal.setZ(normal.z());
                v.color = color;
                allverticesextrema[i] << v;
            }
        }
    }
}

// Function makeSphere: generates vertices for drawing a sphere with radius equal to 1 with a given number of
//      "slices" (longitude) and "stacks" (latitude), to be reescaled afterwards
//      as appropriate. Center at (0,0,0).
//
//  Counterclockwise triangles generated
//

// Function makeSphere
void surface::makeSphere(int slices, int stacks){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines
    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = PI / stacks;

    // Generate vertices coordinates, normal values, and texture coordinates

    int numVertices = slices * (stacks - 1) + 2;
    spherevertices.clear();

    // North pole point
    spherevertices << QVector3D(0.,0.,1.);

    // vertices on the main body
    for (int i = 1; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
            spherevertices << QVector3D(sin(deltaLat * i) * cos(deltaLong * j), sin(deltaLat * i) * sin(deltaLong * j),
                                cos(deltaLat * i));
        }
    }

    // South pole point
    spherevertices << QVector3D(0,0,-1);

    // Generates the indices
//    int numIndices = (stacks - 1) * slices * 6; //why multiply by 6?
    sphereindices.clear();

    //add indices in North Pole region (no. of elements is slices * 3)
    for (int j = 1; j<= slices-1; j++){
        sphereindices << 0;
        sphereindices << j;
        sphereindices << j+1;
    }
    sphereindices << 0;
    sphereindices << slices;
    sphereindices << 1;

    //add indices in South Pole Region (no. of element is slices * 3)
    int temp = numVertices  - 1;
    for (int j  = temp-1; j > temp - slices; j--){
        sphereindices << temp;
        sphereindices << j;
        sphereindices << j - 1;
    }
    sphereindices << temp;
    sphereindices << temp-slices;
    sphereindices << temp-1;

    // add body (no. of element is (stacks - 2) * slices * 6
    for (int i = 1; i < stacks - 1; i++) {
        for (int j = 1; j <= slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            sphereindices << (i - 1) * slices + j;
            sphereindices << i * slices + j;
            sphereindices << i * slices + j + 1;
            // triangle two
            sphereindices << (i - 1) * slices + j;
            sphereindices << i * slices + j + 1;
            sphereindices << (i - 1) * slices + j + 1;


        }
        // triangle one
        sphereindices << i * slices;
        sphereindices << (i+1) * slices;
        sphereindices << i * slices + 1;
        // triangle two
        sphereindices << i * slices;
        sphereindices << i * slices + 1;
        sphereindices << (i - 1) * slices + 1;

    }
}
//  End of function makeSphere
//  ---------------------------------------------------------------------------------------------------------------------------

// Changes precision for extrema coordinates values
void surface::setcoordprecision(int i){
    coordprecision = i;
    emit updatedisplay();
}

void surface::setextremactive(int i, int j, bool a){
    extremactive[i][j] = a;
}

void surface::setfullname(QString a){
    fullname = a;
}

void surface::setinitialposition(QPoint a){
    initialposition = a;
}

void surface::selectfontcolor(){
    QColor colact = fontcolor;
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        fontcolor = col;
        emit fontColor(&fontcolor);
        emit updatedisplay();
    }
}

void surface::selectsurfacecolor(){
    QColor colact = surfacecolor;
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        surfacecolor = col;
        surfacecolor_changed(surfacecolor);
        emit surfaceColor(&surfacecolor);
        emit updatedisplay();
    }
}

void surface::selectmaximacolor(){
    QColor colact = extremacolor[0];
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        extremacolor[0] = col;
        createballs();
        emit maximaColor(&extremacolor[0]);
        emit updatedisplay();
    }
}

void surface::selectminimacolor(){
    QColor colact = extremacolor[1];
    QColor col = QColorDialog::getColor(colact, this);
    if(col.isValid()) {
        extremacolor[1] = col;
        createballs();
        emit minimaColor(&extremacolor[1]);
        emit updatedisplay();
    }
}

void surface::setname(QString a){
    name = a;
}

void surface::setonlyextremactive(bool a){
    onlyextremactive = a;
    emit updatedisplay();
}

void surface::setopacity(float a){
    opacity = a;
}

void surface::setshowextremacoords(bool a){
    showextremacoords = a;
    emit updatedisplay();
}

void surface::setshowextremaindices(bool a){
    showextremaindices = a;
    emit updatedisplay();
}

void surface::setshowextremasymbols(bool a){
    showextremasymbols = a;
    emit updatedisplay();
}

void surface::setshowextremavalues(bool a){
    showextremavalues = a;
    emit updatedisplay();
}

void surface::setshowgridbounds(bool a){
    showgridbounds = a;
    emit updatedisplay();
}

void surface::setshowlocalmax(bool a){
    showlocalmax = a;
    emit updatedisplay();
}

void surface::setshowlocalmin(bool a){
    showlocalmin = a;
    emit updatedisplay();
}

void surface::setsolidsurf(bool a){
    solidsurf = a;
    emit updatedisplay();
}

void surface::settopcolor(float a){
    topcolor = a;
}

void surface::settranslucence(bool a){
    translucence = a;
    emit updatedisplay();
}

void surface::settrianglecolors(){
    if (allvalues.count() != allvertices.count()){
        QMessageBox msgBox;
        msgBox.setText(tr("settrianglecolors"));
        msgBox.setInformativeText(tr("Number of vertices %1 different from number of values %2").arg(allvertices.count())
                .arg(allvalues.count()));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    float tc2i = 1.f / (topcolor*topcolor);
    VertexNormalData v;
    for (int i = 0 ; i < allvalues.count() ; i++){
        v = allvertices.at(i);
        if (allvalues.at(i) >= topcolor){
            v.color.setX(1.);
            v.color.setY(0.);
            v.color.setZ(0.);
        }
        else if(allvalues.at(i) <= -topcolor)
        {
            v.color.setX(0.);
            v.color.setY(0.);
            v.color.setZ(1.);
        }
        else if(allvalues.at(i) >= 0.f){
            v.color.setX(allvalues.at(i)/topcolor);
            v.color.setY(sqrtf(1.f-tc2i*allvalues.at(i)*allvalues.at(i)));
            v.color.setZ(0.);
        }
        else {
            v.color.setX(0.);
            v.color.setY(sqrtf(1.f-tc2i*allvalues.at(i)*allvalues.at(i)));
            v.color.setZ(-allvalues.at(i)/topcolor);
        }
        allvertices.replace(i,v);
    }
}

void surface::setvisible(bool a){
    visible = a;
}


// Changes vertical shift for cps labels in drawing
void surface::setvshift(int i){
    vshift = i;
    emit updatedisplay();
}


// Changes precision for extrema coordinates values
void surface::setvalueprecision(int i){
    valueprecision = i;
    emit updatedisplay();
}

void surface::toggleshowsurf(){
    setvisible(!visible);
    emit updatedisplay();
}
