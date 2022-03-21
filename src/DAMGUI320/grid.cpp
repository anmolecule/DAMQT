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
//	Implementation of class grid
//
//	File:   grid.cpp
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: April 2018
//
#include <QDialog>
#include <QLabel>
#include <QMessageBox>
#include <QHBoxLayout>
#include <QVBoxLayout>

#include <QFile>
#include <QTextStream>

#include <QtDebug>

#include <float.h>     // For FLT_MAX and FLT_MIN  (highest and lowest float numbers
#include <cmath>       // std::abs

#include "grid.h"

grid::grid(QWidget *parent) : QWidget(parent)
{
    initialposition = QPoint(100,500);
    surfaces = new QList<isosurface*>();
    compatderiv = false;
    fun   = nullpointer;
    dxfun = nullpointer;
    dyfun = nullpointer;
    dzfun = nullpointer;
    cisosurface = nullpointer;
    sup_corner_max = QVector3D(0,0,0);
    sup_corner_min = QVector3D(0,0,0);
    nameindex = 0;
    surfcolors << QColor(245,0,0) << QColor(0,0,245) << QColor(0,175,0) <<
                QColor(255,128,0) << QColor(128,0,128) << QColor(0,128,128) <<
                QColor(255,255,0) << QColor(128,128,0) << QColor(255,0,255) <<
                QColor(0,255,255) << QColor(0,170,0) << QColor(0,0,170) <<
                QColor(170,0,170) << QColor(255,220,168) << QColor(192,192,192) <<
                QColor(187,5,27) << QColor(128,0,0) << QColor(0,128,0) <<
                QColor(0,128,128) << QColor(128,0,128) << QColor(128,128,0) <<
                QColor(128,128,128) << QColor(255,255,255) << QColor(0,0,0);
}

grid::~grid(){
    if (fun){
        delete fun;
        fun = nullpointer;
    }
    if (dxfun){
        delete dxfun;
        dxfun = nullpointer;
    }
    if (dyfun){
        delete dyfun;
        dyfun = nullpointer;
    }
    if (dzfun){
        delete dzfun;
        dzfun = nullpointer;
    }
    if (surfaces){
        for (int i = surfaces->count()-1 ; i >= 0  ; i--){
            delete surfaces->at(i);
//            surfaces->replace(i,nullpointer);
        }
        surfaces->clear();
        delete surfaces;
        surfaces = nullpointer;
    }
    if (cisosurface){
        delete cisosurface;
        cisosurface = nullpointer;
    }
}

QString grid::getfullname(){
    return fullname;
}

QString grid::getname(){
    return name;
}

void grid::addisosurf(){
    surfaces->append(new isosurface());
    surfaces->last()->set_ProjectFolder(ProjectFolder);
    surfaces->last()->set_ProjectName(ProjectName);
    surfaces->last()->setname(name.remove(".plt")+QString(tr("_surf_%1")).arg(nameindex));
    surfaces->last()->setfullname(fullname+QString(tr("_surf_%1")).arg(nameindex++));
    surfaces->last()->setmaxcontourvalue(maxcontourvalue);
    surfaces->last()->setmincontourvalue(mincontourvalue);
    surfaces->last()->setinitialposition(QPoint(getinitialposition())+QPoint(20,100)*(surfaces->count()-1)+QPoint(0,10));
    surfaces->last()->setsurfcolor(this->surfcolors.at((nameindex-1)%surfcolors.count()));
    surfaces->last()->setcompatderiv(compatderiv);
    emit surfaceadded();
}

void grid::deletesurf(int i){
    // All surfaces editors must be closed to prevent crash
    for (int j = 0 ; j < surfaces->length() ; j++){
        surfaces->at(j)->closeeditor();
    }
    delete surfaces->at(i);
    surfaces->removeAt(i);   
    emit surfacedeleted();
}


void grid::generatesurf(int i){
    if (cisosurface){
        delete cisosurface;
        cisosurface = nullpointer;
    }
    cisosurface = new CIsoSurface<float>();
    float cell_sizes[3] = {	fun->voxel_x() / (fun->x()-1), fun->voxel_y() / (fun->y()-1), fun->voxel_z() / (fun->z()-1)};
    if (compatderiv && surfaces->at(i)->getnormalgrad()){
        cisosurface->GenerateSurfacewithgrad(fun->data(), dxfun->data(), dyfun->data(), dzfun->data(), surfaces->at(i)->getcontourvalue(),
                    fun->x()-1, fun->y()-1, fun->z()-1, cell_sizes[0], cell_sizes[1], cell_sizes[2]);
    }
    else{
        cisosurface->GenerateSurface(fun->data(), surfaces->at(i)->getcontourvalue(),
                    fun->x()-1, fun->y()-1, fun->z()-1, cell_sizes[0], cell_sizes[1], cell_sizes[2]);
    }
    if (!cisosurface->IsSurfaceValid()){
        QMessageBox msgBox;
        msgBox.setText(tr("generatesurf"));
        msgBox.setInformativeText(tr("Invalid surface"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return;
    }
    else{
//QFile file("surface");
//file.open(QFile::WriteOnly | QFile::Text);
//QTextStream out(&file); // Buffer for writing to file
        float despl[3]={sup_corner_min.x(),sup_corner_min.y(),sup_corner_min.z()};
        if (cisosurface->m_nVertices != cisosurface->m_nNormals){
            QMessageBox msgBox;
            msgBox.setText(tr("generatesurf"));
            msgBox.setInformativeText(tr("Number of vertices (%1)").arg(cisosurface->m_nVertices)
                + QString(tr("\ndoes not coincide with the number of normals (%1).")).arg(cisosurface->m_nNormals));
            msgBox.setIcon(QMessageBox::Warning);
            return;
        }
        surfaces->at(i)->allvertices.clear();
        QColor surfcolor = surfaces->at(i)->getsurfcolor();
        float opacity = surfaces->at(i)->getopacity();
        VertexNormalData v;
        for (unsigned int j = 0 ; j < cisosurface->m_nVertices ; j++){
            v.position.setX(cisosurface->m_ppt3dVertices[j][0]+despl[0]);
            v.position.setY(cisosurface->m_ppt3dVertices[j][1]+despl[1]);
            v.position.setZ(cisosurface->m_ppt3dVertices[j][2]+despl[2]);
            v.normal.setX(cisosurface->m_pvec3dNormals[j][0]);
            v.normal.setY(cisosurface->m_pvec3dNormals[j][1]);
            v.normal.setZ(cisosurface->m_pvec3dNormals[j][2]);
//out << QString("%1  %2  %3  %4  %5  %6\n").arg(v.position.x()).arg(v.position.y()).arg(v.position.z())
//.arg(v.normal.x()).arg(v.normal.y()).arg(v.normal.z());
            v.color.setX(surfcolor.redF());
            v.color.setY(surfcolor.greenF());
            v.color.setZ(surfcolor.blueF());
            v.color.setW(opacity);
            surfaces->at(i)->allvertices.append(v);
        }
//file.close();
        surfaces->at(i)->allindices.clear();
        for (unsigned int j = 0 ; j < 3*cisosurface->m_nTriangles ; j++){
            surfaces->at(i)->allindices.append(cisosurface->m_piTriangleIndices[j]);
        }
        float vaux[6] = {sup_corner_min[0],sup_corner_max[0],sup_corner_min[1],sup_corner_max[1],
                         sup_corner_min[2],sup_corner_max[2]};
        surfaces->at(i)->generategridbounds(vaux);
    }
}

QPoint grid::getinitialposition(){
    return initialposition;
}

float grid::getmaxcontourvalue(){
    return maxcontourvalue;
}

float grid::getmincontourvalue(){
    return mincontourvalue;
}

bool grid::loadnormals(FILE *f1, FILE *f2, FILE *f3, VVBuffer*v1, VVBuffer*v2, VVBuffer*v3,
        int *iref, float *vref, int kntbar, QProgressBar *bar){
    bool fdouble;
    int iaux[3];
//    Reads header of file f1 and checks compatibility
    fread( iaux , sizeof(int) , 2 , f1);     // reads two integer values
    if (iaux[0] == 0){   // if the first one is zero: grid data in double precision
        fdouble = true;
    }
    else{                // elsel grid data: float
        fdouble = false;
    }
    fread( iaux , sizeof(int) , 3 , f1);     // reads nz, ny , nx (in this order)
    if (iaux[0] != iref[0] || iaux[1] != iref[1] || iaux[2] != iref[2])
        return false;
//    Reads header of file f2 and checks compatibility
    fread( iaux , sizeof(int) , 2 , f2);     // reads two integer values
    if ((iaux[0] == 0 && !fdouble) || (iaux[0] != 0 && fdouble)){   // if the first one is zero: grid data in double precision
        return false;
    }
    fread( iaux , sizeof(int) , 3 , f2);     // reads nz, ny , nx (in this order)
    if ((iaux[0] != iref[0]) || (iaux[1] != iref[1]) || (iaux[2] != iref[2]))
        return false;
//    Reads header of file f3 and checks compatibility
    fread( iaux , sizeof(int) , 2 , f3);     // reads two integer values
    if ((iaux[0] == 0 && !fdouble) || (iaux[0] != 0 && fdouble)){   // if the first one is zero: grid data in double precision
        return false;
    }
    fread( iaux , sizeof(int) , 3 , f3);     // reads nz, ny , nx (in this order)
    if ((iaux[0] != iref[0]) || (iaux[1] != iref[1]) || (iaux[2] != iref[2]))
        return false;

    float vaux[6];
    if (fdouble){
        double dvaux[6];
//      Reads grid dimensions from file f1 and checks compatibility
        fread(dvaux , sizeof(double) , 6 , f1);
        for(int i = 0 ; i < 6 ; i++){
            vaux[i] = (float)dvaux[i];
        }
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
//      Reads grid dimensions from file f2 and checks compatibility
        fread(dvaux , sizeof(double) , 6 , f2);
        for(int i = 0 ; i < 6 ; i++){
            vaux[i] = (float)dvaux[i];
        }
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
//      Reads grid dimensions from file f3 and checks compatibility
        fread(dvaux , sizeof(double) , 6 , f3);
        for(int i = 0 ; i < 6 ; i++){
            vaux[i] = (float)dvaux[i];
        }
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
//      Reads gradient components, normalizes gradient and stores normalized components
        float *dp1 = v1->data();
        float *dp2 = v2->data();
        float *dp3 = v3->data();
        double *buff1 = new double[nx];
        double *buff2 = new double[nx];
        double *buff3 = new double[nx];
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                fread(buff1, sizeof(double), nx, f1);
                fread(buff2, sizeof(double), nx, f2);
                fread(buff3, sizeof(double), nx, f3);
                for (int k = 0; k < nx; k++) {
                    double normalization = sqrt(buff1[k]*buff1[k]+buff2[k]*buff2[k]+buff3[k]*buff3[k]);
                    if (normalization > 0.)
                        normalization = 1. / normalization;
                    *dp1++ = (float)(buff1[k] * normalization);
                    *dp2++ = (float)(buff2[k] * normalization);
                    *dp3++ = (float)(buff3[k] * normalization);
                }
                bar->setValue(i*nx+j*nx*ny + kntbar*nx*ny*nz);
            }
        }
        return true;
    }
    else{
//      Reads grid dimensions from file f1 and checks compatibility
        fread(vaux , sizeof(float) , 6 , f1);
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
//      Reads grid dimensions from file f2 and checks compatibility
        fread(vaux , sizeof(float) , 6 , f2);
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
//      Reads grid dimensions from file f3 and checks compatibility
        fread(vaux , sizeof(float) , 6 , f3);
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
//      Reads gradient components, normalizes gradient and stores normalized components         ;
        float *dp1 = v1->data();
        float *dp2 = v2->data();
        float *dp3 = v3->data();
        float *buff1 = new float[nx];
        float *buff2 = new float[nx];
        float *buff3 = new float[nx];
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                fread(buff1, sizeof(float), nx, f1);
                fread(buff2, sizeof(float), nx, f2);
                fread(buff3, sizeof(float), nx, f3);
                for (int k = 0; k < nx; k++) {
                    double normalization = sqrt(buff1[k]*buff1[k]+buff2[k]*buff2[k]+buff3[k]*buff3[k]);
                    if (normalization > 0.)
                        normalization = 1. / normalization;
                    *dp1++ = (float)(buff1[k] * normalization);
                    *dp2++ = (float)(buff2[k] * normalization);
                    *dp3++ = (float)(buff3[k] * normalization);
                }
                bar->setValue(i*nx+j*nx*ny + kntbar*nx*ny*nz);
            }
        }
        return true;
    }
}

bool grid::loadderiv(FILE *f, VVBuffer*v, int *iref, float *vref, int kntbar, QProgressBar *bar){
    bool fdouble;
    int iaux[3];
    fread( iaux , sizeof(int) , 2 , f);     // reads two integer values
    if (iaux[0] == 0){   // if the first one is zero: grid data in double precision
        fdouble = true;
    }
    else{                // elsel grid data: float
        fdouble = false;
    }
    fread( iaux , sizeof(int) , 3 , f);     // reads nz, ny , nx (in this order)
    if (iaux[0] != iref[0] || iaux[1] != iref[1] || iaux[2] != iref[2])
        return false;
    float vaux[6];
    if (fdouble){
        double dvaux[6];
        fread(dvaux , sizeof(double) , 6 , f);
        for(int i = 0 ; i < 6 ; i++){
            vaux[i] = (float)dvaux[i];
        }
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
        float *dp = v->data();
        double *buff=new double[nx];
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                fread(buff, sizeof(double), nx, f);
                for (int k = 0; k < nx; k++) {
                        *dp++ = (float)buff[k];
                }
                bar->setValue(i*nx+j*nx*ny + kntbar*nx*ny*nz);
            }
        }
        return true;
    }
    else{
        fread(vaux , sizeof(float) , 6 , f);
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
        float *dp = v->data();
        float *buff=new float[nx];
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                fread(buff, sizeof(float), nx, f);
                for (int k = 0; k < nx; k++) {
                        *dp++ = (float)buff[k];
                }
                bar->setValue(i*nx+j*nx*ny + kntbar*nx*ny*nz);
            }
        }
        return true;
    }
}

bool grid::loadderivnew(QFile *inputfile, VVBuffer*v, int *iref, float *vref, int kntbar, QProgressBar *bar){
    bool fdouble;
    int iaux[3];
    QDataStream data(inputfile);
    for (int i = 0 ; i < 2 ; i++){
        QByteArray bar;
        bar = inputfile->read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "en loadderivnew: iaux[" << i << "] = " << iaux[i];
    }
    if (iaux[0] == 0){   // if the first one is zero: grid data in double precision
        fdouble = true;
    }
    else{                // elsel grid data: float
        fdouble = false;
    }
    for (int i = 0 ; i < 3 ; i++){
        QByteArray bar;
        bar = inputfile->read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }
    if (iaux[0] != iref[0] || iaux[1] != iref[1] || iaux[2] != iref[2])
        return false;
    float vaux[6];
    QByteArray bytar;
    if (fdouble){
        double dvaux[6];
        for (int i = 0 ; i < 6 ; i++){
            bytar = inputfile->read(sizeof(double));
            memcpy(&dvaux[i], bytar.constData(), sizeof(double));
            vaux[i] = (float)dvaux[i];
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
        float *dp = v->data();
        double val;
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                for (int k = 0; k < nx; k++){
                    bytar = inputfile->read(sizeof(double));
                    memcpy(&val, bytar.constData(), sizeof(double));
                    *dp++ = (float)val;
                }
                bar->setValue(i*nx+j*nx*ny + kntbar*nx*ny*nz);
            }
        }
        return true;
    }
    else{
        for (int i = 0 ; i < 6 ; i++){
            bytar = inputfile->read(sizeof(float));
            memcpy(&vaux[i], bytar.constData(), sizeof(double));
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
        if (std::abs(vaux[0]-vref[0]) + std::abs(vaux[1]-vref[1]) + std::abs(vaux[2]-vref[2]) + std::abs(vaux[3]-vref[3])
                + std::abs(vaux[4]-vref[4]) + std::abs(vaux[5]-vref[5]) > 1.e-5){
            return false;
        }
        float *dp = v->data();
        float val;
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                for (int k = 0; k < nx; k++){
                    bytar = inputfile->read(sizeof(float));
                    memcpy(&val, bytar.constData(), sizeof(float));
                    *dp++ = val;
                }
                bar->setValue(i*nx+j*nx*ny + kntbar*nx*ny*nz);
            }
        }
        return true;
    }
}

bool grid::readplt(QString filename){
    const double factor=0.529177249; // units conversion factor

    int  iaux[3]; //iaux[0]=z, iaux[1]=y, iaux[2]=x
    bool fdouble;
    bool existderivs;
    FILE *f, *fdx, *fdy, *fdz;
    QByteArray ba = filename.toLatin1();
    ba = ba.remove(ba.size()-4,4);
    QByteArray badx = ba+QByteArray("-dx.pltd");
    QByteArray bady = ba+QByteArray("-dy.pltd");
    QByteArray badz = ba+QByteArray("-dz.pltd");
    ba = ba.append(".plt");
    const char *file = ba.data();
    const char *filedx = badx.data();
    const char *filedy = bady.data();
    const char *filedz = badz.data();

//    Open files with function and existderivs
    f   = fopen(file   , "rb" );
    fdx = fopen(filedx , "rb" );
    fdy = fopen(filedy , "rb" );
    fdz = fopen(filedz , "rb" );
    existderivs = true;
    compatderiv = true;
//    Checks if the files exist
    if (fdx==nullpointer || fdy==nullpointer || fdz==nullpointer) {
        existderivs = false;
        compatderiv = false;
        if (fdx != nullpointer) fclose(fdx);
        if (fdy != nullpointer) fclose(fdy);
        if (fdz != nullpointer) fclose(fdz);
    }
//    Read file with function
    fread( iaux , sizeof(int) , 2 , f); //reads two integers
    if (iaux[0] == 0)   // if the first one is zero: grid data in double precision
        fdouble = true;
    else                // elsel grid data: float
        fdouble = false;
    fread( iaux , sizeof(int) , 3 , f); //reads nz, ny , nx (in this order)
    nx=iaux[2]; ny=iaux[1]; nz=iaux[0];
    if ( nx < 0 || ny < 0 || nz < 0) {
        QMessageBox msgBox;
        msgBox.setText(tr("readplt"));
        msgBox.setInformativeText(tr("Error: wrong dimensions"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }

    if (!fun)   fun   = new VVBuffer(file, nx, ny, nz);
    if (!dxfun) dxfun = new VVBuffer(filedx, nx, ny, nz);
    if (!dyfun) dyfun = new VVBuffer(filedy, nx, ny, nz);
    if (!dzfun) dzfun = new VVBuffer(filedz, nx, ny, nz);

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
    if (existderivs)
        bar->setMaximum(4*nx*ny*nz);
    else
        bar->setMaximum(nx*ny*nz);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(bar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    float *dp = fun->data();
    float min = FLT_MAX, max = FLT_MIN;
    float vaux[6];
    if (fdouble){
        double dvaux[6];
        fread(dvaux , sizeof(double) , 6 , f);
        for(int i = 0 ; i < 6 ; i++){
            vaux[i] = (float)dvaux[i];
        }
        fun->set_voxel_size((vaux[5]-vaux[4])/factor, (vaux[3]-vaux[2])/factor, (vaux[1]-vaux[0])/factor);
        double *buff=new double[nx];
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                    fread(buff, sizeof(double), nx, f);
                    for (int k = 0; k < nx; k++) {
                            if (max < buff[k]) max = buff[k];
                            if (min > buff[k]) min = buff[k];
                            *dp++ = (float)buff[k];
                    }
                    bar->setValue(i*nx+j*nx*ny);
            }
        }
    }
    else{
        fread(vaux , sizeof(float) , 6 , f);
        fun->set_voxel_size((vaux[5]-vaux[4])/factor, (vaux[3]-vaux[2])/factor, (vaux[1]-vaux[0])/factor);
        float *buff=new float[nx];
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                fread(buff, sizeof(float), nx, f);
                for (int k = 0; k < nx; k++) {
                        if (max < buff[k]) max = buff[k];
                        if (min > buff[k]) min = buff[k];
                        *dp++ = buff[k];
                }
                bar->setValue(i*nx+j*nx*ny);
            }
        }
    }
    fclose(f);
    sup_corner_min.setX(vaux[4]/factor);
    sup_corner_min.setY(vaux[2]/factor);
    sup_corner_min.setZ(vaux[0]/factor);
    sup_corner_max.setX(vaux[5]/factor);
    sup_corner_max.setY(vaux[3]/factor);
    sup_corner_max.setZ(vaux[1]/factor);
    setmaxcontourvalue(max);
    setmincontourvalue(min);

//    Read files with function derivatives (gradient), and computes and stores them for normals interpolation
    if (existderivs){
        compatderiv = loadderiv(fdx, dxfun, iaux, vaux, 1, bar);
        if (compatderiv) compatderiv = loadderiv(fdy, dyfun, iaux, vaux, 2, bar);
        if (compatderiv) compatderiv = loadderiv(fdz, dzfun, iaux, vaux, 3, bar);

//        compatderiv = loadnormals(fdx, fdy, fdz, dxfun, dyfun, dzfun, iaux, vaux, 1, bar);    // alternative: loads normalized gradient
        fclose(fdx);
        fclose(fdy);
        fclose(fdz);
        if (!compatderiv){
            QMessageBox msgBox;
            msgBox.setText(tr("readplt"));
            msgBox.setInformativeText(tr("Files with derivatives are not compatible with function file.\n")
                                      +tr("Computes gradient numerically."));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
        }
    }
    return true;
}

bool grid::readpltnew(QString fileName){
    const double factor=0.529177249; // units conversion factor

    int  iaux[3]; //iaux[0]=z, iaux[1]=y, iaux[2]=x
    bool fdouble;
    bool existderivs;
    QString filename = fileName;
    QFile inputfile(filename);
    QFile inputfiledx(filename.remove(".plt")+"-dx.pltd");
    QFile inputfiledy(filename+"-dy.pltd");
    QFile inputfiledz(filename+"-dz.pltd");
//    qDebug() << "inputfile = " << filename;
//    qDebug() << "inputfiledx = " << filename+"-dx.pltd";
//    qDebug() << "inputfiledy = " << filename+"-dy.pltd";
//    qDebug() << "inputfiledz = " << filename+"-dz.pltd";

    compatderiv = true;

    if (!inputfile.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename+".plt";
        return false;
    }
//    qDebug() << filename+".plt" << "opened";
    existderivs = true;
    if (!inputfiledx.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename+"-dx.pltd";
        existderivs = false;
    }
//    qDebug() << filename+"-dx.pltd" << "opened";
    if (!inputfiledy.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename+"-dy.pltd";
        existderivs = false;
    }
//    qDebug() << filename+"-dy.pltd" << "opened";
    if (!inputfiledz.open(QIODevice::ReadOnly)){
//        qDebug() << "No puede abrir " << filename+"-dz.pltd";
        existderivs = false;
    }
//    qDebug() << filename+"-dz.pltd" << "opened";

//    Read file with function

    QDataStream data(&inputfile);
    for (int i = 0 ; i < 2 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }
    if (iaux[0] == 0)   // if the first one is zero: grid data in double precision
        fdouble = true;
    else                // elsel grid data: float
        fdouble = false;

    for (int i = 0 ; i < 3 ; i++){
        QByteArray bar;
        bar = inputfile.read(sizeof(int));
        memcpy(&iaux[i], bar.constData(), sizeof(int));
//        qDebug() << "iaux[" << i << "] = " << iaux[i];
    }
    nx=iaux[2]; ny=iaux[1]; nz=iaux[0];
    if ( nx < 0 || ny < 0 || nz < 0) {
        QMessageBox msgBox;
        msgBox.setText(tr("readplt"));
        msgBox.setInformativeText(tr("Error: wrong dimensions"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
//    qDebug() << "nx = " << nx << "ny = " << ny << "nz = " << nz;

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
    if (existderivs)
        bar->setMaximum(4*nx*ny*nz);
    else
        bar->setMaximum(nx*ny*nz);
    layout->addWidget(label,Qt::AlignCenter);
    layout->addWidget(bar,Qt::AlignCenter);
    win->setLayout(layout);
    win->show();
    const char *file = filename.toUtf8()+".plt";
    if (!fun)   fun   = new VVBuffer(file, nx, ny, nz);
    float *dp = fun->data();
    float min = FLT_MAX, max = FLT_MIN;
    float vaux[6];
    QByteArray bytar;
    if (fdouble){
        double dvaux[6];
        for (int i = 0 ; i < 6 ; i++){
            bytar = inputfile.read(sizeof(double));
            memcpy(&dvaux[i], bytar.constData(), sizeof(double));
            vaux[i] = (float)dvaux[i];
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
    }else{
        QByteArray bytar;
        for (int i = 0 ; i < 6 ; i++){
            bytar = inputfile.read(sizeof(float));
            memcpy(&vaux[i], bytar.constData(), sizeof(double));
//            qDebug() << "dvaux[" << i << "] = " << dvaux[i];
        }
    }
    fun->set_voxel_size((vaux[5]-vaux[4])/factor, (vaux[3]-vaux[2])/factor, (vaux[1]-vaux[0])/factor);
    if (fdouble){
        double val;
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                for (int k = 0; k < nx; k++){
                    bytar = inputfile.read(sizeof(double));
                    memcpy(&val, bytar.constData(), sizeof(double));
                    *dp++ = (float)val;
                    if (max < val) max = val;
                    if (min > val) min = val;
                }
                bar->setValue(i*nx+j*nx*ny);
            }
        }
    }
    else{
        float val;
        for (int j = 0; j < nz; j++){
            for (int i = 0; i < ny; i++) {
                for (int k = 0; k < nx; k++){
                    bytar = inputfile.read(sizeof(float));
                    memcpy(&val, bytar.constData(), sizeof(float));
                    *dp++ = val;
                    if (max < val) max = val;
                    if (min > val) min = val;
                }
                bar->setValue(i*nx+j*nx*ny);
            }
        }
    }
//    qDebug() << "fun->data()";
//    for (int i = 0 ; i < 10 ; i++){
//        qDebug() << fun->data()[i];
//    }
    sup_corner_min.setX(vaux[4]/factor);
    sup_corner_min.setY(vaux[2]/factor);
    sup_corner_min.setZ(vaux[0]/factor);
    sup_corner_max.setX(vaux[5]/factor);
    sup_corner_max.setY(vaux[3]/factor);
    sup_corner_max.setZ(vaux[1]/factor);
    setmaxcontourvalue(max);
    setmincontourvalue(min);

    compatderiv = existderivs;

//    Read files with function derivatives (gradient), and computes and stores them for normals interpolation
    if (existderivs){
        const char *filedx = filename.toUtf8()+"-dx.pltd";
        if (!dxfun)   dxfun   = new VVBuffer(filedx, nx, ny, nz);
        compatderiv = loadderivnew(&inputfiledx, dxfun, iaux, vaux, 1, bar);
        if (compatderiv){
            const char *filedy = filename.toUtf8()+"-dy.pltd";
            if (!dyfun)   dyfun   = new VVBuffer(filedy, nx, ny, nz);
            compatderiv = loadderivnew(&inputfiledy, dyfun, iaux, vaux, 1, bar);
        }
        if (compatderiv){
            const char *filedz = filename.toUtf8()+"-dz.pltd";
            if (!dzfun)   dzfun   = new VVBuffer(filedz, nx, ny, nz);
            compatderiv = loadderivnew(&inputfiledz, dzfun, iaux, vaux, 1, bar);
        }
    }

    return true;
}

void grid::setinitialposition(QPoint a){
    initialposition = a;
}

void grid::setmaxcontourvalue(float a){
    maxcontourvalue = a;
}

void grid::setmincontourvalue(float a){
    mincontourvalue = a;
}

void grid::setfullname(QString a){
    fullname = a;
}

void grid::setname(QString a){
    name = a;
}

void grid::set_ProjectFolder(QString name){
    ProjectFolder = name;
}

void grid::set_ProjectName(QString name){
    ProjectName = name;
}

void grid::toggleshowsurf(int i){
    surfaces->at(i)->toggleshowsurf();
}
