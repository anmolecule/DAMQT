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
//  File:   criticalpoints.cpp
//  Description: implements criticalpoints class
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: July 2018
//
#include "criticalpoints.h"
#include <cmath>


criticalpoints::criticalpoints(QString parentpath, QWidget *parent) : QWidget(parent)
{
    connections.clear();
    QVector<QColor> colors;
    colors << QColor(255,0,0) << QColor(170,85,255) << QColor(195,195,195) << QColor(255,128,0);
    for (int i=0;i<MAX_CPS;i++){
        alleigindices[i].clear();
        alleigvertices[i].clear();
        alleigindicesoffset[i].clear();
        allindices[i].clear();
        allvertices[i].clear();
        allindicesoffset[i].clear();
        cpsactive[i].clear();
        cpscolor[i] = colors[i];
        cpseigval[i].clear();
        cpsxyzval[i].clear();
        setvisiblecps(i,true);
        setvisiblecpseigen(false);
        maxindex[i] = 0;
    }
    cpseigcolor[0] = QColor(255,0,0);
    cpseigcolor[1] = QColor(0,0,255);
    cpseigcolor[2] = QColor(0,200,0);
    font = QFont("Noto Sans", 14, QFont::Bold);
    fontcolor = QColor(255, 255, 0, 255);
    ballradius = 10;
    cpcoordprecision = 2;
    cpeigarrowsize = 4;
    cpeigarrowwidth = 10;
    cpeiglength = 10;
    cpeigthickness = 4;
    cpprecision = 2;
    cpvshift = 0;
    cpsfilename = "";
    name = "";
    ProjectFolder = "";
    drawcpscoords = false;
    drawcpsindices = false;
    drawcpssymbols = false;
    drawcpsvalues = false;
    erroreigvec = false;
    onlycpsactive = false;
    setpath(parentpath);
    coneindices.clear();
    conenormals.clear();
    conevertices.clear();
    cylinderindices.clear();
    cylindervertices.clear();
    spherevertices.clear();
    sphereindices.clear();
    makeSphere(15,15);    // computes vertices of a sphere and its indices
    makeCylinder(15,15);    // computes vertices of a cylinder and its indices
    makeCone(15,15,2.);    // computes vertices of a cone and its indices
}
//	Destructor
criticalpoints::~criticalpoints(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

bool criticalpoints::readcpseigen(QString cp_filename)
{	
    QString filename = QFileInfo(cp_filename).absolutePath() + "/" + QFileInfo(cp_filename).completeBaseName() + ".eigv";
    QFile file(filename);
    erroreigvec = false;
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox msgBox;
        msgBox.setText(tr("readcpseigen"));
        msgBox.setInformativeText(tr("File %1 cannot be read. Hessian eigenvectors cannot be visualized")
                                  .arg(filename)+QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        setvisiblecpseigen(false);
        erroreigvec = true;
        return false;
    }
    for (int i=0;i<MAX_CPS;i++){
        cpseigval[i].clear();
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QString line;
    line = in.readLine();      // The first line contains the number of critical points
    int numCPs = line.toInt();
    int naux = 0;
    for(int i = 0 ; i < MAX_CPS ; ++i){
        naux += cpsxyzval[i].count();
    }
    if (naux != numCPs){
        QMessageBox msgBox;
        msgBox.setText(tr("readcpseigen"));
        msgBox.setInformativeText(tr("Number of eigenvectors in file %1").arg(filename)
                + tr("different than three times number of CPs.\nCheck file.. Hessian eigenvectors cannot be visualized"));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        setvisiblecpseigen(false);
        erroreigvec = true;
        QApplication::restoreOverrideCursor();
        return false;
    }
    QVector4D auxvec;
    bool exitloop = false;
    int knt = 0;
    for(int i = 0 ; i < MAX_CPS && !exitloop ; i++ ){
        cpseigval[i].clear();
        for (int j = 0 ; j < cpsxyzval[i].count() && !exitloop; ++j){
            for (int k = 0 ; k < 3 ; k++){
                if (in.atEnd()){     // Never should happen this if file is correct
                    exitloop = true;
                    break;
                }
                line = in.readLine();
                line=line.trimmed();
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
                if (fields.count() < 7){ 
                    k--;
                    continue;
                }
                fields.removeFirst();   // Neglects the first three entries in each line, which correspond to cp coordinates, already stored
                fields.removeFirst();
                fields.removeFirst();
                auxvec.setX(fields.takeFirst().toFloat());    // Reads the components of Hessian eigenvector and sign of eigenvalue
                auxvec.setY(fields.takeFirst().toFloat());
                auxvec.setZ(fields.takeFirst().toFloat());
                auxvec.setW(fields.takeFirst().toFloat());
                cpseigval[i].append(auxvec);
                knt++;
            }
        }
    }

    if (knt != 3*numCPs){
        QMessageBox msgBox;
        msgBox.setText(tr("readcpseigen"));
        msgBox.setInformativeText(tr("Number of eigenvectors read in file %1 different than 3 times the number of critical points. Hessian eigenvectors cannot be visualized")
                                  .arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        setvisiblecpseigen(false);
        erroreigvec = true;
        QApplication::restoreOverrideCursor();
        return false;
    }
    QApplication::restoreOverrideCursor();
    return true;
}

bool criticalpoints::readcpsfile(QString filename)
{	
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox msgBox;
        msgBox.setText(tr("readcpsfile"));
        msgBox.setInformativeText(tr("File %1 cannot be read").arg(filename)+QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    cpsfilename = filename;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QString line;
    line = in.readLine();      // The first line contains the number of critical points
    int numCPs = line.toInt();
    int itype;
    QStringList CPtypes;
    CPtypes << "x" << "y" << "z" << "m";
    
    for (int i=0;i<MAX_CPS;i++){
        cpsxyzval[i].clear();
    }
    QVector4D cpsvec;
    int knt = 0;
    do{
        line = in.readLine();
        line=line.trimmed();
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (fields.count() < 5) continue;
        QString type(fields.takeFirst());
        itype = CPtypes.indexOf(type.toStdString().c_str());
        if (itype < 0 || itype > 3) continue;
        cpsvec.setX(fields.takeFirst().toDouble() * ANGSTROM_TO_BOHR);
        cpsvec.setY(fields.takeFirst().toDouble() * ANGSTROM_TO_BOHR);
        cpsvec.setZ(fields.takeFirst().toDouble() * ANGSTROM_TO_BOHR);
        cpsvec.setW(fields.takeFirst().toDouble());
        cpsxyzval[itype].append(cpsvec);
        cpsactive[itype].append(false);
        knt++;
        if (knt > numCPs){
            QMessageBox msgBox;
            msgBox.setText(tr("readcpsfile"));
            msgBox.setInformativeText(tr("Number of CPs in file %1 higher than declared on top of file ").arg(filename));
            msgBox.setIcon(QMessageBox::Warning);
            msgBox.exec();
            break;
        }
    }while (!in.atEnd());
    if (knt < numCPs){
        QMessageBox msgBox;
        msgBox.setText(tr("readcpsfile"));
        msgBox.setInformativeText(tr("Number of CPs in file %1 lower than declared on top of file ").arg(filename));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
    }
    QApplication::restoreOverrideCursor();
    return true;
}

// -----------------------------------------------------------------------------------
//                      Functions and SLOTS
// -----------------------------------------------------------------------------------

// Function makeCone: generates vertices for drawing a cone with a given number of "slices" (longitude)
//      and "stacks" (latitude), a radius equal to 1 and a "height"
//      Base centered at (0,0,0) top at (0,0,height)
//
//  Counterclockwise triangles generated
//
void criticalpoints::makeCone(int slices, int stacks, float height){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines

    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = height / (stacks-1);

    conevertices.clear();
    conenormals.clear();

    // vertices on the main body
    float aux = 1. / stacks;
    float bux = 1 / height;
    for (int i = 0; i < stacks-1; i++) {
        float rad = 1. - i * aux;
        float h = deltaLat * i;
        for (int j = 0; j < slices; j++) {
            conevertices << QVector3D(rad*cos(deltaLong * j), rad*sin(deltaLong * j), h);
            conenormals << QVector3D(cos(deltaLong * j), sin(deltaLong * j), bux).normalized();
        }
    }
    conevertices << QVector3D(0,0, height);
    conenormals << QVector3D(0,0,1.);

    // vertices on the base
    int indexbase = conevertices.length();
    conevertices << QVector3D(0,0,0);
    conenormals << QVector3D(0,0,-1.);
    for (int j = 0; j < slices; j++) {
        conevertices << QVector3D(cos(deltaLong * j), sin(deltaLong * j), 0);
        conenormals << QVector3D(0,0,-1.);
    }

    // Generates the indices
    coneindices.clear();
    // indices of the body
    for (int i = 1; i < stacks-1 ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            coneindices << (i - 1) * slices + j;
            coneindices << (i - 1) * slices + j + 1;
            coneindices << i * slices + j;
            // triangle two
            coneindices << (i - 1) * slices + j + 1;
            coneindices << i * slices + j + 1;
            coneindices << i * slices + j;

        }
        // triangle one
        coneindices << i * slices-1;
        coneindices << (i-1) * slices;
        coneindices << (i+1) * slices - 1;
        // triangle two
        coneindices << (i-1) * slices;
        coneindices << i * slices;
        coneindices << (i+1) * slices - 1;
    }
//    triangles of cone vertex
    for (int j = 0; j < slices-1; j++) {
        coneindices << (stacks - 2) * slices + j;
        coneindices << (stacks - 2) * slices + j + 1;
        coneindices << (stacks - 1) * slices;
    }
    coneindices << (stacks - 1) * slices - 1;
    coneindices << (stacks - 2) * slices;
    coneindices << (stacks - 1) * slices;
    // indices of the base
    for (int j = 0; j < slices-1; j++) {
        coneindices << indexbase;
        coneindices << j + 2 + indexbase;
        coneindices << j + 1 + indexbase;
    }
    coneindices << indexbase;
    coneindices << 1 + indexbase;
    coneindices << slices + indexbase;
}

//  End of function makeCone
//  ----------------------------------------------------------------------------------------------------------------------------

// Function makeCylinder: generates vertices for drawing a cylinder side surface (neither top nor base) with a given number of
//      "slices" (longitude) and "stacks" (latitude), a radius equal to 1 and a height equal to 1 to be reescaled afterwards
//      as appropriate. Base center at (0,0,0). Orientation: (0,0,1)
//
//  Counterclockwise triangles generated
//
void criticalpoints::makeCylinder(int slices, int stacks){
    // stacks = no. of Latitude lines,
    // slices = no. of Longitude lines
    double PI = 3.14159265358979;
    double deltaLong = PI * 2 / slices;
    double deltaLat = 1. / (stacks-1);

    cylindervertices.clear();

    // vertices on the main body
    for (int i = 0; i < stacks; i++) {
        for (int j = 0; j < slices; j++) {
            cylindervertices << QVector3D(cos(deltaLong * j), sin(deltaLong * j), deltaLat * i);
        }
    }
    // Generates the indices
    cylinderindices.clear();
    // add body (no. of element is (stacks - 2) * slices * 6
    for (int i = 1; i < stacks ; i++) {
        for (int j = 0; j < slices-1; j++) {
            // each quad gives two triangles
            // triangle one
            cylinderindices << (i - 1) * slices + j;
            cylinderindices << i * slices + j + 1;
            cylinderindices << i * slices + j;

            // triangle two
            cylinderindices << (i - 1) * slices + j;
            cylinderindices << (i - 1) * slices + j + 1;
            cylinderindices << i * slices + j + 1;


        }
        // triangle one
        cylinderindices << i * slices - 1;
        cylinderindices << i * slices;
        cylinderindices << (i+1) * slices - 1;

        // triangle two
        cylinderindices << i * slices - 1;
        cylinderindices << (i - 1) * slices;
        cylinderindices << i * slices;

    }
}

//  End of function makeCylinder
//  ----------------------------------------------------------------------------------------------------------------------------

// Function makeSphere: generates vertices for drawing a sphere with radius equal to 1 with a given number of
//      "slices" (longitude) and "stacks" (latitude), to be reescaled afterwards
//      as appropriate. Center at (0,0,0).
//
//  Counterclockwise triangles generated
//

// Function makeSphere
void criticalpoints::makeSphere(int slices, int stacks){
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

//  Creates critical points
void criticalpoints::createCPs(){
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    VertexNormalData v;
    for (int i = 0 ; i < MAX_CPS ; i++){
        allindices[i].clear();
        allvertices[i].clear();
        allindicesoffset[i].clear();
        allindicesoffset[i].append(0);
        color = QVector4D(cpscolor[i].redF(), cpscolor[i].greenF(), cpscolor[i].blueF(), 1.);
        maxindex[i] = 0;
        for (int j = 0 ; j < cpsxyzval[i].length() ; j++ ){
            int kshift = 0;
            if (allindices[i].length() > 0)
                kshift = maxindex[i];
            for (int k = 0 ; k < sphereindices.length() ; k++){
                allindices[i] << sphereindices.at(k) + kshift;
                maxindex[i] = std::max(maxindex[i],allindices[i].last()+1);
            }
            allindicesoffset[i].append(allindices[i].length());
            for (int k = 0 ; k < spherevertices.length() ; k++){
                position = SCALE * ballradius * spherevertices.at(k) + QVector3D(cpsxyzval[i].at(j));
                normal = spherevertices.at(k);
                v.position.setX(position.x());
                v.position.setY(position.y());
                v.position.setZ(position.z());
                v.normal.setX(normal.x());
                v.normal.setY(normal.y());
                v.normal.setZ(normal.z());
                v.color = color;
                allvertices[i] << v;
            }
        }
    }
}

//  Creates critical points eigenvectors arrows
void criticalpoints::createCPeigarrows(){
    if (erroreigvec){
        setvisiblecpseigen(false);
        return;     // If eigenvectors were not correctly read, cannot create
    }
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    QQuaternion q;
    QVector3D scale = SCALE * QVector3D(cpeigarrowwidth,cpeigarrowwidth,CPSCALEHEIGHT*cpeigarrowsize);
    QVector3D shiftorig = SCALE * QVector3D(0,0,CPSCALEHEIGHT*cpeiglength);
    VertexNormalData v;

    for (int i = 0 ; i < MAX_CPS ; i++){
        int kshift = 0;
        alleigindicesoffset[i].append(alleigindices[i].length());
        for (int j = 0, kj = 0 ; j < cpseigval[i].length() ; j+=3, kj++ ){
            for (int jxyz = 0 ; jxyz < 3 ; jxyz++){
                kshift = maxindex[i];
                if (cpseigval[i].at(j+jxyz).w() > 0){
                    for (int k = 0 ; k < coneindices.length() ; k++){
                        alleigindices[i] << coneindices.at(k) + kshift;
                        maxindex[i] = std::max(maxindex[i],alleigindices[i].last()+1);
                    }
                }
                else{
                    for (int k = 0 ; k < sphereindices.length() ; k++){
                        alleigindices[i] << sphereindices.at(k) + kshift;
                        maxindex[i] = std::max(maxindex[i],alleigindices[i].last()+1);
                    }
                }
                alleigindicesoffset[i].append(alleigindices[i].length());
                q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(cpseigval[i].at(j+jxyz)));
                color = QVector4D(cpseigcolor[jxyz].redF(), cpseigcolor[jxyz].greenF(), cpseigcolor[jxyz].blueF(), 1.);
                if (cpseigval[i].at(j+jxyz).w() > 0){
                    for (int k = 0 ; k < conevertices.length() ; k++){
                        position = q.rotatedVector( scale * conevertices.at(k) + shiftorig) + QVector3D(cpsxyzval[i].at(kj));
                        normal = q.rotatedVector(QVector3D(conenormals.at(k)));
                        v.position.setX(position.x());
                        v.position.setY(position.y());
                        v.position.setZ(position.z());
                        v.normal.setX(normal.x());
                        v.normal.setY(normal.y());
                        v.normal.setZ(normal.z());
                        v.color = color;
                        alleigvertices[i] << v;
                    }
                }
                else{
                    for (int k = 0 ; k < spherevertices.length() ; k++){
                        position = q.rotatedVector(SCALE * cpeigarrowwidth * spherevertices.at(k) + shiftorig)
                                + QVector3D(cpsxyzval[i].at(kj));
                        normal = q.rotatedVector(spherevertices.at(k));
                        v.position.setX(position.x());
                        v.position.setY(position.y());
                        v.position.setZ(position.z());
                        v.normal.setX(normal.x());
                        v.normal.setY(normal.y());
                        v.normal.setZ(normal.z());
                        v.color = color;
                        alleigvertices[i] << v;
                    }
                }
            }
        }
    }
}

//  Creates critical points eigenvectors
void criticalpoints::createCPseigen(){
    if (erroreigvec){
        setvisiblecpseigen(false);
        return;     // If eigenvectors were not correctly read, cannot create
    }
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    QQuaternion q;
    QVector3D scale = SCALE * QVector3D(cpeigthickness,cpeigthickness,CPSCALEHEIGHT*cpeiglength);
    VertexNormalData v;
    for (int i = 0 ; i < MAX_CPS ; i++){
        alleigindices[i].clear();
        alleigvertices[i].clear();
        alleigindicesoffset[i].clear();
        alleigindicesoffset[i].append(0);
        maxindex[i] = 0;
        for (int j = 0, kj = 0 ; j < cpseigval[i].length() ; j+=3, kj++ ){
            int kshift = 0;
            for (int jxyz = 0 ; jxyz < 3 ; jxyz++){
                if (alleigindices[i].length() > 0)
                    kshift = maxindex[i];
                for (int k = 0 ; k < cylinderindices.length() ; k++){
                    alleigindices[i] << cylinderindices.at(k) + kshift;
                    maxindex[i] = std::max(maxindex[i],alleigindices[i].last()+1);
                }
                alleigindicesoffset[i].append(alleigindices[i].length());
                q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),QVector3D(cpseigval[i].at(j+jxyz)));
                color = QVector4D(cpseigcolor[jxyz].redF(), cpseigcolor[jxyz].greenF(), cpseigcolor[jxyz].blueF(), 1.);
                for (int k = 0 ; k < cylindervertices.length() ; k++){
                    position = q.rotatedVector( scale * cylindervertices.at(k) ) + QVector3D(cpsxyzval[i].at(kj));
                    normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
                    v.position.setX(position.x());
                    v.position.setY(position.y());
                    v.position.setZ(position.z());
                    v.normal.setX(normal.x());
                    v.normal.setY(normal.y());
                    v.normal.setZ(normal.z());
                    v.color = color;
                    alleigvertices[i] << v;
                }
            }
        }
    }
    createCPeigarrows();
}

//void criticalpoints::emitupdateRightMenu(){
//    emit updateRightMenu();
//}

bool criticalpoints::getcpsactive(int i, int j){
    if (cpsactive[i].isEmpty() || cpsactive[i].length() < j+1) return false;
    return cpsactive[i].at(j);
}

bool criticalpoints::getdrawcpscoords(){
    return drawcpscoords;
}

bool criticalpoints::getdrawcpsindices(){
    return drawcpsindices;
}

bool criticalpoints::getdrawcpssymbols(){
    return drawcpssymbols;
}

bool criticalpoints::getdrawcpsvalues(){
    return drawcpsvalues;
}

bool criticalpoints::geterroreigvec(){
    return erroreigvec;
}

QColor criticalpoints::getcolor(int i){
    return cpscolor[i];
}

QColor criticalpoints::geteigcolor(int i){
    return cpseigcolor[i];
}

QColor criticalpoints::getfontcolor(){
    return fontcolor;
}

QFont criticalpoints::getfont(){
    return font;
}

bool criticalpoints::getonlycpsactive(){
    return onlycpsactive;
}

QString criticalpoints::getcpsfilename(){
    return cpsfilename;
}

QString criticalpoints::getpath(){
    return path;
}

bool criticalpoints::iscpsactive(int i, int j){
    return cpsactive[i][j];
}

bool criticalpoints::isvisiblecps(int i){
    if (i < 0 || i >= MAX_CPS)
        return false;
    return visiblecps[i];
}

bool criticalpoints::isvisiblecpseigen(){
    return visiblecpseigen;
}

int criticalpoints::getcpcoordprecision(){
    return cpcoordprecision;
}

int criticalpoints::getcpprecision(){
    return cpprecision;
}

int criticalpoints::getcpsactivelength(int i){
    return cpsactive[i].length();
}

int criticalpoints::getcpvshift(){
    return cpvshift;
}

int criticalpoints::geteigarrowsize(){
    return cpeigarrowsize;
}

int criticalpoints::geteigarrowwidth(){
    return cpeigarrowwidth;
}

int criticalpoints::geteiglength(){
    return cpeiglength;
}

int criticalpoints::geteigthickness(){
    return cpeigthickness;
}

QString criticalpoints::getname(){
    return name;
}

int criticalpoints::getradius(){
    return ballradius;
}

void criticalpoints::setcolor(int i, QColor a){
    cpscolor[i] = a;
}

void criticalpoints::setcoordprecision(int i){
    cpcoordprecision = i;
}

void criticalpoints::setcpprecision(int i){
    cpprecision = i;
}

void criticalpoints::setcpsactive(int i, int j, bool a){
    cpsactive[i][j] = a;
} 

void criticalpoints::setcpvshift(int i){
    cpvshift = i;
}

void criticalpoints::setdrawcpscoords(bool a){
    drawcpscoords = a;
}

void criticalpoints::setdrawcpsindices(bool a){
    drawcpsindices = a;
}

void criticalpoints::setdrawcpssymbols(bool a){
    drawcpssymbols = a;
}

void criticalpoints::setdrawcpsvalues(bool a){
    drawcpsvalues = a;
}

void criticalpoints::seteigarrowsize(int i){
    cpeigarrowsize = i;
}

void criticalpoints::seteigarrowwidth(int i){
    cpeigarrowwidth = i;
}

void criticalpoints::seteiglength(int i){
    cpeiglength = i;
}

void criticalpoints::seteigthickness(int i){
    cpeigthickness = i;
}

void criticalpoints::seteigcolor(int i, QColor a){
    cpseigcolor[i] = a;
}

void criticalpoints::setfontcolor(QColor a){
    fontcolor = a;
}

void criticalpoints::setfont(QFont a){
    font = a;
}

void criticalpoints::setname(QString a){
    name = a;
}

void criticalpoints::setonlycpsactive(bool a){
    onlycpsactive = a;
}

//  Sets path to file with critical points
void criticalpoints::setpath(QString a){
    path = a;
}

//  Sets Project Folder
void criticalpoints::set_ProjectFolder(QString a){
    ProjectFolder = a;
}

void criticalpoints::setradius(int i){
    ballradius = i;
    createCPs();
}

//  Sets cps show/hide
void criticalpoints::setvisiblecps(int i, bool a){
    if (i < 0 || i >= MAX_CPS)
        return;
    visiblecps[i] = a;
}

//  Sets cps eigenvectors show/hide
void criticalpoints::setvisiblecpseigen(bool a){
    visiblecpseigen = a;
}
