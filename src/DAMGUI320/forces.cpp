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
//  File:   forces.cpp
//  Description: implements forces class
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: February 2019
//
#include "forces.h"
#include <cmath>


forces::forces(QString parentpath, QWidget *parent) : QWidget(parent)
{
    connections.clear();
    QVector<QColor> colors;
    colors << QColor(183,0,0) << QColor(0,0,183) << QColor(0,175,0) <<
            QColor(255,128,0) << QColor(128,0,128) << QColor(0,128,128) <<
            QColor(255,255,0) << QColor(128,128,0) << QColor(255,0,255) <<
            QColor(0,255,255) << QColor(0,170,0) << QColor(0,0,170) <<
            QColor(170,0,170) << QColor(255,220,168) << QColor(192,192,192) <<
            QColor(187,5,27) << QColor(128,0,0) << QColor(0,128,0) <<
            QColor(0,128,128) << QColor(128,0,128) << QColor(128,128,0) <<
            QColor(128,128,128) << QColor(255,255,255) << QColor(0,0,0);
    for (int i=0;i<MAX_FORCES;i++){
        allindices.clear();
        allvertices.clear();
        allindicesoffset.clear();
        forcecolors[i] = colors[i];
        visibleforces[i] = true;
    }
    forcesfilename = "";
    name = "";
    ProjectFolder = "";
    errorforces = false;
    visible = true;
    setpath(parentpath);
    coneindices.clear();
    conenormals.clear();
    conevertices.clear();
    cylinderindices.clear();
    cylindervertices.clear();
    forcesorig->clear();
    forcesxyz->clear();
    arrowlength = 4;
    arrowwidth = 10;
    forceslength = 15;
    forcesthickness = 4;
    makeCylinder(15,15);    // computes vertices of a cylinder and its indices
    makeCone(15,15,2.);    // computes vertices of a cone and its indices
}
//	Destructor
forces::~forces(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
}

bool forces::readforcesfile(QString filename)
{	
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox msgBox;
        msgBox.setText(tr("readforcesfile"));
        msgBox.setInformativeText(tr("File %1 cannot be read").arg(filename)+QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Warning);
        msgBox.exec();
        return false;
    }
    forcesfilename = filename;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QString line;
    line = in.readLine();      // The first line contains the number of critical points
    numcenters = line.toInt();
    QVector3D auxvec;
    for (int i=0 ; i<MAX_FORCES ; i++){
        forcesval[i].clear();
        forcesorig[i].clear();
        forcesxyz[i].clear();
        for (int j=0 ; j < numcenters ; j++){
            line = in.readLine();
            line=line.trimmed();
#if QT_VERSION < 0x050E00
        QStringList fields = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList fields = line.split(' ',Qt::SkipEmptyParts);
#endif
            if (fields.count() != 6){
                QMessageBox msgBox;
                msgBox.setText(tr("readforcesfile"));
                msgBox.setInformativeText(tr("Error reading file %1: insufficient number of entries in forces").arg(filename));
                msgBox.setIcon(QMessageBox::Warning);
                msgBox.exec();
                return false;
            }
            auxvec.setX(fields.takeFirst().toDouble() * ANGSTROM_TO_BOHR);
            auxvec.setY(fields.takeFirst().toDouble() * ANGSTROM_TO_BOHR);
            auxvec.setZ(fields.takeFirst().toDouble() * ANGSTROM_TO_BOHR);
            forcesorig[i].append(auxvec);
            auxvec.setX(fields.takeFirst().toDouble());
            auxvec.setY(fields.takeFirst().toDouble());
            auxvec.setZ(fields.takeFirst().toDouble());
            forcesxyz[i].append(auxvec);
            forcesval[i].append(auxvec.length());
        }
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
void forces::makeCone(int slices, int stacks, float height){
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
void forces::makeCylinder(int slices, int stacks){
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

// Creates forces

void forces::createforces(){
    QVector3D position;
    QVector3D normal;
    QVector4D color;
    QQuaternion q;
    QVector3D scale;
    VertexNormalData v;
    allindices.clear();
    allvertices.clear();
    allindicesoffset.clear();
    allindicesoffset.append(0);
//    Forces bodies
    uint maxindex = 0;
    uint kshift = 0;
    for (int i=0 ; i<MAX_FORCES ; i++){
        if (!visibleforces[i])
            continue;
        for (int j=0 ; j < numcenters ; j++){
            for (int k = 0 ; k < cylinderindices.length() ; k++){
                allindices << cylinderindices.at(k) + kshift;
                maxindex = std::max(maxindex,allindices.last()+1);
            }
            allindicesoffset.append(allindices.length());
            kshift = maxindex;
            QVector3D fxyz = forcesxyz[i].at(j);
            q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),fxyz);
            QColor fcolor = forcecolors[i];
            color = QVector4D(fcolor.redF(), fcolor.greenF(), fcolor.blueF(), 1.);
            scale = SCALE * QVector3D(forcesthickness, forcesthickness, FORCESCALEHEIGHT * fxyz.length() * forceslength);
            for (int k = 0 ; k < cylindervertices.length() ; k++){
                position = q.rotatedVector( scale * cylindervertices.at(k)) + QVector3D(forcesorig[i].at(j));
                normal = q.rotatedVector(QVector3D(QVector2D(cylindervertices.at(k)),0));
                v.position.setX(position.x());
                v.position.setY(position.y());
                v.position.setZ(position.z());
                v.normal.setX(normal.x());
                v.normal.setY(normal.y());
                v.normal.setZ(normal.z());
                v.color = color;
                allvertices << v;
            }
        }
    }
//      Arrows
    scale = SCALE * QVector3D(arrowwidth, arrowwidth, FORCESCALEARROWSHEIGHT*arrowlength);
    QVector3D shiftorig = SCALE * QVector3D(0,0,FORCESCALEHEIGHT*forceslength);
    kshift = maxindex;
    for (int i=0 ; i<MAX_FORCES ; i++){
        if (!visibleforces[i])
            continue;
        for (int j=0 ; j < numcenters ; j++){
            for (int k = 0 ; k < coneindices.length() ; k++){
                allindices << coneindices.at(k) + kshift;
                maxindex = std::max(maxindex,allindices.last()+1);
            }
            allindicesoffset.append(allindices.length());
            kshift = maxindex;
            QVector3D fxyz = forcesxyz[i].at(j);
            q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),fxyz);
            QColor fcolor = forcecolors[i];
            color = QVector4D(fcolor.redF(), fcolor.greenF(), fcolor.blueF(), 1.);
            for (int k = 0 ; k < conevertices.length() ; k++){
                position = q.rotatedVector( scale * conevertices.at(k) + shiftorig * fxyz.length()) + QVector3D(forcesorig[i].at(j));
                normal = q.rotatedVector(QVector3D(conenormals.at(k)));
                v.position.setX(position.x());
                v.position.setY(position.y());
                v.position.setZ(position.z());
                v.normal.setX(normal.x());
                v.normal.setY(normal.y());
                v.normal.setZ(normal.z());
                v.color = color;
                allvertices << v;
            }
        }
    }
    emit updatedisplay();
}

int forces::getarrowlength(){
    return arrowlength;
}

int forces::getarrowwidth(){
    return arrowwidth;
}

int forces::getforceslength(){
    return forceslength;
}

int forces::getforcesthickness(){
    return forcesthickness;
}

QColor forces::getcolor(int i){
    return forcecolors[i];
}

bool forces::getvisibleforces(int i){
    return visibleforces[i];
}

QString forces::getforcesfilename(){
    return forcesfilename;
}

QString forces::getname(){
    return name;
}

QString forces::getpath(){
    return path;
}

bool forces::isvisible(){
    return visible;
}

void forces::setarrowlength(int i){
    arrowlength = i;
}

void forces::setarrowwidth(int i){
    arrowwidth = i;
}

void forces::setcolor(int i, QColor a){
    forcecolors[i] = a;
}

void forces::setforceslength(int i){
    forceslength = i;
}

void forces::setforcesthickness(int i){
    forcesthickness = i;
}

void forces::setname(QString a){
    name = a;
}

//  Sets path to file with critical points
void forces::setpath(QString a){
    path = a;
}

//  Sets Project Folder
void forces::set_ProjectFolder(QString a){
    ProjectFolder = a;
}

//  Sets forces show/hide
void forces::setvisible(bool a){
    visible = a;
}

//  Sets type forces show/hide
void forces::setvisibleforces(int i, bool a){
    if (i < 0 || i >= MAX_FORCES)
        return;
    visibleforces[i] = a;
}
