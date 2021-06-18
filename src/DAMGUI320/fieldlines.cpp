//  Copyright 2008-2019, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//  File:   fieldlines.cpp
//  Description: implements fieldlines class
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//      Last version: March 2019
//
#include "fieldlines.h"
#include <cmath>

fieldlines::fieldlines(QString parentpath, QWidget *parent) : QWidget(parent)
{
    allarrowsindices.clear();
    allarrowsindicesoffset.clear();
    allarrowsvertices.clear();
    allindices.clear();
    allindicesoffset.clear();
    allvertices.clear();
    connections.clear();
    QDLfieldlines = nullpointer;
    SPBlineswidth = nullpointer;
    arrowsseparation = 50;
    arrowssize = 1;
    arrowswidth = 1;
    linescolor = QColor(0,250,0);
    lineswidth = 1.;
    fieldfilename = "";
    name = "";
    ProjectFolder = "";
    showarrows = false;
    visible = true;
    setpath(parentpath);
    makeCircumference();    // computes vertices a polygon to be used in base of arrows
}
//	Destructor
fieldlines::~fieldlines(){
    for (int i = 0 ; i < connections.size() ; i++){
        QObject::disconnect(connections.at(i));
    }
    connections.clear();
    if (QDLfieldlines){
        delete QDLfieldlines;
        QDLfieldlines = nullpointer;
    }
}

//  Reads file with field lines
//
bool fieldlines::readfieldlines(QString filename){
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox msgBox;
        msgBox.setText(tr("readfieldlines"));
        msgBox.setInformativeText(tr("File %1 cannot be read").arg(filename)+QString(":\n%1.").arg(file.errorString()));
        msgBox.setIcon(QMessageBox::Critical);
        msgBox.exec();
        return false;
    }
    fieldfilename = filename;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTextStream in(&file);
    QString line;
    GLuint indices = 0;
    allarrowsindices.clear();
    allarrowsindicesoffset.clear();
    allarrowsvertices.clear();
    allindices.clear();
    allvertices.clear();
    allindicesoffset.clear();
    allindicesoffset.append(0);
    VertexNormalData v;
    bool linestart = true;
    bool previousblank = false;
    while (!in.atEnd()){
        line = in.readLine();
#if QT_VERSION < 0x050E00
        QStringList xy = line.split(' ',QString::SkipEmptyParts);
#else
        QStringList xy = line.split(' ',Qt::SkipEmptyParts);
#endif
        if (xy.count() == 3){
            if (!linestart){
                allindices.append(indices);
            }
            allindices.append(indices++);
            v.position.setX(xy.at(0).toFloat());
            v.position.setY(xy.at(1).toFloat());
            v.position.setZ(xy.at(2).toFloat());
            v.normal.setX(0.5);
            v.normal.setY(0.5);
            v.normal.setZ(0.5);
            v.color = QVector4D(linescolor.redF(), linescolor.greenF(), linescolor.blueF(), 1.);
            allvertices.append(v);
            linestart = false;
            previousblank = false;
        }
        else if(!previousblank){
            allindices.removeLast();
            allindicesoffset.append(allindices.count());
            linestart = true;
            previousblank = true;
        }
    }
    allindicesoffset.append(allindices.count());
    file.close();
    QApplication::restoreOverrideCursor();
    return true;
}

// -----------------------------------------------------------------------------------
//                      Functions and SLOTS
// -----------------------------------------------------------------------------------

bool fieldlines::getshowarrows(){
    return showarrows;
}

bool fieldlines::isvisible(){
    return visible;
}

int fieldlines::getarrowsseparation(){
    return arrowsseparation;
}

int fieldlines::getarrowssize(){
    return arrowssize;
}

int fieldlines::getarrowswidth(){
    return arrowswidth;
}

int fieldlines::getlineswidth(){
    return lineswidth;
}

QString fieldlines::getname(){
    return name;
}

QString fieldlines::getpath(){
    return path;
}

//  Creates arrows on lines
void fieldlines::createarrows(){
    if (allvertices.isEmpty())
        return;
    allarrowsindices.clear();
    allarrowsvertices.clear();
    allarrowsindicesoffset.clear();
    allarrowsindicesoffset.append(0);
    for (int i = 0 ; i < allindicesoffset.count()-1 ; i++){
        if (allindices.length() < allindicesoffset.at(i)+2*(arrowsseparation+arrowssize))
            break;
        int jinf = allindicesoffset.at(i)+2*(arrowsseparation+arrowssize);
        int jsup = allindicesoffset.at(i+1)-(arrowsseparation+arrowssize);
        for (int j = jinf ; j < jsup ; j+=2*arrowsseparation){
            makeArrow(allvertices.at(allindices.at(j-2*arrowssize)),
                  allvertices.at(allindices.at(j+2*arrowssize)));
        }
    }
}

//  Creates an arrow (computes vertices and normals to be used as GL_TRIANGLE_FAN)
void fieldlines::makeArrow(VertexNormalData v1, VertexNormalData v2){
    VertexNormalData vrot;
    QQuaternion q;
    QVector3D v1pos = v1.position;
    QVector3D v2pos = v2.position;
    QVector4D arrowcol = QVector4D(linescolor.redF(), linescolor.greenF(), linescolor.blueF(), 1.);
    double coneheight = (v2pos-v1pos).length();
    double radius = 0.1 * coneheight * arrowswidth;
    QVector3D conevert;
    QVector3D basenormal;

    q = QQuaternion::rotationTo(QVector3D(0.,0.,1.),v2pos-v1pos);
//    Cone vertex
    conevert = q.rotatedVector(QVector3D(0.,0.,1.)); // Rotates vertex
    vrot.normal = conevert;
    vrot.position = coneheight*conevert + v1.position;
    vrot.color = arrowcol;
    int index = allarrowsindices.length();
    allarrowsindices.append(index++);
    allarrowsvertices.append(vrot);
//    Side of cone
    for (int i = 0 ; i < circumference.count() ; i++){
        conevert = q.rotatedVector(circumference.at(i)); // Rotates vertex
        vrot.normal = conevert;
        vrot.position = radius*conevert + v1.position;
        vrot.color = arrowcol;
        allarrowsindices.append(index++);
        allarrowsvertices.append(vrot);
    }
//    End of cone side
//    Base of cone
    conevert = QVector3D(0.,0.,0.);
    basenormal = q.rotatedVector(QVector3D(0.,0.,-1.)); // Rotates base normal (0,0,-1)
    vrot.normal = conevert;
    vrot.position = coneheight*conevert + v1.position;
    vrot.color = arrowcol;
    allarrowsindices.append(index++);
    allarrowsvertices.append(vrot);
    for (int i = 0 ; i < circumference.count() ; i++){
        vrot.normal = basenormal;
        conevert = q.rotatedVector(circumference.at(i)); // Rotates vertex
        vrot.position = radius*conevert + v1.position;
        vrot.color = arrowcol;
        allarrowsindices.append(index++);
        allarrowsvertices.append(vrot);
    }
//    End of cone base
    allarrowsindicesoffset.append(allarrowsindices.length());
}


//  Creates vertices of a circumference on plane XY with center at (0,0,0) and radius 1
void fieldlines::makeCircumference(){
    double PI = 3.14159265358979;
    int numVertices = 32;
    circumference.clear();
    // vertices on circumference
    double deltaLong = 2. * PI / numVertices;
    for (int j = 0; j < numVertices+1 ; j++) {
        circumference.append(QVector3D(cos(deltaLong * j), sin(deltaLong * j),0));
    }
}


QString fieldlines::getfieldfilename(){
    return fieldfilename;
}

QColor fieldlines::getlinescolor(){
    return linescolor;
}

void fieldlines::setarrowsseparation(int i){
    arrowsseparation = i;
}

void fieldlines::setarrowssize(int i){
    arrowssize = i;
}

void fieldlines::setarrowswidth(int i){
    arrowswidth = i;
}

//  Sets lines color
void fieldlines::setlinescolor(QColor a){
    linescolor = a;
    VertexNormalData v;
    for (int i = 0 ; i < allvertices.count() ; i++){
        v = allvertices[i];
        v.color = QVector4D(linescolor.redF(), linescolor.greenF(), linescolor.blueF(), 1.);
        allvertices[i] = v;
    }
    for (int i = 0 ; i < allarrowsvertices.count() ; i++){
        v = allarrowsvertices[i];
        v.color = QVector4D(linescolor.redF(), linescolor.greenF(), linescolor.blueF(), 1.);
        allarrowsvertices[i] = v;
    }
}

void fieldlines::setlineswidth(int i){
    lineswidth = i;
}

void fieldlines::setname(QString a){
    name = a;
}

//  Sets path to file with lines
void fieldlines::setpath(QString a){
    path = a;
}

//  Sets path to file with lines
void fieldlines::set_ProjectFolder(QString a){
    ProjectFolder = a;
}

//  Sets arrows show/hide
void fieldlines::setshowarrows(bool a){
    showarrows = a;
}


//  Sets lines show/hide
void fieldlines::setvisible(bool a){
    visible = a;
}

/*******************************************************************************************************/
/********************************  Class editCpsDialog  implementation  *******************************/
/*******************************************************************************************************/

editFieldDialog::editFieldDialog(QWidget *parent) : QDialog(parent)
{

}

editFieldDialog::~editFieldDialog(){

}

void editFieldDialog::reject(){
    emit closed();
}
void editFieldDialog::closeEvent(QCloseEvent *event){
    event->ignore();
    this->setVisible(false);
    emit closed();
    event->accept();
}
