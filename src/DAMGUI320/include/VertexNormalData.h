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
//	Header file of translations class
//  Description: translations defines a widget that manages all the elements required for
//  translations of lab frame
//
//	File:   VertexNormalData.h
//
//	Author: Rafael Lopez    (rafael.lopez@uam.es)
//
//	Last version: May 2021
//
#ifndef VERTEXNORMALDATA
#define VERTEXNORMALDATA

#include <QVector3D>
#include <QVector4D>
struct VertexNormalData
{
    QVector3D position;
    QVector3D normal;
    QVector4D color;
};
#endif // VERTEXNORMALDATA
