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
//  Class for retrieving chemical data for a given chemical element. 
//  Available data: atomic number (Z); chemical symbol; Spanish, English and 
//  German name; atomic mass, atomic volume; atomic radius; covalent radius;
//  color (RGB); alternative color (R1 G1 B1).
//  Every function receives an argument, either the symbol, name, Z, etc.
//  and returns the data from table.
//
//  File:   Elements.cpp
//
//      Last version: January 2018
//
#include "elements.h"
#include <QtDebug>

Elements::Elements()
{

}

Elements::~Elements()
{

}

//     Gets element atomic number from symbol
int Elements::getZsymbol(QString Symbol)
{
    for (int i=0 ; i < MaxElem ; i++){
        if (Symbol == QString(Properties[i].Symbol)) return Properties[i].AtomicNumber;
     }
     return 0;
}

//     Gets element atomic number from name (only Spanish, English and German)
int Elements::getZname(QString Name)
{
    for (int i=0 ; i < MaxElem ; i++){
        if (Name == QString(Properties[i].SpanishName) || Name == QString(Properties[i].EnglishName) || Name == QString(Properties[i].GermanName))
          return Properties[i].AtomicNumber;
    }
    return 0;
}

//     Gets element symbol from atomic number
QString Elements::getSymbol(int Z)
{
    if (Z >= 0 && Z <= MaxElem) return Properties[Z].Symbol;
    else return QString("");
}

//     Gets element bame from atomic number (Spanish, English, German)
QString Elements::getName(int Z, QString lan)
{
    if (Z >= 0 && Z <= MaxElem){
        if (lan == QString("sp")) return QString(Properties[Z].SpanishName);
        else if (lan == QString("en")) return QString(Properties[Z].EnglishName);
        else if (lan == QString("ge")) return QString(Properties[Z].GermanName);
        else return QString("");
    }
    return QString("");
}

//     Gets element mass from atomic number
double Elements::getMass(int Z)
{
    if (Z >= 0 && Z <= MaxElem) return Properties[Z].AtomicMass;
    else return 0.;
}

//     Gets element volume from atomic number
double Elements::getVolume(int Z)
{
    if (Z >= 0 && Z <= MaxElem) return Properties[Z].AtomicVolume;
    else return 0.;
}

//     Gets element atomic radius from atomic number
double Elements::getAtomicRadius(int Z)
{
    if (Z >= 0 && Z <= MaxElem) return Properties[Z].AtomicRadius;
    else return 0.;
}

//     Gets element covalent radius from atomic number
double Elements::getCovalentRadius(int Z)
{
    if (Z >= 0 && Z <= MaxElem) return Properties[Z].CovalentRadius;
    else return 0.;
}

QColor Elements::getrasmolColor(int Z){
    return QColor(Properties[Z].R,Properties[Z].G,Properties[Z].B);
}

QColor Elements::getjmolColor(int Z){
    return QColor(Properties[Z].R1,Properties[Z].G1,Properties[Z].B1);
}

