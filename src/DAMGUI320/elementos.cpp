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
//  File:   elementos.cpp
//
//      Last version: September 2018
//
#include "elementos.h"
#include <QtDebug>

Elementos::Elementos()
{
	p=&Propiedades[0];
}

Elementos::~Elementos()
{
}
//	Gets element atomic number from symbol
int Elementos::getZsymbol(char *Symbol)
{
	ElementData *q;
	q=p;
	string c = string(Symbol);
	for (int i=0;i<MaxElem;i++){
		string s = string(q->Symbol);
		if (c==s) return q->AtomicNumber;
		q++;
	}
	return 0;
}
//	Gets element atomic number from name (only Spanish, English and German)
int Elementos::getZname(char *Name)
{
	ElementData *q;
	q=p;
	string c = string(Name);
	string s;
	for (int i=0;i<MaxElem;i++){
		s = string(q->SpanishName);
		if (c==s) return q->AtomicNumber;
		s = string(q->EnglishName);
		if (c==s) return q->AtomicNumber;
		s = string(q->GermanName);
		if (c==s) return q->AtomicNumber;
		q++;
	}
	return 0;
}
//	Gets element symbol from atomic number
string Elementos::getSymbol(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s) return string(q->Symbol);
		q++;
	}
	return 0;
}
//	Gets element bame from atomic number (Spanish, English, German)
string Elementos::getName(int Z, string Lan)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			if (Lan=="sp") 
				return string(q->SpanishName);
			else if (Lan=="en")
				return string(q->EnglishName);
			else if (Lan=="ge")
				return string(q->GermanName);
		}
		q++;
	}
	return 0;
}
//	Gets element mass from atomic number
double Elementos::getMass(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			if (q->AtomicMass==NN)
				return 0;
			else
				return double(q->AtomicMass);
		}
		q++;
	}
	return 0;
}
//	Gets element volume from atomic number
double Elementos::getVolume(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			if (q->AtomicVolume==NN)
				return 0;
			else
				return double(q->AtomicVolume);
		}
		q++;
	}
	return 0;
}
//	Gets element atomic radius from atomic number
double Elementos::getARadius(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			if (q->AtomicRadius==NN)
				return 0;
			else
				return double(q->AtomicRadius);
		}
		q++;
	}
	return 0;
}
//	Gets element covalent radius from atomic number
double Elementos::getCRadius(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			if (q->CovalentRadius==NN)
				return 0;
			else
				return double(q->CovalentRadius);
		}
		q++;
	}
	return 0;
}
//	Gets element red component of color from atomic number
float Elementos::getColorR(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			return float(q->R);
		}
		q++;
	}
	return 0.0;
}
//	Gets element green component of color from atomic number
float Elementos::getColorG(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			return float(q->G);
		}
		q++;
	}
	return 0.0;
}
//	Gets element blue component of color from atomic number
float Elementos::getColorB(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			return float(q->B);
		}
		q++;
	}
	return 0.0;
}
//	Gets element alternative red component of color from atomic number
float Elementos::getColorR1(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			return float(q->R1);
		}
		q++;
	}
	return 0.0;
}
//	Gets element alternative green component of color from atomic number
float Elementos::getColorG1(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0;i<MaxElem;i++){
		int s =q->AtomicNumber;
		if (Z==s){
			return float(q->G1);
		}
		q++;
	}
	return 0.0;
}
//	Gets element alternative blue component of color from atomic number
float Elementos::getColorB1(int Z)
{
	ElementData *q;
	q=p;
	for (int i=0 ; i<MaxElem ; i++){
		int s =q->AtomicNumber;
		if (Z==s){
			return float(q->B1);
		}
		q++;
	}
	return 0.0;
}
