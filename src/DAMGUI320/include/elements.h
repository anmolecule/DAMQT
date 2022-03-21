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
//  File:   elements.h
//
//      Last version: January 2018
//
#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <QColor>
#include <QList>
#include <QString>

#include <iostream>
#include <string>
//using namespace std;

enum {MaxElem=119,NameLen=16,NaN=-999};

struct ElementData
{  
     // general data
     int AtomicNumber;          // Z (zero for arbitrary point charge) 
     QString Symbol;          // chemical symbol
     QString SpanishName;     // Spanish name
     QString EnglishName;     // English name
     QString GermanName;     // German name
  
     // atomic data
     double AtomicMass;          // atomic mass [uma]
     double AtomicVolume;           // atomic volume
     double AtomicRadius;          // atomic radius [A]
     double CovalentRadius;          // covalent radius [A]

     // rasmol colors
    int R;
    int G;
    int B;
     // jmol colors
    int R1;
    int G1;
    int B1;
};

static ElementData Properties[MaxElem]=
{
    //------------------------------------------------------------------------------------------------------------------------
    //            Spanish         English        German             mass
    //  Z symbol  name            name           name               [amu]           v.a.   r.a. r.c.(A) R   G   B   R1  G1  B1
    //------------------------------------------------------------------------------------------------------------------------
     {   0, "q" , "Carga "         , "Charge"       , "Charge"          , 0.         , 0.   , 0.5 , 0.3 ,  255,100,  0, 255,100,  0},
     {   1, "H" , "Hidrógeno"      , "Hydrogen"     , "Wasserstoff"     , 1.00794    , 14.4 , 0.79, 0.32,  255,255,255,   0,255,255},
     {   2, "He", "Helio"          , "Helium"       , "Helium"          , 4.002602   , NaN  , 0.49, 0.93,  255,192,203, 217,255,255},
     {   3, "Li", "Litio"          , "Lithium"      , "Lithium"         , 6.941      , 13.10, 2.05, 1.23,  178, 34, 34, 204,128,255},
     {   4, "Be", "Berilio"        , "Beryllium"    , "Beryllium"       , 9.012182   , 5.0  , 1.40, 0.90,  255, 20,147, 194,255,  0},
     {   5, "B" , "Boro"           , "Boron"        , "Bor"             , 10.811     , 4.6  , 1.17, 0.82,    0,255,  0, 255,181,181},
     {   6, "C" , "Carbono"        , "Carbon"       , "Kohlenstoff"     , 12.011     , 4.58 , 0.91, 0.77,  200,200,200, 144,144,144},
     {   7, "N" , "Nitrógeno"      , "Nitrogen"     , "Stickstoff"      , 14.00674   , 17.3 , 0.75, 0.75,  143,143,255,  48, 80,248},
     {   8, "O" , "Oxígeno"        , "Oxygen"       , "Sauerstoff"      , 15.9994    , 14.0 , 0.65, 0.73,  240,  0,  0, 255, 13, 13},
     {   9, "F" , "Fluor"          , "Fluorine"     , "Fluor"           , 18.9984032 , 17.1 , 0.57, 0.72,  218,165, 32, 144,224, 80},
     {  10, "Ne", "Neon"           , "Neon"         , "Neon"            , 20.1797    , 16.7 , 0.51, 0.71,  255, 20,147, 179,227,245},
     {  11, "Na", "Sodio"          , "Sodium"       , "Natrium"         , 22.989768  , 23.7 , 2.23, 1.54,    0,  0,255, 171, 92,242},
     {  12, "Mg", "Magnesio"       , "Magnesium"    , "Magnesium"       , 24.3050    , 13.97, 1.72, 1.36,   34,139, 34, 138,255,  0},
     {  13, "Al", "Aluminio"       , "Aluminum"     , "Aluminium"       , 26.981539  , 10.0 , 1.82, 1.18,  128,128,144, 191,166,166},
     {  14, "Si", "Silicio"        , "Silicon"      , "Silicium"        , 28.0855    , 12.1 , 1.46, 1.11,  218,165, 32, 240,200,160},
     {  15, "P" , "Fósforo"        , "Phosphorus"   , "Phosphor"        , 30.97362   , 17.0 , 1.23, 1.06,  255,165,  0, 255,128,  0},
     {  16, "S" , "Azufre"         , "Sulfur"       , "Schwefel"        , 32.066     , 15.5 , 1.09, 1.02,  255,200, 50, 255,255, 48},
     {  17, "Cl", "Cloro"          , "Chlorine"     , "Chlor"           , 35.4527    , 22.7 , 0.97, 0.99,    0,255,  0,  31,240, 31},
     {  18, "Ar", "Argon"          , "Argon"        , "Argon"           , 39.948     , 28.5 , 0.88, 0.98,  255, 20,147, 128,209,227},
     {  19, "K" , "Potasio"        , "Potassium"    , "Kalium"          , 39.0983    , 45.46, 2.77, 2.03,  255, 20,147, 143, 64,212},
     {  20, "Ca", "Calcio"         , "Calcium"      , "Calcium"         , 40.078     , 29.9 , 2.23, 1.91,  128,128,144,  61,255,  0},
     {  21, "Sc", "Escandio"       , "Scandium"     , "Scandium"        , 44.955910  , 15.0 , 2.09, 1.62,  255, 20,147, 230,230,230},
     {  22, "Ti", "Titanio"        , "Titanium"     , "Titan"           , 47.88      , 10.64, 2.00, 1.45,  128,128,144, 191,194,199},
     {  23, "V" , "Vanadio"        , "Vanadium"     , "Vanadium"        , 50.9415    , 8.78 , 1.92, 1.34,  255, 20,147, 166,166,171},
     {  24, "Cr", "Cromo"          , "Chromium"     , "Chrom"           , 51.9961    , 7.23 , 1.85, 1.18,  128,128,144, 138,153,199},
     {  25, "Mn", "Manganeso"      , "Manganese"    , "Mangan"          , 54.93085   , 1.39 , 1.79, 1.17,  128,128,144, 156,122,199},
     {  26, "Fe", "Hierro"         , "Iron"         , "Eisen"           , 55.847     , 7.1  , 1.72, 1.17,  255,165,  0, 224,102, 51},
     {  27, "Co", "Cobalto"        , "Cobal{t"       , "Cobalt"         , 58.93320   , 6.7  , 1.67, 1.16,  255, 20,147, 240,144,160},
     {  28, "Ni", "Níquel"         , "Nickel"       , "Nickel"          , 58.69      , 6.59 , 1.62, 1.15,  165, 42, 42,  80,208, 80},
     {  29, "Cu", "Cobre"          , "Copper"       , "Kupfer"          , 63.546     , 7.1  , 1.57, 1.17,  165, 42, 42, 200,128, 51},
     {  30, "Zn", "Zinc"           , "Zinc"         , "Zink"            , 65.39      , 9.2  , 1.53, 1.25,  165, 42, 42, 125,128,176},
     {  31, "Ga", "Galio"          , "Gallium"      , "Gallium"         , 69.723     , 11.8 , 1.81, 1.26,  255, 20,147, 194,143,143},
     {  32, "Ge", "Germanio"       , "Germanium"    , "Germanium"       , 72.61      , 13.6 , 1.52, 1.22,  255, 20,147, 102,143,143},
     {  33, "As", "Arsénico"       , "Arsenic"      , "Arsen"           , 74.92159   , 13.1 , 1.33, 1.20,  255, 20,147, 189,128,227},
     {  34, "Se", "Selenio"        , "Selenium"     , "Selen"           , 78.96      , 16.45, 1.22, 1.16,  255, 20,147, 255,161,  0},
     {  35, "Br", "Bromo"          , "Bromine"      , "Brom"            , 79.904     , 23.5 , 1.12, 1.14,  165, 42, 42, 166, 41, 41},
     {  36, "Kr", "Krypton"        , "Krypton"      , "Krypton"         , 83.80      , 38.9 , 1.03, 1.12,  255, 20,147,  92,184,209},
     {  37, "Rb", "Rubidio"        , "Rubidium"     , "Rubidium"        , 85.4678    , 55.9 , 2.98, 2.16,  255, 20,147, 112, 46,176},
     {  38, "Sr", "Estroncio"      , "Strontium"    , "Strontium"       , 87.62      , 33.7 , 2.45, 1.91,  255, 20,147,   0,255,  0},
     {  39, "Y" , "Itrio"          , "Yttrium"      , "Yttrium"         , 88.90585   , 19.8 , 2.27, 1.62,  255, 20,147, 148,255,255},
     {  40, "Zr", "Zirconio"       , "Zirconium"    , "Zirconium"       , 91.224     , 14.1 , 2.16, 1.45,  255, 20,147, 148,224,224},
     {  41, "Nb", "Niobio"         , "Niobium"      , "Niob"            , 92.90638   , 10.87, 2.09, 1.34,  255, 20,147, 115,194,201},
     {  42, "Mo", "Molibdeno"      , "Molybdenum"   , "Molybd�n"        , 95.94      , 9.4  , 2.01, 1.30,  255, 20,147,  84,181,181},
     {  43, "Tc", "Tecnecio"       , "Technetium"   , "Technetium"      , 98.9063    , 8.5  , 1.95, 1.27,  255, 20,147,  59,158,158},
     {  44, "Ru", "Rutenio"        , "Ruthenium"    , "Ruthenium"       , 101.07     , 8.3  , 1.89, 1.25,  255, 20,147,  36,143,143},
     {  45, "Rh", "Rodio"          , "Rhodium"      , "Rhodium"         , 102.90550  , 8.3  , 1.83, 1.25,  255, 20,147,  10,125,140},
     {  46, "Pd", "Paladio"        , "Palladium"    , "Palladium"       , 106.42     , 8.9  , 1.79, 1.28,  255, 20,147,   0,105,133},
     {  47, "Ag", "Plata"          , "Silver"       , "Silber"          , 107.8682   , 10.3 , 1.75, 1.34,  128,128,144, 192,192,192},
     {  48, "Cd", "Cadmio"         , "Cadmium"      , "Cadmium"         , 112.411    , 13.1 , 1.71, 1.48,  255, 20,147, 255,217,143},
     {  49, "In", "Indio"          , "Indium"       , "Indium"          , 114.82     , 15.7 , 2.00, 1.44,  255, 20,147, 166,117,115},
     {  50, "Sn", "Estaño"         , "Tin"          , "Zinn"            , 118.710    , 16.3 , 1.72, 1.41,  165, 42, 42, 102,128,128},
     {  51, "Sb", "Antimonio"      , "Antimony"     , "Antimon"         , 121.75     , 18.23, 1.53, 1.40,  255, 20,147, 158, 99,181},
     {  52, "Te", "Teluro"         , "Tellurium"    , "Tellur"          , 127.60     , 20.5 , 1.42, 1.36,  255, 20,147, 212,122,  0},
     {  53, "I" , "Iodo"           , "Iodine"       , "Iod"             , 126.90447  , 25.74, 1.32, 1.33,  160, 32,240, 148,  0,148},
     {  54, "Xe", "Xenon"          , "Xenon"        , "Xenon"           , 131.29     , 37.3 , 1.24, 1.31,  255,255,255,  66,158,176},
     {  55, "Cs", "Cesio"          , "Cesium"       , "Cesium"          , 132.90543  , 71.07, 3.34, 2.35,  255,255,255,  87, 23,143},
     {  56, "Ba", "Bario"          , "Barium"       , "Barium"          , 137.327    , 39.24, 2.78, 1.98,  255,165,  0,   0,201,  0},
     {  57, "La", "Lantano"        , "Lanthanum"    , "Lanthan"         , 138.9055   , 20.73, 2.74, 1.69,  255, 20,147, 112,212,255},
     {  58, "Ce", "Cerio"          , "Cerium"       , "Cer"             , 140.115    , 20.67, 2.70, 1.65,  255, 20,147, 255,255,199},
     {  59, "Pr", "Praseodimio"    , "Praseodymium" , "Praseodym"       , 140.90765  , 20.8 , 2.67, 1.65,  255, 20,147, 217,255,199},
     {  60, "Nd", "Neodimio"       , "Neodymium"    , "Neodym"          , 144.24     , 20.6 , 2.64, 1.64,  255, 20,147, 199,255,199},
     {  61, "Pm", "Promecio"       , "Promethium"   , "Promethium"      , 146.9151   , 22.39, 2.62, 1.63,  255, 20,147, 163,255,199},
     {  62, "Sm", "Samario"        , "Samarium"     , "Samarium"        , 150.36     , 19.95, 2.59, 1.62,  255, 20,147, 143,255,199},
     {  63, "Eu", "Europio"        , "Europium"     , "Europium"        , 151.965    , 28.9 , 2.56, 1.85,  255, 20,147,  97,255,199},
     {  64, "Gd", "Gadolinio"      , "Gadolinium"   , "Gadolinium"      , 157.25     , 19.9 , 2.54, 1.61,  255, 20,147,  69,255,199},
     {  65, "Tb", "Terbio"         , "Terbium"      , "Terbium"         , 158.92534  , 19.2 , 2.51, 1.59,  255, 20,147,  48,255,199},
     {  66, "Dy", "Disprosio"      , "Dysprosium"   , "Dysprosium"      , 162.50     , 19.0 , 2.49, 1.59,  255, 20,147,  31,255,199},
     {  67, "Ho", "Holmio"         , "Holmium"      , "Holmium"         , 164.93032  , 18.7 , 2.47, 1.58,  255, 20,147,   0,255,156},
     {  68, "Er", "Erbio"          , "Erbium"       , "Erbium"          , 167.26     , 18.4 , 2.45, 1.57,  255, 20,147,   0,230,117},
     {  69, "Tm", "Tulio"          , "Thulium"      , "Thulium"         , 168.93421  , 18.1 , 2.42, 1.56,  255, 20,147,   0,212, 82},
     {  70, "Yb", "Iterbio"        , "Ytterbium"    , "Ytterbium"       , 173.04     , 24.79, 2.40, 1.74,  255, 20,147,   0,191, 56},
     {  71, "Lu", "Lutecio"        , "Lutetium"     , "Lutetium"        , 174.967    , 17.78, 2.25, 1.56,  255, 20,147,   0,171, 36},
     {  72, "Hf", "Hafnio"         , "Hafnium"      , "Hafnium"         , 178.49     , 13.6 , 2.16, 1.44,  255, 20,147,  77,194,255},
     {  73, "Ta", "Tantalio"       , "Tantalum"     , "Tantal"          , 180.9479   , 10.90, 2.09, 1.34,  255, 20,147,  77,166,255},
     {  74, "W" , "Wolframio"      , "Tungsten"     , "Wolfram"         , 183.85     , 9.53 , 2.02, 1.30,  255, 20,147,  33,148,214},
     {  75, "Re", "Renio"          , "Rhenium"      , "Rhenium"         , 186.207    , 8.85 , 1.97, 1.28,  255, 20,147,  38,125,171},
     {  76, "Os", "Osmio"          , "Osmium"       , "Osmium"          , 190.2      , 8.49 , 1.92, 1.26,  255, 20,147,  38,102,150},
     {  77, "Ir", "Iridio"         , "Iridium"      , "Iridium"         , 192.22     , 8.54 , 1.87, 1.27,  255, 20,147,  23, 84,135},
     {  78, "Pt", "Platino"        , "Platinum"     , "Platin"          , 195.08     , 9.10 , 1.83, 1.30,  255, 20,147, 208,208,224},
     {  79, "Au", "Oro"            , "Gold"         , "Gold"            , 196.96654  , 10.2 , 1.79, 1.34,  218,165, 32, 255,209, 35},
     {  80, "Hg", "Mercurio"       , "Mercury"      , "Quecksilber"     , 200.59     , 14.82, 1.76, 1.49,  255, 20,147, 184,184,208},
     {  81, "Tl", "Talio"          , "Thallium"     , "Thallium"        , 204.3833   , 17.2 , 2.08, 1.48,  255, 20,147, 166, 84, 77},
     {  82, "Pb", "Plomo"          , "Lead"         , "Blei"            , 207.2      , 18.17, 1.81, 1.47,  255, 20,147,  87, 89, 97},
     {  83, "Bi", "Bismuto"        , "Bismuth"      , "Bismut"          , 208.98037  , 21.3 , 1.63, 1.46,  255, 20,147, 158, 79,181},
     {  84, "Po", "Polonio"        , "Polonium"     , "Polonium"        , 208.9824   , 22.23, 1.53, 1.46,  255, 20,147, 171, 92,  0},
     {  85, "At", "Astato"         , "Astatine"     , "Astat"           , 209.9871   ,  NaN , 1.43, 1.45,  255, 20,147, 117, 79, 69},
     {  86, "Rn", "Radon"          , "Radon"        , "Radon"           , 222.0176   , 50.5 , 1.34, 1.43,  255, 20,147,  66,130,150},
     {  87, "Fr", "Francio"        , "Francium"     , "Francium"        , 223.0197   ,  NaN , 3.50, 2.50,  255, 20,147,  66,  0,102},
     {  88, "Ra", "Radio"          , "Radium"       , "Radium"          , 226.0254   , 45.20, 3.00, 2.40,  255, 20,147,   0,125,  0},
     {  89, "Ac", "Actinio"        , "Actinium"     , "Actinium"        , 227.0278   , 22.54, 3.20, 2.20,  255, 20,147, 112,171,250},
     {  90, "Th", "Torio"          , "Thorium"      , "Thorium"         , 232.0381   , 19.9 , 3.16, 1.65,  255, 20,147,   0,186,255},
     {  91, "Pa", "Protactinio"    , "Protactinium" , "Protactinium"    , 231.0359   , 15.0 , 3.14,  NaN,  255, 20,147,   0,161,255},
     {  92, "U" , "Uranio"         , "Uranium"      , "Uran"            , 238.0289   , 12.59, 3.11, 1.42,  255, 20,147,   0,143,255},
     {  93, "Np", "Neptunio"       , "Neptunium"    , "Neptunium"       , 237.0482   , 11.62, 3.08,  NaN,  255, 20,147,   0,128,255},
     {  94, "Pu", "Plutonio"       , "Plutonium"    , "Plutonium"       , 244.0642   , 12.32, 3.05,  NaN,  255, 20,147,   0,107,255},
     {  95, "Am", "Americio"       , "Americium"    , "Americium"       , 243.0614   , 17.86, 3.02,  NaN,  255, 20,147,  84, 92,242},
     {  96, "Cm", "Curio"          , "Curium"       , "Curium"          , 247.0703   , 18.28, 2.99,  NaN,  255, 20,147, 120, 92,227},
     {  97, "Bk", "Berkelio"       , "Berkelium"    , "Berkelium"       , 247.0703   ,  NaN , 2.97,  NaN,  255, 20,147, 138, 79,227},
     {  98, "Cf", "Californio"     , "Californium"  , "Californium"     , 251.0796   ,  NaN , 2.95,  NaN,  255, 20,147, 161, 54,212},
     {  99, "Es", "Einstenio"      , "Einsteinium"  , "Einsteinium"     , 252.0829   ,  NaN , 2.92,  NaN,  255, 20,147, 179, 31,212},
     { 100, "Fm", "Fermio"         , "Fermium"      , "Fermium"         , 257.0951   ,  NaN , 2.90,  NaN,  255, 20,147, 179, 31,186},
     { 101, "Md", "Mendelevio"     , "Mendelevium"  , "Mendelevium"     , 258.0986   ,  NaN , 2.87,  NaN,  255, 20,147, 179, 13,166},
     { 102, "No", "Nobelio"        , "Nobelium"     , "Nobelium"        , 259.1009   ,  NaN , 2.85,  NaN,  255, 20,147, 189, 13,135},
     { 103, "Lr", "Lawrencio"      , "Lawrencium"   , "Lawrencium"      , 260.1053   ,  NaN , 2.82,  NaN,  255, 20,147, 199,  0,102},
     { 104, "Rf", "Rutherfordio"   , "Rutherfordium", "Rutherfordium"   , 261.1087   ,  NaN ,  NaN,  NaN,  255, 20,147, 204,  0, 89},
     { 105, "Db", "Dubnio"         , "Dubnium"      , "Dubnium"         , 262.1138   ,  NaN ,  NaN,  NaN,  255, 20,147, 209,  0, 79},
     { 106, "Sg", "Seaborgio"      , "Seaborgium"   , "Seaborgium"      , 263.1182   ,  NaN ,  NaN,  NaN,  255, 20,147, 217,  0, 69},
     { 107, "Bh", "Bohrio"         , "Bohrium"      , "Bohrium"         , 262.1229   ,  NaN ,  NaN,  NaN,  255, 20,147, 224,  0, 56},
     { 108, "Hs", "Hassio"         , "Hassium"      , "Hassium"         , 265.0      ,  NaN ,  NaN,  NaN,  255, 20,147, 230,  0, 46},
     { 109, "Mt", "Meitnerio"      , "Meitnerium"   , "Meitnerium"      , 266.0      ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 110, "Ds", "Darmstadtio"    , "Darmstadtium" , "Darmstadtium"    , 271.0      ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 111, "Uuu", "Unununio"      , "Unnununium"   , "Unnununium"      , 272.0      ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 112, "Uub","Ununbio"        , "Unnunbium"    , "Unnunbium"       , 277.0      ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 113, "Uut","Ununtrio"       , "Unnuntrium"   , "Unnuntrium"      , NaN        ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 114, "Uuq","Ununquadio"     , "Unnunquadium" , "Unnunquadium"    , 289.0      ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 115, "Uup","Ununpentio"     , "Unnunpentium" , "Unnunpentium"    , NaN        ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 116, "Uuh","Ununhexio"      , "Unnunhexium"  , "Unnunhexium"     , NaN        ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 117, "Uus","Ununseptio"     , "Unnunseptium" , "Unnunseptium"    , NaN        ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38},
     { 118, "Uuo","Ununoctio"      , "Unnunoctium"  , "Unnunoctium"     , NaN        ,  NaN ,  NaN,  NaN,  255, 20,147, 235,  0, 38}
};

class Elements
{
public:
     Elements();
     ~Elements();

     int getZsymbol(QString Symbol);
     int getZname(QString Name);

     double getMass(int Z);
     double getVolume(int Z);
     double getAtomicRadius(int Z);
     double getCovalentRadius(int Z);

     QColor getrasmolColor(int Z);
     QColor getjmolColor(int Z);

     QString getSymbol(int Z);
     QString getName(int Z, QString Lan);

private:
     ElementData *p;     
};

#endif
