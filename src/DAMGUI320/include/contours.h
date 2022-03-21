/***************************************************************************
 *   Copyright (C) 2007 by Bjorn Harpe,,,                                  *
 *   bjorn@Ouelong.com                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef CONTOURS_H
#define CONTOURS_H
//#include <QVector>
#include <QPointF>
#include <vector>
#include <algorithm>

#define DIFF 0.0005
#define EQ(_x_,_y_) (((_x_-_y_<DIFF)&&(_y_-_x_<DIFF))?1:0)
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

using namespace std;

struct SVector
{
   double dx,dy;
};

struct SPoint
{
   SPoint(double x,double y){this->x=x;this->y=y;}
   SPoint(){}
   double x,y;
};

struct SPair
{
   SPair(SPoint p1, SPoint p2){this->p1=p1;this->p2=p2;}
   SPair reverse(){return(SPair(p2,p1));}
   SPoint p1,p2;
};

bool operator <(SPoint p1, SPoint p2);
bool operator <(SPair p1,SPair p2);
bool operator ==(SPoint p1,SPoint p2);
bool operator !=(SPoint p1,SPoint p2);
SPoint operator +=(SPoint p, SVector v);

class CContour
{
   public:
      CContour(){contour=NULL;}
      ~CContour();
      int merge(CContour *c);
      int reverse();
      int add_vector(SPoint start,SPoint end);
      int condense(double difference = 0.000000001);
//      int dump();
      bool closed(){return(_start==_end);}
      SPoint start(){return(_start);}
      SPoint end(){return(_end);}
      vector<SVector> *contour;
      SPoint _start,_end;
};

class CRaster
{
   public:
      CRaster(){};
      virtual double value(double, double){return(0);};
      virtual SPoint upper_bound(){return(SPoint(0,0));};
      virtual SPoint lower_bound(){return(SPoint(0,0));};
      virtual ~CRaster(){};
   private:
       
};

class CContourLevel
{
   public:
     CContourLevel(){contour_lines=NULL;raw=NULL;};
//     int dump();
     int merge();
     int consolidate();
     vector<CContour*> *contour_lines;
     vector<SPair> *raw;
     ~CContourLevel();
}; 

class CContourMap
{
   public:
      CContourMap();
      int add_segment(SPair t,int level);
//      int dump();
      int contour(CRaster *r);
      int consolidate();
      CContourLevel* level(int i){return((*contour_level)[i]);}
      double alt(int i){return(levels[i]);}
      ~CContourMap();
      vector<CContourLevel*> *contour_level;
      int n_levels;
      double *levels;      
};

class ToMap:public CRaster
{
   public:
      double value(double x, double y);
      SPoint upper_bound();
      SPoint lower_bound();
      void set_upper_bound(SPoint);
      void set_lower_bound(SPoint);
      QVector<double> *surf2d;
      int ncurves;
   private:
       double nxini;
       double nyini;
       double nxfin;
       double nyfin; 
};
#endif
