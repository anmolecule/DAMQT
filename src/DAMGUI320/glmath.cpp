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
//	Classes for facilitating work with several mathematical objects:
//	enables the definition of the objects and verloads operators.
//	Objects:
//		PUNTO = 3D point
//		VECTOR = 3D vector 3D
//		QUAT = Quaternions
//		MATRIX3 = 3x3 matrices
//		MATRIX4 = 4x4 matrices
//
//  File:   glmath.cpp
//
//      Last version: September 2018
//
#include "glmath.h"


//***************
//Class PUNTO
//***************
PUNTO::PUNTO(float sx, float sy, float sz): x(sx), y(sy), z(sz)
{
}

PUNTO::~PUNTO()
{
}

void PUNTO::Reset()
{
    x = 0.0f; y = 0.0; z = 0.0f;
}

void PUNTO::Set(float sx, float sy, float sz)
{
	x = sx, y = sy, z = sz;
}

PUNTO operator+(const PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2)
{
	PUNTO result;

	result.x = pt3dPoint1.x + pt3dPoint2.x;
	result.y = pt3dPoint1.y + pt3dPoint2.y;
	result.z = pt3dPoint1.z + pt3dPoint2.z;

	return result;
}

PUNTO operator-(const PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2)
{
	PUNTO result;

	result.x = pt3dPoint1.x - pt3dPoint2.x;
	result.y = pt3dPoint1.y - pt3dPoint2.y;
	result.z = pt3dPoint1.z - pt3dPoint2.z;

	return result;
}

PUNTO operator*(const PUNTO& pt3dPoint, float fScale)
{
	PUNTO result;

	result.x = pt3dPoint.x*fScale;
	result.y = pt3dPoint.y*fScale;
	result.z = pt3dPoint.z*fScale;

	return result;
}

PUNTO operator*(float fScale, const PUNTO& pt3dPoint)
{
	PUNTO result;

	result.x = pt3dPoint.x*fScale;
	result.y = pt3dPoint.y*fScale;
	result.z = pt3dPoint.z*fScale;

	return result;
}

PUNTO operator/(const PUNTO& pt3dPoint, float fScale)
{
	PUNTO result;

	result.x = pt3dPoint.x/fScale;
	result.y = pt3dPoint.y/fScale;
	result.z = pt3dPoint.z/fScale;
	
	return result;
}

PUNTO& operator*=(PUNTO& pt3dPoint, float fScale)
{
	pt3dPoint.x *= fScale;
	pt3dPoint.y *= fScale;
	pt3dPoint.z *= fScale;

	return pt3dPoint;
}

PUNTO& operator/=(PUNTO& pt3dPoint, float fScale)
{
	pt3dPoint.x /= fScale;
	pt3dPoint.y /= fScale;
	pt3dPoint.z /= fScale;

	return pt3dPoint;
}

PUNTO& operator+=(PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2)
{
	pt3dPoint1.x += pt3dPoint2.x;
	pt3dPoint1.y += pt3dPoint2.y;
	pt3dPoint1.z += pt3dPoint2.z;

	return pt3dPoint1;
}

PUNTO& operator-=(PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2)
{
	pt3dPoint1.x -= pt3dPoint2.x;
	pt3dPoint1.y -= pt3dPoint2.y;
	pt3dPoint1.z -= pt3dPoint2.z;
	
	return pt3dPoint1;
}


//***************
//Class VECTOR
//***************
VECTOR::VECTOR(float sx, float sy, float sz): x(sx), y(sy), z(sz)
{
}

VECTOR::~VECTOR()
{
}

void VECTOR::Reset()
{
    x = 0.0f; y = 0.0; z = 0.0f;
}

void VECTOR::Set(float sx, float sy, float sz)
{
	x = sx; y = sy; z = sz;
}

float VECTOR::Modulo()
{
	return (float)(sqrt(x*x+y*y+z*z));
}

VECTOR VECTOR::Normalize()
{
	VECTOR res;
	res.x=x; res.y=y; res.z=z;
	float l = res.Modulo();
	if (l == 0.0f){
		res.x=0.0f; res.y=0.0f; res.z=0.0f;
	}else{
		res.x = x / l; res.y = y / l; res.z = z / l;}
	return res;
}

VECTOR VECTOR::CreateVector(PUNTO pt1,PUNTO pt2)
{
	VECTOR temp;
	temp.x=pt2.x-pt1.x;
	temp.y=pt2.y-pt1.y;
	temp.z=pt2.z-pt1.z;
	return (temp);
}

VECTOR VECTOR::operator+ (VECTOR v)
{
	VECTOR temp;
	temp.x = x + v.x;
	temp.y = y + v.y;
	temp.z = z + v.z;
	return (temp);
}
VECTOR VECTOR::operator- (VECTOR v)
{
	VECTOR temp;
	temp.x = x - v.x;
	temp.y = y - v.y;
	temp.z = z - v.z;
	return (temp);
}


VECTOR VECTOR::operator* (float r)
{
	VECTOR temp;
	temp.x = r * x;
	temp.y = r * y;
	temp.z = r * z;
	return (temp);
}

VECTOR VECTOR::operator^ (VECTOR v)
{
	VECTOR res;
	res.x = y*v.z - z*v.y;
	res.y = z*v.x - x*v.z;
	res.z = x*v.y - y*v.x;
	return res;
}

float VECTOR::operator* (VECTOR v)	//dot product
{
	return v.x*x+v.y*y+v.z*z;
}


//***************
//Class QUAT (Quaternions)
//***************
QUAT::QUAT(float sx, float sy, float sz, float sw): x(sx), y(sy), z(sz), w(sw)
{
}

QUAT::~QUAT()
{
}

void QUAT::Reset()
{
    x = 0.0f; y = 0.0;  z = 0.0;  w = 1.0;
}

void QUAT::CopiaQuat(QUAT q)
{
    x = q.x;  y = q.y;  z = q.z;  w = q.w;
}

void QUAT::AxisAngleToQuat(VECTOR axis, float theta)
{
        float halfTheta = theta * 0.5f;
        float cosHalfTheta = cos(halfTheta);
        float sinHalfTheta = sin(halfTheta);
        x = axis.x * sinHalfTheta;
        y = axis.y * sinHalfTheta;
        z = axis.z * sinHalfTheta;
        w = cosHalfTheta;
}

void QUAT::EulerToQuat(float roll, float pitch, float yaw)
{
        float cr, cp, cy, sr, sp, sy, cpcy, spsy;
        cr = cos(roll/2);
        cp = cos(pitch/2);
        cy = cos(yaw/2);
        sr = sin(roll/2);
        sp = sin(pitch/2);
        sy = sin(yaw/2);
        cpcy = cp * cy;
        spsy = sp * sy;
        w = cr * cpcy + sr * spsy;
        x = sr * cpcy - cr * spsy;
        y = cr * sp * cy + sr * cp * sy;
        z = cr * cp * sy - sr * sp * cy;
}

float QUAT::ModuloQuat()
{
      return( sqrt(w*w+x*x+y*y+z*z));
}

void QUAT::NormalizaQuat()
{
      float Mag;
      Mag = ModuloQuat();
      w = w/Mag;
      x = x/Mag;
      y = y/Mag;
      z = z/Mag;
}

QUAT QUAT::operator ^(QUAT &q)
{
	QUAT r;

	r.w = w*q.w - x*q.x - y*q.y - z*q.z;
	r.x = w*q.x + x*q.w + y*q.z - z*q.y;
	r.y = w*q.y + y*q.w + z*q.x - x*q.z;
	r.z = w*q.z + z*q.w + x*q.y - y*q.x;
	
	return(r);
}

float QUAT::operator * (QUAT &q) 
{
	return x * q.x + y * q.y + z * q.z + w * q.w;
}

QUAT QUAT::operator * (float Scalar) 
{
	QUAT r;
	r.x = x * Scalar;
	r.y = y * Scalar;
	r.z = z * Scalar;
	r.w = w * Scalar;
	
	return r;
}

QUAT QUAT::operator + (QUAT &q) 
{
	QUAT r(x + q.x, y + q.y, z + q.z, w + q.w);
	return r;
}

QUAT* QUAT::operator += (QUAT &q)
{
	*this= *this + q;
	return this;
}

QUAT QUAT::operator- (QUAT &q) 
{
	QUAT r(x - q.x, y - q.y, z - q.z, w - q.w);
	return r;
}

QUAT* QUAT::operator -= (QUAT &q)
{
	*this= *this - q;
	return this;
}

bool QUAT::operator == (QUAT &q) 
{
	bool br;
	br= fabs(x - q.x) < epsilon;
	br= (fabs(y - q.y) < epsilon) && br;
	br= (fabs(z - q.z) < epsilon) && br;
	br= (fabs(w - q.w) < epsilon) && br;

	return br;
}

bool QUAT::operator != (QUAT &q) 
{
	return !(*this == q);
}


QUAT QUAT::ConjugaQuat()
{
	return QUAT(-x, -y, -z, w);
}

void QUAT::CreateMatrix(float *pMatrix)
{
	// Make sure the matrix has allocated memory to store the rotation data
	if(!pMatrix) return;

	// First row
	pMatrix[ 0] = 1.0f - 2.0f * ( y * y + z * z ); 
	pMatrix[ 1] = 2.0f * (x * y + z * w);
	pMatrix[ 2] = 2.0f * (x * z - y * w);
	pMatrix[ 3] = 0.0f;  

	// Second row
	pMatrix[ 4] = 2.0f * ( x * y - z * w );  
	pMatrix[ 5] = 1.0f - 2.0f * ( x * x + z * z ); 
	pMatrix[ 6] = 2.0f * (z * y + x * w );  
	pMatrix[ 7] = 0.0f;  

	// Third row
	pMatrix[ 8] = 2.0f * ( x * z + y * w );
	pMatrix[ 9] = 2.0f * ( y * z - x * w );
	pMatrix[10] = 1.0f - 2.0f * ( x * x + y * y );  
	pMatrix[11] = 0.0f;  

	// Fourth row
	pMatrix[12] = 0;  
	pMatrix[13] = 0;  
	pMatrix[14] = 0;  
	pMatrix[15] = 1.0f;

	// Now pMatrix[] is a 4x4 homogeneous matrix that can be applied to an OpenGL Matrix
}

//INTERPOLATION
void QUAT::Slerp(const QUAT &q1, const QUAT &q2, float t)
{
	 float cosTheta = 0.0f;
	 float sinTheta = 0.0f;
	 float beta = 0.0f;
	 float q2Array[4];

	 // Temporary array to hold second quaternion.
	 q2Array[0] = q2.x; q2Array[1] = q2.y; q2Array[2] = q2.z; q2Array[3] = q2.w;

	 cosTheta = q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;

	 if(cosTheta < 0.0f)
		{
		   // Flip sigh if so.
		   q2Array[0] = -q2Array[0]; q2Array[1] = -q2Array[1];
		   q2Array[2] = -q2Array[2]; q2Array[3] = -q2Array[3];
		   cosTheta = -cosTheta;
		}

	 beta = 1.0f - t;

	 if(1.0f - cosTheta > 0.001f)
		{
		   // We are using spherical interpolation.
		   cosTheta = (float)acos(cosTheta);
		   sinTheta = 1.0f / (float)sin(cosTheta);
		   beta = (float)sin(cosTheta * beta) * sinTheta;
		   t = (float)sin(cosTheta * t) * sinTheta;
		}

	 // Interpolation.
	 x = beta * q1.x + t * q2Array[0];
	 y = beta * q1.y + t * q2Array[1];
	 z = beta * q1.z + t * q2Array[2];
 w = beta * q1.w + t * q2Array[3];
}

//rotation respect an axis
void QUAT::Rotatef(float angle, float xAxis, float yAxis, float zAxis)
{
	VECTOR v;
	v.Set(xAxis,yAxis,zAxis);
	v.Normalize(); //Normalizamos
	
	angle = angle * TOradian; //pasamos de grado a radianes

	float sine = (float)sin(angle / 2.0f);

	// Creates a quaternion
	x = xAxis * sine;
	y = yAxis * sine;
	z = zAxis * sine;
	w = (float)cos(angle / 2.0f);
}

//***************
//Class MATRIX 3x3
//***************
MATRIX3::MATRIX3()
{
    LoadIdentity();
}

MATRIX3::MATRIX3(float m11,float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33)
{
    Element[0]=m11;
	Element[1]=m12;
	Element[2]=m13;

	Element[3]=m21;
	Element[4]=m22;
	Element[5]=m23;

	Element[6]=m31;
	Element[7]=m32;
	Element[8]=m33;
}

MATRIX3::~MATRIX3()
{
}

void MATRIX3::LoadNull()
{
	for(int i=0;i<9;i++)
		Element[i]=0.0f;
}

void MATRIX3::LoadIdentity()
{
	LoadNull();
	Element[0]=1.0f;
	Element[4]=1.0f;
	Element[8]=1.0f;
}

MATRIX3 MATRIX3::Rotx(float a)
{
	MATRIX3 r;
	
	r.Element[0]=1.0f;
	r.Element[1]=0.0f;
	r.Element[2]=0.0f;

	r.Element[3]=0.0f;
	r.Element[4]=cos(a);
	r.Element[5]=-sin(a);

	r.Element[6]=0.0f;
	r.Element[7]=sin(a);
	r.Element[8]=cos(a);

	return r;
}

MATRIX3 MATRIX3::Roty(float a)
{
	MATRIX3 r;
	
	r.Element[0]=cos(a);
	r.Element[1]=0.0f;
	r.Element[2]=sin(a);

	r.Element[3]=0.0f;
	r.Element[4]=1.0f;
	r.Element[5]=0.0f;

	r.Element[6]=-sin(a);
	r.Element[7]=0.0f;
	r.Element[8]=cos(a);
	
	return r;
}

MATRIX3 MATRIX3::Rotz(float a)
{
	MATRIX3 r;
	
	r.Element[0]=cos(a);
	r.Element[1]=-sin(a);
	r.Element[2]=0.0f;

	r.Element[3]=sin(a);
	r.Element[4]=cos(a);
	r.Element[5]=0.0f;

	r.Element[6]=0.0f;
	r.Element[7]=0.0f;
	r.Element[8]=1.0f;
	
	return r;
}

void MATRIX3::CopiaMatrix(float m[9])
{
	for(int i=0;i<9;i++)
		Element[i]=m[i];
}

MATRIX3 MATRIX3::operator +(MATRIX3 m)
{
	MATRIX3 temp,r;
	temp.CopiaMatrix(this->Element);

	for (int i=0; i<9; i++)
		r.Element[i]=temp.Element[i]+m.Element[i];
	return r;
}

MATRIX3 MATRIX3::operator -(MATRIX3 m)
{
	MATRIX3 temp,r;
	temp.CopiaMatrix(this->Element);

	for (int i=0; i<9; i++)
		r.Element[i]=temp.Element[i]-m.Element[i];
	return r;
}

MATRIX3 MATRIX3::Transpose(MATRIX3 m)
{
	MATRIX3 temp,r;
	temp.CopiaMatrix(this->Element);

	r.Element[0]=m.Element[0];
	r.Element[1]=m.Element[3];
	r.Element[2]=m.Element[6];
	r.Element[3]=m.Element[1];
	r.Element[4]=m.Element[4];
	r.Element[5]=m.Element[7];
	r.Element[6]=m.Element[2];
	r.Element[7]=m.Element[5];
	r.Element[8]=m.Element[8];

	return r;
}

MATRIX3 MATRIX3::operator *(float m)
{
	MATRIX3 temp,r;
	temp.CopiaMatrix(this->Element);

	for (int i=0; i<9; i++)
		r.Element[i]=temp.Element[i]*m;
	return r;
}

MATRIX3 MATRIX3::operator *(MATRIX3 m)
{
	MATRIX3 temp,r;
	temp.CopiaMatrix(this->Element);
	
	r.Element[0]=temp.Element[0]*m.Element[0]+temp.Element[1]*m.Element[3]+temp.Element[2]*m.Element[6];
	r.Element[1]=temp.Element[0]*m.Element[1]+temp.Element[1]*m.Element[4]+temp.Element[2]*m.Element[7];
	r.Element[2]=temp.Element[0]*m.Element[2]+temp.Element[1]*m.Element[5]+temp.Element[2]*m.Element[8];
	r.Element[3]=temp.Element[3]*m.Element[0]+temp.Element[4]*m.Element[3]+temp.Element[5]*m.Element[6];
	r.Element[4]=temp.Element[3]*m.Element[1]+temp.Element[4]*m.Element[4]+temp.Element[5]*m.Element[7];
	r.Element[5]=temp.Element[3]*m.Element[2]+temp.Element[4]*m.Element[5]+temp.Element[5]*m.Element[8];
	r.Element[6]=temp.Element[6]*m.Element[0]+temp.Element[7]*m.Element[3]+temp.Element[8]*m.Element[6];
	r.Element[7]=temp.Element[6]*m.Element[1]+temp.Element[7]*m.Element[4]+temp.Element[8]*m.Element[7];
	r.Element[8]=temp.Element[6]*m.Element[2]+temp.Element[7]*m.Element[5]+temp.Element[8]*m.Element[8];

	return r;
}


void MATRIX3::MultiplicaMatrix(float m[])
{
	MATRIX3 temp;

	temp.CopiaMatrix(this->Element);

	Element[0]=temp.Element[0]*m[0]+temp.Element[1]*m[3]+temp.Element[2]*m[6];
	Element[1]=temp.Element[0]*m[1]+temp.Element[1]*m[4]+temp.Element[2]*m[7];
	Element[2]=temp.Element[0]*m[2]+temp.Element[1]*m[5]+temp.Element[2]*m[8];
	Element[3]=temp.Element[3]*m[0]+temp.Element[4]*m[3]+temp.Element[5]*m[6];
	Element[4]=temp.Element[3]*m[1]+temp.Element[4]*m[4]+temp.Element[5]*m[7];
	Element[5]=temp.Element[3]*m[2]+temp.Element[4]*m[5]+temp.Element[5]*m[8];
	Element[6]=temp.Element[6]*m[0]+temp.Element[7]*m[3]+temp.Element[8]*m[6];
	Element[7]=temp.Element[6]*m[1]+temp.Element[7]*m[4]+temp.Element[8]*m[7];
	Element[8]=temp.Element[6]*m[2]+temp.Element[7]*m[5]+temp.Element[8]*m[8];
}

VECTOR MATRIX3::operator *(VECTOR m)
{
	VECTOR v;
	MATRIX3 temp;
	temp.CopiaMatrix(this->Element);

	v.x=temp.Element[0]*m.x+temp.Element[1]*m.y+temp.Element[2]*m.z;
	v.y=temp.Element[3]*m.x+temp.Element[4]*m.y+temp.Element[5]*m.z;
	v.z=temp.Element[6]*m.x+temp.Element[7]*m.y+temp.Element[8]*m.z;
	
	return v;
}

//***************
//ClasS MATRIX 4x4
//***************
MATRIX4::MATRIX4()
{
    LoadIdentity();
}

MATRIX4::~MATRIX4()
{
}

void MATRIX4::LoadIdentity()
{
	Element[0]=1.0f;
	Element[1]=0.0f;
	Element[2]=0.0f;
	Element[3]=0.0f;

	Element[4]=0.0f;
	Element[5]=1.0f;
	Element[6]=0.0f;
	Element[7]=0.0f;

	Element[8]=0.0f;
	Element[9]=0.0f;
	Element[10]=1.0f;
	Element[11]=0.0f;

	Element[12]=0.0f;
	Element[13]=0.0f;
	Element[14]=0.0f;
	Element[15]=1.0f;
}

void MATRIX4::CopiaMatrix(float m[16])
{
	Element[0 ] = m[0 ];
	Element[1 ] = m[1 ];
	Element[2 ] = m[2 ];
	Element[3 ] = m[3 ];
	Element[4 ] = m[4 ];
	Element[5 ] = m[5 ];
	Element[6 ] = m[6 ];
	Element[7 ] = m[7 ];
	Element[8 ] = m[8 ];
	Element[9 ] = m[9 ];
	Element[10] = m[10];
	Element[11] = m[11];
	Element[12] = m[12];
	Element[13] = m[13];
	Element[14] = m[14];
	Element[15] = m[15];
}

MATRIX4 MATRIX4::operator *(MATRIX4 m)
{
	MATRIX4 temp,r;
	temp.CopiaMatrix(this->Element);
	
	  r.Element[0] = temp.Element[0 ]*m.Element[0 ] + temp.Element[4 ]*m.Element[1 ]
            + temp.Element[8 ]*m.Element[2 ] + temp.Element[12]*m.Element[3 ];

      r.Element[1] = temp.Element[1 ]*m.Element[0 ] + temp.Element[5 ]*m.Element[1 ]
            + temp.Element[9 ]*m.Element[2 ] + temp.Element[13]*m.Element[3 ];

      r.Element[2] = temp.Element[2 ]*m.Element[0 ] + temp.Element[6 ]*m.Element[1 ]
               + temp.Element[10]*m.Element[2 ] + temp.Element[14]*m.Element[3 ];

      r.Element[3] = temp.Element[3 ]*m.Element[0 ] + temp.Element[7 ]*m.Element[1 ]
               + temp.Element[11]*m.Element[2 ] + temp.Element[15]*m.Element[3 ];

      r.Element[4] = temp.Element[0 ]*m.Element[4 ] + temp.Element[4 ]*m.Element[5 ]
               + temp.Element[8 ]*m.Element[6 ] + temp.Element[12]*m.Element[7 ];

      r.Element[5] = temp.Element[1 ]*m.Element[4 ] + temp.Element[5 ]*m.Element[5 ]
               + temp.Element[9 ]*m.Element[6 ] + temp.Element[13]*m.Element[7 ];

      r.Element[6] = temp.Element[2 ]*m.Element[4 ] + temp.Element[6 ]*m.Element[5 ]
               + temp.Element[10]*m.Element[6 ] + temp.Element[14]*m.Element[7 ];

      r.Element[7] = temp.Element[3 ]*m.Element[4 ] + temp.Element[7 ]*m.Element[5 ]
               + temp.Element[11]*m.Element[6 ] + temp.Element[15]*m.Element[7 ];

      r.Element[8] = temp.Element[0 ]*m.Element[8 ] + temp.Element[4 ]*m.Element[9 ]
               + temp.Element[8 ]*m.Element[10] + temp.Element[12]*m.Element[11];

      r.Element[9] = temp.Element[1 ]*m.Element[8 ] + temp.Element[5 ]*m.Element[9 ]
               + temp.Element[9 ]*m.Element[10] + temp.Element[13]*m.Element[11];

      r.Element[10]= temp.Element[2 ]*m.Element[8 ] + temp.Element[6 ]*m.Element[9 ]
               + temp.Element[10]*m.Element[10] + temp.Element[14]*m.Element[11];

      r.Element[11]= temp.Element[3 ]*m.Element[8 ] + temp.Element[7 ]*m.Element[9 ]
               + temp.Element[11]*m.Element[10] + temp.Element[15]*m.Element[11];

      r.Element[12]= temp.Element[0 ]*m.Element[12] + temp.Element[4 ]*m.Element[13]
               + temp.Element[8 ]*m.Element[14] + temp.Element[12]*m.Element[15];

      r.Element[13]= temp.Element[1 ]*m.Element[12] + temp.Element[5 ]*m.Element[13]
               + temp.Element[9 ]*m.Element[14] + temp.Element[13]*m.Element[15];

      r.Element[14]= temp.Element[2 ]*m.Element[12] + temp.Element[6 ]*m.Element[13]
               + temp.Element[10]*m.Element[14] + temp.Element[14]*m.Element[15];

      r.Element[15]= temp.Element[3 ]*m.Element[12] + temp.Element[7 ]*m.Element[13]
               + temp.Element[11]*m.Element[14] + temp.Element[15]*m.Element[15];
	return r;
}


void MATRIX4::MultiplicaMatrix(float m[])
{
      MATRIX4 temp;

      temp.CopiaMatrix(this->Element);

      Element[0] = temp.Element[0 ]*m[0 ] + temp.Element[4 ]*m[1 ]
            + temp.Element[8 ]*m[2 ] + temp.Element[12]*m[3 ];

      Element[1] = temp.Element[1 ]*m[0 ] + temp.Element[5 ]*m[1 ]
            + temp.Element[9 ]*m[2 ] + temp.Element[13]*m[3 ];

      Element[2] = temp.Element[2 ]*m[0 ] + temp.Element[6 ]*m[1 ]
               + temp.Element[10]*m[2 ] + temp.Element[14]*m[3 ];

      Element[3] = temp.Element[3 ]*m[0 ] + temp.Element[7 ]*m[1 ]
               + temp.Element[11]*m[2 ] + temp.Element[15]*m[3 ];

      Element[4] = temp.Element[0 ]*m[4 ] + temp.Element[4 ]*m[5 ]
               + temp.Element[8 ]*m[6 ] + temp.Element[12]*m[7 ];

      Element[5] = temp.Element[1 ]*m[4 ] + temp.Element[5 ]*m[5 ]
               + temp.Element[9 ]*m[6 ] + temp.Element[13]*m[7 ];

      Element[6] = temp.Element[2 ]*m[4 ] + temp.Element[6 ]*m[5 ]
               + temp.Element[10]*m[6 ] + temp.Element[14]*m[7 ];

      Element[7] = temp.Element[3 ]*m[4 ] + temp.Element[7 ]*m[5 ]
               + temp.Element[11]*m[6 ] + temp.Element[15]*m[7 ];

      Element[8] = temp.Element[0 ]*m[8 ] + temp.Element[4 ]*m[9 ]
               + temp.Element[8 ]*m[10] + temp.Element[12]*m[11];

      Element[9] = temp.Element[1 ]*m[8 ] + temp.Element[5 ]*m[9 ]
               + temp.Element[9 ]*m[10] + temp.Element[13]*m[11];

      Element[10]= temp.Element[2 ]*m[8 ] + temp.Element[6 ]*m[9 ]
               + temp.Element[10]*m[10] + temp.Element[14]*m[11];

      Element[11]= temp.Element[3 ]*m[8 ] + temp.Element[7 ]*m[9 ]
               + temp.Element[11]*m[10] + temp.Element[15]*m[11];

      Element[12]= temp.Element[0 ]*m[12] + temp.Element[4 ]*m[13]
               + temp.Element[8 ]*m[14] + temp.Element[12]*m[15];

      Element[13]= temp.Element[1 ]*m[12] + temp.Element[5 ]*m[13]
               + temp.Element[9 ]*m[14] + temp.Element[13]*m[15];

      Element[14]= temp.Element[2 ]*m[12] + temp.Element[6 ]*m[13]
               + temp.Element[10]*m[14] + temp.Element[14]*m[15];

      Element[15]= temp.Element[3 ]*m[12] + temp.Element[7 ]*m[13]
               + temp.Element[11]*m[14] + temp.Element[15]*m[15];
}

void MATRIX4::MatrixInversa()
 {
      MATRIX4 temp;

      temp.CopiaMatrix(this->Element);

      Element[0 ] = temp.Element[0 ];
      Element[1 ] = temp.Element[4 ];
      Element[2 ] = temp.Element[8 ];

      Element[4 ] = temp.Element[1 ];
      Element[5 ] = temp.Element[5 ];
      Element[6 ] = temp.Element[9 ];

      Element[8 ] = temp.Element[2 ];
      Element[9 ] = temp.Element[6 ];
      Element[10] = temp.Element[10];

      Element[12] *= -1.0f;
      Element[13] *= -1.0f;
      Element[14] *= -1.0f;
}

QUAT MATRIX4::operator * (QUAT &q)
{
	float ValInt;
	int i;
	QUAT r;
	float v[4],vn[4];
	
	v[0]=q.x; v[1]=q.y; v[2]=q.z; v[3]=q.w;

	for (int Index= 0; Index < 4; Index++)
	{
		ValInt= 0;
		for (i= 0; i < 4; i++)
		{
			ValInt+= this->Element[(i * 4) + Index] * v[i];
		}
		vn[Index]= ValInt;
	}
	
	if (vn[3] != 1)
	{
		r.NormalizaQuat();
	}

	r.x=vn[0]; r.y=vn[1]; r.z=vn[2]; r.w=vn[3];
	return r;
}

void MATRIX4::MatrixFromAxisAngle(VECTOR axis, float theta)
{
        QUAT q;
        float halfTheta = theta * 0.5f;
        float cosHalfTheta = cos(halfTheta);
        float sinHalfTheta = sin(halfTheta);
        float xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;
        q.x = axis.x * sinHalfTheta;
        q.y = axis.y * sinHalfTheta;
        q.z = axis.z * sinHalfTheta;
        q.w = cosHalfTheta;
        xs = q.x * 2;  ys = q.y * 2;  zs = q.z * 2;
        wx = q.w * xs; wy = q.w * ys; wz = q.w * zs;
        xx = q.x * xs; xy = q.x * ys; xz = q.x * zs;
        yy = q.y * ys; yz = q.y * zs; zz = q.z * zs;
        Element[0] = 1 - (yy + zz);
        Element[1] = xy - wz;
        Element[2] = xz + wy;
        Element[4] = xy + wz;
        Element[5] = 1 - (xx + zz);
        Element[6] = yz - wx;
        Element[8] = xz - wy;
        Element[9] = yz + wx;
        Element[10] = 1 - (xx + yy);
        Element[12] = Element[13] = Element[14] = Element[3] = Element[7] = Element[11] = 0;
        Element[15] = 1;
}

void MATRIX4::QuatToMatrix(QUAT quat)
{
      float wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
      // calcula coeficientes
      x2 = quat.x + quat.x;
      y2 = quat.y + quat.y;
      z2 = quat.z + quat.z;
      xx = quat.x * x2;
      xy = quat.x * y2;
      xz = quat.x * z2;
      yy = quat.y * y2;
      yz = quat.y * z2;
      zz = quat.z * z2;
      wx = quat.w * x2;
      wy = quat.w * y2;
      wz = quat.w * z2;
      Element[0] = 1.0f - (yy + zz);
      Element[1] = xy - wz;
      Element[2] = xz + wy;
      Element[3] = 0.0f;
      Element[4] = xy + wz;
      Element[5] = 1.0f - (xx + zz);
      Element[6] = yz - wx;
      Element[7] = 0.0f;
      Element[8] = xz - wy;
      Element[9] = yz + wx;
      Element[10] = 1.0f - (xx + yy);
      Element[11] = 0.0f;
      Element[12] = 0.0f;
      Element[13] = 0.0f;
      Element[14] = 0.0f;
      Element[15] = 1.0f;
}
