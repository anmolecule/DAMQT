//  Copyright 2008-2016, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//  File:   glmath.h
//
//      Last version: March 2016
//
#ifndef GLMATH_H
#define GLMATH_H

#include <math.h>

const float m_PI=4.0f*atan(1.0f);
const float TOradian = m_PI / 180.0f;
const float TOdegree = 1.0f/TOradian;
const float epsilon = 0.00001f;

typedef float POINT3D[3];
typedef float VECTOR3D[3];
typedef float NEIGHBOR[3][3][3];

class PUNTO
{
public:
	PUNTO(float sx = 0.0f, float sy = 0.0f, float sz = 0.0f);
    ~PUNTO();

    void Reset();
    void Set(float sx, float sy, float sz);
	float x;
	float y;
	float z;

	friend PUNTO operator+(const PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2);
	friend PUNTO operator-(const PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2);
	friend PUNTO operator*(const PUNTO& pt3dPoint, float fScale);
	friend PUNTO operator*(float fScale, const PUNTO& pt3dPoint);
	friend PUNTO operator/(const PUNTO& pt3dPoint, float fScale);
	friend PUNTO& operator*=(PUNTO& pt3dPoint, float fScale);
	friend PUNTO& operator/=(PUNTO& pt3dPoint, float fScale);
	friend PUNTO& operator+=(PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2);
	friend PUNTO& operator-=(PUNTO& pt3dPoint1, const PUNTO& pt3dPoint2);
};

class VECTOR
{
public:
	VECTOR(float sx = 0.0f, float sy = 0.0f, float sz = 0.0f);
    ~VECTOR();

    void Reset();
    void Set(float sx, float sy, float sz);
    float Modulo();
	VECTOR Normalize();
	VECTOR CreateVector(PUNTO pt1,PUNTO pt2);
	VECTOR operator+ (VECTOR v);
	VECTOR operator- (VECTOR v);
	VECTOR operator* (float r);
	VECTOR operator^ (VECTOR v);
	float operator* (VECTOR v);	//dot product
	
    float x;
    float y;
    float z;
};

class QUAT
{
 public:
	QUAT(float sx = 0.0, float sy = 0.0, float sz = 0.0, float sw = 1.0);
	~QUAT();

	void Reset();
	void CreateMatrix(float *pMatrix);
	void CopiaQuat(QUAT q);
	void Set(float sx, float sy, float sz, float sw) {x = sx, y = sy, z = sz, w = sw;}
	void AxisAngleToQuat(VECTOR axis, float theta);
	void EulerToQuat(float pitch, float yaw, float roll);
	void NormalizaQuat();
	float ModuloQuat();
	QUAT ConjugaQuat();
	void Slerp(const QUAT &q1, const QUAT &q2, float t); //Interpolacion
	void Rotatef(float angle, float xAxis, float yAxis, float zAxis);
	
	QUAT operator + (QUAT &q);
	//QUAT& operator = (QUAT &q);
	QUAT* operator += (QUAT &q);
	QUAT operator- (QUAT &q);
	QUAT* operator -= (QUAT &q);
	float operator * (QUAT &q);
	QUAT operator * (float Scalar);
	QUAT operator ^(QUAT &q);
	bool operator == (QUAT &q);
	bool operator != (QUAT &q);
	

	float x;
	float y;
	float z;
	float w;
};

class MATRIX3
{
public:
	MATRIX3();
	MATRIX3(float m11,float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33);
    ~MATRIX3();

	void LoadNull();
	void LoadIdentity();
	MATRIX3 Rotx(float a);
	MATRIX3 Roty(float a);
	MATRIX3 Rotz(float a);
	void CopiaMatrix(float m[9]);
	void MultiplicaMatrix(float m[9]);
	MATRIX3 operator +(MATRIX3 m);
	MATRIX3 operator -(MATRIX3 m);
	MATRIX3 operator *(float m);
	MATRIX3 operator *(MATRIX3 m);
	MATRIX3 Transpose(MATRIX3 m);
	VECTOR operator *(VECTOR m);

	float Element[9];
};

class MATRIX4
{
public:
	MATRIX4();
    ~MATRIX4();

	void LoadIdentity();
	void CopiaMatrix(float m[16]);
	void MultiplicaMatrix(float m[16]);
	MATRIX4 operator *(MATRIX4 m);
	void MatrixInversa(); 
	void MatrixFromAxisAngle(VECTOR axis, float theta);
	void QuatToMatrix(QUAT quat); 
	QUAT operator * (QUAT &q);

	float Element[16];
};

#endif
