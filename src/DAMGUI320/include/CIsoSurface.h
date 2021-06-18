//  Copyright 2000, Raghavendra Chandrashekara;
//  Copyright 2008-2018, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
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
//------------------------------------------------------------------------
//
// File Name: CIsoSurface.h
// Last Modified: 05/8/2000
// Author: Raghavendra Chandrashekara (based on source code
// provided by Paul Bourke and Cory Gene Bloyd)
// Email: rc99@doc.ic.ac.uk, rchandrashekara@hotmail.com
//
// Description: This is the header file for the CIsoSurface class.
// CIsoSurface can be used to construct an isosurface from a scalar
// field.
//
#ifndef CISOSURFACE_H
#define CISOSURFACE_H


#include <map>
#include <vector>
#include <QFile>
#include <QTextStream>
#include <QVector>
#include <QVector3D>

#if __cplusplus <= 199711L
    #define nullpointer NULL
#else
//  C++11 compliant compiler
    #define nullpointer nullptr
#endif

typedef float POINT3D[3];
typedef float VECTOR3D[3];

struct POINT3DID {
	unsigned int newID;
	float x, y, z;
	float Nx, Ny, Nz;
};

typedef std::map<unsigned int, POINT3DID> ID2POINT3DID;

struct TRIANGLE {
	unsigned int pointID[3];
};

typedef std::vector<TRIANGLE> TRIANGLEVECTOR;

template <class T> class CIsoSurface {
public:
QFile file;
QTextStream *out;

	// Constructor and destructor.
	CIsoSurface();
	~CIsoSurface();

	// Generates the isosurface from the scalar field contained in the
	// buffer ptScalarField[].
	void GenerateSurface(const T* ptScalarField, T tIsoLevel, unsigned int nCellsX, unsigned int nCellsY,  
			unsigned int nCellsZ, float fCellLengthX, float fCellLengthY, float fCellLengthZ);
	
	// Generates the isosurface from the scalar field contained in the buffer ptScalarField[] and its gradient contained in
	// buffers ptScalarField_dx[], ptScalarField_dy[], ptScalarField_dz[].
	void GenerateSurfacewithgrad(const T* ptScalarField, const T* ptScalarField_dx, const T* ptScalarField_dy, const T* ptScalarField_dz,  
			T tIsoLevel, unsigned int nCellsX, unsigned int nCellsY, unsigned int nCellsZ, float fCellLengthX, float fCellLengthY, float fCellLengthZ);
	
	// Returns true if a valid surface has been generated.
	bool IsSurfaceValid();

	// Deletes the isosurface.
	void DeleteSurface();

	// Returns the length, width, and height of the volume in which the
	// isosurface in enclosed in.  Returns -1 if the surface is not valid.
    int GetVolumeLengths(float& fVolLengthX, float& fVolLengthY, float& fVolLengthZ);

//    The following variables and arrays have been moved from protected to public

    // The number of vertices which make up the isosurface.
    unsigned int m_nVertices;

    // The vertices which make up the isosurface.

    POINT3D* m_ppt3dVertices;

    // The number of triangles which make up the isosurface.
    unsigned int m_nTriangles;

    // The indices of the vertices which make up the triangles.

    unsigned int* m_piTriangleIndices;

    // The number of normals.
    unsigned int m_nNormals;

    // The normals.
    VECTOR3D* m_pvec3dNormals;

protected:

    // List of POINT3Ds which form the isosurface.
    ID2POINT3DID m_i2pt3idVertices;

	// List of TRIANGLES which form the triangulation of the isosurface.
	TRIANGLEVECTOR m_trivecTriangles;

	// Returns the edge ID.
	unsigned int GetEdgeID(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);

	// Returns the vertex ID.
	unsigned int GetVertexID(unsigned int nX, unsigned int nY, unsigned int nZ);

	// Calculates the intersection point of the isosurface with an edge.
	POINT3DID CalculateIntersection(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);
	
	// Calculates the intersection point of the isosurface with an edge and interpolates surface gradient.
	POINT3DID CalculateIntersectionfandgrad(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo);

	// Interpolates between two grid points to produce the point at which the isosurface intersects an edge.
	POINT3DID Interpolate(float fX1, float fY1, float fZ1, float fX2, float fY2, float fZ2, T tVal1, T tVal2);
	
	// Interpolates between two grid points to produce the point at which the isosurface intersects an edge. Interpolate normals too.
	POINT3DID Interpolatefandgrad(float fX1, float fY1, float fZ1, float fX2, float fY2, float fZ2, 
	   float fNX1, float fNY1, float fNZ1, float fNX2, float fNY2, float fNZ2, 
	   T tVal1, T tVal2);
 
	// Renames vertices and triangles so that they can be accessed more efficiently.
	void RenameVerticesAndTriangles(float v);

	// Calculates the normals.
	void CalculateNormals(float v);

	// No. of cells in x, y, and z directions.
	unsigned int m_nCellsX, m_nCellsY, m_nCellsZ;

	// Cell length in x, y, and z directions.
	float m_fCellLengthX, m_fCellLengthY, m_fCellLengthZ;

    // The buffer holding the scalar field.
	const T* m_ptScalarField;
	// The buffers holding the gradient of the scalar field
	const T* m_ptScalarField_dx;
	const T* m_ptScalarField_dy;
	const T* m_ptScalarField_dz;
	
	// The isosurface value.
	T m_tIsoLevel;

	// Indicates whether a valid surface is present.
	bool m_bValidSurface;

	// Lookup tables used in the construction of the isosurface.
    static const unsigned int m_edgeTable[256];
    static const int m_triTable[256][16];

};

extern template class CIsoSurface<float>;

#endif // CISOSURFACE_H

