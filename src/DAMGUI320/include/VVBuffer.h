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
//  File:   VVBuffer.h
//
//      Last version: March 2016
//
#ifndef _VVBuffer_H
#define _VVBuffer_H 1

#include <limits.h>
#include <stdlib.h>

class VVBuffer
{
	int size[3];
	float voxel_size[3];
	float *buf;
	const char *bufname;
	float min_value, max_value;

public:
	VVBuffer(const char *name, int _x, int _y, int _z)
	{
		bufname = name;
		size[0] = _x;
		size[1] = _y;
		size[2] = _z;
		buf = new float[_x*_y*_z];
		voxel_size[0] = voxel_size[1] = voxel_size[2] = 1.0;
		min_value = SHRT_MIN;
		max_value = SHRT_MAX;
	}

	~VVBuffer() {
		delete buf;
	}

	float * data() {
		return buf;
	}

	const char *name() const {
		return bufname;
	}

	void set_voxel_size (float vx, float vy, float vz) {
		voxel_size[0] = vx;
		voxel_size[1] = vy;
		voxel_size[2] = vz;
	}

	void set_minmax (float min, float max) {
		min_value = min;
		max_value = max;
	}

	float minn() { return min_value; }
	float maxn() { return max_value; }

	int x() { return size[0]; }
	int y() { return size[1]; }
	int z() { return size[2]; }

	float voxel_x() { return voxel_size[0]; }
	float voxel_y() { return voxel_size[1]; }
	float voxel_z() { return voxel_size[2]; }

	void resize_to(int maxsz);
	bool acceptable_size (int maxsz);
};

#endif
