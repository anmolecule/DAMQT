//  Copyright 2008-2016, Jaime Fernandez Rico, Rafael Lopez, Ignacio Ema,
//  Guillermo Ramirez, David Zorrilla
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
//	Buffer for storing the 3D grid data
//
//  File:   VVBuffer.cpp
//
//      Last version: March 2016
//
#include "VVBuffer.h"

static void
nnresample (float *in, int in_x, int in_y, int in_z,float *out, int out_x, int out_y, int out_z)
{
  for (int oz = 0; oz < out_z; oz++)
    for (int oy = 0; oy < out_y; oy++)
      for (int ox = 0; ox < out_x; ox++)
	{
	  int iidx = (
		      ((oz * in_z / out_z) * in_y +
		       (oy * in_y / out_y)) * in_x +
		      (ox * in_x / out_x) 
		      );
	  *out++ = in[iidx];
	}
}

void VVBuffer::resize_to (int maxsz) 
{
  if (maxsz > size[0] && maxsz > size[1] && maxsz > size[2])
    return;
  double f = maxsz / (double)((size[0] > size[1]) 
			      ? ((size[0] > size[2]) ? size[0] : size[2])
			      : ((size[1] > size[2]) ? size[1] : size[2]));
  int nx = (int)(f * size[0]), ny = (int)(f * size[1]), nz = (int)(f * size[2]);
  float *nbuf = new float[nx * ny * nz];
  nnresample(buf, size[0], size[1], size[2],
	     nbuf, nx, ny, nz);
  delete buf;
  buf = nbuf;
  size[0] = nx; size[1] = ny; size[2] = nz;
}

bool VVBuffer::acceptable_size (int maxsz)
{
	if (maxsz > size[0] && maxsz > size[1] && maxsz > size[2])
		return true;
	else
		return false;
}