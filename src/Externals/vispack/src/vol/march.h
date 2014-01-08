// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: march.h,v 1.2 2003/04/30 01:32:14 whitaker Exp $

#ifndef vis_march_h
#define vis_march_h

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "util/list.h"
#include "vol/volindexlist.h"
#include "vol/volume.h"

struct size {
	int x;
	int y;
	int z;
};

struct location {
	int x;
	int y;
	int z;
};

struct vox_edge {
	int vertex_one;
	int vertex_two;
        int inc_dir;
};

struct triangle {
	float v1;
	float v2;
	float v3;
};

struct tri_edge {
	int e1;
	int e2;
	int e3;
};

struct vox_mesh {
	int num_tri;
	tri_edge *meshes;
};

//Declare Functions

void Set_Stuff(int *bits, vox_edge *voxel_edge, location *rel_vert);		
void Read_Voxel_Mesh_Table(vox_mesh *mesh);
void Marching_Cubes(float *grid, size grid_size, float iso_value, char *nameOut, float x_scale = 1.0, 
		    float y_scale = 1.0, float z_scale = 1.0);
void incVertex (VISTriple<float> vertex, ofstream &output, long &numVertices, 
		VISList <VISTriple<float> > &vertices, long &triVert, 
		VISArray<float> &vectorX, VISArray<float> &vectorY, 
		VISArray<float> &vectorZ, VISVolume<unsigned int> &volN, 
		float x_scale = 1.0, float y_scale = 1.0, float z_scale = 1.0
		);
inline void addNormal (VISTriple<float> normal, long vNum, VISArray<float> &vectorX, VISArray<float> &vectorY, VISArray<float> &vectorZ);

#endif











