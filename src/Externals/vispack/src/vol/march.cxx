// ******************************************************************
// VISPack. Copyright (c) 1994-2000 Ross Whitaker rtw@utk.edu       *
// For conditions of distribution and use, see the file LICENSE.txt *
// accompanying this distribution.                                  *
// ******************************************************************

// $Id: march.cxx,v 1.2 2003/04/30 01:32:14 whitaker Exp $


#include "vol/volume.h"
#include "vol/march.h"
#include "vol/meshtable.h"
#include "util/list.h"
#include "vol/volindexlist.h"
#include "util/array.h"
#include <limits.h>
#include <float.h>
//#include <striostream>
//#define MARCH_DEBUG

//***************************************************************************
//
//  void Marching_Cubes(float *grid, size grid_size, float iso_value, char *nameOut);
//  Use the marching cubes algorithm to extract a surface from the given volume at
//  the value of iso_value.
//
//***************************************************************************

void Marching_Cubes(float *grid, size grid_size, float iso_value, char *nameOut, 
		    float x_scale, float y_scale, float z_scale)
{
    const float normCutoff = 1000.0f*FLT_MIN;  // Discard triangles with norm smaller than this
    int bits[12];            // Holds a mask for the first 12 bits
    int w, x, y;             // For loop variables
    int vert;                // Holds which vertexes are above threshold
			     // and which are below
    vox_edge voxel_edge[12]; // Defines edges by which vertices they
				// connect
    vox_mesh *mesh;          // Pointer to table holding triangles that
				// need to be created for each voxel 
				// configuration
    float vertex[8];         // Value at voxel vertices
    location rel_vert[8];
    ofstream out_file;       // File for output
    int a, b, c;
    int num_triangles;       // Number of triangles for particular
				// voxel during loop
    long count = 0;
    int edge1, edge2, edge3; // Edges where vertices are located for
				// this triangle being constructed
    int vert1_a, vert2_a;
    int vert1_b, vert2_b;
    int vert1_c, vert2_c;
    int x_size, y_size, z_size;
    float between_a;         // Interpolated values where vertices are
    float between_b;         // to be located along given edges
    float between_c;
    float cell_size = 1;     // Relative length of each cell
    mesh = new vox_mesh [256]; // Allocate space for table
    
    // VISArrays to hold vector data, indexed by vertex number, one array each
    // to return x, y, and z data for given normal vector, define length of
    // array to be 4 times number of original volume vertices in one xy slice,
    // array will automatically grow if needed
    long volSlice = grid_size.x * grid_size.y;
    VISArray<float> vectorX(volSlice*4);
    VISArray<float> vectorY(volSlice*4);
    VISArray<float> vectorZ(volSlice*4);

    // Three volumes as indices into vector arrays, used to look for
    // matching vertices, use smallest type possible. Set default (i.e.
    // nothing present) to zero, therefore need to start counting
    // triangles with number 1. Hence the vertex number stored will be one
    // greater than the actual value
    VISVolume<unsigned> volX(grid_size.x, grid_size.y, grid_size.z);
    VISVolume<unsigned> volY(grid_size.x, grid_size.y, grid_size.z);
    VISVolume<unsigned> volZ(grid_size.x, grid_size.y, grid_size.z);
    volX = (unsigned)0;
    volY = (unsigned)0;
    volZ = (unsigned)0;
    // Pointers used to reference along which edge lies a given vertex,
    // i.e. in x, y, or z direction, point these to volX, volY, or volZ as
    // necessary to set edge for vertex 1, 2, or 3
    VISVolume<unsigned> *vol1;
    VISVolume<unsigned> *vol2;
    VISVolume<unsigned> *vol3;

    // VISLists to keep track of location of each vertex, and which
    // vertices compose a given triangle
    VISList <VISTriple<float> > vertices;
    VISList <VISTriple<long> > triangles;
    // Total number of vertices and triangles
    long numVertices = 0, numTriangles = 0;

    // Initialize mesh table etc.
    Set_Stuff(bits, voxel_edge, rel_vert);
    Read_Voxel_Mesh_Table(mesh);

    // Set size of volume
    x_size = grid_size.x;
    y_size = grid_size.y;
    z_size = grid_size.z;
	
    // Temp holders for vertices of current triangle
    VISTriple<float> ver1, ver2, ver3, V1, V2, V3;
    // Norm of cross product
    float V3norm;                        
    // Vertex numbers of current triangles' vertices, and triangle holder
    long triVert1, triVert2, triVert3;
    VISTriple<long> tri;
    //Open iv file that I will write the meshes into
    out_file.open(nameOut);

#ifdef MARCH_DEBUG
    // Open another file to hold vertex positions, display these superimposed
    // over surface for debugging purposes
    ofstream outNorm("Normals.iv");
    outNorm << "#Inventor V2.0 ascii"<< endl<< endl
	    << "Separator {"<<endl
	    << "  DrawStyle {"<<endl
	    << "    style POINTS"<<endl
	    << "    pointSize 3 }"<<endl
	    << "     Coordinate3 {"<< endl
	    << "          point [" << endl;
#endif
    
    //Write the header to the iv file
    out_file << "#Inventor V2.0 ascii"<< endl<< endl
             << "Separator {"<<endl
             << "     ShapeHints {"<< endl
             << "         vertexOrdering CLOCKWISE"<<endl
	     << "         creaseAngle 0 }"<<endl
             << "     Coordinate3 {"<< endl
             << "          point [" << endl;

    // Loop through entire volume, use table to find triangles, vertices etc.
    for (a=0; a<(grid_size.x-1); a++) {
	for (b=0; b<(grid_size.y-1); b++) {
	    for (c=0; c<(grid_size.z-1); c++) {
		vert = 0;
	
		// Assign vertices of voxel to vertex[.] variable.
		vertex[0] = grid[a + b*x_size + c*x_size*y_size];
		vertex[1] = grid[a + (b+1)*x_size + c*x_size*y_size];
		vertex[2] = grid[(a+1) + (b+1)*x_size + c*x_size*y_size];
		vertex[3] = grid[(a+1) + b*x_size + c*x_size*y_size];
		vertex[4] = grid[a + b*x_size + (c+1)*x_size*y_size];
		vertex[5] = grid[a + (b+1)*x_size + (c+1)*x_size*y_size];
		vertex[6] = grid[(a+1) + (b+1)*x_size + (c+1)*x_size*y_size];
		vertex[7] = grid[(a+1) + b*x_size + (c+1)*x_size*y_size];

		// Use bit mask to create binary number that tells which vertices
		// are above the iso_value threshold
		for (w=0; w<8; w++) {
		    if (vertex[w] > iso_value) {
			vert = vert | bits[w];
		    }
		}	
	     
		// Number of triangles for this particular voxel
		num_triangles = mesh[vert].num_tri;
    
		for (x=0; x<num_triangles; x++) {
		    edge1 = mesh[vert].meshes[x].e1;   // Set edges along which
		    edge2 = mesh[vert].meshes[x].e2;   // vertices exist
		    edge3 = mesh[vert].meshes[x].e3;
		    
		    vert1_a = voxel_edge[edge1-1].vertex_one-1;
		    vert2_a = voxel_edge[edge1-1].vertex_two-1;
			  
		    vert1_b = voxel_edge[edge2-1].vertex_one-1;
		    vert2_b = voxel_edge[edge2-1].vertex_two-1;
			  
		    vert1_c = voxel_edge[edge3-1].vertex_one-1;
		    vert2_c = voxel_edge[edge3-1].vertex_two-1;

		    // Interpolate to find vertex locations along edges
		    between_a = cell_size*(iso_value - vertex[vert1_a])/
			(vertex[vert2_a] - vertex[vert1_a]);
		    between_b = cell_size*(iso_value - vertex[vert1_b])/
			(vertex[vert2_b] - vertex[vert1_b]);
		    between_c = cell_size*(iso_value - vertex[vert1_c])/
			(vertex[vert2_c] - vertex[vert1_c]);
		    
		    if (voxel_edge[edge1-1].inc_dir == 1) {
			// Set vertex 1 location
			ver1 = VISTriple<float>(between_a + rel_vert[vert1_a].x + a,
						 rel_vert[vert1_a].y + b,
						 rel_vert[vert1_a].z + c);
			vol1 = &volX;    // Set pointer to X edge
		    } else if (voxel_edge[edge1-1].inc_dir == 2) {
			ver1 = VISTriple<float>(rel_vert[vert1_a].x + a,
						 between_a + rel_vert[vert1_a].y + b,
						 rel_vert[vert1_a].z + c);
			vol1 = &volY;;    // Set pointer to Y edge
		    } else if (voxel_edge[edge1-1].inc_dir == 3) {
			ver1 = VISTriple<float>(rel_vert[vert1_a].x + a,
						 rel_vert[vert1_a].y + b,
						 between_a + rel_vert[vert1_a].z + c);
			vol1 = &volZ;;    // Set pointer to Z edge
		    }
	
		    if (voxel_edge[edge2-1].inc_dir == 1) {
			// Set vertex 2 location
			ver2 = VISTriple<float>(between_b + rel_vert[vert1_b].x + a,
						 rel_vert[vert1_b].y + b,
						 rel_vert[vert1_b].z + c);
			vol2 = &volX;
		    } else if (voxel_edge[edge2-1].inc_dir == 2) {
			ver2 = VISTriple<float>(rel_vert[vert1_b].x + a,
						 between_b + rel_vert[vert1_b].y + b,
						 rel_vert[vert1_b].z + c);
			vol2 = &volY;
		    } else if (voxel_edge[edge2-1].inc_dir == 3) {
			ver2 = VISTriple<float>(rel_vert[vert1_b].x + a,
						 rel_vert[vert1_b].y + b,
						 between_b + rel_vert[vert1_b].z + c);
			vol2 = &volZ;
		    }
	
		    if (voxel_edge[edge3-1].inc_dir == 1) {
			// Set vertex 3 location
			ver3 = VISTriple<float>(between_c + rel_vert[vert1_c].x + a,
						 rel_vert[vert1_c].y + b,
						 rel_vert[vert1_c].z + c);
			vol3 = &volX;
		    } else if (voxel_edge[edge3-1].inc_dir == 2) {
			ver3 = VISTriple<float>(rel_vert[vert1_c].x + a,
						 between_c + rel_vert[vert1_c].y + b,
						 rel_vert[vert1_c].z + c);
			vol3 = &volY;
		    } else if (voxel_edge[edge3-1].inc_dir == 3) {
			ver3 = VISTriple<float>(rel_vert[vert1_c].x + a,
						 rel_vert[vert1_c].y + b,
						 between_c + rel_vert[vert1_c].z + c);
			vol3 = &volZ;
		    }

		    // Calculate triangle normal
		    //Find vectors pointing to vertices, only need 2 of them
		    V1.a(ver1.a()-ver3.a());              
		    V1.b(ver1.b()-ver3.b());
		    V1.c(ver1.c()-ver3.c());
		    V2.a(ver2.a()-ver1.a());
		    V2.b(ver2.b()-ver1.b());
		    V2.c(ver2.c()-ver1.c());
		    // Calculate cross product
		    V3.a(-(V1.b()*V2.c() - V1.c()*V2.b()));  
		    V3.b(-(V1.c()*V2.a() - V1.a()*V2.c()));
		    V3.c(-(V1.a()*V2.b() - V1.b()*V2.a()));
		    // Find norm of the cross product
		    V3norm = sqrt(V3.a()*V3.a() + V3.b()*V3.b() + V3.c()*V3.c()); 
#ifdef MARCH_DEBUG
		    cerr << "Norm: " << V3norm << endl;
#endif

		    // Check value of normal, if near-zero then discard triangle
		    if (V3norm > normCutoff) 
		    {
			// Set vertices (i.e. check if unique etc.)
			incVertex(ver1, out_file, numVertices, vertices, triVert1, vectorX, 
				  vectorY, vectorZ, *vol1, 
				  x_scale, y_scale, z_scale);
			incVertex(ver2, out_file, numVertices, vertices, triVert2, vectorX, 
				  vectorY, vectorZ, *vol2,
				  x_scale, y_scale, z_scale);
			incVertex(ver3, out_file, numVertices, vertices, triVert3, vectorX, 
				  vectorY, vectorZ, *vol3,
				  x_scale, y_scale, z_scale);

			// Set this triangle using vertex numbers found above
			tri.a(triVert1);
			tri.b(triVert2);
			tri.c(triVert3);
#ifdef MARCH_DEBUG
			cerr << "Triangle: "<<triVert1<<" "<<triVert2<<" "<<triVert3<<endl;
#endif
			// Add this triangle to list
			triangles.appendItem(tri);
			numTriangles++;

			// Normalize vector by dividing by norm
			V3.a(V3.a()/V3norm);
			V3.b(V3.b()/V3norm);
			V3.c(V3.c()/V3norm);

#ifdef MARCH_DEBUG
			cerr << "Normal: "<<V3.a()<<" "<<V3.b()<<" "<<V3.c()<<endl;
#endif

			// Add normals to existing normals, new normals already
			// initialized by incVertex
			addNormal(V3, triVert1, vectorX, vectorY, vectorZ);
			addNormal(V3, triVert2, vectorX, vectorY, vectorZ);
			addNormal(V3, triVert3, vectorX, vectorY, vectorZ);
		    }
		    
		} // for (x)
	    } // for (c)
	} // for (b)
    } // for (a)     

    // Send vertex normal output to file
    out_file << "              ]   }" << endl
	     << "NormalBinding {" << endl
	     << "value       PER_VERTEX_INDEXED" << endl
	     << "       }" << endl
	     << "Normal {" << endl
	     << "vector      [" << endl;
	    
    // Loop through normal vector arrays, sending normal corresponding to
    // vertex to output file
#ifdef MARCH_DEBUG
    Link <VISTriple<float> > *vertPtr = vertices.head();
#endif
    
    for (count = 0; count < numVertices; count++) {
	// Normalize vectors
	V3.a(vectorX.itemAt(count)/x_scale);
	V3.b(vectorY.itemAt(count)/y_scale);
	V3.c(vectorZ.itemAt(count)/z_scale);

	

	V3norm = sqrt(V3.a()*V3.a() + V3.b()*V3.b() + V3.c()*V3.c()); 
	V3.a(V3.a()/V3norm);
	V3.b(V3.b()/V3norm);
	V3.c(V3.c()/V3norm);

#ifdef MARCH_DEBUG
	V2 = vertPtr->data();
	outNorm << "        " << V2.a()+V3.a()
		<< " " << V2.b()+V3.b()
		<< " " << V2.c()+V3.c()<<","<<endl;
	vertPtr = vertPtr->next();
#endif
	out_file << "        " << -V3.a()   // Send normal to output file
		 << " " << -V3.b()
		 << " " << -V3.c()
		 << "," << endl;
    }
    
#ifdef MARCH_DEBUG
    outNorm << "          ] }"<<endl
	    << "  FaceSet { }"<<endl
	    << "}" << endl;
    outNorm.close();
#endif

    // Loop through triangles, sending vertices corresponding to triangles to
    // output file
    out_file << "              ]   }"<<endl
             << "     IndexedFaceSet { "<<endl
             << "         coordIndex      [ " << endl;
	     
    Link <VISTriple<long> > *triPtr = triangles.head();
    while (triPtr != NULL) {
      out_file << "        " << (triPtr->data().a())
	       << "," << (triPtr->data().b())
	       << "," << (triPtr->data().c()) << ",-1," << endl;
      triPtr = triPtr->next();
    }

    out_file << "                         ] " << endl
             << "     } " << endl;
#ifdef MARCH_DEBUG
    out_file << "  File { name\"Normals.iv\"}"<<endl;
#endif

    out_file << "}"<<endl;
    
    out_file.close();   
	
    // Delete mesh table
    delete [] mesh;
}


//***************************************************************************
//
// Inc_Vertex: If this vertex already exists then point current triangle's vertex
// to existing vertex number and do not output vertex to file.  (Only loop back
// through 'slice' number of vertices (or until beginning of list), where 'slice'
// represents a guess as to the max number of vertices in one slice of the volume
// (1D slice), around 1.25 times # of vertices in one slice of original volume.
// If the vertex does not exist then output it, add it to list, inc number of
// vertices, and set current triangle's vertex to this new vertex
//
//***************************************************************************

void incVertex (VISTriple<float> vertex, ofstream &output, long &numVertices, 
		VISList <VISTriple<float> > &vertices, long &triVert, 
		VISArray<float> &vectorX, VISArray<float> &vectorY, 
		VISArray<float> &vectorZ, VISVolume<unsigned> &volN, 
		float x_scale, float y_scale, float z_scale) 
{
#ifdef MARCH_DEBUG
    cerr << "Vertex: "<<vertex.a()<<" "<<vertex.b()<<" "<<vertex.c()<<"; ";
#endif

    if (volN.itemAt((int)vertex.a(), (int)vertex.b(), (int)vertex.c()) == 0) {
#ifdef MARCH_DEBUG
	cerr << "Unique" << endl;
#endif
	triVert = numVertices++;     // inc # of vertices
	vertices.appendItem(vertex); // add vertex to list
	output << "            " << 
	  x_scale*vertex.a() << " " << 
	  y_scale*vertex.b() 
	       << " " << z_scale*vertex.c() << ", " << endl; // write output
	vectorX.at(triVert) = 0.0f;       // Set new normal to zero
	vectorY.at(triVert) = 0.0f;
	vectorZ.at(triVert) = 0.0f;
	// Need to add one to vertex number since volume is initialized to zero,
	// default value indicating no vertex present
	volN.at((int)vertex.a(), (int)vertex.b(), (int)vertex.c()) = (unsigned)(triVert + 1);
    }
    else {
#ifdef MARCH_DEBUG
	cerr << "Match" << endl;
#endif
	// Need to subtract one since one added to vertex number when storing it
	// in volume (see above)
	triVert = (long)((int)(volN.itemAt((int)vertex.a(), 
					   (int)vertex.b(), 
					   (int)vertex.c())) - 1);
    }
}


//***************************************************************************
// Add vertex normal to existing normal, if this is a new vertex its normal should
// have been initialized to zero in incVertex function
//***************************************************************************
inline void addNormal (VISTriple<float> normal, long vNum, VISArray<float> &vectorX, VISArray<float> &vectorY, VISArray<float> &vectorZ) {
    vectorX.at(vNum) += normal.a();
    vectorY.at(vNum) += normal.b();
    vectorZ.at(vNum) += normal.c();
}


//***************************************************************************
//
//      void Read_Voxel_Mesh_Table(char *table_name, vox_mesh *mesh);
//
//  Reads the table that tell which triangle to draw for each vertex 
// configuration. Values are stored in mesh[].
//
//***************************************************************************

void Read_Voxel_Mesh_Table(vox_mesh *mesh) {
    
	int count = 0;
	int x;
        int number, bin_num;
	int num_triangle;
	ifstream file_in;

	int datacnt=0;
        while (count < 256) {
	    
	    //Read values from file. "number" and "bin_num" are not
	    //stored. They are used to make finding errors in the
	    //table easier. "num_triangle" tells how many triangles
	    //need to be drawn for each configuration. 256 configurations
	    
	    num_triangle = num_tri[count];
	    mesh[count].meshes = new tri_edge [num_triangle];
	    mesh[count].num_tri = num_triangle;
	
	    //Read in numbers that tell which edges have to be connected
	    //for each triangle.
	    	
	    for (x=0; x<num_triangle; x++) {
		mesh[count].meshes[x].e1 = meshdata[datacnt++];
		//cout << mesh[count].meshes[x].e1 << ",";
		mesh[count].meshes[x].e2 = meshdata[datacnt++];
		//cout << mesh[count].meshes[x].e2 << ",";
		mesh[count].meshes[x].e3 = meshdata[datacnt++];
		//cout << mesh[count].meshes[x].e3 << ",            ";
            }
	    //cout << endl;
	    count++;
	}
	return;    
}


//**********************************************************************
//
//  void Set_Stuff(int *bits, vox_edge *voxel_edge, location *rel_vert)
//
//**********************************************************************

void Set_Stuff(int *bits, vox_edge *voxel_edge, location *rel_vert) {

        //Assign bit masks. Used to create binary numbers that tell
        //which vertices are greater than the iso value and which are
        // below

	bits[0]  = 1;          // Binary 0000 0000 0001
	bits[1]  = 2;          // Binary 0000 0000 0010
	bits[2]  = 4;          // Binary 0000 0000 0100
	bits[3]  = 8;          // Binary 0000 0000 1000
	bits[4]  = 16;         // Binary 0000 0001 0000
	bits[5]  = 32;         // Binary 0000 0010 0000
	bits[6]  = 64;         // Binary 0000 0100 0000
	bits[7]  = 128;        // Binary 0000 1000 0000
	bits[8]  = 256;        // Binary 0001 0000 0000
	bits[9]  = 512;        // Binary 0010 0000 0000
	bits[10] = 1024;       // Binary 0100 0000 0000
	bits[11] = 2048;       // Binary 1000 0000 0000

        // Number the edges of the voxel. Defines which vertices are
	// are connected by which edges. The inc_dir variable tells
	// which direction the edge is going. 1 is x-direction, 
	// 2 is y-direction, 3 is z-direction.

	voxel_edge[0].vertex_one = 1;
	voxel_edge[0].vertex_two = 2;
        voxel_edge[0].inc_dir = 2;
	voxel_edge[1].vertex_one = 2;
	voxel_edge[1].vertex_two = 3;
        voxel_edge[1].inc_dir = 1;
	voxel_edge[2].vertex_one = 4;
	voxel_edge[2].vertex_two = 3;
	voxel_edge[2].inc_dir = 2;
	voxel_edge[3].vertex_one = 1;
	voxel_edge[3].vertex_two = 4;
	voxel_edge[3].inc_dir = 1;
	voxel_edge[4].vertex_one = 5;
	voxel_edge[4].vertex_two = 6;
	voxel_edge[4].inc_dir = 2;
	voxel_edge[5].vertex_one = 6;
	voxel_edge[5].vertex_two = 7;
	voxel_edge[5].inc_dir = 1;
	voxel_edge[6].vertex_one = 8;
	voxel_edge[6].vertex_two = 7;
	voxel_edge[6].inc_dir = 2;
	voxel_edge[7].vertex_one = 5;
	voxel_edge[7].vertex_two = 8;
	voxel_edge[7].inc_dir = 1;
	voxel_edge[8].vertex_one = 1;
	voxel_edge[8].vertex_two = 5;
	voxel_edge[8].inc_dir = 3;
	voxel_edge[9].vertex_one = 2;
	voxel_edge[9].vertex_two = 6;
	voxel_edge[9].inc_dir = 3;
	voxel_edge[10].vertex_one = 4;
	voxel_edge[10].vertex_two = 8;
	voxel_edge[10].inc_dir = 3;
	voxel_edge[11].vertex_one = 3;
	voxel_edge[11].vertex_two = 7;
	voxel_edge[11].inc_dir = 3;
	
	//Defines the relative location of the voxel vertices in
	//relation to the base vertex.

	rel_vert[0].x = 0;
	rel_vert[0].y = 0;
	rel_vert[0].z = 0;
	
	rel_vert[1].x = 0;
	rel_vert[1].y = 1;
	rel_vert[1].z = 0;
	
	rel_vert[2].x = 1;
	rel_vert[2].y = 1;
	rel_vert[2].z = 0;
	
	rel_vert[3].x = 1;
	rel_vert[3].y = 0;
	rel_vert[3].z = 0;
	
	rel_vert[4].x = 0;
	rel_vert[4].y = 0;
	rel_vert[4].z = 1;
	
	rel_vert[5].x = 0;
	rel_vert[5].y = 1;
	rel_vert[5].z = 1;
	
	rel_vert[6].x = 1;
	rel_vert[6].y = 1;
	rel_vert[6].z = 1;
	
	rel_vert[7].x = 1;
	rel_vert[7].y = 0;
	rel_vert[7].z = 1;
	 

}



















