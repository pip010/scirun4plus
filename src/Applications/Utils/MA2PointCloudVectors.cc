/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/



// This program will read in a .pts (specifying the x/y/z coords of each 
// point, one per line, entries separated by white space, file can have 
// an optional one line header specifying number of points... and if it
// doesn't, you have to use the -noPtsCount command-line argument).  
// And the SCIRun output file is written in ASCII, unless you specify 
// -binOutput.

#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/Field.h>
#include <Core/Persistent/Pstreams.h>

#include <Core/Util/FileUtils.h>
#include <Core/Init/init.h>

#include <sci_deprecated.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace SCIRun;

int
main(int argc, char **argv) 
{

  if (argc != 3) 
  {
    std::cerr << "MA2PointCloudVectors <input file (_ma.ptcl)> <output name>" 
         << std::endl;
    return 2;
  }
  SCIRunInit();

  MeshHandle mesh_handle = CreateMesh(POINTCLOUDMESH_E);
  VMesh *pcm = mesh_handle->vmesh();

  //! first line is number of points.
  std::ifstream in(argv[1]);
  if (! in) 
  {
    std::cerr << "could not open " << argv[1] << " exiting." << std::endl;
    return 2;    
  }

  int npts;
  in >> npts;
  
  std::vector<Vector> vecs;
  for (int i = 0; i < npts; i++) 
  {
    double x, y, z, vx, vy, vz, s;
    in >> x >> y >> z >> vx >> vy >> vz >> s;
    pcm->add_point(Point(x,y,z));
    vecs.push_back(Vector(vx, vy, vz));
  }

  //! build field after mesh is done so that its data storage 
  //! gets sized correctly.
  FieldInformation fi(POINTCLOUDMESH_E,CONSTANTDATA_E,VECTOR_E);
  FieldHandle pcH = CreateField(fi,mesh_handle);
  VField* pc = pcH->vfield();
  
  for (int i = 0; i < npts; i++) 
  {
    pc->set_value(vecs[i], VField::index_type(i));
  }

  BinaryPiostream out_stream(argv[2], Piostream::Write);
  Pio(out_stream, pcH);

  return 0;  
}    
