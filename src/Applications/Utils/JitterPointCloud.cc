//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//    File   : JitterPointCloud.cc
//    Author : Jeroen Stinstra and Bill Martin
//    Date   : Wed Dec 2 2009

#include <Core/Algorithms/Converter/ConverterAlgo.h>
#include <Core/Algorithms/Fields/ClipMesh/ClipMeshBySelection.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#ifdef _WIN32
#include <Core/Util/Rand.h>
#endif

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::cout;
using std::ifstream;
using std::endl;

using namespace SCIRun;

Point
compute_jittered_point(const Point& pt, double jitter_rad)
{
	Vector offset_vec;
	do {
		offset_vec = Vector(2*drand48()-1, // uniform random b/t [-1,1]
							2*drand48()-1,
							2*drand48()-1);
	} while (offset_vec.length2() > 1); // eliminate points outside sphere
	cout << "Old pt: " << pt << " New pt: " << pt+offset_vec*jitter_rad << " Offset: " << offset_vec << endl;

	return pt + offset_vec*jitter_rad;
}

int
main(int argc, char *argv[]) {
  if (argc < 4) {
    cerr << "Usage: "<<argv[0]<<" jitterPercentage pointField jitteredPointField\n";
    exit(1);
  }

  double jit_percent = atof(argv[1]);

  PiostreamPtr stream=auto_istream(argv[2]);
  if (!stream) {
    cerr << "Couldn't open file "<<argv[2]<<".  Exiting...\n";
    exit(0);
  }
  FieldHandle pcfld;
  Pio(*stream, pcfld);
  if (!pcfld.get_rep()) {
    cerr << "Error reading surface from file "<<argv[2]<<".  Exiting...\n";
    exit(0);
  }
  VMesh *vmesh = pcfld->vmesh();
  VMesh::Node::size_type nnodes;
  vmesh->size(nnodes);
  Point *pts = vmesh->get_points_pointer();

  BBox pts_bbox = vmesh->get_bounding_box();
  double diag_len = pts_bbox.diagonal().length();
  double jit_rad = diag_len * jit_percent;

  for (VMesh::Node::size_type i=0; i<nnodes; i++)
  {
	  pts[i] = compute_jittered_point(pts[i],jit_rad);
  }

  BinaryPiostream jitPointFieldStream(argv[3], Piostream::Write);
  Pio(jitPointFieldStream, pcfld);
  return 0;


}
