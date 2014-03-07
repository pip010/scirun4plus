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
//    File   : ExtractIsosurface.cc
//    Author : David Weinstein
//    Date   : Mon Jan  6 10:22:37 MST 2009

#include <Core/Algorithms/Fields/MarchingCubes/MarchingCubes.h>
#include <Core/Datatypes/Field.h>
#include <Core/Persistent/Pstreams.h>

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

int
main(int argc, char *argv[]) {
  if (argc < 4) {
    cerr << "Usage: "<<argv[0]<<" inputFieldFile isoSurfaceFieldFile isoVal1 [isoVal2] [isoVal3] [...]\n";
    exit(1);
  }
  FieldHandle inputH;
  PiostreamPtr inputStream = auto_istream(argv[1]);
  if (!inputStream) {
    cerr << "Couldn't open file " << argv[1] << ".  Exiting..." << endl;
    exit(1);
  }
  Pio(*inputStream, inputH);
  if (!inputH.get_rep()) {
    cerr << "Error reading field from file " << argv[1] << ".  Exiting..." 
	 << endl;
    exit(2);
  }
  std::vector<double> isovals;
  for (int i=3; i<argc; i++) {
    double x;
    if (sscanf(argv[i], "%lf", &x) != 1) {
      cerr << "Error: couldn't parse isoval "<<argv[i]<<" as double.\n";
      exit(3);
    }
    isovals.push_back(x);
  }

  SCIRunAlgo::MarchingCubesAlgo algo;
  algo.set_colormap("colormap", 0);
  algo.set_bool("build_geometry", false);
  algo.set_bool("build_field", true);
  algo.set_bool("build_node_interpolant", false);
  algo.set_bool("build_elem_interpolant", false);
  algo.set_int("num_threads", 0); // 0 for multi-threading, 1 for single
  FieldHandle outputH;
  GeomHandle geomH;
  MatrixHandle nodeH, elemH;
  if (!algo.run(inputH, isovals, geomH, outputH, nodeH, elemH)) {
    cerr << "Isosurface Extraction algorithm failed.\n";
    exit(4);
  }
  BinaryPiostream outStream(argv[2], Piostream::Write);
  Pio(outStream, outputH);
  return 0;  
}    
