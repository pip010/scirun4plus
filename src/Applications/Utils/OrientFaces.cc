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
//    File   : OrientFaces.cc
//    Author : Martin Cole
//    Date   : Wed Mar 30 13:02:39 2005

#include <Core/Basis/TriLinearLgn.h>
#include <Core/Datatypes/TriSurfMesh.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Datatypes/Field.h>


#include <iostream>
#include <fstream>


#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;
typedef TriSurfMesh<TriLinearLgn<Point> > TSMesh;

// generate a mesh with bad orientations to test the fix with.
void randomizie_orientation(TSMesh *tsm) {
  TSMesh::Face::iterator fiter, fend;
  tsm->begin(fiter);
  tsm->end(fend);
  srand(69);
  while(fiter != fend) {
    if (rand() % 3) {
      tsm->flip_face(*fiter);
    }
    ++fiter;
  }
}


int
main(int argc, char **argv) {
  if (argc != 3 && argc != 4) {
    cerr << "Usage: OrientFaces input-trisurf.fld [-randomize | -flip] output-trisurf.fld\n";
    exit(1);
  }
  FieldHandle handle;
  PiostreamPtr stream = auto_istream(argv[1]);
  if (!stream) {
    cerr << "Couldn't open file " << argv[1] << ".  Exiting..." << endl;
    exit(1);
  }
  Pio(*stream, handle);
  if (!handle.get_rep()) {
    cerr << "Error reading surface from file " << argv[1] << ".  Exiting..." 
	 << endl;
    exit(2);
  }
  
  MeshHandle mb = handle->mesh();
  TSMesh *tsm = dynamic_cast<TSMesh *>(mb.get_rep());
  if (! tsm) { cerr << "Error: not a TriSurf." << endl; return 99;}


  std::string fout;
  if (argc == 4) {
    if (std::string(argv[2]) == std::string("-randomize")) {
      randomizie_orientation(tsm);
      fout = argv[3];
    } if (std::string(argv[2]) == std::string("-flip")) {
      tsm->flip_faces();
      fout = argv[3];
    } else {
      cerr << "Error: OrientFaces got 3 arguments, but the second one wasn't -randomize or -flip\n";
      exit(3);
    }
  } else {
    fout = argv[2];
    tsm->orient_faces();
  }

  BinaryPiostream out_stream(fout, Piostream::Write);
  Pio(out_stream, handle);
  return 0;  
}    
