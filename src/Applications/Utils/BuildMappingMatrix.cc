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
//    File   : ConvertFieldToNrrd.cc
//    Author : David Weinstein
//    Date   : Mon Jan  5 14:22:37 MST 2009

#include <Core/Algorithms/Fields/Mapping/BuildMappingMatrix.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>
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
  if (argc != 4) {
    cerr << "Usage: "<<argv[0]<<" sourceFieldFile destFieldFile outputMappingMatrix\n";
    exit(-1);
  }
  FieldHandle sourceH;
  PiostreamPtr sourceStream = auto_istream(argv[1]);
  if (!sourceStream) {
    cerr << "Couldn't open file " << argv[1] << ".  Exiting..." << endl;
    exit(-1);
  }
  Pio(*sourceStream, sourceH);
  if (!sourceH.get_rep()) {
    cerr << "Error reading field from file " << argv[1] << ".  Exiting..." 
	 << endl;
    exit(-1);
  }
  FieldHandle destH;
  PiostreamPtr destStream = auto_istream(argv[2]);
  if (!destStream) {
    cerr << "Couldn't open file " << argv[2] << ".  Exiting..." << endl;
    exit(-1);
  }
  Pio(*destStream, destH);
  if (!destH.get_rep()) {
    cerr << "Error reading field from file " << argv[2] << ".  Exiting..." 
	 << endl;
    exit(-1);
  }
  SCIRunAlgo::BuildMappingMatrixAlgo algo;
  algo.set_option("method", "interpolateddata");
  algo.set_scalar("max_distance", -1);
  MatrixHandle matrixH;
  if (!algo.run(sourceH, destH, matrixH)) {
    cerr << "Error generating mapping matrix.\n";
    exit(-1);
  }
  BinaryPiostream matrixStream(argv[3], Piostream::Write);
  Pio(matrixStream, matrixH);
  return 0;
}
