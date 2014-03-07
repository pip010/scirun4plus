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

#include <Core/Algorithms/Converter/ConverterAlgo.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/NrrdData.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

int
main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: "<<argv[0]<<" fieldFile nrrdFile\n";
    exit(1);
  }
  FieldHandle fieldH;
  PiostreamPtr fieldStream = auto_istream(argv[1]);
  if (!fieldStream) {
    cerr << "Couldn't open file " << argv[1] << ".  Exiting..." << endl;
    exit(1);
  }
  Pio(*fieldStream, fieldH);
  if (!fieldH.get_rep()) {
    cerr << "Error reading field from file " << argv[1] << ".  Exiting..." 
	 << endl;
    exit(2);
  }
  ProgressReporter pr;
  SCIRunAlgo::ConverterAlgo algo(&pr);
  NrrdDataHandle nrrdH;
  if (!algo.FieldToNrrd(fieldH, nrrdH)) {
    cerr << "Error converting field to nrrd.\n";
    exit(3);
  }

  NrrdIoState *nio = nrrdIoStateNew();
  // set encoding to be raw
  nio->encoding = nrrdEncodingArray[1];
  // set format to be nrrd
  nio->format = nrrdFormatArray[1];
  // set endian to be endian of machine
  nio->endian = airMyEndian;
  if (AIR_ENDIAN != nio->endian) 
    nrrdSwapEndian(nrrdH->nrrd_);
  if (nrrdSave(argv[2], nrrdH->nrrd_, nio)) {
    char *err = biffGet(NRRD);      
    cerr << "Error writing nrrd " << argv[2] << ". Error: "<< err << endl;
    free(err);
    biffDone(NRRD);
    exit(4);
  }
  return 0;  
}
