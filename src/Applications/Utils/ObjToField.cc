/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2010 Scientific Computing and Imaging Institute,
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


/*
 *  ObjToField.cc
 *
 *  Written by:
 *   Darrell Swenson
 *   Based on PtclToField.cc
 *   University of Utah
 *   Feb 2010
 *
 */

#include <Core/Algorithms/DataIO/ObjToFieldReader.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Datatypes/MatrixTypeConverter.h>

#include <iostream>

using namespace SCIRun;
using namespace SCIRunAlgo;

int
main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " obj_filename field_filename\n";
    exit(-1);
  }

  char *fin, *fout;
  fin = argv[1];
  fout = argv[2];

  FieldHandle fH = 0;
  SCIRunAlgo::ObjToFieldReader reader(0);
  std::string fn(fin);
  if (! reader.read(fn, fH))
  {
    std::cerr << "Convert Obj to field failed." << std::endl;
    exit(-1);
  }

  BinaryPiostream out_stream(fout, Piostream::Write);
  Pio(out_stream, fH);

  return 0;  
}    
