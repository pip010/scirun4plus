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
//    File   : FieldToPtcl.cc
//    Author : David Weinstein, Bill Martin
//    Date   : Thurs Feb  5 7:20:00 MST 2009

#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/PointCloudMesh.h>
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
main(int argc, char* argv[])
{
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " field_filename ptcl_filename\n";
    exit(-1);
  }

  PiostreamPtr stream=auto_istream(argv[1]);
  if (!stream) {
    cerr << "Couldn't open file "<<argv[1]<<".  Exiting...\n";
    exit(0);
  }
  FieldHandle pcfld;
  Pio(*stream, pcfld);
  if (!pcfld.get_rep()) {
    cerr << "Error reading surface from file "<<argv[1]<<".  Exiting...\n";
    exit(0);
  }
  VMesh *vmesh = pcfld->vmesh();
  VMesh::Node::size_type nnodes;
  vmesh->size(nnodes);
  Point *pts = vmesh->get_points_pointer();
  
  // open the ptcl file, write the header, ...
  char *fout = argv[2];
  FILE* fp = fopen(fout,"w");

  if (!fp) {
    cerr << "Error opening output file " << fout << ". Exiting...\n";
    exit(0);
  }

  fprintf(fp,"%d\n",(int)nnodes);

  for(int i=0; i<nnodes; i++)
    fprintf(fp, "%lf %lf %lf 0.0 0.0 0.0 0.0\n", pts[i].x(), pts[i].y(), pts[i].z());

  fclose(fp);

}
