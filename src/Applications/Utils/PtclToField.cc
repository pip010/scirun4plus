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


/*
 *  PtclToField.cc
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2009
 *
 */


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

int readLine(FILE **f, char *buf) {
    char c = 0;
    int cnt=0;
    do {
      if(!feof(*f) && ((c=fgetc(*f))=='#')) {
	while (!feof(*f) && ((c=fgetc(*f))!='\n'));
      }
    } while (c=='\n');
    if (feof(*f)) return 0;
    buf[cnt++]=c;
    while (!feof(*f) && ((c=fgetc(*f))!='\n')) buf[cnt++]=c;
    buf[cnt]=c;
    if (feof(*f)) return 0;
    return 1;
}

int
main(int argc, char **argv) {
  if (argc != 3 && argc != 4) {
    cerr << "Usage: " << argv[0] << " [-nodata] ptcl_filename field_filename\n";
    exit(-1);
  }
  bool nodata=false;
  if (argc == 4) {
    if (std::string(argv[1]) == std::string("-nodata")) {
      nodata=true;
    } else {
      cerr << "Error: found three arguments but first one wasn't '-nodata'\n";
      exit(-1);
    }
  }
  FieldInformation fi("PointCloudMesh", "ConstantBasis", "Vector");
  if (nodata) fi.set_data_basis_type("nodata");
  std::vector<Vector> values;
  MeshHandle mH=CreateMesh(fi);
  char *fin, *fout;
  if (nodata) fin=argv[2]; else fin=argv[1];
  if (nodata) fout=argv[3]; else fout=argv[2];

  FILE *f=fopen(fin, "rt");
  if (!f) {
    cerr << "Error - failed to open input file: "<<fin<<"\n";
    exit(-1);
  }
  char buf[2048];
  readLine(&f, buf);
  int npts;
  sscanf(buf, "%d", &npts);
  float x, y, z, nx, ny, nz, r;
  for (int i=0; i<npts; i++) {
    readLine(&f, buf);
    sscanf(buf, "%f %f %f %f %f %f %f", &x,&y, &z, &nx, &ny, &nz, &r);
    mH->vmesh()->add_point(Point(x,y,z));
    if (!nodata) values.push_back(Vector(nx, ny, nz));
  }
  fclose(f);
  FieldHandle fH=CreateField(fi, mH);
  if (!nodata) fH->vfield()->set_values(values);
  BinaryPiostream out_stream(fout, Piostream::Write);
  Pio(out_stream, fH);
  return 0;  
}    
