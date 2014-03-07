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


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <teem/nrrd.h>

using std::endl;
using std::cerr;

int
main(int argc, char* argv[])
{
  if (argc != 3) {
    cerr << "Usage: "<<argv[0]<<" nrrdFile vffFile\n";
    exit(1);
  }
  Nrrd *nrrd = nrrdNew();
  if (nrrdLoad(nrrd, airStrdup(argv[1]), 0)) {
    char *err = biffGetDone(NRRD);      
    cerr << "Error reading nrrd " << argv[1] << ". Error: "<< err << endl;
    free(err);
    biffDone(NRRD);
    exit(2);
  }
  
  if (nrrd->dim != 3) {
    cerr << "Error -- nrrd has to be 3D.\n";
    exit(3);
  }

  FILE* out_file;
  if ((out_file = fopen(argv[2], "wb")) == NULL) {
    cerr << "Error: output file open failed\n";
    exit(4);
  }
  int nbits=8;
  if (nrrd->type == nrrdTypeShort) nbits=16;
  else if (nrrd->type == nrrdTypeUChar) nbits=8;
  else {
    cerr << "Error: can only convert uchar or short NRRDs to VFF format\n";
    exit(5);
  }
  for (int i=0; i<3; i++) {
    if (airIsNaN(nrrd->axis[i].spacing)) {
      cerr << "WARNING: input NRRD didn't have spacing for axis "<<i<<" (or explicitly had spacing as NaN) -- setting spacing to 1.0 in VFF\n";
      nrrd->axis[i].spacing = 1.0;
    }
  }
  fprintf(out_file, "ncaa\n");
  fprintf(out_file, "rank=3;\n");
  fprintf(out_file, "type=raster;\n");
  fprintf(out_file, "size=%d %d %d;\n", (int)nrrd->axis[0].size, (int)nrrd->axis[1].size, (int)nrrd->axis[2].size);
  fprintf(out_file, "bits=%d;\n", nbits);
  fprintf(out_file, "format=slice;\n");
  fprintf(out_file, "spacing=%lf %lf %lf;\n", nrrd->axis[0].spacing, nrrd->axis[1].spacing, nrrd->axis[2].spacing);
  fprintf(out_file, "origin=0.0 0.0 0.0;\n");
  fprintf(out_file, "elementsize=1.0;\n", nrrd->axis[0].spacing);
  fprintf(out_file, "\n");
  if (nbits == 16 && AIR_ENDIAN == 1234) nrrdSwapEndian(nrrd);
  fwrite(nrrd->data, nbits/8,
	 nrrd->axis[0].size * nrrd->axis[1].size * nrrd->axis[2].size,
	 out_file);
  fclose(out_file);
  return(0);
}
