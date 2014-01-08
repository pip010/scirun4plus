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
    cerr << "Usage: "<<argv[0]<<" nrrdFile volFile\n";
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
  if ((out_file = fopen(argv[2], "w")) == NULL) {
    cerr << "Error: output file open failed\n";
    exit(4);
  }
  if (nrrd->type != nrrdTypeFloat) {
    Nrrd *nrrd2 = nrrdNew();
    if (nrrdConvert(nrrd2, nrrd, nrrdTypeFloat)) {
      char *err = biffGet(NRRD);
      cerr << "Error converting NRRD to float.\n";
      exit(5);
      free(err);
      biffDone(NRRD);
    }
    nrrdNuke(nrrd);
    nrrd=nrrd2;
  }
  fprintf(out_file, "4\n %d %d %d \n", nrrd->axis[0].size,
	  nrrd->axis[1].size, nrrd->axis[2].size);
  
  fwrite(nrrd->data, sizeof(float),
	 nrrd->axis[0].size * nrrd->axis[1].size * nrrd->axis[2].size,
	 out_file);
  
  fclose(out_file);
  return(0);
}
