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
 *  ColorMapToText.cc
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 2003
 *
 */

// This program will read in a SCIRun ColorMap, and will save
// it out to a text version: a .txt file.
// The .txt file will contain one data value per line (r g b a t,
// white-space separated).
// The file will also have a one line header, specifying the number
// of colors, unless the user specifies the -noHeader command-line 
// argument.

#include <Core/Datatypes/ColorMap.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Util/FileUtils.h>
#include <Core/Init/init.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

bool header;

void setDefaults() {
  header=true;
}

int parseArgs(int argc, char *argv[]) {
  int currArg = 3;
  while (currArg < argc) {
    if (!strcmp(argv[currArg],"-noHeader")) {
      header=false;
      currArg++;
    } else {
      cerr << "Error - unrecognized argument: "<<argv[currArg]<<"\n";
      return 0;
    }
  }
  return 1;
}

void printUsageInfo(char *progName) {
  cerr << "\n Usage: "<<progName<<" ColorMap textfile [-noHeader]\n\n";
  cerr << "\t This program will read in a SCIRun ColorMap, and will save \n";
  cerr << "\t it out to a text version: a .txt file. \n";
  cerr << "\t The .txt file will contain one data value per line (r g b a \n";
  cerr << "\t t, white-space separated).\n";
  cerr << "\t The file will also have a one line header, specifying the\n";
  cerr << "\t number of colors, unless the user specifies the -noHeader \n";
  cerr << "\t command-line argument.\n\n";
}

int
main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    printUsageInfo(argv[0]);
    return 2;
  }
  SCIRunInit();
  setDefaults();

  char *colormapName = argv[1];
  char *textfileName = argv[2];
  if (!parseArgs(argc, argv)) {
    printUsageInfo(argv[0]);
    return 2;
  }

  ColorMapHandle handle;
  PiostreamPtr stream=auto_istream(colormapName);
  if (!stream) {
    cerr << "Couldn't open file "<<colormapName<<".  Exiting...\n";
    return 2;
  }
  Pio(*stream, handle);
  if (!handle.get_rep()) {
    cerr << "Error reading colormap from file "<<colormapName<<".  Exiting...\n";
    return 2;
  }

  const unsigned int size = handle->resolution();
  cerr << "Number of colors = "<<size<<"\n";
  FILE *fTxt = fopen(textfileName, "wt");
  if (!fTxt) {
    cerr << "Error -- couldn't open output file "<<textfileName<<"\n";
    return 2;
  }
  if (header) fprintf(fTxt, "%d\n", size);
  //double alpha;
  for (unsigned int c=0; c<size; c++)
  {
    const double t = c / (size - 1.0);
    float r,g,b,a;
    handle->get_color(t, r, g, b, a);
    fprintf(fTxt, "%lf %lf %lf %lf %lf\n", r, g, b, a, t);
  }
  fclose(fTxt);
  return 0;  
}    
