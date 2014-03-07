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
 *  vff2nrrd.cc: create a nrrd header for a vff file
 *
 *  Written by:
 *   Mark Hartner
 *   SCI Institute
 *   University of Utah
 *   December 2004
 *
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <string.h>
#include <stdio.h>

using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::ifstream;
using std::ofstream;

int
main(int argc, char **argv)
{
  int nrrdSkip = 0;
  int dim=-1;
  int i;
  int size[3];
  int nbits = 16;
  char *in = argv[1];
  ifstream vffFileStream(in, std::ios::binary);
  char *out = argv[2];
  ofstream nrrdHeader(out, std::ios::binary);
  float spacing[3];

  if (argc < 3)
  {
    cerr << "Usage: "<<argv[0]<<" inputFile.vff outputFile.nhdr \n";
    return 2;
  }

  if (! vffFileStream)
  {
    cout<<"could not input vff file "<<in<<endl;
    return 1;
  }
  if (! nrrdHeader)
  {
    cout<<"could not create nrrd header file "<<out<<endl;
    return 1;
  }

  spacing[0] = spacing[1] = spacing[2] = 1.0;

  int lineCount = 0;
  int foundFormFeed = 0;	
  char temp[1025];
  while (! vffFileStream.eof() && !foundFormFeed)
  {
    vffFileStream.getline(temp,1024);
    lineCount++;

    if (temp[0] == '\f')
    {
      nrrdSkip = lineCount;
      foundFormFeed = 1;
    }
	  
    if (strncmp(temp,"rank=",5) == 0)
    {
      std::istringstream rankLine(&temp[5]);
      rankLine >> dim;
    }
	  
    if (strncmp(temp,"bits=",5) == 0)
    {
      std::istringstream bitsLine(&temp[5]);
      bitsLine >> nbits;
    }
	  
    if (strncmp(temp,"spacing=",8) == 0)
    {
      std::istringstream spacingLine(&temp[8]);
      spacingLine >> spacing[0] >> spacing[1] >> spacing[2];
    }
	  
    if (strncmp(temp,"size=",5) == 0)
    {
      std::istringstream sizeLine(&temp[5]);
      sizeLine >> size[0] >> size[1] >> size[2];	    
    }
  }
  
  vffFileStream.close();

  string nrrdType="short";
  if (nbits == 8) nrrdType="uchar";

  nrrdHeader << "NRRD0001" << endl
	     << "type: "<< nrrdType << endl
	     << "dimension: " << dim << endl
	     << "sizes:";
  for (i = 0; i < dim; i++)
  {
    nrrdHeader << " " << size[i];
  }
  nrrdHeader << endl
             << "spacings:";
  for (i = 0; i < dim; i++)
  {
    nrrdHeader << " " << spacing[i];
  }
  nrrdHeader << endl
             << "axis mins:";
  for (i = 0; i < dim; i++)
  {
    nrrdHeader << " " << "0.0";
  }
  nrrdHeader << endl
             << "axis maxs:";
  for (i = 0; i < dim; i++)
  {
    nrrdHeader << " " << spacing[i]*(size[i]-1);
  }
  nrrdHeader << endl
             << "centerings:";
  for (i = 0; i < dim; i++)
  {
    nrrdHeader << " node";
  }
  nrrdHeader << endl
             << "data file: " << in << endl
             << "endian: big" << endl
             << "encoding: raw" << endl
             << "line skip: " << nrrdSkip << endl;

  nrrdHeader.close();
  return 0;
}
