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
//    File   : JoinFields.cc
//    Author : David Weinstein
//    Date   : Mon Jan  5 14:22:37 MST 2009

#include <Core/Algorithms/Fields/MergeFields/JoinFields.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Persistent/Pstreams.h>

#include <iostream>
#include <fstream>

#include <stdlib.h>

using std::cerr;
using std::ifstream;
using std::endl;

using namespace SCIRun;

int
main(int argc, char *argv[])
{
  if (argc < 3 || (argc == 3 && std::string(argv[1])==std::string("-indexdata")))
  {
    std::cerr << "Usage: "<<argv[0]<<" [-indexdata] outputFieldFile inputFieldFile1 inputFieldFile2 ... inputFieldFileN\n";
    exit(-1);
  }
  bool indexdata=false;
  if (std::string(argv[1]) == std::string("-indexdata")) indexdata=true;

  SCIRunAlgo::JoinFieldsAlgo algo;
  std::vector<FieldHandle> ffields;

  int j;
  if (indexdata) j=3; else j=2;
  for (; j<argc; j++)
  {
    FieldHandle fldH;
    PiostreamPtr fldStream = auto_istream(argv[j]);
    if (!fldStream)
    {
      std::cerr << "Couldn't open file " << argv[j] << ".  Exiting..." << std::endl;
      exit(-1);
    }
    Pio(*fldStream, fldH);
    if (!fldH.get_rep())
    {
      std::cerr << "Error reading field from file " << argv[1] << ".  Exiting..."
                << std::endl;
      exit(-1);
    }
    FieldInformation fi(fldH);
    if (indexdata)
    {
      fi.set_basis_type(0);
      fi.set_data_type("unsigned_char");
      FieldHandle f2=CreateField(fi, fldH->mesh());
      f2->vfield()->set_all_values(j-3);
      ffields.push_back(f2);
    }
    else
    {
      ffields.push_back(fldH);
    }
  }
  FieldHandle ofield;
  algo.set_bool("merge_nodes", false);
  algo.set_bool("merge_elems", false);
  algo.set_bool("match_node_values", false);

  if (!(algo.run(ffields, ofield)))
  {
    cerr << "Error running JoinFieldsAlgo\n";
    exit(-1);
  }
  char *outname;
  if (indexdata) outname=argv[2]; else outname=argv[1];
  BinaryPiostream outputStream(outname, Piostream::Write);
  Pio(outputStream, ofield);
  return 0;
}
