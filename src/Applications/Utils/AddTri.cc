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
//    File   : AddTri.cc
//    Author : Martin Cole
//    Date   : Thu Mar  3 20:13:23 2005

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Mesh.h>
#include <Core/Persistent/Pstreams.h>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace SCIRun;

int
main(int /*argc*/, char **argv) 
{

  FieldHandle handle;
  PiostreamPtr stream = auto_istream(argv[1]);
  if (!stream) 
  {
    std::cerr << "Couldn't open file " << argv[1] << ".  Exiting..." << std::endl;
    exit(1);
  }
  Pio(*stream, handle);
  if (!handle.get_rep()) 
  {
    std::cerr << "Error reading surface from file " << argv[1] << ".  Exiting..." << std::endl;
    exit(2);
  }
  
  VMesh *tsm = handle->vmesh();
  if (! tsm->is_trisurfmesh()) 
  { 
    std::cerr << "Error: not a TriSurf." << std::endl; return 99;
  }


  VMesh::Node::array_type array(3);
  array[0] = VMesh::Node::index_type(atoi(argv[2]));
  array[1] = VMesh::Node::index_type(atoi(argv[3]));
  array[2] = VMesh::Node::index_type(atoi(argv[4]));
  tsm->add_elem(array);
  
  BinaryPiostream out_stream(argv[5], Piostream::Write);
  Pio(out_stream, handle);
  return 0;  
}    
