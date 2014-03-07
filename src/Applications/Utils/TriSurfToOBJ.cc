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
 *  TriSurfToOBJ.cc
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   May 2001
 *
 */

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Util/StringUtil.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace SCIRun;

int
main(int argc, char **argv) 
{
  if (argc != 7) {
    std::cerr << "Usage: "<<argv[0]<<" scirun_field output_basename color_r color_g color_b Ns\n";
    return 0;
  }

  std::string argv3(argv[3]);
  std::string argv4(argv[4]);
  std::string argv5(argv[5]);
  std::string argv6(argv[6]);

  double red, green, blue, Ns;

  from_string(argv3,red);
  from_string(argv4,green);
  from_string(argv5,blue);
  from_string(argv6,Ns);

  FieldHandle handle;
  PiostreamPtr stream=auto_istream(argv[1]);
  if (!stream) 
  {
    std::cerr << "Couldn't open file "<<argv[1]<<".  Exiting...\n";
    exit(0);
  }
  Pio(*stream, handle);
  
  if (!handle.get_rep()) 
  {
    std::cerr << "Error reading surface from file "<<argv[1]<<".  Exiting...\n";
    exit(0);
  }
  
  FieldInformation fi(handle);
  if (fi.is_trisurfmesh()) 
  {
    std::cerr << "Error -- input field didn't have a TriSurfMesh\n";
    exit(0);
  }

  VMesh *tsm = handle->vmesh();

  // TODO: patch code to modern standards
  // No more fprintf, sprintf fopen, etc
  // These don't manage memory properly
  
  char objname[512];
  sprintf(objname, "%s.obj", argv[2]);
  char mtlname[512];
  sprintf(mtlname, "%s.mtl", argv[2]);
  FILE *fobj = fopen(objname, "wt");
  FILE *fmtl = fopen(mtlname, "wt");
  if (!fobj || !fmtl) 
  {
    std::cerr << "Error - can't open output file "<<objname<<" or "<<mtlname<<"\n";
    exit(0);
  }

  tsm->synchronize(Mesh::NORMALS_E);
  VMesh::Node::iterator niter; 
  VMesh::Node::iterator niter_end; 
  VMesh::Node::size_type nsize; 
  tsm->begin(niter);
  tsm->end(niter_end);
  tsm->size(nsize);
  while(niter != niter_end)
  {
    Point p;
    tsm->get_center(p, *niter);
    fprintf(fobj, "v %lf %lf %lf\n", p.x(), p.y(), p.z());
    Vector n;
    tsm->get_normal(n, *niter);
    fprintf(fobj, "vn %lf %lf %lf\n", n.x(), n.y(), n.z());
    ++niter;
  }
  fprintf(fobj, "usemtl Default\n");  
  
  VMesh::Face::size_type fsize; 
  VMesh::Face::iterator fiter; 
  VMesh::Face::iterator fiter_end; 
  VMesh::Node::array_type fac_nodes(3);
  
  tsm->size(fsize);
  tsm->begin(fiter);
  tsm->end(fiter_end);
  while(fiter != fiter_end)
  {
    tsm->get_nodes(fac_nodes, *fiter);
    fprintf(fobj, "f  %d//%d %d//%d %d//%d\n", 
	    (int)fac_nodes[0]+1, (int)fac_nodes[0]+1,
	    (int)fac_nodes[1]+1, (int)fac_nodes[1]+1,
	    (int)fac_nodes[2]+1, (int)fac_nodes[2]+1);
	    
    ++fiter;
  }
  
  fclose(fobj);
  fprintf(fmtl,"newmtl Default\n");
  fprintf(fmtl," \tKa 0 0 0\n");
  fprintf(fmtl," \tKd %lf %lf %lf\n", red, green, blue);
  fprintf(fmtl," \tKs 1 1 1\n");
  fprintf(fmtl," \tillum 2\n");
  fprintf(fmtl," \tNs %lf\n", Ns);
  fclose(fmtl);
  return 0;  
}    
