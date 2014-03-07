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
 *  TriSurfFieldToVtk.cc
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

#include <sci_deprecated.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace SCIRun;

int
main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: "<<argv[0]<<" scirun_inner_trisurf scirun_outer_trisurf tetgen_basename\n";
    return 0;
  }

  FieldHandle inner_surf;
  PiostreamPtr stream1=auto_istream(argv[1]);
  if (!stream1) 
  {
    std::cerr << "Couldn't open file "<<argv[1]<<".  Exiting...\n";
    exit(0);
  }
  Pio(*stream1, inner_surf);
  
  if (!inner_surf.get_rep()) 
  {
    std::cerr << "Error reading surface from file "<<argv[1]<<".  Exiting...\n";
    exit(0);
  }
  
  FieldInformation ifi(inner_surf);
  if (ifi.is_trisurfmesh()) 
  {
    std::cerr << "Error -- input field wasn't a TriSurfField\n";
    exit(0);
  }
  
  VMesh *inner = inner_surf->vmesh();

  FieldHandle outer_surf;
  PiostreamPtr stream2=auto_istream(argv[2]);
  
  if (!stream2) 
  {
    std::cerr << "Couldn't open file "<<argv[2]<<".  Exiting...\n";
    exit(0);
  }
  
  Pio(*stream2, outer_surf);
  if (!outer_surf.get_rep()) 
  {
    std::cerr << "Error reading surface from file "<<argv[2]<<".  Exiting...\n";
    exit(0);
  }
  
  FieldInformation ofi(outer_surf);
  if (ofi.is_trisurfmesh()) 
  {
    std::cerr << "Error -- input field wasn't a TriSurfField\n";
    exit(0);
  }
  
  VMesh *outer = outer_surf->vmesh();

  // TODO: patch code to modern standards
  // No more fprintf, sprintf fopen, etc
  // These don't manage memory properly
  
  char filename[1000];
  sprintf(filename, "%s.poly", argv[3]);
  FILE *fout = fopen(filename, "wt");
  if (!fout) 
  {
    std::cerr << "Error - can't open input file "<<filename<<"\n";
    exit(0);
  }

  VMesh::Node::iterator niter, niter_end;
  VMesh::Node::size_type ninner;
  inner->size(ninner);
  VMesh::Face::size_type finner;
  inner->size(finner);
  VMesh::Node::size_type nouter;
  outer->size(nouter);
  VMesh::Face::size_type fouter;
  outer->size(fouter);
  
  fprintf(fout, "# Number of nodes, pts/tri, no holes, no boundary markers\n");
  fprintf(fout, "%ld 3 0 0\n", ninner+nouter);
  fprintf(fout, "\n# Inner surface points\n");

  int i;
  Point p;
  inner->begin(niter);
  inner->end(niter_end);
  Point mid_inner;
  for (i=0; niter != niter_end; i++, ++niter) 
  {
    inner->get_center(p, *niter);
    mid_inner += p.asVector();
    fprintf(fout, "%d %lf %lf %lf\n", i+1, p.x(), p.y(), p.z());
  }
  mid_inner /= ninner;

  fprintf(fout, "\n# Outer surface points\n");
  outer->begin(niter);
  outer->end(niter_end);
  for (; niter != niter_end; i++, ++niter) 
  {
    outer->get_center(p, *niter);
    fprintf(fout, "%d %lf %lf %lf\n", i+1, p.x(), p.y(), p.z());
  }

  fprintf(fout, "\n# Number of faces, no boundary markers\n");
  fprintf(fout, "%d 0\n\n", finner+fouter);
  Point mid_outer = AffineCombination(mid_inner, 0.5, p, 0.5);

  VMesh::Face::iterator fiter, fiter_end;
  VMesh::Node::array_type fac_nodes(3);
  inner->begin(fiter);
  inner->end(fiter_end);
  
  fprintf(fout, "# Inner faces\n");
  for (i=0; fiter != fiter_end; i++, ++fiter) 
  {
    inner->get_nodes(fac_nodes, *fiter);
    int i1, i2, i3;
    i1=fac_nodes[0]+1;
    i2=fac_nodes[1]+1;
    i3=fac_nodes[2]+1;
    fprintf(fout, "1\n3 %d %d %d\n", i1, i2, i3);
  }
  outer->begin(fiter);
  outer->end(fiter_end);
  
  fprintf(fout, "\n# Outer faces\n");
  for (; fiter != fiter_end; i++, ++fiter) 
  {
    outer->get_nodes(fac_nodes, *fiter);
    int i1, i2, i3;
    i1=fac_nodes[0]+1+ninner;
    i2=fac_nodes[1]+1+ninner;
    i3=fac_nodes[2]+1+ninner;
    fprintf(fout, "1\n3 %d %d %d\n", i1, i2, i3);
  }
  
  fprintf(fout, "\n# No holes\n0\n\n# Two regions\n2\n\n");
  fprintf(fout, "1 %lf %lf %lf 1\n", mid_inner.x(), mid_inner.y(), mid_inner.z());
  fprintf(fout, "2 %lf %lf %lf 2\n", mid_outer.x(), mid_outer.y(), mid_outer.z());
  
  fclose(fout);
  return 0;
}    
