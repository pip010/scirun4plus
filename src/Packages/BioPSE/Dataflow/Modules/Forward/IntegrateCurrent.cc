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
 *  IntegrateCurrent: Compute current through a surface
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2002
 *
 */

#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/Field.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>

#include <iostream>
#include <sstream>

namespace BioPSE {

using namespace SCIRun;

class IntegrateCurrent : public Module {
  private:
    GuiDouble current_;
  public:
    IntegrateCurrent(GuiContext *context);
    virtual ~IntegrateCurrent() {}
    virtual void execute();
};

DECLARE_MAKER(IntegrateCurrent)


IntegrateCurrent::IntegrateCurrent(GuiContext *context)
  : Module("IntegrateCurrent", context, Filter, "Forward", "BioPSE"),
    current_(context->subVar("current"))
{
}

void
IntegrateCurrent::execute()
{
  FieldHandle efieldH, sigmasH, trisurfH;
  get_input_handle("TetMesh EField", efieldH, true);
  get_input_handle("TetMesh Sigmas", sigmasH, true);
  get_input_handle("TriSurf", trisurfH, true);

  if (efieldH.get_rep() == 0)
  {
    error("Efield input field is empty");
    return;
  }

  if (sigmasH.get_rep() == 0)
  {
    error("Sigmas input field is empty");
    return;
  }

  if (trisurfH.get_rep() == 0)
  {
    error("Trisurf input field is empty");
    return;
  }
  
  FieldInformation tsfi(trisurfH);
  
  if (!(tsfi.is_trisurfmesh()))
  {
    error("Trisurf Field is not a Trisurf mesh.");
    return;
  }
  
  if (efieldH->mesh().get_rep() != sigmasH->mesh().get_rep()) 
  {
    error("EField and Sigma Field need to have the same mesh.");
    return;
  }

  FieldInformation efi(efieldH);

  if (!(efi.is_vector()))
  {
    error("EField is not a vector field.");
    return;
  }
  
  VField* efield = efieldH->vfield();
  VMesh*  emesh  = efieldH->vmesh();
  VMesh*  tris   = trisurfH->vmesh();
  VField* sigmas = sigmasH->vfield();

  std::vector<std::pair<std::string, Tensor> > conds;
  if (!sigmasH->get_property("conductivity_table", conds)) 
  {
    error("No conductivity_table found in Sigmas.");
    return;
  }

  emesh->synchronize(Mesh::ELEM_LOCATE_E);

  // For each face in tris, find its area, centroid, and normal
  // for that centroid, look up its sigma and efield in the tetvol fields
  // compute (sigma * efield * area) and dot it with the face normal
  // sum those up for all tris

  VMesh::Elem::iterator fi, fe;
  tris->begin(fi);
  tris->end(fe);
  double current=0;
  VMesh::Node::array_type nodes;

  VMesh::Elem::size_type nfaces;
  tris->size(nfaces);
  size_t progress = 0;
  size_t progress_max = (size_t)nfaces;
  update_progress(0.0);

  double total_area=0;
  while (fi != fe) 
  {
    Point center;
    tris->get_center(center, *fi);
    double area = tris->get_size(*fi);
    total_area += area;
    tris->get_nodes(nodes, *fi);
    Point p0, p1, p2;
    tris->get_center(p0, nodes[0]);
    tris->get_center(p1, nodes[1]);
    tris->get_center(p2, nodes[2]);
    
    Vector normal(Cross(p2-p1,p2-p0));
    normal.normalize();
    
    VMesh::Elem::index_type tet;
    if (!emesh->locate(tet, center)) 
    {
      error("Trisurf centroid was not located in tetvolmesh.");
      return;
    }
    Vector e;
    efield->get_value(e,tet);
    int sigma_idx;
    sigmas->get_value(sigma_idx,tet);
    Tensor s(conds[sigma_idx].second);
    
    // compute sigma * e
    Vector c(s.mat_[0][0]*e.x()+s.mat_[0][1]*e.y()+s.mat_[0][2]*e.z(),
	     s.mat_[1][0]*e.x()+s.mat_[1][1]*e.y()+s.mat_[1][2]*e.z(),
	     s.mat_[2][0]*e.x()+s.mat_[2][1]*e.y()+s.mat_[2][2]*e.z());
    current += Abs(Dot(c,normal)) * area;

    if ((progress++ & 0xff) == 0) 
    { 
      update_progress(progress, progress_max); 
    }

    ++fi;
  }
  
  update_progress(1.0);
  current_.set(current);
}

} // End namespace BioPSE
