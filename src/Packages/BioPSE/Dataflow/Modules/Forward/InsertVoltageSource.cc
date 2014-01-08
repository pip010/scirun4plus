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
 *  InsertVoltageSource: Insert a voltage source
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2002
 *
 */


#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>
#include <vector>

namespace BioPSE {

using namespace SCIRun;

class InsertVoltageSource : public Module {
  private:
    GuiInt outside_;
    GuiInt groundfirst_;

  public:
    InsertVoltageSource(GuiContext *context);
    virtual ~InsertVoltageSource() {}
    virtual void execute();
};

DECLARE_MAKER(InsertVoltageSource)

InsertVoltageSource::InsertVoltageSource(GuiContext *context)
  : Module("InsertVoltageSource", context, Filter, "Forward", "BioPSE"),
    outside_(context->subVar("outside")),
    groundfirst_(context->subVar("groundfirst"))
{
}

void
InsertVoltageSource::execute()
{
  FieldHandle imeshH;
  get_input_handle("FEMesh", imeshH, true);

  FieldInformation fi(imeshH);
  if (fi.is_pointcloudmesh())
  {
    error("FEMesh is a point cloud mesh, the FE mesh needs to have elements");
    return;
  }
  
  FieldHandle isourceH;
  get_input_handle("VoltageSource", isourceH, true);

  FieldInformation fis(isourceH);
  
  if (fis.is_nodata())
  {
    error("VoltageSource needs to contain data");
    return;
  }

  int groundfirst = groundfirst_.get();
  std::vector<Point> sources;
  std::vector<double> vals;
  
  if (groundfirst) 
  {
    Point pnt;
    if (isourceH->vmesh()->num_nodes() == 0)
    {
      error("VoltageSource field does not have any nodes");
      return;
    }
    isourceH->vmesh()->get_center(pnt,VMesh::Node::index_type(0));

    sources.push_back(pnt);
    vals.push_back(0);
  } 
  else 
  {
  
    VField* field = isourceH->vfield();
    VMesh*  mesh  = isourceH->vmesh();
    VField::size_type num_values = field->num_values();

    Point pnt;
    double val;
    
    for (VField::index_type idx=0; idx<num_values; idx++)
    {
      if (field->basis_order() == 0)
      {
        mesh->get_center(pnt,VMesh::Elem::index_type(idx));
      }
      else
      {
        mesh->get_center(pnt,VMesh::Node::index_type(idx));
      }
      
      field->get_value(val,idx);
      sources.push_back(pnt);
      vals.push_back(val);
    }
  }

  std::vector<std::pair<int, double> > dirichlet;
  imeshH->get_property("dirichlet", dirichlet);

  std::vector<std::pair<std::string, Tensor> > conds;
  imeshH->get_property("conductivity_table", conds);

  // get our own local copy of the Field and mesh
  imeshH.detach();  

  int outside = outside_.get();

  VMesh* vmesh = imeshH->vmesh();
  VMesh::Node::size_type nnodes = vmesh->num_nodes();
  VMesh::Elem::size_type nelems = vmesh->num_elems();
  vmesh->synchronize(Mesh::LOCATE_E);

  VMesh::Node::array_type nbrs;
  Array1<std::vector<std::pair<double, double> > > closest(nnodes); 
                                     // store the dist/val
                                     // to source nodes
  Array1<int> have_some(nnodes);
  have_some.initialize(0);
  Array1<VMesh::Node::index_type> bc_nodes;

  for (size_t di=0; di<dirichlet.size(); di++) 
  {
    int didx=dirichlet[di].first;
    // so other BCs can't over-write this one
    have_some[didx]=1;
    closest[didx].push_back(std::pair<double, double>(0,dirichlet[di].second));
  }
    
  // for each surface data_at position/value...
  for (size_t s=0; s<sources.size(); s++) 
  {
    Point pt=sources[s];
    double val=vals[s];

    // find the tet nodes (nbrs) that it's closest to
    VMesh::Elem::index_type cidx;
    
    if (vmesh->locate(cidx, pt)) 
    {
      vmesh->get_nodes(nbrs, cidx);
    } 
    else if (outside) 
    {
      nbrs.resize(1);
      vmesh->locate(nbrs[0], pt);
    } 
    else
    {
      continue;
    }
    
    // for each nbr, see if this node is the closest of the bc nodes checked
    //   so far -- if so, store it
    double dmin=-1.0;
    VMesh::Node::index_type nmin;

    for (size_t i=0; i<nbrs.size(); i++) 
    {
      Point nbr_pt;
      VMesh::Node::index_type nidx = nbrs[i];
      vmesh->get_center(nbr_pt, nidx);
      double d = (pt - nbr_pt).length();
      if (i==0 || d<dmin)
      {
        nmin=nbrs[i];
        dmin=d;
      }
    }
    if (dmin != -1.0 && !have_some[nmin]) 
    {
      std::pair<double, double> p(dmin, val);
      have_some[nmin]=1;
      bc_nodes.add(nmin);
      closest[nmin].push_back(p);
    }
  }

  for (size_type i=0; i<bc_nodes.size(); i++)
  {
    double val=0;
    int nsrcs=closest[bc_nodes[i]].size();
    for (int j=0; j<nsrcs; j++)
      val+=closest[bc_nodes[i]][j].second/nsrcs;
    dirichlet.push_back(std::pair<int, double>((int)bc_nodes[i], val));
  }
  
  imeshH->set_property("dirichlet", dirichlet, false);
  imeshH->set_property("conductivity_table", conds, false);

  send_output_handle("FEMesh", imeshH);
}


} // End namespace BioPSE
