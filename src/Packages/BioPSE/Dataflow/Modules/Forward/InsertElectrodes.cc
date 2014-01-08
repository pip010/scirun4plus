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
 *  InsertElectrodes: Insert an electrode into a finite element mesh
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2002
 *
 */

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

#include <Core/Geometry/Plane.h>
#include <Core/Geometry/Transform.h>
#include <Core/Util/StringUtil.h>
#include <Core/Containers/Array1.h>
#include <iostream>
#include <sstream>

namespace BioPSE {
using namespace SCIRun;

class InsertElectrodes : public Module {

  bool inside4_p(VMesh* tet_mesh, 
           const VMesh::Node::array_type &nodes,
           const Point &p);
  void insertNodesIntoTetMesh(std::vector<std::pair<int, double> > &dirichlet,
                              VField *tet_fld,
                              Array1<Point> &new_points,
                              Array1<double> &new_values);
  void insertContourIntoTetMesh(std::vector<std::pair<int, double> > &dirichlet,
                                VMesh* mesh,
                                const std::vector<Point> &inner,
                                const std::vector<Point> &outer,
                                const std::string &active_side,
                                double voltage,
                                VMesh* electrodeElements);
  bool point_in_loop(const Point &pt, 
                     const std::vector<Point> &contour, 
                     const Plane &pl);

public:
  InsertElectrodes(GuiContext *context);
  virtual ~InsertElectrodes() {}
  virtual void execute();
};


DECLARE_MAKER(InsertElectrodes)


InsertElectrodes::InsertElectrodes(GuiContext *context)
    : Module("InsertElectrodes", context, Filter, "Forward", "BioPSE")
{
}


bool
InsertElectrodes::point_in_loop(const Point &pt,
				const std::vector<Point> &contour,
				const Plane &pl)
{
  Vector v1, v2, v3;
  unsigned int i;
  double angle, sum=0;
  if (contour.size() == 0) return false;
  for (i=0; i<contour.size()-1; i++)
  {
    v1=pt-contour[i];
    v2=pt-contour[i+1];
    v1.normalize();
    v2.normalize();
    angle=acos(Dot(v1,v2));
    v3=Cross(v1, v2);
    if (pl.eval_point(pt+v3) < 0) angle=-angle;
    sum += angle;
  }
  sum/=(2.*M_PI);
  sum=fabs(sum);
  sum -= (((int)sum)/2)*2;
  // sum is now between 0 and 1.99999
  if (sum > 0.8 && sum < 1.2) return true;
  else if (sum < 0.2 || sum > 1.8) return false;
  else
  {
    error("Point in loop was: " + to_string(sum) + ".");
    return false;
  }
}

bool
InsertElectrodes::inside4_p(VMesh* tet_mesh, 
			    const VMesh::Node::array_type &nodes,
			    const Point &p) 
{
  Point p0, p1, p2, p3;
  tet_mesh->get_center(p0, nodes[0]);
  tet_mesh->get_center(p1, nodes[1]);
  tet_mesh->get_center(p2, nodes[2]);
  tet_mesh->get_center(p3, nodes[3]);
  const double x0 = p0.x();
  const double y0 = p0.y();
  const double z0 = p0.z();
  const double x1 = p1.x();
  const double y1 = p1.y();
  const double z1 = p1.z();
  const double x2 = p2.x();
  const double y2 = p2.y();
  const double z2 = p2.z();
  const double x3 = p3.x();
  const double y3 = p3.y();
  const double z3 = p3.z();

  const double a0 = + x1*(y2*z3-y3*z2) + x2*(y3*z1-y1*z3) + x3*(y1*z2-y2*z1);
  const double a1 = - x2*(y3*z0-y0*z3) - x3*(y0*z2-y2*z0) - x0*(y2*z3-y3*z2);
  const double a2 = + x3*(y0*z1-y1*z0) + x0*(y1*z3-y3*z1) + x1*(y3*z0-y0*z3);
  const double a3 = - x0*(y1*z2-y2*z1) - x1*(y2*z0-y0*z2) - x2*(y0*z1-y1*z0);
  const double iV6 = 1.0 / (a0+a1+a2+a3);

  const double b0 = - (y2*z3-y3*z2) - (y3*z1-y1*z3) - (y1*z2-y2*z1);
  const double c0 = + (x2*z3-x3*z2) + (x3*z1-x1*z3) + (x1*z2-x2*z1);
  const double d0 = - (x2*y3-x3*y2) - (x3*y1-x1*y3) - (x1*y2-x2*y1);
  const double s0 = iV6 * (a0 + b0*p.x() + c0*p.y() + d0*p.z());
  if (s0 < -1.e-6)
    return false;

  const double b1 = + (y3*z0-y0*z3) + (y0*z2-y2*z0) + (y2*z3-y3*z2);
  const double c1 = - (x3*z0-x0*z3) - (x0*z2-x2*z0) - (x2*z3-x3*z2);
  const double d1 = + (x3*y0-x0*y3) + (x0*y2-x2*y0) + (x2*y3-x3*y2);
  const double s1 = iV6 * (a1 + b1*p.x() + c1*p.y() + d1*p.z());
  if (s1 < -1.e-6)
    return false;

  const double b2 = - (y0*z1-y1*z0) - (y1*z3-y3*z1) - (y3*z0-y0*z3);
  const double c2 = + (x0*z1-x1*z0) + (x1*z3-x3*z1) + (x3*z0-x0*z3);
  const double d2 = - (x0*y1-x1*y0) - (x1*y3-x3*y1) - (x3*y0-x0*y3);
  const double s2 = iV6 * (a2 + b2*p.x() + c2*p.y() + d2*p.z());
  if (s2 < -1.e-6)
    return false;

  const double b3 = +(y1*z2-y2*z1) + (y2*z0-y0*z2) + (y0*z1-y1*z0);
  const double c3 = -(x1*z2-x2*z1) - (x2*z0-x0*z2) - (x0*z1-x1*z0);
  const double d3 = +(x1*y2-x2*y1) + (x2*y0-x0*y2) + (x0*y1-x1*y0);
  const double s3 = iV6 * (a3 + b3*p.x() + c3*p.y() + d3*p.z());
  if (s3 < -1.e-6)
    return false;

  return true;
}

void
InsertElectrodes::insertNodesIntoTetMesh(std::vector<std::pair<int, double> > &dirichlet,
					 VField* tet_fld,
					 Array1<Point> &new_points,
					 Array1<double> &new_values) 
{
  VMesh* tet_mesh=tet_fld->vmesh();
  int np=new_points.size();
  tet_mesh->synchronize(Mesh::LOCATE_E);

  Array1<VMesh::Elem::index_type> split_elems;
  Array1<int> new_fdata;
  Array1<Array1<int> > split_elem_nodes;
  VMesh::Elem::index_type cidx;
  Array1<VMesh::Node::index_type> node_map(np); // where do new nodes get
                                                     //  added to TVMesh
  for (int i=0; i<np; i++) 
  {
    if (tet_mesh->locate(cidx, new_points[i])) 
    {
      node_map[i]=tet_mesh->add_point(new_points[i]);
      dirichlet.push_back(std::pair<int,double>((int)node_map[i], new_values[i]));
      int already_split_elem=0;
      for (int j=0; j<split_elems.size() && !already_split_elem; j++) 
      {
        if (cidx == split_elems[j]) 
        {
          already_split_elem=1;
          split_elem_nodes[j].add(i);
        }
      }
      if (!already_split_elem) 
      {
        split_elems.add(cidx);
        int sz=split_elem_nodes.size();
        split_elem_nodes.resize(sz+1);
        split_elem_nodes[sz].add(i);
      }
    } // else outside mesh, so ignore it
  }

  // split_elems now contains the list of all the elements to be split
  // and split_elem_nodes contains the list of new nodes for each element
  // build up a list of new elems -- we'll just keep modifying this list
  Array1<VMesh::Node::array_type> new_tets;
  Array1<int> new_vals;
  VMesh::Node::array_type nodes;
  VMesh::Node::index_type tmp;
  VMesh::Node::index_type nidx;
  for (int e=0; e<split_elems.size(); e++) 
  {
    tet_mesh->get_nodes(nodes, split_elems[e]);
    int val;
    tet_fld->get_value(val,split_elems[e]);
    nidx=node_map[split_elem_nodes[e][0]];
    for (int ii=0; ii<4; ii++) 
    {
      tmp=nodes[ii]; 
      nodes[ii]=nidx;
      new_tets.add(nodes);
      new_vals.add(val);
      nodes[ii]=tmp;
    }
    // find which tet from new_tets this node is in -- then split it
    for (int n=1; n<split_elem_nodes[e].size(); n++) 
    {
      Point p(new_points[split_elem_nodes[e][n]]);
      nidx=node_map[split_elem_nodes[e][n]];
      bool located=false;
      for (int k=0; k<new_tets.size() && !located; k++) 
      {
        if (inside4_p(tet_mesh, new_tets[k], p)) 
        {
          located=true;
          nodes=new_tets[k];
          for (int ii=0; ii<4; ii++) 
          {
            tmp=nodes[ii];
            nodes[ii]=nidx;
            if (ii==0) 
            {
              new_tets[k]=nodes;
            } 
            else 
            {
              new_tets.add(nodes);
              new_vals.add(new_vals[k]);
            }
            nodes[ii]=tmp;
          }
        }
      }
      if (!located) 
      {
        std::cerr << "ERROR -- COULDN'T FIND NODE!!!!\n";
      }
    }
  }

  size_t e;
  for (e=0; e<split_elems.size(); e++) 
  {
    tet_mesh->set_nodes(new_tets[e], split_elems[e]);
  }

  for (; e<new_tets.size(); e++) 
  {
    tet_mesh->add_elem(new_tets[e]);
  }
  
  tet_fld->resize_values();
  
  for (e=0; e<split_elems.size(); e++) 
  {
    tet_fld->set_value(new_vals[e], split_elems[e]);
  }
  
  for (; e<new_tets.size(); e++) 
  {
    tet_fld->set_value(new_vals[e], e);
  }
}
			   

void
InsertElectrodes::insertContourIntoTetMesh(std::vector<std::pair<int, double> > &dirichlet,
					   VMesh* tet_mesh,
					   const std::vector<Point> &inner,
					   const std::vector<Point> &outer,
					   const std::string &active_side,
					   double voltage,
					   VMesh* electrodeElements)
{
  Plane electrode_plane;
  if (Cross(inner[1]-inner[0], inner[1]-inner[2]).length2()>1.e-10)
  {
    electrode_plane=Plane(inner[1], inner[0], inner[2]);
  }
  else if (Cross(inner[1]-inner[0], inner[1]-inner[3]).length2()>1.e-10)
  {
    electrode_plane=Plane(inner[1], inner[0], inner[3]);
  }
  else
  {
    error("First four inner nodes of the elec are colinear.");
    return;
  }
  if (active_side == "back") electrode_plane.flip();

  VMesh::Elem::size_type ncells;
  tet_mesh->size(ncells);
  VMesh::Node::size_type nnodes;
  tet_mesh->size(nnodes);

  VMesh::Node::array_type tet_nodes;
  std::vector<int> is_electrode_node(nnodes, -1); // -1:unknown, 0:no, 1:outer, 2:inner
  std::vector<VMesh::Node::index_type> electrode_nodes;
  std::vector<VMesh::Node::index_type> electrode_node_split_idx(nnodes);

  VMesh::Cell::iterator tet_cell_idx, tet_cell_idx_end;
  tet_mesh->begin(tet_cell_idx); tet_mesh->end(tet_cell_idx_end);
  while (tet_cell_idx != tet_cell_idx_end)
  {
    tet_mesh->get_nodes(tet_nodes, *tet_cell_idx);
    int all_nodes_above_plane=1, all_nodes_below_plane=1;
    std::vector<Point> tet_node_pts(4);
    double signed_distance[4];
    std::vector<VMesh::Node::index_type> projections;
    for (int ii=0; ii<4; ii++)
    {
      tet_mesh->get_center(tet_node_pts[ii], tet_nodes[ii]);
      signed_distance[ii]=electrode_plane.eval_point(tet_node_pts[ii]);
      if (signed_distance[ii] > 0)
      {
        all_nodes_above_plane = 0;
      }
      else
      {
        all_nodes_below_plane = 0;
        projections.push_back(tet_nodes[ii]);
      }
    }
    if (!all_nodes_below_plane && 
        !all_nodes_above_plane)
    { // cell was intersected by plane...
      
      // check to see if projected nodes are inside the electrode
      unsigned int idx;
      for (idx=0; idx<projections.size(); idx++)
      {
        Point p;
        VMesh::Node::index_type nidx(projections[idx]);
        tet_mesh->get_center(p, nidx);
        p = electrode_plane.project(p);
        if (point_in_loop(p, inner, electrode_plane))
        {
          if (is_electrode_node[nidx] == -1)
          {
            is_electrode_node[nidx] = 2;
            electrode_nodes.push_back(nidx);
          }
        }
      }
      for (idx=0; idx<projections.size(); idx++)
      {
        Point p;
        VMesh::Node::index_type nidx(projections[idx]);
        tet_mesh->get_center(p, nidx);
        p = electrode_plane.project(p);
        if (point_in_loop(p, outer, electrode_plane))
        {
          if (is_electrode_node[nidx] == -1)
          {
            is_electrode_node[nidx] = 1;
            electrode_nodes.push_back(nidx);
          }
        }
      }
    }
    ++tet_cell_idx;
  }

  // project electrode nodes to plane, split them, and set Dirichlet
  unsigned int i;
  for (i=0; i<electrode_nodes.size(); i++)
  {
    Point pt;
    VMesh::Node::index_type ni = electrode_nodes[i];
    tet_mesh->get_center(pt, ni);
    Point proj_pt = electrode_plane.project(pt);
    tet_mesh->set_point(proj_pt, ni);
    electrode_node_split_idx[ni] = tet_mesh->add_point(proj_pt);
    if (is_electrode_node[electrode_nodes[i]] == 2)
    {
      dirichlet.push_back(std::pair<int,double>((int)electrode_node_split_idx[ni], 
					   voltage));
      if (active_side == "both") 
      {
        dirichlet.push_back(std::pair<int,double>((int)ni, voltage));
      }
    }
  }

  // set elements above the plane to point to new split node
  VMesh::Elem::iterator ci, ce;
  tet_mesh->begin(ci); tet_mesh->end(ce);
  while (ci != ce)
  {
    tet_mesh->get_nodes(tet_nodes, *ci);
    int have_any =0 ;
    for (i=0; i<4; i++)
    {
      if (is_electrode_node[tet_nodes[i]] > 0) have_any = 1;
    }
    if (have_any)
    {
      Point centroid;
      tet_mesh->get_center(centroid, *ci);
      double centroid_dist = electrode_plane.eval_point(centroid);
      if (centroid_dist>0)
      {
        if (electrodeElements)
        {
          Point pts[4];
          for (i=0; i<4; i++) 
          {
            tet_mesh->get_point(pts[i], tet_nodes[i]);
            electrodeElements->add_point(pts[i]);
          }
          electrodeElements->add_elem(tet_nodes);
        }
        int remap=0;
        for (i=0; i<4; i++)
        {
          if (is_electrode_node[tet_nodes[i]] > 0)
          {
            remap++;
            tet_nodes[i] = electrode_node_split_idx[tet_nodes[i]];
          }
        }
        tet_mesh->set_nodes(tet_nodes, *ci);
      }
    }
    ++ci;
  }
}



void
InsertElectrodes::execute()
{
  FieldHandle imeshH;
  get_input_handle("TetMesh", imeshH,true);

  FieldInformation fi(imeshH);
  if (!(fi.is_tetvolmesh() && fi.is_constantdata()))
  {
    error("Input TetMesh is not a tetvolmesh with data on elements");
    return;
  }


  std::vector<std::pair<int, double> > dirichlet_nodes;
  imeshH->vfield()->get_property("dirichlet", dirichlet_nodes);
  std::vector<std::pair<std::string, Tensor> > conds;
  imeshH->vfield()->get_property("conductivity_table", conds);

  std::vector<FieldHandle> ielecH;
  get_dynamic_input_handles("Electrodes",ielecH,true);

  // get our own local copy of the Field and mesh
  imeshH.detach();
  imeshH->mesh_detach();

  VField* field = imeshH->vfield();
  VMesh*  mesh  = imeshH->vmesh();

  MeshHandle elecElemsH;

  // first add in any wire electrodes
    
  Array1<Point> new_points;
  Array1<double> new_values;
  for (int h=0; h<ielecH.size(); h++) 
  {
    FieldInformation qsfi(ielecH[h]);
    if (!(qsfi.is_quadsurfmesh())) continue;
    
    VField* elecFld = ielecH[h]->vfield();    
    VMesh* elecMesh = ielecH[h]->vmesh();
    
    VMesh::Node::iterator ni, ne;
    
    elecMesh->begin(ni);
    elecMesh->end(ne);
    while(ni != ne) 
    {
      Point p;
      elecMesh->get_center(p, *ni);
      new_points.add(p);
      double val;
      elecFld->get_value(val,*ni);
      new_values.add(val);
      ++ni;
    }
  }
  
  if (new_points.size()) 
    insertNodesIntoTetMesh(dirichlet_nodes, field, new_points, new_values);
  
  // now add in the surface electrodes
  for (int h=0; h<ielecH.size(); h++) 
  {
    FieldInformation cfi(ielecH[h]);
    if (!(cfi.is_curvemesh())) continue;

    VField* elecFld = ielecH[h]->vfield();    
    VMesh* elecMesh = ielecH[h]->vmesh();

    VMesh::Node::iterator ni, ne;
    elecMesh->begin(ni);
    elecMesh->end(ne);
    
    double voltage;
    elecFld->get_value(voltage,*ni);
    std::vector<Point> outer;
    std::vector<Point> inner;
    
    for (int i=0; i<6; i++) ++ni; // the first six nodes are the handle
    while(ni != ne)
    {
      Point p;
      elecMesh->get_center(p, *ni);
      double val;
      elecFld->value(val,*ni);
      if (val == 0.0) inner.push_back(p);
      else outer.push_back(p);
      ++ni;
    }
    inner.push_back(inner[0]);
    outer.push_back(outer[0]);

    VMesh* electrodeElements = 0;
    if (h == 0) 
    {
      elecElemsH = CreateMesh(TETVOLMESH_E); 
      electrodeElements = elecElemsH->vmesh();
    }

    std::string active_side;
    elecFld->get_property("active_side", active_side);
    
    // modify mergedDirichletNodes, and mergedMesh to have new nodes
    insertContourIntoTetMesh(dirichlet_nodes, mesh,
           inner, outer, active_side,
           voltage, electrodeElements);
  }
  
  Transform t;
  t.load_identity();
  mesh->transform(t); // really just want to invalidate mesh->grid_
                      //  because nodes have been moved
  imeshH->vfield()->set_property("dirichlet", dirichlet_nodes, false);
  imeshH->vfield()->set_property("conductivity_table", conds, false);
  
  send_output_handle("TetMesh", imeshH);

  if (elecElemsH.get_rep())
  {
    FieldInformation nfi(TETVOLMESH_E,LINEARDATA_E,DOUBLE_E);
    FieldHandle elecH  = CreateField(nfi,elecElemsH);
    send_output_handle("ElectrodeElements", elecH);
  }
}


} // End namespace BioPSE
