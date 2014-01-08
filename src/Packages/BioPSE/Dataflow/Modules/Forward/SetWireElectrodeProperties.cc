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
 *  SetWireElectrodeProperties: Insert an electrode into a finite element mesh
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
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

#include <Core/Containers/Array1.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <iostream>
#include <sstream>

namespace BioPSE {

using namespace SCIRun;

class SetWireElectrodeProperties : public Module {
  GuiDouble gui_voltage_;
  GuiDouble gui_radius_;
  GuiInt gui_nu_;
public:
  SetWireElectrodeProperties(GuiContext *context);
  virtual ~SetWireElectrodeProperties() {}
  virtual void execute();
};

DECLARE_MAKER(SetWireElectrodeProperties)

SetWireElectrodeProperties::SetWireElectrodeProperties(GuiContext *ctx)
  : Module("SetWireElectrodeProperties", ctx, Filter, "Forward", "BioPSE"),
    gui_voltage_(get_ctx()->subVar("voltage"),5), 
	gui_radius_(get_ctx()->subVar("radius"), 0.1),
    gui_nu_(get_ctx()->subVar("nu"), 5)
{
}


void
SetWireElectrodeProperties::execute()
{
  FieldHandle ielecH;
	MatrixHandle source_matrix;
  get_input_handle("Electrode", ielecH, true);
  const bool input_matrix_p =get_input_handle("Parameters",source_matrix,false);

  FieldInformation fi(ielecH);

  if (!fi.is_curvemesh()) 
  {
    error("Input electrode was not a CurveField.");
    return;
  }
	
	// get info on the parameter matrix
  size_type num_para=0;
  size_type num_col=0;
	
	if (input_matrix_p)
	{
		//cout<< "----matrix input" << endl;
		num_para=source_matrix->nrows();
		num_col=source_matrix->ncols();
	}
	
	//get parameters from parameter matrix if exists. set parameters to ui.
	if (input_matrix_p && num_para == 3 && num_col == 1) 
	{
		double* sm=source_matrix->get_data_pointer();
		
		double temp=sm[0];
		gui_voltage_.set(temp);
		temp=sm[1];
		gui_radius_.set(temp);
		temp=sm[2];
		gui_nu_.set(static_cast<int> (temp));
	}
	
	if (input_matrix_p && (num_para != 3 || num_col != 1))
	{
		error("Parameter matrix needs to be right size (1x3). This input is optional, but remember: voltage, radius, number of faces");
		return;
	}

  VMesh *mesh = ielecH->vmesh();

  double voltage = gui_voltage_.get();
  double radius = gui_radius_.get();
  int nu = gui_nu_.get();

  if (radius < 0) 
  {
    error("Radius can't be negative");
    return;
  }
  if (nu < 3) 
  {
    error("NU can't be less than 3");
    return;
  }
  double du=M_PI*2./nu;

  VMesh::Node::iterator ni, ne;
  VMesh::Node::size_type nn;
  mesh->begin(ni);  
  mesh->end(ne);
  mesh->size(nn);
  
  FieldHandle output = CreateField(QUADSURFMESH_E,LINEARDATA_E,DOUBLE_E);
  
  VMesh* quadmesh = output->vmesh();
  
  Array1<Array1<Point> > pts(nn);
  Array1<Array1<VMesh::Node::index_type> > pts_idx(nn);
  if (nn < 2) 
  {
    error("Need at least two points along Curve");
    return;
  }
  
  Array1<Point> c(nn);
  int idx=0;
  while(ni != ne) 
  {
    mesh->get_center(c[idx], *ni);
    ++ni;
    ++idx;
  }
  
  Vector up, b1, b2;
  for (Mesh::index_type i=0; i<nn; i++) 
  {
    pts[i].resize(nu);
    pts_idx[i].resize(nu);
    if (i==0) 
    {
      up=(c[i+1]-c[i]).safe_normal();
    } 
    else if (i==(nn-1)) 
    {
      up=(c[i]-c[i-1]).safe_normal();
    } 
    else 
    {
      up=(c[i+1]-c[i-1]).safe_normal();
    }
    if (i==0) 
    {
      up.find_orthogonal(b1, b2);
    } 
    else 
    {
      b2=(Cross(up, b1)).safe_normal();
      b1=(Cross(b2, up)).safe_normal();
    }
    for (int u=0; u<nu; u++) 
    {
      pts[i][u]=c[i]+(b1*cos(u*du)+b2*sin(u*du))*radius;
      pts_idx[i][u]=quadmesh->add_point(pts[i][u]);
    }
  }
  VMesh::Node::array_type nodes(4);
  
  for (unsigned int i=0; i<(unsigned int)(nn-1); i++) 
  {
    int u;
    for (u=0; u<nu-1; u++) 
    {
      nodes[0] = pts_idx[i][u];
      nodes[1] = pts_idx[i][u+1];
      nodes[2] = pts_idx[i+1][u+1];
      nodes[3] = pts_idx[i+1][u];
      quadmesh->add_elem(nodes);
    }
    
    nodes[0] = pts_idx[i][u];
    nodes[1] = pts_idx[i][0];
    nodes[2] = pts_idx[i+1][0];
    nodes[3] = pts_idx[i+1][u];
    quadmesh->add_elem(nodes);
  }

  output->vfield()->resize_values();
  output->vfield()->set_all_values(voltage); 
	
	// make parameters output matrix	
	MatrixHandle output_matrix = new DenseMatrix(3,1);
	double* P=output_matrix->get_data_pointer();
	
	double temp=gui_voltage_.get();
	P[0]=temp;
	temp=gui_radius_.get();
	P[1]=temp;
	temp=static_cast<double> (gui_nu_.get());
	P[2]=temp;
	
	send_output_handle("Parameters", output_matrix);

  send_output_handle("Electrode", output);
}

} // End namespace BioPSE
