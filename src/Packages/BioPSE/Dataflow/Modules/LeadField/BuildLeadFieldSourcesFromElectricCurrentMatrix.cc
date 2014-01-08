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
 *  BuildLeadFieldSourcesFromElectricCurrentMatrix: Assign a leadfield solution vector to the field fdata
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   June 2001
 *
 */

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <iostream>
#include <stdio.h>

namespace BioPSE {

using namespace SCIRun;

class BuildLeadFieldSourcesFromElectricCurrentMatrix : public Module
{
  public:
    BuildLeadFieldSourcesFromElectricCurrentMatrix(GuiContext *context);
    virtual ~BuildLeadFieldSourcesFromElectricCurrentMatrix() {}
    virtual void execute();
};


DECLARE_MAKER(BuildLeadFieldSourcesFromElectricCurrentMatrix)


BuildLeadFieldSourcesFromElectricCurrentMatrix::BuildLeadFieldSourcesFromElectricCurrentMatrix(GuiContext *context)
  : Module("BuildLeadFieldSourcesFromElectricCurrentMatrix", context, Filter, "LeadField", "BioPSE")
{
}


void
BuildLeadFieldSourcesFromElectricCurrentMatrix::execute()
{
  // Get input matrix.
  MatrixHandle imatrix;
  get_input_handle("Data", imatrix, true);

  ColumnMatrix *cm = dynamic_cast<ColumnMatrix *>(imatrix.get_rep());
  if (!cm) 
  {
    remark("Matrix was supposed to be a ColumnMatrix.");
    return;
  }

  // Get input field.
  FieldHandle ifield;
  get_input_handle("Mesh", ifield,true);

  MeshHandle mbh = ifield->mesh();
  VMesh *tvm = ifield->vmesh();
  
  FieldInformation fi(ifield);
  if (!fi.is_tetvolmesh()) 
  {
    remark("Field was supposed to be a TetVolField.");
    return;
  }

  VMesh::Elem::size_type csize;  tvm->size(csize);

  if (cm->nrows() != csize * 3) 
  {
    remark("ColumnMatrix should be 3x as big as the number of mesh cells.");
    return;
  }

  // data looks good
  // make a new vector field and copy the matrix data into it
  FieldInformation vfi(TETVOLMESH_E,CONSTANTDATA_E,VECTOR_E);
  FieldHandle ofieldh = CreateField(vfi,ifield->mesh());
  VField* ofield = ofieldh->vfield();

  FieldInformation dfi(TETVOLMESH_E,LINEARDATA_E,DOUBLE_E);
  FieldHandle ofield2h = CreateField(dfi,ifield->mesh());
  VField *ofield2 = ofield2h->vfield();

  VMesh::Node::size_type nsize;  tvm->size(nsize);

  Array1<int> node_refs(nsize);
  Array1<double> node_sums(nsize);
  node_refs.initialize(0);
  node_sums.initialize(0);

  Array1<double> lengths(csize);
  double maxL=0;

  int i;
  for (i=0; i<cm->nrows()/3; i++) 
  {
    Vector v((*cm)[i*3], (*cm)[i*3+1], (*cm)[i*3+2]);
    ofield->set_value(v,i);
    double l=v.length();
    lengths[i]=l;
    if (l>maxL) maxL=l;
    // get all of the nodes that are attached to this cell
    // for each one, increment its count and add this value to its sum
    VMesh::Node::array_type::iterator ni;
    VMesh::Node::array_type na;
    
    tvm->get_nodes(na,VMesh::Elem::index_type(i));
    for (ni = na.begin(); ni != na.end(); ++ni) 
    {
      node_refs[*ni]++;
      node_sums[*ni]+=l;
    }
  }

  std::vector<Vector> vecs;
  std::vector<Vector> vecs2;
  
  FieldInformation pfi(POINTCLOUDMESH_E,CONSTANTDATA_E,VECTOR_E);
  
  FieldHandle pch = CreateField(pfi);  
  FieldHandle pc2h = CreateField(pfi);
  
  VMesh *pcm = pch->vmesh();
  VMesh *pcm2 = pc2h->vmesh();
  
  msg_stream_ << "\n\n\nFocusing ``spikes'':\n";
  double halfMax = maxL/2.;
  
  for (i=0; i<lengths.size(); i++) 
  {
    Point p; Vector v;
    tvm->get_center(p, VMesh::Elem::index_type(i));
    if (lengths[i]>0.000001) 
    {
      pcm->add_point(p);
      ofield->get_value(v,i);
      vecs.push_back(v);
    }
    if (lengths[i]>halfMax) 
    {
      pcm2->add_point(p);
      ofield->get_value(v,i);
      vecs2.push_back(v);
      msg_stream_ << "   Magnitude="<<lengths[i]<<" cellIdx="<<i<<" centroid="<<p<<"\n";
    }
  }
  
  msg_stream_ << "End of focusing spikes.\n";

  pch->vfield()->resize_values();
  pch->vfield()->set_values(vecs);

  pc2h->vfield()->resize_values();
  pc2h->vfield()->set_values(vecs2);
  
  for (Field::index_type ui=0; ui<nsize; ui++)
  {
    double d = node_sums[ui]/node_refs[ui];
    ofield2->set_value(d,ui);
  }
  
  send_output_handle("VectorField", ofieldh);
  send_output_handle("Field(double)", ofield2h);

  send_output_handle("PointCloud(Vector)", pch);
  send_output_handle("PrimaryPointCloud(Vector)", pc2h);
}

} // End namespace BioPSE

