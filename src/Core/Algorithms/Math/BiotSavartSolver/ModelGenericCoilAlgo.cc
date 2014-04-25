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

#include <Core/Algorithms/Math/BiotSavartSolver/ModelGenericCoilAlgo.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Math/MiscMath.h>

#include <vector>

namespace SCIRunAlgo {

using namespace SCIRun;

bool 
ModelGenericCoilAlgo::
run(FieldHandle& meshFieldHandle, MatrixHandle& params)
{
  //MatrixHandle mat = params.dense();
  if (params.get_rep() == 0)
  {
    error("ModelGenericCoilAlgo: Could not convert matrix into dense matrix");
    return (false);    
  } 



  size_type m = 5;//params.ncols();
  size_type n = 4;//params.nrows();
  //double* dataptr = params.get_data_pointer();
  
  FieldInformation fi("CurveMesh",0,"double");
  MeshHandle mesh = CreateMesh(fi,m+1,n+1,Point(0.0,0.0,0.0),Point(static_cast<double>(m+1),static_cast<double>(n+1),0.0));

      //  std::string type = info.get_mesh_type_id();
      //  return (CreateMesh(type,x,y,min,max));

  meshFieldHandle = CreateField(fi,mesh);    
  //output->vfield()->set_values(dataptr,m*n);


  /*
  if (datalocation == "Node")
  {
    FieldInformation fi("ImageMesh",1,"double");
    MeshHandle mesh = CreateMesh(fi,m,n,Point(0.0,0.0,0.0),Point(static_cast<double>(m),static_cast<double>(n),0.0));
    output = CreateField(fi,mesh);    
    output->vfield()->set_values(dataptr,m*n);
  }
  else if (datalocation == "Element")
  {
    FieldInformation fi("ImageMesh",0,"double");
    MeshHandle mesh = CreateMesh(fi,m+1,n+1,Point(0.0,0.0,0.0),Point(static_cast<double>(m+1),static_cast<double>(n+1),0.0));
    output = CreateField(fi,mesh);    
    output->vfield()->set_values(dataptr,m*n);
  }
  else
  {
    pr->error("MatrixToField: Data location information is not recognized");
    return (false);      
  }
  */





  /*
  void
ConvertMatricesToMesh::execute()
{
  MatrixHandle positionshandle;
  MatrixHandle normalshandle;
  
  get_input_handle("Mesh Positions", positionshandle,true);
  if (!get_input_handle("Mesh Normals", normalshandle, false))
  {
    remark("No input normals connected, not used.");
  }

  if (inputs_changed_ || gui_fieldbasetype_.changed() ||
      gui_datatype_.changed() || !oport_cached("Output Field"))
  {
    update_state(Executing);
    
    if (positionshandle->ncols() < 3)
    {
      error("Mesh Positions must contain at least 3 columns for position data.");
      return;
    }
    if (positionshandle->ncols() > 3)
    {
      remark("Mesh Positions contains unused columns, only first three are used.");
    }

    std::string basename = gui_fieldbasetype_.get();
    std::string datatype = gui_datatype_.get();

    FieldInformation fi("CurveMesh",1,datatype);
    if (basename == "Curve") fi.make_curvemesh();
    else if (basename == "HexVol") fi.make_hexvolmesh();
    else if (basename == "PointCloud") fi.make_pointcloudmesh();
    else if (basename == "PrismVol") fi.make_prismvolmesh();
    else if (basename == "QuadSurf") fi.make_quadsurfmesh();
    else if (basename == "TetVol") fi.make_tetvolmesh();
    else if (basename == "TriSurf") fi.make_trisurfmesh();
    
    FieldHandle result_field = CreateField(fi);
    VMesh* mesh = result_field->vmesh();

    index_type i, j;
    const size_type pnrows = positionshandle->nrows();
    for (i = 0; i < pnrows; i++)
    {
      const Point p(positionshandle->get(i, 0),
        positionshandle->get(i, 1),
        positionshandle->get(i, 2));
      mesh->add_point(p);
    }
    
    process_elements(mesh, pnrows, basename != "PointCloud");
    
    result_field->vfield()->resize_values();
    send_output_handle("Output Field", result_field);
  }
}

void 
ConvertMatricesToMesh::process_elements(VMesh* mesh, size_type positionRows, bool required)
{
  MatrixHandle elementshandle;
  if (get_input_handle("Mesh Elements", elementshandle, required))
  {
    index_type ecount = 0;
    const size_type enrows = elementshandle->nrows();
    const size_type encols = elementshandle->ncols();
    VMesh::Node::array_type nodes;

    for (index_type i = 0; i < enrows; i++)
    {
      nodes.clear();
      for (index_type j = 0; j < encols; j++)
      {
        VMesh::Node::index_type index = static_cast<index_type>(elementshandle->get(i, j));
        if (index < 0 || index >= positionRows)
        {
          if (ecount < 10)
          {
            error("Bad index found at " + to_string(i) + ", "+ to_string(j));
          }
          index = 0;
          ecount++;
        }
        nodes.push_back(index);
      }
      mesh->add_elem(nodes);
    }
    if (ecount >= 10)
    {
      error("..." + to_string(ecount-9) + " additional bad indices found.");
    }
  }
}
*/
  
  return (true);
}

std::vector<Vector>
ModelGenericCoilAlgo::
GenerateCircleContour(Vector pos,double r,uint nsegments)
{
  //(x + r*cos(alpha), y + r*sin(alpha)
  double pi = atan(1)*4;
  double anglePerSegment = 2*M_PI/nsegments;
  std::vector<Vector> circlexy(nsegments+1);

  for(int i = 0; i < nsegments; i++)
  {
    //circlexy[i] = new THREE.Vector3(x + r * Math.cos(anglePerSegment*i), y + r * Math.sin(anglePerSegment*i), z);
    circlexy[i].Set(pos.x() + r * cos(anglePerSegment*i), pos.y() + r * sin(anglePerSegment*i), pos.z());
  }

  circlexy[nsegments] = circlexy[0];

  return circlexy;
}

} // end namespace SCIRunAlgo

