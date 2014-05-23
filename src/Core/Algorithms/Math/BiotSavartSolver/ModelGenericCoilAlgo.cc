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
  // if (params.get_rep() == 0)
  // {
  //   error("ModelGenericCoilAlgo: Could not convert matrix into dense matrix");
  //   return (false);    
  // } 

  double R = 15.0;
  uint segments = 33;

  try
  {
      //Vector center(R+5.0, 0.0, 0.0);
      std::vector<Vector> coilPoints;
      std::vector<size_t> coilIndices;

      //this->GenerateCircleContour(coilPoints, coilIndices, center, R, segments);
      this->GenerateFigure8ShapedCoil(coilPoints, coilIndices, R, 2.5, segments);

      //size_type m = 5;//params.ncols();
      //size_type n = 4;//params.nrows();
      //double* dataptr = params.get_data_pointer();
      
      FieldInformation fi("CurveMesh",0,"double");//0 data on elements; 1 data on nodes
      fi.make_curvemesh();

    // ALT ****************
      //MeshHandle meshHandle = CreateMesh(fi,m+1,n+1,Point(0.0,0.0,0.0),Point(static_cast<double>(m+1),static_cast<double>(n+1),0.0));
      //meshFieldHandle = CreateField(fi,meshHandle);

      meshFieldHandle = CreateField(fi);    
      //output->vfield()->set_values(dataptr,m*n);

      VMesh* mesh = meshFieldHandle->vmesh();

      // add nodes
      for(size_t i = 0; i < coilPoints.size(); i++)
      {
        const Point p(coilPoints[i]);
                
        mesh->add_point(p);
      }


      VMesh::Node::array_type edge;

      for(size_t i = 0; i < coilIndices.size(); i++)
      {
          VMesh::Node::index_type p = coilIndices[i];
          edge.push_back(p);

          if(edge.size() == 2)
          {
            mesh->add_elem(edge);
            edge.clear();
          }
      }
      
      std::vector<double> valz(2*segments);

      for(size_t i = 0; i < 2*segments; i++)
      {
        if(i < segments)
        {
          valz[i] = 1;
        }
        else
        {
          valz[i] = -1;
        }
      }

      meshFieldHandle->vfield()->resize_values();
      meshFieldHandle->vfield()->set_values(valz);
  }
  catch (...)
  {
    error("Error alocating output matrix");
    algo_end(); return (false);
  }
  
  return (true);
}

void
ModelGenericCoilAlgo::
GenerateCircleContour(std::vector<Vector> &points, std::vector<size_t> &indices, Vector pos,double r,size_t nsegments)
{
  //(x + r*cos(alpha), y + r*sin(alpha)
  //double pi = atan(1)*4;

  double anglePerSegment = 2*M_PI/nsegments;
  points.resize(nsegments);
  indices.resize(2*nsegments);

  for(size_t i = 0; i < nsegments; i++)
  {
    points[i].Set(pos.x() + r * cos(anglePerSegment*i), pos.y() + r * sin(anglePerSegment*i), pos.z());
  }

  //points[nsegments] = points[0];

  for(size_t i = 0, j = 0; i < nsegments; i++, j+=2)
  {

    indices[j] = i;
    indices[j+1] = i + 1;
  }

  indices[ 2*nsegments - 1 ] = indices[0];
}

std::vector<Vector> 
ModelGenericCoilAlgo::
ConcatPointsForCurve(std::vector<Vector>& points1, std::vector<Vector>& points2)
{
   std::vector<Vector> result(points1);
   result.insert(result.end(), points2.begin(), points2.end());
   return result;
}

std::vector<size_t> 
ModelGenericCoilAlgo::
ConcatIndicesForCurve(std::vector<size_t>& indices1,std::vector<size_t>& indices2)
{
  std::vector<size_t> result(indices1);
  result.insert(result.end(), indices2.begin(), indices2.end());

  size_t offset = indices1.size() / 2;

  for(size_t i = indices1.size(); i < result.size(); i++)
  {
    result[i] += offset;
  }

  return result;
}

void
ModelGenericCoilAlgo::
GenerateFigure8ShapedCoil(std::vector<Vector>& points, std::vector<size_t>& indices, double r, double d, size_t nsegments)
{
  Vector pos_L( -r-(d/2), 0, 0);
  std::vector<Vector> coilPoints_L;
  std::vector<size_t> coilIndices_L;

  GenerateCircleContour(coilPoints_L, coilIndices_L, pos_L, r, nsegments);

  Vector pos_R( r+(d/2), 0, 0);
  std::vector<Vector> coilPoints_R;
  std::vector<size_t> coilIndices_R;

  GenerateCircleContour(coilPoints_R, coilIndices_R, pos_R, r, nsegments);
  
  points = ConcatPointsForCurve(coilPoints_L, coilPoints_R);

  indices = ConcatIndicesForCurve(coilIndices_L, coilIndices_R);
}



} // end namespace SCIRunAlgo

