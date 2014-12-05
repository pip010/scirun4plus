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

    struct Args
    {
      double wireCurrent;//Amps
      double coilRadius;//Meters
      int type;
    };

bool 
ModelGenericCoilAlgo::
run(FieldHandle& meshFieldHandle, MatrixHandle& params, Args& args)
{
	
	try
	{
	  std::vector<Vector> coilPoints;
	  std::vector<size_t> coilIndices;

	  if(args.type == 1)
	  {
		Vector center;//center at origin
		this->GenerateCircleContour(coilPoints, coilIndices, center, args.coilRadius, args.coilSegments);
	  }
	  else if(args.type == 2)
	  {
		this->GenerateFigure8ShapedCoil(coilPoints, coilIndices, args.coilRadius, args.coilDistance, args.coilSegments);
	  }
	  else if(args.type == 3)
	  {
		  this->GenerateFigure8ShapedSpiralCoil(coilPoints, coilIndices, args.coilRadius, 1);
	  }
	  else
	  {
		error("Unknown coil type!");
		algo_end(); return (false);
	  }

	  
	  FieldInformation fi("CurveMesh",0,"double");//0 data on elements; 1 data on nodes
	  fi.make_curvemesh();
	  fi.make_constantdata();
	  fi.make_scalar();

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
		
		//std::cout << " XCXCXCXCXC " << p << std::endl << std::flush;
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
	  

	  if(args.type == 1)
	  {
		std::vector<double> valz(args.coilSegments);

		for(size_t i = 0; i < args.coilSegments; i++)
		{
			valz[i] = args.wireCurrent;

		}

		meshFieldHandle->vfield()->resize_values();
		meshFieldHandle->vfield()->set_values(valz);
	  }
	  else if( args.type == 2 )
	  {
		std::vector<double> valz(2*args.coilSegments);

		for(size_t i = 0; i < 2*args.coilSegments; i++)
		{
		  if(i < args.coilSegments)
		  {
			valz[i] = args.wireCurrent;
		  }
		  else
		  {
			valz[i] = -args.wireCurrent;
		  }
		}

		meshFieldHandle->vfield()->resize_values();
		meshFieldHandle->vfield()->set_values(valz);
	  }
	  else if( args.type == 3 )
	  {
		  
		/*
		std::vector<double> valz(coilPoints.size()-1);

		for(size_t i = 0; i < coilPoints.size()-1; i++)
		{
			valz[i] = args.wireCurrent;
		}

		meshFieldHandle->vfield()->resize_values();
		meshFieldHandle->vfield()->set_values(valz);
		*/
		meshFieldHandle->vfield()->resize_values();
		meshFieldHandle->vfield()->set_all_values(args.wireCurrent);
		
	  }
	  else
	  {
		error("Unknown coil type!");
		algo_end(); return (false);
	  }


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
	double anglePerSegment = 2*M_PI/nsegments;
	points.resize(nsegments);
	indices.resize(2*nsegments);

	for(size_t i = 0; i < nsegments; i++)
	{
	points[i].Set(pos.x() + r * cos(anglePerSegment*i), pos.y() + r * sin(anglePerSegment*i), pos.z());
	}

	for(size_t i = 0, j = 0; i < nsegments; i++, j+=2)
	{
	indices[j] = i;
	indices[j+1] = i + 1;
	}

	indices[ 2*nsegments - 1 ] = indices[0];
}

void
ModelGenericCoilAlgo::
GenerateCircularContour(std::vector<Vector>& points, Vector center, double r, double fromPI, double toPI)
{
	assert(fromPI < toPI);
	double dPI = toPI - fromPI; 

	
	
	//what will be the circumvence of a full 0-2*pi
	//double C = 2*dPI*r;
	//double minSegmentLenght = 0.8; 
	
	// Adaptive LOD for the number of piece-wise segments per full circle
	//int nsegments = Floor( C / minSegmentLenght); 
	
	//double anglePerSegment = (2*M_PI)/nsegments;
	
	double minPI = M_PI / 16;
	assert(dPI > minPI);
	
	size_t nsegments = 2;
	double iPI = dPI / nsegments;
	
	while(iPI > minPI)
	{
		nsegments++;
		iPI = dPI / nsegments;
	}
	
	//points.resize(nsegments);
	
	
	//std::cout << "[dPI:" << dPI << "]" << "[iPI:" << iPI << "]" << "[nsegs:" << nsegments << "]" << std::endl <<std::flush;

	for(size_t i = 0; i < nsegments; i++)
	{
		Vector p(center.x() + r * cos(fromPI + iPI*i), center.y() + r * sin(fromPI + iPI*i), center.z());
		points.push_back(p);
		
		//std::cout << p << std::endl;
		//points[i].Set(center.x() + r * cos(fromPI + iPI*i), center.y() + r * sin(fromPI + iPI*i), center.z());
	}
	
	std::cout << " AAAAAAAA " << std::flush;
}
/*
void
ModelGenericCoilAlgo::
GenerateCircularContour(VMesh* mesh,VField* field, double r, double fromPI, double toPI)
{
	assert(fromPI < toPI);
	double dPI = toPI - fromPI; 

	//what will be the circumvence of a full 0-2*pi
	//double C = 2*dPI*r;
	//double minSegmentLenght = 0.8; 
	
	// Adaptive LOD for the number of piece-wise segments per full circle
	//int nsegments = Floor( C / minSegmentLenght); 
	
	//double anglePerSegment = (2*M_PI)/nsegments;
	
	double minPI = M_PI / 16;
	assert(dPI > minPI);
	
	size_t nsegments = 2;
	double iPI = dPI / nsegments;
	
	while(iPI > minPI)
	{
		nsegments++;
		iPI = dPI / nsegments;
	}
	
	std::cout << "[dPI:" << dPI << "]" << "[iPI:" << iPI << "]" << "[nsegs:" << nsegments << "]" << std::endl <<std::flush;

VMesh::Node::array_type edge;

size_type count = mesh->get_ni();



}
*/

std::vector<Vector> 
ModelGenericCoilAlgo::
ComposePointsForCurve(std::vector<Vector>& points1, std::vector<Vector>& points2)
{
	std::vector<Vector> result(points1);
	result.insert(result.end(), points2.begin(), points2.end());
	return result;
}

std::vector<size_t> 
ModelGenericCoilAlgo::
ComposeIndicesForCurve(std::vector<size_t>& indices1,std::vector<size_t>& indices2)
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

	points = ComposePointsForCurve(coilPoints_L, coilPoints_R);

	indices = ComposeIndicesForCurve(coilIndices_L, coilIndices_R);
}

void
ModelGenericCoilAlgo::
GenerateFigure8ShapedSpiralCoil(std::vector<Vector>& points, std::vector<size_t>& indices, double r, double loops)
{

	Vector origin(0,0,0);

	GenerateCircularContour(points, origin, r, 0, 2*M_PI);
	
	indices.resize(2*points.size());
	
	for(size_t i = 0, j = 0; i < points.size(); i++, j+=2)
	{
	indices[j] = i;
	indices[j+1] = i + 1;
	}
	
	indices[ indices.size() - 1 ] = indices[0];
	
	//std::cout << " BBBB " << points.size() << std::flush;
}

void
ModelGenericCoilAlgo::
GenerateFigure8ShapedSpiralCoil(VMesh* mesh,VField* field, Args args)
{


}

} // end namespace SCIRunAlgo

