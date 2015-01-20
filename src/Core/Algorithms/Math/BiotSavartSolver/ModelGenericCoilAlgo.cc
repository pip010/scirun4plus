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

//! Base class for algorithm
#include <Core/Algorithms/Util/AlgoBase.h>

#include <vector>
//#include <array>
#include <cassert>

#include <boost/lexical_cast.hpp>

namespace SCIRunAlgo {

using namespace SCIRun;



	//! Namespace used for concrete coil implementations
	namespace details
	{
		class BaseCoilgen
		{
			public:
				BaseCoilgen(AlgoBase* algo, ModelGenericCoilAlgo::Args args ) :
				  ref_cnt(0),
				  algo(algo),
				  innerR(args.coilRadiusInner),
				  outerR(args.coilRadiusOuter),
				  outerD(args.coilDistanceOuter),
				  current(args.wireCurrent),
				  windings(args.wireLoops),
				  coilLOD(args.coilLevelDetails)
				{

				}
				
				virtual ~BaseCoilgen()
				{
				}
							
				//! Local entry function, must be implemented by each specific kernel
				virtual void Generate(FieldHandle& meshHandle, MatrixHandle& params) const = 0;
		
				//! Global reference counting
				int ref_cnt;
			
			protected:

				//! ref to the executing algorithm context
				const AlgoBase* algo;
				
				const double innerR;
				const double outerR;
				const double current;
				const double outerD;
				const size_t windings;
				const size_t coilLOD;
				

				
				void GenPointsCircular(
					std::vector<Vector>& points,
					Vector origin, 
					double radius,
					double fromPI, 
					double toPI) const
				{
										
					assert(fromPI < toPI);
					double dPI = toPI - fromPI;
					
					//what will be the circumvence of a full 0-2*pi
					//double C = 2*dPI*r;
					//double minSegmentLenght = 0.8; 
					
					// Adaptive LOD for the number of piece-wise segments per full circle
					//int nsegments = Floor( C / minSegmentLenght); 
					
					//double anglePerSegment = (2*M_PI)/nsegments;
					
					double minPI = M_PI / (8.0d * coilLOD);
					assert(dPI > minPI);
					
					size_t nsegments = 2;
					double iPI = dPI / nsegments;
					
					while(iPI > minPI)
					{
						nsegments++;
						iPI = dPI / nsegments;
					}

					//algo->remark("LOD:  " + boost::lexical_cast<std::string>(this->coilLOD) );

					algo->remark("Segments:  " +  boost::lexical_cast<std::string>(nsegments) );
					
					GenPointsCircular(points, origin, radius, nsegments, fromPI, toPI);

				}
				
				void GenPointsCircular(
					std::vector<Vector>& points,
					Vector origin, 
					double radius,
					double nsegments,
					double fromPI, 
					double toPI) const
				{
										
					assert(fromPI < toPI);
					double dPI = toPI - fromPI;
					
					double iPI = dPI / nsegments;
					
					//std::cout << "[dPI:" << dPI << "]" << "[iPI:" << iPI << "]" << "[nsegs:" << nsegments << "]"  << std::endl << std::flush;

					for(size_t i = 0; i < nsegments; i++)
					{
						Vector p(origin.x() + radius * cos(fromPI + iPI*i), origin.y() + radius * sin(fromPI + iPI*i), origin.z());
						points.push_back(p);
					}
				}
				

				
				void BuildScirunMesh(const std::vector<Vector>& points, 
						const std::vector<size_t>& indices, 
						const std::vector<double>& values,
						FieldHandle& meshHandle) const
				{
					
					VMesh* mesh = meshHandle->vmesh();


					//! add nodes to the new mesh
					for(size_t i = 0; i < points.size(); i++)
					{
						const Point p(points[i]);
						mesh->add_point(p);
					}

					//! add edges to mesh
					VMesh::Node::array_type edge;

					for(size_t i = 0; i < indices.size(); i++)
					{
					  VMesh::Node::index_type p = indices[i];
					  edge.push_back(p);

					  if(edge.size() == 2)
					  {
						mesh->add_elem(edge);
						edge.clear();
					  }
					}

					//! add data to mesh

					VField* field = meshHandle->vfield();

					field->resize_values();
					field->set_values(values);
				}
						
				
		};


		
		//! piece-wise wire discretization
		class SingleloopCoilgen : public BaseCoilgen
		{
			public:
			
				SingleloopCoilgen( 
					AlgoBase* algo, 
					ModelGenericCoilAlgo::Args args)
					: BaseCoilgen( algo,args )
				{
					
				}
				
				~SingleloopCoilgen()
				{
				}
				
				virtual void Generate(FieldHandle& meshHandle, MatrixHandle& params) const
				{
					
					std::vector<Vector> coilPoints;
					std::vector<size_t> coilIndices;
					std::vector<double> coilValues;
	  
					//generate the two coils
					
					double radius = innerR + ((outerR - innerR) / 2.0);
					
					//LEFT
					Vector pos_L( -radius - (outerD/2), 0, 0);
					GenPointsCircular(coilPoints,pos_L,radius,0, 2*M_PI);
					GenSegmentEdges(coilPoints, coilIndices);
					GenSegmentValues(coilPoints, coilValues, current);
					
					//RIGHT
					Vector pos_R( radius + (outerD/2), 0, 0);
					GenPointsCircular(coilPoints,pos_R,radius,0, 2*M_PI);
					GenSegmentEdges(coilPoints, coilIndices);
					GenSegmentValues(coilPoints, coilValues, -current);
					
					//basic topoly assumptions needs to be correct
					assert(coilPoints.size() > 0);
					assert(coilPoints.size() == coilValues.size());
					assert(coilPoints.size()*2 == coilIndices.size());
					
					//SCIrun API creating a new mesh
					//0 data on elements; 1 data on nodes
					FieldInformation fi("CurveMesh",0,"double");
					fi.make_curvemesh();
					fi.make_constantdata();
					fi.make_scalar();

					meshHandle = CreateField(fi);

					BuildScirunMesh(coilPoints,coilIndices,coilValues,meshHandle);
				}
				
			protected:
				void GenSegmentEdges(const std::vector<Vector>& points, std::vector<size_t>& indices) const
				{
					size_t firstIDX = indices.size();
					
					size_t start = firstIDX > 0 ?  firstIDX / 2 : 0 ;
					
					for(size_t i = start; i < points.size() -1; i++)
					{
						indices.push_back(i);
						indices.push_back(i + 1);
					}
					
					size_t lastIDX = indices.size() -1;
					
					indices.push_back(indices[lastIDX]);
					indices.push_back(indices[firstIDX]);
					
					/*
					size_t firstIDX = indices.size();// > 0 ?  indices.size() -1 : 0 ; 
					
					GenEdgesLine(points,indices);
					
					size_t lastIDX = indices.size() -1;

					//indices[ 2*points.size() - 1 ] = indices[firstIDX];
					//close the loop
					indices.push_back(indices[lastIDX]);
					indices.push_back(indices[firstIDX]);
					*/
					
				}
				
				void GenSegmentValues(const std::vector<Vector>& points, std::vector<double>& values, double value) const
				{
					assert(points.size() > 0);
					
					for(size_t i = values.size(); i < points.size(); i++)
					{
						values.push_back(value);
					}
				}

		};


		class TestCoilgen : public SingleloopCoilgen
		{
			public:
				TestCoilgen( 
					AlgoBase* algo, 
					ModelGenericCoilAlgo::Args args)
					: SingleloopCoilgen( algo,args )
				{
					
				}
				
				~TestCoilgen()
				{
				}

				virtual void Generate(FieldHandle& meshHandle, MatrixHandle& params) const
				{
					
					std::vector<Vector> coilPoints;
					std::vector<size_t> coilIndices;
					std::vector<double> coilValues;
	  
					//generate a sible coil single loop
					double radius = innerR + ((outerR - innerR) / 2.0);
					Vector pos( 0, 0, 0);
					GenPointsCircular(coilPoints,pos,radius,0, 2*M_PI);
					GenSegmentEdges(coilPoints, coilIndices);
					GenSegmentValues(coilPoints, coilValues, current);
					
					//basic topoly assumptions needs to be correct
					assert(coilPoints.size() > 0);
					assert(coilPoints.size() == coilValues.size());
					assert(coilPoints.size()*2 == coilIndices.size());
					
					//SCIrun API creating a new mesh
					//0 data on elements; 1 data on nodes
					FieldInformation fi("CurveMesh",0,"double");
					fi.make_curvemesh();
					fi.make_constantdata();
					fi.make_scalar();

					meshHandle = CreateField(fi);

					BuildScirunMesh(coilPoints,coilIndices,coilValues,meshHandle);
				}
		};


		
		//! piece-wise wire discretization
		class MultiloopsCoilgen : public BaseCoilgen
		{
			public:
			
				MultiloopsCoilgen( 
					AlgoBase* algo, 
					ModelGenericCoilAlgo::Args args )
					: BaseCoilgen( algo,args )
				{
					
				}
				
				~MultiloopsCoilgen()
				{
				}
				
				virtual void Generate(FieldHandle& meshHandle, MatrixHandle& params) const
				{
					std::vector<Vector> coilPoints;
					std::vector<size_t> coilIndices;
					std::vector<double> coilValues;
					
					//LEFT
					Vector leftCenter(-outerR - (outerD/2),0,0);
					GenPointsSpiral(coilPoints, leftCenter);
					GenSegmentEdges(coilPoints, coilIndices);
					GenSegmentValues(coilPoints, coilValues, current);
					
					//RIGHT
					Vector rightCenter(outerR + (outerD/2),0,0);
					GenPointsSpiral(coilPoints, rightCenter);
					GenSegmentEdges(coilPoints, coilIndices);
					GenSegmentValues(coilPoints, coilValues, -current);

					FlipX(coilPoints,rightCenter);
					
					//basic topoly assumptions needs to be correct
					assert(coilPoints.size() > 0);
					assert(coilPoints.size() - 2 == coilValues.size());
					assert(coilPoints.size()*2 - 4 == coilIndices.size());
										
					//SCIrun API creating a new mesh
					//0 data on elements; 1 data on nodes
					FieldInformation fi("CurveMesh",0,"double");
					fi.make_curvemesh();
					fi.make_constantdata();
					fi.make_scalar();

					meshHandle = CreateField(fi);
					
					BuildScirunMesh(coilPoints,coilIndices,coilValues,meshHandle);
				}
				
			private:

				void FlipX(std::vector<Vector>& points, Vector origin) const
				{
					for(size_t i = points.size()/2; i < points.size(); i++)
					{
						Vector* v = &points[i];

						//std::cout << "flipping : " << *v ;

						v->x( origin.x() - ( v->x() - origin.x() ) );

						//std::cout << " to " << *v << std::endl;
					}
				}
				
				void GenPointsSpiral(std::vector<Vector>& points, Vector center) const

				{
					double dr = (outerR - innerR) / windings;		
					
					Vector center_offset(center.x()+ dr/2,center.y(),center.z());
					
					for (size_t i = 0; i < windings; i++)
					{
						GenPointsCircular(points, center, innerR + i*dr, 0   , M_PI);
						GenPointsCircular(points, center_offset, innerR + i*dr + dr/2, M_PI, 2*M_PI);	
					}
				
					//TODO refactor to avoid this
					Vector endp(center.x() + outerR * cos(2*M_PI), center.y() + outerR * sin(2*M_PI), center.z());
					points.push_back(endp);

				}
				
				void GenSegmentEdges(const std::vector<Vector>& points, std::vector<size_t>& indices) const
				{
					size_t start = indices.size() > 0 ?  indices.size() / 2 + 1 : 0 ;
					
					for(size_t i = start; i < points.size() -1; i++)
					{
						indices.push_back(i);
						indices.push_back(i + 1);
					}
				}
				
				void GenSegmentValues(const std::vector<Vector>& points, std::vector<double>& values, double value) const
				{
					assert(points.size() > 0);
					
					size_t start = values.size() > 0 ?  values.size() + 1 : 0 ;
					
					for(size_t i = start; i < points.size() - 1; i++)
					{
						values.push_back(value);
					}
				}
				
		};



		//! dipoles domain discretization
		// (replicating paper doi:10.1006/nimg.2002.1282)
		class DipolesCoilgen : public BaseCoilgen
		{
			public:
			
				DipolesCoilgen( 
					AlgoBase* algo, 
					ModelGenericCoilAlgo::Args args )
					: BaseCoilgen( algo,args )
				{
				}
				
				~DipolesCoilgen()
				{
				}
				
				virtual void Generate(FieldHandle& meshHandle, MatrixHandle& params) const
				{
					std::vector<Vector> dipolePoints;
					std::vector<Vector> dipoleValues;
					std::vector<double> radiiInner = preRadiiInner();
					std::vector<double> radiiOuter = preRadiiOuter();
					std::vector<double> numElements = preNumElem();
					std::vector<double> numCoupling = preNumAdjElem();
					
					
					Vector originL(- radiiOuter[15] - outerD / 2.0d, 0, 0);
					Vector originR( radiiOuter[15] + outerD / 2.0d, 0, 0 );

					
					
					
					for (size_t i = 0; i < 16; i++)
					{
						double ringRad = radiiInner[i] + (radiiOuter[i] - radiiInner[i]) / 2.0d;

						double ringArea = M_PI * ( radiiOuter[i] * radiiOuter[i] - radiiInner[i] * radiiInner[i] );

						double dipoleMoment = current * ringArea * numCoupling[i] / numElements[i];
						
						
						// LEFT COIL

						Vector dipoleNormL(0,0,1.0*dipoleMoment);

						GenPointsCircular(dipolePoints, originL, ringRad, numElements[i], 0.0d, 2*M_PI);
						
						GenSegmentValues(dipolePoints, dipoleValues, dipoleNormL );


						// RIGHT COIL
						
						Vector dipoleNormR(0,0,-1.0*dipoleMoment);

						GenPointsCircular(dipolePoints, originR, ringRad, numElements[i], 0.0d, 2*M_PI);

						GenSegmentValues(dipolePoints, dipoleValues, dipoleNormR );
					}


					
					//basic topoly assumptions needs to be correct
					assert(dipolePoints.size() > 0);
					assert(dipolePoints.size() == dipoleValues.size());

										
					//SCIrun API creating a new mesh
					//0 data on elements; 1 data on nodes
					FieldInformation fi("PointCloudMesh",1,"vector");
					fi.make_pointcloudmesh();
					fi.make_lineardata();
					fi.make_vector();

					meshHandle = CreateField(fi);
					
					BuildScirunMesh(dipolePoints,dipoleValues,meshHandle);
				}
				
		private:
				
				const std::vector<double> preRadiiInner() const
				{
					const double vals[16] = {0.00d, 0.003d, 0.007d, 0.011d, 0.015d, 0.019d, 0.023d, 0.026d, 0.028d, 0.030d, 0.032d, 0.034d, 0.036d, 0.038d, 0.040d, 0.042d};
					std::vector<double> preRadii(vals,vals+16);
					return preRadii;
				}
				
				const std::vector<double> preRadiiOuter() const
				{
					const double vals[16] = {0.003d, 0.007d, 0.011d, 0.015d, 0.019d, 0.023d, 0.026d, 0.028d, 0.030d, 0.032d, 0.034d, 0.036d, 0.038d, 0.040d, 0.042d, 0.044d};
					std::vector<double> preRadii(vals,vals+16);
					return preRadii;
				}

				const std::vector<double> preNumElem() const
				{
					const double vals[16] = {3.0d, 9.0d, 12.0d, 16.0d, 20.0d, 24.0d, 28.0d, 30.0d, 32.0d, 34.0d, 36.0d, 38.0d, 40.0d, 42.0d, 44.0d, 44.0d};
					std::vector<double> preNumElem(vals,vals+16);
					return preNumElem;

				}

				const std::vector<double> preNumAdjElem() const
				{
					const double vals[16] = {9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 8.0d, 7.0d, 6.0d, 5.0d, 4.0d, 3.0d, 2.0d, 1.0d};
					std::vector<double> preNumAdjElem(vals,vals+16);
					return preNumAdjElem;
				}
					
				void GenSegmentValues(const std::vector<Vector>& points, std::vector<Vector>& values, Vector val) const
				{
					assert(points.size() > 0);
					
					for(size_t i = values.size(); i < points.size(); i++)
					{
						values.push_back(val);
					}
				}
				
				void BuildScirunMesh(const std::vector<Vector>& points, 
						const std::vector<Vector>& values,
						FieldHandle& meshHandle) const
				{
					
					VMesh* mesh = meshHandle->vmesh();

					//! add nodes to the new mesh
					for(size_t i = 0; i < points.size(); i++)
					{
						const Point p(points[i]);
						mesh->add_point(p);
					}

					//! add data to mesh
					VField* field = meshHandle->vfield();
					field->resize_values();
					field->set_values(values);
				}
		};

	}



	bool 
	ModelGenericCoilAlgo::
	run(FieldHandle& meshFieldHandle, MatrixHandle& params, Args& args)
	{
		
		try
		{
		  std::vector<Vector> coilPoints;
		  std::vector<size_t> coilIndices;
		  
		  using namespace details;
		  
		  Handle<BaseCoilgen> helper;

		  if(args.type == 1)
		  {
			helper = new SingleloopCoilgen(this,args);
		  }
		  else if(args.type == 2)
		  {
		  	helper = new MultiloopsCoilgen(this,args);
		  }
		  else if(args.type == 3)
		  {
		  	helper = new DipolesCoilgen(this,args);
		  }
		  else if(args.type == 4)
		  {
		  	helper = new TestCoilgen(this,args);
		  }
		  else
		  {
			error("Unknown coil type!");
			algo_end(); return (false);
		  }
			
		helper->Generate(meshFieldHandle, params);
		

		}
		catch (...)
		{
		error("Error while running the algorithm ...");
		algo_end(); return (false);
		}

		return (true);
	}


} // end namespace SCIRunAlgo

