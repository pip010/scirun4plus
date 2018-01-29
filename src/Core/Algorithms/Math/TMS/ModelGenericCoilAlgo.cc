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

#include <Core/Algorithms/Math/TMS/ModelGenericCoilAlgo.h>
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Math/MiscMath.h>

//! Base class for algorithm
#include <Core/Algorithms/Util/AlgoBase.h>

#include <vector>
//#include <array>
#include <cassert>

using namespace std;

#include <boost/lexical_cast.hpp>

namespace SCIRunAlgo {

using namespace SCIRun;



	//! Namespace used for concrete coil implementations
	namespace details
	{
		class BaseSegments
		{
			public:
				BaseSegments(				
					std::vector<Vector>& p,
					std::vector<size_t>& i,
					std::vector<double>& v):
					points(p),
					indices(i),
					values(v),
					pc(0)
				{
				}
				virtual ~BaseSegments()
				{
					this->Terminate();
				}

				//TODO consider splitting adding of value
				void AddPoint(Vector point, double value)
				{
					points.push_back(point);

					if( pc > 0)
					{
						size_t psize = points.size();
						//std::cout << "AddPoint-> psize:" << psize << " pc:" << pc << "    " << point << std::endl; 
					
						indices.push_back(psize-2);
						indices.push_back(psize-1);
						values.push_back(value);
					}					

					++pc;
				}

				virtual void Terminate()
				{
					pc = 0;
				}

				std::vector<Vector>& points;
				std::vector<size_t>& indices;
				std::vector<double>& values;
				size_t pc;

			private:

				//! Prevent copying
    			BaseSegments & operator = (const BaseSegments & other);
    			BaseSegments(const BaseSegments & other);
		};

		class ClosedSegments : public BaseSegments
		{
			public:
				ClosedSegments(
					std::vector<Vector>& p,
					std::vector<size_t>& i,
					std::vector<double>& v)
						:BaseSegments(p,i,v)
				{

				}
				virtual ~ClosedSegments()
				{
					if(pc)
					{
						this->Terminate();
					}

					//std::cout << "~ClosedSegments() Points Indices Values: " <<  points.size() << " " << indices.size() << " " << values.size() << std::endl;
					assert(points.size() > 0);
					assert(points.size() == values.size());
					assert(points.size()*2 == indices.size());
				}
				
				void Terminate()
				{
					size_t psize = points.size();

					//close the segments and make a circle
					indices.push_back(psize-1);
					indices.push_back(psize-pc);

					values.push_back(values[values.size()-1]);

					//std::cout << "ClosedSegment-Teminate-> psize:" << psize << " pc:" << pc << "    " << std::endl; 

					pc = 0;
				}
		};

		class OpenSegments : public BaseSegments
		{
			public:
				OpenSegments(
					std::vector<Vector>& p,
					std::vector<size_t>& i,
					std::vector<double>& v)
						:BaseSegments(p,i,v)
				{}
				~OpenSegments()
				{
					assert(points.size() > 0);
					assert(indices.size() == values.size()*2);
					//assert(coilPoints.size()*2 - 2*coilLayers == coilIndices.size());
				}
		};

	



		class BaseCoilgen
		{
			public:
				BaseCoilgen(AlgoBase* algo) :
				  ref_cnt(0),
				  coilLOD(3),
				  coilType(1),
				  coilLayers(1),
				  coilLayersStep(0),
				  algo(algo)
				{

				}
				
				virtual ~BaseCoilgen()
				{
				}
							
				//! Local entry function, must be implemented by each specific kernel
				virtual void Generate(FieldHandle& meshHandle) const = 0;
		
				//! Global reference counting
				int ref_cnt;
			
			protected:

				//! ref to the executing algorithm context
				const AlgoBase* algo;
				size_t coilLOD;
				size_t coilType;
				size_t coilLayers;
				double coilLayersStep;

				void GenPointsCircular(
					BaseSegments& segments,
					Vector origin, 
					double radius,
					double value,
					double fromPI,
					double toPI, 
					double extLOD = 0.0) const
				{

					double dPI = abs(toPI - fromPI);
					
					double minPI = M_PI /  ( 8.0 * coilLOD + coilLOD * extLOD );
					
					assert(dPI > minPI);
					
					size_t nsegments = 2;
					double iPI = dPI / nsegments;
					
					while(iPI > minPI)
					{
						nsegments++;
						iPI = dPI / nsegments;
					}

					algo->remark("#Segments(LOD):  " +  boost::lexical_cast<std::string>(nsegments) );
					
					dPI = toPI - fromPI;
					
					iPI = dPI / nsegments;

					for(size_t i = 0; i < nsegments; i++)
					{
						Vector point(origin.x() + radius * cos(fromPI + iPI*i), origin.y() + radius * sin(fromPI + iPI*i), origin.z());
						segments.AddPoint(point,value);
					}
				}
				
				void GenPointsCircular2(
					BaseSegments& segments,
					Vector origin, 
					double radius,
					double value,
					double fromPI,
					double toPI) const
				{
					cout << endl;

					double dPI = abs(toPI - fromPI);
					
					
					double max_err = (1.0 / (double)coilLOD )  * 0.2d + 0.05d ;

					size_t nsegments = GetDiscretizationParams(radius, max_err);
					
					
					double iPI = dPI / nsegments;

					for(size_t i = 0; i < nsegments; ++i)
					{
						Vector point(origin.x() + radius * cos(fromPI + iPI*i), origin.y() + radius * sin(fromPI + iPI*i), origin.z());
						segments.AddPoint(point,value);
					}
				}
				
				static void CalcAnalyticalBfield(const std::vector<Vector>& points, std::vector<Vector>& data, const double R,const double current)
				{
					/// the analytical B-field: B = MU_0*I / 2 * (  R^2 / ( X^2+ R^2 ) ^ 3/2 )
					/// where R is the radius of the circular coil and X is an offset along its central/middle axis
					for(size_t i =0; i < points.size(); i++)
					{
						Vector v;
						double Z = points[i].z();
						v.z( (current * 2.0 * M_PI * 1E-7) * ( R*R / Pow( (Z*Z + R*R), 3.0 / 2.0 ) ) );
						data.push_back(v);
					}

				}
				
				static double CalcAnalyticalAtZ(const double Z, const double R,const double current = 1.0d)
				{
					/// the analytical B-field: B = MU_0*I / 2 * (  R^2 / ( X^2+ R^2 ) ^ 3/2 )
					/// where R is the radius of the circular coil and X is an offset along its central/middle axis
					return	current * 2.0 * M_PI * 1E-7 * ( R*R / Pow( (Z*Z + R*R), 3.0 / 2.0 ) );
				}
				
				static double CalcDiscreteAtZ(const double Z, const double R, const unsigned int num_segments, const double current = 1.0d)
				{
					Vector F(0.0,0.0,0.0);
					
					Vector ref_node(0,0,Z);
									
					
					double seg_pi = 2 * M_PI / num_segments;

					
					//cout << " Segment : " << seg_begin << " | "  << seg_end << endl;
				
					int nips = R > 0.0 ? 1 / R * num_segments : R * num_segments;
					
					std::vector<Vector> integrPoints;
					
					for (int s = 0; s < num_segments; ++s)
					{
					
						Vector start(cos(seg_pi*s)*R,sin(seg_pi*s)*R,0);
						Vector end(cos(seg_pi*(s+1))*R,sin(seg_pi*(s+1))*R,0);
					
					
						//! curve segment discretization
						for(int iip = 0; iip < nips; iip++)
						{
							double interpolant = static_cast<double>(iip) / static_cast<double>(nips);
							Vector v = Interpolate( start, end, interpolant );
							integrPoints.push_back( v );
							//std::cout << "\t\t integration point: " << integrPoints[iip] << std::endl;//DEBUG
						}


						//! integration step over line segment				
						for(int iip = 0; iip < nips -1; iip++)								
						{

							//! Vector connecting the infinitesimal curve-element			
							Vector Rxyz = (integrPoints[iip] + integrPoints[iip+1] ) / 2  - ref_node;

							//! Infinitesimal curve-element components
							Vector dLxyz = integrPoints[iip+1] - integrPoints[iip];
							
							double Rn = Rxyz.length();


							//! Biot-Savart Magnetic Field
							F +=  1.0e-7 * Cross( Rxyz, dLxyz ) * ( Abs(current) / (Rn*Rn*Rn) );					
						}
						
						integrPoints.clear();
					}
					
					//cout << " _____ " << F  << " _____ " << endl;
					
					return F.z();
				}
				
				// error (0-1) (that is 1-100%)
				unsigned int GetDiscretizationParams(const double radius, const double error, unsigned int num_segments = 8 ) const
				{

					if(num_segments > 1024) return 1024;
					if(num_segments < 8) return 8;
					
					// distance at 10% radius
					//double z_test_distance = 0.1d * radius; 
					
					// 2.5cm
					double z_test_distance = radius > 0.0 ? 25.0 : 0.025; 
									
					
					double B_ana = CalcAnalyticalAtZ( z_test_distance, radius );
					
					double B_dis = CalcDiscreteAtZ( z_test_distance, radius, num_segments );
					
					double B_err = abs( B_dis - B_ana ) / B_ana;
					
					//double T_err = abs( B_err - error) / error;
					
					cout << " Calibration step : " << error << " | "  << B_err << " | " << radius << endl;
					
					if( abs(error - B_err) < 0.01d )
					{
						cout << " Golden step : " << abs(error - B_err) << endl;

					}					
					else if( B_err > error )
					{
						cout << " Improve step : " << num_segments << " to "  <<  2*num_segments << endl;
						
						//! double the number of segments and try again till we reach the relative error required
						return GetDiscretizationParams(radius, error, num_segments * 2 );
					}
					else if( B_err < error )
					{
						cout << " Disimprove step : " << num_segments << " to "  <<  num_segments /  2  + num_segments / 4 << endl;
						
						return GetDiscretizationParams(radius, error,  num_segments /  2  + num_segments / 4 );
					}
					
					return num_segments;

				}
				
				
				
				
				void BuildScirunMesh(
						const std::vector<Vector>& points, 
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
		class CircularWireCoilgen : public BaseCoilgen
		{
			protected:

				const double innerR;
				const double outerR;
				double current;
				const double outerD;
				const size_t windings;
				
			public:
			
				CircularWireCoilgen( 
					AlgoBase* algo, 
					ModelTMSCoilSpiralAlgo::Args args )
					: BaseCoilgen( algo ),
					  	innerR(args.coilRadiusInner),
					  	outerR(args.coilRadiusOuter),
					  	outerD(args.coilDistanceOuter),
					  	current(args.wireCurrent),
					  	windings(args.wireLoops)
				{
					coilLOD = args.coilLevelDetails;
					coilType = args.type;
					coilLayers = args.coilLayers;

					coilLayers = coilLayers == 0 ? 1 : coilLayers;

					coilLayersStep = args.coilLayersStep;
				}
				
				~CircularWireCoilgen()
				{
				}
				
				
				virtual void Generate(FieldHandle& meshHandle) const
				{
					cout << endl << endl;
					
					std::vector<Vector> coilPoints;
					std::vector<size_t> coilIndices;
					std::vector<double> coilValues;

					Vector step(0,0,coilLayersStep);
					
					double dr = (outerR - innerR) / windings;

					if(coilType == 1)
					{
						//Vector origin(0, 0, -0.5*(1.0/coilLayers));
						Vector origin(0, 0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);

						for(size_t l = 0; l < coilLayers; l++)
						{
							ClosedSegments segments(coilPoints,coilIndices,coilValues);
							//GenPointsCircular2(segments, origin, outerR, current, 0.0, 2.0*M_PI);
							GenPointsCircular(segments, origin, outerR, current, 0.0, 2.0*M_PI);
							
							origin += step;
						}

					}
					else if(coilType == 2)
					{
						Vector originLeft ( -outerR - (outerD/2), 0.0, -coilLayersStep*(coilLayers/2)  +coilLayersStep/2);
						Vector originRight(  outerR + (outerD/2), 0.0, -coilLayersStep*(coilLayers/2)  +coilLayersStep/2);

						for(size_t l = 0; l < coilLayers; l++)
						{
							
							for (size_t i = 0; i < windings; i++)
							{
								ClosedSegments segments(coilPoints,coilIndices,coilValues);
								//GenPointsCircular2(segments, originLeft, innerR + dr + i*dr, current, 0.0 , 2*M_PI );
								GenPointsCircular(segments, originLeft, innerR + dr + i*dr, current, 0.0 , 2*M_PI );
							}
							
							originLeft += step;
						}

						for(size_t l = coilLayers; l < 2*coilLayers; l++)
						{	
							ClosedSegments segments(coilPoints,coilIndices,coilValues);

							for (size_t i = 0; i < windings; i++)
							{
								ClosedSegments segments(coilPoints,coilIndices,coilValues);
								//GenPointsCircular2(segments, originRight, innerR + dr + i*dr, -current, 0.0 , 2*M_PI);
								GenPointsCircular(segments, originRight, innerR + dr + i*dr, -current, 0.0 , 2*M_PI);
							}
							
							originRight += step;
						}

					}
					else
					{
						algo->error("coil type value expeced: 1/2 (0-shape/8-shape)");
						return;
					}
				
										
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
		class SingleloopCoilgen : public BaseCoilgen
		{
			public:
			
				SingleloopCoilgen( 
					AlgoBase* algo, 
					ModelTMSCoilSingleAlgo::Args args)
					: BaseCoilgen( algo ),
					  	radius(args.coilRadius),
					  	outerD(args.coilDistanceOuter),
					  	current(args.wireCurrent)
				{
					coilLOD = args.coilLevelDetails;
					coilType = args.type;
					coilLayers = args.coilLayers;

					coilLayers = coilLayers == 0 ? 1 : coilLayers;

					coilLayersStep = args.coilLayersStep;

					current = current / coilLayers;
				}
				
				~SingleloopCoilgen()
				{
				}
				
				virtual void Generate(FieldHandle& meshHandle) const
				{
					Vector vec;
					

					std::vector<Vector> coilPoints;
					std::vector<size_t> coilIndices;
					std::vector<double> coilValues;

					Vector step(0,0,coilLayersStep);
	  
					if(coilType == 1)
					{
						///SINGLE
						//Vector origin(0, 0, -0.5*(1.0/coilLayers));
						Vector origin(0, 0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);
						
						for(size_t l = 0; l < coilLayers; l++)
						{
							ClosedSegments segments(coilPoints,coilIndices,coilValues);

							GenPointsCircular(segments, origin, radius, current, 0.0, 2.0*M_PI);

							origin += step;
						}

						
					}
					else if(coilType == 2)
					{
						Vector originLeft( -radius - (outerD/2), 0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);
						Vector originRight( radius + (outerD/2), 0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);
						
						
						for(size_t l = 0; l < coilLayers; l++)
						{
							ClosedSegments segments(coilPoints,coilIndices,coilValues);

							///LEFT
							GenPointsCircular(segments, originLeft, radius, current, 0.0, 2.0*M_PI);
							
							originLeft += step;
						}

						for(size_t l = 0; l < coilLayers; l++)
						{
							ClosedSegments segments(coilPoints,coilIndices,coilValues);

							///RIGHT
							GenPointsCircular(segments, originRight, radius, -current, 0.0, 2.0*M_PI);
							
							originRight += step;
						}

					}
					else
					{
						algo->error("coil type value expeced: 1/2 (0-shape/8-shape)");
						return;
					}
					
					
					///SCIrun API creating a new mesh
					///0 data on elements; 1 data on nodes
					FieldInformation fi("CurveMesh",0,"double");
					fi.make_curvemesh();
					fi.make_constantdata();
					fi.make_scalar();

					meshHandle = CreateField(fi);

					BuildScirunMesh(coilPoints,coilIndices,coilValues,meshHandle);
				}
				
			protected:

				const double radius;
				double current;
				const double outerD;

		};

		
		//! piece-wise wire discretization
		class MultiloopsCoilgen : public BaseCoilgen
		{
			public:
			
				MultiloopsCoilgen( 
					AlgoBase* algo, 
					ModelTMSCoilSpiralAlgo::Args args )
					: BaseCoilgen( algo ),
					  	innerR(args.coilRadiusInner),
					  	outerR(args.coilRadiusOuter),
					  	outerD(args.coilDistanceOuter),
					  	current(args.wireCurrent),
					  	windings(args.wireLoops)
				{
					coilLOD = args.coilLevelDetails;
					coilType = args.type;
					coilLayers = args.coilLayers;

					coilLayers = coilLayers == 0 ? 1 : coilLayers;

					coilLayersStep = args.coilLayersStep;

					//no auto current adjustment for each layer
					//hidden functionality (rather be explicit up front)
					//leave it up to users judgment

					//current = current / coilLayers;
				}
				
				~MultiloopsCoilgen()
				{
				}
				
				virtual void Generate(FieldHandle& meshHandle) const
				{
					std::vector<Vector> coilPoints;
					std::vector<size_t> coilIndices;
					std::vector<double> coilValues;

					Vector step(0,0,coilLayersStep);



					if(coilType == 1)
					{
						//Vector origin(0, 0, -0.5*(1.0/coilLayers));
						Vector origin(0, 0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);

						for(size_t l = 0; l < coilLayers; l++)
						{
							OpenSegments segments( coilPoints, coilIndices, coilValues );

							///SINGLE coil
							GenPointsSpiralLeft(segments, origin);
							origin += step;
						}

					}
					else if(coilType == 2)
					{
						Vector originLeft ( -outerR - (outerD/2), 0.0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);
						Vector originRight(  outerR + (outerD/2), 0.0, -coilLayersStep*(coilLayers/2) +coilLayersStep/2);
						//Vector originRight ( -outerR - (outerD/2), 0.0, -coilLayersStep*(coilLayers/2) );
						//Vector originLeft (  outerR + (outerD/2), 0.0, -coilLayersStep*(coilLayers/2) );

						for(size_t l = 0; l < coilLayers; l++)
						{
							OpenSegments segments( coilPoints, coilIndices, coilValues );

							//LEFT coil
							GenPointsSpiralLeft(segments, originLeft);

							originLeft += step;
						}

						for(size_t l = coilLayers; l < 2*coilLayers; l++)
						{	
							OpenSegments segments( coilPoints, coilIndices, coilValues );

							//RIGHT coil
							GenPointsSpiralRight(segments, originRight);
							
							originRight += step;
						}

					}
					else
					{
						algo->error("coil type value expeced: 1/2 (0-shape/8-shape)");
						return;
					}

					
										
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

				const double innerR;
				const double outerR;
				double current;
				const double outerD;
				const size_t windings;
				
				void GenPointsSpiralLeft(OpenSegments& segments, Vector center) const

				{
					double dr = (outerR - innerR) / windings;		
					
					Vector center_offset (center.x() + dr/2, center.y(), center.z() );
					
					for (size_t i = 0; i < windings; i++)
					//for (size_t i = windings -1; i > 0; --i)
					{
						GenPointsCircular(segments, center, innerR + i*dr, current, 0   , M_PI, i );

						GenPointsCircular(segments, center_offset, innerR + i*dr + dr/2, current, M_PI, 2*M_PI, i );	
					}
				
					//TODO refactor to avoid this
					Vector endp(center.x() + outerR * cos(2*M_PI), center.y() + outerR * sin(2*M_PI), center.z());
					segments.AddPoint(endp, current);
				}

				void GenPointsSpiralRight(OpenSegments& segments, Vector center) const
				{
					double dr = (outerR - innerR) / windings;		
					
					Vector center_offset( center.x() + dr/2, center.y(), center.z() );

					for (size_t i = windings; i > 0; i--)
					//for (size_t i = windings -1; i > 0; --i)
					{
						GenPointsCircular(segments, center, innerR + i*dr, -current, M_PI, 2*M_PI, i );

						GenPointsCircular(segments, center_offset, innerR + i*dr - dr/2, -current, 0, M_PI, i );	
					}
				
					//TODO refactor to avoid this
					Vector endp(center.x() + innerR * cos(M_PI), center.y() + innerR * sin(M_PI), center.z());
					segments.AddPoint(endp, -current);
				}

				/// this is tricky but doable, the idea is to distribute the current along each coil winding 
				/// so that the top surface flux is not linear but curved (bell shaped like)
				/// this is required since in the AC profile there is inter-winding coupling increasing the resistivity of the net
				/// please see: https://en.wikipedia.org/wiki/Proximity_effect_%28electromagnetism%29
				void AdjustForProximityEffect()
				{
					
				}
				
		};



		//! dipoles domain discretization
		// (replicating paper doi:10.1006/nimg.2002.1282)
		class DipolesCoilgen : public BaseCoilgen
		{
			public:
			
				DipolesCoilgen( 
					AlgoBase* algo, 
					ModelTMSCoilDipoleAlgo::Args args )
					: BaseCoilgen( algo ),
					  	innerR(args.coilRadiusInner),
					  	outerR(args.coilRadiusOuter),
					  	outerD(args.coilDistanceOuter),
					  	current(args.totalCurrent),
					  	segments(args.numberSegments)
					  	
				{
					coilLOD = args.coilLevelDetails;
					coilType = args.type;
					coilLayers = args.coilLayers;

					coilLayers = coilLayers == 0 ? 1 : coilLayers;
					
				}
				
				~DipolesCoilgen()
				{
				}
				
				virtual void Generate(FieldHandle& meshHandle) const
				{
					std::vector<Vector> dipolePoints;
					std::vector<Vector> dipoleValues;
					std::vector<size_t> coilIndices;
					
					std::vector<double> radiiInner = preRadiiInner();
					std::vector<double> radiiOuter = preRadiiOuter();
					//std::vector<double> numElements = preNumElem(radiiInner);
					std::vector<double> numCoupling = preNumAdjElem(radiiInner);

					
					//print_vector(radiiInner);
					//print_vector(radiiOuter);
					
					//print_vector(numElements);
					//print_vector(numCoupling);
					
					assert(radiiInner.size() == radiiOuter.size());
					
					
					//algo->remark("#Rings:  " +  boost::lexical_cast<std::string>(radiiOuter.size()) + " ring-step:" + boost::lexical_cast<std::string>(lod_step_m));
					
					
					if(coilType == 1)
					{
						Vector center(0, 0, 0);

						for (size_t i = 0; i < radiiInner.size(); i++)
						{
							double ringRad = radiiInner[i] + (radiiOuter[i] - radiiInner[i]) / 2.0d;						
							
							/// SINGLE COIL								
							size_t numElements = GenPointsCircular2(dipolePoints, center, ringRad, 0.0d, 2*M_PI, segments);
							double ringArea = M_PI * ( radiiOuter[i] * radiiOuter[i] - radiiInner[i] * radiiInner[i] );							
							double dipoleMoment = (  current * ringArea * numCoupling[i] ) / numElements;
							Vector dipoleNormL(0,0,1.0*dipoleMoment);
							GenSegmentValues(dipolePoints, dipoleValues, dipoleNormL );
						}
					}
					else if(coilType == 2)
					{
						Vector originL( -radiiOuter[radiiOuter.size()-1] - outerD / 2.0d, 0, 0);
						Vector originR( radiiOuter[radiiOuter.size()-1] + outerD / 2.0d, 0, 0 );
						
						for (size_t i = 0; i < radiiInner.size(); i++)
						{
							double ringRad = radiiInner[i] + (radiiOuter[i] - radiiInner[i]) / 2.0d;
							double ringArea = M_PI * ( radiiOuter[i] * radiiOuter[i] - radiiInner[i] * radiiInner[i] );
							
							/// LEFT COIL
							size_t numElementsL = GenPointsCircular2(dipolePoints, originL, ringRad, 0.0d, 2*M_PI, segments);
							
							double dipoleMomentL = ( current * ringArea * numCoupling[i] ) / numElementsL;
							Vector dipoleNormL(0,0,1.0*dipoleMomentL);
							GenSegmentValues(dipolePoints, dipoleValues, dipoleNormL );


							/// RIGHT COIL
							size_t numElementsR = GenPointsCircular2(dipolePoints, originR, ringRad, 0.0d, 2*M_PI, segments);
							double dipoleMomentR = ( current * ringArea * numCoupling[i] ) / numElementsR;
							Vector dipoleNormR(0,0,-1.0*dipoleMomentR);
							GenSegmentValues(dipolePoints, dipoleValues, dipoleNormR );
						}

					}
					else
					{
						algo->error("coil type value expeced: 1/2 (0-shape/8-shape)");
						return;
					}
					
					
					///basic topoly assumptions needs to be correct
					assert(dipolePoints.size() > 0);
					assert(dipolePoints.size() == dipoleValues.size());

										
					///SCIrun API creating a new mesh
					///0 data on elements; 1 data on nodes
					FieldInformation fi("PointCloudMesh",1,"vector");
					fi.make_pointcloudmesh();
					fi.make_lineardata();
					fi.make_vector();

					meshHandle = CreateField(fi);
					
					BuildScirunMesh(dipolePoints,dipoleValues,meshHandle);
				}
				
		protected:

				const double innerR;
				const double outerR;
				const double current;
				const double outerD;
				const size_t segments;
				//const size_t rings;
				//const double lod_step_m;
				
				void print_vector(const std::vector<double>& v) const
				{
					std::cout << std::endl;
					for(int i=0;i<v.size();++i)
					{
						std::cout << v[i] << " "; 
					}
					std::cout << std::endl;
				}
				
				const std::vector<double> preRadiiInner() const
				{
					std::vector<double> preRadii;
					
					double step = (outerR - innerR) / segments;
					
					double d = innerR;
					
					//add first element
					//preRadii.push_back(0.00d);
					
					while( d < outerR)
					{
						preRadii.push_back(d);
						d += step;
					}
					
										
					//const double vals[16] = {0.00d, 0.003d, 0.007d, 0.011d, 0.015d, 0.019d, 0.023d, 0.026d, 0.028d, 0.030d, 0.032d, 0.034d, 0.036d, 0.038d, 0.040d, 0.042d};
					//std::vector<double> preRadii(vals,vals+16);
					return preRadii;
				}
				
				const std::vector<double> preRadiiOuter() const
				{
					std::vector<double> preRadii;
					
					double step = (outerR - innerR) / segments;
					
					double d = innerR;
					
					//add first element
					//preRadii.push_back(d);
					
					while( d < outerR)
					{
						d += step;
						preRadii.push_back(d);
					}
					
					//add last
					//preRadii.push_back(outerR);
					
					//override last to fill to exacct outer radius 
					preRadii[preRadii.size()-1u] = outerR;
					
					
					//const double vals[16] = {0.003d, 0.007d, 0.011d, 0.015d, 0.019d, 0.023d, 0.026d, 0.028d, 0.030d, 0.032d, 0.034d, 0.036d, 0.038d, 0.040d, 0.042d, 0.044d};
					//std::vector<double> preRadii(vals,vals+16);
					return preRadii;
				}

				const std::vector<double> preNumElem(std::vector<double>& radii) const
				{
					std::vector<double> preNumElem;
					
					for(size_t i = 1; i <= radii.size(); ++i)
					{
						
						size_t n = M_PI_2 / (1.0 / ( radii[i]* coilLOD ) );
						
						//size_t n = segments +  ( pow(i,1.5) / coilLOD );
						preNumElem.push_back(n);
					}
					
					//const double vals[16] = {3.0d, 9.0d, 12.0d, 16.0d, 20.0d, 24.0d, 28.0d, 30.0d, 32.0d, 34.0d, 36.0d, 38.0d, 40.0d, 42.0d, 44.0d, 44.0d};
					//std::vector<double> preNumElem(vals,vals+16);
					return preNumElem;

				}

				const std::vector<double> preNumAdjElem(std::vector<double>& radii) const
				{
					std::vector<double> preNumAdjElem;
					
					for(size_t i = 1; i <= radii.size(); ++i)
					{
						preNumAdjElem.push_back(1.0);
					}
					
					//const double vals[16] = {9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 9.0d, 8.0d, 7.0d, 6.0d, 5.0d, 4.0d, 3.0d, 2.0d, 1.0d};
					//std::vector<double> preNumAdjElem(vals,vals+16);
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

				size_t GenPointsCircular2(
					std::vector<Vector>& points,
					Vector origin, 
					double radius,
					double fromPI,
					double toPI, 
					double extLOD) const
				{

					double dPI = abs(toPI - fromPI);
					
					double minPI = M_PI /  ( 8.0 * coilLOD + coilLOD * extLOD );
					
					assert(dPI > minPI);
					
					size_t nsegments = 2;
					double iPI = dPI / nsegments;
					
					while(iPI > minPI)
					{
						nsegments++;
						iPI = dPI / nsegments;
					}

					//algo->remark("#Segments(LOD):  " +  boost::lexical_cast<std::string>(nsegments) );
					
					dPI = toPI - fromPI;
					
					iPI = dPI / nsegments;

					for(size_t i = 0; i < nsegments; i++)
					{
						Vector point(origin.x() + radius * cos(fromPI + iPI*i), origin.y() + radius * sin(fromPI + iPI*i), origin.z());
						//segments.AddPoint(point,value);
						points.push_back(point);
					}
					
					return nsegments;
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

	}// end namespace details


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool 
	ModelTMSCoilSingleAlgo::
	run(FieldHandle& meshFieldHandle, Args& args)
	{
		
		try
		{
			std::vector<Vector> coilPoints;
			std::vector<size_t> coilIndices;

			using namespace details;

			Handle<BaseCoilgen> helper;

			helper = new SingleloopCoilgen(this,args);

			helper->Generate(meshFieldHandle);

		}
		catch (...)
		{
			error("Error while running the algorithm ...");
			algo_end(); return (false);
		}

		return (true);
	}


	bool 
	ModelTMSCoilSpiralAlgo::
	run(FieldHandle& meshFieldHandle, Args& args)
	{
		
		try
		{
			std::vector<Vector> coilPoints;
			std::vector<size_t> coilIndices;

			using namespace details;

			Handle<BaseCoilgen> helper;

			//helper = new MultiloopsCoilgen(this,args);
			helper = new CircularWireCoilgen(this,args);

			helper->Generate(meshFieldHandle);

		}
		catch (...)
		{
			error("Error while running the algorithm ...");
			algo_end(); return (false);
		}

		return (true);
	}


	bool 
	ModelTMSCoilDipoleAlgo::
	run(FieldHandle& meshFieldHandle, Args& args)
	{
		
		try
		{
			std::vector<Vector> coilPoints;
			std::vector<size_t> coilIndices;

			using namespace details;

			Handle<BaseCoilgen> helper;

			helper = new DipolesCoilgen(this,args);

			helper->Generate(meshFieldHandle);

		}
		catch (...)
		{
			error("Error while running the algorithm ...");
			algo_end(); return (false);
		}

		return (true);
	}


} // end namespace SCIRunAlgo

