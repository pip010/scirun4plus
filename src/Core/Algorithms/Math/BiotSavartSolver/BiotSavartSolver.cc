
#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Core/Math/MiscMath.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Core/Thread/Barrier.h>
#include <Core/Thread/Thread.h>

#include <string>
#include <cassert>

#include <boost/lexical_cast.hpp>

//! Namespace used for SCIRun Algorithmic layer
namespace SCIRunAlgo {

	using namespace SCIRun;
	
	//! Namespace used for concrete kernel implementations
	namespace details
	{
		class KernelBase
		{

		public:
			KernelBase(AlgoBase* algo, int t) :
			  ref_cnt(0),
			  algo(algo),
			  barrier("BSV KernelBase Barrier"),
			  mutex("BSV KernelBase Mutex"),
			  typeOut(t),
			  matOut(0)
			{
			}
			
			virtual ~KernelBase()
			{
			}
			
			//! Local entry function, must be implemented by each specific kernel
			virtual bool Integrate(FieldHandle& mesh, FieldHandle& coil, MatrixHandle& outdata) = 0;

	
			//! Global reference counting
			int ref_cnt;
			
		protected:

			//! ref to the executing algorithm context
			AlgoBase* algo;

			//! model miscs.
			VMesh* vmesh;
			VField* vfield;
			size_type modelSize;

			//! coil miscs.
			VMesh* vcoil;
			VField* vcoilField;
			size_type coilSize;

			//! parallel essential primitives 
			Barrier barrier;
			Mutex mutex;
			unsigned int numprocessors;
			std::vector<bool> success;
			
			//! output Field
			int typeOut;
			DenseMatrix *matOut;
			MatrixHandle matOutHandle;

			bool PreIntegration( FieldHandle& mesh, FieldHandle& coil )
			{
					this->vmesh = mesh->vmesh();
					assert(vmesh);

					this->vcoil = coil->vmesh();
					assert(vcoil);

					this->vfield = mesh->vfield();
					assert(vfield);

					this->vcoilField = coil->vfield();
					assert(vcoilField);


					this->numprocessors = Thread::numProcessors();

					int numproc = this->algo->get_int("num_processors");

					if (numproc > 0) 
					{ 
						numprocessors = numproc; 
					}
					
					#ifdef _DEBUG
					//! DEBUG when we want to test with one CPU only
					numprocessors = 1;
					#endif
					
					algo->remark("number of processors:  " + boost::lexical_cast<std::string>(this->numprocessors));
					
					success.resize(numprocessors,true);
					
					//! get number of nodes for the model
					modelSize = vmesh->num_nodes();
					assert(modelSize > 0);
					
					try
					{			
						matOut = new DenseMatrix(static_cast<int>(modelSize),3);
						matOutHandle = matOut;
					}
					catch (...)
					{
						algo->error("Error alocating output matrix");
						return (false);
					}
					
					return (true);
			}
			
			bool PostIntegration( MatrixHandle& outdata )
			{
				//! check for error
				for (size_t j=0; j<success.size(); j++)
				{
					if (success[j] == false) return (false);
				}

				outdata = matOutHandle;
				
				return (true);
			}
		};
		

		class PieceWiseKernel : public KernelBase
		{
			public:
			
				PieceWiseKernel( AlgoBase* algo, int t ) : KernelBase(algo,t)
				{
					//we keep last calculated step
					//however if segments lenght varies,
					//it makes more sense to keep a look-up table of previous steps for given lenght
					autostep = 0.1d;
					extstep = -1.0d;
					
				}
				
				~PieceWiseKernel()
				{
				}
				
				//! Complexity O(M*N) ,where M is the number of nodes of the model and N is the numbder of nodes of the coil
				virtual bool Integrate(FieldHandle& mesh, FieldHandle& coil, MatrixHandle& outdata)
				{

					if(!PreIntegration(mesh,coil))
					{
						return (false);
					}

					vmesh->synchronize(Mesh::NODES_E | Mesh::EDGES_E);
					
					VMesh::Node::array_type enodes;
					Point enode1;
					Point enode2;
					
					//! get number of nodes for the model
					//modelSize = vmesh->num_nodes();

					//! get numbder of nodes for the coil
					coilSize = vcoil->num_nodes();

					//! basic assumption
					assert(modelSize > 0 && coilSize > 1);
					
					coilNodes.clear();
					coilNodes.reserve(coilSize);
					

					for(VMesh::Edge::index_type i = 0; i < vcoil->num_edges(); i++)
					{
						vcoil->get_nodes(enodes,i);
						//std::cout << enodes[0] << "; "<< enodes[1] << "; "<< enodes[2] << "; "<< enodes[3] << "; "<< enodes[4] << "; "<< enodes[5] << "; "<< enodes[6] << "; "<< enodes[7] << "; "<< enodes[8] << "; "<< enodes[9] <<"; "<< enodes[10] << std::endl;
						vcoil->get_point(enode1,enodes[0]);
						vcoil->get_point(enode2,enodes[1]);
						coilNodes.push_back(Vector(enode1));
						coilNodes.push_back(Vector(enode2));
					}

					//! Start the multi threaded
					Thread::parallel(this, &PieceWiseKernel::ParallelKernel, numprocessors);
					
					
					return PostIntegration(outdata);
				}

				void SetIntegrationStep(double step)
				{
					assert(step >= 0.0d);
					extstep = step;
				}

				double GetIntegrationStep() const
				{
					return extstep;
				}
				
				
			private:

				//! integration step, will auto adapt
				double autostep;

				//! integration step, externally provided
				double extstep;

				//! keep nodes on the coil cached
				std::vector<Vector> coilNodes;

				//! buffer of points used for integration
				//std::vector<Vector> integrPoints;
				
				//! execute in parallel
				void ParallelKernel(int proc_num)
				{

					assert(proc_num >= 0);

					int cnt = 0;
					double current = 1.0; 
					Point modelNode;

					const index_type begins = (modelSize * proc_num) / numprocessors;
					const index_type ends  = (modelSize * (proc_num+1)) / numprocessors;

					assert( begins <= ends );

					//! buffer of points used for integration
					std::vector<Vector> integrPoints;
					integrPoints.reserve(256);

					//! keep previous step length
					//! used for optimization purpose
					double prevSegLen = 123456789.12345678;

					//! number of integration points
					int nips = 0;

					try{

						for(index_type iM = begins; 
							iM < ends; 
							iM++)
						{
							vmesh->get_node(modelNode,iM); 

							// result
							Vector F;

							for( size_t iC0 = 0, iC1 =1, iCV = 0; 
								iC0 < coilNodes.size(); 
								iC0+=2, iC1+=2, iCV++)
							{
								vcoilField->get_value(current,iCV);

								current = current == 0.0 ? 1.0 : current;

								Vector coilNodeThis;
								Vector coilNodeNext;

								if(current >= 0.0)
								{
									coilNodeThis = coilNodes[iC0];
									coilNodeNext = coilNodes[iC1];
								}
								else
								{
									coilNodeThis = coilNodes[iC1];
									coilNodeNext = coilNodes[iC0];
								}

								//std::cout << "\t coil_node: THIS[" << coilNodeThis << "] NEXT[" <<  coilNodeNext << "]   E:" << current << std::endl;//DEBUG

								//! Length of the curve element
								Vector diffNodes = coilNodeNext - coilNodeThis;
								double newSegLen = diffNodes.length();

								//first check if externally suplied integration step is available and use it
								if(extstep > 0)
								{
									nips = newSegLen / extstep;
									//std::cout << "external integration" << std::endl;
								}
								else
								{
									//! optimization
									//! only rexompute integration step only if segment length changes
									if( Abs(prevSegLen - newSegLen ) > 0.00000001 )
									{
										prevSegLen = newSegLen;

										//auto adaptive integration step calculation
										nips =  AdjustNumberOfIntegrationPoints(newSegLen);
										
										//std::cout << "\t\t integration step: " << newSegLen / static_cast<double>(nips) << std::endl;//DEBUG
										
										//algo->status( "integration step: " + boost::lexical_cast<std::string>(newSegLen / static_cast<double>(nips)) );
									}
								}


								//assert( nips > 2 );
								if( nips < 3 )
								{
									algo->warning("integration step too big");
								}


								//std::vector<Vector> integrPoints(nips);
								
								integrPoints.clear();

								
								//! curve segment discretization
								for(int iip = 0; iip < nips; iip++)
								{
									double interpolant = static_cast<double>(iip) / static_cast<double>(nips);
									Vector v = Interpolate( coilNodeThis, coilNodeNext, interpolant );
									integrPoints.push_back( v );
									//std::cout << "\t\t integration point: " << integrPoints[iip] << std::endl;//DEBUG
								}


								//! integration step over line segment				
								for(int iip = 0; iip < nips -1; iip++)								
								{
									//double M_MU = 4*M_PI*1.0e-7;
									
									//! Vector connecting the infinitesimal curve-element			
									Vector Rxyz = (integrPoints[iip] + integrPoints[iip+1] ) / 2  - Vector(modelNode);

									//! Infinitesimal curve-element components
									Vector dLxyz = integrPoints[iip+1] - integrPoints[iip];

									//double dLn = dLxyz.length();
									double Rn = Rxyz.length();
									//double Rn = Rxyz.normalize();
									
									//! check for distance between coil and model close to zero
									//! it might cause numerical stability issues with respect to the cross-product
									if(Rn < 0.00001)
									{
										algo->warning("coil<->model distance approaching zero!");
									}

									if(typeOut == 1)
									{
										//! Biot-Savart Magnetic Field
										F +=  1.0e-7 * Cross( Rxyz, dLxyz ) * ( Abs(current) / (Rn*Rn*Rn) );
										//Vector dB = Cross(Rxyz,dLxyz) * ( abs(current)/4/M_PI/Rn/Rn/Rn );	
									
									}	
								
									if(typeOut == 2)
									{
										//! Biot-Savart Magnetic Vector Potential Field
										F += 1.0e-7 * dLxyz * ( Abs(current) / (Rn) );
									}
									
								}

								//std::cout << "\t\t B: " << _results[iM] << std::endl;//DEBUG
							}

							//std::cout << "DEBUG CUR:" << current << std::endl << std::flush;

							matOut->put(iM,0, F[0]);
							matOut->put(iM,1, F[1]);
							matOut->put(iM,2, F[2]);

							//! progress reporter
							if (proc_num == 0) 
							{
								cnt++;
								if (cnt == 200) 
								{ 
									cnt = 0; 
									algo->update_progress(iM,2*(begins-ends)); 
								}
							} 
						}

						success[proc_num] = true;
					}
					catch (...)
					{
						algo->error(std::string("PieceWiseKernel crashed while integrating"));
						success[proc_num] = false;
					}
			  
					//! check point
					barrier.wait(numprocessors);

					// Bail out if one of the processes failed
					for (size_t q=0; q<numprocessors;q++) 
						if (success[q] == false) return;
						
				}
				
				//! Auto adjust accuracy of integration
				int AdjustNumberOfIntegrationPoints(double len)
				{
					//assert(step < len);
					
					int minNP = 100;//more than 1 for sure
					int maxNP = 200;//no more than 1000
					int NP = 0;
					bool over = false;
					bool under = false;

					do
					{
						NP = Ceil( len / autostep );

						under = NP < minNP ? true : false;
						over = NP > maxNP ? true : false; 

						if(under) autostep *= 0.5d;
						if(over) autostep *= 1.5d;
						
						//DEBUG
						//std::cout << "\t integration step : " << step << ",for segment lenght :" << len <<", with number of points per segment :" << NP << std::endl;

					}while( under || over );

					return NP;
				}
			
		};
		
		//! TODO
		class VolumetricKernel : public KernelBase
		{
			public:
			
				VolumetricKernel(AlgoBase* algo, int t) : KernelBase(algo,t)
				{
				}
				
				~VolumetricKernel()
				{
				}
				
				virtual bool Integrate(FieldHandle& mesh, FieldHandle& coil, MatrixHandle& outdata)
				{
					if(!PreIntegration(mesh,coil))
					{
						return (false);
					}
					

					//! get numbder of nodes for the coil
					coilSize = vcoil->num_elems();

					//! basic assumption
					assert(modelSize > 0 && coilSize > 1);
					
					
					
/*					
//RAAANNNNDDDD
					//vcoil->get_element_center(coords_type& coords);
					int basis_order = vcoilField->basis_order();
					
					VMesh::Elem::index_type ei = 1;
					
					double evol = vcoil->get_volume(ei);
					
					Vector eval;
					 //template<class T>  inline void get_value(T& val, VMesh::Elem::index_type idx) const
					 vcoilField->get_value(eval,ei);
					
					Point ecenter;
					Point ecenter2;
					vcoil->get_center(ecenter,ei);		
					vcoilField->get_center(ecenter2, ei);
					
					assert( ecenter == ecenter2 );
					
					double elmsize = vcoil->get_element_size();
					
					int nnodes = vcoil->num_nodes_per_elem();
					
					//tets
					//XXXXXXX  0  -0.00814513  0.166667  [-1.95129 -0.486505 -11.1097]  4
					
					//lattice
					//XXXXXXX  0  0.866878  1  [-1.06667 -6.4 -23.619]  8
					//XXXXXXX  0  0.866878  1  [-1.06667 -6.4 -23.619]  8 [0 0 1]



					std::cout << "XXXXXXX  " << basis_order << "  "<< evol << "  " << elmsize << "  " << ecenter << "  " << nnodes << " " << eval << std::endl;
*/				
					
					
					
					vmesh->synchronize(Mesh::NODES_E | Mesh::EDGES_E);
					
					

					//! Start the multi threaded
					Thread::parallel(this, &VolumetricKernel::ParallelKernel, numprocessors);
					
					
					return PostIntegration(outdata);
				}
				
			private:
				
				//! execute in parallel
				void ParallelKernel(int proc_num)
				{
					assert(proc_num >= 0);

					int cnt = 0;
					Point modelNode;
					Point coilCenter;
					Vector current;
					
					const VMesh::Node::index_type begins = (modelSize * proc_num) / numprocessors;
					const VMesh::Node::index_type ends  = (modelSize * (proc_num+1)) / numprocessors;

					assert( begins <= ends );

					try{

						for(VMesh::Node::index_type iM = begins; iM < ends;	iM++)
						{
							vmesh->get_node(modelNode,iM); 

							//! accumulatedresult
							Vector F;
							
							Vector R;
							
							double evol = 0.0;
							
							double Rl;

							for(VMesh::Elem::index_type  iC = 0; iC < coilSize; iC++)
							{
								vcoilField->get_value(current,iC);
								
								vcoilField->get_center(coilCenter, iC);//auto resolve based on basis_order

								evol = vcoil->get_volume(iC);
								
								R = coilCenter - modelNode;
								
								Rl = R.length();

								if(typeOut == 1)
								{
									//! Biot-Savart Magnetic Field
									//F += Cross( Rxyz, dLxyz ) * ( std::abs(current) / (4.0*M_PI*Rn*Rn*Rn) );	
									F += Cross ( current , R ) * ( evol / (4.0 * M_PI * Rl) );
								}	
							
								if(typeOut == 2)
								{
									//! Biot-Savart Magnetic Vector Potential Field
									//F += dLxyz * ( std::abs(current) / (4.0*M_PI*Rn) );
									F += current * ( evol / (4.0 * M_PI * Rl) );
								}
									
							}

							matOut->put(iM,0, F[0]);
							matOut->put(iM,1, F[1]);
							matOut->put(iM,2, F[2]);

							//! progress reporter
							if (proc_num == 0) 
							{
								cnt++;
								if (cnt == 200) 
								{ 
									cnt = 0; 
									algo->update_progress(iM,2*(begins-ends)); 
								}
							} 
						}

						success[proc_num] = true;
					}
					catch (...)
					{
						algo->error(std::string("VolumetricKernel crashed while integrating"));
						success[proc_num] = false;
					}
			  
					//! check point
					barrier.wait(numprocessors);

					// Bail out if one of the processes failed
					for (size_t q=0; q<numprocessors;q++) 
						if (success[q] == false) return;
						
				}
		};
		
		//! Magnetic Dipoles solver
		class DipolesKernel : public KernelBase
		{
			public:
			
				DipolesKernel(AlgoBase* algo, int t) : KernelBase(algo,t)
				{
				}
				
				~DipolesKernel()
				{
				}
				
				virtual bool Integrate(FieldHandle& mesh, FieldHandle& coil, MatrixHandle& outdata)
				{
					if(!PreIntegration(mesh,coil))
					{
						return (false);
					}

					//algo->remark(std::string("[ Dipole Kernel ]"));
					

					//! get numbder of nodes for the coil
					coilSize = vcoil->num_elems();

					//! basic assumption
					assert(modelSize > 0 && coilSize > 1);
						
					
					//needed?
					vmesh->synchronize(Mesh::NODES_E | Mesh::EDGES_E);
					
					

					//! Start the multi threaded
					Thread::parallel(this, &DipolesKernel::ParallelKernel, numprocessors);
					
					
					return PostIntegration(outdata);
				}
				
			private:
				
				//! execute in parallel
				void ParallelKernel(int proc_num)
				{
					assert(proc_num >= 0);

					int cnt = 0;
					Point modelNode;
					Point dipoleLocation;
					Vector dipoleMoment;
					
					const VMesh::Node::index_type begins = (modelSize * proc_num) / numprocessors;
					const VMesh::Node::index_type ends  = (modelSize * (proc_num+1)) / numprocessors;

					assert( begins <= ends );

					try{


						for(VMesh::Node::index_type iM = begins; iM < ends;	iM++)
						{
							vmesh->get_node(modelNode,iM); 

							//! accumulated result
							Vector F;
							
							Vector R;
							
							double Rl;

							for(VMesh::Elem::index_type  iC = 0; iC < coilSize; iC++)
							{
								vcoilField->get_value(dipoleMoment,iC);
								
								vcoilField->get_center(dipoleLocation, iC);//auto resolve based on basis_order

								
								R = dipoleLocation - modelNode;
								
								Rl = R.length();

								if(typeOut == 1)
								{
									//! Biot-Savart Magnetic Field
									F += 1.0e-7 * ( 3 * R * Dot ( dipoleMoment, R ) / (Rl*Rl*Rl*Rl*Rl) - dipoleMoment / (Rl*Rl*Rl) ) ; 
								}	
							
								if(typeOut == 2)
								{
									//! Biot-Savart Magnetic Vector Potential Field
									F += 1.0e-7 * Cross ( dipoleMoment , R ) / (Rl*Rl*Rl) ;
								}
									
							}

							matOut->put(iM,0, F[0]);
							matOut->put(iM,1, F[1]);
							matOut->put(iM,2, F[2]);

							//! progress reporter
							if (proc_num == 0) 
							{
								cnt++;
								if (cnt == 200) 
								{ 
									cnt = 0; 
									algo->update_progress(iM,2*(begins-ends)); 
								}
							} 
						}

						success[proc_num] = true;
					}
					catch (...)
					{
						algo->error(std::string("DipoleKernel crashed while integrating"));
						success[proc_num] = false;
					}
			  
					//! check point
					barrier.wait(numprocessors);

					// Bail out if one of the processes failed
					for (size_t q=0; q<numprocessors;q++) 
						if (success[q] == false) return;
						
				}
		};
		
	}//! end namespace detail





	//! Run the global algorithm routine
	bool 
	BiotSavartSolverAlgo::
	run(FieldHandle& mesh, FieldHandle& coil, int outtype, MatrixHandle& outdata)
	{
		algo_start("BiotSavartSolver");

		//! Check whether we have domain mesh.
		if (mesh.get_rep() == 0)
		{
			error("No input domain field");
			algo_end(); return (false);
		}
		  
	  	if (coil.get_rep() == 0)
		{
			error("No input coil source field");
			algo_end(); return (false);
		}
		
		if (coil->vfield()->basis_order()  == -1)
		{
			error("Need data on coil mesh.");
			algo_end(); return (false);
		}
	  
		using namespace details;
	  
		Handle<KernelBase> helper;
		

		if( coil->vmesh()->is_curvemesh() )
		{
			if(coil->vfield()->is_constantdata() && coil->vfield()->is_scalar())
			{
				helper = new PieceWiseKernel(this, outtype);

				PieceWiseKernel* pwk = static_cast<PieceWiseKernel*>(helper.get_rep());

				//! deligate the externally defined inteegration step
				pwk->SetIntegrationStep(this->istep);
			}
			else
			{
				error("Curve mesh expected with constant scalar data.");
				algo_end(); return (false);
			}
		}
		else if(coil->vmesh()->is_pointcloudmesh())
		{
			if((coil->vfield()->is_lineardata() || coil->vfield()->is_constantdata() ) && coil->vfield()->is_vector())
			{
				helper = new DipolesKernel(this, outtype);
			}
			else
			{
				error("pointcloud expected with linear vector data.");
				algo_end(); return (false);
			}
		}
		else if( coil->vmesh()->is_volume() )
		{
			//if( (coil->vfield()->is_lineardata() || coil->vfield()->is_constantdata() ) && coil->vfield()->is_vector() )
			if(  coil->vfield()->is_constantdata() && coil->vfield()->is_vector() )
			{
				helper = new VolumetricKernel(this, outtype);
			}
			else
			{
				error("Volumetric mesh expected with constant vector data.");
				algo_end(); return (false);
			}
		}
		else
		{
			error("Unsupported mesh type! Only curve or volumetric.");
			algo_end(); return (false);
		}
		
		//! do the solving
		if( !helper->Integrate(mesh,coil,outdata) )
		{
			error("Aborted during integration");
			algo_end(); return (false);
		}

		algo_end(); return (true);
	}
	
	

} // end namespace SCIRunAlgo
