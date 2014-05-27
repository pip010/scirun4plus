//#include <Core/Algorithms/Fields/FieldData/GetFieldData.h>
#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>

#include <Core/Datatypes/DenseMatrix.h>
//#include <Core/Datatypes/DenseColMajMatrix.h>
#include <Core/Math/MiscMath.h>
//#include <Core/Datatypes/MatrixTypeConverter.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Core/Thread/Barrier.h>
#include <Core/Thread/Thread.h>

#include <string>
#include <cassert>

#include <boost/lexical_cast.hpp>

//! Namespace used for SCIRun Algorithmic layer
namespace SCIRunAlgo {

	using namespace SCIRun;


	// Helper class
	class BSVHelper
	{
	  public:

	    // Constructor needed as Barrier needs to have name
	    BSVHelper(AlgoBase* algo) :
	      ref_cnt(0),
	      algo(algo),
	      barrier("BSVHelper Barrier"),
	      mutex("BSVHelper Mutex"),
	      matBOut(0),
	      matAOut(0)
	    {
	    }

        //! Local entry function, Biot-Savart Contour Piece-wise integration
		bool IntegrateBiotSavart(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& dataB, MatrixHandle& dataA);

	    int ref_cnt;

	  private:

		// ref to the executing algorithm context
	    AlgoBase* algo;

	    // model miscs.
	    VMesh* vmesh;
	    VField* vfield;
	    size_type modelSize;

		// coil miscs.
	    VMesh* vcoil;
	    VField* vcoilField;
	    size_type coilSize;

	  	// parallel essential primitives 
	    Barrier barrier;
	    Mutex mutex;
	    unsigned int numprocessors;
	    std::vector<bool> success;

	    // keep nodes on the coil cached
	    std::vector<Vector> coilNodes;
		
		//output B-Field
		DenseMatrix *matBOut;
		MatrixHandle matBOutHandle;

		//output A-Field
		DenseMatrix *matAOut;
		MatrixHandle matAOutHandle;


		//! TODO
	    int AdjustNumberOfIntegrationPoints(double step, double len);

	    //! Entry point for the parallel version
		void kernel(int proc_num);
	  
	};

	void
	BSVHelper::
	kernel( int proc_num )
	{
		assert(proc_num >= 0);

		int cnt = 0;
		double current = 1.0; 
		Point modelNode;

		const index_type begins = (modelSize * proc_num) / numprocessors;
		const index_type ends  = (modelSize * (proc_num+1)) / numprocessors;

		assert( begins <= ends );

		try{

	  		for(index_type iM = begins; 
	  			iM < ends; 
	  			iM++)
	  		{
	  			vmesh->get_node(modelNode,iM);  			
	  			//std::cout << "model_node: " << modelNode << std::endl;//DEBUG

				Vector A;
				Vector B;

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

					//Length of the curve element
					Vector diffNodes = coilNodeNext - coilNodeThis;
					double lenSegment = diffNodes.length();

					// TODO optimize via promoting to member scope in case a constant is not vaiable for varying line segments lenght
					int numIntegrPoints = AdjustNumberOfIntegrationPoints(0.15,lenSegment);

					std::vector<Vector> integrPoints(numIntegrPoints);
					
					// curve discretization
					for(int iip = 0; iip < numIntegrPoints; iip++)
					{
						integrPoints[iip] = Interpolate( coilNodeThis, coilNodeNext, static_cast<double>(iip) / static_cast<double>(numIntegrPoints) );

						//std::cout << "\t\t integration point: " << integrPoints[iip] << std::endl;//DEBUG
					}



					// integration step over line segment				
					for(int iip = 0; iip < numIntegrPoints -1; iip++)
					{
						//Vector connecting the infinitesimal curve-element			
						Vector Rxyz = integrPoints[iip] - Vector(modelNode);

						//Infinitesimal curve-element components
						Vector dLxyz = integrPoints[iip+1] - integrPoints[iip];

						//Modules
						//double dLn = dLxyz.length();
						double Rn = Rxyz.length();


						//Biot-Savart
						Vector dB = Cross( Rxyz, dLxyz ) * ( std::abs(current) / (4.0*M_PI*Rn*Rn*Rn) );//Vector dB = Cross(Rxyz,dLxyz) * ( abs(current)/4/M_PI/Rn/Rn/Rn );	
						Vector dA = dLxyz * ( std::abs(current) / (4.0*M_PI*Rn) );

						//std::cout << "\t\t dB: " << dB << std::endl;//DEBUG
							
						//Add increment to the B-Field
						B += dB;

						//Add increment to the A-Field
						A += dA;

					
					}

					//std::cout << "\t\t B: " << _results[iM] << std::endl;//DEBUG


				}

				//std::cout << "DEBUG CUR:" << current << std::endl << std::flush;

				matBOut->put(iM,0, B[0]);
				matBOut->put(iM,1, B[1]);
				matBOut->put(iM,2, B[2]);

				matAOut->put(iM,0, A[0]);
				matAOut->put(iM,1, A[1]);
				matAOut->put(iM,2, A[2]);

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
			algo->error(std::string("BuildFEVolRHS crashed while mapping"));
		  	success[proc_num] = false;
		}
  
		//! check point
		barrier.wait(numprocessors);

		// Bail out if one of the processes failed
		for (size_t q=0; q<numprocessors;q++) 
			if (success[q] == false) return;

  	}

	int
	BSVHelper::
	AdjustNumberOfIntegrationPoints(double step, double len)
	{
		assert(step < len);
		int minNP = 3;//more than 1 for sure
		int maxNP = 100;//no more than 1000 i guess
		int NP = 0;
		bool over = false;
		bool under = false;

		do
		{
			under = NP < minNP ? true : false;
			over = NP > maxNP ? true : false; 

			if(under) step*=0.5;
			if(over) step*=1.5;

			NP=Ceil(len/step);

			//std::cout << "\t integration step : " << step << std::endl;//DEBUG TODO TBR

		}while( under || over );

		return NP;
	}


	//! Run the global algorithm routine
	bool 
	BiotSavartSolverAlgo::
	run(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& dataB, MatrixHandle& dataA)
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

		if( !coil->vmesh()->is_curvemesh() )
		{
			error("Only curve mesh type is accepted for coil geometry.");
			return (false);
		}


		Handle<BSVHelper> helper = new BSVHelper(this);

		if( !helper->IntegrateBiotSavart(mesh,coil,outmesh,dataB,dataA) )
		{
			error("Aborted during integration");
			algo_end(); return (false);
		}

		algo_end(); return (true);
	}
	
	bool 
	BSVHelper::
	IntegrateBiotSavart(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& dataB, MatrixHandle& dataA)
	{
		//Complexity O(M*N) ,where M is the number of nodes of the model and N is the numbder of nodes of the coil

		bool status = false;

		this->vmesh = mesh->vmesh();
		assert(vmesh);

		this->vcoil = coil->vmesh();
		assert(vcoil);

		this->vfield = mesh->vfield();
		assert(vfield);

		this->vcoilField = coil->vfield();
		assert(vcoilField);

	  	FieldInformation fimesh(mesh);
	  	FieldInformation ficoil(coil);


	  	this->numprocessors = Thread::numProcessors();

	  	int numproc = this->algo->get_int("num_processors");

	  	if (numproc > 0) 
  		{ 
  			numprocessors = numproc; 
  		}

	  	algo->remark("number of processors:  " + boost::lexical_cast<std::string>(this->numprocessors));

	  	//numprocessors = 1;// DEBUG when want to test with one CPU only


		// get number of nodes for the model
  		modelSize = vmesh->num_nodes();

		// get numbder of nodes for the coil
  		coilSize = vcoil->num_nodes();
  		
  		// basic assumption
  		assert(modelSize > 0 && coilSize > 1);
  
		try
		{
			matBOut = new DenseMatrix((int)modelSize,3);
			matBOutHandle = matBOut;
		}
		catch (...)
		{
			algo->error("Error alocating output matrix");
			return (false);
		}
  
		try
		{
			matAOut = new DenseMatrix((int)modelSize,3);
			matAOutHandle = matAOut;
		}
		catch (...)
		{
			algo->error("Error alocating output matrix");
			return (false);
		}

		vmesh->synchronize(Mesh::NODES_E | Mesh::EDGES_E);

		VMesh::Node::array_type enodes;
		Point enode1;
		Point enode2;

		coilNodes.clear();
		coilNodes.reserve(coilSize);//prenuffer capacity

		for(VMesh::Edge::index_type i = 0; i < vcoil->num_edges(); i++)
		{
			vcoil->get_nodes(enodes,i);
			//std::cout << enodes[0] << "; "<< enodes[1] << "; "<< enodes[2] << "; "<< enodes[3] << "; "<< enodes[4] << "; "<< enodes[5] << "; "<< enodes[6] << "; "<< enodes[7] << "; "<< enodes[8] << "; "<< enodes[9] <<"; "<< enodes[10] << std::endl;
			vcoil->get_point(enode1,enodes[0]);
			vcoil->get_point(enode2,enodes[1]);
			coilNodes.push_back(Vector(enode1));
			coilNodes.push_back(Vector(enode2));
		}

  		//results.resize(modelSize);
		success.resize(numprocessors,true);

		// Start the multi threaded
		Thread::parallel(this, &BSVHelper::kernel,numprocessors);
		for (size_t j=0; j<success.size(); j++)
		{
		if (success[j] == false) return (false);
		}



//WRITE field
  		fimesh.make_vector();
		fimesh.make_lineardata();


	  	outmesh = CreateField(fimesh,mesh->mesh());		


		if (outmesh.get_rep() == 0) 
		{
		algo->error("Could not create output field and output interface");
		return status;
		}  

		VField* ofield = outmesh->vfield();

		for (VMesh::Node::index_type i=0; i< modelSize; i++)
		{
			Vector v( matBOut->get(i,0), matBOut->get(i,1),matBOut->get(i,2) );
			//ofield->set_value(results[i],i);
			//matResults->put(i,0,results[i].x());
			//matResults->put(i,1,results[i].y());
			//matResults->put(i,2,results[i].z());
			ofield->set_value(v,i);
		}

	    //TODO in case we keep out-mesh

    	// copy property manager
		//output->copy_properties(input.get_rep());
//WRITE end

	  success.resize(numprocessors,true);

	  for (size_t j=0; j<success.size(); j++)
	  {
	    if (success[j] == false) return (false);
	  }

	  dataB = matBOutHandle;
	  dataA = matAOutHandle;
    


		

		status = true;
		
		return status;
	}



	
	

} // namespace SCIRunAlgo
