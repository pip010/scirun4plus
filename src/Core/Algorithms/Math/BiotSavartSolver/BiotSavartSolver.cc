//#include <Core/Algorithms/Fields/FieldData/GetFieldData.h>
#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/DenseColMajMatrix.h>
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
	      results()
	    {
	    }

        //! Local entry function, Biot-Savart Contour Piece-wise integration
		bool IntegrateBiotSavart(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& outdata);

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
	    size_type coilSize;

//TBR
	    MatrixHandle rhsmatrixhandle;
	    DenseMatrix* rhsmatrix;
//

	  	// parallel essential primitives 
	    Barrier barrier;
	    Mutex mutex;
	    unsigned int numprocessors;
	    std::vector<bool> success;

	    std::vector<Vector> results;


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


		Point modelNode;
  		Point coilNode;

  		double gamma = 1.0;

		const index_type begins = (modelSize * proc_num) / numprocessors;
		const index_type ends  = (modelSize * (proc_num+1)) / numprocessors;

		try{

  		for(index_type iM = begins; 
  			iM < ends; 
  			iM++)
  		{
  			vmesh->get_node(modelNode,iM);  			
  			//std::cout << "model_node: " << modelNode << std::endl;//DEBUG

			for(index_type iC0 = 0, iC1 =1; 
				iC1 < coilSize; 
				iC0++, iC1++)
			{
				vcoil->get_node(coilNode,iC0);
				//std::cout << "\t coil_node: " << coilNode << std::endl;//DEBUG

				Vector coilNodeThis(coilNode);

				vcoil->get_node(coilNode,iC1);

				Vector coilNodeNext(coilNode);

				//Length of the curve element
				Vector diffNodes = coilNodeNext - coilNodeThis;
				double lenSegment = diffNodes.length();

				// TODO optimize via promoting to member scope in case a constant is not vaiable for varying line segments lenght
				int numIntegrPoints = AdjustNumberOfIntegrationPoints(0.5,lenSegment);

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
					double dLn = dLxyz.length();
					double Rn = Rxyz.length();

					//Biot-Savart
					Vector dB = Cross(Rxyz,dLxyz) * (gamma/4/M_PI/Rn/Rn/Rn);

					//std::cout << "\t\t dB: " << dB << std::endl;//DEBUG
						
					//Add increment to the main field
					//results[i] = numeric.add(results[i], numeric.mul(dB, params.loops) );
					results[iM] += dB;
				
				}

				//std::cout << "\t\t B: " << _results[iM] << std::endl;//DEBUG


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

			std::cout << "\t integration step : " << step << std::endl;//DEBUG TODO TBR

		}while( under || over );

		return NP;
	}

	
	//bool 
	//GetVectorFieldDataV(AlgoBase *algo, FieldHandle& input, MatrixHandle& output);

	//! Run the global algorithm routine
	bool 
	BiotSavartSolverAlgo::
	run(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& outdata)
	{
		algo_start("BiotSavartSolver");

		//! Check whether we have fields.
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

		if( !helper->IntegrateBiotSavart(mesh,coil,outmesh,outdata) )
		{
			error("Aborted during integration");
			algo_end(); return (false);
		}

		algo_end(); return (true);
	}
	
	bool 
	BSVHelper::
	IntegrateBiotSavart(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& outdata)
	{
		//Complexity O(M*N) ,where M is the number of nodes of the model and N is the numbder of nodes of the coil

		bool status = false;

		this->vmesh = mesh->vmesh();
		assert(vmesh);

		this->vcoil = coil->vmesh();
		assert(vcoil);

		this->vfield = mesh->vfield();
		assert(vfield);

	  	FieldInformation fimesh(mesh);
	  	FieldInformation ficoil(coil);
		
		//VMesh::size_type size = vfield->num_values();
		//VMesh::size_type esize = vfield->num_evalues();

		//double gamma = 1.0;


	  	this->numprocessors = Thread::numProcessors();
	  	int numproc = this->algo->get_int("num_processors");
	  	algo->remark("number of processors:  " + boost::lexical_cast<std::string>(numproc));

	  	if (numproc > 0) { numprocessors = numproc; }


	  	numprocessors = 1;
	  
	    /*
		Vector val;
		int k = 0;
		for (VMesh::index_type i=0; i<size; i++)
		{
		vfield->get_value(val,i);
		//dataptr[k] = val.x();
		//dataptr[k+1] = val.y();
		//dataptr[k+2] = val.z();
		algo->status(val.get_string());
		k+=3;
		}
		*/


		// get number of nodes for the model
  		modelSize = vmesh->num_nodes();
  		//vmesh->size( modelSize );

		// get numbder of nodes for the coil
  		coilSize = vcoil->num_nodes();
  		//vcoil->size( coilSize );
  		
  		// basic assumption
  		assert(modelSize > 0 && coilSize > 1);

  		//WRITE Matrix
  		//DenseMatrixGeneric<double> mat((int)modelSize,3);
  		//MatrixHandle matH = mat.dense();
  		MatrixHandle mat = new DenseColMajMatrix(modelSize,3);

  		vmesh->synchronize(Mesh::NODES_E);

//READ
  		results.resize(modelSize);
  		
  		

		success.resize(numprocessors,true);

		// Start the multi threaded
		Thread::parallel(this, &BSVHelper::kernel,numprocessors);
		for (size_t j=0; j<success.size(); j++)
		{
		if (success[j] == false) return (false);
		}



//WRITE
  		fimesh.make_vector();
		fimesh.make_lineardata();


	  	outmesh = CreateField(fimesh,mesh->mesh());		
		//outmesh->copy_properties(mesh.get_rep());




		if (outmesh.get_rep() == 0) 
		{
		algo->error("Could not create output field and output interface");
		//algo_end(); 
		return status;
		}  

		VField* ofield = outmesh->vfield();

		for (VMesh::Node::index_type i=0; i< modelSize; i++)
		{
			//Vector v(0.1,0.2,(double)(i*10));
			ofield->set_value(results[i],i);
		}


	  success.resize(numprocessors,true);

	  // Start the multi threaded FEMVolRHS builder.
	  //Thread::parallel( this, &BSVBuilder::parallel, numprocessors );
	  for (size_t j=0; j<success.size(); j++)
	  {
	    if (success[j] == false) return (false);
	  }

	  outdata = rhsmatrixhandle;

		

		status = true;
		
		return status;
	}



	
/*
	bool
	BSVHelper::build(FieldHandle input, 
	                          MatrixHandle vtable,
	                          MatrixHandle& output)
	{
	  // Get virtual interface to data
	  field = input->vfield();
	  mesh  = input->vmesh();

	  // Determine the number of processors to use:

	  numprocessors = Thread::numProcessors();
	  int numproc = this->algo->get_int("num_processors");
	  if (numproc > 0) { numprocessors = numproc; }
	  
	  
	  // We added a second system of adding a vector table, using a matrix
	  // Convert that matrix into the vector table
	  if (vtable.get_rep())
	  {
	    vectors_.clear();
	    DenseMatrix* mat = vtable->dense();
	    MatrixHandle temphandle = mat;
	    // Only if we can convert it into a dense matrix, otherwise skip it
	    if (mat)
	    {
	      double* data = mat->get_data_pointer();
	      size_type m = mat->nrows();
	      size_type n = mat->ncols();
	      Vector V;

	      // Case the table has isotropic values
	      if (mat->ncols() == 1)
	      {
	        for (size_type p=0; p<m;p++)
	        {
	          V[0] = data[p*n+0];
	          V[1] = data[p*n+0];
	          V[2] = data[p*n+0];
	          
	          vectors_.push_back(std::pair<std::string, Vector>("",V));
	        }
	      }
	      else if (mat->ncols() == 3)
	      {
	        for (size_type p=0; p<m;p++)
	        {
	          V[0] = data[0+p*n];
	          V[1] = data[1+p*n];
	          V[2] = data[2+p*n];

	          vectors_.push_back(std::pair<std::string, Vector>("",V));
	        }
	      }
	      
	    }
	  }
	  
	  success.resize(numprocessors,true);

	  // Start the multi threaded FEMVolRHS builder.
	  //Thread::parallel( this, &BSVBuilder::parallel, numprocessors );
	  for (size_t j=0; j<success.size(); j++)
	  {
	    if (success[j] == false) return (false);
	  }

	  output = rhsmatrixhandle;
	  
	  return (true);
	}
	*/
	

} // namespace SCIRunAlgo
