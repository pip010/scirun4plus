//#include <Core/Algorithms/Fields/FieldData/GetFieldData.h>
#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>

#include <Core/Datatypes/DenseMatrix.h>
//#include <Core/Datatypes/MatrixTypeConverter.h>
#include <Core/Datatypes/FieldInformation.h>
#include <string>
#include <cassert>

//! Namespace used for SCIRun Algorithmic layer
namespace SCIRunAlgo {

	using namespace SCIRun;


	
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
	  
	  //! Construct a class with all the type information of this field
	  //FieldInformation fiMesh(mesh);
	  //FieldInformation fiCoil(coil);

//NEDDED?
	  //! Check whether we have data
	  //if (fiMesh.is_nodata())
	  //{
		//error("Field does not contain any data");
		//algo_end(); return (false);
	  //}
	  
	  //! Depending on the data type select a sub algorithm
	  /*
	  if (fi.is_scalar())
		return(GetScalarFieldDataV(this,input,output));

	  else if (fi.is_vector())
		return(GetVectorFieldDataV(this,input,output));

	  else if (fi.is_tensor())
		return(GetTensorFieldDataV(this,input,output));

	*/

	  if( !IntegrateBiotSavart(mesh,coil,outmesh,outdata) )
	  {
  		error("Aborted during integration");
	  	algo_end(); return (false);
	  }

	  algo_end(); return (true);
	}
	
	bool 
	BiotSavartSolverAlgo::
	IntegrateBiotSavart(FieldHandle& mesh, FieldHandle& coil,FieldHandle& outmesh, MatrixHandle& outdata)
	{

		bool status = false;
		VMesh* vmesh = mesh->vmesh();
		assert(vmesh);
		VMesh* vcoil = coil->vmesh();
		assert(vcoil);
		VField* vfield = mesh->vfield();
		assert(vfield);
	  	FieldInformation fimesh(mesh);
	  	FieldInformation ficoil(coil);
		
		VMesh::size_type size = vfield->num_values();
		VMesh::size_type esize = vfield->num_evalues();

		//algo->remark("REMARK");
		if( !vcoil->is_curvemesh() )
		{
			error("Only curve mesh type is accepted for coil geometry.");
			return status;
		}
	  
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

	    VMesh::Node::size_type nnodes;
  		vmesh->size(nnodes);
  		
  		//just a check
  		int tmp = vmesh->num_nodes();
  		assert(tmp == nnodes);



  		vmesh->synchronize(Mesh::NODES_E);

//READ
  		Point p;
  		for(VMesh::Node::index_type i = 0; i < nnodes; i++)
  		{
  			vmesh->get_node(p,i);
  			std::cout << p << std::endl;
  		}

//WRITE
  		fimesh.make_vector();
		fimesh.make_lineardata();


	  	outmesh = CreateField(fimesh,mesh->mesh());		
		//outmesh->copy_properties(mesh.get_rep());




		if (outmesh.get_rep() == 0) 
		{
		error("Could not create output field and output interface");
		//algo_end(); 
		return status;
		}  

		VField* ofield = outmesh->vfield();

		for (VMesh::Node::index_type i=0; i< nnodes; i++)
		{
			Vector v(0.1,0.2,(double)(i*10));
			ofield->set_value(v,i);
		}

//WRITE Matrix
	  	//MatrixHandle mat = new DenseMatrixGeneric(4,4);
	  	//outdata = mat->dense();

		

		status = true;
		
		return status;
	}
	
	/*
		if (e.data.code == "sim") 
	{
		var params = e.data.params;
		var numPointsDomain = e.data.points.length;
		var numPointsCurv = e.data.curve.length;
		var L = e.data.curve;
		var results = new Array(numPointsDomain);
		
		//for each point in the input calculate induction
		for(i = 0; i < numPointsDomain; i++)
		{
			var point = e.data.points[i];
			results[i] = [0,0,0];
			
			for(c = 0; c < numPointsCurv-1; c++)
			{
				//Length of the curve element
				var diflen = numeric.sub(L[c],L[c+1]);
				var len = numeric.norm2(diflen);

				//Number of points for the curve-element discretization
				//var len_ds = numeric.div(len,params.ds);
				//var Npi = Math.ceil(len/params.ds);
				var Npi = AdjustPointDistribution(len);

				//AVOID ERROR HERE
				//if(Npi < 3){
				//	log("ERROR Integration step is too big!!");
				//}

				//Curve-element discretization
				var Lx = numeric.linspace(L[c][0], L[c+1][0], Npi);
				var Ly = numeric.linspace(L[c][1], L[c+1][1], Npi);
				var Lz = numeric.linspace(L[c][2], L[c+1][2], Npi);

				var Ldiscrete = new Array(Npi);
				for(ci = 0;ci< Npi; ci++)
				{
					Ldiscrete[ci] = [Lx[ci],Ly[ci],Lz[ci]];
				}

				//Integration
				for(s = 0;s<Npi-1;s++)
				{
					//Vector connecting the infinitesimal curve-element			
					var Rxyz = numeric.sub(Ldiscrete[s] , point);

					//Infinitesimal curve-element components
					var dLxyz = numeric.sub(Ldiscrete[s+1] , Ldiscrete[s]);

					//Modules
					var dLn = numeric.norm2(dLxyz);
					var Rn = numeric.norm2(Rxyz);

					//Biot-Savart
					var dB = numeric.mul( crossprod(Rxyz,dLxyz), 
						params.gamma/4/params.pi/Rn/Rn/Rn);

					//Debug
					//log("BS dB: "+dB);
						
					//Add increment to the main field
					results[i] = numeric.add(results[i], numeric.mul(dB, params.loops) );
				}
			}
			
			//Debug
			//log("result [i]:value ::"+i+":"+results[i]);
			
			self.postMessage({code:'progress'});
		}
		
		self.postMessage({ID: e.data.ID, code:'sim', output:results});
		self.postMessage({ID: e.data.ID, code:'finished'});
	}
	*/
	

} // namespace SCIRunAlgo
