#include <Core/Algorithms/Fields/FieldData/GetFieldData.h>
#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/FieldInformation.h>

//! Namespace used for SCIRun Algorithmic layer
namespace SCIRunAlgo {

	using namespace SCIRun;

	bool 
	IntegrateBiotSavart(AlgoBase *algo, FieldHandle& input, FieldHandle& coil, MatrixHandle& output);
	
	bool 
	GetVectorFieldDataV(AlgoBase *algo, FieldHandle& input, MatrixHandle& output);

	//! Function call to convert data from Field into Matrix data
	bool 
	BiotSavartSolverAlgo::
	run(FieldHandle& input, MatrixHandle& output)
	{
	  algo_start("BiotSavartSolver");
	  
		//! Check whether we have a field.
	  if (input.get_rep() == 0)
	  {
		error("No input source field");
		algo_end(); return (false);
	  }
	  
	  //! Construct a class with all the type information of this field
	  FieldInformation fi(input);

	  //! Check whether we have data
	  if (fi.is_nodata())
	  {
		error("Field does not contain any data");
		algo_end(); return (false);
	  }
	  
	  //! Depending on the data type select a sub algorithm
	  /*
	  if (fi.is_scalar())
		return(GetScalarFieldDataV(this,input,output));

	  else if (fi.is_vector())
		return(GetVectorFieldDataV(this,input,output));

	  else if (fi.is_tensor())
		return(GetTensorFieldDataV(this,input,output));
	*/

	  error("Unknown data type");
	  algo_end(); return (false);
	}
	
	bool 
	IntegrateBiotSavart(AlgoBase *algo, FieldHandle& input, FieldHandle& coil, MatrixHandle& output)
	{
		VMesh* vmesh = input->vmesh();
		VField* vfield = input->vfield();
		
		VMesh::size_type size = vfield->num_values();
		VMesh::size_type esize = vfield->num_evalues();
	  
		Vector val;
		int k = 0;
		for (VMesh::index_type i=0; i<size; i++)
		{
		vfield->get_value(val,i);
		//dataptr[k] = val.x();
		//dataptr[k+1] = val.y();
		//dataptr[k+2] = val.z();
		k+=3;
		}
		
	}
	
	bool 
	GetVectorFieldDataV(AlgoBase *algo, FieldHandle& input, MatrixHandle& output)
	{
	  VField* vfield = input->vfield();
	  
	  VMesh::size_type size = vfield->num_values();
	  VMesh::size_type esize = vfield->num_evalues();
	  
	  output = new DenseMatrix(size+esize,3);
	  if (output.get_rep() == 0)
	  {
		algo->error("Could not allocate output matrix");
		algo->algo_end(); return (false);
	  }
	  double* dataptr = output->get_data_pointer();

	  Vector val;
	  int k = 0;
	  for (VMesh::index_type i=0; i<size; i++)
	  {
		vfield->get_value(val,i);
		dataptr[k] = val.x();
		dataptr[k+1] = val.y();
		dataptr[k+2] = val.z();
		k+=3;
	  }
	  
	  if (vfield->basis_order() == 2)
	  {
		vfield->vmesh()->synchronize(Mesh::EDGES_E);

		for (VMesh::index_type i=0; i<esize; i++)
		{
		  vfield->get_evalue(val,i);
		  dataptr[k] = val.x();
		  dataptr[k+1] = val.y();
		  dataptr[k+2] = val.z();
		  k+=3;
		}  
	  }
	  algo->algo_end(); return (true);
	}

} // namespace SCIRunAlgo
