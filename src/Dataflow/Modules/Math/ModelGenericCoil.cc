
/*
 *  ModelGenericCoil.cc:  TODO DEscription
 *
 *  Written by:
 *   Petar Petrov
 *   Utrecht Medical Center
 *   University of Utrecht
 *   January 2014
 *
 */
#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseColMajMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {

	using namespace SCIRun;

	class ModelGenericCoil : public Module
	{
	  public:
		ModelGenericCoil(GuiContext*);
		virtual ~ModelGenericCoil() {}

		virtual void execute();

	  private:
		GuiInt    nrows_;
		GuiInt    ncols_;
		GuiString data_;
		GuiString clabel_;
		GuiString rlabel_;
	};


	DECLARE_MAKER(ModelGenericCoil)

	ModelGenericCoil::ModelGenericCoil(GuiContext* ctx) :
	  Module("ModelGenericCoil", ctx, Source, "Math", "SCIRun"),
	  nrows_(get_ctx()->subVar("rows"), 1),
	  ncols_(get_ctx()->subVar("cols"), 1),
	  data_(get_ctx()->subVar("data"), "{0.0}"),
	  clabel_(get_ctx()->subVar("clabel"), "{0}"),
	  rlabel_(get_ctx()->subVar("rlabel"), "{0}")
	{
	}

	void
	ModelGenericCoil::execute()
	{
	  MatrixHandle handle;
	  TCLInterface::execute(get_id() + " update_matrixdata");
	  
	  size_type nrows = static_cast<size_type>(nrows_.get());
	  size_type ncols = static_cast<size_type>(ncols_.get());
	  std::string data = data_.get();
	  
	  MatrixHandle mat = new DenseColMajMatrix(nrows,ncols);
	  
	  if (mat.get_rep() == 0)
	  {
		error("Could allocate output matrix");
		return;
	  }
	  
	  for (size_t p=0;p<data.size();p++)
	  { 
		if ((data[p] == '}')||(data[p] == '{')) data[p] = ' ';
	  }

	  std::vector<double> nums;
	  multiple_from_string(data,nums);
	  
	  double *ptr = mat->get_data_pointer();
	 
	  //TODO: quick bug fix, this module is due for rewrite.
	  for (index_type p = 0; p < nums.size(); p++)
	  { 
		ptr[p] = nums[p];
	  }

	  handle = mat->dense();
	  send_output_handle("matrix", handle);
	}

} // End namespace SCIRun
