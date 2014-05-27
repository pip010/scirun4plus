
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
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>
#include <Core/Datatypes/ColumnMatrix.h>
#include <Core/Datatypes/DenseColMajMatrix.h>
#include <Core/Datatypes/DenseMatrix.h>

#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

#include <Core/Algorithms/Math/BiotSavartSolver/ModelGenericCoilAlgo.h>

namespace SCIRun {

	using namespace SCIRun;

	class ModelGenericCoil : public Module
	{
	  public:

		ModelGenericCoil(GuiContext*);
		
		virtual ~ModelGenericCoil() {}

		virtual void execute();

	  private:
  	    
  	    GuiDouble wireCurrentTCL;
		GuiDouble coilRadiusTCL;
    	GuiString typeTCL;
		 
		SCIRunAlgo::ModelGenericCoilAlgo algo;
	};


	DECLARE_MAKER(ModelGenericCoil)

	ModelGenericCoil::ModelGenericCoil(GuiContext* ctx) :
		Module("ModelGenericCoil", ctx, Source, "Math", "SCIRun"),
		wireCurrentTCL(ctx->subVar("wireCurrentTCL")),
		coilRadiusTCL(ctx->subVar("coilRadiusTCL")),
		typeTCL(ctx->subVar("typeTCL")) 
	{
		algo.set_progress_reporter(this);
	}

	void
	ModelGenericCoil::execute()
	{
		SCIRunAlgo::ModelGenericCoilAlgo::Args algoArgs;
		MatrixHandle omatrix;
		FieldHandle ofield;

		//TCLInterface::execute(get_id() + " update_matrixdata");

		//size_type nrows = static_cast<size_type>(nrows_.get());
		//size_type ncols = static_cast<size_type>(ncols_.get());

		algoArgs.wireCurrent = static_cast<double>(wireCurrentTCL.get());
		algoArgs.coilRadius = static_cast<double>(coilRadiusTCL.get());
		std::string coilType = static_cast<std::string>(typeTCL.get());

		if (coilType == "0-shaped") algoArgs.type = 1;
		else if (coilType == "8-shaped") algoArgs.type = 2;
		else algoArgs.type = 1;

		//MatrixHandle mat = new DenseColMajMatrix(nrows,ncols);
		//omatrix = mat->dense();

		//if (!oport_cached("Mesh"))
		//	{
		    update_state(Executing);
		    //std::string datalocation = guidatalocation_.get();
		   
		    
		    if (!(algo.run(ofield,omatrix,algoArgs))) return;
		  
		    // send new output if there is any:  
		    send_output_handle("Mesh",ofield);
			
			//
			send_output_handle("Matrix", omatrix);
		//}


	}

} // End namespace SCIRun
