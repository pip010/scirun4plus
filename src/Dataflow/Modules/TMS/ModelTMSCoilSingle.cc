
/*
 *  ModelTMSCoilSingle.cc:  TODO DEscription
 *
 *  Author:
 *   Petar Petrov, MSc
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

#include <Core/Algorithms/Math/TMS/ModelGenericCoilAlgo.h>

namespace SCIRun {

	using namespace SCIRun;

	class ModelTMSCoilSingle : public Module
	{
	  public:

		ModelTMSCoilSingle(GuiContext*);
		
		virtual ~ModelTMSCoilSingle() {}

		virtual void execute();

	  private:
  	    
  	    GuiDouble wireCurrentTCL;
		GuiDouble radiusTCL;
		GuiDouble outerDistanceTCL;
		GuiInt coilDetailsTCL;
    	GuiString typeTCL;
		 
		SCIRunAlgo::ModelTMSCoilSingleAlgo algo;
		SCIRunAlgo::ModelTMSCoilSingleAlgo::Args oldArgs;
	};


	DECLARE_MAKER(ModelTMSCoilSingle)

	ModelTMSCoilSingle::ModelTMSCoilSingle(GuiContext* ctx) :
		Module("ModelTMSCoilSingle", ctx, Source, "TMS", "SCIRun"),
		wireCurrentTCL(ctx->subVar("wireCurrentTCL")),
		radiusTCL(ctx->subVar("coilRadiusTCL")),
		coilDetailsTCL(ctx->subVar("levelDetailTCL")),
		outerDistanceTCL(ctx->subVar("outerDistanceTCL")),
		typeTCL(ctx->subVar("typeTCL")) 
	{
		algo.set_progress_reporter(this);
	}

	void
	ModelTMSCoilSingle::execute()
	{
		SCIRunAlgo::ModelTMSCoilSingleAlgo::Args algoArgs;
		//MatrixHandle omatrix;
		FieldHandle ofield;

		std::string coilType = static_cast<std::string>(typeTCL.get());
		algoArgs.wireCurrent = static_cast<double>(wireCurrentTCL.get());
		algoArgs.coilRadiusInner = static_cast<double>(radiusTCL.get());
		algoArgs.coilRadiusOuter = static_cast<double>(radiusTCL.get());
		algoArgs.coilDistanceOuter = static_cast<double>(outerDistanceTCL.get());
		algoArgs.coilLevelDetails = static_cast<size_t>(coilDetailsTCL.get());
		
		bool need_matrix_data = oport_connected("Matrix");
		bool need_mesh_data = oport_connected("Mesh");

		if (coilType == "single") algoArgs.type = 1;
		else if (coilType == "multi") algoArgs.type = 2;
		else if (coilType == "dipole") algoArgs.type = 3;
		else if (coilType == "test") algoArgs.type = 4;
		else algoArgs.type = 1;
		

		if (!oport_cached("Mesh") || oldArgs != algoArgs)
		{
		    update_state(Executing);
		    
		    if (!(algo.run(ofield,algoArgs))) return;
		  
		    //! send new Mesh output if there is any: 
		    if(need_mesh_data)
		    send_output_handle("Mesh",ofield);
			
			//! send new Matri utput if there is any: 
			// if(need_matrix_data)
			// send_output_handle("Matrix", omatrix);
			
			oldArgs = algoArgs;
		}

	}

} // End namespace SCIRun
