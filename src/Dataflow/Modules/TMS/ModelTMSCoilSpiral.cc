
/*
 *  ModelTMSCoilSpiral.cc:  TODO DEscription
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

	class ModelTMSCoilSpiral : public Module
	{
	  public:

		ModelTMSCoilSpiral(GuiContext*);
		
		virtual ~ModelTMSCoilSpiral() {}

		virtual void execute();

	  private:
  	    
  	    GuiDouble wireCurrentTCL;
  	    GuiInt wireLoopsTCL;
		GuiDouble innerRadiusTCL;
		GuiDouble outerRadiusTCL;
		GuiDouble outerDistanceTCL;
		GuiInt coilDetailsTCL;
    	GuiString typeTCL;
		 
		SCIRunAlgo::ModelTMSCoilSpiralAlgo algo;
		SCIRunAlgo::ModelTMSCoilSpiralAlgo::Args oldArgs;
	};


	DECLARE_MAKER(ModelTMSCoilSpiral)

	ModelTMSCoilSpiral::ModelTMSCoilSpiral(GuiContext* ctx) :
		Module("ModelTMSCoilSpiral", ctx, Source, "TMS", "SCIRun"),
		wireCurrentTCL(ctx->subVar("wireCurrentTCL")),
		wireLoopsTCL(ctx->subVar("wireLoopsTCL")),
		innerRadiusTCL(ctx->subVar("innerRadiusTCL")),
		outerRadiusTCL(ctx->subVar("outerRadiusTCL")),
		outerDistanceTCL(ctx->subVar("outerDistanceTCL")),
		coilDetailsTCL(ctx->subVar("levelDetailTCL")),
		typeTCL(ctx->subVar("typeTCL")) 
	{
		algo.set_progress_reporter(this);
	}

	void
	ModelTMSCoilSpiral::execute()
	{
		SCIRunAlgo::ModelTMSCoilSpiralAlgo::Args algoArgs;
		//MatrixHandle omatrix;
		FieldHandle ofield;

		std::string coilType = static_cast<std::string>(typeTCL.get());
		algoArgs.wireCurrent = static_cast<double>(wireCurrentTCL.get());
		algoArgs.wireLoops = static_cast<size_t>(wireLoopsTCL.get());
		algoArgs.coilRadiusInner = static_cast<double>(innerRadiusTCL.get());
		algoArgs.coilRadiusOuter = static_cast<double>(outerRadiusTCL.get());
		algoArgs.coilDistanceOuter = static_cast<double>(outerDistanceTCL.get());
		algoArgs.coilLevelDetails = static_cast<size_t>(coilDetailsTCL.get());
		
		bool need_matrix_data = oport_connected("Matrix");
		bool need_mesh_data = oport_connected("Mesh");

		if (coilType == "0-shape") algoArgs.type = 1;
		else if (coilType == "8-shape") algoArgs.type = 2;
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
