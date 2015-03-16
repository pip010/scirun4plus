
/*
 *  ModelTMSCoilDipole.cc:  TODO DEscription
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

	class ModelTMSCoilDipole : public Module
	{
	  public:

		ModelTMSCoilDipole(GuiContext*);
		
		virtual ~ModelTMSCoilDipole() {}

		virtual void execute();

	  private:
  	    
  	    GuiDouble totalCurrentTCL;
  	    GuiInt numberSegmentsTCL;
		GuiDouble innerRadiusTCL;
		GuiDouble outerRadiusTCL;
		GuiDouble outerDistanceTCL;
		GuiInt coilDetailsTCL;
		GuiInt coilLayersTCL;
    	GuiString typeTCL;
		 
		SCIRunAlgo::ModelTMSCoilDipoleAlgo algo;
		SCIRunAlgo::ModelTMSCoilDipoleAlgo::Args oldArgs;
	};


	DECLARE_MAKER(ModelTMSCoilDipole)

	ModelTMSCoilDipole::ModelTMSCoilDipole(GuiContext* ctx) :
		Module("ModelTMSCoilDipole", ctx, Source, "TMS", "SCIRun"),
		totalCurrentTCL(ctx->subVar("totalCurrentTCL")),
		numberSegmentsTCL(ctx->subVar("numberSegmentsTCL")),
		innerRadiusTCL(ctx->subVar("innerRadiusTCL")),
		outerRadiusTCL(ctx->subVar("outerRadiusTCL")),
		outerDistanceTCL(ctx->subVar("outerDistanceTCL")),
		coilLayersTCL(ctx->subVar("coilLayersTCL")),
		coilDetailsTCL(ctx->subVar("levelDetailTCL")),
		typeTCL(ctx->subVar("typeTCL")) 
	{
		algo.set_progress_reporter(this);
	}

	void
	ModelTMSCoilDipole::execute()
	{
		SCIRunAlgo::ModelTMSCoilDipoleAlgo::Args algoArgs;

		FieldHandle ofield;

		std::string coilType = static_cast<std::string>(typeTCL.get());
		algoArgs.totalCurrent = static_cast<double>(totalCurrentTCL.get());
		algoArgs.numberSegments = static_cast<size_t>(numberSegmentsTCL.get());
		algoArgs.coilRadiusInner = static_cast<double>(innerRadiusTCL.get());
		algoArgs.coilRadiusOuter = static_cast<double>(outerRadiusTCL.get());
		algoArgs.coilDistanceOuter = static_cast<double>(outerDistanceTCL.get());
		algoArgs.coilLayers = static_cast<size_t>(coilLayersTCL.get());
		algoArgs.coilLevelDetails = static_cast<size_t>(coilDetailsTCL.get());
		
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
			
			oldArgs = algoArgs;
		}

	}

} // End namespace SCIRun
