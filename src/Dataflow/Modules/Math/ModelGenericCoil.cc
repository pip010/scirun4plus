
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
		GuiDouble coilDistanceTCL;
		GuiDouble coilSegmentsTCL;
    	GuiString typeTCL;
		 
		SCIRunAlgo::ModelGenericCoilAlgo algo;
		SCIRunAlgo::ModelGenericCoilAlgo::Args oldArgs;
	};


	DECLARE_MAKER(ModelGenericCoil)

	ModelGenericCoil::ModelGenericCoil(GuiContext* ctx) :
		Module("ModelGenericCoil", ctx, Source, "Math", "SCIRun"),
		wireCurrentTCL(ctx->subVar("wireCurrentTCL")),
		coilRadiusTCL(ctx->subVar("coilRadiusTCL")),
		coilDistanceTCL(ctx->subVar("coilDistanceTCL")),
		coilSegmentsTCL(ctx->subVar("coilSegmentsTCL")),
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

		std::string coilType = static_cast<std::string>(typeTCL.get());
		algoArgs.wireCurrent = static_cast<double>(wireCurrentTCL.get());
		algoArgs.coilRadius = static_cast<double>(coilRadiusTCL.get());
		algoArgs.coilDistance = static_cast<double>(coilDistanceTCL.get());
		algoArgs.coilSegments = static_cast<double>(coilSegmentsTCL.get());
		
		bool need_matrix_data = oport_connected("Matrix");
		bool need_mesh_data = oport_connected("Mesh");

		if (coilType == "0-shaped") algoArgs.type = 1;
		else if (coilType == "8-shaped") algoArgs.type = 2;
		else algoArgs.type = 1;

		if (!oport_cached("Mesh") || oldArgs != algoArgs)
		{
		    update_state(Executing);
		    
		    if (!(algo.run(ofield,omatrix,algoArgs))) return;
		  
		    //! send new Mesh output if there is any: 
		    if(need_mesh_data)
		    send_output_handle("Mesh",ofield);
			
			//! send new Matri utput if there is any: 
			if(need_matrix_data)
			send_output_handle("Matrix", omatrix);
			
			oldArgs = algoArgs;
		}

	}

} // End namespace SCIRun
