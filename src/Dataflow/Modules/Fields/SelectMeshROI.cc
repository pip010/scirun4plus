/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


//    File   : SelectMeshROI.h
//    Author : Petar Petrov
//    Date   : Oct 2017

/*
#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Mesh.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Core/Thread/CrowdMonitor.h>


#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Widgets/BoxWidget.h>
#include <Core/Datatypes/Clipper.h>
#include <Dataflow/GuiInterface/GuiVar.h>


#include <iostream>
#include <stdio.h>
*/

#include <Core/Algorithms/Fields/ROI/MeshROI.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/Field.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

namespace SCIRun {


class SelectMeshROI : public Module
{
  public:
    SelectMeshROI(GuiContext* ctx);
    virtual ~SelectMeshROI() {}
    virtual void execute();

    //! Fix backwards compatibility
    //virtual void post_read();

  private:
    GuiString gui_select_;
    //GuiString gui_method_;
    GuiInt gui_isoval_;
    
    SCIRunAlgo::RoiMeshAlgo algo_;
};


DECLARE_MAKER(SelectMeshROI)


SelectMeshROI::SelectMeshROI(GuiContext* ctx)
  : Module("SelectMeshROI", ctx, Filter, "MiscField", "SCIRun"),
    gui_select_(get_ctx()->subVar("select"),"point"),
    //gui_method_(get_ctx()->subVar("method"),"default"),
    gui_isoval_(get_ctx()->subVar("isoval"), 1)
{
  algo_.set_progress_reporter(this);
}

void
SelectMeshROI::execute()
{
  // Get input field.
  FieldHandle input;
  MatrixHandle roimat;
  MatrixHandle roidist;

  get_input_handle("Input Field", input,true);
  get_input_handle("Input Matrix", roimat, false);
  get_input_handle("Distance Matrix", roidist, false);
  
  if (inputs_changed_ || 
	gui_select_.changed() ||
	gui_isoval_.changed() ||
	!oport_cached("Output Field") )
  {
    update_state(Executing);
      
    FieldHandle output;
    //MatrixHandle mapping;

	int topo_dist = static_cast<int>(gui_isoval_.get());
    //double topo_dist = gui_isoval_.get();
    

    
    if(roidist.get_rep())
    {
	  if (roidist->get_data_size() > 0) 
      {
		  topo_dist = roidist->get(0,0);
	  }
	}
	

    algo_.set_option("select",gui_select_.get());
    algo_.set_int("distance",topo_dist);
    
   
    //if (gui_method_.get() == "convex")
    //{
	//	algo_.set_bool("hex_convex",true);
	//}

    if(!(algo_.run(input,roimat,output)))
    {
      return;
    }
    
    send_output_handle("Output Field",output);
    //send_output_handle("Mapping",mapping);
  }
}


} // End namespace SCIRun
