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

//! Include the algorithm
#include <Core/Algorithms/Fields/FieldData/GetFieldData.h>
#include <Core/Algorithms/Math/BiotSavartSolver/BiotSavartSolver.h>

//! The module class
#include <Dataflow/Network/Module.h>

//! We need to define the ports used
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

namespace SCIRun {

class SolveBiotSavartContour : public Module {
  public:
    //! constructor and execute function
    SolveBiotSavartContour(GuiContext*);
    virtual ~SolveBiotSavartContour() {}
    virtual void execute();
  
  private:
    //! Define algorithms needed
    //SCIRunAlgo::GetFieldDataAlgo algo_;
	SCIRunAlgo::BiotSavartSolverAlgo algo_;
};


DECLARE_MAKER(SolveBiotSavartContour)
SolveBiotSavartContour::SolveBiotSavartContour(GuiContext* ctx)
  : Module("SolveBiotSavartContour", ctx, Source, "Math", "SCIRun")
{
  //! Forward error messages;
  algo_.set_progress_reporter(this);
}

void SolveBiotSavartContour::execute()
{
  //! Define dataflow handles:
  FieldHandle meshField;
  FieldHandle coilField;
  MatrixHandle matrixdata(0);
  
  //! Get data from port:
  if(!(get_input_handle("Mesh",meshField,true))) return;
  if(!(get_input_handle("Coil",coilField,true))) return;

  //! Data is only computed if the output port is connected:
  bool need_matrix_data = oport_connected("BField");

  //! Only do work if needed:
  if (inputs_changed_ ||
      (!oport_cached("BField") && need_matrix_data))
  {    
    update_state(Executing);
    
    if( need_matrix_data) 
    {
      //! Run algorithm
      if(!(algo_.run(meshField,matrixdata))) return;
    }

    //! If port is not connected at time of execute, send down a null handle
    //! send data downstream:
    send_output_handle("BField", matrixdata);

  }
}

} // End namespace SCIRun

