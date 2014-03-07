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

// Include all code for the dynamic engine
#include <Core/Datatypes/String.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Parser/ArrayMathEngine.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/StringPort.h>

namespace SCIRun {

class CreateDataArray : public Module {
public:
  CreateDataArray(GuiContext*);
  ~CreateDataArray() {}
  virtual void execute();
  virtual void tcl_command(GuiArgs&, void*);
  
private:
    GuiString guifunction_;
    GuiString guiformat_;
};


DECLARE_MAKER(CreateDataArray)
CreateDataArray::CreateDataArray(GuiContext* ctx)
  : Module("CreateDataArray", ctx, Source, "DataArrayMath", "SCIRun"),
  guifunction_(get_ctx()->subVar("function")),
  guiformat_(get_ctx()->subVar("format"))
{
}

void CreateDataArray::execute()
{
  // Define handles for input
  MatrixHandle size;
  StringHandle func;
  std::vector<MatrixHandle> matrices;
  
  // Get latest handles from ports
  get_input_handle("Function",func,false);
  get_input_handle("Size",size,false);
  get_dynamic_input_handles("Array",matrices,false);

  TCLInterface::eval(get_id()+" update_text");


  // Do something if data changed
  if (inputs_changed_ || guifunction_.changed() || guiformat_.changed() ||
      !oport_cached("Array"))
  {
    update_state(Executing);
      
    size_type n = 1;

    if (size.get_rep())
    {
      if ((size->ncols() != 1)&&(size->nrows()!=1))
      {
        error("Size input needs to be a 1 by 1 matrix");
        return;
      }
      n = static_cast<int>(size->get(0,0));
      if (n == 0) n = 1;
    }

    // Get number of matrix ports with data (the last one is always empty)
    size_t numinputs = matrices.size();
    if (numinputs > 23)
    {
      error("This module cannot handle more than 23 input matrices.");
      return;
    }
      
    NewArrayMathEngine engine;
    engine.set_progress_reporter(this);

    if (func.get_rep())
    {
      guifunction_.set(func->get());
      get_ctx()->reset();
    }
    
    char mname = 'A';
    std::string matrixname("A");
    
    for (size_t p = 0; p < numinputs; p++)
    {
      if (matrices[p].get_rep() == 0)
      {
        error("No matrix was found on input port.");
        return;      
      }

      matrixname[0] = mname++;
      if (!(engine.add_input_matrix(matrixname,matrices[p]))) return;
    }
    
    if(!(engine.add_output_matrix("RESULT",n))) return;    

    // Add an object for getting the index and size of the array.

    if(!(engine.add_index("INDEX"))) return;
    if(!(engine.add_size("SIZE"))) return;

    std::string function = guifunction_.get();
    if(!(engine.add_expressions(function))) return;

    // Actual engine call, which does the dynamic compilation, the creation of the
    // code for all the objects, as well as inserting the function and looping 
    // over every data point

    if (!(engine.run())) return;

    // Get the result from the engine
    MatrixHandle omatrix;    
    engine.get_matrix("RESULT",omatrix);

    // send new output if there is any: 
    send_output_handle("DataArray", omatrix);
  }
}


void
CreateDataArray::tcl_command(GuiArgs& args, void* userdata)
{
  if(args.count() < 2)
  {
    args.error("CreateScalarArray needs a minor command");
    return;
  }

  if( args[1] == "gethelp" )
  {
    return;
  }
  else
  {
    Module::tcl_command(args, userdata);
  }
}


} // End namespace SCIRun


