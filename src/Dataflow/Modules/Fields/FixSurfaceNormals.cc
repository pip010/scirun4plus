/*
 For more information, please see: http://software.sci.utah.edu
 
 The MIT License
 
 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

#include <Core/Datatypes/Field.h>
#include <Core/Algorithms/Fields/Cleanup/FixSurfaceNormals.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {

using namespace SCIRun;

class FixSurfaceNormals : public Module {
public:
  FixSurfaceNormals(GuiContext*);
  virtual ~FixSurfaceNormals() {}
  virtual void execute();
  
private:
  GuiInt  gui_output_inverted;
  GuiInt  gui_seed;
  
  SCIRunAlgo::FixSurfaceNormalsAlgo algo_;
};


DECLARE_MAKER(FixSurfaceNormals)

FixSurfaceNormals::FixSurfaceNormals(GuiContext* ctx) :
  Module("FixSurfaceNormals", ctx, Source, "ChangeMesh", "SCIRun"),
  gui_output_inverted(get_ctx()->subVar("output_inverted"), 0),
  gui_seed(get_ctx()->subVar("seed"), 0)
{
  algo_.set_progress_reporter(this);
}

void
FixSurfaceNormals::execute()
{
  FieldHandle input, output;
  MatrixHandle seed_matrix;

  get_input_handle("Field",input,true);
  get_input_handle("Seed Matrix",seed_matrix,false);
  
  if (inputs_changed_ ||
      gui_output_inverted.changed() ||
      gui_seed.changed() ||
      ! oport_cached("Field"))
  {
    update_state(Executing);
    
    int seed = static_cast<int>(gui_seed.get());
    
	if(seed_matrix.get_rep())
    {
	  if (seed_matrix->get_data_size() > 0) 
      {
		  seed = seed_matrix->get(0,0);
	  }
	}
    
    algo_.set_bool("gui_output_inverted", gui_output_inverted.get());
    algo_.set_int("gui_seed",seed);

    if (! algo_.run(input, output) ) return;
    
    send_output_handle("Field", output);
    
  }
}

} // End namespace SCIRun


