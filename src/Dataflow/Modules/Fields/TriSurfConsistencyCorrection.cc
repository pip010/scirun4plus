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
#include <Core/Algorithms/Fields/Cleanup/TriSurfConsistencyCorrection.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {

using namespace SCIRun;

class TriSurfConsistencyCorrection : public Module {
public:
  TriSurfConsistencyCorrection(GuiContext*);
  virtual ~TriSurfConsistencyCorrection() {}
  virtual void execute();
  
private:
  GuiInt  output_element_list_;
  
  SCIRunAlgo::TriSurfConsistencyCorrectionAlgo algo_;
};


DECLARE_MAKER(TriSurfConsistencyCorrection)

TriSurfConsistencyCorrection::TriSurfConsistencyCorrection(GuiContext* ctx) :
  Module("TriSurfConsistencyCorrection", ctx, Source, "ChangeMesh", "SCIRun"),
  output_element_list_(get_ctx()->subVar("output_element_list"), 0)
{
  algo_.set_progress_reporter(this);
}

void
TriSurfConsistencyCorrection::execute()
{
  FieldHandle input, output;
  MatrixHandle outputInvertedElementList;

  get_input_handle("Field",input,true);
  
  if (inputs_changed_ ||
      output_element_list_.changed() ||
      ! oport_cached("Field"))
  {
    update_state(Executing);
    
    algo_.set_bool("output_inverted_element_list", output_element_list_.get());

    if (! algo_.run(input, output, outputInvertedElementList) ) return;
    
    send_output_handle("Field", output);
    
    //if (output_element_list_.get() == 0)
    //{
      send_output_handle("InvertedElementsMatrix", outputInvertedElementList);
    //}
  }
}

} // End namespace SCIRun


