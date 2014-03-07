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

#include <Core/Algorithms/FieldArray/MergeMesh/MergeFieldArray.h>
#include <Dataflow/Network/Module.h>

#include <Dataflow/Network/Ports/GeometryPort.h>
#include <Dataflow/Network/Ports/FieldArrayPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {

using namespace SCIRun;

class MergeFieldArray : public Module {
  public:
    MergeFieldArray(GuiContext*);
    virtual ~MergeFieldArray() {}

    virtual void execute();
  
  private:
    GuiDouble guitolerance_;
    GuiInt    guimergenodes_;
    GuiInt    guimatchval_;
    GuiInt    guimeshonly_;
        
    SCIRunAlgo::MergeFieldArrayAlgo algo_;    
        
};


DECLARE_MAKER(MergeFieldArray)

MergeFieldArray::MergeFieldArray(GuiContext* ctx) :
  Module("MergeFieldArray", ctx, Source, "FieldArray", "SCIRun"),
  guitolerance_(get_ctx()->subVar("tolerance"), 0.0001),
  guimergenodes_(get_ctx()->subVar("force-nodemerge"),1),
  guimatchval_(get_ctx()->subVar("matchval"),0),
  guimeshonly_(get_ctx()->subVar("meshonly"),0)  
{
 algo_.set_progress_reporter(this);
}

void
MergeFieldArray::execute()
{
  FieldArrayHandle input;
  FieldHandle output;
  
  get_input_handle("FieldArray",input,true);
  
  if (inputs_changed_ ||  guitolerance_.changed() ||
      guimergenodes_.changed() || guimatchval_.changed() || 
      guimeshonly_.changed() || !oport_cached("Field"))
  {
    algo_.set_bool("merge_nodes",guimergenodes_.get());
    algo_.set_scalar("tolerance",guitolerance_.get());
    algo_.set_bool("match_node_values",guimatchval_.get());
    algo_.set_bool("make_no_data",guimeshonly_.get());
  
    if(!(algo_.run(input,output))) return;
    send_output_handle("Field",output,true);
  }
}

} // End namespace SCIRun


