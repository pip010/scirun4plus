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

#include <Core/Datatypes/FieldArray.h>

#include <Dataflow/Network/Ports/FieldArrayPort.h>
#include <Dataflow/Network/Module.h>


namespace SCIRun {

using namespace SCIRun;

class ExtractFieldArrayRange : public Module {
  public:
    ExtractFieldArrayRange(GuiContext*);

    virtual ~ExtractFieldArrayRange() {}
    virtual void execute();

  private:
    GuiInt  range_start_;
    GuiInt  range_end_;
    GuiInt  range_min_;
    GuiInt  range_max_;
};


DECLARE_MAKER(ExtractFieldArrayRange)

ExtractFieldArrayRange::ExtractFieldArrayRange(GuiContext* ctx) :
  Module("ExtractFieldArrayRange", ctx, Source, "FieldArray", "SCIRun"),
  range_start_(get_ctx()->subVar("range-start"),0),
  range_end_(get_ctx()->subVar("range-end"),10),
  range_min_(get_ctx()->subVar("range-min"),0),
  range_max_(get_ctx()->subVar("range-max"),10)
{
}

void
ExtractFieldArrayRange::execute()
{
  FieldArrayHandle input, output;
  
  get_input_handle("Input",input,true);
  
  if (inputs_changed_ || range_start_.changed() ||
      range_end_.changed() || !oport_cached("Output"))
  {
    std::vector<FieldHandle>& inputs = input->array();
    size_t num_inputs = inputs.size();
    
    range_min_.set(0);
    range_max_.set(num_inputs-1);
    
    int start = range_start_.get();
    int end = range_end_.get();
  
    if (start < 0) start = 0;
    if (start >= num_inputs) start = num_inputs-1;

    if (end < 0) end = 0;
    if (end >= num_inputs) end = num_inputs-1;
    
    if (start > end) { int tmp = start; start = end; end = tmp; }

    range_start_.set(start);
    range_end_.set(end);
    TCLInterface::execute(get_id() + " update_range");
    
    output = new FieldArray;
    std::vector<FieldHandle>& outputs = output->array();
    
    for (size_t j=static_cast<size_t>(start); j<=static_cast<size_t>(end);j++)
    {
      outputs.push_back(inputs[j]);
    }
    
    send_output_handle("Output",output,true);
  }
}

} // End namespace SCIRun


