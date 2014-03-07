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

#include <Core/Algorithms/FieldArray/CreateFieldArray/CreateFieldArrayFromData.h>

#include <Dataflow/Network/Ports/FieldArrayPort.h>
#include <Dataflow/Network/Module.h>

namespace SCIRun {

using namespace SCIRun;

class CreateFieldArrayFromData : public Module {
  public:
    CreateFieldArrayFromData(GuiContext*);
    virtual ~CreateFieldArrayFromData() {}
    virtual void execute();

  private:
    SCIRunAlgo::CreateFieldArrayFromDataAlgo algo_;
    
};


DECLARE_MAKER(CreateFieldArrayFromData)

CreateFieldArrayFromData::CreateFieldArrayFromData(GuiContext* ctx) :
  Module("CreateFieldArrayFromData", ctx, Source, "FieldArray", "SCIRun")
{
  algo_.set_progress_reporter(this);
}

void
CreateFieldArrayFromData::execute()
{
  FieldHandle field;
  MatrixHandle matrix;
  FieldArrayHandle fieldarray;
  
  get_input_handle("Field",field,true);
  get_input_handle("Data",matrix,true);
  
  if (inputs_changed_)
  {
    if(!(algo_.run(field,matrix,fieldarray) )) return;
    send_output_handle("FieldArray",fieldarray);
  }
}

} // End namespace SCIRun


