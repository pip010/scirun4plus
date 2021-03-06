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

#include <Core/Algorithms/Fields/TransformMesh/AlignMeshBoundingBoxes.h>

#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Dataflow/Network/Ports/FieldPort.h>

#include <Dataflow/Network/Module.h>

namespace SCIRun {

using namespace SCIRun;

class AlignMeshBoundingBoxes : public Module {
  public:
    AlignMeshBoundingBoxes(GuiContext*);
    virtual ~AlignMeshBoundingBoxes() {}
    virtual void execute();

  private:
    SCIRunAlgo::AlignMeshBoundingBoxesAlgo algo_;
};


DECLARE_MAKER(AlignMeshBoundingBoxes)

AlignMeshBoundingBoxes::AlignMeshBoundingBoxes(GuiContext* ctx) :
  Module("AlignMeshBoundingBoxes", ctx, Source, "ChangeMesh", "SCIRun")
{
  algo_.set_progress_reporter(this);
}


void
AlignMeshBoundingBoxes::execute()
{
  // Get input field.
  FieldHandle ifield, ofield,objfield;
  MatrixHandle omatrix;

  get_input_handle("Input", ifield,true);
  get_input_handle("AlignmentField", objfield,true);

  if (inputs_changed_ || !oport_cached("Output") || !oport_cached("Transform"))
  {
    // Inform module that execution started
    update_state(Executing);

    // Core algorithm of the module
    if(!(algo_.run(ifield,objfield,ofield,omatrix))) return;    

    // Send output to output ports
    send_output_handle("Output", ofield);
    send_output_handle("Transform", omatrix);
  }
}

} // End namespace SCIRun


