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


/*
 *  UnuConvert: Convert between C types
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2000
 *
 */

#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

#include <iostream>
using std::endl;
#include <stdio.h>

namespace SCITeem {
using namespace SCIRun;

class UnuConvert : public Module {
  private:
    GuiInt type_;
  public:
    UnuConvert(SCIRun::GuiContext *ctx);
    virtual ~UnuConvert() {}
    virtual void execute();
};

} //end namespace SCITeem

using namespace SCITeem;
DECLARE_MAKER(UnuConvert)

UnuConvert::UnuConvert(SCIRun::GuiContext *ctx)
  : Module("UnuConvert", ctx, Filter, "UnuAtoM", "Teem"), 
  type_(get_ctx()->subVar("type"))
{
}

void 
UnuConvert::execute()
{
  update_state(NeedData);

  NrrdDataHandle nrrdH;
  get_input_handle("Nrrd", nrrdH,true);

  if (inputs_changed_ || type_.changed() || !oport_cached("Nrrd"))
  {
    // Force Teem to be locked befoer calling the Teem library
    NrrdGuard nrrd_guard;

    // Inform module that execution started
    update_state(Executing);
      
    Nrrd *nin = nrrdH->nrrd_;
    Nrrd *nout = nrrdNew();

    if (nrrdConvert(nout, nin, type_.get())) 
    {
      char *err = biffGetDone(NRRD);
      error(std::string("Trouble resampling: ") + err);
      msg_stream_ << "  input Nrrd: nin->dim="<<nin->dim<<"\n";
      free(err);
    }

    NrrdDataHandle onrrdH = new NrrdData(nout);
    send_output_handle("Nrrd", onrrdH, true);
  }
}

