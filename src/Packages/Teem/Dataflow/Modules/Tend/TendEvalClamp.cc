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

//    File   : TendEvalClamp.cc
//    Author : Martin Cole
//    Date   : Mon Sep  8 09:46:49 2003

#include <Core/Util/StringUtil.h>

#include <Dataflow/Network/Module.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

#include <teem/ten.h>

#include <sstream>
#include <iostream>
#include <stdio.h>

namespace SCITeem {

using namespace SCIRun;

class TendEvalClamp : public Module {
  public:
    TendEvalClamp(SCIRun::GuiContext *ctx);
    virtual ~TendEvalClamp() {}
    virtual void execute();

  private:
    GuiString       min_;
    GuiString       max_;
};

DECLARE_MAKER(TendEvalClamp)

TendEvalClamp::TendEvalClamp(SCIRun::GuiContext *ctx) : 
  Module("TendEvalClamp", ctx, Filter, "Tend", "Teem"), 
  min_(get_ctx()->subVar("min"), "0.0001"),
  max_(get_ctx()->subVar("max"), "NaN")
{
}


void 
TendEvalClamp::execute()
{
  update_state(NeedData);

  NrrdDataHandle nrrd_handle;
  get_input_handle("nin", nrrd_handle);

  if (inputs_changed_ || min_.changed() || max_.changed() ||
      !oport_cached("nout"))
  {
    // Force Teem to be locked befoer calling the Teem library
    NrrdGuard nrrd_guard;


    // Inform module that execution started
    update_state(Executing);

    Nrrd *nin = nrrd_handle->nrrd_;
    Nrrd *nout = nrrdNew();
    
    float min, max;
    min=max=AIR_NAN;

    if (min_.get() != "NaN" && min_.get() != "nan") 
      from_string(min_.get(),min);
    if (max_.get() != "NaN" && max_.get() != "nan")
      from_string(max_.get(),max);

    if (tenEigenvalueClamp(nout, nin, min, max)) 
    {
      char *err = biffGetDone(TEN);
      error(std::string("Error making tendEvalClamp volume: ") + err);
      free(err);
      return;
    }

    NrrdDataHandle ntmp(new NrrdData(nout));

    send_output_handle("nout", ntmp);
  }
}

} // End namespace SCITeem
