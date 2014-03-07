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
 *  SetWireElectrodePotential: Insert an electrode into a finite element mesh
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   January 2002
 *
 */

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/FieldInformation.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Module.h>

#include <Dataflow/GuiInterface/GuiVar.h>

#include <iostream>
#include <sstream>

namespace BioPSE {

using namespace SCIRun;

class SetWireElectrodePotential : public Module {
  private:
    GuiString active_;
    GuiDouble voltage_;
  public:
    SetWireElectrodePotential(GuiContext *context);
    virtual ~SetWireElectrodePotential() {}
    virtual void execute();
};

DECLARE_MAKER(SetWireElectrodePotential)

SetWireElectrodePotential::SetWireElectrodePotential(GuiContext *context)
  : Module("SetWireElectrodePotential", context, Filter, "Forward", "BioPSE"),
    active_(context->subVar("active")), voltage_(context->subVar("voltage"))
{
}


void
SetWireElectrodePotential::execute()
{
  FieldHandle ielecH;
  get_input_handle("Electrode", ielecH, true);

  FieldHandle elecFld = ielecH->clone();

  FieldInformation fi(elecFld);
  if (!fi.is_curvemesh()) 
  {
    error("Input electrode was not a CurveMesh.");
    return;
  }

  VField *efield = elecFld->vfield();

  double voltage = voltage_.get();
  std::string active = active_.get();
  efield->set_value(voltage,0);
  efield->set_property("active_side", active, false);

  send_output_handle("Electrode", elecFld);
}

} // End namespace BioPSE
