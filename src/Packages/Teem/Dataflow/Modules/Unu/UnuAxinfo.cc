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
 *  UnuAxinfo.cc:
 *
 *  Written by:
 *   Darby Van Uitert
 *   Department of Computer Science
 *   University of Utah
 *   January 2004
 *
 */

#include <Dataflow/GuiInterface/GuiVar.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Network/Ports/NrrdPort.h>

namespace SCITeem {
using namespace SCIRun;

class UnuAxinfo : public Module {
public:
  UnuAxinfo(GuiContext* ctx);
  virtual ~UnuAxinfo() {}

  virtual void execute();

  GuiInt              axis_;
  GuiString           label_;
  GuiString           kind_;
  GuiString           center_;
  GuiDouble           min_;
  GuiDouble           max_;
  GuiDouble           spacing_;
  GuiInt              use_label_;
  GuiInt              use_kind_;
  GuiInt              use_center_;
  GuiInt              use_min_;
  GuiInt              use_max_;
  GuiInt              use_spacing_;
private:
  int                 generation_;
};

DECLARE_MAKER(UnuAxinfo)

UnuAxinfo::UnuAxinfo(GuiContext* ctx)
  : Module("UnuAxinfo", ctx, Source, "UnuAtoM", "Teem"),
    axis_(get_ctx()->subVar("axis"), 0),
    label_(get_ctx()->subVar("label"), "---"),
    kind_(get_ctx()->subVar("kind"), "nrrdKindUnknown"),
    center_(get_ctx()->subVar("center"), "nrrdCenterUnknown"),
    min_(get_ctx()->subVar("min"), 0.0),
    max_(get_ctx()->subVar("max"), 1.0),
    spacing_(get_ctx()->subVar("spacing"), 1.0),
    use_label_(get_ctx()->subVar("use_label"), 1),
    use_kind_(get_ctx()->subVar("use_kind"), 1),
    use_center_(get_ctx()->subVar("use_center"), 1),
    use_min_(get_ctx()->subVar("use_min"), 1),
    use_max_(get_ctx()->subVar("use_max"), 1),
    use_spacing_(get_ctx()->subVar("use_spacing"), 1),
    generation_(-1)
{
}

void
UnuAxinfo::execute()
{
  update_state(NeedData);
  
  // The input port (with data) is required.
  NrrdDataHandle nh;
  get_input_handle("Nrrd", nh, true);
  
  if (inputs_changed_ || axis_.changed() || label_.changed() || kind_.changed() ||
      min_.changed() || max_.changed() || spacing_.changed() || 
      use_label_.changed() || use_kind_.changed() ||
      use_center_.changed() || use_min_.changed() ||
      use_max_.changed() || use_spacing_.changed() ||
      !oport_cached("Nrrd"))
  
  {
    // Force Teem to be locked befoer calling the Teem library
    NrrdGuard nrrd_guard;

    // Inform module that execution started
    update_state(Executing);
    
    unsigned int axis = axis_.get();
    if( !nh.get_rep() || generation_ != nh->generation ) 
    {
      if (axis >= nh->nrrd_->dim) 
      {
        error("Please specify an axis within proper range.");
        return;
      }
      
      generation_ = nh->generation;
    }
    
    reset_vars();

    Nrrd *nin = nh->nrrd_;
    Nrrd *nout = nrrdNew();
    
    // copy input nrrd and modify its label, kind, center, min, max and spacing
    if (nrrdCopy(nout, nin)) 
    {
      char *err = biffGetDone(NRRD);
      error(std::string("Trouble copying input nrrd: ") +  err);
      remark("  Input Nrrd: nin->dim=" + to_string(nin->dim));
      free(err);
      return;
    }
    
    if (use_label_.get() && strlen(label_.get().c_str())) 
    {
      nout->axis[axis].label = (char*)airFree(nout->axis[axis].label);
      nout->axis[axis].label = airStrdup(const_cast<char*>(label_.get().c_str()));
    }
    
    if (use_kind_.get()) 
    {
      std::string kind = kind_.get();
      if (kind == "nrrdKindDomain") 
      {
        nout->axis[axis].kind = nrrdKindDomain;
      } 
      else if (kind == "nrrdKindScalar") 
      {
        nout->axis[axis].kind = nrrdKindScalar;
      } 
      else if (kind == "nrrdKind3Color") 
      {
        nout->axis[axis].kind = nrrdKind3Color;
      } 
      else if (kind == "nrrdKind3Vector") 
      {
        nout->axis[axis].kind = nrrdKind3Vector;
      } 
      else if (kind == "nrrdKind3Normal") 
      {
        nout->axis[axis].kind = nrrdKind3Normal;
      } 
      else if (kind == "nrrdKind3DSymMatrix") 
      {
        nout->axis[axis].kind = nrrdKind3DSymMatrix;
      } 
      else if (kind == "nrrdKind3DMaskedSymMatrix") 
      {
        nout->axis[axis].kind = nrrdKind3DMaskedSymMatrix;
      } 
      else if (kind == "nrrdKind3DMatrix") 
      {
        nout->axis[axis].kind = nrrdKind3DMatrix;
      } 
      else if (kind == "nrrdKindList") 
      {
        nout->axis[axis].kind = nrrdKindList;
      } 
      else if (kind == "nrrdKindStub") 
      {
        nout->axis[axis].kind = nrrdKindStub;
      } 
      else 
      {
        nout->axis[axis].kind = nrrdKindUnknown;
      }
    }
      
    if (use_center_.get()) 
    {
      std::string center = center_.get();
      if (center == "nrrdCenterCell") 
      {
        nout->axis[axis].center = nrrdCenterCell;
      } 
      else if (center == "nrrdCenterNode") 
      {
        nout->axis[axis].center = nrrdCenterNode;
      } 
      else 
      {
        nout->axis[axis].center = nrrdCenterUnknown;
      }
    }
      
    if (use_min_.get() && airExists(min_.get())) 
    {
      nout->axis[axis].min = min_.get();
    }
    if (use_max_.get() && airExists(max_.get())) 
    {
      nout->axis[axis].max = max_.get();
    }
    if (use_spacing_.get() && airExists(spacing_.get())) 
    {
      nout->axis[axis].spacing = spacing_.get();
    }
    
    int centering_offset = (nout->axis[axis].center == nrrdCenterNode);
    double calc_max = nout->axis[axis].min + 
      fabs(nout->axis[axis].spacing) * (nout->axis[axis].size-centering_offset);

    if (airExists(nout->axis[axis].min) && 
        airExists(nout->axis[axis].max) &&
        airExists(nout->axis[axis].spacing) && 
        calc_max != nout->axis[axis].max) 
    {
      if (!use_max_.get()) {
	warning("Warning: Output NRRD's min, max and spacing are not valid. Recalculating max.");
	warning("   old Max = " + to_string(nout->axis[axis].max));
	nout->axis[axis].max = calc_max;
	warning("   new Max = " + to_string(nout->axis[axis].max));
      }
      double calc_min = nout->axis[axis].max - 
	  fabs(nout->axis[axis].spacing * (nout->axis[axis].size-centering_offset));
      if (!use_min_.get() && (calc_min != nout->axis[axis].min)) {
	warning("Warning: Output NRRD's min, max and spacing are not valid. Recalculating min.");
	warning("   old Min = " + to_string(nout->axis[axis].min));
	nout->axis[axis].min = calc_min;
	warning("   new Min = " + to_string(nout->axis[axis].min));
      }
      int flip=1;
      if (nout->axis[axis].spacing < 0) {
	flip=-1;
      }
      double calc_spacing = (nout->axis[axis].max - nout->axis[axis].min) / 
	  (nout->axis[axis].size-centering_offset) * flip;
      if (!use_spacing_.get() && (calc_spacing != nout->axis[axis].spacing)) {
	warning("Warning: Output NRRD's min, max and spacing are not valid. Recalculating spacing.");
	warning("   old Spacing = " + to_string(nout->axis[axis].spacing));
	nout->axis[axis].spacing = calc_spacing;
	warning("   new Spacing = " + to_string(nout->axis[axis].spacing));
      }
    }
	 
    NrrdDataHandle out(new NrrdData(nout));
    
    // Copy the properties.
    out->copy_properties(nh.get_rep());
    
    send_output_handle("Nrrd", out);
  }
}


} // end SCITeem namespace

