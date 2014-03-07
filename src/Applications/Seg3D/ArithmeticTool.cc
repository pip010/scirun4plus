//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : ArithmeticTool.cc
//    Author : McKay Davis
//    Date   : Sun Oct  1 23:22:04 2006



#include <Applications/Seg3D/Painter.h>
#include <Applications/Seg3D/ArithmeticTool.h>
#include <Core/Events/EventManager.h>
#include <Core/Events/LoadVolumeEvent.h>
#include <sci_gl.h>

#include <Applications/Seg3D/GuiCode/arithmeticvolume.h>
#include <Applications/Seg3D/GuiCode/seg3devents.h>

#ifdef _WIN32
#  undef SCISHARE
#  define SCISHARE __declspec(dllimport)
#else
#  define SCISHARE
#endif

namespace SCIRun {


ArithmeticTool::ArithmeticTool(Painter *painter) : 
  BaseTool("Arithmetic"),
  painter_(painter)
{
  update_to_gui();
}


ArithmeticTool::~ArithmeticTool()
{
}


BaseTool::propagation_state_e 
ArithmeticTool::process_event(event_handle_t event)
{
  if (dynamic_cast<FinishEvent *>(event.get_rep())) {
    finish();
  }

  return CONTINUE_E;
}


void
ArithmeticTool::update_to_gui()
{
  vector< std::string > *names = new vector< std::string >;

  painter_->volume_lock_.lock();
  for (size_t i = 0; i < painter_->volumes_.size(); ++i)
  {
    names->push_back( painter_->volumes_[i]->name_ );
  }
  painter_->volume_lock_.unlock();

  wxCommandEvent event(wxEVT_ARITHMETIC_SET_INPUTS, wxID_ANY);
  event.SetClientData((void *)names);
  wxPanel *panel = Painter::global_seg3dframe_pointer_->CurrentToolPanel();
  if (panel)
  {
    wxPostEvent(panel, event);
  }
}


BaseTool::propagation_state_e
ArithmeticTool::finish()
{
  const string input1 = painter_->get_vars()->get_string("Painter::Arithmetic::input1");
  const string input2 = painter_->get_vars()->get_string("Painter::Arithmetic::input2");
  const string operation = painter_->get_vars()->get_string("Painter::Arithmetic::operation");

  bool first_nrrd_ = true;
  bool second_nrrd_ = true;

  NrrdVolumeHandle nrrd_vol_handle1;
  NrrdVolumeHandle nrrd_vol_handle2;

  if (!(first_nrrd_ = painter_->find_volume_by_name(input1).get_rep()))
  {
    painter_->set_status("Arithmetic tool requires two input volumes - the first is invalid.");
    return STOP_E;
  } else
    nrrd_vol_handle1 = painter_->find_volume_by_name(input1);

  if (!(second_nrrd_ = painter_->find_volume_by_name(input2).get_rep()))
  {
    painter_->set_status("Arithmetic tool requires two input volumes - the second is invalid.");
    return STOP_E;
  } else
    nrrd_vol_handle2 = painter_->find_volume_by_name(input2);

  unsigned int nrrdOp;

  if( (nrrdOp = get_op(operation)) == nrrdBinaryOpUnknown )
  {
    painter_->set_status("Arithmetic tool requires an valid arithmetic operation.");
    return STOP_E;
  }

  Nrrd *nin1 = nrrd_vol_handle1->nrrd_handle_->nrrd_;
  Nrrd *nin2 = nrrd_vol_handle2->nrrd_handle_->nrrd_;

  NrrdIter *in1 = nrrdIterNew();
  NrrdIter *in2 = nrrdIterNew();

  nrrdIterSetOwnNrrd(in1, nin1);
  nrrdIterSetOwnNrrd(in2, nin2);

  NrrdDataHandle nout_data_handle = new NrrdData();
  if (nrrdArithIterBinaryOp(nout_data_handle->nrrd_, nrrdOp, in1, in2))
  {
    char *err = biffGetDone(NRRD);
    painter_->set_status( string("Arithmetic tool error performing nrrd operation failed: ") + err);
    free(err);

    cerr << string("Error on line #") 
	 << to_string(__LINE__)
	 << string(" executing nrrd command: nrrdArithIterBinaryOp \n")
	 << string("Message: ") 
	 << err
	 << std::endl;

    return STOP_E;
  }

  // Labels as well as data can be sent but labels are unsigned chars.
  // however, the final output should be float so covert them to
  // floats.
  nout_data_handle =
    VolumeOps::convert_to(nout_data_handle, nrrdTypeFloat );


  nrrdKeyValueCopy(nout_data_handle->nrrd_, nin2);

  // Add the new volume to the layers.
  const string newname =
    painter_->unique_layer_name(input1 + " ArithOP");
  painter_->make_layer(newname, nout_data_handle, 0);

  // Push this undo item.
  UndoHandle undo = new UndoReplaceLayer(painter_, "Undo Arithmetic",
                                         0, painter_->current_volume_, 0);
  painter_->push_undo(undo);

  // Reset the slices to match the new volume.
  painter_->rebuild_layer_buttons();
  painter_->extract_all_window_slices();
  painter_->redraw_all();

  // Autoview to the arithmetic volume.
  event_handle_t empty;
  painter_->Autoview(empty);

  painter_->hide_tool_panel();

  return CONTINUE_E;
}

unsigned int
ArithmeticTool::get_op(const string &op)
{
  if (op == "+") 
    return nrrdBinaryOpAdd;
  else if (op == "-")
    return nrrdBinaryOpSubtract;
  else if (op == "x")
    return nrrdBinaryOpMultiply;
  else if (op == "/")
    return nrrdBinaryOpDivide;
  else if (op == "^")
    return nrrdBinaryOpPow;
  else if (op == "%")
    return nrrdBinaryOpMod;
  else if (op == "fmod")
    return nrrdBinaryOpFmod;
  else if (op == "atan2")
    return nrrdBinaryOpAtan2;
  else if (op == "min")
    return nrrdBinaryOpMin;
  else if (op == "max")
    return nrrdBinaryOpMax;
  else if (op == "lt")
    return nrrdBinaryOpLT;
  else if (op == "lte")
    return nrrdBinaryOpLTE;
  else if (op == "gt")
    return nrrdBinaryOpGT;
  else if (op == "gte")
    return nrrdBinaryOpGTE;
  else if (op == "eq")
    return nrrdBinaryOpEqual;
  else if (op == "neq")
    return nrrdBinaryOpNotEqual;
  else if (op == "comp")
    return nrrdBinaryOpCompare;
  else if (op == "exists")
    return nrrdBinaryOpExists;
  else {
    return nrrdBinaryOpUnknown;
  }
}

}
