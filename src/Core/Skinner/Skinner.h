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
//    File   : Skinner.h
//    Author : McKay Davis
//    Date   : Fri Jun 30 22:06:40 2006

#ifndef SKINNER_SKINNER_H
#define SKINNER_SKINNER_H

#include <Core/Skinner/Parent.h>
#include <Core/Skinner/Signals.h>
#include <string>

#include <Core/Skinner/share.h>

namespace SCIRun {
namespace Skinner {

class Variables;
class GLWindow;
class Skinner;

SCISHARE Skinner *          load_skin(const string&filename);
SCISHARE Skinner *          load_default_skin(const string& default_path="");


class Skinner : public Parent {
public:
  typedef vector<GLWindow *> GLWindows_t;
  Skinner(Variables *);
  virtual ~Skinner();
  void                      spawn_event_processor_threads();
  GLWindows_t &             get_windows() { return windows_; };

private:
  template <class T>
  Drawable *        construct_class_from_maker_signal(event_handle_t &);

  CatcherFunction_t Animation_Maker;
  CatcherFunction_t Arc_Maker;
  CatcherFunction_t Arithmetic_Maker;
  CatcherFunction_t Arrow_Maker;
  CatcherFunction_t Box_Maker;
  CatcherFunction_t Checkerboard_Maker;
  CatcherFunction_t ColorMap1D_Maker;
  CatcherFunction_t ColorPicker_Maker;
  CatcherFunction_t Collection_Maker;
  CatcherFunction_t FileBrowser_Maker;
  CatcherFunction_t Frame_Maker;
  CatcherFunction_t FocusRegion_Maker;
  CatcherFunction_t FocusGrab_Maker;
  CatcherFunction_t GLWindow_Maker;
  CatcherFunction_t GLWindow_Destructor;
  CatcherFunction_t Gradient_Maker;
  CatcherFunction_t Graph2D_Maker;
  CatcherFunction_t Grid_Maker;
  CatcherFunction_t Histogram_Maker;
  CatcherFunction_t Layout_Maker;
  CatcherFunction_t MenuManager_Maker;
  CatcherFunction_t Progress_Maker;
  CatcherFunction_t RenderBuffer_Maker;
  CatcherFunction_t Slider_Maker;
  CatcherFunction_t SceneGraph_Maker;
  CatcherFunction_t Text_Maker;
  CatcherFunction_t TextEntry_Maker;
  CatcherFunction_t Texture_Maker;
  CatcherFunction_t TransferFunction2D_Maker;
  CatcherFunction_t ViewSubRegion_Maker;
  CatcherFunction_t VisibilityGroup_Maker;
  CatcherFunction_t Volume_Maker;

  CatcherFunction_t Quit;
  CatcherFunction_t Redraw;
  CatcherFunction_t Stop;
  CatcherFunction_t Reload_Default_Skin;

  GLWindows_t       windows_;
  set<string>       window_ids_;
  int               running_windows_;
};


template <class T>
Drawable *
Skinner::construct_class_from_maker_signal(event_handle_t &event) {
  MakerSignal *maker_signal =
    dynamic_cast<MakerSignal *>(event.get_rep());
  ASSERT(maker_signal);
  T *obj = new T(maker_signal->get_vars());
  maker_signal->set_signal_thrower(obj);
  maker_signal->set_signal_name(maker_signal->get_signal_name()+"_Done");
  return obj;
}


}
}

#endif
