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
//    File   : Skinner.cc
//    Author : McKay Davis
//    Date   : Fri Jun 30 22:10:07 2006



#include <Core/Skinner/Animation.h>
#include <Core/Skinner/Arithmetic.h>
#include <Core/Skinner/Arrow.h>
#include <Core/Skinner/Arc.h>
#include <Core/Skinner/Box.h>
#include <Core/Skinner/Checkerboard.h>
#include <Core/Skinner/Collection.h>
#include <Core/Skinner/ColorMap1D.h>
#include <Core/Skinner/ColorPicker.h>
#include <Core/Skinner/EventProcessor.h>
#include <Core/Skinner/FileBrowser.h>
#include <Core/Skinner/FocusGrab.h>
#include <Core/Skinner/FocusRegion.h>
#include <Core/Skinner/Frame.h>
#include <Core/Skinner/Gradient.h>
#include <Core/Skinner/Graph2D.h>
#include <Core/Skinner/Grid.h>
#include <Core/Skinner/Histogram.h>
#include <Core/Skinner/Layout.h>
#include <Core/Skinner/MenuManager.h>
#include <Core/Skinner/Progress.h>
#include <Core/Skinner/RenderBuffer.h>
#include <Core/Skinner/SceneGraph.h>
#include <Core/Skinner/Skinner.h>
#include <Core/Skinner/Slider.h>
#include <Core/Skinner/Text.h>
#include <Core/Skinner/TextEntry.h>
#include <Core/Skinner/Texture.h>
#include <Core/Skinner/TransferFunction2D.h>
#include <Core/Skinner/Variables.h>
#include <Core/Skinner/ViewSubRegion.h>
#include <Core/Skinner/VisibilityGroup.h>
#include <Core/Skinner/Volume.h>
#include <Core/Skinner/Window.h>
#include <Core/Skinner/XMLIO.h>

#include <Core/Skinner/FilterRedrawEventsTool.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/Environment.h>
#include <Core/Util/StringUtil.h>
#include <iostream>

#include <algorithm>

using std::cerr;
using std::endl;

namespace SCIRun {
namespace Skinner {


Skinner *
load_skin(const string &filename)
{
#if !defined(DEBUG)
  try {
#endif
    Skinner *skinner = XMLIO::load(filename);
    //        }
    EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
    skinner->spawn_event_processor_threads();
    return skinner;
#if !defined(DEBUG)
  } catch (const string &error) {
    cerr << "Skinner Error: " << error << std::endl;
  } catch (const char *&error) {
    cerr << "Skinner Error: " << error << std::endl;
  } catch (...) {
    cerr << "UNKNOWN Skinner Error" << std::endl;
  }
#endif

  return 0;
}


Skinner *
load_default_skin(const string& default_path /* = "" */)
{
  string scirun_objdir_skin = sci_getenv("SCIRUN_OBJDIR")+string("data");
  string skinner_path;

  if (default_path != "") {
    skinner_path = default_path + ":" + scirun_objdir_skin;
  }
  else {
    skinner_path = scirun_objdir_skin;
  }

  const char *path_ptr = sci_getenv("SKINNER_PATH");
  if (path_ptr) {
    skinner_path = string(path_ptr) + ":" + skinner_path;
  }
  sci_putenv("SKINNER_PATH", skinner_path);
  sci_putenv("SCIRUN_FONT_PATH",skinner_path);
  const char *fn = sci_getenv("SKINNER_FILE");
  string filename = (fn ? fn : "main.skin");

  string path = findFileInPath(filename, skinner_path);
  if (path.empty()) {
    std::cerr << "Cannot find main.skin in SKINNER_PATH.\n";
    std::cerr << "SKINNER_PATH=" << skinner_path << std::endl;;
    return 0;
  }

  return load_skin(path+filename);
}



Skinner::Skinner(Variables *variables) :
  Parent(variables),
  windows_(),
  window_ids_(),
  running_windows_(0)
{
  REGISTER_CATCHER_TARGET(Skinner::Redraw);
  REGISTER_CATCHER_TARGET(Skinner::Animation_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Arc_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Box_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Arithmetic_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Arrow_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Checkerboard_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Collection_Maker);
  REGISTER_CATCHER_TARGET(Skinner::ColorMap1D_Maker);
  REGISTER_CATCHER_TARGET(Skinner::ColorPicker_Maker);
  REGISTER_CATCHER_TARGET(Skinner::GLWindow_Maker);
  REGISTER_CATCHER_TARGET(Skinner::GLWindow_Destructor);
  REGISTER_CATCHER_TARGET(Skinner::Grid_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Gradient_Maker);
  REGISTER_CATCHER_TARGET(Skinner::FileBrowser_Maker);
  REGISTER_CATCHER_TARGET(Skinner::FocusGrab_Maker);
  REGISTER_CATCHER_TARGET(Skinner::FocusRegion_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Frame_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Graph2D_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Histogram_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Layout_Maker);
  REGISTER_CATCHER_TARGET(Skinner::MenuManager_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Progress_Maker);
  REGISTER_CATCHER_TARGET(Skinner::RenderBuffer_Maker);
  REGISTER_CATCHER_TARGET(Skinner::SceneGraph_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Slider_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Text_Maker);
  REGISTER_CATCHER_TARGET(Skinner::TextEntry_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Texture_Maker);
  REGISTER_CATCHER_TARGET(Skinner::TransferFunction2D_Maker);
  REGISTER_CATCHER_TARGET(Skinner::ViewSubRegion_Maker);
  REGISTER_CATCHER_TARGET(Skinner::VisibilityGroup_Maker);
  REGISTER_CATCHER_TARGET(Skinner::Volume_Maker);

  REGISTER_CATCHER_TARGET(Skinner::Stop);
  REGISTER_CATCHER_TARGET(Skinner::Reload_Default_Skin);
  REGISTER_CATCHER_TARGET_BY_NAME(Quit, Skinner::Quit);
}


Skinner::~Skinner() {
}


DECLARE_SKINNER_MAKER(Skinner, Animation);
DECLARE_SKINNER_MAKER(Skinner, Arc);
DECLARE_SKINNER_MAKER(Skinner, Box);
DECLARE_SKINNER_MAKER(Skinner, Arithmetic);
DECLARE_SKINNER_MAKER(Skinner, Arrow);
DECLARE_SKINNER_MAKER(Skinner, Checkerboard);
DECLARE_SKINNER_MAKER(Skinner, Collection);
DECLARE_SKINNER_MAKER(Skinner, ColorMap1D);
DECLARE_SKINNER_MAKER(Skinner, ColorPicker);
DECLARE_SKINNER_MAKER(Skinner, FileBrowser);
DECLARE_SKINNER_MAKER(Skinner, FocusGrab);
DECLARE_SKINNER_MAKER(Skinner, FocusRegion);
DECLARE_SKINNER_MAKER(Skinner, Frame);
DECLARE_SKINNER_MAKER(Skinner, Gradient);
DECLARE_SKINNER_MAKER(Skinner, Graph2D);
DECLARE_SKINNER_MAKER(Skinner, Grid);
DECLARE_SKINNER_MAKER(Skinner, Histogram);
DECLARE_SKINNER_MAKER(Skinner, Layout);
DECLARE_SKINNER_MAKER(Skinner, MenuManager);
DECLARE_SKINNER_MAKER(Skinner, Progress);
DECLARE_SKINNER_MAKER(Skinner, RenderBuffer);
DECLARE_SKINNER_MAKER(Skinner, SceneGraph);
DECLARE_SKINNER_MAKER(Skinner, Slider);
DECLARE_SKINNER_MAKER(Skinner, Text);
DECLARE_SKINNER_MAKER(Skinner, TextEntry);
DECLARE_SKINNER_MAKER(Skinner, Texture);
DECLARE_SKINNER_MAKER(Skinner, TransferFunction2D);
DECLARE_SKINNER_MAKER(Skinner, ViewSubRegion);
DECLARE_SKINNER_MAKER(Skinner, VisibilityGroup);
DECLARE_SKINNER_MAKER(Skinner, Volume);


BaseTool::propagation_state_e
Skinner::GLWindow_Maker(event_handle_t &event)
{
  Variables *vars =
    dynamic_cast<MakerSignal *>(event.get_rep())->get_vars();

  string id = vars->get_id();
  int num = 1;
  while (window_ids_.find(id) != window_ids_.end()) {
    id = vars->get_id() + "-" + to_string(num++);
  }
  window_ids_.insert(id);
  Var<string> win_id(vars,"id");
  win_id = id;


  GLWindow *window = dynamic_cast<GLWindow *>
    (construct_class_from_maker_signal<GLWindow>(event));

  ASSERT(window);
  windows_.push_back(window);
  return STOP_E;
}


BaseTool::propagation_state_e
Skinner::GLWindow_Destructor(event_handle_t &event)
{
  Signal *signal = dynamic_cast<Signal *>(event.get_rep());
  ASSERT(signal);

  GLWindow *window =
    dynamic_cast<GLWindow *>(signal->get_signal_thrower());
  ASSERT(window);

  GLWindows_t::iterator iter =
    find(windows_.begin(), windows_.end(), window);
  ASSERT(iter != windows_.end());

  windows_.erase(iter);

  Drawables_t::iterator citer =
    find(children_.begin(), children_.end(), window);
  if (citer != children_.end()) {
    children_.erase(citer);
  }

  return MODIFIED_E;
}


BaseTool::propagation_state_e
Skinner::Quit(event_handle_t& /*event*/)
{
  EventManager::add_event(new QuitEvent());
  return STOP_E;
}


BaseTool::propagation_state_e
Skinner::Redraw(event_handle_t& /*event*/)
{
  EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E));
  return CONTINUE_E;
}


BaseTool::propagation_state_e
Skinner::Stop(event_handle_t &)
{
  return STOP_E;
}


BaseTool::propagation_state_e
Skinner::Reload_Default_Skin(event_handle_t &)
{
  load_default_skin();
  return CONTINUE_E;
}


void
Skinner::spawn_event_processor_threads()
{
  for (unsigned int w = running_windows_; w < windows_.size(); ++w) {
    GLWindow *window = windows_[w];
    EventProcessor *processor = new EventProcessor(window);
    string name = window->get_id() + " Event Processor";
    Thread *thread = new Thread(processor, name.c_str());
    thread->detach();
  }
  running_windows_ = windows_.size();
}


}
}

