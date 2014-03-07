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
//    CreateAndEditColorMap2D
//    File   : CreateAndEditColorMap2D.cc
//    Author : Milan Ikits
//    Author : Michael Callahan
//    Date   : Thu Jul  8 01:50:58 2004

#include <Dataflow/Network/Module.h>
#include <Dataflow/GuiInterface/GuiVar.h>
#include <Core/Util/StringUtil.h>
#include <Core/Containers/Array3.h>
#include <Core/Persistent/Pstreams.h>
#include <Core/Events/CM2View/CM2View.h>
#include <Core/Events/CursorChangeEvent.h>
#include <Core/Events/ToolManager.h>
#include <Core/Events/BaseTool.h>
#include <Core/Events/CM2View/CM2UndoTool.h>

#include <Core/Volume/CM2Widget.h>
#include <Dataflow/GuiInterface/TkOpenGLContext.h>

#include <Dataflow/Network/Ports/ColorMap2Port.h>
#include <Dataflow/Network/Ports/ColorMapPort.h>
#include <Dataflow/Network/Ports/NrrdPort.h>
#include <Dataflow/Comm/MessageBase.h>
#include <Dataflow/Comm/MessageTypes.h>
#include <Dataflow/Network/Scheduler.h>

#include <set>
#include <sstream>
#include <iostream>

namespace SCIRun {

class CreateAndEditColorMap2D;

class HandleCursorChangeTool : public BaseTool
{
public:
  HandleCursorChangeTool(const std::string& name, CreateAndEditColorMap2D *m) :
    BaseTool(name),
    mod_(m)
  {}
  virtual ~HandleCursorChangeTool() {}

  propagation_state_e process_event(event_handle_t event);

private:
  CreateAndEditColorMap2D *mod_;

};

//! Tools for widget manipulation and event handling.
class WidgetsChangedTool : public BaseTool
{
public:
  WidgetsChangedTool(const std::string& name, CreateAndEditColorMap2D *m) :
    BaseTool(name),
    mod_(m)
  {}
  virtual ~WidgetsChangedTool() {}

  virtual propagation_state_e process_event(event_handle_t event);
private:
    CreateAndEditColorMap2D *mod_;
};

class EventHandler : public Runnable 
{
  public:
    EventHandler(CreateAndEditColorMap2D *m) :
      mod_(m)
    {}
    virtual void          run();
  private:
    CreateAndEditColorMap2D *mod_;
};

class CreateAndEditColorMap2D : public Module {
  friend class EventHandler;
  friend class WidgetsChangedTool;
public:
  CreateAndEditColorMap2D(GuiContext* ctx);
  virtual ~CreateAndEditColorMap2D();
  
  // Virtually Inhereted from Dataflow/Modules/Module.h
  virtual void			execute();
  virtual void			tcl_command(GuiArgs&, void*);
  virtual void			presave();
  // get widgets loaded from an srn set up.
  virtual void                  post_read();
  void				set_window_cursor(const std::string&);

private:
  void        wlock() { wid_mutex_.lock(); }
  void        wunlock() { wid_mutex_.unlock(); }
  void				force_execute();
  void				resize_gui(int n = -1);
  void				update_from_gui();
  void				update_to_gui(bool forward = true);
  void				tcl_unpickle();
  void				save_file(bool save_ppm=false);
  void				load_file();
  void				draw_texture(GLuint &);

  bool				select_widget(int w=-1, int o=-1);
  // functions for changing state via GUI
  void				gui_color_change(GuiArgs &args);
  void				gui_toggle_change(GuiArgs &args);
  void        gui_shade_toggle_change(GuiArgs &args);
  //! send widget state event to view.
  void                         update_view_widgets();
  void                         clear_widgets();
  void                         set_widgets(const std::vector<SLIVR::CM2Widget*> &w, bool updating);
  void destroy_gl_view();

  // Input/Output Ports
  bool dead_;
  CM2View* cm2view_;
  std::string cm2_target_;
  ColorMap2Handle		sent_cmap2_;
  ColorMap2Handle		merge_cmap2_;

  bool				force_execute_;
  std::vector<SLIVR::CM2Widget*>	widgets_;

  int				icmap_generation_;
  int				hist_generation_;

  bool				updating_;
  
  bool				need_to_update_gui_;
  GuiDouble			gui_histo_;
  NrrdDataHandle                histo_;

  GuiInt			gui_selected_widget_;
  // The currently selected widgets selected object
  GuiInt			gui_selected_object_;

  GuiInt			gui_num_entries_;
  std::vector<GuiString *>		gui_name_;
  std::vector<GuiDouble *>		gui_color_r_;
  std::vector<GuiDouble *>		gui_color_g_;
  std::vector<GuiDouble *>		gui_color_b_;
  std::vector<GuiDouble *>		gui_color_a_;
  std::vector<GuiString *>		gui_wstate_;
  std::vector<GuiInt *>		gui_onstate_;
  std::vector<GuiInt *>		gui_shade_type_;

  // variables for file loading and saving
  GuiFilename			filename_;
  GuiString *			end_marker_;

  //! We get events from the CM2View to update widgets.
  ToolManager                      tm_;
  EventManager::event_mailbox_t   *events_;
  RecursiveMutex                   wid_mutex_;
};

void          
EventHandler::run() 
{
  for (;;) 
  {
    event_handle_t e = mod_->events_->receive();
    QuitEvent *quit = dynamic_cast<QuitEvent*>(e.get_rep());
    if (quit)
    {
      quit->signal_done();
      break;
    }
    mod_->tm_.propagate_event(e);
  }
}

BaseTool::propagation_state_e 
HandleCursorChangeTool::process_event(event_handle_t event)
{
  CursorChangeEvent *cce = 
    dynamic_cast<CursorChangeEvent*>(event.get_rep());
  if (! cce) {
    return CONTINUE_E;
  }
  mod_->set_window_cursor(cce->get_cursor_name());
  return STOP_E;
}

BaseTool::propagation_state_e 
WidgetsChangedTool::process_event(event_handle_t event)
{
  SetWidgetsEvent *swe = dynamic_cast<SetWidgetsEvent*>(event.get_rep());
  if (! swe) {
    return CONTINUE_E;
  }
  mod_->wlock();
  mod_->set_widgets(swe->widgets_, swe->updating_);
  mod_->wunlock();
  return STOP_E;
}

DECLARE_MAKER(CreateAndEditColorMap2D)

CreateAndEditColorMap2D::CreateAndEditColorMap2D(GuiContext* ctx) : 
  Module("CreateAndEditColorMap2D", ctx, Filter, "Visualization", "SCIRun"),
  dead_(false),
  cm2view_(0),
  cm2_target_("CM2View:" + get_ctx()->getfullname()),
  sent_cmap2_(0),
  merge_cmap2_(0),
  force_execute_(false),
  widgets_(),
  icmap_generation_(-1),
  hist_generation_(-1),
  updating_(false),
  need_to_update_gui_(true),
  gui_histo_(get_ctx()->subVar("histo"), 0.5),
  histo_(0),
  gui_selected_widget_(get_ctx()->subVar("selected_widget"), -1),
  gui_selected_object_(get_ctx()->subVar("selected_object"), -1),
  gui_num_entries_(get_ctx()->subVar("num-entries"), 1),
  gui_name_(),
  gui_color_r_(),
  gui_color_g_(),
  gui_color_b_(),
  gui_color_a_(),
  gui_wstate_(),
  gui_onstate_(),
  gui_shade_type_(),
  filename_(get_ctx()->subVar("filename"), "MyTransferFunction"),
  end_marker_(0),
  tm_("CAECM2 tool manager"),
  wid_mutex_("CreateAndEditColorMap2D widget mutex")
{
  widgets_.push_back(new RectangleCM2Widget());
  select_widget(widgets_.size()-1, 1);

  // Set initial var state for widget storage
  gui_name_.push_back(new GuiString(get_ctx()->subVar("name-0")));
  gui_color_r_.push_back(new GuiDouble(get_ctx()->subVar("0-color-r")));
  gui_color_g_.push_back(new GuiDouble(get_ctx()->subVar("0-color-g")));
  gui_color_b_.push_back(new GuiDouble(get_ctx()->subVar("0-color-b")));
  gui_color_a_.push_back(new GuiDouble(get_ctx()->subVar("0-color-a")));
  gui_wstate_.push_back(new GuiString(get_ctx()->subVar("state-0")));
  gui_onstate_.push_back(new GuiInt(get_ctx()->subVar("on-0")));
  gui_shade_type_.push_back(new GuiInt(get_ctx()->subVar("shadeType-0")));

  events_ = EventManager::register_event_messages(get_ctx()->getfullname());

  tool_handle_t cct = new HandleCursorChangeTool("HandleCursorChangeTool", 
						 this);
  tm_.add_tool(cct, 100);

  tool_handle_t wct = new WidgetsChangedTool("WidgetsChangedTool", this);
  tm_.add_tool(wct, 101);


  // launch the event handler.
  std::string tid = "EventHandler_" + get_ctx()->getfullname();
  Thread *vt = new Thread(new EventHandler(this), tid.c_str());
  vt->detach(); // runs until thread exits.
}



CreateAndEditColorMap2D::~CreateAndEditColorMap2D()
{
  // destroy the event handler thread.
  event_handle_t event = new QuitEvent(get_ctx()->getfullname());
  EventManager::add_event(event);
  event->wait_signal();

  // unregister.
  EventManager::unregister_event_messages(get_ctx()->getfullname());

  destroy_gl_view();

  // unregister.
  EventManager::unregister_event_messages(cm2_target_);
  
  delete end_marker_;
}

void
CreateAndEditColorMap2D::force_execute()
{
  force_execute_ = true;
  want_to_execute();
}

void
CreateAndEditColorMap2D::tcl_command(GuiArgs& args, void* userdata)
{
  if (args.count() < 2) 
  {
    args.error("No command for EditTransferFunc");
    return;
  }

  if (args[1] == "addtriangle") 
  {
    if (cm2view_)
    {  
      CommandEvent *c = new CommandEvent(cm2_target_);
      c->set_command(CM2View::Add_Triangle_Widget_c);
      event_handle_t event = c;
      EventManager::add_event(event);
    }
  } 
  else if (args[1] == "addrectangle") 
  {
    if (cm2view_)
    {  
      CommandEvent *c = new CommandEvent(cm2_target_);
      c->set_command(CM2View::Add_Rectangle_Widget_c);
      event_handle_t event = c;
      EventManager::add_event(event);
    }
  } 
  else if (args[1] == "addellipsoid") 
  {
    if (cm2view_)
    {  
      CommandEvent *c = new CommandEvent(cm2_target_);
      c->set_command(CM2View::Add_Ellipsoid_Widget_c);
      event_handle_t event = c;
      EventManager::add_event(event);
    }
  } 
  else if (args[1] == "addparaboloid") 
  {
    if (cm2view_)
    {    
      CommandEvent *c = new CommandEvent(cm2_target_);
      c->set_command(CM2View::Add_Paraboloid_Widget_c);
      event_handle_t event = c;
      EventManager::add_event(event);
    }
  } 
  else if (args[1] == "deletewidget") 
  {
    if (cm2view_)
    {  
      CommandEvent *c = new CommandEvent(cm2_target_);
      c->set_command(CM2View::Del_Selected_Widget_c);
      event_handle_t event = c;
      EventManager::add_event(event);
    }
  } 
  else if (args[1] == "undowidget") 
  {
    if (cm2view_)
    {  
      event_handle_t event = new CM2UndoEvent(cm2_target_);
      EventManager::add_event(event);
    }
  }
  else if (args[1] == "unpickle") 
  { 
    // tcl_unpickle(); 
  }
  else if (args[1] == "load") load_file();
  else if (args[1] == "save") save_file((args.count() > 2));
  else if (args[1] == "toggle") gui_toggle_change(args);
  else if (args[1] == "shade_toggle") gui_shade_toggle_change(args);
  else if (args[1] == "color") gui_color_change(args);
  else if (args[1] == "select_widget") 
  {
    select_widget();
    update_view_widgets();
  } 
  else if (args[1] == "mouse") 
  {
    if (cm2view_)
    {  
      PointerEvent *p = new PointerEvent();
      p->set_target(cm2_target_);
      p->set_x(args.get_int(3));
      p->set_y(args.get_int(4));
      unsigned int state = 0;
      unsigned int modifiers = 0;

      if (args[2] == "motion") {
        state = PointerEvent::MOTION_E;
      } else if (args[2] == "widget_motion") {
        state = PointerEvent::MOTION_E | PointerEvent::BUTTON_1_E;
      } else if (args[2] == "push") {
        state = PointerEvent::BUTTON_PRESS_E | PointerEvent::BUTTON_1_E;
      } else if (args[2] == "release") {
        state = PointerEvent::BUTTON_RELEASE_E | PointerEvent::BUTTON_1_E;
      } else if (args[2] == "x_late_start") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::BUTTON_PRESS_E | PointerEvent::BUTTON_1_E;
      } else if (args[2] == "x_late_motion") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::MOTION_E | PointerEvent::BUTTON_1_E;
      } else if (args[2] == "x_late_end") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::BUTTON_RELEASE_E | PointerEvent::BUTTON_1_E;
      } else if (args[2] == "scale_start") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::BUTTON_PRESS_E | PointerEvent::BUTTON_3_E;
      } else if (args[2] == "scale_motion") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::MOTION_E | PointerEvent::BUTTON_3_E;
      } else if (args[2] == "scale_end") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::BUTTON_RELEASE_E | PointerEvent::BUTTON_3_E;
      } else if (args[2] == "reset") {
        modifiers = EventModifiers::SHIFT_E;
        state = PointerEvent::BUTTON_PRESS_E | PointerEvent::BUTTON_2_E;
      }
      p->set_modifiers(modifiers);
      p->set_pointer_state(state);
      event_handle_t event = p;
      EventManager::add_event(event);
    }
  } 
  else if (args[1] == "redraw") 
  {
  } 
  else if (args[1] == "redraw-histo") 
  {
    gui_histo_.reset();
    //! tell the CM2View about the histo
    if (cm2view_)
    {
      if (histo_.get_rep()) 
      {
        event_handle_t event = new SetHistoEvent(histo_->nrrd_, gui_histo_.get(),
                   cm2_target_);
        EventManager::add_event(event);
      }
    }
  } 
  else if (args[1] == "destroygl") 
  { 
    destroy_gl_view();
  } 
  else if (args[1] == "setgl") 
  {
    ASSERT(args.count() == 3);
    if (! cm2view_) 
    {
      TkOpenGLContext *glctx = new TkOpenGLContext(args[2], 0, 512, 256);
      cm2view_ = new CM2View(glctx, get_ctx()->getfullname(), cm2_target_); 
      cm2view_->set_timer_interval(6);
      std::string tid = "CM2View Thread";
      Thread *vt = new Thread(cm2view_, tid.c_str());
      vt->detach(); // runs until thread exits.
      Thread::yield();
      if (widgets_.size()) 
      {
        update_view_widgets();
        update_to_gui();
      } 
      if (histo_.get_rep()) 
      {
        gui_histo_.reset();
        //! tell the CM2View about the histo
        event_handle_t event = new SetHistoEvent(histo_->nrrd_,
                   gui_histo_.get(),
                   cm2_target_);
        EventManager::add_event(event);
      }
      
    }
  } 
  else Module::tcl_command(args, userdata);
}

void
CreateAndEditColorMap2D::post_read()
{
  tcl_unpickle();
  update_view_widgets();
}

void
CreateAndEditColorMap2D::destroy_gl_view()
{
  if (cm2view_)
  {
    cm2view_->set_dead();
    event_handle_t event = new QuitEvent(cm2_target_);
    EventManager::add_event(event);  
    event->wait_signal();
    cm2view_ = 0;
  }
}

void
CreateAndEditColorMap2D::update_view_widgets()
{
  std::vector<SLIVR::CM2Widget*> a;
  if (merge_cmap2_.get_rep()) {
    a = merge_cmap2_->widgets();
  }
  wlock();
  SetWidgetsEvent *e = new SetWidgetsEvent(widgets_, a, cm2_target_);
  wunlock();
  event_handle_t event = e;
  EventManager::add_event(event);
}

void 
CreateAndEditColorMap2D::clear_widgets()
{
  for (size_t i = 0; i < widgets_.size(); ++i) {
    if (widgets_[i]) delete widgets_[i];
  }
  widgets_.clear();
}

void 
CreateAndEditColorMap2D::set_widgets(const std::vector<SLIVR::CM2Widget*> &w,
				     bool updating)
{
  updating_ = updating;
  clear_widgets();
  widgets_ = w;
  need_to_update_gui_ = true;
  force_execute();
}



void
CreateAndEditColorMap2D::save_file(bool save_ppm)
{
  filename_.reset();
  const std::string filename = filename_.get();
  if (filename == "") {
    error("Warning;  No filename provided to CreateAndEditColorMap2D");
    return;
  }
  
  if (save_ppm) {
    CommandEvent *c = new CommandEvent(cm2_target_);
    c->set_command(CM2View::Save_PPM_c + filename);
    event_handle_t e = c;
    EventManager::add_event(e);
  }
  // Open ostream
  PiostreamPtr stream(new BinaryPiostream(filename, Piostream::Write));
  if (stream->error())
    error("Could not open file for writing" + filename);
  else {
    Pio(*stream, sent_cmap2_);
    remark ("Saved ColorMap2 to file: "+filename);
  }
}


void
CreateAndEditColorMap2D::load_file()
{
    // The implementation of this was taken almost directly from
    // NrrdReader Module.  
    filename_.reset();
    std::string fn(filename_.get());
    if(fn == "") {
      error("Please Specify a Transfer Function filename.");
      return;
    }
  
    struct stat buf;
    if(stat(fn.c_str(), &buf) == -1) {
      error(std::string("CreateAndEditColorMap2D error - file not found: '")+fn+"'");
      return;
    }

    const int len = fn.size();
    const std::string suffix(".cmap2");
    // Return if the suffix is wrong
    if (fn.substr(len - suffix.size(), suffix.size()) != suffix) return;

    PiostreamPtr stream = auto_istream(fn, this);
    if (!stream) {
      error("Error reading file '" + fn + "'.");
      return;
    }  
    // read the file.
    ColorMap2Handle icmap = new ColorMap2();
    try {
      Pio(*stream, icmap);
    } catch (...) {
      error("Error loading "+fn);
      icmap = 0;
    }
    stream.reset();
    if (icmap.get_rep()) {
      clear_widgets();
      deep_copy_widgets(icmap->widgets(), widgets_);
      icmap_generation_ = icmap->generation;

      // record the widget state.
      event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::CLEAR_E,
						0, 0, 
						cm2_target_);
      EventManager::add_event(e);
    }
  
    update_to_gui();
    select_widget(-1,0); 
    update_view_widgets();
    force_execute();
}


void
CreateAndEditColorMap2D::presave()
{
  unsigned int i;

  resize_gui();
  update_to_gui(false);

  // Pickle up the tcl states.
  wlock();
  for (i = 0; i < widgets_.size(); i++)
  {
    gui_wstate_[i]->set(widgets_[i]->tcl_pickle());
  }

  const unsigned int ws = widgets_.size();
  wunlock();

  if (ws < gui_name_.size())
  {
    // Delete all of the unused variables.
    for (i = ws; i < gui_name_.size(); i++)
    {
      delete gui_name_[i];
      delete gui_color_r_[i];
      delete gui_color_g_[i];
      delete gui_color_b_[i];
      delete gui_color_a_[i];
      delete gui_wstate_[i];
      delete gui_onstate_[i];
      delete gui_shade_type_[i];
    }

    gui_name_.erase(gui_name_.begin() + ws, gui_name_.end());
    gui_color_r_.erase(gui_color_r_.begin() + ws, gui_color_r_.end());
    gui_color_g_.erase(gui_color_g_.begin() + ws, gui_color_g_.end());
    gui_color_b_.erase(gui_color_b_.begin() + ws, gui_color_b_.end());
    gui_color_a_.erase(gui_color_a_.begin() + ws, gui_color_a_.end());
    gui_wstate_.erase(gui_wstate_.begin() + ws, gui_wstate_.end());
    gui_onstate_.erase(gui_onstate_.begin() + ws, gui_onstate_.end());
    gui_shade_type_.erase(gui_shade_type_.begin() + ws, gui_shade_type_.end());
  }
}


void
CreateAndEditColorMap2D::resize_gui(int n)
{
  gui_num_entries_.reset();
  wlock();
  size_t s = widgets_.size();
  wunlock();
  if (gui_num_entries_.get() == (int)s) return;  
  gui_num_entries_.set(n == -1 ? s : n);
  unsigned int i = 0;
  // Expand the gui elements.
  
  for (i = gui_name_.size(); i < (unsigned int)gui_num_entries_.get(); i++)
  {
    const std::string num = to_string(i);
    gui_name_.push_back(new GuiString(get_ctx()->subVar("name-" + num)));
    gui_color_r_.push_back(new GuiDouble(get_ctx()->subVar(num +"-color-r")));
    gui_color_g_.push_back(new GuiDouble(get_ctx()->subVar(num +"-color-g")));
    gui_color_b_.push_back(new GuiDouble(get_ctx()->subVar(num +"-color-b")));
    gui_color_a_.push_back(new GuiDouble(get_ctx()->subVar(num +"-color-a")));
    gui_wstate_.push_back(new GuiString(get_ctx()->subVar("state-" + num)));
    gui_onstate_.push_back(new GuiInt(get_ctx()->subVar("on-" + num)));
    gui_shade_type_.push_back(new GuiInt(get_ctx()->subVar("shadeType-"+num)));

  }
  // This marker stuff is for TCL, its the last variable created, so
  // its also the last variable written out to the .net script
  // look @ the TCL array ModuleSavedVars(CreateAndEditColorMap2D_0)
  // First: Delete the old variable that marked the end of the variables
  if (end_marker_) 
    delete end_marker_;
  // Second: Create a new marker that marks the end of the variables
  if (i != 0) 
    end_marker_ = new GuiString(get_ctx()->subVar("marker"), "end");
}


void
CreateAndEditColorMap2D::update_to_gui(bool forward)
{
  // Update GUI
  resize_gui();
  wlock();
  gui_selected_widget_.reset();
  
  for (unsigned int i = 0; i < widgets_.size(); i++)
  {
    gui_name_[i]->set(widgets_[i]->get_name());
    Color c(widgets_[i]->get_color());
    gui_color_r_[i]->set(c.r());
    gui_color_g_[i]->set(c.g());
    gui_color_b_[i]->set(c.b());
    gui_color_a_[i]->set(widgets_[i]->get_alpha());
    gui_onstate_[i]->set(widgets_[i]->get_onState());
    gui_shade_type_[i]->set(widgets_[i]->get_shadeType());
    if (widgets_[i]->selected()) gui_selected_widget_.set(i);
  }
  wunlock();

  if (forward) 
  { 
    get_ctx()->reset(); 
    TCLInterface::execute(get_id() + " create_entries"); 
  }
}


void
CreateAndEditColorMap2D::update_from_gui()
{
  get_ctx()->reset();   // Reset GUI vars cache
  resize_gui();   // Make sure we have enough GUI vars to read through
  wlock();
  for (unsigned int i = 0; i < widgets_.size(); i++)
  {
    if (widgets_[i]->get_name() != gui_name_[i]->get()) {
      widgets_[i]->set_name(gui_name_[i]->get());
    }
    Color new_color(gui_color_r_[i]->get(),
		    gui_color_g_[i]->get(),
		    gui_color_b_[i]->get());
    if (widgets_[i]->get_color() != new_color) {
      widgets_[i]->set_color(new_color);
    }
    
    if (fabs(widgets_[i]->get_alpha() - gui_color_a_[i]->get()) > 0.001) {
      widgets_[i]->set_alpha(gui_color_a_[i]->get());
    }

    if (widgets_[i]->get_onState() != gui_onstate_[i]->get()) {
      widgets_[i]->set_onState(gui_onstate_[i]->get());
    }
    if (widgets_[i]->get_shadeType() != gui_shade_type_[i]->get()) {
      widgets_[i]->set_shadeType(gui_shade_type_[i]->get());
    }
  }
  wunlock();
}


void
CreateAndEditColorMap2D::tcl_unpickle()
{
  widgets_.clear();

  gui_num_entries_.reset();
  resize_gui(gui_num_entries_.get());
  for (int i=0; i < gui_num_entries_.get(); i++)
  {
    gui_wstate_[i]->reset();
    if (gui_wstate_[i]->get()[0] == 't')
    {
      widgets_.push_back(new TriangleCM2Widget());
      widgets_[widgets_.size()-1]->tcl_unpickle(gui_wstate_[i]->get());
    }
    else if (gui_wstate_[i]->get()[0] == 'r')
    {
      widgets_.push_back(new RectangleCM2Widget());
      widgets_[widgets_.size()-1]->tcl_unpickle(gui_wstate_[i]->get());
    }
    else if (gui_wstate_[i]->get()[0] == 'i') {
      widgets_.push_back(new ImageCM2Widget());
      widgets_[widgets_.size()-1]->tcl_unpickle(gui_wstate_[i]->get());
    }
    else if (gui_wstate_[i]->get()[0] == 'e') {
      widgets_.push_back(new EllipsoidCM2Widget());
      widgets_[widgets_.size()-1]->tcl_unpickle(gui_wstate_[i]->get());
    }
    else if (gui_wstate_[i]->get()[0] == 'p') {
      widgets_.push_back(new ParaboloidCM2Widget());
      widgets_[widgets_.size()-1]->tcl_unpickle(gui_wstate_[i]->get());
    }
  }

  // Grab colors
  resize_gui();
  update_from_gui();
}

bool
CreateAndEditColorMap2D::select_widget(int widget, int object)
{
  int changed = false;
  if (widget == -1 && object == -1) {
    changed = gui_selected_widget_.changed() || gui_selected_object_.changed();
    widget = gui_selected_widget_.get();
    object = gui_selected_object_.get();
  } else {
    changed = gui_selected_widget_.get() != widget;
    gui_selected_widget_.set(widget);
    gui_selected_object_.set(object);
  }

  for (int i = 0; i < (int)widgets_.size(); i++)
    widgets_[i]->select(i == widget ? object : 0);
  return changed;
}

void
CreateAndEditColorMap2D::set_window_cursor(const std::string& cstr)
{
  TCLInterface::eval(get_ctx()->getfullname() + " set_cursor " + cstr);
}

void
CreateAndEditColorMap2D::execute()
{
  reset_vars();
  ColorMap2Handle icmap = 0;
  NrrdDataHandle h = 0;
  
  if (get_input_handle("Input Colormap", icmap, false)) {
    if (icmap_generation_ != icmap->generation) {
      merge_cmap2_ = icmap;
      icmap_generation_ = icmap->generation;
      update_view_widgets();
    }
  } else if (icmap_generation_ != -1) {
    icmap_generation_ = -1;
    merge_cmap2_ = 0;
    update_view_widgets();
  }
  
  if (get_input_handle("Histogram", h, false)) 
  {
    if (h.get_rep() && h->generation != hist_generation_) 
    {
      hist_generation_ = h->generation;
      if(h->nrrd_->dim != 2 && h->nrrd_->dim != 3) 
      {
        error("Invalid input histogram dimension.");
        return;
      }
      // keep reference to the nrrd.
      histo_ = h; // get_input_handle clears rep to 0 each call.
      
      //! tell the CM2View about the histo
      
      event_handle_t event = new SetHistoEvent(h->nrrd_, gui_histo_.get(),
					       cm2_target_);
      EventManager::add_event(event);
    }
  } 
  else 
  {
    hist_generation_ = -1;
    event_handle_t event = new SetHistoEvent(0, gui_histo_.get(),
					     cm2_target_);
    EventManager::add_event(event);
  }

  // Inform module that execution started
  update_state(Executing);


  wlock();
  if (need_to_update_gui_) {
    update_to_gui();
    need_to_update_gui_ = false;
  }

  //copy the widgets to send.
  std::vector<SLIVR::CM2Widget*>	wds;
  deep_copy_widgets(widgets_, wds);
  std::pair<float, float> vr;
  if (!widgets_.empty())
    vr = widgets_[0]->get_value_range();
  wunlock();

  // If there are input widgets, insert copies in the widget list.
  if (merge_cmap2_.get_rep()) 
  {
    std::vector<SLIVR::CM2Widget*>	mwds;
    deep_copy_widgets(merge_cmap2_->widgets(), mwds);
    std::vector<SLIVR::CM2Widget*>::iterator iter = mwds.begin();
    while (iter != mwds.end()) 
    {
      SLIVR::CM2Widget *w = *iter++;
      wds.push_back(w);
    }
  }

  sent_cmap2_ = new ColorMap2(wds, updating_, 
			      gui_selected_widget_.get(),
			      vr);
  
  send_output_handle("Output Colormap", sent_cmap2_, true);
}





void
CreateAndEditColorMap2D::gui_color_change(GuiArgs &args)
{
  int n = args.get_int(2);
  resize_gui();
  if (n < 0 || n >= gui_num_entries_.get()) return;
  
  gui_color_r_[n]->reset();
  gui_color_g_[n]->reset();
  gui_color_b_[n]->reset();
  gui_color_a_[n]->reset();
  const double a = gui_color_a_[n]->get();
  Color new_color(gui_color_r_[n]->get(), 
		  gui_color_g_[n]->get(),  gui_color_b_[n]->get());
  wlock();
  // record the widget state.
  event_handle_t e = new CM2UndoActionEvent(CM2UndoActionEvent::CHANGE_E, n,
					    widgets_[n]->duplicate(), 
					    cm2_target_);
  EventManager::add_event(e);
  widgets_[n]->set_color(new_color);
  widgets_[n]->set_alpha(a);
  wunlock();
  update_view_widgets();
  force_execute();
}


void
CreateAndEditColorMap2D::gui_toggle_change(GuiArgs &args)
{
  int n = args.get_int(2);
  resize_gui();  // make sure the guivar vector exists
  if (n < 0 || n >= gui_num_entries_.get()) return;
  gui_onstate_[n]->reset();
  int os = gui_onstate_[n]->get();
  wlock();
  widgets_[n]->set_onState(os);  // toggle on/off state.
  wunlock();
  update_view_widgets();
  force_execute();
}

void
CreateAndEditColorMap2D::gui_shade_toggle_change(GuiArgs &args)
{
  int n = args.get_int(2);
  resize_gui();  // make sure the guivar vector exists
  if (n < 0 || n >= gui_num_entries_.get()) return;
  gui_shade_type_[n]->reset();
  int sts = gui_shade_type_[n]->get();
  wlock();
  widgets_[n]->set_shadeType(sts);  // toggle on/off state.
  wunlock();
  update_view_widgets();
  force_execute();
}

} // end namespace SCIRun

