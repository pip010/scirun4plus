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
//
//    File   : Window.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:01:52 2006

#include <Core/Skinner/RenderRegion.h>
#include <Core/Skinner/WindowRedraw.h>
#include <Core/Skinner/Window.h>
#include <Core/Skinner/Variables.h>
#include <Core/Util/ThrottledRunnable.h>
#include <Core/Util/Assert.h>
#include <Core/Util/StringUtil.h>
#include <Core/Thread/Semaphore.h>
#include <Core/Util/Environment.h>
#include <Core/Events/BaseEvent.h>
#include <Core/Math/MinMax.h>
#include <iostream>
#include <sci_gl.h>
#include <sci_defs/x11_defs.h>
#include <sci_defs/bits_defs.h>
#include <Core/Geom/GeomResourceManager.h>

#if defined(_WIN32)
#  include <Core/Geom/Win32OpenGLContext.h>
#  include <Core/Events/Win32EventSpawner.h>
#elif defined(__APPLE__) && !defined(SCI_64BITS)
#  include <Core/Geom/OSXOpenGLContext.h>
#  include <Core/Events/OSXEventSpawner.h>
#else
#  include <Core/Geom/X11OpenGLContext.h>
#  include <Core/Events/X11EventSpawner.h>
#endif


namespace SCIRun {
namespace Skinner {

GLWindow::GLWindow(Variables *vars) :
  Parent(vars),
  lock_(vars->get_id().c_str()),
  width_(vars,"width",640),
  height_(vars,"height",480),
  fps_(vars,"GLWindow::fps",480),
  frames_(0),
  last_time_fps_reset_(0.0),
  posx_(Var<int>(vars,"posx",50)()),
  posy_(Var<int>(vars,"posy",50)()),
  border_(Var<bool>(vars,"border",true)()),
  context_(0),
  own_context_(false),
  spawner_runnable_(0),
  spawner_thread_(0),
  draw_runnable_(0),
  draw_thread_(0),
  redrawables_(),
  force_redraw_(false),
  png_buf_(0),
  png_rows_(0),
  png_num_(0),
  do_png_(sci_getenv("SKINNER_MOVIE_BASE_FILENAME"))
{
  timer_.start();
  REGISTER_CATCHER_TARGET(GLWindow::close);
  REGISTER_CATCHER_TARGET(GLWindow::mark_redraw);
  REGISTER_CATCHER_TARGET(GLWindow::redraw_drawable);
  REGISTER_CATCHER_TARGET(GLWindow::RenderRegion_Maker);
}

bool
GLWindow::create_window(void *data)
{
  if (context_) return true;

#if defined(_WIN32)
  Win32GLContextRunnable* cr =
    new Win32GLContextRunnable(get_id(), 0, posx_, posy_,
                                  (unsigned)width_,(unsigned)height_,
                                  border_);
  spawner_runnable_ = cr;
#elif defined(__APPLE__) && !defined(SCI_64BITS)
  OSXOpenGLContext *context =
    new OSXOpenGLContext(posx_, posy_,
                         (unsigned)width_,(unsigned)height_,
                         border_);
  ASSERT(context);
  //      spawner_runnable_ =
  new OSXEventSpawner(get_id(), context->window_);
  context_ = context;
  own_context_ = true;
#else
  X11OpenGLContext* context =
    new X11OpenGLContext(0, posx_, posy_,
                         (unsigned)width_,(unsigned)height_,
                         border_, 1, data);
  ASSERT(context);
  spawner_runnable_ =
    new X11EventSpawner(get_id(), context->display_, context->window_);
  context_ = context;
  own_context_ = true;
#endif

  const string spawn_thread_name = get_id() + " Event Spawner";
  if (spawner_runnable_) {
    spawner_thread_ = new Thread(spawner_runnable_, spawn_thread_name.c_str());
  }

#ifdef _WIN32
  // this waits for the context to be created,
  // and needs to happen after the thread spawns
  context_ = cr->getContext();;
  own_context_ = true;
  width_ = context_->width();
  height_ = context_->height();
#endif

  const string draw_thread_name = get_id() + " Redraw";
  draw_runnable_ = new WindowRedraw(this, 120.0);
  draw_thread_ = new Thread(draw_runnable_, draw_thread_name.c_str());
  return true;
}


bool
GLWindow::attach_to_context(OpenGLContext *context)
{
  if (context_) return true;

  context_ = context;
  own_context_ = false;
  width_ = context_->width();
  height_ = context_->height();

  const string tname = get_id() + " Redraw";
  draw_runnable_ = new WindowRedraw(this, 120.0);
  draw_thread_ = new Thread(draw_runnable_, tname.c_str(),
                            0, Thread::NotActivated);
  draw_thread_->setStackSize(1024*1024);
  draw_thread_->activate(false);
  return true;
}


GLWindow::~GLWindow()
{
  timer_.stop();
  if (spawner_runnable_) spawner_runnable_->quit();
  if (spawner_thread_) spawner_thread_->join();
  if (draw_runnable_) draw_runnable_->quit();
  if (draw_thread_) draw_thread_->join();

  // Context must be deleted after event spawer is done.
  if (own_context_) delete context_;
  context_ = 0;
  own_context_ = false;

  throw_signal("GLWindow::Destructor");
}


int
GLWindow::get_signal_id(const string &name) const
{
  if (name == "GLWindow::Destructor") return 1;
  if (name == "GLWindow::Destroy") return 2;
  return 0;
}


BaseTool::propagation_state_e
GLWindow::process_event(event_handle_t &event)
{
  WindowEvent *window = dynamic_cast<WindowEvent *>(event.get_rep());
  if (window && (window->get_window_state() & WindowEvent::DESTROY_E))
  {
    throw_signal("GLWindow::Destroy");
    return CONTINUE_E;
  }

  warp_pointer(event);

  // Subdraws with a double buffered context do not work.  1) Drawing
  // to the back results in 1 frame out of date and hairpulling
  // debugging of redraws (because you force redraws and they revert
  // themselves).  2) The back buffer isn't guaranteed to remain
  // intact, so using it as a pristine copy of the front fails.  3)
  // Doing subdraws in the front buffer is better but results in
  // flicker and drawing artifacts.  For now just redraw it all.
  force_redraw_ = true;

  bool redraw = (context_ && window &&
                 window->get_window_state() & WindowEvent::REDRAW_E);
  bool subdraw = redraw && redrawables_.size() && !force_redraw_;

  // If we don't lock the entire redraw then NVidia drivers will die out
  // on 2D colormap framebuffer creation.
  if (redraw)
  {
    ASSERT(context_);

    context_->lock();
    if (!context_->make_current()) {
      context_->unlock();
      return CONTINUE_E;
    }

    GeomResourceManager::delete_pending_objects();

    const int vpw = context_->width();
    const int vph = context_->height();

    width_ = vpw;
    height_ = vph;

    set_region(RectRegion(0.0, 0.0, double(vpw), double(vph)));

    glViewport(0, 0, vpw, vph);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glScaled(2.0, 2.0, 2.0);
    glTranslated(-.5, -.5, -.5);
    glScaled(1.0/double(vpw), 1.0/double(vph), 1.0);
    CHECK_OPENGL_ERROR();

    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDrawBuffer(GL_BACK);
    glClearColor(0,0,0,0);
    glClearDepth(1.0);
    GLbitfield mask = GL_DEPTH_BUFFER_BIT;
    if (!subdraw) {
      mask |= GL_COLOR_BUFFER_BIT;
    }
    glClear(mask);

    CHECK_OPENGL_ERROR();
  }

  if (subdraw) {
    redrawables_t::reverse_iterator iditer = redrawables_.rbegin();
    redrawables_t::reverse_iterator idlast = redrawables_.rend();
    for (; iditer != idlast; ++iditer) {
      iditer->second->process_event(event);
    }
    redrawables_.clear();

  } else {
    if (redraw && force_redraw_) {
      redrawables_.clear();
      force_redraw_ = false;
    }

    for (Drawables_t::iterator child = children_.begin();
         child != children_.end(); ++child) {
      (*child)->set_region(get_region());
      (*child)->process_event(event);
    }

  }

  if (redraw)
  {
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    CHECK_OPENGL_ERROR();

    // The function save_png was removed because the code was commented out,
    // and was not functional.
    // TODO: review this code and determine what else to remove.
    context_->swap();
    //save_png();
    context_->release();

    context_->unlock();

    frames_++;
    double time = timer_.time();
    double delta = time - last_time_fps_reset_;
    if (frames_ > 30 || delta > 1.0) {
      fps_ = frames_ / Max(delta, 0.001);
      last_time_fps_reset_ = time;
      frames_ = 0;
    }
  }

  return CONTINUE_E;
}


BaseTool::propagation_state_e
GLWindow::close(event_handle_t &)
{
  QuitEvent *qe = new QuitEvent();
  qe->set_target(get_id());
  EventManager::add_event(qe);
  return STOP_E;
}


BaseTool::propagation_state_e
GLWindow::mark_redraw(event_handle_t& /*event*/)
{
  force_redraw_ = true;
  EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E, get_name()));
  return CONTINUE_E;
}


BaseTool::propagation_state_e
GLWindow::redraw_drawable(event_handle_t &signalh)
{
  ASSERT(dynamic_cast<Signal *>(signalh.get_rep()));
  Signal *signal = (Signal *)(signalh.get_rep());
  ASSERT(dynamic_cast<Drawable *>(signal->get_signal_thrower()));
  Drawable *drawable = (Drawable *)(signal->get_signal_thrower());

  Var<string> id(signal->get_vars(), "region");
  if (id.exists() && id().length()) {
    for (unsigned int i = 0; i < render_regions_.size(); ++i) {
      if (ends_with(render_regions_[i]->get_id(), id())) {
        drawable = render_regions_[i];
        break;
      }
    }
  }

  Var<bool> visible(drawable->get_vars(), "visible", true);
  if (visible &&
      redrawables_.find(drawable->get_id()) == redrawables_.end()) {
    redrawables_.insert(make_pair(drawable->get_id(), drawable));
    EventManager::add_event(new WindowEvent(WindowEvent::REDRAW_E, get_id()));
  }
  return CONTINUE_E;
}

BaseTool::propagation_state_e
GLWindow::RenderRegion_Maker(event_handle_t &maker_signal)
{
  render_regions_.push_back
    (construct_child_from_maker_signal<RenderRegion>(maker_signal));
  return STOP_E;
}


// This procedure takes care of inverting the y coordinates on X11 systems
BaseTool::propagation_state_e
GLWindow::warp_pointer(event_handle_t &event)
{
  if (!event->is_pointer_event()) {
    return STOP_E;
  }


#if defined(_WIN32)
  // do nothing
#elif defined(__APPLE__)
  // do nothing
#else
  if (dynamic_cast<X11OpenGLContext *>(context_))
  {
    ASSERT(dynamic_cast<PointerEvent *>(event.get_rep()));
    PointerEvent *pointer = (PointerEvent *)(event.get_rep());
    if ((pointer->get_pointer_state() & PointerEvent::MOTION_E) &&
        EventManager::trailfile_is_playing())
    {
      const char *trailmode = sci_getenv("SCIRUN_TRAIL_MODE");
      ASSERT(dynamic_cast<X11OpenGLContext *>(context_));
      X11OpenGLContext *x11 = (X11OpenGLContext *)(context_);
      if (x11 && trailmode) {
        Window src_win = x11->window_;

        // if the second letter of SCIRUN_TRAIL_MODE is G,
        // Grab the pointer from the root window from the user
        if (trailmode[1] == 'G') {
          src_win = None;
        }

        context_->lock();
        XWarpPointer(x11->display_, src_win, x11->window_,
                     0,0, x11->width(), x11->height(),
                     pointer->get_x(), pointer->get_y());
        context_->unlock();
      }
    }

    pointer->set_y(context_->height() - 1 - pointer->get_y());
  }
#endif

  return CONTINUE_E;
}


}
}
