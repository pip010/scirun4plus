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
//    File   : X11OpenGLContext.cc
//    Author : McKay Davis
//    Date   : May 30 2006

#include <sci_defs/x11_defs.h>
#include <sci_defs/opengl_defs.h>

#if defined(HAVE_X11)

#include <Core/Math/MiscMath.h>
#include <Core/Geom/X11OpenGLContext.h>
#include <slivr/ShaderProgramARB.h>
#include <Core/Geom/X11Lock.h>
#include <Core/Util/StringUtil.h>
#include <Core/Datatypes/Color.h>
#include <Core/Exceptions/InternalError.h>

#include <Core/Thread/Mutex.h>
#include <Core/Thread/Thread.h>
#include <Core/Util/Assert.h>
#include <Core/Util/Environment.h>
#include <sci_glx.h>
#include <iostream>
#include <set>

// X11 includes
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xmu/StdCmap.h>


using namespace SCIRun;


extern SCIRun::Mutex dl_lock_;
  
std::vector<int> X11OpenGLContext::valid_visuals_ = std::vector<int>();
#ifdef HAVE_OPENGL
GLXContext X11OpenGLContext::first_context_ = 0;
#endif

X11OpenGLContext::X11OpenGLContext(int visual, 
                                   int x, 
                                   int y,
                                   unsigned int width, 
                                   unsigned int height,
                                   bool border,
                                   bool show,
                                   void *data) : 
  OpenGLContext(),
  mutex_("GL lock")
{
  // Sepeate functions so we can set gdb breakpoints in constructor
  create_context(visual, x, y, width, height, border, show, data);
}



void
X11OpenGLContext::create_context(int visual, 
                                 int x, 
                                 int y,
                                 unsigned int width, 
                                 unsigned int height,
                                 bool border,
                                 bool show,
                                 void *data)
{
#ifdef HAVE_OPENGL
  X11Lock::lock();

  display_ = XOpenDisplay(sci_getenv("DISPLAY"));
  XSync(display_, False);

  screen_ = DefaultScreen(display_);
  root_window_ = DefaultRootWindow(display_);

  if (valid_visuals_.empty())
    listvisuals();
  ASSERT(!valid_visuals_.empty());

  visual = Clamp(visual, 0, (int)valid_visuals_.size()-1);
  visualid_ = valid_visuals_[visual];
  ASSERT(visualid_);

  int n;
  XVisualInfo temp_vi;
  temp_vi.visualid = visualid_;
  vi_ = XGetVisualInfo(display_, VisualIDMask, &temp_vi, &n);
  if(!vi_ || n != 1) {
    throw ("Cannot find Visual ID #" + to_string(visualid_) + 
           std::string(__FILE__)+to_string(__LINE__));
    X11Lock::unlock();
  }
    

  attributes_.colormap = XCreateColormap(display_, root_window_, 
                                         vi_->visual, AllocNone);
  attributes_.override_redirect = border ? false : true;
  //  unsigned int valuemask = (CWX | CWY | CWWidth | CWHeight | 
  //                            CWBorderWidth | CWSibling | CWStackMode);

  if (data) {
    root_window_ = (Window)data;
  }
  
  window_ = XCreateWindow(display_, 
                          root_window_,
                          x, y, 
                          width, height,
                          0,
                          vi_->depth,
                          InputOutput,
                          vi_->visual,
                          //                          0, //valuemask,
                          CWColormap | CWOverrideRedirect,
                          &attributes_);
  
  if (!window_) {
    X11Lock::unlock();
    throw ("Cannot create X11 Window " + 
           std::string(__FILE__)+to_string(__LINE__));
  }

  if (show) {
    XMapRaised(display_, window_);
    XMoveResizeWindow(display_, window_, x, y, width, height);
  }
  XSync(display_, False);

  context_ = glXCreateContext(display_, vi_, first_context_, 1);
  if (!context_) {
    X11Lock::unlock();
    throw ("Cannot create GLX Context" + 
           std::string(__FILE__)+to_string(__LINE__));
  }


  X11Lock::unlock();

  if (!first_context_) {
    first_context_ = context_;
    make_current();
    dl_lock_.lock();
    SLIVR::ShaderProgramARB::init_shaders_supported();
    dl_lock_.unlock();
    release();
  }

  width_ = width;
  height_ = height;
#else
  width_ = width;
  height_ = height;  
#endif
}



X11OpenGLContext::~X11OpenGLContext()
{
#ifdef HAVE_OPENGL
  release();
  X11Lock::lock();
  XSync(display_, False);

  glXDestroyContext(display_, context_);

  XSync(display_, False);

  XDestroyWindow(display_,window_);

  XSync(display_, False);

  // This crashes Seg3D on intel when the splash screen is closed.
  // It seems to blow up the main window's context at the same time.
  //XCloseDisplay(display_);

  X11Lock::unlock();
#endif
}


bool
X11OpenGLContext::make_current()
{
#ifdef HAVE_OPENGL
  ASSERT(context_);
  bool result = true;
  X11Lock::lock();
  result = glXMakeCurrent(display_, window_, context_);
  X11Lock::unlock();
  if (!result)
  {
    std::cerr << "X11 GL context failed make current." << std::endl;
  }

  return result;
#else
  return true;
#endif  
}


void
X11OpenGLContext::release()
{
#ifdef HAVE_OPENGL
  X11Lock::lock();
  glXMakeCurrent(display_, None, NULL);
  X11Lock::unlock();
#endif
}


int
X11OpenGLContext::width()
{
#ifdef HAVE_OPENGL
  X11Lock::lock();

  // TODO: optimize out to configure events
  XWindowAttributes attr;
  XGetWindowAttributes(display_, window_, &attr);
  width_ = attr.width;
  height_ = attr.height;
  X11Lock::unlock();

  return width_;
#else
  return width_;
#endif  
}


int
X11OpenGLContext::height()
{
  return height_;
}


void
X11OpenGLContext::swap()
{  
#ifdef HAVE_OPENGL
  X11Lock::lock();
  glXSwapBuffers(display_, window_);
  X11Lock::unlock();
#endif
}

bool
X11OpenGLContext::has_shaders()
{
#ifdef HAVE_OPENGL
  return(SLIVR::ShaderProgramARB::shaders_supported());
#else
  return (false);
#endif
}

bool
X11OpenGLContext::initialized()
{
#ifdef HAVE_OPENGL
  return(SLIVR::ShaderProgramARB::initialized());
#else
  return (true);
#endif
}

#define GETCONFIG(attrib) \
if(glXGetConfig(display, &vinfo[i], attrib, &value) != 0){\
  std::cerr << "Error getting attribute: " << #attrib << std::endl; \
  return; \
}


void
X11OpenGLContext::listvisuals()
{
#ifdef HAVE_OPENGL
  valid_visuals_.clear();
  std::vector<std::string> visualtags;
  std::vector<int> scores;
  int nvis;
  Display *display;
  display = XOpenDisplay((char *)0);
  int screen = DefaultScreen(display);

  XVisualInfo* vinfo = XGetVisualInfo(display, 0, NULL, &nvis);
  if(!vinfo)
  {
    std::cerr << "XGetVisualInfo failed";
    return;
  }
  for(int i=0;i<nvis;i++)
  {
    int score=0;
    int value;
    GETCONFIG(GLX_USE_GL);
    if(!value)
      continue;
    GETCONFIG(GLX_RGBA);
    if(!value)
      continue;
    GETCONFIG(GLX_LEVEL);
    if(value != 0)
      continue;
    if(vinfo[i].screen != screen)
      continue;
    char buf[20];
    sprintf(buf, "id=%02x, ", (unsigned int)(vinfo[i].visualid));
    valid_visuals_.push_back(vinfo[i].visualid);
    std::string tag(buf);
    GETCONFIG(GLX_DOUBLEBUFFER);
    if(value)
    {
      score+=200;
      tag += "double, ";
    }
    else
    {
      tag += "single, ";
    }
    GETCONFIG(GLX_STEREO);
    if(value)
    {
      score+=1;
      tag += "stereo, ";
    }
    tag += "rgba=";
    GETCONFIG(GLX_RED_SIZE);
    tag+=to_string(value)+":";
    score+=value;
    GETCONFIG(GLX_GREEN_SIZE);
    tag+=to_string(value)+":";
    score+=value;
    GETCONFIG(GLX_BLUE_SIZE);
    tag+=to_string(value)+":";
    score+=value;
    GETCONFIG(GLX_ALPHA_SIZE);
    tag+=to_string(value);
    score+=value;
    GETCONFIG(GLX_DEPTH_SIZE);
    tag += ", depth=" + to_string(value);
    score+=value*5;
    GETCONFIG(GLX_STENCIL_SIZE);
    score += value * 2;
    tag += ", stencil="+to_string(value);
    tag += ", accum=";
    GETCONFIG(GLX_ACCUM_RED_SIZE);
    tag += to_string(value) + ":";
    GETCONFIG(GLX_ACCUM_GREEN_SIZE);
    tag += to_string(value) + ":";
    GETCONFIG(GLX_ACCUM_BLUE_SIZE);
    tag += to_string(value) + ":";
    GETCONFIG(GLX_ACCUM_ALPHA_SIZE);
    tag += to_string(value);

    tag += ", score=" + to_string(score);
    
    visualtags.push_back(tag);
    scores.push_back(score);
  }
  for(int i=0;i < int(scores.size())-1;i++)
  {
    for(int j=i+1;j< int(scores.size());j++)
    {
      if(scores[i] < scores[j])
      {
        std::swap(scores[i], scores[j]);
	      std::swap(visualtags[i], visualtags[j]);
	      std::swap(valid_visuals_[i], valid_visuals_[j]);
      }
    }
  }
#endif
}

#endif
