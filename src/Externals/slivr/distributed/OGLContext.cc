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

#include <OGLContext.h>
#include <slivr/Color.h>
#include <slivr/VolumeRenderer.h>
#include <slivr/VideoCardInfo.h>
#include <slivr/ShaderProgramARB.h>
#include <slivr/BBox.h>
#include <slivr/ColorMap.h>
#include <slivr/ColorMap2.h>
#include <slivr/CM2Widget.h>
#include <slivr/Plane.h>


namespace SLIVR_DIST {

using namespace std;


OGLContext::OGLContext() :
  context_(0),
  display_(0),
  window_(0),
  screen_(0)
{}

OGLContext::~OGLContext()
{}

void 
OGLContext::create() 
{
  display_ = XOpenDisplay(NULL);
  if (display_ == NULL)
  {
    cerr << "error opening display" << endl;
    exit(0);
  }
  screen_ = DefaultScreen(display_);
  XSynchronize(display_, 1);

  int attribs[] = {
    GLX_RGBA,
    GLX_RED_SIZE, 8,
    GLX_GREEN_SIZE, 8,
    GLX_BLUE_SIZE, 8,
    GLX_ALPHA_SIZE, 8,
    GLX_DEPTH_SIZE, 24,
    None};

  XVisualInfo *vi = glXChooseVisual(display_, screen_, attribs);
  if (!vi) {
    cerr << "glXChooseVisual failed" << endl;
    return;
  }
    
  int x = 100;
  int y = 100;
  int width = 256;
  int height = 256;


  XSetWindowAttributes  attributes;
  attributes.colormap = XCreateColormap(display_, 
                                        DefaultRootWindow(display_), 
                                        vi->visual, AllocNone);
  attributes.override_redirect = false;
    
  window_ = XCreateWindow(display_, 
                          DefaultRootWindow(display_),
                          x, y, 
                          width, height,
                          3,
                          vi->depth,
                          InputOutput,
                          vi->visual,
                          CWColormap | CWOverrideRedirect,
                          &attributes);

  context_ = glXCreateContext(display_, vi, 0, GL_TRUE);

  if (!context_) {
    cerr << "glXCreateContext failed" << endl;
    return;
  }

  cerr << "context: " << context_ << endl; 
  cerr << "vi: " << vi << endl; 
}


bool 
OGLContext::make_current() 
{
  if (!context_) {
    cerr << "context not created yet, context = " << context_ << endl;
    return false;
  }

  bool success = glXMakeCurrent(display_, window_, 
                                context_);

  if (!success) {
    cerr << "make current failed, do error check..." << endl;
  }
}

void 
OGLContext::draw() {
  make_current();
  
  glDrawBuffer(GL_FRONT);
  glClearColor(0.5, 0.3, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
}

}
