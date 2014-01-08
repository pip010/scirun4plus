//
//  For more information, please see: http://software.sci.utah.edu
//
//  The MIT License
//
//  Copyright (c) 2009 Scientific Computing and Imaging Institute,
//  University of Utah.
//
//  License for the specific language governing rights and limitations under
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
//    File   : RenderBuffer.cc
//    Author : McKay Davis
//    Date   : Mon Feb 12 00:18:25 2007

#include <Core/Util/StringUtil.h>
#include <Core/Skinner/Variables.h>
#include <Core/Skinner/RenderBuffer.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Thread.h>
#include <Core/Thread/Semaphore.h>
#include <Core/Util/Environment.h>

#include <sci_gl.h>
#include <Core/Geom/OpenGLContext.h>
#include <Core/Geom/X11OpenGLContext.h>


namespace SCIRun {

static  int gcount = 0;
bool
check_fbo_status ()
{
  cerr << gcount++ << " ";
  //----------------------
  // Framebuffer Objects initializations
  //----------------------
  GLuint status = glCheckFramebufferStatusEXT( GL_FRAMEBUFFER_EXT );

  switch( status ) {
  case GL_FRAMEBUFFER_COMPLETE_EXT:
    std::cerr << " GL_FRAMEBUFFER_COMPLETE_EXT \n";
    return true;

  case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
    std::cerr << " GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT \n";
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
    std::cerr << " GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT \n";
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
    std::cerr << " GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT \n";
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
    std::cerr << " GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT \n";
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
    std::cerr << " GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT \n";
    return false;

  case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
    std::cerr << " GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT \n";
    return false;

  case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
    std::cerr << " GL_FRAMEBUFFER_UNSUPPORTED_EXT \n";
    return false;

  default:
    std::cerr << " GL_FRAMEBUFFER_UNKNOWN " << status << std::endl;
    return false;
  }
}


class Skinner::ThreadedRenderBuffer : public Runnable {
public:
  Skinner::RenderBuffer *rb_;
  GLuint fbo_id_[2];
  GLuint depth_id_[2];
  GLuint tex_id_[2];
  unsigned int visible_;
  unsigned int pass_;
  Semaphore sem_;
  OpenGLContext *ctx_;
  ThreadedRenderBuffer(Skinner::RenderBuffer *rb) :
    Runnable(true),
    rb_(rb),
    fbo_id_(),
    depth_id_(),
    tex_id_(),
    visible_(0),
    pass_(0),
    sem_("threadedrunderbuffer",0),
    ctx_(0)
  {
    fbo_id_[0] = fbo_id_[1] = 0;
    depth_id_[0] = depth_id_[1] = 0;
    tex_id_[0] = tex_id_[1] = 0;

    if (!ctx_) {
#ifdef HAVE_X11
      ctx_ = new SCIRun::X11OpenGLContext(0, 0, 0, 1, 1, 0, 0);
#endif
    }
  }

  void
  run() {

    //    cerr << "running\n";
    if (!ctx_) return;
    while (1) {


      sem_.down();
      pass_ = (pass_+1)%2;
      //cerr << "rending id = " << tex_id_ << endl;
      if (!ctx_->make_current()) {
        cerr << "ThreadedRenderBuffer cannot make current\n";
        return;
      }

      int wid = 256;
      int hei = 256;
      gcount = 0;
      if (!fbo_id_[pass_]) {

        // Create binding ids
        glGenFramebuffersEXT( 1, fbo_id_+pass_ );
        glGenRenderbuffersEXT(1, depth_id_+pass_);
        //CHECK_OPENGL_ERROR("FBO bindings");

        // Generate depth buffer storage
        glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depth_id_[pass_]);
        glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT24,
                                 wid, hei);
        //CHECK_OPENGL_ERROR("FBO depth buffer create");
        //        check_fbo_status();

        //        cerr << "  fbo_id_: " << fbo_id_
        //     << "  depth_id_: " << depth_id_
        //     << "  tex_id_: " << tex_id_ << std::endl;


        // Bind Frame Buffer Object
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_id_[pass_]);
        //CHECK_OPENGL_ERROR("FBO bind");

        // Bind Frame Buffer Object to render to Texture Object
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                                  GL_TEXTURE_2D, tex_id_[pass_], 0);
        //CHECK_OPENGL_ERROR("FBO bind to texture");

        // Bind Frame Buffer Object to render to Depth Buffer
        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                                     GL_RENDERBUFFER_EXT, depth_id_[pass_]);
        //CHECK_OPENGL_ERROR("FBO bind to depth");

      }

      // Bind Frame Buffer Object
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_id_[pass_]);
      //CHECK_OPENGL_ERROR("FBO bind");


      // Bind Frame Buffer Object to render to Texture Object
      glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT,
                                GL_TEXTURE_2D, tex_id_[pass_], 0);
      //CHECK_OPENGL_ERROR("FBO bind to texture");

#if 0
      // Bind Frame Buffer Object to render to Depth Buffer
      glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,
                                   GL_RENDERBUFFER_EXT, depth_id_[pass_]);
      //CHECK_OPENGL_ERROR("FBO bind to depth");

#endif

      check_fbo_status();
      glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      //CHECK_OPENGL_ERROR("gldrawbuffer FBO");

      GLint vp[4];
      glGetIntegerv(GL_VIEWPORT, vp);
      glViewport(0,0,wid,hei);
      Skinner::RectRegion reg(0,0,wid,hei);
      glDisable(GL_BLEND);

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();

      glScaled(2.0, 2.0, 2.0);
      glTranslated(-.5, -.5, -.5);
      glScaled(1.0/double(wid), 1.0/double(hei), 1.0);
      //CHECK_OPENGL_ERROR( "glScaled" );

      glDisable(GL_CULL_FACE);
      glEnable(GL_BLEND);
      glDisable(GL_BLEND);
      //      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glClearColor(0.,0.,0.,0.0);
      glClearDepth(1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      event_handle_t event = new WindowEvent(WindowEvent::REDRAW_E);
      for (Drawables_t::iterator citer = rb_->children_.begin();
           citer != rb_->children_.end(); ++citer)
      {
        (*citer)->set_region(reg);
        (*citer)->process_event(event);
      }

      //      rb_->Parent::process_event(event);

      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      while(sem_.tryDown()) {
        cerr << "SEMTRYDOWN\n";
      }


      //      ctx_->release();
    }
  }
};



Skinner::RenderBuffer::RenderBuffer(Variables *vars)
  : Skinner::Parent(vars),
    runner_(0),
    ids_(),
    flipflop_(0)
{
  runner_ = new ThreadedRenderBuffer(this);
  new Thread(runner_, "ThreadedRender");
}


Skinner::RenderBuffer::~RenderBuffer()
{
}


BaseTool::propagation_state_e
Skinner::RenderBuffer::process_event(event_handle_t &event)
{

  WindowEvent *window = dynamic_cast<WindowEvent *>(event.get_rep());
  bool redraw = window && window->get_window_state() == WindowEvent::REDRAW_E;

  if (!redraw) {
    return Parent::process_event(event);
  }

  RectRegion reg = get_region();

  int wid = Ceil(reg.width());
  int hei = Ceil(reg.height());

  wid = 256;
  hei = 256;

  gcount = 1;

  // Generate texture ids
  for (int i = 0; i < 2; ++i) {
    if (!runner_->tex_id_[i]) {
      glGenTextures(1, runner_->tex_id_+i);
      glBindTexture(GL_TEXTURE_2D, runner_->tex_id_[i]);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, wid, hei,
                   0, GL_RGBA, GL_FLOAT, 0);
    }
  }

  runner_->sem_.up();

  // cerr << "drawing id = " << ids_[flipflop_%2] << std::endl;
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, runner_->tex_id_[(runner_->pass_+1)%2]);
  glColor4f(1.0, 0.0, 1.0, 1.0);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  glLoadIdentity();

  glBegin(GL_QUADS);

  glTexCoord2f(0.0, 0.0);
  glVertex3d(reg.x1(), reg.y1(), 0.0);

  glTexCoord2f(1.0, 0.0);
  glVertex3d(reg.x2(), reg.y1(), 0.0);

  glTexCoord2f(1.0, 1.0);
  glVertex3d(reg.x2(), reg.y2(), 0.0);

  glTexCoord2f(0.0, 1.0);
  glVertex3d(reg.x1(), reg.y2(), 0.0);

  glEnd();

  glMatrixMode(GL_TEXTURE);
  glPopMatrix();

  glDisable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //CHECK_OPENGL_ERROR( "Skinner::RenderBuffe::process_event() - end" );

  return CONTINUE_E;
} // end redraw()


} // end namespace SCIRun
